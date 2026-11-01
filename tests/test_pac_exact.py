#!/usr/bin/env python3
"""Validation for card #728: td `.pac` on boost / buck-boost decks with a capacitor at the switch node vs an exact reference.

The six decks of the card #713 arbitration (tools/pac_exact_decks.py: the BOOSTC_580 / BB_583 decks of tests/test_pac_charge.py; boost with
1 nF / 100 pF / 3 nF at the switch node, chargeless or cj 100 pF + tt 2 ns diode; inverting buck-boost with 1 nF, chargeless or charged) are
treated as fast capacitor nodes (card #580: routed jump + continuation).  On s144 they read the phase low against the exact reference
(bst_cl_csw -0.098 deg, bst_q_csw3n -0.156 deg at 100 kHz) and the buck-boost gain high (bb_cl_csw +0.021 dB at 20 kHz).  Card #728 found
six defects in what the continuation gives and fixed them (details in the commit message of the fix); PSS_SALT_NO728=1 restores s144.

    NGSPICE=<prefix>/bin/ngspice python3 tests/test_pac_exact.py [--no-tool]
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
Checks (N/M pass is printed at the end):
  1-6  each deck, 200 points, reltol 1e-4, 20/50/100 kHz: within 0.01 dB / 0.05 deg of the exact values pinned below
  7    the same six decks at 1000 points (the result does not hang on the step): within 0.01 dB / 0.05 deg of the exact values
  8    red proof: PSS_SALT_NO728=1 on the same binary: bst_q_csw3n at 100 kHz more than 0.10 deg and bb_cl_csw at 20 kHz more than
       0.015 dB off the exact values (s144: -0.156 deg, +0.021 dB)
  10   the event buffer never drops a part silently: with PSS_SALT_EVCAP=6 (the buffer lowered from 32 to 6) bst_q_csw3n at 100 kHz,
       1000 points, is REFUSED with the overflow message (its continuation steps carry 6 parts since #728; the s144 buffer of 8 dropped
       anything beyond without a message)
  9    the tool itself (needs numpy + scipy: PAC_EXACT_PYTHON=<python>, else python3 if it imports scipy, else COULD-NOT-RUN; --no-tool
       skips it): tools/pac_exact.py at rtol 1e-9, amplitudes 2 and 1 mV, bst_cl_csw100p at 100 kHz and bb_q_csw at 50 kHz, within
       0.002 dB / 0.01 deg of the pinned values (about 50 s)

Exact values: tools/pac_exact.py (card #713, ~/Work/data/ngspice_augmentation/713b/ex): the circuit as an ODE with ngspice's level-1 diode
equations, the ideal switch's instants found exactly (brentq on the known control), segment-wise Radau at rtol 1e-12, periodic steady
state of the LCM period by shooting, the Fourier coefficient at f as a quadrature state, symmetric difference in the amplitude (4, 2, 1,
0.5 mV) and Richardson in A^2 on 0.5 / 1 mV.  Error bar (spread over the amplitudes and rtol 1e-10 / 1e-12): <= 1.8e-5 dB / 6.3e-4 deg.
An ngspice cross-check with the exact switch instants as breakpoints agreed to 0.0005 dB / 0.0033 deg at 50 / 100 kHz.
Decks and outputs live under ~/Work/data/ngspice_augmentation/728/test (PAC_EXACT_OUT=<dir> redirects them), never the repository.
"""
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.dont_write_bytecode = True   # nothing written into the repository
sys.path.insert(0, os.path.join(ROOT, "tools"))
from pac_exact_decks import DECKS   # noqa: E402  (deck texts only; no numpy needed)

NGSPICE = os.environ.get("NGSPICE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice")
OUT = os.path.abspath(os.path.expanduser(os.environ.get("PAC_EXACT_OUT", "~/Work/data/ngspice_augmentation/728/test")))
BOUND_DB, BOUND_DEG = 0.01, 0.05

# tools/pac_exact.py, rtol 1e-12, Richardson (A^2) on 0.5 / 1 mV: deck -> {f: (dB, deg)} of v(vout) / Vctrl
EXACT = {
    "bst_cl_csw": {20000: (-1.68473, 169.0391), 50000: (-16.52507, 146.3780), 100000: (-25.15105, 130.1531)},
    "bst_q_csw": {20000: (-1.69859, 168.9931), 50000: (-16.52836, 146.2880), 100000: (-25.14476, 130.0699)},
    "bst_cl_csw100p": {20000: (-1.67947, 168.8294), 50000: (-16.54156, 146.5110), 100000: (-25.18426, 130.4368)},
    "bst_q_csw3n": {20000: (-1.71048, 169.4373), 50000: (-16.49096, 145.9882), 100000: (-25.07030, 129.4447)},
    "bb_cl_csw": {20000: (5.62787, 1.0996), 50000: (-10.63227, -11.9666), 100000: (-21.56521, -24.3629)},
    "bb_q_csw": {20000: (5.61128, 1.1092), 50000: (-10.64466, -12.0097), 100000: (-21.57187, -24.4432)},
}

results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"{state:14s} {name}\n{'':14s}   {detail}", flush=True)


def dphase(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def pac(deck, f, points=200, reltol="1e-4", env_extra=None, tag=""):
    """One td .pac point (the deck form of the card #713 arbitration); returns (dB, deg) of v(vout)."""
    os.makedirs(OUT, exist_ok=True)
    name = f"{deck}_{f}_p{points}_r{reltol}{tag}"
    out = os.path.join(OUT, name + ".txt")
    if os.path.exists(out):
        os.remove(out)
    d = (DECKS[deck][4] + ".options reltol=%s\n" % reltol
         + '.pac 500e3 100e-6 vout %d 10 150 1e-6 lin 1 %g %g method="shoot" solver="td"\n' % (points, f, f)
         + ".control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n"
         + f"wrdata {out} v(vout)\nquit\n.endc\n.end\n")
    path = os.path.join(OUT, name + ".cir")
    open(path, "w").write(d)
    env = dict(os.environ, OMP_NUM_THREADS="2")
    if env_extra:
        env.update(env_extra)
    try:
        p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=3G", NGSPICE, "-b", path],
                           capture_output=True, text=True, timeout=900, env=env, cwd=OUT)
    except (OSError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{NGSPICE}: {ex}")
    open(os.path.join(OUT, name + ".log"), "w").write(p.stdout + p.stderr)
    if "REFUS" in p.stdout + p.stderr:
        raise CouldNotRun("refused: " + [l for l in (p.stdout + p.stderr).split("\n") if "REFUS" in l][0][:200])
    if not os.path.exists(out):
        raise CouldNotRun(f"no output for {name} (rc {p.returncode})")
    t = open(out).read().split("\n")[1].split()
    z = complex(float(t[1]), float(t[2]))
    return 20 * math.log10(abs(z)), math.degrees(math.atan2(z.imag, z.real))


def errors(deck, points=200, env_extra=None, tag=""):
    """{f: (dB error, deg error)} of the engine against the pinned exact values"""
    e = {}
    for f, (rdb, rdeg) in EXACT[deck].items():
        d, p = pac(deck, f, points, env_extra=env_extra, tag=tag)
        e[f] = (d - rdb, dphase(p, rdeg))
    return e


def check_deck(deck):
    name = f"{deck}: 200 points, 20/50/100 kHz within {BOUND_DB} dB / {BOUND_DEG} deg of the exact reference"
    t0 = time.time()
    try:
        e = errors(deck)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    wdb, wdeg = max(abs(v[0]) for v in e.values()), max(abs(v[1]) for v in e.values())
    record(name, "pass" if wdb <= BOUND_DB and wdeg <= BOUND_DEG else "FAIL",
           "; ".join(f"{f // 1000}k {v[0]:+.4f} dB {v[1]:+.4f} deg" for f, v in e.items()) + f"; wall {time.time() - t0:.1f}s")


def check_points():
    name = f"all six decks at 1000 points within {BOUND_DB} dB / {BOUND_DEG} deg of the exact reference"
    t0 = time.time()
    wdb = wdeg = 0.0
    where = ""
    try:
        for deck in EXACT:
            for f, v in errors(deck, points=1000).items():
                if abs(v[0]) > wdb:
                    wdb = abs(v[0])
                if abs(v[1]) > wdeg:
                    wdeg, where = abs(v[1]), f"{deck}@{f // 1000}k"
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    record(name, "pass" if wdb <= BOUND_DB and wdeg <= BOUND_DEG else "FAIL",
           f"worst {wdb:.4f} dB, {wdeg:.4f} deg ({where}); wall {time.time() - t0:.1f}s")


def check_red():
    name = "red proof: PSS_SALT_NO728=1 (s144): bst_q_csw3n@100k > 0.10 deg and bb_cl_csw@20k > 0.015 dB off the exact reference"
    try:
        d1, p1 = pac("bst_q_csw3n", 100000, env_extra={"PSS_SALT_NO728": "1"}, tag="_no728")
        d2, p2 = pac("bb_cl_csw", 20000, env_extra={"PSS_SALT_NO728": "1"}, tag="_no728")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    e1 = dphase(p1, EXACT["bst_q_csw3n"][100000][1])
    e2 = d2 - EXACT["bb_cl_csw"][20000][0]
    record(name, "pass" if abs(e1) > 0.10 and abs(e2) > 0.015 else "FAIL",
           f"bst_q_csw3n@100k {e1:+.4f} deg, bb_cl_csw@20k {e2:+.4f} dB (s144: -0.156 deg, +0.021 dB)")


def check_evcap():
    name = "event buffer full is REFUSED, not dropped: PSS_SALT_EVCAP=6, bst_q_csw3n@100k, 1000 points"
    try:
        pac("bst_q_csw3n", 100000, 1000, env_extra={"PSS_SALT_EVCAP": "6"}, tag="_evcap6")
    except CouldNotRun as ex:
        msg = str(ex)
        ok = msg.startswith("refused:") and "the rest would be lost" in msg
        record(name, "pass" if ok else "FAIL", msg[:200])
        return
    record(name, "FAIL", "ran to the end: no refusal with a buffer of 6")


def tool_python():
    py = os.environ.get("PAC_EXACT_PYTHON")
    if py:
        return py
    try:
        if subprocess.run([sys.executable, "-c", "import scipy, numpy"], capture_output=True, timeout=60).returncode == 0:
            return sys.executable
    except (OSError, subprocess.TimeoutExpired):
        pass
    return None


def check_tool():
    name = "tools/pac_exact.py (rtol 1e-9, 2 / 1 mV) reproduces two pinned values within 0.002 dB / 0.01 deg"
    py = tool_python()
    if not py:
        record(name, "COULD-NOT-RUN", "no python with numpy + scipy (set PAC_EXACT_PYTHON)")
        return
    t0 = time.time()
    parts, ok = [], True
    for deck, f in (("bst_cl_csw100p", 100000), ("bb_q_csw", 50000)):
        try:
            p = subprocess.run([py, "-B", os.path.join(ROOT, "tools", "pac_exact.py"), deck, str(f), "2e-3,1e-3", "1e-9"],
                               capture_output=True, text=True, timeout=1800)
            line = [l for l in p.stdout.split("\n") if "Richardson" in l][-1]
            v = line.split(":")[1].split()
            d, ph = float(v[0]), float(v[2])
        except (OSError, subprocess.TimeoutExpired, IndexError, ValueError) as ex:
            record(name, "COULD-NOT-RUN", f"{deck}@{f}: {ex}")
            return
        rdb, rdeg = EXACT[deck][f]
        ok = ok and abs(d - rdb) <= 0.002 and abs(dphase(ph, rdeg)) <= 0.01
        parts.append(f"{deck}@{f // 1000}k {d:.5f} dB {ph:.4f} deg (pinned {rdb:.5f} / {rdeg:.4f})")
    record(name, "pass" if ok else "FAIL", "; ".join(parts) + f"; wall {time.time() - t0:.0f}s")


def main():
    t0 = time.time()
    print(f"binary: {NGSPICE}\noutputs: {OUT}", flush=True)
    for deck in EXACT:
        check_deck(deck)
    check_points()
    check_red()
    check_evcap()
    if "--no-tool" not in sys.argv:
        check_tool()
    npass = sum(1 for r in results if r[1] == "pass")
    nfail = sum(1 for r in results if r[1] == "FAIL")
    ncnr = sum(1 for r in results if r[1] == "COULD-NOT-RUN")
    print(f"\n{npass}/{len(results)} pass, {nfail} FAIL, {ncnr} COULD-NOT-RUN; wall {time.time() - t0:.0f}s", flush=True)
    sys.exit(1 if nfail else (2 if ncnr else 0))


if __name__ == "__main__":
    main()
