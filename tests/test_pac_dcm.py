#!/usr/bin/env python3
"""Validation for card #711: td `.pac` of a buck in discontinuous conduction at the default 200 points per period.

Cause (measured): when the ideal switch turns on in discontinuous conduction the switch node floats (switch off 1 Meg, diode off) and the
pre-event node equation of the inductor-routed jump (pss_salt_ind_jump()) divides by a conductance of a few uS.  It was fed the inductor
current of the event step, which already holds the response to the flip over the event step (V dt / L); the pre-event node voltage came
out as -3.6 V at 200 points, +3.4 V at 800, +5.15 V at 3200 (true 5.741 V), so the jump, and with it the result, was first order in the
step: 0.40 dB / 0.42 deg at 200 points (s140).  The jump now uses the inductor current of the previous accepted point.

    NGSPICE=<prefix>/bin/ngspice python3 tests/test_pac_dcm.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
Checks (N/M pass is printed at the end):
  1  DCM buck, 60 ohm load (tests/fixtures/dcm711/dcm_buck.cir), 200 points, 20/50/100 kHz vs the perturbed-PSS reference: 0.07 dB / 0.35 deg
  2  the same deck at 2000 points: within 0.01 dB / 0.05 deg of the 200-point result (converged in the points)
  3  DCM buck, 150 ohm load (dcm_buck_150.cir), 200 points, 20/50 kHz vs its own perturbed-PSS reference
  4  continuous conduction, four decks (buck, buck at duty 0.7, boost, synchronous buck with dead time), 20 and 100 kHz: within 0.01 dB / 0.05 deg
     of the s140 values (no change); the synchronous buck is compared with its card #727 value instead (s140: -7.7642 / -34.8826 dB, the
     +0.04 dB dead-time offset that #727 removed; now -7.8060 / -34.9243 dB, the perturbed-PSS reference of tests/test_pxf_pstb_switch.py)
  5  red proof: PSS_SALT_NO711=1 (the event step's inductor current again) on this binary: > 0.3 dB off at 20 kHz
  6  red proof: the baseline binary (NGSPICE_BASE, default s140): > 0.3 dB off at 20 kHz

References (dB, deg of v(vout)/v(Vctrl)):
  dcm_buck.cir: tools/pac_perturbed.py sweep --inject source:Vctrl --outputs vout --fs 500e3 --backend pss (800 -> 1600 points/cycle, s119 binary),
    card #581 (~/Work/data/ngspice_augmentation/581), the same deck as the pinned 60 ohm values of tests/test_pac_charge.py (od_cl_dcm).
  dcm_buck_150.cir: the same tool and settings, s140 binary, run for this card (~/Work/data/ngspice_augmentation/711/ref150.log).
  continuous-conduction values: s140, 200 points, reltol 1e-4, run for this card (they equal the new binary's to 4 decimals).
Decks and outputs live under ~/Work/data/ngspice_augmentation/711/test (never the repository).
"""
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures", "dcm711")
NGSPICE = os.environ.get("NGSPICE", "ngspice")
NGSPICE_BASE = os.environ.get(
    "NGSPICE_BASE",
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s140/bin/ngspice"))
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/711/test")

BOUND_DB, BOUND_DEG = 0.07, 0.35

REF_DCM60 = {20000: (-31.1878, -84.6675), 50000: (-38.8166, -77.5846), 100000: (-43.84, -68.0311)}
REF_DCM150 = {}   # filled in from the reference run, see below
# s140, 200 points: deck -> {f: (dB, deg)}
S140_CCM = {
    "ccm_buck.cir": {20000: (-7.4042, -167.379), 100000: (-34.5050, -146.824)},
    "ccm_buck_d70.cir": {20000: (-7.3279, -169.172), 100000: (-34.4432, -147.177)},
    "ccm_boost.cir": {20000: (-1.6781, 168.785), 100000: (-25.1813, 130.431)},
    "ccm_sync.cir": {20000: (-7.8060, -169.820), 100000: (-34.9243, -147.295)},   # card #727 (s140: -7.7642, -34.8826 dB; the +0.04 dB offset)
}

results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"{state:14s} {name}\n{'':14s}   {detail}", flush=True)


def pac(binary, deck, f, points, tag, env_extra=None):
    """One td .pac point of a fixture deck; returns (dB, deg) of v(vout)."""
    os.makedirs(DATA, exist_ok=True)
    body = open(os.path.join(FIX, deck)).read().replace("\n.end\n", "\n")
    out = os.path.join(DATA, f"{tag}.txt")
    if os.path.exists(out):
        os.remove(out)
    d = (body + ".options reltol=1e-4\n"
         + '.pac 500e3 100e-6 vout %d 10 150 1e-6 lin 1 %g %g method="shoot" solver="td"\n' % (points, f, f)
         + ".control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n"
         + f"wrdata {out} v(vout)\nquit\n.endc\n.end\n")
    path = os.path.join(DATA, f"{tag}.cir")
    open(path, "w").write(d)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "2"
    if env_extra:
        env.update(env_extra)
    try:
        p = subprocess.run([binary, "-b", path], capture_output=True, text=True, timeout=600, env=env)
    except (OSError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{binary}: {ex}")
    if "REFUS" in p.stdout + p.stderr:
        raise CouldNotRun("refused: " + [l for l in (p.stdout + p.stderr).split("\n") if "REFUS" in l][0][:200])
    if not os.path.exists(out):
        raise CouldNotRun(f"no output for {tag} (rc {p.returncode})")
    t = open(out).read().split("\n")[1].split()
    z = complex(float(t[1]), float(t[2]))
    return 20 * math.log10(abs(z)), math.degrees(math.atan2(z.imag, z.real))


def dphase(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def check_vs_ref(name, deck, ref, points=200):
    t0 = time.time()
    parts, wdb, wdeg = [], 0.0, 0.0
    try:
        for f, (rdb, rdeg) in ref.items():
            d, p = pac(NGSPICE, deck, f, points, f"{deck[:-4]}_{f}_{points}")
            wdb, wdeg = max(wdb, abs(d - rdb)), max(wdeg, abs(dphase(p, rdeg)))
            parts.append(f"{f / 1e3:g}k {d:.4f}/{p:.3f} vs {rdb:.4f}/{rdeg:.3f}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    ok = wdb <= BOUND_DB and wdeg <= BOUND_DEG
    record(name, "pass" if ok else "FAIL",
           "; ".join(parts) + f"; worst {wdb:.4f} dB / {wdeg:.3f} deg (bound {BOUND_DB} / {BOUND_DEG}), wall {time.time() - t0:.0f}s")


def check_converged():
    name = "DCM buck 60 ohm: 200 points within 0.01 dB / 0.05 deg of 2000 points"
    try:
        a = pac(NGSPICE, "dcm_buck.cir", 20000, 200, "conv_200")
        b = pac(NGSPICE, "dcm_buck.cir", 20000, 2000, "conv_2000")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    dd, dp = abs(a[0] - b[0]), abs(dphase(a[1], b[1]))
    record(name, "pass" if dd <= 0.01 and dp <= 0.05 else "FAIL",
           f"200: {a[0]:.4f}/{a[1]:.3f}, 2000: {b[0]:.4f}/{b[1]:.3f}; diff {dd:.4f} dB / {dp:.3f} deg")


def check_ccm():
    name = "continuous conduction (4 decks x 20/100 kHz): unchanged vs s140 (ccm_sync: vs its card #727 value) within 0.01 dB / 0.05 deg"
    worst_db = worst_deg = 0.0
    parts = []
    try:
        for deck, ref in S140_CCM.items():
            for f, (rdb, rdeg) in ref.items():
                d, p = pac(NGSPICE, deck, f, 200, f"{deck[:-4]}_{f}")
                worst_db, worst_deg = max(worst_db, abs(d - rdb)), max(worst_deg, abs(dphase(p, rdeg)))
                parts.append(f"{deck[:-4]}@{f // 1000}k {d - rdb:+.4f}/{dphase(p, rdeg):+.3f}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    record(name, "pass" if worst_db <= 0.01 and worst_deg <= 0.05 else "FAIL",
           "; ".join(parts) + f"; worst {worst_db:.4f} dB / {worst_deg:.3f} deg")


def check_red(name, binary, env_extra=None):
    try:
        d, p = pac(binary, "dcm_buck.cir", 20000, 200, "red_" + str(abs(hash(name)) % 10000), env_extra)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    rdb = REF_DCM60[20000][0]
    record(name, "pass" if abs(d - rdb) > 0.3 else "FAIL",
           f"20 kHz {d:.4f} dB vs reference {rdb:.4f}: error {abs(d - rdb):.3f} dB (must exceed 0.3)")


def main():
    r150 = os.path.join(HERE, "fixtures", "dcm711", "ref150.txt")
    if os.path.exists(r150):   # "f dB deg" per line, pinned from the reference run
        for l in open(r150):
            t = l.split()
            if len(t) == 3:
                REF_DCM150[int(float(t[0]))] = (float(t[1]), float(t[2]))
    check_vs_ref("DCM buck 60 ohm, 200 points, td .pac vs perturbed PSS", "dcm_buck.cir", REF_DCM60)
    check_converged()
    if REF_DCM150:
        check_vs_ref("DCM buck 150 ohm, 200 points, td .pac vs perturbed PSS", "dcm_buck_150.cir", REF_DCM150)
    else:
        record("DCM buck 150 ohm vs perturbed PSS", "COULD-NOT-RUN", "no pinned reference (tests/fixtures/dcm711/ref150.txt)")
    check_ccm()
    check_red("[RED PROOF] PSS_SALT_NO711=1 on this binary: DCM 60 ohm > 0.3 dB off", NGSPICE, {"PSS_SALT_NO711": "1"})
    if os.path.exists(NGSPICE_BASE):
        check_red("[RED PROOF] baseline binary (s140): DCM 60 ohm > 0.3 dB off", NGSPICE_BASE)
    else:
        record("[RED PROOF] baseline binary", "COULD-NOT-RUN", f"{NGSPICE_BASE} not found")
    n = len(results)
    npass = sum(1 for r in results if r[1] == "pass")
    print(f"\n{npass}/{n} pass")
    if any(r[1] == "FAIL" for r in results):
        return 1
    if any(r[1] == "COULD-NOT-RUN" for r in results):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
