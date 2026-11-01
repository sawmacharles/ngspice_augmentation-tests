#!/usr/bin/env python3
"""Validation for card #710: td `.pac` of the transistor-level (BSIM4) closed-loop buck R7 against a perturbed periodic steady state.

R7 (refs/r7_bsim_buck.cir, card #539) is a 5 V -> 1.46 V, 2 A, 500 kHz synchronous buck in which everything that conducts or decides is
a BSIM4 MOSFET (comparator, NOR dead-time logic, tapered drivers, OTA error amp).  Until #710 its `.pac` was checked only against FRA
curves (0.09-0.12 dB uncertainty).  The reference here is tools/pac_perturbed.py: it runs the real nonlinear circuit to two periodic steady
states (baseline and with a small sine on Vctrl) and reads the same-frequency component of the difference; no small-signal code.

    NGSPICE=<prefix>/bin/ngspice python3 tests/test_r7_perturbed.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
Checks (N/M pass is printed at the end):
  1-3  closed-loop control-to-output v(vout)/v(Vctrl) at 20 / 50 / 100 kHz (tests/fixtures/r7_pac_ctrl2out.cir, 400 points, reltol 1e-4)
       within 0.02 dB / 0.15 deg of the reference (measured on s144: -0.0017 / +0.045, +0.0035 / +0.016, +0.0051 / +0.020)
  4    red proof: PAC_NATIVE_HISTORY=1 (device-model charge history ignored, the card #542 defect) at 20 kHz: > 3 dB off (s144: -119.5 dB)

Reference (dB, deg), 2026-10-19, s144 binary, ~/Work/data/ngspice_augmentation/710/ref.log (61 min, one process):
  python3 tools/pac_perturbed.py sweep --deck r7_ctrl.cir --inject source:Vctrl --outputs vout --fs 500e3 --max-procs 1 \\
      --freqs 20000,50000,100000 --amplitude 0.005
  (r7_ctrl.cir = r7_pac_ctrl2out.cir up to its .options line), backend pss, 200 -> 400 points/cycle:
  amplitude halving +0.0012 / +0.0003 / +0.0005 dB, step residual 0.0001 dB / 0.002 deg, 0.0001 / 0.006, 0.0013 / 0.018.
Decks and outputs live under ~/Work/data/ngspice_augmentation/710/test (never the repository).
"""
import cmath
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE = os.path.join(HERE, "fixtures", "r7_pac_ctrl2out.cir")
NGSPICE = os.environ.get("NGSPICE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice")
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/710/test")

BOUND_DB, BOUND_DEG = 0.02, 0.15
REF = {20000: (3.6001, -7.4101), 50000: (-4.6896, -76.8939), 100000: (-11.8511, -86.4057)}

results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"[{state:13s}] {name}: {detail}", flush=True)


def run_pac(f, tag, env=None):
    """td .pac of the R7 fixture at one frequency; returns (dB, deg) of v(vout)."""
    os.makedirs(DATA, exist_ok=True)
    out = os.path.join(DATA, f"pac_{tag}.txt")
    deck = os.path.join(DATA, f"pac_{tag}.cir")
    with open(FIXTURE) as fh:
        txt = fh.read()
    txt = re.sub(r"(?m)^\.pac .*$", f'.pac 500e3 300e-6 vout 400 10 150 1e-6 lin 1 {f} {f} method="shoot" solver="td"', txt)
    txt = re.sub(r"(?m)^wrdata .*$", f"wrdata {out} v(vout)", txt)
    with open(deck, "w") as fh:
        fh.write(txt)
    if os.path.exists(out):
        os.remove(out)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    if env:
        e.update(env)
    t0 = time.time()
    try:
        p = subprocess.run([NGSPICE, "-b", deck], capture_output=True, text=True, cwd=DATA, env=e, timeout=900)
    except subprocess.TimeoutExpired:
        raise CouldNotRun(f"{tag}: timeout")
    w = time.time() - t0
    if not os.path.exists(out):
        raise CouldNotRun(f"{tag}: no output (rc={p.returncode}): {(p.stdout + p.stderr)[-200:]}")
    rows = [l.split() for l in open(out) if l.strip()]
    v = rows[1]
    H = complex(float(v[1]), float(v[2]))
    return 20 * math.log10(abs(H)), math.degrees(cmath.phase(H)), w


def ddeg(a, b):
    return (a - b + 180) % 360 - 180


def main():
    if not os.path.exists(NGSPICE):
        record("binary", "COULD-NOT-RUN", f"{NGSPICE} not found")
        return 2
    for f, (rdb, rdeg) in REF.items():
        name = f"R7 closed-loop control-to-output at {f // 1000} kHz, td .pac vs perturbed PSS ({BOUND_DB} dB / {BOUND_DEG} deg)"
        try:
            db, deg, w = run_pac(f, str(f))
        except CouldNotRun as ex:
            record(name, "COULD-NOT-RUN", str(ex))
            continue
        e_db, e_deg = db - rdb, ddeg(deg, rdeg)
        ok = abs(e_db) <= BOUND_DB and abs(e_deg) <= BOUND_DEG
        record(name, "pass" if ok else "FAIL",
               f"{db:.4f} dB / {deg:.3f} deg vs {rdb:.4f} / {rdeg:.3f}: {e_db:+.4f} dB / {e_deg:+.3f} deg; wall {w:.1f}s")
    name = "[RED PROOF] PAC_NATIVE_HISTORY=1 (device charge history ignored) at 20 kHz: > 3 dB off"
    try:
        db, deg, w = run_pac(20000, "red_native", {"PAC_NATIVE_HISTORY": "1"})
        e_db = db - REF[20000][0]
        record(name, "pass" if abs(e_db) > 3 else "FAIL", f"{db:.4f} dB / {deg:.3f} deg: {e_db:+.2f} dB; wall {w:.1f}s")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
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
