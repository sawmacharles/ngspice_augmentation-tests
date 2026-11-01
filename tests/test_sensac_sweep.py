#!/usr/bin/env python3
"""Card #672: .sens ac over a sweep of more than one point.

RC divider v1 -> r1 -> out, c1 and r2 to ground; .sens v(out) ac <sweep>, every vector read with wrdata.
Checks (binary from $NGSPICE, default s140; the previous build $NGSPICE_PREV, default s139):
  a. dec 2 1k 1meg (7 points): dV/dR1, dV/dR2, dV/dC1 at every point against the closed form of V = Z/(R1 + Z),
     Z = R2 || 1/(j w C), 1e-5 relative (.sens perturbs by a relative 1e-6, its own accuracy);
  b. lin 4 10k 40k: the frequencies are 10k, 20k, 30k, 40k (like .ac lin: N points from fstart to fstop), and the
     derivatives match at every point;
  c. a single point (lin 1 100k 100k) gives the same numbers as on the previous build (unchanged);
  d. the comparator can fail: on the previous build the dec sweep gives 0 for the resistors after the first point and
     the lin sweep's frequencies are wrong (the defects).

    NGSPICE=<build>/bin/ngspice python3 tests/test_sensac_sweep.py     exit 0 pass, 1 any FAIL, 2 could-not-run
Outputs go to ~/Work/data/ngspice_augmentation/672, never the repo.
"""
import math
import os
import subprocess
import sys

APPS = os.path.expanduser("~/Work/apps/ngspice_augmentation")
NGSPICE = os.environ.get("NGSPICE", os.path.join(APPS, "aaba50ba3f49-s152/bin/ngspice"))
PREV = os.environ.get("NGSPICE_PREV", os.path.join(APPS, "aaba50ba3f49-s139/bin/ngspice"))
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/672")
R1, R2, C = 1e3, 3e3, 1e-9
RESULTS = []


def check(name, ok, detail):
    RESULTS.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


def exact(f):
    w = 2 * math.pi * f
    z = 1 / (1 / R2 + 1j * w * C)
    v = z / (R1 + z)
    dv_dz = R1 / (R1 + z) ** 2
    return {"r1": -v / (R1 + z), "r2": dv_dz * z * z / (R2 * R2), "c1": dv_dz * (-1j * w * z * z)}


def sens(binary, name, sweep):
    """(frequencies, {device: [complex]})"""
    os.makedirs(DATA, exist_ok=True)
    out = os.path.join(DATA, name + ".txt")
    deck = os.path.join(DATA, name + ".cir")
    open(deck, "w").write(f"rc divider\nv1 in 0 dc 0 ac 1\nr1 in out 1k\nc1 out 0 1n\nr2 out 0 3k\n"
                          f".sens v(out) ac {sweep}\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n"
                          f"wrdata {out} r1 c1 r2\n.endc\n.end\n")
    if os.path.exists(out):
        os.remove(out)
    subprocess.run([binary, "-b", deck], capture_output=True, text=True, timeout=300)
    rows = [[float(x) for x in ln.split()] for ln in open(out).read().splitlines()[1:] if ln.strip()]
    return [r[0] for r in rows], {d: [complex(r[1 + 2 * i], r[2 + 2 * i]) for r in rows] for i, d in enumerate(("r1", "c1", "r2"))}


def worst(freqs, got):
    w = 0.0
    for k, f in enumerate(freqs):
        ref = exact(f)
        for d in ("r1", "c1", "r2"):
            w = max(w, abs(got[d][k] - ref[d]) / abs(ref[d]))
    return w


def main():
    for b in (NGSPICE, PREV):
        if not os.path.exists(b):
            print(f"COULD-NOT-RUN: {b} missing")
            return 2
    f, s = sens(NGSPICE, "dec", "dec 2 1k 1meg")
    w = worst(f, s)
    check("a dec sweep: every point vs the closed form", len(f) == 7 and w < 1e-5,
          f"{len(f)} points, max relative error {w:.2e}")
    f, s = sens(NGSPICE, "lin", "lin 4 10k 40k")
    w = worst(f, s)
    check("b lin sweep: N points fstart..fstop, every point vs the closed form",
          [round(x) for x in f] == [10000, 20000, 30000, 40000] and w < 1e-5,
          f"frequencies {[round(x) for x in f]}, max relative error {w:.2e}")
    f1, s1 = sens(NGSPICE, "one", "lin 1 100k 100k")
    f0, s0 = sens(PREV, "one_prev", "lin 1 100k 100k")
    check("c single point unchanged", f1 == f0 == [1e5] and s1 == s0, f"{s1['r1'][0]:.6e} (previous {s0['r1'][0]:.6e})")
    fp, sp = sens(PREV, "dec_prev", "dec 2 1k 1meg")
    fl, _ = sens(PREV, "lin_prev", "lin 4 10k 40k")
    check("d previous build shows the defects (the comparator can fail)",
          all(x == 0 for x in sp["r1"][1:]) and [round(x) for x in fl] != [10000, 20000, 30000, 40000],
          f"r1 after the first point {sorted({abs(x) for x in sp['r1'][1:]})}; lin frequencies {[f'{x:.3g}' for x in fl]}")
    print(f"{RESULTS.count(True)}/{len(RESULTS)} pass")
    return 1 if RESULTS.count(False) else 0


if __name__ == "__main__":
    sys.exit(main())
