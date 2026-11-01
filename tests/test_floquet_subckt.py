#!/usr/bin/env python3
"""Card #682: the shooting PSS Floquet verdict does not depend on how the circuit is netlisted.

The R4d dead-time peak-current-mode buck, once flat (fixtures/floquet682/r4d_flat.cir) and once with every element
wrapped in a .subckt as schematics v4 netlists it (r4d_subckt.cir, devices from tools/s4view/r4d_devices.lib).
On s140 the subckt deck reported spectral_radius=1.266873 UNSTABLE(period-doubling) while the flat one read 0.994
stable; a transient of either shows no period doubling. Cause: the dedicated Floquet re-evaluation differenced the
three logic-like states (latch and dead-time delay capacitors, |x| ~ 0.01 V, below their ripple) with the 1e-6
floor step, so their columns were period-map noise / h; the Newton columns already used 1 % of the ripple (#505).

Checks (fork binary from $NGSPICE):
  1. flat deck: CONVERGED, stable, spectral radius 0.99 .. 1.0 (the slow output-filter mode, 0.994 on s140).
  2. subckt deck: CONVERGED, stable, spectral radius 0.99 .. 1.0.
  3. the two radii agree within 1e-3 (same circuit; their orbits differ by solver tolerance only).
  4. red proof: subckt deck with PSS_SHOOT_NO_BIGH=1 (the old small-step rule) is not 'stable' or its radius
     is off by > 0.05 -- shows check 2 would catch the defect.

    NGSPICE=<build>/bin/ngspice python3 tests/test_floquet_subckt.py    exit 0 pass, 1 any FAIL, 2 could-not-run
"""
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures", "floquet682")
DEVLIB = os.path.join(os.path.dirname(HERE), "tools", "s4view", "r4d_devices.lib")
NGSPICE = os.environ.get("NGSPICE", "ngspice")
PAT = re.compile(r"PSS-SHOOT: CONVERGED .*spectral_radius=([-+0-9.eE]+) (\S+)")


def run(deck, extra_env=None):
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "2")
    if extra_env:
        env.update(extra_env)
    with tempfile.TemporaryDirectory(prefix="floq682_") as td:
        text = open(os.path.join(FIX, deck)).read().replace("@DEVLIB@", DEVLIB)
        path = os.path.join(td, deck)
        open(path, "w").write(text)
        p = subprocess.run([NGSPICE, "-b", path], capture_output=True, text=True, timeout=900, env=env, cwd=td)
    m = PAT.search(p.stdout + p.stderr)
    return (float(m.group(1)), m.group(2)) if m else (None, None)


def main():
    if not os.path.exists(DEVLIB):
        print("COULD-NOT-RUN: missing", DEVLIB)
        return 2
    results = []

    def check(name, ok, detail):
        results.append(ok)
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")

    rf, vf = run("r4d_flat.cir")
    rs, vs = run("r4d_subckt.cir")
    if rf is None or rs is None:
        print(f"COULD-NOT-RUN: no CONVERGED line (flat {rf}, subckt {rs})")
        return 2
    check("flat deck stable", vf == "stable" and 0.99 <= rf <= 1.0, f"spectral_radius={rf} {vf}")
    check("subckt deck stable", vs == "stable" and 0.99 <= rs <= 1.0, f"spectral_radius={rs} {vs}")
    check("same radius", abs(rf - rs) <= 1e-3, f"|{rf} - {rs}| = {abs(rf - rs):.2e} (tol 1e-3)")
    ro, vo = run("r4d_subckt.cir", {"PSS_SHOOT_NO_BIGH": "1"})
    check("red proof: old step rule shows the defect", ro is None or vo != "stable" or abs(ro - rf) > 0.05,
          f"PSS_SHOOT_NO_BIGH=1: spectral_radius={ro} {vo}")
    n = sum(results)
    print(f"\n{n}/{len(results)} pass")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
