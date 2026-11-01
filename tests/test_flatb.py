#!/usr/bin/env python3
"""tests/test_flatb.py -- flat B-source evaluator (card #532, round 15;
src/spicelib/parser/ifeval.c: FBcompile/FBrun).

IFeval compiles the value tree and all derivative trees of a B-source (one DAG)
once into a linear instruction list in the evaluation order of the memoized
tree walk, for ternary-free trees; NGSPICE_FLATB=0 keeps the tree walk.  This
test runs each deck with the evaluator on (default) and off and requires the
numerical output (wrdata, 15 digits) and the ngspice stderr/stdout error text
to be BIT-identical:

  d1  smooth switch / comparator (tanh, + - * /) in .tran and .ac
  d2  functions: pow, pwr, min, max, abs, sqrt, exp, log, atan, u, uramp, -x, ^
      with several controlling variables, in .dc and .tran, plus time/temper/hertz
  d3  ternary (a ? b : c) trees: must take the tree-walk fallback (identical by
      construction; proves the fallback is exercised)
  d4  error path: log of a value that goes <= 0 / division by zero during the
      sweep (HUGE) -- same messages, same rows written before the failure
  red-proof: NGSPICE_FLATB=2 (test-only) scales every flat-evaluated B-source
      value by 1+1e-9; the compare must then FAIL for d1/d2 (the flat program
      is really what evaluated them) and still pass for d3 (ternary: tree walk).

  NGSPICE=<prefix>/bin/ngspice python3 tests/test_flatb.py
  exit 0 pass, 1 fail, 2 could-not-run
"""
import os
import subprocess
import sys

NGSPICE = os.environ.get("NGSPICE") or os.path.expanduser(
    "~/Work/apps/ngspice_augmentation/aaba50ba3f49-r15-4/bin/ngspice")
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/r15/test_flatb")

D1 = """flatb d1
Vin in 0 PULSE(0 1 0 1n 1n 5u 10u) AC 1
Vr ramp 0 PULSE(0 1 0 9u 1u 1n 10u)
Bcmp pwm 0 V = 0.5*(1+tanh(30*(V(in)-V(ramp))))
Bhs in sw I = (V(in)-V(sw)) * (1e-6+(10-1e-6)*0.5*(1+tanh(20*(V(pwm)-0.5))))
Bls sw 0 I = V(sw) * (1e-6+(10-1e-6)*0.5*(1+tanh(20*(0.5-V(pwm)))))
R1 sw o 10
C1 o 0 1n
Rl o 0 1k
.control
set numdgt=15
set wr_singlescale
set wr_vecnames
tran 20n 20u
wrdata t.txt v(o) v(sw) v(pwm)
ac dec 5 1k 10meg
wrdata a.txt v(o)
quit
.endc
.end
"""
D2 = """flatb d2
Vs a 0 DC 0.5
Vb b 0 DC 1.5
Vt t 0 PULSE(0 2 0 1u 1u 3u 10u)
B1 n1 0 V = pow(abs(V(a))+0.1, 1.7) + pwr(V(b), 0.5) + min(V(a), V(b)) * max(V(a), 0.2)
B2 n2 0 V = sqrt(V(a)*V(a)+1) * exp(-V(b)/3) + log(V(b)+2) - atan(V(a)-V(b)) + u(V(a)-0.3) + uramp(V(b)-1)
B3 n3 0 V = -V(n1)/(V(n2)+3) + V(a)^2 + time*1e5 + temper/100 + hertz*0 + V(t)
R1 n1 0 1k
R2 n2 0 1k
R3 n3 0 1k
.control
set numdgt=15
set wr_singlescale
set wr_vecnames
dc Vs 0 2 0.01
wrdata d.txt v(n1) v(n2) v(n3)
tran 50n 20u
wrdata t.txt v(n1) v(n2) v(n3)
quit
.endc
.end
"""
D3 = """flatb d3
Vs a 0 PULSE(-1 1 0 1u 1u 4u 10u)
B1 n1 0 V = V(a) > 0 ? V(a)*2 : -V(a)
B2 n2 0 V = V(a) > 0.5 ? sqrt(V(a)) : 0.1*V(a)
R1 n1 0 1k
R2 n2 0 1k
.control
set numdgt=15
set wr_singlescale
set wr_vecnames
tran 20n 20u
wrdata t.txt v(n1) v(n2)
quit
.endc
.end
"""
D4 = """flatb d4
Vs a 0 DC 1
B1 n1 0 V = log(V(a)) + 1/V(a)
R1 n1 0 1k
.control
set numdgt=15
set wr_singlescale
set wr_vecnames
dc Vs -1 1 0.25
wrdata d.txt v(n1)
quit
.endc
.end
"""
DECKS = {"d1": D1, "d2": D2, "d3": D3, "d4": D4}


def run(tag, name, flatb):
    d = os.path.join(DATA, f"{name}_{tag}")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "deck.cir"), "w") as f:
        f.write(DECKS[name])
    for fn in os.listdir(d):
        if fn.endswith(".txt"):
            os.remove(os.path.join(d, fn))
    env = dict(os.environ)
    env["NGSPICE_FLATB"] = str(flatb)
    p = subprocess.run([NGSPICE, "-b", "deck.cir"], cwd=d, env=env, capture_output=True, text=True, timeout=300)
    files = {fn: open(os.path.join(d, fn)).read() for fn in sorted(os.listdir(d)) if fn.endswith(".txt")}
    # the transcript minus progress/timing lines
    keep = [ln for ln in (p.stdout + "\n" + p.stderr).splitlines()
            if "out of range" in ln or "rror" in ln or "singular" in ln.lower()]
    return p.returncode, files, keep


def main():
    if not os.access(NGSPICE, os.X_OK):
        print(f"COULD-NOT-RUN: ngspice not found at {NGSPICE}")
        return 2
    fail = False
    for name in DECKS:
        try:
            r0, f0, k0 = run("off", name, 0)
            r1, f1, k1 = run("on", name, 1)
        except Exception as e:  # noqa: BLE001
            print(f"COULD-NOT-RUN: {e}")
            return 2
        ok = (r0 == r1) and f0 == f1 and k0 == k1
        fail |= not ok
        rows = {k: v.count("\n") for k, v in f1.items()}
        print(f"[{'pass' if ok else 'FAIL'}] {name}: flat on == off, rc {r0}/{r1}, files {rows}, "
              f"error lines {len(k0)}/{len(k1)}")
        if name in ("d1", "d2", "d3") and not f1:
            print(f"[FAIL] {name}: no output written")
            fail = True
        if name == "d4" and not k1:
            print("[FAIL] d4: expected an error message in the transcript (error path not exercised)")
            fail = True
    for name, must_differ in (("d1", True), ("d2", True), ("d3", False)):
        r0, f0, k0 = run("off", name, 0)
        r2, f2, k2 = run("pert", name, 2)
        differs = f0 != f2
        ok = differs == must_differ
        fail |= not ok
        print(f"[{'pass' if ok else 'FAIL'}] red-proof {name}: FLATB=2 output {'differs' if differs else 'identical'} "
              f"(expected {'differs' if must_differ else 'identical: tree-walk fallback'})")
    print("OVERALL FAIL" if fail else "OVERALL PASS")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
