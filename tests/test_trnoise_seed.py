#!/usr/bin/env python3
"""Card #668: transient noise (trnoise) follows the seed.

Checks (binary from $NGSPICE, default s139; the previous build $NGSPICE_PREV, default s138):
  a. with `.option seed=N` two runs give identical trnoise waveforms (white and 1/f), data blocks byte for byte;
  b. a different seed gives a different waveform;
  c. without a seed two runs still differ, as before (the behaviour without a seed is unchanged);
  d. the comparator can fail: on the previous build two seeded runs differ (the defect this card fixes);
  e. every other seeded draw is unchanged: agauss/gauss/unif parameters with `.option seed` give exactly the previous
     build's values (the pool refill must not move the uniform stream; handbook example 19 caught the first version).

    NGSPICE=<build>/bin/ngspice python3 tests/test_trnoise_seed.py     exit 0 pass, 1 any FAIL, 2 could-not-run
Outputs go to ~/Work/data/ngspice_augmentation/668, never the repo.
"""
import os
import re
import subprocess
import sys

APPS = os.path.expanduser("~/Work/apps/ngspice_augmentation")
NGSPICE = os.environ.get("NGSPICE", os.path.join(APPS, "aaba50ba3f49-s139/bin/ngspice"))
PREV = os.environ.get("NGSPICE_PREV", os.path.join(APPS, "aaba50ba3f49-s138/bin/ngspice"))
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/668")
RESULTS = []


def check(name, ok, detail):
    RESULTS.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


def data_block(binary, name, seed):
    os.makedirs(DATA, exist_ok=True)
    p = os.path.join(DATA, name + ".cir")
    seedline = f".option seed={seed}\n" if seed else ""
    open(p, "w").write(f"trnoise seed\n{seedline}v1 a 0 dc 0 trnoise(1m 1u 1 10u)\nr1 a 0 1k\n"
                       "v2 b 0 dc 0 trnoise(1m 1u 0 0)\nr2 b 0 1k\n.tran 1u 2m\n.end\n")
    raw = os.path.join(DATA, name + ".raw")
    if os.path.exists(raw):
        os.remove(raw)
    subprocess.run([binary, "-b", "-r", raw, p], capture_output=True, text=True, timeout=300)
    b = open(raw, "rb").read()
    return b[b.index(b"Binary:\n"):]


def param_values(binary, name):
    os.makedirs(DATA, exist_ok=True)
    p = os.path.join(DATA, name + ".cir")
    open(p, "w").write("seeded params\n.option seed=1\n.param a = agauss(0,10m,3)\n.param b = {a}\n.param c = {a}\n"
                       ".param u = unif(1, 0.1)\nR1 1 0 '1000*(1+b)'\nR2 2 0 '1000*(1+c)'\nR3 3 0 '1000*u'\n"
                       "V1 1 0 DC 1\nV2 2 0 DC 1\nV3 3 0 DC 1\n.control\nset numdgt=15\nop\n"
                       "print @r1[resistance] @r2[resistance] @r3[resistance]\n.endc\n.end\n")
    out = subprocess.run([binary, "-b", p], capture_output=True, text=True, timeout=300).stdout
    return re.findall(r"resistance\]\s*=\s*(\S+)", out)


def main():
    for b in (NGSPICE, PREV):
        if not os.path.exists(b):
            print(f"COULD-NOT-RUN: {b} missing")
            return 2
    a1, a2 = data_block(NGSPICE, "s5a", 5), data_block(NGSPICE, "s5b", 5)
    check("a same seed, same noise", a1 == a2 and len(a1) > 10000, f"identical: {a1 == a2} ({len(a1)} bytes)")
    a3 = data_block(NGSPICE, "s6", 6)
    check("b another seed, other noise", a3 != a1, "different")
    n1, n2 = data_block(NGSPICE, "n1", None), data_block(NGSPICE, "n2", None)
    check("c no seed: runs still differ (unchanged)", n1 != n2, "different")
    p1, p2 = data_block(PREV, "p5a", 5), data_block(PREV, "p5b", 5)
    check("d previous build: seeded runs differ (the defect; the comparator can fail)", p1 != p2, "different on the previous build")
    v_new, v_old = param_values(NGSPICE, "par_new"), param_values(PREV, "par_old")
    check("e other seeded draws unchanged (agauss, unif)", len(v_new) == 3 and v_new == v_old,
          f"this build {v_new}, previous {v_old}")
    print(f"{RESULTS.count(True)}/{len(RESULTS)} pass")
    return 1 if RESULTS.count(False) else 0


if __name__ == "__main__":
    sys.exit(main())
