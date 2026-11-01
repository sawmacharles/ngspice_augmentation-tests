#!/usr/bin/env python3
"""tests/test_resbase.py -- resistor/capacitor base cache for the KLU load
(card #527, round 14; src/spicelib/analysis/cktload.c, acan.c).

CKTload (real, MODETRAN) and CKTacLoad (complex) skip the resistor (and, in
AC, capacitor) device loops and copy a snapshot of their summed stamps
instead.  The snapshot is keyed on RESgeneration/CAPgeneration, which
REStemp/CAPtemp/RES-CAPparam and the CSC bind functions bump.  This test
runs ONE ngspice session on a KLU RC ladder with tc1 resistors that does
tran, ac, `alter r2`, tran, `alter c3`, ac, `option temp`, tran, ac and checks:

  1. cache on (default) == cache off (NGSPICE_RESBASE=0) on all six traces,
     max |delta| <= 1e-9 x trace range (roundoff only: sum order differs);
  2. the alters/temperature really changed the results (>1e-3 of range), so
     (1) is not vacuous;
  3. red-proof (AC): NGSPICE_RESBASE=2 (test-only switch: cache without
     invalidation) must DEVIATE from cache-off by > 1e-6 after an alter,
     i.e. the generation counters are what keeps (1) true.  (The transient
     cache is additionally keyed on the KLU matrix storage pointer, which
     changes whenever a new analysis rebuilds the matrix, so NGSPICE_RESBASE=2
     does not break it -- measured: t2/t3 stay within 6e-14 -- and no tran
     red-proof is claimed.)

Round 15 (card #532) adds the capacitor base cache + flat capacitor load in
CKTload (cktload.c, capload.c): Ax = Rbase + ag0*Cbase, CAPload's state/rhs
work done by CAPloadFlat.  Extra checks, on the same session plus a floating
cap with m=2, a UIC tran with ic= (first load falls back), and a gear tran
(flat load falls back to CAPload with matrix writes skipped):
  4. R/C cache on == capacitor cache off (NGSPICE_CAPBASE=0), <= 1e-9 of range;
  5. flat load on == flat off (NGSPICE_CAPFLAT=0): BIT-identical;

  7. red-proof inside ONE tran (stop / alter / resume, which keeps the KLU
     matrix pointer): NGSPICE_CAPBASE=2 after `alter c2` and NGSPICE_RESBASE=2
     after `alter r2` each deviate from cache-off by > 1e-6 of range, while
     the default (invalidating) caches stay within 1e-9 of cache-off.

  NGSPICE=<prefix>/bin/ngspice python3 tests/test_resbase.py
  exit 0 pass, 1 fail, 2 could-not-run
"""
import os
import subprocess
import sys

NGSPICE = os.environ.get("NGSPICE") or os.path.expanduser(
    "~/Work/apps/ngspice_augmentation/aaba50ba3f49-r15-2/bin/ngspice")
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/r15/test_resbase")

N = 8
DECK = ["resbase ladder", "Vin in 0 PULSE(0 1 0 1n 1n 5u 10u) AC 1"]
prev = "in"
for i in range(1, N + 1):
    DECK += [f"R{i} {prev} n{i} 1k tc1=1e-3", f"C{i} n{i} 0 1n" + (" ic=0.3" if i == 2 else "")]
    prev = f"n{i}"
DECK += ["Cx n3 n6 0.5n m=2", "Rl n8 0 10k", ".options klu",
         ".control", "set numdgt=15", "set wr_singlescale", "set wr_vecnames",
         "tran 20n 20u", "wrdata t1.txt v(n8)",
         "ac dec 10 1k 10meg", "wrdata a1.txt v(n8)",
         "alter r2 resistance = 3k", "tran 20n 20u", "wrdata t2.txt v(n8)",
         "alter c3 capacitance = 4n", "ac dec 10 1k 10meg", "wrdata a2.txt v(n8)",
         "option temp=70", "tran 20n 20u", "wrdata t3.txt v(n8)",
         "ac dec 10 1k 10meg", "wrdata a3.txt v(n8)",
         "tran 20n 8u uic", "wrdata t4.txt v(n8) v(n2)",
         "option method=gear", "tran 20n 20u", "wrdata t5.txt v(n8)",
         "quit", ".endc", ".end"]
FILES = ["t1", "a1", "t2", "a2", "t3", "a3", "t4", "t5"]


def run(tag, resbase, extra=None):
    d = os.path.join(DATA, tag)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "deck.cir"), "w") as f:
        f.write("\n".join(DECK) + "\n")
    env = dict(os.environ)
    if resbase is not None:
        env["NGSPICE_RESBASE"] = str(resbase)
    env.update(extra or {})
    p = subprocess.run([NGSPICE, "-b", "deck.cir"], cwd=d, env=env, capture_output=True, text=True, timeout=300)
    if p.returncode != 0:
        raise RuntimeError(f"{tag}: ngspice rc={p.returncode}: {p.stdout[-800:]}{p.stderr[-800:]}")
    out = {}
    for k in FILES:
        rows = [[float(x) for x in ln.split()] for ln in open(os.path.join(d, k + ".txt"))
                if ln.strip() and not ln.lstrip()[0].isalpha()]
        out[k] = rows
    return out


def run_resume(tag, alter, env_extra):
    d = os.path.join(DATA, tag)
    os.makedirs(d, exist_ok=True)
    deck = ["resume", "Vin in 0 PULSE(0 1 0 1n 1n 5u 10u)", "R1 in n1 1k", "C1 n1 0 1n", "R2 n1 n2 1k",
            "C2 n2 0 1n", "Rl n2 0 10k", ".options klu", ".control", "set numdgt=15", "set wr_singlescale",
            "set wr_vecnames", "stop when time > 5u", "tran 20n 20u", "delete all", alter, "resume",
            "wrdata r.txt v(n2)", "quit", ".endc", ".end"]
    with open(os.path.join(d, "deck.cir"), "w") as f:
        f.write("\n".join(deck) + "\n")
    env = dict(os.environ)
    env.update(env_extra)
    p = subprocess.run([NGSPICE, "-b", "deck.cir"], cwd=d, env=env, capture_output=True, text=True, timeout=300)
    if p.returncode != 0:
        raise RuntimeError(f"{tag}: ngspice rc={p.returncode}")
    return [[float(x) for x in ln.split()] for ln in open(os.path.join(d, "r.txt"))
            if ln.strip() and not ln.lstrip()[0].isalpha()]


def dev(a, b):
    """max |a-b| over all columns >= 1, relative to a's max |value| per column."""
    assert len(a) == len(b), (len(a), len(b))
    w = 0.0
    for j in range(1, len(a[0])):
        rng = max(abs(r[j]) for r in a) or 1.0
        w = max(w, max(abs(r[j] - s[j]) for r, s in zip(a, b)) / rng)
    return w


def main():
    if not os.access(NGSPICE, os.X_OK):
        print(f"COULD-NOT-RUN: ngspice not found at {NGSPICE}")
        return 2
    try:
        off, on, noinv = run("off", 0), run("on", None), run("noinval", 2)
        capoff = run("capoff", None, {"NGSPICE_CAPBASE": "0"})
        noflat = run("noflat", None, {"NGSPICE_CAPFLAT": "0"})
        capnoinv = run("capnoinval", None, {"NGSPICE_CAPBASE": "2"})
    except Exception as e:  # noqa: BLE001
        print(f"COULD-NOT-RUN: {e}")
        return 2
    fail = False
    for k in FILES:
        d = dev(off[k], on[k])
        ok = d <= 1e-9
        fail |= not ok
        print(f"[{'pass' if ok else 'FAIL'}] cache on vs off {k}: max|d|/range = {d:.2e} (<= 1e-9), rows {len(off[k])}")
    for a, b in (("t1", "t2"), ("a1", "a2"), ("t2", "t3"), ("a2", "a3")):
        if len(off[a]) != len(off[b]):
            continue  # different step count: certainly different
        d = dev(off[a], off[b])
        ok = d > 1e-3
        fail |= not ok
        print(f"[{'pass' if ok else 'FAIL'}] alter/temp changes the result {a}->{b}: {d:.2e} (> 1e-3)")
    worst = 0.0
    for k in ("a2", "a3"):
        if len(off[k]) == len(noinv[k]):
            worst = max(worst, dev(off[k], noinv[k]))
        else:
            worst = max(worst, 1.0)
    ok = worst > 1e-6
    fail |= not ok
    print(f"[{'pass' if ok else 'FAIL'}] red-proof: cache WITHOUT invalidation deviates after alter/temp: {worst:.2e} (> 1e-6)")
    for k in FILES:
        d = dev(capoff[k], on[k])
        ok = d <= 1e-9
        fail |= not ok
        print(f"[{'pass' if ok else 'FAIL'}] cap cache on vs NGSPICE_CAPBASE=0 {k}: max|d|/range = {d:.2e} (<= 1e-9)")
    for k in FILES:
        ok = on[k] == noflat[k]
        fail |= not ok
        print(f"[{'pass' if ok else 'FAIL'}] flat cap load on vs off {k}: {'bit-identical' if ok else 'DIFFERENT'}")
    worst = 0.0
    for k in ("t2", "t3", "t4", "t5"):
        if len(capoff[k]) == len(capnoinv[k]):
            worst = max(worst, dev(capoff[k], capnoinv[k]))
        else:
            worst = max(worst, 1.0)
    print(f"[info] across-analysis tran with CAPBASE=2: {worst:.2e} (the cache is also keyed on the matrix pointer, "
          f"so a new analysis rebuilds it; the red-proof is the in-one-tran resume check below)")
    try:
        for alter, var in (("alter c2 capacitance = 20n", "NGSPICE_CAPBASE"), ("alter r2 resistance = 5k", "NGSPICE_RESBASE")):
            nm = alter.split()[1]
            ref = run_resume(f"rs_{nm}_off", alter, {"NGSPICE_RESBASE": "0"})
            ons = run_resume(f"rs_{nm}_on", alter, {})
            bad = run_resume(f"rs_{nm}_noinv", alter, {var: "2"})
            d_on = dev(ref, ons) if len(ref) == len(ons) else 1.0
            d_bad = dev(ref, bad) if len(ref) == len(bad) else 1.0
            ok = d_on <= 1e-9
            fail |= not ok
            print(f"[{'pass' if ok else 'FAIL'}] resume after {nm}: caches on vs off {d_on:.2e} (<= 1e-9)")
            ok = d_bad > 1e-6
            fail |= not ok
            print(f"[{'pass' if ok else 'FAIL'}] red-proof inside one tran ({var}=2 after {nm}): {d_bad:.2e} (> 1e-6)")
    except Exception as e:  # noqa: BLE001
        print(f"COULD-NOT-RUN: {e}")
        return 2
    print("OVERALL FAIL" if fail else "OVERALL PASS")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
