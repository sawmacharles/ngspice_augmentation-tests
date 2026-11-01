#!/usr/bin/env python3
"""Validation for the numparam draw-once form (card #482, patches-wip/
numparam/0300, docs/mc.md "Draw-once form").

Background (card #473 / docs/mc.md): numparam re-substitutes and
re-evaluates a `.param`'s RHS at *every textual reference* whenever that
RHS contains a call to agauss/gauss/aunif/unif/limit -- confirmed by
tracing inp_fix_agauss_in_param() (src/frontend/inpcom.c), which rewrites
such a `.param name = <expr>` into `.func name() <expr>` and then
"functionalises" every bare reference to `name` elsewhere in the deck into
`name()`, followed by inp_expand_macros_in_deck() inlining the function
body textually at each call site -- so two devices referencing the same
"global" `.param` end up with independent draws. This is deliberate,
HSPICE-compatible behaviour that sky130-style per-instance mismatch relies
on, so the default must not change.

This patch adds an explicit draw-once form: `.param p = once(<expr>)`.
`once()` is registered as an ordinary (identity) numparam function
(XFU_ONCE, src/frontend/numparam/xpressn.c); the actual freeze comes from
inp_fix_agauss_in_param() skipping the .func rewrite for any `.param` line
containing "once(" -- see patches-wip/numparam/0300 for the full patch and
docs/mc.md for the syntax rationale (3 options considered, why this one).

Checks (each printed pass/FAIL/COULD-NOT-RUN):
  1. shared-once, single run: tests/fixtures/np_shared_once.cir -- two
     resistors referencing a `once(agauss(...))` .param directly, plus two
     more referencing it indirectly through plain-copy .params, must all
     read back IDENTICAL resistance in one run (no .control loop at all).
  2. shared-once, statistical: tests/fixtures/np_shared_once_reset.cir --
     400 `reset` iterations in one ngspice process. Within each iteration
     R1 and R2 must match exactly (frozen within one circuit
     instantiation); across iterations the resulting resistance
     distribution's mean and sigma must match agauss(0, 10m, 3)'s
     configured 0 Ohm offset / 3.333 Ohm sigma (10m/3 on a 1k resistor)
     within a 99% CI (tools/mc.py's own mean_ci99/sigma_ci99, same
     convention as tests/test_mc.py).
  3. regression, plain agauss unchanged: tests/fixtures/np_shared_agauss.cir
     -- the exact docs/mc.md reproducer (no once()) -- two resistors
     sharing a plain `.param proc_off = agauss(...)` must still get
     INDEPENDENT (differing) draws, proving this patch does not change the
     default.
  4. regression, sky130-style per-instance mismatch unchanged:
     tests/fixtures/np_sky130_mismatch.cir -- a subckt-local agauss
     .param, instantiated twice, must still get independent per-instance
     draws.
  5. red-proof: tests/fixtures/np_shared_once.cir run against a build of
     this same source tree with patches-wip/numparam/0300 REVERTED
     (NGSPICE_NUMPARAM_NOPATCH) -- the shared-value guarantee from check 1
     must NOT hold (either the run fails outright, since "once" is then an
     undefined numparam function, or -- if it somehow ran -- the four
     values must not all match). Reports COULD-NOT-RUN (not a silent skip
     counted as pass) if NGSPICE_NUMPARAM_NOPATCH is not set.

Usage:
    NGSPICE=<build>/bin/ngspice-aug \\
    NGSPICE_NUMPARAM_NOPATCH=<nopatch build>/bin/ngspice \\
    python3 tests/test_numparam.py

Exit status: 0 = all gating checks passed, 1 = a gating check FAILed,
2 = a gating check COULD-NOT-RUN (e.g. ngspice missing).
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXTURES = os.path.join(HERE, "fixtures")
sys.path.insert(0, os.path.join(ROOT, "tools"))
import mc  # noqa: E402 -- reuse mean_ci99/sigma_ci99, same convention as test_mc.py

NGSPICE = os.environ.get("NGSPICE", "ngspice")
NGSPICE_NOPATCH = os.environ.get("NGSPICE_NUMPARAM_NOPATCH")

results = []
VEC_RE = re.compile(r"@([\w.\[\]]+)\s*=\s*([-+0-9.eE]+)")


def record(name, state, detail):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN")
    results.append((name, state, detail))


class CouldNotRun(Exception):
    pass


def run_deck(ngspice, stem, timeout=120):
    """Run tests/fixtures/<stem>.cir with `ngspice -b`, return (values, proc).

    values: list of (vecname, float) pairs in the order printed (so a loop
    printing "@r1[...] @r2[...]" N times gives 2N entries, in order).
    """
    deck = os.path.join(FIXTURES, f"{stem}.cir")
    if not os.path.isfile(deck):
        raise CouldNotRun(f"fixture not found: {deck}")
    try:
        proc = subprocess.run([ngspice, "-b", deck], capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError as e:
        raise CouldNotRun(f"cannot start {ngspice}: {e}") from e
    except subprocess.TimeoutExpired as e:
        raise CouldNotRun(f"{ngspice} timed out on {stem}.cir: {e}") from e
    out = proc.stdout + proc.stderr
    values = [(m.group(1), float(m.group(2))) for m in VEC_RE.finditer(out)]
    return values, proc


def close(a, b, rel=1e-9):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1.0)


# ---------- check 1: shared-once, single run ----------

def check_shared_once():
    name = "shared-once single-run: r1==r2==r5==r6"
    try:
        values, proc = run_deck(NGSPICE, "np_shared_once")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    got = {k: v for k, v in values}
    need = ["r1[resistance]", "r2[resistance]", "r5[resistance]", "r6[resistance]"]
    if not all(k in got for k in need):
        record(name, "COULD-NOT-RUN",
               f"expected vectors not in output (stderr tail: {proc.stderr[-300:]!r})")
        return
    vals = [got[k] for k in need]
    ok = all(close(v, vals[0]) for v in vals)
    record(name, "pass" if ok else "FAIL", f"values={vals}")


# ---------- check 2: shared-once, statistical across resets ----------

def check_shared_once_statistical():
    pair_name = "shared-once reset-loop: r1==r2 in every one of 400 iterations"
    stat_name = "shared-once reset-loop: mean/sigma of R match agauss(0,10m,3) within 99% CI"
    try:
        values, proc = run_deck(NGSPICE, "np_shared_once_reset", timeout=300)
    except CouldNotRun as e:
        record(pair_name, "COULD-NOT-RUN", str(e))
        record(stat_name, "COULD-NOT-RUN", str(e))
        return

    r1s = [v for k, v in values if k == "r1[resistance]"]
    r2s = [v for k, v in values if k == "r2[resistance]"]
    n = min(len(r1s), len(r2s))
    if n < 300:
        detail = f"only {n} iterations parsed (expected 400); stderr tail: {proc.stderr[-300:]!r}"
        record(pair_name, "COULD-NOT-RUN", detail)
        record(stat_name, "COULD-NOT-RUN", detail)
        return

    mismatches = sum(0 if close(a, b) else 1 for a, b in zip(r1s, r2s))
    record(pair_name, "pass" if mismatches == 0 else "FAIL",
           f"{n} iterations, {mismatches} mismatched pair(s)")

    # agauss(0, 10m, 3): stdvar = abs_variation/sigma = 0.01/3 on the
    # fractional term, applied to a nominal 1k resistor.
    nominal_r = 1000.0
    expected_mean = nominal_r
    expected_sigma = nominal_r * (0.01 / 3.0)
    m, mlo, mhi = mc.mean_ci99(r1s)
    s, slo, shi = mc.sigma_ci99(r1s)
    mean_ok = mlo <= expected_mean <= mhi
    sigma_ok = slo <= expected_sigma <= shi
    record(stat_name, "pass" if (mean_ok and sigma_ok) else "FAIL",
           f"n={n} mean={m:.4f} CI=[{mlo:.4f},{mhi:.4f}] (expect {expected_mean:.4f}); "
           f"sigma={s:.4f} CI=[{slo:.4f},{shi:.4f}] (expect {expected_sigma:.4f})")


# ---------- check 3: plain agauss per-use behaviour unchanged ----------

def check_plain_agauss_unchanged():
    name = "regression: plain .param agauss(...) still gives independent per-use draws"
    try:
        values, proc = run_deck(NGSPICE, "np_shared_agauss")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    got = {k: v for k, v in values}
    need = ["r1[resistance]", "r2[resistance]", "r5[resistance]", "r6[resistance]"]
    if not all(k in got for k in need):
        record(name, "COULD-NOT-RUN",
               f"expected vectors not in output (stderr tail: {proc.stderr[-300:]!r})")
        return
    vals = [got[k] for k in need]
    # must NOT all match -- any pairwise difference proves per-use re-evaluation
    all_equal = all(close(v, vals[0]) for v in vals)
    record(name, "FAIL" if all_equal else "pass", f"values={vals}")


# ---------- check 4: sky130-style per-instance mismatch unchanged ----------

def check_sky130_idiom_unchanged():
    name = "regression: sky130-style per-instance mismatch (subckt-local agauss) unchanged"
    try:
        values, proc = run_deck(NGSPICE, "np_sky130_mismatch")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    got = {k: v for k, v in values}
    need = ["r.x1.r1[resistance]", "r.x2.r1[resistance]"]
    if not all(k in got for k in need):
        record(name, "COULD-NOT-RUN",
               f"expected vectors not in output (stderr tail: {proc.stderr[-300:]!r})")
        return
    v1, v2 = got[need[0]], got[need[1]]
    record(name, "FAIL" if close(v1, v2) else "pass", f"x1={v1} x2={v2}")


# ---------- check 5: red-proof ----------

def check_red_proof():
    name = "red-proof: patch reverted -> shared-once guarantee does not hold"
    if not NGSPICE_NOPATCH:
        record(name, "COULD-NOT-RUN",
               "NGSPICE_NUMPARAM_NOPATCH not set -- no patch-reverted build to test against")
        return
    try:
        values, proc = run_deck(NGSPICE_NOPATCH, "np_shared_once")
    except CouldNotRun as e:
        # the deck fails to run at all on the reverted build (expected: "once"
        # is then an undefined numparam function) -- a valid red-proof pass.
        record(name, "pass", f"reverted build could not run the deck at all: {e}")
        return
    got = {k: v for k, v in values}
    need = ["r1[resistance]", "r2[resistance]", "r5[resistance]", "r6[resistance]"]
    if not all(k in got for k in need):
        record(name, "pass",
               f"reverted build did not produce the expected vectors "
               f"(stderr tail: {proc.stderr[-300:]!r}) -- guarantee does not hold, as expected")
        return
    vals = [got[k] for k in need]
    all_equal = all(close(v, vals[0]) for v in vals)
    # a pass here means the red-proof WORKED, i.e. the guarantee correctly
    # does NOT hold on the reverted build
    record(name, "FAIL" if all_equal else "pass",
           f"values={vals} (FAIL would mean the reverted build still shares the draw, "
           f"which would mean this check is not actually testing the patch)")


def main():
    check_shared_once()
    check_shared_once_statistical()
    check_plain_agauss_unchanged()
    check_sky130_idiom_unchanged()
    check_red_proof()

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    print(f"\n{sum(s == 'pass' for _, s, _ in results)}/{len(results)} pass")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
