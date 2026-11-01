#!/usr/bin/env python3
"""Validation for fra_ref (spec) — verify.md step 6, card #452, 2026-09-27.

Every check is pass / FAIL / COULD-NOT-RUN. tools/fra.py mode "loop" (transient Tian
two-injection FRA) is checked against the same independent references test_stb.py
already uses for stb_script: r1a's closed-form return ratio and r1b's Python nodal
Y-parameter solve (both imported read-only from test_stb, per card #452's instructions;
tests/test_stb.py itself is not modified here and belongs to another worker).

Frequency choice: the two-injection formula T = S/(1-S) (stbcore.tian) is only as
accurate as (1-S) is far from zero. Deep in R1a/R1b's high-gain passband (DC gain
~50 dB for r1a, ~93 dB for r1b) S is within a fraction of a percent of 1, and a
fraction-of-a-percent numerical error in the transient/DFT extraction of V1 becomes a
large error in T — this was measured directly while building fra.py (card #452 report)
and is intrinsic to the two-injection method, not a harness bug (a real network analyzer
has the same conditioning problem). So the frequency lists below span "below and around
crossover" in the sense that matters for stability work (loop gain within a few decades
of 0 dB), not the flat DC-gain plateau; that limit is written up in the card #452 report.

    python3 tests/test_fra.py        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import fra  # noqa: E402
import stb  # noqa: E402
import stbcore  # noqa: E402

sys.path.insert(0, HERE)
import test_stb as ts  # noqa: E402  (read-only reuse: r1a_T, r1b_T, params)

OUT = os.path.join(fra.DATA_ROOT, "tests")
WORKERS = min(4, max(1, (os.cpu_count() or 2) // 2))
DB_TOL = 0.1
DEG_TOL = 1.0
results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))


def db_deg(t):
    return stbcore.db(t), math.degrees(math.atan2(t.imag, t.real))


def compare_point(name, got, ref, db_tol=DB_TOL, deg_tol=DEG_TOL):
    gdb, gph = db_deg(got)
    rdb, rph = db_deg(ref)
    ddb, dph = gdb - rdb, gph - rph
    while dph > 180:
        dph -= 360
    while dph < -180:
        dph += 360
    ok = abs(ddb) <= db_tol and abs(dph) <= deg_tol
    record(name, ok, f"got {gdb:.4f} dB / {gph:.4f} deg, ref {rdb:.4f} dB / {rph:.4f} deg "
                      f"(d={ddb:+.4f} dB / {dph:+.4f} deg, bound {db_tol} dB / {deg_tol} deg)")
    return ddb, dph


def expect_bad_point(name, got, ref, db_tol=DB_TOL, deg_tol=DEG_TOL):
    gdb, gph = db_deg(got)
    rdb, rph = db_deg(ref)
    ddb, dph = gdb - rdb, gph - rph
    while dph > 180:
        dph -= 360
    while dph < -180:
        dph += 360
    bad = abs(ddb) > db_tol or abs(dph) > deg_tol
    record(name, bad, f"got {gdb:.4f} dB / {gph:.4f} deg, ref {rdb:.4f} dB / {rph:.4f} deg "
                       f"(d={ddb:+.4f} dB / {dph:+.4f} deg) must exceed the {db_tol} dB / {deg_tol} deg bound")


# Frequencies "below and around crossover", chosen per the module docstring above.
R1A_DECK = os.path.join(ROOT, "refs", "r1a_unilateral.cir")
R1A_FREQS = [3e5, 1e6, 3e6, 1e7, 3e7]          # r1a crossover ~1.1e7 Hz
R1B_DECK = os.path.join(ROOT, "refs", "r1b_bilateral.cir")
R1B_FREQS = [3e6, 1e7, 3e7, 1e8]               # r1b crossover ~5e7 Hz
PERIODIC_DECK = os.path.join(HERE, "fixtures", "r1a_periodic.cir")

COMMON = dict(settle_time=0.0, settle_periods=20, measure_periods=10, points_per_period=256)


def run_point(deck, f, amplitude=1e-3, outdir_name="pt", **kw):
    lines = stb.prepare_deck(deck)
    d = dict(COMMON)
    d.update(kw)
    return fra.run_loop_point(lines, "XP", f, amplitude, outdir=os.path.join(OUT, outdir_name), **d)


def main():
    t_start = time.time()
    try:
        # ---- R1a: FRA vs closed form (test_stb.r1a_T) ----
        p = ts.params(R1A_DECK)
        pts = fra.sweep_loop(R1A_DECK, "XP", R1A_FREQS, 1e-3, outdir_root=os.path.join(OUT, "r1a_sweep"),
                              workers=WORKERS, **COMMON)
        r1a_T = {}
        for f, (status, val) in zip(R1A_FREQS, pts):
            if status != "ok":
                record(f"R1a FRA @ {f:.4g} Hz", None, val)
                continue
            ref = ts.r1a_T(p, f)
            r1a_T[f] = val["T"]
            compare_point(f"R1a FRA vs closed form @ {f:.4g} Hz", val["T"], ref)

        # ---- combined Tian/Middlebrook estimate in R1a's high-|T| flat region (task T2.3
        #      continuation, card #452): T = S/(1-S) is only as accurate as (1-S) is far
        #      from 0 (module docstring "Ill-conditioning at high loop gain"), and r1a has
        #      no reverse transmission (VCCS input draws no current, refs/r1a_unilateral.cir
        #      header), so Middlebrook and Tian are analytically the same return ratio there
        #      -- any gap between the FRA's numerically-extracted Tian and the closed form is
        #      extraction error, not a different formula. 100 Hz-10 kHz sits on r1a's flat
        #      ~50.46 dB DC-ish plateau (well below the fLC~=159 kHz/fESR~=796 kHz poles),
        #      exactly the region tests/test_fra.py's R1A_FREQS above deliberately avoids.
        #      Uses fra.py's DEFAULT points_per_period (COMMON sets 256, the default; not
        #      overridden here). ----
        R1A_HIGHGAIN_FREQS = [100.0, 300.0, 1e3, 3e3, 1e4]
        pts = fra.sweep_loop(R1A_DECK, "XP", R1A_HIGHGAIN_FREQS, 1e-3,
                              outdir_root=os.path.join(OUT, "r1a_highgain"), workers=WORKERS, **COMMON)
        for f, (status, val) in zip(R1A_HIGHGAIN_FREQS, pts):
            if status != "ok":
                record(f"R1a high-|T| combined FRA @ {f:.4g} Hz", None, val)
                continue
            ref = ts.r1a_T(p, f)
            compare_point(f"R1a high-|T| combined (Tian/Middlebrook) vs closed form @ {f:.4g} Hz",
                           val["combined"], ref)
            record(f"R1a high-|T| combined picks Middlebrook above threshold @ {f:.4g} Hz",
                   val["combined_source"] == "middlebrook",
                   f"|T|={stbcore.db(val['T']):.4f} dB, source={val['combined_source']}, "
                   f"threshold={fra.DEFAULT_COMBINE_THRESHOLD_DB} dB, "
                   f"disagreement(Middlebrook-Tian)={val['disagreement_db']:+.4f} dB / "
                   f"{val['disagreement_deg']:+.4f} deg")
            # red-proof: the *uncombined* Tian estimate from this same run must fail the
            # 0.1 dB / 1 deg bound here -- proves the Middlebrook combination is load-bearing
            # (fixes a real failure), not a no-op that happens to agree either way.
            expect_bad_point(f"red-proof: plain Tian (no Middlebrook) fails the high-|T| bound @ {f:.4g} Hz",
                              val["T"], ref)

        # ---- R1b: FRA vs the Python nodal Y-parameter reference (test_stb.r1b_T) ----
        p = ts.params(R1B_DECK)
        pts = fra.sweep_loop(R1B_DECK, "XP", R1B_FREQS, 1e-3, outdir_root=os.path.join(OUT, "r1b_sweep"),
                              workers=WORKERS, **COMMON)
        for f, (status, val) in zip(R1B_FREQS, pts):
            if status != "ok":
                record(f"R1b FRA @ {f:.4g} Hz", None, val)
                continue
            ref = ts.r1b_T(p, f)
            compare_point(f"R1b FRA vs nodal Y reference @ {f:.4g} Hz", val["T"], ref)

        # ---- amplitude-halving small-signal check (spec fra_ref) ----
        for f in (1e6, 1e7):
            ta, tb, dbdiff = fra.amplitude_halving_check(
                R1A_DECK, "XP", f, 1e-3, outdir_root=os.path.join(OUT, "halving"), **COMMON)
            record(f"amplitude-halving |dB diff| <= 0.1 @ {f:.4g} Hz", abs(dbdiff) <= 0.1,
                   f"A={ta:.6f} A/2={tb:.6f} diff={dbdiff:+.5f} dB")

        # ---- red-proof (c): sign flip must fail the comparison ----
        f = 1e6
        got = r1a_T.get(f)
        if got is not None:
            ref = ts.r1a_T(ts.params(R1A_DECK), f)
            expect_bad_point("red-proof: sign-flipped T must fail the comparison", -got, ref)
        else:
            record("red-proof: sign-flipped T must fail the comparison", None, "no R1a result at 1e6 Hz to flip")

        # ---- red-proof (a): a non-integer number of excitation periods must fail ----
        f = 1e6
        ref = ts.r1a_T(ts.params(R1A_DECK), f)
        bad = run_point(R1A_DECK, f, outdir_name="redproof_a", measure_periods=10.37)
        expect_bad_point("red-proof: non-integer measurement window (10.37 periods) must fail",
                          bad["T"], ref)
        good = run_point(R1A_DECK, f, outdir_name="redproof_a_control", measure_periods=10)
        compare_point("control: integer measurement window (10 periods) still passes", good["T"], ref)

        # ---- red-proof (b1): measuring before the sinusoidal response has settled ----
        bad = run_point(R1A_DECK, f, outdir_name="redproof_b1", settle_periods=0, measure_periods=1)
        expect_bad_point("red-proof: unsettled measurement window (0 settle periods) must fail",
                          bad["T"], ref)

        # ---- red-proof (b2): skipping the perturbed-minus-unperturbed subtraction on a
        #      deck with a large-signal periodic component unrelated to the probe ----
        ref_periodic = ts.r1a_T(ts.params(PERIODIC_DECK), f)
        good = run_point(PERIODIC_DECK, f, outdir_name="redproof_b2_sub", subtract_baseline=True)
        compare_point("control: subtracting the unperturbed run on the periodic-fixture deck passes",
                       good["T"], ref_periodic)
        bad = run_point(PERIODIC_DECK, f, outdir_name="redproof_b2_nosub", subtract_baseline=False)
        expect_bad_point("red-proof: skipping the subtraction on the periodic-fixture deck must fail",
                          bad["T"], ref_periodic)

        # ---- red-proof: could-not-run must never read as a pass ----
        saved = fra.NGSPICE
        stb.NGSPICE = "/nonexistent/ngspice"
        fra.NGSPICE = "/nonexistent/ngspice"
        try:
            run_point(R1A_DECK, f, outdir_name="redproof_bogus")
            record("red-proof: missing ngspice gives COULD-NOT-RUN", False, "run returned a result")
        except stb.CouldNotRun as e:
            record("red-proof: missing ngspice gives COULD-NOT-RUN", True, str(e)[:60])
        finally:
            stb.NGSPICE = saved
            fra.NGSPICE = saved

    except stb.CouldNotRun as e:
        record("run", None, str(e))

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    print(f"\n{sum(s == 'pass' for _, s, _ in results)}/{len(results)} pass"
          f"  ({time.time() - t_start:.1f}s)")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
