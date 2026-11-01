#!/usr/bin/env python3
"""Validation for the converter (switching-deck) FRA fix — card #452, T2.3, 2026-09-28.

Card #452's converter investigation (see the card report) found that tools/fra.py's
amplitude-halving small-signal check failed badly on R3 (refs/r3_vm_buck.cir, voltage-
mode buck) and R4 (refs/r4_pcm_buck.cir, peak-current-mode buck): up to several dB /
tens of degrees of change when halving the injection amplitude, which should leave a
truly small-signal T unchanged to first order. Root cause: `run_loop_point`'s `tran`
call used `tmax = dt = 1/(f*points_per_period)`, which grows with 1/f and, below
roughly 13 kHz on these 500 kHz-switching decks, ends up *larger* than one switching
period — ngspice's solver could step clean over a switching edge, and the resulting
residual did not cancel between the baseline and perturbed transients (see
tools/fra.py's run_loop_point docstring, "tmax"). The fix: `run_loop_point` (and
`sweep_loop`) now take `tmax` and `spice_options` independently of `dt`; a small,
frequency-independent `tmax` (5 ns, well under R3/R4's 2 us switching period) plus a
tighter transient-truncation-error tolerance (`.options trtol=1`, ngspice's default is
7) fix it.

This file is the red-proofed evidence for that fix, at the goal's amplitude ceiling
(<=10 mV, card #452 T2.3): the FIXED settings must pass the halving check (<=0.1 dB /
1 deg) on R3 at three frequencies spanning 500 Hz-fs/4 and on R4 (D~=0.30) at two of
those three (R4's much longer settle_time, 3.5 ms vs R3's 300 us, makes its low-
frequency point expensive — see the card report for the full 3-frequency numbers,
measured but not repeated here to keep this file's runtime down); the OLD settings
(unmodified tmax/spice_options) must FAIL that same check at the low frequency on both
decks, where the old bug was worst.

Frequencies are taken from tools/fra_curves.py's own frequency list (`fra_curves.FREQS`
dec-5, 500 Hz-fs/2, and `halving_freqs`), not retyped as rounded literals: card #452's
investigation found that rounding a frequency to 6 significant figures (as SUMMARY.md's
table prints them) changes the halving numbers by several dB / tens of degrees on its
own — a real window/switching-alignment sensitivity (see the card report, item 5), not
a bug in this file, but it means every comparison here must use the bit-identical float
fra_curves.py's sweep uses.

    python3 tests/test_fra_converter.py     exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN

Runtime: ~3 minutes (measured; kept under the 5 minute budget, card #452 T2.3 instructions).
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
import fra_curves as fc  # noqa: E402  (frequency list + the fixed injection amplitude)

OUT = os.path.join(fra.DATA_ROOT, "tests_converter")
DB_TOL = 0.1
DEG_TOL = 1.0
# card #550 (2026-10-02): on R3 at the default tolerance (reltol 1e-3) the halving difference is tolerance-limited, not a harness property:
# each run's own LTE-driven step sequence differs with the settle time and, in the last bits, between machines. Measured with
# ~/Work/data/ngspice_augmentation/550/fra_spread.py (3 frequencies x settle 300/360/420 us x reltol default/1e-4, second workstation Arch libm and development
# workstation Ubuntu libm). Low f (792 Hz), default reltol: dev. ws. +0.0423 / +0.0952 / +0.0837 dB, second ws. +0.1021 / +0.0604 / +0.0682 dB (settle 300/360/420 us);
# mid f 0.020..0.048 (both); high f dev. ws. -0.0966 / +0.0438 / -0.0301, second ws. -0.0100 / -0.0534 / -0.0285. Same points at reltol=1e-4: low
# +0.0083 / +0.0083 / +0.0082 dB, mid -0.0018, high -0.012 / -0.012 / -0.011, identical on both machines to 3-4 digits and across the three
# settle times (T itself moves 0.17 dB between the two reltols: the default-tolerance value is the tolerance-limited one). So the spread is
# 0.06 dB at equal settle and 0.06 dB across settle times on either machine; the R3 bound is 0.15 dB (largest value seen 0.1021; 3x the 0.06
# spread would be 0.18). R4 measured <= 0.0021 dB on both and keeps 0.1. The red-proofs fail by 1.2 / 4.85 dB, an order of magnitude above
# either bound, and a new check at reltol=1e-4 holds the harness itself to 0.03 dB (measured 0.0083, 3.6x margin).
DB_TOL_R3 = 0.15
DB_TOL_CLEAN = 0.03
results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))


def db_deg(t):
    return stbcore.db(t), math.degrees(math.atan2(t.imag, t.real))


def halving_diff(ra, rb):
    """(dB diff, deg diff) of the combined Tian/Middlebrook estimate, A vs A/2."""
    adb, aph = db_deg(ra["combined"])
    bdb, bph = db_deg(rb["combined"])
    ddb, dph = adb - bdb, aph - bph
    while dph > 180:
        dph -= 360
    while dph < -180:
        dph += 360
    return ddb, dph


def run_halving(name, deck_lines, f, amplitude, settle_time, kw, outdir_root, expect_pass, db_tol=None):
    db_tol = DB_TOL if db_tol is None else db_tol
    ra = fra.run_loop_point(deck_lines, "XP", f, amplitude, settle_time=settle_time,
                             outdir=os.path.join(outdir_root, f"f_{f:.6g}_A"), **kw)
    rb = fra.run_loop_point(deck_lines, "XP", f, amplitude / 2, settle_time=settle_time,
                             outdir=os.path.join(outdir_root, f"f_{f:.6g}_A2"), **kw)
    ddb, dph = halving_diff(ra, rb)
    ok = abs(ddb) <= db_tol and abs(dph) <= DEG_TOL
    detail = (f"T(A)={db_deg(ra['combined'])[0]:.3f}dB T(A/2)={db_deg(rb['combined'])[0]:.3f}dB "
              f"diff={ddb:+.4f}dB / {dph:+.4f}deg (bound {db_tol}dB / {DEG_TOL}deg)")
    if expect_pass:
        record(name, ok, detail)
    else:
        record(name, not ok, detail + "  [red-proof: must exceed the bound]")
    return ddb, dph


# ---------------------------------------------------------------------------
# decks, frequencies, amplitude, settings
# ---------------------------------------------------------------------------
R3_DECK = os.path.join(ROOT, "refs", "r3_vm_buck.cir")
R3_SETTLE = 300e-6
R4_DECK = os.path.join(ROOT, "refs", "r4_pcm_buck.cir")  # D~=0.30 (vin_v=11, deck default)
R4_SETTLE = 3.5e-3

# fc.halving_freqs(fc.FREQS) = [freqs[1], freqs[n//2], freqs[-2]] -- low/mid/high
# representative points from fra_curves.py's own dec-5, 500 Hz-fs/2 sweep list, used
# bit-for-bit (see module docstring above on why rounding these is unsafe).
LOW_F, MID_F, HIGH_F = fc.halving_freqs(fc.FREQS)

AMPLITUDE = fc.AMPLITUDE  # 0.01 V (10 mV) -- the goal's ceiling, card #452 T2.3
assert AMPLITUDE <= 0.01 + 1e-15, f"goal is <=10mV injection, got {AMPLITUDE}"

FIXED_KW = dict(settle_periods=3, measure_periods=5, points_per_period=256, startup="uic",
                tmax=5e-9, spice_options=".options trtol=1\n")
OLD_KW = dict(settle_periods=3, measure_periods=5, points_per_period=256, startup="uic")
# (OLD_KW omits tmax/spice_options entirely -- run_loop_point's pre-card-#452-T2.3-fix
# default, tmax=dt, no extra .options -- exactly what generated the SUMMARY.md numbers
# this fix replaces.)


def main():
    t_start = time.time()
    try:
        r3_lines = stb.prepare_deck(R3_DECK)
        r4_lines = stb.prepare_deck(R4_DECK)

        # ---- FIXED settings must pass on R3 at 3 frequencies (500Hz-fs/4 span) ----
        for f, label in ((LOW_F, "low"), (MID_F, "mid"), (HIGH_F, "high")):
            run_halving(f"R3 (voltage-mode buck) halving <={DB_TOL_R3}dB/1deg @ {label} f={f:.6g} Hz [FIXED]",
                        r3_lines, f, AMPLITUDE, R3_SETTLE, FIXED_KW,
                        os.path.join(OUT, "r3_fixed"), expect_pass=True, db_tol=DB_TOL_R3)
        # tolerance-limited noise removed: the halving difference of the harness itself
        run_halving(f"R3 halving @ low f={LOW_F:.6g} Hz [FIXED, reltol=1e-4]: the harness's own amplitude dependence <= {DB_TOL_CLEAN}dB/1deg",
                    r3_lines, LOW_F, AMPLITUDE, R3_SETTLE, dict(FIXED_KW, spice_options=".options trtol=1 reltol=1e-4\n"),
                    os.path.join(OUT, "r3_fixed_tight"), expect_pass=True, db_tol=DB_TOL_CLEAN)

        # ---- FIXED settings must pass on R4 D~=0.30 at 2 of the 3 (mid, high -- R4's
        #      3.5ms settle_time makes the low-frequency point expensive; measured
        #      separately in the card #452 report at all 3, not repeated here to keep
        #      this file under the 5 minute budget) ----
        for f, label in ((MID_F, "mid"), (HIGH_F, "high")):
            run_halving(f"R4 D~=0.30 (peak-current-mode buck) halving <=0.1dB/1deg @ {label} f={f:.6g} Hz [FIXED]",
                        r4_lines, f, AMPLITUDE, R4_SETTLE, FIXED_KW,
                        os.path.join(OUT, "r4_fixed"), expect_pass=True)

        # ---- red-proof: the OLD settings (no tmax/spice_options fix) must FAIL the same
        #      check, same deck/frequency/amplitude, at the low frequency where the old
        #      bug was worst (tmax=dt >> switching period there) ----
        run_halving("red-proof: R3 OLD settings (tmax=dt, no trtol fix) fails @ low f="
                    f"{LOW_F:.6g} Hz", r3_lines, LOW_F, AMPLITUDE, R3_SETTLE, OLD_KW,
                    os.path.join(OUT, "r3_old_redproof"), expect_pass=False)
        run_halving("red-proof: R4 D~=0.30 OLD settings (tmax=dt, no trtol fix) fails @ low f="
                    f"{LOW_F:.6g} Hz", r4_lines, LOW_F, AMPLITUDE, R4_SETTLE, OLD_KW,
                    os.path.join(OUT, "r4_old_redproof"), expect_pass=False)

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
