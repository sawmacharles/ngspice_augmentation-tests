#!/usr/bin/env python3
"""Validation for tools/fra.py's adaptive settle (card #487, 2026-09-29).

`tools/fra.py` mode "loop"/"transfer" used to settle the UNPERTURBED circuit for a
fixed, caller-guessed `--settle` time before injection turned on. Card #484/#487 found
that guess too short for slow converters (R3's 300us guess left a 6.15% periodicity
residual -- max|v(t)-v(t-T)|/ripple over one switching period T -- biasing R3's FRA
curve by up to 3.7dB/17deg near fs/2; R3 actually needs ~9ms to reach <=1%). This file
checks the new `periodicity_residual`/`periodicity_scan`/`sweep_loop_adaptive`/
`sweep_transfer_adaptive` (opt-in via `--adaptive-settle`; the old fixed-`--settle`
path is completely unchanged and stays every other fra.py test's default) against:

  1. `periodicity_residual`/`_interp_linear` on synthetic, hand-built data (no ngspice,
     fast, deterministic) -- a perfectly periodic signal must read ~0 residual, a
     signal still converging must read a large one, and a flat signal must not divide
     by ~0.
  2. `periodicity_scan` on refs/r2_switched_rc.cir (R2, a genuine LPTV circuit,
     fclk=100kHz/period=10us -- chosen because it is fast, ~0.1s/run, docs/
     checkpoint.md) must find a settle time well under the 300us (30-period) fixed
     settle test_fra_lptv.py already validated as "well converged" for this circuit,
     with residual <= the requested bound; a red-proof with an impossibly tight bound
     must honestly report `capped=True` (never a fabricated pass), and a cap too short
     to cover even 2 periods must raise CouldNotRun rather than silently misbehaving.
  3. `sweep_transfer_adaptive` end to end on R2 vs the SAME fixed-300us-settle FRA
     result (not re-deriving r2_analytic.g0() independently -- that is
     tests/test_fra_lptv.py's job and is reused read-only here): adaptive settle
     (bound=1%) must land close to the well-validated fixed-settle answer, proving the
     new machinery picks a genuinely-equivalent physical state, not just "a time".
  4. The CLI (`--adaptive-settle`/`--settle-period`/`--residual-bound`/`--settle-cap`)
     end to end on refs/r1a_unilateral.cir (mode loop, LTI -- no real periodicity, so
     the scan should settle "immediately" at the first period checked, residual ~0)
     and refs/r2_switched_rc.cir (mode transfer): fra.csv/summary.json must carry the
     new settle_time_s/settle_residual/settle_capped columns; the OLD fixed-settle CLI
     path (no --adaptive-settle) must produce the byte-identical old CSV header (no
     new columns) -- the compatibility claim in write_loop_results' own docstring.

    python3 tests/test_fra_settle.py        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import fra  # noqa: E402
import stb  # noqa: E402
import stbcore  # noqa: E402

OUT = os.path.join(fra.DATA_ROOT, "tests_settle")
results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))
    print(f"{state:13s} {name}  {detail}")


# ---------------------------------------------------------------------------
# 1. periodicity_residual / _interp_linear on synthetic data (no ngspice)
# ---------------------------------------------------------------------------

def check_synthetic():
    period = 1.0
    n = 400
    times = [k * period / n for k in range(3 * n + 1)]

    # perfectly periodic: three identical cycles of a sine -> residual ~0
    perfect = [math.sin(2 * math.pi * t / period) for t in times]
    resid = fra.periodicity_residual(times, perfect, period, 3 * period)
    record("synthetic: perfectly periodic signal -> residual ~0", resid < 1e-9,
           f"residual={resid:.3e}")

    # still converging (exponential approach to a periodic limit): the LAST cycle
    # should be closer to periodic than comparing cycle 0 vs cycle 1 would be
    tau = period  # decay time comparable to one period -> not yet settled by cycle 2-3
    converging = [math.sin(2 * math.pi * t / period) * (1 - math.exp(-t / tau)) + 0.3 * math.exp(-t / tau)
                  for t in times]
    resid_early = fra.periodicity_residual(times, converging, period, 2 * period)
    resid_late = fra.periodicity_residual(times, converging, period, 3 * period)
    record("synthetic: still-converging signal reads a large residual early on",
           resid_early > 0.05, f"residual@2T={resid_early:.4f}")
    record("synthetic: residual shrinks as the signal settles (late < early)",
           resid_late < resid_early, f"residual@2T={resid_early:.4f} residual@3T={resid_late:.4f}")

    # flat signal: ripple ~0, must not raise/produce inf/nan (falls back to abs residual)
    flat = [1.2345 for _ in times]
    resid_flat = fra.periodicity_residual(times, flat, period, 3 * period)
    record("synthetic: flat signal (ripple~0) does not blow up", resid_flat < 1e-9,
           f"residual={resid_flat:.3e}")

    # window before the available trace must raise CouldNotRun, not silently clip
    try:
        fra.periodicity_residual(times, perfect, period, 0.5 * period)
        record("synthetic: window before the trace start raises CouldNotRun", False, "did not raise")
    except stb.CouldNotRun as e:
        record("synthetic: window before the trace start raises CouldNotRun", True, str(e)[:70])


# ---------------------------------------------------------------------------
# 2. periodicity_scan on R2 (real ngspice, but fast: ~0.1s/run per docs/checkpoint.md)
# ---------------------------------------------------------------------------

R2_DECK = os.path.join(ROOT, "refs", "r2_switched_rc.cir")
R2_PERIOD = 10e-6  # 1/fclk, fclk=100kHz


def check_r2_scan():
    deck_lines = stb.prepare_deck(R2_DECK)
    try:
        info = fra.periodicity_scan(deck_lines, ["v(vx)"], R2_PERIOD, residual_bound=0.01,
                                     cap=300e-6, points_per_period=64, startup="uic",
                                     outdir=os.path.join(OUT, "r2_scan"))
    except stb.CouldNotRun as e:
        record("R2 periodicity_scan finds settle_time < 300us at 1% bound", None, str(e))
        return
    record("R2 periodicity_scan finds settle_time < 300us at 1% bound",
           not info["capped"] and info["settle_time"] < 300e-6,
           f"settle_time={info['settle_time']:.3g}s residual={info['residual']:.4%} "
           f"periods_scanned={info['periods_scanned']} capped={info['capped']}")
    record("R2 periodicity_scan residual is at or under the requested bound",
           info["residual"] <= 0.01 + 1e-9 if not info["capped"] else None,
           f"residual={info['residual']:.4%}")

    # red-proof: an impossibly tight bound must honestly report capped, not a fake pass
    try:
        info_tight = fra.periodicity_scan(deck_lines, ["v(vx)"], R2_PERIOD, residual_bound=1e-12,
                                           cap=300e-6, points_per_period=64, startup="uic",
                                           outdir=os.path.join(OUT, "r2_scan_tight"))
        record("red-proof: impossibly tight residual_bound reports capped=True (never a fake pass)",
               info_tight["capped"] and info_tight["residual"] > 1e-12,
               f"capped={info_tight['capped']} residual={info_tight['residual']:.3e} "
               f"settle_time={info_tight['settle_time']:.3g}s (== cap)")
    except stb.CouldNotRun as e:
        record("red-proof: impossibly tight residual_bound reports capped=True (never a fake pass)",
               None, str(e))

    # red-proof: a cap too short to cover 2 periods must raise CouldNotRun, not misbehave
    try:
        fra.periodicity_scan(deck_lines, ["v(vx)"], R2_PERIOD, residual_bound=0.01,
                              cap=1.5e-6, points_per_period=64, startup="uic",
                              outdir=os.path.join(OUT, "r2_scan_shortcap"))
        record("red-proof: a cap under 2 periods raises CouldNotRun", False, "did not raise")
    except stb.CouldNotRun as e:
        record("red-proof: a cap under 2 periods raises CouldNotRun", True, str(e)[:80])


# ---------------------------------------------------------------------------
# 3. sweep_transfer_adaptive vs the SAME fixed-300us-settle result (R2)
# ---------------------------------------------------------------------------

def check_r2_transfer_adaptive_vs_fixed():
    f = 10e3       # fclk/10 -- exact commensurate ratio, matches test_fra_lptv.py's own R2 grid
    measure_periods = 1
    ppp = 1000     # dt=100ns, matches test_fra_lptv.py's DT_TARGET convention

    try:
        fixed = fra.run_transfer_point(stb.prepare_deck(R2_DECK), "vin", "v(vx)", f, 1e-3,
                                        settle_time=300e-6, settle_periods=20.0 * f / 1e5,
                                        measure_periods=measure_periods, points_per_period=ppp,
                                        outdir=os.path.join(OUT, "r2_fixed"))
    except stb.CouldNotRun as e:
        record("R2 adaptive-settle transfer FRA agrees with the fixed-300us-settle result", None, str(e))
        return

    try:
        points, settle_info = fra.sweep_transfer_adaptive(
            R2_DECK, "vin", "v(vx)", [f], 1e-3, R2_PERIOD, residual_bound=0.01, cap=300e-6,
            outdir_root=os.path.join(OUT, "r2_adaptive"), workers=1,
            settle_periods=20.0 * f / 1e5, measure_periods=measure_periods, points_per_period=ppp)
    except stb.CouldNotRun as e:
        record("R2 adaptive-settle transfer FRA agrees with the fixed-300us-settle result", None, str(e))
        return
    status, val = points[0]
    if status != "ok":
        record("R2 adaptive-settle transfer FRA agrees with the fixed-300us-settle result", None, val)
        return
    got_db, got_deg = stbcore.db(val["H"]), math.degrees(__import__("cmath").phase(val["H"]))
    ref_db, ref_deg = stbcore.db(fixed["H"]), math.degrees(__import__("cmath").phase(fixed["H"]))
    ddb, ddeg = got_db - ref_db, got_deg - ref_deg
    while ddeg > 180: ddeg -= 360
    while ddeg < -180: ddeg += 360
    ok = abs(ddb) <= 0.5 and abs(ddeg) <= 3.0
    record("R2 adaptive-settle transfer FRA agrees with the fixed-300us-settle result", ok,
           f"adaptive settle_time={settle_info['settle_time']:.3g}s residual={settle_info['residual']:.4%} "
           f"capped={settle_info['capped']}; got {got_db:.4f}dB/{got_deg:.4f}deg vs fixed-300us "
           f"{ref_db:.4f}dB/{ref_deg:.4f}deg (d={ddb:+.4f}dB/{ddeg:+.4f}deg, bound 0.5dB/3deg)")
    record("R2 adaptive settle used materially less settle time than the fixed 300us guess",
           settle_info["settle_time"] < 300e-6, f"adaptive={settle_info['settle_time']:.3g}s vs fixed=300e-6s")


# ---------------------------------------------------------------------------
# 4. CLI end to end
# ---------------------------------------------------------------------------

R1A_DECK = os.path.join(ROOT, "refs", "r1a_unilateral.cir")


def run_cli(args, timeout=120):
    cmd = [sys.executable, os.path.join(ROOT, "tools", "fra.py")] + args
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return proc.returncode, proc.stdout + proc.stderr


def check_cli():
    # old fixed-settle path: CSV header must be byte-identical to the pre-#487 shape
    outdir_old = os.path.join(OUT, "cli_old")
    rc, out = run_cli([R1A_DECK, "--probe", "XP", "--freqs", "1e6,3e6", "--amplitude", "1e-3",
                        "--out", outdir_old])
    if rc not in (0, 2):
        record("CLI: old fixed-settle path still runs (unchanged)", None, out[-800:])
    else:
        with open(os.path.join(outdir_old, "fra.csv")) as fh:
            header = fh.readline().strip()
        expected = ("freq,status,T_re,T_im,T_db,T_phase_deg,loopGain_db,loopGain_phase_deg,"
                    "middlebrook_db,middlebrook_phase_deg,combined_db,combined_phase_deg,"
                    "combined_source,disagreement_db,disagreement_deg,detail")
        record("CLI: old fixed-settle path's fra.csv header is byte-identical to pre-#487",
               header == expected, header)

    # new adaptive path, mode loop, LTI deck (R1a): should settle "immediately" (residual ~0)
    outdir_new = os.path.join(OUT, "cli_adaptive_loop")
    rc, out = run_cli([R1A_DECK, "--probe", "XP", "--freqs", "1e6,3e6", "--amplitude", "1e-3",
                        "--adaptive-settle", "--settle-period", "1e-6", "--residual-bound", "0.01",
                        "--settle-cap", "1e-4", "--out", outdir_new])
    ok_exit = rc in (0, 2)
    record("CLI: --adaptive-settle mode loop runs and exits 0/2", ok_exit, out[-500:] if not ok_exit else out.strip().splitlines()[-1] if out.strip() else "")
    if ok_exit:
        with open(os.path.join(outdir_new, "fra.csv")) as fh:
            header = fh.readline().strip()
        record("CLI: --adaptive-settle mode loop's fra.csv carries the new settle columns",
               "settle_time_s" in header and "settle_residual" in header and "settle_capped" in header,
               header)
        import json
        with open(os.path.join(outdir_new, "summary.json")) as fh:
            summ = json.load(fh)
        record("CLI: --adaptive-settle mode loop's summary.json carries a 'settle' key",
               "settle" in summ and "settle_time" in summ["settle"], str(summ.get("settle")))

    # new adaptive path, mode transfer, R2
    outdir_tr = os.path.join(OUT, "cli_adaptive_transfer")
    rc, out = run_cli([R2_DECK, "--mode", "transfer", "--source", "vin", "--output", "v(vx)",
                        "--freqs", "10000", "--amplitude", "1e-3", "--adaptive-settle",
                        "--settle-period", str(R2_PERIOD), "--residual-bound", "0.01",
                        "--settle-cap", "3e-4", "--settle-periods", "2", "--measure-periods", "1",
                        "--points-per-period", "1000", "--out", outdir_tr])
    ok_exit = rc in (0, 2)
    record("CLI: --adaptive-settle mode transfer runs and exits 0/2", ok_exit, out[-500:] if not ok_exit else "")
    if ok_exit:
        with open(os.path.join(outdir_tr, "fra.csv")) as fh:
            header = fh.readline().strip()
        record("CLI: --adaptive-settle mode transfer's fra.csv carries the new settle columns",
               "settle_time_s" in header, header)


def main():
    t_start = time.time()
    os.makedirs(OUT, exist_ok=True)
    check_synthetic()
    check_r2_scan()
    check_r2_transfer_adaptive_vs_fixed()
    check_cli()

    width = max(len(n) for n, _, _ in results)
    print("\n" + "=" * 70)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}")
    states = {s for _, s, _ in results}
    n_pass = sum(s == "pass" for _, s, _ in results)
    print(f"\n{n_pass}/{len(results)} pass  ({time.time() - t_start:.1f}s)")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
