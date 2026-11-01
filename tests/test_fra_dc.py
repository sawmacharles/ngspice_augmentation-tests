#!/usr/bin/env python3
"""Validation for tools/fra.py's source-DC-bias fix (card #485, 2026-09-29).

`tools/fra.py`'s mode "transfer" (`run_transfer_point`/`sweep_transfer`, used for
control-to-output, line-to-output and output-impedance measurements, spec.md
`pac_sb0` (c)) injects a SIN onto a named source via `alter @source[sin] = [ VO VA
FREQ TD THETA PHASE ]`. `_sin_spec()` hard-coded VO=0, so injecting on a source that
carries a DC value (e.g. Vin for a line-to-output measurement) silently erased that
bias for BOTH the "base" and "pert" transient runs -- not a small perturbation
around the real operating point, a different circuit entirely (Vin=0V instead of
its real value). Found by the Phase 5 worker: the published R3/R4 line-to-output
curves in ~/Work/data/ngspice_augmentation/fra/curves-fork-s10/SUMMARY.md were wrong
by 30-60 dB.

Fixed: `source_dc_value()` reads the source's own DC token from the deck text (and
refuses -- raises CouldNotRun -- if the source already carries a time function,
rather than silently replacing it); `_resolve_dc_numeric()` resolves a `{param}`
expression token to a plain number via a quick, analysis-free `print @source[dc]`
query (alter's bracket-literal SIN form does not evaluate brace expressions itself,
confirmed empirically); `run_transfer_point` now passes the resolved value through
as `_sin_spec`'s `vo`.

tests/fixtures/fra_dc_tanh.cir is a memoryless nonlinear stage (out=tanh(in)) built
for this file: a LINEAR circuit's small-signal response does not depend on the DC
bias at all (superposition), so it cannot expose this bug even when injected on a
biased source; tanh's local slope does, and has a two-line closed form (see the
fixture's own header), so this is an independent numeric check, not a repeat of
tools/fra.py's own internals.

    python3 tests/test_fra_dc.py        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import math
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import fra  # noqa: E402
import stb  # noqa: E402
import stbcore  # noqa: E402

FIXTURE = os.path.join(HERE, "fixtures", "fra_dc_tanh.cir")
OUT = os.path.join(fra.DATA_ROOT, "tests_fra_dc")
DB_TOL = 0.05   # tanh's slope is exact and the injection amplitude is small (1e-4);
                # first-order error is O(amplitude), so this should be a tight bound
results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))
    print(f"{state:<14}{name:<74}{detail}")


def tanh_gain(vdc):
    """Analytic H = d(tanh(v))/dv at v=vdc = sech^2(vdc) = 1 - tanh(vdc)^2."""
    return 1.0 - math.tanh(vdc) ** 2


def _deck_with_vdc(vdc, outpath):
    with open(FIXTURE) as fh:
        text = fh.read()
    text, n = re.subn(r"^\.param vdc=[-0-9.eE]+$", f".param vdc={vdc}", text, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError("fra_dc_tanh.cir's '.param vdc=' line not found/matched")
    with open(outpath, "w") as fh:
        fh.write(text)
    return outpath


AMPLITUDE = 1e-4
FREQ = 1e3  # arbitrary: H is w-independent for this memoryless circuit (fixture header)


def check_gain_at(name, vdc, run_transfer_fn, tol_db=DB_TOL):
    deck_path = _deck_with_vdc(vdc, os.path.join(OUT, f"vdc_{vdc}.cir"))
    deck_lines = stb.prepare_deck(deck_path)
    try:
        r = run_transfer_fn(deck_lines, "vin", "v(out)", FREQ, AMPLITUDE,
                             settle_time=0.0, settle_periods=2, measure_periods=5,
                             points_per_period=64, outdir=os.path.join(OUT, name.replace(" ", "_")))
    except stb.CouldNotRun as e:
        record(name, None, str(e))
        return None
    h_ref = tanh_gain(vdc)
    got_db, ref_db = stbcore.db(r["H"]), stbcore.db(h_ref)
    ddb = got_db - ref_db
    ok = abs(ddb) <= tol_db
    detail = f"got H={r['H'].real:+.6f}{r['H'].imag:+.6f}j ({got_db:.3f}dB), ref={h_ref:.6f} ({ref_db:.3f}dB), diff={ddb:+.4f}dB (tol {tol_db})"
    record(name, ok, detail)
    return ddb


def _old_buggy_run_transfer_point(deck_lines, source, output_expr, f, amplitude, **kw):
    """Reproduces tools/fra.py's pre-card-#485 run_transfer_point: VO hard-coded to
    0 in the SIN literal (the exact bug), everything else identical. Used only for
    the red-proof below -- not a copy kept as production code anywhere.
    """
    settle_time = kw.get("settle_time", 0.0)
    settle_periods = kw.get("settle_periods", 20)
    measure_periods = kw.get("measure_periods", 10)
    points_per_period = kw.get("points_per_period", 256)
    outdir = kw.get("outdir") or os.path.join(fra.DATA_ROOT, "transfer_point_oldbug")
    ngspice = fra.NGSPICE
    os.makedirs(outdir, exist_ok=True)
    td, dt, t0, twin, tstop = fra.window_times(settle_time, f, settle_periods, measure_periods, points_per_period)
    base_f = os.path.join(outdir, "base.txt")
    pert_f = os.path.join(outdir, "pert.txt")
    for fn in (base_f, pert_f):
        if os.path.exists(fn):
            os.remove(fn)
    off = fra._sin_spec(0.0, f, td)          # the bug: no vo=
    on = fra._sin_spec(amplitude, f, td)     # the bug: no vo=
    control = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
op
alter @{source}[sin] = {off}
tran {dt:.12e} {tstop:.12e} 0 {dt:.12e}
linearize {output_expr}
wrdata {base_f} {output_expr}
alter @{source}[sin] = {on}
tran {dt:.12e} {tstop:.12e} 0 {dt:.12e}
linearize {output_expr}
wrdata {pert_f} {output_expr}
quit
.endc
.end
"""
    fra._run_ngspice(deck_lines, control, outdir, "fra_transfer_oldbug", ngspice)
    tb, (yb,) = fra.read_tran(base_f, 1)
    tp, (yp,) = fra.read_tran(pert_f, 1)
    fra.check_uniform(tb, dt)
    fra.check_uniform(tp, dt)
    dy = fra.subtract(tp, yp, tb, yb, output_expr)
    u0 = fra.injected_reference_phasor(f, td)
    h = fra.extract_component(tp, dy, f, t0, twin) / (amplitude * u0)
    return {"f": f, "H": h, "window": {"td": td, "dt": dt, "t0": t0, "twin": twin, "tstop": tstop}}


def main():
    t_start = time.time()
    os.makedirs(OUT, exist_ok=True)

    # ---- fixed tool: matches the analytic sech^2(Vdc) gain at three bias points ----
    for vdc in (0.0, 0.5, 1.2):
        check_gain_at(f"fixed tool: tanh gain at Vdc={vdc} matches sech^2(Vdc)",
                      vdc, fra.run_transfer_point)

    # ---- red-proof: the OLD (VO=0) code must FAIL at a nonzero bias ----
    for vdc in (0.5, 1.2):
        deck_path = _deck_with_vdc(vdc, os.path.join(OUT, f"redproof_vdc_{vdc}.cir"))
        deck_lines = stb.prepare_deck(deck_path)
        try:
            r = _old_buggy_run_transfer_point(deck_lines, "vin", "v(out)", FREQ, AMPLITUDE,
                                               settle_time=0.0, settle_periods=2, measure_periods=5,
                                               points_per_period=64,
                                               outdir=os.path.join(OUT, f"redproof_{vdc}"))
            h_ref = tanh_gain(vdc)
            got_db, ref_db = stbcore.db(r["H"]), stbcore.db(h_ref)
            ddb = got_db - ref_db
            # red-proof: the bug replaces Vin with a SIN whose baseline is 0V, so the
            # OLD code should measure tanh's slope AT Vdc=0 (gain 1.0, 0dB) regardless
            # of the requested vdc -- i.e. it should disagree with the true sech^2(vdc)
            # answer by roughly ref_db itself once vdc is not near 0.
            record(f"red-proof: OLD (VO=0) tool fails at Vdc={vdc} (measures the Vdc=0 point instead)",
                   abs(ddb) > DB_TOL, f"got {got_db:.3f}dB (~matches sech^2(0)=0dB), "
                                       f"true sech^2({vdc})={ref_db:.3f}dB, diff={ddb:+.3f}dB "
                                       f"(must exceed {DB_TOL}dB tol)")
        except stb.CouldNotRun as e:
            record(f"red-proof: OLD (VO=0) tool fails at Vdc={vdc}", None, str(e))

    # ---- refusal: a source that already carries a time function ----
    deck_lines = stb.prepare_deck(os.path.join(ROOT, "refs", "r3_vm_buck.cir"))
    try:
        fra.source_dc_value(deck_lines, "vramp")
        record("source_dc_value refuses an existing PULSE source (vramp)", False, "did not raise")
    except stb.CouldNotRun as e:
        record("source_dc_value refuses an existing PULSE source (vramp)", True, str(e))

    n_bad = sum(1 for _, s, _ in results if s == "FAIL")
    n_cnr = sum(1 for _, s, _ in results if s == "COULD-NOT-RUN")
    print(f"\n{len(results) - n_bad - n_cnr}/{len(results)} pass"
          + (f", {n_bad} FAIL" if n_bad else "") + (f", {n_cnr} COULD-NOT-RUN" if n_cnr else "")
          + f"  ({time.time() - t_start:.1f}s)")
    if n_bad:
        return 1
    if n_cnr:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
