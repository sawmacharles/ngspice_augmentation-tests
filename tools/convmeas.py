#!/usr/bin/env python3
"""convmeas.py -- switching-converter measurement helpers (card #497,
feature-survey #16, research/feature-survey.md): settled ripple pk-pk,
average efficiency, settling time to a +-band, overshoot/undershoot on a
step, switching frequency, and duty cycle.

Two layers, same split as tools/stb.py / tools/stbcore.py:
  - pure functions (ripple_pp, average_power, efficiency, switching_frequency,
    duty_cycle_threshold, settling_time, overshoot_undershoot) -- stdlib
    only, operate on plain (t, y) sample lists already in memory. Import
    these directly if you already have a trace loaded some other way.
  - a thin CLI driver that runs ONE `.tran` (system ngspice is enough --
    nothing here needs the private fork) via the documented `.control`
    snippet in docs/convmeas.md, reads the resulting wrdata trace, and
    prints/writes every requested measurement.

    tools/convmeas.py DECK --vin Vin --vout-node vout --rload-ohms 3.3 \\
        --period 2e-6 --periods 5 --tstop 1.5e-3 [--pwm-node pwm] \\
        [--target 3.3 --band-pct 1] [--load-step-time 1e-3] [--out DIR]

Input current convention (matches tools/loss.py's p_in = -avg(v*i)): the
native MNA branch current i(Vxxx) of a SPICE voltage source is defined
flowing INTO its + terminal, so the power the source actually DELIVERS to
the rest of the circuit is -v(t)*i(t), not +v(t)*i(t) -- efficiency() here
negates it the same way loss.py does, for the same reason.

Load current is computed as v(vout_node)/rload_ohms (Ohm's law), not read
back as a native i(Rload) vector -- ngspice does not track a resistor's own
branch current as an ordinary circuit vector any more than it does an
`@device[i]` OP parameter without an explicit `save` (tools/loss.py's own
"WHY @device[i]/@device[p] LOOKED FROZEN" writeup, docs/loss.md); computing
it from the known resistor value sidesteps that trap entirely and needs no
`save` line in the deck at all.

Exit status: 0 = ran and computed every requested measurement; 2 = could
not run (ngspice error, missing trace, a window with too few samples, or a
malformed --load-step-time with no step actually crossed).
Outputs (under ~/Work/data/ngspice_augmentation/convmeas/, never the
repository): wrapper.cir, ngspice.log, trace.txt (raw wrdata), summary.json.
"""
import argparse
import csv
import datetime
import json
import math
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stb  # noqa: E402 -- reuses prepare_deck() (include-path fixup, drop
             # existing analyses/control blocks), same convention as stb_tran.py

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/convmeas")
CouldNotRun = stb.CouldNotRun


# ---------------------------------------------------------------------------
# pure measurement functions -- stdlib only, operate on parallel (t, y) lists
# ---------------------------------------------------------------------------

def _window(t, y, t0, t1):
    return [(ti, yi) for ti, yi in zip(t, y) if t0 <= ti <= t1]


def ripple_pp(t, y, t0, t1):
    """max(y) - min(y) over [t0, t1] -- settled ripple peak-to-peak."""
    pts = _window(t, y, t0, t1)
    if not pts:
        raise CouldNotRun(f"ripple_pp: no samples in [{t0}, {t1}]")
    ys = [yi for _, yi in pts]
    return max(ys) - min(ys)


def average_value(t, y, t0, t1):
    """Trapezoidal time-average of y over [t0, t1]. Exact for the AVG-over-a-
    window trick refs/r3_vm_buck.cir's own header documents (average of a
    saturating 0/1 PWM signal over exactly one period = duty cycle, whatever
    the transition shape) -- use directly for that; duty_cycle_threshold()
    below is the more general (any-shaped-signal, explicit threshold)
    alternative."""
    pts = _window(t, y, t0, t1)
    if len(pts) < 2:
        raise CouldNotRun(f"average_value: fewer than 2 samples in [{t0}, {t1}]")
    total = 0.0
    for k in range(len(pts) - 1):
        ta, ya = pts[k]
        tb, yb = pts[k + 1]
        total += 0.5 * (ya + yb) * (tb - ta)
    span = pts[-1][0] - pts[0][0]
    if span <= 0:
        raise CouldNotRun(f"average_value: degenerate window [{t0}, {t1}]")
    return total / span


def average_power(t, v, i, t0, t1):
    """Trapezoidal average of v(t)*i(t) over [t0, t1] -- v and i must already
    share t's grid (one wrdata call, same convention as tools/loss.py)."""
    prod = [vi * ii for vi, ii in zip(v, i)]
    return average_value(t, prod, t0, t1)


def efficiency(t, vin, iin, vout, iout, t0, t1):
    """Pout/Pin over [t0, t1]. iin is the native i(Vxxx) branch current of the
    input source (SPICE sign convention: flows INTO its + terminal) -- see
    module docstring for the -avg(v*i) sign, matching tools/loss.py."""
    p_in = -average_power(t, vin, iin, t0, t1)
    p_out = average_power(t, vout, iout, t0, t1)
    return p_out / p_in if p_in else float("nan")


def _crossings(t, y, threshold, rising=True):
    """Linearly-interpolated times where y crosses threshold, one direction."""
    out = []
    for k in range(len(t) - 1):
        y0, y1 = y[k], y[k + 1]
        hit = (y0 < threshold <= y1) if rising else (y0 >= threshold > y1)
        if hit and y1 != y0:
            frac = (threshold - y0) / (y1 - y0)
            out.append(t[k] + frac * (t[k + 1] - t[k]))
    return out


def switching_frequency(t, y, t0, t1, threshold=0.5):
    """1/period, from the average spacing of RISING threshold crossings within
    [t0, t1] (first crossing to last, divided by crossing count - 1 -- more
    accurate than any single period, same idea as .meas TRIG/TARG RISE=k over
    many periods). None (never a fabricated number) if fewer than 2 crossings
    fall in the window."""
    xs = [c for c in _crossings(t, y, threshold, rising=True) if t0 <= c <= t1]
    if len(xs) < 2:
        return None
    return (len(xs) - 1) / (xs[-1] - xs[0])


def duty_cycle_threshold(t, y, t0, t1, threshold=0.5):
    """Fraction of time y(t) >= threshold within [t0, t1], via linear
    interpolation between samples for the crossing instant within each
    interval (exact for a piecewise-linear y, i.e. the sampled trace taken
    literally). General alternative to average_value() for a signal that
    is not close to a clean 0/1 square wave."""
    pts = _window(t, y, t0, t1)
    if len(pts) < 2:
        raise CouldNotRun(f"duty_cycle_threshold: fewer than 2 samples in [{t0}, {t1}]")
    total = high = 0.0
    for k in range(len(pts) - 1):
        ta, ya = pts[k]
        tb, yb = pts[k + 1]
        dt = tb - ta
        total += dt
        if ya >= threshold and yb >= threshold:
            high += dt
        elif ya < threshold and yb < threshold:
            pass
        else:
            frac = (threshold - ya) / (yb - ya)
            tc = ta + frac * dt
            high += (tc - ta) if ya >= threshold else (tb - tc)
    return high / total if total > 0 else float("nan")


def settling_time(t, y, target, band, t_start=0.0):
    """First time after which y(t) stays within [target-band, target+band] for
    the rest of the given data ("last time it left the band", the standard,
    robust settling-time definition -- NOT "first time it entered the band",
    which a signal can cross through on its way to overshooting). Returns
    t_start if y is already within band for the whole [t_start, end] range;
    returns None (never a fabricated time) if y is still outside the band at
    the very last sample -- "has not settled within the data given", not
    silently reported as the last sample time."""
    pts = [(ti, yi) for ti, yi in zip(t, y) if ti >= t_start]
    if not pts:
        raise CouldNotRun(f"settling_time: no samples at/after t_start={t_start}")
    last_outside = None
    for ti, yi in pts:
        if abs(yi - target) > band:
            last_outside = ti
    if last_outside is None:
        return t_start
    if last_outside >= pts[-1][0]:
        return None
    return last_outside


def overshoot_undershoot(t, y, t_event, y_before, y_after):
    """Peak excursion beyond y_after, as a % of the step size |y_after -
    y_before|, for t >= t_event. "overshoot" = the FIRST, and largest,
    excursion in the SAME direction as the step (the classic control-theory
    definition -- e.g. a rising step overshooting above its new final
    value); "undershoot" = the excursion in the OPPOSITE direction that
    follows it (a dip below the new level on a rising step, or a bump above
    it on a falling one -- what an underdamped step response's second
    half-cycle produces before it finally settles). Undershoot is searched
    only from the overshoot peak's own time onward, not from t_event --
    searching the whole post-event range for "the opposite-direction
    extremum" would instead just find the signal's still-near-y_before
    starting point immediately after t_event on any transient that takes
    time to move at all, which is not what "undershoot" means. Both
    percentages are >= 0 (clamped) and 0.0 exactly if the trace never
    crosses y_after in that direction at all."""
    after = [(ti, yi) for ti, yi in zip(t, y) if ti >= t_event]
    if not after:
        raise CouldNotRun(f"overshoot_undershoot: no samples at/after t_event={t_event}")
    step = y_after - y_before
    if step > 0:
        t_peak, peak_max = max(after, key=lambda p: p[1])
        rest = [yi for ti, yi in after if ti >= t_peak]
        peak_min = min(rest) if rest else peak_max
        overshoot_pct = max(0.0, (peak_max - y_after) / step * 100.0)
        undershoot_pct = max(0.0, (y_after - peak_min) / step * 100.0)
    elif step < 0:
        t_peak, peak_min = min(after, key=lambda p: p[1])
        rest = [yi for ti, yi in after if ti >= t_peak]
        peak_max = max(rest) if rest else peak_min
        overshoot_pct = max(0.0, (y_after - peak_min) / (-step) * 100.0)
        undershoot_pct = max(0.0, (peak_max - y_after) / (-step) * 100.0)
    else:
        peak_max = max(yi for _, yi in after)
        peak_min = min(yi for _, yi in after)
        overshoot_pct = undershoot_pct = float("nan")
    return {"overshoot_pct": overshoot_pct, "undershoot_pct": undershoot_pct,
            "peak_max": peak_max, "peak_min": peak_min}


# ---------------------------------------------------------------------------
# CLI driver: one .tran, one wrdata trace, every requested measurement
# ---------------------------------------------------------------------------

def read_wrdata_real(path, ncols):
    """time + ncols real columns (set wr_singlescale, header line present --
    the same convention tools/stb.py uses for its complex traces, adapted
    for plain real signals; see docs/convmeas.md's documented .control
    snippet)."""
    if not os.path.exists(path):
        raise CouldNotRun(f"ngspice wrote no {path}")
    t, cols = [], [[] for _ in range(ncols)]
    with open(path) as fh:
        next(fh, None)
        for row in fh:
            v = row.split()
            if not v:
                continue
            if len(v) != 1 + ncols:
                raise CouldNotRun(f"unexpected column count in {path}: {row!r}")
            v = [float(x) for x in v]
            t.append(v[0])
            for k in range(ncols):
                cols[k].append(v[1 + k])
    if not t:
        raise CouldNotRun(f"{path} is empty")
    return t, cols


def run_tran(deck, vin, vout_node, rload_ohms, pwm_node, tstep, tstop, outdir, ngspice, timeout=1800):
    os.makedirs(outdir, exist_ok=True)
    lines = stb.prepare_deck(deck)
    vin_node = None
    for line in lines:
        f = line.split()
        if f and f[0].lower() == vin.lower():
            vin_node = f[1]
            break
    if vin_node is None:
        raise CouldNotRun(f"input source {vin!r} not found in deck")
    trace = os.path.join(outdir, "trace.txt")
    if os.path.exists(trace):
        os.remove(trace)
    sigs = [f"v({vin_node})", f"i({vin})", f"v({vout_node})"]
    if pwm_node:
        sigs.append(f"v({pwm_node})")
    body = "\n".join(lines)
    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
tran {tstep} {tstop} 0 uic
wrdata {trace} {' '.join(sigs)}
quit
.endc
.end
"""
    wrapper = os.path.join(outdir, "wrapper.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + ctl)
    try:
        proc = subprocess.run([ngspice, "-b", wrapper], capture_output=True, text=True, timeout=timeout)
    except (FileNotFoundError, PermissionError) as e:
        raise CouldNotRun(f"cannot start {ngspice}: {e}") from e
    except subprocess.TimeoutExpired as e:
        raise CouldNotRun(f"tran timed out after {timeout}s") from e
    with open(os.path.join(outdir, "ngspice.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")
    ncols = 3 + (1 if pwm_node else 0)
    t, cols = read_wrdata_real(trace, ncols)
    result = {"t": t, "vin": cols[0], "iin": cols[1], "vout": cols[2]}
    if pwm_node:
        result["pwm"] = cols[3]
    result["iout"] = [v / rload_ohms for v in cols[2]]
    return result


def compute(trace, period, periods, tstop, target=None, band_pct=1.0,
            load_step_time=None, pre_window=None, post_window=None, threshold=0.5):
    t = trace["t"]
    window_t1 = tstop
    window_t0 = tstop - periods * period
    out = {"window": [window_t0, window_t1], "period": period, "periods": periods}

    out["ripple_pp_v"] = ripple_pp(t, trace["vout"], window_t0, window_t1)
    out["efficiency"] = efficiency(t, trace["vin"], trace["iin"], trace["vout"], trace["iout"],
                                    window_t0, window_t1)
    out["p_in_w"] = -average_power(t, trace["vin"], trace["iin"], window_t0, window_t1)
    out["p_out_w"] = average_power(t, trace["vout"], trace["iout"], window_t0, window_t1)

    if "pwm" in trace:
        out["switching_frequency_hz"] = switching_frequency(t, trace["pwm"], window_t0, window_t1, threshold)
        out["duty_cycle_avg"] = average_value(t, trace["pwm"], window_t0, window_t0 + period)
        out["duty_cycle_threshold"] = duty_cycle_threshold(t, trace["pwm"], window_t0, window_t0 + period, threshold)

    if target is not None:
        band = band_pct / 100.0 * abs(target)
        out["settling_time_s"] = settling_time(t, trace["vout"], target, band)
        out["settling_band_v"] = band

    if load_step_time is not None:
        pw0, pw1 = pre_window or (load_step_time - period, load_step_time)
        qw0, qw1 = post_window or (window_t0, window_t1)
        y_before = average_value(t, trace["vout"], pw0, pw1)
        y_after = average_value(t, trace["vout"], qw0, qw1)
        ou = overshoot_undershoot(t, trace["vout"], load_step_time, y_before, y_after)
        out["load_step"] = {"t_event": load_step_time, "y_before": y_before, "y_after": y_after, **ou}

    return out


def describe(res):
    lines = [f"window: [{res['window'][0]:.6g}, {res['window'][1]:.6g}] s ({res['periods']} periods of {res['period']:.6g} s)",
             f"ripple pk-pk: {res['ripple_pp_v']*1e3:.4f} mV" if abs(res["ripple_pp_v"]) < 1
             else f"ripple pk-pk: {res['ripple_pp_v']:.6g} V",
             f"efficiency: {100*res['efficiency']:.4f} % (Pin={res['p_in_w']:.6g} W, Pout={res['p_out_w']:.6g} W)"]
    if "switching_frequency_hz" in res:
        fsw = res["switching_frequency_hz"]
        lines.append(f"switching frequency: {fsw:.6g} Hz" if fsw is not None else "switching frequency: no 2nd crossing (could not measure)")
        lines.append(f"duty cycle (avg): {res['duty_cycle_avg']:.4f}   (threshold): {res['duty_cycle_threshold']:.4f}")
    if "settling_time_s" in res:
        st = res["settling_time_s"]
        lines.append(f"settling time (+-{100*res['settling_band_v']/abs(res.get('target', res['settling_band_v']) or 1):.3g}% band): "
                      f"{st:.6g} s" if st is not None else "settling time: NOT settled within the given data")
    if "load_step" in res:
        ls = res["load_step"]
        lines.append(f"load step at t={ls['t_event']:.6g}s: {ls['y_before']:.6g} -> {ls['y_after']:.6g} V; "
                      f"overshoot {ls['overshoot_pct']:.3f}%, undershoot {ls['undershoot_pct']:.3f}%")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("--vin", required=True, help="input voltage source instance name")
    ap.add_argument("--vout-node", required=True, help="output node name (no v(...) wrapper)")
    ap.add_argument("--rload-ohms", type=float, required=True, help="load resistor value, Ohms (iout = vout/rload)")
    ap.add_argument("--pwm-node", help="switch-drive/PWM node for switching-frequency + duty-cycle measurement")
    ap.add_argument("--period", type=float, required=True, help="switching period, seconds")
    ap.add_argument("--periods", type=float, default=5, help="number of trailing periods to measure over (default 5)")
    ap.add_argument("--tstop", type=float, required=True, help="total .tran run time, seconds")
    ap.add_argument("--tstep", type=float, default=None, help="fixed max .tran timestep (default period/2000)")
    ap.add_argument("--threshold", type=float, default=0.5, help="pwm threshold for switching-freq/duty-cycle (default 0.5)")
    ap.add_argument("--target", type=float, help="settling-time target value (vout)")
    ap.add_argument("--band-pct", type=float, default=1.0, help="settling-time band, %% of --target (default 1.0)")
    ap.add_argument("--load-step-time", type=float, help="overshoot/undershoot: time of a load/reference step on vout")
    ap.add_argument("--pre-window", nargs=2, type=float, metavar=("T0", "T1"),
                     help="window to average vout over for the PRE-step level (default: [step-period, step])")
    ap.add_argument("--post-window", nargs=2, type=float, metavar=("T0", "T1"),
                     help="window to average vout over for the POST-step level (default: the measurement window)")
    ap.add_argument("--out", help="output directory (default: under ~/Work/data/ngspice_augmentation/convmeas)")
    ap.add_argument("--ngspice", default=NGSPICE, help="ngspice binary (system ngspice is enough)")
    ap.add_argument("--timeout", type=float, default=1800)
    args = ap.parse_args()

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = args.out or os.path.join(DATA_ROOT, f"{os.path.splitext(os.path.basename(args.deck))[0]}-{stamp}")
    tstep = args.tstep or args.period / 2000.0

    try:
        trace = run_tran(args.deck, args.vin, args.vout_node, args.rload_ohms, args.pwm_node,
                          tstep, args.tstop, outdir, args.ngspice, timeout=args.timeout)
        res = compute(trace, args.period, args.periods, args.tstop,
                       target=args.target, band_pct=args.band_pct,
                       load_step_time=args.load_step_time,
                       pre_window=tuple(args.pre_window) if args.pre_window else None,
                       post_window=tuple(args.post_window) if args.post_window else None,
                       threshold=args.threshold)
    except (CouldNotRun, subprocess.TimeoutExpired) as e:
        print(f"COULD NOT RUN: {e}", file=sys.stderr)
        return 2

    if args.target is not None:
        res["target"] = args.target
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print(describe(res))
    print(f"results: {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
