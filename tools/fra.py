#!/usr/bin/env python3
"""fra.py — transient frequency-response reference harness (FRA) on unmodified ngspice-45.2.

Card #452 (plan T2.3), 2026-09-27. This is the independent reference that later phases
(`pac_sb0`, `pstb_loop`, `pac_sampling`) are checked against (verify.md, "Master
acceptance loop"). It never uses a periodic-AC or shooting analysis itself: everything is
plain `.tran`, the PLECS/SIMPLIS-FX "perturbed minus unperturbed" trick (research/
simplis-pop-ac.md Stage 2') applied to Tian's two-injection loop probe (research/
cadence-stb-pstb.md §2.4, path (i) "brute-force transient version").

Reuses tools/stb.py (prepare_deck, probe_nodes, CouldNotRun, NGSPICE, DATA_ROOT) and
tools/stbcore.py (tian, margins, db, unwrapped_phase_deg, max_rel_err) rather than
reimplementing them; those two files and tests/test_stb.py belong to another worker on
this project and are not touched here.

Design choices (see the card #452 report for the numbers behind them):

- **Same starting state, no savestate/loadstate.** ngspice-45.2 exposes no savestate/
  loadstate this harness could find. Instead every run is a deterministic replay from
  t=0 through the same `op`, and the probe's Vi/Ii sources carry a SIN(0 A f TD 0 0)
  waveform: value 0 (a dead short/open) until t=TD, then a clean sine. Because SPICE is
  deterministic and Vi/Ii are 0 before TD in every run, the unperturbed, voltage-injection
  and current-injection runs are bit-for-bit identical up to t=TD (verified empirically,
  card #452 notes) — the "same state" requirement is satisfied by construction, not by
  hoping three separate processes stay in sync.
- **Non-uniform ngspice time steps.** ngspice's adaptive transient stepper does not
  sample uniformly even with a small `tstep` (confirmed empirically: steps as fine as
  5e-8 s next to steps at the nominal 5e-6 s tstep). Every run therefore also passes
  `tmax = tstep` to `tran` and applies the `linearize` control command to the vectors
  before `wrdata`, which resamples them onto an exact uniform grid at multiples of
  tstep — the DFT below assumes that grid is uniform and checks it (`check_uniform`).
- **Sideband-0 extraction.** `extract_component` is a trapezoidal-rule DFT of
  y(t)*e^{-j*2*pi*f*t} over the requested window; for a window that is an exact integer
  number of excitation periods this is the exact (spectrally accurate) sideband-0
  Fourier coefficient of a signal already in sinusoidal steady state, times 2 (peak, not
  RMS). The window is [TD + settle_periods/f, TD + (settle_periods+measure_periods)/f).
- **Phase reference / unit-injection normalisation.** The injected SIN(...) has an
  implicit phase (it starts at value 0, not at its peak) and TD shifts it further, so
  raw extracted phasors are not directly comparable to a `.ac` result with acmag=1
  acphase=0. `injected_reference_phasor` gives the closed-form phasor of sin(w(t-TD))
  under the same extraction convention (exact when the window is an integer number of
  periods, by construction: -1j*exp(-j*w*TD)); every extracted V/I is divided by
  amplitude*that phasor before Tian's formula is applied, which reproduces the
  acmag=1/acphase=0 convention `tools/stb.py` and `.ac` use.
- **Ill-conditioning at high loop gain.** T = S/(1-S) (stbcore.tian) is only as accurate
  as (1-S) is far from 0. Deep in a high-gain passband S -> 1 and a fraction-of-a-percent
  numerical error in the extracted V1 becomes a large error in T; this is intrinsic to
  the two-injection method (a real network analyzer has the same problem), not a bug
  here — see the card #452 report and tests/test_fra.py's frequency choices. Mitigated
  (not eliminated) by also extracting v(a) (v-injection run) and i(v.<probe>.vb)
  (i-injection run) and computing stbcore.middlebrook from the same three transients;
  `combine_tian_middlebrook` reports Tian below a threshold (default 25 dB, the spec's
  20-30 dB ill-conditioning band) and Middlebrook at/above it, plus the disagreement
  between the two at every point, so the two methods' agreement is visible and not just
  asserted.
- **Parallel dispatch.** One ngspice process per frequency (3 sequential `tran` runs
  inside it, sharing one `op`), run under a ProcessPoolExecutor capped at
  max(1, os.cpu_count()//2) workers, each in its own directory under
  ~/Work/data/ngspice_augmentation/fra/... — never in the repository.
- **Adaptive settle (card #487, opt-in via `--adaptive-settle`).** The fixed `--settle`
  above is a human guess at how long the *unforced* circuit needs to reach its own
  periodic steady state before injection turns on; card #484/#487 found this guess too
  short for slow converters (R3's 300us guess left a 6.15% periodicity residual, biasing
  the FRA curve by up to 3.7dB/17deg near fs/2 — R3 actually needs ~9ms to reach <=1%).
  `periodicity_scan` runs the *unperturbed* circuit ONCE from t=0 to a cap, samples the
  same nodes `run_loop_point`/`run_transfer_point` extract from, and finds the earliest
  multiple of the caller-given natural period `T` at which
  `max|v(t)-v(t-T)|/ripple(t-T..t)` (ripple = that window's own peak-to-peak — same
  convention `tools/pss_ref.py`'s `ripple()`/`settle_period_index()` and
  `tools/pssac_ref.py`'s `_measure_h_tran` periodicity check already use elsewhere in
  this project) is at or below `--residual-bound` (default 1%); `sweep_loop_adaptive`/
  `sweep_transfer_adaptive` then reuse that ONE measured settle time for every frequency
  point, exactly like a fixed `--settle` value, so the 3N/N re-settles per sweep are
  unchanged in count — only `td` is measured instead of guessed. Combine with
  `--checkpoint` (mode loop only) to also pay that settle cost once total instead of
  once per sub-transient/frequency (docs/checkpoint.md). The old fixed-settle path
  (`--settle`, no `--adaptive-settle`) is completely unchanged and is still every
  existing test's default.

    tools/fra.py DECK --probe XP --freqs "dec 20 1e5 1e8" --amplitude 1e-3 [--out DIR]
    tools/fra.py DECK --mode transfer --source Vin --output "v(out)" --freqs ... --amplitude ...
    tools/fra.py DECK --probe XP --freqs ... --amplitude ... --adaptive-settle \\
        --settle-period 2e-6 --residual-bound 0.01 --settle-cap 0.05 [--checkpoint]

Exit status: 0 = ran (some points may still be individually could-not-run, see the CSV),
2 = could not run at all (e.g. ngspice missing, deck/probe not found).
"""
import argparse
import bisect
import cmath
import concurrent.futures as cf
import datetime
import json
import math
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stb  # noqa: E402
import stbcore  # noqa: E402

CouldNotRun = stb.CouldNotRun
NGSPICE = stb.NGSPICE
DATA_ROOT = os.path.join(stb.DATA_ROOT, "fra")


# ---------------------------------------------------------------------------
# windowing and the DFT/correlation extraction
# ---------------------------------------------------------------------------

def dt_for(f, points_per_period):
    """Uniform sample spacing for `points_per_period` samples of one cycle at f."""
    return 1.0 / (f * points_per_period)


def window_times(settle_time, f, settle_periods, measure_periods, points_per_period):
    """(td, dt, t0, twin, tstop) for one frequency point.

    td: absolute time the injection SIN turns on (settle_time, rounded up to a multiple
        of dt so it lands on the uniform grid).
    dt: uniform sample spacing (also used as tran's tstep and tmax).
    t0: measurement window start = td + settle_periods excitation periods (lets the
        circuit's own response to the sudden sine ring out).
    twin: measurement window length = measure_periods excitation periods. May be
        deliberately non-integer (red-proof); the window is then not spectrally exact.
    tstop: end of the simulated transient.
    """
    dt = dt_for(f, points_per_period)
    td = math.ceil(settle_time / dt) * dt if settle_time > 0 else 0.0
    period = 1.0 / f
    t0 = td + settle_periods * period
    twin = measure_periods * period
    return td, dt, t0, twin, t0 + twin


def check_uniform(times, dt, tol=0.02):
    """Raise if the (post-`linearize`) time vector is not uniform at spacing dt."""
    if len(times) < 3:
        raise CouldNotRun(f"time series has only {len(times)} points")
    bad = 0
    worst = 0.0
    for t1, t2 in zip(times, times[1:]):
        d = t2 - t1
        err = abs(d - dt)
        if err > tol * dt:
            bad += 1
            worst = max(worst, err)
    if bad:
        raise CouldNotRun(
            f"non-uniform time grid after linearize: {bad}/{len(times) - 1} steps off "
            f"target dt={dt:.6g} by more than {tol * 100:.0f}% (worst {worst:.3g})"
        )


def extract_component(times, values, f, t0, twin):
    """Complex sideband-0 phasor of `values` at frequency f over [t0, t0+twin].

    Trapezoidal-rule DFT: (2/twin) * integral of y(t) * exp(-j*2*pi*f*t) dt, referenced
    to absolute t=0. Exact (to quadrature precision) for a signal in sinusoidal steady
    state when twin is an integer number of periods of f; otherwise this leaks
    (deliberately, for the red-proof window check).
    """
    if twin <= 0:
        raise CouldNotRun(f"non-positive measurement window ({twin:.6g} s)")
    omega = 2 * math.pi * f
    eps = max(1e-15, 1e-12 * twin)
    pts = [(t, y) for t, y in zip(times, values) if t0 - eps <= t <= t0 + twin + eps]
    if len(pts) < 3:
        raise CouldNotRun(
            f"measurement window [{t0:.6g}, {t0 + twin:.6g}] s has only {len(pts)} "
            f"samples; widen the window or shrink points_per_period/dt"
        )
    acc = 0j
    for (t1, y1), (t2, y2) in zip(pts, pts[1:]):
        acc += 0.5 * (y1 * cmath.exp(-1j * omega * t1) + y2 * cmath.exp(-1j * omega * t2)) * (t2 - t1)
    return 2 * acc / twin


def injected_reference_phasor(f, td):
    """Closed-form phasor of sin(2*pi*f*(t-td)), same convention as extract_component,
    exact when the measurement window is an integer number of periods after td:
        integral of sin(w(t-td)) * exp(-jwt) dt over an integer number of periods
        = -j * twin/2 * exp(-j*w*td), so (2/twin)*that = -j*exp(-j*w*td).
    """
    return -1j * cmath.exp(-1j * 2 * math.pi * f * td)


def subtract(times_a, ys_a, times_b, ys_b, what):
    """ys_a - ys_b, after checking the two runs share the same time grid."""
    if len(times_a) != len(times_b):
        raise CouldNotRun(f"{what}: perturbed/unperturbed runs have different lengths")
    for i in (0, len(times_a) // 2, -1):
        if abs(times_a[i] - times_b[i]) > 1e-9 * max(1.0, abs(times_a[i])):
            raise CouldNotRun(f"{what}: perturbed/unperturbed runs are not on the same time grid")
    return [a - b for a, b in zip(ys_a, ys_b)]


# ---------------------------------------------------------------------------
# adaptive settle (card #487): periodicity residual + a single settle-time scan,
# reused by both mode "loop" and mode "transfer" instead of a fixed --settle guess.
# ---------------------------------------------------------------------------

DEFAULT_RESIDUAL_BOUND = 0.01          # 1% of the window's own ripple, card #487's default
DEFAULT_SETTLE_CAP = 0.05              # 50ms; comfortably above R3's measured ~9ms need
DEFAULT_SETTLE_POINTS_PER_PERIOD = 64  # samples/period for the settle-scan trace itself
DEFAULT_SETTLE_MIN_PERIODS = 2         # need at least 2 periods of data to compare 1 vs the next
DEFAULT_SETTLE_N_PHASES = 64           # phase samples/period used BY the residual comparison


def _interp_linear(times, values, t):
    """values(t), linearly interpolated; clamped (not extrapolated) at the ends."""
    if t <= times[0]:
        return values[0]
    if t >= times[-1]:
        return values[-1]
    i = bisect.bisect_right(times, t) - 1
    if i >= len(times) - 1:
        return values[-1]
    t0, t1 = times[i], times[i + 1]
    v0, v1 = values[i], values[i + 1]
    return v0 if t1 == t0 else v0 + (t - t0) * (v1 - v0) / (t1 - t0)


def periodicity_residual(times, values, period, t_end, n_phases=DEFAULT_SETTLE_N_PHASES):
    """max|v(t)-v(t-period)| / ripple, `n_phases` phase samples over the one-period window
    ending at t_end, ripple = that window's own peak-to-peak (same convention as
    tools/pss_ref.py's ripple()/settle_period_index() and tools/pssac_ref.py's
    _measure_h_tran periodicity check — reimplemented here, self-contained, rather than
    imported, to keep fra.py's existing ownership boundary with tools/pss_ref.py/
    tools/pssac_ref.py, card #452's module docstring). Falls back to the plain (non-
    relative) residual if the window's ripple is numerically zero (a flat signal — never
    divide by ~0)."""
    if t_end - 2 * period < times[0] - 1e-15:
        raise CouldNotRun(
            f"periodicity_residual: window [{t_end - 2 * period:.6g}, {t_end:.6g}] starts "
            f"before the available trace ({times[0]:.6g}s)")
    w2 = [_interp_linear(times, values, t_end - period + k * period / n_phases) for k in range(n_phases)]
    w1 = [_interp_linear(times, values, t_end - 2 * period + k * period / n_phases) for k in range(n_phases)]
    ripple = max(w2) - min(w2)
    resid = max(abs(b - a) for a, b in zip(w1, w2))
    return resid / ripple if ripple > 1e-12 else resid


def periodicity_scan(deck_lines, node_exprs, period, residual_bound=DEFAULT_RESIDUAL_BOUND,
                      cap=DEFAULT_SETTLE_CAP, points_per_period=DEFAULT_SETTLE_POINTS_PER_PERIOD,
                      min_periods=DEFAULT_SETTLE_MIN_PERIODS, n_phases=DEFAULT_SETTLE_N_PHASES,
                      outdir=None, ngspice=None, startup="uic", tmax=None, spice_options=""):
    """Find how long `deck_lines` (unperturbed — no injection sources touched) needs to
    settle before its own periodicity residual (max over `node_exprs`, see
    periodicity_residual) is at or below `residual_bound`, by running ONE plain transient
    from t=0 to `cap` and scanning it — not one run per candidate settle time.

    `node_exprs`: ngspice vector expressions to track, e.g. ["v(xp.x)", "v(a)"] (mode
    loop's own two extracted node voltages) or ["v(vout)"] (mode transfer's own
    `output_expr`) — "the measured nodes" card #487 asks the residual to be checked over.
    `period`: the circuit's OWN natural period (e.g. 1/fs for a switching converter) used
    for the periodicity check — NOT an injection period; adaptive settle happens before
    any injection frequency is chosen, and one settle time is reused by every frequency
    point afterward (see sweep_loop_adaptive/sweep_transfer_adaptive).

    Returns a dict: settle_time (the found td, or `cap` if the bound was never reached),
    residual (the periodicity residual actually achieved there), capped (bool — True
    means `cap` was reached without meeting `residual_bound`; the returned settle_time/
    residual are then the BEST seen, never fabricated as if the bound had been met),
    periods_scanned, trace_path, node_exprs.
    """
    ngspice = ngspice or NGSPICE
    outdir = outdir or os.path.join(DATA_ROOT, "settle_scan")
    os.makedirs(outdir, exist_ok=True)
    if period <= 0:
        raise CouldNotRun(f"periodicity_scan: period must be > 0, got {period!r}")
    if cap <= 0:
        raise CouldNotRun(f"periodicity_scan: cap must be > 0, got {cap!r}")
    if not node_exprs:
        raise CouldNotRun("periodicity_scan: node_exprs must be non-empty")
    dt = period / points_per_period
    step_max = tmax if tmax is not None else dt
    if startup not in ("op", "uic"):
        raise CouldNotRun(f"unknown startup {startup!r}, expected 'op' or 'uic'")
    op_line = "op" if startup == "op" else ""
    uic = " uic" if startup == "uic" else ""
    trace_f = os.path.join(outdir, "settle_scan.txt")
    if os.path.exists(trace_f):
        os.remove(trace_f)
    nodes = " ".join(node_exprs)
    control = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
{spice_options}
{op_line}
tran {dt:.12e} {cap:.12e} 0 {step_max:.12e}{uic}
linearize {nodes}
wrdata {trace_f} {nodes}
quit
.endc
.end
"""
    _run_ngspice(deck_lines, control, outdir, "fra_settle_scan", ngspice)
    times, cols = read_tran(trace_f, len(node_exprs))
    check_uniform(times, dt)
    t0 = times[0]
    n_periods_avail = int((times[-1] - t0) / period)
    min_needed = min_periods + 1
    if n_periods_avail < min_needed:
        raise CouldNotRun(
            f"periodicity_scan: cap {cap:.6g}s only covers {n_periods_avail} periods of "
            f"{period:.6g}s (need >= {min_needed}); raise --settle-cap")
    best_k, best_resid = min_periods, None
    for k in range(min_periods, n_periods_avail + 1):
        t_end = t0 + k * period
        resid = max(periodicity_residual(times, col, period, t_end, n_phases) for col in cols)
        if best_resid is None or resid < best_resid:
            best_k, best_resid = k, resid
        if resid <= residual_bound:
            return {"settle_time": t_end, "residual": resid, "capped": False,
                    "periods_scanned": k, "trace_path": trace_f, "node_exprs": list(node_exprs)}
    return {"settle_time": t0 + best_k * period, "residual": best_resid, "capped": True,
            "periods_scanned": best_k, "trace_path": trace_f, "node_exprs": list(node_exprs)}


# ---------------------------------------------------------------------------
# ngspice plumbing
# ---------------------------------------------------------------------------

def read_tran(path, ncols):
    """time list and ncols real columns from a `linearize`d wrdata output."""
    if not os.path.exists(path):
        raise CouldNotRun(f"ngspice wrote no {path}")
    times, cols = [], [[] for _ in range(ncols)]
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if len(parts) != 1 + ncols:
                raise CouldNotRun(f"unexpected column count in {path}: {line!r}")
            vals = [float(x) for x in parts]
            times.append(vals[0])
            for k in range(ncols):
                cols[k].append(vals[1 + k])
    if not times:
        raise CouldNotRun(f"{path} is empty")
    return times, cols


def _sin_spec(amplitude, f, td, vo=0.0):
    """A `[ VO VA FREQ TD THETA PHASE ]` literal for `alter @dev[sin] = ...`.

    `vo` (default 0.0, unchanged behaviour): the source's own DC offset, which
    `alter @dev[sin] = ...` must be given explicitly or it silently zeroes it
    (card #485: `run_transfer_point`'s "off"/"pert" alters were both hard-coding
    VO=0 here, so injecting on any source with a nonzero DC value -- e.g. Vin for
    a line-to-output measurement -- removed its bias entirely, not just added a
    small-signal tone on top of it; wrong by 30-60 dB, since the "base"/"pert"
    runs were both effectively measuring the response with Vin=0V, not Vin=
    <its real DC value>). run_loop_point's callers are unaffected and still pass
    no `vo` (0.0): stb/lgprobe.lib's Vi/Ii probe sources are always `DC 0` by
    construction (their own file, lines 15-16), so 0 was already correct there.

    `vo` may be a plain number (formatted `.12e`, the pre-card-#485 behaviour) or
    a string (source_dc_value's return: a spice-number token like "12" or a
    `{param}` expression, embedded verbatim -- ngspice's `alter` evaluates
    brace expressions per-token, and re-formatting an arbitrary expression string
    through `.12e` would just be wrong).
    """
    vo_text = f"{vo:.12e}" if isinstance(vo, (int, float)) else str(vo)
    return f"[ {vo_text} {amplitude:.12e} {f:.12e} {td:.12e} 0 0 ]"


def _run_ngspice(deck_lines, control, outdir, label, ngspice, cwd=None, allow_patterns=()):
    """cwd (default None, unchanged behaviour): run ngspice with this as its working
    directory. Added for the checkpoint mode (settle_checkpoint/run_loop_point_
    checkpoint), which needs `savestate`/`loadstate`'s filename argument to be a
    bare relative name rather than an absolute path under ~/Work/... -- ngspice's
    control-mode tokenizer lower-cases bare command arguments (confirmed: a
    "ckptB.bin" argument was written/read back as "ckptb.bin", card #472 T10.4),
    which corrupts any path containing upper-case letters (this project's own
    ~/Work/... convention). A relative filename plus an explicit `cwd` sidesteps
    it entirely; existing callers all pass absolute paths already and never set
    this, so they are unaffected.

    allow_patterns (default (), unchanged behaviour): extra known-benign regexes
    (case-insensitive) whose matching lines are excluded from the error scan below
    before deciding could-not-run. Added for the checkpoint mode's netlist-level
    '.tran' card (needed so loadstate's find_tran_job() has a job to read tstop
    from): a `.control` block with no `.plot`/`.print`/`.fourier` line makes
    ngspice print "Error: incomplete or empty netlist / or no '.plot', '.print', "
    "or '.fourier' lines in batch mode" even on an otherwise fully successful run
    (confirmed: the run's own wrdata output was present and correct when this line
    appeared, card #472 T10.4) -- same class of harmless-but-matches-"error" noise
    as R6's XSPICE d_pwm startup messages (refs/README.md, "tools that scan logs
    for 'error' must allow these two lines"), not a new kind of leniency.
    """
    os.makedirs(outdir, exist_ok=True)
    wrapper = os.path.join(outdir, f"{label}.cir")
    with open(wrapper, "w") as fh:
        fh.write("\n".join(deck_lines) + control)
    try:
        # wrapper is always the absolute deck path (a plain argv/getopt argument to
        # the ngspice binary, not a control-mode token, so it is NOT subject to the
        # lower-casing this docstring warns about); `cwd` only affects how bare,
        # relative tokens INSIDE the .control block (e.g. a checkpoint filename)
        # resolve, via ngspice's own process working directory.
        proc = subprocess.run([ngspice, "-b", wrapper], capture_output=True, text=True,
                               timeout=1800, cwd=cwd)
    except (FileNotFoundError, PermissionError) as e:
        raise CouldNotRun(f"cannot start {ngspice}: {e}") from e
    with open(os.path.join(outdir, f"{label}.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    # "error" catches most failures ngspice reports, but a non-convergent transient
    # (e.g. tolerances tightened past what the deck's smooth tanh-transition models can
    # satisfy, card #452 T2.3 converter investigation) exits 0 and never says "error" —
    # it says "aborted" / "Timestep too small" and silently leaves whatever partial,
    # truncated data it had already written, which downstream would look like a
    # normal-but-wrong result (e.g. va1 extracted as exactly 0, a ZeroDivisionError
    # several calls later in stbcore.middlebrook, not a clear CouldNotRun here) instead
    # of the could-not-run it actually is.
    out = proc.stdout + proc.stderr
    errs = [ln for ln in out.splitlines()
            if re.search(r"\berror\b|\baborted\b|timestep too small|singular matrix", ln, re.I)
            and not any(re.search(p, ln, re.I) for p in allow_patterns)]
    # A nonzero returncode with zero real errs left after allow_patterns filtering
    # is tolerated ONLY when the caller opted in (allow_patterns non-empty) -- e.g.
    # the checkpoint-mode "incomplete or empty netlist" notice above sets ngspice's
    # exit code to 1 even on an otherwise fully successful run. Every existing
    # caller passes no allow_patterns, so a nonzero returncode is exactly as fatal
    # for them as before this parameter existed.
    if errs or (proc.returncode != 0 and not allow_patterns):
        raise CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")


# ---------------------------------------------------------------------------
# mode "loop": Tian two-injection loop gain, transient/DFT version of tools/stb.py
# ---------------------------------------------------------------------------

DEFAULT_COMBINE_THRESHOLD_DB = 25.0


def combine_tian_middlebrook(t, mb, threshold_db=DEFAULT_COMBINE_THRESHOLD_DB):
    """Pick Tian's T below `threshold_db`, Middlebrook's estimate at/above it.

    Tian's T = S/(1-S) (stbcore.tian) amplifies extraction error whenever S is close to
    1, i.e. deep in a high-gain passband (module docstring "Ill-conditioning at high
    loop gain"). Middlebrook's result (stbcore.middlebrook) comes from the same two
    transient runs (v(a) in the v-injection run, i(v.<probe>.vb) in the i-injection run)
    and is not subject to the same S->1 blow-up, so it is the better estimate once |T|
    is large; both are always computed and returned, together with their disagreement
    (dB and degrees, got - ref = middlebrook - tian), so the caller can see how the two
    methods track each other around the threshold, not just which one was picked.

    Returns (combined, source, disagreement_db, disagreement_deg); source is "tian" or
    "middlebrook". threshold_db in [20, 30] is the ill-conditioning region the spec
    calls out (spec.md fra_ref / card #452); 25 dB is the midpoint.
    """
    tdb = stbcore.db(t)
    if tdb < threshold_db:
        combined, source = t, "tian"
    else:
        combined, source = mb, "middlebrook"
    ddb = stbcore.db(mb) - stbcore.db(t)
    dph = math.degrees(cmath.phase(mb)) - math.degrees(cmath.phase(t))
    while dph > 180:
        dph -= 360
    while dph < -180:
        dph += 360
    return combined, source, ddb, dph


def run_loop_point(deck_lines, probe, f, amplitude, settle_time=0.0, settle_periods=20,
                    measure_periods=10, points_per_period=256, outdir=None,
                    ngspice=None, subtract_baseline=True, combine_threshold_db=DEFAULT_COMBINE_THRESHOLD_DB,
                    startup="op", tmax=None, spice_options=""):
    """One frequency point of the Tian two-injection loop gain by transient injection.

    Runs three transients from the same t=0 state (see module docstring): unperturbed
    (both probe sources off), voltage-injection (Vi = SIN, Ii off) and current-injection
    (Ii = SIN, Vi off). Returns a dict with T (Bode return ratio), loopGain (-T), the
    Middlebrook estimate from the same runs, a threshold-combined estimate of the two
    (see combine_tian_middlebrook), the normalised probe readings, and the window
    actually used.

    `subtract_baseline=False` skips the perturbed-minus-unperturbed step; it exists only
    so tests/test_fra.py can show that step is load-bearing (verify.md red-proof).

    `startup` ("op", default, or "uic"): how the three transients get their shared t=0
    state. "op" (unchanged default, what R1/R1b's linear decks need) runs ngspice's DC
    operating point first, as before. "uic" skips it and starts the transient from
    device default (zero) initial conditions instead, exactly like refs/r3_vm_buck.cir
    and refs/r4_pcm_buck.cir's own reference `.tran ... uic` (see their .control blocks)
    — needed on these clocked switching converters because "op" evaluates the PWM ramp's
    PULSE source at its t=0 value (0 V) and finds a real but useless DC equilibrium (the
    comparator latched with pwm=0, vout~0, opout railed at its VCCS's negative limit —
    confirmed by inspecting ngspice's "Initial Transient Solution" printout, card #452
    report), not the periodic switched steady state; `settle_time` (below) then has to
    cover the converter's own natural startup transient into that steady state before
    injection turns on, same as the deck headers' own measured settling times. "uic" is
    still deterministic and identical across all three runs (same zero ICs, same
    circuit), so the "bit-for-bit identical up to TD" property this harness relies on
    (module docstring) still holds.

    `tmax` (card #452 T2.3 converter investigation): the maximum internal transient
    step ngspice's adaptive stepper is allowed to take, independent of `dt` (the
    *output* grid spacing used by `linearize`/the DFT). Historically this harness
    passed `tmax = dt`, which is fine for R1/R1b (dt is already small there) but wrong
    for a switching converter at a *low* excitation frequency: dt = 1/(f*points_per_
    period) grows with 1/f, so at e.g. f=792 Hz, points_per_period=256, dt~=4.9us —
    more than double R3/R4's 2us switching period, i.e. the solver was allowed to step
    clean over switching edges. Default `None` keeps the old `tmax = dt` behaviour
    (unchanged for R1/R1b, which never hit this). Pass an explicit value (e.g. a small
    fraction of the deck's switching period) to bound it independently; see
    run_loop_point's caller in fra_curves.py for the converter default.

    `spice_options` (card #452 T2.3): raw text inserted into the `.control` block right
    after `set numdgt=15`, e.g. `".options reltol=1e-6 trtol=1\\n"`, to test/loosen the
    LTE truncation-error tolerances without changing the harness's own defaults.
    """
    ngspice = ngspice or NGSPICE
    outdir = outdir or os.path.join(DATA_ROOT, "loop_point")
    os.makedirs(outdir, exist_ok=True)
    a, _b = stb.probe_nodes(deck_lines, probe)  # validates the probe exists and is an lgprobe
    p = probe.lower()
    td, dt, t0, twin, tstop = window_times(settle_time, f, settle_periods, measure_periods, points_per_period)
    step_max = tmax if tmax is not None else dt

    base_f = os.path.join(outdir, "base.txt")
    vinj_f = os.path.join(outdir, "vinj.txt")
    iinj_f = os.path.join(outdir, "iinj.txt")
    for fn in (base_f, vinj_f, iinj_f):
        if os.path.exists(fn):
            os.remove(fn)

    if startup not in ("op", "uic"):
        raise CouldNotRun(f"unknown startup {startup!r}, expected 'op' or 'uic'")
    op_line = "op" if startup == "op" else ""
    uic = " uic" if startup == "uic" else ""

    # base_f carries all four quantities (it is subtracted from both perturbed runs);
    # vinj_f/iinj_f only need the three quantities Tian+Middlebrook actually use from
    # that run (v(a) is only meaningful/needed from the v-injection run, i(vb) only from
    # the i-injection run — see stbcore.middlebrook's docstring and stb/lgprobe.lib).
    off = _sin_spec(0.0, f, td)
    on = _sin_spec(amplitude, f, td)
    control = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
{spice_options}
{op_line}
alter @v.{p}.vi[sin] = {off}
alter @i.{p}.ii[sin] = {off}
tran {dt:.12e} {tstop:.12e} 0 {step_max:.12e}{uic}
linearize v({p}.x) i(v.{p}.vi) v({a}) i(v.{p}.vb)
wrdata {base_f} v({p}.x) i(v.{p}.vi) v({a}) i(v.{p}.vb)
alter @v.{p}.vi[sin] = {on}
tran {dt:.12e} {tstop:.12e} 0 {step_max:.12e}{uic}
linearize v({p}.x) i(v.{p}.vi) v({a})
wrdata {vinj_f} v({p}.x) i(v.{p}.vi) v({a})
alter @v.{p}.vi[sin] = {off}
alter @i.{p}.ii[sin] = {on}
tran {dt:.12e} {tstop:.12e} 0 {step_max:.12e}{uic}
linearize v({p}.x) i(v.{p}.vi) i(v.{p}.vb)
wrdata {iinj_f} v({p}.x) i(v.{p}.vi) i(v.{p}.vb)
quit
.endc
.end
"""
    _run_ngspice(deck_lines, control, outdir, "fra_loop", ngspice)

    tb, (vxb, ivb, vab, ibb) = read_tran(base_f, 4)
    tv, (vxv, ivv, vav) = read_tran(vinj_f, 3)
    ti, (vxi, ivi, ibi) = read_tran(iinj_f, 3)
    check_uniform(tb, dt)
    check_uniform(tv, dt)
    check_uniform(ti, dt)

    if subtract_baseline:
        dvx_v = subtract(tv, vxv, tb, vxb, "v-injection v(x)")
        div_v = subtract(tv, ivv, tb, ivb, "v-injection i(vi)")
        dva_v = subtract(tv, vav, tb, vab, "v-injection v(a)")
        dvx_i = subtract(ti, vxi, tb, vxb, "i-injection v(x)")
        div_i = subtract(ti, ivi, tb, ivb, "i-injection i(vi)")
        dib_i = subtract(ti, ibi, tb, ibb, "i-injection i(vb)")
    else:
        dvx_v, div_v, dva_v = vxv, ivv, vav
        dvx_i, div_i, dib_i = vxi, ivi, ibi

    u0 = injected_reference_phasor(f, td)
    norm = amplitude * u0
    v1 = extract_component(tv, dvx_v, f, t0, twin) / norm
    i1 = extract_component(tv, div_v, f, t0, twin) / norm
    va1 = extract_component(tv, dva_v, f, t0, twin) / norm
    v2 = extract_component(ti, dvx_i, f, t0, twin) / norm
    i2 = extract_component(ti, div_i, f, t0, twin) / norm
    ib2 = extract_component(ti, dib_i, f, t0, twin) / norm
    t = stbcore.tian(v1, i1, v2, i2)
    mb = stbcore.middlebrook(v1, va1, i2, ib2)
    combined, source, disagree_db, disagree_deg = combine_tian_middlebrook(t, mb, combine_threshold_db)
    return {
        "f": f, "T": t, "loopGain": -t, "middlebrook": mb,
        "combined": combined, "combined_source": source,
        "disagreement_db": disagree_db, "disagreement_deg": disagree_deg,
        "raw": {"v1": v1, "i1": i1, "v2": v2, "i2": i2, "va1": va1, "ib2": ib2},
        "window": {"td": td, "dt": dt, "t0": t0, "twin": twin, "tstop": tstop},
    }


def _loop_worker(args):
    (deck_lines, probe, f, amplitude, settle_time, settle_periods, measure_periods,
     points_per_period, outdir, ngspice, subtract_baseline, combine_threshold_db, startup,
     tmax, spice_options) = args
    try:
        return ("ok", run_loop_point(deck_lines, probe, f, amplitude, settle_time,
                                      settle_periods, measure_periods, points_per_period,
                                      outdir, ngspice, subtract_baseline, combine_threshold_db,
                                      startup, tmax, spice_options))
    except CouldNotRun as e:
        return ("could-not-run", str(e))


def sweep_loop(deck, probe, freqs, amplitude, settle_time=0.0, settle_periods=20,
               measure_periods=10, points_per_period=256, outdir_root=None,
               workers=None, ngspice=None, subtract_baseline=True,
               combine_threshold_db=DEFAULT_COMBINE_THRESHOLD_DB, startup="op",
               tmax=None, spice_options=""):
    """Run mode "loop" over a frequency list, dispatched in parallel ngspice processes.

    Returns an ordered list, one entry per freqs, each ("ok", result_dict) or
    ("could-not-run", message) — pass/could-not-run kept distinct per point.

    `tmax`/`spice_options`: see run_loop_point (card #452 T2.3); passed through unchanged
    to every frequency point.
    """
    deck_lines = stb.prepare_deck(deck)
    stb.probe_nodes(deck_lines, probe)
    outdir_root = outdir_root or _default_outdir(deck, "loop")
    workers = workers or max(1, (os.cpu_count() or 2) // 2)
    jobs = []
    for f in freqs:
        d = os.path.join(outdir_root, f"f_{f:.6g}")
        jobs.append((deck_lines, probe, f, amplitude, settle_time, settle_periods,
                     measure_periods, points_per_period, d, ngspice, subtract_baseline,
                     combine_threshold_db, startup, tmax, spice_options))
    if workers <= 1 or len(jobs) <= 1:
        return [_loop_worker(j) for j in jobs]
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(_loop_worker, jobs))


def amplitude_halving_check(deck, probe, f, amplitude, **kw):
    """Run mode "loop" at `amplitude` and `amplitude/2`; return (T_a, T_a2, db_diff)."""
    deck_lines = stb.prepare_deck(deck)
    outdir_root = kw.pop("outdir_root", None) or _default_outdir(deck, "halving")
    ra = run_loop_point(deck_lines, probe, f, amplitude,
                         outdir=os.path.join(outdir_root, f"f_{f:.6g}_A"), **kw)
    rb = run_loop_point(deck_lines, probe, f, amplitude / 2,
                         outdir=os.path.join(outdir_root, f"f_{f:.6g}_A2"), **kw)
    return ra["T"], rb["T"], stbcore.db(ra["T"]) - stbcore.db(rb["T"])


# ---------------------------------------------------------------------------
# mode "loop", checkpoint variant (card #472 T10.4): settle the converter ONCE,
# savestate, and run every frequency point (all three of its sub-transients: base,
# v-injection, i-injection) as a `loadstate` continuation from that one checkpoint,
# instead of each of the 3*N sub-transients re-running the whole 0..td settle from
# scratch. Requires an ngspice build with the Enhancement-131 `savestate`/`loadstate`
# commands (ngspice-45.2 does not have them -- module docstring); point NGSPICE at
# the ngspice_augmentation fork build. Opt-in only: run_loop_point/sweep_loop/main's
# existing default (`--checkpoint` unset) behaviour is completely unchanged, and
# still runs correctly (and is still what's tested) on plain ngspice-45.2.
#
# Correctness rests on one fact already true of run_loop_point's design (module
# docstring, "Same starting state"): the deck's Vi/Ii injection sources are exactly
# 0 (a dead short/open) until t=td, for every frequency and for all three sub-runs.
# So the trajectory from t=0 to t=td does not depend on f at all -- settling it once
# and checkpointing at td is not an approximation, it is the same state every one of
# the 3*N sub-transients would otherwise have recomputed from scratch.
# ---------------------------------------------------------------------------

def window_times_checkpoint(td, f, settle_periods, measure_periods, points_per_period):
    """Like window_times, but td is the checkpoint's own (fixed, frequency-
    independent) time rather than being derived from `settle_time` and rounded to
    this frequency's own dt grid. injected_reference_phasor's closed form is exact
    for any td (it does not require td to land on the dt grid), so this is not a
    loss of accuracy versus window_times, just a different (checkpoint-driven)
    source for td. Returns (dt, t0, twin, tstop) -- no td in the tuple, since the
    caller already has it (it is the checkpoint's own fixed value, not per-point).
    """
    dt = dt_for(f, points_per_period)
    period = 1.0 / f
    t0 = td + settle_periods * period
    twin = measure_periods * period
    return dt, t0, twin, t0 + twin


def settle_checkpoint(deck_lines, td, tstep, outdir, ngspice=None, startup="uic", spice_options=""):
    """Run the deck from t=0 to t=td once and savestate to <outdir>/settle.ckpt.

    `startup` ("uic", default, or "op"): see run_loop_point's docstring; "uic" is
    the converter default (R3/R4/R3d/R4d all run their reference `.tran` with uic).
    `tstep` is the internal step passed to `tran` for the settle run only (it does
    not need to match any later frequency point's dt -- see the module note above
    on why loadstate does not require it to). The checkpoint's own filename is
    fixed, bare and lower-case ("settle.ckpt") and passed to `savestate` with
    ngspice's cwd set to `outdir` (via _run_ngspice's `cwd`), NOT as an absolute
    path: ngspice's control-mode tokenizer lower-cases bare command arguments
    (confirmed: a "ckptB.bin" argument was written back as "ckptb.bin"), which
    would silently corrupt any absolute path under this project's ~/Work/...
    convention (upper-case W) into a different, likely-nonexistent path (found
    while building this harness, card #472 T10.4; not otherwise documented in
    the ngspice manual). All other files (the deck wrapper, logs) stay absolute,
    unaffected by cwd.
    """
    ngspice = ngspice or NGSPICE
    os.makedirs(outdir, exist_ok=True)
    ckpt_name = "settle.ckpt"
    ckpt_path = os.path.join(outdir, ckpt_name)
    if os.path.exists(ckpt_path):
        os.remove(ckpt_path)
    if startup not in ("op", "uic"):
        raise CouldNotRun(f"unknown startup {startup!r}, expected 'op' or 'uic'")
    op_line = "op" if startup == "op" else ""
    uic = " uic" if startup == "uic" else ""
    control = f"""
.control
set filetype=ascii
set numdgt=15
{spice_options}
{op_line}
tran {tstep:.12e} {td:.12e} 0 {tstep:.12e}{uic}
savestate {ckpt_name}
.endc
.end
"""
    _run_ngspice(deck_lines, control, outdir, "fra_settle", ngspice, cwd=outdir)
    if not os.path.exists(ckpt_path):
        raise CouldNotRun(f"settle_checkpoint: {ngspice} did not write {ckpt_path}")
    return ckpt_path


def run_loop_point_checkpoint(ckpt_path, deck_lines, probe, f, amplitude, td,
                               settle_periods=20, measure_periods=10, points_per_period=256,
                               outdir=None, ngspice=None, subtract_baseline=True,
                               combine_threshold_db=DEFAULT_COMBINE_THRESHOLD_DB,
                               startup="uic", tmax=None, spice_options=""):
    """run_loop_point's checkpoint variant: the same three sub-transients (base,
    v-injection, i-injection) and the same Tian/Middlebrook extraction, but each
    sub-transient is a fresh `loadstate ckpt_path` continuation from the shared
    settle checkpoint instead of a `tran ... uic` from t=0. Returns the same dict
    shape as run_loop_point (same keys), so callers (write_loop_results etc.) do
    not need to know which variant produced a given point.

    The deck needs a top-level '.tran' card (not a control-mode 'tran' command) so
    loadstate's find_tran_job() has a job to read tstop from (com_checkpoint.c);
    that card is added here, once, ahead of the .control block, with tstop = the
    window end time (t0+twin) shared by all three sub-runs.
    """
    ngspice = ngspice or NGSPICE
    outdir = outdir or os.path.join(DATA_ROOT, "loop_point_ckpt")
    os.makedirs(outdir, exist_ok=True)
    a, _b = stb.probe_nodes(deck_lines, probe)
    p = probe.lower()
    dt, t0, twin, tstop = window_times_checkpoint(td, f, settle_periods, measure_periods, points_per_period)
    step_max = tmax if tmax is not None else dt

    if tstop <= td:
        raise CouldNotRun(f"checkpoint window end {tstop:.6g} s is not after td={td:.6g} s")
    if startup not in ("op", "uic"):
        raise CouldNotRun(f"unknown startup {startup!r}, expected 'op' or 'uic'")
    uic = " uic" if startup == "uic" else ""

    base_f = os.path.join(outdir, "base.txt")
    vinj_f = os.path.join(outdir, "vinj.txt")
    iinj_f = os.path.join(outdir, "iinj.txt")
    for fn in (base_f, vinj_f, iinj_f):
        if os.path.exists(fn):
            os.remove(fn)

    off = _sin_spec(0.0, f, td)
    on = _sin_spec(amplitude, f, td)
    body = "\n".join(deck_lines)
    # loadstate's filename argument must be bare/relative and already lower-case
    # (see settle_checkpoint's docstring) -- ngspice's cwd is set to the
    # checkpoint's own directory below (_run_ngspice's cwd=ckpt_dir); every other
    # file reference here (base_f/vinj_f/iinj_f, the deck wrapper) stays absolute.
    ckpt_dir, ckpt_name = os.path.split(ckpt_path)
    control = f"""{body}
.tran {dt:.12e} {tstop:.12e} 0 {step_max:.12e}{uic}
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
{spice_options}
alter @v.{p}.vi[sin] = {off}
alter @i.{p}.ii[sin] = {off}
loadstate {ckpt_name}
linearize v({p}.x) i(v.{p}.vi) v({a}) i(v.{p}.vb)
wrdata {base_f} v({p}.x) i(v.{p}.vi) v({a}) i(v.{p}.vb)
alter @v.{p}.vi[sin] = {on}
loadstate {ckpt_name}
linearize v({p}.x) i(v.{p}.vi) v({a})
wrdata {vinj_f} v({p}.x) i(v.{p}.vi) v({a})
alter @v.{p}.vi[sin] = {off}
alter @i.{p}.ii[sin] = {on}
loadstate {ckpt_name}
linearize v({p}.x) i(v.{p}.vi) i(v.{p}.vb)
wrdata {iinj_f} v({p}.x) i(v.{p}.vi) i(v.{p}.vb)
.endc
.end
"""
    # loadstate's own continuation (com_checkpoint.c's if_run(ckt, "resume", ...))
    # is what actually runs the transient; _run_ngspice's normal "error"/"aborted"
    # log scan (see its docstring) still applies to whatever ngspice prints.
    _run_ngspice([], control, outdir, "fra_loop_ckpt", ngspice, cwd=ckpt_dir,
                 allow_patterns=(r"incomplete or empty netlist",))

    tb, (vxb, ivb, vab, ibb) = read_tran(base_f, 4)
    tv, (vxv, ivv, vav) = read_tran(vinj_f, 3)
    ti, (vxi, ivi, ibi) = read_tran(iinj_f, 3)
    check_uniform(tb, dt)
    check_uniform(tv, dt)
    check_uniform(ti, dt)

    if subtract_baseline:
        dvx_v = subtract(tv, vxv, tb, vxb, "v-injection v(x)")
        div_v = subtract(tv, ivv, tb, ivb, "v-injection i(vi)")
        dva_v = subtract(tv, vav, tb, vab, "v-injection v(a)")
        dvx_i = subtract(ti, vxi, tb, vxb, "i-injection v(x)")
        div_i = subtract(ti, ivi, tb, ivb, "i-injection i(vi)")
        dib_i = subtract(ti, ibi, tb, ibb, "i-injection i(vb)")
    else:
        dvx_v, div_v, dva_v = vxv, ivv, vav
        dvx_i, div_i, dib_i = vxi, ivi, ibi

    u0 = injected_reference_phasor(f, td)
    norm = amplitude * u0
    v1 = extract_component(tv, dvx_v, f, t0, twin) / norm
    i1 = extract_component(tv, div_v, f, t0, twin) / norm
    va1 = extract_component(tv, dva_v, f, t0, twin) / norm
    v2 = extract_component(ti, dvx_i, f, t0, twin) / norm
    i2 = extract_component(ti, div_i, f, t0, twin) / norm
    ib2 = extract_component(ti, dib_i, f, t0, twin) / norm
    t = stbcore.tian(v1, i1, v2, i2)
    mb = stbcore.middlebrook(v1, va1, i2, ib2)
    combined, source, disagree_db, disagree_deg = combine_tian_middlebrook(t, mb, combine_threshold_db)
    return {
        "f": f, "T": t, "loopGain": -t, "middlebrook": mb,
        "combined": combined, "combined_source": source,
        "disagreement_db": disagree_db, "disagreement_deg": disagree_deg,
        "raw": {"v1": v1, "i1": i1, "v2": v2, "i2": i2, "va1": va1, "ib2": ib2},
        "window": {"td": td, "dt": dt, "t0": t0, "twin": twin, "tstop": tstop},
    }


def _loop_ckpt_worker(args):
    (ckpt_path, deck_lines, probe, f, amplitude, td, settle_periods, measure_periods,
     points_per_period, outdir, ngspice, subtract_baseline, combine_threshold_db,
     startup, tmax, spice_options) = args
    try:
        return ("ok", run_loop_point_checkpoint(ckpt_path, deck_lines, probe, f, amplitude, td,
                                                 settle_periods, measure_periods, points_per_period,
                                                 outdir, ngspice, subtract_baseline, combine_threshold_db,
                                                 startup, tmax, spice_options))
    except CouldNotRun as e:
        return ("could-not-run", str(e))


def sweep_loop_checkpoint(deck, probe, freqs, amplitude, td, settle_tstep, settle_periods=20,
                           measure_periods=10, points_per_period=256, outdir_root=None,
                           workers=None, ngspice=None, subtract_baseline=True,
                           combine_threshold_db=DEFAULT_COMBINE_THRESHOLD_DB, startup="uic",
                           tmax=None, spice_options=""):
    """Settle-once-checkpoint-then-per-point variant of sweep_loop.

    `td`: the fixed checkpoint time (settle time), shared by every frequency point
    -- unlike sweep_loop's `settle_time`, this is not rounded to any per-point dt
    grid (window_times_checkpoint's docstring). `settle_tstep`: the internal `tran`
    step used only for the one settle run (settle_checkpoint's docstring).

    Returns (points, settle_ckpt_path): points is the same ordered
    ("ok"|"could-not-run", ...) list sweep_loop returns; settle_ckpt_path is where
    the one-time checkpoint was written, so a caller (e.g. a test) can report its
    own settle-run wall-clock time separately from the per-point time.
    """
    deck_lines = stb.prepare_deck(deck)
    stb.probe_nodes(deck_lines, probe)
    outdir_root = outdir_root or _default_outdir(deck, "loop_ckpt")
    ckpt_path = settle_checkpoint(deck_lines, td, settle_tstep, outdir_root, ngspice, startup, spice_options)
    workers = workers or max(1, (os.cpu_count() or 2) // 2)
    jobs = []
    for f in freqs:
        d = os.path.join(outdir_root, f"f_{f:.6g}")
        jobs.append((ckpt_path, deck_lines, probe, f, amplitude, td, settle_periods,
                     measure_periods, points_per_period, d, ngspice, subtract_baseline,
                     combine_threshold_db, startup, tmax, spice_options))
    if workers <= 1 or len(jobs) <= 1:
        return [_loop_ckpt_worker(j) for j in jobs], ckpt_path
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(_loop_ckpt_worker, jobs)), ckpt_path


def sweep_loop_adaptive(deck, probe, freqs, amplitude, period,
                         residual_bound=DEFAULT_RESIDUAL_BOUND, cap=DEFAULT_SETTLE_CAP,
                         settle_points_per_period=DEFAULT_SETTLE_POINTS_PER_PERIOD,
                         settle_min_periods=DEFAULT_SETTLE_MIN_PERIODS,
                         settle_n_phases=DEFAULT_SETTLE_N_PHASES, checkpoint=False,
                         outdir_root=None, ngspice=None, startup="uic", tmax=None,
                         spice_options="", **kw):
    """Adaptive-settle variant of sweep_loop/sweep_loop_checkpoint (card #487).

    Finds the settle time with ONE periodicity_scan (tracking v(<probe>.x) and v(a),
    the two node voltages run_loop_point's own base sub-transient extracts from —
    "the measured nodes") instead of taking a caller-guessed fixed `--settle`, then
    reuses that ONE settle time for every frequency point exactly like a fixed
    settle_time would be: via sweep_loop_checkpoint (`checkpoint=True`, pays the
    settle cost once total) or sweep_loop (`checkpoint=False`, every point still
    re-settles from t=0 to that time, same cost model as the old fixed-settle path —
    only `td` is measured rather than guessed).

    `**kw`: forwarded to sweep_loop/sweep_loop_checkpoint (settle_periods,
    measure_periods, points_per_period, workers, subtract_baseline,
    combine_threshold_db).

    Returns (points, settle_info, ckpt_path_or_None). settle_info is
    periodicity_scan's own dict (settle_time, residual, capped, periods_scanned,
    trace_path, node_exprs) — attach it to write_loop_results so every CSV row and
    summary.json record what settle time/residual this sweep actually used.
    """
    deck_lines = stb.prepare_deck(deck)
    a, _b = stb.probe_nodes(deck_lines, probe)
    p = probe.lower()
    outdir_root = outdir_root or _default_outdir(deck, "loop_adaptive")
    settle_info = periodicity_scan(
        deck_lines, [f"v({p}.x)", f"v({a})"], period, residual_bound=residual_bound, cap=cap,
        points_per_period=settle_points_per_period, min_periods=settle_min_periods,
        n_phases=settle_n_phases, outdir=os.path.join(outdir_root, "settle_scan"),
        ngspice=ngspice, startup=startup, tmax=tmax, spice_options=spice_options)
    td = settle_info["settle_time"]
    if checkpoint:
        points, ckpt_path = sweep_loop_checkpoint(
            deck, probe, freqs, amplitude, td=td,
            settle_tstep=dt_for(freqs[len(freqs) // 2], kw.get("points_per_period", 256)),
            outdir_root=outdir_root, ngspice=ngspice, startup=startup, tmax=tmax,
            spice_options=spice_options, **kw)
    else:
        points = sweep_loop(deck, probe, freqs, amplitude, settle_time=td,
                             outdir_root=outdir_root, ngspice=ngspice, startup=startup,
                             tmax=tmax, spice_options=spice_options, **kw)
        ckpt_path = None
    return points, settle_info, ckpt_path


# ---------------------------------------------------------------------------
# mode "transfer": sideband-0 transfer from a named source to a named output
# ---------------------------------------------------------------------------

_TRANS_FUNCS = ("sin", "pulse", "pwl", "exp", "sffm", "am", "trnoise", "trrandom")


def source_dc_value(deck_lines, source):
    """Find `source`'s own DC value in the (unmodified) deck text (card #485).

    Returns the DC value as a string (spice-number syntax preserved as-is, e.g.
    "12", "{vin_v}"), or "0" if the device line gives none (matches SPICE's own
    default). Raises CouldNotRun if `source` is not a V/I device found in
    `deck_lines`, or if it already carries a time-function (SIN/PULSE/PWL/EXP/
    SFFM/AM/TRNOISE/TRRANDOM, or a bracketed `[...]` literal) -- callers must not
    silently clobber an existing time-varying source with their own SIN.

    Only matches a TOP-LEVEL device line (this project's fixture/reference decks
    do not nest the sources this is used on inside subcircuits); a dotted name
    like "v.x1.vi" is matched by its last component ("vi") as a best-effort, not
    verified against the actual subcircuit expansion -- see "what was not done".
    Does not follow `+`-continuation lines (none of this project's decks split a
    source's own line that way; not handled if one does).
    """
    last = source.split(".")[-1].lower()
    for line in deck_lines:
        fields = line.split()
        if not fields:
            continue
        devname = fields[0]
        if devname.lower() not in (source.lower(), last):
            continue
        if devname[0].upper() not in ("V", "I"):
            continue
        if len(fields) < 3:
            continue
        rest = fields[3:]
        joined = " ".join(rest).lower()
        for fn in _TRANS_FUNCS:
            if fn + "(" in joined:
                raise CouldNotRun(
                    f"source_dc_value: {source!r} already has a {fn.upper()}(...) time "
                    f"function in the deck ({line.strip()!r}); refusing to replace it "
                    f"with an injected SIN (card #485)")
        if "[" in joined:
            raise CouldNotRun(
                f"source_dc_value: {source!r} already has a bracketed time-function "
                f"literal in the deck ({line.strip()!r}); refusing to replace it")
        i = 0
        if i < len(rest) and rest[i].upper() == "DC":
            i += 1
            return rest[i] if i < len(rest) else "0"
        if i < len(rest) and rest[i].upper() != "AC":
            return rest[i]
        return "0"
    raise CouldNotRun(f"source_dc_value: {source!r} not found as a top-level V/I device")


def _resolve_dc_numeric(deck_lines, source, vo_token, outdir, ngspice):
    """Resolve a brace-expression DC token (e.g. "{vin_v}") to a plain decimal
    string via a quick ngspice `print @source[dc]` query; a plain numeric token
    (e.g. "12") is returned unchanged.

    Needed because `alter @dev[sin] = [ {vin_v} ... ]` (the bracket-literal form
    `_sin_spec` builds) does NOT evaluate a brace expression the way ordinary
    `.param` substitution elsewhere in a deck does -- confirmed empirically
    (card #485): it fails with "Error: cannot evaluate new parameter value."
    Delegating to ngspice's own `@dev[dc]` query is more robust than re-
    implementing `.param` expression evaluation (nested references, spice-suffix
    numbers like "500k") here.

    Deliberately does NOT run `op` first: `@source[dc]` reads the parsed (pre-
    analysis) parameter value regardless (confirmed: identical output with or
    without a preceding `op`), and some decks' DC operating point is itself
    degenerate/non-convergent by design (e.g. refs/r4_pcm_buck.cir's smooth RS
    latch -- the exact reason run_loop_point/run_transfer_point_conv use
    startup="uic" on these decks in the first place, module docstring) --
    running `op` here hit that same "singular matrix: check node q" trap on R4,
    for a value that does not need a circuit solution to read at all.
    """
    if not (vo_token.startswith("{") and vo_token.endswith("}")):
        return vo_token
    control = f"\n.control\nset numdgt=15\nprint @{source}[dc]\n.endc\n.end\n"
    _run_ngspice(deck_lines, control, outdir, "fra_dcquery", ngspice,
                 allow_patterns=(r"incomplete or empty netlist",))
    log_path = os.path.join(outdir, "fra_dcquery.log")
    with open(log_path) as fh:
        out = fh.read()
    m = re.search(rf"@{re.escape(source)}\[dc\]\s*=\s*([-+0-9.eE]+)", out, re.I)
    if not m:
        raise CouldNotRun(f"_resolve_dc_numeric: could not read @{source}[dc] from {log_path}")
    return m.group(1)


def run_transfer_point(deck_lines, source, output_expr, f, amplitude, settle_time=0.0,
                        settle_periods=20, measure_periods=10, points_per_period=256,
                        outdir=None, ngspice=None):
    """Sideband-0 transfer H(f) = Y(f)/U(f) from an injected source to `output_expr`.

    `source` is the ngspice device name to alter, exactly as `listing expand` would show
    it (e.g. "vin" for a top-level source, "v.x1.vi" for one inside a subcircuit); it
    must already carry a DC value in the deck (its SIN function is added/replaced here).
    `output_expr` is any ngspice vector expression ngspice's own `print`/`wrdata` would
    accept, e.g. "v(out)" or "i(rsense)". For control-to-output, line-to-output and
    output-impedance measurements (`pac_sb0` (c) in spec.md); not yet exercised by
    tests/test_fra.py, which only validates mode "loop" on the R1 linear decks.

    card #485: `source`'s own DC value (read from the deck via source_dc_value, which
    raises CouldNotRun if `source` is not found or already carries a time function) is
    now preserved as the injected SIN's VO -- previously hard-coded to 0, which silently
    removed the source's DC bias for both the "base" and "pert" runs (wrong for any
    source with a nonzero DC, e.g. Vin on a line-to-output measurement, by however much
    the circuit's response actually depends on that bias -- 30-60 dB, measured on R3/R4).
    """
    ngspice = ngspice or NGSPICE
    outdir = outdir or os.path.join(DATA_ROOT, "transfer_point")
    os.makedirs(outdir, exist_ok=True)
    td, dt, t0, twin, tstop = window_times(settle_time, f, settle_periods, measure_periods, points_per_period)
    vo = _resolve_dc_numeric(deck_lines, source, source_dc_value(deck_lines, source), outdir, ngspice)

    base_f = os.path.join(outdir, "base.txt")
    pert_f = os.path.join(outdir, "pert.txt")
    for fn in (base_f, pert_f):
        if os.path.exists(fn):
            os.remove(fn)

    off = _sin_spec(0.0, f, td, vo=vo)
    on = _sin_spec(amplitude, f, td, vo=vo)
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
    _run_ngspice(deck_lines, control, outdir, "fra_transfer", ngspice)
    tb, (yb,) = read_tran(base_f, 1)
    tp, (yp,) = read_tran(pert_f, 1)
    check_uniform(tb, dt)
    check_uniform(tp, dt)
    dy = subtract(tp, yp, tb, yb, output_expr)
    u0 = injected_reference_phasor(f, td)
    h = extract_component(tp, dy, f, t0, twin) / (amplitude * u0)
    return {"f": f, "H": h, "window": {"td": td, "dt": dt, "t0": t0, "twin": twin, "tstop": tstop}}


def _transfer_worker(args):
    (deck_lines, source, output_expr, f, amplitude, settle_time, settle_periods,
     measure_periods, points_per_period, outdir, ngspice) = args
    try:
        return ("ok", run_transfer_point(deck_lines, source, output_expr, f, amplitude,
                                          settle_time, settle_periods, measure_periods,
                                          points_per_period, outdir, ngspice))
    except CouldNotRun as e:
        return ("could-not-run", str(e))


def sweep_transfer(deck, source, output_expr, freqs, amplitude, settle_time=0.0,
                    settle_periods=20, measure_periods=10, points_per_period=256,
                    outdir_root=None, workers=None, ngspice=None):
    deck_lines = stb.prepare_deck(deck)
    outdir_root = outdir_root or _default_outdir(deck, "transfer")
    workers = workers or max(1, (os.cpu_count() or 2) // 2)
    jobs = []
    for f in freqs:
        d = os.path.join(outdir_root, f"f_{f:.6g}")
        jobs.append((deck_lines, source, output_expr, f, amplitude, settle_time,
                     settle_periods, measure_periods, points_per_period, d, ngspice))
    if workers <= 1 or len(jobs) <= 1:
        return [_transfer_worker(j) for j in jobs]
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(_transfer_worker, jobs))


def sweep_transfer_adaptive(deck, source, output_expr, freqs, amplitude, period,
                             residual_bound=DEFAULT_RESIDUAL_BOUND, cap=DEFAULT_SETTLE_CAP,
                             settle_points_per_period=DEFAULT_SETTLE_POINTS_PER_PERIOD,
                             settle_min_periods=DEFAULT_SETTLE_MIN_PERIODS,
                             settle_n_phases=DEFAULT_SETTLE_N_PHASES, outdir_root=None,
                             workers=None, ngspice=None, **kw):
    """Adaptive-settle variant of sweep_transfer (card #487): one periodicity_scan on
    `output_expr` (the node this sweep will itself measure) picks the settle time,
    reused as a fixed settle_time for every frequency point (run_transfer_point always
    starts from `op`, never `uic` — see its own docstring — so the settle scan here does
    too, to match exactly what the sweep itself will run).

    Returns (points, settle_info) — see sweep_loop_adaptive's docstring for the shape
    of settle_info; there is no checkpoint variant for mode transfer (run_transfer_point
    has no loadstate/savestate support, unlike run_loop_point/run_loop_point_checkpoint).
    """
    deck_lines = stb.prepare_deck(deck)
    outdir_root = outdir_root or _default_outdir(deck, "transfer_adaptive")
    settle_info = periodicity_scan(
        deck_lines, [output_expr], period, residual_bound=residual_bound, cap=cap,
        points_per_period=settle_points_per_period, min_periods=settle_min_periods,
        n_phases=settle_n_phases, outdir=os.path.join(outdir_root, "settle_scan"),
        ngspice=ngspice, startup="op")
    td = settle_info["settle_time"]
    points = sweep_transfer(deck, source, output_expr, freqs, amplitude, settle_time=td,
                             outdir_root=outdir_root, workers=workers, ngspice=ngspice, **kw)
    return points, settle_info


# ---------------------------------------------------------------------------
# frequency list parsing, output, CLI
# ---------------------------------------------------------------------------

def _spice_num(s):
    m = re.match(r"^([-+0-9.eE]+)(meg|[tgkmunpf])?$", s.strip().lower())
    if not m:
        raise ValueError(f"not a number: {s!r}")
    mult = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "m": 1e-3,
             "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}
    return float(m.group(1)) * (mult[m.group(2)] if m.group(2) else 1.0)


def parse_freqs(spec):
    """"dec N f0 f1" | "oct N f0 f1" | "lin N f0 f1" | a comma list of frequencies."""
    parts = spec.split()
    if len(parts) == 4 and parts[0].lower() in ("dec", "oct", "lin"):
        kind, n, f0, f1 = parts[0].lower(), int(parts[1]), _spice_num(parts[2]), _spice_num(parts[3])
        if kind == "lin":
            return [f0] if n <= 1 else [f0 + (f1 - f0) * k / (n - 1) for k in range(n)]
        base = 10.0 if kind == "dec" else 2.0
        npts = max(1, round(math.log(f1 / f0, base) * n))
        return [f0 * base ** (k / n) for k in range(npts + 1)]
    return [_spice_num(x) for x in spec.split(",") if x.strip()]


def _default_outdir(deck, mode):
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    return os.path.join(DATA_ROOT, mode, f"{os.path.splitext(os.path.basename(deck))[0]}-{stamp}")


def write_loop_results(points, freqs, outdir, settle_info=None):
    """settle_info (default None, unchanged behaviour): periodicity_scan's own dict
    (card #487), attached when the sweep used --adaptive-settle. When given, every CSV
    row (ok or could-not-run alike — "reported per point", card #487's done-when) also
    carries the ONE settle_time/residual/capped this whole sweep used (fixed once,
    before any frequency point ran — see sweep_loop_adaptive), and summary.json gets a
    top-level "settle" key with the full dict. Omitting settle_info (the default)
    reproduces the exact pre-#487 CSV/summary.json shape byte for byte."""
    os.makedirs(outdir, exist_ok=True)
    freqs_ok, t_ok, mb_ok, comb_ok = [], [], [], []
    settle_cols = ",settle_time_s,settle_residual,settle_capped" if settle_info is not None else ""
    settle_vals = (f",{settle_info['settle_time']:.9e},{settle_info['residual']:.6e},"
                    f"{settle_info['capped']}") if settle_info is not None else ""
    with open(os.path.join(outdir, "fra.csv"), "w") as fh:
        fh.write("freq,status,T_re,T_im,T_db,T_phase_deg,loopGain_db,loopGain_phase_deg,"
                  "middlebrook_db,middlebrook_phase_deg,combined_db,combined_phase_deg,"
                  f"combined_source,disagreement_db,disagreement_deg{settle_cols},detail\n")
        for f, (status, val) in zip(freqs, points):
            if status == "ok":
                t, mb, comb = val["T"], val["middlebrook"], val["combined"]
                fh.write(f"{f:.9e},ok,{t.real:.9e},{t.imag:.9e},{stbcore.db(t):.6f},"
                         f"{math.degrees(cmath.phase(t)):.6f},{stbcore.db(-t):.6f},"
                         f"{math.degrees(cmath.phase(-t)):.6f},"
                         f"{stbcore.db(mb):.6f},{math.degrees(cmath.phase(mb)):.6f},"
                         f"{stbcore.db(comb):.6f},{math.degrees(cmath.phase(comb)):.6f},"
                         f"{val['combined_source']},{val['disagreement_db']:.6f},"
                         f"{val['disagreement_deg']:.6f}{settle_vals},\n")
                freqs_ok.append(f)
                t_ok.append(t)
                mb_ok.append(mb)
                comb_ok.append(comb)
            else:
                fh.write(f"{f:.9e},could-not-run,,,,,,,,,,,,,{settle_vals},{val}\n")
    summary = {
        "margins": stbcore.margins(freqs_ok, t_ok) if len(freqs_ok) >= 2 else None,
        "middlebrook_margins": stbcore.margins(freqs_ok, mb_ok) if len(freqs_ok) >= 2 else None,
        "combined_margins": stbcore.margins(freqs_ok, comb_ok) if len(freqs_ok) >= 2 else None,
        "n_ok": len(freqs_ok), "n_total": len(freqs),
    }
    if settle_info is not None:
        summary["settle"] = settle_info
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=lambda z: [z.real, z.imag] if isinstance(z, complex) else str(z))
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("deck")
    ap.add_argument("--mode", choices=["loop", "transfer"], default="loop")
    ap.add_argument("--probe", help="lgprobe instance name (mode loop), e.g. XP")
    ap.add_argument("--source", help="device to inject into (mode transfer), e.g. vin or v.x1.vi")
    ap.add_argument("--output", help="ngspice vector expression to measure (mode transfer), e.g. v(out)")
    ap.add_argument("--freqs", required=True, help='"dec N f0 f1" | "oct N f0 f1" | "lin N f0 f1" | "f1,f2,..."')
    ap.add_argument("--amplitude", type=float, default=1e-3)
    ap.add_argument("--settle", type=float, default=0.0, help="seconds before the injection turns on")
    ap.add_argument("--settle-periods", type=float, default=20)
    ap.add_argument("--measure-periods", type=float, default=10)
    ap.add_argument("--points-per-period", type=int, default=256)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--halving", action="store_true", help="also run at amplitude/2 and report the dB difference")
    ap.add_argument("--out", help="output directory (default: under ~/Work/data/ngspice_augmentation/fra)")
    ap.add_argument("--startup", choices=["op", "uic"], default="op",
                     help="mode loop only: 'op' (default, linear circuits) or 'uic' for clocked switching "
                          "converters whose DC operating point is a degenerate non-switching equilibrium "
                          "(see run_loop_point's docstring); combine with --settle for the converter's own "
                          "startup time")
    ap.add_argument("--checkpoint", action="store_true",
                     help="mode loop only, opt-in (card #472 T10.4): settle once to --settle (used here as "
                          "the fixed checkpoint time td, not per-point-rounded) and savestate, then run every "
                          "frequency point's three sub-transients as a loadstate continuation from that one "
                          "checkpoint instead of each re-running the 0..td settle from scratch. Needs an "
                          "ngspice build with savestate/loadstate (Enhancement-131; ngspice-45.2 does not have "
                          "them -- point NGSPICE at the ngspice_augmentation fork build). --settle must be > 0 "
                          "(or, combined with --adaptive-settle, the settle time it measures).")
    ap.add_argument("--tmax", type=float, default=None,
                     help="max internal transient step (independent of the output dt grid); default None "
                          "keeps tmax=dt (see run_loop_point's docstring, 'tmax')")
    ap.add_argument("--spice-options", default="",
                     help='raw .control-block text inserted after "set numdgt=15", e.g. ".options trtol=1"')
    ap.add_argument("--adaptive-settle", action="store_true",
                     help="opt-in (card #487): ignore --settle and instead find it by settling the "
                          "UNPERTURBED circuit once until its periodicity residual (relative to that "
                          "window's own ripple, over the nodes this run itself measures) is <= "
                          "--residual-bound, capped at --settle-cap; needs --settle-period (the circuit's "
                          "OWN natural period, e.g. 1/fs for a switching converter -- NOT an injection "
                          "period). The old fixed --settle path (this flag unset) is completely unchanged.")
    ap.add_argument("--settle-period", type=float, default=None,
                     help="--adaptive-settle only: the circuit's own natural period")
    ap.add_argument("--residual-bound", type=float, default=DEFAULT_RESIDUAL_BOUND,
                     help=f"--adaptive-settle only: max|v(t)-v(t-period)|/ripple bound "
                          f"(default {DEFAULT_RESIDUAL_BOUND:g} = {DEFAULT_RESIDUAL_BOUND*100:g}%%)")
    ap.add_argument("--settle-cap", type=float, default=DEFAULT_SETTLE_CAP,
                     help=f"--adaptive-settle only: give up (report capped) after this much settle time "
                          f"(default {DEFAULT_SETTLE_CAP:g}s)")
    args = ap.parse_args()

    freqs = parse_freqs(args.freqs)
    kw = dict(settle_time=args.settle, settle_periods=args.settle_periods,
              measure_periods=args.measure_periods, points_per_period=args.points_per_period,
              workers=args.workers)

    if args.mode == "loop":
        if not args.probe:
            ap.error("--probe is required for mode loop")
        if args.adaptive_settle and args.settle_period is None:
            ap.error("--adaptive-settle needs --settle-period")
        if args.checkpoint and not args.adaptive_settle and args.settle <= 0:
            ap.error("--checkpoint needs --settle > 0 (the checkpoint time td), or --adaptive-settle")
        outdir = args.out or _default_outdir(
            args.deck, "loop" + ("_adaptive" if args.adaptive_settle else "") + ("_ckpt" if args.checkpoint else ""))
        settle_info = None
        try:
            if args.adaptive_settle:
                adapt_kw = dict(settle_periods=args.settle_periods, measure_periods=args.measure_periods,
                                 points_per_period=args.points_per_period, workers=args.workers)
                points, settle_info, ckpt_path = sweep_loop_adaptive(
                    args.deck, args.probe, freqs, args.amplitude, args.settle_period,
                    residual_bound=args.residual_bound, cap=args.settle_cap, checkpoint=args.checkpoint,
                    outdir_root=outdir, startup=args.startup, tmax=args.tmax,
                    spice_options=args.spice_options, **adapt_kw)
                print(f"adaptive settle: td={settle_info['settle_time']:.6g}s "
                      f"residual={settle_info['residual']:.4%} "
                      f"({'CAPPED' if settle_info['capped'] else 'bound met'}, "
                      f"{settle_info['periods_scanned']} periods scanned)")
                if ckpt_path:
                    print(f"settle checkpoint: {ckpt_path}")
            elif args.checkpoint:
                ckpt_kw = dict(settle_periods=args.settle_periods, measure_periods=args.measure_periods,
                                points_per_period=args.points_per_period, workers=args.workers)
                points, ckpt_path = sweep_loop_checkpoint(
                    args.deck, args.probe, freqs, args.amplitude, td=args.settle,
                    settle_tstep=dt_for(freqs[len(freqs) // 2], args.points_per_period),
                    outdir_root=outdir, startup=args.startup, tmax=args.tmax,
                    spice_options=args.spice_options, **ckpt_kw)
                print(f"settle checkpoint: {ckpt_path}")
            else:
                points = sweep_loop(args.deck, args.probe, freqs, args.amplitude, outdir_root=outdir,
                                     startup=args.startup, tmax=args.tmax,
                                     spice_options=args.spice_options, **kw)
        except CouldNotRun as e:
            print(f"COULD NOT RUN: {e}", file=sys.stderr)
            return 2
        summary = write_loop_results(points, freqs, outdir, settle_info=settle_info)
        n_bad = sum(1 for s, _ in points if s != "ok")
        if summary["combined_margins"]:
            pm, gm = summary["combined_margins"]["pm"], summary["combined_margins"]["gm"]
            print(f"phase margin (combined Tian/Middlebrook): {pm['pm_deg']:.3f} deg at {pm['freq']:.6g} Hz"
                  if pm else "phase margin: no 0 dB crossing")
            print(f"gain margin (combined Tian/Middlebrook): {gm['gm_db']:.3f} dB at {gm['freq']:.6g} Hz"
                  if gm else "gain margin: no -180 deg crossing")
        print(f"{summary['n_ok']}/{summary['n_total']} frequency points ok"
              + (f", {n_bad} could-not-run" if n_bad else ""))
        if args.halving:
            f0 = freqs[len(freqs) // 2]
            settle_kw = dict(kw)
            if args.adaptive_settle:
                settle_kw["settle_time"] = settle_info["settle_time"]
            _, _, dbdiff = amplitude_halving_check(args.deck, args.probe, f0, args.amplitude,
                                                     outdir_root=os.path.join(outdir, "halving"), **settle_kw)
            print(f"amplitude-halving check at {f0:.6g} Hz: {dbdiff:+.4f} dB")
        print(f"results: {outdir}")
        return 0 if n_bad == 0 else 2

    else:
        if not args.source or not args.output:
            ap.error("--source and --output are required for mode transfer")
        if args.adaptive_settle and args.settle_period is None:
            ap.error("--adaptive-settle needs --settle-period")
        outdir = args.out or _default_outdir(args.deck, "transfer" + ("_adaptive" if args.adaptive_settle else ""))
        settle_info = None
        try:
            if args.adaptive_settle:
                points, settle_info = sweep_transfer_adaptive(
                    args.deck, args.source, args.output, freqs, args.amplitude, args.settle_period,
                    residual_bound=args.residual_bound, cap=args.settle_cap, outdir_root=outdir,
                    workers=args.workers, settle_periods=args.settle_periods,
                    measure_periods=args.measure_periods, points_per_period=args.points_per_period)
                print(f"adaptive settle: td={settle_info['settle_time']:.6g}s "
                      f"residual={settle_info['residual']:.4%} "
                      f"({'CAPPED' if settle_info['capped'] else 'bound met'}, "
                      f"{settle_info['periods_scanned']} periods scanned)")
            else:
                points = sweep_transfer(args.deck, args.source, args.output, freqs, args.amplitude,
                                         outdir_root=outdir, **kw)
        except CouldNotRun as e:
            print(f"COULD NOT RUN: {e}", file=sys.stderr)
            return 2
        os.makedirs(outdir, exist_ok=True)
        n_bad = 0
        settle_cols = ",settle_time_s,settle_residual,settle_capped" if settle_info is not None else ""
        settle_vals = (f",{settle_info['settle_time']:.9e},{settle_info['residual']:.6e},"
                        f"{settle_info['capped']}") if settle_info is not None else ""
        with open(os.path.join(outdir, "fra.csv"), "w") as fh:
            fh.write(f"freq,status,H_re,H_im,H_db,H_phase_deg{settle_cols},detail\n")
            for f, (status, val) in zip(freqs, points):
                if status == "ok":
                    h = val["H"]
                    fh.write(f"{f:.9e},ok,{h.real:.9e},{h.imag:.9e},{stbcore.db(h):.6f},"
                             f"{math.degrees(cmath.phase(h)):.6f}{settle_vals},\n")
                else:
                    n_bad += 1
                    fh.write(f"{f:.9e},could-not-run,,,,{settle_vals},{val}\n")
        if settle_info is not None:
            with open(os.path.join(outdir, "summary.json"), "w") as fh:
                json.dump({"settle": settle_info, "n_ok": len(freqs) - n_bad, "n_total": len(freqs)}, fh, indent=2)
        print(f"{len(freqs) - n_bad}/{len(freqs)} frequency points ok"
              + (f", {n_bad} could-not-run" if n_bad else ""))
        print(f"results: {outdir}")
        return 0 if n_bad == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
