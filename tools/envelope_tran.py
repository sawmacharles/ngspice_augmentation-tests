#!/usr/bin/env python3
"""envelope_tran.py -- envelope-following acceleration for long switching-
converter transients (soft-start, slow load/line ramps) that span thousands
of switching periods (card #512, transient acceleration round 10).

See docs/envelope.md for the full design writeup (literature this is based
on, why mechanism (a) -- a Python driver over savestate/loadstate -- was
chosen over a new C `.tran ... envelope=<M>` engine option, and the
accuracy/speed/refusal results measured on tests/fixtures/env_r3d_softstart*
.cir). Short version of the algorithm (Kundert/White/Sangiovanni-Vincentelli
1988's envelope-following, Petzold 1981's quasi-envelope, and Gear/Kevrekidis
projective integration all describe the same shape: sample the circuit's
FULL state vector once per switching period -- the "stroboscopic map" this
project's own `.pss method=shoot` (pssshoot.c) already uses for a DRIVEN
periodic orbit -- and treat the per-period difference as a slowly varying
"envelope" that can be extrapolated (projected) several periods forward
without re-simulating every one of them):

  1. Simulate 2 real periods from the current checkpoint (or from t=0, uic,
     for the very first window) -- x0 at cycle n, x1 at cycle n+1.
  2. SMOOTHNESS CHECK: compare ||x1-x0|| against the MEDIAN of the last 8
     ACCEPTED per-period state deltas (not a fixed fraction of ||x1||: this
     circuit's Newton-unknown vector always includes the constant 12V input
     rail, so ||x1|| never shrinks even when vout is still near 0V right at
     the start of soft-start, which would make a fixed relative-to-||x1||
     threshold far too loose exactly when it matters most). If this
     period's delta exceeds `--growth-tol` (default 5x) that recent median,
     REFUSE to jump: fall through to plain, period-by-period simulation for
     this window (M=1) and re-sample. tests/fixtures/
     env_r3d_softstart_redproof.cir's 0->4A load step produces a
     single-period delta tens of times its neighbours' -- see docs/
     envelope.md's own measured margin.
  3. Otherwise, take an explicit projective (forward-Euler-in-cycle-index)
     step: x_pred = x1 + M*(x1-x0), landing at cycle n+1+M. `ckpt_io.py`
     patches ONLY `CKTtime` (-> the landing time, so slow-time sources like
     the soft-start ramp itself evaluate correctly there) and `CKTrhsOld`
     (-> x_pred, the full Newton-unknown vector -- see docs/envelope.md,
     "What gets extrapolated, and what does not" for why the per-device
     `CKTstates[]` integration history is deliberately left UNCHANGED
     rather than guessed at). `CKTorder` is forced to 1 and a single
     breakpoint is inserted at the landing time, so the engine re-derives a
     consistent history the same way it already does at any other state
     discontinuity (docs/checkpoint.md's patches 0200/0201 case study).
  4. Simulate `--n-check` (default 4) more real periods from the landing
     point. ERROR CONTROL: compare the LAST TWO check periods' observed
     slope against the pre-jump slope (x1-x0) used to extrapolate -- not
     the first two: measured directly (docs/envelope.md, "the landing
     kick") that the 1-2 periods right after any jump carry their own real
     settling transient (CKTstates[] is deliberately left unmodified by
     the patch, so it is briefly inconsistent with the jumped-forward
     CKTrhsOld/CKTtime -- the compensator's own slow integrator pole, not
     noise, dominates how long this takes to relax). If the last two
     checked periods' slope still disagrees with the pre-jump slope by
     more than `--slope-tol` (default 3.0 = 300% relative -- loose enough
     to accept the landing kick's own characterised, bounded excess,
     tight enough to refuse tests/fixtures/env_r3d_softstart_redproof.cir's
     load step, which is an order of magnitude larger again), the jump
     overshot a smoothly-continuable regime -- quarter M and retry from
     the same x1; if every retry still disagrees, give up on jumping this
     window (M=1 passthrough) rather than accept an unvalidated guess.
  5. Continue from the last check period's own checkpoint; repeat until
     reaching --tstop.

Every window that is actually integrated (the 2-period samples, the
2-period checks, and any M=1 passthrough stretch) is real ngspice output,
concatenated into one combined waveform file with GAPS where periods were
skipped -- the same shape Cadence's own SpectreRF envelope-following user
guide describes for its own `td` (time-domain) output file ("big gaps
between integrated clock cycles").
"""
import argparse
import bisect
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import ckpt_io as ck   # noqa: E402
import stb             # noqa: E402  (prepare_deck)

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation")
HDR = "set filetype=ascii\nset wr_singlescale\nset wr_vecnames\n"


class CouldNotRun(RuntimeError):
    pass


def _run(ngspice, deck_text, outdir, label, timeout=600):
    # `-b <relative name>` with cwd=outdir, NOT `-b <outdir-joined path>`:
    # found the hard way (an interleaved A/B timing run that, for the first
    # time, passed a RELATIVE `outdir`) that joining outdir onto the deck
    # path AND setting cwd=outdir double-nests it (ngspice then looks for
    # "<outdir>/<outdir>/<label>.cir" relative to its own cwd and finds
    # nothing) whenever outdir is itself relative -- absolute outdirs never
    # showed this because the join is idempotent-looking (os.path.join with
    # an absolute second component still works) only by accident of how
    # every prior caller in this project happened to pass paths.
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{label}.cir")
    with open(path, "w") as fh:
        fh.write(deck_text)
    t0 = time.time()
    proc = subprocess.run([ngspice, "-b", f"{label}.cir"], capture_output=True, text=True,
                           timeout=timeout, cwd=outdir)
    wall = time.time() - t0
    with open(os.path.join(outdir, f"{label}.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    return proc, wall


def _danger(proc):
    return re.search(r"\baborted\b|timestep too small|singular matrix|fatal error",
                      proc.stdout + proc.stderr, re.I)


def _load_rows(path, ncols):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.split()
            if len(parts) == 1 + ncols:
                try:
                    rows.append(tuple(float(x) for x in parts))
                except ValueError:
                    continue
    return rows


# ---------------------------------------------------------------------------
# one simulated window: either a fresh uic start, or a loadstate continuation
# (optionally from a PATCHED checkpoint -- the envelope jump itself)
# ---------------------------------------------------------------------------

def run_window(ngspice, body, outdir, label, t_target, tstep, vectors,
               ckpt_in=None, patch=None, uic=False, timeout=600):
    """Run ngspice from t=0 (uic=True, ckpt_in=None) or from `ckpt_in`
    (optionally patched -- envelope jump) to t_target. Always savestates at
    t_target as f"{label}.ckpt" and wrdata's `vectors` as f"{label}.txt".
    Returns a result dict; raises CouldNotRun on any danger string.

    NOTE on checkpoint filenames in the control block: ngspice's
    control-block tokenizer LOWERCASES unquoted (and, found the hard way,
    also quoted -- quoting just leaves literal quote characters in the
    resulting filename instead) bareword arguments, including file paths --
    confirmed directly (`savestate MixedCase.bin` writes `mixedcase.bin`).
    This project's own data root is `~/Work/data/...` (capital W) and the
    scratch/data directories this tool runs under can have mixed-case
    components too, so embedding an ABSOLUTE checkpoint path in the control
    text silently writes/reads the WRONG (lowercased, nonexistent) file.
    Fixed the same way `tran_from_pss.py`/`test_checkpoint.py` already
    avoid it: every checkpoint filename in the control text is a bare
    RELATIVE name (never mixed-case-sensitive on this project's own
    lowercase `h%04d_...` label scheme) with the ngspice subprocess's `cwd`
    set to `outdir` (`_run()`, above) -- the absolute, possibly-mixed-case
    path is used only on the Python side (file existence checks,
    `ckpt_io` reads/writes), never passed to ngspice as text."""
    vecs = " ".join(vectors)
    ckpt_out = os.path.join(outdir, f"{label}.ckpt")
    ckpt_out_rel = os.path.basename(ckpt_out)

    if ckpt_in is None:
        assert uic
        ctl = (f"\n.tran {tstep:g} {t_target:g} 0 uic\n"
               f".control\n{HDR}run\nwrdata {label}.txt {vecs}\n"
               f"savestate {ckpt_out_rel}\n.endc\n.end\n")
    else:
        if patch is not None:
            st = ck.read_checkpoint(ckpt_in)
            st.update(patch)
            ckpt_use = os.path.join(outdir, f"{label}.in.ckpt")
            ck.write_checkpoint(ckpt_use, st)
        else:
            ckpt_use = ckpt_in
        ckpt_use_rel = os.path.basename(ckpt_use)
        ctl = (f"\n.tran {tstep:g} {t_target:g} 0\n"
               f".control\n{HDR}loadstate {ckpt_use_rel}\nwrdata {label}.txt {vecs}\n"
               f"savestate {ckpt_out_rel}\n.endc\n.end\n")

    proc, wall = _run(ngspice, body + ctl, outdir, label, timeout=timeout)
    if _danger(proc):
        raise CouldNotRun(f"{label}: danger string in ngspice output (see {outdir}/{label}.log)")
    if not os.path.exists(ckpt_out):
        raise CouldNotRun(f"{label}: no checkpoint written (see {outdir}/{label}.log)")
    ncols = len(vectors)
    rows = _load_rows(os.path.join(outdir, f"{label}.txt"), ncols)
    return dict(ckpt=ckpt_out, data=os.path.join(outdir, f"{label}.txt"), rows=rows,
                wall=wall, t_target=t_target, label=label)


def _vec_at(rows, t, col):
    """Interpolated value of column `col` at time t (0-based, after time col)."""
    times = [r[0] for r in rows]
    i = bisect.bisect_left(times, t)
    if not rows:
        return None
    if i <= 0:
        return rows[0][1 + col]
    if i >= len(rows):
        return rows[-1][1 + col]
    t0, t1 = times[i - 1], times[i]
    if t1 == t0:
        return rows[i][1 + col]
    frac = (t - t0) / (t1 - t0)
    return rows[i - 1][1 + col] + (rows[i][1 + col] - rows[i - 1][1 + col]) * frac


def _period_stats(rows, t0, t1, col):
    """(avg, max, min) of column `col` over rows with t0 <= t <= t1
    (trapezoidal average)."""
    pts = [(r[0], r[1 + col]) for r in rows if t0 - 1e-15 <= r[0] <= t1 + 1e-15]
    if len(pts) < 2:
        v = _vec_at(rows, (t0 + t1) / 2, col)
        return v, v, v
    area = 0.0
    for (ta, va), (tb, vb) in zip(pts, pts[1:]):
        area += (va + vb) * 0.5 * (tb - ta)
    avg = area / (pts[-1][0] - pts[0][0])
    vals = [v for _, v in pts]
    return avg, max(vals), min(vals)


# ---------------------------------------------------------------------------
# the envelope-following driver
# ---------------------------------------------------------------------------

def envelope_follow(deck_path, outdir, tclk, tstop, tstep, vectors,
                     m_periods=10, m_min=1, n_check=4,
                     growth_tol=5.0, slope_tol=3.0, abs_floor=1e-3,
                     ngspice=None, max_hops=100000):
    """Run the full envelope-following driver from t=0 to tstop on
    `deck_path`. Returns a report dict: `windows` (every real-simulated
    segment, in order, each a dict from run_window() plus 'kind' in
    {'sample','check','passthrough'} and 'jumped_from'/'jumped_to'),
    `jumps` (accepted jump records), `refusals` (windows where a jump was
    attempted and refused, with why), wall totals, and `t_final`."""
    ngspice = ngspice or NGSPICE
    os.makedirs(outdir, exist_ok=True)
    body = "\n".join(stb.prepare_deck(deck_path))

    windows = []   # every real-simulated period, INCLUDING abandoned/refused
                    # jump-check attempts -- for audit/timing only
    path = []      # only the periods on the FINAL accepted trajectory, in
                    # chronological order, each with .t0/.t1 set -- what
                    # validate_against_reference() below compares
    jumps = []
    refusals = []
    hop = [0]  # boxed so the nested helper can increment it

    def norm_diff(a, b):
        return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))

    def norm(a):
        return math.sqrt(sum(p * p for p in a))

    def sample_one_period(ckpt_in, t_from, kind="sample"):
        """Simulate exactly one more period from ckpt_in (t_from -> t_from+tclk).
        Returns (window_dict, full_state_vector, vout_avg_over_that_period).
        Appends to `windows` (audit) but NOT to `path` -- callers add to
        `path` explicitly once they know the period is actually part of the
        final trajectory (a jump-check attempt that gets refused must NOT
        be validated against the reference: it is a discarded branch)."""
        land = t_from + tclk
        w = run_window(ngspice, body, outdir, f"h{hop[0]:04d}_{kind}", land, tstep,
                        vectors, ckpt_in=ckpt_in)
        hop[0] += 1
        w["kind"] = kind
        w["t0"], w["t1"] = t_from, land
        windows.append(w)
        state = ck.read_checkpoint(w["ckpt"])["CKTrhsOld"]
        avg, _, _ = _period_stats(w["rows"], t_from, land, 0)
        return w, state, avg

    # --- bootstrap: 2 real sample periods from t=0, uic ---
    w0 = run_window(ngspice, body, outdir, f"h{hop[0]:04d}_boot", tclk, tstep, vectors,
                     ckpt_in=None, uic=True)
    hop[0] += 1
    w0["kind"] = "sample"
    w0["t0"], w0["t1"] = 0.0, tclk
    windows.append(w0)
    path.append(w0)
    x_prev = ck.read_checkpoint(w0["ckpt"])["CKTrhsOld"]
    avg_prev, _, _ = _period_stats(w0["rows"], 0.0, tclk, 0)
    t = tclk
    ckpt = w0["ckpt"]

    w1, x_cur, avg_cur = sample_one_period(ckpt, t)
    path.append(w1)
    t += tclk
    ckpt = w1["ckpt"]

    # rolling history of recent ACCEPTED (smooth, non-jump) per-period state
    # deltas, used for the smoothness check below in preference to a fixed
    # relative-to-current-norm threshold: this circuit's Newton-unknown
    # vector always contains Vin (a constant 12V rail) at a fixed index, so
    # ||x_cur|| never shrinks even when vout is still near 0V early in the
    # ramp, and a fixed-relative-to-||x_cur|| threshold would be far too
    # loose right when it matters most (the very start of soft-start, where
    # a real per-period change is a LARGE fraction of vout's own tiny
    # current value). Comparing this period's ||dx|| against the recent
    # MEDIAN ||dx|| instead is scale-robust and catches a genuine event
    # (tests/fixtures/env_r3d_softstart_redproof.cir's load step) as an
    # outlier against the smooth ramp's own recent pace, wherever in the
    # ramp it happens.
    recent_dx_norms = []

    def is_smooth(dx_norm):
        if len(recent_dx_norms) < 3:
            return True  # not enough history yet -- accept (bootstrap only)
        hist = sorted(recent_dx_norms[-8:])
        median = hist[len(hist) // 2]
        return dx_norm <= growth_tol * max(median, abs_floor)

    while t < tstop - 1e-15 and hop[0] < max_hops:
        dx = [b - a for a, b in zip(x_prev, x_cur)]
        dx_norm = norm_diff(x_cur, x_prev)
        pre_slope = avg_cur - avg_prev

        smooth = is_smooth(dx_norm)
        if not smooth:
            refusals.append(dict(t=t, reason="growth_tol", dx_norm=dx_norm,
                                  recent=list(recent_dx_norms[-8:])))

        jumped = False
        m_try = m_periods if smooth else m_min
        attempts = 0
        remaining_periods = max(0, round((tstop - t) / tclk))

        while smooth and m_try > m_min and not jumped and remaining_periods > n_check:
            attempts += 1
            m_eff = max(1, min(m_try, remaining_periods - n_check))
            land_t = t + m_eff * tclk
            x_pred = [b + m_eff * d for b, d in zip(x_cur, dx)]
            st = ck.read_checkpoint(ckpt)
            small_delta = min(st["CKTdelta"], tclk / 200.0)
            patch = dict(CKTtime=land_t, CKTrhsOld=x_pred, CKTorder=1,
                         CKTdelta=small_delta, CKTsaveDelta=small_delta,
                         CKTbreaks=[round(land_t, 15)], breakSize=1)
            patched_path = os.path.join(outdir, f"h{hop[0]:04d}_jump.in.ckpt")
            stp = dict(st)
            stp.update(patch)
            ck.write_checkpoint(patched_path, stp)

            # chain n_check REAL, fully-resolved periods from the landing
            # point -- "compare the extrapolated state against a simulated
            # check period" (card #512's own framing)
            try:
                chk_ckpt, chk_t = patched_path, land_t
                chk_avgs = []
                chk_windows = []
                for _ in range(n_check):
                    w_chk, s_chk, a_chk = sample_one_period(chk_ckpt, chk_t, kind="check")
                    chk_ckpt, chk_t = w_chk["ckpt"], chk_t + tclk
                    chk_avgs.append((s_chk, a_chk))
                    chk_windows.append(w_chk)
                danger = False
            except CouldNotRun:
                danger = True

            # Mid-check-window event scan: found by testing on tests/fixtures/
            # env_r3d_softstart_redproof.cir that a big, single-period event
            # landing INSIDE this jump's own n_check window (e.g. exactly at
            # the landing boundary, as the redproof deck's load step does when
            # M happens to land it there) can be invisible to the
            # end-to-end slope_err check below: the load step's own transient
            # decays back to a small per-period delta within 1-2 of the
            # n_check periods (the converter's fast switching dynamics
            # re-settle quickly even though the average is now on a new
            # level), so comparing only the LAST two check periods' slope
            # missed it entirely on the first version of this check (a jump
            # was accepted straight across the event, slope_err=5e-4). Fixed
            # by also running the same history-relative smoothness test
            # (is_smooth(), used for the PRE-jump growth check) across every
            # consecutive pair actually simulated in this check window --
            # not just the summary slope -- so a spike anywhere inside it is
            # caught even if the window's own endpoints look smooth.
            mid_check_event = False
            if not danger:
                chk_states = [x_cur] + [s for s, _ in chk_avgs]
                for a, b in zip(chk_states, chk_states[1:]):
                    if not is_smooth(norm_diff(b, a)):
                        danger = True
                        mid_check_event = True
                        break

            if danger:
                slope_err = float("inf")
            else:
                # Use the LAST two check periods, not the first two: found by
                # measurement (docs/envelope.md, "the landing kick") that the
                # first 1-2 periods right after a jump carry a real, measurable
                # settling transient of their own -- CKTstates[] (the per-device
                # integration history) is deliberately left unmodified by the
                # patch (see this module's own docstring), so it is briefly
                # inconsistent with the jumped-forward CKTrhsOld/CKTtime, the
                # same "state discontinuity" shape docs/checkpoint.md's patches
                # 0200/0201 already characterise at an ordinary breakpoint --
                # and comparing THAT contaminated period's slope against the
                # smooth pre-jump slope refused almost every jump even on the
                # genuinely-smooth ramp. n_check must be >=3 for this to have a
                # settled pair to compare; n_check==2 falls back to the same
                # (contaminated) pair, which is why the default below is 4.
                post_slope = chk_avgs[-1][1] - chk_avgs[-2][1] if n_check >= 2 else 0.0
                denom = max(abs(pre_slope), abs_floor)
                slope_err = abs(post_slope - pre_slope) / denom

            if not danger and slope_err <= slope_tol:
                jumps.append(dict(t_from=t, t_to=land_t, m=m_eff, slope_err=slope_err,
                                   attempts=attempts))
                path.extend(chk_windows)  # this attempt's check periods are now
                                            # part of the final trajectory
                for (s_a, _), (s_b, _) in zip(chk_avgs, chk_avgs[1:]):
                    recent_dx_norms.append(norm_diff(s_b, s_a))
                x_prev, avg_prev = chk_avgs[-2] if n_check >= 2 else (x_cur, avg_cur)
                x_cur, avg_cur = chk_avgs[-1]
                t, ckpt = chk_t, chk_ckpt
                jumped = True
            else:
                reason = "mid_check_event" if mid_check_event else "slope_tol"
                refusals.append(dict(t=t, reason=reason, slope_err=slope_err, m=m_eff,
                                      pre_slope=pre_slope,
                                      post_slope=(post_slope if slope_err != float("inf") else None)))
                m_try = max(m_min, m_try // 4)

        if not jumped:
            # passthrough: plain one-period step (no jump attempted, or every
            # retry was refused) -- re-sample normally and continue
            w_pt, x_new, avg_new = sample_one_period(ckpt, t, kind="passthrough")
            path.append(w_pt)
            recent_dx_norms.append(norm_diff(x_new, x_cur))
            x_prev, avg_prev = x_cur, avg_cur
            x_cur, avg_cur = x_new, avg_new
            t += tclk
            ckpt = w_pt["ckpt"]
        # else: a jump was accepted above -- x_prev/x_cur/avg_prev/avg_cur/t/ckpt
        # are already the last two REAL check-period samples (chk_avgs[-2:]),
        # no extra resample needed before the next outer iteration's dx.

    total_wall = sum(w["wall"] for w in windows)
    path_wall = sum(w["wall"] for w in path)
    return dict(windows=windows, path=path, jumps=jumps, refusals=refusals, t_final=t,
                total_wall=total_wall, path_wall=path_wall, outdir=outdir,
                n_periods_simulated=len(path), tclk=tclk)


# ---------------------------------------------------------------------------
# offline validation against a full plain-.tran reference run
# ---------------------------------------------------------------------------

def validate_against_reference(report, reference_path, ref_vectors, col=0,
                                avg_tol_rel=0.005, avg_tol_floor=1e-3,
                                ripple_tol_rel=0.05, ripple_tol_floor=1e-4):
    """Compare every period on the envelope run's FINAL accepted trajectory
    (`report["path"]`, one window per real-simulated period) against the
    SAME absolute time window in a full plain-`.tran` reference run's
    wrdata output (`reference_path`, same vector order as `ref_vectors`).
    Only periods that were actually simulated are checked -- honest by
    construction, not an interpolation/extrapolation over skipped periods.

    One sequential pass over the (potentially multi-million-row) reference
    file, bucketed against the path windows' own [t0,t1] spans (sorted,
    non-overlapping by construction) -- O(reference rows + path windows),
    never loads the whole reference file into memory at once.

    Returns dict(ok, n_checked, worst_avg_rel, worst_avg_detail,
    worst_ripple_rel, worst_ripple_detail, per_window=[...]).
    """
    spans = sorted(((w["t0"], w["t1"], w) for w in report["path"]), key=lambda s: s[0])
    ncols = len(ref_vectors)
    idx = 0
    n = len(spans)
    # accumulator per active span: trapz area, last (t,v), max, min
    acc = None
    results = {}

    def finalize(t0, t1, w):
        if acc is None or acc["last_t"] is None or acc["first_t"] is None or \
                acc["last_t"] <= acc["first_t"]:
            return None
        ref_avg = acc["area"] / (acc["last_t"] - acc["first_t"])
        ref_max, ref_min = acc["max"], acc["min"]
        my_avg, my_max, my_min = _period_stats(w["rows"], t0, t1, col)
        avg_tol = max(avg_tol_rel * abs(ref_avg), avg_tol_floor)
        my_ripple, ref_ripple = my_max - my_min, ref_max - ref_min
        ripple_tol = max(ripple_tol_rel * abs(ref_ripple), ripple_tol_floor)
        return dict(t0=t0, t1=t1, ref_avg=ref_avg, my_avg=my_avg,
                    avg_err=abs(my_avg - ref_avg), avg_tol=avg_tol,
                    ref_ripple=ref_ripple, my_ripple=my_ripple,
                    ripple_err=abs(my_ripple - ref_ripple), ripple_tol=ripple_tol)

    with open(reference_path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.split()
            if len(parts) != 1 + ncols:
                continue
            try:
                t = float(parts[0])
                v = float(parts[1 + col])
            except ValueError:
                continue
            while idx < n and t > spans[idx][1] + 1e-15:
                t0, t1, w = spans[idx]
                r = finalize(t0, t1, w)
                if r is not None:
                    results[t0] = r
                idx += 1
                acc = None
            if idx >= n:
                break
            t0, t1, w = spans[idx]
            if t0 - 1e-15 <= t <= t1 + 1e-15:
                if acc is None:
                    acc = dict(area=0.0, last_t=None, last_v=None,
                               first_t=None, max=v, min=v)
                if acc["last_t"] is not None:
                    acc["area"] += (acc["last_v"] + v) * 0.5 * (t - acc["last_t"])
                else:
                    acc["first_t"] = t
                acc["last_t"], acc["last_v"] = t, v
                acc["max"] = max(acc["max"], v)
                acc["min"] = min(acc["min"], v)
        if idx < n and acc is not None:
            t0, t1, w = spans[idx]
            r = finalize(t0, t1, w)
            if r is not None:
                results[t0] = r

    per_window = [results[t0] for t0, t1, w in spans if t0 in results]
    if not per_window:
        return dict(ok=False, n_checked=0, worst_avg_rel=None, worst_ripple_rel=None,
                     per_window=[])
    worst_avg = max(per_window, key=lambda r: r["avg_err"] / max(r["avg_tol"], 1e-30))
    worst_ripple = max(per_window, key=lambda r: r["ripple_err"] / max(r["ripple_tol"], 1e-30))
    ok = all(r["avg_err"] <= r["avg_tol"] for r in per_window) and \
        all(r["ripple_err"] <= r["ripple_tol"] for r in per_window)
    return dict(ok=ok, n_checked=len(per_window),
                worst_avg_rel=worst_avg["avg_err"] / max(abs(worst_avg["ref_avg"]), 1e-30),
                worst_avg_detail=worst_avg,
                worst_ripple_rel=worst_ripple["ripple_err"] / max(abs(worst_ripple["ref_ripple"]), 1e-30),
                worst_ripple_detail=worst_ripple,
                per_window=per_window)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("--fguess", type=float, required=True, help="switching frequency (Hz)")
    ap.add_argument("--tstop", type=float, required=True)
    ap.add_argument("--tstep", type=float, default=2e-9)
    ap.add_argument("--m-periods", type=int, default=10)
    ap.add_argument("--n-check", type=int, default=4)
    ap.add_argument("--growth-tol", type=float, default=5.0,
                     help="refuse a jump when this period's state delta exceeds "
                          "this many times the recent median delta (default 5x)")
    ap.add_argument("--slope-tol", type=float, default=3.0)
    ap.add_argument("--abs-floor", type=float, default=1e-3,
                     help="volts: below this, the pre-jump slope is treated as "
                          "'too small to judge tightly' (avoids dividing by "
                          "near-zero right at the very start of a ramp)")
    ap.add_argument("--vector", action="append", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--ngspice", default=None)
    args = ap.parse_args(argv)

    vectors = args.vector or ["v(vout)", "i(L1)"]
    ngspice = args.ngspice or NGSPICE
    name = os.path.splitext(os.path.basename(args.deck))[0]
    outdir = args.out or os.path.join(DATA_ROOT, "envelope", name)
    tclk = 1.0 / args.fguess

    t0 = time.time()
    report = envelope_follow(args.deck, outdir, tclk, args.tstop, args.tstep, vectors,
                              m_periods=args.m_periods, n_check=args.n_check,
                              growth_tol=args.growth_tol, slope_tol=args.slope_tol,
                              abs_floor=args.abs_floor,
                              ngspice=ngspice)
    wall = time.time() - t0

    print(f"envelope_follow: t_final={report['t_final']:.6g} target={args.tstop:.6g} "
          f"windows={len(report['windows'])} jumps={len(report['jumps'])} "
          f"refusals={len(report['refusals'])}")
    print(f"  wall (sum of window subprocess times) = {report['total_wall']:.3f}s, "
          f"driver wall = {wall:.3f}s")
    for j in report["jumps"][:10]:
        print(f"  jump t={j['t_from']:.6g} -> {j['t_to']:.6g} (M={j['m']}, "
              f"slope_err={j['slope_err']:.3g}, attempts={j['attempts']})")
    if len(report["jumps"]) > 10:
        print(f"  ... {len(report['jumps']) - 10} more jumps")
    for r in report["refusals"][:10]:
        print(f"  refusal at t={r['t']:.6g}: {r['reason']} {r}")
    return 0 if report["t_final"] >= args.tstop - tclk else 1


if __name__ == "__main__":
    sys.exit(main())
