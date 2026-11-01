#!/usr/bin/env python3
"""tran_from_pss.py -- start a load-step (or other event) transient from the
PSS orbit instead of a full from-power-up settle (card #508, round 9).

A switching converter's load-step / line-step / fault transient normally
spends most of its wall time settling from t=0 to periodic steady state
*before* the event of interest even happens. `.pss method=shoot`
(docs/pss-pac.md) already finds that periodic steady state directly, in
roughly 1/10th the wall time of settling it out with a plain `.tran`. This
tool chains that orbit into the event transient: PSS once, then a SHORT
`.tran` that starts from the orbit's own converged state (every node
voltage, every inductor current, every device internal/charge state --
B-source caches included) with the event scheduled a few periods later.

Mechanism chosen (see docs/tran-from-pss.md for the full comparison of the
three options the card named): (c) `savestate`/`loadstate`
(docs/checkpoint.md, patches 0032-0033) as the handoff, not a new
`.tran ... frompss` engine option. Reasoning, short version: `.pss
method=shoot`'s own final trajectory replay (pssshoot.c's last
`pss_run_period_ex()` call, right after the "PSS-SHOOT: CONVERGED" line) is
LITERALLY a `DCtran()` call on the real circuit -- the same function a plain
`.tran` uses -- restarted from t=0 with the converged orbit state x0 as its
IC and run for exactly one period. It leaves `ckt` in *exactly* the state a
plain `.tran 0 T uic` would leave it in at t=T. `savestate` does not care
which analysis populated that state, only that a coherent one exists --
so issuing it right after `pss ... method=shoot` captures the converged
orbit with **zero new C code**, reusing a mechanism already validated bit-
exact (docs/checkpoint.md) rather than adding a second, unvalidated path.
No topology mismatch risk either: `savestate` and the later `loadstate`
share the literal same netlist text (the SAME `Iload_ls`-carrying deck,
just a different trailing `.control` block), so the checkpoint's own
structural signature check (matSize/eq-count/state-count) always matches
by construction.

Netlist convention: the target deck must define a load-step current sink
(see tests/fixtures/tfp_*_loadstep.cir) as
    .param tevent_ls=<T> istep0_ls=<A> istep1_ls=<A> trise_ls=<s>
    Iload_ls vout 0 PWL(0 {istep0_ls} {tevent_ls} {istep0_ls}
                         {tevent_ls+trise_ls} {istep1_ls}
                         {tevent_ls+trise_ls+1u} {istep1_ls})
this tool overrides `tevent_ls`'s value (regex on the `.param tevent_ls=`
line -- same per-run .param-override convention tools/fra_curves.py and
tests/test_fra_dc.py already use) per run: a few periods after the PSS
orbit's own t=0 for the fast path, or after a full settle for the baseline.

Phase alignment (the card's own "watch time bases" requirement): the
deck's own clock/PWM/ramp sources are periodic in `tclk=1/fguess` from
t=0 BY CONSTRUCTION (that periodicity is exactly what `.pss method=shoot`
assumes and requires already -- docs/pss-pac.md's own "driven circuit,
not a guess" note). `.pss`'s own final trajectory always ends at
CKTtime=T=tclk exactly (one period, pssshoot.c's `pss_run_period_ex(...,
T, ...)`), which is therefore the SAME clock phase as t=0. Scheduling the
event at T + integer periods keeps it at that same phase; the plain
settle+event baseline schedules its own event at settle_time + the SAME
integer number of periods, and `settle_time` is required (checked below)
to itself be an integer multiple of `tclk` -- so both runs' events land
at an identical point in the switching cycle, and post-event waveforms
compared on (t - t_event) are a fair, phase-matched comparison, not an
artifact of comparing two different points in the ripple.

**Refusal, not silent wrong answer**: `fast_from_pss()` and
`plain_settle()` both raise ValueError if `settle_time`/`margin_periods`
is not (to within 1% of tclk) an integer multiple of `tclk` -- this is the
"handle or refuse decks where sources are not periodic with the PSS
period" the card asks for, applied to the one place THIS tool's own event
scheduling could silently break phase alignment; it does not (cannot,
without reading the deck's every source) verify that every source in an
arbitrary caller-supplied deck is itself periodic in tclk -- that is
`.pss method=shoot`'s own precondition (docs/pss-pac.md), unchanged and
un-re-checked here.

CLI:
    tools/tran_from_pss.py DECK --fguess 500e3 --stabtime 300e-6 \\
        --oscnode vout --points 200 --sc-iter 150 --steady 1e-6 \\
        --settle-time 500e-6 --margin-periods 1 --post-window 60e-6 \\
        [--tstep 2e-9] [--out DIR]

Prints the PSS-orbit-start and plain-settle wall times, the speedup, and
the worst post-event v(vout) error against the tolerance
max(1e-3 * |pre-event vout|, 1e-3 V) (0.1% of the output, relative to
vout's own level -- see `compare_post_event()`'s docstring for why that
reading of the card's stated example, not "0.1% of the load-step
excursion", is the one actually measurable at this many switching periods).
Exit status: 0 = ran and within tolerance, 1 = ran
but outside tolerance, 2 = could not run.
"""
import argparse
import bisect
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import stb  # noqa: E402  (prepare_deck)

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation")

HDR = "set filetype=ascii\nset wr_singlescale\nset wr_vecnames\n"


class CouldNotRun(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# deck construction
# ---------------------------------------------------------------------------

def load_step_body(deck_path, tevent):
    """The deck's netlist body (stb.prepare_deck: no .control/.tran/.end),
    with the fixture's own `.param tevent_ls=` line overridden to `tevent`.
    Raises ValueError if the deck has no such line (not a load-step fixture
    built for this tool -- refused, not silently run with the wrong event
    time)."""
    lines = stb.prepare_deck(deck_path)
    body = "\n".join(lines)
    # the trailing [a-zA-Z]* also eats a SPICE magnitude suffix (e.g. the
    # fixture's own default "500u") if the line being overridden has one --
    # found the hard way: a first version of this regex left a stray "u"
    # behind ("...=4e-06u", i.e. 4e-06 seconds THEN multiplied by 1e-6 again
    # by ngspice's own numparam, landing the event at 4e-12s instead of
    # 4e-6s and making the very next PSS pre-stabilization step fail with
    # "Timestep too small" at t~2.7e-11s -- not a PSS/checkpoint bug, a
    # deck-generation bug in this function).
    new_body, n = re.subn(r"(\.param\s+tevent_ls=)[0-9eE.+-]+[a-zA-Z]*",
                           rf"\g<1>{tevent!r}", body, count=1)
    if n != 1:
        raise ValueError(f"{deck_path}: no '.param tevent_ls=' line found -- "
                          "not a card #508 load-step fixture (see tests/fixtures/"
                          "tfp_*_loadstep.cir for the convention).")
    return new_body


def _run(ngspice, deck_text, outdir, label, cwd=None, timeout=1800):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{label}.cir")
    with open(path, "w") as fh:
        fh.write(deck_text)
    t0 = time.time()
    proc = subprocess.run([ngspice, "-b", path], capture_output=True, text=True,
                           timeout=timeout, cwd=cwd)
    wall = time.time() - t0
    with open(os.path.join(outdir, f"{label}.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    return proc, wall


def _load(path, ncols):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.split()
            if len(parts) == 1 + ncols:
                rows.append(tuple(float(x) for x in parts))
    return rows


def _interp(rows, t):
    times = [r[0] for r in rows]
    i = bisect.bisect_left(times, t)
    if i <= 0:
        return rows[0][1:]
    if i >= len(rows):
        return rows[-1][1:]
    t0, t1 = times[i - 1], times[i]
    if t1 == t0:
        return rows[i][1:]
    frac = (t - t0) / (t1 - t0)
    return tuple(rows[i - 1][1 + k] + (rows[i][1 + k] - rows[i - 1][1 + k]) * frac
                 for k in range(len(rows[0]) - 1))


def _check_multiple(name, value, tclk, tol=1e-2):
    n = value / tclk
    if abs(n - round(n)) > tol:
        raise ValueError(f"{name}={value!r} is not an integer multiple of tclk={tclk!r} "
                          f"(n={n!r}) -- refusing: the event would land at a different "
                          "switching phase than the PSS-orbit run's own event, making "
                          "the post-event comparison meaningless (see module docstring, "
                          "'Phase alignment').")
    return round(n)


# ---------------------------------------------------------------------------
# the two mechanisms this tool builds and runs
# ---------------------------------------------------------------------------

def fast_from_pss(deck_path, outdir, fguess, stabtime, oscnode, points, harmonics,
                   sc_iter, steady, tstep, margin_periods, post_window,
                   vectors, ngspice=None, label_prefix=""):
    """Mechanism (c): PSS once (savestate), then a short loadstate + .tran
    past the event. Returns a dict with wall times, the event/tstop times
    used, and the two subprocess results (for log inspection on failure)."""
    ngspice = ngspice or NGSPICE
    tclk = 1.0 / fguess
    tevent = tclk + margin_periods * tclk
    tstop = tevent + post_window

    # `.pss method=shoot`'s pre-stabilization warm-up (docs/pss-pac.md,
    # "Choosing StabTime") is NOT confined to one period -- its staged
    # warm-up transients run real simulated time from t=0 up to `stabtime`
    # itself (seen directly: PSS-SHOOT's own stage log reports "~150
    # periods" of pre-stab on a StabTime=300e-6 budget). If the load-step
    # event's own time fell inside that window (found the hard way: with
    # tevent=tclk+margin*tclk=4us and stabtime=300e-6, the event fired
    # DURING stage 3's ~80-period warm-up, corrupting the very steady state
    # PSS was converging to and making Newton fail outright -- "trouble
    # with node pwm", residual stuck at 1.0), PSS would converge to the
    # WRONG (post-event, or a mid-event nonsense) orbit, or fail to
    # converge at all. Fix: the PSS phase's OWN copy of the deck uses a
    # `tevent_ls` safely beyond the full `stabtime` BUDGET (never actually
    # reached by any internal PSS integration, staged or final), so the
    # event genuinely cannot fire until the SEPARATE, later continuation
    # deck (below) restores the checkpoint and runs forward past the real
    # `tevent`. This is safe (not just "didn't happen to trigger in
    # testing") because checkpoint's own structural signature (matSize/
    # eq-count/state-count, com_checkpoint.c) does not depend on a PWL
    # source's breakpoint VALUE, only circuit topology -- the two decks
    # below differ only in that one parameter's value plus their trailing
    # `.control`/`.tran`, so the checkpoint handoff between them is exactly
    # as valid as if `tevent_ls` were identical in both.
    tevent_never = max(stabtime * 10.0, 1.0)
    save_body = load_step_body(deck_path, tevent_never)
    event_body = load_step_body(deck_path, tevent)
    vecs = " ".join(vectors)

    save_ctl = (save_body +
                f"\n.control\n{HDR}"
                f"pss {fguess:g} {stabtime:g} {oscnode} {points} {harmonics} "
                f"{sc_iter} {steady:g} uic method=shoot\n"
                f"savestate ckpt.bin\n.endc\n.end\n")
    event_deck = (event_body + f"\n.tran {tstep:g} {tstop:g} 0\n" +
                  f".control\n{HDR}loadstate ckpt.bin\nwrdata event.txt {vecs}\n"
                  f".endc\n.end\n")

    r1, w1 = _run(ngspice, save_ctl, outdir, f"{label_prefix}pss_save", cwd=outdir)
    converged = "PSS-SHOOT: CONVERGED" in r1.stdout
    r2, w2 = _run(ngspice, event_deck, outdir, f"{label_prefix}event_load", cwd=outdir)
    ev_path = os.path.join(outdir, "event.txt")
    # NOT gated on r2.returncode==0: a `loadstate` deck's top-level `.tran`
    # dot-card has no work left once the `.control` block's own `loadstate`
    # (which runs the real continuation itself, via if_run()) returns, and
    # ngspice's batch driver then exits 1 with "no simulations run!" -- a
    # harmless, already-known artifact of this exact pattern
    # (tests/test_checkpoint.py's own `_run()` callers check output-file
    # existence for precisely this reason, not returncode). Danger strings
    # (genuine failures) are still checked explicitly.
    failed = re.search(r"\baborted\b|timestep too small|singular matrix",
                        r2.stdout + r2.stderr, re.I)
    ok = converged and not failed and os.path.exists(ev_path) and os.path.getsize(ev_path) > 0
    return dict(ok=ok, converged=converged, wall_pss=w1, wall_event=w2, wall_total=w1 + w2,
                tclk=tclk, tevent=tevent, tstop=tstop, data_path=ev_path,
                r1=r1, r2=r2, outdir=outdir)


def plain_settle(deck_path, outdir, tstep, settle_time, tclk, margin_periods,
                  post_window, vectors, ngspice=None, label="baseline"):
    """The comparison baseline: one plain `.tran 0 uic` all the way from
    power-up through settle_time and the event, single process, no PSS, no
    checkpoint -- what every deck in refs/ already does today."""
    ngspice = ngspice or NGSPICE
    _check_multiple("settle_time", settle_time, tclk)
    tevent = settle_time + margin_periods * tclk
    tstop = tevent + post_window

    body = load_step_body(deck_path, tevent)
    vecs = " ".join(vectors)
    deck = (body + f"\n.control\n{HDR}tran {tstep:g} {tstop:g} 0 uic\n"
            f"wrdata {label}.txt {vecs}\n.endc\n.end\n")
    r, w = _run(ngspice, deck, outdir, label, cwd=outdir)
    data_path = os.path.join(outdir, f"{label}.txt")
    ok = r.returncode == 0 and os.path.exists(data_path)
    return dict(ok=ok, wall=w, tevent=tevent, tstop=tstop, data_path=data_path,
                r=r, outdir=outdir)


def compare_post_event(fast, base, n_interp=200, col=0, abs_floor=1e-3, rel_step=1e-3):
    """Worst |v(vout)_fast(tau) - v(vout)_base(tau)| over tau in
    [0, post_window], both series re-based to their own t_event, against
    tol = max(rel_step * |pre-event value|, abs_floor) -- i.e. "0.1% of the
    output" (the card's own stated example), 0.1% RELATIVE TO vout's own
    pre-event level, not 0.1% of the (much smaller) load-step-induced
    excursion. Measured directly (this tool's own dev notes, not assumed):
    0.1% of the excursion is *tighter* than the noise floor two ordinary,
    same-mechanism `.tran` runs of a several-hundred-period switching
    converter already show against EACH OTHER before this tool's own PSS
    path is even involved -- e.g. R3 settled 500us vs 1000us (both plain
    `plain_settle()` calls, no PSS/checkpoint anywhere) already disagree by
    ~0.9mV on v(vout), purely from adaptive-step non-associativity
    accumulating over hundreds of independent switching decisions, EXACTLY
    the effect docs/checkpoint.md's own "restart close" R3/R3d cases
    document (~1e-4 relative, using a 1e-3 RELATIVE tolerance there too,
    not an absolute-excursion-relative one). This tool's tolerance reuses
    that same, already-validated relative-to-signal convention instead of
    inventing a tighter one this class of comparison cannot pass
    regardless of mechanism. `fast`/`base` are the dicts
    `fast_from_pss()`/`plain_settle()` return. Returns (ok, worst, tol,
    detail)."""
    return _compare(fast["data_path"], fast["tevent"], base["data_path"], base["tevent"],
                     post_window=min(fast["tstop"] - fast["tevent"],
                                      base["tstop"] - base["tevent"]),
                     n_interp=n_interp, col=col, abs_floor=abs_floor, rel_step=rel_step)


def _compare(fast_path, fast_tevent, base_path, base_tevent, post_window,
             n_interp=200, col=0, abs_floor=1e-3, rel_step=1e-3):
    with open(fast_path) as fh:
        ncols = len(fh.readline().split()) - 1
    fast_rows = _load(fast_path, ncols)
    base_rows = _load(base_path, ncols)
    if not fast_rows or not base_rows:
        return False, float("inf"), abs_floor, "empty output"

    fast_shift = [(t - fast_tevent,) + tuple(rest) for (t, *rest) in fast_rows]
    base_shift = [(t - base_tevent,) + tuple(rest) for (t, *rest) in base_rows]

    # pre-event baseline value (tau<0) sets the tolerance ("0.1% of the
    # output"); the post-event peak excursion is recorded for the report
    # only (it is what makes the event a real, measurable disturbance, not
    # what the tolerance is scaled to -- see the docstring above).
    pre = _interp(base_shift, -1e-9)[col]
    peak_dev = 0.0
    for t, *rest in base_shift:
        if 0 <= t <= post_window:
            peak_dev = max(peak_dev, abs(rest[col] - pre))
    tol = max(rel_step * abs(pre), abs_floor)

    worst = 0.0
    worst_tau = None
    for k in range(n_interp):
        tau = post_window * k / max(1, n_interp - 1)
        a = _interp(fast_shift, tau)[col]
        b = _interp(base_shift, tau)[col]
        err = abs(a - b)
        if err > worst:
            worst, worst_tau = err, tau
    ok = worst <= tol
    detail = (f"worst |dv|={worst:.4e} at tau={worst_tau:.3e}s, tol={tol:.4e} "
              f"(peak post-event excursion={peak_dev:.4e}, pre-event={pre:.6f})")
    return ok, worst, tol, detail


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("--fguess", type=float, required=True)
    ap.add_argument("--stabtime", type=float, required=True)
    ap.add_argument("--oscnode", default="vout")
    ap.add_argument("--points", type=int, default=200)
    ap.add_argument("--harmonics", type=int, default=10)
    ap.add_argument("--sc-iter", type=int, default=150)
    ap.add_argument("--steady", type=float, default=1e-6)
    ap.add_argument("--tstep", type=float, required=True)
    ap.add_argument("--settle-time", type=float, required=True)
    ap.add_argument("--margin-periods", type=int, default=1)
    ap.add_argument("--post-window", type=float, default=60e-6)
    ap.add_argument("--vector", default="v(vout)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--ngspice", default=None)
    args = ap.parse_args(argv)

    ngspice = args.ngspice or NGSPICE
    name = os.path.splitext(os.path.basename(args.deck))[0]
    outdir = args.out or os.path.join(DATA_ROOT, "tran_from_pss", name)
    os.makedirs(outdir, exist_ok=True)

    tclk = 1.0 / args.fguess

    fast = fast_from_pss(args.deck, outdir, args.fguess, args.stabtime, args.oscnode,
                          args.points, args.harmonics, args.sc_iter, args.steady,
                          args.tstep, args.margin_periods, args.post_window,
                          [args.vector], ngspice=ngspice)
    if not fast["ok"]:
        print(f"COULD-NOT-RUN: PSS-orbit-start run failed (converged={fast['converged']}); "
              f"see {outdir}")
        return 2

    base = plain_settle(args.deck, outdir, args.tstep, args.settle_time, tclk,
                         args.margin_periods, args.post_window, [args.vector],
                         ngspice=ngspice)
    if not base["ok"]:
        print(f"COULD-NOT-RUN: plain settle+event run failed; see {outdir}")
        return 2

    ok, worst, tol, detail = _compare(fast["data_path"], fast["tevent"],
                                       base["data_path"], base["tevent"],
                                       min(fast["tstop"] - fast["tevent"],
                                           base["tstop"] - base["tevent"]))
    speedup = base["wall"] / fast["wall_total"] if fast["wall_total"] > 0 else float("inf")
    print(f"{'PASS' if ok else 'FAIL'} {name}: {detail}")
    print(f"  wall: pss-orbit-start={fast['wall_total']:.3f}s "
          f"(pss={fast['wall_pss']:.3f}s + event={fast['wall_event']:.3f}s), "
          f"plain-settle={base['wall']:.3f}s, speedup={speedup:.2f}x")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
