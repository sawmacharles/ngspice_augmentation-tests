#!/usr/bin/env python3
"""Validation for tools/tran_from_pss.py (card #508, transient acceleration
round 9): "start the transient from the PSS orbit" instead of a plain
from-power-up settle before a load-step event.

See docs/tran-from-pss.md for the full design/decision writeup (why
mechanism (c), savestate/loadstate, was chosen over (a) .ic/uic and (b) a
new `.tran ... frompss` engine option -- short version: (c) needs ZERO new
ngspice C code, because `.pss method=shoot`'s own final trajectory replay
is literally a `DCtran()` call, the same function a plain `.tran` uses).

This file checks, on all four load-step fixtures
(tests/fixtures/tfp_{r3,r3d,r4,r4d}_loadstep.cir):

1. Refusal behaviour: a non-integer-multiple-of-tclk `settle_time` is
   refused (ValueError), and a deck with no `.param tevent_ls=` line is
   refused (not silently run with the wrong/no event).
2. PSS-orbit-start runs cleanly (converges, produces output) and its
   post-event v(vout) agrees with a plain settle+event baseline to within
   max(1e-3 * |pre-event vout|, 1e-3 V) -- see docs/tran-from-pss.md's
   "Accuracy" section for why this tolerance (0.1% of the OUTPUT's level,
   not of the load-step excursion) is the one this class of comparison can
   actually pass, and the measured settle times used per deck.

Runtime note: this file runs each deck's PSS-orbit-start and plain-settle
comparison ONCE (not the 5-rep interleaved A/B the round 9 report's own
speed table uses -- that is a timing claim, not a correctness gate, and is
reported in research/transient-acceleration.md instead). R4d's own plain
baseline alone is ~18-19s; total runtime here is ~30s.

    NGSPICE=<fork build> python3 tests/test_tran_from_pss.py
    exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN

Needs a fork build with `.pss method=shoot` and savestate/loadstate (this
card's own `aaba50ba3f49-accel9-1`, or `s55`/`s54` -- unchanged, no new C
this round). Default NGSPICE tried in that order if NGSPICE is unset.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import tran_from_pss as tfp  # noqa: E402

FIXTURES = os.path.join(ROOT, "tests", "fixtures")
OUT = os.path.join(tfp.DATA_ROOT, "tests_tran_from_pss")

_CANDIDATES = [
    os.environ.get("NGSPICE"),
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-accel9-1/bin/ngspice-aug"),
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s55/bin/ngspice-aug"),
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s54/bin/ngspice-aug"),
]
NGSPICE = next((p for p in _CANDIDATES if p and os.path.exists(p) and os.access(p, os.X_OK)), None)

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))
    print(f"{state:<14}{name:<78}{detail}")


DECKS = {
    "r3":  dict(deck=os.path.join(FIXTURES, "tfp_r3_loadstep.cir"),
                stabtime=300e-6, tstep=10e-9, settle=500e-6),
    "r3d": dict(deck=os.path.join(FIXTURES, "tfp_r3d_loadstep.cir"),
                stabtime=700e-6, tstep=2e-9, settle=510e-6),  # 2026-09-30 (card #519): was 500e-6; the r3d shooter no longer converges at 500e-6 with the corrected Goa polarity (converges at 700e-6)
    "r4":  dict(deck=os.path.join(FIXTURES, "tfp_r4_loadstep.cir"),
                stabtime=1500e-6, tstep=10e-9, settle=6000e-6),
    "r4d": dict(deck=os.path.join(FIXTURES, "tfp_r4d_loadstep.cir"),
                stabtime=1500e-6, tstep=2e-9, settle=6688e-6),
}


def check_refusal_non_multiple_settle():
    name = "refuse: settle_time not an integer multiple of tclk"
    tclk = 1.0 / 500e3
    try:
        tfp._check_multiple("settle_time", 500.3e-6, tclk)
        record(name, False, "expected ValueError, none raised")
    except ValueError as e:
        record(name, True, f"raised as expected: {e}")


def check_refusal_no_tevent_ls():
    name = "refuse: deck with no '.param tevent_ls=' line"
    plain_deck = os.path.join(ROOT, "refs", "r3_vm_buck.cir")
    try:
        tfp.load_step_body(plain_deck, 4e-6)
        record(name, False, "expected ValueError, none raised")
    except ValueError as e:
        record(name, True, f"raised as expected: {e}")


def check_deck(name, spec):
    if NGSPICE is None:
        record(name, None, "no fork build found (set NGSPICE); see module docstring")
        return
    tclk = 1.0 / 500e3
    outdir = os.path.join(OUT, name)

    fast = tfp.fast_from_pss(spec["deck"], outdir, 500e3, spec["stabtime"], "vout",
                              200, 10, 150, 1e-6, spec["tstep"], 1, 60e-6, ["v(vout)"],
                              ngspice=NGSPICE, label_prefix="fast_")
    if not fast["ok"]:
        record(name, None, f"PSS-orbit-start run failed (converged={fast['converged']}); "
                            f"see {outdir}")
        return

    base = tfp.plain_settle(spec["deck"], outdir, spec["tstep"], spec["settle"], tclk,
                             1, 60e-6, ["v(vout)"], ngspice=NGSPICE, label="base")
    if not base["ok"]:
        record(name, None, f"plain settle+event run failed; see {outdir}")
        return

    ok, worst, tol, detail = tfp.compare_post_event(fast, base)
    speedup = base["wall"] / fast["wall_total"] if fast["wall_total"] > 0 else float("inf")
    record(name, ok, f"{detail}; wall pss-orbit-start={fast['wall_total']:.3f}s "
                      f"(pss={fast['wall_pss']:.3f}s+event={fast['wall_event']:.3f}s), "
                      f"plain-settle={base['wall']:.3f}s, speedup={speedup:.2f}x "
                      "(single rep -- see research/transient-acceleration.md round 9 "
                      "for the n=5 interleaved A/B timing table)")


def main():
    check_refusal_non_multiple_settle()
    check_refusal_no_tevent_ls()
    for name, spec in DECKS.items():
        check_deck(f"post-event v(vout) accuracy: {name}", spec)

    n_bad = sum(1 for _, s, _ in results if s == "FAIL")
    n_cnr = sum(1 for _, s, _ in results if s == "COULD-NOT-RUN")
    print(f"\n{len(results) - n_bad - n_cnr}/{len(results)} pass"
          + (f", {n_bad} FAIL" if n_bad else "") + (f", {n_cnr} COULD-NOT-RUN" if n_cnr else ""))
    if n_bad:
        return 1
    if n_cnr:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
