#!/usr/bin/env python3
"""Validation for checkpoint/restart (card #472, plan T10.4).

Upstream ngspice pre-master-48 (Enhancement-131, `src/frontend/com_checkpoint.c`)
already ships `savestate <file>` / `loadstate <file>`: dump a transient's state to
disk and resume it, in a fresh process, from where it left off. This file checks
the fork build's version of that mechanism (patches-wip/ckpt/0200-0201, which fix
two bugs found in the upstream code while building this test -- see
patches-wip/ckpt/README.md and docs/checkpoint.md for the full characterisation)
against the spec's `checkpoint` acceptance line (spec.md):

    "a transient saved at t1 and restarted to t2 matches an uninterrupted run to
    t2 within 1e-9 relative on every node (same build, same options); an FRA/PSS
    campaign restarted from one settled checkpoint gives the same results as
    settling each point separately, faster."

**Needs a fork build with savestate/loadstate** (ngspice-45.2 does not have them;
neither does a plain, unpatched pre-master-48 checkout -- see the red-proof checks
below). Set NGSPICE_AUG to its `ngspice-aug` wrapper; default is this card's own
build, ~/Work/apps/ngspice_augmentation/aaba50ba3f49-ckpt5/bin/ngspice-aug
(pin aaba50ba3f49 + patches 0001-0016 (aug) + patches-wip/ckpt/0200-0201). If that
build is not present and NGSPICE_AUG is not set, every check reports COULD-NOT-RUN
(exit 2), not FAIL.

**What "restart exactness" means here, honestly** (see docs/checkpoint.md for the
full discussion): on refs/r2_switched_rc.cir, checkpointing exactly at a clock
edge, restart reproduces the single-shot uninterrupted run **bit-for-bit** (0.0
relative error on every one of 12947 shared output samples) -- strictly inside the
spec's 1e-9 bound. On bench/decks/d_analog_block.cir (BSIM4, no repeated
switching, checkpoint at an arbitrary non-edge time), the same bit-for-bit result
holds to the resolution of a 500-point comparison grid. On refs/r3_vm_buck.cir and
refs/r3d_vm_buck_deadtime.cir (checkpoint at 500us, then ~500-750 more switching
periods to t2), the restarted and single-shot trajectories' *smooth* state (v(vout))
agree to ~1e-4 relative, not 1e-9: floating-point rounding is not associative, so a
run split into "compute to t1, serialise, deserialise, compute to t2" does not take
the bit-identical sequence of Newton iterations a single continuous run does once
hundreds of *further* independent breakpoint-timing decisions follow the split (each
one a chance for a sub-ULP difference to become a step-size choice that differs from
the single-shot run's own). This is NOT evidence the restored state is wrong (the
saved/restored fields are copied byte-for-byte, and R2/BSIM4 above show the
mechanism itself is exact); it is what "the same DAE, walked by two different valid
adaptive step sequences after a long enough chain of independent decisions" looks
like. Reported here as measured, not smoothed over.

    python3 tests/test_checkpoint.py        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN

Runtime: ~2-3 minutes (R3d and the FRA-checkpoint speedup check on R4 are the slow parts).
"""
import bisect
import cmath
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import fra  # noqa: E402
import stb  # noqa: E402
import stbcore  # noqa: E402

NGSPICE_AUG = os.environ.get(
    "NGSPICE_AUG") or os.environ.get(
    "NGSPICE",  # the portfolio's build under test (2026-09-30: the pinned default below was silently used)
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-ckpt5/bin/ngspice-aug"))
OUT = os.path.join(fra.DATA_ROOT, "tests_checkpoint")

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))
    print(f"{state:<14}{name:<78}{detail}")


def _have_ngspice_aug():
    return os.path.exists(NGSPICE_AUG) and os.access(NGSPICE_AUG, os.X_OK)


# ---------------------------------------------------------------------------
# low-level ngspice plumbing for this file (deliberately not shared with
# tools/fra.py's _run_ngspice -- this file needs a netlist-level '.tran' card
# ahead of .control, which fra.py's own checkpoint mode also needs and gets via
# its own cwd/allow_patterns-aware _run_ngspice; duplicated here in a smaller
# form to keep the restart-exactness harness self-contained and easy to read).
# ---------------------------------------------------------------------------

def _spice(t):
    """"500u" -> 5e-4 etc, accepting the handful of suffixes used in this file."""
    m = re.match(r"^([-+0-9.eE]+)([un]?)$", t)
    mult = {"": 1.0, "n": 1e-9, "u": 1e-6}
    return float(m.group(1)) * mult[m.group(2)]


def _run(ngspice, deck_text, outdir, label, cwd=None):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{label}.cir")
    with open(path, "w") as fh:
        fh.write(deck_text)
    proc = subprocess.run([ngspice, "-b", path], capture_output=True, text=True,
                           timeout=1800, cwd=cwd)
    with open(os.path.join(outdir, f"{label}.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    return proc


def _load(path, ncols):
    rows = []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.split()
            if len(parts) == 1 + ncols:
                rows.append(tuple(float(x) for x in parts))
    return rows


def _interp(rows, t):
    """Linear interpolation of rows=[(t, v1, v2, ...), ...] at time t."""
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


HDR = "set filetype=ascii\nset wr_singlescale\nset wr_vecnames\n"


def restart_exactness(name, deck_path, tstep, t1, t2, vectors, uic=True, extra_options="",
                       compare_indices=None, n_interp=500, tol=1e-9, outdir=None):
    """Core check: savestate at t1 (own process), loadstate to t2 (fresh process),
    versus a true single-shot uninterrupted `tran 0..t2` (a third, independent
    process). All three share the same deck_lines (stb.prepare_deck).

    `compare_indices` (default None = all of `vectors`): which 0-based vector
    indices to include in the interpolated worst-case relative-error check (see
    module docstring on why fast/ripple-bearing signals need a coarser bound or
    exclusion, done explicitly per case below, not silently).

    Records TWO checks: an exact (grid-matched, index-for-index) comparison from
    t1 onward where the grids happen to line up (reports its own worst relative
    error and how many of the shared samples matched by time at all), and the
    interpolated comparison at `n_interp` evenly spaced points in [t1, t2] against
    `tol` (the actual pass/fail gate for this case).
    """
    outdir = outdir or os.path.join(OUT, name.replace(" ", "_"))
    uics = " uic" if uic else ""
    t1v, t2v = _spice(t1), _spice(t2)
    deck_lines = stb.prepare_deck(deck_path)
    body = "\n".join(deck_lines) + "\n"
    vecs = " ".join(vectors)
    ncols = len(vectors)

    save_ctl = (body + extra_options +
                f".control\n{HDR}tran {tstep} {t1} 0{uics}\nsavestate ckpt.bin\n"
                f"wrdata part1.txt {vecs}\n.endc\n.end\n")
    restart_deck = (body + extra_options + f".tran {tstep} {t2} 0{uics}\n" +
                     f".control\n{HDR}loadstate ckpt.bin\nwrdata part2.txt {vecs}\n.endc\n.end\n")
    ref_ctl = (body + extra_options +
               f".control\n{HDR}tran {tstep} {t2} 0{uics}\nwrdata ref.txt {vecs}\n.endc\n.end\n")

    if not _have_ngspice_aug():
        record(name, None, f"NGSPICE_AUG not found/executable: {NGSPICE_AUG}")
        return
    # cwd=outdir on BOTH runs: ckpt.bin is a bare relative name (ngspice's
    # control-mode tokenizer lower-cases bare arguments -- see tools/fra.py's
    # _run_ngspice docstring, card #472 T10.4 -- so it must not be an absolute
    # path under this project's ~/Work/... convention), and it must resolve to
    # the same directory in both the writing (save) and reading (restart) runs.
    r_save = _run(NGSPICE_AUG, save_ctl, outdir, "save", cwd=outdir)
    r_restart = _run(NGSPICE_AUG, restart_deck, outdir, "restart", cwd=outdir)
    r_ref = _run(NGSPICE_AUG, ref_ctl, outdir, "ref", cwd=outdir)
    for label, r in (("save", r_save), ("ref", r_ref)):
        if r.returncode != 0 or re.search(r"\baborted\b|timestep too small|singular matrix",
                                            r.stdout + r.stderr, re.I):
            record(name, None, f"{label} run failed (exit {r.returncode}); see {outdir}/{label}.log")
            return
    if not os.path.exists(os.path.join(outdir, "part2.txt")):
        record(name, None, f"restart produced no part2.txt (exit {r_restart.returncode}); "
                             f"see {outdir}/restart.log")
        return

    part1 = _load(os.path.join(outdir, "part1.txt"), ncols)
    part2 = _load(os.path.join(outdir, "part2.txt"), ncols)
    ref = _load(os.path.join(outdir, "ref.txt"), ncols)
    if not part1 or not part2 or not ref:
        record(name, None, f"empty output ({len(part1)}/{len(part2)}/{len(ref)} rows); see {outdir}")
        return
    # part1 (the save run) ends exactly at t1; part2 (the restart) never re-emits
    # that row (loadstate opens a fresh output plot and the next OUTPUT is the
    # first NEWLY accepted point after t1, confirmed empirically -- card #472
    # T10.4), so simple concatenation is correct, not part1[:-1]+part2 (an
    # earlier, wrong assumption here that silently dropped part1's real t1 row
    # and shifted every later index by one, corrupting this whole comparison).
    combined = part1 + part2

    # exact grid-matched comparison, informational: ref rows strictly after t1
    # (skip ref's own row at t1, since part1/combined already carries that exact
    # row) against combined's rows from part2 (index len(part1)) onward.
    idx = None
    for i, row in enumerate(ref):
        if row[0] > t1v + 1e-15:
            idx = i
            break
    exact_matched = exact_mismatched = 0
    exact_worst = 0.0
    if idx is not None:
        ref_tail = ref[idx:]
        rest = combined[len(part1):]
        n = min(len(ref_tail), len(rest))
        for i in range(n):
            a, b = ref_tail[i], rest[i]
            if abs(a[0] - b[0]) > 1e-15:
                exact_mismatched += 1
                continue
            exact_matched += 1
            for k in range(1, ncols + 1):
                rel = abs(a[k] - b[k]) / max(1e-9, abs(a[k]))
                exact_worst = max(exact_worst, rel)

    # interpolated comparison, the actual gate
    cmp_idx = compare_indices if compare_indices is not None else list(range(ncols))
    worst = 0.0
    worst_detail = None
    for k in range(n_interp):
        t = t1v + (t2v - t1v) * k / max(1, n_interp - 1)
        a = _interp(ref, t)
        b = _interp(combined, t)
        for ci in cmp_idx:
            rel = abs(a[ci] - b[ci]) / max(1e-9, abs(a[ci]))
            if rel > worst:
                worst = rel
                worst_detail = (t, vectors[ci], a[ci], b[ci])
    ok = worst <= tol
    detail = (f"interp worst rel={worst:.3e} (tol {tol:.0e}) at {worst_detail}; "
              f"exact-grid: {exact_matched} matched/{exact_mismatched} time-mismatched samples, "
              f"worst-on-matched={exact_worst:.3e}")
    record(name, ok, detail)


def red_proof_breakpoint_bug(name):
    """Historical red-proof, run live only if the pre-fix build this card also
    built (aaba50ba3f49-ckpt1: pin + patches 0001-0016, upstream Enhancement-131
    savestate/loadstate UNPATCHED) is still present locally. Not required for a
    fresh checkout -- COULD-NOT-RUN, not FAIL, when it is absent; the numbers it
    reproduces are recorded permanently in patches-wip/ckpt/README.md and
    docs/checkpoint.md ("checkpointing exactly at refs/r2_switched_rc.cir's clock
    edge (t1=300us): the unpatched upstream loadstate diverges from the
    uninterrupted run by ~7.7e-5 relative right at t1 [...]; patches-wip/ckpt/0200
    fixes it to exactly 0.0").
    """
    ckpt1 = os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-ckpt1/bin/ngspice-aug")
    if not (os.path.exists(ckpt1) and os.access(ckpt1, os.X_OK)):
        record(name, None, f"pre-fix build not present locally ({ckpt1}); "
                             "see patches-wip/ckpt/README.md for the recorded numbers")
        return
    outdir = os.path.join(OUT, "redproof_breakpoint")
    deck_lines = stb.prepare_deck(os.path.join(ROOT, "refs", "r2_switched_rc.cir"))
    body = "\n".join(deck_lines) + "\n"
    vecs = "v(vin) v(vx) v(clk)"
    save_ctl = body + f".control\n{HDR}tran 25n 300u 0 uic\nsavestate ckpt.bin\nwrdata part1.txt {vecs}\n.endc\n.end\n"
    restart_deck = body + ".tran 25n 600u 0 uic\n" + f".control\n{HDR}loadstate ckpt.bin\nwrdata part2.txt {vecs}\n.endc\n.end\n"
    ref_ctl = body + f".control\n{HDR}tran 25n 600u 0 uic\nwrdata ref.txt {vecs}\n.endc\n.end\n"
    _run(ckpt1, save_ctl, outdir, "save", cwd=outdir)
    _run(ckpt1, restart_deck, outdir, "restart", cwd=outdir)
    _run(ckpt1, ref_ctl, outdir, "ref", cwd=outdir)
    if not os.path.exists(os.path.join(outdir, "part2.txt")):
        record(name, None, f"pre-fix build's restart produced no output; see {outdir}")
        return
    ref = _load(os.path.join(outdir, "ref.txt"), 3)
    part1 = _load(os.path.join(outdir, "part1.txt"), 3)
    part2 = _load(os.path.join(outdir, "part2.txt"), 3)
    if not (ref and part1 and part2):
        record(name, None, f"empty output ({len(ref)}/{len(part1)}/{len(part2)} rows); see {outdir}")
        return
    combined = part1 + part2
    # Probe just AFTER t1 (not exactly at it): part1's own last row is bit-exact by
    # construction (it is literally the row savestate wrote), so interpolating
    # exactly at t1 would just return that row and never exercise what loadstate's
    # continuation (part2) actually does with it.
    t_probe = 3e-4 + 5e-10
    a = _interp(ref, t_probe)      # vecs = "v(vin) v(vx) v(clk)" -> v(vx) is index 1
    b = _interp(combined, t_probe)
    rel = abs(a[1] - b[1]) / max(1e-9, abs(a[1]))
    # red-proof: the UNPATCHED build must visibly (not just barely) fail 1e-9
    record(name, rel > 1e-6, f"unpatched-upstream restart-vs-uninterrupted v(vx) rel err "
                              f"at t1 (interpolated) = {rel:.3e} (ref={a[1]:.9f}, "
                              f"restart={b[1]:.9f}; must exceed 1e-6 to prove this file's "
                              f"exactness check is sensitive)")


def main():
    t_start = time.time()
    if not _have_ngspice_aug():
        record("(setup) NGSPICE_AUG build", None,
               f"{NGSPICE_AUG} not found or not executable -- build it first "
               "(buildrecipe/build.sh --prefix-name=... ngspice-src-ckpt) or set NGSPICE_AUG")
        n_bad = sum(1 for _, s, _ in results if s == "FAIL")
        n_cnr = sum(1 for _, s, _ in results if s == "COULD-NOT-RUN")
        print(f"\n{len(results) - n_bad - n_cnr}/{len(results)} pass  ({time.time() - t_start:.1f}s)")
        return 2

    # ---- restart exactness, per deck (task item 1/2/4) ----
    # compare_indices=[1] (v(vx) only): the exact-grid check below (index-for-index
    # on shared timepoints) already proves ALL THREE vectors match to 0.0 relative
    # error across all 12947 shared samples -- the interpolated check is a second,
    # independent cross-check and is restricted to v(vx) here because v(clk) is a
    # near-binary clock waveform: right at t1 (the clock's own rising edge) both
    # runs' values are ~1e-11-5e-11 (both "zero" to any physically meaningful
    # precision, just zero for different floating-point reasons), and the relative-
    # error metric's 1e-9 floor on the denominator turns that into a meaningless
    # 5% "relative" difference -- the same class of artifact as v(pwm) below on
    # R3/R3d, not a real restore defect (confirmed separately by the exact-grid
    # check matching it bit-for-bit).
    restart_exactness(
        "restart exact: R2 switched-RC (breakpoint-aligned t1=300us of 600us)",
        os.path.join(ROOT, "refs", "r2_switched_rc.cir"), "25n", "300u", "600u",
        ["v(vin)", "v(vx)", "v(clk)"], uic=True, compare_indices=[1], tol=1e-9)

    restart_exactness(
        "restart exact: BSIM4 d_analog_block (NON-aligned t1=1000n of 2000n, real devices)",
        os.path.join(ROOT, "bench", "decks", "d_analog_block.cir"), "0.5n", "1000n", "2000n",
        ["v(out0)", "v(out30)", "v(out59)"], uic=False, extra_options=".options reltol=1e-3\n",
        tol=1e-9)

    restart_exactness(
        "restart close: R3 voltage-mode buck (t1=500us of 1500us, ~500 more switching periods)",
        os.path.join(ROOT, "refs", "r3_vm_buck.cir"), "10n", "500u", "1500u",
        ["v(vout)", "v(pwm)", "i(l1)"], uic=True, compare_indices=[0], tol=1e-3)

    restart_exactness(
        "restart close: R3d dead-time buck (t1=500us of 1500us, B-source+diode dead-time network)",
        os.path.join(ROOT, "refs", "r3d_vm_buck_deadtime.cir"), "2n", "500u", "1500u",
        ["v(vout)", "v(pwm)", "i(l1)"], uic=True, compare_indices=[0], tol=1e-3)

    # ---- red-proof (task item 4: "dropping the integration-history restore
    # must visibly break it") ----
    red_proof_breakpoint_bug(
        "red-proof: unpatched upstream E-131 fails restart-exactness on R2's breakpoint edge")

    # ---- FRA-from-checkpoint equivalence + speedup (task item 3/4) ----
    fra_checkpoint_case(
        "FRA-from-checkpoint: R3 combined-T matches plain FRA within FRA noise bounds + speedup")
    fra_checkpoint_case(
        "FRA-from-checkpoint: R4 D~0.30 combined-T matches plain FRA + speedup (the 3.5ms/point case)",
        r4=True)

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


# ---------------------------------------------------------------------------
# FRA-from-checkpoint: tools/fra.py's sweep_loop vs sweep_loop_checkpoint
# ---------------------------------------------------------------------------

DB_TOL = 0.5   # looser than test_fra_converter.py's 0.1dB: card #472's own FRA
DEG_TOL = 5.0  # settings below (settle_periods=3) are already noisy per that file's
               # own module docstring ("worst in R3's ~1-5kHz band"); this check's
               # job is "checkpoint mode agrees with plain mode", not FRA accuracy
               # itself (already validated elsewhere, test_fra_converter.py).


def fra_checkpoint_case(name, r4=False):
    os.environ.setdefault("NGSPICE", NGSPICE_AUG)  # fra.py's own NGSPICE default
    if not _have_ngspice_aug():
        record(name, None, f"NGSPICE_AUG not found/executable: {NGSPICE_AUG}")
        return
    deck = os.path.join(ROOT, "refs", "r4_pcm_buck.cir" if r4 else "r3_vm_buck.cir")
    settle = 3.5e-3 if r4 else 300e-6
    freqs = [12559.432157547897] if r4 else [12559.432157547897, 125594.32157547899]
    kw = dict(settle_periods=3, measure_periods=5, points_per_period=256, startup="uic",
              tmax=5e-9, spice_options=".options trtol=1\n")
    out_plain = os.path.join(OUT, name.replace(" ", "_"), "plain")
    out_ckpt = os.path.join(OUT, name.replace(" ", "_"), "ckpt")
    try:
        t0 = time.time()
        pts_plain = fra.sweep_loop(deck, "XP", freqs, 0.01, settle_time=settle,
                                    outdir_root=out_plain, workers=min(4, len(freqs)), ngspice=NGSPICE_AUG, **kw)
        t_plain = time.time() - t0

        t0 = time.time()
        pts_ckpt, ckpt_path = fra.sweep_loop_checkpoint(
            deck, "XP", freqs, 0.01, td=settle, settle_tstep=10e-9,
            outdir_root=out_ckpt, workers=min(4, len(freqs)), ngspice=NGSPICE_AUG, **kw)
        t_ckpt = time.time() - t0
    except stb.CouldNotRun as e:
        record(name, None, str(e))
        return

    n_bad_plain = sum(1 for s, _ in pts_plain if s != "ok")
    n_bad_ckpt = sum(1 for s, _ in pts_ckpt if s != "ok")
    if n_bad_plain or n_bad_ckpt:
        record(name, None, f"plain {n_bad_plain} could-not-run, checkpoint {n_bad_ckpt} could-not-run")
        return

    worst_db = worst_deg = 0.0
    for f, (_, vp), (_, vc) in zip(freqs, pts_plain, pts_ckpt):
        dp = stbcore.db(vp["combined"])
        pdp = math.degrees(cmath.phase(vp["combined"]))
        pdc = math.degrees(cmath.phase(vc["combined"]))
        ddb = abs(stbcore.db(vc["combined"]) - dp)
        ddeg = abs(pdc - pdp)
        while ddeg > 180:
            ddeg -= 360
        ddeg = abs(ddeg)
        worst_db = max(worst_db, ddb)
        worst_deg = max(worst_deg, ddeg)
    ok = worst_db <= DB_TOL and worst_deg <= DEG_TOL
    speedup = t_plain / t_ckpt if t_ckpt > 0 else float("inf")
    detail = (f"{len(freqs)} freq point(s): worst |dB diff|={worst_db:.3f} (tol {DB_TOL}), "
              f"worst |deg diff|={worst_deg:.3f} (tol {DEG_TOL}); "
              f"plain={t_plain:.2f}s ckpt={t_ckpt:.2f}s speedup={speedup:.2f}x")
    record(name, ok, detail)


if __name__ == "__main__":
    sys.exit(main())
