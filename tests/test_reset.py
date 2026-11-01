#!/usr/bin/env python3
"""tests/test_reset.py -- card #492: state leak across `reset`.

Card #491 (docs/campaign.md, mode (b) section) found that
`campaign_pvt_r3d`'s single-process `alterparam`/`option temp`/`reset`/
`run` loop disagreed with fresh-process results by 2-4% relative on
`duty_final`/`vout_ripple`, concentrated at the `temp=85` corners. Root
cause (research/reset-leak.md, full trace): `reset` (`com_rset()` ->
`com_remcirc()`, `src/frontend/runcoms2.c`) unconditionally frees
`ft_curckt->ci_vars` -- the list an interactive `option name=val` (e.g.
`option temp = 85`) is cached onto so the NEXT circuit (re)load can replay
it onto the fresh circuit (`inp_dodeck()`'s own replay loop, `inp.c`) --
and since `inp_dodeck()`'s `reuse` parameter is, at its one call site,
always `FALSE`, the freshly rebuilt circuit after `reset` always starts
with an EMPTY `ci_vars`: any `option ...=...` set since the circuit was
first loaded is silently discarded, and the new circuit's task falls back
to `CKTnewTask()`'s hard-coded default (27C) instead of the value asked
for -- no error, no warning. Reproduced identically on stock
`/usr/local/bin/ngspice` 45.2 and the unpatched pin (`aaba50ba3f49`):
upstream, not introduced by any patch in this series. `alterparam` is NOT
affected (it edits the deck text itself, which IS correctly replayed by a
reload) -- this is why the bug was masked until a per-point TEMPERATURE
sweep (which goes in only through `option`) specifically exercised it.

Fix (`patches-wip/reset/0800-...patch`, one commit): `com_rset()` detaches
`ci_vars` before `com_remcirc()` and replays each variable through
`cp_vset()` (the same call the interactive `option` command itself makes)
onto the freshly reloaded circuit.

This file checks two things, at two levels:

  1. `check_temp_probe_fixed` / `check_temp_probe_red_proof` -- the direct,
     cheap, deck-independent reproducer (`tests/fixtures/reset_temp_probe.cir`):
     read back a device's OWN operating temperature (`@d1[temp]`) after each
     step of an `option temp=X` / `reset` sequence. Unambiguous: this is not
     inferring the bug from a downstream measure's sensitivity, it is
     reading the exact field the bug corrupts.
  2. `check_r3d_loop_vs_fresh` / `check_r3d_loop_red_proof` -- the
     campaign-level check card #492 asks for: `campaign_pvt_r3d`'s mode (b)
     loop vs fresh-process mode (a), on a sample DELIBERATELY including the
     worst corners found in research/reset-leak.md (temp=85 and temp=-20,
     not just temp=25 -- card #491's own post-mortem is explicit that a
     small RANDOM/leading sample can pass by luck while the full grid does
     not; this sample is chosen to not repeat that mistake), must be
     IDENTICAL (0 relative -- the fix makes this exact, not approximately
     close; see research/reset-leak.md's own 6-point bit-identical result)
     after the fix.

Both levels also carry a red-proof: run against an UNFIXED build (default
aaba50ba3f49-s39, pin + patches/0001-0039, all of which reproduce the bug
per research/reset-leak.md) and assert the check actually FAILS there --
proving the check detects the regression, not just that the fixed build
happens to pass.

Usage:
  NGSPICE_RESET_FIXED=~/Work/apps/ngspice_augmentation/aaba50ba3f49-reset1/bin/ngspice-aug \\
  NGSPICE_RESET_UNFIXED=~/Work/apps/ngspice_augmentation/aaba50ba3f49-s39/bin/ngspice-aug \\
      python3 tests/test_reset.py

Exit code: 0 if every check PASSed, 1 if any FAILed, 2 if any COULD-NOT-RUN
and none FAILed.
"""
import os
import pathlib
import re
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "bench"))

import campaign  # noqa: E402

NGSPICE_FIXED = os.environ.get("NGSPICE_RESET_FIXED") or os.environ.get("NGSPICE") or os.path.expanduser(  # 2026-10-02: the build under test is the "fixed" one
    "~/Work/apps/ngspice_augmentation/aaba50ba3f49-reset1/bin/ngspice-aug")
NGSPICE_UNFIXED = os.environ.get("NGSPICE_RESET_UNFIXED") or os.path.expanduser(
    "~/Work/apps/ngspice_augmentation/aaba50ba3f49-s39/bin/ngspice-aug")
DATA_ROOT = pathlib.Path(os.path.expanduser("~/Work/data/ngspice_augmentation/campaign/tests_reset"))

PROBE_FIXTURE = HERE / "fixtures" / "reset_temp_probe.cir"

# The five '@d1[temp] = ...' prints the fixture produces, in order, and what
# each step SHOULD read back on a correct ngspice. Step indices match the
# fixture's own comments.
PROBE_EXPECTED_FIXED = [27.0, 85.0, 85.0, -20.0, -20.0]
# What the KNOWN-BROKEN mechanism actually produces (stuck at the 27C
# default from the first reset onward) -- used only by the red-proof, and
# only to confirm the check exercises the bug, not to assert this is
# "correct" behaviour anywhere.
PROBE_EXPECTED_UNFIXED = [27.0, 85.0, 27.0, 27.0, 27.0]

_TEMP_LINE_RE = re.compile(r"@d1\[temp\]\s*=\s*([-\d.eE+]+)")


def _ok(path):
    return os.path.exists(path) and os.access(path, os.X_OK)


def _run_probe(ngspice):
    out = subprocess.run([ngspice, "-b", str(PROBE_FIXTURE)],
                          capture_output=True, text=True, timeout=60)
    vals = [float(m.group(1)) for m in _TEMP_LINE_RE.finditer(out.stdout)]
    return vals, out


def check_temp_probe_fixed():
    name = "temp_probe[fixed]"
    if not _ok(NGSPICE_FIXED):
        return "COULD-NOT-RUN", name, f"fixed ngspice not found/executable: {NGSPICE_FIXED}"
    vals, out = _run_probe(NGSPICE_FIXED)
    if len(vals) != 5:
        return "FAIL", name, f"expected 5 @d1[temp] prints, got {len(vals)}: {vals} " \
            f"(stdout tail: {out.stdout[-500:]})"
    if vals != PROBE_EXPECTED_FIXED:
        return "FAIL", name, f"expected {PROBE_EXPECTED_FIXED}, got {vals} -- " \
            "temperature not correctly propagated across reset"
    return "pass", name, f"{vals} matches expected (option temp=X survives reset)"


def check_temp_probe_red_proof():
    name = "temp_probe_red_proof[unfixed]"
    if not _ok(NGSPICE_UNFIXED):
        return "COULD-NOT-RUN", name, f"unfixed ngspice not found/executable: {NGSPICE_UNFIXED}"
    vals, out = _run_probe(NGSPICE_UNFIXED)
    if len(vals) != 5:
        return "COULD-NOT-RUN", name, f"expected 5 @d1[temp] prints, got {len(vals)}: {vals}"
    if vals == PROBE_EXPECTED_FIXED:
        return "FAIL", name, f"red-proof did not reproduce the bug: unfixed build printed " \
            f"{vals} (the CORRECT sequence) -- either the bug is fixed upstream on this " \
            "binary too, or NGSPICE_RESET_UNFIXED points at an already-fixed build"
    if vals != PROBE_EXPECTED_UNFIXED:
        return "FAIL", name, f"unfixed build produced neither the expected-correct nor the " \
            f"expected-broken sequence: {vals} (expected broken: {PROBE_EXPECTED_UNFIXED})"
    return "pass", name, f"{vals} reproduces the known bug (stuck at 27C after the first reset)"


# --------------------------------------------------------------------- #
# Campaign-level check: campaign_pvt_r3d, mode (b) loop vs mode (a) fresh
# --------------------------------------------------------------------- #

# Deliberately not a leading/random slice of POINTS (card #491's own
# post-mortem, docs/campaign.md: "a small-sample identity check is not
# sufficient... the divergence here is corner-dependent, and the sample
# that passed did not happen to include the worst corners"). campaign_pvt_r3d
# POINTS is itertools.product(temps=[-20,25,85], vins=[10,12,14],
# rloads=[1.5,3.3,6.6,10.0]) in that nesting order, i.e. index = temp_i*12 +
# vin_i*4 + rload_i. Indices below: one full (vin, rload) row at temp=-20
# (idx 0-3) and the full temp=85 block (idx 24-35, all 12 combinations --
# the corner research/reset-leak.md found the effect concentrated at) plus
# two temp=25 points (idx 12, 20) as an unaffected control.
R3D_SAMPLE_IDX = list(range(0, 4)) + [12, 20] + list(range(24, 36))

# The fix makes mode (b) match mode (a) EXACTLY on this sample (measured:
# research/reset-leak.md's 6-point subset was bit-identical, 0.000e+00
# relative, all 4 measures) -- not "close within FP noise". A nonzero
# tolerance here would hide a regression of the very bug this file exists
# to catch, so rtol/atol are 0 for every measure, not the looser
# 1e-3/1e-2 tests/test_campaign.py uses for its OWN (pre-fix,
# known-imperfect) identity check.
R3D_TOL_EXACT = {m: (0.0, 0.0) for m in
                 ["t_period", "duty_final", "vout_avg_final", "vout_ripple"]}


def _run_r3d_sample(ngspice, label, timeout_a=60):
    mod = campaign.load_campaign_module(str(ROOT / "bench" / "campaign_pvt_r3d.py"))
    sample_points = [mod.POINTS[i] for i in R3D_SAMPLE_IDX]
    outdir = DATA_ROOT / f"r3d_{label}_{int(time.time())}"
    manifest_a = campaign.run_campaign_a(
        sample_points, mod.build_point_deck, ngspice, str(outdir / "a"),
        timeout=timeout_a, max_procs=6, probe_n=2)
    manifest_b = campaign.run_campaign_b(
        sample_points, mod.build_loop_prelude, mod.build_loop_entry, mod.build_loop_tail,
        ngspice, str(outdir / "b"), timeout=timeout_a * len(sample_points))
    return manifest_a, manifest_b, len(sample_points)


def check_r3d_loop_vs_fresh():
    name = "r3d_loop_vs_fresh[fixed]"
    if not _ok(NGSPICE_FIXED):
        return "COULD-NOT-RUN", name, f"fixed ngspice not found/executable: {NGSPICE_FIXED}"
    try:
        manifest_a, manifest_b, n = _run_r3d_sample(NGSPICE_FIXED, "fixed")
    except Exception as e:
        return "COULD-NOT-RUN", name, f"exception running campaign: {e}"
    if manifest_a["n_ok"] != n or manifest_b["n_ok"] != n:
        return "COULD-NOT-RUN", name, f"not every point completed ok: mode a " \
            f"{manifest_a['n_ok']}/{n}, mode b {manifest_b['n_ok']}/{n}"
    ok, rows = campaign.compare_results(manifest_a["results"], manifest_b["results"],
                                         list(R3D_TOL_EXACT), rtol=0.0, atol=0.0)
    if not ok:
        bad = [r for r in rows if not r["match"]]
        mod = campaign.load_campaign_module(str(ROOT / "bench" / "campaign_pvt_r3d.py"))
        detail = "; ".join(
            f"sample_pos{r['idx']}(point_idx={R3D_SAMPLE_IDX[r['idx']]}, "
            f"temp={mod.POINTS[R3D_SAMPLE_IDX[r['idx']]]['temp']}) "
            f"{r['measure']}: a={r['a']} b={r['b']} {r['detail']}" for r in bad)
        return "FAIL", name, detail
    return "pass", name, f"{n} points ({R3D_SAMPLE_IDX}) x {len(R3D_TOL_EXACT)} measures, " \
        "bit-identical mode (a) vs mode (b)"


def check_r3d_loop_red_proof():
    name = "r3d_loop_red_proof[unfixed]"
    if not _ok(NGSPICE_UNFIXED):
        return "COULD-NOT-RUN", name, f"unfixed ngspice not found/executable: {NGSPICE_UNFIXED}"
    try:
        manifest_a, manifest_b, n = _run_r3d_sample(NGSPICE_UNFIXED, "unfixed")
    except Exception as e:
        return "COULD-NOT-RUN", name, f"exception running campaign: {e}"
    if manifest_a["n_ok"] != n or manifest_b["n_ok"] != n:
        return "COULD-NOT-RUN", name, f"not every point completed ok: mode a " \
            f"{manifest_a['n_ok']}/{n}, mode b {manifest_b['n_ok']}/{n}"
    # loose-ish tolerance matching tests/test_campaign.py's own PVT_TOL for
    # the smooth measures, tight (0) for duty_final/vout_ripple -- the two
    # measures research/reset-leak.md found the bug in. The red-proof must
    # fail on THOSE, specifically, not on generic FP noise elsewhere.
    ok, rows = campaign.compare_results(manifest_a["results"], manifest_b["results"],
                                         ["duty_final", "vout_ripple"], rtol=0.0, atol=0.0)
    if ok:
        return "FAIL", name, "red-proof did not reproduce the bug: unfixed build's mode (b) " \
            "loop matched mode (a) exactly on duty_final/vout_ripple -- either the bug is " \
            "fixed on this binary too, or NGSPICE_RESET_UNFIXED points at an already-fixed build"
    n_mismatch = sum(1 for r in rows if not r["match"])
    return "pass", name, f"{n_mismatch}/{len(rows)} duty_final/vout_ripple mismatches " \
        "reproduce the known bug on the unfixed build"


CHECKS = [
    check_temp_probe_fixed,
    check_temp_probe_red_proof,
    check_r3d_loop_vs_fresh,
    check_r3d_loop_red_proof,
]


def main():
    results = []
    for fn in CHECKS:
        state, name, detail = fn()
        results.append((state, name, detail))
        print(f"[{state:14s}] {name}: {detail}")

    n_fail = sum(1 for s, _, _ in results if s == "FAIL")
    n_cnr = sum(1 for s, _, _ in results if s == "COULD-NOT-RUN")
    n_pass = sum(1 for s, _, _ in results if s == "pass")
    print(f"\n{n_pass} pass, {n_fail} FAIL, {n_cnr} COULD-NOT-RUN (of {len(results)})")
    if n_fail:
        print("OVERALL FAIL")
        return 1
    if n_cnr:
        print("OVERALL COULD-NOT-RUN")
        return 2
    print("OVERALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
