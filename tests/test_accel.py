#!/usr/bin/env python3
"""tests/test_accel.py -- accel_same_results checks for the transient-
ACCELERATION work (ngspice_augmentation card #467, plan T9.1).

For every accepted speedup, this runs each bench deck ONCE under the
baseline prefix and ONCE under the candidate prefix/config, and compares
the deck's own named `.meas`/`print` values within the tolerance stated
for that candidate (chosen BEFORE the candidate was measured -- see
CANDIDATES below and bench/README.md / the card report for where each
number came from). This is a correctness *check* (pass/FAIL/COULD-NOT-RUN),
not the wall-clock benchmark -- that is bench/runner.py's job, run
separately with N>=5 reps under memory caps.

Usage:
  python3 tests/test_accel.py                 # default: ACCEPTED candidates only
                                                # (+ the raw-binary-output regression
                                                # guard), so a plain run never needs
                                                # --candidate-bin and never reports
                                                # COULD-NOT-RUN for a rejected
                                                # candidate that isn't being re-tried.
  python3 tests/test_accel.py --all            # also re-run every REJECTED/pending
                                                # candidate (needs --candidate-bin for
                                                # the ones that were only ever built as
                                                # a separate prefix, e.g.
                                                # omp-schedule-runtime, march-native).
  python3 tests/test_accel.py --candidate klu-everywhere   # one candidate by name,
                                                # regardless of its status.
  python3 tests/test_accel.py --deck d_analog_block --candidate omp-schedule-runtime \
      --candidate-bin ~/Work/apps/ngspice_augmentation/aaba50ba3f49-accel-ompsched/bin/ngspice-aug

Card #478 T9.2 task 0: earlier, the default (no-argument) run walked every
candidate in CANDIDATES regardless of status, so it always reported OVERALL
FAIL -- not because anything was wrong, but because the rejected
omp-schedule-runtime/march-native candidates need --candidate-bin (a
rebuilt prefix) that a default run never supplies. Each candidate now
states its own "status" ("accepted"/"rejected"); the default run only
exercises "accepted" ones (only raw-binary-output today), so OVERALL FAIL
now means an actual regression, not an unmet prerequisite. --all restores
the old "walk everything" behaviour for an explicit full audit.

Exit code: 0 if every check PASSed, 1 if any FAILed or COULD-NOT-RUN.
"""
import argparse
import os
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
DECKS_DIR = REPO / "bench" / "decks"

# card #467's original correctness-check baseline. Kept as s12 (not the
# round-3 A/B wall-clock baseline, s15) because the two prefixes differ
# ONLY in the ngspice-aug wrapper's OMP_WAIT_POLICY default (s15 adds
# passive, a pure runtime-env change with no effect on simulation
# results) -- same code, same .meas/.print output, so s12 remains a valid
# (and already-built) reference for this correctness harness. Wall-clock
# A/B comparisons use s15 (bench/README.md, patches-wip/accel/README.md).
BASELINE_BIN = os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s12/bin/ngspice-aug")

# deck -> (relative .cir path from DECKS_DIR, cwd for the run, memcap)
DECKS = {
    "a_pex_buck": ("a_pex_buck.cir", DECKS_DIR, "6G"),
    "b_long_switch": ("b_long_switch.cir", DECKS_DIR, "2G"),
    "c_pvt_loop": ("c_pvt_loop.cir", DECKS_DIR, "2G"),
    "d_analog_block": ("d_analog_block.cir", DECKS_DIR, "2G"),
    # card #478 T9.2 additions (bench/README.md): same shapes as above, see
    # bench/gen_analog_block.py --tstop-ns / bench/gen_pex_buck.py --topology tree.
    "d_analog_block_long": ("d_analog_block_long.cir", DECKS_DIR, "2G"),
    "e_pex_tree": ("e_pex_tree.cir", DECKS_DIR, "6G"),
    # card #483 T9.3 task 2 additions: same 5T-OTA shape as d_analog_block,
    # BSIM3 (level=49) instead of BSIM4 -- exercises b3ld.c's OMP load loop
    # (candidate 0105) instead of b4ld.c's. See bench/gen_analog_block.py
    # --device bsim3.
    "f_analog_block_bsim3": ("f_analog_block_bsim3.cir", DECKS_DIR, "2G"),
    "f_analog_block_bsim3_long": ("f_analog_block_bsim3_long.cir", DECKS_DIR, "2G"),
    # card #486 T9.4 (round 4) addition: one B-source per PT_*/PTF_* node
    # type/operator (bench/gen_ad_stress.py), argument-swept via PWL
    # sources so every tree is a real PT_VAR-rooted gradient, not something
    # that folds to a constant at parse time. See that generator's header.
    "g_ad_stress": ("g_ad_stress.cir", DECKS_DIR, "2G"),
}

MEAS_LINE_RE = re.compile(r"^(\w[\w#]*)\s*=\s*([-+0-9.eE]+)\s*$")


def parse_meas(text):
    vals = {}
    for line in text.splitlines():
        m = MEAS_LINE_RE.match(line.strip())
        if m:
            try:
                vals[m.group(1)] = float(m.group(2))
            except ValueError:
                pass
    return vals


def run_deck(ngspice_bin, deck, env=None, mem_cap=None):
    cir_rel, cwd, memcap = DECKS[deck]
    cap = mem_cap or memcap
    cmd = ["systemd-run", "--user", "--scope", "--quiet", "-p", f"MemoryMax={cap}",
           "--", ngspice_bin, "-b", cir_rel]
    full_env = dict(os.environ)
    if env:
        full_env.update(env)
    try:
        proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True,
                               timeout=600, env=full_env)
    except subprocess.TimeoutExpired:
        return None, "COULD-NOT-RUN (timeout)"
    if proc.returncode != 0:
        return None, f"COULD-NOT-RUN (rc={proc.returncode}): {proc.stderr[-500:]}"
    meas = parse_meas(proc.stdout)
    if not meas:
        return None, "COULD-NOT-RUN (no .meas/print output parsed)"
    return meas, None


def compare(baseline_meas, cand_meas, tol_rel, tol_abs=0.0):
    problems = []
    for k, bv in baseline_meas.items():
        if k not in cand_meas:
            problems.append(f"{k}: missing in candidate output")
            continue
        cv = cand_meas[k]
        denom = max(abs(bv), abs(cv), 1e-30)
        rel = abs(bv - cv) / denom
        if rel > tol_rel and abs(bv - cv) > tol_abs:
            problems.append(f"{k}: baseline={bv!r} candidate={cv!r} rel_diff={rel:.3e} > tol={tol_rel:.3e}")
    return problems


# --- candidate registry ---------------------------------------------------
# Each candidate: which decks it applies to, how to build the candidate's
# ngspice invocation (binary path and/or env vars), and the tolerance
# chosen BEFORE it was measured (spec.md accel_same_results: default 1e-6
# relative for linear-solver/threading changes; model-bypass/step-control
# changes state and justify their own tolerance).
CANDIDATES = {
    "klu-everywhere": {
        "decks": ["b_long_switch", "c_pvt_loop", "d_analog_block"],
        "bin": BASELINE_BIN,
        "env": {},
        "spiceinit": "option klu\n",
        "tol_rel": 1e-6,
        "status": "rejected",  # card #467: no deck reached the +10% accel_wallclock
                                # bar, d_analog_block regressed 10.7% -- see
                                # bench/README.md's candidates table.
        "note": "linear-solver change (SPARSE->KLU on decks that didn't already select it); "
                "default 1e-6 relative per spec.md.",
    },
    "bypass-1": {
        "decks": ["a_pex_buck", "b_long_switch", "c_pvt_loop", "d_analog_block"],
        "bin": BASELINE_BIN,
        "env": {},
        "spiceinit": "option bypass=1\n",
        "tol_rel": 1e-3,
        "status": "rejected",  # card #467: d_analog_block regressed 33.9%.
        "note": "model-bypass change: BSIM4's own convergence-iteration bypass "
                "(src/spicelib/devices/bsim4/b4ld.c, #ifndef NOBYPASS), gated at "
                "runtime by .options bypass. Tolerance set to 1e-3 relative on named "
                ".meas values BEFORE measuring, justified against each bench deck's "
                "own ripple/settle magnitude (all bench decks' ripple/settle windows "
                "are >>1e-3 relative of the signal, per bench/README.md's baseline "
                "table) -- not the spec's linear-solver default, because bypass "
                "changes which Newton iterations get a full re-evaluation, which can "
                "move the last few significant digits of a converged value.",
    },
    "omp-schedule-runtime": {
        "decks": ["a_pex_buck", "d_analog_block"],
        "bin": None,  # filled in by --candidate-bin at run time (the patched prefix)
        "env": {"OMP_SCHEDULE": "static"},
        "spiceinit": None,
        "tol_rel": 1e-9,
        "status": "rejected",  # card #467: correct but slower everywhere it
                                # matters (schedule-dispatch overhead);
                                # patches-wip/accel/README.md.
        "note": "threading-only change: schedule(runtime) instead of the implied "
                "static default on BSIM4's #pragma omp parallel for, no arithmetic "
                "change, no shared-state write in the parallel region (research/"
                "transient-acceleration.md section 3). Expect near-bit-identical; "
                "1e-9 relative (not 0) only to absorb any non-associative FP summation "
                "order difference in the *unrelated* serial KLU/SPARSE solve if thread "
                "count otherwise perturbs cache/alloc timing -- not the OMP region itself.",
    },
    "march-native": {
        "decks": ["a_pex_buck", "b_long_switch", "c_pvt_loop", "d_analog_block"],
        "bin": None,  # filled in by --candidate-bin
        "env": {},
        "spiceinit": None,
        "tol_rel": 1e-9,
        "status": "rejected",  # card #467: FAILS accel_same_results outright
                                # (~2e-4 rel diff on c_pvt_loop), independent
                                # of the also-not-a-win speed.
        "note": "compiler-flag change (-O3 -march=native). FP reassociation from "
                "different vectorisation/FMA codegen can move the last bit or two; "
                "1e-9 relative chosen before measuring, well inside every bench "
                "deck's own reltol=1e-3 solver tolerance.",
    },
    # card #478 T9.2 (round 2) candidates -- patches-wip/accel/0101-0103,
    # patches-wip/accel/README.md has the full numbers/investigation.
    "accel-0101-0102-bsource-cache": {
        "decks": ["a_pex_buck", "b_long_switch", "c_pvt_loop", "d_analog_block", "e_pex_tree"],
        "bin": None,  # filled in by --candidate-bin (a prefix built from
                      # patches-wip/accel/0101 + 0102 on top of aug)
        "env": {},
        "spiceinit": None,
        "tol_rel": 0.0,
        "status": "rejected",  # card #478: bit-identical (0/0.0 tol) and a
                                # real 22-27% win on b_long_switch, BUT a_pex_buck
                                # regresses 9.6%/6.9% (order-swap-confirmed real,
                                # not noise; mechanism not conclusively identified --
                                # see patches-wip/accel/README.md) -- violates
                                # accel_wallclock's "no deck slower by >5%" bar.
                                # Ranked candidate #1 for a future round (narrow the
                                # cache or root-cause the a_pex_buck cost).
        "note": "PTeval per-IFeval()-call memoization (0101) + self-referential "
                "derivative-node reuse (0102) for B-source expressions. Pure "
                "refactor, no arithmetic change -- tolerance 0 (bit-identical), "
                "confirmed on every deck listed.",
    },
    "accel-0103-omp-schedule-guided": {
        "decks": ["d_analog_block_long", "e_pex_tree", "b_long_switch", "c_pvt_loop"],
        "bin": None,
        "env": {},
        "spiceinit": None,
        "tol_rel": 0.0,
        "status": "folded",    # card #483 T9.3: the parent folded patches-wip/accel/
                                # 0103 into patches/0015 between rounds 2 and 3 --
                                # it is now part of `aug` itself (every prefix built
                                # from `aug`, including BASELINE_BIN and the round-3
                                # A/B wall-clock baseline s15, already has this
                                # change). There is no longer a meaningful separate
                                # "candidate" binary to A/B against: comparing
                                # BASELINE_BIN to itself would trivially pass. Kept
                                # here as a historical record, not a runnable check
                                # -- "folded" is deliberately neither "accepted" NOR
                                # "rejected"/"pending" so neither a default run NOR
                                # --all tries to execute it (both would otherwise
                                # need a --candidate-bin that no longer means
                                # anything). Explicit `--candidate
                                # accel-0103-omp-schedule-guided --candidate-bin
                                # <old prefix>` still works if ever needed.
        "note": "hard-coded schedule(guided) instead of the implied schedule(static) "
                "on BSIM4's OMP load loop -- compile-time-specialised (unlike round "
                "1's rejected schedule(runtime)), no arithmetic change, no shared-"
                "state write in the parallel region -- tolerance 0 (bit-identical). "
                "ACCEPTED card #478, now patches/0015.",
    },
    # card #483 T9.3 (round 3) candidates -- patches-wip/accel/0104-0105,
    # patches-wip/accel/README.md has the full numbers/investigation.
    "accel-0104-bsource-cache-sidetable": {
        "decks": ["a_pex_buck", "b_long_switch", "c_pvt_loop", "d_analog_block", "e_pex_tree"],
        "bin": None,  # filled in by --candidate-bin (a prefix built from
                      # patches-wip/accel/0104 -- and 0102, carried forward
                      # unmodified -- on top of aug/patches/0015)
        "env": {},
        "spiceinit": None,
        "tol_rel": 0.0,
        "status": "folded",  # merged as patches/0030-0031 (2026-09-29)
                                # clears accel_wallclock) but status stays "pending"
                                # not "accepted", same reasoning 0103 documented
                                # last round: bin is still None (needs
                                # --candidate-bin) until the parent folds
                                # patches-wip/accel/0104 into patches/series and a
                                # standing prefix exists -- marking this "accepted"
                                # now would make a plain default `test_accel.py`
                                # run COULD-NOT-RUN every time (task 0's bug,
                                # round 2). bit-identical (0.0 tol) on every deck
                                # listed here plus refs/r3_vm_buck, r4_pcm_buck,
                                # r5_hysteretic_buck (not in DECKS/this table, checked
                                # by hand -- see research/transient-acceleration.md
                                # section 8). Wall-clock, interleaved A/B n=5:
                                # a_pex_buck -3.3% (was -9.6%/-6.9% for round 2's
                                # 0101/0102, the regression this candidate fixes),
                                # b_long_switch +21.2%, c_pvt_loop +19.3% (very
                                # clean), d_analog_block +9.0% (noise, sub-second,
                                # no B-sources), e_pex_tree -3.1%, r3_vm_buck +8.5%,
                                # r4_pcm_buck +26.5%, r5_hysteretic_buck +16.8%.
                                # No deck slower by >5%; several >=10% faster.
                                # ACCEPTED -- clears accel_wallclock cleanly, unlike
                                # round 2's 0101/0102.
        "note": "PTeval per-IFeval()-call memoization, cache moved OUTSIDE "
                "INPparseNode onto INPparseTree (a small per-B-source-instance "
                "open-addressing hash table, not a per-node field) so growing it "
                "no longer perturbs a_pex_buck's cache-sensitive KLU refactor. "
                "Combined with 0102 (self-referential tanh/tan/sqrt/exp node reuse, "
                "carried forward unmodified from round 2 -- it never grew any "
                "struct). Pure refactor, no arithmetic change -- tolerance 0.",
    },
    "accel-0105-omp-schedule-guided-other-loaders": {
        "decks": ["f_analog_block_bsim3_long", "f_analog_block_bsim3", "d_analog_block",
                   "b_long_switch", "c_pvt_loop"],
        "bin": None,  # filled in by --candidate-bin (a prefix built from
                      # patches-wip/accel/0105 on top of aug/patches/0015)
        "env": {},
        "spiceinit": None,
        "tol_rel": 0.0,
        "status": "folded",  # merged as patches/0030-0031 (2026-09-29)
                                # until folded" reasoning as 0104 above and 0103
                                # last round -- bin is None, needs --candidate-bin.
                                # bit-identical everywhere listed (BSIM4
                                # decks d_analog_block/b_long_switch/c_pvt_loop
                                # confirmed no-op, as expected -- this patch never
                                # touches b4ld.c). f_analog_block_bsim3_long (300
                                # BSIM3 instances, 60 cells, the size that matched
                                # d_analog_block_long's convention) gave only
                                # +8.7%/+5.6% (combined ~+7.7%, n=10, order-swap
                                # checked) -- a real, clean, repeatable win but
                                # UNDER the accel_wallclock +10% bar. Doubling to
                                # 120 cells / 600 instances (the deck now checked
                                # into bench/decks/f_analog_block_bsim3_long.cir)
                                # gave a clean +12.7%, n=5, clearing the bar --
                                # BSIM3's simpler model apparently needs more
                                # parallel work per call than BSIM4 did (which
                                # showed 10.4% at only 300 instances, round 2) to
                                # amortise the same fork/join overhead. ACCEPTED at
                                # the 600-instance size; f_analog_block_bsim3
                                # (short, 2us) stays as the bit-identical
                                # regression guard, not a performance claim (too
                                # short for a stable signal, like round 2's
                                # original d_analog_block).
        "note": "hard-coded schedule(guided), same change as 0103/patches/0015, "
                "applied to the other 7 OMP device-loader #pragma omp parallel for "
                "sites found by grep (b3ld.c/b3v32ld.c/hsm2ld.c/b4soild.c/"
                "b4v5ld.c/b4v6ld.c/b4v7ld.c -- 3 more BSIM4-version variants than "
                "the card named, not touched by 0103/0015 which only patched "
                "b4ld.c). OSDI's parallel region (osdiload.c) is a #pragma omp "
                "task structure, not parallel-for, and is NOT touched (no "
                "schedule clause applies to task). Same structure confirmed by "
                "reading each file: flat instance array, per-instance scratch "
                "writes only, deferred serial stamping pass -- no shared "
                "mutable state, iteration-order-independent -- tolerance 0.",
    },
    # card #486 T9.4 (round 4) candidate -- patches-wip/accel/0106,
    # patches-wip/accel/README.md has the full design/coverage/evidence.
    "accel-0106-forward-ad": {
        "decks": ["a_pex_buck", "b_long_switch", "c_pvt_loop", "d_analog_block",
                   "d_analog_block_long", "e_pex_tree", "f_analog_block_bsim3",
                   "f_analog_block_bsim3_long", "g_ad_stress"],
        "bin": None,  # filled in by --candidate-bin (a prefix built from
                      # patches-wip/accel/0106 on top of aug/patches/0031,
                      # i.e. aug-accel4 == pin + patches/0001-0037; card
                      # #486's own build is aaba50ba3f49-accel41)
        "env": {},
        "spiceinit": None,
        "tol_rel": 1e-12,      # stated BEFORE measuring, per this round's brief
                                # ("transient results bit-identical or within
                                # 1e-12 relative (state which)"): 1e-12 is the
                                # number stated here, but what was actually
                                # MEASURED is tighter -- every .meas value on
                                # every deck below matched the baseline
                                # (aaba50ba3f49-s37) at EXACTLY 0.0 relative
                                # (bit-identical `diff`, not just within
                                # tolerance) -- see the note and the round-4
                                # report for the full evidence, including a
                                # separate throwaway self-check binary
                                # (aaba50ba3f49-accel41-selfcheck, never
                                # committed) that compares AD's derivs[]
                                # against the OLD symbolic path's derivs[]
                                # DIRECTLY, call-by-call, inside the SAME
                                # IFeval() invocation: 0 of ~167 million
                                # compared IFeval() calls (every B-source
                                # evaluation across every bench deck, every
                                # B-source ref deck, and
                                # bench/decks/g_ad_stress.cir) exceeded 1e-9
                                # relative -- max relative diff seen across
                                # all of them was 0.0.
        "status": "rejected",  # re-measure 2026-09-29: b_long_switch +3.1 %, r4 +1.5 %, a_pex_buck -0.9 % (no deck >= +10 %); see patches-wip/accel/README.md
                                # deck slower by >5%, worst case e_pex_tree
                                # -2.3%; a_pex_buck +11.9% and b_long_switch
                                # +13.1% both clear the +10% bar -- research/
                                # transient-acceleration.md section 9.6/9.7
                                # has the full table incl. a repeat of
                                # refs/r3d and refs/r4 after their first n=5
                                # pass showed load-noise, not a real
                                # regression). Status stays "pending" not
                                # "accepted" because bin is still None (needs
                                # --candidate-bin) until the parent folds
                                # patches-wip/accel/0106 into patches/series
                                # and a standing prefix exists -- same
                                # reasoning 0103/0104/0105 documented in
                                # earlier rounds.
        "note": "forward-mode (dual-number) AD for B-source expression trees "
                "(src/spicelib/parser/ifeval.c PTevalDual()): value + every "
                "partial derivative computed in ONE recursive walk of the "
                "value tree, replacing the old value-tree-plus-one-symbolic-"
                "derivative-tree-per-variable walks IFeval() used to do. "
                "Covers every PT_*/PTF_* node type/function the parser "
                "supports (38 PTF_* funcnums plus PT_TIME/TEMPERATURE/"
                "FREQUENCY/TERN and the arithmetic operators -- see the "
                "round-4 report's coverage table); falls back to the "
                "existing symbolic + 0104-cache path for any tree with more "
                "than PT_AD_MAX_VARS=16 controlling variables (none observed "
                "in this repo) or when disabled (NGSPICE_ACCEL_NO_AD env var "
                "or -DNGSPICE_ACCEL_DISABLE_AD at compile time). Same "
                "arithmetic as the symbolic path (each dual-number operator "
                "calls the identical low-level PT*() functions the symbolic "
                "derivative tree would have), just fused into one pass -- "
                "measured bit-identical, not merely within tolerance.",
    },
}


def write_spiceinit(text, dirs):
    paths = []
    for d in dirs:
        p = pathlib.Path(d) / ".spiceinit"
        if text is None:
            if p.exists():
                p.unlink()
        else:
            p.write_text(text)
            paths.append(p)
    return paths


# --- raw-binary-output (ACCEPTED, card #467): not a baseline-vs-candidate
# binary/env comparison like the CANDIDATES table above -- it is a change to
# the bench deck itself (bench/decks/b_long_switch.cir now has 'set
# filetype=binary' where it used to say 'ascii'). The original same-results
# evidence was gathered once, by hand, comparing that deck's meas output
# before/after the edit (both runs on the unmodified baseline binary,
# aaba50ba3f49-s12) -- see bench/README.md and patches-wip/accel/README.md
# for the full numbers. This is the regression guard going forward: run
# today's (binary-format) deck and check its .meas values still match the
# reference captured at acceptance time.
RAW_BINARY_REFERENCE_MEAS = {
    # captured 2026-09-28, baseline binary aaba50ba3f49-s12, both the
    # pre-change (ascii) and post-change (binary) deck gave these values
    # bit-identically across 5+5 runs.
    "vout_avg_final": 3.300018e+00,
    "vout_max_final": 3.302050e+00,
    "vout_min_final": 3.297362e+00,
    "vout_ripple": 4.688000e-03,
    "il_max_final": 1.233484e+00,
    "il_min_final": 7.656706e-01,
    "duty_meas_final": 3.035135e-01,
    "fsw_meas": 5.000093e+05,
}


def test_raw_binary_output():
    meas, err = run_deck(BASELINE_BIN, "b_long_switch")
    if err:
        print(f"[raw-binary-output] b_long_switch {err}")
        return False
    problems = compare(RAW_BINARY_REFERENCE_MEAS, meas, tol_rel=0.0)
    if problems:
        print(f"[raw-binary-output] FAIL ({len(problems)} values changed since acceptance):")
        for p in problems[:10]:
            print(f"    {p}")
        return False
    print(f"[raw-binary-output] PASS ({len(RAW_BINARY_REFERENCE_MEAS)}/{len(RAW_BINARY_REFERENCE_MEAS)} reference values match)")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", help="run only this candidate, regardless of its status "
                                         "(accepted or rejected) -- the explicit-selection escape "
                                         "hatch for re-checking a rejected candidate on demand")
    ap.add_argument("--deck", help="run only this deck (within the chosen candidate(s))")
    ap.add_argument("--all", action="store_true",
                     help="also run rejected/pending candidates in the default (no --candidate) run. "
                          "Without this flag, a plain run only re-checks ACCEPTED candidates (card "
                          "#478 T9.2 task 0: a rejected candidate is a closed, documented result -- "
                          "see bench/README.md and patches-wip/accel/README.md for why -- not "
                          "something the default run should re-run every time, and several rejected "
                          "candidates need --candidate-bin, which a default run never has, so "
                          "including them by default made every plain run COULD-NOT-RUN/FAIL for a "
                          "reason unrelated to correctness. Candidates with status \"folded\" (card "
                          "#483 T9.3: accepted work later merged into `aug` itself, e.g. "
                          "accel-0103-omp-schedule-guided -> patches/0015) are excluded from --all "
                          "too, for the same reason -- there is no longer a separate candidate "
                          "binary to compare against; use --candidate <name> --candidate-bin <path> "
                          "explicitly if one is ever needed again.")
    ap.add_argument("--candidate-bin", help="ngspice-aug path for candidates needing a rebuilt binary "
                                             "(omp-schedule-runtime, march-native, and any accel-0100+ "
                                             "candidate whose bin is still None in CANDIDATES)")
    args = ap.parse_args()

    overall_ok = True

    if args.candidate == "raw-binary-output":
        cands = {}
    elif args.candidate:
        cands = {args.candidate: CANDIDATES[args.candidate]}
    elif args.all:
        cands = {name: cfg for name, cfg in CANDIDATES.items() if cfg.get("status") != "folded"}
    else:
        cands = {name: cfg for name, cfg in CANDIDATES.items() if cfg.get("status") == "accepted"}

    for cname, cfg in cands.items():
        cbin = cfg["bin"] or args.candidate_bin
        if cbin is None:
            print(f"[{cname}] COULD-NOT-RUN: needs --candidate-bin (no baseline-binary-only form)")
            overall_ok = False
            continue
        decks = [args.deck] if args.deck else cfg["decks"]
        print(f"=== candidate: {cname} [{cfg.get('status', 'pending')}] "
              f"({cfg['note'][:70]}...) tol_rel={cfg['tol_rel']} ===")
        for deck in decks:
            if deck not in DECKS:
                print(f"  [{deck}] COULD-NOT-RUN: unknown deck")
                overall_ok = False
                continue
            _, cwd, _ = DECKS[deck]
            spiceinit_dirs = [cwd]
            try:
                if cfg["spiceinit"] is not None:
                    write_spiceinit(cfg["spiceinit"], spiceinit_dirs)
                base_meas, base_err = run_deck(BASELINE_BIN, deck)
                if base_err:
                    print(f"  [{deck}] baseline {base_err}")
                    overall_ok = False
                    continue
                # candidate run: spiceinit-only candidates use the baseline binary;
                # binary-only candidates (omp-schedule-runtime, march-native) must
                # NOT pick up a stray spiceinit from a spiceinit-candidate's leftover.
                if cfg["spiceinit"] is None:
                    write_spiceinit(None, spiceinit_dirs)
                cand_meas, cand_err = run_deck(cbin, deck, env=cfg["env"])
            finally:
                write_spiceinit(None, spiceinit_dirs)
            if cand_err:
                print(f"  [{deck}] candidate {cand_err}")
                overall_ok = False
                continue
            problems = compare(base_meas, cand_meas, cfg["tol_rel"])
            if problems:
                print(f"  [{deck}] FAIL ({len(problems)} of {len(base_meas)} values out of tolerance):")
                for p in problems[:10]:
                    print(f"      {p}")
                overall_ok = False
            else:
                print(f"  [{deck}] PASS ({len(base_meas)}/{len(base_meas)} .meas values within tol)")

    if not args.candidate or args.candidate == "raw-binary-output":
        print("=== candidate: raw-binary-output (ACCEPTED, deck-level, see patches-wip/accel/README.md) ===")
        if not test_raw_binary_output():
            overall_ok = False

    print()
    print("OVERALL:", "PASS" if overall_ok else "FAIL")
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
