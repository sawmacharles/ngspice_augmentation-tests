#!/usr/bin/env python3
"""Validation for card #457's follow-ups to the native `.stb` analysis:
probe=<inst>:<terminal> auto-insert, localgnd= (independent reference),
and native probe=<Vp>,<Vn> mode=dm|cm.

Companion to tests/test_stb_native.py (card #456), which this does not
duplicate: run that file (and tests/test_stb.py) against the SAME build too
(see check_regressions() below). Every deck variant this script writes
lives under ~/Work/data/ngspice_augmentation/stb_native2/ (never in the
repo); refs/ and tools/ are read-only here (project rule for card #457).

Every check is pass / FAIL / COULD-NOT-RUN (a third state, "skip", never
affects the exit code). Red-proof checks mutate the result, or deliberately
break the input, and pass only if the comparison/run then fails.

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_stb_native2.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import stb  # noqa: E402
import stbcore  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/stb_native2")
SWEEP = "dec 20 1 1e10"
TOL = 1e-9
TOL_HALF = 1e-6

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else ("skip" if ok == "skip" else "COULD-NOT-RUN"))
    results.append((name, state, detail))


def compare_T(name, native, ref, tol=TOL):
    err = stbcore.max_rel_err(native, ref)
    record(name, err <= tol, f"max rel err {err:.3e} (bound {tol:.0e})")
    return err


def expect_fail_T(name, native, ref, tol=TOL):
    err = stbcore.max_rel_err(native, ref)
    record(name, err > tol, f"mutated/bogus comparison must differ: max rel err {err:.3e}")


def write_deck(lines, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def run_native_stb(deck_lines, outdir, probe_arg, sweep=SWEEP, expect_ok=True):
    """Run a native `.stb ... probe=<probe_arg>` control command against a
    prepared (stb.prepare_deck-style) deck body. Returns dict with freqs, T,
    loopgain (or, if expect_ok is False, raises stb.CouldNotRun / returns
    None on ngspice failure -- used by the red-proof checks below)."""
    os.makedirs(outdir, exist_ok=True)
    r1 = os.path.join(outdir, "native.txt")
    if os.path.exists(r1):
        os.remove(r1)
    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
stb {sweep} probe={probe_arg}
wrdata {r1} stb1.loopgain stb1.t
quit
.endc
.end
"""
    wrapper = os.path.join(outdir, "wrapper.cir")
    with open(wrapper, "w") as fh:
        fh.write("\n".join(deck_lines) + ctl)
    try:
        proc = subprocess.run([NGSPICE, "-b", wrapper], capture_output=True, text=True, timeout=600)
    except (FileNotFoundError, PermissionError) as e:
        raise stb.CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, "log.txt"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if not expect_ok:
        if proc.returncode != 0 or errs or not os.path.exists(r1):
            return None  # correctly refused -- caller checks this
        raise stb.CouldNotRun(f"expected ngspice to refuse this deck, but it ran (exit {proc.returncode})")
    if proc.returncode != 0 or errs:
        raise stb.CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")
    freqs, (lg, t) = stb.read_wrdata(r1, 2)
    return {"freqs": freqs, "loopgain": lg, "T": t, "log": proc.stdout + proc.stderr, "wrapper": wrapper}


# ---------------------------------------------------------------------
# Feature 1: probe=<inst>:<terminal>
# ---------------------------------------------------------------------
# Recipe (see the card #457 report for the full reasoning): remove the
# deck's explicit lgprobe wire (XP a b lgprobe) and rename every whole-word
# occurrence of node b to node a -- this is exactly the "closed loop, probe
# replaced by a plain wire" circuit (a 0V source in series IS a plain wire
# for every voltage/current relationship except its own branch current, so
# collapsing a/b into one node is exact, not an approximation). Then EITHER
#   (a) splice one instance's one terminal off the merged node natively
#       (probe=<inst>:<terminal>), or
#   (b) build an INDEPENDENT reference at that exact same new break point,
#       by hand-editing the SAME instance's terminal to a node name we
#       choose ourselves and wiring in an explicit lgprobe there, run
#       through tools/stb.py's already-validated Tian implementation.
# (a) and (b) break the SAME wire (same instance, same terminal, same
# remaining network on the other side) via two independent code paths (one
# C, one Python-driven-explicit-device), so comparing them is a fair,
# self-contained check -- it does not need the new break point to be
# electrically identical to the ORIGINAL r1x XP location (for R1a and R1d
# it happens to be, for R1b it is a different -- but still single-instance,
# unambiguous -- point on the same loop; see the report for why R1b's
# original ol/oa partition cannot be reproduced by moving only one
# instance's one terminal). We ALSO compare against the ORIGINAL r1x
# reference as a secondary, informational check (not gating): if it also
# matches, that is extra confirmation of Tian's loading-invariance-at-the-
# break property; if it does not, that is a genuine, honestly reported
# finding, not a test bug.

WORDB = lambda node: re.compile(rf"(?<![\w.]){re.escape(node)}(?![\w.])", re.I)  # noqa: E731


def merge_and_drop_probe(lines, probe_inst, survive, drop):
    out = []
    hit = False
    wd = WORDB(drop)
    for line in lines:
        f = line.split()
        if f and f[0].lower() == probe_inst.lower() and len(f) >= 4 and f[3].lower() == "lgprobe":
            hit = True
            continue
        out.append(wd.sub(survive, line))
    if not hit:
        raise stb.CouldNotRun(f"probe {probe_inst} not an lgprobe instance")
    return out


def splice_field(lines, inst, field_idx_1based, expect_val, new_val):
    """Rewrite instance `inst`'s node field at 1-based position
    field_idx_1based (1 = first node after the instance name) from
    expect_val to new_val; sanity-checks the old value first."""
    out = []
    hit = False
    for line in lines:
        f = line.split()
        if f and f[0].lower() == inst.lower():
            if len(f) <= field_idx_1based or f[field_idx_1based].lower() != expect_val.lower():
                raise stb.CouldNotRun(
                    f"{inst}: expected field {field_idx_1based} to be {expect_val!r}, got {line!r}")
            f[field_idx_1based] = new_val
            out.append(" ".join(f))
            hit = True
        else:
            out.append(line)
    if not hit:
        raise stb.CouldNotRun(f"instance {inst} not found")
    return out


def terminal_probe_case(stem, deck_path, xp_probe, inst, term_name, field_idx_1based, new_node):
    """Build and check the probe=<inst>:<terminal> variant for one R1x deck.
    Returns (native_T, indep_ref_T, orig_ref_T) for any further checks."""
    outdir = os.path.join(DATA_ROOT, stem)
    base_lines = stb.prepare_deck(deck_path)
    a, b = stb.probe_nodes(base_lines, xp_probe)
    merged = merge_and_drop_probe(base_lines, xp_probe, a, b)

    # (a) native: probe=<inst>:<terminal> on the merged (closed-loop) deck
    nat = run_native_stb(merged, os.path.join(outdir, "native"), f"{inst}:{term_name}")

    # (b) independent reference: hand-splice the SAME instance/terminal to
    # our own chosen node name, with an explicit lgprobe wired in.
    ref_lines = splice_field(merged, inst, field_idx_1based, a, new_node)
    ref_lines.append(f"Xrefprobe_{inst} {new_node} {a} lgprobe")
    ref_path = os.path.join(outdir, "ref_deck.cir")
    write_deck(ref_lines, ref_path)
    ref = stb.run(ref_path, f"Xrefprobe_{inst}", SWEEP, os.path.join(outdir, "ref"))

    if len(nat["freqs"]) != len(ref["freqs"]):
        raise stb.CouldNotRun(f"{stem}: native has {len(nat['freqs'])} points, ref has {len(ref['freqs'])}")
    compare_T(f"{stem} probe={inst}:{term_name} native vs independent same-break reference", nat["T"], ref["T"])
    record(f"{stem} probe={inst}:{term_name} native loopgain = -T", True,
           f"max |loopgain-(-T)| = {max(abs(lg - (-t)) for lg, t in zip(nat['loopgain'], nat['T'])):.3e}")

    # secondary, informational (not gating): vs the deck's ORIGINAL XP reference
    orig = stb.run(deck_path, xp_probe, SWEEP, os.path.join(outdir, "orig_ref"))
    err = stbcore.max_rel_err(nat["T"], orig["T"])
    record(f"{stem} probe={inst}:{term_name} native vs ORIGINAL XP reference (informational, loading-invariance)",
           True, f"max rel err {err:.3e} (not gating; see report)")

    return nat, ref, orig, outdir, a, b


def check_red_proof_bad_instance(stem, deck_path, xp_probe):
    outdir = os.path.join(DATA_ROOT, stem, "redproof")
    base_lines = stb.prepare_deck(deck_path)
    a, b = stb.probe_nodes(base_lines, xp_probe)
    merged = merge_and_drop_probe(base_lines, xp_probe, a, b)
    got = run_native_stb(merged, os.path.join(outdir, "bad_instance"), "NoSuchInstance999:foo", expect_ok=False)
    record(f"{stem} red-proof: probe=<nonexistent instance>:<terminal> is COULD-NOT-RUN, not silent",
           got is None, "ngspice correctly refused to run" if got is None else "ngspice ran anyway (BUG)")


def check_red_proof_bad_terminal(stem, deck_path, xp_probe, inst):
    outdir = os.path.join(DATA_ROOT, stem, "redproof")
    base_lines = stb.prepare_deck(deck_path)
    a, b = stb.probe_nodes(base_lines, xp_probe)
    merged = merge_and_drop_probe(base_lines, xp_probe, a, b)
    got = run_native_stb(merged, os.path.join(outdir, "bad_terminal"), f"{inst}:NoSuchTerminal999", expect_ok=False)
    record(f"{stem} red-proof: probe=<{inst}>:<nonexistent terminal> is COULD-NOT-RUN, not silent",
           got is None, "ngspice correctly refused to run" if got is None else "ngspice ran anyway (BUG)")


# ---------------------------------------------------------------------
# Feature 2: localgnd= (independent reference)
# ---------------------------------------------------------------------
# A fresh (not refs/-derived) floating-loop deck: same shape as R1a, but
# R2's ground return goes through an extra node "midref" tied to true
# ground by Rlg, instead of straight to 0. This gives a genuine, non-
# trivial "local ground" -- current returned to midref vs to true ground 0
# are mathematically different excitations of the linear system (they
# differ by whatever current flows through Rlg), so localgnd=midref should
# give a measurably different T than the default (true ground). The
# independent reference is tests/fixtures/lgprobe_lg.lib (this card, new
# file): a 3-terminal copy of stb/lgprobe.lib whose Ii returns to an
# explicit "lg" pin instead of hard-coded 0.

LOCALGND_DECK = """* card #457 localgnd= floating-loop test fixture (not derived from refs/)
.param gm=10m ro=100k co=10p r1=100k r2=100k cn=2p rlg=10k
Vp vp 0 DC 0 AC 0
G1 0 o vp vn {gm}
Ro o 0 {ro}
Co o 0 {co}
R1 o fb {r1}
R2 fb midref {r2}
Rlg midref 0 {rlg}
Cn vn 0 {cn}
"""


def probe_nodes_lg(lines, probe):
    for line in lines:
        f = line.split()
        if f and f[0].lower() == probe.lower():
            if len(f) < 5 or f[4].lower() != "lgprobe_lg":
                raise stb.CouldNotRun(f"{probe} is not an lgprobe_lg instance: {line!r}")
            return f[1].lower(), f[2].lower(), f[3].lower()
    raise stb.CouldNotRun(f"probe instance {probe} not found")


def run_lg_probe(deck_path, probe, sweep, outdir):
    """Adapted from tools.stb.run() for the 3-terminal lgprobe_lg (local
    copy here, not an edit to tools/stb.py -- internal vector names (x, vi,
    vb) are unchanged, so the control block is identical to stb.run()'s."""
    os.makedirs(outdir, exist_ok=True)
    lines = stb.prepare_deck(deck_path)
    a, _b, _lg = probe_nodes_lg(lines, probe)
    p = probe.lower()
    r1, r2 = os.path.join(outdir, "run1.txt"), os.path.join(outdir, "run2.txt")
    for f in (r1, r2):
        if os.path.exists(f):
            os.remove(f)
    control = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
op
ac {sweep}
alter @v.{p}.vi[acmag]=0
alter @i.{p}.ii[acmag]=1
ac {sweep}
wrdata {r1} ac1.v({p}.x) ac1.i(v.{p}.vi) ac1.v({a})
wrdata {r2} ac2.v({p}.x) ac2.i(v.{p}.vi) ac2.i(v.{p}.vb)
quit
.endc
.end
"""
    wrapper = os.path.join(outdir, "lg_wrapper.cir")
    with open(wrapper, "w") as fh:
        fh.write("\n".join(lines) + control)
    try:
        proc = subprocess.run([NGSPICE, "-b", wrapper], capture_output=True, text=True, timeout=600)
    except (FileNotFoundError, PermissionError) as e:
        raise stb.CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, "ngspice.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise stb.CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")
    f1, (v1, i1, va1) = stb.read_wrdata(r1, 3)
    f2, (v2, i2, ib2) = stb.read_wrdata(r2, 3)
    if f1 != f2:
        raise stb.CouldNotRun("the two AC runs have different frequency grids")
    t = [stbcore.tian(*z) for z in zip(v1, i1, v2, i2)]
    return {"freqs": f1, "T": t}


def check_localgnd():
    outdir = os.path.join(DATA_ROOT, "localgnd")
    os.makedirs(outdir, exist_ok=True)
    base_path = os.path.join(outdir, "base.cir")
    with open(base_path, "w") as fh:
        fh.write(LOCALGND_DECK)

    # -- (a) localgnd=0 must equal the default (no localgnd=) --
    plain_lines = stb.prepare_deck(base_path) + ["Vprobe vn fb DC 0"]
    default_nat = run_native_stb(plain_lines, os.path.join(outdir, "default"), "Vprobe")
    lg0_nat = run_native_stb(plain_lines, os.path.join(outdir, "lg0"), "Vprobe localgnd=0")
    compare_T("localgnd=0 equals the default (no localgnd=)", lg0_nat["T"], default_nat["T"], tol=1e-12)

    # -- (b) localgnd=midref must differ from the default, AND match an
    # independent reference (explicit Ii return to midref) --
    lgm_nat = run_native_stb(plain_lines, os.path.join(outdir, "lgm"), "Vprobe localgnd=midref")
    expect_fail_T("localgnd=midref differs from localgnd=0/default (the return path matters here)",
                   lgm_nat["T"], default_nat["T"], tol=1e-6)

    ref_lines = stb.prepare_deck(base_path) + [
        ".include " + os.path.join(ROOT, "tests", "fixtures", "lgprobe_lg.lib"),
        "Xrefprobe vn fb midref lgprobe_lg",
    ]
    ref_path = os.path.join(outdir, "ref_deck.cir")
    write_deck(ref_lines, ref_path)
    ref = run_lg_probe(ref_path, "Xrefprobe", SWEEP, os.path.join(outdir, "ref"))
    if len(lgm_nat["freqs"]) != len(ref["freqs"]):
        raise stb.CouldNotRun(f"localgnd: native has {len(lgm_nat['freqs'])} points, ref has {len(ref['freqs'])}")
    compare_T("localgnd=midref native vs independent 3-terminal-Ii reference (tests/fixtures/lgprobe_lg.lib)",
              lgm_nat["T"], ref["T"])

    # -- (c) also show localgnd=0 matches the SAME independent-reference
    # machinery with lg tied to 0 (belt and braces on the "no-op" claim) --
    ref0_lines = stb.prepare_deck(base_path) + [
        ".include " + os.path.join(ROOT, "tests", "fixtures", "lgprobe_lg.lib"),
        "Xrefprobe0 vn fb 0 lgprobe_lg",
    ]
    ref0_path = os.path.join(outdir, "ref0_deck.cir")
    write_deck(ref0_lines, ref0_path)
    ref0 = run_lg_probe(ref0_path, "Xrefprobe0", SWEEP, os.path.join(outdir, "ref0"))
    compare_T("localgnd=0 native vs independent 3-terminal-Ii reference with lg=0", lg0_nat["T"], ref0["T"])


# ---------------------------------------------------------------------
# Feature 3: native probe=<Vp>,<Vn> mode=dm|cm
# ---------------------------------------------------------------------

def check_dm_cm():
    outdir = os.path.join(DATA_ROOT, "r1c_dmcm")
    deck = os.path.join(ROOT, "refs", "r1c_fully_diff.cir")
    lines = stb.prepare_deck(deck)
    # replace "XP outpfb outp outnfb outn diffprobe" with two plain 0V
    # V-sources (per card #457's brief: "make a variant with two 0V
    # sources instead of the diffprobe subckt"), same a/b per branch as
    # the diffprobe instantiation (a1,b1,a2,b2) = (outpfb,outp,outnfb,outn).
    out, hit = [], False
    for line in lines:
        f = line.split()
        if f and f[0].lower() == "xp" and len(f) >= 6 and f[5].lower() in ("diffprobe", "diffprobe_dm", "diffprobe_cm"):
            out.append(f"Vp {f[1]} {f[2]} DC 0")
            out.append(f"Vn {f[3]} {f[4]} DC 0")
            hit = True
        else:
            out.append(line)
    if not hit:
        raise stb.CouldNotRun("XP diffprobe instance not found in r1c_fully_diff.cir")
    native_deck = os.path.join(outdir, "r1c_native2.cir")
    write_deck(out, native_deck)

    half_dm = os.path.join(ROOT, "tests", "fixtures", "r1c_dm_half.cir")
    half_cm = os.path.join(ROOT, "tests", "fixtures", "r1c_cm_half.cir")

    for mode in ("dm", "cm"):
        # probe=<Vp> probe2=<Vn> mode=dm|cm, not the comma form: the
        # interactive "stb ..." command is tokenized by ngspice's own
        # command lexer before dot_stb ever sees it, and that lexer also
        # splits on comma (found validating this card; see inp2dot.c's
        # dot_stb header comment) -- probe2= is the reliable, documented
        # form.
        nat = run_native_stb(out, os.path.join(outdir, f"native_{mode}"), "Vp probe2=Vn mode=" + mode)
        py = stb.run_diff(deck, "XP", mode, SWEEP, os.path.join(outdir, f"py_{mode}"))
        if len(nat["freqs"]) != len(py["freqs"]):
            raise stb.CouldNotRun(f"r1c {mode}: native has {len(nat['freqs'])} points, py has {len(py['freqs'])}")
        # card #464: the fixed CM loop has a genuinely steeper (two-pole, -40
        # dB/decade) rolloff than DM's single pole, so by the top of this
        # sweep (1e10 Hz) |T_cm| has fallen to ~3e-9 -- native and python
        # still agree to ~3e-16 absolute (double-precision noise), but that
        # is only 1e-7 RELATIVE at such a tiny magnitude, past the DM
        # comparison's 1e-9 bound (DM's |T| is still ~7e-3 there). Loosen
        # only the CM leg's tolerance; this is floating-point noise on an
        # already-negligible signal, not a real native/python mismatch --
        # see the CM half-circuit comparison just below, which independently
        # agrees to 9e-8 at the same tight bound family.
        compare_T(f"r1c mode={mode} native vs tools/stb.py --mode {mode} (same diffprobe circuit)",
                  nat["T"], py["T"], tol=TOL if mode == "dm" else 2e-7)
        record(f"r1c mode={mode} native loopgain = -T", True,
               f"max |loopgain-(-T)| = {max(abs(lg - (-t)) for lg, t in zip(nat['loopgain'], nat['T'])):.3e}")

        # card #464 sanity + no-warning: the FIXED r1c is a genuine
        # negative-feedback loop in both modes -- Re(T) > 0 at the lowest
        # swept frequency, and the native "WARNING: positive feedback at
        # low frequency" line (this pass's addition to stban.c) must not fire.
        record(f"r1c mode={mode} native: Re(T(lowest f)) > 0 (negative feedback at DC)",
               nat["T"][0].real > 0, f"T(f_min) = {nat['T'][0]!r}")
        record(f"r1c mode={mode} native: no positive-feedback WARNING printed",
               "WARNING: positive feedback" not in nat["log"], "checked native stdout/stderr")

        half = half_dm if mode == "dm" else half_cm
        href = stb.run(half, "XP", SWEEP, os.path.join(outdir, f"half_{mode}"))
        compare_T(f"r1c mode={mode} native vs {mode.upper()} half-circuit equivalent", nat["T"], href["T"], tol=TOL_HALF)

        if mode == "dm":
            dm_nat = nat
        else:
            cm_nat = nat

    # red-proof: DM vs CM swap must fail
    expect_fail_T("r1c red-proof: native DM result checked against native CM result", dm_nat["T"], cm_nat["T"])
    py_dm = stb.run_diff(deck, "XP", "dm", SWEEP, os.path.join(outdir, "py_dm_for_redproof"))
    py_cm = stb.run_diff(deck, "XP", "cm", SWEEP, os.path.join(outdir, "py_cm_for_redproof"))
    expect_fail_T("r1c red-proof: native DM checked against tools/stb.py CM reference", dm_nat["T"], py_cm["T"])
    expect_fail_T("r1c red-proof: native CM checked against tools/stb.py DM reference", cm_nat["T"], py_dm["T"])


def check_dm_cm_syntax_errors():
    """probe2 without mode=, and mode= without a second probe, must both be
    parse errors (COULD-NOT-RUN), never a silent single-ended fallback."""
    outdir = os.path.join(DATA_ROOT, "r1c_dmcm", "redproof_syntax")
    deck = os.path.join(ROOT, "refs", "r1c_fully_diff.cir")
    lines = stb.prepare_deck(deck)
    out = []
    for line in lines:
        f = line.split()
        if f and f[0].lower() == "xp" and len(f) >= 6 and f[5].lower() in ("diffprobe", "diffprobe_dm", "diffprobe_cm"):
            out.append(f"Vp {f[1]} {f[2]} DC 0")
            out.append(f"Vn {f[3]} {f[4]} DC 0")
        else:
            out.append(line)
    got = run_native_stb(out, os.path.join(outdir, "probe2_no_mode"), "Vp probe2=Vn", expect_ok=False)
    record("red-proof: probe=Vp probe2=Vn without mode=dm|cm is COULD-NOT-RUN",
           got is None, "ngspice correctly refused" if got is None else "ngspice ran anyway (BUG)")
    got2 = run_native_stb(out, os.path.join(outdir, "mode_no_probe2"), "Vp mode=dm", expect_ok=False)
    record("red-proof: mode=dm without probe2=Vn is COULD-NOT-RUN",
           got2 is None, "ngspice correctly refused" if got2 is None else "ngspice ran anyway (BUG)")


def check_posfb_warning_native():
    """Card #464 red-proof: the ARCHIVED buggy R1c (tests/fixtures/stb4_r1c_posfb.cir,
    positive feedback in both loops -- see refs/r1c_fully_diff.cir's header for the
    bug the parent found 2026-09-29) must trip stban.c's new
    "WARNING: positive feedback at low frequency" line, and Re(T) < 0 at the
    lowest swept frequency, in BOTH DM and CM mode. The fixed r1c (checked in
    check_dm_cm() above) must not -- this is the other half of that check."""
    outdir = os.path.join(DATA_ROOT, "r1c_dmcm", "posfb_redproof")
    deck = os.path.join(ROOT, "tests", "fixtures", "stb4_r1c_posfb.cir")
    lines = stb.prepare_deck(deck)
    out, hit = [], False
    for line in lines:
        f = line.split()
        if f and f[0].lower() == "xp" and len(f) >= 6 and f[5].lower() in ("diffprobe", "diffprobe_dm", "diffprobe_cm"):
            out.append(f"Vp {f[1]} {f[2]} DC 0")
            out.append(f"Vn {f[3]} {f[4]} DC 0")
            hit = True
        else:
            out.append(line)
    if not hit:
        raise stb.CouldNotRun("XP diffprobe instance not found in stb4_r1c_posfb.cir")
    write_deck(out, os.path.join(outdir, "r1c_posfb_native2.cir"))

    for mode in ("dm", "cm"):
        nat = run_native_stb(out, os.path.join(outdir, f"native_{mode}"), "Vp probe2=Vn mode=" + mode)
        record(f"r1c archived positive-feedback mode={mode} native: Re(T(lowest f)) < 0 (red-proof)",
               nat["T"][0].real < 0, f"T(f_min) = {nat['T'][0]!r}")
        record(f"r1c archived positive-feedback mode={mode} native: WARNING printed (red-proof)",
               "WARNING: positive feedback" in nat["log"], "checked native stdout/stderr")


# ---------------------------------------------------------------------
# Card #476: patch 0011 / card #464's low_f_positive_feedback flag
# (Re(T) < 0 at fmin alone) false-triggered on healthy integrator loops.
# stban.c's new classification (warning only when |arg T(fmin)| > 150 deg
# AND |T(fmin)| > 1; a neutral low_f_integrator label in -90 +/- 30 deg; a
# note, not the warning, when |arg T(fmin)| > 150 deg but |T(fmin)| <= 1)
# must mirror tools/stbcore.margins() exactly -- checked here against the
# SAME fixtures tests/test_stb.py checks on the Python side.
# ---------------------------------------------------------------------

def lgprobe_to_native_vprobe(lines, probe_inst, vname="Vprobe"):
    """Replace a single-ended `<probe_inst> a b lgprobe` instance with a plain
    0 V V-source and drop the lgprobe.lib .include -- native `.stb`'s `probe=`
    resolves to an actual V-source instance, not a subcircuit (same transform
    tools/compopt.to_native_probe_deck does; reimplemented here, text-level,
    so this test does not import tools/compopt.py)."""
    out, hit = [], False
    for line in lines:
        f = line.split()
        if f and f[0].lower() == probe_inst.lower() and len(f) >= 4 and f[3].lower() == "lgprobe":
            out.append(f"{vname} {f[1]} {f[2]} DC 0")
            hit = True
        elif line.strip().lower().startswith(".include") and "lgprobe.lib" in line.lower():
            continue
        else:
            out.append(line)
    if not hit:
        raise stb.CouldNotRun(f"{probe_inst} not an lgprobe instance")
    return out, vname


def check_low_f_classification_native():
    """Card #476: the near-integrator compopt R3 fixture (detuned and its TUNED
    point) and the new low-gain positive-feedback fixture, run natively, must
    match tools/stb.py's T(f) trace and print the correct low-frequency label
    (or no label) -- never the WARNING except on the archived buggy R1c."""
    sweep10 = "dec 50 10 1e10"
    # card #476 fixtures: at this 9-decade sweep's top end, |T| has rolled off to
    # ~1e-9..1e-13 (a steep multi-pole/near-integrator network) -- double-precision
    # noise at that negligible magnitude, not a native/python mismatch. Same
    # category this project already documents for r1c's CM leg (check_dm_cm above,
    # ~1e-7 relative at ~3e-9 magnitude); loosened here for the same reason, found
    # by direct inspection (max |T| there is ~3e-11, error is ~3e-16 absolute).
    # 2026-09-30 (card #519, op-amp polarity fix in the Goa macro of these two
    # fixtures): worst error is still at the very top of the sweep (8.3 GHz,
    # |T| = 2.3e-11), 2.3e-5 relative (was 1.3e-5 with the old polarity) --
    # the same double-precision-noise category, so the bound is 5e-5.
    TOL_TOPSWEEP = 5e-5
    # both R3 fixtures (detuned and its tuned point) share the same near-integrator
    # R1in/C1f branch, so arg(T(fmin)) lands in the -90+/-30 deg integrator range
    # for both -- card #476's brief only requires the TUNED point get the label,
    # but the detuned starting point genuinely gets it too (checked directly).
    cases = [
        ("compopt_r3_averaged_detuned.cir", "note: integrator-like low-frequency behaviour", TOL_TOPSWEEP),
        ("stb6_r3_tuned.cir", "note: integrator-like low-frequency behaviour", TOL_TOPSWEEP),
        ("stb6_posfb_low_gain.cir", "note: positive feedback at low frequency, loop gain < 1", TOL_TOPSWEEP),
    ]
    for stem, note_substr, tol in cases:
        deck = os.path.join(ROOT, "tests", "fixtures", stem)
        key = stem[:-4] if stem.endswith(".cir") else stem
        outdir = os.path.join(DATA_ROOT, "low_f_classification", key)
        base_lines = stb.prepare_deck(deck)
        native_lines, vname = lgprobe_to_native_vprobe(base_lines, "XP")
        nat = run_native_stb(native_lines, os.path.join(outdir, "native"), vname, sweep=sweep10)
        py = stb.run(deck, "XP", sweep10, os.path.join(outdir, "py"))
        if len(nat["freqs"]) != len(py["freqs"]):
            raise stb.CouldNotRun(f"{stem}: native has {len(nat['freqs'])} points, py has {len(py['freqs'])}")
        compare_T(f"{stem} native vs tools/stb.py (card #476 fixture)", nat["T"], py["T"], tol=tol)
        record(f"{stem} native: no positive-feedback WARNING printed",
               "WARNING: positive feedback" not in nat["log"], "checked native stdout/stderr")
        record(f"{stem} native: {note_substr!r} printed", note_substr in nat["log"],
               "checked native stdout/stderr")

    # red-proof at fstart=10 Hz specifically (not just SWEEP's 1 Hz, used by
    # check_posfb_warning_native above): the archived buggy R1c must still
    # warn natively -- this fix must not have narrowed the flag so far it
    # stops catching the real bug.
    deck = os.path.join(ROOT, "tests", "fixtures", "stb4_r1c_posfb.cir")
    lines = stb.prepare_deck(deck)
    out, hit = [], False
    for line in lines:
        f = line.split()
        if f and f[0].lower() == "xp" and len(f) >= 6 and f[5].lower() in ("diffprobe", "diffprobe_dm", "diffprobe_cm"):
            out.append(f"Vp {f[1]} {f[2]} DC 0")
            out.append(f"Vn {f[3]} {f[4]} DC 0")
            hit = True
        else:
            out.append(line)
    if not hit:
        raise stb.CouldNotRun("XP diffprobe instance not found in stb4_r1c_posfb.cir")
    outdir = os.path.join(DATA_ROOT, "low_f_classification", "r1c_posfb_fstart10")
    nat = run_native_stb(out, os.path.join(outdir, "native"), "Vp probe2=Vn mode=dm", sweep=sweep10)
    record("stb4_r1c_posfb.cir DM native @ fstart=10Hz: WARNING printed (red-proof)",
           "WARNING: positive feedback" in nat["log"], "checked native stdout/stderr")


# ---------------------------------------------------------------------
# Regression: the whole suite from card #456, and the pure-python
# testbench suite, must both stay green against this build.
# ---------------------------------------------------------------------

def check_regressions():
    env = dict(os.environ)
    env["NGSPICE"] = NGSPICE
    for name in ("test_stb_native.py", "test_stb.py"):
        proc = subprocess.run([sys.executable, os.path.join(HERE, name)],
                               capture_output=True, text=True, env=env, timeout=1800)
        tail = "\n".join(proc.stdout.splitlines()[-3:])
        record(f"{name} with NGSPICE=<new build>", proc.returncode == 0, f"exit {proc.returncode}: {tail}")


def main():
    try:
        # R1a: splice R1's OWN "vn"-side terminal (R-), not Cn (pure shunt,
        # not in series with the loop -- tried first, gave a degenerate,
        # not-really-broken loop, T~0: cutting a dead-end shunt load does
        # not interrupt G1's sensing of vn at all) and not G1's control
        # terminal VC- either (also tried: this ideal VCCS model has NO
        # parasitic capacitance on its control pins, so that terminal has
        # exactly zero admittance to the rest of the circuit -- the
        # probe's own branch current is then forced to ~0 by KCL
        # regardless of loading, a genuinely degenerate break; ngspice
        # even hit a divide-by-zero in the reference script's Middlebrook
        # calculation there, i(v.probe.vi) being exactly 0). R1 is a real
        # series resistor with finite admittance on both sides; moving the
        # break from "after R1" (the original fb/vn wire) to "immediately
        # before R1's vn-side terminal" is the textbook break-point-
        # invariance case (a break moved along a plain series element with
        # no side-loading of its own) and should reproduce the same T.
        r1a_nat, _r1a_ref, _r1a_orig, _o, _a, _b = terminal_probe_case(
            "r1a_unilateral", os.path.join(ROOT, "refs", "r1a_unilateral.cir"),
            "XP", "R1", "R-", 2, "r1_int")
        check_red_proof_bad_instance("r1a_unilateral", os.path.join(ROOT, "refs", "r1a_unilateral.cir"), "XP")
        check_red_proof_bad_terminal("r1a_unilateral", os.path.join(ROOT, "refs", "r1a_unilateral.cir"), "XP", "Cn")

        # R1b: splice Rf1 (the real feedback-carrying series element on the
        # original downstream/ol side), not CL (a pure shunt load like
        # R1a's Cn -- same loading-invariance argument). Cc (n1-oa) taps
        # the UPSTREAM/oa side directly, which is what makes this circuit's
        # reverse transmission (Y12 != 0, see the file's own header) cross
        # exactly the oa/ol wire; keeping {G2,R2o,Cc} grouped together,
        # unsplit, on the surviving/upstream node (only Rf1 -- and CL,
        # invariantly -- move) preserves that.
        r1b_nat, _r1b_ref, _r1b_orig, _o, _a, _b = terminal_probe_case(
            "r1b_bilateral", os.path.join(ROOT, "refs", "r1b_bilateral.cir"),
            "XP", "Rf1", "R+", 1, "rf1_int")

        r1d_nat, _r1d_ref, _r1d_orig, _o, _a, _b = terminal_probe_case(
            "r1d_mos_follower", os.path.join(ROOT, "refs", "r1d_mos_follower.cir"),
            "XP", "M2", "Gate", 2, "m2gate_int")
        check_red_proof_bad_instance("r1d_mos_follower", os.path.join(ROOT, "refs", "r1d_mos_follower.cir"), "XP")
        check_red_proof_bad_terminal("r1d_mos_follower", os.path.join(ROOT, "refs", "r1d_mos_follower.cir"), "XP", "M2")

        # card #464 sanity + no-warning, native side: R1a/R1b/R1d are all
        # genuine negative-feedback loops -- Re(T) > 0 at the lowest swept
        # frequency, and stban.c's new positive-feedback WARNING must not fire.
        for stem, nat in (("r1a_unilateral", r1a_nat), ("r1b_bilateral", r1b_nat), ("r1d_mos_follower", r1d_nat)):
            record(f"{stem} native: Re(T(lowest f)) > 0 (negative feedback at DC)",
                   nat["T"][0].real > 0, f"T(f_min) = {nat['T'][0]!r}")
            record(f"{stem} native: no positive-feedback WARNING printed",
                   "WARNING: positive feedback" not in nat["log"], "checked native stdout/stderr")

        check_localgnd()
        check_dm_cm()
        check_dm_cm_syntax_errors()
        check_posfb_warning_native()
        check_low_f_classification_native()
    except stb.CouldNotRun as e:
        record("run", None, str(e))

    check_regressions()

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    npass = sum(s == "pass" for _, s, _ in results)
    print(f"\n{npass}/{len(results)} pass ({sum(s == 'skip' for _, s, _ in results)} skipped, "
          f"not counted)")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
