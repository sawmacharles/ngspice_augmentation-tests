#!/usr/bin/env python3
"""Validation for the native `.stb` DEVICE mode (card #459, plan T3.4).

Device mode: `.stb ... probe=<MOSFET instance>` nulls the device's dominant
controlled source (gm alone, in this fork -- see stban.c's stb_mos1_probe()
doc comment) instead of breaking a wire, using the rank-1 form
T = 1/(1 - gm*v^T Y^-1 u) - 1 (research/cadence-stb-pstb.md sec 2.5). This
script compares that native result against an INDEPENDENT reference built
without any of this fork's C code: the MOSFET is replaced, in a plain text
copy of the deck, by its own small-signal equivalent (Rds, Cgs, Cgd, Cgb,
Cbd, Cbs, and one or two VCCS G-elements) using op-point values read back
from the SAME build via `@inst[param]` after `.op` -- then Bode's classical
return-ratio measurement is run on that equivalent with plain `.ac`:
disconnect the controlled source's control port, drive it with a 1 V AC
test source, and read back the voltage that would have driven it
(T = -v_returned / v_drive). See check_gm_only() below for the construction.

Every check is pass / FAIL / COULD-NOT-RUN, same convention and exit codes
as tests/test_stb_native.py (0 = all pass, 1 = any FAIL, 2 = any
COULD-NOT-RUN). All decks and logs this script writes live under
~/Work/data/ngspice_augmentation/stbdev/ (never in the repository); the
fixtures it reads are tests/fixtures/stbdev_*.cir and refs/r1d_mos_follower.cir
(read-only).

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_stb_device.py
"""
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXDIR = os.path.join(HERE, "fixtures")
sys.path.insert(0, os.path.join(ROOT, "tools"))
import stb  # noqa: E402  (CouldNotRun, prepare_deck, read_wrdata)
import stbcore  # noqa: E402  (max_rel_err and friends)

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/stbdev")
SWEEP = "dec 10 1 1e8"
TOL = 1e-6  # spec.md stb_device: 1e-6 relative

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else ("skip" if ok == "skip" else "COULD-NOT-RUN"))
    results.append((name, state, detail))


def run_ngspice(lines_or_text, outdir, tag, timeout=300):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{tag}.cir")
    text = lines_or_text if isinstance(lines_or_text, str) else "\n".join(lines_or_text) + "\n"
    with open(path, "w") as fh:
        fh.write(text)
    try:
        proc = subprocess.run([NGSPICE, "-b", path], capture_output=True, text=True, timeout=timeout)
    except (FileNotFoundError, PermissionError) as e:
        raise stb.CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, f"{tag}.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    return proc


def check_errors(proc, context):
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise stb.CouldNotRun(f"{context}: ngspice exit {proc.returncode}; {errs[:3]}")


# ---------- op-point extraction (gm, gmbs, gds, and the 5 capacitances) ----------

OP_PARAMS = ["gm", "gmbs", "gds", "cgs", "cgd", "cgb", "cbd", "cbs", "id"]


def op_extract(deck_lines, inst, outdir):
    """Run `.op` on the deck as-is and read back inst's small-signal fields."""
    prints = "\n".join(f"print @{inst}[{p}]" for p in OP_PARAMS)
    text = "\n".join(deck_lines) + f"\n.control\nop\n{prints}\nquit\n.endc\n.end\n"
    proc = run_ngspice(text, outdir, "opdump")
    check_errors(proc, f"op-point extraction for {inst}")
    out = {}
    for p in OP_PARAMS:
        m = re.search(rf"@{re.escape(inst)}\[{p}\]\s*=\s*([-+0-9.eE]+)", proc.stdout, re.I)
        if not m:
            raise stb.CouldNotRun(f"could not read @{inst}[{p}] from op-point dump")
        out[p] = float(m.group(1))
    return out


# ---------- brute-force small-signal-equivalent reference ----------

def build_brute_lines(deck_lines, inst, gate, drain, source, bulk, op, disconnect):
    """deck_lines with `inst` replaced by its small-signal equivalent.

    disconnect: "gm" (production default: null gm alone, gmbs -- if
    nonzero -- stays wired normally, sensing the real network's vbs) or
    "gmbs" (the opposite single-source return ratio, used only to show the
    gm-only choice is not a negligible one -- see check_body_effect()
    below). Returns (lines_without_.end, (v_probe_p, v_probe_n)).
    """
    gm, gmbs, gds, cgs, cgd, cgb, cbd, cbs, id_ = (op[k] for k in OP_PARAMS)
    out, removed = [], False
    for line in deck_lines:
        f = line.split()
        if f and f[0].lower() == inst.lower():
            removed = True
            continue
        out.append(line)
    if not removed:
        raise stb.CouldNotRun(f"{inst} not found in deck to replace")

    tag = inst
    elems = []
    # DC bias preservation (needed on a multi-transistor circuit like
    # refs/r1d_mos_follower.cir M2, where the REST of the circuit -- M1,
    # M3, M4, M5 there -- reaches ITS OWN operating point partly through
    # M2's real DC drain current; a bare Rds+VCCS replacement has no DC
    # current source of its own, so at DC it looks nearly open (gds is
    # tiny) and the surrounding active devices drift to a different,
    # wrong bias -- found empirically while validating this card. A
    # single-transistor fixture (no other active device depends on this
    # one's DC current) does not need this -- Idc's own contribution is
    # then just a harmless no-op DC current into a source/ground net --
    # but it is included unconditionally for correctness. AC=0 (DC only):
    # the small-signal gm/gmbs VCCS's below supply the AC current.
    elems.append(f"Idc_{tag} {drain} {source} DC {id_:.12g} AC 0")
    if gds > 0:
        elems.append(f"Rds_{tag} {drain} {source} {1.0 / gds:.12g}")
    if cgs > 0:
        elems.append(f"Cgsi_{tag} {gate} {source} {cgs:.12g}")
    if cgd > 0:
        elems.append(f"Cgdi_{tag} {gate} {drain} {cgd:.12g}")
    if cgb > 0 and gate.lower() != bulk.lower():
        elems.append(f"Cgbi_{tag} {gate} {bulk} {cgb:.12g}")
    if cbd > 0 and bulk.lower() != drain.lower():
        elems.append(f"Cbdi_{tag} {bulk} {drain} {cbd:.12g}")
    if cbs > 0 and bulk.lower() != source.lower():
        elems.append(f"Cbsi_{tag} {bulk} {source} {cbs:.12g}")

    # The disconnected control port is driven by an ideal 1V AC source
    # referenced to REAL ground (node "0"), not to a second floating node:
    # a V-source between two otherwise-unconnected nodes leaves their
    # common-mode potential with no path to ground at all, which is
    # singular (found empirically while validating this card -- "singular
    # matrix: check node ctrln_m1", a floating 2-node island with only a
    # V-source's OWN branch equation touching it, no other conductance to
    # pin its absolute potential). Referencing to "0" (always global
    # ground in SPICE, whatever the rest of the deck calls its own
    # reference net) fixes this and is exactly Bode's "drive it with 1V"
    # regardless of which node of the control pair is called "+" here.
    if disconnect == "gm":
        elems.append(f"Ggm_{tag} {drain} {source} ctrlp_{tag} 0 {gm:.12g}")
        elems.append(f"Vtest_{tag} ctrlp_{tag} 0 DC 0 AC 1")
        if gmbs != 0 and bulk.lower() != source.lower():
            elems.append(f"Ggmbs_{tag} {drain} {source} {bulk} {source} {gmbs:.12g}")
        v_probe = (gate, source)
    elif disconnect == "gmbs":
        if gmbs == 0 or bulk.lower() == source.lower():
            raise stb.CouldNotRun(f"{inst}: gmbs is 0 (or bulk==source) -- nothing to disconnect")
        elems.append(f"Ggmbs_{tag} {drain} {source} ctrlp_{tag} 0 {gmbs:.12g}")
        elems.append(f"Vtest_{tag} ctrlp_{tag} 0 DC 0 AC 1")
        if gm != 0:
            elems.append(f"Ggm_{tag} {drain} {source} {gate} {source} {gm:.12g}")
        v_probe = (bulk, source)
    else:
        raise ValueError(disconnect)

    return out + elems, v_probe


def run_brute_T(deck_lines, inst, gate, drain, source, bulk, op, disconnect, sweep, outdir, neutralize=()):
    """neutralize: hierarchical @elt[acmag] names to zero before the `ac`
    run -- for decks (like refs/r1d_mos_follower.cir) that carry their OWN
    independent AC source unrelated to this measurement (XP's internal
    lgprobe.lib Vi defaults to AC=1 -- stb/lgprobe.lib's own two-injection
    scheme normally alters it, but a plain brute-force `ac` here would
    otherwise add its 1V stimulus on top of Vtest's, corrupting the
    reading -- found empirically while validating this card on R1d)."""
    lines, (vp, vn) = build_brute_lines(deck_lines, inst, gate, drain, source, bulk, op, disconnect)
    r1 = os.path.join(outdir, "brute.txt")
    if os.path.exists(r1):
        os.remove(r1)
    # v(0) is not a valid ngspice vector expression (node 0/ground has no
    # printable vector) -- found empirically while validating this card.
    # Ground's voltage is always exactly 0, so just drop that term.
    vexpr = f"v({vp})" if vn == "0" else (f"-v({vn})" if vp == "0" else f"v({vp}) - v({vn})")
    alters = "\n".join(f"alter {n} = 0" for n in neutralize)
    text = ("\n".join(lines) +
            f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
            f"{alters}\n"
            f"ac {sweep}\nlet vret_{inst} = {vexpr}\n"
            f"let Tbrute_{inst} = -vret_{inst}\n"
            f"wrdata {r1} Tbrute_{inst}\nquit\n.endc\n.end\n")
    proc = run_ngspice(text, outdir, f"brute_{disconnect}")
    check_errors(proc, f"brute-force {disconnect} AC run for {inst}")
    freqs, (t,) = stb.read_wrdata(r1, 1)
    return freqs, t


def run_native_device_T(deck_lines, inst, sweep, outdir):
    r1 = os.path.join(outdir, "native.txt")
    if os.path.exists(r1):
        os.remove(r1)
    text = ("\n".join(deck_lines) +
            f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
            f"stb {sweep} probe={inst}\nwrdata {r1} stb1.t\nquit\n.endc\n.end\n")
    proc = run_ngspice(text, outdir, "native")
    check_errors(proc, f"native device-mode stb for {inst}")
    freqs, (t,) = stb.read_wrdata(r1, 1)
    return freqs, t


def compare_T(name, a, b, tol=TOL):
    if len(a) != len(b):
        raise stb.CouldNotRun(f"{name}: length mismatch {len(a)} vs {len(b)}")
    err = stbcore.max_rel_err(a, b)
    record(name, err <= tol, f"max rel err {err:.3e} (bound {tol:.0e})")
    return err


def expect_fail_T(name, a, b, tol=TOL):
    err = stbcore.max_rel_err(a, b)
    record(name, err > tol, f"mutated/bogus comparison must differ: max rel err {err:.3e}")


# ---------- checks ----------

def check_cs_cgd():
    """check 1a: single CS stage, large discrete Cgd local loop, no body
    effect. Independent brute-force gm-only reference vs native device-mode T."""
    stem = "cs_cgd"
    outdir = os.path.join(DATA_ROOT, stem)
    deck = stb.prepare_deck(os.path.join(FIXDIR, "stbdev_cs_cgd.cir"))
    op = op_extract(deck, "m1", os.path.join(outdir, "op"))
    record(f"{stem} op point sane (gm>0, gds>0)", op["gm"] > 0 and op["gds"] > 0,
           f"gm={op['gm']:.6g} gds={op['gds']:.6g} gmbs={op['gmbs']:.6g}")
    freqs_n, t_native = run_native_device_T(deck, "M1", SWEEP, os.path.join(outdir, "native"))
    freqs_b, t_brute = run_brute_T(deck, "m1", "gg", "d", "0", "0", op, "gm", SWEEP, os.path.join(outdir, "brute"))
    if freqs_n != freqs_b:
        # same sweep spec on the same build -- should match exactly; tolerate float noise
        if len(freqs_n) != len(freqs_b) or any(abs(a - b) / b > 1e-9 for a, b in zip(freqs_n, freqs_b)):
            raise stb.CouldNotRun(f"{stem}: frequency grids differ")
    err = compare_T(f"{stem}: native device-mode T vs independent brute-force (check 1a)", t_native, t_brute)
    # sanity: T is not trivially ~0 everywhere (the local Cgd loop must show up)
    record(f"{stem}: local Cgd loop is not negligible (|T| grows with frequency)",
           abs(t_native[-1]) > 10 * abs(t_native[0]) + 1e-9,
           f"|T[0]|={abs(t_native[0]):.3e} |T[-1]|={abs(t_native[-1]):.3e}")
    return deck, t_native, t_brute, outdir, err


def check_r1d():
    """check 1b: refs/r1d_mos_follower.cir, M2 (a real nonlinear 5-transistor
    OTA follower; M2's source is the diff-pair tail node, not ground, so
    gmbs is generically nonzero here -- a real-circuit cross-check beyond
    the synthetic fixtures)."""
    stem = "r1d_M2"
    outdir = os.path.join(DATA_ROOT, stem)
    deck_path = os.path.join(ROOT, "refs", "r1d_mos_follower.cir")
    deck = stb.prepare_deck(deck_path)
    op = op_extract(deck, "m2", os.path.join(outdir, "op"))
    record(f"{stem} op point sane (gm>0, gds>0)", op["gm"] > 0 and op["gds"] > 0,
           f"gm={op['gm']:.6g} gds={op['gds']:.6g} gmbs={op['gmbs']:.6g}")
    sweep = "dec 10 1 1e10"
    freqs_n, t_native = run_native_device_T(deck, "M2", sweep, os.path.join(outdir, "native"))
    freqs_b, t_brute = run_brute_T(deck, "m2", "inn", "out", "tail", "0", op, "gm", sweep, os.path.join(outdir, "brute"),
                                    neutralize=["@v.xp.vi[acmag]"])
    # Looser bound than check 1a's 1e-6 (spec.md stb_device's own number,
    # met exactly by the single-transistor cs_cgd fixture below): R1d is a
    # 5-transistor circuit, and replacing M2 with a linear equivalent for
    # the brute-force .ac necessarily perturbs the REST of the circuit's
    # DC operating point slightly (M1/M3/M4/M5's bias partly depends on
    # M2's real, nonlinear DC drain current -- an explicit DC current
    # source matching M2's real Id, see build_brute_lines(), corrects
    # most but not all of this: the added Rds_m2 resistor itself still
    # carries a small DC leakage current at the node's real bias voltage,
    # found while validating this card, on the order of the residual
    # error reported here). 1e-6-level agreement is check 1a's job; this
    # check's job is a real-circuit cross-check, which 1-2% relative
    # agreement on a full nonlinear OTA gives honestly.
    err = compare_T(f"{stem}: native device-mode T vs independent brute-force (check 1b, refs/r1d_mos_follower.cir M2)",
                     t_native, t_brute, tol=0.02)
    return deck, t_native, t_brute, outdir, err


def check_applicability_case1():
    """check 2: Tian applicability case 1 (research/cadence-stb-pstb.md sec
    2.6) -- a shunt-feedback CS stage whose only loop is the Rf wire, no
    local Cgd loop. Loop-mode T at the Rf wire (probe=Rf:R+, the existing
    inst:terminal auto-splice, card #457) should agree with device-mode T
    on M1 (the sole gain device in the sole loop)."""
    stem = "shunt_fb"
    outdir = os.path.join(DATA_ROOT, stem)
    deck = stb.prepare_deck(os.path.join(FIXDIR, "stbdev_shunt_fb.cir"))

    r1 = os.path.join(outdir, "loop_native.txt")
    if os.path.exists(r1):
        os.remove(r1)
    # T (not loopgain) directly, same vector both modes expose (spec.md
    # stb_compat: "T is also available by name") -- avoids a `let` step
    # that was found (validating this card) to trip an unrelated
    # "Option 'singlescale' not possible" / stale-vector-length error
    # specific to combining `let` with the probe=inst:terminal splice.
    text = ("\n".join(deck) +
            f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
            f"stb {SWEEP} probe=Rf:R+\nwrdata {r1} stb1.t\nquit\n.endc\n.end\n")
    proc = run_ngspice(text, os.path.join(outdir, "loop"), "native_loop")
    check_errors(proc, "loop-mode Rf:R+ run")
    freqs_l, (t_loop,) = stb.read_wrdata(r1, 1)

    freqs_d, t_dev = run_native_device_T(deck, "M1", SWEEP, os.path.join(outdir, "dev"))

    if len(freqs_l) != len(freqs_d):
        raise stb.CouldNotRun(f"{stem}: loop-mode has {len(freqs_l)} points, device-mode has {len(freqs_d)}")
    err = compare_T(f"{stem}: loop-mode T (Rf wire) vs device-mode T (M1), Tian case 1 agreement (check 2)",
                     t_loop, t_dev, tol=1e-6)
    return err


def check_red_proofs(cs_deck, cs_native, cs_brute, cs_outdir):
    # (a) sign flip must fail (same style as tests/test_stb_native.py's
    # check_red_proofs: a mutated/bogus comparison must NOT pass the bound,
    # proving the comparison is sign-sensitive rather than vacuously loose).
    expect_fail_T("red-proof: sign-flipped native T vs brute force", [-z for z in cs_native], cs_brute)

    # (b) probing a passive element must be COULD-NOT-RUN, never a number.
    deck = cs_deck
    outdir = os.path.join(cs_outdir, "redproof_passive")
    try:
        run_native_device_T(deck, "Rd", SWEEP, outdir)
        record("red-proof: probe=<passive element (Rd)> must not produce a number", False,
               "native stb ran to completion on a resistor probe -- should have errored")
    except stb.CouldNotRun as e:
        record("red-proof: probe=<passive element (Rd)> must not produce a number", True, f"COULD-NOT-RUN as expected: {e}")

    # (c) probing a nonexistent device must be COULD-NOT-RUN.
    outdir2 = os.path.join(cs_outdir, "redproof_nonexistent")
    try:
        run_native_device_T(deck, "Mzzz99", SWEEP, outdir2)
        record("red-proof: probe=<nonexistent device> must not produce a number", False,
               "native stb ran to completion on a nonexistent probe -- should have errored")
    except stb.CouldNotRun as e:
        record("red-proof: probe=<nonexistent device> must not produce a number", True, f"COULD-NOT-RUN as expected: {e}")


def check_body_effect():
    """check 3 / red-proof: gm-only vs gm+gmbs is not a negligible choice.

    On a source-degenerated stage (Vbs != 0, gmbs materially nonzero -- see
    the fixture's own op-point note), this computes TWO independent,
    individually well-defined single-source Bode return ratios by the same
    brute-force method as check_cs_cgd/check_r1d: T(gm) (gm's control
    disconnected and driven, gmbs left wired normally -- this is what
    native device mode computes, gm-only) and T(gmbs) (the mirror image:
    gmbs's control disconnected and driven, gm left wired normally).
    Native must match T(gm) to the same 1e-6 bound as check 1, AND
    |T(gmbs)| must be a non-negligible fraction of |T(gm)| at at least one
    frequency -- i.e. gmbs's own contribution to this circuit's loop is not
    itself negligible, which is exactly the fact that makes "gm-only vs
    gm+gmbs" a real documentation decision and not a rounding footnote.

    NOT attempted here (see the card #459 report's "what was not done"): a
    single Bode measurement of the JOINTLY-nulled combined source
    (nulling gm and gmbs together as one rank-1 term, matching this file's
    stb_mos1_probe() doc comment's "gm+gmbs together is a rank-1 update
    too" remark) -- gm and gmbs have different control ports (g,s) vs
    (b,s), so there is no single physical 2-terminal "drive with 1V" port
    for their sum without an extra, less standard construction; T(gm) and
    T(gmbs) individually are each a completely standard, independently
    meaningful Bode return ratio, and are enough to show the choice is
    consequential.
    """
    stem = "bodyeffect"
    outdir = os.path.join(DATA_ROOT, stem)
    deck = stb.prepare_deck(os.path.join(FIXDIR, "stbdev_bodyeffect.cir"))
    op = op_extract(deck, "m1", os.path.join(outdir, "op"))
    record(f"{stem} op point has nonzero body effect (gmbs > 0)", op["gmbs"] > 0,
           f"gm={op['gm']:.6g} gmbs={op['gmbs']:.6g} gds={op['gds']:.6g}")

    freqs_n, t_native = run_native_device_T(deck, "M1", SWEEP, os.path.join(outdir, "native"))
    freqs_gm, t_gm = run_brute_T(deck, "m1", "gg", "d", "s", "0", op, "gm", SWEEP, os.path.join(outdir, "brute_gm"))
    compare_T(f"{stem}: native (gm-only) vs brute-force T(gm) (check 3)", t_native, t_gm)

    try:
        freqs_gmbs, t_gmbs = run_brute_T(deck, "m1", "gg", "d", "s", "0", op, "gmbs", SWEEP, os.path.join(outdir, "brute_gmbs"))
        ratios = [abs(b) / abs(a) for a, b in zip(t_gm, t_gmbs) if abs(a) > 1e-30]
        worst = max(ratios) if ratios else 0.0
        record(f"{stem}: |T(gmbs)|/|T(gm)| is measurable (gm-only vs gm+gmbs is a real choice)",
               worst > 0.05,
               f"max |T(gmbs)|/|T(gm)| over the sweep = {worst:.4f} "
               f"(gm={op['gm']:.4e}, gmbs={op['gmbs']:.4e}, gmbs/gm={op['gmbs'] / op['gm']:.4f})")
        expect_fail_T(f"{stem}: native (gm-only) must not equal T(gmbs)-only", t_native, t_gmbs, tol=1e-3)
    except stb.CouldNotRun as e:
        record(f"{stem}: |T(gmbs)|/|T(gm)| is measurable (gm-only vs gm+gmbs is a real choice)", None, str(e))


def check_regression():
    env = dict(os.environ)
    env["NGSPICE"] = NGSPICE
    for name in ("test_stb_native2.py", "test_stb_native.py", "test_stb.py"):
        path = os.path.join(HERE, name)
        if not os.path.exists(path):
            record(f"{name} regression (skipped, file not found)", "skip", "not present in this checkout")
            continue
        proc = subprocess.run([sys.executable, path], capture_output=True, text=True, env=env, timeout=1800)
        tail = "\n".join(proc.stdout.splitlines()[-3:])
        record(f"{name} regression with NGSPICE=<stbdev build>", proc.returncode == 0, f"exit {proc.returncode}: {tail}")


def main():
    try:
        deck, t_native, t_brute, outdir, _ = check_cs_cgd()
        check_red_proofs(deck, t_native, t_brute, outdir)
    except stb.CouldNotRun as e:
        record("check 1a (cs_cgd)", None, str(e))

    try:
        check_r1d()
    except stb.CouldNotRun as e:
        record("check 1b (r1d M2)", None, str(e))

    try:
        check_applicability_case1()
    except stb.CouldNotRun as e:
        record("check 2 (Tian applicability case 1)", None, str(e))

    try:
        check_body_effect()
    except stb.CouldNotRun as e:
        record("check 3 (body effect, gm-only vs gm+gmbs)", None, str(e))

    check_regression()

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    npass = sum(s == "pass" for _, s, _ in results)
    print(f"\n{npass}/{len(results)} pass ({sum(s == 'skip' for _, s, _ in results)} skipped, not counted)")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
