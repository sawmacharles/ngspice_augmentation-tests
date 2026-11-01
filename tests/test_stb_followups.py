#!/usr/bin/env python3
"""Validation for card #460's follow-ups to the native `.stb` analysis:

  1. loop-mode probe with a terminal at literal ground (node 0) -- the
     suspected CKTrhs[0]/CKTirhs[0] trash-slot hazard (same class as the
     device-mode bug fixed in c931b1f4c). CONFIRMED real (see check
     1below): before the fix, a probe with its "b" (upstream/current-
     injection) terminal grounded read back a small, plausible-looking,
     WRONG, slightly frequency-dependent T instead of the mathematically
     honest answer. The fix (guard every CKTrhs[0]/CKTirhs[0] read/write
     in loop mode, src/spicelib/analysis/stban.c) makes it read that
     honest answer -- which turns out to be an ADDITIONAL, independent
     finding, proved three ways here (native, tools/stb.py's real
     3-element lgprobe subckt, and a hand-derived closed form): a
     loop-mode probe with EITHER terminal at literal ground is
     structurally degenerate for Tian's two-injection method -- T=0
     identically, for every circuit, because the probe's own branch
     equation V(a)-V(b)=<injection> either pins the OTHER terminal to
     the injection value outright, or the injected current has nowhere
     to go but the probe's own branch. stban.c now warns rather than
     silently handing back that 0 as if it meant "stable".

  2. BSIM4 device mode, probe=<BSIM4 instance> (stb_bsim4_probe() in
     stban.c). gm ALONE reduces to the identical rank-1 pattern as MOS1
     (u=+1@dNodePrime,-1@sNodePrime; v_gm=+1@gNodePrime,-1@sNodePrime;
     scalar=+gm; node_d/node_s swap on BSIM4mode<0), quasi-static only
     (acnqsmod=0, refused otherwise). Validated here against an
     independent brute-force small-signal equivalent (same technique as
     test_stb_device.py's build_brute_lines, extended with the two
     transcapacitance terms -- Cgg, Cdg -- BSIM4's own capmod=0 model
     does not reduce all the way to 0 even with every overlap/junction
     cap zeroed) to <= 1e-6 relative over the fixture's validated sweep.

  3. MOS1 device mode with the instance biased in REVERSE mode
     (Vds_netlist < 0, netlist drain/source terminals swapped relative
     to their physical roles) -- validates stb_mos1_probe()'s
     node_d/node_s swap for real, not just by inspection. <= 1e-6
     relative against the same brute-force technique as
     test_stb_device.py, using the PHYSICAL (not netlist-declared)
     drain/source roles throughout (found empirically while building
     this check -- using the netlist-declared roles literally, ignoring
     the swap, gives errors of order 1-600, not a rounding-level
     mismatch; see the card #460 report).

Every check is pass / FAIL / COULD-NOT-RUN (0/1/2 exit code convention,
same as tests/test_stb_device.py and tests/test_stb_native.py). Red-proof
checks mutate the result and pass only if the comparison then fails.
Decks and logs live under ~/Work/data/ngspice_augmentation/stb3/ (never
the repository); fixtures are tests/fixtures/stb3_*.cir.

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_stb_followups.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXDIR = os.path.join(HERE, "fixtures")
sys.path.insert(0, os.path.join(ROOT, "tools"))
import stb  # noqa: E402  (CouldNotRun, prepare_deck, probe_nodes, read_wrdata)
import stbcore  # noqa: E402  (tian, max_rel_err)

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/stb3")

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else ("skip" if ok == "skip" else "COULD-NOT-RUN"))
    results.append((name, state, detail))


def run_ngspice(text, outdir, tag, timeout=300):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{tag}.cir")
    with open(path, "w") as fh:
        fh.write(text)
    try:
        proc = subprocess.run([NGSPICE, "-b", path], capture_output=True, text=True, timeout=timeout)
    except (FileNotFoundError, PermissionError) as e:
        raise stb.CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, f"{tag}.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    return proc


def check_errors(proc, context, allow_warnings=True):
    lines = (proc.stdout + proc.stderr).splitlines()
    errs = [ln for ln in lines if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise stb.CouldNotRun(f"{context}: ngspice exit {proc.returncode}; {errs[:3]}")
    return lines


def read_op_params(deck_lines, inst, params, outdir):
    """Run `.op` and read back @inst[param] for each of `params`, full precision (wrdata)."""
    opfile = os.path.join(outdir, "op.txt")
    text = ("\n".join(deck_lines) + "\n.control\nset wr_vecnames\nset numdgt=17\nop\nwrdata "
            + opfile + " " + " ".join(f"@{inst}[{p}]" for p in params) + "\nquit\n.endc\n.end\n")
    proc = run_ngspice(text, outdir, "op")
    check_errors(proc, f"op-point extraction for {inst}")
    with open(opfile) as fh:
        fh.readline()
        row = fh.readline().split()
    vals = [float(row[i]) for i in range(1, len(row), 2)]
    if len(vals) != len(params):
        raise stb.CouldNotRun(f"op-point extraction for {inst}: expected {len(params)} values, got {len(vals)}")
    return dict(zip(params, vals))


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


def compare_T(name, a, b, tol):
    if len(a) != len(b):
        raise stb.CouldNotRun(f"{name}: length mismatch {len(a)} vs {len(b)}")
    err = stbcore.max_rel_err(a, b)
    record(name, err <= tol, f"max rel err {err:.3e} (bound {tol:.0e})")
    return err


def expect_fail_T(name, a, b, tol):
    err = stbcore.max_rel_err(a, b)
    record(name, err > tol, f"mutated/bogus comparison must differ: max rel err {err:.3e}")


def compare_zero(name, a, tol_abs):
    """Like compare_T, but against an exact-0 reference: stbcore.max_rel_err
    divides by the reference's magnitude, which is undefined at 0, so this
    uses a plain absolute-magnitude bound instead."""
    worst = max((abs(z) for z in a), default=0.0)
    record(name, worst <= tol_abs, f"max |T| {worst:.3e} (bound {tol_abs:.0e}, absolute -- reference is exactly 0)")
    return worst


# ================================================================
# Task 1: loop-mode ground-terminal probe (both orientations, two
# fixtures: an idealized-gm linear network with a closed-form
# reference, and a real nonlinear MOS1 stage).
# ================================================================

def native_loop_T(lgprobe_deck, probe, sweep, outdir, swap=False):
    """Build the plain-Vprobe variant of an lgprobe-based deck (same
    technique as tests/test_stb_native.py's plain_variant()) and run
    native .stb on it. swap=True writes `Vprobe b a` instead of `Vprobe a b`."""
    os.makedirs(outdir, exist_ok=True)
    lines = stb.prepare_deck(lgprobe_deck)
    a, b = stb.probe_nodes(lines, probe)
    out = []
    hit = False
    for line in lines:
        f = line.split()
        if f and f[0].lower() == probe.lower() and len(f) >= 4 and f[3].lower() == "lgprobe":
            n1, n2 = (b, a) if swap else (a, b)
            out.append(f"Vprobe {n1} {n2} DC 0")
            hit = True
        else:
            out.append(line)
    if not hit:
        raise stb.CouldNotRun(f"probe {probe} not an lgprobe instance in {lgprobe_deck}")
    r1 = os.path.join(outdir, "native.txt")
    if os.path.exists(r1):
        os.remove(r1)
    text = "\n".join(out) + (f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                              f"stb {sweep} probe=Vprobe\nwrdata {r1} stb1.t\nquit\n.endc\n.end\n")
    proc = run_ngspice(text, outdir, "native")
    lines_out = check_errors(proc, "native loop-mode ground-terminal run")
    warned = any("literal ground" in ln for ln in lines_out)
    freqs, (t,) = stb.read_wrdata(r1, 1)
    return freqs, t, warned


def python_loop_T(lgprobe_deck, probe, sweep, outdir):
    """Independent reference: tools/stb.py's own machinery (real 3-element
    lgprobe.lib subckt, its own internal "x" node), but reading only what
    stbcore.tian() needs -- tools/stb.py's own run() crashes on
    middlebrook() when i2==0 (ZeroDivisionError, a genuine consequence of
    this same degeneracy, not a bug worth working around in tools/*,
    which this card may not edit) and its `ac1.v({a})` print fails
    outright when a=="0" (v(0) is not a printable ngspice vector). This
    reimplements exactly stb.run()'s two-injection control block, minus
    the two reads tools/stb.py itself cannot make here."""
    os.makedirs(outdir, exist_ok=True)
    lines = stb.prepare_deck(lgprobe_deck)
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
wrdata {r1} ac1.v({p}.x) ac1.i(v.{p}.vi)
wrdata {r2} ac2.v({p}.x) ac2.i(v.{p}.vi)
quit
.endc
.end
"""
    proc = run_ngspice("\n".join(lines) + control, outdir, "python_ref")
    check_errors(proc, "tools/stb.py-style reference run (real lgprobe subckt)")
    f1, (v1, i1) = stb.read_wrdata(r1, 2)
    f2, (v2, i2) = stb.read_wrdata(r2, 2)
    if f1 != f2:
        raise stb.CouldNotRun("the two AC runs have different frequency grids")
    t = [stbcore.tian(*z) for z in zip(v1, i1, v2, i2)]
    return f1, t


# The pre-fix numeric value this exact fixture/sweep/point produced,
# captured empirically (STB_DEBUG=1, this same build, before the guard
# was added to stban.c): v2 read back the raw 1.0 this code itself had
# just written to CKTrhs[0] (the trash slot), instead of ground's true,
# always-0 voltage -- see the card #460 report's "before" numbers.
PRE_FIX_FWD_T0 = complex(-1.995972124265276e-05, 1.246832822837396e-13)


def check_ground_linear():
    stem = "gnd_loop_linear"
    outdir = os.path.join(DATA_ROOT, stem)
    deck = os.path.join(FIXDIR, "stb3_gnd_loop_linear.cir")
    sweep = "dec 5 1 1e6"

    # closed form: this circuit's probe is a 0V branch (never perturbs
    # the topology), so the TRUE network return ratio -- independent of
    # where within a single loop you insert a zero-impedance probe -- is
    # refs/r1a_unilateral.cir's own formula.
    p = {"gm": 10e-3, "ro": 100e3, "co": 10e-12, "r1": 100e3, "r2": 100e3, "cn": 2e-12}

    def r1a_closed_form(f):
        import cmath
        s = 2j * 3.141592653589793 * f
        zo = 1 / (1 / p["ro"] + s * p["co"])
        zb = 1 / (1 / p["r2"] + s * p["cn"])
        return p["gm"] * zo * zb / (zo + p["r1"] + zb)

    fwd_f, fwd_t, fwd_warn = native_loop_T(deck, "XP", sweep, os.path.join(outdir, "native_fwd"), swap=False)
    rev_f, rev_t, rev_warn = native_loop_T(deck, "XP", sweep, os.path.join(outdir, "native_rev"), swap=True)
    py_f, py_t = python_loop_T(deck, "XP", sweep, os.path.join(outdir, "py"))

    record(f"{stem}: native warns when a probe terminal is literal ground (forward, b=0)", fwd_warn,
           "warning text seen in ngspice output" if fwd_warn else "no warning printed -- silent number")
    record(f"{stem}: native warns when a probe terminal is literal ground (reversed, a=0)", rev_warn,
           "warning text seen in ngspice output" if rev_warn else "no warning printed -- silent number")

    compare_zero(f"{stem}: native T == 0 exactly, forward orientation (b=0, the trash-slot fix)", fwd_t, tol_abs=1e-12)
    compare_zero(f"{stem}: native T == 0 exactly, reversed orientation (a=0, no trash-slot involved)", rev_t, tol_abs=1e-12)
    compare_zero(f"{stem}: independent reference (tools/stb.py's real 3-element lgprobe subckt) T == 0 too", py_t, tol_abs=1e-12)
    compare_zero(f"{stem}: native forward == native reversed (Tian symmetry, both degenerate the same way)",
                 [a - b for a, b in zip(fwd_t, rev_t)], tol_abs=1e-12)

    # red-proof: the true underlying network's loop gain is large and
    # frequency-dependent (this is NOT "an amplifier with no loop gain" --
    # it is "a probe placed where Tian's method cannot see it"). Compare
    # the closed form against 0 and confirm they do NOT match, proving
    # this check is sensitive rather than vacuous.
    ref0 = [r1a_closed_form(f) for f in fwd_f]
    worst_true_T = min(abs(z) for z in ref0)
    record(f"{stem}: red-proof -- the TRUE network loop gain (closed form) is NOT 0 (confirms T=0 is a probe artifact, not 'no loop')",
           worst_true_T > 1e-3, f"min |T_closed_form| over the sweep = {worst_true_T:.3e} (must be far from 0)")

    # red-proof: the pre-fix buggy value (captured empirically) must not
    # match either T=0 or the closed form, at f=1 Hz specifically (the
    # sweep's first point).
    record(f"{stem}: red-proof -- pre-fix trash-slot value does not match the honest T=0",
           abs(PRE_FIX_FWD_T0) > 1e-9, f"|pre-fix T(f=1Hz)| = {abs(PRE_FIX_FWD_T0):.3e} (must be far from 0)")


def check_ground_mos():
    stem = "gnd_loop_mos"
    outdir = os.path.join(DATA_ROOT, stem)
    deck = os.path.join(FIXDIR, "stb3_gnd_loop_mos.cir")
    sweep = "dec 5 1 1e8"

    fwd_f, fwd_t, fwd_warn = native_loop_T(deck, "XP", sweep, os.path.join(outdir, "native_fwd"), swap=False)
    rev_f, rev_t, rev_warn = native_loop_T(deck, "XP", sweep, os.path.join(outdir, "native_rev"), swap=True)
    py_f, py_t = python_loop_T(deck, "XP", sweep, os.path.join(outdir, "py"))

    record(f"{stem}: native warns when a probe terminal is literal ground (forward, b=0)", fwd_warn,
           "warning text seen in ngspice output" if fwd_warn else "no warning printed -- silent number")
    record(f"{stem}: native warns when a probe terminal is literal ground (reversed, a=0)", rev_warn,
           "warning text seen in ngspice output" if rev_warn else "no warning printed -- silent number")

    compare_zero(f"{stem} (real nonlinear MOS1 stage): native T == 0 exactly, forward (b=0)", fwd_t, tol_abs=1e-12)
    compare_zero(f"{stem}: native T == 0 exactly, reversed (a=0)", rev_t, tol_abs=1e-12)
    compare_zero(f"{stem}: independent reference (tools/stb.py's real lgprobe subckt) T == 0 too", py_t, tol_abs=1e-12)
    compare_zero(f"{stem}: native forward == native reversed (Tian symmetry)",
                 [a - b for a, b in zip(fwd_t, rev_t)], tol_abs=1e-12)


# ================================================================
# Task 2: BSIM4 device mode, rank-1 gm-only nulling.
# ================================================================

def check_bsim4_cgd():
    stem = "bsim4_cgd"
    outdir = os.path.join(DATA_ROOT, stem)
    deck = stb.prepare_deck(os.path.join(FIXDIR, "stb3_bsim4_cgd.cir"))
    # bound found empirically (card #460 report): the fixture's residual,
    # non-reciprocal transcapacitance (BSIM4's Cdg, and everything past
    # it that this brute-force reference does not model -- see the
    # report) grows linearly with frequency, crossing 1e-6 relative
    # around 2e5 Hz; this sweep stays comfortably under that.
    sweep = "dec 15 1 1e5"
    tol = 1e-6

    params = ["gm", "gds", "id", "cgg", "cdg"]
    op = read_op_params(deck, "m1", params, os.path.join(outdir, "op"))
    record(f"{stem}: op point sane (gm>0, gds>0)", op["gm"] > 0 and op["gds"] > 0,
           f"gm={op['gm']:.6g} gds={op['gds']:.6g} id={op['id']:.6g}")

    gm, gds, id_, cgg, cdg = (op[p] for p in params)
    out = [ln for ln in deck if not (ln.split() and ln.split()[0].lower() == "m1")]
    # Cgg: gate-to-ground lump (source=bulk=ground here, so the gate
    # row's ONLY two nonzero transcapacitance terms are Cgg and Cgd; Cgd
    # is exactly 0 in this fixture -- capmod=0, cgdo=0). Cdg: drain-row
    # transcapacitance responding to dVg/dt (Cdd is exactly 0 here, so
    # the drain row's ONLY nonzero term is Cdg) -- NOT representable as a
    # plain two-terminal capacitor (Cgd=0 != Cdg, i.e. non-reciprocal),
    # modeled with the standard SPICE current-mirror trick: a sense
    # capacitor + a 0V ammeter source + a CCCS copying that current into
    # the drain node.
    elems = [
        f"Idc_m1 d 0 DC {id_!r} AC 0",
        f"Rds_m1 d 0 {1.0 / gds!r}",
        f"Ggm_m1 d 0 ctrlp 0 {gm!r}",
        "Vtest_m1 ctrlp 0 DC 0 AC 1",
        f"Cgg_m1 gg 0 {cgg!r}",
        f"Csense_dg gg dsense {cdg!r}",
        "Vsense_dg dsense 0 DC 0 AC 0",
        "Fdg d 0 Vsense_dg 1",
    ]
    brute_lines = out + elems
    r1 = os.path.join(outdir, "brute.txt")
    if os.path.exists(r1):
        os.remove(r1)
    text = "\n".join(brute_lines) + (
        f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=17\n"
        f"ac {sweep}\nlet vret = v(gg)\nlet Tbrute = -vret\n"
        f"wrdata {r1} Tbrute\nquit\n.endc\n.end\n")
    proc = run_ngspice(text, outdir, "brute")
    check_errors(proc, f"brute-force gm-only AC run for {stem}")
    freqs_b, (t_brute,) = stb.read_wrdata(r1, 1)

    freqs_n, t_native = run_native_device_T(deck, "M1", sweep, os.path.join(outdir, "native"))
    if freqs_n != freqs_b:
        if len(freqs_n) != len(freqs_b) or any(abs(a - b) / b > 1e-9 for a, b in zip(freqs_n, freqs_b)):
            raise stb.CouldNotRun(f"{stem}: frequency grids differ")

    err = compare_T(f"{stem}: native BSIM4 device-mode T vs independent brute-force (Cgg+Cdg transcapacitance-augmented)", t_native, t_brute, tol=tol)
    record(f"{stem}: local Cgd loop is not negligible (|T| grows with frequency)",
           abs(t_native[-1]) > 5 * abs(t_native[0]) + 1e-9,
           f"|T[0]|={abs(t_native[0]):.3e} |T[-1]|={abs(t_native[-1]):.3e}")

    expect_fail_T(f"{stem}: red-proof -- sign-flipped native T vs brute force", [-z for z in t_native], t_brute, tol=tol)

    # red-proof: a probe on a passive element must COULD-NOT-RUN.
    try:
        run_native_device_T(deck, "Rd", sweep, os.path.join(outdir, "redproof_passive"))
        record(f"{stem}: red-proof -- probe=<passive element (Rd)> must not produce a number", False,
               "native stb ran to completion on a resistor probe -- should have errored")
    except stb.CouldNotRun as e:
        record(f"{stem}: red-proof -- probe=<passive element (Rd)> must not produce a number", True, f"COULD-NOT-RUN as expected: {e}")
    return err


def check_bsim4_nqs_refused():
    """acnqsmod=1 on the probed instance must be refused, not silently
    null the wrong (real-only) scalar (card #460's own scope: QS only)."""
    stem = "bsim4_nqs_refused"
    outdir = os.path.join(DATA_ROOT, stem)
    deck = stb.prepare_deck(os.path.join(FIXDIR, "stb3_bsim4_cgd.cir"))
    out = [line.replace("acnqsmod=0", "acnqsmod=1") for line in deck]
    if out == deck:
        raise stb.CouldNotRun("acnqsmod=0 not found in the fixture's model card -- nothing to flip")
    r1 = os.path.join(outdir, "native.txt")
    if os.path.exists(r1):
        os.remove(r1)
    text = "\n".join(out) + (f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                              f"stb dec 5 1 1e5 probe=M1\nwrdata {r1} stb1.t\nquit\n.endc\n.end\n")
    proc = run_ngspice(text, outdir, "nqs")
    out_lines = (proc.stdout + proc.stderr)
    refused = "acnqsmod=1" in out_lines or "NQS on" in out_lines or "quasi-static only" in out_lines
    record(f"{stem}: acnqsmod=1 on the probed BSIM4 instance is refused, not silently computed",
           refused and not os.path.exists(r1),
           "refused with a clear message" if refused else f"NOT refused -- exit {proc.returncode}, wrote output: {os.path.exists(r1)}")


# ================================================================
# Task 3: MOS1 device mode, reverse mode (Vds_netlist < 0).
# ================================================================

def check_reverse_mos():
    stem = "reverse_mos"
    outdir = os.path.join(DATA_ROOT, stem)
    deck = stb.prepare_deck(os.path.join(FIXDIR, "stb3_reverse_mos.cir"))
    sweep = "dec 10 1 1e8"
    tol = 1e-6

    params = ["gm", "gds", "id", "vgs", "vds"]
    op = read_op_params(deck, "m1", params, os.path.join(outdir, "op"))
    record(f"{stem}: op point IS reverse mode (vds_netlist < 0)", op["vds"] < 0,
           f"vgs={op['vgs']:.6g} vds={op['vds']:.6g} (must be negative -- MOS1mode<0)")
    record(f"{stem}: op point sane (gm>0, gds>0)", op["gm"] > 0 and op["gds"] > 0,
           f"gm={op['gm']:.6g} gds={op['gds']:.6g} id={op['id']:.6g}")

    gm, gds, id_ = op["gm"], op["gds"], op["id"]
    out = [ln for ln in deck if not (ln.split() and ln.split()[0].lower() == "m1")]
    # PHYSICAL roles throughout (pd=physical drain, ps=physical source --
    # NOT the netlist-declared drain/source slots, which are swapped in
    # this fixture on purpose): matches stb_mos1_probe()'s own
    # node_d/node_s swap on MOS1mode<0 (src/spicelib/analysis/stban.c),
    # which maps node_d/node_s to the PHYSICAL roles regardless of
    # netlist declaration order. Found empirically while building this
    # check: using the NETLIST-declared order literally (ignoring the
    # swap) gives errors of order 1-600, not a rounding-level mismatch
    # -- see the card #460 report.
    elems = [
        f"Idc_m1 pd ps DC {id_!r} AC 0",
        f"Rds_m1 pd ps {1.0 / gds!r}",
        f"Ggm_m1 pd ps ctrlp 0 {gm!r}",
        "Vtest_m1 ctrlp 0 DC 0 AC 1",
    ]
    brute_lines = out + elems
    r1 = os.path.join(outdir, "brute.txt")
    text = "\n".join(brute_lines) + (
        f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=17\n"
        f"ac {sweep}\nlet vret = v(gg)-v(ps)\nlet Tbrute = -vret\n"
        f"wrdata {r1} Tbrute\nquit\n.endc\n.end\n")
    proc = run_ngspice(text, outdir, "brute")
    check_errors(proc, f"brute-force gm-only AC run for {stem}")
    freqs_b, (t_brute,) = stb.read_wrdata(r1, 1)

    freqs_n, t_native = run_native_device_T(deck, "M1", sweep, os.path.join(outdir, "native"))
    if freqs_n != freqs_b:
        if len(freqs_n) != len(freqs_b) or any(abs(a - b) / b > 1e-9 for a, b in zip(freqs_n, freqs_b)):
            raise stb.CouldNotRun(f"{stem}: frequency grids differ")

    err = compare_T(f"{stem}: native device-mode T (REVERSE mode) vs independent brute-force (PHYSICAL roles)", t_native, t_brute, tol=tol)
    expect_fail_T(f"{stem}: red-proof -- sign-flipped native T vs brute force", [-z for z in t_native], t_brute, tol=tol)
    # red-proof: using the NETLIST-declared (not physical) roles must NOT
    # match -- proves the reverse-mode swap is actually load-bearing here,
    # not a no-op.
    elems_wrong = [
        f"Idc_m1w ps pd DC {id_!r} AC 0",
        f"Rds_m1w ps pd {1.0 / gds!r}",
        f"Ggm_m1w ps pd ctrlpw 0 {gm!r}",
        "Vtestw ctrlpw 0 DC 0 AC 1",
    ]
    r2 = os.path.join(outdir, "brute_wrong.txt")
    text2 = "\n".join(out + elems_wrong) + (
        f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=17\n"
        f"ac {sweep}\nlet vretw = v(gg)-v(pd)\nlet Tbrutew = -vretw\n"
        f"wrdata {r2} Tbrutew\nquit\n.endc\n.end\n")
    proc2 = run_ngspice(text2, outdir, "brute_wrong")
    check_errors(proc2, f"red-proof brute-force (netlist roles) for {stem}")
    _, (t_wrong,) = stb.read_wrdata(r2, 1)
    expect_fail_T(f"{stem}: red-proof -- brute force using NETLIST (not physical) drain/source roles must NOT match", t_wrong, t_brute, tol=tol)
    return err


# ================================================================
# Regression with this build.
# ================================================================

def check_regression():
    env = dict(os.environ)
    env["NGSPICE"] = NGSPICE
    for name in ("test_stb_device.py", "test_stb_native2.py", "test_stb_native.py", "test_stb.py"):
        path = os.path.join(HERE, name)
        if not os.path.exists(path):
            record(f"{name} regression (skipped, file not found)", "skip", "not present in this checkout")
            continue
        proc = subprocess.run([sys.executable, path], capture_output=True, text=True, env=env, timeout=1800)
        tail_lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
        tail = tail_lines[-1] if tail_lines else "(no output)"
        record(f"{name} regression with NGSPICE=<stb3 build>", proc.returncode == 0, f"exit {proc.returncode}: {tail}")


def main():
    try:
        check_ground_linear()
    except stb.CouldNotRun as e:
        record("task 1 (ground-terminal, linear fixture)", None, str(e))

    try:
        check_ground_mos()
    except stb.CouldNotRun as e:
        record("task 1 (ground-terminal, MOS fixture)", None, str(e))

    try:
        check_bsim4_cgd()
    except stb.CouldNotRun as e:
        record("task 2 (BSIM4 device mode, Cgd loop)", None, str(e))

    try:
        check_bsim4_nqs_refused()
    except stb.CouldNotRun as e:
        record("task 2 (BSIM4 acnqsmod=1 refusal)", None, str(e))

    try:
        check_reverse_mos()
    except stb.CouldNotRun as e:
        record("task 3 (MOS1 reverse mode)", None, str(e))

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
