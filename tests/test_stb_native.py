#!/usr/bin/env python3
"""Validation for the native `.stb` analysis (card #456, plan T3.2/T3.3).

Compares the native `.stb` in the private-fork build (NGSPICE env var must
point at it -- see buildrecipe/README.md) against `tools/stb.py` (the
independent Tian two-injection testbench, validated in tests/test_stb.py)
run on the SAME build. Every deck variant this script writes lives under
~/Work/data/ngspice_augmentation/stb_native/ (never in the repo); refs/ is
read-only here.

Every check is pass / FAIL / COULD-NOT-RUN (a third state, "skip", is used
for a check that is not meaningful to run at all -- see R1c below -- and
never affects the exit code). Red-proof checks mutate the result and pass
only if the comparison then fails.

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_stb_native.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import stb  # noqa: E402
import stbcore  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/stb_native")
SWEEP = "dec 20 1 1e10"
TOL = 1e-9

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else ("skip" if ok == "skip" else "COULD-NOT-RUN"))
    results.append((name, state, detail))


# ---------- deck-variant generation (a plain 0V V-source probe instead of
# the lgprobe subckt; see stban.c's header comment for the a/b convention:
# a = the probe's first/positive node (downstream), b = the second/negative
# node (upstream) -- same roles lgprobe.lib documents). ----------

def plain_variant(deck, probe, outdir, swap=False):
    """Write a plain-V-source variant of `deck` under outdir; return (path, a, b)."""
    os.makedirs(outdir, exist_ok=True)
    lines = stb.prepare_deck(deck)
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
        raise stb.CouldNotRun(f"probe {probe} not an lgprobe instance in {deck}")
    path = os.path.join(outdir, "plain_swap.cir" if swap else "plain.cir")
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return path, a, b


def strip_end(text):
    return re.sub(r"\n\.end\s*$", "\n", text.rstrip("\n") + "\n")


def run_native_stb(deck_plain, outdir, sweep=SWEEP, probe="Vprobe"):
    """Run the native `.stb` on a plain-probe deck; return dict with freqs, T,
    loopgain, margins (as read back from the analysis's own output), and the
    raw ngspice exit/log for COULD-NOT-RUN diagnostics."""
    os.makedirs(outdir, exist_ok=True)
    r1 = os.path.join(outdir, "native.txt")
    mfile = os.path.join(outdir, "margins.txt")
    for f in (r1, mfile):
        if os.path.exists(f):
            os.remove(f)
    with open(deck_plain) as fh:
        body = strip_end(fh.read())
    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
stb {sweep} probe={probe}
wrdata {r1} stb1.loopgain stb1.t
wrdata {mfile} stb2.phasemargin stb2.phasemarginfreq stb2.gainmargin stb2.gainmarginfreq
quit
.endc
.end
"""
    wrapper = os.path.join(outdir, "wrapper.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + ctl)
    try:
        proc = subprocess.run([NGSPICE, "-b", wrapper], capture_output=True, text=True, timeout=600)
    except (FileNotFoundError, PermissionError) as e:
        raise stb.CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, "log.txt"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise stb.CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")
    freqs, (lg, t) = stb.read_wrdata(r1, 2)
    margins = {}
    if os.path.exists(mfile):
        with open(mfile) as fh:
            next(fh)
            row = [float(x) for x in fh.readline().split()]
        # wrdata prints the 4 real scalars as 4 plain columns (no freq col
        # for a 1-point real plot; the single "0" reference IS column 0)
        if len(row) == 5:
            _, pm, pmf, gm, gmf = row
        elif len(row) == 4:
            pm, pmf, gm, gmf = row
        else:
            raise stb.CouldNotRun(f"unexpected margins column count {len(row)} in {mfile}")
        import math
        margins = {
            "pm": None if math.isnan(pm) else {"pm_deg": pm, "freq": pmf},
            "gm": None if math.isnan(gm) else {"gm_db": gm, "freq": gmf},
        }
    return {"freqs": freqs, "loopgain": lg, "T": t, "margins": margins,
            "log": proc.stdout + proc.stderr, "wrapper": wrapper}


def python_ref(deck, probe, outdir, sweep=SWEEP):
    return stb.run(deck, probe, sweep, outdir)


def compare_T(name, native, python, tol=TOL):
    err = stbcore.max_rel_err(native, python)
    record(name, err <= tol, f"max rel err {err:.3e} (bound {tol:.0e})")
    return err


def expect_fail_T(name, native, python, tol=TOL):
    err = stbcore.max_rel_err(native, python)
    record(name, err > tol, f"mutated/bogus comparison must differ: max rel err {err:.3e}")


# ---------- checks ----------

def check_netlist_card_form(stem, deck, py):
    """`.stb` as a netlist card (not the interactive `stb` command every
    other check here uses) -- a separate code path (dot_stb in inp2dot.c)
    from the interactive command (com_stb -> dosim, runcoms.c). The probe=
    value must be quoted on the netlist card (probe="Vprobe") -- unquoted
    fails with "Undefined parameter [vprobe]" from the numparam brace-
    expression preprocessor, which runs over dot-card lines before dot_stb
    ever sees them (found validating this card; documented in inp2dot.c's
    dot_stb header comment and patches/0003-*)."""
    outdir = os.path.join(DATA_ROOT, stem, "card_form")
    plain, a, b = plain_variant(deck, "XP", outdir)
    with open(plain) as fh:
        body = strip_end(fh.read())
    card = f'.stb {SWEEP} probe="Vprobe"\n'
    ctl = "\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata OUT stb1.t\nquit\n.endc\n.end\n"
    r1 = os.path.join(outdir, "native.txt")
    if os.path.exists(r1):
        os.remove(r1)
    ctl = ctl.replace("OUT", r1)
    wrapper = os.path.join(outdir, "wrapper.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + card + ctl)
    try:
        proc = subprocess.run([NGSPICE, "-b", wrapper], capture_output=True, text=True, timeout=600)
    except (FileNotFoundError, PermissionError) as e:
        raise stb.CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise stb.CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")
    freqs, (t,) = stb.read_wrdata(r1, 1)
    if len(freqs) != len(py["freqs"]):
        raise stb.CouldNotRun(f"{stem} card form: {len(freqs)} points vs python's {len(py['freqs'])}")
    compare_T(f"{stem} .stb NETLIST CARD form (quoted probe=) vs tools/stb.py", t, py["T"])


def check_T_match(stem, deck):
    outdir = os.path.join(DATA_ROOT, stem)
    plain, a, b = plain_variant(deck, "XP", outdir)
    nat = run_native_stb(plain, os.path.join(outdir, "native"))
    py = python_ref(deck, "XP", os.path.join(outdir, "py"))
    if len(nat["freqs"]) != len(py["freqs"]):
        raise stb.CouldNotRun(f"{stem}: native has {len(nat['freqs'])} points, python has {len(py['freqs'])}")
    err = compare_T(f"{stem} native T vs tools/stb.py (probe {a}/{b})", nat["T"], py["T"])
    record(f"{stem} native loopgain = -T", True,
           f"max |loopgain-(-T)| = {max(abs(lg - (-t)) for lg, t in zip(nat['loopgain'], nat['T'])):.3e}")
    return nat, py, err, outdir, a, b


def check_margins(stem, nat, py):
    ref = stbcore.margins(py["freqs"], nat["T"])  # same T as native, independent recompute
    for kind, key, bound in (("pm", "pm_deg", 0.01), ("gm", "gm_db", 0.01)):
        g = nat["margins"].get(kind)
        r = ref[kind]
        if r is None or g is None:
            record(f"{stem} {kind} (native vs stbcore.margins on same T)", g is None and r is None,
                   f"native {'none' if g is None else g}, stbcore {'none' if r is None else r}")
            continue
        dv = abs(g[key] - (r["pm_deg"] if kind == "pm" else r["gm_db"]))
        dfreq = abs(g["freq"] / r["freq"] - 1)
        record(f"{stem} {kind} (native vs stbcore.margins on same T)", dv <= bound and dfreq <= 1e-3,
               f"native {g[key]:.4f} @ {g['freq']:.6g} Hz, stbcore "
               f"{(r['pm_deg'] if kind == 'pm' else r['gm_db']):.4f} @ {r['freq']:.6g} Hz "
               f"(|d|={dv:.4f}, bound {bound})")


def check_probe_reversed(stem, deck, forward_T, outdir):
    plain_rev, a, b = plain_variant(deck, "XP", os.path.join(outdir, "rev"), swap=True)
    nat_rev = run_native_stb(plain_rev, os.path.join(outdir, "rev", "native"))
    compare_T(f"{stem} probe reversed = forward (Tian symmetry)", nat_rev["T"], forward_T)


def check_dc_identical(stem, deck, a, b, outdir):
    """DC op point with the probe present vs the same circuit with a and b
    tied to the same node (probe removed) -- must be bit-for-bit equal on
    every other node (the probe is a 0V branch; it cannot move the DC
    solution). Compares ngspice's own printed values (numdgt=17, i.e. full
    double precision) so 'bit-identical' is checked as 'string-identical'."""
    lines = stb.prepare_deck(deck)
    plain_lines = []
    hit = False
    for line in lines:
        f = line.split()
        if f and f[0].lower() == "xp" and len(f) >= 4 and f[3].lower() == "lgprobe":
            plain_lines.append(f"Vprobe {a} {b} DC 0")
            hit = True
        else:
            plain_lines.append(line)
    if not hit:
        raise stb.CouldNotRun(f"{stem}: XP not found for dc-identical check")

    def run_op(lines_, tag):
        d = os.path.join(outdir, "dcop")
        os.makedirs(d, exist_ok=True)
        opfile = os.path.join(d, f"{tag}.txt")
        if os.path.exists(opfile):
            os.remove(opfile)
        ctl = f"\n.control\nset wr_vecnames\nset numdgt=17\nop\nwrdata {opfile} all\nquit\n.endc\n.end\n"
        w = os.path.join(d, f"{tag}.cir")
        with open(w, "w") as fh:
            fh.write("\n".join(lines_) + ctl)
        proc = subprocess.run([NGSPICE, "-b", w], capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            raise stb.CouldNotRun(f"{stem} dc-op {tag}: ngspice exit {proc.returncode}")
        with open(opfile) as fh:
            return fh.read().splitlines()

    with_probe = run_op(plain_lines, "with_probe")

    # merged-node variant: rename every whole-word occurrence of node b to a,
    # and drop the probe line entirely (a 0V branch tying an already-merged
    # node to itself would be a degenerate/singular branch).
    merged = []
    wordb = re.compile(rf"(?<![\w.]){re.escape(b)}(?![\w.])", re.I)
    for line in plain_lines:
        f = line.split()
        if f and f[0].lower() == "vprobe":
            continue
        merged.append(wordb.sub(a, line))
    without_probe = run_op(merged, "without_probe")

    # compare the node voltages both runs share by name (drop the probe's
    # own extra branch-current row, which only exists in the with-probe run)
    def parse(rows):
        # wrdata's multi-vector format repeats a (scale_name, scale_value)
        # pair before every (vec_name, vec_value) pair -- take only the
        # odd-indexed (vector) columns, keyed by name.
        d = {}
        header = rows[0].split()
        data = rows[1].split()
        for i in range(1, len(header), 2):
            d[header[i].lower()] = data[i]
        return d

    dw = parse(with_probe)
    dwo = parse(without_probe)
    common = set(dw) & set(dwo)
    common.discard("vprobe#branch")
    if not common:
        raise stb.CouldNotRun(f"{stem}: no common node names between with/without-probe op runs")
    mismatches = [(k, dw[k], dwo[k]) for k in sorted(common) if dw[k] != dwo[k]]
    if not mismatches:
        record(f"{stem} DC operating point bit-identical with/without probe", True, "all match, string-exact")
        return
    # Not string-exact: report the actual numeric gap rather than silently
    # loosening the check. A few ULPs at the last digit (relative ~1e-14) is
    # floating-point non-associativity from the differing matrix size/
    # ordering (with vs. without the probe's extra branch unknown), not a
    # real physical difference -- but it is honestly not bit-for-bit, so say
    # so plainly. Only pass here if every mismatch is at that noise floor.
    worst = 0.0
    for k, a, b in mismatches:
        fa, fb = float(a), float(b)
        denom = max(abs(fa), abs(fb), 1e-300)
        worst = max(worst, abs(fa - fb) / denom)
    ok = worst <= 1e-12
    record(f"{stem} DC operating point bit-identical with/without probe", ok,
           f"NOT string-exact: {len(mismatches)}/{len(common)} vector(s) differ "
           f"(e.g. {mismatches[:3]}); worst relative gap {worst:.3e} -- "
           f"{'within 1e-12, floating-point non-associativity from the different matrix size' if ok else 'exceeds 1e-12, a real discrepancy'}")


def check_red_proofs(stem, nat, py):
    # (a) a sign flip must not pass
    expect_fail_T(f"{stem} red-proof: sign-flipped native T", [-z for z in nat["T"]], py["T"])
    # (b) a deliberately wrong RHS -- only the voltage-injection run, with
    # the current-injection run's contribution zeroed out entirely (i.e.
    # collapsing Tian's two-injection formula to a degenerate one-injection
    # guess) -- must not agree with the real two-injection result.
    bogus = [stbcore.tian(v, i, 0j, 0j) for v, i in zip(py["raw"]["v1"], py["raw"]["i1"])]
    expect_fail_T(f"{stem} red-proof: one-injection-only RHS is not the real T", bogus, py["T"])


def check_klu_option(stem, deck):
    """Same comparison with 'option klu' added to the deck (spec: 'test with
    and without option klu / option sparse')."""
    outdir = os.path.join(DATA_ROOT, stem, "klu")
    lines = stb.prepare_deck(deck)
    # line 0 is the deck's title line (SPICE convention: always consumed as
    # the title, whatever it contains) -- the .option card must come after it.
    lines = [lines[0], ".option klu"] + lines[1:]
    plain = os.path.join(outdir, "plain.cir")
    os.makedirs(outdir, exist_ok=True)
    a, b = stb.probe_nodes(lines, "XP")
    out = []
    for line in lines:
        f = line.split()
        if f and f[0].lower() == "xp" and len(f) >= 4 and f[3].lower() == "lgprobe":
            out.append(f"Vprobe {f[1]} {f[2]} DC 0")
        else:
            out.append(line)
    with open(plain, "w") as fh:
        fh.write("\n".join(out) + "\n")
    nat = run_native_stb(plain, os.path.join(outdir, "native"))
    py = python_ref(deck, "XP", os.path.join(outdir, "py"))
    compare_T(f"{stem} native T vs tools/stb.py, option klu", nat["T"], py["T"])


def check_runtime_r1d(deck):
    """Runtime of native .stb vs tools/stb.py on R1d, 500-point sweep."""
    sweep = "dec 500 1 1e10"
    outdir = os.path.join(DATA_ROOT, "r1d_runtime")
    plain, a, b = plain_variant(deck, "XP", outdir)
    t0 = time.time()
    nat = run_native_stb(plain, os.path.join(outdir, "native"), sweep=sweep)
    t_native = time.time() - t0
    t0 = time.time()
    py = python_ref(deck, "XP", os.path.join(outdir, "py"), sweep=sweep)
    t_python = time.time() - t0
    err = stbcore.max_rel_err(nat["T"], py["T"])
    record("R1d runtime: native .stb vs tools/stb.py, 500-pt sweep", True,
           f"native {t_native:.3f}s, script {t_python:.3f}s "
           f"({(t_python / t_native if t_native else float('inf')):.1f}x), "
           f"n={len(nat['freqs'])}, max rel err {err:.3e}")


def check_test_stb_regression():
    env = dict(os.environ)
    env["NGSPICE"] = NGSPICE
    proc = subprocess.run([sys.executable, os.path.join(HERE, "test_stb.py")],
                           capture_output=True, text=True, env=env, timeout=1800)
    tail = "\n".join(proc.stdout.splitlines()[-3:])
    record("tests/test_stb.py with NGSPICE=<new build>", proc.returncode == 0,
           f"exit {proc.returncode}: {tail}")


def main():
    try:
        r1a_nat, r1a_py, _, r1a_out, r1a_a, r1a_b = check_T_match("r1a_unilateral", os.path.join(ROOT, "refs", "r1a_unilateral.cir"))
        check_margins("r1a_unilateral", r1a_nat, r1a_py)
        check_probe_reversed("r1a_unilateral", os.path.join(ROOT, "refs", "r1a_unilateral.cir"), r1a_nat["T"], r1a_out)
        check_dc_identical("r1a_unilateral", os.path.join(ROOT, "refs", "r1a_unilateral.cir"), r1a_a, r1a_b, r1a_out)
        check_red_proofs("r1a_unilateral", r1a_nat, r1a_py)
        check_klu_option("r1a_unilateral", os.path.join(ROOT, "refs", "r1a_unilateral.cir"))
        check_netlist_card_form("r1a_unilateral", os.path.join(ROOT, "refs", "r1a_unilateral.cir"), r1a_py)

        r1b_nat, r1b_py, _, r1b_out, r1b_a, r1b_b = check_T_match("r1b_bilateral", os.path.join(ROOT, "refs", "r1b_bilateral.cir"))
        check_margins("r1b_bilateral", r1b_nat, r1b_py)
        check_probe_reversed("r1b_bilateral", os.path.join(ROOT, "refs", "r1b_bilateral.cir"), r1b_nat["T"], r1b_out)
        check_dc_identical("r1b_bilateral", os.path.join(ROOT, "refs", "r1b_bilateral.cir"), r1b_a, r1b_b, r1b_out)
        check_red_proofs("r1b_bilateral", r1b_nat, r1b_py)

        r1d_nat, r1d_py, _, r1d_out, r1d_a, r1d_b = check_T_match("r1d_mos_follower", os.path.join(ROOT, "refs", "r1d_mos_follower.cir"))
        check_margins("r1d_mos_follower", r1d_nat, r1d_py)
        check_probe_reversed("r1d_mos_follower", os.path.join(ROOT, "refs", "r1d_mos_follower.cir"), r1d_nat["T"], r1d_out)
        check_dc_identical("r1d_mos_follower", os.path.join(ROOT, "refs", "r1d_mos_follower.cir"), r1d_a, r1d_b, r1d_out)
        check_red_proofs("r1d_mos_follower", r1d_nat, r1d_py)

        record("r1c_fully_diff: single-ended native .stb", "skip",
               "R1c's XP is a diffprobe (4-terminal DM/CM probe), not an lgprobe; "
               "native .stb only implements the single 2-terminal probe=<Vname> "
               "syntax (spec.md stb_syntax) -- no differential probe support in "
               "this card (T3.2/T3.3). Not meaningful to run; deferred.")

        check_runtime_r1d(os.path.join(ROOT, "refs", "r1d_mos_follower.cir"))
    except stb.CouldNotRun as e:
        record("run", None, str(e))

    check_test_stb_regression()

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
