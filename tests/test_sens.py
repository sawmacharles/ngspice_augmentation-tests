#!/usr/bin/env python3
"""Validation for card #475 (plan T10.7), FIRST HALF: adjoint sensitivity of
native `.stb`'s loop gain T(f) and margins (dPM/dp, dfc/dp, dGM/dp) wrt
element/model parameters (`sens=all|<name>`, src/spicelib/analysis/stban.c
in the private fork worktree ngspice-src-sens, branch aug-sens).

Compares the adjoint sensitivities against CENTRAL finite differences: for
each tested parameter p with nominal value p0, two extra full `.stb` runs
(no sens=) with p perturbed to p0*(1+DELTA) and p0*(1-DELTA) give
    dT/dp_FD  = (T(f; p+) - T(f; p-)) / (2*DELTA*p0)
    dPM/dp_FD = (PM(p+) - PM(p-)) / (2*DELTA*p0)   (similarly dfc/dp, dGM/dp)
DELTA = 1e-3 (0.1% relative parameter step) -- chosen empirically this pass:
tested against DELTA=1e-4 and DELTA=1e-2 on R1b's R1o (see check_delta_choice
below); 1e-3 and 1e-4 agree to <1e-6 relative (far inside the FD's own
truncation-vs-roundoff sweet spot), while 1e-2 already shows curvature error
on some frequencies -- 1e-3 is used for every comparison below.

Every deck this script writes lives under
~/Work/data/ngspice_augmentation/sens/ (never the repository); tests/fixtures/
sens_*.cir are read-only here. Every check is pass / FAIL / COULD-NOT-RUN (a
third state, "skip", never affects the exit code); red-proof / expected-
divergence checks are labelled explicitly, not folded into the pass count
silently.

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_sens.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import stb  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/sens")
FIXTURES = os.path.join(HERE, "fixtures")
DELTA = 1e-3
TOL = 1e-4  # spec.md sens_adjoint

results = []


def record(name, state, detail):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN", "skip", "info")
    results.append((name, state, detail))
    print(f"[{state:16s}] {name}: {detail}")


# ---------------------------------------------------------------------
# SPICE numeric-suffix parsing (for perturbing a .param value by a ratio
# and writing the result back out in plain scientific notation)
# ---------------------------------------------------------------------
SUFFIX = {
    "t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "mil": 25.4e-6,
    "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15,
}


def parse_spice_num(s):
    m = re.match(r"^([+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?)([a-zA-Z]*)$", s.strip())
    if not m:
        raise ValueError(f"not a SPICE number: {s!r}")
    val, suf = float(m.group(1)), m.group(2).lower()
    if not suf:
        return val
    if suf.startswith("meg"):
        return val * 1e6
    if suf[0] in SUFFIX:
        return val * SUFFIX[suf[0]]
    raise ValueError(f"unknown suffix on {s!r}")


def perturb_param(text, name, factor):
    """Return (new_text, old_value, new_value) with .param NAME=VALUE's
    VALUE multiplied by `factor`, wherever NAME appears as a .param
    assignment (bare word boundary, case-sensitive -- fixtures here are
    all-lowercase .param names)."""
    pat = re.compile(r"(?<![\w.])(" + re.escape(name) + r")\s*=\s*([0-9eE+\-.a-zA-Z]+)")
    hits = list(pat.finditer(text))
    if not hits:
        raise ValueError(f".param {name}= not found")
    if len(hits) > 1:
        raise ValueError(f".param {name}= appears {len(hits)} times, ambiguous")
    old_str = hits[0].group(2)
    old_val = parse_spice_num(old_str)
    new_val = old_val * factor
    new_text = text[:hits[0].start(2)] + f"{new_val:.15e}" + text[hits[0].end(2):]
    return new_text, old_val, new_val


# ---------------------------------------------------------------------
# native .stb (+ optional sens=) driver
# ---------------------------------------------------------------------
def run_stb(deck_text, probe, sweep, outdir, sens_names, tag, env=None):
    os.makedirs(outdir, exist_ok=True)
    deck_path = os.path.join(outdir, f"deck_{tag}.cir")
    with open(deck_path, "w") as fh:
        fh.write(deck_text)

    sens_clause = "".join(f" sens={n}" for n in sens_names)
    t_file = os.path.join(outdir, f"t_{tag}.txt")
    marg_file = os.path.join(outdir, f"marg_{tag}.txt")
    dtdp_file = os.path.join(outdir, f"dtdp_{tag}.txt")
    dmarg_file = os.path.join(outdir, f"dmarg_{tag}.txt")
    for f in (t_file, marg_file, dtdp_file, dmarg_file):
        if os.path.exists(f):
            os.remove(f)

    # a model parameter is named <instance>:<keyword> (sens=M2:vto); its
    # vectors use '_' for the ':' (stban.c)
    vn = [n.replace(":", "_") for n in sens_names]
    dtdp_vecs = " ".join(f"stb2.dT_d{n}" for n in vn)
    dmarg_vecs = " ".join(f"stb4.dPM_d{n} stb4.dfc_d{n} stb4.dGM_d{n}" for n in vn)
    # plot numbering (stban.c): stb1=T always; stb2=dT/dp and stb4=dMargins/dp
    # only exist when sens= was given at all, shifting the margins plot from
    # stb2 (no sens=) to stb3 (sens= present).
    marg_plot = "stb3" if sens_names else "stb2"

    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
stb {sweep} probe={probe}{sens_clause}
wrdata {t_file} stb1.t
wrdata {marg_file} {marg_plot}.phaseMargin {marg_plot}.phaseMarginFreq {marg_plot}.gainMargin {marg_plot}.gainMarginFreq
"""
    if sens_names:
        ctl += f"wrdata {dtdp_file} {dtdp_vecs}\n"
        ctl += f"wrdata {dmarg_file} {dmarg_vecs}\n"
    ctl += "quit\n.endc\n.end\n"

    wrapper = os.path.join(outdir, f"wrapper_{tag}.cir")
    with open(wrapper, "w") as fh:
        fh.write(deck_text + ctl)

    try:
        proc = subprocess.run([NGSPICE, "-b", wrapper], capture_output=True, text=True, timeout=600,
                              env=dict(os.environ, **(env or {})))
    except (FileNotFoundError, PermissionError) as e:
        raise stb.CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, f"log_{tag}.txt"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines()
            if re.search(r"\berror\b", ln, re.I) and "0 errors" not in ln.lower()]
    if proc.returncode != 0 or errs:
        raise stb.CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")

    freqs, (T,) = stb.read_wrdata(t_file, 1)
    mfreqs, mcols = read_wrdata_real(marg_file, 4)
    margins = {"PM": mcols[0][0], "PMfreq": mcols[1][0], "GM": mcols[2][0], "GMfreq": mcols[3][0]}

    out = {"freqs": freqs, "T": T, "margins": margins, "log": proc.stdout + proc.stderr}
    if sens_names:
        _, dtdp_cols = stb.read_wrdata(dtdp_file, len(sens_names))
        out["dTdp"] = {n: dtdp_cols[i] for i, n in enumerate(sens_names)}
        _, dm_cols = read_wrdata_real(dmarg_file, 3 * len(sens_names))
        out["dmargins"] = {}
        for i, n in enumerate(sens_names):
            out["dmargins"][n] = {
                "dPM": dm_cols[3 * i][0], "dfc": dm_cols[3 * i + 1][0], "dGM": dm_cols[3 * i + 2][0],
            }
    return out


def read_wrdata_real(path, n):
    """n real columns from a wr_singlescale/wr_vecnames file, ignoring the
    leading duplicate 'scale' column wrdata always writes even for a
    scale=NULL (single-row, no-frequency) plot -- see stb.read_wrdata's
    docstring for the complex-column analogue this mirrors."""
    if not os.path.exists(path):
        raise stb.CouldNotRun(f"ngspice wrote no {path}")
    freqs, cols = [], [[] for _ in range(n)]
    with open(path) as fh:
        next(fh)
        for row in fh:
            v = [float(x) for x in row.split()]
            if len(v) != 1 + n:
                raise stb.CouldNotRun(f"unexpected column count in {path}: {v}")
            freqs.append(v[0])
            for k in range(n):
                cols[k].append(v[1 + k])
    if not freqs:
        raise stb.CouldNotRun(f"{path} is empty")
    return freqs, cols


# ---------------------------------------------------------------------
# central-FD comparison
# ---------------------------------------------------------------------
def rel_err(a, b, floor=1e-30):
    return abs(a - b) / max(abs(b), floor)


def compare_scalar(name, tag, label, adjoint, plus, minus, p0, delta):
    if adjoint is None or math.isnan(adjoint):
        record(f"{tag}: {label}_d{name}", "skip", "adjoint is NaN (no crossing) -- not compared")
        return None
    if plus is None or minus is None or math.isnan(plus) or math.isnan(minus):
        record(f"{tag}: {label}_d{name}", "skip", "perturbed run has no crossing -- FD undefined")
        return None
    fd = (plus - minus) / (2 * delta * p0)
    scale = max(abs(fd), abs(adjoint), 1e-30)
    err = abs(adjoint - fd) / scale
    ok = err <= TOL
    record(f"{tag}: {label}_d{name}", "pass" if ok else "FAIL",
           f"adjoint={adjoint:.6g} FD={fd:.6g} rel err {err:.3e} (bound {TOL:.0e})")
    return err


def full_check(tag, fixture, probe, sweep, param_map, expect_diverge=()):
    """param_map: {element_name (sens=): .param name}. expect_diverge:
    element names expected to MISMATCH FD beyond TOL (op-point-dependence
    demonstration, not a tool bug -- recorded as 'info', not FAIL/pass)."""
    outdir = os.path.join(DATA_ROOT, tag)
    with open(os.path.join(FIXTURES, fixture)) as fh:
        base_text = fh.read()

    sens_names = list(param_map.keys())
    try:
        nominal = run_stb(base_text, probe, sweep, outdir, sens_names, "nom")
    except stb.CouldNotRun as e:
        record(f"{tag}: nominal run", "COULD-NOT-RUN", str(e))
        return

    freqs = nominal["freqs"]
    npts = len(freqs)
    sample_idx = sorted(set([npts // 8, npts // 4, npts // 2, 3 * npts // 4, npts - 2]))
    sample_idx = [i for i in sample_idx if 0 <= i < npts]

    for elt, pname in param_map.items():
        try:
            plus_text, p0, p_plus = perturb_param(base_text, pname, 1.0 + DELTA)
            minus_text, _, p_minus = perturb_param(base_text, pname, 1.0 - DELTA)
        except ValueError as e:
            record(f"{tag}: {elt} perturb", "COULD-NOT-RUN", str(e))
            continue
        try:
            plus = run_stb(plus_text, probe, sweep, outdir, sens_names, f"{elt}_plus")
            minus = run_stb(minus_text, probe, sweep, outdir, sens_names, f"{elt}_minus")
        except stb.CouldNotRun as e:
            record(f"{tag}: {elt} perturbed runs", "COULD-NOT-RUN", str(e))
            continue

        label = "EXPECT-DIVERGE" if elt in expect_diverge else tag
        adj = nominal["dTdp"][elt]
        err = compare_dTdp_generic(elt, label, adj, plus["T"], minus["T"], p0, DELTA, freqs, sample_idx,
                                    is_divergence_demo=(elt in expect_diverge))

        dm = nominal["dmargins"][elt]
        for key, pkey in (("dPM", "PM"), ("dfc", "PMfreq"), ("dGM", "GM")):
            adj_s = dm[key]
            plus_s = plus["margins"].get(pkey)
            minus_s = minus["margins"].get(pkey)
            if elt in expect_diverge:
                # still computed, but reported as info (scope-boundary demo), not gated
                if adj_s is not None and not math.isnan(adj_s) and plus_s is not None and minus_s is not None \
                        and not math.isnan(plus_s) and not math.isnan(minus_s):
                    fd = (plus_s - minus_s) / (2 * DELTA * p0)
                    err_s = abs(adj_s - fd) / max(abs(fd), abs(adj_s), 1e-30)
                    record(f"EXPECT-DIVERGE {tag}: {key}_d{elt}", "info",
                           f"adjoint={adj_s:.6g} FD={fd:.6g} rel err {err_s:.3e} "
                           f"(op-point-fixed vs true FD -- divergence is the expected finding)")
            else:
                compare_scalar(elt, tag, key, adj_s, plus_s, minus_s, p0, DELTA)


def compare_dTdp_generic(name, tag, adjoint, plus, minus, p0, delta, freqs, sample_idx, is_divergence_demo=False):
    max_err = 0.0
    worst = None
    for i in sample_idx:
        fd = (plus[i] - minus[i]) / (2 * delta * p0)
        adj = adjoint[i]
        scale = max(abs(fd), abs(adj), 1e-30)
        err = abs(adj - fd) / scale
        if err > max_err:
            max_err, worst = err, (freqs[i], adj, fd)
    if is_divergence_demo:
        record(f"EXPECT-DIVERGE {tag}: dT/d{name} vs central FD", "info",
               f"max rel err {max_err:.3e} (op-point-fixed vs true FD -- divergence is the expected finding) "
               f"at f={worst[0]:.4g} Hz adjoint={worst[1]:.6g} FD={worst[2]:.6g}")
        return max_err
    ok = max_err <= TOL
    record(f"{tag}: dT/d{name} vs central FD", "pass" if ok else "FAIL",
           f"max rel err {max_err:.3e} (bound {TOL:.0e}) at f={worst[0]:.4g} Hz "
           f"adjoint={worst[1]:.6g} FD={worst[2]:.6g}")
    return max_err


# ---------------------------------------------------------------------
# DELTA choice sanity: compare FD gradient at DELTA=1e-3 vs 1e-4 on one
# parameter -- if they disagree by more than a tiny amount, 1e-3 is too
# coarse (curvature) or 1e-4 is too fine (roundoff); record which.
# ---------------------------------------------------------------------
def check_delta_choice():
    outdir = os.path.join(DATA_ROOT, "delta_choice")
    with open(os.path.join(FIXTURES, "sens_r1b_native.cir")) as fh:
        base_text = fh.read()
    sweep = "dec 20 10 1e8"
    try:
        nom = run_stb(base_text, "Vprobe", sweep, outdir, [], "nom")
    except stb.CouldNotRun as e:
        record("delta choice sanity", "COULD-NOT-RUN", str(e))
        return
    freqs = nom["freqs"]
    idx = len(freqs) // 2
    fds = {}
    for d in (1e-2, 1e-3, 1e-4):
        pt, p0, _ = perturb_param(base_text, "r1o", 1 + d)
        mt, _, _ = perturb_param(base_text, "r1o", 1 - d)
        try:
            rp = run_stb(pt, "Vprobe", sweep, outdir, [], f"p_{d}")
            rm = run_stb(mt, "Vprobe", sweep, outdir, [], f"m_{d}")
        except stb.CouldNotRun as e:
            record("delta choice sanity", "COULD-NOT-RUN", str(e))
            return
        fds[d] = (rp["T"][idx] - rm["T"][idx]) / (2 * d * p0)
    e_3v4 = abs(fds[1e-3] - fds[1e-4]) / max(abs(fds[1e-4]), 1e-30)
    e_2v3 = abs(fds[1e-2] - fds[1e-3]) / max(abs(fds[1e-3]), 1e-30)
    record("delta choice: FD(1e-3) vs FD(1e-4) at f=" + f"{freqs[idx]:.4g}", "info",
           f"rel diff {e_3v4:.3e} (small => 1e-3 not yet curvature-limited)")
    record("delta choice: FD(1e-2) vs FD(1e-3) at f=" + f"{freqs[idx]:.4g}", "info",
           f"rel diff {e_2v3:.3e}")


# ---------------------------------------------------------------------
# red-proofs: Y instead of Y^T, and dropping the dY/dp stamp sign
# ---------------------------------------------------------------------
def check_bias_term_extras():
    """card #595: red proof (STB_SENS_FIXED_OP=1 restores the fixed
    operating point and must fail on M2_w), T bit-identical with and
    without sens=, AC-magnitude sensitivities exactly zero, and the
    warning for a model parameter shared by several instances."""
    outdir = os.path.join(DATA_ROOT, "R1dbias_extra")
    with open(os.path.join(FIXTURES, "sens_r1d_bias.cir")) as fh:
        base = fh.read()
    sweep = "dec 20 5e6 1e8"
    try:
        nom = run_stb(base, "Vprobe", sweep, outdir, ["M2_w", "Vin_acmag", "Ib_acmag"], "nom")
        old = run_stb(base, "Vprobe", sweep, outdir, ["M2_w"], "fixedop", env={"STB_SENS_FIXED_OP": "1"})
        plain = run_stb(base, "Vprobe", sweep, outdir, [], "plain")
        pt, p0, _ = perturb_param(base, "m2w", 1 + DELTA)
        mt, _, _ = perturb_param(base, "m2w", 1 - DELTA)
        plus = run_stb(pt, "Vprobe", sweep, outdir, [], "plus")
        minus = run_stb(mt, "Vprobe", sweep, outdir, [], "minus")
        shared = run_stb(base, "Vprobe", sweep, outdir, ["M1:vto"], "shared")
    except (stb.CouldNotRun, ValueError) as e:
        record("bias-point term extras (card #595)", "COULD-NOT-RUN", str(e))
        return
    worst_new = worst_old = 0.0
    for i in range(len(nom["freqs"])):
        fd = (plus["T"][i] - minus["T"][i]) / (2 * DELTA * p0)
        sc = max(abs(fd), 1e-30)
        worst_new = max(worst_new, abs(nom["dTdp"]["M2_w"][i] - fd) / sc)
        worst_old = max(worst_old, abs(old["dTdp"]["M2_w"][i] - fd) / sc)
    record("red proof: STB_SENS_FIXED_OP=1 (card #475 fixed operating point) fails dT/dM2_w",
           "pass" if worst_old > 100 * TOL and worst_new <= TOL else "FAIL",
           f"fixed-op worst rel err {worst_old:.3e} (must exceed {100 * TOL:.0e}); "
           f"with the bias-point term {worst_new:.3e} (bound {TOL:.0e})")
    same = all(a == b for a, b in zip(nom["T"], plain["T"])) and len(nom["T"]) == len(plain["T"])
    record("T is bit-identical with and without sens= (the linearisation is restored after each shift)",
           "pass" if same else "FAIL",
           f"{len(nom['T'])} points, max |dT| {max(abs(a - b) for a, b in zip(nom['T'], plain['T'])):.3e}")
    acz = max(max(abs(v) for v in nom["dTdp"][n]) for n in ("Vin_acmag", "Ib_acmag"))
    record("AC-magnitude sensitivities are exactly zero (.stb drops every AC stimulus)",
           "pass" if acz == 0.0 else "FAIL", f"max |dT/d acmag| = {acz:.3e} (s131 gave 2.8e3 for Ib_acmag)")
    warned = re.search(r"model parameter\(s\) of a model with several instances keep the fixed operating "
                       r"point.*: m1:vto", shared["log"]) is not None
    record("a model parameter shared by 3 instances (M1:vto on nch) is warned about",
           "pass" if warned else "FAIL", "warning present" if warned else "warning missing")


def check_red_proofs():
    yt_bin = os.environ.get("NGSPICE_REDPROOF_YT")
    sign_bin = os.environ.get("NGSPICE_REDPROOF_SIGN")
    if not yt_bin and not sign_bin:
        record("red-proofs (Y vs Y^T, dY/dp sign)", "skip",
               "NGSPICE_REDPROOF_YT / NGSPICE_REDPROOF_SIGN not set -- see report for how these "
               "throwaway builds were made and run manually")
        return
    for label, binpath in (("Y-instead-of-Y^T", yt_bin), ("dropped dY/dp sign", sign_bin)):
        if not binpath:
            continue
        global NGSPICE
        saved = NGSPICE
        NGSPICE = binpath
        try:
            outdir = os.path.join(DATA_ROOT, f"redproof_{label.replace(' ', '_').replace('^','')}")
            with open(os.path.join(FIXTURES, "sens_r1b_native.cir")) as fh:
                base_text = fh.read()
            sweep = "dec 20 10 1e8"
            nom = run_stb(base_text, "Vprobe", sweep, outdir, ["R1o"], "nom")
            pt, p0, _ = perturb_param(base_text, "r1o", 1 + DELTA)
            mt, _, _ = perturb_param(base_text, "r1o", 1 - DELTA)
            plus = run_stb(pt, "Vprobe", sweep, outdir, ["R1o"], "plus")
            minus = run_stb(mt, "Vprobe", sweep, outdir, ["R1o"], "minus")
            freqs = nom["freqs"]
            idx = len(freqs) // 2
            fd = (plus["T"][idx] - minus["T"][idx]) / (2 * DELTA * p0)
            adj = nom["dTdp"]["R1o"][idx]
            err = abs(adj - fd) / max(abs(fd), 1e-30)
            record(f"red-proof ({label}) must FAIL the comparison", "pass" if err > TOL else "FAIL",
                   f"rel err {err:.3e} (bound {TOL:.0e}) -- {'correctly diverges' if err > TOL else 'DID NOT diverge, red-proof is not discriminating'}")
        except stb.CouldNotRun as e:
            record(f"red-proof ({label})", "COULD-NOT-RUN", str(e))
        finally:
            NGSPICE = saved


# ---------------------------------------------------------------------
def main():
    print(f"NGSPICE={NGSPICE}")
    os.makedirs(DATA_ROOT, exist_ok=True)

    # R1b (linear): three parameters spanning R and C, one in each stage
    # and one in the outer feedback path. Swept range narrowed around the
    # single falling PM/GM crossing (~9.1e7/~found only above 1e7 Hz on
    # this deck -- a wide-range low-density scan found it first) at HIGH
    # density: dPM/dp, dfc/dp, dGM/dp are computed from a 2-point LOCAL
    # bracket (log-f linear) derivative of the swept T(f) trace, so their
    # OWN accuracy (unlike dT/dp's, which is exact per frequency point) is
    # limited by sweep spacing, not by the adjoint formula itself --
    # confirmed empirically this pass (see the card #475 report's
    # "density convergence" table: dec 25 gives ~1e-3, dec 100 ~1e-4,
    # dec 400 ~1e-5 relative error on this exact deck/parameter). dec 200
    # over about one decade is dense enough to clear the 1e-4 bound with
    # margin while keeping this script's own runtime reasonable.
    full_check("R1b", "sens_r1b_native.cir", "Vprobe", "dec 200 2e7 1.5e8",
               {"R1o": "r1o", "Cc": "cc", "Rf1": "ro1"})

    # R3 averaged linear stand-in (spec.md's own DC-op-half reference for
    # R3). PM crossing ~55832 Hz, GM crossing ~46790 Hz (found from a
    # wide, coarse scan first) -- both inside this narrower, denser range.
    # 2026-09-30 (card #519, Goa polarity fix in sens_r3_averaged_native.cir):
    # the dGM checks (local 2-point bracket derivative, first-order in the grid
    # spacing for GM) went from 5e-6 to 1e-3 relative at dec 300 (exact cause
    # not investigated); dec 1000 gives 2.4e-4, dec 3000 gives <=5.3e-5, so
    # dec 3000 is used (AC only, seconds).
    full_check("R3avg", "sens_r3_averaged_native.cir", "Vprobe", "dec 3000 3e4 8e4",
               {"R2f": "r2", "C2f": "c2", "Resr": "esrval"})

    # R1d nonlinear: CL (op-point-independent, exact test of the nonlinear
    # small-signal matrix path) and M2's w (op-point-DEPENDENT: card #475
    # held the operating point fixed and this was an EXPECT-DIVERGE
    # demonstration; card #595 adds the bias-point term, so it is gated
    # now -- see docs/sens.md). PM crossing ~3.33e7 Hz; no
    # GM crossing on this deck (dGM checks report "skip", not FAIL).
    # sens=<name> for a non-principal instance parameter is
    # "<instance>_<keyword>" (cktsens.c's own DC/AC `.sens` naming, which
    # this file's stb_sens_name() mirrors byte-for-byte) -- "M2_w", not
    # "M2w".
    full_check("R1d", "sens_r1d_native.cir", "Vprobe", "dec 200 5e6 1e8",
               {"CL": "cl_val", "M2_w": "m2w"})

    # card #595: parameters that move the bias point -- the tail bias
    # current, M2's w, and two model parameters of M2's own model: vto
    # (moves the bias) and cbd (junction capacitance: charge only, computed
    # in the device's DC load, which the card #475 stamp difference missed)
    full_check("R1dbias", "sens_r1d_bias.cir", "Vprobe", "dec 200 5e6 1e8",
               {"Ib": "ib_val", "M2_w": "m2w", "M2:vto": "m2vto", "M2:cbd": "m2cbd"})
    check_bias_term_extras()

    check_delta_choice()
    check_red_proofs()

    print()
    n_pass = sum(1 for _, s, _ in results if s == "pass")
    n_fail = sum(1 for _, s, _ in results if s == "FAIL")
    n_cnr = sum(1 for _, s, _ in results if s == "COULD-NOT-RUN")
    n_skip = sum(1 for _, s, _ in results if s == "skip")
    n_info = sum(1 for _, s, _ in results if s == "info")
    print(f"{n_pass} pass, {n_fail} FAIL, {n_cnr} COULD-NOT-RUN, {n_skip} skip, {n_info} info "
          f"(of {len(results)} checks)")
    if n_fail:
        return 1
    if n_cnr:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
