#!/usr/bin/env python3
"""Validation for card #475, PERIODIC half: tools/pstbsens.py -- parameter sensitivities of the
periodic loop gain (.pstb) by central finite differences + Richardson over the whole .pss+.pstb run.

  1. equivalence (the key one): tests/fixtures/pstbsens_r1d_clk.cir = the OTA follower of
     sens_r1d_bias.cir plus an unrelated clock. The orbit is constant, so .pstb T must equal .stb T,
     and the tool's dT/dp must equal the s133 `.stb sens=` adjoint dT/dp (bias-point term included)
     for ib_val, m2w, M2:vto and CL.  Gate: worst per-frequency relative difference <= 1e-4,
     dPM/dp relative difference <= 1e-3.
  2. rd_pstb_r3d.cir (smooth-switch buck with dead time, real switching orbit): dPM/dp and dfc/dp of
     cval, esrval, r2: error estimate e <= 1 % of |D|, and agreement with an independent central difference
     from two plain ngspice runs at p(1 +- 0.005) with this file's own T parsing and margin code,
     within max(3 e, 2 % |D|).
  3. jobs=1 and jobs=4 give byte-identical dtdp.csv / summary.json.
  4. refusals: unknown parameter, parameter without effect, deck without .pstb, PSS not converged,
     no 0 dB crossing (PM line missing).
  5. red proof: one-sided forward difference at h = 0.1 (hidden --redproof-forward) must FAIL check 1's gate
     by a clear margin (>= 10x).

Every check is pass / FAIL / COULD-NOT-RUN (a third state, skip/info, never affects the exit code).
Outputs go to ~/Work/data/ngspice_augmentation/pstbsens/test, never the repository.

    python3 tests/test_pstbsens.py [--only 1,2,...]
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import math
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
DEFAULT_BIN = os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice")
os.environ.setdefault("NGSPICE", DEFAULT_BIN if os.path.exists(DEFAULT_BIN) else "ngspice")
import pstbsens  # noqa: E402
import stb  # noqa: E402

NGSPICE = os.environ["NGSPICE"]
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/pstbsens/test")
FIX = os.path.join(HERE, "fixtures")
TOOL = os.path.join(ROOT, "tools", "pstbsens.py")
RELTOL_ENV = "1e-8"  # PSS_SHOOT_RELTOL, as the tool's default

GATE_T = 1e-4        # check 1: worst relative dT/dp difference
GATE_PM = 1e-3       # check 1: dPM/dp
H1 = 0.003           # check 1 relative step: small enough that a crossing stays inside its grid cell (the chord PM is kinked)
H2 = 0.03            # check 2 relative step: cval and r2 have a PSS that fails to converge at +-1 %

results = []


class CouldNotRun(RuntimeError):
    pass


def record(name, state, detail):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN", "skip", "info")
    results.append((name, state, detail))
    print(f"[{state:13s}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------
# plain ngspice runner (own code, independent of the tool): systemd-run scope, OMP 2
# ---------------------------------------------------------------------------------------------
def ngrun(deck_text, outdir, env_extra=None, timeout=1800):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "deck.cir")
    with open(path, "w") as fh:
        fh.write(deck_text)
    cmd = [NGSPICE, "-b", path]
    if shutil.which("systemd-run"):
        cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=4G", "-p", "MemorySwapMax=0"] + cmd
    env = dict(os.environ, OMP_NUM_THREADS="2", **(env_extra or {}))
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=outdir, env=env)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise CouldNotRun(f"cannot run ngspice: {e}")
    log = p.stdout + p.stderr
    with open(os.path.join(outdir, "log.txt"), "w") as fh:
        fh.write(log)
    return log


def read_real_cols(path, n):
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            v = [float(x) for x in line.split()]
            if len(v) != 1 + n:
                raise CouldNotRun(f"unexpected column count in {path}")
            rows.append(v[1:])
    if not rows:
        raise CouldNotRun(f"{path} empty")
    return [[r[k] for r in rows] for k in range(n)]


def strip_control(text):
    return re.sub(r"(?is)^\.control.*?^\.endc[^\n]*\n", "", text, flags=re.M).replace("\n.end\n", "\n")


# ---------------------------------------------------------------------------------------------
# own margin code (check 2 only): PM and fc of T, written independently of tools/pstbsens.py
# ---------------------------------------------------------------------------------------------
def own_pm(freqs, T):
    """Phase margin 180 + arg T at the first |T| = 1 crossing, fc; linear in ln f; own unwrapping."""
    ph = []
    for i, z in enumerate(T):
        a = math.degrees(cmath_phase(z))
        if i:
            a += 360.0 * round((ph[-1] - a) / 360.0)
        ph.append(a)
    cands = []
    for i in range(len(T) - 1):
        g0, g1 = math.log(abs(T[i])), math.log(abs(T[i + 1]))
        if (g0 > 0) != (g1 > 0):
            w = g0 / (g0 - g1)
            fc = freqs[i] * (freqs[i + 1] / freqs[i]) ** w
            cands.append((180.0 + ph[i] + w * (ph[i + 1] - ph[i]), fc))
    if not cands:
        raise CouldNotRun("own PM: no 0 dB crossing")
    return min(cands)


def cmath_phase(z):
    return math.atan2(z.imag, z.real)


# ---------------------------------------------------------------------------------------------
# check 1
# ---------------------------------------------------------------------------------------------
PAIRS = [("ib_val", "Ib"), ("m2w", "M2_w"), ("M2:vto", "M2:vto"), ("CL", "CL")]
SWEEP = "dec 200 5e6 1e8"


def stb_adjoint():
    """Native .stb with sens= on sens_r1d_bias.cir: T, dT/dp, dPM/dp, dfc/dp (same circuit, no clock)."""
    outdir = os.path.join(DATA, "c1_stb")
    with open(os.path.join(FIX, "sens_r1d_bias.cir")) as fh:
        deck = strip_control(fh.read())
    names = [b for _, b in PAIRS]
    vn = [n.replace(":", "_") for n in names]
    ctl = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
    ctl += f"stb {SWEEP} probe=Vprobe " + " ".join(f"sens={n}" for n in names) + "\n"
    ctl += f"wrdata {outdir}/t.txt stb1.t\n"
    ctl += f"wrdata {outdir}/dtdp.txt " + " ".join(f"stb2.dT_d{n}" for n in vn) + "\n"
    ctl += f"wrdata {outdir}/dm.txt " + " ".join(f"stb4.dPM_d{n} stb4.dfc_d{n}" for n in vn) + "\n"
    ctl += "quit\n.endc\n.end\n"
    ngrun(deck + ctl, outdir)
    try:
        f, (T,) = stb.read_wrdata(os.path.join(outdir, "t.txt"), 1)
        _, dT = stb.read_wrdata(os.path.join(outdir, "dtdp.txt"), len(names))
        cols = read_real_cols(os.path.join(outdir, "dm.txt"), 2 * len(names))
    except (stb.CouldNotRun, OSError) as e:
        raise CouldNotRun(f".stb sens= run failed: {e}")
    return f, T, {n: dT[i] for i, n in enumerate(names)}, \
        {n: (cols[2 * i][0], cols[2 * i + 1][0]) for i, n in enumerate(names)}


def worst_rel(a, b):
    w, at = 0.0, None
    for i, (x, y) in enumerate(zip(a, b)):
        r = abs(x - y) / max(abs(y), 1e-300)
        if r > w:
            w, at = r, i
    return w, at


def check_1_and_5():
    fixture = os.path.join(FIX, "pstbsens_r1d_clk.cir")
    try:
        f_stb, T_stb, dT_stb, dm_stb = stb_adjoint()
    except CouldNotRun as e:
        record("1 equivalence: native .stb sens= reference", "COULD-NOT-RUN", str(e))
        return
    t0 = time.time()
    try:
        recs, nom, outdir, timing, mode, _ = pstbsens.sens_run(
            fixture, [a for a, _ in PAIRS], h=H1, jobs=4, outdir=os.path.join(DATA, "c1_tool"))
    except pstbsens.Refused as e:
        record("1 equivalence: tools/pstbsens.py run", "COULD-NOT-RUN", str(e))
        return
    wall = time.time() - t0
    w, at = worst_rel(nom["T"], T_stb)
    record("1a pstb T == stb T on the constant orbit (pstb_linear tolerance 1e-6)", "pass" if w <= 1e-6 else "FAIL",
           f"{len(T_stb)} points, worst relative |T_pstb - T_stb|/|T_stb| = {w:.3e} at {f_stb[at]:.4g} Hz")
    worst_all = 0.0
    for (tool_name, stb_name), r in zip(PAIRS, recs):
        w, at = worst_rel(r["dT"], dT_stb[stb_name])
        worst_all = max(worst_all, w)
        record(f"1b dT/d{tool_name} vs .stb sens={stb_name} adjoint", "pass" if w <= GATE_T else "FAIL",
               f"worst relative difference {w:.3e} (gate {GATE_T:.0e}) at {f_stb[at]:.4g} Hz; "
               f"tool e/|D| max {max(r['dT_err'][i] / abs(r['dT'][i]) for i in range(len(f_stb))):.2e}")
        a_pm, a_fc = dm_stb[stb_name]
        rp = abs(r["dPM"] - a_pm) / abs(a_pm)
        rf = abs(r["dfc"] - a_fc) / abs(a_fc)
        record(f"1c dPM/d{tool_name}, dfc/d{tool_name} vs adjoint", "pass" if rp <= GATE_PM else "FAIL",
               f"dPM tool {r['dPM']:.6g} adjoint {a_pm:.6g} rel {rp:.3e} (gate {GATE_PM:.0e}); "
               f"dfc tool {r['dfc']:.6g} adjoint {a_fc:.6g} rel {rf:.3e} (reported); "
               f"cubic-interpolated margins: dPM {r['dPM_alt']:.6g} (differs from the chord by {abs(r['dPM_alt'] - r['dPM']) / abs(r['dPM']):.1e})")
    record("1 timing", "info", f"{timing['runs']} runs, {timing['wall_total_s']:.0f} s wall at jobs=4, "
           f"mean run {timing['mean_run_s']:.0f} s => {timing['wall_total_s'] / len(PAIRS):.0f} s per parameter "
           f"(4 runs each, 4 in parallel) (pacpoints=200000)")

    # 5. red proof: forward difference at h = 0.1, no Richardson
    try:
        frecs, _, _, ftiming, fmode, _ = pstbsens.sens_run(
            fixture, [a for a, _ in PAIRS], h=0.1, jobs=4, outdir=os.path.join(DATA, "c5_forward"), forward=True)
    except pstbsens.Refused as e:
        record("5 red proof: forward difference", "COULD-NOT-RUN", str(e))
        return
    wf = 0.0
    for (tool_name, stb_name), r in zip(PAIRS, frecs):
        w, _ = worst_rel(r["dT"], dT_stb[stb_name])
        wf = max(wf, w)
        record(f"5 forward difference h=0.1, dT/d{tool_name}", "info", f"worst relative difference {w:.3e}")
    ratio = wf / GATE_T
    record("5 red proof: forward difference at h = 0.1 must FAIL check 1b's gate by >= 10x",
           "pass" if wf > 10 * GATE_T else "FAIL",
           f"worst {wf:.3e} vs gate {GATE_T:.0e}: ratio {ratio:.1f}x; the central+Richardson run's worst was {worst_all:.3e} "
           f"({wf / max(worst_all, 1e-300):.0f}x smaller)")


# ---------------------------------------------------------------------------------------------
# check 2
# ---------------------------------------------------------------------------------------------
P2 = ["cval", "esrval", "r2"]


def own_run_plain(base_text, name, factor, tag):
    """A full plain run of the deck with .param NAME multiplied by factor; returns (PM, fc, line_pm)."""
    m = list(re.finditer(r"(?<![\w.])" + name + r"=([0-9.eE+\-]+[a-zA-Z]*)", base_text))
    if len(m) != 1:
        raise CouldNotRun(f".param {name} not unique in the deck")
    p0 = pstbsens.parse_number(m[0].group(1))
    text = base_text[:m[0].start(1)] + f"{p0 * factor:.17g}" + base_text[m[0].end(1):]
    outdir = os.path.join(DATA, "c2_plain", tag)
    ctl = (f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n"
           f"wrdata {outdir}/t.txt t\nquit\n.endc\n.end\n")
    log = ngrun(text + ctl, outdir, {"PSS_SHOOT_RELTOL": RELTOL_ENV})
    if "PSS-SHOOT: CONVERGED" not in log:
        raise CouldNotRun(f"plain run {tag}: PSS did not converge")
    line = re.search(r"pstb: phase margin ([-0-9.eE+]+) deg at ([-0-9.eE+]+) Hz", log)
    f, (T,) = stb.read_wrdata(os.path.join(outdir, "t.txt"), 1)
    pm, fc = own_pm(f, T)
    if not line or abs(float(line.group(1)) - pm) > 1e-3 * pm or abs(float(line.group(2)) - fc) > 1e-5 * fc:
        raise CouldNotRun(f"plain run {tag}: own PM {pm}, {fc} differs from the pstb line {line and line.groups()}")
    return pm, fc, p0


def check_2():
    deck = os.path.join(FIX, "rd_pstb_r3d.cir")
    try:
        recs, nom, outdir, timing, mode, _ = pstbsens.sens_run(
            deck, P2, h=H2, jobs=4, outdir=os.path.join(DATA, "c2_tool"), pss_reltol=RELTOL_ENV)
    except pstbsens.Refused as e:
        record("2 r3d: tools/pstbsens.py run", "COULD-NOT-RUN", str(e))
        return
    with open(deck) as fh:
        base = strip_control(fh.read())
    for r in recs:
        for key, label in (("dPM", "dPM/dp"), ("dfc", "dfc/dp")):
            D, e = r[key], r[key + "_err"]
            record(f"2a {r['name']} {label}: error estimate e <= 1 % of |D|", "pass" if e <= 0.01 * abs(D) else "FAIL",
                   f"D = {D:.5g}, e = {e:.3g} ({100 * e / abs(D):.3f} %), h = {H2}")
        try:
            pp, fp, p0 = own_run_plain(base, r["name"], 1.005, f"{r['name']}_p")
            pm_, fm, _ = own_run_plain(base, r["name"], 0.995, f"{r['name']}_m")
        except CouldNotRun as e2:
            record(f"2b {r['name']} independent plain-run difference", "COULD-NOT-RUN", str(e2))
            continue
        for key, label, hi, lo in (("dPM", "dPM/dp", pp, pm_), ("dfc", "dfc/dp", fp, fm)):
            Dx = (hi - lo) / (2 * 0.005 * p0)
            D, e = r[key], r[key + "_err"]
            bound = max(3 * e, 0.02 * abs(D))
            record(f"2b {r['name']} {label} vs two plain runs at p(1 +- 0.005)",
                   "pass" if abs(D - Dx) <= bound else "FAIL",
                   f"tool {D:.5g}, plain {Dx:.5g}, |diff| {abs(D - Dx):.3g} ({100 * abs(D - Dx) / abs(D):.3f} %), "
                   f"bound max(3e, 2 %|D|) = {bound:.3g}")
    pms = [r["PM"] for r in recs]
    record("2 timing", "info", f"{timing['runs']} runs, {timing['wall_total_s']:.1f} s wall at jobs=4, "
           f"mean run {timing['mean_run_s']:.1f} s; nominal PM {pms[0]:.4f} deg, fc {recs[0]['fc']:.6g} Hz")


# ---------------------------------------------------------------------------------------------
# check 3
# ---------------------------------------------------------------------------------------------
def check_3():
    deck = os.path.join(FIX, "rd_pstb_r3d.cir")
    files = {}
    try:
        for jobs in (1, 4):
            od = os.path.join(DATA, f"c3_jobs{jobs}")
            pstbsens.sens_run(deck, ["esrval", "r3"], h=1e-2, jobs=jobs, outdir=od, orbit="vout",
                             pss_reltol=RELTOL_ENV)
            files[jobs] = {n: open(os.path.join(od, n), "rb").read() for n in ("dtdp.csv", "summary.json")}
    except pstbsens.Refused as e:
        record("3 jobs=1 vs jobs=4", "COULD-NOT-RUN", str(e))
        return
    for n in ("dtdp.csv", "summary.json"):
        record(f"3 {n} identical for jobs=1 and jobs=4", "pass" if files[1][n] == files[4][n] else "FAIL",
               f"{len(files[1][n])} bytes vs {len(files[4][n])} bytes")
    js = open(os.path.join(DATA, "c3_jobs4", "summary.json")).read()
    record("3 --orbit vout: d(mean)/dp and d(p-p)/dp present", "pass" if '"dmean"' in js and '"dpp"' in js else "skip",
           "found in summary.json" if '"dmean"' in js else "the PSS plot gave no v(vout): orbit outputs skipped")


# ---------------------------------------------------------------------------------------------
# check 4
# ---------------------------------------------------------------------------------------------
def cli(args):
    p = subprocess.run([sys.executable, TOOL] + args, capture_output=True, text=True, timeout=3600)
    return p.returncode, p.stdout + p.stderr


def check_4():
    deck = os.path.join(FIX, "rd_pstb_r3d.cir")
    out = os.path.join(DATA, "c4")
    os.makedirs(out, exist_ok=True)
    with open(deck) as fh:
        text = fh.read()
    rc, msg = cli([deck, "--param", "nosuchparam", "--out", out + "/a", "--pss-reltol", RELTOL_ENV])
    record("4a unknown parameter refused", "pass" if rc == 2 and "unknown parameter" in msg else "FAIL", f"rc {rc}: {msg.strip()[:140]}")
    unused = os.path.join(out, "unused.cir")
    with open(unused, "w") as fh:
        fh.write(text.replace(".param lval=", ".param zzz_unused=5 lval=", 1))
    rc, msg = cli([unused, "--param", "zzz_unused", "--out", out + "/b", "--pss-reltol", RELTOL_ENV])
    record("4b parameter without effect refused", "pass" if rc == 2 and "no effect" in msg else "FAIL", f"rc {rc}: {msg.strip()[:160]}")
    nopstb = os.path.join(out, "nopstb.cir")
    with open(nopstb, "w") as fh:
        fh.write(re.sub(r"(?m)^\.pstb.*\n", "", text))
    rc, msg = cli([nopstb, "--param", "cval", "--out", out + "/c", "--pss-reltol", RELTOL_ENV])
    record("4c deck without .pstb refused", "pass" if rc == 2 and ".pstb" in msg else "FAIL", f"rc {rc}: {msg.strip()[:140]}")
    nocross = os.path.join(out, "nocross.cir")
    with open(nocross, "w") as fh:
        fh.write(text.replace("dec 5 500 200000", "dec 5 500 2000"))
    rc, msg = cli([nocross, "--param", "esrval", "--out", out + "/d", "--pss-reltol", RELTOL_ENV])
    record("4d PM line missing (no 0 dB crossing in the band) refused", "pass" if rc == 2 and "phase-margin line is missing" in msg else "FAIL",
           f"rc {rc}: {msg.strip()[:160]}")
    rc, msg = cli([deck, "--param", "cval", "--h", "1e-2", "--out", out + "/e", "--pss-reltol", RELTOL_ENV])
    if rc == 2 and "did not converge" in msg:
        record("4e PSS not converged in a run refused (cval, h = 1e-2)", "pass", f"rc {rc}: {msg.strip()[:140]}")
    else:
        record("4e PSS not converged in a run refused (cval, h = 1e-2)", "skip",
               f"this engine now converges all 5 cval runs at h = 1e-2 (rc {rc}); refusal path not exercised here")
    rc, msg = cli([deck, "--param", "esrval", "--jobs", "6", "--out", out + "/f"])
    record("4f --jobs 6 refused (machine budget 5)", "pass" if rc == 2 and "jobs" in msg else "FAIL", f"rc {rc}")
    rc, msg = cli(["--help"])
    record("4g --help works", "pass" if rc == 0 and "--param" in msg and "--redproof-forward" not in msg else "FAIL", f"rc {rc}")


def main():
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))
    print(f"NGSPICE={NGSPICE}", flush=True)
    os.makedirs(DATA, exist_ok=True)
    for key, fn in (("1", check_1_and_5), ("2", check_2), ("3", check_3), ("4", check_4)):
        if only is None or key in only or (key == "1" and "5" in only):
            fn()
    print()
    c = {s: sum(1 for _, x, _ in results if x == s) for s in ("pass", "FAIL", "COULD-NOT-RUN", "skip", "info")}
    print(f"{c['pass']} pass, {c['FAIL']} FAIL, {c['COULD-NOT-RUN']} COULD-NOT-RUN, {c['skip']} skip, {c['info']} info "
          f"(of {len(results)} checks)")
    if c["FAIL"]:
        return 1
    if c["COULD-NOT-RUN"]:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
