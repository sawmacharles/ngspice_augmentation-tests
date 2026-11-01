#!/usr/bin/env python3
"""cards #537, #538: matrix-free Newton-GMRES shooting (`.pss ... method=shoot krylov=1`, or
NGSPICE_PSS_KRYLOV=1) against the dense finite-difference-monodromy path, same binary.

Checks (each writes decks/outputs under ~/Work/data/ngspice_augmentation/538/eq, never the repo):

  equiv <fixture>     dense run vs krylov=1 run (krylov on via the CARD KEY, dense via the
                      environment NGSPICE_PSS_KRYLOV=0, so both entry points are exercised):
                      both converge (or both refuse); the converged orbit agrees -- the
                      orbit's start point (first row of the emitted period; grid independent)
                      and the whole waveform on a common grid, both as max|dev| / peak-to-peak
                      of that vector (the "ripple-relative norm") -- within TOL_ORBIT; the
                      stable/unstable verdict is the same, or krylov says explicitly
                      "UNKNOWN(not-computed)" (never a made-up radius); the Arnoldi spectral
                      radius is within TOL_RHO of the dense one.  Period-run counts of both
                      are reported (info; PSS_SHOOT_COUNT_RUNS=1).
  freerun (card #538) `trigger=` (free-running) shooting is matrix-free as well: R5 hysteretic buck and R5 + extracted RC
                      trees, krylov=1 vs dense: same period (2e-4), same orbit, stable verdict, Arnoldi marker; krylov=0 /
                      AUTO_MIN overrides; automatic from 64 states.  (card #537 had refused it.)
  card_env            krylov=1 on the card, NGSPICE_PSS_KRYLOV=0 in the env -> dense (env wins),
                      and NGSPICE_PSS_KRYLOV=1 with no card key -> krylov.
  pac_after_krylov    `.pac solver=td` after a krylov PSS gives the same transfer curve as after
                      a dense PSS (they need only the orbit state; checked, not assumed).
  pex (--pex)         the PEX-size buck (R3 + 3 RC trees; 490 states, --big 1150): dense vs
                      krylov period runs, wall time, orbit deviation.
  memory (card #553)  the period runs of a `.pss method=shoot` keep no output vectors: peak RSS (getrusage of the child) of the
                      1150-state PEX buck WITHOUT `.save` stays below MEM_BOUND_MB (new build ~40 MB; every run's whole plot
                      used to be retained: 650 MB), the orbit is bit-identical to the run with PSS_KEEP_RUN_DATA=1 (the old
                      retention), and `--redproof` runs the check with that switch and must see it fail.
  identical (NGSPICE_BASE=<baseline binary>)
                      with krylov OFF, this build's stdout (wall_time masked) and wrdata files
                      are byte-identical to the baseline build's on every fixture.
  asan (NGSPICE_ASAN=<asan binary>)
                      two decks with krylov=1 under ASan: exit 0, no AddressSanitizer report.
  auto (card #538)    krylov is chosen automatically from PSS_KRYLOV_AUTO_MIN (default 64) capacitor+inductor
                      states up, with a one-line notice: a 73-state PEX deck without any card runs matrix-free,
                      `krylov=0` on the card or NGSPICE_PSS_KRYLOV=0 forces dense, PSS_KRYLOV_AUTO_MIN=0 turns it
                      off, =73 / =74 are the exact boundary, a 7-state fixture stays dense with no notice, and the
                      automatic orbit equals the dense one.
  warmup (card #538)  from 64 states up the matrix-free run makes no Newton attempt on the first PSS_KRYLOV_WARMUP
                      (default 130) periods of the staged warm-up (the transients and their endpoints are the card
                      #537 ones): same orbit as attempts at every stage (PSS_KRYLOV_WARMUP=0, to 1e-6) and as the
                      dense engine (R3 145 states 1e-5, R4 99 states / R3d 101 states 1e-4), fewer period runs; a
                      deliberately bad start (1 warm-up period, explicit schedule) still converges to the same
                      orbit within the engine's own spread (5e-4).
  --redproof          the krylov deck gets a 1 % different clock guess: the orbit comparison
                      MUST fail, and a deliberately wrong Arnoldi comparison (radius tolerance
                      1e-9) must fail -- proves the checks can fail.

Usage:  NGSPICE=<binary> [NGSPICE_BASE=<baseline>] [NGSPICE_ASAN=<asan>] OMP_NUM_THREADS=2 \
            python3 tests/test_pss_krylov.py [--quick] [--pex] [--big] [--redproof] [--only substr]
Exit 0 pass, 1 fail, 2 could-not-run.
"""
import os, re, sys, subprocess, time, math, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIX = os.path.join(HERE, "fixtures")
NGSPICE = os.environ.get("NGSPICE", "ngspice")
BASE = os.environ.get("NGSPICE_BASE", "")
ASAN = os.environ.get("NGSPICE_ASAN", "")
OUT = os.path.expanduser(os.environ.get("PSS_KRYLOV_TEST_OUT", "~/Work/data/ngspice_augmentation/538/eq"))
RIPPLE_FLOOR = 1e-2   # the engine's own floor under a state's ripple (pssshoot.c ripple_floor, native units)
TOL_ORBIT = 1e-5   # N = max|dev| / max(p-p, RIPPLE_FLOOR): the "ripple-relative norm" of the engine's own
                   # convergence test.  Two runs that both stop at residual 1e-6 of a slow Floquet mode (R4:
                   # rho 0.994, 1/(1-rho) = 166) can differ by ~1e-4 in this norm without either being wrong, so the
                   # pass rule is N(krylov,dense) <= max(TOL_ORBIT, SPREAD_K x spread) where spread is the largest
                   # N(dense, dense variant) over three variants of the DENSE engine itself: FD step (PSS_SHOOT_EPS=2e-5),
                   # no fast-state elimination (PSS_SHOOT_NO_REDUCE=1: all states Newton unknowns, as in krylov mode) and a
                   # different warm-up length (PSS_SHOOT_WARMUP_PERIODS=16).  On a switch deck the period map depends on
                   # the recorded step grid (which stage of the staged warm-up converges): 9e-6 V on R6s, the dense
                   # engine's own spread.  Both numbers are in the report.
SPREAD_K = 3.0
TOL_EDGE = 2e-2    # vectors with edges (pwm, ramp, switch node): the two runs chord the edge on their own grids
TOL_RHO = 5e-3     # |rho_krylov - rho_dense| / rho_dense (both are finite-difference estimates, ~1e-3 apart on R4/R3d)
args = sys.argv[1:]
results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN", "info")
    results.append((name, state, detail))
    print(f"{state:14s} {name}: {detail}", flush=True)


def body(fix):
    with open(os.path.join(FIX, fix)) as fh:
        s = fh.read()
    return s.replace("../../stb/lgprobe.lib", os.path.join(ROOT, "stb", "lgprobe.lib"))


def redirect(deck, d):
    return re.sub(r"(?m)^(\s*(?:wrdata|write))\s+\S*?([^/\s]+)(\s|$)",
                  lambda m: f"{m.group(1)} {d}/{m.group(2)}{m.group(3)}", deck)


def set_krylov_card(deck, on):
    """add `krylov=1` to the (control-command or dot-card) .pss line, or leave the deck alone."""
    if not on:
        return deck
    return re.sub(r"(?m)^(\s*\.?pss\s+.*?)\s*$", lambda m: m.group(1) + " krylov=1", deck, count=1)


def wrdata_files(d):
    return sorted(f for f in os.listdir(d) if f.endswith(".txt") and f not in ("stdout.txt",))


def run(binary, deck_path, env=None, timeout=3600):
    e = dict(os.environ)
    e.setdefault("OMP_NUM_THREADS", "2")
    e["PSS_SHOOT_COUNT_RUNS"] = "1"
    for k in ("NGSPICE_PSS_KRYLOV", "PSS_KRYLOV_ONESIDED", "PSS_KRYLOV_ETA"):
        e.pop(k, None)
    if env:
        e.update(env)
    t0 = time.time()
    try:
        p = subprocess.run([binary, "-b", deck_path], capture_output=True, text=True, timeout=timeout,
                           cwd=os.path.dirname(deck_path), env=e)
    except subprocess.TimeoutExpired:
        raise CouldNotRun(f"{deck_path} timed out")
    except FileNotFoundError as ex:
        raise CouldNotRun(f"cannot start {binary}: {ex}")
    return p.returncode, p.stdout + p.stderr, time.time() - t0


CONV = re.compile(r"PSS-SHOOT: CONVERGED iter=(\d+) residual=(\S+) spectral_radius=(\S+) (\S+) wall_time=(\S+)s")
RUNS = re.compile(r"PSS-SHOOT: period runs=(\d+) \(\+ (\d+) warm-up transient periods\), mode=(\S+)")


def parse(out):
    m = CONV.search(out)
    r = RUNS.search(out)
    if not m:
        return None
    return dict(iter=int(m.group(1)), res=float(m.group(2)), rho=float(m.group(3)), verdict=m.group(4),
                wall=float(m.group(5)), runs=r and int(r.group(1)), warm=r and int(r.group(2)), mode=r and r.group(3))


# ------------------------------------------------------------------ orbit comparison
def read_cols(path):
    """wrdata (no wr_singlescale): rows of (t1 y1 t2 y2 ...) -> [(ts, ys)]"""
    rows = []
    with open(path) as fh:
        for l in fh:
            t = l.split()
            if t:
                try:
                    rows.append([float(v) for v in t])
                except ValueError:
                    pass
    if not rows:
        raise CouldNotRun(f"empty/absent wrdata {path}")
    single = len(rows[0]) % 2 == 1   # wr_singlescale: t y1 y2 ...
    if single:
        return [([r[0] for r in rows], [r[1 + k] for r in rows]) for k in range(len(rows[0]) - 1)]
    return [([r[2 * k] for r in rows], [r[2 * k + 1] for r in rows]) for k in range(len(rows[0]) // 2)]


def interp(xs, ys, x):
    import bisect
    i = bisect.bisect_left(xs, x)
    if i <= 0:
        return ys[0]
    if i >= len(xs):
        return ys[-1]
    x0, x1 = xs[i - 1], xs[i]
    return ys[i - 1] if x1 == x0 else ys[i - 1] + (ys[i] - ys[i - 1]) * (x - x0) / (x1 - x0)


def orbit_dev(pa, pb):
    """per vector: (start-point dev, waveform dev, edge?), smooth vectors normalised by max(p-p, RIPPLE_FLOOR)
    (the engine's own norm), edge vectors by p-p.  Edge vectors (max slope > 15 x p-p/span) are compared
    edge-tolerantly (best of small time shifts), as tests/test_pss_klu.py does."""
    A, B = read_cols(pa), read_cols(pb)
    out = []
    for (xa, ya), (xb, yb) in zip(A, B):
        xa = [t / xa[-1] for t in xa]
        xb = [t / xb[-1] for t in xb]
        pp = (max(ya) - min(ya)) or 1.0
        edge = max(abs(ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]) for i in range(len(ya) - 1) if xa[i + 1] > xa[i]) > 15 * pp
        h = 2 * max(max(b - a for a, b in zip(xa, xa[1:])), max(b - a for a, b in zip(xb, xb[1:]))) if edge else 0.0
        shifts = [h * k / 40 for k in range(-40, 41)] if edge else [0.0]
        dev = 0.0
        n = 4000
        for k in range(n + 1):
            t = xa[0] + (xa[-1] - xa[0]) * k / n
            ta = interp(xa, ya, t)
            dev = max(dev, min(abs(ta - interp(xb, yb, min(max(t + s, xb[0]), xb[-1]))) for s in shifts))
        nrm = pp if edge else max(pp, RIPPLE_FLOOR)
        out.append((abs(ya[0] - yb[0]) / nrm, dev / nrm, edge))
    return out


# ------------------------------------------------------------------ cases
FIXTURES = ["pss_r2_shoot", "pss_r2_q4_shoot", "pss_r3_shoot", "pss_r4_shoot", "pss_r6s_shoot",
            "rd_r3d_shoot", "rd_r4d_shoot", "pss_r3_gainx20_unstable", "pss_r4_d55_unstable"]
QUICK = ["pss_r2_shoot", "pss_r3_shoot", "pss_r6s_shoot", "pss_r4_d55_unstable"]
VECS = {"pss_r6s_shoot": "v(vout) v(opout) i(L1)"}   # i(Shs) is not available in plain wrdata (fixture quirk)


def mkdeck(fix, d, krylov_card=False):
    os.makedirs(d, exist_ok=True)
    s = redirect(body(fix + ".cir"), d)
    if fix in VECS:
        s = re.sub(r"(?m)^(wrdata\s+\S+).*$", lambda m: m.group(1) + " " + VECS[fix], s)
    s = set_krylov_card(s, krylov_card)
    p = os.path.join(d, "deck.cir")
    with open(p, "w") as fh:
        fh.write(s)
    for f in os.listdir(d):
        if f.endswith((".txt", ".raw")):
            os.remove(os.path.join(d, f))
    return p


def first_data_file(d):
    fs = [f for f in os.listdir(d) if f.endswith(".txt") and f != "stdout.txt"]
    if not fs:
        raise CouldNotRun(f"no wrdata output in {d}")
    return os.path.join(d, fs[0])


def verdict_class(v):
    if v.startswith("UNKNOWN"):
        return "unknown"
    return "unstable" if v.startswith("UNSTABLE") else "stable"


def check_equiv(fix, red=False):
    name = f"equiv {fix}" + (" [red-proof]" if red else "")
    try:
        dd = os.path.join(OUT, fix + ("_red" if red else ""), "dense")
        kd = os.path.join(OUT, fix + ("_red" if red else ""), "krylov")
        pd = mkdeck(fix, dd, krylov_card=False)
        pk = mkdeck(fix, kd, krylov_card=True)
        if red:   # deliberate fault: a 1 % different clock guess for the krylov deck only
            s = open(pk).read()
            s = re.sub(r"(?m)^(\s*\.?pss\s+)(\S+)", lambda m: m.group(1) + f"{float(m.group(2)) * 1.01:.6g}", s, count=1)
            open(pk, "w").write(s)
        rcd, od, wd = run(NGSPICE, pd, env={"NGSPICE_PSS_KRYLOV": "0"})
        rck, ok, wk = run(NGSPICE, pk)
        variants = []
        for vn, venv in (("eps", {"PSS_SHOOT_EPS": "2e-5"}), ("noreduce", {"PSS_SHOOT_NO_REDUCE": "1"}),
                         ("warm16", {"PSS_SHOOT_WARMUP_PERIODS": "16"})):
            if red:
                break
            vd = os.path.join(OUT, fix, "dense_" + vn)
            vp = mkdeck(fix, vd, krylov_card=False)
            _, vo, _ = run(NGSPICE, vp, env=dict({"NGSPICE_PSS_KRYLOV": "0"}, **venv))
            if parse(vo):
                variants.append((vn, vd))
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    open(os.path.join(dd, "stdout.txt"), "w").write(od)
    open(os.path.join(kd, "stdout.txt"), "w").write(ok)
    d, k = parse(od), parse(ok)
    if d is None and k is None:
        record(name, "pass", "neither run converges (both report non-convergence)")
        return
    if d is None or k is None:
        record(name, "FAIL", f"convergence differs: dense={'conv' if d else 'no'} krylov={'conv' if k else 'no'}")
        return
    try:
        devs = orbit_dev(first_data_file(dd), first_data_file(kd))
        spreads = [orbit_dev(first_data_file(dd), first_data_file(vd)) for vn, vd in variants]
        nr = [vd for vn, vd in variants if vn == "noreduce"]
        k_nr = max((w for (s_, w, e) in orbit_dev(first_data_file(nr[0]), first_data_file(kd)) if not e), default=0.0) if nr else None
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    bad = []
    worst_s = worst_w = worst_sp = 0.0
    for i, (s_dev, w_dev, edge) in enumerate(devs):
        sp = max((sv[i][1] for sv in spreads), default=0.0)
        tol = TOL_EDGE if edge else max(TOL_ORBIT, SPREAD_K * sp)
        if not edge:
            worst_s, worst_w, worst_sp = max(worst_s, s_dev), max(worst_w, w_dev), max(worst_sp, sp)
        if w_dev > tol:
            bad.append(f"{'edge ' if edge else ''}vec{i} wave {w_dev:.2e}>{tol:.2e}")
    vd, vk = verdict_class(d["verdict"]), verdict_class(k["verdict"])
    if vk == "unknown":
        vtxt = f"stability explicitly not computed ({k['verdict']})"
    elif vd != vk:
        bad.append(f"verdict dense {d['verdict']} vs krylov {k['verdict']}")
        vtxt = "verdict differs"
    else:
        drho = abs(d["rho"] - k["rho"]) / max(d["rho"], 1e-3)
        vtxt = f"{vk} rho {d['rho']:.4f}/{k['rho']:.4f}"
        if drho > TOL_RHO and abs(d["rho"] - k["rho"]) > 1e-3:
            bad.append(f"rho differs by {drho:.2e} (relative)")
    if "[arnoldi]" not in k["verdict"] and k["mode"] == "krylov":
        bad.append("krylov mode without the arnoldi marker in the verdict")
    detail = (f"orbit N start/wave {worst_s:.1e}/{worst_w:.1e} (dense variants' spread {worst_sp:.1e}; krylov vs dense-noreduce {'n/a' if k_nr is None else format(k_nr, '.1e')}); {vtxt}; period runs dense {d['runs']} "
              f"krylov {k['runs']} ({k['mode']}); iter {d['iter']}/{k['iter']}")
    record(name, "FAIL" if bad else "pass", detail + ("; " + "; ".join(bad) if bad else ""))


def check_redproof_arnoldi():
    """a wrong comparison must fail: with a 1e-9 radius tolerance the (finite-difference-limited)
    dense and Arnoldi radii cannot agree."""
    name = "redproof arnoldi radius tolerance 1e-9"
    try:
        dd = os.path.join(OUT, "pss_r3_shoot", "dense")
        kd = os.path.join(OUT, "pss_r3_shoot", "krylov")
        d = parse(open(os.path.join(dd, "stdout.txt")).read())
        k = parse(open(os.path.join(kd, "stdout.txt")).read())
    except OSError:
        record(name, "COULD-NOT-RUN", "run `equiv pss_r3_shoot` first (not in --only?)")
        return
    if d and k and abs(d["rho"] - k["rho"]) > 1e-9:
        record(name, "pass", f"comparison correctly fails at 1e-9 (|drho|={abs(d['rho'] - k['rho']):.2e})")
    else:
        record(name, "FAIL", "the radii agree to 1e-9 -- the equivalence check is vacuous")


FREERUN = re.compile(r"PSS-FREERUN: CONVERGED iter=(\d+) residual=(\S+) period=(\S+) fsw=\S+ spectral_radius=(\S+) (\S+) wall_time=(\S+)s")


def parse_fr(out):
    m = FREERUN.search(out)
    if not m:
        return None
    return dict(iter=int(m.group(1)), res=float(m.group(2)), T=float(m.group(3)), rho=float(m.group(4)), verdict=m.group(5), wall=float(m.group(6)),
                krylov="krylov: GMRES" in out, auto="chosen automatically" in out, fell="continuing with the dense monodromy" in out)


def check_freerun(red=False):
    """card #538: free-running (trigger=) shooting is matrix-free too: GMRES on the Poincare return map (M - I), products by
    pairs of trigger runs, Arnoldi stability estimate.  R5 hysteretic buck fixture and R5 + extracted RC trees."""
    name = "freerun (trigger=) matrix-free" + (" [red-proof]" if red else "")
    fix = "pss_r5_shoot"
    TOL_T = 2e-4 if not red else 1e-9   # relative period difference: the engine stops at steady_coeff 1e-4 of the ripple
    try:
        dd = os.path.join(OUT, fix, "dense")
        kd = os.path.join(OUT, fix, "krylov")
        pd = mkdeck(fix, dd, False)
        pk = mkdeck(fix, kd, True)
        _, od, _ = run(NGSPICE, pd, env={"NGSPICE_PSS_KRYLOV": "0"})
        _, ok_, _ = run(NGSPICE, pk)
        d, k = parse_fr(od), parse_fr(ok_)
        bad = []
        if not (d and k):
            record(name, "FAIL", f"convergence: dense {bool(d)} krylov {bool(k)}")
            return
        dT = abs(d["T"] - k["T"]) / d["T"]
        if dT > TOL_T:
            bad.append(f"period differs by {dT:.2e}")
        if red:
            record(name, "FAIL" if bad else "pass", f"R5 fixture period dense {d['T']:.9g} krylov {k['T']:.9g} (diff {dT:.1e}, tolerance {TOL_T:.0e})")
            return
        if not k["krylov"] or k["fell"] or not k["verdict"].endswith("[arnoldi]"):
            bad.append("not run matrix-free" + (" (fell back to dense)" if k["fell"] else ""))
        if d["krylov"] or d["auto"]:
            bad.append("NGSPICE_PSS_KRYLOV=0 run is not dense")
        if not (d["verdict"] == "stable" and k["verdict"].startswith("stable")) or abs(d["rho"] - k["rho"]) > 0.1:
            bad.append(f"stability: dense {d['verdict']} {d['rho']:.4f} vs krylov {k['verdict']} {k['rho']:.4f}")
        devs = orbit_dev(first_data_file(dd), first_data_file(kd))
        sm = max(w for (s_, w, e) in devs if not e)
        if sm > 2e-4:
            bad.append(f"orbit (smooth vectors) differs by {sm:.2e}")
        detail = (f"R5 fixture (4 states): period dense {d['T']:.9g} krylov {k['T']:.9g} (diff {dT:.1e}), orbit dev smooth {sm:.1e} (edge vectors "
                  f"{max(w for (s_, w, e) in devs if e):.1e} of p-p), rho dense {d['rho']:.4f} krylov {k['rho']:.4f}, wall {d['wall']:.1f}/{k['wall']:.1f}s")
        # PEX-loaded R5: dense vs krylov; override by card and by AUTO_MIN; automatic from 64 states
        a_d = _run_pex("freerun/pex1x2_dense", pexfx_deck(fix, 1, 2, "@D@", False))
        a_k = _run_pex("freerun/pex1x2_krylov", pexfx_deck(fix, 1, 2, "@D@", True))
        a_au = _run_pex("freerun/pex1x2_auto10", pexfx_deck(fix, 1, 2, "@D@", None), {"PSS_KRYLOV_AUTO_MIN": "10"})
        a_ov = _run_pex("freerun/pex1x2_auto10_card0", pexfx_deck(fix, 1, 2, "@D@", False), {"PSS_KRYLOV_AUTO_MIN": "10"})
        a_big = _run_pex("freerun/pex2x4_auto", pexfx_deck(fix, 2, 4, "@D@", None))
        f = {n_: parse_fr(a["o"]) for n_, a in (("d", a_d), ("k", a_k), ("au", a_au), ("ov", a_ov), ("big", a_big))}
        if not all(f.values()):
            record(name, "FAIL", detail + "; PEX R5 run(s) did not converge: " + ",".join(n_ for n_, v in f.items() if not v))
            return
        Td = f["d"]["T"]
        for n_ in ("k", "au", "ov", "big"):
            if abs(f[n_]["T"] - Td) / Td > TOL_T:
                bad.append(f"PEX R5 {n_}: period differs from dense by {abs(f[n_]['T'] - Td) / Td:.2e}")
        if not (f["k"]["krylov"] and f["au"]["krylov"] and f["au"]["auto"]):
            bad.append("PEX R5: krylov=1 / AUTO_MIN=10 did not run matrix-free (with notice)")
        if f["d"]["krylov"] or f["ov"]["krylov"] or f["ov"]["auto"] or f["d"]["auto"]:
            bad.append("PEX R5: krylov=0 override did not give the dense monodromy")
        if not (f["big"]["krylov"] and f["big"]["auto"]):
            bad.append("PEX R5 2x4 trees (96 states): automatic choice did not happen")
        for n_ in ("k", "au", "big"):
            if not f[n_]["verdict"].startswith("stable"):
                bad.append(f"PEX R5 {n_}: verdict {f[n_]['verdict']}")
        detail += (f"; PEX R5 (1x2 trees, 10 states) period dense {Td:.9g} krylov {f['k']['T']:.9g} (diff {abs(f['k']['T'] - Td) / Td:.1e}), wall dense {f['d']['wall']:.1f}s krylov {f['k']['wall']:.1f}s, "
                   f"rho krylov {f['k']['rho']:.4f} (dense {f['d']['rho']:.4f}); overrides ok; automatic (2x4 trees, 96 states): wall {f['big']['wall']:.1f}s rho {f['big']['rho']:.4f}")
        record(name, "FAIL" if bad else "pass", detail + ("; " + "; ".join(bad) if bad else ""))
    except (CouldNotRun, OSError) as e:
        record(name, "COULD-NOT-RUN", str(e))


# ------------------------------------------------------------------ card #540: free-running stability verdict
def _x0_T(out_text, err_text=None):
    """orbit start point (name -> value) and period from a PSS_FREERUN_DEBUG=1 run (x0[k] lines carry the element name)"""
    x0 = {m.group(2).lower(): float(m.group(1)) for m in re.finditer(r"DEBUG x0\[\d+\] \(\w+\) = (\S+) [CL]:(\S+)", out_text)}
    m = FREERUN.search(out_text)
    if not x0 or not m:
        raise CouldNotRun("no PSS_FREERUN_DEBUG orbit in the output")
    return x0, float(m.group(3))


def tran_multiplier(deck_text, x0, T, d, seeds=(1, 2), amp=3e-3, nper=10):
    """INDEPENDENT slowest Floquet multiplier of the free-running orbit: no shooting, no Jacobian.  The settled orbit's
    state (x0, from the engine; only the starting point) is perturbed by 0.3 % x 0.5..1.5 (random signs) in a plain
    .tran (reltol 1e-9, gear), the Poincare section is v(vout) rising through its orbit mean, and the decay per
    period of the deviation of the period (and of i(L1), v(cx) at the section) from the final cycles is measured; the
    median ratio of successive deviations over cycles 2..5 is the slowest multiplier (a heavily damped orbit: the
    deviation is gone after 6 cycles, so cycles 2..5 are the resolvable ones).  Returns (multiplier, T_reference)."""
    import random
    os.makedirs(d, exist_ok=True)
    txt = re.sub(r"(?ims)^\.control.*?^\.endc\s*$", "", deck_text)
    txt = re.sub(r"(?im)^\.(ic|save|end)\b.*$", "", txt)
    ratios, trefs = [], []
    for seed in seeds:
        rng = random.Random(seed)
        lines = []
        for l in txt.splitlines():
            t = l.split()
            if t and t[0][0].lower() in "cl" and t[0].lower() in x0 and not l.startswith("*"):
                v = x0[t[0].lower()] * (1 + amp * rng.choice([-1, 1]) * (0.5 + rng.random()))
                l = re.sub(r"(?i)\s+IC=\S+", "", l) + f" IC={v:.12g}"
            lines.append(l)
        ts = T / 4000
        o = os.path.join(d, f"o{seed}.txt")
        tr = ("\n".join(lines) + f"\n.options reltol=1e-9 method=gear\n.tran {ts:.6g} {nper * T:.6g} 0 {ts:.6g} uic\n"
              f".control\nset wr_singlescale\nset numdgt=15\nrun\nwrdata {o} v(vout) i(L1) v(cx)\nquit\n.endc\n.end\n")
        p = os.path.join(d, f"tran{seed}.cir")
        open(p, "w").write(tr)
        run(NGSPICE, p, timeout=900)
        rows = [[float(x) for x in l.split()] for l in open(o) if l.strip() and not l.lstrip().startswith(("#", "t"))]
        t = [r[0] for r in rows]; vo = [r[1] for r in rows]; il = [r[2] for r in rows]; vc = [r[3] for r in rows]
        sel = [v for tt, v in zip(t, vo) if tt > (nper - 2) * T]
        lvl = sum(sel) / len(sel)
        tc = []
        for i in range(1, len(t)):
            if vo[i - 1] < lvl <= vo[i] and t[i] > 0.3 * T:
                tc.append(t[i - 1] + (lvl - vo[i - 1]) / (vo[i] - vo[i - 1]) * (t[i] - t[i - 1]))
        per = [tc[i + 1] - tc[i] for i in range(len(tc) - 1)]
        if len(per) < 8:
            raise CouldNotRun("too few section crossings in the perturbed .tran")
        tref = (per[-1] + per[-2]) / 2
        dev = [abs(x - tref) / tref for x in per]
        ratios += [dev[k + 1] / dev[k] for k in range(1, 5) if dev[k] > 0]
        trefs.append(tref)
    ratios.sort()
    return ratios[len(ratios) // 2], sum(trefs) / len(trefs)


def check_freerun_stability(red=False):
    """card #540: the free-running Floquet verdict is right, on both engines, against an independent measurement.
    The dense engine used to print radius 24.2 UNSTABLE on the 24-state PEX-loaded R5 (true 0.037) and the matrix-free
    one 4.79 on R5 at the unscaled step; the truth is the decay per period of a perturbation of the settled orbit in a
    plain .tran (tran_multiplier).  Both engines must give stable with the radius within 25 % of it on R5 and on the
    PEX R5 (24 states).  red-proof: tolerance 1e-4 -- must fail."""
    name = "freerun stability vs independent .tran multiplier" + (" [red-proof]" if red else "")
    tolr = 1e-4 if red else 0.25
    try:
        bad, det = [], []
        for label in ("R5 fixture (4 states)", "PEX R5 2x2 (24 states)"):
            tag = "r5" if label.startswith("R5") else "pex24"
            d = os.path.join(OUT, "fr540", tag)
            if tag == "r5":
                p = mkdeck("pss_r5_shoot", os.path.join(d, "dbg"), False)
                deck = open(p).read()
            else:
                a = _run_pex("fr540/pex24_dbg", pexfx_deck("pss_r5_shoot", 2, 2, "@D@", True), {"PSS_FREERUN_DEBUG": "1"})
                deck = open(os.path.join(a["d"], "deck.cir")).read()
            if tag == "r5":
                _, o, _ = run(NGSPICE, p, env={"PSS_FREERUN_DEBUG": "1", "NGSPICE_PSS_KRYLOV": "1"})
                x0, T = _x0_T(o)
            else:
                x0, T = _x0_T(a["o"])
            truth, tref = tran_multiplier(deck, x0, T, os.path.join(d, "tran"))
            det.append(f"{label}: independent multiplier {truth:.4f} (tran period {tref:.8g}s)")
            for mode, k in (("dense", False), ("matrix-free", True)):
                if tag == "r5":
                    pp = mkdeck("pss_r5_shoot", os.path.join(d, mode), k)
                    _, o, _ = run(NGSPICE, pp, env={} if k else {"NGSPICE_PSS_KRYLOV": "0"})
                    r = parse_fr(o)
                else:
                    r = _run_pex(f"fr540/pex24_{mode}", pexfx_deck("pss_r5_shoot", 2, 2, "@D@", k))
                    r = parse_fr(r["o"])
                if not r:
                    bad.append(f"{label} {mode}: did not converge")
                    continue
                err = abs(r["rho"] - truth) / truth
                det.append(f"{mode} {r['rho']:.4f} {r['verdict']} ({r['iter']} it, {r['wall']:.1f}s)")
                if not r["verdict"].startswith("stable") or err > tolr:
                    bad.append(f"{label} {mode}: radius {r['rho']:.4f} {r['verdict']} vs independent {truth:.4f} (error {err:.0%}, tolerance {tolr:.0%})")
        record(name, "FAIL" if bad else "pass", "; ".join(det) + ("; " + "; ".join(bad) if bad else ""))
    except (CouldNotRun, OSError) as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_freerun_guard(red=False):
    """card #540: when the Floquet estimate is not trustworthy the verdict is UNKNOWN(not-computed) (radius -1), never a
    stable/UNSTABLE guess, on both engines.  The stability pass is forced back to the engine's noisy profile (reltol 1e-5,
    steps 1e-6): the step-halving test must disagree and refuse.  Control: the normal profile gives stable.
    red-proof: the check demands a verdict from the forced-noise run -- must fail."""
    name = "freerun stability guard (forced FD noise -> UNKNOWN(not-computed))" + (" [red-proof]" if red else "")
    noisy = {"PSS_FREERUN_STAB_RELTOL": "1e-5", "PSS_FREERUN_STAB_EPS": "1e-6", "PSS_FREERUN_STAB_RIPPLE": "1e-5"}
    try:
        bad, det = [], []
        for mode, k in (("dense", False), ("matrix-free", True)):
            ctl = parse_fr(_run_pex(f"fr540/guard_ctl_{mode}", pexfx_deck("pss_r5_shoot", 1, 2, "@D@", k))["o"])
            a = _run_pex(f"fr540/guard_noisy_{mode}", pexfx_deck("pss_r5_shoot", 1, 2, "@D@", k), noisy)
            r = parse_fr(a["o"])
            if not ctl or not r:
                bad.append(f"{mode}: did not converge")
                continue
            det.append(f"{mode}: control {ctl['rho']:.4f} {ctl['verdict']}; forced noise {r['rho']:.4f} {r['verdict']}")
            if not ctl["verdict"].startswith("stable"):
                bad.append(f"{mode}: control run {ctl['verdict']}")
            if red:
                if not r["verdict"].startswith("stable"):
                    bad.append(f"{mode}: forced-noise run gave no stable verdict ({r['verdict']})")
            else:
                if not (r["verdict"] == "UNKNOWN(not-computed)" and r["rho"] == -1.0 and "stability not computed" in a["o"]):
                    bad.append(f"{mode}: forced noise gave {r['verdict']} radius {r['rho']} (UNKNOWN(not-computed), radius -1 and a reason line expected)")
        record(name, "FAIL" if bad else "pass", "; ".join(det) + ("; " + "; ".join(bad) if bad else ""))
    except (CouldNotRun, OSError) as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_freerun_routing():
    """card #540: free-running shooting is matrix-free by default at any size (notice line, [arnoldi] verdict), also with
    .pac on the orbit; krylov=0 / NGSPICE_PSS_KRYLOV=0 is the dense monodromy (no notice, eigenvalue list); the dense Newton
    converges on the 24-state PEX R5 in <= 12 iterations (17 before, with the 1e-6 step) and finds the matrix-free period."""
    name = "freerun default routing, dense Newton step"
    try:
        bad = []
        d = os.path.join(OUT, "fr540", "route")
        pa = mkdeck("pss_r5_shoot", os.path.join(d, "auto"), False)
        pd = mkdeck("pss_r5_shoot", os.path.join(d, "dense"), False)
        _, oa, _ = run(NGSPICE, pa)
        _, od, _ = run(NGSPICE, pd, env={"NGSPICE_PSS_KRYLOV": "0"})
        a, dn = parse_fr(oa), parse_fr(od)
        if not (a and dn):
            raise CouldNotRun("R5 did not converge")
        if not (a["krylov"] and a["auto"] and a["verdict"].endswith("[arnoldi]")):
            bad.append("R5 fixture (4 states) without a krylov key did not run matrix-free automatically")
        if dn["krylov"] or dn["auto"] or "FLOQUET multiplier[" not in od:
            bad.append("NGSPICE_PSS_KRYLOV=0 did not give the dense monodromy report")
        pk = pexfx_deck("pss_r5_shoot", 2, 2, "@D@", True)
        k24 = parse_fr(_run_pex("fr540/route_k24", pk)["o"])
        d24 = _run_pex("fr540/route_d24", pexfx_deck("pss_r5_shoot", 2, 2, "@D@", False))
        d24p = parse_fr(d24["o"])
        if not (k24 and d24p):
            raise CouldNotRun("PEX R5 2x2 did not converge")
        if d24p["iter"] > 12:
            bad.append(f"dense Newton on the 24-state PEX R5 took {d24p['iter']} iterations (<= 12 expected)")
        if abs(d24p["T"] - k24["T"]) / k24["T"] > 2e-4:
            bad.append(f"dense period {d24p['T']:.9g} vs matrix-free {k24['T']:.9g}")
        pac = os.path.join(d, "pac")
        os.makedirs(pac, exist_ok=True)
        pp = os.path.join(pac, "deck.cir")
        open(pp, "w").write(redirect(body("pac_r5_ctrl2out_20k.cir"), pac))
        _, op, _ = run(NGSPICE, pp)
        if "PAC-TD: done" not in op or "chosen automatically" not in op:
            bad.append("free-running .pac (20 kHz) did not run after the automatic matrix-free PSS")
        else:
            parts = open(os.path.join(pac, "r5_ctrl2out_20k.txt")).read().split()
            db = 20 * math.log10(abs(complex(float(parts[1]), float(parts[2]))))
            if abs(db - 0.0411) > 0.005:   # the dense-path value (tests/test_pac.py, card #503 THIRD pass reference 0.0411 dB)
                bad.append(f".pac after the matrix-free PSS: {db:.4f} dB vs 0.0411 dB (dense path)")
        record(name, "FAIL" if bad else "pass",
               f"R5 auto matrix-free {a['verdict']} rho {a['rho']:.4f}; dense {dn['verdict']} rho {dn['rho']:.4f}; PEX R5 24 states dense {d24p['iter']} it {d24p['wall']:.1f}s "
               f"vs matrix-free {k24['iter']} it {k24['wall']:.1f}s, period diff {abs(d24p['T'] - k24['T']) / k24['T']:.1e}" + ("; " + "; ".join(bad) if bad else ""))
    except (CouldNotRun, OSError) as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_card_env():
    name = "card key vs environment"
    fix = "pss_r2_shoot"
    try:
        d1 = os.path.join(OUT, fix, "ce_card_envoff")
        d2 = os.path.join(OUT, fix, "ce_env_on")
        p1 = mkdeck(fix, d1, True)
        p2 = mkdeck(fix, d2, False)
        _, o1, _ = run(NGSPICE, p1, env={"NGSPICE_PSS_KRYLOV": "0"})
        _, o2, _ = run(NGSPICE, p2, env={"NGSPICE_PSS_KRYLOV": "1"})
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    a, b = parse(o1), parse(o2)
    ok = a and b and a["mode"] == "dense" and b["mode"] == "krylov"
    record(name, "pass" if ok else "FAIL",
           f"card krylov=1 + env 0 -> {a and a['mode']}; no card + env 1 -> {b and b['mode']}")


def _strip_ac(b):
    return re.sub(r"\s+AC\s+\S+", "", b)


R3P = "500e3 300e-6 vout 200 10 150 1e-6"
XFER = {   # name -> (fixture, card, vectors, strip AC, extra env)
    "pac r3 ctrl2out": ("pac_r3_ctrl2out", f'.pac {R3P} dec 5 500 200000 method="shoot" solver="td"', "v(vout)", False, {}),
    "pstb r3": ("pac_r3_ctrl2out", f'.pstb {R3P} dec 5 500 200000 probe="Vprobe" method="shoot"', "t", True, {}),
    "pxf r3 adjoint=0": ("pac_r3_ctrl2out", f'.pxf {R3P} vout dec 5 500 200000 method="shoot" solver="td" '
                         f'sources="Vin" sources="Vctrl" adjoint=0', "pxf_vin pxf_vctrl", False, {"PAC_EXACT_JACOBIAN": "1"}),
    "pnoise r2 smooth": ("pnoise_r2_smooth", '.pnoise 100e3 100e-6 vx 2000 1 50 1e-6 vx Vin dec 5 1k 1meg '
                         'method="shoot" solver="td" dev=1', "onoise_spectrum", False, {}),
}


def check_xfer(label):
    name = f"{label} after krylov PSS"
    fix, card, vecs, strip, xenv = XFER[label]
    safe = label.replace(" ", "_").replace("=", "")
    try:
        res = {}
        for mode in ("dense", "krylov"):
            d = os.path.join(OUT, safe, mode)
            os.makedirs(d, exist_ok=True)
            b = body(fix + ".cir")
            b = "\n".join(l for l in b.splitlines() if l.strip().lower() != ".end")
            b = re.sub(r"(?ms)^\.control.*?^\.endc\s*$", "", b)
            if strip:
                b = _strip_ac(b)
            deck = (b + "\n" + card + (" krylov=1" if mode == "krylov" else "") +
                    "\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata " +
                    os.path.join(d, "out.txt") + " " + vecs + "\nquit\n.endc\n.end\n")
            p = os.path.join(d, "deck.cir")
            if os.path.exists(os.path.join(d, "out.txt")):
                os.remove(os.path.join(d, "out.txt"))
            open(p, "w").write(deck)
            env = dict(xenv)
            env["NGSPICE_PSS_KRYLOV"] = "0" if mode == "dense" else "1"
            rc, o, w = run(NGSPICE, p, env=env)
            open(os.path.join(d, "stdout.txt"), "w").write(o)
            res[mode] = (rc, o, os.path.join(d, "out.txt"))
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return

    def rd(p):
        rows = []
        for l in open(p).read().splitlines()[1:]:
            try:
                rows.append([float(v) for v in l.split()])
            except ValueError:
                pass
        return rows
    try:
        ra, rb = rd(res["dense"][2]), rd(res["krylov"][2])
    except OSError as e:
        record(name, "FAIL", f"no output: {e}")
        return
    if len(ra) != len(rb) or not ra:
        record(name, "FAIL", f"sweep shapes differ ({len(ra)} vs {len(rb)}); krylov tail: {res['krylov'][1][-300:]}")
        return
    ncol = len(ra[0]) - 1
    mx = max(math.sqrt(sum(v * v for v in r[1:])) for r in ra)
    dev = max(math.sqrt(sum((a - b) ** 2 for a, b in zip(r1[1:], r2[1:]))) for r1, r2 in zip(ra, rb)) / mx
    k = parse(res["krylov"][1])
    ok = k is not None and k["mode"] == "krylov" and dev < 1e-3
    record(name, "pass" if ok else "FAIL",
           f"{len(ra)} points x {ncol} columns, max|dH|/max|H| = {dev:.2e} (tol 1e-3); krylov PSS mode: {k and k['mode']}")


def check_identical(red=False):
    name = "krylov OFF byte-identical to baseline" + (" [red-proof]" if red else "")
    if not BASE:
        record(name, "COULD-NOT-RUN", "NGSPICE_BASE not set")
        return
    bad, n = [], 0
    for fix in (["pss_r3_shoot"] if red else FIXTURES + ["pss_r5_shoot", "pss_r4_d45_shoot", "pss_r3_points1000_shoot"]):
        try:
            outs = []
            for tag, binary in (("new", NGSPICE), ("base", BASE)):
                d = os.path.join(OUT, fix, "ident_" + tag)
                p = mkdeck(fix, d, False)
                env = dict(os.environ, OMP_NUM_THREADS="2")
                env.pop("PSS_SHOOT_COUNT_RUNS", None)
                env.pop("NGSPICE_PSS_KRYLOV", None)
                if red and tag == "new":   # deliberate difference: a different FD step in the new build only
                    env["PSS_SHOOT_EPS"] = "2e-5"
                t0 = time.time()
                pr = subprocess.run([binary, "-b", p], capture_output=True, text=True, cwd=d, env=env, timeout=1800)
                txt = re.sub(r"wall_time=[\d.]+s", "wall_time=X", pr.stdout + pr.stderr)
                txt = re.sub(r"(?m)^.*Reference value.*$", "", txt)   # a timing line of the front end (value and presence vary run to run)
                txt = re.sub(r"\n\s*\n", "\n", txt)               # ... and so does the blank line around it
                files = {}
                for f in sorted(os.listdir(d)):
                    if f.endswith((".txt", ".raw")) and f != "stdout.txt":
                        raw = open(os.path.join(d, f), "rb").read()
                        # the .raw header carries the build date and run date: mask those two lines
                        raw = re.sub(rb"(?m)^(Command|Date):.*$", b"", raw)
                        files[f] = raw
                outs.append((txt.replace(d, "D"), files))
        except (CouldNotRun, subprocess.TimeoutExpired, OSError) as e:
            bad.append(f"{fix}: could not run ({e})")
            continue
        n += 1
        (ta, fa), (tb, fb) = outs
        if ta != tb:
            bad.append(f"{fix}: stdout differs")
        elif fa != fb:
            bad.append(f"{fix}: wrdata/raw files differ")
    record(name, "FAIL" if bad else "pass", f"{n} fixtures compared (stdout with wall_time masked, wrdata/raw files)" +
           ("; " + "; ".join(bad) if bad else ""))


def check_asan_pex():
    name = "asan krylov=1 pex b2d3 (73 states)"
    try:
        d = os.path.join(OUT, "pex_b2d3", "asan")
        os.makedirs(d, exist_ok=True)
        deck, nelem = pex_deck(2, 3, d, True)
        p = os.path.join(d, "deck.cir")
        open(p, "w").write(deck)
        rc, o, w = run(ASAN, p, env={"ASAN_OPTIONS": "detect_leaks=0"})
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    k = parse(o)
    bad = rc != 0 or "AddressSanitizer" in o or k is None or k["mode"] != "krylov"
    record(name, "FAIL" if bad else "pass", f"rc={rc}, converged={'yes' if k else 'no'}, mode={k and k['mode']}, wall {w:.0f}s")


def check_asan():
    if not ASAN:
        record("asan", "COULD-NOT-RUN", "NGSPICE_ASAN not set")
        return
    for fix in ("pss_r3_shoot", "rd_r3d_shoot", "pss_r3_gainx20_unstable"):   # last one exercises the dense fallback
        name = f"asan krylov=1 {fix}"
        try:
            d = os.path.join(OUT, fix, "asan")
            p = mkdeck(fix, d, True)
            rc, o, w = run(ASAN, p, env={"ASAN_OPTIONS": "detect_leaks=0"})
        except CouldNotRun as e:
            record(name, "COULD-NOT-RUN", str(e))
            continue
        k = parse(o)
        bad = rc != 0 or "AddressSanitizer" in o or k is None or not k["mode"].startswith("krylov")
        record(name, "FAIL" if bad else "pass", f"rc={rc}, converged={'yes' if k else 'no'}, mode={k and k['mode']}, wall {w:.0f}s")
    # card #538: free-running matrix-free (R5 fixture)
    name = "asan krylov=1 pss_r5_shoot (free-running)"
    try:
        d = os.path.join(OUT, "pss_r5_shoot", "asan")
        p = mkdeck("pss_r5_shoot", d, True)
        rc, o, w = run(ASAN, p, env={"ASAN_OPTIONS": "detect_leaks=0"})
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    k = parse_fr(o)
    bad = rc != 0 or "AddressSanitizer" in o or k is None or not k["krylov"] or k["fell"]
    record(name, "FAIL" if bad else "pass", f"rc={rc}, converged={'yes' if k else 'no'}, matrix-free={k and k['krylov']}, wall {w:.0f}s")


# ------------------------------------------------------------------ PEX-size
def pex_deck(branch, depth, d, krylov):
    sys.path.insert(0, os.path.join(ROOT, "bench"))
    import gen_pex_buck as g
    b = body("pac_r3_ctrl2out.cir")
    b = "\n".join(l for l in b.splitlines() if l.strip().lower() != ".end")
    b = re.sub(r"(?ms)^\.control.*?^\.endc\s*$", "", b)
    b = re.sub(r"(?m)^Cout cx 0 \{cval\} IC=0", "Cout cx gret {cval} IC=0\nRgret gret 0 5m", b)
    lines, nelem = [], 0
    for nm, anchor, rc in (("A", "vout", "10"), ("B", "sw", "0.5"), ("G", "gret", "1")):
        tl, root, far, nn = g.tree_lines(branch, depth, 0.2, 1.6, 20e-15, 0.6, 3, 40, 2e-15)
        lines.append(f"Rcouple_{nm} {anchor} {root.replace('tree_', 'tree' + nm + '_')} {rc}\n")
        for l in tl:
            if l.startswith("*"):
                continue
            t = l.split(None, 1)
            lines.append(t[0] + "_" + nm + " " + t[1].replace("tree_", "tree" + nm + "_"))
        nelem += len([1 for l in tl if not l.startswith("*")])
    pex = "* PEX-size parasitic nets\n" + "".join(lines) + ".save v(vout) v(sw) v(opout) i(L1)\n"
    card = "pss 500e3 300e-6 vout 200 10 150 1e-6 uic method=shoot" + {True: " krylov=1", False: " krylov=0", None: ""}[krylov]
    deck = (b + "\n" + pex + "\n.options klu\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n" + card +
            "\nwrdata " + os.path.join(d, "out.txt") + " v(vout) v(sw) v(opout) i(L1)\nquit\n.endc\n.end\n")
    return deck, nelem


MEM_BOUND_MB = 150.0


def _rss_run(tag, deck_text, env):
    """one ngspice run in its own directory; returns (maxrss_MB, wall, stdout+stderr, out.txt path)"""
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):
        if f.endswith((".txt", ".raw")):
            os.remove(os.path.join(d, f))
    p = os.path.join(d, "deck.cir")
    open(p, "w").write(deck_text.replace("@D@", d))
    e = dict(os.environ, OMP_NUM_THREADS="2", PSS_SHOOT_COUNT_RUNS="1")
    for k in ("NGSPICE_PSS_KRYLOV", "PSS_KEEP_RUN_DATA"):
        e.pop(k, None)
    e.update(env)
    t0 = time.time()
    try:
        with open(os.path.join(d, "stdout.txt"), "w") as fo:
            pr = subprocess.Popen([NGSPICE, "-b", p], stdout=fo, stderr=subprocess.STDOUT, cwd=d, env=e)
            _, status, ru = os.wait4(pr.pid, 0)
    except OSError as ex:
        raise CouldNotRun(f"cannot start {NGSPICE}: {ex}")
    o = open(os.path.join(d, "stdout.txt")).read()
    return ru.ru_maxrss / 1024.0, time.time() - t0, o, os.path.join(d, "out.txt")


def check_memory(red=False):
    name = "memory: peak RSS of a PEX-size .pss (no .save)" + (" [red-proof]" if red else "")
    try:
        deck, nelem = pex_deck(4, 4, "@D@", True)
        deck = "\n".join(l for l in deck.splitlines() if not l.startswith(".save"))
        keep = {"PSS_KEEP_RUN_DATA": "1"}
        rss, wall, o, out = _rss_run("mem553/new", deck, keep if red else {})
        k = parse(o)
        st = re.search(r"(\d+) state\(s\)", o)
        if k is None:
            record(name, "FAIL", "the run did not converge")
            return
        bad = []
        if rss > MEM_BOUND_MB:
            bad.append(f"peak RSS {rss:.0f} MB > {MEM_BOUND_MB:.0f} MB")
        detail = f"{st and st.group(1)} states: peak RSS {rss:.0f} MB (bound {MEM_BOUND_MB:.0f}), wall {wall:.1f}s, {k['runs']} period runs"
        if not red:
            rss_o, wall_o, o_o, out_o = _rss_run("mem553/keep", deck, keep)
            ko = parse(o_o)
            same = ko is not None and ko["runs"] == k["runs"] and open(out, "rb").read() == open(out_o, "rb").read()
            detail += f"; PSS_KEEP_RUN_DATA=1: {rss_o:.0f} MB; orbit {'bit-identical' if same else 'DIFFERS'}, same period-run count {ko and ko['runs'] == k['runs']}"
            if not same:
                bad.append("orbit or period-run count differs from the retaining run")
            if rss_o < 3 * rss:
                bad.append(f"the retaining run ({rss_o:.0f} MB) is not 3x the new one: the check cannot see the fix")
        record(name, "FAIL" if bad else "pass", detail + ("; " + "; ".join(bad) if bad else ""))
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pex(branch, depth, alternations=2):
    name = f"pex b{branch}d{depth}"
    runs = {"dense": [], "krylov": []}
    states = None
    try:
        for r in range(alternations):
            order = ("dense", "krylov") if r % 2 == 0 else ("krylov", "dense")
            for mode in order:
                d = os.path.join(OUT, f"pex_b{branch}d{depth}", mode + str(r))
                os.makedirs(d, exist_ok=True)
                deck, nelem = pex_deck(branch, depth, d, mode == "krylov")
                p = os.path.join(d, "deck.cir")
                open(p, "w").write(deck)
                load = open("/proc/loadavg").read().split()[0]
                cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", NGSPICE]
                t0 = time.time()
                e = dict(os.environ, OMP_NUM_THREADS="2", PSS_SHOOT_COUNT_RUNS="1")
                e.pop("NGSPICE_PSS_KRYLOV", None)
                pr = subprocess.run(cmd[:-1] + [NGSPICE, "-b", p], capture_output=True, text=True, cwd=d, env=e, timeout=3600)
                wall = time.time() - t0
                o = pr.stdout + pr.stderr
                open(os.path.join(d, "stdout.txt"), "w").write(o)
                k = parse(o)
                st = re.search(r"(\d+) state\(s\)", o)
                states = st and int(st.group(1))
                runs[mode].append(dict(k=k, wall=wall, load=load, out=os.path.join(d, "out.txt"), rc=pr.returncode))
                print(f"   {name} {mode} rep{r}: wall {wall:.1f}s load {load} {k and (k['runs'], k['rho'], k['verdict'])}", flush=True)
    except (CouldNotRun, subprocess.TimeoutExpired) as e:
        record(name, "COULD-NOT-RUN", str(e))
        return
    try:
        bad = []
        if any(r["k"] is None for m in runs.values() for r in m):
            record(name, "FAIL", "a run did not converge")
            return
        dev = orbit_dev(runs["dense"][0]["out"], runs["krylov"][0]["out"])
        sm = max(w for (s, w, e) in dev if not e)
        ss = max(s for (s, w, e) in dev if not e)
        if sm > TOL_ORBIT:
            bad.append(f"orbit dev {sm:.2e}")
        dr = runs["dense"][0]["k"]["runs"]
        kr = runs["krylov"][0]["k"]["runs"]
        if dr < 5 * kr:
            bad.append(f"period runs only {dr / kr:.1f}x fewer")
        dw = sorted(r["wall"] for r in runs["dense"])[len(runs["dense"]) // 2]
        kw = sorted(r["wall"] for r in runs["krylov"])[len(runs["krylov"]) // 2]
        record(name, "FAIL" if bad else "pass",
               f"{states} states: period runs dense {dr} krylov {kr} ({dr / kr:.1f}x); wall median dense {dw:.1f}s krylov {kw:.1f}s "
               f"({dw / kw:.1f}x); orbit start/wave dev {ss:.1e}/{sm:.1e} of p-p; rho dense {runs['dense'][0]['k']['rho']:.4f} "
               f"krylov {runs['krylov'][0]['k']['rho']:.4f} ({runs['krylov'][0]['k']['verdict']})" +
               ("; " + "; ".join(bad) if bad else ""))
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))



# ------------------------------------------------------------------ card #538: automatic choice, warm-up schedule
def pexfx_deck(fix, branch, depth, d, krylov):
    """a fixture deck (R4 / R3d / R4d ...) with two extracted RC trees on vout (10 ohm) and sw (0.5 ohm), .options klu"""
    sys.path.insert(0, os.path.join(ROOT, "bench"))
    import gen_pex_buck as g
    b = body(fix + ".cir")
    card = re.search(r"(?m)^pss .*$", b).group(0)
    b = re.sub(r"(?ms)^\.control.*?^\.endc\s*$", "", b)
    b = "\n".join(l for l in b.splitlines() if l.strip().lower() != ".end")
    lines = []
    for nm, anchor, rc in (("A", "vout", "10"), ("B", "sw", "0.5")):
        tl, root, far, nn = g.tree_lines(branch, depth, 0.2, 1.6, 20e-15, 0.6, 3, 40, 2e-15)
        lines.append(f"Rcouple_{nm} {anchor} {root.replace('tree_', 'tree' + nm + '_')} {rc}\n")
        for l in tl:
            if l.startswith("*"):
                continue
            t = l.split(None, 1)
            lines.append(t[0] + "_" + nm + " " + t[1].replace("tree_", "tree" + nm + "_"))
    card += {True: " krylov=1", False: " krylov=0", None: ""}[krylov]
    return (b + "\n* PEX trees\n" + "".join(lines) + ".save v(vout) v(sw) i(L1)\n.options klu\n.control\nset wr_singlescale\n"
            "set wr_vecnames\nset numdgt=15\n" + card + "\nwrdata " + os.path.join(d, "out.txt") + " v(vout) v(sw) i(L1)\nquit\n.endc\n.end\n")


def _run_pex(tag, deck, env=None):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):
        if f.endswith((".txt", ".raw")):
            os.remove(os.path.join(d, f))
    p = os.path.join(d, "deck.cir")
    open(p, "w").write(deck.replace("@D@", d))
    e = dict(os.environ, OMP_NUM_THREADS="2", PSS_SHOOT_COUNT_RUNS="1")
    for k in ("NGSPICE_PSS_KRYLOV", "PSS_KRYLOV_AUTO_MIN", "PSS_SHOOT_WARMUP_PERIODS", "PSS_KRYLOV_WARMUP"):
        e.pop(k, None)
    e.update(env or {})
    t0 = time.time()
    try:
        pr = subprocess.run([NGSPICE, "-b", p], capture_output=True, text=True, cwd=d, env=e, timeout=1800)
    except subprocess.TimeoutExpired:
        raise CouldNotRun(f"{tag} timed out")
    o = pr.stdout + pr.stderr
    open(os.path.join(d, "stdout.txt"), "w").write(o)
    st = re.search(r"(\d+) state\(s\)", o)
    return dict(d=d, out=os.path.join(d, "out.txt"), o=o, k=parse(o), wall=time.time() - t0, states=st and int(st.group(1)),
                auto="chosen automatically" in o, rc=pr.returncode)


def _pex_r3(krylov, d_tag, b=2, dp=3):
    deck, nelem = pex_deck(b, dp, "@D@", krylov)
    return deck


def check_auto(red=False):
    name = "auto threshold picks krylov on a >=64-state deck" + (" [red-proof]" if red else "")
    try:
        r = _run_pex("auto/b2d3_default" + ("_red" if red else ""), _pex_r3(None, ""), {"PSS_KRYLOV_AUTO_MIN": "1000"} if red else None)
        if red:   # the same assertion with the automatic choice disabled by a high threshold: must fail
            bad = [] if (r["k"] and r["k"]["mode"] == "krylov" and r["auto"]) else ["automatic choice did not happen"]
            record(name, "FAIL" if bad else "pass", f"{r['states']} states, mode {r['k'] and r['k']['mode']}, notice {r['auto']}; " + "; ".join(bad))
            return
        bad = []
        if r["states"] is None or r["states"] < 64:
            bad.append(f"deck has {r['states']} states (<64: not a test of the threshold)")
        if not (r["k"] and r["k"]["mode"] == "krylov" and "[arnoldi]" in r["k"]["verdict"]):
            bad.append("not run matrix-free")
        if not r["auto"]:
            bad.append("the one-line notice is missing")
        if len(re.findall(r"chosen automatically", r["o"])) != 1:
            bad.append("notice not exactly once")
        arms = {}
        for tag, deck, env in (("card0", _pex_r3(False, ""), None),
                               ("env0", _pex_r3(None, ""), {"NGSPICE_PSS_KRYLOV": "0"}),
                               ("min0", _pex_r3(None, ""), {"PSS_KRYLOV_AUTO_MIN": "0"}),
                               ("min73", _pex_r3(None, ""), {"PSS_KRYLOV_AUTO_MIN": "73"}),
                               ("min74", _pex_r3(None, ""), {"PSS_KRYLOV_AUTO_MIN": "74"}),
                               ("card1", _pex_r3(True, ""), None)):
            arms[tag] = _run_pex("auto/b2d3_" + tag, deck, env)
        for tag, want_mode, want_notice in (("card0", "dense", False), ("env0", "dense", False), ("min0", "dense", False),
                                            ("min73", "krylov", True), ("min74", "dense", False), ("card1", "krylov", False)):
            a = arms[tag]
            if not (a["k"] and a["k"]["mode"] == want_mode):
                bad.append(f"{tag}: mode {a['k'] and a['k']['mode']}, want {want_mode}")
            if a["auto"] != want_notice:
                bad.append(f"{tag}: notice {'present' if a['auto'] else 'absent'}, want {'present' if want_notice else 'absent'}")
        dev = max(w for (s_, w, e) in orbit_dev(arms["card0"]["out"], r["out"]) if not e)
        if dev > TOL_ORBIT:
            bad.append(f"automatic orbit differs from dense by {dev:.2e}")
        # a small deck: no card, stays dense, no notice
        d7 = os.path.join(OUT, "auto", "r3_fixture")
        p7 = mkdeck("pss_r3_shoot", d7, False)
        rc7, o7, w7 = run(NGSPICE, p7)
        k7 = parse(o7)
        if not (k7 and k7["mode"] == "dense") or "chosen automatically" in o7:
            bad.append(f"7-state fixture: mode {k7 and k7['mode']}, notice {'chosen automatically' in o7}")
        record(name, "FAIL" if bad else "pass",
               f"{r['states']} states: auto -> {r['k']['mode']} ({r['k']['runs']} runs) vs dense {arms['card0']['k']['runs']} runs, orbit dev {dev:.1e}; "
               f"overrides krylov=0/env 0/AUTO_MIN=0/74 -> dense, AUTO_MIN=73 -> krylov; 7-state fixture {k7 and k7['mode']}" + ("; " + "; ".join(bad) if bad else ""))
    except (CouldNotRun, OSError) as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_warmup(red=False):
    name = "matrix-free warm-up: no Newton attempt on the first 130 periods, same orbit" + (" [red-proof]" if red else "")
    # name, deck maker, tolerance vs dense (x TOL_ORBIT)
    cases = [("R3 PEX 145 states", lambda kr: pex_deck(2, 4, "@D@", kr)[0], 1.0, False),
             ("R4 PEX 99 states (rho 0.994)", lambda kr: pexfx_deck("pss_r4_shoot", 2, 4, "@D@", kr), 10.0, True),
             ("R3d PEX 101 states (dead time)", lambda kr: pexfx_deck("rd_r3d_shoot", 2, 4, "@D@", kr), 10.0, False)]
    if red:
        cases = cases[:1]
    try:
        for cname, mk, tolx, slow in cases:
            tag = "warm/" + re.sub(r"\W+", "_", cname) + ("_red" if red else "")
            dense = _run_pex(tag + "_dense", mk(False))
            deck_new = mk(None)
            if red:   # deliberate fault: the krylov deck runs at a 1 % different clock guess
                deck_new = re.sub(r"(?m)^(pss\s+)(\S+)", lambda m: m.group(1) + f"{float(m.group(2)) * 1.01:.6g}", deck_new, count=1)
            new = _run_pex(tag + "_new", deck_new)
            if red:
                dev = max(w for (s_, w, e) in orbit_dev(dense["out"], new["out"]) if not e)
                record(name, "FAIL" if dev > TOL_ORBIT else "pass", f"{cname}: orbit dev {dev:.2e} (tolerance {TOL_ORBIT:.0e})")
                return
            old = _run_pex(tag + "_old", mk(None), {"PSS_KRYLOV_WARMUP": "0"})          # attempts at every stage (card #537 schedule)
            bad1 = _run_pex(tag + "_bad1", mk(None), {"PSS_SHOOT_WARMUP_PERIODS": "1"})  # deliberately bad start, explicit schedule
            bad = []
            for lab, a in (("dense", dense), ("new", new), ("old", old), ("bad start", bad1)):
                if not a["k"]:
                    bad.append(f"{lab}: did not converge")
            if bad:
                record(name, "FAIL", f"{cname}: " + "; ".join(bad))
                continue
            devs = {}
            for lab, a in (("new", new), ("old", old)):
                devs[lab] = max(w for (s_, w, e) in orbit_dev(dense["out"], a["out"]) if not e)
                if devs[lab] > TOL_ORBIT * tolx:
                    bad.append(f"{lab} orbit differs from dense by {devs[lab]:.2e}")
            # the bad start only has to land on the same orbit within the engine's own spread (grid dependence of the period map)
            devs["bad start"] = max(w for (s_, w, e) in orbit_dev(dense["out"], bad1["out"]) if not e)
            if devs["bad start"] > 5e-4:
                bad.append(f"bad start orbit differs from dense by {devs['bad start']:.2e}")
            dno = max(w for (s_, w, e) in orbit_dev(old["out"], new["out"]) if not e)
            # same warm-up endpoint: identical to 1e-6 on R3-class decks; on the slow-mode decks (rho 0.994) a failed
            # Newton attempt leaves noise-level (reltol) differences in the next transient that 1/(1-rho) = 166 amplifies
            if dno > (1e-6 if tolx == 1.0 else TOL_ORBIT * tolx):
                bad.append(f"new vs old schedule orbit differs by {dno:.2e} (same warm-up endpoint expected)")
            if new["k"]["mode"] != "krylov":
                bad.append(f"new schedule ran {new['k']['mode']}")
            if slow:
                # card #550: on the rho 0.994 deck the matrix-free final stage plateaus a few x tol above the target (period-map noise
                # amplified 1/(1-rho)); whether it steps across depends on the last bits of exp/log. Development workstation (AMD Ryzen AI Max+ 395, Linux): 126 runs (old schedule 133-146, dense 273);
                # second workstation (other libm): the same stage stalled at 4.0e-6 and the run fell back to the dense path (405 runs); with the matrix-free replay
                # rung (card #550, branch aug-550) the stall costs one more matrix-free attempt (both workstations over 8 tiny perturbations of eta/maxdim/warm-up: 111-228 runs,
                # old-schedule 133-146). So the saving is not guaranteed here: the claims are "stays matrix-free" (checked above) and "never costlier
                # than the dense engine" (bound = measured worst 228 x 1.2).
                if new["k"]["runs"] > dense["k"]["runs"]:
                    bad.append(f"slow-mode deck: {new['k']['runs']} period runs, more than the dense engine's {dense['k']['runs']}")
            elif not (new["k"]["runs"] < old["k"]["runs"] or new["k"]["runs"] == old["k"]["runs"] and new["k"]["warm"] != 0 and "no Newton attempt" not in new["o"]):
                bad.append(f"no saving: {new['k']['runs']} vs {old['k']['runs']} period runs")
            if "no Newton attempt" not in new["o"] or "no Newton attempt" in old["o"]:
                bad.append("skip notice missing in the default run or present with PSS_KRYLOV_WARMUP=0")
            drho = abs(new["k"]["rho"] - dense["k"]["rho"]) / max(dense["k"]["rho"], 1e-3)
            if drho > TOL_RHO:
                bad.append(f"rho differs by {drho:.2e}")
            record(name, "FAIL" if bad else "pass",
                   f"{cname}: period runs dense {dense['k']['runs']}, new {new['k']['runs']}, attempts at every stage {old['k']['runs']} ({1 - new['k']['runs'] / old['k']['runs']:.0%} fewer), "
                   f"bad start (1 period) {bad1['k']['runs']} ({bad1['k']['mode']}); orbit dev vs dense new/old/bad {devs['new']:.1e}/{devs['old']:.1e}/{devs['bad start']:.1e}, "
                   f"new vs old {dno:.1e}; rho {dense['k']['rho']:.4f}/{new['k']['rho']:.4f}" + ("; " + "; ".join(bad) if bad else ""))
    except (CouldNotRun, OSError) as e:
        record(name, "COULD-NOT-RUN", str(e))


def main():
    only = args[args.index("--only") + 1] if "--only" in args else None
    want = lambda n: only is None or only in n
    fixtures = QUICK if "--quick" in args else FIXTURES
    for fx in fixtures:
        if want("equiv " + fx):
            check_equiv(fx)
    if want("freerun"):
        check_freerun()
    if want("freerun") or want("fr540"):
        check_freerun_stability()
        check_freerun_guard()
        check_freerun_routing()
    if want("card"):
        check_card_env()
    if "--quick" not in args:
        for lab in XFER:
            if want(lab):
                check_xfer(lab)
    if BASE and want("identical"):
        check_identical()
    if ASAN and want("asan"):
        check_asan()
        check_asan_pex()
    if want("memory") and "--quick" not in args:
        check_memory()
    if want("auto"):
        check_auto()
    if want("warm") and "--quick" not in args:
        check_warmup()
    if "--redproof" in args:
        check_auto(red=True)
        check_memory(red=True)
        check_warmup(red=True)
        check_freerun(red=True)
        check_freerun_stability(red=True)
        check_freerun_guard(red=True)
        check_equiv("pss_r3_shoot", red=True)
        check_redproof_arnoldi()
        if BASE:
            check_identical(red=True)
    if "--pex" in args:
        check_pex(3, 4)
        if "--big" in args:
            check_pex(4, 4)
    fails = [r for r in results if r[1] == "FAIL"]
    cnr = [r for r in results if r[1] == "COULD-NOT-RUN"]
    npass = len([r for r in results if r[1] == "pass"])
    print(f"\n{npass} pass, {len(fails)} FAIL, {len(cnr)} could-not-run, of {len(results)}")
    # the red-proof checks are EXPECTED to fail on the deliberately faulted deck: invert them
    if "--redproof" in args:
        red = [r for r in results if "[red-proof]" in r[0]]
        for r in red:
            if r[1] == "FAIL":
                print(f"red-proof OK: {r[0]} fails as it must")
                fails.remove(r)
            elif r[1] == "pass":
                print(f"RED-PROOF BROKEN: {r[0]} passed on a faulted deck")
                fails.append(r)
    if fails:
        sys.exit(1)
    if cnr:
        sys.exit(2)
    sys.exit(0)


if __name__ == "__main__":
    main()
