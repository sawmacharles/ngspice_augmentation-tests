#!/usr/bin/env python3
"""card #533: shooting PSS / td-PAC / PSTB / PXF (forward and adjoint) / PNOISE under
`.options klu` vs the default SPARSE 1.3, same deck, same binary.

Each case writes two decks under ~/Work/data/ngspice_augmentation/533/eq/<case>/{sparse,klu}/
(fixtures hard-code wrdata paths; they are rewritten into that directory), runs both,
and compares:
  * .pss  : the orbit (every wrdata time-series column, KLU linearly interpolated onto the
            SPARSE time grid -- the two runs choose their own steps), deviation / peak-to-peak;
  * .pac/.pstb/.pxf/.pnoise : the complex (or real) transfer curves, max |dH|/|H| and max
            |dH|/max|H| over the sweep.
Pass = deviation <= TOL for the case (orbit TOL_ORBIT, transfer TOL_XFER); the tolerances are
roundoff-to-1e-4 class, see the header of docs/pss-pac.md "KLU" section.

RED-PROOF (--redproof): the same comparison with a deliberate fault -- the KLU deck's
period/clock is perturbed by 1 %% (`fs=495k` is not used; the .pss guess frequency is changed
for the KLU deck only) -- must make the check FAIL, proving the comparison can fail.

Usage:  NGSPICE=<binary> [NGSPICE_ASAN=<asan binary>] OMP_NUM_THREADS=2 python3 tests/test_pss_klu.py [--quick] [--redproof] [--only substr]
Exit 0 pass, 1 fail, 2 could-not-run.   Outputs: ~/Work/data/ngspice_augmentation/533/eq (never the repo).
"""
import os, re, sys, subprocess, time, math, cmath, json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIX = os.path.join(HERE, "fixtures")
NGSPICE = os.environ.get("NGSPICE", "ngspice")
ASAN = os.environ.get("NGSPICE_ASAN", "")
OUT = os.path.expanduser("~/Work/data/ngspice_augmentation/533/eq")
TOL_ORBIT = 1e-4   # smooth vectors: max |dev| / max|y| (the two runs use independent adaptive steps;
                   # a slow Floquet mode, rho~0.994 on R4, turns the 1e-6 residual stop into ~1e-5 of level)
TOL_EDGE = 2e-2    # edge vectors (pwm, ramp, switch node): max |dev| / p-p; both grids chord the corner/tail
TOL_XFER = 1e-4    # relative to the largest |H| of the sweep

CTRL = "set wr_singlescale\nset wr_vecnames\nset numdgt=15\n"


class CouldNotRun(Exception):
    pass


def body(fix):
    with open(os.path.join(FIX, fix)) as fh:
        s = fh.read()
    return s.replace("../../stb/lgprobe.lib", os.path.join(ROOT, "stb", "lgprobe.lib"))


def strip_end(s):
    return "\n".join(l for l in s.splitlines() if l.strip().lower() != ".end")


def strip_ac(s):
    return re.sub(r"\s+AC\s+\S+", "", s)


def with_ac1(s, dev):
    out, hit = [], False
    for ln in strip_ac(s).splitlines():
        if not hit and re.match(r"^\s*" + re.escape(dev) + r"\s", ln):
            ln += " AC 1"
            hit = True
        out.append(ln)
    if not hit:
        raise CouldNotRun(f"device {dev} not found")
    return "\n".join(out)


def add_klu(deck):
    """insert `.options klu` right before the first .control (or before .end)."""
    if re.search(r"^\.control", deck, re.M):
        return re.sub(r"^\.control", ".options klu\n.control", deck, count=1, flags=re.M)
    return re.sub(r"^\.end\s*$", ".options klu\n.end", deck, count=1, flags=re.M)


def redirect(deck, d):
    return re.sub(r"(?m)^(\s*(?:wrdata|write))\s+\S*?([^/\s]+)(\s|$)",
                  lambda m: f"{m.group(1)} {d}/{m.group(2)}{m.group(3)}", deck)


def run(binary, deck_path, env=None, timeout=1800):
    e = dict(os.environ)
    e.setdefault("OMP_NUM_THREADS", "2")
    if env:
        e.update(env)
    t0 = time.time()
    try:
        p = subprocess.run([binary, "-b", deck_path], capture_output=True, text=True,
                           timeout=timeout, cwd=os.path.dirname(deck_path), env=e)
    except subprocess.TimeoutExpired:
        raise CouldNotRun(f"{deck_path} timed out")
    return p.returncode, p.stdout + p.stderr, time.time() - t0


# ------------------------------------------------------------------ cases
def pss_case(fix, vecs=None):
    def mk(d, tag):
        s = redirect(body(fix), d)
        if vecs:  # i(Shs) of the r6s fixture is "not available" in plain wrdata (fixture quirk, same on SPARSE)
            s = re.sub(r"(?m)^(wrdata\s+\S+).*$", lambda m: m.group(1) + " " + vecs, s)
        s = re.sub(r"(?m)^\.control\s*$", ".control\nset numdgt=15", s, count=1)
        return s
    return dict(kind="pss", mk=mk, fix=fix)


def sweep_case(kind, fix, card, vecs, stab_env=None, ac_dev=None, strip=False, env=None, extra_check=None):
    def mk(d, tag):
        b = strip_end(body(fix))
        # drop any .control block the fixture itself carries
        b = re.sub(r"(?ms)^\.control.*?^\.endc\s*$", "", b)
        if ac_dev:
            b = with_ac1(b, ac_dev)
        if strip:
            b = strip_ac(b)
        return b + "\n" + card + "\n.control\n" + CTRL + "run\nwrdata " + os.path.join(d, "out.txt") + " " + vecs + "\nquit\n.endc\n.end\n"
    return dict(kind=kind, mk=mk, fix=fix, env=env)


R3 = ("500e3", "300e-6", "vout", "200", "10", "150", "1e-6")
R4 = ("500e3", "3500e-6", "vout", "200", "10", "150", "1e-6")
SW = "dec 5 500 200000"


def pac_card(p, sweep=SW):
    return f'.pac {" ".join(p)} {sweep} method="shoot" solver="td"'


def pxf_card(p, srcs, adj, out="vout", sweep=SW):
    st = " ".join(f'sources="{s}"' for s in srcs)
    return f'.pxf {" ".join(p)} {out} {sweep} method="shoot" solver="td" {st} adjoint={adj}'


EXACT = {"PAC_EXACT_JACOBIAN": "1"}
CASES = {
    "pss_r3_shoot": pss_case("pss_r3_shoot.cir"),
    "pss_r4_shoot": pss_case("pss_r4_shoot.cir"),
    "pss_r6s_shoot": pss_case("pss_r6s_shoot.cir", "v(vout) v(opout) i(L1)"),
    "pss_r2_shoot": pss_case("pss_r2_shoot.cir"),
    "pss_r2_q4_shoot": pss_case("pss_r2_q4_shoot.cir"),
    "pac_r3_ctrl2out": sweep_case("xfer", "pac_r3_ctrl2out.cir", pac_card(R3), "v(vout)"),
    "pac_r3_line2out": sweep_case("xfer", "pac_r3_line2out.cir", pac_card(R3), "v(vout)"),
    "pac_r3_zout": sweep_case("xfer", "pac_r3_zout.cir", pac_card(R3), "v(vout)"),
    "pac_r4_ctrl2out": sweep_case("xfer", "pac_r4_ctrl2out.cir", pac_card(R4), "v(vout)"),
    "pac_r4_line2out": sweep_case("xfer", "pac_r4_line2out.cir", pac_card(R4), "v(vout)"),
    "pac_r4_zout": sweep_case("xfer", "pac_r4_zout.cir", pac_card(R4), "v(vout)"),
    "pstb_r3": sweep_case("xfer", "pac_r3_ctrl2out.cir",
                          '.pstb ' + " ".join(R3) + ' ' + SW + ' probe="Vprobe" method="shoot"', "t", strip=True),
    "pxf_lti_adj0": sweep_case("xfer", "pxf_lti.cir", pxf_card(("100k", "0", "clkgnd", "2000", "1", "5", "0.001"),
                               ["Vin", "Iinj"], 0, out="out", sweep="dec 10 1k 1meg"), "pxf_vin pxf_iinj", env=EXACT),
    "pxf_lti_adj1": sweep_case("xfer", "pxf_lti.cir", pxf_card(("100k", "0", "clkgnd", "2000", "1", "5", "0.001"),
                               ["Vin", "Iinj"], 1, out="out", sweep="dec 10 1k 1meg"), "pxf_vin pxf_iinj", env=EXACT),
    "pxf_r3_adj0": sweep_case("xfer", "pac_r3_ctrl2out.cir", pxf_card(R3, ["Vin", "Vctrl"], 0), "pxf_vin pxf_vctrl", env=EXACT),
    "pxf_r3_adj1": sweep_case("xfer", "pac_r3_ctrl2out.cir", pxf_card(R3, ["Vin", "Vctrl"], 1), "pxf_vin pxf_vctrl", env=EXACT),
    "pnoise_r2_smooth": sweep_case("xfer", "pnoise_r2_smooth.cir",
                                   '.pnoise 100e3 100e-6 vx 2000 1 50 1e-6 vx Vin dec 5 1k 1meg method="shoot" solver="td" dev=1',
                                   "onoise_spectrum onoise_Rin onoise_Rsw"),
}
CASES["pnoise_r2_smooth"]["real"] = True
# --extra: further engines (free-running trigger= PSS, saltation/ideal-switch decks, PXF with saltation -> adjoint falls back)
EXTRA = {
    "pss_r5_freerun_trig": dict(pss_case("pss_r5_shoot.cir"), tail=0.95),
    "pss_r6o_shoot": pss_case("pss_r6o_shoot.cir"),
    "pxf_r2_adj0": sweep_case("xfer", "pxf_r2_switched_rc.cir", pxf_card(("100e3", "100e-6", "vx", "2000", "1", "25", "1e-6"), ["Vin", "Iinj"], 0, out="vx", sweep="dec 10 1k 40k"), "pxf_vin pxf_iinj"),
    "pxf_r2_adj1": sweep_case("xfer", "pxf_r2_switched_rc.cir", pxf_card(("100e3", "100e-6", "vx", "2000", "1", "25", "1e-6"), ["Vin", "Iinj"], 1, out="vx", sweep="dec 10 1k 40k"), "pxf_vin pxf_iinj"),
}
if "--extra" in sys.argv:
    CASES.update(EXTRA)
_unused = {
}
QUICK = ["pss_r3_shoot", "pss_r2_shoot", "pac_r3_ctrl2out", "pxf_lti_adj1", "pnoise_r2_smooth"]



# ------------------------------------------------------------------ PEX-size converter deck
def pex_deck(branch, depth, d, card, vecs, tag, stab="300e-6"):
    """R3 voltage-mode buck (tests/fixtures/pac_r3_ctrl2out.cir: corrected op-amp polarity
    `Goa 0 opout vrefnode compin`) + three extracted-parasitic RC trees (bench/gen_pex_buck.py
    tree_lines): a decoupling net on vout (10 ohm), a switch-node net on sw (0.5 ohm) and a
    ground-return net (cap return `gret`, 5 mohm to ground, tree on gret via 1 ohm)."""
    sys.path.insert(0, os.path.join(ROOT, "bench"))
    import gen_pex_buck as g
    b = strip_end(body("pac_r3_ctrl2out.cir"))
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
    # .save keeps the per-run transient plot small: the shooting engine runs ~one transient per
    # state (FD monodromy), each leaving a plot; with every parasitic node saved that is tens of GB.
    pex = "* PEX-size parasitic nets (3 RC trees)\n" + "".join(lines) + ".save v(vout) v(sw) v(opout) i(L1)\n"
    if card.startswith("pss"):   # control-command form, as the .pss fixtures
        deck = b + "\n" + pex + "\n.control\n" + CTRL + card + "\nwrdata " + os.path.join(d, "out.txt") + " " + vecs + "\nquit\n.endc\n.end\n"
        return deck, nelem
    deck = b + "\n" + pex + "\n" + card + "\n.control\n" + CTRL + "run\nwrdata " + os.path.join(d, "out.txt") + " " + vecs + "\nquit\n.endc\n.end\n"
    return deck, nelem


# ------------------------------------------------------------------ readers / metrics
def read_pairs(path):
    """wrdata without wr_singlescale: rows of (x1 y1 x2 y2 ...) -> list of (xs, ys)"""
    rows = []
    with open(path) as fh:
        for l in fh:
            t = l.split()
            if t:
                try:
                    rows.append([float(v) for v in t])
                except ValueError:
                    pass
    nc = len(rows[0]) // 2
    return [([r[2 * k] for r in rows], [r[2 * k + 1] for r in rows]) for k in range(nc)]


def interp(xs, ys, x):
    import bisect
    i = bisect.bisect_left(xs, x)
    if i <= 0:
        return ys[0]
    if i >= len(xs):
        return ys[-1]
    x0, x1 = xs[i - 1], xs[i]
    return ys[i - 1] if x1 == x0 else ys[i - 1] + (ys[i] - ys[i - 1]) * (x - x0) / (x1 - x0)


def wrap(t, x):
    """a periodic orbit's time shifted past either end of the record re-enters at the other end"""
    if t > x[-1]:
        t -= x[-1] - x[0]
    elif t < x[0]:
        t += x[-1] - x[0]
    return t


def orbit_dev(pa, pb, tail=1.0):
    """per vector: max |A - B| / p-p on a uniform grid over the common span (KLU curve and SPARSE
    curve both linearly interpolated).  A vector with edges (max slope > 15 x p-p/span, i.e. an edge
    faster than 1/15 of the record: pwm, sawtooth ramp) is compared edge-tolerantly (best of time shifts
    up to 2 x the coarser grid step), because an independent adaptive grid straddles the corner."""
    A, B = read_pairs(pa), read_pairs(pb)
    out = []
    for (xa, ya), (xb, yb) in zip(A, B):
        # time axis in units of the record length (the period): for a free-running (trigger=) orbit
        # the period is itself a Newton unknown and differs by the stop tolerance between runs
        xa = [t / xa[-1] for t in xa]
        xb = [t / xb[-1] for t in xb]
        pp = max(ya) - min(ya) or 1.0
        edge = max(abs(ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]) for i in range(len(ya) - 1) if xa[i + 1] > xa[i]) > 15 * pp / (xa[-1] - xa[0])
        lo, hi = max(xa[0], xb[0]), min(xa[-1], xb[-1])
        # tail < 1 (free-running trigger= decks): the record ends at the trigger event, which each
        # run locates to its own tolerance (steep edge sampled at a different phase) -- leave it out
        hi = lo + tail * (hi - lo)
        h = 2 * max(max(b - a for a, b in zip(xa, xa[1:])), max(b - a for a, b in zip(xb, xb[1:]))) if edge else 0.0
        n = 4000
        shifts = [h * k / 40 for k in range(-40, 41)] if edge else [0.0]
        dev = 0.0
        for k in range(n + 1):
            t = lo + (hi - lo) * k / n
            ta = interp(xa, ya, t)
            dev = max(dev, min(abs(ta - interp(xb, yb, wrap(t + s, xb))) for s in shifts))
        out.append((dev / pp if edge else dev / (max(abs(v) for v in ya) or 1.0), edge))
    return out


def read_sw(path):
    """wr_singlescale/wr_vecnames: header line, then freq then (re,im) per vector -- or real
    columns for noise (re only).  Returns freqs, [[complex-or-real per row] per column]."""
    if not os.path.exists(path):
        raise CouldNotRun(f"no output file {path} (the run did not get as far as wrdata)")
    with open(path) as fh:
        lines = [l.split() for l in fh if l.strip()]
    hdr, rows = lines[0], []
    for t in lines[1:]:
        try:
            rows.append([float(v) for v in t])
        except ValueError:
            pass
    ncol = len(rows[0]) - 1
    return rows, ncol


def xfer_dev(pa, pb, real=False):
    ra, na = read_sw(pa)
    rb, nb = read_sw(pb)
    if len(ra) != len(rb) or na != nb:
        raise CouldNotRun(f"shape mismatch sparse {len(ra)}x{na} vs klu {len(rb)}x{nb}")
    # complex pairs when the column count is even and the sweep is complex; real otherwise.
    # Decide per file: pnoise has an odd data-column count only if real+... -> treat by header
    with open(pa) as fh:
        hdr = fh.readline().split()
    # a complex vector's header carries its name once; real/imag are two numeric columns
    cplx = (not real) and na % 2 == 0   # complex vectors print (re, im) per vector; pnoise spectra are real
    step = 2 if cplx else 1
    worst_rel_pt, worst_rel_max, worst_db, worst_deg = 0.0, 0.0, 0.0, 0.0
    for c in range(0, na, step):
        va = [complex(r[1 + c], r[2 + c]) if cplx else r[1 + c] for r in ra]
        vb = [complex(r[1 + c], r[2 + c]) if cplx else r[1 + c] for r in rb]
        mx = max(abs(v) for v in va) or 1.0
        for x, y in zip(va, vb):
            d = abs(x - y)
            worst_rel_max = max(worst_rel_max, d / mx)
            if abs(x) > 1e-12 * mx:
                worst_rel_pt = max(worst_rel_pt, d / abs(x))
            if cplx and abs(x) > 0 and abs(y) > 0:
                worst_db = max(worst_db, abs(20 * math.log10(abs(y) / abs(x))))
                dp = abs(math.degrees(cmath.phase(y / x)))
                worst_deg = max(worst_deg, dp)
    return dict(rel_max=worst_rel_max, rel_pt=worst_rel_pt, db=worst_db, deg=worst_deg, cplx=cplx)


# ------------------------------------------------------------------ main
def run_case(name, binary, red=False):
    c = CASES[name]
    res = {}
    for mode in ("sparse", "klu"):
        d = os.path.join(OUT, name + ("_red" if red else ""), mode)
        os.makedirs(d, exist_ok=True)
        deck = c["mk"](d, mode)
        if mode == "klu":
            deck = add_klu(deck)
            if red:  # deliberate fault: the KLU deck gets a 1 % different clock guess
                deck = re.sub(r"(?m)^(\.(?:pac|pxf|pnoise|pstb)\s+)(\S+)", lambda m: m.group(1) + f"{float(m.group(2))*1.01:.6g}", deck, count=1) \
                    if c["kind"] == "xfer" else re.sub(r"(?m)^(pss\s+)(\S+)", lambda m: m.group(1) + f"{float(m.group(2))*1.01:.6g}", deck, count=1)
        p = os.path.join(d, "deck.cir")
        for f in os.listdir(d):   # no stale output of an earlier run may satisfy the comparison
            if f.endswith((".txt", ".raw")):
                os.remove(os.path.join(d, f))
        with open(p, "w") as fh:
            fh.write(deck)
        rc, out, wall = run(binary, p, env=c.get("env"))
        with open(os.path.join(d, "log.txt"), "w") as fh:
            fh.write(out)
        if rc != 0 or "not supported with" in out:
            why = " | ".join(l.strip() for l in out.splitlines() if re.search(r"rror|singular|not supported|FAILED", l))[:400]
            raise CouldNotRun(f"{name} {mode}: exit {rc}: {why or out[-300:]}")
        if "FAILED to converge" in out:
            raise CouldNotRun(f"{name} {mode}: PSS did not converge")
        res[mode] = (d, wall, out)
    ds, dk = res["sparse"][0], res["klu"][0]
    if c["kind"] == "pss":
        base = os.path.basename
        fl = [f for f in os.listdir(ds) if f.endswith(".txt") and f != "log.txt"]
        if not fl:
            raise CouldNotRun(f"{name}: no wrdata output ({res['sparse'][2][-300:]})")
        fn = fl[0]
        w = orbit_dev(os.path.join(ds, fn), os.path.join(dk, fn), tail=c.get("tail", 1.0))
        m = max(x / (TOL_EDGE if e else TOL_ORBIT) for x, e in w)   # normalised: pass <= 1
        return dict(name=name, metric=m, tol=1.0, detail="orbit dev per vec (smooth: /max|y|, e: /p-p, edge-tolerant) " + ",".join(f"{x:.2e}{'e' if e else ''}" for x, e in w), walls=(res["sparse"][1], res["klu"][1]))
    x = xfer_dev(os.path.join(ds, "out.txt"), os.path.join(dk, "out.txt"), real=c.get("real", False))
    return dict(name=name, metric=x["rel_max"], tol=TOL_XFER,
                detail=f"max|dH|/max|H|={x['rel_max']:.2e} max|dH|/|H|={x['rel_pt']:.2e}" + (f" dB={x['db']:.2e} deg={x['deg']:.2e}" if x["cplx"] else " (real)"),
                walls=(res["sparse"][1], res["klu"][1]))


def main():
    args = sys.argv[1:]
    names = QUICK if "--quick" in args else list(CASES)
    if "--only" in args:
        s = args[args.index("--only") + 1]
        names = [n for n in CASES if s in n]
    fails = cnr = 0
    rows = []
    for n in names:
        try:
            r = run_case(n, NGSPICE)
            ok = r["metric"] <= r["tol"]
            print(f"{'pass' if ok else 'FAIL'}  {n}: {r['detail']} (norm. metric {r['metric']:.2g}, pass <= {r['tol']:g}) wall sparse {r['walls'][0]:.2f}s klu {r['walls'][1]:.2f}s", flush=True)
            fails += (not ok)
            rows.append(r)
        except CouldNotRun as e:
            print(f"COULD-NOT-RUN  {n}: {e}", flush=True)
            cnr += 1
    if "--redproof" in args or "--quick" in args:
        for n in ("pss_r2_shoot", "pac_r3_ctrl2out"):
            try:
                r = run_case(n, NGSPICE, red=True)
                bad = r["metric"] > r["tol"]
                print(f"{'pass' if bad else 'FAIL'}  redproof {n}: deliberate 1% clock fault is DETECTED={bad} ({r['detail']})", flush=True)
                fails += (not bad)
            except CouldNotRun as e:
                # a deliberately wrong clock may make the engine itself fail: also a detection, but say so
                print(f"pass  redproof {n}: fault made the run fail outright (detected): {str(e)[:160]}", flush=True)
    if ASAN and ("--asan" in args or "--quick" in args):
        for n in ("pss_r3_shoot", "pac_r3_ctrl2out"):
            try:
                c = CASES[n]
                d = os.path.join(OUT, n + "_asan")
                os.makedirs(d, exist_ok=True)
                p = os.path.join(d, "deck.cir")
                with open(p, "w") as fh:
                    fh.write(add_klu(c["mk"](d, "klu")))
                rc, out, wall = run(ASAN, p, env={"ASAN_OPTIONS": "detect_leaks=0"}, timeout=3600)
                bad = rc != 0 or "AddressSanitizer" in out or "FAILED to converge" in out
                print(f"{'FAIL' if bad else 'pass'}  asan klu {n}: rc={rc} wall {wall:.1f}s", flush=True)
                if bad:
                    print(out[-1500:])
                fails += bad
            except CouldNotRun as e:
                print(f"COULD-NOT-RUN  asan {n}: {e}")
                cnr += 1
    print(f"\n{len(rows)} compared, {fails} FAIL, {cnr} could-not-run")
    sys.exit(1 if fails else (2 if cnr else 0))


if __name__ == "__main__":
    main()
