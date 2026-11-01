#!/usr/bin/env python3
"""Card #708: transient-noise ensemble reference for the ABSOLUTE output noise of `.pnoise ... method="shoot" solver="td"` on BSIM4
circuits (fixtures tests/fixtures/pnoise708/*.cir).  The pinned values of tests/test_pnoise_bsim4.py come from this tool.

BSIM4 has no transient noise in ngspice, so the reference replaces each noise generator by an EXPLICIT trnoise-driven source:

  part "res"  (BSIM4 channel noise OFF: ntnoi=0; flicker is off in every fixture): every resistor noisy=0 with a TRNOISE current
              source in parallel, white one-sided PSD 4kT/R.  Checks the engine's absolute scaling with the BSIM4 dynamics (charges,
              time-varying gm/gds) in the loop.  `.pnoise` reference: the same deck with ntnoi=0 and noisy resistors.
  part "mos"  (resistors noisy=0 in both): the BSIM4 channel thermal noise.  Its density S(t) along the periodic orbit is taken from the
              MODEL ITSELF, independently of the `.pnoise` engine: the noiseless steady-state orbit is simulated, and at each of its time
              points the transistor alone, with its four terminals forced to the orbit voltages, gets a plain `.noise` (the short-circuit
              drain-current noise, read through a 0 V source and a 1 ohm CCVS).  In the ensemble the channel noise is a B-source
              I = sqrt(S(t)) * w(t) between drain and source, w unit white TRNOISE, sqrt(S(t)) a periodic PWL voltage of that table
              (quasi-static noise: the density follows the instantaneous bias, which is exactly what BSIM4noise evaluates).
              `.pnoise` reference: the deck with ntnoi=1 and noisy=0 resistors.

trnoise normalisation (tools/pnoise_mc.py, checked again here by the "rc" calibration against plain `.noise`): TRNOISE(NA TS) is iid
N(0,NA^2) every TS, linearly interpolated -> one-sided PSD 2 NA^2 TS sinc^4(f TS); NA = sqrt(S0/(2 TS)) gives S0.
Every injected amplitude is multiplied by K (output noise far above the integ() rounding floor; the circuit stays linear in it) and the
`.pnoise` value is multiplied by K.  Two copies of the circuit (noisy, noiseless) share one transient; the boxcar average of their
difference (integ, sampled at dt) is Welch-averaged (Hann, 50 %), 9-bin bands (pnoise_mc.compare()).

    NGSPICE=<prefix>/bin/ngspice python3 tools/pnoise_mc_bsim4.py <circuit> <part> [nseeds] [K] [TS_ns] [tnoimod]
      circuit: cs | th | lo      part: res | mos | rc (rc = calibration of the chain, circuit ignored)
Outputs under ~/Work/data/ngspice_augmentation/708/mc, or $MC_BSIM4_DATA (never the repository).  At most 2 ngspice at once, MemoryMax=6G,
num_threads=1.

Card #715 (non-quasi-static channel): circuit "lo" = tests/fixtures/pnoise715/lo_nqs_bsim4.cir (BSIM4 trnqsmod=1, 1 MHz LO).  A circuit may
carry `fixdir` (fixture folder) and `scale`: every time of the ensemble (boxcar dt 200 ns, 20 ms per seed, 2 ms dropped) is divided by
`scale` (lo: 1000, i.e. dt 0.2 ns, 20 us, 2 us dropped; bins 610.35 kHz), TS is given in ns as before.  MC_PARAMS="k=v,k=v" rewrites
`.param` values of the fixture in every deck the tool writes (the ensemble, the density table and the .pnoise reference) and is added to the
tag.  Defaults reproduce the #708 runs unchanged.
"""
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pnoise_mc as mc  # noqa: E402

ROOT = os.path.dirname(HERE)
FIX = os.path.join(ROOT, "tests", "fixtures", "pnoise708")
DATA = os.environ.get("MC_BSIM4_DATA", os.path.expanduser("~/Work/data/ngspice_augmentation/708/mc"))
NGSPICE = os.environ.get("NGSPICE", os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice"))
KT4 = 4 * mc.K_BOLTZ * mc.T_KELVIN
# MC_SEED_OFFSET=n: seeds n+1..n+N instead of 1..N (independent ensembles of cases that share a transfer path); the tag gets "_s<n>"
SEED_OFFSET = int(os.environ.get("MC_SEED_OFFSET", "0"))
# card #715: MC_PARAMS="nq=0,lo=0.3" -> .param overrides in every deck (tag suffix "_p<...>")
MC_PARAMS = dict(kv.split("=", 1) for kv in os.environ.get("MC_PARAMS", "").split(",") if "=" in kv)
PTAG = ("_p" + "_".join(f"{k}{v}" for k, v in sorted(MC_PARAMS.items()))) if MC_PARAMS else ""
# card #715: MC_OPTIONS="reltol=1e-6 vntol=1e-9" -> extra .options of the ENSEMBLE transient only (tag suffix "_o<...>")
MC_OPTIONS = os.environ.get("MC_OPTIONS", "").strip()
if MC_OPTIONS:
    PTAG += "_o" + re.sub(r"[^A-Za-z0-9.+-]+", "", MC_OPTIONS.replace("=", ""))

# circuit: fixture, output node, period, settle time for the orbit, the .pnoise head (Fguess StabTime OscNode Points Harmonics SC_iter
# Steady_coeff), the input-source field, and the band centres compared
CIRCUITS = {
    "cs": dict(fix="cs_bsim4.cir", out="d", T=10e-6, settle=200e-6, head="100k 100u gs 2000 1 50 1e-6", insrc="Vg",
               targets=[4000, 10000, 16000, 22000, 28000, 35000, 43000]),
    "th": dict(fix="th_bsim4.cir", out="vout", T=10e-6, settle=200e-6, head="100k 100u clk 2000 1 50 1e-6", insrc="Vin",
               targets=[4000, 10000, 16000, 22000, 28000, 35000, 43000]),
    # card #715: NQS BSIM4 stage, 1 MHz LO; the whole ensemble time base / 1000 (bins 610.35 kHz, bands 4.3 .. 42.7 MHz)
    "lo": dict(fix="lo_nqs_bsim4.cir", fixdir="pnoise715", out="out", T=1e-6, settle=20e-6, head="1meg 20u out 2000 1 50 1e-6",
               insrc="Vin", targets=[4e6, 10e6, 16e6, 22e6, 28e6, 35e6, 43e6], scale=1000.0),
}


def run_ngspice(deck, timeout=3600, env=None):
    e = dict(os.environ if env is None else env)
    e["OMP_NUM_THREADS"] = "1"
    t0 = time.time()
    cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=6G", "-p", "MemorySwapMax=0", NGSPICE, "-b", deck]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=os.path.dirname(deck), env=e)
    except subprocess.TimeoutExpired:
        raise mc.CouldNotRun(f"{deck} timed out after {timeout}s")
    return p.returncode, p.stdout + p.stderr, time.time() - t0


mc.run_ngspice = run_ngspice      # pnoise_reference()/run_batch() of pnoise_mc use this runner


def run_batch2(jobs):
    from concurrent.futures import ThreadPoolExecutor
    def one(j):
        rc, out, wall = run_ngspice(j[0])
        if rc != 0:
            raise mc.CouldNotRun(f"{j[0]} exit {rc}: {out[-1000:]}")
        return wall
    with ThreadPoolExecutor(2) as ex:
        return list(ex.map(one, jobs))


mc.run_batch = lambda jobs, workers=2: run_batch2(jobs)


def fixture_lines(name, fixdir=None):
    with open(os.path.join(os.path.dirname(FIX), fixdir, name) if fixdir else os.path.join(FIX, name)) as fh:
        lines = [l for l in fh.read().splitlines() if l.strip() and not l.startswith("*") and l.strip().lower() != ".end"]
    return set_param(lines, **MC_PARAMS) if MC_PARAMS else lines


def cfix(cname):
    c = CIRCUITS[cname]
    return fixture_lines(c["fix"], c.get("fixdir"))


def set_param(lines, **kw):
    """Rewrite `.param` values (e.g. ntn=0, tnm=1)."""
    out = []
    for l in lines:
        if l.lower().startswith(".param"):
            for k, v in kw.items():
                l = re.sub(rf"\b{k}=\S+", f"{k}={v}", l)
        out.append(l)
    return out


def parse(lines):
    """-> params, models, resistors [(name,n1,n2,valuetext)], mosfets [(name,d,g,s,b,rest)], others (lines)."""
    par, mod, res, mos, oth = [], [], [], [], []
    for l in lines:
        t = l.split()
        c = t[0][0].upper()
        if l.lower().startswith(".param"):
            par.append(l)
        elif l.lower().startswith(".model"):
            mod.append(l)
        elif c == "R":
            res.append((t[0], t[1], t[2], t[3]))
        elif c == "M":
            mos.append((t[0], t[1], t[2], t[3], t[4], " ".join(t[5:])))
        else:
            oth.append(l)
    return par, mod, res, mos, oth


def resolve(val, params):
    v = val.strip("{}")
    if v in params:
        v = params[v]
    m = re.match(r"^([0-9.eE+-]+)([a-zA-Z]*)$", v)
    scale = {"": 1, "k": 1e3, "meg": 1e6, "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12}
    return float(m.group(1)) * scale[m.group(2).lower()]


def param_dict(par):
    d = {}
    for l in par:
        for k, v in re.findall(r"(\w+)=(\S+)", l):
            d[k] = v
    return d


# ------------------------------------------------------------------------------------------------------------------------------------
# orbit and the BSIM4 noise-density table

def orbit(cname, outdir):
    """Noiseless periodic orbit: node voltages of every MOS terminal over the last period (native time points, maxstep T/2000)."""
    c = CIRCUITS[cname]
    lines = cfix(cname)
    par, mod, res, mos, oth = parse(lines)
    nodes = sorted({n for m in mos for n in m[1:5] if n != "0"})
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, "orbit.cir")
    T, ts = c["T"], c["settle"]
    with open(deck, "w") as fh:
        fh.write("* orbit\n" + "\n".join(lines) + "\n")
        fh.write(f".options maxstep={T/2000:.6g} reltol=1e-6\n.tran {T/2000:.6g} {ts + T:.9g}\n.control\nset num_threads=1\nrun\n"
                 "set wr_singlescale\nset wr_vecnames\nset numdgt=15\nwrdata orbit.txt " + " ".join(f"v({n})" for n in nodes)
                 + "\nquit\n.endc\n.end\n")
    rc, out, _ = run_ngspice(deck, timeout=600)
    if rc != 0:
        raise mc.CouldNotRun(out[-1500:])
    rows, prev = [], None
    with open(os.path.join(outdir, "orbit.txt")) as fh:
        next(fh)
        for r in fh:
            v = [float(x) for x in r.split()]
            if v[0] >= ts - 1e-15:
                if not rows and prev is not None and v[0] > ts + 1e-15:
                    # card #715: no accepted point at exactly ts (seen on the 1 MHz lo deck): interpolate the period's start point
                    w = (ts - prev[0]) / (v[0] - prev[0])
                    rows.append([ts] + [a + w * (b - a) for a, b in zip(prev[1:], v[1:])])
                rows.append(v)
            else:
                prev = v
    # phase within the period, dedupe
    pts = {}
    for v in rows:
        ph = v[0] - ts
        if ph > T:
            ph = T
        pts[round(ph, 15)] = dict(zip(nodes, v[1:]))
    tl = sorted(pts)
    if tl[0] > 1e-15 or abs(tl[-1] - T) > 1e-12:
        raise mc.CouldNotRun(f"orbit window {tl[0]}..{tl[-1]} does not cover the period")
    # periodicity of the orbit itself (settled?): first vs last point
    worst = max(abs(pts[tl[0]][n] - pts[tl[-1]][n]) for n in nodes)
    return tl, pts, worst


def density_table(cname, outdir, tnm=0):
    """For every MOS: list of (phase, S_id [A^2/Hz]) from plain .noise at the orbit's terminal voltages (ntnoi=1, given tnoimod)."""
    c = CIRCUITS[cname]
    lines = set_param(cfix(cname), ntn=1, tnm=tnm)
    par, mod, res, mos, oth = parse(lines)
    tl, pts, per_err = orbit(cname, outdir)
    tables = {}
    for (name, d, g, s, b, rest) in mos:
        deck = os.path.join(outdir, f"table_{name}_tnm{tnm}.cir")
        txt = os.path.join(outdir, f"table_{name}_tnm{tnm}.txt")
        if os.path.exists(txt):
            os.remove(txt)
        V = lambda node, k: 0.0 if node == "0" else pts[tl[k]][node]
        with open(deck, "w") as fh:
            fh.write("* BSIM4 drain-current noise at forced terminal voltages\n" + "\n".join(par + mod) + "\n")
            fh.write("Vtd td 0 DC 0\nVtg tg 0 DC 0 AC 1\nVts ts 0 DC 0\nVtb tb 0 DC 0\nVm td tdi DC 0\n")
            fh.write(f"{name} tdi tg ts tb {rest}\nH1 hout 0 Vm 1\n")
            fh.write(".control\nset num_threads=1\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nset appendwrite\n")
            for k in range(len(tl)):
                fh.write(f"alter Vtd dc={V(d, k):.12g}\nalter Vtg dc={V(g, k):.12g}\nalter Vts dc={V(s, k):.12g}\n"
                         f"alter Vtb dc={V(b, k):.12g}\nnoise v(hout) Vtg lin 1 1k 1k\nwrdata {os.path.basename(txt)} onoise_spectrum\n"
                         "destroy all\n")
            fh.write("quit\n.endc\n.end\n")
        rc, out, _ = run_ngspice(deck, timeout=1200)
        if rc != 0:
            raise mc.CouldNotRun(out[-1500:])
        vals = []
        with open(txt) as fh:
            for r in fh:
                t = r.split()
                if len(t) == 2 and not t[0].startswith("f"):
                    vals.append(float(t[1]) ** 2)
        if len(vals) != len(tl):
            raise mc.CouldNotRun(f"table {name}: {len(vals)} values for {len(tl)} points")
        tables[name] = list(zip(tl, vals))
    return tables, per_err


def time_avg(table):
    s = 0.0
    for (t0, a), (t1, b) in zip(table, table[1:]):
        s += 0.5 * (a + b) * (t1 - t0)
    return s / (table[-1][0] - table[0][0])


# ------------------------------------------------------------------------------------------------------------------------------------
# ensemble deck

def mc_deck(cname, part, outdir, seed, ts, duration, dt, K, tables=None):
    c = CIRCUITS[cname]
    lines = set_param(cfix(cname), ntn=0)
    par, mod, res, mos, oth = parse(lines)
    pd = param_dict(par)
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, f"mc_seed{seed}.cir")
    out = c["out"]
    ports = [out] + ([f"sq_{m[0]}" for m in mos] if part == "mos" else [])
    with open(deck, "w") as fh:
        fh.write(f"* card #708 ensemble, {cname} part {part}\n" + "\n".join(par) + "\n")
        if part == "mos":
            T = c["T"]
            for m in mos:
                tab = tables[m[0]]
                pw = " ".join(f"{t:.12g} {math.sqrt(max(s, 0.0)):.9e}" for t, s in tab)
                fh.write(f"Vsq_{m[0]} sq_{m[0]} 0 PWL({pw}) r=0\n")
        fh.write(f".subckt core {' '.join(ports)} params: na=1\n")
        fh.write("\n".join(mod + oth) + "\n")
        for (name, n1, n2, val) in res:
            fh.write(f"{name} {n1} {n2} {val} noisy=0\n")
            if part == "res":
                na = mc.trnoise_na(KT4 / resolve(val, pd), ts)
                fh.write(f"I{name}_nz {n1} {n2} DC 0 TRNOISE({{na*{K * na:.9e}}} {ts:.9e} 0 0)\n")
        for (name, d, g, s, b, rest) in mos:
            fh.write(f"{name} {d} {g} {s} {b} {rest}\n")
            if part == "mos":
                fh.write(f"V{name}_nz nz_{name} 0 DC 0 TRNOISE({{na*{K:.9e}}} {ts:.9e} 0 0)\n")
                fh.write(f"B{name}_nz {d} {s} I = {1 / math.sqrt(2 * ts):.9e}*V(sq_{name})*V(nz_{name})\n")
        fh.write(".ends core\n")
        rest_ports = " ".join(ports[1:])
        fh.write(f"Xn outn {rest_ports} core na=1\nXd outd {rest_ports} core na=0\n")
        fh.write(f".save v(outn) v(outd)\n.options maxstep={ts:.6g} {MC_OPTIONS}\n.tran {dt:.9e} {duration:.9e}\n.control\nset num_threads=1\n"
                 f"setseed {seed + SEED_OFFSET}\nrun\nlet iv = integ(v(outn)-v(outd))\nlinearize iv\n"
                 "set wr_singlescale\nset wr_vecnames\nset numdgt=15\nwrdata mc_out.txt iv\nquit\n.endc\n.end\n")
    return deck, os.path.join(outdir, "mc_out.txt")


# BSIM4 creates the internal S/D nodes that carry the tnoimod=1 rs/rd ("theta") part of the channel noise only when a plain .noise job is in
# the deck (b4set.c noiseAnalGiven, card #708 defect); this inert card makes a .pnoise deck complete.  Fixed by card #729 (builds from w729 /
# the release after s144: CKTnoiseAnalGiven() also counts the .pnoise job), so `wk` is only needed on s144 and older binaries.
DUMMY_NOISE = "Vdum708 ndum708 0 DC 0 AC 1\nRdum708 ndum708 0 1k noisy=0\n.noise v(ndum708) Vdum708 lin 1 1k 1k"


def pnoise_ref_body(cname, part, tnm=0, workaround=False):
    c = CIRCUITS[cname]
    lines = set_param(cfix(cname), ntn=(0 if part == "res" else 1), tnm=tnm)
    if part == "mos":
        lines = [l + " noisy=0" if l[0] in "Rr" else l for l in lines]
    if workaround:
        lines.append(DUMMY_NOISE)
    return f"* card #708 .pnoise reference, {cname} part {part} tnoimod {tnm}\n" + "\n".join(lines)


def run(cname, part, nseeds=16, K=100.0, TS=10e-9, tnm=0, dur=20e-3, drop=2e-3):
    sc = CIRCUITS.get(cname, {}).get("scale", 1.0)
    dt, nfft, dur, drop = 200e-9 / sc, 8192, dur / sc, drop / sc
    if part == "rc":   # scaled circuits (card #715): C / scale, the circuit's own bands, tag rc_<circuit>_ts..
        r, cc = 10e3, 1e-9 / sc
        fc = 1 / (2 * math.pi * r * cc)
        mc.DATA_ROOT = DATA
        an = lambda fl: {f: math.sqrt(KT4 * r / (1 + (f / fc) ** 2)) for f in fl}
        return mc.compare((f"rc_{cname}" if sc != 1.0 else "rc") + f"_ts{TS*1e9:g}", mc.rc_calib_deck, None, "vx", None, an, nseeds, TS,
                          dur, dt, nfft, drop, CIRCUITS[cname if sc != 1.0 else "cs"]["targets"], "analytic", r=r, c=cc)
    c = CIRCUITS[cname]
    tag = f"{cname}_{part}_tnm{tnm}_K{K:g}_ts{TS*1e9:g}_n{nseeds}" + (f"_s{SEED_OFFSET}" if SEED_OFFSET else "") + PTAG
    mc.DATA_ROOT = DATA
    tables = None
    if part == "mos":
        tables, per_err = density_table(cname, os.path.join(DATA, f"{cname}_table"), tnm)
        for n, tab in tables.items():
            print(f"{n}: {len(tab)} table points, <S_id> = {time_avg(tab):.5e} A^2/Hz, max {max(s for _, s in tab):.5e}, "
                  f"orbit periodicity error {per_err:.2e} V")
    builder = lambda outdir, seed, ts, duration, dt_, **_: mc_deck(cname, part, outdir, seed, ts, duration, dt_, K, tables)
    body = pnoise_ref_body(cname, part, tnm)
    ref = lambda fl: {f: K * v for f, v in mc.pnoise_reference(body, os.path.join(DATA, tag, "pnoise_ref"), c["out"], c["insrc"],
                                                                 (c["head"],), fl).items()}
    return mc.compare(tag, builder, body, c["out"], (c["head"],), ref, nseeds, TS, dur, dt, nfft, drop, c["targets"],
                      ".pnoise*K")


def compare_saved(cname, part, nseeds=16, K=100.0, TS=10e-9, tnm=0, nseg=320, workaround=False):
    """Band comparison from a finished ensemble's S_avg.txt (same arithmetic as pnoise_mc.compare())."""
    c = CIRCUITS[cname]
    dt, nfft, NB = 200e-9 / c.get("scale", 1.0), 8192, 9
    tag = f"{cname}_{part}_tnm{tnm}_K{K:g}_ts{TS*1e9:g}_n{nseeds}" + (f"_s{SEED_OFFSET}" if SEED_OFFSET else "") + PTAG
    S = [float(l.split()[1]) for l in open(os.path.join(DATA, tag, "S_avg.txt"))]
    binw = 1.0 / (nfft * dt)
    n_eff = nseg * NB / (1.83 * 1.056)
    _, hib, lob = mc.mc_relative_spread_db(n_eff)
    body = pnoise_ref_body(cname, part, tnm, workaround)
    print(f"=== {tag} (from S_avg.txt, {nseg} segments; 9-bin band 1-sigma +{hib:.3f}/-{lob:.3f} dB) ===")
    print(f"{'centre':>9s} {'ensemble/K':>12s} {'.pnoise':>12s} {'d(dB)':>8s} {'d/sigma':>8s}")
    rows = []
    for t in c["targets"]:
        k0 = max(1, round(t / binw) - NB // 2)
        ks = list(range(k0, k0 + NB))
        fl = [k * binw for k in ks]
        ref = mc.pnoise_reference(body, os.path.join(DATA, tag, "pnoise_ref" + ("_wk" if workaround else "")), c["out"], c["insrc"], (c["head"],), fl)
        mcp = sum(S[k] / mc.sinc_droop(k * binw, dt) for k in ks) / NB / K ** 2
        refp = sum(ref[f] ** 2 for f in fl) / NB
        d = 10 * math.log10(refp / mcp)
        rows.append((fl[NB // 2], math.sqrt(mcp), math.sqrt(refp), d))
        print(f"{fl[NB//2]:9.1f} {math.sqrt(mcp):12.5e} {math.sqrt(refp):12.5e} {d:+8.3f} {d/(hib if d > 0 else lob):+8.2f}")
    print(f"(.pnoise - ensemble) worst |d| {max(abs(r[3]) for r in rows):.3f} dB, mean {sum(r[3] for r in rows)/len(rows):+.3f} dB")
    return rows, (lob, hib)


def main():
    a = sys.argv[1:]
    if a and a[0] == "--saved":
        a = a[1:]
        wk = "wk" in a
        a = [x for x in a if x != "wk"]
        return compare_saved(a[0], a[1], int(a[2]), float(a[3]), float(a[4]) * 1e-9, int(a[5]) if len(a) > 5 else 0, workaround=wk)
    cname, part = a[0], a[1]
    nseeds = int(a[2]) if len(a) > 2 else 16
    K = float(a[3]) if len(a) > 3 else 100.0
    TS = float(a[4]) * 1e-9 if len(a) > 4 else 10e-9
    tnm = int(a[5]) if len(a) > 5 else 0
    run(cname, part, nseeds, K, TS, tnm)


if __name__ == "__main__":
    main()
