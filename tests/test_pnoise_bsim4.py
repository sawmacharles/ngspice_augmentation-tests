#!/usr/bin/env python3
"""Card #708: ABSOLUTE output noise of `.pnoise ... method="shoot" solver="td"` on BSIM4 circuits against an independent
transient-noise ensemble (docs/pss-pac.md, "Absolute .pnoise on BSIM4 vs a transient-noise ensemble (card #708, 2026-10-19)").

Fixtures tests/fixtures/pnoise708/:
  cs_bsim4.cir      BSIM4 common-source stage, gate clocked at 100 kHz through Rg 10k, Rd 20k / Cl 100p: weak conduction <-> triode
                    each period (modulated channel noise AND a periodic output pole).  Output v(d).
  th_bsim4.cir      track-and-hold, BSIM4 switch (W/L 1u/0.5u, Ron ~1k) in series with Rtrack 1k onto 100p, fclk 100 kHz, 20 % track.
  static_bsim4.cir  the cs stage at a DC gate bias, flicker on, PSS period from an unrelated 1 kHz clock (an LTI circuit).

Checks (N/M pass printed at the end; exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN):
  static: `.pnoise dev=1` (td) total and every generator (BSIM4 id, 1/f, rs, rd, corl, Rd, Rg) == plain `.noise` (1e-5 of the total),
       tnoimod 0; tnoimod 1 with the inert .noise card (DUMMY_NOISE); tnoimod 1 WITHOUT it (card #729: the defect found by #708 is
       fixed, the rs/rd part is there); tnoimod 2
  static fd (card #729): the frequency-domain `.pnoise` total == plain `.noise` (1e-5 of the total, 1k/10k/100k) for BSIM4 versions
       4.5 / 4.6.5 / 4.7 / 4.8.2 (the bsim4v5, bsim4v6, bsim4v7, bsim4 copies of the setup), tnoimod 1 and 2, no workaround card
  3-6  ensemble, per 9-bin band (4.3 .. 42.7 kHz, sideband-0 time-averaged PSD), .pnoise vs the PINNED ensemble value, |d| <= 0.5 dB
       per band and |mean d| <= 0.2 dB:
         cs res   BSIM4 channel noise off (ntnoi=0), Rd and Rg thermal noise (explicit white trnoise currents in the ensemble)
         cs mos   resistors noiseless, BSIM4 channel thermal noise (tnoimod 0): in the ensemble a B-source sqrt(S(t)) w(t), S(t) the
                  model's own drain-current noise density at the orbit's terminal voltages (plain .noise per time point)
         th res / th mos   the same split on the track-and-hold
         cs mos tnoimod 1, with the inert .noise card and (card #729) without it
  red  PAC_NOISE_STALE_STATE=1 (the device noise read from the end-of-period state, the pre-#544 bug) on cs mos must FAIL the check;
       a deck with ntnoi=1.25 (+0.97 dB of channel noise) must FAIL the check;
       PSS_SALT_NO729=1 (card #729: BSIM4 creates the tnoimod noise nodes only for a job named "Noise Analysis", the pre-#729 setup)
       must FAIL the static tnoimod 1 check, the static fd 4.8.2 check and the cs mos tnoimod 1 check without the card (-4.4 dB).

The ensembles are too slow for the test (16 seeds x 20 ms, 10-40 min each); they were produced by
    NGSPICE=~/Work/apps/ngspice_augmentation/aaba50ba3f49-s143/bin/ngspice python3 tools/pnoise_mc_bsim4.py <cs|th> <res|mos> 16 100 <TS_ns> [tnoimod]
(th mos with MC_SEED_OFFSET=100; `tools/pnoise_mc_bsim4.py --saved <args> [wk]` redoes the band comparison from a saved spectrum, wk = with the
inert .noise card)
(K = 100 amplitude scale, TS = trnoise step; Hann/Welch 8192 x 200 ns boxcar samples, 320 segments, 9-bin bands; logs and spectra in
~/Work/data/ngspice_augmentation/708/mc).  The values below are the ensemble's band power / K^2 as V/rtHz with the band's 1-sigma in dB
(chi-squared argument of tools/pnoise_mc.py compare(), 320 segments x 9 correlated bins; the chain was calibrated on an RC against
plain .noise at TS = 10 ns: worst band 0.069 dB, mean -0.014 dB).

    NGSPICE=<prefix>/bin/ngspice python3 tests/test_pnoise_bsim4.py      (default: the s152 release; the #729 checks need s145 or later)
Outputs under ~/Work/data/ngspice_augmentation/708/test, or $PNOISE_BSIM4_DATA (never the repository).
"""
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures", "pnoise708")
NGSPICE = os.environ.get("NGSPICE", os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice"))
DATA = os.environ.get("PNOISE_BSIM4_DATA", os.path.expanduser("~/Work/data/ngspice_augmentation/708/test"))

BINW = 1.0 / (8192 * 200e-9)          # Welch bin of the ensembles, 610.35 Hz
NB = 9
TARGETS = [4000, 10000, 16000, 22000, 28000, 35000, 43000]
BAND_DB, MEAN_DB = 0.5, 0.2

CIRC = {
    "cs": dict(fix="cs_bsim4.cir", out="d", head="100k 100u gs 2000 1 50 1e-6", insrc="Vg"),
    "th": dict(fix="th_bsim4.cir", out="vout", head="100k 100u clk 2000 1 50 1e-6", insrc="Vin"),
}

# (circuit, part, tnoimod) -> (sigma_dB of one band, {band centre Hz: ensemble V/rtHz})
# ensemble band centres 4272.5, 9765.6, 15869.1, 21972.7, 28076.2, 34790.0, 42724.6 Hz (keyed by TARGETS); 16 seeds x 20 ms, K = 100
_T = TARGETS
ENS = {
    # TS = 10 ns, wall 382 s (2 ngspice)
    ("cs", "res", 0): (0.113, dict(zip(_T, [1.84330e-08, 1.94041e-08, 1.87399e-08, 1.89165e-08, 1.83087e-08, 1.80862e-08, 1.80278e-08]))),
    # TS = 10 ns, S(t) table 2015 orbit points, <S_id> = 1.98962e-24 A^2/Hz, wall 1218 s
    ("cs", "mos", 0): (0.113, dict(zip(_T, [1.34517e-08, 1.35831e-08, 1.36081e-08, 1.33446e-08, 1.34615e-08, 1.30029e-08, 1.29124e-08]))),
    # tnoimod=1, TS = 10 ns, table <S_id> = 1.87647e-24 A^2/Hz (plain .noise, which includes the rs/rd part), wall 1183 s; .pnoise without
    # the inert .noise card is 4.35..4.55 dB LOW against it (the defect), with it -0.11..+0.07 dB.  Same seeds as cs mos tnoimod 0.
    ("cs", "mos", 1): (0.113, dict(zip(_T, [1.37950e-08, 1.38702e-08, 1.38470e-08, 1.36500e-08, 1.37548e-08, 1.32923e-08, 1.31843e-08]))),
    # TS = 5 ns (tau ~200 ns while tracking), wall 703 s
    ("th", "res", 0): (0.113, dict(zip(_T, [1.34065e-08, 1.26217e-08, 1.30632e-08, 1.21579e-08, 1.19278e-08, 1.12362e-08, 1.01981e-08]))),
    # TS = 5 ns, S(t) table 2104 orbit points, <S_id> = 1.27318e-24 A^2/Hz, wall 1526 s, seeds 101..116 (MC_SEED_OFFSET=100): seeds 1..16
    # share th res's noise sequences and transfer path, so that run (worst 0.283 dB, mean +0.021) was correlated with th res; this one is not
    ("th", "mos", 0): (0.113, dict(zip(_T, [2.15839e-08, 2.13657e-08, 2.04466e-08, 2.01397e-08, 1.91086e-08, 1.85174e-08, 1.72342e-08]))),
}

results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))


def band_bins(t):
    k0 = max(1, round(t / BINW) - NB // 2)
    return [k * BINW for k in range(k0, k0 + NB)]


def fixture_lines(name):
    with open(os.path.join(FIX, name)) as fh:
        return [l for l in fh.read().splitlines() if l.strip().lower() != ".end"]


def set_params(lines, **kw):
    import re
    out = []
    for l in lines:
        if l.lower().startswith(".param"):
            for k, v in kw.items():
                l = re.sub(rf"\b{k}=\S+", f"{k}={v}", l)
        out.append(l)
    return out


def run_ngspice(deck, env=None, timeout=900):
    e = dict(os.environ)
    e["OMP_NUM_THREADS"] = "1"
    e.update(env or {})
    t0 = time.time()
    try:
        p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=6G", "-p", "MemorySwapMax=0", NGSPICE, "-b", deck],
                           capture_output=True, text=True, timeout=timeout, cwd=os.path.dirname(deck), env=e)
    except FileNotFoundError as ex:
        raise CouldNotRun(f"cannot start: {ex}")
    except subprocess.TimeoutExpired:
        raise CouldNotRun(f"{deck} timed out")
    return p.returncode, p.stdout + p.stderr, time.time() - t0


def read_cols(path, n):
    if not os.path.exists(path):
        raise CouldNotRun(f"ngspice wrote no {path}")
    rows = []
    with open(path) as fh:
        next(fh)
        for r in fh:
            v = r.split()
            if v:
                if len(v) != 1 + n:
                    raise CouldNotRun(f"{path}: {len(v)} columns, expected {1 + n}")
                rows.append([float(x) for x in v])
    if not rows:
        raise CouldNotRun(f"{path} is empty")
    return rows


# The defect found by card #708 (docs), fixed by card #729: BSIM4 (b4set.c and the v5/v6/v7 copies) created the internal source/drain nodes
# that carry the tnoimod=1 "theta" part of the channel noise (the source/drain resistance generator rs/rd, conductance reset to 1e3 S) only
# when the job list held a job named "Noise Analysis", i.e. a plain .noise card.  A .pnoise deck has none, so with tnoimod=1 that part was
# silently missing.  Adding this inert .noise card on a separate node makes BSIM4 create the nodes.  Since #729 the setup asks
# CKTnoiseAnalGiven() (cktdojob.c), which also answers yes for a .pnoise (and .sp donoise) job; PSS_SALT_NO729=1 restores the old test.
NO729 = {"PSS_SALT_NO729": "1"}
DUMMY_NOISE = "Vdum708 ndum708 0 DC 0 AC 1\nRdum708 ndum708 0 1k noisy=0\n.noise v(ndum708) Vdum708 lin 1 1k 1k\n"


def pnoise_bands(cname, part, tnm=0, env=None, tag="", ntn_on=1, workaround=False):
    """.pnoise (solver=td) over every Welch bin of the 7 bands in ONE run; returns {centre: band-power-mean V/rtHz}.
    workaround: add DUMMY_NOISE (see above)."""
    c = CIRC[cname]
    lines = set_params(fixture_lines(c["fix"]), ntn=(0 if part == "res" else ntn_on), tnm=tnm)
    if part == "mos":
        lines = [l + " noisy=0" if l[:1] in "Rr" else l for l in lines]
    kfirst = round(band_bins(TARGETS[0])[0] / BINW)
    klast = round(band_bins(TARGETS[-1])[-1] / BINW)
    outdir = os.path.join(DATA, f"{cname}_{part}_tnm{tnm}{tag}")
    os.makedirs(outdir, exist_ok=True)
    txt = os.path.join(outdir, "pn.txt")
    if os.path.exists(txt):
        os.remove(txt)
    deck = os.path.join(outdir, "pn.cir")
    with open(deck, "w") as fh:
        fh.write("\n".join(lines) + "\n" + (DUMMY_NOISE if workaround else ""))
        fh.write(f'.pnoise {c["head"]} {c["out"]} {c["insrc"]} lin {klast - kfirst + 1} {kfirst * BINW:.10g} {klast * BINW:.10g} '
                 'method="shoot" solver="td"\n')
        fh.write(".control\nset num_threads=1\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nwrdata pn.txt onoise_spectrum\n"
                 "quit\n.endc\n.end\n")
    rc, out, wall = run_ngspice(deck, env=env)
    if rc != 0 or "PNOISE-TD: done" not in out:
        raise CouldNotRun(f"{cname} {part} .pnoise exit {rc}: {out[-1500:]}")
    rows = read_cols(txt, 1)
    if len(rows) != klast - kfirst + 1:
        raise CouldNotRun(f"{len(rows)} rows in {txt}, expected {klast - kfirst + 1} (wrong plot written?)")
    byk = {round(r[0] / BINW): r[1] for r in rows}
    res = {}
    for t in TARGETS:
        ks = [round(f / BINW) for f in band_bins(t)]
        res[t] = math.sqrt(sum(byk[k] ** 2 for k in ks) / NB)
    return res, wall


def check_ensemble(cname, part, tnm=0, env=None, red=None, ntn_on=1, workaround=False, known=None):
    """known: a text naming a known defect -- the case is then recorded as info (not counted), with the numbers."""
    key = (cname, part, tnm)
    label = {"res": "resistor noise, BSIM4 channel noise off", "mos": f"BSIM4 channel thermal noise (tnoimod {tnm})"}[part]
    name = f"{cname} {part}: .pnoise vs transient-noise ensemble, {label}"
    if workaround:
        name += " [with the inert .noise card, see DUMMY_NOISE]"
    elif part == "mos" and tnm == 1:
        name += " [no workaround card, #729]"
    if known:
        name = f"[KNOWN DEFECT: {known}] " + name
    if red:
        name = f"[RED PROOF: {red}, must FAIL the check] " + name
    if key not in ENS:
        record(name, "COULD-NOT-RUN", "no pinned ensemble for this case")
        return
    sig, ens = ENS[key]
    try:
        tag = ("_red" + str(abs(hash(red)) % 1000) if red else "") + ("_wk" if workaround else "") + ("_known" if known else "")
        pn, wall = pnoise_bands(cname, part, tnm, env, tag=tag, ntn_on=ntn_on, workaround=workaround)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    ds = [20 * math.log10(pn[t] / ens[t]) for t in TARGETS]
    worst = max(abs(d) for d in ds)
    mean = sum(ds) / len(ds)
    ok = worst <= BAND_DB and abs(mean) <= MEAN_DB
    rows = ", ".join(f"{t / 1e3:g}k {d:+.3f}" for t, d in zip(TARGETS, ds))
    detail = (f".pnoise - ensemble (dB) {rows}; worst {worst:.3f}, mean {mean:+.3f} (limits {BAND_DB}/{MEAN_DB}; band 1-sigma "
              f"{sig:.3f} dB, mean 1-sigma {sig / math.sqrt(len(ds)):.3f}); 4.3k: {pn[TARGETS[0]]:.4e} vs {ens[TARGETS[0]]:.4e} V/rtHz; "
              f"wall {wall:.1f}s")
    if known:
        record(name, "info", detail + ("; the defect is still present" if not ok else "; NOW WITHIN LIMITS: the defect looks fixed, update the test"))
    elif red:
        record(name, "pass" if not ok else "FAIL", detail)
    else:
        record(name, "pass" if ok else "FAIL", detail)


def check_static(tnm, bound=1e-5, workaround=False, known=None, env=None, red=None):
    """DC-biased BSIM4 stage (LTI) with flicker on: .pnoise dev=1 (td) total and every generator vs plain .noise.
    red: a red proof (env applied to both runs) -- recorded pass when the check FAILS."""
    name = f"static: DC-biased BSIM4 stage, tnoimod {tnm}, flicker on, .pnoise dev=1 vs plain .noise (total + every generator)"
    if workaround:
        name += " [with the inert .noise card]"
    elif tnm == 1:
        name += " [no workaround card, #729]"
    if known:
        name = f"[KNOWN DEFECT: {known}] " + name
    if red:
        name = f"[RED PROOF: {red}, must FAIL the check] " + name
    vecs = ["onoise.m1.id", "onoise.m1.1overf", "onoise.m1.rs", "onoise.m1.rd", "onoise.m1.corl", "onoise_rd", "onoise_rg"]
    try:
        body = "\n".join(set_params(fixture_lines("static_bsim4.cir"), tnm=tnm))
        outdir = os.path.join(DATA, f"static_tnm{tnm}" + ("_wk" if workaround else "") + ("_known" if known else "") + ("_red729" if red else ""))
        os.makedirs(outdir, exist_ok=True)
        opts = ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n"
        ctl = ".control\nset num_threads=1\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
        nz = os.path.join(outdir, "nz.cir")
        with open(nz, "w") as fh:
            fh.write(body.replace("Vg gs 0 DC 1.0", "Vg gs 0 DC 1.0 AC 1") + "\n" + opts + ".noise v(d) Vg dec 5 1k 100k 1\n" + ctl
                     + "setplot noise1\nwrdata nz.txt onoise_spectrum " + " ".join(vecs) + "\nquit\n.endc\n.end\n")
        pn = os.path.join(outdir, "pn.cir")
        with open(pn, "w") as fh:
            fh.write(body + "\n" + (DUMMY_NOISE if workaround else "") + opts + '.pnoise 1k 0 clkgnd 2000 1 8 1e-6 d Vg dec 5 1k 100k method="shoot" solver="td" dev=1 '
                     "pacpoints=200000\n" + ctl + "wrdata pn.txt onoise_spectrum " + " ".join(vecs) + "\nquit\n.endc\n.end\n")
        for d, f in ((nz, "nz.txt"), (pn, "pn.txt")):
            p = os.path.join(outdir, f)
            if os.path.exists(p):
                os.remove(p)
            rc, out, wall = run_ngspice(d, env=env)
            if rc != 0 or "REFUSED" in out:
                raise CouldNotRun(f"{os.path.basename(d)} exit {rc}: {out[-1500:]}")
        a = read_cols(os.path.join(outdir, "nz.txt"), 1 + len(vecs))
        b = read_cols(os.path.join(outdir, "pn.txt"), 1 + len(vecs))
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    if len(a) != len(b):
        record(name, "COULD-NOT-RUN", f"{len(a)} .noise points vs {len(b)} .pnoise points (wrong plot written?)")
        return
    names = ["total"] + vecs
    worst, wn, parts = 0.0, "", []
    for j, nm in enumerate(names):
        w = max(abs(b[i][1 + j] - a[i][1 + j]) / a[i][1] for i in range(len(a)))   # relative to the total (no tiny-part blow-up)
        parts.append(f"{nm} {w:.1e}")
        if w > worst:
            worst, wn = w, nm
    share = ", ".join(f"{nm} {100 * (a[0][2 + j] / a[0][1]) ** 2:.1f}%" for j, nm in enumerate(vecs))
    state = "pass" if worst <= bound else "FAIL"
    if red:
        state = "pass" if worst > bound else "FAIL"
    if known:
        state = "info"
        parts.append("the defect is still present" if worst > bound else "NOW WITHIN THE BOUND: the defect looks fixed, update the test")
    record(name, state,
           f"worst {worst:.2e} ({wn}, relative to the total; bound {bound:g}); {', '.join(parts)}; power at 1 kHz: {share}; "
           f"total 1k/100k {a[0][1]:.4e}/{a[-1][1]:.4e} V/rtHz")


def check_static_fd(ver, tnms=(1, 2), bound=1e-5, env=None, red=None):
    """Card #729: frequency-domain .pnoise (default solver; it runs every BSIM4 version, the td engine only bsim4 4.8) total vs plain
    .noise on the static deck, for one BSIM4 version (4.5 -> bsim4v5, 4.6.5 -> bsim4v6, 4.7 -> bsim4v7, 4.8.2 -> bsim4), each tnoimod in
    tnms, no workaround card.  The rs/rd part is ~0.2 % of the total power at 100 kHz (flicker and Rd dominate), far above the bound."""
    name = f"static fd: BSIM4 version {ver}, tnoimod {'/'.join(map(str, tnms))}, frequency-domain .pnoise total vs plain .noise [no workaround card, #729]"
    if red:
        name = f"[RED PROOF: {red}, must FAIL the check] " + name
    worst, parts = 0.0, []
    try:
        for tnm in tnms:
            lines = [l.replace("version=4.8.2", f"version={ver}") for l in set_params(fixture_lines("static_bsim4.cir"), tnm=tnm)]
            body = "\n".join(lines)
            outdir = os.path.join(DATA, f"staticfd_v{ver}_tnm{tnm}" + ("_red729" if red else ""))
            os.makedirs(outdir, exist_ok=True)
            opts = ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n"
            ctl = ".control\nset num_threads=1\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
            nz = os.path.join(outdir, "nz.cir")
            with open(nz, "w") as fh:
                fh.write(body.replace("Vg gs 0 DC 1.0", "Vg gs 0 DC 1.0 AC 1") + "\n" + opts + ".noise v(d) Vg dec 1 1k 100k 1\n" + ctl
                         + "setplot noise1\nwrdata nz.txt onoise_spectrum\nquit\n.endc\n.end\n")
            pn = os.path.join(outdir, "pn.cir")
            with open(pn, "w") as fh:
                fh.write(body + "\n" + opts + ".pnoise 1k 0 clkgnd 64 3 8 1e-6 d Vg dec 1 1k 100k\n" + ctl
                         + "wrdata pn.txt onoise_spectrum\nquit\n.endc\n.end\n")
            for d, f in ((nz, "nz.txt"), (pn, "pn.txt")):
                p = os.path.join(outdir, f)
                if os.path.exists(p):
                    os.remove(p)
                rc, out, wall = run_ngspice(d, env=env)
                if rc != 0 or "REFUSED" in out or "skipped" in out:
                    raise CouldNotRun(f"v{ver} tnoimod {tnm} {os.path.basename(d)} exit {rc}: {out[-1500:]}")
            a = read_cols(os.path.join(outdir, "nz.txt"), 1)
            b = read_cols(os.path.join(outdir, "pn.txt"), 1)
            if len(a) != len(b) or len(a) != 3:
                raise CouldNotRun(f"v{ver} tnoimod {tnm}: {len(a)} .noise points vs {len(b)} .pnoise points, expected 3")
            w = max(abs(b[i][1] - a[i][1]) / a[i][1] for i in range(len(a)))
            parts.append(f"tnoimod {tnm}: worst {w:.1e} (100k .noise {a[-1][1]:.6e}, .pnoise {b[-1][1]:.6e} V/rtHz)")
            worst = max(worst, w)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    ok = worst <= bound
    state = ("pass" if not ok else "FAIL") if red else ("pass" if ok else "FAIL")
    record(name, state, f"worst {worst:.2e} relative to the total (bound {bound:g}); " + "; ".join(parts))


def main():
    if not os.path.exists(NGSPICE):
        print(f"COULD-NOT-RUN: {NGSPICE} not found")
        return 2
    check_static(0)
    check_static(1, workaround=True)
    check_static(1)                                   # card #729: no workaround card
    check_static(2)
    for ver in ("4.5", "4.6.5", "4.7", "4.8.2"):      # card #729: every copy of the BSIM4 setup
        check_static_fd(ver)
    check_ensemble("cs", "res")
    check_ensemble("cs", "mos")
    if ("cs", "mos", 1) in ENS:
        check_ensemble("cs", "mos", tnm=1, workaround=True)
        check_ensemble("cs", "mos", tnm=1)            # card #729: no workaround card
    check_ensemble("th", "res")
    check_ensemble("th", "mos")
    check_ensemble("cs", "mos", env={"PAC_NOISE_STALE_STATE": "1"}, red="PAC_NOISE_STALE_STATE=1")
    check_ensemble("cs", "mos", red="ntnoi=1.25 in the deck", ntn_on=1.25)
    check_static(1, env=NO729, red="PSS_SALT_NO729=1")
    check_static_fd("4.8.2", tnms=(1,), env=NO729, red="PSS_SALT_NO729=1")
    if ("cs", "mos", 1) in ENS:
        check_ensemble("cs", "mos", tnm=1, env=NO729, red="PSS_SALT_NO729=1")
    width = max(len(n) for n, _, _ in results)
    for n, s, d in results:
        print(f"{s:14s} {n:{width}s}  {d}")
    npass = sum(1 for r in results if r[1] == "pass")
    ninfo = sum(1 for r in results if r[1] == "info")
    print(f"\n{npass}/{len(results) - ninfo} pass ({ninfo} info, not counted)")
    if any(r[1] == "FAIL" for r in results):
        return 1
    if any(r[1] == "COULD-NOT-RUN" for r in results):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
