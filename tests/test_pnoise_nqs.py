#!/usr/bin/env python3
"""Card #715: td `.pnoise` (method="shoot" solver="td") on NON-QUASI-STATIC MOS channels (BSIM4 trnqsmod=1, BSIM3v3.3 nqsmod=1), refused
up to s146 ("not validated"), accepted since #715 (docs/pss-pac.md, "Non-quasi-static channels in td .pnoise (card #715, 2026-10-19)").

The NQS models put the channel-charge relaxation in an extra circuit node (qNode) with no noise source of its own: BSIM4noise and BSIM3noise
inject at the drain/source (and gate/bulk) pairs with densities computed from the quasi-static bias.  What NQS changes in .pnoise is the
TRANSFER from those pairs to the output, i.e. the relaxation node as a state of the (transposed) td recursion.

Fixtures tests/fixtures/pnoise715/:
  static_nqs_bsim3.cir  BSIM3v3.3 nqsmod=1 common-source stage at a DC bias (W/L 40u/20u, Rs 1k, Rd 20k, Cl 50f), slow clock elsewhere.
  static_nqs_bsim4.cir  the same with BSIM4 4.8.2 trnqsmod=1 (the card-#559 model).
  lo_nqs_bsim4.cir      the BSIM4 stage with its gate driven by a 1 MHz LO (0.9 V +- 0.35 V through Rs 1k): an LPTV circuit whose device goes
                        from near cut-off into triode every period; analysed at 4.3 .. 42.7 MHz (omega*tau 3..30).

Checks (N/M pass printed at the end; exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN):
  static bsim3   `.pnoise dev=1` total and every generator (M1 id/rd/rs/1overf, Rs, Rd) == plain `.noise`, 0.1 .. 10 MHz, 2e-5 of the total.
                 BSIM3's .noise stamps the linearised TRANSIENT NQS model (b3acld.c), so this reference contains the relaxation state and is
                 independent of the td engine.
  static bsim4   every generator == its `.noise` value x |td .pac transfer| / |.ac transfer| from the same node pair (2e-5 of the total).
                 BSIM4's .noise/.ac are quasi-static with trnqsmod (b4acld.c isolates the charge node), so the density comes from .noise
                 and the NQS transfer from the forward td .pac (validated against a transient FRA by card #559): checks the adjoint
                 noise engine's handling of the relaxation node, not the transfer itself.
  lo res / mos   per 9-bin band (4.3 .. 42.7 MHz, sideband-0 time-averaged PSD), .pnoise vs the PINNED transient-noise ensemble,
                 |d| <= 0.5 dB per band and |mean d| <= 0.2 dB (band 1-sigma 0.113 dB, mean 1-sigma 0.043 dB: 4.4 / 4.7 sigma):
                   res  BSIM4 channel noise off (ntnoi=0), Rs and Rd thermal noise (explicit white trnoise currents in the ensemble)
                   mos  resistors noiseless, BSIM4 channel thermal noise: in the ensemble a B-source sqrt(S(t)) w(t), S(t) the model's own
                        drain-current noise density at the orbit's terminal voltages (plain .noise per time point, card #708 method)
  red            the quasi-static deck (nqsmod=0 / trnqsmod=0) must FAIL the static bsim3 check and the lo res check (the checks see the
                 relaxation state); the BSIM4 static .pnoise must differ from the quasi-static .noise (the Rs generator by 0.13 of the total at 10 MHz);
                 PAC_NQS_PNOISE_REFUSE=1 must restore the refusal; PAC_NQS_REFUSE=1 (pre-#559) must refuse too.

The ensembles are not re-run by the test (16 seeds x 20 us each, 2-8 min); they were produced by
    MC_BSIM4_DATA=<dir> MC_OPTIONS="reltol=1e-7 vntol=1e-10" NGSPICE=<w715>/bin/ngspice python3 tools/pnoise_mc_bsim4.py lo <res|mos> 16 <K> 0.05
(mos also with MC_SEED_OFFSET=100; K = 30 (res) / 100 (mos) amplitude scale, TS = 0.05 ns trnoise step; Hann/Welch 8192 x 0.2 ns boxcar samples, 320 segments, 9-bin bands of 610.35 kHz; the
chain calibrated on R 10k / C 1p against the analytic 4kTR/(1+(f/fc)^2) at this time base: worst band 0.172 dB, mean -0.030 dB).  The
ensemble transient needs a tight Newton tolerance with the NQS model: at the default reltol 1e-3 the noisy-minus-noiseless difference
carries a low-frequency error that is not noise (+13.6 dB in the 4.3 MHz band at K = 100, +8.5 dB at K = 30, +3.1 dB at K = 10; none with
the quasi-static model); at reltol 1e-7 / vntol 1e-10 (1e-8 / 1e-12 gives the same to 0.05 dB) what is left below ~2.5 MHz scales with K^2
(0.61 MHz bin: +1.55 dB at K = 100, +0.37 dB at K = 30), i.e. a second-order (rectified) part of the large-signal ensemble that the linear
.pnoise does not contain; the res pin is therefore the K = 30 run.  Logs in ~/Work/data/ngspice_augmentation/715/mc.

    NGSPICE=<prefix>/bin/ngspice python3 tests/test_pnoise_nqs.py      (default: the s152 release)
Outputs under ~/Work/data/ngspice_augmentation/715/test, or $PNOISE_NQS_DATA (never the repository).
"""
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures", "pnoise715")
NGSPICE = os.environ.get("NGSPICE", os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice"))
DATA = os.environ.get("PNOISE_NQS_DATA", os.path.expanduser("~/Work/data/ngspice_augmentation/715/test"))

BINW = 1.0 / (8192 * 0.2e-9)          # Welch bin of the lo ensembles, 610351.6 Hz
NB = 9
TARGETS = [4e6, 10e6, 16e6, 22e6, 28e6, 35e6, 43e6]
BAND_DB, MEAN_DB = 0.5, 0.2
LO_HEAD = "1meg 20u out 2000 1 50 1e-6"
STATIC_BOUND = 2e-5

# (part) -> (sigma_dB of one band, {band centre Hz: ensemble V/rtHz}); band centres 4.2725, 9.7656, 15.869, 21.973, 28.076, 34.790, 42.725 MHz
# (keyed by TARGETS); 16 seeds x 20 us, K = 100, TS = 0.05 ns, ensemble transient at reltol 1e-7 vntol 1e-10
_T = TARGETS
ENS = {
    # Rs + Rd, K = 30 (see the docstring: at K = 100 a second-order low-frequency part adds +0.06 dB to the 4.3 MHz band; bands from 9.8 MHz
    # agree with K = 100 to 0.006 dB), wall 648 s (2 ngspice); trnqsmod=0 ensemble (reltol 1e-3, K = 100): -0.189 .. +0.149 dB from its .pnoise
    "res": (0.113, dict(zip(_T, [2.17227e-08, 2.05194e-08, 1.98483e-08, 1.92464e-08, 1.87429e-08, 1.83539e-08, 1.80641e-08]))),
    # channel noise, S(t) table 2002 orbit points, <S_id> = 1.84013e-24 A^2/Hz (orbit periodicity 1.8e-7 V); two independent 16-seed
    # ensembles (seeds 1..16: .pnoise - ensemble worst 0.213 dB, mean +0.083; seeds 101..116: worst 0.197, mean +0.033; K = 30 with seeds
    # 1..16 repeats the first to 0.002 dB, so K = 100 is linear), pinned as their power mean: 32 seeds, band 1-sigma 0.080 dB; wall 490 s each
    "mos": (0.080, dict(zip(_T, [2.61013e-08, 2.55825e-08, 2.57774e-08, 2.54158e-08, 2.49785e-08, 2.45845e-08, 2.44627e-08]))),
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
    out = []
    for l in lines:
        if l.lower().startswith(".param"):
            for k, v in kw.items():
                l = re.sub(rf"\b{k}=\S+", f"{k}={v}", l)
        out.append(l)
    return out


def run_ngspice(deck, env=None, timeout=900):
    e = dict(os.environ)
    for k in ("PAC_ALLOW_NQS", "PAC_NQS_REFUSE", "PAC_NQS_PNOISE_REFUSE"):   # only what a check asks for
        e.pop(k, None)
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


OPTS = ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n"
CTL = ".control\nset num_threads=1\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
SWEEP = "dec 5 100k 10meg"   # ends on a decade: .ac then uses the same grid as .noise/.pac/.pnoise (with 20meg it does not)
PNOISE_STATIC = '.pnoise 100k 0 clkgnd 2000 1 8 1e-6 out Vin {sweep} method="shoot" solver="td" dev=1 pacpoints=40000\n'


def run_deck(outdir, name, text, outfile, env=None, need=None):
    deck = os.path.join(outdir, name + ".cir")
    p = os.path.join(outdir, outfile)
    if os.path.exists(p):
        os.remove(p)
    with open(deck, "w") as fh:
        fh.write(text)
    rc, out, wall = run_ngspice(deck, env=env)
    if rc != 0 or "REFUSED" in out or (need and need not in out):
        raise CouldNotRun(f"{name} exit {rc}: " + (" ".join(l for l in out.splitlines() if "REFUSED" in l)[:400] or out[-1200:]))
    return out


def static_runs(fix, tag, vecs, nq_nz=1, nq_pn=1, env=None):
    """.noise (nq_nz) and td .pnoise dev=1 (nq_pn) of a static fixture -> (rows .noise, rows .pnoise); columns f, total, vecs..."""
    outdir = os.path.join(DATA, tag)
    os.makedirs(outdir, exist_ok=True)
    nz = "\n".join(set_params(fixture_lines(fix), nq=nq_nz)) + "\n" + OPTS + f".noise v(out) Vin {SWEEP} 1\n" + CTL \
        + "setplot noise1\nwrdata nz.txt onoise_spectrum " + " ".join(vecs) + "\nquit\n.endc\n.end\n"
    pn = "\n".join(set_params(fixture_lines(fix), nq=nq_pn)) + "\n" + OPTS + PNOISE_STATIC.format(sweep=SWEEP) + CTL \
        + "wrdata pn.txt onoise_spectrum " + " ".join(vecs) + "\nquit\n.endc\n.end\n"
    run_deck(outdir, "nz", nz, "nz.txt")
    run_deck(outdir, "pn", pn, "pn.txt", env=env, need="PNOISE-TD: done")
    a = read_cols(os.path.join(outdir, "nz.txt"), 1 + len(vecs))
    b = read_cols(os.path.join(outdir, "pn.txt"), 1 + len(vecs))
    if len(a) != len(b):
        raise CouldNotRun(f"{len(a)} .noise points vs {len(b)} .pnoise points (wrong plot written?)")
    return a, b


def compare_rows(a, b, names):
    """worst |b - a| relative to the total of a, per column; -> worst, its name, text"""
    worst, wn, parts = 0.0, "", []
    for j, nm in enumerate(names):
        w = max(abs(b[i][1 + j] - a[i][1 + j]) / a[i][1] for i in range(len(a)))
        parts.append(f"{nm} {w:.1e}")
        if w > worst:
            worst, wn = w, nm
    return worst, wn, ", ".join(parts)


def check_static_bsim3(nq_pn=1, red=None):
    name = "static bsim3: BSIM3v3.3 nqsmod=1 stage at a DC bias, td .pnoise dev=1 vs plain .noise (which carries the transient NQS model)"
    if red:
        name = f"[RED PROOF: {red}, must FAIL the check] " + name
    vecs = ["onoise.m1.id", "onoise.m1.rd", "onoise.m1.rs", "onoise.m1.1overf", "onoise_rs", "onoise_rd"]
    try:
        a, b = static_runs("static_nqs_bsim3.cir", f"static_b3_pn{nq_pn}", vecs, nq_nz=1, nq_pn=nq_pn)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    worst, wn, parts = compare_rows(a, b, ["total"] + vecs)
    ok = worst <= STATIC_BOUND
    state = ("pass" if not ok else "FAIL") if red else ("pass" if ok else "FAIL")
    k10 = min(range(len(a)), key=lambda i: abs(a[i][0] - 1e7))
    record(name, state, f"worst {worst:.2e} ({wn}, relative to the total; bound {STATIC_BOUND:g}); {parts}; "
                        f"total at {a[k10][0] / 1e6:.3g} MHz: .noise {a[k10][1]:.5e}, .pnoise {b[k10][1]:.5e} V/rtHz")


def bsim4_transfer(outdir, tag, inj, analysis, nq=1):
    """|v(out)| for a unit AC current between the node pair `inj` (Vin AC removed): analysis 'ac' (quasi-static) or 'pac' (td)."""
    lines = [l.replace("Vin in 0 DC 0.9 AC 1", "Vin in 0 DC 0.9") for l in set_params(fixture_lines("static_nqs_bsim4.cir"), nq=nq)]
    card = (f".ac {SWEEP}\n" if analysis == "ac" else
            f'.pac 100k 0 clkgnd 2000 1 8 1e-6 {SWEEP} method="shoot" solver="td" pacpoints=40000\n')
    fn = f"h_{tag}_{analysis}.txt"
    txt = "\n".join(lines) + f"\nIinj {inj[0]} {inj[1]} DC 0 AC 1\n" + OPTS + card + CTL + f"let hm = mag(v(out))\nwrdata {fn} hm\nquit\n.endc\n.end\n"
    run_deck(outdir, f"h_{tag}_{analysis}", txt, fn)
    return [(r[0], r[1]) for r in read_cols(os.path.join(outdir, fn), 1)]


def check_static_bsim4():
    name = ("static bsim4: BSIM4 trnqsmod=1 stage at a DC bias, td .pnoise dev=1 per generator vs .noise density x |td .pac| / |.ac| "
            "transfer from the generator's node pair")
    red_name = "[RED PROOF: the quasi-static .noise itself, must differ by more than 1e-2 of the total] " + name
    vecs = ["onoise.m1.id", "onoise_rs", "onoise_rd"]
    inj = {"onoise.m1.id": ("out", "0"), "onoise_rd": ("out", "0"), "onoise_rs": ("g", "in")}   # source at ground, rdsmod 0; vdd is AC ground
    try:
        a, b = static_runs("static_nqs_bsim4.cir", "static_b4", vecs)
        outdir = os.path.join(DATA, "static_b4")
        H = {}
        for pair in set(inj.values()):
            t = "_".join(pair)
            H[pair] = (bsim4_transfer(outdir, t, pair, "pac"), bsim4_transfer(outdir, t, pair, "ac"))
            for h in H[pair]:
                if len(h) != len(a) or any(abs(h[i][0] / a[i][0] - 1) > 1e-9 for i in range(len(a))):
                    raise CouldNotRun(f"transfer {pair}: frequency grid differs from .noise ({len(h)} vs {len(a)} points)")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        record(red_name, "COULD-NOT-RUN", str(ex))
        return
    pred = []
    for i in range(len(a)):
        row = [a[i][0]]
        gens = [a[i][2 + j] * H[inj[v]][0][i][1] / H[inj[v]][1][i][1] for j, v in enumerate(vecs)]
        row += [math.sqrt(sum(g * g for g in gens))] + gens
        pred.append(row)
    worst, wn, parts = compare_rows(pred, b, ["total"] + vecs)
    record(name, "pass" if worst <= STATIC_BOUND else "FAIL",
           f"worst {worst:.2e} ({wn}, relative to the predicted total; bound {STATIC_BOUND:g}); {parts}")
    wq, wqn, _ = compare_rows(a, b, ["total"] + vecs)
    record(red_name, "pass" if wq > 1e-2 else "FAIL", f"td .pnoise vs quasi-static .noise: worst {wq:.2e} ({wqn}) (must exceed 1e-2)")


def pnoise_bands(part, nq=1, env=None, tag=""):
    """lo deck: .pnoise (solver=td) over every Welch bin of the 7 bands in ONE run -> {centre: band-power-mean V/rtHz}."""
    lines = set_params(fixture_lines("lo_nqs_bsim4.cir"), ntn=(0 if part == "res" else 1), nq=nq)
    if part == "mos":
        lines = [l + " noisy=0" if l[:1] in "Rr" else l for l in lines]
    kfirst = round(band_bins(TARGETS[0])[0] / BINW)
    klast = round(band_bins(TARGETS[-1])[-1] / BINW)
    outdir = os.path.join(DATA, f"lo_{part}_nq{nq}{tag}")
    os.makedirs(outdir, exist_ok=True)
    txt = ("\n".join(lines) + "\n" + f'.pnoise {LO_HEAD} out Vin lin {klast - kfirst + 1} {kfirst * BINW:.10g} {klast * BINW:.10g} '
           'method="shoot" solver="td"\n' + CTL + "wrdata pn.txt onoise_spectrum\nquit\n.endc\n.end\n")
    t0 = time.time()
    run_deck(outdir, "pn", txt, "pn.txt", env=env, need="PNOISE-TD: done")
    rows = read_cols(os.path.join(outdir, "pn.txt"), 1)
    if len(rows) != klast - kfirst + 1:
        raise CouldNotRun(f"{len(rows)} rows, expected {klast - kfirst + 1} (wrong plot written?)")
    byk = {round(r[0] / BINW): r[1] for r in rows}
    return {t: math.sqrt(sum(byk[round(f / BINW)] ** 2 for f in band_bins(t)) / NB) for t in TARGETS}, time.time() - t0


def check_lo(part, nq=1, red=None):
    label = {"res": "Rs + Rd thermal noise, BSIM4 channel noise off", "mos": "BSIM4 channel thermal noise, resistors noiseless"}[part]
    name = f"lo {part}: 1 MHz LO-driven BSIM4 trnqsmod=1 stage, .pnoise vs transient-noise ensemble, {label}"
    if red:
        name = f"[RED PROOF: {red}, must FAIL the check] " + name
    sig, ens = ENS[part]
    if any(v is None for v in ens.values()):
        record(name, "COULD-NOT-RUN", "no pinned ensemble for this case")
        return
    try:
        pn, wall = pnoise_bands(part, nq=nq, tag=("_red" if red else ""))
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    ds = [20 * math.log10(pn[t] / ens[t]) for t in TARGETS]
    worst = max(abs(d) for d in ds)
    mean = sum(ds) / len(ds)
    ok = worst <= BAND_DB and abs(mean) <= MEAN_DB
    rows = ", ".join(f"{t / 1e6:g}M {d:+.3f}" for t, d in zip(TARGETS, ds))
    detail = (f".pnoise - ensemble (dB) {rows}; worst {worst:.3f}, mean {mean:+.3f} (limits {BAND_DB}/{MEAN_DB}; band 1-sigma "
              f"{sig:.3f} dB, mean 1-sigma {sig / math.sqrt(len(ds)):.3f}); 4.3M: {pn[TARGETS[0]]:.4e} vs {ens[TARGETS[0]]:.4e} V/rtHz; "
              f"wall {wall:.1f}s")
    record(name, ("pass" if not ok else "FAIL") if red else ("pass" if ok else "FAIL"), detail)


def check_refusal(env, label, must):
    name = f"guard: {label} restores the refusal of td .pnoise on the BSIM4 trnqsmod=1 stage"
    outdir = os.path.join(DATA, "refusal_" + "_".join(env))
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, "pn.cir")
    with open(deck, "w") as fh:
        fh.write("\n".join(fixture_lines("static_nqs_bsim4.cir")) + "\n" + PNOISE_STATIC.format(sweep="lin 1 1e7 1e7") + CTL
                 + "quit\n.endc\n.end\n")
    try:
        rc, out, wall = run_ngspice(deck, env=env)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    ref = [l for l in out.splitlines() if "REFUSED" in l]
    ok = bool(ref) and "PNOISE" in ref[0] and must in ref[0]
    record(name, "pass" if ok else "FAIL", (ref[0][:240] if ref else "no REFUSED line: " + out[-300:]))


def main():
    if not os.path.exists(NGSPICE):
        print(f"COULD-NOT-RUN: {NGSPICE} not found")
        return 2
    check_static_bsim3()
    check_static_bsim4()
    check_lo("res")
    check_lo("mos")
    check_static_bsim3(nq_pn=0, red="nqsmod=0 (quasi-static) deck in .pnoise vs the nqsmod=1 .noise")
    check_lo("res", nq=0, red="trnqsmod=0 (quasi-static) deck")
    check_refusal({"PAC_NQS_PNOISE_REFUSE": "1"}, "PAC_NQS_PNOISE_REFUSE=1", "PAC_NQS_PNOISE_REFUSE")
    check_refusal({"PAC_NQS_REFUSE": "1"}, "PAC_NQS_REFUSE=1 (pre-#559)", "PAC_NQS_REFUSE")
    width = max(len(n) for n, _, _ in results)
    for n, s, d in results:
        print(f"{s:14s} {n:{width}s}  {d}")
    npass = sum(1 for r in results if r[1] == "pass")
    print(f"\n{npass}/{len(results)} pass")
    if any(r[1] == "FAIL" for r in results):
        return 1
    if any(r[1] == "COULD-NOT-RUN" for r in results):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
