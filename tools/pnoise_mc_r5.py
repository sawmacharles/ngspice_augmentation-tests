#!/usr/bin/env python3
"""Card #557: transient-noise Monte Carlo for `.pnoise solver=td` on the FREE-RUNNING hysteretic buck R5 (tests/fixtures/r5_freerun_plain.cir
(, `trigger="q"`).  Same construction as pnoise_mc_r3.py (resistors made noisy=0 plus a parallel TRNOISE current source of PSD 4kT/R*K^2;
two copies of the converter in one transient, noisy and noiseless, boxcar-integrated difference, Hann/Welch, 9-bin band averages) except
that the converter is autonomous: both copies start from the same state and their switching phases separate by the phase diffusion of the
noise, which is exactly what `.pnoise` includes (the Floquet phase mode).  K keeps the noise far above the 2 nV/rtHz floor of the integ chain
while the linear (small-signal) assumption holds.

    NGSPICE=<prefix>/bin/ngspice python3 tools/pnoise_mc_r5.py [nseeds] [K] [TS_ns] [targets_kHz,...] [dur_ms]
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pnoise_mc as mc  # noqa: E402

FIX = "r5_freerun_plain.cir"
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/557/mc_r5")
mc.DATA_ROOT = DATA
K_B, TK = mc.K_BOLTZ, mc.T_KELVIN
CARD = '427350 450e-6 vout 200 10 45 1e-4'


def r5_body():
    return "\n".join(mc.body_lines(FIX))


def r5_mc_deck(outdir, seed, ts, duration, dt, K=100.0, **_):
    os.makedirs(outdir, exist_ok=True)
    glob, core, noise = [], [], []
    for l in r5_body().splitlines():
        ls = l.strip()
        if not ls or ls.startswith("*"):
            continue
        if ls.lower().startswith((".param", ".options")):
            glob.append(l)
            continue
        if ls.lower().startswith(".ic"):
            continue
        m = re.match(r"^(R\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(.*)$", ls)
        if m:
            name, n1, n2, val, rest = m.groups()
            core.append(f"{name} {n1} {n2} {val} noisy=0 {rest}")
            rexpr = val[1:-1] if val.startswith("{") and val.endswith("}") else val
            noise.append(f"I_nz_{name} {n1} {n2} DC 0 TRNOISE({{na*{K}*sqrt({2*K_B*TK/ts:.9e}/({rexpr}))}} {ts:.9e} 0 0)")
        elif ls.startswith("."):
            continue
        else:
            core.append(l)
    deck = os.path.join(outdir, f"mc_seed{seed}.cir")
    with open(deck, "w") as fh:
        fh.write("* R5 free-running buck, trnoise replica (card #557)\n")
        fh.write("\n".join(glob) + "\n")
        fh.write(".subckt core vout params: na=1\n" + "\n".join(core) + "\n" + "\n".join(noise) + "\n.ends core\n")
        fh.write("Xn outn core na=1\nXd outd core na=0\n")
        fh.write(".ic v(xn.q)=0.5 v(xd.q)=0.5\n")
        fh.write(f".options reltol={os.environ.get('MC_RELTOL', '1e-4')}\n.save v(outn) v(outd)\n")
        fh.write(f".options maxstep={ts:.6g}\n.tran {dt:.9e} {duration:.9e}\n.control\n"
                 f"setseed {seed}\nrun\nlet iv = integ(v(outn)-v(outd))\nlinearize iv\n"
                 "set wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                 "wrdata mc_out.txt iv\nquit\n.endc\n.end\n")
    return deck, os.path.join(outdir, "mc_out.txt")


def pnoise_ref_r5(freqs, K):
    os.makedirs(os.path.join(DATA, "pnoise_ref"), exist_ok=True)
    d = os.path.join(DATA, "pnoise_ref")
    deck = os.path.join(d, f"pn_{freqs[0]:.6g}.cir")
    with open(deck, "w") as fh:
        fh.write(r5_body() + "\n")
        fh.write(f'.pnoise {CARD} vout Vin lin {len(freqs)} {freqs[0]:.9g} {freqs[-1]:.9g} method="shoot" solver="td" '
                 'trigger="q" trigval=0.5 trigdir="rise"\n')
        fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nwrdata pnout.txt onoise_spectrum\nquit\n.endc\n.end\n")
    rc, o, wall = mc.run_ngspice(deck, timeout=1800)
    if rc != 0:
        raise mc.CouldNotRun(f".pnoise exit {rc}: {o[-1500:]}")
    v = mc._read_wrdata_real(os.path.join(d, "pnout.txt"))
    if len(v) != len(freqs):
        raise mc.CouldNotRun(f"expected {len(freqs)} points, got {len(v)}")
    return {f: K * x for f, x in zip(freqs, v)}


def run(nseeds=16, K=100.0, TS=10e-9, targets=(15000, 40000, 85000), dur=20e-3):
    mc.DATA_ROOT = DATA
    dt, nfft, drop = 200e-9, 8192, 2e-3
    ref = lambda fl: pnoise_ref_r5(fl, K)
    return mc.compare(f"r5_K{K:g}_ts{TS*1e9:g}", r5_mc_deck, r5_body(), "vout", None, ref, nseeds, TS, dur, dt, nfft, drop,
                      list(targets), ".pnoise*K", K=K)


def main():
    a = sys.argv
    nseeds = int(a[1]) if len(a) > 1 else 16
    K = float(a[2]) if len(a) > 2 else 100.0
    TS = float(a[3]) * 1e-9 if len(a) > 3 else 10e-9
    tg = [float(x) * 1e3 for x in a[4].split(",")] if len(a) > 4 else [15000, 40000, 85000]
    dur = float(a[5]) * 1e-3 if len(a) > 5 else 20e-3
    run(nseeds, K, TS, tg, dur)


if __name__ == "__main__":
    main()
