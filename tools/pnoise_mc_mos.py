#!/usr/bin/env python3
"""Card #556: independent transient-noise reference for `.pnoise solver=td` on a MOS switching converter
(tests/fixtures/pnoise_m1buck.cir: MOS1 high-side switch + freewheel diode buck, fs 100 kHz).

Every device noise source is replaced by an explicit TRNOISE-driven current source with the INSTANTANEOUS density the
device model gives at the bias the waveform passes through (MOS1 channel thermal, 4kT*(2/3)*gm(vgs,vds) as
src/spicelib/devices/mos1/mos1noi.c for the default nlev = 2; diode shot 2q|Id|; resistors 4kT/R, noisy=0), written as B-sources of
the node voltages: I = sqrt(S(t)/(2 TS)) * V(nz), V(nz) = unit white TRNOISE (one-sided PSD 2 TS sinc^4).  The sources are scaled by
a common factor K (amplitude) so the output noise is far above the boxcar/difference rounding floor (2 nV/rtHz, see
pnoise_mc.py); the circuit is linear in it (output noise rms << the 5 V swing) and the engine value is scaled by K.
Two copies of the circuit (noisy / noiseless) share every time step and the difference is the noise response (pnoise_mc.py).

    NGSPICE=<prefix>/bin/ngspice-aug python3 tools/pnoise_mc_mos.py [nseeds] [K] [TS_ns] [m1,d1,rl,resr,rload]
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pnoise_mc as mc  # noqa: E402

K_BOLTZ, TK, QE = mc.K_BOLTZ, mc.T_KELVIN, 1.602176634e-19
FIX = "pnoise_m1buck.cir"
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/556/mc_m1buck")
mc.DATA_ROOT = DATA


ALL_SRC = ("m1", "d1", "rl", "resr", "rload")
VEC = {"m1": "onoise_m1_id", "d1": "onoise_d1_id", "rl": "onoise_rl", "resr": "onoise_resr", "rload": "onoise_rload"}


def core_lines(ts, K, only=ALL_SRC):
    """Subckt body: ports vin clk out, params na."""
    beta = 2e-3 * 100.0 / 1.0           # kp * W / L
    vto = 1.0
    kT4 = 4 * K_BOLTZ * TK
    sc = f"{{na*{K}}}"
    def src(name):
        return f"V{name} {name} 0 DC 0 TRNOISE({sc} {ts:.9e} 0 0)\n"
    amp = f"1/sqrt(2*{ts:.9e})"
    s = ""
    s += "M1 vin clk sw 0 nsw w=100u l=1u\n"
    s += ".model nsw nmos level=1 vto=1 kp=2e-3 gamma=0 lambda=0 phi=0.65 tox=5e-8 cgso=1e-10 cgdo=1e-10 cgbo=0 cj=0 kf=0\n"
    s += ".model dm d is=1e-9 n=1 rs=0 cj=0 tt=0\n"
    s += "Vsd 0 da 0\nD1 da sw dm\n"
    s += "Csw sw 0 200p\nL1 sw lx 22u\nRl lx out 50m noisy=0\nCout out cx 10u\nResr cx 0 100m noisy=0\nRload out 0 2.2 noisy=0\n"
    # M1 channel thermal, mos1noi.c with the default nlev = 2 (< 3): Sid = 2/3 |gm|, gm = beta min(vds, vgst) (lambda = gamma = 0)
    vgst = "(V(clk)-V(sw)-%g)" % vto
    vds = "max(V(vin)-V(sw),0)"
    sid = f"((2/3)*{beta}*min({vds},max({vgst},0)))"
    if "m1" in only:
        s += src("nz1")
        s += f"Bm1 vin sw I = {amp}*sqrt({kT4}*{sid}+1e-26)*V(nz1)\n"
    # diode shot, DIOnoise: 2 q |Id|
    if "d1" in only:
        s += src("nz2")
        s += f"Bd1 da sw I = {amp}*sqrt(2*{QE}*(abs(i(Vsd))+1e-6))*V(nz2)\n"
    # resistors: 4kT/R in parallel
    for i, (nm, n1, n2, r) in enumerate((("rl", "lx", "out", 50e-3), ("resr", "cx", "0", 100e-3), ("rload", "out", "0", 2.2))):
        if nm in only:
            s += src(f"nz{3+i}")
            s += f"B{3+i} {n1} {n2} I = {amp}*{math.sqrt(kT4 / r):.9e}*V(nz{3+i})\n"
    return s


def m1buck_mc_deck(outdir, seed, ts, duration, dt, K=1000.0, only=ALL_SRC, **_):
    core = core_lines(ts, K, only)
    top = (".options reltol=1e-4\nVin vin 0 DC 5\nVclk clk 0 PULSE(0 10 0 100n 100n 4.9e-6 10u)\n")
    return mc._deck(outdir, seed, ts, duration, dt, top, core, "out", "m1buck")


def run(nseeds=16, K=300.0, TS=20e-9, only=ALL_SRC):
    """Returns pnoise_mc.compare()'s (rows, (lo, hi sigma in dB), nseg); rows = (freq, MC, ref, d dB, d/sigma)."""
    mc.DATA_ROOT = DATA
    dt, nfft, dur, drop = 200e-9, 8192, 20e-3, 2e-3
    body = "\n".join(mc.body_lines(FIX))
    ca = ("100e3 300e-6 gh 4000 1 50 1e-6",)
    def ref(fl):
        if only == ALL_SRC:
            r = mc.pnoise_reference(body, os.path.join(DATA, "pnoise_ref"), "out", "Vin", ca, fl)
            return {f: K * v for f, v in r.items()}
        # a subset of the sources: sum of squares of the `.pnoise dev=1` per-device vectors
        outdir = os.path.join(DATA, "pnoise_ref_dev")
        os.makedirs(outdir, exist_ok=True)
        deck = os.path.join(outdir, "pn.cir")
        vecs = [VEC[o] for o in only]
        with open(deck, "w") as fh:
            fh.write(body + "\n")
            fh.write(f'.pnoise {ca[0]} out Vin lin {len(fl)} {fl[0]:.9g} {fl[-1]:.9g} method="shoot" solver="td" dev=1\n')
            fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nwrdata pn.txt " + " ".join(vecs) + "\nquit\n.endc\n.end\n")
        rc, o, wall = mc.run_ngspice(deck, timeout=600)
        if rc != 0:
            raise mc.CouldNotRun(o[-1500:])
        cols = [[] for _ in vecs]
        with open(os.path.join(outdir, "pn.txt")) as fh:
            next(fh)
            for row in fh:
                v = row.split()
                for j in range(len(vecs)):
                    cols[j].append(float(v[1 + j]))
        return {f: K * math.sqrt(sum(c[i] ** 2 for c in cols)) for i, f in enumerate(fl)}
    targets = [4000, 10000, 16000, 22000, 28000, 35000, 43000]
    tagn = "" if only == ALL_SRC else "_" + "+".join(only)
    return mc.compare(f"m1buck_K{K:g}_ts{TS*1e9:g}{tagn}", m1buck_mc_deck, body, "out", ca, ref, nseeds, TS, dur, dt, nfft, drop,
                      targets, ".pnoise*K", K=K, only=only)


def main():
    nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    K = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0
    TS = float(sys.argv[3]) * 1e-9 if len(sys.argv) > 3 else 20e-9
    only = tuple(sys.argv[4].split(",")) if len(sys.argv) > 4 else ALL_SRC
    run(nseeds, K, TS, only)


if __name__ == "__main__":
    main()
