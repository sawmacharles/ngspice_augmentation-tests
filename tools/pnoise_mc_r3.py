#!/usr/bin/env python3
"""Card #556: transient-noise Monte Carlo for `.pnoise solver=td` on the CLOSED-LOOP switching converter R3 (tests/fixtures/
pac_r3_ctrl2out.cir: voltage-mode buck, fs 500 kHz, type-III compensation, comparator/PWM, smooth tanh switches).  Every noisy
device is a plain resistor (the switches, op-amp and comparator are noiseless behavioural sources), so each resistor is made
noisy=0 and a TRNOISE current source of one-sided PSD 4kT/R (times K^2, see below) is put in parallel with it, exactly the
pnoise_mc.py construction (two copies of the circuit in one transient share every time step; the difference is the noise
response; boxcar-integrated trace, Hann/Welch, 9-bin band averages, statistics as stated there).  K scales the injected amplitude
so the output noise (1.3e-8 .. 5e-8 V/rtHz) is far above the 2 nV/rtHz difference floor of that chain; the loop is linear in it
(output rms << ripple) and `.pnoise` is scaled by K.

    NGSPICE=<prefix>/bin/ngspice-aug python3 tools/pnoise_mc_r3.py [nseeds] [K] [TS_ns]
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pnoise_mc as mc  # noqa: E402

FIX = "pac_r3_ctrl2out.cir"
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/556/mc_r3")
mc.DATA_ROOT = DATA
K_B, TK = mc.K_BOLTZ, mc.T_KELVIN


def r3_mc_deck(outdir, seed, ts, duration, dt, K=100.0, **_):
    os.makedirs(outdir, exist_ok=True)
    glob, core, noise = [], [], []
    for l in mc.body_lines(FIX):
        ls = l.strip()
        if not ls or ls.startswith("*"):
            continue
        if ls.lower().startswith(".param") or ls.lower().startswith(".options"):
            glob.append(l)
            continue
        m = re.match(r"^(R\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(.*)$", ls)
        if m:
            name, n1, n2, val, rest = m.groups()
            core.append(f"{name} {n1} {n2} {val} noisy=0 {rest}")
            rexpr = val[1:-1] if val.startswith("{") and val.endswith("}") else val
            noise.append(f"I_nz_{name} {n1} {n2} DC 0 TRNOISE({{na*{K}*sqrt({2*K_B*TK/ts:.9e}/({rexpr}))}} {ts:.9e} 0 0)")
        elif ls.lower().startswith("vctrl"):
            core.append(re.sub(r"\s+AC\s+1", "", l))      # no AC source needed in a transient
        elif ls.startswith("."):
            continue
        else:
            core.append(l)
    deck = os.path.join(outdir, f"mc_seed{seed}.cir")
    with open(deck, "w") as fh:
        fh.write("* R3 closed-loop buck, trnoise replica (card #556)\n")
        fh.write("\n".join(glob) + "\n")
        fh.write(".subckt core vout params: na=1\n" + "\n".join(core) + "\n" + "\n".join(noise) + "\n.ends core\n")
        fh.write("Xn outn core na=1\nXd outd core na=0\n")
        # reltol 1e-4: at the default 1e-3 the transient's own truncation error shows in the noise response (+0.14 .. +0.28 dB mean, up to
        # 0.75 dB per band at K = 10..100 in the first runs, gone at 1e-4: mean +0.02 dB)
        fh.write(f".options reltol={os.environ.get('MC_RELTOL', '1e-4')}\n")
        fh.write(f".options maxstep={ts:.6g}\n.tran {dt:.9e} {duration:.9e}\n.control\n"
                 f"setseed {seed}\nrun\nlet iv = integ(v(outn)-v(outd))\nlinearize iv\n"
                 "set wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                 "wrdata mc_out.txt iv\nquit\n.endc\n.end\n")
    return deck, os.path.join(outdir, "mc_out.txt")


def run(nseeds=16, K=100.0, TS=10e-9):
    mc.DATA_ROOT = DATA
    dt, nfft, dur, drop = 200e-9, 8192, 20e-3, 2e-3
    body = "\n".join(mc.body_lines(FIX))
    ca = ("500e3 300e-6 vout 200 10 150 1e-6",)

    def ref(fl):
        r = mc.pnoise_reference(body, os.path.join(DATA, "pnoise_ref"), "vout", "Vin", ca, fl)
        return {f: K * v for f, v in r.items()}
    targets = [4000, 10000, 16000, 22000, 28000, 35000, 43000]
    return mc.compare(f"r3_K{K:g}_ts{TS*1e9:g}", r3_mc_deck, body, "vout", ca, ref, nseeds, TS, dur, dt, nfft, drop,
                      targets, ".pnoise*K", K=K)


def main():
    nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    K = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
    TS = float(sys.argv[3]) * 1e-9 if len(sys.argv) > 3 else 10e-9
    run(nseeds, K, TS)


if __name__ == "__main__":
    main()
