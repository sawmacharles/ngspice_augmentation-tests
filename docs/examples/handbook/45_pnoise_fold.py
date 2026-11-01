#!/usr/bin/env python3
"""45_pnoise_fold.py -- handbook example 45 (card #560): `.pnoise` flicker (1/f) noise folded through a fast periodic modulation (card #556, patch 0110).

    python3 45_pnoise_fold.py NGSPICE OUTDIR

45_pnoise_fold_chop.cir (a resistor's 1/f and thermal noise, chopped by a +-1 square wave of 100 kHz with 1 us edges) is run at 1 / 30 / 90 / 150 kHz,
(a) with the default evaluation, (b) with PNOISE_FLICKER_OLD=1 (the pre-0110 evaluation: the density of every sideband read at the OUTPUT frequency f).
The reference is a closed form that uses no engine output: out = x(t) sq(t), x stationary with S_x(f) = A/f, sq(t) a trapezoid with exact Fourier
coefficients c_q (the second derivative of a piecewise-linear function is a sum of Dirac deltas at the corners), so S_out(f) = Sum_q |c_q|^2 S_x(|f - q fs|)
(thermal: S_x = 4 k T R, flicker: S_x = R^2 kf I^2 / f).  Output lines:
    FOLD <f> Hz  flicker engine <v> exact <v> (<rel>)  old <v> (<rel>)   thermal engine <v> exact <v> (<rel>)
"""
import cmath
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
K_B, T_K = 1.380649e-23, 300.15
FS, TR, R, I, KF = 1e5, 1e-6, 1e3, 1e-3, 4e-12


def trap_coeffs(T, tr, pw, lo=-1.0, hi=1.0, Q=4000):
    ts = [0.0, tr, tr + pw, 2 * tr + pw]
    sl = [(hi - lo) / tr, 0.0, -(hi - lo) / tr, 0.0]
    slp = [0.0, (hi - lo) / tr, 0.0, -(hi - lo) / tr]
    ds = [a - b for a, b in zip(sl, slp)]
    cq = []
    for q in range(1, Q + 1):
        th = 2 * math.pi * q / T
        cq.append(-sum(d * cmath.exp(-1j * th * t) for t, d in zip(ts, ds)) / (th * th * T))
    mean = (lo * (T - 2 * tr - pw) + hi * pw + 0.5 * (lo + hi) * 2 * tr) / T
    return mean, cq


def s_out(f, mean, cq, S):
    tot = mean ** 2 * S(abs(f))
    for i, c in enumerate(cq):
        q = i + 1
        tot += abs(c) ** 2 * (S(abs(f - q * FS)) + S(abs(f + q * FS)))
    return tot


def run(ng, out, f, env=None):
    t = open(os.path.join(HERE, "45_pnoise_fold_chop.cir")).read()
    t = t.replace("lin 1 1k 1k", "lin 1 %g %g" % (f, f))
    tag = "f%g%s" % (f, "_old" if env else "")
    d = os.path.join(out, tag)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "deck.cir")
    open(p, "w").write(t)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", p],
                       cwd=d, env=e, capture_output=True, text=True, timeout=600)
    o = (r.stdout + r.stderr).replace("\r", "\n")
    open(os.path.join(d, "log.txt"), "w").write(o)
    rows = [ln.split() for ln in open(os.path.join(d, "pnoise_fold.txt")).read().splitlines()[1:] if ln.strip()]
    v = [float(x) for x in rows[0]]
    return {"thermal": v[1], "flicker": v[2], "total": v[3], "spectrum": v[4]}


def main(ng, out):
    T = 1 / FS
    mean, cq = trap_coeffs(T, TR, 0.5 * T - TR)
    th = lambda x: R * 4 * K_B * T_K
    fl = lambda x: R * R * KF * I ** 2 / x
    for f in (1e3, 30e3, 90e3, 150e3):
        new, old = run(ng, out, f), run(ng, out, f, {"PNOISE_FLICKER_OLD": "1"})
        rf, rt = math.sqrt(s_out(f, mean, cq, fl)), math.sqrt(s_out(f, mean, cq, th))
        ratio = old["flicker"] / rf
        print("FOLD %6.0f Hz  flicker engine %.4e exact %.4e (%+.1e)  old %.4e (%s)   thermal engine %.4e exact %.4e (%+.1e)"
              % (f, new["flicker"], rf, new["flicker"] / rf - 1, old["flicker"],
                 ("x%.2f" % ratio) if ratio > 2 else ("%+.0f %%" % (100 * (ratio - 1))), new["thermal"], rt, new["thermal"] / rt - 1), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
