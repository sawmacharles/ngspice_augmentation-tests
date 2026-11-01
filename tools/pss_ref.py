#!/usr/bin/env python3
"""Independent references for validating `.pss method=shoot` (card #458
second pass, plan T4.1-T4.3). Two pieces, both used to check the C code in
ngspice-src-pss/src/spicelib/analysis/psseigen.c and pssshoot.c from
OUTSIDE that code, not by re-reading it:

  1. eig_reference(A): eigenvalues of a real n*n matrix via an algorithm
     that shares NO code path with psseigen.c's Hessenberg + double-shift
     QR ("hqr") solver -- Faddeev-LeVerrier for the characteristic
     polynomial (exact, trace-based, no eigendecomposition involved) then
     Durand-Kerner simultaneous iteration for its roots (a different
     numerical method again: polynomial root-finding, not matrix
     iteration). Used (card #458 report) to validate psseigen.c on six
     matrices before it was wired into pssshoot.c; the transcript is in
     the card #458 report, not reproduced here as a fixture since the
     validation was interactive (stdout comparison), not a repo test.

  2. r2_floquet_analytic(...): the exact Floquet multiplier of the R2
     switched-RC circuit (tests/fixtures/pss_r2_shoot.cir / refs/
     r2_switched_rc.cir) from the deck's own R/C/duty values, no
     simulation involved -- R2 is clock-driven (not state-driven)
     switching between two fixed linear-RC time constants each period,
     so the exact one-period state-transition multiplier is the product
     of the two intervals' pure exponential decay:
         lambda = exp(-(T_on/tau_on + T_off/tau_off))
         tau_on  = 1 / ((1/Rin + 1/Rsw) / C)   (switch closed)
         tau_off = Rin * C                       (switch open)
     (This is the homogeneous/zero-input decay of dx/dt = A*x integrated
     over each interval with the fixed A1/A2 from the deck's own header
     comment; matches the g0()/--selftest derivation already in
     refs/r2_analytic.py for the forced response, but this is the
     UNFORCED, one-period multiplier -- the quantity `.pss method=shoot`
     reports as spectral_radius.)

Usage:
    python3 tools/pss_ref.py eig '[[0,-1],[1,0]]'
    python3 tools/pss_ref.py r2 --rin 10e3 --rsw 10e3 --c 1e-9 \\
        --fclk 100e3 --duty 0.5
"""
import argparse
import cmath
import json
import random
import sys


# ---------------------------------------------------------------------
# 1. Independent eigenvalue check: Faddeev-LeVerrier + Durand-Kerner.

def _matmul(A, B, n):
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]


def _trace(A, n):
    return sum(A[i][i] for i in range(n))


def _faddeev_leverrier(A, n):
    """Characteristic polynomial coefficients c0..cn (c0=1, monic) of
    det(x*I - A) = x^n + c1*x^(n-1) + ... + cn. Exact for exact input;
    here run in double precision, which is enough to cross-check the C
    eigensolver to several significant figures."""
    M = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    c = [1.0]
    for k in range(1, n + 1):
        AM = _matmul(A, M, n)
        ck = -_trace(AM, n) / k
        c.append(ck)
        M = [[AM[i][j] + (ck if i == j else 0.0) for j in range(n)]
             for i in range(n)]
    return c


def _poly_eval(c, x):
    r = 0
    for coef in c:
        r = r * x + coef
    return r


def _durand_kerner(c, n, iters=4000, tol=1e-14, seed=0):
    rnd = random.Random(seed)
    roots = [cmath.exp(2j * cmath.pi * k / n) * (0.4 + 0.6 * rnd.random()) + 0.1
             for k in range(n)]
    for _ in range(iters):
        newroots = []
        maxdelta = 0.0
        for i in range(n):
            xi = roots[i]
            denom = 1.0
            for j in range(n):
                if i != j:
                    denom *= (xi - roots[j])
            if abs(denom) < 1e-300:
                denom = 1e-300
            delta = _poly_eval(c, xi) / denom
            newroots.append(xi - delta)
            maxdelta = max(maxdelta, abs(delta))
        roots = newroots
        if maxdelta < tol:
            break
    return roots


def eig_reference(A, iters=4000, seed=0):
    """A: list of lists (n x n, real). Returns n complex eigenvalues,
    sorted by (real, imag) for stable comparison against another method's
    output (which may return them in a different order)."""
    n = len(A)
    c = _faddeev_leverrier(A, n)
    roots = _durand_kerner(c, n, iters=iters, seed=seed)
    return sorted(roots, key=lambda z: (round(z.real, 9), round(z.imag, 9)))


# ---------------------------------------------------------------------
# 2. R2 analytic Floquet multiplier.

def r2_floquet_analytic(rin, rsw, c, fclk, duty):
    T = 1.0 / fclk
    t_on = duty * T
    t_off = T - t_on
    a1 = -(1.0 / rin + 1.0 / rsw) / c   # switch closed
    a2 = -1.0 / (rin * c)                # switch open
    tau_on = 1.0 / abs(a1)
    tau_off = 1.0 / abs(a2)
    lam = -(t_on / tau_on + t_off / tau_off)
    import math
    return math.exp(lam)


# ---------------------------------------------------------------------
# 3. THIRD PASS (card #458, plan T4.1-T4.3, task 1/2): honest settled-.tran
# detection and matched-clock-phase orbit comparison, replacing the second
# pass's max/min comparison (tests/test_pss.py's own docstring, second
# pass, already flagged this as "not yet the tighter target"). Both read
# ngspice `wrdata` ascii output directly (t,val pairs per requested
# signal, column 2*i/2*i+1 for the i-th signal on the wrdata line) -- no
# ngspice internals, pure post-processing, independent of pssshoot.c/
# psseigen.c the same way section 1/2 above are.

import bisect


def read_wrdata_pair(path, signal_index):
    """(times, values) for the signal_index-th `wrdata` column pair (0-based,
    in the order the deck's own `wrdata ... sig0 sig1 ...` line lists them)."""
    times, vals = [], []
    c0, c1 = 2 * signal_index, 2 * signal_index + 1
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) > c1:
                try:
                    t = float(p[c0]); v = float(p[c1])
                except ValueError:
                    continue
                times.append(t)
                vals.append(v)
    return times, vals


def interp(times, vals, t):
    """Linear interpolation; times must be ascending (ngspice wrdata output
    always is, including through breakpoints -- it is the accepted-step
    sequence in order)."""
    if t <= times[0]:
        return vals[0]
    if t >= times[-1]:
        return vals[-1]
    i = bisect.bisect_right(times, t) - 1
    if i >= len(times) - 1:
        return vals[-1]
    t0, t1 = times[i], times[i + 1]
    v0, v1 = vals[i], vals[i + 1]
    return v0 if t1 == t0 else v0 + (t - t0) / (t1 - t0) * (v1 - v0)


def phase_samples(times, vals, period, n_phases, t_start):
    """One period's worth of samples at n_phases evenly spaced clock
    phases starting at t_start (t_start + k*period/n_phases, k=0..n_phases-1)."""
    return [interp(times, vals, t_start + k * period / n_phases)
            for k in range(n_phases)]


def settle_period_index(times, vals, period, abstol_floor, rel_tol,
                         t_start=0.0, phase_frac=0.5):
    """Sample one point per period at a fixed phase (t_start + p*period +
    phase_frac*period); return (settle_period_index, settle_time,
    achieved_rel_diff) for the FIRST period p>=1 where
    |s[p]-s[p-1]| / max(|s[p]|,|s[p-1]|,abstol_floor) <= rel_tol -- the
    same scale/residual convention pssshoot.c's own Newton loop uses, so
    "settled" here means the same thing "converged" means there. If the
    trace never reaches rel_tol, returns (None, None, best_rel_diff_seen)
    -- honest: this is a "not reached within the reference run", not a
    fabricated settle time and not silently the last point."""
    tmax = times[-1]
    n_periods = int((tmax - t_start) / period)
    samples = [interp(times, vals, t_start + p * period + phase_frac * period)
               for p in range(n_periods)]
    best = (None, None)
    for p in range(1, len(samples)):
        scale = max(abs(samples[p]), abs(samples[p - 1]), abstol_floor)
        rel = abs(samples[p] - samples[p - 1]) / scale
        if best[1] is None or rel < best[1]:
            best = (p, rel)
        if rel <= rel_tol:
            return (p, t_start + p * period + phase_frac * period, rel)
    return (None, None, best[1])


def ripple(times, vals, period, t_start):
    n = 200
    s = phase_samples(times, vals, period, n, t_start)
    return max(s) - min(s)


# ---------------------------------------------------------------------

def main(argv):
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_eig = sub.add_parser("eig", help="independent eigenvalues of a JSON matrix")
    p_eig.add_argument("matrix", help='row-major JSON, e.g. "[[0,-1],[1,0]]"')

    p_r2 = sub.add_parser("r2", help="R2 analytic Floquet multiplier")
    p_r2.add_argument("--rin", type=float, default=10e3)
    p_r2.add_argument("--rsw", type=float, default=10e3)
    p_r2.add_argument("--c", type=float, default=1e-9)
    p_r2.add_argument("--fclk", type=float, default=100e3)
    p_r2.add_argument("--duty", type=float, default=0.5)

    p_settle = sub.add_parser("settle", help="first period a wrdata trace settles to rel_tol")
    p_settle.add_argument("path")
    p_settle.add_argument("signal_index", type=int)
    p_settle.add_argument("period", type=float)
    p_settle.add_argument("--rel-tol", type=float, default=1e-6)
    p_settle.add_argument("--abstol-floor", type=float, default=1e-3)
    p_settle.add_argument("--t-start", type=float, default=0.0)

    p_phase = sub.add_parser("phasecmp", help="matched-clock-phase comparison, shoot vs tran")
    p_phase.add_argument("shoot_path")
    p_phase.add_argument("tran_path")
    p_phase.add_argument("signal_index_shoot", type=int)
    p_phase.add_argument("signal_index_tran", type=int)
    p_phase.add_argument("period", type=float)
    p_phase.add_argument("--n-phases", type=int, default=20)
    p_phase.add_argument("--tran-t-start", type=float, required=True,
                          help="start of the tran reference's LAST period to compare")

    args = ap.parse_args(argv)

    if args.cmd == "eig":
        A = json.loads(args.matrix)
        for z in eig_reference(A):
            print(f"{z.real:.8f}{z.imag:+.8f}i")
    elif args.cmd == "r2":
        lam = r2_floquet_analytic(args.rin, args.rsw, args.c, args.fclk, args.duty)
        print(f"{lam:.8f}")
    elif args.cmd == "settle":
        times, vals = read_wrdata_pair(args.path, args.signal_index)
        p, t, rel = settle_period_index(times, vals, args.period,
                                         args.abstol_floor, args.rel_tol,
                                         t_start=args.t_start)
        if p is None:
            print(f"NOT-SETTLED best_rel_diff={rel:.6e} "
                  f"(never reached rel_tol={args.rel_tol:.3e} within the "
                  f"available trace, {times[-1]:.6g}s)")
        else:
            print(f"SETTLED period={p} t={t:.9g} rel_diff={rel:.6e}")
    elif args.cmd == "phasecmp":
        st, sv = read_wrdata_pair(args.shoot_path, args.signal_index_shoot)
        tt, tv = read_wrdata_pair(args.tran_path, args.signal_index_tran)
        shoot_t0 = st[0]
        ssamp = phase_samples(st, sv, args.period, args.n_phases, shoot_t0)
        tsamp = phase_samples(tt, tv, args.period, args.n_phases, args.tran_t_start)
        rip = max(tsamp) - min(tsamp)
        diffs = [abs(a - b) for a, b in zip(ssamp, tsamp)]
        maxdiff = max(diffs)
        meandiff = sum(diffs) / len(diffs)
        pct_max = 100.0 * maxdiff / rip if rip > 0 else float("inf")
        pct_mean = 100.0 * meandiff / rip if rip > 0 else float("inf")
        print(f"ripple={rip:.6g} max_diff={maxdiff:.6g} ({pct_max:.4f}% of "
              f"ripple) mean_diff={meandiff:.6g} ({pct_mean:.4f}% of ripple) "
              f"n_phases={args.n_phases}")


if __name__ == "__main__":
    main(sys.argv[1:])
