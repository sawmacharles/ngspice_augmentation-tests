#!/usr/bin/env python3
"""check_r2_conversion.py -- from-scratch harmonic-conversion-matrix (HTM) solution of
refs/r2_switched_rc.cir's sideband-0 transfer function H0(f) = Vx(f)/Vin(f), card #452.

Independent THIRD method (does not import or reuse any code from refs/r2_analytic.py's
algebra) to arbitrate between refs/r2_analytic.py's g0() (piecewise-LTI eternal-complex-
exponential propagation) and tests/test_fra_lptv.py's transient FRA measurement, which
disagree by up to ~7.8 deg / 3.6% from 10 kHz to 130 kHz.

Method (harmonic conversion matrix / Fourier-domain LPTV solution -- see e.g. the
standard "conversion matrix" treatment of periodically-switched linear circuits):
node equation for x = v(vx), derived from the deck (refs/r2_switched_rc.cir):

    C dx/dt = (u - x)/Rin - x*g_sw(t)
            = -(g_in + g_sw(t)) x + g_in u(t),   g_in = 1/Rin

g_sw(t) is the conductance from vx to ground through Rsw in series with the switch S
(sw1->0), i.e. g_sw(t) = 1/(Rsw + Rswitch(t)) with Rswitch(t) = Ron during the first
duty*T of each clock period T=1/fclk (switch closed, clk high) and Roff during the rest
(switch open, clk low) -- using the deck's actual .model SWCLK SW(...) ron/roff values,
not idealised to 0/infinity (that is refs/r2_analytic.py's idealisation, not made here).

g_sw(t) is a 50%-duty (duty is actually read from the deck, not hardcoded) square wave
between g_hi=1/(Rsw+Ron) and g_lo=1/(Rsw+Roff); its Fourier series g_sw(t) = sum_h G_h
exp(j h w0 t) is derived below in closed form (a sinc-like 1/h decay, odd h only for
50% duty, doesubled to G_0 = duty*g_hi + (1-duty)*g_lo for h=0).

Forced by u(t) = U*exp(j*omega*t) (u(t)=exp(j*omega*t) is the SAME idealisation
refs/r2_analytic.py's g0() uses -- an eternal complex exponential, periodic steady
state), the state is periodic-modulated: x(t) = sum_n V_n exp(j*(omega+n*w0)*t). Matching
Fourier coefficients of the ODE (g_sw(t)*x(t) is a convolution in the harmonic index)
gives, for every sideband n:

    sum_m [ (G_{n-m} + g_in*delta_nm) + j*(omega+n*w0)*C*delta_nm ] V_m = g_in*U*delta_n0

truncated to |n|,|m| <= N and solved by complex Gaussian elimination (own implementation,
standard library only). H0(omega) = V_0/U is what an FRA/PAC measurement at frequency
omega reports (the sideband-0/same-frequency component), directly comparable to
refs/r2_analytic.py's g0() and to tests/test_fra_lptv.py's measured H.

Run:
    python3 tests/check_r2_conversion.py            # full report (convergence, checks, 3-way table)
    python3 tests/check_r2_conversion.py --quick     # skip the N=800 point and the FRA runs (fast, offline)

Exit status: 0 (this file is a report script, not a pass/fail gate; see its printed output).
"""
import argparse
import cmath
import math
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DECK = os.path.join(ROOT, "refs", "r2_switched_rc.cir")

sys.path.insert(0, os.path.join(ROOT, "refs"))
import r2_analytic as ra  # noqa: E402  (read-only reuse of g0() for comparison only; its algebra is NOT reused)

sys.path.insert(0, os.path.join(ROOT, "tools"))
import fra  # noqa: E402  (read-only reuse, for the root-cause demonstration section only)

sys.path.insert(0, HERE)
import test_fra_lptv as tf  # noqa: E402  (read-only reuse of the FRA harness call + frequency list, for an apples-to-apples 3-way compare; this module's own algebra is not touched)


# ---------------------------------------------------------------------------
# deck parsing (own, from-scratch parser -- does not import tools/stb.py's)
# ---------------------------------------------------------------------------

_SUFFIX = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "m": 1e-3,
           "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}


def spice_num(s):
    s = s.strip()
    m = re.match(r"^([-+0-9.eE]+)(meg|[tgkmunpfMEGKUNPF]+)?$", s)
    if not m:
        raise ValueError(f"not a spice number: {s!r}")
    val = float(m.group(1))
    suf = (m.group(2) or "").lower()
    if suf:
        if suf not in _SUFFIX:
            raise ValueError(f"unknown suffix in {s!r}")
        val *= _SUFFIX[suf]
    return val


def parse_deck(path):
    """Pull rin/rsw/c/fclk/duty from .param lines and vt/vh/ron/roff from the SWCLK
    .model line, by direct regex on the deck text -- independent of tools/stb.py and of
    refs/r2_analytic.py's DEFAULTS dict (which is hand-copied from the same deck; this
    reads the deck itself so a future edit to the .param line cannot silently desync the
    two)."""
    text = open(path).read()
    params = {}
    for line in re.finditer(r"^\.param\s+(.+)$", text, re.M | re.I):
        for tok in line.group(1).split():
            if "=" not in tok:
                continue
            k, v = tok.split("=", 1)
            try:
                params[k.strip().lower()] = spice_num(v.strip())
            except ValueError:
                pass  # e.g. per={1/fclk}: an expression, not needed here (T is derived from fclk)
    mm = re.search(r"\.model\s+SWCLK\s+SW\(([^)]*)\)", text, re.I)
    if not mm:
        raise ValueError("SWCLK .model line not found in deck")
    model = {}
    for tok in mm.group(1).split():
        k, v = tok.split("=", 1)
        model[k.strip().lower()] = spice_num(v.strip())
    for need in ("rin", "rsw", "c", "fclk", "duty"):
        if need not in params:
            raise ValueError(f"deck .param missing {need!r}")
    for need in ("ron", "roff"):
        if need not in model:
            raise ValueError(f".model SWCLK missing {need!r}")
    return params, model


# ---------------------------------------------------------------------------
# g_sw(t) Fourier coefficients (50%-or-general-duty square wave between g_hi, g_lo)
# ---------------------------------------------------------------------------

def fourier_coeffs_setup(g_hi, g_lo, duty, T):
    """Returns a function G(h) -> complex Fourier coefficient of g_sw(t), where
    g_sw(t) = g_hi for t in [0, duty*T) and g_lo for t in [duty*T, T), period T.

    Derivation (closed form, from the definition G_h = (1/T) int_0^T g(t) exp(-j h w0 t) dt,
    split at DT=duty*T):
        G_0   = duty*g_hi + (1-duty)*g_lo
        G_h   = (g_hi-g_lo) * (1 - exp(-j h w0 DT)) / (j h w0 T),  h != 0
    For duty=0.5, DT=T/2, h*w0*DT = h*pi, so (1-exp(-j h pi)) is 0 for even h and 2 for
    odd h -- the standard 50%-square-wave "odd harmonics only" result -- but the code
    below does NOT assume duty=0.5; it uses whatever duty parse_deck() read.
    """
    w0 = 2 * math.pi / T
    DT = duty * T
    G0 = duty * g_hi + (1 - duty) * g_lo

    def G(h):
        if h == 0:
            return G0
        return (g_hi - g_lo) * (1 - cmath.exp(-1j * h * w0 * DT)) / (1j * h * w0 * T)

    return G


def fourier_coeffs_selftest(G, g_hi, g_lo, duty, T, hs=(1, 2, 3, 7, 40)):
    """Cross-check the closed-form G(h) above against brute-force trapezoidal quadrature
    of the exact (instantaneous-edge) time-domain square wave, independent of the
    closed-form algebra. Returns list of (h, closed, numeric, reldiff); reldiff is scaled
    by g_hi (not by |numeric|), because even h's true value is ~0 (a null, not just
    small) for this 50%-duty square wave, so a relative-to-itself error is meaningless
    quadrature noise there -- scaling by the signal's own natural size (g_hi) instead
    correctly shows those near-zero values agree with 0 to ~1e-16, i.e. to double
    precision, not "84% off"."""
    n_pts = 200000
    dt = T / n_pts
    DT = duty * T
    out = []
    for h in hs:
        w0 = 2 * math.pi / T
        acc = 0j
        for i in range(n_pts):
            t = (i + 0.5) * dt  # midpoint rule, avoids sampling exactly on the edge
            g = g_hi if (t % T) < DT else g_lo
            acc += g * cmath.exp(-1j * h * w0 * t) * dt
        numeric = acc / T
        closed = G(h)
        reldiff = abs(closed - numeric) / g_hi
        out.append((h, closed, numeric, reldiff))
    return out


# ---------------------------------------------------------------------------
# complex dense linear solve (own Gaussian elimination, partial pivoting; stdlib only)
# ---------------------------------------------------------------------------

def complex_solve(A, b):
    n = len(b)
    M = [row[:] + [bv] for row, bv in zip(A, b)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) == 0.0:
            raise ValueError("singular matrix in complex_solve")
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
        pivrow = M[col]
        pivval = pivrow[col]
        for r in range(col + 1, n):
            row = M[r]
            factor = row[col] / pivval
            if factor != 0:
                for k in range(col, n + 1):
                    row[k] -= factor * pivrow[k]
    x = [0j] * n
    for i in range(n - 1, -1, -1):
        row = M[i]
        s = row[n]
        for k in range(i + 1, n):
            s -= row[k] * x[k]
        x[i] = s / row[i]
    return x


# ---------------------------------------------------------------------------
# HTM assembly and solve
# ---------------------------------------------------------------------------

def htm_solve(omega, N, G, g_in, c, w0, U=1.0 + 0j):
    """Solve the truncated (|n|,|m|<=N) sideband system, return the full V dict {n: V_n}."""
    size = 2 * N + 1
    Gcache = [G(h) for h in range(-2 * N, 2 * N + 1)]  # Gcache[h+2N] = G(h)

    def Gh(h):
        return Gcache[h + 2 * N]

    A = [[0j] * size for _ in range(size)]
    b = [0j] * size
    for i in range(size):
        n = i - N
        row = A[i]
        for jx in range(size):
            m = jx - N
            row[jx] = Gh(n - m)
        row[i] += g_in + 1j * (omega + n * w0) * c
        b[i] = g_in * U if n == 0 else 0j
    V = complex_solve(A, b)
    return {n: V[n + N] for n in range(-N, N + 1)}


def htm_H0(freq_hz, N, G, g_in, c, fclk):
    omega = 2 * math.pi * freq_hz
    w0 = 2 * math.pi * fclk
    Vd = htm_solve(omega, N, G, g_in, c, w0)
    return Vd[0]


# ---------------------------------------------------------------------------
# reporting helpers
# ---------------------------------------------------------------------------

def mag_deg(z):
    return abs(z), math.degrees(cmath.phase(z))


def pct_deg_diff(got, ref):
    gm, gph = mag_deg(got)
    rm, rph = mag_deg(ref)
    pct = (gm - rm) / rm * 100 if rm else float("inf")
    dph = gph - rph
    while dph > 180:
        dph -= 360
    while dph < -180:
        dph += 360
    return pct, dph


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                     help="skip the N=800 convergence point and the FRA (ngspice) runs")
    args = ap.parse_args()

    params, model = parse_deck(DECK)
    rin, rsw, c, fclk, duty = (params[k] for k in ("rin", "rsw", "c", "fclk", "duty"))
    ron, roff = model["ron"], model["roff"]
    g_in = 1.0 / rin
    g_hi = 1.0 / (rsw + ron)
    g_lo = 1.0 / (rsw + roff)
    T = 1.0 / fclk
    w0 = 2 * math.pi * fclk

    print("=" * 100)
    print("check_r2_conversion.py -- harmonic conversion-matrix (HTM) check of R2's sideband-0 H0(f)")
    print("=" * 100)
    print(f"deck params (parsed from {DECK}):")
    print(f"  rin={rin:.6g}  rsw={rsw:.6g}  c={c:.6g}  fclk={fclk:.6g}  duty={duty:.6g}")
    print(f"  SWCLK: ron={ron:.6g}  roff={roff:.6g}")
    print(f"  g_in={g_in:.9e} S   g_hi=1/(rsw+ron)={g_hi:.9e} S   g_lo=1/(rsw+roff)={g_lo:.9e} S")
    print(f"  (idealised g_hi would be 1/rsw={1/rsw:.9e} S, g_lo would be 0 -- "
          f"actual deck values differ by {(g_hi-1/rsw)/(1/rsw)*100:+.4g}% / are {g_lo:.3e} S, not 0)")

    G = fourier_coeffs_setup(g_hi, g_lo, duty, T)

    # ---- Fourier-coefficient self-check (closed form vs brute-force quadrature) ----
    print("\n--- G_h closed-form vs numeric quadrature (independent check of the Fourier-series algebra) ---")
    for h, closed, numeric, reldiff in fourier_coeffs_selftest(G, g_hi, g_lo, duty, T):
        print(f"  h={h:4d}  closed={closed:.6e}  numeric={numeric:.6e}  reldiff={reldiff:.3e}")

    # ---- convergence study ----
    print("\n--- HTM convergence vs N (worst-mismatch region, f=45 kHz = 0.45*fclk) ---")
    f_conv = 45e3
    omega_conv = 2 * math.pi * f_conv
    Ns = [10, 20, 50, 100, 200] if args.quick else [10, 20, 50, 100, 200, 400, 800]
    prev = None
    conv_rows = []
    for N in Ns:
        t0 = time.time()
        Vd = htm_solve(omega_conv, N, G, g_in, c, w0)
        dt_run = time.time() - t0
        v0 = Vd[0]
        m, ph = mag_deg(v0)
        if prev is not None:
            pct, dph = pct_deg_diff(v0, prev)
            print(f"  N={N:5d}  V0={v0.real:+.9e}{v0.imag:+.9e}j  |V0|={m:.9e}  ang={ph:+.6f} deg  "
                  f"dvs_prevN={pct:+.3e}% /{dph:+.3e}deg  ({dt_run:.1f}s)")
        else:
            print(f"  N={N:5d}  V0={v0.real:+.9e}{v0.imag:+.9e}j  |V0|={m:.9e}  ang={ph:+.6f} deg  ({dt_run:.1f}s)")
        conv_rows.append((N, v0, dt_run))
        prev = v0
    N_final = Ns[-1]
    print(f"  -> using N={N_final} as the converged HTM answer for the checks below "
          f"(last-step change {abs(pct_deg_diff(conv_rows[-1][1], conv_rows[-2][1])[0]):.2e}% / "
          f"{abs(pct_deg_diff(conv_rows[-1][1], conv_rows[-2][1])[1]):.2e} deg)")

    # ---- DC limit check ----
    # NOTE: the naive "averaged circuit" DC gain g_in/(g_in+duty*g_hi+(1-duty)*g_lo) is NOT the
    # right ground truth here -- it assumes x(t) is constant over a period, which is false (x(t)
    # has real ripple synchronized with g(t), so <g(t)x(t)> != <g(t)><x(t)>; this circuit is built
    # specifically so the switched and averaged models differ, see r2_switched_rc.cir's header and
    # test_fra_lptv.py's own check 4(a) red-proof). The correct DC reference is the deck's own
    # measured DC-driven .tran average (r2_switched_rc.cir header/comment, vx_avg_final, 2026-09-27
    # run) and r2_analytic.g0()'s own already-recorded near-DC cross-check against it.
    print("\n--- DC limit check: HTM at f->0 vs the deck's measured DC-driven .tran average ---")
    f_dc = 1e-6
    Vd_dc = htm_solve(2 * math.pi * f_dc, 200, G, g_in, c, w0)
    htm_dc = Vd_dc[0]
    measured_dc = 6.699138e-01  # r2_switched_rc.cir header: vx_avg_final, 590-600us window
    g0_dc = ra.g0(f_dc).real
    G0_avg = duty * g_hi + (1 - duty) * g_lo
    naive_averaged_dc = g_in / (g_in + G0_avg)
    print(f"  HTM H0(f={f_dc:.3g} Hz)          = {htm_dc.real:.9e}{htm_dc.imag:+.3e}j")
    print(f"  r2_analytic.g0(f={f_dc:.3g} Hz)  = {g0_dc:.9e}  reldiff HTM vs g0 = "
          f"{abs(htm_dc.real - g0_dc)/g0_dc:.3e}")
    print(f"  deck .tran measured DC average  = {measured_dc:.9e}  reldiff HTM vs measured = "
          f"{abs(htm_dc.real - measured_dc)/measured_dc:.3e}")
    print(f"  (naive averaged-conductance model g_in/(g_in+duty*g_hi+(1-duty)*g_lo) = "
          f"{naive_averaged_dc:.9e} -- WRONG reference, off by "
          f"{abs(htm_dc.real - naive_averaged_dc)/naive_averaged_dc*100:.3f}%, "
          f"confirms the switched/averaged mismatch is real even at DC, not just a high-f effect)")

    # ---- ODE resynthesis check ----
    print("\n--- ODE resynthesis check: does the truncated HTM solution satisfy the true (non-Fourier-truncated) ODE? ---")
    f_ode = 45e3
    omega_ode = 2 * math.pi * f_ode
    N_ode = 200
    Vd_ode = htm_solve(omega_ode, N_ode, G, g_in, c, w0)

    def v_of_t(t):
        return sum(Vd_ode[n] * cmath.exp(1j * (omega_ode + n * w0) * t) for n in range(-N_ode, N_ode + 1))

    def g_exact(t):
        tau = t % T
        return g_hi if tau < duty * T else g_lo

    def rhs(t):
        u = cmath.exp(1j * omega_ode * t)
        return (-(g_in + g_exact(t)) * v_of_t(t) + g_in * u) / c

    eps = T * 1e-6
    test_times = [0.3 * T, 0.7 * T, 1.4 * T, duty * T - 5 * eps, duty * T + 5 * eps, 3.3 * T]
    print(f"  (using N={N_ode} truncation; dv/dt by central difference, eps={eps:.3e} s)")
    for t in test_times:
        dv_num = (v_of_t(t + eps) - v_of_t(t - eps)) / (2 * eps)
        r = rhs(t)
        err = abs(dv_num - r)
        scale = max(abs(dv_num), abs(r), 1e-30)
        print(f"  t={t:.6e}s (t/T={t/T:.4f})  dv/dt_numeric={dv_num:.6e}  RHS(-(g_in+g(t))v+g_in*u)={r:.6e}  "
              f"relerr={err/scale:.3e}")

    # ---- reality condition self-check: G0(-omega) == conj(G0(omega))? (both HTM and g0()) ----
    print("\n--- reality condition H0(-omega) = conj(H0(omega)) (physical necessity for a real system) ---")
    f_re = 30e3
    htm_p = htm_H0(f_re, 200, G, g_in, c, fclk)
    htm_m = htm_H0(-f_re, 200, G, g_in, c, fclk)
    print(f"  HTM:  H0(+f)={htm_p:.6e}  conj(H0(-f))={htm_m.conjugate():.6e}  "
          f"reldiff={abs(htm_p-htm_m.conjugate())/abs(htm_p):.3e}")
    g0_p = ra.g0(f_re)
    g0_m = ra.g0(-f_re)
    print(f"  refs/r2_analytic.g0: g0(+f)={g0_p:.6e}  conj(g0(-f))={g0_m.conjugate():.6e}  "
          f"reldiff={abs(g0_p-g0_m.conjugate())/abs(g0_p):.3e}")

    # ---- 3-way comparison table ----
    print("\n" + "=" * 100)
    print("3-WAY COMPARISON: HTM (this file) vs refs/r2_analytic.g0() vs tests/test_fra_lptv.py FRA measurement")
    print("=" * 100)
    # N=200 is already converged to ~7e-8% / 4e-10 deg relative to N=100 (see the convergence
    # study above, which does go to N=800 to demonstrate that), so N=200 is used here to keep
    # the ngspice sweep below (10 points x 2 tran runs) fast; it is not a precision compromise.
    N_cmp = 100 if args.quick else 200
    header = (f"{'f(Hz)':>10s} {'|HTM|':>12s} {'ang_HTM':>9s} {'|g0|':>12s} {'ang_g0':>9s} "
              f"{'d(g0-HTM)%':>10s} {'d_deg':>7s}")
    if not args.quick:
        header += f" {'|FRA|':>12s} {'ang_FRA':>9s} {'d(FRA-HTM)%':>11s} {'d_deg':>7s}"
    print(header)
    rows = []
    for p, q in tf.R2_RATIOS:
        f, mp, ppp = tf.r2_point_settings(p, q)
        htm_v = htm_H0(f, N_cmp, G, g_in, c, fclk)
        g0_v = ra.g0(f)
        hm, hph = mag_deg(htm_v)
        gm, gph = mag_deg(g0_v)
        pct_g0, dph_g0 = pct_deg_diff(g0_v, htm_v)
        line = f"{f:10.6g} {hm:12.6e} {hph:9.3f} {gm:12.6e} {gph:9.3f} {pct_g0:10.3f} {dph_g0:7.3f}"
        fra_v = None
        if not args.quick:
            r = tf.r2_transfer(f, mp, ppp, outdir_name=f"htmcheck_f_{f:.6g}")
            fra_v = r["H"]
            fm, fph = mag_deg(fra_v)
            pct_fra, dph_fra = pct_deg_diff(fra_v, htm_v)
            line += f" {fm:12.6e} {fph:9.3f} {pct_fra:11.3f} {dph_fra:7.3f}"
        print(line)
        rows.append((f, htm_v, g0_v, fra_v))

    print("\n(d(X-HTM)% / d_deg = percent magnitude / degree phase difference of X relative to HTM, this file's")
    print(" independent third method; HTM does not import refs/r2_analytic.py's or tools/fra.py's algebra.)")

    if not args.quick and any(v is not None for _, _, _, v in rows):
        print("\n" + "=" * 100)
        print("ROOT CAUSE: FRA-vs-(HTM,g0) mismatch traced to tests/test_fra_lptv.py's r2_transfer()")
        print("calling fra.run_transfer_point(..., settle_periods=0): the measurement window opens at")
        print("the exact instant the injected probe sine turns on, so it captures that NEW")
        print("perturbation's own turn-on transient (decaying at ~the switching/Floquet rate, i.e. per")
        print("fclk period, not per excitation period f) -- not just the background switching content,")
        print("which had already settled via settle_time. Demonstration: same points, but with")
        print("settle_periods scaled to a fixed number of fclk periods (20) instead of 0:")
        print(f"{'f(Hz)':>10s} {'settle_periods=0':>18s} {'d_deg(0)':>9s}   {'settle_periods~20fclk':>22s} {'d_deg(fix)':>10s}")
        for p, q in tf.R2_RATIOS:
            f, mp, ppp = tf.r2_point_settings(p, q)
            htm_v = htm_H0(f, N_cmp, G, g_in, c, fclk)
            r0 = tf.r2_transfer(f, mp, ppp, outdir_name=f"rootcause0_f_{f:.6g}")
            sp = 20.0 * f / fclk
            r_fix = fra.run_transfer_point(tf.stb.prepare_deck(tf.R2_DECK), "vin", "v(vx)", f, tf.AMPLITUDE,
                                            settle_time=tf.SETTLE_TIME, settle_periods=sp,
                                            measure_periods=mp, points_per_period=ppp,
                                            outdir=os.path.join(tf.OUT, "r2", f"rootcausefix_f_{f:.6g}"))
            _, dph0 = pct_deg_diff(r0["H"], htm_v)
            _, dphf = pct_deg_diff(r_fix["H"], htm_v)
            print(f"{f:10.6g} {'(as in test file)':>18s} {dph0:9.3f}   {'settle_periods='+format(sp,'.2f'):>22s} {dphf:10.4f}")
        print("\nSuggested fix (in tests/test_fra_lptv.py's r2_transfer(), not this file's to edit):")
        print('  pass settle_periods=20.0*f/FCLK (or similar, scaled to fclk not f) to')
        print('  fra.run_transfer_point instead of the current settle_periods=0.')


if __name__ == "__main__":
    main()
