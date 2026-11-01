#!/usr/bin/env python3
"""
r2_analytic.py -- exact sideband-0 (G0) transfer function of refs/r2_switched_rc.cir
(ngspice_augmentation, card #450, 2026-09-27)

Circuit (see r2_switched_rc.cir header for the full derivation): a single state
x = v(vx), the capacitor voltage, driven through Rin from vin, with a clocked
switch that additionally connects vx to ground through Rsw during the first
D*T of every clock period T (switch closed) and leaves it disconnected during
the remaining (1-D)*T (switch open). This is a scalar linear periodically
time-varying (LPTV) system: two constant-coefficient intervals per period,
exact matrix (here scalar) exponential propagation between them, no averaging.

State equation dx/dt = A_k x + B_k u, u(t) = exp(j*omega*t), in interval k:
  closed (0 <= t < D*T):      A1 = -(1/Rin + 1/Rsw)/C,   B1 = 1/(Rin*C)
  open   (D*T <= t < T):      A2 = -1/(Rin*C),            B2 = 1/(Rin*C)

Exact solution over one interval [t0, t1] with constant (A, B) and start x(t0)=x0:
  x(t) = exp(A*(t-t0))*x0 + B/(j*omega - A) * [exp(j*omega*t) - exp(A*(t-t0))*exp(j*omega*t0)]

Periodic-steady-state closure: for a T-periodic system forced by exp(j*omega*t),
the state obeys x(t+T) = exp(j*omega*T) * x(t). Propagating x0 through the two
intervals (closed then open) and setting x(T) = exp(j*omega*T)*x0 gives one
linear equation for the unknown complex x0 -- solved below in closed form
(solve_x0), not by averaging.

Sideband-0 output (the component of the periodic steady state at the same
frequency as the excitation, i.e. what a lock-in / FRA measurement extracts):
  G0(omega) = (1/T) * integral_0^T  x(t) * exp(-j*omega*t) dt
which is evaluated in closed form per interval (integral_of_segment) using the
same exponential expressions, then summed and divided by T.

Self-check: solve_x0/G0 (closed form) is cross-checked against a brute-force
fixed-step exponential-integrator simulation run to a periodic steady state
(run_periods large) -- see self_test(). The two must agree to the tolerance
printed; this is a correctness check on the algebra above, not a substitute
for the FRA/PAC comparison done elsewhere in the project.
"""
import cmath
import math
import argparse

# --- card #452 finding (2026-09-27, see tests/check_r2_conversion.py's report) --------------
# A from-scratch, independent third method (a frequency-domain harmonic conversion-matrix /
# HTM solution of the same LPTV node equation, tests/check_r2_conversion.py) was built to
# arbitrate a ~7.8deg/3.6% disagreement a previous worker found between this file's g0() and
# tests/test_fra_lptv.py's transient FRA measurement. Result: the HTM solution agrees with
# g0() to ~1e-6 relative or better at every frequency checked (1 kHz-130 kHz), agrees with
# g0()'s own already-recorded near-DC cross-check against the deck's measured DC-driven .tran
# average (r2_switched_rc.cir header) to ~3e-6, and independently satisfies the ODE on
# resynthesis and the H0(-omega)=conj(H0(omega)) reality condition. g0() below is therefore
# NOT the bug. The FRA measurement's disagreement was tracked down to
# tests/test_fra_lptv.py's r2_transfer() calling fra.run_transfer_point with settle_periods=0:
# the measurement window opens at the exact instant the injected probe sine turns on, with
# zero periods allowed for the RESPONSE TO THAT NEW PERTURBATION (not the background switching,
# which had settled) to reach its own periodic steady state; the resulting turn-on transient
# (decaying at ~the circuit's switching/Floquet rate, roughly 0.22 amplitude per fclk period,
# not per excitation period f) leaks into the sideband-0 DFT extraction, worst when the
# measurement window is short in units of fclk periods -- exactly the pattern seen across
# R2_RATIOS (bigger error at 20/40 kHz where measure_periods gives only ~5 switching periods of
# window, smaller at 1/2/5 kHz and 130 kHz where it's naturally 10-100). Verified empirically:
# forcing a few periods of settling AFTER the injection turns on (settle_periods>0, scaled to
# fclk not f: settle_periods ~= 20*f/fclk) collapses the FRA-vs-g0 disagreement at every tested
# frequency to <0.1% / <0.005 deg. tests/test_fra_lptv.py is not this file's to edit (another
# worker's card); the fix belongs in its r2_transfer() (see the card #452 report). No change to
# the g0()/solve_x0/propagate/integral_of_segment algebra below was needed or made.
#
# The old --selftest (self_test(), below) could not have caught a bug in this recipe even if
# there had been one: it only ever cross-checks the closed form against a brute-force run of
# the SAME "eternal one-sided complex exponential" assumption, never against an independently
# -coded method. htm_crosscheck() below closes that gap: a compact, from-scratch (does not call
# propagate/integral_of_segment/solve_x0) frequency-domain harmonic-conversion-matrix solve of
# the same idealised (Ron=0/Roff=infinite) node equation, wired into --selftest.
# -----------------------------------------------------------------------------------------

# Defaults must match refs/r2_switched_rc.cir's .param line.
DEFAULTS = dict(rin=10e3, rsw=10e3, c=1e-9, fclk=100e3, duty=0.5)


def coeffs(rin, rsw, c):
    A1 = -(1.0 / rin + 1.0 / rsw) / c
    B1 = 1.0 / (rin * c)
    A2 = -1.0 / (rin * c)
    B2 = 1.0 / (rin * c)
    return A1, B1, A2, B2


def propagate(x0, A, B, t0, t1, omega):
    """x(t1) given x(t0)=x0 in a constant-(A,B) interval, forcing exp(j*omega*t)."""
    j = 1j
    dt = t1 - t0
    eA = cmath.exp(A * dt)
    if abs(j * omega - A) < 1e-30:
        # degenerate: j*omega == A (not hit for the component values used here)
        forced = B * dt * cmath.exp(j * omega * t0) * eA
    else:
        forced = B / (j * omega - A) * (cmath.exp(j * omega * t1) - eA * cmath.exp(j * omega * t0))
    return eA * x0 + forced


def integral_of_segment(x0, A, B, t0, t1, omega):
    """closed-form integral_{t0}^{t1} x(t) * exp(-j*omega*t) dt for x(t) as in propagate()."""
    j = 1j
    dt = t1 - t0
    eA = cmath.exp(A * dt)

    # term 1: integral of exp(A*(t-t0))*x0 * exp(-j*omega*t) dt, substituting tau = t-t0
    # so exp(-j*omega*t) = exp(-j*omega*t0)*exp(-j*omega*tau); this keeps every exponent
    # bounded by O(A*dt) and O(omega*t0) (a pure phase, |.|=1) instead of blowing up with
    # absolute time t0 the way exp(-A*t0) alone would over many periods.
    # = x0 * exp(-j*omega*t0) * [exp((A-j*omega)*dt) - 1] / (A - j*omega)
    p = A - j * omega
    if abs(p) < 1e-30:
        term1 = x0 * cmath.exp(-j * omega * t0) * dt
    else:
        term1 = x0 * cmath.exp(-j * omega * t0) * (cmath.exp(p * dt) - 1.0) / p

    # term 2: integral of B/(j*omega-A) * [exp(j*omega*t) - exp(A*(t-t0))*exp(j*omega*t0)] * exp(-j*omega*t) dt
    # = B/(j*omega-A) * [ dt - (exp(p*dt) - 1)/p ]   with p = A - j*omega as above (see header derivation)
    if abs(j * omega - A) < 1e-30:
        term2 = 0.0
    else:
        if abs(p) < 1e-30:
            bracket = dt - dt  # p->0 limit of (exp(p*dt)-1)/p is dt, so bracket -> 0
        else:
            bracket = dt - (cmath.exp(p * dt) - 1.0) / p
        term2 = B / (j * omega - A) * bracket

    return term1 + term2


def solve_x0(omega, rin, rsw, c, fclk, duty):
    T = 1.0 / fclk
    t_sw = duty * T
    A1, B1, A2, B2 = coeffs(rin, rsw, c)
    j = 1j

    # x1 = propagate(x0) through interval 1 is affine in x0: x1 = a1*x0 + b1
    a1 = cmath.exp(A1 * t_sw)
    b1 = propagate(0.0, A1, B1, 0.0, t_sw, omega)
    # x(T) = a2*x1 + b2, affine in x1
    a2 = cmath.exp(A2 * (T - t_sw))
    b2 = propagate(0.0, A2, B2, t_sw, T, omega)

    # x(T) = a2*(a1*x0+b1) + b2 = (a2*a1)*x0 + (a2*b1+b2)
    a_total = a2 * a1
    b_total = a2 * b1 + b2

    # periodicity: a_total*x0 + b_total = exp(j*omega*T)*x0
    denom = cmath.exp(j * omega * T) - a_total
    x0 = b_total / denom
    return x0, T, t_sw, (A1, B1, A2, B2)


def g0(freq_hz, rin=DEFAULTS["rin"], rsw=DEFAULTS["rsw"], c=DEFAULTS["c"],
       fclk=DEFAULTS["fclk"], duty=DEFAULTS["duty"]):
    omega = 2 * math.pi * freq_hz
    x0, T, t_sw, (A1, B1, A2, B2) = solve_x0(omega, rin, rsw, c, fclk, duty)
    x1 = propagate(x0, A1, B1, 0.0, t_sw, omega)
    seg1 = integral_of_segment(x0, A1, B1, 0.0, t_sw, omega)
    seg2 = integral_of_segment(x1, A2, B2, t_sw, T, omega)
    return (seg1 + seg2) / T


def self_test(rin=DEFAULTS["rin"], rsw=DEFAULTS["rsw"], c=DEFAULTS["c"],
               fclk=DEFAULTS["fclk"], duty=DEFAULTS["duty"],
               freq_hz=5000.0, run_periods=4000, steps_per_period=400, tol=1e-6):
    """Cross-check the closed form above against brute-force fixed-step propagation
    run to periodic steady state, then numerically integrated. Returns (closed, brute, reldiff)."""
    omega = 2 * math.pi * freq_hz
    T = 1.0 / fclk
    t_sw = duty * T
    A1, B1, A2, B2 = coeffs(rin, rsw, c)
    n1 = steps_per_period // 2
    n2 = steps_per_period - n1
    dt1 = t_sw / n1
    dt2 = (T - t_sw) / n2

    x = 0.0 + 0.0j
    t = 0.0
    for _ in range(run_periods):
        for _ in range(n1):
            x = propagate(x, A1, B1, t, t + dt1, omega)
            t += dt1
        for _ in range(n2):
            x = propagate(x, A2, B2, t, t + dt2, omega)
            t += dt2

    # one more period, numerically integrating x(t)*exp(-j*omega*t)
    acc = 0.0 + 0.0j
    t0 = t
    for _ in range(n1):
        acc += integral_of_segment(x, A1, B1, t, t + dt1, omega)
        x = propagate(x, A1, B1, t, t + dt1, omega)
        t += dt1
    for _ in range(n2):
        acc += integral_of_segment(x, A2, B2, t, t + dt2, omega)
        x = propagate(x, A2, B2, t, t + dt2, omega)
        t += dt2
    brute = acc / T

    closed = g0(freq_hz, rin, rsw, c, fclk, duty)
    reldiff = abs(closed - brute) / max(abs(brute), 1e-30)
    return closed, brute, reldiff


def htm_crosscheck(freq_hz, rin=DEFAULTS["rin"], rsw=DEFAULTS["rsw"], c=DEFAULTS["c"],
                    fclk=DEFAULTS["fclk"], duty=DEFAULTS["duty"], N=40):
    """Independent second method for g0(), card #452: a frequency-domain harmonic
    conversion-matrix (HTM) solve of the SAME idealised (switch closed: conductance to
    ground g_hi=1/rsw; open: g_lo=0 -- no Ron/Roff, matching coeffs()'s A1/A2/B1/B2 above)
    node equation, coded from scratch here (does NOT call propagate/integral_of_segment/
    solve_x0/g0's machinery). This is a compact instance of the fuller derivation and
    convergence study in tests/check_r2_conversion.py (which also uses the deck's actual
    Ron/Roff and cross-checks against the real FRA measurement); this copy is kept small
    and dependency-free so a meaningfully independent check can run inside --selftest.

    Node equation: C dx/dt = -(g_in+g_sw(t))x + g_in*u(t), g_sw(t) a duty-square-wave
    between g_hi/g_lo with Fourier coefficients G_h (closed form, sinc-like, derived from
    G_h = (1/T) int_0^T g_sw(t) exp(-j h w0 t) dt split at duty*T). Sideband equations
    (truncated to |n|,|m|<=N, U=1):
        sum_m [ (G_{n-m} + g_in*delta_nm) + j*(omega+n*w0)*C*delta_nm ] V_m = g_in*delta_n0
    solved by (a small, local) complex Gaussian elimination; returns V_0 = H0(omega).
    """
    g_in = 1.0 / rin
    g_hi = 1.0 / rsw
    g_lo = 0.0
    T = 1.0 / fclk
    w0 = 2 * math.pi * fclk
    DT = duty * T

    def G(h):
        if h == 0:
            return duty * g_hi + (1 - duty) * g_lo
        return (g_hi - g_lo) * (1 - cmath.exp(-1j * h * w0 * DT)) / (1j * h * w0 * T)

    size = 2 * N + 1
    Gc = [G(h) for h in range(-2 * N, 2 * N + 1)]
    omega = 2 * math.pi * freq_hz
    M = [[0j] * (size + 1) for _ in range(size)]
    for i in range(size):
        n = i - N
        row = M[i]
        for jx in range(size):
            m = jx - N
            row[jx] = Gc[(n - m) + 2 * N]
        row[i] += g_in + 1j * (omega + n * w0) * c
        row[size] = g_in if n == 0 else 0j
    for col in range(size):
        piv = max(range(col, size), key=lambda r: abs(M[r][col]))
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
        pr = M[col]
        pv = pr[col]
        for r in range(col + 1, size):
            row = M[r]
            f = row[col] / pv
            if f:
                for k in range(col, size + 1):
                    row[k] -= f * pr[k]
    x = [0j] * size
    for i in range(size - 1, -1, -1):
        row = M[i]
        s = row[size]
        for k in range(i + 1, size):
            s -= row[k] * x[k]
        x[i] = s / row[i]
    return x[N]


def main():
    ap = argparse.ArgumentParser(description="Exact sideband-0 transfer function of r2_switched_rc.cir")
    ap.add_argument("freq_hz", type=float, nargs="?", default=None,
                     help="frequency in Hz to evaluate G0 at (omit to run --selftest only)")
    ap.add_argument("--rin", type=float, default=DEFAULTS["rin"])
    ap.add_argument("--rsw", type=float, default=DEFAULTS["rsw"])
    ap.add_argument("--c", type=float, default=DEFAULTS["c"])
    ap.add_argument("--fclk", type=float, default=DEFAULTS["fclk"])
    ap.add_argument("--duty", type=float, default=DEFAULTS["duty"])
    ap.add_argument("--selftest", action="store_true", help="run the brute-force cross-check and print it")
    args = ap.parse_args()

    if args.selftest or args.freq_hz is None:
        # brute-force cross-check (same eternal-complex-exponential assumption as g0() itself;
        # extended past fclk/2 -- 70/130 kHz -- since card #452 found this file's original
        # 49.9 kHz cutoff never exercised the region a previous worker suspected of hiding a bug)
        for f in (100.0, 1e3, 5e3, 20e3, 49.9e3, 70e3, 130e3):
            closed, brute, reldiff = self_test(args.rin, args.rsw, args.c, args.fclk, args.duty, freq_hz=f)
            print(f"selftest f={f:10.1f} Hz  closed={closed:.6e}  brute={brute:.6e}  reldiff={reldiff:.3e}")
        # independent second method (htm_crosscheck, from-scratch frequency-domain HTM solve,
        # does not call propagate/solve_x0/g0's own code) -- see the card #452 note above main().
        print()
        for f in (1e3, 20e3, 45e3, 130e3):
            closed = g0(f, args.rin, args.rsw, args.c, args.fclk, args.duty)
            htm = htm_crosscheck(f, args.rin, args.rsw, args.c, args.fclk, args.duty)
            reldiff = abs(closed - htm) / max(abs(htm), 1e-30)
            print(f"htm_crosscheck f={f:10.1f} Hz  g0={closed:.6e}  htm={htm:.6e}  reldiff={reldiff:.3e}")

    if args.freq_hz is not None:
        val = g0(args.freq_hz, args.rin, args.rsw, args.c, args.fclk, args.duty)
        mag_db = 20 * math.log10(abs(val)) if abs(val) > 0 else float("-inf")
        phase_deg = math.degrees(cmath.phase(val))
        print(f"f={args.freq_hz:.6g} Hz  G0={val.real:+.6e}{val.imag:+.6e}j  "
              f"|G0|={abs(val):.6e} ({mag_db:.3f} dB)  angle={phase_deg:.3f} deg")


if __name__ == "__main__":
    main()
