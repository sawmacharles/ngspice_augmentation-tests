#!/usr/bin/env python3
"""tools/pac_exact.py (cards #713, #728, 2026-10-19): an ngspice-independent reference for td .pac control-to-output on the boost /
buck-boost capacitor-node decks of tests/test_pac_charge.py (BOOSTC_580, BB_583; deck texts in tools/pac_exact_decks.py).  Needs
numpy + scipy (not in the system python3 here: a venv, e.g. PAC_EXACT_PYTHON for tests/test_pac_exact.py).

The circuit is re-implemented as an ODE with the ngspice level-1 diode equations (dioload.c: Id with gmin, the -3*vte reverse
branch, depletion charge with FC=0.5, diffusion charge TT*Id, RS on the anode side), the ideal switch as a resistor ron/roff
whose state changes when v(ctrl)-v(ramp) crosses vt+vh (on) / vt-vh (off) (swload.c), and the PULSE ramp as vsrcload.c evaluates
it.  The switch control does not depend on the circuit state, so the switch instants are found EXACTLY (brentq to machine
precision on a known function of t) and the ODE is integrated segment by segment between them (scipy Radau, rtol down to
1e-12): no event-location error at all, which is what limits the perturbed-PSS tool.

Small-signal value at f = fs*p/q: periodic steady state over the LCM period Tl = q/fs by shooting (Newton, finite-difference
monodromy), for Vctrl = 0.8 +/- A*sin(2 pi f t); the Fourier coefficient of v(vout) at f is a quadrature state of the ODE; the
symmetric difference (X(+A) - X(-A)) / (-j A) cancels every even power of A, its error is O(A^2), and Richardson in A^2 on
A, A/2 removes that term too.

usage: pac_exact.py <deck> <f> [A1,A2,...] [rtol]      prints one line per amplitude and the Richardson values, JSON to stdout last
       pac_exact.py --dump <deck> <f> <A> <times-file> <out> [phase-rad]
           time-domain small-signal trajectory (x(+A) - x(-A)) / (2A) at the times of column 1 of <times-file> (first period), for the
           input A sin(w t + phase): columns t, iL, v(sw), v(Cout), [junction voltage], v(vout).  phase 0 gives the response to Im of
           ngspice's e^{jwt} small-signal input, pi/2 the response to Re (compare with the td .pac trajectory dumped by PAC_SS_DUMP).
Reference matrix of card #713 (20/50/100 kHz x the 6 decks): rtol 1e-12, amplitudes 4,2,1,0.5 mV, Richardson in A^2 on the two
smallest; its error bar (spread over amplitudes, rtol 1e-10 / 1e-12) is <= 6e-4 deg, 1.8e-5 dB; ~10 min per point.
"""
import sys, os, math, cmath, json, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pac_exact_decks import DECKS

FS = 500e3
T = 1.0 / FS
VT_TH, VH = 0.0, 0.005
RON, ROFF = 10e-3, 1e6
TR, TF, PW = 0.97 * T, 0.03 * T, 1e-9
IS, RS, N, VJ, M, FC = 1e-9, 0.05, 1.05, 0.7, 0.5, 0.5
GMIN = 1e-12
KOQ = 1.38064852e-23 / 1.6021766208e-19
VTE = N * KOQ * (27.0 + 273.15)
F1 = VJ * (1 - (1 - FC) ** (1 - M)) / (1 - M)
F2 = (1 - FC) ** (1 + M)
F3 = 1 - FC * (1 + M)
L, RDCR, COUT, RESR, RLOAD = 10e-6, 10e-3, 22e-6, 10e-3, 10.0
VCTRL = 0.8


def ramp(tau):
    """PULSE(0 2 0 TR TF PW T) at tau in [0, T) (vsrcload.c; TR+PW+TF > PER, so the fall is cut 1 ns short)."""
    if tau <= 0 or tau >= TR + PW + TF:
        return 0.0
    if TR <= tau <= TR + PW:
        return 2.0
    if tau < TR:
        return 2.0 * tau / TR
    return 2.0 - 2.0 * (tau - (TR + PW)) / TF


def idio(vd):
    if vd >= -3 * VTE:
        e = math.exp(vd / VTE)
        return IS * (e - 1) + GMIN * vd, IS * e / VTE + GMIN
    a = 3 * VTE / (vd * math.e)
    a = a * a * a
    return -IS * (1 + a) + GMIN * vd, IS * 3 * a / vd + GMIN


def qcap(vd, cjo, tt, cd, gd):
    """junction charge's capacitance dq/dvd (depletion + diffusion)."""
    if cjo == 0 and tt == 0:
        return 0.0
    if vd < FC * VJ:
        cdep = cjo * (1 - vd / VJ) ** (-M)
    else:
        cdep = cjo / F2 * (F3 + M * vd / VJ)
    return cdep + tt * gd


class Circuit:
    def __init__(self, name):
        self.topo, self.cjo, self.tt, self.csw, _, self.pinned = DECKS[name]
        self.charged = not (self.cjo == 0 and self.tt == 0)
        self.vin = 5.0 if self.topo == "boost" else 12.0
        self.K = 1 + RESR / RS + RESR / RLOAD
        self.vdguess = 0.5

    def solve_vd(self, vsw, vc):
        """chargeless diode: the junction voltage from KCL through RS (monotone in vd)."""
        s = 1.0 if self.topo == "boost" else -1.0
        # boost: i_RS = (vsw - vd - vout)/RS, vout = (vc + RESR*(vsw - vd)/RS)/K
        # bb:    i_RS = (vout - vsw - vd)/RS, vout = (vc + RESR*(vsw + vd)/RS)/K
        def f(vd):
            if s > 0:
                vout = (vc + RESR * (vsw - vd) / RS) / self.K
                irs = (vsw - vd - vout) / RS
                dirs = (-1 + RESR / RS / self.K) / RS
            else:
                vout = (vc + RESR * (vsw + vd) / RS) / self.K
                irs = (vout - vsw - vd) / RS
                dirs = (RESR / RS / self.K - 1) / RS
            i, g = idio(vd)
            return irs - i, dirs - g
        lo, hi = -1e3 - abs(vsw) - abs(vc), 3.0       # f(lo) > 0 > f(hi): safeguarded Newton
        vd = min(max(self.vdguess, lo), hi)
        for _ in range(200):
            r, d = f(vd)
            if r > 0: lo = vd
            else: hi = vd
            nv = vd - r / d
            if not (lo < nv < hi):
                nv = 0.5 * (lo + hi)
            if abs(nv - vd) < 1e-15 * max(1.0, abs(vd)) or hi - lo < 1e-15:
                vd = nv
                break
            vd = nv
        self.vdguess = vd
        return vd

    def outputs(self, x):
        """vout and the diode current through RS from the state."""
        if self.charged:
            iL, vsw, vc, vd = x[0], x[1], x[2], x[3]
        else:
            iL, vsw, vc = x[0], x[1], x[2]
            vd = self.solve_vd(vsw, vc)
        if self.topo == "boost":
            vout = (vc + RESR * (vsw - vd) / RS) / self.K
            irs = (vsw - vd - vout) / RS
        else:
            vout = (vc + RESR * (vsw + vd) / RS) / self.K
            irs = (vout - vsw - vd) / RS
        return vout, irs, vd

    def rhs(self, t, x, gs, w):
        vout, irs, vd = self.outputs(x)
        iL, vsw, vc = x[0], x[1], x[2]
        if self.topo == "boost":
            diL = (self.vin - iL * RDCR - vsw) / L
            dvsw = (iL - gs * vsw - irs) / self.csw
            ic = irs - vout / RLOAD
        else:
            diL = (vsw - iL * RDCR) / L
            dvsw = (gs * (self.vin - vsw) - iL + irs) / self.csw
            ic = -irs - vout / RLOAD
        dvc = ic / COUT
        out = [diL, dvsw, dvc]
        if self.charged:
            i, g = idio(vd)
            out.append((irs - i) / qcap(vd, self.cjo, self.tt, i, g))
        out += [vout * math.cos(w * t), -vout * math.sin(w * t)]
        return out


def switch_events(A, f, q, phase=0.0):
    """[(t_start, t_end, on?)] over [0, q*T] for Vctrl = 0.8 + A sin(2 pi f t + phase): exact crossings of u = vctrl - ramp at +-vh."""
    w = 2 * math.pi * f
    segs, bps = [], []
    for k in range(q):
        t0 = k * T
        u = lambda t, k=k: VCTRL + A * math.sin(w * t + phase) - ramp(t - k * T)
        # turn-off on the rising ramp: u falls through vt - vh
        toff = brentq(lambda t: u(t) - (VT_TH - VH), t0 + 1e-15, t0 + TR, xtol=1e-24, rtol=1e-15, maxiter=500)
        # turn-on on the falling edge: u rises through vt + vh
        ton = brentq(lambda t: u(t) - (VT_TH + VH), t0 + TR + PW + 1e-15, t0 + T - 1e-15, xtol=1e-24, rtol=1e-15, maxiter=500)
        bps += [toff, ton]
    # state at t=0 is ON (u = 0.8 > vh)
    times = [0.0] + bps + [q * T]
    on = True
    for a, b in zip(times[:-1], times[1:]):
        segs.append((a, b, on))
        on = not on
    return segs


def integrate(ckt, x0, segs, w, rtol):
    nphys = len(x0)
    x = list(x0) + [0.0, 0.0]
    atol = [1e-13, 1e-13, 1e-13] + ([1e-14] if ckt.charged else []) + [1e-21, 1e-21]
    for a, b, on in segs:
        gs = 1.0 / RON if on else 1.0 / ROFF
        sol = solve_ivp(ckt.rhs, (a, b), x, method="Radau", rtol=rtol, atol=atol, args=(gs, w), first_step=min(1e-12, (b - a) / 10))
        if not sol.success:
            raise RuntimeError(sol.message)
        x = list(sol.y[:, -1])
    return np.array(x[:nphys]), complex(x[nphys], x[nphys + 1])


def shoot(ckt, x0, segs, w, rtol, J=None, tol=1e-12, it=30):
    x0 = np.array(x0, float)
    scale = np.array([1.0, 10.0, 10.0] + ([1.0] if ckt.charged else []))
    for i in range(it):
        x1, Q = integrate(ckt, x0, segs, w, rtol)
        F = x1 - x0
        if np.max(np.abs(F) / scale) < tol:
            return x0, Q, J, i
        if J is None:
            n = len(x0)
            J = np.zeros((n, n))
            for c in range(n):
                h = 1e-6 * scale[c]
                xp = x0.copy(); xp[c] += h
                J[:, c] = (integrate(ckt, xp, segs, w, rtol)[0] - x1) / h
            J = J - np.eye(n)
        x0 = x0 - np.linalg.solve(J, F)
    raise RuntimeError("shooting did not converge: %s" % F)


def dbdeg(z):
    return 20 * math.log10(abs(z)), math.degrees(cmath.phase(z))


def dump(name, f, A, tfile, out, phase=0.0, rtol=1e-11):
    """card #728: the time-domain small-signal trajectory, see the usage."""
    q = int(round(FS / f))
    w = 2 * math.pi * f
    ts = np.array([float(l.split()[0]) for l in open(tfile) if l.strip()])
    ckt = Circuit(name)
    x0 = [2.0, 0.0, 8.0 if ckt.topo == "boost" else -8.0] + ([0.0] if ckt.charged else [])
    xb, _, _, _ = shoot(ckt, x0, switch_events(0.0, FS, 1), 2 * math.pi * FS, rtol, tol=1e-12, it=60)
    res, J = {}, None
    for sgn in (+1, -1):
        segs = switch_events(sgn * A, f, q, phase)
        xp, Q, J, ni = shoot(ckt, xb, segs, w, rtol, J=J)
        nphys = len(xp)
        x = list(xp) + [0.0, 0.0]
        atol = [1e-13, 1e-13, 1e-13] + ([1e-14] if ckt.charged else []) + [1e-21, 1e-21]
        vals = np.zeros((len(ts), nphys + 1))
        for a, b, on in segs:
            if a >= T:
                break
            gs = 1.0 / RON if on else 1.0 / ROFF
            sol = solve_ivp(ckt.rhs, (a, b), x, method="Radau", rtol=rtol, atol=atol, args=(gs, w), first_step=min(1e-12, (b - a) / 10),
                            dense_output=True)
            for i in np.where((ts > a) & (ts <= b))[0]:
                xi = sol.sol(ts[i])
                vals[i, :nphys] = xi[:nphys]
                vals[i, nphys] = ckt.outputs(xi[:nphys])[0]
            x = list(sol.y[:, -1])
        res[sgn] = vals
    d = (res[1] - res[-1]) / (2 * A)
    with open(out, "w") as fo:
        for i, t in enumerate(ts):
            fo.write("%.17g %s\n" % (t, " ".join("%.10g" % v for v in d[i])))


def main():
    if sys.argv[1] == "--dump":
        a = sys.argv[2:]
        dump(a[0], float(a[1]), float(a[2]), a[3], a[4], float(a[5]) if len(a) > 5 else 0.0)
        return
    name, f = sys.argv[1], float(sys.argv[2])
    amps = [float(a) for a in sys.argv[3].split(",")] if len(sys.argv) > 3 else [4e-3, 2e-3, 1e-3, 5e-4]
    rtol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-11
    q = int(round(FS / f))
    assert abs(q * f - FS) < 1e-6, "f must divide fs"
    w = 2 * math.pi * f
    ckt = Circuit(name)
    t0 = time.time()
    # baseline one-period orbit (A = 0) for the initial guess and the DC check
    x0 = [2.0, 0.0, 8.0 if ckt.topo == "boost" else -8.0] + ([0.0] if ckt.charged else [])
    xb, Qb, Jb, nb = shoot(ckt, x0, switch_events(0.0, FS, 1), 2 * math.pi * FS, rtol, tol=1e-12, it=60)
    # mean vout over one period
    res = {"deck": name, "f": f, "rtol": rtol, "x0_baseline": list(xb), "rows": []}
    J = None
    xs = xb
    H = {}
    for A in amps:
        Xs = []
        for sgn in (+1, -1):
            segs = switch_events(sgn * A, f, q)
            xp, Q, J, ni = shoot(ckt, xs, segs, w, rtol, J=J)
            Xs.append(Q / (q * T))
        h = (Xs[0] - Xs[1]) / (-1j * A)
        H[A] = h
        d, p = dbdeg(h)
        hp = Xs[0] / (-0.5j * A)
        res["rows"].append({"A": A, "re": h.real, "im": h.imag, "db": d, "deg": p, "one_sided_db": dbdeg(hp)[0], "one_sided_deg": dbdeg(hp)[1]})
        print("%-15s f=%g A=%-8g rtol=%g: %.6f dB %.5f deg   (one-sided +A: %.6f dB %.5f deg)  [%.0fs]" % (name, f, A, rtol, d, p, *dbdeg(hp), time.time() - t0), flush=True)
    srt = sorted(H)
    rich = []
    for a1, a2 in zip(srt[:-1], srt[1:]):
        r = a2 / a1   # a1 < a2
        hr = (H[a1] * r * r - H[a2]) / (r * r - 1)
        rich.append({"A_pair": [a1, a2], "re": hr.real, "im": hr.imag, "db": dbdeg(hr)[0], "deg": dbdeg(hr)[1]})
        print("   Richardson(A^2) %g,%g: %.6f dB %.5f deg" % (a1, a2, *dbdeg(hr)), flush=True)
    res["richardson"] = rich
    print("JSON " + json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
