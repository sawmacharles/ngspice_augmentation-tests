#!/usr/bin/env python3
"""Transient-noise Monte Carlo cross-check for card #509's `.pnoise ...
method="shoot" solver="td"` engine (pacshoot.c's pac_engine_run_adjoint_
noise()) -- an INDEPENDENT reference, requested for the card's second pass
because the LTI check (test_pnoise.py's check_pnoise_lti(), against plain
`.noise`) cannot exercise genuine cyclostationary MODULATION or FOLDING (an
LTI circuit's transfer function does not vary with the PSS phase at all),
and the switched-RC fixtures' own `dev=1` per-device-table checks are only
internal self-consistency, not an external reference.

Method: replace each intrinsic device noise source (a plain resistor here)
with an EXPLICIT `TRNOISE` current source (white, NA/TS chosen so the
injected one-sided PSD matches 4kT/R -- see `trnoise_na()` below for the
NA<->PSD relationship, derived from ngspice's own trnoise generator: iid
N(0,NA^2) samples on a TS-spaced grid, linearly interpolated between grid
points), run N independent seeds of a plain `.tran` (long enough to drop
a start-up transient and get many Welch segments), Hann-window and average
the periodogram (implemented from scratch below -- pure Python, this
machine has no numpy), and compare against `.pnoise`'s own onoise_spectrum
at several frequencies, stating the expected Monte-Carlo spread from the
achieved segment count BEFORE comparing (a chi-squared(2 dof)-per-segment
argument, not a number chosen after the fact).

    python3 tools/pnoise_mc.py [--fixture r2|both]
        NGSPICE=<prefix>/bin/ngspice must be set (or edit NGSPICE below).
"""
import cmath
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXTURES = os.path.join(ROOT, "tests", "fixtures")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/pnoise/mc")

NGSPICE = os.environ.get("NGSPICE", "ngspice")

K_BOLTZ = 1.380649e-23
T_KELVIN = 300.15   # ngspice default CKTtemp (27 C = TNOM)


class CouldNotRun(Exception):
    pass


def run_ngspice(deck_path, timeout=300, env=None):
    t0 = time.time()
    try:
        # Resource rules (2026-09-29): OMP_NUM_THREADS=2, real-memory
        # cgroup (NOT ulimit -v), at most 3 ngspice at once (see run_batch()).
        e = dict(os.environ if env is None else env)
        e["OMP_NUM_THREADS"] = "2"
        cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G",
               "-p", "MemorySwapMax=0", NGSPICE, "-b", deck_path]
        proc = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=timeout, cwd=os.path.dirname(deck_path),
                               env=e)
    except FileNotFoundError as ex:
        raise CouldNotRun(f"cannot start {NGSPICE}: {ex}")
    except subprocess.TimeoutExpired:
        raise CouldNotRun(f"{deck_path} timed out after {timeout}s")
    return proc.returncode, proc.stdout + proc.stderr, time.time() - t0


# ---------------------------------------------------------------------
# trnoise NA/TS <-> target one-sided PSD.
#
# ngspice's trnoise (src/frontend/trannoise/1-f-code.c's trnoise_state_gen,
# src/spicelib/devices/isrc/isrcload.c's TRNOISE case) generates iid
# N(0, NA^2) samples on a grid spaced TS apart, and linearly interpolates
# between grid points for any time in between. Treating this as a
# pulse-amplitude-modulated (PAM) process with pulse shape h(t) = the
# normalised (peak 1) triangular kernel of base +-TS (linear interpolation
# between unit impulses IS convolution with this kernel), the standard PAM
# noise-spectrum formula gives a TWO-SIDED PSD
#   S_two(f) = (NA^2/TS) * |H(f)|^2,  H(f) = TS*sinc(f*TS)  (peak-1 triangle)
#            = NA^2 * TS * sinc(f*TS)^4
# flat (sinc^4 ~= 1) for f*TS << 1, and a ONE-SIDED PSD (f>0 convention,
# used throughout this file and by .pnoise/.noise) of
#   S_one(f) = 2 * S_two(f) = 2*NA^2*TS*sinc(f*TS)^4
# so, inverting at f=0 for the target flat-band density S0 [A^2/Hz]:
#   NA = sqrt(S0 / (2*TS))
# TS must be chosen small enough that sinc(f*TS)^4 ~= 1 up to several
# harmonics of the circuit's own clock frequency (the folding bands of
# interest) -- f*TS <= 0.01 keeps the droop under 0.02%.
def trnoise_na(S0, TS):
    return math.sqrt(S0 / (2.0 * TS))


def sinc(x):
    if x == 0.0:
        return 1.0
    px = math.pi * x
    return math.sin(px) / px


# ---------------------------------------------------------------------
# Pure-Python iterative radix-2 FFT (no numpy on this machine).

def fft(a):
    n = len(a)
    if n & (n - 1) != 0:
        raise ValueError("fft: length must be a power of 2")
    a = list(a)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        wlen = cmath.exp(-2j * math.pi / length)
        half = length // 2
        for i in range(0, n, length):
            w = 1 + 0j
            for k in range(half):
                u = a[i + k]
                v = a[i + k + half] * w
                a[i + k] = u + v
                a[i + k + half] = u - v
                w *= wlen
        length <<= 1
    return a


def hann_window(n):
    if n == 1:
        return [1.0]
    return [0.5 * (1.0 - math.cos(2.0 * math.pi * i / (n - 1))) for i in range(n)]


def segment_psd_onesided(x, dt, window, wsum2):
    """One segment's one-sided PSD estimate [unit^2/Hz], index k -> freq
    k/(N*dt), k=0..N/2."""
    N = len(x)
    fs = 1.0 / dt
    xw = [x[i] * window[i] for i in range(N)]
    X = fft([complex(v, 0.0) for v in xw])
    scale = 1.0 / (fs * wsum2)
    S = [0.0] * (N // 2 + 1)
    S[0] = (X[0].real ** 2 + X[0].imag ** 2) * scale
    for k in range(1, N // 2):
        S[k] = 2.0 * (X[k].real ** 2 + X[k].imag ** 2) * scale
    S[N // 2] = (X[N // 2].real ** 2 + X[N // 2].imag ** 2) * scale
    return S


def welch_average(samples_list, dt, nfft, overlap_frac=0.5):
    """samples_list: list of per-seed real sample lists (each a uniformly
    dt-spaced v(t) trace, startup already dropped). Segments each trace
    with a Hann window, nfft points, overlap_frac overlap, across ALL
    seeds, and averages the one-sided periodogram. Returns (freqs, S_avg,
    n_segments)."""
    window = hann_window(nfft)
    wsum2 = sum(w * w for w in window)
    step = max(1, int(nfft * (1.0 - overlap_frac)))
    acc = [0.0] * (nfft // 2 + 1)
    nseg = 0
    for samples in samples_list:
        n = len(samples)
        start = 0
        while start + nfft <= n:
            seg = samples[start:start + nfft]
            S = segment_psd_onesided(seg, dt, window, wsum2)
            for k in range(len(S)):
                acc[k] += S[k]
            nseg += 1
            start += step
    if nseg == 0:
        raise CouldNotRun("no complete FFT segments -- transient too short")
    S_avg = [v / nseg for v in acc]
    freqs = [k / (nfft * dt) for k in range(nfft // 2 + 1)]
    return freqs, S_avg, nseg


def mc_relative_spread_db(nseg):
    """Expected relative standard deviation of an nseg-segment-averaged
    periodogram bin, in dB (power domain): each segment's own periodogram
    value at a fixed bin is (for a Gaussian process) exponentially
    distributed (a chi-squared with 2 degrees of freedom, scaled) --
    relative std = 1 (100%) for ONE segment; averaging nseg such
    (approximately independent, given <=50% overlap and a decorrelating
    Hann window) segments divides the variance by nseg, so the relative
    std of the average is 1/sqrt(nseg). Converted to a +/-1-sigma dB
    spread via 10*log10(1+-1/sqrt(nseg))."""
    rel = 1.0 / math.sqrt(nseg)
    hi_db = 10.0 * math.log10(1.0 + rel)
    lo_db = -10.0 * math.log10(max(1.0 - rel, 1e-9))
    return rel, hi_db, lo_db


def body_lines(fixture):
    with open(os.path.join(FIXTURES, fixture)) as fh:
        return [l for l in fh.read().splitlines() if l.strip().lower() != ".end"]


# ---------------------------------------------------------------------
# MC replica decks.
#
# ANTI-ALIASING (card #509 third pass, root cause of the FIRST run's uniform
# +1.1..+1.7 dB "MC too high" on pnoise_r2_smooth: NOT a .pnoise error): the
# switched-RC output noise extends to ~1.6 MHz (Rsw*C1 pole while ON), and
# point-sampling v(t) at dt=1us (the old `linearize`) folds everything above
# 500 kHz into the band. Fix: the trace handed to the periodogram is the
# BOXCAR AVERAGE of v(t) over each dt sample interval, obtained inside
# ngspice as differences of integ(v-vref) sampled at dt (vref ~ the DC level:
# integrating a ~1 V level over 1e6+ steps leaves a 1/ts^2-scaling rounding floor
# of ~2 nV/rtHz in the differenced trace, found calibrating on the RC deck) (a boxcar's sinc nulls sit
# exactly on the alias frequencies n/dt). The transient itself still runs at
# maxstep = TS = 20 ns (the trnoise grid).

def _read_wrdata_real(path, integ_dt=None):
    if not os.path.exists(path):
        raise CouldNotRun(f"ngspice wrote no {path}")
    vals = []
    with open(path) as fh:
        next(fh)
        for row in fh:
            v = row.split()
            if v:
                vals.append(float(v[1]))
    if not vals:
        raise CouldNotRun(f"{path} is empty")
    if integ_dt:
        vals = [(vals[i + 1] - vals[i]) / integ_dt for i in range(len(vals) - 1)]
    return vals


def _noise_src(name, n1, n2, S0, ts):
    """TRNOISE current source in PARALLEL with the resistor (terminals exactly
    as NevalSrc); amplitude is the subckt parameter {na_<name>} so a second,
    noiseless instance (na=0) shares every time step with the noisy one."""
    return f"{name} {n1} {n2} DC 0 TRNOISE({{na*{trnoise_na(S0, ts):.9e}}} {ts:.9e} 0 0)\n"


def _deck(outdir, seed, ts, duration, dt, top, core, out_node, title):
    """Writes ONE deck holding TWO copies of the circuit (subckt `core`): Xn
    with noise (na=1) and Xd noiseless (na=0). Both live in the same
    transient, so they share every accepted time step and the difference
    v(Xn.out)-v(Xd.out) is exactly the noise response, with the large
    deterministic (periodic-orbit) waveform cancelled to rounding -- this
    replaced an earlier 'separate deterministic run' subtraction whose two
    different time grids leaked edge energy into the high band."""
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, f"mc_seed{seed}.cir")
    with open(deck, "w") as fh:
        fh.write(f"* {title} MC replica (card #509)\n")
        fh.write(top)
        fh.write(".subckt core vin clk out params: na=1\n" + core + ".ends core\n")
        fh.write("Xn vin clk outn core na=1\nXd vin clk outd core na=0\n")
        fh.write(f".options maxstep={ts:.6g}\n.tran {dt:.9e} {duration:.9e}\n.control\n"
                 f"setseed {seed}\nrun\nlet iv = integ(v(outn)-v(outd))\nlinearize iv\n"
                 "set wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                 "wrdata mc_out.txt iv\nquit\n.endc\n.end\n")
    return deck, os.path.join(outdir, "mc_out.txt")


def rc_calib_deck(outdir, seed, ts, duration, dt, r=10e3, c=1e-9, **_):
    """LTI calibration of the whole MC chain: R (noise) into C. Analytic
    4kTR/(1+(f/fc)^2)."""
    core = (f"R1 vin out {r} noisy=0\n"
            + _noise_src("Inz", "vin", "out", 4 * K_BOLTZ * T_KELVIN / r, ts)
            + f"C1 out 0 {c}\n")
    return _deck(outdir, seed, ts, duration, dt, "Vin vin 0 DC 1\nVclk clk 0 DC 0\n",
                 core, "out", "RC calibration")


def r2_smooth_mc_deck(outdir, seed, ts, duration, dt, rin=10e3, rsw=100.0,
                      c=1e-9, fclk=100e3, duty=0.5, gon=1e4, goff=1e-9, ksw=30, **_):
    """pnoise_r2_smooth.cir with noisy=0 resistors + TRNOISE in parallel (Bsw
    is the same, never-noisy, behavioural gate). Output vx."""
    core = (f"Rin vin out {rin} noisy=0\n"
            + _noise_src("Iin_noise", "vin", "out", 4*K_BOLTZ*T_KELVIN/rin, ts)
            + f"C1 out 0 {c}\nRsw out swnode {rsw} noisy=0\n"
            + _noise_src("Isw_noise", "out", "swnode", 4*K_BOLTZ*T_KELVIN/rsw, ts)
            + f"Bsw swnode 0 I = V(swnode) * ({goff}+({gon}-{goff})*0.5*(1+tanh({ksw}*(V(clk)-0.5))))\n")
    top = (f"Vin vin 0 DC 1\nVclk clk 0 PULSE(0 1 0 1n 1n {duty/fclk:.9e} {1/fclk:.9e})\n")
    return _deck(outdir, seed, ts, duration, dt, top, core, "vx", "pnoise_r2_smooth.cir")


def track_hold_mc_deck(outdir, seed, ts, duration, dt, rtrack=25e3, c=10e-12,
                       fclk=100e3, duty=0.2, gon=1.0, goff=1e-9, ksw=20, **_):
    """Track-and-hold (pnoise_track_hold_slow.cir): Rtrack noisy into a tanh
    gate onto C. Output vout."""
    core = (f"Rtrack vin trnode {rtrack} noisy=0\n"
            + _noise_src("Inz", "vin", "trnode", 4*K_BOLTZ*T_KELVIN/rtrack, ts)
            + f"Bsw trnode out I = (V(trnode)-V(out)) * ({goff}+({gon}-{goff})*0.5*(1+tanh({ksw}*(V(clk)-0.5))))\n"
            + f"C1 out 0 {c}\n")
    top = (f"Vin vin 0 DC 1 SIN(1 0.05 {fclk} 0 0)\n"
           f"Vclk clk 0 PULSE(0 1 0 5n 5n {duty/fclk:.9e} {1/fclk:.9e})\n")
    return _deck(outdir, seed, ts, duration, dt, top, core, "vout", "pnoise_track_hold_slow.cir")


def run_batch(jobs, workers=3):
    """jobs: list of (deck,outfile). Runs at most `workers` ngspice at once."""
    from concurrent.futures import ThreadPoolExecutor
    def one(j):
        rc, out, wall = run_ngspice(j[0], timeout=900)
        if rc != 0:
            raise CouldNotRun(f"{j[0]} exit {rc}: {out[-1000:]}")
        return wall
    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(one, jobs))


def run_mc(deck_builder, tag, n_seeds, ts, duration, dt, nfft, drop_start,
           **kwargs):
    """n_seeds noisy/noiseless-pair runs, 3 at a time; each boxcar trace is
    reduced, the start-up dropped, the raw file deleted, and the traces
    Welch-averaged."""
    outdir = os.path.join(DATA_ROOT, tag)
    t0 = time.time()
    jobs = [deck_builder(os.path.join(outdir, f"seed{seed}"), seed, ts,
                         duration, dt, **kwargs) for seed in range(1, n_seeds + 1)]
    run_batch(jobs)
    n_drop = int(round(drop_start / dt))
    samples_list = []
    for deck, outfile in jobs:
        vals = _read_wrdata_real(outfile, dt)
        samples_list.append(vals[n_drop:])
        os.remove(outfile)
    freqs, S_avg, nseg = welch_average(samples_list, dt, nfft, 0.5)
    return freqs, S_avg, nseg, time.time() - t0


# ---------------------------------------------------------------------
# `.pnoise` reference (the engine under test) and analytic references.

def pnoise_reference(body, outdir, out_node, in_src, cardargs, freqs_hz):
    """`.pnoise ... lin n f0 f1` over consecutive Welch bins (freqs_hz is an
    evenly spaced list); one ngspice run. Returns {freq: onoise_spectrum}."""
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, f"pnoise_{freqs_hz[0]:.6g}.cir")
    with open(deck, "w") as fh:
        fh.write(body + "\n")
        fh.write(f'.pnoise {cardargs[0]} {out_node} {in_src} lin {len(freqs_hz)} '
                 f'{freqs_hz[0]:.9g} {freqs_hz[-1]:.9g} method="shoot" solver="td"\n')
        fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\n"
                 "set numdgt=15\nwrdata pnout.txt onoise_spectrum\nquit\n.endc\n.end\n")
    rc, o, wall = run_ngspice(deck, timeout=600)
    if rc != 0:
        raise CouldNotRun(f".pnoise exit {rc}: {o[-1500:]}")
    v = _read_wrdata_real(os.path.join(outdir, "pnout.txt"))
    if len(v) != len(freqs_hz):
        raise CouldNotRun(f"expected {len(freqs_hz)} points, got {len(v)}")
    return dict(zip(freqs_hz, v))


def nearest_bin_freqs(dt, nfft, targets):
    df = 1.0 / (nfft * dt)
    return [max(1, round(t / df)) * df for t in targets]


def sinc_droop(f, dt):
    """Boxcar (dt) magnitude^2 droop in the MC trace, divided out."""
    return sinc(f * dt) ** 2


def compare(tag, deck_builder, body, out_node, cardargs, ref_fn, n_seeds, ts,
            duration, dt, nfft, drop_start, target_freqs, ref_label, **kwargs):
    """Band-averaged comparison: for each target, NB consecutive Welch bins
    are averaged in POWER (MC) and compared with the mean of the reference
    evaluated at exactly those bins. Expected 1-sigma spread of a band
    average = 4.34/sqrt(nseg*NB) dB (independent bins, chi2 argument)."""
    NB = 9
    print(f"\n=== {tag} ===")
    freqs, S_avg, nseg, wall = run_mc(deck_builder, tag, n_seeds, ts, duration,
                                      dt, nfft, drop_start, **kwargs)
    binw = 1.0 / (nfft * dt)
    with open(os.path.join(DATA_ROOT, tag, "S_avg.txt"), "w") as fh:
        for k in range(len(S_avg)):
            fh.write(f"{k*binw:.6f} {S_avg[k]:.9e}\n")
    rel, hi1, lo1 = mc_relative_spread_db(nseg)
    # adjacent Hann bins are correlated (power correlation 0.444 at lag 1,
    # 0.028 at lag 2) and 50%-overlapped segments add ~5.6 % variance, so the
    # NB-bin band average has n_eff = nseg*NB/(1.83*1.056) independent chi2_2
    # samples, not nseg*NB.
    n_eff = nseg * NB / (1.83 * 1.056)
    relb, hib, lob = mc_relative_spread_db(n_eff)
    print(f"MC: {n_seeds} seeds x {duration*1e3:.1f} ms, boxcar dt={dt*1e9:.0f} ns, "
          f"nfft={nfft} ({binw:.1f} Hz bins), {nseg} Hann/50%-overlap Welch segments;"
          f" expected 1-sigma: single bin +{hi1:.2f}/-{lo1:.2f} dB, {NB}-bin "
          f"band +{hib:.3f}/-{lob:.3f} dB; wall {wall:.0f}s")
    print(f"{'band centre':>11s} {'MC V/rtHz':>12s} {ref_label:>12s} {'d(dB)':>8s} {'d/sigma':>8s}   "
          f"| centre bin: MC {'ref':>12s} d(dB)")
    rows = []
    for t in target_freqs:
        k0 = max(1, round(t / binw) - NB // 2)
        ks = list(range(k0, k0 + NB))
        fl = [k * binw for k in ks]
        ref = ref_fn(fl)
        mcp = sum(S_avg[k] / sinc_droop(k * binw, dt) for k in ks) / NB
        refp = sum(ref[f] ** 2 for f in fl) / NB
        d = 10 * math.log10(mcp / refp)
        sig = hib if d > 0 else lob
        kc = ks[NB // 2]
        dc = 10 * math.log10(S_avg[kc] / sinc_droop(kc * binw, dt) / ref[fl[NB // 2]] ** 2)
        rows.append((fl[NB // 2], math.sqrt(mcp), math.sqrt(refp), d, d / sig))
        print(f"{fl[NB//2]:11.1f} {math.sqrt(mcp):12.5e} {math.sqrt(refp):12.5e} "
              f"{d:+8.3f} {d/sig:+8.2f}   | {dc:+.3f}")
    print(f"worst |d| = {max(abs(r[3]) for r in rows):.3f} dB; mean d = "
          f"{sum(r[3] for r in rows)/len(rows):+.3f} dB (mean 1-sigma "
          f"{hib/math.sqrt(len(rows)):.3f} dB)")
    return rows, (lob, hib), nseg


def th_analytic_S0(rtrack, c, T, tt):
    """Exact f->0 one-sided PSD of the held output of the track-and-hold
    (ideal switch limit: Vin noise-free, Rtrack the only noise source, hold
    perfectly isolating). x(t) is an Ornstein-Uhlenbeck process (variance
    kT/C, time constant tau=Rtrack*C) restarted from the held value at each
    track start, run for tt, then frozen for T-tt. Per-cycle integral
    Z_n = A_n + Th*s_n, A_n = int of the OU piece; long-run variance
    Var(Z)+2Cov(Z_n,Z_n+1)+2Cov(Z_n,Z_n+2); S(0)=2*LRV/T."""
    kT = K_BOLTZ * T_KELVIN
    tau = rtrack * c
    a = tt / tau
    Th = T - tt
    q = 1 - math.exp(-a)
    rho = math.exp(-a)
    v = kT / c
    varA = 2 * tau * (tt - tau * q)
    covAs = tau * q
    var0 = varA + Th * Th + 2 * Th * covAs
    # lag 1 (see docstring): A_{n+1} = tau*q*s_n + indep; s_{n+1} = rho*s_n + indep
    c1 = (tau * q) * covAs + Th * rho * covAs + Th * tau * q + Th * Th * rho
    # lag 2: everything picks up another factor rho on the s-chain
    c2 = (tau * q) * covAs * rho + Th * rho * rho * covAs + Th * tau * q * rho + Th * Th * rho * rho
    return 2 * v * (var0 + 2 * c1 + 2 * c2) / T


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    nseeds = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    TS = float(sys.argv[3]) * 1e-9 if len(sys.argv) > 3 else 20e-9
    sfx = "" if TS == 20e-9 else f"_ts{TS*1e9:g}"
    TARGETS = [4000, 10000, 16000, 22000, 28000, 35000, 43000]
    dt, nfft, dur, drop = 200e-9, 8192, 20e-3, 2e-3
    if what in ("rc", "all"):
        r, c = 10e3, 1e-9
        fc = 1 / (2 * math.pi * r * c)
        an = lambda fl: {f: math.sqrt(4*K_BOLTZ*T_KELVIN*r/(1+(f/fc)**2)) for f in fl}
        compare("rc_calib"+sfx, rc_calib_deck, None, "vx", None, an, nseeds, TS,
                dur, dt, nfft, drop, TARGETS, "analytic", r=r, c=c)
    if what in ("r2", "all"):
        body = "\n".join(body_lines("pnoise_r2_smooth.cir"))
        ca = ("100e3 100e-6 vx 2000 1 50 1e-6",)
        ref = lambda fl: pnoise_reference(body, os.path.join(DATA_ROOT, "r2_smooth", "pnoise_ref"), "vx", "Vin", ca, fl)
        compare("r2_smooth"+sfx, r2_smooth_mc_deck, body, "vx", ca, ref, nseeds,
                TS, dur, dt, nfft, drop, TARGETS, ".pnoise")
    if what in ("th", "all"):
        body = "\n".join(body_lines("pnoise_track_hold_slow.cir"))
        ca = ("100e3 100e-6 vout 2000 1 50 1e-6",)
        ref = lambda fl: pnoise_reference(body, os.path.join(DATA_ROOT, "th_slow", "pnoise_ref"), "vout", "Vin", ca, fl)
        compare("track_hold_pnoise"+sfx, track_hold_mc_deck, body, "vout", ca, ref,
                nseeds, TS, dur, dt, nfft, drop, TARGETS, ".pnoise")
    if what in ("tha", "all"):
        body = "\n".join(body_lines("pnoise_track_hold_slow.cir"))
        ca = ("100e3 100e-6 vout 2000 1 50 1e-6",)
        T = 1e-5
        S0 = th_analytic_S0(25e3, 10e-12, T, 0.2 * T)
        naive = 2 * K_BOLTZ * T_KELVIN / 10e-12 * T
        print(f"\nanalytic S(0): exact OU-restart {math.sqrt(S0):.5e} V/rtHz; "
              f"naive ZOH kT/C {math.sqrt(naive):.5e}")
        print("(.pnoise, sinc^2(fT)-shaped; only f<=3kHz compared: shape corrections there < 0.01 dB)")
        fl = [300.0, 1000.0, 2000.0, 3000.0]
        pn = {}
        for f in fl:
            pn.update(pnoise_reference(body, os.path.join(DATA_ROOT, "th_slow", f"pnoise_ana{f:g}"),
                                       "vout", "Vin", ca, [f]))
        for f in fl:
            an = math.sqrt(S0 * sinc(f * T) ** 2)
            print(f"{f:8.0f} Hz  .pnoise {pn[f]:.6e}  analytic {an:.6e}  {20*math.log10(pn[f]/an):+.3f} dB")


if __name__ == "__main__":
    main()
