#!/usr/bin/env python3
"""tools/pssac_ref.py -- THIRD, independent small-signal method for
periodically-driven (clocked) circuits, card #484 ("settle the near-fs/2
disagreement between `.pac solver=td` and the transient-FRA reference with
a THIRD, independent method"). Private fork, no upstream contact, smooth
models only (project rulings #446 B / #447 A, same restriction every other
PSS/PAC/PSTB tool in this project already carries).

METHOD ("perturbed periodic steady state", PLECS/Spectre AC-sweep style --
NOT Spectre's/this fork's own linearised PAC, and NOT FRA's causal
transient-plus-window technique; see the module docstring in
tests/test_pac.py for what those two do):

  1. Pick a target frequency f. Since this fork's `.pss method=shoot`
     solves for a period FIXED at 1/Fguess (a driven circuit, not an
     autonomous one -- docs/pss-pac.md, "Fguess ... the deck's own known
     clock frequency, not a guess to be refined"), f must be commensurate
     with the switching frequency fs for a finite-period shooting solve to
     exist at all: f = fs*p/q for small integers p, q. rational_approx()
     below finds the (p, q) with smallest relative frequency error over
     q = 1..qmax (never assumed exact -- the achieved relative error is
     always reported, and the comparison is made at the ACTUAL achieved
     frequency fs*p/q, not silently mislabelled as the disputed point).

  2. Two shooting-PSS solves, same LCM period T = q/fs (Fguess = fs/q,
     so Points = q * points_per_cycle keeps the same within-a-clock-cycle
     resolution the project's other PSS decks use):
       - baseline: the circuit exactly as in tests/fixtures/pac_*.cir,
         unperturbed (the existing AC-flagged injection element's SIN
         content, if any, removed -- DC only).
       - perturbed: the SAME circuit with a SMALL sinusoidal perturbation
         SIN(offset A f 0 0 0) added to the injection element/source at
         the disputed frequency f = fs*p/q.
     Both converge to their own genuine periodic steady state via this
     fork's shooting-Newton PSS (src/spicelib/analysis/pssshoot.c) -- the
     SAME engine `.pac`/`.pstb` build on, but used here only to reach a
     converged nonlinear orbit, never to linearise one. No perturb-and-
     subtract-around-a-*transient* windowing (FRA's technique) is
     involved: both orbits are each individually converged to a periodic
     steady state before anything is subtracted.

  3. orbit_diff(t) = perturbed(t) - baseline(t) over one LCM period,
     resampled onto a uniform grid (this fork's PSS recording pass uses
     ordinary adaptive DCtran() output -- see pssshoot.c's header, "the
     converged trajectory is emitted as one more DCtran() run" -- so the
     accepted time grid is NOT uniform; resampled here by linear
     interpolation, independent of pacshoot.c's own internal recording-
     pass grid machinery, which this tool never touches or imports).

  4. Sideband-0 (same-frequency) component: a single-bin DFT (Goertzel-
     style, O(N) for the one bin needed) of orbit_diff at exactly bin
     p (frequency p/T = p*fs/q = f by construction), DIVIDED by the same
     single-bin DFT of the analytically-known injected difference
     waveform A*sin(2*pi*f*t) sampled on the identical grid -- this
     ratio is exactly a phasor gain and cancels the sin/cos and 1/2
     transform-scaling conventions without needing to track them by hand.

  5. Linearity check: repeat step 2-4 at amplitude A/2; the two H(f)
     estimates must agree within 0.05 dB (task's own bound) or the
     "perturbation" was not actually small-signal.

Validated (validate_r2, validate_lti below) against refs/r2_analytic.g0()
(exact, closed-form) and against plain `.ac` on an LTI deck before being
trusted on the disputed R3/R4 points -- see the card #484 report for the
numbers actually measured, not asserted here.

Usage:
    python3 tools/pssac_ref.py validate-lti
    python3 tools/pssac_ref.py validate-r2 [--freqs 1000,5000,20000]
    python3 tools/pssac_ref.py point --tag r3_line2out --freq 199054 \\
        [--amplitude 0.01] [--qmax 40] [--vin-v 12]
    python3 tools/pssac_ref.py point --tag r4_line2out --freq 79244.7 --vin-v 11

    NGSPICE=<prefix>/bin/ngspice-aug python3 tools/pssac_ref.py ...
"""
import argparse
import cmath
import math
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES = os.path.join(ROOT, "tests", "fixtures")
REFS = os.path.join(ROOT, "refs")
NGSPICE = os.environ.get("NGSPICE", "ngspice")
# card #484 provenance fix (parent flagged 2026-09-29: a stray helper process
# from the #479/deadtime work ran this tool repeatedly, concurrently,
# and could have written into these same tag/freq/amplitude-keyed output
# directories): every invocation of this file gets its own RUN_ID
# (PID + start time), included in DATA_ROOT, so two processes can never
# share -- let alone silently overwrite or race on -- the same deck/orbit
# files, and a stale or foreign run can never be mistaken for this one's.
# Override with PSSAC_REF_RUN_ID= for a caller that wants a stable,
# reproducible path across repeated invocations of ITS OWN (e.g. a test
# suite); otherwise every run is isolated by construction.
RUN_ID = os.environ.get("PSSAC_REF_RUN_ID") or f"pid{os.getpid()}_{int(time.time())}"
DATA_ROOT = os.path.expanduser(f"~/Work/data/ngspice_augmentation/pssac_ref/{RUN_ID}")


class CouldNotRun(Exception):
    pass


# ---------------------------------------------------------------------
# Generic numeric helpers (self-contained -- no import of tools/fra.py or
# any pacshoot.c-adjacent test helper; see module docstring on why).

def rational_approx(f_target, fs, qmax=64, tol=0.006):
    """(p, q, f_actual, relerr) with f_actual = fs*p/q, 1 <= q <= qmax.

    Prefers the SMALLEST q whose relative error is under `tol` (default
    0.6%, comfortably tighter than every disputed point's own td-PAC-vs-
    FRA gap this card investigates) -- a large q means a q-clock-period
    LCM shooting solve (Points = q*points_per_cycle, and a q-times-longer
    finite-difference monodromy period-map per Newton iteration), so a
    small q that is accurate enough is strictly preferable to the
    global-minimum-error q, which is often needlessly large (e.g. q=63 at
    0.15% error vs q=19 at 0.375% for one of this card's own disputed
    points -- both easily inside FRA/td-PAC's own multi-percent
    disagreement). Falls back to the global minimum-error q found within
    qmax if none clears `tol` (never silently fails; relerr is always
    returned so a caller can judge whether the fallback is good enough)."""
    best = None
    for q in range(1, qmax + 1):
        p = round(f_target / fs * q)
        if p <= 0:
            continue
        f_actual = fs * p / q
        relerr = abs(f_actual - f_target) / f_target
        if relerr < tol:
            return (p, q, f_actual, relerr)
        if best is None or relerr < best[3]:
            best = (p, q, f_actual, relerr)
    return best


def run_ngspice(deck_path, timeout=1800, env=None):
    t0 = time.time()
    try:
        proc = subprocess.run([NGSPICE, "-b", deck_path], capture_output=True,
                               text=True, timeout=timeout,
                               cwd=os.path.dirname(deck_path), env=env)
    except FileNotFoundError as e:
        raise CouldNotRun(f"cannot start {NGSPICE}: {e}")
    except subprocess.TimeoutExpired:
        raise CouldNotRun(f"{deck_path} timed out after {timeout}s")
    return proc.returncode, proc.stdout + proc.stderr, time.time() - t0


def read_wrdata_complex(path, ncols=1):
    """freqs, [ncols complex columns] from a `wr_singlescale`/`wr_vecnames`
    AC-type wrdata file (freq, then real/imag pairs -- NOT the same layout
    as a real time-domain wrdata trace, which is (t, v) pairs per vector;
    see read_wrdata_pair() below for that one). Same format
    tests/test_pac.py's read_wrdata_complex() reads; reimplemented here
    standalone (this file imports nothing from tests/)."""
    freqs, cols = [], [[] for _ in range(ncols)]
    with open(path) as fh:
        next(fh)
        for row in fh:
            v = row.split()
            if not v:
                continue
            v = [float(x) for x in v]
            if len(v) != 1 + 2 * ncols:
                raise CouldNotRun(f"unexpected column count in {path}: {row!r}")
            freqs.append(v[0])
            for k in range(ncols):
                cols[k].append(complex(v[1 + 2 * k], v[2 + 2 * k]))
    if not freqs:
        raise CouldNotRun(f"{path} is empty")
    return freqs, cols


def read_wrdata_pair(path, col_index=0):
    """(t, v) for the col_index-th (0-based) requested vector in a wrdata
    ascii file (each vector contributes a (t, v) column pair)."""
    ts, vs = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            need = 2 * (col_index + 1)
            if len(parts) < need:
                continue
            try:
                ts.append(float(parts[2 * col_index]))
                vs.append(float(parts[2 * col_index + 1]))
            except ValueError:
                continue
    return ts, vs


def resample_uniform(t, v, T, n):
    """Linear-interpolate (t, v) (t increasing, covering >= [0, T]) onto n
    uniform samples t_k = k*T/n, k=0..n-1."""
    if not t or t[0] > 1e-3 * T:
        raise CouldNotRun(f"orbit trace does not start near t=0 (t[0]={t[0] if t else None}, "
                           f"T={T})")
    if t[-1] < T * 0.99:
        raise CouldNotRun(f"orbit trace does not cover the full period "
                           f"(t[-1]={t[-1]}, T={T})")
    out = []
    j = 0
    nt = len(t)
    for k in range(n):
        tk = T * k / n
        if tk <= t[0]:
            out.append(v[0])
            continue
        if tk >= t[-1]:
            # trailing gap before T: interpolate to v[0] assuming exact
            # periodicity v(T) == v(0) (the PSS convergence criterion
            # itself), rather than a flat extrapolation of v[-1].
            if T <= t[-1]:
                out.append(v[-1])
            else:
                frac = (tk - t[-1]) / (T - t[-1])
                out.append(v[-1] + frac * (v[0] - v[-1]))
            continue
        while j + 1 < nt - 1 and t[j + 1] < tk:
            j += 1
        t0, t1 = t[j], t[j + 1]
        v0, v1 = v[j], v[j + 1]
        if t1 <= t0:
            out.append(v0)
        else:
            frac = (tk - t0) / (t1 - t0)
            out.append(v0 + frac * (v1 - v0))
    return out


def resample_window(t, v, t0, T, n):
    """Linear-interpolate (t, v) (t increasing) onto n uniform samples over
    the window [t0, t0+T) -- used by the long-tran fallback backend (see
    run_tran_orbit()) to pull an arbitrary later window out of a long
    trace, unlike resample_uniform() above which always starts at t=0 (a
    single-period `.pss` orbit)."""
    import bisect
    out = []
    for k in range(n):
        tk = t0 + T * k / n
        i = bisect.bisect_right(t, tk) - 1
        if i < 0:
            out.append(v[0])
            continue
        if i >= len(t) - 1:
            out.append(v[-1])
            continue
        t0_, t1_ = t[i], t[i + 1]
        v0_, v1_ = v[i], v[i + 1]
        if t1_ <= t0_:
            out.append(v0_)
        else:
            frac = (tk - t0_) / (t1_ - t0_)
            out.append(v0_ + frac * (v1_ - v0_))
    return out


def dft_bin(y, p):
    """Single-bin complex DFT X_p = (1/N) sum_k y_k * exp(-2j*pi*p*k/N),
    O(N) direct summation (Goertzel-equivalent; N is a few thousand here,
    an FFT would be overkill for one bin)."""
    n = len(y)
    two_pi_p_over_n = 2.0 * math.pi * p / n
    s = 0j
    for k, yk in enumerate(y):
        s += yk * cmath.exp(-1j * two_pi_p_over_n * k)
    return s / n


def db_deg(h):
    return 20.0 * math.log10(abs(h)), math.degrees(cmath.phase(h))


# ---------------------------------------------------------------------
# Deck construction: load a fixture/ref body (everything up to its own
# .control/.end), substitute one injection line for the baseline/perturbed
# variant, append a plain `.pss ... method=shoot` card + a control block
# that writes the one-period orbit via wrdata.

def _load_body(path, stop_at_control=True):
    with open(path) as fh:
        lines = fh.read().splitlines()
    out = []
    for l in lines:
        if stop_at_control and l.strip().lower() == ".control":
            break
        if l.strip().lower() == ".end":
            continue
        out.append(l)
    return out


def _substitute_line(body, pattern, replacement, what):
    out = []
    hit = 0
    rx = re.compile(pattern)
    for l in body:
        if rx.match(l.strip()):
            out.append(replacement)
            hit += 1
        else:
            out.append(l)
    if hit != 1:
        raise CouldNotRun(f"{what}: pattern {pattern!r} matched {hit} lines (expected 1)")
    return out


def _substitute_param(body, param, value, what):
    out = []
    hit = 0
    rx = re.compile(rf"^\.param\s+.*\b{re.escape(param)}=[0-9.eE+-]+")
    for l in body:
        if rx.match(l.strip()):
            new, n = re.subn(rf"({re.escape(param)}=)[0-9.eE+-]+", rf"\g<1>{value}", l, count=1)
            if n != 1:
                raise CouldNotRun(f"{what}: {param} substitution regex did not match inside {l!r}")
            out.append(new)
            hit += 1
        else:
            out.append(l)
    if hit != 1:
        raise CouldNotRun(f"{what}: {param}= line matched {hit} times (expected 1)")
    return out


def build_pss_deck(body, f_lcm, stabtime, oscnode, points, out_expr, outdir, tag,
                    sc_iter=150, steady_coeff=1e-6, harmonics=10):
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, f"{tag}.cir")
    orbit = os.path.join(outdir, f"{tag}_orbit.txt")
    with open(deck, "w") as fh:
        fh.write("\n".join(body) + "\n")
        fh.write(f'.control\nset filetype=ascii\nset wr_singlescale\nset numdgt=17\n'
                  f'pss {f_lcm!r} {stabtime} {oscnode} {points} {harmonics} {sc_iter} '
                  f'{steady_coeff} uic method=shoot\n'
                  f'wrdata {os.path.basename(orbit)} {out_expr}\n'
                  f'.endc\n.end\n')
    return deck, orbit


def build_tran_deck(body, stabtime_s, T, tstep, out_expr, outdir, tag):
    """Long-transient fallback backend (task's own contingency: 'if [.pss]
    doesn't accept a forced period, use a long settled .tran over many LCM
    periods with tight tolerances and a checked periodicity residual').

    Used here for a DIFFERENT, empirically-found reason, not because `.pss
    method=shoot` refuses a forced q-cycle period -- it accepts one and
    Newton genuinely converges (small residual, stable spectral radius) --
    but because this fork's shooting engine has a real bug in its
    TRAJECTORY-EMISSION pass specifically when the requested PSS period is
    a MULTIPLE of the circuit's own fastest periodic source (q > 1): the
    converged x0 is fine, but the final reporting `DCtran()` re-run
    (pssshoot.c's own "the converged trajectory is emitted as one more
    DCtran() run") diverges partway through for some q (e.g. q=5 at
    199054 Hz on R3/R4 measured here -- see the card #484 report), leaving
    physically-impossible values (~1e8 V) in the tail of the emitted
    orbit -- reproducible even on the UNPERTURBED baseline circuit, so
    this is an engine issue with q>1 emission, not an artefact of this
    tool's perturbation. measure_h() below detects this (implausible
    |v| in the orbit) and falls back to this function automatically.

    Runs .tran uic for stabtime_s + 2*T (T = the LCM period, q clock
    cycles), `.options trtol=1` (closer-to-uniform accepted steps, same
    knob tools/fra.py's own transfer-mode settings use), tstep = the same
    per-clock-cycle resolution the .pss backend targets. The CALLER reads
    back two full-T windows ([stabtime, stabtime+T] and [stabtime+T,
    stabtime+2T]) and checks their own periodicity residual before
    trusting either as "the" periodic orbit -- see measure_h()."""
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, f"{tag}.cir")
    orbit = os.path.join(outdir, f"{tag}_orbit.txt")
    tstop = stabtime_s + 2.0 * T
    with open(deck, "w") as fh:
        fh.write("\n".join(body) + "\n")
        fh.write(f".options trtol=1\n")
        fh.write(f".tran {tstep!r} {tstop!r} uic\n")
        fh.write(f'.control\nset filetype=ascii\nset wr_singlescale\nset numdgt=17\n'
                  f'run\nwrdata {os.path.basename(orbit)} {out_expr}\n'
                  f'.endc\n.end\n')
    return deck, orbit, tstop


# ---------------------------------------------------------------------
# Per-tag injection descriptions for tests/fixtures/pac_r{3,4}_*.cir
# (splice points documented in each fixture's own header, transcribed in
# this file's module docstring for r3/r4; identical injection points
# tools/fra_curves_fork.py's instrument_r3()/instrument_r4() use for the
# FRA reference itself and pacshoot.c's own AC-flagged-source convention
# uses for `.pac solver=td`, so this is a genuinely independent RUN of the
# same physical injection point, not a different circuit).
TAGS = {
    "r3_ctrl2out": dict(fixture="pac_r3_ctrl2out.cir",
                         line_pat=r"Vctrl\s+ctrlnode\s+opout\s+DC\s+0(\s+AC\s+1)?$",
                         base="Vctrl ctrlnode opout DC 0",
                         pert=lambda A, f: f"Vctrl ctrlnode opout DC 0 SIN(0 {A:.10g} {f:.10g} 0 0 0)",
                         offset_expr="0", oscnode="vout", out_expr="v(vout)",
                         stabtime="300e-6", fs=500e3, vin_v_param=None),
    "r3_line2out": dict(fixture="pac_r3_line2out.cir",
                         line_pat=r"Vin\s+vin\s+0\s+DC\s+\{vin_v\}(\s+AC\s+1)?$",
                         base="Vin vin 0 DC {vin_v}",
                         pert=lambda A, f: f"Vin vin 0 DC {{vin_v}} SIN({{vin_v}} {A:.10g} {f:.10g} 0 0 0)",
                         offset_expr="{vin_v}", oscnode="vout", out_expr="v(vout)",
                         stabtime="300e-6", fs=500e3, vin_v_param="vin_v", vin_v_default=12),
    "r3_zout": dict(fixture="pac_r3_zout.cir",
                     line_pat=r"Iinj\s+vout\s+0\s+DC\s+0(\s+AC\s+1)?$",
                     base="Iinj vout 0 DC 0",
                     pert=lambda A, f: f"Iinj vout 0 DC 0 SIN(0 {A:.10g} {f:.10g} 0 0 0)",
                     offset_expr="0", oscnode="vout", out_expr="v(vout)",
                     stabtime="300e-6", fs=500e3, vin_v_param=None),
    "r4_ctrl2out": dict(fixture="pac_r4_ctrl2out.cir",
                         line_pat=r"Vctrl\s+ctrlnode\s+vc\s+DC\s+0(\s+AC\s+1)?$",
                         base="Vctrl ctrlnode vc DC 0",
                         pert=lambda A, f: f"Vctrl ctrlnode vc DC 0 SIN(0 {A:.10g} {f:.10g} 0 0 0)",
                         offset_expr="0", oscnode="vout", out_expr="v(vout)",
                         stabtime="3500e-6", fs=500e3, vin_v_param="vin_v", vin_v_default=11),
    "r4_line2out": dict(fixture="pac_r4_line2out.cir",
                         line_pat=r"Vin\s+vin\s+0\s+DC\s+\{vin_v\}(\s+AC\s+1)?$",
                         base="Vin vin 0 DC {vin_v}",
                         pert=lambda A, f: f"Vin vin 0 DC {{vin_v}} SIN({{vin_v}} {A:.10g} {f:.10g} 0 0 0)",
                         offset_expr="{vin_v}", oscnode="vout", out_expr="v(vout)",
                         stabtime="3500e-6", fs=500e3, vin_v_param="vin_v", vin_v_default=11),
    "r4_zout": dict(fixture="pac_r4_zout.cir",
                     line_pat=r"Iinj\s+vout\s+0\s+DC\s+0(\s+AC\s+1)?$",
                     base="Iinj vout 0 DC 0",
                     pert=lambda A, f: f"Iinj vout 0 DC 0 SIN(0 {A:.10g} {f:.10g} 0 0 0)",
                     offset_expr="0", oscnode="vout", out_expr="v(vout)",
                     stabtime="3500e-6", fs=500e3, vin_v_param="vin_v", vin_v_default=11),
    # card #494 SECOND pass (plan T6.1): open-loop ideal-switch buck
    # (Shs/Sls direct differential `sw` control, no compensator), the
    # saltation baseline-failure fixture (tests/fixtures/pac_r6o_ctrl2out.cir
    # / pss_r6o_shoot.cir). Vctrl is the SAME independent source both the
    # `.pac solver=td` run and the DC finite-difference check (card #494
    # report) inject at -- this TAG gives a THIRD, genuinely independent
    # small-signal number (perturbed-PSS + single-bin DFT, not a per-step
    # linearisation and not a DC-only finite difference) at the SAME
    # disputed 1kHz point, per the card's own "build TWO independent
    # references and agree them first" instruction.
    "r6o_ctrl2out": dict(fixture="pac_r6o_ctrl2out.cir",
                          line_pat=r"Vctrl\s+ctrlnode\s+0\s+DC\s+\{vctrl_dc\}(\s+AC\s+1)?$",
                          base="Vctrl ctrlnode 0 DC {vctrl_dc}",
                          pert=lambda A, f: f"Vctrl ctrlnode 0 DC {{vctrl_dc}} SIN({{vctrl_dc}} {A:.10g} {f:.10g} 0 0 0)",
                          offset_expr="{vctrl_dc}", oscnode="vout", out_expr="v(vout)",
                          stabtime="100e-6", fs=500e3, vin_v_param="vin_v", vin_v_default=12),
}


def _tag_body(tag_def, vin_v=None):
    body = _load_body(os.path.join(FIXTURES, tag_def["fixture"]))
    if vin_v is not None:
        if tag_def["vin_v_param"] is None:
            raise CouldNotRun(f"{tag_def['fixture']}: --vin-v given but this fixture has no vin_v param")
        body = _substitute_param(body, "vin_v", vin_v, tag_def["fixture"])
    return body


def _orbit_plausible(v):
    """Heuristic sanity check for a `.pss` emission-pass orbit: this
    project's circuits (R2/R3/R4/LTI) never see node voltages/currents
    outside a low-volt/low-amp range by design, so a multi-order-of-
    magnitude blowup (found empirically at q>1 -- see build_tran_deck()'s
    docstring) is unambiguous corruption, not a legitimate large value."""
    if not v:
        return False
    n = len(v)
    head = v[: max(1, int(n * 0.8))]  # exclude a possible corrupted tail
    med = sorted(abs(x) for x in head)[len(head) // 2]
    worst = max(abs(x) for x in v)
    return worst < 1e4 and worst < 1000.0 * max(med, 1e-9)


def _measure_h_pss(tag, td, f_actual, p, q, amplitude, points_per_cycle, sc_iter,
                    steady_coeff, base_body, pert_body, outdir):
    f_lcm = td["fs"] / q
    points = q * points_per_cycle
    T = 1.0 / f_lcm

    base_deck, base_orbit = build_pss_deck(base_body, f_lcm, td["stabtime"], td["oscnode"],
                                            points, td["out_expr"], outdir, "baseline",
                                            sc_iter=sc_iter, steady_coeff=steady_coeff)
    pert_deck, pert_orbit = build_pss_deck(pert_body, f_lcm, td["stabtime"], td["oscnode"],
                                            points, td["out_expr"], outdir, "perturbed",
                                            sc_iter=sc_iter, steady_coeff=steady_coeff)
    t_wall = 0.0
    for deck, label_ in ((base_deck, "baseline"), (pert_deck, "perturbed")):
        rc, out, wall = run_ngspice(deck, timeout=1800)
        t_wall += wall
        if rc != 0:
            raise CouldNotRun(f"{tag} {label_} .pss exit {rc} (f={f_actual:.3f}Hz q={q}): "
                               f"{out[-2000:]}")
        if "FAILED to converge" in out:
            fline = next((l for l in out.splitlines() if "FAILED to converge" in l), "")
            raise CouldNotRun(f"{tag} {label_} PSS-SHOOT did not converge (f={f_actual:.3f}Hz "
                               f"q={q}): {fline}")

    tb, vb = read_wrdata_pair(base_orbit, 0)
    tp, vp = read_wrdata_pair(pert_orbit, 0)
    if not (_orbit_plausible(vb) and _orbit_plausible(vp)):
        raise CouldNotRun(f"{tag}: .pss trajectory-emission pass produced an implausible "
                           f"orbit (q={q} multi-clock-cycle emission bug -- see "
                           f"build_tran_deck()'s docstring); falling back to the long-tran "
                           f"backend")
    N = points
    vb_u = resample_uniform(tb, vb, T, N)
    vp_u = resample_uniform(tp, vp, T, N)
    diff = [b - a for a, b in zip(vb_u, vp_u)]
    inj = [amplitude * math.sin(2.0 * math.pi * f_actual * (T * k / N)) for k in range(N)]
    return diff, inj, N, t_wall, "pss", None


def _measure_h_tran(tag, td, f_actual, p, q, amplitude, points_per_cycle,
                     base_body, pert_body, outdir, stabtime_mult=10.0):
    """stabtime_mult (default 10x the tag's own `.pss`-tuned stabtime):
    plain `.tran` settling (no Newton acceleration) needs materially more
    real time than shooting-Newton PSS to bring the RIPPLE-RELATIVE
    periodicity residual (not just the DC operating point, which settles
    much faster) down near this backend's own target signal scale --
    measured empirically this card (r3_line2out @200kHz, q=5): 300us (the
    tag's native `.pss` stabtime) left an 11-13% ripple-relative residual,
    1000us ~2%, still not tight enough at the ~1e-4-relative-to-ripple
    level a -70dB sideband needs. See the card #484 report for the actual
    numbers measured and the honesty note on the residual this backend is
    left carrying at its used settings."""
    fs = td["fs"]
    f_lcm = fs / q
    T = 1.0 / f_lcm
    tclk = 1.0 / fs
    tstep = tclk / points_per_cycle
    N = q * points_per_cycle
    stabtime_s = float(td["stabtime"]) * stabtime_mult

    base_deck, base_orbit, tstop = build_tran_deck(base_body, stabtime_s, T, tstep,
                                                     td["out_expr"], outdir, "baseline_tran")
    pert_deck, pert_orbit, _ = build_tran_deck(pert_body, stabtime_s, T, tstep,
                                                td["out_expr"], outdir, "perturbed_tran")
    t_wall = 0.0
    windows = {}
    for deck, orbit, label_ in ((base_deck, base_orbit, "baseline"),
                                 (pert_deck, pert_orbit, "perturbed")):
        rc, out, wall = run_ngspice(deck, timeout=3600)
        t_wall += wall
        if rc != 0:
            raise CouldNotRun(f"{tag} {label_} .tran exit {rc} (f={f_actual:.3f}Hz q={q}): "
                               f"{out[-2000:]}")
        t, v = read_wrdata_pair(orbit, 0)
        if not t or t[-1] < tstop * 0.999:
            raise CouldNotRun(f"{tag} {label_} .tran did not reach tstop "
                               f"(t[-1]={t[-1] if t else None}, tstop={tstop})")
        w1 = resample_window(t, v, stabtime_s, T, N)
        w2 = resample_window(t, v, stabtime_s + T, T, N)
        ripple = max(w2) - min(w2)
        resid = max(abs(b - a) for a, b in zip(w1, w2))
        resid_rel = resid / ripple if ripple > 1e-15 else resid
        windows[label_] = dict(w1=w1, w2=w2, resid_rel=resid_rel)

    diff = [b - a for a, b in zip(windows["baseline"]["w2"], windows["perturbed"]["w2"])]
    t_start = stabtime_s + T
    inj = [amplitude * math.sin(2.0 * math.pi * f_actual * (t_start + T * k / N))
           for k in range(N)]
    periodicity = dict(baseline_resid_rel=windows["baseline"]["resid_rel"],
                        perturbed_resid_rel=windows["perturbed"]["resid_rel"])
    return diff, inj, N, t_wall, "tran", periodicity


def measure_h(tag, freq, amplitude, qmax=64, points_per_cycle=200, sc_iter=200,
              steady_coeff=1e-7, vin_v=None, outdir_root=None, label="", backend="auto",
              tran_stabtime_mult=10.0):
    """One H(freq) estimate for TAGS[tag]: baseline orbit (recomputed each
    call, not cached) + one perturbed orbit at `amplitude`, backend='pss'
    (`.pss method=shoot`, one LCM period) or 'tran' (long transient, see
    build_tran_deck()'s docstring) or 'auto' (try pss first, fall back to
    tran if the emission pass is found corrupted -- see _measure_h_pss()).
    Returns a dict with f_actual, p, q, relerr, H, db, deg, backend_used,
    and (tran backend only) a periodicity residual (relative to that
    window's own ripple) for both baseline and perturbed orbits."""
    td = TAGS[tag]
    fs = td["fs"]
    p, q, f_actual, relerr = rational_approx(freq, fs, qmax=qmax)

    body = _tag_body(td, vin_v=vin_v)
    base_body = _substitute_line(body, td["line_pat"], td["base"], f"{tag} baseline")
    pert_line = td["pert"](amplitude, f_actual)
    pert_body = _substitute_line(body, td["line_pat"], pert_line, f"{tag} perturbed")

    outdir = outdir_root or os.path.join(DATA_ROOT, tag, f"f{freq:.1f}Hz_A{amplitude:g}{label}")

    periodicity = None
    if backend in ("auto", "pss"):
        try:
            diff, inj, N, t_wall, backend_used, periodicity = _measure_h_pss(
                tag, td, f_actual, p, q, amplitude, points_per_cycle, sc_iter,
                steady_coeff, base_body, pert_body, outdir)
        except CouldNotRun as e:
            if backend == "pss":
                raise
            print(f"  [{tag} f={f_actual:.3f}Hz] pss backend COULD-NOT-RUN ({e}); "
                  f"falling back to the long-tran backend", file=sys.stderr)
            diff, inj, N, t_wall, backend_used, periodicity = _measure_h_tran(
                tag, td, f_actual, p, q, amplitude, points_per_cycle,
                base_body, pert_body, outdir, stabtime_mult=tran_stabtime_mult)
    else:
        diff, inj, N, t_wall, backend_used, periodicity = _measure_h_tran(
            tag, td, f_actual, p, q, amplitude, points_per_cycle,
            base_body, pert_body, outdir, stabtime_mult=tran_stabtime_mult)

    Xout = dft_bin(diff, p)
    Xin = dft_bin(inj, p)
    if abs(Xin) < 1e-30:
        raise CouldNotRun(f"{tag}: injected-reference DFT bin is ~0 (p={p}, N={N}) -- bug")
    H = Xout / Xin

    db, deg = db_deg(H)
    return dict(tag=tag, freq_target=freq, p=p, q=q, f_actual=f_actual, relerr=relerr,
                N=N, amplitude=amplitude, H=H, db=db, deg=deg, wall=t_wall,
                backend=backend_used, periodicity=periodicity,
                Xin_check=abs(Xin) / (amplitude / 2.0))


def measure_point(tag, freq, amplitude=0.01, qmax=64, points_per_cycle=200,
                   sc_iter=200, steady_coeff=1e-7, vin_v=None, halve=True, backend="auto",
                   tran_stabtime_mult=10.0):
    """Full disputed-point measurement: H(freq) at `amplitude`, plus a
    halving linearity check (H at amplitude/2, must agree within 0.05 dB
    per the card's own spec -- task text: 'Check linearity by halving the
    amplitude (<=0.05 dB)')."""
    r_full = measure_h(tag, freq, amplitude, qmax=qmax, points_per_cycle=points_per_cycle,
                        sc_iter=sc_iter, steady_coeff=steady_coeff, vin_v=vin_v, label="_A",
                        backend=backend, tran_stabtime_mult=tran_stabtime_mult)
    result = dict(full=r_full)
    if halve:
        r_half = measure_h(tag, freq, amplitude / 2.0, qmax=qmax,
                            points_per_cycle=points_per_cycle, sc_iter=sc_iter,
                            steady_coeff=steady_coeff, vin_v=vin_v, label="_Ahalf",
                            backend=backend, tran_stabtime_mult=tran_stabtime_mult)
        result["half"] = r_half
        result["halving_db_diff"] = r_full["db"] - r_half["db"]
        result["halving_deg_diff"] = r_full["deg"] - r_half["deg"]
    return result


def print_point(result):
    r = result["full"]
    print(f"tag={r['tag']} target={r['freq_target']:.3f}Hz -> actual f=fs*{r['p']}/{r['q']}="
          f"{r['f_actual']:.3f}Hz (relerr {r['relerr']*100:.4f}%), N={r['N']}, "
          f"A={r['amplitude']:g}, backend={r['backend']}, wall={r['wall']:.1f}s, "
          f"|Xin| sanity={r['Xin_check']:.6f} (want ~1)")
    if r["periodicity"]:
        print(f"  periodicity residual (tran backend, relative to that window's own ripple): "
              f"baseline={r['periodicity']['baseline_resid_rel']:.3e}, "
              f"perturbed={r['periodicity']['perturbed_resid_rel']:.3e}")
    print(f"  H = {r['H']:.6g}  -> {r['db']:.4f} dB, {r['deg']:.4f} deg")
    if "half" in result:
        h = result["half"]
        print(f"  halved-amplitude check: A={h['amplitude']:g} (backend={h['backend']}) -> "
              f"{h['db']:.4f} dB, {h['deg']:.4f} deg (diff {result['halving_db_diff']:+.4f} dB, "
              f"{result['halving_deg_diff']:+.4f} deg; target <=0.05 dB)")


# ---------------------------------------------------------------------
# Validation 1: R2 (refs/r2_switched_rc.cir) vs the exact analytic g0().

def validate_r2(freqs, amplitude=0.02, qmax=80):
    sys.path.insert(0, REFS)
    import r2_analytic as ra  # noqa: E402

    fs = 100e3
    body = _load_body(os.path.join(REFS, "r2_switched_rc.cir"))
    line_pat = r"Vin\s+vin\s+0\s+DC\s+1(\s+AC\s+1)?$"

    print("=== validate-r2: pssac_ref.py vs refs/r2_analytic.g0() (exact) ===")
    worst_mag_rel, worst_ph, worst_f = 0.0, 0.0, None
    n_ok = 0
    for freq in freqs:
        p, q, f_actual, relerr = rational_approx(freq, fs, qmax=qmax)
        f_lcm = fs / q
        points = q * 200

        base_body = _substitute_line(body, line_pat, "Vin vin 0 DC 1", "r2 baseline")
        pert_line = f"Vin vin 0 DC 1 SIN(1 {amplitude:.10g} {f_actual:.10g} 0 0 0)"
        pert_body = _substitute_line(body, line_pat, pert_line, "r2 perturbed")

        outdir = os.path.join(DATA_ROOT, "validate_r2", f"f{freq:.1f}Hz")
        base_deck, base_orbit = build_pss_deck(base_body, f_lcm, "100e-6", "vx", points,
                                                "v(vx)", outdir, "baseline", sc_iter=60,
                                                steady_coeff=1e-8)
        pert_deck, pert_orbit = build_pss_deck(pert_body, f_lcm, "100e-6", "vx", points,
                                                "v(vx)", outdir, "perturbed", sc_iter=60,
                                                steady_coeff=1e-8)
        t_wall = 0.0
        for deck, lbl in ((base_deck, "baseline"), (pert_deck, "perturbed")):
            rc, out, wall = run_ngspice(deck, timeout=600)
            t_wall += wall
            if rc != 0:
                print(f"  f={freq}: {lbl} COULD-NOT-RUN exit {rc}: {out[-1000:]}")
                break
        else:
            tb, vb = read_wrdata_pair(base_orbit, 0)
            tp, vp = read_wrdata_pair(pert_orbit, 0)
            N = points
            vb_u = resample_uniform(tb, vb, 1.0 / f_lcm, N)
            vp_u = resample_uniform(tp, vp, 1.0 / f_lcm, N)
            diff = [b - a for a, b in zip(vb_u, vp_u)]
            T = 1.0 / f_lcm
            inj = [amplitude * math.sin(2.0 * math.pi * f_actual * (T * k / N)) for k in range(N)]
            Xout = dft_bin(diff, p)
            Xin = dft_bin(inj, p)
            H = Xout / Xin
            ref = ra.g0(f_actual)
            magrel = abs(abs(H) - abs(ref)) / abs(ref)
            dph = abs(math.degrees(cmath.phase(H)) - math.degrees(cmath.phase(ref)))
            while dph > 180.0:
                dph = abs(dph - 360.0)
            n_ok += 1
            worst_mag_rel = max(worst_mag_rel, magrel)
            if dph > worst_ph:
                worst_ph, worst_f = dph, f_actual
            print(f"  target={freq:.1f}Hz actual={f_actual:.3f}Hz (p/q={p}/{q}, relerr="
                  f"{relerr*100:.3f}%) wall={t_wall:.1f}s: "
                  f"H={H:.6g} ref(g0)={ref:.6g}  |H| relerr={magrel:.3e}  dphase={dph:.4f}deg")
    print(f"\n{n_ok}/{len(freqs)} points ran; worst |H| relative error {worst_mag_rel:.3e} "
          f"(target 1e-3, same as pac_sb0(b)'s own R2 target in tests/test_pac.py), "
          f"worst phase error {worst_ph:.4f} deg at {worst_f} (target 0.05deg's ballpark "
          f"is td-PAC's own target; this independent method is graded here at a looser, "
          f"honestly-stated bound -- see the printed numbers, not restated as pass/FAIL)")
    return worst_mag_rel, worst_ph


# ---------------------------------------------------------------------
# Validation 2: pac_lti.cir vs plain `.ac`.

def validate_lti(freqs, amplitude=0.05, qmax=80):
    fs = 100e3  # Vclk PULSE period 10us
    body = _load_body(os.path.join(FIXTURES, "pac_lti.cir"))
    line_pat = r"Vin\s+in\s+0\s+DC\s+1(\s+AC\s+1)?$"

    ac_outdir = os.path.join(DATA_ROOT, "validate_lti", "ac_ref")
    os.makedirs(ac_outdir, exist_ok=True)

    def ac_at(f_target):
        """Exact single-point `.ac lin 1 f f` at f_target -- no nearest-
        grid-point interpolation error."""
        deck = os.path.join(ac_outdir, f"ac_{f_target:.3f}.cir")
        with open(deck, "w") as fh:
            fh.write("\n".join(_load_body(os.path.join(FIXTURES, "pac_lti.cir"))) + "\n")
            fh.write(f".ac lin 1 {f_target!r} {f_target!r}\n")
            fh.write(".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=17\n"
                      "run\nwrdata ac_out.txt v(out)\nquit\n.endc\n.end\n")
        rc, out, _ = run_ngspice(deck)
        if rc != 0:
            raise CouldNotRun(f".ac reference exit {rc}: {out[-2000:]}")
        f_a, (v_a,) = read_wrdata_complex(os.path.join(ac_outdir, "ac_out.txt"), 1)
        return f_a[0], v_a[0]

    print("=== validate-lti: pssac_ref.py vs plain .ac on pac_lti.cir ===")
    worst_mag_rel, worst_ph, worst_f = 0.0, 0.0, None
    n_ok = 0
    for freq in freqs:
        p, q, f_actual, relerr = rational_approx(freq, fs, qmax=qmax)
        f_lcm = fs / q
        points = q * 200

        base_body = _substitute_line(body, line_pat, "Vin in 0 DC 1", "lti baseline")
        pert_line = f"Vin in 0 DC 1 SIN(1 {amplitude:.10g} {f_actual:.10g} 0 0 0)"
        pert_body = _substitute_line(body, line_pat, pert_line, "lti perturbed")

        outdir = os.path.join(DATA_ROOT, "validate_lti", f"f{freq:.1f}Hz")
        base_deck, base_orbit = build_pss_deck(base_body, f_lcm, "0", "out", points,
                                                "v(out)", outdir, "baseline", sc_iter=60,
                                                steady_coeff=1e-8)
        pert_deck, pert_orbit = build_pss_deck(pert_body, f_lcm, "0", "out", points,
                                                "v(out)", outdir, "perturbed", sc_iter=60,
                                                steady_coeff=1e-8)
        t_wall = 0.0
        ok = True
        for deck, lbl in ((base_deck, "baseline"), (pert_deck, "perturbed")):
            rc, out, wall = run_ngspice(deck, timeout=600)
            t_wall += wall
            if rc != 0:
                print(f"  f={freq}: {lbl} COULD-NOT-RUN exit {rc}: {out[-1000:]}")
                ok = False
                break
        if not ok:
            continue
        tb, vb = read_wrdata_pair(base_orbit, 0)
        tp, vp = read_wrdata_pair(pert_orbit, 0)
        N = points
        vb_u = resample_uniform(tb, vb, 1.0 / f_lcm, N)
        vp_u = resample_uniform(tp, vp, 1.0 / f_lcm, N)
        diff = [b - a for a, b in zip(vb_u, vp_u)]
        T = 1.0 / f_lcm
        inj = [amplitude * math.sin(2.0 * math.pi * f_actual * (T * k / N)) for k in range(N)]
        Xout = dft_bin(diff, p)
        Xin = dft_bin(inj, p)
        H = Xout / Xin
        f_ac_match, v_ac_match = ac_at(f_actual)
        magrel = abs(abs(H) - abs(v_ac_match)) / abs(v_ac_match)
        dph = abs(math.degrees(cmath.phase(H)) - math.degrees(cmath.phase(v_ac_match)))
        while dph > 180.0:
            dph = abs(dph - 360.0)
        n_ok += 1
        worst_mag_rel = max(worst_mag_rel, magrel)
        if dph > worst_ph:
            worst_ph, worst_f = dph, f_actual
        print(f"  target={freq:.1f}Hz actual={f_actual:.3f}Hz (p/q={p}/{q}, relerr="
              f"{relerr*100:.3f}%, nearest .ac point {f_ac_match:.3f}Hz) wall={t_wall:.1f}s: "
              f"H={H:.6g} .ac={v_ac_match:.6g}  |H| relerr={magrel:.3e}  dphase={dph:.4f}deg")
    print(f"\n{n_ok}/{len(freqs)} points ran; worst |H| relative error {worst_mag_rel:.3e}, "
          f"worst phase error {worst_ph:.4f} deg at {worst_f}")
    return worst_mag_rel, worst_ph


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("validate-lti")
    p1.add_argument("--freqs", default="2000,10000,30000,60000,90000")
    p1.add_argument("--amplitude", type=float, default=0.05)
    p1.add_argument("--qmax", type=int, default=80)

    p2 = sub.add_parser("validate-r2")
    p2.add_argument("--freqs", default="1000,10000,49000")
    p2.add_argument("--amplitude", type=float, default=0.02)
    p2.add_argument("--qmax", type=int, default=80)

    p3 = sub.add_parser("point")
    p3.add_argument("--tag", required=True, choices=sorted(TAGS))
    p3.add_argument("--freq", type=float, required=True)
    p3.add_argument("--amplitude", type=float, default=0.01)
    p3.add_argument("--qmax", type=int, default=64)
    p3.add_argument("--points-per-cycle", type=int, default=200)
    p3.add_argument("--sc-iter", type=int, default=200)
    p3.add_argument("--steady-coeff", type=float, default=1e-7)
    p3.add_argument("--vin-v", type=float, default=None)
    p3.add_argument("--no-halve", action="store_true")
    p3.add_argument("--backend", choices=["auto", "pss", "tran"], default="auto")
    p3.add_argument("--tran-stabtime-mult", type=float, default=10.0)

    args = ap.parse_args()
    if args.cmd == "validate-lti":
        freqs = [float(x) for x in args.freqs.split(",")]
        validate_lti(freqs, amplitude=args.amplitude, qmax=args.qmax)
    elif args.cmd == "validate-r2":
        freqs = [float(x) for x in args.freqs.split(",")]
        validate_r2(freqs, amplitude=args.amplitude, qmax=args.qmax)
    elif args.cmd == "point":
        result = measure_point(args.tag, args.freq, amplitude=args.amplitude, qmax=args.qmax,
                                points_per_cycle=args.points_per_cycle, sc_iter=args.sc_iter,
                                steady_coeff=args.steady_coeff, vin_v=args.vin_v,
                                halve=not args.no_halve, backend=args.backend,
                                tran_stabtime_mult=args.tran_stabtime_mult)
        print_point(result)


if __name__ == "__main__":
    try:
        main()
    except CouldNotRun as e:
        print(f"COULD-NOT-RUN: {e}", file=sys.stderr)
        sys.exit(2)
