#!/usr/bin/env python3
"""Validation for periodic (cyclostationary) noise on the shooting-PSS orbit
(`.pnoise ... method="shoot" solver="td" [dev=1]`), card #509, private-fork
branch aug-pnoise, engine in src/spicelib/analysis/pacshoot.c
(pac_td_pnoise()/pac_engine_run_adjoint_noise()), built on the true-adjoint
`.pxf ... adjoint=1` engine (card #504, pac_engine_run_adjoint()) -- see that
function's own header for the shared Phi/G(w)/M(w)/backward-sweep machinery,
copied verbatim (not refactored/shared) into pac_engine_run_adjoint_noise()
so this card touches zero lines of the already-validated .pxf adjoint code.

Read docs/pss-pac.md's `.pnoise` section (added this card) before changing
anything here -- it has the full derivation (including the two WRONG
per-step quadrature-weight candidates this card tried and rejected before
finding the right one: "no weight at all" and "multiply by the output tap's
own w_k/T", both found wrong by check_pnoise_lti() below during development),
the syntax, and the validation numbers this file's checks reproduce.

Design in one line: onoise(f) = T * Sum_k (1/w_k) * dens_k(f), dens_k(f) =
Sum_g S_g(t_k,f)*|nu_k(f).b_g|^2 read straight off the ordinary
NevalSrc()/DEVnoise() machinery by loading ckt->CKTrhs/CKTirhs with nu_k(f)
(the backward adjoint sweep's own per-step SMPtSolve() result, BEFORE the
phase multiply pac_engine_run_adjoint()'s Lambda(w) needs) at the step-k
bias, w_k = 0.5*(h_k+h_next) the SAME periodic-trapezoidal weight the
boundary term's "coef" uses.

Scope (v1, card #509): resistor thermal, diode shot, MOS1/BSIM4 channel
thermal+flicker noise only -- any OTHER noisy device type in the circuit is
REFUSED loudly (pac_noise_check_devices(), matched by DEVnoise function
pointer). Saltation (ideal-switch) events and free-running orbits are
REFUSED too (same PAC_ADJOINT_REFUSED scope cut as `.pxf adjoint=1`), with
NO fallback engine (unlike .pxf's adjoint=0 -- noise has no "one pass per
named source" analogue to fall back to). Input-referred (inoise) is not
computed by this engine.

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_pnoise.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# Reuse test_pac.py's own already-proven helpers (NGSPICE resolution,
# ngspice invocation/timeout handling, the fixture-body reader) -- same
# reason test_pxf.py does this.
import test_pac as _pac  # noqa: E402

run_ngspice = _pac.run_ngspice
body_lines = _pac.body_lines
CouldNotRun = _pac.CouldNotRun

DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/pnoise")

results = []


def record(name, state, detail):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN", "info")
    results.append((name, state, detail))


def read_wrdata_real(path, n):
    """n real columns from a wr_singlescale/wr_vecnames file -- same format
    test_sens.py's own read_wrdata_real() reads, duplicated here (small,
    self-contained) rather than importing that file (avoids pulling in its
    own, unrelated module-level setup)."""
    if not os.path.exists(path):
        raise CouldNotRun(f"ngspice wrote no {path}")
    freqs, cols = [], [[] for _ in range(n)]
    with open(path) as fh:
        next(fh)
        for row in fh:
            v = row.split()
            if not v:
                continue
            v = [float(x) for x in v]
            if len(v) != 1 + n:
                raise CouldNotRun(f"unexpected column count in {path}: {row!r}")
            freqs.append(v[0])
            for k in range(n):
                cols[k].append(v[1 + k])
    if not freqs:
        raise CouldNotRun(f"{path} is empty")
    return freqs, cols


def check_pnoise_lti():
    """LTI circuit (pxf_lti.cir: Rin, Rload, Cout -- resistors and a cap,
    driven by an UNRELATED periodic clock so PSS has a genuine period, same
    fixture .pxf's own card #499/#504 LTI checks use) -- `.pnoise ...
    solver=td` onoise_spectrum must equal plain `.noise`'s onoise_spectrum
    within the card's own 1e-6 relative target (done-when item 1). Checked
    at THREE grid resolutions (Points=200/2000/20000, the SAME PSS points=
    this engine's own recording passes reuse -- no pacpoints= override
    needed, this fixture is cheap) to show real, clean (better than the
    nominal O(h) this engine family's OTHER checks -- test_pac.py's own
    accuracy-floor note -- typically show) convergence, not a lucky single
    point: measured this pass, 200->1.7e-2 max / 8e-3 typical, 2000->8e-5,
    20000->4e-7 -- i.e. the 1e-6 target is reached for real at 20000, not
    silently loosened."""
    name = "pnoise_lti: pxf_lti.cir, .pnoise solver=td onoise_spectrum vs plain .noise (1e-6 target)"
    try:
        outdir = os.path.join(DATA_ROOT, "lti")
        os.makedirs(outdir, exist_ok=True)
        body = "\n".join(body_lines("pxf_lti.cir"))

        # .noise reference -- needs an AC-flagged input source (Vin here;
        # onoise itself does not depend on WHICH source is "input", only
        # inoise would).
        noise_body = body.replace("Vin in 0 DC 1", "Vin in 0 DC 1 AC 1", 1)
        if noise_body == body:
            raise CouldNotRun("Vin DC 1 line not found for AC substitution")
        noise_deck = os.path.join(outdir, "noise_ref.cir")
        with open(noise_deck, "w") as fh:
            fh.write(noise_body + "\n")
            fh.write(".noise v(out) Vin dec 10 1k 1meg\n")
            fh.write(".control\nrun\nsetplot noise1\n"
                      "set wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                      "wrdata noise_ref.txt onoise_spectrum\nquit\n.endc\n.end\n")
        rc, out, _ = run_ngspice(noise_deck, timeout=120)
        if rc != 0:
            raise CouldNotRun(f".noise reference exit {rc}: {out[-2000:]}")
        f_ref, (v_ref,) = read_wrdata_real(os.path.join(outdir, "noise_ref.txt"), 1)

        worst_by_points = {}
        for points in (200, 2000, 20000):
            deck = os.path.join(outdir, f"pnoise_{points}.cir")
            with open(deck, "w") as fh:
                fh.write(body + "\n")
                fh.write(f'.pnoise 100k 0 clkgnd {points} 1 5 0.001 out Vin dec 10 1k 1meg '
                          f'method="shoot" solver="td"\n')
                fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                          "wrdata pnoise_out.txt onoise_spectrum\nquit\n.endc\n.end\n")
            rc, out, wall = run_ngspice(deck, timeout=300)
            if rc != 0:
                raise CouldNotRun(f".pnoise (Points={points}) exit {rc}: {out[-2000:]}")
            f_pn, (v_pn,) = read_wrdata_real(os.path.join(outdir, "pnoise_out.txt"), 1)
            if len(f_pn) != len(f_ref):
                raise CouldNotRun(f"Points={points}: point count mismatch "
                                   f".pnoise {len(f_pn)} vs .noise {len(f_ref)}")
            worst, worst_f = 0.0, None
            for fa, va, fb, vb in zip(f_ref, v_ref, f_pn, v_pn):
                if abs(fa - fb) > 1e-6 * max(fa, fb):
                    raise CouldNotRun(f"frequency mismatch {fa} vs {fb}")
                rel = abs(va - vb) / va if va > 0 else abs(vb)
                if rel > worst:
                    worst, worst_f = rel, fa
            worst_by_points[points] = (worst, worst_f, wall)

        detail = ("; ".join(f"Points={p}: worst rel {w:.3e} at {wf:.4g} Hz ({wall:.2f}s)"
                             for p, (w, wf, wall) in worst_by_points.items()))
        final_worst = worst_by_points[20000][0]
        # Monotonic-ish convergence check: 20000 must be the best of the
        # three (a real accuracy floor, not noise in the measurement).
        converging = worst_by_points[200][0] > worst_by_points[2000][0] > worst_by_points[20000][0]
        if final_worst <= 1e-6 and converging:
            record(name, "pass", detail)
        else:
            record(name, "FAIL", detail + f" (target 1e-6 at Points=20000, converging={converging})")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def _sum_sq_consistency(total, parts_by_name):
    """max relative error between `total` (per-freq list) and
    sqrt(sum of squares of parts_by_name's own per-freq lists) -- the
    internal self-consistency check every .pnoise dev=1 table must satisfy
    exactly (both are read off the SAME `dens` accumulation inside ONE
    engine call, see pac_engine_run_adjoint_noise()'s header: the total is
    NevalSrc's own per-instance sum, the table is the SAME generators'
    individual densities -- P^2 = sum(p_i^2) by construction, no folding
    or cross-terms involved)."""
    worst = 0.0
    for i, t in enumerate(total):
        s = math.sqrt(sum(p[i] ** 2 for p in parts_by_name.values()))
        worst = max(worst, abs(s - t) / t if t > 0 else abs(s))
    return worst


def check_pnoise_r3_buck():
    """R3 voltage-mode buck (pac_r3_ctrl2out.cir, smooth tanh switches,
    project ruling #447 A, the SAME fixture card #504 found its own M(w)-
    convention bug on -- has a capacitor, C1f, whose negNode is NOT ground)
    -- the card's own required "buck output-noise example... with the
    per-device table" fixture (done-when item 3). Every noisy device here
    is a plain resistor (Rdcr/Resr/Rload/Ra/Rb/R1in/R2f/R3in/Roa -- the
    op-amp/switches/comparator are all noiseless behavioural B/E sources),
    squarely in v1 scope. Checks: runs to completion (no REFUSED), every
    onoise value finite and >=0, and the per-device table's sum-of-squares
    equals the reported onoise_spectrum total to near machine precision
    (see _sum_sq_consistency()) -- a real, non-trivial internal-consistency
    check (NOT a tautology: onoise_spectrum and the per-device columns are
    accumulated through TWO DIFFERENT DEVnoise() call sites --
    prtSummary=0 vs prtSummary=1 -- inside the SAME per-step loop, so this
    catches a real class of bug, e.g. a table column silently double- or
    under-counted relative to the scalar total). No independent analytic
    or transient-noise-Monte-Carlo cross-check was run on this fixture
    (see docs/pss-pac.md's `.pnoise` section, "what was skipped") --
    pxf_lti.cir (check_pnoise_lti() above) is this card's own quantitative
    validation against a trusted reference."""
    name = "pnoise_r3: pac_r3_ctrl2out.cir, .pnoise solver=td dev=1 (self-consistency + sanity)"
    try:
        outdir = os.path.join(DATA_ROOT, "r3")
        os.makedirs(outdir, exist_ok=True)
        body = "\n".join(body_lines("pac_r3_ctrl2out.cir"))
        gens = ["Rdcr", "Resr", "Rload", "Ra", "Rb", "R1in", "R2f", "R3in", "Roa"]
        deck = os.path.join(outdir, "pnoise.cir")
        with open(deck, "w") as fh:
            fh.write(body + "\n")
            fh.write('.pnoise 500e3 300e-6 vout 200 10 150 1e-6 vout Vin dec 5 500 50000 '
                      'method="shoot" solver="td" dev=1\n')
            fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                      "wrdata pnoise_out.txt onoise_spectrum " +
                      " ".join(f"onoise_{g}" for g in gens) + "\nquit\n.endc\n.end\n")
        rc, out, wall = run_ngspice(deck, timeout=300)
        if rc != 0:
            raise CouldNotRun(f".pnoise exit {rc}: {out[-2000:]}")
        if "REFUSED" in out:
            raise CouldNotRun(f".pnoise REFUSED unexpectedly: {out[-2000:]}")
        f, cols = read_wrdata_real(os.path.join(outdir, "pnoise_out.txt"), 1 + len(gens))
        total, parts = cols[0], dict(zip(gens, cols[1:]))
        if not all(v >= 0.0 and v == v for v in total):  # v==v: not NaN
            raise CouldNotRun(f"non-finite or negative onoise_spectrum value(s): {total}")
        worst = _sum_sq_consistency(total, parts)
        detail = (f"{len(f)} points 500Hz-50kHz, worst onoise range "
                   f"[{min(total):.3e}, {max(total):.3e}] V/rtHz, "
                   f"per-device sum-of-squares worst rel {worst:.3e}, wall {wall:.2f}s")
        record(name, "pass" if worst <= 1e-9 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pnoise_r2_smooth():
    """Switched-RC LPTV fixture (pnoise_r2_smooth.cir -- SAME core topology
    as pxf_r2_switched_rc.cir but a SMOOTH tanh switch in place of the
    ideal `sw`, since this engine REFUSES saltation with no fallback, see
    that fixture's own header) -- the card's own required "switched-RC"
    fixture (done-when item 2). Same checks as check_pnoise_r3_buck():
    completes, finite/non-negative, per-device table self-consistent.
    This check itself is self-consistency only; the INDEPENDENT checks are
    check_pnoise_track_hold_analytic() and the optional check_pnoise_montecarlo().
    (Original first-pass note: not an independent cross-check against an analytic folded-kT/C result
    or a transient-noise Monte Carlo -- both were judged too expensive/
    error-prone to derive correctly in the time available this pass (a
    correctly-normalised Monte Carlo needs a PWL-modulated trnoise
    injection matched to this engine's own w_k convention, many transient
    runs, and a periodogram; an analytic folded result for a GENUINELY
    LPTV transfer function, as opposed to pxf_lti.cir's constant-density-
    but-LTI-transfer case, needs re-deriving the folding integral by hand)
    -- an honest, stated gap, see docs/pss-pac.md's `.pnoise` section."""
    name = "pnoise_r2_smooth: pnoise_r2_smooth.cir, .pnoise solver=td dev=1 (self-consistency + sanity; independent MC check is the optional PNOISE_MC=1 test)"
    try:
        outdir = os.path.join(DATA_ROOT, "r2_smooth")
        os.makedirs(outdir, exist_ok=True)
        body = "\n".join(body_lines("pnoise_r2_smooth.cir"))
        gens = ["Rin", "Rsw"]
        deck = os.path.join(outdir, "pnoise.cir")
        with open(deck, "w") as fh:
            fh.write(body + "\n")
            fh.write('.pnoise 100e3 100e-6 vx 2000 1 50 1e-6 vx Vin dec 5 1k 1meg '
                      'method="shoot" solver="td" dev=1\n')
            fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                      "wrdata pnoise_out.txt onoise_spectrum " +
                      " ".join(f"onoise_{g}" for g in gens) + "\nquit\n.endc\n.end\n")
        rc, out, wall = run_ngspice(deck, timeout=300)
        if rc != 0:
            raise CouldNotRun(f".pnoise exit {rc}: {out[-2000:]}")
        if "REFUSED" in out:
            raise CouldNotRun(f".pnoise REFUSED unexpectedly: {out[-2000:]}")
        f, cols = read_wrdata_real(os.path.join(outdir, "pnoise_out.txt"), 1 + len(gens))
        total, parts = cols[0], dict(zip(gens, cols[1:]))
        if not all(v >= 0.0 and v == v for v in total):
            raise CouldNotRun(f"non-finite or negative onoise_spectrum value(s): {total}")
        worst = _sum_sq_consistency(total, parts)
        detail = (f"{len(f)} points 1kHz-1MHz, onoise range "
                   f"[{min(total):.3e}, {max(total):.3e}] V/rtHz, "
                   f"per-device sum-of-squares worst rel {worst:.3e}, wall {wall:.2f}s")
        record(name, "pass" if worst <= 1e-9 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pnoise_refuses_unsupported_device():
    """pxf_r2_switched_rc.cir's ideal `sw` (Switch model, ngspice's own
    swnoise.c DOES register a DEVnoise for it) is OUTSIDE this engine's v1
    device scope (R/diode/MOS1/BSIM4 only) -- `.pnoise ... solver=td` must
    REFUSE loudly (pac_noise_check_devices(), pacshoot.c) rather than
    silently running with the Switch's own noise contribution dropped or
    mis-scaled. This fixture's device-scope refusal fires BEFORE this
    engine would even reach its saltation check (PSSsaltEventsSeen(), the
    OTHER refusal this engine shares with `.pxf adjoint=1` -- see
    pac_engine_run_adjoint_noise()'s header) -- that second refusal path is
    a near-verbatim copy of the already-validated .pxf adjoint code and was
    not separately exercised by a fixture this pass (see docs/pss-pac.md's
    `.pnoise` section, "what was skipped"). The done-when here is: ngspice
    itself does not crash/hang (exit 0, same as a `.pac`/`.pstb` post-
    processing failure on an otherwise-converged PSS orbit -- see dcpss.c's
    existing "WARNING -- ... time-domain post-processing failed" idiom,
    reused verbatim by pssshoot.c for this engine too) and prints REFUSED
    naming the actual unsupported device type."""
    name = "pnoise_refuse: pxf_r2_switched_rc.cir (ideal sw), .pnoise solver=td REFUSES (device scope)"
    try:
        outdir = os.path.join(DATA_ROOT, "refuse_device")
        os.makedirs(outdir, exist_ok=True)
        body = "\n".join(body_lines("pxf_r2_switched_rc.cir"))
        deck = os.path.join(outdir, "pnoise.cir")
        with open(deck, "w") as fh:
            fh.write(body + "\n")
            fh.write('.pnoise 100e3 100e-6 vx 2000 1 25 1e-6 vx Vin dec 5 1k 40k '
                      'method="shoot" solver="td"\n')
            fh.write(".control\nrun\nquit\n.endc\n.end\n")
        rc, out, wall = run_ngspice(deck, timeout=120)
        refused = "PNOISE-TD: REFUSED" in out
        named_switch = "'Switch'" in out or '"Switch"' in out
        detail = f"exit {rc}, REFUSED seen={refused}, names 'Switch'={named_switch}, wall {wall:.2f}s"
        if rc == 0 and refused and named_switch:
            record(name, "pass", detail)
        else:
            record(name, "FAIL", detail + f": {out[-1500:]}")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def th_analytic_S0(rtrack, c, T, tt):
    """Exact f->0 one-sided PSD [V^2/Hz] of a track-and-hold's held output
    (ideal-switch limit; Rtrack the only noise source). x(t) is an
    Ornstein-Uhlenbeck process (variance kT/C, tau = Rtrack*C) RESTARTED
    from the held value at each track start, run for tt, then frozen for
    Th = T-tt. Per-cycle integral Z_n = A_n + Th*s_n (A_n = integral of the
    OU piece, s_n = its end value). Long-run variance Var(Z)+2Cov(Z_n,Z_n+1)
    +2Cov(Z_n,Z_n+2), S(0) = 2*LRV/T. (The naive 'ZOH of a kT/C sample'
    answer, 2kT/C*T, is 1.4 dB too high for a 20% track window: the sample
    is correlated with the track-phase samples that precede it.)"""
    kT = 1.380649e-23 * 300.15
    tau = rtrack * c
    a = tt / tau
    Th = T - tt
    q = 1 - math.exp(-a)
    rho = math.exp(-a)
    varA = 2 * tau * (tt - tau * q)
    covAs = tau * q
    var0 = varA + Th * Th + 2 * Th * covAs
    c1 = (tau * q) * covAs + Th * rho * covAs + Th * tau * q + Th * Th * rho
    c2 = ((tau * q) * covAs * rho + Th * rho * rho * covAs
          + Th * tau * q * rho + Th * Th * rho * rho)
    return 2 * (kT / c) * (var0 + 2 * c1 + 2 * c2) / T


def check_pnoise_track_hold_analytic():
    """Independent ANALYTIC cyclostationary check (card #509): a switched-RC
    track-and-hold (pnoise_track_hold_slow.cir: Rtrack=25k, C=10p,
    tau=250 ns, 20 % track window = 8 tau, fclk=100 kHz). Rtrack's noise is
    folded from ~10 harmonics of fs into the baseband, so this exercises the
    modulation + folding + per-step weighting that the LTI check cannot.
    `.pnoise solver=td` at 300 Hz-3 kHz must equal th_analytic_S0() *
    sinc^2(fT) (the sinc shape correction there is < 0.01 dB) within 0.1 dB
    (measured -0.013..-0.017 dB at Points=2000, the residual being the
    finite gate conductance/tanh edges and grid). The naive kT/C ZOH value
    is 1.38 dB away, so this discriminates."""
    name = "pnoise_track_hold_analytic: pnoise_track_hold_slow.cir vs exact folded kT/C (OU restart) at 0.3-3 kHz"
    try:
        outdir = os.path.join(DATA_ROOT, "th_analytic")
        os.makedirs(outdir, exist_ok=True)
        body = "\n".join(body_lines("pnoise_track_hold_slow.cir"))
        T = 1e-5
        S0 = th_analytic_S0(25e3, 10e-12, T, 0.2 * T)
        worst = 0.0
        rows = []
        wall_tot = 0.0
        for f in (300.0, 1000.0, 3000.0):
            deck = os.path.join(outdir, f"pn_{f:g}.cir")
            with open(deck, "w") as fh:
                fh.write(body + "\n")
                fh.write(f'.pnoise 100e3 100e-6 vout 2000 1 50 1e-6 vout Vin lin 1 {f} {f} '
                          'method="shoot" solver="td"\n')
                fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                          "wrdata pn_out.txt onoise_spectrum\nquit\n.endc\n.end\n")
            rc, out, wall = run_ngspice(deck, timeout=300)
            wall_tot += wall
            if rc != 0 or "REFUSED" in out:
                raise CouldNotRun(f".pnoise exit {rc}: {out[-1500:]}")
            _, cols = read_wrdata_real(os.path.join(outdir, "pn_out.txt"), 1)
            x = math.pi * f * T
            an = math.sqrt(S0 * (math.sin(x) / x) ** 2)
            d = 20 * math.log10(cols[0][0] / an)
            worst = max(worst, abs(d))
            rows.append(f"{f:g}Hz {cols[0][0]:.4e} vs {an:.4e} ({d:+.3f} dB)")
        detail = "; ".join(rows) + f"; worst {worst:.3f} dB (limit 0.1), wall {wall_tot:.1f}s"
        record(name, "pass" if worst <= 0.1 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pnoise_track_hold_mos(env=None, red=False):
    """card #544 item 1(b): the same exact folded kT/C closed form (th_analytic_S0, OU restart) with the ideal/tanh switch replaced
    by a MOSFET whose gate is clocked: pnoise_track_hold_mos1.cir (MOS1, Meyer + overlap caps) and pnoise_track_hold_bsim4.cir
    (BSIM4 capmod 2, charge model).  Both go through the generic model-charge path (3 state rows).  The switch's own channel
    thermal noise in triode (0.9 % / 0.24 % of the series resistance, R = Rtrack + Ron) is part of the answer, so the closed
    form uses R = Rtrack + Ron (Ron estimated by hand from kp, W/L and the body-effect threshold, 455 / 120 ohm; the closed form moves by
    < 0.01 dB over that range).  Fixture sizing: Cch/C1 = 0.034 % / 0.09 % (channel charge dumped at turn-off is the one effect
    the textbook formula does not contain: with W=400u, 23 % of C1, the same deck is 0.8 dB off, a physical sampling effect and
    not an engine error), fclk 10 kHz, tau = 2.5 us, track window 20 us = 8 tau.  Bound 0.1 dB at 100 Hz and 1 kHz (measured
    -0.039/-0.000 dB MOS1, +0.006/+0.045 dB BSIM4; the grid residual of the same engine on the tanh switch is -0.016 dB; the naive ZOH
    kT/C value is 1.38 dB away).  red=True: PAC_NATIVE_HISTORY=1 must give > 0.5 dB error (measured -0.97 dB MOS1, -3.45 dB BSIM4)."""
    tag = " [RED PROOF: PAC_NATIVE_HISTORY=1, must exceed 0.5 dB]" if red else ""
    name = "pnoise_track_hold_mos: MOS1 / BSIM4 switch, folded kT/C closed form at 100 Hz and 1 kHz" + tag
    try:
        outdir = os.path.join(DATA_ROOT, "th_mos")
        os.makedirs(outdir, exist_ok=True)
        T = 1e-4
        rows, worst, wall_tot = [], 0.0, 0.0
        for fx, ron in (("pnoise_track_hold_mos1.cir", 455.0), ("pnoise_track_hold_bsim4.cir", 120.0)):
            body = "\n".join(body_lines(fx))
            S0 = th_analytic_S0(50e3 + ron, 50e-12, T, 0.2 * T)
            for f in ((1000.0,) if red else (100.0, 1000.0)):
                deck = os.path.join(outdir, f"pn_{fx[:-4]}_{f:g}.cir")
                with open(deck, "w") as fh:
                    fh.write(body + "\n.options reltol=1e-4\n")
                    fh.write(f'.pnoise 10e3 1e-3 vout 2000 1 50 1e-6 vout Vin lin 1 {f} {f} method="shoot" solver="td"\n')
                    fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                             "wrdata pn_out.txt onoise_spectrum\nquit\n.endc\n.end\n")
                e = dict(os.environ, OMP_NUM_THREADS="2")
                e.update(env or {})
                rc, out, wall = run_ngspice(deck, timeout=600, env=e)
                wall_tot += wall
                if rc != 0 or "REFUSED" in out:
                    raise CouldNotRun(f".pnoise exit {rc}: {out[-1500:]}")
                _, cols = read_wrdata_real(os.path.join(outdir, "pn_out.txt"), 1)
                x = math.pi * f * T
                an = math.sqrt(S0) * (math.sin(x) / x)
                d = 20 * math.log10(cols[0][0] / an)
                worst = max(worst, abs(d))
                rows.append(f"{fx[:-4].split("_")[-1]} {f:g}Hz {cols[0][0]:.4e} vs {an:.4e} ({d:+.3f} dB)")
        detail = "; ".join(rows) + f"; worst {worst:.3f} dB, wall {wall_tot:.1f}s"
        ok = (worst > 0.5) if red else (worst <= 0.1)
        record(name, "pass" if ok else "FAIL", detail + (" (limit 0.1)" if not red else ""))
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


# ---------------------------------------------------------------------------------------------------------------------------
# card #556: flicker (1/f) noise folded by a fast periodic modulation.  A source n(t) = g(t) w(t), w stationary with S_w(f) = A/f^ef,
# gives  S_out(f) = Sum_q |G_q(f)|^2 S_w(|f - q fs|),  G_q = Fourier coefficients of g(t) H(t,f)  -- the density is NOT S_w(f) for
# every sideband.  Closed forms below are exact for the three fixtures (no engine output enters).

K_B, T_K = 1.380649e-23, 300.15


def ref_trap_coeffs(T, tr, pw, lo=-1.0, hi=1.0, tf=None, td=0.0, Q=4000):
    """Mean and Fourier coefficients c_q (q = 1..Q, f(t) = Sum c_q e^{+j q w t}) of PULSE(lo hi td tr tf pw T), exact: the second
    derivative of a piecewise-linear function is a sum of Dirac deltas at its corners."""
    import cmath
    tf = tr if tf is None else tf
    ts = [td, td + tr, td + tr + pw, td + tr + pw + tf]
    sl = [(hi - lo) / tr, 0.0, -(hi - lo) / tf, 0.0]
    slp = [0.0, (hi - lo) / tr, 0.0, -(hi - lo) / tf]
    ds = [a - b for a, b in zip(sl, slp)]
    cq = []
    for q in range(1, Q + 1):
        th = 2 * math.pi * q / T
        cq.append(-sum(d * cmath.exp(-1j * th * t) for t, d in zip(ts, ds)) / (th * th * T))
    mean = (lo * (T - tr - pw - tf) + hi * pw + 0.5 * (lo + hi) * (tr + tf)) / T
    return mean, cq


def ref_S_out(f, fs, mean, cq, S):
    """Sum_q |c_q|^2 S(|f - q fs|) over q = -Q..Q (c_-q = conj c_q)."""
    tot = mean ** 2 * S(abs(f))
    for i, c in enumerate(cq):
        q = i + 1
        tot += abs(c) ** 2 * (S(abs(f - q * fs)) + S(abs(f + q * fs)))
    return tot


def ref_chop_lp(f, kind):
    """pnoise_fold_chop_lp.cir: v_x = Z1 i_n; x sq1 -> LP2 -> x sq2.  Gain from the input at f' = f - m fs to the output at f:
    K_m = Z1(f') Sum_p c1_p H2(f' + p fs) c2_{m-p}; S_out = Sum_m |K_m|^2 S_i(|f'|)."""
    fs = 1e5
    T = 1 / fs
    tr = 1e-6
    Q = 400
    m1, c1p = ref_trap_coeffs(T, tr, 0.5 * T - tr, Q=Q)
    m2, c2p = ref_trap_coeffs(T, tr, 0.3 * T - tr, lo=0.0, hi=1.0, td=0.2 * T, Q=Q)

    def coef(mean, cp, p):
        return mean if p == 0 else (cp[p - 1] if p > 0 else cp[-p - 1].conjugate())
    R1, C1, R2C2 = 1e3, 1.5915e-9, 1e3 * 3.183e-9
    Z1 = lambda x: R1 / (1 + 2j * math.pi * x * R1 * C1)
    H2 = lambda x: 1 / (1 + 2j * math.pi * x * R2C2)
    tot = 0.0
    for m in range(-Q + 1, Q):
        fp = f - m * fs
        Si = 4 * K_B * T_K / R1 if kind == "th" else 4e-12 * 1e-6 / abs(fp)
        K = 0
        for p in range(-Q + 1, Q):
            r = m - p
            if abs(r) <= Q:
                K += coef(m1, c1p, p) * H2(fp + p * fs) * coef(m2, c2p, r)
        tot += abs(K * Z1(fp)) ** 2 * Si
    return tot


def _pn_values(fixture, f, vecs, outdir, node, env=None, extra=""):
    """.pnoise (solver=td, dev=1) of a fold fixture at one frequency; returns {vector: value}."""
    os.makedirs(outdir, exist_ok=True)
    body = "\n".join(body_lines(fixture))
    deck = os.path.join(outdir, f"{fixture[:-4]}_{f:g}.cir")
    with open(deck, "w") as fh:
        fh.write(body + "\n")
        fh.write(f'.pnoise 100k 0 {node[0]} 4000 1 50 1e-6 {node[1]} {node[2]} lin 1 {f:g} {f:g} method="shoot" solver="td" dev=1\n')
        fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                 f"wrdata pn.txt {' '.join(vecs)}\nquit\n.endc\n.end\n")
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    rc, out, wall = run_ngspice(deck, timeout=300, env=e)
    if rc != 0 or "REFUSED" in out or "FAILED" in out:
        raise CouldNotRun(f"{fixture} {f:g} Hz: exit {rc}: {out[-1500:]}")
    _, cols = read_wrdata_real(os.path.join(outdir, "pn.txt"), len(vecs))
    return {v: c[0] for v, c in zip(vecs, cols)}


def check_pnoise_fold(red=False):
    """card #556: folded flicker vs exact closed forms, three fixtures.  (a) pnoise_fold_chop_res: 1/f of a resistor chopped by an ideal
    +-1 modulator (output at 1k/30k/90k/150k Hz, i.e. near DC and near and above fs = 100 kHz); (b) pnoise_fold_chop_lp: LP - chop - LP -
    chop with complex sideband gains (tests the sign of the f -/+ q fs pairing: the wrong sign is 10-13 % off); (c) pnoise_fold_mos:
    MOS1 whose flicker density is modulated by a fast periodic gate drive (g(t) = gm(t)).  Thermal / channel-thermal contributions
    (no folding to evaluate) are checked in the same runs.  Bound 1e-5 relative (measured <= 7e-7).  red=True: PNOISE_FLICKER_OLD=1
    (the pre-#556 evaluation of the density at the output frequency) must be off by > 2 % (measured: +10x chopper at 1 kHz, +50 %
    LP-chop at 20 kHz, -13 % at 70 kHz, +3.4 % / -9 % MOS)."""
    tag = " [RED PROOF: PNOISE_FLICKER_OLD=1, must be off by > 2 %]" if red else ""
    name = "pnoise_fold: flicker folded by a fast periodic modulation vs closed forms (chopper, LP-chop-LP-chop, MOS1 gate-driven)" + tag
    try:
        outdir = os.path.join(DATA_ROOT, "fold")
        env = {"PNOISE_FLICKER_OLD": "1"} if red else None
        fs = 1e5
        T = 1 / fs
        tr = 1e-6
        rows, worst, worst_th, worst_cons = [], 0.0, 0.0, 0.0
        errs_old = []

        # (a) resistor 1/f + thermal, chopped
        mean, cq = ref_trap_coeffs(T, tr, 0.5 * T - tr)
        R, I, kf = 1e3, 1e-3, 4e-12
        th_s = lambda x: R * 4 * K_B * T_K
        fl_s = lambda x: R * R * kf * I ** 2 / x
        for f in (1e3, 30e3, 90e3, 150e3):
            v = _pn_values("pnoise_fold_chop_res.cir", f, ["onoise_r1_thermal", "onoise_r1_1overf", "onoise_r1", "onoise_spectrum"],
                           outdir, ("sq", "out", "Vsq"), env)
            # the corrected flicker slot, the instance total and the scalar spectrum must stay one consistent set
            cons = max(abs(math.hypot(v["onoise_r1_thermal"], v["onoise_r1_1overf"]) / v["onoise_r1"] - 1),
                       abs(v["onoise_spectrum"] / v["onoise_r1"] - 1))
            worst_cons = max(worst_cons, cons)
            rf = math.sqrt(ref_S_out(f, fs, mean, cq, fl_s))
            rt = math.sqrt(ref_S_out(f, fs, mean, cq, th_s))
            e1, e2 = v["onoise_r1_1overf"] / rf - 1, v["onoise_r1_thermal"] / rt - 1
            worst, worst_th = max(worst, abs(e1)), max(worst_th, abs(e2))
            errs_old.append(abs(e1))
            rows.append(f"chop {f / 1e3:g}k {v['onoise_r1_1overf']:.4e} vs {rf:.4e} ({e1:+.1e})")
        # (b) LP - chop - LP - chop
        for f in (20e3, 70e3):
            v = _pn_values("pnoise_fold_chop_lp.cir", f, ["onoise_r1_thermal", "onoise_r1_1overf"], outdir, ("sq1", "out", "Vsq1"), env)
            rf, rt = math.sqrt(ref_chop_lp(f, "fl")), math.sqrt(ref_chop_lp(f, "th"))
            e1, e2 = v["onoise_r1_1overf"] / rf - 1, v["onoise_r1_thermal"] / rt - 1
            worst, worst_th = max(worst, abs(e1)), max(worst_th, abs(e2))
            errs_old.append(abs(e1))
            rows.append(f"lp-chop {f / 1e3:g}k {v['onoise_r1_1overf']:.4e} vs {rf:.4e} ({e1:+.1e})")
        # (c) MOS1 with a modulated bias: g(t) = gm(t) = beta vgst(t), vgst 0.5 <-> 0.9 V
        eps_ox, Cox = 3.9 * 8.854214871e-12, 3.9 * 8.854214871e-12 / 2e-8
        L, W, kp, kfm, Rd = 1e-6, 10e-6, 2e-4, 1e-24, 2e3
        beta = kp * W / L
        mean, cq = ref_trap_coeffs(T, tr, 0.5 * T - tr, lo=0.5, hi=0.9)
        fl_m = lambda x: Rd * Rd * kfm * beta * beta / (W * L * Cox) / x
        for f in (2e3, 90e3):
            v = _pn_values("pnoise_fold_mos.cir", f, ["onoise_m1_id", "onoise_m1_1overf"], outdir, ("g", "d", "Vg"), env)
            rf = math.sqrt(ref_S_out(f, fs, mean, cq, fl_m))
            rt = math.sqrt(Rd * Rd * 4 * K_B * T_K * (2.0 / 3.0) * beta * mean)
            e1, e2 = v["onoise_m1_1overf"] / rf - 1, v["onoise_m1_id"] / rt - 1
            worst, worst_th = max(worst, abs(e1)), max(worst_th, abs(e2))
            errs_old.append(abs(e1))
            rows.append(f"mos {f / 1e3:g}k {v['onoise_m1_1overf']:.4e} vs {rf:.4e} ({e1:+.1e})")
        detail = "; ".join(rows) + f"; worst flicker {worst:.2e}, worst thermal {worst_th:.2e}, slot/total/spectrum consistency {worst_cons:.1e}"
        if red:
            ok = min(errs_old) > 0.02
            record(name, "pass" if ok else "FAIL", detail + " (smallest old-evaluation error must exceed 0.02)")
        else:
            ok = worst <= 1e-5 and worst_th <= 1e-5 and worst_cons <= 1e-9
            record(name, "pass" if ok else "FAIL", detail + " (limit 1e-5)")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pnoise_montecarlo():
    """OPTIONAL (PNOISE_MC=1; ~7 min/case): independent trnoise Monte Carlo
    cross-check via tools/pnoise_mc.py (TRNOISE current sources of PSD 4kT/R
    in parallel with noisy=0 resistors, N seeds, boxcar-anti-aliased Welch
    periodogram) on pnoise_r2_smooth.cir and pnoise_track_hold_slow.cir vs
    `.pnoise`. Pass = every 9-bin band within 0.5 dB (the card's target) AND
    the mean offset within 0.15 dB. Off by default: it is minutes of CPU."""
    name = "pnoise_montecarlo: trnoise MC vs .pnoise (r2_smooth, track_hold_slow, MOS1+diode buck, R3 closed-loop buck)"
    if os.environ.get("PNOISE_MC") != "1":
        record(name, "info", "skipped (set PNOISE_MC=1; see tools/pnoise_mc.py, docs/pss-pac.md)")
        return
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tools"))
        import pnoise_mc as mc
        mc.NGSPICE = _pac.NGSPICE if hasattr(_pac, "NGSPICE") else mc.NGSPICE
        dt, nfft = 200e-9, 8192
        T = [4000, 10000, 16000, 22000, 28000, 35000, 43000]
        worst, means = 0.0, []
        for tag, builder, fx, node, fn in (
                ("r2_smooth", mc.r2_smooth_mc_deck, "pnoise_r2_smooth.cir", "vx", None),
                ("th", mc.track_hold_mc_deck, "pnoise_track_hold_slow.cir", "vout", None)):
            body = "\n".join(mc.body_lines(fx))
            ca = (f"100e3 100e-6 {node} 2000 1 50 1e-6",)
            ref = lambda fl: mc.pnoise_reference(body, os.path.join(mc.DATA_ROOT, "test_" + tag), node, "Vin", ca, fl)
            rows, _, _ = mc.compare("test_mc_" + tag, builder, body, node, ca, ref, 32, 5e-9,
                                    20e-3, dt, nfft, 2e-3, T, ".pnoise")
            worst = max([worst] + [abs(r[3]) for r in rows])
            means.append(sum(r[3] for r in rows) / len(rows))
        # card #556: the same ensemble on a closed-loop PWM buck (R3, resistor noise, tools/pnoise_mc_r3.py) and on a MOS1 + diode buck
        # (tools/pnoise_mc_mos.py: M1 channel thermal and D1 shot noise injected with their instantaneous bias-dependent density)
        import pnoise_mc_mos as mm
        import pnoise_mc_r3 as m3
        for tag, mod, kw in (("m1buck", mm, dict(nseeds=16, K=300.0, TS=20e-9)), ("r3", m3, dict(nseeds=16, K=100.0, TS=10e-9))):
            rows, _, _ = mod.run(**kw)
            worst = max([worst] + [abs(r[3]) for r in rows])
            means.append(sum(r[3] for r in rows) / len(rows))
        ok = worst <= 0.5 and all(abs(m) <= lim for m, lim in zip(means, (0.15, 0.15, 0.25, 0.25)))
        record(name, "pass" if ok else "FAIL",
               f"worst band {worst:.3f} dB, mean offsets {[round(m, 3) for m in means]} dB (r2_smooth, th, m1buck, r3)")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


# --- card #557: `.pnoise` on the FREE-RUNNING hysteretic buck R5 (tests/fixtures/r5_freerun_plain.cir, trigger="q") ---------------------------------
R5_TRIG = 'method="shoot" solver="td" trigger="q" trigval=0.5 trigdir="rise"'
R5_HEAD = "427350 450e-6 vout 200 10 45 1e-4"
# the 9 Welch bins (610.35 Hz) of each band of tools/pnoise_mc_r5.py: first bin, last bin
R5_BANDS = {15258.8: (12817.3828125, 17700.1953125), 40283.2: (37841.796875, 42724.609375), 84838.9: (82397.4609375, 87280.2734375)}
# transient-noise ensemble (tools/pnoise_mc_r5.py: both resistors' noise as TRNOISE currents of PSD 4kT/R x K^2, K = 50, noisy and noiseless copy in one
# transient, reltol 1e-6, Hann/Welch 8192 x 200 ns, 9-bin band power) in V/rtHz, per trnoise step ts (ns), with the 1-sigma of the band average in dB.
# The injected noise is band-limited to 1/(2 ts): .pnoise is white to infinity, so the ensemble approaches it as ts -> 0 (power deficit ~ ts).
R5_ENS = {  # ts_ns: ({band: V/rtHz}, sigma_dB, seeds x ms)
    5.0: ({15258.8: 4.92996e-09, 40283.2: 5.06440e-09, 84838.9: 5.37193e-09}, 0.111, "16 x 20"),
    2.5: ({15258.8: 5.13275e-09, 40283.2: 5.38977e-09, 84838.9: 5.72514e-09}, 0.143, "12 x 16"),
}
R5_K = 50.0


def _r5_pnoise(tag, lo, hi, n, env=None, extra="", head=R5_HEAD, expect_run=True, timeout=300):
    outdir = os.path.join(DATA_ROOT, "r5_" + tag)
    os.makedirs(outdir, exist_ok=True)
    txt = os.path.join(outdir, "pn.txt")
    if os.path.exists(txt):
        os.remove(txt)
    deck = os.path.join(outdir, "pn.cir")
    with open(deck, "w") as fh:
        fh.write("\n".join(body_lines("r5_freerun_plain.cir")) + "\n")
        fh.write(f'.pnoise {head} vout Vin lin {n} {lo:.10g} {hi:.10g} {R5_TRIG}{extra}\n')
        fh.write(".control\nrun\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nwrdata pn.txt onoise_spectrum\nquit\n.endc\n.end\n")
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    rc, out, wall = run_ngspice(deck, timeout=timeout, env=e)
    if not expect_run:
        return rc, out, None
    if rc != 0 or "PNOISE-TD: done" not in out:
        raise CouldNotRun(f"R5 .pnoise {tag}: exit {rc}: {out[-1500:]}")
    _, cols = read_wrdata_real(txt, 1)
    return rc, out, cols[0]


def check_pnoise_freerun_r5():
    """card #557: `.pnoise ... trigger="q"` on the free-running R5 (autonomous hysteretic buck, fsw 422.7 kHz, validated band f <= 0.25 fsw).
    (1) vs the transient-noise ensemble at 15 / 40 / 85 kHz.  The ensemble's injected noise is band-limited to 1/(2 ts) while `.pnoise` is white to infinity
    (R5's comparator folds noise from many harmonics of fsw into the band), so the ensemble approaches `.pnoise` from below as ts -> 0: measured pooled
    deficits -5.3 / -2.8 / -1.1 dB at ts = 40 / 20 / 10 ns (6 seeds, 8 ms), -0.39 dB at 5 ns (16 x 20 ms) and +0.09 dB at 2.5 ns (12 x 16 ms).  Gate: at ts = 2.5 ns
    every band within 0.6 dB (~4 sigma of 0.143 dB) and the pooled mean within 0.3 dB, and the 5 ns ensemble lies below the 2.5 ns one (the trend that explains the
    deficit).  (2) the Floquet phase mode does not diverge in band: flat and finite from 10 Hz to 1.8 kHz (within 10 % of the 10 Hz value).  (3) red proof: the
    gate fails a +2 dB engine."""
    name = "pnoise(R5 free-running, trigger=q): .pnoise vs transient-noise ensemble (ts = 2.5 ns), finite at low f"
    try:
        eng = {}
        for c, (lo, hi) in R5_BANDS.items():
            rc, out, v = _r5_pnoise(f"band{int(c)}", lo, hi, 9)
            eng[c] = math.sqrt(sum(x * x for x in v) / len(v)) * R5_K
        parts, ds = [], []
        for c in R5_BANDS:
            d_s = 10 * math.log10(R5_ENS[2.5][0][c] ** 2 / eng[c] ** 2)
            d_l = 10 * math.log10(R5_ENS[5.0][0][c] ** 2 / eng[c] ** 2)
            ds.append((d_s, d_l))
            parts.append(f"{c / 1e3:.1f} kHz: .pnoise*K {eng[c]:.4e}, MC ts=5 ns {d_l:+.2f} dB, ts=2.5 ns {d_s:+.2f} dB")
        m_s, m_l = sum(d[0] for d in ds) / len(ds), sum(d[1] for d in ds) / len(ds)
        worst = max(abs(d[0]) for d in ds)
        rc, out, lowv = _r5_pnoise("low", 10.0, 1778.0, 6)
        low_ok = all(math.isfinite(x) and abs(x / lowv[0] - 1) < 0.10 for x in lowv)
        red = max(abs(d[0] + 2.0) for d in ds)
        ok = worst <= 0.6 and abs(m_s) <= 0.3 and m_l < m_s and low_ok and red > 0.6
        record(name, "pass" if ok else "FAIL",
               "; ".join(parts) + f"; worst band {worst:.2f} dB (bound 0.6), pooled mean {m_s:+.2f} dB (bound 0.3; ts=5 ns {m_l:+.2f}: lower, as the bandwidth argument says); "
               f"10 Hz..1.8 kHz flat within 10 %: {low_ok} ({lowv[0]:.4e} -> {lowv[-1]:.4e} V/rtHz); red proof: a +2 dB engine would read {red:.2f} dB (> 0.6 fails)")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pnoise_freerun_refusals():
    """card #557: what still refuses on a free-running orbit: the opt-in period-variation projection (PAC_FREERUN_PROJ=1, not transposed) and a sweep above the
    validated band (top frequency > 0.25 fsw = 105.7 kHz: 150 kHz is refused with the band message, no curve)."""
    name = "pnoise(R5 free-running) refusals: PAC_FREERUN_PROJ=1, and a sweep above 0.25 fsw"
    try:
        rc, out, v = _r5_pnoise("proj", 20000.0, 60000.0, 3, env={"PAC_FREERUN_PROJ": "1"}, expect_run=False)
        proj = "PNOISE-TD: REFUSED" in out and "PAC_FREERUN_PROJ" in out and "PNOISE-TD: done" not in out
        rc, out2, v2 = _r5_pnoise("band", 20000.0, 150000.0, 3, expect_run=False)
        band = "VALIDATED BAND" in out2 and "PNOISE-TD: done" not in out2
        record(name, "pass" if proj and band else "FAIL", f"projection refused {proj}; 150 kHz top frequency refused by the band guard {band}")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pnoise_freerun_mc():
    """OPTIONAL (PNOISE_MC=1; ~20 min with 3 ngspice at a time): the R5 ensemble itself, ts = 5 ns, 16 seeds x 20 ms (tools/pnoise_mc_r5.py); the bands must
    reproduce the stored ts = 5 ns values within 3 sigma."""
    name = "pnoise(R5 free-running): live ensemble (ts = 5 ns) reproduces the stored values"
    if os.environ.get("PNOISE_MC") != "1":
        record(name, "info", "skipped (set PNOISE_MC=1; tools/pnoise_mc_r5.py 16 50 5 15,40,85 20 with MC_RELTOL=1e-6)")
        return
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tools"))
        import pnoise_mc as mc
        import pnoise_mc_r5 as m5
        mc.NGSPICE = _pac.NGSPICE
        os.environ["MC_RELTOL"] = "1e-6"
        rows, _, _ = m5.run(nseeds=16, K=R5_K, TS=5e-9, targets=(15000, 40000, 85000), dur=20e-3)
        ref = R5_ENS[5.0][0]
        worst = max(abs(10 * math.log10((r[1] / ref[min(ref, key=lambda c: abs(c - r[0]))]) ** 2)) for r in rows)
        record(name, "pass" if worst <= 3 * R5_ENS[5.0][1] else "FAIL", f"worst band difference to the stored values {worst:.3f} dB (3 sigma = {3 * R5_ENS[5.0][1]:.2f})")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def main():
    check_pnoise_lti()
    check_pnoise_r3_buck()
    check_pnoise_r2_smooth()
    check_pnoise_track_hold_analytic()
    check_pnoise_track_hold_mos()
    check_pnoise_track_hold_mos(env={"PAC_NATIVE_HISTORY": "1"}, red=True)
    check_pnoise_fold()
    check_pnoise_fold(red=True)
    check_pnoise_montecarlo()
    check_pnoise_refuses_unsupported_device()
    check_pnoise_freerun_r5()
    check_pnoise_freerun_refusals()
    check_pnoise_freerun_mc()

    width = max(len(n) for n, _, _ in results)
    for n, state, detail in results:
        print(f"{state:15s} {n:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    npass = sum(s == "pass" for _, s, _ in results)
    ninfo = sum(s == "info" for _, s, _ in results)
    print(f"\n{npass}/{len(results) - ninfo} pass ({ninfo} info, not counted)")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
