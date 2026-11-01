#!/usr/bin/env python3
"""compopt.py — derivative-free loop-compensation optimizer on native `.stb` (card #471, plan T10.3).

Given a deck with a compensator's `.param` values named, bounds for each (log-scale for
R/C, the natural scale a decade sweep works in), a target (crossover f_c, phase margin,
optionally a gain-margin floor and "single falling crossing only"), and a probe, this
tunes the named parameters so the loop hits the target. Each evaluation is one native
`.stb` run (`~/Work/apps/ngspice_augmentation/<prefix>/bin/ngspice-aug`, docs/stb.md);
margins are computed in Python from the raw `T(f)` trace with `tools/stbcore.margins()`
(imported unmodified from this project's already-validated stb work, not reimplemented
here) rather than trusting native `.stb`'s own printed text, so the SAME crossing-finding
code this project's `stb.py`/native `.stb` use is applied to every candidate point,
including the "refuse a degenerate point" checks (a `low_f_positive_feedback` DC sign
flip, or more than one falling 0 dB crossing when the caller has not explicitly allowed
that) — a candidate that trips either is never treated as the accepted optimum, it is
just a heavily-penalised objective value the optimizer learns to move away from.

Optimizer: Nelder-Mead (stdlib only, no numpy/scipy on this machine) in log10-parameter
space for every log-scale parameter, with random restarts (`random.Random(seed)`, so a
run is reproducible) from a perturbation of the best point found so far. This is
derivative-free by construction — appropriate since each "function value" is a subprocess
ngspice run, not a closed-form expression.

Design note on `.pstb` readiness: the only ngspice-facing code is `run_native_stb()`
(build a `.control` block, run one analysis, `wrdata` one vector, parse it) and
`set_params()`/`to_native_probe_deck()` (text-level deck preparation). Swapping `.stb` for
`.pstb` later (plan T5.2) means changing the analysis card written by `run_native_stb()`
and, if `.pstb` needs a `.pss` pre-pass, adding that to the same control block — the
`Param`/`Target`/`evaluate`/`nelder_mead`/`optimize` machinery above it does not know or
care which analysis produced the `T(f)` trace it is handed.

Never touches tools/stb.py, tools/stbcore.py or tools/fra*.py — they are imported,
unmodified, and reused (stbcore.margins for crossing-finding; stb.read_wrdata/prepare_deck
and fra.sweep_loop/write_loop_results for the two independent verification paths).

Outputs (deck copies, ngspice logs, wrdata dumps, JSON reports) never go in the repo:
default output root is ~/Work/data/ngspice_augmentation/compopt/<run-id>/.

    python3 tools/compopt.py DECK --lgprobe-inst XP \\
        --param cc=1e-13:3e-11:log --param c1=1e-14:1e-11:log \\
        --start cc=1e-12 --start c1=0.2e-12 \\
        --target fc=7.6e6 --target pm=65 \\
        --sweep "dec 30 10 1e9" --verify stb

Exit status: 0 = target reached within tolerance and independently confirmed,
1 = ran but the target is reported infeasible (never printed as a success),
2 = could not run at all (ngspice missing, deck error, bad arguments).
"""
import argparse
import cmath
import datetime
import json
import math
import os
import random
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import stb  # noqa: E402  (unmodified — read_wrdata, prepare_deck, CouldNotRun, run())
import stbcore  # noqa: E402  (unmodified — tian, margins, db, unwrapped_phase_deg)

DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/compopt")
DEFAULT_NGSPICE_NATIVE = os.path.expanduser(
    "~/Work/apps/ngspice_augmentation/aaba50ba3f49-s12/bin/ngspice-aug")
NGSPICE_NATIVE = os.environ.get("NGSPICE_AUG") or os.environ.get("NGSPICE") or DEFAULT_NGSPICE_NATIVE  # 2026-10-02: follow the build under test
NGSPICE_STOCK = os.environ.get("NGSPICE", "ngspice")

BIG = 1.0e4          # base cost assigned to a refused/could-not-run evaluation
FC_TOL_REL = 0.05     # spec: crossover within 5%
PM_TOL_DEG = 2.0      # spec: phase margin within 2 deg

LGPROBE_INCLUDE_RE = re.compile(r'(?im)^\s*\.include\s+.*lgprobe\.lib\s*$')

# Standard IEC 60063 E24/E96 decade mantissas, for the optional "round final
# values to a standard series" report step (not used inside the optimizer).
E24 = [1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0,
       3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1]
E96 = [1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21, 1.24, 1.27, 1.30,
       1.33, 1.37, 1.40, 1.43, 1.47, 1.50, 1.54, 1.58, 1.62, 1.65, 1.69, 1.74,
       1.78, 1.82, 1.87, 1.91, 1.96, 2.00, 2.05, 2.10, 2.16, 2.21, 2.26, 2.32,
       2.37, 2.43, 2.49, 2.55, 2.61, 2.67, 2.74, 2.80, 2.87, 2.94, 3.01, 3.09,
       3.16, 3.24, 3.32, 3.40, 3.48, 3.57, 3.65, 3.74, 3.83, 3.92, 4.02, 4.12,
       4.22, 4.32, 4.42, 4.53, 4.64, 4.75, 4.87, 4.99, 5.11, 5.23, 5.36, 5.49,
       5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65, 6.81, 6.98, 7.15, 7.32,
       7.50, 7.68, 7.87, 8.06, 8.25, 8.45, 8.66, 8.87, 9.09, 9.31, 9.53, 9.76]


class CouldNotRun(RuntimeError):
    """ngspice failed to run, or produced no usable T(f) trace."""


# --------------------------------------------------------------------- utils

def round_series(value, series="E96"):
    """Round a positive value to the nearest standard E24/E96 decade value."""
    if value <= 0:
        return value
    table = E24 if series == "E24" else E96
    decade = 10 ** math.floor(math.log10(value))
    mantissa = value / decade
    candidates = [(s, s * decade) for s in table] + [(table[0] * 10, table[0] * decade * 10)]
    _, best = min(candidates, key=lambda sc: abs(math.log(sc[0]) - math.log(mantissa)))
    return best


def set_params(deck_text, values):
    """Return deck_text with every `name=value` token for each name in `values`
    (case-insensitive, whole-token match on any line, typically a `.param` line)
    replaced by the new value, formatted plainly (ngspice accepts scientific
    notation anywhere a suffixed value like "10k" is accepted). Raises
    CouldNotRun if a name is not found at all — silently doing nothing to an
    unmatched parameter would be exactly the kind of silent number this
    project's rules forbid."""
    text = deck_text
    for name, value in values.items():
        pat = re.compile(r'(?i)(?<![A-Za-z0-9_.])(' + re.escape(name) + r')(\s*=\s*)([^\s]+)')
        text, n = pat.subn(lambda m: m.group(1) + m.group(2) + f"{value:.10g}", text)
        if n == 0:
            raise CouldNotRun(f"param {name!r} not found in deck text (no name=value token)")
    return text


def to_native_probe_deck(deck_text, lgprobe_inst, probe_name="Vprobe"):
    """Replace an `<inst> a b lgprobe` instance line with a plain 0 V V-source, and
    drop the lgprobe.lib .include line: native `.stb`'s probe= resolves to an actual
    V-source instance, not a subcircuit (docs/stb.md, "probe=" table). Same
    transform as tests/test_stb_native.py's plain_variant() (a test file, not
    imported here — reimplemented as this module's own small, self-contained
    helper, same as that script documents for itself)."""
    out = []
    hit = False
    for line in deck_text.splitlines():
        if LGPROBE_INCLUDE_RE.match(line):
            continue
        f = line.split()
        if f and f[0].lower() == lgprobe_inst.lower() and len(f) >= 4 and f[3].lower() == "lgprobe":
            n1, n2 = f[1], f[2]
            out.append(f"{probe_name} {n1} {n2} DC 0")
            hit = True
            continue
        out.append(line)
    if not hit:
        raise CouldNotRun(f"lgprobe instance {lgprobe_inst!r} not found in deck")
    return "\n".join(out) + "\n"


def _strip_end(text):
    return re.sub(r"\n\.end\s*$", "\n", text.rstrip("\n") + "\n")


# ---------------------------------------------------------- native `.stb` run

def run_native_stb(deck_text, probe, sweep, outdir, ngspice=NGSPICE_NATIVE, timeout=120):
    """Run native `.stb` on deck_text (already parametrized, plain-V-source probe);
    return (freqs, T, log_text). Raises CouldNotRun on any ngspice error, a missing
    output vector (e.g. `.stb` refused a bad fstart/fstop/n — docs/stb.md's
    "Frequency and point-count validation" table), or too few points."""
    os.makedirs(outdir, exist_ok=True)
    body = _strip_end(deck_text)
    tfile = os.path.join(outdir, "t.txt")
    if os.path.exists(tfile):
        os.remove(tfile)
    ctl = (
        "\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
        f"stb {sweep} probe={probe}\nwrdata {tfile} stb1.t\nquit\n.endc\n.end\n"
    )
    wrapper = os.path.join(outdir, "wrapper.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + ctl)
    try:
        proc = subprocess.run([ngspice, "-b", wrapper], capture_output=True, text=True, timeout=timeout)
    except (FileNotFoundError, PermissionError) as e:
        raise CouldNotRun(f"cannot start {ngspice}: {e}") from e
    except subprocess.TimeoutExpired as e:
        raise CouldNotRun(f"ngspice timed out after {timeout}s: {e}") from e
    log = proc.stdout + proc.stderr
    with open(os.path.join(outdir, "log.txt"), "w") as fh:
        fh.write(log)
    errs = [ln for ln in log.splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")
    if not os.path.exists(tfile):
        raise CouldNotRun("ngspice wrote no loopgain/T vector (stb card refused?)")
    freqs, (t,) = stb.read_wrdata(tfile, 1)
    if len(freqs) < 2:
        raise CouldNotRun(f"only {len(freqs)} frequency point(s) — cannot find a crossing")
    return freqs, t, log


# ------------------------------------------------------------------- Param

class Param:
    """One tunable `.param` name with bounds. log=True (default; the right choice
    for R/C component values) works in log10 space, so a step of the optimizer
    means a fixed *ratio* change, not a fixed absolute change."""

    def __init__(self, name, lo, hi, log=True):
        if lo <= 0 or hi <= 0:
            raise ValueError(f"param {name!r}: bounds must be positive ({lo}, {hi})")
        if not lo < hi:
            raise ValueError(f"param {name!r}: need lo < hi ({lo}, {hi})")
        self.name, self.lo, self.hi, self.log = name, float(lo), float(hi), log

    def to_x(self, value):
        return math.log10(value) if self.log else value

    def from_x(self, x):
        v = (10 ** x) if self.log else x
        return min(max(v, self.lo), self.hi)

    def bounds_x(self):
        return (self.to_x(self.lo), self.to_x(self.hi))

    @staticmethod
    def parse(spec):
        """"name=lo:hi" or "name=lo:hi:lin" (default log)."""
        name, _, rest = spec.partition("=")
        parts = rest.split(":")
        if len(parts) < 2:
            raise ValueError(f"bad --param {spec!r}, want name=lo:hi[:lin|log]")
        lo, hi = float(parts[0]), float(parts[1])
        log = True
        if len(parts) >= 3:
            log = parts[2].lower() != "lin"
        return Param(name, lo, hi, log=log)


# ----------------------------------------------------------------- Evaluation

class Evaluation:
    """The full record of one candidate point: status is one of
    "ok" | "could-not-run" | "positive-feedback" | "multi-falling" | "no-falling-crossing".
    Only "ok" evaluations with within_tol=True are ever reported as a success."""

    def __init__(self, status, cost, values, margins=None, freqs=None, T=None,
                 fc=None, pm=None, rel_fc_err=None, abs_pm_err=None, within_tol=False, detail=""):
        self.status = status
        self.cost = cost
        self.values = dict(values)
        self.margins = margins
        self.freqs = freqs
        self.T = T
        self.fc = fc
        self.pm = pm
        self.rel_fc_err = rel_fc_err
        self.abs_pm_err = abs_pm_err
        self.within_tol = within_tol
        self.detail = detail

    def to_json(self):
        return {
            "status": self.status, "cost": self.cost, "values": self.values,
            "fc": self.fc, "pm": self.pm, "rel_fc_err": self.rel_fc_err,
            "abs_pm_err": self.abs_pm_err, "within_tol": self.within_tol, "detail": self.detail,
            "gm": (self.margins["gm"] if self.margins else None),
            "low_f_positive_feedback": (self.margins["low_f_positive_feedback"] if self.margins else None),
            "n_falling_crossings": (sum(1 for c in self.margins["pm_crossings"] if c["falling"])
                                     if self.margins else None),
        }


def evaluate(base_deck_text, lgprobe_inst, sweep, target, values, outdir,
             ngspice=NGSPICE_NATIVE, allow_multiple_falling=False,
             weight_fc=1.0, weight_pm=1.0, weight_gm=0.3, timeout=120):
    """Run one native `.stb` evaluation of `values` and score it against `target`
    (dict: fc, pm, optional gm_min, optional single_falling — kept for API symmetry
    with allow_multiple_falling, which is the actual switch used below).

    A point that triggers native `.stb`'s own `low_f_positive_feedback` warning, or
    that has more than one falling 0 dB crossing when the caller has not explicitly
    passed allow_multiple_falling=True, is refused: it gets a large, finite cost (so
    Nelder-Mead still has a gradient to follow away from it) and is marked with a
    status other than "ok", so it can never be reported as the accepted optimum."""
    deck_native = to_native_probe_deck(set_params(base_deck_text, values), lgprobe_inst)
    try:
        freqs, T, _log = run_native_stb(deck_native, "Vprobe", sweep, outdir, ngspice=ngspice, timeout=timeout)
    except CouldNotRun as e:
        return Evaluation("could-not-run", BIG * 3, values, detail=str(e))

    m = stbcore.margins(freqs, T)

    if m["low_f_positive_feedback"]:
        return Evaluation("positive-feedback", BIG * 2, values, margins=m, freqs=freqs, T=T, detail=m["note"])

    falling = [c for c in m["pm_crossings"] if c["falling"]]

    if not falling:
        # No valid falling 0 dB crossing at all (rising-only, or none). Give the
        # search something to follow: how far |T| is from 0 dB at the frequency
        # nearest the target crossover (want that to go to 0, creating a crossing).
        idx = min(range(len(freqs)), key=lambda i: abs(math.log10(freqs[i]) - math.log10(target["fc"])))
        db_near = stbcore.db(T[idx])
        return Evaluation("no-falling-crossing", BIG + db_near ** 2, values, margins=m, freqs=freqs, T=T,
                           detail="no falling 0 dB crossing near the target frequency")

    if len(falling) > 1 and not allow_multiple_falling:
        return Evaluation("multi-falling", BIG * 1.5, values, margins=m, freqs=freqs, T=T,
                           detail=f"{len(falling)} falling 0 dB crossings (target asked for a single one)")

    chosen = min(falling, key=lambda c: c["pm_deg"])  # worst-case, same convention as native/stbcore
    fc, pm = chosen["freq"], chosen["pm_deg"]
    err_fc = math.log10(fc / target["fc"])
    err_pm = (pm - target["pm"]) / 60.0
    cost = weight_fc * err_fc ** 2 + weight_pm * err_pm ** 2

    gm_min = target.get("gm_min")
    gm_ok = True
    if gm_min is not None and m["gm"] is not None and m["gm"]["gm_db"] < gm_min:
        gm_ok = False
        cost += weight_gm * ((gm_min - m["gm"]["gm_db"]) / 6.0) ** 2

    rel_fc_err = abs(fc - target["fc"]) / target["fc"]
    abs_pm_err = abs(pm - target["pm"])
    within_tol = (rel_fc_err <= FC_TOL_REL) and (abs_pm_err <= PM_TOL_DEG) and gm_ok
    return Evaluation("ok", cost, values, margins=m, freqs=freqs, T=T, fc=fc, pm=pm,
                       rel_fc_err=rel_fc_err, abs_pm_err=abs_pm_err, within_tol=within_tol)


# ------------------------------------------------------------- Nelder-Mead

def nelder_mead(f, x0, bounds, step=0.3, max_iter=80, tol=1e-9):
    """Minimal stdlib Nelder-Mead, coordinates already in the optimizer's working
    space (log10 for a log-scale Param). f(x) -> float. bounds: list of (lo, hi) in
    that same space, clamped after every simplex move (simple projection, adequate
    for a smooth-ish loop-shaping cost surface). Returns (best_x, best_cost,
    n_evals, history) where history is the best cost seen after each iteration."""
    n = len(x0)

    def clamp(x):
        return [min(max(xi, lo), hi) for xi, (lo, hi) in zip(x, bounds)]

    simplex = [clamp(list(x0))]
    for i in range(n):
        xi = list(x0)
        xi[i] += step
        simplex.append(clamp(xi))
    fvals = [f(x) for x in simplex]
    n_evals = n + 1
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5
    history = []

    for _it in range(max_iter):
        order = sorted(range(len(simplex)), key=lambda i: fvals[i])
        simplex = [simplex[i] for i in order]
        fvals = [fvals[i] for i in order]
        history.append(fvals[0])
        if fvals[-1] - fvals[0] < tol:
            break
        centroid = [sum(s[j] for s in simplex[:-1]) / n for j in range(n)]
        worst = simplex[-1]

        xr = clamp([centroid[j] + alpha * (centroid[j] - worst[j]) for j in range(n)])
        fr = f(xr)
        n_evals += 1
        if fvals[0] <= fr < fvals[-2]:
            simplex[-1], fvals[-1] = xr, fr
        elif fr < fvals[0]:
            xe = clamp([centroid[j] + gamma * (xr[j] - centroid[j]) for j in range(n)])
            fe = f(xe)
            n_evals += 1
            if fe < fr:
                simplex[-1], fvals[-1] = xe, fe
            else:
                simplex[-1], fvals[-1] = xr, fr
        else:
            xc = clamp([centroid[j] + rho * (worst[j] - centroid[j]) for j in range(n)])
            fc_ = f(xc)
            n_evals += 1
            if fc_ < fvals[-1]:
                simplex[-1], fvals[-1] = xc, fc_
            else:
                best = simplex[0]
                for i in range(1, len(simplex)):
                    simplex[i] = clamp([best[j] + sigma * (simplex[i][j] - best[j]) for j in range(n)])
                    fvals[i] = f(simplex[i])
                    n_evals += 1

    order = sorted(range(len(simplex)), key=lambda i: fvals[i])
    return simplex[order[0]], fvals[order[0]], n_evals, history


# ---------------------------------------------------------------- optimize()

class OptimizeResult:
    def __init__(self, feasible, best, all_evals, n_evals, outdir, restarts_done):
        self.feasible = feasible
        self.best = best  # Evaluation
        self.all_evals = all_evals
        self.n_evals = n_evals
        self.outdir = outdir
        self.restarts_done = restarts_done

    def to_json(self):
        return {
            "feasible": self.feasible, "best": self.best.to_json(),
            "n_evals": self.n_evals, "outdir": self.outdir, "restarts_done": self.restarts_done,
        }


def optimize(base_deck_text, lgprobe_inst, params, start, target, sweep, outdir,
             ngspice=NGSPICE_NATIVE, restarts=2, max_iter=80, seed=0,
             allow_multiple_falling=False, step=0.3):
    """Tune `params` (list of Param) of `base_deck_text` (lgprobe-instance form,
    a checked-in tests/fixtures/compopt_*.cir) starting from `start` (dict
    name->value) to reach `target` (dict fc, pm, optional gm_min, optional
    single_falling). Every evaluation's ngspice run lands under
    outdir/eval_<k>/; a run-level trace.csv and result.json are written directly
    under outdir. Returns an OptimizeResult; .feasible is False whenever the
    target was not reached within FC_TOL_REL/PM_TOL_DEG — this function never
    marks a refused (positive-feedback / multi-falling) or could-not-run point
    as the reported optimum, even if random search stumbled onto a low
    residual-looking cost there (both are floored well above any achievable
    "ok" cost by BIG's multipliers)."""
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    counter = {"k": 0}
    all_evals = []
    trace_path = os.path.join(outdir, "trace.csv")
    with open(trace_path, "w") as fh:
        fh.write("eval,restart," + ",".join(p.name for p in params) + ",status,cost,fc,pm,rel_fc_err,abs_pm_err\n")

    def objective(x, restart_idx):
        counter["k"] += 1
        k = counter["k"]
        values = {p.name: p.from_x(xi) for p, xi in zip(params, x)}
        ev_outdir = os.path.join(outdir, f"eval_{k:04d}")
        ev = evaluate(base_deck_text, lgprobe_inst, sweep, target, values, ev_outdir,
                       ngspice=ngspice, allow_multiple_falling=allow_multiple_falling)
        all_evals.append(ev)
        with open(trace_path, "a") as fh:
            vals = ",".join(f"{values[p.name]:.6g}" for p in params)
            fh.write(f"{k},{restart_idx},{vals},{ev.status},{ev.cost:.6g},"
                     f"{ev.fc if ev.fc else ''},{ev.pm if ev.pm else ''},"
                     f"{ev.rel_fc_err if ev.rel_fc_err is not None else ''},"
                     f"{ev.abs_pm_err if ev.abs_pm_err is not None else ''}\n")
        return ev.cost

    bounds = [p.bounds_x() for p in params]
    x0 = [p.to_x(start[p.name]) for p in params]

    best_x, best_cost = None, math.inf
    restarts_done = 0
    for r in range(max(1, restarts)):
        if r == 0:
            x_init = [min(max(xi, lo), hi) for xi, (lo, hi) in zip(x0, bounds)]
        else:
            x_init = [min(max(xi + rng.uniform(-0.6, 0.6), lo), hi)
                      for xi, (lo, hi) in zip(best_x, bounds)]
        xb, cb, _n, _hist = nelder_mead(lambda x: objective(x, r), x_init, bounds,
                                         step=step, max_iter=max_iter)
        restarts_done += 1
        if cb < best_cost:
            best_x, best_cost = xb, cb

    final_values = {p.name: p.from_x(xi) for p, xi in zip(params, best_x)}
    final_outdir = os.path.join(outdir, "final")
    final_eval = evaluate(base_deck_text, lgprobe_inst, sweep, target, final_values, final_outdir,
                           ngspice=ngspice, allow_multiple_falling=allow_multiple_falling)
    all_evals.append(final_eval)

    feasible = (final_eval.status == "ok") and final_eval.within_tol
    result = OptimizeResult(feasible, final_eval, all_evals, counter["k"] + 1, outdir, restarts_done)
    with open(os.path.join(outdir, "result.json"), "w") as fh:
        json.dump(result.to_json(), fh, indent=2)
    return result


# ---------------------------------------------------- independent verification

def verify_with_stb_script(base_deck_text, lgprobe_inst, final_values, sweep, outdir,
                            ngspice=NGSPICE_STOCK, margins_fn=stbcore.margins):
    """Independent confirmation: tools/stb.py (the script-based Tian testbench,
    unmodified, validated in tests/test_stb.py) run on unmodified stock ngspice,
    on the SAME lgprobe-form deck the native optimizer was fed (parametrized to
    the optimizer's final values) — a different code path and, with NGSPICE
    pointed at stock ngspice, a different binary than the one that produced the
    reported optimum. `margins_fn` defaults to the real stbcore.margins; the
    red-proof test in tests/test_compopt.py passes a deliberately corrupted one
    (a sign flip) to show this check is not vacuous."""
    os.makedirs(outdir, exist_ok=True)
    deck_text = set_params(base_deck_text, final_values)
    # the checked-in fixture's .include is relative to tests/fixtures/ (../../stb/...);
    # this copy lives under ~/Work/data/, so make it absolute (same fix as
    # fra_curves.make_r4_variant uses for its own deck copies).
    abs_lib_lg = os.path.join(ROOT, "stb", "lgprobe.lib")
    deck_text = re.sub(r"\.include\s+\S*lgprobe\.lib", f".include {abs_lib_lg}", deck_text)
    deck_path = os.path.join(outdir, "verify_deck.cir")
    with open(deck_path, "w") as fh:
        fh.write(deck_text)
    old_ngspice = stb.NGSPICE
    stb.NGSPICE = ngspice
    try:
        res = stb.run(deck_path, lgprobe_inst, sweep, os.path.join(outdir, "stb_py"))
    finally:
        stb.NGSPICE = old_ngspice
    margins = margins_fn(res["freqs"], res["T"])
    return margins, res


def verify_with_fra(switching_deck_path, final_values, lgprobe_inst, freqs, outdir,
                     settle_time=300e-6, workers=1):
    """Independent confirmation for R3: apply `final_values` to a COPY of the real
    switching deck (refs/r3_vm_buck.cir — read, never written to), run
    tools/fra.py's transient-injection Tian/Middlebrook harness (imported,
    unmodified) at a handful of frequencies around the optimized crossover, using
    the SAME converter settings tools/fra_curves.py documents (POINT_COMMON,
    AMPLITUDE) for this exact deck family. Returns (summary, points, deck_path);
    summary["combined_margins"] is the FRA-measured crossover/PM/GM to compare
    against the averaged-model .stb result — expect a real difference (averaged
    != switching near fs/2, spec's own wording), not agreement to stb's 1e-9."""
    sys.path.insert(0, HERE)
    import fra  # noqa
    import fra_curves  # noqa
    os.makedirs(outdir, exist_ok=True)
    with open(switching_deck_path) as fh:
        text = fh.read()
    text = set_params(text, final_values)
    abs_lib = os.path.join(ROOT, "stb", "lgprobe.lib")
    text, n = re.subn(r"\.include\s+\.\./stb/lgprobe\.lib", f".include {abs_lib}", text, count=1)
    if n != 1:
        raise CouldNotRun("could not rewrite the lgprobe .include path in the switching-deck copy")
    deck_path = os.path.join(outdir, "switching_deck_optimized.cir")
    with open(deck_path, "w") as fh:
        fh.write(text)
    kw = dict(fra_curves.POINT_COMMON)
    points = fra.sweep_loop(deck_path, lgprobe_inst, freqs, fra_curves.AMPLITUDE,
                             settle_time=settle_time, outdir_root=os.path.join(outdir, "sweep"),
                             workers=workers, **kw)
    summary = fra.write_loop_results(points, freqs, outdir)
    return summary, points, deck_path


# --------------------------------------------------------------------- CLI

def _parse_target(pairs):
    target = {}
    for kv in pairs:
        k, _, v = kv.partition("=")
        k = k.strip().lower()
        if k in ("fc", "pm", "gm_min", "gm"):
            target["gm_min" if k == "gm" else k] = float(v)
        elif k == "single_falling":
            target["single_falling"] = v.lower() not in ("0", "false", "no")
        else:
            raise ValueError(f"unknown --target key {k!r}")
    if "fc" not in target or "pm" not in target:
        raise ValueError("--target needs at least fc=<Hz> and pm=<deg>")
    return target


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("deck", help="lgprobe-instance-form deck, e.g. tests/fixtures/compopt_r1b_detuned.cir")
    ap.add_argument("--lgprobe-inst", required=True, help="lgprobe instance name in the deck, e.g. XP")
    ap.add_argument("--param", action="append", required=True, default=[],
                     help='"name=lo:hi[:lin]" (default log scale), repeatable')
    ap.add_argument("--start", action="append", default=[], help='"name=value", repeatable; '
                     "default is the bound's geometric (or linear) midpoint")
    ap.add_argument("--target", action="append", required=True, default=[],
                     help='"fc=<Hz>" | "pm=<deg>" | "gm=<dB>" | "single_falling=1", repeatable')
    ap.add_argument("--sweep", default="dec 30 10 1e9", help='native .stb frequency sweep, e.g. "dec 30 10 1e9"')
    ap.add_argument("--restarts", type=int, default=2)
    ap.add_argument("--max-iter", type=int, default=80)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--allow-multiple-falling", action="store_true")
    ap.add_argument("--ngspice", default=NGSPICE_NATIVE)
    ap.add_argument("--verify", choices=("none", "stb"), default="none",
                     help="independent confirmation after optimizing: stb = tools/stb.py on stock ngspice")
    ap.add_argument("--out", help="output directory (default: under ~/Work/data/ngspice_augmentation/compopt)")
    args = ap.parse_args()

    try:
        params = [Param.parse(s) for s in args.param]
        target = _parse_target(args.target)
        start = {}
        for s in args.start:
            k, _, v = s.partition("=")
            start[k] = float(v)
        for p in params:
            if p.name not in start:
                start[p.name] = math.sqrt(p.lo * p.hi) if p.log else (p.lo + p.hi) / 2
    except ValueError as e:
        print(f"COULD NOT RUN: bad arguments: {e}", file=sys.stderr)
        return 2

    with open(args.deck) as fh:
        base_deck_text = fh.read()

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = args.out or os.path.join(DATA_ROOT, f"{os.path.splitext(os.path.basename(args.deck))[0]}-{stamp}")

    try:
        result = optimize(base_deck_text, args.lgprobe_inst, params, start, target, args.sweep, outdir,
                           ngspice=args.ngspice, restarts=args.restarts, max_iter=args.max_iter,
                           seed=args.seed, allow_multiple_falling=args.allow_multiple_falling)
    except CouldNotRun as e:
        print(f"COULD NOT RUN: {e}", file=sys.stderr)
        return 2

    b = result.best
    print(f"path: {result.outdir}")
    print(f"evaluations: {result.n_evals} ({result.restarts_done} restart(s))")
    print("final values:")
    for p in params:
        v = b.values[p.name]
        print(f"  {p.name} = {v:.6g}  (E96 ~= {round_series(v):.6g}, E24 ~= {round_series(v, 'E24'):.6g})")
    if b.status == "ok":
        print(f"fc = {b.fc:.6g} Hz (target {target['fc']:.6g}, rel err {b.rel_fc_err*100:.2f}%)")
        print(f"PM = {b.pm:.3f} deg (target {target['pm']:.3f}, abs err {b.abs_pm_err:.3f} deg)")
    else:
        print(f"status: {b.status} -- {b.detail}")

    if result.feasible:
        print("TARGET REACHED (within 5% fc / 2 deg PM)")
    else:
        print("TARGET NOT REACHED -- reporting INFEASIBLE (best point found, not a success)")

    if args.verify == "stb" and result.feasible:
        vmargins, _res = verify_with_stb_script(base_deck_text, args.lgprobe_inst, b.values,
                                                  args.sweep, os.path.join(outdir, "verify"))
        pm = vmargins["pm"]
        print("independent check (tools/stb.py, stock ngspice):")
        print(f"  PM = {pm['pm_deg']:.3f} deg @ {pm['freq']:.6g} Hz" if pm else "  PM: no 0 dB crossing")

    return 0 if result.feasible else 1


if __name__ == "__main__":
    sys.exit(main())
