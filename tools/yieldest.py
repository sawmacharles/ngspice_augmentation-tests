#!/usr/bin/env python3
"""yieldest.py -- fast / high-sigma yield (failure-probability) estimation
(card #474, plan T10.6). Companion to tools/mc.py (imported, not modified --
this card's hard rule): mc.py draws N independent samples and computes plain
statistics; at 4-6 sigma targets plain Monte Carlo needs millions of ngspice
runs to bound Pf honestly (at Pf~3e-5, a 20%-relative-CI 95% estimate needs
roughly N ~ (1.96/0.20)^2 / 3e-5 ~= 3.2e6 runs -- see docs/yield.md for the
worked number), which is not affordable. This tool estimates Pf with far
fewer runs, using the standard variance-reduction techniques from
research/feature-survey.md item #6 (Jonsson & Lelong 2021's importance
sampling; the industry "FMC"/"HSMC" family's worst-case-distance search),
not an ML surrogate.

## Method 1: mean-shift importance sampling around the worst-case-distance
   (most-probable-failure) point

All sampling happens in a STANDARDIZED parameter space: every declared
Gaussian parameter x_i (nominal m_i, sigma s_i) is represented as
u_i = (x_i - m_i) / s_i, so the *nominal* (design) sampling distribution is
always the standard multivariate normal N(0, I) in u-space, regardless of
each parameter's physical units. A "failure" is a single scalar measure
(read back via mc.py's own MCMEAS convention) crossing a threshold.

Step 1 -- worst-case-distance (WCD) search, gradient-free: a deterministic
set of directions in u-space (coordinate axes; for 2 parameters a dense
15-degree angular grid; for 3+ parameters the coordinate axes plus the
all-sign unit-diagonal directions) is searched outward from the origin by
expand-then-bisect line search along each ray, giving a boundary radius per
direction. The smallest such radius is refined further by a gradient-free
directional pattern (Hooke-Jeeves-style) search on the direction vector
itself. This gives beta = ||u*|| (the "worst-case distance", exactly what
FORM/HSMC-style tools call it) and u* (the most-probable-failure point,
MPP) -- see `search_mpp`.

Step 2 -- mean-shift importance sampling: draw N samples from N(u*, I) (same
covariance, shifted mean -- "mean-shift" IS, the simplest, most standard
variant), evaluate each on the real circuit, and weight every sample by the
likelihood ratio w(u) = f0(u)/f_IS(u) = exp(-u.u* + 0.5*||u*||^2) (f0 is the
true N(0,I) density, f_IS is N(u*,I)). Pf_hat = mean(w * 1{fail}) over the N
IS draws -- a standard, unbiased importance-sampling estimator (see
`importance_sampling`). Most of the probability mass under N(u*, I) sits
right where the failure region is, so most draws contribute a nonzero
weighted term, instead of almost all draws being an uninformative "pass"
the way plain MC's draws from N(0,I) are at 4-6 sigma.

## Method 2 (cross-check): scaled-sigma sampling (SSS) extrapolation

Independent of the MPP search's direction/weights (it does reuse the found
beta only to pick a sensible scale range, see `default_sss_scales`): sample
plain MC at several *enlarged* standard deviations (u ~ N(0, c^2 I) for a
few c > 1), where the failure probability Pf_c is large enough to measure
directly and cheaply (thousands, not millions, of runs). Assuming (as
WCD/FORM does) the failure boundary is approximately a hyperplane at
distance beta from the origin, Pf_c ~= Phi(-beta/c) exactly for a true
hyperplane, so beta_hat_c = -c * Phi^-1(Pf_c_hat) is recovered from each
scale's plain-MC estimate; the per-scale estimates are combined by weighted
least squares (regressing Phi^-1(Pf_c) against 1/c through the origin,
weighted by the delta-method variance of each scale's own binomial
proportion -- found empirically to matter: a scale with only a handful of
failures is much noisier than one with dozens, and a plain unweighted
average of beta_hat_c let a noisy scale swing the answer). This needs no
gradient and no importance weights at all -- a genuinely independent second
method, per this card's brief ("scaled-sigma sampling (SSS) extrapolation,
as a cross-check") -- see `scaled_sigma_sampling`.

## Honesty check: the two methods should roughly agree

`estimate_yield` compares the IS estimate against the SSS estimate; a large
disagreement (>3x either way) is reported as `limitation_warning`, not
silently absorbed into "the answer" -- this is this tool's only defence
against a failure region the direction-limited WCD search cannot see (a
narrow region between two search rays, or a second disjoint region far from
the one direction search found) -- see docs/yield.md Limits and
tests/test_yield.py's disjoint-region red-proof, which is built to trigger
exactly this.

## What this tool does NOT do

- Only `dist: "gauss"` parameters (the whole likelihood-ratio/beta machinery
  is Gaussian-specific); `mc.py`'s "unif" dist is not supported here.
- Only a single scalar pass/fail measure per estimate (`measure` + `fail`
  {"op": "lt"|"gt", "threshold": v}) -- no multi-measure combined yield yet.
- The WCD search assumes each searched ray crosses pass->fail at most once
  (standard FORM assumption); a non-monotonic ray can mislead the bisection.
- No ML surrogate (deliberately -- research/feature-survey.md item #6's own
  ranking prefers a WCD/importance-sampling wrapper over a surrogate for a
  first version; a surrogate is future work, not attempted here).

Every ngspice evaluation is a separate process (mc.run_one, same as
tools/mc.py), run through a bounded ThreadPoolExecutor (<=4 concurrent,
enforced regardless of what a spec asks for -- this card's hard rule).
Outputs (decks, manifest.json) go under `~/Work/data/ngspice_augmentation/
yield/<run>/`, never the repository.
"""
import argparse
import concurrent.futures
import hashlib
import itertools
import json
import math
import os
import random
import statistics
import sys
import threading
import time
from statistics import NormalDist

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import mc  # noqa: E402 -- imported, not modified (card #474 hard rule): reuses
# mc.render_deck (the @@TOKEN@@ substitution) and mc.run_one (one ngspice
# subprocess per evaluation, MCMEAS stdout parsing).

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/yield")

Z95 = NormalDist().inv_cdf(0.975)  # 95% two-sided normal quantile, ~1.95996
# (95%, not mc.py's 99%: yield/Pf reporting conventionally uses 95% CIs --
# stated here, not silently different from mc.py's own convention)


# --------------------------------------------------------------------- #
# Standardized <-> physical parameter values
# --------------------------------------------------------------------- #

def _sigma_of(param):
    if "pelgrom" in param:
        p = param["pelgrom"]
        return p["coeff"] / math.sqrt(p["w"] * p["l"])
    return param["sigma"]


def _seed_int(*parts):
    h = hashlib.sha256("|".join(str(p) for p in parts).encode()).digest()
    return int.from_bytes(h[:8], "big")


def check_params_supported(params):
    for p in params:
        if p.get("dist", "gauss") != "gauss":
            raise ValueError(
                f"yieldest.py only supports dist='gauss' params (param "
                f"{p['name']!r} has dist={p.get('dist')!r}) -- the "
                f"likelihood-ratio weights and worst-case-distance beta are "
                f"Gaussian-specific; see docs/yield.md Limits")


def to_values(params, u):
    """standardized u (dict name -> value, missing == 0) -> physical values."""
    values = {}
    for p in params:
        nominal = p.get("nominal", 0.0)
        sigma = _sigma_of(p)
        values[p["name"]] = nominal + sigma * u.get(p["name"], 0.0)
    return values


def build_deck(tag, template, values):
    header = f"* yieldest.py eval {tag} (card #474)\n* sampled params:\n"
    for name, v in values.items():
        header += f"*   {name} = {v:.10g}\n"
    return header + mc.render_deck(template, values)


# --------------------------------------------------------------------- #
# Evaluator: one ngspice process per standardized point
# --------------------------------------------------------------------- #

class Evaluator:
    """Evaluates the circuit at a standardized point u, returns pass/fail
    against the declared measure/threshold. Every call is counted (honest
    `n_sims` reporting) and every deck is written under `workdir` for
    provenance (no simulation output goes in the repository -- outdir is
    always under ~/Work/data/ngspice_augmentation/yield/)."""

    def __init__(self, spec, workdir, ngspice=None, timeout=None):
        self.spec = spec
        self.workdir = workdir
        self.ngspice = ngspice or spec.get("ngspice", NGSPICE)
        self.timeout = timeout or spec.get("timeout", 20)
        self.params = spec["params"]
        self.template = spec["template"]
        self.measure = spec["measure"]
        self.fail_op = spec["fail"]["op"]
        self.fail_threshold = spec["fail"]["threshold"]
        self.n_calls = 0
        self.n_ok = 0
        self.n_not_ok = 0
        self._lock = threading.Lock()

    def eval(self, u, tag=None):
        with self._lock:
            self.n_calls += 1
            call_id = self.n_calls
        tag = tag or f"e{call_id:06d}"
        values = to_values(self.params, u)
        deck_text = build_deck(tag, self.template, values)
        path = os.path.join(self.workdir, f"{tag}.cir")
        with open(path, "w") as f:
            f.write(deck_text)
        result = mc.run_one(self.ngspice, path, self.timeout)
        ok = bool(result["ok"]) and self.measure in result.get("measures", {})
        val = result["measures"].get(self.measure) if ok else None
        fail = None
        if ok:
            if self.fail_op == "lt":
                fail = val < self.fail_threshold
            elif self.fail_op == "gt":
                fail = val > self.fail_threshold
            else:
                raise ValueError(f"unknown fail op {self.fail_op!r} (only 'lt'/'gt' supported)")
        with self._lock:
            if ok:
                self.n_ok += 1
            else:
                self.n_not_ok += 1
        return {"ok": ok, "fail": fail, "value": val, "u": dict(u), "tag": tag,
                "error": result.get("error"), "returncode": result.get("returncode")}


def _pool(max_workers):
    return concurrent.futures.ThreadPoolExecutor(max_workers=min(max_workers, 4))


# --------------------------------------------------------------------- #
# Step 1: worst-case-distance (MPP) search, gradient-free
# --------------------------------------------------------------------- #

def _grid_directions(dim_names):
    """Deterministic direction set -- no RNG, so the search itself is exactly
    reproducible run to run: coordinate unit vectors for every dimension,
    plus (2 dims) a dense 15-degree angular grid (24 directions), plus
    (3+ dims) the all-sign unit-diagonal directions (up to 32), which a
    coordinate-only search would miss entirely."""
    d = len(dim_names)
    dirs = []
    for name in dim_names:
        for s in (1.0, -1.0):
            u = {n: 0.0 for n in dim_names}
            u[name] = s
            dirs.append(u)
    if d == 2:
        a, b = dim_names
        n_dir = 24
        for k in range(n_dir):
            theta = 2 * math.pi * k / n_dir
            dirs.append({a: math.cos(theta), b: math.sin(theta)})
    elif d >= 3:
        cap = 32
        for signs in itertools.islice(itertools.product((1.0, -1.0), repeat=d), cap):
            norm = math.sqrt(d)
            dirs.append({n: s / norm for n, s in zip(dim_names, signs)})
    seen = set()
    uniq = []
    for u in dirs:
        key = tuple(round(u.get(n, 0.0), 9) for n in dim_names)
        if key not in seen:
            seen.add(key)
            uniq.append(u)
    return uniq


def _find_crossing(evaluator, direction, dim_names, r_max=15.0, expand_factor=1.7, bisect_iters=14):
    """Find the smallest r>=0 at which u=r*direction fails, assuming (the
    standard WCD/FORM assumption) the ray crosses pass->fail at most once in
    [0, r_max]. A non-ok (non-converging) ngspice point aborts that ray
    conservatively rather than guessing past it."""
    n_evals = 0
    r_lo, r_hi = 0.0, None
    r = 1.0
    while r <= r_max:
        u = {n: r * direction.get(n, 0.0) for n in dim_names}
        res = evaluator.eval(u)
        n_evals += 1
        if not res["ok"]:
            return {"found": False, "n_evals": n_evals, "reason": "ngspice-not-ok"}
        if res["fail"]:
            r_hi = r
            break
        r_lo = r
        r *= expand_factor
    if r_hi is None:
        return {"found": False, "n_evals": n_evals, "reason": "no-crossing-within-r_max"}
    for _ in range(bisect_iters):
        r_mid = 0.5 * (r_lo + r_hi)
        u = {n: r_mid * direction.get(n, 0.0) for n in dim_names}
        res = evaluator.eval(u)
        n_evals += 1
        if not res["ok"]:
            r_hi = r_mid  # conservative: treat non-convergence as inside the failing side
            continue
        if res["fail"]:
            r_hi = r_mid
        else:
            r_lo = r_mid
    return {"found": True, "r": r_hi, "n_evals": n_evals}


def _bisect_around(evaluator, direction, dim_names, r_guess, bracket_frac=0.5,
                    bisect_iters=12, expand_factor=1.5, r_max=30.0):
    """Like `_find_crossing`, but starts its bracket near `r_guess` -- far
    fewer evaluations. Used only by the local refinement step below, where
    the previous iterate's beta is already a good starting guess."""
    n_evals = 0
    r_lo = max(0.0, r_guess * (1 - bracket_frac))
    r_hi = r_guess * (1 + bracket_frac)
    u_lo = {n: r_lo * direction.get(n, 0.0) for n in dim_names}
    res_lo = evaluator.eval(u_lo)
    n_evals += 1
    if res_lo["ok"] and res_lo["fail"]:
        full = _find_crossing(evaluator, direction, dim_names, r_max=r_hi)
        full["n_evals"] += n_evals
        return full
    u_hi = {n: r_hi * direction.get(n, 0.0) for n in dim_names}
    res_hi = evaluator.eval(u_hi)
    n_evals += 1
    tries = 0
    while (not res_hi["ok"] or not res_hi["fail"]) and tries < 6 and r_hi < r_max:
        r_hi *= expand_factor
        u_hi = {n: r_hi * direction.get(n, 0.0) for n in dim_names}
        res_hi = evaluator.eval(u_hi)
        n_evals += 1
        tries += 1
    if not res_hi["ok"] or not res_hi["fail"]:
        return {"found": False, "n_evals": n_evals, "reason": "no-crossing-near-guess"}
    for _ in range(bisect_iters):
        r_mid = 0.5 * (r_lo + r_hi)
        u_mid = {n: r_mid * direction.get(n, 0.0) for n in dim_names}
        res = evaluator.eval(u_mid)
        n_evals += 1
        if not res["ok"]:
            r_hi = r_mid
            continue
        if res["fail"]:
            r_hi = r_mid
        else:
            r_lo = r_mid
    return {"found": True, "r": r_hi, "n_evals": n_evals}


def search_mpp(evaluator, dim_names, max_workers=4, refine_iters=6, refine_step0=0.35):
    """Worst-case-distance search: deterministic direction grid (parallel
    across <=4 ngspice processes) + gradient-free directional pattern-search
    refinement (Hooke-Jeeves-style on the unit direction vector)."""
    directions = _grid_directions(dim_names)
    results = []
    with _pool(max_workers) as pool:
        futs = {pool.submit(_find_crossing, evaluator, d, dim_names): d for d in directions}
        for fut in concurrent.futures.as_completed(futs):
            r = fut.result()
            r["direction"] = futs[fut]
            results.append(r)
    n_evals = sum(r["n_evals"] for r in results)
    found = [r for r in results if r["found"]]
    if not found:
        return {"ok": False, "reason": "no direction found a failure within r_max",
                "n_directions_tried": len(directions), "n_evals": n_evals}
    best = min(found, key=lambda r: r["r"])
    direction = dict(best["direction"])
    beta = best["r"]

    step = refine_step0
    for _ in range(refine_iters):
        candidates = []
        for name in dim_names:
            for s in (step, -step):
                cand = dict(direction)
                cand[name] = cand.get(name, 0.0) + s
                norm = math.sqrt(sum(v * v for v in cand.values()))
                if norm > 0:
                    candidates.append({n: v / norm for n, v in cand.items()})
        with _pool(max_workers) as pool:
            futs = {pool.submit(_bisect_around, evaluator, c, dim_names, beta): c for c in candidates}
            crossed = []
            for fut in concurrent.futures.as_completed(futs):
                cross = fut.result()
                n_evals += cross["n_evals"]
                cross["direction"] = futs[fut]
                crossed.append(cross)
        improving = [c for c in crossed if c["found"] and c["r"] < beta - 1e-9]
        if improving:
            b = min(improving, key=lambda c: c["r"])
            beta = b["r"]
            direction = b["direction"]
        else:
            step *= 0.5
        if step < 1e-3:
            break

    u_star = {n: beta * direction.get(n, 0.0) for n in dim_names}
    return {"ok": True, "beta": beta, "direction": direction, "u_star": u_star,
            "n_directions_tried": len(directions), "n_evals": n_evals}


# --------------------------------------------------------------------- #
# Step 2a: mean-shift importance sampling
# --------------------------------------------------------------------- #

def importance_sampling(evaluator, dim_names, u_star, n, seed, max_workers=4):
    rng = random.Random(seed)
    samples = [{name: u_star.get(name, 0.0) + rng.gauss(0.0, 1.0) for name in dim_names}
               for _ in range(n)]
    recs = [None] * n
    with _pool(max_workers) as pool:
        futs = {pool.submit(evaluator.eval, u, f"is{i:05d}"): i for i, u in enumerate(samples)}
        for fut in concurrent.futures.as_completed(futs):
            recs[futs[fut]] = fut.result()

    beta2 = sum(v * v for v in u_star.values())
    weighted = []  # (weight, fail) for every OK sample
    for r in recs:
        if not r["ok"]:
            continue
        dot = sum(r["u"].get(name, 0.0) * u_star.get(name, 0.0) for name in dim_names)
        w = math.exp(-dot + 0.5 * beta2)
        weighted.append((w, r["fail"]))

    n_ok = len(weighted)
    if n_ok == 0:
        return {"ok": False, "reason": "no ok runs", "n": n, "n_ok": 0}

    contrib = [w if fail else 0.0 for w, fail in weighted]
    pf_hat = sum(contrib) / n_ok
    clipped = pf_hat > 1.0
    pf_hat = min(pf_hat, 1.0)
    n_fail = sum(1 for _, fail in weighted if fail)

    sum_w = sum(w for w, _ in weighted)
    sum_w2 = sum(w * w for w, _ in weighted)
    ess = (sum_w * sum_w / sum_w2) if sum_w2 > 0 else 0.0

    if n_ok >= 2 and any(c != contrib[0] for c in contrib):
        se = math.sqrt(statistics.variance(contrib) / n_ok)
    else:
        se = 0.0

    if pf_hat > 0 and se > 0:
        cv = se / pf_hat
        lo = pf_hat * math.exp(-Z95 * cv)
        hi = min(1.0, pf_hat * math.exp(Z95 * cv))
    elif pf_hat > 0:
        lo = hi = pf_hat
    else:
        # no weighted failures at all -- report a rule-of-three-style upper
        # bound rather than fabricating a point estimate
        lo, hi = 0.0, 3.0 / n_ok

    equiv_sigma = -NormalDist().inv_cdf(pf_hat) if 0 < pf_hat < 1 else None

    out = {
        "ok": True, "method": "mean_shift_importance_sampling", "n": n, "n_ok": n_ok,
        "n_fail_weighted_contrib": n_fail, "pf": pf_hat, "pf_ci95": [lo, hi],
        "equivalent_sigma": equiv_sigma, "ess": ess, "u_star": dict(u_star),
        # the UNweighted fail rate under the shifted sampling distribution --
        # not a valid Pf estimator, exposed only as a diagnostic/red-proof
        # tool (see docs/yield.md and tests/test_yield.py's dropped-weight
        # red-proof, which uses exactly this field).
        "naive_unweighted_fail_rate_NOT_A_PF_ESTIMATE": n_fail / n_ok,
    }
    if clipped:
        out["note"] = "pf_hat clipped from >1 -- high importance-weight variance at this n; increase n_is"
    return out


# --------------------------------------------------------------------- #
# Step 2b: scaled-sigma sampling (SSS) extrapolation -- independent cross-check
# --------------------------------------------------------------------- #

def default_sss_scales(beta):
    """Scale factors chosen (and checked, see docs/yield.md) so each scale's
    Pf_c lands in a moderate, cheaply-measurable range (roughly 1%-60%)
    instead of being tiny (needs huge N) or saturated (needs huge N the
    other way) -- spanning below and up to beta itself. Clipped at 1.2 so a
    "scale" is never actually a shrink."""
    factors = (0.35, 0.5, 0.65, 0.8, 1.0, 1.2)
    return sorted({round(max(1.2, beta * f), 4) for f in factors})


def scaled_sigma_sampling(evaluator, dim_names, n_per_scale, scales, seed, max_workers=4):
    """Independent cross-check (no WCD search, no importance weights): plain
    MC at several enlarged standard deviations c, giving a directly-measured
    Pf_c at each; beta_hat = -c*Phi^-1(Pf_c) recovers the worst-case
    distance under the same "boundary is approximately a hyperplane at
    distance beta" assumption WCD/FORM already makes. The per-scale
    estimates are combined by WEIGHTED least squares (regressing
    z_c=Phi^-1(Pf_c) against 1/c through the origin, slope=-beta), weight =
    1/Var(z_c) via the delta method on a binomial proportion -- this
    naturally downweights a scale whose failure count was too small to be
    reliable, rather than needing a separate escalation/retry heuristic.
    See docs/yield.md for why simple unweighted averaging across scales was
    tried first and found noisier at a matched run budget."""
    scale_results = {}
    for c in scales:
        rng = random.Random(_seed_int(seed, "sss", c))
        samples = [{name: c * rng.gauss(0.0, 1.0) for name in dim_names}
                   for _ in range(n_per_scale)]
        recs = [None] * n_per_scale
        with _pool(max_workers) as pool:
            futs = {pool.submit(evaluator.eval, u, f"sss{c:g}_{i:05d}"): i
                    for i, u in enumerate(samples)}
            for fut in concurrent.futures.as_completed(futs):
                recs[futs[fut]] = fut.result()
        ok_recs = [r for r in recs if r["ok"]]
        n_ok = len(ok_recs)
        k = sum(1 for r in ok_recs if r["fail"])
        scale_results[c] = {"n": n_per_scale, "n_ok": n_ok, "k": k,
                             "pf_c": (k / n_ok) if n_ok else None}

    nd = NormalDist()
    num, den = 0.0, 0.0
    used = []
    for c, r in scale_results.items():
        pf_c, n_ok, k = r["pf_c"], r["n_ok"], r["k"]
        if pf_c is None or k == 0 or k == n_ok:
            continue  # 0/n or n/n: Phi^-1 undefined / infinite variance -- excluded, not guessed
        z = nd.inv_cdf(pf_c)
        phi_z = math.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)
        var_z = (pf_c * (1 - pf_c)) / (n_ok * phi_z * phi_z)
        if var_z <= 0:
            continue
        w = 1.0 / var_z
        x = 1.0 / c
        r["beta_est"] = -c * z
        r["weight"] = w
        num += w * x * z
        den += w * x * x
        used.append(c)

    n_total = sum(r["n"] for r in scale_results.values())
    if den == 0:
        return {"ok": False, "reason": "no scale gave a usable (0<k<n) failure count",
                "scales": scale_results, "n_total": n_total}

    beta_hat = -num / den
    se_beta = math.sqrt(1.0 / den)
    pf_sss = NormalDist().cdf(-beta_hat)
    lo_beta, hi_beta = beta_hat - Z95 * se_beta, beta_hat + Z95 * se_beta
    pf_lo, pf_hi = NormalDist().cdf(-hi_beta), NormalDist().cdf(-lo_beta)
    return {
        "ok": True, "method": "scaled_sigma_sampling", "pf": pf_sss, "pf_ci95": [pf_lo, pf_hi],
        "equivalent_sigma": beta_hat, "beta_mean": beta_hat, "beta_se": se_beta,
        "scales_used": used, "scales": scale_results, "n_total": n_total,
    }


# --------------------------------------------------------------------- #
# Plain Monte Carlo (sanity check at low sigma, where it is affordable)
# --------------------------------------------------------------------- #

def plain_mc(evaluator, dim_names, n, seed, max_workers=4):
    rng = random.Random(seed)
    samples = [{name: rng.gauss(0.0, 1.0) for name in dim_names} for _ in range(n)]
    recs = [None] * n
    with _pool(max_workers) as pool:
        futs = {pool.submit(evaluator.eval, u, f"pmc{i:05d}"): i for i, u in enumerate(samples)}
        for fut in concurrent.futures.as_completed(futs):
            recs[futs[fut]] = fut.result()
    ok_recs = [r for r in recs if r["ok"]]
    n_ok = len(ok_recs)
    if n_ok == 0:
        return {"ok": False, "reason": "no ok runs", "n": n}
    k = sum(1 for r in ok_recs if r["fail"])
    p = k / n_ok
    se = math.sqrt(max(p * (1 - p), 0.0) / n_ok)
    lo, hi = max(0.0, p - Z95 * se), min(1.0, p + Z95 * se)
    return {"ok": True, "n": n, "n_ok": n_ok, "k": k, "pf": p, "pf_ci95": [lo, hi]}


# --------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------- #

def estimate_yield(spec, outdir, progress=None):
    spec = dict(spec)
    if "template" not in spec:
        with open(spec["template_path"]) as f:
            spec["template"] = f.read()
    check_params_supported(spec["params"])
    dim_names = [p["name"] for p in spec["params"]]
    max_workers = min(spec.get("max_workers", 4), 4)

    decks_dir = os.path.join(outdir, "decks")
    os.makedirs(decks_dir, exist_ok=True)
    evaluator = Evaluator(spec, decks_dir)

    nominal = evaluator.eval({n: 0.0 for n in dim_names}, tag="nominal")
    result = {"ok": True, "seed": spec.get("seed", 0), "params": spec["params"],
              "measure": spec["measure"], "fail": spec["fail"], "nominal_check": nominal}
    if not nominal["ok"]:
        result.update(ok=False, stage="nominal", reason="nominal point did not converge")
        result["n_sims"] = evaluator.n_calls
        _write_manifest(result, outdir)
        return result
    if nominal["fail"]:
        result.update(ok=False, stage="nominal",
                       reason="nominal (0-sigma) point already fails the given threshold -- "
                              "not a yield problem in the usual sense")
        result["n_sims"] = evaluator.n_calls
        _write_manifest(result, outdir)
        return result

    mpp = search_mpp(evaluator, dim_names, max_workers=max_workers,
                      refine_iters=spec.get("refine_iters", 6))
    result["mpp_search"] = mpp
    if not mpp["ok"]:
        result.update(ok=False,
                       reason="worst-case-distance search found no failure within range -- "
                              "cannot build an importance-sampling shift; see docs/yield.md Limits")
        result["n_sims"] = evaluator.n_calls
        _write_manifest(result, outdir)
        return result

    n_is = spec.get("n_is", 300)
    is_res = importance_sampling(evaluator, dim_names, mpp["u_star"], n_is,
                                  seed=_seed_int(spec.get("seed", 0), "is"),
                                  max_workers=max_workers)
    result["importance_sampling"] = is_res

    # SSS's scale range is chosen from the WCD search's own beta (see
    # default_sss_scales) -- a fixed scale list would put Pf_c in a
    # nonsensical range (too tiny, or saturated) for a target far from what
    # it was tuned for; an explicit spec["sss_scales"] overrides this.
    sss_scales = spec.get("sss_scales") or default_sss_scales(mpp["beta"])
    n_sss = spec.get("n_sss", 3000)
    sss_res = scaled_sigma_sampling(evaluator, dim_names, n_sss, sss_scales,
                                     seed=_seed_int(spec.get("seed", 0), "sss"),
                                     max_workers=max_workers)
    result["scaled_sigma_sampling"] = sss_res

    disagreement = None
    if is_res.get("ok") and sss_res.get("ok") and is_res.get("pf", 0) > 0:
        disagreement = sss_res["pf"] / is_res["pf"]
        if disagreement > 3.0 or disagreement < 1.0 / 3.0:
            result["limitation_warning"] = (
                f"importance-sampling Pf ({is_res['pf']:.3e}) and scaled-sigma-sampling Pf "
                f"({sss_res['pf']:.3e}) disagree by {disagreement:.1f}x. Two known causes (see "
                f"docs/yield.md Limits, which this tool cannot itself tell apart): (1) a second "
                f"failure region the worst-case-distance search's direction set did not find, "
                f"or (2) SSS's hyperplane-boundary assumption breaking down on a strongly curved "
                f"failure boundary even with a single region (importance sampling does not share "
                f"this assumption and is the more robust of the two in that case). Do not trust "
                f"either point estimate alone here."
            )
    result["method_disagreement_ratio"] = disagreement

    result["pf_estimate"] = is_res.get("pf") if is_res.get("ok") else sss_res.get("pf")
    result["pf_ci95"] = is_res.get("pf_ci95") if is_res.get("ok") else sss_res.get("pf_ci95")
    result["equivalent_sigma"] = (is_res.get("equivalent_sigma") if is_res.get("ok")
                                   else sss_res.get("equivalent_sigma"))
    beta = mpp["beta"] or 1.0
    result["most_probable_failure_point"] = {
        "u_star": mpp["u_star"], "beta": mpp["beta"],
        "values": to_values(spec["params"], mpp["u_star"]),
        "contribution_fraction": {n: (mpp["u_star"].get(n, 0.0) ** 2) / (beta * beta)
                                   for n in dim_names},
    }
    result["n_sims"] = evaluator.n_calls
    result["n_sims_breakdown"] = {
        "nominal_check": 1, "mpp_search": mpp.get("n_evals", 0),
        "importance_sampling": is_res.get("n", 0),
        "scaled_sigma_sampling": sss_res.get("n_total", 0),
    }
    _write_manifest(result, outdir)
    return result


def _write_manifest(result, outdir):
    path = os.path.join(outdir, "manifest.json")
    with open(path, "w") as f:
        json.dump(result, f, indent=1, default=str)
    return path


# --------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------- #

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="estimate a failure probability from a spec JSON file")
    r.add_argument("spec")
    r.add_argument("--out", default=None)

    args = ap.parse_args(argv)
    if args.cmd == "run":
        with open(args.spec) as f:
            spec = json.load(f)
        outdir = args.out or os.path.join(DATA_ROOT, f"run_{int(time.time())}")
        os.makedirs(outdir, exist_ok=True)
        result = estimate_yield(spec, outdir)
        printable = {k: v for k, v in result.items() if k != "mpp_search"}
        print(json.dumps(printable, indent=1, default=str))
        print(f"manifest: {os.path.join(outdir, 'manifest.json')}", file=sys.stderr)
        return 0 if result.get("ok") else 2


if __name__ == "__main__":
    sys.exit(main())
