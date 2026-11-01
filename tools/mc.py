#!/usr/bin/env python3
"""mc.py -- Monte Carlo runner with correct global (process) vs local
(mismatch) variation (card #473, plan T10.5).

## Why this exists (what plain ngspice does today, with evidence)

ngspice already has the random primitives (`agauss`/`gauss`/`aunif`/`unif`
inside `.param`, `setseed`, `mc_source`, the `.control` `dowhile`/`alter`
loop idiom shown in `ngspice-src/examples/Monte_Carlo/*.sp`). What it does
NOT have is a way to make a *shared* `.param` behave like a single shared
random draw. Probed directly on the pinned system ngspice-45.2 (reproducer
below, and see docs/mc.md for the transcript):

    .param proc_off = agauss(0, 10m, 3)
    R1 1 0 R='1k*(1+proc_off)'
    R2 2 0 R='1k*(1+proc_off)'

Both R1 and R2 reference the *same* `.param` by name. If `.param` bound a
value once at declaration (the behaviour anyone would assume from the
manual's wording and from every other SPICE dialect's `.param`), R1 and R2
would always read back equal. They do not -- on a single parse, with no
`.control` loop at all:

    @r1[resistance] = 1.010073e+03
    @r2[resistance] = 1.000933e+03

Indirecting through two more `.param`s that are plain copies (`a = proc_off`,
`b = proc_off`, no `agauss` retyped) still does not converge them. ngspice's
numparam treats a `.param` whose RHS contains a random function as macro
text that is re-substituted and re-evaluated at *every* point of use, not a
value bound once. This is a real, reproducible defect for Monte Carlo:
the natural idiom for "shared/global process variation" -- one `.param`,
referenced by several devices -- silently gives *independent*, not
*correlated*, draws. (`setseed`+`mc_source` inside a `.control dowhile` loop
does not fix this either -- it re-parses the whole deck each iteration, so
the same per-reference re-evaluation happens every run.) This matches
research/feature-survey.md item #7's flag ("no distinct local-mismatch
primitive separate from a plain global-parameter MC") and goes a step
further: even the *plain global* case is broken by this re-evaluation
behaviour, not just "missing a mismatch primitive".

A second, independently-reported pitfall (in this fork's own copy of
`examples/Monte_Carlo/MC_2_control.sp`, a comment in the file itself):
looping `tran`/`op` inside one ngspice process via `mc_source` without a
`reset` shows "some strange drift towards lower and lower frequencies" --
i.e. the in-process reload path is not trusted even by its own example.

Given both, this tool does **not** rely on ngspice's `agauss`/`.param`
machinery for the actual random draws at all -- it treats that as broken
for the "global" case and awkward for "local" (no way to say "N independent
per-instance draws, Pelgrom-scaled" in one `.param`). Instead:

  - all sampling happens here, in Python, with an explicit, documented,
    hash-derived seed per (run, parameter) -- reproducible without needing
    any other run's state (so replaying run k alone reproduces run k alone);
  - each Monte Carlo point is a *separate* ngspice process on a deck with
    the sampled values substituted in as plain numeric literals (via
    `@@TOKEN@@` placeholders -- never `.param ... agauss(...)`), so no
    ngspice numparam behaviour (correct or broken) is depended on;
  - "global" vs "local" is purely a matter of whether two devices' tokens
    share a name (global: same token, same per-run value everywhere) or
    have distinct names (local: independent per-run, per-instance value,
    optionally Pelgrom-scaled by 1/sqrt(W*L)) -- the same convention the
    sky130 PDK's public documentation describes for its own
    `mc_mm_switch` (mismatch, per-instance)/`mc_pr_switch` (process,
    shared) idiom, but enforced correctly here rather than through the
    `.param`-sharing idiom just shown to be broken.

This is the "tools-level solution" the card's brief prefers when the core
already has the primitives it needs (a deterministic RNG, subprocess
ngspice, and plain numeric substitution) -- no C change was made or was
found necessary. See the bottom of this file / docs/mc.md for what was
checked and what a native-side fix would look like, left for a future card.

## Design

A "spec" (plain dict, or a JSON file with the same shape) declares:

  - `template`: deck text (or `template_path`: a file) containing
    `@@NAME@@` placeholders -- never inside a SPICE string in a way that
    collides with `{}`/`'...'`/`$` (ngspice's own expression/comment
    syntax), so `@@...@@` was chosen deliberately: ngspice never emits or
    parses that token shape itself.
  - `seed`: master seed (int).
  - `n`: number of Monte Carlo points.
  - `params`: list of dicts, each:
      `name`   -- placeholder name (matches `@@name@@` in the template,
                  case-sensitive).
      `scope`  -- "global" (this exact param dict is normally referenced,
                  by using the SAME `name`, in more than one place in the
                  template -- one shared draw per run) or "local" (an
                  independent per-run draw; give each mismatched instance
                  its own `params` entry with a different `name`).
      `dist`   -- "gauss" (default) or "unif".
      `nominal`-- centre value (added to the sampled deviation).
      `sigma`  -- 1-sigma (gauss) or half-range (unif) deviation, OR:
      `pelgrom`-- {"coeff": A, "w": W, "l": L} -> sigma = A / sqrt(W*L)
                  (Pelgrom's law; W, L in whatever consistent units A was
                  characterised in, e.g. um and V*um).
  - `measures`: list of names the deck's own `.control` block must report
      back via a line of the exact form `MCMEAS <name>=<value>` on stdout
      (see docs/mc.md for the idiom: `let x=...; set xs=$&x;
      echo MCMEAS name=$xs`). This is how "sampled parameters are read
      back per run ... via .meas of the device values or currents" is
      satisfied without depending on ngspice's own `.meas` (which pairs
      awkwardly with a bare `.op`) or a rawfile parser.
  - `checks` (optional): list of assertions evaluated against the 99% CI
      of the resulting sample statistics -- see `evaluate_checks`.
  - `ngspice`, `max_workers` (<=4, enforced), `timeout` (seconds/run),
    `omp_threads` ("auto" default = clamp(cores // max_workers, 1, 8), N, or 0
    to keep the build's spinit value; card #527, tools/ompthreads.py).

Runner: `run_campaign(spec, outdir)` draws N independent samples (one
Python `random.Random` stream per (run, param-name), see `_rng_for`),
writes one deck per run under `outdir/decks/`, runs them with bounded
concurrency (<=4 concurrent ngspice processes -- a `ThreadPoolExecutor` of
that size; each call is a `subprocess.run`, so the 4-process cap is a real
OS-process cap, not just a thread cap), collects the `MCMEAS` lines,
computes statistics, evaluates any `checks`, and writes `manifest.json`
(full provenance: every per-run sampled value and its own seed derivation,
so ANY single run can be replayed with `replay_run` with no dependency on
having run the rest of the campaign).

## Statistics (stdlib only -- no numpy/scipy on this box, confirmed)

- mean: `statistics.mean`, sample stdev: `statistics.stdev` (ddof=1).
- 99% CI for a mean: normal approximation via `statistics.NormalDist`
  (z_0.995 = 2.5758). At the N ~= 500 this tool's validation uses, the
  Student-t correction (df=499) changes z by < 0.4% -- normal approx is
  used and this approximation is stated here, not hidden.
- 99% CI for a sample stdev: chi-square quantiles via the Wilson-Hilferty
  cube-root normal approximation to the chi-square distribution (a
  standard, well-known approximation; accurate to a fraction of a percent
  at N ~= 500 degrees of freedom, which is the only regime this tool
  currently uses it in -- not claimed accurate at small N).
- 99% CI for a Pearson correlation: Fisher z-transform, standard-normal CI
  in z-space, inverse-transformed back (asymptotically normal, standard
  practice, accurate at N ~= 500).
"""
import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import random
import re
import statistics
import subprocess
import sys
import time
from statistics import NormalDist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ompthreads  # noqa: E402  (card #527)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/mc")

TOKEN_RE = re.compile(r"@@(\w+)@@")
MCMEAS_RE = re.compile(r"^MCMEAS\s+(\w+)\s*=\s*([-+0-9.eE]+)\s*$")

Z995 = NormalDist().inv_cdf(0.995)  # 99% two-sided normal quantile, ~2.5758


# --------------------------------------------------------------------- #
# Seeding / sampling
# --------------------------------------------------------------------- #

def _seed_int(*parts):
    """Deterministic 64-bit seed from an arbitrary tuple of identifiers.

    Hash-derived (not "master_seed + run_idx" arithmetic) so every
    (run, param) stream is independent of every other -- no carried state,
    no risk of accidental correlation between nearby run indices or
    between two differently-named params that happen to add up the same.
    """
    h = hashlib.sha256("|".join(str(p) for p in parts).encode()).digest()
    return int.from_bytes(h[:8], "big")


def _rng_for(seed, run_idx, param_name):
    return random.Random(_seed_int(seed, run_idx, param_name))


def _sigma_of(param):
    if "pelgrom" in param:
        p = param["pelgrom"]
        return p["coeff"] / math.sqrt(p["w"] * p["l"])
    return param["sigma"]


def sample_param(param, seed, run_idx):
    rng = _rng_for(seed, run_idx, param["name"])
    nominal = param.get("nominal", 0.0)
    sigma = _sigma_of(param)
    dist = param.get("dist", "gauss")
    if dist == "gauss":
        return nominal + sigma * rng.gauss(0.0, 1.0)
    elif dist == "unif":
        return nominal + sigma * rng.uniform(-1.0, 1.0)
    else:
        raise ValueError(f"unknown dist {dist!r} for param {param['name']!r}")


# --------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------- #

def mean_ci99(values):
    n = len(values)
    m = statistics.mean(values)
    if n < 2:
        return m, m, m
    s = statistics.stdev(values)
    half = Z995 * s / math.sqrt(n)
    return m, m - half, m + half


def _chi2_quantile(p, k):
    """Wilson-Hilferty cube-root approximation to the chi-square quantile
    with k degrees of freedom at cumulative probability p. See module
    docstring for accuracy caveat (large-k regime only, used here at
    k ~= N-1 ~= 499)."""
    z = NormalDist().inv_cdf(p)
    val = k * (1.0 - 2.0 / (9.0 * k) + z * math.sqrt(2.0 / (9.0 * k))) ** 3
    return max(val, 1e-12)


def sigma_ci99(values):
    n = len(values)
    s = statistics.stdev(values)
    if n < 3:
        return s, s, s
    k = n - 1
    lo_chi2 = _chi2_quantile(0.005, k)
    hi_chi2 = _chi2_quantile(0.995, k)
    lo = s * math.sqrt(k / hi_chi2)
    hi = s * math.sqrt(k / lo_chi2)
    return s, lo, hi


def pearson(xs, ys):
    n = len(xs)
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    denom = math.sqrt(sxx * syy)
    if denom == 0:
        return 0.0
    r = sxy / denom
    return max(-1.0, min(1.0, r))


def corr_ci99(xs, ys):
    n = len(xs)
    r = pearson(xs, ys)
    if n < 4 or abs(r) >= 1.0 - 1e-12:
        return r, r, r
    z = math.atanh(r)
    se = 1.0 / math.sqrt(n - 3)
    half = Z995 * se
    lo = math.tanh(z - half)
    hi = math.tanh(z + half)
    return r, lo, hi


# --------------------------------------------------------------------- #
# Deck generation / running
# --------------------------------------------------------------------- #

def render_deck(template, values):
    def sub(m):
        name = m.group(1)
        if name not in values:
            raise KeyError(f"template placeholder @@{name}@@ has no sampled value")
        return format(values[name], ".10g")
    return TOKEN_RE.sub(sub, template)


def _parse_mcmeas(stdout_text):
    out = {}
    for line in stdout_text.splitlines():
        m = MCMEAS_RE.match(line.strip())
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def run_one(ngspice, deck_path, timeout, extra_env=None):
    t0 = time.time()
    try:
        proc = subprocess.run(
            [ngspice, "-b", deck_path],
            capture_output=True, text=True, timeout=timeout,
            env=dict(os.environ, **extra_env) if extra_env else None)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "timeout", "wall_s": time.time() - t0}
    except FileNotFoundError as e:
        return {"ok": False, "error": f"ngspice not found: {e}", "wall_s": time.time() - t0}
    wall = time.time() - t0
    measures = _parse_mcmeas(proc.stdout)
    ok = (proc.returncode == 0)
    return {
        "ok": ok, "returncode": proc.returncode, "measures": measures,
        "wall_s": wall, "stdout_tail": proc.stdout[-2000:], "stderr_tail": proc.stderr[-2000:],
    }


def sample_run_params(spec, run_idx):
    """Sample every declared param for one run index. Depends only on
    spec['seed'], run_idx, and each param's own name -- no other run's
    state -- so this is exactly what a lone replay recomputes."""
    values = {}
    for p in spec["params"]:
        values[p["name"]] = sample_param(p, spec["seed"], run_idx)
    return values


def build_run_deck(spec, run_idx, values):
    header = (
        f"* mc.py run {run_idx} (card #473) -- master seed {spec['seed']}\n"
        f"* sampled params:\n"
    )
    for p in spec["params"]:
        header += f"*   {p['name']} ({p.get('scope', '?')}, {p.get('dist', 'gauss')}) = {values[p['name']]:.10g}\n"
    return header + render_deck(spec["template"], values)


def run_campaign(spec, outdir, progress=None):
    spec = dict(spec)
    if "template" not in spec:
        with open(spec["template_path"]) as f:
            spec["template"] = f.read()
    n = spec["n"]
    max_workers = min(spec.get("max_workers", 4), 4)
    timeout = spec.get("timeout", 30)
    ngspice = spec.get("ngspice", NGSPICE)
    # card #527: OpenMP threads per process = clamp(cores // max_workers, 1, 8)
    # ("omp_threads": "auto" default | N | 0 = leave the build's spinit value);
    # a no-op for max_workers <= cores // 8 (<= 4 workers on 32 cores).
    omp_n = ompthreads.resolve_setting(spec.get("omp_threads", "auto"), max_workers) if max_workers > 1 else None
    omp_env = ompthreads.spinit_env(ngspice, omp_n, outdir) or None

    decks_dir = os.path.join(outdir, "decks")
    logs_dir = os.path.join(outdir, "logs")
    os.makedirs(decks_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    runs = [None] * n

    def _do(i):
        values = sample_run_params(spec, i)
        deck_text = build_run_deck(spec, i, values)
        deck_path = os.path.join(decks_dir, f"run_{i:05d}.cir")
        with open(deck_path, "w") as f:
            f.write(deck_text)
        result = run_one(ngspice, deck_path, timeout, omp_env)
        log_path = os.path.join(logs_dir, f"run_{i:05d}.log")
        with open(log_path, "w") as f:
            f.write(result.get("stdout_tail", "") + "\n---stderr---\n" + result.get("stderr_tail", ""))
        record = {"idx": i, "params": values, "deck": deck_path, "log": log_path}
        record.update({k: v for k, v in result.items() if k not in ("stdout_tail", "stderr_tail")})
        return record

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        futs = {pool.submit(_do, i): i for i in range(n)}
        done = 0
        for fut in concurrent.futures.as_completed(futs):
            i = futs[fut]
            runs[i] = fut.result()
            done += 1
            if progress:
                progress(done, n)

    n_ok = sum(1 for r in runs if r["ok"])
    stats = compute_stats(spec, runs)
    checks = evaluate_checks(spec.get("checks", []), stats)

    manifest = {
        "seed": spec["seed"], "n": n, "n_ok": n_ok,
        "omp_threads_applied": omp_n if omp_env else None,
        "params": spec["params"], "measures": spec.get("measures", []),
        "runs": runs, "stats": stats, "checks": checks,
    }
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    return manifest


def compute_stats(spec, runs):
    ok_runs = [r for r in runs if r["ok"]]
    stats = {"n_ok": len(ok_runs), "n_total": len(runs), "params": {}, "measures": {}}

    for pname in [p["name"] for p in spec["params"]]:
        vals = [r["params"][pname] for r in ok_runs]
        if len(vals) >= 2:
            m, mlo, mhi = mean_ci99(vals)
            s, slo, shi = sigma_ci99(vals)
            stats["params"][pname] = {
                "mean": m, "mean_ci99": [mlo, mhi],
                "sigma": s, "sigma_ci99": [slo, shi], "n": len(vals),
            }

    for mname in spec.get("measures", []):
        vals = [r["measures"][mname] for r in ok_runs if mname in r["measures"]]
        if len(vals) >= 2:
            m, mlo, mhi = mean_ci99(vals)
            s, slo, shi = sigma_ci99(vals)
            stats["measures"][mname] = {
                "mean": m, "mean_ci99": [mlo, mhi],
                "sigma": s, "sigma_ci99": [slo, shi], "n": len(vals),
            }

    # pairwise correlations between every two declared measures (small N of
    # measures expected -- a matched pair or two -- so no combinatorial
    # blow-up)
    mnames = spec.get("measures", [])
    stats["correlations"] = {}
    for i, a in enumerate(mnames):
        for b in mnames[i + 1:]:
            xs, ys = [], []
            for r in ok_runs:
                if a in r["measures"] and b in r["measures"]:
                    xs.append(r["measures"][a])
                    ys.append(r["measures"][b])
            if len(xs) >= 4:
                r_, lo, hi = corr_ci99(xs, ys)
                stats["correlations"][f"{a}:{b}"] = {"r": r_, "ci99": [lo, hi], "n": len(xs)}
    return stats


def evaluate_checks(checks, stats):
    out = []
    for c in checks:
        ctype = c["type"]
        result = {"check": c, "state": "FAIL", "detail": ""}
        try:
            if ctype == "correlation":
                key = f"{c['a']}:{c['b']}"
                key_alt = f"{c['b']}:{c['a']}"
                entry = stats["correlations"].get(key) or stats["correlations"].get(key_alt)
                if entry is None:
                    result["state"] = "COULD-NOT-RUN"
                    result["detail"] = f"no correlation computed for {c['a']}/{c['b']}"
                else:
                    lo, hi = entry["ci99"]
                    expect = c["expect"]
                    ok = (lo - c.get("tol", 0.0)) <= expect <= (hi + c.get("tol", 0.0))
                    result["state"] = "pass" if ok else "FAIL"
                    result["detail"] = f"r={entry['r']:.4f} ci99={entry['ci99']} expect={expect} n={entry['n']}"
            elif ctype == "mean":
                entry = stats["measures"].get(c["measure"]) or stats["params"].get(c["measure"])
                if entry is None:
                    result["state"] = "COULD-NOT-RUN"
                    result["detail"] = f"no stats for {c['measure']}"
                else:
                    lo, hi = entry["mean_ci99"]
                    expect = c["expect"]
                    ok = lo <= expect <= hi
                    result["state"] = "pass" if ok else "FAIL"
                    result["detail"] = f"mean={entry['mean']:.6g} ci99={entry['mean_ci99']} expect={expect}"
            elif ctype == "sigma":
                entry = stats["measures"].get(c["measure"]) or stats["params"].get(c["measure"])
                if entry is None:
                    result["state"] = "COULD-NOT-RUN"
                    result["detail"] = f"no stats for {c['measure']}"
                else:
                    lo, hi = entry["sigma_ci99"]
                    expect = c["expect"]
                    ok = lo <= expect <= hi
                    result["state"] = "pass" if ok else "FAIL"
                    result["detail"] = f"sigma={entry['sigma']:.6g} ci99={entry['sigma_ci99']} expect={expect}"
            else:
                result["state"] = "COULD-NOT-RUN"
                result["detail"] = f"unknown check type {ctype!r}"
        except Exception as e:
            result["state"] = "COULD-NOT-RUN"
            result["detail"] = f"exception: {e}"
        out.append(result)
    return out


# --------------------------------------------------------------------- #
# Replay
# --------------------------------------------------------------------- #

def replay_run(manifest_path, run_idx, spec=None):
    """Recompute run `run_idx`'s sampled params purely from (seed, params,
    run_idx) -- NOT by reading the manifest's recorded values -- rerun it,
    and compare against the manifest's recorded measures. This is the
    actual "replaying run k alone reproduces its result" proof: it does
    not reuse anything from the original campaign except the spec.
    """
    with open(manifest_path) as f:
        manifest = json.load(f)
    if spec is None:
        raise ValueError("replay_run needs the original spec (params/template/seed) -- "
                          "the manifest stores results, not the template text, by design "
                          "(templates can be large; pass the same spec dict/file used for the run)")
    spec = dict(spec)
    if "template" not in spec:
        with open(spec["template_path"]) as f:
            spec["template"] = f.read()
    spec["seed"] = manifest["seed"]

    values = sample_run_params(spec, run_idx)
    orig = manifest["runs"][run_idx]
    for pname, v in orig["params"].items():
        if abs(values.get(pname, float("nan")) - v) > 1e-9 * max(1.0, abs(v)):
            return {"ok": False, "reason": f"param {pname} mismatch on resample: {values.get(pname)} != {v}"}

    deck_text = build_run_deck(spec, run_idx, values)
    tmp_path = orig["deck"] + ".replay"
    with open(tmp_path, "w") as f:
        f.write(deck_text)
    result = run_one(spec.get("ngspice", NGSPICE), tmp_path, spec.get("timeout", 30))
    match = result["ok"] == orig["ok"]
    if result["ok"] and orig["ok"]:
        for k, v in orig["measures"].items():
            rv = result["measures"].get(k)
            if rv is None or abs(rv - v) > 1e-9 * max(1.0, abs(v)):
                match = False
    return {"ok": match, "replayed": result, "original": orig}


# --------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------- #

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run a Monte Carlo campaign from a spec JSON file")
    r.add_argument("spec", help="path to a spec JSON file")
    r.add_argument("--out", default=None)

    p = sub.add_parser("replay", help="replay a single run from a previous campaign")
    p.add_argument("manifest")
    p.add_argument("--run", type=int, required=True)
    p.add_argument("--spec", required=True)

    args = ap.parse_args(argv)

    if args.cmd == "run":
        with open(args.spec) as f:
            spec = json.load(f)
        outdir = args.out or os.path.join(DATA_ROOT, f"run_{int(time.time())}")
        os.makedirs(outdir, exist_ok=True)

        def progress(done, total):
            print(f"\r{done}/{total} runs", end="", file=sys.stderr)

        manifest = run_campaign(spec, outdir, progress=progress)
        print(file=sys.stderr)
        print(json.dumps({"n_ok": manifest["n_ok"], "n": manifest["n"],
                           "stats": manifest["stats"], "checks": manifest["checks"]}, indent=1))
        print(f"manifest: {os.path.join(outdir, 'manifest.json')}", file=sys.stderr)
        any_fail = any(c["state"] == "FAIL" for c in manifest["checks"])
        return 1 if any_fail else 0

    elif args.cmd == "replay":
        with open(args.spec) as f:
            spec = json.load(f)
        res = replay_run(args.manifest, args.run, spec=spec)
        print(json.dumps(res, indent=1, default=str))
        return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
