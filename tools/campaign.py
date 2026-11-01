#!/usr/bin/env python3
"""tools/campaign.py -- campaign-level scheduler (card #491, plan T9.6,
transient-ACCELERATION round 6).

Five rounds of prior work (research/transient-acceleration.md sections 1-10,
patches-wip/accel/README.md) squeezed single-run speedups out of the fork
itself. The remaining big, multiplicative lever is at CAMPAIGN level: many
PVT/MC/sweep points run one after another. Per the project's earlier
spice-research (cited in this card, not re-derived here): corner-count
reduction, single-process V/T sweeping (`alterparam`/`option temp` + reset,
already proven in this repo's own `bench/gen_pvt.py` -- ~30% off just from
not re-paying parse/setup per point), and parallel dispatch, with memory as
the ceiling on how far parallel dispatch can go (`min(cores, RAM /
per-sim RSS)`, not cores alone).

## What this module is, and is not

This is a **generic scheduler/runner engine**, not a circuit generator. A
"campaign" is a small Python module (see `bench/campaign_pvt_r3d.py`,
`bench/campaign_mc_diffpair.py` for the two worked examples this card
built) that supplies:

  - `POINTS`: a list of dicts, one per campaign point (sampled/enumerated
    however that campaign likes -- a PVT grid via `itertools.product`, an
    MC draw via `tools/mc.py`'s own seeded sampler, ...).
  - `MEASURES`: the bare names (no per-point prefix) every point's run is
    expected to report.
  - `build_point_deck(point, idx) -> str`: a complete, standalone deck for
    mode (a) -- one ngspice process per point. Must end its own `.control`
    block with the measured values auto-printed (ngspice does this itself
    after every named `.meas`) or an explicit `print`/`echo` line in the
    `name = value` / `MCMEAS name=value` shape this module's parser reads.
  - `build_loop_prelude(points) -> str` / `build_loop_entry(point, idx) ->
    str` / `build_loop_tail() -> str`: the three pieces of ONE multi-point
    deck for mode (b) -- a single ngspice session, `alterparam`/`option
    temp`/`reset`/`run`/`meas` per point (the exact idiom `bench/gen_pvt.py`
    already validated, generalised here so any campaign can use it, not
    just the one deck it was written for). Each entry's `.meas` names MUST
    be tagged `p{idx:03d}_<measure>` so the parser can attribute them.

Three run modes (the build task's (a)/(b)/(c)):

  (a) one process per point, scheduled by `run_campaign_a` with a
      memory-aware concurrency picker (see below).
  (b) one ngspice session looping over every point via
      `alterparam`/`option temp`/`reset`/`run` -- `run_campaign_b`. Note
      (card #473's own finding, cited in this card's brief): a `.param`
      whose RHS calls a random function (`agauss`/...) is re-evaluated at
      every point of *use*, not bound once -- irrelevant here because every
      value a campaign loop `alterparam`s in is a plain literal (sampled in
      Python, exactly `tools/mc.py`'s own design for the same reason), never
      an in-deck random function. `once()` (patch 0016,
      `.param p = once(agauss(...))`) is the fix for the case where a deck
      DOES want an in-ngspice random draw frozen at one shared value; not
      needed by anything in this module, cited here because the build task
      asks for the citation.
  (c) checkpoint-shared corner settling -- `run_campaign_c` -- for the
      narrower case where several DIFFERENT analyses share the same
      (T, Vin, Rload, ...) corner (not several DIFFERENT corners): settle
      once with `savestate` (docs/checkpoint.md), then `loadstate` per
      downstream view instead of re-running the shared settle prefix each
      time. This is the same mechanism `tools/fra.py --checkpoint` already
      uses in production (1.4-3x measured there); this module only adds the
      generic settle-once/branch-many plumbing, it does not re-derive
      checkpoint/restart correctness (docs/checkpoint.md owns that).

## Memory-aware scheduling

`pick_concurrency()` implements `min(cores, RAM / per-sim RSS)` (the ceiling
named in this card's brief) from a real measurement, not a guess: the first
`--probe-n` points of a mode-(a) run are executed serially under
`/usr/bin/time -v`, their peak RSS is taken (with a safety margin -- real
points can vary, and the probe sample is small), and the remaining points
run at `min(--max-procs, cores, floor((mem_budget - headroom) /
per_point_rss))` concurrency. `--mem-budget-mb` defaults to the box's
currently-free memory (`/proc/meminfo MemAvailable`) at call time, not total
RAM, since other users' work already resident is not this campaign's to
evict. `--force-concurrency` bypasses the picker entirely (kept for the
red-proof test in `tests/test_campaign.py`, which needs to show that
*skipping* the memory-aware picker on a memory-hungry deck under a real
`systemd-run --user --scope -p MemoryMax=...` cap gets processes killed,
while the picker's own choice does not).

`--per-point-mem-cap-mb` additionally wraps EACH point's own process in its
own `systemd-run --user --scope -p MemoryMax=<cap>` (the same idiom
`bench/runner.py` already uses) as a per-point safety net: a single point
that balloons past the probe's estimate is killed cleanly (recorded as
`"oom": true`, not a silent hang), without taking down its concurrent
siblings.

## Solver choice (opt-in)

`probe_solver()` / `solver_table()` implement the build task's "automatic
solver choice": run a deck's own first (short) point twice, once as-is and
once with `.option klu` injected, measure wall-clock, and report which is
faster. This does not change any deck by default -- `--solver-probe` on a
campaign run only *measures and records* the recommendation per campaign
(campaigns here are all PVT/MC point shapes, not PEX-mesh shapes, so round
5's own finding -- SPARSE default already right for small/well-ordered
matrices, KLU only needed for `a_pex_buck`-class high-fill-in grids -- is
expected to reproduce, and is checked against, not assumed; see
`research/transient-acceleration.md` section 11).

## Output

Every run writes a JSON manifest and a CSV results table under `outdir`
(caller-supplied; this project's own convention keeps these under
`~/Work/data/ngspice_augmentation/campaign/...`, never the repo -- see
`bench/campaign_pvt_r3d.py`/`bench/campaign_mc_diffpair.py`'s own CLI
wrappers for the convention this module itself does not hardcode).
"""
import argparse
import concurrent.futures
import csv
import importlib.util
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import ompthreads  # noqa: E402  (card #527 thread-count rule)

NGSPICE_DEFAULT = os.environ.get("NGSPICE", "ngspice")

MEAS_LINE_RE = re.compile(r"^(\w[\w#]*)\s*=\s*([-+0-9.eE]+)\s*$")
LOOP_MEAS_RE = re.compile(r"^p(\d+)_(\w[\w#]*)\s*=\s*([-+0-9.eE]+)\s*$")
# tools/mc.py's own tagged-print convention ("echo MCMEAS name=value") for
# analyses (like a bare .op) where a plain .meas card isn't reliable --
# campaigns that reuse tools/mc.py's sampling scheme (bench/campaign_mc_*)
# reuse this convention too instead of inventing a second one.
MCMEAS_LINE_RE = re.compile(r"^MCMEAS\s+(\w[\w#]*)\s*=\s*([-+0-9.eE]+)\s*$")
MAXRSS_RE = re.compile(r"Maximum resident set size \(kbytes\): (\d+)")
WALL_RE = re.compile(r"Elapsed \(wall clock\) time.*?: (?:(\d+):)?(\d+):(\d+(?:\.\d+)?)")
USER_RE = re.compile(r"User time \(seconds\): (\d+(?:\.\d+)?)")
SYS_RE = re.compile(r"System time \(seconds\): (\d+(?:\.\d+)?)")


# --------------------------------------------------------------------- #
# Low-level process running
# --------------------------------------------------------------------- #

def load_avg():
    try:
        return os.getloadavg()
    except OSError:
        return (None, None, None)


def mem_available_mb():
    """Currently-free memory (/proc/meminfo MemAvailable), not total RAM --
    other users' resident work is not this campaign's to evict."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    kb = int(line.split()[1])
                    return kb / 1024.0
    except OSError:
        pass
    return None


def _parse_time_v(text):
    out = {}
    m = WALL_RE.search(text)
    if m:
        h = int(m.group(1)) if m.group(1) else 0
        mi = int(m.group(2))
        s = float(m.group(3))
        out["wall_s"] = h * 3600 + mi * 60 + s
    m = USER_RE.search(text)
    if m:
        out["user_s"] = float(m.group(1))
    m = SYS_RE.search(text)
    if m:
        out["sys_s"] = float(m.group(1))
    m = MAXRSS_RE.search(text)
    if m:
        out["maxrss_kb"] = int(m.group(1))
    return out


def parse_meas(text, regex=MEAS_LINE_RE):
    """Parse per-line `name = value` (ngspice's own auto-printed `.meas`
    result, or an explicit `print`) AND `MCMEAS name=value` (tools/mc.py's
    tagged-print convention, for analyses where a plain `.meas` isn't
    reliable, e.g. a bare `.op`) into one flat dict -- a loop-mode deck's
    `pNNN_<name>` tags work identically whichever line shape produced them,
    see `_split_loop_measures`."""
    vals = {}
    for line in text.splitlines():
        line = line.strip()
        m = regex.match(line)
        if m:
            try:
                if regex is LOOP_MEAS_RE:
                    idx, name, val = m.groups()
                    vals.setdefault(int(idx), {})[name] = float(val)
                else:
                    vals[m.group(1)] = float(m.group(2))
            except ValueError:
                pass
            continue
        m2 = MCMEAS_LINE_RE.match(line)
        if m2:
            try:
                vals[m2.group(1)] = float(m2.group(2))
            except ValueError:
                pass
    return vals


GNU_TIME = "/usr/bin/time" if os.access("/usr/bin/time", os.X_OK) else None


class _Proc:
    def __init__(self, returncode, stdout, stderr):
        self.returncode, self.stdout, self.stderr = returncode, stdout, stderr


def _run_wait4(cmd, cwd, timeout, env=None):
    """Run cmd, return (_Proc, {"user_s","sys_s","maxrss_kb"}) from the
    child's own rusage (os.wait4), or (None, {}) on timeout."""
    import tempfile
    with tempfile.TemporaryFile() as fo, tempfile.TemporaryFile() as fe:
        p = subprocess.Popen(cmd, cwd=cwd, stdout=fo, stderr=fe, env=env)
        deadline = time.time() + timeout
        while True:
            pid, status, ru = os.wait4(p.pid, os.WNOHANG)
            if pid:
                break
            if time.time() > deadline:
                p.kill()
                os.wait4(p.pid, 0)
                return None, {}
            time.sleep(0.02)
        p.returncode = os.waitstatus_to_exitcode(status)
        fo.seek(0); fe.seek(0)
        out = fo.read().decode(errors="replace"); err = fe.read().decode(errors="replace")
    return (_Proc(p.returncode, out, err),
            {"user_s": ru.ru_utime, "sys_s": ru.ru_stime, "maxrss_kb": ru.ru_maxrss})


def run_ngspice(ngspice, deck_path, timeout, mem_cap_mb=None, cwd=None, extra_env=None):
    """Run one ngspice batch process under /usr/bin/time -v (always -- this
    is how maxrss/wall/user/sys are measured for every point, not just the
    probe), optionally further wrapped in its own systemd-run memory-capped
    scope (per-point safety net, `--per-point-mem-cap-mb`)."""
    # Normalise to absolute BEFORE any cwd= override is applied to the
    # child -- deck_path is built by callers relative to the *parent*
    # process's cwd (os.path.join(outdir, ...) with a caller-supplied,
    # possibly-relative outdir); resolving it here with the parent's own
    # os.path.abspath (not the child's, post-chdir, cwd) avoids a
    # doubled-path FileNotFound whenever both a relative outdir and an
    # explicit cwd= are in play at once (mode (a)'s per-point cwd, mode
    # (c)'s per-corner cwd).
    deck_path = os.path.abspath(str(deck_path))
    # GNU time is not installed everywhere (a second workstation has none, found 2026-09-29);
    # without it the child's rusage comes from os.wait4 instead.
    inner = ([GNU_TIME, "-v"] if GNU_TIME else []) + [ngspice, "-b", deck_path]
    if mem_cap_mb:
        cmd = ["systemd-run", "--user", "--scope", "--quiet",
               "-p", f"MemoryMax={int(mem_cap_mb)}M",
               # swap would hide the cap (that workstation has swap, found 2026-09-29)
               "-p", "MemorySwapMax=0", "--"] + inner
    else:
        cmd = inner
    t0 = time.time()
    la_before = load_avg()
    env = dict(os.environ, **extra_env) if extra_env else None
    if GNU_TIME:
        try:
            proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=env)
        except subprocess.TimeoutExpired:
            return {"ok": False, "error": "timeout", "wall_s": time.time() - t0,
                    "load_before": la_before, "measures": {}}
        tv = _parse_time_v(proc.stderr)
    else:
        proc, tv = _run_wait4(cmd, cwd, timeout, env=env)
        if proc is None:
            return {"ok": False, "error": "timeout", "wall_s": time.time() - t0,
                    "load_before": la_before, "measures": {}}
    wall_measured = time.time() - t0
    oom = (proc.returncode in (-9, 137)) or "Killed" in proc.stderr
    measures = parse_meas(proc.stdout)
    result = {
        "ok": proc.returncode == 0,
        "returncode": proc.returncode,
        "oom": oom,
        "wall_s": tv.get("wall_s", wall_measured),
        "wall_measured_s": wall_measured,
        "user_s": tv.get("user_s"),
        "sys_s": tv.get("sys_s"),
        "maxrss_kb": tv.get("maxrss_kb"),
        "measures": measures,
        "load_before": la_before,
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-2000:],
    }
    return result


# --------------------------------------------------------------------- #
# Memory-aware concurrency picker
# --------------------------------------------------------------------- #

def pick_concurrency(per_point_rss_kb, mem_budget_mb, cores, max_procs=None,
                      headroom_mb=512):
    """min(cores, RAM / per-sim RSS), from a real per-point RSS measurement.
    headroom_mb is subtracted from the budget before dividing (leave room
    for the scheduler process itself, the OS, other users) -- default 512MB,
    overridable. Always returns >= 1 (never refuses to run at all)."""
    usable_mb = max(0.0, mem_budget_mb - headroom_mb)
    per_point_mb = max(per_point_rss_kb / 1024.0, 1.0)
    by_mem = max(1, int(usable_mb // per_point_mb))
    n = min(cores, by_mem)
    if max_procs:
        n = min(n, max_procs)
    return max(1, n)


# --------------------------------------------------------------------- #
# Mode (a): one process per point, memory-aware concurrency
# --------------------------------------------------------------------- #

def run_campaign_a(points, build_point_deck, ngspice, outdir, timeout=300,
                    max_procs=None, mem_budget_mb=None, headroom_mb=512,
                    per_point_mem_cap_mb=None, probe_n=3, retries=1,
                    force_concurrency=None, cores=None, progress=None,
                    keep_decks=True, omp_threads="auto"):
    os.makedirs(outdir, exist_ok=True)
    decks_dir = os.path.join(outdir, "decks")
    os.makedirs(decks_dir, exist_ok=True)
    n = len(points)
    cores = cores or os.cpu_count() or 1
    if mem_budget_mb is None:
        avail = mem_available_mb()
        mem_budget_mb = avail if avail else 2048.0

    def _write_deck(idx):
        text = build_point_deck(points[idx], idx)
        path = os.path.join(decks_dir, f"pt_{idx:04d}.cir")
        with open(path, "w") as f:
            f.write(text)
        return path

    omp_env = [None]   # set after the probe phase, once the concurrency is known

    def _run_one(idx, mem_cap=per_point_mem_cap_mb):
        path = _write_deck(idx)
        last = None
        for attempt in range(max(1, retries)):
            last = run_ngspice(ngspice, path, timeout, mem_cap_mb=mem_cap, cwd=decks_dir,
                               extra_env=omp_env[0])
            if last["ok"]:
                break
        last["idx"] = idx
        last["point"] = points[idx]
        last["deck"] = path
        last["attempts"] = attempt + 1
        return last

    t_campaign0 = time.time()

    # --- probe phase: run the first probe_n points serially, measure RSS ---
    probe_n = max(0, min(probe_n, n))
    probe_results = []
    for i in range(probe_n):
        r = _run_one(i)
        probe_results.append(r)
        if progress:
            progress("probe", i + 1, probe_n)

    per_point_rss_kb = 0
    for r in probe_results:
        if r.get("maxrss_kb"):
            per_point_rss_kb = max(per_point_rss_kb, r["maxrss_kb"])
    if per_point_rss_kb == 0:
        per_point_rss_kb = 262144  # 256MB fallback if /usr/bin/time gave nothing
    per_point_rss_kb = int(per_point_rss_kb * 1.25)  # safety margin over a small probe sample

    if force_concurrency:
        concurrency = force_concurrency
    else:
        concurrency = pick_concurrency(per_point_rss_kb, mem_budget_mb, cores,
                                        max_procs=max_procs, headroom_mb=headroom_mb)

    # --- OpenMP threads per process (card #527): with `concurrency` processes
    # at once, ngspice's default 8 threads each oversubscribe the machine;
    # threads = clamp(cores // concurrency, 1, 8) (measured +10..+20 % points
    # per hour at 8..12 concurrent BSIM4 points on 32 cores; flat for decks
    # without OpenMP devices).  omp_threads: "auto" (rule), N, or 0/"default"
    # (leave the build's spinit alone).  The probe phase above is serial and
    # keeps the default.
    omp_n = ompthreads.resolve_setting(omp_threads, concurrency, cores) if concurrency > 1 else None
    omp_env[0] = ompthreads.spinit_env(ngspice, omp_n, outdir) or None
    omp_applied = omp_n if omp_env[0] else None

    # --- scheduled phase: remaining points, at the chosen concurrency ---
    remaining = list(range(probe_n, n))
    results = [None] * n
    for r in probe_results:
        results[r["idx"]] = r

    done_count = [probe_n]
    if remaining:
        with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as pool:
            futs = {pool.submit(_run_one, i): i for i in remaining}
            for fut in concurrent.futures.as_completed(futs):
                i = futs[fut]
                results[i] = fut.result()
                done_count[0] += 1
                if progress:
                    progress("run", done_count[0], n)

    wall_total = time.time() - t_campaign0
    manifest = {
        "mode": "a", "n": n, "ngspice": ngspice,
        "cores": cores, "mem_budget_mb": mem_budget_mb, "headroom_mb": headroom_mb,
        "per_point_rss_kb_probed": per_point_rss_kb, "probe_n": probe_n,
        "concurrency": concurrency, "force_concurrency": force_concurrency,
        "omp_threads_setting": omp_threads, "omp_threads_applied": omp_applied,
        "per_point_mem_cap_mb": per_point_mem_cap_mb,
        "wall_total_s": wall_total,
        "n_ok": sum(1 for r in results if r and r["ok"]),
        "n_oom": sum(1 for r in results if r and r.get("oom")),
        "results": results,
    }
    _write_manifest(outdir, manifest)
    return manifest


# --------------------------------------------------------------------- #
# Mode (b): single-process loop (alterparam/option temp/reset/run)
# --------------------------------------------------------------------- #

def build_loop_deck(points, build_loop_prelude, build_loop_entry, build_loop_tail):
    parts = [build_loop_prelude(points)]
    for idx, p in enumerate(points):
        parts.append(build_loop_entry(p, idx))
    parts.append(build_loop_tail())
    return "".join(parts)


def run_campaign_b(points, build_loop_prelude, build_loop_entry, build_loop_tail,
                    ngspice, outdir, timeout=1800, mem_cap_mb=None):
    os.makedirs(outdir, exist_ok=True)
    deck_text = build_loop_deck(points, build_loop_prelude, build_loop_entry, build_loop_tail)
    deck_path = os.path.join(outdir, "loop.cir")
    with open(deck_path, "w") as f:
        f.write(deck_text)

    t0 = time.time()
    r = run_ngspice(ngspice, deck_path, timeout, mem_cap_mb=mem_cap_mb, cwd=outdir)
    wall_total = time.time() - t0

    # run_ngspice's own `measures` field is parsed from the FULL stdout (not
    # the truncated tail kept for the manifest), via MEAS_LINE_RE -- which
    # already matches "pNNN_name = value" lines (\w includes digits/
    # underscore), just flatly keyed by the whole tagged name. Split that
    # flat dict into {idx: {measure: value}} here -- no second ngspice run
    # needed (a second run would silently double this mode's own wall-clock
    # cost, which would corrupt the very number this mode exists to report).
    per_point = _split_loop_measures(r["measures"])

    n = len(points)
    results = [None] * n
    for idx in range(n):
        vals = per_point.get(idx, {})
        results[idx] = {
            "idx": idx, "point": points[idx], "ok": r["ok"] and bool(vals),
            "measures": vals,
        }

    manifest = {
        "mode": "b", "n": n, "ngspice": ngspice, "deck": deck_path,
        "wall_total_s": wall_total, "run": {k: v for k, v in r.items()
                                             if k not in ("stdout_tail",)},
        "n_ok": sum(1 for x in results if x["ok"]),
        "results": results,
    }
    _write_manifest(outdir, manifest)
    return manifest


_LOOP_TAG_RE = re.compile(r"^p(\d+)_(\w+)$")


def _split_loop_measures(flat_measures):
    per_point = {}
    for k, v in flat_measures.items():
        m = _LOOP_TAG_RE.match(k)
        if m:
            idx, name = int(m.group(1)), m.group(2)
            per_point.setdefault(idx, {})[name] = v
    return per_point


# --------------------------------------------------------------------- #
# Mode (c): checkpoint-shared corner settling
# --------------------------------------------------------------------- #

def run_campaign_c(corners, build_settle_deck, build_view_deck, view_names,
                    ngspice, outdir, timeout=300):
    """For each corner (a dict, e.g. {"vin_v":12,"rload":3.3,"temp":27}):
      1. build_settle_deck(corner) -> a deck that runs a `.tran ... uic` up
         to the shared settle time and ends its .control with `savestate
         ckpt.bin` (lower-case filename -- docs/checkpoint.md's own
         lower-casing trap).
      2. for each name in view_names: build_view_deck(corner, view_name) ->
         a deck with the SAME topology/`.tran <step> <t2>` job and a
         `.control` that does `loadstate ckpt.bin` then whatever that
         view's own further run/meas is.
    Every corner gets its own subdirectory (own cwd) so `ckpt.bin` never
    collides between corners run concurrently by a caller that parallelises
    across corners (this function itself runs one corner at a time -- mode
    (a)'s scheduler is the place to add cross-corner parallelism if wanted;
    not done here, this mode's whole point is settle-once-per-corner, not
    max concurrency)."""
    os.makedirs(outdir, exist_ok=True)
    results = []
    for ci, corner in enumerate(corners):
        cdir = os.path.join(outdir, f"corner_{ci:03d}")
        os.makedirs(cdir, exist_ok=True)
        settle_path = os.path.join(cdir, "settle.cir")
        with open(settle_path, "w") as f:
            f.write(build_settle_deck(corner))
        t0 = time.time()
        settle_r = run_ngspice(ngspice, settle_path, timeout, cwd=cdir)
        settle_wall = time.time() - t0

        views = {}
        for vname in view_names:
            view_path = os.path.join(cdir, f"view_{vname}.cir")
            with open(view_path, "w") as f:
                f.write(build_view_deck(corner, vname))
            t0 = time.time()
            vr = run_ngspice(ngspice, view_path, timeout, cwd=cdir)
            vr["wall_s_own"] = time.time() - t0
            views[vname] = vr

        results.append({
            "idx": ci, "corner": corner,
            "settle": settle_r, "settle_wall_s": settle_wall,
            "views": views,
            "ok": settle_r["ok"] and all(v["ok"] for v in views.values()),
        })

    manifest = {
        "mode": "c", "n_corners": len(corners), "view_names": view_names,
        "ngspice": ngspice, "results": results,
        "n_ok": sum(1 for r in results if r["ok"]),
    }
    _write_manifest(outdir, manifest)
    return manifest


# --------------------------------------------------------------------- #
# Solver choice (opt-in)
# --------------------------------------------------------------------- #

def _inject_option(deck_text, option_line):
    """Insert a `.option ...` card right before the first `.control` line
    (top-level dot-cards must precede .control in ngspice's own grammar)."""
    lines = deck_text.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.strip().lower().startswith(".control"):
            lines.insert(i, option_line.rstrip("\n") + "\n")
            return "".join(lines)
    return deck_text + "\n" + option_line + "\n"  # no .control found -- append


_OPTIONS_KLU_TOKEN_RE = re.compile(r"(^[ \t]*\.options\b[^\n]*?)[ \t]*\bklu\b[ \t]*", re.IGNORECASE | re.MULTILINE)


def _force_sparse(deck_text):
    """Some of this repo's own bench decks bake `klu` in as a bare flag on
    their own `.options` line (e.g. `bench/gen_pex_buck.py`'s
    `.options reltol=... klu`, unconditional in every topology it
    generates -- found while building this prober: an earlier version that
    only ever ADDED `.option klu` degenerated to a no-op "klu vs klu"
    comparison on exactly the decks (a_pex_buck/e_pex_tree) this table
    most needs to say something about). Strip a bare `klu` token off any
    `.options` line so a genuine SPARSE-forced variant can be built and
    timed, symmetric with `_inject_option(..., ".option klu")` for the
    KLU-forced side -- reproducing round 5's actual "force klu vs force
    sparse" comparison generically, not just "as-is vs +klu"."""
    return _OPTIONS_KLU_TOKEN_RE.sub(r"\1", deck_text)


def probe_solver(ngspice, deck_path, timeout=300, cwd=None):
    """Run a deck twice -- once with SPARSE forced, once with KLU forced
    (both derived from the SAME base deck text, so this is meaningful even
    if the deck's own default already bakes one or the other in -- see
    `_force_sparse`) -- and report which is faster. Correctness note: KLU
    vs SPARSE are NOT bit-identical in general (different pivoting/
    round-off order -- expected, documented in
    research/transient-acceleration.md section 10, not a bug), so this
    function reports timing only; a campaign opting into the recommendation
    should treat it as a performance choice, not assume identical numbers
    across the switch."""
    with open(deck_path) as f:
        base_text = f.read()
    sparse_text = _force_sparse(base_text)
    klu_text = _inject_option(_force_sparse(base_text), ".option klu")
    tmp_dir = cwd or os.path.dirname(deck_path) or "."
    sparse_path = os.path.join(tmp_dir, "._solver_probe_sparse_" + os.path.basename(deck_path))
    klu_path = os.path.join(tmp_dir, "._solver_probe_klu_" + os.path.basename(deck_path))
    with open(sparse_path, "w") as f:
        f.write(sparse_text)
    with open(klu_path, "w") as f:
        f.write(klu_text)
    try:
        default_r = run_ngspice(ngspice, sparse_path, timeout, cwd=cwd)
        klu_r = run_ngspice(ngspice, klu_path, timeout, cwd=cwd)
    finally:
        for p in (sparse_path, klu_path):
            try:
                os.remove(p)
            except OSError:
                pass
    default_wall = default_r.get("wall_s") or float("inf")
    klu_wall = klu_r.get("wall_s") or float("inf")
    chosen = "klu" if (klu_r["ok"] and klu_wall < default_wall) else "sparse"
    if not default_r["ok"] and not klu_r["ok"]:
        chosen = "COULD-NOT-RUN"
    elif not default_r["ok"]:
        chosen = "klu"  # SPARSE didn't even finish -- KLU wins by default (a_pex_buck-class deck)
    elif not klu_r["ok"]:
        chosen = "sparse"
    speedup = (default_wall / klu_wall) if (klu_r["ok"] and default_r["ok"] and klu_wall > 0) else None
    return {
        "deck": deck_path, "sparse_wall_s": default_wall, "klu_wall_s": klu_wall,
        "sparse_ok": default_r["ok"], "klu_ok": klu_r["ok"],
        "chosen": chosen, "speedup_klu_over_sparse": speedup,
    }


def solver_table(ngspice, deck_paths, timeout=300):
    rows = []
    for d in deck_paths:
        rows.append(probe_solver(ngspice, d, timeout=timeout))
    return rows


# --------------------------------------------------------------------- #
# Results comparison (identity check between two manifests, or two modes)
# --------------------------------------------------------------------- #

def compare_results(results_a, results_b, measures, rtol=0.0, atol=0.0):
    """Compare two per-point results lists (same shape as a manifest's
    "results", indexed by point idx) on the given measure names. Returns
    (ok: bool, rows: list) -- rows has one entry per (idx, measure) with the
    two values and whether they matched within tolerance. rtol=atol=0 (the
    default) requires bit-identical float parsing, i.e. exact string->float
    agreement -- the strictest useful bar, matching this project's own
    accel_same_results convention of stating any nonzero tolerance and why,
    rather than defaulting to one."""
    rows = []
    ok = True
    n = min(len(results_a), len(results_b))
    for i in range(n):
        ra = results_a[i]
        rb = results_b[i]
        for m in measures:
            va = (ra or {}).get("measures", {}).get(m)
            vb = (rb or {}).get("measures", {}).get(m)
            if va is None or vb is None:
                match = False
                detail = "missing"
            else:
                tol = atol + rtol * max(abs(va), abs(vb))
                match = abs(va - vb) <= tol
                detail = f"diff={abs(va - vb):.6g} tol={tol:.6g}"
            if not match:
                ok = False
            rows.append({"idx": i, "measure": m, "a": va, "b": vb, "match": match, "detail": detail})
    return ok, rows


# --------------------------------------------------------------------- #
# I/O helpers
# --------------------------------------------------------------------- #

def _write_manifest(outdir, manifest):
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1, default=str)
    csv_path = os.path.join(outdir, "results.csv")
    results = manifest.get("results", [])
    if results and isinstance(results, list) and results and "measures" in (results[0] or {}):
        measure_names = sorted({k for r in results if r for k in r.get("measures", {})})
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["idx", "ok", "wall_s", "maxrss_kb", "returncode"] + measure_names)
            for r in results:
                if r is None:
                    continue
                w.writerow([r.get("idx"), r.get("ok"), r.get("wall_s"), r.get("maxrss_kb"),
                            r.get("returncode")] + [r.get("measures", {}).get(m) for m in measure_names])


def load_campaign_module(path):
    spec = importlib.util.spec_from_file_location("campaign_module", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------- #

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common_run_args(p):
        p.add_argument("--module", required=True, help="campaign module .py (POINTS/build_* functions)")
        p.add_argument("--ngspice", default=NGSPICE_DEFAULT)
        p.add_argument("--out", required=True)
        p.add_argument("--timeout", type=float, default=300)

    a = sub.add_parser("run-a", help="mode (a): one process per point, memory-aware concurrency")
    common_run_args(a)
    a.add_argument("--max-procs", type=int, default=None)
    a.add_argument("--mem-budget-mb", type=float, default=None)
    a.add_argument("--headroom-mb", type=float, default=512)
    a.add_argument("--per-point-mem-cap-mb", type=float, default=None)
    a.add_argument("--probe-n", type=int, default=3)
    a.add_argument("--retries", type=int, default=1)
    a.add_argument("--force-concurrency", type=int, default=None)
    a.add_argument("--omp-threads", default="auto",
                   help="OpenMP threads per ngspice process: auto (= clamp(cores//concurrency,1,8), default), N, or 0/default to leave the build's spinit value")

    b = sub.add_parser("run-b", help="mode (b): single-process alterparam loop")
    common_run_args(b)
    b.add_argument("--mem-cap-mb", type=float, default=None)

    ss = sub.add_parser("solver-table", help="probe KLU vs default solver on a list of decks")
    ss.add_argument("--ngspice", default=NGSPICE_DEFAULT)
    ss.add_argument("--decks", nargs="+", required=True)
    ss.add_argument("--timeout", type=float, default=300)
    ss.add_argument("--out", default=None)

    cp = sub.add_parser("compare", help="compare two manifest.json result sets")
    cp.add_argument("--a", required=True)
    cp.add_argument("--b", required=True)
    cp.add_argument("--measures", nargs="+", required=True)
    cp.add_argument("--rtol", type=float, default=0.0)
    cp.add_argument("--atol", type=float, default=0.0)

    args = ap.parse_args(argv)

    def progress(phase, done, total):
        print(f"\r[{phase}] {done}/{total}", end="", file=sys.stderr)

    if args.cmd == "run-a":
        mod = load_campaign_module(args.module)
        manifest = run_campaign_a(
            mod.POINTS, mod.build_point_deck, args.ngspice, args.out,
            timeout=args.timeout, max_procs=args.max_procs,
            mem_budget_mb=args.mem_budget_mb, headroom_mb=args.headroom_mb,
            per_point_mem_cap_mb=args.per_point_mem_cap_mb, probe_n=args.probe_n,
            retries=args.retries, force_concurrency=args.force_concurrency,
            progress=progress, omp_threads=args.omp_threads)
        print(file=sys.stderr)
        print(json.dumps({"n": manifest["n"], "n_ok": manifest["n_ok"], "n_oom": manifest["n_oom"],
                           "concurrency": manifest["concurrency"], "wall_total_s": manifest["wall_total_s"]}, indent=1))
        return 0 if manifest["n_ok"] == manifest["n"] else 1

    elif args.cmd == "run-b":
        mod = load_campaign_module(args.module)
        manifest = run_campaign_b(
            mod.POINTS, mod.build_loop_prelude, mod.build_loop_entry, mod.build_loop_tail,
            args.ngspice, args.out, timeout=args.timeout, mem_cap_mb=args.mem_cap_mb)
        print(json.dumps({"n": manifest["n"], "n_ok": manifest["n_ok"],
                           "wall_total_s": manifest["wall_total_s"]}, indent=1))
        return 0 if manifest["n_ok"] == manifest["n"] else 1

    elif args.cmd == "solver-table":
        rows = solver_table(args.ngspice, args.decks, timeout=args.timeout)
        print(json.dumps(rows, indent=1))
        if args.out:
            with open(args.out, "w") as f:
                json.dump(rows, f, indent=1)
        return 0

    elif args.cmd == "compare":
        with open(args.a) as f:
            ma = json.load(f)
        with open(args.b) as f:
            mb = json.load(f)
        ok, rows = compare_results(ma["results"], mb["results"], args.measures, rtol=args.rtol, atol=args.atol)
        n_mismatch = sum(1 for r in rows if not r["match"])
        print(json.dumps({"ok": ok, "n_rows": len(rows), "n_mismatch": n_mismatch}, indent=1))
        if not ok:
            for r in rows:
                if not r["match"]:
                    print(json.dumps(r), file=sys.stderr)
        return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
