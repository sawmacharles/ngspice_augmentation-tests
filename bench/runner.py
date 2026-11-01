#!/usr/bin/env python3
"""bench/runner.py -- benchmark runner for the transient-ACCELERATION work
(ngspice_augmentation card #467, plan T9.1).

Runs each bench deck N times, memory-capped via `systemd-run --user --scope`,
wrapped in `/usr/bin/time -v` for wall/CPU/RSS, and records:
  - wall clock (from /usr/bin/time -v "Elapsed (wall clock) time")
  - user/sys CPU seconds, max RSS
  - ngspice's own `rusage cputime` phase split (parse+setup vs the tran loop,
    from the "elapsed time since last call" deltas already in each deck's
    .control block)
  - named .meas values printed via the deck's `print` line (for
    accel_same_results comparison against a reference JSON)
  - `uptime` load average at the moment of each run (other users' jobs
    may be running concurrently -- this is recorded, not corrected
    for)

Usage:
  bench/runner.py --ngspice /path/to/ngspice-aug --deck a_pex_buck --n 5
  bench/runner.py --ngspice /path/to/ngspice-aug --all --n 5 --label baseline-s12
  bench/runner.py --ngspice /path/to/ngspice-aug --deck c_pvt_separate --n 3

Output: JSON results under
  ~/Work/data/ngspice_augmentation/accel/bench/<deck>/results/<label>_<ts>.json
raw per-run logs under .../<deck>/logs/. Never under the repo.
"""
import argparse
import glob
import json
import os
import pathlib
import re
import statistics
import subprocess
import time

REPO = pathlib.Path(__file__).resolve().parent.parent
DECKS_DIR = REPO / "bench" / "decks"
DATA_ROOT = pathlib.Path(os.path.expanduser("~/Work/data/ngspice_augmentation/accel/bench"))

# deck registry: name -> (relative .cir path from DECKS_DIR, memory cap, kind)
DECKS = {
    "a_pex_buck": ("a_pex_buck.cir", "6G", "single"),
    "b_long_switch": ("b_long_switch.cir", "2G", "single"),
    "c_pvt_loop": ("c_pvt_loop.cir", "2G", "single"),
    "c_pvt_separate": (None, "4G", "pvt_separate"),
    "d_analog_block": ("d_analog_block.cir", "2G", "single"),
    # card #478 T9.2 task 2/3 additions:
    "d_analog_block_long": ("d_analog_block_long.cir", "2G", "single"),  # same 300
        # BSIM4 instances as d_analog_block, tstop stretched to 50us (>=10s
        # wall) so an OMP-schedule change has enough Newton iterations to
        # show a stable signal -- see bench/gen_analog_block.py --tstop-ns.
    "e_pex_tree": ("e_pex_tree.cir", "6G", "single"),  # a_pex_buck's same
        # power stage/60 BSIM4 fingers, but the NxN grid mesh replaced by a
        # tree-shaped synthetic PEX net (branches/local meshes/coupling
        # caps) -- bench/gen_pex_buck.py --topology tree; see bench/README.md.
}

TIME_RE = {
    "wall": re.compile(r"Elapsed \(wall clock\) time.*?: (?:(\d+):)?(\d+):(\d+(?:\.\d+)?)"),
    "user": re.compile(r"User time \(seconds\): (\d+(?:\.\d+)?)"),
    "sys": re.compile(r"System time \(seconds\): (\d+(?:\.\d+)?)"),
    "maxrss": re.compile(r"Maximum resident set size \(kbytes\): (\d+)"),
}
RUSAGE_RE = re.compile(r"(elapsed|CPU) time since last call \(seconds\) = ([\d.]+)")
MEAS_LINE_RE = re.compile(r"^(\w[\w#]*)\s*=\s*([-+0-9.eE]+)\s*$")


def parse_time_v(text):
    out = {}
    m = TIME_RE["wall"].search(text)
    if m:
        h = int(m.group(1)) if m.group(1) else 0
        mi = int(m.group(2))
        s = float(m.group(3))
        out["wall_s"] = h * 3600 + mi * 60 + s
    for k in ("user", "sys"):
        m = TIME_RE[k].search(text)
        if m:
            out[k + "_s"] = float(m.group(1))
    m = TIME_RE["maxrss"].search(text)
    if m:
        out["maxrss_kb"] = int(m.group(1))
    return out


def parse_rusage_phases(text):
    return [float(m.group(2)) for m in RUSAGE_RE.finditer(text)]


def parse_meas(text):
    vals = {}
    for line in text.splitlines():
        m = MEAS_LINE_RE.match(line.strip())
        if m:
            try:
                vals[m.group(1)] = float(m.group(2))
            except ValueError:
                pass
    return vals


def load_avg():
    try:
        return os.getloadavg()
    except OSError:
        return (None, None, None)


def run_one(ngspice_bin, cir_path, mem_cap, workdir, timeout=1200):
    cmd = [
        "systemd-run", "--user", "--scope", "--quiet",
        f"-p", f"MemoryMax={mem_cap}",
        "--", "/usr/bin/time", "-v", ngspice_bin, "-b", str(cir_path),
    ]
    t0 = time.time()
    la_before = load_avg()
    proc = subprocess.run(cmd, cwd=str(workdir), capture_output=True, text=True, timeout=timeout)
    t1 = time.time()
    la_after = load_avg()
    return {
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "wall_measured_s": t1 - t0,
        "load_before": la_before,
        "load_after": la_after,
    }


def bench_single(ngspice_bin, deck, n, label):
    cir_rel, mem_cap, kind = DECKS[deck]
    cir_path = DECKS_DIR / cir_rel
    workdir = DECKS_DIR
    outdata = DATA_ROOT / deck
    (outdata / "logs").mkdir(parents=True, exist_ok=True)
    (outdata / "results").mkdir(parents=True, exist_ok=True)

    runs = []
    for i in range(n):
        r = run_one(ngspice_bin, cir_path, mem_cap, workdir)
        tv = parse_time_v(r["stderr"])
        phases = parse_rusage_phases(r["stdout"])
        meas = parse_meas(r["stdout"])
        logf = outdata / "logs" / f"{label}_{deck}_run{i}.log"
        logf.write_text(
            f"=== cmd stdout ===\n{r['stdout']}\n=== cmd stderr (/usr/bin/time -v) ===\n{r['stderr']}\n"
        )
        run_rec = {
            "run": i,
            "returncode": r["returncode"],
            "wall_measured_s": r["wall_measured_s"],
            "load_before": r["load_before"],
            "load_after": r["load_after"],
            "rusage_phases_s": phases,
            "meas": meas,
            **tv,
        }
        runs.append(run_rec)
        print(f"  [{deck}] run {i+1}/{n}: wall={run_rec.get('wall_s', run_rec['wall_measured_s']):.3f}s "
              f"rc={r['returncode']} load={r['load_before']}")

    walls = [r.get("wall_s", r["wall_measured_s"]) for r in runs if r["returncode"] == 0]
    result = {
        "deck": deck,
        "label": label,
        "ngspice_bin": ngspice_bin,
        "n": n,
        "mem_cap": mem_cap,
        "median_wall_s": statistics.median(walls) if walls else None,
        "all_wall_s": walls,
        "runs": runs,
        "ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    ts = time.strftime("%Y%m%d-%H%M%S")
    outpath = outdata / "results" / f"{label}_{ts}.json"
    outpath.write_text(json.dumps(result, indent=2))
    print(f"[{deck}] median wall = {result['median_wall_s']} s over {len(walls)}/{n} ok runs -> {outpath}")
    return result


def bench_pvt_separate(ngspice_bin, n, label, parallel=1):
    """Run the 50 pvt_points/pt_NNN.cir files as separate processes; the
    metric is the TOTAL wall clock to finish all 50 (representing the
    many-short-points workload). One systemd-run scope wraps each process;
    with parallel>1 several run concurrently (still each individually
    memory-capped)."""
    pts_dir = DECKS_DIR / "pvt_points"
    pts = sorted(glob.glob(str(pts_dir / "pt_*.cir")))
    assert pts, "no pvt point files -- run bench/gen_pvt.py first"
    outdata = DATA_ROOT / "c_pvt_separate"
    (outdata / "logs").mkdir(parents=True, exist_ok=True)
    (outdata / "results").mkdir(parents=True, exist_ok=True)

    runs = []
    for i in range(n):
        t0 = time.time()
        la_before = load_avg()
        per_point = []
        if parallel <= 1:
            for p in pts:
                r = run_one(ngspice_bin, pathlib.Path(p), "1G", pts_dir, timeout=120)
                per_point.append({"file": os.path.basename(p), "rc": r["returncode"],
                                   "wall": parse_time_v(r["stderr"]).get("wall_s")})
        else:
            import concurrent.futures as cf
            with cf.ThreadPoolExecutor(max_workers=parallel) as ex:
                futs = {ex.submit(run_one, ngspice_bin, pathlib.Path(p), "1G", pts_dir, 120): p for p in pts}
                for fut in cf.as_completed(futs):
                    p = futs[fut]
                    r = fut.result()
                    per_point.append({"file": os.path.basename(p), "rc": r["returncode"],
                                       "wall": parse_time_v(r["stderr"]).get("wall_s")})
        t1 = time.time()
        la_after = load_avg()
        total_wall = t1 - t0
        nfail = sum(1 for pp in per_point if pp["rc"] != 0)
        runs.append({"run": i, "total_wall_s": total_wall, "n_points": len(pts),
                     "n_fail": nfail, "load_before": la_before, "load_after": la_after,
                     "per_point": per_point})
        print(f"  [c_pvt_separate] run {i+1}/{n}: total_wall={total_wall:.3f}s "
              f"({len(pts)} points, parallel={parallel}, fail={nfail})")

    walls = [r["total_wall_s"] for r in runs if r["n_fail"] == 0]
    result = {
        "deck": "c_pvt_separate", "label": label, "ngspice_bin": ngspice_bin, "n": n,
        "parallel": parallel, "median_wall_s": statistics.median(walls) if walls else None,
        "all_wall_s": walls, "runs": runs, "ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    ts = time.strftime("%Y%m%d-%H%M%S")
    outpath = outdata / "results" / f"{label}_par{parallel}_{ts}.json"
    outpath.write_text(json.dumps(result, indent=2))
    print(f"[c_pvt_separate] median total wall (parallel={parallel}) = {result['median_wall_s']} s -> {outpath}")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ngspice", required=True, help="path to ngspice-aug binary/wrapper")
    ap.add_argument("--deck", help="one deck name from the registry")
    ap.add_argument("--all", action="store_true", help="run every deck")
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--label", default="run")
    ap.add_argument("--pvt-parallel", type=int, default=1)
    args = ap.parse_args()

    print(f"uptime: {subprocess.run(['uptime'], capture_output=True, text=True).stdout.strip()}")
    print(f"nproc: {os.cpu_count()}")

    deck_names = list(DECKS.keys()) if args.all else [args.deck]
    for deck in deck_names:
        if deck == "c_pvt_separate":
            bench_pvt_separate(args.ngspice, args.n, args.label, parallel=args.pvt_parallel)
        else:
            bench_single(args.ngspice, deck, args.n, args.label)


if __name__ == "__main__":
    main()
