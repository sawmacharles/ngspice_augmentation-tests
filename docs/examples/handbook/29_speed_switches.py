#!/usr/bin/env python3
"""29_speed_switches.py -- handbook example (card #534): the rounds 13-15 speed work, measured on THIS machine, with the
off-switches. Same binary, same deck, interleaved runs (start position rotated per repetition), median wall time, and an
output-identity check against the switch-off run.

  KLU PEX tran deck  (bench/decks/e_pex_tree.cir, ~11k capacitors, `.options klu`):
      default | NGSPICE_CAPBASE=0 (capacitor cache off) | NGSPICE_RESBASE=0 (R and C caches off)
      | NGSPICE_KLU_REFSCALE=1 (scaled KLU refactor, the pre-round-13 behaviour)
  PSS deck with tanh B-source switches  (26_pss_r3d_deadtime.cir):
      default | NGSPICE_FLATB=0 (flat B-source evaluator off)

Usage: 29_speed_switches.py NGSPICE_BINARY OUTDIR [REPS]
Each ngspice runs with `set num_threads=2` (the project's 2-thread rule; ngspice reads num_threads, not OMP_NUM_THREADS).
One ngspice at a time. Prints lines starting `SPEED` and `IDENT` (the handbook quotes them).
"""
import os, re, statistics, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
binary, outdir = sys.argv[1], sys.argv[2]
reps = int(sys.argv[3]) if len(sys.argv) > 3 else 5
os.makedirs(outdir, exist_ok=True)


def make_pex(path):
    s = open(os.path.join(ROOT, "bench", "decks", "e_pex_tree.cir")).read()
    s = s.replace("set filetype=ascii\n", "set filetype=ascii\nset num_threads=2\nset wr_singlescale\nset numdgt=15\n", 1)
    s = re.sub(r"wrdata \S+ ", "wrdata wave.txt ", s)
    open(path, "w").write(s)


def make_pss(path):
    s = open(os.path.join(HERE, "26_pss_r3d_deadtime.cir")).read()
    s = s.replace("set filetype=ascii\n", "set filetype=ascii\nset num_threads=2\n", 1)
    open(path, "w").write(s)


def run_once(deck, cfg_env, tag):
    d = os.path.join(outdir, tag)
    os.makedirs(d, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="2", OMP_WAIT_POLICY="passive")
    for k in ("NGSPICE_CAPBASE", "NGSPICE_RESBASE", "NGSPICE_KLU_REFSCALE", "NGSPICE_FLATB", "NGSPICE_CAPFLAT"):
        env.pop(k, None)
    env.update(cfg_env)
    t0 = time.time()
    p = subprocess.run([binary, "-b", deck], cwd=d, env=env, capture_output=True, text=True)
    return time.time() - t0, p.returncode, p.stdout + p.stderr


def read_cols(path):
    rows = []
    for ln in open(path):
        try:
            rows.append([float(x) for x in ln.split()])
        except ValueError:
            pass
    return rows


def maxrel(a, b):
    """max |delta| / column range over all columns; None if the row counts differ"""
    if len(a) != len(b):
        return None
    worst = 0.0
    for c in range(len(a[0])):
        ca, cb = [r[c] for r in a], [r[c] for r in b]
        rng = (max(ca) - min(ca)) or 1.0
        worst = max(worst, max(abs(x - y) for x, y in zip(ca, cb)) / rng)
    return worst


def campaign(name, deck, cfgs, outfile, reps):
    times = {k: [] for k in cfgs}
    keys = list(cfgs)
    for rep in range(reps):
        order = keys[rep % len(keys):] + keys[:rep % len(keys)]
        for k in order:
            w, rc, out = run_once(deck, cfgs[k], f"{name}_{k}")
            if rc != 0:
                print(f"SPEED {name} {k}: ngspice exit {rc}")
                sys.exit(1)
            times[k].append(w)
    base = statistics.median(times["default"])
    for k in keys:
        med = statistics.median(times[k])
        print(f"SPEED {name:8s} {k:14s} median {med:7.3f} s  [min {min(times[k]):.3f} max {max(times[k]):.3f}]  "
              f"{100 * (med / base - 1):+6.1f} % vs default  (n={reps})")
    ref = read_cols(os.path.join(outdir, f"{name}_default", outfile))
    for k in keys:
        if k == "default":
            continue
        o = read_cols(os.path.join(outdir, f"{name}_{k}", outfile))
        m = maxrel(ref, o)
        if m is None:
            print(f"IDENT {name} {k}: row counts differ ({len(ref)} vs {len(o)})")
        elif m == 0.0:
            print(f"IDENT {name} {k}: output bit-identical to default")
        else:
            print(f"IDENT {name} {k}: max |delta| / column range = {m:.1e} vs default")


pex = os.path.join(outdir, "e_pex_tree_2thr.cir")
make_pex(pex)
campaign("pex_tran", pex, {"default": {}, "CAPBASE=0": {"NGSPICE_CAPBASE": "0"},
                          "RESBASE=0": {"NGSPICE_RESBASE": "0"},
                          "KLU_REFSCALE=1": {"NGSPICE_KLU_REFSCALE": "1"}}, "wave.txt", reps)
pss = os.path.join(outdir, "26_2thr.cir")
make_pss(pss)
campaign("pss_r3d", pss, {"default": {}, "FLATB=0": {"NGSPICE_FLATB": "0"}}, "pss_r3d.txt", max(reps, 9))
