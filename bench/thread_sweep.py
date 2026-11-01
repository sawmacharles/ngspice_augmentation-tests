#!/usr/bin/env python3
"""bench/thread_sweep.py -- interleaved wall-clock sweep of ngspice builds x
OpenMP thread counts x decks (card #520, round 13).

Each deck is copied to a work dir with `set num_threads=N` injected as the
first line of its .control block (spinit's num_threads=8 is overridden; the
OMP_NUM_THREADS env var is NOT what ngspice reads: cktsetup.c calls
omp_set_num_threads(num_threads)) and `set numdgt=15` when --ident. A
"config" is  label=/path/to/ngspice[:N]  (N absent = leave spinit default 8).
Configs are rotated in start position each repetition; one ngspice at a time.
Records wall, user+sys CPU, load1 before the run, stdout digest.
Config may carry per-run env: label=/path/ngspice:N@K=V,K2=V2
Usage: thread_sweep.py --out X.json --n 3 --decks a,b --config base=/p/ngspice:1 ...
"""
import argparse, hashlib, json, os, re, shutil, subprocess, sys, time, statistics

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DECKS = {n: os.path.join(REPO, "bench/decks", n + ".cir") for n in
         "a_pex_buck b_long_switch c_pvt_loop d_analog_block d_analog_block_long e_pex_tree f_analog_block_bsim3 f_analog_block_bsim3_long g_ad_stress h_pex_tree_sparse".split()}
DECKS["r3d"] = os.path.join(REPO, "refs/r3d_vm_buck_deadtime.cir")
DECKS["r4d"] = os.path.join(REPO, "refs/r4d_pcm_buck_deadtime.cir")

def load1():
    return float(open("/proc/loadavg").read().split()[0])

def prep(deck, nthreads, work, ident, extra):
    src = open(DECKS[deck]).read().replace("../../stb/", os.path.join(REPO, "stb") + "/").replace("../stb/", os.path.join(REPO, "stb") + "/").splitlines()
    out, done = [], False
    for ln in src:
        out.append(ln)
        if not done and ln.strip().lower() == ".control":
            if ident: out.append("set numdgt=15")
            if nthreads: out.append("set num_threads=%d" % nthreads)
            for e in extra: out.append(e)
            done = True
    assert done, deck
    p = os.path.join(work, deck + ("_t%s" % nthreads) + ".cir")
    open(p, "w").write("\n".join(out) + "\n")
    return p

def run1(exe, deck_path, cwd, env):
    t0 = time.perf_counter()
    p = subprocess.Popen([exe, "-b", deck_path], cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = p.stdout.read()
    _, st, ru = os.wait4(p.pid, 0)
    wall = time.perf_counter() - t0
    return wall, ru.ru_utime + ru.ru_stime, os.waitstatus_to_exitcode(st), out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--decks", required=True); ap.add_argument("--config", action="append", required=True)
    ap.add_argument("--ident", action="store_true"); ap.add_argument("--maxload", type=float, default=12)
    ap.add_argument("--extra", action="append", default=[], help="extra control line injected")
    ap.add_argument("--env", action="append", default=[], help="K=V")
    a = ap.parse_args()
    work = os.path.expanduser("~/Work/data/ngspice_augmentation/r13/work"); os.makedirs(work, exist_ok=True)
    env = dict(os.environ); env.pop("OMP_NUM_THREADS", None)
    for kv in a.env: k, v = kv.split("=", 1); env[k] = v
    cfgs = []
    for c in a.config:
        lab, rest = c.split("=", 1)
        rest, _, cenv = rest.partition("@")   # optional @K=V,K=V per-config env
        exe, _, th = rest.partition(":")
        cfgs.append((lab, exe, int(th) if th else 0, dict(kv.split("=", 1) for kv in cenv.split(",") if kv)))
    res = []
    for deck in a.decks.split(","):
        paths = {c[0]: prep(deck, c[2], work, a.ident, a.extra) for c in cfgs}
        for rep in range(a.n):
            k = rep % len(cfgs)
            for lab, exe, th, cenv in cfgs[k:] + cfgs[:k]:
                while load1() > a.maxload: time.sleep(10)
                l = load1()
                w, cpu, rc, out = run1(exe, paths[lab], work, dict(env, **cenv))
                meas = b"\n".join(x for x in out.splitlines() if re.match(rb"^[A-Za-z_][\w.]*\s*=\s*[-+0-9.eE]+\s*$", x) and not re.search(rb"time|elapsed|cpu|memory|speed", x, re.I))
                r = dict(deck=deck, cfg=lab, rep=rep, wall=w, cpu=cpu, rc=rc, load1=l, digest=hashlib.md5(meas).hexdigest())
                res.append(r); print(json.dumps(r), flush=True)
                json.dump(res, open(a.out, "w"))
if __name__ == "__main__": main()
