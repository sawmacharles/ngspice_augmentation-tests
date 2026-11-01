#!/usr/bin/env python3
"""49_omp_threads.py -- handbook example 49 (card #560): OMP_NUM_THREADS does NOT limit ngspice's OpenMP threads; `set num_threads=N` does.

    python3 49_omp_threads.py NGSPICE OUTDIR

ngspice reads its OpenMP thread count from the variable `num_threads`, which the build's `spinit` sets to 8 at start-up; the environment variable
OMP_NUM_THREADS is overruled by that.  Deck: bench/decks/d_analog_block.cir (60 five-transistor OTA cells = 300 BSIM4 instances, 2 us transient), its
`wrdata` redirected into OUTDIR.  Arms, one ngspice at a time (raw bin/ngspice: its idle OpenMP threads spin, so CPU time shows the thread count):
    default                 (spinit: 8 threads)
    OMP_NUM_THREADS=2       environment only
    set num_threads=2       in the deck's .control, before the first analysis command
    ompthreads.py 2         tools/ompthreads.py: SPICE_SCRIPTS -> a copy of the scripts directory whose spinit says 2 (no deck change)
    set num_threads=1
Each arm prints  OMPRUN <arm>: wall <s> s, CPU <s> s, cores used (CPU/wall) <x>, peak RSS ...  and the arms' output vectors are compared (wrdata md5).
"""
import hashlib
import os
import re
import resource
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import ompthreads  # noqa: E402


def deck(out, set_threads):
    t = open(os.path.join(ROOT, "bench", "decks", "d_analog_block.cir")).read()
    t = re.sub(r"(?m)^wrdata \S+", "wrdata " + os.path.join(out, "wave.txt"), t)
    t = t.replace("rusage cputime\n", "")
    if set_threads:
        t = t.replace(".control\n", ".control\nset num_threads=%d\n" % set_threads, 1)
    return t


def run(ng, out, name, set_threads=None, env=None):
    d = os.path.join(out, name.replace(" ", "_").replace("=", ""))
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "deck.cir"), "w").write(deck(d, set_threads))
    e = dict(os.environ)
    e.pop("OMP_NUM_THREADS", None)
    e.update(env or {})
    t0 = time.time()
    p = subprocess.Popen(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", "deck.cir"],
                         cwd=d, env=e, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    o = p.communicate()[0]
    wall = time.time() - t0
    return d, wall, o


def main(ng, out):
    arms = [("default", None, {}), ("OMP_NUM_THREADS=2", None, {"OMP_NUM_THREADS": "2"}), ("set num_threads=2", 2, {}),
            ("ompthreads.py 2", None, ompthreads.spinit_env(ng, 2, out)), ("set num_threads=1", 1, {})]
    sums = {}
    for name, st, env in arms:
        r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
        d, wall, o = run(ng, out, name, st, env)
        r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime)
        w = open(os.path.join(d, "wave.txt"), "rb").read() if os.path.exists(os.path.join(d, "wave.txt")) else b""
        sums[name] = hashlib.md5(w).hexdigest()[:8]
        print("OMPRUN %s: wall %.2f s, CPU %.2f s, cores used (CPU/wall) %.2f, wrdata md5 %s" % (name, wall, cpu, cpu / wall, sums[name]), flush=True)
    print("OMPSAME: the wrdata file of every arm byte-identical: %s" % (len(set(sums.values())) == 1), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
