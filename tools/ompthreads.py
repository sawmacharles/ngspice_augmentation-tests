"""tools/ompthreads.py -- OpenMP thread-count rule for concurrent ngspice
processes (card #527, round 14; research/transient-acceleration.md section 18).

ngspice reads its OpenMP thread count from the `num_threads` variable (set to
8 by the build's `spinit`), NOT from OMP_NUM_THREADS.  N concurrent processes
each with 8 threads oversubscribe a machine once N*8 > cores; measured on a second workstation
(32 cores, 12 concurrent BSIM4 points): 2 threads/process = +19.7 % points per
hour over 8, 4 threads +13 %, (8 procs x 4 threads +10.5 %, 6 x 5 +0.6 %).
Rule: threads = clamp(cores // n_procs, 1, 8).  Decks without OpenMP devices
(switching/B-source decks) do not care (measured flat).

The override is delivered without touching the deck: SPICE_SCRIPTS points at a
copy of the build's scripts directory whose `spinit` has the new num_threads.
"""
import os
import re
import shutil

DEFAULT_CAP = 8      # spinit's value; the measured optimum for one process


def auto_threads(n_procs, cores=None, cap=DEFAULT_CAP):
    """clamp(cores // n_procs, 1, cap)."""
    cores = cores or os.cpu_count() or 1
    return max(1, min(cap, cores // max(1, int(n_procs))))


def resolve_setting(setting, n_procs, cores=None):
    """setting: None/'auto' -> rule; positive int (or digit string) -> that
    many; 0/'default'/'off' -> None (leave the build's spinit alone)."""
    if setting is None or str(setting).lower() == "auto":
        return auto_threads(n_procs, cores)
    if str(setting).lower() in ("default", "off", "0"):
        return None
    n = int(setting)
    if n < 1:
        raise ValueError(f"omp threads must be >= 1, got {setting}")
    return n


def find_scripts_dir(ngspice):
    """The build's share/ngspice/scripts directory (holds spinit), or None."""
    env = os.environ.get("SPICE_SCRIPTS")
    if env and os.path.isfile(os.path.join(env, "spinit")):
        return env
    exe = ngspice if os.path.sep in ngspice else (shutil.which(ngspice) or ngspice)
    bindir = os.path.dirname(os.path.realpath(exe))
    for cand in (os.path.join(bindir, "..", "share", "ngspice", "scripts"),
                 os.path.join(os.path.dirname(os.path.abspath(exe)), "..", "share", "ngspice", "scripts")):
        if os.path.isfile(os.path.join(cand, "spinit")):
            return os.path.normpath(cand)
    return None


_NT_RE = re.compile(r"(?m)^\s*set\s+num_threads\s*=\s*(\d+)")


def spinit_env(ngspice, threads, workdir):
    """Extra environment ({'SPICE_SCRIPTS': dir}) giving `threads` OpenMP
    threads per process, or {} if threads is None, already what spinit sets,
    or the scripts directory cannot be found (then ngspice keeps its default;
    the caller records that via `threads_applied`)."""
    if threads is None:
        return {}
    src = find_scripts_dir(ngspice)
    if not src:
        return {}
    text = open(os.path.join(src, "spinit")).read()
    m = _NT_RE.search(text)
    if m and int(m.group(1)) == threads:
        return {}
    dst = os.path.join(os.path.abspath(workdir), f".spice_scripts_t{threads}")
    if not os.path.isfile(os.path.join(dst, "spinit")):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copytree(src, dst, dirs_exist_ok=True)
        text = _NT_RE.sub(f"set num_threads={threads}", text) if m else text + f"\nset num_threads={threads}\n"
        with open(os.path.join(dst, "spinit"), "w") as f:
            f.write(text)
    return {"SPICE_SCRIPTS": dst}
