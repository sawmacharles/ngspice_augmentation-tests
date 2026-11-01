#!/usr/bin/env python3
"""opbench -- DC operating-point robustness benchmark harness (card #637).

Measures how reliably (and how cheaply) an ngspice binary finds the DC operating
point of a corpus of hard decks (tests/fixtures/op/*.cir plus the refs/*.cir
decks turned into .op decks).  stdlib only.

  python3 tools/opbench.py --ngspice BIN --label NAME [--env K=V ...]
        [--options "gminsteps=0 srcsteps=0"] [--nodeset-seeds N] [--jobs 4]
        [--out DIR] [--decks GLOB] [--no-refs] [--seed 637] [--refcache DIR]
  python3 tools/opbench.py compare DIR_A DIR_B

Per run (deck x variant; variant 'plain' or 'nsK') it writes a row of
DIR/results.csv, DIR/solutions/<id>.json, DIR/logs/<id>.log, DIR/rundecks/<id>.cir;
plus DIR/summary.json and a printed table.

How a run is built
  run deck = deck body (analysis cards .op/.tran/... and any .control block removed,
  .end removed) + `.options ...` (user) + optional `.nodeset` + a .control block:
      set numdgt=15 / op / rusage totiter / rusage equations / print all / quit
  `print all` lists every node voltage and every branch current of the OP plot.

What is parsed (see docs in the final report of card #637 for caveats)
  converged   : at least one 'name = value' line, all values finite, a
                'Total iterations' line, and none of: 'no convergence',
                'could not be simulated', 'Error', 'timed out', 'aborted'.
                ('Warning: singular matrix' lines are NOT fatal: they are printed
                by every failed intermediate attempt and also by runs that go on
                to converge via gmin stepping.)
  stage       : which algorithm produced the answer, from messages ngspice prints
                at default verbosity (no `set ngdebug` needed):
                  newton  - no fallback message
                  dyngmin - 'Dynamic gmin stepping completed'
                  gmin    - 'Gmin stepping completed' (classic gmin stepping)
                  dynsrc/src - 'Source stepping completed' (dynamic/classic variants)
                  optran  - 'Transient op finished successfully' (pseudo-transient)
                'tried' holds all stages entered ('Starting/Trying ...' lines).
  iters       : 'Total iterations' from `rusage totiter`: Newton iterations of the
                whole OP, summed over every fallback attempt (not just the last).
  unknowns    : 'Circuit Equations' from `rusage equations`.
Accuracy
  Each deck is also run once with reltol=1e-9 abstol=1e-15 vntol=1e-12 itl1=2000
  (default fallbacks) -> 'ref'.  A run's solution is compared on node voltages
  (not branch currents).  Beyond 1e-3*max(1,|v|) it is WRONG, unless the
  solution is a genuine other root ('other'): the solution is fed back as a full
  .nodeset (tight tolerances, 15 digits) and must converge to the same point
  (movement < 1e-2*max(1,|v|); a default-tolerance solution of a stiff circuit can be a few 1e-3 off its root) in <= OTHER_MAX_ITERS iterations.
"""
import argparse
import concurrent.futures as cf
import csv
import fnmatch
import glob
import hashlib
import json
import math
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join(ROOT, "tests", "fixtures", "op")
REFS = os.path.join(ROOT, "refs")
TIMEOUT = 60
OTHER_MAX_ITERS = 25    # sanity cap on the nodeset re-run; NB a nodeset run costs >= 5 iterations even from an exact root (pinned phase + release), so the spec's "<= 3" is not reachable; the test is "Newton does not leave the point"
REF_OPTIONS = "reltol=1e-9 abstol=1e-15 vntol=1e-12 itl1=2000"
WRONG_TOL = 1e-3
ABSURD_FACTOR = 1e6   # |node voltage| above this x max(1, largest source) counts as 'absurd' = not converged
OTHER_TOL = 1e-2   # a loosely converged (default reltol) point of a stiff circuit can sit ~reltol*cond from its root

CSV_FIELDS = ["deck", "variant", "seed", "converged", "status", "stage", "tried", "iters",
              "wall_s", "unknowns", "flag", "kcl_iters", "maxerr", "maxerr_node", "ref_ok", "ref_method"]

# analysis / output cards removed when turning a reference deck into an OP deck
DROP_CARDS = re.compile(r"^\s*\.(op|tran|pss|ac|pac|pstb|pxf|pnoise|stb|pz|noise|dc|four|meas|measure|fra|sens|tf|"
                        r"disto|print|plot|probe|save|ic|end)\b", re.I)
CTRL_START = re.compile(r"^\s*\.control\b", re.I)
CTRL_END = re.compile(r"^\s*\.endc\b", re.I)
VAL_LINE = re.compile(r"^([^\s=]+)\s*=\s*(\S+)\s*$")
FATAL = re.compile(r"no convergence|could not be simulated|^\s*error|timed out|aborted|"
                   r"doAnalyses: .*(failed|error)|simulation interrupted", re.I | re.M)

# message -> (stage, kind); kind 'try' = entered, 'ok' = completed
MSG = [
    (re.compile(r"starting dynamic gmin stepping", re.I), ("dyngmin", "try")),
    (re.compile(r"dynamic gmin stepping completed", re.I), ("dyngmin", "ok")),
    (re.compile(r"trying (true )?gmin stepping|starting (true )?gmin stepping", re.I), ("gmin", "try")),
    (re.compile(r"(?<!dynamic )gmin stepping completed", re.I), ("gmin", "ok")),
    (re.compile(r"starting dynamic source stepping|trying dynamic source", re.I), ("dynsrc", "try")),
    (re.compile(r"dynamic source stepping completed", re.I), ("dynsrc", "ok")),
    (re.compile(r"trying source stepping|starting source stepping", re.I), ("src", "try")),
    (re.compile(r"(?<!dynamic )source stepping completed", re.I), ("src", "ok")),
    (re.compile(r"transient op started|starting transient op|optran.*start", re.I), ("optran", "try")),
    (re.compile(r"transient op finished successfully|transient op.*completed", re.I), ("optran", "ok")),
]


# ----------------------------------------------------------------- deck handling
def read_text(path):
    with open(path, errors="replace") as f:
        return f.read()


def strip_deck(text, for_ref=False, base_dir=None):
    """Return deck lines without analysis cards, .control blocks and .end."""
    out, in_ctrl, dropping = [], False, False
    for ln in text.splitlines():
        if in_ctrl:
            if CTRL_END.match(ln):
                in_ctrl = False
            continue
        if CTRL_START.match(ln):
            in_ctrl = True
            continue
        if ln.startswith("+") and dropping:
            continue
        dropping = False
        if DROP_CARDS.match(ln):
            dropping = True
            continue
        if base_dir:
            m = re.match(r"^(\s*\.(?:include|inc|lib)\s+)(\S+)(.*)$", ln, re.I)
            if m and not os.path.isabs(m.group(2).strip('"')):
                ln = m.group(1) + os.path.normpath(os.path.join(base_dir, m.group(2).strip('"'))) + m.group(3)
        if base_dir and re.match(r"^\s*x\S*\s.*\sdiffprobe\s*$", ln, re.I):
            ln = ln.rstrip() + "_dm"   # refs/r1c uses the generic name; the lib defines _dm/_cm
        out.append(ln)
    return out


def corpus(decks_glob, use_refs, corpus_dir=None):
    items = []  # (name, path, is_ref)
    for p in sorted(glob.glob(os.path.join(corpus_dir or CORPUS, "*.cir"))):
        items.append((os.path.basename(p)[:-4], p, False))
    if use_refs:
        for p in sorted(glob.glob(os.path.join(REFS, "*.cir"))):
            items.append(("ref:" + os.path.basename(p)[:-4], p, True))
    if decks_glob:
        pats = decks_glob.split(",")
        items = [i for i in items if any(fnmatch.fnmatch(i[0], q) for q in pats)]
    return items


def source_range(lines):
    """min/max of DC voltage source values (ground, 0, always included)."""
    params = {}
    for ln in lines:
        m = re.match(r"^\s*\.param\s+(.*)$", ln, re.I)
        if m:
            for k, v in re.findall(r"(\w+)\s*=\s*([-+0-9.eE]+[a-zA-Z]*)\b", m.group(1)):
                params[k.lower()] = v
    vals = [0.0]
    for ln in lines:
        if not re.match(r"^\s*v", ln, re.I):
            continue
        toks = ln.replace("(", " ( ").replace(")", " ) ").split()
        rest = toks[3:]
        if not rest:
            continue
        cand = None
        low = [t.lower() for t in rest]
        if "dc" in low and low.index("dc") + 1 < len(rest):
            cand = rest[low.index("dc") + 1]
        elif low[0] in ("pulse", "sin", "pwl", "exp", "sffm", "ac"):
            nxt = [t for t in rest[1:] if t not in "()"]
            cand = nxt[0] if (nxt and low[0] in ("pulse", "sin", "exp")) else None
        else:
            cand = rest[0]
        if cand is None:
            continue
        v = spice_num(cand.strip("{}"), params)
        if v is not None:
            vals.append(v)
    return min(vals), max(vals)


SUFFIX = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}


def spice_num(tok, params=None):
    tok = tok.strip().lower()
    if params and tok in params:
        tok = params[tok].lower()
    m = re.match(r"^([-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?)([a-z]*)$", tok)
    if not m:
        return None
    x = float(m.group(1))
    suf = m.group(2)
    if suf.startswith("meg"):
        return x * 1e6
    return x * SUFFIX.get(suf[:1], 1.0) if suf else x


def optran_off_line(options):
    """`optran` command that deselects the pseudo-transient fallback (step size 0) and
    re-states gminsteps/srcsteps (the command overwrites them in the task)."""
    g = re.search(r"\bgminsteps\s*=\s*(-?\d+)", options or "")
    sr = re.search(r"\bsrcsteps\s*=\s*(-?\d+)", options or "")
    return "optran 1 %s %s 0 0 0" % (g.group(1) if g else "1", sr.group(1) if sr else "1")


def build_run_deck(body_lines, options, nodeset, no_optran=False):
    lines = list(body_lines)
    if options:
        lines.append(".options " + options)
    if nodeset:
        items = " ".join(f"v({n})={v:.15g}" for n, v in nodeset.items())
        # keep lines short: ngspice has an input line length limit
        chunk, cur = [], ""
        for it in items.split(" "):
            if len(cur) + len(it) > 200:
                chunk.append(cur)
                cur = ""
            cur += (" " if cur else "") + it
        chunk.append(cur)
        lines.append(".nodeset " + chunk[0])
        lines += ["+ " + c for c in chunk[1:]]
    lines += [".control", "set numdgt=15"] + ([optran_off_line(options)] if no_optran else []) + ["op", "rusage totiter", "rusage equations",
              "print all", "quit", ".endc", ".end"]
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------- running
def launch(ngspice, deck_path, env_extra, limit=True):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_THREAD_LIMIT"] = "1"   # ngspice calls omp_set_num_threads itself; only the limit keeps it to one thread
    env.update(env_extra)
    cmd = [ngspice, "-b", deck_path]
    if limit and shutil.which("systemd-run"):
        cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=4G", "-p", "MemorySwapMax=0"] + cmd
    t0 = time.time()
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                             errors="replace", env=env, start_new_session=True, cwd=os.path.dirname(deck_path))
    except OSError as e:
        return "", time.time() - t0, f"cannot start: {e}"
    try:
        out, _ = p.communicate(timeout=TIMEOUT)
        rc = p.returncode
        if rc < 0 or rc >= 128:
            return out, time.time() - t0, f"crash: exit status {rc}"
        return out, time.time() - t0, None
    except subprocess.TimeoutExpired:
        try:
            os.killpg(p.pid, 9)
        except OSError:
            pass
        out, _ = p.communicate()
        return out, time.time() - t0, "timeout"


def parse_log(log):
    """-> dict(converged,status,stage,tried,iters,unknowns,nodes,branches)"""
    nodes, branches = {}, {}
    for ln in log.splitlines():
        m = VAL_LINE.match(ln.strip())
        if not m:
            continue
        k = m.group(1).lower()
        if k in ("total", "circuit"):
            continue
        try:
            v = float(m.group(2))
        except ValueError:
            continue
        (branches if ("#" in k or k.startswith("@") or k.endswith("branch")) else nodes)[k] = v
    iters = unk = None
    m = re.search(r"^Total iterations\s*=\s*(\d+)", log, re.M)
    if m:
        iters = int(m.group(1))
    m = re.search(r"^Circuit Equations\s*=\s*(\d+)", log, re.M)
    if m:
        unk = int(m.group(1))
    tried, final = [], None
    for ln in log.splitlines():
        for rx, (st, kind) in MSG:
            if rx.search(ln):
                if kind == "try" and st not in tried:
                    tried.append(st)
                if kind == "ok":
                    final = st
                break
    allv = list(nodes.values()) + list(branches.values())
    finite = all(math.isfinite(x) for x in allv)
    fatal = FATAL.search(log)
    nonfinite_text = re.search(r"=\s*-?(nan|inf)\b", log, re.I)
    ok = bool(nodes) and finite and iters is not None and not fatal and not nonfinite_text
    status = "ok" if ok else ("noconv" if iters is not None or fatal else "nodata")
    return dict(converged=ok, status=status, stage=(final or "newton") if ok else "",
                tried=",".join(tried), iters=iters, unknowns=unk, nodes=nodes, branches=branches)


def run_one(ctx, deck_id, body, options, nodeset, tag, with_env=True, no_optran=None):
    """Run a built deck; return parse dict + wall + log."""
    wd = os.path.join(ctx["out"], "rundecks")
    os.makedirs(wd, exist_ok=True)
    path = os.path.join(wd, f"{deck_id}__{tag}.cir".replace(":", "_"))
    with open(path, "w") as f:
        f.write(build_run_deck(body, options, nodeset, ctx["no_optran"] if no_optran is None else no_optran))
    log, wall, err = launch(ctx["ngspice"], path, ctx["env"] if with_env else {}, ctx["limit"])
    r = parse_log(log)
    r["wall"] = wall
    r["log"] = log
    if r["converged"]:
        lo, hi = source_range(body)
        bound = ABSURD_FACTOR * max(1.0, abs(lo), abs(hi))
        big = max(abs(v) for v in r["nodes"].values())
        if big > bound:
            r["converged"], r["status"], r["stage"] = False, "absurd", ""
            log_note = f"\n[opbench] absurd solution: max |v| = {big:.3g} > {bound:.3g}\n"
            r["log"] += log_note
    if err:
        r["converged"], r["status"], r["stage"] = False, err.split(":")[0], ""
        r["log"] += f"\n[opbench] {err}\n"
    return r


def get_ref(ctx, name, body):
    key = hashlib.sha1((ctx["ngspice"] + "\0".join(body) + REF_OPTIONS + "v2" + json.dumps({}, sort_keys=True)).encode()).hexdigest()[:16]   # references never carry the variant's --env
    cdir = ctx["refcache"]
    os.makedirs(cdir, exist_ok=True)
    cpath = os.path.join(cdir, name.replace(":", "_") + ".json")
    if os.path.exists(cpath):
        try:
            c = json.load(open(cpath))
            if c.get("key") == key:
                return c
        except (OSError, ValueError):
            pass
    # reference ladder: tight tolerances with the default fallbacks; if that fails or is absurd,
    # retry without gmin stepping (source stepping first), then Newton + optran only.
    method, r = "", None
    for method, extra in (("default", ""), ("no-gmin", " gminsteps=0"), ("no-gmin-no-src", " gminsteps=0 srcsteps=0")):
        r = run_one(ctx, name, body, REF_OPTIONS + extra, None, "ref" if not extra else "ref" + str(len(extra)), with_env=False, no_optran=False)
        if r["converged"]:
            break
    if not r["converged"]:
        # last resort: any sane solution at default tolerances, polished with tight tolerances via .nodeset
        for i, extra in enumerate(("", "gminsteps=0", "srcsteps=0", "gminsteps=0 srcsteps=0")):
            q = run_one(ctx, name, body, extra, None, f"refloose{i}", with_env=False, no_optran=False)
            if q["converged"]:
                pol = run_one(ctx, name, body, REF_OPTIONS, dict(q["nodes"]), f"refpolish{i}", with_env=False, no_optran=False)
                r = pol if pol["converged"] else q
                method = "polished:" + (extra or "default") if pol["converged"] else "loose:" + (extra or "default")
                break
    c = dict(key=key, ok=r["converged"], nodes=r["nodes"], iters=r["iters"], wall=r["wall"], stage=r["stage"], method=method)
    with open(cpath, "w") as f:
        json.dump(c, f)
    with open(os.path.join(ctx["out"], "logs", name.replace(":", "_") + "__ref.log"), "w") as f:
        f.write(r["log"])
    return c


def compare_sol(sol, ref, tol=None):
    """max |v - vref| over common nodes, relative-scaled for the WRONG test."""
    worst, wn, bad = 0.0, "", False
    for k, vr in ref.items():
        if k not in sol:
            continue
        d = abs(sol[k] - vr)
        if d > worst:
            worst, wn = d, k
        if d > (tol or WRONG_TOL) * max(1.0, abs(vr)):
            bad = True
    return worst, wn, bad


def job(ctx, name, body, variant, seed, nodeset, ref):
    tag = variant
    r = run_one(ctx, name, body, ctx["options"], nodeset, tag)
    row = dict(deck=name, variant=variant, seed=seed if seed is not None else "", converged=int(r["converged"]),
               status=r["status"], stage=r["stage"], tried=r["tried"],
               iters=r["iters"] if r["iters"] is not None else "", wall_s=f"{r['wall']:.3f}",
               unknowns=r["unknowns"] if r["unknowns"] is not None else "", flag="", kcl_iters="", maxerr="",
               maxerr_node="", ref_ok=int(bool(ref and ref["ok"])), ref_method=(ref or {}).get("method", ""))
    with open(os.path.join(ctx["out"], "logs", f"{name}__{tag}.log".replace(":", "_")), "w") as f:
        f.write(r["log"])
    if r["converged"]:
        with open(os.path.join(ctx["out"], "solutions", f"{name}__{tag}.json".replace(":", "_")), "w") as f:
            json.dump(dict(nodes=r["nodes"], branches=r["branches"]), f)
        if ref and ref["ok"]:
            worst, wn, bad = compare_sol(r["nodes"], ref["nodes"])
            row["maxerr"], row["maxerr_node"] = f"{worst:.3e}", wn
            if not bad:
                row["flag"] = "ok"
            else:
                # is it a genuine other root?  Feed the whole solution back as a .nodeset under
                # tight tolerances: from a true root Newton stays (moves < WRONG_TOL), from garbage it moves.
                chk = run_one(ctx, name, body, REF_OPTIONS, dict(r["nodes"]), tag + "_kcl", with_env=False, no_optran=False)
                row["kcl_iters"] = chk["iters"] if chk["iters"] is not None else ""
                if chk["converged"] and chk["iters"] is not None and chk["iters"] <= OTHER_MAX_ITERS:
                    _, _, moved = compare_sol(chk["nodes"], r["nodes"], OTHER_TOL)
                    row["flag"] = "WRONG" if moved else "other"
                else:
                    row["flag"] = "WRONG"
                with open(os.path.join(ctx["out"], "logs", f"{name}__{tag}_kcl.log".replace(":", "_")), "w") as f:
                    f.write(chk["log"])
        else:
            row["flag"] = "noref"
    return row


def seeds_nodeset(base_seed, name, k, nodes, lo, hi):
    rng = random.Random(f"{base_seed}:{name}:{k}")
    return {n: rng.uniform(lo - 1.0, hi + 1.0) for n in sorted(nodes)}


# ----------------------------------------------------------------- driver
def run_bench(a):
    out = os.path.abspath(a.out)
    for sub in ("logs", "solutions", "rundecks"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    env = dict(kv.split("=", 1) for kv in (a.env or []))
    ctx = dict(ngspice=a.ngspice, env=env, out=out, options=(a.options or "").strip(),
               refcache=os.path.abspath(a.refcache) if a.refcache else os.path.join(out, "ref"),
               limit=not a.no_limit, no_optran=a.no_optran)
    items = corpus(a.decks, not a.no_refs, a.corpus_dir)
    if not items:
        sys.exit("no decks selected")
    bodies, refs, tasks = {}, {}, []
    t_all = time.time()
    # 1. references (parallel)
    for name, path, is_ref in items:
        bodies[name] = strip_deck(read_text(path), base_dir=os.path.dirname(path) if is_ref else None)
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(get_ref, ctx, n, bodies[n]): n for n in bodies}
        for f in cf.as_completed(futs):
            refs[futs[f]] = f.result()
    print(f"references: {sum(r['ok'] for r in refs.values())}/{len(refs)} converged", flush=True)
    # 2. runs
    for name, path, is_ref in items:
        body = bodies[name]
        tasks.append((name, body, "plain", None, None))
        if a.nodeset_seeds:
            lo, hi = source_range(body)
            nodes = [n for n in refs[name]["nodes"] if n != "0"]
            for k in range(1, a.nodeset_seeds + 1):
                tasks.append((name, body, f"ns{k}", f"{a.seed}:{name}:{k}",
                              seeds_nodeset(a.seed, name, k, nodes, lo, hi) if nodes else None))
    rows = []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(job, ctx, n, b, v, s, ns, refs[n]) for (n, b, v, s, ns) in tasks]
        for f in futs:
            rows.append(f.result())
    rows.sort(key=lambda r: (r["deck"], r["variant"]))
    with open(os.path.join(out, "results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)
    summ = summarize(rows)
    summ["label"] = a.label
    summ["ngspice"] = a.ngspice
    summ["options"] = ctx["options"]
    summ["env"] = env
    summ["nodeset_seeds"] = a.nodeset_seeds
    summ["seed"] = a.seed
    summ["elapsed_s"] = round(time.time() - t_all, 1)
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print_summary(summ, rows)
    return 0


def summarize(rows):
    def stats(rs):
        conv = [r for r in rs if int(r["converged"])]
        its = [int(r["iters"]) for r in conv if r["iters"] != ""]
        stages = {}
        for r in conv:
            stages[r["stage"]] = stages.get(r["stage"], 0) + 1
        return dict(runs=len(rs), converged=len(conv), conv_pct=round(100.0 * len(conv) / max(1, len(rs)), 1),
                    stages=stages, iters_median=statistics.median(its) if its else None, iters_total=sum(its),
                    wall_total_s=round(sum(float(r["wall_s"]) for r in rs), 2),
                    wrong=sum(r["flag"] == "WRONG" for r in rs), other=sum(r["flag"] == "other" for r in rs),
                    failed_decks=sorted({r["deck"] for r in rs if not int(r["converged"])}))
    plain = [r for r in rows if r["variant"] == "plain"]
    seeded = [r for r in rows if r["variant"] != "plain"]
    s = dict(decks=len({r["deck"] for r in rows}), all=stats(rows), plain=stats(plain))
    if seeded:
        s["seeded"] = stats(seeded)
    return s


def print_summary(s, rows):
    print(f"\n== {s.get('label','')}  decks={s['decks']}  elapsed={s.get('elapsed_s','?')}s ==")
    print(f"{'set':8} {'runs':>5} {'conv%':>6} {'iters med':>9} {'iters tot':>9} {'wall s':>8} {'WRONG':>5} {'other':>5}  stages")
    for k in ("plain", "seeded", "all"):
        if k in s:
            t = s[k]
            print(f"{k:8} {t['runs']:5d} {t['conv_pct']:6.1f} {str(t['iters_median']):>9} {t['iters_total']:9d} "
                  f"{t['wall_total_s']:8.1f} {t['wrong']:5d} {t['other']:5d}  {t['stages']}")
    f = s["plain"]["failed_decks"]
    print("plain-variant failures:", ", ".join(f) if f else "none")


# ----------------------------------------------------------------- compare
def load(d):
    rows = {}
    with open(os.path.join(d, "results.csv")) as f:
        for r in csv.DictReader(f):
            rows[(r["deck"], r["variant"])] = r
    return rows


def compare(da, db):
    A, B = load(da), load(db)
    keys = sorted(set(A) & set(B))
    flips_lost, flips_gain, stage_chg, ratios, acc = [], [], [], [], []
    for k in keys:
        a, b = A[k], B[k]
        ca, cb = int(a["converged"]), int(b["converged"])
        if ca and not cb:
            flips_lost.append(k)
        if cb and not ca:
            flips_gain.append(k)
        if ca and cb:
            if a["stage"] != b["stage"]:
                stage_chg.append((k, a["stage"], b["stage"]))
            if a["iters"] and b["iters"] and int(a["iters"]) > 0:
                ratios.append((int(b["iters"]) / int(a["iters"]), k))
            if a["flag"] != b["flag"]:
                acc.append((k, a["flag"], b["flag"], a["maxerr"], b["maxerr"]))
    only = sorted(set(A) ^ set(B))
    print(f"A={da}\nB={db}\ncommon runs: {len(keys)}  (only in one: {len(only)})")
    print(f"converged  A={sum(int(A[k]['converged']) for k in keys)}  B={sum(int(B[k]['converged']) for k in keys)}")
    print(f"flips converged->failed (A->B): {len(flips_lost)}")
    for k in flips_lost:
        print("   LOST ", "/".join(k), B[k]["status"])
    print(f"flips failed->converged (A->B): {len(flips_gain)}")
    for k in flips_gain:
        print("   GAINED", "/".join(k), "stage", B[k]["stage"])
    print(f"stage changes among both-converged: {len(stage_chg)}")
    for k, sa, sb in stage_chg:
        print(f"   {'/'.join(k)}: {sa} -> {sb}")
    if ratios:
        rs = [r for r, _ in ratios]
        gm = math.exp(sum(math.log(max(r, 1e-9)) for r in rs) / len(rs))
        print(f"iteration ratio B/A over {len(rs)} both-converged runs: geomean {gm:.3f}, median {statistics.median(rs):.3f}, "
              f"total A={sum(int(A[k]['iters']) for _, k in ratios)} B={sum(int(B[k]['iters']) for _, k in ratios)}")
        big = sorted(ratios, reverse=True)
        print("   worst B/A:", ", ".join(f"{'/'.join(k)}={r:.2f}" for r, k in big[:5] if r > 1.2) or "none > 1.2")
        print("   best  B/A:", ", ".join(f"{'/'.join(k)}={r:.2f}" for r, k in big[::-1][:5] if r < 0.8) or "none < 0.8")
    print(f"accuracy-flag changes among both-converged: {len(acc)}")
    for k, fa, fb, ea, eb in acc:
        print(f"   {'/'.join(k)}: {fa or '-'}({ea}) -> {fb or '-'}({eb})")
    print(f"WRONG  A={sum(A[k]['flag']=='WRONG' for k in keys)} B={sum(B[k]['flag']=='WRONG' for k in keys)};  "
          f"other  A={sum(A[k]['flag']=='other' for k in keys)} B={sum(B[k]['flag']=='other' for k in keys)}")
    wa = sum(float(A[k]['wall_s']) for k in keys)
    wb = sum(float(B[k]['wall_s']) for k in keys)
    print(f"wall total A={wa:.1f}s B={wb:.1f}s")
    return 0


def main(argv):
    if len(argv) > 1 and argv[1] == "compare":
        if len(argv) != 4:
            sys.exit("usage: opbench.py compare DIR_A DIR_B")
        return compare(argv[2], argv[3])
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ngspice", required=True)
    ap.add_argument("--env", action="append", metavar="K=V")
    ap.add_argument("--options", default="", help='e.g. "gminsteps=0 srcsteps=0"')
    ap.add_argument("--nodeset-seeds", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", default=None)
    ap.add_argument("--decks", default=None, help="comma-separated globs on deck names ('ref:*' for refs)")
    ap.add_argument("--corpus-dir", default=None, help="directory of *.cir decks (default tests/fixtures/op)")
    ap.add_argument("--no-refs", action="store_true", help="skip the refs/*.cir entries")
    ap.add_argument("--label", required=True)
    ap.add_argument("--seed", type=int, default=637)
    ap.add_argument("--refcache", default=None, help="directory for reference solutions (shared between labels of one binary)")
    ap.add_argument("--no-optran", action="store_true",
                    help="deselect the pseudo-transient (optran) fallback in the measured runs (not in the reference run)")
    ap.add_argument("--no-limit", action="store_true", help="do not wrap runs in systemd-run")
    a = ap.parse_args(argv[1:])
    if a.jobs > 4:
        sys.exit("--jobs must be <= 4 (project rule)")
    if a.out is None:
        a.out = os.path.expanduser(f"~/Work/data/ngspice_augmentation/637/{a.label}")
    return run_bench(a)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
