#!/usr/bin/env python3
"""tranbench -- transient-analysis robustness / cost / accuracy benchmark (card #641).

Transient analogue of tools/opbench.py.  For every deck of the corpus (tests/fixtures/tran/*.cir
plus a few converter fixtures, ids 'fx:<name>') it runs the deck's .tran with the binary under
test and records whether it completed, why not, how many accepted / rejected time points and
Newton iterations it needed, the wall time, how often the dev-build warning
'time point skipped' appeared on stderr, and how far the probed waveforms (and .meas values)
are from a REFERENCE run of the same deck by the same binary.  stdlib only.

  python3 tools/tranbench.py --ngspice BIN --label NAME [--env K=V ...] [--options "..."]
        [--decks a,b | 'fx:*' | 'tr_*'] [--jobs 4] [--refcache DIR] [--out DIR]
        [--timeout 300] [--grid 4000] [--no-fixtures] [--no-limit] [--list]
  python3 tools/tranbench.py compare DIR_A DIR_B

Deck conventions (corpus decks)
  * exactly one .tran card (a `tran` line in a .control block is accepted too and turned into
    a .tran card, as are `meas tran ...` lines; the rest of a .control block is dropped);
  * a header comment line   * probe: v(a) v(b) i(vsense) l1#branch   names the compared vectors;
  * optional   * refopts: itl4=100 ...   extra .options for the reference run only.
  .print/.plot/.probe/.save/.end cards are dropped; the harness adds its own .control block:
      set filetype=ascii / set numdgt=15 / run / rusage tranpoints trantime traniter rejected /
      wrdata <DIR/data/id.dat> <probes> / quit
  ASCII waveform data goes to OUT/data (or REFCACHE) only, never next to the deck.

Measured run   = deck as written + `.options <--options>` + the --env variables in the
                 environment (any inherited TR_* variable is removed first).
Reference run  = same binary, NO --env, `.options reltol=1e-6 abstol=1e-15 vntol=1e-9 trtol=1`,
                 .tran with tmax = the deck's own tmax if it has one, else tstep/4 (but never
                 fewer than span/50000 -> at most ~50000 forced points).  If that rung fails the
                 ladder tries (2) itl4=100 and tmax/4, (3) additionally method=gear, (4) reltol=1e-5
                 abstol=1e-12 vntol=1e-6, (5) reltol=1e-4 (same abstol/vntol); the rung used
                 is stored as ref_method.  Cache: REFCACHE/<id>.json (+ .dat), keyed on the
                 binary (path, size, mtime), the reference deck text, the reference settings and the
                 probes -- NEVER on --env / --options of the measured run.

Accuracy: both waveforms are linearly interpolated onto a uniform grid of --grid points over
[tstart, tstop] (only the part both runs cover is compared).  Per signal: max |err|, rms err;
normalised = err / max(peak-to-peak of the reference signal, 1e-3).  Deck figure = worst signal.
.meas values: relative error |a-b| / max(|ref|, 1e-12); a value missing in the run (failed) is
counted in meas_missing.  Statistics (median / max normalised error) are over decks that
completed and have a reference.
"""
import argparse
import bisect
import concurrent.futures as cf
import csv
import fnmatch
import glob
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join(ROOT, "tests", "fixtures", "tran")
FIXDIR = os.path.join(ROOT, "tests", "fixtures")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/641")
TIMEOUT = 300
GRID = 4000
NORM_FLOOR = 1e-3
MEAS_FLOOR = 1e-12
REF_OPTIONS = "reltol=1e-6 abstol=1e-15 vntol=1e-9 trtol=1"
REF_MAX_FORCED = 50000      # reference tmax is never smaller than span/REF_MAX_FORCED
REF_VERSION = "v2"
REF_METHODS = {1: "tight", 2: "tight+itl4+tmax/4", 3: "tight+gear", 4: "relaxed1e-5", 5: "relaxed1e-4"}

# existing converter decks used by path (their own .control block is not used: the tran / meas
# lines are lifted from it).  id -> (file, probes)
FIXTURES = {
    "fx:loss_diode_buck": ("loss_diode_buck.cir", "v(vout) l1#branch"),
    "fx:loss_mos1_buck": ("loss_mos1_buck.cir", "v(vout) l1#branch"),
    "fx:loss_bsim4_buck": ("loss_bsim4_buck.cir", "v(vout) l1#branch"),
    "fx:pss_r2_tran_settled": ("pss_r2_tran_settled.cir", "v(vx)"),
    "fx:stbtran_ldo_softstart": ("stbtran_ldo_softstart.cir", "v(vout) v(ea)"),
    "fx:convmeas_rlc_step": ("convmeas_rlc_step.cir", "v(vout)"),
}

CSV_FIELDS = ["deck", "completed", "fail_kind", "points", "rejected", "iters", "wall_s", "tran_s", "skipped",
              "reached", "nerr_max", "nerr_rms", "abs_max", "abs_rms", "worst_sig", "meas_err", "meas_missing",
              "ref_ok", "ref_method"]

SKIP_CARDS = re.compile(r"^\s*\.(print|plot|probe|save|end|tran)\b", re.I)
CTRL_START = re.compile(r"^\s*\.control\b", re.I)
CTRL_END = re.compile(r"^\s*\.endc\b", re.I)
SUFFIX = {"t": 1e12, "g": 1e9, "k": 1e3, "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}


# ----------------------------------------------------------------- numbers / deck parsing
def spice_num(tok):
    """'10n' -> 1e-8, '1meg' -> 1e6, plain floats; None when not a plain number."""
    m = re.match(r"^([-+]?(?:\d+\.?\d*|\.\d+)(?:e[-+]?\d+)?)([a-z]*)$", tok.strip().lower())
    if not m:
        return None
    x, suf = float(m.group(1)), m.group(2)
    if suf.startswith("meg"):
        return x * 1e6
    if suf.startswith("mil"):
        return x * 25.4e-6
    return x * SUFFIX.get(suf[:1], 1.0) if suf else x


def join_continuations(lines):
    out = []
    for ln in lines:
        if ln.startswith("+") and out:
            out[-1] += " " + ln[1:].strip()
        else:
            out.append(ln)
    return out


def params_of(lines):
    p = {}
    for ln in lines:
        m = re.match(r"^\s*\.param\s+(.*)$", ln, re.I)
        if m:
            for k, v in re.findall(r"(\w+)\s*=\s*(\{[^}]*\}|[^\s]+)", m.group(1)):
                p[k.lower()] = v
    return p


def eval_value(tok, params, depth=0):
    """Number from a token: plain number, .param name or {expression} of those (best effort)."""
    tok = tok.strip()
    v = spice_num(tok)
    if v is not None:
        return v
    if depth > 6:
        return None
    expr = tok.strip("{}").lower()

    def sub(m):
        w = m.group(0)
        if w in params:
            r = eval_value(params[w], params, depth + 1)
            return repr(r) if r is not None else w
        return w
    expr = re.sub(r"[a-z_]\w*", sub, expr)
    expr = re.sub(r"(\d(?:\.\d*)?(?:e[-+]?\d+)?)(meg|mil|[tgkmunpf])\b",
                  lambda m: repr(spice_num(m.group(1) + m.group(2))), expr)
    if not re.match(r"^[0-9eE+\-*/().\s]+$", expr):
        return None
    try:
        return float(eval(expr, {"__builtins__": {}}, {}))
    except Exception:
        return None


def parse_tran_card(text, params=None):
    """'.tran 1n 10u 0 5n uic' -> dict(tstep,tstop,tstart,tmax or None,uic).  Raises ValueError."""
    toks = text.replace(",", " ").split()
    if toks and toks[0].lower() in (".tran", "tran"):
        toks = toks[1:]
    uic = any(t.lower() == "uic" for t in toks)
    toks = [t for t in toks if t.lower() != "uic"]
    vals = [eval_value(t, params or {}) for t in toks]
    if len(vals) < 2 or vals[0] is None or vals[1] is None:
        raise ValueError(f"cannot read .tran card: {text!r}")
    return dict(tstep=vals[0], tstop=vals[1], tstart=(vals[2] or 0.0) if len(vals) > 2 else 0.0,
                tmax=vals[3] if len(vals) > 3 else None, uic=uic)


def parse_deck(text):
    """Split a deck into body lines (title.. without analysis/output cards and .control), the tran
    dict, the meas names, probes and refopts."""
    raw = join_continuations(text.splitlines())
    probes, refopts = None, ""
    for ln in raw:
        m = re.match(r"^\s*\*\s*probe:\s*(.*)$", ln, re.I)
        if m and probes is None:
            probes = m.group(1).strip()
        m = re.match(r"^\s*\*\s*refopts:\s*(.*)$", ln, re.I)
        if m:
            refopts = m.group(1).strip()
    params = params_of(raw)
    body, tran_line, ctrl_tran, ctrl_meas, in_ctrl = [], None, None, [], False
    first = True
    for ln in raw:
        if first:                      # title line, untouched
            body.append(ln)
            first = False
            continue
        if in_ctrl:
            if CTRL_END.match(ln):
                in_ctrl = False
            elif re.match(r"^\s*tran\s", ln, re.I) and ctrl_tran is None:
                ctrl_tran = ln
            elif re.match(r"^\s*meas(ure)?\s+tran\s", ln, re.I):
                ctrl_meas.append("." + ln.strip())
            continue
        if CTRL_START.match(ln):
            in_ctrl = True
            continue
        if re.match(r"^\s*\.tran\b", ln, re.I):
            if tran_line is None:
                tran_line = ln
            continue
        if SKIP_CARDS.match(ln):
            continue
        body.append(ln)
    if tran_line is None and ctrl_tran is not None:
        tran_line = ctrl_tran
        body += ctrl_meas
    if tran_line is None:
        raise ValueError("deck has no .tran card")
    tran = parse_tran_card(tran_line, params)
    meas_lines = [re.sub(r"^\s*\.", "", ln).strip() for ln in body if re.match(r"^\s*\.meas(?:ure)?\s+tran\s", ln, re.I)]
    body = [ln for ln in body if not re.match(r"^\s*\.meas(?:ure)?\s+tran\s", ln, re.I)]
    meas = [re.match(r"^meas(?:ure)?\s+tran\s+(\S+)", ln, re.I).group(1).lower() for ln in meas_lines]
    return dict(body=body, tran=tran, tran_line=tran_line, meas=meas, meas_lines=meas_lines, probes=probes, refopts=refopts)


def tran_card(t, tmax=None):
    s = f".tran {t['tstep']:.12g} {t['tstop']:.12g} {t['tstart']:.12g}"
    tm = t["tmax"] if tmax is None else tmax
    if tm is not None:
        s += f" {tm:.12g}"
    return s + (" uic" if t["uic"] else "")


def ref_tmax(t):
    """Reference max step: the deck's own tmax if it fixes one, else tstep/4 limited to span/REF_MAX_FORCED."""
    if t["tmax"] is not None:
        return t["tmax"]
    return max(t["tstep"] / 4.0, (t["tstop"] - t["tstart"]) / REF_MAX_FORCED)


def build_run_deck(d, options, datapath, probes, ref_rung=None):
    """Deck text.  ref_rung None = measured run; 1..3 = reference ladder rung."""
    lines = list(d["body"])
    t = d["tran"]
    if ref_rung is None:
        if options:
            lines.append(".options " + options)
        lines.append(tran_card(t))
    else:
        opts = REF_OPTIONS
        tm = ref_tmax(t)
        if ref_rung >= 2:
            opts += " itl4=100"
            tm = tm / 4.0 if t["tmax"] is None else tm
        if ref_rung == 3:
            opts += " method=gear"
        if ref_rung == 4:       # relaxed: some decks (floating gate drive, 1 uOhm switches) cannot run at 1e-6
            opts = "reltol=1e-5 abstol=1e-12 vntol=1e-6 trtol=1 itl4=100"
        if ref_rung == 5:
            opts = "reltol=1e-4 abstol=1e-12 vntol=1e-6 trtol=1 itl4=100"
        if d["refopts"]:
            opts += " " + d["refopts"]
        lines.append(".options " + opts)
        lines.append(tran_card(t, tm))
    lines += [".control", "set filetype=ascii", "set numdgt=15", "run",
              "rusage tranpoints trantime traniter rejected"]
    if probes:
        lines.append(f"wrdata {datapath} {probes}")
    # .meas cards are run as control-mode `meas` after the run so that `print` gives full precision
    lines += d["meas_lines"]
    for m in d["meas"]:
        lines.append(f"print {m}")
    lines += ["quit", ".endc", ".end"]
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------- log parsing
RX = {
    "points": re.compile(r"^Transient timepoints\s*=\s*(\d+)", re.M),
    "iters": re.compile(r"^Transient iterations\s*=\s*(\d+)", re.M),
    "rejected": re.compile(r"^Rejected timepoints\s*=\s*(\d+)", re.M),
    "tran_s": re.compile(r"^Transient analysis time\s*=\s*([-+0-9.eE]+)", re.M),
}
SKIPPED_RX = re.compile(r"time point skipped", re.I)
KIND_RX = [
    ("timestep too small", re.compile(r"timestep too small|time step too small", re.I)),
    ("no convergence", re.compile(r"no convergence|singular matrix|doAnalyses: .*(failed|error)|"
                                  r"gmin stepping failed|source stepping failed|iteration limit", re.I)),
]
OTHER_FATAL = re.compile(r"^\s*error|simulation\(s\) aborted|simulation interrupted|could not be simulated|"
                         r"unknown (model|subckt)|no such (vector|command)", re.I | re.M)


def parse_log(out, err):
    """-> dict(points,iters,rejected,tran_s,skipped,fail_kind)  (fail_kind '' when no error text)."""
    r = {}
    for k, rx in RX.items():
        m = rx.search(out)
        r[k] = (float(m.group(1)) if k == "tran_s" else int(m.group(1))) if m else None
    r["skipped"] = len(SKIPPED_RX.findall(err)) + len(SKIPPED_RX.findall(out))
    both = err + "\n" + out
    kind = ""
    for name, rx in KIND_RX:
        if rx.search(both):
            kind = name
            break
    if not kind and OTHER_FATAL.search(both):
        kind = "other"
    r["fail_kind"] = kind
    return r


def parse_meas(out, names):
    """Values of the .meas names from 'name = value' lines; None for 'failed' or absent."""
    vals = {}
    wanted = set(names)
    for ln in out.splitlines():
        m = re.match(r"^\s*([A-Za-z_][\w.]*)\s*=\s*(\S+)", ln)
        if not m or m.group(1).lower() not in wanted:
            continue
        k = m.group(1).lower()
        try:
            vals[k] = float(m.group(2))
        except ValueError:
            vals[k] = None
    return {k: vals.get(k) for k in names}


# ----------------------------------------------------------------- waveforms / error math
def load_wrdata(path, nsig):
    """ASCII wrdata file -> [(t[], y[]) per signal].  Columns: (time, value) pairs, or time + values."""
    cols = None
    try:
        f = open(path, errors="replace")
    except OSError:
        return None
    with f:
        for ln in f:
            parts = ln.split()
            if not parts:
                continue
            try:
                row = [float(x) for x in parts]
            except ValueError:
                continue
            if cols is None:
                cols = [[] for _ in row]
            if len(row) != len(cols):
                continue
            for c, x in zip(cols, row):
                c.append(x)
    if not cols or not cols[0]:
        return None
    if len(cols) >= 2 * nsig:
        return [(cols[2 * i], cols[2 * i + 1]) for i in range(nsig)]
    if len(cols) == nsig + 1:
        return [(cols[0], cols[i + 1]) for i in range(nsig)]
    return None


def interp(ts, ys, x):
    """Linear interpolation of the piecewise-linear waveform (ts ascending); clamps outside."""
    n = len(ts)
    if n == 1 or x <= ts[0]:
        return ys[0]
    if x >= ts[-1]:
        return ys[-1]
    j = bisect.bisect_right(ts, x)
    t0, t1 = ts[j - 1], ts[j]
    if t1 == t0:
        return ys[j]
    return ys[j - 1] + (ys[j] - ys[j - 1]) * (x - t0) / (t1 - t0)


def uniform_grid(t0, t1, n):
    if n < 2 or t1 <= t0:
        return [t0]
    return [t0 + (t1 - t0) * i / (n - 1) for i in range(n)]


def signal_error(run, ref, grid, floor=NORM_FLOOR):
    """run, ref = (t[], y[]).  Error over the grid points both waveforms cover.
    -> dict(abs_max, abs_rms, nerr_max, nerr_rms, pp, covered, n)."""
    rt, ry = run
    qt, qy = ref
    lim, lo = min(rt[-1], qt[-1]), max(rt[0], qt[0])
    eps = 1e-12 * max(abs(grid[-1]), 1e-30)
    sel = [g for g in grid if lo - eps <= g <= lim + eps]
    covered = len(sel) / len(grid)
    pp = max(qy) - min(qy)
    norm = max(pp, floor)
    if not sel:
        return dict(abs_max=None, abs_rms=None, nerr_max=None, nerr_rms=None, pp=pp, covered=0.0, n=0)
    mx, ss = 0.0, 0.0
    for g in sel:
        a, b = interp(rt, ry, g), interp(qt, qy, g)
        d = abs(a - b)
        if not math.isfinite(d):
            d = float("inf")
        mx = max(mx, d)
        ss += d * d
    rms = math.sqrt(ss / len(sel)) if math.isfinite(ss) else float("inf")
    return dict(abs_max=mx, abs_rms=rms, nerr_max=mx / norm, nerr_rms=rms / norm, pp=pp, covered=covered, n=len(sel))


def meas_errors(run, ref):
    """-> (max relative error or None, number missing, {name: rel err or None})."""
    per, worst, missing = {}, None, 0
    for k, rv in ref.items():
        if rv is None:
            continue
        a = run.get(k)
        if a is None or not math.isfinite(a):
            per[k] = None
            missing += 1
            continue
        e = abs(a - rv) / max(abs(rv), MEAS_FLOOR)
        per[k] = e
        worst = e if worst is None else max(worst, e)
    return worst, missing, per


def compare_waveforms(run_sigs, ref_sigs, names, tstart, tstop, ngrid=GRID):
    """-> (deck-level dict, per-signal dict)."""
    grid = uniform_grid(tstart, tstop, ngrid)
    per, worst = {}, None
    for i, nm in enumerate(names):
        if run_sigs is None or ref_sigs is None:
            per[nm] = None
            continue
        e = signal_error(run_sigs[i], ref_sigs[i], grid)
        per[nm] = e
        if e["nerr_max"] is not None and (worst is None or e["nerr_max"] > per[worst]["nerr_max"]):
            worst = nm
    deck = dict(nerr_max=None, nerr_rms=None, abs_max=None, abs_rms=None, worst_sig="", reached=None)
    if worst is not None:
        w = per[worst]
        deck.update(nerr_max=w["nerr_max"], abs_max=w["abs_max"], worst_sig=worst,
                    nerr_rms=max(e["nerr_rms"] for e in per.values() if e and e["nerr_rms"] is not None),
                    abs_rms=max(e["abs_rms"] for e in per.values() if e and e["abs_rms"] is not None),
                    reached=min(e["covered"] for e in per.values() if e))
    return deck, per


def reached_fraction(sigs, tstart, tstop):
    if not sigs:
        return 0.0
    tend = min(s[0][-1] for s in sigs)
    return max(0.0, min(1.0, (tend - tstart) / (tstop - tstart))) if tstop > tstart else 1.0


# ----------------------------------------------------------------- running
def clean_env(extra):
    env = {k: v for k, v in os.environ.items() if not k.startswith("TR_")}   # step-control switches only via --env
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_THREAD_LIMIT"] = "1"
    env.update(extra)
    return env


def launch(ngspice, deck_path, env_extra, limit, timeout):
    """-> (stdout, stderr, wall, error-or-None)."""
    cmd = [ngspice, "-b", deck_path]
    if limit and shutil.which("systemd-run"):
        cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=4G", "-p", "MemorySwapMax=0"] + cmd
    t0 = time.time()
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, errors="replace",
                             env=clean_env(env_extra), start_new_session=True, cwd=os.path.dirname(deck_path))
    except OSError as e:
        return "", "", time.time() - t0, f"cannot start: {e}"
    try:
        out, err = p.communicate(timeout=timeout)
        rc = p.returncode
        if rc < 0 or rc >= 128:
            return out, err, time.time() - t0, f"crash: exit status {rc}"
        return out, err, time.time() - t0, None
    except subprocess.TimeoutExpired:
        try:
            os.killpg(p.pid, 9)
        except OSError:
            pass
        out, err = p.communicate()
        return out or "", err or "", time.time() - t0, "timeout"


def clip(text, head=60000, tail=20000):
    """Keep logs small: a runaway warning stream can be millions of lines."""
    if len(text) <= head + tail:
        return text
    return text[:head] + f"\n[tranbench: {len(text) - head - tail} characters clipped]\n" + text[-tail:]


def one_run(ctx, deck_id, d, probes, options, tag, env, datadir, rundir, ref_rung=None):
    """Run a built deck once.  -> result dict (incl. completed/fail_kind/signals/meas)."""
    safe = deck_id.replace(":", "_")
    os.makedirs(datadir, exist_ok=True)
    os.makedirs(rundir, exist_ok=True)
    dat = os.path.join(datadir, f"{safe}{tag}.dat")
    open(dat, "w").close()      # truncate: a stale file can never be read back
    path = os.path.join(rundir, f"{safe}{tag}.cir")
    with open(path, "w") as f:
        f.write(build_run_deck(d, options, dat, probes, ref_rung))
    out, err, wall, lerr = launch(ctx["ngspice"], path, env, ctx["limit"], ctx["timeout"])
    r = parse_log(out, err)
    r.update(wall=wall, out=out, err=err, dat=dat)
    t = d["tran"]
    names = probes.split() if probes else []
    sigs = load_wrdata(dat, len(names)) if names else None
    r["sigs"] = sigs
    r["names"] = names
    r["meas"] = parse_meas(out, d["meas"])
    r["reached"] = reached_fraction(sigs, t["tstart"], t["tstop"]) if sigs else (1.0 if (r["points"] and not r["fail_kind"] and not lerr) else 0.0)
    kind = ""
    if lerr:
        kind = lerr.split(":")[0]
    elif r["fail_kind"]:
        kind = r["fail_kind"]
    elif r["points"] is None:
        kind = "other"
    elif r["reached"] < 1.0 - 1e-6:
        kind = "other"
    r["completed"] = (kind == "")
    r["fail_kind"] = kind
    if lerr:
        r["err"] += f"\n[tranbench] {lerr}\n"
    return r


def ref_key(ctx, d, probes):
    st = os.stat(ctx["ngspice"])
    blob = "\0".join([os.path.realpath(ctx["ngspice"]), str(st.st_size), str(int(st.st_mtime)), REF_VERSION,
                      REF_OPTIONS, str(REF_MAX_FORCED), "\n".join(d["body"]),
                      tran_card(d["tran"]), d["refopts"], probes or "", json.dumps({}, sort_keys=True)])
    return hashlib.sha1(blob.encode()).hexdigest()[:16]    # the measured run's --env / --options are NOT part of the key


def get_ref(ctx, deck_id, d, probes):
    """Cached reference.  -> dict(ok,method,key,names,dat,meas,points,wall,rung)  (waveforms are re-read from dat)."""
    key = ref_key(ctx, d, probes)
    cdir = ctx["refcache"]
    os.makedirs(cdir, exist_ok=True)
    safe = deck_id.replace(":", "_")
    cpath = os.path.join(cdir, safe + ".json")
    if os.path.exists(cpath):
        try:
            c = json.load(open(cpath))
            if c.get("key") == key and (not c["ok"] or os.path.exists(c["dat"])):
                c["cached"] = True
                return c
        except (OSError, ValueError, KeyError):
            pass
    c = dict(key=key, ok=False, method="", names=probes.split() if probes else [], dat="", meas={}, points=None,
             wall=0.0, kinds=[], cached=False)
    for rung in (1, 2, 3, 4, 5):
        r = one_run(ctx, deck_id, d, probes, None, f"__ref{rung}", {}, cdir, os.path.join(cdir, "rundecks"), ref_rung=rung)
        c["wall"] += r["wall"]
        with open(os.path.join(cdir, f"{safe}__ref{rung}.log"), "w") as f:
            f.write(clip(r["out"]) + "\n---- stderr ----\n" + clip(r["err"]))
        c["kinds"].append(r["fail_kind"])
        if r["completed"] and (r["sigs"] or not probes):
            c.update(ok=True, method=REF_METHODS[rung],
                     dat=r["dat"], meas=r["meas"], points=r["points"])
            break
        if r["fail_kind"] == "timeout":
            break
    with open(cpath, "w") as f:
        json.dump(c, f)
    return c


def job(ctx, deck_id, d, probes, ref):
    safe = deck_id.replace(":", "_")
    out = ctx["out"]
    r = one_run(ctx, deck_id, d, probes, ctx["options"], "", ctx["env"], os.path.join(out, "data"),
                os.path.join(out, "rundecks"))
    with open(os.path.join(out, "logs", safe + ".log"), "w") as f:
        f.write(clip(r["out"]) + "\n---- stderr ----\n" + clip(r["err"]))
    t = d["tran"]
    row = dict(deck=deck_id, completed=int(r["completed"]), fail_kind=r["fail_kind"],
               points=r["points"] if r["points"] is not None else "", rejected=r["rejected"] if r["rejected"] is not None else "",
               iters=r["iters"] if r["iters"] is not None else "", wall_s=f"{r['wall']:.3f}",
               tran_s=f"{r['tran_s']:.4f}" if r["tran_s"] is not None else "", skipped=r["skipped"],
               reached=f"{r['reached']:.4f}", nerr_max="", nerr_rms="", abs_max="", abs_rms="", worst_sig="",
               meas_err="", meas_missing="", ref_ok=int(bool(ref and ref["ok"])), ref_method=(ref or {}).get("method", ""))
    detail = dict(deck=deck_id, probes=probes, meas=r["meas"], signals={}, meas_err={})
    if ref and ref["ok"]:
        refsigs = load_wrdata(ref["dat"], len(ref["names"])) if ref["names"] else None
        if refsigs and r["sigs"]:
            deck_e, per = compare_waveforms(r["sigs"], refsigs, r["names"], t["tstart"], t["tstop"], ctx["grid"])
            for k in ("nerr_max", "nerr_rms", "abs_max", "abs_rms"):
                if deck_e[k] is not None:
                    row[k] = f"{deck_e[k]:.4e}"
            row["worst_sig"] = deck_e["worst_sig"]
            detail["signals"] = per
        if d["meas"]:
            worst, missing, per = meas_errors(r["meas"], ref["meas"])
            row["meas_err"] = "" if worst is None else f"{worst:.4e}"
            row["meas_missing"] = missing
            detail["meas_err"] = per
            detail["ref_meas"] = ref["meas"]
    with open(os.path.join(out, "results", safe + ".json"), "w") as f:
        json.dump(detail, f, indent=1)
    return row


# ----------------------------------------------------------------- corpus
def corpus(decks, use_fixtures, corpus_dir=None):
    """-> [(id, path, probes_override)]"""
    items = [(os.path.basename(p)[:-4], p, None) for p in sorted(glob.glob(os.path.join(corpus_dir or CORPUS, "*.cir")))]
    if use_fixtures:
        for k, (fn, pr) in FIXTURES.items():
            p = os.path.join(FIXDIR, fn)
            if os.path.exists(p):
                items.append((k, p, pr))
    if decks:
        pats = [x for x in decks.split(",") if x]
        items = [i for i in items if any(fnmatch.fnmatch(i[0], q) or i[0] == q for q in pats)]
    return items


def load_deck(path, probes_override):
    d = parse_deck(open(path, errors="replace").read())
    probes = probes_override or d["probes"]
    if probes is None:
        raise ValueError(f"{path}: no '* probe:' line")
    if os.path.dirname(os.path.abspath(path)) != os.path.abspath(CORPUS):
        base = os.path.dirname(os.path.abspath(path))   # fixture decks: make relative includes absolute
        d["body"] = [re.sub(r"^(\s*\.(?:include|inc|lib)\s+)(?!/)(\S+)", lambda m: m.group(1) + os.path.normpath(os.path.join(base, m.group(2))), ln, flags=re.I)
                     for ln in d["body"]]
    return d, probes


# ----------------------------------------------------------------- driver
def run_bench(a):
    out = os.path.abspath(a.out)
    for sub in ("logs", "data", "rundecks", "results"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    env = dict(kv.split("=", 1) for kv in (a.env or []))
    ctx = dict(ngspice=os.path.abspath(a.ngspice), env=env, out=out, options=(a.options or "").strip(),
               refcache=os.path.abspath(a.refcache) if a.refcache else os.path.join(os.path.dirname(out), "refcache"),
               limit=not a.no_limit, timeout=a.timeout, grid=a.grid)
    items = corpus(a.decks, not a.no_fixtures, a.corpus_dir)
    if a.list:
        for i in items:
            print(i[0], i[1])
        return 0
    if not items:
        sys.exit("no decks selected")
    decks, t_all = {}, time.time()
    for name, path, pr in items:
        decks[name] = load_deck(path, pr)
    refs = {}
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(get_ref, ctx, n, decks[n][0], decks[n][1]): n for n in decks}
        for f in cf.as_completed(futs):
            refs[futs[f]] = f.result()
    print(f"references: {sum(r['ok'] for r in refs.values())}/{len(refs)} ok "
          f"({sum(r['cached'] for r in refs.values())} cached)", flush=True)
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(job, ctx, n, decks[n][0], decks[n][1], refs[n]) for n in decks]
        rows = [f.result() for f in futs]
    rows.sort(key=lambda r: r["deck"])
    with open(os.path.join(out, "results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)
    summ = summarize(rows)
    summ.update(label=a.label, ngspice=a.ngspice, options=ctx["options"], env=env, grid=a.grid, timeout=a.timeout,
                refcache=ctx["refcache"], elapsed_s=round(time.time() - t_all, 1))
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print_summary(summ)
    return 0


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def summarize(rows):
    done = [r for r in rows if int(r["completed"])]
    nerrs = [fnum(r["nerr_max"]) for r in done if fnum(r["nerr_max"]) is not None]
    rmss = [fnum(r["nerr_rms"]) for r in done if fnum(r["nerr_rms"]) is not None]
    merrs = [fnum(r["meas_err"]) for r in done if fnum(r["meas_err"]) is not None]

    def tot(k):
        return sum(int(r[k]) for r in rows if r[k] != "")
    kinds = {}
    for r in rows:
        if not int(r["completed"]):
            kinds[r["fail_kind"]] = kinds.get(r["fail_kind"], 0) + 1
    return dict(decks=len(rows), completed=len(done), completed_pct=round(100.0 * len(done) / max(1, len(rows)), 1),
                points=tot("points"), rejected=tot("rejected"), iters=tot("iters"), skipped=tot("skipped"),
                wall_s=round(sum(float(r["wall_s"]) for r in rows), 2),
                nerr_median=statistics.median(nerrs) if nerrs else None, nerr_max=max(nerrs) if nerrs else None,
                nerr_rms_median=statistics.median(rmss) if rmss else None,
                meas_err_max=max(merrs) if merrs else None, nerr_n=len(nerrs),
                ref_methods={m: sum(1 for r in rows if r["ref_method"] == m) for m in sorted({r["ref_method"] for r in rows if r["ref_method"]})},
                no_ref=sorted(r["deck"] for r in rows if not int(r["ref_ok"])),
                fail_kinds=kinds,
                failed=sorted(f"{r['deck']} ({r['fail_kind']})" for r in rows if not int(r["completed"])))


def _g(x):
    return "n/a" if x is None else f"{x:.3g}"


def print_summary(s):
    print(f"\n== {s.get('label', '')}  decks={s['decks']}  elapsed={s.get('elapsed_s', '?')}s ==")
    print(f"decks {s['decks']}  completed {s['completed']} ({s['completed_pct']}%)  points {s['points']}  rejected {s['rejected']}  "
          f"iterations {s['iters']}  skipped-warnings {s['skipped']}  wall {s['wall_s']} s")
    print(f"normalised error (completed decks with reference, n={s['nerr_n']}): median {_g(s['nerr_median'])}  "
          f"max {_g(s['nerr_max'])}  (rms median {_g(s['nerr_rms_median'])});  .meas rel. error max {_g(s['meas_err_max'])}")
    print("failure kinds:", s["fail_kinds"] or "none")
    print("failed decks:", ", ".join(s["failed"]) if s["failed"] else "none")
    if s.get("ref_methods"):
        print("reference rungs used:", s["ref_methods"])
    if s["no_ref"]:
        print("decks without a reference:", ", ".join(s["no_ref"]))


# ----------------------------------------------------------------- compare
def load(d):
    with open(os.path.join(d, "results.csv")) as f:
        return {r["deck"]: r for r in csv.DictReader(f)}


def _ratio(a, b):
    a, b = fnum(a), fnum(b)
    return None if not a or b is None else b / a


def compare(da, db, out=None):
    out = out or sys.stdout
    A, B = load(da), load(db)
    keys = sorted(set(A) & set(B))
    P = lambda *x: print(*x, file=out)
    P(f"A={da}\nB={db}\ncommon decks: {len(keys)}  (only in one: {len(set(A) ^ set(B))})")
    ca = [k for k in keys if int(A[k]["completed"])]
    cb = [k for k in keys if int(B[k]["completed"])]
    both = [k for k in keys if int(A[k]["completed"]) and int(B[k]["completed"])]
    P(f"completed  A={len(ca)}  B={len(cb)}  both={len(both)}")
    P(f"{'deck':32} {'done':>7} {'points A':>9} {'B':>9} {'rej A':>7} {'B':>7} {'skip A':>6} {'B':>6} {'iter B/A':>8} {'wall B/A':>8} {'nerr A':>9} {'B':>9}  note")
    lost, gained, worse, better = [], [], [], []
    for k in keys:
        a, b = A[k], B[k]
        da_, db_ = int(a["completed"]), int(b["completed"])
        note = []
        if da_ and not db_:
            lost.append(k)
            note.append("LOST:" + b["fail_kind"])
        if db_ and not da_:
            gained.append(k)
            note.append("GAINED(was " + a["fail_kind"] + ")")
        if not da_ and not db_ and a["fail_kind"] != b["fail_kind"]:
            note.append(f"{a['fail_kind']}->{b['fail_kind']}")
        ea, eb = fnum(a["nerr_max"]), fnum(b["nerr_max"])
        if da_ and db_ and ea is not None and eb is not None:
            if eb > 2 * ea and eb > 1e-4:
                worse.append(k)
                note.append("ERR UP")
            elif ea > 2 * eb and ea > 1e-4:
                better.append(k)
                note.append("err down")
        ir = _ratio(a["iters"], b["iters"]) if da_ and db_ else None
        wr = _ratio(a["wall_s"], b["wall_s"]) if da_ and db_ else None
        P(f"{k:32} {('Y' if da_ else 'n') + '>' + ('Y' if db_ else 'n'):>7} {a['points']!s:>9} {b['points']!s:>9} "
          f"{a['rejected']!s:>7} {b['rejected']!s:>7} {a['skipped']!s:>6} {b['skipped']!s:>6} {'-' if ir is None else format(ir, '.2f'):>8} "
          f"{'-' if wr is None else format(wr, '.2f'):>8} {a['nerr_max'] or '-':>9} {b['nerr_max'] or '-':>9}  {' '.join(note)}")

    def tot(rs, ks, f):
        return sum(int(rs[k][f]) for k in ks if rs[k][f] != "")
    P(f"\nflips completed->failed (A->B): {len(lost)} {lost}")
    P(f"flips failed->completed (A->B): {len(gained)} {gained}")
    P(f"error up >2x: {len(worse)} {worse};  error down >2x: {len(better)} {better}")
    if both:
        P(f"over the {len(both)} decks completed by both:  points A={tot(A, both, 'points')} B={tot(B, both, 'points')}  "
          f"rejected A={tot(A, both, 'rejected')} B={tot(B, both, 'rejected')}  "
          f"iterations A={tot(A, both, 'iters')} B={tot(B, both, 'iters')}  "
          f"wall A={sum(float(A[k]['wall_s']) for k in both):.1f}s B={sum(float(B[k]['wall_s']) for k in both):.1f}s")
        rs = [r for r in (_ratio(A[k]["iters"], B[k]["iters"]) for k in both) if r]
        if rs:
            gm = math.exp(sum(math.log(max(r, 1e-9)) for r in rs) / len(rs))
            P(f"iteration ratio B/A: geomean {gm:.3f}, median {statistics.median(rs):.3f}")
        for lab, R in (("A", A), ("B", B)):
            e = [fnum(R[k]["nerr_max"]) for k in both if fnum(R[k]["nerr_max"]) is not None]
            if e:
                P(f"normalised error {lab}: median {statistics.median(e):.3g} max {max(e):.3g}")
    return 0


def main(argv):
    if len(argv) > 1 and argv[1] == "compare":
        if len(argv) != 4:
            sys.exit("usage: tranbench.py compare DIR_A DIR_B")
        return compare(argv[2], argv[3])
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ngspice", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--env", action="append", metavar="K=V", help="environment variable for the measured runs only (repeatable)")
    ap.add_argument("--options", default="", help='extra .options for the measured runs, e.g. "reltol=1e-2"')
    ap.add_argument("--decks", default=None, help="comma-separated names / globs ('fx:*' for fixtures)")
    ap.add_argument("--corpus-dir", default=None, help="directory of *.cir decks (default tests/fixtures/tran)")
    ap.add_argument("--no-fixtures", action="store_true", help="skip the fx:* fixture decks")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", default=None)
    ap.add_argument("--refcache", default=None, help="directory for reference runs (shared between labels of one binary)")
    ap.add_argument("--timeout", type=float, default=TIMEOUT)
    ap.add_argument("--grid", type=int, default=GRID, help="points of the uniform comparison grid")
    ap.add_argument("--no-limit", action="store_true", help="do not wrap runs in systemd-run")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args(argv[1:])
    if a.jobs > 4:
        sys.exit("--jobs must be <= 4 (project rule)")
    if a.out is None:
        a.out = os.path.join(DATA_ROOT, a.label)
    for kv in a.env or []:
        if "=" not in kv:
            sys.exit(f"--env expects K=V, got {kv!r}")
    return run_bench(a)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
