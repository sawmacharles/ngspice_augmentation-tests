#!/usr/bin/env python3
"""pstbsens.py -- parameter sensitivities of the periodic loop gain (.pstb) by central finite
differences over the complete .pss + .pstb run.  Card #475, periodic half (2026-10-14).

    tools/pstbsens.py DECK --param cl_val [--param M2:w ...] [--h 1e-2] [--abs-step A]
                      [--jobs 4] [--orbit NODE] [--out DIR]

For each parameter p the deck is run at p(1+h), p(1-h), p(1+h/2), p(1-h/2) (relative steps;
with --abs-step A the steps are A, A/2 -- required when p = 0), plus once at nominal.
    D(s)  = (f(p+s) - f(p-s)) / (2 s)
    D     = (4 D(h/2) - D(h)) / 3          Richardson: the s^2 term cancels
    e     = |D(h/2) - D(h)|                error estimate, per output
Outputs: complex dT/dp(f) at every sweep frequency (T = the `t` vector .pstb writes), dPM/dp, dfc/dp,
dGM/dp (margins computed from T with the definition of pacshoot.c, see margins()), normalised
sensitivities p/PM dPM/dp, and with --orbit NODE d(mean)/dp and d(peak-to-peak)/dp of v(NODE)
over one PSS period.  Why finite differences and not an adjoint: docs/pstb-sens.md.

Parameter names (as sens= in docs/sens.md):  a .param name; an element name for its principal value
(R, C, L, V/I dc); <inst>_<kw> or <inst>:<kw> for an instance parameter (M2:w), and <inst>:<kw> for a
parameter of the instance's model (M2:vto; a model shared by several instances moves all of them).

Exit status: 0 = ran, 2 = refused or could not run (message on stderr).  Outputs go to --out
(default ~/Work/data/ngspice_augmentation/pstbsens/<deck>-<time>/), never the repository.
"""
import argparse
import cmath
import concurrent.futures
import csv
import datetime
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/pstbsens")
DEFAULT_BIN = os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice")


def ngspice_bin():
    env = os.environ.get("NGSPICE")
    if env:
        return env
    return DEFAULT_BIN if os.path.exists(DEFAULT_BIN) else "ngspice"


class Refused(RuntimeError):
    """The request or a run cannot be honoured; reported, exit 2."""


# ---------------------------------------------------------------------------------------------
# SPICE numbers and a small expression evaluator (only to learn the nominal value of p)
# ---------------------------------------------------------------------------------------------
SUFFIX = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "mil": 25.4e-6, "m": 1e-3, "u": 1e-6,
          "n": 1e-9, "p": 1e-12, "f": 1e-15, "a": 1e-18}
NUM_RE = r"(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"


def parse_number(s):
    m = re.fullmatch(r"([+-]?" + NUM_RE + r")([a-zA-Z]*)", s.strip())
    if not m:
        raise ValueError(f"not a SPICE number: {s!r}")
    val, suf = float(m.group(1)), m.group(2).lower()
    if not suf:
        return val
    if suf.startswith("meg"):
        return val * 1e6
    if suf.startswith("mil"):
        return val * 25.4e-6
    if suf[0] in SUFFIX:
        return val * SUFFIX[suf[0]]
    return val  # unit letters (V, ohm, F ...) as ngspice does


def eval_expr(expr, params, _depth=0):
    """Evaluate a SPICE expression with .param names; raises ValueError if it cannot."""
    if _depth > 20:
        raise ValueError("recursive .param")
    e = expr.strip().strip("{}'\"")

    def num(m):
        return repr(parse_number(m.group(0)))
    e = re.sub(r"(?<![\w.])" + NUM_RE + r"(?:meg|mil|[tgkmunpfa])?[a-zA-Z]*(?![\w.(])", num, e, flags=re.I)
    e = e.replace("^", "**")
    env = {k: getattr(math, k) for k in ("sqrt", "exp", "log", "log10", "sin", "cos", "tan", "pow", "pi")}
    env.update({"ln": math.log, "abs": abs, "min": min, "max": max})
    for name in set(re.findall(r"(?<![\w.])[A-Za-z_]\w*(?!\w*\()", e)):
        if name in env:
            continue
        raw = params.get(name.lower())
        if raw is None:
            raise ValueError(f"unknown name {name!r}")
        env[name] = eval_expr(raw, params, _depth + 1)
    try:
        return float(eval(e, {"__builtins__": {}}, env))  # noqa: S307 -- own deck text, no builtins
    except Exception as exc:
        raise ValueError(f"cannot evaluate {expr!r}: {exc}") from exc


# ---------------------------------------------------------------------------------------------
# deck handling
# ---------------------------------------------------------------------------------------------
VAL = r"(\{[^}]*\}|'[^']*'|[^\s)]+)"


def load_deck(path):
    """Logical lines (continuations joined), absolute include paths, control blocks and .end removed.
    Returns (lines, removed_control_lines)."""
    base = os.path.dirname(os.path.abspath(path))
    with open(path) as fh:
        raw = fh.read().splitlines()
    lines = []
    for n, ln in enumerate(raw):
        if n and ln.startswith("+") and lines:
            lines[-1] += " " + ln[1:].strip()
        else:
            lines.append(ln)
    out, removed, in_ctl = [], [], False
    for n, ln in enumerate(lines):
        low = ln.strip().lower()
        if n and low.startswith(".control"):
            in_ctl = True
            continue
        if in_ctl:
            if low.startswith(".endc"):
                in_ctl = False
            else:
                removed.append(ln.strip())
            continue
        if n and low == ".end":
            continue
        m = re.match(r"^(\s*\.(?:include|inc|lib)\s+)(['\"]?)([^'\"\s]+)(['\"]?)(.*)$", ln, re.I)
        if n and m and not os.path.isabs(m.group(3)):
            ln = f"{m.group(1)}{m.group(2)}{os.path.join(base, m.group(3))}{m.group(4)}{m.group(5)}"
        out.append(ln)
    return out, removed


def find_pstb(lines):
    for n, ln in enumerate(lines):
        if n and re.match(r"^\s*\.pstb\b", ln, re.I):
            return n
    return None


def param_table(lines):
    tab = {}
    for n, ln in enumerate(lines):
        if n and re.match(r"^\s*\.param\b", ln, re.I):
            for m in re.finditer(r"(?<![\w.])([A-Za-z_]\w*)\s*=\s*" + VAL, ln[len(ln) - len(ln.lstrip()) + 6:]):
                tab[m.group(1).lower()] = m.group(2)
    return tab


class Target:
    """One perturbable number in the deck: where it is, its nominal value and raw text."""

    def __init__(self, name, line_no, span, raw, p0, kind, note=""):
        self.name, self.line_no, self.span, self.raw, self.p0, self.kind, self.note = \
            name, line_no, span, raw, p0, kind, note

    def new_text(self, delta):
        """Replacement for the raw value token with p0 + delta."""
        if self.p0 is not None and not re.match(r"^[{'\"]", self.raw):
            return f"{self.p0 + delta:.17g}"
        inner = self.raw.strip("{}'\"")
        return "{(" + inner + ")+(" + f"{delta:.17g}" + ")}"

    def apply(self, lines, delta):
        out = list(lines)
        ln = out[self.line_no]
        a, b = self.span
        out[self.line_no] = ln[:a] + self.new_text(delta) + ln[b:]
        return out


def _value_target(name, lines, n, regex, params, kind, note=""):
    hits = list(re.finditer(regex, lines[n]))
    if len(hits) != 1:
        return None
    raw = hits[0].group(1)
    try:
        p0 = eval_expr(raw, params)
    except ValueError as e:
        raise Refused(f"parameter {name}: cannot learn its nominal value from {raw!r} ({e})")
    return Target(name, n, hits[0].span(1), raw, p0, kind, note)


def resolve(spec, lines):
    """Find the number the user means. Order: .param, element principal value, inst:kw / inst_kw."""
    params = param_table(lines)
    key = spec.lower()
    elems = {}
    models = {}
    for n, ln in enumerate(lines):
        if not n or not ln.strip() or ln.lstrip()[0] in "*.":
            m = re.match(r"^\s*\.model\s+(\S+)", ln, re.I)
            if n and m:
                models[m.group(1).lower()] = n
            continue
        elems.setdefault(ln.split()[0].lower(), n)

    def instparam(inst, kw):
        n = elems.get(inst)
        if n is None:
            return None
        return _value_target(spec, lines, n, r"(?i)(?<![\w.])" + re.escape(kw) + r"\s*=\s*" + VAL,
                             params, "instance parameter")

    def modparam(inst, kw):
        n = elems.get(inst)
        if n is None:
            return None
        toks = lines[n].split()
        mname = next((t.lower() for t in toks[2:] if t.lower() in models), None)
        if mname is None:
            return None
        mn = models[mname]
        t = _value_target(spec, lines, mn, r"(?i)(?<![\w.])" + re.escape(kw) + r"\s*=\s*" + VAL,
                          params, "model parameter")
        if t is not None:
            users = [x for x, nn in elems.items()
                     if mname in [q.lower() for q in lines[nn].split()[2:]]]
            if len(users) > 1:
                t.note = f"model {mname} is shared by {len(users)} instances ({', '.join(sorted(users))}): all move"
        return t

    if ":" in key:
        inst, kw = key.split(":", 1)
        t = instparam(inst, kw) or modparam(inst, kw)
        if t:
            return t
        raise Refused(f"unknown parameter {spec!r}: no instance parameter or model parameter {kw!r} of {inst!r}")
    if key in params:
        hits = [n for n, ln in enumerate(lines) if n and re.match(r"^\s*\.param\b", ln, re.I)
                and re.search(r"(?i)(?<![\w.])" + re.escape(key) + r"\s*=", ln)]
        if len(hits) != 1:
            raise Refused(f"parameter {spec!r} is defined {len(hits)} times in .param lines: ambiguous")
        n = hits[0]
        regex = r"(?i)(?<![\w.])" + re.escape(key) + r"\s*=\s*" + VAL
        t = _value_target(spec, lines, n, regex, params, ".param")
        if t is None:
            raise Refused(f"parameter {spec!r} appears twice on its .param line: ambiguous")
        return t
    if key in elems:
        n = elems[key]
        toks = lines[n].split()
        kind = key[0]
        if kind in "rcl" and len(toks) > 3 and "=" not in toks[3]:
            start = len(lines[n]) - len(lines[n].lstrip())
            pos = 0
            for _ in range(3):
                m = re.compile(r"\s*\S+").match(lines[n], pos if pos else start)
                pos = m.end()
            m = re.compile(r"\s*" + VAL).match(lines[n], pos)
            raw = m.group(1)
            try:
                p0 = eval_expr(raw, params)
            except ValueError as e:
                raise Refused(f"parameter {spec}: cannot learn the nominal value of {raw!r} ({e})")
            return Target(spec, n, m.span(1), raw, p0, "element value")
        if kind in "vi":
            m = re.search(r"(?i)\bdc\s+" + VAL, lines[n])
            if m:
                return Target(spec, n, m.span(1), m.group(1), eval_expr(m.group(1), params), "element dc value")
        raise Refused(f"parameter {spec!r}: element {toks[0]} has no principal value this tool can rewrite "
                      f"(use {toks[0]}:<keyword>)")
    if "_" in key:
        parts = key.split("_")
        for i in range(1, len(parts)):
            inst, kw = "_".join(parts[:i]), "_".join(parts[i:])
            if inst in elems:
                t = instparam(inst, kw) or modparam(inst, kw)
                if t:
                    return t
    raise Refused(f"unknown parameter {spec!r}: not a .param, an element, or <inst>:<keyword> of this deck")


# ---------------------------------------------------------------------------------------------
# one run
# ---------------------------------------------------------------------------------------------
def control_block(outdir, orbit):
    c = ["", ".control", "set wr_singlescale", "set wr_vecnames", "set numdgt=15", "run",
         f"wrdata {outdir}/t.txt t"]
    if orbit:
        c += ["foreach p $plots", "setplot $p", f"wrdata {outdir}/orb.txt v({orbit})", "end"]
    c += ["quit", ".endc", ".end", ""]
    return "\n".join(c)


def read_complex(path):
    """(freqs, [complex]) from wrdata of one complex vector with wr_singlescale + wr_vecnames."""
    f, t = [], []
    with open(path) as fh:
        next(fh)
        for row in fh:
            v = [float(x) for x in row.split()]
            if len(v) != 3:
                raise Refused(f"unexpected column count in {path}")
            f.append(v[0])
            t.append(complex(v[1], v[2]))
    if not f:
        raise Refused(f"{path} is empty")
    return f, t


def margins(freqs, T):
    """PM, fc, GM, fGM of T with the definition of pacshoot.c's .pstb (verbatim logic): phase unwrapped
    from the first point, |T| in dB; the 0 dB crossing and the -180 deg crossings are found by linear
    interpolation in log f (frequency = f0 (f1/f0)^t of the interpolated quantity); PM = 180 + phase there
    (phase linear in the dB fraction), GM = -|T|dB at a -180 deg crossing; the lowest of several
    crossings is reported.  Returns dict with None for a missing crossing."""
    n = len(freqs)
    ph, prev = [], 0.0
    for i, t in enumerate(T):
        p = math.degrees(math.atan2(t.imag, t.real))
        if i:
            while p - prev > 180.0:
                p -= 360.0
            while p - prev < -180.0:
                p += 360.0
        ph.append(p)
        prev = p
    db = [20.0 * math.log10(abs(t)) for t in T]
    pms, gms = [], []
    for i in range(n - 1):
        f0, f1, m0, m1, p0, p1 = freqs[i], freqs[i + 1], db[i], db[i + 1], ph[i], ph[i + 1]
        if m0 * m1 < 0 or (m1 == 0.0 and m0 != 0.0):
            t = (0.0 - m0) / (m1 - m0)
            pms.append((180.0 + p0 + t * (p1 - p0), f0 * (f1 / f0) ** t))
        lo, hi = min(p0, p1), max(p0, p1)
        k = math.ceil((lo + 180.0) / 360.0)
        while -180.0 + 360.0 * k <= hi:
            target = -180.0 + 360.0 * k
            if target != p0 and p1 != p0:
                t = (target - p0) / (p1 - p0)
                gms.append((-(m0 + t * (m1 - m0)), f0 * (f1 / f0) ** t))
            k += 1
    pm = min(pms) if pms else None
    gm = min(gms) if gms else None
    return {"PM": pm[0] if pm else None, "fc": pm[1] if pm else None,
            "GM": gm[0] if gm else None, "fGM": gm[1] if gm else None, "npm": len(pms)}


def _lagrange(xs, ys, x):
    tot = 0.0
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        w = yi
        for j, xj in enumerate(xs):
            if j != i:
                w *= (x - xj) / (xi - xj)
        tot += w
    return tot


def margins_smooth(freqs, T):
    """Same definition as margins() (0 dB crossing of |T|, -180 deg crossings, PM = 180 + phase there, lowest of
    several), but |T| in dB and the unwrapped phase are interpolated by the cubic through the four nearest
    sweep points (in ln f) instead of the chord.  The chord has a kink every time a crossing moves into the next
    grid cell, which makes a finite difference over a step that spans cells wrong by the interpolation error
    (measured: dPM/dIb 0.5 % off at 200 points/decade).  Falls back to the chord when the cubic has no root in
    the cell."""
    n = len(freqs)
    if n < 4:
        return margins(freqs, T)
    x = [math.log(f) for f in freqs]
    ph, prev = [], 0.0
    for i, t in enumerate(T):
        p = math.degrees(math.atan2(t.imag, t.real))
        if i:
            while p - prev > 180.0:
                p -= 360.0
            while p - prev < -180.0:
                p += 360.0
        ph.append(p)
        prev = p
    db = [20.0 * math.log10(abs(t)) for t in T]

    def root(i, ys, target):
        k = min(max(i - 1, 0), n - 4)
        xs4, y4 = x[k:k + 4], ys[k:k + 4]
        a, b = x[i], x[i + 1]
        fa, fb = _lagrange(xs4, y4, a) - target, _lagrange(xs4, y4, b) - target
        if fa * fb > 0:
            return None
        for _ in range(100):
            m = 0.5 * (a + b)
            fm = _lagrange(xs4, y4, m) - target
            if fa * fm <= 0:
                b, fb = m, fm
            else:
                a, fa = m, fm
            if b - a < 1e-15:
                break
        return 0.5 * (a + b), xs4, k

    pms, gms = [], []
    for i in range(n - 1):
        m0, m1, p0, p1 = db[i], db[i + 1], ph[i], ph[i + 1]
        if m0 * m1 < 0 or (m1 == 0.0 and m0 != 0.0):
            r = root(i, db, 0.0)
            if r:
                xr, xs4, k = r
                pms.append((180.0 + _lagrange(xs4, ph[k:k + 4], xr), math.exp(xr)))
            else:
                t = (0.0 - m0) / (m1 - m0)
                pms.append((180.0 + p0 + t * (p1 - p0), freqs[i] * (freqs[i + 1] / freqs[i]) ** t))
        lo, hi = min(p0, p1), max(p0, p1)
        kk = math.ceil((lo + 180.0) / 360.0)
        while -180.0 + 360.0 * kk <= hi:
            target = -180.0 + 360.0 * kk
            if target != p0 and p1 != p0:
                r = root(i, ph, target)
                if r:
                    xr, xs4, k = r
                    gms.append((-_lagrange(xs4, db[k:k + 4], xr), math.exp(xr)))
                else:
                    t = (target - p0) / (p1 - p0)
                    gms.append((-(m0 + t * (m1 - m0)), freqs[i] * (freqs[i + 1] / freqs[i]) ** t))
            kk += 1
    pm = min(pms) if pms else None
    gm = min(gms) if gms else None
    return {"PM": pm[0] if pm else None, "fc": pm[1] if pm else None,
            "GM": gm[0] if gm else None, "fGM": gm[1] if gm else None, "npm": len(pms)}


def run_one(cfg, key, deck_lines):
    """Run one deck variant. Returns dict(freqs, T, margins, line_pm, line_fc, orbit, wall)."""
    outdir = os.path.join(cfg["workdir"], key)
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, "deck.cir")
    with open(deck, "w") as fh:
        fh.write("\n".join(deck_lines) + "\n" + control_block(outdir, cfg["orbit"]))
    cmd = [cfg["bin"], "-b", deck]
    if cfg["systemd"]:
        cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=4G", "-p", "MemorySwapMax=0"] + cmd
    env = dict(os.environ, OMP_NUM_THREADS="2")
    if cfg["pss_reltol"]:
        env["PSS_SHOOT_RELTOL"] = cfg["pss_reltol"]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=cfg["timeout"], cwd=outdir, env=env)
    except subprocess.TimeoutExpired:
        raise Refused(f"run {key}: timeout after {cfg['timeout']} s")
    except OSError as e:
        raise Refused(f"cannot start {cmd[0]}: {e}")
    wall = time.time() - t0
    log = proc.stdout + proc.stderr
    with open(os.path.join(outdir, "log.txt"), "w") as fh:
        fh.write(log)
    if "PSS-SHOOT: CONVERGED" not in log:
        raise Refused(f"run {key}: PSS did not converge (no 'PSS-SHOOT: CONVERGED' line; log {outdir}/log.txt)")
    pm_line = re.search(r"pstb: phase margin ([-+0-9.eE]+) deg at ([-+0-9.eE]+) Hz", log)
    if not pm_line:
        why = re.search(r"pstb: phase margin: no 0 dB crossing", log)
        raise Refused(f"run {key}: the pstb phase-margin line is missing"
                      + (" (no 0 dB crossing of T in the swept band)" if why else f" (log {outdir}/log.txt)"))
    tfile = os.path.join(outdir, "t.txt")
    if not os.path.exists(tfile):
        raise Refused(f"run {key}: ngspice wrote no T (exit {proc.returncode}; log {outdir}/log.txt)")
    freqs, T = read_complex(tfile)
    mg = margins(freqs, T)
    ms = margins_smooth(freqs, T)
    lpm, lfc = float(pm_line.group(1)), float(pm_line.group(2))
    if mg["PM"] is None or abs(mg["PM"] - lpm) > 1e-3 * max(abs(lpm), 1.0) \
            or abs(mg["fc"] - lfc) > 1e-5 * lfc:
        raise Refused(f"run {key}: PM computed from T ({mg['PM']}, {mg['fc']}) differs from the pstb line "
                      f"({lpm}, {lfc})")
    cross = re.search(r"phase margin .*\((\d+) crossings found", log)
    res = {"freqs": freqs, "T": T, "m": ms if cfg["interp"] == "cubic" else mg, "m_chord": mg,
           "m_alt": mg if cfg["interp"] == "cubic" else ms, "line_pm": lpm, "line_fc": lfc, "wall": wall,
           "ncross": int(cross.group(1)) if cross else 1, "orbit": None}
    if cfg["orbit"]:
        of = os.path.join(outdir, "orb.txt")
        if os.path.exists(of):
            tt, vv = [], []
            with open(of) as fh:
                next(fh)
                for row in fh:
                    a = row.split()
                    tt.append(float(a[0]))
                    vv.append(float(a[1]))
            if len(tt) > 3:
                span = tt[-1]
                # periodic rectangle/trapezoid: v(0) = v(T); mean over (0, T]
                integ = 0.0
                prev_t, prev_v = 0.0, vv[-1]
                for a, b in zip(tt, vv):
                    integ += 0.5 * (prev_v + b) * (a - prev_t)
                    prev_t, prev_v = a, b
                res["orbit"] = {"mean": integ / span, "pp": max(vv) - min(vv), "period": span, "n": len(tt)}
    return res


# ---------------------------------------------------------------------------------------------
# differences
# ---------------------------------------------------------------------------------------------
def combine(vals, s, mode):
    """vals: dict '+h','-h','+h2','-h2','0' -> value (float/complex/None). Returns (D, e) or (None, None)."""
    if mode == "forward":
        if vals["+h"] is None or vals["0"] is None:
            return None, None
        return (vals["+h"] - vals["0"]) / s, 0.0
    need = ("+h", "-h", "+h2", "-h2")
    if any(vals[k] is None for k in need):
        return None, None
    d1 = (vals["+h"] - vals["-h"]) / (2 * s)
    d2 = (vals["+h2"] - vals["-h2"]) / s
    return (4 * d2 - d1) / 3, abs(d2 - d1)


def sensitivities(tgt, s, runs, mode):
    """runs: dict label -> run result. Returns the result record of one parameter."""
    labels = ["0", "+h"] if mode == "forward" else ["0", "+h", "-h", "+h2", "-h2"]
    pick = lambda f: {k: (f(runs[k]) if k in runs else None) for k in ["0", "+h", "-h", "+h2", "-h2"]}  # noqa: E731
    rec = {"name": tgt.name, "kind": tgt.kind, "value": tgt.p0, "step": s, "note": tgt.note}
    nom = runs["0"]
    nf = len(nom["freqs"])
    dT, eT = [], []
    for i in range(nf):
        D, e = combine(pick(lambda r: r["T"][i]), s, mode)
        dT.append(D)
        eT.append(e)
    rec["dT"], rec["dT_err"] = dT, eT
    for key, out in (("PM", "dPM"), ("fc", "dfc"), ("GM", "dGM"), ("fGM", "dfGM")):
        D, e = combine(pick(lambda r: r["m"][key]), s, mode)
        rec[out], rec[out + "_err"] = D, e
    for key, out in (("PM", "dPM"), ("fc", "dfc")):
        D, e = combine(pick(lambda r: r["m_alt"][key]), s, mode)
        rec[out + "_alt"], rec[out + "_alt_err"] = D, e
    rec["PM"], rec["fc"] = nom["m"]["PM"], nom["m"]["fc"]
    p0 = tgt.p0
    rec["nPM"] = p0 / rec["PM"] * rec["dPM"] if rec["dPM"] is not None and p0 != 0 else None
    rec["nPM_err"] = p0 / rec["PM"] * rec["dPM_err"] if rec["dPM"] is not None and p0 != 0 else None
    rec["nfc"] = p0 / rec["fc"] * rec["dfc"] if rec["dfc"] is not None and p0 != 0 else None
    if all(runs[k]["orbit"] for k in labels):
        for key in ("mean", "pp"):
            D, e = combine(pick(lambda r: r["orbit"][key]), s, mode)
            rec["d" + key], rec["d" + key + "_err"] = D, e
    return rec


def fmt(x, w=11):
    return f"{'n/a':>{w}}" if x is None else f"{x:>{w}.4e}"


def print_table(recs, mode, orbit):
    print(f"{'param':<14}{'value':>11}{'dPM/dp':>12}{'+-err':>11}{'p/PM dPM/dp':>13}{'dfc/dp':>12}{'+-err':>11}{'dGM/dp':>12}")
    for r in recs:
        print(f"{r['name']:<14}{fmt(r['value'])}{fmt(r['dPM'], 12)}{fmt(r['dPM_err'])}{fmt(r['nPM'], 13)}"
              f"{fmt(r['dfc'], 12)}{fmt(r['dfc_err'])}{fmt(r['dGM'], 12)}")
    print(f"(PM = {recs[0]['PM']:.6g} deg at fc = {recs[0]['fc']:.6g} Hz, nominal; "
          f"{'one-sided forward difference [RED PROOF MODE]' if mode == 'forward' else 'central difference + Richardson'})")
    if orbit:
        print(f"\n{'param':<14}{'d mean/dp':>12}{'+-err':>11}{'d p-p/dp':>12}{'+-err':>11}   (v({orbit}) over one PSS period)")
        for r in recs:
            print(f"{r['name']:<14}{fmt(r.get('dmean'), 12)}{fmt(r.get('dmean_err'))}{fmt(r.get('dpp'), 12)}{fmt(r.get('dpp_err'))}")
    for r in recs:
        if r["note"]:
            print(f"note: {r['name']}: {r['note']}")


def zlist(z, part):
    return None if z is None else (z.real if part == "re" else z.imag)


def write_outputs(outdir, freqs, nomT, recs, meta):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "dtdp.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        head = ["freq", "T_re", "T_im"]
        for r in recs:
            head += [f"dT_d{r['name']}_re", f"dT_d{r['name']}_im", f"dT_d{r['name']}_err"]
        w.writerow(head)
        for i, f in enumerate(freqs):
            row = [repr(f), repr(nomT[i].real), repr(nomT[i].imag)]
            for r in recs:
                d, e = r["dT"][i], r["dT_err"][i]
                row += [repr(d.real), repr(d.imag), repr(e)] if d is not None else ["", "", ""]
            w.writerow(row)
    js = {"meta": meta, "params": []}
    for r in recs:
        q = {k: v for k, v in r.items() if k not in ("dT", "dT_err")}
        js["params"].append(q)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(js, fh, indent=1, sort_keys=True)
        fh.write("\n")


def sens_run(deck, params, h=1e-2, abs_step=None, jobs=4, orbit=None, outdir=None, forward=False,
             timeout=3600, systemd=True, pss_reltol=None, interp="chord"):
    """The tool as a function. Returns (records, freqs, nominalT, outdir, timing)."""
    if not 1 <= jobs <= 5:
        raise Refused("--jobs must be between 1 and 5 (machine budget: 5 ngspice processes)")
    if not params:
        raise Refused("give at least one --param")
    lines, removed = load_deck(deck)
    if find_pstb(lines) is None:
        raise Refused("the deck has no .pstb card (nothing to differentiate)")
    mode = "forward" if forward else "richardson"
    if forward and h == 1e-2:
        h = 0.1
    targets = [resolve(p, lines) for p in params]
    steps = []
    for t in targets:
        if t.p0 == 0:
            if abs_step is None:
                raise Refused(f"parameter {t.name} is 0: give --abs-step")
            steps.append(abs_step)
        else:
            steps.append(h * t.p0)
    stem = os.path.splitext(os.path.basename(deck))[0]
    outdir = outdir or os.path.join(DATA_ROOT, f"{stem}-{datetime.datetime.now():%Y%m%d-%H%M%S}")
    outdir = os.path.abspath(outdir)
    workdir = os.path.join(outdir, "runs")
    if os.path.isdir(workdir):
        for sub in os.listdir(workdir):
            shutil.rmtree(os.path.join(workdir, sub))
    cfg = {"bin": ngspice_bin(), "workdir": workdir, "orbit": orbit, "timeout": timeout,
           "systemd": systemd and shutil.which("systemd-run") is not None,
           "pss_reltol": None if pss_reltol in (None, "0") else str(pss_reltol), "interp": interp}
    todo = [("nom", "0", None, lines)]
    for i, (t, s) in enumerate(zip(targets, steps)):
        variants = [("+h", s), ("-h", -s), ("+h2", s / 2), ("-h2", -s / 2)]
        if forward:
            variants = variants[:1]
        for lab, d in variants:
            todo.append((f"p{i}_{lab}", lab, i, t.apply(lines, d)))
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {k: ex.submit(run_one, cfg, k, dl) for k, lab, i, dl in todo}
        results = {}
        err = None
        for k, f in futs.items():
            try:
                results[k] = f.result()
            except Refused as e:
                err = err or e
        if err:
            raise err
    wall = time.time() - t0
    nom = results["nom"]
    recs = []
    for i, t in enumerate(targets):
        runs = {"0": nom}
        for k, lab, ii, _ in todo:
            if ii == i:
                runs[lab] = results[k]
        if runs["+h"]["T"] == nom["T"]:
            raise Refused(f"parameter {t.name!r} had no effect (T bit-identical to nominal at +{h:g}): "
                          f"not found, or not used by the loop/orbit")
        if any(r["ncross"] != nom["ncross"] for r in runs.values()):
            print(f"warning: {t.name}: the number of 0 dB crossings changes between runs; margin derivatives "
                  f"are not smooth there", file=sys.stderr)
        recs.append(sensitivities(t, steps[i], runs, mode))
    if orbit and nom["orbit"] is None:
        print(f"warning: --orbit {orbit}: no v({orbit}) vector found in the PSS plot; orbit outputs skipped",
              file=sys.stderr)
    meta = {"deck": os.path.abspath(deck), "h": h, "abs_step": abs_step, "mode": mode, "orbit": orbit,
            "nominal_PM": nom["m"]["PM"], "nominal_fc": nom["m"]["fc"], "nominal_GM": nom["m"]["GM"],
            "runs": len(todo), "control_block_replaced": bool(removed), "pss_reltol": cfg["pss_reltol"],
            "margin_interp": interp, "nominal_PM_chord": nom["m_chord"]["PM"], "nominal_fc_chord": nom["m_chord"]["fc"]}
    write_outputs(outdir, nom["freqs"], nom["T"], recs, meta)
    timing = {"wall_total_s": wall, "runs": len(todo), "jobs": jobs,
              "mean_run_s": sum(r["wall"] for r in results.values()) / len(results)}
    with open(os.path.join(outdir, "timing.json"), "w") as fh:
        json.dump(timing, fh, indent=1)
        fh.write("\n")
    return recs, nom, outdir, timing, mode, removed


def main(argv=None):
    ap = argparse.ArgumentParser(description="Parameter sensitivities of the periodic loop gain (.pstb) "
                                 "by central finite differences + Richardson over the full .pss+.pstb run.")
    ap.add_argument("deck")
    ap.add_argument("--param", action="append", default=[], metavar="P",
                    help=".param name, element name, <inst>:<kw> / <inst>_<kw> (repeatable)")
    ap.add_argument("--h", type=float, default=1e-2, help="relative step (default 1e-2)")
    ap.add_argument("--abs-step", type=float, default=None, help="absolute step, needed when p = 0")
    ap.add_argument("--jobs", type=int, default=4, help="parallel ngspice processes (1..5, default 4)")
    ap.add_argument("--orbit", metavar="NODE", help="also d(mean)/dp, d(peak-to-peak)/dp of v(NODE) over one PSS period")
    ap.add_argument("--out", metavar="DIR", help=f"output directory (default {DATA_ROOT}/<deck>-<time>)")
    ap.add_argument("--timeout", type=int, default=3600, help="per-run timeout in seconds")
    ap.add_argument("--pss-reltol", default=None, metavar="X",
                    help="set PSS_SHOOT_RELTOL=X for every run (engine default 1e-5). On a switching deck 1e-8 cuts the "
                         "run-to-run noise of T (grid selection, ~1e-3) by ~1000x, which a 1%% parameter step needs; "
                         "not the default because it makes a linear deck fail (timestep too small)")
    ap.add_argument("--margin-interp", choices=("cubic", "chord"), default="chord",
                    help="interpolation of T between sweep points when locating the crossings for the margin "
                         "derivatives (default cubic; chord = the pstb line's own, kinked in p)")
    ap.add_argument("--no-systemd", action="store_true", help="do not wrap ngspice in systemd-run --scope")
    ap.add_argument("--redproof-forward", action="store_true", help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    try:
        recs, nom, outdir, timing, mode, removed = sens_run(
            a.deck, a.param, h=a.h, abs_step=a.abs_step, jobs=a.jobs, orbit=a.orbit, outdir=a.out,
            forward=a.redproof_forward, timeout=a.timeout, systemd=not a.no_systemd,
            pss_reltol=a.pss_reltol, interp=a.margin_interp)
    except Refused as e:
        print(f"pstbsens: refused: {e}", file=sys.stderr)
        return 2
    print_table(recs, mode, a.orbit)
    if removed:
        print(f"note: the deck's own .control block ({len(removed)} lines) was replaced by run + wrdata t")
    print(f"\nwrote {outdir}/dtdp.csv, summary.json, timing.json; {timing['runs']} runs, "
          f"{timing['wall_total_s']:.1f} s wall, jobs={timing['jobs']} (mean run {timing['mean_run_s']:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
