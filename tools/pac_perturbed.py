#!/usr/bin/env python3
"""tools/pac_perturbed.py -- correct-by-construction small-signal AC for
switching converters built with IDEAL switches (sw/csw, XSPICE PWM), card
#494 parallel track B.

WHY THIS EXISTS: native `.pac solver=td` linearises the circuit per
accepted transient step, which is exactly right for a SMOOTH switch model
(a continuous B-source conductance) but gives an ideal `sw`/`csw`/XSPICE
switch's own control-to-power-port gain as *identically zero* in every step
except the ones where the switch actually flips -- the whole PWM modulator
gain lives in a saltation (event) correction this project's fork does not
yet have right (`docs/pss-pac.md`'s "Saltation" section, card #494 main
track: 13.19 V/V vs a ~6 V/V reference, not resolved). This tool sidesteps
the problem entirely by never linearising anything: it runs the REAL
nonlinear circuit (ideal switches included) to two separate converged
periodic steady states (baseline, and the same circuit with a small
sinusoidal perturbation added at the injection point), subtracts them, and
reads the sideband-0 (same-frequency) component of the difference by a
single-bin DFT -- the "perturbed periodic steady state" method
`tools/pssac_ref.py` (card #484, `research/near-fs2.md`) already validated
against `.ac`, an exact analytic LPTV reference, and (for R6o at 1kHz) the
saltation investigation's own two other independent references. This file
reuses `pssac_ref.py`'s numerical core (rational frequency matching, DFT,
resampling, orbit-plausibility check, `CouldNotRun`) but is a GENERIC,
user-facing tool: any deck, an injection point chosen from the command
line (not a hardcoded per-fixture TAGS table), any output(s), and a
frequency list -- plus an optional Tian two-injection loop-gain probe
(reusing `tools/stbcore.py`'s `tian()`) and parallel dispatch of the many
independent ngspice runs a sweep needs via `tools/campaign.py`'s
memory-aware scheduler.

METHOD (see tools/pssac_ref.py's own docstring for the full derivation;
summarised here):

  1. A target frequency f must be commensurate with the deck's switching
     frequency `fs` (--fs) for a finite LCM-period orbit to exist:
     f_actual = fs*p/q, found by `rational_approx()` (imported from
     pssac_ref.py, unmodified).
  2. Two orbits are brought to their own periodic steady state over the
     LCM period T=q/fs: baseline (injection source at its plain DC value,
     no perturbation) and perturbed (the SAME deck with a small
     SIN(offset A f_actual 0 0 0) added at the injection point). Backend:
     `.pss method=shoot` (fast, one LCM-period solve) when it converges
     and produces a plausible orbit (q small enough); a long, settled
     `.tran` with a periodicity-residual check otherwise (automatic
     fallback -- q>1 shooting-PSS emission is known-buggy on this fork,
     research/near-fs2.md's "Engine finding", and this tool's whole point
     is ideal-switch decks, which are also the ones least likely to
     shoot-converge at all -- docs/pss-pac.md's saltation section reports
     R6/R6s closed-loop `.pss method=shoot` non-convergence and an XSPICE
     `d_pwm` segfault in the SAME shooting engine's period-restart, for
     reasons unrelated to saltation; the long-`.tran` fallback here is a
     PLAIN, non-restarting transient and does not hit either bug).
  3. orbit_diff(t) = perturbed(t) - baseline(t), single-bin DFT at bin p
     (Goertzel-style, `dft_bin()`), divided by the same DFT of the
     analytically-known injected waveform -- a phasor gain, sidesteps
     phase-convention bugs.
  4. Linearity: repeat at amplitude/2 (Amplitude halving), report the
     dB/deg spread (default target <=0.05dB, informational, not a hard
     gate -- see tools/pssac_ref.py's own honesty note on this: a couple
     of validated points there did not fully reach 0.05dB either).
  5. Loop-gain mode (--loop-gain): a flat Tian two-injection probe (same
     three-device topology as stb/lgprobe.lib: Ii (shunt current into an
     internal node x), Vi (series voltage, x to the downstream node a),
     Vb (x to the upstream node b), instantiated directly in the deck
     text rather than via the shared .lib -- this tool never touches or
     .includes stb/lgprobe.lib) -- two runs (Vi perturbed / Ii perturbed),
     sharing one baseline (both DC 0), sideband-0 v(x) and i(Vi) from
     each, combined via `tools/stbcore.py`'s `tian()`, EXACTLY the
     formula `tools/stb.py`'s native-.ac loop-gain tool already uses,
     just fed perturbed-PSS sideband-0 numbers instead of linear-AC ones.

BACKEND SELECTION AND PARALLEL DISPATCH: every ngspice invocation a sweep
needs (one baseline per distinct q, one full-amplitude and one
half-amplitude perturbed run per requested frequency, doubled again for
loop-gain mode's two injections) is built as a `tools/campaign.py`
"point" (a self-contained deck string) and dispatched through
`campaign.run_campaign_a()`'s memory-aware, one-process-per-point
scheduler -- the SAME scheduler this project's own PVT/MC campaigns use.
Backend='auto' runs a first PSS-attempt BATCH (all frequencies in
parallel), then a SECOND batch, tran-fallback only, for whichever
frequencies' PSS attempt failed to converge or produced an implausible
(corrupted-emission) orbit -- two batches, not one ngspice call blocking
on another's fallback decision, so the parallelism is not serialised by
the auto-fallback logic.

Usage:
    python3 tools/pac_perturbed.py sweep \\
        --deck tests/fixtures/pac_r3_ctrl2out.cir --fs 500e3 \\
        --inject source:Vctrl --outputs vout \\
        --freqs 1000,20000,199054 --amplitude 0.01 \\
        --ngspice ~/Work/apps/ngspice_augmentation/aaba50ba3f49-s41/bin/ngspice-aug \\
        --out ~/Work/data/ngspice_augmentation/pac_perturbed/run1 --max-procs 4

    python3 tools/pac_perturbed.py sweep --deck refs/r6_vm_buck_ideal.cir --fs 500e3 \\
        --inject node:opout_pac_unused --loop-gain --probe-nodes fbtop:vout \\
        --sweep-spec "dec 5 500 200000" --backend tran ...

    NGSPICE=<prefix>/bin/ngspice-aug python3 tools/pac_perturbed.py sweep ...

Exit status: 0 if the sweep ran (individual points may still be reported
COULD-NOT-RUN in the CSV/JSON -- see the "status" column/field), 2 if the
whole invocation could not run at all (bad deck, bad injection spec, no
commensurate frequency found for ANY requested point).
"""
import argparse
import cmath
import csv
import json
import math
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import pssac_ref as ref  # noqa: E402  -- reused, not modified (numeric core)
import stbcore  # noqa: E402 -- reused, not modified (tian())
import campaign  # noqa: E402 -- reused, not modified (run_campaign_a scheduler)

CouldNotRun = ref.CouldNotRun

RUN_ID = os.environ.get("PACP_RUN_ID") or f"pid{os.getpid()}_{int(time.time())}"
DATA_ROOT = os.path.expanduser(f"~/Work/data/ngspice_augmentation/pac_perturbed/{RUN_ID}")


# --------------------------------------------------------------------- #
# Small parsing helpers (deck-agnostic; no TAGS table, unlike pssac_ref.py)
# --------------------------------------------------------------------- #

_SI = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "mil": 25.4e-6, "m": 1e-3,
       "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}
_SI_RE = re.compile(r"^([+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?)([a-zA-Z]*)$")


def parse_num(s):
    """A float, or a SPICE-suffixed number (500k, 1meg, 10u, ...)."""
    s = s.strip()
    m = _SI_RE.match(s)
    if not m:
        raise ValueError(f"not a number: {s!r}")
    val = float(m.group(1))
    suf = m.group(2).lower()
    if not suf:
        return val
    if suf in _SI:
        return val * _SI[suf]
    # "meg" is matched whole above; a bare leading "m" followed by other
    # letters (e.g. stray unit text) still resolves via the 1-char table
    for k in ("meg", "mil", "t", "g", "k", "m", "u", "n", "p", "f"):
        if suf.startswith(k):
            return val * _SI[k]
    raise ValueError(f"unknown suffix on {s!r}")


def parse_freqs(freqs_arg, sweep_spec_arg):
    if freqs_arg:
        return [parse_num(x) for x in freqs_arg.split(",") if x.strip()]
    if sweep_spec_arg:
        toks = sweep_spec_arg.split()
        if len(toks) != 4:
            raise SystemExit(f"--sweep-spec expects '{{lin|dec|oct}} N fstart fstop', got {sweep_spec_arg!r}")
        kind, n, fstart, fstop = toks[0].lower(), int(toks[1]), parse_num(toks[2]), parse_num(toks[3])
        out = []
        if kind == "lin":
            if n <= 1:
                out = [fstart]
            else:
                for k in range(n):
                    out.append(fstart + (fstop - fstart) * k / (n - 1))
        elif kind in ("dec", "oct"):
            ratio = 10.0 ** (1.0 / n) if kind == "dec" else 2.0 ** (1.0 / n)
            f = fstart
            while f <= fstop * 1.0000001:
                out.append(f)
                f *= ratio
        else:
            raise SystemExit(f"unknown sweep kind {kind!r} (want lin/dec/oct)")
        return out
    raise SystemExit("one of --freqs or --sweep-spec is required")


def parse_outputs(outputs_arg):
    out = []
    for tok in outputs_arg.split(","):
        tok = tok.strip()
        if not tok:
            continue
        out.append(tok if "(" in tok else f"v({tok})")
    if not out:
        raise SystemExit("--outputs gave no expressions")
    return out


def guess_oscnode(out_exprs, fallback="0"):
    m = re.match(r"v\(([^)]+)\)", out_exprs[0], re.I)
    return m.group(1) if m else fallback


def parse_injection(spec):
    parts = spec.split(":")
    kind = parts[0].lower()
    if kind == "source":
        if len(parts) != 2:
            raise SystemExit(f"--inject source:<InstName>, got {spec!r}")
        return {"kind": "source", "name": parts[1]}
    if kind == "node":
        if len(parts) not in (2, 3):
            raise SystemExit(f"--inject node:<node>[:<ref>], got {spec!r}")
        return {"kind": "node", "node": parts[1], "ref": parts[2] if len(parts) == 3 else "0"}
    if kind == "series":
        if len(parts) != 3:
            raise SystemExit(f"--inject series:<a>:<b>, got {spec!r}")
        return {"kind": "series", "a": parts[1], "b": parts[2]}
    raise SystemExit(f"--inject must start with source:/node:/series:, got {spec!r}")


# card #541 (bug 2 of #539): a PULSE whose TR+PW+TF exceeds PER (R3/R4/R7 ramps: tr=0.97*T, tf=0.03*T, pw=1n -> sum T+1n)
# makes the next rising edge start before the previous falling edge ends. ngspice (stock 45.2 too) then sets breakpoints
# 1 ns apart at every period boundary and aborts `.tran` "Timestep too small ... trouble with node vramp#branch" for many
# tight option sets (reltol 1e-6/1e-7, vntol 1e-7, ...: deterministic, at a multiple of the period). The waveform is the
# same with TF shortened by the overlap, so the tool does that for the runs it makes (the user's deck file is untouched).
PULSE_FIXES = []   # one note per source changed in the last load_deck_body() call
_PULSE_RE = re.compile(r"^(\s*[vi]\S*\s+\S+\s+\S+\s+(?:dc\s+\S+\s+)?pulse\s*\()([^)]*)(\).*)$", re.I)


def _eval_param_expr(txt, env):
    import ast
    txt = txt.strip()
    if txt.startswith("{") and txt.endswith("}"):
        txt = txt[1:-1]
    try:
        return parse_num(txt)
    except ValueError:
        pass
    def ev(n):
        if isinstance(n, ast.Expression):
            return ev(n.body)
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
            return float(n.value)
        if isinstance(n, ast.Name):
            return env[n.id.lower()]
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            v = ev(n.operand)
            return -v if isinstance(n.op, ast.USub) else v
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            ops = {ast.Add: lambda: a + b, ast.Sub: lambda: a - b, ast.Mult: lambda: a * b,
                   ast.Div: lambda: a / b, ast.Pow: lambda: a ** b}
            return ops[type(n.op)]()
        raise ValueError("unsupported")
    # SPICE suffixes inside expressions (1n, 2u): turn a number+suffix token into its value
    txt = re.sub(r"(?<![\w.])(\d+\.?\d*(?:e[+-]?\d+)?)(meg|mil|[tgkmunpf])(?![\w])",
                 lambda m: repr(parse_num(m.group(1) + m.group(2))), txt, flags=re.I)
    return ev(ast.parse(txt, mode="eval"))


def fix_pulse_overlap(body):
    """Shorten TF of every PULSE with TR+PW+TF > PER by the excess (see PULSE_FIXES above). Skips anything it cannot
    evaluate (unknown parameter, expression it does not understand): the line is then left as the deck wrote it."""
    env = {}
    for line in body:
        if re.match(r"\s*\.param\b", line, re.I):
            for m in re.finditer(r"(\w+)\s*=\s*(\{[^}]*\}|[^\s]+)", line.split(None, 1)[1] if len(line.split(None, 1)) > 1 else ""):
                try:
                    env[m.group(1).lower()] = _eval_param_expr(m.group(2), env)
                except Exception:  # noqa: BLE001
                    pass
    out = []
    for line in body:
        m = _PULSE_RE.match(line)
        if m:
            args = m.group(2).split()
            try:
                if len(args) >= 7:
                    tr, tf, pw, per = (_eval_param_expr(args[i], env) for i in (3, 4, 5, 6))
                    ex = tr + pw + tf - per
                    if per > 0 and ex > 1e-12 * per and ex < tf:
                        args[4] = repr(tf - ex)
                        PULSE_FIXES.append(f"{line.split()[0]}: PULSE TR+PW+TF exceeds PER by {ex:.3g} s; TF shortened to {tf - ex:.6g}")
                        line = m.group(1) + " ".join(args) + m.group(3)
            except Exception:  # noqa: BLE001
                pass
        out.append(line)
    return out


def load_deck_body(deck_path, params=None):
    """Body of `deck_path` up to .control/.end, with optional .param
    substitutions (reuses pssac_ref._substitute_param, unmodified)."""
    body = ref._load_body(deck_path)
    for name, value in (params or {}):
        body = ref._substitute_param(body, name, value, deck_path)
    # card #536 defect 1: a fixture's own analysis card (`.pac`, `.pss`, `.tran`
    # ...) left in the body competes with the `run`/`pss` this tool appends (a
    # converging `.pac` consumes `run`, no `.tran` ever runs, wrdata writes the
    # .pac card's own orbit, H=0 -> math domain error). Strip them ALL, here,
    # so every backend (pss, tran, freerun, loop-gain) is covered, and record
    # what was ignored so the summary can say so.
    body, stripped = _strip_dot_cards(body, ANALYSIS_CARDS, report=True)
    STRIPPED_CARDS[:] = stripped
    PULSE_FIXES[:] = []
    fixed = fix_pulse_overlap(body)
    # Opt-in (PACP_FIX_PULSE=1): shortening TF moves the falling edge, which is
    # not neutral on ideal-switch decks (handbook example 10 / pac_r6o: the
    # reference run stopped settling, 2026-10-02). Default: only warn.
    if os.environ.get("PACP_FIX_PULSE") == "1":
        body = fixed
    else:
        PULSE_FIXES[:] = [n.replace("TF shortened to", "NOT changed (set PACP_FIX_PULSE=1 to shorten TF to") + ")" for n in PULSE_FIXES]
    return body


_SOURCE_LINE_TMPL = r"^{name}\s+(\S+)\s+(\S+)\s+DC\s+(\S+)(?:\s+AC\s+\S+)?\s*$"


def find_source_line(body, name):
    rx = re.compile(_SOURCE_LINE_TMPL.format(name=re.escape(name)), re.I)
    hits = [(i, rx.match(l.strip())) for i, l in enumerate(body) if rx.match(l.strip())]
    if len(hits) != 1:
        raise CouldNotRun(f"source:{name}: found {len(hits)} matching independent-source lines "
                           f"(expected exactly 1; only a plain 'Name n+ n- DC <val>[ AC <val>]' "
                           f"line is supported)")
    idx, m = hits[0]
    return idx, m.group(1), m.group(2), m.group(3)


def _device_name_present(body, name):
    for l in body:
        f = l.strip().split()
        if f and f[0].lower() == name.lower():
            return True
    return False


def apply_injection(body, injection, amplitude, f_actual):
    """A NEW list (does not mutate body): the injection device set to its
    plain DC value (amplitude==0, baseline) or DC + SIN(...) (perturbed)."""
    body = list(body)
    kind = injection["kind"]
    if kind == "source":
        idx, n1, n2, dcexpr = find_source_line(body, injection["name"])
        if amplitude == 0.0:
            body[idx] = f"{injection['name']} {n1} {n2} DC {dcexpr}"
        else:
            body[idx] = (f"{injection['name']} {n1} {n2} DC {dcexpr} "
                         f"SIN({dcexpr} {amplitude:.10g} {f_actual:.10g} 0 0 0)")
        return body
    if kind == "node":
        name = "Iinject_pac"
        if _device_name_present(body, name):
            raise CouldNotRun(f"deck already has a device named {name}; --inject node: cannot reuse it")
        if amplitude == 0.0:
            body.append(f"{name} {injection['node']} {injection['ref']} DC 0")
        else:
            body.append(f"{name} {injection['node']} {injection['ref']} DC 0 "
                         f"SIN(0 {amplitude:.10g} {f_actual:.10g} 0 0 0)")
        return body
    if kind == "series":
        name = "Vinject_pac"
        if _device_name_present(body, name):
            raise CouldNotRun(f"deck already has a device named {name}; --inject series: cannot reuse it")
        if amplitude == 0.0:
            body.append(f"{name} {injection['a']} {injection['b']} DC 0")
        else:
            body.append(f"{name} {injection['a']} {injection['b']} DC 0 "
                         f"SIN(0 {amplitude:.10g} {f_actual:.10g} 0 0 0)")
        return body
    raise CouldNotRun(f"unknown injection kind {kind!r}")


LOOP_PROBE_NODE = "xprobe_pac"


def apply_loop_injection(body, a, b, which, amplitude, f_actual):
    """which in ('none','v','i'): the flat Tian probe (Ii shunt-into-x,
    Vi series x->a, Vb x->b -- same topology as stb/lgprobe.lib, built
    inline so this tool never .includes or edits that shared file)."""
    body = list(body)
    x = LOOP_PROBE_NODE
    if _device_name_present(body, "Ii_pac") or _device_name_present(body, "Vi_pac"):
        raise CouldNotRun("deck already has Ii_pac/Vi_pac devices; --loop-gain cannot reuse them")
    ii_dc = f"Ii_pac 0 {x} DC 0"
    vi_dc = f"Vi_pac {x} {a} DC 0"
    if which == "v" and amplitude != 0.0:
        vi = f"{vi_dc} SIN(0 {amplitude:.10g} {f_actual:.10g} 0 0 0)"
        ii = ii_dc
    elif which == "i" and amplitude != 0.0:
        vi = vi_dc
        ii = f"{ii_dc} SIN(0 {amplitude:.10g} {f_actual:.10g} 0 0 0)"
    else:
        vi, ii = vi_dc, ii_dc
    body.append(ii)
    body.append(vi)
    body.append(f"Vb_pac {b} {x} DC 0")
    return body


# --------------------------------------------------------------------- #
# Deck-TEXT builders (return a string, do NOT write files -- campaign.py's
# run_campaign_a owns deck-file writing; see build_point_deck below).
# --------------------------------------------------------------------- #

def pss_deck_text(body, f_lcm, stabtime, oscnode, points, out_expr, orbit_name,
                   sc_iter=200, steady_coeff=1e-7, harmonics=10):
    lines = list(body)
    lines.append(".control")
    lines.append("set filetype=ascii")
    lines.append("set wr_singlescale")
    lines.append("set numdgt=17")
    lines.append(f"pss {f_lcm!r} {stabtime} {oscnode} {points} {harmonics} {sc_iter} "
                 f"{steady_coeff} uic method=shoot")
    lines.append(f"wrdata {orbit_name} {out_expr}")
    lines.append(".endc")
    lines.append(".end")
    return "\n".join(lines) + "\n"


def _tran_options(reltol):
    # card #536 defect 2: `reltol` (None = ngspice default 1e-3) is the knob that
    # sets the .tran truncation-noise floor of perturbed-minus-baseline.
    return ".options trtol=1" + (f" reltol={reltol:g}" if reltol is not None else "")


def tran_deck_text(body, stabtime_s, T, tstep, out_expr, orbit_name, reltol=None):
    tstop = stabtime_s + 2.0 * T
    lines = list(body)
    lines.append(_tran_options(reltol))
    lines.append(f".tran {tstep!r} {tstop!r} uic")
    lines.append(".control")
    lines.append("set filetype=ascii")
    lines.append("set wr_singlescale")
    lines.append("set numdgt=17")
    lines.append("run")
    lines.append(f"wrdata {orbit_name} {out_expr}")
    lines.append(".endc")
    lines.append(".end")
    return "\n".join(lines) + "\n", tstop


def tran_deck_text_window(body, tstep, tstop, out_expr, orbit_name):
    """card #493: same shape as tran_deck_text() above, but takes tstop
    DIRECTLY (a single settle+window run) instead of computing it from a
    period T -- used by the free-running ("--freerun") backend below, where
    there is no fixed period to double. Kept as a separate function (not a
    tstop=None branch bolted onto tran_deck_text) so the existing driven-mode
    callers/tests are byte-for-byte unaffected."""
    lines = list(body)
    lines.append(".options trtol=1")
    lines.append(f".tran {tstep!r} {tstop!r} uic")
    lines.append(".control")
    lines.append("set filetype=ascii")
    lines.append("set wr_singlescale")
    lines.append("set numdgt=17")
    lines.append("run")
    lines.append(f"wrdata {orbit_name} {out_expr}")
    lines.append(".endc")
    lines.append(".end")
    return "\n".join(lines) + "\n"


def read_wrdata_columns(path, n):
    """One pass: n value columns (plus a SHARED time column) from an ascii
    wrdata file written with `set wr_singlescale` (every deck this file
    builds sets it) and n vectors requested in ONE `wrdata` call. FOUND
    WHILE VALIDATING THIS TOOL (loop-gain mode's 2-output v(x)/i(Vi) read,
    the first caller here to ever request n>1 vectors in one wrdata call):
    `wr_singlescale` writes the time column ONCE per row, not once per
    vector (`t v1 v2 ... vn`, n+1 columns) -- NOT n repeated (t, v) pairs
    (2n columns), which is what tools/pssac_ref.py's own read_wrdata_pair
    assumes but never actually exercises past n=1 (every existing caller
    in this project only ever wrdata's a single vector, where the two
    layouts are indistinguishable: 2 columns either way). Returns n COPIES
    of the same shared time array (never mutated) so every existing call
    site here that indexes ts[k] per output column keeps working
    unchanged."""
    if not os.path.exists(path):
        raise CouldNotRun(f"ngspice wrote no {path}")
    ts, vs = [], [[] for _ in range(n)]
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            need = 1 + n
            if len(parts) < need:
                continue
            try:
                vals = [float(x) for x in parts[:need]]
            except ValueError:
                continue
            ts.append(vals[0])
            for k in range(n):
                vs[k].append(vals[1 + k])
    if not ts:
        raise CouldNotRun(f"{path} is empty or malformed")
    return [ts] * n, vs


# --------------------------------------------------------------------- #
# Frequency planning
# --------------------------------------------------------------------- #

EXACT_FREQ_REL = 1e-9
EXACT_QMAX_CAP = 5000


def exact_pq(ftar, fs):
    """card #536 defect 4: the (p, q) with fs*p/q == ftar exactly (to the
    decimal spelling of the two numbers), or None. Fraction(repr(x)) keeps the
    decimal the user typed (1000.0 / 500000.0 -> 1/500) instead of the binary
    approximation."""
    from fractions import Fraction
    r = Fraction(repr(float(ftar))) / Fraction(repr(float(fs)))
    if r <= 0:
        return None
    return r.numerator, r.denominator


def compute_freq_plan(freqs, fs, qmax, tol, strict_freq=False, exact=False):
    """rational_approx() (pssac_ref.py, unmodified) always returns SOME
    (p, q) -- exact by construction (f_actual = fs*p/q is always hit
    exactly, so the DFT bin is always an integer number of cycles over the
    LCM window; this is what makes the method correct, not merely
    approximate) -- but the achieved f_actual can be far from the
    REQUESTED f_target when qmax is too small to find a close rational
    approximation (pssac_ref.py's own validated table has exactly this
    case: target 1000Hz -> achieved 1250Hz on R2 at qmax=80, still trusted
    because 3.9e-6 relative error against the ANALYTIC reference at the
    ACHIEVED 1250Hz, not the requested 1000Hz). Default here (strict_freq=
    False, matching that precedent): accept the achieved frequency and
    report it honestly (relerr, f_actual) even if far from the request --
    never silently mislabel it as the requested frequency. --strict-freq
    (the RED-PROOF mode this tool's own test suite exercises): refuse
    (COULD-NOT-RUN) instead when relerr exceeds --tol, for a caller who
    needs the SPECIFIC requested frequency, not merely A commensurate one."""
    plan = []
    for ftar in freqs:
        if exact:
            # --exact-freq (card #536 defect 4): the requested frequency itself,
            # q as large as it takes (<= EXACT_QMAX_CAP), instead of the
            # smallest-q approximation within --tol (1000 Hz at fs=500 kHz is
            # 1004.02 Hz at q=498 otherwise).
            pq = exact_pq(ftar, fs)
            if pq is None or pq[1] > EXACT_QMAX_CAP:
                plan.append(dict(freq_target=ftar, p=None, q=None, f_actual=None, relerr=None, ok=False,
                                 reason=f"--exact-freq: {ftar:.10g}Hz at fs={fs:.10g}Hz needs "
                                        f"q={pq[1] if pq else '?'} > cap {EXACT_QMAX_CAP} (or is not a "
                                        f"finite ratio) -- refused, pick a frequency fs*p/q with small q"))
                continue
            p_, q_ = pq
            plan.append(dict(freq_target=ftar, p=p_, q=q_, f_actual=fs * p_ / q_,
                             relerr=abs(fs * p_ / q_ - ftar) / ftar, ok=True))
            continue
        got = ref.rational_approx(ftar, fs, qmax=qmax, tol=tol)
        if got is None:
            # rational_approx() (pssac_ref.py) returns None when not even ONE
            # q in [1, qmax] gives an integer p>=1 (e.g. ftar << fs/qmax) --
            # refuse cleanly rather than crash on the unpack.
            plan.append(dict(freq_target=ftar, p=None, q=None, f_actual=None, relerr=None, ok=False,
                             reason=f"no integer p>=1 for any q in [1,{qmax}] at fs={fs:.6g}Hz -- "
                                    f"{ftar:.6g}Hz is too far below fs/qmax; raise --qmax"))
            continue
        p, q, f_actual, relerr = got
        ok = (not strict_freq) or relerr <= tol
        entry = dict(freq_target=ftar, p=p, q=q, f_actual=f_actual, relerr=relerr, ok=ok)
        if not ok:
            entry["reason"] = (f"RED-PROOF (--strict-freq): no commensurate f=fs*p/q within "
                               f"tol={tol*100:.3f}% up to qmax={qmax} (best available: p={p} q={q} "
                               f"f_actual={f_actual:.3f}Hz, relerr={relerr*100:.3f}%) -- refused, "
                               f"not silently reported at the wrong frequency")
        elif relerr > tol:
            entry["note"] = (f"f_actual={f_actual:.3f}Hz is {relerr*100:.3f}% away from the "
                             f"requested {ftar:.3f}Hz (no closer commensurate point found within "
                             f"qmax={qmax}) -- reported honestly, not silently substituted")
        plan.append(entry)
    return plan


# --------------------------------------------------------------------- #
# Point construction + campaign dispatch, single-injection mode
# --------------------------------------------------------------------- #

def _mkpoint(deck_text, orbit_name, role, **meta):
    d = dict(deck_text=deck_text, orbit_name=orbit_name, role=role)
    d.update(meta)
    return d


def _build_point_deck(point, idx):
    return point["deck_text"]


def plan_single_pss(freq_ok, body, injection, amplitude, fs, ppc, stabtime, oscnode,
                     sc_iter, steady_coeff, out_expr, halve):
    points, base_by_q, full_by_f, half_by_f = [], {}, {}, {}
    for q in sorted({e["q"] for e in freq_ok}):
        f_lcm = fs / q
        base_body = apply_injection(body, injection, 0.0, 0.0)
        idx = len(points)
        orbit = f"pt_{idx:04d}_orbit.txt"
        text = pss_deck_text(base_body, f_lcm, stabtime, oscnode, q * ppc, out_expr, orbit,
                             sc_iter=sc_iter, steady_coeff=steady_coeff)
        points.append(_mkpoint(text, orbit, "base", q=q, T=1.0 / f_lcm, N=q * ppc))
        base_by_q[q] = idx
    for e in freq_ok:
        q, f_actual = e["q"], e["f_actual"]
        f_lcm = fs / q
        for role, frac, dest in (("full", 1.0, full_by_f), ("half", 0.5, half_by_f)):
            if role == "half" and not halve:
                continue
            pb = apply_injection(body, injection, amplitude * frac, f_actual)
            idx = len(points)
            orbit = f"pt_{idx:04d}_orbit.txt"
            text = pss_deck_text(pb, f_lcm, stabtime, oscnode, q * ppc, out_expr, orbit,
                                 sc_iter=sc_iter, steady_coeff=steady_coeff)
            points.append(_mkpoint(text, orbit, role, q=q, freq_target=e["freq_target"],
                                    f_actual=f_actual, T=1.0 / f_lcm, N=q * ppc,
                                    amplitude=amplitude * frac))
            dest[e["freq_target"]] = idx
    return points, base_by_q, full_by_f, half_by_f


def plan_single_tran(freq_ok, body, injection, amplitude, fs, ppc, stabtime_s, out_expr, halve,
                     reltol=None):
    points, base_by_q, full_by_f, half_by_f = [], {}, {}, {}
    tclk = 1.0 / fs
    tstep = tclk / ppc
    for q in sorted({e["q"] for e in freq_ok}):
        T = q / fs
        base_body = apply_injection(body, injection, 0.0, 0.0)
        idx = len(points)
        orbit = f"pt_{idx:04d}_orbit.txt"
        text, tstop = tran_deck_text(base_body, stabtime_s, T, tstep, out_expr, orbit, reltol=reltol)
        points.append(_mkpoint(text, orbit, "base", q=q, T=T, N=q * ppc, tstop=tstop, stabtime_s=stabtime_s))
        base_by_q[q] = idx
    for e in freq_ok:
        q, f_actual = e["q"], e["f_actual"]
        T = q / fs
        for role, frac, dest in (("full", 1.0, full_by_f), ("half", 0.5, half_by_f)):
            if role == "half" and not halve:
                continue
            pb = apply_injection(body, injection, amplitude * frac, f_actual)
            idx = len(points)
            orbit = f"pt_{idx:04d}_orbit.txt"
            text, tstop = tran_deck_text(pb, stabtime_s, T, tstep, out_expr, orbit, reltol=reltol)
            points.append(_mkpoint(text, orbit, role, q=q, freq_target=e["freq_target"],
                                    f_actual=f_actual, T=T, N=q * ppc, amplitude=amplitude * frac,
                                    tstop=tstop, stabtime_s=stabtime_s))
            dest[e["freq_target"]] = idx
    return points, base_by_q, full_by_f, half_by_f


# --------------------------------------------------------------------- #
# card #493: free-running (autonomous) backend.
#
# Every backend above needs a KNOWN, FIXED switching frequency `--fs`,
# either to build a `.pss method=shoot` deck whose Fguess/period IS fs
# (driven engine, pssac_ref.py/docs/pss-pac.md's own restriction) or, in
# the driven `.tran` fallback, to pick a commensurate f_actual=fs*p/q so
# `perturbed(t) - baseline(t)` is EXACTLY periodic over one LCM period and
# a plain rectangular-window single-bin DFT is exact. Neither holds for a
# free-running (hysteretic/COT/QR) converter: `--fs` is at best a rough
# ESTIMATE (the true switching frequency is an output, docs/pss-pac.md's
# free-running section), and -- the harder problem -- a perturbation
# generally SHIFTS that frequency too (this is exactly the period-variation
# effect `docs/pss-pac.md`'s pacshoot.c correction derives), so the
# perturbed and baseline orbits do not share a common period at all: there
# is no T for which perturbed(t)-baseline(t) is periodic, and picking a
# commensurate f_actual against fs (the driven backends' whole trick) does
# not fix that.
#
# card #489 SECOND pass's own validation table sidestepped this with an ad
# hoc, uncommitted script: a long (~500+ switching-period), Hann-windowed,
# single-frequency DFT on independent baseline/perturbed `.tran` runs (NOT
# a time-domain subtraction of two orbits assumed to share a period) --
# "docs/pss-pac.md ... What was not done" explicitly flagged that script as
# not promoted into a real tool. This is that promotion: `run_sweep_freerun()`
# below is the SAME technique, generalised (any deck/injection/output, same
# as every other backend here) and reusing this file's own DFT/resample/
# injection machinery (_extract_h(hann=True), apply_injection(), dft_bin()).
#
# COMMENSURATE-FREQUENCY HANDLING, the free-running analogue: instead of
# picking f_actual=fs*p/q (against the switching frequency, which is not
# fixed here), this backend picks f_actual=M/W for an INTEGER cycle count M
# against the analysis WINDOW length W (shared across baseline and every
# perturbed run at every requested frequency) -- exact by construction,
# same "always exact, report relerr honestly" principle
# compute_freq_plan()/rational_approx() use for the driven backends, just
# anchored to the window instead of to fs. W itself is sized generously
# (>= `--freerun-min-switch-periods` switching periods at the `--fs`
# ESTIMATE, and >= `--freerun-cycles` cycles of the LOWEST requested
# frequency) so the Hann taper's leakage suppression has enough switching
# periods AND enough injected-signal cycles to work with -- both floors are
# reported in each row's own note.
# --------------------------------------------------------------------- #

def _freerun_window_len(freq_targets, fs_est, cycles, min_switch_periods):
    w_cycles = cycles / min(freq_targets)
    w_floor = min_switch_periods / fs_est
    return max(w_cycles, w_floor)


def _freerun_window_n(W, fs_est, ppc):
    return max(int(round(W * fs_est * ppc)), 16)


def plan_freerun_tran(freq_targets, body, injection, amplitude, fs_est, ppc, stabtime_s, W,
                      out_expr, halve):
    tstep = 1.0 / (fs_est * ppc)
    tstop = stabtime_s + W
    points = []
    base_body = apply_injection(body, injection, 0.0, 0.0)
    base_idx = 0
    orbit = "pt_0000_orbit.txt"
    text = tran_deck_text_window(base_body, tstep, tstop, out_expr, orbit)
    points.append(_mkpoint(text, orbit, "base", tstop=tstop, stabtime_s=stabtime_s, W=W))
    full_by_f, half_by_f, freq_entries = {}, {}, []
    for ft in freq_targets:
        M = max(1, int(round(ft * W)))
        f_actual = M / W
        relerr = abs(f_actual - ft) / ft if ft else 0.0
        freq_entries.append(dict(freq_target=ft, f_actual=f_actual, M=M, relerr=relerr))
        for role, frac, dest in (("full", 1.0, full_by_f), ("half", 0.5, half_by_f)):
            if role == "half" and not halve:
                continue
            pb = apply_injection(body, injection, amplitude * frac, f_actual)
            idx = len(points)
            orbit = f"pt_{idx:04d}_orbit.txt"
            text = tran_deck_text_window(pb, tstep, tstop, out_expr, orbit)
            points.append(_mkpoint(text, orbit, role, freq_target=ft, f_actual=f_actual, M=M,
                                   amplitude=amplitude * frac, tstop=tstop, stabtime_s=stabtime_s, W=W))
            dest[ft] = idx
    return points, base_idx, full_by_f, half_by_f, freq_entries


ANALYSIS_CARDS = (".pac", ".pstb", ".pss", ".tran", ".ac", ".dc", ".noise", ".pxf", ".pnoise",
                  ".op", ".tf", ".sens", ".disto", ".pz", ".fourier")
STRIPPED_CARDS = []  # analysis cards ignored in the last load_deck_body() call (for the summary)

_IDEAL_SWITCH_MODELS = ("sw", "csw", "d_pwm")


def detect_ideal_switches(body):
    """card #536 defect 3: names of ideal-switch devices in `body` -- `S`/`W`
    devices (sw/csw) or XSPICE `A` devices whose .model type is sw, csw or
    d_pwm. Their switching INSTANTS move with the perturbation by a time
    comparable to a coarse step, so the fixed grid quantises the event shift
    (handbook example 10: 3 % error at 200 points/cycle, 50 mV)."""
    models = {}
    for l in body:
        t = l.strip().split()
        if len(t) >= 3 and t[0].lower() == ".model":
            models[t[1].lower()] = re.split(r"[(\s=]", t[2].lower())[0]
    found = []
    for l in body:
        t = l.strip().split()
        if len(t) < 5 or t[0][0].lower() not in "swa":
            continue
        if any(models.get(x.lower()) in _IDEAL_SWITCH_MODELS for x in t[3:]):
            found.append(t[0])
    return found


def _strip_dot_cards(body, names, report=False):
    """card #493: remove top-level dot-cards (and their '+' continuation
    lines) whose keyword (case-insensitive) is in `names`. Used ONLY by the
    free-running backend below -- every OTHER backend in this file is
    already validated (tests/test_pac_perturbed.py) with a fixture's own
    leftover analysis card (e.g. `.pac ... method="shoot" solver="td"`)
    left in body, and is left byte-for-byte alone here.

    WHY: found investigating this card's own R5 sweep -- a fixture's
    `.pac ... trigger="q"` card, when its underlying PSS-FREERUN shoot
    CONVERGES (which it reliably does on R5), consumes the appended
    `.control` block's `run` command ENTIRELY: ngspice's batch dispatch
    ran the .pac card's own PSS-FREERUN job (converging, emitting its own
    single-period "Time Domain Periodic Steady State Analysis" trajectory,
    then refusing the .pac post-processing itself, as designed) and the
    freshly-appended `.tran` card, later in the same deck, never ran at
    all -- no error, no warning, `run`'s exit code was still 0. wrdata then
    wrote the .pac card's OWN ~2.4us one-period orbit instead of this
    tool's intended multi-millisecond Hann-DFT window, giving a near-
    identical baseline/perturbed pair (both are that same PSS orbit,
    independent of injection amplitude, since the .pac card's own
    perturbation machinery is a completely different, unrelated code path
    from the plain `.tran` this tool actually wants) and a bogus H around
    -200dB (numerical noise, not a real near-zero gain). The driven-mode
    tran-fallback backends elsewhere in this file (plan_single_tran(), the
    existing, validated path) are not affected the same way in practice
    only because a q too large for shooting-PSS to converge is EXACTLY
    the condition that sends a point down that fallback in the first
    place -- the leftover .pss/.pac card's own attempt FAILS there and
    does not produce a competing plot -- not because the collision itself
    is absent; this is a latent risk in those paths too, flagged here,
    NOT fixed there (out of scope for this card -- changing an already-
    validated, already-tested code path without a card asking for it)."""
    out = []
    removed = []
    skip_continuation = False
    names = {nm.lower() for nm in names}
    for line in body:
        stripped = line.strip()
        if skip_continuation and stripped.startswith("+"):
            continue
        skip_continuation = False
        first = stripped.split()[0].lower() if stripped.split() else ""
        if first in names:
            skip_continuation = True
            removed.append(first)
            continue
        out.append(line)
    return (out, removed) if report else out


def run_sweep_freerun(args):
    body = load_deck_body(os.path.abspath(args.deck), params=[p.split("=", 1) for p in (args.param or [])])
    run_notes = resolve_defaults(args, body, driven=False)
    injection = parse_injection(args.inject)
    out_exprs = parse_outputs(args.outputs)
    out_expr_joined = " ".join(out_exprs)
    freq_targets = parse_freqs(args.freqs, args.sweep_spec)
    fs_est = args.fs
    halve = not args.no_halve
    n = len(out_exprs)

    stabtime_s = args.stabtime if args.stabtime is not None else max(500.0 / fs_est, 500e-6)
    W = _freerun_window_len(freq_targets, fs_est, args.freerun_cycles, args.freerun_min_switch_periods)
    N = _freerun_window_n(W, fs_est, args.points_per_cycle)

    os.makedirs(args.out, exist_ok=True)
    points, base_idx, full_by_f, half_by_f, freq_entries = plan_freerun_tran(
        freq_targets, body, injection, args.amplitude, fs_est, args.points_per_cycle,
        stabtime_s, W, out_expr_joined, halve)
    results, decks_dir = _run_batch(points, args.ngspice, os.path.join(args.out, "freerun_tran"),
                                    args.timeout, args.max_procs, probe_n=1)

    note_common = (f"freerun Hann-DFT: stabtime={stabtime_s:.6g}s W={W:.6g}s "
                  f"(~{W * fs_est:.0f} switching periods @ fs_est={fs_est:.6g}Hz) N={N}")
    rows = []
    base_r = results[base_idx]
    if not base_r.get("ok"):
        for e in freq_entries:
            entry = dict(freq_target=e["freq_target"], p=e["M"], q=1, f_actual=e["f_actual"],
                        relerr=e["relerr"])
            rows.append(_row(entry, None, "COULD-NOT-RUN", out_exprs,
                             note=f"freerun baseline .tran failed: "
                                  f"{(base_r or {}).get('stderr_tail', '')[-300:]}"))
        return rows, out_exprs

    ts_b, vs_b = read_wrdata_columns(os.path.join(decks_dir, points[base_idx]["orbit_name"]), n)
    for e in freq_entries:
        ft = e["freq_target"]
        entry = dict(freq_target=ft, p=e["M"], q=1, f_actual=e["f_actual"], relerr=e["relerr"])
        full_idx = full_by_f[ft]
        full_r = results[full_idx]
        if not full_r.get("ok"):
            rows.append(_row(entry, None, "COULD-NOT-RUN", out_exprs,
                             note=f"freerun perturbed .tran failed: {full_r.get('stderr_tail', '')[-300:]}"))
            continue
        try:
            ts_f, vs_f = read_wrdata_columns(os.path.join(decks_dir, points[full_idx]["orbit_name"]), n)
            H_full = _extract_h(ts_b, vs_b, ts_f, vs_f, W, N, e["M"], e["f_actual"], args.amplitude,
                                window0=stabtime_s, hann=True)
            H_half = None
            half_idx = half_by_f.get(ft)
            if half_idx is not None and results[half_idx].get("ok"):
                ts_h, vs_h = read_wrdata_columns(os.path.join(decks_dir, points[half_idx]["orbit_name"]), n)
                H_half = _extract_h(ts_b, vs_b, ts_h, vs_h, W, N, e["M"], e["f_actual"],
                                    args.amplitude / 2.0, window0=stabtime_s, hann=True)
            wall = ((base_r.get("wall_s") or 0) + (full_r.get("wall_s") or 0)
                    + ((results[half_idx].get("wall_s") or 0) if half_idx is not None else 0))
            r = dict(entry=entry, backend="tran-freerun-hann", H=H_full, H_half=H_half, wall=wall,
                     periodicity=None)
            rows.append(_row(entry, r, "ok", out_exprs, note=note_common))
        except CouldNotRun as exc:
            rows.append(_row(entry, None, "COULD-NOT-RUN", out_exprs, note=str(exc)))
    return rows, out_exprs


def _run_batch(points, ngspice, outdir, timeout, max_procs, probe_n=2):
    if not points:
        return [], os.path.join(outdir, "decks")
    manifest = campaign.run_campaign_a(points, _build_point_deck, ngspice, outdir,
                                       timeout=timeout, max_procs=max_procs, probe_n=min(probe_n, len(points)))
    return manifest["results"], os.path.join(outdir, "decks")


def _pss_point_ok(result):
    if not result or not result.get("ok"):
        return False, (result or {}).get("stderr_tail", "no result")[-400:]
    text = result.get("stdout_tail", "") + result.get("stderr_tail", "")
    if "FAILED to converge" in text:
        line = next((l for l in text.splitlines() if "FAILED to converge" in l), "")
        return False, f"PSS-SHOOT did not converge: {line}"
    if "Timestep too small" in text or "timestep too small" in text.lower():
        return False, "transient/monodromy step failure (timestep too small)"
    return True, None


def _orbit_ok(decks_dir, orbit_name, n):
    path = os.path.join(decks_dir, orbit_name)
    ts, vs = read_wrdata_columns(path, n)
    for col in vs:
        if not ref._orbit_plausible(col):
            return False, "implausible orbit values (q>1 shooting-PSS emission corruption, see " \
                          "research/near-fs2.md's 'Engine finding')", (ts, vs)
    return True, None, (ts, vs)


def hann_window(n):
    """card #493: a plain Hann taper, n samples, 1.0 at the centre tapering
    to 0.0 at both edges (n==1 -> [1.0], matching numpy.hanning's own
    n<=1 convention). Used only by the free-running ("--freerun") backend
    below -- see _extract_h()'s own `hann=` argument."""
    if n <= 1:
        return [1.0] * n
    return [0.5 - 0.5 * math.cos(2.0 * math.pi * k / (n - 1)) for k in range(n)]


def _extract_h(ts_base, vs_base, ts_pert, vs_pert, T, N, p, f_actual, amplitude, window0=0.0,
               hann=False):
    """`hann=True` (card #493, the free-running backend): apply a Hann taper
    to the resampled baseline/perturbed/injected-reference windows before
    the single-bin DFT, INSTEAD OF the driven-mode assumption this function
    otherwise relies on (that baseline and perturbed orbits share EXACTLY
    the same period T, so a plain rectangular window of one period is
    already periodic and leakage-free). A free-running orbit's own period
    is not fixed -- it is what card #493 is investigating -- so `diff =
    perturbed - baseline` is not periodic over an arbitrary window and a
    rectangular DFT leaks the (incommensurate) switching content into the
    target bin; the Hann taper suppresses that leakage instead (same
    technique docs/pss-pac.md's free-running section's own ad hoc reference
    script used: "a long, Hann-windowed, single-frequency DFT on a plain
    cold-start .tran"). The taper is applied identically to baseline,
    perturbed AND the known-injected reference waveform, so its own
    amplitude scaling cancels exactly in the Xout/Xin ratio below -- no
    window-gain correction needed. Default False: existing driven-mode
    callers (rectangular window, exact periodicity) are byte-for-byte
    unaffected."""
    out = []
    win = hann_window(N) if hann else None
    for base_col, pert_col in zip(vs_base, vs_pert):
        vb_u = ref.resample_uniform(ts_base[0], base_col, T, N) if window0 == 0.0 else \
            ref.resample_window(ts_base[0], base_col, window0, T, N)
        vp_u = ref.resample_uniform(ts_pert[0], pert_col, T, N) if window0 == 0.0 else \
            ref.resample_window(ts_pert[0], pert_col, window0, T, N)
        inj = [amplitude * math.sin(2.0 * math.pi * f_actual * (window0 + T * k / N)) for k in range(N)]
        if hann:
            vb_u = [a * w for a, w in zip(vb_u, win)]
            vp_u = [a * w for a, w in zip(vp_u, win)]
            inj = [a * w for a, w in zip(inj, win)]
        diff = [b - a for a, b in zip(vb_u, vp_u)]
        Xout = ref.dft_bin(diff, p)
        Xin = ref.dft_bin(inj, p)
        if abs(Xin) < 1e-30:
            raise CouldNotRun("injected-reference DFT bin is ~0 -- bug")
        out.append(Xout / Xin)
    return out


def _tran_residual(ts, vcol, stabtime_s, T, N):
    """Ripple-relative periodicity residual between two successive LCM-period
    windows (same idea as pssac_ref.py's tran-backend residual). FIX vs that
    file's own metric (found while validating THIS tool on pac_lti.cir,
    which has no real ripple on its output -- the unrelated clock network
    never reaches it): dividing by ripple alone blows up on a near-flat
    signal (residual/~0 => a spuriously huge number even though the
    absolute mismatch is numerical noise). Floor the denominator at the
    window's own signal SCALE too, not just its ripple, so a genuinely flat
    output does not get spuriously refused."""
    w1 = ref.resample_window(ts, vcol, stabtime_s, T, N)
    w2 = ref.resample_window(ts, vcol, stabtime_s + T, T, N)
    ripple = max(w2) - min(w2)
    resid = max(abs(b - a) for a, b in zip(w1, w2))
    scale = max(ripple, 1e-6 * max((abs(x) for x in w2), default=0.0), 1e-12)
    return resid / scale


# --------------------------------------------------------------------- #
# Single-injection sweep orchestration
# --------------------------------------------------------------------- #

def resolve_defaults(args, body, driven=True):
    """card #536 defect 3: --amplitude / --points-per-cycle default to None and
    are resolved here. Plain decks: 0.01 V, 200 points per cycle (unchanged).
    Ideal-switch decks (detect_ideal_switches): 0.05 V and 800 points per cycle
    -- measured on handbook example 10 (exact LPTV model 15.8556 dB at
    1004.016 Hz): 200 ppc/10 mV 16.12 dB (+0.26 dB), 800 ppc/10 mV 15.964,
    800 ppc/50 mV 15.8549, 3200 ppc/50 mV 15.8567. An explicit CLI value always
    wins. Returns human-readable notes for the summary."""
    notes = []
    sw = detect_ideal_switches(body) if driven else []
    args.ideal_switches = sw
    if args.amplitude is None:
        args.amplitude = 0.05 if sw else 0.01
        if sw:
            notes.append("ideal switches detected (" + ", ".join(sw[:4]) + (" ..." if len(sw) > 4 else "")
                         + f"): default amplitude {args.amplitude:g} V")
    if args.points_per_cycle is None:
        args.points_per_cycle = 800 if sw else 200
        if sw:
            notes.append(f"ideal switches: default grid {args.points_per_cycle} points/cycle "
                         f"(the perturbation-induced event shift must be >> the step)")
    return notes


def _wrap_deg(d):
    while d > 180.0:
        d -= 360.0
    while d < -180.0:
        d += 360.0
    return d


def _pss_stage(args, entries, body, injection, out_exprs, oscnode, stabtime, ppc, halve, outdir):
    """One `.pss method=shoot` pass over `entries` at `ppc` points per cycle.
    Returns (good {freq_target: result}, failed [(entry, reason)])."""
    out_expr_joined = " ".join(out_exprs)
    pts, base_by_q, full_by_f, half_by_f = plan_single_pss(
        entries, body, injection, args.amplitude, args.fs, ppc, stabtime, oscnode,
        args.sc_iter, args.steady_coeff, out_expr_joined, halve)
    results, decks_dir = _run_batch(pts, args.ngspice, outdir, args.timeout, args.max_procs)
    good, failed = {}, []
    n = len(out_exprs)
    for e in entries:
        ft, q = e["freq_target"], e["q"]
        try:
            base_r = results[base_by_q[q]]
            base_ok, base_reason = _pss_point_ok(base_r)
            full_ok_run, full_reason = _pss_point_ok(results[full_by_f[ft]]) if ft in full_by_f else (False, "not run")
            half_idx = half_by_f.get(ft)
            half_ok_run, half_reason = (_pss_point_ok(results[half_idx]) if half_idx is not None else (True, None))
            reason = base_reason or full_reason or half_reason
            ok = base_ok and full_ok_run and half_ok_run
            if ok:
                base_col_ok, r1, base_orbit = _orbit_ok(decks_dir, pts[base_by_q[q]]["orbit_name"], n)
                full_col_ok, r2, full_orbit = _orbit_ok(decks_dir, pts[full_by_f[ft]]["orbit_name"], n)
                ok = base_col_ok and full_col_ok
                reason = reason or r1 or r2
                if ok and half_idx is not None:
                    half_col_ok, r3, half_orbit = _orbit_ok(decks_dir, pts[half_idx]["orbit_name"], n)
                    ok = ok and half_col_ok
                    reason = reason or r3
            if ok:
                T, N, p = pts[full_by_f[ft]]["T"], pts[full_by_f[ft]]["N"], e["p"]
                H_full = _extract_h(base_orbit[0], base_orbit[1], full_orbit[0], full_orbit[1],
                                    T, N, p, e["f_actual"], args.amplitude)
                H_half = None
                if half_idx is not None:
                    H_half = _extract_h(base_orbit[0], base_orbit[1], half_orbit[0], half_orbit[1],
                                        T, N, p, e["f_actual"], args.amplitude / 2.0)
                wall = (base_r.get("wall_s") or 0) + (results[full_by_f[ft]].get("wall_s") or 0) \
                    + ((results[half_idx].get("wall_s") or 0) if half_idx is not None else 0)
                good[ft] = dict(entry=e, backend="pss", H=H_full, H_half=H_half, wall=wall,
                                periodicity=None, ppc=ppc)
            else:
                failed.append((e, reason))
        except CouldNotRun as exc:
            failed.append((e, str(exc)))
    return good, failed


def _tran_stage(args, entries, body, injection, out_exprs, tran_stabtime, ppc, halve, outdir, notes):
    """One long-`.tran` pass over `entries` at `ppc` points per cycle. Returns
    {freq_target: result}; a result without "H" carries an "error"."""
    out_expr_joined = " ".join(out_exprs)
    pts, base_by_q, full_by_f, half_by_f = plan_single_tran(
        entries, body, injection, args.amplitude, args.fs, ppc, tran_stabtime, out_expr_joined,
        halve, reltol=args.reltol)
    results, decks_dir = _run_batch(pts, args.ngspice, outdir, args.timeout, args.max_procs)
    n = len(out_exprs)
    out = {}
    for e in entries:
        ft, q = e["freq_target"], e["q"]
        note = notes.get(ft)
        sfx = f" (pss attempt: {note})" if note else ""
        try:
            base_r = results[base_by_q[q]]
            base_ok, base_reason = base_r["ok"], (None if base_r["ok"] else base_r.get("stderr_tail", "")[-300:])
            full_r = results[full_by_f[ft]]
            half_idx = half_by_f.get(ft)
            half_ok = results[half_idx]["ok"] if half_idx is not None else True
            if not (base_ok and full_r["ok"] and half_ok):
                reason = base_reason or full_r.get("stderr_tail", "")[-300:]
                out[ft] = dict(entry=e, backend="tran", H=None, H_half=None, wall=0.0, periodicity=None,
                               ppc=ppc, error=f".tran run failed: {reason}" + sfx)
                continue
            ts_b, vs_b = read_wrdata_columns(os.path.join(decks_dir, pts[base_by_q[q]]["orbit_name"]), n)
            ts_f, vs_f = read_wrdata_columns(os.path.join(decks_dir, pts[full_by_f[ft]]["orbit_name"]), n)
            T, N = pts[full_by_f[ft]]["T"], pts[full_by_f[ft]]["N"]
            resid_base = max(_tran_residual(ts_b[k], vs_b[k], tran_stabtime, T, N) for k in range(n))
            resid_full = max(_tran_residual(ts_f[k], vs_f[k], tran_stabtime, T, N) for k in range(n))
            periodicity = dict(baseline_resid_rel=resid_base, perturbed_resid_rel=resid_full)
            if resid_base > args.residual_bound or resid_full > args.residual_bound:
                out[ft] = dict(
                    entry=e, backend="tran", H=None, H_half=None, wall=0.0, periodicity=periodicity, ppc=ppc,
                    error=f"RED-PROOF: periodicity residual ({resid_base:.3g}/{resid_full:.3g}, "
                          f"ripple-relative) exceeds --residual-bound={args.residual_bound:.3g} -- "
                          f"orbit not settled, refused rather than reporting an untrustworthy H")
                continue
            H_full = _extract_h(ts_b, vs_b, ts_f, vs_f, T, N, e["p"], e["f_actual"], args.amplitude,
                                window0=tran_stabtime + T)
            H_half = None
            if half_idx is not None:
                ts_h, vs_h = read_wrdata_columns(os.path.join(decks_dir, pts[half_idx]["orbit_name"]), n)
                H_half = _extract_h(ts_b, vs_b, ts_h, vs_h, T, N, e["p"], e["f_actual"], args.amplitude / 2.0,
                                    window0=tran_stabtime + T)
            wall = (base_r.get("wall_s") or 0) + (full_r.get("wall_s") or 0) \
                + ((results[half_idx].get("wall_s") or 0) if half_idx is not None else 0)
            out[ft] = dict(entry=e, backend="tran", H=H_full, H_half=H_half, wall=wall,
                           periodicity=periodicity, ppc=ppc)
        except CouldNotRun as exc:
            out[ft] = dict(entry=e, backend="tran", H=None, H_half=None, wall=0.0, periodicity=None,
                           ppc=ppc, error=f"{exc}" + sfx)
    return out


def _step_residual(H, H_ref):
    """Per-output (|dB diff|, |deg diff|) between the result at the working grid
    and at the refined grid."""
    out = []
    for h, hr in zip(H, H_ref):
        db1, d1 = ref.db_deg(h)
        db2, d2 = ref.db_deg(hr)
        out.append((abs(db1 - db2), abs(_wrap_deg(d1 - d2))))
    return out


def run_sweep(args):
    body = load_deck_body(os.path.abspath(args.deck), params=[p.split("=", 1) for p in (args.param or [])])
    run_notes = resolve_defaults(args, body)
    injection = parse_injection(args.inject)
    out_exprs = parse_outputs(args.outputs)
    oscnode = args.osc_node or guess_oscnode(out_exprs)
    freqs = parse_freqs(args.freqs, args.sweep_spec)
    plan = compute_freq_plan(freqs, args.fs, args.qmax, args.tol, strict_freq=args.strict_freq,
                             exact=args.exact_freq)

    os.makedirs(args.out, exist_ok=True)
    rows = []
    remaining_ok = [e for e in plan if e["ok"]]
    results_by_freq = {}
    halve = not args.no_halve
    ppc = args.points_per_cycle
    stabtime = args.stabtime if args.stabtime is not None else max(100.0 / args.fs, 100e-6)
    tran_stabtime = stabtime * args.tran_stabtime_mult
    pss_reasons = {}

    if args.backend == "auto":
        # card #536: a shooting solve over q*ppc points is both slow (q=64 at 800 ppc: 140 s per
        # point) and unreliable at large q (research/near-fs2.md); go straight to .tran above the cap.
        big = [e for e in remaining_ok if e["q"] * ppc > args.pss_max_points]
        pss_reasons.update({e["freq_target"]: f"q*ppc={e['q'] * ppc} > --pss-max-points={args.pss_max_points} "
                                              f"(PSS attempt skipped)" for e in big})
        remaining_ok = [e for e in remaining_ok if e not in big]
        tran_only = big
    else:
        tran_only = []
    if args.backend in ("auto", "pss") and remaining_ok:
        good, failed = _pss_stage(args, remaining_ok, body, injection, out_exprs, oscnode, stabtime,
                                  ppc, halve, os.path.join(args.out, "pss_attempt"))
        results_by_freq.update(good)
        if args.backend == "pss":
            for e, reason in failed:
                results_by_freq[e["freq_target"]] = dict(entry=e, backend="pss", H=None, H_half=None,
                                                          wall=0.0, periodicity=None, error=reason, ppc=ppc)
            remaining_ok = []
        else:
            remaining_ok = [e for e, _ in failed]
            pss_reasons.update({e["freq_target"]: reason for e, reason in failed})
    remaining_ok = remaining_ok + tran_only
    if args.backend == "tran" or (args.backend == "auto" and remaining_ok):
        results_by_freq.update(_tran_stage(args, remaining_ok, body, injection, out_exprs, tran_stabtime,
                                           ppc, halve, os.path.join(args.out, "tran_fallback"), pss_reasons))

    # card #536 defect 3: step-refinement self-check. Every point that produced a
    # value is re-measured with the grid refined (default x2) -- baseline and
    # perturbed both, same backend, no halving run -- and the dB/deg difference
    # is reported as step_resid_*; run apply_gates() gates on it.
    if args.step_check:
        ppc2 = ppc * args.step_refine
        for backend, stage in (("pss", "pss_refined"), ("tran", "tran_refined")):
            ents = [r["entry"] for r in results_by_freq.values() if r.get("H") is not None and r["backend"] == backend]
            if not ents:
                continue
            if backend == "pss":
                good2, failed2 = _pss_stage(args, ents, body, injection, out_exprs, oscnode, stabtime,
                                            ppc2, False, os.path.join(args.out, stage))
                res2 = dict(good2)
                for e, reason in failed2:
                    res2[e["freq_target"]] = dict(H=None, error=reason)
            else:
                res2 = _tran_stage(args, ents, body, injection, out_exprs, tran_stabtime, ppc2, False,
                                   os.path.join(args.out, stage), {})
            for e in ents:
                r, r2 = results_by_freq[e["freq_target"]], res2.get(e["freq_target"], {})
                r["ppc_refined"] = ppc2
                if r2.get("H") is not None:
                    r["step_resid"] = _step_residual(r["H"], r2["H"])
                    r["wall"] = (r.get("wall") or 0) + (r2.get("wall") or 0)
                else:
                    r["step_error"] = r2.get("error", "refined run missing")

    for e in plan:
        ft = e["freq_target"]
        if not e["ok"]:
            rows.append(_row(e, None, "COULD-NOT-RUN", out_exprs, note=e["reason"]))
            continue
        r = results_by_freq.get(ft)
        if r is None or r.get("H") is None:
            rows.append(_row(e, r, "COULD-NOT-RUN", out_exprs, note=(r or {}).get("error", "no result")))
            continue
        rows.append(_row(e, r, "ok", out_exprs, note=e.get("note")))
    for row in rows:
        row["run_notes"] = run_notes
    return rows, out_exprs


FREQ_MISMATCH_REL = 1e-6


def _freq_fields(entry):
    """card #536 defect 4: the frequency the numbers are AT, stated per row."""
    ft, fa, rel = entry.get("freq_target"), entry.get("f_actual"), entry.get("relerr")
    mism = bool(fa is not None and ft and abs(fa - ft) / abs(ft) > FREQ_MISMATCH_REL)
    d = dict(f_requested=ft, f_mismatch=mism)
    if mism:
        d["freq_warning"] = (f"values are at f_actual={fa:.6g} Hz, NOT the requested {ft:.6g} Hz "
                             f"({(fa - ft) / ft * 100:+.3f} %); use --exact-freq for the requested frequency")
    return d


def _row(entry, r, status, out_exprs, note=None):
    row = dict(freq_target=entry["freq_target"], p=entry.get("p"), q=entry.get("q"),
              f_actual=entry.get("f_actual"), relerr=entry.get("relerr"), status=status,
              backend=(r or {}).get("backend"), wall_s=(r or {}).get("wall"), note=note)
    row.update(_freq_fields(entry))
    row["periodicity"] = (r or {}).get("periodicity")
    row["ppc"] = (r or {}).get("ppc")
    row["ppc_refined"] = (r or {}).get("ppc_refined")
    row["step_error"] = (r or {}).get("step_error")
    row["outputs"] = {}
    if r and r.get("H") is not None:
        step = (r.get("step_resid") or [None] * len(out_exprs))
        for expr, h, hh, sr in zip(out_exprs, r["H"], r["H_half"] or [None] * len(out_exprs), step):
            if abs(h) < 1e-30 or not math.isfinite(abs(h)):
                row["status"] = "COULD-NOT-RUN"
                row["note"] = (f"H=0 / non-finite for {expr} (perturbed == baseline: the injection never "
                               f"reached the output, or the run did not execute the intended analysis)")
                row["outputs"] = {}
                return row
            db, deg = ref.db_deg(h)
            entry_o = dict(H_re=h.real, H_im=h.imag, H_db=db, H_deg=deg)
            if hh is not None and abs(hh) > 1e-30:
                dbh, degh = ref.db_deg(hh)
                entry_o["halving_db_diff"] = db - dbh
                entry_o["halving_deg_diff"] = _wrap_deg(deg - degh)
            if sr is not None:
                entry_o["step_resid_db"], entry_o["step_resid_deg"] = sr
            row["outputs"][expr] = entry_o
    return row


def apply_gates(rows, args, loopgain=False):
    """card #536 defect 2/3: a point that fails its own self-checks is not
    reported as a value. Checks: amplitude halving (|dB| <= --halving-bound-db,
    |deg| <= --halving-bound-deg) and, when run, the step-refinement residual
    (--step-bound-db/-deg). A failing row becomes status REFUSED, `outputs` is
    emptied (the measured numbers move to `refused_outputs`, for diagnosis
    only) and `note` states the measured errors. --no-refuse keeps status ok
    but sets reliable=false and the same note."""
    for row in rows:
        if row["status"] != "ok":
            row.setdefault("reliable", None)
            continue
        reasons = []
        if loopgain:
            hd, hg = row.get("halving_db_diff"), row.get("halving_deg_diff")
            if hd is not None and (abs(hd) > args.halving_bound_db or abs(_wrap_deg(hg or 0.0)) > args.halving_bound_deg):
                reasons.append(f"amplitude-halving error T: {hd:+.3f} dB / {hg:+.2f} deg "
                               f"(bound {args.halving_bound_db:g} dB / {args.halving_bound_deg:g} deg)")
        else:
            for expr, o in row["outputs"].items():
                hd, hg = o.get("halving_db_diff"), o.get("halving_deg_diff")
                if hd is not None and (abs(hd) > args.halving_bound_db or abs(hg) > args.halving_bound_deg):
                    reasons.append(f"amplitude-halving error on {expr}: {hd:+.3f} dB / {hg:+.2f} deg "
                                   f"(bound {args.halving_bound_db:g} dB / {args.halving_bound_deg:g} deg)")
                sd, sg = o.get("step_resid_db"), o.get("step_resid_deg")
                if sd is not None and (sd > args.step_bound_db or sg > args.step_bound_deg):
                    reasons.append(f"step-refinement residual on {expr}: {sd:.3f} dB / {sg:.2f} deg between "
                                   f"{row.get('ppc')} and {row.get('ppc_refined')} points/cycle "
                                   f"(bound {args.step_bound_db:g} dB / {args.step_bound_deg:g} deg)")
            if row.get("step_error"):
                reasons.append(f"step self-check could not run: {row['step_error'][:200]}")
        if not reasons:
            row["reliable"] = True
            continue
        msg = "UNRELIABLE: " + "; ".join(reasons)
        row["note"] = msg + (f" | {row['note']}" if row.get("note") else "")
        if args.no_refuse:
            row["reliable"] = False
            continue
        row["status"] = "REFUSED"
        row["reliable"] = False
        if loopgain:
            row["refused_values"] = {k: row.get(k) for k in ("T_re", "T_im", "T_db", "T_deg",
                                                              "loopGain_db", "loopGain_deg")}
            for k in row["refused_values"]:
                row[k] = None
        else:
            row["refused_outputs"] = row["outputs"]
            row["outputs"] = {}


# --------------------------------------------------------------------- #
# Loop-gain (Tian double injection) sweep orchestration
# --------------------------------------------------------------------- #

def plan_loop_pss(freq_ok, body, a, b, amplitude, fs, ppc, stabtime, sc_iter, steady_coeff, halve):
    out_expr = f"v({LOOP_PROBE_NODE}) i(Vi_pac)"
    points, base_by_q = [], {}
    variant_idx = {}  # (freq_target, which, role) -> idx
    for q in sorted({e["q"] for e in freq_ok}):
        f_lcm = fs / q
        bb = apply_loop_injection(body, a, b, "none", 0.0, 0.0)
        idx = len(points)
        orbit = f"pt_{idx:04d}_orbit.txt"
        text = pss_deck_text(bb, f_lcm, stabtime, a, q * ppc, out_expr, orbit,
                             sc_iter=sc_iter, steady_coeff=steady_coeff)
        points.append(_mkpoint(text, orbit, "base", q=q, T=1.0 / f_lcm, N=q * ppc))
        base_by_q[q] = idx
    for e in freq_ok:
        q, f_actual = e["q"], e["f_actual"]
        f_lcm = fs / q
        for which in ("v", "i"):
            for role, frac in (("full", 1.0), ("half", 0.5)):
                if role == "half" and not halve:
                    continue
                pb = apply_loop_injection(body, a, b, which, amplitude * frac, f_actual)
                idx = len(points)
                orbit = f"pt_{idx:04d}_orbit.txt"
                text = pss_deck_text(pb, f_lcm, stabtime, a, q * ppc, out_expr, orbit,
                                     sc_iter=sc_iter, steady_coeff=steady_coeff)
                points.append(_mkpoint(text, orbit, role, q=q, freq_target=e["freq_target"],
                                       f_actual=f_actual, T=1.0 / f_lcm, N=q * ppc,
                                       amplitude=amplitude * frac))
                variant_idx[(e["freq_target"], which, role)] = idx
    return points, base_by_q, variant_idx


def plan_loop_tran(freq_ok, body, a, b, amplitude, fs, ppc, stabtime_s, halve, reltol=None):
    out_expr = f"v({LOOP_PROBE_NODE}) i(Vi_pac)"
    tclk = 1.0 / fs
    tstep = tclk / ppc
    points, base_by_q = [], {}
    variant_idx = {}
    for q in sorted({e["q"] for e in freq_ok}):
        T = q / fs
        bb = apply_loop_injection(body, a, b, "none", 0.0, 0.0)
        idx = len(points)
        orbit = f"pt_{idx:04d}_orbit.txt"
        text, tstop = tran_deck_text(bb, stabtime_s, T, tstep, out_expr, orbit, reltol=reltol)
        points.append(_mkpoint(text, orbit, "base", q=q, T=T, N=q * ppc, tstop=tstop, stabtime_s=stabtime_s))
        base_by_q[q] = idx
    for e in freq_ok:
        q, f_actual = e["q"], e["f_actual"]
        T = q / fs
        for which in ("v", "i"):
            for role, frac in (("full", 1.0), ("half", 0.5)):
                if role == "half" and not halve:
                    continue
                pb = apply_loop_injection(body, a, b, which, amplitude * frac, f_actual)
                idx = len(points)
                orbit = f"pt_{idx:04d}_orbit.txt"
                text, tstop = tran_deck_text(pb, stabtime_s, T, tstep, out_expr, orbit, reltol=reltol)
                points.append(_mkpoint(text, orbit, role, q=q, freq_target=e["freq_target"],
                                       f_actual=f_actual, T=T, N=q * ppc, amplitude=amplitude * frac,
                                       tstop=tstop, stabtime_s=stabtime_s))
                variant_idx[(e["freq_target"], which, role)] = idx
    return points, base_by_q, variant_idx


def _tian_from_orbits(decks_dir, base_orbit, v_full_orbit, i_full_orbit, T, N, p, f_actual, amplitude,
                      window0=0.0):
    ts_b, vs_b = read_wrdata_columns(os.path.join(decks_dir, base_orbit), 2)
    ts_v, vs_v = read_wrdata_columns(os.path.join(decks_dir, v_full_orbit), 2)
    ts_i, vs_i = read_wrdata_columns(os.path.join(decks_dir, i_full_orbit), 2)
    Hv = _extract_h(ts_b, vs_b, ts_v, vs_v, T, N, p, f_actual, amplitude, window0=window0)
    Hi = _extract_h(ts_b, vs_b, ts_i, vs_i, T, N, p, f_actual, amplitude, window0=window0)
    v1, i1 = Hv  # v(x), i(Vi) response to the voltage-injection run
    v2, i2 = Hi  # v(x), i(Vi) response to the current-injection run
    T_loop = stbcore.tian(v1, i1, v2, i2)
    return T_loop, (v1, i1, v2, i2)


def run_sweep_loopgain(args):
    body = load_deck_body(os.path.abspath(args.deck), params=[p.split("=", 1) for p in (args.param or [])])
    a, b = args.probe_nodes.split(":")
    run_notes = resolve_defaults(args, body)
    freqs = parse_freqs(args.freqs, args.sweep_spec)
    plan = compute_freq_plan(freqs, args.fs, args.qmax, args.tol, strict_freq=args.strict_freq,
                             exact=args.exact_freq)
    os.makedirs(args.out, exist_ok=True)
    ok_entries = [e for e in plan if e["ok"]]
    remaining_ok = list(ok_entries)
    results_by_freq = {}
    stabtime = args.stabtime if args.stabtime is not None else max(100.0 / args.fs, 100e-6)
    halve = not args.no_halve

    if args.backend in ("auto", "pss") and remaining_ok:
        pts, base_by_q, vidx = plan_loop_pss(remaining_ok, body, a, b, args.amplitude, args.fs,
                                             args.points_per_cycle, stabtime, args.sc_iter,
                                             args.steady_coeff, halve)
        results, decks_dir = _run_batch(pts, args.ngspice, os.path.join(args.out, "pss_attempt"),
                                        args.timeout, args.max_procs)
        fallback = []
        for e in remaining_ok:
            ft, q = e["freq_target"], e["q"]
            try:
                needed = [base_by_q[q], vidx[(ft, "v", "full")], vidx[(ft, "i", "full")]]
                if halve:
                    needed += [vidx[(ft, "v", "half")], vidx[(ft, "i", "half")]]
                run_ok, reason = True, None
                for idx in needed:
                    ok_i, r_i = _pss_point_ok(results[idx])
                    run_ok = run_ok and ok_i
                    reason = reason or r_i
                plaus_ok, plaus_reason = True, None
                if run_ok:
                    for idx in needed:
                        ok_i, r_i, _ = _orbit_ok(decks_dir, pts[idx]["orbit_name"], 2)
                        plaus_ok = plaus_ok and ok_i
                        plaus_reason = plaus_reason or r_i
                good = run_ok and plaus_ok
                if good:
                    T, N, p = pts[vidx[(ft, "v", "full")]]["T"], pts[vidx[(ft, "v", "full")]]["N"], e["p"]
                    T_loop, raw = _tian_from_orbits(decks_dir, pts[base_by_q[q]]["orbit_name"],
                                                    pts[vidx[(ft, "v", "full")]]["orbit_name"],
                                                    pts[vidx[(ft, "i", "full")]]["orbit_name"],
                                                    T, N, p, e["f_actual"], args.amplitude)
                    T_half = raw_half = None
                    if halve:
                        T_half, raw_half = _tian_from_orbits(decks_dir, pts[base_by_q[q]]["orbit_name"],
                                                              pts[vidx[(ft, "v", "half")]]["orbit_name"],
                                                              pts[vidx[(ft, "i", "half")]]["orbit_name"],
                                                              T, N, p, e["f_actual"], args.amplitude / 2.0)
                    wall = sum((results[i].get("wall_s") or 0) for i in needed)
                    results_by_freq[ft] = dict(entry=e, backend="pss", T=T_loop, T_half=T_half, wall=wall,
                                               periodicity=None, raw=raw, raw_half=raw_half)
                else:
                    fallback.append((e, reason or plaus_reason))
            except CouldNotRun as exc:
                fallback.append((e, str(exc)))
        if args.backend == "pss":
            for e, reason in fallback:
                results_by_freq[e["freq_target"]] = dict(entry=e, backend="pss", T=None, wall=0.0,
                                                          periodicity=None, error=reason)
        else:
            remaining_ok = [e for e, _ in fallback]
            for e, reason in fallback:
                results_by_freq.setdefault(e["freq_target"], {})["_pss_fallback_reason"] = reason

    if args.backend == "tran" or (args.backend == "auto" and remaining_ok):
        tran_stabtime = stabtime * args.tran_stabtime_mult
        pts, base_by_q, vidx = plan_loop_tran(remaining_ok, body, a, b, args.amplitude, args.fs,
                                              args.points_per_cycle, tran_stabtime, halve, reltol=args.reltol)
        results, decks_dir = _run_batch(pts, args.ngspice, os.path.join(args.out, "tran_fallback"),
                                        args.timeout, args.max_procs)
        for e in remaining_ok:
            ft, q = e["freq_target"], e["q"]
            note = results_by_freq.get(ft, {}).pop("_pss_fallback_reason", None)
            try:
                needed = [base_by_q[q], vidx[(ft, "v", "full")], vidx[(ft, "i", "full")]]
                if halve:
                    needed += [vidx[(ft, "v", "half")], vidx[(ft, "i", "half")]]
                if not all(results[i]["ok"] for i in needed):
                    results_by_freq[ft] = dict(entry=e, backend="tran", T=None, wall=0.0, periodicity=None,
                                               error=".tran run failed" + (f" (pss attempt: {note})" if note else ""))
                    continue
                T, N = pts[vidx[(ft, "v", "full")]]["T"], pts[vidx[(ft, "v", "full")]]["N"]
                ts_b, vs_b = read_wrdata_columns(os.path.join(decks_dir, pts[base_by_q[q]]["orbit_name"]), 2)
                resid_base = max(_tran_residual(ts_b[k], vs_b[k], tran_stabtime, T, N) for k in range(2))
                ts_v, vs_v = read_wrdata_columns(os.path.join(decks_dir, pts[vidx[(ft, "v", "full")]]["orbit_name"]), 2)
                resid_v = max(_tran_residual(ts_v[k], vs_v[k], tran_stabtime, T, N) for k in range(2))
                ts_i, vs_i = read_wrdata_columns(os.path.join(decks_dir, pts[vidx[(ft, "i", "full")]]["orbit_name"]), 2)
                resid_i = max(_tran_residual(ts_i[k], vs_i[k], tran_stabtime, T, N) for k in range(2))
                periodicity = dict(baseline_resid_rel=resid_base, v_run_resid_rel=resid_v, i_run_resid_rel=resid_i)
                if max(resid_base, resid_v, resid_i) > args.residual_bound:
                    results_by_freq[ft] = dict(entry=e, backend="tran", T=None, wall=0.0, periodicity=periodicity,
                                               error=f"RED-PROOF: periodicity residual exceeds "
                                                     f"--residual-bound={args.residual_bound:.3g}")
                    continue
                window0 = tran_stabtime + T
                T_loop, raw = _tian_from_orbits(decks_dir, pts[base_by_q[q]]["orbit_name"],
                                                pts[vidx[(ft, "v", "full")]]["orbit_name"],
                                                pts[vidx[(ft, "i", "full")]]["orbit_name"],
                                                T, N, e["p"], e["f_actual"], args.amplitude, window0=window0)
                T_half = None
                if halve:
                    T_half, _ = _tian_from_orbits(decks_dir, pts[base_by_q[q]]["orbit_name"],
                                                  pts[vidx[(ft, "v", "half")]]["orbit_name"],
                                                  pts[vidx[(ft, "i", "half")]]["orbit_name"],
                                                  T, N, e["p"], e["f_actual"], args.amplitude / 2.0, window0=window0)
                wall = sum((results[i].get("wall_s") or 0) for i in needed)
                results_by_freq[ft] = dict(entry=e, backend="tran", T=T_loop, T_half=T_half, wall=wall,
                                           periodicity=periodicity, raw=raw)
            except CouldNotRun as exc:
                results_by_freq[ft] = dict(entry=e, backend="tran", T=None, wall=0.0, periodicity=None,
                                           error=f"{exc}" + (f" (pss attempt: {note})" if note else ""))

    rows = []
    for e in plan:
        ft = e["freq_target"]
        if not e["ok"]:
            rows.append(dict(freq_target=ft, p=e.get("p"), q=e.get("q"), f_actual=e.get("f_actual"),
                             relerr=e.get("relerr"), status="COULD-NOT-RUN", backend=None, wall_s=None,
                             note=e["reason"], T_re=None, T_im=None, T_db=None, T_deg=None,
                             loopGain_db=None, loopGain_deg=None, halving_db_diff=None, halving_deg_diff=None,
                             periodicity=None))
            continue
        r = results_by_freq.get(ft)
        if r is None or r.get("T") is None:
            rows.append(dict(freq_target=ft, p=e.get("p"), q=e.get("q"), f_actual=e.get("f_actual"),
                             relerr=e.get("relerr"), status="COULD-NOT-RUN",
                             backend=(r or {}).get("backend"), wall_s=(r or {}).get("wall"),
                             note=(r or {}).get("error", "no result"), T_re=None, T_im=None, T_db=None,
                             T_deg=None, loopGain_db=None, loopGain_deg=None, halving_db_diff=None,
                             halving_deg_diff=None, periodicity=(r or {}).get("periodicity")))
            continue
        Tl = r["T"]
        db, deg = ref.db_deg(Tl)
        lg = -Tl
        lgdb, lgdeg = ref.db_deg(lg)
        row = dict(freq_target=ft, p=e["p"], q=e["q"], f_actual=e["f_actual"], relerr=e["relerr"],
                  status="ok", backend=r["backend"], wall_s=r["wall"], note=None,
                  T_re=Tl.real, T_im=Tl.imag, T_db=db, T_deg=deg, loopGain_db=lgdb, loopGain_deg=lgdeg,
                  periodicity=r["periodicity"])
        if r.get("T_half") is not None:
            dbh, _ = ref.db_deg(r["T_half"])
            row["halving_db_diff"] = db - dbh
            row["halving_deg_diff"] = deg - ref.db_deg(r["T_half"])[1]
        else:
            row["halving_db_diff"] = row["halving_deg_diff"] = None
        rows.append(row)
    for row, e in zip(rows, plan):
        row.update(_freq_fields(e))
        row["run_notes"] = run_notes
        row["ppc"] = args.points_per_cycle
    return rows


# --------------------------------------------------------------------- #
# Output writers
# --------------------------------------------------------------------- #

def write_single_outputs(rows, out_exprs, outdir):
    json_path = os.path.join(outdir, "pac_perturbed.json")
    csv_path = os.path.join(outdir, "pac_perturbed.csv")
    with open(json_path, "w") as fh:
        json.dump(rows, fh, indent=2, default=str)
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["freq_target", "f_actual", "p", "q", "relerr", "backend", "status", "output",
                   "H_re", "H_im", "H_db", "H_deg", "halving_db_diff", "halving_deg_diff",
                   "baseline_resid_rel", "perturbed_resid_rel", "wall_s", "note",
                   "f_requested", "f_mismatch", "reliable", "step_resid_db", "step_resid_deg", "ppc"])
        for row in rows:
            per = row.get("periodicity") or {}
            tail = [row.get("f_requested"), row.get("f_mismatch"), row.get("reliable")]
            if row["status"] != "ok" or not row["outputs"]:
                w.writerow([row["freq_target"], row.get("f_actual"), row.get("p"), row.get("q"),
                           row.get("relerr"), row.get("backend"), row["status"], "", "", "", "", "",
                           "", "", per.get("baseline_resid_rel"), per.get("perturbed_resid_rel"),
                           row.get("wall_s"), row.get("note")] + tail + ["", "", row.get("ppc")])
                continue
            for expr, o in row["outputs"].items():
                w.writerow([row["freq_target"], row["f_actual"], row["p"], row["q"], row["relerr"],
                           row["backend"], row["status"], expr, o["H_re"], o["H_im"], o["H_db"], o["H_deg"],
                           o.get("halving_db_diff"), o.get("halving_deg_diff"),
                           per.get("baseline_resid_rel"), per.get("perturbed_resid_rel"),
                           row.get("wall_s"), row.get("note")] + tail
                          + [o.get("step_resid_db"), o.get("step_resid_deg"), row.get("ppc")])
    return json_path, csv_path


def write_loop_outputs(rows, outdir):
    json_path = os.path.join(outdir, "pac_perturbed_loopgain.json")
    csv_path = os.path.join(outdir, "pac_perturbed_loopgain.csv")
    with open(json_path, "w") as fh:
        json.dump(rows, fh, indent=2, default=str)
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["freq_target", "f_actual", "p", "q", "relerr", "backend", "status",
                   "T_re", "T_im", "T_db", "T_deg", "loopGain_db", "loopGain_deg",
                   "halving_db_diff", "halving_deg_diff", "baseline_resid_rel", "v_run_resid_rel",
                   "i_run_resid_rel", "wall_s", "note", "f_requested", "f_mismatch", "reliable"])
        for row in rows:
            per = row.get("periodicity") or {}
            w.writerow([row["freq_target"], row.get("f_actual"), row.get("p"), row.get("q"),
                       row.get("relerr"), row.get("backend"), row["status"], row.get("T_re"),
                       row.get("T_im"), row.get("T_db"), row.get("T_deg"), row.get("loopGain_db"),
                       row.get("loopGain_deg"), row.get("halving_db_diff"), row.get("halving_deg_diff"),
                       per.get("baseline_resid_rel"), per.get("v_run_resid_rel"), per.get("i_run_resid_rel"),
                       row.get("wall_s"), row.get("note"), row.get("f_requested"), row.get("f_mismatch"),
                       row.get("reliable")])
    return json_path, csv_path


def print_summary(rows, loopgain=False):
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_ref = sum(1 for r in rows if r["status"] == "REFUSED")
    print(f"{n_ok}/{len(rows)} points ok" + (f", {n_ref} REFUSED (failed their own self-checks, no value reported)"
                                           if n_ref else ""))
    if STRIPPED_CARDS:
        print(f"note: ignored analysis card(s) in the deck: {', '.join(sorted(set(STRIPPED_CARDS)))} "
              f"(this tool runs its own analyses)")
    for n in PULSE_FIXES:
        print(f"note: {n} (stock ngspice can abort .tran at tight tolerances on such a source; card #541)")
    for n in (rows[0].get("run_notes") if rows else None) or []:
        print(f"note: {n}")
    for r in rows:
        ft = r["freq_target"]
        if r["status"] == "REFUSED":
            print(f"  f_target={ft:.6g}Hz: REFUSED -- {r['note']}")
            continue
        if r["status"] != "ok":
            print(f"  f_target={ft:.6g}Hz: COULD-NOT-RUN ({r['note']})")
            continue
        if r.get("f_mismatch"):
            print(f"  *** f_target={ft:.6g}Hz was MEASURED AT f_actual={r['f_actual']:.6g}Hz "
                  f"({(r['f_actual'] - ft) / ft * 100:+.3f} %) -- use --exact-freq for the requested frequency")
        head = (f"  f_target={ft:.6g}Hz -> f_actual={r['f_actual']:.6g}Hz "
                f"(p/q={r['p']}/{r['q']}, relerr={r['relerr']*100:.3f}%) backend={r['backend']} "
                f"wall={r['wall_s']:.1f}s")
        flag = "" if r.get("reliable") is not False else " [UNRELIABLE, see note]"
        if loopgain:
            print(f"{head}: T={r['T_db']:.4f}dB/{r['T_deg']:.4f}deg "
                  f"loopGain={r['loopGain_db']:.4f}dB/{r['loopGain_deg']:.4f}deg"
                  + (f" halving={r['halving_db_diff']:+.4f}dB" if r.get("halving_db_diff") is not None else "")
                  + flag)
        else:
            def one(k, v):
                t = f"{k}={v['H_db']:.4f}dB/{v['H_deg']:.4f}deg"
                ex = []
                if v.get("halving_db_diff") is not None:
                    ex.append(f"halving {v['halving_db_diff']:+.4f}dB")
                if v.get("step_resid_db") is not None:
                    ex.append(f"step-resid {v['step_resid_db']:.4f}dB/{v['step_resid_deg']:.3f}deg "
                              f"@{r.get('ppc')}->{r.get('ppc_refined')}ppc")
                return t + (f" ({'; '.join(ex)})" if ex else "")
            print(f"{head}: " + ", ".join(one(k, v) for k, v in r["outputs"].items()) + flag)
        if r.get("note") and r.get("reliable") is False:
            print(f"    {r['note']}")


# --------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("sweep", help="control-to-output/line-to-output/output-impedance sweep, "
                                     "or (--loop-gain) a Tian two-injection loop-gain sweep")
    p.add_argument("--deck", required=True)
    p.add_argument("--fs", type=parse_num, required=True,
                   help="switching/clock frequency, Hz -- EXACT for a driven deck; with --freerun "
                        "(card #493) this is only a rough ESTIMATE, used to size the settle time "
                        "and analysis window, never for commensurate-frequency matching")
    p.add_argument("--freerun", action="store_true",
                   help="card #493: autonomous (hysteretic/COT/QR) deck -- switching frequency is "
                        "an OUTPUT, not fixed, and a perturbation generally shifts it, so baseline "
                        "and perturbed orbits share no common period. Uses a long, Hann-windowed, "
                        "single-bin DFT on independent .tran runs instead of periodic subtraction "
                        "(see the 'free-running (autonomous) backend' comment above "
                        "run_sweep_freerun()); .tran only, no .pss attempt, --backend/--qmax/--tol/"
                        "--strict-freq/--tran-stabtime-mult/--residual-bound do not apply; not "
                        "combinable with --loop-gain (not attempted this card, see docs)")
    p.add_argument("--freerun-cycles", type=float, default=100.0,
                   help="--freerun: minimum cycles of the LOWEST requested frequency inside the "
                        "shared analysis window W (W=max(this/f_min, --freerun-min-switch-periods/fs))")
    p.add_argument("--freerun-min-switch-periods", type=float, default=300.0,
                   help="--freerun: minimum switching periods (at the --fs ESTIMATE) inside the "
                        "shared analysis window W -- the Hann taper's own leakage suppression needs "
                        "this to be generous (card #489 second pass's own ad hoc reference used "
                        "~500; not reduced without re-checking leakage)")
    p.add_argument("--inject", help="source:<InstName> | node:<node>[:<ref>] | series:<a>:<b>")
    p.add_argument("--outputs", help="comma-separated node/branch expressions, e.g. vout,i(L1)")
    p.add_argument("--loop-gain", action="store_true")
    p.add_argument("--probe-nodes", help="a:b for --loop-gain (a=downstream, b=upstream, same "
                                        "convention as stb/lgprobe.lib)")
    p.add_argument("--freqs", help="comma-separated target frequencies")
    p.add_argument("--sweep-spec", help="'{lin|dec|oct} N fstart fstop'")
    p.add_argument("--amplitude", type=float, default=None,
                   help="injection amplitude (default 0.01; 0.05 when ideal sw/csw/d_pwm switches are "
                        "detected in the deck, card #536)")
    p.add_argument("--qmax", type=int, default=64)
    p.add_argument("--tol", type=float, default=0.006)
    p.add_argument("--strict-freq", action="store_true", help="refuse (COULD-NOT-RUN) a point "
                   "whose achieved commensurate frequency misses --tol, instead of accepting "
                   "the closest achievable one and reporting it honestly (default)")
    p.add_argument("--points-per-cycle", type=int, default=None,
                   help="grid points per switching cycle (default 200; 800 when ideal switches are "
                        "detected, card #536)")
    p.add_argument("--exact-freq", action="store_true",
                   help="card #536: measure at EXACTLY the requested frequency (q = denominator of "
                        "f/fs, up to %d) instead of the smallest-q commensurate approximation within "
                        "--tol (1 kHz at fs=500 kHz is 1004.02 Hz otherwise); refuses if q is too large"
                        % EXACT_QMAX_CAP)
    p.add_argument("--reltol", type=float, default=1e-5,
                   help="card #536: .options reltol for the .tran backend decks (default 1e-5; ngspice's "
                        "own 1e-3 leaves a truncation-noise floor in perturbed-minus-baseline that "
                        "swamps small-gain points: r3d 100/200 kHz and r4d 20 kHz missed the "
                        "amplitude-halving check by 0.5-17 dB at 1e-3 and by <=0.02 dB at 1e-5; "
                        "costs ~2.8x wall per point). 1e-3 reproduces the old behaviour.")
    p.add_argument("--pss-max-points", type=int, default=60000,
                   help="--backend auto: skip the .pss attempt (go straight to .tran) when q*points-per-cycle "
                        "exceeds this (card #536)")
    p.add_argument("--no-step-check", dest="step_check", action="store_false",
                   help="skip the step-refinement self-check (default on: every point is re-measured "
                        "on a grid refined by --step-refine and the dB/deg difference is reported "
                        "and gated). Single-injection driven mode only.")
    p.add_argument("--step-refine", type=int, default=2, help="grid refinement factor for the step self-check")
    p.add_argument("--halving-bound-db", type=float, default=0.1,
                   help="REFUSE a point whose amplitude-halving dB difference exceeds this")
    p.add_argument("--halving-bound-deg", type=float, default=1.0)
    p.add_argument("--step-bound-db", type=float, default=0.1,
                   help="REFUSE a point whose step-refinement dB residual exceeds this")
    p.add_argument("--step-bound-deg", type=float, default=1.0)
    p.add_argument("--no-refuse", action="store_true",
                   help="report a point that fails its self-checks with its value, status ok, "
                        "reliable=false and the UNRELIABLE note (default: status REFUSED, no value)")
    p.add_argument("--sc-iter", type=int, default=200)
    p.add_argument("--steady-coeff", type=float, default=1e-7)
    p.add_argument("--stabtime", type=float, default=None, help="pss StabTime, seconds (default 100/fs)")
    p.add_argument("--tran-stabtime-mult", type=float, default=10.0)
    p.add_argument("--residual-bound", type=float, default=0.05, help="tran-backend periodicity "
                   "residual bound (ripple-relative); exceeded => COULD-NOT-RUN")
    p.add_argument("--osc-node", default=None)
    p.add_argument("--backend", choices=["auto", "pss", "tran"], default="auto")
    p.add_argument("--no-halve", action="store_true")
    p.add_argument("--param", action="append", help="NAME=VALUE, repeatable, substituted into the "
                   "deck's own .param line before injection")
    p.add_argument("--ngspice", default=os.environ.get("NGSPICE", "ngspice"))
    p.add_argument("--out", default=None)
    p.add_argument("--max-procs", type=int, default=4)
    p.add_argument("--timeout", type=int, default=1800)
    args = ap.parse_args()
    args.out = os.path.abspath(args.out) if args.out else DATA_ROOT

    try:
        if args.loop_gain and args.freerun:
            raise SystemExit("--freerun --loop-gain: not implemented (card #493 did not attempt "
                             "loop-gain on a free-running orbit -- see tools/pac_perturbed.py's "
                             "own free-running backend comment)")
        if args.loop_gain:
            if not args.probe_nodes:
                raise SystemExit("--loop-gain needs --probe-nodes a:b")
            rows = run_sweep_loopgain(args)
            apply_gates(rows, args, loopgain=True)
            print_summary(rows, loopgain=True)
            jp, cp = write_loop_outputs(rows, args.out)
        else:
            if not args.inject or not args.outputs:
                raise SystemExit("--inject and --outputs are required (or use --loop-gain)")
            rows, out_exprs = run_sweep_freerun(args) if args.freerun else run_sweep(args)
            apply_gates(rows, args)
            print_summary(rows, loopgain=False)
            jp, cp = write_single_outputs(rows, out_exprs, args.out)
        print(f"wrote {jp}\nwrote {cp}")
    except CouldNotRun as e:
        print(f"COULD-NOT-RUN: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
