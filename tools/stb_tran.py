#!/usr/bin/env python3
"""stb_tran.py -- native `.stb` loop gain / phase / gain margin AT A SERIES OF
TIME POINTS along a large-signal transient (LDO startup, soft-start, a load
step), card #497 (feature-survey #12, research/feature-survey.md).

    NGSPICE_AUG=<fork prefix>/bin/ngspice-aug python3 tools/stb_tran.py DECK \\
        --probe Vprobe --freeze Vss --times "10u,50u,100u,250u,500u,1000u" \\
        --sweep "dec 30 10 10meg" [--out DIR]

WHAT "freeze the operating point at t_k" ACTUALLY MEANS HERE, AND WHY (read
this before trusting a number -- also see docs/stb-tran.md for the full
writeup and the reproducer commands this section summarises):

Native `.stb` (like `.ac`) always finds its own linearisation point by
calling the engine's ordinary DC operating-point solver (`CKTop()` ->
`NIiter()`, full Newton) immediately before loading the small-signal
matrix -- `stban.c`/`acan.c`, both unconditionally, unless `.option noopac`
is set for a linear circuit. This Newton solve converges to the circuit's
TRUE nonlinear DC equilibrium given whatever independent-source values are
present *at that moment* -- it does NOT stop early at, or get pulled toward,
any previously-loaded state. Empirically verified three independent ways
(reproducers under this docstring's own `_verify_freeze_evidence()`, also
run standalone and captured in docs/stb-tran.md):
  1. a capacitor's own per-instance `IC=` value is discarded by a plain `op`
     (Newton fully converges past it to the true equilibrium);
  2. `tran ... uic` (which DOES respect `IC=` for its own first timestep)
     followed immediately by `op`, in the SAME session, still discards it --
     the subsequent `op` re-solves from scratch;
  3. `savestate` a genuinely non-equilibrium transient state, `loadstate` it
     into a FRESH process, then run `.stb`/`ac` immediately after: the
     diode-current discriminator in docs/stb-tran.md shows the post-stb bias
     exactly matches a plain `op`'s true equilibrium, not the loaded
     checkpoint's frozen state.
So neither `savestate`/`loadstate` nor `.ic`/per-device `IC=`/`.nodeset` can
make native `.stb` linearise about an arbitrary NON-equilibrium point --
this is a real engine property (script-layer code cannot route around it;
card #497 is script-layer only, no fork/C changes).

THE APPROACH THIS TOOL ACTUALLY USES (documented trade-off, per the card's
own "document the trade-off" instruction): freeze only the deck's genuinely
EXOGENOUS, explicitly time-varying independent sources (`--freeze NAME`,
repeatable -- e.g. a soft-start reference ramp, a supply ramp, a PWL load
step) at their t_k instantaneous value (evaluated directly from the source's
own PULSE/PWL/SIN spec, in Python -- no simulation needed for this step),
rewrite them as plain `DC <value>` sources, and run a completely ordinary
`op` + native `.stb` on that frozen deck. Every OTHER node (including every
energy-storage state -- output capacitor voltage, inductor current, the
error amplifier's own compensation cap) is then whatever the circuit's own
Newton solve converges to GIVEN those frozen exogenous inputs -- this is
the standard "quasi-static" / "frozen-time" linearisation used throughout
analog and power-electronics practice for exactly this question (loop
margin during a slow soft-start or supply ramp), and it is EXACT provided
the loop's own natural dynamics settle fast compared to how fast the frozen
source itself is moving. It is NOT exact, and should not be trusted, at a
t_k where the circuit's fast states have not yet caught up with a frozen
source that is itself moving quickly (a ripple-phase snapshot mid-switching-
cycle, for instance) -- see "Validity check" below for a diagnostic this
tool prints for exactly that case, and "Known limits" in docs/stb-tran.md.

This explains the two validation identities in this card's own brief, both
now exact BY CONSTRUCTION and confirmed against an independently hand-built
reference deck in tests/test_stb_tran.py:
  - late t_k (the transient has settled): the frozen source's value has
    ALSO settled to its own final value by then, so the frozen deck IS the
    steady-state deck -- equals a plain DC-op `.stb` of the unmodified deck.
  - early t_k: equals a DC-op `.stb` of "the same circuit with its [frozen]
    sources ... set to the frozen values" -- exactly what is built and run.

Validity check (optional, `--check-tran "TSTEP TSTOP"` + `--check-node
NODE`, repeatable): runs the deck's OWN real transient once (unmodified,
same binary), and at every t_k reports max relative deviation between each
`--check-node`'s REAL transient voltage and its value in the frozen op's own
solution -- a direct, honest measure of how well the quasi-static assumption
held at that t_k (small = trustworthy freeze; large = the loop had not
caught up with the frozen source yet, read that t_k's margin with caution).

Output: a table PM(t_k)/GM(t_k)/fc(t_k) (fc = the 0 dB crossing frequency
the reported PM was measured at, `NaN`/None where there is no crossing --
same convention as tools/stbcore.margins(), never a fake number), printed to
stdout and written as times.csv + summary.json under
~/Work/data/ngspice_augmentation/stb_tran/<deck>-<timestamp>/ (never the
repository) unless --out is given. Optional plot data: --plot-loopgain also
writes each t_k's full loopgain.csv (freq, T re/im, dB, phase) for anyone
who wants to plot the Bode trace at a specific instant.

Exit status: 0 = ran (every t_k's op+stb completed -- "no crossing" is a
valid, non-failing result, not an error); 2 = could not run (ngspice-aug
missing, a t_k's run errored, or a --freeze source's waveform could not be
parsed/evaluated).
"""
import argparse
import csv
import datetime
import json
import math
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stb  # noqa: E402 -- reuses prepare_deck() (include-path fixup, drop
             # existing analyses/control blocks) and CouldNotRun

NGSPICE_AUG = os.environ.get(
    "NGSPICE_AUG") or os.environ.get(
    "NGSPICE",
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s46/bin/ngspice"))
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation")

CouldNotRun = stb.CouldNotRun

_SUFFIX = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "mil": 25.4e-6,
           "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}


def spice_num(tok):
    """"500u" -> 5e-4, "1meg" -> 1e6, plain floats/ints too. Case-insensitive,
    "meg" checked before the bare "m" (milli) suffix, same convention this
    project's other tools (tests/test_stb.py's spice_num, tools/loss.py's)
    already use."""
    s = tok.strip()
    m = re.match(r"^([-+0-9.eE]+)(meg|mil|[tgkmunpf])?[a-zA-Z]*$", s, re.I)
    if not m:
        raise CouldNotRun(f"cannot parse SPICE number {tok!r}")
    mult = _SUFFIX.get(m.group(2).lower(), 1.0) if m.group(2) else 1.0
    return float(m.group(1)) * mult


# ---------------------------------------------------------------------------
# waveform evaluation -- enough of PULSE/PWL/SIN to freeze a soft-start ramp,
# a supply ramp, or a PWL load-step source at an arbitrary time t. Anything
# else (EXP, SFFM, AM, a PWL with r=/td= repeat options) is refused, not
# silently mishandled.
# ---------------------------------------------------------------------------

_PARAM_LINE = re.compile(r"(\w+)\s*=\s*(\{[^{}]*\}|\S+)")
_NUM_TOKEN = re.compile(r"(?<![\w.])[0-9]+\.?[0-9]*(?:[eE][-+]?[0-9]+)?(?:meg|mil|[tgkmunpf])?(?![\w.])",
                         re.I)


def resolve_params(lines):
    """{param name (lowercase): float value}, resolved from every top-level
    `.param name=expr[,...]` line in deck order (a param's own expr may
    reference any EARLIER param, brace-wrapped or bare -- e.g.
    `.param tclk={1/fs}` after `.param fs=500k`, same numparam convention
    every fixture in this project already relies on). Needed because this
    tool parses the raw deck TEXT (it does not run ngspice's own numparam
    preprocessor) to evaluate a --freeze source's PULSE/PWL/SIN args, which
    are routinely written as `{param}` brace expressions, not literal
    numbers (see tests/fixtures/stbtran_ldo_softstart.cir's Vss line)."""
    params = {}
    for line in lines:
        s = line.strip()
        if not s or s[0] in "*$":
            continue
        if s.lower().startswith(".param"):
            for m in _PARAM_LINE.finditer(s[len(".param"):]):
                params[m.group(1).lower()] = _eval_expr(m.group(2), params)
    return params


def _eval_expr(expr, params):
    """A `{...}` (braces optional) numparam-style expression -> float, given
    already-resolved params. Substitutes known param names (longest first,
    so e.g. "vref_final" isn't partially shadowed by a shorter "vref"), then
    every remaining SPICE-suffixed numeric literal (500k, 1meg, 10u, ...),
    then evaluates with plain Python arithmetic (+-*/**, parens) and no
    builtins -- the same small, safe-eval convention this project's other
    scripts (e.g. tools/loss.py's Bcmp/Bhs-style expression handling) use."""
    e = expr.strip()
    if e.startswith("{") and e.endswith("}"):
        e = e[1:-1]
    for name in sorted(params, key=len, reverse=True):
        e = re.sub(rf"(?<![\w.]){re.escape(name)}(?![\w])", repr(params[name]), e, flags=re.I)
    e = _NUM_TOKEN.sub(lambda m: repr(spice_num(m.group(0))), e)
    try:
        return float(eval(e, {"__builtins__": {}}, {}))
    except Exception as ex:
        raise CouldNotRun(f"cannot evaluate expression {expr!r} (-> {e!r}): {ex}") from ex


def eval_token(tok, params):
    """One PULSE/PWL/SIN argument -> float: a `{expr}` brace expression, a
    bare param name, or a plain SPICE-suffixed number."""
    tok = tok.strip()
    if tok.startswith("{") and tok.endswith("}"):
        return _eval_expr(tok, params)
    try:
        return spice_num(tok)
    except CouldNotRun:
        return _eval_expr(tok, params)


def _args(inner, params):
    return [eval_token(t, params) for t in re.split(r"[,\s]+", inner.strip()) if t]


def eval_pulse(a, t):
    if len(a) < 2:
        raise CouldNotRun(f"PULSE needs at least V1,V2: {a}")
    v1, v2 = a[0], a[1]
    td = a[2] if len(a) > 2 else 0.0
    tr = a[3] if len(a) > 3 else 1e-9
    tf = a[4] if len(a) > 4 else 1e-9
    pw = a[5] if len(a) > 5 else 1e30
    per = a[6] if len(a) > 6 and a[6] > 0 else float("inf")
    if t < td:
        return v1
    tt = (t - td) % per if math.isfinite(per) else (t - td)
    if tt < tr:
        return v1 + (v2 - v1) * (tt / tr if tr > 0 else 1.0)
    tt -= tr
    if tt < pw:
        return v2
    tt -= pw
    if tt < tf:
        return v2 + (v1 - v2) * (tt / tf if tf > 0 else 1.0)
    return v1


def eval_pwl(a, t):
    if len(a) < 4 or len(a) % 2:
        raise CouldNotRun(f"PWL needs an even list of >=2 (time,value) pairs: {a}")
    pts = list(zip(a[0::2], a[1::2]))
    for k in range(1, len(pts)):
        if pts[k][0] < pts[k - 1][0]:
            raise CouldNotRun(f"PWL times must be non-decreasing: {a}")
    if t <= pts[0][0]:
        return pts[0][1]
    if t >= pts[-1][0]:
        return pts[-1][1]
    for k in range(len(pts) - 1):
        t0, v0 = pts[k]
        t1, v1 = pts[k + 1]
        if t0 <= t <= t1:
            return v0 if t1 == t0 else v0 + (v1 - v0) * (t - t0) / (t1 - t0)
    raise CouldNotRun(f"PWL evaluation fell through: {a} at t={t}")  # unreachable


def eval_sin(a, t):
    if len(a) < 3:
        raise CouldNotRun(f"SIN needs at least VO,VA,FREQ: {a}")
    vo, va, freq = a[0], a[1], a[2]
    td = a[3] if len(a) > 3 else 0.0
    theta = a[4] if len(a) > 4 else 0.0
    phase = a[5] if len(a) > 5 else 0.0
    if t < td:
        return vo + va * math.sin(math.radians(phase))
    damp = math.exp(-theta * (t - td)) if theta else 1.0
    return vo + va * damp * math.sin(2 * math.pi * freq * (t - td) + math.radians(phase))


WAVEFORMS = {"pulse": eval_pulse, "pwl": eval_pwl, "sin": eval_sin}
_FUNC_RE = re.compile(r"\b(pulse|pwl|sin|exp|sffm|am)\s*\(([^)]*)\)", re.I)


def freeze_value(line, t, params):
    """Instantaneous value of a "Vxxx n+ n- ... TRAN_FUNC(...)" source line at
    time t. A line with no recognised TRAN_FUNC is already a plain DC source
    -- returns its DC value unchanged (freezing it is then a documented no-op,
    not an error, since --freeze naming a source that is already constant is
    harmless)."""
    m = _FUNC_RE.search(line)
    if not m:
        dm = re.search(r"\bDC\s+(\{[^{}]*\}|[-+0-9.eE]+[a-zA-Z]*)", line, re.I)
        if dm:
            return eval_token(dm.group(1), params)
        toks = line.split()
        if len(toks) >= 4:
            try:
                return eval_token(toks[3], params)
            except CouldNotRun:
                pass
        raise CouldNotRun(f"no DC value or recognised TRAN_FUNC on: {line!r}")
    kind = m.group(1).lower()
    if kind not in WAVEFORMS:
        raise CouldNotRun(f"{kind.upper()}(...) sources are not supported by "
                           f"--freeze (only PULSE/PWL/SIN): {line!r}")
    return WAVEFORMS[kind](_args(m.group(2), params), t)


def freeze_deck(lines, freeze_names, t):
    """Deck lines with every --freeze source rewritten to a plain DC source at
    its value(t). Returns (new_lines, {name: value})."""
    params = resolve_params(lines)
    wanted = {n.lower() for n in freeze_names}
    seen = set()
    out, values = [], {}
    for line in lines:
        f = line.split()
        if f and f[0].lower() in wanted and (f[0][0].lower() in ("v", "i")):
            name = f[0]
            val = freeze_value(line, t, params)
            values[name] = val
            out.append(f"{f[0]} {f[1]} {f[2]} DC {val:.15g}")
            seen.add(f[0].lower())
        else:
            out.append(line)
    missing = wanted - seen
    if missing:
        raise CouldNotRun(f"--freeze source(s) not found in deck: {sorted(missing)}")
    return out, values


# ---------------------------------------------------------------------------
# one op + native .stb run, at a single frozen time point
# ---------------------------------------------------------------------------

def run_ngspice(wrapper, outdir, label, timeout=600):
    try:
        proc = subprocess.run([NGSPICE_AUG, "-b", wrapper], capture_output=True,
                               text=True, timeout=timeout)
    except (FileNotFoundError, PermissionError) as e:
        raise CouldNotRun(f"cannot start {NGSPICE_AUG}: {e}") from e
    except subprocess.TimeoutExpired as e:
        raise CouldNotRun(f"{label} timed out after {timeout}s") from e
    with open(os.path.join(outdir, f"{label}.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise CouldNotRun(f"ngspice-aug exit {proc.returncode}; {errs[:3]}")
    return proc


def _read_margins(mfile):
    """{pm_deg, fc_hz, gm_db, gm_freq_hz} from a wrdata'd
    stb2.phasemargin/phasemarginfreq/gainmargin/gainmarginfreq file --
    shared by both stb_at_point() (quasistatic) and stb_at_point_instant()
    (instant), same wrdata column layout either way."""
    if not os.path.exists(mfile):
        raise CouldNotRun(f"ngspice wrote no {mfile}")
    with open(mfile) as fh:
        next(fh, None)  # wr_vecnames header line
        row = [float(x) for x in fh.readline().split()]
    if len(row) == 5:
        _, pm, pmf, gm, gmf = row
    elif len(row) == 4:
        pm, pmf, gm, gmf = row
    else:
        raise CouldNotRun(f"unexpected margins column count {len(row)} in {mfile}")
    return {
        "pm_deg": None if math.isnan(pm) else pm,
        "fc_hz": None if math.isnan(pmf) else pmf,
        "gm_db": None if math.isnan(gm) else gm,
        "gm_freq_hz": None if math.isnan(gmf) else gmf,
    }


def _read_loopgain(lgfile):
    freqs, t = [], []
    if not os.path.exists(lgfile):
        return None
    with open(lgfile) as fh:
        next(fh, None)
        for row in fh:
            v = row.split()
            if len(v) == 5:
                freqs.append(float(v[0]))
                t.append(complex(float(v[3]), float(v[4])))
    return (freqs, t)


def stb_at_point(frozen_lines, probe, sweep, outdir, label, want_loopgain=False):
    """op + native .stb on an already-frozen deck; returns {pm, gm, fc,
    low_f_positive_feedback, note, loopgain: (freqs,T) or None}."""
    os.makedirs(outdir, exist_ok=True)
    mfile = os.path.join(outdir, f"{label}_margins.txt")
    lgfile = os.path.join(outdir, f"{label}_loopgain.txt")
    for f in (mfile, lgfile):
        if os.path.exists(f):
            os.remove(f)
    body = "\n".join(frozen_lines)
    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
op
stb {sweep} probe={probe}
wrdata {mfile} stb2.phasemargin stb2.phasemarginfreq stb2.gainmargin stb2.gainmarginfreq
"""
    if want_loopgain:
        ctl += f"wrdata {lgfile} stb1.loopgain stb1.t\n"
    ctl += "quit\n.endc\n.end\n"
    wrapper = os.path.join(outdir, f"{label}.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + ctl)
    run_ngspice(wrapper, outdir, label)
    result = _read_margins(mfile)
    if want_loopgain:
        lg = _read_loopgain(lgfile)
        if lg is not None:
            result["loopgain"] = lg
    return result


# ---------------------------------------------------------------------------
# --mode=instant: op-point-free linearisation at the REAL transient state at
# t_k, card #500. See this module's docstring addendum below and
# docs/stb-tran.md's "instant mode" section for the full design and why this
# needs BOTH `.tran` and `.stb ... instant` as TOP-LEVEL NETLIST CARDS run
# together via a single `run` (not as two sequential `.control` commands --
# ngspice's own per-command dispatch (if_run()/doAnalyses(ckt,RESTART,...))
# resets ckt->CKTtime to 0 and re-runs CKTsetup/CKTunsetup for EVERY
# interactive analysis command, discarding a preceding command's ending
# state regardless of whether the new command's own analysis calls CKTop;
# only jobs processed within the SAME doAnalyses(ckt, RESTART, task) call --
# i.e. netlist-level jobs run together via one `run`/batch invocation --
# share state, and only because CKTdoJob's job loop iterates by ANALYSIS
# TYPE INDEX (AC, ..., TRAN, ..., STB last), which happens to run `.tran`
# before `.stb` regardless of netlist declaration order -- found and
# verified empirically while validating this card, see docs/stb-tran.md).
# ---------------------------------------------------------------------------

def stb_at_point_instant(base_lines, probe, sweep, tstep, t, outdir, label, want_loopgain=False):
    """`.tran {tstep} {t} 0 uic` + `.stb {sweep} probe={probe} instant`, as
    top-level netlist cards run together via a single `run` -- linearises
    about the REAL (possibly non-equilibrium) transient state `.tran` left
    at t=t, no frozen-source/re-solved-equilibrium approximation at all
    (contrast stb_at_point()'s quasistatic mode, which freezes only the
    --freeze source and re-solves every OTHER state to a fresh DC
    equilibrium). base_lines must NOT already contain a `.tran`/`.stb`/
    `.control` block -- pass stb.prepare_deck()'s output, same as
    quasistatic mode."""
    os.makedirs(outdir, exist_ok=True)
    mfile = os.path.join(outdir, f"{label}_margins.txt")
    lgfile = os.path.join(outdir, f"{label}_loopgain.txt")
    for f in (mfile, lgfile):
        if os.path.exists(f):
            os.remove(f)
    body = "\n".join(base_lines)
    jobs = f'.tran {tstep} {t:.15g} 0 uic\n.stb {sweep} probe="{probe}" instant\n'
    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
run
wrdata {mfile} stb2.phasemargin stb2.phasemarginfreq stb2.gainmargin stb2.gainmarginfreq
"""
    if want_loopgain:
        ctl += f"wrdata {lgfile} stb1.loopgain stb1.t\n"
    ctl += "quit\n.endc\n.end\n"
    wrapper = os.path.join(outdir, f"{label}.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + jobs + ctl)
    run_ngspice(wrapper, outdir, label, timeout=1800)
    result = _read_margins(mfile)
    if want_loopgain:
        lg = _read_loopgain(lgfile)
        if lg is not None:
            result["loopgain"] = lg
    return result


# ---------------------------------------------------------------------------
# validity check: does the frozen op's --check-node value actually resemble
# the deck's own real .tran at t_k? (see module docstring, "Validity check")
# ---------------------------------------------------------------------------

def run_real_tran(lines, check_tran, check_nodes, outdir):
    tstep, tstop = check_tran
    trace = os.path.join(outdir, "check_tran.txt")
    if os.path.exists(trace):
        os.remove(trace)
    body = "\n".join(lines)
    sig = " ".join(f"v({n})" for n in check_nodes)
    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
tran {tstep} {tstop} 0 uic
wrdata {trace} {sig}
quit
.endc
.end
"""
    wrapper = os.path.join(outdir, "check_tran.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + ctl)
    run_ngspice(wrapper, outdir, "check_tran", timeout=1800)
    if not os.path.exists(trace):
        raise CouldNotRun(f"ngspice wrote no {trace}")
    rows = []
    with open(trace) as fh:
        next(fh, None)
        for line in fh:
            v = line.split()
            if len(v) == 1 + len(check_nodes):
                rows.append([float(x) for x in v])
    if not rows:
        raise CouldNotRun(f"{trace} is empty")
    return rows  # [t, node1, node2, ...]


def nearest_row(rows, t):
    return min(rows, key=lambda r: abs(r[0] - t))


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def sweep_times(deck, probe, freeze, times, sweep, outdir, check_tran=None,
                 check_nodes=(), plot_loopgain=False, mode="quasistatic", tstep=None):
    os.makedirs(outdir, exist_ok=True)
    base_lines = stb.prepare_deck(deck)

    if mode == "instant":
        out = []
        for t in times:
            label = f"t_{t:.9g}".replace("+", "").replace(".", "p")
            res = stb_at_point_instant(base_lines, probe, sweep, tstep, t, outdir, label,
                                        want_loopgain=plot_loopgain)
            out.append({"t": t, "frozen": {}, **res})
        return out

    tran_rows = None
    if check_tran:
        tran_rows = run_real_tran(base_lines, check_tran, check_nodes, outdir)

    out = []
    for t in times:
        label = f"t_{t:.9g}".replace("+", "").replace(".", "p")
        frozen_lines, froze = freeze_deck(base_lines, freeze, t)
        res = stb_at_point(frozen_lines, probe, sweep, outdir, label, want_loopgain=plot_loopgain)
        row = {"t": t, "frozen": froze, **res}
        if tran_rows is not None:
            nr = nearest_row(tran_rows, t)
            deltas = {}
            if abs(nr[0] - t) > max(1e-12, 0.05 * (tran_rows[1][0] - tran_rows[0][0])) and len(tran_rows) > 1:
                pass  # nearest sample not exactly at t; still report it, it's the closest we have
            frozen_op = frozen_op_values(frozen_lines, check_nodes, outdir, label)
            for k, node in enumerate(check_nodes):
                real_v = nr[1 + k]
                frozen_v = frozen_op.get(node)
                if frozen_v is not None:
                    denom = abs(real_v) if real_v else 1.0
                    deltas[node] = abs(real_v - frozen_v) / denom
            row["quasistatic_check"] = {"real_t": nr[0], "rel_delta": deltas}
        out.append(row)
    return out


def frozen_op_values(frozen_lines, nodes, outdir, label):
    """v(node) from the frozen deck's own `op` (reused for the validity check
    so it is exactly the same operating point .stb itself linearised about --
    a second, tiny ngspice call, not parsed back out of the .stb run, since
    .stb's own wrdata output above only carries the margins/loopgain, not the
    op node voltages)."""
    if not nodes:
        return {}
    opfile = os.path.join(outdir, f"{label}_op.txt")
    if os.path.exists(opfile):
        os.remove(opfile)
    body = "\n".join(frozen_lines)
    sig = " ".join(f"v({n})" for n in nodes)
    ctl = f"\n.control\nop\nprint {sig}\n.endc\n.end\n"
    wrapper = os.path.join(outdir, f"{label}_op.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + ctl)
    proc = run_ngspice(wrapper, outdir, f"{label}_op")
    vals = {}
    for m in re.finditer(r"v\((\w+)\)\s*=\s*([-+0-9.eE]+)", proc.stdout, re.I):
        vals[m.group(1).lower()] = float(m.group(2))
    return {n: vals.get(n.lower()) for n in nodes}


def write_results(rows, outdir):
    csv_path = os.path.join(outdir, "times.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t", "pm_deg", "fc_hz", "gm_db", "gm_freq_hz"])
        for r in rows:
            w.writerow([r["t"], r["pm_deg"], r["fc_hz"], r["gm_db"], r["gm_freq_hz"]])
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(rows, fh, indent=2, default=lambda o: None)


def describe(rows):
    lines = [f"{'t':>12s}  {'PM (deg)':>10s}  {'fc (Hz)':>12s}  {'GM (dB)':>10s}"]
    for r in rows:
        pm = f"{r['pm_deg']:.3f}" if r["pm_deg"] is not None else "no crossing"
        fc = f"{r['fc_hz']:.6g}" if r["fc_hz"] is not None else "--"
        gm = f"{r['gm_db']:.3f}" if r["gm_db"] is not None else "no crossing"
        lines.append(f"{r['t']:12.6g}  {pm:>10s}  {fc:>12s}  {gm:>10s}")
        qc = r.get("quasistatic_check")
        if qc and qc["rel_delta"]:
            worst = max(qc["rel_delta"].values())
            lines.append(f"{'':12s}  quasi-static check at real t={qc['real_t']:.6g}: "
                          f"worst node rel delta {worst:.3e}")
    return "\n".join(lines)


def main():
    global NGSPICE_AUG
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("--probe", required=True, help="native .stb probe V-source instance")
    ap.add_argument("--mode", choices=("quasistatic", "instant"), default="quasistatic",
                     help="quasistatic (default): freeze --freeze source(s), re-solve every other "
                          "state to a fresh DC equilibrium (docs/stb-tran.md). instant (card #500): "
                          "no freezing/re-solving at all -- run the deck's own REAL `.tran` up to "
                          "each t_k and linearise about whatever non-equilibrium state it left "
                          "(needs --tstep, not --freeze; see docs/stb.md's 'instant mode' section)")
    ap.add_argument("--freeze", action="append", default=[],
                     help="quasistatic mode only: independent source instance to freeze at each "
                          "t_k (evaluated from its own PULSE/PWL/SIN spec); repeatable; required "
                          "(and only meaningful) for --mode=quasistatic")
    ap.add_argument("--tstep", help="instant mode only: the .tran step size used to run the deck's "
                     "own real transient up to each t_k (e.g. \"50n\"); required for --mode=instant")
    ap.add_argument("--times", required=True, help='comma list of t_k, e.g. "10u,50u,500u"')
    ap.add_argument("--sweep", default="dec 30 10 10meg", help='native .stb sweep, e.g. "dec 30 10 10meg"')
    ap.add_argument("--out", help="output directory (default: under ~/Work/data/ngspice_augmentation/stb_tran)")
    ap.add_argument("--ngspice-aug", default=NGSPICE_AUG, help="fork ngspice-aug/ngspice binary (native .stb)")
    ap.add_argument("--check-tran", help='"TSTEP TSTOP" -- quasistatic mode only -- also run the '
                     "deck's own real transient once and report each t_k's quasi-static validity "
                     "(needs --check-node)")
    ap.add_argument("--check-node", action="append", default=[],
                     help="node to compare real-transient vs frozen-op voltage at each t_k; repeatable")
    ap.add_argument("--plot-loopgain", action="store_true",
                     help="also write each t_k's full loopgain.csv (freq,T) for plotting")
    args = ap.parse_args()
    NGSPICE_AUG = args.ngspice_aug

    if args.mode == "instant":
        if not args.tstep:
            print("COULD NOT RUN: --mode=instant needs --tstep", file=sys.stderr)
            return 2
        if args.freeze:
            print("COULD NOT RUN: --freeze is quasistatic-mode only (instant mode runs the deck's "
                  "own real .tran unmodified, nothing to freeze) -- do not combine with --mode=instant",
                  file=sys.stderr)
            return 2
        if args.check_tran:
            print("COULD NOT RUN: --check-tran is quasistatic-mode only (it measures how well the "
                  "frozen approximation held vs the real transient -- instant mode already IS the "
                  "real transient, nothing to check it against) -- do not combine with --mode=instant",
                  file=sys.stderr)
            return 2
    else:
        if not args.freeze:
            print("COULD NOT RUN: --mode=quasistatic (the default) needs at least one --freeze",
                  file=sys.stderr)
            return 2

    try:
        times = sorted(spice_num(t) for t in args.times.split(","))
    except CouldNotRun as e:
        print(f"COULD NOT RUN: {e}", file=sys.stderr)
        return 2

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = args.out or os.path.join(DATA_ROOT, "stb_tran",
                                       f"{os.path.splitext(os.path.basename(args.deck))[0]}-{stamp}")

    check_tran = None
    if args.check_tran:
        toks = args.check_tran.split()
        if len(toks) != 2:
            print('COULD NOT RUN: --check-tran needs "TSTEP TSTOP"', file=sys.stderr)
            return 2
        check_tran = (toks[0], toks[1])

    try:
        rows = sweep_times(args.deck, args.probe, args.freeze, times, args.sweep, outdir,
                            check_tran=check_tran, check_nodes=args.check_node,
                            plot_loopgain=args.plot_loopgain, mode=args.mode, tstep=args.tstep)
    except (CouldNotRun, subprocess.TimeoutExpired) as e:
        print(f"COULD NOT RUN: {e}", file=sys.stderr)
        return 2

    if args.plot_loopgain:
        for r in rows:
            lg = r.pop("loopgain", None)
            if lg:
                freqs, t = lg
                label = f"t_{r['t']:.9g}".replace("+", "").replace(".", "p")
                with open(os.path.join(outdir, f"{label}_loopgain.csv"), "w") as fh:
                    fh.write("freq,T_re,T_im,T_db,T_phase_deg\n")
                    for f, z in zip(freqs, t):
                        db = 20 * math.log10(abs(z)) if z else float("-inf")
                        ph = math.degrees(math.atan2(z.imag, z.real))
                        fh.write(f"{f:.9e},{z.real:.9e},{z.imag:.9e},{db:.6f},{ph:.6f}\n")

    write_results(rows, outdir)
    print(describe(rows))
    print(f"results: {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
