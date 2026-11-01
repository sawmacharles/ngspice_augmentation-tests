#!/usr/bin/env python3
"""loss.py -- per-element efficiency / loss breakdown for a switching converter
deck (card #470, plan T10.2).

Given a deck and a time window (the last N periods of a settled `.tran` on any
ngspice, or one period from a `.pss ... method=shoot` orbit on the fork build),
computes: input power (from a named source), output power (a named load), the
average dissipated power of every OTHER element grouped into categories
(switch conduction, switch transition/overlap, gate drive, diode conduction
(incl. reverse-recovery), inductor DCR, capacitor ESR, resistive/other), the
energy-balance residual, and efficiency. Device coverage: R, L, C, V, I
(constant-0), B/G/E (simple linear/expression forms), D (native diode, any
model level) and M (native MOSFET, any model level, e.g. level=1 or BSIM4
level=54) -- see "Device coverage" below. Anything else is never silently
dropped: it is listed under "excluded" and blocks a "pass" (exit 2).

    tools/loss.py DECK --vin Vin --load Rload --period 2e-6 --periods 5 \\
        --tstop 1.5e-3 [--dcr Rdcr] [--esr Resr] [--switch Bhs] [--switch Bls] \\
        [--diode Dbody]

    tools/loss.py DECK --pss --vin Vin --load Rload --fguess 500e3 \\
        --stabtime 300e-6 --oscnode vout [...]

WHY `@device[i]` / `@device[p]` LOOKED FROZEN (card #470), AND THE FIX (card
#481): on this ngspice-45.2 build, `@name[i]`/`@name[p]` read back through
`wrdata` inside a `tran ... uic` run appeared to NOT update per accepted
timestep -- frozen at their t=0 value, for every device type tried. Card #481
found the actual cause: those vectors were never `.save`d / `save`d. By
default (no `.save`/`save` at all) ngspice keeps every ordinary circuit
vector -- node voltages and the native MNA branch currents of V/L devices --
but an `@device[param]` reference is NOT an ordinary circuit vector; it is
only accumulated across the run if explicitly named in a `save` command
*issued before the analysis starts*. Once that `save @name[param]` line is
added, `@name[id]` (diodes), `@name[id]`/`[is]`/`[ig]`/`[ib]`/`[p]` (MOS1) and
`@name[p]` (inductors, etc.) all track correctly per timestep -- confirmed
against ngspice's own terminal-voltage-times-current product at several
sample points (docs/loss.md has the reproducers and the numbers). SECOND trap
found the same day: `save` is NOT additive to the "no save at all -> save
everything ordinary" default -- issuing any `save` command REPLACES that
default with exactly the named list, so a `save @dev[id]` line alone silently
stopped even `i(l1)`/`v(node)` (normally free) from being recorded, breaking
`wrdata` for everything else the moment one `@device[...]` vector was needed.
build_wrapper() therefore saves its FULL signal request (every `i_sig` and
`v(node)` too, not just the `@device[...]` ones) whenever any `@device[...]`
vector is needed at all. This tool
still prefers node voltages (`v(node)`) and native MNA branch currents
(`i(Vxxx)`/`i(Lxxx)`, including `i(<type>.<inst>.<name>)`) wherever they
apply (R, L, C, V, B, G, E), since those never needed `save` and were already
validated; `save`d `@device[...]` vectors are now used for the two device
kinds that have no such native branch current and whose own model equations
are impractical to re-derive in Python: diodes (`@name[id]`, `save`d) and
MOSFETs. MOSFETs get a THIRD primitive instead of `@name[id/is/ig/ib]`: BSIM4
(level=54) does not expose `is`/`ib`/`p` in its own OP-parameter table at all
(only `id`, `ibs`, `ibd`, checked directly against ngspice's own device
source, `src/spicelib/devices/bsim4/b4.c`), so relying on `@device[...]`
would silently miss gate/bulk/body-diode current for that model. Instead,
every top-level MOSFET's 4 terminals (drain, gate, source, bulk) are each
wired through an inserted 0V sense voltage source in a COPY of the deck
(`inject_mosfet_sense`, never the file on disk), and each terminal's power is
`V(original node) * i(that sense source)` -- the native branch-current
primitive again, just applied per-terminal instead of per-2-terminal-device,
and correct for ANY MOSFET compact model since it never depends on what that
model's own OP table happens to expose. See docs/loss.md for the full
writeup, the reproducers and the validation numbers.

Exit status: 0 = ran and computed a complete, trustworthy breakdown; 2 = could
not run/compute at all, OR ran but at least one device in the deck could not
be accounted (kind not handled, unresolved subcircuit, MOSFET inside a
.subckt, ...) -- report.txt/summary.json still list what WAS computed and,
under "excluded", every device that was not, but the run is never reported as
a "pass" in that case (card #481: an incomplete energy-balance sum used to be
silently reported as if it were complete).
Outputs (under ~/Work/data/ngspice_augmentation/loss/, never the repository):
  wrapper.cir, ngspice.log, trace.txt (raw wrdata), summary.json, report.txt
"""
import argparse
import datetime
import json
import math
import os
import re
import subprocess
import sys

NGSPICE = os.environ.get("NGSPICE", "ngspice")
NGSPICE_AUG = os.environ.get(
    "NGSPICE_AUG",
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-pss9/bin/ngspice-aug"),
)
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/loss")

DROP_DOT = (".ac", ".tran", ".op", ".dc", ".noise", ".tf", ".pz", ".sens", ".save",
            ".end", ".four", ".meas", ".measure", ".pss", ".stb", ".pac", ".pstb")

SUFFIX = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "m": 1e-3,
          "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}

TYPEMAP = {"r": "R", "l": "L", "c": "C", "v": "V", "i": "I", "b": "B", "g": "G", "e": "E",
           "d": "D", "m": "M"}
# device kinds whose OWN branch current ngspice tracks as a native MNA unknown
# (reliable, updates every accepted timestep, accessed via i(name)).
NATIVE_I_TYPES = ("V", "L")


class CouldNotRun(RuntimeError):
    pass


# ---------------------------------------------------------------------
# deck text: strip the original control/analysis (like tools/stb.py's
# prepare_deck), make include paths absolute. Self-contained copy, not an
# import of tools/stb.py, so this file has no dependency on it.

def prepare_deck(path):
    base = os.path.dirname(os.path.abspath(path))
    with open(path) as fh:
        lines = fh.read().splitlines()
    out, in_control = [], False
    for n, line in enumerate(lines):
        low = line.strip().lower()
        if n and low.startswith(".control"):
            in_control = True
            continue
        if in_control:
            if low.startswith(".endc"):
                in_control = False
            continue
        if n and low.split(" ")[0] in DROP_DOT:
            continue
        m = re.match(r"^(\s*\.(?:include|inc|lib)\s+)(['\"]?)([^'\"\s]+)(['\"]?)(.*)$", line, re.I)
        if n and m and not os.path.isabs(m.group(3)):
            line = f"{m.group(1)}{os.path.join(base, m.group(3))}{m.group(5)}"
        out.append(line)
    return out


def merge_continuations(lines):
    """Merge SPICE line continuations (a physical line starting with `+`) into
    the previous line, so every later pass sees one logical line per card.
    Needed for any multi-line `.model` card (this project's own BSIM4
    fixtures use these, card #481) -- without this, a continuation line's own
    first token ("+") used to be silently dropped as "not a recognized device
    kind"; since enumerate_devices() no longer silently drops ANY unrecognized
    line (card #481, "never silently skip"), it would otherwise be wrongly
    reported as 8 separate excluded devices named "+"."""
    out = []
    for line in lines:
        if line.strip().startswith("+") and out:
            out[-1] = out[-1] + " " + line.strip()[1:].strip()
        else:
            out.append(line)
    return out


def included_files(lines):
    for line in lines:
        m = re.match(r"^\s*\.(?:include|inc|lib)\s+['\"]?([^'\"\s]+)['\"]?", line, re.I)
        if m and os.path.exists(m.group(1)):
            yield m.group(1)


# ---------------------------------------------------------------------
# .param resolution: SPICE-suffix numbers, {ref} substitution, memoized.

def spice_num(tok):
    m = re.match(r"^\s*([-+0-9.eE]+)(meg|[tgkmunpf])?\s*$", tok.strip().lower())
    if not m:
        return None
    val = float(m.group(1))
    if m.group(2):
        val *= SUFFIX[m.group(2)]
    return val


def parse_params(lines):
    """{name: raw_expr_text} from every `.param name=expr ...` line, braces stripped."""
    params = {}
    for line in lines:
        s = line.strip()
        if not re.match(r"^\.param\b", s, re.I):
            continue
        body = s[len(".param"):].strip()
        for name, expr in re.findall(r"(\w+)\s*=\s*(\{[^{}]*\}|\S+)", body, re.I):
            if expr.startswith("{") and expr.endswith("}"):
                expr = expr[1:-1]
            params[name.lower()] = expr
    return params


def resolve_param(name, params, cache, stack=()):
    name = name.lower()
    if name in cache:
        return cache[name]
    if name not in params:
        raise CouldNotRun(f".param {name} not found")
    if name in stack:
        raise CouldNotRun(f".param cycle: {' -> '.join(stack)} -> {name}")
    val = resolve_expr(params[name], params, cache, stack + (name,))
    cache[name] = val
    return val


def resolve_expr(expr, params, cache, stack=()):
    """A .param value or a {...}-braced expression -> float, resolving nested {refs}."""
    def repl(m):
        return f"({resolve_param(m.group(1), params, cache, stack)!r})"
    text = re.sub(r"\{([A-Za-z_]\w*)\}", repl, expr)
    text = re.sub(r"\{([A-Za-z_]\w*)\}", repl, text)  # second pass for refs revealed by the first
    # bare param-name tokens (no braces) inside an expr, e.g. "1/fs"
    def repl_bare(m):
        tok = m.group(0)
        if tok.lower() in params:
            return f"({resolve_param(tok, params, cache, stack)!r})"
        return tok
    text = re.sub(r"[A-Za-z_]\w*", repl_bare, text)
    text = re.sub(r"(?<![0-9eE.])(\d+\.?\d*(?:[eE][-+]?\d+)?)(meg|[tgkmunpf])\b",
                  lambda m: repr(spice_num(m.group(0))), text, flags=re.I)
    try:
        return float(eval(text, {"__builtins__": {}}, {}))
    except Exception as e:
        raise CouldNotRun(f"cannot evaluate expr {expr!r} (-> {text!r}): {e}") from e


# ---------------------------------------------------------------------
# device enumeration (top-level + one level of .subckt expansion). Every
# device line is enumerated as a Device, REGARDLESS of whether this tool
# knows how to compute its power (card #481: enumerate_devices used to
# `continue` past any line whose type letter wasn't in TYPEMAP -- M, D, J,
# an XSPICE `A` block, etc. were invisible to the tool, not even reported as
# excluded. Now every line becomes a Device; build_plan() is what decides
# "unsupported", and main()/render_report() refuse to call the run a pass
# when anything ends up there -- see the module docstring's Exit status).

class Device:
    def __init__(self, name, kind, n1, n2, fields, qual=None):
        self.name = name          # display name (original case)
        self.kind = kind          # R L C V I B G E D M, or any other SPICE type letter
        self.n1, self.n2 = n1, n2
        self.fields = fields      # full split netlist line (subckt-local for expanded devices)
        self.qual = qual or name.lower()  # ngspice query name: plain, or "t.inst.name"

    def __repr__(self):
        return f"<{self.kind} {self.name} {self.n1}-{self.n2}>"


def load_subckts(text):
    subs, cur, curname = {}, None, None
    for line in text.splitlines():
        s = line.strip()
        m = re.match(r"^\.subckt\s+(\S+)\s+(.*)$", s, re.I)
        if m:
            curname = m.group(1).lower()
            cur = {"ports": m.group(2).split(), "devs": []}
            subs[curname] = cur
            continue
        if re.match(r"^\.ends\b", s, re.I):
            cur = None
            continue
        if cur is not None and s and not s.startswith("*"):
            cur["devs"].append(s)
    return subs


def enumerate_devices(lines, subs):
    devices = []
    for lineno, line in enumerate(lines):
        if lineno == 0:
            # card #511: the first line is the title (ngspice skips it); a title
            # such as "Buck ..." used to become a phantom B device, silently.
            continue
        s = line.strip()
        if not s or s.startswith("*") or s.startswith("."):
            continue
        f = s.split()
        name = f[0]
        t = name[0].lower()
        if t == "x":
            subname = f[-1].lower()
            if subname not in subs:
                # unresolved subckt: still record it (never silently drop it) so it
                # shows up under "excluded" and blocks a pass, rather than vanishing.
                devices.append(Device(name, "X", f[1] if len(f) > 1 else "0",
                                       f[2] if len(f) > 2 else "0", f, qual=name.lower()))
                continue
            sub = subs[subname]
            nodes = f[1:-1]
            portmap = dict(zip([p.lower() for p in sub["ports"]], nodes))

            def mapnode(nd, inst=name):
                low = nd.lower()
                if low == "0":
                    return "0"
                if low in portmap:
                    return portmap[low]
                return f"{inst.lower()}.{low}"

            for dl in sub["devs"]:
                df = dl.split()
                dname = df[0]
                dt = dname[0].lower()
                if len(df) < 3:
                    continue  # not even 2 nodes to parse; not a real device line
                qual = f"{dt}.{name.lower()}.{dname.lower()}"
                kind = TYPEMAP.get(dt, dname[0].upper())
                devices.append(Device(f"{name}.{dname}", kind,
                                       mapnode(df[1]), mapnode(df[2]), df, qual=qual))
            continue
        if len(f) < 3:
            continue  # not even 2 nodes to parse; not a real device line
        kind = TYPEMAP.get(t, name[0].upper())
        n1, n2 = f[1], f[2]
        if kind == "M" and len(f) >= 5:
            n1, n2 = f[1], f[3]  # drain, source (skip the gate) -- the pair this
            # tool treats as "the blocking terminals" for classification and the
            # switch conduction/transition split (Vds), same role n1/n2 play for
            # every other 2-terminal device kind.
        devices.append(Device(name, kind, n1, n2, f, qual=name.lower()))
    return devices


def fanout(devices):
    """node (lowercase) -> number of device terminal slots (n1/n2) touching it."""
    fo = {}
    for d in devices:
        for nd in (d.n1, d.n2):
            if nd == "0":
                continue
            fo[nd] = fo.get(nd, 0) + 1
    return fo


# ---------------------------------------------------------------------
# MOSFET terminal-current instrumentation: rewrite every TOP-LEVEL `M` device
# line so each of its 4 terminals (drain, gate, source, bulk) is measured by
# an inserted 0V sense voltage source, in a COPY of the deck text -- never
# the file on disk, same "wrapper, not an edit" pattern build_wrapper() uses
# for the control block. See the module docstring for why: BSIM4 (level=54)
# has no is/ib/p in its own OP-parameter table (checked directly against
# ngspice's device source), so @name[id/is/ig/ib] cannot give a trustworthy
# total device power for every MOSFET level this tool is asked to support --
# but i(Vsense), the native MNA branch current of an inserted 0V source, is
# exactly the same reliable primitive this tool already trusts for V/L
# devices, just applied once per terminal instead of once per 2-terminal
# device. P(t) for the MOSFET is then sum over d/g/s/b of
# V(original node) * i(sense source at that terminal) -- valid for ANY
# compact model (it never reads a single model-specific OP parameter).

def inject_mosfet_sense(lines):
    """Return (new_lines, sense). sense[Mname] = {"d": (orig_node, isig), "g": ...,
    "s": ..., "b": ...}. Only rewrites M lines OUTSIDE any .subckt block (depth
    tracked here); an M line inside a .subckt is left untouched and build_plan()
    reports it as unsupported rather than silently measuring it wrong."""
    out = []
    sense = {}
    depth = 0
    for lineno, line in enumerate(lines):
        s = line.strip()
        low = s.lower()
        if lineno == 0:
            # the netlist's first line is its title, never a card (card #511)
            out.append(line)
            continue
        if low.startswith(".subckt"):
            depth += 1
            out.append(line)
            continue
        if low.startswith(".ends"):
            depth = max(0, depth - 1)
            out.append(line)
            continue
        if depth == 0 and s and not s.startswith("*") and not s.startswith("."):
            f = s.split()
            if f[0][0].lower() == "m" and len(f) >= 5:
                name = f[0]
                tag = name.lower()
                orig = {"d": f[1], "g": f[2], "s": f[3], "b": f[4]}
                new_node = {t: f"__lsn{t}_{tag}" for t in ("d", "g", "s", "b")}
                rest = f[5:]
                out.append(" ".join([name, new_node["d"], new_node["g"],
                                      new_node["s"], new_node["b"]] + rest))
                terms = {}
                for t in ("d", "g", "s", "b"):
                    vname = f"V__ls{t}_{name}"
                    out.append(f"{vname} {orig[t]} {new_node[t]} DC 0")
                    terms[t] = (orig[t], f"i({vname.lower()})")
                sense[name] = terms
                continue
        out.append(line)
    return out, sense


# ---------------------------------------------------------------------
# B/G controlling expressions -> a Python-evaluable callable of node voltages.

def b_expr_text(fields):
    """The RHS expression text of `Bxxx n1 n2 I=expr` / `V=expr`, or None."""
    joined = " ".join(fields[3:])
    m = re.match(r"^\s*[IV]\s*=\s*(.*)$", joined, re.I)
    if not m:
        return (None, None)
    expr = m.group(1).strip()
    # card #511: ngspice accepts I='expr' / I="expr"; the quotes are not part of it
    if len(expr) >= 2 and expr[0] == expr[-1] and expr[0] in "'\"":
        expr = expr[1:-1].strip()
    return (expr, joined.strip().upper().startswith("V"))


def expr_nodes(expr):
    """distinct node names referenced as V(node) or V(n1,n2) in expr, lowercased."""
    nodes = set()
    for m in re.finditer(r"V\(\s*([^,()\s]+)\s*(?:,\s*([^()\s]+)\s*)?\)", expr, re.I):
        nodes.add(m.group(1).lower())
        if m.group(2):
            nodes.add(m.group(2).lower())
    return nodes


def compile_expr(expr, params, cache):
    """expr text (with {param} refs and V(...)/V(...,...) reads) -> python source
    string using a local dict `V` keyed by lowercase node name, and math fns."""
    text = expr

    def repl_param(m):
        return f"({resolve_param(m.group(1), params, cache)!r})"
    text = re.sub(r"\{([A-Za-z_]\w*)\}", repl_param, text)

    def repl_v(m):
        n1 = m.group(1).lower()
        if m.group(2):
            n2 = m.group(2).lower()
            return f"(V[{n1!r}]-V[{n2!r}])"
        return f"V[{n1!r}]"
    text = re.sub(r"V\(\s*([^,()\s]+)\s*(?:,\s*([^()\s]+)\s*)?\)", repl_v, text, flags=re.I)
    text = text.replace("**", "@@POW@@").replace("^", "**").replace("@@POW@@", "**")
    # card #511: ngspice B-source logic -- &&, ||, ! and the ternary c ? a : b
    text = text.replace("&&", " and ").replace("||", " or ")
    text = re.sub(r"!(?!=)", " not ", text)
    return _ternary_to_py(text)


def _ternary_to_py(text):
    """c ? a : b (right-associative, any nesting, parentheses respected)
    -> (a if c else b)."""
    depth, q = 0, -1
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "?" and depth == 0:
            q = i
            break
    if q < 0:
        # no top-level ternary: recurse into parenthesised groups
        out, i = [], 0
        while i < len(text):
            if text[i] == "(":
                d, j = 1, i + 1
                while j < len(text) and d:
                    d += {"(": 1, ")": -1}.get(text[j], 0)
                    j += 1
                out.append("(" + _ternary_to_py(text[i + 1:j - 1]) + ")")
                i = j
            else:
                out.append(text[i])
                i += 1
        return "".join(out)
    depth, nest, c = 0, 0, -1
    for i in range(q + 1, len(text)):
        ch = text[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif depth == 0 and ch == "?":
            nest += 1
        elif depth == 0 and ch == ":":
            if nest == 0:
                c = i
                break
            nest -= 1
    if c < 0:
        raise SyntaxError(f"ternary '?' without ':' in B expression: {text!r}")
    cond = _ternary_to_py(text[:q])
    a = _ternary_to_py(text[q + 1:c])
    b = _ternary_to_py(text[c + 1:])
    return f"(({a}) if ({cond}) else ({b}))"


MATH_NS = {k: getattr(math, k) for k in ("tanh", "exp", "sqrt", "sin", "cos", "log", "atan", "pi")}
MATH_NS["abs"] = abs


def make_evaluator(text):
    code = compile(text, "<b_expr>", "eval")

    def f(V):
        ns = dict(MATH_NS)
        ns["V"] = V
        return eval(code, {"__builtins__": {}}, ns)
    return f


# ---------------------------------------------------------------------
# wrdata reader: N pairs (t, val) per row, in request order.

def read_wrdata(path, n):
    rows = []
    with open(path) as fh:
        for line in fh:
            v = line.split()
            if len(v) != 2 * n:
                continue
            rows.append([float(x) for x in v])
    if not rows:
        raise CouldNotRun(f"{path} has no usable rows")
    return rows


def read_wrdata_window(path, n, t0, keep_margin=1):
    """Stream `path`, keeping only rows with time >= t0 (plus up to `keep_margin`
    rows just before it, for a clean first interval). Keeps memory bounded to the
    window, not the whole tran run -- important for R4-length (6ms) decks (found
    2026-09-29, card #470: loading the full ~6ms trace -- 1.5M+ rows -- into a
    Python list used several GB and made post-processing itself the dominant
    cost, not ngspice; this streaming filter fixes it)."""
    import collections
    pre = collections.deque(maxlen=keep_margin)
    kept = []
    saw_any = False
    with open(path) as fh:
        for line in fh:
            v = line.split()
            if len(v) != 2 * n:
                continue
            row = [float(x) for x in v]
            saw_any = True
            if row[0] >= t0:
                kept.append(row)
            else:
                pre.append(row)
    if not saw_any:
        raise CouldNotRun(f"{path} has no usable rows")
    return list(pre) + kept


# ---------------------------------------------------------------------
# signal-plan: decide, per device, what to fetch and how to turn it into P(t).

class Plan:
    def __init__(self):
        self.i_sigs = []     # native i(name) branch currents -- never need `save` (card #481)
        self.at_sigs = []    # @dev[param] device-internal vectors -- DO need an explicit
                              # `save` issued before the analysis, or they freeze at t=0
                              # (card #481; see the module docstring)
        self.node_sigs = []  # node names (lowercase) to request as v(node)
        self._seen_i, self._seen_at, self._seen_v = set(), set(), set()

    def want_i(self, qual):
        key = f"i({qual})"
        if key not in self._seen_i:
            self._seen_i.add(key)
            self.i_sigs.append(key)
        return key

    def want_i_raw(self, sig):
        """Register an already-formatted `i(...)` signal (e.g. from a sense source
        this tool inserted itself, `inject_mosfet_sense`) -- same native-current
        guarantee as want_i(), just not built from a Device's own qual."""
        if sig not in self._seen_i:
            self._seen_i.add(sig)
            self.i_sigs.append(sig)
        return sig

    def want_at(self, sig):
        """Register an `@dev[param]` vector. Callers must NOT rely on this tracking
        per-timestep unless it ends up in a `save` line before the analysis --
        build_wrapper() does that for every at_sig automatically."""
        if sig not in self._seen_at:
            self._seen_at.add(sig)
            self.at_sigs.append(sig)
        return sig

    def want_v(self, node):
        if node != "0" and node not in self._seen_v:
            self._seen_v.add(node)
            self.node_sigs.append(node)
        return node


def build_plan(devices, params, cache, fan, mos_sense=None):
    """Return (plan, per_device_info). per_device_info[dev] = dict describing
    how to compute P(t) for that device from the fetched signals. mos_sense
    (from inject_mosfet_sense) maps a top-level M device's display name to its
    per-terminal (orig_node, i(sense_source)) pairs."""
    mos_sense = mos_sense or {}
    plan = Plan()
    info = {}
    for d in devices:
        if d.kind in NATIVE_I_TYPES:
            key = plan.want_i(d.qual)
            plan.want_v(d.n1); plan.want_v(d.n2)
            info[d] = {"method": "native_i", "isig": key}
        elif d.kind == "R":
            # R's value is fetched from ngspice itself (@name[resistance], a static
            # elaboration-time parameter, queried once via `print`), not parsed from
            # the netlist text, so a resistor line with extra fields (tc1=, etc.)
            # never trips this up.
            plan.want_v(d.n1); plan.want_v(d.n2)
            info[d] = {"method": "ohm", "R": None}
        elif d.kind == "C":
            plan.want_v(d.n1); plan.want_v(d.n2)
            info[d] = {"method": "cap", "C": None}
        elif d.kind in ("B", "G", "E"):
            if d.kind == "B":
                expr, is_v = b_expr_text(d.fields)
                if expr is None:
                    info[d] = {"method": "zero"}
                    continue
                if is_v:
                    if fan.get(d.n1, 0) + fan.get(d.n2, 0) <= 2:
                        # only this device's own two terminals touch these nodes
                        info[d] = {"method": "zero"}
                    else:
                        # externally loaded: its current is not in its own expression,
                        # but a V-defining B-source has a branch, so ngspice's own
                        # i(name) gives it exactly, same convention as a V source
                        info[d] = {"method": "native_i", "isig": plan.want_i(d.qual)}
                    plan.want_v(d.n1); plan.want_v(d.n2)
                    continue
                nodes = expr_nodes(expr) | {d.n1, d.n2}
                for nd in nodes:
                    plan.want_v(nd)
                fn = make_evaluator(compile_expr(expr, params, cache))
                info[d] = {"method": "bexpr", "fn": fn}
            else:
                # G/E: Gname n1 n2 nc1 nc2 gain  (simple linear VCCS/VCVS only)
                if len(d.fields) < 6:
                    info[d] = {"method": "unsupported", "why": "non-linear/short G-E line"}
                    continue
                nc1, nc2, gaintok = d.fields[3], d.fields[4], d.fields[5]
                gain = spice_num(gaintok)
                if gain is None:
                    gain = resolve_param(gaintok, params, cache) if gaintok.lower() in params \
                        else resolve_expr(gaintok.strip("{}"), params, cache)
                plan.want_v(d.n1); plan.want_v(d.n2); plan.want_v(nc1); plan.want_v(nc2)
                info[d] = {"method": "vccs", "gain": gain, "nc1": nc1.lower(), "nc2": nc2.lower()}
        elif d.kind == "I":
            # independent current source: only a constant DC 0 is supported generically
            # (lgprobe.lib's Ii), matching this project's only actual use.
            joined = " ".join(d.fields[3:]).upper()
            m = re.match(r"^\s*DC\s+([-+0-9.eE]+)", joined)
            val = float(m.group(1)) if m else None
            if val == 0.0:
                info[d] = {"method": "const", "I": 0.0}
                plan.want_v(d.n1); plan.want_v(d.n2)
            else:
                info[d] = {"method": "unsupported", "why": "non-zero/time-varying I source"}
        elif d.kind == "D":
            # native diode: no MNA branch-current unknown of its own (unlike V/L), and
            # re-deriving the exponential/junction-cap model in Python (the way B/G
            # expressions are re-derived) is impractical -- so this uses ngspice's own
            # computed terminal current, @name[id], which (once `save`d -- card #481)
            # tracks correctly per timestep and was validated to equal V(anode,cathode)
            # * id to 5+ significant figures against ngspice's own @name[p] at several
            # sampled points (docs/loss.md). Works at any .subckt nesting depth, same
            # as the static @name[resistance]/[capacitance] queries above.
            sig = plan.want_at(f"@{d.qual}[id]")
            plan.want_v(d.n1); plan.want_v(d.n2)
            info[d] = {"method": "native_i", "isig": sig}
        elif d.kind == "M":
            # MOSFET: see the module docstring for why this does NOT use @name[id/is/
            # ig/ib] (BSIM4 doesn't expose is/ib/p at all) -- instead every top-level
            # M device gets 4 inserted 0V sense sources (inject_mosfet_sense, run
            # before this function, on a COPY of the deck), one per terminal, and this
            # sums V(original terminal node) * i(that sense source) over d/g/s/b. A
            # MOSFET inside a .subckt instance is NOT covered (sense-source injection
            # only rewrites top-level lines) and is reported unsupported rather than
            # silently measured wrong or dropped.
            if "." in d.qual:
                info[d] = {"method": "unsupported",
                           "why": "MOSFET inside a .subckt instance; this tool's "
                                  "sense-source terminal-current method only rewrites "
                                  "top-level M lines"}
            elif d.name not in mos_sense:
                info[d] = {"method": "unsupported",
                           "why": "MOSFET line was not recognized for sense-source "
                                  "injection (unexpected field count)"}
            else:
                terms = mos_sense[d.name]
                for orig_node, isig in terms.values():
                    plan.want_v(orig_node)
                    plan.want_i_raw(isig)
                info[d] = {"method": "mosfet_sense", "terms": terms}
        else:
            info[d] = {"method": "unsupported", "why": f"unhandled device kind {d.kind}"}
    return plan, info


# ---------------------------------------------------------------------
# run ngspice, build the wrdata request, parse the static R/C queries.

def run_ngspice(binary, wrapper, timeout):
    try:
        proc = subprocess.run([binary, "-b", wrapper], capture_output=True, text=True, timeout=timeout)
    except (FileNotFoundError, PermissionError) as e:
        raise CouldNotRun(f"cannot start {binary}: {e}") from e
    except subprocess.TimeoutExpired as e:
        raise CouldNotRun(f"{binary} timed out after {timeout}s") from e
    return proc.returncode, proc.stdout + proc.stderr


def build_wrapper(lines, plan, static_r, static_c, trace_path, analysis_ctl):
    static_prints = [f"print @{name}[resistance]" for name in static_r] + \
                    [f"print @{name}[capacitance]" for name in static_c]
    sigs = list(plan.i_sigs) + list(plan.at_sigs) + [f"v({nd})" for nd in plan.node_sigs]
    ctl = "\n.control\nset filetype=ascii\nset numdgt=17\n"
    ctl += "\n".join(static_prints) + "\n" if static_prints else ""
    if plan.at_sigs:
        # MUST come before the analysis line: `@dev[param]` vectors only get their
        # per-timestep history recorded if `save`d before the run starts (card #481,
        # see the module docstring) -- issuing this after `tran`/`pss` would be too
        # late, same freeze bug this line exists to avoid. Crucially, `save` is NOT
        # additive to ngspice's own default "no .save/save at all -> save every
        # ordinary vector" behaviour -- issuing ANY `save` command REPLACES that
        # default with exactly the named list (found empirically, card #481: a
        # `save @dev[id]` line alone made even i(l1)/v(node) -- normally saved for
        # free -- unavailable to the later wrdata). So the save list here is the
        # FULL signal request (i_sigs + at_sigs + node voltages), not just at_sigs.
        ctl += "save " + " ".join(sigs) + "\n"
    ctl += analysis_ctl
    ctl += f"\nwrdata {trace_path} " + " ".join(sigs) + "\n.endc\n.end\n"
    return "\n".join(lines) + ctl, sigs


def parse_static(out, names, param):
    vals = {}
    for name in names:
        m = re.search(rf"@{re.escape(name)}\[{param}\]\s*=\s*([-+0-9.eE]+)", out, re.I)
        if not m:
            raise CouldNotRun(f"ngspice did not report @{name}[{param}]")
        vals[name] = float(m.group(1))
    return vals


# ---------------------------------------------------------------------
# category classification (deck-agnostic defaults, overridable via CLI).

CATEGORIES = ("switch_conduction", "switch_transition", "gate_drive",
              "diode_conduction", "inductor_dcr", "capacitor_esr", "resistive_other")


def classify(devices, vin_name, load_name, dcr_names, esr_names, switch_names, gate_names,
             diode_names=()):
    lower = lambda s: s.lower()
    dcr, esr, sw, gate, diode = (set(map(lower, x)) for x in
                                  (dcr_names, esr_names, switch_names, gate_names, diode_names))
    vin_dev = load_dev = None
    cat = {}
    for d in devices:
        nl = d.name.lower()
        if nl == vin_name.lower() and vin_dev is None and "." not in d.name:
            vin_dev = d
            continue
        if nl == load_name.lower() and load_dev is None and "." not in d.name:
            load_dev = d
            continue
        if nl in sw:
            cat[d] = "switch"
        elif nl in gate:
            cat[d] = "gate_drive"
        elif nl in dcr:
            cat[d] = "inductor_dcr"
        elif nl in esr:
            cat[d] = "capacitor_esr"
        elif nl in diode:
            cat[d] = "diode_conduction"
        else:
            cat[d] = "resistive_other"
    if vin_dev is None:
        raise CouldNotRun(f"input source {vin_name!r} not found in deck")
    if load_dev is None:
        raise CouldNotRun(f"load {load_name!r} not found in deck")
    return vin_dev, load_dev, cat


# ---------------------------------------------------------------------
# power(t) evaluation and trapezoidal/interval integration over a window.

def device_power_series(d, dinfo, rows, idx, t):
    """[(t_i, P_i)] pointwise for native/ohm/bexpr/vccs/const/zero devices, using
    node-voltage columns already resolved via idx (name -> column offset into
    each row's value list, values at rows[k][idx[name]])."""
    method = dinfo["method"]
    out = []
    for r in rows:
        if method == "native_i":
            I = r[idx[dinfo["isig"]]]
            v1 = 0.0 if d.n1 == "0" else r[idx[d.n1]]
            v2 = 0.0 if d.n2 == "0" else r[idx[d.n2]]
            out.append((v1 - v2) * I)
        elif method == "ohm":
            v1 = 0.0 if d.n1 == "0" else r[idx[d.n1]]
            v2 = 0.0 if d.n2 == "0" else r[idx[d.n2]]
            R = dinfo["R"]
            out.append((v1 - v2) ** 2 / R if R else 0.0)
        elif method == "bexpr":
            V = {nd: (0.0 if nd == "0" else r[idx[nd]]) for nd in idx
                 if not nd.startswith("i(") and not nd.startswith("@")}
            I = dinfo["fn"](V)
            v1 = 0.0 if d.n1 == "0" else r[idx[d.n1]]
            v2 = 0.0 if d.n2 == "0" else r[idx[d.n2]]
            out.append((v1 - v2) * I)
        elif method == "mosfet_sense":
            p = 0.0
            for orig_node, isig in dinfo["terms"].values():
                v = 0.0 if orig_node == "0" else r[idx[orig_node]]
                I = r[idx[isig]]
                p += v * I
            out.append(p)
        elif method == "vccs":
            v1 = 0.0 if d.n1 == "0" else r[idx[d.n1]]
            v2 = 0.0 if d.n2 == "0" else r[idx[d.n2]]
            vc1 = 0.0 if dinfo["nc1"] == "0" else r[idx[dinfo["nc1"]]]
            vc2 = 0.0 if dinfo["nc2"] == "0" else r[idx[dinfo["nc2"]]]
            I = dinfo["gain"] * (vc1 - vc2)
            out.append((v1 - v2) * I)
        elif method == "const":
            v1 = 0.0 if d.n1 == "0" else r[idx[d.n1]]
            v2 = 0.0 if d.n2 == "0" else r[idx[d.n2]]
            out.append((v1 - v2) * dinfo["I"])
        else:  # zero / unsupported (unsupported handled by caller, excluded upstream)
            out.append(0.0)
    return out


def cap_energy(d, dinfo, rows, idx):
    """Capacitor energy over the window via Q=C*V, dE = V_avg * dQ per interval
    (exact for the trapezoidal V(t) this implies; avoids needing dV/dt directly)."""
    C = dinfo["C"]
    if not C:
        return 0.0
    vs = []
    ts = []
    for r in rows:
        v1 = 0.0 if d.n1 == "0" else r[idx[d.n1]]
        v2 = 0.0 if d.n2 == "0" else r[idx[d.n2]]
        vs.append(v1 - v2)
        ts.append(r[0])
    e = 0.0
    for k in range(len(vs) - 1):
        dq = C * (vs[k + 1] - vs[k])
        e += 0.5 * (vs[k] + vs[k + 1]) * dq
    return e


def trapz_energy(series_t, series_p):
    e = 0.0
    for k in range(len(series_t) - 1):
        e += 0.5 * (series_p[k] + series_p[k + 1]) * (series_t[k + 1] - series_t[k])
    return e


def switch_split_energy(d, dinfo, rows, idx, low_frac, high_frac):
    """(conduction_energy, transition_energy) for a switch device by the
    documented criterion: at each interval, classify by the *midpoint*
    fraction of |V(n1)-V(n2))| relative to that switch's own max over the
    window ("the rail" it blocks). Below low_frac*rail or above
    high_frac*rail -> conduction (on, or off with negligible current either
    way); in between -> transition/overlap."""
    method = dinfo["method"]
    vdrop, p = [], []
    for r in rows:
        v1 = 0.0 if d.n1 == "0" else r[idx[d.n1]]
        v2 = 0.0 if d.n2 == "0" else r[idx[d.n2]]
        vdrop.append(v1 - v2)
    series = device_power_series(d, dinfo, rows, idx, None)
    rail = max(abs(x) for x in vdrop) or 1.0
    ts = [r[0] for r in rows]
    e_cond = e_trans = 0.0
    for k in range(len(ts) - 1):
        mid = 0.5 * (abs(vdrop[k]) + abs(vdrop[k + 1])) / rail
        dt = ts[k + 1] - ts[k]
        seg_e = 0.5 * (series[k] + series[k + 1]) * dt
        if low_frac <= mid <= high_frac:
            e_trans += seg_e
        else:
            e_cond += seg_e
    return e_cond, e_trans


# ---------------------------------------------------------------------

def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("--vin", required=True, help="input voltage source instance name")
    ap.add_argument("--load", required=True, help="output load resistor instance name")
    ap.add_argument("--dcr", action="append", default=[], help="resistor name = inductor DCR (repeatable)")
    ap.add_argument("--esr", action="append", default=[], help="resistor name = capacitor ESR (repeatable)")
    ap.add_argument("--switch", action="append", default=[],
                     help="B-source name = a power switch, split into conduction/transition (repeatable)")
    ap.add_argument("--gate-drive", action="append", default=[], dest="gate_drive",
                     help="device name = gate-drive loss element, if modelled (repeatable)")
    ap.add_argument("--diode", action="append", default=[],
                     help="D device name = a power diode, bucketed into diode_conduction "
                          "(includes reverse-recovery -- see docs/loss.md) instead of "
                          "resistive_other (repeatable)")
    ap.add_argument("--switch-low", type=float, default=0.10,
                     help="below this fraction of the switch's own max |v| -> conduction (default 0.10)")
    ap.add_argument("--switch-high", type=float, default=0.90,
                     help="above this fraction of the switch's own max |v| -> conduction, i.e. off (default 0.90)")
    ap.add_argument("--period", type=float, required=True, help="switching period, seconds")

    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--tstop", type=float, help="tran path: total run time, seconds")
    mode.add_argument("--pss", action="store_true", help="pss path: use `.pss ... method=shoot` on NGSPICE_AUG")
    ap.add_argument("--periods", type=float, default=5,
                     help="tran path: window = last N periods (default 5; N need not be an "
                          "integer -- used by the red-proof test to show a non-integer-period "
                          "window breaks the energy balance)")
    ap.add_argument("--tstep", default=None, help="tran path: fixed max timestep (default period/2000)")
    ap.add_argument("--fguess", type=float, help="pss path: Fguess (default 1/period)")
    ap.add_argument("--stabtime", type=float, help="pss path: StabTime (default 150*period)")
    ap.add_argument("--oscnode", help="pss path: OscNode")
    ap.add_argument("--points", type=int, default=200, help="pss path: Points (default 200)")
    ap.add_argument("--harmonics", type=int, default=10, help="pss path: Harmonics (default 10)")
    ap.add_argument("--sciter", type=int, default=150, help="pss path: SC_iter (default 150)")
    ap.add_argument("--steady", type=float, default=1e-6, help="pss path: Steady_coeff (default 1e-6)")
    ap.add_argument("--ngspice-aug", default=NGSPICE_AUG, help="fork ngspice-aug binary for --pss")

    ap.add_argument("--out", help="output directory (default: under ~/Work/data/ngspice_augmentation/loss)")
    ap.add_argument("--timeout", type=float, default=1800)
    args = ap.parse_args(argv)

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = args.out or os.path.join(
        DATA_ROOT, f"{os.path.splitext(os.path.basename(args.deck))[0]}-{stamp}")
    os.makedirs(outdir, exist_ok=True)

    try:
        result = run_breakdown(args, outdir)
    except (CouldNotRun, subprocess.TimeoutExpired) as e:
        print(f"COULD NOT RUN: {e}", file=sys.stderr)
        return 2

    write_outputs(result, outdir)
    print(render_report(result))
    print(f"results: {outdir}")
    if result["excluded"]:
        names = [f"{e['name']}({e['kind']})" for e in result["excluded"]]
        print(f"COULD NOT RUN: energy balance refused -- {len(names)} unsupported "
              f"device(s) not included in the sum: {names}", file=sys.stderr)
        return 2
    return 0


def run_breakdown(args, outdir):
    lines = prepare_deck(args.deck)
    lines = merge_continuations(lines)
    # MOSFET terminal-current instrumentation (card #481): rewrites top-level M
    # lines in THIS in-memory copy only, before anything else reads `lines` --
    # params/subckt parsing below see only ordinary R/V/M lines either way, so
    # this has to happen first, not as a build_plan-time afterthought.
    lines, mos_sense = inject_mosfet_sense(lines)
    params = parse_params(lines)
    cache = {}
    subs = {}
    for inc in included_files(lines):
        subs.update(load_subckts("\n".join(merge_continuations(open(inc).read().splitlines()))))
    devices = enumerate_devices(lines, subs)
    # fanout() must see the REWRITTEN M lines' internal node names here (a top-
    # level M's own n1/n2 no longer touch the original nodes post-injection --
    # only its 4 sense sources do, one terminal-slot each, exactly replacing
    # the M device's own former contribution there; computing fanout after the
    # node fix-up below would double-count each MOSFET terminal node).
    fan = fanout(devices)
    # NOW fix up a top-level M device's own n1/n2 (drain, source) from
    # inject_mosfet_sense's internal __lsnd_*/__lsns_*/... junction names back
    # to the ORIGINAL circuit nodes -- what the rest of this tool (e.g.
    # switch_split_energy's Vds computation) actually has voltage data for,
    # since build_plan's M dispatch calls want_v() on the original nodes only.
    for d in devices:
        if d.kind == "M" and d.name in mos_sense:
            d.n1 = mos_sense[d.name]["d"][0]
            d.n2 = mos_sense[d.name]["s"][0]

    vin_dev, load_dev, cat = classify(devices, args.vin, args.load, args.dcr, args.esr,
                                       args.switch, args.gate_drive, args.diode)
    # the sense sources inject_mosfet_sense() added are internal instrumentation,
    # not real deck elements -- they measure exactly 0 W themselves (ideal 0V
    # sources) and are folded into their MOSFET's own "mosfet_sense" power, so
    # they are excluded here rather than double-counted or listed as noise.
    other_devices = [d for d in devices if d not in (vin_dev, load_dev)
                      and not d.name.lower().startswith("v__ls")]
    all_for_plan = [vin_dev, load_dev] + other_devices
    plan, info = build_plan(all_for_plan, params, cache, fan, mos_sense)

    trace_path = os.path.join(outdir, "trace.txt")
    static_r = sorted({d.qual for d, dinfo in info.items() if dinfo["method"] == "ohm"})
    static_c = sorted({d.qual for d, dinfo in info.items() if dinfo["method"] == "cap"})

    if args.pss:
        fguess = args.fguess or (1.0 / args.period)
        stabtime = args.stabtime if args.stabtime is not None else 150 * args.period
        oscnode = args.oscnode or args.load
        analysis_ctl = (f"\npss {fguess:.9g} {stabtime:.9g} {oscnode} {args.points} "
                         f"{args.harmonics} {args.sciter} {args.steady:.3g} uic method=shoot\n")
        binary = args.ngspice_aug
    else:
        if args.tstop is None:
            raise CouldNotRun("--tstop is required for the tran path (or pass --pss)")
        tstep = args.tstep or (args.period / 1000.0)
        analysis_ctl = f"\ntran {tstep} {args.tstop} 0 uic\n"
        binary = NGSPICE

    wrapper_text, sigs = build_wrapper(lines, plan, static_r, static_c, trace_path, analysis_ctl)
    wrapper_path = os.path.join(outdir, "wrapper.cir")
    with open(wrapper_path, "w") as fh:
        fh.write(wrapper_text)

    rc, out = run_ngspice(binary, wrapper_path, args.timeout)
    with open(os.path.join(outdir, "ngspice.log"), "w") as fh:
        fh.write(out)
    errs = [ln for ln in out.splitlines() if re.search(r"\berror\b", ln, re.I)]
    if args.pss:
        # card #505: the shooting engine's own "PSS-SHOOT: ..." diagnostics quote ngspice
        # status codes ("... did not complete (ngspice error 106)") for RECOVERABLE
        # events -- a probe warm-up stage that was not enough, a monodromy column that
        # fell back to a one-sided difference -- and those used to be misread here as a
        # fatal ngspice error (loss.py --pss then reported COULD NOT RUN on a run whose
        # PSS-SHOOT actually converged: R3d/R4d, card #505 part d). In --pss mode the
        # verdict is the engine's own: a "PSS-SHOOT: CONVERGED" line, or the FAILED line.
        errs = [ln for ln in errs if not ln.lstrip().startswith("PSS-SHOOT:")]
        if rc == 0 and not errs and "PSS-SHOOT: CONVERGED" not in out:
            failed = [ln for ln in out.splitlines() if "PSS-SHOOT: FAILED" in ln]
            raise CouldNotRun("PSS-SHOOT did not converge: " + (failed[0] if failed else "no CONVERGED line"))
    # R6-style XSPICE duty-clamp lines are known-benign elsewhere in this project;
    # here we have no XSPICE devices, so any "error" line is treated as fatal.
    if rc != 0 or errs:
        raise CouldNotRun(f"ngspice exit {rc}; {errs[:5]}")

    static_r_vals = parse_static(out, static_r, "resistance")
    static_c_vals = parse_static(out, static_c, "capacitance")
    for d, dinfo in info.items():
        if dinfo["method"] == "ohm":
            dinfo["R"] = static_r_vals[d.qual]
        elif dinfo["method"] == "cap":
            dinfo["C"] = static_c_vals[d.qual]

    n_sig = len(sigs)
    idx = {}
    for k, s in enumerate(sigs):
        key = s[2:-1] if s.startswith("v(") else s  # "v(node)" -> node name; "i(...)" kept whole
        idx[key] = 2 * k + 1
    # time column: use the first requested signal's own time column (all share one grid)
    t_idx = 0

    if args.pss:
        rows = read_wrdata(trace_path, n_sig)
        window_rows = rows
        window_t0, window_t1 = rows[0][t_idx], rows[-1][t_idx]
        n_periods = 1
    else:
        # stream the file, keeping only the trailing window -- the tran run can be
        # ms-long at ns-scale steps (millions of rows); reading it all into memory
        # made post-processing the dominant cost, not ngspice (card #470 finding).
        window_t0 = args.tstop - args.periods * args.period
        window_rows = read_wrdata_window(trace_path, n_sig, window_t0)
        window_t1 = window_rows[-1][t_idx]
        if len(window_rows) < 4:
            raise CouldNotRun(f"window [{window_t0},{window_t1}] has only {len(window_rows)} samples")
        n_periods = args.periods

    tw0, tw1 = window_rows[0][t_idx], window_rows[-1][t_idx]
    T = tw1 - tw0
    if T <= 0:
        raise CouldNotRun("empty/degenerate window")

    # rebuild each row as (t, value_at_idx[name]) keyed lookup done by position;
    # `idx` maps name -> column offset directly into each row list.
    def col(r, name):
        return r[idx[name]]

    # device power series
    per_device_avg = {}
    switch_split = {}
    for d, dinfo in info.items():
        if dinfo["method"] == "unsupported":
            continue
        if dinfo["method"] == "cap":
            e = cap_energy(d, dinfo, window_rows, idx)
            per_device_avg[d] = e / T
            continue
        series = device_power_series(d, dinfo, window_rows, idx, None)
        ts = [r[t_idx] for r in window_rows]
        e = trapz_energy(ts, series)
        per_device_avg[d] = e / T
        if cat.get(d) == "switch":
            e_cond, e_trans = switch_split_energy(d, dinfo, window_rows, idx,
                                                    args.switch_low, args.switch_high)
            switch_split[d] = (e_cond / T, e_trans / T)

    p_in = -per_device_avg[vin_dev]
    p_out = per_device_avg[load_dev]

    totals = {c: 0.0 for c in CATEGORIES}
    per_device_report = []
    for d in other_devices:
        dinfo = info[d]
        if dinfo["method"] == "unsupported":
            per_device_report.append({"name": d.name, "kind": d.kind, "category": "EXCLUDED",
                                       "avg_power_w": None, "note": dinfo["why"]})
            continue
        p = per_device_avg[d]
        c = cat[d]
        if c == "switch":
            e_cond, e_trans = switch_split.get(d, (p, 0.0))
            totals["switch_conduction"] += e_cond
            totals["switch_transition"] += e_trans
            per_device_report.append({"name": d.name, "kind": d.kind, "category": "switch",
                                       "avg_power_w": p, "conduction_w": e_cond, "transition_w": e_trans})
        else:
            totals[c] += p
            per_device_report.append({"name": d.name, "kind": d.kind, "category": c, "avg_power_w": p})

    sum_losses = sum(totals.values())
    residual = p_in - (p_out + sum_losses)
    residual_pct = 100.0 * residual / p_in if p_in else float("nan")
    efficiency = p_out / p_in if p_in else float("nan")

    excluded = [r for r in per_device_report if r["category"] == "EXCLUDED"]

    return {
        "deck": args.deck, "mode": "pss" if args.pss else "tran",
        "period": args.period, "n_periods": n_periods,
        "window": [tw0, tw1],
        "vin": args.vin, "load": args.load,
        "p_in_w": p_in, "p_out_w": p_out,
        "categories_w": totals, "sum_losses_w": sum_losses,
        "residual_w": residual, "residual_pct_of_pin": residual_pct,
        "efficiency": efficiency,
        "per_device": per_device_report,
        "excluded": excluded,
        "gate_drive_modelled": bool(args.gate_drive),
        "switch_low": args.switch_low, "switch_high": args.switch_high,
        # kept for post-hoc red-proof recomputation by tests (drop one device's avg power
        # and recheck the balance) without re-running ngspice:
        "_per_device_avg": {d.name: v for d, v in per_device_avg.items()},
        "_vin_name": vin_dev.name, "_load_name": load_dev.name,
        "_switch_names": [d.name for d in switch_split],
        "_switch_split": {d.name: v for d, v in switch_split.items()},
        "_category_of": {d.name: c for d, c in cat.items()},
    }


def write_outputs(result, outdir):
    pub = {k: v for k, v in result.items() if not k.startswith("_")}
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(pub, fh, indent=2)
    with open(os.path.join(outdir, "report.txt"), "w") as fh:
        fh.write(render_report(result) + "\n")


def render_report(r):
    lines = [
        f"deck: {r['deck']}  mode: {r['mode']}  window: [{r['window'][0]:.6g}, {r['window'][1]:.6g}]s "
        f"({r['n_periods']} period(s) of {r['period']:.6g}s)",
        f"P_in  ({r['vin']}) = {r['p_in_w']:.6f} W",
        f"P_out ({r['load']}) = {r['p_out_w']:.6f} W",
        f"efficiency = {100*r['efficiency']:.4f} %",
        "loss categories (average W over the window):",
    ]
    for c in CATEGORIES:
        note = "" if (c != "gate_drive" or r["gate_drive_modelled"]) else "  (not modelled in this deck)"
        lines.append(f"  {c:20s} {r['categories_w'][c]: .6f}{note}")
    lines.append(f"  {'sum of losses':20s} {r['sum_losses_w']: .6f}")
    if r["excluded"]:
        lines.append(f"ENERGY BALANCE: COULD NOT RUN -- {len(r['excluded'])} unsupported "
                     f"device(s) are NOT included in the sum below (see 'excluded' list); "
                     f"this is not a valid pass/fail result (card #481: an incomplete sum "
                     f"is refused, never silently reported as passing):")
        lines.append(f"  residual (untrusted, incomplete sum) = P_in - (P_out + sum_losses) "
                     f"= {r['residual_w']: .6e} W ({r['residual_pct_of_pin']:+.4f}% of P_in)")
        lines.append("excluded (unsupported) devices:")
        for e in r["excluded"]:
            lines.append(f"  {e['name']} (kind {e['kind']}): {e['note']}")
    else:
        lines.append(f"residual = P_in - (P_out + sum_losses) = {r['residual_w']: .6e} W "
                     f"({r['residual_pct_of_pin']:+.4f}% of P_in)")
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
