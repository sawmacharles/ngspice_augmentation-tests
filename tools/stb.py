#!/usr/bin/env python3
"""stb.py — Cadence-style stb (Tian two-injection loop gain) on an unmodified ngspice.

Card #451 (plan T2.2), 2026-09-27. The deck carries an lgprobe instance
(stb/lgprobe.lib) in series with a wire that breaks every loop. This script adds a
control block that runs op, then two AC sweeps (voltage injection, then current
injection via alter), and applies Tian's formula.

    tools/stb.py DECK --probe XP --sweep "dec 50 1 1e9" [--out DIR]

Subtask 2.2.4: with --mode dm|cm, the deck instead carries a diffprobe instance
(stb/diffprobe.lib) and this drives the correlated differential/common-mode
two-injection probe instead, applying the same Tian formula to branch 1's
readings. Default (no --mode) is the single-ended lgprobe behaviour above,
unchanged.

    tools/stb.py DECK --probe XP --mode dm|cm --sweep "dec 50 1 1e9" [--out DIR]

Prints phase and gain margins; writes stb.csv (T and Spectre-sign loopGain per
frequency) and summary.json in the output directory (default under
~/Work/data/ngspice_augmentation/stb/, never the repository).
Exit status: 0 = ran, 2 = could not run (ngspice error, missing output).
"""
import argparse
import datetime
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stbcore  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation")
DROP_DOT = (".ac", ".tran", ".op", ".dc", ".noise", ".tf", ".pz", ".sens", ".save", ".end", ".four", ".meas", ".measure")


class CouldNotRun(RuntimeError):
    pass


def prepare_deck(path):
    """Deck text with absolute include paths, and without .end, analyses or control blocks."""
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


def probe_nodes(deck_lines, probe):
    """(a, b) nodes of the lgprobe instance named probe."""
    for line in deck_lines:
        f = line.split()
        if f and f[0].lower() == probe.lower():
            if len(f) < 4 or f[3].lower() != "lgprobe":
                raise CouldNotRun(f"{probe} is not an lgprobe instance: {line!r}")
            return f[1].lower(), f[2].lower()
    raise CouldNotRun(f"probe instance {probe} not found in deck")


DIFFPROBE_NAMES = ("diffprobe", "diffprobe_dm", "diffprobe_cm")


def probe_nodes_diff(deck_lines, probe):
    """(a1, b1, a2, b2) nodes of the diffprobe instance named probe.

    Accepts the deck's generic marker "diffprobe" as well as either variant
    name, so this can be called both before and after select_diff_mode().
    """
    for line in deck_lines:
        f = line.split()
        if f and f[0].lower() == probe.lower():
            if len(f) < 6 or f[5].lower() not in DIFFPROBE_NAMES:
                raise CouldNotRun(f"{probe} is not a diffprobe instance: {line!r}")
            return f[1].lower(), f[2].lower(), f[3].lower(), f[4].lower()
    raise CouldNotRun(f"probe instance {probe} not found in deck")


def select_diff_mode(deck_lines, probe, mode):
    """Rewrite the probe instance's subckt name to diffprobe_dm or diffprobe_cm.

    The deck names the generic marker "diffprobe" (or either variant already);
    this is a text-level rewrite of that one field, same style as prepare_deck's
    include-path rewrite. No `alter`/`alterparam` of a live gain is needed.
    """
    if mode not in ("dm", "cm"):
        raise CouldNotRun(f"unknown diffprobe mode {mode!r}")
    variant = f"diffprobe_{mode}"
    out, hit = [], False
    for line in deck_lines:
        f = line.split()
        if f and f[0].lower() == probe.lower() and len(f) >= 6 and f[5].lower() in DIFFPROBE_NAMES:
            f[5] = variant
            out.append(" ".join(f))
            hit = True
        else:
            out.append(line)
    if not hit:
        raise CouldNotRun(f"probe instance {probe} not found in deck")
    return out


def read_wrdata(path, ncomplex):
    """frequency list and ncomplex complex columns from a wr_singlescale/wr_vecnames file."""
    if not os.path.exists(path):
        raise CouldNotRun(f"ngspice wrote no {path}")
    freqs, cols = [], [[] for _ in range(ncomplex)]
    with open(path) as fh:
        next(fh)
        for row in fh:
            v = [float(x) for x in row.split()]
            if len(v) != 1 + 2 * ncomplex:
                raise CouldNotRun(f"unexpected column count in {path}")
            freqs.append(v[0])
            for k in range(ncomplex):
                cols[k].append(complex(v[1 + 2 * k], v[2 + 2 * k]))
    if not freqs:
        raise CouldNotRun(f"{path} is empty")
    return freqs, cols


def run(deck, probe, sweep, outdir):
    """Run the two injections; return dict with freqs, T, loopGain, raw readings, margins."""
    os.makedirs(outdir, exist_ok=True)
    lines = prepare_deck(deck)
    a, _b = probe_nodes(lines, probe)
    p = probe.lower()
    r1, r2 = os.path.join(outdir, "run1.txt"), os.path.join(outdir, "run2.txt")
    for f in (r1, r2):
        if os.path.exists(f):
            os.remove(f)
    control = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
op
ac {sweep}
alter @v.{p}.vi[acmag]=0
alter @i.{p}.ii[acmag]=1
ac {sweep}
wrdata {r1} ac1.v({p}.x) ac1.i(v.{p}.vi) ac1.v({a})
wrdata {r2} ac2.v({p}.x) ac2.i(v.{p}.vi) ac2.i(v.{p}.vb)
quit
.endc
.end
"""
    wrapper = os.path.join(outdir, "stb_wrapper.cir")
    with open(wrapper, "w") as fh:
        fh.write("\n".join(lines) + control)
    try:
        proc = subprocess.run([NGSPICE, "-b", wrapper], capture_output=True, text=True, timeout=3600)
    except (FileNotFoundError, PermissionError) as e:
        raise CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, "ngspice.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")
    f1, (v1, i1, va1) = read_wrdata(r1, 3)
    f2, (v2, i2, ib2) = read_wrdata(r2, 3)
    if f1 != f2:
        raise CouldNotRun("the two AC runs have different frequency grids")
    t = [stbcore.tian(*z) for z in zip(v1, i1, v2, i2)]
    mb = [stbcore.middlebrook(*z) for z in zip(v1, va1, i2, ib2)]
    return {"freqs": f1, "T": t, "loopGain": [-z for z in t], "middlebrook": mb,
            "raw": {"v1": v1, "i1": i1, "v2": v2, "i2": i2, "va1": va1, "ib2": ib2},
            "margins": stbcore.margins(f1, t)}


def run_diff(deck, probe, mode, sweep, outdir):
    """Run the two correlated injections for --mode dm|cm; same Tian formula as run().

    Reads and drives only branch 1 (a1/b1) of the diffprobe instance; branch 2
    (a2/b2) is driven in lockstep (same acmag) at the fixed phase baked into the
    diffprobe_dm/diffprobe_cm variant (stb/diffprobe.lib). See that file for why
    branch-1-only readback gives the DM/CM loop gain of the whole coupled network.
    """
    os.makedirs(outdir, exist_ok=True)
    lines = prepare_deck(deck)
    a1, _b1, _a2, _b2 = probe_nodes_diff(lines, probe)
    lines = select_diff_mode(lines, probe, mode)
    p = probe.lower()
    r1, r2 = os.path.join(outdir, "run1.txt"), os.path.join(outdir, "run2.txt")
    for f in (r1, r2):
        if os.path.exists(f):
            os.remove(f)
    control = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
op
ac {sweep}
alter @v.{p}.vi1[acmag]=0
alter @v.{p}.vi2[acmag]=0
alter @i.{p}.ii1[acmag]=1
alter @i.{p}.ii2[acmag]=1
ac {sweep}
wrdata {r1} ac1.v({p}.x1) ac1.i(v.{p}.vi1) ac1.v({a1})
wrdata {r2} ac2.v({p}.x1) ac2.i(v.{p}.vi1) ac2.i(v.{p}.vb1)
quit
.endc
.end
"""
    wrapper = os.path.join(outdir, f"stb_diff_{mode}_wrapper.cir")
    with open(wrapper, "w") as fh:
        fh.write("\n".join(lines) + control)
    try:
        proc = subprocess.run([NGSPICE, "-b", wrapper], capture_output=True, text=True, timeout=3600)
    except (FileNotFoundError, PermissionError) as e:
        raise CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, "ngspice.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    errs = [ln for ln in (proc.stdout + proc.stderr).splitlines() if re.search(r"\berror\b", ln, re.I)]
    if proc.returncode != 0 or errs:
        raise CouldNotRun(f"ngspice exit {proc.returncode}; {errs[:3]}")
    f1, (v1, i1, va1) = read_wrdata(r1, 3)
    f2, (v2, i2, ib2) = read_wrdata(r2, 3)
    if f1 != f2:
        raise CouldNotRun("the two AC runs have different frequency grids")
    t = [stbcore.tian(*z) for z in zip(v1, i1, v2, i2)]
    mb = [stbcore.middlebrook(*z) for z in zip(v1, va1, i2, ib2)]
    return {"freqs": f1, "T": t, "loopGain": [-z for z in t], "middlebrook": mb,
            "raw": {"v1": v1, "i1": i1, "v2": v2, "i2": i2, "va1": va1, "ib2": ib2},
            "margins": stbcore.margins(f1, t)}


def write_results(res, outdir):
    ph = stbcore.unwrapped_phase_deg(res["T"])
    phl = stbcore.unwrapped_phase_deg(res["loopGain"])
    with open(os.path.join(outdir, "stb.csv"), "w") as fh:
        fh.write("freq,T_re,T_im,T_db,T_phase_deg,loopGain_db,loopGain_phase_deg\n")
        for f, z, p, pl in zip(res["freqs"], res["T"], ph, phl):
            fh.write(f"{f:.9e},{z.real:.9e},{z.imag:.9e},{stbcore.db(z):.6f},{p:.6f},{stbcore.db(z):.6f},{pl:.6f}\n")
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(res["margins"], fh, indent=2)


def describe(m):
    pm, gm = m["pm"], m["gm"]
    if pm:
        kind = "falling" if pm["falling"] else "rising"
        s = [f"phase margin: {pm['pm_deg']:.3f} deg at {pm['freq']:.6g} Hz ({kind} 0 dB crossing)"]
        if not pm["falling"]:
            s.append("note: PM refers to a RISING |T| crossing (|T| increasing through 0 dB), "
                     "not the usual falling-crossing stability margin"
                     + (" -- a value over 180 deg here is not a real margin" if pm["pm_deg"] > 180 else ""))
    else:
        s = ["phase margin: no 0 dB crossing"]
    s.append(f"gain margin: {gm['gm_db']:.3f} dB at {gm['freq']:.6g} Hz" if gm else "gain margin: no -180 deg crossing")
    if len(m["pm_crossings"]) > 1:
        s.append(f"note: {len(m['pm_crossings'])} 0 dB crossings; the lowest margin is reported")
    if len(m["gm_crossings"]) > 1:
        s.append(f"note: {len(m['gm_crossings'])} -180 deg crossings; the lowest margin is reported")
    # card #476: exactly one of these three fires (disjoint arg(T(fmin)) ranges) --
    # the warning (positive feedback WITH loop gain), a note (positive feedback but
    # loop gain < 1 there), or a neutral integrator-like label. See
    # tools/stbcore.margins()'s docstring for the |arg|/|T| thresholds.
    if m.get("low_f_positive_feedback"):
        s.append(f"WARNING: {m['note']}")
    elif m.get("low_f_positive_feedback_low_gain"):
        s.append(f"note: {m['low_f_positive_feedback_low_gain_note']}")
    elif m.get("low_f_integrator"):
        s.append(f"note: {m['low_f_integrator_note']}")
    return "\n".join(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("deck")
    ap.add_argument("--probe", required=True, help="lgprobe (or diffprobe) instance name, e.g. XP")
    ap.add_argument("--sweep", default="dec 50 1 1e9", help='AC sweep, e.g. "dec 50 1 1e9"')
    ap.add_argument("--out", help="output directory (default: under ~/Work/data/ngspice_augmentation/stb)")
    ap.add_argument("--mode", choices=("dm", "cm"), help="differential/common-mode loop gain via a diffprobe "
                     "instance (stb/diffprobe.lib); default: single-ended lgprobe behaviour, unchanged")
    args = ap.parse_args()
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = args.out or os.path.join(DATA_ROOT, "stb", f"{os.path.splitext(os.path.basename(args.deck))[0]}-{stamp}")
    try:
        if args.mode:
            res = run_diff(args.deck, args.probe, args.mode, args.sweep, outdir)
        else:
            res = run(args.deck, args.probe, args.sweep, outdir)
    except (CouldNotRun, subprocess.TimeoutExpired) as e:
        print(f"COULD NOT RUN: {e}", file=sys.stderr)
        return 2
    write_results(res, outdir)
    print(describe(res["margins"]))
    print(f"results: {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
