#!/usr/bin/env python3
"""Validation for refs/r3d_vm_buck_deadtime.cir and refs/r4d_pcm_buck_deadtime.cir
(card #479): realistic (dead-time, finite switch transition, realistic Ron,
low-side body/freewheel diode) counterparts of refs/r3_vm_buck.cir and
refs/r4_pcm_buck.cir, which tools/loss.py (card #470, docs/loss.md) showed
have NO dead-time and run at ~29% efficiency from genuine shoot-through.

Checks, per deck, pass/FAIL/COULD-NOT-RUN:
  - runs to steady state on ngspice-45.2 (/usr/local/bin/ngspice, or $NGSPICE)
  - runs to steady state on the fork build ($NGSPICE_AUG, default
    ~/Work/apps/ngspice_augmentation/aaba50ba3f49-s14/bin/ngspice-aug)
  - Vout within 1% of 3.3V (both binaries)
  - signal-live / V-G7 sanity: the deck's own .meas windows are not
    silently inverted/garbage (t_period near the nominal switching period,
    0 < duty < 1, ripple >= 0) -- the same trap refs/README.md documents
    for R3/R4's TRIG/TARG FALL-vs-RISE desync
  - no shoot-through: peak SIMULTANEOUS high-side/low-side conduction
    current, recomputed independently from the recorded v(sw)/v(dly_hs)/
    v(dly_ls) trace using the decks' own Bhs/Bls controlling expressions
    (the same technique tools/loss.py uses for B-source devices, per
    docs/loss.md) -- not just a gate-signal-crossing proxy, which a first
    trial (card #479) showed was NOT tight enough: gon=70 is large enough
    that a shallow tanh tail still gives amp-scale conductance*V even
    where the two gate signals' own 0.5-crossings are cleanly separated
  - efficiency > 85% and energy balance <= 0.5%, via `python3 tools/loss.py`
    (its documented CLI, tran path, last 5 periods, same as tests/test_loss.py
    uses for R3/R4)

Usage:
    python3 tests/test_refs_d.py
    NGSPICE=<path> NGSPICE_AUG=<path> python3 tests/test_refs_d.py

Exit status: 0 = all gating checks passed, 1 = a gating check FAILed,
2 = a gating check COULD-NOT-RUN (e.g. a binary missing) and nothing FAILed.
"""
import json
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REFS = os.path.join(ROOT, "refs")
LOSS_PY = os.path.join(ROOT, "tools", "loss.py")

NGSPICE = os.environ.get("NGSPICE", "/usr/local/bin/ngspice")
NGSPICE_AUG = os.environ.get(
    "NGSPICE_AUG") or os.environ.get(
    "NGSPICE",  # the portfolio's build under test (2026-09-30: the pinned default below was silently used)
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s14/bin/ngspice-aug"))

DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/refs")
LOSS_OUT_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/loss/test_refs_d")

VOUT_TARGET = 3.3
VOUT_TOL_PCT = 1.0
OVERLAP_BOUND_FRAC = 0.01  # simultaneous conduction current < 1% of load current
EFF_MIN_PCT = 85.0
BALANCE_BOUND_PCT = 0.5

# common switch-model constants shared by both decks (kept in sync with the
# .param lines in refs/r3d_vm_buck_deadtime.cir / refs/r4d_pcm_buck_deadtime.cir
# by hand -- not re-parsed from the deck, see "what was not done")
GON = 70.0
GOFF = 1e-6
KSWD = 150.0

DECKS = {
    "r3d": dict(
        cir="r3d_vm_buck_deadtime.cir",
        outdir="r3d_vm_buck_deadtime",
        vin=12.0, rload=3.3, period=2e-6, tstop=1.5e-3,
        window=(1498e-6, 1500e-6),
    ),
    "r4d": dict(
        cir="r4d_pcm_buck_deadtime.cir",
        outdir="r4d_pcm_buck_deadtime",
        vin=11.0, rload=3.3, period=2e-6, tstop=6e-3,
        window=(5998e-6, 6000e-6),
    ),
}
# wrdata column layout is identical for both decks:
# v(vout) <ctrl> v(dly_hs) v(dly_ls) v(sw) i(L1) i(vin), each preceded by its
# own time column (wrdata repeats t per vector)
COL_VOUT, COL_DLY_HS, COL_DLY_LS, COL_SW, COL_IL1, COL_IVIN = 1, 5, 7, 9, 11, 13

results = []
GATING = set()


def record(name, state, detail, gating=True):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN", "info")
    results.append((name, state, detail))
    if gating:
        GATING.add(name)
    print(f"[{state:15s}] {name}: {detail}")


class CouldNotRun(Exception):
    pass


def run_ngspice(binary, cir, timeout=200):
    if not os.path.exists(binary) and "/" in binary:
        raise CouldNotRun(f"binary not found: {binary}")
    # the refs' absolute wrdata directory may not exist yet on a fresh
    # machine (found on a second workstation, 2026-09-29)
    try:
        with open(os.path.join(REFS, cir)) as fh:
            for line in fh:
                m = re.match(r"\s*(?:wrdata|write)\s+([~/]\S+)", line)
                if m:
                    os.makedirs(os.path.dirname(os.path.expanduser(m.group(1))), exist_ok=True)
    except OSError:
        pass
    t0 = time.time()
    try:
        proc = subprocess.run([binary, "-b", cir], cwd=REFS, capture_output=True,
                               text=True, timeout=timeout)
    except FileNotFoundError as e:
        raise CouldNotRun(f"cannot start {binary}: {e}") from e
    except subprocess.TimeoutExpired as e:
        raise CouldNotRun(f"{binary} timed out after {timeout}s") from e
    wall = time.time() - t0
    out = proc.stdout + proc.stderr
    if proc.returncode != 0:
        raise CouldNotRun(f"{binary} exit {proc.returncode}: {out[-800:]}")
    return out, wall


def parse_meas(out):
    """Parse the deck's own `print` lines: `name = 1.234e+00` (ngspice's
    control-section print format), the last occurrence of each name wins
    (the .control block also echoes .meas results earlier in a slightly
    different format ("name  = value  from= ... to= ..."); the plain
    `print` block afterwards is what this regex is tuned to)."""
    vals = {}
    for m in re.finditer(r"^(\w+)\s*=\s*([-+0-9.eE]+)\s*$", out, re.MULTILINE):
        name, v = m.group(1), m.group(2)
        try:
            vals[name] = float(v)
        except ValueError:
            pass
    return vals


def steady_state_check(tag, deck_name, binary, cfg):
    check = f"{deck_name}/{tag}: steady-state run + Vout/signal-live"
    try:
        out, wall = run_ngspice(binary, cfg["cir"])
    except CouldNotRun as e:
        record(check, "COULD-NOT-RUN", str(e))
        return None
    vals = parse_meas(out)
    needed = ["t_period", "duty_meas_final", "vout_avg_final", "vout_ripple"]
    missing = [n for n in needed if n not in vals]
    if missing:
        record(check, "FAIL", f"missing measurements in output: {missing}")
        return None
    problems = []
    t_period = vals["t_period"]
    if not (0.9 * cfg["period"] <= t_period <= 1.5 * cfg["period"]):
        problems.append(f"t_period={t_period:.4e}s far from nominal {cfg['period']:.4e}s "
                         f"(V-G7: a desynced TRIG/TARG window silently returns garbage)")
    duty = vals["duty_meas_final"]
    if not (0.0 < duty < 1.0):
        problems.append(f"duty_meas_final={duty} not in (0,1) -- inverted/garbage window")
    ripple = vals["vout_ripple"]
    if ripple < 0:
        problems.append(f"vout_ripple={ripple} < 0 -- inverted MAX/MIN window")
    vout = vals["vout_avg_final"]
    err_pct = abs(vout - VOUT_TARGET) / VOUT_TARGET * 100.0
    if err_pct > VOUT_TOL_PCT:
        problems.append(f"Vout={vout:.6f}V is {err_pct:.3f}% off {VOUT_TARGET}V (bound {VOUT_TOL_PCT}%)")
    if problems:
        record(check, "FAIL", "; ".join(problems))
        return None
    record(check, "pass",
           f"Vout={vout:.6f}V ({err_pct:.4f}% off target), fsw={1/t_period:.1f}Hz, "
           f"duty={duty:.4f}, ripple={ripple*1e3:.4f}mV, wall={wall:.1f}s")
    return vals


def read_wrdata(path):
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    return rows


def switch_cond(dly):
    g = 0.5 * (1.0 + math.tanh(KSWD * (dly - 0.5)))
    return GOFF + (GON - GOFF) * g


def shoot_through_check(deck_name, cfg, vout_final):
    check = f"{deck_name}: no-shoot-through (simultaneous hs/ls current)"
    trace_path = os.path.join(DATA_ROOT, cfg["outdir"], cfg["outdir"] + ".txt")
    if not os.path.exists(trace_path):
        record(check, "COULD-NOT-RUN", f"no trace at {trace_path} (steady-state run must "
               "have written it -- run that check first)")
        return
    rows = read_wrdata(trace_path)
    if not rows:
        record(check, "COULD-NOT-RUN", f"{trace_path} is empty")
        return
    lo, hi = cfg["window"]
    vin = cfg["vin"]
    iload = vout_final / cfg["rload"] if vout_final else cfg["vin"] * 0.3 / cfg["rload"]
    max_overlap = 0.0
    t_max = None
    n_in_window = 0
    for r in rows:
        t = r[0]
        if not (lo <= t <= hi):
            continue
        n_in_window += 1
        sw = r[COL_SW]
        dly_hs = r[COL_DLY_HS]
        dly_ls = r[COL_DLY_LS]
        i_hs = (vin - sw) * switch_cond(dly_hs)
        i_ls = sw * switch_cond(dly_ls)
        overlap = min(abs(i_hs), abs(i_ls))
        if overlap > max_overlap:
            max_overlap = overlap
            t_max = t
    if n_in_window == 0:
        record(check, "COULD-NOT-RUN", f"no trace samples in window [{lo},{hi}]")
        return
    bound = OVERLAP_BOUND_FRAC * iload
    detail = (f"peak simultaneous |I_hs|,|I_ls| overlap = {max_overlap:.6e}A at t={t_max}, "
              f"bound = {OVERLAP_BOUND_FRAC*100:.0f}% of load current ({iload:.4f}A) = {bound:.6e}A, "
              f"n_samples={n_in_window}")
    if max_overlap <= bound:
        record(check, "pass", detail)
    else:
        record(check, "FAIL", detail)


def efficiency_check(deck_name, cfg):
    check = f"{deck_name}: efficiency > {EFF_MIN_PCT:.0f}% and energy balance <= {BALANCE_BOUND_PCT}%"
    deck_path = os.path.join(REFS, cfg["cir"])
    outdir = os.path.join(LOSS_OUT_ROOT, deck_name)
    cmd = [sys.executable, LOSS_PY, deck_path,
           "--vin", "Vin", "--load", "Rload", "--dcr", "Rdcr", "--esr", "Resr",
           "--switch", "Bhs", "--switch", "Bls",
           "--period", str(cfg["period"]), "--tstop", str(cfg["tstop"]), "--periods", "5",
           "--out", outdir]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, cwd=REFS, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        record(check, "COULD-NOT-RUN", "tools/loss.py timed out after 300s")
        return
    wall = time.time() - t0
    if proc.returncode != 0:
        record(check, "COULD-NOT-RUN", f"tools/loss.py exit {proc.returncode}: "
               f"{(proc.stdout+proc.stderr)[-500:]}")
        return
    summ_path = os.path.join(outdir, "summary.json")
    if not os.path.exists(summ_path):
        record(check, "COULD-NOT-RUN", "tools/loss.py exited 0 but wrote no summary.json")
        return
    with open(summ_path) as f:
        summary = json.load(f)
    eff_pct = summary["efficiency"] * 100.0
    residual_pct = abs(summary["residual_pct_of_pin"])
    detail = (f"efficiency={eff_pct:.4f}% (bound >{EFF_MIN_PCT}%), "
              f"|residual|={residual_pct:.4f}% of Pin (bound <={BALANCE_BOUND_PCT}%), "
              f"P_in={summary['p_in_w']:.4f}W P_out={summary['p_out_w']:.4f}W, wall={wall:.1f}s")
    if eff_pct > EFF_MIN_PCT and residual_pct <= BALANCE_BOUND_PCT:
        record(check, "pass", detail)
    else:
        record(check, "FAIL", detail)


def main():
    for deck_name, cfg in DECKS.items():
        cir_path = os.path.join(REFS, cfg["cir"])
        if not os.path.exists(cir_path):
            record(f"{deck_name}: deck exists", "COULD-NOT-RUN", f"missing {cir_path}")
            continue

        # fork build first (so the ngspice-45.2 run below leaves the fresher
        # wrdata trace on disk for the shoot-through check, which is defined
        # against the required/deterministic ngspice-45.2 path)
        steady_state_check("fork", deck_name, NGSPICE_AUG, cfg)
        vals = steady_state_check("ngspice-45.2", deck_name, NGSPICE, cfg)

        vout_final = vals["vout_avg_final"] if vals else None
        shoot_through_check(deck_name, cfg, vout_final)
        efficiency_check(deck_name, cfg)

    print()
    print("=" * 70)
    n_fail = sum(1 for n, s, d in results if s == "FAIL" and n in GATING)
    n_cnr = sum(1 for n, s, d in results if s == "COULD-NOT-RUN" and n in GATING)
    n_pass = sum(1 for n, s, d in results if s == "pass")
    print(f"{n_pass} pass, {n_fail} FAIL, {n_cnr} COULD-NOT-RUN (of {len(GATING)} gating checks)")
    if n_fail:
        return 1
    if n_cnr:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
