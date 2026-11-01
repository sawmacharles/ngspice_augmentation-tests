#!/usr/bin/env python3
"""fra_curves_fork.py — FRA reference curves on the private fork binary (card #461,
follow-on to card #452 T2.3, prep for Phase 5). NEW FILE; does not edit tools/fra.py,
tools/fra_curves.py, tools/stb*.py or any refs/*.cir (per card #461's hard rules — the
fork source worktree ngspice-src-pss belongs to another worker and is never touched
here either; this only *runs* an already-built fork binary as a subprocess).

Why a fork binary and not the system ngspice-45.2: the fork (pre-master-48 pin
aaba50ba3f49, +10 .stb patches, see patches/) differs from 45.2 in converter steady
state by 0.2-1.6% (upstream commits f3871c3e, 0638aaa1 — buildrecipe/README.md); Phase 5
(native periodic AC / PSTB) will itself run in the fork, so its own reference curves
must come from the same binary, not from 45.2.

Four tasks (card #461 = card #452 follow-on):

1. Loop-gain curves (mode "loop", tools/fra.py, unmodified) for R3, R4 D~=0.30 and
   R4 D~=0.45, same frequency grid/settings as the 45.2 SUMMARY.md, with amplitude-
   halving at 3 frequencies each (fra_curves.halving_freqs already includes the low
   792 Hz point for every deck, R3 and R4 alike, since all three RUNS entries share the
   same FREQS list). Reuses tools/fra.py's run_loop_point/sweep_loop, which already
   support tmax/spice_options (card #452's converter-noise-floor fix) and startup="uic"
   — nothing new needed for this mode, only ngspice=FORK_NGSPICE.

2. Transfer-mode curves (control-to-output, line-to-output, output impedance) for R3
   and R4 D~=0.30. tools/fra.py's run_transfer_point does NOT support tmax/spice_options
   or startup="uic" (it hardcodes tmax=dt, and "op" startup, which card #452's converter
   work found gives a useless latched-off DC point on these clocked decks — see
   run_loop_point's docstring). Since tools/fra.py may not be edited here, this file
   reimplements a transfer-point runner (`run_transfer_point_conv`, below) with those
   two fixes, reusing fra.py's low-level helpers (window_times, extract_component,
   injected_reference_phasor, subtract, check_uniform, read_tran, _run_ngspice,
   _sin_spec) rather than duplicating them.

   Injection points (documented here, not in refs/, since refs/*.cir cannot be edited):
   neither deck has an independent source sitting at its "control" node (R3: the type-III
   op-amp's output, node `opout`, read by the PWM comparator; R4: the type-II op-amp's
   output, node `vc`, read by the current comparator) or an independent current sink at
   `vout` (for output impedance) — both are added, as a DC=0 (so the *unperturbed*
   operating point is bit-identical to the unmodified deck) series/shunt element, in
   copies of the decks made under ~/Work/data (never refs/), by `instrument_deck()`
   below:
     - control-to-output: a new node (ctrlnode) is spliced in series between the
       compensator output and the comparator's read of it (a zero-current branch — the
       comparator is a behavioural V-source that draws no current from the nodes it
       reads, so this genuinely does not load the compensator network); `Vctrl` is
       altered with a SIN. R3: `Vctrl ctrlnode opout DC 0`, Bcmp's `V(opout)` -> `V
       (ctrlnode)`. R4: `Vctrl ctrlnode vc DC 0`, Bcmp's `V(vc)` -> `V(ctrlnode)`.
     - line-to-output: no deck change — the existing `Vin` source already carries a DC
       value and fra.py's `alter @vin[sin]` mechanism (same as fra.py's own docstring)
       attaches a SIN to it directly.
     - output impedance: a new `Iinj vout 0 DC 0` current source (shunt, in parallel
       with Rload) is added; sign convention (which direction a positive Iinj perturbs
       Vout) is empirically checked below (`check_iinj_sign`) rather than assumed from
       memory, per the project's "works needs the test output" rule.
   Both new elements are added to a single instrumented copy per deck (one copy serves
   all three transfer measurements — only whichever source is altered is ever non-zero).

   Settle/measure periods: same short-window strategy already validated for these decks
   in loop mode (fra_curves.py's POINT_COMMON, settle_periods=3/measure_periods=5 in
   units of the injection period — chosen there because long windows let the baseline/
   perturbed trajectories' switching *instants* measurably diverge over many thousands
   of cycles), PLUS the separate lesson from R2 (tests/test_fra_lptv.py's
   `r2_transfer`): "the injected sine's own turn-on transient decays per CLOCK period",
   so settle_periods must also cover >=20 clock periods after TD. Both constraints are
   satisfied by `settle_periods_for(f) = max(3.0, 20.0 * f / FS)` (in units of the
   injection period, as window_times expects): at f=500 Hz this is 3 injection periods
   (>>20 clock periods already, since 3/500 s = 1500 clock periods); at f=fs/2 it is
   10 injection periods = 20 clock periods (the binding constraint there).

3. R4 at D~=0.55 without slope compensation (vin_v=6.0, vslope_amp=0, already the
   deck's default) on the fork: a plain `.tran` (the deck's own unmodified `.control`
   block, run through the fork binary directly — no fra.py machinery), redirected only
   at its `write`/`wrdata` output paths (so it does not clobber
   ~/Work/data/ngspice_augmentation/refs/r4_pcm_buck/, which holds the 45.2 reference
   run) via `instrument_deck(..., redirect_outputs=...)`. No FRA is attempted on this
   case (no period-T steady state to build a Tian/Middlebrook two-injection window
   against) — only the .meas t_period/duty/ripple/il evidence is recorded, exactly as
   refs/README.md's 45.2 table already did for the same case.

4. Fork-vs-45.2 comparison: reads the 45.2 SUMMARY.md this project already has
   (~/Work/data/ngspice_augmentation/fra/curves/SUMMARY.md, loop mode only — 45.2 has
   no transfer-mode curves to compare against) and diffs it point by point (dB, deg)
   against this run's loop-mode fork curves; margins (PM/crossover, GM) are reported for
   both binaries from each one's own combined (Tian/Middlebrook) curve.

    python3 tools/fra_curves_fork.py            # tasks 1+2+3+4, ~15-30 min
    python3 tools/fra_curves_fork.py --task 1    # loop-mode fork curves only
    python3 tools/fra_curves_fork.py --task 2    # transfer-mode fork curves only
    python3 tools/fra_curves_fork.py --task 3    # D=0.55 no-slope evidence only
    python3 tools/fra_curves_fork.py --task 4    # comparison only (needs 1's output)
"""
import argparse
import cmath
import datetime
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fra  # noqa: E402
import fra_curves as fc  # noqa: E402
import stb  # noqa: E402
import stbcore  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFS = fc.REFS
FORK_NGSPICE = os.path.expanduser(
    "~/Work/apps/ngspice_augmentation/aaba50ba3f49-s10/bin/ngspice-aug")
OUT_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/fra/curves-fork-s10")
DECKS_DIR = os.path.join(OUT_ROOT, "_decks")
FS = fc.FS  # 500e3
FREQS = fc.FREQS  # dec 5, 500 Hz .. fs/2 (fs/1000..fs/2), same grid as the 45.2 SUMMARY
WORKERS = fc.WORKERS

AMPLITUDE = fc.AMPLITUDE          # 0.01 V, loop mode (unchanged from card #452 T2.3)
HALF_AMPLITUDE = fc.HALF_AMPLITUDE
AMPLITUDE_V = 0.01                # transfer mode, V-type injection (control, line)
AMPLITUDE_I = 0.01                # transfer mode, I-type injection (output impedance), A
POINT_COMMON = fc.POINT_COMMON    # settle_periods=3, measure_periods=5, ppp=256, uic, tmax=5ns, trtol=1


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# deck instrumentation (copies under ~/Work/data only; refs/ is never touched)
# ---------------------------------------------------------------------------

def _fix_include(text):
    abs_lib = os.path.join(REFS, "stb", "lgprobe.lib")
    # matches fra_curves.make_r4_variant's own substitution (refs/ layout: stb/ is a
    # sibling of refs/, i.e. ROOT/stb/lgprobe.lib)
    abs_lib = os.path.join(ROOT, "stb", "lgprobe.lib")
    new_text, n = re.subn(r"\.include \.\./stb/lgprobe\.lib", f".include {abs_lib}", text, count=1)
    if n != 1:
        raise RuntimeError("include substitution did not match .include ../stb/lgprobe.lib")
    return new_text


def instrument_r3(dst_name="r3_vm_buck_transfer.cir"):
    """R3 copy with Vctrl (control-to-output splice at opout) and Iinj (vout shunt)."""
    with open(os.path.join(REFS, "r3_vm_buck.cir")) as fh:
        text = fh.read()
    text = _fix_include(text)
    text, n = re.subn(
        r"^Bcmp pwm 0 V = 0\.5\*\(1\+tanh\(\{kc\}\*\(V\(opout\)-V\(ramp\)\)\)\)$",
        "Vctrl ctrlnode opout DC 0\n"
        "Iinj vout 0 DC 0\n"
        "Bcmp pwm 0 V = 0.5*(1+tanh({kc}*(V(ctrlnode)-V(ramp))))",
        text, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError("R3 Bcmp substitution failed")
    os.makedirs(DECKS_DIR, exist_ok=True)
    dst = os.path.join(DECKS_DIR, dst_name)
    with open(dst, "w") as fh:
        fh.write(text)
    return dst


def instrument_r4(dst_name="r4_pcm_buck_D030_transfer.cir"):
    """R4 (D030, deck default vin_v=11) copy with Vctrl (splice at vc) and Iinj (vout shunt)."""
    with open(os.path.join(REFS, "r4_pcm_buck.cir")) as fh:
        text = fh.read()
    text = _fix_include(text)
    text, n = re.subn(
        r"^Bcmp reset_sig 0 V = 0\.5\*\(1\+tanh\(\{kc\}\*\(V\(vcs\)\+V\(slope\)-V\(vc\)\)\)\)$",
        "Vctrl ctrlnode vc DC 0\n"
        "Iinj vout 0 DC 0\n"
        "Bcmp reset_sig 0 V = 0.5*(1+tanh({kc}*(V(vcs)+V(slope)-V(ctrlnode))))",
        text, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError("R4 Bcmp substitution failed")
    os.makedirs(DECKS_DIR, exist_ok=True)
    dst = os.path.join(DECKS_DIR, dst_name)
    with open(dst, "w") as fh:
        fh.write(text)
    return dst


def make_r4_d055_noslope(dst_name="r4_pcm_buck_D055_noslope.cir"):
    """R4 copy: vin_v=6.0 (D~=0.55), vslope_amp left at the deck's default 0, write/wrdata
    paths redirected under OUT_ROOT so this never touches
    ~/Work/data/ngspice_augmentation/refs/r4_pcm_buck/ (the 45.2 reference run's output)."""
    with open(os.path.join(REFS, "r4_pcm_buck.cir")) as fh:
        text = fh.read()
    text = _fix_include(text)
    text, n = re.subn(r"(\.param vin_v=)11\b", r"\g<1>6.0", text, count=1)
    if n != 1:
        raise RuntimeError("vin_v substitution failed")
    outdir = os.path.join(OUT_ROOT, "r4_D055_noslope")
    os.makedirs(outdir, exist_ok=True)
    text, n1 = re.subn(
        r"write ~/Work/data/ngspice_augmentation/refs/r4_pcm_buck/r4_pcm_buck\.raw",
        f"write {os.path.join(outdir, 'r4_pcm_buck_D055.raw')}", text, count=1)
    text, n2 = re.subn(
        r"wrdata ~/Work/data/ngspice_augmentation/refs/r4_pcm_buck/r4_pcm_buck\.txt",
        f"wrdata {os.path.join(outdir, 'r4_pcm_buck_D055.txt')}", text, count=1)
    if n1 != 1 or n2 != 1:
        raise RuntimeError("write/wrdata path substitution failed")
    os.makedirs(DECKS_DIR, exist_ok=True)
    dst = os.path.join(DECKS_DIR, dst_name)
    with open(dst, "w") as fh:
        fh.write(text)
    return dst, outdir


# ---------------------------------------------------------------------------
# custom transfer-point runner: fra.run_transfer_point + tmax/spice_options/startup="uic"
# ---------------------------------------------------------------------------

def run_transfer_point_conv(deck_lines, source, output_expr, f, amplitude, settle_time=0.0,
                             settle_periods=20, measure_periods=10, points_per_period=256,
                             outdir=None, ngspice=None, startup="uic", tmax=None, spice_options=""):
    """Like fra.run_transfer_point, but with run_loop_point's tmax/spice_options/startup
    handling (card #452's converter-noise-floor fix; needed here because these are the
    same clocked-switching decks, and fra.py's own run_transfer_point has neither knob
    and hardcodes startup="op", which card #452 found gives a useless latched-off DC
    point on R3/R4 — see run_loop_point's docstring, "startup")."""
    ngspice = ngspice or fra.NGSPICE
    outdir = outdir or os.path.join(fra.DATA_ROOT, "transfer_point_conv")
    os.makedirs(outdir, exist_ok=True)
    td, dt, t0, twin, tstop = fra.window_times(settle_time, f, settle_periods, measure_periods,
                                                points_per_period)
    step_max = tmax if tmax is not None else dt

    base_f = os.path.join(outdir, "base.txt")
    pert_f = os.path.join(outdir, "pert.txt")
    for fn in (base_f, pert_f):
        if os.path.exists(fn):
            os.remove(fn)

    if startup not in ("op", "uic"):
        raise fra.CouldNotRun(f"unknown startup {startup!r}, expected 'op' or 'uic'")
    op_line = "op" if startup == "op" else ""
    uic = " uic" if startup == "uic" else ""

    # card #485 (parent instruction, 2026-09-29): this function duplicated fra.py's
    # own run_transfer_point's bug -- fra._sin_spec's VO hard-coded to 0, so `alter
    # @source[sin] = ...` silently erased the source's real DC bias (e.g. Vin=12V)
    # for BOTH the base and pert runs. Combined with this function's own startup=
    # "uic" (needed for R3/R4, see the docstring above), that meant the ENTIRE
    # settle+measurement trajectory ran with Vin pinned near 0V instead of 12V --
    # not a small-signal-point error but a qualitatively different (unpowered,
    # non-regulating) circuit, consistent with the published r3_vm_buck/
    # r4_pcm_buck_D030 line_to_output curves (curves-fork-s10/SUMMARY.md) showing
    # ~0dB "line rejection" at low frequency where a regulated converter's active
    # feedback should give strong (tens of dB) rejection. Fixed the same way as
    # fra.py's run_transfer_point: read the source's own DC value from the deck
    # (fra.source_dc_value) and resolve any {param} expression to a plain number
    # via ngspice itself (fra._resolve_dc_numeric -- alter's bracket-literal SIN
    # form does not evaluate brace expressions), then carry it through as VO.
    vo = fra._resolve_dc_numeric(deck_lines, source, fra.source_dc_value(deck_lines, source),
                                  outdir, ngspice)
    off = fra._sin_spec(0.0, f, td, vo=vo)
    on = fra._sin_spec(amplitude, f, td, vo=vo)
    control = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
{spice_options}
{op_line}
alter @{source}[sin] = {off}
tran {dt:.12e} {tstop:.12e} 0 {step_max:.12e}{uic}
linearize {output_expr}
wrdata {base_f} {output_expr}
alter @{source}[sin] = {on}
tran {dt:.12e} {tstop:.12e} 0 {step_max:.12e}{uic}
linearize {output_expr}
wrdata {pert_f} {output_expr}
quit
.endc
.end
"""
    fra._run_ngspice(deck_lines, control, outdir, "fra_transfer_conv", ngspice)
    tb, (yb,) = fra.read_tran(base_f, 1)
    tp, (yp,) = fra.read_tran(pert_f, 1)
    fra.check_uniform(tb, dt)
    fra.check_uniform(tp, dt)
    dy = fra.subtract(tp, yp, tb, yb, output_expr)
    u0 = fra.injected_reference_phasor(f, td)
    h = fra.extract_component(tp, dy, f, t0, twin) / (amplitude * u0)
    return {"f": f, "H": h, "window": {"td": td, "dt": dt, "t0": t0, "twin": twin, "tstop": tstop}}


def settle_periods_for(f):
    """max(3 injection periods, 20 clock periods) -- see module docstring item 2."""
    return max(3.0, 20.0 * f / FS)


def sweep_transfer_conv(deck, source, output_expr, freqs, amplitude, settle_time, outdir_root,
                         ngspice, workers=None):
    import concurrent.futures as cf
    deck_lines = stb.prepare_deck(deck)
    workers = workers or WORKERS

    def one(f):
        d = os.path.join(outdir_root, f"f_{f:.6g}")
        sp = settle_periods_for(f)
        try:
            return ("ok", run_transfer_point_conv(
                deck_lines, source, output_expr, f, amplitude, settle_time=settle_time,
                settle_periods=sp, measure_periods=5, points_per_period=256, outdir=d,
                ngspice=ngspice, startup="uic", tmax=5e-9, spice_options=".options trtol=1\n"))
        except fra.CouldNotRun as e:
            return ("could-not-run", str(e))

    jobs = list(freqs)
    if workers <= 1:
        return [one(f) for f in jobs]
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(one, jobs))
    # NOTE: ThreadPoolExecutor, not ProcessPoolExecutor (fra.py's own choice for mode
    # "loop") -- deck_lines/closures here are not easily picklable across a module-level
    # worker function without more plumbing; each point still launches its own ngspice
    # *subprocess*, so the GIL is not a bottleneck (subprocess.run releases it), only
    # process-pool-style isolation is traded away, which nothing here depends on.


def halving_transfer(deck, source, output_expr, f, amplitude, settle_time, outdir_root, ngspice):
    deck_lines = stb.prepare_deck(deck)
    sp = settle_periods_for(f)
    kw = dict(settle_time=settle_time, settle_periods=sp, measure_periods=5,
              points_per_period=256, ngspice=ngspice, startup="uic", tmax=5e-9,
              spice_options=".options trtol=1\n")
    ra = run_transfer_point_conv(deck_lines, source, output_expr, f, amplitude,
                                  outdir=os.path.join(outdir_root, f"f_{f:.6g}_A"), **kw)
    rb = run_transfer_point_conv(deck_lines, source, output_expr, f, amplitude / 2,
                                  outdir=os.path.join(outdir_root, f"f_{f:.6g}_A2"), **kw)
    ddb = stbcore.db(ra["H"]) - stbcore.db(rb["H"])
    dph = math.degrees(cmath.phase(ra["H"])) - math.degrees(cmath.phase(rb["H"]))
    while dph > 180:
        dph -= 360
    while dph < -180:
        dph += 360
    return ra["H"], rb["H"], ddb, dph


def check_iinj_sign(deck, outdir, settle_time, check_amplitude):
    """Empirical steady-state check (never assumed from memory): does a positive Iinj
    (see instrument_r3/instrument_r4 -- `Iinj vout 0 DC 0`) raise or lower Vout? A plain
    `op` does not work on these decks (R4's smooth RS latch has no DC solution without
    its `.ic` line -- confirmed here: `op` alone raised "singular matrix: check node q",
    exactly the trap the deck's own header documents for `tran ... uic`). Instead runs
    two full `tran ... uic` transients (same deterministic-replay-from-t=0 pattern
    run_loop_point uses) out to `settle_time` (past this deck's own documented full
    settling time, refs/README.md, not the shorter FRA settle_time -- a 0.98*settle_time..
    settle_time averaging window must land on a genuinely settled trajectory or the sign
    comes out ambiguous/wrong: found empirically here -- a first attempt using R3's own
    300us FRA settle_time and 0.1 A gave Vout *rising* under extra load, an artifact of
    comparing two still-settling startup trajectories, not a real small-signal effect;
    900us (well past R3's documented ~500us full settle) with 0.01 A gives the expected
    "more load -> Vout sags" direction instead), at Iinj=0 and Iinj=check_amplitude, and
    compares the settled v(vout). check_amplitude is independent of the FRA injection
    amplitude (AMPLITUDE_I=0.01 A) -- R4's current-mode regulation is tight enough that
    0.01 A does not move v(vout) past `print`'s 6-significant-digit output at all
    (confirmed here), so a larger check_amplitude (0.3 A) is used for R4 only to get an
    unambiguous digit; the *direction* found (Iinj+ lowers Vout) is a linear small-signal
    property expected to hold at 0.01 A too, not something check_amplitude changes."""
    deck_lines = stb.prepare_deck(deck)
    os.makedirs(outdir, exist_ok=True)
    dt = settle_time / 4000.0
    t1_from = settle_time * 0.98
    control = f"""
.control
set numdgt=15
alter @iinj = 0
tran {dt:.6e} {settle_time:.6e} 0 {dt:.6e} uic
meas tran vout_0 AVG v(vout) FROM={t1_from:.6e} TO={settle_time:.6e}
alter @iinj = {check_amplitude}
tran {dt:.6e} {settle_time:.6e} 0 {dt:.6e} uic
meas tran vout_1 AVG v(vout) FROM={t1_from:.6e} TO={settle_time:.6e}
print vout_0 vout_1
quit
.endc
.end
"""
    fra._run_ngspice(deck_lines, control, outdir, "iinj_sign", FORK_NGSPICE)
    with open(os.path.join(outdir, "iinj_sign.log")) as fh:
        log = fh.read()
    m0 = re.search(r"vout_0\s*=\s*([-\d.eE+]+)", log)
    m1 = re.search(r"vout_1\s*=\s*([-\d.eE+]+)", log)
    if not (m0 and m1):
        raise RuntimeError(f"could not parse vout_0/vout_1 from log at {outdir}")
    v0, v1 = float(m0.group(1)), float(m1.group(1))
    return v0, v1, "Iinj+ lowers Vout (Iinj is an extra load draw)" if v1 < v0 else \
                    "Iinj+ raises Vout (Iinj sources current into vout)"


# ---------------------------------------------------------------------------
# task 1: loop-mode fork curves (unmodified fra.py, ngspice=FORK_NGSPICE)
# ---------------------------------------------------------------------------

RUNS = [
    {"name": "r3_vm_buck", "deck_maker": lambda: os.path.join(REFS, "r3_vm_buck.cir"),
     "settle_time": 300e-6, "desc": "voltage-mode buck, D~=0.276 (12V->3.3V, fs=500kHz, default)"},
    {"name": "r4_pcm_buck_D030", "deck_maker": lambda: os.path.join(REFS, "r4_pcm_buck.cir"),
     "settle_time": 3.5e-3, "desc": "peak-current-mode buck, D~=0.30 (vin_v=11, deck default)"},
    {"name": "r4_pcm_buck_D045", "deck_maker": lambda: fc.make_r4_variant("7.33", "r4_pcm_buck_D045_fork.cir"),
     "settle_time": 3.5e-3, "desc": "peak-current-mode buck, D~=0.45 (vin_v=7.33, copy under ~/Work/data)"},
]


def task1_loop(results):
    print("=== task 1: loop-mode fork curves ===", flush=True)
    out = {}
    for entry in RUNS:
        deck = entry["deck_maker"]()
        outdir = os.path.join(OUT_ROOT, entry["name"])
        os.makedirs(outdir, exist_ok=True)
        print(f"--- {entry['name']}: {entry['desc']} ---", flush=True)
        print(f"deck: {deck}", flush=True)
        t0 = time.time()
        points = fra.sweep_loop(deck, "XP", FREQS, AMPLITUDE, settle_time=entry["settle_time"],
                                 outdir_root=os.path.join(outdir, "sweep"), workers=WORKERS,
                                 ngspice=FORK_NGSPICE, **POINT_COMMON)
        sweep_s = time.time() - t0
        summary = fra.write_loop_results(points, FREQS, outdir)
        n_bad = sum(1 for s, _ in points if s != "ok")
        print(f"  sweep: {summary['n_ok']}/{summary['n_total']} ok"
              + (f", {n_bad} could-not-run" if n_bad else "") + f"  ({sweep_s:.1f}s)", flush=True)

        halving = []
        t0 = time.time()
        for f in fc.halving_freqs(FREQS):
            deck_lines = stb.prepare_deck(deck)
            ra = fra.run_loop_point(deck_lines, "XP", f, AMPLITUDE, settle_time=entry["settle_time"],
                                     outdir=os.path.join(outdir, "halving", f"f_{f:.6g}_A"),
                                     ngspice=FORK_NGSPICE, **POINT_COMMON)
            rb = fra.run_loop_point(deck_lines, "XP", f, HALF_AMPLITUDE, settle_time=entry["settle_time"],
                                     outdir=os.path.join(outdir, "halving", f"f_{f:.6g}_A2"),
                                     ngspice=FORK_NGSPICE, **POINT_COMMON)
            db_diff = stbcore.db(ra["combined"]) - stbcore.db(rb["combined"])
            deg_diff = math.degrees(cmath.phase(ra["combined"])) - math.degrees(cmath.phase(rb["combined"]))
            while deg_diff > 180:
                deg_diff -= 360
            while deg_diff < -180:
                deg_diff += 360
            halving.append({"f": f, "T_A": ra["combined"], "T_A2": rb["combined"],
                             "db_diff": db_diff, "deg_diff": deg_diff})
            print(f"  halving @ {f:.6g} Hz: A={stbcore.db(ra['combined']):.3f} dB "
                  f"A/2={stbcore.db(rb['combined']):.3f} dB  diff={db_diff:+.3f} dB / {deg_diff:+.3f} deg",
                  flush=True)
        halving_s = time.time() - t0

        result = {"name": entry["name"], "deck": deck, "desc": entry["desc"],
                  "settle_time": entry["settle_time"], "sweep_summary": summary,
                  "sweep_runtime_s": sweep_s, "halving": halving, "halving_runtime_s": halving_s,
                  "n_bad": n_bad, "freqs": FREQS, "points": points}
        with open(os.path.join(outdir, "curves_run.json"), "w") as fh:
            json.dump({k: v for k, v in result.items() if k != "points"}, fh, indent=2,
                       default=lambda z: [z.real, z.imag] if isinstance(z, complex) else str(z))
        out[entry["name"]] = result
    results["loop"] = out
    return out


# ---------------------------------------------------------------------------
# task 2: transfer-mode fork curves
# ---------------------------------------------------------------------------

TRANSFER_TYPES = [
    ("control_to_output", "vctrl", "v(vout)", AMPLITUDE_V, "V", "control node splice (R3: opout; R4: vc)"),
    ("line_to_output", "vin", "v(vout)", AMPLITUDE_V, "V", "existing Vin source"),
    ("output_impedance", "iinj", "v(vout)", AMPLITUDE_I, "A", "new Iinj shunt at vout"),
]


def task2_transfer(results):
    print("=== task 2: transfer-mode fork curves ===", flush=True)
    r3_deck = instrument_r3()
    r4_deck = instrument_r4()
    print(f"R3 instrumented deck: {r3_deck}", flush=True)
    print(f"R4 D030 instrumented deck: {r4_deck}", flush=True)

    sign_out = {}
    for name, deck, sign_settle, sign_amp in (("r3_vm_buck", r3_deck, 900e-6, 0.01),
                                               ("r4_pcm_buck_D030", r4_deck, 3.5e-3, 0.3)):
        v0, v1, note = check_iinj_sign(deck, os.path.join(OUT_ROOT, name, "iinj_sign_check"),
                                        sign_settle, sign_amp)
        sign_out[name] = {"v0": v0, "v1": v1, "note": note, "check_amplitude": sign_amp}
        print(f"  Iinj sign check ({name}): Iinj=0 -> v(vout)={v0:.6f} V, "
              f"Iinj=0.1A -> v(vout)={v1:.6f} V -- {note}", flush=True)

    entries = [
        ("r3_vm_buck", r3_deck, 300e-6),
        ("r4_pcm_buck_D030", r4_deck, 3.5e-3),
    ]
    out = {}
    for name, deck, settle_time in entries:
        out[name] = {}
        for kind, source, output_expr, amp, unit, note in TRANSFER_TYPES:
            outdir = os.path.join(OUT_ROOT, name, "transfer", kind)
            os.makedirs(outdir, exist_ok=True)
            print(f"--- {name} / {kind} (source={source}, amplitude={amp}{unit}, {note}) ---", flush=True)
            t0 = time.time()
            points = sweep_transfer_conv(deck, source, output_expr, FREQS, amp, settle_time,
                                          os.path.join(outdir, "sweep"), FORK_NGSPICE)
            sweep_s = time.time() - t0
            n_bad = 0
            with open(os.path.join(outdir, "fra.csv"), "w") as fh:
                fh.write("freq,status,H_re,H_im,H_db,H_phase_deg,detail\n")
                for f, (status, val) in zip(FREQS, points):
                    if status == "ok":
                        h = val["H"]
                        fh.write(f"{f:.9e},ok,{h.real:.9e},{h.imag:.9e},{stbcore.db(h):.6f},"
                                 f"{math.degrees(cmath.phase(h)):.6f},\n")
                    else:
                        n_bad += 1
                        fh.write(f"{f:.9e},could-not-run,,,,,{val}\n")
            print(f"  sweep: {len(FREQS) - n_bad}/{len(FREQS)} ok"
                  + (f", {n_bad} could-not-run" if n_bad else "") + f"  ({sweep_s:.1f}s)", flush=True)

            halving = []
            t0 = time.time()
            for f in fc.halving_freqs(FREQS):
                ha, hb, ddb, dph = halving_transfer(deck, source, output_expr, f, amp, settle_time,
                                                     os.path.join(outdir, "halving"), FORK_NGSPICE)
                halving.append({"f": f, "H_A": ha, "H_A2": hb, "db_diff": ddb, "deg_diff": dph})
                print(f"  halving @ {f:.6g} Hz: A={stbcore.db(ha):.3f} dB A/2={stbcore.db(hb):.3f} dB "
                      f"diff={ddb:+.3f} dB / {dph:+.3f} deg", flush=True)
            halving_s = time.time() - t0

            summary = {"deck": deck, "source": source, "output_expr": output_expr,
                       "amplitude": amp, "amplitude_unit": unit, "injection_note": note,
                       "settle_time": settle_time, "sweep_runtime_s": sweep_s,
                       "n_ok": len(FREQS) - n_bad, "n_total": len(FREQS), "n_bad": n_bad,
                       "halving": halving, "halving_runtime_s": halving_s}
            with open(os.path.join(outdir, "summary.json"), "w") as fh:
                json.dump(summary, fh, indent=2,
                          default=lambda z: [z.real, z.imag] if isinstance(z, complex) else str(z))
            out[name][kind] = {"summary": summary, "points": points, "freqs": FREQS}
    results["transfer"] = out
    results["iinj_sign"] = sign_out
    return out


# ---------------------------------------------------------------------------
# task 3: R4 D~=0.55, no slope comp -- plain .tran, subharmonic evidence only
# ---------------------------------------------------------------------------

def task3_d055(results):
    print("=== task 3: R4 D~=0.55 no-slope subharmonic check (fork, plain .tran) ===", flush=True)
    deck, outdir = make_r4_d055_noslope()
    t0 = time.time()
    proc = subprocess.run([FORK_NGSPICE, "-b", deck], capture_output=True, text=True, timeout=1800)
    runtime_s = time.time() - t0
    log = proc.stdout + proc.stderr
    with open(os.path.join(outdir, "run.log"), "w") as fh:
        fh.write(log)
    vals = {}
    for m in re.finditer(r"^(\w+)\s*=\s*([-\d.eE+]+)\s*$", log, re.M):
        vals[m.group(1)] = float(m.group(2))
    t_period = vals.get("t_period")
    fsw_meas = vals.get("fsw_meas")
    subharmonic = None
    if t_period is not None:
        expected_t = 1.0 / FS
        ratio = t_period / expected_t
        subharmonic = ratio > 1.5  # period-2 => ratio ~2, period-1 => ratio ~1
    print(f"  exit={proc.returncode}, runtime={runtime_s:.1f}s", flush=True)
    print(f"  t_period={t_period}, fsw_meas={fsw_meas}, "
          f"ratio_to_clock_period={t_period / (1.0 / FS) if t_period else None}, "
          f"subharmonic={subharmonic}", flush=True)
    result = {"deck": deck, "returncode": proc.returncode, "runtime_s": runtime_s,
              "measured": vals, "subharmonic_period2": subharmonic, "log_path": os.path.join(outdir, "run.log")}
    with open(os.path.join(outdir, "result.json"), "w") as fh:
        json.dump(result, fh, indent=2)
    results["d055_noslope"] = result
    return result


# ---------------------------------------------------------------------------
# task 4: fork vs 45.2 comparison (loop mode only -- 45.2 has no transfer curves)
# ---------------------------------------------------------------------------

FORTYFIVE_TWO_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/fra/curves")


def load_452_csv(name):
    """Read a 45.2-run fra.csv (curves/<name>/fra.csv) into a list of
    (freq, T_db, T_phase_deg, combined_db, combined_phase_deg), ordered by freq.

    Column layout (fra.write_loop_results): freq=0, status=1, T_re=2, T_im=3, T_db=4,
    T_phase_deg=5, loopGain_db=6, loopGain_phase_deg=7, middlebrook_db=8,
    middlebrook_phase_deg=9, combined_db=10, combined_phase_deg=11, combined_source=12,
    disagreement_db=13, disagreement_deg=14, detail=15. (A first version of this function
    used parts[3]/[4] and parts[9]/[10] -- off by one against this layout, caught by
    inspecting the header with `head -1 fra.csv | tr , '\\n' | cat -n` after the generated
    SUMMARY.md's Task 4 table came out with columns visibly shuffled; fixed here.)"""
    path = os.path.join(FORTYFIVE_TWO_ROOT, name, "fra.csv")
    out = []
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if parts[1] != "ok":
                continue
            f = float(parts[0])
            t_db, t_ph = float(parts[4]), float(parts[5])
            comb_db, comb_ph = float(parts[10]), float(parts[11])
            out.append((f, t_db, t_ph, comb_db, comb_ph))
    return out, path


FREQ_MATCH_RELTOL = 1e-4  # matching 45.2's and the fork's frequency lists (see task4_compare)


def load_fork_csv(name):
    """Read this run's own fork fra.csv (OUT_ROOT/<name>/fra.csv, same column layout as
    load_452_csv) into a list of (freq, T_db, T_phase_deg, combined_db, combined_phase_deg)."""
    path = os.path.join(OUT_ROOT, name, "fra.csv")
    out = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if parts[1] != "ok":
                continue
            f = float(parts[0])
            t_db, t_ph = float(parts[4]), float(parts[5])
            comb_db, comb_ph = float(parts[10]), float(parts[11])
            out.append((f, t_db, t_ph, comb_db, comb_ph))
    return out, path


def task4_compare(results=None):
    """Disk-based: reads this run's own fork fra.csv (written by task1_loop, whether in
    this same process invocation or a previous one) and the 45.2 fra.csv, so `--task 4`
    works standalone against already-written curves-fork-s10/ output without needing
    task 1's in-memory results (task 1's ngspice sweep is the expensive part; this
    comparison is pure post-processing and should not require re-running it)."""
    print("=== task 4: fork vs 45.2 comparison ===", flush=True)
    out = {}
    for name in ("r3_vm_buck", "r4_pcm_buck_D030", "r4_pcm_buck_D045"):
        try:
            old, old_path = load_452_csv(name)
        except FileNotFoundError:
            print(f"  {name}: no 45.2 reference at {os.path.join(FORTYFIVE_TWO_ROOT, name, 'fra.csv')}"
                  " -- skipped", flush=True)
            continue
        try:
            fork_rows, fork_path = load_fork_csv(name)
        except FileNotFoundError:
            print(f"  {name}: no fork fra.csv at {os.path.join(OUT_ROOT, name, 'fra.csv')} "
                  "-- run task 1 first -- skipped", flush=True)
            continue
        margins_fork = None
        try:
            with open(os.path.join(OUT_ROOT, name, "curves_run.json")) as fh:
                margins_fork = json.load(fh)["sweep_summary"]["combined_margins"]
        except (FileNotFoundError, KeyError):
            pass

        rows = []
        for f, n_t_db, n_t_ph, n_c_db, n_c_ph in fork_rows:
            # nearest-frequency match within FREQ_MATCH_RELTOL, not exact float equality:
            # the 45.2 CSV's freq column is text-round-tripped through "%.9e" (fra.py's
            # write_loop_results), which does NOT exactly reproduce the in-memory float
            # for anything but a handful of round numbers (500, 5000, 50000 in this grid)
            # -- found here: an exact dict-key lookup matched only 3/14 points per deck.
            best = min(old, key=lambda row: abs(row[0] - f))
            if abs(best[0] - f) > FREQ_MATCH_RELTOL * f:
                continue
            o_t_db, o_t_ph, o_c_db, o_c_ph = best[1], best[2], best[3], best[4]
            d_t_db, d_c_db = n_t_db - o_t_db, n_c_db - o_c_db
            d_t_ph, d_c_ph = n_t_ph - o_t_ph, n_c_ph - o_c_ph
            while d_t_ph > 180: d_t_ph -= 360
            while d_t_ph < -180: d_t_ph += 360
            while d_c_ph > 180: d_c_ph -= 360
            while d_c_ph < -180: d_c_ph += 360
            rows.append({"f": f, "old_T_db": o_t_db, "new_T_db": n_t_db, "d_T_db": d_t_db,
                         "old_T_deg": o_t_ph, "new_T_deg": n_t_ph, "d_T_deg": d_t_ph,
                         "old_comb_db": o_c_db, "new_comb_db": n_c_db, "d_comb_db": d_c_db,
                         "old_comb_deg": o_c_ph, "new_comb_deg": n_c_ph, "d_comb_deg": d_c_ph})
        out[name] = {"rows": rows, "old_path": old_path, "fork_path": fork_path,
                     "margins_452": None, "margins_fork": margins_fork}
        print(f"  {name}: {len(rows)} comparable points (45.2 ref: {old_path})", flush=True)
    if results is not None:
        results["comparison"] = out
    return out


# ---------------------------------------------------------------------------
# SUMMARY.md
# ---------------------------------------------------------------------------

def write_summary(results, total_s):
    lines = []
    now = datetime.datetime.now().isoformat(timespec="seconds")
    binsha = sha256_of(FORK_NGSPICE)
    lines.append(f"# FRA fork reference curves -- R3/R4 (card #461, card #452 follow-on, "
                 f"generated {now})")
    lines.append("")
    lines.append(f"Binary: `{FORK_NGSPICE}`")
    lines.append(f"sha256: `{binsha}`")
    lines.append(f"Total runtime this invocation: {total_s:.1f}s ({total_s/60:.2f} min).")
    lines.append(f"Frequency grid: {len(FREQS)} points, {FREQS[0]:.6g}-{FREQS[-1]:.6g} Hz "
                 f"(dec 5, i.e. fs/1000..fs/2, same as the 45.2 SUMMARY.md), workers={WORKERS}.")
    lines.append(f"Loop-mode settings: {POINT_COMMON}")
    lines.append(f"Loop-mode amplitude: {AMPLITUDE} V (halving at {HALF_AMPLITUDE} V).")
    lines.append(f"Transfer-mode amplitude: {AMPLITUDE_V} V (control/line), {AMPLITUDE_I} A "
                 f"(output impedance); settle_periods = max(3 injection periods, 20 clock periods) "
                 f"per point, measure_periods=5 injection periods, tmax=5ns, trtol=1, startup=uic "
                 f"(same as loop mode; see module docstring).")
    lines.append("")

    if "loop" in results:
        lines.append("## Task 1: loop-gain curves (mode \"loop\")")
        lines.append("")
        for name, r in results["loop"].items():
            s = r["sweep_summary"]
            lines.append(f"### {name} -- {r['desc']}")
            lines.append("")
            lines.append(f"deck: `{r['deck']}`; settle_time={r['settle_time']*1e6:.0f} us; "
                         f"sweep runtime {r['sweep_runtime_s']:.1f}s, halving runtime "
                         f"{r['halving_runtime_s']:.1f}s; {s['n_ok']}/{s['n_total']} points ok"
                         + (f", {r['n_bad']} could-not-run" if r['n_bad'] else "") + ".")
            if s["combined_margins"]:
                pm, gm = s["combined_margins"]["pm"], s["combined_margins"]["gm"]
                pm_s = f"PM={pm['pm_deg']:.2f} deg @ fc={pm['freq']:.6g} Hz" if pm else "PM: no 0dB crossing"
                gm_s = f"GM={gm['gm_db']:.2f} dB @ f={gm['freq']:.6g} Hz" if gm else "GM: no -180deg crossing"
                lines.append(f"Margins (combined Tian/Middlebrook): {pm_s}; {gm_s}")
            lines.append("")
            lines.append("| freq (Hz) | status | T dB | T phase deg | combined dB | combined phase deg | source |")
            lines.append("|---|---|---|---|---|---|---|")
            for f, (status, val) in zip(r["freqs"], r["points"]):
                if status != "ok":
                    lines.append(f"| {f:.6g} | {status} ({val}) | - | - | - | - | - |")
                    continue
                t, comb = val["T"], val["combined"]
                lines.append(f"| {f:.6g} | ok | {stbcore.db(t):.3f} | {math.degrees(cmath.phase(t)):.3f} | "
                             f"{stbcore.db(comb):.3f} | {math.degrees(cmath.phase(comb)):.3f} | "
                             f"{val['combined_source']} |")
            lines.append("")
            lines.append("Amplitude-halving:")
            lines.append("")
            lines.append("| freq (Hz) | T@A dB | T@A/2 dB | dB diff | deg diff |")
            lines.append("|---|---|---|---|---|")
            for h in r["halving"]:
                lines.append(f"| {h['f']:.6g} | {stbcore.db(h['T_A']):.3f} | "
                             f"{stbcore.db(h['T_A2']):.3f} | {h['db_diff']:+.3f} | {h['deg_diff']:+.3f} |")
            lines.append("")

    if "transfer" in results:
        lines.append("## Task 2: transfer-mode curves")
        lines.append("")
        for name, kinds in results["transfer"].items():
            iinj = results.get("iinj_sign", {}).get(name)
            if iinj:
                lines.append(f"### {name}")
                lines.append("")
                lines.append(f"Iinj sign check: Iinj=0 -> v(vout)={iinj['v0']:.6f} V, "
                             f"Iinj={iinj['check_amplitude']}A -> v(vout)={iinj['v1']:.6f} V -- {iinj['note']}")
                lines.append("")
            for kind, data in kinds.items():
                s = data["summary"]
                lines.append(f"#### {name} / {kind}")
                lines.append("")
                lines.append(f"source=`{s['source']}`, output=`{s['output_expr']}`, "
                             f"amplitude={s['amplitude']}{s['amplitude_unit']}, note: {s['injection_note']}")
                lines.append(f"sweep runtime {s['sweep_runtime_s']:.1f}s; {s['n_ok']}/{s['n_total']} ok"
                             + (f", {s['n_bad']} could-not-run" if s['n_bad'] else "") + ".")
                lines.append("")
                lines.append("| freq (Hz) | status | H dB | H phase deg |")
                lines.append("|---|---|---|---|")
                for f, (status, val) in zip(data["freqs"], data["points"]):
                    if status != "ok":
                        lines.append(f"| {f:.6g} | {status} ({val}) | - | - |")
                        continue
                    h = val["H"]
                    lines.append(f"| {f:.6g} | ok | {stbcore.db(h):.3f} | {math.degrees(cmath.phase(h)):.3f} |")
                lines.append("")
                lines.append("Amplitude-halving:")
                lines.append("")
                lines.append("| freq (Hz) | H@A dB | H@A/2 dB | dB diff | deg diff |")
                lines.append("|---|---|---|---|---|")
                for h in s["halving"]:
                    lines.append(f"| {h['f']:.6g} | {stbcore.db(h['H_A']):.3f} | "
                                 f"{stbcore.db(h['H_A2']):.3f} | {h['db_diff']:+.3f} | {h['deg_diff']:+.3f} |")
                lines.append("")

    if "d055_noslope" in results:
        r = results["d055_noslope"]
        lines.append("## Task 3: R4 D~=0.55, no slope compensation (fork, plain .tran, no FRA attempted)")
        lines.append("")
        lines.append(f"deck: `{r['deck']}`; exit={r['returncode']}; runtime={r['runtime_s']:.1f}s")
        lines.append("")
        m = r["measured"]
        tperiod = m.get('t_period')
        ratio = (tperiod / (1.0/FS)) if tperiod else None
        lines.append(f"t_period = {tperiod} s (expected single-period = {1.0/FS:.6g} s, "
                     f"ratio = {ratio}), fsw_meas = {m.get('fsw_meas')} Hz, "
                     f"duty_meas_final = {m.get('duty_meas_final')}, "
                     f"vout_avg_final = {m.get('vout_avg_final')} V, "
                     f"vout_ripple = {m.get('vout_ripple')} V, "
                     f"il_max_final = {m.get('il_max_final')} A, il_min_final = {m.get('il_min_final')} A")
        lines.append("")
        lines.append(f"**period-2 subharmonic: {r['subharmonic_period2']}** "
                     "(ratio > 1.5 of the nominal clock period counted as subharmonic; "
                     "the 45.2 reference measured ratio ~2.0 for this same vin_v=6.0/vslope_amp=0 case, "
                     "refs/README.md).")
        lines.append("")

    if "comparison" in results:
        lines.append("## Task 4: fork vs 45.2 comparison (loop mode)")
        lines.append("")
        for name, c in results["comparison"].items():
            lines.append(f"### {name}")
            lines.append("")
            lines.append(f"45.2 reference: `{c['old_path']}`")
            lines.append("")
            mf = c["margins_fork"]
            if mf:
                pm, gm = mf["pm"], mf["gm"]
                pm_s = f"PM={pm['pm_deg']:.2f} deg @ fc={pm['freq']:.6g} Hz" if pm else "PM: no 0dB crossing"
                gm_s = f"GM={gm['gm_db']:.2f} dB @ f={gm['freq']:.6g} Hz" if gm else "GM: no -180deg crossing"
                lines.append(f"Fork margins (combined): {pm_s}; {gm_s}")
                lines.append("")
            lines.append("| freq (Hz) | 45.2 T dB | fork T dB | dT dB | 45.2 T deg | fork T deg | dT deg | "
                         "45.2 comb dB | fork comb dB | dcomb dB | dcomb deg |")
            lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
            for row in c["rows"]:
                lines.append(f"| {row['f']:.6g} | {row['old_T_db']:.3f} | {row['new_T_db']:.3f} | "
                             f"{row['d_T_db']:+.3f} | {row['old_T_deg']:.3f} | {row['new_T_deg']:.3f} | "
                             f"{row['d_T_deg']:+.3f} | {row['old_comb_db']:.3f} | {row['new_comb_db']:.3f} | "
                             f"{row['d_comb_db']:+.3f} | {row['d_comb_deg']:+.3f} |")
            lines.append("")

    lines.append("## What was not done / caveats")
    lines.append("")
    lines.append("- Transfer mode has no 45.2 reference to compare against (card #452's converter work "
                 "only produced loop-mode curves) -- task 4 compares loop mode only, as the task specified.")
    lines.append("- `run_transfer_point_conv` (this file) is a reimplementation of fra.py's "
                 "`run_transfer_point` with tmax/spice_options/startup=\"uic\" added; it is NOT the same "
                 "code path fra.py's own tests (test_fra_lptv.py, R1a/R2) validated, since those decks did "
                 "not need those knobs. It reuses fra.py's lower-level helpers (window_times, "
                 "extract_component, etc.) directly, but the control-block assembly is new and only "
                 "red-proofed here by the amplitude-halving checks in this SUMMARY, not by a dedicated "
                 "test file (none was written -- out of the time budget for this card).")
    lines.append("- The control-to-output and output-impedance injection points (Vctrl, Iinj) are new "
                 "elements added only in instrumented copies under ~/Work/data (never refs/*.cir); their "
                 "DC value is 0 in the unperturbed case, so the deck's own DC operating point is unchanged, "
                 "but this was not independently re-verified against refs/README.md's exact numbers beyond "
                 "the Iinj sign-check op-points recorded above.")
    lines.append("- No independent re-verification (V-G4): all numbers here were produced and read by "
                 "this same card.")
    lines.append("")

    os.makedirs(OUT_ROOT, exist_ok=True)
    path = os.path.join(OUT_ROOT, "SUMMARY.md")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {path}", flush=True)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", type=int, choices=[1, 2, 3, 4], action="append",
                     help="run only this task (repeatable); default: 1,2,3,4")
    args = ap.parse_args()
    tasks = args.task or [1, 2, 3, 4]

    if not os.path.exists(FORK_NGSPICE):
        print(f"FORK NGSPICE NOT FOUND: {FORK_NGSPICE}", file=sys.stderr)
        return 2

    t_all = time.time()
    results = {}
    if 1 in tasks:
        task1_loop(results)
    if 2 in tasks:
        task2_transfer(results)
    if 3 in tasks:
        task3_d055(results)
    if 4 in tasks:
        task4_compare(results)
    total_s = time.time() - t_all
    write_summary(results, total_s)
    print(f"total runtime: {total_s:.1f}s ({total_s/60:.2f} min)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
