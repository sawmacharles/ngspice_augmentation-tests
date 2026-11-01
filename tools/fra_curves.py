#!/usr/bin/env python3
"""fra_curves.py — converter loop-gain reference curves (card #452, T2.3 continuation, 2026-09-27).

Runs tools/fra.py mode "loop" (transient Tian/Middlebrook FRA) across refs/r3_vm_buck.cir
(voltage-mode buck, fs=500kHz) and refs/r4_pcm_buck.cir (peak-current-mode buck, fs=500kHz)
from fs/1000 (500 Hz) to fs/2 (250 kHz), ~5 points/decade, at D~=0.30 (R4's default) and
D~=0.45 (a `.param vin_v=7.33` copy, made under ~/Work/data — refs/ is never edited).

**Startup**: both decks' own native `.control` blocks use `tran ... uic` (device default/
`.ic`-line initial conditions), not ngspice's DC operating point — because `op` evaluates
the PWM ramp/slope PULSE sources at their t=0 value (0 V) and finds a real but useless
"converter latched off" DC equilibrium (pwm=0, vout~0, the compensator op-amp railed at
its VCCS's negative limit; confirmed by inspecting ngspice's "Initial Transient Solution"
printout while building this). `run_loop_point(..., startup="uic")` (tools/fra.py) mirrors
that: the three transients bootstrap from zero ICs (or the deck's own `.ic` line, which
r4_pcm_buck.cir carries) instead, and `settle_time` below covers the converter's own
startup transient (refs/README.md: R3 <1% by ~220us, R4 <1% by a few ms) before injection
turns on.

**Measurement window**: R1/R1b's linear-circuit defaults (settle_periods=20,
measure_periods=10) make the "perturbed minus unperturbed" trick break down here — over
enough switching cycles (tens of thousands, at these low injection frequencies with that
many periods) the two deterministic-but-different trajectories' switching *instants*
measurably diverge (found empirically while building this: T at a fixed frequency swung by
tens of dB and the phase by ~180 deg between settle/measure_periods choices of (5,8),
(20,10) and (40,20) — see the card #452 report for the numbers), so the baseline
subtraction stops cancelling the switching ripple cleanly and the residual leaks into the
sideband-0 extraction. Short windows (settle_periods=3, measure_periods=5, ~4-6 switching-
period-independent injection periods) plus a larger-than-R1/R1b injection amplitude (0.1 V,
chosen below) keep this under control, at the cost of visible point-to-point scatter,
worst in R3's ~1-5kHz band (right at its type-III compensator's double zero, fLC/2~=2.5kHz
per the deck header) — documented, not hidden, in SUMMARY.md.

    python3 tools/fra_curves.py
"""
import cmath
import json
import math
import os
import re
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fra  # noqa: E402
import stb  # noqa: E402
import stbcore  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFS = os.path.join(ROOT, "refs")
OUT_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/fra/curves")
DECKS_DIR = os.path.join(OUT_ROOT, "_decks")

WORKERS = max(1, (os.cpu_count() or 2) // 2)
FS = 500e3
FREQS = fra.parse_freqs(f"dec 5 {FS/1000:.0f} {FS/2:.0f}")  # fs/1000 .. fs/2, ~5 pts/decade

AMPLITUDE = 0.01      # V (10 mV), injected in series in the feedback path (goal: <=10mV, card #452 T2.3)
HALF_AMPLITUDE = AMPLITUDE / 2
# card #452 T2.3 converter-noise fix: the old tmax (= dt = 1/(f*points_per_period)) grows
# with 1/f and, at low excitation frequencies, ends up *larger* than the deck's 2us
# switching period -- the solver could step clean over a switching edge, and the residual
# didn't cancel in the baseline subtraction (root cause; see tools/fra.py's run_loop_point
# docstring and the card #452 report). Forcing a small, frequency-independent tmax plus a
# tighter transient-truncation-error tolerance (trtol=1, vs ngspice's default 7) fixed it:
# halving error dropped from several dB/tens of degrees to <0.1dB/1deg at 10mV on both R3
# and R4 D~=0.30, at three frequencies each spanning 500Hz-fs/4 (report, item 2). tmax=5ns
# was picked empirically: 10ns left an occasional marginal (~0.1-0.2dB) miss at 10mV,
# tmax=1-2ns and/or tighter reltol risked "Timestep too small" non-convergence on these
# smooth-tanh switch models (same trap the deck headers already document for ksw/kc=50).
POINT_COMMON = dict(settle_periods=3, measure_periods=5, points_per_period=256, startup="uic",
                     tmax=5e-9, spice_options=".options trtol=1\n")
COMMON = dict(POINT_COMMON, workers=WORKERS)  # sweep_loop also takes workers; run_loop_point does not

R3_DECK = os.path.join(REFS, "r3_vm_buck.cir")
R4_DECK = os.path.join(REFS, "r4_pcm_buck.cir")


def make_r4_variant(vin_v, name):
    """A `.param vin_v=...` copy of r4_pcm_buck.cir under ~/Work/data (refs/ untouched)."""
    with open(R4_DECK) as fh:
        text = fh.read()
    new_text, n = re.subn(r"(\.param vin_v=)11\b", rf"\g<1>{vin_v}", text, count=1)
    if n != 1:
        raise RuntimeError("vin_v substitution did not match r4_pcm_buck.cir's .param line")
    # this copy no longer lives next to refs/../stb/, so make the include absolute
    abs_lib = os.path.join(REFS, "..", "stb", "lgprobe.lib")
    new_text, n = re.subn(r"\.include \.\./stb/lgprobe\.lib", f".include {abs_lib}", new_text, count=1)
    if n != 1:
        raise RuntimeError("include substitution did not match r4_pcm_buck.cir's .include line")
    os.makedirs(DECKS_DIR, exist_ok=True)
    dst = os.path.join(DECKS_DIR, name)
    with open(dst, "w") as fh:
        fh.write(new_text)
    return dst


RUNS = [
    {"name": "r3_vm_buck", "deck": R3_DECK, "settle_time": 300e-6,
     "desc": "voltage-mode buck, D~=0.276 (12V->3.3V, fs=500kHz, default)"},
    {"name": "r4_pcm_buck_D030", "deck": R4_DECK, "settle_time": 3.5e-3,
     "desc": "peak-current-mode buck, D~=0.30 (vin_v=11, deck default)"},
    {"name": "r4_pcm_buck_D045", "deck": None, "settle_time": 3.5e-3,
     "desc": "peak-current-mode buck, D~=0.45 (vin_v=7.33, copy under ~/Work/data)"},
]


def halving_freqs(freqs):
    """3 representative frequencies for the small-signal amplitude-halving check."""
    n = len(freqs)
    return [freqs[1], freqs[n // 2], freqs[-2]]  # avoid the very first/last (fs/2) points


def run_one(entry):
    deck = entry["deck"]
    if deck is None:
        deck = make_r4_variant("7.33", "r4_pcm_buck_D045.cir")
    outdir = os.path.join(OUT_ROOT, entry["name"])
    os.makedirs(outdir, exist_ok=True)
    print(f"=== {entry['name']}: {entry['desc']} ===", flush=True)
    print(f"deck: {deck}", flush=True)

    t0 = time.time()
    points = fra.sweep_loop(deck, "XP", FREQS, AMPLITUDE, settle_time=entry["settle_time"],
                             outdir_root=os.path.join(outdir, "sweep"), **COMMON)
    sweep_s = time.time() - t0
    summary = fra.write_loop_results(points, FREQS, outdir)
    n_bad = sum(1 for s, _ in points if s != "ok")
    print(f"  sweep: {summary['n_ok']}/{summary['n_total']} ok"
          + (f", {n_bad} could-not-run" if n_bad else "") + f"  ({sweep_s:.1f}s)", flush=True)

    # amplitude-halving small-signal check at 2-3 representative frequencies
    halving = []
    t0 = time.time()
    for f in halving_freqs(FREQS):
        deck_lines = stb.prepare_deck(deck)
        ra = fra.run_loop_point(deck_lines, "XP", f, AMPLITUDE, settle_time=entry["settle_time"],
                                 outdir=os.path.join(outdir, "halving", f"f_{f:.6g}_A"), **POINT_COMMON)
        rb = fra.run_loop_point(deck_lines, "XP", f, HALF_AMPLITUDE, settle_time=entry["settle_time"],
                                 outdir=os.path.join(outdir, "halving", f"f_{f:.6g}_A2"), **POINT_COMMON)
        db_diff = stbcore.db(ra["combined"]) - stbcore.db(rb["combined"])
        deg_diff = math.degrees(cmath.phase(ra["combined"])) - math.degrees(cmath.phase(rb["combined"]))
        while deg_diff > 180:
            deg_diff -= 360
        while deg_diff < -180:
            deg_diff += 360
        halving.append({"f": f, "T_A": ra["combined"], "T_A2": rb["combined"],
                         "db_diff": db_diff, "deg_diff": deg_diff})
        print(f"  halving @ {f:.6g} Hz: A={stbcore.db(ra['combined']):.3f} dB "
              f"A/2={stbcore.db(rb['combined']):.3f} dB  diff={db_diff:+.3f} dB / {deg_diff:+.3f} deg", flush=True)
    halving_s = time.time() - t0
    print(f"  halving checks: ({halving_s:.1f}s)", flush=True)

    result = {
        "name": entry["name"], "deck": deck, "desc": entry["desc"],
        "amplitude": AMPLITUDE, "half_amplitude": HALF_AMPLITUDE,
        "settle_time": entry["settle_time"], "common": COMMON,
        "sweep_summary": summary, "sweep_runtime_s": sweep_s,
        "halving": halving, "halving_runtime_s": halving_s,
        "n_bad": n_bad, "freqs": FREQS, "points": points,
    }
    with open(os.path.join(outdir, "curves_run.json"), "w") as fh:
        json.dump(
            {k: v for k, v in result.items() if k != "points"}, fh, indent=2,
            default=lambda z: [z.real, z.imag] if isinstance(z, complex) else str(z))
    return result


def fmt_margins(m, label):
    if not m:
        return f"{label}: no margins (fewer than 2 ok points)"
    pm, gm = m["pm"], m["gm"]
    pm_s = f"PM={pm['pm_deg']:.2f} deg @ fc={pm['freq']:.6g} Hz" if pm else "PM: no 0 dB crossing"
    gm_s = f"GM={gm['gm_db']:.2f} dB @ f={gm['freq']:.6g} Hz" if gm else "GM: no -180 deg crossing"
    return f"{label}: {pm_s}; {gm_s}"


def main():
    t_all = time.time()
    results = [run_one(e) for e in RUNS]
    total_s = time.time() - t_all

    lines = []
    lines.append("# FRA converter reference curves — R3/R4 (card #452, T2.3 continuation, "
                  f"generated {__import__('datetime').datetime.now().isoformat(timespec='seconds')})")
    lines.append("")
    lines.append(f"Total runtime: {total_s:.1f}s ({total_s/60:.2f} min), "
                  f"{len(FREQS)} frequency points/sweep ({FREQS[0]:.6g}-{FREQS[-1]:.6g} Hz, "
                  f"dec 5 i.e. ~5 pts/decade), workers={WORKERS}.")
    lines.append(f"Injection amplitude: {AMPLITUDE} V (halving check at {HALF_AMPLITUDE} V) -- the "
                  "goal's <=10mV ceiling (card #452 T2.3). The original attempt at this card used 0.1V "
                  "because the *old* tmax handling made 1mV/10mV too noisy to use (order-10dB/many-degree "
                  "scatter); with the tmax/trtol fix below, 10mV passes the halving check cleanly (report) "
                  "and no longer needs the extra headroom.")
    lines.append("")
    lines.append("**Timestep-control fix (card #452 T2.3, see the card report for the full investigation):** "
                  "`tmax` used to equal `dt = 1/(f*points_per_period)`, which grows with 1/f and, below "
                  "~13kHz on these 500kHz-switching decks, is *larger* than one switching period -- the "
                  "solver could step clean over a switching edge and the residual leaked into the sideband-0 "
                  "DFT, not cancelling between the baseline and perturbed runs. `tools/fra.py`'s "
                  "`run_loop_point` now takes `tmax` and `spice_options` independently of `dt`; this sweep "
                  "uses `tmax=5ns` (frequency-independent, well under the 2us switching period) and "
                  "`.options trtol=1` (tighter transient-truncation-error tolerance than ngspice's default "
                  "7). Halving error dropped from several dB / tens of degrees to comfortably under "
                  "0.1dB/1deg at 10mV on both decks (see the halving table below and tests/test_fra_converter.py).")
    lines.append("")
    lines.append("## Margins: FRA curve vs averaged-model estimate")
    lines.append("")
    lines.append("| deck | FRA T (Tian) crossover/PM/GM | FRA combined (Tian/Middlebrook) crossover/PM/GM "
                  "| averaged-model estimate (deck header) |")
    lines.append("|---|---|---|---|")
    header_est = {
        "r3_vm_buck": "fc=50.000 kHz, PM~=73.5 deg (estimate, python cmath calc, not measured)",
        "r4_pcm_buck_D030": "fc=20 kHz, PM~=88.1 deg (estimate, python cmath calc, not measured)",
        "r4_pcm_buck_D045": "fc=20 kHz, PM~=88.1 deg (estimate; same averaged model, D not a factor "
                             "in the current-mode approximation used for it)",
    }
    for r in results:
        s = r["sweep_summary"]
        t_m = fmt_margins(s["margins"], "T")
        c_m = fmt_margins(s["combined_margins"], "combined")
        lines.append(f"| {r['name']} | {t_m} | {c_m} | {header_est.get(r['name'], '?')} |")
    lines.append("")

    lines.append("## Amplitude-halving small-signal check (combined estimate)")
    lines.append("")
    lines.append("| deck | freq (Hz) | T @ A dB | T @ A/2 dB | dB diff | deg diff |")
    lines.append("|---|---|---|---|---|---|")
    for r in results:
        for h in r["halving"]:
            lines.append(f"| {r['name']} | {h['f']:.6g} | {stbcore.db(h['T_A']):.3f} | "
                          f"{stbcore.db(h['T_A2']):.3f} | {h['db_diff']:+.3f} | {h['deg_diff']:+.3f} |")
    lines.append("")

    lines.append("## Full curves")
    lines.append("")
    for r in results:
        lines.append(f"### {r['name']} — {r['desc']}")
        lines.append("")
        lines.append(f"deck: `{r['deck']}`; settle_time={r['settle_time']*1e6:.0f} us; "
                      f"sweep runtime {r['sweep_runtime_s']:.1f}s, halving runtime {r['halving_runtime_s']:.1f}s; "
                      f"{r['sweep_summary']['n_ok']}/{r['sweep_summary']['n_total']} points ok"
                      + (f", {r['n_bad']} could-not-run" if r['n_bad'] else "") + ".")
        lines.append("")
        lines.append("| freq (Hz) | status | T dB | T phase deg | Middlebrook dB | combined dB | "
                      "combined phase deg | combined source | disagree dB | disagree deg |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|")
        for f, (status, val) in zip(r["freqs"], r["points"]):
            if status != "ok":
                lines.append(f"| {f:.6g} | {status} ({val}) | - | - | - | - | - | - | - | - |")
                continue
            t, mb, comb = val["T"], val["middlebrook"], val["combined"]
            lines.append(
                f"| {f:.6g} | ok | {stbcore.db(t):.3f} | {math.degrees(cmath.phase(t)):.3f} | "
                f"{stbcore.db(mb):.3f} | {stbcore.db(comb):.3f} | {math.degrees(cmath.phase(comb)):.3f} | "
                f"{val['combined_source']} | {val['disagreement_db']:+.3f} | {val['disagreement_deg']:+.3f} |")
        lines.append("")

    lines.append("## Suspicious / caveats")
    lines.append("")
    lines.append("- **Not `.ac`-validated**: unlike R1/R1b, no independent closed-form or `.ac` reference "
                  "exists for a switching converter's real (non-averaged) loop gain, so these curves have "
                  "no ground truth to check against beyond the averaged-model estimate in each deck's "
                  "header, which is itself an approximation (ignores switching/sampling dynamics, per "
                  "T5.3 pac_sampling's job) — the task explicitly says differences are expected.")
    lines.append("- **R3's ~1-5kHz region** (near the type-III compensator's double zero at fLC/2~=2516Hz, "
                  "deck header) shows visible run-to-run scatter in both magnitude (a few dB) and phase "
                  "(tens of degrees) between amplitude/window choices tried while tuning this — plausibly "
                  "a real feature (the double zero's fast phase transition) partly obscured by residual "
                  "switching-ripple leakage, not cleanly separable with this brute-force technique at "
                  "reasonable runtime. Treat R3's crossover/PM number with more caution than R4's.")
    lines.append("- **Points at/near fs/2 (250 kHz)**: this is exactly the ambiguous region the task "
                  "warned about — the injection frequency approaches the switching Nyquist, sideband-0 "
                  "is no longer cleanly separated from switching-frequency content, and R4's 250kHz point "
                  "showed an odd phase jump during tuning (not necessarily present at the final settings; "
                  "check the table above for this run).")
    lines.append("- **R4 D=0.55 (with/without slope compensation)** was not run — the task asked for D~=0.30 "
                  "and, if time allowed, D~=0.45; it did, but D=0.55 (the subharmonic-instability case, "
                  "refs/README.md) was out of scope for this card and was not attempted.")
    lines.append("- **Short measurement windows** (settle_periods=3, measure_periods=5, vs R1/R1b's 20/10): "
                  "chosen because longer windows made things *worse*, not better, on these switching decks "
                  "(baseline/perturbed trajectory divergence over many thousands of switching cycles — see "
                  "the module docstring). This trades averaging for reduced divergence risk; it is not the "
                  "R1/R1b-validated default, but the amplitude-halving accuracy at this window IS now "
                  "red-proofed by tests/test_fra_converter.py (R3 at 3 frequencies, R4 D~=0.30 at 2), which "
                  "also red-proofs that the *old* tmax handling (before the card #452 T2.3 fix above) fails "
                  "the same check on the same deck/frequency/amplitude.")
    lines.append("- **10ns tmax and reltol=1e-4 together caused ngspice non-convergence** (\"Timestep too "
                  "small ... trouble with node sw\", exit 0) at one of the three R3 frequencies while tuning "
                  "the fix above — ngspice reports this without the word \"error\", so `tools/fra.py`'s "
                  "`_run_ngspice` used to treat it as a successful run with truncated data (found via a "
                  "downstream ZeroDivisionError in stbcore.middlebrook, not a clear failure at the point of "
                  "the bug); `_run_ngspice` now also treats \"aborted\" / \"timestep too small\" / \"singular "
                  "matrix\" as could-not-run. tmax=5ns + trtol=1 (no reltol change) avoided the trap and is "
                  "what this sweep uses.")
    lines.append("- **R4 D=0.55 fix not re-checked**: the tmax/trtol fix above was validated on R3 and R4 "
                  "D~=0.30/D~=0.45 (the decks the task named); it was not re-run on the D=0.55 "
                  "subharmonic-instability case (out of scope here, same as before).")
    lines.append("- **No independent re-verification** (V-G4): all numbers here were produced and read by "
                  "this same card.")
    lines.append("")

    summary_path = os.path.join(OUT_ROOT, "SUMMARY.md")
    os.makedirs(OUT_ROOT, exist_ok=True)
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {summary_path}")
    print(f"total runtime: {total_s:.1f}s ({total_s/60:.2f} min)")


if __name__ == "__main__":
    main()
