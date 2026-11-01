#!/usr/bin/env python3
"""Validation for tools/compopt.py, the loop-compensation optimizer (card #471, plan T10.3).

Every check is pass / FAIL / COULD-NOT-RUN, same convention as this project's other
test_*.py scripts (see tests/test_stb_native.py). Fixtures: tests/fixtures/
compopt_r1b_detuned.cir (a real, poorly-margined starting point on refs/r1b_bilateral.cir's
topology) and tests/fixtures/compopt_r3_averaged_detuned.cir (an averaged linear stand-in
for refs/r3_vm_buck.cir's power stage + PWM gain + type-III compensator, unstable at its
detuned starting values) — see each fixture's own header comment for the numbers behind
its "detuned" starting point and chosen target.

    NGSPICE_AUG=<fork prefix>/bin/ngspice-aug python3 tests/test_compopt.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN

All ngspice runs are sequential (one process at a time — this project's CPU budget while
another worker benchmarks), and every output lands under
~/Work/data/ngspice_augmentation/compopt/ (never the repository).
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import compopt  # noqa: E402
import stbcore  # noqa: E402

DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/compopt/tests")
FIXTURES = os.path.join(HERE, "fixtures")

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else ("skip" if ok == "skip" else "COULD-NOT-RUN"))
    results.append((name, state, detail))


def read_fixture(name):
    with open(os.path.join(FIXTURES, name)) as fh:
        return fh.read()


# ---------------------------------------------------------------- R1b (T10.3)

R1B_PARAMS = [compopt.Param("cc", 1e-13, 3e-11, log=True), compopt.Param("c1", 1e-14, 1e-11, log=True)]
R1B_START = {"cc": 1e-12, "c1": 0.2e-12}
R1B_SWEEP = "dec 30 10 1e9"


def check_r1b_feasible():
    """Detuned R1b (PM=10.3 deg @ 44.5 MHz, fixture header) -> target fc=7.6 MHz,
    PM=65 deg, a point a hand sweep during this card's design found reachable
    (c1~=0.2p, cc~=10p gives PM=66.4 deg @ 7.64 MHz -- see the fixture header).
    Confirmed independently by tools/stb.py on stock ngspice (a different code
    path, and with NGSPICE unset from the fork build, a different binary)."""
    base = read_fixture("compopt_r1b_detuned.cir")
    target = {"fc": 7.6e6, "pm": 65.0}
    outdir = os.path.join(DATA_ROOT, "r1b_feasible")
    t0 = time.time()
    result = compopt.optimize(base, "XP", R1B_PARAMS, R1B_START, target, R1B_SWEEP, outdir,
                               restarts=2, max_iter=80, seed=0)
    elapsed = time.time() - t0
    b = result.best
    record("r1b feasible: optimizer reports feasible=True (never a false success)",
           result.feasible, f"{result.n_evals} evals, {elapsed:.2f}s, status={b.status}")
    if not result.feasible:
        return None
    record("r1b feasible: crossover within 5%", b.rel_fc_err is not None and b.rel_fc_err <= 0.05,
           f"fc={b.fc:.6g} Hz vs target {target['fc']:.6g} Hz, rel err {b.rel_fc_err*100 if b.rel_fc_err else float('nan'):.3f}%")
    record("r1b feasible: phase margin within 2 deg", b.abs_pm_err is not None and b.abs_pm_err <= 2.0,
           f"PM={b.pm:.4f} deg vs target {target['pm']:.3f} deg, abs err {b.abs_pm_err:.4f} deg")
    record("r1b feasible: final values inside the declared bounds",
           all(p.lo <= b.values[p.name] <= p.hi for p in R1B_PARAMS),
           f"values={b.values}")

    # independent confirmation: tools/stb.py on stock ngspice, same deck text
    try:
        vmargins, _res = compopt.verify_with_stb_script(base, "XP", b.values, R1B_SWEEP,
                                                          os.path.join(outdir, "verify"))
    except compopt.CouldNotRun as e:
        record("r1b feasible: independent tools/stb.py confirmation ran", None, str(e))
        return b
    vpm = vmargins["pm"]
    record("r1b feasible: independent tools/stb.py confirmation found a falling 0 dB crossing",
           vpm is not None and vpm["falling"], f"{vpm}")
    if vpm is not None:
        record("r1b feasible: independent PM within 2 deg of target",
               abs(vpm["pm_deg"] - target["pm"]) <= 2.0,
               f"independent PM={vpm['pm_deg']:.4f} deg vs target {target['pm']:.3f} deg")
        record("r1b feasible: independent PM matches native .stb's PM tightly "
               "(script vs native, same values -- this project's documented 1e-9-ish agreement)",
               abs(vpm["pm_deg"] - b.pm) < 1e-3,
               f"independent {vpm['pm_deg']:.6f} deg vs native {b.pm:.6f} deg")
    return b


def check_r1b_infeasible():
    """PM=89 deg at fc=1e8 Hz (above R1b's second pole region -- a hand sweep during
    this card's design pass found only PM in the low-tens-to-negative range or
    ill-defined rising crossings anywhere near that frequency, never close to 89
    deg on a falling crossing) must be reported infeasible, never as success."""
    base = read_fixture("compopt_r1b_detuned.cir")
    target = {"fc": 1e8, "pm": 89.0}
    outdir = os.path.join(DATA_ROOT, "r1b_infeasible")
    t0 = time.time()
    result = compopt.optimize(base, "XP", R1B_PARAMS, R1B_START, target, R1B_SWEEP, outdir,
                               restarts=3, max_iter=80, seed=0)
    elapsed = time.time() - t0
    record("r1b infeasible: optimizer reports feasible=False",
           result.feasible is False,
           f"{result.n_evals} evals, {elapsed:.2f}s, best fc={result.best.fc}, best pm={result.best.pm}")
    record("r1b infeasible: best point still misses tolerance (not silently close)",
           (result.best.rel_fc_err is None or result.best.rel_fc_err > 0.05
            or result.best.abs_pm_err is None or result.best.abs_pm_err > 2.0),
           f"rel_fc_err={result.best.rel_fc_err}, abs_pm_err={result.best.abs_pm_err}")


# -------------------------------------------------------------- R3 (T10.3)

R3_PARAMS = [compopt.Param("r2", 1e3, 1e6, log=True), compopt.Param("c2", 1e-11, 1e-7, log=True)]
R3_START = {"r2": 100e3, "c2": 0.3e-9}
R3_SWEEP = "dec 40 500 1e7"  # fstart=500Hz: see fixture header, avoids the near-integrator
                             # low_f_positive_feedback false positive at fstart=10Hz


def check_r3_feasible():
    """Detuned R3-averaged (PM=-5.4 deg @ 55.8 kHz, unstable -- fixture header) ->
    target fc=34 kHz, PM=38 deg, close to the real deck's own nominal design point
    (r2=27.51k, c2=2.299n gives PM=38.10 @ 33660 Hz on this averaged model, which
    itself matches the real switching deck's ALREADY-MEASURED FRA result from card
    #452, ~/Work/data/ngspice_augmentation/fra/curves/SUMMARY.md, to 0.6%/0.35 deg --
    see the fixture header for the full validation). Confirmed two ways: tools/stb.py
    on the SAME averaged linear circuit (bonus sanity -- this is still a linear stand-in,
    not the switching deck), and tools/fra.py's transient-injection harness on the REAL
    switching deck (refs/r3_vm_buck.cir, read-only, never edited) with the optimized r2/c2
    substituted in and the converter settings tools/fra_curves.py documents for this deck
    family -- the spec's own required "for R3 also compare ... report the difference"."""
    base = read_fixture("compopt_r3_averaged_detuned.cir")
    target = {"fc": 34000.0, "pm": 38.0}
    outdir = os.path.join(DATA_ROOT, "r3_feasible")
    t0 = time.time()
    result = compopt.optimize(base, "XP", R3_PARAMS, R3_START, target, R3_SWEEP, outdir,
                               restarts=2, max_iter=80, seed=0)
    elapsed = time.time() - t0
    b = result.best
    record("r3 feasible: optimizer reports feasible=True (never a false success)",
           result.feasible, f"{result.n_evals} evals, {elapsed:.2f}s, status={b.status}")
    if not result.feasible:
        return None
    record("r3 feasible: crossover within 5%", b.rel_fc_err is not None and b.rel_fc_err <= 0.05,
           f"fc={b.fc:.6g} Hz vs target {target['fc']:.6g} Hz, rel err {b.rel_fc_err*100 if b.rel_fc_err else float('nan'):.3f}%")
    record("r3 feasible: phase margin within 2 deg", b.abs_pm_err is not None and b.abs_pm_err <= 2.0,
           f"PM={b.pm:.4f} deg vs target {target['pm']:.3f} deg, abs err {b.abs_pm_err:.4f} deg")
    record("r3 feasible: final values inside the declared bounds",
           all(p.lo <= b.values[p.name] <= p.hi for p in R3_PARAMS),
           f"values={b.values}")

    # independent confirmation (1): tools/stb.py on the same averaged linear circuit
    try:
        vmargins, _res = compopt.verify_with_stb_script(base, "XP", b.values, R3_SWEEP,
                                                          os.path.join(outdir, "verify_stb"))
        vpm = vmargins["pm"]
        record("r3 feasible: independent tools/stb.py (averaged model) PM within 2 deg of target",
               vpm is not None and abs(vpm["pm_deg"] - target["pm"]) <= 2.0, f"{vpm}")
    except compopt.CouldNotRun as e:
        record("r3 feasible: independent tools/stb.py (averaged model) confirmation ran", None, str(e))

    # independent confirmation (2): FRA on the real switching deck (required by spec)
    fc = b.fc
    freqs = sorted({round(fc * m, -1) for m in (0.5, 0.75, 1.0, 1.5, 2.0)})
    try:
        summary, _points, deck_path = compopt.verify_with_fra(
            os.path.join(ROOT, "refs", "r3_vm_buck.cir"), b.values, "XP", freqs,
            os.path.join(outdir, "verify_fra"), settle_time=300e-6, workers=1)
    except compopt.CouldNotRun as e:
        record("r3 feasible: independent FRA (real switching deck) confirmation ran", None, str(e))
        return b
    cm = summary.get("combined_margins")
    record("r3 feasible: independent FRA confirmation ran and found a crossing",
           cm is not None and cm["pm"] is not None, f"deck={deck_path}, summary={cm}")
    if cm and cm["pm"]:
        fra_fc, fra_pm = cm["pm"]["freq"], cm["pm"]["pm_deg"]
        fc_diff_pct = abs(fra_fc - fc) / fc * 100
        pm_diff_deg = abs(fra_pm - b.pm)
        # NOT required to match .stb tightly -- averaged != switching near fs/2 (spec's own
        # wording); this just checks the two are in the same ballpark (catches a gross
        # error) and reports the real difference, honestly, rather than asserting agreement.
        record("r3 feasible: FRA (real switching deck) crossover/PM reasonably close to the "
               "averaged-model .stb result (difference reported, not required to vanish)",
               fc_diff_pct <= 25.0 and pm_diff_deg <= 15.0,
               f"FRA: fc={fra_fc:.6g} Hz, PM={fra_pm:.3f} deg  |  averaged .stb: fc={fc:.6g} Hz, "
               f"PM={b.pm:.3f} deg  |  diff: {fc_diff_pct:.2f}% fc, {pm_diff_deg:.3f} deg PM")
    return b


def check_r3_infeasible():
    """fc=150 kHz, PM=80 deg: a hand sweep during this card's design pass found
    {r2,c2} alone (with everything else fixed) cannot push this compensator's
    crossover much past ~60 kHz regardless of gain (R3in/C3in's fixed pole and the
    ESR zero dominate above that) -- must be reported infeasible."""
    base = read_fixture("compopt_r3_averaged_detuned.cir")
    target = {"fc": 150000.0, "pm": 80.0}
    outdir = os.path.join(DATA_ROOT, "r3_infeasible")
    t0 = time.time()
    result = compopt.optimize(base, "XP", R3_PARAMS, R3_START, target, R3_SWEEP, outdir,
                               restarts=3, max_iter=80, seed=0)
    elapsed = time.time() - t0
    record("r3 infeasible: optimizer reports feasible=False",
           result.feasible is False,
           f"{result.n_evals} evals, {elapsed:.2f}s, best fc={result.best.fc}, best pm={result.best.pm}")
    record("r3 infeasible: best point still misses tolerance (not silently close)",
           (result.best.rel_fc_err is None or result.best.rel_fc_err > 0.05
            or result.best.abs_pm_err is None or result.best.abs_pm_err > 2.0),
           f"rel_fc_err={result.best.rel_fc_err}, abs_pm_err={result.best.abs_pm_err}")


# ------------------------------------------------------ red-proof (V-G6)

def check_red_proof_corrupted_parser(r1b_best_values):
    """A deliberately corrupted margin "parser" (sign-flipped T, which shifts every
    phase by 180 deg) fed into verify_with_stb_script's margins_fn must make the
    independent confirmation DISAGREE with the target/native result -- proving the
    confirmation step is load-bearing, not a vacuous check that would pass no
    matter what tools/compopt.py handed it."""
    if r1b_best_values is None:
        record("red-proof: corrupted margin parser makes confirmation fail", None,
               "r1b feasible check did not produce values to corrupt-check against")
        return
    base = read_fixture("compopt_r1b_detuned.cir")
    target_pm = 65.0

    def corrupted_margins(freqs, T):
        return stbcore.margins(freqs, [-t for t in T])

    outdir = os.path.join(DATA_ROOT, "red_proof_corrupted_parser")
    try:
        vmargins, _res = compopt.verify_with_stb_script(base, "XP", r1b_best_values, R1B_SWEEP,
                                                          outdir, margins_fn=corrupted_margins)
    except compopt.CouldNotRun as e:
        record("red-proof: corrupted margin parser makes confirmation fail", None, str(e))
        return
    vpm = vmargins["pm"]
    disagrees = (vpm is None) or (abs(vpm["pm_deg"] - target_pm) > 2.0)
    record("red-proof: corrupted (sign-flipped) margin parser makes the independent "
           "confirmation disagree with the target (never passes silently)",
           disagrees, f"corrupted independent PM = {vpm}")

    # and the SAME values with the REAL parser must still agree (proves the failure
    # above is caused by the corruption, not by r1b_best_values itself being bad)
    vmargins_real, _res2 = compopt.verify_with_stb_script(base, "XP", r1b_best_values, R1B_SWEEP,
                                                            os.path.join(outdir, "control"))
    vpm_real = vmargins_real["pm"]
    record("red-proof control: the SAME values with the real (uncorrupted) parser agree",
           vpm_real is not None and abs(vpm_real["pm_deg"] - target_pm) <= 2.0,
           f"real independent PM = {vpm_real}")


def main():
    r1b_best = check_r1b_feasible()
    check_r1b_infeasible()
    r3_best = check_r3_feasible()
    check_r3_infeasible()
    check_red_proof_corrupted_parser(r1b_best.values if r1b_best else None)

    if r1b_best is None:
        pass  # already recorded above
    if r3_best is None:
        pass

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    npass = sum(s == "pass" for _, s, _ in results)
    print(f"\n{npass}/{len(results)} pass ({sum(s == 'skip' for _, s, _ in results)} skipped, not counted)")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
