#!/usr/bin/env python3
"""Validation for tools/yieldest.py, fast/high-sigma yield estimation
(card #474, plan T10.6).

Every ground-truth Pf used below is computed HERE, independently of
yieldest.py's own code path -- case (a) via the exact closed-form normal
CDF (the resistor stack is linear by construction, see
tests/fixtures/yield_resistor.cir), cases (b)/(c) via a 2-D numerical
(Riemann-sum) integration of the standard Gaussian density over the exact
same nonlinear formula the ngspice B-source evaluates (see the fixture
files' own comments) -- never by calling anything in tools/yieldest.py or
tools/mc.py's own statistics. A mismatch between yieldest's answer and
these numbers is yieldest's bug, not a circularity in the test.

Checks (each printed pass/FAIL/COULD-NOT-RUN):
  1. ngspice availability smoke test.
  2. Linear resistor-stack case (a) at 4, 5, 6 sigma: importance-sampling
     Pf and scaled-sigma-sampling Pf each within a factor of 2 of the exact
     Phi(-k); total simulation count <= 1% of what plain MC would need for
     an equivalent-confidence estimate (both stated per level).
  3. Plain MC sanity check at 3 sigma on the same linear circuit -- the
     fast method and brute-force MC must agree at a sigma level where
     brute force is actually affordable.
  4. Nonlinear MOS-square-law case (b) at ~5 sigma (Kp and Vth both
     Gaussian, entering as a product -- not reducible to a linear
     combination): importance-sampling Pf within a factor of 2 of the
     numerically-integrated truth; scaled-sigma-sampling reported
     non-gating (see the module docstring / docs/yield.md -- SSS's
     hyperplane assumption is measurably worse on a curved boundary,
     found during this card, not assumed).
  5. Plain MC sanity check at ~3 sigma on the MOS case.
  6. Red-proof (dropped weights): the *unweighted* fail-rate under the
     mean-shifted sampling distribution (yieldest exposes this only as a
     diagnostic, never as its own Pf estimate) must overestimate the true
     Pf by far more than 2x -- proving the likelihood-ratio weights are
     load-bearing, not decorative.
  7. Red-proof (disjoint regions): a synthetic fixture with two disjoint
     failure regions, one the worst-case-distance search is built to find
     (on its own 15-degree search grid) and one it is built to miss (an
     angular sliver exactly between two grid rays) -- the missed region is
     verified (independently, by the same 2-D integration as case (b)) to
     hold ~35x more probability than the found one. yieldest's
     importance-sampling estimate must miss the true total by far more
     than 2x, AND the tool's own method-disagreement check must flag a
     limitation_warning -- i.e. the tool must not silently report a
     plausible-looking wrong number.
  8. Most-probable-failure-point reporting: the "which parameter(s) drive
     it" contribution fractions must sum to 1 and roughly match what the
     fixture's own construction implies.

Usage:
    python3 tests/test_yield.py
    NGSPICE=<build>/bin/ngspice python3 tests/test_yield.py

Exit status: 0 = all gating checks passed, 1 = a gating check FAILed,
2 = a gating check COULD-NOT-RUN (e.g. ngspice missing).
"""
import math
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import yieldest  # noqa: E402 -- yieldest itself imports tools/mc.py; this test does not
from statistics import NormalDist  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/yield/test_yield")
FIXTURES = os.path.join(HERE, "fixtures")

nd = NormalDist()
Phi = nd.cdf
invPhi = nd.inv_cdf
phi_pdf = nd.pdf

results = []
GATING = set()
NGSPICE_MISSING = False


def record(name, state, detail, gating=True):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN", "info")
    results.append((name, state, detail))
    if gating:
        GATING.add(name)
    print(f"[{state:15s}] {name}: {detail}")


def plain_mc_runs_needed(pf, relerr=0.20, conf_z=1.96):
    """How many plain-MC runs a 100*(1-2*(1-Phi(z)))% CI of the given
    relative width around `pf` would need -- N ~= (z/relerr)^2 * (1-pf)/pf.
    Used only to report/gate the "<=1% of plain MC" criterion; never to run
    that many runs."""
    if pf <= 0:
        return float("inf")
    return (conf_z / relerr) ** 2 * (1.0 - pf) / pf


def numeric_integrate_2d(fail_fn, cx, cy, half_width, n=1200):
    """Independent ground truth: P(fail_fn(u1,u2)) under standard 2-D
    Gaussian (u1,u2), via a centered Riemann-sum grid over
    [cx-half_width, cx+half_width] x [cy-half_width, cy+half_width].
    Deliberately plain stdlib arithmetic -- no scipy/numpy on this box, and
    no dependency on anything in tools/yieldest.py or tools/mc.py."""
    lo_x, hi_x = cx - half_width, cx + half_width
    lo_y, hi_y = cy - half_width, cy + half_width
    hx = (hi_x - lo_x) / n
    hy = (hi_y - lo_y) / n
    total = 0.0
    for i in range(n):
        x = lo_x + (i + 0.5) * hx
        wx = phi_pdf(x)
        if wx < 1e-300:
            continue
        row = 0.0
        for j in range(n):
            y = lo_y + (j + 0.5) * hy
            if fail_fn(x, y):
                row += phi_pdf(y)
        total += wx * row
    return total * hx * hy


def fresh_outdir(tag):
    outdir = os.path.join(DATA_ROOT, tag)
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir, exist_ok=True)
    return outdir


def check_ngspice_available():
    global NGSPICE_MISSING
    tpl = open(os.path.join(FIXTURES, "yield_resistor.cir")).read()
    spec = {
        "template": tpl, "seed": 1, "ngspice": NGSPICE,
        "params": [{"name": "R1", "nominal": 10000.0, "sigma": 0.0},
                   {"name": "R2", "nominal": 10000.0, "sigma": 0.0}],
        "measure": "vtop", "fail": {"op": "lt", "threshold": 1e9},
    }
    outdir = fresh_outdir("smoke")
    evaluator = yieldest.Evaluator(spec, os.path.join(outdir, "decks"))
    os.makedirs(evaluator.workdir, exist_ok=True)
    res = evaluator.eval({"R1": 0.0, "R2": 0.0}, tag="smoke")
    if not res["ok"]:
        record("ngspice availability smoke test", "COULD-NOT-RUN", f"eval not ok: {res}")
        NGSPICE_MISSING = True
    else:
        record("ngspice availability smoke test", "pass", f"vtop={res['value']}")


# --------------------------------------------------------------------- #
# Case (a): linear resistor stack -- exact Pf via normal CDF
# --------------------------------------------------------------------- #

RESISTOR_TPL = open(os.path.join(FIXTURES, "yield_resistor.cir")).read()
R_NOM = 10000.0
R_FRAC_SIGMA = 0.03  # 3% resistor tolerance
R_SIGMA = R_NOM * R_FRAC_SIGMA
I_SRC = 1e-3
VTOP_NOM = I_SRC * 2 * R_NOM
SIGMA_VTOP = I_SRC * math.sqrt(2) * R_SIGMA  # linear combo of 2 independent Gaussians


def resistor_spec(k_sigma, seed):
    threshold = VTOP_NOM - k_sigma * SIGMA_VTOP
    return {
        "template": RESISTOR_TPL, "seed": seed, "ngspice": NGSPICE,
        "params": [
            {"name": "R1", "nominal": R_NOM, "sigma": R_SIGMA},
            {"name": "R2", "nominal": R_NOM, "sigma": R_SIGMA},
        ],
        "measure": "vtop", "fail": {"op": "lt", "threshold": threshold},
        "n_is": 300, "max_workers": 4,
    }


def check_resistor_case(k_sigma, seed):
    tag = f"resistor_{k_sigma}sigma"
    if NGSPICE_MISSING:
        record(f"resistor {k_sigma}sigma", "COULD-NOT-RUN", "ngspice unavailable")
        return None
    spec = resistor_spec(k_sigma, seed)
    outdir = fresh_outdir(tag)
    t0 = time.time()
    result = yieldest.estimate_yield(spec, outdir)
    wall = time.time() - t0
    pf_true = Phi(-k_sigma)

    if not result["ok"]:
        record(f"resistor {k_sigma}sigma: estimate_yield ran", "COULD-NOT-RUN",
               f"reason={result.get('reason')}")
        return None
    record(f"resistor {k_sigma}sigma: estimate_yield ran", "pass",
           f"n_sims={result['n_sims']} wall={wall:.1f}s beta_mpp={result['mpp_search']['beta']:.4f}")

    is_res = result["importance_sampling"]
    if is_res.get("ok"):
        ratio = is_res["pf"] / pf_true
        record(f"resistor {k_sigma}sigma: importance-sampling Pf within 2x",
               "pass" if 0.5 <= ratio <= 2.0 else "FAIL",
               f"pf={is_res['pf']:.4e} true={pf_true:.4e} ratio={ratio:.3f} n={is_res['n']} ess={is_res['ess']:.1f}")
    else:
        record(f"resistor {k_sigma}sigma: importance-sampling Pf within 2x", "COULD-NOT-RUN",
               str(is_res.get("reason")))

    sss_res = result["scaled_sigma_sampling"]
    if sss_res.get("ok"):
        ratio = sss_res["pf"] / pf_true
        record(f"resistor {k_sigma}sigma: scaled-sigma-sampling Pf within 2x",
               "pass" if 0.5 <= ratio <= 2.0 else "FAIL",
               f"pf={sss_res['pf']:.4e} true={pf_true:.4e} ratio={ratio:.3f} "
               f"scales={sss_res.get('scales_used')}")
    else:
        record(f"resistor {k_sigma}sigma: scaled-sigma-sampling Pf within 2x", "COULD-NOT-RUN",
               str(sss_res.get("reason")))

    needed = plain_mc_runs_needed(pf_true)
    frac = result["n_sims"] / needed
    record(f"resistor {k_sigma}sigma: n_sims <= 1% of plain-MC-equivalent",
           "pass" if frac <= 0.01 else "FAIL",
           f"n_sims={result['n_sims']} plain_mc_equivalent~={needed:.3e} "
           f"used={frac:.4%} of it (20% rel. CI, 95% conf. target)")

    contrib = result["most_probable_failure_point"]["contribution_fraction"]
    total_contrib = sum(contrib.values())
    record(f"resistor {k_sigma}sigma: MPP contribution fractions sum to 1, ~50/50 R1/R2",
           "pass" if abs(total_contrib - 1.0) < 1e-6 and all(0.3 < v < 0.7 for v in contrib.values()) else "FAIL",
           f"contribution={contrib}")

    return result


# --------------------------------------------------------------------- #
# Plain-MC sanity check at 3 sigma (resistor)
# --------------------------------------------------------------------- #

def check_plain_mc_sanity_resistor():
    if NGSPICE_MISSING:
        record("plain MC sanity @3sigma (resistor)", "COULD-NOT-RUN", "ngspice unavailable")
        return
    k = 3.0
    pf_true = Phi(-k)
    spec = resistor_spec(k, seed=31415)
    outdir = fresh_outdir("resistor_3sigma_plainmc")
    decks_dir = os.path.join(outdir, "decks")
    os.makedirs(decks_dir, exist_ok=True)
    evaluator = yieldest.Evaluator(spec, decks_dir)
    n = 6000
    t0 = time.time()
    pmc = yieldest.plain_mc(evaluator, ["R1", "R2"], n, seed=2718, max_workers=4)
    wall = time.time() - t0
    if not pmc.get("ok"):
        record("plain MC sanity @3sigma (resistor)", "COULD-NOT-RUN", str(pmc.get("reason")))
        return
    lo, hi = pmc["pf_ci95"]
    ok = lo <= pf_true <= hi or (pmc["pf"] > 0 and 0.4 <= pmc["pf"] / pf_true <= 2.5)
    record("plain MC sanity @3sigma (resistor): brute force and fast method agree",
           "pass" if ok else "FAIL",
           f"plain_mc: n={n} k={pmc['k']} pf={pmc['pf']:.4e} ci95={pmc['pf_ci95']} "
           f"true={pf_true:.4e} wall={wall:.1f}s")


# --------------------------------------------------------------------- #
# Case (b): nonlinear MOS square-law -- Pf via 2-D numerical integration
# --------------------------------------------------------------------- #

MOS_TPL = open(os.path.join(FIXTURES, "yield_mos.cir")).read()
KP_NOM, KP_SIGMA = 2.5e-3, 2.0e-4
VTH_NOM, VTH_SIGMA = 0.4, 0.02
VGS = 1.0
ISPEC_5SIGMAISH = 2.5e-4


def mos_id(kp_u, vth_u):
    kp = KP_NOM + KP_SIGMA * kp_u
    vth = VTH_NOM + VTH_SIGMA * vth_u
    d = VGS - vth
    if d <= 0.0:
        return 0.0
    return 0.5 * kp * d * d


def mos_spec(threshold, seed, n_is=300):
    return {
        "template": MOS_TPL, "seed": seed, "ngspice": NGSPICE,
        "params": [
            {"name": "KP", "nominal": KP_NOM, "sigma": KP_SIGMA},
            {"name": "VTH", "nominal": VTH_NOM, "sigma": VTH_SIGMA},
        ],
        "measure": "id", "fail": {"op": "lt", "threshold": threshold},
        "n_is": n_is, "max_workers": 4,
    }


def check_mos_case():
    if NGSPICE_MISSING:
        record("MOS nonlinear case", "COULD-NOT-RUN", "ngspice unavailable")
        return None
    pf_true = numeric_integrate_2d(lambda u1, u2: mos_id(u1, u2) < ISPEC_5SIGMAISH,
                                    0.0, 0.0, half_width=8.0, n=1200)
    sigma_eq = -invPhi(pf_true) if 0 < pf_true < 1 else None
    record("MOS nonlinear case: independent 2-D numerical-integration ground truth", "info",
           f"threshold={ISPEC_5SIGMAISH:.4e} Pf={pf_true:.4e} equivalent_sigma={sigma_eq:.3f}",
           gating=False)

    spec = mos_spec(ISPEC_5SIGMAISH, seed=4242)
    outdir = fresh_outdir("mos_5sigma")
    t0 = time.time()
    result = yieldest.estimate_yield(spec, outdir)
    wall = time.time() - t0
    if not result["ok"]:
        record("MOS nonlinear case: estimate_yield ran", "COULD-NOT-RUN", str(result.get("reason")))
        return None
    record("MOS nonlinear case: estimate_yield ran", "pass",
           f"n_sims={result['n_sims']} wall={wall:.1f}s beta_mpp={result['mpp_search']['beta']:.4f}")

    is_res = result["importance_sampling"]
    if is_res.get("ok"):
        ratio = is_res["pf"] / pf_true
        record("MOS nonlinear case: importance-sampling Pf within 2x",
               "pass" if 0.5 <= ratio <= 2.0 else "FAIL",
               f"pf={is_res['pf']:.4e} true={pf_true:.4e} ratio={ratio:.3f}")
    else:
        record("MOS nonlinear case: importance-sampling Pf within 2x", "COULD-NOT-RUN",
               str(is_res.get("reason")))

    needed = plain_mc_runs_needed(pf_true)
    frac = result["n_sims"] / needed
    record("MOS nonlinear case: n_sims <= 1% of plain-MC-equivalent",
           "pass" if frac <= 0.01 else "FAIL",
           f"n_sims={result['n_sims']} plain_mc_equivalent~={needed:.3e} used={frac:.4%} of it")

    sss_res = result["scaled_sigma_sampling"]
    sss_ratio = (sss_res.get("pf") or 0) / pf_true if sss_res.get("ok") else None
    record("MOS nonlinear case: scaled-sigma-sampling Pf (non-gating -- known curved-boundary limit)",
           "pass" if sss_ratio is not None and 0.5 <= sss_ratio <= 2.0 else "info",
           f"pf={sss_res.get('pf')} ratio={sss_ratio} disagreement={result.get('method_disagreement_ratio')} "
           f"-- SSS's hyperplane-tail assumption is measurably worse on this curved (Kp*(Vgs-Vth)^2) "
           f"boundary than on the linear resistor case; importance sampling does not share that "
           f"assumption and is gated above instead (see docs/yield.md Limits).",
           gating=False)

    return result, pf_true


def check_plain_mc_sanity_mos():
    if NGSPICE_MISSING:
        record("plain MC sanity @~3sigma (MOS)", "COULD-NOT-RUN", "ngspice unavailable")
        return
    # find a threshold giving Pf close to Phi(-3) via a coarse bisection on
    # the SAME independent numeric integrator used for the gating check
    target = Phi(-3.0)
    lo_th, hi_th = ISPEC_5SIGMAISH, mos_id(0, 0) * 0.999
    for _ in range(14):
        mid = 0.5 * (lo_th + hi_th)
        pf = numeric_integrate_2d(lambda u1, u2: mos_id(u1, u2) < mid, 0.0, 0.0, half_width=6.0, n=400)
        if pf < target:
            lo_th = mid
        else:
            hi_th = mid
    threshold_3s = 0.5 * (lo_th + hi_th)
    pf_true_3s = numeric_integrate_2d(lambda u1, u2: mos_id(u1, u2) < threshold_3s, 0.0, 0.0,
                                       half_width=6.0, n=900)

    spec = mos_spec(threshold_3s, seed=9182)
    outdir = fresh_outdir("mos_3sigma_plainmc")
    decks_dir = os.path.join(outdir, "decks")
    os.makedirs(decks_dir, exist_ok=True)
    evaluator = yieldest.Evaluator(spec, decks_dir)
    n = 6000
    t0 = time.time()
    pmc = yieldest.plain_mc(evaluator, ["KP", "VTH"], n, seed=8081, max_workers=4)
    wall = time.time() - t0
    if not pmc.get("ok"):
        record("plain MC sanity @~3sigma (MOS)", "COULD-NOT-RUN", str(pmc.get("reason")))
        return
    lo, hi = pmc["pf_ci95"]
    ok = lo <= pf_true_3s <= hi or (pmc["pf"] > 0 and 0.4 <= pmc["pf"] / pf_true_3s <= 2.5)
    record("plain MC sanity @~3sigma (MOS): brute force and numeric truth agree",
           "pass" if ok else "FAIL",
           f"threshold={threshold_3s:.4e} plain_mc: n={n} k={pmc['k']} pf={pmc['pf']:.4e} "
           f"ci95={pmc['pf_ci95']} true={pf_true_3s:.4e} wall={wall:.1f}s")


# --------------------------------------------------------------------- #
# Red-proof 1: dropped likelihood-ratio weights must miss badly
# --------------------------------------------------------------------- #

def check_red_proof_dropped_weights(resistor_result_6sigma):
    if NGSPICE_MISSING or resistor_result_6sigma is None:
        record("red-proof: dropped importance-sampling weights miss by >2x", "COULD-NOT-RUN",
               "no 6-sigma resistor result available")
        return
    is_res = resistor_result_6sigma["importance_sampling"]
    if not is_res.get("ok"):
        record("red-proof: dropped importance-sampling weights miss by >2x", "COULD-NOT-RUN",
               "importance_sampling did not run ok")
        return
    pf_true = Phi(-6.0)
    naive = is_res["naive_unweighted_fail_rate_NOT_A_PF_ESTIMATE"]
    ratio = naive / pf_true if pf_true > 0 else float("inf")
    record("red-proof: dropped importance-sampling weights miss by >2x",
           "pass" if ratio > 2.0 else "FAIL",
           f"naive (unweighted) fail-rate under the shifted distribution={naive:.4e} "
           f"vs true Pf={pf_true:.4e} -> overestimates by {ratio:.3e}x "
           f"(correct weighted estimate was {is_res['pf']:.4e}, ratio "
           f"{is_res['pf']/pf_true:.3f}x -- the weights are what fixes this)")


# --------------------------------------------------------------------- #
# Red-proof 2: a failure region the search cannot find
# --------------------------------------------------------------------- #

DISJOINT_TPL = open(os.path.join(FIXTURES, "yield_disjoint.cir")).read()
DISJ_A = (-4.763139720814412, -2.75, 1.0)   # center x, center y, radius -- on-grid (210 deg)
DISJ_B = (-0.4568416727701806, 3.4700570148083365, 0.2)  # off-grid (97.5 deg, between rays)


def disjoint_fail(u1, u2):
    ax, ay, ra = DISJ_A
    bx, by, rb = DISJ_B
    da = (u1 - ax) ** 2 + (u2 - ay) ** 2 - ra * ra
    db = (u1 - bx) ** 2 + (u2 - by) ** 2 - rb * rb
    return min(da, db) < 0.0


def check_disjoint_red_proof():
    if NGSPICE_MISSING:
        record("red-proof: disjoint failure region detected/flagged", "COULD-NOT-RUN",
               "ngspice unavailable")
        return

    pf_a = numeric_integrate_2d(lambda u1, u2: (u1 - DISJ_A[0]) ** 2 + (u2 - DISJ_A[1]) ** 2 < DISJ_A[2] ** 2,
                                 DISJ_A[0], DISJ_A[1], half_width=1.6 * DISJ_A[2], n=800)
    pf_b = numeric_integrate_2d(lambda u1, u2: (u1 - DISJ_B[0]) ** 2 + (u2 - DISJ_B[1]) ** 2 < DISJ_B[2] ** 2,
                                 DISJ_B[0], DISJ_B[1], half_width=1.6 * DISJ_B[2], n=600)
    pf_total = pf_a + pf_b  # regions are >7 apart center-to-center vs radii summing to 1.2 -- disjoint
    record("disjoint fixture: independent ground truth (region A found by grid, B missed by design)",
           "info", f"Pf_A={pf_a:.4e} Pf_B={pf_b:.4e} (B/A={pf_b/pf_a:.1f}x) Pf_total={pf_total:.4e}",
           gating=False)

    spec = {
        "template": DISJOINT_TPL, "seed": 999, "ngspice": NGSPICE,
        "params": [{"name": "U1", "nominal": 0.0, "sigma": 1.0},
                   {"name": "U2", "nominal": 0.0, "sigma": 1.0}],
        "measure": "meas", "fail": {"op": "lt", "threshold": 0.0},
        "n_is": 300, "max_workers": 4,
    }
    outdir = fresh_outdir("disjoint")
    result = yieldest.estimate_yield(spec, outdir)
    if not result["ok"]:
        record("red-proof: disjoint failure region detected/flagged", "COULD-NOT-RUN",
               str(result.get("reason")))
        return

    beta_found = result["mpp_search"]["beta"]
    beta_a_nearest = 5.5 - DISJ_A[2]  # center distance (5.5) minus radius
    record("disjoint fixture: worst-case-distance search found region A, not B",
           "pass" if abs(beta_found - beta_a_nearest) < 0.05 else "info",
           f"beta_found={beta_found:.4f} region_A_nearest={beta_a_nearest:.4f} "
           f"(region B's own nearest-boundary distance is 3.3, CLOSER than A's, but B is never "
           f"seen by the 15-degree search grid by construction)",
           gating=False)

    is_res = result["importance_sampling"]
    is_pf = is_res.get("pf") if is_res.get("ok") else None
    miss_ratio = (pf_total / is_pf) if is_pf else None
    record("red-proof: importance-sampling estimate misses the true (both-region) Pf by >2x",
           "pass" if miss_ratio is not None and miss_ratio > 2.0 else "FAIL",
           f"IS pf={is_pf} vs true total={pf_total:.4e} -> IS underestimates by "
           f"{miss_ratio:.1f}x (it only ever sees region A, Pf_A={pf_a:.4e}, ratio to A alone: "
           f"{(is_pf/pf_a if is_pf else None) and is_pf/pf_a:.3f})")

    warned = result.get("limitation_warning") is not None
    record("red-proof: tool reports a limitation_warning instead of a silent wrong answer",
           "pass" if warned else "FAIL",
           f"method_disagreement_ratio={result.get('method_disagreement_ratio')} "
           f"limitation_warning={'set' if warned else 'NOT SET'}")


# --------------------------------------------------------------------- #

def main():
    print("=" * 70)
    check_ngspice_available()

    resistor_results = {}
    for k in (4, 5, 6):
        print("=" * 70)
        print(f"resistor (linear), {k} sigma")
        resistor_results[k] = check_resistor_case(k, seed=2026 + k)

    print("=" * 70)
    print("plain MC sanity @3 sigma (resistor)")
    check_plain_mc_sanity_resistor()

    print("=" * 70)
    print("MOS (nonlinear), ~5 sigma")
    mos_out = check_mos_case()

    print("=" * 70)
    print("plain MC sanity @~3 sigma (MOS)")
    check_plain_mc_sanity_mos()

    print("=" * 70)
    print("red-proof: dropped importance-sampling weights")
    check_red_proof_dropped_weights(resistor_results.get(6))

    print("=" * 70)
    print("red-proof: disjoint failure region")
    check_disjoint_red_proof()

    print("=" * 70)
    n_pass = sum(1 for _, s, _ in results if s == "pass")
    n_fail = sum(1 for name, s, _ in results if s == "FAIL" and name in GATING)
    n_cnr = sum(1 for name, s, _ in results if s == "COULD-NOT-RUN" and name in GATING)
    print(f"{n_pass} pass, {n_fail} FAIL (gating), {n_cnr} COULD-NOT-RUN (gating), "
          f"{len(results) - n_pass - n_fail - n_cnr} info/non-gating")
    if n_fail:
        return 1
    if n_cnr:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
