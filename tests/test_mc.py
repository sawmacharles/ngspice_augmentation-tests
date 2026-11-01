#!/usr/bin/env python3
"""Validation for tools/mc.py, the Monte Carlo global/local(mismatch) runner
(card #473, plan T10.5).

A matched resistor pair (spec's own "matched MOS or resistor pair" option),
read back per run via `@r1[resistance]`/`@r2[resistance]` (a device-value
readback, per the spec's "e.g. via .meas of the device values or currents").
For all statistical checks N ~= 500, matching the spec's stated N.

Checks (each printed pass/FAIL/COULD-NOT-RUN):
  1. global-only: sampled global parameter moves both devices together
     (correlation == 1, exactly, since with local sigma=0 the two devices'
     values are literally the same formula).
  2. local-only: independent (correlation ~= 0 within its 99% CI), each
     measured sigma matches the configured sigma within its 99% CI.
  3. both: correct mixture -- measured correlation matches the analytic
     sigma_g^2/(sigma_g^2+sigma_l^2) prediction within its 99% CI.
  4. Pelgrom scaling: sigma measured at two different W*L values scales
     as 1/sqrt(W*L) (ratio check, tolerance-based -- see note at the check).
  5. reproducibility: same seed -> identical per-run sampled values AND
     measured results across two independent campaign runs; replaying a
     single run (mc.replay_run) alone reproduces that run's own result,
     using nothing but the spec + seed + run index (no other run's state).
  6. red-proof: a deliberately-shared "local" parameter (two device tokens
     both bound to the SAME param name, i.e. the bug this tool exists to
     avoid) must FAIL the independence check -- proving the checker
     actually distinguishes correct local (independent) from broken local
     (accidentally shared), not just that it always reports "pass".

Usage:
    python3 tests/test_mc.py
    NGSPICE=<build>/bin/ngspice python3 tests/test_mc.py

Exit status: 0 = all gating checks passed, 1 = a gating check FAILed,
2 = a gating check COULD-NOT-RUN (e.g. ngspice missing).
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import mc  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/mc/test_mc")

results = []
GATING = set()
NGSPICE_MISSING = False


def record(name, state, detail, gating=True):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN", "info")
    results.append((name, state, detail))
    if gating:
        GATING.add(name)
    print(f"[{state:15s}] {name}: {detail}")


PAIR_TEMPLATE = """\
* matched resistor pair -- global (shared @@GA@@) + local (@@LA@@/@@LB@@)
R1 1 0 '1000*(1+@@GA@@+@@LA@@)'
R2 2 0 '1000*(1+@@GA@@+@@LB@@)'
V1 1 0 DC 1
V2 2 0 DC 1
.control
op
let r1v = @r1[resistance]
let r2v = @r2[resistance]
set r1s = $&r1v
set r2s = $&r2v
echo MCMEAS r1=$r1s
echo MCMEAS r2=$r2s
.endc
.end
"""

# same as PAIR_TEMPLATE but R2 reuses @@LA@@ instead of its own @@LB@@ --
# the deliberate "shared local" bug for the red-proof.
PAIR_TEMPLATE_SHARED_BUG = PAIR_TEMPLATE.replace("@@LB@@", "@@LA@@")

SINGLE_TEMPLATE = """\
* single resistor, local/mismatch param only -- for the Pelgrom check
R1 1 0 '1000*(1+@@LOC@@)'
V1 1 0 DC 1
.control
op
let r1v = @r1[resistance]
set r1s = $&r1v
echo MCMEAS r1=$r1s
.endc
.end
"""

N = 500
SIGMA_G = 0.02
SIGMA_L = 0.01


def base_params(sigma_g, sigma_l_a, sigma_l_b):
    return [
        {"name": "GA", "scope": "global", "dist": "gauss", "nominal": 0.0, "sigma": sigma_g},
        {"name": "LA", "scope": "local", "dist": "gauss", "nominal": 0.0, "sigma": sigma_l_a},
        {"name": "LB", "scope": "local", "dist": "gauss", "nominal": 0.0, "sigma": sigma_l_b},
    ]


def check_ngspice_available():
    global NGSPICE_MISSING
    outdir = tempfile.mkdtemp(prefix="mc_smoke_")
    try:
        spec = {
            "template": SINGLE_TEMPLATE, "seed": 1, "n": 2, "ngspice": NGSPICE,
            "params": [{"name": "LOC", "scope": "local", "dist": "gauss", "nominal": 0.0, "sigma": 0.0}],
            "measures": ["r1"],
        }
        manifest = mc.run_campaign(spec, outdir)
        if manifest["n_ok"] == 0:
            errs = {r.get("error") for r in manifest["runs"]}
            record("ngspice availability smoke test", "COULD-NOT-RUN", f"0/{manifest['n']} runs ok: {errs}")
            NGSPICE_MISSING = True
        else:
            record("ngspice availability smoke test", "pass", f"{manifest['n_ok']}/{manifest['n']} runs ok")
    finally:
        shutil.rmtree(outdir, ignore_errors=True)


def run_campaign(tag, template, params, measures, checks, n=N, seed=1):
    outdir = os.path.join(DATA_ROOT, tag)
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir, exist_ok=True)
    spec = {
        "template": template, "seed": seed, "n": n, "ngspice": NGSPICE,
        "params": params, "measures": measures, "checks": checks, "max_workers": 4,
    }
    manifest = mc.run_campaign(spec, outdir)
    return spec, manifest


def check_global_only():
    if NGSPICE_MISSING:
        record("global-only: r1/r2 correlation == 1", "COULD-NOT-RUN", "ngspice unavailable")
        return
    params = base_params(SIGMA_G, 0.0, 0.0)
    checks = [
        {"type": "correlation", "a": "r1", "b": "r2", "expect": 1.0},
        {"type": "mean", "measure": "r1", "expect": 1000.0},
        {"type": "sigma", "measure": "r1", "expect": 1000.0 * SIGMA_G},
    ]
    spec, m = run_campaign("global_only", PAIR_TEMPLATE, params, ["r1", "r2"], checks)
    if m["n_ok"] < n_needed(m):
        record("global-only campaign completeness", "COULD-NOT-RUN", f"only {m['n_ok']}/{m['n']} runs ok")
        return
    corr = m["stats"]["correlations"].get("r1:r2", {})
    for c in m["checks"]:
        label = f"global-only: {c['check']['type']} {c['check'].get('a', c['check'].get('measure'))}"
        record(label, c["state"], c["detail"])
    record("global-only: correlation exactly 1 (identical formula, not just correlated)",
           "pass" if corr.get("r") == 1.0 else "FAIL", f"r={corr.get('r')}")


def n_needed(manifest):
    return int(0.95 * manifest["n"])


def check_local_only():
    if NGSPICE_MISSING:
        record("local-only: r1/r2 correlation ~= 0", "COULD-NOT-RUN", "ngspice unavailable")
        return
    params = base_params(0.0, SIGMA_L, SIGMA_L)
    checks = [
        {"type": "correlation", "a": "r1", "b": "r2", "expect": 0.0},
        {"type": "sigma", "measure": "r1", "expect": 1000.0 * SIGMA_L},
        {"type": "sigma", "measure": "r2", "expect": 1000.0 * SIGMA_L},
    ]
    spec, m = run_campaign("local_only", PAIR_TEMPLATE, params, ["r1", "r2"], checks)
    if m["n_ok"] < n_needed(m):
        record("local-only campaign completeness", "COULD-NOT-RUN", f"only {m['n_ok']}/{m['n']} runs ok")
        return
    for c in m["checks"]:
        label = f"local-only: {c['check']['type']} {c['check'].get('a', c['check'].get('measure'))}"
        record(label, c["state"], c["detail"])


def check_both_mixture():
    if NGSPICE_MISSING:
        record("both: mixture correlation matches analytic prediction", "COULD-NOT-RUN", "ngspice unavailable")
        return
    params = base_params(SIGMA_G, SIGMA_L, SIGMA_L)
    expect_r = SIGMA_G ** 2 / (SIGMA_G ** 2 + SIGMA_L ** 2)
    checks = [
        {"type": "correlation", "a": "r1", "b": "r2", "expect": expect_r},
        {"type": "sigma", "measure": "r1", "expect": 1000.0 * (SIGMA_G ** 2 + SIGMA_L ** 2) ** 0.5},
    ]
    spec, m = run_campaign("both_mixture", PAIR_TEMPLATE, params, ["r1", "r2"], checks)
    if m["n_ok"] < n_needed(m):
        record("both-mixture campaign completeness", "COULD-NOT-RUN", f"only {m['n_ok']}/{m['n']} runs ok")
        return
    for c in m["checks"]:
        label = f"both: {c['check']['type']} {c['check'].get('a', c['check'].get('measure'))} (analytic expect={c['check']['expect']:.4f})"
        record(label, c["state"], c["detail"])


def check_pelgrom():
    if NGSPICE_MISSING:
        record("Pelgrom scaling: sigma ~ 1/sqrt(W*L)", "COULD-NOT-RUN", "ngspice unavailable")
        return
    A = 0.02  # coefficient, arbitrary consistent units for this synthetic check
    small = {"name": "LOC", "scope": "local", "dist": "gauss", "nominal": 0.0,
              "pelgrom": {"coeff": A, "w": 1.0, "l": 1.0}}
    big = {"name": "LOC", "scope": "local", "dist": "gauss", "nominal": 0.0,
            "pelgrom": {"coeff": A, "w": 4.0, "l": 4.0}}
    spec_s, m_s = run_campaign("pelgrom_small", SINGLE_TEMPLATE, [small], ["r1"], [])
    spec_b, m_b = run_campaign("pelgrom_big", SINGLE_TEMPLATE, [big], ["r1"], [])
    if m_s["n_ok"] < n_needed(m_s) or m_b["n_ok"] < n_needed(m_b):
        record("Pelgrom scaling campaign completeness", "COULD-NOT-RUN",
               f"small {m_s['n_ok']}/{m_s['n']}, big {m_b['n_ok']}/{m_b['n']}")
        return
    sigma_small = m_s["stats"]["measures"]["r1"]["sigma"]
    sigma_big = m_b["stats"]["measures"]["r1"]["sigma"]
    ratio = sigma_small / sigma_big if sigma_big else float("inf")
    expect_ratio = (4.0 * 4.0 / (1.0 * 1.0)) ** 0.5  # sqrt(WL_big/WL_small) = sqrt(16) = 4
    rel_err = abs(ratio - expect_ratio) / expect_ratio
    # tolerance-based, not CI-combined (see docs/mc.md) -- 500-point sigma CI
    # width at this N is a few percent, 20% is a generous, honestly-stated
    # margin, not a tuned-to-pass number
    record("Pelgrom scaling: sigma(small WL)/sigma(big WL) ~= sqrt(WL_big/WL_small)",
           "pass" if rel_err <= 0.20 else "FAIL",
           f"sigma_small={sigma_small:.4f} sigma_big={sigma_big:.4f} ratio={ratio:.3f} "
           f"expect={expect_ratio:.3f} rel_err={rel_err:.1%}")


def check_reproducibility():
    if NGSPICE_MISSING:
        record("reproducibility: same seed -> identical results", "COULD-NOT-RUN", "ngspice unavailable")
        return
    params = base_params(SIGMA_G, SIGMA_L, SIGMA_L)
    n_small = 50
    spec_a, m_a = run_campaign("repro_a", PAIR_TEMPLATE, params, ["r1", "r2"], [], n=n_small, seed=777)
    spec_b, m_b = run_campaign("repro_b", PAIR_TEMPLATE, params, ["r1", "r2"], [], n=n_small, seed=777)
    if m_a["n_ok"] < n_small or m_b["n_ok"] < n_small:
        record("reproducibility campaign completeness", "COULD-NOT-RUN",
               f"a={m_a['n_ok']}/{n_small} b={m_b['n_ok']}/{n_small}")
        return
    mismatches = []
    for i in range(n_small):
        ra, rb = m_a["runs"][i], m_b["runs"][i]
        for pname in ("GA", "LA", "LB"):
            if ra["params"][pname] != rb["params"][pname]:
                mismatches.append(f"run {i} param {pname}: {ra['params'][pname]} != {rb['params'][pname]}")
        for mname in ("r1", "r2"):
            va, vb = ra["measures"].get(mname), rb["measures"].get(mname)
            if va != vb:
                mismatches.append(f"run {i} measure {mname}: {va} != {vb}")
    record("reproducibility: two campaigns, same seed -> bit-identical sampled params and measures",
           "pass" if not mismatches else "FAIL",
           f"{n_small} runs compared, {len(mismatches)} mismatches" + (f": {mismatches[:3]}" if mismatches else ""))

    # lone replay of one run, using only the spec + seed + run index -- not
    # the manifest's recorded values, and not any neighbouring run's state
    k = 17
    replay_dir = os.path.join(DATA_ROOT, "repro_a")
    manifest_path = os.path.join(replay_dir, "manifest.json")
    res = mc.replay_run(manifest_path, k, spec=spec_a)
    record(f"reproducibility: replaying run {k} alone reproduces its recorded result",
           "pass" if res["ok"] else "FAIL", str(res.get("reason", "matched")))


def check_red_proof():
    if NGSPICE_MISSING:
        record("red-proof: shared 'local' param must FAIL independence check", "COULD-NOT-RUN",
               "ngspice unavailable")
        return
    # GA sigma 0 so only the (buggy, shared) "local" contribution shows;
    # LA is used for BOTH r1 and r2 via PAIR_TEMPLATE_SHARED_BUG, i.e. this
    # is a deck that CLAIMS local/independent mismatch but actually wires
    # both devices to the same draw -- exactly the bug this tool's whole
    # design (distinct @@name@@ per instance) exists to make impossible to
    # fall into silently. The check below must catch it.
    params = [
        {"name": "GA", "scope": "global", "dist": "gauss", "nominal": 0.0, "sigma": 0.0},
        {"name": "LA", "scope": "local", "dist": "gauss", "nominal": 0.0, "sigma": SIGMA_L},
    ]
    checks = [{"type": "correlation", "a": "r1", "b": "r2", "expect": 0.0}]
    spec, m = run_campaign("red_proof", PAIR_TEMPLATE_SHARED_BUG, params, ["r1", "r2"], checks)
    if m["n_ok"] < n_needed(m):
        record("red-proof campaign completeness", "COULD-NOT-RUN", f"only {m['n_ok']}/{m['n']} runs ok")
        return
    c = m["checks"][0]
    corr = m["stats"]["correlations"].get("r1:r2", {})
    record("red-proof: independence check correctly reports FAIL on a shared 'local' param",
           "pass" if c["state"] == "FAIL" else "FAIL",
           f"checker state={c['state']} (expected FAIL) r={corr.get('r')} detail={c['detail']}")


DIFFPAIR_TEMPLATE_PATH = os.path.join(ROOT, "docs", "examples", "mc_diffpair_offset.cir")
DIFFPAIR_GAIN = 3.78598  # V/V, Vout_diff per V of injected gate offset, measured by
# finite difference at W=8u/L=1u, Vtail bias 100uA/pair -- see docs/mc.md
# "How the gain figure was measured". Non-gating: a real BSIM4 differential
# pair (docs/mc.md's worked example), included here as extra regression
# coverage, not the pass/FAIL gate (that's the resistor checks above).
DIFFPAIR_PELGROM_A = 4e-3  # V*um, illustrative coefficient for this synthetic model


def check_diffpair_worked_example():
    if NGSPICE_MISSING:
        record("worked example: MOS diff-pair offset vs Pelgrom prediction (non-gating)",
               "COULD-NOT-RUN", "ngspice unavailable", gating=False)
        return
    if not os.path.exists(DIFFPAIR_TEMPLATE_PATH):
        record("worked example: MOS diff-pair offset vs Pelgrom prediction (non-gating)",
               "COULD-NOT-RUN", f"missing {DIFFPAIR_TEMPLATE_PATH}", gating=False)
        return
    with open(DIFFPAIR_TEMPLATE_PATH) as f:
        template = f.read()

    def one(tag, w, l, n=300):
        sigma = DIFFPAIR_PELGROM_A / (w * l) ** 0.5
        params = [
            {"name": "MW", "scope": "fixed", "dist": "gauss", "nominal": w, "sigma": 0.0},
            {"name": "ML", "scope": "fixed", "dist": "gauss", "nominal": l, "sigma": 0.0},
            {"name": "VOFF", "scope": "local", "dist": "gauss", "nominal": 0.0, "sigma": sigma},
        ]
        spec, m = run_campaign(f"diffpair_{tag}", template, params, ["vdiff"], [], n=n, seed=99)
        return sigma, m

    try:
        sigma_small, m_small = one("small", 8, 1)
        sigma_big, m_big = one("big", 32, 4)
    except Exception as e:
        record("worked example: MOS diff-pair offset vs Pelgrom prediction (non-gating)",
               "COULD-NOT-RUN", f"exception: {e}", gating=False)
        return
    if m_small["n_ok"] < n_needed(m_small) or m_big["n_ok"] < n_needed(m_big):
        record("worked example: MOS diff-pair offset vs Pelgrom prediction (non-gating)",
               "COULD-NOT-RUN", f"small {m_small['n_ok']}/{m_small['n']} big {m_big['n_ok']}/{m_big['n']}",
               gating=False)
        return
    back_small = m_small["stats"]["measures"]["vdiff"]["sigma"] / DIFFPAIR_GAIN
    back_big = m_big["stats"]["measures"]["vdiff"]["sigma"] / DIFFPAIR_GAIN
    rel_err_small = abs(back_small - sigma_small) / sigma_small
    rel_err_big = abs(back_big - sigma_big) / sigma_big
    ok = rel_err_small <= 0.15 and rel_err_big <= 0.15
    record("worked example: MOS diff-pair offset vs Pelgrom prediction (non-gating)",
           "pass" if ok else "FAIL",
           f"W=8u/L=1u: predicted={sigma_small*1e3:.4f}mV back-converted={back_small*1e3:.4f}mV "
           f"(rel_err={rel_err_small:.1%}); W=32u/L=4u: predicted={sigma_big*1e3:.4f}mV "
           f"back-converted={back_big*1e3:.4f}mV (rel_err={rel_err_big:.1%})",
           gating=False)


def main():
    print("=" * 70)
    check_ngspice_available()

    print("=" * 70)
    print("1. global-only")
    check_global_only()

    print("=" * 70)
    print("2. local-only")
    check_local_only()

    print("=" * 70)
    print("3. both (mixture)")
    check_both_mixture()

    print("=" * 70)
    print("4. Pelgrom scaling")
    check_pelgrom()

    print("=" * 70)
    print("5. reproducibility")
    check_reproducibility()

    print("=" * 70)
    print("6. red-proof")
    check_red_proof()

    print("=" * 70)
    print("7. worked example (non-gating, docs/mc.md)")
    check_diffpair_worked_example()

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
