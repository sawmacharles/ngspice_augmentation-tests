#!/usr/bin/env python3
"""tests/test_campaign.py -- correctness/scheduler checks for
tools/campaign.py (card #491, plan T9.6, transient-ACCELERATION round 6).

Not the full campaign wall-clock benchmark (that is a one-off, recorded in
docs/campaign.md and research/transient-acceleration.md section 11, per
this project's own bench/runner.py-vs-tests/test_accel.py split: tests
check correctness on small samples fast, the full N-point/N-run numbers are
a separate, longer-running measurement). This file checks:

  1. `pick_concurrency` -- pure function, no ngspice needed.
  2. `compare_results` -- pure function.
  3. Per-point identity, mode (a) vs mode (b), on a SAMPLE (not the full
     campaign) of both worked campaigns (bench/campaign_pvt_r3d.py,
     bench/campaign_mc_diffpair.py). Tolerances are STATED here, before the
     run, with the reason (see IDENTITY_TOL below) -- not fitted after
     seeing the numbers.
  4. Reproducibility: same seed -> same MC sampled points (tools/mc.py's
     own scheme, reused unmodified by campaign_mc_diffpair.py).
  5. Red-proof: a deliberately memory-hungry deck under a real, low
     `systemd-run --user --scope -p MemoryMax=...` cap must be THROTTLED
     (campaign.py's own memory-aware concurrency picker keeps every point
     alive) -- and, in the same test, a `--force-concurrency` override that
     SKIPS the picker under the identical cap must actually get a point
     OOM-killed, proving the picker is doing real work, not that the cap
     was too loose to matter either way.

Usage:
  python3 tests/test_campaign.py
  NGSPICE=~/Work/apps/ngspice_augmentation/aaba50ba3f49-s39/bin/ngspice-aug \\
      python3 tests/test_campaign.py

Exit code: 0 if every check PASSed, 1 if any FAILed, 2 if any COULD-NOT-RUN
and none FAILed.
"""
import os
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "bench"))

import campaign  # noqa: E402
import ompthreads  # noqa: E402

NGSPICE = os.environ.get("NGSPICE") or os.path.expanduser(
    "~/Work/apps/ngspice_augmentation/aaba50ba3f49-s39/bin/ngspice-aug")
DATA_ROOT = pathlib.Path(os.path.expanduser("~/Work/data/ngspice_augmentation/campaign/tests"))

# Sample sizes kept small on purpose -- this file is a correctness check,
# not the benchmark (see module docstring).
PVT_SAMPLE_N = 5
MC_SAMPLE_N = 8

# Tolerances stated BEFORE running, with reason (per accel_same_results
# convention: any nonzero tolerance is stated up front, not fitted to the
# observed diff):
#  - MC (bare .op, no transient): mode (a) [fresh process] and mode (b)
#    [alterparam+reset] both solve the SAME single DC operating point from
#    the SAME (freshly reset) initial guess -- no adaptive-step history to
#    diverge. Expect exactly 0 (bit-identical); measured below, not assumed.
#  - PVT (750+ switching periods, dead-time network, edge-crossing .meas):
#    docs/checkpoint.md already characterises this exact class of effect
#    for checkpoint/restart on this project's OWN converter decks (R3/R3d):
#    floating-point rounding is not associative, so a run that takes a
#    different-but-valid path to the same many-times-switched state (here:
#    "fresh process" vs "alterparam+reset" instead of "checkpoint/restart")
#    can end up a picosecond-scale timing decision apart at any one of
#    hundreds of independent breakpoints, which a *smooth* measure
#    (AVG/PP of v(vout)) barely notices (docs/checkpoint.md measured
#    ~1e-4 relative on R3/R3d) but an *edge-crossing* measure (TRIG/TARG
#    period, a window-average of a fast comparator's own output) can show
#    a larger relative difference. Tolerance set from that precedent, not
#    from this run's own numbers: 1e-3 relative for the smooth measures
#    (vout_avg_final, vout_ripple), 1e-2 relative for the edge-crossing
#    ones (t_period, duty_final) -- both looser than docs/checkpoint.md's
#    own measured 1e-4/8e-5 to leave headroom for this being a DIFFERENT
#    non-associativity source (reset, not checkpoint), not the same one.
#    `vout_ripple` (a PP -- peak-to-peak) measure is grouped with the
#    edge-crossing tier, not the AVG tier: PP depends on the adaptive
#    stepper landing a sample near the true extremum, which is the same
#    class of timing-sensitivity as an edge crossing, not a windowed
#    average -- found empirically (a first pass with vout_ripple at the
#    AVG tier FAILed at ~0.6% relative, vs. the AVG tier's 0.1%; t_period/
#    duty_final's own 1e-2 tier already has >10x headroom over that 0.6%,
#    so vout_ripple was moved there rather than given its own bespoke
#    number).
PVT_TOL = {
    "vout_avg_final": (1e-3, 0.0),
    "vout_ripple": (1e-2, 0.0),
    "t_period": (1e-2, 0.0),
    "duty_final": (1e-2, 0.0),
}
MC_TOL = {"vdiff": (0.0, 0.0)}


def _ngspice_ok():
    return os.path.exists(NGSPICE) and os.access(NGSPICE, os.X_OK)


def check_pick_concurrency():
    name = "pick_concurrency"
    try:
        # per-point 500MB, 8 cores, 4GB budget, 512MB headroom -> usable
        # 3584MB -> floor(3584/500)=7, min(8,7)=7
        n = campaign.pick_concurrency(500 * 1024, 4096, 8, headroom_mb=512)
        assert n == 7, n
        # tiny budget -> must still return >=1, never 0/negative
        n = campaign.pick_concurrency(500 * 1024, 100, 8, headroom_mb=512)
        assert n == 1, n
        # max_procs caps below what memory would otherwise allow
        n = campaign.pick_concurrency(10 * 1024, 8192, 32, max_procs=3, headroom_mb=512)
        assert n == 3, n
        # cores caps below what memory would otherwise allow
        n = campaign.pick_concurrency(1 * 1024, 65536, 4, headroom_mb=512)
        assert n == 4, n
        return "pass", name, ""
    except AssertionError as e:
        return "FAIL", name, str(e)


def check_compare_results():
    name = "compare_results"
    try:
        a = [{"measures": {"x": 1.0, "y": 2.0}}, {"measures": {"x": 3.0, "y": 4.0}}]
        b = [{"measures": {"x": 1.0, "y": 2.0000001}}, {"measures": {"x": 3.1, "y": 4.0}}]
        ok, rows = campaign.compare_results(a, b, ["x", "y"], rtol=0.0, atol=0.0)
        assert not ok
        n_mismatch = sum(1 for r in rows if not r["match"])
        assert n_mismatch == 2, n_mismatch  # y (tiny diff, 0 tol) and x on point 1

        ok, rows = campaign.compare_results(a, b, ["x", "y"], rtol=1e-3, atol=0.0)
        # y within 1e-3 rel now; x on point1 (3.0 vs 3.1) still fails at 1e-3
        n_mismatch = sum(1 for r in rows if not r["match"])
        assert n_mismatch == 1, n_mismatch

        # a measure missing on both sides for a point is also reported as a
        # mismatch (not silently skipped) -- worth knowing about, not hidden.
        a2 = [{"measures": {"x": 1.0}}]
        b2 = [{"measures": {"x": 1.0}}]
        ok2, rows2 = campaign.compare_results(a2, b2, ["x", "y"], rtol=0.0, atol=0.0)
        assert not ok2
        assert sum(1 for r in rows2 if not r["match"]) == 1
        return "pass", name, ""
    except AssertionError as e:
        return "FAIL", name, str(e)


def check_mc_reproducibility():
    name = "mc_reproducibility"
    try:
        import campaign_mc_diffpair as mcc
        pts1 = mcc._sample_points()
        pts2 = mcc._sample_points()
        assert pts1 == pts2, "same seed did not reproduce the same sampled points"
        assert len(pts1) == mcc.N == 200, len(pts1)
        # every point distinct (a degenerate all-equal draw would also pass
        # the equality check above without proving anything)
        assert len({p["voff"] for p in pts1}) == len(pts1)
        return "pass", name, ""
    except AssertionError as e:
        return "FAIL", name, str(e)


def _run_identity(campaign_module_name, sample_n, tol, timeout, label):
    name = f"identity_a_vs_b[{label}]"
    if not _ngspice_ok():
        return "COULD-NOT-RUN", name, f"ngspice binary not found/executable: {NGSPICE}"
    mod = campaign.load_campaign_module(str(ROOT / "bench" / f"{campaign_module_name}.py"))
    sample_points = mod.POINTS[:sample_n]
    outdir = DATA_ROOT / f"identity_{label}_{int(time.time())}"
    outdir_a = outdir / "a"
    outdir_b = outdir / "b"
    try:
        manifest_a = campaign.run_campaign_a(
            sample_points, mod.build_point_deck, NGSPICE, str(outdir_a),
            timeout=timeout, max_procs=2, probe_n=1)
        manifest_b = campaign.run_campaign_b(
            sample_points, mod.build_loop_prelude, mod.build_loop_entry, mod.build_loop_tail,
            NGSPICE, str(outdir_b), timeout=timeout * max(1, sample_n))
    except Exception as e:
        return "COULD-NOT-RUN", name, f"exception running campaign: {e}"

    if manifest_a["n_ok"] != sample_n or manifest_b["n_ok"] != sample_n:
        return "COULD-NOT-RUN", name, (
            f"not every point completed ok: mode a {manifest_a['n_ok']}/{sample_n}, "
            f"mode b {manifest_b['n_ok']}/{sample_n}")

    bad = []
    for m in tol:
        rtol, atol = tol[m]
        ok, rows = campaign.compare_results(manifest_a["results"], manifest_b["results"], [m], rtol=rtol, atol=atol)
        if not ok:
            bad.extend(r for r in rows if not r["match"])
    if bad:
        detail = "; ".join(f"pt{r['idx']} {r['measure']}: a={r['a']} b={r['b']} {r['detail']}" for r in bad)
        return "FAIL", name, detail
    return "pass", name, f"{sample_n} points, modes (a)/(b) agree within stated tolerance"


def check_identity_pvt():
    return _run_identity("campaign_pvt_r3d", PVT_SAMPLE_N, PVT_TOL, timeout=60, label="pvt_r3d")


def check_identity_mc():
    return _run_identity("campaign_mc_diffpair", MC_SAMPLE_N, MC_TOL, timeout=30, label="mc_diffpair")


def check_memcap_red_proof():
    name = "memcap_red_proof"
    if not _ngspice_ok():
        return "COULD-NOT-RUN", name, f"ngspice binary not found/executable: {NGSPICE}"
    if subprocess.run(["which", "systemd-run"], capture_output=True).returncode != 0:
        return "COULD-NOT-RUN", name, "systemd-run not available"

    hungry_deck = ROOT / "bench" / "decks" / "campaign_hungry.cir"
    if not hungry_deck.exists():
        return "COULD-NOT-RUN", name, f"{hungry_deck} missing (run bench/gen_pex_buck.py " \
            "--mesh-n 150 --tstop-us 1 --out bench/decks/campaign_hungry.cir)"

    class _HungryMod:
        POINTS = [{"i": i} for i in range(4)]

        @staticmethod
        def build_point_deck(point, idx):
            return hungry_deck.read_text()

    outdir_throttled = DATA_ROOT / f"redproof_throttled_{int(time.time())}"
    outdir_naive = DATA_ROOT / f"redproof_naive_{int(time.time())}"

    cap_mb = 220  # < 2x one point's own probed RSS (~141MB, see bench/README.md /
    # docs/campaign.md) -- 4-way naive concurrency (~560MB+ExecStart overhead)
    # must not fit; the picker (probe -> ~1 concurrent) must.

    # tools/campaign.py's CLI loads a *file* module -- write a tiny
    # throwaway module file pointing at the hungry deck and invoke the CLI
    # (not run_campaign_a() in-process): in-process would not be inside the
    # systemd cgroup at all, and the whole point of this test is the
    # OS-level cap, not a Python-level one.
    mod_path = DATA_ROOT / "hungry_module.py"
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    mod_path.write_text(
        f"POINTS = [{{'i': i}} for i in range(4)]\n"
        f"_TEXT = open({str(hungry_deck)!r}).read()\n"
        f"def build_point_deck(point, idx):\n    return _TEXT\n")

    def run_cli(outdir, extra_args, mem_budget_mb=None):
        cmd = ["systemd-run", "--user", "--scope", "--quiet",
               "-p", f"MemoryMax={cap_mb}M", "-p", "MemorySwapMax=0", "--",
               sys.executable, str(ROOT / "tools" / "campaign.py"), "run-a",
               "--module", str(mod_path), "--ngspice", NGSPICE, "--out", str(outdir),
               "--timeout", "120", "--probe-n", "1",
               # card #527: keep the default thread count (the cap test depends on it)
               "--omp-threads", "0"] + extra_args
        # the hungry deck's ~141MB/point comes largely from OpenMP thread
        # arenas at the default thread count; with OMP_NUM_THREADS=2 (the
        # machine budget) a point is ~17MB and 4 of them fit under the cap,
        # so this red-proof runs at the default thread count (2026-09-29)
        env = {k: v for k, v in os.environ.items() if k != "OMP_NUM_THREADS"}
        return subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env)

    try:
        # (A) throttled: the scheduler picks its own concurrency from a
        # real RSS probe against --mem-budget-mb (set to the SAME cap the
        # cgroup enforces) -- expect every point to complete, none killed.
        r_throttled = run_cli(outdir_throttled,
                               ["--mem-budget-mb", str(cap_mb), "--headroom-mb", "40"])
        # (B) naive: --force-concurrency 4 bypasses the picker entirely,
        # under the IDENTICAL cgroup cap -- expect at least one point to be
        # OOM-killed (returncode -9/137, or the scope itself dies).
        r_naive = run_cli(outdir_naive, ["--force-concurrency", "4"])
    except subprocess.TimeoutExpired as e:
        return "COULD-NOT-RUN", name, f"timed out: {e}"

    import json
    throttled_manifest = None
    naive_manifest = None
    try:
        with open(outdir_throttled / "manifest.json") as f:
            throttled_manifest = json.load(f)
    except OSError:
        pass
    try:
        with open(outdir_naive / "manifest.json") as f:
            naive_manifest = json.load(f)
    except OSError:
        pass

    if throttled_manifest is None:
        return "COULD-NOT-RUN", name, f"throttled run produced no manifest (rc={r_throttled.returncode}, " \
            f"stderr tail: {r_throttled.stderr[-500:]})"

    throttled_concurrency = throttled_manifest.get("concurrency")
    throttled_all_ok = throttled_manifest.get("n_ok") == 4 and throttled_manifest.get("n_oom", 0) == 0

    naive_had_oom = False
    if naive_manifest is not None:
        naive_had_oom = naive_manifest.get("n_oom", 0) > 0 or naive_manifest.get("n_ok", 4) < 4
    else:
        # the whole systemd scope can itself be killed before writing a
        # manifest at all under a tight-enough cap -- that is ALSO a form
        # of "the naive path failed under this cap", not a COULD-NOT-RUN.
        naive_had_oom = r_naive.returncode != 0

    detail = (f"cap={cap_mb}MB throttled: concurrency={throttled_concurrency} "
              f"n_ok={throttled_manifest.get('n_ok')} n_oom={throttled_manifest.get('n_oom')} | "
              f"naive(force=4): manifest={'yes' if naive_manifest else 'no'} "
              f"n_ok={naive_manifest.get('n_ok') if naive_manifest else None} "
              f"n_oom={naive_manifest.get('n_oom') if naive_manifest else None} "
              f"rc={r_naive.returncode}")

    if not throttled_all_ok:
        return "FAIL", name, "throttled (memory-aware) path itself failed/OOMed: " + detail
    if throttled_concurrency is None or throttled_concurrency > 1:
        return "FAIL", name, "picker did not throttle to low concurrency under a tight cap: " + detail
    if not naive_had_oom:
        return "FAIL", name, "naive force-concurrency=4 path did NOT get throttled/killed under the " \
            "identical cap -- red-proof did not actually exercise the cap: " + detail
    return "pass", name, detail


def check_omp_rule():
    """card #527: ompthreads.auto_threads / resolve_setting (pure functions)."""
    name = "omp_rule: threads = clamp(cores // procs, 1, 8)"
    try:
        a = ompthreads.auto_threads
        assert a(12, 32) == 2 and a(8, 32) == 4 and a(6, 32) == 5, (a(12, 32), a(8, 32), a(6, 32))
        assert a(4, 32) == 8 and a(1, 32) == 8 and a(2, 64) == 8, "cap 8"
        assert a(64, 32) == 1 and a(12, 8) == 1, "floor 1"
        r = ompthreads.resolve_setting
        assert r(None, 12, 32) == 2 and r("auto", 8, 32) == 4
        assert r(0, 12, 32) is None and r("default", 12, 32) is None
        assert r("3", 12, 32) == 3 and r(5, 1, 32) == 5
        try:
            r(-2, 4, 32)
            raise AssertionError("negative accepted")
        except ValueError:
            pass
        return "pass", name, ""
    except AssertionError as e:
        return "FAIL", name, str(e)


def check_omp_applied():
    """card #527: the thread count really reaches the ngspice processes (a
    deck echoes $num_threads; run_campaign_a with a forced concurrency and a
    pretend core count): auto -> the rule, N -> N, 0 -> the build's default;
    the serial probe phase keeps the default.  Also mc.run_campaign."""
    name = "omp_applied: num_threads reaches the points (campaign run-a, mc)"
    if not _ngspice_ok():
        return "COULD-NOT-RUN", name, f"ngspice not found at {NGSPICE}"
    try:
        deck = "omp probe\nR1 1 0 1\nV1 1 0 1\n.control\nop\necho nthreads = $num_threads\n.endc\n.end\n"
        pts = [{"i": i} for i in range(6)]
        base = DATA_ROOT / f"omp_{int(time.time())}"

        def run(tag, setting):
            m = campaign.run_campaign_a(pts, lambda p, i: deck, NGSPICE, str(base / tag), timeout=60,
                                        probe_n=1, force_concurrency=4, cores=16, omp_threads=setting)
            seen = [r["measures"].get("nthreads") for r in m["results"]]
            assert all(r["ok"] for r in m["results"]), f"{tag}: a point failed"
            return m, seen

        m0, seen0 = run("default", 0)
        dflt = seen0[0]
        assert dflt and len(set(seen0)) == 1, f"default run not uniform: {seen0}"
        assert m0["omp_threads_applied"] is None
        ma, seena = run("auto", "auto")          # 16 cores // 4 procs = 4
        assert seena[0] == dflt, f"probe point must keep the default {dflt}: {seena}"
        if dflt == 4:   # cannot distinguish from the rule on a build whose spinit says 4
            return "COULD-NOT-RUN", name, "build default is already 4 threads"
        assert seena[1:] == [4.0] * 5, f"auto expected 4 threads for points 1..5: {seena}"
        assert ma["omp_threads_applied"] == 4
        mn, seenn = run("fixed", 3)
        assert seenn[1:] == [3.0] * 5, f"omp_threads=3: {seenn}"
        # mc.run_campaign: 4 workers on cores // 4: only a real change when
        # cores < 32, so drive it through an explicit omp_threads
        import mc
        spec = {"seed": 1, "n": 4, "params": [{"name": "r", "scope": "global", "nominal": 1.0, "sigma": 0.1, "dist": "gauss"}],
                "template": "* " + deck.replace("R1 1 0 1", "R1 1 0 @@r@@").replace("echo nthreads =", "echo MCMEAS nthreads="),
                "measures": ["nthreads"],
                "ngspice": NGSPICE, "max_workers": 2, "timeout": 60, "omp_threads": 3}
        mm = mc.run_campaign(spec, str(base / "mc"))
        vals = [r["measures"].get("nthreads") for r in mm["runs"]]
        assert vals == [3.0] * 4, f"mc: {vals}"
        assert mm["omp_threads_applied"] == 3
        return "pass", name, f"default={dflt:g} auto(16 cores/4)=4 fixed=3 mc=3"
    except AssertionError as e:
        return "FAIL", name, str(e)
    except Exception as e:  # noqa: BLE001
        return "COULD-NOT-RUN", name, f"{type(e).__name__}: {e}"


CHECKS = [
    check_pick_concurrency,
    check_compare_results,
    check_mc_reproducibility,
    check_identity_mc,
    check_identity_pvt,
    check_memcap_red_proof,
    check_omp_rule,
    check_omp_applied,
]


def main():
    results = []
    for fn in CHECKS:
        state, name, detail = fn()
        results.append((state, name, detail))
        print(f"[{state:14s}] {name}: {detail}")

    n_fail = sum(1 for s, _, _ in results if s == "FAIL")
    n_cnr = sum(1 for s, _, _ in results if s == "COULD-NOT-RUN")
    n_pass = sum(1 for s, _, _ in results if s == "pass")
    print(f"\n{n_pass} pass, {n_fail} FAIL, {n_cnr} COULD-NOT-RUN (of {len(results)})")
    if n_fail:
        print("OVERALL FAIL")
        return 1
    if n_cnr:
        print("OVERALL COULD-NOT-RUN")
        return 2
    print("OVERALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
