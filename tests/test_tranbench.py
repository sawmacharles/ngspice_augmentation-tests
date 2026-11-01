#!/usr/bin/env python3
"""Unit tests for tools/tranbench.py (card #641, transient benchmark harness).

Most checks need no ngspice: number/.tran parsing, deck splitting (control block, probe /
refopts headers, tran + meas lifted from a .control block), reference-deck construction
(tight options, tmax rule), rusage / warning / failure-kind parsing, .meas parsing, wrdata
reading, linear interpolation and error math (hand-computed cases), reference-cache key
independence from --env, summary statistics and the `compare` output.

One smoke run (RC deck through run_bench, reference + measured run, --env passed to the measured
run only) happens when TRANBENCH_NGSPICE names an ngspice binary; otherwise it is skipped.

Usage:  [TRANBENCH_NGSPICE=<binary>] python3 tests/test_tranbench.py
Smoke outputs go to ~/Work/data/ngspice_augmentation/641/test_tranbench/, never the repo.
Exit status: 0 = all pass (skips allowed), 1 = a check failed.
"""
import argparse
import csv
import io
import math
import os
import sys
import tempfile
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import tranbench as tb  # noqa: E402

results = []


def check(name):
    def deco(fn):
        try:
            r = fn()
            results.append((name, "skip" if r == "skip" else "pass", ""))
        except AssertionError as e:
            results.append((name, "FAIL", str(e) or traceback.format_exc(limit=1)))
        except Exception:
            results.append((name, "FAIL", traceback.format_exc(limit=3)))
        return fn
    return deco


def near(a, b, rel=1e-9, abs_=1e-12):
    return abs(a - b) <= max(abs_, rel * abs(b))


DECK = """* title line -- test deck
* probe: v(out) i(v1)
* refopts: itl4=100
.param tstop=2m
.model d1 D(IS=1e-14)
V1 in 0 PULSE(0 1 0 1n 1n 1m 2m)
R1 in out 1k
C1 out 0 1u
.options reltol=1e-3
.print tran v(out)
.save v(out)
.meas tran v1 FIND v(out) AT=1m
.tran 1u {tstop}
.control
run
.endc
.end
"""


@check("spice_num suffixes")
def _():
    assert near(tb.spice_num("10n"), 1e-8) and near(tb.spice_num("1meg"), 1e6) and near(tb.spice_num("2.5k"), 2500)
    assert near(tb.spice_num("1e-3"), 1e-3) and tb.spice_num("abc") is None and near(tb.spice_num("5"), 5)


@check("eval_value params and expressions")
def _():
    p = {"fs": "500k", "tclk": "{1/fs}"}
    assert near(tb.eval_value("{tclk*10}", p), 2e-5)
    assert tb.eval_value("{nosuch}", p) is None


@check("parse_tran_card: tstep tstop tstart tmax uic, params")
def _():
    t = tb.parse_tran_card(".tran 1n 10u 2u 5n uic")
    assert t["uic"] and all(near(a, b) for a, b in ((t["tstep"], 1e-9), (t["tstop"], 1e-5), (t["tstart"], 2e-6), (t["tmax"], 5e-9)))
    t = tb.parse_tran_card(".tran 1u {tstop}", {"tstop": "2m"})
    assert near(t["tstop"], 2e-3) and t["tmax"] is None and not t["uic"] and t["tstart"] == 0.0
    try:
        tb.parse_tran_card(".tran 1u")
    except ValueError:
        return
    raise AssertionError("no ValueError")


@check("parse_deck: headers, dropped cards, meas conversion")
def _():
    d = tb.parse_deck(DECK)
    assert d["probes"] == "v(out) i(v1)" and d["refopts"] == "itl4=100"
    text = "\n".join(d["body"])
    for gone in (".print", ".save", ".tran", ".control", ".end", ".meas"):
        assert gone not in text, gone
    assert ".options reltol=1e-3" in text and ".model d1" in text and d["body"][0].startswith("* title")
    assert d["meas"] == ["v1"] and d["meas_lines"] == ["meas tran v1 FIND v(out) AT=1m"]
    assert near(d["tran"]["tstop"], 2e-3)


@check("parse_deck: tran and meas lifted from a .control block")
def _():
    d = tb.parse_deck("* t\nV1 a 0 1\nR1 a 0 1k\n.control\nset filetype=ascii\ntran 2n 1.5m 0 uic\n"
                      "meas tran vavg AVG v(a) FROM=1m TO=1.5m\nprint vavg\nwrite x.raw\n.endc\n.end\n")
    assert near(d["tran"]["tstop"], 1.5e-3) and d["tran"]["uic"] and d["meas"] == ["vavg"]
    assert not any("write" in ln or "print" in ln for ln in d["body"])


@check("parse_deck: no .tran is an error")
def _():
    try:
        tb.parse_deck("* t\nR1 a 0 1\n.end\n")
    except ValueError:
        return
    raise AssertionError("expected ValueError")


@check("reference deck: tight options, tmax = tstep/4, span cap, deck tmax kept, no env anywhere")
def _():
    d = tb.parse_deck(DECK)
    ref = tb.build_run_deck(d, "reltol=0.1", "/x/ref.dat", d["probes"], ref_rung=1)
    assert ".options reltol=1e-6 abstol=1e-15 vntol=1e-9 trtol=1 itl4=100" in ref
    assert "reltol=0.1" not in ref                      # the measured run's options never reach the reference
    assert ".tran 1e-06 0.002 0 2.5e-07" in ref, ref
    assert "wrdata /x/ref.dat v(out) i(v1)" in ref and "rusage tranpoints trantime traniter rejected" in ref
    run = tb.build_run_deck(d, "reltol=0.1", "/x/run.dat", d["probes"])
    assert ".options reltol=0.1" in run and ".tran 1e-06 0.002 0" in run and "reltol=1e-6" not in run
    big = tb.parse_tran_card(".tran 1n 1m")
    assert near(tb.ref_tmax(big), 1e-3 / 50000)           # 2.5e-10 would be 4M points: capped
    fixed = tb.parse_tran_card(".tran 1n 1m 0 3n")
    assert near(tb.ref_tmax(fixed), 3e-9)
    assert "method=gear" in tb.build_run_deck(d, "", "/x", d["probes"], ref_rung=3)
    assert "reltol=1e-5" in tb.build_run_deck(d, "", "/x", d["probes"], ref_rung=4)


@check("parse_log: rusage counts, skipped warnings, failure kinds")
def _():
    out = "Transient timepoints = 610\n\nTransient analysis time = 0.0006\n\nTransient iterations = 1218\n\nRejected timepoints = 7\n"
    r = tb.parse_log(out, "Warning: time point skipped at t = 1e-9 (step 1e-21, 3 in a row)\n" * 3)
    assert (r["points"], r["iters"], r["rejected"], r["skipped"], r["fail_kind"]) == (610, 1218, 7, 3, "")
    assert near(r["tran_s"], 0.0006)
    assert tb.parse_log("", "doAnalyses: TRAN:  Timestep too small; time = 1e-9")["fail_kind"] == "timestep too small"
    assert tb.parse_log("", "Error: Transient op failed, timestep too small")["fail_kind"] == "timestep too small"
    assert tb.parse_log("", "Warning: singular matrix:  check nodes x and y")["fail_kind"] == "no convergence"
    assert tb.parse_log("", "Error: no such vector sw")["fail_kind"] == "other"
    assert tb.parse_log(out, "")["fail_kind"] == ""


@check("parse_meas: values, failed, missing; last assignment wins")
def _():
    out = "v1                  =  6.32120e-01\nbad = failed\nv1 = 6.321204000000000e-01\n"
    m = tb.parse_meas(out, ["v1", "bad", "gone"])
    assert m == {"v1": 0.6321204, "bad": None, "gone": None}


@check("interp: linear, clamping, duplicate time stamps")
def _():
    ts, ys = [0.0, 1.0, 1.0, 3.0], [0.0, 2.0, 5.0, 1.0]
    assert near(tb.interp(ts, ys, 0.5), 1.0) and near(tb.interp(ts, ys, 2.0), 3.0)
    assert tb.interp(ts, ys, -1) == 0.0 and tb.interp(ts, ys, 9) == 1.0 and tb.interp([2.0], [7.0], 5) == 7.0


@check("signal_error: constant offset, normalisation, pp floor")
def _():
    t = [i / 100.0 for i in range(101)]
    ref = (t, [math.sin(2 * math.pi * x) for x in t])        # pp = 2
    run = (t, [math.sin(2 * math.pi * x) + 0.01 for x in t])
    e = tb.signal_error(run, ref, tb.uniform_grid(0, 1, 400))
    assert near(e["abs_max"], 0.01, 1e-9) and near(e["abs_rms"], 0.01, 1e-9)
    assert near(e["nerr_max"], 0.005, 1e-3) and near(e["pp"], 2.0, 1e-3) and e["covered"] == 1.0
    flat = ([0, 1], [0.0, 0.0])
    e2 = tb.signal_error(([0, 1], [1e-3, 1e-3]), flat, tb.uniform_grid(0, 1, 11))
    assert near(e2["nerr_max"], 1.0)                          # 1 mV error / floor 1 mV


@check("signal_error: resampling a coarser waveform, partial coverage, non-finite")
def _():
    ref = ([0.0, 1.0], [0.0, 1.0])                            # ramp, pp 1
    run = ([0.0, 0.5], [0.0, 0.5])                            # same ramp, stops at 0.5
    e = tb.signal_error(run, ref, tb.uniform_grid(0, 1, 101))
    assert near(e["abs_max"], 0.0, abs_=1e-12) and near(e["covered"], 51 / 101, 1e-9) and e["n"] == 51
    coarse = ([0.0, 1.0], [0.0, 1.0])                         # triangle sampled at the corners only: error at the peak
    tri = ([0.0, 0.5, 1.0], [0.0, 1.0, 0.0])
    e = tb.signal_error(coarse, tri, tb.uniform_grid(0, 1, 3))
    assert near(e["abs_max"], 1.0) and near(e["abs_rms"], math.sqrt(1.25 / 3))   # errors 0, 0.5, 1 at t = 0, 0.5, 1
    bad = tb.signal_error(([0, 1], [float("nan"), 0.0]), ref, tb.uniform_grid(0, 1, 5))
    assert math.isinf(bad["abs_max"]) and math.isinf(bad["nerr_max"])
    none = tb.signal_error(([5.0, 6.0], [0, 0]), ref, tb.uniform_grid(0, 1, 5))
    assert none["covered"] == 0.0 and none["abs_max"] is None


@check("meas_errors: relative error, missing counted")
def _():
    worst, missing, per = tb.meas_errors({"a": 1.1, "b": None, "c": 2.0}, {"a": 1.0, "b": 3.0, "c": 2.0, "d": None})
    assert near(worst, 0.1) and missing == 1 and per["b"] is None and per["c"] == 0.0 and "d" not in per


@check("load_wrdata: (time,value) column pairs and shared-time layout")
def _():
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "a.dat")
        open(p, "w").write(" 0 1 0 10 \n 1 2 1 20 \n bad line\n 2 3 2 30\n")
        s = tb.load_wrdata(p, 2)
        assert s[0] == ([0, 1, 2], [1, 2, 3]) and s[1] == ([0, 1, 2], [10, 20, 30])
        open(p, "w").write("0 1 10\n1 2 20\n")
        s = tb.load_wrdata(p, 2)
        assert s[1] == ([0, 1], [10, 20])
        assert tb.load_wrdata(os.path.join(td, "none.dat"), 1) is None
        open(p, "w").close()
        assert tb.load_wrdata(p, 1) is None


@check("compare_waveforms picks the worst signal")
def _():
    t = [i / 10.0 for i in range(11)]
    ref = [(t, [x for x in t]), (t, [2 * x for x in t])]
    run = [(t, [x + 0.001 for x in t]), (t, [2 * x + 0.1 for x in t])]
    deck, per = tb.compare_waveforms(run, ref, ["v(a)", "v(b)"], 0.0, 1.0, 50)
    assert deck["worst_sig"] == "v(b)" and near(deck["nerr_max"], 0.05, 1e-6) and near(per["v(a)"]["nerr_max"], 0.001, 1e-6)


@check("clean_env drops inherited TR_* and keeps explicit --env")
def _():
    os.environ["TR_STEP"] = "iter"
    try:
        e = tb.clean_env({})
        assert "TR_STEP" not in e and e["OMP_NUM_THREADS"] == "1"
        assert tb.clean_env({"TR_SKIP": "1"})["TR_SKIP"] == "1"
    finally:
        del os.environ["TR_STEP"]


@check("reference cache key ignores --env / --options, follows the deck and probes")
def _():
    d = tb.parse_deck(DECK)
    base = dict(ngspice=sys.executable)
    k1 = tb.ref_key(dict(base, env={"TR_STEP": "iter"}, options="reltol=0.1"), d, "v(out)")
    k2 = tb.ref_key(dict(base, env={}, options=""), d, "v(out)")
    assert k1 == k2
    assert tb.ref_key(base, d, "v(out) v(in)") != k2
    d2 = tb.parse_deck(DECK.replace("1k", "2k"))
    assert tb.ref_key(base, d2, "v(out)") != k2


def fake_rows():
    def row(deck, completed, kind, pts, rej, it, wall, nerr, ref_ok=1, skipped=0):
        return dict(deck=deck, completed=completed, fail_kind=kind, points=pts, rejected=rej, iters=it, wall_s=wall,
                    tran_s="0", skipped=skipped, reached="1", nerr_max=nerr, nerr_rms=nerr, abs_max="", abs_rms="",
                    worst_sig="v(a)", meas_err="", meas_missing="", ref_ok=ref_ok, ref_method="tight" if ref_ok else "")
    A = [row("d1", 1, "", 100, 5, 300, "1.0", "1e-3"), row("d2", 1, "", 200, 0, 500, "2.0", "1e-2"),
         row("d3", 0, "timestep too small", 50, 9, 120, "0.5", "", 0), row("d4", 1, "", 80, 1, 160, "0.4", "2e-4")]
    B = [row("d1", 1, "", 90, 2, 270, "0.9", "1.2e-3"), row("d2", 0, "timeout", 20, 0, 40, "300", "", 1, 5),
         row("d3", 1, "", 70, 3, 200, "0.6", "5e-3", 0), row("d4", 1, "", 80, 1, 160, "0.4", "9e-4")]
    return A, B


def write_dir(td, name, rows):
    d = os.path.join(td, name)
    os.makedirs(d)
    with open(os.path.join(d, "results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, tb.CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return d


@check("summarize: totals, median / max normalised error, failed list")
def _():
    A, _b = fake_rows()
    s = tb.summarize(A)
    assert s["decks"] == 4 and s["completed"] == 3 and s["completed_pct"] == 75.0
    assert s["points"] == 430 and s["rejected"] == 15 and s["iters"] == 1080
    assert near(s["nerr_median"], 1e-3) and near(s["nerr_max"], 1e-2) and s["nerr_n"] == 3
    assert s["failed"] == ["d3 (timestep too small)"] and s["no_ref"] == ["d3"] and s["fail_kinds"] == {"timestep too small": 1}
    buf = io.StringIO()
    old, sys.stdout = sys.stdout, buf
    try:
        tb.print_summary(dict(s, label="x"))
    finally:
        sys.stdout = old
    assert "completed 3 (75.0%)" in buf.getvalue() and "d3 (timestep too small)" in buf.getvalue()


@check("compare: flips, error regressions, per-deck table")
def _():
    A, B = fake_rows()
    with tempfile.TemporaryDirectory() as td:
        da, db = write_dir(td, "A", A), write_dir(td, "B", B)
        buf = io.StringIO()
        assert tb.compare(da, db, out=buf) == 0
    txt = buf.getvalue()
    assert "completed  A=3  B=3  both=2" in txt
    assert "LOST:timeout" in txt and "GAINED(was timestep too small)" in txt
    assert "flips completed->failed (A->B): 1 ['d2']" in txt and "flips failed->completed (A->B): 1 ['d3']" in txt
    assert "error up >2x: 1 ['d4']" in txt and "ERR UP" in txt
    assert "iterations A=460 B=430" in txt          # over d1 + d4 only (completed by both)


def smoke():
    nsp = os.environ.get("TRANBENCH_NGSPICE")
    if not nsp:
        return "skip"
    assert os.access(nsp, os.X_OK), f"TRANBENCH_NGSPICE not executable: {nsp}"
    out = os.path.expanduser("~/Work/data/ngspice_augmentation/641/test_tranbench")
    os.makedirs(out, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out) as td:
        cdir = os.path.join(td, "corpus")
        os.makedirs(cdir)
        open(os.path.join(cdir, "rc.cir"), "w").write(
            "* rc smoke\n* probe: v(out)\nV1 in 0 PULSE(0 1 100u 1n 1n 2m 4m)\nR1 in out 1k\nC1 out 0 1u\n"
            ".meas tran v1 FIND v(out) AT=1.1m\n.tran 5u 5m\n.end\n")
        base = dict(ngspice=nsp, label="smoke", env=None, options="", decks=None, corpus_dir=cdir, no_fixtures=True,
                    jobs=2, refcache=os.path.join(td, "refcache"), timeout=60.0, grid=2000, no_limit=True, list=False)
        o1 = os.path.join(td, "o1")
        assert tb.run_bench(argparse.Namespace(out=o1, **base)) == 0
        rows = tb.load(o1)
        r = rows["rc"]
        assert r["completed"] == "1" and r["ref_ok"] == "1" and r["ref_method"] == "tight", r
        assert int(r["points"]) > 100 and int(r["iters"]) >= int(r["points"]) and float(r["nerr_max"]) < 1e-3, r
        assert float(r["meas_err"]) < 1e-3
        # measured run with an env var and a sloppy tolerance: reference is reused (cached), env reaches the run only
        o2 = os.path.join(td, "o2")
        assert tb.run_bench(argparse.Namespace(out=o2, **dict(base, env=["TR_TESTONLY=1"], options="reltol=0.05"))) == 0
        r2 = tb.load(o2)["rc"]
        assert r2["completed"] == "1" and int(r2["points"]) <= int(r["points"]), (r2, r)
        assert os.path.exists(os.path.join(td, "refcache", "rc.json"))
        buf = io.StringIO()
        tb.compare(o1, o2, out=buf)
        assert "rc" in buf.getvalue()


check("smoke run on an RC deck (needs TRANBENCH_NGSPICE)")(smoke)


def main():
    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:5s} {name:{width}s}  {detail}")
    npass = sum(s == "pass" for _, s, _ in results)
    print(f"\n{npass}/{len(results)} pass, {sum(s == 'skip' for _, s, _ in results)} skipped, "
          f"{sum(s == 'FAIL' for _, s, _ in results)} FAIL")
    return 1 if any(s == "FAIL" for _, s, _ in results) else 0


if __name__ == "__main__":
    sys.exit(main())
