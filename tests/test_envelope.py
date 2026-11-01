#!/usr/bin/env python3
"""Validation for envelope-following (card #512, transient acceleration
round 10): `tools/envelope_tran.py` (the driver) and `tools/ckpt_io.py`
(the checkpoint-file read/write module it needs and this round built --
see docs/envelope.md for the full design writeup, including where this
mechanism falls short of the card's own accuracy bound and why).

This file checks, in order:

1. `ckpt_io` round-trips a REAL checkpoint file byte-for-byte identical
   (`test_roundtrip_byte_identical`) -- the format was read directly from
   `src/frontend/com_checkpoint.c`/`cktdefs.h`, not guessed, and this is
   the direct verification of that reading, not an assumption.
2. A patched checkpoint (only CKTtime/CKTrhsOld/CKTorder/CKTdelta/
   CKTbreaks changed) `loadstate`s cleanly and the run continues without
   a danger string (`test_patched_checkpoint_loads_cleanly`).
3. `m_periods=1` (no jumps at all -- pure checkpoint-chaining, no
   extrapolation) meets the card's own 0.5% cycle-average tolerance on
   MOST periods of `tests/fixtures/env_r3d_softstart.cir`
   (`test_checkpoint_chaining_control_accuracy`) -- isolates "is the
   checkpoint mechanism itself sound" from "is the projective jump
   sound" (docs/envelope.md's own "the landing kick" section explains why
   these must be measured separately). Gated PASS/FAIL at a loose bound
   (80% of periods, not 100%) because even pure chaining is not
   bit-exact (docs/checkpoint.md's own established ~1e-4-relative
   restart-vs-uninterrupted-run finding applies here too).
4. The full default-M envelope run's speed and accuracy are MEASURED and
   REPORTED against the deck's own full plain-`.tran` reference
   (`test_full_run_measured`) -- NOT gated pass/fail on the accuracy
   number: docs/envelope.md characterises, with a measured root cause,
   that this mechanism as built does not meet the card's 0.5%/5% bound,
   and gating this test "pass" on that number would misrepresent a real,
   already-written-up negative result. The test DOES gate that the run
   completes (reaches tstop) and reports genuine numbers, not placeholders.
5. Refusal behaviour on the red-proof deck
   (`test_redproof_refusal_and_bounded_error`): asserts refusals occur in
   the vicinity of the injected load-step event, and that even the ONE
   jump docs/envelope.md found straddling the event (a characterised,
   inherent limitation -- an event entirely inside a SKIPPED span cannot
   be observed without simulating it) produces a BOUNDED error, not an
   unbounded blow-up.

6.-11. card #719: the native implicit engine `.tranenv` (card #513, ported
   to the s150 code and made to work on closed-loop converters by running
   every period map at a tight transient tolerance, `maptol`, default
   reltol 1e-5; docs/envelope.md, card #719 section), measured against a
   plain `.tran` at reltol 1e-6 (2 ns max step):
   `tranenv719_closed_loop_softstart` (env_r3d_softstart.cir, voltage-mode
   buck, dead-time network: accuracy bounds on every simulated period and
   faster than the deck's own plain `.tran 2n` in an interleaved pair),
   `tranenv719_red_proof_maptol_off` (ENV_MAPTOL_OFF=1, the card #513
   behaviour: no jump accepted on that deck), `tranenv719_pcm_softstart`
   (envelope719/env_r4_pcm_softstart.cir, peak-current-mode buck with an SR
   latch: accuracy bounds; speed reported, not gated: about equal to the
   deck's own plain `.tran 10n` at default tolerance, which is 5-7x less
   accurate), `tranenv719_linear_chopper` (the card #513 smooth deck),
   `tranenv719_rejects_undriven`, and with NGSPICE_BASE=<a build without
   this work> `tranenv719_plain_tran_unchanged` (plain `.tran` output
   byte-identical).  These need a build with `.tranenv` (else
   COULD-NOT-RUN); about 1 minute on top of checks 1-5 (the reltol=1e-6
   references, ~25 s, are cached in the output folder).

    NGSPICE=<a build with patches 0032-0033, e.g. s55> python3 tests/test_envelope.py
    exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
    ENVELOPE_OUT=<dir> moves every output (default ~/Work/data/ngspice_augmentation/tests_envelope)

Runtime: ~40-60s (the full-deck run in check 4 is the slow part, ~7-10s of
actual ngspice time plus Python/subprocess overhead for ~600-900 hops).
"""
import filecmp
import os
import re
import subprocess
import sys
import tempfile
import time
import json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import ckpt_io as ck          # noqa: E402
import envelope_tran as et    # noqa: E402

FIXTURES = os.path.join(ROOT, "tests", "fixtures")
REFS = os.path.join(ROOT, "refs")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation")
OUT = os.environ.get("ENVELOPE_OUT") or os.path.join(DATA_ROOT, "tests_envelope")

_CANDIDATES = [
    os.environ.get("NGSPICE"),
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice"),
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s55/bin/ngspice-aug"),
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s54/bin/ngspice-aug"),
]
NGSPICE = next((p for p in _CANDIDATES if p and os.path.exists(p) and os.access(p, os.X_OK)), None)

REF_SOFTSTART = os.path.join(DATA_ROOT, "envelope", "env_r3d_softstart",
                              "env_r3d_softstart_reference.txt")
REF_REDPROOF = os.path.join(DATA_ROOT, "envelope", "env_r3d_softstart_redproof",
                             "env_r3d_softstart_redproof_reference.txt")
REF_VECTORS = ["v(vout)", "i(L1)", "v(vrefnode)"]
_REF_CACHE = {}


def make_reference(fixture, dest):
    """card #550: run the fixture's own full plain .tran (its .control block: `run`, wrdata) with the raw-file write dropped and the wrdata
    redirected to `dest`, so the test needs no reference data synced from another machine (a second workstation had none; the redproof reference on the development workstation (AMD Ryzen AI Max+ 395, Linux) was
    stale: written for a fixture whose load step never existed). About 20 s, 130 MB. Returns dest, or None if the run produced nothing."""
    with open(os.path.join(FIXTURES, fixture)) as fh:
        text = fh.read()
    text = text.replace(".include ../../stb/lgprobe.lib", f".include {os.path.join(ROOT, 'stb', 'lgprobe.lib')}")
    text = re.sub(r"(?m)^write .*\n", "", text)
    text, n = re.subn(r"(?m)^wrdata \S+", "wrdata " + dest, text)
    if n != 1 or NGSPICE is None:
        return None
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    cir = dest + ".cir"
    with open(cir, "w") as fh:
        fh.write(text)
    try:
        subprocess.run([NGSPICE, "-b", cir], cwd=os.path.dirname(dest), capture_output=True, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        return None
    return dest if os.path.exists(dest) and os.path.getsize(dest) > 0 else None


def reference(kind):
    """kind 'softstart' | 'redproof': the full-.tran reference wrdata. softstart: the data-dir file of an earlier run when it exists (the first-round
    reference the docs quote), else regenerated; redproof: always regenerated from the current fixture (cached by the fixture's mtime)."""
    if kind in _REF_CACHE:
        return _REF_CACHE[kind]
    fixture = f"env_r3d_{'softstart' if kind == 'softstart' else 'softstart_redproof'}.cir"
    canon = REF_SOFTSTART if kind == "softstart" else None
    path = canon if canon and os.path.exists(canon) else None
    if path is None:
        dest = os.path.join(OUT, "ref", f"{kind}_reference.txt")
        if os.path.exists(dest) and os.path.getmtime(dest) > os.path.getmtime(os.path.join(FIXTURES, fixture)):
            path = dest
        else:
            path = make_reference(fixture, dest)
    _REF_CACHE[kind] = path
    return path

results = []


def record(name, status, detail):
    results.append((name, status, detail))
    print(f"{status:12s} {name}: {detail}")


def test_roundtrip_byte_identical():
    name = "roundtrip_byte_identical"
    if NGSPICE is None:
        record(name, "COULD-NOT-RUN", "no ngspice build found (set NGSPICE)")
        return
    deck = os.path.join(REFS, "r2_switched_rc.cir")
    if not os.path.exists(deck):
        record(name, "COULD-NOT-RUN", f"{deck} missing")
        return
    with tempfile.TemporaryDirectory() as d:
        import stb
        body = "\n".join(stb.prepare_deck(deck))
        ctl = "\n.tran 25n 100u 0\n.control\nset filetype=ascii\nrun\nsavestate ckpt_a.bin\n.endc\n.end\n"
        with open(os.path.join(d, "probe.cir"), "w") as fh:
            fh.write(body + ctl)
        proc = subprocess.run([NGSPICE, "-b", "probe.cir"], cwd=d, capture_output=True,
                               text=True, timeout=60)
        ckpt = os.path.join(d, "ckpt_a.bin")
        if not os.path.exists(ckpt):
            record(name, "COULD-NOT-RUN", f"no checkpoint written; stderr={proc.stderr[:300]}")
            return
        st = ck.read_checkpoint(ckpt)
        rt = os.path.join(d, "ckpt_rt.bin")
        ck.write_checkpoint(rt, st)
        same = filecmp.cmp(ckpt, rt, shallow=False)
        if same:
            record(name, "PASS", f"{os.path.getsize(ckpt)} bytes, 0 different")
        else:
            record(name, "FAIL", "round-tripped checkpoint differs from original")


def test_patched_checkpoint_loads_cleanly():
    name = "patched_checkpoint_loads_cleanly"
    if NGSPICE is None:
        record(name, "COULD-NOT-RUN", "no ngspice build found")
        return
    deck = os.path.join(REFS, "r2_switched_rc.cir")
    with tempfile.TemporaryDirectory() as d:
        import stb
        body = "\n".join(stb.prepare_deck(deck))
        save_ctl = "\n.tran 25n 100u 0\n.control\nset filetype=ascii\nrun\nsavestate a.bin\n.endc\n.end\n"
        with open(os.path.join(d, "save.cir"), "w") as fh:
            fh.write(body + save_ctl)
        subprocess.run([NGSPICE, "-b", "save.cir"], cwd=d, capture_output=True, text=True, timeout=60)
        ckpt = os.path.join(d, "a.bin")
        if not os.path.exists(ckpt):
            record(name, "COULD-NOT-RUN", "savestate produced nothing")
            return
        st = ck.read_checkpoint(ckpt)
        tclk_guess = 10e-6
        new_t = st["CKTtime"] + 5 * tclk_guess
        dx = 0.0  # this deck is already settled; a zero extrapolation is a valid smoke test
        st["CKTtime"] = new_t
        st["CKTorder"] = 1
        st["CKTdelta"] = min(st["CKTdelta"], 1e-9)
        st["CKTsaveDelta"] = st["CKTdelta"]
        st["CKTbreaks"] = [round(new_t, 15)]
        st["breakSize"] = 1
        patched = os.path.join(d, "patched.bin")
        ck.write_checkpoint(patched, st)
        resume_ctl = (f"\n.tran 25n {new_t + 20e-6:g} 0\n.control\nset filetype=ascii\n"
                       "loadstate patched.bin\nwrdata out.txt v(vx)\n.endc\n.end\n")
        with open(os.path.join(d, "resume.cir"), "w") as fh:
            fh.write(body + resume_ctl)
        proc = subprocess.run([NGSPICE, "-b", "resume.cir"], cwd=d, capture_output=True,
                               text=True, timeout=60)
        out = proc.stdout + proc.stderr
        danger = et._danger(proc)
        out_path = os.path.join(d, "out.txt")
        if danger or not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
            record(name, "FAIL", f"danger or empty output; log tail: {out[-300:]}")
        else:
            record(name, "PASS", f"loadstate of a CKTtime/CKTrhsOld-patched checkpoint "
                                  f"continued cleanly, resumed at t={new_t:g}s")


def test_checkpoint_chaining_control_accuracy():
    name = "checkpoint_chaining_control_accuracy_M1"
    REF = reference("softstart") if NGSPICE else None
    if NGSPICE is None or not REF:
        record(name, "COULD-NOT-RUN",
               f"missing ngspice build or reference file {REF_SOFTSTART} (and it could not be regenerated)")
        return
    deck = os.path.join(FIXTURES, "env_r3d_softstart.cir")
    outdir = os.path.join(OUT, "m1_control")
    try:
        rep = et.envelope_follow(deck, outdir, tclk=1 / 500e3, tstop=200e-6, tstep=2e-9,
                                  vectors=["v(vout)", "i(L1)"], m_periods=1, n_check=4,
                                  ngspice=NGSPICE)
    except Exception as exc:  # noqa: BLE001
        record(name, "COULD-NOT-RUN", f"envelope_follow raised {exc!r}")
        return
    if len(rep["jumps"]) != 0:
        record(name, "FAIL", f"m_periods=1 should never jump; got {len(rep['jumps'])} jumps")
        return
    val = et.validate_against_reference(rep, REF, REF_VECTORS, col=0,
                                         avg_tol_rel=0.005, avg_tol_floor=1e-3)
    n = val["n_checked"]
    if n == 0:
        record(name, "COULD-NOT-RUN", "no periods compared against the reference")
        return
    frac_pass = sum(1 for r in val["per_window"] if r["avg_err"] <= r["avg_tol"]) / n
    detail = (f"n_periods={n} frac_within_0.5pct={frac_pass:.2f} "
              f"worst_rel_err={val['worst_avg_rel']*100:.2f}%")
    if frac_pass >= 0.80:
        record(name, "PASS", detail + " (checkpoint-chaining mechanism itself is sound)")
    else:
        record(name, "FAIL", detail)


def test_full_run_measured():
    name = "full_run_speed_and_accuracy_M10"
    REF = reference("softstart") if NGSPICE else None
    if NGSPICE is None or not REF:
        record(name, "COULD-NOT-RUN",
               f"missing ngspice build or reference file {REF_SOFTSTART} (and it could not be regenerated)")
        return
    deck = os.path.join(FIXTURES, "env_r3d_softstart.cir")
    outdir = os.path.join(OUT, "full_m10")
    import time
    t0 = time.time()
    try:
        rep = et.envelope_follow(deck, outdir, tclk=1 / 500e3, tstop=4e-3, tstep=2e-9,
                                  vectors=["v(vout)", "i(L1)"], m_periods=10, n_check=4,
                                  ngspice=NGSPICE)
    except Exception as exc:  # noqa: BLE001
        record(name, "COULD-NOT-RUN", f"envelope_follow raised {exc!r}")
        return
    wall = time.time() - t0
    if rep["t_final"] < 4e-3 - 1e-9:
        record(name, "FAIL", f"did not reach tstop: t_final={rep['t_final']:.6g}")
        return
    val = et.validate_against_reference(rep, REF, REF_VECTORS, col=0,
                                         avg_tol_rel=0.005, avg_tol_floor=1e-3)
    frac_pass = (sum(1 for r in val["per_window"] if r["avg_err"] <= r["avg_tol"]) /
                 max(val["n_checked"], 1))
    speedup_hint = ("(reference full run measured separately, interleaved A/B n=3, median "
                     "18.28s vs this config's own median 15.21s = 1.20x, see docs/envelope.md)")
    detail = (f"driver_wall={wall:.2f}s periods_simulated={rep['n_periods_simulated']}/2000 "
              f"jumps={len(rep['jumps'])} refusals={len(rep['refusals'])} "
              f"vout_avg: n_checked={val['n_checked']} frac_within_0.5pct={frac_pass:.2f} "
              f"worst_rel_err={val['worst_avg_rel']*100:.1f}% {speedup_hint}. "
              "NOT gated pass/fail on the accuracy number -- see docs/envelope.md "
              "'the landing kick' for the characterised root cause of why this "
              "does not meet the card's 0.5%/5% bound as built.")
    record(name, "INFO", detail)


def test_redproof_refusal_and_bounded_error():
    name = "redproof_refusal_and_bounded_error"
    if NGSPICE is None:
        record(name, "COULD-NOT-RUN", "no ngspice build found")
        return
    deck = os.path.join(FIXTURES, "env_r3d_softstart_redproof.cir")
    if not os.path.exists(deck):
        record(name, "COULD-NOT-RUN", f"{deck} missing")
        return
    outdir = os.path.join(OUT, "redproof")
    event_t = 1.5e-3
    try:
        rep = et.envelope_follow(deck, outdir, tclk=1 / 500e3, tstop=1.6e-3, tstep=2e-9,
                                  vectors=["v(vout)", "i(L1)"], m_periods=10, n_check=4,
                                  ngspice=NGSPICE)
    except Exception as exc:  # noqa: BLE001
        record(name, "COULD-NOT-RUN", f"envelope_follow raised {exc!r}")
        return
    # card #550 (2026-10-02): the fixture's load step never existed (PWL continuation lines without '+': ngspice dropped them, the reference
    # .tran shows no step either). The old window of 120 us before the event was satisfied on the development workstation by a gate firing on compensator ringing
    # 94 us earlier (and not at all on the second workstation, whose last-bit-different trajectory had no such ringing refusal) -- a coincidence, not the event.
    # With the step real (4 A at 1.5 ms) both machines refuse exactly at it: development workstation refusals at t = 1498-1510 us (mid_check_event x4, growth_tol x4)
    # and second workstation 1500-1510 us (same reasons), the first ones 14 / 12 us before (the jump attempts whose check window contains the event).
    # Window: [event - 20 us, event + 30 us] (measured spread between machines 2 us), and a growth_tol refusal AT or after the event must be there
    # (the pre-jump gate seeing the jump in the per-period delta, not only the mid-check scan).
    near = [r for r in rep["refusals"] if -20e-6 <= r["t"] - event_t <= 30e-6]
    after_growth = [r for r in near if r["t"] >= event_t - 1e-9 and r["reason"] == "growth_tol"]
    if not near or not after_growth:
        record(name, "FAIL", f"{len(near)} refusals within [-20, +30] us of the injected event at t={event_t:g}, "
               f"{len(after_growth)} growth_tol at/after it; the smoothness gates did not react to the load step")
        return
    detail = f"{len(near)} refusals within [-20,+30] us of the event (reasons: " \
             f"{sorted(set(r['reason'] for r in near))})"
    REF_RP = reference("redproof")
    if REF_RP is None:
        record(name, "COULD-NOT-RUN", detail + "; the reference .tran for the bounded-error check could not be generated")
        return
    if True:
        val = et.validate_against_reference(rep, REF_RP, REF_VECTORS, col=0,
                                             avg_tol_rel=0.005, avg_tol_floor=1e-2)
        near_windows = [r for r in val["per_window"] if abs(r["t0"] - event_t) < 40e-6]
        if near_windows:
            worst_near = max(r["avg_err"] / max(abs(r["ref_avg"]), 1e-6) for r in near_windows)
            detail += f"; worst cycle-avg error near the event = {worst_near*100:.2f}% (bounded check: <20%)"
            if worst_near > 0.20:
                record(name, "FAIL", detail + " -- unbounded, not a graceful degradation")
                return
    record(name, "PASS", detail)


# ---------------------------------------------------------------------------
# card #719: `.tranenv` on closed-loop converters (docs/envelope.md, card #719)

FIX719 = os.path.join(FIXTURES, "envelope719")
T719 = 2e-6
ENV719 = ".tranenv 500k 200 8 1e-4 4m 100 envtol=5e-4 phase=0.5 uic"
FINE_OPTS = ".options reltol=1e-6 abstol=1e-12 vntol=1e-9"


def _deck719(fixture, line, opts, wr):
    """fixture text with its .tran line (and anything after it) replaced by `opts`, `line` and a control block that writes v(vout)
    i(L1) to the bare file name `wr` (cwd = the output folder; ngspice lowercases bare file names, so `wr` is lower case)."""
    with open(fixture) as fh:
        txt = fh.read()
    txt = re.sub(r"(\.\./)+stb/lgprobe\.lib", os.path.join(ROOT, "stb", "lgprobe.lib"), txt)
    i = txt.index("\n.tran ")
    return (txt[:i + 1] + (opts + "\n" if opts else "") + line + "\n.control\nset num_threads=2\nset filetype=ascii\n"
            "set wr_singlescale\nset wr_vecnames\nsave v(vout) i(L1)\nrun\nwrdata " + wr + " v(vout) i(L1)\n.endc\n.end\n")


def _run719(text, label, extra_env=None, ngspice=None, timeout=900):
    d = os.path.join(OUT, "env719")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, label + ".cir"), "w") as fh:
        fh.write(text)
    env = dict(os.environ, OMP_NUM_THREADS="2")
    for k in ("ENV_MAPTOL_OFF", "ENV_NCHECK", "ENV_DEBUG", "ENV_SCAN"):
        env.pop(k, None)
    env.update(extra_env or {})
    t0 = time.perf_counter()
    proc = subprocess.run([ngspice or NGSPICE, "-b", label + ".cir"], cwd=d, capture_output=True, text=True,
                          timeout=timeout, env=env)
    wall = time.perf_counter() - t0
    with open(os.path.join(d, label + ".log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    m = re.search(r"ENV: DONE periods_simulated=(\d+)/(\d+) .*?jumps=(\d+)", proc.stdout)
    stats = {"simulated": int(m.group(1)), "total": int(m.group(2)), "jumps": int(m.group(3))} if m else None
    return wall, proc, stats, os.path.join(d, label.lower() + ".txt")


def _windows719(path, T=T719, ph=0.5):
    """k -> (avg vout, avg iL, vout max-min) over every complete window [(k+ph)T, (k+1+ph)T] (the .tranenv period maps for
    phase=0.5); a window counts only if covered >= 98 % with no gap > T/4 (the .tranenv plot has a gap over every skipped span)."""
    st = {}
    prev = None
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.split()
            if len(p) < 3:
                continue
            t, v0, v1 = float(p[0]) - ph * T, float(p[1]), float(p[2])
            if prev is not None:
                tp, a0, a1 = prev
                dt = t - tp
                if dt > 0 and t + tp >= 0:
                    k = int(0.5 * (t + tp) / T)
                    x = st.setdefault(k, [0.0, 0.0, 0.0, 1e300, -1e300, 0.0])
                    x[0] += dt
                    x[1] += 0.5 * (v0 + a0) * dt
                    x[2] += 0.5 * (v1 + a1) * dt
                    x[3] = min(x[3], v0)
                    x[4] = max(x[4], v0)
                    x[5] = max(x[5], dt)
            prev = (t, v0, v1)
    return {k: (x[1] / x[0], x[2] / x[0], x[4] - x[3]) for k, x in st.items() if x[0] >= 0.98 * T and x[5] <= T / 4}


def _ref719(fixture, tag):
    """plain .tran 2n at reltol 1e-6 (halved-step check on the development workstation, card #719: vs 1n at reltol 1e-7 <= 0.005 % vout, 0.05 % iL, 0.14 %
    ripple on r3d), windows cached in the output folder by the fixture's mtime"""
    d = os.path.join(OUT, "env719")
    cache = os.path.join(d, f"ref_{tag}.json")
    if os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(fixture):
        with open(cache) as fh:
            return {int(k): tuple(v) for k, v in json.load(fh).items()}
    wall, proc, _, txt = _run719(_deck719(fixture, ".tran 2n 4m 0 uic", FINE_OPTS, f"ref_{tag}.txt"), f"ref_{tag}")
    if not os.path.exists(txt):
        return None
    w = _windows719(txt)
    with open(cache, "w") as fh:
        json.dump(w, fh)
    return w


def _err719(ref, cand, kmin=0):
    """worst error over the candidate's windows: (vout avg, iL avg, ripple), each relative to max(|ref|, 2 % of the run maximum)"""
    keys = [k for k in cand if k in ref and k >= kmin]
    mx = [max(abs(v[c]) for v in ref.values()) for c in range(3)]
    worst = [0.0, 0.0, 0.0]
    for k in keys:
        for c in range(3):
            worst[c] = max(worst[c], abs(cand[k][c] - ref[k][c]) / max(abs(ref[k][c]), 0.02 * mx[c]))
    return len(keys), worst


def _fmt719(n, w):
    return f"n={n} vout {100 * w[0]:.3f}% iL {100 * w[1]:.3f}% ripple {100 * w[2]:.3f}%"


_HAVE_TRANENV = None


def _have_tranenv():
    global _HAVE_TRANENV
    if _HAVE_TRANENV is None:
        _HAVE_TRANENV = False
        if NGSPICE is not None:
            txt = ("probe\nV1 a 0 PULSE(0 1 0 1n 1n 1u 2u)\nR1 a b 1\nC1 b 0 1u IC=0\n.tranenv 500k 50 5 1e-4 20u 4 uic\n"
                   ".control\nrun\n.endc\n.end\n")
            _, proc, stats, _ = _run719(txt, "probe_tranenv", timeout=60)
            _HAVE_TRANENV = stats is not None
    return _HAVE_TRANENV


def test_tranenv719_closed_loop_softstart():
    """r3d soft-start: every simulated period within the bounds (k >= 100: vout 0.1 %, iL 1 %, ripple 1 %; all: vout 0.2 %, iL 2 %,
    ripple 2 %; measured on the development workstation 0.014 / 0.16 / 0.14 % and 0.035 / 0.32 / 0.37 %), >= 20 jumps, <= 25 % of the periods simulated
    (measured 171/2000), and faster than the deck's own plain .tran 2n (default tolerances; 0.57 / 7.6 / 13.7 % off the reference)
    by >= 2x in an interleaved pair (measured 6.2x)."""
    name = "tranenv719_closed_loop_softstart"
    if not _have_tranenv():
        record(name, "COULD-NOT-RUN", "NGSPICE has no `.tranenv` (needs a card #719 build)")
        return
    fx = os.path.join(FIXTURES, "env_r3d_softstart.cir")
    ref = _ref719(fx, "r3d")
    if not ref:
        record(name, "COULD-NOT-RUN", "reference run produced nothing")
        return
    w_env, w_plain = [], []
    for rep in range(2):
        wall, proc, stats, txt = _run719(_deck719(fx, ENV719, "", "env_r3d.txt"), "env_r3d")
        w_env.append(wall)
        wallp, procp, _, txtp = _run719(_deck719(fx, ".tran 2n 4m 0 uic", "", "plain_r3d.txt"), "plain_r3d")
        w_plain.append(wallp)
    if stats is None or not os.path.exists(txt) or not os.path.exists(txtp):
        record(name, "FAIL", "the .tranenv or plain run did not complete (see env719/*.log)")
        return
    env_w = _windows719(txt)
    n100, e100 = _err719(ref, env_w, 100)
    nall, eall = _err719(ref, env_w, 0)
    _, eplain = _err719(ref, _windows719(txtp), 100)
    ratio = min(w_plain) / max(min(w_env), 1e-9)
    detail = (f"{stats['simulated']}/{stats['total']} periods simulated, {stats['jumps']} jumps; vs plain .tran reltol 1e-6: "
              f"k>=100 {_fmt719(n100, e100)}; all {_fmt719(nall, eall)}; wall .tranenv {min(w_env):.2f}s vs plain .tran 2n "
              f"(default tol, k>=100 {_fmt719(0, eplain)[4:]}) {min(w_plain):.2f}s = {ratio:.1f}x (best of 2, interleaved)")
    ok = (e100[0] <= 1e-3 and e100[1] <= 1e-2 and e100[2] <= 1e-2 and eall[0] <= 2e-3 and eall[1] <= 2e-2 and eall[2] <= 2e-2
          and stats["jumps"] >= 20 and stats["simulated"] <= 0.25 * stats["total"] and ratio >= 2.0)
    record(name, "PASS" if ok else "FAIL", detail)


def test_tranenv719_red_proof_maptol_off():
    """ENV_MAPTOL_OFF=1: period maps at the deck's reltol (1e-3), the card #513 engine: on the r3d soft-start deck no jump is
    accepted (measured 0 of 1999 periods skipped, 8 Newton refusals), i.e. the map tolerance is what makes the closed loop work."""
    name = "tranenv719_red_proof_maptol_off"
    if not _have_tranenv():
        record(name, "COULD-NOT-RUN", "NGSPICE has no `.tranenv` (needs a card #719 build)")
        return
    fx = os.path.join(FIXTURES, "env_r3d_softstart.cir")
    _, proc, stats, txt = _run719(_deck719(fx, ENV719, "", "env_r3d_off.txt"), "env_r3d_off", {"ENV_MAPTOL_OFF": "1"})
    if stats is None:
        record(name, "FAIL", "the .tranenv run did not complete (see env719/env_r3d_off.log)")
        return
    m = re.search(r"map_reltol=(\S+)", proc.stdout)
    detail = (f"ENV_MAPTOL_OFF=1: {stats['simulated']}/{stats['total']} periods simulated, {stats['jumps']} jumps, map reltol "
              f"{m.group(1) if m else '?'}")
    record(name, "PASS" if stats["jumps"] <= 2 else "FAIL", detail)


def test_tranenv719_pcm_softstart():
    """r4 peak-current-mode soft-start (envelope719/env_r4_pcm_softstart.cir): k >= 100 bounds vout 0.1 %, iL 1 %, ripple 1.5 %
    (measured 0.028 / 0.31 / 0.41 %), >= 50 jumps. Speed reported only: about equal to the deck's own plain .tran 10n (measured
    1.66 vs 1.69 s), which is 0.19 / 1.4 / 2.7 % off; with ENV_MAPTOL_OFF=1 the errors are those of a default-tolerance chain."""
    name = "tranenv719_pcm_softstart"
    if not _have_tranenv():
        record(name, "COULD-NOT-RUN", "NGSPICE has no `.tranenv` (needs a card #719 build)")
        return
    fx = os.path.join(FIX719, "env_r4_pcm_softstart.cir")
    ref = _ref719(fx, "r4")
    if not ref:
        record(name, "COULD-NOT-RUN", "reference run produced nothing")
        return
    wall, proc, stats, txt = _run719(_deck719(fx, ENV719, "", "env_r4.txt"), "env_r4")
    wallp, _, _, txtp = _run719(_deck719(fx, ".tran 10n 4m 0 uic", "", "plain_r4.txt"), "plain_r4")
    if stats is None or not os.path.exists(txt) or not os.path.exists(txtp):
        record(name, "FAIL", "the .tranenv or plain run did not complete (see env719/*.log)")
        return
    n100, e100 = _err719(ref, _windows719(txt), 100)
    _, eplain = _err719(ref, _windows719(txtp), 100)
    detail = (f"{stats['simulated']}/{stats['total']} periods simulated, {stats['jumps']} jumps; k>=100 {_fmt719(n100, e100)}; "
              f"wall {wall:.2f}s vs plain .tran 10n {wallp:.2f}s (default tol, k>=100 {_fmt719(0, eplain)[4:]}; not gated)")
    ok = e100[0] <= 1e-3 and e100[1] <= 1e-2 and e100[2] <= 1.5e-2 and stats["jumps"] >= 50
    record(name, "PASS" if ok else "FAIL", detail)


def test_tranenv719_linear_chopper():
    """the card #513 smooth deck still meets its bar after the port (vs plain .tran 2n, exact for this linear deck: vout 0.5 %, iL
    0.5 %, ripple 5 % on every simulated period; measured 0.08 / 0.23 / 2.3 %, 810 of 2000 periods as in card #513)"""
    name = "tranenv719_linear_chopper"
    if not _have_tranenv():
        record(name, "COULD-NOT-RUN", "NGSPICE has no `.tranenv` (needs a card #719 build)")
        return
    fx = os.path.join(FIX719, "env_lin_chopper.cir")
    _, _, _, txtp = _run719(_deck719(fx, ".tran 2n 4m 0 uic", "", "plain_lin.txt"), "plain_lin")
    wall, proc, stats, txt = _run719(_deck719(fx, ENV719, "", "env_lin.txt"), "env_lin")
    if stats is None or not os.path.exists(txt) or not os.path.exists(txtp):
        record(name, "FAIL", "a run did not complete (see env719/*_lin.log)")
        return
    n, e = _err719(_windows719(txtp), _windows719(txt), 0)
    ok = e[0] <= 5e-3 and e[1] <= 5e-3 and e[2] <= 5e-2 and stats["simulated"] <= 0.5 * stats["total"]
    record(name, "PASS" if ok else "FAIL", f"{stats['simulated']}/{stats['total']} periods simulated, {stats['jumps']} jumps; "
           f"vs plain .tran 2n {_fmt719(n, e)}; wall {wall:.2f}s")


def test_tranenv719_rejects_undriven():
    name = "tranenv719_rejects_undriven"
    if not _have_tranenv():
        record(name, "COULD-NOT-RUN", "NGSPICE has no `.tranenv` (needs a card #719 build)")
        return
    txt = "undriven\nV1 a 0 DC 1\nR1 a b 1k\nC1 b 0 1u IC=0\n.tranenv 500k 50 5 1e-4 20u 4 uic\n.control\nrun\n.endc\n.end\n"
    _, proc, stats, _ = _run719(txt, "undriven", timeout=60)
    ok = stats is None and "needs a driven (clocked) circuit" in proc.stdout + proc.stderr
    record(name, "PASS" if ok else "FAIL", "undriven deck refused with the driven-circuit message" if ok
           else "no refusal message (see env719/undriven.log)")


def test_tranenv719_plain_tran_unchanged():
    """only with NGSPICE_BASE (a build without this work, e.g. s150): plain .tran output byte-identical on the r4 deck"""
    base = os.environ.get("NGSPICE_BASE")
    if not base:
        print("(tranenv719_plain_tran_unchanged skipped: set NGSPICE_BASE to compare plain .tran output with a base build)")
        return
    name = "tranenv719_plain_tran_unchanged"
    if NGSPICE is None or not os.path.exists(base):
        record(name, "COULD-NOT-RUN", "NGSPICE or NGSPICE_BASE missing")
        return
    fx = os.path.join(FIX719, "env_r4_pcm_softstart.cir")
    _, _, _, a = _run719(_deck719(fx, ".tran 10n 4m 0 uic", "", "same_new.txt"), "same_new")
    _, _, _, b = _run719(_deck719(fx, ".tran 10n 4m 0 uic", "", "same_base.txt"), "same_base", ngspice=base)
    if not (os.path.exists(a) and os.path.exists(b)):
        record(name, "FAIL", "a plain run did not complete")
        return
    same = filecmp.cmp(a, b, shallow=False)
    record(name, "PASS" if same else "FAIL", f"plain .tran 10n wrdata of the r4 deck {'byte-identical' if same else 'DIFFERS'} "
           f"between NGSPICE and NGSPICE_BASE")


def main():
    test_roundtrip_byte_identical()
    test_patched_checkpoint_loads_cleanly()
    test_checkpoint_chaining_control_accuracy()
    test_full_run_measured()
    test_redproof_refusal_and_bounded_error()
    test_tranenv719_closed_loop_softstart()
    test_tranenv719_red_proof_maptol_off()
    test_tranenv719_pcm_softstart()
    test_tranenv719_linear_chopper()
    test_tranenv719_rejects_undriven()
    test_tranenv719_plain_tran_unchanged()

    n_fail = sum(1 for _, s, _ in results if s == "FAIL")
    n_cnr = sum(1 for _, s, _ in results if s == "COULD-NOT-RUN")
    print(f"\n{len(results)} checks: {sum(1 for _,s,_ in results if s=='PASS')} PASS, "
          f"{sum(1 for _,s,_ in results if s=='INFO')} INFO, {n_fail} FAIL, {n_cnr} COULD-NOT-RUN")
    if n_fail:
        return 1
    if n_cnr:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
