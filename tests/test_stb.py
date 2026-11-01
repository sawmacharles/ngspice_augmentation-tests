#!/usr/bin/env python3
"""Validation for stb_script (spec) — verify.md step 5, card #451, 2026-09-27.

Every check is pass / FAIL / COULD-NOT-RUN. The references are independent of the
Tian probe maths:
  R1a  closed-form return ratio from the deck's .param values;
  R1b  a plain-Python nodal solve of the broken loop's two-port Y parameters;
  R1b, R1d  Y parameters measured in ngspice with the probe replaced by two port
       sources held at the closed-loop DC voltage (a different circuit and formula).
Red-proof checks mutate the result (sign flip, swapped injections, bogus ngspice)
and pass only if the comparison then fails.

    python3 tests/test_stb.py        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import cmath
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import stb  # noqa: E402
import stbcore  # noqa: E402

OUT = os.path.join(stb.DATA_ROOT, "tests", "stb")
SWEEP = "dec 50 1 1e10"
TOL = 1e-6
results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))


def params(deck):
    vals = {}
    with open(deck) as fh:
        for line in fh:
            if line.lower().startswith(".param"):
                for k, v in re.findall(r"(\w+)=(\S+)", line):
                    vals[k.lower()] = spice_num(v)
    return vals


def spice_num(s):
    m = re.match(r"^([-+0-9.eE]+)(meg|[tgkmunpf])?$", s.lower())
    mult = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}
    return float(m.group(1)) * (mult[m.group(2)] if m.group(2) else 1.0)


# ---------- independent references ----------

def r1a_T(p, f):
    s = 2j * math.pi * f
    zo = 1 / (1 / p["ro"] + s * p["co"])
    zb = 1 / (1 / p["r2"] + s * p["cn"])
    return p["gm"] * zo * zb / (zo + p["r1"] + zb)


def solve(a, b):
    """Complex Gaussian elimination with partial pivoting, a x = b."""
    n = len(a)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[piv] = m[piv], m[c]
        for r in range(c + 1, n):
            k = m[r][c] / m[c][c]
            for j in range(c, n + 1):
                m[r][j] -= k * m[c][j]
    x = [0j] * n
    for r in range(n - 1, -1, -1):
        x[r] = (m[r][n] - sum(m[r][j] * x[j] for j in range(r + 1, n))) / m[r][r]
    return x


def r1b_T(p, f):
    """Two-port Y of R1b broken at oa | ol by nodal analysis; T = (Y12+Y21)/(Y11+Y22)."""
    s = 2j * math.pi * f
    n1, oa, ol, vn = range(4)
    y = [[0j] * 4 for _ in range(4)]

    def adm(i, j, g):
        for a in (i, j):
            if a is not None:
                y[a][a] += g
        if i is not None and j is not None:
            y[i][j] -= g
            y[j][i] -= g

    adm(n1, None, 1 / p["r1o"] + s * p["c1"])
    adm(oa, None, 1 / p["r2o"])
    adm(n1, oa, s * p["cc"])
    adm(ol, None, s * p["cl"])
    adm(ol, vn, 1 / p["ro1"])
    adm(vn, None, 1 / p["ro2"] + s * p["cin"])
    adm(vn, n1, s * p["cf"])
    y[n1][vn] += -p["gm1"]  # G1 n1 0 vp vn: current gm1*(vp-vn) leaves n1, vp = 0
    y[oa][n1] += p["gm2"]   # G2 oa 0 n1 0: current gm2*v(n1) leaves oa
    ports, internal = [oa, ol], [n1, vn]
    ymat = [[0j, 0j], [0j, 0j]]
    for c, pc in enumerate(ports):
        rhs = [-y[i][pc] for i in internal]
        vi = solve([[y[i][j] for j in internal] for i in internal], rhs)
        v = [0j] * 4
        v[pc] = 1
        for k, i in enumerate(internal):
            v[i] = vi[k]
        for r, pr in enumerate(ports):
            ymat[r][c] = sum(y[pr][j] * v[j] for j in range(4))
    return stbcore.tian_from_y(ymat[0][0], ymat[0][1], ymat[1][0], ymat[1][1])


def ngspice_y_reference(deck, probe, outdir):
    """T from ngspice-measured two-port Y parameters, ports held at the closed-loop DC voltage."""
    os.makedirs(outdir, exist_ok=True)
    lines = stb.prepare_deck(deck)
    a, b = stb.probe_nodes(lines, probe)
    # closed-loop DC voltage at the probe
    op_file = os.path.join(outdir, "op.txt")
    ctl = f"\n.control\nset wr_vecnames\nset numdgt=15\nop\nwrdata {op_file} v({a})\nquit\n.endc\n.end\n"
    w = os.path.join(outdir, "op.cir")
    with open(w, "w") as fh:
        fh.write("\n".join(lines) + ctl)
    subprocess.run([stb.NGSPICE, "-b", w], capture_output=True, text=True, timeout=600)
    with open(op_file) as fh:
        next(fh)
        vdc = float(fh.readline().split()[-1])
    opened = []
    for line in lines:
        f = line.split()
        if f and f[0].lower() == probe.lower():
            opened += [f"Vpa {a} 0 DC {vdc!r} AC 1", f"Vpb {b} 0 DC {vdc!r} AC 0"]
        else:
            opened.append(line)
    r1, r2 = os.path.join(outdir, "y1.txt"), os.path.join(outdir, "y2.txt")
    ctl = (f"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nop\nac {SWEEP}\n"
           f"alter vpa ac=0\nalter vpb ac=1\nac {SWEEP}\n"
           f"wrdata {r1} ac1.i(vpa) ac1.i(vpb)\nwrdata {r2} ac2.i(vpa) ac2.i(vpb)\nquit\n.endc\n.end\n")
    w = os.path.join(outdir, "yref.cir")
    with open(w, "w") as fh:
        fh.write("\n".join(opened) + ctl)
    proc = subprocess.run([stb.NGSPICE, "-b", w], capture_output=True, text=True, timeout=600)
    with open(os.path.join(outdir, "yref.log"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    f1, (ia1, ib1) = stb.read_wrdata(r1, 2)
    _f2, (ia2, ib2) = stb.read_wrdata(r2, 2)
    # current into the network = -(SPICE source current); unit port voltage
    return f1, [stbcore.tian_from_y(-y11, -y12, -y21, -y22)
                for y11, y21, y12, y22 in zip(ia1, ib1, ia2, ib2)]


def continuous_margins(fn, fmin, fmax):
    """PM/GM of a continuous T(f) by a dense scan (2000/decade) plus bisection."""
    n = int(2000 * math.log10(fmax / fmin))
    fs = [fmin * (fmax / fmin) ** (k / n) for k in range(n + 1)]
    ts = [fn(f) for f in fs]
    ph = stbcore.unwrapped_phase_deg(ts)
    out = {"pm": None, "gm": None}

    def bisect(g, lo, hi):
        glo = g(lo)
        for _ in range(80):
            mid = math.sqrt(lo * hi)
            if (g(mid) > 0) == (glo > 0):
                lo, glo = mid, g(mid)
            else:
                hi = mid
        return math.sqrt(lo * hi)

    for k in range(n):
        if (abs(ts[k]) - 1) * (abs(ts[k + 1]) - 1) < 0 and out["pm"] is None:
            fc = bisect(lambda f: abs(fn(f)) - 1, fs[k], fs[k + 1])
            base = ph[k] - math.degrees(cmath.phase(ts[k]))
            p = math.degrees(cmath.phase(fn(fc))) + base
            p += 360 * round((ph[k] - p) / 360)
            out["pm"] = (fc, 180 + p)
        if (ph[k] + 180) * (ph[k + 1] + 180) < 0 and out["gm"] is None:
            fc = bisect(lambda f: math.degrees(cmath.phase(-fn(f))), fs[k], fs[k + 1])
            out["gm"] = (fc, -20 * math.log10(abs(fn(fc))))
    return out


# ---------- checks ----------

def compare(name, got, ref, tol=TOL):
    err = stbcore.max_rel_err(got, ref)
    record(name, err <= tol, f"max rel err {err:.2e} (bound {tol:.0e})")
    return err


def expect_fail(name, got, ref, tol=TOL):
    err = stbcore.max_rel_err(got, ref)
    record(name, err > tol, f"mutated result must differ: max rel err {err:.2e}")


def margin_check(name, got, ref):
    for kind, key in (("pm", "pm_deg"), ("gm", "gm_db")):
        g, r = got[kind], ref[kind]
        if r is None or g is None:
            record(f"{name} {kind}", g is None and r is None,
                   f"reference {'none' if r is None else r}, script {'none' if g is None else g[key]}")
            continue
        bound = 0.1 if kind == "pm" else 0.05
        dv, dfreq = abs(g[key] - r[1]), abs(g["freq"] / r[0] - 1)
        record(f"{name} {kind}", dv <= bound and dfreq <= 0.01,
               f"script {g[key]:.4f} @ {g['freq']:.5g} Hz, reference {r[1]:.4f} @ {r[0]:.5g} Hz "
               f"(|d|={dv:.4f}, bound {bound}; freq {dfreq:.1e}, bound 1e-2)")


def run_deck(stem):
    deck = os.path.join(ROOT, "refs", f"{stem}.cir")
    return deck, stb.run(deck, "XP", SWEEP, os.path.join(OUT, stem))


def check_low_f_sign(name, margins, t0):
    """Card #464 sanity check: a genuine negative-feedback reference loop must
    have Re(T) > 0 at the lowest swept frequency, and stbcore.margins() must
    not raise the low_f_positive_feedback flag for it."""
    record(f"{name}: Re(T(lowest f)) > 0 (negative feedback at DC)", t0.real > 0, f"T(f_min) = {t0!r}")
    record(f"{name}: low_f_positive_feedback is False", margins["low_f_positive_feedback"] is False,
           f"note={margins['note']!r}")


def check_low_f_warns(name, margins):
    """Red-proof: a genuine positive-feedback loop must raise the flag."""
    record(f"{name}: low_f_positive_feedback is True (red-proof)", margins["low_f_positive_feedback"] is True,
           f"note={margins['note']!r}")


def check_low_f_integrator(name, margins):
    """Card #476: a healthy near-integrator loop (arg T(fmin) in -90+/-30 deg)
    must get the neutral low_f_integrator label, and must NOT raise the
    positive-feedback warning or its low-gain note."""
    record(f"{name}: low_f_integrator is True", margins["low_f_integrator"] is True,
           f"low_f_integrator_note={margins['low_f_integrator_note']!r}")
    record(f"{name}: low_f_positive_feedback is False (not the warning)",
           margins["low_f_positive_feedback"] is False, f"note={margins['note']!r}")
    record(f"{name}: low_f_positive_feedback_low_gain is False",
           margins["low_f_positive_feedback_low_gain"] is False,
           f"low_gain_note={margins['low_f_positive_feedback_low_gain_note']!r}")


def check_low_f_low_gain_note(name, margins):
    """Card #476: positive feedback (|arg T(fmin)| > 150 deg) with |T(fmin)| <= 1
    must raise low_f_positive_feedback_low_gain (a note), never the warning
    (low_f_positive_feedback) or the integrator label."""
    record(f"{name}: low_f_positive_feedback_low_gain is True",
           margins["low_f_positive_feedback_low_gain"] is True,
           f"note={margins['low_f_positive_feedback_low_gain_note']!r}")
    record(f"{name}: low_f_positive_feedback is False (not the warning)",
           margins["low_f_positive_feedback"] is False, f"note={margins['note']!r}")
    record(f"{name}: low_f_integrator is False",
           margins["low_f_integrator"] is False, f"integrator_note={margins['low_f_integrator_note']!r}")


def main():
    try:
        # R1a: closed form
        deck, res = run_deck("r1a_unilateral")
        p = params(deck)
        ref = [r1a_T(p, f) for f in res["freqs"]]
        compare("R1a T vs closed form", res["T"], ref)
        compare("R1a Middlebrook = Tian (no reverse transmission)", res["middlebrook"], res["T"])
        record("R1a loopGain = -T", all(lg == -t for lg, t in zip(res["loopGain"], res["T"])), "Spectre sign")
        margin_check("R1a margins vs bisection on closed form", res["margins"],
                     continuous_margins(lambda f: r1a_T(p, f), 1, 1e10))
        expect_fail("R1a red-proof: sign-flipped T", [-z for z in res["T"]], ref)
        r = res["raw"]
        swapped = [stbcore.tian(*z) for z in zip(r["v2"], r["i2"], r["v1"], r["i1"])]
        expect_fail("R1a red-proof: injections swapped", swapped, ref)
        check_low_f_sign("R1a", res["margins"], res["T"][0])

        # R1a probe orientation: T must not change when a and b are swapped
        text = open(deck).read().replace("XP vn fb lgprobe", "XP fb vn lgprobe")
        rev = os.path.join(OUT, "r1a_reversed.cir")
        os.makedirs(OUT, exist_ok=True)
        with open(rev, "w") as fh:
            fh.write(text.replace(".include ../stb/lgprobe.lib", f".include {ROOT}/stb/lgprobe.lib"))
        res_rev = stb.run(rev, "XP", SWEEP, os.path.join(OUT, "r1a_reversed"))
        compare("R1a probe reversed = forward", res_rev["T"], res["T"])

        # R1b: Python nodal reference and ngspice Y reference
        deck, res = run_deck("r1b_bilateral")
        p = params(deck)
        ref = [r1b_T(p, f) for f in res["freqs"]]
        compare("R1b T vs Python nodal Y-parameter reference", res["T"], ref)
        fy, refy = ngspice_y_reference(deck, "XP", os.path.join(OUT, "r1b_yref"))
        compare("R1b T vs ngspice Y-parameter reference", res["T"], refy)
        mb_err = stbcore.max_rel_err(res["middlebrook"], res["T"])
        record("R1b Middlebrook differs from Tian (reverse transmission present)", mb_err > 1e-3,
               f"max rel diff {mb_err:.2e}")
        margin_check("R1b margins vs bisection on Python reference", res["margins"],
                     continuous_margins(lambda f: r1b_T(p, f), 1, 1e10))
        check_low_f_sign("R1b", res["margins"], res["T"][0])

        # R1d: transistor level, ngspice Y reference with ports held at the DC voltage
        deck, res = run_deck("r1d_mos_follower")
        fy, refy = ngspice_y_reference(deck, "XP", os.path.join(OUT, "r1d_yref"))
        compare("R1d T vs ngspice Y-parameter reference (DC-held ports)", res["T"], refy)
        live = abs(res["T"][0])
        record("R1d signal live: low-frequency |T| > 10", live > 10, f"|T(1 Hz)| = {live:.4g}")
        check_low_f_sign("R1d", res["margins"], res["T"][0])

        # R1c: differential and common-mode loop gain via the diffprobe (card #451
        # subtask 2.2.4), each checked against an independent half-circuit
        # equivalent built directly by the DM/CM half-circuit theorem (exact for
        # this exactly symmetric linear network -- see stb/diffprobe.lib for why
        # reading branch 1 alone gives the whole coupled network's DM/CM T).
        deck = os.path.join(ROOT, "refs", "r1c_fully_diff.cir")
        res_dm = stb.run_diff(deck, "XP", "dm", SWEEP, os.path.join(OUT, "r1c_dm"))
        res_cm = stb.run_diff(deck, "XP", "cm", SWEEP, os.path.join(OUT, "r1c_cm"))
        half_dm = os.path.join(ROOT, "tests", "fixtures", "r1c_dm_half.cir")
        half_cm = os.path.join(ROOT, "tests", "fixtures", "r1c_cm_half.cir")
        ref_dm = stb.run(half_dm, "XP", SWEEP, os.path.join(OUT, "r1c_dm_half"))
        ref_cm = stb.run(half_cm, "XP", SWEEP, os.path.join(OUT, "r1c_cm_half"))
        compare("R1c DM loop gain vs DM half-circuit equivalent", res_dm["T"], ref_dm["T"])
        compare("R1c CM loop gain vs CM half-circuit equivalent", res_cm["T"], ref_cm["T"])
        # red-proof: the DM and CM loop gains are genuinely different circuits (different
        # DC gain and bandwidth -- see refs/r1c_fully_diff.cir); comparing one mode's
        # result against the other mode's reference (as if the diffprobe orientation or
        # E/F sign were wrong) must fail, not pass by coincidence.
        expect_fail("R1c red-proof: DM result checked against CM reference", res_dm["T"], ref_cm["T"])
        expect_fail("R1c red-proof: CM result checked against DM reference", res_cm["T"], ref_dm["T"])
        check_low_f_sign("R1c DM", res_dm["margins"], res_dm["T"][0])
        check_low_f_sign("R1c CM", res_cm["margins"], res_cm["T"][0])

        # card #464 red-proof: the ARCHIVED buggy R1c (positive feedback in both
        # loops, found by the parent 2026-09-29 -- see refs/r1c_fully_diff.cir's
        # header) must trip the low_f_positive_feedback warning in both modes;
        # the fixed R1c above (checked with check_low_f_sign, just above) must not.
        posfb_deck = os.path.join(ROOT, "tests", "fixtures", "stb4_r1c_posfb.cir")
        posfb_dm = stb.run_diff(posfb_deck, "XP", "dm", SWEEP, os.path.join(OUT, "r1c_posfb_dm"))
        posfb_cm = stb.run_diff(posfb_deck, "XP", "cm", SWEEP, os.path.join(OUT, "r1c_posfb_cm"))
        check_low_f_warns("R1c archived positive-feedback DM", posfb_dm["margins"])
        check_low_f_warns("R1c archived positive-feedback CM", posfb_cm["margins"])

        # -----------------------------------------------------------------
        # card #476: patch 0011 / card #464's low_f_positive_feedback flag
        # (Re(T) < 0 at fmin alone) false-triggered on healthy integrator
        # loops -- found on tests/fixtures/compopt_r3_averaged_detuned.cir
        # (a type-III compensated buck, near-integrator R1in/C1f branch) at
        # fstart=10 Hz: T is genuinely near -90 deg there, and a little extra
        # lag makes Re(T) slightly negative with no real instability.
        # DECISION (option A): warn only when |arg T(fmin)| > 150 deg
        # AND |T(fmin)| > 1; a neutral "integrator-like low-frequency
        # behaviour" label when arg T(fmin) is within -90 +/- 30 deg; and a
        # note (not the warning) when |arg T(fmin)| > 150 deg but
        # |T(fmin)| <= 1. Re-run every reference loop above (and the
        # compopt R3 fixture, detuned and tuned) at fstart=10 Hz specifically
        # -- the exact frequency the false positive was reported at -- rather
        # than relying only on SWEEP's fstart=1 Hz above.
        SWEEP10 = "dec 50 10 1e10"

        r1a_10 = stb.run(os.path.join(ROOT, "refs", "r1a_unilateral.cir"), "XP", SWEEP10,
                          os.path.join(OUT, "r1a_fstart10"))
        check_low_f_sign("R1a @ fstart=10Hz", r1a_10["margins"], r1a_10["T"][0])

        r1b_10 = stb.run(os.path.join(ROOT, "refs", "r1b_bilateral.cir"), "XP", SWEEP10,
                          os.path.join(OUT, "r1b_fstart10"))
        check_low_f_sign("R1b @ fstart=10Hz", r1b_10["margins"], r1b_10["T"][0])

        r1d_10 = stb.run(os.path.join(ROOT, "refs", "r1d_mos_follower.cir"), "XP", SWEEP10,
                          os.path.join(OUT, "r1d_fstart10"))
        check_low_f_sign("R1d @ fstart=10Hz", r1d_10["margins"], r1d_10["T"][0])

        r1c_dm_10 = stb.run_diff(deck, "XP", "dm", SWEEP10, os.path.join(OUT, "r1c_dm_fstart10"))
        r1c_cm_10 = stb.run_diff(deck, "XP", "cm", SWEEP10, os.path.join(OUT, "r1c_cm_fstart10"))
        check_low_f_sign("R1c DM (corrected) @ fstart=10Hz", r1c_dm_10["margins"], r1c_dm_10["T"][0])
        check_low_f_sign("R1c CM (corrected) @ fstart=10Hz", r1c_cm_10["margins"], r1c_cm_10["T"][0])

        # the ARCHIVED buggy R1c must still warn at fstart=10 Hz too (not
        # just at SWEEP's 1 Hz starting point) -- red-proof this fix did not
        # narrow the flag so far it stops catching the real bug.
        posfb_dm_10 = stb.run_diff(posfb_deck, "XP", "dm", SWEEP10, os.path.join(OUT, "r1c_posfb_dm_fstart10"))
        posfb_cm_10 = stb.run_diff(posfb_deck, "XP", "cm", SWEEP10, os.path.join(OUT, "r1c_posfb_cm_fstart10"))
        check_low_f_warns("R1c archived positive-feedback DM @ fstart=10Hz", posfb_dm_10["margins"])
        check_low_f_warns("R1c archived positive-feedback CM @ fstart=10Hz", posfb_cm_10["margins"])

        # compopt_r3_averaged_detuned.cir: the near-integrator type-III
        # compensator that motivated this card. Its own DETUNED starting
        # point (the fixture's own default r2=100k/c2=0.3n) must NOT raise
        # low_f_positive_feedback at fstart=10 Hz -- this was the reported
        # false positive (fixture header, docs/compopt.md's old "Sweep
        # note", now corrected).
        r3_detuned_deck = os.path.join(ROOT, "tests", "fixtures", "compopt_r3_averaged_detuned.cir")
        r3_detuned = stb.run(r3_detuned_deck, "XP", SWEEP10, os.path.join(OUT, "r3_detuned_fstart10"))
        m3d = r3_detuned["margins"]
        record("R3 detuned @ fstart=10Hz: low_f_positive_feedback is False (card #476 fix)",
               m3d["low_f_positive_feedback"] is False,
               f"T(10Hz)={r3_detuned['T'][0]!r}, note={m3d['note']!r}")

        # the TUNED point (docs/compopt.md's worked example: r2=27950.5,
        # c2=2.4425e-9, PM=37.998 deg @ 34000 Hz -- a genuinely healthy,
        # well-margined loop of the same near-integrator family) must get
        # the neutral low_f_integrator label, not a warning or note.
        r3_tuned_deck = os.path.join(ROOT, "tests", "fixtures", "stb6_r3_tuned.cir")
        r3_tuned = stb.run(r3_tuned_deck, "XP", SWEEP10, os.path.join(OUT, "r3_tuned_fstart10"))
        check_low_f_integrator("R3 tuned @ fstart=10Hz", r3_tuned["margins"])

        # new fixture (card #476): genuine positive feedback (T(fmin) real
        # and negative, arg=180 deg) but |T(fmin)|<1 there -- the note, not
        # the warning.
        low_gain_deck = os.path.join(ROOT, "tests", "fixtures", "stb6_posfb_low_gain.cir")
        low_gain = stb.run(low_gain_deck, "XP", SWEEP, os.path.join(OUT, "posfb_low_gain"))
        check_low_f_low_gain_note("posfb low-gain fixture", low_gain["margins"])
        t0_lg = low_gain["T"][0]
        record("posfb low-gain fixture: T(fmin) matches closed form (Re<0, |T|<1)",
               abs(t0_lg - (-0.333333)) < 1e-3, f"T(f_min)={t0_lg!r}")

        # margins(): crossings bookkeeping on synthetic data
        fs = [10 ** (k / 10) for k in range(0, 101)]
        flat = [0.5 + 0j for _ in fs]
        m = stbcore.margins(fs, flat)
        record("margins: |T| < 1 everywhere reports no crossing", m["pm"] is None and m["gm"] is None, str(m["pm"]))
        bump = [complex(3 if 1e3 < f < 1e5 else 0.5) * cmath.exp(-1j * math.radians(90)) for f in fs]
        m = stbcore.margins(fs, bump)
        record("margins: two 0 dB crossings both reported", len(m["pm_crossings"]) == 2, f"{len(m['pm_crossings'])} found")

        # margins(): card #476 low-frequency classification, synthetic (a constant
        # T(f) trace so T(fmin) is exactly the value under test, independent of
        # any ngspice run).
        fs2 = [10 ** (k / 10) for k in range(0, 51)]

        def const(z0):
            return [z0] * len(fs2)

        m = stbcore.margins(fs2, const(cmath.rect(2.0, math.radians(170))))
        record("margins synthetic: |arg|=170deg, |T|=2 -> low_f_positive_feedback True (warning)",
               m["low_f_positive_feedback"] is True and m["low_f_positive_feedback_low_gain"] is False
               and m["low_f_integrator"] is False, str(m))

        m = stbcore.margins(fs2, const(cmath.rect(0.5, math.radians(170))))
        record("margins synthetic: |arg|=170deg, |T|=0.5 -> low_f_positive_feedback_low_gain True (note)",
               m["low_f_positive_feedback"] is False and m["low_f_positive_feedback_low_gain"] is True
               and m["low_f_integrator"] is False, str(m))

        m = stbcore.margins(fs2, const(cmath.rect(5.0, math.radians(-90))))
        record("margins synthetic: arg=-90deg -> low_f_integrator True",
               m["low_f_positive_feedback"] is False and m["low_f_positive_feedback_low_gain"] is False
               and m["low_f_integrator"] is True, str(m))

        m = stbcore.margins(fs2, const(cmath.rect(5.0, math.radians(-119))))
        record("margins synthetic: arg=-119deg (inside -90+/-30) -> low_f_integrator True",
               m["low_f_integrator"] is True, str(m))

        m = stbcore.margins(fs2, const(cmath.rect(5.0, math.radians(-121))))
        record("margins synthetic: arg=-121deg (just outside -90+/-30 AND |arg|<150) -> none of the three fire",
               m["low_f_positive_feedback"] is False and m["low_f_positive_feedback_low_gain"] is False
               and m["low_f_integrator"] is False, str(m))

        m = stbcore.margins(fs2, const(5.0 + 0j))
        record("margins synthetic: arg=0deg, ordinary negative feedback -> none of the three fire",
               m["low_f_positive_feedback"] is False and m["low_f_positive_feedback_low_gain"] is False
               and m["low_f_integrator"] is False, str(m))
    except stb.CouldNotRun as e:
        record("run", None, str(e))

    # red-proof for could-not-run: a missing ngspice must never read as a pass
    saved = stb.NGSPICE
    stb.NGSPICE = "/nonexistent/ngspice"
    try:
        stb.run(os.path.join(ROOT, "refs", "r1a_unilateral.cir"), "XP", SWEEP, os.path.join(OUT, "bogus"))
        record("red-proof: missing ngspice gives COULD-NOT-RUN", False, "run returned a result")
    except stb.CouldNotRun as e:
        record("red-proof: missing ngspice gives COULD-NOT-RUN", True, str(e)[:60])
    finally:
        stb.NGSPICE = saved

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    print(f"\n{sum(s == 'pass' for _, s, _ in results)}/{len(results)} pass")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
