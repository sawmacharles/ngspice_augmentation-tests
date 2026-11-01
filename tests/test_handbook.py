#!/usr/bin/env python3
"""Validation for the user handbook (card #507; the handbook HTML page
itself is not published).

Runs docs/examples/handbook/run_all.sh (every example deck/command the
handbook page quotes output from), then checks that the specific numbers
printed in the page actually appear in what that run produced -- within a
small rounding tolerance for floating-point values, exact for integer
counts/labels. This is a red-proof against the page drifting from the
build it claims to document, not a re-derivation of any tool's own
correctness (that is each tool's own docs/*.md and tests/test_*.py's job).

Needs the private-fork build the handbook was written against
(NGSPICE_PREFIX, default matches run_all.sh's own default:
~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152, or the prefix of $NGSPICE) and system ngspice on
PATH (tools/convmeas.py, tools/mc.py, tools/yieldest.py's own fixture
smoke test use it). If run_all.sh itself fails to complete, every check
in this file reports COULD-NOT-RUN, never a silent pass.

Card #543 (s96): examples 32-37 (BSIM4 buck R7 against the converged FRA, forward vs adjoint, the s91 red-proof, the PULSE-overlap trap,
pac_perturbed self-checks, the free-running stability verdict, the 490-state PEX deck, the PWL clock) and an ASan-twin run of two of them.
The ASan build is NGSPICE_ASAN (its bin/ngspice or its prefix; default <prefix>-asan); the s91 red-proofs need the s91 build next to the prefix
(otherwise those checks report COULD-NOT-RUN). The converged FRA curves of R7 are read from
~/Work/data/ngspice_augmentation/fra/curves-r7-20261001/ (not regenerated: 16-24 min per point).

Card #552 (s108): examples 38-44 (.pnoise at the per-step operating point, the MOS track-and-hold, hidden-state continuation / pacreltol / gear on R7,
the free-running R8, ideal-switch events with device charges, XSPICE under .pss and td .pac, the digitally controlled buck R9). Every new example has
a red proof inside the same run (PAC_NOISE_STALE_STATE=1, PAC_NATIVE_HISTORY=1, PSS_NO_HIDDEN=1, PSS_XS_NO_STATE=1, the s96/s105 builds next to the
prefix, a limit-cycling variant). Machine-dependent quantities (iteration counts, wall times) are not checked; bounds sit at >= 3x the measured
spread so the suite passes on both workstations (card #550 rule).

Card #560 (s115): examples 45-50 (.pnoise folded flicker against the closed form of a chopper, PSS memory before/after, speed rounds 16-17 with their
off-switches, free-running .pxf adjoint=1 / .pnoise on R5, OMP_NUM_THREADS vs `set num_threads`, td .pac through an s_xfer compensator).  Red proofs inside
the same run: PNOISE_FLICKER_OLD=1, PSS_KEEP_RUN_DATA=1 and the real s108 build, PAC_FREERUN_PROJ=1, PAC_HID_RED=1, NGSPICE_PSS_AUTOKLU=0 / `.options sparse`;
the cited references (FRA of the R3 + s_xfer loop, the transient-noise ensembles of R5) are constants in the example scripts and the checks.

Usage:
    NGSPICE=<prefix>/bin/ngspice NGSPICE_ASAN=<asan prefix> OMP_NUM_THREADS=2 python3 tests/test_handbook.py     (about 70 min on a loaded machine; examples 45-50 add about 12 min)
    python3 tests/test_handbook.py --prefix ~/Work/apps/ngspice_augmentation/<name>

Exit status: 0 = every check passed, 1 = a check FAILed (page/output
mismatch), 2 = could not run (run_all.sh itself did not complete, or a
log file this script needs is missing).
"""
import argparse
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
HANDBOOK_DIR = os.path.join(ROOT, "docs", "examples", "handbook")
def _default_prefix():
    """NGSPICE=<prefix>/bin/ngspice (or ngspice-aug) selects the build; else s115 (card #560)."""
    env = os.environ.get("NGSPICE") or os.environ.get("NGSPICE_AUG")
    if env:
        return os.path.dirname(os.path.dirname(os.path.abspath(os.path.expanduser(env))))
    return os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152")


DEFAULT_PREFIX = _default_prefix()
OUT = os.path.expanduser("~/Work/data/ngspice_augmentation/handbook/run_all")

results = []
RAN_OK = True


def record(name, state, detail):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN")
    results.append((name, state, detail))
    print(f"[{state:15s}] {name}: {detail}")


def read_log(name):
    path = os.path.join(OUT, f"{name}.log")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return fh.read()


def read_file(*parts):
    path = os.path.join(OUT, *parts)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return fh.read()


def check_contains(name, label, text, needle):
    if text is None:
        record(f"{name}: {label}", "COULD-NOT-RUN", "log file missing")
        return
    ok = needle in text
    record(f"{name}: {label}", "pass" if ok else "FAIL",
           needle if ok else f"expected substring not found: {needle!r}")


def check_float_near(name, label, text, pattern, expect, rel=1e-3, abs_tol=0.0):
    """Find `pattern` (one capture group, a float) in text and compare to expect."""
    if text is None:
        record(f"{name}: {label}", "COULD-NOT-RUN", "log file missing")
        return
    m = re.search(pattern, text)
    if not m:
        record(f"{name}: {label}", "FAIL", f"pattern not found: {pattern!r}")
        return
    got = float(m.group(1))
    tol = max(abs_tol, rel * abs(expect))
    ok = abs(got - expect) <= tol
    record(f"{name}: {label}", "pass" if ok else "FAIL",
           f"got={got!r} expect={expect!r} tol={tol!r}")


def check_converged(name, label, text, tag, it_ref, resid_cap, rho_ref=None, rho_tol=1e-3, verdict="stable", period_ref=None):
    """Cross-machine robust version of the exact 'PSS-SHOOT: CONVERGED iter=.. residual=.. spectral_radius=..' string match (card #550).
    The handbook HTML quotes the digits of the development workstation (AMD Ryzen AI Max+ 395, Linux, Ubuntu glibc); a second workstation (Arch, other libm) gives the same orbit with last-bit differences that the
    shooting residual (a noise-limited quantity just under the convergence tolerance) and the iteration count show. Accepted: the line is
    there, iteration count within +-5 of the development workstation's (second workstation/s103 measured: 8/8, 5/5, free-running 17 vs 20, i.e. spread 0/0/3), residual <= resid_cap (the engine's
    own convergence tolerance, not the digits), spectral radius within rho_tol (measured spread 3e-4, x3), verdict token equal, and for
    the free-running case the period within 2e-5 relative (measured 5.6e-6, x3)."""
    if text is None:
        record(f"{name}: {label}", "COULD-NOT-RUN", "log file missing")
        return
    m = re.search(tag + r": CONVERGED iter=(\d+) residual=(\S+)(?: period=(\S+) fsw=\S+)? spectral_radius=(\S+) (\S+)", text)
    if not m:
        record(f"{name}: {label}", "FAIL", f"no {tag} CONVERGED line")
        return
    it, res, rho, vd = int(m.group(1)), float(m.group(2)), float(m.group(4)), m.group(5)
    ok = abs(it - it_ref) <= 5 and res <= resid_cap and vd == verdict and (rho_ref is None or abs(rho - rho_ref) <= rho_tol)
    det = f"iter {it} (ref {it_ref} +-5), residual {res:.3e} (<= {resid_cap:g}), rho {rho:.6f}" + (f" (ref {rho_ref} +-{rho_tol:g})" if rho_ref is not None else "") + f" {vd}"
    if period_ref is not None:
        per = float(m.group(3))
        ok = ok and abs(per / period_ref - 1) <= 2e-5
        det += f", period {per:.9g} (ref {period_ref:g} +-2e-5 rel)"
    record(f"{name}: {label}", "pass" if ok else "FAIL", det)


def check_complex_near(name, label, text, freq, mag_expect, deg_expect, rel, deg_tol):
    """wrdata complex line `freq re im`: compare |H| (relative) and, if deg_expect is given,
    the phase (absolute degrees) -- card #535."""
    import cmath, math
    if text is None:
        record(f"{name}: {label}", "COULD-NOT-RUN", "log file missing")
        return
    pat = re.escape(f"{freq:.15e}") + r"\s+([-\d.eE+]+)\s+([-\d.eE+]+)"
    m = re.search(pat, text)
    if not m:
        record(f"{name}: {label}", "FAIL", f"line for {freq} Hz not found")
        return
    h = complex(float(m.group(1)), float(m.group(2)))
    mag, deg = abs(h), math.degrees(cmath.phase(h))
    ok = abs(mag - mag_expect) <= rel * mag_expect
    if deg_expect is not None:
        ok = ok and abs(deg - deg_expect) <= deg_tol
    record(f"{name}: {label}", "pass" if ok else "FAIL",
           f"|H|={mag:.6g} (expect {mag_expect:.6g} +-{rel:.1%}), phase={deg:.4f} deg"
           + (f" (expect {deg_expect} +-{deg_tol})" if deg_expect is not None else ""))


def run_all_sh(prefix):
    global RAN_OK
    script = os.path.join(HANDBOOK_DIR, "run_all.sh")
    proc = subprocess.run([script, prefix], cwd=HANDBOOK_DIR,
                           capture_output=True, text=True, timeout=10800)
    if proc.returncode != 0:
        RAN_OK = False
        record("run_all.sh", "COULD-NOT-RUN",
               f"exit {proc.returncode}; tail: {proc.stdout[-1500:]}\n{proc.stderr[-1500:]}")
    else:
        record("run_all.sh", "pass", "completed, exit 0")



# ======================================================================================================= card #543 (s96)
FRA_R7 = os.path.expanduser("~/Work/data/ngspice_augmentation/fra/curves-r7-20261001")


def fra_r7(case):
    """{freq_rounded: (dB, deg, unc_dB, unc_deg)} from the converged transient FRA of the BSIM4 buck (an independent reference:
    tools/fra.py perturb-and-DFT on stock transient simulation, settle x2 / amplitude/2 / reltol variants give the uncertainty)."""
    import csv
    p = os.path.join(FRA_R7, case, "fra.csv")
    if not os.path.exists(p):
        return None
    out = {}
    for r in csv.DictReader(open(p)):
        if r["status"] != "ok":
            continue
        if "combined_db" in r:
            out[round(float(r["freq"]))] = (float(r["combined_db"]), float(r["combined_phase_deg"]), float(r["unc_db"]), float(r["unc_deg"]))
        else:
            out[round(float(r["freq"]))] = (float(r["H_db"]), float(r["H_phase_deg"]), float(r["unc_db"]), float(r["unc_deg"]))
    return out


def rd_cplx(text):
    """wrdata (wr_singlescale, wr_vecnames) of complex vectors: rows (f, [complex per vector])."""
    rows = []
    for ln in (text or "").splitlines()[1:]:
        v = ln.split()
        if len(v) >= 3:
            x = [float(a) for a in v]
            rows.append((x[0], [complex(x[1 + 2 * k], x[2 + 2 * k]) for k in range((len(x) - 1) // 2)]))
    return rows


def vs_fra(name, rows, col, case, bounds, label):
    """compare vector `col` of rows to the FRA table of `case`; |dB| <= max(bounds[0], 3 unc), |deg| <= max(bounds[1], 3 unc)."""
    import cmath
    ref = fra_r7(case)
    if ref is None or not rows:
        record(f"{name}: {label}", "COULD-NOT-RUN", "FRA reference or engine output missing")
        return
    worst_db = worst_deg = 0.0
    worst_ok = True
    n = 0
    det = []
    for f, vec in rows:
        k = round(f)
        if k not in ref:
            continue
        n += 1
        h = vec[col]
        ddb = 20 * math.log10(abs(h)) - ref[k][0]
        ddeg = (math.degrees(cmath.phase(h)) - ref[k][1] + 180.0) % 360.0 - 180.0
        ok = abs(ddb) <= max(bounds[0], 3 * ref[k][2]) and abs(ddeg) <= max(bounds[1], 3 * ref[k][3])
        worst_ok &= ok
        worst_db, worst_deg = max(worst_db, abs(ddb)), max(worst_deg, abs(ddeg))
        if not ok:
            det.append(f"{k} Hz off by {ddb:+.3f} dB / {ddeg:+.2f} deg")
    record(f"{name}: {label}", "pass" if (n >= 8 and worst_ok) else "FAIL",
           f"{n} FRA frequencies, worst {worst_db:.3f} dB / {worst_deg:.2f} deg (bounds max({bounds[0]} dB, 3 unc) / max({bounds[1]} deg, 3 unc))"
           + ("; " + "; ".join(det) if det else ""))


def orbit_rows(path):
    t = read_file(*path) if isinstance(path, tuple) else None
    if t is None:
        return None
    return [[float(x) for x in ln.split()] for ln in t.splitlines()[1:] if ln.strip()]


def orbit_interp(rows, t, c):
    import bisect
    ts = [r[0] for r in rows]
    i = bisect.bisect_left(ts, t)
    if i <= 0:
        return rows[0][c]
    if i >= len(rows):
        return rows[-1][c]
    a = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
    return rows[i - 1][c] + a * (rows[i][c] - rows[i - 1][c])


def orbit_mean(rows, c):
    s = 0.0
    for a, b in zip(rows, rows[1:]):
        s += (b[0] - a[0]) * (a[c] + b[c]) / 2
    return s / (rows[-1][0] - rows[0][0])


def orbit_dev(A, B, c, T=2e-6, n=2000):
    return max(abs(orbit_interp(A, (j + 0.5) * T / n, c) - orbit_interp(B, (j + 0.5) * T / n, c)) for j in range(n))


def check_r7():
    N = "32 R7 BSIM4 buck"
    conv = r"PSS-SHOOT: CONVERGED iter=(\d+) residual=(\S+) spectral_radius=(\S+) (\S+) wall_time"
    # ---- .pss, three tolerances
    base = orbit_rows(("32run", "orbit.txt"))
    for tag, fname, ofile, itmax in (("reltol 1e-4", "32_r7_pss", "orbit.txt", 6), ("reltol 3e-4", "32_r7_pss_reltol3e-4", "orbit_3e-4.txt", 8),
                                      ("reltol 1e-3 (the ngspice default)", "32_r7_pss_reltol1e-3", "orbit_1e-3.txt", 12)):
        t = read_log(fname)
        if t is None or base is None:
            record(f"{N} .pss {tag}", "COULD-NOT-RUN", "log or orbit missing")
            continue
        m = re.search(conv, t)
        if not m:
            record(f"{N} .pss {tag}", "FAIL", "no CONVERGED line (the page says it converges)")
            continue
        it, rho, verdict = int(m.group(1)), float(m.group(3)), m.group(4)
        rows = orbit_rows(("32run", ofile))
        dev = orbit_dev(base, rows, 1)
        ok = it <= itmax and 0.88 < rho < 0.905 and verdict == "stable" and dev < 1e-5
        record(f"{N} .pss {tag}: converges, stable" + ("" if fname == "32_r7_pss" else ", same orbit as reltol 1e-4"), "pass" if ok else "FAIL",
               f"iter {it} (<= {itmax}), rho {rho:.6f} {verdict} (0.88..0.905; the settled .tran envelope decays at 0.899), "
               f"max |vout - vout(1e-4)| {dev:.1e} V (< 1e-5)")
    if base:
        vo = orbit_mean(base, 1)
        vv = [r[1] for r in base]
        il = orbit_mean(base, 2)
        pp = (max(vv) - min(vv)) * 1e3
        # reference: the settled .tran at reltol 1e-5 (refs/README.md R7, 2026-10-01, tests/test_r7.py): Vout 1.45940 V, ripple 70.15 mV,
        # il 1.9464 A; the .pss orbit carries the documented restart offset (+0.36 mV vout, +0.48 mA il)
        ok = abs(vo - 1.45940) < 1.0e-3 and abs(pp - 70.15) < 1.0 and abs(il - 1.9464) < 2e-3
        record(f"{N} .pss orbit vs the settled .tran (cited, tests/test_r7.py)", "pass" if ok else "FAIL",
               f"vout mean {vo:.5f} V (ref 1.45940 +-1 mV), ripple {pp:.2f} mV (ref 70.15 +-1), il mean {il:.4f} A (ref 1.9464 +-2 mA)")
    # ---- .pstb vs FRA
    t = read_file("32run", "pstb.txt")
    lg = read_log("32_r7_pstb")
    rows = rd_cplx(t)
    vs_fra(N + " .pstb", rows, 0, "loop", (0.25, 0.6), "loop gain T vs the converged FRA, every FRA frequency (500 Hz grid)")
    if lg:
        m = re.search(r"pstb: phase margin ([\d.]+) deg at ([\d.eE+]+) Hz", lg)
        ok = bool(m) and abs(float(m.group(1)) - 53.52) < 1.0 and abs(float(m.group(2)) / 26190.8 - 1) < 0.02
        record(f"{N} .pstb: phase margin / crossover vs FRA (53.52 deg, 26.19 kHz)", "pass" if ok else "FAIL", m.group(0) if m else "no PM line")
    else:
        record(f"{N} .pstb: phase margin / crossover vs FRA", "COULD-NOT-RUN", "log missing")
    rows3 = rd_cplx(read_file("32run", "pstb_1e-3.txt"))
    vs_fra(N + " .pstb at the ngspice default reltol 1e-3", rows3, 0, "loop", (0.25, 0.6), "loop gain vs the converged FRA (the page says the 1e-3 run is as good)")
    lg3 = read_log("32_r7_pstb_reltol1e-3")
    m3 = re.search(r"pstb: phase margin ([\d.]+) deg at ([\d.eE+]+) Hz", lg3 or "")
    record(f"{N} .pstb at reltol 1e-3: phase margin / crossover vs FRA (53.52 deg, 26.19 kHz)", "pass" if (m3 and abs(float(m3.group(1)) - 53.52) < 1.0 and abs(float(m3.group(2)) / 26190.8 - 1) < 0.02) else "FAIL", m3.group(0) if m3 else "no PM line")
    # ---- .pac, .pxf
    pac = rd_cplx(read_file("32run", "pac.txt"))
    vs_fra(N + " .pac", pac, 0, "ctrl2out", (0.35, 1.2), "control-to-output vs the converged FRA")
    fwd = rd_cplx(read_file("32run", "pxf_fwd.txt"))
    adj = rd_cplx(read_file("32run", "pxf_adj.txt"))
    if fwd and adj and pac:
        w_adj = max(abs(a - b) / abs(a) for (_, x), (_, y) in zip(fwd, adj) for a, b in zip(x, y))
        record(f"{N} .pxf adjoint=1 equals adjoint=0 (all 3 sources, {len(fwd)} frequencies)", "pass" if w_adj < 1e-9 else "FAIL", f"worst relative difference {w_adj:.1e} (docs: 1e-12)")
        w_pac = max(abs(a[1] - b[0]) / abs(b[0]) for (_, a), (_, b) in zip(fwd, pac))
        record(f"{N} .pxf Vctrl path equals the .pac control-to-output", "pass" if w_pac < 1e-9 else "FAIL", f"worst relative difference {w_pac:.1e}")
    else:
        record(f"{N} .pxf", "COULD-NOT-RUN", "pxf/pac output missing")
    vs_fra(N + " .pxf forward", fwd, 0, "line2out", (0.35, 1.2), "Vin -> vout (line-to-output) vs the converged FRA")
    vs_fra(N + " .pxf forward", fwd, 2, "zout", (0.45, 1.2), "Iinj -> vout (output impedance) vs the converged FRA")
    # ---- which build: generic charge history line present on s96; absent and silently wrong on s91
    lp = read_log("32_r7_pac")
    check_contains(N + " build check", "s96 prints the generic charge-history line (45 state rows)", lp,
                   "PAC-TD: generic charge-based history (card #542): 14 native C/L + model-charge states -> 45 state rows")
    old = read_log("32_r7_pstb_s91")
    t_old = read_file("32old", "pstb.txt")
    if old is None or "build not found" in old:
        record(f"{N} red-proof on the pre-fix build s91", "COULD-NOT-RUN", "s91 build missing")
    else:
        r_old = rd_cplx(t_old)
        mx_old = max((abs(v[0]) for _, v in r_old), default=float("nan"))
        mx_new = max((abs(v[0]) for _, v in rows), default=float("nan"))
        ok = ("generic charge-based history" not in old and mx_old < 1e-6 and mx_new > 1.0 and "no 0 dB crossing" in old)
        record(f"{N} red-proof: s91 .pstb is silently wrong (|T| < 1e-6 everywhere, no crossing, no charge-history line); s96 |T| > 1", "pass" if ok else "FAIL",
               f"s91 max |T| {mx_old:.2e}, s96 max |T| {mx_new:.3f}; s91 says: " + ", ".join(re.findall(r"pstb: [^\n]*", old)[:1]))
    # ---- the ramp-PULSE overlap trap
    res = {}
    for v in ("old", "fixed"):
        for who in ("stock", "fork"):
            tt = read_log(f"32tran_{v}_{who}")
            if tt is None:
                res[(v, who)] = None
            else:
                res[(v, who)] = "aborts" if "Timestep too small" in tt else ("completes" if re.search(r"^vo\s*=", tt, re.M) else "?")
    ok = res == {("old", "stock"): "aborts", ("old", "fork"): "aborts", ("fixed", "stock"): "completes", ("fixed", "fork"): "completes"}
    record(f"{N} ramp PULSE overlap (tr+pw+tf > period): old line aborts .tran at reltol 1e-6 on stock 45.2 AND the fork, fixed line completes on both",
           "pass" if ok else "FAIL", str({f"{k[0]}/{k[1]}": v for k, v in res.items()}))


def check_pac_perturbed_s96():
    N = "33 pac_perturbed"
    t = read_log("33")
    if t is None:
        record(f"{N} ideal-switch deck at exactly 1 kHz", "COULD-NOT-RUN", "log missing")
    else:
        # exact model of this deck (matrix exponentials, event-time perturbation; see example 10): 6.20355 V/V, -1.8874 deg
        m = re.search(r"f_target=1000Hz -> f_actual=1000Hz .*?v\(vout\)=([-\d.]+)dB/([-\d.]+)deg", t)
        if not m:
            record(f"{N} ideal-switch deck at exactly 1 kHz", "FAIL", "no f_actual=1000Hz result line: " + t[-400:])
        else:
            ddb = float(m.group(1)) - 20 * math.log10(6.20355)
            ddeg = float(m.group(2)) - (-1.8874)
            record(f"{N} ideal-switch deck at exactly 1 kHz vs the exact LPTV model", "pass" if abs(ddb) <= 0.02 and abs(ddeg) <= 0.1 else "FAIL",
                   f"{m.group(1)} dB / {m.group(2)} deg vs model {20*math.log10(6.20355):.4f} dB / -1.8874 deg: {ddb:+.4f} dB / {ddeg:+.4f} deg (bounds 0.02 / 0.1)")
    o = read_log("33old")
    check_contains(N + " pre-#536 settings (10 mV, 200 points/cycle, reltol 1e-3)", "point is REFUSED, no value printed", o, "REFUSED")
    if o is not None:
        record(f"{N} pre-#536 settings: no result line with a value for 1000 Hz", "pass" if "v(vout)=" not in o else "FAIL", "no `v(vout)=` line" if "v(vout)=" not in o else "a value was printed")
    # example 34: requested vs measured frequency, PULSE overlap note
    d = read_log("34")
    f = read_log("34fix")
    N = "34 pac_perturbed (clock PULSE overlaps its period)"
    def closed(freq):
        # Rin 1k, Cout 10n, Rload 10k: H = Rload/(Rin+Rload) / (1 + j w Cout (Rin||Rload))
        import cmath
        w = 2 * math.pi * freq
        h = (10e3 / 11e3) / (1 + 1j * w * 10e-9 * (1e3 * 10e3 / 11e3))
        return 20 * math.log10(abs(h)), math.degrees(cmath.phase(h))
    for tag, txt in (("default", d), ("PACP_FIX_PULSE=1 --exact-freq", f)):
        if txt is None:
            record(f"{N} {tag}", "COULD-NOT-RUN", "log missing")
            continue
        m = re.search(r"f_actual=([\d.eE+-]+)Hz .*?v\(out\)=([-\d.]+)dB/([-\d.]+)deg", txt)
        if not m:
            record(f"{N} {tag}", "FAIL", "no result line: " + txt[-300:])
            continue
        fa = float(m.group(1))
        cdb, cdeg = closed(fa)
        ok = abs(float(m.group(2)) - cdb) < 0.01 and abs(float(m.group(3)) - cdeg) < 0.05
        record(f"{N} {tag}: result equals the closed-form RC response at the frequency actually measured ({fa:.6g} Hz)", "pass" if ok else "FAIL",
               f"{m.group(2)} dB / {m.group(3)} deg vs closed form {cdb:.4f} / {cdeg:.4f}")
    check_contains(N, "default run: f_actual differs from the request and the tool says so", d, "was MEASURED AT f_actual=")
    check_contains(N, "default run: PULSE overlap noted, TF NOT changed", d, "NOT changed (set PACP_FIX_PULSE=1")
    check_contains(N, "PACP_FIX_PULSE=1 run: TF shortened", f, "TF shortened to")
    if f is not None:
        record(f"{N} --exact-freq: no `MEASURED AT` warning", "pass" if "MEASURED AT" not in f else "FAIL", "")


def check_freerun_verdict():
    N = "35 free-running verdict"
    # independent truth (card #540, docs/pss-pac.md, cited): the slowest multiplier of the R5 orbit (and of the PEX-loaded R5), measured by perturbing
    # the settled orbit in a plain .tran, is 0.036-0.043 (cycles 3-5, two seeds). The pre-#540 dense engine printed 24.2 UNSTABLE on the 24-state deck
    # and 0.888 on a sibling fixture: far away from it.
    t = read_log("35")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    r = {}
    for m in re.finditer(r"FRVERDICT (\w+): (.*?), iter (\d+), period (\S+), radius (\S+) (\S+), wall (\S+) s(.*)", t):
        r[m.group(1)] = dict(notice=m.group(2), iters=int(m.group(3)), period=float(m.group(4)), rho=float(m.group(5)), verdict=m.group(6), why=m.group(8))
    if set(r) != {"matrixfree_normal", "dense_normal", "matrixfree_forced", "dense_forced"}:
        record(N, "FAIL", "a run is missing or failed: " + t[-400:])
        return
    for tag in ("matrixfree_normal", "dense_normal"):
        x = r[tag]
        ok = x["verdict"].startswith("stable") and 0.025 < x["rho"] < 0.05
        record(f"{N} (R5 + 2 RC trees, 10 states) {tag}: stable, radius within the independent .tran multiplier (0.036-0.043) +-35 %", "pass" if ok else "FAIL",
               f"radius {x['rho']} {x['verdict']}, period {x['period']}, {x['iters']} iterations")
    ok = (r["matrixfree_normal"]["notice"].startswith("notice") and r["dense_normal"]["notice"] == "no notice"
          and r["matrixfree_normal"]["verdict"].endswith("[arnoldi]") and abs(r["matrixfree_normal"]["period"] / r["dense_normal"]["period"] - 1) < 2e-4)
    record(f"{N}: matrix-free is the default (notice, [arnoldi] marker), NGSPICE_PSS_KRYLOV=0 is dense, both find the same period (2e-4)", "pass" if ok else "FAIL",
           f"{r['matrixfree_normal']['notice']} / {r['dense_normal']['notice']}; periods {r['matrixfree_normal']['period']} / {r['dense_normal']['period']}")
    for tag in ("matrixfree_forced", "dense_forced"):
        x = r[tag]
        ok = x["verdict"] == "UNKNOWN(not-computed)" and x["rho"] == -1.0 and "UNSTABLE" not in x["verdict"]
        record(f"{N} {tag}: the s90 finite-difference noise gives UNKNOWN(not-computed) and radius -1, never UNSTABLE", "pass" if ok else "FAIL",
               f"{x['verdict']} radius {x['rho']}{x['why'][:140]}")
    check_contains(N, "dense forced run says why (estimates at h and h/2 disagree)", t, "the estimates at perturbation h and h/2 disagree")
    d = read_log("35d")
    base = read_log("30")
    for tag, tt in (("example 30 default (matrix-free)", base), ("example 30 dense (NGSPICE_PSS_KRYLOV=0)", d)):
        if tt is None:
            record(f"{N} {tag}", "COULD-NOT-RUN", "log missing")
            continue
        m = re.search(r"PSS-FREERUN: CONVERGED iter=(\d+) residual=\S+ period=(\S+) fsw=\S+ spectral_radius=(\S+) (\S+)", tt)
        if not m:
            record(f"{N} {tag}", "FAIL", "no CONVERGED line")
            continue
        ok = m.group(4).startswith("stable") and 0.025 < float(m.group(3)) < 0.05 and abs(float(m.group(2)) / 2.36584e-6 - 1) < 1e-3
        record(f"{N} {tag}: stable, radius 0.036-0.043 +-35 %, period 2.36584 us (1e-3)", "pass" if ok else "FAIL",
               f"iter {m.group(1)}, period {m.group(2)}, radius {m.group(3)} {m.group(4)}")


def check_pwl():
    N = "37 .pss on PWL(...) r=0"
    t = read_log("37")
    o = read_file("37run", "pwl_pss.txt")
    if t is None or o is None:
        record(N, "COULD-NOT-RUN", "log or orbit missing")
        return
    mx = re.search(r"^vmax\s*=\s*(\S+)", t, re.M)
    mn = re.search(r"^vmin\s*=\s*(\S+)", t, re.M)
    v = [float(l.split()[1]) for l in o.splitlines()[1:] if l.strip()]
    ok = ("PSS-SHOOT: CONVERGED" in t and mx and mn and abs(max(v) - float(mx.group(1))) < 1e-4 and abs(min(v) - float(mn.group(1))) < 1e-4)
    record(f"{N}: converges; orbit extremes equal the settled .tran's (1e-4)", "pass" if ok else "FAIL",
           f"orbit {max(v):.6f}/{min(v):.6f}, tran {mx.group(1) if mx else None}/{mn.group(1) if mn else None}")
    old = read_log("37_s91")
    if old is None:
        record(f"{N}: red-proof on s91 (refused)", "COULD-NOT-RUN", "s91 log missing (build absent?)")
    else:
        record(f"{N}: red-proof, the same deck on s91 is refused ('no V/I source with a time function')", "pass" if "no V/I source with a time function" in old else "FAIL", "")


def check_pex_krylov():
    N = "36 PEX-size .pss (490 states)"
    t = read_log("36")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    runs = {}
    for m in re.finditer(r"PEXRUN (\w+) rep(\d+): states (\d+), wall ([\d.]+) s \(PSS-SHOOT ([\d.]+) s\), load (\S+), iter (\d+), period runs (\d+) \((\S+)\), rho (\S+) (\S+), notice lines (\d+)", t):
        runs.setdefault(m.group(1), []).append(dict(states=int(m.group(3)), wall=float(m.group(4)), iters=int(m.group(7)), nruns=int(m.group(8)),
                                                     mode=m.group(9), rho=float(m.group(10)), verdict=m.group(11), notice=int(m.group(12))))
    if set(runs) != {"sparse_dense", "klu_dense", "klu_auto"} or any("FAILED" in l for l in t.splitlines() if l.startswith("PEXRUN")):
        record(N, "FAIL", "an arm is missing or failed: " + t[-400:])
        return
    allr = [x for a in runs.values() for x in a]
    ok = all(x["states"] == 490 and x["verdict"].startswith("stable") and abs(x["rho"] - 0.9803) < 0.002 for x in allr)
    record(f"{N}: every run has 490 states and converges stable at the R3 multiplier 0.9803", "pass" if ok else "FAIL",
           "; ".join(f"{a}: {x['states']} states, rho {x['rho']}, {x['iters']} it" for a, v in runs.items() for x in v[:1]))
    ok = (all(x["mode"] == "dense" and x["notice"] == 0 for k in ("sparse_dense", "klu_dense") for x in runs[k])
          and all(x["mode"] == "krylov" and x["notice"] == 1 for x in runs["klu_auto"]))
    record(f"{N}: krylov=0 stays dense (no notice); with no krylov key the matrix-free path is chosen automatically, notice printed exactly once", "pass" if ok else "FAIL",
           "modes " + str({a: sorted({x['mode'] for x in v}) for a, v in runs.items()}) + ", notices " + str({a: [x['notice'] for x in v] for a, v in runs.items()}))
    nd, na = min(x["nruns"] for x in runs["klu_dense"]), max(x["nruns"] for x in runs["klu_auto"])
    record(f"{N}: matrix-free needs at least 10x fewer period runs than the dense monodromy", "pass" if nd >= 10 * na else "FAIL", f"dense {nd}, matrix-free {na}")
    mn = lambda a: min(x["wall"] for x in runs[a])
    record(f"{N}: matrix-free is faster than dense even on the best-of-n minimum (shared machine)", "pass" if mn("klu_auto") < mn("klu_dense") and mn("klu_auto") < mn("sparse_dense") else "FAIL",
           f"min wall: sparse dense {mn('sparse_dense'):.1f} s, KLU dense {mn('klu_dense'):.1f} s, KLU matrix-free {mn('klu_auto'):.1f} s")
    ok = True
    det = []
    n = 0
    for m in re.finditer(r"PEXORBIT (\w+) vs sparse_dense .*?: vout (\S+), opout (\S+), il (\S+), sw (\S+) \(edge", t):
        n += 1
        v, o, i, w = (float(m.group(k)) for k in (2, 3, 4, 5))
        # smooth vectors agree to the engine's own grid-noise level (README 537: 2.2e-6); the opout and the switch-node edge are looser
        ok &= v < 1e-4 and i < 1e-4 and o < 3e-4 and w < 3e-3
        det.append(f"{m.group(1)}: vout {v:.1e}, opout {o:.1e}, il {i:.1e}, sw edge {w:.1e}")
    record(f"{N}: orbit of KLU-dense and of matrix-free equals the SPARSE-dense orbit (vout, il < 1e-4; opout < 3e-4; switch-node edge < 3e-3 of ripple)", "pass" if (ok and n == 2) else "FAIL", "; ".join(det) or "no PEXORBIT lines")


def check_asan(prefix):
    """Run a few of the new examples under the ASan twin of the build (NGSPICE_ASAN = its bin/ngspice or its prefix; default <prefix>-asan):
    no sanitizer report, and the same converged iteration count as the normal build."""
    N = "ASan twin"
    exe = os.environ.get("NGSPICE_ASAN") or (prefix.rstrip("/") + "-asan")
    if os.path.isdir(exe):
        exe = os.path.join(exe, "bin", "ngspice")
    if not os.path.exists(exe):
        record(N, "COULD-NOT-RUN", f"no ASan build at {exe}")
        return
    scratch = os.path.join(os.path.dirname(OUT), "asan")
    os.makedirs(scratch, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="2", ASAN_OPTIONS="detect_leaks=0")
    cases = [("free-running matrix-free + .pac (example 30)", os.path.join(HANDBOOK_DIR, "30_pac_freerun_r5.cir"), r"PSS-FREERUN: CONVERGED iter=(\d+)", "30"),
             ("PEX-size matrix-free .pss (example 36, klu_auto)", os.path.join(OUT, "36run", "klu_auto_rep0", "deck.cir"), r"PSS-SHOOT: CONVERGED iter=(\d+)", "36")]
    for label, deck, pat, ref in cases:
        if not os.path.exists(deck):
            record(f"{N}: {label}", "COULD-NOT-RUN", f"{deck} missing")
            continue
        d = os.path.join(scratch, ref)
        os.makedirs(d, exist_ok=True)
        try:
            p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", exe, "-b", deck],
                               capture_output=True, text=True, cwd=d, env=env, timeout=1500)
        except (FileNotFoundError, subprocess.TimeoutExpired) as ex:
            record(f"{N}: {label}", "COULD-NOT-RUN", str(ex))
            continue
        o = (p.stdout + p.stderr).replace("\r", "\n")
        a = re.search(pat, o)
        n = re.search(pat, read_log(ref) or read_log("36") or "") if ref == "30" else None
        if ref == "36":
            n = re.search(pat, read_file("36run", "klu_auto_rep0", "stdout.txt") or "")
        ok = ("AddressSanitizer" not in o and "runtime error" not in o and a is not None and n is not None and a.group(1) == n.group(1))
        record(f"{N}: {label}: no sanitizer report, same iteration count as the normal build", "pass" if ok else "FAIL",
               f"rc {p.returncode}; ASan iter {a.group(1) if a else None}, normal iter {n.group(1) if n else None}; sanitizer text: {'yes' if 'AddressSanitizer' in o else 'no'}")


# ======================================================================================================= card #552 (s108)
def th_closed_form_s0(rtrack, c, T, tt):
    """exact one-sided low-frequency PSD of the track-and-hold output (Ornstein-Uhlenbeck restart process), tests/test_pnoise.py th_analytic_S0"""
    kT = 1.380649e-23 * 300.15
    tau = rtrack * c
    a = tt / tau
    Th = T - tt
    q = 1 - math.exp(-a)
    rho = math.exp(-a)
    varA = 2 * tau * (tt - tau * q)
    covAs = tau * q
    var0 = varA + Th * Th + 2 * Th * covAs
    c1 = (tau * q) * covAs + Th * rho * covAs + Th * tau * q + Th * Th * rho
    c2 = ((tau * q) * covAs * rho + Th * rho * rho * covAs + Th * tau * q * rho + Th * Th * rho * rho)
    return 2 * (kT / c) * (var0 + 2 * c1 + 2 * c2) / T


def lines_of(name, prefix):
    t = read_log(name)
    return None if t is None else [l for l in t.splitlines() if l.startswith(prefix)]


def check_pnoise_opoint():
    N = "38 .pnoise operating point"
    t = read_log("38")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    for kind in ("diode", "bjt"):
        m = re.search(rf"PNLTI {kind}: .*worst relative difference (\S+) ", t)
        # measured 6.5e-7 (diode), 8.8e-7 (BJT); bound 3e-6
        record(f"{N}: LTI bias, {kind}: .pnoise dev=1 equals plain .noise, total and every device contribution (3e-6)",
               "pass" if m and float(m.group(1)) <= 3e-6 else "FAIL", m.group(0)[-120:] if m else "no PNLTI line")
    m = re.search(r"PNQS diode: total at 1 kHz: quasi-static \(96 x \.noise\) (\S+), \.pnoise (\S+) \(\S+\), \.pnoise PAC_NOISE_STALE_STATE=1 (\S+) \(x(\S+)\)", t)
    if not m:
        record(f"{N}: diode under a 100 kHz LO", "FAIL", "no PNQS diode line")
    else:
        ref, new, stale = float(m.group(1)), float(m.group(2)), float(m.group(3))
        record(f"{N}: diode under a 100 kHz LO: .pnoise equals the quasi-static average of 96 x .noise (1e-4; folded-flicker reference, card #556)", "pass" if abs(new / ref - 1) <= 1e-4 else "FAIL",
               f"{new:.5e} vs {ref:.5e} ({new / ref - 1:+.1e}); measured 4e-7")
        record(f"{N}: red proof PAC_NOISE_STALE_STATE=1 (the pre-#544 operating point) is > 2x too high on the diode", "pass" if stale / ref > 2 else "FAIL", f"x{stale / ref:.2f} (measured 3.01 against the folded-flicker reference of card #556; 4.31 against the plain time average it replaced)")
    m = re.search(r"PNQS diode on the s96 build \(before patch 0097\): total (\S+) \(x(\S+) of quasi-static\)", t)
    if not m:
        record(f"{N}: the real pre-#544 build (s96) on the diode", "COULD-NOT-RUN", "s96 build missing or its line absent")
    else:
        record(f"{N}: the s96 build (no switch) is > 3x too high on the diode", "pass" if float(m.group(2)) > 3 else "FAIL", f"x{m.group(2)}")
    m = re.search(r"PNQS bjt: total at 1 kHz: quasi-static \(96 x \.noise\) (\S+), \.pnoise (\S+) \(\S+\), \.pnoise PAC_NOISE_STALE_STATE=1 (\S+) \(x(\S+)\)", t)
    if not m:
        record(f"{N}: BJT under a 100 kHz LO", "FAIL", "no PNQS bjt line")
    else:
        ref, new, stale = float(m.group(1)), float(m.group(2)), float(m.group(3))
        record(f"{N}: BJT under a 100 kHz LO: .pnoise equals the quasi-static average (1e-3; folded-flicker reference, card #556)", "pass" if abs(new / ref - 1) <= 1e-3 else "FAIL", f"{new:.5e} vs {ref:.5e} ({new / ref - 1:+.1e}); measured 4.5e-4")
        record(f"{N}: red proof PAC_NOISE_STALE_STATE=1 moves the BJT total by more than 10 %", "pass" if abs(stale / ref - 1) > 0.1 else "FAIL", f"x{stale / ref:.2f} (measured 0.80)")


def check_th_mos():
    N = "39 .pnoise MOS track-and-hold"
    T = 1e-4
    got = {}
    for k, ron in (("mos1", 455.0), ("bsim4", 120.0)):
        for v in ("normal", "native"):
            txt = read_file(f"39_{k}_{v}", f"pnoise_th_{k}.txt")
            if txt is None:
                record(f"{N} {k} {v}", "COULD-NOT-RUN", "output missing")
                continue
            rows = [[float(x) for x in l.split()] for l in txt.splitlines()[1:] if l.strip()]
            S0 = th_closed_form_s0(50e3 + ron, 50e-12, T, 0.2 * T)
            ds = {}
            for f, val in rows:
                x = math.pi * f * T
                ds[round(f)] = 20 * math.log10(val / (math.sqrt(S0) * math.sin(x) / x))
            got[(k, v)] = ds
    for k in ("mos1", "bsim4"):
        if (k, "normal") in got:
            ds = got[(k, "normal")]
            w = max(abs(d) for d in ds.values())
            record(f"{N} {k}: folded kT/C closed form at 100 Hz and 1 kHz (0.1 dB)", "pass" if w <= 0.1 and len(ds) == 2 else "FAIL",
                   ", ".join(f"{f} Hz {d:+.3f} dB" for f, d in sorted(ds.items())) + " (measured: mos1 -0.039/-0.000, bsim4 +0.006/+0.045)")
        if (k, "native") in got:
            d = abs(got[(k, "native")][1000])
            record(f"{N} {k}: red proof PAC_NATIVE_HISTORY=1 is more than 0.5 dB off at 1 kHz", "pass" if d > 0.5 else "FAIL", f"{d:.2f} dB (measured 0.97 / 3.45)")
    naive = 20 * math.log10(math.sqrt(2 * 1.380649e-23 * 300.15 / 50e-12 * T) / math.sqrt(th_closed_form_s0(50e3 + 455.0, 50e-12, T, 0.2 * T)))
    record(f"{N}: the check discriminates (the naive kT/C zero-order hold is > 1 dB away)", "pass" if abs(naive) > 1.0 else "FAIL", f"{naive:+.2f} dB")


def tran_rows(d, fn, ncol):
    t = read_file(d, fn)
    return None if t is None else [[float(x) for x in l.split()][:1 + ncol] for l in t.splitlines()[1:] if l.strip()]


def check_r7_552():
    N = "40 R7 after #546/#549"
    tr = tran_rows("40run", "tran.txt", 2)
    hid, nohid = tran_rows("40run", "orbit_hid.txt", 2), tran_rows("40run", "orbit_nohid.txt", 2)
    if not (tr and hid and nohid):
        record(N, "COULD-NOT-RUN", "tran / orbit files missing")
        return
    ref = [[r[0] - 298e-6] + r[1:] for r in tr if 298e-6 - 1e-12 <= r[0] <= 300e-6 + 1e-12]
    prev = [[r[0] - 296e-6] + r[1:] for r in tr if 296e-6 - 1e-12 <= r[0] <= 298e-6 + 1e-12]
    rip = max(r[1] for r in ref) - min(r[1] for r in ref)
    per = max(abs(orbit_interp(prev, r[0], 1) - r[1]) for r in ref) / rip
    record(f"{N}: the settled .tran (reltol 1e-5) is periodic to 1e-5 of the ripple (the reference is converged)", "pass" if per < 1e-5 else "FAIL", f"{per:.1e} of ripple {rip * 1e3:.2f} mV")
    res = {}
    for tag, orb in (("hid", hid), ("nohid", nohid)):
        dm = (orbit_mean(orb, 1) - orbit_mean(ref, 1)) * 1e3
        pt = max(abs(orbit_interp(orb, (j + 0.5) * 2e-6 / 2000, 1) - orbit_interp(ref, (j + 0.5) * 2e-6 / 2000, 1)) for j in range(2000) if (j + 0.5) * 2e-6 / 2000 >= 2e-9) / rip * 100
        res[tag] = (dm, pt)
    record(f"{N}: .pss orbit with the hidden-state continuation equals the settled .tran (mean 0.02 mV, pointwise 0.03 % of ripple)",
           "pass" if abs(res["hid"][0]) < 0.02 and res["hid"][1] < 0.03 else "FAIL", f"mean offset {res['hid'][0]:+.4f} mV, max pointwise {res['hid'][1]:.4f} % (measured -0.004 mV, 0.015 %)")
    record(f"{N}: red proof PSS_NO_HIDDEN=1 (restart from the native states) is off by > 0.2 mV and > 0.3 % of ripple",
           "pass" if res["nohid"][0] > 0.2 and res["nohid"][1] > 0.3 else "FAIL", f"mean offset {res['nohid'][0]:+.4f} mV, max pointwise {res['nohid'][1]:.4f} % (measured +0.362 mV, 0.567 %)")
    fra = fra_r7("loop")
    if not fra:
        record(f"{N}: pacreltol", "COULD-NOT-RUN", "FRA reference missing")
        return
    import cmath

    def dd(tag):
        rows = rd_cplx(read_file("40run", f"pstb_{tag}.txt"))
        out = {}
        for f, v in rows:
            k = round(f)
            if k in fra:
                out[k] = (20 * math.log10(abs(v[0])) - fra[k][0], (math.degrees(cmath.phase(v[0])) - fra[k][1] + 180) % 360 - 180)
        return out, rows
    try:
        dflt, _ = dd("default")
        legacy, _ = dd("legacy")
        pacr, _ = dd("pacreltol")
        gear, rg = dd("gear")
    except (TypeError, AttributeError):
        record(f"{N}: pacreltol", "COULD-NOT-RUN", "pstb output missing")
        return
    sel = (20000, 40000, 80000, 200000)
    if not all(f in d for d in (dflt, legacy, pacr, gear) for f in sel):
        record(f"{N}: pacreltol", "COULD-NOT-RUN", "a .pstb output is missing or lacks the FRA frequencies")
        return
    sel = (20000, 40000, 80000, 200000)
    record(f"{N}: legacy recording tolerance (PAC_RECORD_RELTOL_AUTO=0, the pre-#625 default) leaves the loop gain +0.09..0.13 dB above the converged FRA at 20-200 kHz",
           "pass" if all(0.07 < legacy[f][0] < 0.16 for f in sel) else "FAIL", " ".join(f"{f // 1000}k {legacy[f][0]:+.3f}" for f in sel) + " (measured +0.115..+0.122)")
    record(f"{N}: the default (card #625: recording reltol 1e-8) is within 0.05 dB of the FRA at 20-200 kHz (red: the legacy run above is > 0.07)",
           "pass" if all(abs(dflt[f][0]) < 0.05 for f in sel) else "FAIL", " ".join(f"{f // 1000}k {dflt[f][0]:+.3f}" for f in sel) + " (measured on s136 +0.014..+0.021)")
    record(f"{N}: pacreltol=1e-7 brings it to within 0.05 dB of the FRA at 20-200 kHz (red: the legacy default is > 0.07)", "pass" if all(abs(pacr[f][0]) < 0.05 for f in sel) else "FAIL",
           " ".join(f"{f // 1000}k {pacr[f][0]:+.3f}" for f in sel) + " (measured +0.014..+0.021)")
    record(f"{N}: pacreltol=1e-7 phase within 0.5 deg of the FRA at 20-200 kHz", "pass" if all(abs(pacr[f][1]) < 0.5 for f in sel) else "FAIL", " ".join(f"{f // 1000}k {pacr[f][1]:+.2f}" for f in sel))
    tr_rows = rd_cplx(read_file("40run", "pstb_pacreltol.txt"))
    wd = wp = 0.0
    for (f, a), (_, b) in zip(tr_rows, rg):
        wd = max(wd, abs(20 * math.log10(abs(a[0]) / abs(b[0]))))
        wp = max(wp, abs((math.degrees(cmath.phase(a[0])) - math.degrees(cmath.phase(b[0])) + 180) % 360 - 180))
    record(f"{N}: method=gear equals method=trap (both pacreltol=1e-7) on all {len(rg)} frequencies (0.02 dB / 0.05 deg)", "pass" if wd <= 0.02 and wp <= 0.05 and len(rg) == 400 else "FAIL",
           f"worst {wd:.4f} dB / {wp:.4f} deg (measured 0.0047 / 0.0016)")
    check_contains(N, "method=gear .pnoise is refused with a message", read_log("40"), "method=gear is not supported by the adjoint/noise recursion")


def check_r8():
    N = "41 R8 free-running"
    lg = {k: read_file("41run", f"r8_pss_{k}.log") for k in ("hid", "nohid")}
    tr = tran_rows("41run", "tran.txt", 4)
    if not tr or not all(lg.values()):
        record(N, "COULD-NOT-RUN", "tran / pss logs missing")
        return
    ts = [r[0] for r in tr]
    cx = [ts[i - 1] + (1.7055 - tr[i - 1][4]) * (ts[i] - ts[i - 1]) / (tr[i][4] - tr[i - 1][4]) for i in range(1, len(tr)) if tr[i - 1][4] < 1.7055 <= tr[i][4]]
    per = [b - a for a, b in zip(cx, cx[1:])]
    tmean = sum(per) / len(per)
    record(f"{N}: the settled .tran period (reltol 1e-6) is steady to 2e-5 over {len(per)} periods", "pass" if len(per) >= 15 and (max(per) - min(per)) / tmean < 2e-5 else "FAIL",
           f"mean {tmean * 1e6:.6f} us, spread {(max(per) - min(per)) / tmean:.1e}")
    p = {}
    for k in ("hid", "nohid"):
        m = re.search(r"PSS-FREERUN: CONVERGED iter=(\d+) residual=\S+ period=(\S+) fsw=\S+ spectral_radius=(\S+) (\S+)", lg[k])
        p[k] = float(m.group(2)) / tmean - 1 if m else None
        if not m:
            record(f"{N} {k}: converges", "FAIL", "no CONVERGED line")
    if p["hid"] is not None:
        record(f"{N}: hidden-state continuation: .pss period equals the settled .tran period at equal tolerance (1e-5 relative)", "pass" if abs(p["hid"]) <= 1e-5 else "FAIL", f"{p['hid']:+.1e} (measured +2.0e-6)")
    if p["nohid"] is not None:
        record(f"{N}: red proof PSS_NO_HIDDEN=1: the restarted map's period is off by more than 5e-4 relative", "pass" if abs(p["nohid"]) > 5e-4 else "FAIL", f"{p['nohid']:+.1e} (measured -9.8e-4)")
    t = read_log("41")
    for k, lim, red in (("hid", (0.05, 0.05, 0.1), False), ("nohid", (0.5, 0.5, 1.0), True)):
        m = re.search(rf"R8ORB {k}: .*?vout ripple \S+ mV, mean offset \S+ %, max pointwise (\S+) % of ripple; il ripple \S+ A, mean offset \S+ %, max pointwise (\S+) % of ripple; cx ripple \S+ mV, mean offset \S+ %, max pointwise (\S+) %", t or "")
        if not m:
            record(f"{N} orbit {k}", "FAIL", "no R8ORB line")
            continue
        v = [float(m.group(i)) for i in (1, 2, 3)]
        ok = all(x > l for x, l in zip(v, lim)) if red else all(x <= l for x, l in zip(v, lim))
        record(f"{N}: orbit vs the settled .tran, max pointwise vout/il/cx % of ripple " + ("[RED PROOF PSS_NO_HIDDEN=1: each above " if red else "(each below ") + f"{lim})",
               "pass" if ok else "FAIL", f"{v} (measured hid 0.012/0.013/0.034, nohid 0.816/0.896/2.237)")
    check_contains(N, ".pac runs on the generic charge-history path", t, "PAC-TD: generic charge-based history")
    ref = {20000: (7.8977, -0.2728), 50000: (7.9969, -0.8614), 80000: (8.1740, -1.4905)}
    got = {int(m.group(1)): (float(m.group(2)), float(m.group(3))) for m in re.finditer(r"R8PAC f=(\d+) Hz: ([+-][\d.]+) dB / ([+-][\d.]+) deg", t or "")}
    ok = set(got) == set(ref) and all(abs(got[f][0] - ref[f][0]) <= 0.05 and abs(got[f][1] - ref[f][1]) <= 0.2 for f in ref)
    record(f"{N}: .pac control-to-output at 20/50/80 kHz vs pac_perturbed --freerun (cited, not rerun; 0.05 dB / 0.2 deg)", "pass" if ok else "FAIL",
           "; ".join(f"{f // 1000}k {got[f][0]:+.4f} dB/{got[f][1]:+.3f} deg vs {ref[f][0]:+.4f}/{ref[f][1]:+.4f}" for f in sorted(got)) + " (measured worst 0.0044 dB / 0.052 deg)")


def check_salt_charge_552():
    N = "42 ideal switch + device charges"
    t = read_log("42")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    for f in (20000, 100000):
        for tag in ("default", "PSS_NO_HIDDEN=1", "PAC_NATIVE_HISTORY=1"):
            m = re.search(rf"SALTID {f} Hz {re.escape(tag)}: .*? = (\S+); R6ob \S+ dB / \S+ deg, F x pac_perturbed\(R6o\) \S+ / \S+: ([+-][\d.]+) dB / ([+-][\d.]+) deg", t)
            if not m:
                record(f"{N} {f} Hz {tag}", "FAIL", "no SALTID line")
                continue
            e, db, dg = float(m.group(1)), float(m.group(2)), float(m.group(3))
            if tag == "default":
                ok, lab = e <= 2e-3 and abs(db) <= 0.01 and abs(dg) <= 0.05, "identity H(R6ob) = F x H(R6o) (2e-3) and F x pac_perturbed(R6o) (0.01 dB / 0.05 deg)"
            elif tag == "PSS_NO_HIDDEN=1":
                ok, lab = e > 2e-3 and abs(db) >= 0.02, "red proof: restart without hidden state breaks the identity (> 2e-3, > 0.02 dB)"
            else:
                ok, lab = abs(db) > 3, "red proof: native-only history is > 3 dB off"
            record(f"{N} {f // 1000} kHz {tag}: {lab}", "pass" if ok else "FAIL", f"identity {e:.2e}, {db:+.3f} dB / {dg:+.2f} deg (measured: default 2.9e-4 / -0.001; no-hidden 4.8e-3 / +0.043; native -27 dB)")
    # s118 (card #575): the charged diode at a node without capacitor runs (fast-node path) and is checked against the cited perturbed PSS;
    # with an RC snubber at the node it was refused (s118-s143) and runs since s144 (card #713), checked against its own perturbed PSS
    m1 = re.search(r"SALTDIODE charged diode \(CJO=50p TT=1n\): runs, fast-node path \(inductor-routed jump: yes\): \S+ dB / \S+ deg at 20 kHz, pac_perturbed \S+ / \S+: (\S+) dB / (\S+) deg", t)
    record(f"{N}: a charged diode at a switch node without capacitor runs on the fast-node path, within 0.1 dB / 0.5 deg of the cited pac_perturbed (s118; refused up to s117)",
           "pass" if m1 and abs(float(m1.group(1))) <= 0.1 and abs(float(m1.group(2))) <= 0.5 else "FAIL", (f"{m1.group(1)} dB / {m1.group(2)} deg" if m1 else "no line"))
    m1b = re.search(r"SALTDIODE charged diode with an RC snubber: runs, fast-node path \(inductor-routed jump: yes\): \S+ dB / \S+ deg at 20 kHz, pac_perturbed \S+ / \S+: (\S+) dB / (\S+) deg", t)
    record(f"{N}: a charged diode at a switch node with an RC snubber and no node capacitor runs, within 0.02 dB / 0.1 deg of its pac_perturbed (s144; refused s118-s143)",
           "pass" if m1b and abs(float(m1b.group(1))) <= 0.02 and abs(float(m1b.group(2))) <= 0.1 else "FAIL", (f"{m1b.group(1)} dB / {m1b.group(2)} deg (measured s144: -0.004 / -0.01)" if m1b else "no line"))
    m2 = re.search(r"SALTDIODE chargeless diode: (.*)", t)
    record(f"{N}: a chargeless diode beside the switch runs without refusal", "pass" if m2 and m2.group(1).startswith("no refusal") else "FAIL", (m2.group(1)[:100] if m2 else "no line"))


def check_xspice_552():
    N = "43 XSPICE"
    t = read_log("43")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return

    def val(pat):
        m = re.search(pat, t)
        return m.groups() if m else None
    g = val(r"XSPAC s(?!105)\d+: .*worst relative difference (\S+) ")
    record(f"{N}: td .pac of gain(out_offset=0.1) equals plain .ac (1e-9)", "pass" if g and float(g[0]) <= 1e-9 else "FAIL", f"{g[0] if g else 'no line'} (measured 6e-13)")
    g = val(r"XSPAC s105: .*worst relative difference (\S+) ")
    record(f"{N}: the s105 build (before) is off by more than 1 % (the model's DC rhs in the small-signal source)", "pass" if g and float(g[0]) > 0.01 else ("COULD-NOT-RUN" if not g else "FAIL"),
           f"{g[0] if g else 's105 build missing'} (measured 4.0e-2)")
    g = val(r"XSPSS s(?!105)\d+ s_xfer: orbit vs settled tran: vout (\S+) %, o1 (\S+) %")
    record(f"{N}: .pss of s_xfer equals the settled .tran (1.5 % of ripple; the model's own first-order integrator limits it)", "pass" if g and max(float(x) for x in g) <= 1.5 else "FAIL", f"{g} (measured 0.62)")
    g = val(r"XSPSS s105 s_xfer: orbit vs settled tran: vout (\S+) %, o1 (\S+) %")
    record(f"{N}: the s105 build converges to a WRONG s_xfer orbit (> 3 % of ripple)", "pass" if g and min(float(x) for x in g) > 3 else ("COULD-NOT-RUN" if not g else "FAIL"), f"{g} (measured 4.9)")
    g = val(r"XSPSS s(?!105)\d+ s_xfer PSS_XS_NO_STATE=1: orbit vs settled tran: vout (\S+) %, o1 (\S+) %")
    record(f"{N}: red proof PSS_XS_NO_STATE=1 on s_xfer is > 3 % off", "pass" if g and min(float(x) for x in g) > 3 else "FAIL", f"{g} (measured 4.9)")
    g = val(r"XSPSS s(?!105)\d+ d_tff divider .*?: orbit vs settled tran (\S+) % of ripple; .*?digital node mismatches (\d+), analog node drift (\S+)")
    record(f"{N}: .pss of the d_tff divider by 4 (adc_bridge, d_tff, dac_bridge) equals the settled .tran (0.4 %), no digital mismatch", "pass" if g and float(g[0]) <= 0.4 and g[1] == "0" else "FAIL",
           f"{g} (measured 0.110 %, 0 mismatches)")
    g = re.search(r"XSPSS s105 d_tff divider[^\n]*", t)
    record(f"{N}: the s105 build crashes (segfault) on the same digital deck", "pass" if g and "segfault" in g.group(0) else ("COULD-NOT-RUN" if not g else "FAIL"), (g.group(0)[-60:] if g else "s105 build missing"))
    g = val(r"XSPSS s(?!105)\d+ d_dff sample-and-hold .*?default (\S+) % of ripple; PSS_XS_NO_STATE=1 (\S+) %")
    record(f"{N}: d_dff sample-and-hold: .pss equals the settled .tran (0.1 %); red proof PSS_XS_NO_STATE=1 is > 1 %", "pass" if g and float(g[0]) <= 0.1 and float(g[1]) > 1 else "FAIL", f"{g} (measured 0.023 / 8.55)")
    m = re.search(r"XSTD \.pac on s_xfer: (.*)", t)
    record(f"{N}: td .pac now runs THROUGH s_xfer (card #557: hidden states as Floquet states)", "pass" if m and m.group(1).startswith("runs, 2 hidden states") else "FAIL", (m.group(1)[:110] if m else "no line"))
    m = re.search(r"XSTD \.pac on the digital deck: (.*)", t)
    record(f"{N}: td .pac runs through the digital deck (card #714: edge shifts)", "pass" if m and m.group(1).startswith("runs") else "FAIL", (m.group(1)[:110] if m else "no line"))
    for key, pat in ((".pac on oneshot", "keeps integration state"), (".pac on the digital deck with PAC_DIG_REFUSE=1", "'adc_bridge'"), (".pss on an unlisted model (zener)", "'zener'"),
                     ("free-running (trigger=) with an XSPICE model", "free-running (trigger=)")):
        m = re.search(r"XSREF " + re.escape(key) + r": (.*)", t)
        record(f"{N}: refused: {key}", "pass" if m and "REFUSED" in m.group(1) and pat in m.group(1) else "FAIL", (m.group(1)[:110] if m else "no line"))
    m = re.search(r"XSREF d_tff toggling at the clock period, \.pss at 1 MHz: (.*)", t)
    record(f"{N}: a digital state not periodic with the analysis period (d_tff toggling at 1 MHz) is refused", "pass" if m and m.group(1).startswith("REFUSED") and "not periodic" in m.group(1) else "FAIL", (m.group(1)[:110] if m else "no line"))
    m = re.search(r"XSREF d_tff toggling at the clock period, \.pss at 0\.5 MHz: converges, orbit vs settled tran (\S+) %", t)
    record(f"{N}: the same deck at twice the period converges and equals the settled .tran (0.4 %)", "pass" if m and float(m.group(1)) <= 0.4 else "FAIL", m.group(0)[-60:] if m else "no line")


def check_r9_552():
    N = "44 R9 digital buck"
    t = read_log("44")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    m = re.search(r"R9TRAN: .*DPWM code \[(\d+)\] over", t)
    record(f"{N}: settled .tran: constant DPWM code 556 from 6 ms (no limit cycle)", "pass" if m and m.group(1) == "556" else "FAIL", m.group(0) if m else "no line")
    m = re.search(r"R9PSS: iter (\d+), residual (\S+), radius (\S+) (\S+), \S+s; orbit vs the last tran period \(58.8\.\.60 ms\): ripple (\S+) mV, max pointwise (\S+) mV = (\S+) % of ripple, mean offset (\S+) mV; DPWM code \[(\d+)\]; .*?digital node mismatches (\d+)", t)
    if not m:
        record(f"{N}: .pss at T = 1.2 ms", "FAIL", "no converged R9PSS line: " + "; ".join(l for l in t.splitlines() if l.startswith("R9PSS"))[:200])
    else:
        rip, pct, mo = float(m.group(5)), float(m.group(7)), float(m.group(8))
        ok = m.group(4) == "stable" and pct <= 1.0 and abs(mo) <= 0.1 * rip and m.group(9) == "556" and m.group(10) == "0"
        record(f"{N}: .pss of the whole digital buck (600 switching cycles per period) equals the settled .tran (1 % of ripple), code 556, no digital mismatch", "pass" if ok else "FAIL",
               f"{pct:.2f} % of the {rip:.3f} mV ripple, mean offset {mo:+.4f} mV, radius {m.group(3)} {m.group(4)}, mismatches {m.group(10)} (measured 0.33 %, +0.005 mV)")
    m = re.search(r"R9PAC: (.*?)( \(script wall)", t)
    record(f"{N}: .pstb is refused, naming the XSPICE models", "pass" if m and "REFUSED" in m.group(1) and re.search(r"adc_bridge|d_tff|dac_bridge|d_dff", m.group(1)) else "FAIL", (m.group(1)[:120] if m else "no line"))
    m = re.search(r"R97BIT tran: .*differ by up to (\d+) mV; DPWM codes visited \[([\d, ]+)\]", t)
    record(f"{N}: 7-bit variant: the settled .tran hunts (consecutive 1.2 ms windows differ by > 50 mV, at least 2 codes)", "pass" if m and int(m.group(1)) > 50 and len(m.group(2).split(",")) >= 2 else "FAIL",
           m.group(0)[-90:] if m else "no line")
    m = re.search(r"R97BIT pss at T: (.*?)( \(script wall)", t)
    record(f"{N}: 7-bit variant: .pss at T is refused (no CONVERGED)", "pass" if m and not m.group(1).startswith("CONVERGED") and re.search(r"REFUSED|FAILED|not periodic", m.group(1)) else "FAIL", (m.group(1)[:140] if m else "no line"))


# ======================================================================================================= card #560 (s115)
def check_fold_560():
    N = "45 .pnoise folded flicker"
    t = read_log("45")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    rows = [r for r in re.findall(r"FOLD\s+(\d+) Hz\s+flicker engine (\S+) exact (\S+) \((\S+)\)\s+old (\S+) \(([^)]*)\)\s+thermal engine (\S+) exact (\S+) \((\S+)\)", t)]
    if len(rows) != 4:
        record(N, "FAIL", f"{len(rows)} FOLD lines (expected 4)")
        return
    worst = max(abs(float(r[3])) for r in rows)
    worst_th = max(abs(float(r[8])) for r in rows)
    record(f"{N}: chopped resistor 1/f at 1 / 30 / 90 / 150 kHz equals the closed form (1e-5; thermal too)", "pass" if worst <= 1e-5 and worst_th <= 1e-5 else "FAIL",
           f"worst flicker {worst:.1e}, worst thermal {worst_th:.1e} (measured 9e-8 / 2.2e-7)")
    old = {int(r[0]): (float(r[4]) / float(r[2]) - 1) for r in rows}
    ok = old[1000] > 5 and old[30000] > 0.3 and old[90000] < -0.3 and old[150000] < -0.1 and min(abs(v) for v in old.values()) > 0.1
    record(f"{N}: red proof PNOISE_FLICKER_OLD=1 is far off at every frequency (x10 / +80 % / -52 % / -23 %)", "pass" if ok else "FAIL",
           ", ".join(f"{f // 1000} kHz {v:+.2f}" for f, v in sorted(old.items())))


def check_memory_560():
    N = "46 PSS memory"
    t = read_log("46")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    runs = {}
    for m in re.finditer(r"MEMRUN (.*?), (before \(s108\)|after \(s115\)|after, PSS_KEEP_RUN_DATA=1): peak RSS (\d+) MB, wall (\S+) s, period runs (\d+), rho (\S+ \S+), orbit md5 (\w+)", t):
        runs[(m.group(1), m.group(2).split()[0])] = (int(m.group(3)), float(m.group(4)), int(m.group(5)), m.group(6), m.group(7))
    for deck, bound_after, min_ratio in (("490 states", 60, 5.0), ("1150 states", 90, 8.0)):
        a, k, b = runs.get((deck, "after")), runs.get((deck, "after,")), runs.get((deck, "before"))
        if a is None or k is None:
            record(f"{N}: {deck}", "FAIL", "MEMRUN lines missing")
            continue
        record(f"{N}: {deck} without .save: peak RSS after (s115) <= {bound_after} MB and PSS_KEEP_RUN_DATA=1 (the old retention) is >= {min_ratio:g}x larger (red proof)",
               "pass" if a[0] <= bound_after and k[0] >= min_ratio * a[0] else "FAIL", f"after {a[0]} MB, switch {k[0]} MB (x{k[0] / max(a[0], 1):.1f}); measured 17/243 and 26/622 (s115 / switch, 490 / 1150 states)")
        if b is None:
            record(f"{N}: {deck}: the real s108 build is >= {min_ratio:g}x larger", "COULD-NOT-RUN", "s108 build missing")
        else:
            record(f"{N}: {deck}: the real s108 build (no patch 0109) is >= {min_ratio:g}x larger", "pass" if b[0] >= min_ratio * a[0] else "FAIL", f"s108 {b[0]} MB, s115 {a[0]} MB")
    m = re.search(r"MEMSAME (\S+ states): orbit files byte-identical across (\d+) arms: (\w+)", t)
    same = re.findall(r"MEMSAME .*?: orbit files byte-identical across \d+ arms: (\w+)", t)
    record(f"{N}: the wrdata orbit is byte-identical before / after / with the switch (all three decks)", "pass" if len(same) == 3 and all(x == "True" for x in same) else "FAIL", str(same))
    sv = runs.get(("490 states, .save", "after"))
    record(f"{N}: with .save the peak is small either way (<= 40 MB after)", "pass" if sv and sv[0] <= 40 else "FAIL", str(sv))
    a, b = runs.get(("490 states", "after")), runs.get(("490 states", "before"))
    if a and b:
        record(f"{N}: period runs 45 (s115) against 47 (s108): the Arnoldi stop of round 16, same radius", "pass" if a[2] < b[2] and a[3] == b[3] else "FAIL", f"{a[2]} vs {b[2]}, rho {a[3]} / {b[3]}")
    m8 = re.search(r"MEMPLOTS s108: (.*)", t)
    m15 = re.search(r"MEMPLOTS s115: (.*)", t)
    record(f"{N}: after .pss method=shoot only the emitted trajectory is listed (2 plots: const + pss<N>); s108 listed >= 10", "pass" if m15 and len(m15.group(1).split()) == 2 and (not m8 or len(m8.group(1).split()) >= 10) else "FAIL",
           f"s115: {m15.group(1) if m15 else None}; s108: {(m8.group(1)[:60] + '...') if m8 else 'build missing'}")


def check_speed_560():
    N = "47 speed rounds 16-17"
    t = read_log("47")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    for tag, nodes, wall_bound in (("b2d4", None, None), ("b3d3", 136, 1.0), ("b3d4", 379, 0.9)):
        m = re.search(r"AKLU %s \((\d+) states, ([^)]*)\): notice in default (\d+), AUTOKLU=0 (\d+), options sparse (\d+); median wall default (\S+) s, AUTOKLU=0 (\S+) s \(([^)]*)\), options sparse (\S+) s; output byte-identical to AUTOKLU=0: (\w+); orbit deviation default vs SPARSE (\S+) of ripple" % tag, t)
        if not m:
            record(f"{N}: automatic KLU {tag}", "FAIL", "no AKLU line")
            continue
        nd, noff, nsp = int(m.group(3)), int(m.group(4)), int(m.group(5))
        wd, wo = float(m.group(6)), float(m.group(7))
        if nodes is None:
            ok = nd == 0 and noff == 0 and nsp == 0 and m.group(10) == "True" and "below 128" in m.group(2)
            record(f"{N}: automatic KLU, {tag} ({m.group(1)} states, 109 nodes, below 128): no notice, output byte-identical to NGSPICE_PSS_AUTOKLU=0", "pass" if ok else "FAIL", m.group(0)[:200])
        else:
            ok = (nd >= 1 and noff == 0 and nsp == 0 and ("%d nodes" % nodes) in m.group(2) and float(m.group(11)) <= 1e-3 and wd / wo < wall_bound)
            record(f"{N}: automatic KLU, {tag} ({nodes} nodes): notice printed once per run, absent with AUTOKLU=0 and with `.options sparse`; orbit within 1e-3 of ripple of the SPARSE run; faster (x{wall_bound:g} bound)",
                   "pass" if ok else "FAIL", f"{m.group(2)}, wall {wd:.2f} vs {wo:.2f} s ({m.group(8)}), orbit dev {m.group(11)} (measured 2.8e-6 / 6.2e-5; -34 % / -27 %)")
    for case, ident_expected, bound in (("R3d (example 26)", True, 1.1), ("R4d (example 27)", False, 1.1), ("PEX 73 states, krylov=0", True, 0.9)):
        m = re.search(r"FORK %s: sequential median (\S+) s .*?PSS_FORK_COLS=2 median (\S+) s .*?; wrdata byte-identical: (\w+); iterations (\d+) / (\d+)" % re.escape(case), t)
        if not m:
            record(f"{N}: PSS_FORK_COLS=2, {case}", "FAIL", "no FORK line")
            continue
        r = float(m.group(2)) / float(m.group(1))
        # not byte-identical (R4d): the iterates differ, so the iteration count may too (second workstation: 5 against 8) -- only a converged result is required
        ok = (m.group(3) == str(ident_expected)) and r < bound and (abs(int(m.group(4)) - int(m.group(5))) <= 1 or (not ident_expected and int(m.group(5)) < 20))
        record(f"{N}: PSS_FORK_COLS=2, {case}: wrdata {'byte-identical' if ident_expected else 'NOT byte-identical (why it is not the default)'}; wall x{bound:g} bound",
               "pass" if ok else "FAIL", f"identical {m.group(3)}, wall ratio {r:.2f}, iterations {m.group(4)} / {m.group(5)} (measured 0.91 / 0.93 / 0.70 on a quiet development workstation, 0.75 / 0.79 / 0.63 on a loaded one)")
    m = re.search(r"MRHS R3d .pstb .*?n=(\d+): default median (\S+) s .*?PAC_MULTIRHS=1 (\S+) s .*?; output byte-identical: (\w+)", t)
    if not m:
        record(f"{N}: PAC_MULTIRHS=1", "FAIL", "no MRHS line")
    else:
        r = float(m.group(3)) / float(m.group(2))
        record(f"{N}: PAC_MULTIRHS=1 on a td .pstb: output byte-identical, wall within 20 % (no gain, hence opt-in)", "pass" if m.group(4) == "True" and 0.8 < r < 1.2 else "FAIL", f"ratio {r:.3f}, identical {m.group(4)} (measured 0.99)")
    for case, dev_bound, need_gain, it_tol in (("R3d (example 26)", 1e-9, False, 1), ("R4d (example 27)", 1e-3, False, 1), ("R7 BSIM4 (example 32)", 1e-5, True, 1)):
        m = re.search(r"SEED %s: default iter (\d+) residual (\S+) radius (\S+) (\S+), median (\S+) s .*?; PSS_SEED_STAGE1=1 iter (\d+) residual (\S+) radius (\S+) (\S+), median (\S+) s .*?; orbit deviation (\S+) of ripple" % re.escape(case), t)
        if not m:
            record(f"{N}: PSS_SEED_STAGE1=1, {case}", "FAIL", "no SEED line")
            continue
        # R4d: the Floquet multiplier is 0.994 (0.6 % from 1) and the seeded run's finite-difference estimate is not stable across machines (second workstation printed
        # 1.52 UNSTABLE(period-doubling) for the seeded arm, development workstation 0.994057 stable): the orbit and the residual are gated, the radius and verdict are only reported
        r4d = case.startswith("R4d")
        ok = (abs(int(m.group(1)) - int(m.group(6))) <= it_tol and float(m.group(7)) <= 1e-6 and (r4d or (abs(float(m.group(3)) - float(m.group(8))) <= 2e-3 and m.group(9).startswith("stable")))
              and float(m.group(11)) <= dev_bound and ((float(m.group(10)) / float(m.group(5)) < 0.95) if need_gain else True))
        record(f"{N}: PSS_SEED_STAGE1=1, {case}: converges to the same orbit (<= {dev_bound:g} of ripple{'' if case.startswith('R4d') else ', radius within 2e-3, stable'}){' and is faster than the default (-22 % measured, bound x0.95)' if need_gain else ''}",
               "pass" if ok else "FAIL", f"iterations {m.group(1)} / {m.group(6)}, residual {m.group(2)} / {m.group(7)}, radius {m.group(3)} {m.group(4)} / {m.group(8)} {m.group(9)}, wall {m.group(5)} / {m.group(10)} s, orbit dev {m.group(11)}")
    m = re.search(r"ARN PEX 490 states matrix-free: default (\d+) period runs, radius (\S+), median (\S+) s; PSS_KRYLOV_ARN_EARLY=0 (\d+) period runs, radius (\S+), median (\S+) s .*?output byte-identical: (\w+)", t)
    record(f"{N}: Arnoldi early stop: fewer period runs than PSS_KRYLOV_ARN_EARLY=0, same radius, byte-identical output", "pass" if m and int(m.group(1)) < int(m.group(4)) and m.group(2) == m.group(5) and m.group(7) == "True" else "FAIL",
           m.group(0)[:160] if m else "no ARN line")
    m = re.search(r"ARN R3d dense: default (\d+) period runs, PSS_FLOQUET_BASELINE=1 (\d+) period runs; output byte-identical: (\w+)", t)
    record(f"{N}: dense Floquet re-evaluation without its baseline period: one period run fewer, byte-identical output", "pass" if m and int(m.group(2)) - int(m.group(1)) == 1 and m.group(3) == "True" else "FAIL", m.group(0) if m else "no line")


def check_freerun_noise_560():
    N = "48 free-running .pxf / .pnoise"
    t = read_log("48")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    m = re.search(r"ADJ \(PAC_EXACT_JACOBIAN=1\) forward: .pxf adjoint=1 (ran|DID NOT RUN) the adjoint engine; .*?: (\S+)$", t, re.M)
    record(f"{N}: .pxf adjoint=1 equals the forward per-source engine (PAC_EXACT_JACOBIAN=1) to 1e-6 on R5, Vin and Vrefsrc at 20/55/90 kHz", "pass" if m and m.group(1) == "ran" and float(m.group(2)) <= 1e-6 else "FAIL", m.group(0)[-90:] if m else "no line")
    m = re.search(r"ADJ \(default forward\) forward: .*?: (\S+)$", t, re.M)
    record(f"{N}: the same against the default forward engine (the O(1/Points) floor of card #504) <= 1e-4", "pass" if m and float(m.group(1)) <= 1e-4 else "FAIL", m.group(1) if m else "no line")
    m = re.search(r"ADJ red proof PAC_FREERUN_PROJ=1: the forward answer moves by a relative (\S+); adjoint=1 with it (REFUSES)", t)
    record(f"{N}: red proof PAC_FREERUN_PROJ=1: the forward operator is different (> 1e-3) and the adjoint refuses instead of transposing it", "pass" if m and float(m.group(1)) > 1e-3 else "FAIL", m.group(0)[-110:] if m else "no line")
    m = re.search(r"PN low frequency, 10 Hz .. 1.8 kHz: .*? V/rtHz; spread (\S+) %", t)
    vals = [float(x) for x in re.findall(r"Hz (\d\.\d+e-\d+)", (re.search(r"PN low frequency[^\n]*", t) or [""])[0] if re.search(r"PN low frequency[^\n]*", t) else "")]
    record(f"{N}: .pnoise is finite and flat in band (10 Hz .. 1.8 kHz within 10 %; no phase-diffusion divergence)", "pass" if m and float(m.group(1)) <= 10 and len(vals) == 5 and all(1e-10 < v < 2e-10 for v in vals) else "FAIL", f"spread {m.group(1) if m else '?'} %, values {vals}")
    bands = re.findall(r"PN band (\S+) kHz \((\S+) fsw\): .pnoise (\S+) V/rtHz .*? ensemble (\S+) \(ts 5 ns, (\S+) dB\) and (\S+) \(ts 2.5 ns, (\S+) dB\)", t)
    if len(bands) != 3:
        record(f"{N}: .pnoise vs the cited transient-noise ensemble", "FAIL", f"{len(bands)} PN band lines")
    else:
        d25 = [float(b[6]) for b in bands]
        d5 = [float(b[4]) for b in bands]
        pooled = 10 * math.log10(sum(10 ** (d / 10) for d in d25) / 3)
        # the cited ensemble (ts 2.5 ns) is within 0.6 dB per band and 0.3 dB pooled of the engine; the 5 ns ensemble lies below the 2.5 ns one (finite noise bandwidth)
        ok = max(abs(d) for d in d25) <= 0.6 and abs(pooled) <= 0.3 and all(a < b for a, b in zip(d5, d25))
        record(f"{N}: .pnoise at 15 / 40 / 85 kHz against the CITED transient-noise ensemble (ts 2.5 ns): each band within 0.6 dB, pooled within 0.3 dB, the 5 ns ensemble below it",
               "pass" if ok else "FAIL", f"2.5 ns: {d25} dB, pooled {pooled:+.2f}; 5 ns: {d5} dB (cited: -0.39/+0.30/+0.37, pooled +0.09)")
        # discriminates: a +2 dB engine would fail the same gate
        record(f"{N}: the ensemble gate can fail (a +2 dB engine is outside it)", "pass" if max(abs(d - 2.0) for d in d25) > 0.6 else "FAIL", "")
    m = re.search(r"GUARD .pnoise up to 200 kHz \(0.47 fsw\): (.*)", t)
    record(f"{N}: a .pnoise sweep above 0.25 fsw is refused (band guard)", "pass" if m and "REFUSED" in m.group(1) and "0.25" in m.group(1) else "FAIL", (m.group(1)[:100] if m else "no line"))


def check_omp_560():
    N = "49 OpenMP threads"
    t = read_log("49")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return
    cores = {}
    for m in re.finditer(r"OMPRUN (.*?): wall (\S+) s, CPU (\S+) s, cores used \(CPU/wall\) (\S+), wrdata md5 (\w+)", t):
        cores[m.group(1)] = float(m.group(4))
    need = ["default", "OMP_NUM_THREADS=2", "set num_threads=2", "ompthreads.py 2", "set num_threads=1"]
    if any(k not in cores for k in need):
        record(N, "FAIL", f"missing arms: {sorted(cores)}")
        return
    d = cores["default"]
    record(f"{N}: the spinit default uses many threads (CPU/wall >= 3 on a 300-BSIM4 deck)", "pass" if d >= 3.0 else "FAIL", f"{d:.2f} (measured 7.6)")
    record(f"{N}: OMP_NUM_THREADS=2 does NOT limit it (CPU/wall within 25 % of the default)", "pass" if cores["OMP_NUM_THREADS=2"] >= 0.75 * d else "FAIL", f"{cores['OMP_NUM_THREADS=2']:.2f} vs {d:.2f}")
    record(f"{N}: `set num_threads=2` and tools/ompthreads.py (SPICE_SCRIPTS) limit it to about 2 (1.2 .. 2.6)", "pass" if all(1.2 <= cores[k] <= 2.6 for k in ("set num_threads=2", "ompthreads.py 2")) else "FAIL",
           f"{cores['set num_threads=2']:.2f} / {cores['ompthreads.py 2']:.2f}")
    record(f"{N}: `set num_threads=1` runs on one core (<= 1.3)", "pass" if cores["set num_threads=1"] <= 1.3 else "FAIL", f"{cores['set num_threads=1']:.2f}")
    record(f"{N}: the output of every arm is byte-identical", "pass" if "OMPSAME: the wrdata file of every arm byte-identical: True" in t else "FAIL", "")


def check_xspice_pac_560():
    N = "50 td .pac through s_xfer"
    t = read_log("50")
    if t is None:
        record(N, "COULD-NOT-RUN", "log missing")
        return

    def line(tag):
        m = re.search(r"R3XS %s: (\d+) hidden states; difference to the cited FRA: (.*?); wall" % re.escape(tag), t)
        if not m:
            return None
        v = [(float(a), float(b)) for a, b in re.findall(r"kHz (\S+) dB / (\S+) deg", m.group(2))]
        return int(m.group(1)), v
    a = line("pp2000")
    if not a:
        record(f"{N}: R3 + s_xfer", "FAIL", "no R3XS pp2000 line")
    else:
        w1, w2 = max(abs(x) for x, y in a[1]), max(abs(y) for x, y in a[1])
        record(f"{N}: R3 buck with an s_xfer compensator, td .pac (pacpoints 2000, 2 hidden states) equals the cited FRA at 4 / 20 / 50 kHz (0.03 dB / 0.15 deg)", "pass" if a[0] == 2 and len(a[1]) == 3 and w1 <= 0.03 and w2 <= 0.15 else "FAIL",
               f"worst {w1:.4f} dB / {w2:.3f} deg (measured 0.0004 / 0.002)")
    b = line("default grid")
    record(f"{N}: on the default recording grid it is within 0.3 dB / 1 deg of the FRA (the grid is the error)", "pass" if b and max(abs(x) for x, y in b[1]) <= 0.3 and max(abs(y) for x, y in b[1]) <= 1.0 else "FAIL", str(b))
    c = line("PAC_HID_RED=1")
    record(f"{N}: red proof PAC_HID_RED=1 (hidden states do not act on the nodes) is more than 1 dB off the FRA", "pass" if c and max(abs(x) for x, y in c[1]) > 1.0 else "FAIL", str(c))
    m = re.findall(r"XSLTI s_xfer td .pac at pacpoints=(\d+): worst relative difference to plain .ac (\S+) over", t)
    rich = re.search(r"XSLTI Richardson .*?: (\S+)$", t, re.M)
    d = dict((int(p), float(v)) for p, v in m)
    ok = 100000 in d and 50000 in d and d[100000] <= 2e-2 and d[50000] > 1.5 * d[100000] and rich and float(rich.group(1)) <= 5e-4
    record(f"{N}: a linear s_xfer on a periodic orbit: td vs .ac is first order in the step (the error at 50000 points is > 1.5x that at 100000, <= 2e-2) and Richardson extrapolation reaches 5e-4", "pass" if ok else "FAIL",
           f"{d}, Richardson {rich.group(1) if rich else None} (measured 1.9e-2 / 9.4e-3 / 8.6e-5)")
    for key, pat in (("slew model that slews on the orbit", "slew-rate limiting of a 'slew' code model is active"), ("hyst model that flips on the orbit", "Boolean state of a 'hyst' code model changes on the orbit"),
                     (".pnoise through s_xfer", "'s_xfer'")):
        m = re.search(r"XSREF %s: (.*)" % re.escape(key), t)
        record(f"{N}: refused: {key}", "pass" if m and "REFUSED" in m.group(1) and pat in m.group(1) else "FAIL", (m.group(1)[:100] if m else "no line"))
    m = re.search(r"XSREF .pxf adjoint=1 through s_xfer: falls back to the per-source engine \(message: (\w+)\), differs from .pac by (\S+)", t)
    record(f"{N}: .pxf adjoint=1 through s_xfer says it cannot transpose the hidden states and falls back; the result equals .pac (1e-9)", "pass" if m and m.group(1) == "True" and float(m.group(2)) <= 1e-9 else "FAIL", m.group(0)[-80:] if m else "no line")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default=DEFAULT_PREFIX)
    ap.add_argument("--skip-run", action="store_true",
                     help="reuse the existing run_all/ output instead of re-running "
                          "(for iterating on this script only)")
    args = ap.parse_args()

    if not args.skip_run:
        run_all_sh(args.prefix)
    elif not os.path.isdir(OUT):
        record("run_all.sh", "COULD-NOT-RUN", "--skip-run given but no prior output dir")

    # 01: stb LDO loop
    t = read_log("01")
    check_contains("01 stb LDO loop", "phase margin line", t,
                    "phase margin 35.82 deg at 43488.8 Hz (falling 0 dB crossing)")
    check_contains("01 stb LDO loop", "no gain-margin crossing", t,
                    "gain margin: no -180 deg crossing of T")

    # 02: stb DM/CM diffamp
    t = read_log("02")
    check_contains("02 stb DM/CM", "DM phase margin", t, "phase margin 94.45 deg")
    check_contains("02 stb DM/CM", "CM phase margin", t, "phase margin 56.88 deg")

    # 03: stb device mode
    t = read_log("03")
    check_contains("03 stb device mode", "rising 0 dB crossing, PM=258.3", t,
                    "phase margin 258.3 deg at 3288.26 Hz (rising 0 dB crossing)")

    # 04: stb sens= -- card #507 second-reader correction: the plot prefix
    # (stb2/stb4/...) is not predictable, so this deck uses `display` to
    # show the current plot's own name, and references the vectors it just
    # produced WITHOUT any plot prefix (they are already the current plot).
    t = read_log("04")
    check_contains("04 stb sens=", "display shows the current plot's own name", t,
                    "Name: stb4 (STB dMargins/dp)")
    check_float_near("04 stb sens=", "dPM/dR1o (no plot prefix)", t,
                      r"dpm_dr1o\s*=\s*(-?[\d.eE+-]+)", -6.28493e-08, rel=1e-3)
    check_float_near("04 stb sens=", "dPM/dCc (no plot prefix)", t,
                      r"dpm_dcc\s*=\s*(-?[\d.eE+-]+)", 3.24916e13, rel=1e-3)

    # 05: stb instant
    t = read_log("05")
    check_contains("05 stb instant", "mid-ramp phase margin line", t,
                    "phase margin 18.36 deg at 69243.1 Hz (falling 0 dB crossing)")
    check_float_near("05 stb instant", "stb2.phasemargin", t,
                      r"stb2\.phasemargin\s*=\s*([\d.eE+-]+)", 18.36451, rel=1e-3)

    # 06: pss method=shoot (s85: iter=3 -> iter=1, spectral radius 0.223059 -> 0.223111; same orbit to 5 digits)
    t = read_log("06")
    check_contains("06 pss shoot", "CONVERGED line", t,
                    "PSS-SHOOT: CONVERGED iter=1 residual=2.358290e-14 spectral_radius=0.223111 stable")

    # 07: pac solver=td (s136, card #625: recording reltol 1e-8 moved the 1 kHz value by 1e-12 relative; PAC_RECORD_RELTOL=1e-10 gives
    # 9.061344830139665e-01 -5.175828063966958e-02, the s133 value was 9.061344830130268e-01 -5.175828064053489e-02)
    t = read_file("07", "pac_td_out.txt")
    check_contains("07 pac solver=td", "1kHz real part", t, "9.061344830139495e-01")
    check_contains("07 pac solver=td", "1kHz imag part", t, "-5.175828064063284e-02")

    # 08: pstb
    t = read_log("08")
    check_contains("08 pstb", "phase margin line", t,
                    "pstb: phase margin 9.778 deg at 1.10724e+07 Hz")   # s133: 1.10723e+07 (card #625)

    # 09: pxf
    t = read_file("09", "pxf_out.txt")
    check_contains("09 pxf", "pxf_vin 1kHz real part matches pac", t, "9.061344830139")
    check_contains("09 pxf", "pxf_iinj 1kHz real part (1000x pxf_vin)", t, "9.061344830139")

    # 10: saltation / ideal switch
    t_out = read_file("10", "pac_r6o_1k.txt")
    # card #535: the s85 drop (6.204 -> 6.157 V/V at 1 kHz, 248.7 kHz point +1.85 %) was a spurious
    # t=0 saltation jump, fixed in aug-535 (a496a1bac). The pinned reference is the exact LPTV
    # model of this deck (piecewise-linear switched circuit, matrix exponentials, event-time
    # perturbation; cross-checked by an event-exact nonlinear time-domain run to 5e-9):
    # 1 kHz 6.20355 V/V, -1.8874 deg; 248744 Hz 4.5456e-3 V/V. Series s63 and aug-535 sit within
    # +0.03 % / -0.04 % of it (pac_saltation(b) reltol sweep spread 0.014 %); the bad s85 sat at
    # -0.75 % / +1.85 %, so 0.2 % / 0.2 deg / 0.3 % separates them with margin.
    check_complex_near("10 saltation ideal sw", "1kHz |H| and phase vs the exact LPTV model",
                       t_out, 1000.0, 6.20355, -1.8874, rel=2e-3, deg_tol=0.2)
    check_complex_near("10 saltation ideal sw", "248.744kHz |H| vs the exact LPTV model",
                       t_out, 248744.0, 4.5456e-3, None, rel=3e-3, deg_tol=None)
    t = read_log("10")
    check_contains("10 saltation ideal sw", "fs/2 unvalidated-margin warning fires", t,
                    "UNVALIDATED margin near the switching Nyquist frequency")

    # 11: pac_perturbed
    t = read_log("11")
    check_contains("11 pac_perturbed", "10kHz result", t, "-2.0542dB/-29.7349deg")
    check_contains("11 pac_perturbed", "30kHz result", t, "-6.7789dB/-59.7329deg")

    # 12: fra --adaptive-settle
    t = read_log("12")
    check_contains("12 fra adaptive-settle", "settle line", t, "td=6e-05s")
    check_contains("12 fra adaptive-settle", "1/1 points ok", t, "1/1 frequency points ok")

    # 13: convmeas
    t = read_log("13")
    check_contains("13 convmeas", "ripple pk-pk", t, "0.5068 mV")
    check_contains("13 convmeas", "overshoot/undershoot", t,
                    "overshoot 37.242%, undershoot 13.857%")

    # 14: loss.py
    t = read_log("14")
    check_contains("14 loss.py", "efficiency line", t, "efficiency = 84.6933 %")
    check_contains("14 loss.py", "diode_conduction", t, "diode_conduction      0.574108")
    check_contains("14 loss.py", "residual", t, "-0.1373% of P_in")

    # 15: mc.py
    t = read_log("15")
    check_contains("15 mc.py", "correlation check pass", t, '"state": "pass"')
    check_contains("15 mc.py", "measured correlation ~0.81", t, "r=0.8082")

    # 16: yieldest.py
    t = read_log("16")
    check_float_near("16 yieldest.py", "beta_mpp", t,
                      r"beta_mpp=([\d.]+)", 4.0, abs_tol=0.01)
    check_float_near("16 yieldest.py", "importance-sampling/true ratio", t,
                      r"ratio=([\d.]+)", 0.888, rel=0.05)

    # 17: compopt.py
    t = read_log("17")
    check_contains("17 compopt.py", "target reached", t, "TARGET REACHED")
    check_contains("17 compopt.py", "independent stb.py check", t,
                    "PM = 65.000 deg @ 7.60034e+06 Hz")

    # 18: checkpoint save/restore
    t = read_log("18")
    if t is not None:
        n = t.count("vx_at_600u = 7.531775e-01")
        record("18 checkpoint", "pass" if n >= 2 else "FAIL",
               f"vx_at_600u = 7.531775e-01 appears {n} time(s) "
               "(expect >=2: uninterrupted run and restored run agree)")
    else:
        record("18 checkpoint", "COULD-NOT-RUN", "log file missing")

    # 19: once() vs without -- patch 0016, IN patches/series, confirmed present
    # in this build (card #507 correction: the page originally, wrongly,
    # claimed once() was not in the aaba50ba3f49-s54 build).
    t = read_log("19")
    check_float_near("19 once() vs without", "without once(): r1", t,
                      r"resistance\]\s*=\s*([\d.eE+-]+)", 996.2011, rel=1e-3)
    if t is not None:
        # vals: [without-r1, without-r2, with-r1, with-r2]
        vals = re.findall(r"resistance\]\s*=\s*([\d.eE+-]+)", t)
        ok = (len(vals) == 4 and vals[0] != vals[1] and vals[2] == vals[3])
        record("19 once() vs without: without-once r1!=r2, with-once r1==r2",
               "pass" if ok else "FAIL", f"parsed values={vals}")
    else:
        record("19 once() vs without: without-once r1!=r2, with-once r1==r2",
               "COULD-NOT-RUN", "log file missing")

    # 20: .pss method=shoot driven-detection red-proof (card #507 second
    # reader): a B-source clock using bare `time` is NOT detected as a
    # driven source -- expected to be REFUSED, not silently misread.
    t = read_log("20")
    check_contains("20 pss driven-detection red-proof", "documented refusal message", t,
                    "no V/I source with a time function was found")

    # 21: saltation on a switch node without a capacitor -- was exactly 0
    # (found by the #507 reader), fixed by card #510 (0060-0063): nonzero,
    # 13.28 V/V at 1 kHz, and the inductor-routed warning is printed.
    t = read_log("21")
    check_contains("21 saltation no-cap switch node (card #510 fixed)", "1kHz result nonzero", t,
                    "1.000000e+03\t1.321289e+01")   # 13.28437 up to s117; 13.32035 s118; 13.32032 s119; 13.21289 from s123 (routed jump after T, generic recursion, #583): 22.43 dB; a hand estimate (Vin + Vd) / Vramp x 1 / (1 - (f/f0)^2) gives 22.43 dB; no pac_perturbed reference at 1 kHz (not commensurate within its limits)
    check_contains("21 saltation no-cap switch node (card #510 fixed)", "warning printed", t,
                    "saltation via inductor-routed jump at node sw")


    # ---------------------------------------------------------------- card #534 (s85) additions
    # 22: closed-loop ideal-switch buck: .pss converges (card #498) -- and fails without exact event location
    t = read_log("22")
    check_contains("22 pss closed-loop sw", "CONVERGED line", t,
                    # spectral radius to 5 decimals: 0.980206 before s141, 0.980205 since #682 (Floquet columns use the Newton step rule)
                    "PSS-SHOOT: CONVERGED iter=4 residual=1.827417e-07 spectral_radius=0.98020")
    t = read_log("22off")
    check_contains("22 pss closed-loop sw", "red-proof: PSS_NO_EVENT_LOCATE=1 does NOT converge", t,
                    "FAILED to converge within 150 Newton iterations (final residual 2.167824e-05")

    # 23: .pstb on that orbit vs the independent perturbed-nonlinear Tian reference (pac_perturbed --loop-gain)
    t = read_file("23", "pstb_r6s.txt")
    ref = read_log("23ref")
    if t is None or ref is None:
        record("23 pstb closed-loop sw", "COULD-NOT-RUN", "pstb output or pac_perturbed log missing")
    else:
        rows = [ln.split() for ln in t.splitlines() if ln.strip() and not ln.strip().startswith("freq")]
        mine = {round(float(r[0])): complex(float(r[1]), float(r[2])) for r in rows}
        refs = {}
        for m in re.finditer(r"f_actual=(\d+)Hz.*?: T=([-\d.]+)dB/([-\d.]+)deg", ref):
            refs[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
        worst_db = worst_deg = 0.0
        det = []
        for f in (30000, 100000):
            if f not in mine or f not in refs:
                record("23 pstb closed-loop sw", "COULD-NOT-RUN", f"missing point {f}: {sorted(mine)} {sorted(refs)}")
                break
            db = 20 * math.log10(abs(mine[f]))
            deg = math.degrees(math.atan2(mine[f].imag, mine[f].real))
            ddeg = abs(deg - refs[f][1])
            ddeg = min(ddeg, 360 - ddeg)
            worst_db, worst_deg = max(worst_db, abs(db - refs[f][0])), max(worst_deg, ddeg)
            det.append(f"{f/1e3:g}k: pstb {db:.3f} dB/{deg:.2f} deg vs ref {refs[f][0]:.3f}/{refs[f][1]:.2f}")
        else:
            record("23 pstb closed-loop sw", "pass" if (worst_db <= 0.5 and worst_deg <= 3.0) else "FAIL",
                   f"worst {worst_db:.3f} dB / {worst_deg:.2f} deg (target 0.5/3); " + "; ".join(det))
    t = read_log("23")
    check_contains("23 pstb closed-loop sw", "PSS converged under .pstb", t, "PSS-SHOOT: CONVERGED iter=3 residual=5.936077e-07")

    # 24: .pnoise on the switched RC: values, and the per-device table must sum (in power) to the total
    t = read_file("24", "pnoise_r2.txt")
    if t is None:
        record("24 pnoise switched RC", "COULD-NOT-RUN", "pnoise_r2.txt missing")
    else:
        rows = [[float(x) for x in ln.split()] for ln in t.splitlines()[1:] if ln.strip()]
        r1k, r1m = rows[0], rows[-1]
        # card #625 (s136): recording reltol 1e-8; PAC_RECORD_RELTOL=1e-10 gives 4.443138e-09, s136 4.443128e-09 (2e-6 off), s133 4.444167e-09 (2.3e-4 off)
        record("24 pnoise switched RC: 1 kHz onoise_spectrum", "pass" if abs(r1k[1] / 4.443138e-9 - 1) < 1e-4 else "FAIL",
               f"{r1k[1]:.6e} V/rtHz (expect 4.443138e-09, the PAC_RECORD_RELTOL=1e-10 value; bound 1e-4)")
        record("24 pnoise switched RC: 1 MHz onoise_spectrum", "pass" if abs(r1m[1] / 7.86190e-10 - 1) < 1e-4 else "FAIL",
               f"{r1m[1]:.6e} V/rtHz (expect 7.86190e-10)")
        worst = max(abs(r[2] ** 2 + r[3] ** 2 - r[1] ** 2) / r[1] ** 2 for r in rows)
        record("24 pnoise switched RC: per-device sum of squares == total", "pass" if worst < 1e-9 else "FAIL",
               f"{len(rows)} points, worst relative {worst:.2e}")
        record("24 pnoise switched RC: all finite and positive", "pass" if all(r[1] > 0 and r[1] == r[1] for r in rows) else "FAIL", "")

    # 25: .pnoise track-and-hold vs the EXACT folded kT/C (OU restart) -- th_analytic_S0 re-implemented from tests/test_pnoise.py
    def th_s0(rtrack, c, T, tt):
        kT = 1.380649e-23 * 300.15
        tau = rtrack * c
        a = tt / tau
        Th = T - tt
        q = 1 - math.exp(-a)
        rho = math.exp(-a)
        varA = 2 * tau * (tt - tau * q)
        covAs = tau * q
        var0 = varA + Th * Th + 2 * Th * covAs
        c1 = (tau * q) * covAs + Th * rho * covAs + Th * tau * q + Th * Th * rho
        c2 = ((tau * q) * covAs * rho + Th * rho * rho * covAs + Th * tau * q * rho + Th * Th * rho * rho)
        return 2 * (kT / c) * (var0 + 2 * c1 + 2 * c2) / T
    t = read_file("25", "pnoise_th.txt")
    if t is None:
        record("25 pnoise track-and-hold", "COULD-NOT-RUN", "pnoise_th.txt missing")
    else:
        T = 1e-5
        S0 = th_s0(25e3, 10e-12, T, 0.2 * T)
        naive_db = 20 * math.log10(math.sqrt(2 * 1.380649e-23 * 300.15 / 10e-12 * T) / math.sqrt(S0))
        rows = [[float(x) for x in ln.split()] for ln in t.splitlines()[1:] if ln.strip()]
        worst, det = 0.0, []
        for f, v in rows:
            x = math.pi * f * T
            an = math.sqrt(S0 * (math.sin(x) / x) ** 2)
            d = 20 * math.log10(v / an)
            worst = max(worst, abs(d))
            det.append(f"{f:g} Hz {v:.4e} vs {an:.4e} ({d:+.3f} dB)")
        record("25 pnoise track-and-hold vs exact folded kT/C", "pass" if (len(rows) == 3 and worst <= 0.1) else "FAIL",
               "; ".join(det) + f"; worst {worst:.3f} dB (limit 0.1); naive kT/C ZOH is {naive_db:+.2f} dB away")
        record("25 pnoise track-and-hold: the check discriminates (naive ZOH is > 1 dB off)",
               "pass" if naive_db > 1.0 else "FAIL", f"{naive_db:.2f} dB")

    # 26/27/28: dead-time converters
    t = read_log("26")
    check_converged("26 pss R3d (R3 default 300us/200/150)", "CONVERGED line", t, "PSS-SHOOT", 8, 1e-6, 0.980658)
    t = read_log("27")
    check_converged("27 pss R4d", "CONVERGED line", t, "PSS-SHOOT", 5, 1e-6, 0.994028)
    t = read_log("28")
    # development workstation 38.6 deg at 33473.1 Hz, second workstation 38.59 deg at 33491.7 Hz (card #550): 0.01 deg / 5.6e-4 rel apart; bounds are 3x the spread
    check_float_near("28 pstb R3d", "phase margin 38.6 deg (dev. ws.; second ws. 38.59)", t, r"pstb: phase margin ([\d.]+) deg", 38.6, abs_tol=0.03)
    check_float_near("28 pstb R3d", "crossover 33473.1 Hz (dev. ws.; second ws. 33491.7)", t, r"pstb: phase margin [\d.]+ deg at ([\d.eE+]+) Hz", 33473.1, rel=2e-3)
    check_float_near("28 pstb R3d", "PM within 1 deg of R3's FRA reference 38.46 deg (plausibility only)", t,
                      r"pstb: phase margin ([\d.]+) deg", 38.46, abs_tol=1.0)

    # 30: free-running .pac vs pac_perturbed --freerun reference values (tests/test_pac.py, card #503 third pass table)
    t = read_file("30", "pac_r5.txt")
    lg = read_log("30")
    if t is None or lg is None:
        record("30 pac free-running R5", "COULD-NOT-RUN", "output missing")
    else:
        refs30 = {20000: (-0.009, -0.09), 60000: (0.279, -0.45), 100000: (0.769, -2.13)}
        rows = [ln.split() for ln in t.splitlines()[1:] if ln.strip()]
        worst_db = worst_deg = 0.0
        det = []
        for r in rows:
            f = round(float(r[0]))
            h = complex(float(r[1]), float(r[2]))
            db, deg = 20 * math.log10(abs(h)), math.degrees(math.atan2(h.imag, h.real))
            worst_db = max(worst_db, abs(db - refs30[f][0]))
            worst_deg = max(worst_deg, abs(deg - refs30[f][1]))
            det.append(f"{f/1e3:g}k {db:+.3f} dB/{deg:+.2f} deg (ref {refs30[f][0]:+.3f}/{refs30[f][1]:+.2f})")
        ok = len(rows) == 3 and worst_db <= 0.5 and worst_deg <= 3.0 and "PSS-FREERUN: running .pac/.pstb" in lg
        record("30 pac free-running R5 vs pac_perturbed --freerun", "pass" if ok else "FAIL",
               f"worst {worst_db:.3f} dB / {worst_deg:.2f} deg (target 0.5/3); " + "; ".join(det))
    check_converged("30 pac free-running R5", "free-running orbit converged", lg, "PSS-FREERUN", 17, 2e-5, 0.036816, rho_tol=1e-3,
                    verdict="stable[arnoldi]", period_ref=2.36584808e-06)
    check_contains("30 pac free-running R5", "10 % band warning still printed at 60 kHz", lg,
                    "free-running (trigger=) .pac/.pstb result at f=60000 Hz is inside the UNVALIDATED margin")
    t1 = read_file("30t1", "pac_r5.txt")
    if t is not None and t1 is not None:
        record("30 pac free-running R5: 1-thread result equals 2-thread result", "pass" if t == t1 else "FAIL",
               "identical wrdata" if t == t1 else "wrdata differs")
    else:
        record("30 pac free-running R5: 1-thread result equals 2-thread result", "COULD-NOT-RUN", "30t1 output missing")

    # 29: speed switches (rounds 13-15): identity must hold; timing only checked where the margin is large
    t = read_log("29")
    if t is None:
        record("29 speed switches", "COULD-NOT-RUN", "log missing")
    else:
        med, mn = {}, {}
        for m in re.finditer(r"SPEED (\w+)\s+(\S+)\s+median\s+([\d.]+) s\s+\[min ([\d.]+)", t):
            med[(m.group(1), m.group(2))] = float(m.group(3))
            mn[(m.group(1), m.group(2))] = float(m.group(4))
        for k in (("pex_tran", "default"), ("pex_tran", "RESBASE=0"), ("pss_r3d", "default"), ("pss_r3d", "FLATB=0")):
            if k not in med:
                record("29 speed switches", "COULD-NOT-RUN", f"missing timing {k}")
        if ("pss_r3d", "FLATB=0") in med:
            check_contains("29 speed: FLATB=0 output", "bit-identical to default", t, "IDENT pss_r3d FLATB=0: output bit-identical to default")
            r = med[("pss_r3d", "default")] / med[("pss_r3d", "FLATB=0")]
            record("29 speed: flat B-source evaluator is not slower (PSS R3d)", "pass" if r < 1.0 else "FAIL",
                   f"default/FLATB=0 = {r:.3f}")
        if ("pex_tran", "RESBASE=0") in med:
            # the machine is shared and loaded: compare best-of-n (minimum), the median is too noisy for a hard check
            r = mn[("pex_tran", "default")] / mn[("pex_tran", "RESBASE=0")]
            # card #550: record() only knows pass/FAIL/COULD-NOT-RUN; the "info" state asserted (AssertionError, whole suite aborted) on a
            # loaded development-workstation run (r >= 1), so the non-gating case is printed, not recorded.
            (record if r < 1.0 else (lambda n, st, d: print(f"[info           ] {n}: {d}")))("29 speed: R/C base caches are faster (KLU PEX tran, e_pex_tree; best-of-n)", "pass",
                   f"min default/min RESBASE=0 = {r:.3f}" + ("" if r < 1.0 else
                   " -- timing on a shared machine; not a gate (2026-10-02: 1.475 on the development workstation at load >20 while the same"
                   " builds read 0.66 on the idle second workstation, research section 20). The identity checks below are the gate."))
            for sw in ("CAPBASE=0", "RESBASE=0", "KLU_REFSCALE=1"):
                m = re.search(rf"IDENT pex_tran {re.escape(sw)}: max \|delta\| / column range = ([\d.eE+-]+)", t)
                if not m:
                    record(f"29 speed: {sw} output identity", "FAIL", "no IDENT line (or row counts differ)")
                else:
                    v = float(m.group(1))
                    record(f"29 speed: {sw} output within 1e-6 of default", "pass" if v <= 1e-6 else "FAIL", f"{v:.1e}")


    # 31: campaign/mc thread rule (pure function of tools/ompthreads.py)
    t = read_log("31")
    check_contains("31 thread rule", "32 cores", t,
                    "THREADRULE cores=32: N=1:8  N=2:8  N=4:8  N=6:5  N=8:4  N=12:2  N=16:2  N=32:1")
    check_contains("31 thread rule", "8 cores", t,
                    "THREADRULE cores=8: N=1:8  N=2:4  N=4:2  N=6:1  N=8:1  N=12:1  N=16:1  N=32:1")
    check_contains("31 thread rule", "settings", t,
                    "THREADRULE settings: 'auto'->2  3->3  'default'->None  0->None")


    # ---------------------------------------------------------------- card #543 (s96) additions
    check_r7()
    check_pac_perturbed_s96()
    check_freerun_verdict()
    check_pex_krylov()
    check_pwl()
    check_asan(args.prefix)

    # ---------------------------------------------------------------- card #552 (s108) additions
    check_pnoise_opoint()
    check_th_mos()
    check_r7_552()
    check_r8()
    check_salt_charge_552()
    check_xspice_552()
    check_r9_552()

    # ---------------------------------------------------------------- card #560 (s115) additions
    check_fold_560()
    check_memory_560()
    check_speed_560()
    check_freerun_noise_560()
    check_omp_560()
    check_xspice_pac_560()

    n_pass = sum(1 for _, s, _ in results if s == "pass")
    n_fail = sum(1 for _, s, _ in results if s == "FAIL")
    n_cnr = sum(1 for _, s, _ in results if s == "COULD-NOT-RUN")
    print(f"\n{n_pass} pass, {n_fail} FAIL, {n_cnr} COULD-NOT-RUN")
    if not RAN_OK or n_cnr:
        return 2
    if n_fail:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
