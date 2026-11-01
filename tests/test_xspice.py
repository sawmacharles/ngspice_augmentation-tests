#!/usr/bin/env python3
"""XSPICE code models under `.pss method=shoot` and td `.pac/.pstb/.pxf/.pnoise` (card #551, 2026-10-02).

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_xspice.py [--quick] [--no-r9] [--r9-only]
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN (a FAIL wins)

Outputs under ~/Work/data/ngspice_augmentation/551/test (never the repo).  Fixtures: tests/fixtures/xspice_*.cir, r9_*.cir.

What is checked (every check can fail; red proofs at the end make the bounds fail on purpose)
  1. td `.pac` of STATELESS analog code models (gain summer mult divide limit pwl aswitch pswitch ...) vs plain `.ac` of the same deck
     (a slow clock on an unrelated net gives `.pss` a period).  Before #551 the response carried the model's DC rhs (gain=2.5, out_offset=0.1
     gave 2.6): pac_capture_b0() called the XSPICE AC load with the analysis type still DC.  Bound 1e-6.
  2. td `.pac` REFUSES the stateful analog models without a td model (oneshot, clock sources) with a message naming the model; the digital
     divider deck runs since card #714 (edge shifts, tests/test_digital_pac.py checks values) and PAC_DIG_REFUSE=1 restores its refusal.
  3. `.pss method=shoot` of stateless analog models vs the settled `.tran` (0.02 % of ripple), of stateful analog models (state carried by
     XSsnap; their own first-order cm_analog_integrate error bounds the agreement: s_xfer/int/slew/inductoric <= 1.5 %), and of digital /
     hybrid models (adc_bridge dac_bridge d_tff d_dff d_jkff d_srff d_dlatch d_srlatch gates d_fdiv d_lut d_state d_osc d_pwm; Points 4000:
     <= 0.4 %).  Before #551: stateful models converged to a WRONG orbit (s_xfer 4.9 %, int 153 %, d_dt 250 % of ripple), digital models crashed.
  4. unvalidated models are REFUSED by name; a digital loop whose state is not periodic with the analysis period (d_tff toggling at the
     clock period) is REFUSED, at twice the period it converges.
  5. red proofs: PSS_XS_NO_STATE=1 (restart without the carried XSPICE state) must exceed the bounds on a held-state deck (comparator +
     d_dff sample-and-hold) and on s_xfer.
  7. (card #557) td `.pac/.pstb/.pxf` THROUGH the stateful analog models s_xfer, int, d_dt (and slew / hyst on an orbit that stays on one branch):
     the models' cm_analog_alloc state slots are Floquet states, the per-step map is found by finite differences of the model (pacshoot.c
     pac_hid_*).  Checked on a periodic LTI orbit (sine source, nonconstant state) against (a) the exact backward-Euler discrete transfer of the
     model chain (the transient's own cm_analog_integrate is first order: that is what the linearisation must reproduce, 2e-5), (b) plain `.ac`
     by Richardson extrapolation of two grids (the O(h) difference of the models' own integral vanishes), (c) the continuous hand calculation of
     `.ac` itself; closed loops (s_xfer compensator in the R3 buck, tests/fixtures/xspice_r3_sxfer.cir) vs the independent perturbed-periodic-
     steady-state reference (tools/pac_perturbed.py); red proofs PAC_HID_RED=1 (hidden states do not act) and PAC_HID_NO_JCOL=1 (declared PARTIALs
     instead of finite differences: s_xfer's PARTIAL is not the derivative) must fail the bound; remaining refusals (pnoise, method=gear, slew with
     slew-rate limiting on the orbit, hyst that flips, oneshot) still refuse.
  6. R9 (refs/r9_digital_buck.cir, digitally controlled buck): orbit vs settled tran, red proof, limit-cycling variant (7-bit DPWM) refused at T and
     converging at 2 T.  Slow (minutes): skipped with --no-r9 / --quick.
"""
import bisect
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
NG = os.environ.get("NGSPICE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-557-4/bin/ngspice-aug")
OUT = os.path.expanduser("~/Work/data/ngspice_augmentation/557/test")
QUICK = "--quick" in sys.argv
NOR9 = "--no-r9" in sys.argv or QUICK
R9ONLY = "--r9-only" in sys.argv
results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"[{state:13s}] {name}: {detail}", flush=True)


def run(deck, name, env=None, timeout=1800):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name + ".cir")
    with open(p, "w") as fh:
        fh.write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    try:
        r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", NG, "-b", p],
                           capture_output=True, text=True, cwd=OUT, env=e, timeout=timeout)
    except (FileNotFoundError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{name}: {ex}")
    return (r.stdout + r.stderr).replace("\r", "\n"), r.returncode


def fixture(name):
    with open(os.path.join(FIX, name)) as fh:
        txt = fh.read()
    txt = txt.replace("STATEFILE", os.path.join(FIX, "xspice_state_cnt.in"))
    txt = re.sub(r"^\.include\s+\S*lgprobe\.lib", ".include " + os.path.join(os.path.dirname(HERE), "stb", "lgprobe.lib"), txt, flags=re.M)
    return re.sub(r"\n\.end\s*$", "\n", txt.rstrip() + "\n")


def rd(path):
    if not os.path.exists(path):
        raise CouldNotRun(f"missing {path}")
    rows = []
    for l in open(path):
        try:
            rows.append([float(x) for x in l.split()])
        except ValueError:
            pass
    if not rows:
        raise CouldNotRun(f"empty {path}")
    return rows


def interp(ts, vs, t):
    i = bisect.bisect_left(ts, t)
    if i <= 0:
        return vs[0]
    if i >= len(ts):
        return vs[-1]
    a = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
    return vs[i - 1] + a * (vs[i] - vs[i - 1])


WR = "set filetype=ascii\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"


def pss_vs_tran(fx, tag, f, node="out", vecs="v(out)", pts=400, tstab=None, tstop=None, extra="", step=None, env=None, tmin=1e-9, cols=(1,)):
    """max |orbit - settled tran| over one period as % of the tran ripple, per column; returns (list, log, converged)."""
    T = 1.0 / f
    tstab = tstab or 40 * T
    tstop = tstop or 60 * T
    step = step or T / 400
    base = fixture(fx) + ".options reltol=1e-5\n"
    for suf in ("pss", "tran"):
        pth = os.path.join(OUT, f"{tag}_{suf}.txt")
        if os.path.exists(pth):
            os.remove(pth)
    log, _ = run(base + ".control\n" + WR + f"pss {f:.15g} {tstab:g} {node} {pts} 10 150 1e-6 uic method=shoot {extra}\n"
                 f"wrdata {OUT}/{tag}_pss.txt {vecs}\nquit\n.endc\n.end\n", f"{tag}_pss", env=env)
    conv = "PSS-SHOOT: CONVERGED" in log
    if not conv:
        return None, log, False
    run(base + ".control\n" + WR + f"tran {step:g} {tstop:g} uic\nwrdata {OUT}/{tag}_tran.txt {vecs}\nquit\n.endc\n.end\n", f"{tag}_tran")
    P, R = rd(os.path.join(OUT, f"{tag}_pss.txt")), rd(os.path.join(OUT, f"{tag}_tran.txt"))
    rt = [r[0] for r in R]
    t0 = rt[-1] - T
    pt = [0.0] + [r[0] for r in P]
    Tp = pt[-1]
    out = []
    N = 800
    for c in cols:
        pv = [P[0][c]] + [r[c] for r in P]
        rv = [r[c] for r in R]
        rp = [interp(rt, rv, t0 + T * j / N) for j in range(N + 1)]
        pp = [interp(pt, pv, Tp * j / N) for j in range(N + 1)]
        rip = max(rp) - min(rp)
        out.append(100 * max(abs(a - b) for j, (a, b) in enumerate(zip(pp, rp)) if T * j / N >= tmin) / rip)
    return out, log, True


def tag_of(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s)


# ---------------------------------------------------------------------------------------------------------------------
STATELESS = ["gain", "summer", "mult", "divide", "limit", "climit", "ilimit", "pwl", "aswitch", "pswitch"]
# bound (% of ripple, column list) per stateful class-2 model; the first-order cm_analog_integrate error dominates (docs/pss-pac.md #551)
STATEFUL = {"s_xfer": (1.5, (1, 2)), "int": (1.0, (1, 2)), "d_dt": (0.5, (1,)), "slew": (1.0, (1, 2)), "hyst": (0.05, (1, 2)),
            "oneshot": (0.5, (1,)), "capacitor": (0.5, (1, 2)), "inductor": (1.5, (1,)), "square": (0.5, (1, 2)),
            "sine": (0.05, (1, 2)), "triangle": (0.1, (1, 2))}
DIGITAL = {"tff": 4e-6, "logic": 4e-6, "seq": 4e-6, "fdiv_lut": 4e-6, "osc_pwm": 4e-6, "pull": 4e-6, "dsample": 4e-6}


def check_stateless_pss():
    worst, parts = 0.0, []
    for m in STATELESS:
        try:
            e, log, ok = pss_vs_tran(f"xspice_an_{m}.cir", "sl_" + m, 100e3, vecs="v(out) v(o1)", step=2e-9, cols=(1, 2))
        except CouldNotRun as ex:
            record(f"PSS stateless {m}", "COULD-NOT-RUN", str(ex))
            continue
        if not ok:
            record(f"PSS stateless {m}", "FAIL", "did not converge")
            continue
        parts.append(f"{m} {max(e):.4f}")
        worst = max(worst, max(e))
    record("PSS orbit vs settled tran: stateless analog code models (10)", "pass" if worst <= 0.02 and len(parts) == len(STATELESS) else "FAIL",
           f"worst {worst:.4f} % of ripple (bound 0.02); " + ", ".join(parts))


def check_stateful_pss():
    for m, (bound, cols) in STATEFUL.items():
        extra = "driven=1" if m in ("square", "sine", "triangle") else ""
        try:
            e, log, ok = pss_vs_tran(f"xspice_an_{m}.cir", "sf_" + m, 100e3, vecs="v(out) v(o1)", step=2e-9, extra=extra, cols=cols)
        except CouldNotRun as ex:
            record(f"PSS stateful {m}", "COULD-NOT-RUN", str(ex))
            continue
        if not ok:
            record(f"PSS stateful {m}", "FAIL", "did not converge: " + " | ".join(l for l in log.split("\n") if "FAILED" in l or "REFUSED" in l)[:200])
            continue
        w = max(e)
        record(f"PSS orbit vs settled tran: stateful analog {m}", "pass" if w <= bound else "FAIL",
               f"max {w:.4f} % of ripple per column {['%.4f' % x for x in e]} (bound {bound}; XSPICE integrators are first order, "
               f"the error shrinks with Points)")


def check_digital_pss():
    for n, T in DIGITAL.items():
        try:
            e, log, ok = pss_vs_tran(f"xspice_dig_{n}.cir", "dg_" + n, 1.0 / T, vecs="v(out)", pts=4000, tstab=40e-6, tstop=240e-6,
                                     extra="driven=1", step=0.25e-9)
        except CouldNotRun as ex:
            record(f"PSS digital {n}", "COULD-NOT-RUN", str(ex))
            continue
        if not ok:
            record(f"PSS digital {n}", "FAIL", "did not converge: " + " | ".join(l for l in log.split("\n") if "FAILED" in l or "REFUSED" in l)[:300])
            continue
        m = re.search(r"digital node mismatches (\d+), analog node drift (\S+)", log)
        record(f"PSS orbit vs settled tran: digital/hybrid deck {n}", "pass" if e[0] <= 0.4 and m and m.group(1) == "0" else "FAIL",
               f"{e[0]:.4f} % of ripple (bound 0.4, Points 4000, tran step 0.25 ns); hidden-state check: {m.group(0) if m else 'n/a'}")


# --- td PAC --------------------------------------------------------------------------------------------------------
CLK = "Vclk clk 0 PULSE(0 1 0 1n 1n 500u 1m)\nRclk clk clkgnd 100k\nCclk clkgnd 0 10n\n"
OPTS = ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n"


def pac_fixture(m):
    txt = fixture(f"xspice_an_{m}.cir")
    lines = [l for l in txt.split("\n") if not l.lower().startswith(("rl ", "co ", "vin "))]
    # LTI-biased deck: DC bias + AC 1 source, resistive load on o1
    return "\n".join(lines) + "\nVin in 0 DC 0.3 AC 1\nRl o1 0 1k\n"


def pac_run(m, kind):
    body = "* xspice pac\n" + CLK + pac_fixture(m) + OPTS
    wr = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata %s/%s_%s.txt o1\nquit\n.endc\n.end\n"
    if kind == "ac":
        return run(body + ".ac dec 5 1k 100k\n" + wr % (OUT, "pac_ac", m), f"pac_ac_{m}")
    return run(body + '.pac 1k 0 clkgnd 2000 1 8 1e-6 dec 5 1k 100k method="shoot" solver="td" pacpoints=100000\n' + wr % (OUT, "pac_td", m),
               f"pac_td_{m}")


def rdc(path):
    if not os.path.exists(path):
        raise CouldNotRun(f"missing {path}")
    rows = []
    for l in open(path).read().split("\n")[1:]:
        t = l.split()
        if len(t) >= 3:
            rows.append((float(t[0]), complex(float(t[1]), float(t[2]))))
    if not rows:
        raise CouldNotRun(f"empty {path}")
    return rows


def check_pac_stateless():
    worst, parts = 0.0, []
    for m in ["gain", "summer", "mult", "divide", "limit", "pwl", "aswitch", "pswitch", "climit", "ilimit"]:
        try:
            for f in (f"pac_td_{m}.txt", f"pac_ac_{m}.txt"):
                if os.path.exists(os.path.join(OUT, f)):
                    os.remove(os.path.join(OUT, f))
            pac_run(m, "ac")
            log, _ = pac_run(m, "td")
            a, t = rdc(os.path.join(OUT, f"pac_ac_{m}.txt")), rdc(os.path.join(OUT, f"pac_td_{m}.txt"))
        except CouldNotRun as ex:
            record(f"td PAC vs .ac: {m}", "COULD-NOT-RUN", str(ex))
            continue
        if len(a) != len(t):
            record(f"td PAC vs .ac: {m}", "FAIL", "point count differs")
            continue
        if m in ("aswitch", "pswitch"):
            # plain .ac of the XSPICE conductance-port switches is itself wrong here (1.111 instead of the divider 1k/(1k+r_on)=0.90909);
            # the td result is compared with the exact divider instead
            w = max(abs(y[1] - 1000.0 / 1100.0) / (1000.0 / 1100.0) for y in t)
            parts.append(f"{m} {w:.1e} (vs divider; .ac says {abs(a[0][1]):.4f})")
        else:
            w = max(abs(x[1] - y[1]) / abs(x[1]) for x, y in zip(a, t))
            parts.append(f"{m} {w:.1e}")
        worst = max(worst, w)
    record("td .pac vs plain .ac: stateless analog code models", "pass" if worst <= 1e-6 and len(parts) == 10 else "FAIL",
           f"worst rel {worst:.2e} (bound 1e-6; s103 gave 0.04..0.57 for gain/summer/mult/divide); " + ", ".join(parts))


def check_pac_refusals():
    cases = [("oneshot", "keeps integration state"), ("square", "keeps integration state")]   # card #557: s_xfer int d_dt slew hyst now have a td model
    ok_all = True
    parts = []
    for m, why in cases:
        try:
            log, rc = pac_run(m, "td")
        except CouldNotRun as ex:
            record(f"td PAC refuses {m}", "COULD-NOT-RUN", str(ex))
            continue
        mm = re.search(r"REFUSED -- XSPICE code model '(\w+)'", log)
        ok = bool(mm) and why in log
        ok_all &= ok
        parts.append(f"{m}: {'refused' if ok else 'NOT refused'}")
    # a digital deck: .pac solver=td on the divider deck
    try:
        deck = fixture("xspice_dig_tff.cir") + ".options reltol=1e-5\n.pac 250k 40u out 400 1 8 1e-6 dec 3 1k 10k method=\"shoot\" solver=\"td\" driven=1\n.control\nrun\nquit\n.endc\n.end\n"
        log, rc = run(deck, "pac_dig")
    except CouldNotRun as ex:
        record("td PAC refuses digital", "COULD-NOT-RUN", str(ex))
        return
    ok = "PAC-TD: done" in log and "REFUSED" not in log and "event-driven XSPICE blocks (card #714)" in log
    ok_all &= ok
    parts.append(f"adc_bridge/d_tff/dac_bridge: {'runs (card #714)' if ok else 'did NOT run'}")
    try:
        logr, rc = run(deck, "pac_dig_refuse", env={"PAC_DIG_REFUSE": "1"})
        okr = re.search(r"REFUSED -- XSPICE code model '(adc_bridge|d_tff|dac_bridge)'", logr) is not None and "event-driven" in logr
    except CouldNotRun:
        okr = False
    ok_all &= okr
    parts.append(f"PAC_DIG_REFUSE=1: {'refused' if okr else 'NOT refused'}")
    record("td .pac/.pstb/.pxf/.pnoise REFUSE the remaining stateful analog (oneshot, clock sources), naming the model; the digital divider runs, PAC_DIG_REFUSE=1 refuses it", "pass" if ok_all else "FAIL", "; ".join(parts))


# --- card #557: td small-signal through stateful analog code models ----------------------------------------------------------
TD_MODELS = {
    "s_xfer": (".model m s_xfer(num_coeff=[1] den_coeff=[2.5e-12 2.2e-6 1] int_ic=[0 0])", "A1 in o1 m", "SIN(0.3 0.2 1k)"),
    "int": (".model m int(in_offset=0.0 gain=1e5 out_lower_limit=-1e9 out_upper_limit=1e9 limit_range=1e-3 out_ic=0)\n.model sm summer(in_gain=[1 -1])",
            "A0 [in o1] e1 sm\nA1 e1 o1 m", "SIN(0.3 0.2 1k)"),
    "d_dt": (".model m d_dt(gain=1e-6 out_offset=0 out_lower_limit=-1e9 out_upper_limit=1e9 limit_range=1e-3)", "A1 in o1 m", "SIN(0.3 0.2 1k)"),
    "slew": (".model m slew(rise_slope=1e9 fall_slope=1e9)", "A1 in o1 m", "SIN(0.5 0.05 1k)"),
    "hyst": (".model m hyst(in_low=0 in_high=1 hyst=0.1 out_lower_limit=0 out_upper_limit=1 input_domain=0.05 fraction=FALSE)", "A1 in o1 m",
             "SIN(0.5 0.05 1k)"),
}


def td_disc(model, f, h):
    """exact transfer of the backward-Euler chain cm_analog_integrate() builds (uniform step h)"""
    z = complex(math.cos(2 * math.pi * f * h), math.sin(2 * math.pi * f * h))
    zi = 1 / z
    if model == "int":
        return 1e5 * h / (1 - zi + 1e5 * h)
    if model == "d_dt":
        return 1e-6 * (1 - zi) / h
    if model == "s_xfer":
        G, d0, d1 = 1 / 2.5e-12, 1.0 / 2.5e-12, 2.2e-6 / 2.5e-12
        k = 1 - zi
        return h * h * (G / (1 + d0 * zi * h * h / k ** 2 + d1 * zi * h / k)) / k ** 2
    return 1.0


def td_cont(model, f):
    s = complex(0, 2 * math.pi * f)
    if model == "int":
        return 1e5 / (s + 1e5)
    if model == "d_dt":
        return 1e-6 * s
    if model == "s_xfer":
        return 1 / (2.5e-12 * s * s + 2.2e-6 * s + 1)
    return 1.0


def td_deck(model, kind, pp, extra="", amp=None, card="pac", pxf_src=False):
    mdl, inst, src = TD_MODELS[model]
    if amp is not None:
        src = re.sub(r"SIN\(([\d.]+) [\d.]+ 1k\)", lambda m: f"SIN({m.group(1)} {amp} 1k)", src)
    d = (f"* stateful xspice {model}\nVin in 0 {src} AC 1\nRdum in dum 1meg\nCdum dum 0 1p\n{inst}\n{mdl}\nRl o1 0 1k\n"
         ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n")
    wr = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata %s/td_%s_%s_%s.txt %s\nquit\n.endc\n.end\n"
    tag = f"{model}_{kind}_{pp}{extra and '_' + re.sub(r'[^a-z0-9]', '', extra.lower())}"
    if kind == "ac":
        return d + ".ac dec 5 1k 100k\n" + wr % (OUT, model, "ac", 0, "o1"), tag
    return (d + f'.pac 1k 0 o1 2000 1 8 1e-6 dec 5 1k 100k method="shoot" solver="td" pacpoints={pp}\n' + wr % (OUT, model, "td", pp, "o1"), tag)


def td_run(model, kind, pp, env=None):
    deck, tag = td_deck(model, kind, pp)
    out = os.path.join(OUT, f"td_{model}_{'ac' if kind == 'ac' else 'td'}_{0 if kind == 'ac' else pp}.txt")
    if os.path.exists(out):
        os.remove(out)
    log, rc = run(deck, "td_" + tag + ("_env" if env else ""), env=env)
    return rdc(out), log


def check_td_stateful_lti():
    parts, worst_d, worst_r, worst_ac = [], 0.0, 0.0, 0.0
    ok_n = 0
    for m in ("s_xfer", "int", "d_dt", "slew", "hyst"):
        try:
            ac, _ = td_run(m, "ac", 0)
            t1, log1 = td_run(m, "td", 100000)
            if m in ("s_xfer", "int", "d_dt"):
                t2, _ = td_run(m, "td", 200000)
        except CouldNotRun as ex:
            record(f"td .pac through stateful {m}", "COULD-NOT-RUN", str(ex))
            continue
        h = 1e-3 / 100000
        wd = max(abs(t[1] - td_disc(m, t[0], h)) / abs(td_disc(m, t[0], h)) for t in t1)
        wa = max(abs(a[1] - td_cont(m, a[0])) / abs(td_cont(m, a[0])) for a in ac)
        worst_d, worst_ac = max(worst_d, wd), max(worst_ac, wa)
        txt = f"{m}: vs BE-discrete {wd:.1e}, .ac vs hand {wa:.1e}"
        if m in ("s_xfer", "int", "d_dt"):
            rich = max(abs(2 * y[1] - x[1] - a[1]) / abs(a[1]) for x, y, a in zip(t1, t2, ac))
            worst_r = max(worst_r, rich)
            txt += f", Richardson(1e5,2e5 pts) vs .ac {rich:.1e}, hidden states: " + (re.search(r"(\d+) hidden states", log1) or [None, "?"])[1]
        else:
            txt += f", td vs .ac {max(abs(a[1] - t[1]) / abs(a[1]) for a, t in zip(ac, t1)):.1e}"
        parts.append(txt)
        ok_n += 1
    record("td .pac THROUGH stateful analog code models s_xfer / int / d_dt / slew / hyst (periodic LTI orbit, sine source) = exact BE-discrete transfer "
           "and, extrapolated in h, plain .ac",
           "pass" if ok_n == 5 and worst_d <= 2e-5 and worst_r <= 1e-4 and worst_ac <= 1e-9 else "FAIL",
           f"bounds: discrete 2e-5, Richardson 1e-4, .ac vs hand 1e-9; worst {worst_d:.1e} / {worst_r:.1e} / {worst_ac:.1e}; " + "; ".join(parts))


def check_td_stateful_red():
    try:
        ref, _ = td_run("s_xfer", "td", 100000)
        h = 1e-3 / 100000
        out = {}
        for env in ({"PAC_HID_RED": "1"}, {"PAC_HID_NO_JCOL": "1"}):
            r, _ = td_run("s_xfer", "td", 100000, env=env)
            out[list(env)[0]] = max(abs(t[1] - td_disc("s_xfer", t[0], h)) / abs(td_disc("s_xfer", t[0], h)) for t in r)
    except CouldNotRun as ex:
        record("red proofs PAC_HID_RED / PAC_HID_NO_JCOL", "COULD-NOT-RUN", str(ex))
        return
    good = max(abs(t[1] - td_disc("s_xfer", t[0], h)) / abs(td_disc("s_xfer", t[0], h)) for t in ref)
    record("red proofs: without the hidden states (PAC_HID_RED=1) or with the declared PARTIALs (PAC_HID_NO_JCOL=1) the s_xfer check must fail its 2e-5 bound",
           "pass" if good <= 2e-5 and out["PAC_HID_RED"] > 1e-2 and out["PAC_HID_NO_JCOL"] > 1e-3 else "FAIL",
           f"normal {good:.1e}; PAC_HID_RED {out['PAC_HID_RED']:.1e} (must exceed 1e-2); PAC_HID_NO_JCOL {out['PAC_HID_NO_JCOL']:.1e} (must exceed 1e-3)")


def td_other(model, tail, name, env=None, amp=None, src=None):
    mdl, inst, s0 = TD_MODELS[model]
    if src:
        s0 = src
    deck = (f"* stateful xspice {model}\nVin in 0 {s0} AC 1\nRdum in dum 1meg\nCdum dum 0 1p\n{inst}\n{mdl}\nRl o1 0 1k\n"
            ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n" + tail)
    return run(deck, name, env=env)


def check_td_stateful_refusals_and_pxf():
    parts, ok_all = [], True
    opts = 'method="shoot" solver="td" pacpoints=100000'
    wrp = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata %s/%s.txt %s\nquit\n.endc\n.end\n"
    try:
        # pnoise
        log, _ = td_other("s_xfer", f".pnoise 1k 0 o1 2000 1 8 1e-6 o1 Vin lin 3 1k 50k {opts}\n.control\nrun\nquit\n.endc\n.end\n", "td_pnoise")
        ok = "PNOISE-TD: REFUSED" in log
        ok_all &= ok; parts.append(f"pnoise through s_xfer refused: {ok}")
        # gear
        log, _ = td_other("s_xfer", f".options method=gear\n.pac 1k 0 o1 2000 1 8 1e-6 lin 3 1k 50k {opts}\n.control\nrun\nquit\n.endc\n.end\n", "td_gear")
        ok = "REFUSED -- a stateful XSPICE code model" in log and "method=gear" in log
        ok_all &= ok; parts.append(f"method=gear with s_xfer refused: {ok}")
        # slew with slew-rate limiting on the orbit (slews near the zero crossings only: a regime change)
        mdl = TD_MODELS["slew"][0].replace("1e9", "3000")
        deck = ("* slew slewing\nVin in 0 SIN(0 0.5 1k) AC 1\nRdum in dum 1meg\nCdum dum 0 1p\nA1 in o1 m\n" + mdl + "\nRl o1 0 1k\n.options reltol=1e-6\n"
                ".pac 1k 0 o1 2000 1 8 1e-6 lin 3 1k 50k " + opts + "\n.control\nrun\nquit\n.endc\n.end\n")
        log, _ = run(deck, "td_slew_slewing")
        ok = "slew-rate limiting of a 'slew' code model is active" in log and "REFUSED" in log
        ok_all &= ok; parts.append(f"slew that slews on the orbit refused: {ok}")
        # hyst that flips
        mdl = TD_MODELS["hyst"][0]
        deck = ("* hyst flipping\nVin in 0 SIN(0.5 1.5 1k) AC 1\nRdum in dum 1meg\nCdum dum 0 1p\nA1 in o1 m\n" + mdl + "\nRl o1 0 1k\n.options reltol=1e-6\n"
                ".pac 1k 0 o1 2000 1 8 1e-6 lin 3 1k 50k " + opts + "\n.control\nrun\nquit\n.endc\n.end\n")
        log, _ = run(deck, "td_hyst_flip")
        ok = "Boolean state of a 'hyst' code model changes on the orbit" in log
        ok_all &= ok; parts.append(f"hyst that flips refused: {ok}")
        # pxf: adjoint=1 falls back, equals adjoint=0 and .pac
        res = {}
        for tag, tail in (("a1", f".pxf 1k 0 o1 2000 1 8 1e-6 o1 lin 3 1k 50k sources=\"Vin\" adjoint=1 {opts}\n"),
                          ("a0", f".pxf 1k 0 o1 2000 1 8 1e-6 o1 lin 3 1k 50k sources=\"Vin\" adjoint=0 {opts}\n"),
                          ("pac", f".pac 1k 0 o1 2000 1 8 1e-6 lin 3 1k 50k {opts}\n")):
            out = os.path.join(OUT, f"td_pxf_{tag}.txt")
            if os.path.exists(out):
                os.remove(out)
            vec = "v(o1)" if tag == "pac" else "pxf_Vin"
            log, _ = td_other("s_xfer", tail + wrp % (OUT, "td_pxf_" + tag, vec), "td_pxf_" + tag)
            res[tag] = (rdc(out), log)
        fb = "adjoint engine REFUSED" in res["a1"][1] and "does not transpose its hidden states" in res["a1"][1]
        d = max(abs(a[1] - b[1]) / abs(b[1]) for a, b in zip(res["a1"][0], res["pac"][0]))
        d0 = max(abs(a[1] - b[1]) / abs(b[1]) for a, b in zip(res["a0"][0], res["pac"][0]))
        ok = fb and d <= 1e-9 and d0 <= 1e-9
        ok_all &= ok; parts.append(f".pxf adjoint=1 falls back (message {fb}), equals .pac to {d:.1e} / adjoint=0 {d0:.1e}")
    except CouldNotRun as ex:
        record("stateful XSPICE refusals and .pxf", "COULD-NOT-RUN", str(ex))
        return
    record("remaining refusals through stateful XSPICE models (pnoise, method=gear, slew-rate limiting, hyst flips) and .pxf adjoint fallback",
           "pass" if ok_all else "FAIL", "; ".join(parts))


# R3 voltage-mode buck with the compensator replaced by an XSPICE s_xfer, control-to-output (series injection Vctrl).  Reference: a frequency-response
# measurement on the real nonlinear circuit (sine of amplitude 0.01 V in series at Vctrl, .tran reltol 1e-5 maxstep 1 ns, single-bin DFT over the last
# 2 ms = whole numbers of the signal and of the 500 kHz switching periods, no baseline subtraction needed); amplitude 0.0025 repeats it.
# f -> (dB, deg) of vout / Vctrl; ~/Work/data/ngspice_augmentation/557/dev/r3xs/fra.py
R3XS_FRA = {4000.0: (-6.1833, 12.682), 20000.0: (-11.7160, -95.560), 50000.0: (-25.2293, -131.765)}


def check_td_r3_sxfer():
    deck = fixture("xspice_r3_sxfer.cir") + ".options reltol=1e-4\n"
    res = {}
    for pp, env in ((2000, None), (2000, {"PAC_HID_RED": "1"})):
        out = os.path.join(OUT, f"td_r3xs_{pp}{'_red' if env else ''}.txt")
        if os.path.exists(out):
            os.remove(out)
        log, rc = run(deck + f'.pac 500e3 1e-3 vout 200 10 150 1e-6 lin 24 4000 50000 method="shoot" solver="td" pacpoints={pp}\n'
                      f".control\n{WR}run\nwrdata {out} v(vout)\nquit\n.endc\n.end\n", "td_r3xs" + ("_red" if env else ""), env=env)
        rows = dict((r[0], complex(r[1], r[2])) for r in rd(out))
        res[bool(env)] = (rows, log)
    worst_db = worst_deg = 0.0
    parts = []
    rows, log = res[False]
    for f, (db, deg) in R3XS_FRA.items():
        H = rows[f]
        d_db = 20 * math.log10(abs(H)) - db
        d_deg = (math.degrees(math.atan2(H.imag, H.real)) - deg + 180) % 360 - 180
        worst_db, worst_deg = max(worst_db, abs(d_db)), max(worst_deg, abs(d_deg))
        parts.append(f"{f / 1e3:g} kHz {d_db:+.4f} dB / {d_deg:+.3f} deg")
    nh = re.search(r"(\d+) hidden states", log)
    red = max(abs(20 * math.log10(abs(res[True][0][f])) - db) for f, (db, deg) in R3XS_FRA.items())
    record("td .pac through an s_xfer compensator in the R3 buck (periodic closed loop, 2 hidden states) vs FRA on the nonlinear circuit",
           "pass" if worst_db <= 0.03 and worst_deg <= 0.15 and nh and nh.group(1) == "2" and red > 1.0 else "FAIL",
           f"worst {worst_db:.4f} dB / {worst_deg:.3f} deg (bound 0.03 dB / 0.15 deg), pacpoints 2000; deltas vs FRA: " + "; ".join(parts) +
           f"; red proof PAC_HID_RED=1 off by {red:.1f} dB (must exceed 1)")


def check_td_pstb_r3_sxfer():
    """`.pstb` (Tian, two forcing patterns on the same hidden-state Phi) of the R3 + s_xfer loop, probe Vprobe, against the closed-loop identity
    H_cl = G / (1 + T) with H_cl the series-injection control-to-output of check_td_r3_sxfer and G the same transfer with the loop opened (compensator
    replaced by a DC source at the regulated control value).  The opened-loop operating point is not exactly the closed one, so G is only known up to a
    constant factor (measured 1.0926): the check is that G/(1+T)/H_cl is CONSTANT in magnitude (spread <= 2e-3) and drifts <= 2 deg over 4..50 kHz,
    i.e. the frequency dependence of T through the 2 hidden states is that of the loop; a T without them is off by tens of dB."""
    name = "td .pstb through an s_xfer compensator in the R3 buck: G/(1+T)/H_cl constant over 4..50 kHz (closed-loop identity)"
    base = fixture("xspice_r3_sxfer.cir")
    cl = base.replace("Vctrl ctrlnode opout DC 0 AC 1", "Vctrl ctrlnode opout DC 0 AC 1")
    ol = base.replace("A1 e1 opout mc", "* open loop").replace("Rop opout 0 1meg", "Vop opout 0 DC 0.6405").replace("A0 [vrefnode fbnode] e1 sm", "* open loop")
    tail = f".control\n{WR}run\nwrdata %s %s\nquit\n.endc\n.end\n"
    try:
        res = {}
        for tag, deck, card, vec in (
                ("h", cl, '.pac 500e3 1e-3 vout 200 10 150 1e-6 lin 24 4000 50000 method="shoot" solver="td" pacpoints=2000', "v(vout)"),
                ("g", ol, '.pac 500e3 1e-3 vout 200 10 150 1e-6 lin 24 4000 50000 method="shoot" solver="td" pacpoints=2000', "v(vout)"),
                ("t", cl.replace(" AC 1", ""), '.pstb 500e3 1e-3 vout 200 10 150 1e-6 lin 24 4000 50000 probe="Vprobe" method="shoot" pacpoints=2000', "loopgain t")):
            out = os.path.join(OUT, f"td_pstb_{tag}.txt")
            if os.path.exists(out):
                os.remove(out)
            log, rc = run(deck + ".options reltol=1e-4\n" + card + "\n" + tail % (out, vec), "td_pstb_" + tag)
            res[tag] = (dict((r[0], r[1:]) for r in rd(out)), log)
        ok_h = re.search(r"(\d+) hidden states", res["t"][1]) and re.search(r"(\d+) hidden states", res["t"][1]).group(1) == "2"
        ratios, phs = [], []
        for f in (4000.0, 10000.0, 20000.0, 30000.0, 50000.0):
            t = complex(*res["t"][0][f][2:4])
            g = complex(*res["g"][0][f][0:2])
            h = complex(*res["h"][0][f][0:2])
            pred = g / (1 + t)
            ratios.append(abs(pred) / abs(h))
            phs.append(math.degrees(math.atan2((pred / h).imag, (pred / h).real)))
        spread = (max(ratios) - min(ratios)) / (sum(ratios) / len(ratios))
        drift = max(phs) - min(phs)
        record(name, "pass" if ok_h and spread <= 2e-3 and drift <= 2.0 else "FAIL",
               f"|G/(1+T)|/|H_cl| = {ratios[0]:.4f} .. {ratios[-1]:.4f} (spread {spread:.1e}, bound 2e-3), phase of the ratio {min(phs):+.2f} .. {max(phs):+.2f} deg (drift {drift:.2f}, "
               f"bound 2), hidden states 2: {bool(ok_h)}")
    except (CouldNotRun, KeyError) as ex:
        record(name, "COULD-NOT-RUN", str(ex))


def check_asan_stateful():
    """card #557: ASan twin (NGSPICE_ASAN) on the new hidden-state paths: R3 + s_xfer `.pac` and `.pstb`, the LTI s_xfer `.pac`, `.pxf` (adjoint=1 fallback)."""
    name = "ASan: td .pac/.pstb/.pxf through stateful analog code models (R3 + s_xfer, LTI s_xfer)"
    asan = os.environ.get("NGSPICE_ASAN") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-557-asan/bin/ngspice")
    if not os.path.exists(asan):
        record(name, "COULD-NOT-RUN", f"{asan} not found (set NGSPICE_ASAN=)")
        return
    global NG
    saved = NG
    NG = asan
    try:
        base = fixture("xspice_r3_sxfer.cir") + ".options reltol=1e-4\n"
        wr = f".control\n{WR}run\nquit\n.endc\n.end\n"
        logs = []
        for tag, card in (("pac", '.pac 500e3 1e-3 vout 200 10 150 1e-6 lin 6 4000 50000 method="shoot" solver="td" pacpoints=1000'),
                          ("pstb", '.pstb 500e3 1e-3 vout 200 10 150 1e-6 lin 6 4000 50000 probe="Vprobe" method="shoot" pacpoints=1000')):
            deck = base.replace(" AC 1", "") if tag == "pstb" else base
            log, rc = run(deck + card + "\n" + wr, "td_asan_r3_" + tag, env={"ASAN_OPTIONS": "detect_leaks=0"})
            logs.append((tag, log, rc))
        mdl, inst, src = TD_MODELS["s_xfer"]
        deck = (f"* asan s_xfer\nVin in 0 {src} AC 1\nRdum in dum 1meg\nCdum dum 0 1p\n{inst}\n{mdl}\nRl o1 0 1k\n.options reltol=1e-6\n"
                '.pxf 1k 0 o1 2000 1 8 1e-6 o1 lin 3 1k 50k sources="Vin" adjoint=1 method="shoot" solver="td" pacpoints=20000\n' + wr)
        log, rc = run(deck, "td_asan_pxf", env={"ASAN_OPTIONS": "detect_leaks=0"})
        logs.append(("pxf", log, rc))
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    finally:
        NG = saved
    hit = [t for t, l, r in logs if "ERROR: AddressSanitizer" in l or "runtime error" in l]
    done = all(("PAC-TD: done" in l) if t == "pac" else ("PSTB: done" in l) if t == "pstb" else ("PXF-TD: done" in l) for t, l, r in logs)
    record(name, "pass" if not hit and done else "FAIL", f"ASan hits {hit or 'none'}; all three analyses completed: {done}")


def check_pss_refusals():
    parts, ok_all = [], True
    # unvalidated model by name
    deck = ("* t\nVin in 0 SIN(0 2 100k)\nRs in o1 100\nA1 %gd(o1 0) m\n.model m zener(v_breakdown=1.5 i_breakdown=1m r_breakdown=10 i_sat=1e-12 n_forward=1)\n"
            "Co o1 0 100p\n.options reltol=1e-5\n.control\npss 100e3 40e-6 o1 400 10 150 1e-6 uic method=shoot\nquit\n.endc\n.end\n")
    try:
        log, rc = run(deck, "ref_zener")
        ok = "REFUSED -- XSPICE code model 'zener'" in log and "CONVERGED" not in log
        ok_all &= ok
        parts.append(f"zener refused by name: {ok}")
        # non-periodic digital state
        for f, want in ((1e6, False), (5e5, True)):
            e, log, conv = pss_vs_tran("xspice_dig_toggle.cir", f"tog_{int(f)}", f, vecs="v(out)", pts=4000, tstab=40e-6, tstop=240e-6,
                                       extra="driven=1", step=0.25e-9)
            ok = (conv == want) and (want or "not periodic with the analysis period" in log) and (not want or e[0] <= 0.4)
            ok_all &= ok
            parts.append(f"d_tff toggling, PSS at {f / 1e6:g} MHz {'converges (%.3f %%)' % e[0] if want and conv else 'refused' if not conv else 'CONVERGED?'}: {ok}")
        # free-running with XSPICE
        deck = fixture("xspice_dig_tff.cir") + ".options reltol=1e-5\n.control\npss 250k 40u out 400 10 150 1e-6 uic method=shoot trigger=out trigval=0.5 trigdir=rise\nquit\n.endc\n.end\n"
        log, rc = run(deck, "ref_freerun")
        ok = "free-running (trigger=)" in log and "CONVERGED" not in log
        ok_all &= ok
        parts.append(f"free-running with XSPICE refused: {ok}")
    except CouldNotRun as ex:
        record("PSS refusals", "COULD-NOT-RUN", str(ex))
        return
    record("PSS refuses unvalidated models, non-periodic digital state and free-running XSPICE decks", "pass" if ok_all else "FAIL", "; ".join(parts))


def check_red_proofs():
    try:
        e, log, ok = pss_vs_tran("xspice_dig_dsample.cir", "red_dsample", 250e3, vecs="v(out)", pts=4000, tstab=40e-6, tstop=240e-6,
                                 extra="driven=1", step=0.25e-9, env={"PSS_XS_NO_STATE": "1"})
        e2, log2, ok2 = pss_vs_tran("xspice_an_s_xfer.cir", "red_sxfer", 100e3, vecs="v(out) v(o1)", step=2e-9, env={"PSS_XS_NO_STATE": "1"}, cols=(1, 2))
    except CouldNotRun as ex:
        record("red proof PSS_XS_NO_STATE=1", "COULD-NOT-RUN", str(ex))
        return
    d1 = e[0] if ok else float("nan")
    d2 = max(e2) if ok2 else float("nan")
    record("red proof PSS_XS_NO_STATE=1 (restart without the carried XSPICE state) must exceed the bounds",
           "pass" if (not ok or d1 > 1.0) and (not ok2 or d2 > 1.5) else "FAIL",
           f"held-state deck (bound 0.4 %): {d1:.3f} %; s_xfer (bound 1.5 %): {d2:.3f} %")


# --- R9 ----------------------------------------------------------------------------------------------------------------
def r9_pss(fx, tag, f, tstab, pts, env=None, timeout=3000, sciter=150):
    base = fixture(fx) + ".options reltol=1e-4\n"
    pth = os.path.join(OUT, f"{tag}_pss.txt")
    if os.path.exists(pth):
        os.remove(pth)
    t0 = time.time()
    log, _ = run(base + ".control\n" + WR.replace("numdgt=15", "numdgt=12") +
                 f"pss {f:.15g} {tstab:g} vout {pts} 10 {sciter} 1e-6 uic method=shoot\nwrdata {pth} v(vout) v(code)\nquit\n.endc\n.end\n",
                 f"{tag}_pss", env=env, timeout=timeout)
    return log, time.time() - t0


def check_r9():
    T = 1.2e-3
    f = 1.0 / T
    try:
        log, w = r9_pss("r9_fine.cir", "r9_fine", f, 12e-3, 24000)
        if "PSS-SHOOT: CONVERGED" not in log:
            record("R9 PSS", "FAIL", "did not converge: " + " | ".join(l for l in log.split("\n") if "FAILED" in l)[:300])
            return
        run(fixture("r9_fine.cir") + ".options reltol=1e-4\n.control\n" + WR.replace("numdgt=15", "numdgt=12") +
            f"tran 100n 60m uic\nwrdata {OUT}/r9_fine_tran.txt v(vout) v(code)\nquit\n.endc\n.end\n", "r9_fine_tran", timeout=3000)
        P, R = rd(os.path.join(OUT, "r9_fine_pss.txt")), rd(os.path.join(OUT, "r9_fine_tran.txt"))
    except CouldNotRun as ex:
        record("R9 PSS vs settled tran", "COULD-NOT-RUN", str(ex))
        return
    rt = [r[0] for r in R]
    rv = [r[1] for r in R]
    t0 = rt[-1] - T
    pt = [0.0] + [r[0] for r in P]
    pv = [P[0][1]] + [r[1] for r in P]
    N = 3000
    c = [interp(rt, rv, t0 + T * j / N) for j in range(N + 1)]
    p = [interp(pt, pv, pt[-1] * j / N) for j in range(N + 1)]
    rip = max(c) - min(c)
    d = max(abs(a - b) for a, b in zip(c, p))
    mo = sum(a - b for a, b in zip(p, c)) / len(c)
    codes = sorted(set(round(r[2] * 2048) for r in P))
    m = re.search(r"digital node mismatches (\d+), analog node drift (\S+)", log)
    record("R9 (digital buck, 11-bit DPWM) PSS orbit vs settled tran (period 1.2 ms = 600 switching cycles)",
           "pass" if d / rip <= 0.01 and abs(mo) <= 0.1 * rip and codes == [556] else "FAIL",
           f"ripple {rip * 1e3:.3f} mV, max |d| {d * 1e3:.4f} mV ({100 * d / rip:.2f} % of ripple; bound 1 %), mean offset {mo * 1e3:+.4f} mV, "
           f"DPWM code {codes} (settled tran: 556), {m.group(0) if m else ''}, PSS wall {w:.0f} s")
    # red proof
    try:
        log2, w2 = r9_pss("r9_fine.cir", "r9_ns", f, 6e-3, 6000, env={"PSS_XS_NO_STATE": "1"})
        P2 = rd(os.path.join(OUT, "r9_ns_pss.txt"))
    except CouldNotRun as ex:
        record("R9 red proof", "COULD-NOT-RUN", str(ex))
        return
    mean2 = sum(r[1] for r in P2) / len(P2)
    mean1 = sum(r[1] for r in P) / len(P)
    codes2 = sorted(set(round(r[2] * 2048) for r in P2))
    record("R9 red proof PSS_XS_NO_STATE=1: the orbit of the restarted map (DPWM register back at its initial code) is off by more than the ripple",
           "pass" if abs(mean2 - mean1) > rip else "FAIL", f"mean Vout {mean2:.5f} vs {mean1:.5f} V ({(mean2 - mean1) * 1e3:+.2f} mV, ripple {rip * 1e3:.2f} mV), codes {codes2}")
    # limit-cycling variant
    try:
        log3, w3 = r9_pss("r9_lc.cir", "r9_lc", f, 12e-3, 6000, sciter=12)
        run(fixture("r9_lc.cir") + ".options reltol=1e-4\n.control\n" + WR.replace("numdgt=15", "numdgt=12") +
            f"tran 100n 30m uic\nwrdata {OUT}/r9_lc_tran.txt v(vout) v(code)\nquit\n.endc\n.end\n", "r9_lc_tran", timeout=3000)
        R = rd(os.path.join(OUT, "r9_lc_tran.txt"))
    except CouldNotRun as ex:
        record("R9 limit-cycling variant", "COULD-NOT-RUN", str(ex))
        return
    ok3 = "CONVERGED" not in log3 and ("not periodic with the analysis period" in log3 or "FAILED to converge" in log3)
    wins = {}
    for r in R:
        if r[0] > 12e-3:
            wins.setdefault(int(r[0] / T), []).append(r)
    means = [sum(x[1] for x in v) / len(v) for k, v in sorted(wins.items())]
    codes = sorted(set(round(r[2] * 128) for r in R if r[0] > 12e-3))
    alt = max(abs(means[i] - means[i + 1]) for i in range(len(means) - 1)) if len(means) > 2 else 0.0
    ok4 = alt > 0.05 and len(codes) >= 2
    record("R9 limit-cycling variant (7-bit DPWM: LSB 92 mV >> the 6 mV zero-error bin): the settled tran hunts with period 2T and PSS at T is REFUSED",
           "pass" if ok3 and ok4 else "FAIL",
           f"tran: window-mean Vout alternates by {alt * 1e3:.0f} mV between consecutive {T * 1e3:g} ms windows (bound 50 mV), codes {codes}; "
           f"PSS at T: {'refused' if ok3 else 'NOT refused'} ({w3:.0f} s)")


def main():
    if not os.path.exists(NG):
        print("ngspice build not found:", NG)
        return 2
    for fn in (() if R9ONLY else (check_pac_stateless, check_td_stateful_lti, check_td_stateful_red, check_td_stateful_refusals_and_pxf, check_td_r3_sxfer, check_td_pstb_r3_sxfer, check_asan_stateful, check_pac_refusals, check_stateless_pss, check_stateful_pss, check_digital_pss,
                                  check_pss_refusals, check_red_proofs)):
        t0 = time.time()
        try:
            fn()
        except CouldNotRun as ex:
            record(fn.__name__, "COULD-NOT-RUN", str(ex))
        if not QUICK or fn in (check_pac_stateless, check_pac_refusals, check_td_stateful_lti):
            print(f"   ({fn.__name__}: {time.time() - t0:.0f} s)", flush=True)
    if not NOR9:
        try:
            check_r9()
        except CouldNotRun as ex:
            record("R9", "COULD-NOT-RUN", str(ex))
    nf = sum(1 for r in results if r[1] == "FAIL")
    nc = sum(1 for r in results if r[1] == "COULD-NOT-RUN")
    print(f"\n{len(results) - nf - nc}/{len(results)} pass, {nf} FAIL, {nc} COULD-NOT-RUN")
    return 1 if nf else (2 if nc else 0)


if __name__ == "__main__":
    sys.exit(main())
