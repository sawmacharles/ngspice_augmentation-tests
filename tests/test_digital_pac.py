#!/usr/bin/env python3
"""td small signal through digital XSPICE blocks -- card #714.

    NGSPICE=<prefix>/bin/ngspice [NGSPICE_BASE=<s149 prefix>/bin/ngspice] [DIG714_OUT=<dir>] python3 tests/test_digital_pac.py [--quick]

Exit 0 pass, 1 fail, 2 could-not-run.  Outputs under ~/Work/data/ngspice_augmentation/714/test (DIG714_OUT overrides), never the repo.
Every ngspice runs under systemd-run (MemoryMax=3G, no swap), one at a time, `set num_threads=2`.  About 3 minutes (--quick: 1.5).

What is linearised (docs/pss-pac.md, card #714): a digital edge whose time is set by an adc_bridge threshold crossing moves by
dtau = -dv/v' (v the comparator input, located exactly in the periodic runs), every digital block with fixed delays passes the shift of the
input edge that triggered an output on to that output, a d_pwm transition moves by the shift of its previous transition +- T d'(c) dc (that
shift is a Floquet state: it carries over the period), and a moved dac_bridge edge moves the analog ramp after it ((s_before - s_after)
dtau injected on the dac output while it ramps).  Edges of fixed time (clocks) are fixed waveforms.

References: tests/fixtures/digital714/refs.json, brute-force FRA (plain .tran, uic, one sine added to a DC source or to the lgprobe Vi / Ii
sources, phasor of the difference to the base run by exact integration of the piecewise-linear waveform over whole periods of f and of the
1/fs switching period; adc_bridge decks with XS_ADC_LOCATE=1 so the plain transient places the comparator edges at the crossings, decks
without adc_bridge on the unmodified s149 binary), with an amplitude-halving and a tmax-halving variant pinned next to each value: the
bounds below are 10x or more the spread of those variants and of the measured td-vs-FRA difference (see the per-check details).

Checks
  1  open loop, comparator PWM (a_ol: PULSE ramp vs control -> adc_bridge -> d_buffer -> dac_bridge -> tanh switches): td .pac 2/25/50 kHz
  2  open loop, latch PWM (a_ff: d_osc sets a d_dff, the comparator resets it asynchronously, two d_inverters)
  3  the same with IDEAL sw switches driven by the dac output (a_ffsw: the dac ramp moves the switch event, saltation reads it)
  4  open loop, d_pwm with an analog duty (b_ol2), control-to-output
  5  open loop, clock-only digital path (c_ol: d_pwm at fixed duty), line-to-output
  6  closed loop, comparator PWM (a_cl, type-III compensator): .pstb at 2/10/20/31.25/35.7/50 kHz vs the FRA loop gain (card done-when),
     and the phase margin / crossover printed by .pstb vs the FRA interpolated between 31.25 and 35.7 kHz
  7  closed loop, d_osc-reset ramp, comparator with a +-1 mV band (two crossings per edge, dac through out_undef) (a_cl2): .pstb
  8  closed loop, clock-only digital path, analog comparator (c_cl): .pstb
  9  a_cl .pss: converges, orbit averages = the settled located .tran (pinned)
 10  .pxf on a_ol: forward Vc column == .pac; adjoint=1 falls back to the forward engine with a message
 11  refusals naming the instance: a tie (r_tie), a multi-bit adc (r_mbit), a VCO (r_vco), a moved ramp in flight at t = T (a_ol with
     t_rise 20 ns), a d_pwm transition whose dac ramp crosses t = T (b_ol2 with init_phase 0), .pnoise
 12  red proofs: PAC_DIG_REFUSE=1 (the s149 refusal), PAC_DIG_NOSHIFT=1 (no edge shift: zero modulator gain, closed loop open),
     PAC_DIG_PWM_NOSTATE=1 (d_pwm shift without its Floquet state), PSS_ADC_LOCATE=0 (closed-loop .pss does not converge)
 13  (NGSPICE_BASE set) the baseline refuses 1 and does not converge 9
"""
import cmath
import json
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIX = os.path.join(HERE, "fixtures", "digital714")
NG = os.environ.get("NGSPICE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice")
BASE = os.environ.get("NGSPICE_BASE")
OUT = os.environ.get("DIG714_OUT") or os.path.expanduser("~/Work/data/ngspice_augmentation/714/test")
QUICK = "--quick" in sys.argv
results = []
WR = "set num_threads=2\nset filetype=ascii\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
PTS = 4000


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"[{state:13s}] {name}: {detail}", flush=True)


def run(deck, name, env=None, ng=None, timeout=1800):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name + ".cir")
    with open(p, "w") as fh:
        fh.write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    for k in ("PAC_DIG_REFUSE", "PAC_DIG_NOSHIFT", "PAC_DIG_PWM_NOSTATE", "PSS_ADC_LOCATE", "XS_ADC_LOCATE", "PAC_DIG_COARSE"):
        e.pop(k, None)
    if env:
        e.update(env)
    try:
        r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=3G", "-p", "MemorySwapMax=0", ng or NG, "-b", p],
                           capture_output=True, text=True, cwd=OUT, env=e, timeout=timeout)
    except (FileNotFoundError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{name}: {ex}")
    out = (r.stdout + r.stderr).replace("\r", "\n")
    with open(os.path.join(OUT, name + ".log"), "w") as fh:
        fh.write(out)
    return out


def fixture(name):
    with open(os.path.join(FIX, name + ".cir")) as fh:
        txt = fh.read()
    txt = txt.replace("../../../stb/lgprobe.lib", os.path.join(ROOT, "stb", "lgprobe.lib"))
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


def refused(log):
    return [l.strip() for l in log.split("\n") if "REFUSED" in l]


def td_pac(deck, src, f, name, env=None, extra=""):
    """complex v(vout) per unit AC on src, or (None, log)"""
    d = re.sub(r"^(%s\s+\S+\s+\S+\s+.*)$" % re.escape(src), lambda m: m.group(1) + " AC 1", deck, count=1, flags=re.M | re.I)
    fn = os.path.join(OUT, name + ".txt")
    if os.path.exists(fn):
        os.replace(fn, fn + ".old")
    log = run(d + f".options reltol=1e-5\n.pac 500k 100u vout {PTS} 1 50 1e-6 lin 1 {f} {f} method=\"shoot\" solver=\"td\" {extra}\n"
              f".control\n{WR}run\nwrdata {fn} v(vout)\nquit\n.endc\n.end\n", name, env)
    if refused(log) or not os.path.exists(fn) or "PAC-TD: done" not in log:
        return None, log
    r = rd(fn)[0]
    return complex(r[1], r[2]), log


def td_pstb(deck, freqs, name, env=None, sweep=None):
    fn = os.path.join(OUT, name + ".txt")
    if os.path.exists(fn):
        os.replace(fn, fn + ".old")
    sw = sweep or f"lin 1 {freqs[0]} {freqs[0]}"
    log = run(deck + f".options reltol=1e-5\n.pstb 500k 1m vout {PTS} 10 150 1e-6 {sw} method=\"shoot\" probe=\"XP\"\n"
              f".control\n{WR}run\nwrdata {fn} T\nquit\n.endc\n.end\n", name, env)
    if refused(log) or not os.path.exists(fn) or "PSTB: done" not in log:
        return None, log
    return [(r[0], complex(r[1], r[2])) for r in rd(fn)], log


def db(z):
    return 20 * math.log10(abs(z)) if abs(z) > 0 else -999.0


def dphase(a, b):
    return (math.degrees(cmath.phase(a / b)) + 180.0) % 360.0 - 180.0


REFS = None


def refs():
    global REFS
    if REFS is None:
        with open(os.path.join(FIX, "refs.json")) as fh:
            REFS = json.load(fh)
    return REFS


def ref_tf(deck, f):
    v = refs()[deck]["variants"]
    a = refs()[deck]["amps"][0]
    main = complex(*v["main"]["H"]["%g/%g" % (f, a)])
    spread = 0.0
    for key, var in v.items():
        for k, z in var.get("H", {}).items():
            if k.split("/")[0] == "%g" % f:
                spread = max(spread, abs(db(complex(*z)) - db(main)))
    return main, spread


def ref_loop(deck, f):
    v = refs()[deck]["variants"]
    main = complex(*v["main"]["T"]["%g" % f])
    spread = 0.0
    for var in v.values():
        z = var.get("T", {}).get("%g" % f)
        if z:
            spread = max(spread, abs(db(complex(*z)) - db(main)))
    return main, spread


def check_tf(label, deck, src, bound_db, bound_deg, freqs=(2000.0, 25000.0, 50000.0)):
    try:
        d = fixture(deck)
        parts, ok = [], True
        for f in freqs:
            z, log = td_pac(d, src, f, f"{deck}_pac_{f:g}")
            if z is None:
                record(label, "FAIL", f"td .pac at {f:g} Hz: " + (" | ".join(refused(log))[:300] or log[-300:]))
                return
            r, spread = ref_tf(deck, f)
            ddb, ddeg = db(z) - db(r), dphase(z, r)
            ok = ok and abs(ddb) <= bound_db and abs(ddeg) <= bound_deg
            parts.append(f"{f / 1e3:g}k {db(z):.4f} dB/{math.degrees(cmath.phase(z)):.3f} deg vs FRA {db(r):.4f}/{math.degrees(cmath.phase(r)):.3f} "
                         f"({ddb:+.4f} dB {ddeg:+.4f} deg; FRA variants {spread:.4f} dB)")
        record(label, "pass" if ok else "FAIL", f"bound {bound_db} dB / {bound_deg} deg: " + "; ".join(parts))
    except (CouldNotRun, KeyError, OSError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])


def check_loop(label, deck, freqs, bound_db, bound_deg):
    try:
        d = fixture(deck)
        parts, ok = [], True
        for f in freqs:
            res, log = td_pstb(d, [f], f"{deck}_pstb_{f:g}")
            if res is None:
                record(label, "FAIL", f".pstb at {f:g} Hz: " + (" | ".join(refused(log))[:300] or log[-300:]))
                return
            z = res[0][1]
            r, spread = ref_loop(deck, f)
            ddb, ddeg = db(z) - db(r), dphase(z, r)
            ok = ok and abs(ddb) <= bound_db and abs(ddeg) <= bound_deg
            parts.append(f"{f / 1e3:.4g}k {db(z):.4f} dB/{math.degrees(cmath.phase(z)):.3f} deg ({ddb:+.4f} dB {ddeg:+.4f} deg; FRA variants "
                         f"{spread:.4f} dB)")
        record(label, "pass" if ok else "FAIL", f"vs FRA loop gain, bound {bound_db} dB / {bound_deg} deg: " + "; ".join(parts))
    except (CouldNotRun, KeyError, OSError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])


def check_margin():
    label = "6b. closed loop a_cl: .pstb phase margin / crossover vs the FRA (interpolated in dB and deg against ln f, 31.25..35.7 kHz)"
    try:
        res, log = td_pstb(fixture("a_cl"), [20000.0], "a_cl_margin", sweep="lin 13 20k 50k")
        m = re.search(r"phase margin ([-0-9.]+) deg at ([0-9.]+) Hz", log)
        if res is None or not m:
            record(label, "FAIL", " | ".join(refused(log))[:300] or "no phase margin line")
            return
        pm, fc = float(m.group(1)), float(m.group(2))
        f1, f2 = 31250.0, 500e3 / 14
        (z1, _), (z2, _) = ref_loop("a_cl", f1), ref_loop("a_cl", f2)
        x = db(z1) / (db(z1) - db(z2))
        fcr = f1 * math.exp(x * math.log(f2 / f1))
        p1, p2 = math.degrees(cmath.phase(z1)), math.degrees(cmath.phase(z2))
        pmr = 180.0 + p1 + x * (p2 - p1)
        ok = abs(pm - pmr) <= 0.1 and abs(fc / fcr - 1) <= 0.005
        record(label, "pass" if ok else "FAIL", f".pstb PM {pm:.3f} deg at {fc:.1f} Hz; FRA {pmr:.3f} deg at {fcr:.1f} Hz (bounds 0.1 deg, 0.5 %)")
    except (CouldNotRun, KeyError, OSError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])


def check_orbit():
    label = "9. closed loop a_cl: .pss converges (adc_bridge crossings located); orbit averages = settled located .tran (pinned)"
    try:
        fn = os.path.join(OUT, "a_cl_pss.txt")
        log = run(fixture("a_cl") + f".options reltol=1e-5\n.control\n{WR}pss 500k 1m vout {PTS} 10 150 1e-6 uic method=shoot\n"
                  f"wrdata {fn} v(vout) v(opout)\nquit\n.endc\n.end\n", "a_cl_pss")
        m = re.search(r"PSS-SHOOT: CONVERGED iter=(\d+) residual=(\S+) spectral_radius=(\S+)", log)
        if not m:
            record(label, "FAIL", " | ".join(l for l in log.split("\n") if "FAILED" in l or "REFUSED" in l)[:300])
            return
        P = rd(fn)
        ts = [r[0] for r in P]

        def uavg(col, n=20000):   # average on a uniform grid over the orbit (the output points are not uniform)
            import bisect
            vs = [r[col] for r in P]
            acc = 0.0
            for j in range(n):
                t = ts[0] + (ts[-1] - ts[0]) * j / n
                i = min(max(bisect.bisect_left(ts, t), 1), len(ts) - 1)
                acc += vs[i - 1] + (vs[i] - vs[i - 1]) * (t - ts[i - 1]) / (ts[i] - ts[i - 1])
            return acc / n
        av, ao = uavg(1), uavg(2)
        tr = refs()["a_cl_tran"]
        ok = abs(av / tr["vout_avg"] - 1) <= 1e-6 and abs(ao / tr["opout_avg"] - 1) <= 1e-5 and int(m.group(1)) <= 10
        record(label, "pass" if ok else "FAIL", f"{m.group(0)}; avg vout {av:.8f} V (settled tran {tr['vout_avg']:.8f}), avg opout {ao:.8f} V "
               f"(tran {tr['opout_avg']:.8f}); bounds 1e-6 / 1e-5 relative, at most 10 Newton iterations")
    except (CouldNotRun, KeyError, OSError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])


def check_pxf():
    label = "10. .pxf on a_ol: forward Vc column == td .pac; adjoint=1 falls back to the forward engine"
    try:
        d = fixture("a_ol")
        z, _ = td_pac(d, "Vc", 25000.0, "a_ol_pac_pxfref")
        vals, msgs = [], []
        for adj in (0, 1):
            fn = os.path.join(OUT, f"a_ol_pxf{adj}.txt")
            log = run(d + f".options reltol=1e-5\n.pxf 500k 100u vout {PTS} 1 50 1e-6 vout lin 1 25000 25000 method=\"shoot\" solver=\"td\" "
                      f"sources=\"Vc\" sources=\"Vin\" adjoint={adj}\n.control\n{WR}run\nwrdata {fn} pxf_vc pxf_vin\nquit\n.endc\n.end\n", f"a_ol_pxf{adj}")
            r = rd(fn)[0]
            vals.append(complex(r[1], r[2]))
            msgs.append("adjoint engine REFUSED (event-driven XSPICE blocks" in log)
        ok = z is not None and all(abs(v - z) <= 1e-9 * abs(z) for v in vals) and msgs[1] and not msgs[0]
        record(label, "pass" if ok else "FAIL", f".pac {z}; .pxf adjoint=0 {vals[0]}, adjoint=1 {vals[1]} (fallback message: {msgs[1]})")
    except (CouldNotRun, KeyError, OSError, IndexError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])


def check_refusals():
    cases = [
        ("r_tie", "Vc", None, "", "simultaneous input edges of different timing origin at XSPICE instance aand"),
        ("r_mbit", "Vc", None, "", "multi-bit adc_bridge acmp"),
        ("r_vco", "Vc", None, "driven=1", "d_osc aosc: its control input has a small-signal part"),
        ("a_ol", "Vc", lambda d: d.replace("t_rise=5n t_fall=5n", "t_rise=20n t_fall=20n"), "", "still in flight at the end of the period"),
        ("b_ol2", "Vc", lambda d: d.replace("init_phase=90", "init_phase=0"), "", "still in flight at the end of the period"),
    ]
    for deck, src, edit, extra, want in cases:
        label = f"11. refused: {deck}{' (edited)' if edit else ''} -> '{want[:60]}'"
        try:
            d = fixture(deck)
            if edit:
                d = edit(d)
            z, log = td_pac(d, src, 25000.0, f"refuse_{deck}{'_e' if edit else ''}", extra=extra)
            msg = [l for l in refused(log) if "card #714" in l]
            ok = z is None and any(want in l for l in msg)
            record(label, "pass" if ok else "FAIL", (msg[0][:260] if msg else "no card #714 refusal: " + " | ".join(refused(log))[:200]))
        except (CouldNotRun, OSError) as ex:
            record(label, "COULD-NOT-RUN", str(ex)[:300])
    label = "11. refused: .pnoise through digital blocks"
    try:
        log = run(fixture("a_ol") + f".options reltol=1e-5\n.pnoise 500k 100u vout {PTS} 1 50 1e-6 vout Vc lin 1 25000 25000 method=\"shoot\" "
                  f"solver=\"td\"\n.control\n{WR}run\nquit\n.endc\n.end\n", "refuse_pnoise")
        msg = refused(log)
        record(label, "pass" if msg and "PNOISE-TD" in msg[0] else "FAIL", (msg[0][:200] if msg else "not refused"))
    except (CouldNotRun, OSError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])


def check_red():
    d = fixture("a_ol")
    try:
        z, log = td_pac(d, "Vc", 25000.0, "red_refuse", env={"PAC_DIG_REFUSE": "1"})
        ok = z is None and any("event-driven / digital" in l for l in refused(log))
        record("12a. red proof PAC_DIG_REFUSE=1: the s149 refusal comes back", "pass" if ok else "FAIL", (refused(log) or ["not refused"])[0][:200])
    except (CouldNotRun, OSError) as ex:
        record("12a. red proof PAC_DIG_REFUSE", "COULD-NOT-RUN", str(ex)[:300])
    try:
        z, log = td_pac(d, "Vc", 25000.0, "red_noshift", env={"PAC_DIG_NOSHIFT": "1"})
        r, _ = ref_tf("a_ol", 25000.0)
        ok = z is not None and abs(z) <= 1e-6 * abs(r)
        record("12b. red proof PAC_DIG_NOSHIFT=1: no edge shift, the comparator path has zero gain (a_ol 25 kHz)", "pass" if ok else "FAIL",
               f"|H| = {abs(z) if z is not None else 'n/a'} vs FRA {abs(r):.4f}")
        res, log = td_pstb(fixture("a_cl"), [10000.0], "red_noshift_cl", env={"PAC_DIG_NOSHIFT": "1"})
        r, _ = ref_loop("a_cl", 10000.0)
        ok = res is not None and abs(db(res[0][1]) - db(r)) > 20
        record("12c. red proof PAC_DIG_NOSHIFT=1 on the closed loop (a_cl 10 kHz): the loop is open", "pass" if ok else "FAIL",
               f"|T| {db(res[0][1]) if res else 'n/a'} dB vs FRA {db(r):.3f} dB")
    except (CouldNotRun, OSError, KeyError) as ex:
        record("12b/c. red proof PAC_DIG_NOSHIFT", "COULD-NOT-RUN", str(ex)[:300])
    try:
        z, log = td_pac(fixture("b_ol2"), "Vc", 2000.0, "red_pwm_nostate", env={"PAC_DIG_PWM_NOSTATE": "1"})
        r, _ = ref_tf("b_ol2", 2000.0)
        ok = z is not None and abs(db(z) - db(r)) > 1.0
        record("12d. red proof PAC_DIG_PWM_NOSTATE=1: d_pwm transitions without the carried shift (b_ol2 2 kHz)", "pass" if ok else "FAIL",
               f"{db(z) if z is not None else 'n/a'} dB vs FRA {db(r):.4f} dB")
    except (CouldNotRun, OSError, KeyError) as ex:
        record("12d. red proof PAC_DIG_PWM_NOSTATE", "COULD-NOT-RUN", str(ex)[:300])
    try:
        log = run(fixture("a_cl") + f".options reltol=1e-5\n.control\n{WR}pss 500k 1m vout {PTS} 10 150 1e-6 uic method=shoot\nquit\n.endc\n.end\n",
                  "red_noloc_pss", env={"PSS_ADC_LOCATE": "0"})
        ok = "PSS-SHOOT: CONVERGED" not in log
        record("12e. red proof PSS_ADC_LOCATE=0: the closed comparator loop's .pss does not converge (edges on the grid)", "pass" if ok else "FAIL",
               " | ".join(l for l in log.split("\n") if "FAILED" in l or "CONVERGED" in l)[:220])
    except (CouldNotRun, OSError) as ex:
        record("12e. red proof PSS_ADC_LOCATE", "COULD-NOT-RUN", str(ex)[:300])


def check_base():
    if not BASE:
        return
    try:
        d = re.sub(r"^(Vc\s+\S+\s+\S+\s+.*)$", lambda m: m.group(1) + " AC 1", fixture("a_ol"), count=1, flags=re.M)
        log = run(d + f".options reltol=1e-5\n.pac 500k 100u vout {PTS} 1 50 1e-6 lin 1 25000 25000 method=\"shoot\" solver=\"td\"\n"
                  f".control\n{WR}run\nquit\n.endc\n.end\n", "base_a_ol", ng=BASE)
        ok1 = any("event-driven / digital" in l for l in refused(log))
        log2 = run(fixture("a_cl") + f".options reltol=1e-5\n.control\n{WR}pss 500k 1m vout {PTS} 10 150 1e-6 uic method=shoot\nquit\n.endc\n.end\n",
                   "base_a_cl_pss", ng=BASE)
        ok2 = "PSS-SHOOT: CONVERGED" not in log2
        record("13. baseline build: td .pac refused on a_ol, .pss of a_cl not converged", "pass" if ok1 and ok2 else "FAIL",
               f"refused: {ok1}; a_cl .pss not converged: {ok2}")
    except (CouldNotRun, OSError) as ex:
        record("13. baseline", "COULD-NOT-RUN", str(ex)[:300])


def main():
    if not os.path.exists(NG):
        print("ngspice build not found:", NG)
        return 2
    if not os.path.exists(os.path.join(FIX, "refs.json")):
        print("missing", os.path.join(FIX, "refs.json"))
        return 2
    check_tf("1. open loop, comparator PWM (adc_bridge -> d_buffer -> dac_bridge), control-to-output", "a_ol", "Vc", 0.005, 0.02)
    check_tf("2. open loop, latch PWM (d_osc sets a d_dff, comparator resets it; d_inverter pair)", "a_ff", "Vc", 0.005, 0.02)
    if not QUICK:
        check_tf("3. open loop, latch PWM driving IDEAL sw switches", "a_ffsw", "Vc", 0.01, 0.05)
    check_tf("4. open loop, d_pwm with an analog duty (Floquet shift state)", "b_ol2", "Vc", 0.005, 0.02)
    check_tf("5. open loop, clock-only digital path, line-to-output", "c_ol", "Vin", 0.005, 0.02)
    check_loop("6. closed loop, comparator PWM (a_cl): .pstb", "a_cl", (2000.0, 10000.0, 20000.0, 31250.0, 500e3 / 14, 50000.0), 0.01, 0.05)
    check_margin()
    if not QUICK:
        check_loop("7. closed loop, d_osc ramp, comparator with a +-1 mV band (a_cl2): .pstb (2 kHz left out: the FRA spreads 0.07 dB there)", "a_cl2", (10000.0, 20000.0, 50000.0), 0.01, 0.05)
        check_loop("8. closed loop, clock-only digital path, analog comparator (c_cl): .pstb (2 kHz left out: the FRA spreads 0.017 dB there)", "c_cl", (10000.0, 20000.0, 50000.0), 0.01, 0.05)
    check_orbit()
    check_pxf()
    check_refusals()
    check_red()
    check_base()
    npass = sum(1 for r in results if r[1] == "pass")
    print(f"{npass}/{len(results)} pass")
    if any(r[1] == "COULD-NOT-RUN" for r in results):
        return 2
    return 0 if npass == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
