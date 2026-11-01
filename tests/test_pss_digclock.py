#!/usr/bin/env python3
"""`.pss method=shoot` of a converter whose ONLY clock is an XSPICE digital block (d_osc / d_pwm) -- card #712.

    NGSPICE=<prefix>/bin/ngspice [NGSPICE_BASE=<s140 prefix>/bin/ngspice] python3 tests/test_pss_digclock.py [--quick] [--r9]

Exit 0 pass, 1 fail, 2 could-not-run.  Outputs under ~/Work/data/ngspice_augmentation/712/test (never the repo).

Checks (every ngspice runs under systemd-run MemoryMax=3G, num_threads=2):
  1. open-loop buck, d_pwm clock (500 kHz, duty 0.275), `.pss 500k ...` WITHOUT driven=1: the period is taken from Fguess because the d_pwm
     runs at 1 x Fguess; the orbit equals the settled .tran of the same deck (ripple waveform within 1 % of ripple, average within 1e-4).
  2. the same with a d_osc clock.
  3. d_pwm at duty 0.55: the pending clock edge lands on the period boundary (before #712: "Timestep too small" at t = T).
  4. a clock that is not a whole multiple of Fguess is not taken as the period: refused with a message naming the clock.
  5. closed-loop buck with the d_pwm (compensator output as control): the d_pwm decides its next edge time from the analog signal of the
     previous period, a quantity that is not a shooting unknown; the hidden state does not settle and the run is REFUSED (not a wrong orbit).
  6. td `.pac` on the converged circuit runs through the digital blocks since card #714 (a clock-only digital path: its edges are fixed
     waveforms; values checked against FRA in tests/test_digital_pac.py); PAC_DIG_REFUSE=1 restores the s149 refusal (up to s149 it refused).
  7. (--r9) the digitally controlled R9 buck with a d_osc as its sampling clock (period 600 clock cycles).
  8. red proof on the baseline build NGSPICE_BASE (s140): decks 1-3 are refused without driven=1; the two d_pwm decks (1, 3) also do not
     converge with driven=1 (digital state not periodic / "Timestep too small" at t = T); the d_osc deck does converge there.
  9. card #741, latch-PWM buck (fixtures/dosc741: d_osc sets a d_dff, an adc_bridge comparator resets it), d_osc init_phase 306 / 270 / 0:
     a. 306 and 270 (the clock edge 0.05 T / 0.15 T into the period: the d_osc toggles its output in the t = 0 operating point of every
        restarted period run): CONVERGED, no crash, orbit vs pinned settled-.tran numbers (plain .tran 5n 2m with XS_ADC_LOCATE=1, last
        period): cycle-average v(vout) within 0.1 % of the reference ripple (peak-peak), ripple within 1 %.  Up to s151: "Too many
        analog/event-driven solution alternations" in the restarted runs (the limit was compared with the job-total count), then a
        segfault (the op-failure path of DCtran left a freed plot pointer that pss_run_period ended again).
     b. 0 (edge at 0.9 T, the reset is high: the latch never sets, the converter is off; plain .tran: vout 1.2e-7 V): CONVERGED onto the
        all-off orbit, accepted because the pre-stabilization transient itself ends periodic there (up to s151: rejected as the
        degenerate all-off point and reported "FAILED to converge within 50 Newton iterations" with a residual below tol).
     c. red proofs on the same build: PSS_XS_ALTCOUNT_OFF=1 brings the alternation failure back on 306; PSS_ALLOFF_ACCEPT_OFF=1 the
        rejection on 0.
"""
import bisect
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures", "digclock")
FIX741 = os.path.join(HERE, "fixtures", "dosc741")
# card #741: settled plain .tran 5n 2m uic, XS_ADC_LOCATE=1 (comparator edges located exactly, as in the shooting runs), w741 build,
# last period: (cycle-average v(vout), peak-peak v(vout)) in V
REF741 = {306: (4.92017961, 0.00791231), 270: (3.72740854, 0.00716408)}
NG = os.environ.get("NGSPICE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-w712/bin/ngspice")
BASE = os.environ.get("NGSPICE_BASE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s140/bin/ngspice")
OUT = os.path.expanduser("~/Work/data/ngspice_augmentation/712/test")
QUICK = "--quick" in sys.argv
R9 = "--r9" in sys.argv
results = []
T = 2e-6


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"[{state:13s}] {name}: {detail}", flush=True)


def run(deck, name, ng=None, timeout=1800, env=None):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name + ".cir")
    with open(p, "w") as fh:
        fh.write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    try:
        r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=3G", "-p", "MemorySwapMax=0", ng or NG, "-b", p],
                           capture_output=True, text=True, cwd=OUT, env=e, timeout=timeout)
    except (FileNotFoundError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{name}: {ex}")
    return (r.stdout + r.stderr).replace("\r", "\n")


def fixture(name):
    with open(os.path.join(FIX, name)) as fh:
        txt = fh.read()
    txt = re.sub(r"^\.include\s+\S*lgprobe\.lib", ".include " + os.path.join(os.path.dirname(HERE), "stb", "lgprobe.lib"), txt, flags=re.M)
    return re.sub(r"\n\.end\s*$", "\n", txt.rstrip() + "\n")


WR = "set num_threads=2\nset filetype=ascii\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"


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


def pss_deck(deck, tag, f=500e3, extra="", pts=2000):
    pth = os.path.join(OUT, f"{tag}_pss.txt")
    if os.path.exists(pth):
        os.remove(pth)
    return deck + ".options reltol=1e-5\n.control\n" + WR + f"pss {f:.15g} 100u vout {pts} 10 150 1e-6 uic method=shoot {extra}\n" \
        f"wrdata {pth} v(vout)\nquit\n.endc\n.end\n"


def orbit_vs_tran(fx, tag, text_edit=None):
    """returns (ripple error % of ripple, average relative error, log); raises CouldNotRun"""
    deck = fixture(fx)
    if text_edit:
        deck = text_edit(deck)
    log = run(pss_deck(deck, tag), tag + "_pss")
    if "PSS-SHOOT: CONVERGED" not in log:
        return None, None, log
    tp = os.path.join(OUT, f"{tag}_tran.txt")
    if os.path.exists(tp):
        os.remove(tp)
    run(deck + ".options reltol=1e-5\n.control\n" + WR + f"tran 5n 5m 4.99m uic\nwrdata {tp} v(vout)\nquit\n.endc\n.end\n", tag + "_tran")
    P, R = rd(os.path.join(OUT, f"{tag}_pss.txt")), rd(tp)
    rt, rv = [r[0] for r in R], [r[1] for r in R]
    pt, pv = [0.0] + [r[0] for r in P], [P[0][1]] + [r[1] for r in P]
    t0 = rt[-1] - T
    n = 800
    rp = [interp(rt, rv, t0 + T * j / n) for j in range(n + 1)]
    pp = [interp(pt, pv, pt[-1] * j / n) for j in range(n + 1)]
    rip = max(rp) - min(rp)
    err = 100 * max(abs(a - b) for a, b in zip(pp, rp)) / rip
    ap, ar = sum(pp[:-1]) / n, sum(rp[:-1]) / n
    return err, abs(ap - ar) / abs(ar), log


def check_orbit(fx, tag, label, edit=None):
    try:
        err, avg, log = orbit_vs_tran(fx, tag, edit)
    except CouldNotRun as ex:
        record(label, "COULD-NOT-RUN", str(ex))
        return
    if err is None:
        record(label, "FAIL", "did not converge: " + " | ".join(l for l in log.split("\n") if "FAILED" in l or "REFUSED" in l or "Error" in l)[:300])
        return
    m = re.search(r"digital node mismatches (\d+), analog node drift (\S+)", log)
    clk = "XSPICE digital clock" in log and "driven=1" not in log.split("PSS-SHOOT")[0]
    ok = err <= 1.0 and avg <= 1e-4 and clk and m and m.group(1) == "0"
    record(label, "pass" if ok else "FAIL",
           f"orbit vs settled tran: max |d| {err:.4f} % of ripple (bound 1), average rel {avg:.2e} (bound 1e-4), period taken from the digital clock "
           f"without driven=1: {clk}, {m.group(0) if m else 'hidden-state check n/a'}")


def deck741(phase):
    with open(os.path.join(FIX741, "ff_dosc.cir")) as fh:
        txt = fh.read().replace("init_phase=PHASE", f"init_phase={phase}")
    return re.sub(r"\n\.end\s*$", "\n", txt.rstrip() + "\n")


def vout_stats(path):
    """cycle average and peak-peak of v(vout) over the written PSS period (column 1 of wrdata)"""
    P = rd(path)
    t, v = [r[0] for r in P], [r[1] for r in P]
    n = 2000
    s = [interp(t, v, t[-1] * j / n) for j in range(n + 1)]
    return sum(s[:-1]) / n, max(s) - min(s)


def check_741():
    finished = lambda log: re.search(r"ngspice\S* done", log) is not None
    for ph in (306, 270):
        label = f"9a. #741 latch-PWM buck, d_osc init_phase={ph}: CONVERGED, no crash, orbit vs settled .tran (average within 0.1 % of ripple)"
        try:
            tag = f"d741_{ph}"
            log = run(pss_deck(deck741(ph), tag, pts=4000), tag + "_pss")
            if "PSS-SHOOT: CONVERGED" not in log:
                bad = [l for l in log.split("\n") if "FAILED" in l or "Too many" in l or "Error" in l]
                record(label, "FAIL", ("no crash" if finished(log) else "CRASHED (no 'done' line)") + "; not converged: " + " | ".join(bad)[:240])
                continue
            avg, pp = vout_stats(os.path.join(OUT, f"{tag}_pss.txt"))
            ra, rp = REF741[ph]
            da, dr = abs(avg - ra) / rp, abs(pp - rp) / rp
            ok = finished(log) and da <= 1e-3 and dr <= 1e-2
            record(label, "pass" if ok else "FAIL",
                   f"average {avg:.8f} V vs {ra:.8f} V: {100 * da:.4f} % of ripple (bound 0.1); ripple {pp * 1e3:.4f} mV vs {rp * 1e3:.4f} mV: "
                   f"{100 * dr:.3f} % (bound 1); finished: {finished(log)}")
        except CouldNotRun as ex:
            record(label, "COULD-NOT-RUN", str(ex))
    label = "9b. #741 init_phase=0 (latch never set, converter off): CONVERGED onto the all-off orbit, |average v(vout)| < 1 uV, accepted as such"
    try:
        log = run(pss_deck(deck741(0), "d741_0", pts=4000), "d741_0_pss")
        conv = "PSS-SHOOT: CONVERGED" in log
        msg = "does not switch in steady state" in log
        avg = vout_stats(os.path.join(OUT, "d741_0_pss.txt"))[0] if conv else None
        ok = conv and msg and finished(log) and abs(avg) < 1e-6
        record(label, "pass" if ok else "FAIL",
               f"converged: {conv}, all-off acceptance message: {msg}, average v(vout) {avg if avg is None else f'{avg:.3e} V'}"
               + ("" if conv else " | " + " | ".join(l for l in log.split("\n") if "FAILED" in l)[:200]))
    except CouldNotRun as ex:
        record(label, "COULD-NOT-RUN", str(ex))
    label = "9c. #741 red proofs on this build: PSS_XS_ALTCOUNT_OFF=1 -> alternation failure on 306; PSS_ALLOFF_ACCEPT_OFF=1 -> all-off rejected on 0"
    try:
        la = run(pss_deck(deck741(306), "d741_red_alt", pts=4000), "d741_red_alt_pss", env={"PSS_XS_ALTCOUNT_OFF": "1"})
        lo = run(pss_deck(deck741(0), "d741_red_off", pts=4000), "d741_red_off_pss", env={"PSS_ALLOFF_ACCEPT_OFF": "1"})
        oka = "Too many analog/event-driven solution alternations" in la and "PSS-SHOOT: CONVERGED" not in la and finished(la)
        oko = "degenerate all-off fixed point" in lo and "PSS-SHOOT: CONVERGED" not in lo
        record(label, "pass" if oka and oko else "FAIL", f"alternation failure back (and no crash): {oka}; all-off rejected: {oko}")
    except CouldNotRun as ex:
        record(label, "COULD-NOT-RUN", str(ex))


def main():
    for p in (NG, BASE):
        if not os.path.exists(p):
            print("ngspice build not found:", p)
            return 2
    check_orbit("buck_ol_dpwm.cir", "dp", "1. open-loop buck, d_pwm clock: PSS orbit (period from Fguess) vs settled tran")
    check_orbit("buck_ol_dosc.cir", "do", "2. open-loop buck, d_osc clock: PSS orbit vs settled tran")
    if not QUICK:
        check_orbit("buck_ol_dpwm.cir", "dp55", "3. d_pwm duty 0.55 (clock edge on the period boundary): PSS orbit vs settled tran",
                    edit=lambda d: d.replace("Vc c 0 DC 0.275", "Vc c 0 DC 0.55"))
    # 4. non-multiple
    try:
        log = run(pss_deck(fixture("buck_ol_dpwm.cir"), "nm", f=300e3), "nm_pss")
        ok = "none runs at a whole multiple of Fguess" in log and "needs a driven" in log and "CONVERGED" not in log
        record("4. d_pwm at 500 kHz is not the clock of a 300 kHz analysis: refused naming the clock", "pass" if ok else "FAIL",
               " | ".join(l for l in log.split("\n") if "whole multiple" in l or "needs a driven" in l)[:260])
    except CouldNotRun as ex:
        record("4. non-multiple clock", "COULD-NOT-RUN", str(ex))
    # 5. closed loop
    try:
        log = run(pss_deck(fixture("buck_cl_dpwm.cir"), "cl"), "cl_pss")
        ok = "CONVERGED" not in log and "does not settle" in log
        record("5. closed-loop buck, d_pwm (edge time decided from the analog signal): refused, not a wrong orbit", "pass" if ok else "FAIL",
               " | ".join(l for l in log.split("\n") if "does not settle" in l or "CONVERGED" in l)[:260] or "no message")
    except CouldNotRun as ex:
        record("5. closed-loop refusal", "COULD-NOT-RUN", str(ex))
    # 6. pac
    try:
        d = fixture("buck_ol_dpwm.cir").replace("Vin vin 0 DC {vin_v}", "Vin vin 0 DC {vin_v} AC 1")
        d += ".options reltol=1e-5\n.pac 500k 0 vout 2000 1 8 1e-6 dec 3 1k 10k method=\"shoot\" solver=\"td\"\n.control\nset num_threads=2\nrun\nquit\n.endc\n.end\n"
        log = run(d, "pac")
        ok = "PSS-SHOOT: CONVERGED" in log and "PAC-TD: done" in log and "REFUSED" not in log and "event-driven XSPICE blocks (card #714)" in log
        logr = run(d, "pac_refuse", env={"PAC_DIG_REFUSE": "1"})
        okr = "PAC-TD: REFUSED" in logr and "event-driven / digital" in logr
        record("6. td .pac through the clock-only digital blocks runs (card #714); PAC_DIG_REFUSE=1 restores the refusal", "pass" if ok and okr else "FAIL",
               f"runs: {ok}; PAC_DIG_REFUSE=1 refused: {okr}")
    except CouldNotRun as ex:
        record("6. td .pac refusal", "COULD-NOT-RUN", str(ex))
    # 7. R9
    if R9:
        try:
            fxd = fixture("r9_dosc.cir").replace("../../stb/", os.path.join(os.path.dirname(HERE), "stb") + "/")
            pth = os.path.join(OUT, "r9d_pss.txt")
            log = run(fxd + ".options reltol=1e-4\n.control\n" + WR + f"pss {1 / 1.2e-3:.15g} 12e-3 vout 6000 10 150 1e-6 uic method=shoot\nwrdata {pth} v(vout) v(code)\nquit\n.endc\n.end\n",
                      "r9d_pss", timeout=3000)
            tp = os.path.join(OUT, "r9d_tran.txt")
            run(fxd + ".options reltol=1e-4\n.control\n" + WR + f"tran 100n 60m uic\nwrdata {tp} v(vout) v(code)\nquit\n.endc\n.end\n", "r9d_tran", timeout=3000)
            if "PSS-SHOOT: CONVERGED" not in log:
                record("7. R9 digital buck with d_osc sampling clock", "FAIL", " | ".join(l for l in log.split("\n") if "FAILED" in l or "REFUSED" in l)[:260])
            else:
                P, R = rd(pth), rd(tp)
                Tr = 1.2e-3
                rt, rv = [r[0] for r in R], [r[1] for r in R]
                pt, pv = [0.0] + [r[0] for r in P], [P[0][1]] + [r[1] for r in P]
                t0 = rt[-1] - Tr
                n = 3000
                c = [interp(rt, rv, t0 + Tr * j / n) for j in range(n + 1)]
                p = [interp(pt, pv, pt[-1] * j / n) for j in range(n + 1)]
                rip = max(c) - min(c)
                d = max(abs(a - b) for a, b in zip(c, p))
                am = abs(sum(p[:-1]) - sum(c[:-1])) / n / abs(sum(c[:-1]) / n)
                record("7. R9 digital buck, d_osc sampling clock, period 600 clock cycles: PSS orbit vs settled tran",
                       "pass" if d / rip <= 0.01 and am <= 1e-4 else "FAIL", f"ripple {rip * 1e3:.3f} mV, max |d| {100 * d / rip:.3f} % of ripple (bound 1), average rel {am:.2e} (bound 1e-4)")
            lb = run(fxd + ".options reltol=1e-4\n.control\n" + WR + f"pss {1 / 1.2e-3:.15g} 12e-3 vout 6000 10 150 1e-6 uic method=shoot\nquit\n.endc\n.end\n",
                     "r9d_base_pss", ng=BASE)
            record("7b. red proof, baseline build: R9 with the d_osc clock is refused without driven=1", "pass" if "needs a driven" in lb and "CONVERGED" not in lb else "FAIL",
                   "refused" if "needs a driven" in lb else "NOT refused")
        except CouldNotRun as ex:
            record("7. R9 with d_osc", "COULD-NOT-RUN", str(ex))
    # 8. red proof on the baseline (a d_osc deck converges on the baseline with driven=1: only the missing period detection is red there)
    for fx, tag, edit, must_fail in (("buck_ol_dpwm.cir", "dp", None, True), ("buck_ol_dosc.cir", "do", None, False),
                                     ("buck_ol_dpwm.cir", "dp55", lambda d: d.replace("Vc c 0 DC 0.275", "Vc c 0 DC 0.55"), True)):
        try:
            deck = fixture(fx)
            if edit:
                deck = edit(deck)
            l0 = run(pss_deck(deck, tag + "_base"), tag + "_base_pss", ng=BASE)
            l1 = run(pss_deck(deck, tag + "_base1", extra="driven=1"), tag + "_base1_pss", ng=BASE)
            conv1 = "PSS-SHOOT: CONVERGED" in l1
            ok = "needs a driven" in l0 and "CONVERGED" not in l0 and (not conv1 or not must_fail)
            record(f"8. red proof, baseline build, {fx}{' (duty 0.55)' if edit else ''}: refused without driven=1" +
                   (", no converged orbit with driven=1" if must_fail else ""), "pass" if ok else "FAIL",
                   f"no driven=1: {'refused' if 'needs a driven' in l0 else 'NOT refused'}; driven=1: "
                   f"{'converged' if conv1 else ' | '.join(x[:100] for x in l1.split(chr(10)) if 'FAILED' in x or 'doAnalyses' in x)[:200]}")
        except CouldNotRun as ex:
            record(f"8. red proof {fx}", "COULD-NOT-RUN", str(ex))
    # 9. card #741
    check_741()
    npass = sum(1 for r in results if r[1] == "pass")
    print(f"{npass}/{len(results)} pass")
    if any(r[1] == "COULD-NOT-RUN" for r in results):
        return 2
    return 0 if npass == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
