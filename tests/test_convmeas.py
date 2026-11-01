#!/usr/bin/env python3
"""Validation for tools/convmeas.py (card #497, feature-survey #16).

Three layers:
  1. Pure-Python synthetic checks (hand-built (t,y) arrays with exactly
     known answers) for every function -- ripple_pp, average_value,
     average_power/efficiency, switching_frequency, duty_cycle_threshold,
     settling_time. No ngspice needed, always run.
  2. A closed-form underdamped 2nd-order step response (textbook formula,
     built as a dense Python array -- still no ngspice) for
     overshoot_undershoot(), since that one genuinely needs a non-trivial
     waveform shape to exercise meaningfully.
  3. ngspice-driven fixtures: tests/fixtures/convmeas_sine.cir (ripple_pp
     against a hand-computed 2*VA) and tests/fixtures/convmeas_rlc_step.cir
     (the SAME closed-form 2nd-order response, but now the real simulated
     trace of an actual R-L-C circuit built for exactly that zeta -- checks
     overshoot/undershoot/settling_time against the closed-form numbers on
     REAL ngspice output, not just a hand-built array); plus
     refs/r3_vm_buck.cir (a real, independently-documented switching
     converter, refs/README.md's own MEASURED table) for ripple/efficiency/
     switching-frequency/duty-cycle on a genuine converter trace.

Needs system ngspice (`ngspice -b`) for layer 3 -- nothing here needs the
private fork. If ngspice is not on PATH, layer-3 checks report
COULD-NOT-RUN (exit 2), never a silent pass; layers 1/2 still run.

    python3 tests/test_convmeas.py        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXTURES = os.path.join(HERE, "fixtures")
REFS = os.path.join(ROOT, "refs")
sys.path.insert(0, os.path.join(ROOT, "tools"))
import convmeas  # noqa: E402
import stb  # noqa: E402

DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/tests_convmeas")
NGSPICE = os.environ.get("NGSPICE", "ngspice")

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))


def close(a, b, tol):
    return abs(a - b) <= tol * max(1.0, abs(b))


def _have_ngspice():
    return shutil.which(NGSPICE) is not None or os.path.exists(NGSPICE)


# ---------------------------------------------------------------------------
# layer 1: pure-Python synthetic checks
# ---------------------------------------------------------------------------

def check_ripple_and_average():
    n = 2001
    t = [k * 1e-6 for k in range(n)]  # 0 .. 2ms, 1us steps
    y = [2.5 + 0.3 * math.sin(2 * math.pi * 1000 * ti) for ti in t]  # 1kHz sine on 2.5V offset
    rp = convmeas.ripple_pp(t, y, 0.5e-3, 1.5e-3)  # exactly one period window
    record("ripple_pp on a hand-built 1kHz/0.3V sine = 0.6V", close(rp, 0.6, 1e-3), f"got {rp:.6f}")

    # average_value of a symmetric sine over an integer number of periods = the DC offset
    av = convmeas.average_value(t, y, 0.0, 1.999e-3)  # ~2 periods
    record("average_value of a sine over ~2 periods = its DC offset (2.5V)", close(av, 2.5, 2e-3),
           f"got {av:.6f}")

    # a plain 0/1 square wave, hand-built at fine resolution: period 10us, duty 0.3
    n2 = 100000
    t2 = [k * 0.1e-6 for k in range(n2)]  # 0.1us steps, 10ms span
    y2 = [1.0 if (ti % 10e-6) < 3e-6 else 0.0 for ti in t2]
    av2 = convmeas.average_value(t2, y2, 0.0, 1e-3)  # exactly 100 periods
    record("average_value of a hand-built 30% duty square wave (period=10us) = 0.30", close(av2, 0.30, 1e-3),
           f"got {av2:.6f}")


def check_power_efficiency():
    t = [0.0, 1.0, 2.0, 3.0]
    vin = [10.0] * 4
    iin = [-2.0] * 4  # native i(Vin) sign: current INTO the + terminal is negative when the source is DELIVERING 2A
    vout = [5.0] * 4
    iout = [3.0] * 4
    pin = -convmeas.average_power(t, vin, iin, 0.0, 3.0)
    pout = convmeas.average_power(t, vout, iout, 0.0, 3.0)
    record("average_power: constant V*I = 10*2=20W in, 5*3=15W out", close(pin, 20.0, 1e-9) and close(pout, 15.0, 1e-9),
           f"pin={pin} pout={pout}")
    eff = convmeas.efficiency(t, vin, iin, vout, iout, 0.0, 3.0)
    record("efficiency = 15/20 = 75% on constant DC power", close(eff, 0.75, 1e-9), f"got {eff}")


def check_switching_and_duty():
    # hand-built ideal 25kHz, 40% duty square wave, 20ns resolution, 12 periods
    period = 1 / 25e3
    duty = 0.40
    dt = 20e-9
    n = int(12 * period / dt)
    t = [k * dt for k in range(n)]
    y = [1.0 if (ti % period) < duty * period else 0.0 for ti in t]
    t0, t1 = 2 * period, t[-1] - period  # a clean interior window, several periods
    fsw = convmeas.switching_frequency(t, y, t0, t1, threshold=0.5)
    record("switching_frequency on a hand-built ideal 25kHz square wave", close(fsw, 25e3, 1e-3), f"got {fsw:.3f}")
    dc = convmeas.duty_cycle_threshold(t, y, t0, t0 + period, threshold=0.5)
    record("duty_cycle_threshold on a hand-built 40% duty square wave", close(dc, 0.40, 2e-2), f"got {dc:.4f}")
    dc_avg = convmeas.average_value(t, y, t0, t0 + period)
    record("average_value duty (AVG-over-one-period trick) also = 0.40", close(dc_avg, 0.40, 2e-2), f"got {dc_avg:.4f}")

    # red-proof: fewer than 2 crossings -> None, not a fabricated number
    flat = [0.0] * len(t)
    fsw_none = convmeas.switching_frequency(t, flat, t0, t1, threshold=0.5)
    record("red-proof: switching_frequency on a flat (never-crossing) signal returns None", fsw_none is None,
           f"got {fsw_none}")


def check_settling_time():
    # y = 2*band, then drops to target and stays -- clean, discrete, hand-picked
    t = list(range(11))
    target, band = 0.0, 0.5
    y = [5.0] * 5 + [0.0] * 6  # samples t=0..4 outside band, t=5..10 inside
    st = convmeas.settling_time(t, y, target, band)
    record("settling_time: last-outside-band sample (t=4) on a clean step", st == 4, f"got {st}")

    st0 = convmeas.settling_time(t, [0.0] * 11, target, band, t_start=3.0)
    record("settling_time: already within band for all of [t_start,end] returns t_start", st0 == 3.0, f"got {st0}")

    st_never = convmeas.settling_time(t, [5.0] * 11, target, band)
    record("red-proof: settling_time returns None (never a fabricated time) when never settled", st_never is None,
           f"got {st_never}")


# ---------------------------------------------------------------------------
# layer 2: closed-form 2nd-order step response, pure Python
# ---------------------------------------------------------------------------

ZETA, LVAL, CVAL = 0.3, 1e-3, 1e-6
RVAL = 2 * ZETA * math.sqrt(LVAL / CVAL)   # must match tests/fixtures/convmeas_rlc_step.cir's .param rval
WN = 1 / math.sqrt(LVAL * CVAL)
WD = WN * math.sqrt(1 - ZETA ** 2)
MP_FRAC = math.exp(-ZETA * math.pi / math.sqrt(1 - ZETA ** 2))


def closed_form(t, td=10e-6):
    if t < td:
        return 0.0
    tt = t - td
    phi = math.acos(ZETA)
    return 1.0 - math.exp(-ZETA * WN * tt) / math.sqrt(1 - ZETA ** 2) * math.sin(WD * tt + phi)


def check_closed_form_overshoot():
    t = [k * 2e-7 for k in range(6000)]  # 0..1.2ms, 0.2us steps -- resolves the ~104us peak well
    y = [closed_form(ti) for ti in t]
    res = convmeas.overshoot_undershoot(t, y, 10e-6, 0.0, 1.0)
    record("overshoot_undershoot (closed-form array) matches Mp=37.2326% for zeta=0.3",
           close(res["overshoot_pct"], MP_FRAC * 100, 5e-3), f"got {res['overshoot_pct']:.4f}%")
    record("overshoot_undershoot (closed-form array) matches Mp^2=13.8627% undershoot",
           close(res["undershoot_pct"], MP_FRAC ** 2 * 100, 1e-2), f"got {res['undershoot_pct']:.4f}%")

    # red-proof: a deliberately wrong zeta's Mp must NOT match
    wrong_mp = math.exp(-0.6 * math.pi / math.sqrt(1 - 0.6 ** 2)) * 100  # zeta=0.6's Mp, not 0.3's
    record("red-proof: zeta=0.6's Mp does not match this zeta=0.3 circuit's measured overshoot",
           not close(res["overshoot_pct"], wrong_mp, 5e-3), f"measured {res['overshoot_pct']:.3f}% vs wrong-zeta {wrong_mp:.3f}%")


# ---------------------------------------------------------------------------
# layer 3: ngspice-driven fixtures
# ---------------------------------------------------------------------------

def check_sine_fixture():
    deck = os.path.join(FIXTURES, "convmeas_sine.cir")
    outdir = os.path.join(DATA_ROOT, "sine")
    trace = convmeas.run_tran(deck, "Vtest", "sig", 1e9, None, "20n", "500u", outdir, NGSPICE)
    t = trace["t"]
    rp = convmeas.ripple_pp(t, trace["vout"], 100e-6, 400e-6)
    record("ripple_pp on the real ngspice SIN(2.5 0.6 10k) trace = 1.2V (hand-computed 2*VA)",
           close(rp, 1.2, 5e-3), f"got {rp:.6f}")


def check_rlc_step_fixture():
    deck = os.path.join(FIXTURES, "convmeas_rlc_step.cir")
    outdir = os.path.join(DATA_ROOT, "rlc_step")
    trace = convmeas.run_tran(deck, "Vin", "vout", 1e9, None, "20n", "1.01m", outdir, NGSPICE)
    t = trace["t"]
    vout = trace["vout"]
    y_before = convmeas.average_value(t, vout, 0.0, 9e-6)
    y_after = convmeas.average_value(t, vout, 9e-4, 1.01e-3)
    res = convmeas.overshoot_undershoot(t, vout, 10e-6, y_before, y_after)
    record("real ngspice R-L-C step trace: overshoot matches closed-form Mp=37.2326% (zeta=0.3)",
           close(res["overshoot_pct"], MP_FRAC * 100, 1e-2), f"got {res['overshoot_pct']:.4f}%")
    record("real ngspice R-L-C step trace: undershoot matches closed-form Mp^2=13.8627%",
           close(res["undershoot_pct"], MP_FRAC ** 2 * 100, 2e-2), f"got {res['undershoot_pct']:.4f}%")

    # settling_time to a +-2% band, cross-checked against a dense closed-form scan
    band = 0.02
    st = convmeas.settling_time(t, vout, 1.0, band, t_start=10e-6)
    fine = [k * 1e-7 for k in range(int(10e-3 / 1e-7))]
    fine = [ti for ti in fine if ti <= t[-1]]
    last_outside_cf = max((ti for ti in fine if abs(closed_form(ti) - 1.0) > band), default=None)
    record("settling_time (+-2%) on the real trace matches a dense closed-form scan within 5us",
           st is not None and last_outside_cf is not None and abs(st - last_outside_cf) < 5e-6,
           f"trace: {st}  closed-form scan: {last_outside_cf}")


def check_r3_buck_fixture():
    """refs/r3_vm_buck.cir -- a real, independently-documented (refs/README.md)
    switching converter, not built for this card. Compared against its own
    already-published MEASURED numbers, with a documented, honest tolerance:
    ripple/duty use a loose (~15%) bound (this tool re-runs its OWN .tran at
    its OWN tstep, not the deck's original 10n/.meas windows, so some
    numeric difference from a different-but-equally-valid sampling is
    expected, not a bug); switching frequency is compared against the IDEAL
    500kHz clock (not the deck's own single-period 500.30kHz .meas number)
    because switching_frequency() here averages over several periods
    end-to-end (first crossing to last), which converges to the clock's true
    long-run rate -- a materially different, and for this purpose more
    precise, method than a single TRIG/TARG period measurement; see
    docs/convmeas.md for the full discussion."""
    deck = os.path.join(REFS, "r3_vm_buck.cir")
    outdir = os.path.join(DATA_ROOT, "r3_vm_buck")
    trace = convmeas.run_tran(deck, "Vin", "vout", 3.3, "pwm", "10n", "1.5e-3", outdir, NGSPICE, timeout=120)
    res = convmeas.compute(trace, period=2e-6, periods=5, tstop=1.5e-3, threshold=0.5)

    record("r3_vm_buck.cir: ripple_pp within 15% of the deck's own documented 4.836mV",
           close(res["ripple_pp_v"], 4.836e-3, 0.15), f"got {res['ripple_pp_v']*1e3:.4f} mV")
    record("r3_vm_buck.cir: duty (avg) within 10% of the deck's own documented 0.2777",
           close(res["duty_cycle_avg"], 0.2777, 0.10), f"got {res['duty_cycle_avg']:.4f}")
    record("r3_vm_buck.cir: switching_frequency (multi-period average) within 0.5% of the ideal 500kHz clock",
           close(res["switching_frequency_hz"], 500e3, 5e-3), f"got {res['switching_frequency_hz']:.1f} Hz")
    record("r3_vm_buck.cir: efficiency is a sane 0-100%, non-nan (no independently-documented reference number)",
           res["efficiency"] is not None and 0.0 < res["efficiency"] < 1.0,
           f"got {100*res['efficiency']:.2f}% (Pin={res['p_in_w']:.4g}W Pout={res['p_out_w']:.4g}W)")


def check_cli_could_not_run():
    saved = convmeas.NGSPICE if hasattr(convmeas, "NGSPICE") else None
    try:
        convmeas.run_tran(os.path.join(FIXTURES, "convmeas_sine.cir"), "Vtest", "sig", 1e9, None,
                           "20n", "500u", os.path.join(DATA_ROOT, "bogus"), "/nonexistent/ngspice")
        record("red-proof: missing ngspice binary gives CouldNotRun", False, "did not raise")
    except stb.CouldNotRun as e:
        record("red-proof: missing ngspice binary gives CouldNotRun", True, str(e)[:60])


def check_cli_subprocess():
    """One real end-to-end CLI invocation (argparse, main(), summary.json) --
    the library-level checks above cover the maths; this covers the CLI
    plumbing itself actually working."""
    outdir = os.path.join(DATA_ROOT, "cli_sine")
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    cmd = [sys.executable, os.path.join(ROOT, "tools", "convmeas.py"),
           os.path.join(FIXTURES, "convmeas_sine.cir"), "--vin", "Vtest", "--vout-node", "sig",
           "--rload-ohms", "1e9", "--period", "100e-6", "--periods", "3", "--tstop", "500e-6",
           "--out", outdir, "--ngspice", NGSPICE]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    ok = proc.returncode == 0 and os.path.exists(os.path.join(outdir, "summary.json"))
    record("CLI subprocess run (convmeas_sine.cir) exits 0 and writes summary.json", ok,
           f"exit={proc.returncode} stderr_tail={proc.stderr[-200:]}")


def main():
    check_ripple_and_average()
    check_power_efficiency()
    check_switching_and_duty()
    check_settling_time()
    check_closed_form_overshoot()

    if _have_ngspice():
        try:
            check_sine_fixture()
            check_rlc_step_fixture()
            check_r3_buck_fixture()
            check_cli_could_not_run()
            check_cli_subprocess()
        except stb.CouldNotRun as e:
            record("run", None, str(e))
    else:
        record("ngspice on PATH", None, f"{NGSPICE!r} not found -- layer-3 (ngspice-dependent) checks "
                "skipped as COULD-NOT-RUN, not silently passed")

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    npass = sum(s == "pass" for _, s, _ in results)
    print(f"\n{npass}/{len(results)} pass")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
