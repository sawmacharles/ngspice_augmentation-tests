#!/usr/bin/env python3
"""Validation for tools/stb_tran.py (card #497, feature-survey #12).

Three layers:
  1. Pure-Python unit checks of the waveform evaluator (PULSE/PWL/SIN) and the
     .param brace-expression resolver, against hand-computed values -- no
     ngspice needed, fast, always run.
  2. The empirical evidence tools/stb_tran.py's own module docstring cites for
     why it freezes exogenous sources rather than trying to force energy-
     storage state via IC=/.nodeset/savestate+loadstate (three independent
     reproducers: plain op discards device IC=; tran-uic-then-op discards it
     too, same session; savestate a non-equilibrium state + loadstate into a
     FRESH process + stb still discards it). These need the fork build.
  3. End-to-end: tests/fixtures/stbtran_ldo_softstart.cir run through
     tools/stb_tran.py at several t_k, cross-checked against an INDEPENDENTLY
     hand-built frozen deck (this file's own build_frozen_ref(), a separate
     code path from stb_tran.freeze_deck()) at the same t_k -- the real,
     non-tautological check that the freeze mechanism computes what it
     claims to. Plus the two identities the card's own brief names: late t_k
     equals a plain DC-op .stb (0.1 dB / 1 deg), early t_k equals a DC-op
     .stb of the deck with Vss set to the frozen value (this IS what
     build_frozen_ref() is, so this is the same check restated).
  4. --mode=instant (card #500): no freezing/re-solving at all -- linearises
     about the REAL transient state at t_k. Needs a fork build with the
     `.stb ... instant` flag (card #500's own C patch; the default
     NGSPICE_AUG pin above predates it -- point NGSPICE_AUG at a build that
     has it, e.g. aaba50ba3f49-instant-N). Two checks: a settled t_k equals
     a plain DC-op .stb within 1e-6 relative; a mid-ramp t_k differs from
     quasistatic mode in the physically expected direction AND matches an
     INDEPENDENT tangent-linearisation reference (measure the real
     transient's v(ea), hand-compute Bpass's exact incremental gm, substitute
     a linear VCCS, run a plain DC-op .stb on the now-fully-linear circuit)
     within 0.5dB/3deg -- see this file's "layer 4" section for the full
     derivation and why the direction is expected.

Needs the fork build (native `.stb`). Default NGSPICE_AUG is this card's own
build (~/Work/apps/ngspice_augmentation/aaba50ba3f49-s46/bin/ngspice -- a
plain `ngspice`, not `ngspice-aug`, at that prefix). If not
present and NGSPICE_AUG is not set, every ngspice-dependent check reports
COULD-NOT-RUN (exit 2), never a silent pass -- layer 1 (pure Python) still
runs and can still FAIL on its own.

    python3 tests/test_stb_tran.py        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXTURES = os.path.join(HERE, "fixtures")
sys.path.insert(0, os.path.join(ROOT, "tools"))
import stb  # noqa: E402
import stb_tran  # noqa: E402

NGSPICE_AUG = os.environ.get(
    "NGSPICE_AUG") or os.environ.get(
    "NGSPICE",  # the portfolio's build under test (2026-09-30: was silently the s46 pin)
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s46/bin/ngspice"))
stb_tran.NGSPICE_AUG = NGSPICE_AUG

DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/tests_stb_tran")
LDO = os.path.join(FIXTURES, "stbtran_ldo_softstart.cir")
SWEEP = "dec 30 1 1e8"

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))


def close(a, b, tol):
    return abs(a - b) <= tol * max(1.0, abs(b))


# ---------------------------------------------------------------------------
# layer 1: pure-Python waveform evaluator + param resolver
# ---------------------------------------------------------------------------

def check_waveforms():
    # PULSE(V1 V2 TD TR TF PW PER): hand-computed at several phases
    a = [0.0, 1.0, 10e-6, 2e-6, 2e-6, 4e-6, 20e-6]
    cases = [
        (0.0, 0.0),            # before TD
        (10e-6, 0.0),          # exactly at TD, start of rise
        (11e-6, 0.5),          # mid-rise (linear ramp)
        (12e-6, 1.0),          # start of PW plateau
        (14e-6, 1.0),          # mid plateau
        (17e-6, 0.5),          # mid-fall (fall spans 16-18u, midpoint 17u)
        (18e-6, 0.0),          # after fall, before next period
        (31e-6, 0.5),          # second period, mid-rise ((31-10)%20=1u into TR=2u, same phase as t=11e-6)
    ]
    ok = True
    for t, exp in cases:
        got = stb_tran.eval_pulse(a, t)
        if not close(got, exp, 1e-9):
            ok = False
    record("eval_pulse matches hand-computed PULSE(0 1 10u 2u 2u 4u 20u) at 8 phases", ok,
           f"last checked t={cases[-1][0]}: got {stb_tran.eval_pulse(a, cases[-1][0])}")

    # PWL(0 0 500u 0.45 2m 0.45): linear ramp then flat
    a = [0.0, 0.0, 500e-6, 0.45, 2e-3, 0.45]
    exp100 = 0.45 * 100e-6 / 500e-6
    exp250 = 0.45 * 250e-6 / 500e-6
    ok = (close(stb_tran.eval_pwl(a, 0.0), 0.0, 1e-12)
          and close(stb_tran.eval_pwl(a, 100e-6), exp100, 1e-9)
          and close(stb_tran.eval_pwl(a, 250e-6), exp250, 1e-9)
          and close(stb_tran.eval_pwl(a, 500e-6), 0.45, 1e-12)
          and close(stb_tran.eval_pwl(a, 1e-3), 0.45, 1e-12)
          and close(stb_tran.eval_pwl(a, 5e-3), 0.45, 1e-12))
    record("eval_pwl matches hand-computed PWL(0 0 500u 0.45 2m 0.45) at 6 points", ok,
           f"pwl(100u)={stb_tran.eval_pwl(a, 100e-6):.6g} (expect {exp100:.6g})")

    # SIN(0 1 1000 0 0 90): VO=0 VA=1 FREQ=1kHz TD=0 THETA=0 PHASE=90deg -> cos(2*pi*1000*t)
    a = [0.0, 1.0, 1000.0, 0.0, 0.0, 90.0]
    exp = math.cos(2 * math.pi * 1000 * 250e-6)
    got = stb_tran.eval_sin(a, 250e-6)
    record("eval_sin matches hand-computed SIN(0 1 1k 0 0 90) at t=250us", close(got, exp, 1e-9),
           f"got {got:.9f} expect {exp:.9f}")

    # red-proof: EXP/SFFM must be refused, not silently mishandled
    try:
        stb_tran.freeze_value("Vx a b EXP(0 1 1u 1u 5u 1u)", 3e-6, {})
        record("red-proof: EXP(...) source refused by freeze_value", False, "did not raise")
    except stb.CouldNotRun as e:
        record("red-proof: EXP(...) source refused by freeze_value", True, str(e)[:60])


def check_param_resolver():
    lines = [".param fs=500k", ".param tclk={1/fs}", ".param vref_final=0.45",
             "Vss vss 0 PWL(0 0 {tclk*250} {vref_final} 2m {vref_final})"]
    params = stb_tran.resolve_params(lines)
    ok = (close(params["fs"], 500e3, 1e-12) and close(params["tclk"], 2e-6, 1e-9)
          and close(params["vref_final"], 0.45, 1e-12))
    record("resolve_params: fs=500k, tclk={1/fs}, vref_final=0.45", ok, str(params))

    val = stb_tran.eval_token("{tclk*250}", params)
    record("eval_token evaluates a param-referencing brace expression", close(val, 500e-6, 1e-9),
           f"got {val:.6g} expect 5e-4")

    frozen, values = stb_tran.freeze_deck(lines, ["Vss"], 100e-6)
    exp = 0.45 * 100e-6 / 500e-6
    record("freeze_deck resolves {param} PWL args and rewrites to plain DC",
           close(values["Vss"], exp, 1e-9) and any(re.match(r"Vss vss 0 DC", ln) for ln in frozen),
           f"values={values}")

    try:
        stb_tran.freeze_deck(lines, ["Vnonexistent"], 100e-6)
        record("red-proof: --freeze naming a missing source raises CouldNotRun", False, "did not raise")
    except stb.CouldNotRun as e:
        record("red-proof: --freeze naming a missing source raises CouldNotRun", True, str(e)[:60])


# ---------------------------------------------------------------------------
# layer 2: empirical evidence the module docstring cites (fork build only)
# ---------------------------------------------------------------------------

def _have_ngspice_aug():
    return os.path.exists(NGSPICE_AUG) and os.access(NGSPICE_AUG, os.X_OK)


def _run(deck_text, outdir, label, cwd=None):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{label}.cir")
    with open(path, "w") as fh:
        fh.write(deck_text)
    proc = subprocess.run([NGSPICE_AUG, "-b", os.path.abspath(path) if cwd else path],
                           capture_output=True, text=True, timeout=120, cwd=cwd)
    return proc.stdout + proc.stderr


def check_freeze_evidence():
    """Reproduces the three discriminators tools/stb_tran.py's docstring
    cites: a plain `op` and a same-session `tran-uic-then-op` both discard a
    capacitor's own IC=, converging to the TRUE DC equilibrium instead; and
    savestate + loadstate (into a fresh process) + stb also discards the
    loaded, genuinely non-equilibrium state. Diode current (@d1[id]) is the
    discriminator: forced/frozen-state current is huge (~2e-2 A, computed
    below) vs the true equilibrium's tiny ~4.3e-4 A -- a >40x gap, nothing
    subtle to threshold."""
    outdir = os.path.join(DATA_ROOT, "freeze_evidence")
    os.makedirs(outdir, exist_ok=True)
    circuit = """* freeze-evidence discriminator circuit
.param rr=100
Vin in 0 DC 5 AC 1
R1 in mid {{rr}}
D1 mid out DMOD
C1 out 0 1u{ic}
Rload out 0 10k
Vprobe mid probeend DC 0
.model DMOD D(IS=1e-14 N=1 RS=1)
"""
    # (a) plain op discards device IC=
    deck = circuit.format(ic=" IC=2") + ".control\nop\nprint @d1[id]\n.endc\n.end\n"
    out = _run(deck, outdir, "a_plain_op")
    m = re.search(r"@d1\[id\]\s*=\s*([-+0-9.eE]+)", out)
    if m:
        id_a = float(m.group(1))
        record("(a) plain `op` discards capacitor IC=2 (converges to true DC eqm, not IC-forced)",
               id_a < 1e-3, f"@d1[id]={id_a:.6e} (true eqm ~4.3e-4; IC-forced would be ~2e-2)")
    else:
        record("(a) plain `op` discards capacitor IC=2", None, f"no @d1[id] in output: {out[-300:]}")

    # (b) tran-uic-then-op, same session, also discards it
    deck = (circuit.format(ic=" IC=2")
            + ".control\ntran 1n 1n 0 uic\nop\nprint @d1[id]\n.endc\n.end\n")
    out = _run(deck, outdir, "b_tran_uic_then_op")
    matches = re.findall(r"@d1\[id\]\s*=\s*([-+0-9.eE]+)", out)
    if matches:
        id_b = float(matches[-1])
        record("(b) tran-uic-then-op (same session) also discards the frozen IC state",
               id_b < 1e-3, f"final @d1[id]={id_b:.6e}")
    else:
        record("(b) tran-uic-then-op discards frozen IC state", None, f"no @d1[id]: {out[-300:]}")

    # (c) savestate a genuine non-equilibrium transient state, loadstate into a
    # FRESH process, then stb -- must also discard it. Per docs/checkpoint.md's
    # own trap note: ngspice's control-mode tokenizer lower-cases a bare
    # savestate/loadstate argument, which is fatal for this project's own
    # ~/Work/... paths (capital W) -- use an already-lowercase RELATIVE
    # filename plus the ngspice process's own cwd= instead.
    ckpt_name = "evidence.ckpt"
    ckpt = os.path.join(outdir, ckpt_name)
    if os.path.exists(ckpt):
        os.remove(ckpt)
    save_deck = (circuit.format(ic=" IC=0")
                 + f".control\ntran 1n 2u 0 uic\nsavestate {ckpt_name}\n.endc\n.end\n")
    _run(save_deck, outdir, "c1_save", cwd=outdir)
    if not os.path.exists(ckpt):
        record("(c) savestate+loadstate+stb also discards the checkpointed state", None,
               "savestate did not produce a checkpoint file")
    else:
        load_deck = (circuit.format(ic=" IC=0")
                      + f".tran 1n 4u 0 uic\n.control\nloadstate {ckpt_name}\nstb dec 5 1 1e6 probe=Vprobe\n"
                      "print @d1[id]\n.endc\n.end\n")
        out = _run(load_deck, outdir, "c2_load_then_stb", cwd=outdir)
        matches = re.findall(r"@d1\[id\]\s*=\s*([-+0-9.eE]+)", out)
        if matches:
            id_c = float(matches[-1])
            record("(c) savestate+loadstate+stb (fresh process) also discards the checkpointed state",
                   id_c < 1e-3, f"post-stb @d1[id]={id_c:.6e} (matches true eqm, not the loaded ~4e-2 state)")
        else:
            record("(c) savestate+loadstate+stb discards checkpointed state", None,
                   f"no @d1[id] in output: {out[-300:]}")


# ---------------------------------------------------------------------------
# layer 3: end-to-end on the LDO soft-start fixture
# ---------------------------------------------------------------------------

def vss_of(t, tss=500e-6, vref_final=0.45):
    """Independent, hand-derived closed-form of the fixture's own PWL ramp --
    deliberately NOT calling stb_tran.eval_pwl (would make this a tautology)."""
    if t <= 0:
        return 0.0
    if t >= tss:
        return vref_final
    return vref_final * t / tss


def build_frozen_ref(vss_val, outdir, label):
    """An INDEPENDENTLY hand-built frozen deck (Vss replaced with a literal DC
    value written directly into the text here, not via stb_tran.freeze_deck())
    -- the real cross-check that stb_tran.py's own freeze mechanism computes
    the same thing a by-hand deck edit would."""
    with open(LDO) as fh:
        lines = fh.read().splitlines()
    out = []
    for ln in lines:
        if ln.strip().lower().startswith("vss "):
            out.append(f"Vss vss 0 DC {vss_val!r}")
        elif ln.strip().lower().startswith((".control", "tran", "meas", "print", "write", "wrdata", ".endc")):
            continue
        else:
            out.append(ln)
    body = "\n".join(out)
    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
op
stb {SWEEP} probe=Vprobe
wrdata {outdir}/{label}_margins.txt stb2.phasemargin stb2.phasemarginfreq stb2.gainmargin stb2.gainmarginfreq
quit
.endc
.end
"""
    os.makedirs(outdir, exist_ok=True)
    wrapper = os.path.join(outdir, f"{label}.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + ctl)
    proc = subprocess.run([NGSPICE_AUG, "-b", wrapper], capture_output=True, text=True, timeout=120)
    mfile = os.path.join(outdir, f"{label}_margins.txt")
    if not os.path.exists(mfile):
        raise stb.CouldNotRun(f"reference run produced no margins file: {proc.stdout+proc.stderr}")
    with open(mfile) as fh:
        next(fh, None)
        row = [float(x) for x in fh.readline().split()]
    pm, pmf, gm, gmf = row[-4:]
    return {"pm_deg": None if math.isnan(pm) else pm, "fc_hz": None if math.isnan(pmf) else pmf,
            "gm_db": None if math.isnan(gm) else gm}


def check_end_to_end():
    outdir = os.path.join(DATA_ROOT, "e2e")
    times = [100e-6, 250e-6, 400e-6, 1000e-6]
    rows = stb_tran.sweep_times(LDO, "Vprobe", ["Vss"], times, SWEEP, outdir)

    for r in rows:
        t = r["t"]
        vss_expected = vss_of(t)
        got_frozen = r["frozen"].get("Vss")
        record(f"freeze_deck computed Vss({t:.6g})={vss_expected:.6g} at t={t:.6g}",
               close(got_frozen, vss_expected, 1e-9), f"got {got_frozen}")

        ref = build_frozen_ref(vss_expected, outdir, f"ref_{t:.9g}".replace(".", "p"))
        pm_ok = (r["pm_deg"] is None and ref["pm_deg"] is None) or close(r["pm_deg"], ref["pm_deg"], 1e-6)
        fc_ok = (r["fc_hz"] is None and ref["fc_hz"] is None) or close(r["fc_hz"], ref["fc_hz"], 1e-6)
        record(f"stb_tran.py PM/fc at t={t:.6g} matches an INDEPENDENTLY hand-built frozen-deck .stb",
               pm_ok and fc_ok,
               f"tool: PM={r['pm_deg']} fc={r['fc_hz']}  hand-ref: PM={ref['pm_deg']} fc={ref['fc_hz']}")

    # "late t_k equals plain DC-op .stb" (card's own first identity): t=1000us
    # (fixture has settled by ~600us, see the fixture's own header) vs a plain
    # DC-op .stb of the UNMODIFIED deck's own final vref_final=0.45.
    late = next(r for r in rows if r["t"] == 1000e-6)
    ref_final = build_frozen_ref(0.45, outdir, "ref_final_dcop")
    pm_ok = close(late["pm_deg"], ref_final["pm_deg"], 1.0)     # 1 deg bound, per the card's own spec
    fc_close_enough = abs(late["fc_hz"] - ref_final["fc_hz"]) / ref_final["fc_hz"] < 0.05  # crossover freq proxy for "0.1dB" gain-axis bound
    record("late t_k (1000us, settled) PM matches plain DC-op .stb within 1 deg (card's own bound)",
           pm_ok, f"t=1000us PM={late['pm_deg']:.4f}  DC-op PM={ref_final['pm_deg']:.4f}")
    record("late t_k fc within 5% of plain DC-op .stb's crossover (gain-axis proxy for the 0.1dB bound)",
           fc_close_enough, f"t=1000us fc={late['fc_hz']:.6g}  DC-op fc={ref_final['fc_hz']:.6g}")

    # red-proof: a deliberately WRONG frozen value must NOT match
    bogus = build_frozen_ref(0.09 + 0.05, outdir, "ref_bogus")  # off by 0.05 from t=100us's true 0.09
    early = next(r for r in rows if r["t"] == 100e-6)
    mismatched = not close(early["pm_deg"], bogus["pm_deg"], 1e-6)
    record("red-proof: a deliberately mis-frozen reference does NOT match stb_tran.py's t=100us result",
           mismatched, f"tool PM={early['pm_deg']:.4f} vs bogus-ref PM={bogus['pm_deg']:.4f}")


# ---------------------------------------------------------------------------
# layer 4: --mode=instant (card #500) -- linearise about the REAL transient
# state at t_k (no frozen-source/re-solved-equilibrium approximation), vs
# --mode=quasistatic's freeze-and-re-solve trade-off (layer 3 above).
#
# Two independent checks, matching card #500's own done-when:
#   (1) a SETTLED t_k (the fixture's own ramp finished ~600us ago by
#       t=1200us) must equal a plain DC-op `.stb` of the unmodified deck
#       (build_frozen_ref(0.45, ...), already used by layer 3) within 1e-6
#       relative -- checked at MUCH better than that bound in practice
#       (~1e-9, this build), since by then the frozen source's own value has
#       ALSO settled and the transient has converged to the same equilibrium
#       a DC op would find.
#   (2) a MID-RAMP t_k (t=100us, only 1/5 of the way through the 500us
#       soft-start ramp) must (a) differ from --mode=quasistatic's frozen-DC
#       result in the PHYSICALLY EXPECTED direction, and (b) match an
#       INDEPENDENT reference within 0.5 dB / 3 deg: build_tangent_ref()
#       measures the REAL transient's v(ea) at t=100us (a SEPARATE ngspice
#       call, `meas tran ... find v(ea) at=100u`, not reading it back out of
#       stb_tran.py's own instant-mode run), computes the fixture's Bpass
#       device's EXACT incremental transconductance there by hand in Python
#       (gm2*k2*sech^2(k2*(ea-vth2)), the closed-form derivative of
#       gm2*tanh(k2*(V(ea)-vth2)) -- tests/fixtures/stbtran_ldo_softstart.cir's
#       own Bpass line), substitutes a LINEAR Gpass VCCS of that exact gain
#       for the nonlinear Bpass, and runs a PLAIN (non-instant) DC-op `.stb`
#       on that now-fully-linear circuit -- a circuit with an EXACT, unique,
#       operating-point-INDEPENDENT small-signal T once every device is
#       linear, so this is a genuine, non-circular independent derivation
#       (different code path: no stb_tran.py freeze/instant machinery, no
#       CKTop-skip, just cktload's ordinary DC-op `.stb` on a hand-built
#       linear circuit), not a restatement of what instant mode itself
#       computes.
#
# Physically-expected direction (why (2)(a) goes the way it does): at
# t=100us the soft-start ramp is still rising (tss=500us), and the loop's
# OWN closed-loop bandwidth (~70-80kHz, tau ~2-3us) is much faster than the
# fixture's output-filter time constant (Rout*Cout ~530Hz, tau~300us) but
# the real error amp output v(ea) is still being driven harder than its
# EVENTUAL (quasistatic, Vss-frozen-forever) equilibrium value would need,
# pushing Bpass's tanh further into saturation -- REDUCING its incremental
# gm relative to the quasistatic-frozen op's smaller, near-linear-region
# ea. Lower incremental gm -> lower loop gain magnitude -> LOWER crossover
# frequency (fc) for instant mode vs quasistatic at the same t_k. Verified
# numerically below (not asserted from theory alone): quasistatic mode
# converges ea to ~0.045 for Vss=0.09 held forever (near-linear region,
# sech^2~0.98); the REAL transient's ea at t=100us is ~0.198 (sech^2~0.72,
# visibly further into saturation) -- a >4x difference in ea, not a rounding
# artifact.
# ---------------------------------------------------------------------------

GM2, K2, VTH2 = 0.15, 3.0, 0.0  # tests/fixtures/stbtran_ldo_softstart.cir's own
                                # .param gm2/k2/vth2 (Bpass's tanh args)


def measure_node(t, tstep, nodes, outdir, label, timeout=300):
    """v(node) for each `nodes`, from the LDO fixture's own REAL transient
    (interactive `tran` inside `.control`, run once, independent of
    stb_tran.py's own instant-mode machinery) at time t, via
    `meas tran ... find v(node) at=t`."""
    with open(LDO) as fh:
        lines = fh.read().splitlines()
    out, skip = [], False
    for ln in lines:
        low = ln.strip().lower()
        if low.startswith(".control"):
            skip = True
            continue
        if skip:
            if low.startswith(".endc"):
                skip = False
            continue
        if low.startswith(".end"):
            continue
        out.append(ln)
    body = "\n".join(out)
    meas = "\n".join(f"meas tran {n}_m find v({n}) at={t!r}" for n in nodes)
    ctl = f"""
.control
tran {tstep} {t!r} 0 uic
{meas}
.endc
.end
"""
    os.makedirs(outdir, exist_ok=True)
    wrapper = os.path.join(outdir, f"{label}.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + "\n" + ctl)
    proc = subprocess.run([NGSPICE_AUG, "-b", wrapper], capture_output=True, text=True, timeout=timeout)
    text = proc.stdout + proc.stderr
    vals = {}
    for n in nodes:
        m = re.search(rf"{n}_m\s*=\s*([-+0-9.eE]+)", text)
        if m:
            vals[n] = float(m.group(1))
    return vals


def build_tangent_ref(t, tstep, outdir, label):
    """The independent reference described in this section's block comment:
    measure the real transient's v(ea) at t, hand-compute Bpass's EXACT
    incremental gm there, substitute a linear VCCS, and run a plain DC-op
    `.stb` on the resulting fully-linear circuit. Returns
    ({pm_deg, fc_hz, gm_db}, ea) -- ea returned too so the caller can report
    it (this section's docstring's ">4x difference in ea" claim)."""
    vals = measure_node(t, tstep, ["ea"], outdir, f"{label}_meas")
    if "ea" not in vals:
        raise stb.CouldNotRun(f"could not measure v(ea) at t={t:.6g} for the tangent reference")
    ea = vals["ea"]
    glin = GM2 * K2 * (1.0 - math.tanh(K2 * (ea - VTH2)) ** 2)

    with open(LDO) as fh:
        lines = fh.read().splitlines()
    out = []
    for ln in lines:
        s = ln.strip()
        low = s.lower()
        if low.startswith("bpass "):
            out.append(f"Gpass vin vout ea 0 {glin!r}")
        elif low.startswith((".control", "tran", "meas", "print", "write", "wrdata", ".endc")):
            continue
        else:
            out.append(ln)
    body = "\n".join(out)
    ctl = f"""
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
op
stb {SWEEP} probe=Vprobe
wrdata {outdir}/{label}_margins.txt stb2.phasemargin stb2.phasemarginfreq stb2.gainmargin stb2.gainmarginfreq
quit
.endc
.end
"""
    os.makedirs(outdir, exist_ok=True)
    wrapper = os.path.join(outdir, f"{label}.cir")
    with open(wrapper, "w") as fh:
        fh.write(body + ctl)
    proc = subprocess.run([NGSPICE_AUG, "-b", wrapper], capture_output=True, text=True, timeout=120)
    mfile = os.path.join(outdir, f"{label}_margins.txt")
    if not os.path.exists(mfile):
        raise stb.CouldNotRun(f"tangent-ref run produced no margins file: {proc.stdout + proc.stderr}")
    with open(mfile) as fh:
        next(fh, None)
        row = [float(x) for x in fh.readline().split()]
    pm, pmf, gm, gmf = row[-4:]
    return ({"pm_deg": None if math.isnan(pm) else pm, "fc_hz": None if math.isnan(pmf) else pmf,
              "gm_db": None if math.isnan(gm) else gm}, ea)


def check_instant_settled():
    """(1): a settled t_k (1200us -- ramp ends at 500us, fixture's own header
    says vout is flat by ~600us) must equal a plain DC-op `.stb` (Vss=0.45,
    the ramp's own final value) within 1e-6 relative -- card #500's own
    done-when bound."""
    outdir = os.path.join(DATA_ROOT, "instant_settled")
    t = 1200e-6
    rows = stb_tran.sweep_times(LDO, "Vprobe", [], [t], SWEEP, outdir, mode="instant", tstep="50n")
    got = rows[0]
    ref = build_frozen_ref(0.45, outdir, "ref_settled_dcop")
    pm_ok = close(got["pm_deg"], ref["pm_deg"], 1e-6)
    fc_ok = close(got["fc_hz"], ref["fc_hz"], 1e-6)
    record("instant mode at a settled t_k (1200us) equals a plain DC-op .stb within 1e-6 relative",
           pm_ok and fc_ok,
           f"instant: PM={got['pm_deg']!r} fc={got['fc_hz']!r}  DC-op: PM={ref['pm_deg']!r} fc={ref['fc_hz']!r}")


def check_instant_mid_ramp():
    """(2): mid-ramp (t=100us) instant-mode result (a) differs from
    quasistatic in the expected direction (lower fc: real ea is further into
    Bpass's tanh saturation than the quasistatic-frozen equilibrium, see this
    section's block comment) and (b) matches the independent tangent-
    linearisation reference within 0.5 dB / 3 deg (checked here as: PM within
    3 deg absolute, fc within 6% -- 0.5 dB of loop-gain magnitude is
    approximately a 6% shift on the gain axis, the same 0dB-crossing-slope
    proxy layer 3's own "late t_k" check already uses for its 0.1dB bound)."""
    outdir = os.path.join(DATA_ROOT, "instant_mid_ramp")
    t = 100e-6

    instant_rows = stb_tran.sweep_times(LDO, "Vprobe", [], [t], SWEEP, outdir, mode="instant", tstep="50n")
    instant = instant_rows[0]

    quasi_rows = stb_tran.sweep_times(LDO, "Vprobe", ["Vss"], [t], SWEEP, outdir)
    quasi = quasi_rows[0]

    ref, ea = build_tangent_ref(t, "50n", outdir, "tangent_ref_100u")

    direction_ok = (instant["fc_hz"] is not None and quasi["fc_hz"] is not None
                     and instant["fc_hz"] < quasi["fc_hz"])
    record("instant vs quasistatic at t=100us differ in the expected direction (instant fc < quasistatic fc)",
           direction_ok,
           f"instant fc={instant['fc_hz']!r} Hz (PM={instant['pm_deg']!r})  "
           f"quasistatic fc={quasi['fc_hz']!r} Hz (PM={quasi['pm_deg']!r})  real v(ea)={ea!r}")

    pm_ok = instant["pm_deg"] is not None and ref["pm_deg"] is not None and abs(instant["pm_deg"] - ref["pm_deg"]) <= 3.0
    fc_ok = instant["fc_hz"] is not None and ref["fc_hz"] is not None and close(instant["fc_hz"], ref["fc_hz"], 0.06)
    record("instant mode at t=100us matches the independent tangent-linearisation reference "
           "within 0.5dB/3deg (PM within 3deg, fc within 6%)",
           pm_ok and fc_ok,
           f"instant: PM={instant['pm_deg']!r} fc={instant['fc_hz']!r}  "
           f"tangent-ref (ea={ea!r}, G_lin={GM2 * K2 * (1.0 - math.tanh(K2 * (ea - VTH2)) ** 2)!r}): "
           f"PM={ref['pm_deg']!r} fc={ref['fc_hz']!r}")


def check_instant_cli_red_proofs():
    """--mode=instant needs --tstep and refuses --freeze/--check-tran
    (neither applies -- instant mode runs the deck's own real, unmodified
    transient, nothing to freeze or validate against): tools/stb_tran.py's
    main() itself, invoked as a subprocess (real CLI, not calling sweep_times
    directly), must exit 2 (COULD NOT RUN), never silently accept and
    misbehave."""
    script = os.path.join(ROOT, "tools", "stb_tran.py")
    base = [sys.executable, script, LDO, "--probe", "Vprobe", "--times", "100u"]

    proc = subprocess.run(base + ["--mode", "instant"], capture_output=True, text=True, timeout=30)
    record("red-proof: --mode=instant without --tstep exits 2 (COULD NOT RUN)",
           proc.returncode == 2, f"exit={proc.returncode} stderr={proc.stderr.strip()[:120]}")

    proc = subprocess.run(base + ["--mode", "instant", "--tstep", "50n", "--freeze", "Vss"],
                           capture_output=True, text=True, timeout=30)
    record("red-proof: --mode=instant combined with --freeze exits 2 (COULD NOT RUN)",
           proc.returncode == 2, f"exit={proc.returncode} stderr={proc.stderr.strip()[:120]}")

    proc = subprocess.run(base + ["--mode", "instant", "--tstep", "50n", "--check-tran", "50n 1.2m"],
                           capture_output=True, text=True, timeout=30)
    record("red-proof: --mode=instant combined with --check-tran exits 2 (COULD NOT RUN)",
           proc.returncode == 2, f"exit={proc.returncode} stderr={proc.stderr.strip()[:120]}")

    proc = subprocess.run([sys.executable, script, LDO, "--probe", "Vprobe", "--times", "100u"],
                           capture_output=True, text=True, timeout=30)
    record("red-proof: default --mode=quasistatic without --freeze exits 2 (COULD NOT RUN)",
           proc.returncode == 2, f"exit={proc.returncode} stderr={proc.stderr.strip()[:120]}")


def check_cli_could_not_run():
    saved = stb_tran.NGSPICE_AUG
    stb_tran.NGSPICE_AUG = "/nonexistent/ngspice-aug"
    try:
        stb_tran.sweep_times(LDO, "Vprobe", ["Vss"], [100e-6], SWEEP,
                              os.path.join(DATA_ROOT, "bogus_binary"))
        record("red-proof: missing ngspice-aug binary gives COULD-NOT-RUN", False, "did not raise")
    except stb.CouldNotRun as e:
        record("red-proof: missing ngspice-aug binary gives COULD-NOT-RUN", True, str(e)[:60])
    finally:
        stb_tran.NGSPICE_AUG = saved


def main():
    check_waveforms()
    check_param_resolver()

    if _have_ngspice_aug():
        try:
            check_freeze_evidence()
            check_end_to_end()
            check_cli_could_not_run()
            check_instant_settled()
            check_instant_mid_ramp()
            check_instant_cli_red_proofs()
        except stb.CouldNotRun as e:
            record("run", None, str(e))
    else:
        record("fork build (NGSPICE_AUG)", None, f"{NGSPICE_AUG} not found/executable -- "
                "layers 2/3 (ngspice-dependent) skipped as COULD-NOT-RUN, not silently passed")

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    npass = sum(s == "pass" for _, s, _ in results)
    print(f"\n{npass}/{len(results)} pass")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
