#!/usr/bin/env python3
"""Validation of the generic charge-based history of the time-domain PAC engine (card #542, 2026-10-02):
`.pac/.pstb/.pxf/.pnoise ... solver=td` on circuits whose devices carry MODEL charges (MOS1 Meyer/overlap/junction,
BSIM4 intrinsic/overlap/junction, diode cj/tt, BJT, mutual inductors).  Engine: src/spicelib/analysis/pacshoot.c
(pac_generic_prepare(), pac_gen_loads(), pac_step_generic()); design: docs/pss-pac.md "Generic charge history (card #542)".

    NGSPICE=<prefix>/bin/ngspice-aug [NGSPICE_ASAN=<asan prefix>/bin/ngspice] python3 tests/test_pac_charge.py [--quick]
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN (a FAIL wins)

Outputs under ~/Work/data/ngspice_augmentation/544/test (never the repo; PAC_CHARGE_OUT=<dir> redirects them, card #713).  Every check can fail: red proofs at the end set
PAC_NATIVE_HISTORY=1 (the pre-#542 native-only history) and require the checks' own metric to FAIL.

Checks
  1. LTI vs plain `.ac` (the gold reference for any LTI circuit): each deck has a slow periodic clock on an unrelated net so
     `.pss` has a period, the bias is fixed, so the periodic small-signal response of the DUT IS its `.ac`.  Bound 1e-6 relative
     (worst point of 1 kHz..100 kHz, pacpoints 400000).  MOS1 gate-capacitance ladder, MOS1 amplifier (Meyer+overlap+junction),
     diode (cj+tt), BJT amplifier, JFET, MOS3, BSIM3v3.3, BSIM4 common source capmod 0 / capmod 2 / rgatemod 1, transformer with a K statement.
  2. bug 1 reproducer (tests/bugs/539/bug1_pac_model_charge): chain of N inverters, control-to-average-output gain vs the static
     derivative of the periodic steady state (two `.tran`, reltol 1e-5), N = 1..3, MOS1 + overlap caps and BSIM4.  Bound 1 %.
  3. `.pxf` adjoint=1 vs adjoint=0 (forward) and `.pnoise` vs `.noise` on the same decks (the transposed recursion must use
     exactly the matrices the forward one used): 1e-8 / 1e-5 relative.
  4. Honesty guard: PAC_GEN_TEST_NONAFFINE=BSIM4 makes that type look non-affine in ag0 -> must REFUSE naming the type; BSIM4
     trnqsmod=1 (card #559: now accepted, see 12.); free-running + charged diode must REFUSE; an ideal
     switch (saltation) next to a charged diode must REFUSE; a diode without any charge next to the same switch must still run
     on the native path (no refusal).
  5. Cross-checks: native-only deck forced through the generic path (PAC_FORCE_GENERIC=1) equals the native path (1e-6); `.options
     klu` equals sparse (1e-9) on a MOS deck.
  6. Red proofs (PAC_NATIVE_HISTORY=1 must FAIL checks 1 and 2).
  7. ASan (only if NGSPICE_ASAN is set): generic `.pac` + `.pxf adjoint=1` + `.pnoise` on a BSIM4 chain, no sanitizer report.

Card #544 (2026-10-02) additions
  8. `.pnoise dev=1` (total and every per-device contribution) vs plain `.noise` on LTI-biased decks: MOS1 (with kf), BSIM4, diode, BJT, JFET,
     MOS3, BSIM3 (check_pnoise_dev; bound 1e-5, measured <= 9e-7).
  9. `.pnoise` under a TIME-VARYING bias vs the quasi-static time average of `.noise` over LO phases (MOS1, BSIM4, diode, BJT; check_pnoise_qs;
     bound 5e-4, measured 1.5e-7..3e-4).  Red proof PAC_NOISE_STALE_STATE=1 (the pre-#544 device operating point: stale CKTstate0) must fail.
 9b. card #556: the quasi-static reference of 9 is the white-noise answer; a flicker source whose density is modulated at the 100 kHz LO
     folds as  Sum_q |a_q|^2 S_w(f - q fs)  (a_q: DFT of the per-phase amplitude), which the check now uses for the flicker part (the engine
     before #556 evaluated S_w(f) for every sideband and agreed with the plain time average: 3-31 % off the folded value for the diode,
     BJT and BSIM4 decks, which have kf).  Red proof PNOISE_FLICKER_OLD=1 (the old evaluation) must fail the diode check.
 10. free-running transistor-level (BSIM4) hysteretic buck R8: `.pac solver=td` control-to-output and line-to-output at 20/50/80 kHz vs
     tools/pac_perturbed.py --freerun (0.1 dB / 0.5 deg; measured <= 0.021 dB / 0.052 deg); R5 + charged diode runs on the generic path and equals
     the native R5 value.  Red proof PAC_NATIVE_HISTORY=1.
 11. ideal-switch (saltation) events with model charges on the orbit: R6o with a BSIM4 source follower in the control path vs F(.ac) x pac_perturbed(R6o)
     and the exact factorisation identity; red proof PAC_NATIVE_HISTORY=1 (28 dB off).  A charged diode at the switched node is still refused.

Card #559 (2026-10-03) additions (checks_559(), red_proofs_559(); `--only 559`)
 12. ideal-switch (saltation) events at a switch node with device-model charge AND a native capacitor (R6od: freewheel diode cj+tt; R6odm: MOS body
     diode + cgdo with a floating-through-cgdo gate): .pac td ctrl2out at 20/50/100 kHz vs tools/pac_perturbed (0.1 dB / 0.5 deg; measured <= 0.021 dB /
     0.16 deg); native decks: see card #568 below; still refused: model charge only (no native capacitor) at the node, mutual inductance on the routed
     branch.  Red proofs PAC_SALT_NATIVEJUMP=1 (the old jump: 0.6 dB / 4.0 dB off; 4.8 / 0.3 dB before card #568 shortened the settle window) and PAC_SALT_NOQ=1 (the refusal back).
 13. non-quasi-static MOS channel (BSIM4 trnqsmod=1, BSIM3 nqsmod=1) in td PAC: LTI bias vs a converged transient FRA (1e-4 relative at 1..30 MHz; .ac is
     not the reference: acnqsmod is a different approximation), an LO-switched BSIM4 stage vs transient FRA (0.01 dB), .pxf adjoint == forward, .pnoise on
     NQS accepted since card #715 (tests/test_pnoise_nqs.py holds its accuracy checks).  An NQS flag on the device line is found too (card #730; BSIM4 4.7 / BSIM3 3.2.4 refused, BSIM4 4.8 accepted).  Red proofs: PAC_NQS_REFUSE=1, PAC_NQS_PNOISE_REFUSE=1, PAC_NQS_INSTANCE_OFF=1; a quasi-static deck fails the NQS references (40 % / 3 dB).

Card #568 (2026-10-06)
 14. native saltation with a chargeless diode at the switch node (R6od, cj = tt = 0): the order-1 settle window after an event is bounded to 2 accepted
     steps (it ran 10..92 steps and the result was 5.3 dB off, silently, up to s115); vs tools/pac_perturbed at 20/50/100 kHz (0.1 dB / 0.5 deg; measured
     0.026 dB / 0.16 deg); R6o within 2e-4 of its s115 value.  Red proof PSS_SALT_MAXSTEPS=-1 (the old window).

Card #575 (2026-10-06) (checks_575(), red_proofs_575(); `--only 575`)
 15. ideal-switch events at an inductor node without a dominant native capacitor: fast-node treatment (routed jump continued until the node
     has settled) on 6 decks, charge jump where the native capacitor dominates or sits with a snubber, chargeless nodes with the exact 2-step
     settle window and the snubber-capacitor jump; each vs its own perturbed-PSS reference at 20/50/100 kHz.  A charge-only node with an RC
     snubber stayed refused (until card #713, item 21).  Red proofs: PSS_SALT_NOCONT, PAC_SALT_QONLY, PSS_SALT_MAXSTEPS=-1, PSS_SALT_NOSNUB.

Card #579 (2026-10-06) (checks_579(), red_proofs_579(); `--only 579`)
 16. the same on other circuits: boost (chargeless, cj, cj + tt), synchronous buck with dead time (chargeless and charged body diode), line-to-
     output, duty 0.7; each vs its own perturbed-PSS reference.  Red proofs: PSS_SALT_NO579, PSS_SALT_NODIOFAR, PAC_NO_IMPULSE, PSS_SALT_NOCONT.

Cards #580 / #581 (2026-10-07) (checks_580_581(), red_proofs_580_581(); `--only 580`)
 17. boost with a capacitor at the switch node (fast-node treatment, red proof PSS_SALT_NOFARFAST); discontinuous conduction: converged at
     2000 points per period, the DCM note printed (and absent in continuous conduction); since #711 (2026-10-19) 200 points are within
     0.06 dB, red proof: PSS_SALT_NO711=1 at 200 points is > 0.3 dB off (it was the plain 200-point run before #711).

Card #582 (2026-10-07) (checks_582(), red_proofs_582(); `--only 582`)
 18. half bridge without dead time and without node capacitor (two switches flipping in one step = one event; generic recursion for capacitor-
     less switch nodes): open loop with a chargeless / charged / no body diode, closed-loop .pstb vs the loop-gain reference at 30 kHz, and a
     closed-loop fast-node deck.  Red proofs: PSS_SALT_NO582, PAC_NO582_GENERIC.

Card #583 (2026-10-07) (checks_583(), red_proofs_583(); `--only 583`)
 19. routed jumps added after T (no order-1 window): inverting buck-boost with and without a node capacitor, boost with 3 nF, and tighter phase
     bounds on three fast-node buck decks.  Red proof PAC_SALT_SMOOTHT=0 (phase).

Card #586 (2026-10-07) (checks_586(), red_proofs_586(); `--only 586`)
 20. a source that is not periodic from t = 0 (delayed sawtooth) is refused by the shooting engine; a delayed pulse that fits in its period is
     not; a two-phase interleaved buck with ramps periodic from 0 vs the perturbed PSS.  Red proof PSS_ALLOW_DELAYED_SOURCE=1.

Card #713 (2026-10-19) (checks_713(), red_proofs_713(); `--only 713`)
 21. a switch node whose only charge is the devices' (diode cj / tt, MOS body diode) TOGETHER with an RC snubber is no longer refused: fast
     node (routed jump + continuation) and the snubber capacitor's jump G / C_far x (D_end - D_before), D = v_node - v_far, when the node's
     transition has ended; every continuation jump takes the event coefficient of its switch flip (not of the step it is given in).
     9 open-loop decks (snubber 10..50 ohm, 100 pF..10 nF, duty 0.275 / 0.7, MOS body diode, cj 50 pF + tt 5 ns) and a DCM deck (2000 points)
     vs their perturbed-PSS references at 20/50/100 kHz, two at reltol 1e-6, and a closed-loop .pstb (asynchronous buck, charged diode,
     snubber) vs pac_perturbed --loop-gain at 10/30/100 kHz; the s123 bar (0.067 dB / 0.32 deg).  Red proofs: PSS_SALT_NO713=1 (the refusal
     back), PAC_SALT_FASTSNUB=1 (no snubber jump: the gain error of the slow-snubber decks grows 3-10x), PAC_NO713_COEF=1 (continuation
     coefficient at its own step: the DCM deck 0.47 deg off at 100 kHz).
"""
import cmath
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NG = os.environ.get("NGSPICE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice-aug")
ASAN = os.environ.get("NGSPICE_ASAN")
OUT = os.path.expanduser("~/Work/data/ngspice_augmentation/544/test")   # the r8_* fixtures write here too; card #559 outputs share it
OUT_FIXED = OUT                                    # the folder the fixtures' own wrdata lines name
if os.environ.get("PAC_CHARGE_OUT"):               # card #713: a second run beside another one (PAC_CHARGE_OUT=<dir>; the fixtures are rewritten into it)
    OUT = os.path.abspath(os.path.expanduser(os.environ["PAC_CHARGE_OUT"]))
QUICK = "--quick" in sys.argv
results = []


class CouldNotRun(Exception):
    pass


QUIET = False


def record(name, state, detail):
    results.append((name, state, detail))
    if not QUIET:
        print(f"[{state:13s}] {name}: {detail}", flush=True)


def run(deck, name, env=None, timeout=1500, exe=None):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name + ".cir")
    if OUT != OUT_FIXED:
        deck = deck.replace(OUT_FIXED, OUT)   # card #713: PAC_CHARGE_OUT
    with open(p, "w") as fh:
        fh.write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    try:
        r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0",
                            exe or NG, "-b", p], capture_output=True, text=True, cwd=OUT, env=e, timeout=timeout)
    except (FileNotFoundError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{name}: {ex}")
    return (r.stdout + r.stderr).replace("\r", "\n"), r.returncode


def rd(fn, cplx=True):
    path = os.path.join(OUT, fn)
    if not os.path.exists(path):
        raise CouldNotRun(f"missing {fn}")
    rows = []
    for l in open(path).read().split("\n")[1:]:
        t = l.split()
        if len(t) >= 2:
            rows.append((float(t[0]), complex(float(t[1]), float(t[2])) if cplx else float(t[1])))
    if not rows:
        raise CouldNotRun(f"empty {fn}")
    return rows


def worst_rel(a, b):
    if len(a) != len(b) or any(abs(x[0] - y[0]) > 1e-6 * x[0] for x, y in zip(a, b)):
        raise CouldNotRun("frequency lists differ (a sweep compared point by point must have identical frequencies)")
    return max(abs(x[1] - y[1]) / abs(x[1]) for x, y in zip(a, b))


# slow clock on an unrelated net: gives .pss a period (1 ms) without touching the bias of the circuit under test.
CLK = ("* pac charge lti deck\nVclk clk 0 PULSE(0 1 0 1n 1n 500u 1m)\nRclk clk clkgnd 100k\nCclk clkgnd 0 10n\n")
OPTS = ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n"
BSIM = ("nch4 nmos level=54 version=4.8.2 mobmod=0 capmod=%d igcmod=0 igbmod=0 rdsmod=0 rbodymod=0 rgatemod=%d acnqsmod=0 "
        "toxe=3e-9 toxp=3e-9 toxm=3e-9 vth0=0.5 k1=0.4 k2=0 k3=0 w0=2.5e-6 u0=300 vsat=1.1e5 cgso=2e-10 cgdo=2e-10 cgbo=1e-11 "
        "cgsl=1e-10 cgdl=1e-10 cf=1e-10 cjs=1e-3 cjd=1e-3 cjswd=1e-10 cjswgd=1e-10 xpart=0 %s")
BSIM_AMP = ("Vdd vdd 0 3\nVin in 0 DC 0.9 AC 1\nRs in g 20k\nM1 out g 0 0 nch4 w=4u l=0.5u as=4p ad=4p ps=10u pd=10u\n"
            "Rd vdd out 20k\nCl out 0 50f\n")
DECKS = {
    # MOS1 gate capacitances only (cutoff): Meyer cgb + overlap, an RC ladder made of MOS gates
    "mos1_ladder": CLK + ".model nm nmos level=1 vto=0.7 kp=2e-4 lambda=0.02 tox=2e-8 cgso=5e-10 cgdo=5e-10 cgbo=2e-10 cj=1e-4\n"
                    "Vin in 0 DC 0.2 AC 1\nR1 in g1 10Meg\nM1 g1 0 0 0 nm w=100u l=2u\nR2 g1 g2 10Meg\nM2 g2 0 0 0 nm w=100u l=2u\n"
                    "R3 g2 out 10Meg\nM3 out 0 0 0 nm w=100u l=2u\n",
    "mos1_amp": CLK + ".model nm nmos level=1 vto=0.7 kp=2e-4 lambda=0.02 tox=2e-8 cgso=5e-10 cgdo=5e-10 cgbo=2e-10 cj=1e-4 "
                "cjsw=1e-10 pb=0.8\nVdd vdd 0 5\nVin in 0 DC 1.2 AC 1\nRs in g 100k\n"
                "M1 out g 0 0 nm w=40u l=2u as=1e-10 ad=1e-10 ps=1e-4 pd=1e-4\nRd vdd out 20k\nCl out 0 10p\n",
    "diode": CLK + ".model dm d is=1e-14 cjo=2p vj=0.7 m=0.5 tt=5n\nVin in 0 DC 2 AC 1\nRs in a 10k\nD1 a 0 dm\nRl a out 1k\nCl out 0 1n\n",
    "bjt_amp": CLK + ".model qm npn is=1e-15 bf=100 vaf=80 cje=1p cjc=0.5p tf=0.2n tr=5n\nVcc vcc 0 5\nVin in 0 DC 1 AC 1\n"
                "Rb in b 200k\nRc vcc out 4k\nQ1 out b 0 0 qm\nCl out 0 1p\n",
    # same with a collector-substrate junction (cjs): its charge history has no device-level ic to re-derive (see docs), so a
    # startup remainder of ~7e-6 stays; documented, bounded at 2e-5 here
    "bjt_amp_cjs": CLK + ".model qm npn is=1e-15 bf=100 vaf=80 cje=1p cjc=0.5p tf=0.2n tr=5n cjs=1p\nVcc vcc 0 5\nVin in 0 DC 1 AC 1\n"
                "Rb in b 200k\nRc vcc out 4k\nQ1 out b 0 0 qm\nCl out 0 1p\n",
    "bsim4_cs_capmod0": CLK + ".model " + BSIM % (0, 0, "") + "\n" + BSIM_AMP,
    "bsim4_cs_capmod2": CLK + ".model " + BSIM % (2, 0, "") + "\n" + BSIM_AMP,
    "bsim4_cs_rgate": CLK + ".model " + BSIM % (2, 1, "rshg=50 ngcon=1") + "\n" + BSIM_AMP,
    "jfet_amp": CLK + ".model jm njf vto=-2 beta=1e-3 lambda=0.01 cgs=2p cgd=1p pb=0.8\nVdd vdd 0 10\nVin in 0 DC -1 AC 1\nRs in g 100k\n"
                "J1 out g 0 jm\nRd vdd out 5k\nCl out 0 10p\n",
    "mos3_amp": CLK + ".model nm3 nmos level=3 vto=0.7 kp=2e-4 tox=2e-8 cgso=5e-10 cgdo=5e-10 cj=1e-4 theta=0.05 eta=0.1 kappa=0.2\n"
                "Vdd vdd 0 5\nVin in 0 DC 1.2 AC 1\nRs in g 100k\nM1 out g 0 0 nm3 w=40u l=2u as=1e-10 ad=1e-10 ps=1e-4 pd=1e-4\n"
                "Rd vdd out 20k\nCl out 0 10p\n",
    "bsim3_amp": CLK + ".model nb3 nmos level=8 version=3.3.0 tox=5e-9 vth0=0.5 cgso=2e-10 cgdo=2e-10 cj=1e-3 capmod=3\nVdd vdd 0 3\n"
                 "Vin in 0 DC 0.9 AC 1\nRs in g 20k\nM1 out g 0 0 nb3 w=4u l=0.5u as=4p ad=4p ps=10u pd=10u\nRd vdd out 20k\nCl out 0 50f\n",
    "mutual_xfmr": CLK + "Vin in 0 DC 0 AC 1\nRs in p 100\nL1 p 0 1m\nL2 s 0 4m\nK1 L1 L2 0.9\nRl s out 50\nCl out 0 100n\n",
}
SWEEP = "dec 5 1k 100k"
BOUNDS = {"bjt_amp_cjs": 2e-5, "jfet_amp": 1e-5}


def wrline(name, vec):
    return (".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata %s %s\nquit\n.endc\n.end\n" % (name, vec))


def ac_deck(name):
    return DECKS[name] + OPTS + ".ac %s\n" % SWEEP + wrline("ac_%s.txt" % name, "v(out)")


def pac_deck(name, pts, extra=""):
    return (DECKS[name] + OPTS + ".pac 1k 0 clkgnd 2000 1 8 1e-6 %s method=\"shoot\" solver=\"td\" pacpoints=%d %s\n" % (SWEEP, pts, extra)
            + wrline("td_%s.txt" % name, "out"))


# ---------------------------------------------------------------------------------------------------------------------
def check_lti(name, bound=1e-6, pts=400000, env=None, tag=""):
    t0 = time.time()
    try:
        run(ac_deck(name), "ac_" + name)
        txt, rc = run(pac_deck(name, pts), "td_" + name, env=env)
        a, t = rd("ac_%s.txt" % name), rd("td_%s.txt" % name)
        w = worst_rel(a, t)
    except CouldNotRun as ex:
        record(f"LTI vs .ac: {name}{tag}", "COULD-NOT-RUN", str(ex))
        return None
    m = re.search(r"-> (\d+) state rows", txt)
    ok = w <= bound
    record(f"LTI vs .ac: {name}{tag}", "pass" if ok else "FAIL",
           f"worst rel {w:.3e} (bound {bound:g}), {len(a)} points 1k..100k, state rows {m.group(1) if m else 'native path'}, "
           f"|H|(1k)={abs(a[0][1]):.4g} |H|(100k)={abs(a[-1][1]):.4g}, wall {time.time() - t0:.0f}s")
    return w


def load_repro():
    import importlib.util
    p = os.path.join(HERE, "bugs", "539", "bug1_pac_model_charge", "run_repro.py")
    spec = importlib.util.spec_from_file_location("run_repro", p)
    m = importlib.util.module_from_spec(spec)
    os.environ["NGSPICE"] = NG
    spec.loader.exec_module(m)
    m.NG = NG
    m.W = os.path.join(OUT, "repro")
    os.makedirs(m.W, exist_ok=True)
    return m


def static_gain(r, stages, models, dv=0.01):
    av = []
    for vc in (2.5 + dv, 2.5 - dv):
        o = r.run("st", r.deck(stages, models, vc="%g" % vc, analysis=".options reltol=1e-5\n.control\ntran 1n 60u 0 uic\n"
                                "meas tran av AVG v(out) FROM=40u TO=60u\n.endc\n.end\n"))
        mm = re.search(r"^av\s*=\s*([-0-9.e+]+)", o, re.M)
        if not mm:
            raise CouldNotRun("static tran gave no av")
        av.append(float(mm.group(1)))
    return (av[0] - av[1]) / (2 * dv)


def check_chain(models=("level1_cg1e-8", "bsim4"), stages=(1, 2, 3), bound=0.01, env=None, tag="", pts=40000):
    try:
        r = load_repro()
    except Exception as ex:
        record("bug1 inverter chain", "COULD-NOT-RUN", str(ex))
        return None
    worst = 0.0
    parts = []
    for mn in models:
        for st in stages:
            try:
                s = abs(static_gain(r, st, r.MODELS[mn]))
                old = dict(os.environ)
                os.environ.update(env or {})
                try:
                    g = abs(r.pac_gain(st, r.MODELS[mn], pts, env=env))
                finally:
                    os.environ.clear(); os.environ.update(old)
            except Exception as ex:
                record(f"bug1 chain {mn} N={st}{tag}", "COULD-NOT-RUN", str(ex))
                continue
            e = abs(g - s) / s
            worst = max(worst, e)
            parts.append(f"{mn} N={st}: static {s:.4f} td {g:.4f} ({100 * e:.2f}%)")
    record(f"bug1 chain: td .pac vs static d<v>/dVc{tag}", "pass" if worst <= bound else "FAIL",
           f"worst {100 * worst:.2f}% (bound {100 * bound:g}%); " + "; ".join(parts))
    return worst


def check_pxf_adjoint(name, bound=1e-8):
    try:
        res = []
        for adj in (0, 1):
            d = (DECKS[name] + OPTS + '.pxf 1k 0 clkgnd 2000 1 8 1e-6 out %s sources="Vin" method="shoot" solver="td" adjoint=%d '
                 'pacpoints=100000\n' % (SWEEP, adj) + wrline("px%d_%s.txt" % (adj, name), "pxf_vin"))
            run(d, "px%d_%s" % (adj, name))
            res.append(rd("px%d_%s.txt" % (adj, name)))
        w = worst_rel(res[0], res[1])
    except CouldNotRun as ex:
        record(f".pxf adjoint vs forward: {name}", "COULD-NOT-RUN", str(ex))
        return
    record(f".pxf adjoint=1 vs adjoint=0: {name}", "pass" if w <= bound else "FAIL", f"worst rel {w:.3e} (bound {bound:g})")


def check_pnoise_lti(name="mos1_ladder", bound=1e-5):
    try:
        nz = DECKS[name] + OPTS + ".noise v(out) Vin %s\n" % SWEEP + wrline("nz_%s.txt" % name, "onoise_spectrum").replace("run\n", "run\nsetplot noise1\n")
        pn = (DECKS[name] + OPTS + '.pnoise 1k 0 clkgnd 2000 1 8 1e-6 out Vin %s cyclo method="shoot" solver="td" pacpoints=200000\n' % SWEEP
              + wrline("pn_%s.txt" % name, "onoise_spectrum"))
        run(nz, "nz_" + name)
        run(pn, "pn_" + name)
        a, t = rd("nz_%s.txt" % name, False), rd("pn_%s.txt" % name, False)
        w = max(abs(x[1] - y[1]) / abs(x[1]) for x, y in zip(a, t))
        assert len(a) == len(t)
    except (CouldNotRun, AssertionError) as ex:
        record(f".pnoise vs .noise: {name}", "COULD-NOT-RUN", str(ex))
        return
    record(f".pnoise td vs plain .noise (LTI): {name}", "pass" if w <= bound else "FAIL", f"worst rel {w:.3e} (bound {bound:g})")


# --- card #544 item 1: .pnoise absolute level with model-charge devices vs plain .noise (LTI bias) -----------------------------
NDECKS = {
    "mos1_amp_kf": DECKS["mos1_amp"].replace("cjsw=1e-10 pb=0.8", "cjsw=1e-10 pb=0.8 kf=1e-24 af=1"),
    "bsim4_cs_noise": CLK + ".model " + BSIM % (2, 0, "") + "\n" + BSIM_AMP,
    "diode_noise": CLK + ".model dm d is=1e-14 cjo=2p vj=0.7 m=0.5 tt=5n rs=50 kf=1e-12 af=1\nVin in 0 DC 2 AC 1\nRs in a 10k\nD1 a 0 dm\nRl a out 1k\nCl out 0 1n\n",
    "bjt_amp_noise": CLK + ".model qm npn is=1e-15 bf=100 vaf=80 cje=1p cjc=0.5p tf=0.2n tr=5n rb=100 re=2 kf=1e-14 af=1\nVcc vcc 0 5\nVin in 0 DC 1 AC 1\n"
                    "Rb in b 200k\nRc vcc out 4k\nQ1 out b 0 0 qm\nCl out 0 1p\n",
    "jfet_amp_noise": DECKS["jfet_amp"].replace("pb=0.8", "pb=0.8 kf=1e-18 af=1"),
    "mos3_amp_noise": DECKS["mos3_amp"].replace("kappa=0.2", "kappa=0.2 kf=1e-24 af=1"),
    "bsim3_amp_noise": DECKS["bsim3_amp"],
}


def listvecs(txt):
    return [m.group(1) for m in re.finditer(r"^\s+(onoise[\w.]+)\s*:", txt, re.M) if m.group(1) != "onoise_spectrum"]


def rdn(fn, ncol):
    path = os.path.join(OUT, fn)
    if not os.path.exists(path):
        raise CouldNotRun(f"missing {fn}")
    rows = []
    for l in open(path).read().split("\n")[1:]:
        t = l.split()
        if len(t) >= 1 + ncol:
            rows.append([float(x) for x in t[:1 + ncol]])
    if not rows:
        raise CouldNotRun(f"empty {fn}")
    return rows


def check_pnoise_dev(name, bound=1e-5, pts=200000, floor=1e-6, env=None, tag=""):
    """total and every per-device contribution of `.pnoise dev=1` (td, generic path) vs plain `.noise` (summary table on)
    on an LTI-biased deck.  A contribution below `floor` of the total is compared in absolute terms (error / (floor*total) * floor)."""
    t0 = time.time()
    try:
        deck = NDECKS[name]
        nz = deck + OPTS + ".noise v(out) Vin %s 1\n" % SWEEP + ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nsetplot noise1\ndisplay\n.endc\n.end\n"
        txt, _ = run(nz, "nzl_" + name)
        vecs = listvecs(txt)
        if not vecs:
            raise CouldNotRun("no per-device vectors from .noise")
        wr = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nsetplot noise1\nwrdata nzl_%s.txt onoise_spectrum %s\nquit\n.endc\n.end\n" % (name, " ".join(vecs))
        run(deck + OPTS + ".noise v(out) Vin %s 1\n" % SWEEP + wr, "nzw_" + name)
        pn = (deck + OPTS + '.pnoise 1k 0 clkgnd 2000 1 8 1e-6 out Vin %s cyclo method="shoot" solver="td" dev=1 pacpoints=%d\n' % (SWEEP, pts)
              + ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata pnl_%s.txt onoise_spectrum %s\nquit\n.endc\n.end\n" % (name, " ".join(vecs)))
        otxt, _ = run(pn, "pnl_" + name, env=env)
        if "REFUSED" in otxt:
            raise CouldNotRun("REFUSED: " + [l for l in otxt.split("\n") if "REFUSED" in l][0][:200])
        a = rdn("nzl_%s.txt" % name, 1 + len(vecs))
        b = rdn("pnl_%s.txt" % name, 1 + len(vecs))
    except CouldNotRun as ex:
        record(f".pnoise dev=1 vs .noise: {name}{tag}", "COULD-NOT-RUN", str(ex))
        return None
    names = ["total"] + vecs
    worst, wn, parts = 0.0, "", []
    for j, nm in enumerate(names):
        ws = 0.0
        for i in range(len(a)):
            ref, tot = a[i][1 + j], a[i][1]
            e = abs(b[i][1 + j] - ref) / (ref if ref >= floor * tot else tot)
            ws = max(ws, e)
        parts.append(f"{nm} {ws:.1e}" + (f" [{100 * (a[0][1 + j] / a[0][1]) ** 2:.1f}% of power at 1k]" if j and a[0][1 + j] ** 2 > 0.02 * a[0][1] ** 2 else ""))
        if ws > worst:
            worst, wn = ws, nm
    ok = worst <= bound
    record(f".pnoise dev=1 vs plain .noise (LTI bias): {name}{tag}", "pass" if ok else "FAIL",
           f"worst rel {worst:.3e} at {wn} (bound {bound:g}); total at 1k/100k {a[0][1]:.4g}/{a[-1][1]:.4g} V/rtHz; " + ", ".join(parts) + f"; wall {time.time() - t0:.0f}s")
    return worst


# --- card #544 items 3 and 4: free-running orbit and ideal-switch events with model charges -------------------------------------
FIXD = os.path.join(HERE, "fixtures")


def run_fx(fixture, name, env=None, timeout=3000):
    """run tests/fixtures/<fixture> (decks write their outputs under ~/Work/data/ngspice_augmentation/544/test)."""
    os.makedirs(OUT, exist_ok=True)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    deck, cwd = os.path.join(FIXD, fixture), FIXD
    if OUT != OUT_FIXED:   # card #713: PAC_CHARGE_OUT -- a copy of the fixture whose output paths name OUT
        txt = open(deck).read()
        if "\n.include" in txt.lower() or "\n.lib" in txt.lower():
            raise CouldNotRun(f"{fixture}: includes files; cannot be redirected with PAC_CHARGE_OUT")
        deck = os.path.join(OUT, "fx_" + fixture)
        with open(deck, "w") as fh:
            fh.write(txt.replace(OUT_FIXED, OUT))
        cwd = OUT
    try:
        r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", NG, "-b",
                            deck], capture_output=True, text=True, cwd=cwd, env=e, timeout=timeout)
    except (FileNotFoundError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{fixture}: {ex}")
    return (r.stdout + r.stderr).replace("\r", "\n"), r.returncode


def rd_rows(fn):
    return rd(fn)


def dbdeg(z):
    return 20 * math.log10(abs(z)), math.degrees(math.atan2(z.imag, z.real))


# pac_perturbed.py --freerun (independent, nonlinear, Hann-windowed DFT of two long .tran runs; no shared code with the engine), run on
# tests/fixtures/r8_plain (= r8_pac_ctrl2out without analysis cards): amplitude 0.001 V at Vrefsrc / 0.01 V at Vin, stabtime 300 us,
# reltol 1e-4; "halving" = the amplitude-halving spread the tool reports.  ~/Work/data/ngspice_augmentation/544/r8/pp, pp2.
R8_REF = {   # case -> {f: (dB, deg, halving dB)}
    "ctrl2out": {20000.0: (7.8977, -0.2728, 0.009), 50000.0: (7.9969, -0.8614, 0.009), 80000.0: (8.1740, -1.4905, 0.008)},
    "line2out": {20000.0: (-21.0533, -0.1610, 0.003), 50000.0: (-20.9651, -0.4295, 0.028), 80000.0: (-20.7652, -0.7912, 0.007)},
}


def check_r8_freerun(case="ctrl2out", env=None, tag="", bound_db=0.05, bound_deg=0.2):
    """item 3: .pac solver=td on a FREE-RUNNING transistor-level (BSIM4) hysteretic buck (fixtures/r8_pac_<case>.cir: R7's power stage,
    drivers, dead time, comparator; hysteresis by resistive positive feedback; trigger=cx rise) through the generic model-charge path.
    Bound: 0.05 dB / 0.2 deg against the reference above (card #549: was 0.1 / 0.5; project target 0.5 dB / 3 deg, inside f <= 0.25 fsw = 99 kHz;
    reference uncertainty = halving spread <= 0.028 dB; measured engine - reference with the hidden-state free-running orbit: ctrl2out 0.0018 dB /
    0.050 deg, line2out 0.0179 dB / 0.054 deg; before: 0.0044 / 0.052 and 0.0205 / 0.051)."""
    fx = "r8_pac_%s.cir" % case
    name = f"R8 free-running BSIM4 hysteretic buck, .pac td {case} 20/50/80 kHz vs pac_perturbed --freerun{tag}"
    ref = R8_REF[case]
    if not ref:
        record(name, "COULD-NOT-RUN", "no reference recorded")
        return
    try:
        t0 = time.time()
        txt, rc = run_fx(fx, "r8_" + case, env=env)
        if "REFUSED" in txt or ".pac/.pstb skipped" in txt:
            raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l or ".pac/.pstb skipped" in l][:1] or txt[-300:])
        rows = rd_rows("r8_pac_%s.txt" % case)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    parts, worst_db, worst_deg = [], 0.0, 0.0
    for fq, z in rows:
        if fq not in ref:
            continue
        d, p = dbdeg(z)
        rd_, rp = ref[fq][0], ref[fq][1]
        worst_db = max(worst_db, abs(d - rd_)); worst_deg = max(worst_deg, abs(p - rp))
        parts.append(f"{fq / 1e3:g}k {d:.4f} dB/{p:.3f} deg vs {rd_:.4f}/{rp:.3f}")
    ok = len(parts) == len(ref) and worst_db <= bound_db and worst_deg <= bound_deg
    record(name, "pass" if ok else "FAIL", "; ".join(parts) + f"; worst {worst_db:.4f} dB / {worst_deg:.3f} deg (bound {bound_db} dB / {bound_deg} deg), wall {time.time() - t0:.0f}s")
    return worst_db


# tools/pac_perturbed.py (independent, nonlinear perturbed periodic steady state, .pss backend, amplitude halving <= 2e-4 dB, step
# refinement 800 -> 1600 points/cycle <= 2e-4 dB) on tests/fixtures/pac_r6o_ctrl2out.cir with vctrl_dc = 0.5285212 (the control DC of
# pac_r6ob_ctrl2out.cir), `--inject source:Vctrl --outputs vout`: (dB, deg) of vout.  ~/Work/data/ngspice_augmentation/544/r6om/ppd.
PP_R6O = {20000: (-7.8049, -169.9099), 50000: (-23.8535, -161.4504), 100000: (-34.9237, -147.3063)}


def check_salt_charge(env=None, tag="", bound_db=0.01, bound_deg=0.05, red=False):
    """item 4: ideal-switch (saltation) events with model charges on the orbit.  fixtures/pac_r6ob_ctrl2out.cir = R6o (ideal sw buck, 500 kHz) with a
    BSIM4 (capmod 2, W 100u) source follower and a 2 pF gate capacitor in the control path, so the control perturbation is low-passed by the gate
    network (Rgs 1 Meg, pole ~55 kHz: F = 0.81 / -14.8 deg at 20 kHz, 0.51 / -52.9 deg at 100 kHz) and model charges are on an orbit with events.
    The follower is LTI at its constant bias and decoupled from the converter, hence exactly H(R6ob) = F(f) * H(R6o at the same control DC) with F from
    `.ac`; the reference is F(f) * the independent pac_perturbed value of R6o above (an engine-free reference; the .pss backend of the tool cannot be
    used on R6ob itself: it restarts the hidden follower state every period, and its .tran backend does not settle).
    Measured: +0.042..0.043 dB (constant in frequency) before card #546, now 0.002 dB / 0.00 deg: the offset was the PSS period restart (the
    follower's source node and BSIM4 charges restarted at zero every period, which shifted the periodic gate state and so the follower's bias; the
    hidden-state continuation of pssshoot.c removes it).  Bounds 0.01 dB / 0.05 deg (were 0.1 / 0.5)."""
    name = "saltation + BSIM4 control-path charges, .pac td ctrl2out 20/50/100 kHz vs F(.ac) x pac_perturbed(R6o)" + (" [RED PROOF: PAC_NATIVE_HISTORY=1, must exceed 3 dB]" if red else tag)
    try:
        t0 = time.time()
        res, Fs = {}, {}
        b4 = open(os.path.join(FIXD, "pac_r6ob_ctrl2out.cir")).read().replace("\n.end\n", "\n")
        for f in ((20000,) if red else (20000, 50000, 100000)):
            d = b4 + (".pac 500e3 100e-6 vout 200 10 150 1e-6 lin 1 %d %d method=\"shoot\" solver=\"td\"\n" % (f, f) + wrline("r6ob_%d.txt" % f, "v(vout)"))
            txt, rc = run(d, "r6ob_%d" % f, env=env)
            if "REFUSED" in txt:
                raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l][0][:200])
            res[f] = rd("r6ob_%d.txt" % f)[0][1]
            run(b4 + ".ac lin 1 %d %d\n" % (f, f) + wrline("r6ob_ac_%d.txt" % f, "v(ctrlnode)"), "r6ob_ac_%d" % f)
            Fs[f] = rd("r6ob_ac_%d.txt" % f)[0][1]
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    parts, wdb, wdeg = [], 0.0, 0.0
    for f in res:
        h = 10 ** (PP_R6O[f][0] / 20) * cmath.exp(1j * math.radians(PP_R6O[f][1]))
        d, p = dbdeg(res[f])
        rdb, rdeg = dbdeg(Fs[f] * h)
        dp = (p - rdeg + 180.0) % 360.0 - 180.0
        wdb, wdeg = max(wdb, abs(d - rdb)), max(wdeg, abs(dp))
        parts.append(f"{f / 1e3:g}k {d:.3f} dB/{p:.2f} deg vs {rdb:.3f}/{rdeg:.2f}")
    ok = (wdb > 3.0) if red else (wdb <= bound_db and wdeg <= bound_deg)
    record(name, "pass" if ok else "FAIL", "; ".join(parts) + f"; worst {wdb:.3f} dB / {wdeg:.2f} deg" + ("" if red else f" (bound {bound_db} / {bound_deg})") + f", wall {time.time() - t0:.0f}s")
    return wdb


def check_salt_charge_identity(bound=2e-3, env=None, tag="", fixture="pac_r6ob_ctrl2out.cir", kind="BSIM4", measured="4.9e-4"):
    name = f"saltation + {kind} control path: H({fixture[4:8]}) == F(.ac) * H(R6o, same control DC), 100 Hz..100 kHz" + tag
    try:
        b4 = open(os.path.join(FIXD, fixture)).read().replace("\n.end\n", "\n")
        r6o = open(os.path.join(FIXD, "pac_r6o_ctrl2out.cir")).read().replace("\n.end\n", "\n")
        op = b4 + ".control\nop\nprint v(ctrlnode)\n.endc\n.end\n"
        txt, _ = run(op, "r6ob_op")
        m = re.search(r"v\(ctrlnode\) = ([-0-9.e+]+)", txt)
        if not m:
            raise CouldNotRun("no operating point")
        r6o = r6o.replace("vctrl_dc=0.55", "vctrl_dc=%s" % m.group(1))
        worst, parts = 0.0, []
        for f in (100, 5000, 100000):
            for nm, dk in (("r6ob", b4), ("r6o", r6o)):
                d = dk + ('.pac 500e3 100e-6 vout 200 10 150 1e-6 lin 1 %d %d method="shoot" solver="td"\n' % (f, f) + wrline("id_%s_%d.txt" % (nm, f), "v(vout)"))
                run(d, "id_%s_%d" % (nm, f), env=env if nm == "r6ob" else None)
            run(b4 + ".ac lin 1 %d %d\n" % (f, f) + wrline("id_ac_%d.txt" % f, "v(ctrlnode)"), "id_ac_%d" % f)
            hb, ho, F = rd("id_r6ob_%d.txt" % f)[0][1], rd("id_r6o_%d.txt" % f)[0][1], rd("id_ac_%d.txt" % f)[0][1]
            e = abs(hb - F * ho) / abs(F * ho)
            worst = max(worst, e)
            parts.append(f"{f}Hz {e:.2e}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    record(name, "pass" if worst <= bound else "FAIL", f"worst rel {worst:.2e} (bound {bound:g}; measured {measured}; was 5.4e-3 = 0.05 dB before card #546's hidden-state continuation of the period restart); " + ", ".join(parts))
    return worst


# --- card #544 item 1: time-varying bias, absolute noise level against a quasi-static brute force ---------------------------------
# LO sine at 100 kHz on the gate of a common-source stage whose time constant (Rd*C ~ 1-3 ns) is 3-4 orders below the LO period and whose
# device stays in saturation while it conducts (smooth gm(t)).  Then the output PSD at 1 kHz is the time average of the instantaneous
# PSD: <S_inst(t)> over one LO period, S_inst from plain `.noise` at the DC bias of N equally spaced LO phases (error O(tau/T) ~ 1e-4).
# (A first version of this check with Rd = 100k, a conduction sliver and 400 phases: the sliver is narrower than the phase grid
# and the FD engine's 20 harmonics, so neither the quasi-static nor the FD value is a reference there; the check below avoids it.)
QS_DECKS = {
    "mos1": ("* qs mos1\n.model nm nmos level=1 vto=0.7 kp=2e-4 lambda=0.02 tox=2e-8 cgso=5e-10 cgdo=5e-10 cgbo=2e-10 cj=1e-4 cjsw=1e-10 pb=0.8\n"
             "Vdd vdd 0 3\n%s\nM1 out g 0 0 nm w=10u l=2u as=1e-10 ad=1e-10 ps=4e-5 pd=4e-5\nRd vdd out 2k\nCl out 0 0.5p\nVin in 0 DC 0 AC 1\nRin in 0 1k\n", 0.9, 0.6),
    "diode": ("* qs diode\n.model dm d is=1e-14 cjo=0.2p vj=0.7 m=0.5 tt=0.1n rs=10 kf=1e-14 af=1\n%s\nRs g a 1k\nD1 a out dm\nRl out 0 1k\nCl out 0 1p\nVin in 0 DC 0 AC 1\nRin in 0 1k\n", 0.9, 0.5),
    "bjt": ("* qs bjt\n.model qm npn is=1e-15 bf=100 vaf=80 cje=0.2p cjc=0.1p tf=0.1n tr=1n rb=100 re=2 kf=1e-14 af=1\nVcc vcc 0 3\n%s\nRb g b 10k\nRc vcc out 2k\nQ1 out b 0 0 qm\nCl out 0 0.1p\nVin in 0 DC 0 AC 1\nRin in 0 1k\n", 0.75, 0.15),
    "bsim4": ("* qs bsim4\n.model " + BSIM % (2, 0, "") + "\nVdd vdd 0 3\n%s\nM1 out g 0 0 nch4 w=1u l=0.5u as=1p ad=1p ps=4u pd=4u\nRd vdd out 5k\nCl out 0 50f\nVin in 0 DC 0 AC 1\nRin in 0 1k\n", 0.8, 0.5),
}


def check_pnoise_qs(kind, N=96, bound=5e-4, env=None, tag=""):
    name = f".pnoise td vs quasi-static time average of .noise (LO 100 kHz, 1 kHz): {kind}{tag}"
    if kind == "bjt":
        bound = max(bound, 1e-3)   # card #556: its tau/T quasi-static error (3e-4 white, 4.7e-4 flicker part) is the slowest of the four decks
    t0 = time.time()
    try:
        deck, vdc, vac = QS_DECKS[kind]
        lo_tran = "Vlo g 0 SIN(%g %g 100k)" % (vdc, vac)
        lo_dc = "Vlo g 0 DC %s"
        opts = ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n"
        txt, _ = run(deck % (lo_dc % "0.9") + opts + ".noise v(out) Vin dec 1 1000 1000 1\n.control\nrun\nsetplot noise1\ndisplay\n.endc\n.end\n", "qs_list_" + kind)
        vecs = listvecs(txt)
        wr = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=12\nrun\nsetplot noise1\nwrdata qs_%s.txt onoise_spectrum %s\nquit\n.endc\n.end\n" % (kind, " ".join(vecs))
        acc = None
        phases = []
        for k in range(N):
            vg = vdc + vac * math.sin(2 * math.pi * (k + 0.5) / N)
            run(deck % (lo_dc % ("%.9g" % vg)) + opts + ".noise v(out) Vin dec 1 1000 1000 1\n" + wr, "qs_%s" % kind)
            row = rdn("qs_%s.txt" % kind, 1 + len(vecs))[-1]
            phases.append(row[1:])
            acc = [x * x for x in row[1:]] if acc is None else [a + x * x for a, x in zip(acc, row[1:])]
        qs = [math.sqrt(a / N) for a in acc]
        # card #556: the quasi-static TIME AVERAGE is exact for white sources only.  A flicker source n = g(t) w(t) whose g(t) is
        # modulated at the 100 kHz LO puts its modulation harmonics q at S_w(f - q fs), not S_w(f):  S_out(f) = Sum_q |a_q|^2 f/|f - q fs|
        # (ef = 1 for these decks), a_q = DFT of the per-phase amplitude a(t) = onoise_1overf(phase) (positive: the same sign at every
        # phase).  The reference is therefore the time average for everything that is not flicker and this sum for the flicker part
        # (a composite -- device total and spectrum -- keeps its white part and swaps its flicker part).
        names0 = ["total"] + vecs
        fol = {}
        for j, nm in enumerate(names0):
            if "1overf" in nm and "1overfsw" not in nm:
                a = [ph[j] for ph in phases]
                tot = 0.0
                for q in range(-(N // 2) + 1, N // 2 + 1):
                    aq = sum(a[i] * complex(math.cos(2 * math.pi * q * i / N), -math.sin(2 * math.pi * q * i / N)) for i in range(N)) / N
                    rho = 1.0 if q == 0 else 1000.0 / abs(1000.0 - q * 100e3)
                    tot += abs(aq) ** 2 * rho
                fol[nm] = math.sqrt(tot)
        fl_of = {}
        for nm in fol:
            base = nm[:nm.rindex("1overf") - 1]
            fl_of.setdefault(base, []).append(nm)
        qs_old = list(qs)
        for j, nm in enumerate(names0):
            if nm == "total":
                parts = list(fol)
            else:
                parts = fl_of.get(nm, [])
            if parts:
                white = qs_old[j] ** 2 - sum(qs_old[names0.index(x)] ** 2 for x in parts)
                qs[j] = math.sqrt(max(white, 0.0) + sum(fol[x] ** 2 for x in parts))
            elif nm in fol:
                qs[j] = fol[nm]
        pn = (deck % lo_tran + opts + '.pnoise 100k 50u out 400 1 20 1e-6 out Vin lin 1 1000 1000 cyclo method="shoot" solver="td" dev=1 pacpoints=8000\n'
              + ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=12\nrun\nwrdata qsp_%s.txt onoise_spectrum %s\nquit\n.endc\n.end\n" % (kind, " ".join(vecs)))
        otxt, _ = run(pn, "qsp_" + kind, env=env)
        if "REFUSED" in otxt:
            raise CouldNotRun("REFUSED")
        td = rdn("qsp_%s.txt" % kind, 1 + len(vecs))[-1][1:]
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return None
    names = ["total"] + vecs
    worst, parts = 0.0, []
    for j, nm in enumerate(names):
        if qs[j] < 0.01 * qs[0] and j:
            continue
        e = abs(td[j] - qs[j]) / qs[j]
        worst = max(worst, e)
        parts.append(f"{nm} td {td[j]:.5e} qs {qs[j]:.5e} ({e:.1e})")
    record(name, "pass" if worst <= bound else "FAIL", f"worst rel {worst:.2e} (bound {bound:g}, N={N} phases); " + "; ".join(parts) + f"; wall {time.time() - t0:.0f}s")
    return worst


def check_refusals():
    # (a) non-affine load guard (test hook names the device type)
    try:
        txt, rc = run(pac_deck("bsim4_cs_capmod0", 20000), "refuse_nonaffine", env={"PAC_GEN_TEST_NONAFFINE": "BSIM4"})
        ok = ("REFUSED" in txt and "BSIM4" in txt and "not affine" in txt)
        record("guard: non-affine device load is refused, naming the type", "pass" if ok else "FAIL",
               [l for l in txt.split("\n") if "REFUSED" in l][:1] or txt[-300:])
    except CouldNotRun as ex:
        record("guard: non-affine device load is refused", "COULD-NOT-RUN", str(ex))
    # (a2) BSIM4 non-quasi-static: since card #559 the relaxation node is a state of the generic recursion (accepted, validated vs a transient FRA in
    # checks_559); the pre-#559 refusal is restored by PAC_NQS_REFUSE=1 (check_nqs_refusals)
    try:
        d = pac_deck("bsim4_cs_capmod2", 20000).replace("acnqsmod=0", "acnqsmod=0 trnqsmod=1")
        assert "trnqsmod=1" in d
        txt, rc = run(d, "refuse_nqs")
        record("BSIM4 trnqsmod=1 is accepted since card #559 (non-quasi-static note, no refusal)", "pass" if ("REFUSED" not in txt and "non-quasi-static MOS model" in txt) else "FAIL",
               [l for l in txt.split("\n") if "REFUSED" in l or "non-quasi" in l][:1] or txt[-300:])
    except (CouldNotRun, AssertionError) as ex:
        record("BSIM4 trnqsmod=1 accepted", "COULD-NOT-RUN", str(ex))
    # (b) saltation + charged diode is refused, saltation + chargeless diode still runs (native path)
    base = ("* sw + diode\nVin vin 0 DC 12\nVramp vramp 0 PULSE(0 1 0 2u 1n 1n 2u)\nVctrl vctrl 0 DC 0.5 AC 1\n"
            ".model SWMOD SW(vt=0 vh=0 ron=0.05 roff=1e6)\nS1 vin sw vctrl vramp SWMOD\nD1 0 sw DMOD\n.model DMOD D(IS=1e-14 RS=0.02 N=1.05 %s)\n"
            "L1 sw vout 10u\nResr vout voutc 0.01\nCout voutc 0 100u IC=6\nRload vout 0 5\n"
            '.pac 500k 0 vout 200 1 25 1e-6 lin 2 1000 248744 method="shoot" solver="td"\n.control\nrun\nwrdata sw.txt vout\n.endc\n.end\n')
    try:
        txt1, _ = run(base % "CJO=50p TT=1n", "sw_charged")
        txt0, _ = run(base % "", "sw_nocharge")
        # card #575: this deck (charge-only switch node, no snubber) is no longer refused: the routed jump with its continuation
        # carries it (accuracy vs the perturbed PSS: checks_575); with a snubber too since card #713 (checks_713)
        ok1 = "REFUSED" not in txt1 and "generic charge-based" in txt1 and "saltation via inductor-routed jump" in txt1
        # (card #583: a switch node with a diode uses the generic recursion, with or without charge)
        ok0 = "REFUSED" not in txt0 and "PAC-TD: time-domain" in txt0 and "PAC-TD: done" in txt0
        record("saltation + device-model charge at an inductor node without a snubber runs on the fast-node path since card #575", "pass" if ok1 else "FAIL",
               [l for l in txt1.split("\n") if "REFUSED" in l or "inductor-routed" in l][:1] or txt1[-300:])
        record("saltation + chargeless diode runs without refusal (generic recursion since card #583)", "pass" if ok0 else "FAIL",
               "ran" if ok0 else txt0[-300:])
    except CouldNotRun as ex:
        record("guard: saltation + charged diode", "COULD-NOT-RUN", str(ex))


def check_freerun_refusal():
    """card #544: a free-running (trigger=) orbit with device-model charge (R5 + a charged diode behind 1 Meg) is no longer refused: it runs on
    the generic path.  The diode is far from the loop, so the control-to-output value must equal the native R5 deck's (cjo=0 -> native path,
    bit-identical to s96) within 0.05 dB / 0.5 deg (measured 0.014 dB / 0.03 deg: the charged deck uses the interpolated crossing and the
    1e-6 polish profile)."""
    fx = os.path.join(HERE, "fixtures", "pac_r5_ctrl2out_20k.cir")
    out = os.path.join(OUT, "r5_diode_%s.txt")
    try:
        base = open(fx).read()
        assert "wrdata ~/Work/data" in base
        base = re.sub(r"wrdata \S+ ", "wrdata " + out + " ", base)
        dd = base.replace("\n.end\n", "\nDx 0 xd DX\nRx xd vout 1meg\n.model DX D(cjo=1p)\n.end\n").replace("%s", "d")
        d0 = base.replace("\n.end\n", "\nDx 0 xd DX\nRx xd vout 1meg\n.model DX D(cjo=0)\n.end\n").replace("%s", "n")
        assert "cjo=1p" in dd and "cjo=0" in d0
        out1, _ = run(dd, "freerun_charged", env={"PSS_FREERUN_PAC_EXPERIMENTAL": "1"})
        out0, _ = run(d0, "freerun_native", env={"PSS_FREERUN_PAC_EXPERIMENTAL": "1"})
        def last(fn):   # the R5 fixture's wrdata has no header line
            t = open(os.path.join(OUT, fn)).read().strip().split("\n")[-1].split()
            return complex(float(t[1]), float(t[2]))
        zc, zn = last("r5_diode_d.txt"), last("r5_diode_n.txt")
    except (CouldNotRun, AssertionError) as ex:
        record("free-running + device-model charge runs on the generic path", "COULD-NOT-RUN", str(ex))
        return
    ok_g = "generic charge-based" in out1 and "REFUSED" not in out1 and "PAC-TD: done" in out1
    ok_n = "generic charge-based" not in out0 and "PAC-TD: done" in out0
    (dc, pc), (dn, pn) = dbdeg(zc), dbdeg(zn)
    ok = ok_g and ok_n and abs(dc - dn) <= 0.05 and abs(pc - pn) <= 0.5
    record("free-running (trigger=) R5 + charged diode: generic path, equals the native R5 value (0.05 dB / 0.5 deg)", "pass" if ok else "FAIL",
           f"charged {dc:.4f} dB/{pc:.3f} deg (generic={ok_g}), native {dn:.4f} dB/{pn:.3f} deg (native path={ok_n})")


def check_force_generic():
    fx = os.path.join(HERE, "fixtures", "pac_lti.cir")
    body = [l for l in open(fx).read().split("\n") if not l.startswith(".end")]
    deck = "\n".join(body) + '\n.pac 100k 0 clkgnd 2000 1 5 0.001 dec 5 1k 1meg method="shoot" solver="td" pacpoints=200000\n' + wrline("fg.txt", "out")
    try:
        run(deck, "fg_native")
        a = rd("fg.txt")
        txt, _ = run(deck, "fg_generic", env={"PAC_FORCE_GENERIC": "1"})
        b = rd("fg.txt")
        w = worst_rel(a, b)
        assert "generic charge-based" in txt
    except (CouldNotRun, AssertionError) as ex:
        record("native path == generic path on a native-only deck", "COULD-NOT-RUN", str(ex))
        return
    record("native C/L path == generic path on a native-only deck (pac_lti.cir, PAC_FORCE_GENERIC=1)", "pass" if w <= 1e-6 else "FAIL",
           f"worst rel {w:.3e} (bound 1e-6)")


def check_klu():
    try:
        d = pac_deck("mos1_amp", 100000)
        run(d, "klu_off")
        a = rd("td_mos1_amp.txt")
        run(d.replace(OPTS, OPTS + ".options klu\n"), "klu_on")
        b = rd("td_mos1_amp.txt")
        w = worst_rel(a, b)
    except CouldNotRun as ex:
        record("klu == sparse (generic path)", "COULD-NOT-RUN", str(ex))
        return
    record("`.options klu` == sparse on a MOS1 deck (generic path)", "pass" if w <= 1e-9 else "FAIL", f"worst rel {w:.3e} (bound 1e-9)")


def check_asan():
    if not ASAN:
        record("ASan generic .pac/.pxf adjoint/.pnoise", "COULD-NOT-RUN", "NGSPICE_ASAN not set")
        return
    try:
        r = load_repro()
        deck = r.deck(2, r.MODELS["bsim4"], analysis="")
        outs = []
        for card in ('.pac 500e3 20e-6 out 400 10 50 1e-7 lin 3 1000 3000 method="shoot" solver="td" pacpoints=400',
                     '.pxf 500e3 20e-6 out 400 10 50 1e-7 out lin 2 1000 3000 sources="Vc" method="shoot" solver="td" adjoint=1 pacpoints=400',
                     '.pnoise 500e3 20e-6 out 400 10 50 1e-7 out Vc lin 2 1000 3000 cyclo method="shoot" solver="td" pacpoints=400'):
            d = deck + card + "\n.options reltol=1e-4\n.control\nrun\n.endc\n.end\n"
            txt, rc = run(d, "asan_run", env={"ASAN_OPTIONS": "detect_leaks=0"}, exe=ASAN, timeout=3000)
            outs.append((card.split()[0], rc, "AddressSanitizer" in txt, "generic charge-based" in txt))
        bad = [o for o in outs if o[2] or not o[3]]
        record("ASan: generic .pac / .pxf adjoint=1 / .pnoise on a BSIM4 chain", "pass" if not bad else "FAIL", str(outs))
    except Exception as ex:
        record("ASan generic paths", "COULD-NOT-RUN", f"{type(ex).__name__}: {ex}")


def red_proofs():
    global QUIET
    red = {"PAC_NATIVE_HISTORY": "1"}
    QUIET = True
    w = check_lti("mos1_amp", bound=1e-6, pts=100000, env=red, tag=" [RED PROOF: PAC_NATIVE_HISTORY=1, must FAIL]")
    last = results[-1] if results else None
    if last and last[0].startswith("LTI vs .ac: mos1_amp [RED"):
        results.pop()
        QUIET = False
        ok = (w is not None and w > 1e-2)
        record("red proof: PAC_NATIVE_HISTORY=1 fails the LTI check (old behaviour reproduced)", "pass" if ok else "FAIL",
               f"worst rel {w:.3e} (must be > 1e-2); the same deck passes at 3.6e-9 without the switch" if w is not None else last[2])
    QUIET = True
    w2 = check_chain(models=("level1_cg1e-8",), stages=(2,), bound=0.01, env=red, tag=" [RED PROOF]")
    last = results[-1]
    QUIET = False
    if last[0].startswith("bug1 chain") and "[RED PROOF]" in last[0]:
        results.pop()
        ok = (w2 is not None and w2 > 0.2)
        record("red proof: PAC_NATIVE_HISTORY=1 fails the inverter-chain check (bug 1 reproduced)", "pass" if ok else "FAIL",
               f"error {100 * w2:.1f}% (must be > 20%)" if w2 is not None else last[2])


def red_proofs_556():
    """card #556 red proof: the pre-#556 flicker evaluation (density at the output frequency for every sideband) must fail the folded
    quasi-static check of the diode deck."""
    global QUIET
    QUIET = True
    w = check_pnoise_qs("diode", env={"PNOISE_FLICKER_OLD": "1"}, tag=" [RED PROOF #556]")
    r = results.pop() if results and "RED PROOF" in results[-1][0] else None
    QUIET = False
    record("red proof: PNOISE_FLICKER_OLD=1 (flicker density at the output frequency) fails the folded-flicker quasi-static check",
           "pass" if (w is not None and w > 1e-2) else "FAIL",
           f"diode mismatch {w:.3e} (must be > 1e-2)" if w is not None else (r[2] if r else "no result"))


def red_proofs_544():
    """card #544 red proofs: each new check must be able to fail."""
    global QUIET
    QUIET = True
    w = check_pnoise_qs("diode", env={"PAC_NOISE_STALE_STATE": "1"}, tag=" [RED PROOF]")
    r = results.pop() if results and "RED PROOF" in results[-1][0] else None
    QUIET = False
    record("red proof: PAC_NOISE_STALE_STATE=1 (pre-#544 noise operating point) fails the time-varying-bias noise check",
           "pass" if (w is not None and w > 1e-2) else "FAIL", f"diode quasi-static mismatch {w:.3e} (must be > 1e-2); the same deck passes at 1.6e-6" if w is not None else (r[2] if r else "no result"))
    QUIET = True
    w = check_r8_freerun("ctrl2out", env={"PAC_NATIVE_HISTORY": "1"}, tag=" [RED PROOF]")
    r = results.pop() if results and "RED PROOF" in results[-1][0] else None
    QUIET = False
    record("red proof: PAC_NATIVE_HISTORY=1 fails the free-running R8 check",
           "pass" if (w is not None and w > 1.0) else "FAIL", f"worst {w:.2f} dB (must be > 1 dB)" if w is not None else (r[2] if r else "no result"))
    check_salt_charge(env={"PAC_NATIVE_HISTORY": "1"}, red=True)
    # card #546: the hidden-state continuation of the PSS period restart is what makes the follower identity exact
    QUIET = True
    w = check_salt_charge_identity(env={"PSS_NO_HIDDEN": "1"}, tag=" [RED PROOF]")
    r = results.pop() if results and "RED PROOF" in results[-1][0] else None
    QUIET = False
    record("red proof: PSS_NO_HIDDEN=1 (pre-#546 period restart) fails the follower identity check",
           "pass" if (w is not None and w > 3e-3) else "FAIL",
           f"worst rel {w:.2e} (must be > 3e-3; the same deck passes at 4.9e-4)" if w is not None else (r[2] if r else "no result"))


# ---------------------------------------------------------------------------------------------------------------------------------
# card #549: `.options method=gear` for the td recording pass (it failed with "error 1" since #465: PSSpacRecordStep refused every
# non-trapezoidal step).  The generic charge recursion carries a BDF/Gear step of any order (J(ag0) y_k = -sum_j ag_j Q_{k-j}, the charge
# history Q_{k-1..k-6} instead of the trapezoidal derivative history), so gear is routed through it for native decks too.
OUT549 = os.path.join(OUT, "549") if OUT != OUT_FIXED else os.path.expanduser("~/Work/data/ngspice_augmentation/549/test")


def _gear_run(tag, meth, fixture, card, sig, env=None, outdir=None):
    """fixture body without its own analysis card / control block + `.options method=<meth>` + <card> + wrdata <sig>."""
    global OUT
    body = open(os.path.join(FIXD, fixture)).read().rsplit(".end", 1)[0]
    body = re.sub(r"(?ms)^\.control.*?^\.endc\s*", "", body)
    body = re.sub(r"(?m)^\.(pac|pstb|pss|pnoise|pxf)\b.*\n", "", body)
    deck = (body + f".options method={meth}\n" + card + ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n"
            f"wrdata g_{tag}.txt {sig}\nquit\n.endc\n.end\n")
    save = OUT
    OUT = outdir or OUT549
    try:
        txt, rc = run(deck, "g_" + tag, env=env)
        rows = None
        if os.path.exists(os.path.join(OUT, f"g_{tag}.txt")) and not re.search(r"FAILED|REFUSED", txt):
            rows = rd(f"g_{tag}.txt")
        return txt, rows
    finally:
        OUT = save


def _gear_diff(a, b):
    """worst (dB, deg) difference of two complex sweeps with equal frequencies."""
    if len(a) != len(b) or any(abs(x[0] - y[0]) > 1e-6 * x[0] for x, y in zip(a, b)):
        raise CouldNotRun("gear/trap frequency lists differ")
    wd = wp = 0.0
    for x, y in zip(a, b):
        if not (math.isfinite(x[1].real) and math.isfinite(x[1].imag) and math.isfinite(y[1].real) and math.isfinite(y[1].imag)):
            return float("inf"), float("inf")   # a NaN curve is never an agreement (the red proof ends here: order-1 history diverges)
        d1, p1 = dbdeg(x[1]); d2, p2 = dbdeg(y[1])
        dp = abs(p1 - p2)
        wd = max(wd, abs(d1 - d2)); wp = max(wp, min(dp, 360 - dp))
    return wd, wp


def check_gear(red=False):
    """method=gear vs method=trap on R3 control-to-output (smooth buck, native C/L: the legacy trap path is replaced by the generic
    recursion for gear) at 20..180 kHz, five points, and on the R7 `.pstb` (BSIM4, 10 points 20..200 kHz).  The two integrators discretise
    the base trajectory differently, so they agree only to the tolerance of the recording pass: with `pacreltol=1e-7` R3 agrees to 0.0012 dB /
    0.004 deg (at the default 1e-4 deck tolerance the spread is 0.014 dB, trap itself moving 0.04 dB between 1e-4 and 1e-7), bounds 0.01 dB / 0.05 deg.
    Red proof PAC_GEAR_TEST_ORD1=1 (every gear step treated as order 1: history beyond Q_{k-1} dropped) must break the agreement."""
    tag = "red" if red else ""
    name = "method=gear vs method=trap, .pac td R3 ctrl2out 20..180 kHz, pacreltol=1e-7" + (" [RED PROOF: PAC_GEAR_TEST_ORD1=1, must exceed the bound]" if red else "")
    card = '.pac 500e3 300e-6 vout 200 10 150 1e-6 lin 5 20k 180k method="shoot" solver="td" pacreltol=1e-7\n'
    try:
        _, rt = _gear_run("r3_trap", "trap", "pac_r3_ctrl2out.cir", card, "v(vout)")
        txt, rg = _gear_run("r3_gear" + tag, "gear", "pac_r3_ctrl2out.cir", card, "v(vout)", env={"PAC_GEAR_TEST_ORD1": "1"} if red else None)
        if rt is None or rg is None:
            raise CouldNotRun("gear/trap run produced no sweep: " + "; ".join(l for l in txt.split("\n") if "PAC-TD" in l)[:300])
        wd, wp = _gear_diff(rt, rg)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    ok = wd <= 0.01 and wp <= 0.05
    if red:
        record(name, "pass" if not ok else "FAIL", f"worst {wd:.4f} dB / {wp:.3f} deg (the bound 0.01 dB / 0.05 deg must be exceeded)")
    else:
        record(name, "pass" if ok else "FAIL", f"worst {wd:.5f} dB / {wp:.4f} deg (bound 0.01 dB / 0.05 deg)")


def check_gear_r7():
    name = "method=gear vs method=trap, R7 .pstb (BSIM4, generic path) 20..200 kHz, pacreltol=1e-7"
    card = '.pstb 500e3 300e-6 vout 400 10 150 1e-6 lin 10 20000 200000 probe="Vprobe" method="shoot" pacreltol=1e-7\n'
    try:
        _, rt = _gear_run("r7_trap", "trap", "r7_pstb_tight.cir", card, "t")
        txt, rg = _gear_run("r7_gear", "gear", "r7_pstb_tight.cir", card, "t")
        if rt is None or rg is None:
            raise CouldNotRun("R7 gear/trap run produced no sweep")
        wd, wp = _gear_diff(rt, rg)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    ok = wd <= GEAR_R7_DB and wp <= GEAR_R7_DEG
    record(name, "pass" if ok else "FAIL", f"worst {wd:.4f} dB / {wp:.3f} deg (bound {GEAR_R7_DB} dB / {GEAR_R7_DEG} deg)")


GEAR_R7_DB, GEAR_R7_DEG = 0.02, 0.05   # measured 0.0047 dB (constant over 20..200 kHz: the gear PSS orbit differs by ~0.05 % of Vout ripple) / 0.002 deg


def check_gear_refusals():
    """gear is not carried by the adjoint/noise recursion (the transposed sweep replays the trapezoidal derivative history) nor by the saltation
    jump: both must say so, not return a curve."""
    name = "method=gear refuses .pnoise and an ideal-switch (saltation) event with a clear message"
    try:
        body = open(os.path.join(FIXD, "pxf_lti.cir")).read().rsplit(".end", 1)[0]
        deck = (body + ".options method=gear\n.pnoise 100k 0 clkgnd 200 1 5 0.001 out Vin dec 3 1k 1meg method=\"shoot\" solver=\"td\"\n"
                ".control\nset wr_singlescale\nset wr_vecnames\nrun\nwrdata g_pn.txt onoise_spectrum\nquit\n.endc\n.end\n")
        global OUT
        save = OUT
        OUT = OUT549
        try:
            t1, _ = run(deck, "g_pn")
        finally:
            OUT = save
        t2, _ = _gear_run("salt", "gear", "pac_r6o_ctrl2out.cir", '.pac 500e3 300e-6 vout 200 10 150 1e-6 lin 1 20k 20k method="shoot" solver="td"\n', "v(vout)")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    m1 = "method=gear is not supported by the adjoint/noise recursion" in t1
    m2 = "ideal-switch (saltation) event with method=gear" in t2
    record(name, "pass" if (m1 and m2) else "FAIL", f".pnoise message {'yes' if m1 else 'NO'}, saltation message {'yes' if m2 else 'NO'}")


# ---------------------------------------------------------------------------------------------------------------------------------
# card #559 (2026-10-03): (1) ideal-switch (saltation) events at a switch node that carries device-model charge on top of a native
# capacitor; (2) non-quasi-static MOS channel charge (BSIM4 trnqsmod, BSIM3 nqsmod) in td PAC.
#
# (1) tests/fixtures/pac_r6od_ctrl2out.cir: open-loop asynchronous buck, ideal `sw` (500 kHz) + freewheel diode with cj (300 pF) and tt (2 ns),
#     1 nF at the switch node; pac_r6odm_ctrl2out.cir: the freewheel path is the body diode (is, cbd) of an off NMOS whose gate sits behind 1 kOhm
#     and moves with the node through cgdo (capacitive divider).  Reference: tools/pac_perturbed.py (perturbed periodic steady state, .pss backend,
#     amplitude halving <= 0.004 dB, step refinement 800 -> 1600 points/cycle <= 1e-4 dB), `--inject source:Vctrl --outputs vout --fs 500e3`,
#     s115 binary, ~/Work/data/ngspice_augmentation/559/pp/{od,odm}.  A converged transient FRA (tools/fra_curves_fork.run_transfer_point_conv, 3 ms settle,
#     amplitude 0.1 V, tmax 2 ns, two settings) agrees with the perturbed PSS to 0.05 dB at 20 kHz and 0.3 dB at 50 kHz (its own spread; 100 kHz is
#     noise-limited): 559/fra/res_odL.json.  The dB/deg (z) pairs are the reference values.
PP_OD = {20000: (-7.4671, -164.6137), 50000: (-23.4620, -159.6350), 100000: (-34.5263, -146.6731)}
PP_ODM = {20000: (-7.4858, -166.4806), 50000: (-23.5063, -160.3457), 100000: (-34.5735, -147.0254)}


def _salt_deck(fixture, f, repl=(), opts=".options reltol=1e-4\n", out="o.txt"):
    b = open(os.path.join(FIXD, fixture)).read().replace("\n.end\n", "\n")
    for a, c in repl:
        assert a in b, a
        b = b.replace(a, c)
    return b + opts + '.pac 500e3 100e-6 vout 200 10 150 1e-6 lin 1 %g %g method="shoot" solver="td"\n' % (f, f) + wrline(out, "v(vout)")


def check_salt_model_charge(fixture, ref, label, env=None, tag="", bound_db=0.1, bound_deg=0.5, red_min=None):
    """Engine td .pac ctrl2out at 20/50/100 kHz vs the perturbed-PSS reference.  red_min: red proof, the worst error must EXCEED it (20 kHz only)."""
    name = f"saltation + device-model charge at the switch node ({label}), .pac td ctrl2out vs perturbed PSS" + tag
    t0 = time.time()
    parts, wdb, wdeg = [], 0.0, 0.0
    try:
        for f in ((20000,) if red_min is not None else (20000, 50000, 100000)):
            txt, rc = run(_salt_deck(fixture, f, out="sm_%s_%d.txt" % (label[:4].replace(" ", "_"), f)), "sm_%s_%d" % (label[:4].replace(" ", "_"), f), env=env)
            if "REFUSED" in txt:
                raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l][0][:240])
            z = rd("sm_%s_%d.txt" % (label[:4].replace(" ", "_"), f))[0][1]
            d, p = dbdeg(z)
            rdb, rdeg = ref[f]
            dp = (p - rdeg + 180.0) % 360.0 - 180.0
            wdb, wdeg = max(wdb, abs(d - rdb)), max(wdeg, abs(dp))
            parts.append(f"{f / 1e3:g}k {d:.4f} dB/{p:.3f} deg vs {rdb:.4f}/{rdeg:.3f}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return None
    if red_min is not None:
        record(name, "pass" if wdb > red_min else "FAIL", "; ".join(parts) + f"; worst {wdb:.3f} dB (must exceed {red_min} dB)")
    else:
        record(name, "pass" if (wdb <= bound_db and wdeg <= bound_deg) else "FAIL",
               "; ".join(parts) + f"; worst {wdb:.4f} dB / {wdeg:.3f} deg (bound {bound_db} / {bound_deg}), wall {time.time() - t0:.0f}s")
    return wdb


# card #568: native decks (no model charge at the switch node).  Until s115 the order-1 settle window after an event lasted 30 x C/G of the node with
# G counted from resistors and switches only: with a (chargeless) freewheel diode at the node the window ran for 10..92 accepted steps and the
# answer was 5 dB off, silently.  The window now also ends after 2 accepted steps.  (a) R6od with cj = tt = 0 vs tools/pac_perturbed (same settings
# as PP_OD; ~/Work/data/ngspice_augmentation/559/pp/odn and 568/pp/odn); (b) R6o, whose window was 14..66 steps on a node that stays stiff, moves by
# 5e-5 relative: pinned to the s115 value within 2e-4 (15-digit wrdata values at 50 kHz, `.options reltol=1e-4`, s115 binary).
PP_ODN = {20000: (-7.4629, -164.7743), 50000: (-23.4604, -159.6840), 100000: (-34.5250, -146.6774)}
CHARGELESS = (("cjo=300p tt=2n", "cjo=0 tt=0"),)
S115_R6O_50K = complex(-6.083823069644149e-02, -2.041346736289501e-02)


def check_salt_native_window(env=None, red_min=None):
    name = "saltation, chargeless diode at the switch node with a 1 nF capacitor (R6od, cj = tt = 0; capacitor jump), .pac td ctrl2out vs perturbed PSS"
    if red_min is not None:
        name += " [RED PROOF: PSS_SALT_MAXSTEPS=-1 (the s115 time-only window), must exceed %g dB]" % red_min
    parts, wdb, wdeg = [], 0.0, 0.0
    try:
        for f in ((50000,) if red_min is not None else (20000, 50000, 100000)):
            txt, rc = run(_salt_deck("pac_r6od_ctrl2out.cir", f, CHARGELESS, out="nw_%d.txt" % f), "nw_%d" % f, env=env)
            if "REFUSED" in txt:
                raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l][0][:240])
            d, p = dbdeg(rd("nw_%d.txt" % f)[0][1])   # (generic recursion since card #583; the capacitor jump keeps its window)
            rdb, rdeg = PP_ODN[f]
            dp = (p - rdeg + 180.0) % 360.0 - 180.0
            wdb, wdeg = max(wdb, abs(d - rdb)), max(wdeg, abs(dp))
            parts.append(f"{f / 1e3:g}k {d:.4f} dB/{p:.3f} deg vs {rdb:.4f}/{rdeg:.3f}")
    except (CouldNotRun, OSError, IndexError) as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    if red_min is not None:
        record(name, "pass" if wdb > red_min else "FAIL", "; ".join(parts) + f"; worst {wdb:.3f} dB")
    else:
        record(name, "pass" if (wdb <= 0.1 and wdeg <= 0.5) else "FAIL", "; ".join(parts) + f"; worst {wdb:.4f} dB / {wdeg:.3f} deg (bound 0.1 / 0.5)")


def check_salt_native_r6o():
    name = "saltation, R6o (native capacitor, no diode) stays within 2e-4 of its s115 value at 50 kHz (card #568 moved it 5e-5)"
    try:
        txt, rc = run(_salt_deck("pac_r6o_ctrl2out.cir", 50000, out="bi_0.txt"), "bi_0")
        z = rd("bi_0.txt")[0][1]   # (generic recursion for every ideal-switch circuit since card #585)
    except (CouldNotRun, OSError, IndexError) as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    e = abs(z - S115_R6O_50K) / abs(S115_R6O_50K)
    record(name, "pass" if e <= 2e-4 else "FAIL", f"rel {e:.2e}")


def check_salt_refusals():
    """What is still refused after #559, and why (docs/pss-pac.md 'Saltation with model charges, continued')."""
    # (a) model charge only at the switch node together with an RC snubber: refused since card #575 (0.07-0.11 dB / 0.57 deg on 4 decks then);
    #     card #713: runs (fast node + the snubber jump at the end of the transition; accuracy in checks_713), PSS_SALT_NO713=1 refuses as before
    #     (red_proofs_713), PSS_SALT_NO713=1 PAC_SALT_FASTSNUB=1 runs the s143 fast node without any snubber jump
    try:
        for f_ in ("rf_a.txt", "rf_a2.txt"):
            if os.path.exists(os.path.join(OUT, f_)):
                os.remove(os.path.join(OUT, f_))
        txt, rc = run(_salt_deck("pac_r6od_ctrl2out.cir", 20000, (("csw=1n", "csw=0"),), out="rf_a.txt"), "rf_a")
        ok = "REFUSED" not in txt and "card #713" in txt and os.path.exists(os.path.join(OUT, "rf_a.txt"))
        txt2, rc2 = run(_salt_deck("pac_r6od_ctrl2out.cir", 20000, (("csw=1n", "csw=0"),), out="rf_a2.txt"), "rf_a2", env={"PSS_SALT_NO713": "1", "PAC_SALT_FASTSNUB": "1"})
        ok2 = "REFUSED" not in txt2 and "card #713" not in txt2 and os.path.exists(os.path.join(OUT, "rf_a2.txt"))
        record("saltation at a node with device-model charge, no native capacitor and an RC snubber runs since card #713 (note printed; "
               "PSS_SALT_NO713=1 PAC_SALT_FASTSNUB=1 runs it without the snubber jump)",
               "pass" if (ok and ok2) else "FAIL", [l for l in (txt + txt2).split("\n") if "REFUSED" in l or "card #713" in l][:1] or txt[-300:])
    except CouldNotRun as ex:
        record("saltation at a charge-only node with a snubber", "COULD-NOT-RUN", str(ex))
    # (b) mutual inductance on the inductor-routed branch (chargeless diode): refused as before
    base = ("* sw + mutual\nVin vin 0 DC 12\nVramp vramp 0 PULSE(0 1 0 2u 1n 1n 2u)\nVctrl vctrl 0 DC 0.5 AC 1\n"
            ".model SWMOD SW(vt=0 vh=0 ron=0.05 roff=1e6)\nS1 vin sw vctrl vramp SWMOD\nD1 0 sw DMOD\n.model DMOD D(IS=1e-14 RS=0.02 N=1.05)\n"
            "L1 sw vout 10u\nResr vout voutc 0.01\nCout voutc 0 100u IC=6\nRload vout 0 5\nL2 x 0 10u\nK1 L1 L2 0.5\nRx x 0 1\n"
            '.pac 500k 0 vout 200 1 25 1e-6 lin 1 20000 20000 method="shoot" solver="td"\n.control\nrun\nquit\n.endc\n.end\n')
    try:
        txt, rc = run(base, "rf_b")
        record("guard: saltation + mutual inductance on the routed branch stays refused",
               "pass" if ("REFUSED" in txt and "mutual inductance" in txt) else "FAIL", [l for l in txt.split("\n") if "REFUSED" in l][:1] or txt[-300:])
    except CouldNotRun as ex:
        record("guard: saltation + mutual inductance", "COULD-NOT-RUN", str(ex))


# --- (2) non-quasi-static channel ----------------------------------------------------------------------------------------------
# BSIM4 trnqsmod=1 / BSIM3 nqsmod=1: the charge-relaxation node (qNode) is an ordinary circuit node whose charge terms are proportional to ag0 (the
# affine guard passes), i.e. an ordinary state of the generic recursion.  LTI-biased common-source stage (L = 20 um: relaxation time ~0.1 us),
# a slow clock on an unrelated net, `.pac` at 1..30 MHz.  `.ac` is NOT the reference: its acnqsmod is a single-pole approximation that is 20-40 %
# away from the transient model at omega*tau ~ 1 (BSIM4 at 10 MHz: .ac -1.05+3.37j, transient FRA -1.64+2.34j, quasi-static -2.56+2.95j).  The reference
# is a converged transient FRA of the same model (tools/fra_curves_fork.run_transfer_point_conv, startup op, 1 mV, 256 points/cycle, reltol 1e-6, 20 cycles
# after 20 periods of settling; a second run at 0.5 mV, 512 ppc, 1e-7, 40/40 cycles differs by 3e-5..8e-5): the values below are the second run.
NQS_BSIM4 = ("nch4 nmos level=54 version=4.8.2 mobmod=0 capmod=2 igcmod=0 igbmod=0 rdsmod=0 rbodymod=0 rgatemod=0 acnqsmod=%d trnqsmod=%d "
             "toxe=3e-9 toxp=3e-9 toxm=3e-9 vth0=0.5 k1=0.4 k2=0 k3=0 w0=2.5e-6 u0=300 vsat=1.1e5 cgso=2e-10 cgdo=2e-10 cgbo=1e-11 "
             "cgsl=1e-10 cgdl=1e-10 cf=1e-10 cjs=1e-3 cjd=1e-3 cjswd=1e-10 cjswgd=1e-10 xpart=0")
NQS_B3 = ("nb3 nmos level=8 version=3.3.0 tox=5e-9 vth0=0.5 cgso=2e-10 cgdo=2e-10 cj=1e-3 capmod=3 nqsmod=%d")
NQS_CLK = "Vclk clk 0 PULSE(0 1 0 1n 1n 5u 10u)\nRclk clk clkgnd 1k\nCclk clkgnd 0 1p\n"
FRA_NQS4 = {1e6: -2.898523389085374 + 0.296953610800028j, 3e6: -2.7745136207151977 + 0.8722809488476134j,
            1e7: -1.6410731017874856 + 2.339914629245763j, 3e7: 1.4337340350807914 + 2.247982102656766j}
FRA_QS4 = {1e7: -2.563147248531186 + 2.954627026804471j}
FRA_NQS3 = {3e6: -3.0640949832671494 + 0.5610854937346226j, 1e7: -2.2582797389356863 + 1.1301732727880023j, 3e7: -1.3232109932669514 + 0.9691447920035382j}


def _nqs_deck(kind, ac, tr, f):
    if kind == "bsim4":
        mod, dev = ".model " + NQS_BSIM4 % (ac, tr), "nch4"
    else:
        mod, dev = ".model " + NQS_B3 % tr, "nb3"
    return ("* nqs cs\n" + mod + "\n" + NQS_CLK + "Vdd vdd 0 3\nVin in 0 DC 0.9 AC 1\nRs in g 100\nM1 out g 0 0 %s w=40u l=20u as=4p ad=4p ps=10u pd=10u\n"
            "Rd vdd out 20k\nCl out 0 50f\n.options reltol=1e-6\n" % dev
            + '.pac 100k 0 clkgnd 2000 1 8 1e-6 lin 1 %g %g method="shoot" solver="td" pacpoints=40000\n' % (f, f)
            + wrline("nq_%s_%d%d_%g.txt" % (kind, ac, tr, f), "v(out)"))


def check_nqs_lti(kind, ac, tr, ref, label, bound=5e-4, env=None, red_min=None):
    name = f"NQS channel, td .pac vs transient FRA: {label}" + (" [RED PROOF]" if red_min else "")
    t0 = time.time()
    worst, parts = 0.0, []
    try:
        for f, h in ref.items():
            txt, rc = run(_nqs_deck(kind, ac, tr, f), "nq_%s_%d%d_%g" % (kind, ac, tr, f), env=env)
            if "REFUSED" in txt:
                raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l][0][:160])
            z = rd("nq_%s_%d%d_%g.txt" % (kind, ac, tr, f))[0][1]
            e = abs(z - h) / abs(h)
            worst = max(worst, e)
            parts.append(f"{f / 1e6:g}M {e:.1e}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return None
    ok = (worst > red_min) if red_min else (worst <= bound)
    record(name, "pass" if ok else "FAIL", ", ".join(parts) + f"; worst rel {worst:.2e}" + (f" (must exceed {red_min})" if red_min else f" (bound {bound:g})") + f", wall {time.time() - t0:.0f}s")
    return worst


# switching deck: 1 MHz LO 0.35 V on the gate of an NQS BSIM4 stage (the channel swings from cut-off into strong conduction), small signal at Vin.
# Reference: converged transient FRA of the same deck (two settings, 2 mV/256 ppc/1e-6/20 cycles and 1 mV/512 ppc/1e-7/40 cycles; the 3 MHz point is
# the noisier one: 8.6027/162.22 vs 8.6192/162.69), tools/pac_perturbed.py's .pss backend gives 8.6127/162.76 and 8.0874/127.88 (0.011 dB / 0.03 dB
# and 0.3 / 0.5 deg away): it restarts the relaxation node's history every period.
FRA_LPTV = {3e6: (8.6110, 162.45), 1e7: (8.1214, 127.347)}   # mean of the two FRA settings (dB, deg)


def check_nqs_lptv(bound_db=0.05, bound_deg=0.5, env=None, nq=1, red_min=None):
    name = "NQS channel in an LO-switched stage (1 MHz LO), td .pac vs transient FRA" + (" [RED PROOF: trnqsmod=0 (quasi-static) deck, must exceed %g dB]" % red_min if red_min else "")
    try:
        parts, wdb, wdeg = [], 0.0, 0.0
        for f, (rdb, rdeg) in FRA_LPTV.items():
            d = open(os.path.join(FIXD, "pac_nqs_lptv.cir")).read().replace("\n.end\n", "\n").replace(".param nq=1", ".param nq=%d" % nq)
            d += ".options reltol=1e-6\n.pac 1meg 0 out 2000 1 8 1e-6 lin 1 %g %g method=\"shoot\" solver=\"td\" pacpoints=4000\n" % (f, f) + wrline("lp_%g.txt" % f, "v(out)")
            txt, rc = run(d, "lp_%g" % f, env=env)
            if "REFUSED" in txt:
                raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l][0][:160])
            z = rd("lp_%g.txt" % f)[0][1]
            dd, pp = dbdeg(z)
            dp = (pp - rdeg + 180.0) % 360.0 - 180.0
            wdb, wdeg = max(wdb, abs(dd - rdb)), max(wdeg, abs(dp))
            parts.append(f"{f / 1e6:g}M {dd:.4f} dB/{pp:.2f} deg vs {rdb:.4f}/{rdeg:.2f}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return None
    ok = (wdb > red_min) if red_min else (wdb <= bound_db and wdeg <= bound_deg)
    record(name, "pass" if ok else "FAIL", "; ".join(parts) + f"; worst {wdb:.4f} dB / {wdeg:.2f} deg" + (f" (must exceed {red_min} dB)" if red_min else f" (bound {bound_db} / {bound_deg})"))
    return wdb


def check_nqs_refusals():
    # PAC_NQS_REFUSE=1 restores the pre-#559 refusal; `.pnoise` accepts NQS models since #715 (PAC_NQS_PNOISE_REFUSE=1 refuses); acnqsmod alone is accepted (td linearises the transient model)
    try:
        txt, rc = run(_nqs_deck("bsim4", 1, 1, 1e7), "nq_ref1", env={"PAC_NQS_REFUSE": "1"})
        record("guard: PAC_NQS_REFUSE=1 restores the refusal of an NQS model", "pass" if ("REFUSED" in txt and "trnqsmod" in txt) else "FAIL",
               [l for l in txt.split("\n") if "REFUSED" in l][:1] or txt[-200:])
        d = _nqs_deck("bsim4", 1, 1, 1e7).split(".pac")[0] + '.pnoise 100k 0 clkgnd 2000 1 8 1e-6 out Vin lin 1 1e6 1e6 cyclo method="shoot" solver="td"\n.control\nrun\nquit\n.endc\n.end\n'
        # card #715: accepted (no noise source on the relaxation node; accuracy in tests/test_pnoise_nqs.py); PAC_NQS_PNOISE_REFUSE=1 refuses again
        txt, rc = run(d, "nq_ref2")
        record(".pnoise on an NQS model is accepted since card #715 (runs to done, no refusal)", "pass" if ("REFUSED" not in txt and "PNOISE-TD: done" in txt) else "FAIL",
               [l for l in txt.split("\n") if "REFUSED" in l or "PNOISE-TD: done" in l][:1] or txt[-200:])
        txt, rc = run(d, "nq_ref3", env={"PAC_NQS_PNOISE_REFUSE": "1"})
        record("guard: PAC_NQS_PNOISE_REFUSE=1 restores the .pnoise refusal of an NQS model", "pass" if ("REFUSED" in txt and "PNOISE" in txt and "PAC_NQS_PNOISE_REFUSE" in txt) else "FAIL",
               [l for l in txt.split("\n") if "REFUSED" in l][:1] or txt[-200:])
    except CouldNotRun as ex:
        record("guard: NQS refusals", "COULD-NOT-RUN", str(ex))


def check_nqs_instance():
    # card #730: an NQS flag on the DEVICE line is found like one on the .model line (up to s147 only the model was asked: a refused family
    # with trnqsmod/nqsmod on the device line ran unvalidated without a message).  PAC_NQS_INSTANCE_OFF=1 restores the model-only detection.
    # .pnoise is not checked here: it refuses BSIM4 4.7 / BSIM3 3.2 noise devices anyway (outside its device scope, card #509/#544).
    def inst(d, flag):
        return re.sub(r"(?m)^(M1 .*)$", r"\1 " + flag, d)
    try:
        d4 = _nqs_deck("bsim4", 0, 0, 1e7)
        d47 = inst(d4.replace("version=4.8.2", "version=4.7.0"), "trnqsmod=1")
        txt, rc = run(d47, "nq_inst47")
        record("BSIM4 4.7 with trnqsmod=1 on the device line is refused (card #730)", "pass" if ("REFUSED" in txt and "on instance m1" in txt) else "FAIL",
               [l for l in txt.split("\n") if "REFUSED" in l][:1] or txt[-200:])
        txt, rc = run(d47, "nq_inst47off", env={"PAC_NQS_INSTANCE_OFF": "1"})
        record("[RED PROOF] PAC_NQS_INSTANCE_OFF=1 (s147 model-only detection): the same deck runs without a refusal",
               "pass" if ("REFUSED" not in txt and "PAC-TD: done" in txt) else "FAIL", [l for l in txt.split("\n") if "REFUSED" in l or "PAC-TD: done" in l][:1])
        d3 = re.sub(r"(?m)^(M1 .*?) as=.*$", r"\1", _nqs_deck("bsim3", 0, 0, 1e7).replace("version=3.3.0", "version=3.2.4"))
        txt, rc = run(inst(d3, "nqsmod=1"), "nq_inst32")
        record("BSIM3 3.2.4 with nqsmod=1 on the device line is refused (card #730)", "pass" if ("REFUSED" in txt and "BSIM3v32" in txt and "on instance m1" in txt) else "FAIL",
               [l for l in txt.split("\n") if "REFUSED" in l][:1] or txt[-200:])
        txt, rc = run(inst(d4, "trnqsmod=1"), "nq_inst48")
        record("BSIM4 4.8 with trnqsmod=1 on the device line is accepted with the NQS note naming the instance", "pass" if ("REFUSED" not in txt and "trnqsmod=1 on instance m1" in txt and "PAC-TD: done" in txt) else "FAIL",
               [l for l in txt.split("\n") if "non-quasi-static" in l or "REFUSED" in l][:1] or txt[-200:])
    except CouldNotRun as ex:
        record("NQS flag on the device line", "COULD-NOT-RUN", str(ex))


def check_nqs_pxf():
    """`.pxf adjoint=1` == forward on an NQS deck (the transposed recursion uses the stored J_k, C_k including the relaxation-node rows)."""
    try:
        res = []
        for adj in (0, 1):
            d = (_nqs_deck("bsim4", 1, 1, 1e7).split(".pac")[0] + '.pxf 100k 0 clkgnd 2000 1 8 1e-6 out lin 1 1e7 3e7 sources="Vin" method="shoot" solver="td" adjoint=%d pacpoints=40000\n' % adj
                 + wrline("nqpx%d.txt" % adj, "pxf_vin"))
            txt, rc = run(d, "nqpx%d" % adj)
            if "REFUSED" in txt:
                raise CouldNotRun("REFUSED")
            res.append(rd("nqpx%d.txt" % adj))
        w = worst_rel(res[0], res[1])
    except CouldNotRun as ex:
        record(".pxf adjoint vs forward on an NQS deck", "COULD-NOT-RUN", str(ex))
        return
    record(".pxf adjoint=1 vs adjoint=0 on a BSIM4 trnqsmod=1 deck", "pass" if w <= 1e-8 else "FAIL", f"worst rel {w:.2e} (bound 1e-8)")


# ---------------------------------------------------------------------------------------------------------------------------------
# card #575 (2026-10-06): ideal-switch events at an inductor node WITHOUT a dominant native capacitor.
#   fast node  : the node's only (or dominant) charge is a device's (diode cj / tt, MOS body diode and overlap), no RC snubber: the jump goes
#                to the inductor and is continued until the node voltage has settled (PSSsaltCont);
#   chargeless : the settle window is the event step + exactly 2 accepted steps (the time rule gave none after a located event: 3.0x);
#                an RC snubber's capacitor gets its own jump;
#   charge jump (#559) only where the native capacitor dominates (>= 3x the smallest device capacitance on the orbit) or sits with a snubber.
# References: tools/pac_perturbed.py, same settings as PP_OD, s117 binary, one run per deck (~/Work/data/ngspice_augmentation/575/m/ref_*; the
# 18-deck table is in docs/pss-pac.md).  OD = pac_r6od_ctrl2out.cir, ODM = pac_r6odm_ctrl2out.cir, SWC = the pac_sw_ctrl_direct topology.
SWC_575 = ("* sw + diode, no capacitor at the switch node\nVin vin 0 DC 12\nVramp vramp 0 PULSE(0 1 0 2u 1n 1n 2u)\nVctrl vctrl 0 DC 0.5 AC 1\n"
           ".model SWMOD SW(vt=0 vh=0 ron=0.05 roff=1e6)\nS1 vin sw vctrl vramp SWMOD\nD1 0 sw DMOD\n.model DMOD D(IS=1e-14 RS=0.02 N=1.05 CJO=%s TT=%s)\n"
           "L1 sw vout 10u\nResr vout voutc 0.01\nCout voutc 0 100u IC=6\nRload vout 0 5\n")
NOSNUB = ("Rsnub sw nsnub 10", "Rsnub sw nsnub 1e9")
D575 = {   # name: (fixture or None (SWC), replacements or (cjo, tt), reference {f: (dB, deg)}, what)
    "swc_50p_1n": (None, ("50p", "1n"), {20000: (-1.1907, -168.9907), 50000: (-17.2338, -161.1212), 100000: (-28.302, -147.1468)}, "fast node, diode cj 50 pF tt 1 ns"),
    "swc_300p_5n": (None, ("300p", "5n"), {20000: (-1.1914, -168.8931), 50000: (-17.2339, -161.0941), 100000: (-28.3019, -147.1462)}, "fast node, diode cj 300 pF tt 5 ns"),
    "swc_1n_0": (None, ("1n", "0"), {20000: (-1.194, -168.6234), 50000: (-17.2349, -161.0227), 100000: (-28.3024, -147.148)}, "fast node, diode cj 1 nF"),
    "od_nosnub": ("pac_r6od_ctrl2out.cir", (("csw=1n", "csw=0"), NOSNUB), {20000: (-7.4068, -167.2059), 50000: (-23.4372, -160.4421), 100000: (-34.5054, -146.8161)}, "fast node, R6od without node capacitor and snubber"),
    "odm_nosnub": ("pac_r6odm_ctrl2out.cir", (("csw=1n", "csw=0"), NOSNUB), {20000: (-7.4392, -169.03), 50000: (-23.4833, -161.1359), 100000: (-34.5532, -147.1668)}, "fast node, MOS body diode as the freewheel device"),
    "od_csw10p_nosnub": ("pac_r6od_ctrl2out.cir", (("csw=1n", "csw=10p"), NOSNUB), {20000: (-7.4071, -167.1918), 50000: (-23.4373, -160.4379), 100000: (-34.5056, -146.8159)}, "fast node, 10 pF native beside cj 300 pF (not dominant)"),
    "od_csw1n_nosnub": ("pac_r6od_ctrl2out.cir", (NOSNUB,), {20000: (-7.4389, -165.7329), 50000: (-23.4512, -159.9983), 100000: (-34.5173, -146.7707)}, "charge jump, 1 nF native dominates, no snubber"),
    "od_csw10p": ("pac_r6od_ctrl2out.cir", (("csw=1n", "csw=10p"),), {20000: (-7.4197, -166.578), 50000: (-23.4426, -160.2455), 100000: (-34.5102, -146.7508)}, "charge jump, 10 pF native with the snubber"),
    "od_cl_nosnub": ("pac_r6od_ctrl2out.cir", (("csw=1n cjo=300p tt=2n", "csw=0 cjo=0 tt=0"), NOSNUB), {20000: (-7.4041, -167.3459), 50000: (-23.4359, -160.4843), 100000: (-34.5043, -146.8206)}, "chargeless diode, vh > 0, no capacitor (was 3.0x off)"),
    "od_cl_snub": ("pac_r6od_ctrl2out.cir", (("csw=1n cjo=300p tt=2n", "csw=0 cjo=0 tt=0"),), {20000: (-7.4118, -166.9616), 50000: (-23.4388, -160.3736), 100000: (-34.5069, -146.7763)}, "chargeless diode with an RC snubber (was 8.4 dB off)"),
}


def _deck_575(name, f, reltol="1e-4", out=None, points=200):
    fx, repl, ref, what = D575[name]
    if fx is None:
        b = SWC_575 % repl
    elif fx == "@":
        b = repl
    else:
        b = open(os.path.join(FIXD, fx)).read().replace("\n.end\n", "\n")
        for a, c in repl:
            assert a in b, a
            b = b.replace(a, c)
    return b + ".options reltol=%s\n" % reltol + '.pac 500e3 100e-6 vout %d 10 150 1e-6 lin 1 %g %g method="shoot" solver="td"\n' % (points, f, f) + wrline(out or ("f5_%s_%d.txt" % (name, f)), "v(vout)")


def check_575(name, reltol="1e-4", bound_db=0.1, bound_deg=0.5, env=None, red_min=None, tag="", freqs=(20000, 50000, 100000), expect=None, red_freq=20000, points=200, absent=None, red_deg=None):
    """td .pac ctrl2out vs the perturbed-PSS reference of that deck.  red_min: red proof (the worst dB error must exceed it, at red_freq only).
    expect: a substring that must be in the log (which treatment ran)."""
    fx, repl, ref, what = D575[name]
    label = f"saltation without a dominant node capacitor: {what} [{name}, reltol {reltol}" + ("" if points == 200 else f", {points} points") + "]" + tag
    parts, wdb, wdeg = [], 0.0, 0.0
    try:
        for f in ((red_freq,) if red_min is not None else freqs):
            o = "f5_%s_%d.txt" % (name, f)
            txt, rc = run(_deck_575(name, f, reltol, o, points), "f5_%s_%d" % (name, f), env=env)
            if "REFUSED" in txt:
                raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l][0][:200])
            if expect and expect not in txt:
                record(label, "FAIL", f"expected '{expect}' in the log")
                return None
            if absent and absent in txt:
                record(label, "FAIL", f"'{absent}' must not be in the log")
                return None
            d, p_ = dbdeg(rd(o)[0][1])
            rdb, rdeg = ref[f]
            dp = (p_ - rdeg + 180.0) % 360.0 - 180.0
            wdb, wdeg = max(wdb, abs(d - rdb)), max(wdeg, abs(dp))
            parts.append(f"{f / 1e3:g}k {d - rdb:+.3f} dB {dp:+.2f} deg")
    except (CouldNotRun, OSError, IndexError) as ex:
        record(label, "COULD-NOT-RUN", str(ex))
        return None
    if red_min is not None and red_deg is not None:
        record(label, "pass" if wdeg > red_deg else "FAIL", "; ".join(parts) + f"; worst {wdeg:.2f} deg (must exceed {red_deg})")
    elif red_min is not None:
        record(label, "pass" if wdb > red_min else "FAIL", "; ".join(parts) + f"; worst {wdb:.3f} dB (must exceed {red_min})")
    else:
        record(label, "pass" if (wdb <= bound_db and wdeg <= bound_deg) else "FAIL", "; ".join(parts) + f"; worst {wdb:.3f} dB / {wdeg:.2f} deg (bound {bound_db} / {bound_deg})")
    return wdb


def checks_575():
    for n in ("swc_50p_1n", "swc_300p_5n", "swc_1n_0", "od_nosnub", "odm_nosnub", "od_csw10p_nosnub"):
        check_575(n, expect="saltation via inductor-routed jump")
    check_575("swc_50p_1n", reltol="1e-6")          # the charge jump was 3.6 dB off here at this tolerance
    check_575("od_csw10p_nosnub", reltol="1e-6")
    check_575("od_csw1n_nosnub", expect="formed on the charge vector")
    check_575("od_csw10p", expect="formed on the charge vector")
    check_575("od_cl_nosnub", bound_db=0.02, bound_deg=0.1)
    # chargeless + snubber: limited by the trajectory tolerance over the snubber's 10 ns (0.135 dB at 1e-4, 0.043 at 1e-6, 0.001 at 1e-7)
    check_575("od_cl_snub", bound_db=0.2, bound_deg=0.2)
    check_575("od_cl_snub", reltol="1e-6", bound_db=0.1, bound_deg=0.2)


def red_proofs_575():
    check_575("swc_50p_1n", env={"PSS_SALT_NOCONT": "1"}, red_min=10.0, tag=" [RED PROOF: PSS_SALT_NOCONT=1 (the single-step routed jump)]")
    check_575("swc_50p_1n", env={"PAC_SALT_QONLY": "1"}, red_min=0.2, tag=" [RED PROOF: PAC_SALT_QONLY=1 (the charge jump of #559 at a charge-only node)]")
    check_575("od_csw10p_nosnub", env={"PAC_SALT_QONLY": "1"}, red_min=0.2, tag=" [RED PROOF: PAC_SALT_QONLY=1 (charge jump with a non-dominant native capacitor: ran unrefused on s117)]")
    check_575("odm_nosnub", env={"PSS_SALT_NOCONT": "1"}, red_min=10.0, tag=" [RED PROOF: PSS_SALT_NOCONT=1]")
    check_575("od_cl_nosnub", env={"PSS_SALT_MAXSTEPS": "-1", "PAC_SALT_SMOOTHT": "0"}, red_min=5.0, tag=" [RED PROOF: PSS_SALT_MAXSTEPS=-1 + PAC_SALT_SMOOTHT=0 (jump before T with the time-only settle window of s117: 3.0x)]")
    check_575("od_cl_snub", env={"PSS_SALT_NOSNUB": "1"}, red_min=0.5, tag=" [RED PROOF: PSS_SALT_NOSNUB=1 (no jump on the snubber capacitor)]")


# ---------------------------------------------------------------------------------------------------------------------------------
# card #579 (2026-10-06): the same matrix method on other circuits.  Boost (low-side ideal switch, output diode, 22 uF + 10 mohm ESR), a
# synchronous buck with two ideal switches and dead time (body diode), line-to-output and a second duty cycle on the R6od buck.  On s118 the
# boost was 71 dB off (the result ~0) and the charged synchronous buck 0.25 dB, unrefused.  References: tools/pac_perturbed.py, s118 binary,
# ~/Work/data/ngspice_augmentation/579/m/ref_*.
BOOST_579 = ("* boost, low-side ideal switch, output diode\n.param fs=500k tclk={1/fs}\nVin vin 0 DC 5\nRdcr vin vl 10m\nL1 vl sw 10u IC=0\n"
             "S1 sw 0 ctrlnode ramp SWMOD\n.model SWMOD SW(vt=0 vh=0.005 ron=10m roff=1Meg)\nD1 sw vout DMOD\n"
             ".model DMOD D(IS=1e-9 RS=0.05 N=1.05 CJO=%s VJ=0.7 M=0.5 TT=%s)\nCout vout cx 22u IC=0\nResr cx 0 10m\nRload vout 0 10\n"
             "Vramp ramp 0 PULSE(0 2 0 {0.97*tclk} {0.03*tclk} 1n {tclk})\nVctrl ctrlnode 0 DC 0.8 AC 1\n")
SYNC_579 = ("* synchronous buck, two ideal switches with dead time, body diode\n.param fs=500k tclk={1/fs}\nVin vin 0 DC 12\n"
            "Shs vin sw ctrlnode ramp SWMOD\nSls sw 0 ramp ctrlp SWMOD\nVdt ctrlp ctrlnode DC 0.06\n.model SWMOD SW(vt=0 vh=0.005 ron=10m roff=1Meg)\n"
            "D1 0 sw DMOD\n.model DMOD D(IS=1e-9 RS=0.05 N=1.05 CJO=%s VJ=0.7 M=0.5 TT=%s)\nL1 sw swl 10u IC=0\nRdcr swl vout 10m\nResr vout cx 10m\n"
            "Cout cx 0 100u IC=0\nRload vout 0 3.3\nVramp ramp 0 PULSE(0 2 0 {0.97*tclk} {0.03*tclk} 1n {tclk})\nVctrl ctrlnode 0 DC 0.55 AC 1\n")
LINE_579 = (("Vin vin 0 DC {vin_v}", "Vin vin 0 DC {vin_v} AC 1"), ("Vctrl ctrlnode 0 DC {vctrl_dc} AC 1", "Vctrl ctrlnode 0 DC {vctrl_dc}"))
D575.update({
    "bst_cl": ("@", BOOST_579 % ("0", "0"), {20000: (-1.68, 168.8104), 50000: (-16.544, 146.5291), 100000: (-25.1976, 130.4729)}, "boost, chargeless diode (was 71 dB off)"),
    "bst_cj": ("@", BOOST_579 % ("100p", "0"), {20000: (-1.6801, 168.8189), 50000: (-16.5423, 146.5187), 100000: (-25.1942, 130.4604)}, "boost, diode cj 100 pF"),
    "bst_cj1n": ("@", BOOST_579 % ("1n", "0"), {20000: (-1.6808, 168.8891), 50000: (-16.5313, 146.4667), 100000: (-25.1729, 130.3961)}, "boost, diode cj 1 nF"),
    "bst_q": ("@", BOOST_579 % ("100p", "2n"), {20000: (-1.694, 168.7673), 50000: (-16.5474, 146.4393), 100000: (-25.1912, 130.3889)}, "boost, diode cj 100 pF tt 2 ns (reverse recovery through the switch)"),
    "syn_cl": ("@", SYNC_579 % ("0", "0"), {20000: (-7.806, -169.8188), 50000: (-23.8541, -161.4173), 100000: (-34.9243, -147.2917)}, "synchronous buck with dead time, chargeless body diode"),
    "syn_q": ("@", SYNC_579 % ("300p", "2n"), {20000: (-7.8074, -169.7001), 50000: (-23.8548, -161.3819), 100000: (-34.9253, -147.291)}, "synchronous buck with dead time, charged body diode (was 0.25 dB off)"),
    "od_cl_line": ("pac_r6od_ctrl2out.cir", (("csw=1n cjo=300p tt=2n", "csw=0 cjo=0 tt=0"), NOSNUB) + LINE_579, {20000: (-34.5405, -167.3775), 50000: (-50.5727, -160.4842), 100000: (-61.6409, -146.8183)}, "line-to-output, chargeless"),
    "od_q_d30": ("pac_r6od_ctrl2out.cir", (("csw=1n", "csw=0"), NOSNUB, ("vctrl_dc=0.55", "vctrl_dc=1.4")), {20000: (-7.3256, -169.1488), 50000: (-23.3707, -161.1689), 100000: (-34.4408, -147.1773)}, "duty 0.7, fast node"),
})


def checks_579():
    for n in ("bst_cl", "bst_cj", "bst_cj1n", "bst_q", "syn_cl", "syn_q", "od_cl_line", "od_q_d30"):
        check_575(n)
    check_575("bst_q", reltol="1e-6")
    check_575("syn_q", reltol="1e-6")


def red_proofs_579():
    check_575("bst_cl", env={"PSS_SALT_NO579": "1"}, red_min=10.0, tag=" [RED PROOF: PSS_SALT_NO579=1 (the pre-event node equation without the anode-side diode)]")
    check_575("bst_cl", env={"PSS_SALT_NODIOFAR": "1"}, red_min=0.2, tag=" [RED PROOF: PSS_SALT_NODIOFAR=1 (no jump on the capacitor the diode feeds)]")
    check_575("bst_cl", env={"PAC_NO_IMPULSE": "1"}, red_min=0.5, red_freq=100000, tag=" [RED PROOF: PAC_NO_IMPULSE=1 (no output impulse), 100 kHz]")
    check_575("bst_q", env={"PSS_SALT_NOCONT": "1"}, red_min=1.0, tag=" [RED PROOF: PSS_SALT_NOCONT=1]")


# ---------------------------------------------------------------------------------------------------------------------------------
# cards #580 / #581 (2026-10-07).  #580: a capacitor at a boost's switch node (the diode feeds the output capacitor): fast-node treatment
# instead of the capacitor jump (0.14-0.18 dB / 0.55-0.74 deg before).  #581: discontinuous conduction is first order in the step (the diode's
# turn-off instant is not located): 0.40 dB at 200 points per period, 0.06 dB at 2000; a note is printed once.
BOOSTC_580 = BOOST_579.replace("Cout vout cx 22u IC=0", "Csw sw 0 %s\nCout vout cx 22u IC=0")
D575.update({
    "bst_cl_csw": ("@", BOOSTC_580 % ("0", "0", "1n"), {20000: (-1.6851, 169.0429), 50000: (-16.532, 146.3875), 100000: (-25.1428, 130.1333)}, "boost, chargeless diode, 1 nF at the switch node"),
    "bst_q_csw": ("@", BOOSTC_580 % ("100p", "2n", "1n"), {20000: (-1.6989, 169.0004), 50000: (-16.5319, 146.3029), 100000: (-25.1456, 130.1234)}, "boost, cj 100 pF tt 2 ns, 1 nF at the switch node"),
    "bst_cl_csw100p": ("@", BOOSTC_580 % ("0", "0", "100p"), {20000: (-1.6795, 168.8301), 50000: (-16.5465, 146.5302), 100000: (-25.1741, 130.3579)}, "boost, chargeless diode, 100 pF at the switch node"),
    "od_cl_dcm": ("pac_r6od_ctrl2out.cir", (("csw=1n cjo=300p tt=2n", "csw=0 cjo=0 tt=0"), NOSNUB, ("rload=3.3", "rload=60")), {20000: (-31.1878, -84.6675), 50000: (-38.8166, -77.5846), 100000: (-43.84, -68.0311)}, "buck in discontinuous conduction (60 ohm)"),
})
DCM_NOTE = "discontinuous conduction"


def checks_580_581():
    check_575("bst_cl_csw")
    check_575("bst_q_csw")
    check_575("bst_cl_csw100p")
    check_575("bst_q_csw", reltol="1e-6")
    # #581: converged with the points, and the note is there; no note in continuous conduction; PAC_NO_DCM_NOTE=1 silences it
    check_575("od_cl_dcm", points=2000, expect=DCM_NOTE)
    check_575("od_cl_dcm", points=200, bound_db=0.06, bound_deg=0.2, expect=DCM_NOTE, tag=" (200 points: the note must be printed; within 0.03 dB since #711, was 0.40 dB)", freqs=(20000,))
    check_575("od_cl_nosnub", absent=DCM_NOTE, tag=" (continuous conduction: no DCM note)", freqs=(20000,))
    check_575("bst_cl", absent=DCM_NOTE, tag=" (continuous conduction: no DCM note)", freqs=(20000,))
    check_575("od_cl_dcm", points=2000, env={"PAC_NO_DCM_NOTE": "1"}, absent=DCM_NOTE, tag=" (PAC_NO_DCM_NOTE=1)", freqs=(20000,))


def red_proofs_580_581():
    check_575("bst_cl_csw", env={"PSS_SALT_NOFARFAST": "1"}, red_min=0.12, red_freq=100000, tag=" [RED PROOF: PSS_SALT_NOFARFAST=1 (the capacitor jump at a boost node), 100 kHz]")
    check_575("od_cl_dcm", points=200, env={"PSS_SALT_NO711": "1"}, red_min=0.3, tag=" [RED PROOF: PSS_SALT_NO711=1 (switch turn-on reads the inductor current of the trial point), 200 points in DCM must be > 0.3 dB off]")


# ---------------------------------------------------------------------------------------------------------------------------------
# card #582 (2026-10-07): (a) a half bridge whose two ideal switches flip at the same instant (no dead time), no capacitor at the switch node:
# one event per node per step (td .pac was 94 dB off with a chargeless body diode, 12.6 dB with a charged one, .pstb 42 dB closed loop);
# without any diode at the node the legacy per-C/L recursion stays 12 dB off, so circuits with a capacitor-less switch node use the generic
# recursion; (b) closed-loop .pstb on the fast-node and routed treatments vs tools/pac_perturbed.py --loop-gain (tran backend, 0.01 V;
# ~/Work/data/ngspice_augmentation/582/m; the tool refused 10 and 100 kHz on these decks, 30 kHz only).
SYNC0_582 = SYNC_579.replace("Vdt ctrlp ctrlnode DC 0.06", "Vdt ctrlp ctrlnode DC 0")
D575.update({
    "syn0_cl": ("@", SYNC0_582 % ("0", "0"), {20000: (-7.805, -169.9099), 50000: (-23.8534, -161.4507), 100000: (-34.9237, -147.3069)}, "half bridge without dead time, chargeless body diode (was 94 dB off)"),
    "syn0_nodiode": ("@", (SYNC0_582 % ("0", "0")).replace("D1 0 sw DMOD\n", ""), {20000: (-7.805, -169.9099), 50000: (-23.8534, -161.4507), 100000: (-34.9237, -147.3069)}, "half bridge without dead time, no diode"),
    "syn0_q": ("@", SYNC0_582 % ("300p", "2n"), {20000: (-7.805, -169.9108), 50000: (-23.8534, -161.4501), 100000: (-34.9237, -147.3069)}, "half bridge without dead time, charged body diode (was 12.6 dB off)"),
})
CL_582 = {   # name: (replacements on the handbook's closed-loop R6s deck, reference T at 30 kHz (dB, deg), what)
    "cl_sync_nocap": ((("Csw sw 0 1n\n", ""), ("Rsnub sw nsnub 10\n", ""), ("Csnub nsnub 0 1n\n", "")), (1.463, -139.18), "closed-loop half bridge, no node capacitor (was 42 dB off)"),
    "cl_async_q": ((("Sls sw 0   ramp opout SWMOD\n", "D1 0 sw DMOD\n.model DMOD D(IS=1e-9 RS=0.05 N=1.05 CJO=300p VJ=0.7 M=0.5 TT=2n)\n"),
                    ("Csw sw 0 1n\n", ""), ("Rsnub sw nsnub 10\n", ""), ("Csnub nsnub 0 1n\n", "")), (1.910, -137.76), "closed-loop asynchronous buck, charged diode, no node capacitor (fast node)"),
}


def check_pstb_582(name, env=None, red_min=None, bound_db=0.2, bound_deg=1.0):
    repl, (rdb, rdeg), what = CL_582[name]
    label = f"closed-loop .pstb vs pac_perturbed --loop-gain at 30 kHz: {what} [{name}]" + (" [RED PROOF: must exceed %g dB]" % red_min if red_min else "")
    try:
        b = open(os.path.join(ROOT, "docs", "examples", "handbook", "23_pstb_closed_loop_sw_ref_deck.cir")).read()
        for a, c in repl:
            assert a in b, a
            b = b.replace(a, c)
        b = b.replace("\n.end\n", "\n") if "\n.end\n" in b else b
        b += 'Vprobe fbtop vout DC 0\n.pstb 500e3 300e-6 vout 200 10 150 1e-6 lin 1 30000 30000 probe="Vprobe" method="shoot"\n' + wrline("cl_%s.txt" % name, "t")
        txt, rc = run(b, "cl_%s" % name, env=env)
        if "REFUSED" in txt:
            raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l][0][:200])
        d, p_ = dbdeg(rd("cl_%s.txt" % name)[0][1])
    except (CouldNotRun, OSError, IndexError, AssertionError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])
        return
    dp = (p_ - rdeg + 180.0) % 360.0 - 180.0
    if red_min:
        record(label, "pass" if abs(d - rdb) > red_min else "FAIL", f"{d - rdb:+.3f} dB {dp:+.2f} deg")
    else:
        record(label, "pass" if (abs(d - rdb) <= bound_db and abs(dp) <= bound_deg) else "FAIL", f"T {d:.3f} dB / {p_:.2f} deg vs {rdb} / {rdeg}: {d - rdb:+.3f} dB {dp:+.2f} deg (bound {bound_db} / {bound_deg})")


def checks_582():
    for n in ("syn0_cl", "syn0_nodiode", "syn0_q"):
        check_575(n, bound_db=0.05, bound_deg=0.2)
    check_pstb_582("cl_sync_nocap")
    check_pstb_582("cl_async_q")


def red_proofs_582():
    check_575("syn0_cl", env={"PSS_SALT_NO582": "1"}, red_min=10.0, tag=" [RED PROOF: PSS_SALT_NO582=1 (one routed jump per switch)]")
    check_575("syn0_nodiode", env={"PAC_NO582_GENERIC": "1"}, red_min=3.0, tag=" [RED PROOF: PAC_NO582_GENERIC=1 (the legacy per-C/L recursion)]")
    check_pstb_582("cl_sync_nocap", env={"PSS_SALT_NO582": "1"}, red_min=10.0)


# ---------------------------------------------------------------------------------------------------------------------------------
# card #583 (2026-10-07): routed jumps are added after T is formed (no order-1 window for them): the fast-node phase error no longer
# grows with the node capacitance.  Inverting buck-boost (ideal high-side switch, inductor to ground, diode from the negative output to
# the switch node, 20 V swing at the node), with and without 1 nF at the node; references ~/Work/data/ngspice_augmentation/579/m/ref_bb_*
# (the tool refused 20 kHz on the cj + tt deck without node capacitor: its own amplitude-halving check, 1.35 deg).
BB_583 = ("* inverting buck-boost\n.param fs=500k tclk={1/fs}\nVin vin 0 DC 12\nS1 vin sw ctrlnode ramp SWMOD\n.model SWMOD SW(vt=0 vh=0.005 ron=10m roff=1Meg)\n"
          "L1 sw lg 10u IC=0\nRdcr lg 0 10m\nD1 vout sw DMOD\n.model DMOD D(IS=1e-9 RS=0.05 N=1.05 CJO=%s VJ=0.7 M=0.5 TT=%s)\n%s"
          "Cout vout cx 22u IC=0\nResr cx 0 10m\nRload vout 0 10\nVramp ramp 0 PULSE(0 2 0 {0.97*tclk} {0.03*tclk} 1n {tclk})\nVctrl ctrlnode 0 DC 0.8 AC 1\n")
D575.update({
    "bb_cl": ("@", BB_583 % ("0", "0", ""), {20000: (5.6601, -0.2818), 50000: (-10.6352, -12.0204), 100000: (-21.6039, -23.6838)}, "buck-boost, chargeless diode"),
    "bb_q": ("@", BB_583 % ("100p", "2n", ""), {50000: (-10.6479, -12.0557), 100000: (-21.6116, -23.7504)}, "buck-boost, cj 100 pF tt 2 ns"),
    # card #728: the 20 kHz points of the two buck-boost decks with 1 nF at the node are the exact values of tools/pac_exact.py (the
    # perturbed-PSS references were 0.105 / 0.108 deg above them; at 50 / 100 kHz they are within 0.03 / 0.04 deg of the exact values)
    "bb_cl_csw": ("@", BB_583 % ("0", "0", "Csw sw 0 1n\n"), {20000: (5.6279, 1.0996), 50000: (-10.6306, -11.9392), 100000: (-21.5722, -24.3383)}, "buck-boost, chargeless diode, 1 nF at the node (was 1.06 deg off)"),
    "bb_q_csw": ("@", BB_583 % ("100p", "2n", "Csw sw 0 1n\n"), {20000: (5.6113, 1.1092), 50000: (-10.644, -12.002), 100000: (-21.5716, -24.4017)}, "buck-boost, cj + tt, 1 nF at the node (was 1.38 deg off)"),
    "bst_q_csw3n": ("@", BOOSTC_580 % ("100p", "2n", "3n"), {20000: (-1.7107, 169.4522), 50000: (-16.4893, 145.9996), 100000: (-25.0794, 129.4944)}, "boost, cj + tt, 3 nF at the node (was 1.04 deg off)"),
})


def checks_583():
    check_575("bb_cl")
    check_575("bb_q", freqs=(50000, 100000))
    check_575("bb_cl_csw", bound_deg=0.4)
    check_575("bb_q_csw", bound_deg=0.4)
    check_575("bst_q_csw3n", bound_deg=0.45)
    # the decks whose phase was at 0.14-0.34 deg with the order-1 window: tighter now
    check_575("swc_1n_0", bound_db=0.05, bound_deg=0.2, tag=" (phase 0.34 deg before #583)")
    check_575("odm_nosnub", bound_db=0.05, bound_deg=0.1, tag=" (phase 0.14 deg before #583)")
    check_575("od_nosnub", bound_db=0.05, bound_deg=0.1, tag=" (phase 0.15 deg before #583)")


def red_proofs_583():
    check_575("bb_q_csw", env={"PAC_SALT_SMOOTHT": "0"}, red_min=0.0, red_deg=1.0, tag=" [RED PROOF: PAC_SALT_SMOOTHT=0 (jumps before T, order-1 window at every continuation step): phase]")
    check_575("swc_1n_0", env={"PAC_SALT_SMOOTHT": "0"}, red_min=0.0, red_deg=0.25, tag=" [RED PROOF: PAC_SALT_SMOOTHT=0: phase]")


# ---------------------------------------------------------------------------------------------------------------------------------
# card #586 (2026-10-07): a source that is not periodic from t = 0 (a PULSE sawtooth with a delay: the natural way to write the second ramp
# of an interleaved converter) made the shooting engine converge to the orbit of another circuit, silently (td .pac 6.5 dB off).  It is
# refused now; the same two-phase interleaved buck with a ramp that IS periodic from 0 (unshifted sawtooth + square source in series) is
# checked against tools/pac_perturbed.py (~/Work/data/ngspice_augmentation/579/m/ref_il2c_cl).
IL2_586 = ("* two-phase interleaved asynchronous buck\n.param fs=500k tclk={1/fs}\nVin vin 0 DC 12\nS1 vin sw1 ctrlnode ramp1 SWMOD\nS2 vin sw2 ctrlnode ramp2 SWMOD\n"
           ".model SWMOD SW(vt=0 vh=0.005 ron=10m roff=1Meg)\nD1 0 sw1 DMOD\nD2 0 sw2 DMOD\n.model DMOD D(IS=1e-9 RS=0.05 N=1.05 CJO=0 VJ=0.7 M=0.5 TT=0)\n"
           "L1 sw1 l1x 10u IC=0\nR1 l1x vout 100m\nL2 sw2 l2x 10u IC=0\nR2 l2x vout 100m\nResr vout cx 10m\nCout cx 0 100u IC=0\nRload vout 0 1\n"
           "Vramp1 ramp1 0 PULSE(0 2 0 {0.97*tclk} {0.03*tclk} 1n {tclk})\n%s\nVctrl ctrlnode 0 DC 0.55 AC 1\n")
IL2_GOOD = "Vramp2 ramp2 rsq PULSE(0 2 0 {0.97*tclk} {0.03*tclk} 1n {tclk})\nVsq rsq 0 PULSE(1 -1 {0.485*tclk} 1n 1n {0.5*tclk} {tclk})"
IL2_DELAYED = "Vramp2 ramp2 0 PULSE(0 2 {tclk/2} {0.97*tclk} {0.03*tclk} 1n {tclk})"
REF_IL2 = {20000: (-0.9619, -158.762), 50000: (-17.3783, -157.492), 100000: (-28.5003, -145.356)}
D575.update({
    "il2c_cl": ("@", IL2_586 % IL2_GOOD, REF_IL2, "two-phase interleaved buck, ramps periodic from t = 0"),
    "il2_delayed": ("@", IL2_586 % IL2_DELAYED, REF_IL2, "two-phase interleaved buck, second ramp as a delayed PULSE"),
})


def checks_586():
    check_575("il2c_cl", bound_db=0.05, bound_deg=0.2)
    name = "guard: a PULSE that is not periodic from t = 0 (delayed sawtooth) is REFUSED by the shooting engine"
    try:
        txt, rc = run(_deck_575("il2_delayed", 20000, "1e-4", "il2d.txt"), "il2d")
        record(name, "pass" if ("REFUSED" in txt and "not periodic from t = 0" in txt and "vramp2" in txt.lower()) else "FAIL",
               [l for l in txt.split("\n") if "REFUSED" in l][:1] or txt[-300:])
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
    # a delayed gate pulse that fits in its period (TD + TR + PW + TF <= PER) is periodic from 0: must not be refused
    name = "a delayed PULSE that fits in its period is not refused"
    try:
        d = (IL2_586 % IL2_GOOD).replace("Rload vout 0 1\n", "Rload vout 0 1\nVaux aux 0 PULSE(0 1 {0.2*tclk} 1n 1n {0.3*tclk} {tclk})\nRaux aux 0 1k\n")
        d += '.options reltol=1e-4\n.pac 500e3 100e-6 vout 200 10 150 1e-6 lin 1 20000 20000 method="shoot" solver="td"\n' + wrline("il2aux.txt", "v(vout)")
        txt, rc = run(d, "il2aux")
        record(name, "pass" if ("REFUSED" not in txt and "PAC-TD: done" in txt) else "FAIL", [l for l in txt.split("\n") if "REFUSED" in l][:1] or "ran")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))


def red_proofs_586():
    check_575("il2_delayed", env={"PSS_ALLOW_DELAYED_SOURCE": "1"}, red_min=3.0, tag=" [RED PROOF: PSS_ALLOW_DELAYED_SOURCE=1 (the delayed ramp accepted as before): must be > 3 dB off]")


# ---------------------------------------------------------------------------------------------------------------------------------
# card #713 (2026-10-19): a switch node whose only charge is the devices' together with an RC snubber (refused since #575).  Fast node (routed
# jump continued until the node has settled) plus one jump of the snubber capacitor when the continuation ends: G / C_far x (D_end - D_before),
# D = v_node - v_far (the saltation rule across the fast transition; a chargeless node gives G / C_far x dV, the jump of #575).
# References: tools/pac_perturbed.py, same settings as the #575 matrix (pss backend, 0.05 V, 800 points per cycle, halving and step checks),
# s143 binary, one run per deck (~/Work/data/ngspice_augmentation/713/ref_*; od_snub, odm_snub, od_snub_r50 and od_snub_c100p are identical
# to the s117 references of #575, 575/m/ref_*).  Closed loop: pac_perturbed --loop-gain --probe-nodes fbtop:vout --amplitude 0.01
# --backend pss (713/ref_cl_async_q_snub_pss; halving <= 0.0003 dB; the tran backend's halving is 0.05 dB there and it refused 10 kHz).
SNUB = (("csw=1n", "csw=0"),)
D575.update({
    "od_snub": ("pac_r6od_ctrl2out.cir", SNUB, {20000: (-7.4191, -166.6055), 50000: (-23.4423, -160.2547), 100000: (-34.5099, -146.7526)}, "charge-only node + RC snubber 10 ohm / 1 nF (refused #575..#712)"),
    "od_snub_r50": ("pac_r6od_ctrl2out.cir", SNUB + (("Rsnub sw nsnub 10", "Rsnub sw nsnub 50"),), {20000: (-7.4076, -167.167), 50000: (-23.4375, -160.4308), 100000: (-34.5058, -146.8146)}, "charge-only node + snubber 50 ohm / 1 nF"),
    "od_snub_c100p": ("pac_r6od_ctrl2out.cir", SNUB + (("Csnub nsnub 0 1n", "Csnub nsnub 0 100p"),), {20000: (-7.4091, -167.0899), 50000: (-23.4382, -160.405), 100000: (-34.5064, -146.8065)}, "charge-only node + snubber 10 ohm / 100 pF"),
    "od_snub_slow": ("pac_r6od_ctrl2out.cir", SNUB + (("Rsnub sw nsnub 10", "Rsnub sw nsnub 50"), ("Csnub nsnub 0 1n", "Csnub nsnub 0 10n")), {20000: (-7.4035, -167.1392), 50000: (-23.4334, -160.4298), 100000: (-34.5024, -146.8288)}, "charge-only node + snubber 50 ohm / 10 nF (tau 0.5 us)"),
    "od_snub_d30": ("pac_r6od_ctrl2out.cir", SNUB + (("vctrl_dc=0.55", "vctrl_dc=1.4"),), {20000: (-7.3256, -169.1341), 50000: (-23.3706, -161.1644), 100000: (-34.4406, -147.1804)}, "charge-only node + snubber, duty 0.7"),
    "od_snub_q5n": ("pac_r6od_ctrl2out.cir", (("csw=1n cjo=300p tt=2n", "csw=0 cjo=50p tt=5n"),), {20000: (-7.4132, -166.8925), 50000: (-23.4394, -160.3547), 100000: (-34.5076, -146.7742)}, "diode cj 50 pF tt 5 ns + snubber 10 ohm / 1 nF"),
    "odm_snub": ("pac_r6odm_ctrl2out.cir", SNUB, {20000: (-7.4484, -168.4387), 50000: (-23.4882, -160.9517), 100000: (-34.5574, -147.1024)}, "MOS body diode + overlap, snubber 10 ohm / 1 nF"),
    "odm_snub_r50": ("pac_r6odm_ctrl2out.cir", SNUB + (("Rsnub sw nsnub 10", "Rsnub sw nsnub 50"),), {20000: (-7.4398, -168.9895), 50000: (-23.4836, -161.125), 100000: (-34.5536, -147.1637)}, "MOS body diode + overlap, snubber 50 ohm / 1 nF"),
    # 10 nF snubber: at 0.05 V the tool's 20 kHz point is amplitude-limited (halving 0.75 deg); 20 / 50 kHz from a 0.0125 V run (halving 0.047 /
    # 0.003 deg, 713/ref_od_snub_c10n_a0125), 100 kHz from the 0.05 V run (halving 0.002 deg)
    "od_snub_c10n": ("pac_r6od_ctrl2out.cir", SNUB + (("Csnub nsnub 0 1n", "Csnub nsnub 0 10n"),), {20000: (-7.473, -164.1987), 50000: (-23.4628, -159.3648), 100000: (-34.5265, -146.3928)}, "charge-only node + snubber 10 ohm / 10 nF"),
    # discontinuous conduction (60 ohm), reference at 0.005 V (the tool refused 0.05 V: halving 0.05-0.75 dB), 713/ref_od_snub_dcm_a005; run at 2000
    # points per period (at 200: 0.072 dB / 0.21 deg, like the other DCM + snubber decks: charge jump 0.047 dB / 0.39 deg, chargeless 0.070 / 0.25)
    "od_snub_dcm": ("pac_r6od_ctrl2out.cir", SNUB + (("rload=3.3", "rload=60"),), {20000: (-32.2538, -82.5734), 50000: (-39.6844, -72.5753), 100000: (-44.0374, -59.5356)}, "charge-only node + snubber, discontinuous conduction (60 ohm)"),
})
DECKS_713 = ("od_snub", "od_snub_r50", "od_snub_c100p", "od_snub_slow", "od_snub_c10n", "od_snub_d30", "od_snub_q5n", "odm_snub", "odm_snub_r50")
CL_713 = {   # closed loop (handbook R6s): asynchronous, diode cj 300 pF tt 2 ns, no node capacitor, the deck's 10 ohm / 1 nF snubber kept
    "cl_async_q_snub": ((("Sls sw 0   ramp opout SWMOD\n", "D1 0 sw DMOD\n.model DMOD D(IS=1e-9 RS=0.05 N=1.05 CJO=300p VJ=0.7 M=0.5 TT=2n)\n"), ("Csw sw 0 1n\n", "")),
                        {10000: (16.1583, -120.757), 30000: (1.884, -137.428), 100000: (-16.1719, -162.832)},
                        "closed-loop asynchronous buck, charged diode, RC snubber, no node capacitor"),
}


def check_pstb_713(name, bound_db=0.067, bound_deg=0.32, env=None):
    repl, ref, what = CL_713[name]
    label = f"closed-loop .pstb vs pac_perturbed --loop-gain (pss backend): {what} [{name}]"
    parts, wdb, wdeg = [], 0.0, 0.0
    try:
        b = open(os.path.join(ROOT, "docs", "examples", "handbook", "23_pstb_closed_loop_sw_ref_deck.cir")).read()
        for a, c in repl:
            assert a in b, a
            b = b.replace(a, c)
        b = b.replace("\n.end\n", "\n")
        for f in sorted(ref):
            o = "cl713_%s_%d.txt" % (name, f)
            d_ = b + 'Vprobe fbtop vout DC 0\n.pstb 500e3 300e-6 vout 200 10 150 1e-6 lin 1 %g %g probe="Vprobe" method="shoot"\n' % (f, f) + wrline(o, "t")
            txt, rc = run(d_, "cl713_%s_%d" % (name, f), env=env)
            if "REFUSED" in txt:
                raise CouldNotRun([l for l in txt.split("\n") if "REFUSED" in l][0][:200])
            d, p_ = dbdeg(rd(o)[0][1])
            dp = (p_ - ref[f][1] + 180.0) % 360.0 - 180.0
            wdb, wdeg = max(wdb, abs(d - ref[f][0])), max(wdeg, abs(dp))
            parts.append(f"{f / 1e3:g}k {d - ref[f][0]:+.4f} dB {dp:+.3f} deg")
    except (CouldNotRun, OSError, IndexError, AssertionError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])
        return
    record(label, "pass" if (wdb <= bound_db and wdeg <= bound_deg) else "FAIL", "; ".join(parts) + f"; worst {wdb:.4f} dB / {wdeg:.3f} deg (bound {bound_db} / {bound_deg})")


def checks_713():
    for n in DECKS_713:
        check_575(n, bound_db=0.067, bound_deg=0.32, expect="card #713")
    check_575("od_snub_dcm", points=2000, bound_db=0.067, bound_deg=0.32, expect="card #713")
    check_575("od_snub", reltol="1e-6", bound_db=0.067, bound_deg=0.32)
    check_575("odm_snub", reltol="1e-6", bound_db=0.067, bound_deg=0.32)
    check_pstb_713("cl_async_q_snub")


def red_proofs_713():
    # PSS_SALT_NO713=1: the #575 refusal is back (the accuracy check cannot run without the feature)
    global QUIET
    QUIET = True
    check_575("od_snub", env={"PSS_SALT_NO713": "1"}, tag=" [RED PROOF NO713]", freqs=(20000,))
    r = results.pop() if results and "RED PROOF NO713" in results[-1][0] else None
    QUIET = False
    record("red proof: PSS_SALT_NO713=1 restores the refusal of a charge-only switch node with an RC snubber",
           "pass" if (r and r[1] == "COULD-NOT-RUN" and "REFUSED" in r[2] and "RC snubber" in r[2]) else "FAIL", r[2][:200] if r else "no result")
    # PAC_SALT_FASTSNUB=1: the fast node without the snubber jump (what PAC_SALT_FASTSNUB ran before #713)
    check_575("od_snub_c10n", env={"PAC_SALT_FASTSNUB": "1"}, red_min=0.045, red_freq=50000,
              tag=" [RED PROOF: PAC_SALT_FASTSNUB=1 (no snubber jump), 50 kHz: must exceed 0.045 dB (0.020 with it)]")
    check_575("od_snub_slow", env={"PAC_SALT_FASTSNUB": "1"}, red_min=0.008, red_freq=50000,
              tag=" [RED PROOF: PAC_SALT_FASTSNUB=1, 50 kHz: must exceed 0.008 dB (0.0004 with the jump)]")
    # PAC_NO713_COEF=1: continuation jumps take the event coefficient of the step they are given in (before #713): the slow DCM turn-off edge
    check_575("od_snub_dcm", points=2000, env={"PAC_NO713_COEF": "1"}, red_min=0.0, red_deg=0.3, red_freq=100000,
              tag=" [RED PROOF: PAC_NO713_COEF=1 (continuation coefficient at its own step), 100 kHz phase must exceed 0.3 deg (0.02 with the flip's)]")


def checks_559():
    check_salt_model_charge("pac_r6od_ctrl2out.cir", PP_OD, "diode cj+tt, 1 nF")
    check_salt_model_charge("pac_r6odm_ctrl2out.cir", PP_ODM, "MOS body diode + cgdo, 1 nF")
    check_salt_native_window()
    check_salt_native_r6o()
    check_salt_refusals()
    check_nqs_lti("bsim4", 1, 1, FRA_NQS4, "BSIM4 trnqsmod=1 acnqsmod=1, LTI bias, 1/3/10/30 MHz")
    check_nqs_lti("bsim4", 1, 0, FRA_QS4, "BSIM4 acnqsmod=1 only (td = transient model, quasi-static), 10 MHz")
    check_nqs_lti("bsim3", 0, 1, FRA_NQS3, "BSIM3v3.3 nqsmod=1, LTI bias, 3/10/30 MHz")
    check_nqs_lptv()
    check_nqs_pxf()
    check_nqs_refusals()
    check_nqs_instance()


def red_proofs_559():
    """card #559 red proofs: each new check must be able to fail."""
    check_salt_model_charge("pac_r6od_ctrl2out.cir", PP_OD, "diode cj+tt, 1 nF", env={"PAC_SALT_NATIVEJUMP": "1"}, tag=" [RED PROOF: PAC_SALT_NATIVEJUMP=1 (the jump divided by the native capacitance only, no refusal)]", red_min=0.3)
    check_salt_model_charge("pac_r6odm_ctrl2out.cir", PP_ODM, "MOS body diode + cgdo, 1 nF", env={"PAC_SALT_NATIVEJUMP": "1"}, tag=" [RED PROOF: PAC_SALT_NATIVEJUMP=1]", red_min=1.0)
    # PAC_SALT_NOQ=1: the pre-#559 refusal is back (check must report COULD-NOT-RUN with REFUSED, i.e. the live check cannot pass without the feature)
    global QUIET
    QUIET = True
    check_salt_model_charge("pac_r6od_ctrl2out.cir", PP_OD, "diode cj+tt, 1 nF", env={"PAC_SALT_NOQ": "1"}, tag=" [RED PROOF NOQ]")
    r = results.pop() if results and "RED PROOF NOQ" in results[-1][0] else None
    QUIET = False
    record("red proof: PAC_SALT_NOQ=1 restores the refusal (the saltation check cannot run without the feature)",
           "pass" if (r and r[1] == "COULD-NOT-RUN" and "REFUSED" in r[2]) else "FAIL", r[2][:160] if r else "no result")
    check_salt_native_window(env={"PSS_SALT_MAXSTEPS": "-1"}, red_min=3.0)
    # a quasi-static deck against the NQS reference must fail (the NQS checks discriminate): BSIM4 trnqsmod=0 at 10 MHz is 40 % off the transient NQS FRA
    check_nqs_lti("bsim4", 0, 0, {1e7: FRA_NQS4[1e7]}, "trnqsmod=0 (quasi-static) deck vs the NQS FRA, 10 MHz", red_min=0.2)
    check_nqs_lptv(nq=0, red_min=0.3)




def main():
    if "--only" in sys.argv:
        for n in sys.argv[sys.argv.index("--only") + 1].split(","):
            if n in NDECKS:
                check_pnoise_dev(n)
            elif n.startswith("qs_"):
                check_pnoise_qs(n[3:])
            elif n == "r8":
                check_r8_freerun("ctrl2out")
                check_r8_freerun("line2out")
            elif n == "freerun":
                check_freerun_refusal()
            elif n == "gear":
                check_gear()
                check_gear_refusals()
                check_gear_r7()
                check_gear(red=True)
            elif n == "559":
                checks_559()
                red_proofs_559()
            elif n == "575":
                checks_575()
                red_proofs_575()
            elif n == "579":
                checks_579()
                red_proofs_579()
            elif n == "580":
                checks_580_581()
                red_proofs_580_581()
            elif n == "582":
                checks_582()
                red_proofs_582()
            elif n == "583":
                checks_583()
                red_proofs_583()
            elif n == "586":
                checks_586()
                red_proofs_586()
            elif n == "713":
                checks_713()
                red_proofs_713()
            elif n == "salt":
                check_salt_charge_identity()
                check_salt_charge_identity(fixture="pac_r6om_ctrl2out.cir", kind="MOS1", measured="2.0e-4")
                check_salt_charge()
        n = {k: sum(1 for r in results if r[1] == k) for k in ("pass", "FAIL", "COULD-NOT-RUN")}
        print(f"\n{n['pass']} pass, {n['FAIL']} FAIL, {n['COULD-NOT-RUN']} could-not-run (of {len(results)})")
        sys.exit(1 if n["FAIL"] else (2 if n["COULD-NOT-RUN"] else 0))
    for n in DECKS:
        if QUICK and n not in ("mos1_amp", "bsim4_cs_capmod0", "mutual_xfmr"):
            continue
        check_lti(n, bound=BOUNDS.get(n, 1e-6))
    check_chain()
    for n in ("mos1_ladder", "bsim4_cs_capmod0", "mos1_amp"):
        check_pxf_adjoint(n)
    check_pnoise_lti()
    check_refusals()
    check_freerun_refusal()
    check_force_generic()
    check_klu()
    for n in NDECKS:
        check_pnoise_dev(n)
    for k in QS_DECKS:
        check_pnoise_qs(k)
    check_r8_freerun("ctrl2out")
    if not QUICK:
        check_r8_freerun("line2out")
    check_salt_charge_identity()
    check_salt_charge_identity(fixture="pac_r6om_ctrl2out.cir", kind="MOS1", measured="2.0e-4")
    check_salt_charge()
    check_gear()
    check_gear_refusals()
    if not QUICK:
        check_gear_r7()
    check_gear(red=True)
    red_proofs()
    red_proofs_544()
    red_proofs_556()
    checks_559()
    red_proofs_559()
    checks_575()
    red_proofs_575()
    checks_579()
    red_proofs_579()
    checks_580_581()
    red_proofs_580_581()
    checks_582()
    red_proofs_582()
    checks_583()
    red_proofs_583()
    checks_586()
    red_proofs_586()
    checks_713()
    red_proofs_713()
    check_asan()
    n = {k: sum(1 for r in results if r[1] == k) for k in ("pass", "FAIL", "COULD-NOT-RUN")}
    print(f"\n{n['pass']} pass, {n['FAIL']} FAIL, {n['COULD-NOT-RUN']} could-not-run (of {len(results)})")
    sys.exit(1 if n["FAIL"] else (2 if n["COULD-NOT-RUN"] else 0))


if __name__ == "__main__":
    main()
