#!/usr/bin/env python3
"""Validation of td `.pxf` and closed-loop `.pstb` on ideal-switch converters (card #709, 2026-10-19; binary s141).

After the saltation work (cards #568..#586) td `.pac` was validated on 47 open-loop decks against the perturbed periodic steady state
(tools/pac_perturbed.py, tests/test_pac_charge.py).  This suite closes the two remaining gaps:

  1. `.pxf` (periodic transfer function, per-source engine and `adjoint=1`, solver="td") on ideal-switch converters: buck (duty 0.275 and
     0.7), boost (chargeless and charged diode), half bridge with dead time (chargeless and charged body diode), synchronous buck without
     dead time.  Three independent sources per deck: Vin (line), Vctrl (control, the PWM comparator input), Iinj (output current sink).
  2. closed-loop `.pstb` (periodic loop gain, Tian probe `lgprobe`) on a BOOST converter with an ideal switch.

    NGSPICE=<prefix>/bin/ngspice python3 tests/test_pxf_pstb_switch.py [--fra-live]
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN (a FAIL wins).  N/M pass is printed at the end.

Decks: tests/fixtures/pxfpstb709/*.cir (copied from tests/fixtures/dcm711 and refs with an output current sink Iinj added).  Outputs under
~/Work/data/ngspice_augmentation/709/test, never the repository.  Every ngspice runs under systemd-run MemoryMax=3G, OMP_NUM_THREADS=2,
`set num_threads=2`; one at a time.

Checks
  1  `.pxf` (adjoint=0) vs `.pac` driven at the same source, read at the same output: 7 decks x 3 sources x 9 frequencies (20..100 kHz),
     1e-6 relative.  Measured 0.0 (bit-identical: the same engine, `.pxf` reuses pac_engine_run()).
  2  `.pxf adjoint=1`: on every deck with ideal-switch events the adjoint engine REFUSES ("saltation events seen") and falls back to the
     per-source engine; checked: the message is printed and the numbers equal the `.pac` ones (1e-6).  The transposed recursion itself is
     NOT exercised by an ideal-switch deck (docs/pss-pac.md, "Saltation ... adjoint"), it is covered by tests/test_pxf.py on smooth decks.
  3  `.pxf` vs the perturbed periodic steady state (tools/pac_perturbed.py, pinned below): 7 decks x 3 sources x {20, 50, 100} kHz,
     0.07 dB / 0.35 deg.  Measured worst 0.042 dB / 0.09 deg on s141/s142 (the dead-time half bridges, see 4); after card #727 0.034 dB / 0.09 deg (boost_q).
  4  the half bridge with dead time, control-to-output within 0.005 dB of the reference at 20/50/100 kHz (Vin, Iinj too).  Card #709 found
     it +0.0407 dB (chargeless body diode) / +0.035 dB (charged) high at every frequency; card #727 (2026-10-19) traced it to the pre-event
     node equation of the routed jump, which left out the body diode's series resistance (the pre-event voltage of the low-side and
     high-side turn-on after a dead time short by Id x RS), and on the charged diode a continuation cut before the diffusion-charge tail
     (docs/pss-pac.md, card #727).  Measured on the #727 build: syn_dt +0.0000..+0.0001 dB, syn_dt_q -0.0003..+0.0005 dB.  Red proof
     PSS_SALT_NO727=1 brings the +0.04 dB back (in 6).  syn0 (no dead time) stays within 0.005 dB as before.
  5  closed-loop `.pstb` (tests/fixtures/pxfpstb709/cl_boost.cir: boost, ideal low-side `sw`, chargeless diode, integrating compensator,
     fs 500 kHz, Vout 7.7 V, crossover 1.35 kHz, phase margin 85.3 deg), 9 frequencies 1 kHz..fs/4:
     a  vs the FRA (tools/fra.py loop mode, injected-tone transients, pinned below): 0.25 dB / 2.5 deg below fs/10;
     b  vs the perturbed-PSS loop gain (tools/pac_perturbed.py --loop-gain) at 10 and 20 kHz: 0.05 dB / 0.3 deg;
     c  vs the analytic composition T = Zf/R1 x beta x Gvc, Gvc being the `.pxf` control-to-output of the open-loop boost: 0.08 dB / 0.3 deg
        (independent of the FRA at every frequency, including where T is -65 .. -85 dB);
     d  phase margin and crossover frequency printed by `.pstb` vs the FRA interpolated at the crossover.
  6  red proofs (each must break the agreement): PSS_SALTATION=0, PSS_SALT_NO579, PSS_SALT_NODIOFAR, PAC_NO_IMPULSE, PSS_SALT_NO582,
     PSS_SALT_NOCONT, PSS_SALT_NO727 (syn_dt and syn_dt_q off by more than 0.025 dB) on `.pxf`; PSS_SALTATION=0 and PSS_SALT_NO579 on `.pstb`.

References (all computed once, s141, 2026-10-19, ~/Work/data/ngspice_augmentation/709):
  PERTURBED  tools/pac_perturbed.py sweep --deck tests/fixtures/pxfpstb709/<deck>.cir --fs 500e3 --inject source:<S> --outputs vout
             --freqs 20000,50000,100000 --backend pss --amplitude 0.01 (Vctrl, 20/50 kHz; 0.05 at 100 kHz, halving bound 0.3 dB) / 0.2 (Vin, Iinj)
             Default grid 800 points per cycle, halving and 800->1600 step self-checks inside the tool (halving <= 0.018 dB, step <= 0.011 dB).
             Logs: 709/refall.log, 709/ref709.json.
  PSTB LG    tools/pac_perturbed.py sweep --deck <cl_boost without XP> --loop-gain --probe-nodes fbtop:vout --inject node:opout_pac_unused
             --freqs 10000 | 20000 --amplitude 0.2 --backend pss --stabtime 2e-3.
  FRA        tools/fra.py run_loop_point (probe XP), startup uic, .options trtol=1 reltol=1e-4, tmax and amplitude per row below, see FRA below.
"""
import cmath
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FX = os.path.join(HERE, "fixtures", "pxfpstb709")
NG = os.environ.get("NGSPICE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice")
OUT = os.path.expanduser("~/Work/data/ngspice_augmentation/709/test")
LIVE_FRA = "--fra-live" in sys.argv
results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"[{state:13s}] {name}: {detail}", flush=True)


def run(deck, name, env=None, timeout=900):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name + ".cir")
    with open(p, "w") as fh:
        fh.write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    try:
        r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=3G", "-p", "MemorySwapMax=0", NG, "-b", p],
                           capture_output=True, text=True, cwd=OUT, env=e, timeout=timeout)
    except (FileNotFoundError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{name}: {ex}")
    return (r.stdout + r.stderr).replace("\r", "\n")


def rd(fn):
    """wrdata file with complex columns -> [(f, [z...])]"""
    path = os.path.join(OUT, fn)
    if not os.path.exists(path):
        raise CouldNotRun(f"missing {fn}")
    rows = []
    for l in open(path).read().split("\n")[1:]:
        t = l.split()
        if len(t) >= 3:
            v = [float(x) for x in t]
            rows.append((v[0], [complex(v[1 + 2 * i], v[2 + 2 * i]) for i in range((len(v) - 1) // 2)]))
    if not rows:
        raise CouldNotRun(f"empty {fn}")
    return rows


def dbdeg(z):
    return 20 * math.log10(abs(z)), math.degrees(math.atan2(z.imag, z.real))


def dphase(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def fixture(fx):
    return open(os.path.join(FX, fx)).read()


WR = ".control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata %s %s\nquit\n.endc\n.end\n"
SOURCES = ("Vin", "Vctrl", "Iinj")
SWEEP9 = "lin 9 20000 100000"
DECKS = ("boost", "boost_q", "buck", "buck_d70", "syn_dt", "syn_dt_q", "syn0")
WHAT = {"boost": "boost, chargeless diode", "boost_q": "boost, diode cj 100 pF tt 2 ns", "buck": "buck (R6od family), duty 0.275",
        "buck_d70": "buck, duty 0.7", "syn_dt": "half bridge with dead time, chargeless body diode",
        "syn_dt_q": "half bridge with dead time, charged body diode (cj 300 pF, tt 2 ns)", "syn0": "synchronous buck without dead time"}


def pxf_deck(fx, out, adjoint=0, sweep=SWEEP9, srcs=SOURCES, vec=None):
    b = fixture(fx + ".cir")
    card = (".options reltol=1e-4\n.pxf 500e3 100e-6 vout 200 10 150 1e-6 vout %s method=\"shoot\" solver=\"td\" %s adjoint=%d\n"
            % (sweep, " ".join('sources="%s"' % s for s in srcs), adjoint))
    return b + card + WR % (out, " ".join("pxf_" + s.lower() for s in srcs))


def pac_deck(fx, src, out, sweep=SWEEP9):
    import re
    b = re.sub(r"^(%s .*)$" % src, r"\1 AC 1", fixture(fx + ".cir"), flags=re.M)
    card = ".options reltol=1e-4\n.pac 500e3 100e-6 vout 200 10 150 1e-6 %s method=\"shoot\" solver=\"td\"\n" % sweep
    return b + card + WR % (out, "v(vout)")


# ------------------------------------------------------------------------------------------------------------------------------
# perturbed-PSS references (dB, deg of v(vout) / source), tools/pac_perturbed.py, s141, 2026-10-19 (see the docstring for the command)
# ------------------------------------------------------------------------------------------------------------------------------
PERT = {
    "boost": {"Vin": {20000: (-14.3670, -172.9314), 50000: (-31.0455, -174.0511), 100000: (-43.1315, -171.1572)},
              "Vctrl": {20000: (-1.6773, 168.7883), 50000: (-16.5579, 146.5238), 100000: (-25.1976, 130.4729)},
              "Iinj": {20000: (-7.9096, 94.3770), 50000: (-16.6374, 94.8795), 100000: (-22.7038, 98.3202)}},
    "boost_q": {"Vin": {20000: (-14.3833, -172.9167), 50000: (-31.0604, -174.0452), 100000: (-43.1463, -171.1594)},
                "Vctrl": {20000: (-1.6913, 168.7451), 50000: (-16.5617, 146.4336), 100000: (-25.1912, 130.3889)},
                "Iinj": {20000: (-7.9113, 94.3772), 50000: (-16.6376, 94.8793), 100000: (-22.7038, 98.3201)}},
    "buck": {"Vin": {20000: (-34.5405, -167.3766), 50000: (-50.5727, -160.4842), 100000: (-61.6409, -146.8183)},
             "Vctrl": {20000: (-7.4036, -167.3770), 50000: (-23.4359, -160.4859), 100000: (-34.5043, -146.8206)},
             "Iinj": {20000: (-21.3884, 99.3347), 50000: (-29.4737, 108.1942), 100000: (-34.5232, 122.5110)}},
    "buck_d70": {"Vin": {20000: (-26.4468, -169.1699), 50000: (-42.4920, -161.1684), 100000: (-53.5619, -147.1647)},
                 "Vctrl": {20000: (-7.3256, -169.1718), 50000: (-23.3709, -161.1729), 100000: (-34.4409, -147.1752)},
                 "Iinj": {20000: (-21.3841, 99.2242), 50000: (-29.4734, 108.1874), 100000: (-34.5231, 122.5102)}},
    "syn_dt": {"Vin": {20000: (-34.5245, -169.8190), 50000: (-50.5726, -161.4177), 100000: (-61.6429, -147.2911)},
               "Vctrl": {20000: (-7.8060, -169.8180), 50000: (-23.8541, -161.4184), 100000: (-34.9243, -147.2917)},
               "Iinj": {20000: (-21.3828, 99.1839), 50000: (-29.4732, 108.1850), 100000: (-34.5231, 122.5099)}},
    "syn_dt_q": {"Vin": {20000: (-34.5150, -169.6950), 50000: (-50.5628, -161.3778), 100000: (-61.6343, -147.2716)},
                 "Vctrl": {20000: (-7.8073, -169.6995), 50000: (-23.8552, -161.3839), 100000: (-34.9253, -147.2910)},
                 "Iinj": {20000: (-21.3832, 99.1921), 50000: (-29.4732, 108.1855), 100000: (-34.5231, 122.5099)}},
    "syn0": {"Vin": {20000: (-34.5231, -169.9100), 50000: (-50.5716, -161.4504), 100000: (-61.6417, -147.3080)},
             "Vctrl": {20000: (-7.8049, -169.9100), 50000: (-23.8534, -161.4506), 100000: (-34.9237, -147.3069)},
             "Iinj": {20000: (-21.3827, 99.1782), 50000: (-29.4732, 108.1846), 100000: (-34.5231, 122.5098)}},
}
PXF_DB, PXF_DEG = 0.07, 0.35


def check_pxf_vs_pac(d):
    label = f"pxf vs pac (same source, same output) 1e-6: {WHAT[d]} [{d}]"
    try:
        o0 = run(pxf_deck(d, d + "_pxf0.txt", 0), d + "_pxf0")
        p0 = rd(d + "_pxf0.txt")
        worst = 0.0
        for i, s in enumerate(SOURCES):
            run(pac_deck(d, s, d + "_pac_%s.txt" % s), d + "_pac_" + s)
            pc = rd(d + "_pac_%s.txt" % s)
            if len(pc) != len(p0):
                raise CouldNotRun("point counts differ")
            for a, b in zip(p0, pc):
                worst = max(worst, abs(a[1][i] - b[1][0]) / abs(b[1][0]))
        record(label, "pass" if worst <= 1e-6 else "FAIL", f"{len(p0)} points x 3 sources, worst relative difference {worst:.2e} (bound 1e-6)")
        # adjoint=1
        o1 = run(pxf_deck(d, d + "_pxf1.txt", 1), d + "_pxf1")
        p1 = rd(d + "_pxf1.txt")
        w1 = max(abs(a[1][i] - b[1][i]) / abs(b[1][i]) for a, b in zip(p1, p0) for i in range(3))
        refused = "adjoint engine REFUSED" in o1 and "adjoint done" not in o1
        record(f"pxf adjoint=1 refuses (saltation events) and falls back to the per-source engine: {d}", "pass" if refused and w1 <= 1e-6 else "FAIL",
               f"refusal message {'present' if refused else 'ABSENT'}; fallback vs adjoint=0 {w1:.2e} (bound 1e-6)")
    except (CouldNotRun, OSError, IndexError, ZeroDivisionError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])


def pxf_points(d, env=None, srcs=SOURCES, tag=""):
    run(pxf_deck(d, d + "_pxf9%s.txt" % tag, 0, srcs=srcs), d + "_pxf9" + tag, env=env)
    return {int(round(f)): v for f, v in rd(d + "_pxf9%s.txt" % tag)}


def pert_errors(d, byf, srcs=SOURCES):
    """list of (source, f, dB error, deg error) of the .pxf values against PERT"""
    errs = []
    for i, s in enumerate(srcs):
        for f, (rdb, rdeg) in sorted(PERT[d][s].items()):
            z = byf[f][i]
            if abs(z) == 0.0:
                errs.append((s, f, -999.0, 180.0))
                continue
            db, dg = dbdeg(z)
            errs.append((s, f, db - rdb, dphase(dg, rdeg)))
    return errs


def check_pxf_vs_pert(d):
    label = f"pxf vs the perturbed PSS (3 sources x 20/50/100 kHz) {PXF_DB} dB / {PXF_DEG} deg: {WHAT[d]} [{d}]"
    try:
        errs = pert_errors(d, pxf_points(d))
    except (CouldNotRun, OSError, IndexError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:300])
        return
    wdb = max(abs(e[2]) for e in errs)
    wdg = max(abs(e[3]) for e in errs)
    per = "; ".join("%s %+.3f/%+.2f" % (s, max((e[2] for e in errs if e[0] == s), key=abs), max((e[3] for e in errs if e[0] == s), key=abs)) for s in SOURCES)
    record(label, "pass" if wdb <= PXF_DB and wdg <= PXF_DEG else "FAIL", f"worst {wdb:.3f} dB / {wdg:.3f} deg (per source, dB/deg: {per})")


# check 4: the dead-time half bridges (see docstring); card #709 pinned the +0.04 dB residual here as a window, card #727 removed it
DEADTIME_WINDOW = {"syn_dt": (-0.005, 0.005), "syn_dt_q": (-0.005, 0.005)}


def check_deadtime_residual():
    try:
        for d, (lo, hi) in DEADTIME_WINDOW.items():
            errs = pert_errors(d, pxf_points(d), srcs=SOURCES)
            c = [e[2] for e in errs if e[0] == "Vctrl"]
            others = max(abs(e[2]) for e in errs if e[0] != "Vctrl")
            ok = all(lo <= x <= hi for x in c) and others <= 0.005
            record(f"dead-time half bridge control gain vs the perturbed PSS, window [{lo}, {hi}] dB (card #727; was +0.04 dB) [{d}]",
                   "pass" if ok else "FAIL",
                   "Vctrl error %s dB at 20/50/100 kHz; Vin/Iinj worst %.4f dB (<= 0.005)%s" % (", ".join("%+.4f" % x for x in c), others,
                   "" if ok else "  (outside the window: the card #727 offset is back, or worse; re-read docs/pss-pac.md)"))
        # no dead time: the same deck without the offset
        d = "syn0"
        errs = pert_errors(d, pxf_points(d))
        c = max(abs(e[2]) for e in errs if e[0] == "Vctrl")
        record("same circuit without dead time has no such offset [syn0]", "pass" if c <= 0.005 else "FAIL", f"Vctrl worst {c:.4f} dB (<= 0.005)")
    except (CouldNotRun, OSError, IndexError) as ex:
        record("dead-time half bridge", "COULD-NOT-RUN", str(ex)[:300])


# ------------------------------------------------------------------------------------------------------------------------------
# closed loop: boost, ideal low-side sw, type-I compensator (tests/fixtures/pxfpstb709/cl_boost.cir)
# ------------------------------------------------------------------------------------------------------------------------------
CL_FREQS = (1000, 1500, 2000, 5000, 10000, 20000, 50000, 100000, 125000)
# FRA: tools/fra.py run_loop_point(deck lines of cl_boost.cir + `.save v(xp.x) v(fbtop) i(v.xp.vi) i(v.xp.vb)`, "XP", f, A, settle_time=1.5e-3,
# settle_periods=spp, measure_periods=mpp, points_per_period=256, startup="uic", tmax=tmax, spice_options="options trtol=1 reltol=1e-4\nset
# num_threads=2\n"), NGSPICE = the s141 binary; Tian combined estimate (Middlebrook agrees to 1e-8 dB: no reverse transmission).  Rows: f: (dB, deg) of T
# [A, tmax, spp/mpp].  Logs: ~/Work/data/ngspice_augmentation/709/fra_fb.jsonl, fra_cl_a.jsonl.
#   1000 [0.2, 5e-10, 5/5]  1500 [0.2, 5e-10, 5/5]  2000 [0.2, 5e-10, 5/5]  5000 [0.2, 2.5e-10, 5/5]  10000 [0.2, 2.5e-10, 20/10]
#   20000 [0.2, 2.5e-10, 20/10]  50000 [0.2, 2.5e-10, 20/10]  100000 [0.2, 2.5e-10, 20/10]
# Halving checks (same deck): 2000: A 0.1 gives (-3.001, -96.81); 5000: A 0.1 gives (-5.005, -122.85), tmax 5e-10 gives (-4.675, -122.54): the FRA's own spread
# at 5 kHz is 0.33 dB / 0.3 deg.  10 kHz: tmax 5e-10 / 1e-9 give (-20.25 .. -20.07 dB, 102.7 .. 102.3 deg) with 3 settle periods: the injection transient has
# not died (spp 20 needed, hence the 20/10 setting: at spp 3 the phase is 4 deg off).
# Above 20 kHz the FRA is useless on this loop: T is -65 dB (50 kHz) and the FRA returns -58.6 dB / 19.7 deg, -59.0 / 23.8 with the loop broken at the
# comparator input instead, for A = 0.05 .. 0.3 and tmax 1.25e-10 .. 5e-10 alike (an additive error of about -59 dB that does not depend on the amplitude or
# the step).  The perturbed PSS (PSTB_LG) gives -64.94 dB / 56.1 deg there; they are kept as information, not as a bar.
FRA = {1000: (2.378, -93.22), 1500: (-0.914, -94.98), 2000: (-3.014, -96.85), 5000: (-4.695, -122.36), 10000: (-20.750, 98.02),
       20000: (-41.747, 78.46), 50000: (-58.590, 19.65), 100000: (-67.635, 4.98)}
FRA_BAR = {1000: (0.25, 2.5), 1500: (0.25, 2.5), 2000: (0.25, 2.5), 5000: (0.25, 2.5), 10000: (0.30, 2.5), 20000: (0.50, 2.5)}   # 20 kHz: FRA floor
# perturbed-PSS loop gain (tools/pac_perturbed.py sweep --deck <cl_boost without the XP line> --loop-gain --probe-nodes fbtop:vout --inject
# node:opout_pac_unused --freqs f --backend pss --stabtime 2e-3 --amplitude A [--halving-bound-db 3 --halving-bound-deg 20]); (dB, deg) of T, [A, halving]
PSTB_LG = {7812.5: ((-10.6731, 118.2290), (0.05, 0.2)), 10000: ((-20.9842, 97.8105), (0.2, 0.2)), 20000: ((-42.1193, 78.7174), (1.0, 0.2)),
           50000: ((-64.9395, 56.1154), (1.0, 0.2)), 100000: ((-79.7999, 40.0575), (1.0, 0.2)), 125000: ((-83.9763, 39.4029), (1.0, 0.2))}
# bars: 0.05 dB / 0.4 deg up to 50 kHz; the tool's own halving error at 100 / 125 kHz is 0.57 / 0.97 dB (T = -80 / -84 dB): 0.5 dB / 3 deg there
PSTB_LG_BAR = {7812.5: (0.05, 0.4), 10000: (0.05, 0.4), 20000: (0.05, 0.4), 50000: (0.05, 0.4), 100000: (0.5, 3.0), 125000: (0.5, 3.0)}
RA, RB, R1, CF = 6700.0, 1000.0, 10e3, 10e-9


def cl_deck(cf="10n", mode="fb"):
    b = fixture("cl_boost.cir").replace("../../../stb/lgprobe.lib", os.path.join(ROOT, "stb", "lgprobe.lib")).replace("cf=10n", "cf=" + cf)
    if mode == "comp":   # loop broken at the comparator input instead of at the feedback divider
        b = (b.replace("S1 sw 0 opout ramp SWMOD", "S1 sw 0 opc ramp SWMOD").replace("XP fbtop vout lgprobe", "Vfb fbtop vout DC 0\nXP opc opout lgprobe")
             .replace("v(opout)=0.8", "v(opout)=0.8 v(opc)=0.8"))
    return b


def pstb_at(f, env=None, tag="", mode="fb"):
    name = "cl_%s_%g%s" % (mode, f, tag)
    d = (cl_deck(mode=mode) + ".options reltol=1e-4\n.pstb 500e3 2e-3 vout 200 10 150 1e-6 lin 1 %g %g method=\"shoot\" probe=\"XP\"\n" % (f, f)
         + WR % (name + ".txt", "t"))
    out = run(d, name, env=env)
    return rd(name + ".txt")[0][1][0], out


def pstb_margin():
    """phase margin / crossover as printed by .pstb on a lin 11 1000..2000 Hz sweep"""
    import re
    d = (cl_deck() + ".options reltol=1e-4\n.pstb 500e3 2e-3 vout 200 10 150 1e-6 lin 11 1000 2000 method=\"shoot\" probe=\"XP\"\n" + WR % ("cl_margin.txt", "t"))
    out = run(d, "cl_margin")
    m = re.search(r"phase margin ([-0-9.]+) deg at ([0-9.]+) Hz", out)
    if not m:
        raise CouldNotRun("no phase margin line: " + out[-300:])
    return float(m.group(1)), float(m.group(2))


def analytic_T(f, gvc):
    beta = (RB * R1 / (RB + R1)) / (RA + RB * R1 / (RB + R1))
    return gvc * beta / (2j * math.pi * f * R1 * CF)


def fra_crossover():
    """crossover frequency and phase margin from the pinned FRA at 1 and 1.5 kHz (linear in dB and degrees against ln f)"""
    (d1, p1), (d2, p2) = FRA[1000], FRA[1500]
    x = d1 / (d1 - d2)
    return 1000.0 * math.exp(x * math.log(1.5)), 180.0 + p1 + x * (p2 - p1)


def check_pstb():
    try:
        T = {f: pstb_at(f)[0] for f in CL_FREQS}
        T["7812.5"] = pstb_at(7812.5)[0]
    except (CouldNotRun, OSError, IndexError) as ex:
        record("pstb closed loop", "COULD-NOT-RUN", str(ex)[:300])
        return
    # a. FRA
    parts, bad, info = [], [], []
    for f in CL_FREQS:
        if f not in FRA:
            continue
        db, dg = dbdeg(T[f])
        e1, e2 = db - FRA[f][0], dphase(dg, FRA[f][1])
        if f in FRA_BAR:
            parts.append("%g %+.3f dB %+.2f deg" % (f, e1, e2))
            if abs(e1) > FRA_BAR[f][0] or abs(e2) > FRA_BAR[f][1]:
                bad.append(f)
        else:
            info.append("%g %+.1f dB %+.1f deg" % (f, e1, e2))
    record("closed-loop boost .pstb vs FRA (tools/fra.py loop mode), 1 kHz .. 20 kHz: 0.25 dB / 2.5 deg (0.30 at 10 kHz, 0.50 at 20 kHz: the FRA's own floor)",
           "pass" if not bad else "FAIL", "; ".join(parts) + ("; FRA unusable above 20 kHz (not a bar): " + "; ".join(info)))
    # a2. phase margin and crossover
    try:
        pm, fc = pstb_margin()
        fcf, pmf = fra_crossover()
        record("closed-loop boost phase margin / crossover: .pstb printed vs FRA interpolated at 1 / 1.5 kHz (2.5 deg, 3 %)",
               "pass" if abs(pm - pmf) <= 2.5 and abs(fc / fcf - 1) <= 0.03 else "FAIL",
               ".pstb %.2f deg at %.1f Hz; FRA %.2f deg at %.1f Hz" % (pm, fc, pmf, fcf))
    except (CouldNotRun, OSError) as ex:
        record("closed-loop boost phase margin", "COULD-NOT-RUN", str(ex)[:300])
    # b. perturbed-PSS loop gain
    parts, bad = [], []
    for f, (ref, _) in sorted(PSTB_LG.items()):
        z = T[f] if f in T else T[str(f)]
        db, dg = dbdeg(z)
        e1, e2 = db - ref[0], dphase(dg, ref[1])
        parts.append("%g %+.3f dB %+.2f deg" % (f, e1, e2))
        if abs(e1) > PSTB_LG_BAR[f][0] or abs(e2) > PSTB_LG_BAR[f][1]:
            bad.append(f)
    record("closed-loop boost .pstb vs the perturbed-PSS loop gain (pac_perturbed --loop-gain), 7.8 kHz .. fs/4: 0.05 dB / 0.4 deg to 50 kHz, 0.5 / 3 above (reference halving noise)",
           "pass" if not bad else "FAIL", "; ".join(parts))
    # c. analytic composition from .pxf
    try:
        gv = {}
        for f in CL_FREQS:
            run(pxf_deck("boost", "cl_gvc_%g.txt" % f, 0, sweep="lin 1 %g %g" % (f, f), srcs=("Vctrl",)), "cl_gvc_%g" % f)
            gv[f] = rd("cl_gvc_%g.txt" % f)[0][1][0]
    except (CouldNotRun, OSError, IndexError) as ex:
        record("closed-loop boost .pstb vs the analytic composition", "COULD-NOT-RUN", str(ex)[:300])
        return
    parts, wdb, wdg = [], 0.0, 0.0
    for f in CL_FREQS:
        adb, adg = dbdeg(analytic_T(f, gv[f]))
        db, dg = dbdeg(T[f])
        parts.append("%g %+.3f/%+.2f" % (f, db - adb, dphase(dg, adg)))
        wdb, wdg = max(wdb, abs(db - adb)), max(wdg, abs(dphase(dg, adg)))
    record("closed-loop boost .pstb vs T = Zf/R1 x beta x Gvc (open-loop .pxf control-to-output), 1 kHz .. fs/4: 0.08 dB / 0.3 deg",
           "pass" if wdb <= 0.08 and wdg <= 0.3 else "FAIL", "dB/deg: " + "; ".join(parts) + f"; worst {wdb:.3f} dB / {wdg:.2f} deg")
    # d. probe location: comparator input instead of the divider
    try:
        parts, bad = [], []
        for f, (bdb, bdg) in ((1000, (0.05, 0.3)), (10000, (0.05, 0.3)), (50000, (0.5, 2.5))):
            z = pstb_at(f, mode="comp")[0]
            db, dg = dbdeg(z)
            d1, d2 = dbdeg(T[f])
            parts.append("%g %+.3f dB %+.2f deg" % (f, db - d1, dphase(dg, d2)))
            if abs(db - d1) > bdb or abs(dphase(dg, d2)) > bdg:
                bad.append(f)
        record("closed-loop boost .pstb: loop broken at the comparator input vs at the divider (same loop), 1 / 10 kHz 0.05 dB / 0.3 deg, 50 kHz 0.5 dB / 2.5 deg "
               "(at fs/4 the two differ by 0.86 dB / 8.9 deg: sideband-0 only, informational)", "pass" if not bad else "FAIL", "; ".join(parts))
    except (CouldNotRun, OSError, IndexError) as ex:
        record("closed-loop boost probe location", "COULD-NOT-RUN", str(ex)[:300])


# ------------------------------------------------------------------------------------------------------------------------------
# red proofs: each toggle removes one piece of the saltation treatment; the agreement the checks above demand must break
# ------------------------------------------------------------------------------------------------------------------------------
def red_pxf(d, env, min_db, freq=20000, what=""):
    label = f"RED PROOF {env}: .pxf {d} control-to-output vs the perturbed PSS must be off by more than {min_db} dB at {freq / 1e3:g} kHz {what}"
    try:
        byf = pxf_points(d, env=env, srcs=("Vctrl",), tag="_red")
        z = byf[freq][0]
        ref = PERT[d]["Vctrl"][freq][0]
        err = (20 * math.log10(abs(z)) - ref) if abs(z) > 0 else -999.0
    except (CouldNotRun, OSError, IndexError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:200])
        return
    record(label, "pass" if abs(err) > min_db else "FAIL", f"error {err:+.3f} dB")


def red_pstb(env, min_db, f, ref, what=""):
    label = f"RED PROOF {env}: .pstb closed-loop boost at {f / 1e3:g} kHz must be off by more than {min_db} dB {what}"
    try:
        z = pstb_at(f, env=env, tag="_red")[0]
        err = (20 * math.log10(abs(z)) - ref) if abs(z) > 0 else -999.0
    except (CouldNotRun, OSError, IndexError) as ex:
        record(label, "COULD-NOT-RUN", str(ex)[:200])
        return
    record(label, "pass" if abs(err) > min_db else "FAIL", f"error {err:+.3f} dB")


def red_proofs():
    red_pxf("boost", {"PSS_SALTATION": "0"}, 100.0, what="(no saltation at all: the control gain is 0)")
    red_pxf("buck", {"PSS_SALTATION": "0"}, 100.0)
    red_pxf("syn_dt", {"PSS_SALTATION": "0"}, 100.0)
    red_pxf("boost", {"PSS_SALT_NO579": "1"}, 10.0, what="(the pre-event node equation without the anode-side diode: 71 dB)")
    red_pxf("boost", {"PSS_SALT_NODIOFAR": "1"}, 0.15, what="(no jump on the capacitor the diode feeds)")
    red_pxf("boost", {"PAC_NO_IMPULSE": "1"}, 0.5, freq=100000, what="(no output impulse)")
    red_pxf("syn0", {"PSS_SALT_NO582": "1"}, 10.0, what="(one routed jump per switch, switches flipping in one step)")
    red_pxf("syn_dt_q", {"PSS_SALT_NOCONT": "1"}, 10.0, what="(single-step routed jump)")
    red_pxf("boost_q", {"PSS_SALT_NOCONT": "1"}, 5.0, freq=100000)
    red_pxf("syn_dt", {"PSS_SALT_NO727": "1"}, 0.025, what="(the body diode's series resistance left out of the pre-event node equation: +0.04 dB)")
    red_pxf("syn_dt_q", {"PSS_SALT_NO727": "1"}, 0.025, what="(the same, and the continuation cut at 1e-3 of the peak rate: +0.035 dB)")
    # the .pstb: against the pinned FRA at 1 kHz (T = +2.4 dB) and the perturbed-PSS loop gain at 10 kHz / 50 kHz
    red_pstb({"PSS_SALTATION": "0"}, 50.0, 1000, FRA[1000][0], "(no saltation: the loop has no modulator gain)")
    red_pstb({"PSS_SALT_NO579": "1"}, 50.0, 10000, PSTB_LG[10000][0][0])
    red_pstb({"PSS_SALT_NODIOFAR": "1"}, 0.8, 50000, PSTB_LG[50000][0][0], "(1.5 dB and 46 deg)")
    red_pstb({"PAC_NO_IMPULSE": "1"}, 0.2, 50000, PSTB_LG[50000][0][0], "(0.36 dB)")


# ------------------------------------------------------------------------------------------------------------------------------
# optional: recompute one FRA point (about 2 minutes) and compare it with the pinned value
# ------------------------------------------------------------------------------------------------------------------------------
def fra_live():
    try:
        sys.path.insert(0, os.path.join(ROOT, "tools"))
        os.environ["NGSPICE"] = NG
        import fra
        import stb
        lines = stb.prepare_deck(os.path.join(FX, "cl_boost.cir")) + [".save v(xp.x) v(fbtop) i(v.xp.vi) i(v.xp.vb)"]
        r = fra.run_loop_point(lines, "XP", 2000.0, 0.2, settle_time=1.5e-3, settle_periods=5, measure_periods=5, points_per_period=256, startup="uic",
                               tmax=5e-10, spice_options="options trtol=1 reltol=1e-4\nset num_threads=2\n", outdir=os.path.join(OUT, "fra_live"), ngspice=NG)
        db, dg = dbdeg(r["combined"])
        e1, e2 = db - FRA[2000][0], dphase(dg, FRA[2000][1])
        record("FRA live at 2 kHz reproduces the pinned value (0.05 dB / 0.3 deg)", "pass" if abs(e1) <= 0.05 and abs(e2) <= 0.3 else "FAIL",
               "%.3f dB %.2f deg vs pinned %.3f %.2f" % (db, dg, FRA[2000][0], FRA[2000][1]))
    except Exception as ex:   # noqa: BLE001
        record("FRA live at 2 kHz", "COULD-NOT-RUN", str(ex)[:300])


def main():
    if not os.path.exists(NG):
        print("ngspice not found:", NG)
        return 2
    for d in DECKS:
        check_pxf_vs_pac(d)
    for d in DECKS:
        check_pxf_vs_pert(d)
    check_deadtime_residual()
    check_pstb()
    red_proofs()
    if LIVE_FRA:
        fra_live()
    n = len(results)
    npass = sum(1 for r in results if r[1] == "pass")
    print(f"\n{npass}/{n} pass")
    if any(r[1] == "FAIL" for r in results):
        return 1
    if any(r[1] == "COULD-NOT-RUN" for r in results):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
