#!/usr/bin/env python3
"""Validation for card #716: lossless transmission lines (TRA, the T element) in the shooting `.pss method=shoot` and in the
time-domain small-signal analyses (`.pac/.pstb/.pxf solver="td"`).

Before (s145, measured for this card): a TRA restarted every shooting period from the constant history of its initial conditions;
`.pss` reported CONVERGED with an orbit 23 % to 100 % of the peak-peak off the settled transient on RC-terminated lines and up to
80 times the ripple on a buck output line; td `.pac/.pstb/.pxf` refused the circuit (card #542).  Now the line's delay table is
part of the period map (pssshoot.c: pss_tl_latch(), pss_tl_newton(); tra/traload.c: TRApssHist) and of the small-signal
recursion and its Floquet boundary problem (pacshoot.c: pac_tl_*).  LTRA (lossy), TXL and CPL lines stay refused, so does a TRA
in `.pnoise` and in a free-running (trigger=) analysis.

    NGSPICE=<prefix>/bin/ngspice python3 tests/test_tline_pss.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN (and no FAIL)
    outputs: ~/Work/data/ngspice_augmentation/716/test (TLINE_PSS_OUT overrides), never the repository

Checks (decks in tests/fixtures/tline716; references: the orbit checks run their own settled `.tran` with this binary, the
transmission line's transient model being unchanged by this card; the small-signal LTI checks run `.ac`/`.stb` (frequency
domain, the line's exact e^{-jwTD})):
  orbit   1-5  .pss orbit vs a settled .tran (400 periods, reltol 1e-8), pointwise over the period, as a fraction of the reference's
               peak-peak, allowing a timing shift of up to T/1000 (the reference's output step: an edge one output step off is
               not an amplitude error); decks: RC-terminated line TD 0.3 T, TD 1.5 T (longer than the period), a strongly
               mismatched RC line, an LC filter feeding a 5 ohm line into 3.3 ohm || 100 nF (high-Q line), an ideal-switch buck
               with that line in its output path.  Also: the joint Newton converged.  Bounds: ORBIT below.
  red     6    PSS_TL_OFF=1 (constant history, the s145 behaviour): the RC orbit > 10 % of the peak-peak off.
  pac     7-10 td .pac == .ac on three linear decks (a source driving a line into an LTI load: sideband 0 of an LTI circuit is its
               .ac response), 3 frequencies each, and on one of them through the generic (charge) recursion.
  pxf     11-12 td .pxf (forward engine) == .ac, TD 0.3 T and 1.5 T.
  pstb    13   .pstb == native .stb on a transconductance loop whose output reaches the divider through a line, 1 kHz - 1 MHz
               (above that this deck and the same deck without the line have the same O(h^2) time-step error of the engine).
  sw      14-15 td .pac on the switching decks vs pinned perturbed-PSS references (tests/fixtures/tline716/refs.txt: tools/
               pac_perturbed.py sweep --backend pss --no-step-check, card #716 build), 20 / 50 / 100 kHz.
  refuse  16-19 TXL and LTRA (lossy RLC) in .pss: REFUSED; TRA in .pnoise: REFUSED; TRA in a free-running .pss: REFUSED.
  red     20-21 PAC_TL_NO_DELAY=1 (the small-signal waves absorbed by the line) puts td .pac > 1 dB off .ac; PSS_TL_OFF=1 brings back
               the td .pac refusal.
"""
import cmath
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures", "tline716")
NGSPICE = os.environ.get("NGSPICE", os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s152/bin/ngspice"))
DATA = os.environ.get("TLINE_PSS_OUT", os.path.expanduser("~/Work/data/ngspice_augmentation/716/test"))
T = 2e-6   # every deck but tl_pstb switches at 500 kHz

# ---- bounds (measured on the card #716 build, see the card report and docs/pss-pac.md; roughly 2-3x the measured value) ----
ORBIT = {   # deck tag -> (bound as fraction of peak-peak, timing-tolerant)
    "rc": 0.008, "rc_td3": 0.025, "rc_hiq": 0.006, "lcp": 0.002, "buckp": 0.02,
}
# measured (w716, Points 400): rc 0.37 %, rc_td3 1.22 % (a history corner: the reflected source edge arrives with the next source
# edge; the sampled history smears it over ~2 samples of T/Points, O(dt^2): 3.4 / 1.2 / 0.24 % at Points 200 / 400 / 800), rc_hiq
# 0.28 %, lcp 0.026 %, buckp 0.80 % (the switching deck's transient reference itself moves by 1.7 % of the ripple between reltol
# 1e-8 at T/1000 and T/10000 maximum step: the bound is that spread)
PAC_DB, PAC_DEG = 2e-3, 0.02          # LTI td .pac vs .ac (measured <= 1.3e-4 dB / 3.2e-4 deg)
PXF_REL = 2e-5                        # measured <= 1.9e-6
PSTB_REL = 3e-4                       # up to 1 MHz; measured 6.3e-5 (the same deck without the line: 1.1e-4)
SW_DB, SW_DEG = 0.01, 0.1             # switching decks vs perturbed PSS: measured <= 0.001 dB / 0.016 deg; the references
                                      # amplitude halving <= 0.0008 dB / 0.012 deg (refs.txt)

# pinned references: tools/pac_perturbed.py sweep --backend pss (card #716 build), ~/Work/data/ngspice_augmentation/716/ref
REF_BUCKP = {}   # f -> (dB, deg) of v(vl)/v(Vin), filled below
REF_BUCKC = {}   # f -> (dB, deg) of v(vl)/v(Vctrl)

results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"{state:14s} {name}\n{'':14s}   {detail}", flush=True)


def body(deck, params=None):
    t = open(os.path.join(FIX, deck)).read()
    for k, v in (params or {}).items():
        t, n = re.subn(r"(?m)(^\.param\b.*?\b%s=)(\S+)" % re.escape(k), lambda m: m.group(1) + str(v), t)
        if not n:
            t += f".param {k}={v}\n"
    return t.replace("\n.end\n", "\n").rstrip("\n") + "\n"


def run(text, tag, env_extra=None, timeout=1800):
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, tag + ".cir")
    open(path, "w").write(text)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "2"
    for k in ("PSS_TL_OFF", "PAC_TL_NO_DELAY"):
        env.pop(k, None)
    if env_extra:
        env.update(env_extra)
    try:
        p = subprocess.run([NGSPICE, "-b", path], capture_output=True, text=True, timeout=timeout, env=env)
    except (OSError, subprocess.TimeoutExpired) as ex:
        raise CouldNotRun(f"{tag}: {ex}")
    log = p.stdout + "\n" + p.stderr
    open(os.path.join(DATA, tag + ".log"), "w").write(log)
    return log


def rows(path):
    out = []
    if not os.path.exists(path):
        return out
    for l in open(path):
        try:
            out.append([float(x) for x in l.split()])
        except ValueError:
            pass
    return out


def interp(r, col, t):
    lo, hi = 0, len(r) - 1
    if t <= r[0][0]:
        return r[0][col]
    if t >= r[-1][0]:
        return r[-1][col]
    while hi - lo > 1:
        m = (lo + hi) // 2
        if r[m][0] <= t:
            lo = m
        else:
            hi = m
    a, b = r[lo], r[hi]
    return a[col] if b[0] == a[0] else a[col] + (b[col] - a[col]) * (t - a[0]) / (b[0] - a[0])


def orbit_error(orb, ref, col, n=2000):
    """(timing-tolerant max error, plain max error, peak-peak of the reference) over one period"""
    d = T / 1000
    vals = [interp(ref, col, T * k / n) for k in range(n)]
    pp = max(vals) - min(vals)
    w1 = w2 = 0.0
    for k in range(n):
        t = T * k / n
        o = interp(orb, col, t)
        w1 = max(w1, abs(o - vals[k]))
        win = [interp(ref, col, (t + s * d / 4) % T) for s in range(-4, 5)]
        lo, hi = min(win), max(win)
        w2 = max(w2, 0.0 if lo <= o <= hi else min(abs(o - lo), abs(o - hi)))
    return w2, w1, pp


def tran_ref(deck, params, nodes, tag):
    out = os.path.join(DATA, tag + ".txt")
    if os.path.exists(out):
        os.remove(out)
    run(body(deck, params) + ".options reltol=1e-8 abstol=1e-12 vntol=1e-9\n.control\nset wr_singlescale\nset numdgt=15\n"
        f"tran {T / 1000:g} {400 * T:g} {399 * T:g} {T / 1000:g}\nwrdata {out} " + " ".join(f"v({x})" for x in nodes)
        + "\nquit\n.endc\n.end\n", tag)
    r = rows(out)
    if not r:
        raise CouldNotRun(f"no .tran output for {tag}")
    t0 = r[0][0]
    return [[x[0] - t0] + x[1:] for x in r]


def pss_orbit(deck, params, nodes, tag, env=None):
    out = os.path.join(DATA, tag + ".txt")
    if os.path.exists(out):
        os.remove(out)
    log = run(body(deck, params) + ".options reltol=1e-5\n.control\nset num_threads=2\nset wr_singlescale\nset numdgt=15\n"
              f"pss 500e3 100e-6 {nodes[0]} 400 10 150 1e-6 method=shoot\nwrdata {out} " + " ".join(f"v({x})" for x in nodes)
              + "\nquit\n.endc\n.end\n", tag, env)
    return rows(out), log


def check_orbit(name, key, deck, params, nodes, env=None, red=False):
    t0 = time.time()
    try:
        ref = tran_ref(deck, params, nodes, key + "_tran")
        orb, log = pss_orbit(deck, params, nodes, key + ("_red" if red else "_pss"), env)
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    if not orb:
        record(name, "COULD-NOT-RUN" if not red else "FAIL", "no orbit: " + " | ".join(
            l for l in log.split("\n") if "FAILED" in l or "REFUS" in l)[:300])
        return
    worst, parts = 0.0, []
    for c, nd in enumerate(nodes, 1):
        w2, w1, pp = orbit_error(orb, ref, c)
        worst = max(worst, w2 / pp)
        parts.append(f"v({nd}) {100 * w2 / pp:.3f} % (plain {100 * w1 / pp:.3f} %) of p-p {pp:.4g} V")
    jn = [l for l in log.split("\n") if "joint Newton" in l]
    if red:
        ok = worst > 0.10
        record(name, "pass" if ok else "FAIL", "; ".join(parts) + f"; worst {100 * worst:.2f} % (must exceed 10 %)")
        return
    conv = bool(jn) and "converged" in jn[0] and "NOT converged" not in jn[0] and "CONVERGED" in log
    ok = conv and worst <= ORBIT[key]
    m = re.search(r"of (\d+) line\(s\): [^,]*, residual ([^ ]+) on its step grid; a grid recorded at the solution moves it by "
                  r"([^ ]+) \(native, of the ripple\) / ([^ ]+) ", jn[0]) if jn else None
    extra = f"; grid check {m.group(3)} / {m.group(4)}" if m else ""
    record(name, "pass" if ok else "FAIL", "; ".join(parts) + f"; worst {100 * worst:.3f} % (bound {100 * ORBIT[key]:.2f} %)"
           + extra + ("" if conv else "; joint Newton NOT converged") + f"; wall {time.time() - t0:.0f}s")


def pac_point(deck, params, node, f, tag, env=None):
    out = os.path.join(DATA, tag + ".txt")
    if os.path.exists(out):
        os.remove(out)
    log = run(body(deck, params) + ".options reltol=1e-5\n"
              f'.pac 500e3 100e-6 {node} 400 10 150 1e-6 lin 1 {f} {f} method="shoot" solver="td"\n'
              ".control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n"
              f"wrdata {out} v({node})\nquit\n.endc\n.end\n", tag, env)
    if "REFUS" in log:
        raise CouldNotRun("REFUSED: " + [l for l in log.split("\n") if "REFUS" in l][0][:200])
    r = rows(out)
    if not r or len(r[0]) < 3:
        raise CouldNotRun(f"no .pac output ({tag})")
    return complex(r[0][1], r[0][2])


def ac_point(deck, params, node, f, tag):
    out = os.path.join(DATA, tag + ".txt")
    run(body(deck, params) + f".control\nset wr_singlescale\nset numdgt=15\nac lin 1 {f} {f}\nwrdata {out} v({node})\nquit\n"
        ".endc\n.end\n", tag)
    r = rows(out)
    if not r:
        raise CouldNotRun(f"no .ac output ({tag})")
    return complex(r[0][1], r[0][2])


def db(z):
    return 20 * math.log10(abs(z)) if abs(z) > 0 else -400.0


def dph(a, b):
    return (math.degrees(cmath.phase(a)) - math.degrees(cmath.phase(b)) + 180.0) % 360.0 - 180.0


def check_pac_lti(name, key, deck, params, node, freqs, env=None):
    wdb = wdeg = 0.0
    parts = []
    try:
        for f in freqs:
            z = pac_point(deck, params, node, f, f"{key}_pac_{int(f)}", env)
            a = ac_point(deck, params, node, f, f"{key}_ac_{int(f)}")
            wdb, wdeg = max(wdb, abs(db(z) - db(a))), max(wdeg, abs(dph(z, a)))
            parts.append(f"{f / 1e3:g}k {db(z):.4f}/{math.degrees(cmath.phase(z)):.3f} vs .ac {db(a):.4f}/{math.degrees(cmath.phase(a)):.3f}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    record(name, "pass" if wdb <= PAC_DB and wdeg <= PAC_DEG else "FAIL",
           "; ".join(parts) + f"; worst {wdb:.2e} dB / {wdeg:.2e} deg (bound {PAC_DB} / {PAC_DEG})")


def check_pxf(name, key, deck, params):
    out = os.path.join(DATA, key + "_pxf.txt")
    if os.path.exists(out):
        os.remove(out)
    try:
        log = run(body(deck, params) + '.options reltol=1e-5\n.pxf 500e3 100e-6 b 400 10 150 1e-6 b dec 2 30k 300k method="shoot" '
                  'solver="td" sources="Vs"\n.control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n'
                  f"wrdata {out} pxf_vs\nquit\n.endc\n.end\n", key + "_pxf")
        r = rows(out)
        if not r:
            raise CouldNotRun("no .pxf output: " + " | ".join(l for l in log.split("\n") if "REFUS" in l or "rror" in l)[:200])
        worst, parts = 0.0, []
        for x in r:
            z = complex(x[1], x[2])
            a = ac_point(deck, params, "b", x[0], f"{key}_pxfac_{int(x[0])}")
            e = abs(z - a) / abs(a)
            worst = max(worst, e)
            parts.append(f"{x[0] / 1e3:.4g}k rel {e:.1e}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    record(name, "pass" if worst <= PXF_REL else "FAIL", "; ".join(parts) + f"; worst {worst:.2e} (bound {PXF_REL})")


def check_pstb(name):
    stb = os.path.join(DATA, "pstb_stb.txt")
    pst = os.path.join(DATA, "pstb_pstb.txt")
    for f in (stb, pst):
        if os.path.exists(f):
            os.remove(f)
    b = body("tl_pstb.cir")
    try:
        run(b + '.stb dec 5 1k 1Meg probe="Vprobe"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n'
            f"wrdata {stb} stb1.loopgain\nquit\n.endc\n.end\n", "pstb_stb")
        log = run(b + '.pstb 100k 0 clkgnd 2000 1 5 0.001 dec 5 1k 1Meg probe="Vprobe" method="shoot"\n.control\n'
                  f"set wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata {pst} loopgain\nquit\n.endc\n.end\n", "pstb_pstb")
        a, p = rows(stb), rows(pst)
        if not a or not p or len(a) != len(p):
            raise CouldNotRun("no .stb/.pstb output: " + " | ".join(l for l in log.split("\n") if "REFUS" in l or "FAILED" in l)[:200])
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    worst, wf = 0.0, 0.0
    for x, y in zip(a, p):
        za, zp = complex(x[1], x[2]), complex(y[1], y[2])
        e = abs(zp - za) / abs(za)
        if e > worst:
            worst, wf = e, x[0]
    record(name, "pass" if worst <= PSTB_REL else "FAIL",
           f"{len(a)} points 1k..1Meg, worst relative loop-gain difference {worst:.2e} at {wf:.4g} Hz (bound {PSTB_REL})")


def check_switching(name, key, deck, node, ref):
    if not ref:
        record(name, "COULD-NOT-RUN", "no pinned reference")
        return
    t0 = time.time()
    wdb = wdeg = 0.0
    parts = []
    try:
        for f, (rdb, rdeg) in sorted(ref.items()):
            z = pac_point(deck, {}, node, f, f"{key}_pac_{int(f)}")
            d1 = db(z) - rdb
            d2 = (math.degrees(cmath.phase(z)) - rdeg + 180.0) % 360.0 - 180.0
            wdb, wdeg = max(wdb, abs(d1)), max(wdeg, abs(d2))
            parts.append(f"{f / 1e3:g}k {db(z):.4f}/{math.degrees(cmath.phase(z)):.3f} vs {rdb:.4f}/{rdeg:.3f}")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
        return
    record(name, "pass" if wdb <= SW_DB and wdeg <= SW_DEG else "FAIL",
           "; ".join(parts) + f"; worst {wdb:.4f} dB / {wdeg:.3f} deg (bound {SW_DB} / {SW_DEG}); wall {time.time() - t0:.0f}s")


LTRA_DECK = """* tline716 refusal: lossy line (LTRA, RLC) between a pulse source and an RC load
Vs src 0 PULSE(0 1 0 20n 20n 0.98u 2u) AC 1
Rs src a 20
O1 a 0 b 0 lline
.model lline ltra r=0.5 l=250n g=0 c=100p len=60
Rl b 0 100
Cl b 0 1n
"""


TXL_DECK = """* tline716 refusal: lossy line (TXL, the Y element) between a pulse source and an RC load
Vs src 0 PULSE(0 1 0 20n 20n 0.98u 2u) AC 1
Rs src a 20
Y1 a 0 b 0 ymod LEN=1
.model ymod txl R=0.5 L=250n G=0 C=100p length=60
Rl b 0 100
Cl b 0 1n
"""


def check_refusals():
    name = "TXL (Y element) in .pss method=shoot: REFUSED"
    try:
        log = run(TXL_DECK + ".control\npss 500e3 100e-6 b 400 10 150 1e-6 method=shoot\nquit\n.endc\n.end\n", "refuse_txl")
        ok = "REFUSED" in log and "TransLine" in log
        record(name, "pass" if ok else "FAIL", ([l for l in log.split("\n") if "REFUSED" in l] or ["no REFUSED line"])[0][:220])
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
    name = "LTRA (lossy RLC line) in .pss method=shoot: REFUSED"
    try:
        log = run(LTRA_DECK + ".control\npss 500e3 100e-6 b 400 10 150 1e-6 method=shoot\nquit\n.endc\n.end\n", "refuse_ltra")
        ok = "REFUSED" in log and "LTRA" in log
        record(name, "pass" if ok else "FAIL", ([l for l in log.split("\n") if "REFUSED" in l] or ["no REFUSED line"])[0][:220])
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
    name = "TRA in .pnoise solver=td: REFUSED"
    try:
        log = run(body("tra_rc.cir", {"td": "0.6u"}) + '.pnoise 500e3 100e-6 b 400 10 150 1e-6 b Vs dec 2 30k 300k method="shoot" '
                  'solver="td"\n.control\nrun\nquit\n.endc\n.end\n', "refuse_pnoise")
        ok = "REFUSED" in log and "pnoise" in log
        record(name, "pass" if ok else "FAIL", ([l for l in log.split("\n") if "REFUSED" in l] or ["no REFUSED line"])[0][:220])
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
    name = "TRA in a free-running (trigger=) .pss: REFUSED"
    try:
        log = run(body("tra_rc.cir", {"td": "0.6u"}) + ".control\npss 500e3 100e-6 b 400 10 150 1e-6 method=shoot trigger=a "
                  "trigval=0.5 trigdir=rise\nquit\n.endc\n.end\n", "refuse_freerun")
        ok = "REFUSED" in log and "free-running" in log
        record(name, "pass" if ok else "FAIL", ([l for l in log.split("\n") if "REFUSED" in l] or ["no REFUSED line"])[0][:220])
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))


def check_red_pac():
    name = "[RED PROOF] PAC_TL_NO_DELAY=1: td .pac > 1 dB off .ac (tra_rc, 100 kHz)"
    try:
        z = pac_point("tra_rc.cir", {"td": "0.6u"}, "b", 100e3, "red_nodelay", {"PAC_TL_NO_DELAY": "1"})
        a = ac_point("tra_rc.cir", {"td": "0.6u"}, "b", 100e3, "red_nodelay_ac")
        e = abs(db(z) - db(a))
        record(name, "pass" if e > 1.0 else "FAIL", f"{db(z):.3f} dB vs .ac {db(a):.3f} dB: {e:.2f} dB off")
    except CouldNotRun as ex:
        record(name, "COULD-NOT-RUN", str(ex))
    name = "[RED PROOF] PSS_TL_OFF=1: td .pac refused again (the s145 behaviour)"
    try:
        pac_point("tra_rc.cir", {"td": "0.6u"}, "b", 100e3, "red_tloff_pac", {"PSS_TL_OFF": "1"})
        record(name, "FAIL", "td .pac ran with PSS_TL_OFF=1")
    except CouldNotRun as ex:
        record(name, "pass" if "REFUSED" in str(ex) else "COULD-NOT-RUN", str(ex)[:220])


def main():
    if not os.path.exists(NGSPICE):
        print(f"COULD-NOT-RUN: no binary {NGSPICE} (set NGSPICE)")
        return 2
    rf = os.path.join(FIX, "refs.txt")   # "deck f dB deg" per line: the pinned perturbed-PSS references
    if os.path.exists(rf):
        for l in open(rf):
            t = l.split()
            if len(t) == 4 and not l.startswith("#"):
                (REF_BUCKP if t[0] == "tra_buckp" else REF_BUCKC)[float(t[1])] = (float(t[2]), float(t[3]))
    check_orbit("orbit: RC-terminated line, TD = 0.3 T (tra_rc)", "rc", "tra_rc.cir", {"td": "0.6u"}, ["b", "a"])
    check_orbit("orbit: RC-terminated line, TD = 1.5 T, longer than the period", "rc_td3", "tra_rc.cir", {"td": "3u"}, ["b", "a"])
    check_orbit("orbit: strongly mismatched RC line (Rs 5, Z0 50, 1k || 1n)", "rc_hiq", "tra_rc.cir",
                {"td": "0.6u", "rs": 5, "rl": "1k"}, ["b", "a"])
    check_orbit("orbit: LC filter -> 5 ohm line -> 3.3 ohm || 100 nF (high-Q line, tra_lcp)", "lcp", "tra_lcp.cir", {}, ["vl", "vout"])
    check_orbit("orbit: ideal-switch buck with the line in its output path (tra_buckp)", "buckp", "tra_buckp.cir", {}, ["vl", "vout"])
    check_orbit("[RED PROOF] PSS_TL_OFF=1: RC orbit > 10 % of the peak-peak off", "rc", "tra_rc.cir", {"td": "0.6u"}, ["b", "a"],
                env={"PSS_TL_OFF": "1"}, red=True)
    check_pac_lti("td .pac == .ac: tra_rc TD 0.3 T", "rc", "tra_rc.cir", {"td": "0.6u"}, "b", [20e3, 100e3, 300e3])
    check_pac_lti("td .pac == .ac: tra_rc TD 1.5 T", "rc_td3", "tra_rc.cir", {"td": "3u"}, "b", [20e3, 100e3, 300e3])
    check_pac_lti("td .pac == .ac: tra_lcp (high-Q line)", "lcp", "tra_lcp.cir", {}, "vl", [20e3, 100e3, 300e3])
    check_pac_lti("td .pac == .ac: tra_rc TD 0.3 T on the generic (charge) recursion (PAC_FORCE_GENERIC=1)", "rc_gen", "tra_rc.cir",
                  {"td": "0.6u"}, "b", [20e3, 300e3], env={"PAC_FORCE_GENERIC": "1"})
    check_pxf("td .pxf (forward) == .ac: tra_rc TD 0.3 T", "rc", "tra_rc.cir", {"td": "0.6u"})
    check_pxf("td .pxf (forward) == .ac: tra_rc TD 1.5 T", "rc_td3", "tra_rc.cir", {"td": "3u"})
    check_pstb(".pstb == .stb: loop with a line (tl_pstb), 1 kHz..1 MHz")
    check_switching("switching: tra_buckp line-to-output td .pac vs perturbed PSS", "buckp", "tra_buckp.cir", "vl", REF_BUCKP)
    check_switching("switching: tra_buckc control-to-output td .pac vs perturbed PSS", "buckc", "tra_buckc.cir", "vl", REF_BUCKC)
    check_refusals()
    check_red_pac()
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
