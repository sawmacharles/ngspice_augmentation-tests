#!/usr/bin/env python3
"""38_pnoise_opoint.py -- handbook example 38 (card #552): `.pnoise` evaluates every device's noise at the operating point of EACH STEP
(card #544, patch 0097; before it the diode / BJT / MOS3 / BSIM4 noise densities were read from the state of the END of the period).

    python3 38_pnoise_opoint.py NGSPICE OUTDIR

Two checks, both against plain `.noise` (the tool every analog engineer trusts), small decks of tests/test_pac_charge.py:

  1. LTI bias (a slow clock on an unrelated net gives `.pss` a period; the circuit under test sits at a constant bias): `.pnoise dev=1` must
     equal plain `.noise ... 1` (total and every per-device contribution) -- worst relative difference printed (PNLTI lines).
  2. TIME-VARYING bias: a 100 kHz sine LO on the gate / anode of a diode, BJT or MOS1 stage whose time constant is 3-4 orders below the LO
     period, so that the PSD at 1 kHz is the time average of the instantaneous PSD, which plain `.noise` gives at the DC bias of N = 96 equally
     spaced LO phases (quasi-static brute force, error O(tau/T) ~ 1e-4). `.pnoise` (default) against that, and the same `.pnoise` with
     PAC_NOISE_STALE_STATE=1 (the pre-#544 operating point: the density read from the state of the end of the period) -- PNQS lines.
     Card #556: the plain time average is exact for white noise only; a flicker (1/f) source modulated by the 100 kHz LO folds from f - q fs, so the
     reference for the flicker vectors is Sum_q |a_q|^2 S_w(f - q fs) (a_q = DFT over the 96 phases of the per-phase flicker amplitude).
"""
import math
import os
import re
import subprocess
import sys

OPTS = ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n"
CLK = "* pnoise lti deck\nVclk clk 0 PULSE(0 1 0 1n 1n 500u 1m)\nRclk clk clkgnd 100k\nCclk clkgnd 0 10n\n"
LTI = {
    "diode": CLK + ".model dm d is=1e-14 cjo=2p vj=0.7 m=0.5 tt=5n rs=50 kf=1e-12 af=1\nVin in 0 DC 2 AC 1\nRs in a 10k\nD1 a 0 dm\nRl a out 1k\nCl out 0 1n\n",
    "bjt": CLK + ".model qm npn is=1e-15 bf=100 vaf=80 cje=1p cjc=0.5p tf=0.2n tr=5n rb=100 re=2 kf=1e-14 af=1\nVcc vcc 0 5\nVin in 0 DC 1 AC 1\n"
          "Rb in b 200k\nRc vcc out 4k\nQ1 out b 0 0 qm\nCl out 0 1p\n",
}
QS = {   # (deck with %s for the LO source line, DC, amplitude)
    "diode": ("* qs diode\n.model dm d is=1e-14 cjo=0.2p vj=0.7 m=0.5 tt=0.1n rs=10 kf=1e-14 af=1\n%s\nRs g a 1k\nD1 a out dm\nRl out 0 1k\nCl out 0 1p\nVin in 0 DC 0 AC 1\nRin in 0 1k\n", 0.9, 0.5),
    "bjt": ("* qs bjt\n.model qm npn is=1e-15 bf=100 vaf=80 cje=0.2p cjc=0.1p tf=0.1n tr=1n rb=100 re=2 kf=1e-14 af=1\nVcc vcc 0 3\n%s\nRb g b 10k\nRc vcc out 2k\nQ1 out b 0 0 qm\nCl out 0 0.1p\nVin in 0 DC 0 AC 1\nRin in 0 1k\n", 0.75, 0.15),
}
WR = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=12\nrun\n%s\nwrdata %s %s\nquit\n.endc\n.end\n"


def run(ngspice, out, name, deck, env=None):
    p = os.path.join(out, name + ".cir")
    open(p, "w").write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ngspice, "-b", p],
                       cwd=out, env=e, capture_output=True, text=True, timeout=1800)
    return (r.stdout + r.stderr).replace("\r", "\n")


def rd(out, fn, ncol):
    rows = []
    for ln in open(os.path.join(out, fn)).read().splitlines()[1:]:
        t = ln.split()
        if len(t) >= 1 + ncol:
            rows.append([float(x) for x in t[:1 + ncol]])
    return rows


def listvecs(txt):
    return [m.group(1) for m in re.finditer(r"^\s+(onoise[\w.]+)\s*:", txt, re.M) if m.group(1) != "onoise_spectrum"]


def lti(ngspice, out, name):
    deck = LTI[name]
    sweep = "dec 5 1k 100k"
    txt = run(ngspice, out, "nzl_" + name, deck + OPTS + ".noise v(out) Vin %s 1\n.control\nset wr_singlescale\nrun\nsetplot noise1\ndisplay\n.endc\n.end\n" % sweep)
    vecs = listvecs(txt)
    run(ngspice, out, "nzw_" + name, deck + OPTS + ".noise v(out) Vin %s 1\n" % sweep + WR % ("setplot noise1", "nzl_%s.txt" % name, "onoise_spectrum " + " ".join(vecs)))
    pn = deck + OPTS + '.pnoise 1k 0 clkgnd 2000 1 8 1e-6 out Vin %s cyclo method="shoot" solver="td" dev=1 pacpoints=200000\n' % sweep
    run(ngspice, out, "pnl_" + name, pn + WR % ("", "pnl_%s.txt" % name, "onoise_spectrum " + " ".join(vecs)))
    a, b = rd(out, "nzl_%s.txt" % name, 1 + len(vecs)), rd(out, "pnl_%s.txt" % name, 1 + len(vecs))
    worst, wn = 0.0, ""
    floor = 1e-6
    for j, nm in enumerate(["total"] + vecs):
        for x, y in zip(a, b):
            ref, tot = x[1 + j], x[1]
            e = abs(y[1 + j] - ref) / (ref if ref >= floor * tot else tot)
            if e > worst:
                worst, wn = e, nm
    print("PNLTI %s: .pnoise dev=1 vs .noise, %d points 1k..100k, total + %d device contributions: worst relative difference %.2e (at %s); total %.4e V/rtHz at 1 kHz"
          % (name, len(a), len(vecs), worst, wn, a[0][1]), flush=True)


def qs(ngspice, out, name, N=96):
    deck, vdc, vac = QS[name]
    dc = "Vlo g 0 DC %s"
    sine = "Vlo g 0 SIN(%g %g 100k)" % (vdc, vac)
    txt = run(ngspice, out, "qs_list_" + name, deck % (dc % "0.9") + OPTS + ".noise v(out) Vin dec 1 1000 1000 1\n.control\nrun\nsetplot noise1\ndisplay\n.endc\n.end\n")
    vecs = listvecs(txt)
    cols = "onoise_spectrum " + " ".join(vecs)
    acc = None
    phases = []
    for k in range(N):
        vg = vdc + vac * math.sin(2 * math.pi * (k + 0.5) / N)
        run(ngspice, out, "qs_" + name, deck % (dc % ("%.9g" % vg)) + OPTS + ".noise v(out) Vin dec 1 1000 1000 1\n"
            + WR % ("setplot noise1", "qs_%s.txt" % name, cols))
        row = rd(out, "qs_%s.txt" % name, 1 + len(vecs))[-1][1:]
        phases.append(row)
        acc = [x * x for x in row] if acc is None else [a + x * x for a, x in zip(acc, row)]
    ref = [math.sqrt(a / N) for a in acc]
    # card #556: the time average is the WHITE-noise answer.  A flicker source whose density is modulated by the 100 kHz LO keeps S_w(f) only for the
    # q = 0 term of the modulation; its harmonics q sit at S_w(f - q fs):  S = Sum_q |a_q|^2 f/|f - q fs|  (a_q = DFT of the per-phase amplitude).  The
    # reference is the time average for the white part and this sum for every flicker vector (instance totals and the spectrum recomposed).
    names0 = ["total"] + vecs
    fol = {}
    for j, nm in enumerate(names0):
        if "1overf" in nm and "1overfsw" not in nm:
            a = [ph[j] for ph in phases]
            tot = 0.0
            for q in range(-(N // 2) + 1, N // 2 + 1):
                aq = sum(a[i] * complex(math.cos(2 * math.pi * q * i / N), -math.sin(2 * math.pi * q * i / N)) for i in range(N)) / N
                tot += abs(aq) ** 2 * (1.0 if q == 0 else 1000.0 / abs(1000.0 - q * 100e3))
            fol[nm] = math.sqrt(tot)
    qs_white = list(ref)
    fl_of = {}
    for nm in fol:
        fl_of.setdefault(nm[:nm.rindex("1overf") - 1], []).append(nm)
    for j, nm in enumerate(names0):
        parts = list(fol) if nm == "total" else fl_of.get(nm, [])
        if parts:
            white = qs_white[j] ** 2 - sum(qs_white[names0.index(x)] ** 2 for x in parts)
            ref[j] = math.sqrt(max(white, 0.0) + sum(fol[x] ** 2 for x in parts))
        elif nm in fol:
            ref[j] = fol[nm]
    pn = (deck % sine + OPTS + '.pnoise 100k 50u out 400 1 20 1e-6 out Vin lin 1 1000 1000 cyclo method="shoot" solver="td" dev=1 pacpoints=8000\n'
          + WR % ("", "qsp_%s.txt" % name, cols))
    res = {}
    for tag, env in (("new", None), ("old", {"PAC_NOISE_STALE_STATE": "1"})):
        run(ngspice, out, "qsp_" + name, pn, env=env)
        res[tag] = rd(out, "qsp_%s.txt" % name, 1 + len(vecs))[-1][1:]
    s96 = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(ngspice))), "aaba50ba3f49-s96", "bin", "ngspice")
    if name == "diode" and os.path.exists(s96):      # the real pre-#544 build (series s96), no switch
        run(s96, out, "qsp96_" + name, pn.replace("qsp_", "qsp96_"))
        old96 = rd(out, "qsp96_%s.txt" % name, 1 + len(vecs))[-1][1:]
        print("PNQS %s on the s96 build (before patch 0097): total %.5e (x%.2f of quasi-static)" % (name, old96[0], old96[0] / ref[0]), flush=True)
    names = ["total"] + vecs
    keep = [j for j in range(len(names)) if j == 0 or ref[j] > 0]   # every contribution
    wn = max(abs(res["new"][j] / ref[j] - 1) for j in keep)
    jo = max(keep, key=lambda j: abs(math.log(res["old"][j] / ref[j])))
    print("PNQS %s: total at 1 kHz: quasi-static (%d x .noise) %.5e, .pnoise %.5e (%+.1e), .pnoise PAC_NOISE_STALE_STATE=1 %.5e (x%.2f); "
          "%d contributions: .pnoise worst %.1e, stale-state worst x%.1f (%s: %.4e vs %.4e)"
          % (name, N, ref[0], res["new"][0], res["new"][0] / ref[0] - 1, res["old"][0], res["old"][0] / ref[0], len(keep), wn,
             res["old"][jo] / ref[jo], names[jo], res["old"][jo], ref[jo]), flush=True)


if __name__ == "__main__":
    ng, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for n in ("diode", "bjt"):
        lti(ng, out, n)
    for n in ("diode", "bjt"):
        qs(ng, out, n)
