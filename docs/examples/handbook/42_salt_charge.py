#!/usr/bin/env python3
"""42_salt_charge.py -- handbook example 42 (card #552): ideal-switch (saltation) events on an orbit that also carries device-model charges
(card #544, patch 0098) and the hidden-state continuation that makes it exact (card #546, patch 0100).

    python3 42_salt_charge.py NGSPICE OUTDIR

42_r6o_ctrl.cir = R6o, the open-loop ideal-`sw` buck of example 10 (12 V -> about 3.2 V, 500 kHz, `vctrl_dc` the control DC).
42_r6ob_ctrl_follower.cir = the same buck with a BSIM4 (capmod 2, W 100 u) source follower and a 2 pF gate capacitor in its control path:
Vctrl -> Rgs 1 Meg -> gate (2 pF) -> follower -> `ctrlnode`. The follower is linear at its constant bias and decoupled from the converter, so
the exact result is H(R6ob) = F(f) * H(R6o at the same control DC) with F(f) = v(ctrlnode)/Vctrl from a plain `.ac`. Prints, at 20 and 100 kHz:
  SALTID f: |H_b - F H_o| / |F H_o|                  (default; PSS_NO_HIDDEN=1 = the s99 restart; PAC_NATIVE_HISTORY=1 = the pre-#542 history)
  SALTREF f: .pac of R6ob against F * the independent pac_perturbed value of R6o (tests/test_pac_charge.py PP_R6O, cited, not rerun)
and the two refusal/no-regression checks of an ideal switch beside a charged / a chargeless diode (SALTDIODE lines).
"""
import cmath
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# tools/pac_perturbed.py (.pss backend, amplitude halving and step refinement <= 2e-4 dB) on 42_r6o_ctrl.cir with vctrl_dc = 0.5285212: dB, deg of vout
PP_R6O = {20000: (-7.8049, -169.9099), 100000: (-34.9237, -147.3063)}
WR = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata %s %s\nquit\n.endc\n.end\n"


def run(ngspice, out, name, deck, env=None):
    p = os.path.join(out, name + ".cir")
    open(p, "w").write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ngspice, "-b", p],
                       cwd=out, env=e, capture_output=True, text=True, timeout=1800)
    return (r.stdout + r.stderr).replace("\r", "\n")


def rd1(out, fn):
    ln = [l for l in open(os.path.join(out, fn)).read().splitlines()[1:] if l.strip()][0].split()
    return complex(float(ln[1]), float(ln[2]))


def body(name):
    return open(os.path.join(HERE, name)).read().replace("\n.end\n", "\n")


def dbdeg(z):
    return 20 * math.log10(abs(z)), math.degrees(cmath.phase(z))


def pac(f):
    return '.pac 500e3 100e-6 vout 200 10 150 1e-6 lin 1 %d %d method="shoot" solver="td"\n' % (f, f)


def main(ng, out):
    os.makedirs(out, exist_ok=True)
    b4, o6 = body("42_r6ob_ctrl_follower.cir"), body("42_r6o_ctrl.cir")
    txt = run(ng, out, "op", b4 + ".control\nop\nprint v(ctrlnode)\n.endc\n.end\n")
    vc = re.search(r"v\(ctrlnode\) = ([-0-9.e+]+)", txt).group(1)
    o6 = o6.replace("vctrl_dc=0.55", "vctrl_dc=%s" % vc)
    print("SALT control DC of the follower (v(ctrlnode) at the bias): %s V -> R6o is run at vctrl_dc=%s" % (vc, vc), flush=True)
    cases = (("default", None), ("PSS_NO_HIDDEN=1", {"PSS_NO_HIDDEN": "1"}), ("PAC_NATIVE_HISTORY=1", {"PAC_NATIVE_HISTORY": "1"}))
    for f in (20000, 100000):
        run(ng, out, "ac_%d" % f, b4 + ".ac lin 1 %d %d\n" % (f, f) + WR % ("ac_%d.txt" % f, "v(ctrlnode)"))
        F = rd1(out, "ac_%d.txt" % f)
        run(ng, out, "o_%d" % f, o6 + pac(f) + WR % ("o_%d.txt" % f, "v(vout)"))
        Ho = rd1(out, "o_%d.txt" % f)
        h = 10 ** (PP_R6O[f][0] / 20) * cmath.exp(1j * math.radians(PP_R6O[f][1]))
        rdb, rdeg = dbdeg(F * h)
        for tag, env in cases:
            t = run(ng, out, "b_%d" % f, b4 + pac(f) + WR % ("b_%d.txt" % f, "v(vout)"), env)
            gl = "generic charge-based history" in t
            Hb = rd1(out, "b_%d.txt" % f)
            e = abs(Hb - F * Ho) / abs(F * Ho)
            d, p = dbdeg(Hb)
            print("SALTID %d Hz %s: |H(R6ob) - F H(R6o)| / |F H(R6o)| = %.2e; R6ob %.3f dB / %.2f deg, F x pac_perturbed(R6o) %.3f / %.2f: %+.3f dB / %+.2f deg%s"
                  % (f, tag, e, d, p, rdb, rdeg, d - rdb, (p - rdeg + 180) % 360 - 180, "" if gl else "  [no charge-history line: native path]"), flush=True)
    # an ideal switch beside a diode, no capacitor at the switch node (s118, card #575): a charged diode runs on the fast-node path (it was
    # refused up to s117), the same with an RC snubber at the node runs since s144 (card #713; refused s118-s143), a chargeless diode keeps
    # the native path
    base = ("* sw + diode\nVin vin 0 DC 12\nVramp vramp 0 PULSE(0 1 0 2u 1n 1n 2u)\nVctrl vctrl 0 DC 0.5 AC 1\n"
            ".model SWMOD SW(vt=0 vh=0 ron=0.05 roff=1e6)\nS1 vin sw vctrl vramp SWMOD\nD1 0 sw DMOD\n.model DMOD D(IS=1e-14 RS=0.02 N=1.05 %s)\n"
            "L1 sw vout 10u\nResr vout voutc 0.01\nCout voutc 0 100u IC=6\nRload vout 0 5\n%s"
            '.pac 500e3 100e-6 vout 200 10 150 1e-6 lin 1 20000 20000 method="shoot" solver="td"\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata sw.txt v(vout)\n.endc\n.end\n')
    # tools/pac_perturbed.py on the charged deck (cited: tests/test_pac_charge.py D575 swc_50p_1n, 2026-10-06, s117 binary): -1.1907 dB / -168.9907 deg at 20 kHz
    # the snubber deck (10 ohm / 1 nF): tools/pac_perturbed.py, s144 binary, 2026-10-19 (halving 0.0006 dB, step residual 0.0001 dB / 0.002 deg):
    # -1.1902 dB / -168.9719 deg at 20 kHz
    refs = {"": (-1.1907, -168.9907), "Rsn sw ns 10\nCsn ns 0 1n\n": (-1.1902, -168.9719)}
    for tag, mod, extra in (("charged diode (CJO=50p TT=1n)", "CJO=50p TT=1n", ""), ("charged diode with an RC snubber", "CJO=50p TT=1n", "Rsn sw ns 10\nCsn ns 0 1n\n"),
                            ("chargeless diode", "", "")):
        rdb, rdeg = refs[extra]
        try: os.remove(os.path.join(out, "sw.txt"))
        except OSError: pass
        t = run(ng, out, "swd", base % (mod, extra))
        ref = [l for l in t.splitlines() if "REFUSED" in l]
        if ref:
            res = ref[0].strip()[:230]
        elif "PAC-TD: done" not in t:
            res = "UNEXPECTED: " + t[-200:]
        elif not mod:
            res = "no refusal; " + ("generic recursion (any switch node with a diode, s123)" if "generic charge-based" in t else "native path")
        elif "generic charge-based" in t:
            d, p = dbdeg(rd1(out, "sw.txt"))
            res = "runs, fast-node path (inductor-routed jump: %s): %.4f dB / %.3f deg at 20 kHz, pac_perturbed %.4f / %.3f: %+.3f dB / %+.2f deg" % (
                "yes" if "saltation via inductor-routed jump" in t else "NO", d, p, rdb, rdeg, d - rdb, (p - rdeg + 180) % 360 - 180)
        else:
            res = "no refusal; native path"
        print("SALTDIODE %s: %s" % (tag, res), flush=True)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
