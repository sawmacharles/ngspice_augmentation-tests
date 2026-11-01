#!/usr/bin/env python3
"""50_xspice_pac_through.py -- handbook example 50 (card #560): td `.pac` / `.pstb` / `.pxf` THROUGH stateful analog XSPICE models (card #557, patch 0115).

    python3 50_xspice_pac_through.py NGSPICE OUTDIR

  1. R3 voltage-mode buck whose type-III compensator is one `s_xfer` (50_xspice_r3_sxfer.cir): `.pac` control-to-output, 24 points 4..50 kHz, pacpoints=2000,
     against a frequency-response measurement (FRA) of the real nonlinear circuit -- the FRA values are CITED (docs/pss-pac.md card #557: one 0.01 V sine in
     series at Vctrl, .tran reltol 1e-5 maxstep 1 ns, single-bin DFT over the last 2 ms, amplitude 0.0025 repeats it to 0.0004 dB / 0.004 deg).
     Also the default recording grid, and PAC_HID_RED=1 (the hidden states do not act on the nodes: red proof).                           (R3XS lines)
  2. a linear s_xfer on a periodic orbit: td `.pac` against plain `.ac` (the exact backward-Euler discrete transfer of the model is the
     reference of tests/test_xspice.py; here the first-order convergence to `.ac` and its Richardson extrapolation), at pacpoints 50000 / 100000 (td vs `.ac` is first order in h)                               (XSLTI lines)
  3. refusals: a `slew` model that slews on the orbit, a `hyst` that flips, `.pnoise` through s_xfer, `.pxf adjoint=1` (falls back to the per-source
     engine, bit-identical to `.pac`).                                                                                                   (XSREF lines)
"""
import cmath
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FRA = {4000.0: (-6.1833, 12.682), 20000.0: (-11.7160, -95.560), 50000.0: (-25.2293, -131.765)}   # cited: dB, deg of vout / Vctrl
WR = ".control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata %s %s\nquit\n.endc\n.end\n"


def run(ng, out, name, deck, env=None):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "deck.cir"), "w").write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.pop("PAC_HID_RED", None)
    e.update(env or {})
    t0 = time.time()
    p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", "deck.cir"],
                       cwd=d, env=e, capture_output=True, text=True, timeout=1800)
    o = (p.stdout + p.stderr).replace("\r", "\n")
    open(os.path.join(d, "stdout.txt"), "w").write(o)
    return d, o, time.time() - t0


def rdc(path):
    rows = {}
    for ln in open(path).read().splitlines()[1:]:
        t = [float(x) for x in ln.split()]
        if len(t) >= 3:
            rows[t[0]] = [complex(t[k], t[k + 1]) for k in range(1, len(t) - 1, 2)]
    return rows


def db_deg(h):
    return 20 * math.log10(abs(h)), math.degrees(cmath.phase(h))


def main(ng, out):
    base = open(os.path.join(HERE, "50_xspice_r3_sxfer.cir")).read()
    deck = base.rstrip() + "\n.options reltol=1e-4\n"
    card = '.pac 500e3 1e-3 vout 200 10 150 1e-6 lin 24 4000 50000 method="shoot" solver="td"'
    for tag, extra, env in (("pp2000", " pacpoints=2000", None), ("default grid", "", None), ("PAC_HID_RED=1", " pacpoints=2000", {"PAC_HID_RED": "1"})):
        d, o, wall = run(ng, out, "r3xs_" + tag.replace(" ", "_").replace("=", ""), deck + card + extra + "\n" + WR % ("out.txt", "v(vout)"), env)
        rows = rdc(os.path.join(d, "out.txt"))
        nh = re.search(r"(\d+) hidden states", o)
        parts = []
        for f, (db, deg) in FRA.items():
            x, y = db_deg(rows[f][0])
            parts.append("%g kHz %+.4f dB / %+.3f deg" % (f / 1e3, x - db, (y - deg + 180) % 360 - 180))
        print("R3XS %s: %s hidden states; difference to the cited FRA: %s; wall %.1f s" % (tag, nh.group(1) if nh else "no", "; ".join(parts), wall), flush=True)
    # ---- linear s_xfer on a periodic orbit
    src = ("* stateful xspice s_xfer\nVin in 0 SIN(0.3 0.2 1k) AC 1\nRdum in dum 1meg\nCdum dum 0 1p\nA1 in o1 m\n"
           ".model m s_xfer(num_coeff=[1] den_coeff=[2.5e-12 2.2e-6 1] int_ic=[0 0])\nRl o1 0 1k\n.options reltol=1e-6 abstol=1e-15 vntol=1e-9\n")
    d, o, wall = run(ng, out, "lti_ac", src + ".ac dec 5 1k 100k\n" + WR % ("out.txt", "o1"))
    ac = rdc(os.path.join(d, "out.txt"))
    res = {}
    for pp in (50000, 100000):
        d, o, wall = run(ng, out, "lti_pp%d" % pp, src + '.pac 1k 0 o1 2000 1 8 1e-6 dec 5 1k 100k method="shoot" solver="td" pacpoints=%d\n' % pp + WR % ("out.txt", "v(o1)"))
        res[pp] = rdc(os.path.join(d, "out.txt"))

    for pp, rows in res.items():
        w_ac = max(abs(rows[f][0] - ac[f][0]) / abs(ac[f][0]) for f in rows if f in ac)
        print("XSLTI s_xfer td .pac at pacpoints=%d: worst relative difference to plain .ac %.2e over %d points" % (pp, w_ac, len(rows)), flush=True)
    rich = max(abs(2 * res[100000][f][0] - res[50000][f][0] - ac[f][0]) / abs(ac[f][0]) for f in res[50000] if f in ac)
    print("XSLTI Richardson 2 td(pacpoints 100000) - td(pacpoints 50000) vs .ac: %.2e" % rich, flush=True)
    # ---- refusals
    opts = 'method="shoot" solver="td" pacpoints=100000'
    mod = ".model m slew(rise_slope=3000 fall_slope=3000)"
    d, o, w = run(ng, out, "ref_slew", "* slew slewing\nVin in 0 SIN(0 0.5 1k) AC 1\nRdum in dum 1meg\nCdum dum 0 1p\nA1 in o1 m\n" + mod + "\nRl o1 0 1k\n.options reltol=1e-6\n"
                  ".pac 1k 0 o1 2000 1 8 1e-6 lin 3 1k 50k " + opts + "\n.control\nrun\nquit\n.endc\n.end\n")
    l = [x.strip() for x in o.splitlines() if "REFUSED" in x]
    print("XSREF slew model that slews on the orbit: " + (l[0][:240] if l else "NOT REFUSED"), flush=True)
    mod = ".model m hyst(in_low=0 in_high=1 hyst=0.1 out_lower_limit=0 out_upper_limit=1 input_domain=0.05 fraction=FALSE)"
    d, o, w = run(ng, out, "ref_hyst", "* hyst flipping\nVin in 0 SIN(0.5 1.5 1k) AC 1\nRdum in dum 1meg\nCdum dum 0 1p\nA1 in o1 m\n" + mod + "\nRl o1 0 1k\n.options reltol=1e-6\n"
                  ".pac 1k 0 o1 2000 1 8 1e-6 lin 3 1k 50k " + opts + "\n.control\nrun\nquit\n.endc\n.end\n")
    l = [x.strip() for x in o.splitlines() if "REFUSED" in x]
    print("XSREF hyst model that flips on the orbit: " + (l[0][:240] if l else "NOT REFUSED"), flush=True)
    d, o, w = run(ng, out, "ref_pnoise", src + ".pnoise 1k 0 o1 2000 1 8 1e-6 o1 Vin lin 3 1k 50k " + opts + "\n.control\nrun\nquit\n.endc\n.end\n")
    l = [x.strip() for x in o.splitlines() if "REFUSED" in x]
    print("XSREF .pnoise through s_xfer: " + (l[0][:240] if l else "NOT REFUSED"), flush=True)
    outs = {}
    for tag, tail, vec in (("a1", '.pxf 1k 0 o1 2000 1 8 1e-6 o1 lin 3 1k 50k sources="Vin" adjoint=1 ' + opts, "pxf_Vin"), ("a0", '.pxf 1k 0 o1 2000 1 8 1e-6 o1 lin 3 1k 50k sources="Vin" adjoint=0 ' + opts, "pxf_Vin"),
                           ("pac", '.pac 1k 0 o1 2000 1 8 1e-6 lin 3 1k 50k ' + opts, "v(o1)")):
        d, o, w = run(ng, out, "pxf_" + tag, src + tail + "\n" + WR % ("out.txt", vec))
        outs[tag] = (rdc(os.path.join(d, "out.txt")), o)
    fb = "adjoint engine REFUSED" in outs["a1"][1] and "does not transpose its hidden states" in outs["a1"][1]
    w1 = max(abs(outs["a1"][0][f][0] - outs["pac"][0][f][0]) / abs(outs["pac"][0][f][0]) for f in outs["pac"][0])
    print("XSREF .pxf adjoint=1 through s_xfer: falls back to the per-source engine (message: %s), differs from .pac by %.1e" % (fb, w1), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
