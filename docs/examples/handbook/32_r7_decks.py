#!/usr/bin/env python3
"""32_r7_decks.py -- writes the handbook's transistor-level (BSIM4) buck decks (card #543, section 6b).

    python3 32_r7_decks.py OUTDIR

Every deck is refs/r7_bsim_buck.cir (card #539) with its own `.control` block cut off and ONLY these edits, nothing else:
  * the `../stb/lgprobe.lib` include path made absolute (the decks are written outside refs/);
  * `.options reltol=1e-4` (see the page: `.pss` on BSIM4 decks needs reltol <= 3e-4);
  * the analysis card (`.pss` / `.pstb` / `.pac` / `.pxf`) and a `.control` that writes the result with `wrdata`;
  * for `.pstb`, `.pac` and `.pxf`: the Tian `lgprobe` subcircuit replaced by a plain 0 V source `Vprobe fbtop vout DC 0`
    (a `.pstb` probe is a plain V source; leaving `lgprobe` in a deck used for `.pac` adds a second, unwanted AC injection);
  * for `.pac` / `.pxf`: a 0 V source `Vctrl` spliced into the comparator's opout input (control-to-output) and a current
    source `Iinj` at vout (output impedance); `Vin` is the line input. `.pac` takes the one that carries `AC 1`.

Decks written (names are what run_all.sh runs):
  r7_pss.cir          .pss method=shoot, reltol 1e-4, Points 400           -> orbit.txt (t vout il)
  r7_pss_reltol3e-4.cir  same at reltol 3e-4                                  (the documented limit)
  r7_pss_reltol1e-3.cir  same at the ngspice DEFAULT reltol 1e-3, at most 12 Newton iterations (refs/README.md measured a stall here
                      on s88; on s96 it converges -- the page reports what the run shows)
  r7_pstb.cir         .pstb lin 400 500 200000 (500 Hz grid: every FRA point is on it)  -> pstb.txt (f t)
  r7_pstb_reltol1e-3.cir  the same .pstb at the ngspice DEFAULT reltol 1e-3                -> pstb_1e-3.txt
  r7_pac.cir          .pac  ctrl-to-output, same grid                           -> pac.txt
  r7_pxf_fwd.cir / r7_pxf_adj.cir   .pxf sources=all, adjoint=0 / adjoint=1, same grid -> pxf_fwd.txt / pxf_adj.txt
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
REF = os.path.join(ROOT, "refs", "r7_bsim_buck.cir")
LIB = os.path.join(ROOT, "stb", "lgprobe.lib")

PSS = "500e3 300e-6 vout 400 10 {iters} 1e-6"
GRID = "lin 400 500 200000"
CTL = "set filetype=ascii\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"


def body():
    s = open(REF).read()
    b = s[:s.index(".control")]
    assert "../stb/lgprobe.lib" in b and "XP fbtop vout lgprobe" in b
    return b.replace("../stb/lgprobe.lib", LIB)


def probe_body(b):
    b = b.replace("XP fbtop vout lgprobe", "Vprobe fbtop vout DC 0")
    return b.replace(".include " + LIB + "\n", "")


def splice(b, ctrl=False, line=False, zout=False):
    if ctrl:
        assert "Mc1 cd1 opout ctail" in b
        b = b.replace("Mc1 cd1 opout ctail", "Vctrl ctrlnode opout DC 0%s\nMc1 cd1 ctrlnode ctail" % (" AC 1" if ctrl == "ac" else ""))
    if line:
        b = b.replace("Vin vin 0 DC {vin_v}", "Vin vin 0 DC {vin_v} AC 1")
    if zout:
        b = b.replace("Rload vout 0 {rload}", "Rload vout 0 {rload}\nIinj vout 0 DC 0%s" % (" AC 1" if zout == "ac" else ""))
    return b


def deck(b, opts, card, wr):
    return b + ".options %s\n%s\n.control\n%s run\n%s\n.endc\n.end\n" % (opts, card, CTL, wr)


def main(out):
    os.makedirs(out, exist_ok=True)
    b = body()
    pb = probe_body(b)
    w = lambda name, txt: open(os.path.join(out, name), "w").write(txt)
    pss = lambda it: ".pss " + PSS.format(iters=it) + " uic method=shoot"
    wr_orbit = "wrdata orbit.txt v(vout) i(L1)"
    w("r7_pss.cir", deck(b, "reltol=1e-4", pss(150), wr_orbit))
    w("r7_pss_reltol3e-4.cir", deck(b, "reltol=3e-4", pss(150), "wrdata orbit_3e-4.txt v(vout) i(L1)"))
    w("r7_pss_reltol1e-3.cir", deck(b, "reltol=1e-3", pss(12), "wrdata orbit_1e-3.txt v(vout) i(L1)"))
    w("r7_pstb.cir", deck(pb, "reltol=1e-4", '.pstb ' + PSS.format(iters=150) + ' ' + GRID + ' probe="Vprobe" method="shoot"',
                          "wrdata pstb.txt t"))
    w("r7_pstb_reltol1e-3.cir", deck(pb, "reltol=1e-3", '.pstb ' + PSS.format(iters=150) + ' ' + GRID + ' probe="Vprobe" method="shoot"',
                                     "wrdata pstb_1e-3.txt t"))
    w("r7_pac.cir", deck(splice(pb, ctrl="ac"), "reltol=1e-4", '.pac ' + PSS.format(iters=150) + ' ' + GRID + ' method="shoot" solver="td"',
                         "wrdata pac.txt v(vout)"))
    pxb = splice(pb, ctrl="dc", zout="dc")
    for tag, adj in (("fwd", 0), ("adj", 1)):
        w("r7_pxf_%s.cir" % tag, deck(pxb, "reltol=1e-4",
                                      '.pxf ' + PSS.format(iters=150) + ' vout ' + GRID + ' method="shoot" solver="td" sources="all" adjoint=%d' % adj,
                                      "wrdata pxf_%s.txt pxf_vin pxf_vctrl pxf_iinj" % tag))


def tran_trap(out):
    """r7_tran_old.cir: the ramp line as it was before card #541 (tr + pw + tf = period + 1 ns); r7_tran_fixed.cir: the shipped
    line (tf shortened by 1 ns). Both: .options reltol=1e-6 abstol=1e-9, .tran 1n 30u uic."""
    os.makedirs(out, exist_ok=True)
    b = body()
    new = "Vramp ramp 0 PULSE({vrlo} {vrhi} 0 {0.97*tclk} {0.03*tclk-1n} 1n {tclk})"
    old = "Vramp ramp 0 PULSE({vrlo} {vrhi} 0 {0.97*tclk} {0.03*tclk} 1n {tclk})"
    assert new in b
    for tag, line in (("old", old), ("fixed", new)):
        t = b.replace(new, line)
        open(os.path.join(out, "r7_tran_%s.cir" % tag), "w").write(
            t + ".options reltol=1e-6 abstol=1e-9\n.control\ntran 1n 30u 0 uic\nmeas tran vo avg v(vout) from=20u to=30u\n.endc\n.end\n")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--tran-trap":
        tran_trap(sys.argv[2])
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else ".")
