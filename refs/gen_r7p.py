#!/usr/bin/env python3
"""gen_r7p.py -- build refs/r7p_bsim_buck_pex.cir (card #539) from refs/r7_bsim_buck.cir.

R7p = R7 plus three synthetic extracted-parasitic RC trees (bench/gen_pex_buck.py's tree_lines(), the same generator
the PEX-sized bench decks use; deterministic, no RNG): one hung on `vout` (output-plane parasitics), one on the switch
node `sw`, one on the high-side gate net `ghs` (distributed gate-routing R/C). branch=3, depth=4: 121 nodes, 120 trace
resistors, 121 capacitors, ~40 sibling/coupling elements each: ~375 extra explicit capacitors (states for the shooter)
and ~900 extra elements. Not derived from any real layout. Run: python3 refs/gen_r7p.py (rewrites the deck).
"""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench"))
from gen_pex_buck import tree_lines

TREES = [  # tag, attach net, series R root coupling, rtrace0, rgrow, cplane0, cshrink, coupling_c
    ("tv", "vout", "20m", 0.05, 1.5, 200e-12, 0.6, 20e-12),
    ("ts", "sw", "5m", 0.02, 1.5, 20e-12, 0.6, 2e-12),
    ("tg", "ghs", "0.5", 0.5, 1.5, 2e-12, 0.6, 0.2e-12),
]

def build():
    base = open(os.path.join(HERE, "r7_bsim_buck.cir")).read()
    head, ctl = base.split(".control")[0], ".control" + base.split(".control", 1)[1]
    out = ["\n* ---- R7p: synthetic extracted parasitic RC trees (refs/gen_r7p.py, bench/gen_pex_buck.py tree_lines) ----\n"]
    ncap = nres = 0
    for tag, net, rc, r0, rg, c0, cs, cc in TREES:
        lines, root, far, n = tree_lines(3, 4, r0, rg, c0, cs, 3, 12, cc)
        out.append(f"Rcpl_{tag} {net} {tag}_0 {rc}\n")
        for ln in lines:
            if ln.startswith("*"):
                out.append(ln); continue
            toks = ln.split()
            toks[0] += "_" + tag
            toks[1] = toks[1].replace("tree_", tag + "_"); toks[2] = toks[2].replace("tree_", tag + "_")
            if toks[0].startswith("C"): ncap += 1
            else: nres += 1
            out.append(" ".join(toks) + "\n")
    out.append(f"* {ncap} parasitic capacitors, {nres} parasitic resistors\n\n")
    txt = head + "".join(out) + ctl
    txt = txt.replace("R7 transistor-level synchronous buck", "R7p (R7 + parasitic RC trees) transistor-level synchronous buck", 1)
    open(os.path.join(HERE, "r7p_bsim_buck_pex.cir"), "w").write(txt)
    return ncap, nres

if __name__ == "__main__":
    print("caps %d res %d" % build())
