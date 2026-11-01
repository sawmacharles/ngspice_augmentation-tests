#!/usr/bin/env python3
"""Bug 3 (card #539): `.pxf ... solver=td adjoint=1` disagrees with adjoint=0 (and with .pac) for a source that sets a MOS
threshold-crossing time, even with level-1 MOSFETs without any capacitance. Needs helpers.py (copy of bug 1's run_repro.py)."""
import os
import helpers as R
def pxf(st, models, adj, srcs):
    s = " ".join('sources="%s"' % x for x in srcs)
    an = ('.pxf 500e3 20e-6 out 400 10 50 1e-7 out lin 1 1000 1000 method="shoot" solver="td" %s adjoint=%d\n.options reltol=1e-4\n'
          '.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata pxf.txt %s\n.endc\n.end\n') % (s, adj, " ".join("pxf_" + x.lower() for x in srcs))
    R.run("pxf", R.deck(st, models, vc="2.5", analysis=an))
    r = open(os.path.join(R.W, "pxf.txt")).read().split("\n")[1].split()
    return [abs(complex(float(r[1 + 2 * k]), float(r[2 + 2 * k]))) for k in range(len(srcs))]
for m in ("level1_cg0", "level1_cg1e-8", "bsim4"):
    for st in (1, 2):
        print("%-14s N=%d  |pxf_vc| forward %.4f  adjoint %.4f   (static d<v>/dVc: ~0.5)" % (m, st, pxf(st, R.MODELS[m], 0, ("Vc",))[0], pxf(st, R.MODELS[m], 1, ("Vc",))[0]), flush=True)
