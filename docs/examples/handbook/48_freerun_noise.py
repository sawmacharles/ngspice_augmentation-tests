#!/usr/bin/env python3
"""48_freerun_noise.py -- handbook example 48 (card #560): `.pxf adjoint=1` and `.pnoise` on a FREE-RUNNING (trigger=) orbit (card #557, patches 0114, 0116).

    python3 48_freerun_noise.py NGSPICE OUTDIR

48_freerun_noise.cir is the R5 hysteretic buck (example 30) without analysis cards.  Everything is appended here, `trigger="q" trigval=0.5 trigdir="rise"`:
  1. `.pxf 427350 450e-6 vout 200 10 45 1e-4 vout lin 3 20000 90000 sources="Vin" sources="Vrefsrc" adjoint=1` against adjoint=0 (the per-source forward
     engine; with PAC_EXACT_JACOBIAN=1 the fair comparison of card #504, and with its default): worst relative difference.         (ADJ lines)
     Red proof PAC_FREERUN_PROJ=1: the forward operator with the period-variation projection is a different operator; the adjoint refuses.
  2. `.pnoise ... vout Vin lin N f1 f2` at 10 Hz .. 1.8 kHz (the phase mode does not diverge in band), and in three 9-bin bands (610.35 Hz bins) at
     15.3 / 40.3 / 84.8 kHz = 0.036 / 0.095 / 0.20 fsw, in V/rtHz (the values of the transient-noise ensemble of tools/pnoise_mc_r5.py are CITED).   (PN lines)
  3. the band guard: a sweep whose top frequency is above 0.25 fsw is refused.                                                              (GUARD line)
"""
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TRIG = 'method="shoot" solver="td" trigger="q" trigval=0.5 trigdir="rise"'
HEAD = "427350 450e-6 vout 200 10 45 1e-4"
FSW = 422681.0   # the orbit frequency of example 30 (R5), Hz
# CITED: transient-noise ensemble of tools/pnoise_mc_r5.py (docs/pss-pac.md card #557): K = 50, TRNOISE step ts, V/rtHz per band (16 seeds x 20 ms at 5 ns, 12 x 16 ms at 2.5 ns)
ENS = {5.0: {15258.8: 4.92996e-09, 40283.2: 5.06440e-09, 84838.9: 5.37193e-09}, 2.5: {15258.8: 5.13275e-09, 40283.2: 5.38977e-09, 84838.9: 5.72514e-09}}
BANDS = {15258.8: (12817.3828125, 17700.1953125), 40283.2: (37841.796875, 42724.609375), 84838.9: (82397.4609375, 87280.2734375)}


def run(ng, out, name, card, vecs, env=None, timeout=900):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    body = open(os.path.join(HERE, "48_freerun_noise.cir")).read()
    deck = body.rstrip() + "\n" + card + "\n.control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\n"
    deck += ("wrdata out.txt " + " ".join(vecs) + "\n" if vecs else "") + "quit\n.endc\n.end\n"
    open(os.path.join(d, "deck.cir"), "w").write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.pop("PAC_EXACT_JACOBIAN", None)
    e.pop("PAC_FREERUN_PROJ", None)
    e.update(env or {})
    p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", "deck.cir"],
                       cwd=d, env=e, capture_output=True, text=True, timeout=timeout)
    o = (p.stdout + p.stderr).replace("\r", "\n")
    open(os.path.join(d, "stdout.txt"), "w").write(o)
    return d, o


def rows(d, ncomplex=0, nreal=0):
    out = []
    for ln in open(os.path.join(d, "out.txt")).read().splitlines()[1:]:
        t = [float(x) for x in ln.split()]
        if len(t) > 1:
            out.append(t)
    return out


def worst_rel(a, b):
    w = 0.0
    for ra, rb in zip(a, b):
        for k in range(1, len(ra) - 1, 2):
            x, y = complex(ra[k], ra[k + 1]), complex(rb[k], rb[k + 1])
            w = max(w, abs(x - y) / abs(y))
    return w


def main(ng, out):
    # ---- 1. adjoint = forward
    vecs = ["pxf_Vin", "pxf_Vrefsrc"]
    sw = 'vout lin 3 20000 90000 sources="Vin" sources="Vrefsrc"'
    res = {}
    for tag, env in (("exact", {"PAC_EXACT_JACOBIAN": "1"}), ("default", {})):
        for adj in (1, 0):
            d, o = run(ng, out, "pxf_%s_a%d" % (tag, adj), ".pxf %s %s adjoint=%d %s" % (HEAD, sw, adj, TRIG), vecs, env)
            res[(tag, adj)] = (rows(d), o)
    for tag in ("exact", "default"):
        a, f = res[(tag, 1)][0], res[(tag, 0)][0]
        did = "adjoint done" in res[(tag, 1)][1]
        print("ADJ %s forward: .pxf adjoint=1 %s the adjoint engine; worst relative difference to adjoint=0 over Vin and Vrefsrc at 20/55/90 kHz: %.2e"
              % ("(PAC_EXACT_JACOBIAN=1)" if tag == "exact" else "(default forward)", "ran" if did else "DID NOT RUN", worst_rel(a, f)), flush=True)
    h = res[("exact", 0)][0]
    d, o = run(ng, out, "pxf_proj_a0", ".pxf %s %s adjoint=0 %s" % (HEAD, sw, TRIG), vecs, {"PAC_EXACT_JACOBIAN": "1", "PAC_FREERUN_PROJ": "1"})
    fp = rows(d)
    d, o2 = run(ng, out, "pxf_proj_a1", ".pxf %s %s adjoint=1 %s" % (HEAD, sw, TRIG), vecs, {"PAC_EXACT_JACOBIAN": "1", "PAC_FREERUN_PROJ": "1"})
    print("ADJ red proof PAC_FREERUN_PROJ=1: the forward answer moves by a relative %.2f; adjoint=1 with it %s" %
          (worst_rel(fp, h), "REFUSES and falls back to the per-source engine" if ("adjoint engine REFUSED" in o2 and "adjoint done" not in o2) else "DID NOT REFUSE"), flush=True)
    # ---- 2. pnoise
    d, o = run(ng, out, "pn_low", ".pnoise %s vout Vin lin 5 10 1800 %s" % (HEAD, TRIG), ["onoise_spectrum"])
    low = rows(d)
    print("PN low frequency, 10 Hz .. 1.8 kHz: " + ", ".join("%.0f Hz %.4e" % (r[0], r[1]) for r in low) + " V/rtHz; spread %.1f %%" % (100 * (max(r[1] for r in low) / min(r[1] for r in low) - 1)), flush=True)
    pw = {}
    for c, (lo, hi) in BANDS.items():
        d, o = run(ng, out, "pn_band%d" % int(c), ".pnoise %s vout Vin lin 9 %.10g %.10g %s" % (HEAD, lo, hi, TRIG), ["onoise_spectrum"])
        v = [r[1] for r in rows(d)]
        pw[c] = math.sqrt(sum(x * x for x in v) / len(v))
        e5, e25 = ENS[5.0][c], ENS[2.5][c]
        print("PN band %.1f kHz (%.3f fsw): .pnoise %.4e V/rtHz (9-bin rms) x K=50 = %.4e; cited ensemble %.4e (ts 5 ns, %+.2f dB) and %.4e (ts 2.5 ns, %+.2f dB)"
              % (c / 1e3, c / FSW, pw[c], pw[c] * 50, e5, 20 * math.log10(e5 / (pw[c] * 50)), e25, 20 * math.log10(e25 / (pw[c] * 50))), flush=True)
    # ---- 3. guard
    d, o = run(ng, out, "pn_guard", ".pnoise %s vout Vin lin 3 50000 200000 %s" % (HEAD, TRIG), ["onoise_spectrum"])
    l = [x.strip() for x in o.splitlines() if "REFUSED" in x]
    print("GUARD .pnoise up to 200 kHz (0.47 fsw): " + (l[0][:230] if l else "NOT REFUSED"), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
