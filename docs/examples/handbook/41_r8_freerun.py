#!/usr/bin/env python3
"""41_r8_freerun.py -- handbook example 41 (card #552): R8, the transistor-level (BSIM4) FREE-RUNNING hysteretic buck (cards #544, #549).

    python3 41_r8_freerun.py NGSPICE OUTDIR

Deck: refs/r8_bsim_hysteretic_buck.cir (R7's BSIM4 power stage, drivers, dead-time generator and comparator; hysteresis by resistive positive
feedback; 5 V -> 1.7055 V, fsw about 395 kHz) with its own `.control` cut off and only the analysis cards / `.options reltol=1e-6` added.
Runs, one ngspice at a time, OMP_NUM_THREADS=2:
  T  settled `.tran` (reltol 1e-6, 2 ns grid, 300 us; the last 20 periods are used): period = mean spacing of the rising crossings of cx
     through 1.7055 V (the trigger of the `.pss`), orbit = the last period between two crossings
  A  `.pss ... method=shoot trigger=cx trigval=1.7055 trigdir=rise` (default: the hidden-state continuation of the trigger runs, card #549)
  B  the same with PSS_NO_HIDDEN=1 (every trigger run restarts from the native C/L states only: s101 behaviour, "before")
  C  `.pac solver=td` control-to-output at 20/50/80 kHz on the free-running orbit (generic model-charge path, card #544): printed against the
     reference of tests/test_pac_charge.py (tools/pac_perturbed.py --freerun, two long nonlinear transients per point, 8-11 min a point: cited, NOT rerun)
Prints lines R8TRAN / R8PSS / R8ORB / R8PAC.
"""
import cmath
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
REF = os.path.join(ROOT, "refs", "r8_bsim_hysteretic_buck.cir")
CTL = "set filetype=ascii\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
# tools/pac_perturbed.py --freerun on the same circuit (tests/test_pac_charge.py R8_REF): dB, deg, amplitude-halving spread dB
PP = {20000: (7.8977, -0.2728, 0.009), 50000: (7.9969, -0.8614, 0.009), 80000: (8.1740, -1.4905, 0.008)}
LEVEL = 1.7055


def sh(ngspice, deck, cwd, env=None, log=None):
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    t0 = time.time()
    p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ngspice, "-b", deck],
                       cwd=cwd, env=e, capture_output=True, text=True, timeout=3600)
    out = (p.stdout + p.stderr).replace("\r", "\n")
    if log:
        open(log, "w").write(out)
    return out, time.time() - t0


def rd(path):
    rows = []
    for ln in open(path).read().splitlines()[1:]:
        v = ln.split()
        if len(v) >= 3:
            rows.append([float(x) for x in v])
    return rows


def cross(ts, vs, lev):
    out = []
    for i in range(1, len(ts)):
        if vs[i - 1] < lev <= vs[i]:
            out.append(ts[i - 1] + (lev - vs[i - 1]) * (ts[i] - ts[i - 1]) / (vs[i] - vs[i - 1]))
    return out


def interp(ts, vs, t):
    import bisect
    i = bisect.bisect_left(ts, t)
    if i <= 0:
        return vs[0]
    if i >= len(ts):
        return vs[-1]
    a = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
    return vs[i - 1] + a * (vs[i] - vs[i - 1])


def main(ngspice, out):
    os.makedirs(out, exist_ok=True)
    s = open(REF).read()
    b = s[:s.index(".control")]
    w = lambda n, t: open(os.path.join(out, n), "w").write(t)
    # columns of every wrdata file below: t vout q il cx
    sigs = "v(vout) v(q) i(L1) v(cx)"
    w("r8_tran.cir", b + ".options reltol=1e-6\n.control\n%ssave %s\ntran 2n 300u 250u uic\nwrdata tran.txt %s\n.endc\n.end\n" % (CTL, sigs, sigs))
    card = "pss 395e3 100u q 400 10 150 1e-6 uic method=shoot trigger=cx trigval=%g trigdir=rise" % LEVEL
    for tag in ("hid", "nohid"):
        w("r8_pss_%s.cir" % tag, b + ".options reltol=1e-6\n.control\n%s%s\nwrdata orbit_%s.txt %s\n.endc\n.end\n" % (CTL, card, tag, sigs))
    w("r8_pac.cir", b.replace("Vrefsrc vrefnode 0 DC {vref}", "Vrefsrc vrefnode 0 DC {vref} AC 1") +
      '.options reltol=1e-4 klu\n.pac 395e3 100u q 400 10 150 1e-6 lin 3 20e3 80e3 method="shoot" solver="td" trigger="cx" trigval=%g trigdir="rise"\n'
      ".control\n%srun\nwrdata pac.txt v(vout)\n.endc\n.end\n" % (LEVEL, CTL))
    txt, wl = sh(ngspice, "r8_tran.cir", out, None, os.path.join(out, "r8_tran.log"))
    tr = rd(os.path.join(out, "tran.txt"))
    tt = [r[0] for r in tr]
    cxc = cross(tt, [r[4] for r in tr], LEVEL)
    per = [y - x for x, y in zip(cxc, cxc[1:])]
    tmean = sum(per) / len(per)
    t0, t1 = cxc[-2], cxc[-1]
    print("R8TRAN: settled tran reltol 1e-6: %d periods (rising crossings of cx through %g V), mean period %.6f us, last %.6f us, spread %.1e relative (script wall %.0f s)"
          % (len(per), LEVEL, tmean * 1e6, per[-1] * 1e6, (max(per) - min(per)) / tmean, wl), flush=True)
    for tag, env in (("hid", None), ("nohid", {"PSS_NO_HIDDEN": "1"})):
        txt, wl = sh(ngspice, "r8_pss_%s.cir" % tag, out, env, os.path.join(out, "r8_pss_%s.log" % tag))
        m = re.search(r"PSS-FREERUN: CONVERGED iter=(\d+) residual=(\S+) period=(\S+) fsw=(\S+) spectral_radius=(\S+) (\S+) wall_time=(\S+)s", txt)
        if not m:
            print("R8PSS %s: NO CONVERGED LINE" % tag, flush=True)
            continue
        T = float(m.group(3))
        print("R8PSS %s: iter %s, period %.6f us (settled tran mean %.6f us: %+.1e relative), %s, script wall %.0f s"
              % (tag, m.group(1), T * 1e6, tmean * 1e6, T / tmean - 1, "radius " + m.group(5) + " " + m.group(6), wl), flush=True)
        P = rd(os.path.join(out, "orbit_%s.txt" % tag))
        pt = [0.0] + [r[0] for r in P]
        n = 2000
        parts = []
        for c, nm, scale, unit in ((1, "vout", 1e3, "mV"), (3, "il", 1.0, "A"), (4, "cx", 1e3, "mV")):
            pv = [P[0][c]] + [r[c] for r in P]
            rv = [r[c] for r in tr]
            rp = [interp(tt, rv, t0 + (t1 - t0) * j / n) for j in range(n + 1)]
            pp = [interp(pt, pv, pt[-1] * j / n) for j in range(n + 1)]
            rip = max(rp) - min(rp)
            d = [a - b2 for a, b2 in zip(pp, rp)]
            parts.append("%s ripple %.3f %s, mean offset %+.3f %%, max pointwise %.3f %% of ripple" % (nm, rip * scale, unit, 100 * sum(d[:-1]) / n / rip, 100 * max(abs(x) for x in d) / rip))
        print("R8ORB %s: orbit vs the last tran period, phase 0 = the cx crossing: %s" % (tag, "; ".join(parts)), flush=True)
    txt, wl = sh(ngspice, "r8_pac.cir", out, None, os.path.join(out, "r8_pac.log"))
    gl = re.search(r"PAC-TD: generic charge-based history[^\n]*", txt)
    print("R8PAC: " + (gl.group(0)[:170] if gl else "no generic-history line") + " (script wall %.0f s)" % wl, flush=True)
    for ln in open(os.path.join(out, "pac.txt")).read().splitlines()[1:]:
        v = ln.split()
        if len(v) >= 3:
            f = round(float(v[0]))
            h = complex(float(v[1]), float(v[2]))
            db, dg = 20 * math.log10(abs(h)), math.degrees(cmath.phase(h))
            r = PP.get(f)
            print("R8PAC f=%d Hz: %+.4f dB / %+.3f deg   pac_perturbed --freerun (cited): %+.4f / %+.4f (halving %.3f)   difference %+.4f dB / %+.3f deg"
                  % (f, db, dg, r[0], r[1], r[2], db - r[0], dg - r[1]), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
