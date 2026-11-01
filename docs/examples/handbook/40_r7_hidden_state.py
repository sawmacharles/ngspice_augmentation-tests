#!/usr/bin/env python3
"""40_r7_hidden_state.py -- handbook example 40 (card #552): the transistor-level buck R7 after cards #546 and #549.

    python3 40_r7_hidden_state.py NGSPICE OUTDIR

Decks are the ones of example 32 (32_r7_decks.py: refs/r7_bsim_buck.cir with only the analysis card, `.options reltol=1e-4` and the plain
`Vprobe` changed). One ngspice at a time, OMP_NUM_THREADS=2. Runs:

  A. `.pss method=shoot` (reltol 1e-4)                          default (hidden-state continuation of the period restart, card #546)
  B. the same with PSS_NO_HIDDEN=1                              the pre-#546 restart from the native C/L states only ("before")
  C. a settled `.tran` of the same circuit, reltol 1e-5, 1 ns grid, 300 us (150 periods; the slowest decay is 0.9 per period): the reference.
     Orbit = last period (298..300 us); phase 0 = clock phase 0 for the orbit and for the transient, so no alignment is needed.
  D. `.pstb` (500 Hz grid) default / legacy (PAC_RECORD_RELTOL_AUTO=0, the pre-#625 default) / `pacreltol=1e-7` (card #546) / `pacreltol=1e-7` with `.options method=gear` (card #549), each against
     the converged FRA loop gain of ~/Work/data/ngspice_augmentation/fra/curves-r7-20261001/loop/fra.csv (not regenerated here: 16-24 min a point)
     and the gear run against the trapezoidal one.

Prints lines starting R7ORB / R7PSTB / R7GEAR (parsed by tests/test_handbook.py, which recomputes the comparisons itself).
"""
import csv
import importlib.util
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("d32", os.path.join(HERE, "32_r7_decks.py"))
d32 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d32)

FRA = os.path.expanduser("~/Work/data/ngspice_augmentation/fra/curves-r7-20261001/loop/fra.csv")
T = 2e-6


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


def interp(ts, vs, t):
    import bisect
    i = bisect.bisect_left(ts, t)
    if i <= 0:
        return vs[0]
    if i >= len(ts):
        return vs[-1]
    a = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
    return vs[i - 1] + a * (vs[i] - vs[i - 1])


def orbit_vs_tran(orb, ref):
    """orb rows (t, vout, il) over one period starting at 0; ref rows (t, vout, il) = the last period of the settled tran (t - 298 us).
    Returns {signal: (mean offset in mV or mA, % of ripple, max pointwise % of ripple for t >= 2 ns)}."""
    rt = [r[0] for r in ref]
    ot = [r[0] for r in orb]
    out = {}
    for c, nm, unit in ((1, "vout", 1e3), (2, "il", 1e3)):
        rv = [r[c] for r in ref]
        ov = [r[c] for r in orb]
        n = 2000
        ripple = max(rv) - min(rv)
        d = [interp(ot, ov, (j + 0.5) * T / n) - interp(rt, rv, (j + 0.5) * T / n) for j in range(n)]
        dpt = [abs(x) for j, x in enumerate(d) if (j + 0.5) * T / n >= 2e-9]
        out[nm] = (sum(d) / n * unit, 100 * sum(d) / n / ripple, 100 * max(dpt) / ripple, ripple)
    return out


def fra_loop():
    if not os.path.exists(FRA):
        return {}
    r = {}
    for row in csv.DictReader(open(FRA)):
        if row["status"] == "ok":
            r[round(float(row["freq"]))] = (float(row["combined_db"]), float(row["combined_phase_deg"]), float(row["unc_db"]))
    return r


def pstb_vs_fra(path, fra):
    import cmath
    out = {}
    for ln in open(path).read().splitlines()[1:]:
        v = ln.split()
        if len(v) < 3:
            continue
        f = round(float(v[0]))
        if f in fra:
            h = complex(float(v[1]), float(v[2]))
            out[f] = (20 * math.log10(abs(h)) - fra[f][0], (math.degrees(cmath.phase(h)) - fra[f][1] + 180) % 360 - 180)
    return out


def pstb_vs_pstb(a, b):
    import cmath
    ra = [[float(x) for x in ln.split()] for ln in open(a).read().splitlines()[1:] if ln.strip()]
    rb = [[float(x) for x in ln.split()] for ln in open(b).read().splitlines()[1:] if ln.strip()]
    wd = wp = 0.0
    for x, y in zip(ra, rb):
        h1, h2 = complex(x[1], x[2]), complex(y[1], y[2])
        wd = max(wd, abs(20 * math.log10(abs(h1) / abs(h2))))
        wp = max(wp, abs((math.degrees(cmath.phase(h1)) - math.degrees(cmath.phase(h2)) + 180) % 360 - 180))
    return wd, wp


def main(ngspice, out):
    os.makedirs(out, exist_ok=True)
    b = d32.body()
    pb = d32.probe_body(b)
    pss = ".pss " + d32.PSS.format(iters=150) + " uic method=shoot"
    # --- A, B: .pss with and without the hidden-state continuation
    for tag, env in (("hid", None), ("nohid", {"PSS_NO_HIDDEN": "1"})):
        open(os.path.join(out, "r7_pss_%s.cir" % tag), "w").write(d32.deck(b, "reltol=1e-4", pss, "wrdata orbit_%s.txt v(vout) i(L1)" % tag))
        txt, w = sh(ngspice, "r7_pss_%s.cir" % tag, out, env, os.path.join(out, "r7_pss_%s.log" % tag))
        m = re.search(r"PSS-SHOOT: CONVERGED iter=(\d+) residual=(\S+) spectral_radius=(\S+) (\S+) wall_time=(\S+)s", txt)
        print("R7PSS %s: %s (script wall %.0f s)" % (tag, m.group(0) if m else "NO CONVERGED LINE", w), flush=True)
    # --- C: settled transient
    open(os.path.join(out, "r7_tran.cir"), "w").write(
        b + ".options reltol=1e-5\n.control\nset filetype=ascii\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
        "tran 1n 300u 296u uic\nlinearize v(vout) i(L1)\nwrdata tran.txt v(vout) i(L1)\n.endc\n.end\n")
    txt, w = sh(ngspice, "r7_tran.cir", out, None, os.path.join(out, "r7_tran.log"))
    tr = rd(os.path.join(out, "tran.txt"))
    ref = [[r[0] - 298e-6] + r[1:] for r in tr if 298e-6 - 1e-12 <= r[0] <= 300e-6 + 1e-12]
    prev = [[r[0] - 296e-6] + r[1:] for r in tr if 296e-6 - 1e-12 <= r[0] <= 298e-6 + 1e-12]
    rv, pv = [r[1] for r in ref], [r[1] for r in prev]
    per = max(abs(interp([p[0] for p in prev], pv, r[0]) - r[1]) for r in ref) / (max(rv) - min(rv))
    print("R7TRAN: settled tran reltol 1e-5, %d rows in the last period, vout mean %.5f V, ripple %.2f mV, last period vs the one before %.1e of ripple (script wall %.0f s)"
          % (len(ref), sum(rv) / len(rv), (max(rv) - min(rv)) * 1e3, per, w), flush=True)
    for tag in ("nohid", "hid"):
        orb = rd(os.path.join(out, "orbit_%s.txt" % tag))
        r = orbit_vs_tran(orb, ref)
        print("R7ORB %s: vout mean offset %+.4f mV (%+.3f %% of ripple), max pointwise %.4f %%; il mean offset %+.4f mA (%+.3f %%), max pointwise %.4f %%"
              % (tag, r["vout"][0], r["vout"][1], r["vout"][2], r["il"][0], r["il"][1], r["il"][2]), flush=True)
    # --- D: .pstb, recording-pass tolerance and integration method
    fra = fra_loop()
    card = '.pstb ' + d32.PSS.format(iters=150) + ' ' + d32.GRID + ' probe="Vprobe" method="shoot"'
    # card #625: "default" records at reltol 1e-8 (switch-free circuit); "legacy" is the pre-#625 default (recording at the inherited reltol)
    variants = (("default", "", "reltol=1e-4", None), ("legacy", "", "reltol=1e-4", {"PAC_RECORD_RELTOL_AUTO": "0"}),
                ("pacreltol", " pacreltol=1e-7", "reltol=1e-4", None), ("gear", " pacreltol=1e-7", "reltol=1e-4\n.options method=gear", None))
    for tag, extra, opts, env in variants:
        open(os.path.join(out, "r7_pstb_%s.cir" % tag), "w").write(d32.deck(pb, opts, card + extra, "wrdata pstb_%s.txt t" % tag))
        txt, w = sh(ngspice, "r7_pstb_%s.cir" % tag, out, env, os.path.join(out, "r7_pstb_%s.log" % tag))
        pm = re.search(r"pstb: phase margin [\d.]+ deg at [\d.eE+]+ Hz", txt)
        eng = re.search(r"PSTB: done, wall_time=([\d.]+)s", txt)
        d = pstb_vs_fra(os.path.join(out, "pstb_%s.txt" % tag), fra)
        sel = " ".join("%gk:%+.3f" % (f / 1e3, d[f][0]) for f in (20000, 40000, 80000, 200000) if f in d)
        wdb = max((abs(v[0]) for v in d.values()), default=float("nan"))
        print("R7PSTB %s: %s; engine-FRA dB at %s; worst over %d FRA frequencies %.3f dB; recording pass %ss, script wall %.0f s"
              % (tag, pm.group(0) if pm else "no PM line", sel, len(d), wdb, eng.group(1) if eng else "?", w), flush=True)
    wd, wp = pstb_vs_pstb(os.path.join(out, "pstb_gear.txt"), os.path.join(out, "pstb_pacreltol.txt"))
    print("R7GEAR: gear vs trapezoidal (both pacreltol=1e-7), all 400 frequencies: worst %.4f dB / %.4f deg" % (wd, wp), flush=True)
    # the refusal of gear in .pnoise (message only: the deck stops at the card)
    pn = d32.deck(b, "reltol=1e-4\n.options method=gear", '.pnoise ' + d32.PSS.format(iters=150) + ' vout Vin lin 2 1000 2000 method="shoot" solver="td"',
                  "wrdata pn_gear.txt onoise_spectrum")
    open(os.path.join(out, "r7_pnoise_gear.cir"), "w").write(pn)
    txt, w = sh(ngspice, "r7_pnoise_gear.cir", out, None, os.path.join(out, "r7_pnoise_gear.log"))
    ln = [l for l in txt.splitlines() if "method=gear" in l]
    print("R7GEAR pnoise: " + (ln[0].strip()[:200] if ln else "NO gear refusal line"), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
