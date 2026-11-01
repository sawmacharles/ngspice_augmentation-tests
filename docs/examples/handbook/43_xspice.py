#!/usr/bin/env python3
"""43_xspice.py -- handbook example 43 (card #552): XSPICE code models under `.pss method=shoot` and td `.pac` (card #551, patches 0106-0108).

    python3 43_xspice.py NGSPICE OUTDIR

Small decks (copies of tests/fixtures/xspice_*.cir), one ngspice at a time, OMP_NUM_THREADS=2. NGSPICE is the build under test (s108 for card #552, s115 for card #560; its label is printed); the "before" build is
the series s105 build next to it (<apps>/aaba50ba3f49-s105, the series just before patch 0106) when present.
  1. td `.pac` of a stateless analog model (`gain`, gain=2.5 out_offset=0.1) against plain `.ac` of the same deck     (XSPAC lines)
  2. `.pss` orbit against the settled `.tran`, as % of the transient's ripple (max pointwise over one period):
       s_xfer (a stateful analog model), the d_tff divider (adc_bridge -> 2 x d_tff -> dac_bridge: a digital/hybrid chain, `driven=1`)    (XSPSS lines)
       and the red proof PSS_XS_NO_STATE=1 (no carried XSPICE state) on s_xfer and on the d_dff sample-and-hold deck (43_xspice_dig_dsample.cir)
  3. refusals, each by name or reason: td `.pac` on oneshot and on the digital deck with PAC_DIG_REFUSE=1 (it runs since card #714), `.pss` on an unlisted model (zener), a digital
     state that is not periodic with the analysis period (d_tff toggling at the clock period: refused at T, converges at 2 T), a free-running
     (trigger=) analysis with an XSPICE model                                                                        (XSREF lines)
"""
import bisect
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WR = "set filetype=ascii\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"


def run(ng, out, name, deck, env=None):
    p = os.path.join(out, name + ".cir")
    open(p, "w").write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2")
    e.update(env or {})
    r = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", p],
                       cwd=out, env=e, capture_output=True, text=True, timeout=1800)
    o = (r.stdout + r.stderr).replace("\r", "\n")
    open(os.path.join(out, name + ".log"), "w").write(o)
    return o, r.returncode


def fixture(name):
    t = open(os.path.join(HERE, name)).read()
    return re.sub(r"\n\.end\s*$", "\n", t.rstrip() + "\n")


def rd(path):
    rows = []
    for l in open(path):
        try:
            rows.append([float(x) for x in l.split()])
        except ValueError:
            pass
    return rows


def interp(ts, vs, t):
    i = bisect.bisect_left(ts, t)
    if i <= 0:
        return vs[0]
    if i >= len(ts):
        return vs[-1]
    a = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
    return vs[i - 1] + a * (vs[i] - vs[i - 1])


def pss_vs_tran(ng, out, fx, tag, f, vecs, pts, tstab, tstop, step, extra="", env=None, cols=(1,), tmin=1e-9):
    """(list of max pointwise % of ripple per column | None, pss log, converged, returncode)"""
    T = 1.0 / f
    base = fixture(fx) + ".options reltol=1e-5\n"
    for suf in ("pss", "tran"):
        p = os.path.join(out, "%s_%s.txt" % (tag, suf))
        if os.path.exists(p):
            os.remove(p)
    log, rc = run(ng, out, tag + "_pss", base + ".control\n" + WR + "pss %.15g %g %s %d 10 150 1e-6 uic method=shoot %s\nwrdata %s_pss.txt %s\nquit\n.endc\n.end\n"
                  % (f, tstab, "out", pts, extra, tag, vecs), env)
    if "PSS-SHOOT: CONVERGED" not in log:
        return None, log, False, rc
    run(ng, out, tag + "_tran", base + ".control\n" + WR + "tran %g %g uic\nwrdata %s_tran.txt %s\nquit\n.endc\n.end\n" % (step, tstop, tag, vecs))
    P, R = rd(os.path.join(out, tag + "_pss.txt")), rd(os.path.join(out, tag + "_tran.txt"))
    rt = [r[0] for r in R]
    t0 = rt[-1] - T
    pt = [0.0] + [r[0] for r in P]
    res = []
    N = 800
    for c in cols:
        pv = [P[0][c]] + [r[c] for r in P]
        rv = [r[c] for r in R]
        rp = [interp(rt, rv, t0 + T * j / N) for j in range(N + 1)]
        pp = [interp(pt, pv, pt[-1] * j / N) for j in range(N + 1)]
        rip = max(rp) - min(rp)
        res.append(100 * max(abs(a - b) for j, (a, b) in enumerate(zip(pp, rp)) if T * j / N >= tmin) / rip)
    return res, log, True, rc


def pac_gain(ng, out, tag):
    clk = "Vclk clk 0 PULSE(0 1 0 1n 1n 500u 1m)\nRclk clk clkgnd 100k\nCclk clkgnd 0 10n\n"
    body = fixture("43_xspice_an_gain.cir")
    body = "\n".join(l for l in body.split("\n") if not l.lower().startswith(("rl ", "co ", "vin "))) + "\nVin in 0 DC 0.3 AC 1\nRl o1 0 1k\n"
    deck = "* xspice pac\n" + clk + body + ".options reltol=1e-6 abstol=1e-15 vntol=1e-9\n"
    wr = ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata %s_%%s.txt o1\nquit\n.endc\n.end\n" % tag
    ta, _ = run(ng, out, tag + "_ac", deck + ".ac dec 5 1k 100k\n" + wr % "ac")
    tt, _ = run(ng, out, tag + "_td", deck + '.pac 1k 0 clkgnd 2000 1 8 1e-6 dec 5 1k 100k method="shoot" solver="td" pacpoints=100000\n' + wr % "td")

    def rdc(fn):
        r = []
        for l in open(os.path.join(out, fn)).read().split("\n")[1:]:
            t = l.split()
            if len(t) >= 3:
                r.append(complex(float(t[1]), float(t[2])))
        return r
    a, t = rdc(tag + "_ac.txt"), rdc(tag + "_td.txt")
    return max(abs(x - y) / abs(x) for x, y in zip(a, t)), abs(a[0]), abs(t[0])


def main(ng, out):
    os.makedirs(out, exist_ok=True)
    apps = os.path.dirname(os.path.dirname(os.path.dirname(ng)))
    old = os.path.join(apps, "aaba50ba3f49-s105", "bin", "ngspice")
    m = re.search(r"-(s\d+)(?:/|$)", os.path.dirname(os.path.dirname(ng)))
    cur = m.group(1) if m else "current"   # the label of the build under test (s115 for card #560)
    builds = [(cur, ng)] + ([("s105", old)] if os.path.exists(old) else [])
    for tag, b in builds:
        w, a0, t0 = pac_gain(b, out, "pg_" + tag)
        print("XSPAC %s: td .pac vs .ac of gain(gain=2.5 out_offset=0.1), 5 points 1k..100k: worst relative difference %.2e (|H| at 1 kHz: .ac %.6f, td %.6f)" % (tag, w, a0, t0), flush=True)
    for tag, b in builds:
        e, log, ok, rc = pss_vs_tran(b, out, "43_xspice_an_s_xfer.cir", "sx_" + tag, 100e3, "v(out) v(o1)", 400, 4e-4, 6e-4, 2e-9, cols=(1, 2))
        print("XSPSS %s s_xfer: %s" % (tag, ("orbit vs settled tran: vout %.3f %%, o1 %.3f %% of ripple (max pointwise)" % tuple(e)) if ok else "no convergence (rc %d)" % rc), flush=True)
    e, log, ok, rc = pss_vs_tran(ng, out, "43_xspice_an_s_xfer.cir", "sx_red", 100e3, "v(out) v(o1)", 400, 4e-4, 6e-4, 2e-9, env={"PSS_XS_NO_STATE": "1"}, cols=(1, 2))
    print("XSPSS %s s_xfer PSS_XS_NO_STATE=1: %s" % (cur, ("orbit vs settled tran: vout %.3f %%, o1 %.3f %%" % tuple(e)) if ok else "no convergence"), flush=True)
    for tag, b in builds:
        e, log, ok, rc = pss_vs_tran(b, out, "43_xspice_dig_tff.cir", "tf_" + tag, 250e3, "v(out)", 4000, 40e-6, 240e-6, 0.25e-9, extra="driven=1")
        hs = re.search(r"XSPICE hidden state[^\n]*", log)
        crash = ("signal/segfault (rc %d)" % rc) if (not ok and (rc < 0 or rc > 128 or "Segmentation" in log or "SIGSEGV" in log)) else "no convergence (rc %d)" % rc
        print("XSPSS %s d_tff divider by 4 (adc_bridge, 2 x d_tff, dac_bridge; 1 MHz clock divided by 4 = 250 kHz, Points 4000): %s%s"
              % (tag, ("orbit vs settled tran %.4f %% of ripple" % e[0]) if ok else crash, ("; " + hs.group(0)[:150]) if hs else ""), flush=True)
    e, log, ok, rc = pss_vs_tran(ng, out, "43_xspice_dig_dsample.cir", "ds", 250e3, "v(out)", 4000, 40e-6, 240e-6, 0.25e-9, extra="driven=1")
    e2, log2, ok2, rc2 = pss_vs_tran(ng, out, "43_xspice_dig_dsample.cir", "ds_red", 250e3, "v(out)", 4000, 40e-6, 240e-6, 0.25e-9, extra="driven=1", env={"PSS_XS_NO_STATE": "1"})
    print("XSPSS %s d_dff sample-and-hold (comparator + d_dff of a 250 kHz sine at 1 MHz): default %s; PSS_XS_NO_STATE=1 %s"
          % (cur, ("%.4f %% of ripple" % e[0]) if ok else "no convergence", ("%.2f %%" % e2[0]) if ok2 else "no convergence"), flush=True)
    # ---- refusals (s108)
    clk = fixture("43_xspice_dig_tff.cir")
    sxdeck = ("* xspice pac\nVclk clk 0 PULSE(0 1 0 1n 1n 500u 1m)\nRclk clk clkgnd 100k\nCclk clkgnd 0 10n\n"
              + "\n".join(l for l in fixture("43_xspice_an_s_xfer.cir").split("\n") if not l.lower().startswith(("rl ", "co ", "vin ")))
              + "\nVin in 0 DC 0.3 AC 1\nRl o1 0 1k\n.options reltol=1e-6\n")
    # card #557: td .pac/.pstb/.pxf now run THROUGH s_xfer / int / d_dt (the model's state slots are Floquet states); oneshot (event) is still refused
    log, rc = run(ng, out, "ref_pac_sx", sxdeck + '.pac 1k 0 clkgnd 2000 1 8 1e-6 dec 3 1k 100k method="shoot" solver="td" pacpoints=100000\n.control\nrun\nquit\n.endc\n.end\n')
    hs = re.search(r"stateful XSPICE code models \(card #557\): (\d+) hidden states", log)
    print("XSTD .pac on s_xfer: " + (("runs, %s hidden states" % hs.group(1)) if hs and "PAC-TD: done" in log else "NOT RUN: " + "; ".join(l.strip()[:120] for l in log.splitlines() if "REFUSED" in l)), flush=True)
    osdeck = ("* xspice pac oneshot\nVclk clk 0 PULSE(0 1 0 1n 1n 500u 1m)\nRclk clk clkgnd 100k\nCclk clkgnd 0 10n\n"
              + "\n".join(l for l in fixture("43_xspice_an_s_xfer.cir").split("\n") if not l.lower().startswith(("rl ", "co ", "vin ", "a1 ", ".model")))
              + "\nVin in 0 DC 0.3 AC 1\nVc c 0 DC 1.0\nVclr z 0 DC 0\nA1 in c z o1 m\n"
              ".model m oneshot(clk_trig=0.5 pos_edge_trig=TRUE cntl_array=[0 1] pw_array=[1u 2u] out_low=0 out_high=1)\nRl o1 0 1k\n.options reltol=1e-6\n")
    log, rc = run(ng, out, "ref_pac_os", osdeck + '.pac 1k 0 clkgnd 2000 1 8 1e-6 dec 3 1k 100k method="shoot" solver="td"\n.control\nrun\nquit\n.endc\n.end\n')
    l1 = [l for l in log.splitlines() if "REFUSED" in l]
    print("XSREF .pac on oneshot: " + (l1[0].strip()[:260] if l1 else "NOT REFUSED"), flush=True)
    # card #714 (s150): td .pac runs through digital blocks (edge shifts); PAC_DIG_REFUSE=1 restores the s149 refusal
    digdeck = clk + '.options reltol=1e-5\n.pac 250k 40u out 400 1 8 1e-6 dec 3 1k 10k method="shoot" solver="td" driven=1\n.control\nrun\nquit\n.endc\n.end\n'
    log, rc = run(ng, out, "ref_pac_dig", digdeck)
    l2 = [l for l in log.splitlines() if "REFUSED" in l]
    print("XSTD .pac on the digital deck: " + ("runs" if "PAC-TD: done" in log and not l2 else "NOT RUN: " + "; ".join(l.strip()[:120] for l in l2)), flush=True)
    log, rc = run(ng, out, "ref_pac_dig_refuse", digdeck, env={"PAC_DIG_REFUSE": "1"})
    l2 = [l for l in log.splitlines() if "REFUSED" in l]
    print("XSREF .pac on the digital deck with PAC_DIG_REFUSE=1: " + (l2[0].strip()[:260] if l2 else "NOT REFUSED"), flush=True)
    zd = ("* zener\nVin in 0 SIN(0 2 100k)\nRs in o1 100\nA1 %gd(o1 0) m\n.model m zener(v_breakdown=1.5 i_breakdown=1m r_breakdown=10 i_sat=1e-12 n_forward=1)\n"
          "Co o1 0 100p\n.options reltol=1e-5\n.control\npss 100e3 40e-6 o1 400 10 150 1e-6 uic method=shoot\nquit\n.endc\n.end\n")
    log, rc = run(ng, out, "ref_zener", zd)
    l3 = [l for l in log.splitlines() if "REFUSED" in l]
    print("XSREF .pss on an unlisted model (zener): " + (l3[0].strip()[:260] if l3 else "NOT REFUSED"), flush=True)
    for f in (1e6, 5e5):
        e, log, ok, rc = pss_vs_tran(ng, out, "43_xspice_dig_toggle.cir", "tog_%d" % f, f, "v(out)", 4000, 40e-6, 240e-6, 0.25e-9, extra="driven=1")
        l4 = [l for l in log.splitlines() if "not periodic with the analysis period" in l]
        print("XSREF d_tff toggling at the clock period, .pss at %g MHz: %s" % (f / 1e6, ("converges, orbit vs settled tran %.4f %%" % e[0]) if ok else ("REFUSED: " + (l4[0].strip()[:200] if l4 else "no convergence"))), flush=True)
    log, rc = run(ng, out, "ref_free", clk + ".options reltol=1e-5\n.control\npss 250k 40u out 400 10 150 1e-6 uic method=shoot trigger=out trigval=0.5 trigdir=rise\nquit\n.endc\n.end\n")
    l5 = [l for l in log.splitlines() if "free-running (trigger=)" in l]
    print("XSREF free-running (trigger=) with an XSPICE model: " + (l5[0].strip()[:260] if l5 else "NOT REFUSED"), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
