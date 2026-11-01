#!/usr/bin/env python3
"""44_r9.py -- handbook example 44 (card #552): R9, the digitally controlled synchronous buck made of XSPICE primitives (card #551).

    python3 44_r9.py NGSPICE OUTDIR

Deck: refs/r9_digital_buck.cir (an 11-bit up/down counter integrator, adc_bridge comparators with a zero-error bin, d_pwm; R3's power stage with
smooth tanh switches) with its own `.control` cut off and only the analysis added. One ngspice at a time, OMP_NUM_THREADS=2. Runs (minutes):
  T  settled `.tran` 60 ms (100 ns grid): the DPWM code is constant from 6 ms; the orbit reference is its last period, 58.8..60 ms (T = 1.2 ms = 600 cycles)
  P  `.pss 833.33 Hz` (T = 1.2 ms, the enable-clock period of the integrator), stabilisation 12 ms, Points 24000, method=shoot
  Q  td `.pac` on the same deck: the refusal (a quantiser has no small-signal gain)
  L  the 7-bit variant (tools/gen_r9.py --nbits 7: DPWM step 92 mV against a 6 mV zero-error bin): its settled `.tran` (30 ms) hunts with period 2 T
     (window means alternate), and `.pss` at T is REFUSED
Prints lines R9TRAN / R9PSS / R9PAC / R97BIT. The loop-gain FRA of this deck (amplitude-dependent by up to 100 dB, docs/pss-pac.md card #551) is
cited on the page, not rerun: it needs 24 ms of settling per point.
"""
import bisect
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CTL = "set filetype=ascii\nset wr_singlescale\nset wr_vecnames\nset numdgt=12\n"
T = 1.2e-3


def sh(ng, deck, cwd, log):
    e = dict(os.environ, OMP_NUM_THREADS="2")
    t0 = time.time()
    p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", deck],
                       cwd=cwd, env=e, capture_output=True, text=True, timeout=3600)
    o = (p.stdout + p.stderr).replace("\r", "\n")
    open(log, "w").write(o)
    return o, time.time() - t0


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


def body(path):
    s = open(path).read()
    b = s[:s.index(".control")]
    return b.replace("../stb/lgprobe.lib", os.path.join(ROOT, "stb", "lgprobe.lib"))


def main(ng, out):
    os.makedirs(out, exist_ok=True)
    b = body(os.path.join(ROOT, "refs", "r9_digital_buck.cir"))
    w = lambda n, t: open(os.path.join(out, n), "w").write(t)
    pss = "pss %.15g 12m vout 24000 10 150 1e-6 uic method=shoot" % (1 / T)
    w("r9_tran.cir", b + ".options reltol=1e-4\n.control\n%stran 100n 60m uic\nwrdata r9_tran.txt v(vout) v(code)\n.endc\n.end\n" % CTL)
    w("r9_pss.cir", b + ".options reltol=1e-4\n.control\n%s%s\nwrdata r9_pss.txt v(vout) v(code)\n.endc\n.end\n" % (CTL, pss))
    w("r9_pac.cir", b.replace("XP fbtop vout lgprobe", "Vprobe fbtop vout DC 0") + ".options reltol=1e-4\n"
      '.pstb %.15g 12m vout 24000 10 150 1e-6 lin 2 300 1000 probe="Vprobe" method="shoot"\n.control\nrun\n.endc\n.end\n' % (1 / T))
    o, wl = sh(ng, "r9_tran.cir", out, os.path.join(out, "r9_tran.log"))
    R = rd(os.path.join(out, "r9_tran.txt"))
    rt = [r[0] for r in R]
    codes_t = sorted({round(r[2] * 2048) for r in R if r[0] > 6e-3})
    print("R9TRAN: settled tran 60 ms: DPWM code %s over 6..60 ms, script wall %.0f s" % (codes_t, wl), flush=True)
    o, wp = sh(ng, "r9_pss.cir", out, os.path.join(out, "r9_pss.log"))
    m = re.search(r"PSS-SHOOT: CONVERGED iter=(\d+) residual=(\S+) spectral_radius=(\S+) (\S+) wall_time=(\S+)s", o)
    hs = re.search(r"XSPICE hidden state over one period[^\n]*", o)
    if not m:
        print("R9PSS: NO CONVERGED LINE: " + " | ".join(l for l in o.splitlines() if "FAILED" in l or "REFUSED" in l)[:300], flush=True)
    else:
        P = rd(os.path.join(out, "r9_pss.txt"))
        pt = [0.0] + [r[0] for r in P]
        pv = [P[0][1]] + [r[1] for r in P]
        rv = [r[1] for r in R]
        t0 = rt[-1] - T
        n = 3000
        c = [interp(rt, rv, t0 + T * j / n) for j in range(n + 1)]
        p = [interp(pt, pv, pt[-1] * j / n) for j in range(n + 1)]
        rip = max(c) - min(c)
        d = max(abs(a - x) for a, x in zip(c, p))
        mo = sum(a - x for a, x in zip(p, c)) / len(c)
        codes = sorted({round(r[2] * 2048) for r in P})
        print("R9PSS: iter %s, residual %s, radius %s %s, %ss; orbit vs the last tran period (58.8..60 ms): ripple %.3f mV, max pointwise %.4f mV = %.2f %% of ripple, "
              "mean offset %+.4f mV; DPWM code %s; %s (script wall %.0f s)"
              % (m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), rip * 1e3, d * 1e3, 100 * d / rip, mo * 1e3, codes,
                 hs.group(0)[:260] if hs else "no hidden-state line", wp), flush=True)
    o, wq = sh(ng, "r9_pac.cir", out, os.path.join(out, "r9_pac.log"))
    ln = [l for l in o.splitlines() if "REFUSED" in l]
    print("R9PAC: " + (ln[0].strip()[:300] if ln else "NOT REFUSED") + " (script wall %.0f s)" % wq, flush=True)
    # ---- the limit-cycling 7-bit variant
    subprocess.run([sys.executable, os.path.join(ROOT, "tools", "gen_r9.py"), "--nbits", "7", "--out", os.path.join(out, "r9_7bit_gen.cir"),
                    "--lgprobe", os.path.join(ROOT, "stb", "lgprobe.lib")], check=True, capture_output=True)
    b7 = body(os.path.join(out, "r9_7bit_gen.cir"))
    w("r9_7bit_tran.cir", b7 + ".options reltol=1e-4\n.control\n%stran 100n 30m uic\nwrdata r9_7bit_tran.txt v(vout) v(code)\n.endc\n.end\n" % CTL)
    w("r9_7bit_pss.cir", b7 + ".options reltol=1e-4\n.control\n%spss %.15g 12m vout 6000 10 6 1e-6 uic method=shoot\nwrdata r9_7bit_pss.txt v(vout) v(code)\n.endc\n.end\n" % (CTL, 1 / T))
    o, wl = sh(ng, "r9_7bit_tran.cir", out, os.path.join(out, "r9_7bit_tran.log"))
    R7 = rd(os.path.join(out, "r9_7bit_tran.txt"))
    wins = {}
    for r in R7:
        if r[0] > 12e-3:
            wins.setdefault(int(r[0] / T), []).append(r)
    means = [sum(x[1] for x in v) / len(v) for k, v in sorted(wins.items())]
    codes = sorted({round(r[2] * 128) for r in R7 if r[0] > 12e-3})
    alt = max(abs(means[i] - means[i + 1]) for i in range(len(means) - 1))
    print("R97BIT tran: window means of Vout (12..30 ms, windows of T = 1.2 ms): %s V; consecutive windows differ by up to %.0f mV; DPWM codes visited %s (script wall %.0f s)"
          % (" ".join("%.3f" % x for x in means[:8]), alt * 1e3, codes, wl), flush=True)
    o, wp = sh(ng, "r9_7bit_pss.cir", out, os.path.join(out, "r9_7bit_pss.log"))
    ln = [l for l in o.splitlines() if "REFUSED" in l or "FAILED" in l or "not periodic" in l]
    print("R97BIT pss at T: %s%s (script wall %.0f s)" % ("CONVERGED?! " if "PSS-SHOOT: CONVERGED" in o else "", (ln[0].strip()[:260] if ln else "no refusal line"), wp), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
