#!/usr/bin/env python3
"""Reproducer for bug 1 (card #539): .pac/.pstb/.pxf/.pnoise solver=td lose the history of device-model charge storage.

A chain of N MOS inverters (level 1, optionally with overlap caps cgso/cgdo; also BSIM4 below) is driven by ramp+Vc
(period 2 us). Its low-frequency control-to-average-output gain is the static derivative d<v>/dVc of the periodic
steady state, which two plain .tran runs give (central difference). A correct .pac solver=td must agree at 1 kHz.

usage: NGSPICE=<ngspice-aug> python3 run_repro.py      (writes only under ./work)
"""
import os, re, subprocess, sys
NG = os.environ.get("NGSPICE") or os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s88/bin/ngspice-aug")
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(HERE, "work"); os.makedirs(W, exist_ok=True)
MODELS = {
 "level1_cg0":  ".model nch4 nmos level=1 vto=0.5 kp=3e-4 lambda=0.05\n.model pch4 pmos level=1 vto=-0.5 kp=1e-4 lambda=0.05\n",
 "level1_cg1e-8": ".model nch4 nmos level=1 vto=0.5 kp=3e-4 lambda=0.05 cgso=1e-8 cgdo=1e-8\n.model pch4 pmos level=1 vto=-0.5 kp=1e-4 lambda=0.05 cgso=1e-8 cgdo=1e-8\n",
 "bsim4": """.model nch4 nmos level=54 version=4.8.2
+ mobmod=0 capmod=0 igcmod=0 igbmod=0 rdsmod=0 rbodymod=0 rgatemod=0 acnqsmod=0
+ toxe=3e-9 toxp=3e-9 toxm=3e-9
+ vth0=0.5 k1=0.4 k2=0 k3=0 w0=2.5e-6
+ u0=300 vsat=1.1e5
+ cgso=2e-10 cgdo=2e-10 cgbo=0 cgsl=0 cgdl=0
+ cf=0 cjs=1e-4 cjd=1e-4 cjswd=1e-10 cjswgd=1e-10
.model pch4 pmos level=54 version=4.8.2
+ mobmod=0 capmod=0 igcmod=0 igbmod=0 rdsmod=0 rbodymod=0 rgatemod=0 acnqsmod=0
+ toxe=3e-9 toxp=3e-9 toxm=3e-9
+ vth0=-0.5 k1=0.4 k2=0 k3=0 w0=2.5e-6
+ u0=100 vsat=1.0e5
+ cgso=2e-10 cgdo=2e-10 cgbo=0 cgsl=0 cgdl=0
+ cf=0 cjs=1e-4 cjd=1e-4 cjswd=1e-10 cjswgd=1e-10
""",
}
def deck(stages, models, vc="2.5 AC 1", analysis="", extra=""):
    d = "* bug1 reproducer: %d inverters\n" % stages + models + "Vdd vdd 0 5\nVc ctl 0 DC %s\nVr inp ctl PULSE(-2.5 2.5 0 0.9u 0.1u 1n 2u)\n" % vc
    prev = "inp"
    for i in range(stages):
        nxt = "n%d" % i if i < stages - 1 else "out"
        d += ("Mn%d %s %s 0 0 nch4 w=4u l=0.5u as=4p ad=4p ps=10u pd=10u\n"
              "Mp%d %s %s vdd vdd pch4 w=10u l=0.5u as=10p ad=10p ps=22u pd=22u\n") % (i, nxt, prev, i, nxt, prev)
        prev = nxt
    return d + "Rf out vo 1k\nCf vo 0 1n\n" + analysis
def run(name, text, env=None):
    p = os.path.join(W, name + ".cir"); open(p, "w").write(text)
    e = dict(os.environ, OMP_NUM_THREADS="2"); e.update(env or {})
    return subprocess.run([NG, "-b", p], capture_output=True, text=True, env=e, cwd=W).stdout
def static_gain(stages, models):
    av = []
    for vc in (2.51, 2.49):
        o = run("st", deck(stages, models, vc="%g" % vc, analysis=".control\ntran 1n 60u 0 uic\nmeas tran av AVG v(out) FROM=40u TO=60u\n.endc\n.end\n"))
        av.append(float(re.search(r"^av\s*=\s*([-0-9.e+]+)", o, re.M).group(1)))
    return (av[0] - av[1]) / 0.02
def pac_gain(stages, models, pts=400, env=None):
    an = (".pac 500e3 20e-6 vo 400 10 50 1e-7 dec 1 1000 1000 method=\"shoot\" solver=\"td\" pacpoints=%d\n.options reltol=1e-4\n"
          ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\nrun\nwrdata pac.txt v(out)\n.endc\n.end\n") % pts
    run("pac", deck(stages, models, analysis=an), env)
    r = open(os.path.join(W, "pac.txt")).read().split("\n")[1].split()
    return complex(float(r[1]), float(r[2]))
if __name__ == "__main__":
    print("%-14s %-7s %9s %12s %12s" % ("models", "stages", "static", "pac(400)", "pac(40000)"))
    for mname in ("level1_cg0", "level1_cg1e-8", "bsim4"):
        for st in (1, 2, 3):
            s = static_gain(st, MODELS[mname])
            a = pac_gain(st, MODELS[mname], 400); b = pac_gain(st, MODELS[mname], 40000)
            print("%-14s %-7d %9.4f %12.4f %12.4f" % (mname, st, abs(s), abs(a), abs(b)), flush=True)
