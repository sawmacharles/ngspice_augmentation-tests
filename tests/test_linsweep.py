#!/usr/bin/env python3
"""Card #651: a linear frequency sweep of 2 points gives both points (.ac, .sp, .stb).

Checks (fork binary from $NGSPICE; stock comparison from $NGSPICE_STOCK, default s137):
  a. `lin 2 f1 f2` gives exactly [f1, f2] in .ac, .sp and .stb.
  b. lin 1, lin 3, lin 7, dec 5 and oct 3 give the same frequencies as the stock build (only N = 2 changed),
     and the stock build gives one point for lin 2 (the test can fail: the comparator sees the old behaviour).
  c. .noise lin 2 (which was already right) is unchanged and has 2 points.

    NGSPICE=<build>/bin/ngspice python3 tests/test_linsweep.py     exit 0 pass, 1 any FAIL, 2 could-not-run
Outputs go to ~/Work/data/ngspice_augmentation/651, never the repo.
"""
import os
import subprocess
import sys

APPS = os.path.expanduser("~/Work/apps/ngspice_augmentation")
NGSPICE = os.environ.get("NGSPICE", os.path.join(APPS, "aaba50ba3f49-s152/bin/ngspice"))
STOCK = os.environ.get("NGSPICE_STOCK", os.path.join(APPS, "aaba50ba3f49-s137/bin/ngspice"))
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/651/test")
RESULTS = []

CIRCUITS = {
    "ac": "v1 a 0 dc 0 ac 1\nr1 a b 1k\nc1 b 0 1n\n.ac {sweep}\n",
    "sp": "v1 a 0 dc 0 ac 1 portnum 1 z0 50\nr1 a b 75\nr2 b 0 40\nv2 b 0 dc 0 ac 0 portnum 2 z0 50\n.sp {sweep}\n",
    "stb": ("vdd vdd 0 5\nvin inp 0 dc 1\ne1 out 0 inp fb 1000\nvprobe out fbx 0\nr1 fbx fb 10k\nr2 fb 0 10k\n"
            "c1 fbx 0 1n\n.stb {sweep} probe=vprobe\n"),
    "noise": "v1 a 0 dc 0 ac 1\nr1 a b 1k\nr2 b 0 1k\n.noise v(b) v1 {sweep}\n",
}


def check(name, ok, detail):
    RESULTS.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


def freqs(binary, kind, sweep, tag):
    os.makedirs(DATA, exist_ok=True)
    deck = os.path.join(DATA, f"{tag}_{kind}_{sweep.replace(' ', '_')}.cir")
    body = CIRCUITS[kind].format(sweep=sweep)
    ctl = ".control\nrun\n" + ("setplot noise1\n" if kind == "noise" else "setplot stb1\n" if kind == "stb" else "") + "print frequency\n.endc\n.end\n"
    with open(deck, "w") as fh:
        fh.write(f"lin sweep {kind}\n" + body + ctl)
    p = subprocess.run([binary, "-b", deck], capture_output=True, text=True, timeout=120)
    out = p.stdout + p.stderr
    vals = []
    single = None
    for ln in out.splitlines():
        t = ln.split()
        if ln.startswith("frequency = "):
            single = float(ln.split("=")[1].split(",")[0])
        elif len(t) >= 2 and t[0].isdigit():
            try:
                vals.append(float(t[1].split(",")[0]))
            except ValueError:
                pass
    if single is not None and not vals:
        vals = [single]
    return vals


def main():
    for b in (NGSPICE, STOCK):
        if not os.path.exists(b):
            print(f"COULD-NOT-RUN: {b} missing")
            return 2
    for kind in ("ac", "sp", "stb"):
        f = freqs(NGSPICE, kind, "lin 2 1k 1meg", "new")
        check(f"a {kind} lin 2", f == [1e3, 1e6], f"{f}")
    for kind in ("ac", "sp", "stb"):
        same = []
        for sweep in ("lin 1 1k 1meg", "lin 3 1k 1meg", "lin 7 10 1meg", "dec 5 1k 1meg", "oct 3 1k 64k"):
            a, b = freqs(NGSPICE, kind, sweep, "new"), freqs(STOCK, kind, sweep, "stock")
            same.append((sweep, a == b and len(a) > 0, len(a)))
        check(f"b {kind} other sweeps as stock", all(s[1] for s in same),
              ", ".join(f"{s[0]}: {s[2]} pts {'=' if s[1] else '!='}" for s in same))
        old = freqs(STOCK, kind, "lin 2 1k 1meg", "stock")
        check(f"b {kind} stock lin 2 had one point (comparator sees it)", old == [1e3], f"{old}")
    a, b = freqs(NGSPICE, "noise", "lin 2 1k 1meg", "new"), freqs(STOCK, "noise", "lin 2 1k 1meg", "stock")
    check("c noise lin 2 unchanged", a == b == [1e3, 1e6], f"new {a}, stock {b}")
    print(f"{RESULTS.count(True)}/{len(RESULTS)} pass")
    return 1 if RESULTS.count(False) else 0


if __name__ == "__main__":
    sys.exit(main())
