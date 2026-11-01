#!/usr/bin/env python3
"""46_pss_memory.py -- handbook example 46 (card #560): `.pss method=shoot` no longer keeps the output plot of every period run (card #553, patch 0109).

    python3 46_pss_memory.py NGSPICE OUTDIR

Deck: the PEX-size R3 buck of example 36 (R3 voltage-mode buck + three synthetic extracted RC trees, `.options klu`, matrix-free `.pss`), with the
`.save` line of example 36 REMOVED (`.save` hides the effect: four vectors of 8 B per row), at two sizes: branch 3 depth 4 = 490 capacitor/inductor
states and branch 4 depth 4 = 1150.  Each deck is run (one ngspice at a time, OMP_NUM_THREADS=2, MemoryMax 8G) with the peak RSS read from the child's rusage (ru_maxrss):
    before   the s108 build next to the prefix (no patch 0109), when present
    after    NGSPICE (s115)
    switch   NGSPICE with PSS_KEEP_RUN_DATA=1 (the old retention, the red proof of tests/test_pss_krylov.py)
and the 490-state deck with `.save` on the s108 and the s115 build.  Lines:
    MEMRUN <deck> <arm>: peak RSS <MB> MB, wall <s> s, period runs <n>, rho <r>, orbit md5 <8 hex>
    MEMPLOTS <build>: <plots listed after .pss of example 06>   (the intermediate "PSS shoot" plots, `echo $plots`)
and a summary MEMSAME line: the wrdata orbit of every arm of a deck is byte-identical.
"""
import hashlib
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bench"))
import gen_pex_buck as g  # noqa: E402


def build(branch, depth, save):
    s = open(os.path.join(ROOT, "tests", "fixtures", "pac_r3_ctrl2out.cir")).read()
    b = s.replace("../../stb/lgprobe.lib", os.path.join(ROOT, "stb", "lgprobe.lib"))
    b = "\n".join(l for l in b.splitlines() if l.strip().lower() != ".end")
    b = re.sub(r"(?ms)^\.control.*?^\.endc\s*$", "", b)
    b = b.replace("Cout cx 0 {cval} IC=0", "Cout cx gret {cval} IC=0\nRgret gret 0 5m")
    lines = []
    for nm, anchor, rc in (("A", "vout", "10"), ("B", "sw", "0.5"), ("G", "gret", "1")):
        tl, root, far, nn = g.tree_lines(branch, depth, 0.2, 1.6, 20e-15, 0.6, 3, 40, 2e-15)
        lines.append("Rcouple_%s %s %s %s\n" % (nm, anchor, root.replace("tree_", "tree" + nm + "_"), rc))
        for l in tl:
            if l.startswith("*"):
                continue
            t = l.split(None, 1)
            lines.append(t[0] + "_" + nm + " " + t[1].replace("tree_", "tree" + nm + "_"))
    pex = "* PEX-size parasitic nets (handbook example 46)\n" + "".join(lines) + (".save v(vout) v(sw) v(opout) i(L1)\n" if save else "")
    return (b + "\n" + pex + "\n.options klu\n.control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
            "pss 500e3 300e-6 vout 200 10 150 1e-6 uic method=shoot\nwrdata orbit.txt v(vout) v(sw) v(opout) i(L1)\nquit\n.endc\n.end\n")


def run(ng, out, name, deck, env=None):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "deck.cir"), "w").write(deck)
    e = dict(os.environ, OMP_NUM_THREADS="2", PSS_SHOOT_COUNT_RUNS="1")
    e.pop("PSS_KEEP_RUN_DATA", None)
    e.update(env or {})
    t0 = time.time()
    # peak RSS = ru_maxrss of the waited-for child (systemd-run --scope execs the program in the same process; /usr/bin/time is not on every machine)
    p = subprocess.Popen(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", "deck.cir"],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=d, env=e)
    outb = p.stdout.read()
    _, status, ru = os.wait4(p.pid, 0)
    p.returncode = 0
    wall = time.time() - t0
    o = outb.decode(errors="replace").replace("\r", "\n")
    open(os.path.join(d, "stdout.txt"), "w").write(o)
    m = re.search(r"PEAKRSS_KB (\d+)", "PEAKRSS_KB %d" % ru.ru_maxrss)
    rn = re.search(r"PSS-SHOOT: period runs=(\d+)", o)
    rho = re.search(r"CONVERGED iter=\d+ residual=\S+ spectral_radius=(\S+) (\S+)", o)
    md5 = hashlib.md5(open(os.path.join(d, "orbit.txt"), "rb").read()).hexdigest()[:8] if os.path.exists(os.path.join(d, "orbit.txt")) else "none"
    return dict(rss=int(m.group(1)) / 1024.0 if m else -1, wall=wall, runs=int(rn.group(1)) if rn else -1, rho=(rho.group(1) + " " + rho.group(2)) if rho else "none", md5=md5)


def main(ng, out):
    old = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(ng))), "aaba50ba3f49-s108", "bin", "ngspice")
    have_old = os.path.exists(old)
    cases = [("490 states", 3, 4, False), ("1150 states", 4, 4, False), ("490 states, .save", 3, 4, True)]
    for label, br, dp, save in cases:
        arms = ([("before (s108)", old, None)] if have_old else []) + [("after (s115)", ng, None)]
        if not save:
            arms.append(("after, PSS_KEEP_RUN_DATA=1", ng, {"PSS_KEEP_RUN_DATA": "1"}))
        res = []
        for an, exe, env in arms:
            r = run(exe, out, ("%d_%d_%d_%s" % (br, dp, save, an.split()[0])), build(br, dp, save), env)
            res.append(r)
            print("MEMRUN %s, %s: peak RSS %.0f MB, wall %.1f s, period runs %d, rho %s, orbit md5 %s" % (label, an, r["rss"], r["wall"], r["runs"], r["rho"], r["md5"]), flush=True)
        print("MEMSAME %s: orbit files byte-identical across %d arms: %s" % (label, len(res), len({r["md5"] for r in res}) == 1), flush=True)
    deck = open(os.path.join(HERE, "06_pss_shoot.cir")).read()
    deck = deck.replace("print vx", "").replace(".endc", "echo plots: $plots\n.endc")
    for tag, exe in ([("s108", old)] if have_old else []) + [("s115", ng)]:
        d = os.path.join(out, "plots_" + tag)
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "deck.cir"), "w").write(deck)
        p = subprocess.run([exe, "-b", "deck.cir"], capture_output=True, text=True, cwd=d, timeout=300)
        m = re.search(r"^plots: (.*)$", p.stdout.replace("\r", "\n"), re.M)
        print("MEMPLOTS %s: %s" % (tag, m.group(1) if m else "no plots line"), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
