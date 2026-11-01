#!/usr/bin/env python3
"""35_freerun_verdict.py -- handbook example 35 (card #543): the stability verdict of a free-running (trigger=) .pss, card #540.

    python3 35_freerun_verdict.py NGSPICE OUTDIR

Deck: tests/fixtures/pss_r5_shoot.cir (R5 hysteretic buck, fsw 422.7 kHz) + two synthetic extracted RC trees (bench/gen_pex_buck.py
tree_lines, branch 1, depth 2: on vout through 10 ohm, on the switch node sw through 0.5 ohm) = 10 states, `.options klu`.
Four runs, one ngspice at a time, OMP_NUM_THREADS=2:
    matrix-free (the default) / dense (NGSPICE_PSS_KRYLOV=0)  x  normal stability pass / forced noisy pass
    (PSS_FREERUN_STAB_RELTOL=1e-5 PSS_FREERUN_STAB_EPS=1e-6 PSS_FREERUN_STAB_RIPPLE=1e-5: the finite-difference noise of s90 and earlier)
The independent truth is not computed here: card #540 measured the slowest multiplier of this orbit by perturbing the settled orbit in a
plain .tran (0.037; 0.036-0.043 over cycles 3-5, two seeds; docs/pss-pac.md) -- cited, tests/test_handbook.py checks against it.
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bench"))
import gen_pex_buck as g  # noqa: E402

NOISY = {"PSS_FREERUN_STAB_RELTOL": "1e-5", "PSS_FREERUN_STAB_EPS": "1e-6", "PSS_FREERUN_STAB_RIPPLE": "1e-5"}


def build():
    s = open(os.path.join(ROOT, "tests", "fixtures", "pss_r5_shoot.cir")).read()
    s = s.replace("../../stb/lgprobe.lib", os.path.join(ROOT, "stb", "lgprobe.lib"))
    card = re.search(r"(?m)^pss .*$", s).group(0)
    b = re.sub(r"(?ms)^\.control.*?^\.endc\s*$", "", s)
    b = "\n".join(l for l in b.splitlines() if l.strip().lower() != ".end")
    lines = []
    for nm, anchor, rc in (("A", "vout", "10"), ("B", "sw", "0.5")):
        tl, root, far, nn = g.tree_lines(1, 2, 0.2, 1.6, 20e-15, 0.6, 3, 40, 2e-15)
        lines.append("Rcouple_%s %s %s %s\n" % (nm, anchor, root.replace("tree_", "tree" + nm + "_"), rc))
        for l in tl:
            if l.startswith("*"):
                continue
            t = l.split(None, 1)
            lines.append(t[0] + "_" + nm + " " + t[1].replace("tree_", "tree" + nm + "_"))
    return (b + "\n* PEX trees (handbook example 35)\n" + "".join(lines) + ".save v(vout) v(sw) i(L1)\n.options klu\n.control\nset wr_singlescale\n"
            "set wr_vecnames\nset numdgt=15\n" + card + "\nwrdata orbit.txt v(vout) v(sw) i(L1)\nquit\n.endc\n.end\n")


def main():
    ng, out = sys.argv[1], sys.argv[2]
    for tag, env in (("matrixfree_normal", {}), ("dense_normal", {"NGSPICE_PSS_KRYLOV": "0"}),
                     ("matrixfree_forced", NOISY), ("dense_forced", dict(NOISY, NGSPICE_PSS_KRYLOV="0"))):
        d = os.path.join(out, tag)
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "deck.cir"), "w").write(build())
        e = dict(os.environ, OMP_NUM_THREADS="2")
        for k in ("NGSPICE_PSS_KRYLOV", "PSS_KRYLOV_AUTO_MIN", "PSS_FREERUN_AUTO_MIN"):
            e.pop(k, None)
        e.update(env)
        p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", "deck.cir"],
                           capture_output=True, text=True, cwd=d, env=e, timeout=1800)
        o = (p.stdout + p.stderr).replace("\r", "\n")
        open(os.path.join(d, "stdout.txt"), "w").write(o)
        st = re.search(r"PSS-FREERUN: (\d+) states >= 1: matrix-free", o)
        notice = "notice" if st else "no notice"
        c = re.search(r"PSS-FREERUN: CONVERGED iter=(\d+) residual=\S+ period=(\S+) fsw=\S+ spectral_radius=(\S+) (\S+) wall_time=(\S+)s", o)
        why = re.search(r"PSS-FREERUN: stability not computed: ([^\n]*)", o)
        if not c:
            print("FRVERDICT %s FAILED rc=%d: %s" % (tag, p.returncode, o[-300:]))
            continue
        print("FRVERDICT %s: %s%s, iter %s, period %s, radius %s %s, wall %s s%s" % (
            tag, notice, (" (%s states)" % st.group(1)) if st else "", c.group(1), c.group(2), c.group(3), c.group(4), c.group(5),
            "; stability not computed: " + why.group(1)[:160] if why else ""), flush=True)


if __name__ == "__main__":
    main()
