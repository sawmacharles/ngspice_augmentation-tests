#!/usr/bin/env python3
"""36_pex_krylov.py -- handbook example 36 (card #543): `.pss` on a converter with extracted parasitics, dense monodromy vs
matrix-free Newton-GMRES, SPARSE vs KLU.  The deck is built here, not stored:

    refs-style R3 voltage-mode buck (tests/fixtures/pac_r3_ctrl2out.cir body, its .control cut off, the output capacitor returned
    through a 5 mOhm "ground return" node) + three synthetic extracted RC trees (bench/gen_pex_buck.py tree_lines, branch 3,
    depth 4: on vout through 10 ohm, on the switch node sw through 0.5 ohm, on the capacitor's ground return through 1 ohm)
    = 490 capacitor/inductor states.   .pss 500e3 300e-6 vout 200 10 150 1e-6 uic method=shoot

    python3 36_pex_krylov.py NGSPICE OUTDIR [REPS]

Arms (one ngspice at a time, OMP_NUM_THREADS=2, MemoryMax 8G, order rotated each repetition, `uptime` load recorded):
    sparse_dense  `.options sparse` (since s112 a shooting `.pss` of >= 128 nodes takes KLU by itself unless the deck names a solver), `krylov=0` on the card
    klu_dense     `.options klu`, `krylov=0`
    klu_auto      `.options klu`, no krylov key at all  (the default: matrix-free from 64 states)
Prints one `PEXRUN` line per run and, per arm, the median/min/max wall; then the orbit deviation of each arm against sparse_dense
(max |dev| / max(peak-to-peak, 1e-2) over v(vout) v(sw) v(opout) i(L1), the engine's own ripple-relative norm), the Floquet
radius, the period-run counts and the `PSS-SHOOT: 490 states >= 64` notice (it must appear exactly once, in klu_auto only).
"""
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bench"))
import gen_pex_buck as g  # noqa: E402

BRANCH, DEPTH = 3, 4


def build(klu, krylov_card):
    s = open(os.path.join(ROOT, "tests", "fixtures", "pac_r3_ctrl2out.cir")).read()
    b = s.replace("../../stb/lgprobe.lib", os.path.join(ROOT, "stb", "lgprobe.lib"))
    b = "\n".join(l for l in b.splitlines() if l.strip().lower() != ".end")
    b = re.sub(r"(?ms)^\.control.*?^\.endc\s*$", "", b)
    b = b.replace("Cout cx 0 {cval} IC=0", "Cout cx gret {cval} IC=0\nRgret gret 0 5m")
    assert "Rgret" in b
    lines = []
    for nm, anchor, rc in (("A", "vout", "10"), ("B", "sw", "0.5"), ("G", "gret", "1")):
        tl, root, far, nn = g.tree_lines(BRANCH, DEPTH, 0.2, 1.6, 20e-15, 0.6, 3, 40, 2e-15)
        lines.append("Rcouple_%s %s %s %s\n" % (nm, anchor, root.replace("tree_", "tree" + nm + "_"), rc))
        for l in tl:
            if l.startswith("*"):
                continue
            t = l.split(None, 1)
            lines.append(t[0] + "_" + nm + " " + t[1].replace("tree_", "tree" + nm + "_"))
    pex = "* PEX-size parasitic nets (handbook example 36)\n" + "".join(lines) + ".save v(vout) v(sw) v(opout) i(L1)\n"
    card = "pss 500e3 300e-6 vout 200 10 150 1e-6 uic method=shoot" + (" krylov=0" if krylov_card is False else "")
    return (b + "\n" + pex + ("\n.options klu\n" if klu else "\n.options sparse\n") + ".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n" +
            card + "\nwrdata orbit.txt v(vout) v(sw) v(opout) i(L1)\nquit\n.endc\n.end\n")


ARMS = {"sparse_dense": (False, False), "klu_dense": (True, False), "klu_auto": (True, None)}


def orbit(path):
    rows = []
    for l in open(path):
        try:
            rows.append([float(x) for x in l.split()])
        except ValueError:
            pass
    return rows


def main():
    ng, out = sys.argv[1], sys.argv[2]
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    res = {a: [] for a in ARMS}
    names = list(ARMS)
    for r in range(reps):
        order = names[r % len(names):] + names[:r % len(names)]
        for a in order:
            d = os.path.join(out, "%s_rep%d" % (a, r))
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, "deck.cir"), "w").write(build(*ARMS[a]))
            load = open("/proc/loadavg").read().split()[0]
            env = dict(os.environ, OMP_NUM_THREADS="2", PSS_SHOOT_COUNT_RUNS="1")
            env.pop("NGSPICE_PSS_KRYLOV", None)
            env.pop("PSS_KRYLOV_AUTO_MIN", None)
            t0 = time.time()
            p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", ng, "-b", "deck.cir"],
                               capture_output=True, text=True, cwd=d, env=env, timeout=1800)
            wall = time.time() - t0
            o = (p.stdout + p.stderr).replace("\r", "\n")
            open(os.path.join(d, "stdout.txt"), "w").write(o)
            m = re.search(r"PSS-SHOOT: CONVERGED iter=(\d+) residual=(\S+) spectral_radius=(\S+) (\S+) wall_time=(\S+)s", o)
            st = re.search(r"(\d+) state\(s\)", o)
            rn = re.search(r"PSS-SHOOT: period runs=(\d+) .*mode=(\S+)", o)
            notice = len(re.findall(r"states >= 64: matrix-free Newton-GMRES", o))
            if not m or not st:
                print("PEXRUN %s rep%d FAILED rc=%d: %s" % (a, r, p.returncode, o[-300:]))
                continue
            res[a].append(dict(wall=wall, load=load, iters=int(m.group(1)), rho=float(m.group(3)), verdict=m.group(4), pss=float(m.group(5)),
                               states=int(st.group(1)), runs=int(rn.group(1)) if rn else -1, mode=rn.group(2) if rn else "?", notice=notice,
                               orbit=os.path.join(d, "orbit.txt")))
            x = res[a][-1]
            print("PEXRUN %s rep%d: states %d, wall %.1f s (PSS-SHOOT %.1f s), load %s, iter %d, period runs %d (%s), rho %.6f %s, notice lines %d"
                  % (a, r, x["states"], wall, x["pss"], load, x["iters"], x["runs"], x["mode"], x["rho"], x["verdict"], x["notice"]), flush=True)
    for a in names:
        w = sorted(x["wall"] for x in res[a])
        if w:
            print("PEXWALL %s n=%d median %.1f s [min %.1f .. max %.1f]" % (a, len(w), w[len(w) // 2], w[0], w[-1]))
    if res["sparse_dense"] and all(res[a] for a in names):
        import bisect
        ref = orbit(res["sparse_dense"][0]["orbit"])
        T = 2e-6

        def at(rows, t, c):
            ts = [r[0] for r in rows]
            i = bisect.bisect_left(ts, t)
            if i <= 0:
                return rows[0][c]
            if i >= len(rows):
                return rows[-1][c]
            f = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
            return rows[i - 1][c] + f * (rows[i][c] - rows[i - 1][c])
        for a in ("klu_dense", "klu_auto"):
            o = orbit(res[a][0]["orbit"])
            dev = {}
            for c, nm in ((1, "vout"), (2, "sw"), (3, "opout"), (4, "il")):
                pp = max(max(r[c] for r in ref) - min(r[c] for r in ref), 1e-2)
                dev[nm] = max(abs(at(ref, (j + 0.5) * T / 2000, c) - at(o, (j + 0.5) * T / 2000, c)) for j in range(2000)) / pp
            print("PEXORBIT %s vs sparse_dense (orbits on a 2000-point grid, max |dev|/max(p-p,1e-2)): vout %.2e, opout %.2e, il %.2e, sw %.2e (edge vector)"
                  % (a, dev["vout"], dev["opout"], dev["il"], dev["sw"]))


if __name__ == "__main__":
    main()
