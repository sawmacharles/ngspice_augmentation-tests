#!/usr/bin/env python3
"""47_speed_r16_r17.py -- handbook example 47 (card #560): the speed rounds 16 and 17 (cards #555, #558) on THIS machine, same binary, off-switches,
interleaved (start position rotated per repetition), one ngspice at a time (PSS_FORK_COLS arms: the parent plus at most 2 forked workers), `set num_threads=2`
in every deck, `uptime` load recorded.  Lines printed (the handbook quotes them):

    AKLU   automatic KLU for `.pss method=shoot` (>= 128 nodes): three PEX bucks (R3 + 3 extracted RC trees, no `.options klu` in the deck), default vs
           NGSPICE_PSS_AUTOKLU=0 vs `.options sparse` written in the deck; the notice line, the median wall, the orbit deviation vs the SPARSE run
    FORK   PSS_FORK_COLS=2 on dense `.pss` (krylov=0): identity of the wrdata file against the sequential run, median wall
    MRHS   PAC_MULTIRHS=1 on a td `.pstb` under `.options klu`: identity and median wall
    SEED   PSS_SEED_STAGE1=1 on dense `.pss`: iterations, residual, radius, orbit deviation vs default, median wall
    ARN    PSS_KRYLOV_ARN_EARLY=0 (Arnoldi stop) on the 490-state PEX deck: period runs, radius;  PSS_FLOQUET_BASELINE=1 on the dense R3d deck: period runs

Usage: 47_speed_r16_r17.py NGSPICE OUTDIR [REPS]
"""
import bisect
import hashlib
import os
import re
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bench"))
import gen_pex_buck as g  # noqa: E402

NG, OUT = sys.argv[1], sys.argv[2]
REPS = int(sys.argv[3]) if len(sys.argv) > 3 else 3
os.makedirs(OUT, exist_ok=True)
ENVOFF = ("NGSPICE_PSS_AUTOKLU", "NGSPICE_PSS_AUTOKLU_MIN", "PSS_FORK_COLS", "PAC_MULTIRHS", "PSS_SEED_STAGE1", "PSS_KRYLOV_ARN_EARLY", "PSS_FLOQUET_BASELINE",
          "NGSPICE_PSS_KRYLOV", "PSS_KRYLOV_AUTO_MIN", "PSS_FORK_TEST_FAIL")


def pex_deck(branch, depth, solver="", krylov0=False, wr="orbit.txt"):
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
    pex = "* PEX-size parasitic nets (handbook example 47)\n" + "".join(lines) + ".save v(vout) v(sw) v(opout) i(L1)\n"
    return (b + "\n" + pex + solver + ".control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
            "pss 500e3 300e-6 vout 200 10 150 1e-6 uic method=shoot" + (" krylov=0" if krylov0 else "") + "\nwrdata " + wr + " v(vout) v(sw) v(opout) i(L1)\nquit\n.endc\n.end\n")


def file_deck(name):
    t = open(os.path.join(HERE, name)).read()
    return t.replace("set filetype=ascii\n", "set filetype=ascii\nset num_threads=2\n", 1)


def run(deck, env, tag):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "deck.cir"), "w").write(deck)
    e = dict(os.environ, OMP_WAIT_POLICY="passive", PSS_SHOOT_COUNT_RUNS="1")
    for k in ENVOFF:
        e.pop(k, None)
    e.update(env)
    load = open("/proc/loadavg").read().split()[0]
    t0 = time.time()
    p = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", NG, "-b", "deck.cir"],
                       cwd=d, env=e, capture_output=True, text=True, timeout=1800)
    wall = time.time() - t0
    o = (p.stdout + p.stderr).replace("\r", "\n")
    open(os.path.join(d, "stdout.txt"), "w").write(o)
    m = re.search(r"PSS-SHOOT: CONVERGED iter=(\d+) residual=(\S+) spectral_radius=(\S+) (\S+)", o)
    rn = re.search(r"PSS-SHOOT: period runs=(\d+)", o)
    outs = [f for f in os.listdir(d) if f.endswith(".txt") and f != "stdout.txt"]
    md5 = hashlib.md5(b"".join(open(os.path.join(d, f), "rb").read() for f in sorted(outs))).hexdigest()[:8]
    return dict(wall=wall, load=load, out=o, dir=d, md5=md5, rc=p.returncode,
                iters=int(m.group(1)) if m else -1, resid=m.group(2) if m else "-", rho=m.group(3) if m else "-", verdict=m.group(4) if m else "-",
                runs=int(rn.group(1)) if rn else -1, rcode=p.returncode)


def interleaved(arms, reps, pre="x"):
    """arms: {name: (deck, env)}; returns {name: [result, ...]}"""
    names = list(arms)
    res = {n: [] for n in names}
    for r in range(reps):
        for n in names[r % len(names):] + names[:r % len(names)]:
            res[n].append(run(arms[n][0], arms[n][1], "%s_%s_r%d" % (pre, n.replace(" ", "_").replace("=", ""), r)))
    return res


def med(rs):
    w = sorted(x["wall"] for x in rs)
    return statistics.median(w), w[0], w[-1]


def orbit(path):
    rows = []
    for l in open(path):
        try:
            rows.append([float(x) for x in l.split()])
        except ValueError:
            pass
    return rows


def at(rows, t, c):
    ts = [r[0] for r in rows]
    i = bisect.bisect_left(ts, t)
    if i <= 0:
        return rows[0][c]
    if i >= len(rows):
        return rows[-1][c]
    f = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
    return rows[i - 1][c] + f * (rows[i][c] - rows[i - 1][c])


def dev(pa, pb, cols, T=2e-6):
    a, b = orbit(pa), orbit(pb)
    w = 0.0
    for c in cols:
        pp = max(max(r[c] for r in a) - min(r[c] for r in a), 1e-2)
        w = max(w, max(abs(at(a, (j + 0.5) * T / 2000, c) - at(b, (j + 0.5) * T / 2000, c)) for j in range(2000)) / pp)
    return w


SECT = os.environ.get("SPEED_SECTIONS", "ABCDE")


def main():
    # ---------------------------------------------------------------- A: automatic KLU
    for br, dp in (((2, 4), (3, 3), (3, 4)) if "A" in SECT else ()):
        arms = {"default": (pex_deck(br, dp), {}), "AUTOKLU=0": (pex_deck(br, dp), {"NGSPICE_PSS_AUTOKLU": "0"}),
                "options sparse": (pex_deck(br, dp, solver="\n.options sparse\n"), {})}
        res = interleaved(arms, REPS, "aklu%d%d" % (br, dp))
        d0 = res["default"][0]
        note = re.search(r"Note: periodic analysis switching to KLU \((\d+) nodes >= 128", d0["out"])
        n_notes = {k: sum(1 for x in v if "switching to KLU" in x["out"]) for k, v in res.items()}
        dv = dev(res["AUTOKLU=0"][0]["dir"] + "/orbit.txt", d0["dir"] + "/orbit.txt", (1, 3, 4))
        ident = d0["md5"] == res["AUTOKLU=0"][0]["md5"]
        st = re.search(r"(\d+) state\(s\)", d0["out"])
        w = {k: med(v) for k, v in res.items()}
        print("AKLU b%dd%d (%s states, %s): notice in default %s, AUTOKLU=0 %s, options sparse %s; median wall default %.2f s, AUTOKLU=0 %.2f s (%+.1f %%), options sparse %.2f s; "
              "output byte-identical to AUTOKLU=0: %s; orbit deviation default vs SPARSE %.1e of ripple; load %s"
              % (br, dp, st.group(1) if st else "?", ("%s nodes" % note.group(1)) if note else "below 128 nodes", n_notes["default"], n_notes["AUTOKLU=0"], n_notes["options sparse"],
                 w["default"][0], w["AUTOKLU=0"][0], 100 * (w["default"][0] / w["AUTOKLU=0"][0] - 1), w["options sparse"][0], ident, dv, d0["load"]), flush=True)

    # ---------------------------------------------------------------- B: fork-parallel monodromy columns
    cases = [("R3d (example 26)", file_deck("26_pss_r3d_deadtime.cir"), 2), ("R4d (example 27)", file_deck("27_pss_r4d_deadtime.cir"), 2),
             ("PEX 73 states, krylov=0", pex_deck(2, 3, krylov0=True), 2)]
    for name, deck, k in (cases if "B" in SECT else ()):
        res = interleaved({"sequential": (deck, {}), "FORK_COLS=%d" % k: (deck, {"PSS_FORK_COLS": str(k)})}, REPS, "fork" + name.split()[0])
        a, b = res["sequential"], res["FORK_COLS=%d" % k]
        ident = len({x["md5"] for x in a + b}) == 1
        wa, wb = med(a), med(b)
        print("FORK %s: sequential median %.2f s [%.2f .. %.2f], PSS_FORK_COLS=%d median %.2f s [%.2f .. %.2f] (%+.1f %% median, %+.1f %% min); wrdata byte-identical: %s; iterations %d / %d, residual %s / %s; load %s"
              % (name, wa[0], wa[1], wa[2], k, wb[0], wb[1], wb[2], 100 * (wb[0] / wa[0] - 1), 100 * (wb[1] / wa[1] - 1), ident, a[0]["iters"], b[0]["iters"], a[0]["resid"], b[0]["resid"], a[0]["load"]), flush=True)

    # ---------------------------------------------------------------- C: multi-RHS
    if "C" in SECT:
        deck = file_deck("28_pstb_r3d_deadtime.cir").replace(".control", ".options klu\n.control", 1)
        res = interleaved({"column by column": (deck, {}), "PAC_MULTIRHS=1": (deck, {"PAC_MULTIRHS": "1"})}, max(REPS, 5), "mrhs")
        a, b = res["column by column"], res["PAC_MULTIRHS=1"]
        wa, wb = med(a), med(b)
        print("MRHS R3d .pstb (.options klu), n=%d: default median %.3f s [%.3f .. %.3f], PAC_MULTIRHS=1 %.3f s [%.3f .. %.3f] (%+.1f %% median, %+.1f %% min); output byte-identical: %s; load %s"
              % (len(a), wa[0], wa[1], wa[2], wb[0], wb[1], wb[2], 100 * (wb[0] / wa[0] - 1), 100 * (wb[1] / wa[1] - 1), len({x["md5"] for x in a + b}) == 1, a[0]["load"]), flush=True)

    # ---------------------------------------------------------------- D: stage-1 seed
    r7dir = os.path.join(OUT, "r7decks")
    os.makedirs(r7dir, exist_ok=True)
    subprocess.run([sys.executable, os.path.join(HERE, "32_r7_decks.py"), r7dir], check=True, capture_output=True)
    r7 = open(os.path.join(r7dir, "r7_pss.cir")).read().replace(".control\n", ".control\nset num_threads=2\n", 1)
    seed_cases = [("R3d (example 26)", file_deck("26_pss_r3d_deadtime.cir"), (1, 2, 3), REPS), ("R4d (example 27)", file_deck("27_pss_r4d_deadtime.cir"), (1, 2, 3), REPS),
                  ("R7 BSIM4 (example 32)", r7, (1, 2), 3)]
    for name, deck, cols, reps in (seed_cases if "D" in SECT else []):
        res = interleaved({"default": (deck, {}), "SEED_STAGE1=1": (deck, {"PSS_SEED_STAGE1": "1"})}, reps, "seed" + name.split()[0])
        a, b = res["default"], res["SEED_STAGE1=1"]
        wa, wb = med(a), med(b)
        fn = [f for f in os.listdir(a[0]["dir"]) if f.endswith(".txt") and f != "stdout.txt"][0]
        print("SEED %s: default iter %d residual %s radius %s %s, median %.2f s [%.2f .. %.2f]; PSS_SEED_STAGE1=1 iter %d residual %s radius %s %s, median %.2f s [%.2f .. %.2f] (%+.1f %% median, %+.1f %% min); "
              "orbit deviation %.1e of ripple; load %s"
              % (name, a[0]["iters"], a[0]["resid"], a[0]["rho"], a[0]["verdict"], wa[0], wa[1], wa[2], b[0]["iters"], b[0]["resid"], b[0]["rho"], b[0]["verdict"], wb[0], wb[1], wb[2],
                 100 * (wb[0] / wa[0] - 1), 100 * (wb[1] / wa[1] - 1), dev(a[0]["dir"] + "/" + fn, b[0]["dir"] + "/" + fn, [c for c in cols if c <= 3]), a[0]["load"]), flush=True)

    # ---------------------------------------------------------------- E: Arnoldi stop, Floquet baseline
    deck = pex_deck(3, 4)
    res = interleaved({"default": (deck, {}), "ARN_EARLY=0": (deck, {"PSS_KRYLOV_ARN_EARLY": "0"})}, REPS, "arn")
    a, b = res["default"], res["ARN_EARLY=0"]
    wa, wb = med(a), med(b)
    print("ARN PEX 490 states matrix-free: default %d period runs, radius %s, median %.2f s; PSS_KRYLOV_ARN_EARLY=0 %d period runs, radius %s, median %.2f s (%+.1f %% for the default); output byte-identical: %s"
          % (a[0]["runs"], a[0]["rho"], wa[0], b[0]["runs"], b[0]["rho"], wb[0], 100 * (wa[0] / wb[0] - 1), a[0]["md5"] == b[0]["md5"]), flush=True)
    deck = file_deck("26_pss_r3d_deadtime.cir")
    res = interleaved({"default": (deck, {}), "FLOQUET_BASELINE=1": (deck, {"PSS_FLOQUET_BASELINE": "1"})}, REPS, "flq")
    a, b = res["default"], res["FLOQUET_BASELINE=1"]
    print("ARN R3d dense: default %d period runs, PSS_FLOQUET_BASELINE=1 %d period runs; output byte-identical: %s" % (a[0]["runs"], b[0]["runs"], a[0]["md5"] == b[0]["md5"]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
