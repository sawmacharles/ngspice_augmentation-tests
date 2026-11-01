#!/usr/bin/env python3
"""card #558, acceleration round 17: automatic KLU for shooting periodic analyses, fork()-parallel dense monodromy columns,
multi-right-hand-side KLU solves in the time-domain PAC recursion (opt-in, measured: no gain).

  autoklu_*  the size rule: a `.pss method=shoot` deck with >= 128 nodes switches to KLU by itself (one note line),
             a smaller deck, an explicit `.options klu` / `.options sparse` and NGSPICE_PSS_AUTOKLU=0 do not;
             NGSPICE_PSS_AUTOKLU_MIN moves the threshold; the orbit under the automatic KLU equals the SPARSE orbit to the
             KLU-vs-SPARSE tolerance of tests/test_pss_klu.py.
  fork_*     PSS_FORK_COLS=<n> (opt-in): R3, R3d and the dense PEX deck give a byte-identical orbit file against the sequential
             loop; for every deck the result is byte-identical for every worker count (a column starts from the process as
             it was on entry); R4d (algebraic latch: the sequential chain carries the previous run's node vector) agrees
             with the sequential orbit to the orbit tolerance only; a dying worker (PSS_FORK_TEST_FAIL=1) falls back to the
             sequential loop with the identical result and one message.
  multirhs_* the td PAC / .pstb recursion with the columns solved in one klu_solve (default) vs column by column
             (opt-in PAC_MULTIRHS=1; default column by column): byte-identical curves on a native-capacitor deck (SPARSE and KLU) and on the BSIM4 R7 deck.
  red proofs: the identity checks are run against a deliberately perturbed deck (reltol 1e-3 for the second run) and must FAIL.

Usage:  NGSPICE=<binary> [NGSPICE_ASAN=<asan binary>] python3 tests/test_accel17.py [--quick] [--only substr]
Exit 0 pass, 1 fail, 2 could-not-run.   Outputs: ~/Work/data/ngspice_augmentation/558/test (never the repo).
Every run sets `set num_threads=2` in the deck (the spinit default of 8 would ignore OMP_NUM_THREADS)."""
import os, re, sys, hashlib, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bench"))
import test_pss_klu as K   # body(), redirect(), orbit_dev(), pex_deck(), strip_end()

NGSPICE = os.environ.get("NGSPICE", "ngspice")
ASAN = os.environ.get("NGSPICE_ASAN", "")
OUT = os.path.expanduser("~/Work/data/ngspice_augmentation/558/test")
QUICK = "--quick" in sys.argv
ONLY = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else ""
TOL_ORBIT, TOL_EDGE = 1e-4, 2e-2
rows = []


def record(name, status, detail=""):
    if ONLY and ONLY not in name:
        return
    rows.append((name, status))
    print(f"{status:<13} {name}: {detail}", flush=True)


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def prep(deck, d, nt=2):
    """redirect wrdata/write into d, drop raw writes, cap the OpenMP threads"""
    deck = K.redirect(deck, d)
    deck = re.sub(r"(?mi)^\s*write\s+.*\n", "", deck)
    deck = re.sub(r"(?m)^\.control\s*$", ".control\nset num_threads=%d\nset numdgt=15" % nt, deck, count=1)
    if "\nquit" not in deck:
        deck = deck.replace("\n.endc", "\nquit\n.endc", 1)
    return deck


def run(tag, deck, env=None, binary=None, timeout=1800):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):
        if f.endswith(".txt"):
            os.unlink(os.path.join(d, f))
    p = os.path.join(d, "deck.cir")
    open(p, "w").write(prep(deck, d))
    e = dict(os.environ, OMP_NUM_THREADS="2")
    for k in ("PSS_FORK_COLS", "PSS_FORK_TEST_FAIL", "NGSPICE_PSS_AUTOKLU", "NGSPICE_PSS_AUTOKLU_MIN", "PAC_MULTIRHS"):
        e.pop(k, None)
    if env:
        e.update(env)
    t0 = time.time()
    cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=8G", "-p", "MemorySwapMax=0", binary or NGSPICE, "-b", p]
    pr = subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=e, timeout=timeout)
    o = pr.stdout + pr.stderr
    open(os.path.join(d, "log.txt"), "w").write(o)
    outs = sorted(f for f in os.listdir(d) if f.endswith(".txt") and f != "log.txt")
    return dict(rc=pr.returncode, log=o, files=[os.path.join(d, f) for f in outs], wall=time.time() - t0, d=d)


def same_files(a, b):
    return len(a["files"]) > 0 and [os.path.basename(x) for x in a["files"]] == [os.path.basename(x) for x in b["files"]] and \
        all(md5(x) == md5(y) for x, y in zip(a["files"], b["files"]))


def orbit_ok(a, b):
    worst = 0.0
    for fa, fb in zip(a["files"], b["files"]):
        for dev, edge in K.orbit_dev(fa, fb):
            worst = max(worst, dev / (TOL_EDGE if edge else TOL_ORBIT))
    return worst <= 1.0, worst


def fix(name):
    return K.body(name)


def sparse_deck(deck):
    return re.sub(r"(?mi)^\.options\s+klu\s*$", "", deck)


def klu_deck(deck):
    return re.sub(r"^\.control", ".options klu\n.control", sparse_deck(deck), count=1, flags=re.M)


def with_opt(deck, opt):
    return re.sub(r"^\.control", ".options %s\n.control" % opt, sparse_deck(deck), count=1, flags=re.M)


NOTE = "periodic analysis switching to KLU"


def pex_deck(branch, depth):
    """R3 + three RC trees (as tests/test_pss_krylov.py pex_deck), no solver option; the matrix-free engine takes over at >= 64 states"""
    import test_pss_krylov as T
    deck, nelem = T.pex_deck(branch, depth, "@D@", None)
    return sparse_deck(deck)


# ------------------------------------------------------------------ checks
def check_autoklu():
    big = pex_deck(3, 3)            # 136 nodes
    small = fix("pss_r3_shoot.cir")  # 15 nodes
    a = run("autoklu/big_auto", big)
    note = a["log"].count(NOTE)
    record("autoklu: 136-node .pss method=shoot switches to KLU, one note", "pass" if a["rc"] == 0 and note == 1 else "FAIL",
           f"rc {a['rc']}, notes {note}, wall {a['wall']:.1f}s")
    off = run("autoklu/big_off", big, {"NGSPICE_PSS_AUTOKLU": "0"})
    record("autoklu: NGSPICE_PSS_AUTOKLU=0 keeps SPARSE", "pass" if off["rc"] == 0 and NOTE not in off["log"] else "FAIL",
           f"notes {off['log'].count(NOTE)}")
    sp = run("autoklu/big_sparse_opt", with_opt(big, "sparse"))
    record("autoklu: explicit .options sparse wins", "pass" if sp["rc"] == 0 and NOTE not in sp["log"] else "FAIL", f"notes {sp['log'].count(NOTE)}")
    kl = run("autoklu/big_klu_opt", with_opt(big, "klu"))
    record("autoklu: explicit .options klu: no automatic note", "pass" if kl["rc"] == 0 and NOTE not in kl["log"] else "FAIL", f"notes {kl['log'].count(NOTE)}")
    s = run("autoklu/small", small)
    record("autoklu: 15-node deck stays SPARSE", "pass" if s["rc"] == 0 and NOTE not in s["log"] else "FAIL", f"notes {s['log'].count(NOTE)}")
    t = run("autoklu/small_min", small, {"NGSPICE_PSS_AUTOKLU_MIN": "10"})
    record("autoklu: NGSPICE_PSS_AUTOKLU_MIN=10 moves the threshold", "pass" if t["rc"] == 0 and t["log"].count(NOTE) == 1 else "FAIL", f"notes {t['log'].count(NOTE)}")
    tr = run("autoklu/tran_only", fix("pss_r3_shoot.cir").replace("pss 500e3 300e-6 vout 200 10 150 1e-6 uic method=shoot", "tran 1u 20u"),
             {"NGSPICE_PSS_AUTOKLU_MIN": "1"})
    record("autoklu: a .tran job never takes the periodic rule", "pass" if NOTE not in tr["log"] else "FAIL", f"notes {tr['log'].count(NOTE)}")
    # orbit equality, automatic KLU vs SPARSE (same binary, same deck)
    ok, w = orbit_ok(a, sp) if a["files"] and sp["files"] else (False, -1)
    record("autoklu: orbit under the automatic KLU equals the SPARSE orbit (tolerance of test_pss_klu)", "pass" if ok else "FAIL", f"worst deviation / tolerance {w:.3g}")
    # red proof: the same comparison against a deck with a changed period must fail
    red = run("autoklu/big_red", big.replace("pss 500e3", "pss 505e3"))
    ok2, w2 = orbit_ok(a, red) if red["files"] else (False, -1)
    record("autoklu: red proof (period moved by 1 %) must fail the orbit comparison", "pass" if not ok2 else "FAIL", f"worst/tol {w2:.3g}")


def check_fork():
    for name, fx in (("r3", "pss_r3_shoot.cir"), ("r3d", "rd_pss_r3d_shoot.cir")):
        deck = fix(fx)
        base = run(f"fork/{name}_base", deck)
        f1 = run(f"fork/{name}_k1", deck, {"PSS_FORK_COLS": "1"})
        f3 = run(f"fork/{name}_k3", deck, {"PSS_FORK_COLS": "3"})
        tr = f3["log"].count("worker failed")
        record(f"fork: {name} dense .pss, PSS_FORK_COLS=1 and 3 byte-identical to the sequential loop",
               "pass" if base["rc"] == f1["rc"] == f3["rc"] == 0 and same_files(base, f1) and same_files(base, f3) and tr == 0 else "FAIL",
               f"rc {base['rc']}/{f1['rc']}/{f3['rc']}, wall {base['wall']:.1f}/{f1['wall']:.1f}/{f3['wall']:.1f}s")
    deck = fix("rd_pss_r4d_shoot.cir")
    base = run("fork/r4d_base", deck)
    f1 = run("fork/r4d_k1", deck, {"PSS_FORK_COLS": "1"})
    f3 = run("fork/r4d_k3", deck, {"PSS_FORK_COLS": "3"})
    f5 = run("fork/r4d_k5", deck, {"PSS_FORK_COLS": "5"})
    ok, w = orbit_ok(base, f3) if f3["files"] and base["files"] else (False, -1)
    record("fork: r4d (algebraic latch) identical for every worker count", "pass" if same_files(f1, f3) and same_files(f3, f5) else "FAIL",
           f"K=1/3/5 md5 equal {same_files(f1, f3)}/{same_files(f3, f5)}")
    record("fork: r4d agrees with the sequential orbit to the orbit tolerance", "pass" if ok else "FAIL", f"worst deviation / tolerance {w:.3g}; byte-identical {same_files(base, f3)}")
    # dense PEX deck (73 states, krylov=0): the large-n case
    px = pex_deck(2, 3)
    pb = run("fork/pex_base", px, {"NGSPICE_PSS_KRYLOV": "0"})
    p4 = run("fork/pex_k4", px, {"NGSPICE_PSS_KRYLOV": "0", "PSS_FORK_COLS": "4"})
    record("fork: dense 73-state PEX .pss byte-identical to the sequential loop", "pass" if pb["rc"] == p4["rc"] == 0 and same_files(pb, p4) else "FAIL",
           f"wall {pb['wall']:.1f}s sequential, {p4['wall']:.1f}s with 4 workers")
    # a dying worker
    fl = run("fork/r3d_fail", fix("rd_pss_r3d_shoot.cir"), {"PSS_FORK_COLS": "3", "PSS_FORK_TEST_FAIL": "1"})
    b3 = run("fork/r3d_base2", fix("rd_pss_r3d_shoot.cir"))
    nm = fl["log"].count("a forked column worker failed")
    record("fork: dying workers fall back to the sequential loop (identical result, message per phase)",
           "pass" if fl["rc"] == 0 and nm >= 3 and same_files(fl, b3) else "FAIL", f"rc {fl['rc']}, messages {nm}")
    # red proof: identity against a deck at another tolerance must fail
    rd = run("fork/r3d_red", fix("rd_pss_r3d_shoot.cir").replace("pss 500e3", "pss 501e3"), {"PSS_FORK_COLS": "3"})
    record("fork: red proof (period moved by 0.2 %) must break the byte identity", "pass" if rd["files"] and not same_files(rd, b3) else "FAIL", "")
    if ASAN and not QUICK:
        a = run("fork/asan_r4d", deck, {"PSS_FORK_COLS": "3", "ASAN_OPTIONS": "detect_leaks=0"}, binary=ASAN)
        record("fork: ASan run of r4d with 3 workers", "pass" if a["rc"] == 0 and "AddressSanitizer" not in a["log"] else "FAIL", f"rc {a['rc']}")


def check_multirhs():
    for name, deck in (("r3d_pstb sparse", sparse_deck(fix("rd_pstb_r3d.cir"))), ("r3d_pstb klu", klu_deck(fix("rd_pstb_r3d.cir")))):
        a = run("multirhs/" + name.replace(" ", "_") + "_off", deck)
        b = run("multirhs/" + name.replace(" ", "_") + "_on", deck, {"PAC_MULTIRHS": "1"})
        record(f"multirhs: {name} curves byte-identical, column-by-column vs one multi-RHS solve",
               "pass" if a["rc"] == b["rc"] == 0 and same_files(a, b) else "FAIL", f"files {len(a['files'])}")
    if not QUICK:
        deck = klu_deck(fix("r7_pstb.cir"))
        a = run("multirhs/r7_pstb_klu_off", deck)
        b = run("multirhs/r7_pstb_klu_on", deck, {"PAC_MULTIRHS": "1"})
        record("multirhs: BSIM4 R7 .pstb under KLU byte-identical", "pass" if a["rc"] == b["rc"] == 0 and same_files(a, b) else "FAIL",
               f"wall {a['wall']:.0f}s off, {b['wall']:.0f}s on")
    # red proof: the second run uses another probe frequency
    d2 = klu_deck(fix("rd_pstb_r3d.cir")).replace("dec 5 500 200000", "dec 5 500 150000")
    a = run("multirhs/r3d_pstb_klu_ref", klu_deck(fix("rd_pstb_r3d.cir")), {"PAC_MULTIRHS": "1"})
    c = run("multirhs/red", d2, {"PAC_MULTIRHS": "1"})
    record("multirhs: red proof (other sweep) must break the byte identity", "pass" if c["files"] and not same_files(a, c) else "FAIL", "")


def main():
    os.makedirs(OUT, exist_ok=True)
    try:
        check_autoklu()
        check_fork()
        check_multirhs()
    except subprocess.TimeoutExpired as e:
        print("COULD-NOT-RUN: timeout", e)
        sys.exit(2)
    fails = [n for n, s in rows if s == "FAIL"]
    print(f"\n{len(rows)} checks, {len(fails)} FAIL")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
