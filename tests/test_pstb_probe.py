#!/usr/bin/env python3
"""Validation for card #717: `.pstb` probe forms at parity with `.stb`.

`.pstb` used to take only `probe="Vname"` (a plain top-level 0 V source). It now also takes

  probe="<inst>:<terminal>"   auto-insert a 0 V source in series with that instance terminal
                              (the same parse-time splice `.stb` uses)
  probe="<Xname>"             an instance of the lgprobe subcircuit (stb/lgprobe.lib)

Checks (every comparison is against the SAME circuit's plain-probe `.pstb` run):
  1  linear loop (tests/fixtures/pstb_linear.cir), probe=inst:terminal  == plain probe, 1e-9 relative
  2  linear loop, probe=<lgprobe instance>                              == plain probe, 1e-9 relative
  3  linear loop, probe=inst:terminal on .pstb == the same probe on native .stb, 1e-6 relative
  4  switching converter (refs/r3_vm_buck.cir, closed loop), probe=inst:terminal == plain probe, 1e-9
  5  switching converter, the deck's own lgprobe instance (XP)                   == plain probe, 1e-9
  6  refusals: unknown instance, unknown terminal, non-V-source probe, comma (differential) form,
     a name that is neither a V-source nor an lgprobe instance -- each with a clear message, no data
  7  red proof: the same inst:terminal and lgprobe decks on the baseline binary (s140, NGSPICE_BASE)
     do NOT produce a loop gain (refused); and the sign-flip switch makes the lgprobe form disagree

    NGSPICE=<prefix>/bin/ngspice python3 tests/test_pstb_probe.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
Decks and outputs live under ~/Work/data/ngspice_augmentation/717 (never the repository).
"""
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXTURES = os.path.join(HERE, "fixtures")
NGSPICE = os.environ.get("NGSPICE", "ngspice")
NGSPICE_BASE = os.environ.get(
    "NGSPICE_BASE",
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-s140/bin/ngspice"))
DATA = os.path.expanduser("~/Work/data/ngspice_augmentation/717/test")
LGLIB = os.path.join(ROOT, "stb", "lgprobe.lib")

results = []


class CouldNotRun(Exception):
    pass


def record(name, state, detail):
    results.append((name, state, detail))
    print(f"{state:14s} {name}\n{'':14s}   {detail}", flush=True)


def run(binary, deck, env_extra=None, timeout=900):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "2"
    if env_extra:
        env.update(env_extra)
    t0 = time.time()
    try:
        p = subprocess.run([binary, "-b", deck], capture_output=True, text=True, timeout=timeout,
                           cwd=os.path.dirname(deck), env=env)
    except FileNotFoundError as e:
        raise CouldNotRun(f"cannot start {binary}: {e}")
    except subprocess.TimeoutExpired:
        raise CouldNotRun(f"{deck} timed out after {timeout}s")
    return p.returncode, p.stdout + p.stderr, time.time() - t0


def read_wr(path):
    """freqs, loopgain (complex list) from a wr_singlescale/wr_vecnames file, 'loopgain' first."""
    if not os.path.exists(path):
        return None, None
    f, lg = [], []
    with open(path) as fh:
        next(fh)
        for row in fh:
            v = row.split()
            if len(v) >= 3:
                f.append(float(v[0]))
                lg.append(complex(float(v[1]), float(v[2])))
    return f, lg


def write_deck(name, body, card, vec="loopgain"):
    os.makedirs(DATA, exist_ok=True)
    deck = os.path.join(DATA, name + ".cir")
    out = os.path.join(DATA, name + ".txt")
    if os.path.exists(out):
        os.remove(out)
    with open(deck, "w") as fh:
        fh.write(body.rstrip("\n") + "\n" + card + "\n")
        fh.write(".control\nset num_threads=2\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                 f"run\nwrdata {name}.txt {vec}\nquit\n.endc\n.end\n")
    return deck, out


def worst_rel(a, b):
    if a is None or b is None or len(a) != len(b):
        raise CouldNotRun("missing output or point-count mismatch")
    return max(abs(x - y) / abs(x) for x, y in zip(a, b))


# ---------------------------------------------------------------- decks
def linear_bodies():
    with open(os.path.join(FIXTURES, "pstb_linear.cir")) as fh:
        src = "".join(l for l in fh if l.strip().lower() != ".end")
    assert "Vprobe fb vn DC 0\n" in src
    plain = src
    # the probe wire merged away (fb == vn): the loop is then broken by probe=R1:R- instead
    inst = (src.replace("Vprobe fb vn DC 0\n", "").replace("G1 0 o vp vn", "G1 0 o vp fb")
            .replace("Cn vn 0", "Cn fb 0"))
    lg = src.replace("Vprobe fb vn DC 0", "XP fb vn lgprobe") + f".include {LGLIB}\n"
    return plain, inst, lg, plain


def r3_bodies():
    with open(os.path.join(ROOT, "refs", "r3_vm_buck.cir")) as fh:
        lines = fh.read().splitlines()
    keep = []
    for l in lines:
        if l.strip().lower() == ".control":
            break
        keep.append(l)
    src = "\n".join(keep) + "\n"
    xp = "XP fbtop vout lgprobe"
    assert xp in src and "Ra fbtop fbnode {ra}" in src
    inc = ".include ../stb/lgprobe.lib"
    assert inc in src
    # The plain-probe reference writes its Vprobe LAST in the deck, the order in which the parser
    # creates a spliced probe source: moving the probe line elsewhere in an otherwise identical deck
    # changes the PSS orbit/loop gain by ~1e-7 (node and pivot order, shooting tolerance), so the
    # 1e-9 comparison is only meaningful at equal element order (check_r3_order_floor reports that
    # floor).
    plain = src.replace(inc + "\n", "").replace(xp + "\n", "") + "Vprobe fbtop vout DC 0\n"
    inplace = src.replace(inc + "\n", "").replace(xp, "Vprobe fbtop vout DC 0")
    inst = src.replace(inc + "\n", "").replace(xp + "\n", "").replace("Ra fbtop fbnode {ra}",
                                                                      "Ra vout fbnode {ra}")
    lg = src.replace(inc, f".include {LGLIB}")
    return plain, inst, lg, inplace


LIN_CARD = ".pstb 100k 0 clkgnd 2000 1 5 0.001 dec 3 1k 200k probe=\"{p}\" method=\"shoot\""
R3_CARD = ".pstb 500e3 300e-6 vout 200 10 150 1e-6 dec 3 1k 100k probe=\"{p}\" method=\"shoot\""


def run_form(binary, tag, body, card, probe, env_extra=None):
    deck, out = write_deck(tag, body, card.format(p=probe))
    rc, text, wall = run(binary, deck, env_extra)
    f, lg = read_wr(out)
    return rc, text, wall, f, lg


# ---------------------------------------------------------------- checks
def check_forms(prefix, label, bodies, card, plain_probe, forms):
    """forms: (key, body, probe, reference body). The reference is the plain-probe deck with the
    probe source at the same place in the element order as the form under test (bodies[0] for the
    spliced probe, which the parser creates last; bodies[3] for a probe replacing XP in place)."""
    refs = {}
    for key, body, probe, refbody in forms:
        nm = f"{label}: probe={probe} vs plain probe={plain_probe} (1e-9 relative)"
        try:
            rk = "ref_" + str(id(refbody))
            if rk not in refs:
                _, text, w0, f0, lg0 = run_form(NGSPICE, f"{prefix}_plain_{len(refs)}", refbody, card,
                                                plain_probe)
                if lg0 is None:
                    raise CouldNotRun(f"plain-probe reference run gave no data: {text[-1500:]}")
                refs[rk] = (f0, lg0, w0)
            f0, lg0, w0 = refs[rk]
            rc, text, w, f, v = run_form(NGSPICE, f"{prefix}_{key}", body, card, probe)
            if v is None:
                record(nm, "FAIL", f"no loop gain produced; last output: {text.strip()[-300:]!r}")
                continue
            if len(f) != len(f0) or any(abs(a - b) > 1e-9 * a for a, b in zip(f0, f)):
                record(nm, "FAIL", "frequency grid differs")
                continue
            e = worst_rel(lg0, v)
            record(nm, "pass" if e <= 1e-9 else "FAIL",
                   f"{len(f)} points, worst relative |dT| = {e:.3e}, wall {w:.1f}s (plain {w0:.1f}s)")
        except CouldNotRun as ex:
            record(nm, "COULD-NOT-RUN", str(ex))


def check_linear_vs_stb():
    nm = "linear: .pstb probe=R1:R- vs native .stb probe=R1:R- (1e-6 relative)"
    try:
        _, inst, _, _ = linear_bodies()
        deck, out = write_deck("lin_stb", inst, '.stb dec 3 1k 200k probe="R1:R-"',
                               vec="stb1.loopgain")
        rc, text, _ = run(NGSPICE, deck)
        f1, v1 = read_wr(out)
        if v1 is None:
            raise CouldNotRun(f".stb reference gave no data: {text[-800:]}")
        rc, text, w, f2, v2 = run_form(NGSPICE, "lin_stb_p", inst, LIN_CARD, "R1:R-")
        if v2 is None:
            record(nm, "FAIL", "no .pstb loop gain")
            return
        e = worst_rel(v1, v2)
        record(nm, "pass" if e <= 1e-6 else "FAIL", f"{len(f1)} points, worst relative = {e:.3e}")
    except CouldNotRun as ex:
        record(nm, "COULD-NOT-RUN", str(ex))


def check_r3_order_floor():
    nm = "info: r3 plain probe, Vprobe line moved within the deck (reordering floor of the PSS orbit)"
    try:
        b = r3_bodies()
        _, _, _, f0, v0 = run_form(NGSPICE, "r3_order_a", b[0], R3_CARD, "Vprobe")
        _, _, _, f1, v1 = run_form(NGSPICE, "r3_order_b", b[3], R3_CARD, "Vprobe")
        record(nm, "info", f"worst relative |dT| = {worst_rel(v0, v1):.3e} between two orderings of the same circuit")
    except CouldNotRun as ex:
        record(nm, "COULD-NOT-RUN", str(ex))


def check_refusals():
    plain, inst, lg, _ = linear_bodies()
    cases = [
        ("unknown instance", plain, "Nosuch:R+", r"not in circuit"),
        ("unknown terminal", plain, "R1:zz", r"has no terminal"),
        ("not a V-source and not an lgprobe", plain, "R2", r"not a V-source|neither"),
        ("unknown name", plain, "Qnone", r"neither|not a V-source"),
        ("differential comma form", plain, "Vprobe,Vprobe2", r"not supported"),
    ]
    for what, body, probe, pat in cases:
        nm = f"refused with a message: probe={probe} ({what})"
        try:
            rc, text, w, f, v = run_form(NGSPICE, "ref_" + re.sub(r"\W", "_", probe), body, LIN_CARD, probe)
            if v is None and re.search(pat, text, re.I):
                m = [l for l in text.splitlines() if re.search(pat, l, re.I)][0].strip()
                record(nm, "pass", m[:200])
            else:
                record(nm, "FAIL", f"data produced={v is not None}; no matching message ({pat}); "
                       f"tail {text.strip()[-300:]!r}")
        except CouldNotRun as ex:
            record(nm, "COULD-NOT-RUN", str(ex))


def check_redproof():
    plain, inst, lg, _ = linear_bodies()
    # (a) the baseline binary cannot run the new forms (a refusal, not a number)
    for key, body, probe in (("inst", inst, "R1:R-"), ("lg", lg, "XP")):
        nm = f"red proof: baseline binary does not run probe={probe}"
        try:
            if not os.path.exists(NGSPICE_BASE):
                raise CouldNotRun(f"baseline binary {NGSPICE_BASE} not found")
            if os.path.realpath(NGSPICE_BASE) == os.path.realpath(NGSPICE):
                record(nm, "info", "NGSPICE is the baseline itself; see the other checks (they FAIL)")
                continue
            rc, text, w, f, v = run_form(NGSPICE_BASE, f"red_{key}", body, LIN_CARD, probe)
            if v is None:
                last = [l for l in text.splitlines() if "probe" in l.lower()][-1:]
                record(nm, "pass", f"baseline refuses ({last[0].strip()[:150] if last else 'no data'})")
            else:
                record(nm, "FAIL", "baseline produced data; the new checks would not be red")
        except CouldNotRun as ex:
            record(nm, "COULD-NOT-RUN", str(ex))
    # (b) a sign error in the lgprobe wiring is detected by the comparison
    nm = "red proof: PAC_REDPROOF_SIGNFLIP makes the lgprobe form disagree with the plain probe"
    try:
        rc, text, w, f0, v0 = run_form(NGSPICE, "red_plain", plain, LIN_CARD, "Vprobe")
        rc, text, w, f1, v1 = run_form(NGSPICE, "red_lgflip", lg, LIN_CARD, "XP",
                                       {"PAC_REDPROOF_SIGNFLIP": "1"})
        if v0 is None or v1 is None:
            raise CouldNotRun("no data")
        e = worst_rel(v0, v1)
        record(nm, "pass" if e > 1e-3 else "FAIL", f"worst relative difference {e:.3e} (must be > 1e-3)")
    except CouldNotRun as ex:
        record(nm, "COULD-NOT-RUN", str(ex))


def main():
    os.makedirs(DATA, exist_ok=True)
    lin = linear_bodies()
    check_forms("lin", "linear", lin, LIN_CARD, "Vprobe",
                [("inst", lin[1], "R1:R-", lin[0]), ("lg", lin[2], "XP", lin[3])])
    check_linear_vs_stb()
    r3 = r3_bodies()
    check_forms("r3", "r3_vm_buck", r3, R3_CARD, "Vprobe",
                [("inst", r3[1], "Ra:R+", r3[0]), ("lg", r3[2], "XP", r3[3])])
    check_r3_order_floor()
    check_refusals()
    check_redproof()

    npass = sum(1 for r in results if r[1] == "pass")
    nfail = sum(1 for r in results if r[1] == "FAIL")
    ncnr = sum(1 for r in results if r[1] == "COULD-NOT-RUN")
    ntot = sum(1 for r in results if r[1] != "info")
    print(f"\n{npass}/{ntot} pass" + (f", {nfail} FAIL" if nfail else "") +
          (f", {ncnr} COULD-NOT-RUN" if ncnr else ""))
    return 1 if nfail else (2 if ncnr else 0)


if __name__ == "__main__":
    sys.exit(main())
