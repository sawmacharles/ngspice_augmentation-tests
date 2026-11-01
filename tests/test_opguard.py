#!/usr/bin/env python3
"""Card #642: the OP ladder sets aside a 'converged' operating point whose node voltages are absurd
(beyond 1e6 x max(1 V, largest |V source|), or not finite) and goes on to the next stage; when nothing
else converges, the set-aside result is restored, so no deck ends worse than before.

Checks (fork binary from $NGSPICE):
  1. fixtures/opguard_exp_chain.cir: stock gmin stepping lands near 1e84 V; with the guard the OP is the
     true one (every node within 0..30 V, the chain bottom n5 near one diode drop) and the warning names gmin stepping.
  2. same deck with OP_ABSURD_GUARD=0: the old absurd result comes back (|n0| > 1e30), i.e. the switch restores stock.
  3. fixtures/opguard_highv.cir: a legitimate 1e9 V operating point (current source only): the same values as
     without the guard (OP_ABSURD_GUARD=0), bit for bit in the printed digits, plus the 'keeping' warning.

    NGSPICE=<build>/bin/ngspice python3 tests/test_opguard.py    exit 0 pass, 1 any FAIL, 2 could-not-run
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
NGSPICE = os.environ.get("NGSPICE", "ngspice")


def run(deck, guard=True):
    env = dict(os.environ)
    if not guard:
        env["OP_ABSURD_GUARD"] = "0"
    p = subprocess.run([NGSPICE, "-b", os.path.join(FIX, deck)], capture_output=True, text=True, timeout=300, env=env)
    out = p.stdout + p.stderr
    vals = {m.group(1): float(m.group(2)) for m in re.finditer(r"^(\w+) = ([-+0-9.eE]+)\s*$", out, re.M)}
    return out, vals


def main():
    results = []

    def check(name, ok, detail):
        results.append(ok)
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")

    try:
        out, v = run("opguard_exp_chain.cir")
    except (OSError, subprocess.TimeoutExpired) as e:
        print(f"COULD-NOT-RUN: {NGSPICE}: {e}")
        return 2
    nodes = [k for k in v if re.fullmatch(r"[nt]\d", k)]
    sane = nodes and all(-1e-6 <= v[k] <= 30.0 for k in nodes)
    check("exp chain solved", bool(sane) and 0.5 < v.get("n5", -1) < 1.5,
          f"{len(nodes)} nodes, n0={v.get('n0')}, n5={v.get('n5')}")
    check("exp chain warning", "operating point from gmin stepping has node voltages beyond" in out,
          "warning names gmin stepping")

    out0, v0 = run("opguard_exp_chain.cir", guard=False)
    check("OP_ABSURD_GUARD=0 restores stock", abs(v0.get("n0", 0.0)) > 1e30, f"n0={v0.get('n0')}")

    outh, vh = run("opguard_highv.cir")
    _, vh0 = run("opguard_highv.cir", guard=False)
    same = vh0 and all(vh.get(k) == vh0[k] for k in ("a", "b"))
    check("high-voltage OP unchanged", bool(same) and vh0.get("a", 0) > 1e8, f"a={vh.get('a')} vs {vh0.get('a')}")
    check("high-voltage keep warning", "keeping the operating point" in outh, "restore warning printed")

    n_fail = results.count(False)
    print(f"{results.count(True)}/{len(results)} pass")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
