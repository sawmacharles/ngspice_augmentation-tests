#!/usr/bin/env python3
"""Smoke test for tools/opbench.py (card #637, DC OP benchmark harness).

Runs the harness on three small corpus decks (latch_nmos_res, jfet_isrc,
casc_selfbias) plus a deliberately singular deck (two voltage sources across
one node), written to a temp corpus dir.  Checks (pass / FAIL / COULD-NOT-RUN):

  1. parse_log() on canned logs: fallback stage detection (newton, dyngmin,
     src, optran), 'Warning: singular matrix' alone is not a failure,
     'could not be simulated' is, non-finite values are.
  2. default run (+ 2 nodeset seeds): good decks converged, iterations > 0,
     solution file present with the expected nodes, reference comparison
     filled in and 'ok' (maxerr < 1e-6), casc_selfbias reports a fallback
     stage, nodeset variants carry their seed and are ok/other (not WRONG);
     the singular deck is not converged in every variant and has no solution.
  3. summary.json counts agree with results.csv.
  4. Newton-only run (--options "gminsteps=0 srcsteps=0" --no-optran):
     casc_selfbias no longer converges, latch still does; `compare` prints
     the lost deck.

Usage:  NGSPICE=<fork binary> python3 tests/test_opbench.py
Outputs go to ~/Work/data/ngspice_augmentation/637/test_opbench/, never the repo.
Exit status: 0 = all pass, 1 = a check FAILed, 2 = could not run.
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import opbench  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
OUT = os.path.expanduser("~/Work/data/ngspice_augmentation/637/test_opbench")
results = []


def record(name, state, detail=""):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN")
    results.append((name, state, detail))


def check(name, cond, detail=""):
    record(name, "pass" if cond else "FAIL", detail)


def rows_of(d):
    with open(os.path.join(d, "results.csv")) as f:
        return {(r["deck"], r["variant"]): r for r in csv.DictReader(f)}


def run_harness(corpus_dir, label, extra):
    cmd = [sys.executable, os.path.join(ROOT, "tools", "opbench.py"), "--ngspice", NGSPICE, "--label", label,
           "--no-refs", "--corpus-dir", corpus_dir, "--out", os.path.join(OUT, label),
           "--refcache", os.path.join(OUT, "refcache")] + extra
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    return p


def test_parser():
    good = "Total iterations = 4\nCircuit Equations = 3\na = 1.0e+00\nb = 5e-01\n"
    r = opbench.parse_log(good)
    check("parse: plain Newton success", r["converged"] and r["stage"] == "newton" and r["iters"] == 4
          and r["unknowns"] == 3 and r["nodes"] == {"a": 1.0, "b": 0.5}, str(r["stage"]))
    dg = "Warning: singular matrix:  check node x\nNote: Starting dynamic gmin stepping\nNote: Dynamic gmin stepping completed\n" + good
    r = opbench.parse_log(dg)
    check("parse: dynamic gmin, singular warning not fatal", r["converged"] and r["stage"] == "dyngmin", f"{r['stage']} {r['tried']}")
    sr = ("Note: Starting dynamic gmin stepping\nWarning: Dynamic gmin stepping failed\nNote: Starting true gmin stepping\n"
          "Warning: True gmin stepping failed\nNote: Starting source stepping\nNote: Source stepping completed\n" + good)
    r = opbench.parse_log(sr)
    check("parse: source stepping after failed gmin stages", r["stage"] == "src" and r["tried"] == "dyngmin,gmin,src", f"{r['stage']} {r['tried']}")
    ot = sr.replace("Note: Source stepping completed", "Warning: source stepping failed\nNote: Transient op started\nNote: Transient op finished successfully")
    r = opbench.parse_log(ot)
    check("parse: optran", r["stage"] == "optran" and r["tried"].endswith("optran"), f"{r['stage']} {r['tried']}")
    bad = "Error: The operating point could not be simulated successfully.\nTotal iterations = 101\n"
    r = opbench.parse_log(bad)
    check("parse: failed OP", not r["converged"] and r["status"] == "noconv", r["status"])
    r = opbench.parse_log(good.replace("5e-01", "nan"))
    check("parse: non-finite value is not converged", not r["converged"])


def main():
    try:
        subprocess.run([NGSPICE, "--version"], capture_output=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as e:
        record("ngspice available", "COULD-NOT-RUN", f"{NGSPICE}: {e}")
        return finish()
    if not os.path.isdir(os.path.join(ROOT, "tests", "fixtures", "op")):
        record("corpus present", "COULD-NOT-RUN", "tests/fixtures/op missing")
        return finish()
    os.makedirs(OUT, exist_ok=True)
    test_parser()

    tmp = tempfile.mkdtemp(prefix="opbench_corpus_", dir=OUT)
    for d in ("latch_nmos_res", "jfet_isrc", "casc_selfbias"):
        shutil.copy(os.path.join(ROOT, "tests", "fixtures", "op", d + ".cir"), tmp)
    with open(os.path.join(tmp, "singular.cir"), "w") as f:
        f.write("* two ideal sources across one node: singular matrix\nV1 a 0 1\nV2 a 0 2\nR1 a 0 1k\n.op\n.end\n")

    p = run_harness(tmp, "default", ["--nodeset-seeds", "2", "--jobs", "2"])
    if p.returncode != 0:
        record("harness default run", "COULD-NOT-RUN", (p.stdout + p.stderr)[-400:])
        return finish()
    R = rows_of(os.path.join(OUT, "default"))
    check("rows: 4 decks x (plain + 2 seeds)", len(R) == 12, str(len(R)))
    for d, node in (("latch_nmos_res", "a"), ("jfet_isrc", "s")):
        r = R[(d, "plain")]
        sol_path = os.path.join(OUT, "default", "solutions", f"{d}__plain.json")
        sol = json.load(open(sol_path))["nodes"] if os.path.exists(sol_path) else {}
        check(f"{d}: converged, iterations > 0, solution present",
              r["converged"] == "1" and r["iters"] != "" and int(r["iters"]) > 0 and node in sol and int(r["unknowns"]) > 0,
              f"conv={r['converged']} iters={r['iters']} nodes={sorted(sol)[:4]}")
        check(f"{d}: reference comparison", r["ref_ok"] == "1" and r["flag"] == "ok" and r["maxerr"] != "" and float(r["maxerr"]) < 1e-6,
              f"flag={r['flag']} maxerr={r['maxerr']}")
    r = R[("casc_selfbias", "plain")]
    check("casc_selfbias: converged through a fallback stage",
          r["converged"] == "1" and r["stage"] in ("dyngmin", "gmin", "src", "optran") and r["tried"] != "",
          f"stage={r['stage']} tried={r['tried']}")
    ns = R[("latch_nmos_res", "ns1")]
    check("nodeset variant: seed recorded, not WRONG", ns["seed"] == "637:latch_nmos_res:1" and ns["flag"] in ("ok", "other")
          and ns["converged"] == "1", f"seed={ns['seed']} flag={ns['flag']}")
    check("singular deck: never converged, no solution",
          all(R[("singular", v)]["converged"] == "0" for v in ("plain", "ns1", "ns2"))
          and not os.path.exists(os.path.join(OUT, "default", "solutions", "singular__plain.json")),
          str({v: R[("singular", v)]["status"] for v in ("plain", "ns1", "ns2")}))
    summ = json.load(open(os.path.join(OUT, "default", "summary.json")))
    nconv = sum(r["converged"] == "1" for r in R.values() if r["variant"] == "plain")
    check("summary.json agrees with results.csv", summ["plain"]["runs"] == 4 and summ["plain"]["converged"] == nconv == 3,
          f"summary plain={summ['plain']['converged']}/{summ['plain']['runs']} csv={nconv}")

    p = run_harness(tmp, "newton", ["--options", "gminsteps=0 srcsteps=0", "--no-optran", "--jobs", "2"])
    if p.returncode != 0:
        record("harness newton run", "COULD-NOT-RUN", (p.stdout + p.stderr)[-400:])
        return finish()
    N = rows_of(os.path.join(OUT, "newton"))
    check("newton-only: latch converges, casc_selfbias does not (no optran either)",
          N[("latch_nmos_res", "plain")]["converged"] == "1" and N[("casc_selfbias", "plain")]["converged"] == "0"
          and N[("latch_nmos_res", "plain")]["stage"] == "newton",
          f"latch={N[('latch_nmos_res','plain')]['stage']} casc={N[('casc_selfbias','plain')]['status']}")
    c = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "opbench.py"), "compare", os.path.join(OUT, "default"),
                        os.path.join(OUT, "newton")], capture_output=True, text=True, timeout=60)
    check("compare: reports the lost deck", c.returncode == 0 and "LOST  casc_selfbias/plain" in c.stdout, c.stdout[-300:])
    shutil.rmtree(tmp, ignore_errors=True)
    return finish()


def finish():
    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:13s} {name:{width}s}  {detail}")
    print(f"\n{sum(s == 'pass' for _, s, _ in results)}/{len(results)} pass")
    states = {s for _, s, _ in results}
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
