#!/usr/bin/env python3
"""Card #466: native `.stb` frequency-sweep parameter validation, both forms.

An independent reader (verify-v-g4, 2026-09-29) found that native `.stb` /
the interactive `stb` command accepted fstart<=0, fstop<=0, fstart>fstop, and
n<=0 with only a generic "assumes default parameter(s)" warning, even though
stbsetp.c already had an E_PARMVAL check ("Frequency of <= 0 is invalid for
STB start/stop") that could never fire -- dot_stb (inp2dot.c) was silently
substituting a valid default value for anything it judged invalid BEFORE ever
calling the analysis's own setter, so the setter's check never saw the bad
value. Fixed by only defaulting a field that is truly OMITTED (no token
there at all) and otherwise passing the user's actual value straight through
to STBsetParm/STBaskQuest (stbsetp.c), which now also validates fstop<fstart
and numsteps<=0 (neither had a check before). This script is the
COULD-NOT-RUN-proof for that fix: every bad-value case below must produce NO
loopgain vector at all (not a wrong number), in BOTH the interactive `stb`
command and the netlist-card `.stb` line.

While investigating the netlist-card form, also found (and fixed,
src/frontend/inpcom.c's inp_stb_autoquote()) that probe=/probe2=/mode=/
localgnd= all needed manual quoting on a netlist card (numparam tries to
evaluate the unquoted bareword RHS as a number and fails with "Undefined
parameter [...]", reproduced identically on a bogus, unrecognized dot-card --
not something specific to `.stb`) -- dot_stb now auto-quotes these before
numparam ever sees the line, so probe=Vprobe works unquoted on a netlist
card too, matching the interactive command (which was never run through
numparam to begin with). Checked here against the quoted-equivalent form,
which must give an identical T (both forms parse to the same C string once
past this preprocessing step).

Every check is pass / FAIL / COULD-NOT-RUN (never affects the exit code
except as a real failure -- there is no "skip" state used in this file).

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_stb_cli.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import stb  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/stb_cli")
TOL = 1e-12  # same-text-after-preprocessing comparisons (quoted vs auto-quoted)

results = []


def record(name, ok, detail):
    state = "pass" if ok is True else ("FAIL" if ok is False else "COULD-NOT-RUN")
    results.append((name, state, detail))


# A minimal deck, not a real feedback loop (Vin drives E1 directly, nothing
# feeds back into "in") -- fine for every check here, since none of them
# depend on T having any particular closed-form value except the
# quoted-vs-unquoted-text-equivalence ones, which only need the SAME T on
# both sides of that comparison, whatever it is.
CIRCUIT_BODY = """* card #466 stb cli test fixture
Vin in 0 DC 0 AC 1
E1 mid 0 in 0 100
Vprobe mid fb DC 0
R1 fb 0 1k
"""

# A second, differential-capable fixture for probe2=/mode= checks.
CIRCUIT_BODY_DM = """* card #466 stb cli test fixture (dm/cm)
Vin in 0 DC 0 AC 1
Vp in a DC 0
Vn in b DC 0
R1 a 0 1k
R2 b 0 1k
"""

# A third fixture with a reference node ("lref") distinct from the probe's
# own two nodes, for localgnd= -- reusing one of the probe's own nodes
# (e.g. "fb") is a degenerate localgnd choice (found while writing this
# check: it produces NaN in the loopgain), not a meaningful comparison.
CIRCUIT_BODY_LOCALGND = """* card #466 stb cli test fixture (localgnd)
Vin in 0 DC 0 AC 1
E1 mid 0 in 0 100
Vprobe mid fb DC 0
R1 fb lref 1k
Rlg lref 0 1k
"""


def _write(text, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def _run(path, outdir):
    try:
        proc = subprocess.run([NGSPICE, "-b", path], capture_output=True, text=True, timeout=120)
    except (FileNotFoundError, PermissionError) as e:
        raise stb.CouldNotRun(f"cannot start {NGSPICE}: {e}") from e
    with open(os.path.join(outdir, "log.txt"), "w") as fh:
        fh.write(proc.stdout + proc.stderr)
    return proc


def run_interactive(stb_args, outdir, body=CIRCUIT_BODY):
    """`.control ... stb <args> ...` form. Returns dict on success; raises
    stb.CouldNotRun (caller decides whether that is expected) on any
    ngspice-level error or missing output."""
    os.makedirs(outdir, exist_ok=True)
    r1 = os.path.join(outdir, "out.txt")
    if os.path.exists(r1):
        os.remove(r1)
    text = f"""{body}.control
set wr_singlescale
set wr_vecnames
set numdgt=15
stb {stb_args}
wrdata {r1} stb1.loopgain
quit
.endc
.end
"""
    path = _write(text, outdir, "interactive.cir")
    proc = _run(path, outdir)
    log = proc.stdout + proc.stderr
    errs = [ln for ln in log.splitlines() if re.search(r"\berror\b", ln, re.I)]
    produced = os.path.exists(r1)
    if proc.returncode != 0 or errs or not produced:
        raise stb.CouldNotRun(f"ngspice exit {proc.returncode}; produced={produced}; {errs[:3]}")
    freqs, (lg,) = stb.read_wrdata(r1, 1)
    return {"freqs": freqs, "loopgain": lg, "log": log}


def run_netlistcard(stb_card, outdir, body=CIRCUIT_BODY):
    """`.stb ...` as a netlist card, `.control/run/wrdata` for output --
    same recipe test_stb_native.py's check_netlist_card_form() uses."""
    os.makedirs(outdir, exist_ok=True)
    r1 = os.path.join(outdir, "out.txt")
    if os.path.exists(r1):
        os.remove(r1)
    text = f"""{body}{stb_card}
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
run
wrdata {r1} stb1.loopgain
quit
.endc
.end
"""
    path = _write(text, outdir, "netlistcard.cir")
    proc = _run(path, outdir)
    log = proc.stdout + proc.stderr
    errs = [ln for ln in log.splitlines() if re.search(r"\berror\b", ln, re.I)]
    produced = os.path.exists(r1)
    if proc.returncode != 0 or errs or not produced:
        raise stb.CouldNotRun(f"ngspice exit {proc.returncode}; produced={produced}; {errs[:3]}")
    freqs, (lg,) = stb.read_wrdata(r1, 1)
    return {"freqs": freqs, "loopgain": lg, "log": log}


def expect_refused(run_fn, args_or_card, outdir, must_contain):
    """Call run_fn(args_or_card, outdir); pass only if it raised
    CouldNotRun (no vector produced) AND the ngspice log (read back from
    outdir/log.txt even though the exception itself does not carry it)
    contains `must_contain` -- the COULD-NOT-RUN-proof: a wrong number
    would make run_fn return normally (FAIL), not raise."""
    try:
        got = run_fn(args_or_card, outdir)
    except stb.CouldNotRun:
        pass
    else:
        record(f"{outdir}: refused, no loopgain vector produced", False,
               f"ngspice ran anyway and returned data (BUG): loopgain[0]={got['loopgain'][0]!r}")
        return
    logpath = os.path.join(outdir, "log.txt")
    log = open(logpath).read() if os.path.exists(logpath) else ""
    record(f"{outdir}: refused, no loopgain vector produced", True, "COULD-NOT-RUN as expected")
    record(f"{outdir}: refusal message contains {must_contain!r}", must_contain in log,
           f"log tail: {log.strip().splitlines()[-4:]}")


# ---------------------------------------------------------------------
# (3, CODE BUG) fstart<=0, fstop<=0, fstart>fstop, n<=0 must all REFUSE
# (error, no analysis, no vectors), in both forms.
# ---------------------------------------------------------------------

BAD_CASES = [
    # (label, dec_args, expected E_PARMVAL substring)
    ("fstart=0", "dec 10 0 1e6", "Frequency of <= 0 is invalid for STB start"),
    ("fstart=-1", "dec 10 -1 1e6", "Frequency of <= 0 is invalid for STB start"),
    ("fstop=0", "dec 10 1 0", "Frequency of <= 0 is invalid for STB stop"),
    ("fstart>fstop", "dec 10 5000 1", "STB stop frequency must not be less than STB start frequency"),
    ("n=0", "dec 0 1 1e6", "Number of points must be > 0 for STB numsteps"),
    ("n=-5", "dec -5 1 1e6", "Number of points must be > 0 for STB numsteps"),
]


def check_bad_values_refused():
    for label, args, msg in BAD_CASES:
        outdir = os.path.join(DATA_ROOT, "refuse_interactive", label.replace(" ", "_").replace(">", "gt").replace("=", "eq"))
        expect_refused(lambda a, o: run_interactive(a + " probe=Vprobe", o), args, outdir, msg)

        outdir2 = os.path.join(DATA_ROOT, "refuse_netlistcard", label.replace(" ", "_").replace(">", "gt").replace("=", "eq"))
        card = f'.stb {args} probe="Vprobe"\n'
        expect_refused(lambda c, o: run_netlistcard(c, o), card, outdir2, msg)


# ---------------------------------------------------------------------
# Good values must be unaffected by the fix (regression), and omitted
# values must still fall back to defaults with the generic warning
# (unchanged behaviour -- only an explicitly-bad value now refuses).
# ---------------------------------------------------------------------

def check_good_values_unchanged():
    # dec N FSTART FSTOP sweeps N points PER DECADE (like .ac): 1 Hz to
    # 1e6 Hz is 6 decades, so 10*6+1 = 61 points total -- not 10.
    outdir = os.path.join(DATA_ROOT, "good_interactive")
    got = run_interactive("dec 10 1 1e6 probe=Vprobe", outdir)
    record("good values (interactive): 61 points (10/decade x 6 decades + 1)",
           len(got["freqs"]) == 61, f"got {len(got['freqs'])}")
    record("good values (interactive): fstart=1", abs(got["freqs"][0] - 1.0) < 1e-9, f"got {got['freqs'][0]!r}")
    record("good values (interactive): fstop=1e6", abs(got["freqs"][-1] - 1e6) / 1e6 < 1e-9, f"got {got['freqs'][-1]!r}")
    record("good values (interactive): no default-parameter warning",
           "assumes default parameter" not in got["log"], "checked log")

    outdir2 = os.path.join(DATA_ROOT, "good_netlistcard")
    got2 = run_netlistcard('.stb dec 10 1 1e6 probe="Vprobe"\n', outdir2)
    record("good values (netlist card): 61 points (10/decade x 6 decades + 1)",
           len(got2["freqs"]) == 61, f"got {len(got2['freqs'])}")
    record("good values (netlist card): no default-parameter warning",
           "assumes default parameter" not in got2["log"], "checked log")


def check_missing_still_defaults():
    """n, fstart AND fstop all omitted: still default (n->10/decade,
    fstart->1, fstop->1000x fstart = 1000, i.e. 3 decades -> 10*3+1 = 31
    points) with the generic warning, and the run still completes -- this
    is NOT the "invalid explicit value" case the fix changed; genuinely
    omitted fields are unaffected."""
    outdir = os.path.join(DATA_ROOT, "missing_n_interactive")
    got = run_interactive("dec probe=Vprobe", outdir)
    record("n/fstart/fstop omitted (interactive): still runs (31 default points)",
           len(got["freqs"]) == 31, f"got {len(got['freqs'])}")
    record("n/fstart/fstop omitted (interactive): default-parameter warning IS printed",
           "assumes default parameter" in got["log"], "checked log")

    outdir2 = os.path.join(DATA_ROOT, "missing_n_netlistcard")
    got2 = run_netlistcard('.stb dec probe="Vprobe"\n', outdir2)
    record("n/fstart/fstop omitted (netlist card): still runs (31 default points)",
           len(got2["freqs"]) == 31, f"got {len(got2['freqs'])}")
    record("n/fstart/fstop omitted (netlist card): default-parameter warning IS printed",
           "assumes default parameter" in got2["log"], "checked log")


# ---------------------------------------------------------------------
# (2) netlist-card auto-quoting: probe=/probe2=+mode=/localgnd= must now
# work UNQUOTED, matching the already-quoted form exactly (T identical).
# ---------------------------------------------------------------------

def compare_loopgain(name, a, b, tol=TOL):
    """Both sides here are the SAME computation on two netlist texts that
    the fix makes byte-identical after preprocessing (an auto-quoted
    bareword and a hand-quoted one), so an absolute-or-relative bound
    (not stbcore.max_rel_err's pure-relative one, which divides by zero
    when this fixture's open-loop T is exactly 0 at some frequencies) is
    the right comparator: both values should match near machine
    precision, whatever their magnitude, including exactly 0."""
    worst = max(abs(x - y) / max(abs(y), 1e-300) if abs(y) > 1e-12 else abs(x - y)
                 for x, y in zip(a, b, strict=True))
    record(name, worst <= tol, f"max abs-or-rel diff {worst:.3e} (bound {tol:.0e})")


def check_autoquote_probe():
    quoted = run_netlistcard('.stb dec 10 1 1e6 probe="Vprobe"\n', os.path.join(DATA_ROOT, "autoquote_probe_quoted"))
    unquoted = run_netlistcard('.stb dec 10 1 1e6 probe=Vprobe\n', os.path.join(DATA_ROOT, "autoquote_probe_unquoted"))
    compare_loopgain("netlist card probe=Vprobe unquoted == probe=\"Vprobe\" quoted",
                      unquoted["loopgain"], quoted["loopgain"])


def check_autoquote_probe2_mode():
    quoted = run_netlistcard('.stb dec 10 1 1e6 probe="Vp" probe2="Vn" mode="dm"\n',
                              os.path.join(DATA_ROOT, "autoquote_dmcm_quoted"), body=CIRCUIT_BODY_DM)
    unquoted = run_netlistcard('.stb dec 10 1 1e6 probe=Vp probe2=Vn mode=dm\n',
                                os.path.join(DATA_ROOT, "autoquote_dmcm_unquoted"), body=CIRCUIT_BODY_DM)
    compare_loopgain("netlist card probe=/probe2=/mode= unquoted == all quoted",
                      unquoted["loopgain"], quoted["loopgain"])


def check_autoquote_localgnd():
    body = CIRCUIT_BODY_LOCALGND
    quoted = run_netlistcard('.stb dec 10 1 1e6 probe="Vprobe" localgnd="lref"\n',
                              os.path.join(DATA_ROOT, "autoquote_localgnd_quoted"), body=body)
    unquoted = run_netlistcard('.stb dec 10 1 1e6 probe=Vprobe localgnd=lref\n',
                                os.path.join(DATA_ROOT, "autoquote_localgnd_unquoted"), body=body)
    compare_loopgain("netlist card localgnd=lref unquoted == localgnd=\"lref\" quoted",
                      unquoted["loopgain"], quoted["loopgain"])


# ---------------------------------------------------------------------
# Regression: the rest of the suite must stay green against this build.
# ---------------------------------------------------------------------

def check_regressions():
    env = dict(os.environ)
    env["NGSPICE"] = NGSPICE
    for name in ("test_stb_native2.py", "test_stb_native.py", "test_stb_device.py",
                 "test_stb_followups.py", "test_stb.py"):
        proc = subprocess.run([sys.executable, os.path.join(HERE, name)],
                               capture_output=True, text=True, env=env, timeout=1800)
        tail = "\n".join(proc.stdout.splitlines()[-3:])
        record(f"{name} with NGSPICE=<new build>", proc.returncode == 0, f"exit {proc.returncode}: {tail}")


def main():
    try:
        check_bad_values_refused()
        check_good_values_unchanged()
        check_missing_still_defaults()
        check_autoquote_probe()
        check_autoquote_probe2_mode()
        check_autoquote_localgnd()
    except stb.CouldNotRun as e:
        record("run", None, str(e))

    check_regressions()

    width = max(len(n) for n, _, _ in results)
    for name, state, detail in results:
        print(f"{state:16s} {name:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    npass = sum(s == "pass" for _, s, _ in results)
    print(f"\n{npass}/{len(results)} pass")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
