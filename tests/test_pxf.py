#!/usr/bin/env python3
"""Validation for the time-domain shooting-based periodic transfer function
(`.pxf ... method="shoot" solver="td"`), card #499, private-fork branch
aug-pxf, engine in src/spicelib/analysis/pacshoot.c (pac_td_pxf(), reusing
pac_engine_run() -- see that function's own design comment in pacshoot.c
for the option this pass chose and why: a batched-forward per-source sweep
sharing ONE PSS convergence, not a from-scratch backward/adjoint sweep of
the per-step recursion).

Read docs/pss-pac.md's `.pxf` section (added this card) before changing
anything here -- it has the syntax, the design decision (and the
alternative considered), and the accuracy numbers this file's checks
reproduce.

Scope: `.pxf` vs `.pac ... solver=td` (SAME source, isolated) on three
fixtures (an LTI network, a switched-RC LPTV network, and the R3
voltage-mode buck, smooth switches only per project ruling #447 A) plus a
direct `.pxf` vs plain `.ac` cross-check on the LTI fixture (the card's
own "plus an LTI deck where it equals plain .ac/.tf" requirement). Every
check compares COMPLEX values at every swept frequency, gated at the
card's own 1e-6 relative target.

card #504 (`.pxf ... adjoint=1`, branch `aug-adj`, `pac_engine_run_adjoint()`
in `pacshoot.c`): the from-scratch true backward/transposed engine card
#499's own design comment above deferred -- `check_pxf_adjoint_*()` below
cross-check `adjoint=1` against `adjoint=0` (this file's ORIGINAL v1
engine, still the default for a small source count) on the SAME three
fixtures (both under `PAC_EXACT_JACOBIAN=1` -- see docs/pss-pac.md's
adjoint section for why that is the fair, floating-point-exact comparison)
plus a wall-clock timing sweep (`tests/fixtures/pxf_manysrc.cir`, N=2/10/50
sources) confirming the "K sources for ~1 sweep" design claim, and an ASan
check of the adjoint engine's own new code (`ADJ_ASAN_NGSPICE`).

    NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_pxf.py
        exit 0 all pass, 1 any FAIL, 2 any COULD-NOT-RUN
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# Reuse test_pac.py's own already-proven helpers (NGSPICE resolution,
# ngspice invocation/timeout handling, the wrdata-complex reader, the
# fixture-body reader, the R3 `.pac` fixture-deck builder) rather than
# duplicating them -- both scripts read the SAME NGSPICE env var, at
# import time, from the SAME os.environ, so this is equivalent to having
# defined them locally.
import test_pac as _pac  # noqa: E402

run_ngspice = _pac.run_ngspice
read_wrdata_complex = _pac.read_wrdata_complex
body_lines = _pac.body_lines
CouldNotRun = _pac.CouldNotRun
FIXTURES = _pac.FIXTURES

DATA_ROOT = os.path.expanduser("~/Work/data/ngspice_augmentation/pxf")

# card #499: ASan build of THIS card's own worktree (ngspice-src-pxf/asan/),
# same buildrecipe/README.md "AddressSanitizer build" recipe test_pac.py's
# own ASAN_NGSPICE already uses (CFLAGS -O1 -g -fsanitize=address
# -fno-omit-frame-pointer, ASAN_OPTIONS=detect_leaks=0 during make/install
# and here at runtime -- a build-time helper leaks and unrelated
# XSPICE/cmpp-style leaks in one-shot host tooling are not this card's
# concern).
ASAN_NGSPICE = os.environ.get(
    "ASAN_NGSPICE") or os.environ.get(
    "NGSPICE_ASAN",  # the portfolio's build under test (2026-09-30: the pinned default below was silently used)
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-pxf-1-asan/bin/ngspice"))

# card #504: this card's OWN ASan build (ngspice-src-adj/, branch aug-adj --
# see buildrecipe/README.md's "AddressSanitizer build" recipe, applied to
# THIS worktree so the adjoint engine's own new code -- SMPtSolve, the
# trajectory recorder, pac_engine_run_adjoint() -- is what is under test,
# not an older worktree's snapshot of pacshoot.c).
ADJ_ASAN_NGSPICE = os.environ.get(
    "ADJ_ASAN_NGSPICE") or os.environ.get(
    "NGSPICE_ASAN",  # the portfolio's build under test (2026-09-30: the pinned default below was silently used)
    os.path.expanduser("~/Work/apps/ngspice_augmentation/aaba50ba3f49-adj-asan/bin/ngspice"))

results = []


def record(name, state, detail):
    assert state in ("pass", "FAIL", "COULD-NOT-RUN", "info")
    results.append((name, state, detail))


def _strip_ac(body):
    """Remove every existing ` AC <value>` token from every line first --
    some fixtures this file reuses (pac_r3_ctrl2out.cir) are built for
    `.pac` and already carry their OWN `AC 1` on a device (Vctrl) unrelated
    to whichever ONE source a given single-source reference deck here is
    isolating; leaving it in place would silently sum two active AC
    sources (exactly the trap docs/pss-pac.md's ".pac ... solver=td"
    section documents -- "Exactly one AC-flagged device should be active
    per run"), which is what a first version of this file's helper did
    (found validating this card: an apparent ~96% .pxf-vs-.pac mismatch on
    Vin that was actually a Vin+Vctrl sum on the .pac side, not a .pxf
    engine bug -- see the card #499 report)."""
    return re.sub(r"\s+AC\s+\S+", "", body)


def _with_ac1(body, devname):
    """body (a "\\n".join(body_lines(...))-style string), with every
    existing AC flag stripped (see _strip_ac()) and then " AC 1" added to
    the first line whose FIRST token is exactly devname (case-sensitive,
    matching this project's own device-name convention) -- the same
    "exactly one AC-flagged device active per run" idiom pacshoot.c's own
    header documents for `.pac`'s netlist-AC-source forcing. Raises
    CouldNotRun if devname is not found."""
    lines = _strip_ac(body).splitlines()
    out = []
    hit = False
    for ln in lines:
        if not hit and re.match(r"^\s*" + re.escape(devname) + r"\s", ln):
            out.append(ln.rstrip("\n") + " AC 1")
            hit = True
        else:
            out.append(ln)
    if not hit:
        raise CouldNotRun(f"device '{devname}' not found in fixture body")
    return "\n".join(out)


def _pxf_deck(fixture, outdir, fname, out_node, sources, fguess, stabtime, oscnode,
              points, harmonics, sc_iter, steady_coeff, sweep, out_vecs, adjoint=None):
    """adjoint=None (default): no `adjoint=` token at all (auto-select --
    card #504's own PXF_ADJOINT_AUTO_THRESHOLD decides). adjoint=1/0:
    force the true-adjoint engine on/off explicitly (card #504)."""
    os.makedirs(outdir, exist_ok=True)
    body = "\n".join(body_lines(fixture))
    src_tokens = " ".join(f'sources="{s}"' for s in sources)
    adj_token = "" if adjoint is None else f" adjoint={int(adjoint)}"
    deck = os.path.join(outdir, fname)
    with open(deck, "w") as fh:
        fh.write(body + "\n")
        fh.write(f'.pxf {fguess} {stabtime} {oscnode} {points} {harmonics} {sc_iter} '
                  f'{steady_coeff} {out_node} {sweep} method="shoot" solver="td" '
                  f'{src_tokens}{adj_token}\n')
        fh.write(".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                  f"run\nwrdata pxf_out.txt {' '.join(out_vecs)}\nquit\n.endc\n.end\n")
    return deck


def _pac_single_source_deck(fixture, outdir, fname, devname, out_expr, fguess, stabtime,
                             oscnode, points, harmonics, sc_iter, steady_coeff, sweep):
    os.makedirs(outdir, exist_ok=True)
    body = _with_ac1("\n".join(body_lines(fixture)), devname)
    deck = os.path.join(outdir, fname)
    with open(deck, "w") as fh:
        fh.write(body + "\n")
        fh.write(f'.pac {fguess} {stabtime} {oscnode} {points} {harmonics} {sc_iter} '
                  f'{steady_coeff} {sweep} method="shoot" solver="td"\n')
        fh.write(".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                  f"run\nwrdata pac_out.txt {out_expr}\nquit\n.endc\n.end\n")
    return deck


def _ac_single_source_deck(fixture, outdir, fname, devname, out_expr, sweep):
    os.makedirs(outdir, exist_ok=True)
    body = _with_ac1("\n".join(body_lines(fixture)), devname)
    deck = os.path.join(outdir, fname)
    with open(deck, "w") as fh:
        fh.write(body + "\n")
        fh.write(f".ac {sweep}\n")
        fh.write(".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                  f"run\nwrdata ac_out.txt {out_expr}\nquit\n.endc\n.end\n")
    return deck


def _worst_rel(a, b):
    worst = 0.0
    for x, y in zip(a, b):
        rel = abs(x - y) / abs(x) if abs(x) > 0 else abs(y)
        worst = max(worst, rel)
    return worst


def check_pxf_lti():
    """.pxf vs .pac(one source) AND vs plain .ac(one source), on
    pxf_lti.cir (LTI network, two independent sources Vin/Iinj, unrelated
    periodic clock for PSS to shoot on) -- the card's own "LTI deck where
    it equals plain .ac/.tf" requirement, and one of the three .pxf-vs-.pac
    fixtures, both satisfied by this one fixture (see pxf_lti.cir's own
    header)."""
    name = "pxf(lti): pxf_lti.cir, sources=Vin,Iinj vs .pac(one source) and vs .ac (1e-6 target)"
    try:
        outdir = os.path.join(DATA_ROOT, "lti")
        sweep = "dec 10 1k 1meg"
        fguess, stabtime, oscnode = "100k", "0", "clkgnd"
        points, harmonics, sc_iter, steady_coeff = "2000", "1", "5", "0.001"

        # card #465's own accuracy note (test_pac.py's check_pac_lti()) applies
        # unchanged here: PSSpacRecordStep()'s per-step J_k reuse needs a fine
        # recording-pass grid (PAC_RECORD_POINTS, decoupled from PSS's own
        # coarse convergence Points= above) to reach 1e-6 against the EXACT
        # `.ac` reference -- .pxf reuses the SAME pac_engine_run(), so the SAME
        # override applies. Not needed for the .pxf-vs-.pac cross-check itself
        # (both sides read the SAME engine at whatever grid, so they agree to
        # machine precision regardless -- confirmed FIRST, at the default
        # points= above, before this env was added).
        fine_env = _pac.pac_record_points_env(200000)

        pxf_deck = _pxf_deck("pxf_lti.cir", outdir, "pxf.cir", "out", ["Vin", "Iinj"],
                              fguess, stabtime, oscnode, points, harmonics, sc_iter,
                              steady_coeff, sweep, ["pxf_vin", "pxf_iinj"])
        rc, out, wall = run_ngspice(pxf_deck, timeout=300, env=fine_env)
        if rc != 0:
            raise CouldNotRun(f".pxf exit {rc}: {out[-2000:]}")
        f_pxf, (h_vin, h_iinj) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        worst = 0.0
        detail_parts = []
        for src, h_pxf in (("Vin", h_vin), ("Iinj", h_iinj)):
            pac_deck = _pac_single_source_deck("pxf_lti.cir", outdir, f"pac_{src}.cir", src,
                                                "out", fguess, stabtime, oscnode, points,
                                                harmonics, sc_iter, steady_coeff, sweep)
            rc, out2, _ = run_ngspice(pac_deck, timeout=300, env=fine_env)
            if rc != 0:
                raise CouldNotRun(f".pac(only {src}) exit {rc}: {out2[-2000:]}")
            f_pac, (h_pac,) = read_wrdata_complex(os.path.join(outdir, "pac_out.txt"), 1)
            if len(f_pac) != len(f_pxf):
                raise CouldNotRun(f"{src}: point count mismatch pxf={len(f_pxf)} pac={len(f_pac)}")
            rel_pac = _worst_rel(h_pac, h_pxf)
            worst = max(worst, rel_pac)
            detail_parts.append(f"{src} vs .pac: {rel_pac:.3e}")

            ac_deck = _ac_single_source_deck("pxf_lti.cir", outdir, f"ac_{src}.cir", src,
                                              "out", sweep)
            rc, out3, _ = run_ngspice(ac_deck, timeout=60)
            if rc != 0:
                raise CouldNotRun(f".ac(only {src}) exit {rc}: {out3[-2000:]}")
            f_ac, (h_ac,) = read_wrdata_complex(os.path.join(outdir, "ac_out.txt"), 1)
            if len(f_ac) != len(f_pxf):
                raise CouldNotRun(f"{src}: point count mismatch pxf={len(f_pxf)} ac={len(f_ac)}")
            rel_ac = _worst_rel(h_ac, h_pxf)
            worst = max(worst, rel_ac)
            detail_parts.append(f"{src} vs .ac: {rel_ac:.3e}")

        detail = (f"{len(f_pxf)} points 1kHz-1MHz, worst relative error " + ", ".join(detail_parts) +
                  f" (target 1e-6 each), wall {wall:.1f}s")
        record(name, "pass" if worst <= 1e-6 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pxf_r2():
    """.pxf vs .pac(one source), on pxf_r2_switched_rc.cir (switched-RC
    LPTV network, ideal ROM `sw` -- but its control node is driven only by
    the fixed clock, never by Vin/Iinj, so saltation is a confirmed no-op
    here exactly as it is on the SAME topology's existing pac_sb0(b) check
    -- see the fixture's own header)."""
    name = "pxf(r2): pxf_r2_switched_rc.cir, sources=Vin,Iinj vs .pac(one source) (1e-6 target)"
    try:
        outdir = os.path.join(DATA_ROOT, "r2")
        sweep = "dec 10 1k 40k"
        fguess, stabtime, oscnode = "100e3", "100e-6", "vx"
        points, harmonics, sc_iter, steady_coeff = "2000", "1", "25", "1e-6"

        pxf_deck = _pxf_deck("pxf_r2_switched_rc.cir", outdir, "pxf.cir", "vx", ["Vin", "Iinj"],
                              fguess, stabtime, oscnode, points, harmonics, sc_iter,
                              steady_coeff, sweep, ["pxf_vin", "pxf_iinj"])
        rc, out, wall = run_ngspice(pxf_deck, timeout=300)
        if rc != 0:
            raise CouldNotRun(f".pxf exit {rc}: {out[-2000:]}")
        f_pxf, (h_vin, h_iinj) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        worst = 0.0
        detail_parts = []
        for src, h_pxf in (("Vin", h_vin), ("Iinj", h_iinj)):
            pac_deck = _pac_single_source_deck("pxf_r2_switched_rc.cir", outdir, f"pac_{src}.cir",
                                                src, "vx", fguess, stabtime, oscnode, points,
                                                harmonics, sc_iter, steady_coeff, sweep)
            rc, out2, _ = run_ngspice(pac_deck, timeout=300)
            if rc != 0:
                raise CouldNotRun(f".pac(only {src}) exit {rc}: {out2[-2000:]}")
            f_pac, (h_pac,) = read_wrdata_complex(os.path.join(outdir, "pac_out.txt"), 1)
            if len(f_pac) != len(f_pxf):
                raise CouldNotRun(f"{src}: point count mismatch pxf={len(f_pxf)} pac={len(f_pac)}")
            rel = _worst_rel(h_pac, h_pxf)
            worst = max(worst, rel)
            detail_parts.append(f"{src}: {rel:.3e}")

        detail = (f"{len(f_pxf)} points 1kHz-40kHz, worst relative error " + ", ".join(detail_parts) +
                  f" (target 1e-6 each), wall {wall:.1f}s")
        record(name, "pass" if worst <= 1e-6 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pxf_r3_buck():
    """.pxf vs .pac(one source), on pac_r3_ctrl2out.cir (R3 voltage-mode
    buck, smooth tanh-type switches only -- project ruling #447 A), sources
    Vin (line/supply) and Vctrl (control, already spliced into the
    compensator output for ctrl2out use) -- the card's own required
    "buck" fixture of the three .pxf-vs-.pac cross-checks."""
    name = "pxf(r3 buck): pac_r3_ctrl2out.cir, sources=Vin,Vctrl vs .pac(one source) (1e-6 target)"
    try:
        outdir = os.path.join(DATA_ROOT, "r3")
        sweep = "dec 5 500 200000"
        fguess, stabtime, oscnode = "500e3", "300e-6", "vout"
        points, harmonics, sc_iter, steady_coeff = "200", "10", "150", "1e-6"

        pxf_deck = _pxf_deck("pac_r3_ctrl2out.cir", outdir, "pxf.cir", "vout", ["Vin", "Vctrl"],
                              fguess, stabtime, oscnode, points, harmonics, sc_iter,
                              steady_coeff, sweep, ["pxf_vin", "pxf_vctrl"])
        rc, out, wall = run_ngspice(pxf_deck, timeout=300)
        if rc != 0:
            raise CouldNotRun(f".pxf exit {rc}: {out[-2000:]}")
        if "FAILED to converge" in out:
            raise CouldNotRun("PSS-SHOOT did not converge")
        f_pxf, (h_vin, h_vctrl) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        worst = 0.0
        detail_parts = []
        for src, h_pxf in (("Vin", h_vin), ("Vctrl", h_vctrl)):
            pac_deck = _pac_single_source_deck("pac_r3_ctrl2out.cir", outdir, f"pac_{src}.cir",
                                                src, "vout", fguess, stabtime, oscnode, points,
                                                harmonics, sc_iter, steady_coeff, sweep)
            rc, out2, _ = run_ngspice(pac_deck, timeout=300)
            if rc != 0:
                raise CouldNotRun(f".pac(only {src}) exit {rc}: {out2[-2000:]}")
            f_pac, (h_pac,) = read_wrdata_complex(os.path.join(outdir, "pac_out.txt"), 1)
            if len(f_pac) != len(f_pxf):
                raise CouldNotRun(f"{src}: point count mismatch pxf={len(f_pxf)} pac={len(f_pac)}")
            rel = _worst_rel(h_pac, h_pxf)
            worst = max(worst, rel)
            detail_parts.append(f"{src}: {rel:.3e}")

        detail = (f"{len(f_pxf)} points 500Hz-200kHz, worst relative error " + ", ".join(detail_parts) +
                  f" (target 1e-6 each), wall {wall:.1f}s")
        record(name, "pass" if worst <= 1e-6 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pxf_sources_default_all():
    """sources= omitted entirely defaults to "all" -- pxf_lti.cir has
    three independent V/I sources (Vclk, Vin, Iinj); a bare `.pxf` with no
    sources= token must therefore emit three vectors (pxf_vclk, pxf_vin,
    pxf_iinj), and the vin/iinj columns must match the explicit
    sources="Vin" sources="Iinj" run exactly (same engine, same
    isolation -- this is testing the PARSER's default, not new numerics)."""
    name = "pxf(default sources=all): pxf_lti.cir, bare .pxf finds all 3 sources"
    try:
        outdir = os.path.join(DATA_ROOT, "default_all")
        sweep = "dec 5 1k 1meg"
        fguess, stabtime, oscnode = "100k", "0", "clkgnd"
        points, harmonics, sc_iter, steady_coeff = "2000", "1", "5", "0.001"

        os.makedirs(outdir, exist_ok=True)
        body = "\n".join(body_lines("pxf_lti.cir"))
        deck = os.path.join(outdir, "pxf_all.cir")
        with open(deck, "w") as fh:
            fh.write(body + "\n")
            fh.write(f'.pxf {fguess} {stabtime} {oscnode} {points} {harmonics} {sc_iter} '
                      f'{steady_coeff} out {sweep} method="shoot" solver="td"\n')
            fh.write(".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                      "run\nwrdata pxf_all_out.txt pxf_vclk pxf_vin pxf_iinj\nquit\n.endc\n.end\n")
        rc, out, wall = run_ngspice(deck, timeout=300)
        if rc != 0:
            raise CouldNotRun(f".pxf (default sources) exit {rc}: {out[-2000:]}")
        f_all, (h_vclk, h_vin, h_iinj) = read_wrdata_complex(
            os.path.join(outdir, "pxf_all_out.txt"), 3)

        explicit_deck = _pxf_deck("pxf_lti.cir", outdir, "pxf_explicit.cir", "out",
                                   ["Vin", "Iinj"], fguess, stabtime, oscnode, points,
                                   harmonics, sc_iter, steady_coeff, sweep,
                                   ["pxf_vin", "pxf_iinj"])
        rc, out2, _ = run_ngspice(explicit_deck, timeout=300)
        if rc != 0:
            raise CouldNotRun(f".pxf (explicit sources) exit {rc}: {out2[-2000:]}")
        f_exp, (h_vin_e, h_iinj_e) = read_wrdata_complex(
            os.path.join(outdir, "pxf_out.txt"), 2)

        if len(f_all) != len(f_exp):
            raise CouldNotRun(f"point count mismatch default={len(f_all)} explicit={len(f_exp)}")
        worst = max(_worst_rel(h_vin, h_vin_e), _worst_rel(h_iinj, h_iinj_e))
        detail = (f"3 vectors found (pxf_vclk, pxf_vin, pxf_iinj), default-vs-explicit "
                  f"worst relative error {worst:.3e} (target 1e-6), wall {wall:.1f}s")
        record(name, "pass" if worst <= 1e-6 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_asan_pxf():
    """card #499 (report requirement): `.pxf` under AddressSanitizer, a
    multi-source, multi-point sweep (this card's own worktree ASan build,
    see ASAN_NGSPICE above) -- clean exit is the gate, same convention
    test_pac.py's check_asan_multipoint() already uses for `.pac`/`.pstb`."""
    name = "ASan: .pxf multi-source, multi-point (pxf_lti.cir, sources=Vin,Iinj)"
    if not os.path.exists(ASAN_NGSPICE):
        record(name, "COULD-NOT-RUN",
               f"{ASAN_NGSPICE} not found -- set ASAN_NGSPICE= to the fork's ASan build "
               f"(see buildrecipe/README.md's 'AddressSanitizer build' section)")
        return
    try:
        outdir = os.path.join(DATA_ROOT, "asan")
        deck = _pxf_deck("pxf_lti.cir", outdir, "asan.cir", "out", ["Vin", "Iinj"],
                          "100k", "0", "clkgnd", "2000", "1", "5", "0.001",
                          "dec 10 1k 1meg", ["pxf_vin", "pxf_iinj"])
        env = dict(os.environ)
        env["ASAN_OPTIONS"] = "detect_leaks=0"
        saved_ngspice = _pac.NGSPICE
        _pac.NGSPICE = ASAN_NGSPICE
        try:
            rc, out, wall = run_ngspice(deck, timeout=600, env=env)
        finally:
            _pac.NGSPICE = saved_ngspice
        asan_hit = "ERROR: AddressSanitizer" in out or "runtime error" in out
        if rc != 0 or asan_hit:
            record(name, "FAIL", f"exit {rc}, ASan hit={asan_hit}: {out[-3000:]}")
        else:
            record(name, "pass", f"clean exit, 31-point sweep, 2 sources, wall {wall:.1f}s")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def exact_jac_env(base_env=None):
    """card #504: the true-adjoint engine's own forward Phi/G passes force
    the exact-Jacobian convention unconditionally (see pacshoot.c's
    pac_engine_run_adjoint() design comment -- its backward sweep
    RECOMPUTES J_k at the true accepted x_k, so its own forward halves
    must use the SAME convention, not v1's default stale-x_{k-1} one).
    Comparing adjoint=1 against adjoint=0's DEFAULT output is therefore
    only accurate to the O(1/Points) staleness error this project's own
    accuracy section already documents, not floating-point noise -- every
    adjoint-vs-v1 check below runs the v1 (adjoint=0) side ALSO under
    PAC_EXACT_JACOBIAN=1 so both sides are mathematically exact and agree
    to roundoff. base_env lets pac_record_points_env()'s own dict be
    extended rather than starting fresh."""
    env = dict(base_env) if base_env is not None else dict(os.environ)
    env["PAC_EXACT_JACOBIAN"] = "1"
    return env


def check_pxf_adjoint_lti():
    """card #504: `.pxf ... adjoint=1` (true backward/transposed engine)
    vs `adjoint=0` (v1, card #499's batched-forward-per-source engine),
    BOTH under PAC_EXACT_JACOBIAN=1 (see exact_jac_env()), on pxf_lti.cir
    -- an LTI network, so this ALSO exercises "sources with no AC value"
    (neither Vin nor Iinj carries one, same as every other check in this
    file) and the trivial free-of-saltation case (no switches at all)."""
    name = "pxf-adjoint(lti): adjoint=1 vs adjoint=0 (both PAC_EXACT_JACOBIAN=1), pxf_lti.cir (1e-8 target)"
    try:
        outdir = os.path.join(DATA_ROOT, "adj_lti")
        sweep = "dec 10 1k 1meg"
        fguess, stabtime, oscnode = "100k", "0", "clkgnd"
        points, harmonics, sc_iter, steady_coeff = "2000", "1", "5", "0.001"
        env = exact_jac_env()

        adj_deck = _pxf_deck("pxf_lti.cir", outdir, "adj.cir", "out", ["Vin", "Iinj"],
                              fguess, stabtime, oscnode, points, harmonics, sc_iter,
                              steady_coeff, sweep, ["pxf_vin", "pxf_iinj"], adjoint=1)
        rc, out, wall_adj = run_ngspice(adj_deck, timeout=300, env=env)
        if rc != 0:
            raise CouldNotRun(f"adjoint=1 exit {rc}: {out[-2000:]}")
        if "adjoint done" not in out:
            raise CouldNotRun(f"adjoint=1 did not report its own 'adjoint done' line "
                               f"(engine not actually exercised?): {out[-1500:]}")
        f_adj, (h_vin_a, h_iinj_a) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        v1_deck = _pxf_deck("pxf_lti.cir", outdir, "v1.cir", "out", ["Vin", "Iinj"],
                             fguess, stabtime, oscnode, points, harmonics, sc_iter,
                             steady_coeff, sweep, ["pxf_vin", "pxf_iinj"], adjoint=0)
        rc, out2, wall_v1 = run_ngspice(v1_deck, timeout=300, env=env)
        if rc != 0:
            raise CouldNotRun(f"adjoint=0 exit {rc}: {out2[-2000:]}")
        # both decks' .control block wrdata's to the SAME hardcoded
        # "pxf_out.txt" (_pxf_deck's own convention) -- h_vin_a/h_iinj_a
        # above were already read from the adjoint run's copy before this
        # v1 run overwrites it.
        f_v1, (h_vin_v1, h_iinj_v1) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        if len(f_adj) != len(f_v1):
            raise CouldNotRun(f"point count mismatch adjoint={len(f_adj)} v1={len(f_v1)}")
        worst = max(_worst_rel(h_vin_a, h_vin_v1), _worst_rel(h_iinj_a, h_iinj_v1))
        detail = (f"{len(f_adj)} points 1kHz-1MHz, worst relative error {worst:.3e} "
                  f"(target 1e-8), wall adjoint={wall_adj:.2f}s v1={wall_v1:.2f}s")
        record(name, "pass" if worst <= 1e-8 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pxf_adjoint_r2_falls_back():
    """card #504: pxf_r2_switched_rc.cir has a real ideal-`sw` saltation
    event (Sswitch flips twice a period) -- pac_engine_run_adjoint() must
    DETECT this (PSSsaltEventsSeen()>0) and REFUSE, falling back to
    adjoint=0 transparently (see that function's own "scope cut" comment
    and docs/pss-pac.md's matching section). Checks BOTH the fallback
    message appears AND that the numbers are byte-for-byte what adjoint=0
    alone gives (trivially true by construction if the fallback actually
    ran, but checked, not assumed)."""
    name = "pxf-adjoint(r2, saltation): adjoint=1 REFUSES and falls back to v1 (exact match)"
    try:
        outdir = os.path.join(DATA_ROOT, "adj_r2")
        sweep = "dec 10 1k 40k"
        fguess, stabtime, oscnode = "100e3", "100e-6", "vx"
        points, harmonics, sc_iter, steady_coeff = "2000", "1", "25", "1e-6"

        adj_deck = _pxf_deck("pxf_r2_switched_rc.cir", outdir, "adj.cir", "vx", ["Vin", "Iinj"],
                              fguess, stabtime, oscnode, points, harmonics, sc_iter,
                              steady_coeff, sweep, ["pxf_vin", "pxf_iinj"], adjoint=1)
        rc, out, _ = run_ngspice(adj_deck, timeout=300)
        if rc != 0:
            raise CouldNotRun(f"adjoint=1 exit {rc}: {out[-2000:]}")
        if "REFUSED" not in out or "falling back" not in out:
            raise CouldNotRun(f"adjoint=1 did not report the expected saltation refusal/"
                               f"fallback message: {out[-2000:]}")
        f_adj, (h_vin_a, h_iinj_a) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        v1_deck = _pxf_deck("pxf_r2_switched_rc.cir", outdir, "v1.cir", "vx", ["Vin", "Iinj"],
                             fguess, stabtime, oscnode, points, harmonics, sc_iter,
                             steady_coeff, sweep, ["pxf_vin", "pxf_iinj"], adjoint=0)
        rc, out2, _ = run_ngspice(v1_deck, timeout=300)
        if rc != 0:
            raise CouldNotRun(f"adjoint=0 exit {rc}: {out2[-2000:]}")
        f_v1, (h_vin_v1, h_iinj_v1) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        worst = max(_worst_rel(h_vin_a, h_vin_v1), _worst_rel(h_iinj_a, h_iinj_v1))
        detail = f"fallback message present, worst relative diff vs standalone adjoint=0 run {worst:.3e}"
        record(name, "pass" if worst <= 1e-9 else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pxf_adjoint_r3_buck():
    """card #504: adjoint=1 vs adjoint=0 (both PAC_EXACT_JACOBIAN=1) on
    pac_r3_ctrl2out.cir -- smooth tanh-type switches (project ruling #447 A),
    a REAL nonlinear small-signal matrix (not just linear R/C/L like
    pxf_lti.cir), sources Vin (line) and Vctrl (control)."""
    # card #625 (s136): the default recording now runs at reltol 1e-8 on this switch-free deck, with steps down to 1.4e-13 s (2.6e-12 at the
    # inherited 1e-5), and the two engines agree to 6.4e-8 at 199 kHz (7e-10 before): a roundoff floor of the finer grid; the answer itself
    # moves 1.2e-4 between recording reltol 1e-7 and 1e-8.  Default grid: 1e-7 target; the legacy grid (PAC_RECORD_RELTOL_AUTO=0) keeps 1e-8.
    for tag, tgt, extra in (("", 1e-7, None), (", legacy grid PAC_RECORD_RELTOL_AUTO=0", 1e-8, {"PAC_RECORD_RELTOL_AUTO": "0"})):
        _check_pxf_adjoint_r3_one(tag, tgt, extra)


def _check_pxf_adjoint_r3_one(tag, tgt, extra):
    name = f"pxf-adjoint(r3 buck){tag}: adjoint=1 vs adjoint=0 (both PAC_EXACT_JACOBIAN=1) ({tgt:g} target)"
    try:
        outdir = os.path.join(DATA_ROOT, "adj_r3")
        sweep = "dec 5 500 200000"
        fguess, stabtime, oscnode = "500e3", "300e-6", "vout"
        points, harmonics, sc_iter, steady_coeff = "200", "10", "150", "1e-6"
        env = exact_jac_env()
        env.update(extra or {})

        adj_deck = _pxf_deck("pac_r3_ctrl2out.cir", outdir, "adj.cir", "vout", ["Vin", "Vctrl"],
                              fguess, stabtime, oscnode, points, harmonics, sc_iter,
                              steady_coeff, sweep, ["pxf_vin", "pxf_vctrl"], adjoint=1)
        rc, out, wall_adj = run_ngspice(adj_deck, timeout=300, env=env)
        if rc != 0:
            raise CouldNotRun(f"adjoint=1 exit {rc}: {out[-2000:]}")
        if "FAILED to converge" in out:
            raise CouldNotRun("PSS-SHOOT did not converge")
        if "adjoint done" not in out:
            raise CouldNotRun(f"adjoint=1 did not report 'adjoint done' (fell back?): {out[-2000:]}")
        f_adj, (h_vin_a, h_vctrl_a) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        os.rename(os.path.join(outdir, "pxf_out.txt"), os.path.join(outdir, "adj_out.txt"))
        v1_deck = _pxf_deck("pac_r3_ctrl2out.cir", outdir, "v1.cir", "vout", ["Vin", "Vctrl"],
                             fguess, stabtime, oscnode, points, harmonics, sc_iter,
                             steady_coeff, sweep, ["pxf_vin", "pxf_vctrl"], adjoint=0)
        rc, out2, wall_v1 = run_ngspice(v1_deck, timeout=300, env=env)
        if rc != 0:
            raise CouldNotRun(f"adjoint=0 exit {rc}: {out2[-2000:]}")
        f_v1, (h_vin_v1, h_vctrl_v1) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 2)

        if len(f_adj) != len(f_v1):
            raise CouldNotRun(f"point count mismatch adjoint={len(f_adj)} v1={len(f_v1)}")
        worst = max(_worst_rel(h_vin_a, h_vin_v1), _worst_rel(h_vctrl_a, h_vctrl_v1))
        detail = (f"{len(f_adj)} points 500Hz-200kHz, worst relative error {worst:.3e} "
                  f"(target {tgt:g}), wall adjoint={wall_adj:.2f}s v1={wall_v1:.2f}s")
        record(name, "pass" if worst <= tgt else "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


MOS_CHAIN_DECK = """* card #541: two cap-less MOS1 inverters (NO model charge -> native td path) driven by a ramp + DC Vc: Vc sets a MOS threshold-crossing time
.model nch4 nmos level=1 vto=0.5 kp=3e-4 lambda=0.05
.model pch4 pmos level=1 vto=-0.5 kp=1e-4 lambda=0.05
Vdd vdd 0 5
Vc ctl 0 DC 2.5
Vr inp ctl PULSE(-2.5 2.5 0 0.9u 0.1u 1n 2u)
Mn0 n0 inp 0 0 nch4 w=4u l=0.5u as=4p ad=4p ps=10u pd=10u
Mp0 n0 inp vdd vdd pch4 w=10u l=0.5u as=10p ad=10p ps=22u pd=22u
Mn1 out n0 0 0 nch4 w=4u l=0.5u as=4p ad=4p ps=10u pd=10u
Mp1 out n0 vdd vdd pch4 w=10u l=0.5u as=10p ad=10p ps=22u pd=22u
Rf out vo 1k
Cf vo 0 1n
.pxf 500e3 20e-6 out 400 10 50 1e-7 out lin 1 1000 1000 method="shoot" solver="td" sources="Vc" adjoint=%d
.options reltol=1e-4
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
run
wrdata pxf_out.txt pxf_vc
quit
.endc
.end
"""


def _mos_chain_pxf(outdir, tag, adjoint, env):
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, tag + ".cir")
    with open(deck, "w") as fh:
        fh.write(MOS_CHAIN_DECK % adjoint)
    rc, out, _ = run_ngspice(deck, timeout=300, env=env)
    if rc != 0:
        raise CouldNotRun(f"{tag} exit {rc}: {out[-1500:]}")
    if adjoint and "adjoint done" not in out:
        raise CouldNotRun(f"{tag}: no 'adjoint done' (fell back?): {out[-800:]}")
    _, (h,) = read_wrdata_complex(os.path.join(outdir, "pxf_out.txt"), 1)
    return h


def check_pxf_adjoint_mos_gate():
    """card #541 (#539 bug 3): adjoint=1 on a cap-less MOS1 chain (native td path) where the source sets a MOS
    threshold crossing. Before: adjoint ~0.02 vs forward 0.53 (the backward sweep reloaded J_k in reverse order
    from device state; the forward 'exact' J_k was also loaded in MODEINITPRED, i.e. at a PREDICTED point).
    (a) adjoint vs forward, both PAC_EXACT_JACOBIAN=1 (1e-6); (b) default forward vs exact forward (1e-3: the exact
    load must be at x_k); (c) magnitude sanity (static d<v>/dVc ~ 0.5, so > 0.3)."""
    name = "pxf-adjoint(mos gate, native td): adjoint==forward (1e-6), exact==default J (1e-3), |H|>0.3"
    try:
        outdir = os.path.join(DATA_ROOT, "adj_mosgate")
        env_x = exact_jac_env()
        h_adj = _mos_chain_pxf(outdir, "adj", 1, env_x)
        h_fx = _mos_chain_pxf(outdir, "fwd_exact", 0, env_x)
        h_f0 = _mos_chain_pxf(outdir, "fwd_default", 0, dict(os.environ))
        e_a = _worst_rel(h_adj, h_fx)
        e_j = _worst_rel(h_fx, h_f0)
        mag = min(abs(complex(z)) for z in h_fx)
        ok = e_a <= 1e-6 and e_j <= 1e-3 and mag > 0.3
        record(name, "pass" if ok else "FAIL",
               f"adjoint vs forward {e_a:.3e} (1e-6), exact vs default-J forward {e_j:.3e} (1e-3), "
               f"min |H_fwd| {mag:.4f}, min |H_adj| {min(abs(complex(z)) for z in h_adj):.4f}")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pxf_adjoint_timing():
    """card #504 done-when: "wall time for N sources is ~1 sweep: measure
    N=2, 10, 50 ... vs v1" -- tests/fixtures/pxf_manysrc.cir (Vin + 50
    independent current injectors I1..I50, all LTI, no switches). Forces
    BOTH engines explicitly (adjoint=1 / adjoint=0) for a clean timing
    comparison (no auto-selection ambiguity) at N in {2,10,50}. This is a
    wall-clock MEASUREMENT (recorded as "info"), plus one real pass/FAIL
    gate: adjoint's own N=50 wall time must not exceed its N=2 wall time
    by more than a generous margin (5x) -- the qualitative "~flat in N"
    claim, checked, not just asserted -- while v1's is expected (and
    allowed) to grow close to linearly."""
    name_prefix = "pxf-adjoint timing"
    try:
        outdir = os.path.join(DATA_ROOT, "adj_timing")
        sweep = "dec 5 1k 100k"
        fguess, stabtime, oscnode = "100k", "0", "clkgnd"
        points, harmonics, sc_iter, steady_coeff = "2000", "1", "5", "0.001"

        def sources_for(nsrc):
            srcs = ["Vin"]
            for k in range(1, nsrc):
                srcs.append(f"I{k}")
            return srcs

        table = []
        for nsrc in (2, 10, 50):
            srcs = sources_for(nsrc)
            out_vecs = [f"pxf_{s.lower()}" for s in srcs]
            row = {"n": nsrc}
            for label, adj in (("adjoint", 1), ("v1", 0)):
                deck = _pxf_deck("pxf_manysrc.cir", outdir, f"{label}_{nsrc}.cir", "out", srcs,
                                  fguess, stabtime, oscnode, points, harmonics, sc_iter,
                                  steady_coeff, sweep, out_vecs, adjoint=adj)
                rc, out, wall = run_ngspice(deck, timeout=600)
                if rc != 0:
                    raise CouldNotRun(f"{label} N={nsrc} exit {rc}: {out[-2000:]}")
                if label == "adjoint" and "adjoint done" not in out:
                    raise CouldNotRun(f"adjoint N={nsrc} fell back unexpectedly: {out[-2000:]}")
                row[label] = wall
            table.append(row)
            record(f"{name_prefix}: N={nsrc}", "info",
                   f"adjoint={row['adjoint']:.3f}s v1={row['v1']:.3f}s "
                   f"(v1/adjoint={row['v1']/row['adjoint']:.2f}x)")

        adj2 = table[0]["adjoint"]
        adj50 = table[2]["adjoint"]
        v1_2 = table[0]["v1"]
        v1_50 = table[2]["v1"]
        ratio_adj = adj50 / adj2 if adj2 > 0 else float("inf")
        ratio_v1 = v1_50 / v1_2 if v1_2 > 0 else float("inf")
        detail = (f"table: " + "; ".join(f"N={r['n']} adjoint={r['adjoint']:.3f}s v1={r['v1']:.3f}s"
                                          for r in table) +
                  f" -- adjoint(N=50)/adjoint(N=2)={ratio_adj:.2f}x, "
                  f"v1(N=50)/v1(N=2)={ratio_v1:.2f}x (target: adjoint stays roughly flat, "
                  f"<=5x; v1 is allowed and expected to scale up close to linearly)")
        record(f"{name_prefix}: adjoint stays ~flat in N (<=5x from N=2 to N=50)",
               "pass" if ratio_adj <= 5.0 else "FAIL", detail)
    except CouldNotRun as e:
        record(f"{name_prefix}: table", "COULD-NOT-RUN", str(e))


def check_asan_pxf_adjoint():
    """card #504 done-when: ASan clean on the adjoint engine itself (not
    just v1) -- a multi-source, multi-point sweep with adjoint=1 forced,
    on this card's own ASan build (ADJ_ASAN_NGSPICE, buildrecipe/README.md
    recipe, --prefix-name=aaba50ba3f49-adj-asan on ngspice-src-adj)."""
    name = "ASan: .pxf adjoint=1 multi-source (pxf_manysrc.cir, N=10)"
    if not os.path.exists(ADJ_ASAN_NGSPICE):
        record(name, "COULD-NOT-RUN",
               f"{ADJ_ASAN_NGSPICE} not found -- set ADJ_ASAN_NGSPICE= to this card's ASan "
               f"build (see buildrecipe/README.md's 'AddressSanitizer build' section)")
        return
    try:
        outdir = os.path.join(DATA_ROOT, "adj_asan")
        srcs = ["Vin"] + [f"I{k}" for k in range(1, 10)]
        out_vecs = [f"pxf_{s.lower()}" for s in srcs]
        deck = _pxf_deck("pxf_manysrc.cir", outdir, "asan.cir", "out", srcs,
                          "100k", "0", "clkgnd", "2000", "1", "5", "0.001",
                          "dec 5 1k 100k", out_vecs, adjoint=1)
        env = dict(os.environ)
        env["ASAN_OPTIONS"] = "detect_leaks=0"
        saved_ngspice = _pac.NGSPICE
        _pac.NGSPICE = ADJ_ASAN_NGSPICE
        try:
            rc, out, wall = run_ngspice(deck, timeout=600, env=env)
        finally:
            _pac.NGSPICE = saved_ngspice
        asan_hit = "ERROR: AddressSanitizer" in out or "runtime error" in out
        fell_back = "adjoint done" not in out
        if rc != 0 or asan_hit:
            record(name, "FAIL", f"exit {rc}, ASan hit={asan_hit}: {out[-3000:]}")
        elif fell_back:
            record(name, "FAIL", f"adjoint engine fell back unexpectedly: {out[-2000:]}")
        else:
            record(name, "pass", f"clean exit, 10 sources, wall {wall:.1f}s")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_asan_harmonics1():
    """card #502: the pin's frequency-domain PSS engine (dcpss.c, used by plain
    `.pss` and by `.pxf` without method="shoot") read pssResults[msize] and
    pssfreqs[1] with harmonics=1, one past a CKTharms-entry buffer (heap
    overflow at dcpss.c:4477 under ASan, found reviewing #499). Red-proof: the
    pre-fix s47 ASan build reports heap-buffer-overflow on both decks."""
    for kind in ("pss", "pxf"):
        name = f"ASan: frequency-domain .{kind} with harmonics=1 (card #502)"
        if not os.path.exists(ASAN_NGSPICE):
            record(name, "COULD-NOT-RUN", f"{ASAN_NGSPICE} not found")
            continue
        try:
            outdir = os.path.join(DATA_ROOT, "asan_h1")
            os.makedirs(outdir, exist_ok=True)
            body = "\n".join(body_lines("pxf_lti.cir"))
            line = (".pss 100k 0 clkgnd 2000 1 5 0.001" if kind == "pss" else
                    ".pxf 100k 0 clkgnd 2000 1 5 0.001 out dec 10 1k 1meg")
            deck = os.path.join(outdir, f"h1_{kind}.cir")
            with open(deck, "w") as fh:
                fh.write(body + "\n" + line + "\n.control\nrun\nquit\n.endc\n.end\n")
            env = dict(os.environ)
            env["ASAN_OPTIONS"] = "detect_leaks=0"
            saved_ngspice = _pac.NGSPICE
            _pac.NGSPICE = ASAN_NGSPICE
            try:
                rc, out, wall = run_ngspice(deck, timeout=600, env=env)
            finally:
                _pac.NGSPICE = saved_ngspice
            asan_hit = "ERROR: AddressSanitizer" in out or "runtime error" in out
            if rc != 0 or asan_hit:
                record(name, "FAIL", f"exit {rc}, ASan hit={asan_hit}: {out[-2000:]}")
            else:
                record(name, "pass", f"clean exit, wall {wall:.1f}s")
        except CouldNotRun as e:
            record(name, "COULD-NOT-RUN", str(e))


def _sw_ctrl_direct_pxf_point(f, env_extra=None):
    """One native `.pxf sources="all"` frequency point on the committed
    tests/fixtures/pac_sw_ctrl_direct.cir (card #510) -- strips the
    fixture's OWN `.pac`/`.control` block first (see test_pac.py's
    `_ctrl_direct_circuit_only_lines()` docstring for why: reusing it
    verbatim, the way `.pac` fixtures with no block of their own do
    elsewhere in this file, silently re-runs the FIXTURE's own .pac point
    instead of this function's requested one -- found while writing this
    check). Returns (pxf_vctrl, pxf_vramp, pxf_vin) as complex."""
    tag = f"{f:g}"
    outdir = os.path.join(DATA_ROOT, "sw_ctrl_direct", tag)
    os.makedirs(outdir, exist_ok=True)
    body = "\n".join(_pac._ctrl_direct_circuit_only_lines())
    deck = os.path.join(outdir, "pxf.cir")
    with open(deck, "w") as fh:
        fh.write(body + "\n")
        fh.write(f'.pxf 500e3 0 vout 200 1 25 1e-6 vout lin 1 {f!r} {f!r} '
                  f'method="shoot" solver="td"\n')
        fh.write(".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                  "run\nwrdata pxf_out.txt pxf_vctrl pxf_vramp pxf_vin\nquit\n.endc\n.end\n")
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    rc, out, wall = run_ngspice(deck, timeout=120, env=env)
    if rc != 0:
        raise CouldNotRun(f".pxf exit {rc}: {out[-1500:]}")
    if "FAILED to converge" in out:
        raise CouldNotRun("PSS-SHOOT did not converge")
    freqs, (h_vctrl, h_vramp, h_vin) = read_wrdata_complex(
        os.path.join(outdir, "pxf_out.txt"), 3)
    return h_vctrl[0], h_vramp[0], h_vin[0]


def check_pxf_sw_ctrl_direct_nonzero():
    """card #510: the same silently-wrong-zero bug test_pac.py's
    check_pac_sw_ctrl_direct_saltation()/_redproof() cover for `.pac`
    also hit `.pxf sources="all"` on the identical fixture -- pxf_vctrl
    and pxf_vramp were EXACTLY 0+0j at every frequency (the card's own
    bug report), while pxf_vin was already sensible (Vin does not drive
    the switch's control port, so it was never on the affected path).
    `.pxf`'s v1 (non-adjoint) engine is pacshoot.c's SAME
    PSSpacRecordStep()/g_pac_Cnew machinery `.pac` uses (card #504's true-
    adjoint engine REFUSES a saltation-affected deck and falls back to
    v1 -- see pxf-adjoint(r2, saltation) above -- so this fixture, which
    DOES have saltation events, always takes that path), so the fix is
    the SAME pssshoot.c change; this check confirms it end-to-end through
    `.pxf`'s own command/vector-naming path rather than assuming.
    Checks (1kHz): pxf_vctrl and pxf_vramp both nonzero; pxf_vctrl ==
    -pxf_vramp (h=v_ctrl-v_ramp -- a unit perturbation on either terminal
    of the SAME differential control must give exactly opposite
    sensitivities, a free internal-consistency check, not just "nonzero");
    pxf_vctrl matches test_pac.py's own native `.pac` vout gain at 1kHz
    (same engine, Vctrl AC=1 in `.pac` IS exactly pxf_vctrl's own
    decomposition) to a tight numerical tolerance; red-proof,
    PSS_SALTATION=0 -> pxf_vctrl exactly 0."""
    name = "pxf(sw_ctrl_direct, card #510): pxf_vctrl/pxf_vramp nonzero, mirror-symmetric, matches .pac"
    try:
        h_vctrl, h_vramp, h_vin = _sw_ctrl_direct_pxf_point(1000)
        if h_vctrl == 0j or h_vramp == 0j:
            record(name, "FAIL",
                   f"pxf_vctrl={h_vctrl} pxf_vramp={h_vramp} -- expected both nonzero "
                   "(the original card #510 bug is back)")
            return
        mirror_rel = abs(h_vctrl + h_vramp) / abs(h_vctrl)
        h_pac = _pac._ctrl_direct_native_pac_point(1000)
        pac_rel = abs(h_vctrl - h_pac) / abs(h_pac) if abs(h_pac) > 0 else abs(h_vctrl)
        detail = (f"pxf_vctrl={h_vctrl:.6g} pxf_vramp={h_vramp:.6g} pxf_vin={h_vin:.6g} "
                  f"mirror |vctrl+vramp|/|vctrl|={mirror_rel:.3e} (target 1e-6) "
                  f".pac vout={h_pac:.6g} vs pxf_vctrl relerr={pac_rel:.3e} (target 1e-6)")
        if mirror_rel <= 1e-6 and pac_rel <= 1e-6:
            record(name, "pass", detail)
        else:
            record(name, "FAIL", detail)
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
    except Exception as e:  # noqa: BLE001
        record(name, "COULD-NOT-RUN", f"{type(e).__name__}: {e}")


def check_pxf_sw_ctrl_direct_redproof():
    """card #510 red-proof: PSS_SALTATION=0 -> pxf_vctrl exactly 0+0j at
    1kHz on pac_sw_ctrl_direct.cir (mirrors test_pac.py's
    check_pac_sw_ctrl_direct_redproof(), through `.pxf` instead of
    `.pac`)."""
    name = "pxf(sw_ctrl_direct) red-proof: PSS_SALTATION=0 -> pxf_vctrl exactly 0"
    try:
        h_vctrl, h_vramp, h_vin = _sw_ctrl_direct_pxf_point(1000, env_extra={"PSS_SALTATION": "0"})
        detail = f"PSS_SALTATION=0: pxf_vctrl={h_vctrl} pxf_vramp={h_vramp}"
        if h_vctrl == 0j and h_vramp == 0j:
            record(name, "pass", detail)
        else:
            record(name, "FAIL", detail + " (expected both exactly 0)")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))
    except Exception as e:  # noqa: BLE001
        record(name, "COULD-NOT-RUN", f"{type(e).__name__}: {e}")


# --- card #557: `.pxf adjoint=1` on FREE-RUNNING (trigger=) orbits ----------------------------------------------------------
R5_TRIG = 'method="shoot" solver="td" trigger="q" trigval=0.5 trigdir="rise"'
R8_TRIG = 'method="shoot" solver="td" trigger="cx" trigval=1.7055 trigdir="rise"'


def _plain_body(fixture):
    """fixture body without analysis cards / control block / AC flags (r8_pac_*.cir carry them)"""
    out = []
    for l in body_lines(fixture):
        t = l.strip().lower()
        if t.startswith((".pac", ".control", ".endc", "set ", "wrdata")) or t == "run":
            continue
        out.append(l)
    return _strip_ac("\n".join(out))


def _fr_deck(outdir, fname, body, card, outvecs, outtxt):
    os.makedirs(outdir, exist_ok=True)
    deck = os.path.join(outdir, fname)
    with open(deck, "w") as fh:
        fh.write(body + "\n" + card + "\n")
        fh.write(".control\nset wr_singlescale\nset wr_vecnames\nset numdgt=15\n"
                  f"run\nwrdata {outtxt} {' '.join(outvecs)}\nquit\n.endc\n.end\n")
    return deck


def _freerun_adjoint_case(tag, body, head, trig, sweep, env_fwd=None, env_adj=None, timeout=900):
    """returns (adjoint=1 curves, adjoint=0 curves, logs) for the two sources Vin / Vrefsrc"""
    outdir = os.path.join(DATA_ROOT, "adj_" + tag)
    res = {}
    for adj in (1, 0):
        txt = os.path.join(outdir, f"pxf{adj}.txt")
        if os.path.exists(txt):
            os.remove(txt)
        deck = _fr_deck(outdir, f"pxf{adj}.cir", body,
                        f'.pxf {head} vout {sweep} sources="Vin" sources="Vrefsrc" adjoint={adj} {trig}',
                        ["pxf_Vin", "pxf_Vrefsrc"], txt)
        rc, out, _ = run_ngspice(deck, timeout=timeout, env=(env_adj if adj else env_fwd))
        if rc != 0:
            raise CouldNotRun(f"{tag} adjoint={adj} exit {rc}: {out[-1500:]}")
        f, (h1, h2) = read_wrdata_complex(txt, 2)
        res[adj] = ((h1, h2), out)
    return res[1][0], res[0][0], (res[1][1], res[0][1])


def check_pxf_adjoint_freerun_r5():
    """card #557: `.pxf adjoint=1` on the free-running R5 (hysteretic buck, tests/fixtures/r5_freerun_plain.cir, trigger="q") equals the forward
    per-source engine (adjoint=0, PAC_EXACT_JACOBIAN=1) at 20/55/90 kHz for Vin and the hysteresis setpoint Vrefsrc.  The forward free-running
    operator is the fixed-period LPTV engine on the polished orbit (the period-variation projection is opt-in, PAC_FREERUN_PROJ), whose exact transpose
    the adjoint sweep is.  Red proof: the forward engine WITH the projection (PAC_FREERUN_PROJ=1) is a different operator and must differ; and with the
    projection requested the adjoint refuses (and falls back) instead of transposing the wrong operator."""
    name = "pxf-adjoint(R5 free-running, trigger=q): adjoint=1 == forward per-source to 1e-6 (and the opt-in projection still refuses)"
    try:
        body = "\n".join(body_lines("r5_freerun_plain.cir"))
        head, sweep = "427350 450e-6 vout 200 10 45 1e-4", "lin 3 20000 90000"
        env_x = exact_jac_env()
        a, f, (la, lf) = _freerun_adjoint_case("r5", body, head, R5_TRIG, sweep, env_fwd=env_x, env_adj=env_x)
        if "adjoint done" not in la:
            raise CouldNotRun("the adjoint engine did not run on the free-running orbit: " + la[-800:])
        w = max(_worst_rel(a[0], f[0]), _worst_rel(a[1], f[1]))
        env_p = exact_jac_env()
        env_p["PAC_FREERUN_PROJ"] = "1"
        ap, fp, (lap, lfp) = _freerun_adjoint_case("r5proj", body, head, R5_TRIG, sweep, env_fwd=env_p, env_adj=env_p)
        refused = "adjoint engine REFUSED" in lap and "PAC_FREERUN_PROJ" in lap and "adjoint done" not in lap
        red = _worst_rel(f[0], fp[0])
        ok = w <= 1e-6 and refused and red > 1e-3
        record(name, "pass" if ok else "FAIL",
               f"worst relative difference adjoint vs forward {w:.2e} (bound 1e-6); with PAC_FREERUN_PROJ=1 the forward answer moves by {red:.2e} "
               f"(red proof, must exceed 1e-3) and the adjoint refuses+falls back: {refused}")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def check_pxf_adjoint_freerun_r8():
    """card #557: same on the free-running transistor-level R8 (BSIM4 hysteretic buck, trigger=cx, generic charge path, hidden-state continuation)."""
    name = "pxf-adjoint(R8 free-running BSIM4, trigger=cx): adjoint=1 == forward per-source to 1e-6"
    try:
        body = _plain_body("r8_pac_ctrl2out.cir")
        a, f, (la, lf) = _freerun_adjoint_case("r8", body, "395e3 100u q 400 10 150 1e-6", R8_TRIG, "lin 3 20000 80000", timeout=1500)
        if "adjoint done" not in la:
            raise CouldNotRun("the adjoint engine did not run on the free-running orbit: " + la[-800:])
        w = max(_worst_rel(a[0], f[0]), _worst_rel(a[1], f[1]))
        record(name, "pass" if w <= 1e-6 else "FAIL", f"worst relative difference adjoint vs forward {w:.2e} (bound 1e-6; measured 3e-14), 20/50/80 kHz, Vin and Vrefsrc")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))



def check_asan_pxf_freerun():
    """card #557: ASan clean on `.pxf adjoint=1` through a FREE-RUNNING orbit (R5, trigger q, two sources, three frequencies)."""
    name = "ASan: .pxf adjoint=1 on the free-running R5 (2 sources, 3 frequencies)"
    if not os.path.exists(ADJ_ASAN_NGSPICE):
        record(name, "COULD-NOT-RUN", f"{ADJ_ASAN_NGSPICE} not found -- set NGSPICE_ASAN=")
        return
    try:
        outdir = os.path.join(DATA_ROOT, "adj_r5asan")
        txt = os.path.join(outdir, "pxf1.txt")
        deck = _fr_deck(outdir, "asan.cir", "\n".join(body_lines("r5_freerun_plain.cir")),
                        f'.pxf 427350 450e-6 vout 200 10 45 1e-4 vout lin 3 20000 90000 sources="Vin" sources="Vrefsrc" adjoint=1 {R5_TRIG}',
                        ["pxf_Vin", "pxf_Vrefsrc"], txt)
        env = dict(os.environ)
        env["ASAN_OPTIONS"] = "detect_leaks=0"
        saved = _pac.NGSPICE
        _pac.NGSPICE = ADJ_ASAN_NGSPICE
        try:
            rc, out, wall = run_ngspice(deck, timeout=900, env=env)
        finally:
            _pac.NGSPICE = saved
        hit = "ERROR: AddressSanitizer" in out or "runtime error" in out
        record(name, "pass" if rc == 0 and not hit and "adjoint done" in out else "FAIL",
               f"exit {rc}, ASan hit {hit}, adjoint ran {'adjoint done' in out}, wall {wall:.1f}s")
    except CouldNotRun as e:
        record(name, "COULD-NOT-RUN", str(e))


def main():
    check_pxf_lti()
    check_pxf_r2()
    check_pxf_r3_buck()
    check_pxf_sources_default_all()
    check_asan_pxf()
    check_asan_harmonics1()

    # card #510: silently-wrong-zero saltation forcing-term bug, .pxf side.
    check_pxf_sw_ctrl_direct_nonzero()
    check_pxf_sw_ctrl_direct_redproof()

    # card #504: true-adjoint (`adjoint=1`) engine.
    check_pxf_adjoint_lti()
    check_pxf_adjoint_r2_falls_back()
    check_pxf_adjoint_r3_buck()
    check_pxf_adjoint_mos_gate()
    check_pxf_adjoint_timing()
    check_asan_pxf_adjoint()

    # card #557: free-running orbits.
    check_pxf_adjoint_freerun_r5()
    check_pxf_adjoint_freerun_r8()
    check_asan_pxf_freerun()

    width = max(len(n) for n, _, _ in results)
    for n, state, detail in results:
        print(f"{state:15s} {n:{width}s}  {detail}")
    states = {s for _, s, _ in results}
    npass = sum(s == "pass" for _, s, _ in results)
    ninfo = sum(s == "info" for _, s, _ in results)
    print(f"\n{npass}/{len(results) - ninfo} pass ({ninfo} info, not counted)")
    return 1 if "FAIL" in states else (2 if "COULD-NOT-RUN" in states else 0)


if __name__ == "__main__":
    sys.exit(main())
