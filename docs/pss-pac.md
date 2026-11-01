# Periodic steady state, periodic AC and periodic loop gain (`.pss method=shoot`, `.pac solver=td`, `.pstb`) in ngspice_augmentation

Written for analog designers who know Cadence Spectre's `pss`/`pac`/`pstb`
but not this project's code. Card #465 (plan T5.1-T5.3), private fork only
(project ruling #446 B, 2026-09-27 -- **never** contact upstream ngspice about
any of this). All example output below was captured on this card's FIFTH
pass build, prefix `aaba50ba3f49-pac5`
(`~/Work/apps/ngspice_augmentation/aaba50ba3f49-pac5/bin/ngspice-aug`, branch
`aug-pac` = pin + `patches/0001-0014` + `patches-wip/pss/0015-0020` +
`patches-wip/pac/0021-0026`), 2026-09-29.

Companion tools in this project: `tools/stb.py`/native `.stb` (loop gain of
a circuit already at a **DC** operating point -- see `docs/stb.md`) and
`tools/fra.py` (frequency response by finite-difference perturbation of a
**switching** circuit's periodic orbit, no custom build needed). Everything
in this document is instead a **native, in-simulator** analysis on this
project's private fork: no finite differencing, no perturb-and-subtract --
the engine linearises the periodic orbit itself (LPTV: linear
periodically-time-varying) and solves for the small-signal response
directly.

## What these three do

- **`.pss method=shoot`**: finds the periodic steady state (one switching
  period) of a **driven** (clocked) switching converter by shooting-Newton
  iteration on the one-period map, instead of the existing (default, still
  unmodified) `dcpss.c` algorithm. This is the circuit every `.pac`/`.pstb`
  card below needs first -- either run explicitly, or implicitly (a `.pac`/
  `.pstb` card on its own also runs the PSS first).
- **`.pac ... solver=td`**: sideband-0 periodic AC. Given the converged
  periodic orbit, linearises the circuit's own time-varying G(t)/C(t) around
  it and computes the small-signal transfer function, at each requested
  frequency, from **the netlist's own AC-flagged source(s)** (same `AC 1`
  convention `.ac` uses) to **every node/branch** in the circuit -- read
  back with an ordinary `wrdata`/`print`. This is a control-to-output,
  line-to-output or output-impedance measurement of a switching converter,
  done natively (no perturb-and-subtract, no injected-amplitude choice to
  get wrong).
- **`.pstb`**: sideband-0 **periodic loop gain**, using Tian's two-injection
  method (same convention `.stb`/`tools/stb.py` use: `loopGain = -T`) at a
  plain top-level probe V-source spliced into the loop. This is the
  switching-converter analogue of `.stb` -- Spectre's `pstb`.

All three share one engine (`src/spicelib/analysis/pssshoot.c` for the PSS,
`pacshoot.c` for the PAC/PSTB post-processing) and one restriction, **project
ruling #447 A (2026-09-27): smooth switch models only.** Every device in the
circuit must be a continuous function of its inputs (behavioural `tanh`-type
switches/comparators, as every reference deck in `refs/` uses) -- no ideal
SPICE `sw`, no XSPICE digital. Saltation correction (the extra term an
event-driven switch's small-signal model needs at the instant it trips) is
out of scope and not implemented; running these cards against an ideal-`sw`
deck will not be refused, but the result is not trustworthy at the switching
instants.

**Free-running (autonomous) oscillators are not supported** (project ruling
#446 B/#447 A) -- `method=shoot` only handles a **driven** circuit (a fixed,
known clock), the same restriction the PSS engine's caller already enforces
before this code is ever entered. `fguess` here is the deck's own known
clock frequency, not a guess to be refined.

## `.pss method=shoot`

```
.pss Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff [uic] [method=shoot]
```

or, in a `.control` block: `pss Fguess StabTime OscNode Points Harmonics
SC_iter Steady_coeff [uic] method=shoot`.

| field | meaning |
|---|---|
| `Fguess` | the clock frequency (Hz) -- exact, this is a driven circuit, not a guess |
| `StabTime` | the pre-stabilisation budget (seconds) -- see "Choosing StabTime" below. `0` skips straight to Newton from the deck's own `.ic` (only safe on a circuit with no degenerate all-off fixed point, e.g. `pac_lti.cir`/`r2_switched_rc.cir`) |
| `OscNode` | a node name -- for a driven PSS this is only used for reporting, not for period detection (the period is `1/Fguess`, fixed) |
| `Points` | integration resolution within one period (`CKTmaxStep = T/Points`) -- also the default recording-pass resolution for `.pac`/`.pstb` unless `pacpoints=` overrides it (see below) |
| `Harmonics` | accepted for compatibility with the existing (unmodified) `.pss` card grammar; **unused** by `method=shoot` (no Fourier-harmonic machinery in this path) |
| `SC_iter` | maximum Newton iterations for the **final** stage (probe stages get a smaller, fixed budget -- see "Choosing StabTime") |
| `Steady_coeff` | Newton convergence tolerance, **ripple-relative**: converged when every Newton-unknown state's own residual is below `Steady_coeff * max(that state's own switching ripple, a floor)` -- not an absolute or magnitude-relative test, so it does not falsely "converge" on a state that is barely moving relative to its own DC value |
| `uic` | flag, no value -- always effectively on for this engine (every period run restarts from IC) |
| `method=shoot` | selects this engine; any other value (or unset) keeps the existing, unmodified `dcpss.c` path |

On success, stdout carries the summary line and a Floquet report:

```
PSS-SHOOT: CONVERGED iter=4 residual=1.218278e-07 spectral_radius=0.993978 stable wall_time=1.346s
PSS-SHOOT: FLOQUET dominant_multiplier=0.993978+0.000000i angle_deg=0.000
```

`spectral_radius` is the modulus of the dominant Floquet multiplier of the
one-period map (restricted to the Newton-unknown/"slow" states -- see fast-
state reduction below): **`< 1` is stable, `>= 1` is not a real steady
state** (the engine reports non-convergence instead of a false "converged"
answer in that case -- see the red-proof decks in `tests/fixtures/` and
"Interpreting non-convergence" below). A **negative real** dominant
multiplier near `-1` is the textbook signature of **period-doubling**
(subharmonic) instability -- the peak-current-mode-without-slope-compensation
case at D>0.5 (`refs/r4_pcm_buck.cir`, `vin_v=6.0`).

On the converged orbit, the trajectory is also emitted as an ordinary
time-domain plot ("Time Domain Periodic Steady State Analysis"), same as the
existing algorithm -- `wrdata`/`plot`/`print` work exactly as after any
other analysis.

### Choosing StabTime and Points

A raw `.ic`/IC=0 start is very often itself a degenerate fixed point (a
converter sitting powered-off, comparator/compensator settled at whatever a
cold `op` finds, is self-consistent and non-switching) -- Newton would
"converge" to that in a couple of iterations without ever seeing the real
orbit. `StabTime` exists to escape that trap: a **staged warm-up**, trying a
short pre-stabilisation transient first (tens of periods) and only paying
for a longer one if a fast Newton probe on top of it still fails, capped by
`StabTime` itself in the worst case -- this makes `StabTime` a **budget**,
not a value the deck author needs to tune precisely: pass whatever a plain
`.tran` needs to fully settle (the same number `refs/*.cir`'s own header
comments already give), and the shooter will usually converge well before
spending all of it.

**How long a real converter needs**: the dominant Floquet multiplier
(reported above) tells you directly -- a slow mode decays by `multiplier^N`
over `N` periods, so a multiplier of `0.98` needs roughly `N ~
ln(tol)/ln(0.98)` periods to settle to a given tolerance. **The closer the
multiplier is to the unit circle (stable side), the longer the warm-up (and
the more Newton iterations) actually needed** -- this project's own R4
(peak-current-mode buck) at D~=0.45 has spectral radius 0.994 (a genuinely
slow, but stable, mode -- see "R4 D~=0.45: a near-unity-but-stable Floquet
mode" below for the full story, including a real engine bug this pass found
and fixed there), needing `StabTime=3500e-6` (1750 periods) where the same
circuit at D~=0.30 (spectral radius 0.994 too, coincidentally similar,
but converges faster in practice) is comfortable at `1500e-6`.

`Points` sets the integration grid within one period; it is a genuine
accuracy/speed trade-off for the FD monodromy and for `.pac`/`.pstb`'s own
recording-pass demodulation (see `pacpoints=` below to decouple the two).
200 is this project's own default for R3/R4 (7/14-state converters); 2000+
for a simple 1-2 state circuit costs little and buys more digits (see
`pac_sb0(a)` below, Points=200000 reaching 8.2e-9 relative error against
`.ac` on a trivial RC).

### Interpreting non-convergence

A circuit that is genuinely unstable, or genuinely period-doubling, will
**correctly fail to converge** -- this is not a bug to work around. Two
worked examples from this project's own regression suite:

- R4 at D~=0.55 with no slope compensation (`vin_v=6.0`,
  `tests/fixtures/pss_r4_d55_unstable.cir`): Newton residual plateaus around
  `~1.3e-3` -- consistent with **oscillating around, not converging to**, a
  period-1 fixed point that does not exist there (the real orbit is
  period-2). `.tran` independently confirms `t_period` measuring almost
  exactly `2x` the nominal clock period at this operating point.
- R3 with the feedback gain multiplied 20x (`pss_r3_gainx20_unstable.cir`):
  a genuinely unstable loop. An earlier version of this engine used to
  report a spurious `CONVERGED specrad~1.000001` here -- traced to Newton
  finding the **degenerate all-off trivial fixed point**, not the true
  marginal orbit; a guard added since correctly turns that into an honest
  non-convergence instead.

`PSS_SHOOT_DEBUG=1` in the environment prints the residual at every Newton
iteration (including every warm-up stage) and, on convergence, the Floquet
eigenvalues themselves -- useful for telling "genuinely not a steady state"
apart from "needs more StabTime/Points" (see the case study below for
exactly this diagnosis).

### R4 D~=0.45: a near-unity-but-stable Floquet mode, and the fix (card #465 FIFTH pass)

**Symptom** (found on `refs/r4_pcm_buck.cir`, `vin_v=7.33`, D~=0.45): PSS
shooting reliably converged at D~=0.30 (this same deck's default) but
**failed to converge at D~=0.45**, even though the fork's own FRA loop-gain
measurement showed a perfectly stable, well-damped orbit there (PM=89.37deg)
and `StabTime` was already generous (3500us = 1750 periods). The reported
final residual plateaued around `3e-4` to `6e-5` depending on exactly which
fixture/parameters were used -- close to, but never reaching, the `1e-6`
target, and the run reported `FAILED to converge` rather than a wrong
answer.

**Diagnosis** (`PSS_SHOOT_DEBUG=1`, reproduced on a local repro deck): the
dominant Floquet multiplier at D~=0.45 is **0.994** -- large (slow-decaying)
but genuinely `< 1` (stable), consistent with the FRA loop gain's own high
phase margin there. This makes `(Phi - I)` (the matrix Newton's linear
solve uses) have an eigenvalue of only about `0.006` along that direction,
amplifying whatever noise is in the Newton step roughly `166x` along it.
Tracing the actual failure with `PSS_SHOOT_DEBUG=1`: on the final Newton
stage, a freshly-recomputed finite-difference Jacobian produced a direction
that, even after 24 backtracking halvings down to `alpha~1e-7` (a step size
that should leave the trial point numerically indistinguishable from the
current one), still reported a **residual 3-5x worse than the current
baseline** -- and did not shrink further as `alpha` continued to shrink,
which is not possible for a genuine small perturbation. Root cause: the
damped-step **trial** period run in `pssshoot.c`'s Newton loop
(`pss_newton_run()`) was run on its own **unforced** adaptive time grid,
unlike every other period run in this engine (the baseline, every
finite-difference perturbation column, and the Floquet re-evaluation are all
explicitly forced onto the same accepted-time grid -- see the file's own
"SUSPECTED CAUSE 2" design note). With the Newton step's own amplification
factor (~166x here) turned way up by the near-unity multiplier, the small,
otherwise-harmless grid-selection noise between two *separately* adaptively-
stepped runs of a circuit with steep `tanh` switch transitions became large
enough to swamp the (small, genuinely converging) true residual improvement
-- so the line search could never find a trial point the comparison
recognised as better, at any step size.

**Fix**: force the damped-step trial run onto the same accepted grid too
(`pss_run_period_ex(..., force=1)`, was `force=0`) -- removing that
comparison noise. Confirmed on the repro case: converges in 4 iterations to
`residual=1.2e-7` (was `FAILED`, plateaued `~3e-4`-`6e-5`), spectral radius
`0.994` (correctly reported stable). Re-ran the **entire existing** PSS
regression suite (`tests/test_pss.py`, R2/R3/R4 nominal, R3 long-settling,
R3 Points=1000, both red-proofs, all four determinism checks) with the fix:
**no regression** -- every other check's residual/spectral-radius/iteration
count is unchanged to the significant figures already reported, and both
red-proof decks (R4 D~=0.55, R3 x20-gain) still correctly refuse to
converge. Gated, like every other tunable in this file, by the existing
`PSS_SHOOT_NO_FORCE` debug/isolation env var (forces the old, unforced
behaviour back on for comparison -- no new env var was needed for this
fix).

## `.pac ... solver=td`

```
.pac Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff \
     {lin|dec|oct} SweepPoints Fstart Fstop [method=shoot] [solver=td] [pacpoints=N]
```

The first seven fields are exactly the `.pss` fields above (a `.pac`/`.pstb`
card with no separate `.pss` card runs the same shooting-Newton PSS first,
using these same values). The sweep spec is the same `{lin|dec|oct} N
Fstart Fstop` convention `.ac` uses. `solver="td"` (quoted, this project's
own convention throughout -- see "Netlist-card quoting" below) selects this
time-domain shooting-based engine; the default (`fd`, or unset) is the
existing, unmodified frequency-domain PAC. `method="shoot"` selects the
shooting PSS underneath (same as `.pss`).

**Forcing**: whichever device(s) in the netlist carry an `AC` magnitude
(`Vin ... DC {vin_v} AC 1`, `Iinj ... DC 0 AC 1`, the same convention `.ac`
itself uses) are the excitation -- captured once via the standard
"`CKTacLoad` at omega=1" idiom (same one the existing frequency-domain PAC
and `.ac` itself use), then propagated through the LPTV orbit at every
sweep frequency. **Exactly one AC-flagged device should be active per run**
-- every one found is summed into a single excitation vector, so two active
at once silently measures their sum, not either one alone (see "Trap:
`lgprobe`'s own hidden `AC 1`" below for a real example this pass hit).

**Output**: sideband-0 response at **every** node/branch in the circuit,
read back the ordinary way:

```
.control
run
wrdata out.txt v(vout)
.endc
```

with `v(vout)` (or any other node/branch) giving the complex small-signal
transfer function **from the AC-flagged source to that node, at unit
input** -- divide by the source's own AC magnitude if it is not 1.

### `pacpoints=N`

The shooting PSS's own convergence grid (`Points` above) and the `.pac`/
`.pstb` **recording pass**'s own demodulation grid are two different
concerns paying the same `Points` value by default -- often wastefully,
since PSS convergence usually needs far fewer points than an accurate
sideband-0 extraction does. `pacpoints=N` decouples them: PSS converges at
its own (coarse, fast) `Points`, and the recording pass re-runs at `N`
instead. Measured on `pac_lti.cir` (task 2's 1e-6 target): PSS `Points=2000`
+ `pacpoints=200000` reaches the same 8.2e-9 relative-error result PSS
`Points=200000` alone did, in **1.3s instead of 113s** (~87x) -- the
shooting-Newton convergence phase was the entire cost, not this engine's own
sweep (which is `O(2 period reruns)` total regardless of sweep length or
resolution, by design). `PAC_RECORD_POINTS` (an env var) is the deprecated
predecessor, still read if `pacpoints=` is unset.

### Accuracy

Sideband-0 accuracy is governed by two **separate** effects, and raising
`pacpoints=`/`Points` only helps with the first:

1. **Trapezoidal discretisation, `O(h^2)`** -- shrinks quickly with finer
   `Points`: confirmed `~1.6e-2` at Points=200, `~1.6e-5` at Points=20000,
   `~1.6e-7` at Points=200000 on a trivial RC (`pac_lti.cir`).
2. **A separate, non-vanishing floor**, believed to come from this engine's
   own documented trade-off (`pacshoot.c`'s header): the per-step recursion
   reuses the real transient's own last-Newton-iterate factorization
   (`ckt->CKTmatrix` as `NIiter()` left it) rather than refactoring exactly
   at the accepted point -- an opt-in fix exists
   (`PAC_EXACT_JACOBIAN`, env var, default off, "not proven safe... only
   exercised end to end... not by the same kind of targeted, single-device
   validation" as the rest of this file -- see `pacshoot.c`'s own comment).
   Measured on `pstb_linear.cir`: refining Points 50000->300000 (6x) only
   cut the loop-gain relative error ~6x, not the ~36x an `O(h^2)` law would
   give -- an `O(1/Points)` signature, not `O(h^2)`.

Both checks in `tests/test_pac.py` (`pac_sb0(a)`/`pac_sb0(b)`) therefore
assert the accuracy **actually measured** at a grid this project's own test
suite can afford routinely, not a blind spec number -- see the worked
examples below for the real figures.

**A separate floor was found this pass, correlated with how deeply the
true signal is attenuated, not with frequency as such** -- confirmed **not**
reduced by `pacpoints=` (tried 200 to 200000, no measurable change) nor by
`PAC_EXACT_JACOBIAN` (also no measurable change, on a synthetic reproduction
of the same symptom, an LC filter behind a switch with no feedback at all --
see the card #465 report). For control-to-output/output-impedance (true
signal never drops much past -40dB across this project's own 500Hz-~200kHz
sweep), this only ever shows up at the single grid point closest to fs/2
(199054 Hz on a 500kHz-switching converter), 0.1-0.6dB/1-2.4deg there.
**Line-to-output is a much harder case for this floor**: closed-loop line
rejection pushes the true signal down to -45..-73dB across most of the same
sweep (strong loop gain at low/mid frequency, see the worked example below),
deep enough that 3-4 of the 14 points (roughly 50kHz and up, on both R3 and
R4) sit in this floor's own range and need the same bounded-exclusion
treatment `tests/test_pac.py` already applies to the single near-fs/2 point
elsewhere. Every point NOT in this floor's range, across every quantity this
project's own regression checks, matches its reference within 0.5dB and
typically within 0.1-0.3dB.

### Trap: `lgprobe`'s own hidden `AC 1`

`stb/lgprobe.lib` (this project's Tian loop-gain probe subcircuit, used by
`tools/stb.py`/native `.stb`/`.pstb`) is:

```
.subckt lgprobe a b
Ii 0 x DC 0 AC 0
Vi x a DC 0 AC 1
Vb b x DC 0
.ends lgprobe
```

`Vi` is **unconditionally** `AC 1` -- by design, for `.stb`'s own
alter-driven two-run protocol (Vi AC 1/Ii AC 0 on one run, the reverse on
the other). If a deck's `XP fbtop vout lgprobe` instance is left in place
for a **`.pac`** run (not `.pstb`), this hidden `AC 1` is picked up **in
addition to** whatever source the deck author actually intended to drive --
silently summing a Tian-style series injection, right inside the feedback
path, on top of the real measurement. This pass hit exactly this: an early
version of the `pac_sb0(c)` control-to-output/line-to-output/output-
impedance fixtures below reused `refs/r3_vm_buck.cir`'s/`r4_pcm_buck.cir`'s
`XP ... lgprobe` unmodified and read back **8-58 dB** off the FRA reference
before this was found and fixed (swap `XP fbtop vout lgprobe` for a plain
`Vprobe fbtop vout DC 0` -- same DC behaviour, no hidden AC flag -- exactly
the substitution `.pstb`'s own `probe=` mechanism already needed for the
unrelated reason that native `probe=<Vname>` does not support `lgprobe`'s
3-device topology at all). **If a `.pac` result looks like an open-loop
response despite a working feedback loop, check for a stray AC-flagged
device first.**

### Trap: a perturbation source with a nonzero DC bias needs its own DC in the SIN offset

Not a bug in this engine -- a bug **this pass found in the finite-
difference FRA reference tooling** (`tools/fra.py`'s `_sin_spec()`), while
cross-checking `.pac`'s own line-to-output measurement against
`~/Work/data/ngspice_augmentation/fra/curves-fork-s10/SUMMARY.md`'s
published "line-to-output" curve and finding an **enormous** (30-60dB)
disagreement that did **not** shrink with `pacpoints=`/`Points` the way a
genuine engine floor would, and that `.pac`'s own value matched simple
feedback theory (`H_closed = H_open/(1+T)`, using this project's own
already-validated loop-gain measurement) far better than the published FRA
number did.

`tools/fra.py`'s `_sin_spec(amplitude, f, td)` builds a SPICE
`SIN(0 amplitude f td 0 0)` literal -- **VO (the DC offset) is hardcoded to
0**. `alter @dev[sin] = ...` **replaces** the device's transient function
entirely, so for a source whose own DC value is legitimately 0 (`Vctrl`,
`Iinj`, `lgprobe`'s own `Vi`) this is correct. **`Vin` is not one of those**
-- it carries the converter's real 12V (R3) / 11V (R4) supply. `alter
@vin[sin] = [0 0.01 500 <td> 0 0]` silently drops that 12V bias, turning the
"perturbed" transient into a converter whose input rail toggles between
-0.01V and +0.01V -- not a small-signal perturbation of the real operating
point at all, but an entirely different (and, for a power converter,
nonsensical) large-signal condition.

**Confirmed the fix**: reproducing the SAME measurement with `VO=vin_v`
(the deck's own real DC input) instead of `VO=0` gives `-38.597dB` at 500Hz
on R3 (vs. the published `+0.032dB`) -- matching both `.pac`'s own
`-38.87dB` and the independent feedback-theory estimate
(`H_open(500Hz)~=D~=0.275 (-11.2dB)`, loop gain `T(500Hz)=27.7dB` already
validated by this project's own `.pstb` check, `H_closed ~= H_open/(1+T) ~=
-38.9dB`) to within a few tenths of a dB. **`tools/fra.py` was not edited**
(this project's `tools/fra_curves_fork.py` establishes the "read-only,
never edit the shared FRA tooling" convention this pass follows too) -- a
corrected reference curve was generated with a small, separate script
reusing `tools/fra.py`'s own lower-level helpers (`window_times`,
`extract_component`, `injected_reference_phasor`, `subtract`,
`check_uniform`, `read_tran`) with an explicit, correct `VO`, the same
pattern `tools/fra_curves_fork.py` itself already uses. **Any future
line-to-output (or other nonzero-DC-source) FRA measurement using
`tools/fra.py`'s `alter @dev[sin]` convention should pass an explicit VO,
not rely on the default.**

## `.pstb`

```
.pstb Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff \
      {lin|dec|oct} SweepPoints Fstart Fstop probe="Vname" [method=shoot]
```

Same first-seven-fields/sweep-spec convention as `.pac` above. Uses Tian's
two-injection method (same maths as `.stb`, `research/cadence-stb-pstb.md`
sec 2.4) on the LPTV orbit instead of a single DC operating point.

### `probe=` forms (card #717: parity with `.stb`)

| Form | Meaning |
|---|---|
| `probe="Vname"` | an existing plain top-level `DC 0` V-source already in series with the loop-break wire |
| `probe="<inst>:<terminal>"` | auto-insert: at netlist-parse time the named terminal of an existing instance is moved onto a new internal node, with a new `DC 0` V-source between that node and the original net (`a` = the instance side, `b` = the original net). This is the same splice (`stb_probe_splice()`, `inp2dot.c`) as `.stb`, so terminal names are the device's own (`R+`/`R-` for a resistor, `C+`/`C-`, ...), the instance must appear before the `.pstb` card, and a probe already spliced for `.stb` on the same terminal is reused (from the code; not tested) |
| `probe="<Xname>"` | an instance of the `lgprobe` subcircuit (`stb/lgprobe.lib`, `Xname a b lgprobe`): the loop is read at its `Vi` (voltage injection, branch current) and `x` (current injection, voltage) exactly as `tools/stb.py` does; `a` is downstream, `b` upstream. The subcircuit's own `AC 1` flag plays no role (`.pstb` drives only its own injections) |

All three give the same loop gain for the same circuit. Checked by
`tests/test_pstb_probe.py` (relative difference of `T` against the plain
probe: linear loop 2.4e-11 / 0 for `inst:terminal` / `lgprobe`; the
`refs/r3_vm_buck.cir` closed-loop buck 0 / 1.8e-11, 7 points 1 kHz-100 kHz).
Caveat for switching circuits: the comparison is only at the level of the
PSS orbit noise when the two decks list their elements in the same order --
moving the probe line elsewhere in an otherwise identical r3 deck changes `T`
by about 1.3e-7, not because of the probe but because node/pivot order moves
the shooting solution within its tolerance. Compare like with like.

Refused with a message (no loop gain is produced): an unknown instance or
terminal in `inst:terminal`; a probe name that is neither a V-source nor a
complete `lgprobe` instance (`v.<X>.vi`, `v.<X>.vb`, `i.<X>.ii` with the
`Ii 0 x; Vi x a; Vb b x` topology); the differential `probe=Vp,Vn` form; and
device-based (gain-nulling) probes. A probe with a terminal at ground is
accepted but reports `T = 0`, as in `.stb`.

**Output**: `loopgain` and `t` (Spectre's own naming, `loopGain = -T`, same
sign convention native `.stb` and `tools/stb.py` both already use) -- read
back with `wrdata out.txt t` (or `loopgain`). A phase-margin line is also
printed to stdout:

```
pstb: phase margin 89.38 deg at 5646.41 Hz
```

## `.pxf` (periodic transfer function, card #499)

```
.pxf Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff OutNode \
     {lin|dec|oct} SweepPoints Fstart Fstop [maxsideband] \
     [method="shoot"] [solver="td"] [sources="all"|"Vname"|"Iname"] [pacpoints=N]
```

Spectre's `pxf`: the sideband-0 transfer function from **every** requested
independent source to **one fixed output**, on the shooting PSS orbit --
the PAC/PSTB engine above answers "what does source X do to every node",
`.pxf` answers "what does every source do to node Y", one card, one run.
Card #499 (plan T5.4), private fork only, same project ruling #447 A
restriction (smooth switch models only) as `.pac`/`.pstb` above, since
it's built on the exact same engine (`pacshoot.c`'s `pac_engine_run()`,
completely unchanged -- see "Design decision" below).

**This `.pxf` dot-card already existed** in this fork's base (Enhancement-
125, `dcpss.c`'s `pxf_sweep()`, the frequency-domain harmonic-conversion-
matrix engine's own adjoint) -- **unchanged, and still the default**
(`solver="fd"`, or `method`/`solver` omitted entirely): it solves one
combined transfer function from whichever source(s) the deck itself
AC-flags, summed, across every conversion sideband. `method="shoot"
solver="td"` (both required, same gating `.pac`/`.pstb` already use)
selects the time-domain shooting-based engine documented here instead --
additive, same pattern as `.pac`'s own `solver=td` on top of the existing
frequency-domain `.pac`. **Field order was kept from the pre-existing
card** (`OutNode` before the sweep spec, matching `.pss`'s own `OscNode`
placement) rather than adopting a new, incompatible grammar for this
pass's own suggested syntax (`v(out[,ref])` after the sweep) -- see
"Syntax decisions" below for the full reasoning.

**`OutNode`**: a bare node name, **ground-referenced only** (no
`v(out,ref)` differential output in this pass -- matches the pre-existing
`PxOutNode` field's own "(reference = ground)" comment; a future pass
could add a second node and difference the two rows, not attempted here).

**`sources=`**: which independent V/I source instance(s) to sweep,
`solver="td"` only. Repeatable (`sources="V1" sources="I2"`), **not** a
comma list -- same reason `.stb`'s `sens=` isn't one (docs/sens.md,
card #457's tokenizer comma-splitting hazard). `sources="all"` (also the
default when `sources=` is omitted entirely) sweeps **every** native
`Vsource`/`Isource` instance in the circuit, regardless of whether the
deck itself gives it a netlist `AC` value -- confirmed against Spectre's
own `pxf` semantics (WebSearch, card #499): `.pxf` does not need a
source's own `AC` flag, unlike `.pac`/`.ac`.

**Output**: one complex vector per swept source, named `pxf_<source>`
(e.g. `pxf_v1`, `pxf_iinj`), in one plot -- read back the ordinary way:

```
.control
run
wrdata out.txt pxf_vin pxf_vctrl
.endc
```

### Syntax decisions (what was decided, and why)

1. **`OutNode` position and form kept from the pre-existing dot-card**
   (bare node, before the sweep, ground reference only) instead of this
   pass's own initially-suggested `v(out[,ref])`-after-the-sweep form --
   the project rule "one name per thing... never rename" is written for
   projects/repos, not dot-cards, but the same spirit applies: `.pxf`
   already existed and worked (for the frequency-domain engine); inventing
   a second, incompatible `.pxf` grammar for the SAME card name would have
   broken every existing `.pxf` deck's field order for no reason a v1
   differential-output feature needs. Extending it additively (new
   trailing tokens only) was the lower-risk, equally expressive choice.
2. **`sources=` repeated, not comma-separated** -- established project
   convention (`sens=`, card #457/#475), same tokenizer hazard.
3. **`sources="all"` default** -- matches Spectre's own "every independent
   source" `pxf` default; explicit `sources=` narrows it.
4. **One adjoint (transposed) solve, or a batched forward sweep?** See
   `pacshoot.c`'s `pac_td_pxf()` design comment for the full derivation of
   both options and which one shipped -- short version: a genuine
   backward/transposed sweep of this file's own per-step recursion was
   derived (needs a new `SMPtSolve`-style real-transpose wrapper around
   Sparse's `spSolveTransposed`, confirmed available in this build, plus a
   transpose of the per-step "history" operator complicated by order-2
   trapezoidal's own 2-term recursion, plus rebuilding every `J_k` out of
   order) but judged too much NEW, unvalidated adjoint machinery for the
   time available, given this project's OWN saltation work is a direct
   demonstration of what an under-derived transposed correction costs here
   (months, several passes, a never-root-caused ~2.2x scaling bug). This
   pass instead reuses `pac_engine_run()` (the SAME already-validated
   engine `.pac solver=td` uses) once per requested source -- still
   eliminates the dominant cost "one `.pac` per source" would pay (a full
   PSS shooting-Newton RECONVERGENCE per source; PSS converges exactly
   ONCE per `.pxf` job here) at the cost of not sharing the source-
   independent monodromy (`Phi`) computation across sources. A true
   backward-adjoint engine remains valid future work if source count ever
   makes this pass's `O(K)` cost the bottleneck.

### Validation (card #499)

Every `.pxf`-vs-`.pac` comparison below reuses the literally IDENTICAL
`pac_engine_run()` code path (`pac_td_pxf()` just supplies a different,
per-source isolated forcing pattern), so agreement is at floating-point
precision, not merely inside the card's 1e-6 target -- confirming the
per-source isolation (`pac_capture_source_b0()`, zero every OTHER
source's AC stamp, capture the one requested via the SAME
`CKTload+CKTacLoad` idiom `.pac`'s own forcing capture uses) exactly
reproduces what a real single-source `AC 1` netlist would give `.pac`:

| fixture | sources | vs `.pac` (one source) | vs `.ac`/`.tf` |
|---|---|---|---|
| `pxf_lti.cir` (LTI, unrelated clock for PSS) | Vin, Iinj | 0.0 (exact) | 8.225e-09 (target 1e-6) |
| `pxf_r2_switched_rc.cir` (switched-RC LPTV, ideal `sw`) | Vin, Iinj | 0.0 (exact) | -- |
| `pac_r3_ctrl2out.cir` (R3 voltage-mode buck, smooth switches) | Vin, Vctrl | 0.0 (exact) | -- |

(`pxf_r2_switched_rc.cir`'s ideal `sw` is safe here for the same reason
`pac_sb0(b)`'s existing check already relies on: its control node is
driven only by the fixed clock, never by Vin/Iinj, so saltation's own
correction coefficient is identically zero -- confirmed no-op on this
topology already, saltation section above.) The `.ac`/`.tf` cross-check
(the card's separate "an LTI deck where it equals plain .ac/.tf"
requirement) reuses `pxf_lti.cir`; its 8.225e-09 figure is the SAME number
`pac_sb0(a)` already reports for the underlying `.pac solver=td` engine at
`PAC_RECORD_POINTS=200000` (`test_pac.py`'s own accuracy-floor note
applies unchanged -- `.pxf` inherits it, does not add to it).

`tests/test_pxf.py` also checks the `sources=` default (omitted entirely
= `"all"`, confirmed to find all three independent sources on
`pxf_lti.cir` -- Vclk, Vin, Iinj -- and to match the explicit
`sources="Vin" sources="Iinj"` run exactly).

Full regression (`test_pac.py`, `test_pss.py`, `test_stb_native.py`,
`test_checkpoint.py`, `test_sens.py`) re-run against this card's build:
no change from the documented baselines (see the card #499 report for the
exact counts) -- this feature is purely additive (a new dot-card branch
and a new static function in `pacshoot.c`; zero lines changed in
`pac_engine_run()`, `PSSpacRecordStep()`, `pac_reset()`, or any other
already-shipped `.pac`/`.pstb` code path).

### Independent cyclostationary validation (card #509, second pass, 2026-09-29)

Result: **no formulation change needed** -- `T * Sum_k dens_k / w_k` (the 1/w_k weighting, modulation and folding) agrees with two independent references on genuinely cyclostationary (folding) circuits. A first Monte Carlo attempt disagreed by +1.1..+1.7 dB, and that was entirely a fault of the Monte Carlo, found and removed below.

**(a) Exact analytic, track-and-hold** (`tests/fixtures/pnoise_track_hold_slow.cir`: Rtrack=25k, C=10p, tau=250 ns, fclk=100 kHz, 20 % track window = 8 tau; Rtrack's noise, bandwidth 640 kHz, folds from ~6 harmonics of fs). The held output is an Ornstein-Uhlenbeck process (variance kT/C) restarted from the held value at each track start. Its low-frequency one-sided PSD from the per-cycle integral Z_n = A_n + T_h*s_n and its lag-1/2 covariances is 7.7806e-8 V/rtHz (`th_analytic_S0()` in tests/test_pnoise.py). The naive "zero-order hold of a kT/C sample" value 9.104e-8 is 1.38 dB HIGHER and is wrong here: the sample is correlated with the track-phase samples before it (an analytic trap, not an engine error). `.pnoise` at 300 Hz/1 kHz/2 kHz/3 kHz (Points=2000): -0.017/-0.016/-0.015/-0.013 dB from the exact value (test `pnoise_track_hold_analytic`, limit 0.1 dB).

**(b) trnoise Monte Carlo** (`tools/pnoise_mc.py`, run with `NGSPICE=... python3 tools/pnoise_mc.py all 32 5`): each resistor is made `noisy=0` and a `TRNOISE` current source is placed in parallel with it (terminals as NevalSrc), NA = sqrt(S0/(2 TS)) for a one-sided flat PSD S0 = 4kT/R (iid N(0,NA^2) samples every TS, linearly interpolated: PSD 2*NA^2*TS*sinc^4(f TS)); T = 300.15 K in both. 32 seeds x 20 ms, Welch (Hann, 50 % overlap, 8192 points of 200 ns = 610 Hz bins), 640 segments. Expected spread stated before comparing: single bin +0.17/-0.18 dB (1 sigma, chi-squared with 2 dof averaged over 640 segments); a 9-bin band average has n_eff = 640*9/(1.83*1.056) (adjacent Hann bins correlate 0.44) so **1 sigma = 0.08 dB**. Each row averages the MC power over 9 consecutive bins and compares with the mean of `.pnoise` evaluated at exactly those 9 bins:

| centre (Hz) | RC calibration (analytic) | r2_smooth (`.pnoise`) | track_hold_slow (`.pnoise`) |
|---|---|---|---|
| 4272 | +0.12 | -0.22 | +0.03 |
| 9766 | -0.23 | +0.02 | -0.10 |
| 15869 | +0.04 | -0.09 | -0.03 |
| 21973 | -0.15 | +0.05 | -0.10 |
| 28076 | -0.06 | +0.03 | -0.09 |
| 34790 | +0.09 | -0.09 | -0.10 |
| 42725 | +0.17 | -0.03 | +0.08 |
| worst / mean (dB) | 0.23 / -0.003 | 0.22 / -0.05 | 0.10 / -0.05 |

(MC minus reference, dB in power; 20*log10 of the amplitude ratio is the same number.) Every band is inside 3 sigma of zero and all worst values are inside the card's 0.5 dB target; the RC row is the calibration of the whole MC chain (NA convention, boxcar, Welch, normalisation) against a known LTI answer, so the MC chain itself is good to about +-0.2 dB per band, mean -0.003. Data: `~/Work/data/ngspice_augmentation/pnoise/mc_all_ts5_v2.txt`. 32 seeds cost 4.7/8.3/8.6 min (RC/r2/T&H) at 3 ngspice with OMP_NUM_THREADS=2.

Three artefacts in the Monte Carlo itself, each found by the RC calibration or a discriminating rerun, worth knowing before trusting any trnoise-vs-`.pnoise` comparison:
1. **Aliasing.** Point-sampling v(t) (`linearize`, dt=1 us) folded the switched-RC's >500 kHz noise (Rsw*C1 pole at 1.6 MHz) into the band: +1.1..+1.7 dB, flat. Fix: the periodogram input is the boxcar average over each 200 ns sample (differences of `integ(v)`), whose sinc nulls sit on the alias frequencies.
2. **Time-step bias.** With maxstep = TS = 20 ns (tau 99-250 ns) the transient integration itself biased the noise by -0.3 dB (both circuits, same sign); TS = 5 ns removes it (T&H mean -0.26 -> -0.05 dB). Use TS <= tau/20.
3. **Numerical floor.** `integ` of a ~1 V level over >1e6 steps leaves a rounding floor ~2 nV/rtHz (RC calibration +4.5 dB at 43 kHz). Fix: integrate v(noisy) - v(noiseless replica); the two copies of the circuit share one transient (so identical time steps), which also cancels the ~1 V periodic orbit without a second run. (A separate deterministic run with a different time grid leaked edge energy into the high band.)

Not covered by this second pass: MC and analytic checks use resistor thermal noise only (constant S; folding through the LPTV transfer, not a time-varying S); MC needs `PNOISE_MC=1` (about 30 min) and is off in the default test run.

### What was skipped

- **`pnoise`** (periodic noise) is unrelated scope, not attempted.
- **Sidebands other than 0** -- same sideband-0-only scope as `.pac`
  `solver=td` above; Spectre's own `pxf` supports sweeping the OTHER
  sidebands of a source's contribution too (WebSearch, card #499), not
  implemented here.
- **A true backward-adjoint engine** -- derived (see "Syntax decisions"
  item 4 above) and NOW IMPLEMENTED, `.pxf ... adjoint=1` (card #504) --
  see "`adjoint=1`: a true backward-adjoint engine (card #504)" below.
  The batched-forward-per-source design above (`adjoint=0`, still the
  default for a small source count) stays available and is what the
  adjoint engine itself falls back to for saltation/free-running
  circuits (see that section's own scope cut).
- **Differential output** (`v(out,ref)`) -- ground-referenced only.
- **No new fixture with more than 2 independent sources was built** --
  every validation fixture above sweeps exactly 2 (plus the bare-`.pxf`
  default-`sources=all` check, which finds 3 on `pxf_lti.cir` including
  the clock source, but is not separately cross-checked against `.pac`
  for the clock source itself). Card #504's own `pxf_manysrc.cir` (below)
  is the exception, but it exists for TIMING only (50 near-identical
  current injectors, not a topology cross-check).

## `adjoint=1`: a true backward-adjoint engine (card #504)

```
.pxf ... method="shoot" solver="td" sources=... [adjoint=1|0]
```

`adjoint=1` selects `pacshoot.c`'s `pac_engine_run_adjoint()` in place of
the `adjoint=0` batched-forward-per-source loop above (`pac_td_pxf()`'s
original card #499 design, still the default for a small source count --
see `PXF_ADJOINT_AUTO_THRESHOLD` in `pacshoot.c`, currently 4: more than 4
requested sources auto-selects the adjoint engine unless `adjoint=` says
otherwise explicitly). Both give the SAME answer (this section's whole
validation is proving that); the difference is cost: `adjoint=0` pays one
extra pair of period reruns PER requested source (`O(K)`), `adjoint=1`
pays a small, FIXED number of extra passes REGARDLESS of source count
(`O(1)` in K) -- see "Design" and "Timing" below.

### Design

For one fixed output row and frequency `w`, `adjoint=0`'s own two-pass
engine (`pac_engine_run()`, this file's header) computes, for a
source's own isolated forcing pattern `b`:

```
Y0(w,b) = xi^T G(w) + h(w,b)
xi       = (e^{jwT}I - Phi)^{-1} z(w)
z(w)     = E_read(u_N(b))              (pass 1's own b-forced, zero-IC endpoint)
h(w,b)   = demod(u_k(b)[out])          (pass 2's own b-forced, zero-IC demodulation)
```

`z(w)` and `h(w,b)` are BOTH read off the literal SAME trajectory (the
b-driven, zero-initial-condition particular solution of the per-step
recursion -- pass 1's forced column and pass 2's own forcing use the
identical `b`, same per-step map, and re-running from the converged `x0`
reproduces the transient bit for bit, this file's own header). Substituting
gives one single linear functional of `b`:

```
Y0(w,b) = [E_read(u_N(b))]^T M(w)  +  (1/T) sum_k alpha_k[out] u_k(b)[out]
M(w) := ((e^{jwT}I - Phi)^T)^{-1} G(w)          (the TASK's own "transpose
                                                  of the periodic boundary
                                                  condition")
```

`G(w)[d]` is "the `Y0` a homogeneous-IC=`e_d`, `b=0` pass 2 would give" --
computed ONCE (batched over every frequency, `M=n*nfreq` columns, one
period run) alongside `Phi` (`M=n` columns, one more period run) -- NEITHER
depends on the source, so both are shared across every requested source.
`alpha_k[out]` is the raw trapezoidal-plus-wraparound demodulation weight
`PSSpacRecordStep()`/`pac_finish_demod_wraparound()` already compute,
re-derived in closed periodic form: `alpha_k = ph2_k * 0.5*(h_k +
h_{k+1 mod N})` (`h_{N+1} := h_1`) -- the expected periodic-trapezoidal
quadrature weight (a genuine cross-check on the derivation: it fell out of
re-deriving the existing code's own formulas, not assumed from scratch).

Since `u_k(b)` is generated by ONE linear per-step recursion with `b`
entering identically (phase-modulated) at every step, `Y0(w,.)` is a single
COMPLEX-LINEAR functional of `b` -- representable as `Y0(w,b) = b^T
Lambda(w)` for ONE `Lambda(w)` in `C^neq` that does not depend on `b` (or
on which source) at all. `Lambda(w)` is computed by a BACKWARD sweep of the
per-step recursion's adjoint (`pac_engine_run_adjoint()`'s own long design
comment in `pacshoot.c` has the full per-step derivation: `SMPtSolve()`,
the new REAL transposed-solve wrapper this card added around Sparse's
`spSolveTransposed`/KLU's `klu_tsolve` -- dispatched exactly like the
existing complex `SMPcaSolve()` docs/sens.md's card #475 adjoint already
uses; the trapezoidal order-2 history term's OWN two-term recursion,
transposed into `beta`/`gamma` "gather" formulas that are the direct
mirror of `PSSpacRecordStep()`'s existing `dceq`/`dveq`/`T_k` "scatter"
formulas, re-derived FROM that code rather than from physics intuition, to
avoid a docs/sens.md-style transpose-or-sign bug). Folding `M(w)` in only
at the LAST step and `alpha_k[out]` at EVERY step gives `Lambda(w)` in
ONE traversal; every requested source's PXF is then a single dot product
`b_source^T Lambda(w)` -- the "K sources for the cost of ~1 sweep" payoff.

`J_k` (the per-step Jacobian, needed in REVERSE order): the backward sweep
does not store or clone factorizations. It RECOMPUTES `J_k` on demand from
a recorded trajectory (`g_pac_traj[]`: `x_k` plus
`t_k`/`h_k`/`ag0_k`/`ag1_k`/`order_k`, captured during the engine's OWN
"`Phi` + trajectory" forward pass) via a fresh `CKTload()`+`SMPluFac()` per
backward step -- mirroring `pac_refactor_at_accepted()`'s existing,
validated "put `x_k` where devices read node voltages from, reload,
refactor" recipe, but standalone (no live transient to protect; the
handful of `ckt`-global scalars this touches are saved/restored around the
whole backward sweep). See "Memory" below for why recompute was chosen
over storing factorizations, and what it costs.

Because this recompute-based `J_k` is internally consistent only with
itself (`x_k`, not the live transient's possibly-stale `x_{k-1}` the
DEFAULT, non-`PAC_EXACT_JACOBIAN` `adjoint=0` path uses -- see
`pac_refactor_at_accepted()`'s own header in `pacshoot.c`), the adjoint
engine's OWN forward `Phi`/`G` passes ALSO force the exact-Jacobian
convention (independent of the `PAC_EXACT_JACOBIAN` env var), so its own
forward and backward halves use the IDENTICAL `J_k` sequence -- a
well-defined adjoint of ONE map, not an accidental mix of two.
**Consequence**: comparing `adjoint=1` against `adjoint=0`'s own DEFAULT
(stale-Jacobian, `O(1/Points)`-accurate per this file's earlier accuracy
section) output is only accurate to `O(1/Points)`, not floating-point
noise -- `tests/test_pxf.py`'s adjoint checks therefore run `adjoint=0`
ALSO under `PAC_EXACT_JACOBIAN=1` for the comparison (both sides then
mathematically exact, agreement limited only by roundoff). See
"Validation" below for the numbers both ways.

### Scope cut (refused, not silently wrong)

Saltation (ideal-switch/event) corrections and the free-running
(`trigger=`) period-variation projection are NOT transposed by this
engine. Detected via `PSSsaltEventsSeen()>0` (checked right after the
`Phi` pass) or `trig_idx>=0` (the free-running caller) and returned as
`PAC_ADJOINT_REFUSED`, so `pac_td_pxf()` falls back to the `adjoint=0`
per-source engine transparently, with a one-line stderr note -- correctness
is preserved on every fixture: `pxf_r2_switched_rc.cir` (the one
`test_pxf.py` fixture with real saltation events -- its `Sswitch` flips
twice a period) is answered by `adjoint=0` even when `adjoint=1`/auto picks
this engine, so the done-when's "== v1" is trivially true there by
construction, not because saltation was actually transposed. A saltation-
or free-running-aware adjoint remains valid future work.

### Memory

The recorded trajectory is `O(N*(neq+2))` doubles (`N` = accepted steps in
one recording-pass period, `neq` = full circuit unknown count) -- the
`Phi`/`G` forward pass's own `g_pac_traj_append()` call per accepted step.
Chosen over deep-cloning/storing `N` factorizations (the "task's own"
other listed option) because: (1) it needed no new understanding of
Sparse's/KLU's internal factored-matrix representation (`CKTload()`+
`SMPluFac()` is an already-existing, already-validated recipe --
`pac_refactor_at_accepted()`, unchanged, just invoked out of the live
transient's own time order); (2) for every fixture this card measured
(`neq` in the tens, `N` from a few hundred to a few thousand at this
project's usual `pacpoints=`/`Points=`), the trajectory is at most a few
MB -- see the numbers table below. A much larger circuit or a much finer
recording grid (e.g. `pac_lti.cir`'s own `Points=200000` accuracy-floor
fixture, this file's earlier "Accuracy" section) would grow this linearly;
storing factorizations instead remains valid future work if per-step
REFACTOR cost (not memory) ever dominates for a circuit far larger than
this project's own fixtures.

### Validation (card #504)

`tests/test_pxf.py`'s `check_pxf_adjoint_*()` checks, against this card's
own build (`aaba50ba3f49-adj-1`, branch `aug-adj` = pin + `patches/0001`-
`0053` + this card's own diff), `adjoint=1` vs `adjoint=0` (both under
`PAC_EXACT_JACOBIAN=1` -- see "Design" above for why that's the fair
comparison, not `adjoint=0`'s own default):

| fixture | sources | worst relative error | target |
|---|---|---|---|
| `pxf_lti.cir` (LTI, unrelated clock for PSS) | Vin, Iinj | **2.380e-14** | 1e-8 |
| `pac_r3_ctrl2out.cir` (R3 voltage-mode buck, smooth switches, has an inductor AND a capacitor -- `C1f` -- whose negNode is NOT ground) | Vin, Vctrl | **4.014e-09** | 1e-8 |
| `pxf_r2_switched_rc.cir` (switched-RC LPTV, ideal `sw` -- real saltation events, `Sswitch` flips twice a period) | Vin, Iinj | **0.0 (exact)** | -- (adjoint REFUSES and falls back to `adjoint=0` -- see "Scope cut" above; the fallback message is checked present, and the numbers are then trivially the SAME run, not a separate cross-check) |

Both floating-point-exact fixtures (`pxf_lti.cir`, well inside 1e-8;
`pac_r3_ctrl2out.cir`, also inside 1e-8, a real nonlinear small-signal
matrix from BSIM/behavioural-switch stamps, not just linear R/L/C) --
"or explain a tighter floor" was not needed. Every source on every
fixture carries no netlist `AC` value (same as every other `.pxf` fixture
in this file), confirming "sources with no AC value" is exercised, not
just claimed.

**A real bug found and fixed validating this** (both are in this card's
own diff, not shipped separately as their own patch): (1) the per-step
demodulation seed in the backward sweep was missing the `(1/T)`
normalisation `G(w)` already carries internally -- caught immediately (a
gross, ~40-100% error, not a subtle one) by the `pxf_lti.cir` check before
this section existed. (2) `M(w)`'s own scatter into the backward sweep's
LAST step used the per-step recursion's internal "posNode minus negNode"
convention instead of `z(w)`/`Phi`/`G(w)`'s own "posNode only"
(`pac_state_row()`) convention -- invisible on `pxf_lti.cir`/
`pxf_r2_switched_rc.cir` (every capacitor there has `negNode=0`, so the
two conventions coincide) but a **106% error on `pac_r3_ctrl2out.cir`**
(`C1f` is between `compin` and `opout`, neither ground) -- exactly the
kind of fixture-coverage argument this project's own docs/sens.md report
makes for why more than one topology needs checking. See `pac_engine_run_adjoint()`'s own code comments in `pacshoot.c` for both
fixes in place -- the `(1/T) normalisation` comment right before the
per-step demod `coef` seed (fix 1), and the `card #504 bug fix` comment on
the `is_last` block's `M(w)` scatter (fix 2).

### Native path: J_k stored, loaded at x_k (card #541, 2026-10-02)

Circuits without model charges (cap-less MOS1, ideal switches, native C/L) take the native path. Until #541 its backward sweep re-loaded
`J_k` from the devices in reverse order, and the forward "exact" `J_k` (`PAC_EXACT_JACOBIAN`, forced by the adjoint engine) was loaded in
`MODETRAN|MODEINITPRED`, the mode `dctran` has already switched to for the next step: MOS1 and similar devices then evaluate at the
voltages PREDICTED from `CKTstate1/2`, not at `x_k`. Symptoms (#539 bug 3): `adjoint=1` ~0.02 vs forward 0.53 on an inverter chain whose
source sets a MOS threshold crossing, and the exact-J forward 0.4462 vs the default-J forward 0.4677 (5 %). Now `pac_refactor_at_accepted()` loads
in `MODEINITFLOAT`, bypass off (as `pac_gen_loads()`), and the adjoint recording keeps the loaded `J_k` per step (16 B per stored entry,
same `PAC_TRAJ_MAX_MB` cap) for the transpose. Result: adjoint == forward to 1e-16, exact-J == default-J forward to 2e-6 on that chain.
`PAC_EXACT_JACOBIAN=1` on other native decks (R2/R3/R4 switched fixtures) goes through the same load: see tests/test_pxf.py, test_pac.py.

### Timing (card #504)

`check_pxf_adjoint_timing()`, `tests/fixtures/pxf_manysrc.cir` (LTI, no
switches -- isolates the engine's own per-source cost from any
switching/saltation concern), `sweep=dec 5 1kHz-100kHz` (11 points),
`adjoint=1`/`adjoint=0` both FORCED explicitly (no auto-selection
ambiguity), wall time in seconds (whole-process `ngspice -b`, includes the
PSS convergence both engines pay identically):

| N sources | `adjoint=1` | `adjoint=0` (v1) | v1 / adjoint |
|---|---|---|---|
| 2  | 0.172s | 0.180s | 1.05x |
| 10 | 0.176s | 0.265s | 1.50x |
| 50 | 0.180s | 0.631s | 3.51x |

`adjoint=1`'s own wall time is flat within measurement noise from N=2 to
N=50 (0.172s -> 0.180s, 1.05x) -- exactly the design's O(1)-in-K claim,
not just an architectural assertion: the backward sweep and its two
forward passes (`Phi`, `G(w)`) never re-run per source, only the final
dot product does, and that is `O(neq)` per source, too cheap to show up
at this circuit's size. `adjoint=0` grows from 0.180s to 0.631s (3.51x)
over the SAME range -- sub-linear here (not the full ~25x a naive `O(K)`
model would predict) because `pac_engine_run()`'s own PER-CALL fixed
costs (PSS convergence is NOT re-run per source -- only `pac_engine_run()`
itself, i.e. the `Phi`+forced-column and re-propagate/demodulate passes,
re-runs) dominate at this small a circuit/sweep; the gap would widen
further at a longer sweep or a bigger circuit (not separately measured
this pass). The crossover point (`adjoint=1` cheaper than `adjoint=0`)
is already at N=10 here, well below `PXF_ADJOINT_AUTO_THRESHOLD`'s
default of 4 sources -- i.e. `pac_td_pxf()`'s own auto-selection is, if
anything, conservative (selects the cheaper engine even earlier than
strictly needed by this one fixture's own numbers).

**Memory**: not separately instrumented with a byte count this pass (a
real gap, not silently skipped) -- what WAS measured: the recorded
trajectory's own length (`g_pac_traj_n`, printed via `PXF_ADJ_DEBUG=1`,
the engine's own opt-in diagnostic env var, same convention as
`PAC_ORDER_TRACE`/etc elsewhere in this file) was 2015 steps on
`pxf_lti.cir` (`Points=2000`, `neq=6`: nodes `clk`/`clkgnd`/`in`/`out` +
branch equations for `Vin`/`Vclk`, hand-counted from the fixture) and 369
steps on `pac_r3_ctrl2out.cir` (`Points=200`). Both are `O(N*(neq+2))`
doubles per the design, and at these sizes (`neq` in the single/low-double
digits, `N` in the hundreds to low thousands) the trajectory is
unambiguously sub-megabyte -- ASan (which would have caught any gross
over-allocation as a build/runtime issue, not silently) raised nothing.
A byte-precise number, and a check of a MUCH finer recording grid (e.g.
this project's own `Points=200000` `pac_lti.cir` accuracy-floor fixture,
which would put `N` around 2e5 and the trajectory around a few tens of
MB by the same formula) were not measured this pass.

## Worked examples

### R3 loop gain vs. the fork's own FRA reference

`refs/r3_vm_buck.cir` (voltage-mode buck, D~=0.276, fs=500kHz) with `XP
fbtop vout lgprobe` swapped for a plain `Vprobe fbtop vout DC 0` (see the
trap above):

```
.pstb 500e3 300e-6 vout 200 10 150 1e-6 dec 5 500 200000 probe="Vprobe" method="shoot"
.control
run
wrdata pstb_out.txt t
.endc
```

Result: 14 points, 500Hz-~200kHz, worst `|T|` error **0.273 dB**, worst
phase error **1.812 deg** against the fork's own FRA loop-gain curve
(`~/Work/data/ngspice_augmentation/fra/curves-fork-s10/SUMMARY.md`, same
fork family, card #461); phase margin **38.4 deg @ 33554 Hz** (FRA: **38.46
deg @ 33463 Hz**). Both comfortably inside a 0.5dB/3deg target.

### R2 vs. the exact analytic sideband-0 value

`refs/r2_switched_rc.cir` (a genuinely switched, but clock-driven not
state-driven, single-state LPTV RC network) has a closed-form exact
sideband-0 transfer function (`refs/r2_analytic.py`'s `g0()`, independently
self-checked against a brute-force fixed-step propagation to ~3e-11
relative):

```
.pac 100e3 100e-6 vx 20000 1 25 1e-6 dec 10 1k 130k method="shoot" solver="td"
.control
run
wrdata pac_td_out.txt vx
.endc
```

Result: 22 points, 1kHz-130kHz, worst `|g0|` relative error **5.884e-05**
(target 1e-3), worst phase error **0.0017 deg** (target 0.05) -- and, at
every one of the 22 points, **at least 10x closer** to the exact `g0()` than
the standard duty-averaged small-signal model is (the averaged model is a
reasonable low-frequency approximation but visibly wrong as frequency
approaches `fclk`; this engine tracks the true LPTV response instead).

### R4, D~=0.30 vs. D~=0.45: half-switching-frequency peaking grows with duty

Peak-current-mode control without slope compensation has well-known
sampled-data dynamics that get closer to the (period-doubling) stability
boundary as `D -> 0.5` -- visible here as **growing peaking in the
control-to-output response near half the switching frequency**, confirmed
in **both** the FRA reference and this engine's own `.pac` measurement (so
this is a real physical effect the fixture's FRA data already shows, not an
artefact of either measurement method):

| | D~=0.30 @ ~199kHz | D~=0.45 @ ~199kHz |
|---|---|---|
| FRA (control-to-output) | -26.138 dB | -21.890 dB |
| `.pac solver=td` | -25.723 dB | -22.406 dB |

Both show the same direction and a similar magnitude of change (+4.2dB FRA,
+3.3dB `.pac`) as D goes from 0.30 to 0.45.

## Free-running (autonomous) shooting (`.pss ... trigger=`, card #489)

Card #489 (plan T6.2), private fork only, same project ruling #447 A
restriction (smooth switch models -- see "Known limits" below): a
`trigger=<node>` on a `.pss ... method=shoot` card selects a SECOND
shooting engine, `PSSshootFreerun()` (`src/spicelib/analysis/pssshoot.c`,
appended after the driven engine `PSSshootDriven()`), for converters whose
switching frequency is an OUTPUT, not an input -- hysteretic/bang-bang,
constant-on-time (COT), quasi-resonant (QR). Where the driven engine's
`Fguess` is the deck's own exact, known clock, here it is only a rough
period ESTIMATE (used to size the search window and the per-period
integration grid); the true period is measured.

```
.pss Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff \
     [uic] method=shoot trigger=<node> [trigval=<v>] [trigdir=rise|fall]
```

`trigger=` is the only new required field; `trigval` defaults to 0.5,
`trigdir` to `rise`. `StabTime` matters MORE here than in the driven
engine: this engine has no staged-warm-up/Broyden escalation to fall back
on, so it runs ONE plain warm-up transient of `StabTime` seconds from the
deck's own IC before Newton starts (mirroring what the field already means
on a `.pss` card), and needs to be long enough to clear a hysteretic
converter's own cold-start transient (`refs/r5_hysteretic_buck.cir`'s own
header: ~400-450us before its q latch settles into clean switching).

### The algorithm (Poincare-section shooting, SIMPLIS-POP style)

The trigger node's own crossing of `trigval` in the `trigdir` direction
defines the Poincare section `g(x,t) := v(trigger) - trigval = 0`. Unlike
a bordered-Newton design (period T as an explicit second unknown with a
separate phase condition -- considered and NOT used here), this engine
follows SIMPLIS POP / the Wong shooting patent (`research/simplis-pop-
ac.md` Sec.2, Sec.6 Stage 1) and lets T fall out of EVENT LOCATION: the
period map integrates from a given `x0` until the FIRST qualifying
trigger crossing, however long that takes (up to `MAX_PERIOD = 10 x
Fguess-estimate`, the same "free-running gets 10x, driven gets >=3x"
convention SIMPLIS's own reference manual gives), and returns both the
endpoint state and the measured period. Newton's unknown stays exactly
`x0` (n-dimensional, the same native-capacitor/inductor state vector the
driven engine uses), with the ordinary shooting residual
`F(x0) = x(T(x0)) - x0`.

Monodromy: `pss_run_until_trigger()` locates the exact crossing time by
regula-falsi/bisection refinement (restarting the sub-period from the
bracket's pre-crossing state, a handful of cheap short re-runs, each at
most one of the outer probe run's own accepted steps), then a plain
central finite difference over this SAME period map gives M:
`M[:,j] = (P(x0+h e_j) - P(x0-h e_j)) / (2h)`, with `P(x0) := x(T(x0))`
run independently for every column -- each perturbed trajectory finds ITS
OWN crossing time. This is the reason no separate terminal-projection
formula (`M = (I - f(T) grad(g)^T / (grad(g).f(T))) Psi(T)`, the analytic
alternative `research/simplis-pop-ac.md` Sec.6 Stage 1 gives) is
implemented: because every column's OUTPUT is, by construction, always
exactly ON the section (`g=0`) -- the map's image is always constrained to
the section regardless of where a perturbed `x0` itself lands -- the FD
difference quotient already estimates the derivative of the section-to-
section return map directly, mechanically reproducing what the projection
formula computes analytically. `M` is consequently RANK n-1, not n: it has
exactly one eigenvalue collapsed to (near) ZERO in the direction the
projection removes -- this is the expected signature of a working
terminal-projected monodromy, verified empirically (R5: trivial eigenvalue
3.6e-6, all others real and small, see the card #489 report).

Saltation: NONE needed, for the same reason the driven engine needs none
-- project ruling #447 A restricts both engines to smooth (continuous)
switch models, so "discretise-then-linearise is consistent" (`research/
ngspice-internals-pss-pac.md` Sec.4.5) applies to every "event" here too;
it is a smooth threshold crossing, not a true discontinuity.

Floquet: eigenvalues of M (`psseigen.c`, unchanged). Because M is already
the projected monodromy, its trivial direction sits near 0 (not 1, where
it would sit in the un-projected `Psi(T)`); this engine finds it as the
smallest-modulus eigenvalue and reports it separately, labelled "trivial
(projected flow direction)". The other n-1 eigenvalues are the genuine
Floquet multipliers; `max|lambda| < 1` over those n-1 is the stability
test (same convention the driven engine's `spectral_radius` line uses,
just excluding a different eigenvalue).

Honest failure: `pss_run_until_trigger()` returning "no crossing found
within MAX_PERIOD" is reported as FAILED, never a steady state -- this is
what makes a collapsed-hysteresis-band or stable-DC-point deck come back
"no periodic orbit" rather than silently reporting a non-oscillating DC
point as if it were one (verified, `tests/test_pss.py`'s
`check_r5_freerun_redproofs()`: a `trigval` outside the trigger node's own
achievable range, and a `trigger=` naming a nonexistent node, both refuse
cleanly with no special-case code -- both just never produce a crossing).

### Validation (R5, `refs/r5_hysteretic_buck.cir`, card #489)

`tests/fixtures/pss_r5_shoot.cir`: converges in 31 Newton iterations
(residual 6.3e-6, tol 1e-4), period 2365.85ns (fsw=422681Hz),
spectral_radius=0.035 (stable; all 3 non-trivial multipliers real and
`<0.04` in magnitude), trivial multiplier 3.6e-6 (correctly near 0).

Period vs a settled `.tran`: this deck has a fast internal state (`q`, the
latch, relaxation time constant `cq/kdrive = 1ns`) that a COARSE reference
`.tran`'s own output-grid resolution aliases badly -- measured this card:
`refs/r5_hysteretic_buck.cir`'s own 10ns tstep gives a period estimate
with +-3% cycle-to-cycle jitter (average over 432 cycles in [2,3]ms:
432.3kHz) that is a measurement artifact, not the true period. Re-running
the same deck at 2ns tstep (`tests/fixtures/pss_r5_tran_settled.cir`)
gives a much more stable 423.6kHz (0.23% from PSS-FREERUN's 422.68kHz);
re-running it AGAIN starting from PSS-FREERUN's own converged `x0` (device
`IC=` set to the converged state) at 10ps tstep gives a rock-stable
(<1e-6 relative cycle-to-cycle) 422.46kHz -- 0.05% from PSS-FREERUN's own
answer, WITHIN the 0.1% target. This monotonic convergence (10ns: 432.3kHz
w/ heavy jitter -> 2ns: 423.6kHz -> 10ps from the shooter's own x0:
422.5kHz) is consistent with PSS-FREERUN's answer being correct and the
coarser references carrying their own resolution bias, not a single clean
apples-to-apples proof at a fixed, practical reference resolution -- the
2ns comparison is reported as "info" in `tests/test_pss.py`
(`check_r5_freerun()`), not gated at 0.1%, for exactly this reason. NOT
independently re-measured: the 20-phase orbit-shape check the card's own
validate section asks for (R5 orbit vs the settled `.tran` at 20 phases)
-- only the period and the from-x0 stability cross-check were run this
pass; see the card #489 report's "what was not done".

### `.pac`/`.pstb` on a free-running orbit (card #489 SECOND pass)

**Derivation.** `pacshoot.c`'s existing LPTV linearisation assumes a FIXED
period T; a free-running orbit's small-signal response needs a period-
VARIATION term, because perturbing e.g. the hysteresis window's own
reference shifts not just the state trajectory but the trigger crossing
TIME itself. Project ruling #447 A (smooth switch models only, same
restriction as everywhere else in this file) turns out to simplify this a
lot: with NO true event/jump, `f+ = f-` at the crossing, so the saltation
term in SIMPLIS's own Stage 2 event correction (`research/simplis-pop-
ac.md` Sec.6 Stage 2, second bullet) is identically ZERO -- the only
correction needed is the SAME terminal projection Stage 1 already uses for
the large-signal problem, applied here to the complex boundary condition
instead of to the real monodromy.

Setting: `x*(t)` the converged orbit, trigger guard `g(x)=x_trig-trigval`
(no explicit t-dependence for a plain threshold), period T. A forcing
`b e^{jwt}` produces a state response `y(t)` (pacshoot.c's own per-step
propagated variable) satisfying, to FIRST order, that the perturbed
trajectory's own crossing happens at `t*=T+dT` where (Taylor-expanding the
guard around the nominal crossing, `x*_trig(T)=trigval` by construction):

```
trigval + f_trig(T)*dT + y_trig(T) = trigval   =>   dT = -y_trig(T)/f_trig(T)
```

(`f_trig(T)` is the trigger STATE's own slope on the nominal orbit at the
crossing -- a real, non-infinitesimal number, since the physical dynamics
there are NOT small even though y is). The state re-anchored to the ACTUAL
(moved) crossing, in the SAME co-moving-frame sense Stage 1's projection
uses, is `y(T) + f(T)*dT`:

```
y0_eff = y(T) - f(T)*y_trig(T)/f_trig(T) = Proj . y(T),
Proj = I - f(T) e_trig^T / (e_trig . f(T))
```

exactly SIMPLIS's Stage 1 terminal-projection MATRIX, reused here on
pacshoot.c's complex `y(T) = Phi*xi + z(w)` instead of on the real
monodromy `Psi(T)`. The `e^{jwT}` phase factor in the periodicity boundary
condition is UNCHANGED to first order (correcting it too would need
`dT` multiplied by an already-first-order quantity, i.e. second order --
negligible for linear AC): the corrected boundary-value solve is simply

```
(e^{jwT} I - Proj.Phi) xi = Proj . z(w)
```

in place of the driven engine's `(e^{jwT} I - Phi) xi = z(w)`; the
re-propagation/demodulation pass (`pacshoot.c`'s "pass 2") is UNCHANGED
(same second-order argument -- see `pacshoot.c`'s own comment at the
correction site). `Proj` is a rank-1 update (only ROW `trig_idx` of
Phi/z is ever subtracted, weighted by `f(T)[d]/f(T)[trig_idx]` per row
`d`) -- never built as a dense matrix.

**Implementation** (`pacshoot.c`'s `pac_engine_run()`/`pac_pstb_run_ex()`,
`pssshoot.c`'s `PSSshootFreerun()`): `f(T)` (all n components) is a cheap
backward-difference estimate, `(x_T - x_prev)/(T_meas - t_prev)`, using
`pss_run_until_trigger()`'s OWN pre-crossing sample (no new run -- that
function already has to locate this point to refine the crossing time;
`x_prev_out`/`t_prev_out` were added as optional outputs). `e_trig` is
resolved by matching `trigger=`'s raw node number against a native
capacitor's own positive terminal among pacshoot.c's reduced-state map
(the paradigm case for every hysteretic/COT/QR reference this project
has: the trigger IS the latch/comparator capacitor itself, e.g. R5's `q`)
-- any other case (an inductor-current trigger, or an algebraic node not
directly a reduced state) refuses cleanly (`E_UNSUPP`, clear message)
rather than running uncorrected. `trigger=`/`trigval=`/`trigdir=` needed
adding to `.pac`/`.pstb`'s own dot-card parser (`inp2dot.c`'s `dot_pac`/
`dot_pstb`) -- `.pss` already had this (used only via the CONTROL-MODE
`pss` command in every existing R5 fixture, which bypasses numparam
entirely), but no fixture had ever exercised `trigger=` on a TOP-LEVEL dot
card before this pass; hit the SAME numparam bareword-substitution bug
`.stb`'s `probe=`/`mode=`/`localgnd=` already has a documented, quoted-
value workaround for (card #466, `inp_stb_autoquote()`) -- `trigger="q"`/
`trigdir="rise"` (quoted, matching this project's `probe="Vname"`
convention) works; unquoted does not, on a dot card (the CONTROL-MODE
`pss` command is unaffected either way, per that command's own numparam
bypass).

**Validation, and an honest limit found.** New fixtures
`tests/fixtures/pac_r5_ctrl2out.cir` (a `Vrefsrc` DC=vref/AC=1 source
broken out of the hysteresis threshold expression, giving R5 -- which has
no separate compensator loop to break, unlike R3/R4's `XP lgprobe` splice
-- somewhere for a "control" perturbation to enter) and
`pac_r5_line2out.cir` (`Vin ... AC 1`, no splice needed). Compared against
a NEW, independent reference: a long (1.2ms), Hann-windowed, single-
frequency DFT on a plain cold-start `.tran` (baseline vs. perturbed,
`(Y_pert-Y_base)/U_pert` -- deliberately NOT a perturbed-minus-baseline
subtraction of two long free trajectories, since their switching instants
drift apart over a long window; the Hann window plus a long absolute
window (~500 switching periods) suppresses the incommensurate-switching-
content leakage into the excitation frequency's own DFT bin instead),
cross-checked by amplitude halving as the task asked:

| freq | reference (halving spread) | td-PAC, corrected | delta (corrected) | td-PAC, UNcorrected | delta (uncorrected) |
|---|---|---|---|---|---|
| 20kHz (4.7% fsw) | -0.0045..-0.0012dB / -0.36..-0.65deg (spread 0.003dB/0.29deg) | 0.063dB / -0.015deg | 0.07dB / 0.34deg | 0.067dB / -0.042deg | 0.07dB / 0.32deg |
| 60kHz (14% fsw) | 0.318..0.329dB / -0.49..-0.66deg (spread 0.01dB/0.17deg) | 0.627dB / 0.092deg | 0.30dB / 0.59deg | 0.647dB / -0.657deg | 0.32dB / -0.16deg |
| 150kHz (35% fsw) | 1.465dB / -5.78..-5.96deg (spread 0.0002dB/0.17deg) | 5.197dB / 3.820deg | **3.73dB / 9.60deg** | 4.513dB / -12.99deg | **3.05dB / -7.20deg** |

At 20-60kHz the correction is inside this project's usual 0.5dB/3deg
target (and consistently a little CLOSER to the reference than leaving it
off). At 150kHz BOTH the corrected and uncorrected engine disagree with
the reference by 3+dB/7-10deg -- the reference's own halving spread there
is 0.0002dB, so this is not reference noise. **Root cause not identified
this pass.** Two candidates, not distinguished: (1) the single-local-step
backward-difference `f(T)` estimate is imprecise at R5's OWN extremely
stiff switch-node time constant (`Csw`'s RC ~= 100pF/100S = 1ps -- two
orders of magnitude faster than the ~1ns comparator/latch constant this
file's other timing choices are tuned around, and exactly at the trigger
crossing `sw`'s own `tanh(ksw*(q-0.5))` is at ITS steepest point too,
since it shares the same argument zero-crossing as the trigger guard --
`PSS_FREERUN_FT_DEBUG=1` shows `f(T)`'s `Csw`-row component at ~4.5e10
V/s there); the correction's row weights `f(T)[d]/f(T)[trig_idx]` could be
amplifying that imprecision, worse at higher frequency. (2) `pacshoot.c`'s
own already-documented ("Accuracy" above) Jacobian-reuse floor, never
before characterised on a circuit this stiff (R3/R4 have nothing like
`Csw`'s 1ps constant) -- possibly just now exposed because `.pac` had
never run on R5 at all before this pass, corrected or not (the UNcorrected
column above is not close to the reference either).

**Because of this, `.pac`/`.pstb` on a free-running orbit STAYS REFUSED BY
DEFAULT** (same honest-failure-over-wrong-number principle as before this
pass) -- `PSSshootFreerun()` prints the clear stderr refusal and skips the
post-processing. The implementation above is real and wired up (not a
placeholder), reachable for further investigation via
`PSS_FREERUN_PAC_EXPERIMENTAL=1` (env var, explicitly NOT for production
use, printed as a warning on every use) -- see that env var's own comment
in `pssshoot.c`'s `PSSshootFreerun()` for exactly what is/isn't trusted.
`PSS_FREERUN_FT_DEBUG=1` additionally prints the `f(T)` estimate itself
(the `T_meas`/`t_prev` backward-difference bracket and all n components)
for anyone chasing the two candidate root causes above.

### card #493: root-cause investigation using a THIRD, independent reference

Card #493 revisited the 150kHz disagreement above using
`tools/pac_perturbed.py --freerun` (new this card -- a Hann-windowed,
single-bin DFT on independent long `.tran` runs, promoted from the SECOND
pass's own ad hoc, uncommitted reference script; see `docs/pac-perturbed.md`'s
"Free-running (autonomous) mode" section for the method and the full
numbers table). Fork work: worktree `ngspice-src-fr150`, branch
`aug-fr150`, based on the series head (`8cd652490`, patches 0001-0046).

**The independent reference confirmed the SECOND pass's own ad hoc
reference was right, not wrong** -- `--freerun` reproduces it closely at
20/60/150kHz (150kHz: 1.412dB/-5.892deg here vs 1.465dB/-5.78..-5.96deg
there). Unlike card #494's SIXTH pass (saltation), where the third
reference found the REFERENCE wrong and the C code right, here it is the
native engine that disagrees with two independent methods, not the other
way round. New this card: the disagreement does not start at 150kHz --
it is already outside the 0.5dB/3deg target by 100kHz (23.7% of R5's
422.68kHz fsw), well below where the driven engine's own near-fs/2
accuracy floor usually bites:

| freq (% of fsw) | pac_perturbed --freerun | native corrected | native uncorrected | delta corrected |
|---|---|---|---|---|
| 20kHz (4.7%) | -0.009dB/-0.09deg | 0.061dB/-0.01deg | 0.067dB/-0.04deg | +0.07dB/+0.08deg |
| 60kHz (14.2%) | 0.279dB/-0.45deg | 0.615dB/0.16deg | 0.647dB/-0.66deg | +0.34dB/+0.61deg |
| 100kHz (23.7%) | 0.769dB/-2.13deg | 1.875dB/0.92deg | 1.866dB/-3.11deg | +1.11dB/+3.05deg |
| 150kHz (35.5%) | 1.412dB/-5.89deg | 5.049dB/5.23deg | 4.513dB/-12.99deg | +3.64dB/+11.13deg |
| 200kHz (47.3%) | 2.329dB/-11.33deg | 11.588dB/41.95deg | 7.707dB/-47.73deg | +9.26dB/+53.28deg |

("native corrected" here is card #493's OWN improved `f(T)` estimate, see
below -- the old backward-difference estimate gives 5.197dB/3.820deg at
150kHz, i.e. no meaningfully different verdict.)

**Three candidates tried, all ruled out as the (sole) explanation**:

1. **`f(T)` precision** (the SECOND pass's leading suspect): replaced the
   backward difference (whatever width the outer probe run's own
   bracketing step happened to give -- measured 15.8ps on R5, already
   fairly fine, not the ~11.8ns outer Points-grid width originally
   suspected) with a controlled, fixed-width FORWARD difference from the
   converged crossing state `x_T` itself (`f+ = f-`, project ruling #447 A,
   so this side is exactly as valid), default 1e-13s
   (`PSS_FREERUN_FT_DT`, `PSS_FREERUN_FT_BACKWARD=1` keeps the old
   behaviour for comparison; `pssshoot.c`'s `PSSshootFreerun()`). Result
   at 150kHz barely moves (5.05dB/5.23deg vs the old 5.20dB/3.82deg) --
   both equally far from the 1.41dB/-5.89deg reference.
2. **`pacshoot.c`'s own Jacobian-reuse floor** (`PAC_EXACT_JACOBIAN=1`,
   card #465 THIRD pass, unrelated to free-running specifically): changes
   the 150kHz result by <1e-4 relative -- ruled out cleanly.
3. **Recording-pass grid resolution** (`pacpoints=4000` vs the fixture's
   default 200): changes the WRONG answer to a DIFFERENT wrong answer
   (5.22dB/-5.46deg, phase moves toward the reference but magnitude moves
   further away) -- consistent with this project's own prior experience
   with `pacpoints=` on the driven engine's near-fs/2 floor ("tried 200
   to 200000, no measurable change" above), not a fix.

**Also new this card**: both the corrected AND the UNCORRECTED
(`PAC_FREERUN_NO_PROJ=1`) native results grow wrong TOGETHER, tracking
each other closely and diverging from the reference starting at the same
~100kHz point -- the period-variation projection correction is not
uniquely at fault. This is consistent with a general per-step-
linearisation accuracy floor on R5's own extremely stiff switch node
(`Csw`, ~1ps local RC -- two to three orders of magnitude faster than
any state in R3/R4) that the driven engine's existing validation never
stress-tested, rather than a bug specific to the free-running correction
term. Not root-caused further this card (time budget).

**Disposition: documented and enforced, not fixed.** `.pac`/`.pstb` on a
free-running orbit remains REFUSED BY DEFAULT (`PSS_FREERUN_PAC_EXPERIMENTAL=1`
unchanged). New: a stderr warning (`pacshoot.c`'s per-frequency boundary-
solve loop, same style as patch 0046's saltation band guard) whenever an
experimental free-running `.pac`/`.pstb` frequency reaches 10% of the
MEASURED free-running switching frequency (there is no fixed `fs` here,
unlike the saltation guard's `fs/2`) -- deliberately more conservative
than the ~14% point actually measured clean, since only one fixture/
output was checked (same reasoning patch 0046 used for its own 0.90
threshold against a ~98.77%-clean measurement). Committed on
`ngspice-src-fr150`, branch `aug-fr150`: the improved `f(T)` estimate (a
genuine, if inconclusive, numerical improvement -- it rules out candidate
1 rather than just asserting it) and the new band guard.

### card #503: driven-clocked-replica experiment rules out Csw stiffness alone; a near-degenerate Floquet pair found; disposition unchanged

Card #503 revisited the card #493 disagreement (still open at 100-200kHz)
with three experiments the task specifically asked for: soften the switch
node, replace the hysteretic comparator by a clock at the same switching
instants and compare the (well-validated) DRIVEN engine, and sanity-check
the sideband-0 demodulation. Build: `ngspice-src-fr2`/`aug-fr2` (pin + 53
patches, `aaba50ba3f49-fr2-N`). Numeric baseline reproduced BIT-EXACT
against card #493's own table at 150kHz/200kHz before any change (confirms
the build; see the report table below) -- the new C change below only adds
a diagnostic `fprintf`, no numeric-path code changed, so every other card
#493 number is unchanged and re-used here rather than re-run.

**Experiment 1 (decisive): a DRIVEN clocked replica of R5's own switch
node.** `refs/r5_hysteretic_buck.cir`'s power stage (`Bhs`/`Bls`/`Csw`
100pF||100S = the SAME ~1ps local RC/`L1`/`Rdcr`/`Resr`/`Cout`/`Rload`)
unchanged; the hysteretic self-referencing latch (`Cq`/`Bq`) REMOVED and
node `q` driven directly by an ideal clock (`PULSE(0 1 0 4n 4n 0.65061e-6
2.36585204e-6)`, period/duty and 0-to-100% transition speed all matched to
the free-running orbit's own converged values measured this card:
T=2.36585204e-6s from `PSS-FREERUN: CONVERGED`, on-time 0.65461us from the
emitted orbit's own `v(q)` 0.5-crossings, transition slope ~2.5e8 V/s from
the same trace, matching the docs' own `f(T)[q]` estimate). This makes
period/duty EXACTLY known (driven, no trigger=/Newton-over-T needed) while
keeping the identical stiff switch-node dynamics. Native (ordinary, no
trigger=) `.pac solver=td` line-to-output (`Vin AC 1`) vs
`tools/pac_perturbed.py` (plain driven backend, NOT `--freerun`, `--fs`
now exact) at the same 20-200kHz range:

| f (Hz) | native driven | pac_perturbed (driven) | delta |
|---|---|---|---|
| 20000 (4.9% "fsw") | -34.5281dB/-169.9101deg | -34.6325dB/-169.8706deg (f_actual=20127.7) | +0.104dB/-0.040deg |
| 65000 (15.4%) | -54.9088dB/-156.9352deg | -54.9052dB/-156.9405deg (65027.8) | -0.004dB/+0.005deg |
| 110000 (26.0%) | -63.0556dB/-144.8497deg | -63.0920dB/-144.7846deg (110265) | +0.036dB/-0.065deg |
| 155000 (36.7%) | -67.8210dB/-135.4031deg | -67.7083dB/-135.6420deg (153702) | -0.113dB/+0.239deg |
| 200000 (47.3%) | -71.0348dB/-128.2369deg | -70.9660dB/-128.3837deg (198909) | -0.069dB/+0.147deg |

Every point inside 0.11dB/0.24deg -- an order of magnitude inside the
project's usual 0.5dB/3deg target, all the way to 47% of the free-running
`fsw`, the SAME frequency where the free-running engine itself is off by
9.26dB/53deg (see the card #493 table). This clocked replica's own driven
monodromy is `spectral_radius=0.9929` (`PSS-SHOOT: CONVERGED ... stable`)
-- MORE marginal (closer to the unit circle) than the real free-running
orbit's own projected monodromy (`spectral_radius=0.0507`), i.e. a HARDER
case by the metric this project already uses to flag near-unity-multiplier
trouble (see the R4 D~=0.45 section above) -- and the driven engine still
handled it cleanly. **Conclusion: R5's ~1ps Csw stiffness, in isolation, is
NOT sufficient to explain the free-running engine's breakdown.** Caveat:
this replica has `n=3` reduced states (no `Cq` latch state at all), so it
tests Csw's stiffness alone, not the Csw+Cq COMBINED fast pair discussed
next -- a fully faithful "clocked-but-with-latch" replica (PWL-replicated
`v(vout)` feeding `Bq`'s own tanh argument, breaking just the feedback
loop while keeping `Cq`'s 1ns relaxation dynamics intact) was designed but
NOT built this pass (time budget) -- see "what was not done".

**Experiment 2: tightening PSSshootFreerun()'s own Newton tolerance.**
`.pac` card field 7 (`Steady_coeff`)/field 6 (`SC_iter`) on
`pac_r5_ctrl2out_150k.cir`/`_200k.cir`: default (1e-4, 45 iters) converges
at iter=41, residual 3.05e-5, `spectral_radius=0.0507 stable` (matches
card #489/#493 exactly). Tightened to 1e-9/80 iters: **FAILS to converge**
at all, plateauing at residual 5.17e-5 regardless of iteration budget -- a
genuine numerical floor in the free-running Newton loop on this circuit,
not simply "needs more iterations". Tightened to 1e-5/60 iters: reports
CONVERGED (iter=42, residual 7.05e-6, tighter than default) but
`spectral_radius=1.125196 UNSTABLE` -- directly contradicting the
well-established stable, well-damped orbit (0.0507 at default tol, cross-
checked against a settled `.tran` in card #489). This is a red flag that
the finite-difference monodromy (`eps=1e-6` relative FD step, `pssshoot.c`
`PSSshootFreerun()`) is noise-dominated at this precision, most plausibly
in whichever column probes the `Csw`/`Cq` fast pair: `h = eps*|x0[j]|` is
only ~1e-6-1e-7 (V) for these states, and the resulting endpoint-state
perturbation needs to survive `pss_run_until_trigger()`'s own
regula-falsi/bisection crossing relocation (with an absolute crossing
tolerance around `1e-12*(1+trigval)` on the TRIGGER STATE, but no
equivalent floor on the OTHER n-1 states' own endpoint precision) through
a ~1ps-scale switch transition -- plausibly below the transient
integrator's own achievable precision there.

**Experiment 3 (no code change needed -- re-reading the ALREADY-PRINTED
Floquet report at the default, working tolerance): a near-degenerate
trivial eigenvalue.** `PSS-FREERUN: FLOQUET trivial(... index 2)=
-7.92329e-07-9.77806e-06i` and `multiplier[3]=-7.92329e-07+9.77806e-06i`
are EXACT complex conjugates, `|.|=9.81011e-06` for BOTH -- not a single,
well-separated REAL near-zero trivial eigenvalue (the theory's prediction:
ONE trivial direction, the real flow-direction eigenvector `x'(t)`, whose
eigenvalue must itself be real). The "smallest-modulus = trivial" pick
(`PSSshootFreerun()`'s own Floquet-report code) has, on R5, landed on one
member of a near-degenerate PAIR arbitrarily -- consistent with R5 having
TWO independent fast/stiff directions (`Csw`~1ps, `Cq`~1ns, both several
orders of magnitude faster than the ~2.365us period) whose Floquet
multipliers both collapse close to zero, not well separated from each
other the way a single-fast-state system's would be. Reproduced
identically across repeated runs (bit-exact). **This does not directly
touch `pacshoot.c`'s trigger-node-based projection** (`pac_find_trig_idx()`
resolves `trig_idx` from the `trigger=` NODE NAME, never from this
eigenvalue pick) -- it is corroborating evidence, not a mechanism, for
Experiment 2's own noise-floor finding.

**C change this card**: `pssshoot.c`'s `PSSshootFreerun()` Floquet-report
block only -- a diagnostic warning (`PSS-FREERUN: WARNING -- the "trivial"
... eigenvalue ... is NOT well separated from the next-smallest genuine
Floquet multiplier ...`) whenever the two smallest-magnitude eigenvalues
are within 5x of each other (R5: identical to 6 significant figures, far
inside this threshold). No numeric-path code changed -- every existing
`.pac`/`.pstb`/Floquet NUMBER this file or `pacshoot.c` produces is
byte-for-byte unchanged (confirmed: 150kHz/200kHz native-corrected/
uncorrected values reproduce card #493's table exactly on this build, see
the report). `pacshoot.c` itself is UNCHANGED this card -- the new
evidence points upstream, into `PSSshootFreerun()`'s own monodromy/Newton
conditioning, not into `pacshoot.c`'s period-variation formula (already
independently ruled out as the sole cause by card #493's own
corrected-vs-uncorrected-track-together finding, reaffirmed here).

**Root cause: a structural numerical-conditioning limit, not a single
targeted coding bug.** Taken together -- (a) Csw's stiffness ALONE,
isolated on a driven/no-trigger deck (even a harder, more-marginal one by
spectral radius), is handled cleanly to <0.11dB/0.24deg up to 47% of fsw;
(b) the free-running Newton solve's own finite-difference monodromy sits
at or past its numerical noise floor well before the tolerance the default
fixtures actually use (fails outright at 1e-9, reports a false
"UNSTABLE" at 1e-5); (c) a genuinely near-degenerate pair of near-zero
Floquet multipliers, consistent with TWO widely-separated fast time
constants both far faster than the switching period -- the evidence points
to the shooting-Newton + finite-difference-monodromy + backward/forward-
difference-`f(T)` approach itself being numerically fragile specifically
on circuits with MORE THAN ONE time constant several orders of magnitude
faster than the switching period, not a single fixable bug in the
period-variation projection term. R5 is, so far, the only fixture in this
project with two such disparate fast states (R3/R4 have nothing
comparable) -- whether a single very-fast state alone (no latch) would
also break free-running `.pac` is untested, since R5's own hysteretic
control fundamentally needs the latch state to define the trigger= node
at all; no free-running fixture without a latch-like second fast state
exists to check this independently.

**Disposition: still documented and enforced, not fixed.** `.pac`/`.pstb`
on a free-running orbit remains REFUSED BY DEFAULT
(`PSS_FREERUN_PAC_EXPERIMENTAL=1` unchanged, same 10%-of-measured-fsw band
guard from card #493). The new diagnostic makes a real, previously-silent
numerical red flag loud (same "honest failure over wrong number"
convention as this file's other warnings) but does not itself change any
number or close the gap -- not attempted as a "fix" this pass given the
structural nature of the finding and the remaining time budget.

### card #503 SECOND pass: the autonomous formulation, not the latch, is the dominant cause -- narrowed further, still not fixed

The coordinator did not accept "structural" from the first pass and asked
for four things, in order: (1) a clocked replica that keeps the LATCH
(`Cq`/`Bq`, not just `Csw`) so the first pass's replica (which removed the
latch) does not get credit for isolating the autonomous part; (2) a check
of whether finer steps through the crossing window closes the gap (the
first pass only checked general Newton-tolerance tightening, not targeted
step refinement); (3) a non-transient reference, since both existing
references are transient DFTs; (4) a fix if a cause is found, targeting
0.5dB/3deg. All four were attempted.

**Task 1: clocked-WITH-latch replica.** `Cq`/`Bq` kept EXACTLY as in
`refs/r5_hysteretic_buck.cir` (same 1ns relaxation, same self-referencing
`tanh` comparator); ONLY `Bq`'s `V(vout)` argument is replaced by a FIXED
PWL source (`voutrepl`, one period of the real converged orbit's own
`v(vout)`, tiled since this build's PWL parser does not support the `r=0`
repeat flag) -- this breaks the vout->Bq feedback path (the one closed
loop that makes the system autonomous) while leaving the latch's exact
nonlinear form/time-constant, and the whole circuit's switching instants
in place. The result is a genuinely DRIVEN system (period T known exactly
= the recorded orbit's own T=2.36585204e-6s) -- `.pss method=shoot` (no
`trigger=`) converges in 4 iterations, `spectral_radius=0.9929` (the SAME
marginal value the first pass's latch-FREE replica also reported -- q's
own dynamics are now open-loop/decoupled from the L-C-Resr subsystem's
own monodromy, a block-triangular structure, so this is a consistency
check, not a coincidence). **Caveat found while building this**: removing
the vout feedback removes the dominant NEGATIVE feedback too, so the
open-loop gain is much larger (~40dB at 20kHz, not ~0dB) -- this changes
the magnitude scale from the closed-loop measurement but does not affect
the comparison to `pac_perturbed.py` on the SAME (also open-loop) deck,
which is the actual test. `pac_perturbed.py`'s own amplitude had to be
cut from the project's usual 0.001 to 2e-6 to stay in the linear regime
at this gain (0.001 gave a 5.8dB amplitude-halving spread -- badly
nonlinear; 2e-6 gives 0.006-0.037dB, trustworthy):

| f (Hz) | native driven (Points=200 default) | pac_perturbed (amplitude=2e-6) | delta |
|---|---|---|---|
| 20000 | 40.328dB/-169.998deg | 38.723dB/-170.001deg (f=20127.7) | +1.604dB/+0.003deg |
| 65000 | 19.915dB/-157.001deg | 18.317dB/-157.551deg (f=65027.8) | +1.598dB/+0.550deg |
| 110000 | 11.771dB/-144.560deg | 10.206dB/-145.110deg (f=110264.6) | +1.565dB/+0.550deg |
| 155000 | 7.078dB/-134.785deg | 5.559dB/-136.621deg (f=153702.2) | +1.519dB/+1.836deg |
| 200000 | 3.984dB/-127.670deg | 2.385dB/-129.343deg (f=198908.7) | +1.599dB/+1.674deg |

**Not exact** (per the coordinator's own framing) -- but a roughly
CONSTANT ~1.5-1.6dB/<2deg gap, not the free-running engine's 9.26dB/53deg
at the same 200kHz point, and critically not GROWING/diverging the way
the free-running result does. With `pacpoints=200000` (finer recording-
pass grid, task 2's own lever, applied here to the DRIVEN+latch deck):

| f (Hz) | native driven (pacpoints=200000) | pac_perturbed | delta |
|---|---|---|---|
| 20000 | 38.517dB/-170.007deg | 38.723dB/-170.001deg | -0.206dB/-0.006deg |
| 65000 | 18.137dB/-157.249deg | 18.317dB/-157.551deg | -0.180dB/+0.302deg |
| 110000 | 9.990dB/-145.381deg | 10.206dB/-145.110deg | -0.216dB/-0.271deg |
| 155000 | 5.224dB/-136.152deg | 5.559dB/-136.621deg | -0.335dB/+0.469deg |
| 200000 | 2.010dB/-129.204deg | 2.385dB/-129.343deg | -0.375dB/+0.139deg |

**Every point now inside the project's usual 0.5dB/3deg target, up to 47%
of fsw, with the EXACT SAME latch dynamics as the real R5 orbit.** This
answers task 1's own framing directly: with the period fixed externally
(driven) and enough recording-pass resolution, the latch's steepness is
NOT a barrier -- the residual it leaves (the first, coarse-grid table
above) is a genuine but SECONDARY, grid-fixable discretization effect, not
the dominant cause. The dominant cause is specifically what changes
between "driven, T known exactly" and "free-running, T found by
trigger-relocation Newton" -- the autonomous formulation itself.

**Task 2: finer steps through the crossing window.** Two things were
tried. First, the SAME `pacpoints=200000` lever that just cleanly closed
the driven+latch replica's gap (table above), applied instead to the REAL
free-running fixtures (`pac_r5_ctrl2out_150k.cir`/`_200k.cir`,
`PSS_FREERUN_PAC_EXPERIMENTAL=1`): 150kHz goes 5.049dB/5.234deg (Points=
200 default) -> 5.216dB/-5.459deg (pacpoints=4000, card #493) ->
5.343dB/-0.475deg (pacpoints=200000) -- NOT monotonic, NOT converging
toward the reference's 1.412dB/-5.89deg (the phase actually swings past
the reference's sign and back). 200kHz: 11.588dB (default) -> 10.937dB
(pacpoints=200000) vs reference 2.329dB -- still 8.6dB/44deg off after a
1000x finer recording grid. **The SAME lever that cleanly fixed the
driven+latch replica does nothing systematic here** -- direct evidence
the free-running gap is not a recording-pass discretization floor, even
though discretization DOES explain the driven+latch replica's own
(much smaller) residual. Second, `PSS_FREERUN_FD_EPS` (new env var, this
pass, `pssshoot.c`'s `PSSshootFreerun()`) overrides the Newton monodromy's
relative FD step (default 1e-6, unchanged unless set) to test whether a
LARGER step (less FD noise, more truncation bias -- the opposite trade-off
from tightening `steady_coeff`, which the first pass already showed makes
things worse) helps: 1e-5 gives 150kHz 4.951dB/2.451deg, 200kHz
11.183dB/38.94deg (barely moved); 1e-4 FAILS outright (loses the trigger
crossing entirely, iteration 15); 1e-3 "converges" to a DIFFERENT orbit
(`spectral_radius=0.1555`, not the true 0.0507) giving 150kHz
5.393dB/-6.399deg, 200kHz 17.352dB/-22.35deg -- WORSE at 200kHz than the
default. Non-monotonic, sign-flipping, sometimes worse -- ruled out as a
usable lever, kept only as a gated diagnostic (default behaviour
unchanged; no test sets this env var).

**Task 3: a non-transient reference.** Two, since neither alone was fully
satisfying on its own:
1. **Internal consistency, no external reference needed.** The SAME
   `PSS-FREERUN` run that produces the 200kHz `.pac` answer ALSO prints
   this orbit's own large-signal Floquet report: `spectral_radius=0.0507`,
   every non-trivial multiplier real and small (-0.0507, 0.00257, plus the
   near-degenerate pair near 0 from the first pass) -- a heavily damped
   orbit with nothing resembling a high-Q or near-marginal pole pair. A
   linearised small-signal response of such an orbit cannot physically
   show sharp gain PEAKING (native: +11.6dB, ~3.8x) combined with LEADING
   phase (native: +42deg) at any frequency -- that combination needs a
   resonant pole pair close to the unit circle, which this SAME run's own
   stability analysis explicitly rules out. The native engine's own two
   outputs (large-signal Floquet, small-signal `.pac`) are mutually
   inconsistent on this point -- a self-contained discriminator, no
   external model needed, and the free-running engine's own instability
   (task 2) is a plausible mechanism for how an inconsistent answer could
   arise.
2. **A literature-informed order-of-magnitude bound.** R5 has no separate
   compensator -- the only dynamics between the hysteretic threshold and
   vout are the latch (`Cq`/`Bq`) and the L-C-Resr output filter, so the
   general shape of its small-signal response is bounded by the same
   "effective sampling" argument the sampled-data literature uses for
   current-mode/ripple-based control (Ridley's sampled-data poles for
   peak-current-mode control; Redl & Sun's ripple-based-control describing-
   function analysis) -- a zero-order-hold envelope at the switching rate
   (`H(f) = sinc(fT)*e^{-j*pi*f*T}`, T and T/2 bracketing the once- vs
   twice-per-period correction this deck's symmetric threshold shift could
   give) bounds `|H|` to roughly 0 to -3.5dB and the phase to a LAG of
   order `-40` to `-85` degrees at 200kHz (NOT derived to be exact --
   this project's own hysteretic loop reacts continuously, not on a fixed
   clock, so the true response should track TIGHTER than this envelope,
   which is exactly what the reference shows: -11.33deg lag, well inside
   the envelope's lag-only sign and comfortably smaller in magnitude).
   `pac_perturbed --freerun`'s own +2.329dB/-11.33deg sits inside this
   qualitative envelope (mild peaking is plausible for a bang-bang loop,
   nowhere near the envelope's own outer bound); native free-running's
   +11.588dB/+41.95deg sits far outside it in BOTH magnitude and phase
   SIGN (a lag-dominated envelope cannot produce a leading-phase point).

**Task 4: fix.** Not achieved for the free-running engine itself. Task 1's
result shows precisely where a fix would need to land -- the autonomous
orbit's own `x0`/`T` (and by extension the recording pass's `Phi`/`z`,
computed by integrating from that `x0` for that `T`) must be imprecise or
inconsistent in a way a driven deck's `x0`/`T` (found without any
trigger-relocation search) is not, since IDENTICAL code
(`pac_engine_run()`) driven by those two different-quality inputs gives a
grid-fixable ~1.5dB residual in one case and a grid-IMMUNE 9dB+ divergence
in the other. `PSS_FREERUN_FD_EPS` was the one concrete lever this pass
tried at the suspected mechanism (the trigger-crossing finite-difference
monodromy) and it made things worse or unstable, not better -- ruled out,
not a fix. No other candidate was narrowed precisely enough this pass to
implement and VERIFY a corrective change within budget; per this file's
own "honest failure over wrong number" convention, none is shipped.
`pacshoot.c` remains unchanged both passes -- all new evidence continues
to point at `pssshoot.c`'s `PSSshootFreerun()` (the orbit/trigger-
relocation machinery itself), not at the period-variation projection
formula.

**Disposition: still refused by default, cannot be ungated.**
`PSS_FREERUN_PAC_EXPERIMENTAL=1` unchanged. What is now known with much
higher confidence than the first pass: the dominant cause is the
autonomous/trigger-relocation formulation's own orbit precision, not
general stiffness and not the latch's steepness (both now demonstrated,
not just argued) -- ruled out as fixes this pass: recording-grid
refinement (task 2, works for driven+latch, not for free-running) and the
Newton FD-step `eps` (task 2, unstable/non-monotonic). The natural next
step (not attempted this pass, time budget): instrument the free-running
`x0`/`T`/`f(T)` precision directly against a from-scratch high-precision
reference (e.g. a very-long-settled `.tran` restarted from the converged
`x0`, many more digits than the 150k/200k comparison already in this
file), to quantify exactly how much precision is lost in
`pss_run_until_trigger()`'s own bisection/regula-falsi relocation and
whether feeding a deliberately higher-precision `x0`/`T`/`f(T)` (computed
OUTSIDE the main, necessarily-loose-for-robustness Newton loop) into
`pac_pstb_run_ex()` closes the gap.

### card #503 THIRD pass: hand-off precision confirmed as the dominant cause, fixed for 20-100kHz (<=0.4fsw); a smaller residual remains at 150-200kHz

Rebased `ngspice-src-fr2`/`aug-fr2` onto `aug` (pin + 55 patches) first,
clean rebase, no conflicts; build `aaba50ba3f49-fr2-N`, baseline
`aaba50ba3f49-s55`. Reproduced the SECOND pass's own numbers bit-exact
before any change (confirms the rebase changed nothing on this path).

**Task 1: an independent high-precision free-running orbit.** Seeded a
plain `.tran` (`.ic v(q)=0.5 v(cx)=3.300028428 v(sw)=5.947131838`, `L1
IC=0.7162297427` -- exactly `PSS_FREERUN_DEBUG=1`'s own converged x0, so no
cold-start transient is needed) at `reltol=1e-9`/`trtol=1`/`Tmax=20ps`,
measured 15 consecutive cycles by interpolated `v(q)=0.5` crossings.
Result: period and all four states are stable to display precision
(`numdgt=17`) from cycle 2 onward -- T=2.367108e-6s (vs PSS-FREERUN's own
2.365852e-6s, a **530ppm** difference), state at the crossing
[q, cx, sw, iL1] = [0.5, 3.300029, 5.947131, 0.7160859] (vs PSS-FREERUN's
[0.5, 3.300028, 5.947132, 0.7162297] -- q/cx/sw agree to <2e-7 relative,
but **iL1 differs by ~0.02%**). Re-run at `reltol=1e-11`/`Tmax=5ps`
reproduced the SAME numbers to displayed precision -- this reference is
genuinely converged, not limited by its own tolerance choice.

**Task 2: hand-off comparison.** The gap above (530ppm in T, 0.02% in
iL1, the other three states already agreeing to <2e-7) is what
`PSSshootFreerun()` hands to `pac_pstb_run_ex()` as `x0`/`T`.

**Task 3: inject the high-precision x0/T (decisive).** New env-var
override this pass, `PSS_FREERUN_X0_OVERRIDE`/`PSS_FREERUN_T_OVERRIDE`
(`pssshoot.c`'s `PSSshootFreerun()`, gated behind
`PSS_FREERUN_PAC_EXPERIMENTAL=1`, replaces the Newton loop's own x0/T for
the `.pac`/`.pstb` hand-off only -- the CONVERGED/FLOQUET report and the
emitted trajectory keep using the unmodified Newton x0/T). Feeding in
task 1's independent x0/T:

| f (Hz) | native (x0/T injected) | pac_perturbed --freerun | delta | native (unfixed, same build) | delta (unfixed) |
|---|---|---|---|---|---|
| 60000 | 0.0147dB/-1.479deg | 0.279dB/-0.45deg | -0.264dB/-1.029deg | 0.615dB/0.155deg | +0.336dB/+0.605deg |
| 100000 | 0.2140dB/-3.191deg | 0.769dB/-2.13deg | -0.555dB/-1.061deg | 1.875dB/0.920deg | +1.106dB/+3.050deg |
| 150000 | 0.6004dB/-5.460deg | 1.412dB/-5.89deg | -0.812dB/+0.430deg | 5.049dB/5.234deg | +3.637dB/+11.124deg |
| 200000 | 1.1331dB/-7.520deg | 2.329dB/-11.33deg | -1.196dB/+3.810deg | 11.588dB/41.950deg | +9.259dB/+53.280deg |

**Answers the decision rule directly: injecting the independent x0/T cuts
the error 5-8x at every frequency (e.g. 200kHz: 9.26dB/53.3deg ->
1.20dB/3.81deg) -- inside or very close to 0.5dB/3deg at every point,
decisively confirming hand-off precision as the DOMINANT cause**, though
not perfectly zero everywhere (see "residual" below). `pacpoints=200000`
stacked on top of the injected x0/T made 60-200kHz MUCH WORSE (200kHz:
13.7dB/-62.3deg) -- an unexplained adverse interaction between a
precision-corrected x0 and a finer recording grid, found but not
understood this pass (unlike the SECOND pass's driven+latch replica,
where finer grid cleanly helped); flagged, not chased further (time
budget) -- default `Points`/`pacpoints` is used throughout the fix below.

**Fix implemented and shipped.** Root cause: `PSSshootFreerun()` tightens
`ckt->CKTreltol` to AT MOST 1e-5 near the top of the function (`if (1e-5 <
save_reltol) ckt->CKTreltol = 1e-5;`, unchanged since card #489) -- a
robustness measure for circuits in general, but not tight enough to
resolve R5's own ~1e6:1 stiffness ratio (`Csw`~1ps RC vs the ~2.4us
period) to the precision `.pac`/`.pstb` needs at higher fraction-of-fsw.
Simply running the WHOLE shooting (Newton AND the recording pass) at a
tight reltol throughout (`.options reltol=1e-9` on the deck) was tried and
measured to help less, and sometimes hurt at high frequency -- probably
the recording pass's own already-documented per-step-Jacobian-reuse floor
(this file's "Accuracy" section) behaving differently at a very different
step-size regime is a confound. The shipped fix (`PSSshootFreerun()`, this
pass) instead POLISHES `x0`/`T` for the `.pac`/`.pstb` hand-off ONLY, with
a SCOPED reltol tightening that reverts before the recording pass runs:
after Newton converges (at its normal, robust tolerance), a few (default
5, `PSS_FREERUN_POLISH_ITERS`) extra PLAIN period-map evaluations
(successive substitution -- `pss_run_until_trigger()`, the SAME machinery
Newton's own probe steps already use, no new numerics) are run under a
temporarily tightened `ckt->CKTreltol` (default 1e-9,
`PSS_FREERUN_POLISH_RELTOL`) that also covers the immediately-following
`f(T)` re-estimation, then reltol is restored before
`pac_pstb_run_ex()`. This converges FAST: since the orbit is heavily
damped (`spectral_radius~=0.05`), the period map is a contraction at
roughly that rate, so 5 extra iterations already reach a fixed point
matching task 1's independent reference to 4-5 significant figures (T:
2.36585e-6 -> 2.36701e-6, a 530ppm -> ~40ppm improvement); 15 iterations
or a further-tightened `PSS_FREERUN_POLISH_RELTOL=1e-11` move the answer
by <1e-5 relative -- already converged, not iteration-starved. Extending
the tight-reltol scope to the `f(T)` estimation (this pass's own
refinement over an earlier version of this fix) changed nothing further
(<1e-4 relative) -- `f(T)`'s own precision was not the limiting factor
either.

**Numbers with the fix (default: polish on, 5 iterations, reltol=1e-9,
default Points/pacpoints), same build/deck/env as the "unfixed" column
above, `PSS_FREERUN_PAC_EXPERIMENTAL=1` still required:**

| f (Hz, %fsw) | native (fixed) | pac_perturbed --freerun | delta | delta (unfixed, for comparison) |
|---|---|---|---|---|
| 20000 (4.7%) | -0.005dB/-0.111deg | -0.009dB/-0.09deg | +0.004dB/-0.021deg | +0.070dB/+0.080deg |
| 60000 (14.2%) | 0.098dB/-1.416deg | 0.279dB/-0.45deg | -0.181dB/-0.966deg | +0.336dB/+0.605deg |
| 100000 (23.7%) | 0.334dB/-3.507deg | 0.769dB/-2.13deg | -0.435dB/-1.377deg | +1.106dB/+3.050deg |
| 150000 (35.5%) | 0.675dB/-6.430deg | 1.412dB/-5.89deg | -0.737dB/-0.540deg | +3.637dB/+11.124deg |
| 200000 (47.3%) | 1.076dB/-8.740deg | 2.329dB/-11.33deg | -1.253dB/+2.590deg | +9.259dB/+53.280deg |

**Target (0.5dB/3deg up to 0.4fsw=169072Hz): MET at 20kHz and 60kHz and
100kHz** (all three comfortably inside on both magnitude and phase).
**150kHz (still <=0.4fsw) is CLOSE but not fully inside** -- phase is
well inside (0.54deg of 3deg budget), magnitude misses by 0.24dB (0.74dB
vs the 0.5dB budget). 200kHz (outside the required 0.4fsw range) improved
by 7-8x in dB and ~20x in phase but is further from the target (0.75dB
over budget in magnitude, phase now inside at 2.59deg of 3deg).

**The remaining 150-200kHz residual is a genuine, separate, smaller
effect -- not further hand-off imprecision.** Confirmed by pushing
precision well past what the fix needs: `PSS_FREERUN_POLISH_RELTOL=1e-11`
and `PSS_FREERUN_POLISH_ITERS=15` both left the 150kHz/200kHz numbers
unchanged to 3-4 significant figures (150kHz delta stayed -0.74 to
-0.81dB across every precision level tried); extending tight reltol to
`f(T)` also changed nothing. Per the coordinator's own decision rule,
since hand-off precision demonstrably explains the DOMINANT part (5-8x
reduction) but a smaller, precision-independent residual remains at the
higher frequencies, the honest conclusion is a MIX: hand-off precision
was the dominant bug (now fixed, verified below), and a smaller,
genuine gap in how the free-running small-signal problem is posed likely
explains the rest. **The explicit bordered-system derivation (unknowns
y0 and dT, phase condition `e_trig . y(T) + f_trig(T) dT = 0` at the
trigger crossing) was NOT written out and compared symbolically against
`pacshoot.c`'s own rank-1-projection implementation this pass** (time
budget, after the numerical fix above consumed most of it) -- `pacshoot.c`
remains textually unchanged all three passes; docs/pss-pac.md's existing
derivation (card #489 SECOND pass section above) already states this same
phase condition and claims the rank-1 projection reproduces it exactly to
first order -- verifying that algebra symbolically, and checking the
"second-order, negligible" claim for the `e^{jwT}` phase factor
quantitatively at R5's own `dT` scale, is the natural next step if the
150-200kHz residual needs closing further.

**C changes**: `pssshoot.c`'s `PSSshootFreerun()` only, three additions
--(1) `PSS_FREERUN_X0_OVERRIDE`/`PSS_FREERUN_T_OVERRIDE` (test-only,
used to decide the hypothesis above, gated behind
`PSS_FREERUN_PAC_EXPERIMENTAL=1`); (2) the polish loop itself (default
ON whenever `PSS_FREERUN_PAC_EXPERIMENTAL=1` and no X0_OVERRIDE is given
-- `PSS_FREERUN_POLISH_ITERS=0` disables); (3) `PSS_FREERUN_POLISH_RELTOL`
(default 1e-9). Nothing outside the already-gated `.pac`/`.pstb` block
changed -- the CONVERGED/FLOQUET report, the emitted trajectory, and
every existing `.pss`/driven-`.pac` code path are byte-for-byte
unaffected (confirmed: `test_pss.py`'s R3/R4 determinism hashes,
unaffected by construction since they never reach this code path, plus
full regression below). `pacshoot.c` is unchanged all three passes.

**Disposition: PARTIAL ungate proposal.** `.pac`/`.pstb` on a
free-running orbit is still gated behind `PSS_FREERUN_PAC_EXPERIMENTAL=1`
by default -- NOT proposing removing the gate entirely, since 150kHz (a
point inside the required <=0.4fsw range) does not yet cleanly meet the
0.5dB target. Proposing instead: (a) keep the experimental gate, but (b)
note in this file and in the env var's own stderr message that the FIX
above (on by default whenever the gate is opted into) narrows the
validated-clean range to approximately 20-100kHz / <=0.24fsw on R5 (was
<=0.14fsw before this pass) -- a caller opting in for THAT range now gets
a materially more trustworthy answer, not yet a fully validated one
above it. The 10%-of-fsw warning band guard (card #493, `pacshoot.c`)
is UNCHANGED (still fires above 10% of measured fsw, i.e. still fires
inside the now-good 14-24% range too) -- not tightened or loosened this
pass, since it is pacshoot.c's own guard and this pass's fix lives
entirely in pssshoot.c.

**Tests** (NGSPICE=`aaba50ba3f49-fr2-N`, vs `aaba50ba3f49-s55`): `test_pac.py`
25/25 (the R5 free-running experimental info-check improved: 0.064dB/
0.49deg -> -0.002dB/0.39deg delta at 20kHz vs the Hann-DFT reference),
`test_pss.py` 40/40 (chains `test_stb_native.py`/`test_stb_native2.py`/
`test_stb_device.py`/`test_stb.py`, all clean), `test_pac_perturbed.py`
8/8 (including its own R5 `--freerun` 20/60kHz check, unaffected/still
passing), `test_checkpoint.py` 7/7, `test_stb_native.py` (standalone)
28/29 (1 pre-existing skip, unchanged) -- no regressions.

### card #503 FOURTH pass: ungated inside the validated band (project decision)

The THIRD pass's own numbers (above) showed the x0/T polish fix meets
this project's usual 0.5dB/3deg target against `tools/pac_perturbed.py
--freerun` at 20/60/100kHz (<=23.7% of R5's fsw) but not yet at 150kHz+
(35%+). Three options were weighed (keep gated / ungate in the
validated band and refuse outside it / ungate everywhere with a warning)
and decided: **UNGATE INSIDE THE VALIDATED BAND.**

**Behaviour, `pssshoot.c`'s `PSSshootFreerun()`, this pass:**

- `.pac`/`.pstb` on a free-running orbit now **runs BY DEFAULT** (no
  `PSS_FREERUN_PAC_EXPERIMENTAL` needed) whenever the requested sweep's
  own top frequency is **<= 0.25 x fsw** (`fsw = 1/T`, `T` the orbit's
  own POLISHED period -- the THIRD pass's fix, still on by default in
  this path). 0.25 is the chosen cutoff: comfortably above
  the 23.7% point actually validated clean, comfortably below the 35%
  point that was not (same margin-not-exact-boundary convention this
  file's other band guards already use -- patch 0046's saltation guard,
  card #493's now-superseded 10% guard).
- Above that band, the WHOLE analysis is **refused by default** (a
  clear stderr message naming the band, the measured `fsw`, and the
  escape hatch -- checked against the sweep's own top frequency,
  `job->PACfStop`; refusing the whole sweep rather than skipping
  individual points is the simpler of the two options the task allowed,
  chosen because a deck wanting some points inside and some outside the
  band can already split into two `.pac` cards, the same way this
  project's other sweep-range constraints already work) -- unless
  `PSS_FREERUN_PAC_EXPERIMENTAL=1` is set, in which case it runs anyway
  with a loud "ABOVE the validated ... band" / UNVALIDATED warning
  (same wording style as this file's other opt-in-but-loud env vars).
  `pacshoot.c`'s own existing per-frequency 10%-of-fsw warning (card
  #493) is UNCHANGED and still fires independently inside the now-good
  14-24% range too -- not tightened or loosened this pass.
- The x0/T polish (THIRD pass) stays on by default in this now-default
  path, unchanged.
- Never a silently-unvalidated number either way -- every code path
  prints exactly what it did and why.

**Validated-clean range, stated plainly: ~20-100kHz (<=23.7% of fsw) on
R5** -- the same range the THIRD pass's own numbers support, now the
actual default behaviour instead of an opt-in. 150kHz+ (>=35% of fsw)
needs `PSS_FREERUN_PAC_EXPERIMENTAL=1` and is not yet within target.

**Tests**: `tests/test_pac.py`'s R5 free-running checks were replaced
(the old "always refused by default" pair no longer matches reality) with
three checks matching the new behaviour: `check_r5_freerun_pac_default_band()`
(20/60/100kHz run WITHOUT the env var and meet 0.5dB/3deg against
`pac_perturbed --freerun`'s own reference numbers, hardcoded from the
THIRD pass table above -- worst 0.435dB/1.377deg, both inside target),
`check_r5_freerun_pac_band_refusal()` (150kHz WITHOUT the env var is
REFUSED, red-proof), `check_r5_freerun_pac_experimental_above_band()`
(150kHz WITH the env var runs and prints the ABOVE-band warning -- checks
the gate mechanism, not the still-unvalidated number at that point).

**C changes**: `pssshoot.c`'s `PSSshootFreerun()` only -- the outer gate
on the `.pac`/`.pstb` block changed from `&& getenv(...)` to always-enter,
with the band-check/refusal logic added inside (where the old "REFUSED by
default, no exceptions" `else if` branch used to be, now removed
entirely). `pacshoot.c` remains unchanged all four passes.

## Known limits

- **Smooth switch models only** (project ruling #447 A) -- no ideal `sw`, no
  XSPICE digital; not enforced by a refusal, just not validated or
  supported.
- **Driven (clocked) circuits**: `.pss method=shoot` (no `trigger=`).
  **Free-running (autonomous) circuits**: `.pss method=shoot trigger=...`
  (card #489, `PSSshootFreerun()`) -- see the section above; `.pac`/`.pstb`
  on a free-running orbit are IMPLEMENTED (card #489 SECOND pass, the
  period-variation projection; card #503 THIRD pass, an x0/T hand-off-
  precision fix) and, as of card #503 FOURTH pass, **run BY DEFAULT for
  frequencies <= 0.25 x fsw** (validated to 0.5dB/3deg against an
  independent reference, `tools/pac_perturbed.py --freerun`, up to ~23.7%
  of fsw on R5) -- refused by default above that band (a clear message,
  never a silently-unvalidated number), reachable there only with
  `PSS_FREERUN_PAC_EXPERIMENTAL=1` (root cause of the 150kHz+/35%+
  disagreement not fully identified; see the card #503 sections above).
- **No saltation/event handling** -- follows directly from "smooth models
  only" (true of both shooting engines).
- **`probe=`/`pstb_probe=` takes a plain top-level `DC 0` V-source only** --
  no `probe=<inst>:<terminal>` auto-splice (native `.stb` has this; `.pstb`
  does not, yet), no `lgprobe`-style multi-device probe subcircuits (see the
  trap above), no device-mode (gain-nulling) probing.
- **Above fs/2, only the scalar sideband-0 loop gain/transfer function is
  meaningful** -- this engine computes exactly one sideband (`k=0`); a real
  switching converter's response above half the switching frequency
  involves real coupling to the neighbouring sidebands (`k=+-1`, ..., not
  computed here) that a sideband-0-only measurement cannot see. Even within
  that scope, a deep-attenuation-correlated accuracy floor (see "Accuracy"
  above) affects at least the single grid point closest to fs/2 on every
  quantity this project checks, and several more points (down to ~50kHz on
  a 500kHz-switching converter) on line-to-output specifically, where
  closed-loop rejection pushes the true signal deep enough to reach it
  sooner -- not reduced by `pacpoints=`/`Points`, bounded and reported, not
  hidden, in every check that hits it.
- **`pnoise`/`psp`, sidebands != 0, `pac_events` (event-driven small-signal
  excitation) are all out of scope** -- see `pacshoot.c`'s own "what this
  does not do" header section. **`.pxf` (time-domain, `method="shoot"
  solver="td"`) is now implemented** (card #499, sideband 0 only, see its
  own section above) -- a batched-forward per-source sweep on top of the
  SAME `pac_engine_run()` `.pac`/`.pstb` already use, not the from-scratch
  backward-adjoint engine that section's "Syntax decisions" discusses and
  did not attempt.
- **Private fork only** -- none of `.pss method=shoot`, `.pac solver=td`, or
  `.pstb` exist in the system `/usr/local/bin/ngspice` (45.2); a build from
  `buildrecipe/` (or a prefix under `~/Work/apps/ngspice_augmentation/`)
  is required.

## Netlist-card quoting

`method="shoot"`/`solver="td"`/`probe="Vname"` (quoted) is this project's
own convention throughout its `.pac`/`.pstb` decks -- a real, since-fixed
upstream parser bug (`inp_casefix()`, patch 0014, card #477) used to corrupt
memory when a dot-card line's **last** character before the newline closed a
quoted argument; every example above keeps a token after the last quoted
field for exactly this historical reason, though the branch now includes the
fix and no longer strictly needs the workaround. `.pss`'s own
`method=shoot` (unquoted) is unaffected either way (never the last quoted
token on its line in this project's own decks).

## Saltation (event) corrections for ideal `sw`/`csw` and XSPICE (card #494, plan T6.1)

Second half of project ruling #447 A ("smooth models first, saltation for
sw/csw and XSPICE later"). Series `patches-wip/salt/` (0900-0901; branch
`aug-salt`, build `aaba50ba3f49-salt1`). Full derivation, exact numbers and
honest status: `patches-wip/salt/README.md` and the card #494 report. This
section is the theory summary + what shipped.

### Why smooth-model PAC/PSTB gets ideal switches wrong

Every `.pac`/`.pstb` linearisation above (both this file's `solver=td` and
the frequency-domain default) propagates a small-signal perturbation
through a SEQUENCE of per-step LINEAR solves, each using the accepted
step's own fixed topology (J_k = G(x_k)+ag0*C(x_k)). For a SMOOTH switch
(a B-source `tanh`/sigmoid conductance), G(x_k) is a genuine, differentiable
function of the node voltages, so a perturbation reaching the switch's own
control input changes G(x_k) itself, continuously, and the chain rule
carries that dependence through every step. For an IDEAL `sw`/`csw`
device, G(x_k) is a DISCONTINUOUS function of the control voltage (exactly
`SWonConduct` or exactly `SWoffConduct`, decided once per accepted step by
`swload.c`) -- an infinitesimal perturbation that does not flip which side
of the threshold the step lands on changes NOTHING in that step's J_k, so
the small-signal path from an ideal switch's OWN control port to its power
port is exactly zero in this per-step-linear scheme, at every step except
the ones where the switch actually changes state. The entire modulator
gain (∂t_switch/∂control -- the classic PWM small-signal gain Fm) lives
ONLY in those event steps. Measured this card (`tests/fixtures/
pac_r6o_ctrl2out.cir`, open-loop ideal-switch buck, Vctrl AC-injected
directly at the switch's own control node): baseline `.pac solver=td`
control-to-output gain is **exactly 0** at 1kHz/25.5kHz/50kHz -- not small,
identically zero, confirming the mechanism above.

### The correction: Fang's saltation matrix, mapped onto this file's recursion

At a switching event (surface h(x,t)=0, e.g. `v_ctrl - vt` for `sw`), any
small perturbation delta_x jumps:

```
delta_x(t_e+) = S * delta_x(t_e-)
S = I + (f+ - f-) grad(h)^T / (grad(h).f- + dh/dt)
```

(f+/f- = the vector field dx/dt just after/before the event; see
`research/simplis-pop-ac.md` S4/S6). For `sw`, grad(h) is nonzero only at
the switch's own two CONTROL nodes (an ideal switch's control port has no
electrical coupling to its power port, so h depends on nothing else); the
`grad(h).f-` term is folded into an empirically-measured `dh/dt` (the
control voltage's own finite-difference slope across the accepted step,
from `SWctrlvalue`'s existing `CKTstate0`/`CKTstate1` history -- reused,
not recomputed) rather than computed separately, since on ngspice's own
discretised grid that IS the total derivative along the trajectory.

Mapped onto this project's own discretise-then-linearise recursion (the
SAME per-step machinery `PSSsensRecordStep`/`PSSpacRecordStep` already use
for smooth devices' G(t)/C(t) history terms -- see those functions' own
header comments), `(f+-f-)` never needs to be built as an explicit
state-space vector: instead, one extra current-injection RHS term at the
switch's own POWER terminals, solved through the SAME per-step
factorization every other history term already uses, gives the network's
correct response directly:

```
coeff = (y_prev[posCntrl] - y_prev[negCntrl]) / h_dot
dI    = ag0 * dG * v_switch * coeff
rhs[posNode] -= dI;  rhs[negNode] += dI
```

`y_prev` is THIS column's own propagated sensitivity from the END of the
PREVIOUS accepted step (delta_x(t_e-), what the formula calls for -- not
the step currently being built). `dG` = conductance jump (g_after-g_before,
from the model's `SWonConduct`/`SWoffConduct`), `v_switch` = the switch's
own terminal voltage at the accepted point. `ag0` converts the resulting
charge-like quantity (dG*v_switch*coeff has units of Coulombs -- an Amp
times a "coeff" that works out to seconds) into the companion-model's
expected Amps for this RHS row -- the SAME role `ag0` plays in the
existing `ceq`/`veq` history terms just above it in both functions.
`PSSsaltFindEvents()` (`pssshoot.c`, shared via `pssshoot.h`) locates the
event each accepted step, post hoc, from `SWload.c`'s own
`CKTstate0[SWswitchstate] != CKTstate1[SWswitchstate]` (no new device
callback, no exact event-time location -- see "what remains" below).
`PSS_SALTATION=0` disables the whole thing (the red-proof toggle).

### Status (honest)

**Confirmed no-op on every existing fixture** (R1-R5, including R2's own
ideal `sw` -- its control node is driven only by a fixed clock, never by
any perturbation these tests inject, so the correction's own coefficient
is identically zero there): full regression on the `salt1` build --
`test_pac.py` 21/21 (3 info), `test_pss.py` 40/44 (4 info),
`test_checkpoint.py` 7/7, `test_sens.py` 27/27, and the `stb*` suites
(60/60, 28/29, 54/54, 16/16, 29/29, 42/42) all unchanged from the
pre-saltation baseline.

**Non-trivial but NOT yet quantitatively correct** on
`pac_r6o_ctrl2out.cir` (open-loop ideal-switch buck: `sw` compared
directly against a ramp via its own differential control pair, no
compensator loop, no XSPICE): with the patch, the same 1kHz sweep that
gave exactly 0 before now gives **13.19 V/V** (phase -1.9deg), the right
sign and right order of magnitude and falling off with frequency as
expected, vs an independent reference of **~5.8-6.0 V/V** (a
finite-difference of two converged `.pss method=shoot` orbits at
`vctrl_dc` = 0.55V and 0.551V -- see the card #494 report for the run) --
roughly **2.2x too high**, not resolved this pass. Finding and fixing the
missing `ag0` factor (without it the result was ~1e-12, i.e. still
effectively zero) was this pass's main progress; at least one more scaling
or double-counting issue remains (two switch instances, Shs and Sls,
change state at the same instant and both contribute a correction term --
summing both is believed physically correct but was not proven so this
pass; that is the leading suspect).

### What remains

- **CSW** (current-controlled switch) is not implemented -- `SW` only.
- **The closed-loop case** (a compensator/op-amp feedback loop driving the
  switch, i.e. R6/R6s as originally specified): `.pss method=shoot` itself
  was found NOT to converge on an ideal-`sw`-in-feedback-loop buck this
  card, independent of saltation (see `tests/fixtures/pss_r6s_shoot.cir`'s
  header for the full trail -- several convergence aids tried: softer
  comparator gain, a gate-drive RC filter, a direct differential-control
  comparator instead of a B-source tanh). A genuine, separate robustness
  gap in the shooting engine's forced-tight (reltol=1e-5) Newton/period-
  integration pass when combined with a hard conductance jump in a stiff
  RLC loop -- not chased further this pass.
- **XSPICE digital PWM** (`d_pwm`, R6 as originally specified in
  `refs/r6_vm_buck_ideal.cir`): `.pss method=shoot` SEGFAULTS in
  `cm_d_pwm` during the shooting engine's pre-stabilization transient.
  Root-caused via `gdb` (crash inside `cm_d_pwm+297`, a NULL/stale
  `cm_event_get_ptr(0,0)` dereference at TIME==0.0): XSPICE's MIF
  event-queue storage (`cm_event_alloc`/`cm_event_get_ptr`) is private
  state OUTSIDE `ckt->CKTstate0/1`, and the shooting engine's period-
  restart (re-entering a fresh transient at t=0 from a saved state, used
  for both pre-stabilization and every Newton trial period) does not
  know about it and cannot reset/preserve it consistently. `driven=1`
  (patch 0900) only fixes the SEPARATE "not recognised as driven"
  refusal; the segfault is unrelated and unfixed. Not chased further this
  pass -- fixing it would need either a new XSPICE-aware save/restore hook
  in the shooting engine, or accepting that XSPICE digital models are out
  of scope for shooting PSS entirely (a project-level scope decision, not
  made here).
- **Exact event-time location**: the event is detected AFTER the fact,
  from whichever accepted step `SWload.c`'s own state machine flipped on
  -- no `CKTsetBreak`-based precise crossing time (unlike `cm_d_pwm`'s own
  `cm_analog_set_temp_bkpt`, the idiomatic mechanism for this in the
  XSPICE side of this codebase). The error this introduces is O(h) in the
  step size, believed comparable to the existing discretise-then-linearise
  recursion's own O(h) treatment of a smooth device's G(t)/C(t), but not
  measured separately from the ~2.2x scaling gap above.
- **The 0.5dB/3deg-vs-FRA validation table** the brief asks for (R6
  control-to-output/line-to-output/output-impedance/.pstb loop gain from
  fs/1000 to fs/2) was NOT produced -- blocked on the closed-loop
  convergence gap above (no converged R6/R6s closed-loop orbit exists to
  sweep against FRA). Only the open-loop, single-frequency-point,
  approximately-DC comparison above was completed.
- `tests/test_pss.py`/`tests/test_pac.py` were NOT extended with a new,
  committed automated check for this feature this pass (the manual runs
  above were not turned into a `check_*()` function) -- the brief's
  "add checks; keep existing green" was only half done (existing suites
  confirmed green; no new check added).

### SECOND pass (card #494 second pass, plan T6.1): root-cause narrowing, still UNRESOLVED

Re-ran the exact `salt1` build (no C changes made this pass -- see "why no
fix" below) and reproduced the first pass's numbers exactly: `.pac
solver=td` on `pac_r6o_ctrl2out.cir` at 1kHz gives **13.1885 V/V, -1.888
deg** (13.19 V/V/-1.9deg as reported). Three independent references were
built and cross-checked against each other before trusting the "~6 V/V"
figure:

1. **Analytic averaged-model DC anchor**: Vin/Vramp = 12/2 = **6.00 V/V**
   (the textbook voltage-mode PWM modulator DC gain, valid at 1kHz since
   1kHz << the LC corner and << fs/2=250kHz).
2. **Finite difference of two converged `.pss method=shoot` orbits**
   (freshly re-run this pass, `vctrl_dc`=0.549/0.551, both CONVERGED
   cleanly: residual ~1e-7, spectral_radius=0.994, stable): trapezoidal
   cycle-average V(vout) gives **5.934 V/V** -- confirms the first pass's
   "~5.8-6.0" figure and agrees with the analytic anchor to ~1%.
3. **`tools/pssac_ref.py`'s perturbed-PSS method** (new `r6o_ctrl2out` TAG
   added this pass, additive): at 1kHz the exact-commensurate `.pss`
   backend needs q~498 (1kHz is 1/500th of fs=500kHz) and does NOT
   converge at that period count (`Timestep too small ... trouble with
   swmod-instance shs`) -- an LCM-period practicality limit, not a
   saltation bug. The tool's automatic long-`.tran` fallback then gives
   **6.39 V/V, -1.94 deg**, in the same ballpark as (1)/(2) but with a
   periodicity residual (8.2x/1.0x the window's own ripple) too large to
   trust rigorously as a third precise reference -- reported here as
   corroborating, not dispositive.

All three land at 5.9-6.4 V/V; the code's 13.19 V/V is a real, reproducible
~2.1-2.2x-too-high discrepancy, not a reference-agreement problem.

**Investigated and RULED OUT this pass** (with the `PSS_SALT_TRACE=1`
runtime trace, `pssshoot.c`'s `PSSsaltFindEvents()`):

- *Double-counting of the two complementary switch events* (Shs
  turning off, Sls turning on, at the same instant -- the leading suspect
  named at the end of the first pass): traced and found that BOTH events
  are indeed detected together (`n_events=2` at the same accepted step,
  `t~5.38e-7s`, one period's real duty edge). But `Shs`'s and `Sls`'s own
  `coeff = (y_prev[posCntrl]-y_prev[negCntrl])/h_dot` are **algebraically
  identical** (Sls's control pins are the exact swap of Shs's, on the same
  model/vt/vh, so both its numerator and `h_dot` flip sign together,
  cancelling) -- and the two current injections land at DIFFERENT node
  pairs (Shs: vin/sw; Sls: sw/gnd). Worked through by hand: the combined
  RHS contribution reduces to one clean term, `coeff * dG_shs * Vin`,
  injected at node `sw` -- a single physically-sensible combined
  correction (proportional to Vin, independent of the intermediate switch
  voltage), not two redundant applications of the same correction. This is
  NOT the 2.2x.
- *A missing sign / `dh/du` "input term"* (`research/simplis-pop-ac.md`
  S6 Stage 2's `-(f+-f-)(dh/du . b e^{jwt})/denom`, which `pss_r6o_shoot.cir`'s
  own header correctly flags as the relevant case -- Vctrl injects
  directly at the switch's own control node): worked through the reduced-
  state derivation and found the code's homogeneous `S*delta_x-` term,
  evaluated with `y_prev` at a DIRECTLY-DRIVEN node, algebraically
  reduces to exactly the same quantity the separate input-term formula
  calls for, but with what looked like the wrong sign on paper. Checked
  against the MEASURED output instead of trusting the derivation: the
  code's result is **+13.19 V/V** (same sign as the expected **+6 V/V**,
  small near-zero phase, physically the right shape for a buck's
  control-to-output response well below resonance) -- a pure sign error
  would flip this 180 deg, not scale it, so this is NOT the 2.2x either
  (the paper derivation above is unresolved/inconclusive -- flagged, not
  trusted over the measurement).

**Investigated, NOT resolved** (leading remaining suspect): the `ag0`
factor in `dI = ag0 * dG * v_switch * coeff`. A rough local approximation
(treating the affected node as capacitance-dominated, `J_k ~= ag0*C`)
shows `ag0`'s inclusion is dimensionally consistent with Fang's `(f+-f-)*
delta_t_e` state jump (the `ag0` in the numerator cancels the `ag0` inside
`J_k^-1`, leaving `dG*v_switch*coeff/C`, matching a rate times a time-shift
-- so simply DELETING `ag0` is not obviously correct either, and was not
tried as a blind guess). But the trace shows the real duty-cycle event is
detected at an extremely small accepted step (`delta=1.087e-12s`, giving
`ag0~9.2e11`) -- and AT THAT SAME NODE the switch's own instantaneous
conductance (~100 S, the `sw` model's `Ron`) is only ~9x smaller than
`ag0*Csw` (~920), i.e. comparable, not negligible -- so the simple
"`ag0` cancels to leave `X/C`" approximation used above is rough at best
at the exact step size these events actually land on, and a real scaling
sensitivity to that (incidental, LTE-driven) tiny bracketing step size
cannot be ruled out without a deeper, exact re-derivation matched to this
circuit's actual `J_k` (mixed L/C/switch-conductance at node `sw`) --
NOT completed this pass.

**Also found, separate, smaller bug, NOT fixed**: `PSSsaltFindEvents()`
fires a spurious event at `t=1e-10` (`n_events=1`, `shs` only) on EVERY
period-restart run (every Newton iteration, every PAC column) -- `h_dot`
there is anomalously huge (~5.5e9, ~6600x the real edge's), consistent
with `CKTstate1` reflecting a period-restart reset rather than genuine
end-of-previous-period history, not a real mid-period transition (Shs is
ON on both sides of the period boundary in the converged orbit). Its
`1/h_dot` scaling likely makes its contribution small relative to the
real edges, but it is a real correctness gap (events should be filtered
against the period boundary, or `CKTstate1` should carry over the
converged orbit's own wrap-around state) not chased further this pass.

**Why no code fix was attempted this pass**: the leading remaining
suspect (`ag0`/tiny-step scaling) needs a from-scratch, exact
re-derivation matched to this fixture's own mixed L/C/switch-conductance
node, not a plausible-looking one-line change -- and the double-counting
and sign hypotheses that WERE cheap to test both checked out negative
against the measurements above. Given the time available this pass, a
guessed edit that happened to move 13.19 closer to 6 without a real
derivation behind it would not have been trustworthy, so none was made;
the C source tree (`ngspice-src-salt`, branch `aug-salt`) is unchanged
from the first pass (`HEAD=7879490d8`), and `patches-wip/salt/0901-...` is
unchanged. New this pass (additive, uncommitted, for the parent to
review): the three-reference agreement above, and `tools/pssac_ref.py`'s
new `r6o_ctrl2out` TAG.

**Third pass should**: (1) do the exact re-derivation of what `ag0`
(or the correct step-independent replacement) should be at a node with
comparable G/ag0*C, rather than the rough capacitive-node approximation
here; (2) filter/fix the spurious `t~0` period-restart event; (3) only
then attempt the closed-loop (R6s), XSPICE (`d_pwm` segfault), red-proof
and full-regression tasks, all NOT reached this pass either.

### THIRD pass (card #494 third pass, plan T6.1): state-space jump implemented (no `ag0`), still 13.19 -- root cause moved to event-time location

This pass implemented the parent's suggested approach directly: compute
Fang's saltation correction as a genuine STATE-SPACE jump,
`delta_x(t_e+) = delta_x(t_e-) - (f+-f-)*dt_e`, `dt_e =
-(grad(h).delta_x(t_e-))/h_dot`, instead of routing an `ag0`-scaled
companion-model current through the event step's own `G+ag0*C` solve (the
first two passes' approach, and the leading suspect they left
unresolved).

**The `(f+-f-)` re-derivation (closed form, no `ag0`).** For an ideal `sw`
event, `grad(h)` is nonzero only at the switch's own two control nodes
(unchanged from the first two passes). The new part: `(f+-f-)` at the
switch's own POWER nodes is not routed through any matrix solve at all --
it is computed directly from KCL. At the event, the state itself (every
capacitor voltage, every inductor current) is CONTINUOUS -- only its time
derivative jumps. Every branch current touching a switch's power node
OTHER than the switch's own current injection therefore does not jump
either (an inductor's current cannot jump; a resistor's current is a
continuous function of two continuous node voltages) -- so by KCL, the
switch's ENTIRE current-injection jump `dI = dG*v_switch` (`dG` =
`g_after-g_before`, `v_switch` = the switch's own terminal voltage, held
fixed as the continuous state it is) must be absorbed by that node's own
directly-attached capacitance alone:

```
f+_node - f-_node = -+ dI / C_local(node)
```

`C_local(node)` = the sum of `CAPcapac` over every native capacitor
instance with a terminal at that node (`pss_salt_local_cap()`,
`pssshoot.c`; walks `CKThead[CAPtype]` the same way `pss_build_state_map()`
already does) -- 0, correctly, for a node with no capacitor directly on
it (an ideal voltage source's own node has no state to jump). This is a
genuine Farads-based physical quantity: no `ag0`, no dependence anywhere
on the step size the event happens to land on. Hand-verified against this
same pass-2 section's own "double-counting ruled out" manual algebra: for
`pac_r6o_ctrl2out.cir`'s complementary Shs/Sls pair, summing both
switches' `(f+-f-)` contributions at node `sw` reduces to `-dG_shs*Vin /
Csw` exactly -- matching the earlier "`coeff * dG_shs * Vin`, a single
physically-sensible combined correction" finding above, now derived from
first principles rather than observed algebraically.

**Where it is applied: after the step's own solve, not through it.** The
jump `dfdt_pos/neg * coeff` (`coeff = (y_prev[posCntrl]-y_prev[negCntrl])
/h_dot`, unchanged from the first two passes) is added DIRECTLY onto the
already-solved `g_sens_Snew`/`g_pac_Cnew` for this column, AFTER the
step's own smooth dceq/dveq history solve -- not injected into that
solve's RHS. An earlier version of this pass tried the more "obvious"
placement instead -- applying the jump to `g_sens_Sprev`/`g_pac_Cprev`
BEFORE the solve, letting the existing per-step recursion carry it
through the event step's own (tiny) `G+ag0*C` matrix -- and measured
**7.73 V/V**, closer to the ~6 V/V references than either this pass's
final number or the first two passes' 13.19. Investigating WHY before
accepting it: routing the jump through that matrix divides it by
`(ag0*C_local + G_new_switch + ...)` instead of leaving it as the pure
value Fang's formula calls for, and the just-flipped switch's own new
conductance is directly comparable to `ag0*C_local` at that node (the
same "not negligible" finding the second pass already flagged) -- so the
"before" placement's closeness to 6 V/V is a COINCIDENTAL, step-size-
dependent partial cancellation, not a correct result. Confirmed by
re-deriving the algebra: in the `ag0 -> infinity` (`h -> 0`) limit, the
OLD ag0-routed formula and this pass's closed-form `(f+-f-)` converge to
the exact same quantity (the `ag0*C/(ag0*C+G)` factor -> 1) -- meaning
the first two passes' approach was ALREADY, at leading order, computing
something close to this pass's own "after" formula, not a fundamentally
different (wrong) one; the "before" placement's 7.73 V/V is the outlier
that needed explaining, not the 13.19 V/V both other placements agree on.
The "after" placement is what shipped (`patches-wip/salt/0902-...patch`).

**Result: unchanged from the first two passes.** `pac_r6o_ctrl2out.cir`
at 1kHz: **13.1885 V/V, -1.888 deg** (13.19/-1.9, matching the first two
passes almost to the digit) -- despite arriving via a completely
different, from-scratch, `ag0`-free derivation. The three references
carried over from the second pass (6.00 analytic, 5.93 FD-of-converged-
orbits, 6.39 `pssac_ref.py`) were not independently re-run this pass
(only the trivial analytic Vin/Vramp=6.00 anchor needs no simulation and
was re-confirmed); still 2.1-2.2x below the measured 13.19.

**New finding this pass (the actual result, not a fix): `h_dot` is not
converging.** `PSSsaltFindEvents()`'s `h_dot` is a post-hoc finite
difference of the switch's own control voltage
(`(vctrl_now-vctrl_prev)/CKTdelta`, from `SWctrlvalue`'s `CKTstate0/1`
history) across the accepted step SWload.c's own hysteresis-driven
step-halving happened to land on. On this fixture's real duty-cycle edge,
`Vramp` is a PULSE source rising perfectly linearly (zero curvature) at
exactly `vramp/(0.97*tclk)` = 1.0309e6 V/s over the whole window in
question -- so an exact finite difference over ANY sub-interval of that
segment should reproduce 1.0309e6 V/s to full floating-point precision,
independent of window size. Measured instead, varying only `.options
reltol`:

| reltol | event step size | h_dot (V/s) | `.pac` gain (V/V) |
|---|---|---|---|
| 1e-3 | 1.087e-12 s | 825443 | 13.19 |
| 1e-4 | (same order) | -- | 13.38 |
| 1e-5 | 1.087e-12 s | 825443 | 13.38 |
| 1e-6 | 6.978e-13 s | 521443 | 18.02 |
| 1e-9 | -- | -- | PSS-SHOOT itself fails to converge (ngspice error 106 in the pre-stabilization transient) -- plausibly a generic reltol=1e-9 practicality limit unrelated to saltation, not chased further |

`h_dot` should be a CONSTANT and is not: it drifts from 825443 to 521443
(AWAY from the true 1.0309e6, not toward it) between reltol=1e-5 and
1e-6, and the final `.pac` gain drifts correspondingly (13.19 -> 18.02)
rather than settling to a fixed value. **This means the card's validation
step (1) -- "show the result no longer changes when you change
reltol/tmax" -- is NOT met**; the opposite is measured and reported here
honestly. Since this pass's `(f+-f-)` formula is itself provably free of
`ag0`/step-size dependence (a closed-form Farads-based quantity), the
remaining step-size sensitivity traces to something else that both this
pass and the first two share unmodified: the EVENT LOCATION AND `h_dot`
ESTIMATE ITSELF, both derived post hoc from whichever accepted step
`SWload.c`'s internal state-machine/step-halving logic happens to
converge to (not from an exact, breakpoint-located crossing time) --
exactly the gap the first two passes' README already flagged as "NOT
reached" but had not identified as the LIKELY dominant error source until
this pass's reltol sweep.

**Confirmed unaffected (full regression, identical to the first two
passes' baseline numbers)**: `test_pac.py` 21/21 (3 info, including ASan
under a new `aaba50ba3f49-salt3-asan` build built this pass since C
changed), `test_pss.py` 40/44 (4 info), `test_checkpoint.py` 7/7,
`test_sens.py` 27/27, `test_stb_native.py` 28/29, `test_stb.py` 60/60,
`test_stb_native2.py` 54/54, `test_stb_device.py` 16/16,
`test_stb_followups.py` 29/29, `test_stb_cli.py` 42/42, `test_compopt.py`
20/20. `PSS_SALTATION=0` still gives exactly 0 gain on
`pac_r6o_ctrl2out.cir` (red-proof).

**NOT attempted this pass, same gaps as passes 1/2**: CSW; the closed-loop
R6/R6s case (`.pss method=shoot` still does not converge there, not
investigated further); the XSPICE `d_pwm` segfault; the 0.5dB/3deg-vs-FRA
validation table (blocked on the closed-loop gap, and now also on this
pass's own event-location finding above -- a table built on an
event-location-sensitive number would not be trustworthy). The spurious
`t~0` restart event the second pass found IS fixed this pass
(`PSSsaltReset()`, a one-shot latch called from `pss_sens_reset()`/
`pac_reset()` that makes the next `PSSsaltFindEvents()` call return 0).

**Why the `ag0`-magnitude question is now considered answered, not still
open**: the second pass's own "leading remaining suspect" was `ag0`'s
scaling at the tiny event step. This pass shows the `ag0`-free reformulation
gives the SAME 13.19 V/V the `ag0`-routed formula gave -- and explains
why (the `ag0 -> infinity` limit argument above) -- so `ag0` itself was
never the bug; it was, at leading order, computing close to the right
thing all along. The real lever is event-time/`h_dot` precision.

**Fourth pass should**: implement exact event-time location for `sw`
crossings (a `CKTsetBreak`-based bisection to the true `h=0` instant,
analogous to XSPICE's own `cm_analog_set_temp_bkpt` idiom in this
codebase) instead of the post-hoc `CKTstate0/1` finite difference this
and the first two passes all share -- the reltol-sweep table above is
reasonably strong evidence this, not a further scaling guess at the jump
formula, is what closing the remaining ~2.2x gap and achieving genuine
step-size independence needs. Only after that should the closed-loop
(R6s), XSPICE, red-proof-table and full 0.5dB/3deg-vs-FRA validation
tasks be attempted -- an event-location fix could plausibly change the
converged R6 numbers enough to be worth doing first.

### FOURTH pass (card #494 fourth pass, plan T6.1): arithmetic check done FIRST, per the brief -- exact-`h_dot` extrapolation still ~2x off; NEW root cause found, NOT the one the third pass predicted; no C change made

Per the brief's explicit task order ("before coding, check the arithmetic"),
this pass did NOT touch `ngspice-src-salt` (no C change, `HEAD` still
`9a21d9eb5` after a clean `git rebase aug` onto the current pin+0001-0041,
no conflicts). Instead it built an independent, from-scratch Python model
of `pac_r6o_ctrl2out.cir` (no numpy/scipy available in this sandbox --
pure-python complex arithmetic + a hand-rolled dense 4x4 real-matrix
exponential via scaling-and-squaring/Taylor series, unconditionally stable
regardless of stiffness, unlike the first attempt with explicit RK4, which
blew up to NaN -- see below) to answer task item (1) directly:
`patches-wip/salt/check_r6o_arithmetic.py` (+ `linalg4.py`), reproduce with
`python3 patches-wip/salt/check_r6o_arithmetic.py`.

**Model**: reduced 4-state LTI system `x = [V(sw), V(nsnub), I(L1), V(cx)]`
(the fixture's only capacitor/inductor states; `swl`/`vout` are algebraic,
eliminated by hand), two fixed topologies (`Shs` ON/`Sls` OFF and vice
versa) switching once per period at the exact analytic event time
`t_e = vctrl_dc/ramp_rate = 5.335e-7 s`. Hand-derived KCL confirms the
THIRD pass's closed-form result exactly: summing both complementary
switches' `(f+-f-)` contributions at node `sw`, every `V(sw)`/`V(nsnub)`/
`I(L1)`-dependent term cancels algebraically between the two switches,
leaving `f+_sw - f-_sw = -(Gon-Goff)*Vin/Csw = -1.2e12 V/s` exactly,
independent of the operating point -- the jump FORMULA itself checks out
by independent derivation, confirming the third pass's own finding.

**Task (1) result -- the exact-`h_dot` arithmetic check**: propagating
that jump with the ANALYTICALLY EXACT ramp slope (`ramp_rate = 1.030928e6
V/s`, not a simulated/drifting finite difference) through an EXACT
(matrix-exponential) integration of the reduced model's own stiff
dynamics gives **6.0173 V/V, -1.887 deg at 1kHz** -- matching all three
independent references (6.00 analytic, 5.93 FD-of-orbits, 6.39
`pssac_ref.py`) to within 0.3-2%, and matching the ACTUAL CODE's reported
PHASE (-1.888/-1.89 deg, all three passes) to three significant figures.
This proves the saltation jump formula is correct, both algebraically and
numerically, when correctly propagated -- exactly as the third pass
suspected. **But** cross-checking with the SAME (drifted, buggy) `h_dot`
values the actual code measured (825443 and 521443 V/s, from the third
pass's own reltol table) -- still using this pass's exact, correctly-
relaxed propagation -- gives **7.52 V/V and 11.90 V/V**, NOT the code's
actually-reported 13.19 V/V and 18.02 V/V. **This is the decisive
result**: even granting the code its own (wrong) `h_dot`, a correctly-
propagated model does not reproduce what the C code outputs. There are
therefore (at least) TWO separate, compounding defects, not one -- so per
the brief's own conditional ("if it is still ~2x off... look for a
missing/extra factor"), implementing exact event-time location (task 2)
alone would NOT be expected to close the gap to the <0.5dB/3deg tolerance
task (3) asks for; it would plausibly move the code from 13.19 toward
something in the 7.5-12 V/V neighbourhood, still 1.25-2x off.

**The second, newly-identified defect**: node `V(sw)`'s own local
dynamics are governed by `Csw=1nF` against the switch's own ~100S ON
conductance -- a time constant of about `1e-11s`, roughly FIVE ORDERS OF
MAGNITUDE faster than everything else in the circuit (the L1/Cout corner
is ~5kHz) and faster than the `coeff`/`dt_e` scale (~1e-6s) the jump is
multiplied by. The jump itself is a legitimate SENSITIVITY (not a literal
voltage), so a huge raw value at that node is not wrong in itself -- but
`PSSsensRecordStep`/`PSSpacRecordStep` propagate the rest of the period
using TRAPEZOIDAL integration (both gate explicitly on
`CKTintegrateMethod != TRAPEZOIDAL`) at accepted step sizes chosen by LTE
control on the SLOW states, nowhere near fine enough to resolve a
`1e-11s` transient. Trapezoidal integration is A-stable but NOT L-stable:
for `|lambda*h| -> infinity` (exactly this node's regime), its
amplification factor `(1+lambda*h/2)/(1-lambda*h/2) -> -1`, not `0` --
instead of damping the huge post-jump spike to near zero within the next
accepted step (what the exact solution, or an L-stable method, does), it
can under-damp/ring, leaking an uncontrolled fraction of it into the
genuinely slow `I(L1)`/`V(cx)` states over the following steps. This is a
textbook trapezoidal-stiffness artifact, and a well-fitting explanation
for both the residual excess gain (13.19/7.52 = 1.75x, 18.02/11.90 =
1.51x, on top of the already-accounted-for `h_dot` drift) and for why the
result does not stabilize with reltol (a ringing artifact's residual is
controlled by essentially arbitrary step-size boundaries, not by a
converging quadrature).

**This also re-opens the third pass's "before vs after placement"
finding**: the third pass measured 7.73 V/V routing the jump through the
event step's OWN `G+ag0*C` solve ("before" placement) and REJECTED it as
"a coincidental, step-size-dependent partial cancellation" in favour of
the 13.19 V/V "after" (undiluted) placement that shipped, reasoning from
the `ag0 -> infinity` limit. This pass's independently-computed
drifted-`h_dot`-but-correctly-relaxed number, **7.52 V/V**, lands within
3% of that rejected 7.73 V/V figure -- using a completely different
method (an exact matrix exponential, no `ag0` or companion model
anywhere) with the SAME `h_dot`. That is a strong, independent
corroboration that the "before" placement's proximity to the ~6 V/V
references was NOT coincidental: it happens to supply (imperfectly) the
same physical damping this pass's exact model supplies deliberately. The
third pass's rejection reasoning examined the wrong limit (`ag0 ->
infinity`, i.e. `h -> 0`) -- the relevant physical criterion is whether
the accepted step is large compared to the node's own `~1e-11s` time
constant, which it overwhelmingly is; some real damping there is what
correctness requires, not an artifact to route around.

**RK4 blow-up, itself informative**: an initial attempt using explicit
RK4 (fine substeps, 4000/segment) produced NaN outright -- `|lambda*h| ~
1e11 * 1.3e-10 ~ 13`, far outside RK4's stability region. This is an
independent numerical confirmation that node `sw` is exactly as stiff as
the analysis above says; switching to an unconditionally-stable matrix
exponential (pure-python scaling-and-squaring, `patches-wip/salt/
linalg4.py`, no numpy/scipy in this sandbox) fixed it, and the 1kHz
result is stable to 6 significant figures across Taylor order (24 vs 30)
and substep count (200 vs 800) -- not a numerical-precision artifact.

**Recommendation for a fifth pass** (no C change attempted this pass,
per the time-box, once the *original* plan -- "just fix `h_dot`" -- was
shown insufficient by the arithmetic check): do NOT spend the pass on
`CKTsetBreak` bisection alone. Either (a) treat node `sw` as algebraic/
quasi-static in the "fuller reduced-Jacobian solve" the third pass's own
derivation comment already flagged as not implemented ("a circuit where a
capacitor bridges two simultaneously-switching nodes would need..." --
the same idea generalises to any node whose local RC time constant is
far below the step grid), eliminating the fast state from the jump
target entirely rather than carrying a raw value through it; or (b) force
one or more L-stable (backward-Euler) sub-steps immediately after a salt
event before resuming the normal trapezoidal grid, so the fast mode is
damped the way the true physics (and this pass's model) says it should
be. Exact event-time location (the third pass's own recommendation) is
still worth doing -- the `h_dot` drift is real and costs ~25-60% on its
own -- but should be paired with, not substituted for, a fix to this
stiffness/placement issue.

**Confirmed unchanged**: no C source touched, so no new build and no new
regression run this pass -- the existing `aaba50ba3f49-salt3`/
`salt3-asan` builds and the third pass's own regression numbers (test_pac
21/21+3info, test_pss 40/44+4info, test_checkpoint 7/7, test_sens 27/27,
all `stb*` suites, test_compopt 20/20) are unchanged and were not
re-run (nothing that would affect them changed).

**NOT reached this pass**: any C code change; CSW; the closed-loop R6/R6s
case; the XSPICE `d_pwm` segfault; the 0.5dB/3deg-vs-FRA validation table
(same gaps as passes 1-3, now also blocked on the fifth-pass fix above);
red-proof re-verification (unchanged from third pass, not re-run since no
C source changed).

### FIFTH pass (card #494 fifth pass, plan T6.1): both diagnosed bugs fixed -- 13.19 V/V -> 6.20 V/V

Implemented the fourth pass's own two-bug diagnosis (`patches-wip/salt/
0903-...patch`, HEAD `ebafb2db8` on `aug-salt`; build `aaba50ba3f49-
salt5`/`-salt5-asan`). Full derivation in `pssshoot.c`'s comment block
above `PSSsaltActive`'s definition and `patches-wip/salt/README.md`'s own
"FIFTH pass" section; summary here.

**Task 1 (exact event location/`h_dot`).** Root-caused (not guessed) by
reading `dctran.c`'s actual control flow: `ckt->CKTdelta`, by the time
`PSSsaltFindEvents()` runs (right after `CKTaccept()`), has ALREADY been
overwritten with `CKTtrunc()`'s suggestion for the NEXT trial step
(`ckt->CKTdelta = newdelta;`, dctran.c ~line 993) -- not the size of the
step that produced the `CKTstate0`/`CKTstate1` transition being
differenced. The exact same "read after accept but already advanced"
trap `PSSsensRecordStep()` already documents for `CKTorder` (card #458
fourth pass), never previously traced through for `CKTdelta`. Fix:
`PSSsaltFindEvents()` now tracks its own history of `ckt->CKTtime`
across calls and uses the TRUE elapsed time between consecutive
accepted points as `h_dot`'s denominator -- self-contained, no
dependence on `dctran.c`'s internal `CKTdeltaOld[]`/`CKTorder`
bookkeeping. Measured: `h_dot` now matches the exact ramp slope
(1.030928e6 V/s) to 6 significant figures, independent of reltol.

**Task 2 (L-stable propagation after an event).** New
`pss_salt_local_cond()` (mirrors `pss_salt_local_cap()`, sums native
Resistor `RESconduct` + Switch `SWcond` on a node) estimates the
triggering node's own `tau = C_local/G_local`. A settle window
(`PSS_SALT_SETTLE_TAUS=30` local time constants of real elapsed time)
opens on every real event; `PSSsaltForceOrder1()` overrides
`PSSsensRecordStep()`/`PSSpacRecordStep()`'s own order inference to 1
(backward Euler, L-stable) for every step inside it -- option (a) from
the brief's three choices, chosen and justified over (b)/(c): the raw
Fang jump, EXACTLY propagated (`check_r6o_arithmetic.py`'s matrix-
exponential model), already gives the right answer (6.017 V/V) -- the
bug was this file's DISCRETE (trapezoidal) propagation of an
already-correct jump, not the jump target itself, so a reduced-Jacobian
quasi-static elimination would have discarded a correct formula to fix a
propagation defect. This is the same "order 1 right after a breakpoint"
idiom `dctran.c`'s own XSPICE path already uses for the real trajectory,
applied here to this file's separate sensitivity/phasor recursion (never
touching `ckt->CKTorder` itself, so the real trajectory is unaffected).

**Result**: `pac_r6o_ctrl2out.cir` at 1kHz: **6.2044 V/V, -1.888 deg**
(13.1885/-1.888 unchanged through the first three passes). Matches a
from-scratch, two-saltation-event-per-period exact model (this pass's
extension to `check_r6o_arithmetic.py` -- see below) to **0.001dB/
0.001deg**, and the second pass's three independent references (6.00
analytic, 5.93 FD-of-orbits, 6.39 `pssac_ref.py`) within their own
mutual spread. **Step-size independence** (reltol 1e-3/1e-4/1e-5/1e-6,
the exact table the third pass showed drifting 13.19->13.38->18.02
V/V): now 6.2044/6.2045/6.2045/6.2044 V/V, spread **0.0002dB** -- the
validation task the brief's step (1) asks for, actually met.
**Red-proof**: `PSS_SALTATION=0` still gives exactly 0.

**Reference model extended, and why**: the fourth pass's single-event
Python model matched well at 1kHz (~3%) but the gap did not shrink with
frequency the way discretisation error would -- investigated and
explained, not left as a mystery: `Vramp`'s PULSE source has a SECOND
real crossing near period end (the fast fall-back edge), confirmed
present in the actual C code via `PSS_SALT_TRACE=1` (`h_dot=
-3.33333e+07` at `t=1.98465e-6s`, matching `-vramp/(0.03*tclk)` exactly)
but not modelled by the fourth pass's single-event reference. Adding a
second event to the Python model (same combined-KCL derivation, mirror
sign) closed the gap from ~3% to 0.001dB at 1kHz and to <=0.15dB from
100Hz-25.5kHz -- the C code was right; the earlier REFERENCE was
incomplete, not a second code bug.

**Known residual, honestly not root-caused this pass**: code-vs-model
gap grows above 25.5kHz -- 0.55dB at 50kHz, 2.0dB/2.8deg at 100kHz
(fs/5, fs=500kHz), OUT of the 0.5dB/3deg target at those two points
only. Ruled out with evidence: a PSS/PAC grid-resolution artifact
(`points` 200->2000, `PAC_RECORD_POINTS`=20000: <0.001dB change); the
pre-existing stale-Jacobian accuracy floor this file's own PAC section
above already documents for other fixtures (`PAC_EXACT_JACOBIAN=1`:
<0.001dB change). Leading unconfirmed suspect: a sampled-data/aliasing
effect as `w*T` approaches O(1) (1.257 rad at 100kHz), not chased
further this pass. The 0.5dB/3deg target IS met from 100Hz to 25.5kHz --
over a decade below fs/2, covering this converter's actual control
bandwidth (LC corner ~5kHz).

**Regression**: unchanged from the third/fourth passes' baseline, plus
three new automated checks the brief asked for (`check_pac_r6o_
saltation_gain/_stepsize/_redproof()`, `tests/test_pac.py`).
`test_pac.py` 21/21+3info -> **24/24+3info** (also reconfirmed clean
under a new `aaba50ba3f49-salt5-asan` build, `ASAN_OPTIONS=
detect_leaks=0` for the XSPICE `cmpp` code-model preprocessor's own
known unrelated build-time leak); `test_pss.py` 40/44+4info;
`test_checkpoint.py` 7/7; `test_sens.py` 27/27 gating; `test_reset.py`
OVERALL PASS (4/4); `test_compopt.py` 20/20; `test_stb_native.py`
28/29; `test_stb.py` 60/60; `test_stb_native2.py` 54/54;
`test_stb_device.py` 16/16; `test_stb_followups.py` 29/29;
`test_stb_cli.py` 42/42 -- every number matches the pre-this-pass
baseline exactly.

**Closed-loop R6s re-tried, confirmed still failing, confirmed
orthogonal**: `pss_r6s_shoot.cir` on this pass's build still fails
identically (`Timestep too small ... trouble with swmod-instance shs`,
error 106). This is expected, not a new gap: plain `.pss method=shoot`
convergence uses its own default finite-difference monodromy, a path
this pass's changes (both gated behind `PSSsensActive`/`PSSpacActive`,
analytic-Jacobian-only) never touch.

**NOT reached this pass**: CSW; the closed-loop R6/R6s case (confirmed
orthogonal to this pass, above, not further investigated); the XSPICE
`d_pwm` segfault (lowest priority in the brief, explicitly "only if time
remains" -- it did not); the full 0.5dB/3deg-vs-FRA table across a
closed-loop orbit (still blocked on the convergence gap); the >=50kHz
residual above (found, bounded, three hypotheses ruled out, not
root-caused).

### SIXTH pass (card #494 sixth pass, plan T6.1): >=50kHz residual ROOT-CAUSED -- it was the Python model, not the code; validated-band guard shipped; closed-loop diagnosed further, not fixed

HEAD is now `8cd652490` (one new commit on `aug-salt`,
`patches-wip/salt/0904-...patch`, series file updated; build
`aaba50ba3f49-salt6`). Full derivation and every number in the table
below: this pass's own working notes and `patches-wip/salt/README.md`'s
mirrored section; summary here.

**Task 1 (HIGH-FREQUENCY TRUTH): resolved.** Ran `tools/pac_perturbed.py`
(a third, independent method -- nonlinear perturb-and-DFT, shares no code
with `pacshoot.c` at all) on `pac_r6o_ctrl2out.cir` at 25/50/100/150/
200kHz plus seven more points chosen to approach fs/2=250kHz closely
(218750, 235294, 243590, 246914, 248744, 249000-snapped-to-250000Hz), and
cross-checked every one against a matching native `.pac solver=td` run at
the EXACT SAME achieved frequency (not the nearest round number). Result:

| f (Hz) | f/(fs/2) | native .pac (dB/deg) | pac_perturbed (dB/deg) | delta |
|---|---|---|---|---|
| 25000 | 0.100 | -11.850/-168.786 | -11.847/-168.798 | -0.003dB/+0.012deg |
| 50000 | 0.200 | -23.859/-161.450 | -23.859/-161.472 | +0.004dB/+0.022deg |
| 100000 | 0.400 | -34.926/-147.307 | -34.926/-147.349 | +0.001dB/+0.042deg |
| 150000 | 0.600 | -40.662/-136.330 | -40.662/-136.392 | -0.002dB/+0.062deg |
| 200000 | 0.800 | -44.312/-128.238 | -44.312/-128.395 | -0.000dB/+0.157deg |
| 218750 | 0.875 | -45.377/-125.789 | -45.382/-125.790 | +0.006dB/+0.001deg |
| 235294 | 0.941 | -46.219/-123.844 | -46.212/-123.866 | -0.007dB/+0.022deg |
| 243590 | 0.974 | -46.613/-122.936 | -46.615/-122.966 | +0.003dB/+0.029deg |
| 246914 | 0.988 | -46.766/-122.585 | -46.763/-122.598 | -0.003dB/+0.013deg |
| 248744 | 0.995 | -46.848/-122.394 | -46.560/-122.311 | -0.289dB/-0.083deg |
| 250000 | 1.000 | -47.358/-122.265* | -43.682/-80.801 | -3.223dB/-41.464deg |

(*250000Hz native magnitude/phase computed from the raw `v(vout)` complex
pair the same way as every other row, shown for completeness -- this row
is the failure case, not a passing reference point.)

**Verdict: the C CODE was right all along; the two-event Python toy model
(`check_r6o_arithmetic.py`, used only as a hand-derivation sanity check,
never part of the shipped engine) was the one that drifted at high
frequency.** Native `.pac` agrees with the independent `pac_perturbed.py`
reference to <=0.03dB/0.16deg from 25kHz all the way to ~247kHz (98.8% of
fs/2) -- FAR better than the 0.5dB/3deg target, and far better than the
"0.55dB@50kHz, 2.0dB@100kHz vs the Python model" residual the third
through fifth passes chased without success. Only in the last ~1-2% below
fs/2 does real, structural disagreement appear (0.289dB at 99.5%,
3.2dB/41deg at fs/2 exactly) -- and BOTH methods break down there for the
same reason: a genuine sideband-0/sideband-(-1) aliasing degeneracy for a
REAL system probed at exactly half the switching frequency (at w=ws/2,
the k=0 and k=-1 sidebands coincide in frequency for a real-valued
response -- neither this engine nor `pac_perturbed.py` computes sideband
coupling, see this file's own "Known limits" section, which already
documented this in general for smooth-switch decks; this pass shows it
applies to the saltation-corrected case too, and pins down almost exactly
WHERE it starts to matter). This is a single, measure-zero singular point
in principle, not a growing discretisation error -- empirically confirmed
by how flat the delta stays right up to 98.8% of fs/2 and then how
sharply it turns over in the last percent or two.

**Why the Python model drifted (found, not fixed).**
`check_r6o_arithmetic.py`'s `gain_at_freq_2event()` computes the
saltation-corrected PAC gain by propagating a complex "envelope" state
through each smooth interval using the REAL matrix exponential
(`expm(A, dt)`, no `-jw` term) and deferring the `e^{-jwT}` periodicity
phase factor to a single multiply at the very end of each period -- a
valid algebraic trick for the STATE recursion itself (`jwI` commutes with
every real `A`, verified by hand, term-by-term, this pass) but NOT valid
for `avg_over()`, the helper that cycle-averages the envelope to extract
the reported sideband-0 gain: integrating the UNPHASED trajectory and
correcting only the final sum is an approximation whose own error grows
with `w*T` (the original single-event function's OWN code comment already
said this integration is "valid...since w*T<<1 at 1kHz" -- true at 1kHz
(`w*T`=0.013rad) but not at 100kHz (`w*T`=1.257rad), and the fix was never
applied to `avg_over()` when the second event was added in the fifth
pass). An attempted fix this pass (weighting each sub-step's contribution
by `exp(-jw*(t0+s))`, `t0` the interval's own absolute start time) was
implemented and tested but made the LOW-frequency agreement (which was
already excellent, 0.001dB/0.001deg at 1kHz) measurably WORSE rather than
better, meaning the correction's sign/frame convention was not correctly
derived this pass either -- reported honestly rather than shipped broken.
`check_r6o_arithmetic.py` is UNCHANGED this pass (no functional edit); the
practical fix applied instead was to stop using it as ground truth above
25kHz and use `pac_perturbed.py`'s own measured values there (see
`tests/test_pac.py`'s updated `R6O_CTRL2OUT_MODEL`, no longer excludes
50/100/150/200kHz -- all four now pass against the corrected reference,
worst case 0.157deg).

**Task 2 (VALIDATED-BAND GUARD): shipped.** `patches-wip/salt/
0904-...patch` adds `PSSsaltEventsSeen()` (pssshoot.c/.h, a counter of
real saltation events found since the last `PSSsaltReset()`, 0 for every
existing no-`sw`-event fixture -- exact no-op there) and a new stderr
WARNING in `pacshoot.c`'s per-frequency boundary-solve loop
(`pac_engine_run()`, shared by both `.pac` and `.pstb`): whenever a
saltation-affected deck's requested frequency reaches 90% of the
switching Nyquist frequency (`f >= 0.9*fs/2`), one `PAC-TD: WARNING --
saltation-corrected result at f=... is inside the UNVALIDATED margin...`
line is printed, naming the exact ratio and pointing at this section.
0.90 was chosen deliberately conservative relative to the ~0.995 actually
measured clean above (only two points near the edge were checked, and
only on one fixture -- R6o) -- never silent, and gated strictly to
decks with real saltation events (confirmed: no warning on
`pac_r3_ctrl2out.cir`, a smooth-switch deck, even swept right up to ITS
OWN fs/2). New automated check: `tests/test_pac.py`'s
`check_pac_r6o_nyquist_guard()` (three sub-checks: fires at fs/2 on R6o,
silent at 40%-of-fs/2 on R6o, silent on R3 near ITS fs/2) -- passes.

**Task 3 (closed-loop R6s convergence): diagnosed further, still not
fixed.** Reproduced the pre-existing failure exactly (`pss_r6s_shoot.cir`,
`ngspice error 106`, "trouble with swmod-instance shs", timestep
collapsing to ~1e-18s at t=1.98492e-6s -- suspiciously close to R6o's own
independently-measured fall-edge saltation event time, t_e2=1.9846e-6s).
New this pass: (a) a PLAIN `.tran uic` of the identical circuit (no PSS
machinery at all) runs the full 50us requested window CLEANLY through and
past that exact instant with no timestep trouble -- this rules out a
general ideal-switch/feedback-loop stiffness problem in the transient
engine itself; the failure is specific to how `.pss method=shoot`'s own
warm-up ("stage 3/probe", "stage 4/final") period-integration calls this
circuit, not a property of the circuit's transient response in isolation.
(b) An experiment shifting `Vramp`'s phase by half a period (`TD=0.5*tclk`
added to its `PULSE(...)`, tested in a scratch copy, NOT committed --
still does not converge, see below) got PAST the warm-up-stage failure
that has blocked every prior pass (stage 3/4 now report a real Newton
residual instead of failing to even complete their own probe transient),
but surfaced a SECOND, independent failure of the exact same signature
(`ngspice error 106`, same switch) inside the DEFAULT finite-difference
monodromy's own per-state `+h` perturbation column (state `C:c2f`, one of
the compensator's own capacitors), now at t~2e-11s -- i.e. essentially
the very first step of that perturbed period run. This is consistent with
the brief's own leading hypothesis ("FD monodromy on the forced grid with
events that move between perturbed runs") but pins it down further: the
`+h` perturbation on a compensator state evidently moves the very FIRST
switch decision of that column's period run into the same numerically
pathological regime the unperturbed warm-up hits near the ramp's fast
reset edge. **Not fixed**: a real fix needs either (i) precise event-time
location in the DEFAULT (non-saltation) FD-monodromy path too (the same
"no `CKTsetBreak`-based crossing time" gap this card's OWN saltation
machinery was built to address, but for the *unperturbed* trajectory's own
step control, not the sensitivity/phasor recursion), or (ii) a smaller/
adaptive FD step `h` for columns whose perturbation is found to flip a
switch decision, neither attempted (both are changes to the
DEFAULT/shared FD-monodromy code path every existing fixture depends on --
too invasive to risk in the time remaining this pass without a much
larger regression budget). The phase-shifted fixture variant is NOT
committed (it does not converge either, so it is not a net improvement,
only a diagnostic data point); `tests/fixtures/pss_r6s_shoot.cir` is
UNCHANGED this pass.

**Task 4 (XSPICE `d_pwm` segfault): not attempted**, per the brief's own
explicit priority ("only if time remains" -- it did not, given tasks 1-3
above).

**Regression on this pass's build (`aaba50ba3f49-salt6`)**: `test_pac.py`
24/24+3info -> **25/25+3info** (the new `check_pac_r6o_nyquist_guard()`
check, plus `R6O_CTRL2OUT_MODEL`'s 50/100/150/200kHz entries switched from
the flawed Python model to `pac_perturbed.py`'s measured values and
un-excluded -- all pass, worst case 0.157deg); `test_pss.py` 40/44+4info;
`test_checkpoint.py` 7/7; `test_sens.py` 27/27 gating; `test_reset.py`
OVERALL PASS (4/4); `test_compopt.py` 20/20; `test_stb_native.py` 28/29;
`test_stb.py` 60/60; `test_stb_native2.py` 54/54; `test_stb_device.py`
16/16; `test_stb_followups.py` 29/29; `test_stb_cli.py` 42/42 -- every
number matches the fifth pass's baseline exactly except the intentional
`test_pac.py` improvements above (patch 0904 is warning-only, zero
behavioural change to any computed result, confirmed by every other
number being byte-identical to the fifth pass).

### card #510 (high priority, silently-wrong result): the input-forcing/inductor-routed case -- a second, distinct saltation gap found and fixed for the capacitor-free switching node

Card #510: `tests/fixtures/pac_sw_ctrl_direct.cir` (open-loop ideal-switch
buck, single `sw` S1 straight into an inductor L1 with a freewheeling
diode D1, no snubber capacitor anywhere on the switch's own power node)
gave `.pac`/`.pxf` control-to-output **exactly 0+0j at every frequency**,
with no warning, while a `.tran` shows a real ~13 V/V DC gain and the
saltation machinery's own event-detection (fs/2-band warning, `PSS_SALT_
TRACE`) confirms it DOES see the switching event. Two fixture bugs fixed
first, per the card's own instruction ("fix the fixture... and check
both"): the model line used `vt1`/`vt2` (not `SW`'s real parameter names,
`vt`/`vh` -- ngspice printed `unrecognized parameter, ignored` and fell
back to the SAME defaults, `VT=0 VH=0`, confirmed bit-identical either
way, so this was never the cause); and `print out` (no such vector) ->
`print vout`.

**Cause, confirmed with an experiment.** The card's own hypothesis --
"the AC source enters the switch's event function DIRECTLY (control
voltage = input u, no state in between)... the event-time shift dtau...
is missing" -- was tested directly and REFUTED: inserting a tiny RC
between `Vctrl` and the switch's control node, at RC time constants from
1ns down to 1ps, did NOT change the zero result AT ALL, at any RC value
(a state-driven control path should have made the gain approach the
right answer as RC->0; it stayed exactly 0+0j throughout). `PSS_SALT_
TRACE=1` showed why: the saltation machinery DOES locate the event
(`n_events=1`) and DOES compute a nonzero `coeff` (the event-time shift,
`dtau`, from `y_prev[posCntrl]-y_prev[negCntrl]` -- this already reaches
a directly-forced control node correctly, via the SAME full-circuit
per-step `SMPsolve()` that already handles it for the validated
R6o fixture, so no change was needed there), but `dfdt_pos=0 dfdt_neg=0`
-- the (f+-f-) JUMP TARGET itself is zero. Root cause: `pss_salt_local_
cap()` (THIRD pass) returns 0 whenever a switch's power node has no
DIRECTLY-ATTACHED capacitor -- correct for an exogenous node (an
independent source's own terminal genuinely cannot jump, its comment
already says so) but WRONG for a node like this fixture's `sw` (S1's
negNode), which has no capacitor but IS the near terminal of inductor
L1 -- the one genuine, persistent STATE adjacent to the event. Confirmed
by the complementary experiment: adding a capacitor directly at node
`sw` (1nF down to 100fF) DOES make the result nonzero (though not a
clean answer on its own -- see below).

A second, deeper reason this matters, found while deriving the fix:
`PSSsensRecordStep()`/`PSSpacRecordStep()`'s own smooth per-step
recursion only ever READS a state's history (`g_sens_Sprev`/`g_pac_
Cprev`) at STATE rows -- capacitor node pairs and inductor branch-
current rows, the `for (d=0; d<n; d++)` loop in both files. A jump
written onto a NON-state node's own row (as the code did before this
fix, always targeting `posNode`/`negNode` directly) is therefore INERT
for propagation whenever that node is not itself a state: it is
overwritten from scratch, correctly, by the very next accepted step's
own `SMPsolve()`. This is exactly why R6o's own "0 at node `vin`" case
was harmless (nothing forward-propagating was ever lost there -- `vin`
is an exogenous source node) but is a genuine bug at an inductor-
adjacent node with no cap: the ONE state that DOES persist there (the
inductor's own current) was never given a jump at all, at any
capacitance value -- which is also WHY the "add a capacitor" experiment
above does not, by itself, cleanly fix anything: it changes the state
topology (adds a new fast state, interacting with `.pss`'s own fast-
state elimination and the L-stable settle window) rather than modelling
the real circuit.

**Fix, `pss_salt_ind_jump()` (`pssshoot.c`, shared by `PSSsensRecordStep()`
and `PSSpacRecordStep()`, hence `.pss` monodromy, `.pac`, and `.pxf`'s v1
fallback -- card #504's true-adjoint engine already REFUSES any
saltation-affected deck and falls back to v1, so `.pxf` on this fixture
always takes the fixed path; confirmed end to end, not assumed, by a new
`.pxf` check, see Validation below).** When a switch's power node has no
local capacitor, look for a UNIQUE adjacent inductor and route the
`(f+-f-)*dt_e` jump onto ITS OWN branch-current state instead, in
inductor-current-derivative units (`f_L=(V(pos)-V(neg))/L`).
`f_L+-f_L- = (V(node)+ - V(node)-)/L`; `V(node)+` is simply the ACTUAL,
already-solved accepted-point value (`ckt->CKTrhsOld[node]`) -- no
approximation needed there. Only `V(node)-` (the COUNTERFACTUAL "what
would this node's voltage be at this same instant, holding the adjacent
inductor's current `I_L` at its continuous, unable-to-jump actual value,
if the triggering switch had NOT flipped") needs deriving, from a LOCAL
one-node KCL solve at that node: the triggering switch (conductance
`g_prev`, to a verified-stable `otherNode`) plus every OTHER device
actually touching the node (resistor / other switch / diode) plus the
fixed `I_L`, using PRE-EVENT (`CKTstate1`, NOT `CKTstate0`) small-signal
Norton equivalents for every nonlinear one of those:

```
V(node)- = [g_prev*V(otherNode) + K_other_prev - s*I_L] / (g_prev + G_other_prev)
```

(`s = +1` if `node` is the found inductor's own `posNode`, else `-1`;
`K_other_prev`/`G_other_prev` = the sum of every other device's own
PRE-event Norton current/conductance at `node`, from its `CKTstate1`
history -- for a diode, `K = I_D_prev + G_D_prev*(V_far - V_D_prev)` on
the cathode side, `K = -I_D_prev + G_D_prev*(V_far + V_D_prev)` on an
`RS=0` anode; full sign derivation for both cases in `pssshoot.c`'s own
comment above `pss_salt_ind_jump()`). No `ag0`, no step size, anywhere --
the same step-size-independence property the capacitor case has.

**A second bug found while validating, not guessed at and shipped
blind.** A SIMPLER first version of this fix used only the POST-event
`G_other` (`dV=dI/(g_prev+G_other_now)`, the direct analogue of the
capacitor formula with `C` replaced by `G`) -- this gave a WILDLY wrong
answer, many orders of magnitude off, specifically on this fixture's
turn-ON edge (`g_prev` there is `SWoffConduct`~1e-6, so the result is
almost entirely determined by `G_other`), traced to the freewheeling
diode at the SAME node effectively co-switching in sympathy with S1
(forward-biased, carrying the full inductor current, right before the
edge; reverse-biased, ~0 conductance, right after) -- the "conductance
stays fixed across the event" assumption that formula makes for every
OTHER device at the node is badly violated exactly when that other
device is itself a synchronous/freewheeling element, an entirely
ordinary situation in a switching converter. Caught by cross-validating
against `tools/pac_perturbed.py` (see Validation below) on a purpose-
built diagnostic variant (a symmetric triangular-ramp version of this
fixture, both switching edges placed safely mid-period, away from any
period-restart artifact -- not committed, a scratch validation aid) --
NOT by trusting the first formula that merely stopped returning exactly
0. The PRE-event (`CKTstate1`) Norton formula above does not have this
failure mode: it solves the counterfactual PRE-event topology's own node
equation using THAT topology's own (correct) diode operating point, not
one borrowed from the post-event topology.

**A third, separate gap found and fixed the same pass: the FIFTH pass's
L-stable settle window never activated for this case.** `PSSsaltFindEvents()`'s
`tau_max` (feeding `PSSsaltForceOrder1()`'s post-event backward-Euler
window, the FIFTH pass's fix for TRAPEZOIDAL ringing right after a jump
at a stiff node) was computed ONLY from `pss_salt_local_cap()`/
`pss_salt_local_cond()` -- both 0 whenever there is no local capacitor,
so this fixture's inductor-routed jump got NO settle window at all,
letting the surrounding order-2 recursion ring into the slow states the
same way the FOURTH/FIFTH pass found for R6o's `Csw`. Fixed by having
`pss_salt_ind_jump()` also return an L/R-based tau estimate
(`L/(g_prev+G_other_prev)`), wired into the SAME `tau_max`/`g_salt_
settle_until` logic the capacitor case already uses. Measured effect:
cut the worst-case error against `tools/pac_perturbed.py` from ~6dB to
~2-3.6dB (see the validation table below) -- a real, large improvement,
not fully closing the gap.

**Validation (`tools/pac_perturbed.py`, independent, nonlinear,
shares no code with `pacshoot.c`, s55 baseline binary), 1k/10k/50k/100k/
200kHz, all below 0.9*fs/2 (fs=500kHz):**

| f (Hz) | code (dB/deg) | `pac_perturbed.py` (dB/deg) | delta (dB/deg) |
|---|---|---|---|
| 1000   | 25.197 / -3.589   | 22.405 / -2.604   | +2.79 / -0.99 |
| 10000  | 9.048 / -169.595  | 12.613 / -166.481 | -3.57 / -3.11 |
| 50000  | -19.624 / -158.938 | -17.446 / -161.514 | -2.18 / +2.58 |
| 100000 | -30.375 / -143.695 | -28.276 / -145.079 | -2.10 / +1.38 |
| 200000 | -39.398 / -125.594 | -38.156 / -126.844 | -1.24 / +1.25 |

Worst case 3.57dB / 3.11deg at 10kHz -- **does NOT meet the card's own
0.5dB/3deg target** (reported honestly, as this file's own established
convention for an unmet target requires -- see the module docstring of
`tests/test_pac.py`). What IS achieved and validated: the silent
exact-zero (the card's own stated priority, "a silently wrong result")
is eliminated; the result is correctly SIGNED, correctly PHASED (within
a few degrees, closer than the dB error), of the correct ORDER OF
MAGNITUDE, and correctly DECAYS with frequency, matching the physical
`.tran` DC-gain estimate's order of magnitude. Line-to-output (`Vin`,
via `.pxf`'s `pxf_vin`) and the R6o capacitor-case fixture (`pac_
saltation(a-d)` in `tests/test_pac.py`) are confirmed byte-identical to
the pre-this-pass baseline (this fix's inductor-routed path never
triggers for either: `vin` has no local cap AND no adjacent inductor,
correctly a no-op per the "purely algebraic node" case; R6o's switch
nodes both have local capacitors, so the THIRD pass's original formula
is used, unchanged).

**Not root-caused further this pass (time budget).** Leading, NOT
confirmed, suspect for the remaining ~2-3.6dB: the settle window's own
size. `PSS_SALT_SETTLE_TAUS=30` was tuned (FIFTH pass) against R6o's
`Csw`-based tau (nanoseconds, far shorter than one switching period);
this fixture's L/R-based tau (`L1=10uH` against the switch/diode's own
conductances) can be hundreds of nanoseconds to microseconds -- 30x that
can extend the forced-backward-Euler window across MULTIPLE subsequent
switching periods, an untested regime the FIFTH pass's own tuning never
exercised. A quick, uncommitted diagnostic (`PSS_SALT_NOFILTER=1`, a
one-line opt-in bypass of the THIRD pass's restart-event filter, added
this pass and left in as a diagnostic aid, default-off, zero behaviour
change otherwise) also surfaces a related, SEPARATE, unresolved
question: this fixture's OWN fast-fall ramp edge (the "turn back ON"
transition) lands essentially exactly at the period boundary (`t=T`,
by construction -- `Vramp`'s `PULSE(0 1 0 2u 1n 1n 2u)` reaches its peak
at exactly `t=2us=T` and resets within ~2ns), and the THIRD pass's
restart-event filter (which discards the first `PSSsaltFindEvents()`
call after every period reset, to suppress a known stale-`CKTstate1`
artifact) currently discards THAT edge along with the artifact,
confirmed via the diagnostic (with the filter bypassed, an event IS
found there, but its own `h_dot` is NOT reproducible between calls --
4.9995e9 vs 9.9995e9 V/s for what should be the same instant -- itself
consistent with the FIFTH pass's own root-caused "`CKTdelta` freshly
after a reset has no valid true-elapsed-time denominator yet" gap, never
fully closed for the FIRST call after a restart specifically). Whether
fixing THAT gap (not attempted this pass -- touches the shared restart-
filter/event-location machinery every existing saltation fixture
depends on, judged too invasive to risk without a much larger regression
budget than this pass had) would also close the remaining dB gap is
unknown -- a genuinely open question for a future pass, not conflated
with the settle-window-size suspicion above.

**Regression.** `test_pac.py`: R6o's four `pac_saltation(a-d)` checks
byte-identical to the pre-this-pass baseline (0.145dB/0.157deg worst
case, unchanged); the two new checks this pass adds
(`check_pac_sw_ctrl_direct_saltation()`, `check_pac_sw_ctrl_direct_
redproof()`) report accurately (an honest FAIL on the 0.5dB/3deg target,
a passing red-proof both via `PSS_SALTATION=0` and via literally running
the named pre-fix binary, `aaba50ba3f49-s55`, confirming it gives
exactly 0). `test_pxf.py`: two new checks
(`check_pxf_sw_ctrl_direct_nonzero()`, `check_pxf_sw_ctrl_direct_
redproof()`) pass -- `pxf_vctrl`/`pxf_vramp` nonzero and exactly mirror-
symmetric (`|pxf_vctrl+pxf_vramp|/|pxf_vctrl|=0`, `h=v_ctrl-v_ramp`'s own
differential-control symmetry), matching `.pac`'s own vout gain to
machine precision (same engine). `test_pss.py` 40/45+5info,
`test_pac_perturbed.py` (unrelated pre-existing `campaign.py` manifest-
directory bug on `pac_r6o_ctrl2out.cir`'s own tran-fallback path, `(d)`,
not touched by or related to this card -- COULD-NOT-RUN, not a
regression), `test_stb_native.py` 28/29, `test_checkpoint.py` 7/7 --
every number matches the pre-this-pass baseline exactly except the
intentional `test_pac.py`/`test_pxf.py` additions above. The two
existing `test_pac.py` R5 free-running FAILs and the two existing
`test_pxf.py` card #502 harmonics=1 ASan FAILs are unchanged,
pre-existing, and unrelated to this card (confirmed identical on the
s55 baseline).

**What was skipped.** CSW (current-controlled switch), XSPICE: not
reached (this pass's scope was the `sw` case the card's own fixture
needs). A resistor/other-switch's own far terminal, or a diode's own
far terminal, that is not itself ground/a source/capacitor-anchored is
REFUSED, not modelled (a genuinely coupled multi-node event this local
solve does not attempt); a `RS!=0` diode anode is REFUSED for the same
reason (only the cathode side, and an `RS=0` anode, are exact -- see
`pss_salt_other_cond()`... i.e. the inline diode loop in `pss_salt_ind_
jump()`'s own comment). A non-triggering switch's own PRE-event
conductance is approximated by its CURRENT `SWcond` (no separate
`CKTstate` history kept for that field, unlike the diode's) -- fine
unless that other switch is ALSO flipping at the exact same instant, not
exercised by any fixture checked. The remaining ~2-3.6dB gap above is
not closed. The restart-filter/period-boundary interaction is diagnosed
(`PSS_SALT_NOFILTER=1`) but not fixed.

**Step 2 result, 2026-09-29 (card #510, aug-salt2 commit 251926ea4, to be
patch 0062): accuracy target MET.** The remaining ~2-3.6dB gap above was
the settle window's sizing, as suspected. Step 1 sized the L-stable
(backward-Euler) post-event window like the capacitor case, 30 time
constants with tau=L/(g_prev+Gsum), about 15us here = 7.5 switching
periods of forced order 1. That is wrong for an inductor: its L/R tau is
large by design, so |lambda*h|<<1 and there is no stiff ring to damp; what
the window protects against is reusing the pre-event order-2 history term
for the first step or two after the jump. The window is now sized in
steps: `PSS_SALT_IND_STEPS` (diagnostic env, default 2) x this event's
accepted step `dt_actual` (`pss_salt_ind_jump()` gained a `dt_actual`
argument). Per the commit's sweep: N=1 is a cliff (~6.5dB), N>=1.5 plateaus
at ~0.45-0.48dB; N=2 shipped. The restart-filter/period-boundary
suspect (`PSS_SALT_NOFILTER`) was not needed and was not re-examined.

Final accuracy, `pac_sw_ctrl_direct.cir` vs `tools/pac_perturbed.py`
(build aaba50ba3f49-salt2-1; code minus reference):

| f (Hz) | code dB / deg | ref dB / deg | error dB / deg |
|---|---|---|---|
| 1000 | 22.475 / -2.540 | 22.405 / -2.604 | +0.070 / +0.064 |
| 10000 | 12.618 / -167.088 | 12.613 / -166.481 | +0.005 / -0.607 |
| 50000 | -17.269 / -161.246 | -17.446 / -161.514 | +0.177 / +0.268 |
| 100000 | -28.319 / -147.452 | -28.276 / -145.079 | -0.043 / -2.373 |
| 200000 | -37.672 / -129.041 | -38.156 / -126.844 | +0.484 / -2.197 |

Worst 0.484dB / 2.373deg (target 0.5 / 3); the 200kHz point has only 0.016dB
of margin. The superseded ~3.6dB/3.1deg figures and the "remaining gap not
closed" statements above describe step 1. `tests/test_pac.py`'s accuracy
check for this fixture is now gating (FAIL above 0.5dB/3deg).

The inductor-routed WARNING (once per analysis) is deliberately still
printed: only this one fixture validates the inductor-routed path (one
inductor, one freewheeling diode, one switch), so the path stays flagged as
less validated than the capacitor case; its test still gates. Regression:
R6o results byte-identical between aaba50ba3f49-s61 and -salt2-1 (all eight
`pac_r6o_ctrl2out.cir` frequencies 100/1k/10k/25.5k/50k/100k/150k/200kHz
`pac_out.txt`, and `pss_r6o_shoot.cir`'s `.txt` and `.raw`, the latter
identical except the Date/Command header lines); test_pac 31/31,
test_pss 40/45, test_pxf 12/14 (2 ASan card #502 FAILs, same on s61),
test_pac_perturbed 8/8, test_stb_native 28/29, test_checkpoint 7/7. Not done:
other topologies (several inductors, an RS!=0 diode anode) are still refused
or unvalidated; N was not re-swept in this pass.

## card #498: closed-loop ideal-switch PSS (exact `sw` event location)

**Result (2026-09-29, second pass).** `tests/fixtures/pss_r6s_shoot.cir` (the
closed-loop voltage-mode buck built from two native `sw` devices whose
control is `opout` vs the ramp) now converges with `.pss method=shoot` at the
deck's own tolerance: iter 3, residual 8.3e-9, spectral radius 0.980
(stable), 0.4 s (Points 200; Points 100 and 400 also converge). Before this
card it plateaued at 1.8e-3 (Points 200) / 9.3e-4 (Points 2000) or aborted.
The orbit agrees with the last period of a 1.6 ms settled `.tran`: vout
1.7e-4 V, i(L1) 1.7e-3 A (on ~0.48 A ripple), v(opout) 1.0e-3 V (max
difference). `.pstb` on the same orbit agrees with `tools/pac_perturbed.py
--loop-gain` to 0.13 dB / 0.26 deg over 10k-200k (`test_pac.py`,
`check_pstb_r6s_closed_loop`).

**Root cause.** A native `sw` flips at the first accepted time point that sees
its control past the threshold, so the switching instant is quantised to the
adaptive grid. When the control is itself a state (a comparator in a closed
loop), the period map x0 -> x(T) is piecewise constant in the event time and
every finite-difference monodromy column sees either zero or a whole-step jump
of d(t_event)/d(x0): the Newton residual plateaus at a Points-dependent floor.
Two aggravating effects were found on the way: pass 1's looser-reltol retry
produced Jacobian columns and line-search trials computed at a different
tolerance from the baseline (garbage FD columns), and after a flip the
trapezoidal LTE history straddles a ~1 kA current jump in the snubber/`Csw`
caps, so the step collapsed to `CKTdelmin` ("Timestep too small ... swmod").

**Fix** (`swload.c`, `dctran.c`, `pssshoot.c`, `pacshoot.c`; only for the
FD-Jacobian path and only when the deck has a native `sw`):

1. *Lagged switch decision* (`PSSswLagActive`): during a period run a switch
   changes state only if the control value at the previous **accepted** point
   was already beyond the threshold, so the state is constant over a step and
   the control voltage is a smooth function of time.
2. *Exact event location* (`PSSeventLocate`, called in `dctran.c` after a
   converged step): if a control crossed its flip threshold within the step,
   the step is rejected and re-run to the interpolated crossing time (Illinois
   regula falsi on the control value, bracket kept across trials, accepted when
   the overshoot is <= `PSS_EVENT_TTOL` seconds, default 1e-16). The next step
   starts exactly at the crossing and the flip happens there.
3. *Flat LTE history after a flip*: the first step after a real switch flip
   restarts with `CKTstate2/3 = CKTstate1` (as the first transient step does),
   otherwise the LTE sees the current jump as a huge second difference.
4. *Located crossings are not forced as breakpoints* in perturbed runs
   (each run locates its own crossing, so `d t_event / d x0` is kept).
5. *Honest verdict*: a converged verdict is refused if its baseline period run
   only completed at the looser retry reltol.
6. `.pac`/`.pstb` recording passes follow the same treatment
   (`PSSshootLagEnabled()`), so the small-signal orbit is the PSS orbit.
   `PSS_NO_EVENT_LOCATE=1` restores stock switching everywhere (used by the
   red-proof in `test_pss.py`). The analytic Jacobian
   (`PSS_SHOOT_JACOBIAN=analytic`) keeps stock switching.

**Convergence history (r6s, Points 200).** pass 1: plateau 1.8e-3 -> event
location alone (no lag/flat history): 2e-5 -> lag + flat history: converged
at iter 3 (8.3e-9). `PSS_NO_EVENT_LOCATE=1`: 6.8e-4 plateau, then a period
run aborts (red-proof in `test_pss.py`). `pss_r4_d55_unstable` still refuses.

**Tests.** `test_pss.py`: `check_r6s_closed_loop_converges` (converged,
orbit vs settled `.tran` bounds 1e-3 V / 6e-3 A / 4e-3 V) and
`check_r6s_event_location_redproof`. `test_pac.py`:
`check_pstb_r6s_closed_loop` (reference values stored in the test:
10k 15.86 dB/-127.48 deg, 20k 6.64/-132.03, 30k 1.50/-139.45, 100k
-16.63/-163.28, 200k -28.32/-171.42, from `pac_perturbed.py --backend tran`).
Regression versus s63: test_pss 40/45+, test_pac 31/31, test_pxf 14/14,
test_sens 27, test_checkpoint 7/7, test_stb_native 28/29 unchanged; ASan
clean.

**Not done / limits.** The analytic Jacobian still diverges on r6s, also with
exact event location forced on (experiment, reverted), so grid quantisation is
not its cause; suspects are the simultaneous complementary-switch pair and
`v_switch` sampled at the end of the flip step. Points 2000 does not converge.
The loose-reltol retry still fires for some warm-up/Jacobian-column/trial runs
(the verdict no longer depends on it). `csw`, XSPICE `d_pwm` and
negative-hysteresis `sw` models are not handled. The 5 kHz `.pstb` point was
not used (its `pac_perturbed` reference failed its own amplitude-halving check
by 3.5 dB); 1 kHz/3 kHz are not commensurate at qmax 100.

## card #505: realistic dead-time bucks (R3d/R4d) under shooting-Newton -- root causes and fixes

Symptoms (card #479 part 2 / #487): (a) the R4d `.pac`/`.pstb` fixtures stalled at residual 1.0 unless a
DC=0 series source was spliced at the comparator input; (b) R3d `.pss` plateaued at 3e-4 (and the warm-up
was non-monotonic in StabTime); (c) R4d vc orbit 0.21 % vs the 0.1 % target; (d) `tools/loss.py --pss`
failed on R3d/R4d.

**(a) is not a topology or state-selection effect.** State vectors, fast-state detection and the first
Newton step are identical with and without the splice to 4-6 digits (same 9 states, same reductions). The
first perturbed period run that differs is one FD monodromy column (`-h` on `cdly_ls`): the corrector
collapses at t = 39.25 ns ("Timestep too small ... node reset_sig", step 1e-19 s). Trace of that run
(temporary per-step dump, not committed): the node solution runs away exponentially (vc, vout ->
1e13 V) once ngspice's step has been cut below ~1 ps by corrector failures on the deck's very sharp
dead-time turn-on edge (`kswd=150`, `gon=70`). The runaway is a genuine property of the deck: **the op-amp
`Goa vc 0 vrefnode compin` (r3/r3d/r4/r4d/r6 all use it) has the wrong polarity -- it is a positive-feedback
amplifier through `C1f`/`C2f`, with a right-half-plane pole at ~1e12 s^-1.** Verified in stock
ngspice-45.2 with no PSS at all: `.tran 1e-14 4e-11 uic` on the deck gives v(vc) = -1e168 at 40 ps; with
the two Goa control nodes swapped it stays at 0.6 V. Also stock: the same deck at `reltol=1e-6 vntol=1e-9`
aborts with "Timestep too small" at 13 us, swapped it runs; R3d's `tran uic` aborts at every reltol <= 1e-4
(vout 1e274 at t = 8e-11). Coarse implicit steps hide the pole (the slow dynamics are those of an ideal
virtual-short op-amp either way -- identical vout/vc after 3 ms), so nothing looked wrong until an
integration step landed below ~1 ps. The splice merely changed the corrector's path (an extra branch
unknown) enough to miss that one collapse. The engine cannot remove a model instability; it can stop it
killing the Newton run.

**Engine fixes (src/spicelib/analysis/pssshoot.c, dctran.c, pssshoot.h):**
1. *One-sided FD rung.* If one side of a monodromy column fails on a switch-free deck, the surviving side
   is differenced against a **forced run at x0 itself** (same grid; differencing against the free baseline
   would divide the free-vs-forced mismatch by h) and the column is used; message
   "used a one-sided difference for this column". Only reached on failure. `PSS_NO_ONESIDED=1` disables.
2. *Noise-robust fast-state detection.* States whose magnitude is below their own swing (the R4 latch
   `cq`, the dead-time RC nodes) have period-map noise ~1e-5..3e-5 (their end value is exp(-K) in the
   comparator on-time) against an FD step of 1e-6, so their columns read 5-30 instead of ~0, they stayed
   Newton unknowns and their garbage diagonal produced Floquet multipliers of 2..14 "UNSTABLE" on an
   orbit whose multiplier is 0.994. Detection now re-measures every above-threshold column with a step of
   5 % of the state's ripple and keeps the smaller value. The Newton FD step for a state with
   |x| < ripple is raised to 1 % of the ripple. `PSS_SHOOT_NO_BIGH=1` disables both.
3. *Exact step/order replay* (`PSSreplayStep`). The classic "force every recorded time as a breakpoint"
   scheme restarts backward Euler at each forced point and re-derives the sub-steps by LTE, so a forced
   run at the SAME x0 is a different function from the free baseline (R3d: i(L1) differs by 1e-4 relative
   at alpha -> 0, trial residual 1.7e-3 vs baseline 3.2e-4 -- that mismatch IS the R3d plateau, not
   grid-quantised events). Replay integrates the baseline's recorded accepted times with the recorded
   integration order and no LTE rejection, so P(x0) reproduces the baseline. It is a *fallback rung*: only
   when the final classic Newton attempt fails, on switch-free decks (`sw` decks keep the #498 machinery):
   classic -> replay from the classic iterate -> replay from the warm-up endpoint. `PSS_REPLAY=1/0`
   forces it on/off. Not the default because on R4d it converges to the right residual but its FD
   Floquet is noisy (multiplier 11-14 spurious) and the emitted orbit was 1-2 % of ripple off.
4. `tools/loss.py --pss`: the word "error" in the engine's own recoverable diagnostics ("...(ngspice error
   106)") was treated as fatal; in `--pss` mode the verdict is now the `PSS-SHOOT: CONVERGED` line.
5. Fixtures: the card #487 `Vdummy` splice is removed from `rd_pac_r4d_line2out/zout`, `rd_pstb_r4d`.
   (`Vctrl` in `rd_pac_r4d_ctrl2out` is the real injection and stays.)

**Not done / not fixed:** the op-amp polarity in `refs/*.cir` (a project decision: the swap gives the
physically correct inverting amplifier, makes plain `.tran` robust at reltol 1e-5, but shifts R3d's
Floquet multiplier 0.973 -> 0.980 so the FRA/PAC references would need regenerating). A tightened `vntol`
for the Newton runs (`PSS_SHOOT_VNTOL`, off by default) lowered the R3d plateau (2e-5 -> 5e-6) but never
converged it and broke another deck at 1e-9. R3d's settled-`.tran` orbit reference cannot be generated from
the deck as written (see above), so the R3d orbit-vs-tran checks report COULD-NOT-RUN with that reason;
on the polarity-corrected twin the PSS orbit matches the settled reltol=1e-5 `.tran` to 0.021 % (vout),
0.019 % (i(L1)), 0.029 % (opout) of ripple.

**Fixed 2026-09-30 (card #519):** the Goa polarity is corrected in all five `refs/*.cir` and the 60 test
fixtures (`Goa 0 opout vrefnode compin`, `Goa 0 vc ...` in R4/R4d); `bench/decks` keep the old line so the
stored speed baselines stay comparable. Vout, ripple, fsw, PM and crossover move within the FRA noise; the
FRA/envelope references were regenerated (old copies kept as `*oldpol*` under ~/Work/data). The R3d
stabtime in the fixtures was raised (300 -> 500 us) to converge the shooter; that is a workaround, not an
engine fix. See `refs/README.md`.

## What was not done this pass (honest accounting)

- `pnoise`/`pxf`/`psp`, sidebands != 0, PXF/adjoint PAC, device-mode `.pstb`
  probing, GEAR integration: out of scope (see "Known limits" above).
- R4's line-to-output/output-impedance FRA references only exist at D~=0.30
  in this project's own data (`SUMMARY.md`); the D~=0.45 pac_sampling check
  above covers control-to-output only, per this card's own task scope.
- The corrected line-to-output FRA reference (VO-fixed) was generated for
  R3 and R4 D~=0.30 only, with the one-off script described above (not
  committed as a new `tools/` file this pass -- reusing
  `tools/fra_curves_fork.py`'s already-committed `instrument_r3()`/
  `instrument_r4()` plus a small local helper); a proper fix belongs in
  `tools/fra.py`'s own `_sin_spec()`/callers, which is out of scope for
  this project/card to edit (see the trap above).

## `.pnoise` on the shooting-PSS orbit: periodic (cyclostationary) noise (card #509)

```
.pnoise Fguess StabTime OscNode Points Harmonics SC_iter Steady_coeff \
        OutNode InSrc {dec|oct|lin} NumPts Fstart Fstop [cyclo] \
        [method="shoot" solver="td" [dev=1] [pacpoints=N]]
```

Spectre's `pnoise ... noisetype=timeaverage`: the time-averaged output
noise PSD at sideband 0 of one output node, folding every device's
cyclostationary noise (density evaluated at the operating point EACH
instant of the periodic orbit passes through, not just one DC bias) through
the LPTV system's own periodic Green's function. `.pnoise` **already
existed** in this fork's base (Enhancement-124/126/138/139/177/178,
`dcpss.c`'s `pnoise_sweep()`, the frequency-domain harmonic-conversion-
matrix engine, including its own already-sophisticated
Enhancement-178 "exact separable cyclostationary folding" for a colored
(non-white) source) -- **unchanged, and still the default** (`solver="fd"`
or `method`/`solver` omitted). `method="shoot" solver="td"` (both required,
same gating `.pac`/`.pxf`/`.pstb` already use) selects the NEW time-domain
shooting-based engine documented here instead, built on the `.pxf
adjoint=1` true-adjoint engine (card #504, `pacshoot.c`'s
`pac_engine_run_adjoint()`) -- additive, same pattern `.pxf`'s own
`solver=td` used on top of the pre-existing frequency-domain `.pxf`. Field
order/positions were kept from the pre-existing card (same reasoning as
`.pxf`'s own "Syntax decisions" item 1: extend additively with trailing
keyword tokens, don't invent a second, incompatible grammar for the same
dot-card name) -- the card's own suggested `.pnoise v(out[,ref]) dec N
fstart fstop` syntax was NOT adopted for this reason.

**`dev=1`**: additionally emits a per-device noise-generator contribution
table (`onoise_<inst><suffix>` vectors, e.g. `onoise_Rin_thermal`,
`onoise_Rin_1overf`, `onoise_Rin` (the per-instance total) -- the SAME
naming `NOISE_ADD_OUTVAR` already gives plain `.noise`'s own summary
report). Default (`dev=` omitted) emits `onoise_spectrum` only, cheaper (no
per-generator naming/summary pass). **`InSrc`** (input-referred spectrum)
is accepted syntactically (an old deck's field is not broken) but IGNORED
by `solver=td` -- no input-referred (`inoise`) output in this v1 engine;
give any valid source name, or reuse the deck's own probe/reference source.
**`pacpoints=N`**: same recording-pass grid override `.pac`/`.pstb`/`.pxf`
already have (`pacshoot.c`'s `pac_run_period()`), solver=td only.

### Design

`pac_engine_run_adjoint()`'s own backward sweep (card #504) computes, for
one fixed output row and frequency `w`, `Lambda(w) = Sum_k nu_k(w) *
e^{+j*w*t_k}`, where `nu_k(w)` is the per-step `SMPtSolve()` result at
backward step `k`, BEFORE that phase multiply. By that function's own
derivation, `Y0(w,b) = b^T Lambda(w)` for a source forced IDENTICALLY (the
same spatial pattern `b`, phase-modulated `e^{jwt_k}`) at every step.
Substituting `Lambda(w) = Sum_k nu_k*e^{jwt_k}` and regrouping:

```
Y0(w,b) = Sum_k (b*e^{jwt_k})^T nu_k
```

`b*e^{jwt_k}` IS this engine family's own per-step forced-column convention
(the same discretisation a stationary source's own `br`/`bi` get phase-
modulated through at every step) -- so `nu_k` itself is exactly
`dY0(w)/d(beta_k)`, the sensitivity of the output sideband-0 response to a
UNIT, UNMODULATED forcing injected into the per-step RHS at step `k` ONLY:
the discrete-time periodic transfer function `H(t_k,f)` this card's own
brief calls for ("a noise current injected over step k has the same stamp
as a small-signal source at step k") -- the SAME per-step RHS convention as
a stationary AC source, just not phase-modulated across steps (a
genuinely LOCALISED injection at one step, not a continuous sinusoid).

For a noise generator `g` with a FIXED (bias-independent) terminal
incidence pattern `b_g` (e.g. drain-minus-source for a MOS channel current
-- exactly `NevalSrc()`'s own posNode/negNode convention) and a
cyclostationary density `S_g(t_k,f)` (the device's own `DEVnoise` routine,
evaluated at the step-`k` operating point), the Telichevesky/Kundert
(1996, "efficient AC and noise analysis of two-tone RF circuits")/Demir-
Roychowdhury cyclostationary-noise integral

```
onoise(f) = (1/T) INT_0^T |H(t,f)|^2 S(t,f) dt
```

discretises, on this engine's own periodic-trapezoidal grid, to

```
onoise(f) = T * Sum_k (1/w_k) * Sum_g |b_g . nu_k(f)|^2 S_g(t_k,f)
```

`w_k = 0.5*(h_k+h_{k+1 mod N})` is the SAME periodic-trapezoidal weight the
boundary term's own `coef` uses for the OUTPUT side in `pac_engine_run_
adjoint()` -- used here, DIVIDED rather than multiplied, for the INJECTION
side. Continuous white noise cannot be "sampled" as `beta_k = i_n(t_k)`
directly (a Dirac-delta autocorrelation has no well-defined pointwise
sample); the physically consistent discretisation (the same one ngspice's
own `trnoise` white generator and the TK96/Demir-Roychowdhury discrete-
noise literature use) represents `i_n(t)` over step `k` as a ZERO-ORDER-
HOLD value `n_k` entering the per-step recursion exactly like a
deterministic `beta_k`, with `Var(n_k) = S_k/w_k` so that the CHARGE
`n_k*w_k` over the step reproduces the true continuous process's own charge
variance `S_k*w_k` (a white process integrated over an interval of length
`w_k` has variance `S_k*w_k`). Since the `n_k` are independent across `k`
(white/uncorrelated), the output sideband-0 response `Y0_random(f) =
Sum_k nu_k(f)*n_k` has variance `Sum_k |nu_k(f)|^2 * Var(n_k) = Sum_k
S_k*|nu_k(f)|^2 / w_k` -- `S_g*|nu_k.b_g|^2` is `dens_k` below, and the
remaining `T` is this card's own empirically-confirmed normalisation (see
"A derivation false start" below for how it was found).

`Sum_g |b_g.nu_k(f)|^2 S_g(t_k,f)` is EXACTLY what `NevalSrc()`/`DEVnoise()`
already compute when `ckt->CKTrhs`/`CKTirhs` are loaded with `nu_k(f)` at
the step-`k` bias (`NevalSrc`: `gain = |CKTrhs[node1]-CKTrhs[node2]|^2`,
`*noise = gain*S(bias)`) -- so this engine needs **no new device-noise
physics at all**: it reuses `pac_engine_run_adjoint()`'s own backward sweep
(`Phi`/`G(w)`/`M(w)` passes, and the per-step `abarS`/`abarT`/`Cdiag`/
`beta`/`gamma` recursion -- copied VERBATIM into a new function,
`pac_engine_run_adjoint_noise()`, not shared by a refactor, so this card
touches ZERO lines of the already-validated `.pxf` adjoint engine) and, at
each `(k,f)`, loads `nu_k` into `ckt->CKTrhs`/`CKTirhs` with the bias
ALREADY at `x_k` (`pac_adj_refactor_at()` -- the SAME `CKTload()` that
rebuilds `J_k` also refreshes every device's own bias-dependent noise
parameters, e.g. `RESconduct`/`REScurrent`/`gm`) and calls the ordinary
`DEVnoise(N_DENS,N_CALC,...)` loop -- exactly the idiom `dcpss.c`'s own
frequency-domain `pnoise_sweep()`/`QPnoiseAnalyze()` already use for their
own (different) adjoint vector.

**Scope (v1)**: restricted to resistor thermal, diode shot, and MOS1/BSIM4
channel thermal+flicker noise -- the only `DEVnoise` routines this
derivation has been checked against. `pac_noise_check_devices()`
(`pacshoot.c`) scans the circuit for any OTHER device type with
noise-producing instances (matched by `DEVnoise` FUNCTION POINTER, e.g.
`RESnoise`/`DIOnoise`/`MOS1noise`/`BSIM4noise`, not by device-name string)
and REFUSES loudly (`E_UNSUPP`, naming the actual device type) rather than
silently omitting or mis-scaling its contribution. Saltation (ideal-switch)
events and free-running orbits are refused too (`PAC_ADJOINT_REFUSED`,
checked identically to `.pxf adjoint=1` -- `PSSsaltEventsSeen()>0` after
the `Phi` pass, or `trig_idx>=0`) with **no fallback engine** -- unlike
`.pxf`'s `adjoint=0`, there is no meaningful "one pass per named source"
analogue for noise, which comes continuously from every device, not a
short list of named sources. Input-referred (`inoise`) is not computed.

### A derivation false start (found, and fixed, validating this card)

Two OTHER candidate per-step weightings were tried, in order, before the
one above, both found WRONG by `check_pnoise_lti()` below (the LTI-vs-
plain-`.noise` check) during this card's own development:

1. **No extra weight at all** (`onoise(f) = Sum_k dens_k`, `dens_k =
   S_k*|nu_k|^2`, unweighted): under-reported by almost exactly a factor of
   `N` (the accepted-step count) at three grid resolutions
   (Points=200/2000/20000) -- the value SHRANK roughly as `1/N` as the grid
   was refined, the opposite of convergence, immediately flagging it wrong
   (a correct discretisation of a fixed integral must NOT vanish as the
   grid refines).
2. **Multiply by `w_k/T`** (mirroring the output tap's own `coef` weight,
   by direct but wrong analogy -- "the input side should look like the
   output side"): under-reported by almost exactly `N^2`.

The correct formula, `T * dens_k/w_k` (divide, not multiply, by the
quadrature weight; multiply, not divide, by the period) was found by fitting
several candidate weight functions against the SAME per-step
`(t_k,h_k,dens_k)` data (a debug dump, `PNOISE_DEBUG_STEPS=1`, still in
`pacshoot.c`, gated off by default) directly in Python until one matched
`.noise`'s own reference value to `<1e-5` relative -- confirmed
analytically afterward (see "Design" above, the `Var(n_k)=S_k/w_k`
white-noise-sample argument) once the right answer was known, not derived
purely on paper first. This project's OWN working rule ("done needs the
read-back... works needs the test output") is exactly why this was caught
before being shipped wrong: the LTI-vs-`.noise` check was written and run
BEFORE any of the other fixtures, specifically because it has an
independent, trusted, already-validated reference (plain `.noise`) to catch
exactly this class of error.

### Validation (card #509)

`tests/test_pnoise.py`:

| check | fixture | result |
|---|---|---|
| `pnoise_lti` | `pxf_lti.cir` (LTI: Rin/Rload/Cout, unrelated clock for PSS) vs plain `.noise`, `onoise_spectrum`, 1kHz-1MHz | Points=200: 7.9e-3 worst rel; Points=2000: 7.9e-5; **Points=20000: 4.3e-7** (target 1e-6, met) -- convergence is clean and monotonic across all three, not a lucky single point |
| `pnoise_r3` | `pac_r3_ctrl2out.cir` (R3 voltage-mode buck, smooth switches, `dev=1`) | runs to completion, every `onoise_spectrum` value finite and >=0 (range 1.3e-8 to 5.3e-8 V/rtHz, 500Hz-50kHz), per-device sum-of-squares vs the reported total: worst relative **8.2e-16** (machine precision) |
| `pnoise_r2_smooth` | `pnoise_r2_smooth.cir` (switched-RC LPTV, SMOOTH tanh switch -- see below) `dev=1` | (Monte Carlo cross-check in the next section) runs to completion, finite/>=0 (range 7.9e-10 to 4.4e-9 V/rtHz, 1kHz-1MHz), per-device sum-of-squares worst relative **1.4e-14** |
| `pnoise_refuse` | `pxf_r2_switched_rc.cir` (ideal `sw`) | REFUSED loudly (device-scope gate: the `Switch` model's own `swnoise.c` DEVnoise is outside v1 scope), exit 0, no crash/hang |

`test_pxf.py`, `test_pac.py`, `test_pss.py`, `test_stb_native.py`,
`test_checkpoint.py`, `test_sens.py` re-run against this card's build: no
change from the s55 baseline (`test_pxf.py`'s two `ASan: frequency-domain
... harmonics=1` FAILs and `test_pac.py`'s R5-free-running-above-band
FAIL+COULD-NOT-RUN are ALL pre-existing -- confirmed identical, byte-for-
byte same failure lines, running the SAME checks against the unmodified
`aaba50ba3f49-s55` baseline binary) -- this feature is purely additive (two
new functions in `pacshoot.c`, one new field on `PSSan`, new trailing
dot-card tokens; zero lines changed in `pac_engine_run_adjoint()`,
`PSSpacRecordStep()`, or any other already-shipped `.pac`/`.pstb`/`.pxf`
code path). ASan build (`buildrecipe/README.md`'s recipe,
`aaba50ba3f49-pnoise-asan`): all four `test_pnoise.py` checks clean, no
AddressSanitizer report.

**Switched-RC fixture note**: `pxf_r2_switched_rc.cir` (the existing
`.pxf`/`.pac` switched-RC fixture) uses an IDEAL `sw` -- this engine
REFUSES saltation with NO fallback (see "Scope" above), so
`pnoise_r2_smooth.cir` (new this card) replaces the ideal switch with a
SMOOTH tanh-based behavioural conductance (same idiom `pac_r3_ctrl2out.cir`
already uses for its own power-stage switches), same `Rin`-then-switched-
`Rsw`-to-ground topology otherwise. Component values needed retuning from
`pxf_r2_switched_rc.cir`'s own (`rin=rsw=10k`) -- PSS-SHOOT rejected those
as a "degenerate all-off fixed point" (not a real orbit) with the smooth
switch; `rsw=100` (well below `rin=10k`, giving a clean, well-converging
ON/OFF swing) converges.

### What was skipped

- **(Resolved, second pass -- see the next section.)** The independent analytic and transient-noise Monte Carlo cross-checks that the first pass skipped were done in the second pass; no formulation change was needed.
- **The saltation-refusal path itself** (`PSSsaltEventsSeen()>0`, as
  opposed to the device-scope refusal `pnoise_refuse` above actually hits
  first on `pxf_r2_switched_rc.cir`) was not separately exercised by a
  fixture this pass -- it is a near-verbatim copy of the already-validated
  `.pxf adjoint=1` code (unchanged logic, same check), so risk is judged
  low, but it is untested BY THIS CARD's own suite.
- **Sidebands other than 0** -- same sideband-0-only scope as `.pxf`/`.pac
  solver=td`.
- **Input-referred (`inoise`) spectrum** -- not computed by `solver=td`;
  `InSrc` is accepted syntactically but ignored. The existing frequency-
  domain `.pnoise` still reports it.
- **A device-scope-refusal fixture with a device OTHER than `Switch`** (a
  BJT, JFET, or MESFET circuit, say) was not built -- `pxf_r2_switched_rc.
  cir`'s ideal `sw` happens to ALSO be outside scope (via `swnoise.c`'s own
  `DEVnoise`), so `pnoise_refuse` exercises the device-scope gate at all,
  just not on a device this project's other fixtures don't already avoid
  for unrelated reasons.
- **Flicker (1/f) noise's own frequency dependence was not separately
  isolated** -- the LTI check (thermal-noise-only resistors) does not
  exercise `data->freq`-dependent densities at all; `pac_r3_ctrl2out.cir`'s
  resistors are also purely thermal (no `Kf`/`Af` given). The mechanism
  (`nd.freq = freqs[f]`, set once per frequency, identically to
  `pnoise_sweep()`'s own convention) is the SAME code path plain `.noise`
  already validates for flicker noise generally, so risk is judged low but
  untested BY THIS CARD's own suite specifically.
- **Cyclostationary flicker/shot noise where `S` itself has REAL frequency
  structure that interacts with the folding** (Enhancement-177/178's own
  "noise originates at a DIFFERENT physical frequency than the output"
  concern, which the existing frequency-domain engine's `cyclo` mode
  handles explicitly) -- this engine evaluates `S_g(t_k,f)` at the OUTPUT
  frequency `f` directly (matching a "modulated-white" source exactly, the
  same simplifying assumption `pnoise_sweep()`'s OWN non-`cyclo` default
  makes), not the sideband-folded physical frequency E-177 derived for the
  frequency-domain engine's OWN harmonic sidebands. This time-domain engine
  has no discrete "sideband" concept to fold across (it IS the time
  domain), so E-177's specific correction does not directly translate --
  whether an analogous correction is needed here (and what it would look
  like) was not investigated.

## `.options klu` with the shooting engines (card #533, 2026-10-01)

`.pss ... method=shoot` (driven and `trigger=` free-running), `.pac`/`.pstb`/`.pxf`/`.pnoise ... solver="td"`
(including `.pxf adjoint=1`, saltation decks, the transposed solve) run on the KLU matrix when `.options klu`
is set. Before this card the whole family refused it ("periodic steady state analysis (.pss) is not supported
with 'option KLU'"), so converters with extracted parasitics (thousands of R/C) were stuck on SPARSE 1.3 and got
none of the KLU work of rounds 13-15 (unscaled refactor, R/C base cache, flat capacitor load).

**Still refused under KLU:** only the Lannutti PSS (`.pss` without `method=shoot`, `method=gear` etc.), with the
message "... not supported with 'option KLU' unless method=shoot is given ...". Nothing else in the family
refuses. (Not tested under KLU: `.pnoise` with `dev=`/per-device tables beyond `pnoise_r2_smooth`, the
`tools/pac_perturbed.py` flow, XSPICE code models, BSIM4-heavy PSS decks.)

### What had to change (3 things; the engine itself used only SMP* calls)

1. `dcpss.c`: the refusal moved after the `method=shoot` dispatch (it was in front of it).
2. `pacshoot.c` `pac_capture_b0()`: it puts the matrix in complex stamp mode for `CKTacLoad(w=1)` (to capture the
   AC source vector B0) and never put it back. SPARSE does not care; under KLU every device's stamp pointers are
   re-bound to the complex CSC array (`DEVbindCSCComplex`), so the first real transient after it stamped into the
   wrong array ("singular matrix: check node vin#branch", PAC-TD FAILED). It now converts back
   (`DEVbindCSCComplexToReal`, `KLUmatrixIsComplex = KLUmatrixReal`) on every exit path, as `acan.c` does at the
   end of `.ac`.
3. `klusmp.c` `SMPluFac()` + `dcpss.c pss_shoot_klu_guard()`: a pivot-growth guard while a shooting engine runs.
   `klu_refactor` reuses the pivot order of the last `klu_factor` with no threshold test; SPARSE's `spFactor`
   re-pivots a small pivot. On the stiff R4 peak-current-mode deck (latch node `q`) the reused order degrades to a
   reciprocal pivot growth of 1e-18 and the shooting warm-up stops with "Timestep too small ... node q" (plain
   `.tran` of the same deck is fine under KLU; which of the R4 `.pac` fixtures failed changed with the rounding
   of the refactor variant, i.e. it is marginal, not deck-specific). Every 8th refactor (`NGSPICE_KLU_RGEVERY`)
   `klu_rgrowth` is evaluated; below 1e-10 (`NGSPICE_KLU_RGTHR`; 0 disables) `SMPluFac` returns `E_SINGULAR`,
   which NIiter already answers with a full re-pivoting `klu_factor` of the same matrix. The count is printed:
   `PSS-SHOOT: KLU pivot-growth guard re-pivoted N refactorization(s)` (R4 fixtures: 16-30, R3 and the PEX deck: 0).
   Cost: checking every refactor was ~10 % of a PEX-size PSS run, every 8th is ~1.4 %.
   `NGSPICE_PSS_KLU_SCALED=1` additionally keeps KLU's row-scaled refactor (round-13 unscaled refactor is the
   default); the scaled refactor costs +15 % on the PEX deck and is not needed once the guard is on.

The transposed real solve (`SMPtSolve`, card #504) already had a KLU path (`klu_tsolve`); `.pxf adjoint=1` under KLU
uses it unchanged. SPARSE runs are byte-identical to the s85 build (22 decks, `cmp` of the wrdata files).

### KLU vs SPARSE, same deck, same binary (`tests/test_pss_klu.py`, outputs in `~/Work/data/ngspice_augmentation/533/eq`)

| case | KLU vs SPARSE max deviation |
|---|---|
| pss_r3_shoot / r4 / r6s / r2 / r2_q4 orbit (smooth vectors, worst, /max abs) | 1.4e-6 / 3.6e-6 / 1.6e-9 / 1.8e-13 / 1.8e-13 (R4 v(vc): slow Floquet mode rho 0.994 x residual stop 1e-6) |
| orbit edge vectors (pwm, ramp, switch node; edge-tolerant, /p-p) | 7.6e-3 (R3 pwm), 1.2e-3 (R4 sw) -- both runs chord the edge on their own adaptive grids (R5 free-running: edge crossing times agree to 3 ps) |
| pac_r3 ctrl2out / line2out / zout (max abs dH / max abs H) | 2.9e-6 / 7.2e-6 / 8.7e-6 (0.074 mdB, 0.44 mdeg worst point) |
| pac_r4 ctrl2out / line2out / zout | 7.8e-8 / 4.6e-8 / 1.3e-10 |
| pstb_r3 (T) | 5.8e-6 |
| pxf_lti adjoint=0 / adjoint=1 | 2e-13 / 2e-13 |
| pxf_r3 (Vin, Vctrl) adjoint=0 / adjoint=1 | 7.2e-6 / 7.2e-6 (adjoint engine confirmed running under KLU) |
| pxf_r2 (saltation; adjoint=1 falls back to adjoint=0) | 5e-14 |
| pnoise_r2_smooth | 6e-14 |
| pss_r5 free-running (trigger=; last 5 % of the record left out, it ends at the trigger event), pss_r6o (ideal switches) | 2.2e-6 smooth (edge 6e-4) / 2e-12 |

The ~1e-5 level of the converter decks is the PSS stop criterion (residual 1e-6, the two solvers stop at
different Newton iterates), not a solver difference; the LTI and switched-RC decks (exact orbits) agree to 1e-13.

### Scale: R3 buck + three extracted RC trees, 2418 parasitic R/C (1150 capacitor states), `.pac` after PSS

R3 voltage-mode buck (`pac_r3_ctrl2out.cir`, corrected `Goa 0 opout vrefnode compin`) + `bench/gen_pex_buck.py`
trees: decoupling net on vout (10 ohm), switch-node net on sw (0.5 ohm), ground-return net (cap return `gret`,
1 ohm); `.save` limits the per-run transient plots (else every parasitic node of the 1150 FD-monodromy runs is
kept: >8 GB). One ngspice at a time, interleaved, start position rotated, OMP_NUM_THREADS=2, development workstation (AMD Ryzen AI Max+ 395, Linux) loadavg 5.1-7.3.

| arm (n=4) | wall median s | PSS-SHOOT s | PAC-TD s | min-max wall |
|---|---|---|---|---|
| SPARSE | 242.9 | 229 | 12.6-15.8 | 228.9-245.7 |
| KLU (this card) | 173.7 (-28.5 %) | 158.7 (-31 %) | 13.5-15.8 | 169.9-177.1 |

PSS and the one `.pac` sweep converge in 3 Newton iterations on both; transfer curves agree to 3.8e-6 (max/max) /
1e-4 dB / 0.6 mdeg. PSS gets the same -25 % that plain `.tran` gets on this deck (0.48 s vs 0.37 s per 100 us);
the PAC recording pass is not faster (about +2 s under KLU, not understood, within the noise of the load).
Earlier arms of the same A/B (build with guard on every refactor, 3 reps): KLU 192.3 s, KLU + scaled refactor
222.5 s, SPARSE 239.0 s.

Limits found: the shooting FD monodromy runs one transient per capacitor/inductor state, so wall time and, without
`.save`, memory grow ~quadratically with the number of parasitic capacitors (b3d4: 490 states, SPARSE 49 s / KLU
35 s; 1150 states: 243 s / 174 s; 5700 states not attempted). Fast-state elimination (1150 -> 5 Newton unknowns) is what
keeps the Newton phase cheap; the one-time detection sweep is the cost.

## Matrix-free Newton-GMRES shooting, `krylov=1` (card #537, 2026-10-01)

`.pss ... method=shoot` builds the monodromy M of the one-period map column by column: two perturbed
period runs (centred difference) per retained state, per Newton iteration, plus the one-time fast-state
detection sweep of 2n+ period runs. That is fine for 4-20 states and quadratic in nothing but painful at
PEX scale (R3 buck + 2418 parasitic R/C = 1150 capacitor states: 2348 period runs, 152 s on KLU). `krylov=1`
replaces it with the matrix-free Newton-GMRES method of Telichevesky, Kundert and White (DAC 1995): nothing of
M is assembled, no fast-state detection runs, every state is a Newton unknown, and the Newton step
`(M - I) dx = -(P(x0) - x0)` is solved by GMRES where each Krylov vector costs one finite-difference
directional derivative, i.e. one pair of perturbed period runs. Parasitic states decay inside the period (their
Floquet multipliers are about 0), so `M - I` has a handful of eigenvalues away from -1 and GMRES needs about
(number of slow modes + 1) products: 5-6 products on every deck measured, from 7 to 2470 states.

### Interface

| | |
|---|---|
| card key | `krylov=1` on `.pss`, and on `.pac`, `.pstb`, `.pxf`, `.pnoise` (their embedded PSS); control form: `pss ... method=shoot krylov=1` |
| environment | `NGSPICE_PSS_KRYLOV=1` forces it on, `=0` forces it off (the environment wins over the card) |
| default | off. With it off the engine is byte-identical to the previous build (12 fixtures: stdout with `wall_time` masked, `wrdata` and `.raw`) |
| `krylov=1` without `method=shoot` | ignored, with a note on stderr |
| free-running (`trigger=`) | not supported (the Poincare-section map carries the period as an unknown): the run uses the dense monodromy and says so on stderr (`PSS-FREERUN: krylov=1 is not supported ...`) |
| `PSS_SHOOT_JACOBIAN=analytic` | wins over `krylov=1` (a note is printed) |

Tunables (debug/experiment environment variables, same style as the other `PSS_SHOOT_*`): `PSS_KRYLOV_ETA` (inexact-Newton
tolerance, default 1e-3), `PSS_KRYLOV_MAXDIM` (GMRES dimension limit, default 40, no restart), `PSS_KRYLOV_ARNOLDI` (largest
Arnoldi dimension for the stability estimate, default 40), `PSS_KRYLOV_ONESIDED=1` (one-sided products, see below),
`PSS_KRYLOV_DEBUG=1` (GMRES residual history), `PSS_SHOOT_COUNT_RUNS=1` (print the period-run count also in dense mode).

### What is shared with the dense engine (so F and the products are one map)

Only the linear solve for the Newton direction is replaced. The baseline `P(x0)` block (recorded step grid, forced rerun on
switch decks, ripple measurement), every perturbed run (`pss_run_period_ex_safe(..., force=1)`: same recorded grid or exact
replay as the dense FD columns and the damped-step trials), the ripple-relative `scale[k] = max(ripple_k, 1e-2)` that
defines the residual, the trust-radius cap, the damped line search on dx, the degenerate-all-off-fixed-point guard, the staged
warm-up and the replay/nudge ladder are the existing code. The line search and the convergence test always use the true
residual `max_k |P(x0)_k - x0_k| / scale_k <= steady_coeff`, so an inexact Krylov direction cannot produce a false
`CONVERGED`.

### The linear solve

* GMRES in the weighted variables `y = x / scale` (the weights of the residual test), x0 = 0, modified Gram-Schmidt with one
  re-orthogonalisation, Givens least squares. Stops at relative residual `eta = 1e-3` (inexact Newton), at the dimension limit, or
  at breakdown. Relative residual still above 0.5 at the limit = stagnation.
* Product `M v`: v (unit norm in the weighted variables) maps to the physical direction `d_k = scale_k v_k`; the perturbation
  `p = sigma d` has sigma chosen so that `max_k |p_k| / h_k = 1`, where `h_k` is the dense engine's own per-state FD step (`|x|*1e-5`,
  floor 1e-6, 1 % of the ripple for a state smaller than its ripple): no state is moved by more than its dense FD step.
  Centred difference `(P(x0+p) - P(x0-p)) / 2`. A product whose period run fails retries at sigma/4 (up to 5 times).
* Failure handling: a product that cannot be run, GMRES stagnation, or (final stage only) a Newton attempt that does not
  converge makes the analysis continue **on the dense path** from the current x0 (fast-state detection, dense FD monodromy), exactly as
  without `krylov=1`; the output says so (`... continuing on the dense path ...`, and the run-count line ends `mode=krylov-fell-back-to-dense`).
  A probe stage that merely runs out of its short iteration budget stays on the Krylov path (the normal warm-up escalation).
* One-sided products (`PSS_KRYLOV_ONESIDED=1`, half the runs) are not the default: on R3 they converge in 43 instead of 116 runs, but
  on R4/R3d (rho 0.994) the direction error is amplified 1/(1-rho) = 166 times and the Newton fails or the Arnoldi radius is wrong
  (R3d: reported "UNSTABLE 1.275" for a stable orbit), the same reason the dense engine uses centred differences.

### Stability report (no M exists)

After convergence an Arnoldi process (fixed pseudo-random start vector, grown in steps of 4/8 until the spectral radius settles,
at most 40 products) gives the Hessenberg matrix whose eigenvalues approximate the leading Floquet multipliers:

```
PSS-SHOOT: Arnoldi estimate of the Floquet multipliers: dimension 12 of 1150, 24 period runs; spectral radius 0.980318 (0.980318 at dimension 8)
PSS-SHOOT: CONVERGED iter=2 residual=6.099072e-08 spectral_radius=0.980318 stable[arnoldi] wall_time=7.977s
PSS-SHOOT: FLOQUET dominant_multiplier=0.980318+0i angle_deg=0.000 (Arnoldi estimate)
PSS-SHOOT: period runs=75 (+ 140 warm-up transient periods), mode=krylov
```

The verdict token carries the marker (`stable[arnoldi]`, `UNSTABLE(period-doubling)[arnoldi]`), so a parser that matches
`spectral_radius=(\S+) (\S+) wall_time=` still works and a verdict check sees which engine measured it. If the radius has not settled to 1 % at the
maximum dimension, or the interval between the last two estimates straddles 1, or the eigen-solve fails, the verdict is
`UNKNOWN(not-computed)`, the radius is printed as -1 and the line above says why. When the Arnoldi dimension reaches the
number of states (every deck up to 8 states) the Hessenberg spectrum is the whole spectrum. Accuracy against the dense FD
Floquet (both are finite-difference estimates, about 1e-3 apart on R3/R4/R3d): 0.980210/0.980199 (R3), 0.994/0.994 (R4),
0.980278/0.980205 (R6s), 0.9801/0.9795 (R3d), 0.9944/0.9940 (R4d), 1.1910/1.1910 (R3 x20, via dense fallback).

### `.pac` / `.pstb` / `.pxf` / `.pnoise` after a Krylov PSS

They take only the orbit state x0 (`pac_pstb_run(ckt, job, x0, n, T)`) and build their own monodromy with the per-step
propagated recording pass; they do not use anything from the PSS Newton matrix. Checked: `.pac`, `.pstb`, `.pxf adjoint=0`
and `.pnoise` on R3 / R2 after `krylov=1` agree with the dense-PSS curves to 1.6e-7, 2.5e-7, 1.6e-7 and 0 (max |dH| / max |H|).
The PEX-scale `.pac` itself (PAC-TD, 1150 states, ~14 s) is unchanged.

### Results (KLU, OMP_NUM_THREADS=2, one ngspice at a time, development workstation load 5.5-6.0, alternating order, `tests/test_pss_krylov.py`)

PEX-size converter, `.pss` only, R3 buck + three extracted RC trees (`bench/gen_pex_buck.py`, `.save`, `.options klu`), orbit
deviation = max |dev| / max(peak-peak, 1e-2) over v(vout), v(sw), v(opout), i(L1):

| states | dense runs / wall | krylov runs / wall | runs ratio | wall ratio | orbit dev | rho dense / Arnoldi |
|---|---|---|---|---|---|---|
| 37 | 122 / 0.7 s | 90 / 0.7 s | 1.4x | 1.0x | | 0.9805 / 0.9803 |
| 73 | 194 / 1.5 s | 90 / 1.1 s | 2.2x | 1.4x | | 0.9798 / 0.9803 |
| 145 | 338 / 4.0 s | 88 / 1.7 s | 3.8x | 2.4x | | 0.9810 / 0.9810 |
| 205 | 458 / 6.5 s | 77 / 1.7 s | 5.9x | 3.8x | | 0.9803 / 0.9804 |
| 490 | 1028 / 31.8 s | 77 / 4.1 s | 13.4x | 7.8x | 2.2e-6 | 0.9803 / 0.9803 |
| 1150 | 2348 / 152.1 s | 75 / 8.0 s | 31.3x | 19.0x | 5.3e-7 | 0.9803 / 0.9803 |
| 2470 | not run (extrapolated > 10 min) | 86 / 21.0 s | | | | / 0.9803 |

(1150 states, final build -4, loadavg 4.5-5.2: dense wall 152.1 s and 151.4 s, Krylov 8.0 s and 7.9 s, two alternations (dense-krylov, krylov-dense); the dense 2348 runs are the 2300-run
detection sweep plus the Newton stages. The Krylov run is 140 periods of warm-up transient + 75 period runs; the warm-up is
now the larger part of the 8 s.)

Small fixtures (`tests/fixtures`), period runs dense / krylov and verdict:

| deck | dense runs | krylov runs | verdict dense / krylov | note |
|---|---|---|---|---|
| pss_r2_shoot (1 state) | 13 | 7 | stable 0.2231 / stable[arnoldi] 0.2231 | |
| pss_r2_q4_shoot | 13 | 7 | stable 0.0025 / 0.0025 | |
| pss_r3_shoot (7) | 62 | 116 | stable 0.9802 / 0.9802 | small deck: Krylov costs more |
| pss_r4_shoot (7) | 89 | 127 | stable 0.9940 / 0.9940 | |
| pss_r6s_shoot (8, ideal sw) | 86 | 68 | stable 0.9802 / 0.9803 | |
| rd_r3d_shoot (9) | 276 | 95 | stable 0.9795 / 0.9801 | |
| rd_r4d_shoot (9) | 132 | 195 | stable 0.9940 / 0.9944 | |
| pss_r3_gainx20_unstable | 148 | 281 | UNSTABLE 1.191 / same | Krylov attempt did not converge, dense fallback |
| pss_r4_d55_unstable | 152 | 228 | UNSTABLE(period-doubling) 1.3136 / 1.3084 | same, dense fallback |
| pss_r5_shoot (free-running) | | | | refused, dense, identical result |

Orbit equivalence (`N` = max |dev| / max(p-p, 1e-2), the engine's own ripple-relative norm): against the **dense run with no fast-state elimination**
(the dense variant with the same unknowns as Krylov) the orbit agrees to 1.0e-6 (R3), 2.4e-5 (R4, rho 0.994), 1.8e-6 (R6s), 1.0e-6 (R3d), 8.2e-5 (R4d); against the
default dense run the difference is up to 7.2e-4 (R6s) and 2e-2 (R3d), and it is the **dense engine's own spread**: its variants (different FD step, no
elimination, warm-up 16 periods) differ from it by 7.2e-4 and 2.1e-2 on the same two decks. The cause is the period map's dependence on the recorded
step grid on switch / dead-time decks (R6s: converging at warm-up stage 3 or stage 4 moves v(opout) by 9e-6 V; a converged residual of 1e-6 of a slow
Floquet mode, rho 0.994, is worth 1/(1-rho) = 166 in the orbit). It is not a Krylov effect; the equivalence test uses the dense variants' spread as the
tolerance floor.

### Automatic threshold (proposed, NOT enabled)

The Krylov cost is flat (75-90 period runs from 37 to 2470 states: warm-up probes, 3-4 Newton iterations of 10-12 product runs, and the Arnoldi estimate) while the
dense cost is 2n+ detection runs plus Newton: the run counts cross at about 25-30 states (R3 7 states: 62 / 116; 37 states: 122 / 90), the wall time at
about 40 (0.7 s each). Proposal: enable `krylov` automatically when the number of capacitor + inductor states (`PSS-SHOOT: ... N state(s)`, known before any
period run) is at least 64 (2.2x fewer runs, 1.4x less wall at 73 states; 5x at 205). Not enabled because (a) the unstable/non-converging decks fall back to dense after a failed
Krylov attempt (cost: the failed attempt, e.g. 281 vs 148 runs on R3 x20), (b) a deck with many slow modes (several tens of non-negligible
multipliers) has a larger GMRES dimension and a smaller advantage (not measured: every deck here has 5-6 slow modes), (c) changing the defaults is a separate project decision.

### Not done / limits (2026-10-01)

* Free-running (`trigger=`) shooting is not matrix-free (refused, runs dense; its states are few).
* The warm-up transient (140 periods here) and the baseline/forced reruns are unchanged and are now the bulk of the 8 s on the 1150-state deck.
* No restart in GMRES (limit 40, then dense fallback), no preconditioner, no Eisenstat-Walker forcing term (fixed eta 1e-3), no Krylov-subspace recycling between Newton iterations.
* `krylov=1` with `PSS_SHOOT_JACOBIAN=analytic` is not combined (analytic wins).
* The Arnoldi estimate converges to the leading multipliers; an isolated eigenvalue near the unit circle that the pseudo-random start vector excites
  weakly could be missed in principle (the settle test compares dimension m and m-4 only).
* Not tested: Krylov on BSIM4-heavy decks, XSPICE-clocked (`driven=1`) decks, `.pac` with `pacpoints`, PEX decks with a slow parasitic mode (rho of a parasitic near 1), more than 2470 states.
* Dense baseline for 2470 states not run (would be > 10 min, memory not checked).

## Automatic matrix-free shooting, warm-up cost, free-running `trigger=` (card #538, 2026-10-02)

Series s88 + branch `aug-538` (commits d19ee52d4 and 6285624ca on top of aug-537); build `aaba50ba3f49-538-10`
(ASan: `aaba50ba3f49-538-7-asan`, rebuilt from the final source). All numbers: KLU, OMP_NUM_THREADS=2, development workstation shared (loadavg 18-23
during the timed runs; arms interleaved, CPU time recorded next to wall time), decks and logs in
`~/Work/data/ngspice_augmentation/538/`.

### 1. Automatic `krylov` (default on from 64 states)

`krylov` is chosen automatically when the number of capacitor + inductor states (`PSS-SHOOT: ... N state(s)`, known before any
period run) is **>= 64**, for `.pss method=shoot` and for the PSS embedded in `.pac/.pstb/.pxf/.pnoise`, and for free-running `trigger=` (below).
One line on stdout when auto picks it:

```
PSS-SHOOT: 73 states >= 64: matrix-free Newton-GMRES (krylov) chosen automatically (krylov=0 or NGSPICE_PSS_KRYLOV=0 forces the dense monodromy)
```

| | |
|---|---|
| force dense | `krylov=0` on the card, or `NGSPICE_PSS_KRYLOV=0` (the environment wins over the card) |
| force matrix-free | `krylov=1` on the card, or `NGSPICE_PSS_KRYLOV=1`; both work below 64 states as in card #537 |
| threshold | `PSS_KRYLOV_AUTO_MIN` (default 64; `0` switches the automatic choice off; the boundary is `n >= AUTO_MIN`) |
| not automatic | `PSS_SHOOT_JACOBIAN=analytic`; free-running with `.pac/.pstb` on the orbit; whenever the card or the environment name a choice |
| below 64 states | untouched: all 12 fixtures byte-identical to s88 (`tests/test_pss_krylov.py identical`, stdout with `wall_time` masked, `wrdata` and `.raw`) |

**State-count sweep** (R3 buck + extracted RC trees, `bench/gen_pex_buck.py` trees of branch x depth, 3 trees, `.save`, 2 interleaved
repetitions, runs are deterministic; "krylov" below 64 = `krylov=1` forced, from 64 = automatic):

| states | dense runs | dense wall s | krylov runs | krylov wall s | wall ratio |
|---|---|---|---|---|---|
| 16 | 80 | 0.58 | 84 | 0.82 | 0.7x |
| 19 | 86 | 0.61 | 86 | 0.88 | 0.7x |
| 22 | 92 | 0.70 | 86 | 0.95 | 0.7x |
| 25 | 98 | 0.72 | 71 | 0.77 | 0.9x |
| 31 | 110 | 0.83 | 73 | 0.84 | 1.0x |
| 34 | 116 | 0.98 | 86 | 1.06 | 0.9x |
| 37 | 122 | 0.99 | 86 | 1.09 | 0.9x |
| 46 | 140 | 1.31 | 73 | 1.04 | 1.3x |
| 58 | 164 | 1.72 | 86 | 1.41 | 1.2x |
| **64** | 176 | 1.92 | 47 | 0.90 | 2.1x |
| 70 | 188 | 2.08 | 45 | 0.83 | 2.5x |
| 82 | 212 | 2.71 | 60 | 1.31 | 2.1x |
| 100 | 248 | 3.56 | 47 | 1.20 | 3.0x |
| 115 | 278 | 3.98 | 47 | 1.19 | 3.3x |
| 145 | 338 | 6.04 | 58 | 1.94 | 3.1x |
| 172 | 392 | 7.04 | 60 | 1.99 | 3.5x |
| 205 | 458 | 10.08 | 47 | 1.88 | 5.4x |

Run counts cross at about 20-25 states, wall time at about 40 (R3 class). The slow-mode / dead-time classes (PEX trees on the R4, R3d, R4d, R3 fixtures,
`krylov=1` forced, 21-53 states) are slower matrix-free than dense in wall time almost everywhere (0.5-1.1x: R4 21/27/51 states 0.70/0.73/0.85x, R3d 23/29/53 states
0.96/0.74/1.03x, R4d 23/29/53 states 0.81/0.47/1.08x -- at 29 states the Krylov attempt fails and falls back to dense, 348 vs 181 runs -- R3 fixture 21/27 states 0.6x), and win at
99-101 states (R4 1.3x, R3d 1.2x, R4d 1.4x). The threshold stays at **64**: below it the dense engine costs under 4 s and has no Krylov-failure mode; from 64 up the advantage is
2x on R3-class decks and growing. `PSS_KRYLOV_AUTO_MIN` moves it.

### 2. Where the 8 s of the 1150-state deck go, and what was changed

Profile (`PSS_SHOOT_PROFILE=1`, s88 build, 1150 states, loadavg ~20, total 11.6 s; 8.0 s at lighter load):

| part | s | share |
|---|---|---|
| warm-up transients (140 periods, probe stages, 10x looser step) | 0.45 | 4 % |
| stage 1 + 2 Newton attempts (both abandoned: residual 1.0, step clipped to alpha 7e-6 / 0.19; 10 + 12 product runs + baselines) | 2.9 | 25 % |
| stage 3 Newton (2 iterations of 12 product runs + baselines/trials) | 4.3 | 37 % |
| Arnoldi stability estimate (24 runs) | ~3.4 | 29 % |
| setup (netlist, KLU symbolic) | ~0.2 | 2 % |

So **the warm-up is not the cost** (a warm-up period is 0.0036 s, a period run 0.165 s: 45x): the cost is period runs, and 35 % of the runs go
into two Newton attempts that cannot succeed. Measured options (490- and 1150-state R3 PEX, 99/101-state R4/R3d/R4d PEX, `krylov=1`):

| option | result |
|---|---|
| (a1) shorter warm-up (`PSS_SHOOT_WARMUP_PERIODS` 1 / 5; deliberately bad start) | still converges (R3 490: 129 / 154 runs vs 77; R4d fell back to dense at 1 period, 458 runs; R4d at 5: 12 iterations): a Newton iteration costs ~13 runs, a warm-up period 1/45 of a run, so a shorter warm-up loses |
| (a2) longer first stage (50 .. 150 periods as the first stage) | 20-45 % fewer runs, **but the orbit moves by 2e-4 against the dense orbit** (s88 and the dense engine agree to 5e-7 at 1150 states): the period map depends on the recorded step grid, the dense and old-krylov x0 share the 20/40/80 warm-up history, any other history lands at the grid-noise level (1e-5 .. 3e-4, amplified 1/(1-rho) = 50). Rejected: misses the 1e-5 target |
| (a3) **adopted**: same 20/40/80 warm-up transients and endpoints, **no Newton attempt while cumulative warm-up < 130 periods** (`PSS_KRYLOV_WARMUP`, `0` = attempt at every stage; from `PSS_KRYLOV_WARMUP_MIN` = 64 states; an explicit `PSS_SHOOT_WARMUP_PERIODS` or the final stage always attempts) | 1150 states: 75 -> 47 runs, orbit **bit-identical** to s88's; R3 PEX 99-205 states -30..-53 % runs; R3d -16 %; R4 -5 % (its probes fail at 140, 300 and 620 periods too) |
| (b) cheaper warm-up | nothing to win (4 %) |
| (c) GMRES recycling (GCRO-style deflation by stored `(v, (M-I)v)` pairs, within a Newton run and across stages; code kept at `538/pssshoot_with_recycle.c`, not in the tree) | stale images, no check: -12..-50 % runs, but 2 of 16 decks (R4 PEX, R4d PEX; slow mode rho 0.994) stagnated at residual 1.9e-5 and fell back to dense (a 1e-3 error of the stored images is a 17 % error of the slow eigenvalue of M-I); with one true product after the deflation (true residual check) the guess is off by 3e-3 .. 0.2 (FD noise floor of the period map) and the run needs as many products as before: **no gain, not adopted** |
| Eisenstat-Walker forcing term (eta = 0.1 .. 1e-3 by residual) | worse on every deck, R4 fell back to dense: **not adopted** |
| Arnoldi checkpoint every 2 dimensions from 8 (was 4; settle test, 2e-3, unchanged) | stops at 10 instead of 12 on every deck: -4 runs. Adopted |

**Proof** (s88 `krylov=1` vs this build default, 4 alternations, rotated order; wall median, min-max; loadavg of the runs):

| states | s88 runs | new runs | s88 wall s | new wall s | orbit dev vs dense | new vs s88 orbit | Floquet radius (Arnoldi) |
|---|---|---|---|---|---|---|---|
| 490 | 77 | 47 | 6.16 (6.09-6.21) | 4.23 (4.16-4.30) | 2.2e-6 | 0 (identical files) | 0.980301 stable[arnoldi] |
| 1150 | 75 | 47 | 12.60 (11.76-12.74) | 8.53 (8.18-8.83) | 5.3e-7 | 0 | 0.980318 stable[arnoldi] |
| 2470 | 86 | 58 | 32.47 (30.51-32.76) | 23.77 (22.83-23.96) | (no dense run; = s88) | 0 | 0.980297 stable[arnoldi] |

loadavg 19-23 throughout. Orbit deviation = max |dev| / max(peak-to-peak, 1e-2) over v(vout), v(sw), v(opout), i(L1), against the dense orbit
(`krylov=0`, 1028 / 2348 runs, 49.6 s / 248.9 s at this load); dense Floquet radius 0.9803. Other classes (PEX R4/R3d/R4d/R6s 99-101 states):
runs 137 -> 126, 162 -> 136, 150 -> 148, 159 -> 129; orbit identical to s88's, 1e-5 .. 3e-5 of the dense orbit (the dense engine's own variant spread there is 1e-4 .. 2e-2).

### 3. Free-running (`trigger=`) matrix-free: done, no bordered system needed

The free-running engine shoots on the Poincare section: `P(x)` = the state at the next trigger crossing (`pss_run_until_trigger`), whose return time *is* the
period. The dense engine's `M` is the terminal-projected monodromy of that map, so `(M - I) dx = -(P(x0) - x0)` already contains the phase condition and `T`
needs no extra unknown. The matrix-free solve is the card #537 GMRES with `P` = the trigger map (one product = a centred pair of trigger runs),
the Arnoldi estimate for the stability verdict. Same interface as section 1 (`krylov=1`, `NGSPICE_PSS_KRYLOV`, automatic from 64 states, notice line,
`krylov=0` forces dense). On a failed product or GMRES stagnation it continues on the dense monodromy and says so. Output additions:
`PSS-FREERUN: iter=N krylov: GMRES dim D relres R, K trigger runs in the solve`, the CONVERGED line with `stable[arnoldi]`.

**FD step.** The dense free-running engine steps each state by `1e-6 |x|` (floor 1e-9): at the noise floor of the trigger map (card #503). Used for the
matrix-free products it gives a wrong Arnoldi radius on R5 (4.79, "UNSTABLE" for a stable orbit); `PSS_FREERUN_KRY_HMUL` multiplies the step, default 1000
(1e-3 relative; floor 1e-6): radius at x1 / x10 / x100 / x1000 / x10000 = 4.79 / 0.0728 / 0.136 / 0.0301 / 0.0446 (dense 0.0351). The stable/unstable class is
right from x10 up, the value is +-0.05 on R5 (its dense value is itself noisy: the dense run on the 24-state PEX R5 reports 24.2 "UNSTABLE", the matrix-free run 0.0346).

| deck | dense | matrix-free |
|---|---|---|
| R5 hysteretic buck fixture (4 states) | 31 Newton iterations (the cold start wanders), 12.4 s, T 2.36584963 us, rho 0.0351 | 18 iterations, 6.2 s, T 2.36584692 us (1.1e-6), rho 0.0301; orbit dev 9e-6 (smooth vectors; edge vectors 5e-2 of p-p as in the KLU comparison) |
| R5 + 2 RC trees, 10 states | 15 it, 10.9 s, T 2.36594009 us, rho 0.0215 | 4 it, 1.1 s, T 2.36581657 us (5.2e-5; steady_coeff is 1e-4), rho 0.0113 |
| 24 states | 17 it, 40.4 s, rho 24.2 (wrong) | 4 it, 1.6 s, T 2.36601796 us (vs dense 6.8e-6), rho 0.0346 |
| 48 states | **FAILED** to converge in 45 iterations, 636 s | 4 it, 2.9 s, rho 0.0270 |
| 136 states | not run (dense hopeless) | 7 it, 7.2 s, rho 0.0329 |
| 96 states, automatic (test) | | 4.6 s, rho 0.0361, notice present |

### Tests (`tests/test_pss_krylov.py`, edited in the working tree, not committed)

New: `auto` (73-state PEX deck without any card runs matrix-free with exactly one notice; `krylov=0` card, `NGSPICE_PSS_KRYLOV=0`, `PSS_KRYLOV_AUTO_MIN=0` and `=74` -> dense;
`=73` -> krylov; 7-state fixture dense with no notice; automatic orbit equals dense, 1.2e-6), `warmup` (R3 145 states / R4 99 / R3d 101: no-attempt schedule equals attempts-at-every-stage
(`PSS_KRYLOV_WARMUP=0`) and the dense orbit; a deliberately bad start of 1 warm-up period lands within 5e-4), `freerun` (replaces "freerun refused cleanly": R5 and PEX R5 matrix-free vs dense,
overrides, automatic at 96 states), ASan of the free-running krylov run. Red-proofs (must fail, and do): automatic choice disabled by `PSS_KRYLOV_AUTO_MIN=1000`; krylov deck at a 1 % different
clock guess; free-running period tolerance 1e-9; plus the card #537 red-proofs. `pex_deck()` now writes `krylov=0` into its dense arm (it would otherwise go automatic).
Result with the final build: 26 pass + 5 red-proofs failing as they must, 0 could-not-run; test_pss_klu 17/0, test_pss 43/47 (4 info), test_pac 32/32, test_r3d_r4d 21/25, test_tran_from_pss 6/6 -- the s88 numbers.

### Not done / limits (2026-10-02)

* Threshold stays 64 although wall time crosses at ~40-46 states on R3-class decks (kept conservative; `PSS_KRYLOV_AUTO_MIN`).
* Skipping Newton attempts during the first 130 periods helps R4-class decks only 2-17 %: their probes also fail at 140/300/620 periods. A deck-dependent rule (from the stabtime, or from a cheap drift measure) was not found; 130 is tuned on R3/R3d/R4/R4d PEX decks.
* The Arnoldi estimate is now ~40 % of the remaining runs (20 of 47); settle test compares dimension m with m-4 (the radius is already right at 8); a cheaper rule or reuse of the last GMRES basis was not tried.
* Free-running: the radius value is noisy (+-0.05 on R5; verdict class right), `PSS_FREERUN_KRY_HMUL` (1000) tuned on R5 and PEX R5 only, no dense reference above 48 states (dense fails there), `.pac/.pstb` on a free-running orbit stays refused (and disables the automatic choice); the dense run with the same fixture needs 31 iterations because of the cold start, matrix-free needed 18: warm-up of the free-running engine (one StabTime transient) untouched.
* Orbit of the dense engine is the reference only to the grid-noise level; the 1e-5 target holds for the s88 warm-up history, not for any other.
* Not tested: BSIM4-heavy decks, XSPICE-clocked (`driven=1`) decks, `.pac` with `pacpoints` after an automatic PSS, more than 2470 states, dense 2470.

## Free-running Floquet verdict, dense Newton step, matrix-free by default (card #540, 2026-10-02)

Defect: the dense free-running shooter (`.pss ... trigger=`, below 64 states or `krylov=0`) printed radius 24.2 `UNSTABLE` on the stable 24-state
PEX-loaded R5, did not converge on the 48-state one (636 s, 45 iterations), and on the `pac_r5_line2out` fixture printed radius 0.888 `stable` where the
truth is 0.037 (right class, wrong value, and 0.888 is close enough to 1 to mislead). The matrix-free path read 4.79 on R5 until its step was scaled (card #538).
A simulator must never print a wrong stability verdict.

### Root cause

The monodromy is differenced on the Poincare return map `P(x)` = state at the next trigger crossing. `P` is an adaptive-step integration at the engine's
`reltol` 1e-5 (the cap at the top of `PSSshootFreerun`): two runs from states 1e-6 apart take different step sequences and agree only to ~1e-5 of the state.
Step 1e-6 relative therefore differences noise: `noise/h` ~ 10 per matrix entry, and the Floquet radius of that matrix is whatever the noise makes it
(0.02 ... 24 on the same orbit). It also made the Newton Jacobian poor (31 iterations on R5, 17 on the 24-state PEX R5, no convergence at 48 states).
Evidence: with `.options reltol=1e-9` on the deck the dense radius is 0.0368 at steps 1e-6, 1e-5 and 1e-4 (all three equal), 4 Newton iterations instead of 15
on the 10-state deck; at the default reltol the radius vs step is 0.035 / 0.017 / 0.013 / 0.154 / 0.011 at 1e-6 / 1e-5 / 1e-4 / 1e-3 / 1e-2 (R5).

### Independent truth (no shooting, no Jacobian)

The settled orbit's state is perturbed by 0.3 % x 0.5..1.5 (random signs) in a plain `.tran` (reltol 1e-9, gear), the Poincare section is `v(vout)` rising
through its orbit mean, and the decay per period of the deviation (period, `i(L1)`, `v(cx)` at the section) from the final cycles is measured
(`~/Work/data/ngspice_augmentation/540/truth.py`, and `tran_multiplier()` in `tests/test_pss_krylov.py`). Ratios of successive deviations settle at
0.036 - 0.043 (cycles 3-5) on R5 and on the PEX R5 at 10 / 24 / 48 states, two seeds each: slowest multiplier **0.037** (the RC trees are passive, the slow mode is R5's own).
The same `.tran` gives the period 2.36710 us; both engines print 2.3658 - 2.3661 us (5e-4 low, the reltol 1e-5 floor of card #503; the polish at card #503/#517 only runs for `.pac`).
This card did not change the period or the orbit.

### What changed

1. **Stability pass (both engines)**, after convergence, under a scoped tight profile (`reltol` 1e-8, `PSS_FREERUN_STAB_RELTOL`, only ever tightens, restored before anything
   else runs, so the `.pac`/`.pstb` hand-off is unchanged): one tight base trigger run, then one-sided differences. Dense: `n` columns at step 1e-5 |x|
   (`PSS_FREERUN_STAB_EPS`), eigenvalues, trivial (smallest) multiplier excluded. Matrix-free: the Arnoldi process of card #538, with per-state step ceilings
   proportional to the state's ripple (1e-3 of it, `PSS_FREERUN_STAB_RIPPLE`; |x| 1e-6 and 1e-9 as floors). With ceilings proportional to |x| the stiffest state (Csw, ripple ~ |x|)
   fixes the product size and vout is moved by 1e-8 of itself, 10 % of the tight noise floor (first try: 0.12 instead of 0.037 at 48 states, consistently, so the guard below did not see it).
2. **Guard.** The pass is made twice, at h and h/2. Noise-dominated differences scale like 1/h, truncation-dominated like h: they disagree exactly when the
   estimate is not trustworthy. The verdict is trusted only if both radii are on the same side of 1 and agree to 10 % (0.02 absolute below 0.2); otherwise the
   CONVERGED line reads `spectral_radius=-1.000000 UNKNOWN(not-computed)` preceded by `PSS-FREERUN: stability not computed: <reason>` (also when the Arnoldi estimate is
   not settled, a base/column run fails, or the eigen-solve fails). The reported radius is the larger of the two. Same token as card #537.
   Forcing the old noise (`PSS_FREERUN_STAB_RELTOL=1e-5 PSS_FREERUN_STAB_EPS=1e-6 PSS_FREERUN_STAB_RIPPLE=1e-5`) gives `UNKNOWN(not-computed)` on both engines
   (radii 22.7 vs 0.07 dense, 1.3 vs 2.7 matrix-free on the 10-state deck; 7.0 vs 13.9 on 24 states), never `UNSTABLE`.
3. **Dense Newton step** 1e-4 relative (floor 1e-7; was 1e-6, floor 1e-9), 100x above the return-map noise. `PSS_FREERUN_FD_EPS` still overrides (and applies to the matrix-free
   products as eps x `PSS_FREERUN_KRY_HMUL`, unchanged 1e-3). With `.pac`/`.pstb` on the orbit and the dense engine (explicit `krylov=0`) the step stays 1e-6: the cold start of the
   pac_r5 fixtures (residual 0.7) takes 24-41 iterations at 1e-6 .. 1e-5 and wanders into a lost crossing at 3e-5 and 1e-4.
4. **Free-running is matrix-free from the first state** (`PSS_FREERUN_AUTO_MIN`, default 1; `krylov=0` or `NGSPICE_PSS_KRYLOV=0` = dense; `PSS_KRYLOV_AUTO_MIN` stays the driven engine's 64),
   **also with `.pac`/`.pstb`** on the orbit (was dense only).

### Choice: matrix-free by default, dense kept as the explicit fallback

(development workstation, load 5-18, OMP 2, `.options klu`, builds `aaba50ba3f49-540-N`; wall in s incl. the 450 us warm-up; one run each, matrix-free = automatic)

| deck | independent multiplier | dense before (s90) | dense after | matrix-free before (#538) | matrix-free after | period dense / m-f |
|---|---|---|---|---|---|---|
| R5 fixture, 4 states | 0.037-0.038 | 0.0351 stable, 31 it, 7.8 s | 0.0384 stable, 15 it, 4.7 s | 0.0301 (4.79 unscaled), 18 it, 4.5 s | 0.0368 stable[arnoldi], 18 it, 4.9 s | 2.36584256 / 2.36584692 us |
| PEX R5 10 states | 0.037 | 0.0215, 15 it, 10.9 s | 0.0384, 4 it, 3.8 s | 0.0113, 4 it, 1.1 s | 0.0367, 4 it, 2.1 s | 2.36593766 / 2.36581657 us |
| PEX R5 24 states | 0.037-0.041 | **24.24 UNSTABLE**, 17 it, 40 s | 0.0384 stable, 8 it, 17 s | 0.0346, 4 it, 1.6 s | 0.0368 stable[arnoldi], 4 it, 2.8 s | 2.36601967 / 2.36601796 us |
| PEX R5 48 states | 0.037 | **FAILED**, 45 it, 636 s | 0.0384 stable, 9 it, 57 s | 0.0270, 4 it, 2.9 s | 0.0369 stable[arnoldi], 4 it, 4.7 s | 2.36610072 / 2.36609845 us |
| PEX R5 136 states | not run | - | not run | 0.0329, 7 it, 7.2 s | 0.0369 stable[arnoldi], 7 it, 9.9 s | - / 2.36605571 us |
| pac_r5_line2out (`.pac`) | - | **0.888 stable**, 29 it, 12.1 s | | - | 0.0368, 20 it, 9.7 s (`.pac` ran) | 2.36585312 / 2.36583518 us |
| pac_r5_ctrl2out 20k (`.pac`) | - | 0.0507, 41 it, 15.6 s | | - | 0.0368, 17 it, 8.7 s (`.pac` ran) | 2.36585204 / 2.36584808 us |

Matrix-free is equal at 4 states (wall, iterations) and better everywhere above (period within 2e-6 of dense up to 24 states, radius within 0.5 % of the independent value, 3x - 12x
faster from 10 to 48 states), and it converges the pac_r5 cold starts in 17-20 iterations instead of 29-41. No deck where dense wins by more than the run-to-run spread; so by the
evidence matrix-free is the default at all sizes. The cost of the corrected stability pass is 17 tight trigger runs (+1 - 3 s at 24 - 48 states, tight runs are ~10x a normal one),
the old noisy estimate was cheaper and wrong. The dense radius is 0.0384 on every deck (one-sided at 1e-5, not the central 0.0368 of the tight-profile run: 4 %, the
truncation of a one-sided quotient on a switching map; the two agree to 10 % by the guard's own rule).

### `.pac` / `.pstb` after a matrix-free free-running PSS: identical, enabled

`pac_r5_ctrl2out_{20k,60k,100k}` (the tests/test_pac.py R5 checks): dense-path numbers (s90) vs after a matrix-free PSS (this build), `|H|` and phase:
20k 1.00474447 / 7.32212e-4 vs 1.00474439 / 7.32104e-4; 60k 1.03900182 / -6.36473e-3 vs 1.03900170 / -6.36519e-3; 100k 1.09744530 / -3.64659e-2 vs 1.09744513 / -3.64667e-2
(relative differences <= 3.4e-7 in magnitude, <= 8e-7 rad). Against `pac_perturbed --freerun`: 0.0501 / 0.0535 / 0.0435 dB and 0.1317 / 0.0990 / 0.2268 deg, same as before (worst 0.0535 dB / 0.2268 deg;
the 0.435 dB / 1.38 deg quoted in the card are the older numbers, the s90 log has 0.0535 / 0.2269). The 150 kHz refusal and the experimental override are untouched
(tests/test_pac.py). So `.pac`/`.pstb` is no longer a reason for the dense engine.

### Tests

`tests/test_pss_krylov.py` (working tree, not committed): `fr540` group, `--only fr540`:
`freerun stability vs independent .tran multiplier` (R5 4 states and PEX R5 24 states, both engines must be `stable` within 25 % of a multiplier measured by perturbing the orbit in a `.tran`
inside the test; the s90 binary prints 24.24 UNSTABLE on the 24-state deck and fails it), `freerun stability guard` (forced noise gives `UNKNOWN(not-computed)`, radius -1, a reason line, on both engines;
control stable), `freerun default routing, dense Newton step` (automatic matrix-free at 4 states with notice, `NGSPICE_PSS_KRYLOV=0` dense with the eigenvalue list, dense Newton on the 24-state deck <= 12 iterations
and the matrix-free period, `.pac` runs after the automatic PSS with the dense-path value). Red-proofs (must fail, and do): stability tolerance 1e-4; demanding a verdict from the forced-noise run.
`tests/test_pss.py`: R5 free-running now also runs `NGSPICE_PSS_KRYLOV=0` and checks the dense radius (0.02 - 0.06) and the trivial-multiplier line (it was silently skipped once the default became matrix-free).

### Not done / limits (2026-10-02)

* The period/orbit of both engines is 5e-4 below the tight-profile value (2.36584 vs 2.36710 us); only the `.pac` hand-off polishes. Not changed here.
* The independent multiplier is a median of period-deviation ratios over cycles 2-5 of a heavily damped orbit (0.037, +-10 % between seeds); it resolves "0.04 not 24", not 0.037 vs 0.039.
* Only R5 (a single heavily damped orbit, radius 0.037) was measured; the guard's bounds (10 % / 0.02, same side of 1) were not exercised on a free-running orbit with a radius near 1
  (no such fixture); in that region a disagreement straddling 1 gives UNKNOWN by construction. Cost of the pass on very large free-running decks (> 1000 states) not measured.
* Dense above 48 states, and the dense stability pass for `n` in the hundreds (2n tight runs), were not run; matrix-free is the default there.
* BSIM4-heavy free-running decks, XSPICE-clocked free-running decks: none exist here, untested.

## Generic charge history for `solver=td` (card #542, 2026-10-02)

Fixes the silent error of card #539 "bug 1": `.pac/.pstb/.pxf/.pnoise ... solver=td` carried the history of native
capacitors and inductors only, while the recorded Jacobian `J_k = G_k + ag0*C_k` also contains the model charges
(MOS1 overlap/Meyer/junction, BSIM4 intrinsic/overlap/junction, diode cj/tt, BJT, mutual inductance). Every model
capacitance then acted as a conductance `ag0*C` to a perturbation-free reference (an inverter chain lost 0.55, 0.31, 0.04,
0.001 of the timing per stage, worse on a finer grid; the R7 BSIM4 buck had a loop gain of 2e-10).

### Recursion (every charge state, native and model)

Linearised trapezoidal / BE step, `Q` = perturbation of the charge entering each equation row, `T` = its time derivative:

    (G_k + ag0 C_k) y_k = b_k + ag0 Q_{k-1} + a1 T_{k-1}
    Q_k = C_k y_k            T_k = ag0 (Q_k - Q_{k-1}) - a1 T_{k-1}        (a1 = ag1 for order 2, 0 for order 1)

`C_k = dq/dx` is a full sparse matrix and is obtained with **no per-device code**: every charge-conserving device
stamps its capacitive part through `NIintegrate` as `geq = ag0*cap`, so the loaded matrix is affine in `ag0` and
`C_k = (J(ag0) - J(ag0/2)) / (ag0/2)` entry by entry (`pac_gen_loads()`, two loads at the true accepted point `x_k`,
state/rhs/mode/bypass saved and restored; new `SMPnumElements/SMPcopyElements/SMPsetElements` read and rewrite the loaded
matrix under sparse1.3 and KLU). The native C/L recursion is the special case `C` = diag(C, -L) and is kept bit-identical
for circuits without model charges (it is also what `PAC_FORCE_GENERIC=1` is cross-checked against: 1e-10).

### Which path a circuit takes

1. Only Resistor, Capacitor, Inductor, V/I sources, controlled sources, Switch, CSwitch, ASRC: native path, nothing new runs
   (zero cost, results bit-identical to s90; verified on test_pac/test_pxf/test_pnoise).
2. Another device type present: one **detection period** loads each such type alone at ag0 and ag0/2; no difference
   anywhere on the orbit (e.g. a diode with cj = tt = 0 beside an ideal switch) -> native path (saltation still works). The
   detection period snapshots and restores CKTrhs/CKTrhsOld/CKTstates and the integration scalars, because the legacy flow's
   first recorded period starts from whatever those hold: without the restore R4d's `.pstb` moved by 1.8e-3 (it contains chargeless
   diodes), with it the output is bit-identical to s90.
3. Mutual inductor, or detection finds a charge: **generic path**.
4. `Tranline/TransLine/CplLines/LTRA` (delay history, not a charge): refused.

### State set, prepass, first step

The Floquet state (Phi, xi, z) of the generic path is the set of matrix rows that are a **column of C_k anywhere on the
orbit** (node voltages and inductor/flux branch currents; 3-4 rows on the LTI decks, 44-45 on R7 against 14 native). It is found by
a prepass period (`pac_generic_prepare()`), which also runs the **affinity guard**: a third load at 2*ag0 must satisfy
`J(2a) = 3 J(a) - 2 J(a/2)` entrywise to 1e-7; otherwise the analysis is REFUSED naming the device type (each type is loaded
alone to find it). A device whose load is not affine in ag0 cannot silently produce a curve (none was found among MOS1/MOS3/BSIM3v3.3/BSIM4
incl. capmod 0/2/3 and rgatemod=1, JFET, diode, BJT, mutual, with `.options bypass=1` too: the loads are forced to bypass=0).
Non-quasi-static BSIM3/BSIM4 (`nqsmod`, `trnqsmod`, `acnqsmod` != 0) pass the affinity test but carry an internal relaxation
state (qcheq) that is a CKTstate history of its own, not a node charge: REFUSED (`PAC_ALLOW_NQS=1` overrides; measured 3e-5 off .ac on a
ps-class NQS deck). The first accepted step seeds `Q_0 = C_1 y_0` on the state rows and `T_0 = 0`
(the recorded period restarts with a first order-1 step, as before).

### Start of the recorded period (restart artefact, measured and mostly removed)

The recording runs start with UIC from the PSS orbit's native states only. The history of every model charge comes from
the INITJCT load, i.e. from the device-level ICs (MOS icVDS/icVGS/icVBS, diode icVD, BJT icVBE/icVCE) that `CKTic()` derived
once from zero node ICs; so nodes carried only by model charges start every period at 0 V and relax over the first
ns..us, a transient the orbit does not have (refs/README.md R7: "first ~15 ns of the orbit wrong"). The linearisation along it is off
by about (startup duration)/T: MOS1 amplifier DC gain 2.7e-4 (T = 1 ms), diode 1.9e-5. The generic path therefore (a) runs
two prepass periods and keeps the end-of-period node vector, (b) restores it into `ckt->CKTrhs` before every recording run
and re-derives the device ICs with `DEVsetic` (native C/L excluded: `pac_state_set_ic` sets those from x0; ICs the deck gave explicitly
are honoured). Result: 2.7e-4 -> 3.6e-9 (MOS1), BSIM4 3e-7 -> 2e-10, diode 1.9e-5 -> 1.2e-7. Residual: a charge whose
history has no device-level ic (BJT collector-substrate cjs: 7e-6, bounded 2e-5 in the test). `PAC_NO_IC_INJECT=1` switches (b) off.

### Adjoint (`.pxf adjoint=1`) and `.pnoise`

The backward sweep transposes exactly the forward map. With `nu_k = J_k^{-T} z_k`, `z_k = w_k + C_k^T (alpha_k + ag0_k beta_k)`,
`gamma_k = nu_k - beta_k`, `alpha_{k-1} = ag0_k gamma_k`, `beta_{k-1} = a1_k gamma_k` (alpha, beta: adjoints of Q, T; the
periodic closure M(w) enters `w_N` on the state rows). `J_k` and `C_k` are **stored** during the forward trajectory pass (16 B per
stored matrix entry, cap `PAC_TRAJ_MAX_MB`, default 3000, beyond it `.pxf` falls back to the per-source forward engine and
`.pnoise` reports it) instead of being re-loaded from the devices in reverse order: a reload is not exact for devices whose load
reads `CKTstate1` (MOS1 Meyer averages the previous step's capacitance; first measured symptom: adjoint 1e57 on a MOS1 chain)
or limiting state. With stored matrices adjoint == forward to 1e-12 on LTI decks and 1e-14 on R7 (3 sources). `DEVnoise` still needs
the device instance values at `x_k` (one limiting-warmed load per step, then the stored `J_k` replaces the matrix).

### What happens at an event, free-running, delay lines

- Saltation (ideal `sw/csw` event) with model charges on the orbit: REFUSED (Fang's jump is derived for native C/L states;
  dq continuity through a switched topology is not derived). Native circuits and chargeless diodes beside switches are unchanged.
- `trigger=` free-running orbit with model charges: REFUSED (period-variation projection needs f(T) at every charge state).
- KLU: same results (9e-15 on a MOS1 deck, 2e-9 on the BSIM4 chain at pacpoints 4000). On R7 at Points 400 KLU and sparse differ by
  1.4e-3 (imag part 0.185 vs 0.183), the same spread as a 1e-4 relative change of reltol (1.5016..1.5018 / -0.1851..-0.1807): the base
  trajectory's adaptive step sequence at the comparator/driver edges sets a ~1e-3 (0.2 deg, 0.02 dB) repeatability floor of R7's td PAC at reltol 1e-4.

### Numbers (aaba50ba3f49-542-3, development workstation, OMP 2)

LTI decks vs plain `.ac` (slow clock on an unrelated net, bias fixed, pacpoints 400000, 1 kHz..100 kHz, worst point), `tests/test_pac_charge.py`:

| deck | worst rel | | deck | worst rel |
|---|---|---|---|---|
| MOS1 gate-cap ladder | 2.2e-7 | | BSIM4 CS capmod 0 | 1.9e-10 |
| MOS1 amp (Meyer+overlap+junction) | 3.6e-9 | | BSIM4 CS capmod 2 | 2.0e-10 |
| diode (cj, tt) | 1.2e-7 | | BSIM4 CS rgatemod 1 | 2.0e-10 |
| BJT amp (cje cjc tf tr) | 2.0e-8 | | BSIM3v3.3 amp | 1.9e-10 |
| BJT + cjs (substrate) | 7.3e-6 (bounded 2e-5) | | MOS3 amp / JFET amp | 3.6e-9 / 6.9e-8 |
| transformer (K) | 3.4e-7 | | BSIM4 capmod 3 (pacpoints 1e5) | 3.2e-9 |

bug 1 chain (control to average output, `.pac` pacpoints 40000 vs static `.tran` derivative, reltol 1e-5): MOS1+cgso/cgdo 1e-8 N=1,2,3:
0.5355/0.4984/0.5045 vs 0.5355/0.4985/0.5045; BSIM4 N=1,2,3: 0.5648/0.5020/0.5151 vs 0.5650/0.5020/0.5155 (worst 0.07 %; before: 0.2507/0.0600/0.0138 and
0.0214/0.0201/0.0028). The static reference needs reltol 1e-5 (default reltol gives 0.5249/0.5575 for the same chain depending on the step).

R7 (BSIM4 buck, converged FRA `fra/curves-r7-20261001`, reltol 1e-4, Points 400, 8 frequencies 500 Hz..200 kHz), worst point:
ctrl-to-out 0.105 dB / 0.93 deg, line-to-out 0.173 dB / 0.79 deg, Zout 0.243 dB / 0.62 deg, loop gain 0.149 dB / 0.31 deg (largest at 500 Hz, where
the FRA uncertainty is 0.09-0.12 dB; at 40-200 kHz ctrl-to-out and loop gain sit 0.10-0.12 dB above the FRA, 0.07 dB with reltol 2e-5, 0.015 dB less with
pacpoints 4000). `.pxf` adjoint == forward to <= 1e-12 at 500 Hz / 20 kHz / 200 kHz for Vin, Vctrl, Iinj.

Cost (R7 `.pac`, 1 point): PAC-TD engine 0.49 s generic vs 0.34 s with `PAC_NATIVE_HISTORY=1`, plus the prepass (2 periods + detection) 0.68 s; the whole job
is dominated by `.pss` (22 s): +1.2 s (5 %). `pac_r3_ctrl2out` (native devices only): no extra work, output bit-identical to s90, PAC-TD 0.006 s both.

### Switches (env, test only)

`PAC_ALLOW_NQS=1` accept NQS models (approximate),
`PAC_NATIVE_HISTORY=1` old native-only history (red proof: reproduces bug 1, MOS1 amp error 0.81, chain 88 %),
`PAC_FORCE_GENERIC=1` generic path on a native deck, `PAC_NO_IC_INJECT=1`, `PAC_GEN_TRACE=1`, `PAC_TRAJ_MAX_MB`,
`PAC_GEN_TEST_NONAFFINE=<type>` (makes a type look non-affine: guard test).

### Limits and what is not claimed

- MOS1 Meyer and BSIM4 capmod 0 charges are not exact functions of x (the model integrates increments): the generic
  recursion uses the stamped C (what Newton uses), exact for constant bias (LTI proofs) and for slowly varying C; the extra
  `v dC/dx` term of a nonconservative charge is not in the linearisation.
- The trapezoidal ring mode (`T_k = ... - T_{k-1}`) is neutral and not part of the Floquet closure (restart per period with T_0 = 0, as
  for native C/L); error O(w h) of one step, invisible on the LTI decks (1e-10..1e-7).
- `.pnoise` on R7: Rload/Resr thermal values now include the closed loop and the sideband folding (+9.6 % over the old passive-LTI value);
  not validated against an independent reference (no MC for BSIM4); MOS channel noise not validated.
- Memory: the generic engine keeps `(neq+1) x M` complex Q, T, y arrays plus the stored trajectory matrices; R7p (413 native states,
  dense Phi of (state rows)^2) was not run on this path.

## Model-charge td PAC, continued: `.pnoise` level, R7p, free-running orbits, saltation (card #544, 2026-10-02)

Branch aug-544 on aug-541 (s96), three commits: `e2ef937f7` pnoise (item 1), `240254301` saltation (item 4), `6e88e4d8e` free-running (item 3); final build
`aaba50ba3f49-544-13` (= HEAD), ASan `aaba50ba3f49-544-asan`. Four jobs, in the order of the card; each number below is from a run in
`~/Work/data/ngspice_augmentation/544/` and each check can fail (red proofs listed with it). Item 2 needed no code (a trajectory report line only).

### 1. `.pnoise` absolute level with device-model noise

**LTI-biased decks vs plain `.noise`** (`tests/test_pac_charge.py check_pnoise_dev`: slow 1 ms clock on an unrelated net, bias fixed, `dev=1`, `.noise ... 1`
for the per-device vectors, 5 points 1k..100k, pacpoints 200000). Total and every device contribution (thermal, shot, flicker, `rd/rs/rg`, ...):

| deck | worst rel (total and per device) |
|---|---|
| MOS1 CS amp (kf, Meyer+overlap+junction) | 4.3e-8 |
| BSIM4 CS amp (capmod 2, flicker 99.9 % of the power at 1 kHz, thermal above) | 3.8e-8 |
| diode (cj, tt, rs, kf) | 6.5e-7 |
| BJT amp (cje cjc tf tr rb re kf) | 8.8e-7 |
| JFET / MOS3 / BSIM3v3.3 amp | 4.5e-8 / 4.3e-8 / 3.7e-8 |

The engine's whitelist (`pac_noise_check_devices`) had only resistor, diode, MOS1, BSIM4: **BJT, JFET, MOS3, BSIM3 are now admitted** (they use the same
`NevalSrc` gain idiom; the table is their validation). Other noisy types (Switch, ...) are still refused.

*Flicker under a periodic bias.* The engine evaluates `S_g(t_k, f)` at the OUTPUT frequency `f` (docs "What was skipped" of card #509): for a
periodic bias it treats flicker as a modulated source that is white across the sidebands, i.e. a source at `f + k fs` gets the density it has at `f`,
not at `f + k fs`. That is exact for LTI bias (table above) and for a slowly varying bias (item below), and it OVERESTIMATES the folded flicker of a
fast-switching circuit (the true density at `f + k fs` is lower). Not validated for a switched circuit (no reference; the FD engine's sideband
folding is not one either, see below); thermal and shot noise do not have the issue.

**A real bug, found by the time-varying check, fixed.** `DIOnoise`, `BJTnoise`, `MOS3noise` and `BSIM4noise` read the operating point from `CKTstate0`
(diode and BJT currents, `vds/vgs`). The backward adjoint sweep never advances device state, so at every step k they saw the state of the END of the
period (the generic path even restored `CKTstate0` after its warm loads). Constant bias hides it, so every LTI check passed. With a time-varying bias a
switching diode's shot/flicker density, or a BSIM4's `vds` term, was wrong. Fix (`pac_noise_warm_load`): the noise pass of the forward recording keeps the
converged `CKTstate0` of every step (`PACtrajStep.st0`, `nstates x 8 B` per step, counted in `PAC_TRAJ_MAX_MB`); the backward sweep restores it and does
one load at `x_k` with the integration coefficients ZERO, so the capacitive current the transient load adds to the junction-current state (computed
from stale `CKTstate1`: 10^6-fold currents on the first steps, ~1 % elsewhere) vanishes and the density is the DC junction current plain `.noise` uses.
The native-path adjoint did not reload the devices at all (stored `J_k` only), same fix.

**Time-varying bias vs a quasi-static brute force** (`check_pnoise_qs`): a 100 kHz LO sine on the gate of a common-source stage whose time constant
(1-3 ns) is 3-4 orders below the LO period and whose device stays in saturation while it conducts (smooth `gm(t)`). Then the PSD at 1 kHz is the
time average of the instantaneous PSD, which plain `.noise` gives at the DC bias of N = 96 equally spaced LO phases (error O(tau/T) ~ 1e-4; N = 96,
192, 384 identical to 6 digits):

| deck | engine - quasi-static (total and per device) | before the fix |
|---|---|---|
| MOS1 (Rd 2k, channel thermal 3.7e-9 V/rtHz, Rd 5.7e-9) | 3.0e-5 | 3.0e-5 |
| BSIM4 (capmod 2, flicker dominated) | 1.5e-7 | 5.1e-3 |
| diode (cj, tt, rs, kf) | 1.6e-6 | 4.3 x |
| BJT (cje cjc tf tr rb re kf) | 3.0e-4 | 170 x |

Red proof: `PAC_NOISE_STALE_STATE=1` (skips the warm load: the pre-#544 operating point) fails the diode check. A first version of this check (Rd 100k,
a conduction sliver of 0.3-3 us per period) is NOT a reference: the sliver is narrower than any 96-point phase grid, and the FD `.pnoise` (20
harmonics) smooths it (FD 3.15e-8, TD 4.73e-8, quasi-static 4.46e-8 at 100 kHz LO; the TD value moves from 4.1e-8 to 3.9e-8 as the LO goes
10 -> 100 kHz because tau/T is then not small). Kept in the notes so nobody repeats the FD comparison there.

**Switched MOS sampler vs the kT/C closed form** (`tests/test_pnoise.py check_pnoise_track_hold_mos`; fixtures `pnoise_track_hold_mos1.cir`,
`..._bsim4.cir`): `pnoise_track_hold_slow.cir` with the tanh switch replaced by a gate-clocked MOSFET (MOS1 with Meyer + overlap, BSIM4 capmod 2),
fclk 10 kHz, Rtrack 50k, C 50p (tau 2.5 us, 20 us = 8 tau track window), compared with `th_analytic_S0` (exact folded kT/C of the Ornstein-Uhlenbeck
restart process) with R = Rtrack + Ron (Ron 455 / 120 ohm estimated by hand, moves the closed form by < 0.01 dB):

| switch | 100 Hz | 1 kHz | PAC_NATIVE_HISTORY=1 (red proof) |
|---|---|---|---|
| MOS1 (W/L 10/1) | -0.039 dB | -0.0003 dB | -0.97 dB |
| BSIM4 (8/0.5) | +0.006 dB | +0.045 dB | -3.45 dB |

The tanh-switch fixture is -0.016 dB; the naive kT/C value is 1.38 dB away. The fixture is sized so the textbook formula applies: with W/L = 400/0.5
(channel capacitance 23 % of C1) the same deck is 0.8 dB off, which is the charge-sharing of the channel charge at turn-off, a physical effect not in
the formula and not an engine error.

**R7 (BSIM4 buck): absolute level NOT validated.** A brute-force estimate is not feasible: BSIM4 has no transient noise in ngspice, `tools/pnoise_mc.py`
handles resistors only, and a trnoise Monte Carlo of the resistors alone (20 ms x 32 seeds at 2 ns steps of a 400-device deck) is hours. What is
checked: the two resistor contributions bound from below by the FRA `Zout` (test_r7, +8 %), R7p == R7 for `onoise_rload` at 200 kHz to 5 digits
(1.18508e-11 vs 1.1851e-11), and the generic machinery for MOS channel noise under a time-varying bias and with charges: the quasi-static and kT/C
checks above. The loop-shaped part and the sideband folding of the MOS noise of the converter itself have no independent reference.

### 2. R7p (413 states) on the generic td PAC path

`tests/fixtures/r7p_{pstb,pac_ctrl2out}_{sparse,klu}.cir`, `r7p_pnoise.cir`, `r7p_pxf_adj.cir` (+ `r7_pstb_sparse`, `r7_pac_ctrl2out_sparse`), made by
`~/Work/data/ngspice_augmentation/544/mkfix544.py`; `tests/test_r7.py check_r7p_td`. `.pss` matrix-free (krylov chosen automatically at 413 states),
then the generic path: 413 native C/L + model charges -> **407-408 state rows**, 2192-2199 recorded steps, charge-state prepass 1.7-1.8 s.

| job (development workstation, load 15-23, OMP 2) | wall | of which PSS | PAC engine | peak RSS |
|---|---|---|---|---|
| `.pstb` 2 points, sparse / klu | 80 / 97 s | ~60 s | 17.9 / 21.1 s (4.5-5.3 s per point) | 1.96 / 1.94 GB |
| `.pac` 2 points, klu / sparse | 81 / 80 s | ~60 s | 11.5 / 9.2 s (4.6-5.7 s per point) | 2.0 / 1.56 GB |
| `.pnoise dev=1` 2 points, klu | 94 s | ~60 s | 32.3 s | 2.05 GB |
| `.pxf adjoint=1` 1 point, klu | 92 s | ~60 s | 22.9 s | 2.02 GB |
| R7 `.pstb` / `.pac` 2 points (reference) | 34 / 29 s | | 1.8 s / 0.8 s | 0.22 GB |

*Memory.* `PAC_TRAJ_REPORT=1` (or `PAC_GEN_TRACE=1`) prints what the adjoint trajectory costs: R7p **77.0 MB of `J_k`/`C_k` + 6.9 MB of `x_k`** for 2194 steps (about
2200 stored entries per step, 16 B each). That is 2.8 % of the default `PAC_TRAJ_MAX_MB = 3000`, which fits ~36 R7p; `C_k` already shares `J_k`'s pattern
only by value, so sharing the pattern across steps would halve it (to ~40 MB): not done, it is not wasteful. The forward `.pac/.pstb` engines store
nothing. The ~2 GB peak RSS is the matrix-free `.pss` (454 MB at the start of the Newton phase, 1.85 GB at its end, +40 MB for the PAC), not the PAC
engine; a plain `r7p_pss_shoot_both` (with its wrdata of 7 vectors) peaks at 3.4 GB. The Phi/G/z arrays are `(neq+1) x (n + nfreq)` complex with
neq ~ 500: a few MB. Not investigated: why the krylov PSS grows by 1.4 GB.

*Agreement with R7* (the parasitics are small; expected within ~0.5 dB of R7's FRA; measured):

| point | R7 td | R7p sparse | R7p klu | FRA (R7, unc) |
|---|---|---|---|---|
| loop gain 20 kHz | 4.380 dB | 4.345 dB | 4.365 dB | 4.290 dB (0.009) |
| loop gain 80 kHz | -12.734 dB | -12.758 dB | -12.745 dB | -12.831 dB (0.012) |
| phase margin | 56.09 deg at 28.52 kHz | 56.15 deg at 28.44 kHz | 56.18 deg at 28.49 kHz | 53.5 deg (docs; different definition at the crossing) |
| ctrl-to-out 20 kHz | 3.594 dB / -6.94 deg | 3.580 / -7.32 | 3.580 / -7.18 | 3.627 / -7.79 (0.01 / 0.2) |
| ctrl-to-out 80 kHz | -9.608 / -84.28 | -9.668 / -84.36 | -9.641 / -84.33 | -9.690 / -84.41 |

R7p - R7 <= 0.06 dB / 0.4 deg, KLU - sparse <= 0.03 dB / 0.15 deg (the known repeatability floor of the td PAC at reltol 1e-4), R7p - FRA <= 0.09 dB
on the loop gain and ctrl-to-out (the R7 td engine itself is 0.10 dB above the FRA at 20-200 kHz). `check_r7p_td` bounds: R7p vs R7 0.15 dB / 1 deg, vs
FRA the project target `max(0.5 dB, 3 unc)` / 3 deg, KLU vs sparse 0.1 dB, RSS < 3000 MB.

### 3. Free-running (`trigger=`) orbits with model charges

**Refusal lifted** in `pac_generic_prepare`; `pac_find_trig_idx` accepts a generic state row; the period-variation projection stays what it has been
since card #517: opt-in (`PAC_FREERUN_PROJ`, "not needed for w != k ws") and is **not formed on a generic state set** (f(T) over the charge rows is not
available; the opt-in prints a note and is ignored). So the generic free-running `.pac/.pstb` is the fixed-period engine on the polished orbit, as the
native one is by default. `.pnoise` and `.pxf adjoint=1` on a free-running orbit stay refused (the adjoint is not transposed through the trigger).

**What had to change below it (pssshoot.c): the free-running PSS did not work on a transistor-level circuit at all.** `pss_run_until_trigger` refines
the crossing time by restarting from the bracketing step's native states a few picoseconds before the crossing; every non-native node restarts at zero,
so the trigger node follows the restarted, wrong, derivative ("Timestep too small" at t = 3e-16 s on the first try, or a refinement that never
converges). For circuits with hidden dynamic state (any device other than R, C, L, V, I, controlled sources, switches, mutual, ASRC, or a diode with
cjo/tt/cjsw: `pss_trig_hidden_state`) the crossing, its time and the native states are now **interpolated** linearly from the two accepted steps that
bracket it. Native decks (R5) are bit-identical (`PSS_TRIG_RESTART_REFINE=1` / `PSS_TRIG_INTERP=1` force either). Two defaults for such circuits: the
polish and the stability estimate use reltol 1e-6, not 1e-9 / 1e-8 (1e-9 is below the evaluation noise of a charge-conserving BSIM4: "polish iteration lost
the trigger crossing", then `.pac` skipped), and no 1e-13 s forward-difference probe for f(T) (it restarts too).

**R8**, `refs/r8_bsim_hysteretic_buck.cir` (made by `mkr8.py`): R7's BSIM4 power stage, drivers, dead-time generator and comparator, no ramp, no error
amplifier; hysteresis by resistive positive feedback from q (Rh1 100k / Rh2 1k: 50 mV window at the comparator); the comparator output drives a native
`Cq` through `Rq` so the hysteresis memory is a native state. fsw 395.5 kHz (settled `.tran` 2.5266 us; `.pss` 2.5284 us, the restart map's orbit is
0.07 % longer), Vout 1.7055 V. **Choice of the trigger**: `trigger=cx` (Cout's positive node, +-1.9 mV ripple, one rising crossing of 1.7055 V per
period), not `q`: the comparator output is hidden state, a native-state restart resets it, and a crossing of `q` then re-occurs 20 ns after the
restart (residual 0.3-0.9, no convergence in 21 iterations); `cx` converges in 10 (90 s, matrix-free).

`.pac` (generic, 39 state rows from 13 native, 7271 steps, prepass 3.6 s, PAC 1.8 s for 3 points; whole job 2 min) vs `tools/pac_perturbed.py --freerun`
(independent: two long nonlinear `.tran`, Hann-windowed single-bin DFT, amplitude halving; 0.001 V at Vrefsrc / 0.01 V at Vin, stabtime 300 us, reltol
1e-4; 8-11 min per point; FRA does not apply to a free-running converter: no fixed clock to lock to). f/fsw = 0.05, 0.13, 0.20:

| case | f | engine | pac_perturbed (halving) | difference |
|---|---|---|---|---|
| ctrl-to-out | 20 kHz | 7.8984 dB / -0.266 deg | 7.8977 / -0.2728 (0.009) | +0.0007 dB / +0.007 deg |
| | 50 kHz | 8.0013 / -0.856 | 7.9969 / -0.8614 (0.009) | +0.0044 / +0.006 |
| | 80 kHz | 8.1760 / -1.439 | 8.1740 / -1.4905 (0.008) | +0.0020 / +0.052 |
| line-to-out | 20 kHz | -21.0471 / -0.153 | -21.0533 / -0.1610 (0.003) | +0.0062 / +0.008 |
| | 50 kHz | -20.9446 / -0.481 | -20.9651 / -0.4295 (0.028) | +0.0205 / -0.051 |
| | 80 kHz | -20.7708 / -0.818 | -20.7652 / -0.7912 (0.007) | -0.0056 / -0.027 |

Bound in `check_r8_freerun`: 0.1 dB / 0.5 deg. Red proof `PAC_NATIVE_HISTORY=1`. R5 + a charged diode behind 1 Meg (generic path) vs the native R5 deck: 0.014
dB / 0.03 deg (bound 0.05 / 0.5; the old "free-running + charged diode must REFUSE" test became this). The `PAC-TD` warning "f/fsw >= 0.10 UNVALIDATED"
still prints: its band is R5's, this card adds R8 up to 0.20.

### 4. Ideal-switch (saltation) events with model charges on the orbit

The refusal is replaced by the jump itself in `pac_step_generic`: Fang's jump (`delta_x+ = delta_x- - (f+ - f-) dt_e`, event row / coefficient exactly as the native
branch: `coeff = (y(posCntrl) - y(negCntrl)) / h_dot`, `dfdt_pos/neg = -dG v_sw / C_local`) is added to the solved `y_k` at the event's capacitor-node or
inductor-branch rows BEFORE `Q_k = C_k y_k` is formed, so `Q_k` carries `C_k * jump` and `T_k = ag0 (Q_k - Q_{k-1}) - a1 T_{k-1}` its derivative
history, as in the native branch. It divides by the NATIVE local capacitance, so it is only right where the loaded charge matrix at the event rows
is the native capacitors/inductors alone: `pac_gen_salt_guard` compares the |C| row and column sums at `posRow/negRow` with the native expectation
and REFUSES otherwise (charged diode at the switched node: still refused, test kept; mutual inductance on the branch). The trajectory/adjoint is still
refused when events were seen (`.pxf` falls back to the forward engine, `.pnoise` refuses). Native circuits (including chargeless diodes next to
switches) are untouched.

Checks (`check_salt_charge*`): `pac_r6ob_ctrl2out.cir` = R6o (ideal-sw buck, 500 kHz) with a BSIM4 (capmod 2, W 100u) source follower plus 2 pF on its
gate in the control path (Rgs 1 Meg: F = 0.814 / -14.8 deg at 20 kHz, 0.508 / -52.9 deg at 100 kHz). The follower is LTI at its constant bias and
decoupled from the converter, so H(R6ob) = F(f) H(R6o) exactly; the reference is `F(.ac)` x the independent `pac_perturbed` value of R6o at the same
control DC (0.5285212 V; `.pss` backend, halving and step refinement <= 2e-4 dB): 20 kHz -9.548 dB / 175.27 deg vs -9.591 / 175.28, 50 kHz -26.878 /
165.07 vs -26.921 / 165.08, 100 kHz -40.765 / 159.78 vs -40.807 / 159.79: **+0.043 dB, <= 0.01 deg** (bound 0.1 dB / 0.5 deg). Red proof
`PAC_NATIVE_HISTORY=1`: -37.6 dB vs -9.6 dB (28 dB off). A native RC control path through the SAME generic code (`PAC_FORCE_GENERIC=1`) matches the
exact factorisation to 0.002 dB, native R6o forced through the generic path to 3e-4 relative; the MOS1 follower deck is +0.078 dB, a BSIM4 follower +0.034
/ +0.042 / +0.054 dB at Ibias 30 / 100 / 300 uA. So there is a small constant (in frequency) gain offset specific to model charges in the control path,
unexplained (not step size: reltol 1e-3..1e-6 and pacpoints 400..4000 do not move it; not the hidden-state restart: a native 2 pF on the gate does not
change it); it is 1/10 of the project target and is bounded, not root-caused. `pac_perturbed` cannot referee R6ob itself: its `.pss` backend restarts the
hidden follower state each period (phase 6.9 deg off at 100 kHz, which the exact factorisation proves wrong) and its `.tran` backend does not settle
inside its periodicity bound.

### Not done / limits (2026-10-02)

- `.pnoise` of R7 (and any switched circuit) absolute level, sideband folding of flicker: no independent reference.
- Free-running: the period-variation projection on a generic state set; `.pnoise`/`.pxf adjoint` on a free-running orbit; `.pstb` on R8 (no probe in the deck);
  more than the 20-80 kHz points (0.25 fsw = 99 kHz not run); the PSS stability verdict of R8 at the tight profile; the hidden-state restart artefact of
  the PSS orbit itself (0.07 % in the period, polished away for the PAC but not removed).
- Saltation: events at a node that also carries model charge (refused), adjoint/noise with events, the unexplained +0.04-0.08 dB offset.
- The R7p memory of the krylov PSS (1.4 GB growth) was measured, not investigated.
- FD `.pnoise` is not a reference for narrow conduction windows (harmonic truncation).

### Regression (OMP 2, one suite at a time; logs `~/Work/data/ngspice_augmentation/544/port-544/`, baselines `port-s96/`)

| suite | s96 | 544 | note |
|---|---|---|---|
| test_pac.py | 32/32 | 32/32 | |
| test_pxf.py | 15/15 | 15/15 | |
| test_pnoise.py | 5/5 | 7/7 | + MOS1/BSIM4 track-and-hold, + its red proof |
| test_pac_charge.py | 28/28 | 45/45 | + 7 LTI noise decks, 4 quasi-static noise decks, 2 R8, 2 saltation + R5 charged-diode, 3 red proofs; ASan clean |
| test_r7.py | 16/16 | 18/18 | + R7p `.pstb` and `.pac` (sparse, KLU, vs R7, vs FRA, RSS) |
| test_pss.py | 44/48 (4 info) | 44/48 | |
| test_pss_krylov.py | 27/27 | 27/27 | |
| test_r3d_r4d.py | 21/25 (4 info) | 21/25 | |

pac/pxf/pnoise/pac_charge(first run) used build 544-12; the one failure of that run (the charged-diode-at-the-switch-node refusal: the inductor-routed
jump let a diode cj at the switch node through the guard, which only looked at the jump rows) was fixed by also checking the switch's power nodes
(544-13) and all suites listed after it, and pac_charge again, ran on 544-13 (pac_charge 45/45). ASan: the existing chain check plus, on the new
paths, the four quasi-static noise decks, MOS1/BSIM4 track-and-hold, R6ob saltation, R5 + charged diode (free-running) and the R8 `.pac`: no report.

## Period restart on transistor-level circuits, and the three "small residuals" (card #546, 2026-10-02)

Branch aug-546 on aug-544 (s99). Build `aaba50ba3f49-546-3` (ASan `aaba50ba3f49-546-asan`). Two engine changes, both small:

1. **Hidden-state continuation of the `.pss method=shoot` period restart** (`pssshoot.c`, hook in `dctran.c`; `PSS_NO_HIDDEN=1` switches it off).
2. **`pacreltol=r`** on `.pac/.pstb/.pxf/.pnoise` (RELTOL of the td recording passes only; `inp2dot.c`, `psssetp.c`, `pssdefs.h`, `pacshoot.c`). Opt-in: unset = byte-identical.

### 1. Cause common to residuals 1 and 3: the shooting engine restarted every period run from the native C/L states only

`pss_run_period()` starts each period run with UIC from the native states. Every other node (gate chains, logic nodes, the follower's source
node, nodes carried only by device-model charge) is re-initialised (the node vector the run inherits, device ics derived once from zero node
ics), so every run begins with a startup transient that the periodic orbit does not have, and the Newton fixed point is the fixed point of the
RESTARTED map, not of the circuit. Card #542 had fixed this only for the recording runs of `.pac/.pstb` (`g_pac_xinit`); the `.pss` orbit
itself (Newton, Floquet, the emitted period, and the x0 handed to `.pac`) still had it.

Fix: the full node vector at the end of the last *latching* run (the warm-up transient stages and every Newton baseline period, `record=1`) is
put back into `CKTrhs` before the next run and the device ics are re-derived from it with `DEVsetic` (native C/L excluded: their ics come from the
shooting state), exactly what `pac_run_period()` does. The finite-difference / Krylov / trial runs of one Newton iteration do not latch, so all
perturbations of an iteration start from the same hidden vector and the Jacobian stays self-consistent; the hidden vector follows the baseline
(a lagged-hidden block iteration, no extra Newton iterations: R7 iter 4 as before). Only circuits with a non-native device that is not a
diode take it (anything but R, C, L, V, I, controlled sources, sw/csw, ASRC, mutual, diode); native and diode decks run the old code bit for
bit (checked: `pss_r6o_shoot`, `pss_r3_shoot`, `pss_r4_shoot`, `pss_r2_shoot`, `pac_r3_ctrl2out`: every output file and the log identical
to s99 except the raw-file Date/Command header lines and the output path in the log; `~/Work/data/ngspice_augmentation/546/bitid.sh`).
Chargeless/charged diodes keep the old restart (not investigated: out of the card, and the R4d/R5 numbers stay bit-identical).
`.pss trigger=` free-running orbits (`PSSshootFreerun`) are not changed.

#### Residual 1: R7 orbit offset (+0.36 mV, 0.57 % of ripple): ROOT-CAUSED AND FIXED

| orbit (`r7_pss_shoot.cir`, reltol 1e-4) vs settled tran (1 ns grid, reltol 1e-5, `tests/fixtures/r7_tran.cir`) | s99 | 546-3 |
|---|---|---|
| mean(vout) offset | +0.362 mV (+0.516 % of ripple 70.15 mV) | -0.004 mV (-0.006 %) |
| mean(il) offset | +0.484 mA (+0.050 %) | -0.005 mA (-0.001 %) |
| max pointwise vout, t >= 50 ns | 0.567 % | 0.007 % |
| max pointwise il, t >= 50 ns | 0.102 % | 0.002 % |
| first 15 ns, vout / il | 0.65 % / 0.13 % | 0.01 % / 0.00 % |
| opout (restart spike 47 % at 1 ns, 8 % at 50 ns) | not comparable | <= 0.33 % everywhere (1 ns reference grid) |
| `.pss` iterations / rho / wall (OMP 2) | 3-4 / 0.893 / 25-30 s | 4 / 0.8929 / 18-28 s |

Tolerance sweep on both sides (ramp repaired with `tf = 0.03*tclk-1n`, the #541 trap, since 1e-6 aborts the old PULSE; mean(vout) offset of the
orbit minus the last period of the tran, `~/Work/data/ngspice_augmentation/546/cmp.py`):

| orbit \ tran | 1e-5 | 1e-6 | 1e-7 |
|---|---|---|---|
| 1e-5 | -0.0041 mV | -0.0042 | -0.0042 |
| 1e-6 | -0.0040 | -0.0041 | -0.0042 |
(1e-4 on the old PULSE vs the old tran: -0.0040 mV.) Every cell: max pointwise 0.007 % of ripple (vout), 0.002 % (il). Neither side moves, so the
reference was NOT the limit: the whole 0.5 % was the restart. New test bounds (`tests/test_r7.py check_pss_base`): vout/il pointwise 0.03 % (was 1 %),
mean offset 0.02 mV (was 1 mV), opout 1 % (was "not compared"), t >= 2 ns (was 50 ns); red proof `PSS_NO_HIDDEN=1` must give > 0.2 mV.

#### Residual 3: BSIM4/MOS1 source follower in an ideal-switch control path (+0.043 / +0.078 dB): same cause, FIXED

The earlier reading "not the hidden-state restart (a native 2 pF on the gate does not change it)" was wrong about what the restart does. With s99 the
`.pss` orbit of `pac_r6ob_ctrl2out.cir` was not an orbit at all: the follower's gate state at the fixed point of the restarted map was 1.1499 V instead
of its bias 1.2163 V, ctrlnode swung 0.395..0.511 V instead of sitting at 0.52852 V, and vout was 2.959 V instead of 3.171 V (-7 %). `.pac` then
linearised along a recorded period in which the follower was relaxing, which gave a constant gain error small enough to look like 0.04 dB.

| `pac_r6ob_ctrl2out.cir` | s99 | 546-3 |
|---|---|---|
| follower orbit: v(g) / v(ctrlnode) | 1.150..1.195 / 0.395..0.511 V | 1.216300000 / 0.528523247 V, constant |
| vout at the orbit start | 2.9591 V | 3.1715 V |
| `.pac` v(ctrlnode) / `.ac` v(ctrlnode), 100 Hz / 20 kHz | 0.99927 / 0.99925 | 1.0000000 / 0.99999 |
| ctrl-to-out 20 / 50 / 100 kHz vs F(.ac) x pac_perturbed(R6o) | +0.043 dB | -0.0015 dB (-9.593 vs -9.591 at 20 kHz; worst 0.002 dB, <= 0.00 deg) |
| identity H(R6ob) = F(.ac) H(R6o, same control DC), 100 Hz..100 kHz | 5.4e-3 | 4.9e-4 |
| MOS1 follower identity (`pac_r6om_ctrl2out.cir`) | 9e-3 (0.078 dB) | 2.0e-4 |

Red proof (`tests/test_pac_charge.py`): `PSS_NO_HIDDEN=1` gives 4.98e-3 again. New bounds: `check_salt_charge` 0.01 dB / 0.05 deg (was 0.1 / 0.5),
identity 2e-3 (was 1.5e-2), and a MOS1 identity check is added.

### 2. Residual 2: R7 td loop gain +0.09 dB over the FRA from 20 kHz up: it is the RELTOL of the recording pass, not the orbit, Points or the FRA

Controlled experiments (`.pstb`, `lin 10 20k 200k`, FRA points 20/40/80/200 kHz; diff = engine - FRA in dB, FRA unc 0.009 / 0.017 / 0.011 / 0.032):

| experiment | 20k | 40k | 80k | 200k |
|---|---|---|---|---|
| s99 (orbit offset present), default, old PULSE deck | +0.089 | +0.096 | +0.096 | +0.093 |
| 546-1 (orbit fixed), same deck | +0.0885 | +0.095 | +0.095 | +0.092 |
| ... Points 800 | +0.087 | +0.094 | +0.093 | +0.090 |
| ... pacpoints 4000 | +0.113 | +0.120 | +0.119 | +0.116 |
| gear instead of trap | not supported by the td recording pass (error 1), not investigated | | | |
| deck reltol 1e-4 (recording runs at the PSS tight 1e-5), repaired PULSE | +0.115 | +0.122 | +0.121 | +0.118 |
| + `pacreltol=` 1e-4 / 3e-5 | -0.298 / +0.185 | -0.292 / +0.192 | -0.292 / +0.191 | -0.295 / +0.188 |
| + `pacreltol=` 1e-5 / 3e-6 | +0.115 / +0.055 | +0.122 / +0.061 | +0.121 / +0.061 | +0.118 / +0.058 |
| + `pacreltol=` 1e-6 / 3e-7 / **1e-7** | +0.028 / +0.019 / **+0.014** | +0.035 / +0.026 / **+0.021** | +0.034 / +0.025 / **+0.020** | +0.031 / +0.022 / **+0.017** |
| + `pacreltol=` 5e-8 (floor; 2e-8 and 1e-8 abort the recording pass, error 1) | +0.013 | +0.020 | +0.019 | +0.016 |
| LTE only tight (trtol 0.07, reltol 1e-4) | +0.015 | +0.022 | +0.021 | +0.018 |
| NR only tight (reltol 1e-6, trtol 700) | -0.020 | -0.014 | -0.014 | -0.017 |

(not monotone at loose tolerance: the dependence is on where the LTE-chosen steps fall on the comparator/driver edges.)

- Not the orbit offset (1): +0.089..0.096 dB on s99 and +0.0885..0.095 dB with the orbit fixed (0.001 dB apart), and recording-only tightening (orbit stays at 1e-4) gives the same numbers as tightening everything.
- Not the Points grid: Points 400 -> 800 moves 0.002 dB, pacpoints 4000 moves it the wrong way (the grid is set by LTE at the edges, not by maxstep).
- Not the FRA: `tools/pac_perturbed.py` (real nonlinear transients, repaired deck, reltol 1e-4, 80 kHz, amplitude 0.02, halving 0.012 dB) gives -12.8295 dB vs
  FRA -12.8312 dB, and -12.8295 dB at pac_perturbed reltol 1e-5 and 1e-6 (tran backend 23 / 44 / 74 min at 1e-4 / 1e-5 / 1e-6): the reference is converged to 4 digits. The engine at pacreltol 1e-7: -12.8108 dB (+0.019 dB).
- It is the base trajectory's step grid: the recursion linearises exactly *that discrete trajectory*, whose edge resolution is set by the LTE
  of the base run, not by the accuracy of the sensitivity along the edge. Smallest deck: the 3-stage BSIM4 inverter chain of card #539 bug 1 (`tests/bugs/539/bug1_pac_model_charge`),
  `.pac pacpoints=400`, static reference = `.tran` derivative at reltol 1e-7 (0.51500), td gain 0.52087 / 0.51657 / 0.51536 / 0.51510 at recording reltol
  1e-4 / 1e-5 / 1e-6 / 1e-7 (+1.14 / +0.31 / +0.07 / +0.02 %; N = 1 chain: 0.05 % at any reltol): error ~ reltol^0.6, no model involved but the
  charges' edge response (`546/chain.py`).
- Fix: `pacreltol=1e-7` on the `.pstb`/`.pac` card (cost: the recording pass only, R7 `.pstb` job 22.6 s -> 26.9 s on a loaded development workstation). R7 loop gain with it: +0.014 .. +0.021 dB over the converged FRA (inside its
  uncertainty at 80 kHz +/- 0.011 and 200 kHz +/- 0.032, 0.005-0.01 dB above at 20/40 kHz). The default stays the deck's tight-PSS reltol (1e-5 for a 1e-4 deck): a deck
  that needs < 0.1 dB should set `pacreltol`. New test `check_pstb_tight` (fixture `r7_pstb_tight.cir`): |diff| <= 0.05 dB (was max(0.5, 3 unc)), red proof: without `pacreltol` > 0.05 dB.

### Not done / open (2026-10-02)

- Free-running (`trigger=`) orbits still restart from native states (0.07 % period error of R8, polished away for `.pac`, not removed): `PSSshootFreerun` was not touched.
- Diode (charged or not) decks keep the old restart: they were excluded so that R4d/R5 outputs stay bit-identical; not measured whether they carry an offset.
- `pacreltol` is not the default; the remaining +0.014..+0.021 dB over the FRA at 1e-7 is within the reference uncertainty but was not driven lower (1e-8 aborts "timestep too small").
- No FRA was re-run (16-24 min per point): the FRA-amplitude hypothesis was judged by the existing halving data and by `pac_perturbed` at three reltols (80 kHz only).
- gear integration for `.pac`/`.pstb` (fails with error 1), and the per-device `C_k = dq/dx` vs secant question: not needed for any of the three residuals, not tested.

### Regression (OMP 2, one suite at a time, redirected copy of tests/ under `546/run`, outputs `546/out`; logs `~/Work/data/ngspice_augmentation/546/port-546/`, baselines `port-s99/`)

| suite | s99 | 546-3 | note |
|---|---|---|---|
| test_r7.py | 18/18 | 21/21 | + orbit bounds 0.03 %/0.02 mV (red proof PSS_NO_HIDDEN), + `pacreltol` loop gain 0.05 dB (+ red proof) |
| test_pac_charge.py | 45/45 | 47/47 | + MOS1 follower identity, + red proof PSS_NO_HIDDEN; salt/identity bounds tightened 10x / 7.5x; ASan chain clean |
| test_pss.py | 43/48 (5 info) | 43/48 (5 info) | |
| test_pac.py | 32/32 | 32/32 | |
| test_pxf.py | 15/15 | 15/15 | |
| test_pnoise.py | 7/7 | 7/7 | |
| test_pss_krylov.py | 27/27 | 27/27 | |
| test_r3d_r4d.py | 21/25 (4 info) | 21/25 (4 info) | the two orbit-vs-tran rows need the standing r3d/r4d references (copied) and shoot output (regenerated with 546-3) |
| test_handbook.py | 109 pass, 1 FAIL (R/C-cache speed ratio 1.475 under load) | 110 pass, 0 FAIL | ASan twins clean |

Native / diode decks bit-identical: `546/bitid.sh pss_r6o_shoot.cir pss_r3_shoot.cir pss_r4_shoot.cir pac_r3_ctrl2out.cir pss_r2_shoot.cir` (s99 vs 546-1 and 546-3): every output file and
log equal (raw-file Date/Command lines and the output path filtered). ASan (`546-asan`): R6ob `.pac` (same -0.33029 as the normal build), R7 `.pstb pacreltol=1e-7`, R7 `.pss krylov`: no report.

### refs/README.md edits to make with this card (R7 / R7p section)

- table row `.pss method=shoot`: replace "Orbit vs settled tran: vout +0.36 mV mean offset (0.57 % of ripple max), il +0.48 mA (0.10 %) ... first ~15 ns of the orbit wrong (opout spike 47 % ...)"
  by "(s99) ... since card #546: vout -0.004 mV (0.007 % of ripple max), il 0.002 %, first 15 ns 0.01 %; the restart artefact is gone (hidden-state continuation)".
- table row `.pstb`: add "td loop gain +0.09 dB above the FRA at the default recording reltol; `pacreltol=1e-7` +0.014..0.021 dB (card #546)".
- "What was not done": delete "The orbit mean offset (+0.36 mV) is attributed to the period restart; the attribution was not isolated" (isolated and fixed in #546: PSS_NO_HIDDEN=1 brings it back).

## Hidden-state continuation for free-running orbits and charged diodes, and `method=gear` for the td recording pass (card #549, 2026-10-02)

Branch aug-549 on aug-546 (s101). Builds `aaba50ba3f49-549-5` (ASan `aaba50ba3f49-549-asan`). Two source files: `pssshoot.c` (items 1 and 2), `pacshoot.c` (item 3).
Switches: `PSS_NO_HIDDEN=1` (the #546 switch, now also for free-running orbits), `PSS_HID_DIODE=1` (force the continuation on every diode deck), `PSS_HID_NODIODE=1` (never on
diodes), `PAC_GEAR_TEST_ORD1=1` (red proof of the gear recursion). Not done in any item: see the end.

### 1. Free-running (`trigger=`) orbits carry the hidden state (R8)

Cause, as in #546: `PSSshootFreerun` restarted every trigger run from the NATIVE states only (every BSIM4 gate/driver/comparator node and the device charge history
re-initialised), so the Newton fixed point was that of the restarted map. Fix (`pssshoot.c`): the hidden vector belongs to the Poincare section, so the trigger probe keeps the
full node vector at the two accepted steps that bracket the crossing (same call site as the native snapshots) and interpolates it with the weight that interpolates the native states
(`g_hid_T`); every trigger run restarts from the section's vector (restored into `CKTrhs`, device ics re-derived by `DEVsetic`, native C/L excluded: `pss_run_period`, unchanged #546
code); a trigger run never latches by itself (its end is not the section); the callers latch the interpolated vector explicitly: the Newton baseline `P(x0)` (the hidden state of the
next iterate: a lagged block iteration, no extra runs), an accepted damped trial, each polish iteration of the `.pac` hand-off, and the end of the warm-up. Perturbation / Krylov /
Floquet runs all start from the same vector (non-latching), so the finite-difference Jacobian stays self-consistent. Taken only on the interpolating path (`pss_trig_hidden_state`,
i.e. exactly the circuits #544 already interpolates) with a charge-carrying non-native device; the restart-refine path (R5 and other native decks) is untouched and bit-identical
(R5 native: `pss_r5*`/`pac_r5*` untouched code path). Dense and matrix-free use the same `pss_run_until_trigger`.

R8 (`tests/fixtures/r8_pss_1e6.cir`, `r8_tran_1e6.cir`; helper scripts `549/r8/gen.py cmp8.py`), `trigger=cx`. Reference: the settled `.tran` (2 ns grid, last 20 periods of 250..300 us),
period = mean spacing of the rising crossings of cx through 1.7055 V; the tran period depends on its own reltol (the reference is the limit below ~2e-5):

| settled .tran reltol | 1e-4 | 1e-5 | 1e-6 | 1e-7 | 1e-8 |
|---|---|---|---|---|---|
| period | 2.530351 us | 2.530948 | 2.531095 | 2.531153 | aborts ("timestep too small" at 2.2 us) |

(the reltol-1e-3 `.tran` of card #544, 2.5266 us, was a different tolerance, not a different circuit.)

| `.pss trigger=cx` (deck reltol; the engine caps at 1e-5) | s101 (native restart) | 549 (hidden continuation) |
|---|---|---|
| deck 1e-4 / 1e-6 / 1e-7: period | 2.528398 / 2.528620 / n.r. us | 2.530886 / 2.531100 / 2.531151 us |
| period vs the `.tran` at the SAME reltol (1e-5 / 1e-6 / 1e-7) | -9.2e-4 / -9.8e-4 / n.r. | **-2.0e-5 / +1.3e-6 / -6.1e-7** |
| orbit vs `.tran` 1e-6, max pointwise, % of ripple: vout (95.6 mV) / il (1.32 A) / cx (3.75 mV) | 0.816 / 0.896 / 2.237 | **0.012 / 0.013 / 0.034** (deck 1e-6) |
| same at the deck / tran 1e-7 | n.r. | 0.016 / 0.017 / 0.034 |
| mean(vout) offset | -0.041 % of ripple | +0.000 % |
| Newton iterations / wall (krylov, loaded development workstation, OMP 2) | 10 / 141 s (1e-4); 4 / 111 s (1e-6) | 8 / 111 s (1e-4); 6 / 134 s (1e-6); 3 / 139 s (1e-7) |
| dense (`NGSPICE_PSS_KRYLOV=0`, deck 1e-6): period / orbit % / iterations / wall | not run (interrupted to stay within the ngspice limit) | 2.531102 us (+2.0e-6) / 0.011, 0.012, 0.034 / 3 / 391 s |
| Floquet radius (Arnoldi, whole spectrum, 5 of 13) | 0.7235 (1e-4), 0.7285 (1e-6) | 0.7319 (1e-4), 0.7316 (1e-6), 0.7327 (1e-7) |

Target met: orbit <= 0.1 % of ripple (0.012 .. 0.034 %), period <= 1e-5 relative at equal tolerance (1.3e-6, 6e-7); at unequal tolerance the tran reference is the limit (its period moves
2.3e-5 relative from 1e-6 to 1e-7, the `.pss` 1e-6 deck is 2.1e-5 from the 1e-7 tran). The remaining 0.012-0.017 % is probably the linear interpolation of the 2 ns reference grid (not
isolated). The Floquet radius moved 0.005-0.009 (no independent reference was built; the dense estimate is "not computed: estimates at h and h/2 disagree" (0.574 / 0.410) on the new
build).

`.pac solver=td` on the new R8 orbit vs `tools/pac_perturbed.py --freerun` (`check_r8_freerun`, the #544 reference table; the reference does not depend on the engine's orbit):
ctrl2out worst **0.0018 dB / 0.050 deg** (was 0.0044 / 0.052), line2out **0.0179 dB / 0.054 deg** (was 0.0205 / 0.051); bound tightened 0.1 dB / 0.5 deg -> 0.05 dB / 0.2 deg (reference halving spread <= 0.028 dB).

Test: `tests/test_r7.py check_r8_orbit` (period 5e-6, orbit 0.05 % vout/il, 0.1 % cx, vs the reltol-1e-6 tran) and `check_r8_orbit_red` (`PSS_NO_HIDDEN=1` must give > 5e-4 and > 0.3 %).

### 2. Diode decks

What was measured: the discharge diode of the dead-time network (`.model ddt D(IS=1e-6 N=1.0)`, R3d/R4d as shipped) has NO charge (no CJO, TT): there is no hidden diode state to carry. Two experiments:

**(a) as shipped (chargeless diode), `PSS_HID_DIODE=1` (hidden path) vs default (restart)**, shoot reltol 1e-5 vs settled `.tran` reltol 1e-6 (last period, 2000 phase samples; `549/rd/cmprd.py`):

| | R3d restart / hidden | R4d restart / hidden |
|---|---|---|
| vout (4.8 / 4.7 mV ripple), max % of ripple | 0.0040 / 0.0040 | 0.0026 / 0.0018 |
| i(L1) | 0.0010 / 0.0011 | 0.0024 / 0.0015 |
| compensator (opout 12.5 mV / vc 1.33 mV) | 0.041 / 0.041 | 0.089 / 0.057 |
| Floquet radius at deck 1e-5 | 0.980438 / 0.980408 | 0.994031 / **-1.1649 "UNSTABLE(period-doubling)"** (0.99403 at deck 1e-6 and 1e-7: FD noise of the 1e-5 monodromy, not the orbit; the settled tran is periodic to 1e-7 of ripple) |
| `.pac/.pstb` vs converged FRA (`tests/test_r3d_r4d.py` cases; worst dB) | R3d ctrl2out .044, line2out .077, zout .022, pstb .360; R4d ctrl2out .071, line2out **.168**, zout .022, pstb .123 | R3d .029, .076, .030, .345; R4d .074, **.319**, .022, .113 |

The orbit is equal within the reference. The hidden-vector orbit is WORSE for R4d line2out: a constant +0.23 .. +0.26 dB offset over 500 Hz .. 200 kHz (restart: 0.00-0.02 dB; FRA unc 0.005 dB at 500 Hz).
Reason: a deck without model charge takes the NATIVE `.pac` recording path (`pac_circuit_has_model_charges`: the detection run finds no charge), which restarts from the native states only (hidden nodes
re-initialised, bit-identical to before); it closes the orbit that the RESTARTED `.pss` map converged to, not the hidden-vector one (a -57 dB cancellation quantity is sensitive to the 1e-5 x ripple
differences). With a charge anywhere on the orbit `.pac` takes the generic path with the periodic-continuation start (`g_pac_xinit`), consistent with the hidden-vector orbit. Hence the rule below.

**(b) the same decks with a charged diode (`CJO=1p TT=1n`)** (restart = `PSS_NO_HIDDEN=1`; hidden = new default), max % of ripple vs `.tran` 1e-6 / 1e-7:

| | R3d restart -> hidden | R4d restart -> hidden |
|---|---|---|
| vout | 0.0875 -> 0.0037 / 0.0037 | 0.0022 -> 0.0014 / 0.0013 |
| i(L1) | 0.0693 -> 0.0009 / 0.0009 | 0.0021 -> 0.0005 / 0.0007 |
| compensator | **3.24 -> 0.031** / 0.018 | 0.113 -> 0.025 / 0.019 |
| Floquet radius | 0.980344 -> 0.980325 | 0.994036 -> 0.994031 |

(`.pac`/`.pstb` on the charged decks: 0.04-0.34 dB (hidden) vs 0.02-0.31 dB (restart) against the FRA of the UNcharged decks, which is not a reference for them: not decisive, no charged FRA exists.)

**Decision.** A diode with charge (CJO, CJSW or TT nonzero in any diode model) is a charge-carrying device exactly like a MOSFET and takes the continuation by default (`pss_diode_charged`, the test
`pss_trig_hidden_state` already used for free-running orbits). A diode without charge keeps the old restart bit for bit: its "hidden state" is empty, the hidden-vector orbit is no better, and it
breaks the self-consistency with the native `.pac` pass (R4d line2out +0.24 dB). Bit identity (s101 vs 549-4, all output files and the log): `pss_r6o_shoot`, `pss_r3_shoot`, `pss_r4_shoot`, `rd_r3d_shoot`,
`rd_r4d_shoot`, `pac_r3_ctrl2out`, `rd_pac_r4d_line2out` (`549/bitid.sh`). No existing pinned number changes for that reason. Test changes (`tests/test_r3d_r4d.py`): new `check_charged_diode`
(R3d, R4d; bounds 0.01 % vout/i(L1), 0.1 % compensator; red proof `PSS_NO_HIDDEN=1` must exceed the compensator bound: 3.2 %); equal-tolerance orbit check bounds tightened per signal from 0.1 % to
0.01 % (vout, i(L1); measured <= 0.0042 %), 0.1 % (opout) and 0.15 % (R4d vc, measured 0.089 %: tolerance-limited, the 1e-5 / 1e-6 / 1e-7 trans differ by 0.23 / 0.089 / 0.02 % of its 1.3 mV ripple), and
added for R4d (it had none); the settled-tran decks now write only the last 25 us (the full file was 0.8 GB for R3d).

### 3. `method=gear` for the td recording pass

Cause: `PSSpacRecordStep` refused every step whose integration method was not trapezoidal order 1-2 ("analytic propagation only implements TRAPEZOIDAL ... abandoning this pass", then the PSS driver's
generic "error 1"); the recursions (native per-C/L and generic) were written for the trapezoidal derivative history `T`. Fix (`pacshoot.c`): gear is carried by the GENERIC charge recursion, which has no
integration-specific history other than `Q`: a BDF/Gear step of order p gives `dq/dt|k = ag[0] Q_k + sum_j ag[j] Q_{k-j}`, so `J(ag0) y_k = -sum_{j=1..p} ag[j] Q_{k-j}` (+ the forcing); the charge history
`Q_{k-1..k-6}` is kept (`g_pac_Qprev`, `g_pac_Qh[2..6]`), the order of the step is the highest nonzero `CKTag[j]`. Order 1 is backward Euler (same formula); trap keeps its `T` recursion bit for bit.
Native-only decks are routed through the generic path when the method is gear (`force_gen` in `pac_generic_prepare`: prepass + affine-in-ag0 check cost 7 ms on R3). Refused with a message (not an error 1): `.pnoise`
and `.pxf adjoint=1` (the transposed sweep replays the trapezoidal history) and an ideal-switch (saltation) event (the jump is derived for the trapezoidal history). `.pxf adjoint=0`, `.pac`, `.pstb` work.

Validation (gear vs trap, `pacreltol=1e-7`): R3 ctrl2out 20..180 kHz, 5 points: worst **0.0013 dB / 0.0035 deg** (at the default tolerance 0.014 dB; trap itself moves 0.04 dB between the default and 1e-7:
the spread is the recording pass's tolerance). R7 `.pstb` 20..200 kHz, 10 points: **0.0047 dB (constant) / 0.002 deg**; the gear loop gain is +0.019 .. +0.026 dB over the converged FRA (trap +0.014 .. +0.021).
Tests (`tests/test_pac_charge.py`): `check_gear` (0.01 dB / 0.05 deg), `check_gear_r7` (0.02 dB / 0.05 deg), `check_gear_refusals`, red proof `PAC_GEAR_TEST_ORD1=1` (history beyond Q_{k-1} dropped -> NaN curve). 

### Regression and what was not done

(OMP 2, one suite at a time from the main repo `tests/`, build 549-5 with ASan twin 549-asan; logs `~/Work/data/ngspice_augmentation/549/port-549/`, baselines `port-s101/`, handbook `port-s99/`)

| suite | s101 | 549-5 | note |
|---|---|---|---|
| test_pss.py | 43/48 (5 info) | 43/48 (5 info) | |
| test_pac.py | 32/32 | 32/32 (2 info) | |
| test_pxf.py | 15/15 | 15/15 | |
| test_pnoise.py | 7/7 | 7/7 | |
| test_pac_charge.py | 47/47 | 51/51 | + gear (R3, R7 pstb, refusals, red proof); R8 bound 0.05 dB / 0.2 deg; ASan chain clean |
| test_r7.py | 21/21 | 23/23 | + R8 orbit vs tran at equal tolerance, + its red proof (`PSS_NO_HIDDEN=1`) |
| test_pss_krylov.py | 27/27 | 27/27 | |
| test_r3d_r4d.py | 21/25 (4 info) | 27/31 (4 info) | + charged-diode orbit R3d, R4d + red proof, + R4d equal-tolerance orbit (3 signals), per-signal bounds tightened |
| test_tran_from_pss.py | 6/6 | 6/6 | |
| test_handbook.py | 110/110 (s99 log) | 110/110 | ASan twins clean |
| test_loss.py | 34 pass (1 info) | 34 pass (1 info) | |

Bit identity s101 vs 549-4 (outputs + log): pss_r6o_shoot, pss_r3_shoot, pss_r4_shoot, rd_r3d_shoot, rd_r4d_shoot, pac_r3_ctrl2out, rd_pac_r4d_line2out. ASan (549-asan): R8 `.pss` (same iterations/period as the normal build), gear R3 `.pac`, gear R7 `.pstb`,
charged-diode R3d `.pss`: no report. The one test bug found on the way: `_variant_deck`'s `save` list lacked `v(vc)`, so the R4d settled tran wrote nothing (fixed).

Not done: (1) the dense free-running Floquet estimate (disagrees at h / h/2, "not computed") was not investigated; the hidden vector is frozen in every perturbation run, so the monodromy omits the
sensitivity of the hidden state at the end of the run to the start state (it is exact when the hidden state relaxes within the period: R7 0.8929 vs 0.893 before and after) -- a verdict for a deck with
slowly relaxing hidden state (diode TT comparable to the period) is not independently validated; (2) no independent reference for `.pac`/`.pstb` on the charged-diode decks (FRA 16-24 min per point not run);
(3) gear: orders above 2 were not seen in the runs (ngspice's default maxord is 2; `.options maxord=` untested), no gear `.pnoise`/`.pxf adjoint=1`, no gear with saltation; (4) R8 `.pss` period below ~2e-5 relative
is unprovable against the tran reference (its own tolerance dependence); the 0.012 % residual orbit error was not attributed.

### refs/README.md edits to make with this card

- R8 section: replace "(settled `.tran` period 2.5266 us; `.pss trigger=cx` 2.5284 us)" by "(settled `.tran` period 2.5311 us at reltol 1e-6, 2.5312 at 1e-7 -- the 2.5266 of card #544 was reltol 1e-3; `.pss trigger=cx` 2.5311 us, hidden-state continuation of the trigger runs, card #549; before: 2.5284 us)" and "the restart map's orbit is 0.07 % longer than the settled tran" -> "was 0.1 % shorter".
- the `.pac` line: "<= 0.021 dB / 0.052 deg" -> "<= 0.018 dB / 0.054 deg (card #549)".

## XSPICE code models under `.pss method=shoot` and td small-signal analyses (card #551, 2026-10-02)

Branch aug-551 on aug-549 (s103), commits `4fdd73f20`, `3ff606356`, `4220f339f`. Builds `aaba50ba3f49-551-6` (ASan twin `aaba50ba3f49-551-asan`). Suite `tests/test_xspice.py`, fixtures `tests/fixtures/xspice_*.cir`, `r9_*.cir`, reference `refs/r9_digital_buck.cir` (generator `tools/gen_r9.py`). Switch (red proof only): `PSS_XS_NO_STATE=1`; `PSS_XS_TRACE=1` prints every period run's restore/latch; `PSS_XSPICE_ANY=1` lifts the model whitelist (experiments only).

### What was wrong on s103 (measured, `~/Work/data/ngspice_augmentation/551/dev/inv_s103.txt`)

Minimal decks, `.pss method=shoot` orbit vs the settled `.tran` (% of the tran ripple), and td `.pac` vs plain `.ac`:

| model class | s103 `.pss` | s103 td `.pac` |
|---|---|---|
| stateless analog (gain, summer, mult, divide, limit, climit, ilimit, pwl, aswitch, pswitch) | right (<= 0.003 %) | **wrong by the model's DC rhs**: gain=2.5 out_offset=0.1 gave 2.6 (4 %), summer 4.4 %, mult 22 %, divide 57 % |
| s_xfer, int, d_dt, slew, hyst, d_dt | **converged and wrong**: s_xfer 4.9 % (271 % with a fast pole), int 153 %, d_dt 250 %, square (clock) 8 % / 100 % | s_xfer `inf`, d_dt 1.6e4 x, int "REFUSED: not affine in ag0", hyst 90 % |
| digital / hybrid (adc_bridge, dac_bridge, d_tff, ...) | **SEGV in the 2nd period run** (`cm_adc_bridge`, `OUTPUT_STATE` through a stale output pointer) | not reached |
| delay (analog delay line) | 86 % wrong | - |

Root causes (all three had been guessed at in the #494 notes, none located):

1. **td PAC `B0` capture** (`pacshoot.c pac_capture_b0`): `CKTload(DCOP|INITSMSIG)` then `CKTacLoad`. `MIFload` (the XSPICE `DEVacLoad`) takes its analysis type from `g_mif_info.circuit.anal_type`, which only `CKTload` and `acan.c` set; the DC value was still there, so every analog code model stamped its DC rhs (`out_offset`, `output - partial*input`) into the small-signal source vector. Fix: set `MIF_AC` around `CKTacLoad`, restore after. Stateless analog models are then exact (1e-12 vs `.ac`).
2. **Hidden XSPICE state at the period restart**: a shooting period run is a fresh `DCtran(restart=1)`. (a) Analog code-model state (`cm_analog_alloc`/`cm_analog_integrate` ranges and the old input of every port in `CKTstates`) was whatever the previous run (trial, perturbation) ended with: the Newton fixed point was that of the restarted map. (b) The event system is only built once per job (`EVTsetup` from `CKTdoJob`): a second `DCtran` ran with output pointers still aimed at the previous run's event structs (the crash), stale queues, stale node lists, `inst->initialized` still true.
3. The `TIME == 0.0` initialisation branches of the models (d_dt, slew, d_pwm, d_tff, dac_bridge, ...) reset the state at every restart even if it had been carried.

### What changed

* `evtsetup.c`: `XSsnap` (take at the end of a latching run, put before every restarted run). Contents: analog `CKTstate0` ranges of every MIF instance, resolved + staged value of every digital node, the tail state block of every event instance (the previous-state chain is rebuilt by `EVTcreate_state`), pending instance and output events (times relative to the end of the run), the model-time offset. `put` rebuilds the event framework exactly as `EVTsetup` does for the queues and output pointers but keeps `initialized`, empties the node history lists (reseeded by `EVTop_save`) and compacts the state lists as `EVTaccept` does (first version put the single remaining state into the free list: models then saw `clk_old == clk`, no edge ever, caught by the divider deck).
* **Continuing model time** (`g_xs_toff`, `mif.h`): the models see `TIME = CKTtime + toff`; `toff` advances by the run length at every latching run. Their stored change times / oscillator phases stay valid, the `TIME == 0` branches run only in the very first run, `cm_event_queue`, `cm_analog_set_temp_bkpt`, `cm_analog_set_perm_bkpt` subtract it (a rounding error of the difference of two large numbers landed 1e-20 below the current event time and the event was refused: the adc_bridge never fired in restarted runs). The t=0 pass of a restarted run (DC mode, models still see time 0 there) lets models re-initialise from their inputs; `dctran.c` calls `PSSxsPostInitHook` right after it to put the carried state back (`pssshoot.c pss_xs_post_init`).
* `pssshoot.c`: `pss_run_period` restores before / takes after (a latching run = the Newton baseline and the warm-up stages, like the #546 node-vector hidden state; perturbation, Krylov and trial runs never latch, so a Jacobian is self-consistent). **Hidden-state periodicity check** after convergence: up to 12 latching baselines from the converged native state; digital node values must be equal at the start and the end of the period (else REFUSED: limit cycle / subharmonic of the period) and the analog node vector must have returned (> 1e-2 refused, > 1e-3 warned); the line `XSPICE hidden state over one period ... digital node mismatches N, analog node drift X` is printed.
* **Model policy** (`PSSxspicePolicy`, key = `Spice_Model_Name`, which differs from the directory for `capacitoric`, `inductoric`): class 1 stateless analog (orbit + small signal), class 2 analog with state (orbit only), class 3 digital/hybrid (orbit, driven only), everything else REFUSED naming the model, so are class 2/3 models in td `.pac/.pstb/.pxf/.pnoise` (the per-step linearisation has no recursion for `cm_analog_integrate` state, events are discontinuities, a quantiser has no small-signal gain) and every XSPICE model in a free-running `trigger=` analysis (no defined hidden state at the interpolated crossing).

### Results (build 551-x, development workstation load 10-14, OMP 2; `tests/test_xspice.py`)

Inventory after (PSS = orbit vs settled tran % of ripple; PAC = td vs `.ac`):

| model | `.pss shoot` | td `.pac` |
|---|---|---|
| gain summer mult divide limit climit ilimit pwl | works, <= 0.003 % | works, <= 2e-12 vs `.ac` (s103: wrong -> fixed) |
| aswitch pswitch | works, 0.0002 % | works (0.909 = the exact divider; plain `.ac` gives 1.111, the XSPICE ac model is wrong, not compared) |
| s_xfer / int / slew / d_dt / hyst | works: 0.62 / 0.20 / 0.34 / 0.06 / 0.003 % (s103: 4.9 / 153 / 0.34 / 250 / 0.003; wrong -> fixed) | refused (names the model, "keeps integration state") |
| oneshot, capacitoric, inductoric | works 0.16 / 0.09 / 0.29 % | refused |
| sine, square, triangle (clock sources, `driven=1`) | works 0.0002 / 0.014 / 0.004 % (square s103: 8 %) | refused |
| adc_bridge dac_bridge d_tff d_dff d_jkff d_srff d_dlatch d_srlatch d_and d_nand d_or d_nor d_xor d_xnor d_inverter d_buffer d_pullup d_pulldown d_fdiv d_lut d_state d_osc d_pwm | works (7 decks, Points 4000: 0.001 .. 0.16 %; s103: SEGV) | refused ("event-driven / digital") |
| delay, zener, d_ram, d_genlut, d_tristate, d_open_*, all unlisted | REFUSED by name ("not validated") | REFUSED |
| any XSPICE model, `trigger=` | REFUSED | - |

The 0.1-0.6 % of the integrating models is their own first-order `cm_analog_integrate` (`cm.c`: at trapezoidal order 2 `geq` is halved and the history term replaced, i.e. backward Euler, "WARNING - This code needs to be redone"): the orbit equals a settled tran of the same step; s_xfer 0.35 / 0.14 / 0.05 % at Points 400 / 1600 / 6400 against a 0.2 ns tran. Digital decks need resolution too: the adc_bridge detects a threshold crossing at the next accepted time point (1.3 ns apart on the edges at Points 400, 0.13 % at 4000).

Red proofs (`PSS_XS_NO_STATE=1`): comparator + d_dff sample-and-hold (held bit at t=0 is not the ic) 0.02 % -> 8.5 %; s_xfer back to the restart orbit; R9 orbit 17 mV (4.6x the ripple) off with the DPWM register back at its initial code.

Unfixed / not validated: see the end of this section.

### R9: digitally controlled synchronous buck (`refs/r9_digital_buck.cir`)

See `refs/README.md` "R9" for the deck description. Summary of the validation:

| check | result |
|---|---|
| settled tran | code 556 constant from 6 ms, Vout 3.2926 V (inside the +-9.9 mV bin), ripple 3.72 mV |
| `.pss method=shoot` T = 1.2 ms (600 cycles), Points 24000, tstab 12 ms | 3 iterations, 152-250 s, rho 0.029; orbit vs the last 1.2 ms of the settled tran (100 ns): max 0.012 mV = 0.33 % of ripple, mean offset +0.005 mV, code 556 both; hidden state "digital node mismatches 0, analog drift 7e-5" |
| red proof `PSS_XS_NO_STATE=1` | register back at ic 553 each period: orbit -17.15 mV (4.6 x ripple) |
| 7-bit DPWM variant (LSB 92 mV >> 6 mV bin) | settled tran hunts with period 2 T (Vout 3.225 / 3.315 V in alternate windows); `.pss` at T does not converge (residual 0.55) -> refused |
| td `.pac/.pstb` | REFUSED, names adc_bridge / d_tff / dac_bridge |
| FRA control-to-output (fra.py transfer) | vs averaged LC model <= 0.17 dB / 3.7 deg 2-50 kHz (-1.57 dB at the 5 kHz resonance, Ron not modelled); amplitude spread 2.5..40 mV <= 0.11 dB |
| FRA loop gain (fra.py loop, 300/1000/3000 Hz, 0.5/2/5/20 mV) | amplitude-dependent by up to 100 dB (zero-error bin + 1.2 ms sampling): no small-signal loop gain exists; only 20 mV (> bin) shows the integrator describing function |
| pac_perturbed | killed after 40 min, no point: not validated |

Design rule that made the limit cycle go away (Peterchev/Sanders): DPWM step (11 bit: 5.7 mV at Vout) below the ADC zero-error bin (6 mV wide), plus an integrator step spaced longer than the LC ring (0.66 ms) -- with a step every 96 us the 11-bit variant also hunted +-1 LSB (ringing across the bin).


### Not done / limits (2026-10-02)

* td small-signal analysis through XSPICE state, events or quantisers: refused, not implemented. A stateful-analog recursion would need the model's per-step map (the integral is not a node charge; `J(ag0)` is `A + B/ag0`); digital: saltation of an event time that depends on an analog input (adc_bridge threshold, d_pwm duty) plus a zero gain through every quantiser. Use a converged FRA.
* Pending events at the period boundary are carried (output queue validated with a clock edge straddling `t = T`), but a source that is not T-periodic at the boundary (PULSE with `td` > 0: `v(0) != v(T)`) is a trap of `.pss` in general (found while testing phase shifts), not XSPICE specific.
* Models with state outside `CKTstates`/the event blocks (analog `delay`: a malloc'd buffer) are refused by name; the periodicity check caught it (drift 2.7 %) before the whitelist was cut.
* `d_source`, `d_cosim`, `d_process`, `filesource`, `pwlts`, tlines, table models, the xtraevt real/int nodes: untested, refused.
* `.pss` Krylov/matrix-free on an XSPICE deck: not run (the hooks are in `pss_run_period`, which both use); BSIM4 + XSPICE mixed decks not run.
* The hidden analog state follows the baseline by a lagged block iteration (no extra Newton iterations): a hidden state with a time constant of many periods (a digital integrator) needs a stabilisation time that settles it; the periodicity check measures the drift and refuses > 1e-2.

## PSS memory: period runs no longer keep their output (card #553, 2026-10-02)

Branch `aug-553` (one commit on `aug-s108`, patch file `~/Work/data/ngspice_augmentation/553/patch/0001-*.patch`, not yet in `patches/series`);
build `aaba50ba3f49-553-5`, ASan `aaba50ba3f49-553-asan`. Measurements: `~/Work/data/ngspice_augmentation/553/m/` (`measure.py`, one directory per run with deck, stdout, time line).
Development workstation shared (load 13-25), OMP_NUM_THREADS=2, every run under `systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=0`, `/usr/bin/time -f %M`.

### Where the memory went

Not the Krylov basis, not the hidden-state snapshots (#546/#551), not the replay records (#505), not KLU copies: after the fix the whole heap of the
R7p matrix-free job is 7 MB. It was the **output plots of the period runs**. Every period run (baseline, product, trial, detection column, PAC recording
pass) is a `DCtran()` restart, and `DCtran` opens an ordinary plot with one vector per node/branch (`.save` absent = all of them) and fills it with every
accepted step; the front end keeps a finished plot in its plot list until `destroy`, and nothing destroyed these. R7p: 7293 rows x ~430 vectors x 8 B = 25 MB
per run (`PSS_MEM_TRACE=1`: +25.05 MB RSS and malloc in-use per run, flat otherwise), 43 runs = 1.45 GB. The dense monodromy makes 2 runs per state, so
413 states = 866 runs x 25 MB: OOM at 8 GB after 167 s (the "n x n is 1.4 MB, so something per column is retained" of the card: one 25 MB plot per column).
A `.save` list of four vectors hides it (4 x 8 B x rows), which is why the earlier PEX measurements with `.save` showed 17-24 MB at any size.

### The fix (pure memory management)

`PSSdiscardRunData` (pssshoot.c) is set around every period run (`pss_run_period`, and the PAC recording pass `pac_run_period`):
* `outitf.c`: while it is set `AddRealValueToVector`/`plotAddComplexValue` keep one point per vector (the existing `savenone` mechanism, no allocation of the
  run length up front), and `OUTendPlot` drops the finished plot (`plot_discard` in postcoms.c = `destroy <plot>`) and advances `plot_num`. The plot shells were the
  second leak: ~90 KB per run at 490 states (2348 runs of a dense 1150-state analysis = 0.6 GB even with one-point vectors).
* The run that emits the converged trajectory ("Time Domain Periodic Steady State Analysis") keeps its data: it is the plot `wrdata`/`plot`/`write` read.
* Debug: `PSS_KEEP_RUN_DATA=1` restores the old retention (used by the red-proof), `PSS_MEM_TRACE=1` prints `PSS-MEM run #n <name> rss= inuse= arena= mmap=` after each period run.

User-visible change: `display`/`setplot` after a `.pss method=shoot` no longer lists the 60-2300 intermediate plots (`tran1`, `unknown2`, ...: every one of them was a
"PSS shoot (...)" run), only the emitted trajectory; its number is not always the old one (pss62 -> pss66 on pss_r3_shoot). `wrdata`, `write` (.raw), stdout (wall_time masked) are unchanged.

### Result: peak RSS (MB), `.pss` only, matrix-free (krylov=1), before = s108, after = 553-5

Decks: R3 buck + 3 extracted RC trees (`tests/test_pss_krylov.py pex_deck`, b3d4 / b4d4 / b5d4), "no .save" = the same deck without the `.save` line, R7p = `tests/fixtures/r7p_pss_shoot_krylov.cir` (413 states, no `.save`, `wrdata` of 7 vectors).

| deck | states | solver | before | after | period runs (before = after) |
|---|---|---|---|---|---|
| R3 PEX, .save | 490 | KLU / sparse | 17 / 17 | 14 / 13 | 47 |
| R3 PEX, no .save | 490 | KLU / sparse | 254 / 255 | 19 / 19 | 47 |
| R3 PEX, .save | 1150 | KLU / sparse | 19 / 19 | 16 / 14 | 47 / 60 |
| R3 PEX, no .save | 1150 | KLU / sparse | 649 / 854 | 28 / 28 | 47 / 60 |
| R3 PEX, .save | 2470 | KLU / sparse | 24 / 23 | 21 / 20 | 58 / 60 |
| R3 PEX, no .save | 2470 | KLU / sparse | 1806 / 1930 | 49 / 49 | 58 / 60 |
| R7p | 413 | KLU / sparse | 1447 / 1449 | 39 / 40 | 43 |

R7p krylov target <= 500 MB: 39 MB. Memory is now flat in the number of period runs (the 2470-state no-.save deck: 49 MB, of which the emitted trajectory is the part that grows with the states).

Dense monodromy (`krylov=0`, KLU):

| deck | states | before | after | period runs | wall before / after |
|---|---|---|---|---|---|
| R3 PEX, no .save | 490 | 3579 | 18 | 1028 | 52 / 36 s |
| R3 PEX, .save | 490 | 63 | 17 | 1028 | 50 / 36 s |
| R3 PEX, no .save | 1150 | OOM-killed at 8 GB after 69 s | 34 | 2348 | - / 179 s |
| R7p | 413 | OOM-killed at 8 GB after 167 s | 40 | 866 | - / 303 s (rho 0.8914; the Arnoldi estimate of the krylov run: 0.8935 klu / 0.8927 sparse) |

(Wall times are from a shared machine at load 13-25 and are only a check that nothing got slower; the arms were not interleaved. R3 PEX 2470 dense not run.)

### Identity

`out.txt` (wrdata of v(vout) v(sw) v(opout) i(L1), 15 digits) of every run above is **byte-identical** before/after (14 krylov runs incl. R7p KLU and sparse, 2 dense runs), the period-run
counts and rho are the same. The change touches no number the engine reads (the engine reads the circuit's own state, never the plot vectors).

### Test

`tests/test_pss_krylov.py` `memory` check: 1150-state PEX buck without `.save`, peak RSS (getrusage of the child) <= 150 MB, bit-identical orbit and period-run count versus the
run with `PSS_KEEP_RUN_DATA=1` (which must show >= 3x the RSS, else the check cannot see the fix); `--redproof` runs it with the switch and must FAIL
(measured: 648 MB > 150 MB, red-proof OK). Regression on 553-5 vs `port-s108`: test_pss 44/48 (4 info) same; test_pss_krylov 28/28 (27 + the new check); test_pss_klu 17 compared 0 FAIL; test_pac 32/32;
test_r7 23/23; test_pac_charge 51/51; test_xspice 26/26; test_r3d_r4d 27/31 same lines as s108. ASan (553-asan): the pss_krylov ASan checks (krylov on R3, R3d, unstable fallback, R5 free-running, PEX 73 states) and the test_pac multi-point ASan runs clean.

### Not done / limits

* XSPICE circuits keep per-run event-system job data (`Evt_Job_t`, one job per DCtran restart, `EVTsetup_plot`): not measured, not changed; the plot part is fixed for them as well.
* Plain `.tran`/`.ac` etc. still keep their plots (that is their output). Only period runs of `.pss method=shoot` and the td `.pac` recording passes discard.
* The 40 MB left on R7p is the emitted trajectory (7298 rows x 430 vectors = 25 MB) + base; `.save` of the needed vectors reduces it (the deck's `wrdata` list does not).
* Dense 2470 states not run; second workstation not used; no interleaved wall-time A/B (memory change only, 47-60 runs unchanged).
* Not added to `patches/series` (that is done at merge).

## `.pnoise` on switching circuits: folded flicker noise and an independent absolute-level check (card #556, 2026-10-03)

Branch aug-556 on aug-553 (s109), one engine commit `70856795e`; build `aaba50ba3f49-556-2` (ASan `aaba50ba3f49-556-1-asan`, same source). Two jobs: the flicker (1/f)
density under a fast periodic bias, which card #544 listed as "not validated, probably overestimated", and an absolute-level reference for the loop-shaped and folded part
of `.pnoise`. Everything below was run in `~/Work/data/ngspice_augmentation/556/`.

### 1. Flicker noise under a fast periodic bias

**The error.** A device noise source is `n(t) = g(t) w(t)` with `w` stationary and `S_w(f) = A/f^ef`. On a periodic orbit the time-averaged output density at `f` is

    S_out(f) = Sum_q |G_q(f)|^2 S_w(|f - q fs|),   a(t) = g(t) H(t,f) = Sum_q G_q e^{j q ws t}

(`H(t,f)` = response at `f` to a unit source at time `t`; q = 0 is the baseband term, q != 0 the sidebands that fold from `f - q fs`). The `.pnoise ... solver=td` recursion
(`pac_engine_run_adjoint_noise`) evaluated `S_g(t_k, f)` at the OUTPUT frequency for every sideband, i.e. `S_out = S_w(f) Sum_q |G_q|^2`. That is exact for white noise
(Parseval), exact when `a(t)` has no harmonics (LTI bias, which is what the #544 LTI table tests), and wrong when the orbit modulates the source or the transfer fast: a chopper
or mixer puts the flicker of a baseband source at `fs`, where `S_w(fs) = S_w(f) f/fs`, so the old value is too high by `fs/f`; the error has either sign for other `f`
(measured +10x ... -52 %, table below).

**The fix** (`pacshoot.c`, `nevalsrc.c`): for every flicker generator the sweep keeps, per output frequency and step, `z_k = u_k sqrt(S_k) e^{+j w t_k}`, `u_k` = the complex gain of the
generator's two nodes in `nu_k` (a capture hook in `NevalSrc` for the `N_GAIN` call each device's flicker term makes, one call per `1overf` slot in the same instance order;
a count mismatch falls back to the old evaluation with a warning) and `S_k` its density at `f`. The exponent `ef` comes from a second `DEVnoise` call at `2f` (a flicker density
that is not a power law in `f` also falls back with a warning; that branch is not exercised by any deck here - BSIM4 `fnoimod=1` is still a power law, `Ssi` and `Swi` both scale as `f^-ef`). The correction added to the accumulated density of that generator, to its instance total and to the
scalar spectrum is `Sum_q |G_q|^2 (rho_q - 1)`, `rho_q = (f/|f - q fs|)^ef`, where `G_q = Sum_k z_k sinc(q ws w_k/2) e^{-j q ws tc_k}` are the Fourier coefficients of the zero-order-hold
function the engine's own discretisation represents (`w_k` the trapezoid weight, `tc_k` the centre of the interval), so that `Sum_q |G_q|^2 = T Sum_k |z_k|^2 / w_k` is exactly what
the sweep summed (Parseval holds for the piecewise-constant function; printed as `head/E` by `PNOISE_FOLD_DEBUG=1`). The sum runs over `|q| <= 2048` (`PNOISE_FOLD_QMAX`) or until
the Parseval remainder is below 1e-9; the rest is added with the density of a `q^-2` tail, `(f/(Q fs))^ef/(1+ef)` (a 10 ns-edge chopper agrees with its closed form to 2.9e-8 at
`Q = 2048`, 4e-7 at 256, 4e-6 at 64). The pairing sign (`+q` with `f - q fs`) is checked by a deck whose sideband gains are asymmetric (the wrong sign is 10-13 % off).
A sideband at `|f - q fs| < 1e-3 fs` takes the density of `1e-3 fs` (S(0) is infinite; a warning is printed when such a sideband carries more than 1 % of the value): this is what
keeps an LTI deck at `f` = a multiple of the PSS period exact (fs = 1 kHz, f = 1/10/100 kHz in the #544 table: 6e-7 after, 2-3.6 % before the floor existed).

White, thermal and shot sources and every deck without a flicker source are untouched: `.pnoise` outputs are **bit-identical** (m1buck, all vectors, `cmp`; R7 non-flicker slots; the
test_pnoise / test_pac / test_pxf / test_r7 / test_xspice / test_handbook logs match the s109 logs, see below). A deck with flicker but no modulation (LTI bias) changes by 1e-7 or less.
`PNOISE_FLICKER_OLD=1` restores the old evaluation (red proof). Cost: one extra `DEVnoise` call per step and frequency for the flicker slots, `Q * N * (generators x frequencies)`
for the sums; R7 (2161 steps, 70 flicker generators, 3 frequencies): 19.9 s with the fold, 20.9 s with `PNOISE_FLICKER_OLD=1`.

**Closed-form references** (`tests/test_pnoise.py check_pnoise_fold`, fixtures `pnoise_fold_*.cir`, no engine output enters the references; exact Fourier coefficients of the PULSE
trapezoid; T = 300.15 K; tolerance 1e-5, measured <= 9e-8):

| deck | what is folded | f | closed form (V/rtHz) | engine now | old evaluation |
|---|---|---|---|---|---|
| `pnoise_fold_chop_res` (R 1k, kf af=2 ef=1, 1 mA, x +-1 square wave, 100 kHz, 1 us edges) | resistor 1/f chopped | 1 kHz | 5.6902e-9 | +8.9e-8 | x 10.3 |
| | | 30 kHz | 5.9571e-9 | +8.7e-8 | +80 % |
| | | 90 kHz | 1.2892e-8 | +7.2e-8 | -52 % |
| | | 150 kHz | 6.2414e-9 | +9.0e-8 | -23 % |
| `pnoise_fold_chop_lp` (R1 \|\| C1 -> chop -> LP -> chop, delayed 0.2 T, 30 %: asymmetric sideband gains) | complex H(t,f) | 20 kHz | 2.1947e-9 | +2.7e-8 | +51 % |
| | | 70 kHz | 1.8591e-9 | +2.3e-8 | -13 % |
| `pnoise_fold_mos` (MOS1 nlev 2, gate trapezoid 1.2 <-> 1.6 V, saturated, no caps) | flicker `kf gm^2/(f W L Cox)` with `gm(t)` modulated | 2 kHz | 4.7680e-7 | -8.7e-9 | +3.4 % |
| | | 90 kHz | 8.1151e-8 | +9.5e-9 | -9.4 % |

The thermal / channel-thermal vector of the same runs agrees with its exact value (`Sum_q |c_q|^2 4kTR`, `4kT (2/3) beta <vgst> Rd^2`) to 1.8e-6 and is bit-identical to the old evaluation.

**Time-varying bias with real device models** (`tests/test_pac_charge.py check_pnoise_qs`). The #544 check compared `.pnoise` at 1 kHz with the quasi-static time average of plain
`.noise` over 96 phases of a 100 kHz LO, and the old evaluation agreed with it. That reference is the WHITE answer: for a flicker source modulated at the LO, only the q = 0 term
carries `S_w(f)`, the harmonics of `a(t)` sit at `S_w(f - q fs)`. The check now takes the time average for everything that is not flicker and `Sum_q |a_q|^2 f/|f - q fs|` (`a_q` =
DFT of the per-phase flicker amplitude) for the flicker part, with the instance totals and the spectrum recomposed. Engine vs reference (worst of total and every per-device vector):
MOS1 3.0e-5 (no kf: all white), diode 1.6e-6, BSIM4 7.0e-7 (flicker part 1.1e-9), BJT 4.7e-4 (its tau/T quasi-static error: 3e-4 for the white part; bound for the BJT 1e-3). The old
evaluation was 9.2e-2 / 3.1e-1 / 3.0e-1 away from the folded reference (diode / BJT / BSIM4: the unfolded time average overestimated by the variance of `a(t)`). Red proof
`PNOISE_FLICKER_OLD=1` fails the diode check.

**Effect on R7** (`.pnoise dev=1`, `r7_pnoise.cir` with 3 points, new / old):

| f | total new / old | OTA pair `m1..m4` flicker | `mls` / `mhs` flicker | gate-driver and comparator flicker |
|---|---|---|---|---|
| 1 kHz | 2.705e-7 / 4.224e-7 (0.64) | 0.87 | 0.80 / 0.51 | 0.01 - 0.08 |
| 10 kHz | 1.292e-7 / 5.784e-7 (0.22) | 0.33 | 0.80 / 0.52 | 0.02 - 0.11 |
| 100 kHz | 2.579e-8 / 5.826e-8 (0.44) | 0.30 | 0.83 / 0.59 | 0.03 - 0.13 |

The dominant contributors at 1 kHz are the OTA devices (flicker, 87 % of the power) whose bias is constant; for them the change is the sideband energy of the loop-shaped `H(t,f)` (PWM), not a modulated
source. The gate drivers and the comparator switch fast and their flicker is almost entirely removed, as the chopper tables predict. R7's values from #542/#544 for Rload/Resr thermal
are unchanged.

### 2. Absolute level against transient-noise ensembles

Method (`tools/pnoise_mc.py` chain, validated on an LTI RC to 0.2 dB per band): every noise generator is made `noisy=0` and replaced by a TRNOISE-driven current source with the same
one-sided density; the noisy circuit and a noiseless copy share every time step (two subcircuit instances), the difference is the noise response, boxcar-integrated, Hann/Welch
(nfft 8192, dt 200 ns, 610 Hz bins, 50 % overlap), 9-bin band averages compared with `.pnoise` at the same bins. Statistics stated before comparing: 1 sigma of a band = 0.065 dB for
48 seeds x 20 ms (960 segments), 0.08 dB for 32, 0.11 dB for 16. A common amplitude factor K lifts the output above the 2 nV/rtHz floor of the differencing chain (the engine value is scaled by K; the
loop is linear in it, checked by running two K).

**R7 itself is not covered by an ensemble.** BSIM4 has no transient noise, a time-varying BSIM4 `gm(t), gds(t)` density cannot be written as a B-source without the device's internal
quantities, and the 2161-step-per-period transistor-level deck costs hours per seed (R7: 20 ms x 32 seeds). R7 is covered by analogy with two decks that contain the same machinery:
a closed-loop PWM converter and a MOS switch + diode converter. What is checked on R7 is the closed-form part: the folded flicker (above), and the Rload/Resr thermal values (test_r7).

**(a) R3, closed-loop voltage-mode buck** (`tools/pnoise_mc_r3.py`, `pac_r3_ctrl2out.cir`: 500 kHz, type-III loop, tanh switches and comparator, all noise from nine resistors, K = 100,
TS = 10 ns, `.options reltol=1e-4` in the transient, 32 seeds x 20 ms, 1 sigma per band 0.08 dB): per band (4.3 .. 42.7 kHz) MC/engine in dB -0.30, +0.08, -0.17, -0.05, +0.17, +0.26, +0.14; worst 0.30 dB,
mean +0.02 dB. The same ensemble at the default transient `reltol=1e-3` was +0.14 (K 100, 48 seeds) / +0.28 (K 30) / +0.14 (K 10) dB mean with bands up to +0.65 / -0.75 dB: that offset is the
transient's own truncation error in the small noise response (it is gone at 1e-4), not the engine. The engine value is converged in the grid (Points 200 / 800 / 3200: 3.1013 / 3.1016 / 3.1020e-8 at 9.8 kHz).

**(b) MOS1 high-side switch + freewheel diode buck** (`tests/fixtures/pnoise_m1buck.cir`, `tools/pnoise_mc_mos.py`: 100 kHz, D 0.5, 100 ns edges, CCM, 4kT-noise of Rl, Resr, Rload, M1
channel thermal (`2/3 gm`, the default nlev 2) and D1 shot noise 2q|Id| injected as B-sources with the instantaneous density `sqrt(S(v(t)))`; M1 24 %, Rload 56 %, Rl 17 %, D1 3 % of the power at 1 kHz,
M1 6 % at 10 kHz; `reltol=1e-4`): K = 300, TS = 20 ns, 48 seeds x 20 ms, bands 4.3 .. 42.7 kHz: +0.21, -0.29, +0.02, -0.08, +0.06, +0.03, -0.04 dB; worst 0.29 dB, mean -0.01 dB (1 sigma 0.065 dB per band).
K = 100, 24 seeds: worst 0.24 dB, mean +0.04 dB (the circuit is linear in K). Per source (12-16 seeds, K 300): Rload alone -0.01 dB mean (worst 0.16), M1 alone -0.15 (worst 0.30; TS 10 ns: -0.12, worst 0.23),
D1 alone -0.13 (worst 0.37).

Things that went wrong on the way and are kept here so the next reference does not repeat them: (i) the first MC injected the nlev 3 MOS1 channel formula while the default model (nlev 2)
uses `2/3 gm` (4.9x higher in power, a pure scale factor that showed up as +7.8 dB flat in frequency; found by comparing the B-source with plain `.noise` on an LTI deck, not by the engine);
(iii) the first R3 ensemble used the default transient `reltol` and sat 0.1-0.3 dB high (see (a)); (iv) a discontinuous-conduction version of the converter (L 10 uH, 5 ohm load) rings at 3.5 MHz on a floating switching node after the diode turns off and gave +5.8 dB for the diode
alone; the CCM version above does not.

The residual +-0.3 dB per band, against a 0.065-0.08 dB statistical error, is the calibration accuracy of the Monte-Carlo chain stated in #509 (+-0.2 dB per band on a known LTI answer) plus the
band limit of the injected white noise (TS); the means over the bands are -0.01 / +0.02 dB.

### Tests (all in the working tree, not committed) and regression

Regression with `aaba50ba3f49-556-2` / ASan `-556-1-asan`, `OMP_NUM_THREADS=2`, logs `~/Work/data/ngspice_augmentation/556/reg/` against `port-s109/`: test_pnoise 9/9 (was 7/7: + the fold check and its red proof), test_pac 32/32,
test_pxf 15/15, test_r7 24/24, test_xspice 26/26, test_pac_charge 52/52 (was 51: + `red_proofs_556`; the LTI-bias numbers move by <= 1e-7), test_handbook 165/165 (second run after the example-38 change; the first run had 3 FAIL, all in example 38): example 38 compared the `.pnoise` of a flicker-bearing diode / BJT under a 100 kHz LO with
the plain quasi-static average (the white answer) and failed 3 checks until its reference was replaced by the folded one (`docs/examples/handbook/38_pnoise_opoint.py`, `tests/test_handbook.py`: diode 4.1e-7 (bound 1e-4), BJT
4.5e-4 (bound 1e-3), stale-state red proof x3.01 against the folded reference (was x4.31 against the white one; bound now > 2); the handbook HTML on the shelf quotes the old 4.47449e-09 / x4.31 and needs the new numbers 4.28631e-09 / x3.01). Everything else in the logs
is identical to s109 apart from wall times.

`tests/test_pnoise.py`: `check_pnoise_fold` (+ red proof `PNOISE_FLICKER_OLD=1`: every case must be off by more than 2 %, measured 3.4 % .. 10x), fixtures `pnoise_fold_chop_res/chop_lp/mos.cir`;
`PNOISE_MC=1` adds the two ensembles above (16 seeds each, bound 0.5 dB per band, 0.25 dB mean). `tests/test_pac_charge.py`: `check_pnoise_qs` folded reference, `red_proofs_556`.
`tools/pnoise_mc_mos.py`, `tools/pnoise_mc_r3.py`, `tests/fixtures/pnoise_m1buck.cir`.

### Statements of earlier sections this replaces

#544 item 1 ("Flicker under a periodic bias ... overestimates the folded flicker of a fast-switching circuit ... not validated") and its "Not done" line ("sideband folding of flicker: no independent reference"), and the #509 skip
list entry on cyclostationary flicker (docs "What was skipped"): the evaluation at the output frequency was wrong in both directions, is replaced, and is validated by three closed forms and four device-model
quasi-static decks above. The #544 `check_pnoise_qs` numbers "measured 1.5e-7..3e-4" were against the white reference and the engine's old evaluation; the table above is against the folded reference.

### Not done / limits

- R7 absolute level: no ensemble (see above); covered by analogy (R3 closed loop + MOS1/diode converter, mean offsets +0.02 / -0.01 dB).
- Flicker only: the engine has no other frequency-dependent density; a user model with its own `NevalSrc N_GAIN` flicker call is matched by slot count, a mismatch falls back with a warning.
- A flicker density that is not a power law in `f` (none among the device models here) is not folded (warning, old evaluation; the branch is untested). `1overfsw` diode sidewall slots are not folded.
- Tail beyond `Q` = 2048 sidebands is a `q^-2` model, not summed; the floor at `1e-3 fs` is a regularisation of an infinite density.
- Extra ensemble caveat: the injected white noise is band-limited by TS (sinc^4); the engine's is not.

## Stateful analog XSPICE models in td small-signal analyses, and `.pxf adjoint=1` / `.pnoise` on free-running orbits (card #557, 2026-10-03)

Branch `aug-557` on `aug-s111` (s111), commits `ec6baf66e` (free-running adjoint/noise), `4c023b825` (stateful code models). Builds `aaba50ba3f49-557-6`
(ASan twin `aaba50ba3f49-557-asan`). Outputs `~/Work/data/ngspice_augmentation/557/`. Tests: `tests/test_xspice.py` (7 new checks),
`tests/test_pxf.py` (2), `tests/test_pnoise.py` (see below); fixtures `tests/fixtures/xspice_r3_sxfer.cir`, `r5_freerun_plain.cir`; tool
`tools/pnoise_mc_r5.py`. Switches (test/experiments only): `PAC_HID_RED=1`, `PAC_HID_NO_JCOL=1` (red proofs), `PAC_HID_ALLOW_SLEW=1`, `PAC_HID_ALLOW_FLIP=1`
(lift the two refusals below), `PAC_HID_OFF=1` (back to the #551 refusal), `PAC_HID_TRACE=1`.

Two refusals lifted, one of them in part:

| refusal (#551 / #504 / #544) | now |
|---|---|
| td `.pac/.pstb/.pxf` through s_xfer, int, d_dt | **lifted** (hidden-state recursion, below) |
| td through slew, hyst | lifted only where the model stays on one branch (never slews / never flips); otherwise **refused with the reason, measured** |
| `.pxf adjoint=1` on a `trigger=` orbit | **lifted** (it is the exact transpose of the default forward operator) |
| `.pnoise` on a `trigger=` orbit | **lifted**, validated against a transient-noise ensemble |
| `.pnoise` / `.pxf adjoint=1` through stateful code models, `method=gear` with them, saltation with them | still refused |

### 1. td small signal through stateful analog code models

**Model.** A stateful analog code model keeps its state in `CKTstate` slots (`cm_analog_alloc`: the integrator chain of s_xfer, the integral of int, the old
input of d_dt, ...). One accepted step of the real transient is a map `(x_k, h_k) = S(x_{k-1}, h_{k-1})` of the node vector and of those slots.
Linearised at the accepted point:

    J'_k y_k + B_k dh_{k-1} = (forcing + capacitive history of the generic recursion of #542)        dh_k = P_k y_k + R_k dh_{k-1}
    B = dF/dh1,  P = dh0/dx,  R = dh0/dh1,   F = J x - rhs  (the load residual at the converged point, so a model's output current / branch equation)

The slots are extra Floquet states: virtual rows `neq+1 .. neq+nh` of the generic state vector (`Cprev/Cnew`, `srows`), so the boundary solve
`(e^{jwT} I - Phi) xi = z`, pass 2, `.pstb` and `.pxf adjoint=0` carry them without change. B, P, R and the model's input columns of J are **central
differences of the model itself** (`pac_hid_stage`: full `CKTload()` at x_k with one input node, or one `CKTstate1` slot, moved; the step is enlarged until
the response is resolved to 1e-7 of its entry, at most three times, so a slot whose natural scale is 1e5 times its value on the orbit is still seen): no
per-model code. The first prepass samples the dependence pattern of every candidate slot; **active** = reachable from the inputs (P or an R chain) AND read by
a later step (B or an R chain). s_xfer 2 of 20 slots, int 1 of 4, d_dt 1 of 2; coefficients copied state1 -> state0, time stamps and write-only
copies (`in`, `out`) are dropped. J's model input columns are replaced by the finite-difference dF/dx before the factorisation: **s_xfer's declared
`PARTIAL(out,in)` is the value of an integrator, not a derivative** (it is only a Newton iteration matrix); the charge matrix C_k = (J(ag0)-J(ag0/2))/(ag0/2)
is taken with `g_xs_ag0_fix` so the model integral does not enter it (see 3. below). Between recording runs the converged orbit's carried XSPICE state
(`XSsnap`, #551) is put back (`PSSxsPacBegin`), so each recording period starts where the orbit's last one ended.

**What the linearisation has to reproduce.** `cm_analog_integrate` is first order: at trapezoidal order 1 and 2 it is backward Euler, `I_k = I_{k-1} + h u_k`
(order 2 halves `geq` and drops the history, `cm.c`: "this code needs to be redone"). The transient of a model chain is therefore the backward-Euler
discretisation of its transfer function, and td small signal, which linearises the transient, must equal **that**, not plain `.ac`; it approaches `.ac`
as O(h). Reference chain, uniform step h (exact z-domain transfer of the code, `tests/test_xspice.py td_disc`): int `g h /(1-z^-1)` (closed loop `g h/(1-z^-1+g h)`),
d_dt `g (1-z^-1)/h`, s_xfer `I1 (1-z^-1) = h v2, I0 (1-z^-1) = h I1, v2 = G u - d0 z^-1 I0 - d1 z^-1 I1`. `.ac` itself of these models was checked against the
continuous hand calculation first: 5e-16 (s_xfer), 2e-16 (int), 7e-16 (d_dt).

**Results**, periodic LTI orbit (sine source 0.3 V + 0.2 V at 1 kHz on the model input, so the hidden states are not constant), `pacpoints` 1e5 (h = 10 ns), 1 kHz..100 kHz,
worst of 11 points, `tests/test_xspice.py check_td_stateful_lti`:

| model | td vs exact BE-discrete | td vs `.ac`, pacpoints 1e5 / 2e5 / 4e5 | Richardson (2 td(h/2) - td(h)) vs `.ac`, 1e5/2e5 | hidden states |
|---|---|---|---|---|
| s_xfer (1/(1+2.2e-6 s+2.5e-12 s^2)) | 1.7e-7 | 9.4e-3 / 4.7e-3 / 2.4e-3 | 2.1e-5 | 2 |
| int in a closed loop (summer, gain 1e5) | 5.7e-8 | 3.1e-3 / 1.6e-3 / 7.8e-4 | 1.7e-6 | 1 |
| d_dt (gain 1e-6) | 1.8e-7 | 1.3e-2 (2.5e4) / 6.3e-3 (5e4) / 1.6e-3 (2e5) / 7.9e-4 (4e5) | 3.3e-6 | 1 |
| slew (rise/fall 1e9 V/s, never slews) / hyst (stays on one branch) | 0 / 0 | 0 / 0 | - | 0 (J exact) |

(the constant-state decks of #551's clock layout give the same: 3e-6 against the discrete reference at 25 000 points). The error is first order in h and the
extrapolated value is `.ac` to 2e-5, as it must be. The scaling of the closed loop with the number of states was checked on an s_xfer loop (C(s) = 1e5 (1+6.4e-5 s)/(s (1+s/1.26e6)), two
hidden states, `.ac` cannot run: no DC solution of an integrating s_xfer): 2e-3 .. 4e-3 against the hand calculation at pacpoints 2e5 (O(h) again).

**Closed periodic loop**, `tests/fixtures/xspice_r3_sxfer.cir`: R3 with its type-III op-amp compensator replaced by an s_xfer, C(s) = 1e5 (1+s/wz)^2/(s (1+s/wp)),
wz = 2 pi 5 kHz, wp = 2 pi 200 kHz; regulates 3.3 V at 500 kHz; 3 state rows + **2 hidden states**; control-to-output (series injection Vctrl) vs a frequency-response
measurement of the real nonlinear circuit (single-bin DFT of the response to a 0.01 V sine, reltol 1e-5, maxstep 1 ns, last 2 ms = whole numbers of both periods;
amplitude 0.0025 repeats it to <= 0.0004 dB / 0.004 deg):

| f | FRA | td `pacpoints=2000` | delta | `pacpoints` 4000 | 8000 | default grid (1727 steps) |
|---|---|---|---|---|---|---|
| 4 kHz | -6.1833 dB / 12.682 deg | -6.1830 / 12.682 | +0.0003 dB / -0.001 deg | +0.0012 / -0.012 | +0.0011 / -0.019 | -0.008 / +0.215 |
| 20 kHz | -11.7160 / -95.560 | -11.7164 / -95.559 | -0.0004 / +0.001 | -0.0027 / +0.075 | -0.0074 / +0.022 | +0.084 / -0.196 |
| 50 kHz | -25.2293 / -131.765 | -25.2292 / -131.763 | +0.0001 / +0.002 | +0.0111 / +0.063 | +0.0001 / +0.030 | +0.009 / -0.288 |

(`tools/pac_perturbed.py` refused all three points of this deck, "UNRELIABLE": amplitude-halving 0.23 / 0.48 dB and step-refinement 0.5 dB: the models' own
first-order integrals make the nonlinear reference grid dependent; its raw values, 20 kHz -11.82 dB / -95.1 deg and 50 kHz -25.07 / -131.9, lie within those spreads of the table.)
Red proofs (`check_td_stateful_red`, `check_td_r3_sxfer`): `PAC_HID_RED=1` (the hidden states do not act on the nodes) the s_xfer deck is wrong by 1.0 (relative) and
the R3 curve by 25.6 dB; `PAC_HID_NO_JCOL=1` (declared PARTIALs) 6e-3; both must exceed their bounds.

**Four things the periodic orbit found that a constant-state deck hides** (all were error sources of 1e-2 ... 40x before they were fixed; the first R3 curve was 0.012 where the FRA says 0.49):
1. `cm_static_integrate` at trapezoidal order 2 uses `CKTag[0]` also for the history term `cur = -0.5 ag0 I_old`; the C_k extraction varies CKTag[0], so the model's
   integrals (and s_xfer's value-of-integrator "partial") moved between the two loads and became a bogus charge matrix entry, **proportional to the orbit amplitude**
   (sine amplitude 0.2 / 0.02 / 0.002 gave 1.3e-2 / 1.3e-3 / 1.3e-4). `cm.c` now uses `g_xs_ag0_fix` for every use of ag0 (set by `pac_jload`).
2. the model's declared PARTIAL (s_xfer) is not dF/dx (above): finite-difference columns.
3. `ckt->CKTorder` at the recording hook is already the next step's (as `PSSpacRecordStep` notes); the loads set the order recovered from ag (a precaution: it did not change a number on these decks).
4. the hidden state must be carried between recording runs, and a slot whose natural scale (1e5 V/s) is invisible on a DC orbit needs an adaptive finite-difference step.

**Refusals kept (each checked, `check_td_stateful_refusals_and_pxf`)**: `.pnoise` through any stateful model (`pac_noise_check_devices`, and the adjoint sweep does not carry the
hidden states); `.pxf adjoint=1` -> "adjoint sweep does not transpose its hidden states", falls back to the per-source engine, **bit-identical to `.pac`** (adjoint=0 too);
`.pstb` runs; `method=gear` with them (the models' integral is read at gear order p from `CKTstates[2..]`, not carried); saltation events with them; oneshot, the clock sources
(sine, square, triangle), capacitoric / inductoric and every digital / hybrid model as before.

**slew and hyst, measured, not guessed.** A slew model that slews part of the period, or a hyst that flips, is a regime change: the small-signal response is not defined.
`PAC_HID_ALLOW_SLEW=1 / ALLOW_FLIP=1` with a 3 kHz perturbation on a 1 kHz orbit, nonlinear transient reference, amplitude 0.01 / 0.0025: slew (rate 3000 V/s, 0.5 V sine,
slews near the zero crossings) |H| = 0.59 / 1.08, hyst (in 0.5 +- 1.5 V, flips every half period) |H| = 17.9 / 71.1 (a relay: gain ~ 1/amplitude), td (saltation missing) 1.89 / 0.217.
No reference exists to validate against, so both are refused when the prepass sees slewing (an active slew slot) or a state flip (`hyst`'s Boolean slot differing from the previous step's).

**Cost.** One accepted step loads the circuit 2 (n_active + n_inputs) + 3 times more than the generic recursion (the prepass samples the dependence pattern of all candidate slots at every step, or every k-th step beyond 20 000 steps per
period). R3 + s_xfer (24 candidate slots): prepass 0.044 s at the default grid (1727 steps), 0.32 s at `pacpoints=2000` (3524 steps); PAC 0.006 s per frequency; `.pstb` 24 points 0.3 s.

### 2. `.pxf adjoint=1` and `.pnoise` on a free-running (`trigger=`) orbit

**What the forward operator is.** The free-running forward `.pac/.pstb/.pxf` is, since #517, the fixed-period LPTV engine on the polished orbit
(`(e^{jwT} I - Phi) xi = z`): the period-variation projection of #489 (`Proj = I - f(T) e_trig^T / e_trig.f(T)`, bordered with the phase condition) is opt-in
(`PAC_FREERUN_PROJ=1`, "not needed for w != k ws") and is not formed on a generic state set. For w != k ws the phase mode (the unit Floquet multiplier) is not excited
singularly, and the fixed-period solution is the exact response of the linearised orbit; that is what was validated against `pac_perturbed --freerun` in #503/#544.
So the adjoint does not need a bordered system: it is the transpose of **that** operator, and the existing sweep (`pac_engine_run_adjoint`, `_adjoint_noise`) already is.
The refusals were `if (trig_idx >= 0)` guards. They are now taken only with the opt-in projection (`PAC_FREERUN_PROJ`, adjoint refuses, `.pxf` falls back to per-source, `.pnoise` refuses).
Three more places had to change for `.pxf`/`.pnoise` to run on a trigger orbit at all (neither had ever been run there): `inp2dot.c` `dot_pxf`/`dot_pnoise` did not parse `trigger= trigval= trigdir=`
("unknown trailing token ... ignored", then "doAnalyses: operation not supported"), and `PSSshootFreerun` only post-processed `.pac/.pstb` cards (`PSSdoPAC || PSSdoPSTB`).

**`.pxf adjoint=1` = forward per source** (`tests/test_pxf.py check_pxf_adjoint_freerun_r5/_r8`; forward under `PAC_EXACT_JACOBIAN=1`, the fair comparison of #504):

| deck | sources | frequencies | adjoint=1 vs adjoint=0 | adjoint=0 vs `.pac` |
|---|---|---|---|---|
| R5 (`r5_freerun_plain.cir`, trigger q, native path) | Vin, Vrefsrc | 20 / 55 / 90 kHz | 2.4e-9 / 2.2e-9 / 2.0e-9 (Vin), 1.2e-10 (Vrefsrc) | 8e-6 (default stale-Jacobian forward, the #504 O(1/Points) floor) |
| R8 (BSIM4 hysteretic buck, trigger cx, generic path, hidden-state continuation) | Vin, Vrefsrc | 20 / 50 / 80 kHz | 3e-14 (worst) | 0 (bit-identical) |

Red proof: the forward engine with the projection (`PAC_FREERUN_PROJ=1`) is a different operator and moves the answer by a relative 3.5 (R5 Vin curve), and with it the
adjoint refuses ("PXF-TD: adjoint engine REFUSED ... PAC_FREERUN_PROJ ... using the per-source engine") instead of transposing the wrong operator.

**What `.pnoise` reports on a free-running orbit.** The exact output noise PSD of the linearised cyclostationary orbit at the (baseband, sideband-0) frequency f: every device noise
source folded through all sidebands (white noise exactly by Parseval, flicker by #556), **including the Floquet phase mode** (the unit multiplier: the response of the orbit's phase to the
noise). Its sideband-0 projection is the mean of x*'(t), which is zero, so the phase-diffusion term of the baseband density does not diverge as f -> 0: R5, 10 Hz .. 1.8 kHz, 1.383e-10 -> 1.357e-10 V/rtHz
(finite, flat within 2 %). The phase diffusion shows up where it physically does, as the 1/(f - k fsw)^2 skirt of the carrier harmonics (k >= 1), outside the validated band: the linearised
result diverges at f -> k fsw and means nothing inside the Lorentzian line width. The existing band guard of #503 keeps those frequencies out: a sweep whose top frequency exceeds 0.25 fsw
is REFUSED (`PSS_FREERUN_PAC_EXPERIMENTAL=1` runs it, unvalidated). So nothing is excluded or refused in band, and nothing near a carrier is claimed.

**`.pnoise` vs a transient-noise ensemble, R5** (`tools/pnoise_mc_r5.py`, `tests/test_pnoise.py check_pnoise_freerun_r5`). Method of #556 on the autonomous converter: Rdcr, Resr, Rload made
`noisy=0` and shunted by TRNOISE current sources of one-sided PSD 4kT/R x K^2 (K = 50 keeps the output noise 5e-9 V/rtHz and the loop linear; K = 100 leaves the linear regime after ~3 ms),
two copies of the converter in one transient (noisy, noiseless), boxcar-integrated difference, Hann/Welch 8192 x 200 ns (610 Hz bins), 9-bin bands at 15.3 / 40.3 / 84.8 kHz (0.036 / 0.095 / 0.20 fsw).
Both copies are autonomous: their switching phases decorrelate under the noise (the phase diffusion), which is exactly what the ensemble has to contain and `.pnoise` has to predict.

Two things the free-running case needed that R3's ensemble did not:
1. **A numerical floor, K-independent.** At `reltol` 1e-4 (R3's setting) the two copies of the autonomous converter differ by 90 .. 170 microvolt rms even at K = 0 (two copies solved in one matrix converge
   to the Newton tolerance independently, the free-running phase has no restoring force). That white-ish floor read 3.9e-9 / 5.9e-9 / 14.5e-9 V/rtHz at K = 15 / 50 / 150 (`.pnoise` x K: 1.6e-9 / 5.4e-9 / 16e-9)
   and made the first ensemble 0.8 .. 2.4 dB "above" the engine. `MC_RELTOL=1e-6` removes it (the tool's default stays 1e-4, R3 needs nothing more).
2. **The injected noise bandwidth matters.** TRNOISE of step ts is white only to 1/(2 ts); `.pnoise` is white to infinity and R5's comparator folds the noise of many harmonics of fsw into the
   band. The ensemble approaches the engine from below as ts shrinks (K = 50, `MC_RELTOL=1e-6`, pooled over the three bands; the table is MC / `.pnoise` x K in dB per band):

| ts | seeds x ms | 15.3 kHz | 40.3 kHz | 84.8 kHz | pooled (1-sigma) |
|---|---|---|---|---|---|
| 40 ns | 6 x 8 | -5.54 | -5.13 | -5.16 | -5.27 (0.19) |
| 20 ns | 6 x 8 | -3.23 | -2.40 | -2.91 | -2.85 (0.19) |
| 10 ns | 6 x 8 | -1.78 | -0.78 | -0.79 | -1.12 (0.19) |
| 5 ns | 16 x 20 | -0.74 | -0.25 | -0.19 | -0.39 (0.064) |
| 2.5 ns | 12 x 16 | -0.39 | +0.30 | +0.37 | +0.09 (0.082) |

(the power deficit is linear in ts, 0.70 / 0.48 / 0.23 / 0.07 of the power missing at 40 / 20 / 10 / 5 ns; per band 1-sigma 0.33 (6 seeds) / 0.11 (16 x 20) / 0.14 dB.) At the finest step the engine is
within +-0.4 dB per band (2.7 sigma at worst, 15 kHz) and +0.09 +- 0.08 dB pooled: agreement at the 0.1-0.4 dB level, short of #556's 0.07 dB because the finite-bandwidth deficit has to be driven
down by a step that makes each seed 6.4 M steps (about 3.5 min of wall time per 16 ms seed at 2.5 ns, 3 at a time; a 0.07 dB ensemble was not run). The test stores the ts = 5 and 2.5 ns values, gates the 2.5 ns bands at 0.6 dB and the pooled
mean at 0.3 dB, requires the 5 ns ensemble to lie below the 2.5 ns one, and fails a +2 dB engine (red proof); `PNOISE_MC=1` re-runs the 5 ns ensemble live.

R8 (BSIM4, trigger cx): `.pnoise` runs (20 / 50 / 80 kHz: 1.160e-7 / 1.082e-7 / 1.080e-7 V/rtHz, generic path, 65 s PSS); **no independent reference** (no BSIM4 ensemble exists for a free-running orbit).

### Statements of earlier sections this replaces

* #551 "Model policy" and "Not done": "s_xfer int d_dt slew hyst ... REFUSED in td `.pac/.pstb/.pxf/.pnoise`" and "td small-signal analysis through XSPICE state ... refused, not implemented. A stateful-analog recursion would need the
  model's per-step map (the integral is not a node charge; `J(ag0)` is `A + B/ag0`)" -> implemented as the finite-difference per-step map above; still refused: `.pnoise`, gear, saltation, oneshot, clock sources, capacitoric/inductoric, digital.
* #504 "Scope cut" and #509: "the free-running (`trigger=`) period-variation projection are NOT transposed by this engine ... `trig_idx>=0` -> PAC_ADJOINT_REFUSED" / "`.pnoise`/`.pxf` on a free-running orbit refused" and #544 §3 "`.pnoise` and
  `.pxf adjoint=1` on a free-running orbit stay refused" and its "Not done" line -> only the opt-in projection refuses.
* #503 FOURTH pass / Known limits "`.pac`/`.pstb` on a free-running orbit run by default <= 0.25 fsw": the same band now applies to `.pxf` and `.pnoise`.

### Not done / limits (2026-10-03)

* Stateful models: `.pnoise` (no transposed hidden-state sweep; the model's noise is in any case outside the engine's device scope), gear / order > 2 (`CKTstates[2..]` history), saltation with them, a slew model that
  slews and a hyst that flips (no small-signal response exists, measured), oneshot, sine/square/triangle clock models, capacitoric/inductoric, every digital/hybrid model; free-running + XSPICE still refused (#551);
  untested: models in a Krylov / matrix-free `.pss`, many instances (the cost is 2 (n_active + n_inputs) + 3 loads per step, all steps in the prepass).
* td vs plain `.ac` is O(h) for the integrating models (their own backward-Euler integral), not 1e-6: the exact reference is the discrete transfer; at the default recording grid of an R3-class deck (1727 steps per 2 us) the R3 + s_xfer curve is
  0.08 dB / 0.29 deg from the FRA, 0.0004 dB / 0.002 deg at `pacpoints=2000`: use `pacpoints` for a stateful analog model in a loop.
* `.pstb` through the models runs and satisfies the closed-loop identity H_cl = G/(1+T) to a constant factor (spread 3.5e-4 over 4 .. 50 kHz, phase drift 1 deg; the constant 1.0926 is the operating-point difference of the opened-loop G), it was not compared with an FRA loop-gain measurement.
* The hidden analog state of the PSS orbit (#551: lagged block iteration) is only as periodic as that iteration made it: the R3 + s_xfer orbit's "analog node drift" was 8e-4 relative at 12 latching baselines (|rho| 0.87 of the slowest mode), the
  td curve is unaffected at the 0.0004 dB level but a stateful loop with a slower hidden mode needs a longer `stabtime`.
* Free-running: the adjoint of the opt-in projection (`PAC_FREERUN_PROJ`), `.pstb` on R8 (no probe), the R8 `.pnoise` absolute level (no reference), anything above 0.25 fsw, `.pnoise` on R5 at 0.07 dB (see above).
* `tools/pac_perturbed.py` is not a usable reference for a deck with first-order XSPICE integrals (its self-checks refuse; spread 0.2-0.5 dB); the FRA (one sine, tight reltol) was used.

### Regression (OMP 2, one suite at a time; logs `~/Work/data/ngspice_augmentation/557/port-557/`, baselines `port-s111/`)

| suite | s111 (`port-s111`) | 557 | note |
|---|---|---|---|
| test_xspice.py | 26/26 | 32/32 (build 7) | + 7 checks (LTI stateful, red proofs, R3 + s_xfer vs FRA, `.pstb` identity, ASan, refusals + `.pxf`), refusal list updated; R9 included |
| test_pxf.py | 15/15 (3 info) | 18/18 (3 info) (build 7) | + R5 and R8 free-running adjoint == forward (+ red proof / projection refusal), ASan on the free-running adjoint |
| test_pnoise.py | 9/9 (1 info) | 11/11 (2 info) (build 7) | + R5 free-running vs ensemble, refusals (projection, band), optional live ensemble |
| test_pac.py | 32/32 (2 info) | 32/32 (2 info) (build 6) | |
| test_pac_charge.py | 52/52 | 52/52 (build 6) | ASan chain clean |
| test_r7.py | 24/24 | 24/24 (build 6) | |
| test_pss.py | 44/48 (4 info) | 44/48 (4 info) (build 6) | |
| test_handbook.py | 165 | 165 (43.log regenerated on build 7, the rest from the build-6 run) | the first run had the one intended failure ("refused: .pac on s_xfer": the example now shows `.pac` running through s_xfer and refuses oneshot instead: `docs/examples/handbook/43_xspice.py`, `tests/test_handbook.py`) |

Build 6 and build 7 differ only by one stderr line ("PXF-TD: adjoint engine REFUSED ... PAC_FREERUN_PROJ", commit `ca28e3277`); test_xspice/test_pxf/test_pnoise were (re)run on build 7, the others on build 6.
ASan (`aaba50ba3f49-557-asan`): the new paths run clean: R3 + s_xfer `.pac` / `.pstb`, LTI s_xfer `.pxf adjoint=1` fallback (`check_asan_stateful`), R5 `.pxf adjoint=1` (`check_asan_pxf_freerun`), R5 `.pnoise`, R8 `.pnoise` (manual runs), test_pac_charge's chain.

Not run: test_pss_krylov, test_r3d_r4d, test_tran_from_pss, test_loss (the change is in `pacshoot.c` td paths, the free-running post-processing gate and `.pxf/.pnoise` card parsing; the PSS engines are untouched).

## Saltation at charged switch nodes and non-quasi-static MOS channels in td PAC (card #559, 2026-10-03)

Patch 0117 (branch `aug-559`, commit c880d9c21 on the s115 tip). Two refusals of #542/#544 looked at; one lifted
inside a validated envelope, one lifted for `.pac/.pstb/.pxf` and kept for `.pnoise`.

### 1. Ideal-switch (sw) events where the switch node carries device-model charge

**Jump on the charge vector** (`pac_salt_qjump_prepare/finish`, pacshoot.c). For a switch power node whose loaded charge row is not the native
capacitor/inductor one (prepass `pac_salq_note`: diode cj/tt, MOS junction/overlap, ...):

    dQ(row)   = coeff * (-dI at pos, +dI at neg),   coeff = (y(posCntrl) - y(negCntrl)) / h_dot,   dI = dG v_sw
    dy        = C_ss^-1 dQ over the charge-coupled block of state rows (capacitive divider to e.g. the gate; singular block: dQ / diag C)
    Q_k       = Q_smooth + dQ  (exactly),   y_k = y_smooth + dy,   T_k = ag0 (Q_smooth - Q_{k-1}) - a1 T_{k-1}   (T from the SMOOTH part)

Nothing is divided by a capacitance (the native jump divided by the native local C). T carries no jump spike, so the trapezoidal recursion does not
re-inject it into the next step and no order-1 settle window is needed for these nodes (`PSS_SALT_QSTEPS=N` forces N steps, test only: N = 0..2
agree to 1e-3 dB; the first scheme tried, a tau = C/G window, moved the answer by 3 dB with the window length: the freewheel diode's conductance is
not on at the accepted point of the event step, and the diode's stored charge is not a capacitor behind a conductance).

**Validated envelope** (guard `pac_gen_salt_guard`): the node has the model charge ON TOP OF a native capacitor. (Card #575 narrowed this: the native capacitor must dominate the node or sit with an RC snubber; see the section of that card.) Evidence, `.pac` ctrl2out vs
tools/pac_perturbed (perturbed periodic steady state, `.pss` backend, halving <= 0.004 dB, step refinement <= 1e-4 dB), 20/50/100 kHz:

| deck | worst | 
|---|---|
| `pac_r6od_ctrl2out.cir`: ideal sw + freewheel diode cj 300 pF tt 2 ns + 1 nF at the node | 0.017 dB / 0.16 deg |
| `pac_r6odm_ctrl2out.cir`: body diode (is, cbd) + cgdo of an off NMOS, gate behind 1 kOhm moves with the node | 0.021 dB / 0.15 deg |

Result is independent of reltol (1e-4 .. 1e-6: -7.4797 / -7.4795 dB at 20 kHz, reference -7.4671). A converged transient FRA
(`fra_curves_fork.run_transfer_point_conv`, 3 ms settle, 0.1 V, tmax 2 ns) gives -7.415 / -23.775 dB at 20 / 50 kHz, i.e. 0.05 dB / 0.3 dB from the
perturbed PSS (the FRA's own spread; a 0.05 V run 0.9 dB, 100 kHz noise-limited): the FRA cannot referee below ~0.3 dB on an ideal-switch buck, the
perturbed PSS (self-checking, refuses unreliable points) is the reference. Red proofs: `PAC_SALT_NATIVEJUMP=1` (the old jump without the refusal) 0.61 dB
and 4.0 dB off (4.8 / 0.32 dB before card #568 shortened the settle window); `PAC_SALT_NOQ=1` restores the refusal. The charged decks do not use the order-1 settle window, so card #568 (next section) leaves them where they were.

**Kept refused, and exactly why**
- A switch node whose only charge is device-model charge (no native capacitor at the node; `PAC_SALT_QONLY=1` runs it): R6od with the node
  capacitor removed is 0.06-0.10 dB off and drifts with reltol (0.073 dB at 1e-4, 0.099 at 1e-6); the diode-only inductor-adjacent
  `pac_sw_ctrl_direct` + cj 50 pF/tt 1 ns is 0.30 dB off at 20/100 kHz and 3 dB off at reltol 1e-6 (the node is ~12 pF reverse biased: the
  post-event transient lasts ~100 ps, about one step, so the impulse is not resolved), while the quasi-static inductor-routed jump (#510) of the
  same deck without charge is robust (-1.237 .. -1.213 vs -1.19). A node whose charge vanishes at the event (tt without cj, reverse biased) gives a
  nonsense curve (-42 dB) without the check.
- Mutual inductance on the inductor-routed branch (unchanged).
- Adjoint/noise with events (`.pxf` falls back to the forward engine, `.pnoise` refuses), `method=gear` with events (unchanged).

### 2. Non-quasi-static MOS channel (BSIM4 trnqsmod, BSIM3 nqsmod)

The relaxation "charge-deficit" node qNode of the transient NQS model is an ordinary circuit node: its charge terms are proportional to ag0 (the
third-load affinity guard passes) and the rest is conductive, so the generic recursion carries it as a state (no hidden-state machinery as in #557
is needed; the #542 refusal rested on a measurement against `.ac`, which is the wrong reference). `.ac` is NOT a reference: acnqsmod is a single-pole
approximation (gm, gds, C scaled by 1/(1+j w tau)); BSIM4 at 10 MHz: `.ac` -1.05+3.37j, transient model -1.64+2.34j, quasi-static -2.56+2.95j.
Reference: converged transient FRA of the same model (LTI bias, L = 20 um, tau ~ 0.1 us; two settings differ 3e-5..8e-5).

| deck | td PAC vs FRA (relative) |
|---|---|
| BSIM4 trnqsmod=1, 1 / 3 / 10 / 30 MHz | 1.4e-5 / 1.8e-5 / 1.8e-5 / 2.0e-4 |
| BSIM3v3.3 nqsmod=1, 3 / 10 / 30 MHz | 1.8e-5 / 1.2e-5 / 4.8e-5 |
| BSIM4 acnqsmod=1 only (transient quasi-static), 10 MHz | 1.8e-5 (td linearises the transient model; acnqsmod affects `.ac` only: note printed) |
| BSIM4 NQS stage, 1 MHz LO 0.35 V, 3 / 10 MHz | 0.009 dB / 0.02 deg (FRA spread 0.01 dB) |

`.pxf adjoint=1` == forward to 5e-13 on the NQS deck. tools/pac_perturbed (.pss backend) is 0.011 / 0.034 dB and 0.3 / 0.5 deg from the FRA on the
LO deck: its period restart drops the relaxation node's derivative history, which the FRA does not; the FRA is the reference. Red proofs: the same
decks with trnqsmod=0 are 39 % (LTI) and 3.3 dB (LO) from the NQS references; `PAC_NQS_REFUSE=1` restores the refusal. Still refused: `.pnoise` on
any NQS model (noise of the relaxation node not validated), NQS-capable devices of other families (`PAC_ALLOW_NQS=1` overrides).

### Not done
- NQS device in the control path of an ideal-switch event (the combination was not run separately); NQS with a switching converter at power level.
- Q-mode saltation with `.pnoise` / adjoint; model-charge-only switch nodes (see above): a proper treatment of the impulse across an unresolved
  fast transient (apply the event at t_e inside the step and propagate the relaxation analytically) was not attempted.
- The native path's settle window (chargeless diode at the switch node, 5 dB off) was found here and is fixed under card #568, next section.

## Native saltation: the order-1 settle window is bounded in steps (card #568, 2026-10-06)

**The defect (present up to s115, no message).** After an ideal-switch event the S/T recursion (`.pss` sensitivities) and the C/T recursion (td
`.pac/.pstb/.pxf`) are forced to order 1 for a settle window, so that the jump written onto a stiff node is damped instead of ringing under the
trapezoidal rule (card #494). The window was 30 time constants, tau = C/G of the switch node, with G counted from the resistors and switches at
the node only (`pss_salt_local_cond`). When the conductance that really discharges the node belongs to another device, tau is far too long. On
`pac_r6od_ctrl2out.cir` with cj = tt = 0 (ideal high-side `sw`, a chargeless freewheel diode, 1 nF and a 10 ohm snubber at the switch node) G was the
snubber's 0.1 S, the window 300 ns of a 2 us period, 10..92 accepted steps per event. Order 1 on the recursion while the trajectory itself was
integrated at order 2 drops the derivative history of the slow states, and every forced step costs accuracy:

| forced order-1 steps after the event step | 0 | 1 | 2 | 3 | 4 | 6 | 10 | s115 (10..92) |
|---|---|---|---|---|---|---|---|---|
| error vs perturbed PSS, 20 kHz, dB | 9.4 | 0.014 | 0.021 | 0.032 | 0.052 | 0.18 | 1.3 | 5.36 |

**The fix (patch 0118).** The window also ends after 2 accepted steps (`PSS_SALT_SETTLE_MAXSTEPS`; `PSS_SALT_MAXSTEPS=n` overrides, `-1` restores
the time-only rule). One backward-Euler step already crushes a spike on a node that is stiff against the step, and a node that is not stiff does not
ring at order 2; zero steps is not enough (9 dB), which the table shows.

| deck | s115 | now | reference (tools/pac_perturbed, `.pss` backend, halving <= 0.004 dB) |
|---|---|---|---|
| R6od, cj = tt = 0, 20 kHz | -12.820 dB / -167.07 deg | -7.4842 / -164.938 | -7.4629 / -164.774 |
| 50 kHz | -28.757 / -160.24 | -23.4860 / -159.695 | -23.4604 / -159.684 |
| 100 kHz | -39.760 / -146.74 | -34.5507 / -146.679 | -34.5250 / -146.677 |

Independent of reltol (1e-3, 1e-4: -7.4842; 1e-6: -7.4788 dB at 20 kHz). Decks whose node stays stiff through the old window barely move: R6o 20 /
100 kHz -7.8052 / -34.9231 -> -7.8052 / -34.9238 dB; R6ob -9.5925 / -40.8084 -> -9.5913 / -40.8071; R6om -9.3284 / -37.6747 -> -9.3273 / -37.6738
(<= 0.0013 dB; they are no longer bit-identical to s115). The charged decks of card #559 (`pac_r6od`, `pac_r6odm` as shipped) do not use the window
and are unchanged (-7.4797 / -34.5434 and -7.5029 / -34.5942 dB). `pac_sw_ctrl_direct` (inductor-routed jump) already had a one-step window.

Tests: `test_pac_charge.py` `check_salt_native_window` (0.1 dB / 0.5 deg at 20/50/100 kHz; red proof `PSS_SALT_MAXSTEPS=-1` must exceed 3 dB) and
`check_salt_native_r6o` (R6o within 2e-4 of its s115 value).

### Not done
- The remaining 0.02 dB / 0.16 deg to the perturbed PSS is the same on the charged decks (which force no order-1 step at all), so it is not the
  window; it was not traced further (the reference's own amplitude is 0.05 V on a 2 V ramp).
- The time estimate itself still ignores diodes and transistors; with the step bound it can only shorten the window, never lengthen it.
- A node that needs more than two order-1 steps (stiff, with steps much shorter than its tau right after the event) was not constructed.

## Ideal-switch events at an inductor node without a dominant capacitor (card #575, 2026-10-06)

Patch 0119, build `aaba50ba3f49-s118`. Card #559 left a switch node whose only charge is a device's refused. Building the treatment for it on a
matrix of 20 asynchronous-buck decks, each with its own `tools/pac_perturbed.py` reference at 20/50/100 kHz (same settings as card #559; data
`~/Work/data/ngspice_augmentation/575/m`, `matrix-*.log`), also showed that s117 was silently wrong on six of them. Decks: SWC = the
`pac_sw_ctrl_direct` topology (ideal `sw`, vh = 0, diode, 10 uH, no capacitor at the node); OD = `pac_r6od_ctrl2out.cir` (vh = 5 mV, diode cj 300 pF
tt 2 ns unless said, optional node capacitor, 10 ohm / 1 nF snubber); ODM = `pac_r6odm_ctrl2out.cir` (MOS body diode + overlap charge).

| deck | treatment on s118 | worst dB / deg, reltol 1e-4 | reltol 1e-6 | s117 (worst dB, 1e-4) |
|---|---|---|---|---|
| SWC, diode cj 50 pF, tt 1 ns | fast node | 0.029 / 0.07 | 0.028 / 0.07 | refused |
| SWC, cj 50 pF | fast node | 0.002 / 0.07 | 0.002 / 0.07 | refused |
| SWC, cj 300 pF, tt 5 ns | fast node | 0.032 / 0.15 | 0.034 / 0.15 | refused |
| SWC, cj 1 nF | fast node | 0.008 / 0.34 | 0.008 / 0.34 | refused |
| OD, no node capacitor, no snubber | fast node | 0.020 / 0.15 | 0.020 / 0.14 | refused |
| ODM (MOS body diode), no node capacitor, no snubber | fast node | 0.005 / 0.14 | 0.002 / 0.14 | refused |
| OD, 10 pF native, no snubber | fast node | 0.019 / 0.15 | 0.019 / 0.15 | 0.429 dB, no message |
| ODM, 100 pF native, no snubber | fast node | 0.003 / 0.25 | 0.003 / 0.25 | 0.177 dB, no message |
| OD, 1 nF native, no snubber | charge jump | 0.033 / 0.10 | 0.032 / 0.10 | 0.033 |
| OD, 10 pF native + snubber | charge jump | 0.055 / 0.11 | 0.036 / 0.11 | 0.055 |
| OD, 100 pF native + snubber | charge jump | 0.040 / 0.12 | 0.032 / 0.11 | 0.040 |
| ODM, 100 pF native + snubber | charge jump | 0.041 / 0.11 | 0.022 / 0.11 | 0.041 |
| OD, chargeless diode, no capacitor, no snubber | chargeless, routed | 0.001 / 0.03 | 0.001 / 0.03 | 9.54 dB, no message |
| OD, chargeless diode, 10 pF native + snubber | native capacitor jump | 0.024 / 0.12 | 0.015 / 0.12 | 0.421 dB, no message |
| OD, chargeless diode, 100 pF native + snubber | native capacitor jump | 0.044 / 0.13 | 0.027 / 0.13 | 0.249 dB, no message |
| OD, chargeless diode, snubber only | chargeless, routed + snubber jump | 0.135 / 0.05 | 0.043 / 0.11 | 8.39 dB, no message |
| OD, charge-only node + snubber 10 ohm / 1 nF | refused | refused | refused | refused |
| same, snubber 50 ohm / 1 nF | refused | refused | refused | refused |
| same, snubber 10 ohm / 100 pF | refused | refused | refused | refused |
| ODM, charge-only node + snubber | refused | refused | refused | refused |

### What was wrong in s117

1. **The settle window could have no step in it.** The order-1 window after an event was a time (30 tau of the node, or 2 event-step widths
   for the inductor-routed jump). After a located event (vh > 0: the event step is ~30 ps and the next step is twice as long) the next
   accepted point lands exactly on the end of a 2-step-width window, so no order-1 step followed the event; the derivative history of a jump
   taken over a 30 ps step is then carried on by the trapezoidal rule: i_L gets the jump three times (9.54 dB = 3.00x). With a small node
   capacitor 30 tau is below the next step: 0.42 dB at 10 pF, 0.25 dB at 100 pF. **Now:** the window is the event step and exactly the next 2
   accepted steps (`PSS_SALT_SETTLE_MAXSTEPS`), no time rule. This completes card #568, which only stopped the window being too long.
   `PSS_SALT_MAXSTEPS=-1` restores the time rule (red proof).
2. **An RC snubber at a chargeless routed node had no jump of its own.** Its capacitor's rate (v_node - v_far) G / C changes with the node voltage
   exactly as the inductor's does, so it gets G / C x dV at the event (8.4 dB -> 0.135 dB; `PSS_SALT_NOSNUB=1`: 0.96 dB). What remains converges
   with the trajectory tolerance over the snubber's 10 ns: 0.135 dB at reltol 1e-4, 0.043 at 1e-6, 0.001 at 1e-7.
3. **The charge jump of card #559 ran outside what it is good for.** Its guard asked for "a native capacitor at the node"; with 10 pF native
   beside a 300 pF junction and no snubber it ran and was 0.43 dB off. The charge jump is right where the linear capacitor governs the node
   (1 nF native: 0.033 dB with or without snubber) or a snubber damps it; without either it is 0.13-3.6 dB off on every deck tried, including a
   node with 1 nF of pure junction capacitance whose transient the steps resolve well (so "unresolved transient", the reason given in the #559
   notes, is not it: the error is made where the diode takes over the inductor current within a step or two, and the linearisation along the
   coarse trajectory loses part of the injected charge there).

### The fast-node treatment

For a node whose charge is (mostly) a device's, next to exactly one inductor, the slow state the event changes is the inductor current, by the
integral of the time-shifted node voltage, (V_before - V_settled) / L x dt_e. The existing routed jump uses the node voltage at the accepted
event point for V_settled, which is right for a chargeless node only. A charged node is still on its way: each following accepted step adds
(V_j - V_{j-1}) / L as a further routed jump with the same event function (`PSSsaltCont` in pssshoot.c), while the node moves faster than
1e-3 of the fastest rate seen since the event (`PSS_SALT_CONT_FRAC`), until two quiet steps. A chargeless node is quiet at once and nothing
is added. V_before is the node equation of the pre-event topology (resistors, switches, diodes) as before; with any other device at the node
(the MOSFET of ODM) it is the node voltage of the previous accepted point. `PSS_SALT_NOCONT=1` (red proof): 32 dB and 27 dB off.

Which treatment an event gets (generic td path, i.e. device charges present):

| at the switch's power node | treatment |
|---|---|
| native capacitor >= 3x the smallest device capacitance of the node on the orbit, or native capacitor + RC snubber | charge jump (#559) |
| device charge, one inductor at the node, no RC snubber | fast node (routed jump + continuation) |
| device charge only, RC snubber at the node | refused (`PAC_SALT_FASTSNUB=1` runs it: 0.07-0.11 dB / 0.57 deg on 4 decks) |
| device charge, no inductor at the node, no native capacitor | refused as before (`PAC_SALT_QONLY=1`) |

### Limits and what was not done
- The fast-node phase error grows with the node's capacitance: 0.07 deg at 50 pF, 0.15 at 300 pF, 0.34 deg at 1 nF of junction capacitance at
  20 kHz (12 V, 1.2 A, 10 uH, 500 kHz); beyond about 1.5 nF of device capacitance without a native capacitor it would pass 0.5 deg. Not guarded.
- All 20 decks are one topology family (asynchronous buck, one switch, 500 kHz, control-to-output). Synchronous half-bridges with two ideal
  switches and device charge, `.pstb` and `.pxf` on these decks, and closed-loop decks were not run.
- A snubber at a charge-only node: the fast transient and the snubber's do not separate; adding the snubber jump to the continuation made it
  worse (+0.4-0.75 dB). Refused.
- The remaining 0.02-0.03 dB common to the fast-node decks was not traced.
- The `.pss` sensitivity recursion uses the same events; in a `.pss` pass the node is not known as charged, so it gets the plain routed jump
  (convergence only, the orbit is the same); Floquet multipliers of such decks were not checked against an independent value.

Tests: `tests/test_pac_charge.py --only 575`: 13 checks on 10 decks and 6 red proofs (`PSS_SALT_NOCONT`, `PAC_SALT_QONLY` x2,
`PSS_SALT_MAXSTEPS=-1`, `PSS_SALT_NOSNUB`), plus the refusal with its override in `check_salt_refusals`.

## Boost and half-bridge switch events: diode-fed capacitor, output impulse, hand-over (card #579, 2026-10-06)

Patch 0120, build `aaba50ba3f49-s119`. The matrix method of card #575 on other circuits: a boost (ideal low-side `sw`, output diode, 22 uF with
10 mohm ESR, 5 V in, duty 0.39), a synchronous buck with two ideal switches and dead time (58 ns / 1.8 ns) and a body diode, and the R6od buck
line-to-output, at duty 0.7 and in discontinuous conduction. Each deck has its own `tools/pac_perturbed.py` reference at 20/50/100 kHz
(`~/Work/data/ngspice_augmentation/579/m`, `matrix2-*.log`). The middle columns were produced on the candidate build of the same source
(`aaba50ba3f49-579-a`); ten of the decks are re-run on s119 by `tests/test_pac_charge.py --only 579`.

| deck | s119: worst dB / deg, reltol 1e-4 | reltol 1e-6 | s118, reltol 1e-4 |
|---|---|---|---|
| boost, chargeless diode (rs 50 mohm), vh 5 mV | 0.018 / 0.22 | 0.018 / 0.22 | 71.555 / 43.93 |
| boost, chargeless diode, rs = 0 | 0.017 / 0.22 | 0.017 / 0.22 | 5.292 / 72.05 |
| boost, chargeless diode, vh = 0 (event not located) | 0.064 / 0.74 | 0.021 / 0.42 | 71.450 / 43.70 |
| boost, diode cj 100 pF | 0.015 / 0.04 | 0.015 / 0.04 | not run |
| boost, diode cj 1 nF | 0.013 / 0.14 | 0.014 / 0.14 | not run |
| boost, diode cj 100 pF + tt 2 ns | 0.033 / 0.18 | 0.030 / 0.23 | 71.571 / 44.00 |
| boost, cj + tt, rs = 0 | 0.071 / 0.28 | 0.080 / 0.32 | not run |
| boost, chargeless diode, 1 nF at the switch node | 0.183 / 0.55 | 0.175 / 0.67 | 0.183 / 0.55 |
| boost, cj + tt, 1 nF at the switch node | 0.136 / 0.74 | 0.105 / 0.61 | 0.136 / 0.74 |
| synchronous buck, dead time, chargeless body diode | 0.040 / 0.00 | 0.040 / 0.00 | 0.040 / 0.00 |
| synchronous buck, charged body diode (cj 300 pF, tt 2 ns) | 0.024 / 0.10 | 0.027 / 0.09 | 0.253 / 0.10 |
| synchronous buck, chargeless, 1 nF at the node | 0.064 / 0.05 | 0.039 / 0.07 | 0.064 / 0.05 |
| synchronous buck, charged, 1 nF at the node | 0.039 / 0.04 | 0.037 / 0.04 | 0.039 / 0.04 |
| R6od buck, chargeless, line-to-output | 0.004 / 0.01 | 0.004 / 0.01 | 0.004 / 0.01 |
| R6od buck, charged, line-to-output | 0.016 / 0.11 | 0.013 / 0.11 | 0.014 / 0.11 |
| R6od buck, chargeless, duty 0.7 | 0.003 / 0.00 | 0.003 / 0.00 | 0.003 / 0.00 |
| R6od buck, charged, duty 0.7 | 0.024 / 0.02 | 0.026 / 0.02 | 0.031 / 0.02 |
| R6od buck, chargeless, discontinuous conduction (60 ohm load) | 0.399 / 0.40 | 0.400 / 0.40 | 0.399 / 0.40 |

**On s118 and every build before it a boost with an ideal switch and no capacitor at the switch node returned a td `.pac` of essentially zero**
(71 dB below the reference) with only the generic "inductor-routed jump" warning. Four things were missing:

1. **The pre-event node equation did not contain the output diode.** `pss_salt_ind_jump()` solves the node's equation in the pre-event topology
   to get the voltage the node would have had without the flip. A diode is taken into it when its cathode is at the node, or its anode with
   rs = 0; with the anode at the node behind a series resistance (internal node) it was skipped, V_prev came out as -I_L / g_off and the jump
   1e5 times too large on the turn-on edge. For such a diode, for any device the equation does not model (#575), and when the node's previous
   event is still being continued, V_prev is now the node voltage of the previous accepted point. `PSS_SALT_NO579=1`: 71.6 dB.
2. **The capacitor a diode at the switch node feeds had no jump.** In a buck the output capacitor's current is continuous at the event; in a
   boost the diode current into it switches between 0 and I_L, so its rate jumps by (i_D after - i_D before) / C (`pss_salt_diode_far()`).
   `PSS_SALT_NODIOFAR=1`: 0.22 dB and 21 deg at 20 kHz.
3. **The output itself is discontinuous.** v_out = v_C + ESR x i_C steps by ESR x I_L (17 mV) at each event. When the event moves by dt_e the
   output changes by a pulse of area (y+ - y-) x (-dt_e): an impulse that no state carries, so the recursion never shows it. Its sideband-0
   weight J x coeff x exp(-j w t_e), J being the large-signal change of each row across the event step (and across each continuation step), is
   now added to the demodulation integral for events at an algebraic or fast node (`pac_imp_*` in pacshoot.c); with a capacitor at the switch
   node nothing but the switch current is discontinuous and nothing is added. It grows with frequency against a falling state response:
   `PAC_NO_IMPULSE=1` is 0.08 dB at 20 kHz, 0.36 dB at 50 kHz, 0.72 dB / 6.4 deg at 100 kHz.
4. **The continuation of card #575** gave per-step differences and dropped the quiet steps. It now gives the difference to the last GIVEN value
   (node voltage for the inductor, diode current for the fed capacitor), flushes once at the end, and when the node flips again (dead time)
   it is flushed up to the previous accepted point and ended before the new event is formed (hand-over: 0.25 -> 0.027 dB on the charged
   synchronous buck). It also runs for a charged diode whose junction sits behind its series resistance (no charge on the node's own row).
   The diode-current part ends when the current moves no faster than 10x the inductor current it follows (with a rate threshold relative
   to its own peak it never ended, and the whole off-time ripple was given to the output capacitor: 6 dB).

### Outside the bounds, not fixed
- **A capacitor at the boost's switch node** (1 nF): 0.14-0.18 dB / 0.55-0.74 deg, with and without diode charge, same as s118. The node is
  then a state (capacitor / charge jump); the error is made where the diode takes over within a step or two after the slew.
- **Discontinuous conduction**: 0.40 dB / 0.40 deg, independent of reltol. No switch event is involved in the diode's turn-off at zero current;
  the inductor current is then a stiff state under the trapezoidal rule (L x g_off). Not investigated; neither side of the comparison has been
  checked against a third method there. No message is printed: treat td small-signal results of DCM decks as unvalidated.
- **vh = 0 on the boost** (event not located, 10 ns event steps): 0.064 dB / 0.74 deg at reltol 1e-4, 0.021 / 0.42 at 1e-6: J and the diode-current
  difference then include the smooth change over the step. Use vh > 0 or a tight reltol.
- Not run: boost charged in DCM (the reference tool returned no reliable point), `.pstb` / `.pxf` / closed loop on any of these decks, buck-boost,
  isolated topologies, more than one inductor at the switch node (refused as before).

Tests: `tests/test_pac_charge.py --only 579`: 10 checks on 8 decks, 4 red proofs.

## Boost with a capacitor at the switch node; discontinuous conduction (cards #580, #581, 2026-10-07)

Patches 0121 (#581) and 0122 (#580), first built as `aaba50ba3f49-s121`, merged with patches 0123-0124 as `aaba50ba3f49-s123`. The two decks card #579 left outside the bounds.

### A capacitor at a boost's switch node (#580)

A switch node with a capacitor is a state: the event gives it the capacitor jump (or the charge jump of #559) and the recursion carries the
node from there. On a boost that was 0.14-0.18 dB / 0.55-0.74 deg off the perturbed PSS, with the same error at 200 and 2000 points per
period and with 1 or 4 order-1 steps, so neither a step nor a window effect. What the capacitor path does not have is what card #579
added to the routed path for a boost: the jump of the capacitor the diode feeds, and the output impulse. A capacitor node **whose diode
feeds another capacitor, and that has no RC snubber**, is now treated as a fast node (routed jump + continuation):

| deck (boost of #579 + capacitor at the switch node) | s119 worst dB / deg | s121 / s122, reltol 1e-4 | reltol 1e-6 |
|---|---|---|---|
| chargeless diode, 1 nF | 0.183 / 0.55 | 0.022 / 0.10 | 0.021 / 0.09 |
| cj 100 pF + tt 2 ns, 1 nF | 0.136 / 0.74 | 0.033 / 0.44 | 0.028 / 0.43 |
| chargeless diode, 100 pF | not run | 0.006 / 0.19 | 0.008 / 0.20 |
| cj + tt, 3 nF | not run | 0.053 / 1.04 | 0.041 / 1.02 |

The last row is the limit of the fast-node treatment: its phase error grows with the capacitance of the node (the slew takes C dV / I;
3 nF x 8.3 V / 1.7 A = 15 ns of a 2 us period). Not guarded. A buck node (the diode goes to ground or to a source: nothing is fed) keeps the
capacitor jump, and so does any node with an RC snubber; all 40 other matrix decks are unchanged (`579/matrix1-f.log`, `matrix2-f.log`).
Forcing the fast-node treatment on every capacitor node (`PSS_SALT_FORCEFAST=1`) is wrong: 4.7 dB on the R6o half bridge, 1.9 deg on a
buck with a snubber. `PSS_SALT_NOFARFAST=1` restores the capacitor jump (red proof: 0.183 dB at 100 kHz).

### Discontinuous conduction (#581)

The 0.40 dB / 0.40 deg of the DCM buck is step resolution. When the inductor current runs out, the instant the diode stops conducting falls
inside an accepted step; no switch flips, the step control does not look for it, and the orbit (and with it the small-signal result) is first
order in the step, whatever reltol is. The reference (`tools/pac_perturbed.py`, 800 -> 1600 points per cycle, step-converged) is right:

| points per period (4th field of `.pac`) | 200 | 2000 | 20000 |
|---|---|---|---|
| R6od DCM (60 ohm), error at 20 kHz, dB | 0.400 | 0.061 | 0.027 |

At 2000 points the error is 0.061 dB at 20, 50 and 100 kHz alike. The engine now says so once per run: when an inductor's current stays
below 1e-3 of its largest value for more than 1 % of the period it prints "PAC-TD: note -- the current of inductor L1 is zero for part of
the period (discontinuous conduction) ... Raise the points of the .pss/.pac line." It fires on the DCM deck and is silent on the
continuous-conduction decks tried (9 by hand, 2 in the test suite). `PAC_NO_DCM_NOTE=1` silences it.

### Not done
- Locating the diode's turn-off in the transient itself (it would remove the first-order error at 200 points); the note is the whole fix.
- The note is printed by the td small-signal engines only, not by a bare `.pss`.
- DCM with device charge has no reference (the reference tool returned no reliable point); DCM on a boost was not run.
- `.pstb` / `.pxf` / closed loop on these decks.

Tests: `tests/test_pac_charge.py --only 580`: 9 checks, 2 red proofs.

## Half bridges without dead time; closed-loop `.pstb` on the switch-event treatments (card #582, 2026-10-07)

Patch 0123 (first built as `aaba50ba3f49-s122`, merged with patch 0124 as `aaba50ba3f49-s123`). Cards #575-#580 validated td `.pac` open loop. `.pstb` goes through the same events with the loop
closed, and had only ever been checked on the R6s half bridge (node capacitor + snubber). Five closed-loop bucks of the R6s family against
`tools/pac_perturbed.py --loop-gain` (tran backend, 0.01 V; `~/Work/data/ngspice_augmentation/582/m`):

| closed-loop deck | engine - reference, T at 30 kHz | 10 / 100 kHz |
|---|---|---|
| synchronous, 1 nF + snubber at the node (the handbook's example 23) | -0.001 dB / 0.00 deg | +0.006 dB / -0.07 deg; +0.009 dB / -0.02 deg |
| synchronous, no capacitor at the node | +0.024 dB / -0.02 deg (s121: +41.9 dB) | no reference |
| asynchronous, chargeless diode, no node capacitor | +0.009 dB / +0.08 deg | 100 kHz: -0.08 dB / -0.15 deg vs a 0.03 V reference |
| asynchronous, charged diode (cj 300 pF, tt 2 ns), no node capacitor: fast node | -0.037 dB / +0.02 deg | no reference |
| the same with 1 nF + snubber | +0.008 dB / +0.06 deg | no reference |

The reference tool refused the 10 and 100 kHz points of the four new decks (its amplitude-halving check: 0.1-7 dB between amplitudes at
0.003, 0.01 and 0.03 V on the tran backend); only 30 kHz, and one 100 kHz point at 0.03 V, are checked.

**Found: two switches flipping in the same step.** In a half bridge driven without dead time both switches flip at the same instant. With
a capacitor at the node each event is a capacitor jump and the two add up correctly. Without one, each switch formed its own routed jump,
taking the other switch's conductance as already flipped and the pre-event node equation as if only this switch had moved. On s121 and
every build before it, with no message beyond the routed-jump warning:

| deck (no capacitor at the switch node, no dead time) | s121 | s122 / s123 |
|---|---|---|
| open loop, chargeless body diode, td `.pac` 20/50/100 kHz | 93.8 dB off | 0.001 dB / 0.00 deg |
| open loop, charged body diode | 12.6 dB off | 0.003 dB / 0.00 deg |
| open loop, no diode at the node | 95.5 dB off | 0.001 dB / 0.00 deg |
| closed loop, `.pstb` at 30 kHz | 41.9 dB off | 0.024 dB / 0.02 deg |

Now the first flipping switch of a node (list order) gives the node's jump, with the node voltage of the previous accepted point as the
pre-event value; the other gives nothing at that node (no routed jump, no charge jump), and the output impulse treats the pair as one
event coefficient. `PSS_SALT_NO582=1` restores the per-switch jumps.

**Not understood: the legacy recursion on a node with nothing but switches.** With the events corrected, the deck without any diode at the
node was still 12.4 dB / 94 deg off on the legacy per-C/L recursion, although its events and its orbit are identical to those of the deck
with a (never conducting) chargeless diode, which is right. Forced through the generic charge-history recursion (`PAC_FORCE_GENERIC=1`)
it is right, and the 15 chargeless matrix decks of #575/#579 are equal or better there (boost phase 0.22 -> 0.04 deg). The generic recursion
is therefore selected for every circuit that has an ideal-switch power node with neither a native capacitor nor a voltage source on it
(`pac_switch_node_without_cap()`; `PAC_NO582_GENERIC=1` gives the legacy path back: red proof, 12.4 dB). The defect in the legacy path
itself was not found; decks whose switch nodes all carry a capacitor still use it and are validated as before (R6 family).

### Not done
- Reference points at 10 and 100 kHz for the closed-loop decks; `.pxf` on any of the matrix decks; a closed-loop boost.
- The cause of the legacy recursion's error (above).
- Two switches at one node that flip in the same step but belong to different control functions (not a complementary pair): the first
  one's coefficient is used for the node's jump; not constructed, not tested.

Tests: `tests/test_pac_charge.py --only 582`: 5 checks (3 open loop, 2 closed loop), 3 red proofs.

## Routed jumps after the derivative history; the fast-node phase error (card #583, 2026-10-07)

Patch 0124, build `aaba50ba3f49-s123`. After card #580 the fast-node treatment was right in gain everywhere but its phase error at low
frequency grew with the capacitance of the switch node and the voltage it swings: 0.07 deg at 50 pF and 0.34 deg at 1 nF of junction
capacitance on the buck, 1.04 deg on the boost with 3 nF, 1.06-1.38 deg at 20 kHz on an inverting buck-boost (20 V swing) with 1 nF.

**Cause.** A jump added to the solution before the charge and its derivative history T are formed puts dy / h into T, and the trapezoidal
rule carries that spike into the next step; that is why every jump was followed by order-1 steps (cards #494, #568, #575). The
continuation of a fast node gives a jump at every step of the node's transition, so the whole transition (20-40 steps at 1 nF) ran at
order 1 - and card #568 had already measured that every order-1 step on this recursion costs accuracy.

**Fix.** In the generic recursion a routed jump (inductor branch row, snubber- or diode-fed capacitor, continuation step) is now added
after T_k has been formed from the smooth part (y += dy, Q += C dy), exactly as the charge jump of card #559 is; T carries no spike and no
order-1 window is opened for these events (`PSSsaltNoRoutedWindow`). Capacitor jumps at the switch's own node keep the 2-step window. The
generic recursion is selected for every circuit with a diode at a switch node as well (#582 selected it for capacitor-less switch nodes).
`PAC_SALT_SMOOTHT=0` restores jumps-before-T with the window.

| deck | worst phase at 20/50/100 kHz, s122 | s123 | worst dB, s123 |
|---|---|---|---|
| buck SWC, cj 1 nF | 0.34 deg | 0.10 | 0.005 |
| buck SWC, cj 300 pF + tt 5 ns | 0.15 | 0.07 | 0.016 |
| buck R6od, no node capacitor | 0.15 | 0.04 | 0.007 |
| buck, MOS body diode | 0.14 | 0.02 | 0.002 |
| buck R6od, 10 pF native, no snubber | 0.15 | 0.04 | 0.006 |
| synchronous buck, charged body diode | 0.10 | 0.01 | 0.035 |
| boost, cj + tt, 1 nF at the node | 0.44 | 0.06 | 0.029 |
| boost, cj + tt, 3 nF at the node | 1.04 | 0.32 | 0.047 |
| inverting buck-boost, chargeless | 0.05 | 0.05 | 0.008 |
| inverting buck-boost, cj 100 pF + tt 2 ns | 0.19 | 0.09 | 0.016 |
| inverting buck-boost, chargeless, 1 nF at the node | 1.06 | 0.14 | 0.039 |
| inverting buck-boost, cj + tt, 1 nF at the node | 1.38 | 0.25 | 0.014 |

The references' own amplitude-halving spread in phase is 0.03-0.09 deg at 20 kHz on these decks. All 47 open-loop matrix decks and the 5
closed-loop ones were re-run on the candidate build of this source (`aaba50ba3f49-583-a`; `579/matrix1-h.log`, `matrix2-h.log`,
`582/matrix3-h.log`): gain within 0.067 dB everywhere except the two known cases (chargeless buck with a snubber only: 0.149 dB at reltol
1e-4, 0.035 at 1e-6; discontinuous conduction at 200 points: 0.40 dB with the note), phase within 0.32 deg except the boost with vh = 0
(0.55 deg at reltol 1e-4).

### Not done
- The legacy per-C/L recursion is unchanged (it is still used when every switch node has a capacitor and no diode: the R6 family).
- The remaining phase on the largest node capacitors (0.25-0.32 deg) was not traced further.
- The inverting buck-boost with cj + tt and no node capacitor has no 20 kHz reference (the tool's amplitude-halving check refused it).

Tests: `tests/test_pac_charge.py --only 583`: 8 checks, 2 red proofs (phase).

## The legacy recursion is no longer used with ideal switches (card #585, 2026-10-07)

Patch 0125 (first built as `aaba50ba3f49-s124`, merged with patch 0126 as `aaba50ba3f49-s125`). Card #582 left one defect avoided but unexplained: on a half bridge whose switch node has neither a
capacitor nor a diode the legacy per-C/L recursion is 12.4 dB / 94 deg off (reference -7.8050 dB / -169.910 deg, legacy -20.2470 / -76.338,
generic -7.8049 / -169.910) although its events and its orbit are those of the deck with a chargeless diode at the node, which is right.
Narrowed with the legacy recursion forced (`PAC_NO582_GENERIC=1` on s123):

| change to the deck or the run | legacy result |
|---|---|
| a chargeless diode elsewhere in the circuit (across the load) | still wrong |
| a diode at the switch node, or 1 pF at the switch node (capacitor jump) | right |
| no output impulse; native history forced | still wrong |
| matrix reloaded and refactored before the legacy solve of every step | still wrong (so not a stale matrix) |

The cause was not found. Nothing ties it to the node being capacitor-less, so the generic recursion is now selected for **every** circuit
that contains an ideal switch (`pac_switch_node_without_cap()` returns 1 as soon as there is a Switch instance; `PAC_LEGACY_SW=1` gives the
narrower rule of #582/#583 back, `PAC_NO582_GENERIC=1` the legacy recursion). On the R6 family, where the legacy recursion was validated,
the two give the same values to 4 decimals at 20/50/100 kHz (R6o, R6ob, R6om, R6od, R6odm) at the same run time (0.29-2.06 s per deck).
`.pxf` with `adjoint=0` and `adjoint=1` equals `.pac` exactly on five matrix decks (with events the adjoint falls back to the forward
engine, so this is a consistency check).

Not done: the cause; circuits without ideal switches still use the legacy recursion when they have no device charges (smooth B-source
switches, linear circuits) - validated there since card #465.

## Sources must be periodic from t = 0 (card #586, 2026-10-07)

Patch 0126, build `aaba50ba3f49-s125`. A two-phase interleaved asynchronous buck was added to the matrix with its second ramp written the
natural way, `PULSE(0 2 {T/2} 0.97T 0.03T 1n T)`. td `.pac` came out 6.5 dB below the perturbed PSS at 20 kHz (the gain of one phase),
with the DCM note of card #581 as the only hint. Cause: every shooting period run starts at t = 0, and a PULSE with a delay sits at V1
for t < TD - so within each period run the second ramp is 0 for the first half period, the second switch is on for most of it, and the
engine converges, correctly, to the periodic orbit of that other circuit (one phase at zero current). A long `.tran` of the same deck
shows both phases at 1.36 A.

A delayed PULSE is periodic from t = 0 only if the pulse has ended PER - TD after its start (TR + PW + TF <= PER - TD): true for a delayed
gate pulse, false for a delayed sawtooth. The shooting engine now checks every independent V and I source before it starts and refuses

    PSS-SHOOT: REFUSED -- source vramp2 is not periodic from t = 0 (PULSE with a delay TD: before TD it sits at V1, its periodic
    continuation does not (TR+PW+TF > PER-TD); TD = 1e-06 s). ...

and likewise a SIN with TD > 0. `PSS_ALLOW_DELAYED_SOURCE=1` runs anyway (red proof: 6.6 dB off). **How to write a phase-shifted
sawtooth that is periodic from 0:** the unshifted sawtooth in series with a square source,

    Vramp2 ramp2 rsq PULSE(0 2 0 {0.97*tclk} {0.03*tclk} 1n {tclk})
    Vsq    rsq   0   PULSE(1 -1 {0.485*tclk} 1n 1n {0.5*tclk} {tclk})

With that the interleaved buck (chargeless diodes, 100 mohm phase resistors) is within 0.001 dB / 0.01 deg of the perturbed PSS at
20/50/100 kHz, at reltol 1e-4 and 1e-6. With charged diodes (cj 300 pF, tt 2 ns) it is within 0.037 dB / 0.02 deg, and with 1 nF at each
switch node as well within 0.047 dB / 0.04 deg (`579/matrix2-il2c.log`, s125).

Not done: PWL, EXP, SFFM, AM and behavioural (B) sources are not checked; a PULSE whose PER is not the analysis period or a divisor of it is
not checked here either; running the period from a time offset where every source is periodic (which would make the delayed form work)
was not attempted.

Tests: `tests/test_pac_charge.py --only 586`: 3 checks, 1 red proof.


## Upstream fixes since the base; VDMOS quasi-saturation in the transient Jacobian (card #590, 2026-10-08)

Upstream `pre-master-48` moved from the base `aaba50ba3f49` (2026-09-27) to `f961c98d0` (2026-10-07): 11 commits. The base pin is
unchanged; five fixes are taken on top of the series as cherry-picks with upstream authorship (`git cherry-pick -x`):

| Patch | Upstream | What |
|---|---|---|
| 0127 | fc601dcab | diode: instance parameter description text |
| 0128 | 797409238 | `s_xfer`: crash when `den_coeff` has length 1 |
| 0129 | 42f32880f | VDMOS `.ac`/`.pz`: the quasi-saturation term of the drain resistor was missing |
| 0130 | dfbbc64ca | MOS2: derivatives (gds, threshold dependence on vds) |
| 0131 | f961c98d0 | BSIM3v3 and BSIM4v5: capacitance derivative `dVdseffCV_dVb` |

**Not taken: 059066b53 ("Fix bug 873", random seed).** It makes `inp_getopts` hoist every `option ...` line of the deck, including the
ones inside `.control`, into the circuit's option list. A loop that sets `option temp = 0` and `option temp = 85` per point then runs
every point at 85 C (measured on `bench/decks/c_pvt_loop`: all 50 points at 85 C with the commit, 25 at 0 C and 25 at 85 C without).
That is a regression of the fix of card #474. Under project ruling #446 B (no upstream contact) it is not reported upstream. Also not
taken: 4 Visual Studio project-file commits and a merge commit.

**The same missing term in the transient load (patch 0132, ours).** With quasi-saturation (`rq`, `vq`) the VDMOS drain resistor depends
on vds. Upstream's `VDMOSload` puts the conductance in the matrix but not its derivative with respect to vds, so the transient Jacobian
is not the derivative of the residual. Newton still converges to the right waveform (the residual is right), but everything that reads
the Jacobian as the linearisation is wrong: td `.pac`, `.pstb`, `.pxf`, the shooting monodromy. Patch 0132 stamps the term `gdq` at
(D,d)+, (D,s)-, (DP,d)-, (DP,s)+ with the matching right-hand side. `NGSPICE_VDMOS_QSJAC=0` restores upstream's load.

Measured on `tests/fixtures/pac_vdmos_qs.cir` (a VDMOS in quasi-saturation at a DC operating point under a periodic source, so td `.pac`
must equal `.ac`): td `.pac` gave -0.02737 before, -0.04297 with the patch; `.ac` (with 0129) gives -0.04298 and a finite difference
of two operating points agrees. In the suite: worst relative error 1.4e-4 (target 2e-3); red proof with `NGSPICE_VDMOS_QSJAC=0`:
0.363 relative, 3.92 dB low.

No pinned value of the 41 suites or of the handbook moved on s131.

Not done: the derivative of the self-heating power with respect to the quasi-saturation term (`tnodeout` decks) is not stamped; no
VDMOS switching converter was compared against the perturbed PSS, only the fixture above; the transient waveform itself was not
expected to move and none of the suites shows a change, but no VDMOS transient deck is in the speed set, so the effect of the better
Jacobian on iteration counts was not measured.

Tests: `tests/test_pac.py` check `pac_vdmos_qs`: 1 check, 1 red proof.

## The td recording passes run at reltol 1e-8 on circuits without switch events (card #625, 2026-10-14)

Problem. Without `pacreltol=`, the recording passes of td `.pac/.pstb/.pxf/.pnoise` inherited the PSS tight reltol (1e-5, `pssshoot.c` tight_reltol). The recursion is
built on the transient's own adaptive grid, so its first-order error depended on how many steps the LTE controller happened to put on each fast edge. That grid
changes with any parameter perturbation, down to p*(1+1e-10). On `tests/fixtures/rd_pstb_r3d.cir`, |T| scattered by 4.0e-3 rms (up to 0.15 dB, PM 0.056 deg rms)
between neighbouring parameter values. That breaks finite-difference sensitivities (#475, `tools/pstbsens.py`) and any sweep. Newton tolerance and the PSS orbit
were measured and ruled out (`~/Work/data/ngspice_augmentation/625/findings.md`).

Change (patch 0134, `pacshoot.c`):
- When the deck sets no `pacreltol=` and there is no `PAC_RECORD_RELTOL`, the recording passes run at reltol **1e-8**.
- If a recording pass stops with "timestep too small", the analysis is redone at 1e-7, then 1e-6, then the inherited value, with a `Note:` line each time.
  rd_pstb_r4d (peak-current buck with a dead-time network) lands on 1e-6: its turn-on edge collapses the step below 1e-7, cause not established.
- Circuits with ideal switches (SW, CSW) or any XSPICE code model keep the inherited value. Their event corrections and hidden-state recursions were validated on
  that grid, and with 1e-8 test_pac_charge lost 5 checks (swc decks 0.018 -> 0.12 dB).
- `PAC_RECORD_RELTOL_AUTO=0` restores the old behaviour everywhere; `pacreltol=` still wins.

Results (s136 vs s133):

| deck | s133 (inherited 1e-5) | s136 (1e-8 default) |
|---|---|---|
| rd_pstb_r3d, scatter of abs(T) over 41 neighbouring cval | 4.0e-3 rms | 7.3e-7 rms |
| r4_pstb_d030 vs the perturbed-PSS reference (`pac_perturbed.py --loop-gain`), 50-200 kHz | 0.0017-0.080 dB / 0.08-0.46 deg | 0.0003-0.025 dB / 0.001-0.11 deg |
| R7 BSIM buck vs the converged FRA, 20-200 kHz, no `pacreltol` | +0.115..+0.122 dB | +0.014..+0.021 dB (= `pacreltol=1e-7`) |
| handbook 24 (.pnoise, smooth switch) 1 kHz vs PAC_RECORD_RELTOL=1e-10 | 2.3e-4 relative | 2.2e-6 |
| handbook 07 (.pac td) 1 kHz | 1e-12 relative change | |

Tests: test_pac `pstb_r4_d030_perturbed` (0.04 dB / 0.2 deg vs the perturbed-PSS table, red proof with PAC_RECORD_RELTOL_AUTO=0). The FRA comparison at
199 kHz is now excluded from `pstb_r4_d030`, because the FRA there is 0.32 dB / 3.4 deg off the perturbed reference. test_r7 and handbook 40: the plain deck is
within 0.03 / 0.05 dB of the FRA; the red proof is the legacy default. Handbook 07/08/09/24 snapshots were re-pinned.
test_pxf `pxf-adjoint(r3 buck)`: adjoint and direct agree to 6.4e-8 at 199 kHz (7e-10 on the legacy grid). The fine grid takes steps down to 1.4e-13 s, so this is
a roundoff floor; the answer itself moves 1.2e-4 between recording reltol 1e-7 and 1e-8. The target is 1e-7 on the default grid, and 1e-8 is kept on the legacy grid.
The first build of this change (s135, discarded) detected code models by `DEVload == MIFload`. That never matches a code model loaded from a `.cm` library (its DEVload is the library's wrapper),
so XSPICE decks got the fine grid and test_xspice's s_xfer check went from 1.7e-7 to 2.0e-3. s136 uses `DEVpublic.cm_func`, as the hidden-state code does.

Not done:
- the cause of the r4d edge collapse below 1e-7;
- tighter recording for ideal-switch and XSPICE decks (their event corrections were not re-derived for a fine grid);
- the wall cost on 261-frequency sweeps was not measured (r3d 14 points: +14 %);
- no run on a second workstation.

## Floquet re-evaluation step size (card #682, 2026-10-19)

The R4d buck wrapped in subcircuits reported `spectral_radius=1.266873 UNSTABLE(period-doubling)` on s140 while the same circuit written flat
read 0.994. Cause: the Floquet re-evaluation built its finite-difference columns with `h = max(|x| * eps, 1e-6)`, so a state whose magnitude
is below its own ripple (the subckt deck's logic-state nodes) got a column that was orbit noise divided by `h`; that rank-one block gave the
-1.27 multiplier. The columns now use the step rule of the Newton columns (`pss_newton_fd_h()`, 1 % of the state's ripple, card #505).
Result on s141: flat 0.994039, subckt 0.994025, both stable and repeatable. `PSS_SHOOT_NO_BIGH=1` restores the old step (the red proof in
`tests/test_floquet_subckt.py`, fixtures `tests/fixtures/floquet682/`).

Not done: no other subckt-wrapped converter was checked beyond R4d; no run on the second workstation.

## Half bridge with dead time: the body diode's series resistance; the fast-node tail (card #727, 2026-10-19)

Build `aaba50ba3f49-w727` (source branch w727 = s142 + one commit, `pssshoot.c`). Card #709 found the td `.pac/.pxf` control-to-output gain of
a synchronous buck with dead time (two ideal switches, body diode, no capacitor at the switch node; `tests/fixtures/pxfpstb709/syn_dt*.cir`)
+0.042 dB (chargeless diode) / +0.035 dB (cj 300 pF, tt 2 ns) above the perturbed periodic steady state at every frequency, phase right,
absent without dead time. It was the "0.040 / 0.00" synchronous-buck row of card #579.

**Cause 1 (both decks).** The inductor-routed jump takes the pre-event node voltage from the node equation of the pre-event topology
(`pss_salt_ind_jump()`). A diode with its cathode at the node entered it with the junction's voltage, current and conductance (CKTstate1),
as if the junction sat directly between the anode terminal and the node; the series resistance (between the anode terminal and the
internal node) was left out. For a conducting diode the pre-event voltage came out short by Id x RS. After a dead time both turn-on events
start from the conducting diode. `PSS_SALT_TRACE` on syn_dt (RS 50 mohm, s142):

| event (one period) | node voltage before (previous accepted point) | jump dV used, s142 | w727 |
|---|---|---|---|
| high side off, rising ramp | +11.988 V | -12.6181 | -12.6181 |
| low side on, rising ramp (after 58 ns dead time, Id 1.21 A) | -0.6288 V | +0.5559 | +0.6167 |
| low side off, falling ramp | -0.0076 V | -0.5855 | -0.5855 |
| high side on, falling ramp (after 1.8 ns, Id 0.76 A) | -0.593 V | +12.5476 | +12.5854 |

The two events of one ramp edge share one event-time sensitivity, so what counts is their sum: 12.062 V instead of 12.000 V on the rising
edge (+0.51 %), 11.962 V on the falling edge (-0.31 %); weighted by the edges' share of the duty sensitivity (97 % / 3 %) that is the +0.04 dB.
This explains every observation of #709: proportional to Id x RS / Vin, independent of the dead time and of the diode's Is (the error is
Id x RS, not Vd), dependent on the ramp's fall fraction (0.0407 dB at 3 %, 0.0081 dB for a 50/50 triangle where the two edges nearly cancel),
absent without dead time (the switches then flip in one step and the jump uses the previous accepted node voltage, card #582). Fix: the
diode branch is linearised as a whole, G = Gd / (1 + Gd RS) about the branch voltage Vd + Id RS. An anode at the node behind RS already
used the previous accepted voltage (card #579).

**Cause 2 (charged diode only).** After cause 1 the charged deck was -0.007 dB. The fast-node continuation (card #575) hands the node's
further movement to the inductor until the node moves slower than 1e-3 of its fastest rate since the event. A charged diode keeps moving the
node for some tt after the slew (its diffusion charge fills) at a rate far below 1e-3 of the slew's (about 1e-5 on average): the continuation of the high-side turn-off ended
at -0.618 V while the node went on to -0.629 V before the low side turned on; those 11 mV never reached the inductor. The default fraction
(`PSS_SALT_CONT_FRAC`) is now 1e-4 at nodes without a native capacitor (1e-5 and 1e-6 measured: syn_dt_q within 0.0013 / 0.0023 dB, no
better). Capacitor nodes treated as fast (boost, buck-boost: card #580) keep 1e-3: with 1e-4 there the phase got worse (boost, cj + tt,
3 nF at the node: 0.32 -> 0.46 deg, over its 0.4-deg bound; bb_q_csw 0.25 -> 0.37 deg; bst_q_csw 0.06 -> 0.17 deg). `PSS_SALT_CONT_FRAC`
sets both.

`PSS_SALT_NO727=1` restores both old choices (red proof: +0.0418 / +0.0348 dB again, identical to s142).

| deck, Vctrl -> v(vout), error vs the perturbed PSS (dB) at 20 / 50 / 100 kHz | s142 | w727 |
|---|---|---|
| syn_dt (chargeless body diode) | +0.0418 / +0.0418 / +0.0417 | +0.0000 / +0.0001 / +0.0000 |
| syn_dt_q (cj 300 pF, tt 2 ns) | +0.0348 / +0.0351 / +0.0352 | +0.0005 / +0.0005 / -0.0003 |
| syn0 (no dead time) | 0.0000 | 0.0000 |
| buck_d70 / buck (asynchronous, same diode: the high-side turn-on) | -0.002 / -0.001 | 0.000 / 0.000 |

Phase unchanged (worst 0.012 deg on these decks). On the fast-node decks of `tests/test_pac_charge.py` the longer continuation also
removes most of the 0.01-0.02 dB residual that card #575 left "not traced" (worst dB / deg, s142 -> w727): swc_50p_1n 0.018 -> 0.006
(reltol 1e-6: 0.021 -> 0.006), swc_300p_5n 0.016 -> 0.009, od_nosnub 0.007 -> 0.002, od_q_d30 0.019 -> 0.008, syn0_q 0.004 -> 0.000,
syn_q at reltol 1e-6 0.014 -> 0.001, bst_q at reltol 1e-6 0.026 / 0.21 deg -> 0.028 / 0.09 deg; every other check within 0.001 dB and
0.01 deg of s142.

Tests (w727): `tests/test_pxf_pstb_switch.py` 44/44 (check 4 re-pinned from the +0.04 window to +-0.005 dB, two new red proofs
PSS_SALT_NO727 on syn_dt / syn_dt_q; s142 fails exactly those two windows, 42/44); `tests/test_pac_dcm.py` 6/6 (ccm_sync compared with its
#727 value; ccm_buck_d70 +0.0023 dB, ccm_buck +0.0007 dB vs s140, both toward the reference); `tests/test_pac.py` 35/35, output identical
to s142; `tests/test_r3d_r4d.py` 27/31 (4 info), numbers identical to s141; `tests/test_pac_charge.py` 136/137, the one not run is the
ASan check (NGSPICE_ASAN not set; s142 137/137 with it).

Not done:
- the remaining 0.034 dB / 0.09 deg of the charged boost (boost_q) and 0.018 dB of the chargeless boost in the same suite were not traced;
- a cathode-at-node diode with sidewall resistance (DIOtConductanceSW, a separate sidewall branch) is still taken without it; no deck has one;
- the continuation fraction was swept only on the half bridge (1e-4..1e-6); why the longer continuation hurts the phase at capacitor
  nodes was not investigated (they keep 1e-3); no ASan build, no closed-loop deck beyond the suites, no run on the second workstation;
- the doc rows of cards #579/#709 above keep their historical numbers.

## `.pss method=shoot` with a digital clock as the only clock (card #712)

A converter whose only clock is an XSPICE `d_osc` or `d_pwm` has no V/I source with a time function to take the period from; `s140` refused it
("needs a driven circuit"), and `driven=1` (card #494) only lifted that refusal.

* **Period.** `.pss Fguess ... method=shoot` (no `trigger=`) is admitted as driven, without `driven=1`, when a fixed-frequency `d_pwm`
  (`frequency`) or `d_osc` (all `freq_array` entries equal) runs at a whole multiple (1e-6 relative) of Fguess: the period is `1/Fguess` and the
  message `PSS: no V/I source with a time function; the XSPICE digital clock ... sets the period` is printed. A clock at another frequency is
  not taken (message names it); `driven=1` still forces. A control-dependent `d_osc` frequency is not recognised as a clock (give `driven=1`).
* **Clock edge on the period boundary.** A `d_pwm`/`d_osc` edge lands on `t = T` every period. (a) Whether it is processed at the end of one run
  or the start of the next depends on rounding of the continuing model time (1e-20 s); the digital node then differed between the two
  snapshots of the hidden-state periodicity check although the orbit is the same. The check now counts a pending output event due within 1e-15 s
  of the snapshot as already applied. (b) The same rounding put the edge 9e-19 s before the final time ("Timestep too small; time = T,
  timestep = 0"); event times, output events and breakpoints within 64 eps of the final time are snapped onto it (only when the continuing
  offset is non-zero, i.e. in the periodic analyses).
* **Hidden state that does not settle.** If repeated periods from the same native state do not bring the hidden node vector back (drift above
  1e-3 and not decaying after 4+ periods) in a circuit with digital (event) nodes, the run is now REFUSED ("does not settle"); a circuit
  with only analog code models keeps the #551 warning (the first build refused the XSPICE `capacitor` deck of test_xspice, whose hidden node
  sits near 0 V: 8e-12 V of drift, 1.6e-3 relative on the 1e-9 floor; `PSS_XS_SETTLE_ANY=1` restores that rule). Measured case: closed loop with the compensator output
  driving a `d_pwm` directly (`tests/fixtures/digclock/buck_cl_dpwm.cir`): the `d_pwm` decides the length of its next half period from the
  analog signal about 0.8 of a half period after the previous edge, and the edge time is state the shooting unknowns do not contain (the
  lagged hidden state multiplies by about 1.1 per period); before the refusal the run reported "converged" with the orbit 140 % of the ripple
  off the settled `.tran` (average 8e-4 off). Needs the pending edge times as Newton unknowns: not done.
* **Comparators through an `adc_bridge` in the loop** (analog ramp vs compensator, `buck_cl_dosc.cir`-style, not kept as a fixture): the
  edge time is quantised to the step grid (no breakpoint at the crossing), the finite-difference monodromy is wrong and the Newton step does not
  reduce the residual (0.75 after 150 iterations). Not solved.
* **Results** (`tests/test_pss_digclock.py`, build `aaba50ba3f49-w712`): open-loop 12 V -> 3.3 V buck, 500 kHz, `d_pwm` clock (duty 0.275 and 0.55)
  and `d_osc` clock: orbit equals the settled `.tran` to 0.04-0.05 % of the ripple (4.8 mV), average 4e-7..6e-7 relative. R9 (digitally controlled
  buck, 11-bit DPWM, `tests/fixtures/digclock/r9_dosc.cir`: `r9_fine.cir` with the PULSE sampling clock replaced by a `d_osc`), period 1.2 ms
  = 600 cycles: 0.415 % of the ripple (3.72 mV), average 5.4e-7 relative, 31 s.
* **Small signal.** td `.pac/.pstb/.pxf/.pnoise` through the digital blocks keep refusing (`PAC-TD: REFUSED ... event-driven / digital`), card #714.

## td `.pac` in discontinuous conduction at 200 points: the 0.40 dB removed (card #711, 2026-10-19)

The 0.40 dB / 0.42 deg of the DCM buck (card #581, "step resolution") was not the diode's turn-off at zero current. It was the **turn-on of the
ideal switch** at the end of the period, in the inductor-routed jump (`pss_salt_ind_jump()`, pssshoot.c). In discontinuous conduction the switch
node floats before the flip (switch off 1 Mohm, diode off), so the pre-event node equation V_prev = (g V_in + K - I_L) / (g + G) divides by a
conductance of 1 uS. It was given I_L from `CKTrhsOld`, the event step's value, which already holds the response to the flip over the event step
(V dt / L: 12 V x 30 ps / 10 uH = 36 uA at 200 points). Measured on `tests/fixtures/dcm711/dcm_buck.cir` (event at 1.9847 us):

| points per period | event step | I_L used | V_prev (true node voltage 5.7408 V) | 20 kHz, dB |
|---|---|---|---|---|
| 200 | 30 ps | 15.6 uA | -3.64 V | -30.784 (0.40 off) |
| 800 | 7.5 ps | 8.6 uA | 3.40 V | -31.07 |
| 3200 | 1.9 ps | 6.8 uA | 5.15 V | -31.14 |

Found by narrowing the maximum step to 1/16 in windows of the period (only 1.98-1.9847 us, the approach to the event, moved the result from -30.784
to -31.142 dB) and printing the equation's terms at the event. A tighter `pacreltol` (1e-9) hid it by shrinking the event step. The jump now uses
the inductor currents of the previous accepted point (rolled with the switch-terminal voltages in `PSSsaltFindEvents()`): V_prev = 5.7408 V at
every resolution. `PSS_SALT_NO711=1` reads the event step's current again (red proof).

| deck, 200 points | s140 error dB / deg | this build |
|---|---|---|
| R6od DCM 60 ohm, 20 / 50 / 100 kHz | 0.40 / 0.42 (worst of the three) | 0.022 / 0.009 (worst of the three) |
| same, 150 ohm, 20 / 50 kHz (new reference) | 0.41 / 0.21 | 0.024 / 0.003 |

The remaining 0.022 dB is the reference's own amplitude-halving residual (0.019 dB on the 150 ohm run). The result no longer depends on the points
(100 ... 20000: -31.167 ... -31.165 dB). Four continuous-conduction decks (buck, buck at duty 0.7, boost, synchronous buck with dead time) are
unchanged to 4 decimals at 20 and 100 kHz: with a conducting diode or second switch G is 10 S and I_L's step-to-step change does not matter.
The DCM note of card #581 is kept with new wording (no longer "first order in the step"). Test: `tests/test_pac_dcm.py` (6 checks, 2 red proofs).

Not done: DCM with device charge and a DCM boost (no reference, as in #581); `.pstb` / `.pxf` on DCM decks; a DCM deck whose switch node carries
a capacitor (it takes the capacitor-jump path, which does not use I_L).

## Charge-only switch node with an RC snubber; continuation jumps at the event coefficient (card #713, 2026-10-19)

Build `aaba50ba3f49-w713` (source branch w713 = s143 + one commit, `pssshoot.c`, `pssshoot.h`, `pacshoot.c`). Since card #575 a switch node
whose only charge is the devices' (diode cj / tt, MOS body diode, no native capacitor) together with an RC snubber was refused: the fast-node
treatment gave the snubber capacitor no jump and was 0.07-0.11 dB / 0.57 deg off then (before #583 and #727; on s143 with the override
`PAC_SALT_FASTSNUB=1` the same decks were 0.002-0.07 dB off, growing with the snubber's time constant).

**1. The snubber of a fast node.** The snubber capacitor's rate G (v_node - v_far) / C_far changes across the node's transition. Taking the
saltation rule across the whole (unresolved) transition, it gets one jump when the continuation ends: G / C_far x (D_end - D_before),
D = v_node - v_far, D_before from the pre-event node voltage of the routed jump, D_end at the last continuation point (the previous accepted
point on a hand-over). A chargeless node has no transition and this is the G / C_far x dV jump of #575; a snubber that settles together with the
node gets ~0 (turn-off of the 10 ohm / 1 nF deck: D -1.3 V at the end; turn-on: 12.4 V, the node is held by the switch while the snubber
charges). Driving the snubber with the node increments at every continuation step instead (what #575 tried) is +1.2 dB (10 ohm / 1 nF),
+0.8 dB (10 nF) and +8.4 dB (DCM) off: the snubber's own relaxation is counted twice, once in the increments and once by the recursion.

**2. Continuation jumps take the coefficient of their event.** Every part of a continued routed jump (inductor row, diode-fed and snubber
capacitors, the output impulse of card #579) belongs to the switch flip, i.e. to the event time, but its coefficient (y(posCntrl) -
y(negCntrl)) / h_dot was read from the small-signal solution at the step it was given in: the control's phase at a later time, a delay of the
order of the node's transition. Found on the DCM snubber deck: its turn-off edge is slow (0.33 A into cj + 1 nF, 47 ns) and the error grew
linearly with frequency (0.13 / 0.31 / 0.47 deg at 20 / 50 / 100 kHz), the same with or without the snubber jump; the charge jump of #559 on
the same deck (`PAC_SALT_QONLY=1`, applied at the event only) was right in phase. The flip's coefficient is now kept per frequency column under
an event serial (`PSSsaltEvent.serial`, `.cont` for continuation jumps) and reused by its continuation: 0.02 deg. `PAC_NO713_COEF=1` reads every
coefficient at its own step again (red proof).

| deck (R6od family unless said; all without a native capacitor at the switch node), worst of 20/50/100 kHz | s143 | s143 `PAC_SALT_FASTSNUB=1` | w713 |
|---|---|---|---|
| snubber 10 ohm / 1 nF | refused | 0.007 dB / 0.11 deg | 0.004 dB / 0.12 deg |
| snubber 50 ohm / 1 nF | refused | 0.010 / 0.04 | 0.002 / 0.04 |
| snubber 10 ohm / 100 pF | refused | 0.002 / 0.04 | 0.002 / 0.05 |
| snubber 50 ohm / 10 nF | refused | 0.013 / 0.04 | 0.001 / 0.05 |
| snubber 10 ohm / 10 nF | refused | 0.072 / 0.10 | 0.024 / 0.07 |
| duty 0.7 | refused | 0.031 / 0.06 | 0.014 / 0.05 |
| diode cj 50 pF, tt 5 ns | refused | 0.032 / 0.13 | 0.010 / 0.13 |
| MOS body diode + overlap (R6odm), 10 ohm / 1 nF | refused | 0.002 / 0.10 | 0.004 / 0.11 |
| R6odm, 50 ohm / 1 nF | refused | 0.003 / 0.05 | 0.002 / 0.05 |
| DCM (60 ohm), 2000 points | refused | 0.028 / 0.47 | 0.033 / 0.08 |
| 10 ohm / 1 nF and R6odm at reltol 1e-6 | refused | 0.008 / 0.11, 0.005 / 0.10 | 0.004 / 0.12, 0.004 / 0.11 |
| closed-loop .pstb (handbook 23 deck, asynchronous, cj 300 pF tt 2 ns, snubber), 10/30/100 kHz | refused | 0.009 / 0.04 | 0.002 / 0.045 |

References: `tools/pac_perturbed.py`, pss backend, the #575 settings (0.05 V, 800 points per cycle, halving and step checks), s143 binary
(`~/Work/data/ngspice_augmentation/713/ref_*`); the four decks #575 had are identical to its s117 references. The 10 nF deck's 20 kHz point is
amplitude-limited at 0.05 V (halving 0.75 deg; the engine was 1.0 deg off it): 20 / 50 kHz from a 0.0125 V run. DCM at 0.005 V (the tool refused
0.05 V). Closed loop: `--loop-gain --amplitude 0.01 --backend pss` (halving <= 0.0003 dB; the tran backend's halving is 0.05 dB there). The
remaining 20 kHz phase (0.04-0.13 deg on every snubber deck) is the reference's: it is 4/3 of the tool's own amplitude-halving difference at
that point on each deck (quadratic in the amplitude), e.g. 10 ohm / 1 nF halving 0.087 deg, engine -0.117; 10 nF halving 0.745, engine -1.006
against the 0.05 V value.

**What moved on the other decks (fix 2 only; fix 1 touches nothing that was not refused).** All continued (fast-node) decks: bucks by at most
0.05 deg and 0.001 dB; the boost / buck-boost decks with a capacitor at the switch node (card #580; slower edges) by more, all inside their bounds:

| deck, error vs its pinned reference at 50 / 100 kHz | s143 dB / deg | w713 dB / deg |
|---|---|---|
| buck-boost, cj + tt, 1 nF (bb_q_csw) | +0.006 / +0.12, +0.014 / +0.25 | -0.000 / +0.01, -0.006 / -0.02 |
| buck-boost, chargeless, 1 nF (bb_cl_csw) | +0.024 / +0.03, +0.040 / +0.14 | +0.016 / -0.06, +0.017 / -0.12 |
| boost, cj + tt, 3 nF (bst_q_csw3n) | +0.026 / +0.11, +0.047 / +0.32 | +0.013 / -0.10, +0.024 / -0.21 |
| boost, cj + tt, 1 nF (bst_q_csw) | +0.020 / -0.02, +0.029 / +0.06 | +0.016 / -0.09, +0.021 / -0.12 |
| boost, chargeless, 1 nF (bst_cl_csw) | +0.018 / +0.00, +0.005 / +0.10 | +0.012 / -0.06, -0.004 / -0.08 |
| boost, diode cj 1 nF, no capacitor (bst_cj1n) | +0.004 / +0.01, +0.015 / +0.03 | +0.006 / -0.00, +0.018 / +0.03 |

The gains move toward the references. The boost phases change sign; whether they got better is below what the references resolve at 100 kHz:
the 0.05 V references have an amplitude-halving difference of -0.39 / -0.57 deg there (bst_q_csw / bst_cl_csw), and 0.0125 V runs
(`713/refa_*`) a step residual of 0.08-0.29 deg; against the 0.0125 V values w713 / s143 are +0.07 / +0.24 deg (bst_q_csw), +0.21 / +0.39
(bst_cl_csw), -0.34 / +0.19 (bst_q_csw3n), -0.20 / +0.06 (bb_cl_csw), -0.09 / +0.19 (bb_q_csw) at 100 kHz, and at 50 kHz (step residual
<= 0.055 deg) -0.12 / -0.05, -0.11 / -0.05, -0.15 / +0.06, -0.05 / +0.05, -0.00 / +0.11. The output-impulse part of fix 2 is <= 0.05 deg of it.

Tests: `tests/test_pac_charge.py --only 713` (checks_713, red_proofs_713): 13 checks and 4 red proofs (`PSS_SALT_NO713=1` refusal back,
`PAC_SALT_FASTSNUB=1` on the 10 nF and the 50 ohm / 10 nF decks > 0.045 / 0.008 dB, `PAC_NO713_COEF=1` on the DCM deck > 0.3 deg at 100 kHz);
the old refusal check in check_salt_refusals now asserts that the deck runs with the note and that `PSS_SALT_NO713=1 PAC_SALT_FASTSNUB=1` runs
it. `PAC_CHARGE_OUT=<dir>` redirects the suite's outputs (fixtures rewritten into it). Full suite on w713: 153/154 (the ASan check not run,
NGSPICE_ASAN unset); s143 on the same suite: 138 pass, 1 FAIL (the refusal check), 15 could-not-run (the refused decks, ASan); every other
check identical except the continued decks above. `tests/test_pac.py` 35/35, `tests/test_pac_dcm.py` 6/6, `tests/test_r3d_r4d.py` 27/31
(4 info): output identical to s143; `tests/test_pxf_pstb_switch.py` 44/44, boost_q +0.001 dB, otherwise identical.

Not done:
- DCM with a snubber at 200 points per period is 0.072 dB / 0.21 deg off (2000 points: 0.033 / 0.08); the existing treatments on the same
  circuit are no better at 200 points (charge jump with 10 pF: 0.047 dB / 0.39 deg; chargeless node: 0.070 dB / 0.25 deg), and the DCM note
  of #581 is not printed for these decks (the ringing inductor current never stays below 1e-3 of its peak). Not traced.
- whether the boost phases at 100 kHz are now better or worse is not resolved (references, above); no FRA cross-check;
- only asynchronous bucks (and one closed loop) with a snubber at a charge-only node; a half bridge with a snubber at a capacitor-less node,
  a boost with a snubber, two snubbers at one node (handled, up to 4, not run), `.pxf` (same engine as `.pac`, not run separately);
- the `.pss` sensitivity recursion does not know the node is charged and gives the chargeless snubber jump (convergence only, as #575 noted);
- no ASan build, no run on the second workstation.

## Absolute .pnoise on BSIM4 vs a transient-noise ensemble (card #708, 2026-10-19)

Question: is the ABSOLUTE output noise of `.pnoise ... method="shoot" solver="td"` right on BSIM4 circuits, where the device noise comes
from the model's own noise equations? Binary s143 (`aaba50ba3f49-s143`). Test `tests/test_pnoise_bsim4.py` (9/9 pass, 2 info lines),
fixtures `tests/fixtures/pnoise708/`, ensemble tool `tools/pnoise_mc_bsim4.py`, all runs and logs in `~/Work/data/ngspice_augmentation/708/`.

**Circuits.** `cs_bsim4.cir`: BSIM4 4.8.2 common-source stage (the model card of `pnoise_track_hold_bsim4.cir`, capmod 2), W/L 1u/0.5u,
gate driven at 100 kHz through Rg 10k (0.6 <-> 1.4 V trapezoid, 1 us edges), Rd 20k, Cl 100p: each period the device goes from weak
conduction (vd 1.63 V, id 8 uA) to triode (vd 0.32 V, id 74 uA), so both the channel-noise density and the output pole are periodic.
`th_bsim4.cir`: track-and-hold, Rtrack 1k + BSIM4 switch W/L 1u/0.5u (Ron ~1k) onto 100p, fclk 100 kHz, 20 % track window; the switch
carries 73 % of the output noise power. Flicker is off in both (noia=noib=noic=0); `ntn` (ntnoi) and `tnm` (tnoimod) are parameters.

**Method.** BSIM4 has no transient noise, so each generator is replaced by an explicit trnoise-driven source, in two parts:
- *res*: BSIM4 channel noise off (ntnoi=0); each resistor noisy=0 with a TRNOISE current in parallel, white one-sided PSD 4kT/R.
  `.pnoise` reference: the same deck with noisy resistors. This checks the engine's absolute scaling with the BSIM4 dynamics in the loop.
- *mos*: resistors noisy=0 in both; the BSIM4 channel thermal noise. Its density along the orbit is taken from the model itself, without
  the `.pnoise` engine: the noiseless periodic orbit is simulated (maxstep T/2000, reltol 1e-6; settled, first-last point 0 / 1.6e-7 V),
  and at each of its ~2000 time points the transistor alone, its four terminals forced to the orbit voltages, gets a plain `.noise`
  (short-circuit drain-current noise through a 0 V source and a 1 ohm CCVS). In the ensemble the channel noise is a B-source
  I = sqrt(S(t)) w(t) between drain and source, w unit white TRNOISE, sqrt(S(t)) a periodic PWL (`r=0`) of that table: a quasi-static
  modulated source, which is what BSIM4noise evaluates at each step of the orbit. (A B-source reading `@m1[gm]` was not tried: for
  tnoimod 0 the density is 4kT ueff |qinv| / Leff^2 * ntnoi, not a function of gm, so a gm-based source would be the wrong density anyway.)
- Ensemble: 16 seeds x 20 ms, first 2 ms dropped, noisy and noiseless copies in one transient, boxcar average of their difference over
  200 ns (integ), Hann/Welch 8192 points, 50 % overlap, 320 segments; 9-bin bands (610.35 Hz bins) at 4.3 .. 42.7 kHz, power averaged and
  compared with `.pnoise` at exactly those bins (same arithmetic as `tools/pnoise_mc.py compare()`). Every amplitude is scaled by K = 100
  (outputs ~2 uV/rtHz, rms a few mV, far above the integ rounding floor) and the `.pnoise` value by K. trnoise step TS = 10 ns (cs,
  min tau ~200 ns) and 5 ns (th). Band 1-sigma, stated before comparing: +0.111/-0.114 dB (chi-squared, 320 segments x 9 bins correlated
  as in #509: n_eff = 320*9/1.93); the mean over 7 bands 0.043 dB.
- The ensemble averages the periodogram over many periods per segment (163.8), i.e. it estimates the time-averaged PSD, which is what
  `.pnoise` sideband 0 reports (this section's first part: Spectre `noisetype=timeaverage`); both one-sided, V/rtHz.

**Calibration of the chain** (trnoise normalisation NA = sqrt(S0/(2 TS)), boxcar, Welch): R 10k into C 1n, TS = 10 ns, against the
analytic 4kTR/(1+(f/fc)^2), which equals plain `.noise` on that deck to 1.7e-7: worst band 0.069 dB, mean -0.014 dB.

**Results** (`.pnoise` minus ensemble, dB in power, per band centre in kHz):

| case | 4.27 | 9.77 | 15.87 | 21.97 | 28.08 | 34.79 | 42.72 | worst | mean |
|---|---|---|---|---|---|---|---|---|---|
| cs res (Rd + Rg, BSIM4 noise off) | +0.225 | -0.239 | +0.028 | -0.106 | +0.110 | +0.124 | +0.020 | 0.239 | +0.023 |
| cs mos, tnoimod 0 | +0.104 | -0.001 | -0.059 | +0.050 | -0.106 | +0.087 | -0.006 | 0.106 | +0.010 |
| cs mos, tnoimod 1, with the inert `.noise` card (below) | +0.069 | +0.000 | -0.028 | +0.034 | -0.114 | +0.071 | -0.016 | 0.114 | +0.002 |
| cs mos, tnoimod 1, plain deck (**defect**) | -4.346 | -4.417 | -4.450 | -4.395 | -4.552 | -4.378 | -4.483 | 4.552 | -4.432 |
| th res (Rtrack, BSIM4 noise off) | -0.150 | +0.285 | -0.196 | +0.158 | -0.033 | -0.012 | +0.095 | 0.285 | +0.021 |
| th mos, tnoimod 0 (seeds 101..116) | +0.030 | +0.029 | +0.228 | +0.091 | +0.190 | -0.035 | -0.146 | 0.228 | +0.055 |

Absolute levels at 4.27 kHz (`.pnoise`, V/rtHz): cs res 1.892e-8, cs mos 1.361e-8, th res 1.318e-8, th mos 2.166e-8. Every band of the
valid cases is within 2.6 sigma and inside the 0.5 dB target; every mean is within 1.3 sigma of zero. Ensembles that share seeds and a
transfer path are correlated: th mos with seeds 1..16 gave the same pattern as th res (worst 0.283 dB, mean +0.021), so the pinned th mos
value is a second run with seeds 101..116; cs mos tnoimod 0 and 1 share seeds too (similar patterns). Ensemble wall time (2 ngspice):
cs res 6 min, cs mos 20 min, th res 12 min, th mos 25 min.

**Static check** (`static_bsim4.cir`: the cs stage at a DC gate bias, flicker on (fnoimod 1), the period from an unrelated 1 kHz clock):
`.pnoise dev=1` total and every generator (BSIM4 id, 1/f, rs, rd, Rd, Rg) vs plain `.noise`, 1k..100k: tnoimod 0 worst 4.7e-7 of the
total; tnoimod 1 with the inert `.noise` card 4.7e-7.

**Defect found (not fixed): BSIM4 tnoimod=1 under `.pnoise` loses the source/drain-resistance part of the channel noise.**
`b4set.c` (and `b4v5set.c`, `b4v6set.c`, `b4v7set.c`) creates the internal source/drain nodes only if `rsh`/geometry asks for them, or if
`tnoiMod == 1 && noiseAnalGiven` (v5/v6: `tnoiMod != 0`), and `noiseAnalGiven` is set only when the job list holds a job named
`"Noise Analysis"`, i.e. a plain `.noise` card. With the nodes, the tnoimod=1 model splits the channel noise into the `id` generator
(T2 - igsquare) and the `rs` (vds >= 0) / `rd` generator (conductance reset to 1e3 S, scaled by 1 + theta^2 gspr / IdovVds) that carries the
igsquare part; without them `rs`/`rd` are evaluated between identical nodes and contribute nothing. A `.pnoise` deck has no
"Noise Analysis" job, so the part is silently missing. Numbers: static deck at 1 kHz, `.noise` rs 2.2787e-8 V/rtHz, `.pnoise` (td) 0,
`id` equal (1.23154e-8); on cs mos the transient ensemble (whose table comes from `.noise`) is 4.43 dB above `.pnoise` (the missing part is
1.8x the remaining one). The frequency-domain `.pnoise` (default solver) has the same loss (static deck 1 kHz: 3.567042e-6 without, 3.567113e-6
with the card; `.noise` 3.567113e-6): the cause is in the device setup, not in `pacshoot.c`. Workaround: add an inert `.noise` card on a
separate node (`Vdum ndum 0 DC 0 AC 1`, `Rdum ndum 0 1k noisy=0`, `.noise v(ndum) Vdum lin 1 1k 1k`); then `.pnoise` is complete (table
rows above). tnoimod 0 and 2 are not affected in b4set.c (v4.8); v5/v6 also for tnoimod 2 (not run). A fix would set `noiseAnalGiven`
for the PSS/`.pnoise` job too (or always create the nodes for tnoimod=1); not done here. The test records this as two info lines
("KNOWN DEFECT"), which say "now within limits" once it is fixed.

**Red proofs** (both must fail the cs mos check, and do): `PAC_NOISE_STALE_STATE=1` (device noise from the end-of-period state, the
pre-#544 bug): -6.8 dB; a deck with ntnoi = 1.25: +0.98 dB.

Not done:
- flicker under a periodic bias is not compared with an ensemble (only the static check has flicker on): a modulated 1/f source needs a
  1/f trnoise per table point or a filter bank, and the engine's treatment of folded flicker (card #556) is a separate question;
- tnoimod 2 (correlated gate noise, `NevalSrc2`) and induced gate noise, igcmod/igbmod shot noise, rdsmod/rgatemod/rbodymod resistor
  noise: not in these decks;
- the trnoise step TS and the amplitude scale K were not varied for the BSIM4 decks (TS = tau/20 follows #509's measured rule; the 0.1 dB
  agreement of cs mos bounds both effects there);
- the ensembles are not re-run by the test (pinned values; `tools/pnoise_mc_bsim4.py` reproduces them in 6-25 min each);
- R7/R7p absolute level (card #544) still has no reference; the defect above is not fixed (src untouched); no ASan run or run on the second workstation.

## BSIM4 tnoimod noise nodes under .pnoise (card #729, 2026-10-19)

Fixes the defect found by card #708 (previous section). Source branch w729 (on s144 = 1cb1e1566), commit 0efdc3905; test build
`aaba50ba3f49-w729`.

**Mechanism.** `BSIM4setup` (`b4set.c`; the same block in `b4v5set.c`, `b4v6set.c`, `b4v7set.c`) set a local `noiseAnalGiven` by walking
`ft_curckt->ci_curTask->jobs` for a job whose `JOBname` is the string `"Noise Analysis"`, the name `dot_noise()` (inp2dot.c) gives a `.noise`
card. The internal drain/source nodes are then created if `rdsmod != 0`, or `tnoiMod == 1 && noiseAnalGiven` (v4.7/v4.8; v4.5/v4.6:
`tnoiMod != 0`), or `rsh`/geometry asks for them. `.pnoise` is a PSS job (`dot_pnoise()`, analysis "PSS", `JOBname` "Periodic Noise
Analysis", `PSSdoPnoise = 1`), so the nodes were not created and the rs/rd generator (`NevalSrc(..., sNodePrime, sNode, ...)` in b4noi.c)
sat between identical nodes. Both `.pnoise` engines (frequency-domain `pnoise_sweep()` and td `pac_td_pnoise()`) were affected; so was
`.sp ... donoise=1` (analysis "SP", `SPdoNoise`). No other device uses this pattern (`grep -rn "Noise Analysis\|JOBname\|ci_curTask"
src/spicelib/devices`: only the four BSIM4 setups; bsimsoi has a tnoimod but no noise-conditional node).

**Fix.** New `CKTnoiseAnalGiven(TSKtask *)` in `cktdojob.c` (declared in `cktdefs.h`): true if the task holds a job named "Noise Analysis"
(the old test), any NOISE-type job, a PSS job with `PSSdoPnoise`, or an SP job with `SPdoNoise`. The four BSIM4 setups call it instead of
their own loop. A transient, .ac, .op, .pss, .pac, .pxf, .pstb task answers no, so its node count is unchanged (byte-identical outputs,
below). `PSS_SALT_NO729=1` restores the name-only test (red proof).

**Numbers** (static deck, `tests/fixtures/pnoise708/static_bsim4.cir`, tnoimod 1, no workaround card):
- td `.pnoise dev=1` at 1 kHz: rs 0 (s144) -> 2.27867e-8 V/rtHz (w729), `.noise` 2.27867e-8; total 3.567042e-6 -> 3.567113e-6 (`.noise`
  3.567113e-6); every generator within 4.7e-7 of the total (the same as the workaround deck).
- frequency-domain `.pnoise` (harmonics 3) total vs `.noise`, 1k/10k/100k: v4.6.5, v4.7, v4.8.2 tnoimod 1: 2.0e-3 (s144) -> 6.8e-9 (w729);
  v4.5 tnoimod 1: 3.2e-7 -> 3.6e-8 (in 4.5 the rs part is ~1e-14 V/rtHz even under `.noise`); tnoimod 2 (all four versions): 1e-15 on
  both binaries. v4.5/v4.6 tnoimod 2 do create the nodes under `.noise` (condition `!= 0`), but the rs/rd generator is 0 there with this
  deck (rdsmod 0, no rsh), so no visible change; not exercised with rsh > 0.
- cs mos ensemble, tnoimod 1 without the card: -4.432 dB mean (s144) -> +0.002 dB, worst 0.114 dB (identical to the workaround deck).
- the td `.pnoise` engine refuses BSIM4v5/v6/v7 instances (card #509 scope), so those versions are covered only by the frequency-domain
  engine.

**Test** `tests/test_pnoise_bsim4.py`: the two KNOWN DEFECT info lines are now checks; added static tnoimod 2 (td, per generator incl.
`corl`), four `static fd` checks (one per BSIM4 version, tnoimod 1 and 2), and three red proofs with `PSS_SALT_NO729=1` (static tnoimod 1:
rs off by 6.3e-2 of the total; static fd 4.8.2: 2.0e-3; cs mos tnoimod 1: -4.43 dB), all failing their checks as required.
`PNOISE_BSIM4_DATA` overrides the output directory. w729: 19/19 pass. s144: 14/19 (the five #729 checks fail, as expected).

**No regression** (w729 vs s144): `tests/test_pnoise.py` 11/11 on both, all 11 output spectra byte-identical; cs_bsim4 at versions
4.5/4.6.5/4.8.2 x tnoimod 1/2: `.tran`, `.ac`, `.noise` and `.pss` outputs byte-identical (24 + 6 files). No other suite or deck in the
repository sets tnoimod, and with the default tnoimod 0 the setup is unchanged.

`.sp donoise=1` on a tnoimod 1 BSIM4 stage now includes the rs/rd part too (NF 27.52 -> 29.22 dB on a smoke deck); this is not checked
against an independent reference.

Not done: no ASan run or run on the second workstation; the handbook suite was not run (no handbook deck uses tnoimod); `.sp` noise not validated; the interactive
`qpnoise`/`phasenoise` commands evaluate noise on a circuit set up by `qpss`/`hbosc` with the deck's own task, so a BSIM4 tnoimod 1 device
there still lacks the nodes unless the deck has a `.noise`/`.pnoise` card (not fixed, not run).

## R7 (BSIM4 closed-loop buck) td .pac vs a perturbed periodic steady state (card #710, 2026-10-19)

R7 (`refs/r7_bsim_buck.cir`, card #539: 5 V -> 1.46 V, 2 A, 500 kHz, every switch, comparator, gate and amplifier a BSIM4 MOSFET) had
only FRA curves as a reference for its `.pac` (0.09-0.12 dB uncertainty; `.pac` sat 0.10-0.12 dB above them at 40-200 kHz). A first
`tools/pac_perturbed.py` run was stopped after 75 min on #587. Run to the end on the development workstation (one process, s144 binary, 61 min):
`pac_perturbed.py sweep --deck r7_ctrl.cir --inject source:Vctrl --outputs vout --fs 500e3 --freqs 20000,50000,100000 --amplitude 0.005`
(the deck is `tests/fixtures/r7_pac_ctrl2out.cir` up to its `.options reltol=1e-4`; backend pss, 200 -> 400 points/cycle).

| f | perturbed PSS (dB / deg) | halving / step residual | td `.pac` s144 (400 points) | difference |
|---|---|---|---|---|
| 20 kHz | 3.6001 / -7.4101 | +0.0012 dB; 0.0001 dB / 0.002 deg | 3.5984 / -7.3652 | -0.0017 dB / +0.045 deg |
| 50 kHz | -4.6896 / -76.8939 | +0.0003 dB; 0.0001 / 0.006 | -4.6861 / -76.8784 | +0.0035 / +0.016 |
| 100 kHz | -11.8511 / -86.4057 | +0.0005 dB; 0.0013 / 0.018 | -11.8460 / -86.3857 | +0.0051 / +0.020 |

The 20 kHz FRA value the R7 fixtures cite (1.502 V/V = 3.533 dB) is 0.067 dB below this reference, inside the FRA's own uncertainty;
the `.pac` / FRA offset noted above was the FRA's. `tests/test_r7_perturbed.py` pins the three points (bound 0.02 dB / 0.15 deg) with a
red proof: `PAC_NATIVE_HISTORY=1` (the device charge history of #542 off) reads -119.5 dB at 20 kHz. 4/4 on s144, about 95 s.
Not done: the R7 loop gain (`.pstb`) and line-to-output / output impedance against the same kind of reference; R7p (the extracted
variant); `.pnoise` absolute level of R7 still has no ensemble reference (the #708 method could give one).

## Boost / buck-boost with a capacitor at the switch node vs an exact reference (card #728, 2026-10-19)

Patches 0146-0147 (source commits `pac td (#728)` and `pss salt (#728)`). The #713 arbitration had built an exact reference for six
decks whose switch node carries a native capacitor: boost with 1 nF / 100 pF / 3 nF at the node, chargeless or cj 100 pF + tt 2 ns diode,
and the inverting buck-boost with 1 nF, chargeless or charged. On s144 they read the phase low (bst_cl_csw -0.098 deg, bst_q_csw3n
-0.156 deg at 100 kHz) and the buck-boost gain high (bb_cl_csw +0.021 dB at 20 kHz).

**Exact reference, `tools/pac_exact.py`** (deck texts in `tools/pac_exact_decks.py`; needs numpy + scipy): the circuit as an ODE with
ngspice's level-1 diode equations, the ideal switch's instants found exactly (brentq on the known control), segment-wise Radau at rtol
1e-12, periodic steady state of the LCM period by shooting, the Fourier coefficient at f as a quadrature state, symmetric difference in
the amplitude (4, 2, 1, 0.5 mV) and Richardson in A^2. Error bar <= 1.8e-5 dB / 6.3e-4 deg. `--dump` writes the time-domain small-signal
trajectory.

**Cause.** Not the #713 event coefficient (with the parts below fixed, the flip's coefficient is right for every continuation part; the
own-step coefficient of s143 gave +0.08..+0.37 deg). Found by dumping the engine's small-signal trajectory (`PAC_SS_DUMP=<file>`) and the
exact one, and comparing every state jump of the continuation edge by edge. Six defects, all in what the continuation gives:
1. far-capacitor share: a diode-fed capacitor's charge impulse divides between the capacitor branch (with its ESR) and the load; all of it
   went into the capacitor (0.12 % too large, -0.06 deg at 100 kHz above the boost's right-half-plane zero);
2. node correction: with the inductor current changing during the slew the shifted-slew model is off by dt_e (I(t) - I(t_e)) / C
   (inductor jump 4.4e-4 too large -> 4e-5);
3. output impulse with drift: the impulse held the output's own smooth drift (2-3 % of the ESR step; bst_cl_csw 100 kHz -0.070 ->
   -0.005 deg). Only for capacitor nodes: on the charge-only fast nodes of the bucks (#575) the drift-free impulse moved results AWAY
   from their perturbed-PSS references (R6od 100 kHz -0.014 dB / +0.06 deg, DCM with a snubber 1.8 deg), not understood; they keep the
   old impulse;
4. far-capacitor lump: during the diode's own fast transient (reverse recovery: -106 A and back on the charged boost) the parts are held
   back and given as one lump after it (was 30 % too large step by step; a lump at the far node alone made the clamped node ring);
5. pre-event voltage of a capacitor node: the routed jump's node equation was 0.23 V off at the charged buck-boost's turn-on (jump
   1.15 % short); a node with a native capacitor takes its previous accepted voltage;
6. delay terms: parts accrue over their step but were given at its end (0.2 % of the output capacitor's charge at 200 points).

**Event buffer.** A continuation step now gives up to 6 parts against the old buffer of 8, beyond which parts were dropped without a
message. The buffer is 32 and a step that fills it is REFUSED ("... the rest would be lost"); `PSS_SALT_EVCAP=<n>` lowers it for a test
(with 6 the 3 nF boost is refused at 100 kHz / 1000 points).

| deck, worst of 20/50/100 kHz vs exact, 200 points, reltol 1e-4 | s144 dB / deg | s146 dB / deg |
|---|---|---|
| bst_cl_csw (boost, 1 nF, chargeless diode) | 0.0055 / 0.098 | 0.0023 / 0.004 |
| bst_q_csw (boost, 1 nF, cj + tt) | 0.0205 / 0.075 | 0.0006 / 0.001 |
| bst_cl_csw100p (boost, 100 pF) | 0.0061 / 0.043 | 0.0002 / 0.001 |
| bst_q_csw3n (boost, 3 nF, cj + tt) | 0.0147 / 0.156 | 0.0021 / 0.005 |
| bb_cl_csw (buck-boost, 1 nF, chargeless) | 0.0210 / 0.093 | 0.0070 / 0.004 |
| bb_q_csw (buck-boost, 1 nF, cj + tt) | 0.0061 / 0.019 | 0.0034 / 0.002 |

1000 points: worst 0.0058 dB / 0.008 deg; 4000 points 0.0045 / 0.014; reltol 1e-6 0.0071 / 0.015; 125 kHz within 0.007 dB / 0.005 deg
(s144 up to 0.19 deg). The two buck-boost 20 kHz pins of `tests/test_pac_charge.py` are now the exact values (the perturbed-PSS
references were 0.105 / 0.108 deg high). Moved inside their bounds, not explained: od_csw10p_nosnub at reltol 1e-6, 100 kHz -0.007 ->
-0.013 dB; bst_cj1n 100 kHz +0.03 -> +0.06 deg (no exact reference); bst_cl_csw100p 100 kHz moved from its perturbed-PSS pin (-0.004 ->
-0.010 dB) while the pin itself is 0.010 dB off exact.

`tests/test_pac_exact.py` pins the 18 exact values and checks the six decks at 200 and 1000 points (0.01 dB / 0.05 deg), the red proof
`PSS_SALT_NO728=1` (restores s144 exactly), the overflow refusal, and re-runs the tool on two points (PAC_EXACT_PYTHON=<python with
scipy>). `PSS_728_OFF=<bits>` turns single parts off (1 share, 2 correction, 4 impulse drift, 8 clamped rows with the lump, 16 lump,
32 pre-event voltage, 128 correction always into the node, 256 delay terms).
Not done: the residual step-related +0.006-0.007 dB of bb_cl_csw at 200 points (0.002 dB at 4000); why the drift-free impulse is wrong
on charge-only bucks; no exact reference for the buck and charge-only decks (they stay validated against the perturbed PSS); no run on the second workstation.

## Non-quasi-static channels in td .pnoise (card #715, 2026-10-19)

Up to s146 td `.pnoise` refused BSIM4 `trnqsmod=1` and BSIM3v3.3 `nqsmod=1` ("not validated", card #559 kept the refusal for `.pnoise`
only). Measured here and lifted: source branch w715 (on s146 = 5d296208c), test build `aaba50ba3f49-w715`; test
`tests/test_pnoise_nqs.py` (9/9 on w715), fixtures `tests/fixtures/pnoise715/`, ensembles `tools/pnoise_mc_bsim4.py lo ...`; runs and logs in
`~/Work/data/ngspice_augmentation/715/`.

**What NQS adds to noise.** Neither model has a noise source on the relaxation (charge-deficit) node: `BSIM4noise` and `BSIM3noise` inject at
the drain/source (and gate/bulk) node pairs, with densities from `qinv`, `ueff`, `gm`, ... of the quasi-static bias (b4ld.c line 3023,
b3ld.c `qinv = -(qgate + qbulk)`); the NQS block only adds the `qcheq` state and the qNode stamps. So NQS changes the TRANSFER from the
generators to the output, i.e. the qNode as a state row of the transposed td recursion, which #559 already carries (`.pxf` adjoint ==
forward to 5e-13 on an NQS deck). `.noise`/`.ac` differ between the two models: BSIM4's AC load isolates the charge node with trnqsmod
(b4acld.c, "not doing a real job"), so its `.noise` is quasi-static (`.noise` with trnqsmod 1 == trnqsmod 0 to all digits); BSIM3's
`b3acld.c` stamps the linearised transient NQS model (gtau, j omega on the charge node) when `nqsmod=1, acnqsmod=0`, so its `.noise` is an
independent reference that contains the relaxation state.

**Static (LTI) stages** (`static_nqs_bsim3.cir`, `static_nqs_bsim4.cir`: W/L 40u/20u, relaxation time ~0.1 us, Vin 0.9 V through Rs 1k,
Rd 20k, Cl 50f; PSS period from a 100 kHz clock on an unrelated net; `.pnoise dev=1 pacpoints=40000`, 0.1 .. 10 MHz):

| case | worst, relative to the total | the quasi-static deck / reference instead |
|---|---|---|
| BSIM3 nqsmod 1: td `.pnoise` vs `.noise`, total and every generator (M1 id/rd/rs/1overf, Rs, Rd) | 2.3e-6 (Rs) | nqsmod 0 `.pnoise` vs the nqsmod 1 `.noise`: 0.137 (Rs), total 3.719e-8 vs 3.564e-8 V/rtHz at 10 MHz |
| BSIM4 trnqsmod 1: td `.pnoise` per generator vs its `.noise` value x \|td `.pac`\| / \|`.ac`\| from the generator's node pair | 6.4e-8 | the quasi-static `.noise` itself: 0.134 (Rs) |

The BSIM4 row checks the adjoint noise engine against the forward td `.pac` (validated against a transient FRA by #559), not the transfer
itself; at 25 .. 100 MHz both models leave the 2e-5 class through the time step (40000 points per 10 us; the trnqsmod 0 / nqsmod 0 decks
do the same against `.noise`: 7e-4 at 100 MHz).

**LPTV stage vs transient-noise ensembles** (`lo_nqs_bsim4.cir`: the #559 LO deck with Rs 100 -> 1k: gate 0.9 V +- 0.35 V at 1 MHz,
the device goes from near cut-off into triode; bands 4.3 .. 42.7 MHz, omega tau 3..30). The #708 method on a time base divided by 1000
(`CIRCUITS["lo"]["scale"]`): boxcar 0.2 ns, 16 seeds x 20 us, 2 us dropped, Hann/Welch 8192, 320 segments, 9-bin bands of 610.35 kHz,
TS = 0.05 ns, K = 100 (res: 30, below), Newton tolerance reltol 1e-7 / vntol 1e-10 (below); the RC calibration at this time base (10k / 1p vs the analytic 4kTR/(1+(f/fc)^2)): worst band 0.172 dB, mean
-0.030 dB. *res*: BSIM4 channel noise off, Rs and Rd as white trnoise currents; *mos*: resistors noiseless, the channel noise as
sqrt(S(t)) w(t) with S(t) from BSIM4's own `.noise` at the orbit's terminal voltages (2002 points, orbit periodicity 1.8e-7 V). Through
the NQS stage the Rs noise is 0.2 .. 3.9 dB lower than with trnqsmod 0; the channel noise moves only 0.01 .. 0.1 dB with NQS (it is injected at the drain).

`.pnoise` minus ensemble, dB in power, per band centre in MHz:

| case | 4.27 | 9.77 | 15.87 | 21.97 | 28.08 | 34.79 | 42.72 | worst | mean |
|---|---|---|---|---|---|---|---|---|---|
| lo res (Rs + Rd), trnqsmod 1, 16 seeds, K = 30 (pinned) | -0.148 | +0.058 | +0.022 | +0.032 | +0.056 | +0.040 | -0.046 | 0.148 | +0.002 |
| lo res, same seeds, K = 100 | -0.210 | +0.062 | +0.023 | +0.035 | +0.062 | +0.044 | -0.045 | 0.210 | -0.004 |
| lo mos (channel), trnqsmod 1, seeds 1..16 | -0.119 | +0.165 | +0.213 | -0.001 | +0.106 | +0.067 | +0.152 | 0.213 | +0.083 |
| lo mos, seeds 101..116 | +0.157 | +0.119 | -0.197 | +0.085 | +0.085 | +0.155 | -0.176 | 0.197 | +0.033 |
| lo mos, both (32 seeds, pinned) | +0.017 | +0.142 | +0.003 | +0.042 | +0.096 | +0.111 | -0.015 | 0.142 | +0.057 |
| red: trnqsmod 0 deck vs the lo res ensemble | +0.102 | +1.019 | +1.874 | +2.620 | +3.188 | +3.598 | +3.860 | 3.860 | +2.323 |
| trnqsmod 0 deck vs its own ensemble (reltol 1e-3, K = 100) | -0.149 | +0.061 | -0.029 | +0.069 | +0.091 | +0.189 | -0.093 | 0.189 | +0.020 |

Band 1-sigma 0.113 dB (16 seeds), 0.080 dB (32); every band of the pinned NQS cases within 1.8 sigma, means within 1.9 sigma of zero
(res +0.05, mos 32 seeds +1.9). The mos ensemble at K = 30 (same seeds) repeats K = 100 to 0.002 dB; the res ensemble does not in its
lowest band (next paragraph). Absolute levels at 4.27 MHz:
res 2.136e-8, mos 2.615e-8 V/rtHz. `.pnoise` cost: 0.5 s for the 72 bins.

**The ensemble needs a tight Newton tolerance with the NQS model.** At the default reltol 1e-3 the noisy-minus-noiseless difference of
the two circuit copies carries a low-frequency component that is not noise: it depends on the tolerance and on K, and the quasi-static
model does not show it. Ensemble minus `.pnoise`, lo res:

| ensemble setting | bin 0.61 MHz | bin 1.22 MHz | bin 2.44 MHz | band 4.27 MHz | bands 9.8 .. 42.7 MHz |
|---|---|---|---|---|---|
| reltol 1e-3, K = 100 | +21.9 dB | +21.3 | +17.1 | +13.60 | within 0.09 dB |
| reltol 1e-3, K = 30 / K = 10 | | | | +8.48 / +3.12 | within 0.07 dB |
| reltol 1e-6, K = 100 | +2.03 | +1.18 | +0.36 | +0.259 | within 0.06 dB |
| reltol 1e-7, vntol 1e-10, K = 100 | +1.55 | +0.79 | +0.23 | +0.210 | within 0.06 dB |
| reltol 1e-8, vntol 1e-12, K = 100 | +1.50 | +0.75 | +0.22 | +0.207 | within 0.06 dB |
| reltol 1e-7, vntol 1e-10, K = 30 (pinned) | +0.37 | +0.20 | +0.12 | +0.148 | within 0.06 dB |

Two effects. (1) The tolerance: the part that goes away between reltol 1e-3 and 1e-7 (1e-8 changes nothing more: bins above 3 MHz
agree to 0.02 %). Not traced further; a likely cause is the charge-deficit node's Newton tolerance (its "voltage" is a charge times 1e9, so
vntol is a large relative error on it), which makes the two copies' Newton errors differ and accumulate in the slow relaxation near
cut-off. It is a property of transient simulation with trnqsmod, not of `.pnoise`, and matters for anyone running transient noise on NQS
devices. (2) What is left at a tight tolerance below ~2.5 MHz scales with K^2 in the K-referred power (0.61 MHz bin: excess power 0.43 at
K = 100, 0.09 at K = 30; a K-independent part of ~0.05, i.e. the bin's own 1-sigma, remains): a second-order (rectified) part of the large-
signal ensemble, which the linear `.pnoise` does not contain by construction. It reaches the 4.3 MHz band (bins 1.8 .. 6.7 MHz) through
its lowest bins, hence the pinned res ensemble is the K = 30 run; the bands from 9.8 MHz are the same at K = 30 and 100 to 0.006 dB. The
quasi-static deck shows neither effect at reltol 1e-3 and K = 100 (bins 0.6 .. 3.7 MHz within +-0.17 dB of its `.pnoise`).

**Change** (pacshoot.c `pac_find_nqs_device()`): in `.pnoise` (strict) BSIM4 trnqsmod and BSIM3 nqsmod are accepted like in
`.pac/.pstb/.pxf`, with the note "... (card #559; .pnoise: no noise source on it, card #715)"; an NQS model of any other family is still
refused (BSIM3v0/v1/v32, BSIM4v5/v6/v7: the td noise engine refuses their noise routines anyway). `PAC_NQS_PNOISE_REFUSE=1` restores the
.pnoise refusal (red proof), `PAC_NQS_REFUSE=1` the pre-#559 refusal in every td analysis; the refusal text names the variable that caused
it. The note buffer grew 160 -> 256 bytes (the longer note was cut). No other code path changes: decks without an NQS model run the same
code.

**Test** `tests/test_pnoise_nqs.py` (exit 0/1/2; outputs `~/Work/data/ngspice_augmentation/715/test` or `$PNOISE_NQS_DATA`; default binary
s147, `NGSPICE` overrides): static bsim3 (2e-5 of the total), static bsim4 (2e-5) plus its red proof (the quasi-static `.noise` must differ
by > 1e-2), lo res and lo mos against the pinned ensembles (0.5 dB per band, 0.2 dB mean: 4.4 / 4.7 sigma of 16 seeds), red proofs nqsmod 0
in the static bsim3 check and trnqsmod 0 in the lo res check (both must FAIL the check, and do: 0.137 and +2.3 dB mean), and both refusal
switches. w715: 9/9, 10 s. s146: 2/9 (every NQS `.pnoise` refused: 5 COULD-NOT-RUN; the two refusal checks FAIL because the s146 text
does not name the variable). `tools/pnoise_mc_bsim4.py`: circuit `lo` (fixture folder and time scale per circuit), `MC_BSIM4_DATA`,
`MC_PARAMS` (.param overrides), `MC_OPTIONS` (ensemble .options), `rc` calibration at a circuit's time scale, and an orbit start point
interpolated when the transient has no point exactly at the settle time; the #708 invocations are unchanged.

**No regression** (w715 vs the s146 logs in `port-s146/`, copies of the suites with every output redirected under 715/suites, fixtures
copied with their hard-wired output paths rewritten): `test_pnoise.py` 11/11, `test_pnoise_bsim4.py` 19/19, `test_pac.py` 35/35, every
reported number identical to s146.

**Needs a follow-up outside this card:** `tests/test_pac_charge.py` `check_nqs_refusals()` asserts "`.pnoise` on an NQS model stays refused";
on w715 that deck now runs (PNOISE-TD: done), so the check FAILs until it either expects acceptance or sets `PAC_NQS_PNOISE_REFUSE=1`.

Not done:
- BSIM3 nqsmod in an LPTV circuit against an ensemble (the static BSIM3 check is exact and independent; the LPTV aspect is validated on
  BSIM4, with the same recursion);
- tnoimod 1/2 and flicker with NQS (only tnoimod 0, flicker off, in the lo deck; the static decks have flicker off);
- NQS devices with ideal-switch events in `.pnoise` (refused anyway: `.pnoise` with saltation events), free-running orbits;
- an NQS stage whose own channel noise is strongly NQS-sensitive (here 0.1 dB): the channel-noise path through the relaxation state is
  checked, but weakly discriminated; the Rs path carries the discrimination (3.9 dB);
- instance-level `nqsmod`/`trnqsmod` (given on the device line, not the model) is not seen by `pac_find_nqs_device()` (model parameters
  only), in `.pac` as before and now in `.pnoise`: accepted silently for BSIM3/BSIM4 (validated), and for other families not caught (a
  pre-existing gap of #559, not changed);
- the cause of the ensemble's low-frequency tolerance artefact; no ASan run or run on the second workstation; test_pac_charge.py not run.


## Transmission lines in the shooting PSS and the td small-signal analyses (card #716, 2026-10-22)

**Before** (s145, measured for this card; decks in `tests/fixtures/tline716/`, data in `~/Work/data/ngspice_augmentation/716/`): the
lossless line (TRA, the `T` element) restarted every shooting period from the constant history of its initial conditions (IC=, default
0). `.pss method=shoot` reported CONVERGED with a wrong orbit; td `.pac/.pstb/.pxf` refused the circuit ("delay history that is not a
charge state", card #542).

| deck (500 kHz) | line | s145 .pss orbit vs settled .tran (max, of the peak-peak) |
|---|---|---|
| tra_rc: pulse, Rs 20, line, 100 ohm \|\| 1 nF | Z0 50, TD 0.6 us (0.3 T) | v(b) 22.9 %, v(a) 29.0 % |
| tra_rc, TD 3 us (1.5 T) | | v(b) 100 %, v(a) 58 % |
| tra_rc, Rs 5, 1 k \|\| 1 nF | TD 0.6 us | v(b) 52.9 %, v(a) 20.5 % |
| tra_lcp: pulse 12 V, LC filter, line, 3.3 ohm \|\| 100 nF | Z0 5, TD 0.3 us | v(vl) 4.22 V (79 x the 54 mV ripple) |
| tra_buckp: ideal-switch buck (fixed duty), same output line | Z0 5, TD 0.3 us | v(vl) 3.94 V (70 x the ripple) |

**Change.** `tra/traload.c`: the three-point interpolation of the delayed excitation is now `TRAinterp()` (the former MODEINITPRED code,
unchanged) and records the points and weights of the current step (`TRAipT/TRAipF`); a run started while `TRApssRestore` is set takes the
delay table from `TRApssHist` (times <= 0) at MODEINITTRAN and sets again the breakpoints `TRAaccept()` set when those points were accepted
(`TRApssBreaks`). `pssshoot.c`:
- `pss_tl_latch()`: the table at the end of every latching run (warm-up, Newton baselines), shifted to the next run's start. Within a
  periodic orbit the waves over [-TD, 0] are those of [T - TD, T] (TD > T: the periodic extension, kept by the shifted table), so this
  lagged history is right at the fixed point and carries the warm-up and the native Newton near it.
- `pss_tl_newton()`: the lagged history alone does not settle where the line is lightly damped (tra_lcp: a 5 ohm line between a 20 mOhm
  ESR and 100 nF keeps its high-frequency waves for tens of periods; every native Newton step launches a new mismatch; the history moved by
  10-80 % of its peak-peak from round to round for hundreds of rounds, the orbit stayed 58 % of the ripple off). The history becomes a
  Newton unknown: both excitations of every line sampled at fixed times tau_j over [-TD - 2 dt, 0], dt = T/Points (PSS_TL_OVERSAMPLE
  divides it); the restored table IS the samples; Phi(x, h) = (native states at T, the run's own table at tau_j + T). Dense FD Newton
  on one discrete map at a time: a free run records the grid, F, the columns and the trials REPLAY it (PSSreplayActive; with an ideal switch
  the recorded times are breakpoints, as in the engine's own Newton). On tra_lcp (linear) one Newton step reaches 1e-8; breakpoints alone
  stalled at 0.2, re-recording without replay wandered at 1e-3. Converged on a grid (1e-7, or stalled below 1e-4), a new grid is recorded
  at the solution: its F is the "grid check" (how far the discretisation moves the fixed point; PSS_TL_TOL_X / _H 1e-2). History FD step
  1e-3 of max(|h|, peak-peak) (the excitation is the size of the port voltage, whose solve noise swamps a smaller step). The joint
  Floquet spectral radius is printed (the CONVERGED line's is the native states' with the history frozen: tra_lcp native 0.936, joint
  e.g. 0.888 on tra_buckc).
- LTRA other than R-G (convolution over the whole past on the transient's time-point list), TXL (`TransLine`) and CplLines are REFUSED in
  the shooting PSS (before: silently wrong); a TRA in a free-running (trigger=) analysis is REFUSED. `PSS_SHOOT_JACOBIAN=analytic` is
  ignored with a line (the sensitivity recursion has no delay terms).

`pacshoot.c` (`pac_tl_*`): the linearised line is the line: the excitation of port 1 at step k is the column's own wave y(p2) - y(n2) +
Z0 y(br2) at the large-signal step's interpolation points with its weights (TRAip*), entering the branch-equation right-hand side like
TRAinput (native, generic and gear recursions); the waves of every column are stored per accepted step. A wave at t <= 0 is a history
sample of the column; y(tau) = e^{-jwT} y(tau + T), read at the end of the run with the large-signal table's own quadratic interpolation
(PSStlTableAt, the operator of the PSS map). The Floquet problem is solved on (native states, history samples): pass 1 carries one
homogeneous column per sample, (e^{jwT} I - Phi_aug) xi = z_aug, pass 2 starts every frequency from its xi (both parts). `.pstb`'s Phi
reuse is skipped with a line. `.pnoise` with a line is REFUSED (no transpose of the delay in the adjoint sweep); `.pxf adjoint=1` falls back
to the forward engine. `PAC_TL_NO_DELAY=1` (red proof) drops the delayed terms; `PSS_TL_OFF=1` restores the whole old behaviour.

**Measured** (w716 = s146 candidate + this card, commits e18c59508 c8119851b; `tests/test_tline_pss.py` 21/21 in 141 s, s145
2/21 (5 orbit FAIL, 4 refusal FAIL, every td analysis COULD-NOT-RUN: refused); Points 400):

| check | result |
|---|---|
| tra_rc TD 0.3 T orbit vs .tran (reltol 1e-8) | 0.37 % of p-p at v(a) (timing-tolerant T/1000; plain 1.33 % at an edge), v(b) 0.03 % |
| tra_rc TD 1.5 T | 1.22 % (a history corner, O(dt^2): 3.4 / 1.2 / 0.24 % at Points 200 / 400 / 800; PSS_TL_OVERSAMPLE=2: 0.55 %) |
| tra_rc Rs 5, 1 k | 0.28 % |
| tra_lcp (high-Q line) | 0.026 % (14 uV); grid check 1.5e-3 / 3.7e-3 |
| tra_buckp (ideal switch) | 0.80 % of the 56 mV ripple (0.46 mV); the transient reference itself moves 0.97 mV (1.7 %) between max step T/1000 and T/10000 at reltol 1e-8; without the line (buckp_noline) 0.20 % |
| td .pac vs .ac, tra_rc TD 0.3 T / 1.5 T / tra_lcp, 20-300 kHz | <= 3.7e-6 / 1.7e-6 / 1.3e-4 dB, <= 1.1e-4 / 1.6e-5 / 3.2e-4 deg; generic recursion (PAC_FORCE_GENERIC=1) identical |
| td .pxf vs .ac | <= 1.9e-6 / 3.3e-7 relative |
| .pstb vs .stb (tl_pstb, line TD 20 ns in the loop), 1 kHz - 1 MHz | 6.3e-5 relative at 1 MHz (no line: 1.1e-4); above, both grow as the engine's O(h^2) |
| tra_buckp line-to-output (v(vl)/Vin), 20 / 50 / 100 kHz vs perturbed PSS | <= 0.0001 dB / 0.003 deg (reference halving <= 1e-4 dB); a transient-backend reference (independent of the large-signal line history, amplitude 0.3 V, reltol 1e-7) agrees at 20 / 50 kHz within its own uncertainty (0.006 dB / 0.10 deg, 0.016 dB / 0.13 deg; 100 kHz refused by its halving check) |
| tra_buckc control-to-output (ramp comparator, saltation), 20 / 50 / 100 kHz | <= 0.0010 dB / 0.016 deg (reference halving <= 0.0008 dB / 0.012 deg) |
| red proofs | PSS_TL_OFF=1: orbit 27 % off, .pac refused; PAC_TL_NO_DELAY=1: .pac output 0 (398 dB off) |

Cost: tra_rc 0.5 s (248 samples, 506 runs), TD 1.5 T 1.9 s (1206 samples), tra_lcp 1.5 s, tra_buckp 13-15 s (128 samples, ~1900 runs; no
replay with an ideal switch, the Newton converges linearly), against 0.1 s for the wrong s145 orbit; td .pac per frequency point 0.4-19 s.

**No regression** (copies of the suites under 716/regress, outputs and fixtures' hard-wired paths redirected; vs `port-s146/`):
`test_pac.py` 35/35, every line identical; `test_pss.py` 39/44 + 5 info: every value identical to s146; differences: the copy runs
the determinism checks at 2 threads instead of 8 (memory rule; pass), check_regressions() (which runs other suites) is disabled, and
r3_longsettle's speed ratio read 9.92x (info) vs 10.91x; interleaved A/B on that fixture, 8 runs each: s146 median 0.383 s, w716 0.374 s
(no slowdown: machine load during that run).

Not done:
- LTRA (lossy), TXL and CPL lines: refused, not supported (their history is a convolution since t = 0; would need a periodic
  convolution); LTRA R-G is accepted (a deck converges; its orbit was not compared);
- TRA in `.pnoise` (refused) and in the adjoint `.pxf` (forward fallback); TRA with stateful XSPICE models in td small-signal (refused);
  free-running orbits with a line (refused);
- the joint Newton is dense: (n + 2 Ns) runs per iteration and an O(N^3) solve; a long line against T/Points (PSS_TL_MAXSAMPLES 4000) is
  refused; no Krylov variant;
- the history representation (samples at T/Points) limits the orbit at corners of the waves (O(dt^2), measured above);
- switching decks: the pinned reference is the perturbed PSS of this same build (it shares the large-signal line history, not the
  small-signal code); the transient-backend reference is independent but noisy (above); no step-refinement self-check (1600 points per
  cycle on the 25-period LCM orbits: hours); one control-to-output and one line-to-output deck; no closed-loop converter with a line;
- multi-line decks, coupled lines, a line inside a loop with an ideal switch: not run.

## Small signal through digital XSPICE blocks: td `.pac` / `.pstb` / `.pxf` (card #714, 2026-10-23)

Up to s149 every event-driven XSPICE model (adc_bridge, dac_bridge, the d_ gates, flip-flops, latches, d_osc, d_pwm, ...) made td
`.pac/.pstb/.pxf` refuse ("event-driven / digital", card #551). With the refusal lifted (`PSS_XSPICE_ANY=1`) the comparator-PWM buck
below gave a control-to-output of exactly 0 (the quantiser has no gain at any accepted step) and the d_pwm decks crashed the recording
pass (SEGV: the recording runs did not put the carried event-system state back). And the closed comparator loop had no orbit:
`.pss` failed after 150 iterations (residual 4e-3), or, with the d_osc-reset ramp, "converged" with a spectral radius 1.087.

### The linearisation of a digital path

A digital block has no small-signal gain. What moves is the TIME of a digital edge, when an analog quantity sets it:

* **adc_bridge threshold crossing.** The input v crosses the level L at which the bridge's state changes (in_high / in_low; with
  in_low > in_high the Schmitt levels; with a band the 0 -> X and X -> 1 crossings are two edges). The edge moves by
  `dtau = -dv(t_c) / v'(t_c)`: dv(t_c) is the small-signal input at the crossing (linear interpolation of the recursion's columns between
  the accepted points around it), v' the orbit's slope there.
* **Blocks with fixed delays** (gates, buffers, inverters, flip-flops, latches, dividers, d_lut, d_state) pass the shift of the input edge
  that triggered an output change on to that output: the delay is a constant. An output triggered by an edge of fixed time (a clock)
  does not move, whatever its data inputs did (the zero gain of a quantiser: a d_dff clocked by a fixed clock passes no small signal).
  This is followed per event at run time, not from the netlist: every queued output event carries a tag (0: fixed time, else a timing
  record of the analysis), every event node keeps the tag of its last change, and an instance called because inputs changed passes the
  tag of those inputs to its outputs (`evtload.c`, `evtqueue.c`, `evtdeque.c`; `g_xs_dig_*` in `mif.h`).
* **dac_bridge.** An input edge moved by dtau moves the ramp that follows it: from the edge until the output is flat again
  `dout(t) = (s_before - s_after) dtau` (s_before / s_after the orbit's output slope just before / after the edge; from a flat level
  this is minus the ramp slope times the shift; a ramp interrupted by another edge keeps the offset until it reaches a level). It is
  injected on the dac output's branch row (voltage outputs) at every accepted point inside that window, with the left half of the
  point's trapezoid weight, which gives the exact area of the shifted ramp when the window starts and ends on accepted points. Two grid
  details had to be fixed for that: cm_dac_bridge sets its ramp-end breakpoint only when an edge falls inside a step, never when the edge
  is on a time point (the usual case), so the recording run sets one at the end of every moved ramp; and in the periodic analyses the
  breakpoint a model sets at the current event time came back a few 1e-19 s off it (continuing model time, #551) and dctran.c refused
  it ("bp >= evt_step"): the edge then did not land on a time point. `cm.c` now snaps such a breakpoint onto the event time
  (`XS_NO_EVT_SNAP=1`: old behaviour; an edge inside a step is handled too, with the partial weight of its part of the step; measured:
  same loop gain to 1e-8). The analog network then carries the injection like any source: through a smooth switch model, or through an
  ideal sw whose control it moves (the switch's own saltation term reads it at its event).
* **d_pwm with an analog duty.** cm_d_pwm decides each transition at one call, from the control it reads then:
  `when = last_time + T d(c)` (falling) or `+ T (1 - d(c))` (rising). The shift is the shift of last_time plus `+- T d'(c) dc` (d' the
  slope of its PWL table, as clamped at INIT; 0 where the duty is clamped). last_time carries over the period boundary (a transition can
  be in flight across t = T), so its shift is an extra Floquet state of the boundary solve, per d_pwm whose control can carry a small
  signal (not one held by a V source this analysis does not excite), like the transmission-line history of #716. A transition carried in
  from the previous period gets the state's value as its tag (`g_xs_dig_carry_tag`, read by XSsnapPut). The orbit control must be
  constant: a control that moves along the orbit makes the PWM period depend on it.
* **Clocks** (d_osc, d_pwm whose control is a DC source): fixed waveforms; nothing is injected (case (c), no refusal needed).

### Exact crossings: adc_bridge location (`evtsetup.c` XSadcLocate)

The stock adc_bridge flips at the first accepted point past its level: the edge is quantised to the step grid. In a closed loop the
shooting map is then piecewise constant in the edge time and its finite-difference Jacobian is wrong (s149: no convergence, or a false
radius > 1). A converged step in which a single-bit adc_bridge input crossed its level is now rejected and redone to end at the linearly
interpolated crossing (plus half the tolerance), until the overshoot is below `XS_ADC_TTOL` (1e-16 s), like the ideal-switch location of
#498. Active in the shooting and td recording runs (`PSS_ADC_LOCATE=0`: off) and in any transient with `XS_ADC_LOCATE=1` (so a settled
`.tran` or an FRA reference places the edges the same way). Multi-bit adc_bridges (one input, several outputs) are not located.

The td recording runs of circuits with digital blocks take the fine default recording reltol (1e-8, card #625) instead of the inherited
one (`PAC_DIG_COARSE=1`: inherited): at the inherited 1e-5 the comparator buck was 0.006 dB off the FRA, at 1e-7 0.0004 dB.

### Refused, with the instance named (the recursion would be wrong without notice)

* simultaneous input edges of different timing origin at one instance (a tie: which input decides is not differentiable);
* a moved edge reaching the analog side through a multi-bit adc_bridge, an inverted adc output port, an adc input of type vnam, a
  dac_bridge in multi-input mode or with a current / other output;
* a d_osc whose control has a small-signal part (a VCO), a d_pwm whose control moves along the orbit;
* a moved dac ramp or a moved edge (other than a d_pwm transition) still in flight at t = T (shift the period phase);
* `.pnoise` (the transposed sweep does not carry the edge shifts); `.pxf adjoint=1` falls back to the forward engine with a message;
* free-running (`trigger=`) analyses with XSPICE (unchanged, #551).

A closed loop through a d_pwm (compensator output as its duty) stays refused by `.pss` ("does not settle", #712), and that is right: the
settled transient of `tests/fixtures/digclock/buck_cl_dpwm.cir` switches at 2.00876 us, 0.44 % longer than 1/fs, because cm_d_pwm times
each half period from the previous transition with the control it reads then; the steady state is a free-running oscillation, not periodic
with the analysis period.

### Measured (build `aaba50ba3f49-w714`, `tests/test_digital_pac.py` 23/23 in 173 s; decks `tests/fixtures/digital714/`)

Facts on s149 first (`~/Work/data/ngspice_augmentation/714/facts/`): every deck below refused in td small signal; with
`PSS_XSPICE_ANY=1` the comparator buck gave exactly 0 and the d_pwm decks crashed (SEGV); `.pss` of the closed comparator loop with a PULSE
ramp failed (residual 4.1e-3 after 150 iterations); with the d_osc-reset ramp it "converged" at iteration 27 with a spectral radius
1.087 (UNSTABLE) in one run and failed in another; the closed d_pwm loop is refused ("does not settle").

References: brute-force FRA (one sine on a DC source, or on the lgprobe Vi / Ii, base run subtracted, exact integral of the
piecewise-linear waveform over whole periods of f and of 1/fs), plain transient with crossings located (`XS_ADC_LOCATE=1`) on the
adc decks, the unmodified s149 binary on the others. "variants": the largest move of the reference under amplitude halving and tmax
halving (1 -> 0.5 ns); FRA with reltol 1e-7 instead of 1e-5 moved < 1e-5 dB (a_ol). Each td point is a few seconds.

| deck (all 12 V -> buck, 500 kHz, tanh switches unless noted) | f | td vs FRA | FRA variants |
|---|---|---|---|
| a_ol: PULSE ramp vs control -> adc_bridge -> d_buffer -> dac_bridge, control-to-output | 2 / 25 / 50 kHz | +0.0001 / -0.0004 / -0.0002 dB, <= 0.0017 deg | <= 1e-5 dB |
| a_ff: d_osc sets a d_dff, comparator resets it, two d_inverters | 2 / 25 / 50 kHz | +0.0001 / -0.0004 / -0.0002 dB, <= 0.0017 deg | <= 1e-5 dB |
| a_ffsw: a_ff driving IDEAL sw switches (saltation reads the moved dac ramp) | 2 / 25 / 50 kHz | -0.0000 / -0.0009 / -0.0004 dB, <= 0.0044 deg | 0.0015-0.0019 dB (sw edges on the grid in the plain transient) |
| b_ol2: d_pwm with an analog duty (Floquet shift state) | 2 / 25 / 50 kHz | +0.0001 / -0.0005 / -0.0003 dB, <= 0.0014 deg | <= 1e-4 dB |
| c_ol: clock-only digital path (d_pwm fixed duty), line-to-output | 2 / 25 / 50 kHz | +0.0000 / -0.0005 / -0.0003 dB, <= 0.0016 deg | <= 1e-5 dB |
| a_cl: closed loop, type III, PULSE ramp + adc comparator, `.pstb` (card done-when) | 2 / 10 / 20 / 31.25 / 35.7 / 50 kHz | +0.0002 / +0.0002 / +0.0001 / -0.0005 / -0.0009 / -0.0012 dB, <= 0.0002 deg | <= 0.0010 dB |
| a_cl phase margin / crossover | | .pstb 38.350 deg at 33547 Hz, FRA (interpolated) 38.345 deg at 33536 Hz | |
| a_cl2: closed loop, d_osc-reset ramp, adc with a +-1 mV band (two edges per crossing, dac via out_undef) | 10 / 20 / 50 kHz | -0.0004 / -0.0019 / -0.0042 dB, <= 0.0033 deg | 0.0030 / 0.0022 / 0.0028 dB (2 kHz: 0.07 dB, not used) |
| c_cl: closed loop, clock-only digital path, analog tanh comparator | 10 / 20 / 50 kHz | -0.0000 / -0.0001 / -0.0037 dB, <= 0.0025 deg | 0.0015 / 0.0006 / 0.0029 dB (2 kHz: 0.017 dB, not used) |

Orbits (`.pss` vs a settled transient with located crossings, % of the vout ripple): a_cl 3 iterations, rho 0.980, vout 0.000 %
(0.01 uV), opout 0.001 %; open-loop decks 0.012-0.042 % (the 3 ms transient's own settling). a_cl2: vout 0.011 %, opout 0.45 % of its
ripple (an average offset of 55 uV, 1e-4 relative, the same with a 14 ms transient): the #551 lagged hidden-state iteration
(hidden-state drift 5e-5, the dac_bridge ramp state of the clock path), not this card.

Red proofs (each breaks its check): `PAC_DIG_REFUSE=1` (the s149 refusal), `PAC_DIG_NOSHIFT=1` (a_ol control-to-output 0, a_cl loop
gain -328 dB: the loop is open), `PAC_DIG_PWM_NOSTATE=1` (b_ol2 27.67 dB instead of 21.65), `PSS_ADC_LOCATE=0` (a_cl `.pss` fails after
150 iterations). Without the event-time snap and without the edge-inside-step weight (first build) the moved ramp started one step late
and a_cl's loop gain came out 2.02 times the FRA (s_before read mid-ramp); with only the snap off (`XS_NO_EVT_SNAP=1`) the result is the
same to 1e-8. Without the fine recording reltol (`PAC_DIG_COARSE=1`): +0.006 dB.

### Regression (copies of the suites under `714/regr`, outputs and fixtures' hard-wired paths redirected; vs `port-s149/`)

- `test_pac.py` 32/33: every value identical to s149; the ASan check could not run (no ASan build of w714).
- `test_pac_charge.py` 158/159: every value identical; the ASan check fails because no ASan binary was given.
- `test_xspice.py` 30/32: the refusal check "adc_bridge/d_tff/dac_bridge refused in td .pac" now FAILS by design (accepted; that deck
  has no AC source, the answer is 0); ASan not run. The digital PSS orbits vs the settled PLAIN transient moved from 0.014-0.155 % to
  0.035-0.21 % of the ripple (bound 0.4): the orbit now places the adc_bridge edges at the crossings, the plain transient on its 0.25 ns
  grid. Against a transient with `XS_ADC_LOCATE=1` the same orbits agree to 0.0008-0.0014 % (tff, logic, seq, fdiv_lut, pull, dsample;
  osc_pwm 0.0013 % either way). R9 identical within 0.0001 mV (hidden-state drift 7e-5 -> 2e-7).
- `test_pss_digclock.py` 8/9: check 6 ("td .pac through the digital blocks keeps refusing") now FAILS by design: that deck (d_pwm clock at
  fixed duty, AC on Vin) is case (c) and runs; check 3 (clock edge on the period boundary) 0.0486 -> 0.0314 % of the ripple (event-time
  snap). Both suites' refusal checks need updating to the new behaviour (not done here: those files are outside this card's scope).

### Not done / limits

- d_pwm in a closed loop (free-running, see above): refused; no free-running PSS with XSPICE.
- d_osc as a VCO, multi-bit adc_bridge codes reaching the analog side, dac_bridge multi-input mode and current outputs, inverted adc
  ports, vnam inputs: refused (not linearised). A multi-bit adc whose code is only sampled by clocked flip-flops passes (its edges are
  then fixed: zero gain), not run on R9.
- `.pnoise` with digital blocks (refused), adjoint `.pxf` (forward fallback).
- A moved edge in flight at t = T (other than a d_pwm transition): refused; the shift would have to become a Floquet state like the d_pwm one.
- Ties are refused even when the two edges have the same shift (two identical comparators).
- The adc crossing of the very last step of a period (event processed after the run ends) is not caught by the in-flight check.
- The dac injection window is exact only when its ends are accepted points (a ramp interrupted before reaching a level is weighted from
  the step's slope); measured cases all have breakpoints at both ends.
- Not run: R9 (multi-bit DPWM + digital compensator), BSIM4 + digital decks, Krylov PSS on these decks, ASan build of w714.
- The adc location changes the PSS grid of every deck with a single-bit adc_bridge (`PSS_ADC_LOCATE=0` restores the old grid); the
  event-time snap changes the grid of periodic runs with dac_bridge edges.
