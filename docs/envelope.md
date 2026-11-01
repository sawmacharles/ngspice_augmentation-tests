# Envelope following for long switching-converter transients (card #512, round 10)

Card #512 (transient acceleration round 10), 2026-09-29. Goal: accelerate
soft-start / slow load-line-ramp transients that span thousands of
switching periods, the way round 9's `tools/tran_from_pss.py`
(docs/tran-from-pss.md) accelerated the OTHER long-transient shape
(power-up settle before a single event) by starting from `.pss
method=shoot`'s converged orbit instead of simulating the settle. This
document covers the mechanism built (`tools/envelope_tran.py` +
`tools/ckpt_io.py`), what it measures, and -- stated as plainly as
docs/tran-from-pss.md's own "Why not (a)" section -- where it falls short
of the card's own accuracy bar and why, with a characterised root cause
and a concrete next step.

## Mechanism chosen: (a), a Python driver over savestate/loadstate

The card named two shapes: (a) a Python driver using `savestate`/
`loadstate` + short `.tran` segments, needing a way to WRITE a modified
checkpoint state (`loadstate` only restores exactly, as round 9 already
found); or (b) a new C engine option (`.tran ... envelope=<M>` or a new
`.envtran`) doing the same thing natively.

**(a) was chosen.** Reasoning, in the same shape round 9's own mechanism
choice used:

- The checkpoint file format (`src/frontend/com_checkpoint.c`,
  `com_savestate`/`com_loadstate`) is a small, fully characterised
  sequence of raw `fwrite`/`fread` calls -- read directly from the C
  source, not guessed -- so writing a MODIFIED checkpoint is tractable in
  pure Python with zero new C. `tools/ckpt_io.py` implements this and is
  verified (not assumed) to round-trip an unmodified file byte-for-byte
  identical to what `savestate` itself wrote (`tests/test_envelope.py`'s
  `test_roundtrip_byte_identical`), and a patched file that only changes
  `CKTtime`/`CKTrhsOld`/`CKTorder`/`CKTdelta`/`CKTbreaks` `loadstate`s
  cleanly and continues the simulation (checked by hand first, then
  automated -- see "The checkpoint file format" below).
- This project's build/test machine, the development workstation (AMD Ryzen AI Max+ 395, Linux), was, for this whole round, a
  genuinely tight-RAM shared box: `free -h` showed 2.4-2.5 GiB free and
  swap at 7.9/8.0 GiB used throughout (other users' builds/tests running
  concurrently), and `uptime` load average 1.0-6.3. A new C engine option
  needs a fresh `buildrecipe/build.sh` pass (8-way parallel `make`) on top
  of that load -- risk round 9 already weighed the same way for its own,
  simpler single-hop handoff and decided against; round 10's repeated,
  many-hop driver has MORE moving parts to validate, not fewer, so the
  same reasoning applies at least as strongly. `ngspice-src-env` (the
  worktree card #512 sets up for this) was consequently NOT built this
  round -- `aug-env` is untouched, no C commit exists for this card.
- Mechanism (a) is directly testable per-hop (a single Python process
  drives many small, independently-inspectable ngspice subprocess calls),
  which is how the accuracy shortfall below was actually FOUND -- an
  isolating control run (M=1, no jumps at all, pure checkpoint chaining)
  separated "checkpoint format/chaining error" from "projective-jump
  error" cleanly, in a few lines of Python. That kind of ablation is much
  harder to get inside a single C engine pass.

## The checkpoint file format (verified, not assumed)

Read directly from `ngspice-src-env`'s (== `aug`'s, patches 0032-0033,
unchanged) `src/frontend/com_checkpoint.c` and `src/include/ngspice/
cktdefs.h` for exact C types (found the hard way: `CKTmode` is `long`,
8 bytes on this build's x86-64 Linux -- not `int` like its neighbours;
guessing `int` there would have silently misaligned every field after
it). A sequence of raw, unpadded `fwrite`/`fread` calls:

```
magic            8 bytes literal b"NGCKPT01"
maxEqNum, matSize, numStates, maxOrder, breakSize    5 x int32 (header)
CKTorder                                              int32
CKTmode                                                int64
CKTtime, CKTdelta, CKTsaveDelta, CKTdelmin, CKTminBreak   5 x float64
CKTdeltaOld[7]                                        7 x float64
CKTbreaks[]                                           breakSize x float64 (if breakSize>0)
rhsOld_present (int32) + CKTrhsOld[]                  (matSize+1) x float64 (if present)
irhsOld_present (int32) + CKTirhsOld[]                (matSize+1) x float64 (if present)
8x: state_present (int32) + CKTstates[i][]            numStates x float64 (if present)
```

`tools/ckpt_io.py` implements `read_checkpoint()`/`write_checkpoint()`
against exactly this layout. Verified directly (`tests/test_envelope.py`):
`savestate` a real deck's mid-run state, `read_checkpoint()` it,
`write_checkpoint()` it straight back out, `filecmp` the two files --
**0 bytes different**. A further check confirms `loadstate`-ing the
round-tripped file produces **byte-identical continuation output** to
`loadstate`-ing the original (both compared against the same further
`.tran`).

**A trap found and worked around, not in the format itself**: ngspice's
control-block tokenizer LOWERCASES bareword file arguments (confirmed
directly: `savestate MixedCase.bin` writes `mixedcase.bin`; quoting does
NOT help -- `savestate "MixedCase.bin"` writes a file literally named
`"mixedcase.bin"`, quote characters included). This project's own data
root is `~/Work/data/...` (capital W), so embedding an ABSOLUTE checkpoint
path in a control block silently reads/writes the WRONG, lowercased,
usually-nonexistent file. Fixed the same way `tools/tran_from_pss.py`/
`tests/test_checkpoint.py` already avoid it (their own decks never
reference an absolute path either): every checkpoint filename in the
control text is a bare relative name, with the ngspice subprocess's `cwd`
set to the run's own output directory (`tools/envelope_tran.py`'s
`run_window()`); the absolute, possibly-mixed-case path is used only on
the Python side.

## The algorithm

Sample the circuit's full Newton-unknown vector (`CKTrhsOld`) once per
switching period -- the same "stroboscopic map" this project's own `.pss
method=shoot` (`pssshoot.c`) already uses for a driven periodic orbit,
and exactly what Petzold 1981's "quasi-envelope" and Kundert/White/
Sangiovanni-Vincentelli 1988's envelope-following both sample. Per
`research/transient-acceleration.md`'s round 10 section for the literature
detail; short version of what's implemented here:

1. **Bootstrap**: 2 real periods from t=0 (`uic`) give x0, x1.
2. **Smoothness gate** (before attempting a jump): compare `||x1-x0||`
   against the MEDIAN of the last 8 ACCEPTED per-period deltas -- not a
   fixed fraction of `||x1||`. This circuit's state vector always
   contains `Vin` (a constant 12V rail) at a fixed index, so `||x1||`
   never shrinks even when vout is still near 0V right at the start of
   soft-start; a naive relative-to-`||x1||` threshold is far too loose
   exactly when it matters most (checked directly: replacing the
   history-relative test with a fixed-relative one let early-ramp
   compensator-startup transients through undetected in a first version
   of this tool). `--growth-tol` (default 5x) is the allowed multiple of
   that recent median.
3. **Explicit projective jump** (Gear/Kevrekidis' "projective
   integration", forward-Euler in cycle-index n; Kundert's own algorithm
   is the IMPLICIT/Newton-corrected version of the same idea -- see "What
   was not built" below): `x_pred = x1 + M*(x1-x0)`. `ckpt_io.py` patches
   ONLY `CKTtime` (so slow-time sources -- the soft-start ramp itself --
   evaluate correctly at the landing time) and `CKTrhsOld` (-> x_pred).
   `CKTorder` is forced to 1 and a single breakpoint inserted at the
   landing time -- the same "state discontinuity -> low-order restart"
   mechanism `docs/checkpoint.md`'s patches 0200/0201 already
   characterise at an ordinary source-edge breakpoint, reused here rather
   than inventing a second one.
4. **Check periods and error control**: `--n-check` (default 4) real,
   fully-resolved periods are simulated from the landing point. The slope
   between the LAST two of those (not the first two -- see "The landing
   kick" below) is compared against the pre-jump slope; disagreement past
   `--slope-tol` (default 3.0 = 300% relative) shrinks M by 4x and
   retries from the same x1. A second, independent scan checks EVERY
   consecutive pair actually simulated inside the check window (not just
   its endpoints) against the same history-relative smoothness test --
   added after testing on the red-proof deck showed a single event that
   fires and re-settles entirely inside one jump's own check window can
   look smooth end-to-end while hiding a real spike in the middle (see
   "Refusal behaviour" below).
5. Every window that is actually integrated (bootstrap, check periods,
   and any passthrough fallback) is real ngspice output; the algorithm
   never simulates a period it does not also validate the smoothness of.

## What gets extrapolated, and what does not

Only `CKTrhsOld` (the full Newton-unknown vector: every node voltage and
branch current) is extrapolated. `CKTstates[]` -- the roughly
per-device-and-per-order integration HISTORY every multistep/trapezoidal
formula's predictor and LTE estimate reads -- is deliberately left
UNTOUCHED, copied verbatim from the pre-jump checkpoint. This was a
scope decision, not an oversight: `CKTstates[]`'s layout is private to
each device model (a capacitor's charge/current history sits at different
offsets than a BSIM4's dozen internal capacitances), undocumented at the
file-format level, and mapping "extrapolate this device's charge state by
M periods' worth of drift" generically, for every device model this
project's reference decks use, was judged out of scope for one round --
flagged explicitly as the clearest next step below.

## The landing kick (measured, root-caused)

Leaving `CKTstates[]` unmodified means that right after a jump, the
compensator's own slow-pole integrator state (this deck's dominant pole
is `R2f*C2f = 27.51k * 2.299n ~= 63us`, i.e. ~31 switching periods at
500kHz -- far longer than the 2-4 check periods originally used to judge
smoothness) is briefly INCONSISTENT with the jumped-forward
`CKTrhsOld`/`CKTtime`. Newton's algebraic reconciliation at the landing
instant produces a measurable, systematic excess in the immediately
following periods' own per-period slope -- not noise, and not something a
looser random-noise tolerance can distinguish from a genuine event.
Measured directly (a debug instrumentation pass, `pre_slope`/`post_slope`
printed per jump attempt, `env_r3d_softstart.cir`, first 200us): the
pre-jump slope stays a stable ~2.4-2.7 mV/period throughout (matching the
expected ~3.3V/1500 periods ~= 2.2 mV/period ramp rate closely), while the
FIRST TWO check periods after landing show `post_slope` consistently
2-6x larger, scaling with M (M=5: ~1.2-1.7x excess; M=20: ~5-6x excess).
Comparing the LAST two of `n_check=4` check periods instead (giving the
kick 1-2 periods to relax) reduces but does not eliminate this -- accepted
jumps still show `slope_err` in the low-percent-to-tens-of-percent range,
not the sub-percent a genuinely smooth continuation would show.

**This directly explains the accuracy shortfall below**: an isolating
control run with `m_periods=1` (the inner jump loop's own `m_try > m_min`
guard means NO jump is ever attempted -- pure checkpoint-chaining, one
real period at a time, no extrapolation at all) passes the 0.5%
cycle-average tolerance on 93% of its periods (median error ratio 0.29 of
tolerance, worst case only 1.4x over). The SAME deck with jumping enabled
(M=3, 5, or 10) fails far more often and far worse (worst-case error
ratios 9-13x tolerance, scaling up with M) -- isolating the error source
to the projective jump's own landing kick, not to `savestate`/`loadstate`
chaining, subprocess-per-window overhead, or the checkpoint format.

## Accuracy: measured, against a full plain `.tran` reference

`tests/fixtures/env_r3d_softstart.cir`: `refs/r3d_vm_buck_deadtime.cir`'s
own power stage/compensator/dead-time network, unchanged, with `Vrefsrc`
replaced by a `PWL(0 0 {tss} {vref} {ttot} {vref})` ramp, `tss=3ms`
(1500 periods at fs=500kHz), `ttot=4ms` (2000 periods total: the ramp plus
500 periods of post-ramp settle). Reference: the same deck's own
standalone `.control` block, a single plain `tran 2n 4e-3 0 uic` run
(`~/Work/data/ngspice_augmentation/envelope/env_r3d_softstart/
env_r3d_softstart_reference.txt`, 2036858 rows, `s55`, 18.34s wall).

Validated with `envelope_follow()`'s own `validate_against_reference()`
(a single streaming pass over the reference file, bucketed against the
envelope run's own real-simulated periods -- never interpolates or
compares a period neither run actually produced): every real-simulated
period's cycle-averaged `v(vout)` and `i(L1)`, and its ripple
(max-min), against `max(0.5% relative, 1mV floor)` and `max(5% relative,
0.1mV floor)` respectively (the card's own stated bounds).

Full 4ms run, every real-simulated period checked (not a slice -- the
final, authoritative sweep; an earlier partial-slice measurement made
while root-causing "the landing kick" above is superseded by this one).
"worst rel err" is the single worst period's `|error|/|reference value|`
(a plain relative error, e.g. 0.055 = 5.5%), NOT normalised by the
tolerance -- read against the stated 0.5%/5% bounds directly:

| config | periods simulated (of 2000) | vout avg: frac within 0.5% | vout avg worst rel err | vout ripple: frac within 5% | iL avg: frac within 0.5% |
|---|---|---|---|---|---|
| M=1 (no jumps -- checkpoint-chaining control) | 2000 (all) | 99.2% | 0.7% | 90.2% | 73.6% |
| M=5 | 1070 | 52.7% | 2.5% | 34.5% | 20.7% |
| M=10 (this tool's shipped default) | 790 | 25.1% | 5.5% | 15.8% | 7.8% |
| M=20 | 630 | 35.9% | 3.0% | 24.1% | 7.1% |

(M=20's better vout-avg fraction than M=10's, despite a larger nominal M,
is real and explained by its own much higher refusal count in this table
-- 560 vs 141, see the speed table below -- the mid-check-window scan
shrinks most M=20 attempts down to something smaller before they are
ever accepted, so the EFFECTIVE jump size actually taken is not simply
"20"; not chased further this round.)

**Verdict, stated plainly: the mechanism as built does NOT meet the
card's own 0.5%/5% accuracy bound at any M>1.** The M=1 control shows the
checkpoint-chaining machinery itself is sound (99.2% pass, worst case only
0.7% relative error, close to the 0.5% target and consistent with
docs/checkpoint.md's own established ~1e-4-relative restart-vs-
uninterrupted-run finding once compounded over 2000 chained hops); every
M>1 configuration's shortfall traces directly to the landing-kick effect
characterised above, which gets WORSE with larger nominal M (the very
thing that buys the speedup) -- `i(L1)`'s own average is hit even harder
than `v(vout)`'s (7-21% pass vs 25-53%), consistent with inductor current
being the state most directly perturbed by a duty-cycle guess that Newton
then has to reconcile against stale charge/flux history at landing. This
is a genuine, characterised negative result on the accuracy axis, not a
hidden one -- see "What a future round should do differently" below for
the concrete, literature-backed fix (an implicit/Newton-corrected
envelope step, which this round did not have time to build).

## Speed: measured, full run, interleaved A/B

`m_periods=10, n_check=4, growth_tol=5.0x, slope_tol=3.0` (this tool's
shipped defaults), full `env_r3d_softstart.cir` (0 to 4ms), `s55`, one
ngspice process at a time (subprocess calls are strictly sequential in
this driver -- never more than 1 concurrent ngspice process, well inside
this card's 2-process limit). **Interleaved A/B, n=3** (per this project's
own "re-measure speed claims" convention -- an earlier, non-interleaved
single-shot measurement made before the mid-check-window scan (see
"Refusal behaviour" below) was added found a misleadingly larger ~2.6x;
that scan's extra retries cost real wall time and the number below
supersedes it):

| rep | uptime before A | A (reference) wall | uptime before B | B (envelope, M=10) wall |
|---|---|---|---|---|
| 1 | load 2.95,3.33,3.36 | 18.35s | load 3.18,3.36,3.37 | 15.91s |
| 2 | load 3.20,3.35,3.36 | 18.28s | load 3.08,3.32,3.35 | 15.05s |
| 3 | load 2.76,3.24,3.32 | 18.12s | load 2.76,3.20,3.31 | 15.21s |

Medians: reference **18.28s**, envelope **15.21s** -- **1.20x speedup**,
well under the card's informal ">=5x" framing. Deterministic across reps
(same 790/2000 periods simulated, 189 jumps, 141 refusals every time --
no run-to-run variation, as expected for a fixed algorithm over a fixed
deck with no randomness). **Root cause of the shortfall, stated
plainly**: 790 real periods out of 2000 is a genuine **60% reduction** in
simulated periods, which on a per-period-cost basis alone would suggest
something closer to a 2.5x speedup -- the actual 1.2x is far LOWER
because (a) every one of those 790 periods, AND every one of the 141
refused/abandoned jump-check attempts' own check periods, pays its OWN
ngspice process-launch and netlist-parse overhead (this driver launches
one ngspice subprocess per single-period hop, following the same
per-window-subprocess pattern `tools/tran_from_pss.py` uses for its own
two-phase handoff), and (b) the mid-check-window smoothness scan added to
fix the red-proof gap below (see "Refusal behaviour") makes MORE jump
attempts fail after already paying for their own check periods, without
ever landing a skip -- a direct, measured accuracy-vs-speed trade both
ways: the same fix that (partially) improves safety also taxes the speed
this round's own earlier, less-safe measurement had claimed. At this
converter's own per-period cost (a 2us window of a 7-state deadtime-buck
circuit is cheap to integrate), the FIXED per-process overhead is a much
larger fraction of each hop's wall time than it would be for a longer or
heavier per-period simulation -- a real, measured, and reasoned-about
limitation, not a hidden one. Chaining more real work per ngspice
invocation (e.g. running all `n_check` check periods in ONE `.tran` call
instead of `n_check` separate loadstate/savestate hops, or batching
several accepted jumps' worth of passthrough periods together) is the
clear next lever, not attempted this round.

## Refusal behaviour, on the red-proof deck

`tests/fixtures/env_r3d_softstart_redproof.cir`: the same soft-start deck
plus `Iload_step`, an abrupt 0->4A current-sink load step (`PWL`, 100ns
rise + 1us hold) at `t=1.5ms`, mid-ramp -- substituting for the card's own
"e.g. a DCM/CCM transition" example, since this deck's low-side switch
(`Bls`, a smooth bidirectional B-source, no diode-emulation/forced-CCM
distinction -- see `refs/r3d_vm_buck_deadtime.cir`'s own header) has no
DCM mode to transition into without a separate device model this round
did not build; stated here as a substitution, not hidden.

**Correction, 2026-10-02 (card #550).** The two continuation lines of the `Iload_step` PWL in this fixture had no leading `+`
(a tab-less indent is not a continuation in ngspice 45.2, the pre-master-48 base or this fork), so ngspice dropped them and the source was an
unterminated `PWL(0 0 1.5m 0` that never stepped: the reference `.tran` has no load step either, and every measurement in the bullets below
(the straddling jump, "120 refusals", the 0.6-2.9 % error "near the event") described ordinary compensator ringing and landing kicks, not the
step. The old test passed on the development workstation only by a gate firing 94 us before a non-existent event, and not at all on a second workstation (different last bits of
exp/log, a different ringing). Fixed (`+` added); with the 4 A step real, both workstations refuse exactly at it (refusals at 1498-1510 us on
the development workstation, 1500-1510 us on the second workstation: `mid_check_event` x4 and `growth_tol` x4 within [-20, +30] us of the event; no jump straddles it in either run) and
the worst cycle-average error near the event is 5.17 % (development workstation) / 5.12 % (second workstation), bound 20 %. The test now generates its reference .tran itself.
The qualitative statements below about skipped spans remain true of the method, but the numbers are not measurements of the step.

- The pre-jump growth-tol gate and the last-two-check-periods slope-tol
  gate, working alone, were NOT sufficient: a first pass found one jump
  (M=10) landing so that the 4A step fired and largely re-settled entirely
  INSIDE that jump's own 4-period check window -- `slope_err=0.0008`,
  comfortably under `slope_tol`, because the converter's fast switching
  dynamics re-settle to a new ripple pattern within 1-2 periods even
  though the compensator's slow pole means the AVERAGE keeps moving for
  much longer. **A red-proof event can be invisible to an end-to-end
  before/after slope comparison if it starts and mostly finishes inside
  the very check window meant to validate the jump.**
- Fixed by adding a scan across EVERY consecutive pair of real periods
  actually simulated inside each check window (the same history-relative
  `is_smooth()` test the pre-jump gate already uses), not just its
  endpoints. After this fix: of 120 total refusals in a 0-1.6ms run
  spanning the event, 114 are `reason=mid_check_event` -- the added scan
  is doing real, load-bearing work, not a no-op.
- **Even so, one jump (M=10, landing at t=1.506ms, i.e. 6 periods AFTER
  the event at t=1.500ms) still straddled the event UNDETECTED** --
  because the event fell entirely inside the SKIPPED span itself (the 10
  periods between t=1.486ms and t=1.506ms that are never simulated at
  all when a jump succeeds), not inside a check window this tool ever
  looks at. **This is an inherent limitation of envelope-following as a
  class of method, not specific to this implementation**: by
  construction, nothing that happens strictly inside a skipped span can
  be observed without simulating it, which is the entire point of
  skipping it. Cadence's own SpectreRF user guide states the same
  trade-off plainly: "if the analysis can not find cycles to skip, it
  reduces to transient analysis" -- the method's soundness rests on the
  ASSUMPTION that skipped spans are smooth, not a guarantee of catching
  every possible violation. Measured impact here, not just described:
  the missed straddle produced a **bounded ~0.6-2.9% cycle-average error**
  in `v(vout)` over the next several periods (`t=1.480ms..1.512ms`,
  against the deck's own full reference run) -- a real, out-of-tolerance
  error, but a bounded, graceful one, not a blow-up, because M is capped
  at 10 periods (20us) by this tool's own default rather than the
  hundreds the card's own framing discusses. Keeping M modest is this
  round's own practical mitigation for a limitation that cannot be fully
  engineered away from outside the skipped span.

## Tests

`tests/test_envelope.py` (new, uncommitted, per this card's own
instructions): `ckpt_io` round-trip byte-identity; a loadstate-of-a-patch
sanity check; the M=1 checkpoint-chaining-only accuracy control; the
full-deck accuracy/speed measurement (reported, not gated pass/fail, given
the characterised shortfall above -- gating it "pass" would misrepresent
a real, measured result); and the red-proof refusal check (asserts
`mid_check_event`/`growth_tol`/`slope_tol` refusals occur near the event,
and reports -- does not silently hide -- the one straddling jump found).
Run with `NGSPICE=<a build with patches 0032-0033, e.g. s55> python3
tests/test_envelope.py`.

## What a future round should do differently (ranked)

1. **An implicit, Newton-corrected envelope step** -- literally Kundert's
   own algorithm (`docs`/this round's own literature read of Cadence's "A
   User's Guide to Envelope Following Analysis": a 3-point quadratic
   interpolation through `x(t_(n-1))`, `x(t_n)`, and the unknown
   `x(t_(n+1))`, solved by shooting-Newton using `∂φ/∂x(t_(n+1))` -- the
   ONE-PERIOD MAP'S OWN JACOBIAN, i.e. this project's own monodromy
   machinery from `pssshoot.c`/`docs/pss-pac.md`, already built and
   validated for the PSS case) instead of this round's explicit
   forward-Euler-in-cycle-index step. This directly targets the landing
   kick: an implicit step is SOLVED to be self-consistent with the true
   one-period dynamics (using the monodromy Jacobian's own `(Phi-I)`
   structure, the same singular-near-a-slow-pole matrix `docs/
   pss-pac.md`'s own "R4 D~=0.45" case study already characterises), not
   guessed and then algebraically reconciled after the fact. This is a
   substantially bigger build (reusing, not rewriting, the existing FD
   monodromy code, but needing new glue to drive it from a NON-periodic,
   time-varying-source starting point) -- ranked first because it
   addresses the actual measured root cause, not a symptom.
2. **Chain more real work per ngspice invocation** to close the
   speed gap between "70% fewer periods simulated" and "2.6x wall
   speedup" -- run all `n_check` check periods (and/or several
   passthrough periods) in ONE `.tran` call instead of one subprocess per
   period, amortising process-launch/parse overhead the way a real
   production tool would.
3. **A physically-motivated partial `CKTstates[]` extrapolation** -- even
   without the full implicit step above, extrapolating just the
   compensator's own dominant slow-pole state (identifiable from the
   netlist: the largest R*C time constant among the feedback network's
   own capacitors) alongside `CKTrhsOld` would likely shrink the landing
   kick substantially, as a cheaper intermediate step before (1).
4. **PCM buck (r4d) and line-ramp variants** -- this round only built the
   VM buck (r3d) soft-start case; the card also named line ramps and the
   peak-current-mode deadtime deck.

---

# Implicit envelope following inside the fork (card #513, round 11, 2026-09-29)

> Status 2026-09-29: parked in `patches-wip/env/` (not in the series); see its README.

Round 10's explicit projective jump (above) failed on accuracy (stale
`CKTstates`) and speed (1.2x, process overhead). Round 11 builds
Kundert/White/Sangiovanni-Vincentelli's implicit envelope step natively:
`src/spicelib/analysis/pssenv.c` on branch `aug-env2` (not in the `patches/`
series yet). **Result in one paragraph:** on a deck whose one-period map is
smooth (`tests/fixtures/env_lin_chopper.cir`, linear chopper + LC) it meets
the card's bar -- every simulated period's cycle-average vout / iL within
0.5 % and ripple within 5 % of full `.tran 2n`, wall 8.3x faster than that
full `.tran` (1.8x faster than an equal-resolution period-by-period chain);
a declared load step in the middle of the run is handled (jumps never land
within 2 periods of a PWL corner) and a hidden one is caught by the error
estimate. On the card's own target deck, the closed-loop soft-start buck
`env_r3d_softstart.cir`, **the bar is not met and cannot be met by this
method as it stands**: the stroboscopic one-period map of that deck is not
smooth at the level Newton needs (evidence below), so the engine refuses
every jump, backs off, and degrades to a plain period-by-period chain
(no speed-up, accurate only as far as that chain and the deck's own
discretisation sensitivity allow).

## Syntax

    .tranenv Fguess Points SC_iter Steady_coeff Tstop Mmax [envtol=<tol>] [phase=<0..1>] [uic]

A new dot-card (own parser function, own analysis path) -- not
`.tran ... envelope=`: `.tran`'s grammar and engine are shared by every
deck, a new directive provably cannot change any of them (and
`tranenv_plain_tran_unchanged` checks the plain `.tran` output is
byte-identical to the s61 baseline). `Fguess` clock frequency, `Points`
integration points per period (max step T/Points), `SC_iter` Newton budget
(capped at 8), `Steady_coeff` Newton tolerance (scaled step or residual,
default 1e-4), `Tstop`, `Mmax` largest jump in periods, `envtol` local
error tolerance relative to the state scale (default 5e-4), `phase`
sampling phase as a fraction of a period (period maps run
[kT+phase*T, (k+1)T+phase*T]; default 0; 0.5 avoids restarting in the
middle of a hard clock edge). Needs a driven (clocked) circuit -- else the
error "`.tranenv needs a driven (clocked) circuit`" -- and >= 1 native C or
L. Output: **one** ordinary time-domain plot holding every real simulated
period (gaps in time across the skipped spans; `wrdata`/`plot`/`meas`
work), plus a small "Envelope Samples" plot with the states at each real
period boundary. `ENV_DEBUG=1` prints every decision on stderr.

## Design

* **State / period map.** The reduced state of `.pss method=shoot` (all
  native C voltages and L currents, n values). Every period map
  F(x; t0) is a fresh `DCtran()` restart (uic) from those values at an
  ABSOLUTE start time t0. Nothing from a previous map survives (no
  `CKTstates` history), so round 10's "stale state after the jump" cannot
  occur. Needed changes to run at t0 != 0: `dctran.c` (restart branch,
  `envT0`, incl. the uic first-step breakpoint `CKTsetBreak(ckt, CKTstep)`
  that was absolute -- the cause of a "breakpoint in the past" panic found
  and fixed this round), and the per-instance `VSRCbreak_time`/
  `ISRCbreak_time` bookkeeping is reset before every map (a stale value
  suppresses PULSE/PWL edges or panics).
* **Implicit step.** Solve y - x_n = M (F(y; t_n+M T) - y) by Newton,
  J = (1+M) I - M Phi. Phi by forward FD of the free map (step 2e-3 of the
  state scale), cached across steps and across M: chord Newton, one
  period run per iteration, refreshed when 12 solves old, when an
  iteration fails to halve the residual, or when a step is rejected.
* **Error control on the slow states only** (Phi_kk >= 0.9; the fast
  re-slaved states are handled exactly by the implicit step and their
  jitter would otherwise pin M at 2-3): M <= sqrt(0.6 envtol scale/|c|)
  (c = last change of the per-period slope), post-solve check
  |y - explicit projection| <= envtol scale (else M shrinks by
  sqrt(envtol/err) and the step is refused), doubling ladder to Mmax,
  history-relative smoothness gate (delta vs median of last 8 deltas).
* **Events.** A jump never lands within 2 periods of a corner of any
  native PWL V/I source (soft-start knees, load/line steps); real periods
  cross the corner and the post-corner transient keeps M small through the
  curvature predictor. A hidden persistent event inside a jump shows up in
  F(y) at the landing and is refused by the error estimate. A transient
  wholly inside a skipped span is invisible to any envelope method (same
  limitation as round 10).
* **Fall back to M=1.** Newton failure, singular Jacobian, period-map
  error, repeated error refusals: plain real period, then exponential
  back-off of 4, 8 ... 128 real periods before the next attempt, so a
  deck Newton cannot handle costs ~1.3x a plain chain instead of paying a
  failed solve per period.
* **Memory.** Scratch runs (Newton/FD) produce no plot at all
  (`PSSenvOutMode` 1). The first version made one ngspice plot per scratch
  run and reached 19.8 GB RSS in 86 s (thousands of plots); now 30 MB
  (linear deck) and 235 MB (closed-loop deck, 2000 periods x 200 points
  stored in the one kept plot). `test_envelope.py` checks max RSS < 1 GB.

## Measured (interleaved A/B, load average recorded; OMP_NUM_THREADS=2)

Reference for accuracy = full `.tran 2n` (the fixture's own tstep); the
`.tranenv` runs use 200 points/period (10 ns max step), `envtol=5e-4`,
`Steady_coeff=1e-4`, `phase=0.5`, `Mmax=100`. "Simulated periods" = periods
actually integrated and emitted (the rest are skipped, not interpolated).
Errors are worst cases over all simulated periods, relative with a floor of
2 % of the run maximum. Run of `tests/test_envelope.py` at load1 ~6.

| deck | simulated / 2000 | jumps | vout avg | iL avg | ripple | verdict |
|---|---|---|---|---|---|---|
| env_lin_chopper (smooth) | 810 | 188 | 0.08 % | 0.27 % | 2.3 % | meets 0.5 % / 5 % |
| + PWL load step at 1.451 ms | 951 | 199 | 0.12 % | 0.28 % | 3.2 % | meets; +-100 periods around the event: 0.07 / 0.13 / 2.0 % |
| + same step as time expression (hidden) | 951 | 199 | 0.12 % | 0.28 % | 3.2 % | meets (error estimate refuses the straddling jump) |
| env_r3d_softstart (closed loop) | 1999 | 0 | 0.38 % (periods 300+) | 9.3 % | 29 % | NOT met; steady periods 1600+: 0.03 % / 0.98 % |

| deck | full `.tran 2n` | `.tranenv` | speed-up | vs equal-resolution chain (Mmax=1) |
|---|---|---|---|---|
| env_lin_chopper | 3.61 s (median of 3) | 0.43 s (median of 3) | **8.3x** (>= 5x) | 1.8x (0.80 s) |
| env_r3d_softstart | 8.6 s (single run, load1 3.9) | 3.0 s (period chain, no jump accepted) | 2.8x, but ALL of it is the coarser 10 ns step, none is envelope following | ~1.0x (about 1.1-1.3x overhead from the refused Newton attempts) |

(The 8.3x is mostly two things: fewer periods simulated -- 810 of 2000 --
and the coarser 10 ns step cap versus 2 ns; the linear deck's answer is
insensitive to that step, the 1.8x against the equal-resolution chain is
the part due to envelope jumps alone. State which one you mean when
quoting the number.)

## Why the closed-loop soft-start deck does not work (measured)

1. **The deck is discretisation-sensitive at the level of the bar.** Plain
   full `.tran` with tstep 5 ns vs the 2 ns reference: iL cycle-average
   -6 % (period 150-200), 33 ns: -27 %..+5 %; on the open-loop comparator
   deck (`env_openloop_buck.cir`, kept as a fixture) 2 ns vs 5 ns already
   differ by 8 % in iL and 20 % in ripple. Any restarted or coarser
   simulation, envelope or not, misses the 0.5 % / 5 % bound there.
2. **The stroboscopic map is not smooth at Newton's scale.** FD columns of
   Phi at relative steps 1e-5 and 1e-3 differ by 10-100x for the
   compensator states (e.g. dF_iL/d c3in = 458 vs 12 A/V): a 1e-5 V change
   of the compensator's c3in state changes the next-period iL by ~5 mA
   (0.4 %). The sampled c3in/coa/c1f sequences show sustained
   period-to-period jitter of that size (median |second difference| of
   coa ~1e-3 V at steady state), matching the +-0.3 % iL jitter visible in
   the full `.tran` reference itself. Newton on that map fails (20 of ~30
   attempts refused in the test run) and the error estimate refuses the
   rest.
3. Two Jacobian routes were tried and dropped: FD on a forced time grid
   (tiny steps: the transient solver's convergence-tolerance noise in F,
   divided by 1e-9..1e-7 V, produced entries of 1e2..1e8 -- the forced and
   free runs also differ by ~1e-3 because a forced run restarts at order 1
   at every recorded time) and the analytic (discretise-then-linearise)
   monodromy of `pssshoot.c` (documented inaccurate near high-gain
   comparators in docs/pss-pac.md; it also aborted in `spSolve` on the
   first sensitivity run at t0 > 0 -- not chased).
4. Next steps if this deck matters: a slow-manifold reduction (Newton only
   on the integrator states with the fast states slaved by a short
   settling run, `docs/pss-pac.md`'s fast-state reduction), a smoother
   comparator (kc, kswd) model, or accepting the chain result there.

## Tests / how to run

`NGSPICE=<aug-env2 build> python3 tests/test_envelope.py` (about 1 min):
the round-10 checks plus `tranenv_linear_accuracy_speed`,
`tranenv_red_proof_pwl_step`, `tranenv_hidden_event_bounded`,
`tranenv_closed_loop_degrades`, `tranenv_rejects_undriven`,
`tranenv_plain_tran_unchanged` (10/10 PASS + 1 INFO on the 2026-09-29
run). `tools/envelope_cmp.py` computes per-period window statistics from
`wrdata` files and compares runs.

## What was not done (round 11)

Broyden/analytic Jacobian, slow-manifold reduction of the Newton system,
non-PWL slow-source corner detection (a PULSE/SIN/expression event is only
caught by the error estimate), a non-uic start, autonomous circuits,
`patches/` series entry (the work is commits on `aug-env2`), the r4d (PCM)
and line-ramp variants, and any run of the closed-loop deck with a tighter
transient tolerance.

---

# Closed-loop soft-start with `.tranenv` (card #719, 2026-10-23)

**Result.** The closed-loop soft-start decks now run faster than plain
`.tran` at a stated, measured accuracy. On `env_r3d_softstart.cir`
(voltage-mode buck, dead-time network, 2000 periods) `.tranenv` simulates
171 of 2000 periods (60 jumps) in 1.09 s; every simulated period's cycle
average of vout / iL and its vout ripple is within 0.014 % / 0.16 % / 0.14 %
(periods 100+) of a plain `.tran` at reltol 1e-7. The deck's own plain
`.tran 2n` takes 6.81 s and is 0.57 % / 7.6 % / 13.7 % off (6.2x slower and
10-100x less accurate); the cheapest plain `.tran` at least as accurate as
`envtol=1e-4` (10 ns, reltol 1e-5) takes 6.06 s against 1.60 s (3.8x). On the
second closed-loop deck, the peak-current-mode buck
`tests/fixtures/envelope719/env_r4_pcm_softstart.cir`, `.tranenv` only
draws level: about the same wall time as a plain `.tran` of the same accuracy
(1.0-1.05x), or 5-7x more accurate at the same wall time as the deck's plain
`.tran 10n`. What made the difference was not a new map or Jacobian. Every
period map now runs at a tight transient tolerance (`maptol`, default reltol
1e-5). Card #513's engine (`patches-wip/env`) is ported to the s150 code with
that change and three smaller ones (below).

## Facts on s150 (before)

* Round 10's Python driver (`tools/envelope_tran.py`, the only envelope
  option in s150): 806/2000 periods, 185 jumps, worst vout cycle-average
  error 5.2 % against the deck's own `.tran 2n` (port-s150 log, rerun here
  with the same numbers).
* Round 11's `.tranenv` ported unchanged onto s150 (`ENV_MAPTOL_OFF=1` on
  the new build gives exactly that engine): 1999/2000 periods simulated,
  **0 jumps**, 8 Newton refusals, 2 error refusals, as card #513 found.
* The reference. Plain `.tran` 1 ns at reltol 1e-7 (abstol 1e-13,
  vntol 1e-10) against 2 ns at reltol 1e-6 (the halved-step and tighter-tol
  check). The worst difference over all 2000 periods is 0.005 % vout, 0.05 %
  iL and 0.14 % ripple, in periods < 30; from period 100 on it is 0.001 % /
  0.011 % / 0.11 %. r4: 0.0015 % / 0.015 % / 0.13 %. All errors below are
  against the 1 ns / 1e-7 run. Errors are worst cases over every simulated
  period, relative to max(|ref|, 2 % of the run maximum), on windows
  [(k+0.5)T, (k+1.5)T] (the `.tranenv` period maps at phase=0.5).
* **The deck is tolerance-sensitive, not step-sensitive.** The plain
  `.tran 2n` reference that rounds 10 and 11 measured against is itself
  0.57 % / 7.6 % / 13.7 % off (periods 100+; 15 % iL and 24 % ripple in the
  first periods). At reltol 1e-5, 10 ns gives 0.005 % / 0.05 % / 0.7 %.
  Round 11's "5 ns vs 2 ns differ by 6 % in iL" was the reltol 1e-3
  noise of both runs.

## Why no jump was accepted (measured)

The cycle map is rough at the default reltol 1e-3 because of the step
control: the switch edges are 1-2 ns wide (tanh, kswd=150) and the LTE
control places them differently from one period to the next. Event
location is not involved, since these decks have B-source switches and no
`sw` instances (PSSeventLocate handles native `sw` only).

Period-to-period jitter of the stroboscopic samples (phase 0.5), median
|second difference| of iL over periods 700-900. The true sequence is smooth:
the ramp is linear.

| run | iL median abs d2 | vout median abs d2 |
|---|---|---|
| plain .tran 1n reltol 1e-7 (reference) | 9.2e-7 A | 4.3e-9 V |
| plain .tran 2n reltol 1e-3 (deck default) | 5.1e-4 A | 1.8e-5 V |
| period chain (.tranenv Mmax=1) reltol 1e-3 | 4.7e-4 A | 5.1e-5 V |
| period chain, fixed 10 ns grid (trial, dropped) | 4.7e-4 A (p90 1.2e-2: limit cycle) | 6.1e-5 V |
| period chain reltol 1e-5 | 1.1e-6 A | 1e-8 V |
| period chain reltol 1e-6 | 1.4e-7 A | 1e-8 V |

`ENV_SCAN=1.5e-3` gives the finite-difference column dF_iL/dv(Cout) of the
one-period map for relative steps 1e-8 ... 1e-2 of the state:

| map | 1e-8 | 1e-7 | 1e-6 | 1e-5 | 1e-4 | 1e-3 | 1e-2 |
|---|---|---|---|---|---|---|---|
| reltol 1e-3 (round 11) | -0.19 | -0.19 | -0.19 | -0.19 | -2.28 | -7.43 | -7.47 |
| fixed 10 ns grid, reltol 1e-3 | 0.055 | 0.055 | 0.049 | 0.0002 | 0.70 | -14.5 | -8.1 |
| reltol 1e-5 (new default) | 6.9 | -5.8 | -7.72 | -7.80 | -7.85 | -7.83 | -7.68 |
| reltol 1e-6 | -19.4 | -9.0 | -7.80 | -7.85 | -7.85 | -7.84 | -7.68 |

At reltol 1e-3 a small perturbation never moves the step grid, so the
finite difference sees the map on a frozen grid. That is the duty-cycle
sensitivity without the edge-timing part, 40x too small, so Newton gets a
wrong Jacobian. A large perturbation crosses grid changes. A line scan of
v(Cout) over +-1e-4 of its value has a residual of 5.8e-5 A around a straight
line at reltol 1e-3 (3e-3 A with the fixed grid, and a jump), against
6.9e-7 A at 1e-5 and 1.3e-7 A at 1e-6. The curvature predictor
(M <= sqrt(0.6 envtol scale / |second difference|)) sees the 5e-4 A jitter
and holds M below 1 everywhere. ENV_DEBUG=4 on the fixed-grid chain shows
Mpred 0.02-0.7 over 0-2 ms. The reltol 1e-3 adaptive chain has the same
jitter and accepted no jump either.

## Levers tried, ranked by measured smoothness

1. **Tight period-map tolerance (kept).** Every period map (bootstrap,
   Newton, finite differences, landing, plain periods) runs at reltol =
   min(deck reltol, maptol), with maptol 1e-5 by default. `maptol=` on the
   card overrides it, and `maptol=0` keeps the deck's reltol. r3d,
   envtol 5e-4, maptol swept (an intermediate build with n_check=2 and plain
   doubling): reltol 1e-3 0 jumps; 1e-4 828 periods / 90
   jumps / 4.1 s; 1e-5 237 / 63 / 1.27 s; 1e-6 238 / 64 / 1.88 s;
   1e-7 241 / 66 / 3.27 s, all with the same accuracy from 1e-5 down.
   abstol and vntol do not matter (reltol alone gives the same accuracy
   as reltol + abstol 1e-12 + vntol 1e-9).
2. **Fixed step grid (tried, dropped).** No LTE control; every map takes
   the same steps (the maximum step, clipped by source breakpoints). The
   map no longer depends on step decisions, but a 10 ns grid on a 1 ns
   edge is a steep staircase: the chain falls into a quantisation limit
   cycle (iL sawtooth of +-7 mA with a period of about 7 periods at steady
   state), the columns are worse than the adaptive ones, and no jump is
   accepted.
3. **Grid replay, as in the PSS F evaluations (card #505), or a smoothed or
   least-squares Jacobian: not needed.** With lever 1 the free adaptive map
   gives finite-difference columns that agree from 1e-6 to 1e-3, and Newton
   converges with the cached chord Jacobian (0 Newton refusals on both
   decks).
4. **PSS-based quasi-steady orbit tracking the reference: not built.** A
   frozen-reference PSS has vout = vref (1 + Ra/Rb) exactly. It misses the
   loop's ramp-tracking lag and the C dv/dt current, which the implicit
   envelope step does capture.

Smaller changes made while measuring (all in pssenv.c):

* After a landing the curvature is taken between the jump's secant slope
  (y - x_n)/M and the landing period's slope, so the next jump needs no
  extra real period. The confirmation period becomes optional:
  `ENV_NCHECK`, default now 1, was 2. r3d at envtol 1e-4: 2.87 s ->
  1.55 s, interleaved (555 -> 279 periods).
* The next M is M min(2, 0.9 sqrt(envtol / err_est)) instead of 2M. Plain
  doubling after a step that used 60 % of envtol alternated accept/refuse
  on the PCM deck (85 error refusals, now 1).
* The Jacobian seed run (`ENV Phi`) was written to the output plot, which
  put one duplicated period (time going backwards) into every `.tranenv`
  output. It is now a scratch run.
* `ENV_SCAN=<t>` prints the column and line-scan diagnostic above.
  `ENV_DEBUG=4` names the state that limits M. The DONE line reports
  `rejected_steps` and `map_reltol`. The `.tranenv` parser no longer
  leaks its tokens.

## Accuracy (measured)

Plain `.tran` at 1 ns, reltol 1e-7 is the reference. `.tranenv 500k 200 8
1e-4 4m 100 envtol=<e> phase=0.5 uic` uses the default maptol 1e-5. The
columns give the worst error over every simulated period with k >= 100,
then over all k, then at the end of the soft-start (periods 1490-1510).

r3d (`env_r3d_softstart.cir`):

| run | periods | vout k>=100 | iL k>=100 | ripple k>=100 | vout / iL / ripple all k | end of ramp vout / iL |
|---|---|---|---|---|---|---|
| plain .tran 2n default (deck) | 2000 | 0.57 % | 7.6 % | 13.7 % | 0.65 / 15.4 / 23.9 % | 0.008 / 0.20 % |
| plain .tran 10n default | 2000 | 1.7 % | 18 % | 40 % | 5.9 / 38 / 48 % | 0.028 / 1.1 % |
| plain .tran 10n reltol 1e-4 | 2000 | 0.098 % | 0.86 % | 1.96 % | 0.11 / 1.4 / 2.1 % | 0.0006 / 0.039 % |
| plain .tran 10n reltol 3e-5 | 2000 | 0.013 % | 0.047 % | 0.70 % | 0.072 / 0.70 / 1.6 % | 0.0002 / 0.041 % |
| plain .tran 10n reltol 1e-5 | 2000 | 0.0049 % | 0.047 % | 0.68 % | 0.034 / 0.31 / 0.89 % | 0.0002 / 0.043 % |
| .tranenv envtol 2e-3 | 127 | 0.056 % | 0.55 % | 0.49 % | 0.056 / 0.55 / 0.49 % | 0.0002 / 0.007 % |
| **.tranenv envtol 5e-4 (default)** | 171 | **0.014 %** | **0.16 %** | **0.14 %** | 0.035 / 0.32 / 0.37 % | 0.0002 / 0.006 % |
| .tranenv envtol 1e-4 | 283 | 0.0046 % | 0.035 % | 0.054 % | 0.035 / 0.32 / 0.37 % | 0.0001 / 0.005 % |
| .tranenv, ENV_MAPTOL_OFF=1 (red proof) | 1999 (0 jumps) | 1.97 % | 20.6 % | 46 % | 4.8 / 37 / 46 % | 0.033 / 1.5 % |

r4 PCM (`envelope719/env_r4_pcm_softstart.cir`, Vin 11 V, smooth SR latch,
type-II loop):

| run | periods | vout k>=100 | iL k>=100 | ripple k>=100 | end of ramp vout / iL |
|---|---|---|---|---|---|
| plain .tran 10n default (deck) | 2000 | 0.19 % | 1.36 % | 2.74 % | 0.0025 / 0.17 % |
| plain .tran 2n default | 2000 | 0.0014 % | 0.015 % | 0.13 % | 0.0000 / 0.008 % |
| plain .tran 10n reltol 1e-4 | 2000 | 0.018 % | 0.089 % | 0.78 % | 0.0003 / 0.048 % |
| plain .tran 10n reltol 1e-5 | 2000 | 0.0048 % | 0.052 % | 0.74 % | 0.0002 / 0.033 % |
| .tranenv envtol 2e-3 | 190 | 0.069 % | 0.67 % | 0.87 % | 0.035 / 0.32 % |
| **.tranenv envtol 5e-4 (default)** | 369 | **0.028 %** | **0.31 %** | **0.41 %** | 0.020 / 0.18 % |
| .tranenv envtol 1e-4 | 744 | 0.0070 % | 0.081 % | 0.11 % | 0.0015 / 0.006 % |
| .tranenv, ENV_MAPTOL_OFF=1 | 1610 (99 jumps) | 0.19 % | 1.37 % | 2.46 % | 0.0024 / 0.14 % |

**Error bound and why.** envtol bounds the local error of each jump
relative to each slow state's scale. Measured over the whole soft-start,
the global error on simulated periods 100+ is at most 0.03 % vout,
0.35 % iL and 0.45 % ripple at envtol 5e-4, and 0.007 % / 0.08 % / 0.11 %
at envtol 1e-4, on both decks. That is roughly envtol x 0.6 in vout and
envtol x 7 in iL. iL is the small difference of two large slow states
(vout/R + C dv/dt). The worst iL error of every run lies in periods
1512-1539, right after the ramp knee, where the knee transient moves the
accumulated envelope error. The worst vout error on r3d lies in periods
100-160 (the first jumps). The first ~20 periods (startup from rest) are
real periods in every run, at the map tolerance. Their error is that of a
reltol=1e-5 transient (0.32 % iL, the same as a plain reltol 1e-5 run).
The reference contributes <= 0.011 % iL, 0.11 % ripple. Skipped periods
are not output and not interpolated; the bound covers simulated periods
only. A transient that starts and ends inside a skipped span is invisible
(card #513; PWL corners are never jumped over).

## Speed (interleaved A/B)

Three alternations in the order shown, one ngspice at a time,
OMP_NUM_THREADS=2, all on the same w719 build (plain `.tran` output is
byte-identical to s150). Load average at the start 1.73-1.93 (one-minute
value; another user was simulating).

| deck / run | wall (3 runs) | median | accuracy (k>=100 vout / iL / ripple) |
|---|---|---|---|
| r3d plain .tran 2n default (deck) | 6.88 6.81 6.79 | 6.81 s | 0.57 / 7.6 / 13.7 % |
| r3d plain .tran 10n reltol 1e-5 | 6.03 6.10 6.06 | 6.06 s | 0.0049 / 0.047 / 0.68 % |
| r3d plain .tran 10n reltol 1e-4 | 3.54 3.54 3.56 | 3.54 s | 0.098 / 0.86 / 1.96 % |
| r3d .tranenv envtol 5e-4 | 1.13 1.09 1.07 | **1.09 s** | 0.014 / 0.16 / 0.14 % |
| r3d .tranenv envtol 1e-4 | 1.60 1.59 1.63 | **1.60 s** | 0.0046 / 0.035 / 0.054 % |
| r4 plain .tran 10n default (deck) | 1.71 1.69 1.69 | 1.69 s | 0.19 / 1.36 / 2.74 % |
| r4 plain .tran 10n reltol 1e-4 | 2.26 2.27 2.28 | 2.27 s | 0.018 / 0.089 / 0.78 % |
| r4 .tranenv envtol 5e-4 | 1.66 1.66 1.64 | 1.66 s | 0.028 / 0.31 / 0.41 % |
| r4 .tranenv envtol 1e-4 | 2.16 2.15 2.21 | 2.16 s | 0.0070 / 0.081 / 0.11 % |

r3d: 6.2x faster than the deck's own plain `.tran` and more accurate on
every measure. At matched accuracy (envtol 1e-4 against plain 10 ns
reltol 1e-5, at least as accurate in vout, iL and ripple) it is 3.8x; at
envtol 5e-4 against plain reltol 1e-4 (less accurate on all three) it is
3.2x. r4: 1.05x at matched accuracy (envtol 1e-4 against plain 10 ns
reltol 1e-4, which is less accurate). It is 1.0x against the deck's plain
`.tran 10n` while 5-7x more accurate.

The PCM deck gains less for two reasons. Its loop is still settling
during the first millisecond (c2f, the 330 us type-II integrator, limits M
to 1-3 for 500 periods), and a 10 ns default-tolerance plain run of it is
cheap. On the linear chopper (card #513's smooth deck) the port gives
card #513's numbers with n_check=2 (810/2000 periods, 0.08 / 0.23 / 2.3 %).
With the new n_check=1 it simulates 644 periods and stays within
0.10 / 0.30 / 1.9 %.

## Red proof

`ENV_MAPTOL_OFF=1` restores the card #513 behaviour: period maps at the
deck's reltol. On r3d that gives 0 jumps in 2000 periods and the errors of a
default-tolerance chain (20.6 % iL). On r4 it gives 99 jumps but 1610
periods simulated, 3.7 s, and default-tolerance errors (1.37 % iL).
`tests/test_envelope.py` checks the r3d case
(`tranenv719_red_proof_maptol_off`).

## Tests

`tests/test_envelope.py` adds `tranenv719_closed_loop_softstart` (r3d
accuracy bounds, >= 20 jumps, <= 25 % of the periods simulated, >= 2x
faster than the deck's plain `.tran 2n` in an interleaved pair),
`tranenv719_red_proof_maptol_off`, `tranenv719_pcm_softstart` (r4
accuracy bounds, speed reported only), `tranenv719_linear_chopper`,
`tranenv719_rejects_undriven`, and `tranenv719_plain_tran_unchanged`, which
runs only with NGSPICE_BASE and compares byte-identical plain `.tran`
output. The references (2 ns, reltol 1e-6) are cached under
`$ENVELOPE_OUT/env719`, a new override; the default output folder is
unchanged. The default binary is s151, with NGSPICE overriding it; a build
without `.tranenv` gives COULD-NOT-RUN. On the w719 build: 11 checks,
10 PASS, 1 INFO (the round 10 measurement), 1 min 43 s including the
references. On s150: 4 PASS, 1 INFO, 5 COULD-NOT-RUN, exit 2. Copies of
test_pss / test_pac / test_tran_from_pss with every data path redirected
match the port-s150 logs (44/48 + 4 info, 35/35 + 2 info, 6/6), with
identical result lines apart from wall times.

## What was not done (card #719)

* `.tranenv` is still a separate directive, not `.tran ... envelope=`.
  Plain `.tran` runs unchanged code: the new branches are inactive unless a
  `.tranenv` job is running, and byte-identical output is checked on three
  decks.
* No `patches/` series entry or s151 build was made: the work is on branch
  w719. A third closed-loop deck (r3, no dead time) was not run.
* The slow-state error scale is still |x| floored at 1 % of the largest
  state. A per-state range would let the PCM deck's small c2f state take
  larger M; this was not tried.
* No line-ramp or load-ramp variant was tried.
* The Envelope Samples plot is not created when the deck's `save` list
  leaves none of its vectors. The test decks have such a list, so the
  samples were not checked.
* No ASan run of the new code was made: there is no ASan twin of w719.
* Memory was not gated in the new tests. ngspice memory was not measured;
  all runs stayed under the 3 GB scope cap.
