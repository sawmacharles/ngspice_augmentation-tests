# checkpoint/restart -- transient state save/resume (card #472, plan T10.4)

Adopts and fixes upstream ngspice pre-master-48's `savestate`/`loadstate`
commands (Enhancement-131, `src/frontend/com_checkpoint.c`): dump a running
`.tran`'s complete state to a file, and resume it later -- in a fresh
process, possibly on a different invocation entirely -- from exactly where
it left off. Two bugs in the upstream mechanism are fixed here
(first written as `patches-wip/ckpt/0200`, `0201`; promoted into the series
as `patches/0032`, `0033` on 2026-09-29 and present in every build from s33
on -- see the patch commit messages for the full characterisation).

## What is saved, and what is NOT (characterisation)

Read directly from `src/frontend/com_checkpoint.c` and
`src/spicelib/analysis/dctran.c` (pre-master-48 pin `aaba50ba`, this
project's `aug` branch = pin + `patches/` 0001-0016). `savestate <file>`
writes, and `loadstate <file>` restores:

- the solution vector (`CKTrhsOld`, `CKTirhsOld` if AC-analysis-shaped
  imaginary state is present -- not normally, for a plain `.tran`)
- the integration history: all (up to eight) `CKTstates[]` device-state
  vectors -- the past-step charge/flux/current history every integrator
  (trapezoidal, Gear 2-6) needs
- the scalar integration state: `CKTorder`, `CKTmode`, `CKTtime`,
  `CKTdelta`, `CKTsaveDelta`, `CKTdelmin`, `CKTminBreak`, `CKTdeltaOld[7]`
- the pending breakpoint table, `CKTbreaks[]`

**Not saved, and confirmed not needed** (characterised by reading
`src/maths/ni/nicomcof.c`, not assumed): `CKTag[]`/`CKTagp[]` (the
Gear/trapezoidal coefficient matrix and its predictor) are a *pure
function* of `CKTdelta`, `CKTdeltaOld[]`, `CKTorder` and the integration
method/`CKTxmu` (all either restored or fixed by the build/options) --
`NIcomCof()` recomputes them from scratch before every step, called
unconditionally at the top of the stepping loop (`dctran.c:755`) regardless
of whether the run started fresh, resumed in-memory, or resumed from a
checkpoint. There is nothing to save here.

**Not saved, and a real limitation** (confirmed, not just inferred):
XSPICE event-driven state (`g_mif_info`, the event queue) and digital/hybrid
code-model internal state are outside `CKTstates[]`/`CKTrhsOld` entirely
and are not touched by `savestate`/`loadstate` at all. A checkpoint taken
mid-run on a deck using XSPICE digital/event code models (e.g. refs/
r6_vm_buck_ideal.cir's `d_pwm`) would restore the analog solution correctly
but silently lose the event-driven side's state. **Not tested here** (out
of this card's four target decks) -- flagged as a real gap for any future
work that wants to checkpoint an XSPICE-event-heavy deck, not verified to
be safe or unsafe.

**Solver restriction** (in the code, both directions): `savestate` and
`loadstate` both refuse to run under `.option klu` ("checkpoint/restart is
only supported with the default Sparse solver, not KLU"). None of this
card's four target decks use `.option klu`, so this was not exercised
further.

**File format**: a small binary header (an 8-byte magic + a topology
signature: equation count, matrix size, state count, max order, breakpoint
count) followed by the raw arrays above, written with `fwrite`/native type
sizes -- not portable across machines/builds with different type sizes, and
rejected via the signature if the circuit it is loaded into does not match
structurally (same node/device count) the one that wrote it. Confirmed:
loading a checkpoint that does not structurally match the live circuit is
rejected with a clear error, not silently misapplied.

## The two bugs found and fixed (patches-wip/ckpt/0200, 0201)

Found by characterising restart *exactness* against an uninterrupted run on
`refs/r2_switched_rc.cir` (a switched-RC LPTV circuit with a clocked
`SW` switch -- chosen first because it is fast, 0.13s/run, and has a clean,
literal source-edge breakpoint to check against).

**0200 -- the checkpoint-restore path dropped the "we just landed on a
breakpoint" breakpoint.** `dctran.c`'s checkpoint branch filters the
restored `CKTbreaks[]` table with `CKTbreaks[nb] > CKTtime + CKTdelmin`,
intending to drop stale past breakpoints. When the checkpoint is taken
*exactly at* a breakpoint (e.g. a clocked source's edge -- exactly what
happens if you checkpoint a converter right after it settles onto its own
switching grid), this also drops the *current* breakpoint, which an
in-memory continuation would still see as `CKTbreaks[0]` (`CKTaccept` only
clears a breakpoint once `CKTtime > CKTbreaks[0]`, not `==`). Losing it
means the `AlmostEqualUlps(CKTtime, CKTbreaks[0])` check that cuts
integration order to 1 and shrinks the next step right at an edge never
fires after a `loadstate`, so the restarted run takes one large step
straight through the edge instead of the many small ones an uninterrupted
run takes. **Measured before the fix**: checkpointing R2 at t1=300us (a
multiple of its 10us clock period, landing exactly on an edge), `v(vx)`
diverged from the uninterrupted run by **~7.7e-5 relative right at t1**
(interpolated probe 0.5ns after resume; a directly-sampled 300us+1ns point
showed the restarted run had already skipped straight to `clk=1.0` while
the uninterrupted run was still resolving the edge's 1ns rise time in
several sub-nanosecond steps).

**0201 -- fixing 0200 naively created a divide-by-zero-shaped bug.** The
straightforward fix (keep breakpoints `>= CKTtime - CKTdelmin` instead of
strictly excluding the current one) exposed a second issue: the saved
breakpoint table can hold *two* entries that both equal (or are within
`CKTdelmin` of) the resume time -- e.g. a clocked source's own edge
breakpoint *and* the checkpointed sub-run's own `.tran`-command final-time
breakpoint, which coincide exactly whenever the checkpoint happens to be
taken at a source edge (the case 0200 was built to fix). Keeping both
verbatim left `CKTbreaks[1] == CKTbreaks[0]`, so `.1 * MIN(CKTsaveDelta,
CKTbreaks[1]-CKTbreaks[0])` evaluated to exactly 0 and the very first step
after resume failed outright: `doAnalyses: TRAN: Timestep too small; time =
0.0003, timestep = 0`. Fixed by de-duplicating adjacent kept breakpoints
within `CKTdelmin` of each other.

**After both fixes**, the same R2 checkpoint-at-an-edge case reproduces the
uninterrupted run **bit-for-bit**: 0.0 relative error on `v(vin)`, `v(vx)`
and `v(clk)` across all 12947 samples the two runs share from t1 to t2 (an
exact, index-for-index, time-for-time match -- not a tolerance band).

## Restart exactness, per deck (`tests/test_checkpoint.py`)

Ran with the fork build `aaba50ba3f49-ckpt5` (pin `aaba50ba` + `aug`
patches 0001-0016 + `patches-wip/ckpt/0200-0201`). Every case: `savestate`
at t1 (its own process), `loadstate` to t2 (a fresh process), compared
against a **true single-shot uninterrupted `tran 0..t2`** (a third,
independent process) -- the literal reading of the spec's "matches an
uninterrupted run".

| deck | t1 | t2 | checkpoint alignment | result |
|---|---|---|---|---|
| `refs/r2_switched_rc.cir` | 300us | 600us | exactly on a clock edge (300us = 30 x 10us period) | **bit-exact, 0.0 relative error**, 12947/12947 samples grid-matched |
| `bench/decks/d_analog_block.cir` (BSIM4 x300, no `uic`) | 1000ns | 2000ns | arbitrary (mid-plateau, no source edge there) | **bit-exact to the resolution of a 500-point comparison grid**, 0.0 relative error |
| `refs/r3_vm_buck.cir` (voltage-mode buck, tanh switches) | 500us | 1500us | on a PWM ramp reset (500us = 250 x 2us period) | `v(vout)` (the smooth, regulated state) agrees to **~1.4e-4 relative** over the ~500 further switching periods to t2 |
| `refs/r3d_vm_buck_deadtime.cir` (+ dead-time network, body diode) | 500us | 1500us | same | `v(vout)` agrees to **~8.3e-5 relative** over ~500 further periods |

**Why R3/R3d are ~1e-4, not 1e-9, and why that is not the same bug as
0200/0201**: floating-point rounding is not associative. A run split into
"integrate to t1, serialise every byte, deserialise, integrate to t2" does
not walk the *bit-identical* sequence of Newton iterations a single
continuous run does, once hundreds of *further, independent* breakpoint-
timing decisions follow the split (R3/R3d switch every 2us; 500 more edges
happen between t1 and t2). A sub-ULP difference at any one of those 500
independent decisions can tip a step-size choice a few picoseconds either
way, and that is enough to shift a fast comparator's crossing instant by a
similarly tiny amount -- not a wrong answer, a *different but equally
valid* adaptive step sequence through the same DAE. This is not evidence
the restored state is wrong: R2 and BSIM4 above restore the *identical*
byte patterns and, precisely because neither has hundreds of *further*
independent switching decisions after t1 amplifying that non-associativity,
reproduce the uninterrupted run bit-for-bit. `v(pwm)`/`i(L1)` (fast,
ripple-bearing signals) are excluded from the tight R3/R3d comparison for
the same reason interpolating a near-binary or ripple-phase-sensitive
signal at an arbitrary sample time inflates a picosecond-scale timing
difference into an apparently large "relative error" at that one sample --
confirmed a comparison-methodology artifact, not a real divergence, by
checking it disappears when comparing at the *same* accepted timepoints
instead of interpolating (done for R2's `v(clk)`, which showed the same
symptom and the same 0.0 truth once compared correctly).

## Red-proof

`tests/test_checkpoint.py`'s `red_proof_breakpoint_bug` re-runs the exact
R2-at-an-edge case above against the **unpatched** upstream build
(`aaba50ba3f49-ckpt1`, kept locally alongside `-ckpt5` specifically so this
check can run live) and requires it to show more than 1e-6 relative error
-- it measures **2.3e-5**, comfortably over the bar, confirming this file's
own exactness check is sensitive and not a tautology.

A second, deliberately-broken build was made (not kept, not patched in --
`com_checkpoint.c`'s `loadstate` edited to read `CKTorder`/`CKTdeltaOld`
from the file but then force `CKTorder=1` and reset `CKTdeltaOld[]` to
`CKTmaxStep` instead of applying them -- i.e. "drop the integration-history
restore" per this card's own instructions) to check that dropping the
history restore is *also* visibly caught, separately from the breakpoint
bug 0200/0201 fixes. Result: on R2, this broke restart **outright and
loudly** -- `doAnalyses: TRAN: Timestep too small; time = 0.0003, timestep
= 0` on the very first step after `loadstate`, not a subtle wrong number.
Not kept as a patch (it was never meant to be one); recorded here as the
red-proof evidence the task asked for.

## Using it directly

```
savestate <file>      # in a running .tran, after e.g. `tran 10n 500u 0 uic`
```
then, in a **fresh** ngspice process, with the *same* deck (same build,
same options) and a top-level `.tran <tstep> <t2> ...` card (its `tstop`
becomes the new end time -- `loadstate` needs a `.tran` job to read that
from, per `com_checkpoint.c`'s `find_tran_job()`):
```
loadstate <file>       # runs the continuation to <t2> and returns
```
**Trap, found while building the FRA integration below and not documented
anywhere in the ngspice manual**: ngspice's control-mode tokenizer
lower-cases bare (unquoted) command arguments, so `savestate Foo.bin` /
`loadstate Foo.bin` silently read and write `foo.bin`. This is invisible
until the file "isn't found" -- and it is *fatal* for this project's own
`~/Work/data/ngspice_augmentation/...` convention (a capital `W`) if a
checkpoint path is built from it directly. Use a bare, already-lower-case
relative filename plus the ngspice process's own working directory
(`cwd=`) instead of an absolute path under `~/Work/...`; never build the
argument from a mixed-case path.

## FRA/PSS-from-checkpoint (task item 3): `tools/fra.py --checkpoint`

`tools/fra.py`'s mode "loop" (the transient Tian/Middlebrook FRA harness,
card #452) re-settles the converter from t=0 through its own `op`/`uic`
startup **for every one of its three sub-transients (unperturbed,
v-injection, i-injection), for every frequency point** -- 3N full settles
for an N-point sweep, even though the pre-injection trajectory (t=0 to
`td`, the injection turn-on time) is *identical* for every frequency and
every sub-run by construction (the injection SIN sources are 0 before
`td`, module docstring). Opt-in checkpoint mode (`--checkpoint`, needs
`NGSPICE` pointed at a `savestate`/`loadstate`-capable fork build):
settle **once** to `td`, `savestate`, then run every sub-transient of every
frequency point as a `loadstate` continuation from that one file instead
of a fresh `tran ... uic` from t=0. New functions: `settle_checkpoint`,
`run_loop_point_checkpoint`, `sweep_loop_checkpoint`,
`window_times_checkpoint` (all additive; `run_loop_point`/`sweep_loop`/
every existing default are completely unchanged -- `tests/test_fra.py`,
33/33, and `tests/test_fra_converter.py` still pass unmodified).

Measured (`tests/test_checkpoint.py`'s `fra_checkpoint_case`, same fork
build, `.options trtol=1`, `settle_periods=3`/`measure_periods=5`):

| deck | settle time | freq points | plain `sweep_loop` | `--checkpoint` | speedup | agreement |
|---|---|---|---|---|---|---|
| R3 (voltage-mode buck) | 300us | 2 | 2.85-2.96s | 1.99-2.04s | **~1.4-1.5x** | worst 0.069 dB / 0.772 deg |
| R4 D~=0.30 (peak-current-mode buck, the task's own "3.5ms/point" case) | 3.5ms | 1 | 16.79-17.05s | 5.46-5.56s | **~3.0-3.1x** | worst 0.004 dB / 0.055 deg |

Agreement is well inside the FRA harness's own noise bounds at these
settings (`test_fra_converter.py` uses 0.1dB/1deg at `settle_periods=3` on
a much larger frequency set; this check's looser 0.5dB/5deg bound exists
only to separate "checkpoint mode agrees with plain mode" from "the FRA
itself is noisy at short settle_periods", already characterised
separately). Not bit-identical (same non-associativity discussion as
R3/R3d above -- both modes see hundreds of switching periods before the
measurement window), but well within the ill-conditioning bounds
`tools/fra.py`'s own `combine_tian_middlebrook` already manages.

The speedup at 1-2 frequency points already exceeds the settle-time
fraction alone would predict (R4: 3.5ms settle vs. a ~2ms measurement
window per point should give at most ~1.7x on the *dominant* cost, not
3x) -- consistent with the settle phase being redone 3 times per point in
plain mode (base/v-injection/i-injection) and once total in checkpoint
mode; a full N-point sweep (this card did not run one to keep total
runtime bounded, see "what was not done") should show a larger win still,
approaching the theoretical `(1 + 3N)/(1 + N)` factor as N grows, capped
by the per-frequency measurement cost itself becoming the floor.

## Adaptive settle (card #487): measuring `td` instead of guessing it

`tools/fra.py`'s fixed `--settle` (above, and the FRA/PSS-from-checkpoint section)
is a human guess at how long the UNPERTURBED circuit needs to reach its own
periodic steady state before injection turns on. Card #484 found the guesses this
project had been using too short for slow converters: R3's 300us guess left a
6.15% periodicity residual (`max|v(t)-v(t-T)|` over one switching period `T`,
relative to that window's own peak-to-peak ripple — same convention
`tools/pss_ref.py`'s `ripple()`/`settle_period_index()` and `tools/pssac_ref.py`'s
`_measure_h_tran` periodicity check already use elsewhere in this project), biasing
R3's own FRA curve by up to 3.7dB/17deg near fs/2; 1ms left ~0.6%; R3 actually
needs roughly 9ms to reach <=1%.

`--adaptive-settle` (opt-in; the old fixed `--settle` path is completely unchanged
and stays every existing test's default) replaces the guess with a measurement:
`periodicity_scan` runs the circuit ONCE, unperturbed, from t=0 to a cap, samples
the same nodes the sweep itself will measure (mode loop: `v(<probe>.x)` and the
probe's `a` node; mode transfer: `output_expr`), and finds the earliest multiple of
the caller-given natural period `--settle-period` (e.g. `1/fs` for a switching
converter — NOT an injection period) at which the periodicity residual is at or
below `--residual-bound` (default 1%). That ONE measured settle time is then reused
for every frequency point in the sweep, exactly like a fixed `--settle` value would
be — the number of re-settles per sweep is unchanged, only `td` is measured rather
than guessed. `--settle-cap` bounds the search (default 50ms); if the bound is
never reached, the scan honestly reports `capped=True` with the best residual
actually seen — never a fabricated pass. Combine with `--checkpoint` (mode loop
only) to also pay that settle cost once total instead of once per sub-transient/
frequency point (same mechanism as this file's "FRA/PSS-from-checkpoint" section
above — `sweep_loop_adaptive(..., checkpoint=True)` calls `sweep_loop_checkpoint`
with the measured `td`).

Every CSV row (`fra.csv`, ok or could-not-run alike) and `summary.json` carry the
one settle_time/residual/capped this sweep used (`settle_time_s`, `settle_residual`,
`settle_capped` columns; a `settle` key in the JSON) whenever `--adaptive-settle`
was used; the old fixed-settle CLI path produces the byte-identical pre-#487
CSV/JSON shape (checked directly, `tests/test_fra_settle.py`).

    tools/fra.py DECK --probe XP --freqs "dec 5 500 250000" --amplitude 0.01 \
        --adaptive-settle --settle-period 2e-6 --residual-bound 0.01 \
        --settle-cap 0.05 --startup uic --tmax 5e-9 \
        --spice-options ".options trtol=1" [--checkpoint]

Measured (`tests/test_fra_settle.py`, R2 switched-RC, fclk=100kHz/period=10us, fast
deck): adaptive settle at a 1% bound found `td=60us` (6 periods, residual 0.91%),
well under the 300us (30-period) fixed settle `tests/test_fra_lptv.py` already
validated as converged for this circuit — the resulting FRA point (10kHz,
control-to-`v(vx)`) agrees with the fixed-300us-settle result to <0.0001dB/0.0001deg
(same trajectory, just measured instead of over-settled).

R3d/R4d (card #479, the slow-converter case this card was written for — R3 itself
needs ~9ms, card #487's own report): tracking `v(vout)` PLUS the loop probe's own
`v(xp.x)`/`v(fbtop)` PLUS the compensator's own output node (`v(opout)` for R3d,
`v(vc)` for R4d) at a 1% bound, `tools/fra.py --adaptive-settle` found R3d
`td=510us` (255 periods, residual 0.97%) and R4d `td=1672us` (836 periods, residual
0.98%) — both comfortably reused as one fixed settle time for every point of their
own 14-point loop-gain and ctrl2out/line2out/zout FRA curves (regenerated under
`~/Work/data/ngspice_augmentation/fra/curves-fork-s37-rd/`, card #479 part 2). The
resulting R3d/R4d loop-gain curves agree closely with R3/R4's own long-published
tables (e.g. R4d @ 50kHz: -18.988dB/-101.43deg vs R4's -18.987dB/-101.317deg) — a
strong physical sanity check, since both share the same R1..C3/type-II compensator
values and dead-time is expected to be a second-order effect at these frequencies
(the decks' own header comments). **Tracking `v(vout)` alone is not safe**: a first
attempt on R3d tracking only `v(vout)` found an implausibly short `td=256us`
(0.89% residual) — root cause, confirmed by comparison: `vout` is tightly regulated
by the loop and reads periodic almost immediately, even while the compensator's own
slower internal state is still settling (correcting exactly that drift is what a
loop does), so `vout` alone is not a safe proxy for "has the WHOLE orbit reached
periodic steady state" — the node list above (adding the loop-probe and compensator
nodes) is what this project's own `tools/fra_curves_fork.py`-style regeneration
scripts should track for a converter with feedback, not `vout` alone.

New functions in `tools/fra.py` (all additive; `run_loop_point`/`sweep_loop`/
`run_transfer_point`/`sweep_transfer` and every existing default are unchanged):
`_interp_linear`, `periodicity_residual`, `periodicity_scan`, `sweep_loop_adaptive`,
`sweep_transfer_adaptive`. `write_loop_results` gained an optional `settle_info`
parameter (default `None`, reproduces the old CSV/JSON shape exactly).

### What adaptive settle does NOT do

- It does not change the FRA measurement itself (`run_loop_point`/
  `run_transfer_point`, the DFT extraction, Tian/Middlebrook combination) — only
  which `settle_time`/`td` value feeds into it.
- It does not retry with a longer cap automatically if `capped=True` — the caller
  (or the CLI's own printed line, `adaptive settle: ... CAPPED`) has to notice and
  raise `--settle-cap`; a silently-widening cap risks masking a circuit that
  genuinely never reaches a periodic steady state (e.g. a subharmonic-unstable
  converter, refs/README.md's R4 D=0.55 case) as "just needs more time".
- The periodicity-residual scan's own sampling (`--settle-period` /
  `points_per_period`/`tmax`/`spice_options`, all forwarded to the scan) should
  match the accuracy the caller actually needs from `td` — a converter deck should
  pass the same `tmax`/`.options trtol=...` the measurement sweep itself uses (see
  the "Timestep-control fix" note in `tools/fra_curves.py`), or the scan can settle
  on a `td` that looks periodic at a coarser resolution than the measurement will
  actually run at.

## What was not done

- XSPICE event-driven / digital code-model state across a checkpoint is
  identified as a real gap (com_checkpoint.c saves none of it) but not
  tested -- none of this card's four target decks use XSPICE event models.
- No dedicated red-proof for `CKTbreaks` corruption beyond the two bugs
  found; other possible ways to break restart (e.g. dropping `CKTrhsOld`
  entirely) were not separately red-proofed -- the two real bugs found by
  characterisation were fixed and red-proofed; a from-scratch fault-
  injection sweep over every saved field was not attempted.
- `tools/fra.py --checkpoint` was validated at 1-2 frequency points per
  deck (kept small to bound `tests/test_checkpoint.py`'s own runtime,
  ~60s total); a full N-point sweep (`tools/fra_curves.py`'s own 14-point
  list) was not re-run in checkpoint mode to get a first-hand full-sweep
  wall-clock number -- the 1-3 point numbers above are real, measured, not
  extrapolated, but a full-sweep number is extrapolated from them, not
  measured.
- `run_loop_point_checkpoint`'s `linearize` call pads its output from
  absolute t=0 rather than from the checkpoint's own first real data point
  (confirmed: the resulting file has `tstop/dt + 1` rows regardless of how
  large `td` is, with garbage/extrapolated values before the real data
  starts) -- correctness is unaffected (the extraction window is always
  well after `td`, so the garbage rows are never read), but it is a real,
  unoptimised inefficiency on a checkpoint mode whose whole point is
  speed; not fixed here (would need either a smaller `linearize` range,
  which ngspice's `linearize` does not appear to support, or dropping
  `linearize`/`check_uniform` for the checkpoint path specifically and
  relying on `extract_component`'s own already-non-uniform-safe
  trapezoidal integration instead -- a design change to fra.py beyond
  what "opt-in mode, keep defaults passing" called for here).
- No PSS-restarted-from-checkpoint campaign was built (the task's item 3
  names "an FRA/PSS campaign" in the spec line; only the FRA half was
  built, since `patches-wip/pss`'s shooting-Newton PSS is a separate,
  in-progress line of work this card does not touch).
- Checkpoint files were not tested across machines or ngspice builds
  (the format is explicitly native/non-portable per the code and this
  doc's own "File format" section; not independently re-verified here
  beyond reading the code).
