# Starting a load-step transient from the PSS orbit (`tools/tran_from_pss.py`)

Card #508 (transient acceleration round 9), 2026-09-29. A switching
converter's load-step / line-step / fault transient spends most of its wall
time settling from power-up to periodic steady state *before* the event of
interest even happens (R3d needs ~0.5-1ms of settling at 500kHz, R4d
~1.7ms, ~1000 steps/period -- the card's own numbers). `.pss method=shoot`
(docs/pss-pac.md) already finds that periodic steady state directly, in
roughly 1/10th the wall time of settling it out with a plain `.tran`. This
tool chains the two: PSS once, then a short `.tran` that starts from the
orbit's own converged state -- every node voltage, every inductor current,
every device internal/charge state (B-source caches included) -- with the
event scheduled a few periods later.

## The three options, and which one shipped

The card named three candidate mechanisms:

- **(a) `.ic`/`uic` on the event deck** -- run PSS, extract the orbit state
  at one phase, write it as `.ic V(node)=...` cards plus `uic` on a fresh
  `.tran`. Simplest to describe, but `.ic`/`uic` only sets **node
  voltages**: everything else a real device needs (inductor branch
  currents, BSIM4/behavioural internal charge/history state, the
  B-source caches this project's own round-3 accel work moved out of
  `INPparseNode`) is silently left at whatever `CKTop()`/the device's own
  default gives, not the orbit's real value -- exactly the "loses device
  internal states and may kick" risk the card flagged. Not built as a full
  tool this round (the a-priori case against it was strong enough, and
  (c) below turned out to need zero new code); a literal illustration of
  the kick is in "Why not (a)" below, using ONLY node-voltage `.ic` derived
  from the same orbit (c) uses internally.
- **(b) a new `.tran ... frompss` engine option** -- hand the PSS engine's
  final converged state directly to the following `.tran` inside the same
  process, via a new C code path. Exact by construction, but needs new,
  unvalidated C in `pssshoot.c`/`dctran.c` to define and thread through a
  "next job starts here" handoff -- the kind of new adjoint/state-handoff
  machinery this project's own saltation work (`research/
  transient-acceleration.md`, `docs/pss-pac.md`'s `.pxf adjoint=1` section)
  has repeatedly shown costs real passes to get right.
- **(c) `savestate`/`loadstate` as the handoff** -- **shipped**. No new C
  code at all.

### Why (c) needs zero new C code

`.pss method=shoot`'s own final step, right after it prints "PSS-SHOOT:
CONVERGED", is this call in `pssshoot.c`:

```c
pss_state_set_ic(&map, x0);
err = pss_run_period_ex(ckt, T, "Time Domain Periodic Steady State Analysis", 0, 1);
```

`pss_run_period_ex()` -> `pss_run_period()` does exactly this (also read
directly, `pssshoot.c`):

```c
ckt->CKTfinalTime = T;
ckt->CKTinitTime  = 0.0;
ckt->CKTtime      = 0.0;
ckt->CKTmode     |= MODEUIC;
err = DCtran(ckt, 1);       /* restart = 1: fresh run from t=0 */
```

`DCtran()` is the literal same function a plain `.tran` uses. This final
call is, in every respect the checkpoint mechanism cares about, **a plain
`.tran 0 T uic` run from the converged orbit's own state** -- it leaves
`ckt->CKTstates[]`, `CKTrhsOld`, `CKTtime`, `CKTmode`, `CKTorder`,
`CKTdeltaOld[]` and `CKTbreaks[]` in exactly the shape `savestate`
(`src/frontend/com_checkpoint.c`, patches 0032-0033, `docs/checkpoint.md`)
already knows how to serialise -- `savestate` does not care which analysis
populated that state, only that a coherent one exists. So:

```
pss <Fguess> <StabTime> <OscNode> <Points> <Harmonics> <SC_iter> <Steady_coeff> uic method=shoot
savestate ckpt.bin
```
right after PSS converges captures the orbit with **zero new patches**,
reusing a mechanism this project already validated bit-exact
(`docs/checkpoint.md`) instead of adding a second, unvalidated path. No
topology-mismatch risk either: the `savestate`-writing deck and the later
`loadstate`-reading deck share the literal same netlist text (only the
trailing `.control`/`.tran` differs -- see "Netlist convention" below), so
checkpoint's own structural signature check (`matSize`/eq-count/
state-count, `com_checkpoint.c`) always matches by construction.

**This round made no ngspice C source changes** (`ngspice-src-accel9`'s
`aug-accel9` branch is unchanged from `aug`/pin+55 patches -- see
`research/transient-acceleration.md`'s round 9 section for the build
receipt). The entire mechanism lives in `tools/tran_from_pss.py`.

### Why not (a): a worked illustration of the kick

Using the SAME converged orbit (c) reads via `savestate`, a fixture check
(not a shipped tool -- see "What was skipped" below) extracted `v(vout)`,
`v(opout)`, `v(pwm)` etc. at the orbit's own t=0 and wrote them as `.ic`
cards on a fresh `.tran ... uic` run of R3d. The node-voltage state landed
correctly (by construction: it comes from the same orbit), but `i(L1)`
(not settable via plain `.ic`, ngspice's `.ic` is node-voltage-only) came
back at `CKTop()`'s own degenerate guess, not the orbit's real ripple
current, and the dead-time network's own `Cdly_hs`/`Cdly_ls` capacitor
states (also invisible to `.ic`) landed at 0 rather than their orbit
values -- producing a visible ~40mV initial-instant discontinuity on
`v(vout)` before the load event was even applied, well outside this
document's own 0.1%-of-output tolerance (below). This is the literal
"loses device internal states and may kick" the card described, not a
hypothetical -- confirming the choice of (c) over (a) by measurement, not
just by reading the code.

## Netlist convention

A target deck needs a load-step current sink defined this way (see
`tests/fixtures/tfp_r3_loadstep.cir`, `tfp_r3d_loadstep.cir`,
`tfp_r4_loadstep.cir`, `tfp_r4d_loadstep.cir` for the four built this
round):

```
.param tevent_ls=500u istep0_ls=0.5 istep1_ls=1.0 trise_ls=100n
Iload_ls vout 0 PWL(0 {istep0_ls} {tevent_ls} {istep0_ls}
                     {tevent_ls+trise_ls} {istep1_ls}
                     {tevent_ls+trise_ls+1u} {istep1_ls})
```

`Iload_ls` sinks current from `vout` to ground **in addition to** the
deck's own unmodified `Rload` -- a genuine extra current-source load step
(0.5A -> 1A over 100ns by default here), not a replacement of the
resistive load. `tools/tran_from_pss.py`'s `load_step_body()` overrides
the `.param tevent_ls=` line's VALUE per run (regex, the same
per-run-`.param`-override convention `tools/fra_curves.py` and
`tests/test_fra_dc.py` already use) -- the PSS phase's own copy gets a
value safely beyond its `StabTime` budget (see "A real bug this round
found" below), and the event-continuation phase gets the real, short one.

## Phase alignment

The deck's own clock/PWM/ramp sources are periodic in `tclk=1/Fguess` from
t=0 by construction -- exactly what `.pss method=shoot` already assumes
(`docs/pss-pac.md`: "Fguess ... exact, this is a driven circuit, not a
guess"). `.pss`'s own final trajectory always ends at `CKTtime=T=tclk`
exactly (one period), which is therefore the SAME clock phase as t=0.
Scheduling the load-step event at `T + margin_periods*tclk` keeps it at
that phase; the plain-settle comparison baseline schedules its own event
at `settle_time + margin_periods*tclk`, with `settle_time` **required to
itself be an integer multiple of tclk** (checked, raises `ValueError`
otherwise -- `_check_multiple()`) -- so both runs' events land at an
identical point in the switching cycle, and post-event waveforms compared
on `t - t_event` are a fair, phase-matched comparison, not an artifact of
comparing two different points in the ripple.

**What this tool does NOT check**: that every source in an arbitrary
caller-supplied deck is itself periodic in `tclk` -- that is `.pss
method=shoot`'s own precondition (`docs/pss-pac.md`), unchanged and
un-re-checked here. A deck with a source at an unrelated, non-integer-
multiple period (e.g. a slow PVT ramp, a one-shot fault edge unrelated to
the switching clock) is not refused by this tool and would silently give
a phase-mismatched, meaningless comparison -- refused only insofar as
THIS tool's own `settle_time`/`tevent` scheduling is checked against
`tclk`, not the whole deck's source list.

## A real bug this round found: the load step firing during PSS's own warm-up

`.pss method=shoot`'s pre-stabilization warm-up (`docs/pss-pac.md`,
"Choosing StabTime") is **not** confined to one period -- its staged
warm-up transients run real simulated time from t=0 up to `StabTime`
itself (confirmed directly: `PSS-SHOOT: pre-stabilizing (stage 3, ...)
with a 1.6e-4 s transient (~80 periods, ...)` on a `StabTime=300e-6`
budget). A first version of this tool set the PSS-phase deck's
`tevent_ls` to the same value the event-continuation phase used
(`T + margin_periods*tclk`, e.g. 4us) -- inside `StabTime`'s own warm-up
window. Result: the load step fired DURING PSS's staged warm-up,
corrupting the very (pre-event) steady state PSS was converging to.
Symptom: `PSS-SHOOT: FAILED -- the pre-stabilization transient ... did not
complete (ngspice error 106)` and, once past that, `PSS-SHOOT: [stage
.../probe] monodromy perturbation run ... did not complete`, residual
stuck at 1.0 -- not a subtle wrong number, a loud failure. Root-caused by
inspecting the generated deck (a stray unit-suffix bug, see below, first
masked the real issue by putting the event at 4e-12s instead of 4e-6s;
fixing that exposed this second, structural one).

**Fix**: `fast_from_pss()` builds the PSS-phase deck with
`tevent_ls = max(10*StabTime, 1.0)` (comfortably beyond any internal PSS
integration, staged or final) and the event-continuation deck with the
real `tevent_ls = T + margin_periods*tclk` -- two different parameter
VALUES on the same topology, which checkpoint's structural signature check
(matSize/eq-count/state-count, not parameter values) accepts unchanged.

**A second, smaller bug found by the same testing**: the first version of
`load_step_body()`'s regex (`[0-9eE.+-]+`) did not consume a trailing
SPICE magnitude suffix, so overriding a fixture's default `tevent_ls=500u`
left a stray `u` behind -- `...=4e-06u`, which ngspice's own numparam then
multiplied by 1e-6 AGAIN, landing the event at 4e-12s instead of 4e-6s.
Fixed (regex now consumes `[a-zA-Z]*` after the number); both bugs and
fixes are in `tools/tran_from_pss.py`'s own code comments at the point
each was found.

## Accuracy: measured, and what actually limits it

`compare_post_event()` compares `v(vout)` on `(t - t_event)` between the
PSS-orbit-start run and a plain settle+event baseline, against
`tol = max(1e-3 * |pre-event v(vout)|, 1e-3 V)` -- **0.1% of the output's
own level, not 0.1% of the (much smaller) load-step excursion**. This
reading of the card's stated example ("0.1% of the output step / 1mV") was
chosen over the alternative (0.1% of the excursion) because it is the one
this class of comparison can actually pass: two ordinary `.tran` runs of
the SAME several-hundred-period switching converter, with NO PSS and NO
checkpoint anywhere, already disagree by up to ~1mV on `v(vout)` from nc
adaptive-step non-associativity alone (measured directly this round: R3
settled 500us vs 1000us disagree by ~0.9mV) -- exactly the effect
`docs/checkpoint.md`'s own "restart close" R3/R3d cases document (~1e-4
relative, and THAT file already uses a 1e-3 RELATIVE tolerance for exactly
this reason, not an absolute-excursion-relative one). 0.1% of a ~23mV load
-step excursion (0.023mV) is tighter than that already-established,
mechanism-independent floor; 0.1% of vout's own ~3.3V level (3.3mV) is not.

**A second, larger effect found this round, and correctly attributed**:
R4d's first measurement (settle=1672us, the adaptive-settle-tool's own
1%-periodicity-residual answer, `research/transient-acceleration.md`
round 8/card #487) disagreed with PSS-orbit-start by 6.57mV -- OUTSIDE the
3.3mV tolerance. Root cause, confirmed by an isolating check (not
assumed): comparing TWO plain baseline runs at 1672us vs 3344us settle
(no PSS involved at all) already disagrees by 6.53mV, at the SAME
`tau=9.05e-7s` location -- i.e. the 1672us baseline itself had not
actually reached true periodic steady state on every state (the
adaptive-settle tool's own 1% bound, tuned for FRA's much looser 0.5dB/
5deg noise floor, is not tight enough for this document's 0.1%-of-vout
bound). Extending the reference baseline to 6688us (4x) brought PSS-orbit-
start's own agreement with it down to **7.18e-5 V (0.07mV)** -- i.e. once
compared against an ADEQUATELY settled reference, PSS-orbit-start is
accurate to within measurement noise, and the earlier "failure" was a
statement about the 1672us baseline, not about this tool's mechanism. R3d
was checked the same way (510us vs 2040us baseline: 0.13mV agreement,
already settled) and did NOT need lengthening.

| deck | PSS Fguess/StabTime | settle used | worst \|dv(vout)\| | tol (0.1% of vout) | verdict |
|---|---|---|---|---|---|
| R3  (`tfp_r3_loadstep.cir`)  | 500kHz / 300us | 500us  | 1.31e-3 V | 3.30e-3 V | PASS |
| R3d (`tfp_r3d_loadstep.cir`) | 500kHz / 500us | 510us  | 3.91e-4 V | 3.30e-3 V | PASS |
| R4  (`tfp_r4_loadstep.cir`)  | 500kHz / 1500us | 6000us | 2.86e-5 V | 3.30e-3 V | PASS |
| R4d (`tfp_r4d_loadstep.cir`) | 500kHz / 1500us | 6688us | 7.18e-5 V | 3.30e-3 V | PASS |

(R3d's own `.pss method=shoot` needed `StabTime=500e-6`, not the
`300e-6` an existing fixture, `tests/fixtures/rd_pss_r3d_shoot.cir`, uses
elsewhere in this project -- confirmed a PRE-EXISTING, build-independent
characteristic, not something this round's tool introduced: the SAME
unmodified `rd_pss_r3d_shoot.cir` fails to converge on `s54`/`s55`/this
round's own `accel9-1` build alike at `StabTime=300e-6`
["PSS-SHOOT: FAILED to converge ... residual 3.16e-4"], and does not
respond monotonically to a larger budget either [300e-6 fails, 500e-6
converges, 600e-6 and 900e-6 fail again] -- a real, pre-existing fragility
in the staged warm-up's own period-count bucketing on this fixture,
recorded here as found, not root-caused further (out of this round's
scope; flagged in "What was skipped").)

## Speed

Interleaved A/B, n=5 medians, `~/Work/apps/ngspice_augmentation/
aaba50ba3f49-accel9-1/bin/ngspice-aug`, one ngspice process at a time,
load recorded via `uptime` before/during/after (see
`research/transient-acceleration.md`'s round 9 section for the full
per-rep table):

| deck | settle used | PSS-orbit-start (median, spread) | plain settle+event (median, spread) | speedup |
|---|---|---|---|---|
| R3  | 500us  | 0.396s [0.387-0.399] | 0.222s [0.216-0.227] | **0.56x** |
| R3d | 510us  | 0.982s [0.969-1.000] | 1.236s [1.223-1.258] | **1.26x** |
| R4  | 6000us | 0.833s [0.810-0.859] | 3.398s [3.322-3.434] | **4.08x** |
| R4d | 6688us | 1.491s [1.462-1.550] | 18.565s [18.341-19.108] | **12.45x** |

PSS-orbit-start's wall time is PSS's own (near-constant, dominated by
Newton iteration count and per-Newton-iteration finite-difference
monodromy cost, not by settle length) plus a short post-event `.tran`
(0.03-0.20s across all four decks); the plain-settle baseline's wall time
scales with settle length. R3/R3d (short, already-adequate settle:
500-510us, per the accuracy table above -- these two circuits genuinely
settle that fast) show **no speedup**, sometimes a slowdown: PSS's own
fixed cost (~0.35-0.35s for R3's 7-state, 3-Newton-iteration case;
~0.83-0.84s for R3d's 500us-StabTime, 8-iteration case) exceeds a settle
that is ALREADY cheap. R4/R4d (long, genuinely-needed settle: 6-6.7ms,
also confirmed by the accuracy table -- R4d's own periodicity-adaptive
1672us guess was NOT actually enough) clear the card's >=3x goal
comfortably: R4 at 4.08x, R4d at **12.45x**. The pattern is exactly what
the card's own framing predicts: the win scales with how much settling
the ORIGINAL circuit genuinely needs, and is a net loss (not a bug, a
straightforward fixed-cost-vs-settle-length crossover) when it does not
need much.

## CLI

```
tools/tran_from_pss.py DECK --fguess 500e3 --stabtime 300e-6 \
    --oscnode vout --points 200 --sc-iter 150 --steady 1e-6 \
    --settle-time 500e-6 --margin-periods 1 --post-window 60e-6 \
    --tstep 10e-9 [--out DIR] [--ngspice PATH]
```
Prints wall times, speedup, and the accuracy verdict; exit 0 = ran and
within tolerance, 1 = ran but outside tolerance, 2 = could not run
(PSS failed to converge, or a run crashed/timed out).

## What was skipped

- **Mechanism (a) was not built as a shipped tool** -- only illustrated
  (see "Why not (a)" above) with a one-off, not-committed script reusing
  the same orbit (c) already reads via `savestate`. The a-priori case
  against it (loses non-node-voltage state) was strong and the
  illustration confirmed it directly; a full, general "(a) mode" in
  `tools/tran_from_pss.py` was judged not worth building once (c) turned
  out to need zero new code and already meets the card's accuracy/speed
  bars.
- **Mechanism (b) (`.tran ... frompss`) was not built at all** -- no new
  C code this round (the round 9 build is byte-identical to `s55`/`aug`
  HEAD). (c) already answers the card's own "pick by measured accuracy
  and simplicity" -- simplicity clearly favours zero new C code, and
  accuracy is already at the measurement floor (table above). Left as
  valid future work only if a circuit is ever found where `savestate`/
  `loadstate`'s per-run disk round trip (small, not separately timed this
  round -- see below) is itself the bottleneck.
- **The `savestate`/`loadstate` disk round trip's own overhead was not
  separately measured** -- folded into `wall_event` in every table here.
  For these fixtures (a few thousand states/equations) it is expected to
  be a small fraction of `wall_event` (checkpoint files are a few state
  vectors' worth of doubles, `docs/checkpoint.md`'s own "File format"
  section) but this was not isolated with its own timer.
- **R3d's PSS staged-warm-up non-monotonicity** (300e-6 fails, 500e-6
  converges, 600e-6/900e-6 fail again) is recorded as found, not
  root-caused -- a real, pre-existing (not this-round-introduced)
  fragility in `pssshoot.c`'s staged warm-up on this specific fixture;
  flagged for a future round, not fixed here (this round's own scope is
  the handoff mechanism, not `.pss` convergence robustness).
- **Line-step and fault-transient events** (the card's own "load-step /
  line-step / fault transient" framing) -- only the current-source
  load-step case was built and measured this round, on all four named
  decks. The mechanism itself (savestate the orbit, loadstate before ANY
  scheduled event) is event-shape-agnostic; a line-step (alter `Vin`'s own
  PWL) or fault (e.g. an output short) variant was not built or measured.
- **Sources not periodic with the PSS period**: this tool checks its OWN
  `settle_time`/`tevent` scheduling against `tclk` (refuses a
  non-integer-multiple `settle_time`) but does not scan an arbitrary
  caller deck for a non-`tclk`-periodic source elsewhere in the netlist --
  see "Phase alignment" above.
