# Loop stability at transient time points (`stb_tran`)

Card #497 (feature-survey #12, `research/feature-survey.md`), 2026-09-29;
**updated by card #500, 2026-09-29: `tools/stb_tran.py` now has a second
mode, `--mode=instant`, that removes this doc's own original limitation**
(the section below, "What 'freeze the operating point at t_k' means here,
and why", proved that native `.stb` could not be made to literally
linearise about a non-equilibrium transient point **without a fork
change** -- card #500 made that fork change, a new `.stb ... instant`
flag; `docs/stb.md`'s "Instant mode" section has the full C-side design
and validation). That original analysis is kept below **unedited** because
it is still exactly correct for `--mode=quasistatic` (now the explicit,
still-default mode name -- what this doc otherwise calls "this tool" below
refers to quasistatic mode unless a passage says otherwise) and because
the reasoning (`CKTop()` always re-solves to the true DC equilibrium) is
what makes `--mode=instant`'s design necessary in the first place. See
"`--mode=instant`: no freezing at all" below for the new mode.

`tools/stb_tran.py`: given a deck with a large-signal transient (LDO
startup, soft-start, a load step) and a native `.stb` probe, compute loop
gain / phase margin / gain margin / crossover frequency at a list of time
points `t_k` along that transient. Needs the private fork build (native
`.stb`, `docs/stb.md`) -- there is no `tools/stb.py` (system-ngspice)
equivalent, since this reuses native `.stb` directly at each frozen point.

```sh
# quasistatic (default): freeze the exogenous --freeze source(s), re-solve
# every OTHER state to a fresh DC equilibrium (see below for why/when).
NGSPICE_AUG=<fork prefix>/bin/ngspice python3 tools/stb_tran.py DECK \
    --probe Vprobe --freeze Vss --times "10u,100u,250u,400u,1000u" \
    --sweep "dec 30 1 1e8" [--check-tran "TSTEP TSTOP" --check-node NODE]

# instant (card #500): no freezing/re-solving -- run the deck's own REAL
# .tran up to each t_k and linearise about whatever it actually converged
# to there. Needs --tstep (the .tran step size); --freeze/--check-tran do
# not apply (refused if given) -- see docs/stb.md's "Instant mode" section.
NGSPICE_AUG=<fork prefix>/bin/ngspice python3 tools/stb_tran.py DECK \
    --probe Vprobe --mode instant --tstep 50n \
    --times "10u,100u,250u,400u,1000u" --sweep "dec 30 1 1e8"
```

Prints a `PM(t) / fc(t) / GM(t)` table; writes `times.csv` and
`summary.json` under `~/Work/data/ngspice_augmentation/stb_tran/<deck>-
<timestamp>/` (never the repository) unless `--out` is given. Optional
`--plot-loopgain` also writes each `t_k`'s full `loopgain.csv` (freq, `T`
re/im, dB, phase) for anyone who wants the Bode trace at one instant, not
just its margin. Exit status: 0 = ran (a `t_k` with no 0 dB crossing is
a legitimate, non-failing result, same convention as `docs/stb.md`'s
`tools/stb.py`/native `.stb`); 2 = could not run.

## `--mode=instant`: no freezing at all

Runs the deck completely unmodified (no `--freeze`, no source rewriting,
no re-solved equilibrium for any state) via `.tran <tstep> {t_k} 0 uic` +
`.stb ... instant` as netlist cards executed together in one process per
`t_k` (`stb_at_point_instant()`) -- see `docs/stb.md`'s "Instant mode"
section for exactly why both must be netlist cards run together (not two
sequential `.control` commands) and the full C-side design/validation.

This is a **strictly different, and for most purposes more faithful,
answer** than quasistatic mode: it reports the small-signal behaviour of
the circuit's OWN real, possibly non-equilibrium, transient state at
`t_k` -- not the behaviour of a hypothetical circuit whose exogenous
source(s) were held at their `t_k` value forever. The two agree once the
transient has settled (both converge to the same equilibrium, see
`docs/stb.md`'s "settled point equals DC `.stb`" check) and can disagree
substantially during a fast ramp/transient relative to the loop's own
bandwidth (`docs/stb.md`'s mid-ramp LDO example: ~15% crossover-frequency
shift, ~1 deg PM, from a real tracking-error effect quasistatic mode
cannot see since it always re-converges to a zero-tracking-error
equilibrium for whatever source value is frozen).

When to use which:

| | `--mode=quasistatic` (default) | `--mode=instant` |
|---|---|---|
| What it computes | small-signal behaviour of a DIFFERENT (hypothetical, source-frozen-forever) circuit at each t_k | small-signal behaviour of the ACTUAL circuit's own real transient state at t_k |
| Needs | `--freeze <source>` (repeatable) | `--tstep` |
| Cost per t_k | one `op` (fast, no time integration) | one `.tran` from 0 to `t_k` (cost grows with `t_k`) |
| Exact when | the loop settles fast vs. the frozen source's own rate of change (checkable via `--check-tran`) | always -- it is the real transient, not an approximation of it (mod the usual "one LTI snapshot" caveat any single-frequency-sweep Bode analysis has, see `docs/stb.md`) |
| A "how well did the approximation hold" diagnostic | yes, `--check-tran`/`--check-node` | not applicable (nothing is being approximated -- refused if given, see below) |

`--freeze` and `--check-tran`/`--check-node` are refused (COULD NOT RUN,
exit 2) if combined with `--mode=instant`: neither applies (there is
nothing to freeze, and nothing to validate the freeze against), and
silently ignoring them would risk the reader thinking a validity check
ran when it did not.

## What "freeze the operating point at t_k" means here, and why

This is the part of the card's own brief this doc has to answer honestly:
**native `.stb` cannot be made to linearise about an arbitrary
NON-equilibrium point using stock, script-layer techniques** (no fork/C
change was made or is needed to reach this conclusion -- it is a real,
verified property of the existing engine). Read `stban.c`/`acan.c`: both
`.ac` and `.stb` call `CKTop()` (the engine's ordinary DC operating-point
solver, full Newton via `NIiter()`) immediately before loading the
small-signal matrix, unconditionally (unless `.option noopac` short-circuits
a genuinely linear circuit). Newton always converges to the circuit's TRUE
nonlinear DC equilibrium given whatever independent-source values are
present at that moment -- nothing about a prior `.tran`, a `savestate`d
checkpoint, a per-device `IC=`, or a `.nodeset` changes what it converges
*to*, only (at most) how many iterations it takes to get there.

**Three independent, reproducible checks confirm this** (all live in
`tests/test_stb_tran.py`'s `check_freeze_evidence()`, all re-run for this
doc against `~/Work/apps/ngspice_augmentation/aaba50ba3f49-s46/bin/ngspice`;
each uses a diode (`D1`, `IS=1e-14 N=1 RS=1`) whose DC operating current is
a sharp, unambiguous discriminator between "the forced/frozen state" and
"the true equilibrium" -- `~2e-2 A` forced vs `4.329920e-04 A` true, a >40x
gap, nothing subtle to threshold):

1. **A plain `op` discards a capacitor's own instance `IC=`.** `C1 out 0 1u
   IC=2` (`out` pinned to 2 V, far from equilibrium) followed by a bare
   `op`: the diode's operating current comes back `4.329920e-04 A` -- the
   TRUE equilibrium, not the current a 2 V forced bias would produce.
2. **`tran ... uic` then `op`, same session, ALSO discards it.** The first
   transient timestep genuinely does honour `IC=2` (that is what UIC means),
   but the very next `op` call in the same control block re-solves from
   scratch and gives the identical `4.329920e-04 A` -- confirming this is
   not a savestate/loadstate-specific gap, it is `CKTop()` itself.
3. **`savestate` a genuinely non-equilibrium transient state, `loadstate`
   it into a FRESH process, then `stb` immediately after: still discarded.**
   `tran 1n 2u 0 uic` (no `IC=`, so `out` is still mid-charge, far from
   equilibrium) then `savestate`; a fresh process's `loadstate` (which
   itself *runs the continuation to the deck's own `.tran` tstop* per
   `docs/checkpoint.md`) leaves the diode at `@d1[id]=4.04e-2 A` -- but the
   `stb` command issued right after that `loadstate` reports the SAME
   `4.329920e-04 A` a plain `op` on the unmodified circuit gives, discarding
   the loaded, genuinely non-equilibrium state entirely.

None of `.ic`, per-device `IC=`, `.nodeset`, or `savestate`/`loadstate` can
therefore make native `.stb` linearise about an arbitrary time-domain
snapshot -- they are all consumed (at most) by the Newton solve's *starting
guess*, and Newton walks past a starting guess to the unique fixed point
regardless (barring genuine circuit bistability, out of scope here).

### The approach this tool actually uses, and its validity condition

Freeze only the deck's genuinely EXOGENOUS, explicitly time-varying
independent sources (`--freeze NAME`, repeatable -- a soft-start reference
ramp, a supply ramp, a PWL load step) at their `t_k` instantaneous value,
evaluated directly from the source's own `PULSE`/`PWL`/`SIN` spec in Python
(`tools/stb_tran.py`'s own small waveform evaluator + `.param` brace-
expression resolver -- no simulation needed for this step), rewritten as a
plain `DC <value>` source; then run an entirely ordinary `op` + native
`.stb` on that frozen deck. Every OTHER node -- including every real
energy-storage state, exactly the states the card's brief asks to be "fixed
at their t_k values" -- is then whatever the circuit's own Newton solve
converges to GIVEN those frozen exogenous inputs, not literally pinned to
its real transient value at `t_k`.

This is the standard **quasi-static / frozen-time linearisation** used
throughout analog and power-electronics practice for exactly this question
(loop margin during a slow soft-start or supply ramp), and it is exact
provided the loop's own natural dynamics settle fast compared to how fast
the frozen source itself is moving. **It is a documented approximation, not
a literal state freeze**, at any `t_k` where the circuit has not yet caught
up with a frozen source that is itself still moving -- see "Validity check"
below for a diagnostic this tool computes for exactly that case, and
"Known limits" below for when to distrust a number.

This also explains, and is exactly why, the two validation identities the
card's own brief names are true **by construction**:
- **late `t_k`** (the transient has settled): the frozen source has ALSO
  settled to its own final value by then, so the frozen deck IS the
  steady-state deck -- equals a plain DC-op `.stb` of the unmodified deck.
- **early `t_k`**: equals "a DC-op `.stb` of the same circuit with its
  sources/ICs set to the frozen values" -- literally what is built and run.

### Validity check: `--check-tran "TSTEP TSTOP" --check-node NODE`

Runs the deck's own real, unmodified `.tran` once (same binary), and at
every `t_k` reports the worst relative deviation between each
`--check-node`'s REAL transient voltage (nearest sampled point to `t_k`)
and that same node's value in the frozen deck's own `op` solution -- a
direct, honest measure of how well the quasi-static assumption held at that
`t_k`. Real captured numbers, `tests/fixtures/stbtran_ldo_softstart.cir`
(`--check-node vout --check-node ea`):

| `t_k` | worst node rel. delta (real tran vs. frozen op) |
|---|---|
| 100 us (mid-ramp, `tss`=500 us) | **0.774** (77%) |
| 1000 us (settled, ~500 us after the ramp ends) | **9.5e-08** |

This fixture's error dynamic pole (`Rout*Cout` ~ 530 Hz, tau ~300 us) is
**deliberately slower** than its own 500 us soft-start ramp, so the loop
genuinely lags the frozen-time assumption during the ramp -- chosen
specifically so this diagnostic has a real, non-trivial case to catch (a
fixture where quasi-static held everywhere would only ever print a
rubber-stamp near-zero number and prove nothing about the check itself).
**Read `t_k=100us`'s PM/GM numbers on this fixture as "what the loop's own
margin would be if the ramp stopped moving at that instant," not as "what
the real circuit's instantaneous small-signal behaviour is at t=100us"** --
the two are different questions, and only the first one is what `.stb`
(frozen or otherwise) can ever answer. On a real design where the loop
bandwidth is comfortably above the frozen source's own rate of change (the
common, well-designed case), this delta is expected to be small at every
`t_k`, and the frozen-op PM/GM IS a trustworthy proxy for the instantaneous
margin -- always check this diagnostic before trusting a number from a fast
ramp / a loop with marginal bandwidth relative to it.

## Fixture: `tests/fixtures/stbtran_ldo_softstart.cir`

An LDO-style loop (single-pole error amp `G1`/`Rp1`/`Cp1`, nonlinear pass
device `Bpass` -- a saturating `tanh` VCCS so its incremental gm genuinely
depends on the DC operating point, not a linear-everywhere LDO which would
give the same PM/GM at every `t_k` regardless of whether the freeze
mechanism does anything at all -- output filter `Rout`/`Cout`/`Resr`,
feedback divider `Rf1`/`Rf2`) driven by a soft-start reference ramp `Vss`
(`PWL(0 0 {tss} {vref_final} 2m {vref_final})`, `tss`=500 us). `Vprobe` is
a native `.stb` probe (plain 0 V source) in the feedback wire, same
convention as `docs/examples/b_ldo_loop_native.cir`.

Measured, this build: `tran 50n 1.2m 0 uic` shows `v(vout)` rising with the
ramp, a small resonant overshoot to 3.0265 V at `t=514us` (PM is a real,
moderate 17-20 deg through the ramp, not a textbook-huge margin), settling
to 2.995537 V by `t~600us` and flat thereafter (2.995537 V at
600/800/1000/1100/1200 us). Frozen DC-op `.stb` (native, `dec 30 1 1e8`) at
several fixed `Vss` values:

| `Vss` (fixed) | PM | `fc` | GM |
|---|---|---|---|
| 0.09 (= `Vss(100us)`) | 17.288 deg | 81378.9 Hz | no crossing |
| 0.225 (= `Vss(250us)`) | 17.584 deg | 77327.2 Hz | no crossing |
| 0.36 (= `Vss(400us)`) | 18.370 deg | 69192.2 Hz | no crossing |
| 0.45 (final / settled) | 19.560 deg | 60744.0 Hz | no crossing |

A real, ~34% shift in crossover frequency and ~2.3 deg of PM across the
ramp -- confirms this fixture's nonlinearity is not negligible and the tool
is measuring something real, not a constant restated four times.

## Validation (`tests/test_stb_tran.py`)

Three layers, this build (`aaba50ba3f49-s46`), 2026-09-29:

1. **Pure-Python waveform evaluator + `.param` resolver** (no ngspice): 8
   hand-computed `PULSE` phases, 6 `PWL` points, 1 `SIN` phase-offset case,
   a `{param}`-referencing brace expression, and a red-proof (`EXP(...)` is
   refused, not silently mishandled) -- all independent of `tools/stb_tran.py`'s
   own `eval_pwl` (the `PWL` hand-check would be circular otherwise).
2. **The three freeze-evidence reproducers above**, run live and asserted
   (diode current within 1e-3 A of the true `4.33e-4 A` equilibrium in all
   three cases, not the forced/loaded value).
3. **End-to-end**, `tests/fixtures/stbtran_ldo_softstart.cir` at
   `t_k in {100,250,400,1000} us`: `tools/stb_tran.py`'s own computed
   `Vss(t_k)` matches an independent closed-form `vref_final*t/tss`; its
   PM/`fc` at each `t_k` match an **independently hand-built** frozen deck
   (a separate code path in the test file, not `stb_tran.freeze_deck()`,
   editing the `Vss` line directly) to **1e-6 relative**; `t=1000us` (settled)
   matches a plain DC-op `.stb` of the unmodified deck within **1 deg
   PM / 5% of `fc`** (well inside the card's own 0.1 dB/1 deg bound -- this
   fixture's late point already needed no rounding, both sides agree to
   4-5 significant figures); a deliberately mis-frozen reference (`Vss`
   off by 0.05 from the true `t=100us` value) is confirmed to NOT match, so
   the comparison is not vacuously loose; a missing `NGSPICE_AUG` binary is
   confirmed COULD-NOT-RUN, not a silent pass.

**23/23 pass** (0 skipped). Every ngspice-dependent check reports
COULD-NOT-RUN, not a silent pass, when `NGSPICE_AUG` is not set and the
default build is not present -- the 8 pure-Python checks (layer 1) still
run and can still FAIL on their own regardless.

## Known limits and caveats

- **Not a literal energy-storage-state freeze** -- see "What freeze means
  here" above. This is a documented, load-bearing trade-off, not an
  oversight: the engine's own `CKTop()` makes a literal freeze impossible
  at the script layer (verified, not assumed).
- **`--check-tran`'s validity diagnostic is the only defence against
  trusting a bad quasi-static point** -- without it, a `t_k` where the loop
  has not caught up with its own frozen source prints a PM/GM number with
  no visible warning that it may not represent the real instant well. Run
  it whenever the frozen source's own rate of change is not obviously much
  slower than the loop's bandwidth.
- **Only `PULSE`/`PWL`/`SIN`** are supported by `--freeze`'s waveform
  evaluator; `EXP`/`SFFM`/`AM`/a `PWL` using the `r=`/`td=` repeat options
  are refused outright (`CouldNotRun`), never silently mis-evaluated.
- **A `--freeze` source that is switching (a PWM ramp, a clock)**, not
  slowly varying, is a fundamentally different (and NOT quasi-static)
  situation this tool does not address -- freezing a switching source at an
  arbitrary phase gives the small-signal loop gain of a *fictitious static
  circuit held at that instantaneous switch position*, not a meaningful
  stability statement about the switching converter itself (that question
  is `pstb`'s job -- periodic AC about a PSS orbit -- which this project
  does not implement; `docs/stb.md`'s "Known limits" section).
- **Only one probe/sweep per run** -- no DM/CM (`mode=dm|cm`) or device-mode
  (`probe=<MOS instance>`) support in `stb_tran.py` itself, though nothing
  about the freeze mechanism is specific to loop mode; not exercised here
  (out of this card's four-fixture time budget).
- `--check-tran`'s "nearest sampled point to `t_k`" can be off by up to half
  a `TSTEP` if the transient's own adaptive step lands unevenly around
  `t_k` -- not a fixed-grid interpolation; fine for this doc's coarse
  diagnostic, not intended as a precise co-simulation.

## What was not done

- No `docs/examples/` walkthrough deck was added (unlike `docs/stb.md`'s
  four lettered examples) -- `tests/fixtures/stbtran_ldo_softstart.cir`
  plus this doc's own tables serve that role for this card; a dedicated
  `docs/examples/e_*` pair was not built, time budget.
- Only one fixture (LDO soft-start) was built and validated; the card's own
  wording also names "a load step" as a motivating case -- `--freeze` on a
  `PULSE`-shaped load-step source is exercised by the pure-Python `eval_pulse`
  unit checks, but no end-to-end load-step FIXTURE (deck + tool run +
  independent cross-check) was built, only the soft-start-ramp one.
- The `--check-tran` validity diagnostic was validated on exactly one
  fixture at exactly two `t_k` (a bad case and a good case) -- not swept
  across a range of ramp-speed-vs-loop-bandwidth ratios to characterise
  where the "small delta" threshold should actually be drawn; no numeric
  threshold is asserted by the tool itself (it reports the number, the
  reader judges it), which is a deliberate choice (see `docs/stb.md`'s own
  "never a silent number" convention) but means there is no tested
  pass/fail line for "quasi-static enough."
- `--plot-loopgain` was implemented and runs without error but its output
  file was not independently cross-checked against a separately-run native
  `.stb` loopgain trace at the same frozen point (the PM/`fc`/GM scalars
  already are, via the margins-file comparison above, and `stban.c`'s own
  margin/loopgain post-processing is the same code path already validated
  project-wide in `docs/stb.md` -- so this is a low-risk, not a zero-risk,
  gap).
- No independent re-verification by a second reader -- this doc, its
  fixture, and its numbers are all from a single session.
