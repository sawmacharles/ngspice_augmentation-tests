# Loop-compensation optimizer (`tools/compopt.py`)

Card #471 (plan T10.3), 2026-09-28. Tunes a compensator's named `.param` values so a
loop (probed with native `.stb`, docs/stb.md) hits a requested crossover frequency and
phase margin, starting from a detuned/bad starting point. Built on native `.stb` (each
evaluation is one `.stb` run); designed so the analysis card it writes can be swapped for
`.pstb` later (plan T5.2) without touching the optimizer, `Param`/`Target` or scoring code
above it — see `run_native_stb()`'s docstring in `tools/compopt.py`.

Never edits `tools/stb.py`, `tools/stbcore.py` or `tools/fra*.py` — they are imported,
unmodified, and reused: `stbcore.margins()` for crossing-finding (the SAME code this
project's native `.stb` and `tools/stb.py` already use and this project has already
validated — docs/stb.md, "How this was validated"), `stb.read_wrdata`/`stb.run` for
reading native `.stb`'s output and for the independent script-based check, and
`fra.sweep_loop`/`fra_curves.POINT_COMMON` for the independent switching-deck check on R3.

## Quick start

```sh
python3 tools/compopt.py tests/fixtures/compopt_r1b_detuned.cir \
    --lgprobe-inst XP \
    --param cc=1e-13:3e-11:log --param c1=1e-14:1e-11:log \
    --start cc=1e-12 --start c1=0.2e-12 \
    --target fc=7.6e6 --target pm=65 \
    --sweep "dec 30 10 1e9" --restarts 2 --max-iter 80 --verify stb
```

- `deck`: an **lgprobe-instance-form** deck (`XP a b lgprobe` + `.include .../lgprobe.lib`
  — the same form `refs/*.cir` and `tools/stb.py` use). `compopt.py` derives the
  native-`.stb`-compatible plain-V-source variant on the fly
  (`to_native_probe_deck()`) for every evaluation; the checked-in fixture stays in
  lgprobe form so the independent `tools/stb.py` check can run on the exact same deck
  text, unmodified.
- `--param name=lo:hi[:lin]`: a `.param` name to tune, log-scale by default (right for
  R/C values — a step means a fixed *ratio*, not a fixed absolute change); `:lin` for a
  linear-scale parameter.
- `--start name=value`: starting value (defaults to the bounds' geometric, or for `:lin`
  arithmetic, midpoint if omitted).
- `--target fc=<Hz> pm=<deg> [gm=<dB>] [single_falling=1]`: `fc`/`pm` are required. `gm`
  is an optional gain-margin floor (penalized, not gating, if a `-180°` crossing exists
  and reads below it — a loop with no such crossing is not penalized: "no crossing" is
  not the same as "fails the floor", docs/stb.md's own margin convention).
  `single_falling` is accepted for symmetry with the spec's wording but is **informational
  only**: the optimizer already refuses (heavily penalizes, never accepts as the reported
  optimum) any point with more than one falling 0 dB crossing **by default** — pass
  `--allow-multiple-falling` to relax that, not `single_falling=0`.
- `--verify stb`: after optimizing, also runs the independent `tools/stb.py` check
  (stock ngspice) and prints its margin. `verify_with_fra()` (R3's independent switching-
  deck check) is Python-API-only, since it needs a real switching deck path and a
  frequency list, not a single flag — see the R3 example below.

Exit status: 0 = target reached and (if `--verify` given) independently confirmed within
tolerance; 1 = ran, but the target is reported **infeasible** — never printed as success;
2 = could not run at all.

## How an evaluation works

1. `set_params()` regex-substitutes every named `.param name=value` token in the base
   deck text (whole-token, case-insensitive match; raises if a name is not found at all
   — a silently-ignored parameter would be exactly the kind of silent number this
   project's rules forbid).
2. `to_native_probe_deck()` replaces the lgprobe instance with a plain `Vprobe a b DC 0`
   and drops the `.include .../lgprobe.lib` line — native `.stb`'s `probe=` resolves to
   an actual V-source instance, not a subcircuit (docs/stb.md's `probe=` table).
3. `run_native_stb()` runs `~/Work/apps/ngspice_augmentation/<prefix>/bin/ngspice-aug -b`
   on a `.control` wrapper (`stb <sweep> probe=Vprobe`, `wrdata` the `T` vector), and
   raises `CouldNotRun` on any ngspice error or a missing/short output (e.g. `.stb`
   refused a bad `fstart`/`fstop`/`n` — docs/stb.md's "Frequency and point-count
   validation" table applies unchanged here).
4. `stbcore.margins(freqs, T)` (imported, unmodified) finds every 0 dB and -180° crossing,
   labels each falling/rising, and flags `low_f_positive_feedback` (`Re(T) < 0` at the
   lowest swept frequency) — exactly as native `.stb`'s own text output and `tools/stb.py`
   do it.
5. Scoring (`evaluate()`):
   - `low_f_positive_feedback` fires -> refused (status `"positive-feedback"`, cost
     floored high) — never the reported optimum.
   - No falling 0 dB crossing at all -> status `"no-falling-crossing"`, cost still
     finite and gradient-bearing (`|T|` in dB at the frequency nearest the target,
     squared — pushes the search toward creating a crossing there).
   - More than one falling crossing, `allow_multiple_falling` not set -> refused
     (status `"multi-falling"`).
   - Otherwise: the worst (lowest-margin) falling crossing is scored,
     `cost = w_fc*(log10(fc/fc_target))^2 + w_pm*((pm-pm_target)/60)^2`, plus a gain-margin
     penalty if `gm_min` given and the found gain margin (if any) is below it. `within_tol`
     is `True` only when `|fc-fc_target|/fc_target <= 5%` **and** `|pm-pm_target| <= 2°`
     **and** the gain-margin floor (if any) is met.

## The optimizer

Nelder-Mead (stdlib only — no numpy/scipy on this machine), in log10-parameter space for
every log-scale `Param`, with bound-clamping after every simplex move and random restarts
(`random.Random(seed)`, so a run is reproducible) perturbing the best point found so far.
Derivative-free by construction: each "function value" is one ngspice subprocess run, not
a closed-form expression a gradient could be taken of. All evaluations run **sequentially**
(one ngspice process at a time, well inside this project's CPU budget while another
worker benchmarks) — a 2-parameter, 2-restart optimization typically takes 100-250
evaluations and, on the fixtures below (small linear/averaged circuits), wall-clock well
under 2 seconds.

After the restarts, one more evaluation at the reported-best point is run cleanly
(`final/`) so the reported margins come from a fresh, isolated run, not whichever simplex
vertex happened to be evaluated last.

## Independent verification

- **`verify_with_stb_script()`**: re-parametrizes the SAME lgprobe-form deck text with the
  optimizer's final values and runs `tools/stb.py` (the script-based Tian testbench,
  unmodified) — by default against `NGSPICE` (stock system ngspice, a different binary
  from the fork build the optimizer itself used, and a different code path: two ordinary
  `ac` analyses plus Python's Tian formula, vs native `.stb`'s in-simulator two-solve).
  Its `margins_fn` parameter defaults to the real `stbcore.margins`; passing a corrupted
  one is this doc's own red-proof (below) and `tests/test_compopt.py`'s.
- **`verify_with_fra()`** (R3 only, spec's own requirement): applies the optimizer's final
  values to a **copy** of the real switching deck (`refs/r3_vm_buck.cir`, read-only, never
  edited in place) and runs `tools/fra.py`'s transient-injection Tian/Middlebrook harness
  (imported, unmodified) at a handful of frequencies around the found crossover, using the
  SAME converter settings `tools/fra_curves.py` documents for this exact deck family
  (`fra_curves.POINT_COMMON`: `settle_periods=3, measure_periods=5, points_per_period=256,
  startup="uic", tmax=5e-9, spice_options=".options trtol=1"`, and `fra_curves.AMPLITUDE`
  = 10 mV). This is **not** expected to match the averaged-model `.stb` result to `.stb`'s
  own precision — averaged != switching, especially near `fs/2` (spec's own wording) — so
  the check reports the *difference*, honestly, rather than asserting near-equality.

## Worked example: R1b, detuned two-stage amplifier

`tests/fixtures/compopt_r1b_detuned.cir` is `refs/r1b_bilateral.cir`'s topology (a
two-stage transconductance amplifier with a reverse path through `Cf`/`Cc`, closed-loop
gain of 2) with its own as-given `cc`/`c1` (the Miller/dominant-pole caps) as the
"detuned" starting point: native `.stb` on this deck gives **PM = 10.3 deg at
fc = 4.45e7 Hz** — barely stable. Tuning `cc` (bounds 0.1p-30p) and `c1` (bounds
0.01p-10p) to reach **fc = 7.6 MHz, PM = 65 deg**:

```sh
python3 tools/compopt.py tests/fixtures/compopt_r1b_detuned.cir \
    --lgprobe-inst XP \
    --param cc=1e-13:3e-11:log --param c1=1e-14:1e-11:log \
    --start cc=1e-12 --start c1=0.2e-12 \
    --target fc=7.6e6 --target pm=65 \
    --sweep "dec 30 10 1e9" --restarts 2 --max-iter 80 --verify stb
```

Real captured output (2026-09-28, this card):

```
path: ~/Work/data/ngspice_augmentation/compopt/_cli_sanity
evaluations: 129 (2 restart(s))
final values:
  cc = 1.00262e-11  (E96 ~= 1e-11, E24 ~= 1e-11)
  c1 = 2.02597e-12  (E96 ~= 2.05e-12, E24 ~= 2e-12)
fc = 7.60034e+06 Hz (target 7.6e+06, rel err 0.00%)
PM = 65.000 deg (target 65.000, abs err 0.000 deg)
TARGET REACHED (within 5% fc / 2 deg PM)
independent check (tools/stb.py, stock ngspice):
  PM = 65.000 deg @ 7.60034e+06 Hz
```

129 evaluations, well under a second of wall-clock. The independent `tools/stb.py` run
(stock ngspice, script-based Tian formula) agrees with native `.stb`'s own PM to 6
significant figures (65.000282 deg both ways) — the same near-exact native-vs-script
agreement docs/stb.md documents project-wide for loop-mode probes.

**Infeasible case** (same fixture, same bounds): target `fc=1e8, pm=89` (89 deg at a
crossover above R1b's second-pole region). The optimizer's best point after 243
evaluations across 3 restarts was `fc=1.71e7 Hz, PM=48.2 deg` — 83% off in frequency and
41 deg off in phase margin — reported as `TARGET NOT REACHED -- reporting INFEASIBLE`,
never as a success.

## Worked example: R3, averaged voltage-mode buck compensator

`tests/fixtures/compopt_r3_averaged_detuned.cir` replaces `refs/r3_vm_buck.cir`'s switch
network and PWM comparator with one linear VCVS implementing the standard state-space-
averaged PWM-switch model (`v(sw)_avg = Vin/Vramp * v(opout)`); everything downstream (the
LC output filter with ESR) and the type-III compensator are unchanged linear elements from
the real deck. See the fixture's own header comment for the full derivation.

**Validation of the averaging itself** (found during this card's design pass, not assumed):
at the real deck's own nominal compensator values, this averaged model's native `.stb`
gives **PM = 38.10 deg @ fc = 33660 Hz**. Card #452's own FRA reference curves
(`~/Work/data/ngspice_augmentation/fra/curves/SUMMARY.md`, generated independently, before
this card existed) measured the REAL switching deck's loop gain at **PM = 38.45 deg @
fc = 33463 Hz** — a 0.6%/0.35 deg difference. That is far closer than the real deck's own
header comment's hand-derived estimate (PM~=73.5 deg @ fc=50 kHz, explicitly labelled
there as "ESTIMATE ... not measured"), which the same #452 FRA sweep had already shown to
be inaccurate. This averaged model is a good linear stand-in for compensator tuning near
this crossover, even though — by construction — it drops every switching/sampling effect
(no half-`fs` peaking, no ripple) and should not be trusted near `fs/2` (250 kHz).

**Sweep note (card #476, 2026-09-29: fixed, no longer applies)**: this compensator's
`R1in`/`C1f` branch is a near-integrator (a pole close to the origin). Sweeping from
`fstart=10 Hz` used to give spurious `low_f_positive_feedback` false positives for
several otherwise-reasonable `(r2, c2)` values, because `Re(T)` genuinely straddles 0
that close to the near-integrator pole — not a real instability. Card #476 narrowed
`stbcore.margins()`'s `low_f_positive_feedback` trigger to require BOTH
`|arg T(fmin)| > 150 deg` AND `|T(fmin)| > 1` (positive feedback WITH loop gain — see
docs/stb.md's "Interpreting margins"), and added a neutral `low_f_integrator` label for
`arg T(fmin)` within -90 +/- 30 deg; a near-integrator branch like this one's `Re(T)`
straddling 0 near DC no longer trips the flag at all, at `fstart=10 Hz` or otherwise
— checked directly: `tests/fixtures/compopt_r3_averaged_detuned.cir` (this file's own
detuned starting point) and its tuned point (`tests/fixtures/stb6_r3_tuned.cir`) both run
clean, with `low_f_integrator` true, at `"dec 50 10 1e10"` in both `tools/stb.py` and
native `.stb` (`tests/test_stb.py`, `tests/test_stb_native2.py`). `fstart=500 Hz`
(matching `tools/fra_curves.py`'s own `fs/1000` lower bound for this deck family) is
still used below and by the fixture/`tests/test_compopt.py` — not because it is needed
to avoid the old false positive any more, but because it was already the chosen sweep
before this fix and changing it is outside this card's scope.

Detuned start `r2=100k, c2=0.3n` gives **PM = -5.4 deg @ fc = 55.8 kHz** — genuinely
unstable, not merely off-target. Tuning `r2` (bounds 1k-1M) and `c2` (bounds 10p-100n) to
reach **fc = 34 kHz, PM = 38 deg** (close to the real design's own nominal point — "recover
a working design"), real captured output from `tests/test_compopt.py` (2026-09-28):

```
final values: r2=27950.5, c2=2.4425e-09
fc=34000 Hz (target 34000, rel err 0.001%)
PM=37.998 deg (target 38.000, abs err 0.0016 deg)
151 evaluations, 0.61s
```

Independent confirmations:

- `tools/stb.py` on the same averaged linear circuit: PM = 37.998 deg @ 33999.6 Hz —
  matches the native result to 4 decimal places (same near-exact agreement as R1b).
- `tools/fra.py` on the **real switching deck** (`refs/r3_vm_buck.cir`, with `r2`/`c2`
  substituted to the optimized values), 5 frequency points around the found crossover,
  `fra_curves.POINT_COMMON` settings: **FRA measured PM = 38.38 deg @ fc = 33878 Hz** —
  a difference of 0.36% in frequency and 0.385 deg in phase margin from the averaged-model
  `.stb` result. This is well inside what the averaged-model-vs-switching agreement at
  card #452's own reference point suggested to expect, and is reported here as a genuine
  independent measurement, not assumed.

  Note (card #476, 2026-09-29: no longer misfires here, re-checked directly): the raw
  `T(f)` at these 2-5 points straddling the crossover has `Re(T) < 0` at the lowest of
  them (phase is already past -90 deg approaching the crossover, which is ordinary
  there). Before card #476, `stbcore.margins()`'s `low_f_positive_feedback` flag (which
  assumes the lowest swept point is near true DC) read `True` on this narrow FRA
  verification window from that sign alone — not a real low-frequency instability
  finding, just a known limitation of checking that flag over a window that does not
  start near DC. Card #476's narrower trigger (`|arg T(fmin)| > 150 deg` AND
  `|T(fmin)| > 1`, docs/stb.md's "Interpreting margins") does not fire here: re-run
  directly (`tools/compopt.verify_with_fra`'s own summary,
  `tests/test_compopt.py`'s `check_r3_feasible`), `low_f_positive_feedback` now reads
  `False` on this same window (`arg T(fmin)` is short of -150 deg approaching the
  crossover, not past it). `tests/test_compopt.py` still does not gate on this flag for
  the FRA window either way — it only compares the PM/fc numbers — so this was always
  informational, never a masked failure.

**Infeasible case** (same fixture, same bounds): target `fc=150000, pm=80`. A hand sweep
during this card's design pass found `{r2, c2}` alone (everything else fixed) cannot push
this compensator's crossover much past roughly 60 kHz regardless of gain (the fixed
`R3in`/`C3in` pole and the output-filter ESR zero dominate above that). The optimizer's
best point after 168 evaluations was `fc=29.1 kHz, PM=48.0 deg` — 81% off in frequency,
32 deg off in phase margin — reported `TARGET NOT REACHED`.

## Limits and what was not done

- **Not wired to `.pstb`.** `.pstb` does not exist in this fork yet (plan T5.2, not
  reached) — the module is structured so swapping it in later only touches
  `run_native_stb()` (see that function's docstring), but this was not tested against a
  real `.pstb` build because none exists.
- **`single_falling=0`/`1` in `--target` is parsed but not yet wired to override
  `--allow-multiple-falling`** — see "Quick start" above. Low risk (the default already
  matches the spec's "unless allowed" wording), but a CLI user asking for
  `single_falling=0` alone, without also passing `--allow-multiple-falling`, gets the
  strict (single-falling) behavior anyway. Not caught by any test here.
- **The E24/E96 rounding report** (`round_series()`) is a plain nearest-mantissa lookup
  against the standard IEC 60063 IEC tables, offered as a convenience in the CLI's printed
  output; it is not applied to the values used for the independent verification runs (those
  use the optimizer's raw, unrounded final values) and its own accuracy (does the rounded
  component actually still meet the target) was not separately checked.
- **Gain-margin floor (`gm=<dB>` target) and `--allow-multiple-falling`** are implemented
  (see "How an evaluation works") but not exercised by any fixture in
  `tests/test_compopt.py` — only `fc`/`pm` targets are tested end-to-end.
  `stbcore.margins()`'s own "no crossing = None, never a number" convention is preserved
  (a missing gain-margin crossing is not treated as "fails the floor").
  `no-falling-crossing` scoring (the "no valid falling crossing at all" branch of
  `evaluate()`) is exercised implicitly by the infeasible-target searches above (the
  optimizer visits many such points while searching) but has no dedicated, isolated test.
- **Only two fixtures, one target each (plus one infeasible target each)** — R1b and the
  R3 averaged model, 2 tunable parameters apiece. Nothing here tests 3+ parameter
  simplices, a `:lin` (non-log) parameter, or a deck with a device-mode probe.
- **The R3 FRA verification's frequency list is narrow** (5 points bracketing the found
  crossover, not a full `fs/1000`-to-`fs/2` sweep like `fra_curves.py`'s own reference
  curves) — chosen to keep runtime and CPU load modest (this card's hard limit: at most 2
  concurrent ngspice processes, another worker timing benchmarks). A wider sweep might
  show a larger averaged-vs-switching difference near `fs/2`; this was not checked here.
- **No independent re-verification by a second reader** — this doc, its fixtures,
  and every number in it are from this one card/session.
- **CPU discipline**: every ngspice run in `tools/compopt.py`, `tests/test_compopt.py` and
  this doc's examples runs one process at a time (`workers=1` where the underlying tool
  supports concurrency, e.g. `verify_with_fra`'s `fra.sweep_loop`); the full test suite
  (`tests/test_compopt.py`, 20 checks, both fixtures, both feasible and infeasible targets,
  both independent-verification paths, plus the red-proof) completes in about 12 seconds
  wall-clock on this machine.
