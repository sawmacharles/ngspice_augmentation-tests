# Converter measurement helpers (`convmeas`)

Card #497 (feature-survey #16, `research/feature-survey.md`), 2026-09-29.
`tools/convmeas.py`: settled ripple pk-pk, average efficiency, settling
time to a +-band, overshoot/undershoot on a step, switching-frequency
estimate, and duty cycle -- for any deck with a `.tran`. Needs only system
ngspice (`/usr/local/bin/ngspice`); nothing here uses the private fork.

Two layers, same split as `tools/stb.py`/`tools/stbcore.py`:
- **pure functions** (`ripple_pp`, `average_value`, `average_power`,
  `efficiency`, `switching_frequency`, `duty_cycle_threshold`,
  `settling_time`, `overshoot_undershoot`) -- stdlib only, operate on plain
  `(t, y)` sample lists already in memory. Import these directly if you
  already have a trace loaded some other way (a `wrdata` file you read
  yourself, a `.raw` file read with another tool, ...).
- **a thin CLI driver** that runs one `.tran` (the documented `.control`
  snippet below, or read `tools/convmeas.py`'s own `run_tran()`), reads the
  resulting trace, and prints/writes every requested measurement.

```sh
python3 tools/convmeas.py DECK --vin Vin --vout-node vout --rload-ohms 3.3 \
    --pwm-node pwm --period 2e-6 --periods 5 --tstop 1.5e-3 \
    [--target 3.3 --band-pct 1] [--load-step-time 1e-3]
```

Prints a report; writes `wrapper.cir`, `ngspice.log`, `trace.txt` (raw
`wrdata`) and `summary.json` under
`~/Work/data/ngspice_augmentation/convmeas/<deck>-<timestamp>/` (never the
repository) unless `--out` is given. Exit status: 0 = ran and computed
every requested measurement; 2 = could not run.

## The documented `.control` snippet set

If you would rather instrument a deck yourself (or already have your own
`.tran`/`wrdata` machinery, as most of this project's `refs/*.cir` decks
do) and just want the pure functions, this is the snippet `tools/
convmeas.py`'s own `run_tran()` builds, reusable directly in any deck:

```spice
.control
set wr_singlescale
set wr_vecnames
set numdgt=15
tran {tstep} {tstop} 0 uic
wrdata trace.txt v(vin_node) i(Vin) v(vout_node) v(pwm_node)
.endc
.end
```

`set wr_singlescale` + `set wr_vecnames` is what makes `wrdata` emit ONE
shared time column (with a header line) instead of a separate time column
per signal -- `tools/convmeas.py`'s own `read_wrdata_real()` (and this doc's
own hand-parsing, if you write your own reader) expects exactly this
layout. Load current is **not** requested as `i(Rload)` -- ngspice does not
track a plain resistor's branch current as an ordinary circuit vector any
more than it tracks an `@device[i]` OP parameter without an explicit `save`
(`tools/loss.py`'s own "WHY `@device[i]`/`@device[p]` LOOKED FROZEN" finding,
`docs/loss.md`) -- `iout = v(vout_node) / rload_ohms` (Ohm's law, using the
already-known resistor value) sidesteps that trap entirely and needs no
`save` line at all. Input current, `i(Vin)`, IS a native MNA branch current
(every V-source's is) and needs no `save`.

**Input-current sign convention** (matches `tools/loss.py`'s
`p_in = -avg(v*i)`, same reasoning, same numbers on the same circuit): a
SPICE voltage source's own branch current `i(Vxxx)` is defined flowing INTO
its `+` terminal, so the power it actually DELIVERS to the rest of the
circuit is `-v(t)*i(t)`, not `+v(t)*i(t)`. `efficiency()` negates it
internally -- pass `i(Vin)` straight through, do not pre-negate it yourself.

## Functions

| Function | Signature | What it computes |
|---|---|---|
| `ripple_pp(t, y, t0, t1)` | -> float | `max(y) - min(y)` over `[t0, t1]` -- settled ripple peak-to-peak |
| `average_value(t, y, t0, t1)` | -> float | trapezoidal time-average of `y` over `[t0, t1]` -- exact for the AVG-over-one-period duty-cycle trick (see below) |
| `average_power(t, v, i, t0, t1)` | -> float | trapezoidal average of `v(t)*i(t)` -- `v`/`i` must share `t`'s grid |
| `efficiency(t, vin, iin, vout, iout, t0, t1)` | -> float | `Pout/Pin`, `iin` in the native `i(Vxxx)` sign convention (see above) |
| `switching_frequency(t, y, t0, t1, threshold=0.5)` | -> float\|None | `1/period`, from the average spacing of rising `threshold` crossings (first to last, divided by count-1) -- `None` (never a fabricated number) if fewer than 2 crossings fall in the window |
| `duty_cycle_threshold(t, y, t0, t1, threshold=0.5)` | -> float | fraction of time `y(t) >= threshold` in `[t0, t1]`, via linear interpolation for each crossing instant -- general, works for any waveform shape |
| `settling_time(t, y, target, band, t_start=0.0)` | -> float\|None | the LAST time `y` was outside `[target-band, target+band]`, for `t >= t_start` (the standard, robust definition -- not "first time it entered the band," which a signal can cross on its way to overshooting); `t_start` if always within band; `None` (never a fabricated time) if still outside the band at the very last sample |
| `overshoot_undershoot(t, y, t_event, y_before, y_after)` | -> dict | `overshoot_pct`/`undershoot_pct` (both >= 0, as a % of `abs(y_after - y_before)`), `peak_max`/`peak_min` |

### Duty cycle: two methods, on purpose

`average_value()` over exactly one period IS the duty cycle for a signal
that saturates close to a clean 0/1 square wave -- this is
`refs/r3_vm_buck.cir`'s own already-validated convention (its header:
*"average `pwm` (a smooth 0/1 signal) over a plain FROM/TO window exactly
one period wide ... a periodic signal's average over any exactly-one-period
window equals the duty cycle regardless of phase"*), and it needs no
threshold choice at all. `duty_cycle_threshold()` is the more general
alternative -- any waveform shape, an explicit threshold, does not need the
window to be exactly one period wide. Both are provided; `tools/
convmeas.py`'s own CLI computes and reports both (`duty_cycle_avg` /
`duty_cycle_threshold`) so a caller can sanity-check one against the other.

### `overshoot_undershoot`'s undershoot is searched from the peak onward, not from `t_event`

A real transient does not jump instantly -- for a short while after
`t_event` it is still close to `y_before`. Searching the WHOLE post-event
range for "the extremum opposite the step direction" would then just find
that still-near-`y_before` starting point, not the genuine undershoot dip
that follows the first overshoot. `overshoot_undershoot()` therefore finds
the overshoot peak first (the global extremum in the step's own direction,
over the whole post-event range), then searches only from THAT peak's own
time onward for the undershoot extremum -- correct for the common
decaying-oscillation case (each successive peak/trough is smaller in
magnitude than the last, so the first trough after the overshoot peak IS
the global minimum of the remaining trace). Found and fixed during this
card's own validation (`tests/test_convmeas.py`'s closed-form check
originally reported ~100% undershoot on a clean, textbook step response,
because the naive "global min over the whole post-event window" picked up
the pre-rise starting point at `t_event+0.2us`, not the real trough at
`t_event+208us` -- see that file's `check_closed_form_overshoot()`).

## Validation (`tests/test_convmeas.py`)

Three layers, 25/25 pass (system `ngspice`, this build's default):

1. **Pure-Python synthetic checks**, hand-built `(t, y)` arrays with exactly
   known answers, no ngspice: a 1 kHz/0.3 V sine's `ripple_pp` = 0.6 V; a
   sine's `average_value` over ~2 periods = its DC offset; a hand-built 30%
   duty square wave's `average_value` = 0.30; constant `V*I` gives an exact
   `average_power`/`efficiency`; a hand-built ideal 25 kHz/40%-duty square
   wave's `switching_frequency`/`duty_cycle_threshold`/`average_value`
   (duty, cross-checked two ways) all match; `settling_time` on a clean,
   hand-picked step (last-outside-band sample), an already-settled case, and
   a red-proof (never settles -> `None`, not a fabricated time); a flat
   (never-crossing) signal's `switching_frequency` -> `None`, not a
   fabricated number.
2. **A closed-form underdamped 2nd-order step response** (textbook formula,
   `y(t) = 1 - exp(-zeta*wn*t)/sqrt(1-zeta^2) * sin(wd*t + arccos(zeta))`,
   `zeta=0.3`), built as a dense Python array (no ngspice): `overshoot_pct`
   matches `Mp = exp(-zeta*pi/sqrt(1-zeta^2)) = 37.2326%` and
   `undershoot_pct` matches the second peak's `Mp^2 = 13.8627%` (successive-
   peak deviations of this canonical system decay by the factor `Mp` each
   half-cycle, so the `k`-th peak's deviation fraction is `Mp^k`); a
   red-proof (a different `zeta`'s `Mp` does NOT match this `zeta`'s
   measured overshoot).
3. **ngspice-driven fixtures**:
   - `tests/fixtures/convmeas_sine.cir` (`SIN(2.5 0.6 10k)`): real simulated
     `ripple_pp` = 1.2 V, matching the hand-computed `2*VA` to 5e-3 relative.
   - `tests/fixtures/convmeas_rlc_step.cir` -- the SAME closed-form response
     as layer 2, but now an ACTUAL series R-L-C circuit built for exactly
     `zeta=0.3` (`R=18.9736659610... Ohms` computed in Python to 10 sig figs
     for `L=1mH, C=1uF`), and its real, ngspice-simulated `v(vout)` trace:
     overshoot 37.2216% (closed-form 37.2326%, matches to 1e-2 relative),
     undershoot 13.8696% (closed-form 13.8627%, matches to 2e-2 relative),
     `settling_time` (+-2% band) matches a dense closed-form scan to within
     5 us. The deck's own `.meas` statements independently confirm the
     circuit itself is correctly built: `peak1=1.372326` (= `1+Mp` exactly)
     at `t=114.17us` (closed-form predicts `td+pi/wd=114.14us`),
     `peak2=0.861373` (= `1-Mp^2` exactly) -- see the fixture's own header.
   - `refs/r3_vm_buck.cir` (a real, independently-documented switching
     converter, `refs/README.md`'s own MEASURED table -- not built for this
     card): `ripple_pp` 4.8358 mV vs. the deck's documented 4.836 mV (within
     15%, a deliberately loose bound -- this tool re-runs its OWN `.tran` at
     its OWN `tstep`, not the deck's original 10 ns/`.meas` windows, so some
     numeric difference from a different-but-equally-valid sampling is
     expected, not a bug); `duty_cycle_avg` 0.2748 vs. documented 0.2777
     (within 10%); `switching_frequency` 500110.9 Hz, within 0.5% of the
     IDEAL 500 kHz clock (see "Why `switching_frequency` doesn't match the
     deck's own 500.30 kHz number" below); `efficiency` 21.85% (Pin=15.1 W,
     Pout=3.299 W) -- sane and non-`NaN`, no independently-documented
     reference number exists for this exact quantity on this deck (`tools/
     loss.py`'s own ~29% figure on R3/R4 is a DIFFERENT, device-level
     energy-balance measurement, not directly comparable number-for-number
     -- see below).
   - A CLI red-proof (missing ngspice binary -> `CouldNotRun`, exit 2) and
     one real subprocess invocation of the CLI end-to-end (exit 0,
     `summary.json` written) -- the library checks above cover the maths,
     this covers the CLI plumbing actually working.

### Why `switching_frequency` doesn't match the deck's own 500.30 kHz number

`refs/r3_vm_buck.cir`'s own header reports `fsw measured 500.30kHz` from a
SINGLE-period `.meas tran t_period TRIG v(pwm) VAL=0.5 RISE=250 TARG v(pwm)
VAL=0.5 RISE=251` -- one specific period's timing, which is measurably
perturbed by that period's own loop-dynamics-driven duty jitter.
`tools/convmeas.py`'s `switching_frequency()` instead averages over SEVERAL
periods end-to-end (first crossing to last, divided by the crossing count
minus 1), which converges to the driving PWM ramp's true long-run clock
rate (`Vramp`'s own `PULSE(0 {vramp} 0 {0.97*tclk} {0.03*tclk} 1n {tclk})`,
exactly 500 kHz by construction) -- a materially different, and for a
"what frequency is this converter really switching at" question, more
precise method, not a bug or a mismatch to be explained away. `duty_cycle`
has the analogous, smaller effect (this project's own single-window
`.meas` numbers, 0.2762/0.2777 at two different times, already differ from
each other by more than `convmeas.py`'s multi-period average differs from
either).

### Why `efficiency` has no independently-documented cross-check number

`tools/loss.py` (card #470) already computes a full per-device
energy-balance efficiency for R3/R4 (~29%, `refs/README.md`'s own
"WHY R3/R3d's efficiency..." discussion, shoot-through between the smooth
`Bhs`/`Bls` switches). `tools/convmeas.py`'s `efficiency()` is a much
simpler `Pout/Pin` from just the input source and the load resistor's own
Ohm's-law current -- it should be in the same ballpark (confirmed: 21.85%
here vs. `loss.py`'s ~29%, same order of magnitude, same "this deck's
shoot-through is real and dominant" conclusion) but is not expected to
match `loss.py`'s number exactly, since `loss.py` integrates over a
different, `.pss`-aware window and accounts for every device individually
rather than just the two terminals. Read `convmeas.py`'s `efficiency()` as
a fast, single-number sanity check, not a replacement for `tools/loss.py`'s
full breakdown when the loss BREAKDOWN (not just the total) is what is
needed.

## Known limits and caveats

- **`switching_frequency`/`duty_cycle_threshold` assume a clean, mostly-
  saturating signal** (a comparator or PWM node) -- neither has been
  validated against something like a switching converter's own inductor
  current ripple (a triangular, not square, waveform) where "duty cycle"
  is not naturally a 0/1 threshold question.
- **`overshoot_undershoot` assumes the classic decaying-oscillation shape**
  (each successive extremum smaller than the last) -- see "undershoot is
  searched from the peak onward" above; a pathological waveform whose
  SECOND excursion is larger than its first (not a standard linear-system
  step response) would not be handled correctly by the "first trough after
  the peak is the global minimum of the rest" assumption.
- **`settling_time`'s band is an absolute value**, not a percentage --
  the CLI's `--band-pct` converts `band_pct/100 * abs(target)` to an
  absolute band before calling the pure function; a caller using the pure
  function directly must do this conversion themselves if a relative band
  is wanted.
- **No device-level loss breakdown** -- `tools/loss.py` already does that;
  `convmeas.py`'s `efficiency()` is deliberately a simpler two-terminal
  `Pout/Pin`, not a replacement.
- **`run_tran()`'s own `.tran` uses `uic`** (matching every fixture and
  `refs/*.cir` deck in this project) -- a deck relying on a real DC
  operating-point solve before its transient (no `IC=`/`.ic` anywhere) would
  start from all-zero initial conditions instead; not exercised here since
  every fixture and `refs/` deck this project has already uses `uic`.

## What was not done

- No dedicated fixture for switching-frequency/duty-cycle measured
  purely synthetically through ngspice (a `PULSE`-driven square wave deck)
  -- the hand-built-array checks (layer 1) cover the maths exactly; the
  real-converter fixture (`refs/r3_vm_buck.cir`, layer 3) covers "does this
  work on a real, non-ideal PWM node," but nothing in between (an ngspice-
  simulated, but still clean/synthetic, square wave) was built, time budget.
- `--load-step-time`/`overshoot_undershoot`'s CLI path was smoke-tested
  directly (a manual run against `convmeas_rlc_step.cir`, printed in this
  doc's own development notes) but not added as a CHECKED-IN test case in
  `tests/test_convmeas.py` -- the equivalent library-level call IS checked
  in (`check_rlc_step_fixture`), covering the same maths through the same
  code path `compute()` calls, just not through `main()`'s own argument
  parsing for this specific flag combination.
- `duty_cycle_avg`/`duty_cycle_threshold`'s cross-check against each other
  was only done on the hand-built synthetic square wave (layer 1), not on
  `refs/r3_vm_buck.cir`'s real, smoothly-transitioning `pwm` node -- the
  CLI computes and reports both there too, but the test only asserts
  `duty_cycle_avg` against the deck's own documented number.
- No independent re-verification by a second reader -- this doc, its
  fixtures, and its numbers are all from a single session.
