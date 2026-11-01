# Reference circuit set (R2-R6) — ngspice_augmentation, card #450, plan T2.1

All decks run with `ngspice -b <deck>.cir` from this `refs/` directory (they
`.include ../stb/lgprobe.lib` by relative path, and `.control` writes a raw
file plus a `wrdata` text trace under `~/Work/data/ngspice_augmentation/refs/<deck>/`
— never into the repo). Header comments in each `.cir` carry the full design
derivation, the measured numbers, and every trap hit while building it;
this table is a summary, not a replacement for reading them.

R1 (op-amp follower / LDO / fully differential, card #451) is a separate
worker's deck (`refs/r1*.cir`) and is not covered here. `stb/`, `tools/fra.py`
and `tests/` belong to other workers in progress and were not touched.

| deck | purpose | key values | run | measured steady state | runtime | plan role |
|---|---|---|---|---|---|---|
| `r2_switched_rc.cir` | switched-RC LPTV circuit with an **exact** (not averaged) sideband-0 transfer function, derived by piecewise-LTI state propagation and given in closed form in `r2_analytic.py` | Rin=Rsw=10k, C=1nF, fclk=100kHz, duty=0.5, Vin=DC 1V (AC 1 marked, unused by this `.tran`) | `ngspice -b r2_switched_rc.cir` | fclk 100.000kHz (live); duty 0.5001; v(vx) avg (last period) 6.699138e-01V, ripple 1.600631e-01V; settled to ~1e-5 relative by 100us; near-DC cross-check against `r2_analytic.py` agrees to ~7e-5 (residual from finite switch Ron/Roff and 25ns tran step, not an algebra error — the script's own `--selftest` closed-form-vs-brute-force check agrees to ~3e-11) | 0.13s | T5.1 `pac_sb0` exact reference (R2 vs analytic) |
| `r3_vm_buck.cir` | voltage-mode synchronous buck, fixed 500kHz, type-III compensator (ideal-ish op-amp + R/C network), smooth (tanh) switches and PWM comparator, `XP` loop-gain probe in the feedback divider path | 12V->3.3V, L=10uH C=100uF Esr=10mOhm Rload=3.3, fs=500kHz, R1=10k R2=27.51k R3=101.7 C1=36.35pF C2=2.299nF C3=6.261nF | `ngspice -b r3_vm_buck.cir` | fsw 500.30kHz; duty 0.2762(500us)/0.2777(1500us) vs ideal 0.275; Vout 3.300360/3.299680/3.299578 V at 500/1000/1500us (settles to <1% by ~220us, monotonic no ringing); ripple 4.836mV | 2.1s | T2.3 FRA reference, T5.1/T5.2 `pac_sb0`/`pstb_loop` target; averaged-model PM estimate 73.5 deg at fc=50kHz (estimate, not measured) |
| `r4_pcm_buck.cir` | peak-current-mode buck, fixed 500kHz clock, smooth RS latch + smooth current comparator, slope compensation as `.param vslope_amp` (0 by default), type-II voltage loop, duty set by `.param vin_v` | Same power stage as R3; Rsense=0.05 Ohm amplified `rigain`=8 (Ri_eff=0.4); R1=10k R2=50.51k C1=25.21pF C2=6.534nF | `ngspice -b r4_pcm_buck.cir` (default vin_v=11, D~0.30, vslope_amp=0) | **D~0.30** (vin_v=11): fsw 500.04kHz, duty 0.3055(final), Vout settles to 3.300018V by ~3ms, ripple 4.651mV. **D~0.45** (vin_v=7.33, verification run, not checked in separately): duty 0.4583, Vout 3.300017V, ripple 3.619mV, stable. **D~0.55, no slope** (vin_v=6.0): **unstable** — period-2 subharmonic (t_period doubles to 4.000027us / 250kHz), ripple rises to 7.409mV — the textbook peak-current-mode instability, and exactly verify.md's T4.3 red-proof case. **D~0.55, vslope_amp=0.15**: stable again, duty 0.5600, Vout 3.300020V, ripple 2.935mV (lower than the D<0.5 cases) | 7.7s (default 6ms run); ~7.2-7.6s each for the D=0.45/0.55 verification variants | T5.3 `pac_sampling` (sampling/subharmonic effects across duty, with the D>0.5-no-slope-comp instability as the required red-proof case) |
| `r5_hysteretic_buck.cir` | free-running (self-oscillating) hysteretic buck, no clock, no compensator — direct voltage hysteresis via a self-referencing (regenerative) smooth Schmitt trigger | 12V->3.3V target, hys=20mV band, L=10uH C=100uF Esr=10mOhm Rload=3.3 (same power stage as R3/R4/R6) | `ngspice -b r5_hysteretic_buck.cir` | Large cold-start overshoot (Vout peaks ~7.6V, no soft-start/current-limit), settles to steady free-running oscillation by ~400us. **fsw ~427.2-427.6kHz** (an output of the simulation, not a parameter — see the "what was not done" note below on this number). Vout 3.3004V(1000us)/3.3005V(3000us final); ripple 5.637mV (0.17%); duty(final) 0.2789 vs ideal 0.275 | 2.3s | T6.2 `pss_freerun` (free-running/variable-frequency converter, needs the period-as-unknown shooting the plan describes) |
| `r6_vm_buck_ideal.cir` | R3 rebuilt with ideal `sw` switches (native ngspice, not XSPICE) and the XSPICE `d_pwm` code model (+ `dac_bridge`) instead of smooth devices; same power stage and compensator values as R3 for direct comparison | Identical to R3 (12V->3.3V, fs=500kHz, same L/C/Esr/Rload, same R1..C3, same op-amp) | `ngspice -b r6_vm_buck_ideal.cir` | fsw 497.85kHz; duty 0.2783(500us)/0.2712(1500us) vs ideal 0.275, comparable to R3's 0.2762/0.2777; Vout 3.294988/3.299978/3.300021 V at 500/1000/1500us; ripple 4.811mV, essentially matching R3's 4.836mV | 3.1s | T2.1 R6 target; T6.1 `pac_events` red-proof — "an ideal-switch deck (R6) must fail [PAC/PSTB] before Phase 6" per verify.md, because the abrupt `sw` events have no saltation correction yet. This deck's *large-signal* (`.tran`) behaviour matches R3 closely, as expected; only the future small-signal/PAC comparison is expected to diverge |

## Traps found and fixed while building these (also in each deck's header)

- **ngspice `meas ... param='expr referencing another meas'`** fails in this
  build in both control-mode `meas` and netlist `.meas` ("no such function
  as 'param=...'"), and separately **`AVG` combined with `TRIG`/`TARG`**
  fails too ("bad syntax, equal sign missing"). Worked around with `let` for
  derived quantities and plain `FROM=`/`TO=` windows for `AVG`.
- **`.meas ... TRIG ... RISE=n TARG ... FALL=n`** desynchronises during a
  noisy startup transient (extra comparator chatter can offset the RISE and
  FALL counters relative to each other), silently returning an inverted
  (negative) window and garbage — exactly the V-G7 trap the task called out.
  Fixed by measuring duty as `AVG` of the actual switch-node voltage
  (`v(sw)`, which saturates close to a true 0/Vin square wave) over a
  plain one-period-wide window, not by counting edges of a soft internal
  gate signal.
- **B-source current sign**: `Bxxx n1 n2 I=expr` delivers `expr` as current
  flowing n1->n2 (resistor-like), so a capacitor at n1 referenced to n2 sees
  `-expr`, not `+expr`. Writing a capacitor-charging ODE with the "obvious"
  sign made R4's smooth RS latch diverge to `-1e28` within a few
  microseconds; found with an isolated 6-line reproducer and fixed by
  negating the current expression.
- **Op-amp compensator local instability**: a "realistic" 100Hz-dominant-pole
  op-amp (GBW=10MHz) wrapped in R3/R4's type-II/III feedback network was
  *not* locally stable in isolation (no power stage at all) — `v(opout)`
  diverged to `1e98`. The feedback network's own noise gain intersects the
  amplifier's rolloff with too little phase margin at the amplifier's own
  (very high, since A0=1e5) local crossover. Fixed by pushing the op-amp's
  own pole out to ~16MHz (`co_oa=100p`), which is not a realistic op-amp
  GBW but is a reasonable "ideal-ish" approximation given the task's wording.
- **Ideal-switch (R6) convergence**: the same power stage that is fine with
  smooth switches blew up (`v(vout)` past `1e270`) around t=369us with ideal
  `sw` switches driven by independently-bridged complementary digital gate
  signals. Fixed with a standard RC snubber (10 Ohm + 1nF) at the switch
  node — a legitimate, common technique, not a change to the circuit's
  intent, but required here and noted for honesty.
- **R4 current-sense loading**: an initial design used a single 0.4 Ohm
  sense resistor directly in the power path (chosen for a convenient
  comparator voltage scale) and it dropped enough voltage at ~1A to
  visibly distort Vout and the effective duty cycle. Fixed by using a small
  0.05 Ohm physical sense resistor plus a separate x8 gain block, giving the
  same comparator-side signal scale without loading the power path.
- **R4 smooth-latch floor/ceiling**: the smooth RS latch's `q` did not
  swing cleanly 0..1 with the initial gain constants — it plateaued
  around 0.13/0.89 because the reset comparator is only briefly active
  right at the current crossing, then its own pull-down effectively
  vanishes (found via `avg(v(q))` giving duty=0.41 against an ideal
  Vout/Vin=0.30, even though Vout itself regulated correctly). The
  switches still worked correctly (they saturate well before q reaches
  either end), but `avg(v(q))` is not a usable duty proxy here; `avg(v(sw))`
  is used instead. Tightened anyway (kc 30->100, greset 1e-4->1) to narrow
  the floor/ceiling gap.
- **R5's two-comparator latch had a real dead zone**, distinct from R4's
  floor/ceiling: once the startup ring-down decayed enough, `q` drifted
  into the region where both the set and reset comparators are near-zero
  and *froze* at a non-binary value (q=0.484), settling the whole converter
  into a static linear-regulator fixed point (constant Vout, zero
  switching) instead of continuing to oscillate. Confirmed by running out
  to 800us with no change. Fixed by rebuilding the latch as a genuine
  self-referencing (regenerative) Schmitt trigger, whose threshold depends
  on its own output, removing the fixed point entirely. A pure-algebraic
  (zero-capacitance) version of that Schmitt trigger hit "timestep too
  small" exactly at the switching instant (a real jump discontinuity with
  nothing to regularise it in time); fixed by driving it through a small
  capacitor (finite, fast relaxation dynamics) instead.
- **R5 switching-frequency misread**: a first debug look at `v(q)` sampled
  every ~6.6us appeared to show a ~47us (~21kHz) oscillation envelope. That
  was aliasing — the real period, confirmed by exact `TRIG`/`TARG`
  crossing-time measurements (not sampled plotting), is ~2.34us (~427kHz).

## What was not done

- No FRA (transient-injection frequency response) was run against any of
  these decks — that is T2.3's job (`fra_ref`), a separate card/worker
  building `tools/fra.py` in parallel; these decks only had their `.tran`
  run to steady state per T2.1's exit condition.
- No `stb`/loop-gain probe run was exercised through `XP` on R3/R4 (only
  placed in the netlist) — `tools/stb.py`/`stbcore.py` belong to another
  worker; not touched here.
- R3's and R4's compensator designs are averaged-model *estimates*, computed
  with a standalone Python scratch calculation (not delivered as a repo
  file — only R2's analytic recipe was asked to ship as a script), and are
  explicitly labelled as such in the deck headers; no PAC/PSTB or `.ac`
  loop-gain measurement exists yet to check them against.
- R5's hysteresis band (20mV) and R4's slope-compensation amplitude (0.15V
  for the D>0.5 stability check) were picked by a rough back-of-envelope
  estimate and one round of "does it work", not finely tuned or optimised.
- R4's D=0.45 and D=0.55 (both with and without slope compensation) cases
  were run as one-off verification copies (sed-substituted `.param`
  overrides) to characterise the deck's behaviour across duty cycles, per
  the task; they are documented in this README and in `r4_pcm_buck.cir`'s
  header but are not separate checked-in `.cir` files — only the
  `vin_v=11, vslope_amp=0` (D~0.30) case is.
- No dead-time was modelled between high-side/low-side switches in any
  deck (R3/R4/R5's smooth switches naturally overlap near their crossover,
  which is why they did not need R6's snubber).
- The private #448 target decks are not part of this set (per plan T2.1,
  item 7) — #448 was still undecided as of this card.
- No independent re-verification (V-G4) has happened yet: all six numbers
  above were produced and read by this same card.

## Review notes (2026-09-27)

- All five decks re-run on ngspice-45.2: exit 0, numbers above reproduced (R3 Vout 3.299578 V, R4 3.300018 V, R5 3.300484 V at 427 kHz, R6 3.300021 V; periods 2.0 us, R5 2.34 us).
- Probe orientation corrected in R3, R4, R6 to `XP fbtop vout` (a = downstream, b = upstream, per stb/lgprobe.lib). Tian's T is unaffected; Middlebrook comparison output needs this order.
- R6: XSPICE `d_pwm` prints "Error: duty cycle 0 is not positive" / "1 is 1 or larger" during start-up and clamps duty to 0.01..0.99. Steady-state duty (0.27) is inside that range, so results are unaffected; tools that scan logs for "error" must allow these two lines for R6.
- `duty_meas` in R4/R5 `.meas` output is avg(v(sw)) in volts; the `let` lines after it divide by Vin to give the duty cycle.

## R3d/R4d — realistic (dead-time) counterparts of R3/R4, card #479, 2026-09-29

`tools/loss.py` (card #470, `docs/loss.md`) measured R3 and R4 at ~29%
efficiency: their smooth Bhs/Bls switches are both driven from the same
`pwm`/`1-pwm` signal with no separation in time, so every transition is
genuine shoot-through (hundreds of amps, confirmed against ngspice's own
solved `i(Vin)`). R3 and R4 themselves are UNCHANGED (still the math
references Phase 4-5 validations use, per this file's own table above).
`r3d_vm_buck_deadtime.cir` and `r4d_pcm_buck_deadtime.cir` are new, separate
decks alongside them: same power stage, same compensator/current-loop
component values, same probe placement and orientation, plus: non-
overlapping gate drive with a dead-time (a smooth RC-charge / diode-fast-
discharge steering network per side — break-before-make, not a symmetric
delay), realistic switch Ron (14.3 mOhm), a smooth low-side body/freewheel
diode for the dead-time interval, and gate-drive loss stated as excluded
(not modelled), all built from smooth devices only (no ideal `sw`, no
XSPICE, per #447 A). Full design derivation, every trap hit, and the
component-value reasoning are in each deck's own header comment; this table
is a summary.

| deck | purpose | key values (vs R3/R4) | run | measured steady state | efficiency (`tools/loss.py`) | runtime | test |
|---|---|---|---|---|---|---|---|
| `r3d_vm_buck_deadtime.cir` | R3 (voltage-mode buck) rebuilt with dead-time | gon=70 (Ron=14.3mOhm, R3 was 10mOhm); rdt=430 rdisch=60 cdt=100p kswd=150 (dead-time network); kc=100 (R3 was 30, raised so the PWM comparator's own edge doesn't dominate the dead-time delay); Bbody vfbody=0.6 gbody=50 (body/freewheel diode, as a B-source — see below) | `ngspice -b r3d_vm_buck_deadtime.cir` | ngspice-45.2: Vout 3.300007V, fsw 499.95kHz, duty 0.2870, ripple 4.931mV. Fork (aaba50ba3f49-s14): Vout 3.299975V, fsw 500.00kHz — agrees to ~3e-5V. Peak simultaneous hs/ls conduction current (independently re-derived from the trace, not just the gate-crossing times): 1.26e-05 A, vs a 1%-of-load bound of 1.00e-02 A | P_in=3.3558W P_out=3.3000W, **98.34%** (vs R3's 28.65%), energy-balance residual -0.0023% of P_in | 8.8s (ngspice-45.2, 1.5ms/750 periods) | `tests/test_refs_d.py`: 4/4 checks pass |
| `r4d_pcm_buck_deadtime.cir` | R4 (peak-current-mode buck) rebuilt with dead-time | same gon/dead-time-network values as R3d; `vslope_amp`/`vin_v` unchanged from R4 (D~0.30 default, slope-comp param still there for the D>0.5 case); kc LEFT at R4's own 100 (unlike R3d — see header, a first trial at 260 was unstable before a loading bug was found) | `ngspice -b r4d_pcm_buck_deadtime.cir` | ngspice-45.2: Vout 3.300026V, fsw 500.00kHz, duty 0.3048, ripple 4.690mV, il 0.766-1.236A (matches R4's own 0.768-1.233A). Fork: Vout 3.300017V, fsw 500.00kHz. Peak simultaneous hs/ls conduction current: 1.16e-05 A, vs bound 1.00e-02 A | P_in=3.3956W P_out=3.3001W, **97.18%** (vs R4's 28.96%), energy-balance residual -0.0034% of P_in | 38.8s (ngspice-45.2, 6ms/3000 periods) | `tests/test_refs_d.py`: 4/4 checks pass |

**Design notes and traps (full detail in each deck's own header):**

- **Break-before-make, not a symmetric delay.** An RC delay applied
  identically to both a signal and its complement shifts both edges
  together and does NOT create a dead zone — a first R3d trial made
  exactly this mistake (both `pwm`'s complement generated as a plain
  `1-pwm` fed through the same symmetric RC) and `i(Vin)` showed a -210A
  spike, i.e. no real dead-time at all. Fixed with an asymmetric network:
  a slow RC charge path for turn-ON, a fast diode-steered discharge path
  for turn-OFF, on each side.
- **A "normal" ~0.6V-knee diode barely conducts in a 0-1V control signal.**
  The dead-time network's fast-discharge diode (`ddt` model) needed
  IS=1e-6 (Schottky-like, low knee), not a first trial's IS=1e-14 (a
  normal Si diode knee, enormous relative to the available ~1V swing) —
  with the normal-knee diode the "fast" discharge path barely conducted
  and the two switches' gate nodes still crossed 0.5 almost together
  (-159A spike).
- **A clean gate-crossing-time gap is not the same as a clean current
  overlap bound.** With the dead-time network's first working kswd=15
  (chosen from a hand estimate targeting a ~15ns switch transition), the
  two gate nodes' 0.5-crossings were genuinely ~20-25ns apart, but
  gon=70 is large enough that the tanh's residual tail still gave
  amp-scale simultaneous conduction current (up to ~0.78A, comparable to
  the ~1A load) — found only by independently re-deriving I_hs(t)/I_ls(t)
  from the trace using the switches' own controlling expressions (the
  same technique `tools/loss.py` uses for B-source devices), not by
  checking the gate signals' crossing times alone. Fixed with kswd=150,
  which makes the switch transition itself faster than the task's 10-20ns
  suggestion (see each deck's header for the trade-off this created).
- **A B-source current-source node (real capacitor + controlled current,
  high impedance) is not the same as a B-source voltage-defined node
  (zero impedance).** R3's `pwm` (from `Bcmp`, `V=expr`) can be loaded by
  the dead-time network for free. R4's `q` (a genuine RS-latch state, a
  real capacitor `Cq` driven only by a small controlled current `Bq`) is
  NOT — tapping the dead-time network directly onto `q` in a first R4d
  trial loaded it enough that the latch stopped oscillating at 500kHz and
  instead crept up over tens of microseconds, never switching. Fixed by
  buffering `q` through its own `V=expr` B-source (`Bqbuf`) first, and
  building the dead-time network off that zero-impedance copy instead.
- **`tools/loss.py` has no handling for ngspice's native `D` device.**
  It silently excludes an unhandled device kind from its sum (by design —
  see `docs/loss.md`), which is honest but meant the low-side body diode's
  real dissipated power (it carries the inductor current during every
  dead-time interval) was missing from the energy balance: a native-`D`
  version of the body diode gave a 0.5249% residual, just over the 0.5%
  bound. Fixed by rebuilding the body diode as a smooth B-source
  (soft-knee, `Bbody`, vfbody=0.6V, an effective 20 mOhm "on-resistance"
  above that) instead — `tools/loss.py` DOES handle `B, I=expr` devices
  (it re-derives their current from the parsed expression, the same way
  it validates Bhs/Bls, per `docs/loss.md`), so this is counted correctly
  without changing the diode's physical (smooth, continuous) behaviour.
  The dead-time steering diodes (`Ddis_hs`/`Ddis_ls`) are still native `D`
  devices — their own power turned out small enough not to matter.

**What was not done:**

- No independent re-verification (V-G4) — this card wrote and ran
  everything above.
- `tools/stb.py` (XP probe, Tian two-injection loop gain) was attempted on
  both decks and returned "no 0 dB crossing" — it runs its own DC
  operating point + `.ac`, and neither deck's DC fixed point is the
  periodic operating point the averaged small-signal loop-gain estimate
  assumes. This is not new to R3d/R4d: R3/R4's own table above already
  says no `stb`/loop-gain probe run was ever exercised on THEM either.
  FRA crossover/PM for R3d/R4d is therefore still only the unchanged
  averaged-model estimate from R3/R4 (same L/C/Esr/Rload, same
  compensator component values), with a hand argument (in each deck's
  header) that dead-time's effect on that estimate is second-order at the
  design crossover (a sub-degree extra phase lag from the ~15-25ns delay
  itself) — not a new measurement.
- Dead-time (~15-25ns, from the trace) and switch transition time (~1-2ns,
  see the kswd=150 trade-off note above) came out narrower than the
  task's 20-40ns / 10-20ns examples — correctness (no shoot-through,
  >85% efficiency) was prioritised over hitting those example numbers
  exactly; both are reported honestly rather than adjusted to match.
- Dead-time's own duty-cycle/gain-loss effect (a well-known real
  phenomenon in hardware buck converters) was not separately quantified
  here beyond the phase-lag argument above.
- Slope-compensation/duty-sweep verification (R4's own D=0.45/D=0.55
  cases) was not repeated on R4d — out of scope for this card, which kept
  `vslope_amp`/`vin_v` at R4's own default (D~0.30) values unchanged.
- `refs/r3_vm_buck.cir` and `refs/r4_pcm_buck.cir` themselves were not
  touched (per the task) and still show ~29% efficiency / genuine
  shoot-through — they remain the intended math references for Phase 4-5.

## Op-amp polarity fix (card #519, 2026-09-30)

Found in #505: the error-amplifier macro `Goa opout 0 vrefnode compin {gm_oa}` (R4/R4d: `Goa vc 0 ...`) in R3, R3d,
R4, R4d and R6 (and in every `tests/fixtures/*.cir` derived from them) had inverted polarity. With `C1f`/`C2f` from
`compin` to the op-amp output that is LOCAL POSITIVE feedback with a fast unstable pole (~gm/C ~ 1e12 1/s). Coarse
implicit steps land on the correct slow manifold, so every earlier FRA/PAC/PSTB/`.tran` number agreed with every other
one, but any small step blew up (fixed-step `.tran` -> 1e168 V, `reltol=1e-6` aborts, shooting monodromy runs died with
"Timestep too small"). Fixed: `Goa 0 opout vrefnode compin {gm_oa}` (`Goa 0 vc ...` in R4/R4d): current flows INTO the
output when vref > compin, i.e. an inverting amplifier with negative feedback through `C1f`/`C2f`. 65 files changed
(5 refs, 60 fixtures); `bench/` generators and decks carry the same old line and were NOT changed (speed baselines).

Confirmation on R3 (s69 build, `tmax` = step where stated): `.tran 1e-14 4e-11 uic` old v(opout) = -7.2e173 V at 40 ps,
new 1.0 V; `.tran 1p 5u 0 1p uic` old aborts at 27.8 ps ("Timestep too small", 4.6e277 V), new runs to 5 us;
`reltol=1e-6 vntol=1e-9` full 1.5 ms run old aborts at 102 ps, new completes (Vout 3.29998 V over 1490-1500 us);
`tran 1n 1.5m 0 1n uic` and the default step both run for old and new (Vout 3.30002/3.30012 old, 3.29998/3.29995 new).
The converter still regulates at the same Vout.

Steady-state numbers, old -> new (each deck run from `refs/` with stock 45.2 and with the s69 fork build; last period /
last 10 us windows as in the decks' own `.meas`):

| deck | binary | Vout final avg (V) | ripple (mV) | duty (final) | fsw (kHz) |
|---|---|---|---|---|---|
| R3 | 45.2 | 3.29958 -> 3.29947 | 4.836 -> 5.000 | 0.2777 -> 0.2762 | 500.30 -> 499.25 |
| R3 | s69 | 3.30012 -> 3.29995 | 4.881 -> 4.941 | 0.2782 -> 0.2764 | 499.85 -> 500.08 |
| R3d | 45.2 | 3.30001 -> 3.29983 | 4.931 -> 4.929 | 0.2870 -> 0.2872 | 499.95 -> 499.90 |
| R3d | s69 | 3.29997 -> 3.29988 | 4.994 -> 4.869 | 0.2867 -> 0.2874 | 500.00 -> 500.20 |
| R4 D~0.30 | 45.2 | 3.30002 -> 3.29998 | 4.651 -> 4.651 | 0.3055 -> 0.3055 | 500.04 -> 500.04 |
| R4 D~0.45 | 45.2 | 3.30002 -> 3.29998 | 3.619 -> 3.619 | 0.4583 -> 0.4584 | 500.00 -> 500.00 |
| R4 D~0.55 + slope 0.15 | 45.2 | 3.30002 -> 3.29998 | 2.935 -> 2.935 | 0.5600 -> 0.5600 | 500.00 -> 500.00 |
| R4d | 45.2 | 3.30003 -> 3.29997 | 4.690 -> 4.693 | 0.3048 -> 0.3058 | 500.00 -> 500.00 |
| R4d | s69 | 3.30002 -> 3.29998 | 4.692 -> 4.695 | 0.3049 -> 0.3062 | 500.00 -> 500.00 |
| R6 | 45.2 | 3.30002 -> 3.29998 | 4.811 -> 4.811 | 0.2712 -> 0.2711 | 497.85 -> 497.85 |

The R3/R3d scatter (duty 0.3-0.6 %, ripple 1-3 %) is the same size as the old binary-to-binary scatter of the same old
deck (45.2 vs s69: duty 0.2777 vs 0.2782, ripple 4.836 vs 4.881 mV, fsw 500.30 vs 499.85 kHz): it is the smooth-switch
converter's sampling/step-selection noise, not a change of the circuit. R4 is identical to 1e-5. R4d's duty +1e-3 is the
sharp dead-time edge (kswd=150) landing on a different time step. R4 D~0.55 without slope compensation is period-2
unstable both before and after (the per-clock-period on-fraction of `q` alternates between ~0.1-0.3 and ~0.95 in both; the deck's `t_period` meas, which
pairs `q` rises 250 and 251, reads 4.00 us old and 2.00 us new only because the pairing phase differs -- it is a
measurement artefact, the orbit is the same). What DID change for the better: the shooter now converges to that
unstable period-1 orbit and flags it `UNSTABLE(period-doubling)` (multiplier -1.316); before, it died in the monodromy
runs. FRA reference curves regenerated (s10 fork, `fra/curves-fork-s10-polfix/`, and s37 for R3d/R4d,
`fra/curves-fork-s37-rd-polfix/`; the R3d/R4d set replaced `fra/curves-fork-s37-rd/`, the old files are in
`fra/curves-fork-s37-rd-oldpol/`): loop PM/crossover R3 38.46 deg @ 33463 -> 38.46 @ 33457 Hz, R4 D~0.30 89.25 @ 5631 ->
89.25 @ 5631, R4 D~0.45 89.37 @ 5646 -> 89.37 @ 5645, R3d 38.74 @ 33489 -> 38.72 @ 33451, R4d 89.26 @ 5625 -> 89.25 @
5625. R4/R4d loop, control-to-output and output-impedance curves agree to <=0.09 dB / 0.6 deg at every point (R4
line-to-output 0.04 dB; R4d line-to-output differs up to 5 dB at -75 dB magnitude, where that FRA reference is unreliable
anyway). R3/R3d differ only within the settle-window noise of the method: R3 loop worst 0.15 dB / 1.2 deg at the 5 kHz
resonance, the same spread as re-running the OLD deck with settle 300/450/600 us (35.4-35.7 dB, -55 to -57 deg); R3d at 5 kHz
swings by +-3 dB / +-20 deg with settle time for old and new decks alike. The stored tables in `tests/` (FRA_R3_*, FRA_R4_*)
were left as they are (the new values are inside that noise; td-PAC/PSTB still agree with them, `tests/test_pac.py` green).
R3d/R4d were never checkable against `.pac`/`.pstb` before (the shooter could not converge on r3d or crashed in its
monodromy runs); now they run, and R3d's PSTB phase margin (38.6 deg @ 33.48 kHz) matches FRA and R3.

Shooter side effects of the fix: (1) `.pss` warm-up probe stages that used to abort instantly (blown-up monodromy runs) now
do real Newton work, so `r3_longsettle` shoot is 30 % slower (0.45 -> 0.59 s, .tran unchanged; ratio 11.1x -> 8.5x, below
the 10x gate in `tests/test_pss.py`); (2) r3d shooting convergence depends on the `.pss` stabtime (converges at 400-600 us,
fails at 300 and >=700 us for `rd_r3d_shoot.cir`; `tfp_r3d_loadstep` converges at 700 us, not 500 us) -- the r3d fixtures use
500 us (`rd_r3d_shoot.cir`) / 700 us (`test_tran_from_pss.py`) now.

## R7 / R7p -- transistor-level (BSIM4) synchronous buck, card #539, 2026-10-01

`r7_bsim_buck.cir` (run `ngspice -b r7_bsim_buck.cir` from `refs/`) is the first reference deck in which everything that
conducts or decides is a MOSFET: BSIM4 (level 54, v4.8.2) power FETs, BSIM4 gate drivers with a gate-delay dead-time generator, a
BSIM4 diff-pair PWM comparator, a BSIM4 5T OTA error amplifier with type-II compensation. Only Vin, vref, the OTA bias current, the
ramp, the passives and the Tian probe are ideal. All model cards are copied from cards already in the repo
(`tests/fixtures/loss_bsim4_buck.cir` for the power FETs, `bench/decks/d_analog_block.cir` for logic/analog; capmod=0 as there).
`r7p_bsim_buck_pex.cir` (made by `refs/gen_r7p.py`, which reuses `bench/gen_pex_buck.py`'s `tree_lines`) adds three synthetic
RC trees (on vout, sw and the HS gate net): 399 capacitors + 480 resistors, 413 explicit states in all.

| item | value |
|---|---|
| power stage | 5 V -> 1.46 V, 2 A (Rload 0.75), L 2.2 uH + 10 mOhm, Cout 100 uF + ESR 80 mOhm (ESR zero ~20 kHz), fs 500 kHz |
| switches | HS PMOS 60 mm/0.35 um, LS NMOS 40 mm/0.35 um (junction as/ad/ps/pd set so body diodes and Coss exist) |
| dead time | NOR pair with 6-inverter 3 pF-loaded delays; measured at the gates 18.76 ns (HS off -> LS on), 18.84 ns (LS off -> HS on) |
| control | NMOS-pair comparator vs 1.2..2.2 V ramp; OTA gm 1 mS, ro 49 kOhm; Rc 2k, Cc 6n, Cp 20p; divider 1.5k/1k, vref 0.6 (OTA inn = fb, inp = vref: opout falls when Vout rises, #519) |
| steady state (`.tran`, reltol 1e-5) | Vout 1.45940 V (DC gain error of the 34 dB OTA; target 1.5), ripple 70.15 mV, il 1.9464 A (pp 0.970), fsw 500.000 kHz, gate duty 0.311, periodic to 3e-7 of ripple; efficiency 92.9 % (own Pin/Pout over 10 periods, drivers and control included; `tools/loss.py` does NOT apply: it cannot instrument MOSFETs inside `.subckt`, "no such device c1") |
| loop (converged FRA, `fra/curves-r7-20261001/loop`) | crossover 26.19 kHz, PM 53.5 deg (11 points, unc <= 0.09 dB / 0.1 deg per point; the averaged-model design target was ~55 deg) |

### What was run and the verdict (engine `aaba50ba3f49-s88`, development workstation (AMD Ryzen AI Max+ 395, Linux) load 5-19, OMP 2, <= 4 ngspice)

| analysis | result | verdict |
|---|---|---|
| `.tran` | table above; reltol 1e-3 .. 1e-5 and tmax 2n .. 0.5n move Vout by <= 1e-6 V | pass |
| `.pss method=shoot` | needs `.options reltol <= 3e-4` at Points 400 (default 1e-3 stalls at residual 4e-4; Points 1000 also converges, 85 s). iter 3-4, rho 0.893 (tran envelope decay 0.899), 25-30 s. Orbit vs settled tran: vout +0.36 mV mean offset (0.57 % of ripple max), il +0.48 mA (0.10 %), edge times agree to 0.3 ns; the period restart re-initialises every non-state node (gate chain, logic, MOS charges): first ~15 ns of the orbit wrong (opout spike 47 % of its ripple at 1 ns, 8 % at 50 ns) | pass, with the restart artefact stated (PSS state = native C/L only: 14 states here) |
| `.pss` klu / krylov=1 / both | same orbit to <= 5e-6 of ripple, rho 0.8917 .. 0.8930 (finite-difference spread) | pass |
| R7p (413 states) | `krylov=1 klu`: 36-48 s, rho 0.8933; krylov alone 63 s; **dense monodromy (plain and klu) is OOM-killed at the 8 GB cap** (not run) | pass for krylov; dense not runnable under the cap |
| OpenMP | `num_threads` 1 / 2 / 8: tran byte-identical, `.pss` orbit + iterations + residual byte-identical (1 thread 36 s, 8 threads 19 s) | pass |
| `.pstb` (loop gain) | |T| ~ 2e-10 at all of 500 Hz .. 200 kHz ("no 0 dB crossing") vs FRA 30 .. -21 dB | **FAIL, engine bug 1** |
| `.pac` ctrl-to-out / line-to-out / Zout | ctrl: -129 .. -121 dB vs FRA -17.6 .. +3.6 .. -18.2 dB; line: +25 .. -1 .. -3.4 dB error; Zout: +30 dB error at 500 Hz, <= 2 dB above 20 kHz (full table `bugs/bug1_pac_model_charge/r7_td_vs_fra.txt`) | **FAIL, engine bug 1** |
| `.pxf` forward vs adjoint | forward == `.pac` exactly; adjoint differs: Iinj 0.01-0.27 dB, Vin constant -24 dB | FAIL (forward is itself bug 1; adjoint additionally bug 3) |
| `.pnoise solver=td dev=1` | Rload and Resr thermal contributions at 100/200 kHz match the exact LTI value to 0.14-0.2 % (they reach vout through passives only). MOS channel noise, loop-shaped and folded parts NOT validated: no MC reference (`tools/pnoise_mc.py` handles resistors only, BSIM4 has no transient noise) and the engine's PAC recursion is bug 1 | sanity only, not validated |
| `.stb ... instant` at 300.5/300.9 us | "no 0 dB crossing" (mid-ramp comparator saturated, frozen-instant loop gain ~0): not applicable to a switching loop | not applicable |

Independent references (all agree where they run): converged transient FRA (`tools/fra.py`, loop and transfer, per point settle x2 /
amplitude/2 (/ reltol 1e-5 at 3 points) variants: unc <= 0.12 dB / 0.34 deg, 500 Hz worst) and `tools/pac_perturbed.py` (reltol 1e-4):
ctrl-to-out 20 kHz 3.623 dB / -7.78 deg vs FRA 3.627 / -7.79; line-to-out 20 kHz -16.202 / -8.73 vs -16.193 / -8.77; loop gain
20 kHz 4.2910 dB / -129.61 deg vs 4.2904 / -129.63, 80 kHz -12.828 / -107.12 vs -12.831 / -107.11; ctrl-to-out 80 kHz -9.687 / -84.41 vs -9.690 / -84.39, 200 kHz -18.166 / -90.87 vs -18.161 / -90.88 (each point 16-24 min, so only a few points).

### Engine findings (reproducers + diagnosis in `~/Work/data/ngspice_augmentation/539/bugs/`)

1. **`solver=td` .pac/.pstb/.pxf/.pnoise drop the history of device-model charge storage** (`bug1_pac_model_charge/`). The reduced state is native
   C/L only (`pacshoot.c`: history terms built from `CAPcapac`/`INDinduct`), while `J_k = G + ag0*C` contains the model capacitances. A chain of N
   inverters loses its timing sensitivity per stage (0.55, 0.31, 0.04, 0.001 for BSIM4; the finer the grid the worse); level-1 MOS without caps is
   right, level-1 with `cgso/cgdo` is wrong. Silent: no warning. Not BSIM4-specific; hits every transistor-level converter.
2. **FIXED (card #541, 2026-10-02)** `.tran` aborted "Timestep too small ... vramp#branch" (`bug2_tight_reltol/`) for many tight option sets (reltol 1e-6 / 1e-7,
   vntol 1e-7, ... deterministic, always at a multiple of the period; the node named varies: vramp, vrefsrc, vin). Cause: the DECK, not an engine patch. The
   ramp `PULSE(.. td=0 tr={0.97*tclk} tf={0.03*tclk} pw=1n per={tclk})` has tr+pw+tf = tclk + 1n, so the next rising edge starts 1 ns before the falling edge
   ends; ngspice then places breakpoints 1 ns apart at every period boundary and the step collapses (1.25e-21 s). Stock `/usr/local/bin/ngspice` 45.2 aborts
   identically at 4 us, and so does `aaba50ba3f49-np1-nopatch` (no patches); with tf shortened by 1 ns every option set tried runs (reltol 1e-5/1e-6/1e-7,
   vntol 1e-7/1e-9, abstol 1e-9), so the NGSPICE_RESBASE/CAPBASE/FLATB/KLU_REFSCALE switches were not needed. `refs/r7_bsim_buck.cir` now uses
   `tf = {0.03*tclk-1n}` (rise identical, fall 59 ns instead of 60); the fixtures `tests/fixtures/r7_*.cir`, R3/R4 decks and handbook examples keep the old line
   (same trap, only hit at tight options; golden numbers unchanged). `tools/pac_perturbed.py` repairs such PULSE lines itself (`fix_pulse_overlap`, prints a
   note), so its default `--reltol 1e-5` is safe again. Rule for any deck: PULSE needs tr + pw + tf <= per.
3. **FIXED (card #541)** `.pxf adjoint=1` ~0 across a MOS gate (`bug3_pxf_adjoint_mos_gate/`): the native-path backward sweep re-loaded J_k in reverse order from
   device state, and the forward "exact" J_k (`PAC_EXACT_JACOBIAN`, forced by the adjoint engine) was loaded in `MODETRAN|MODEINITPRED`, i.e. at the PREDICTED point;
   now loaded in INITFLOAT at x_k and stored per step for the transpose. See docs/pss-pac.md "Native adjoint" (card #541).
4. **FIXED (card #541)** `.pss` refusing a PWL clock: `VSRCparam()` clears `VSRCfuncTGiven` on every parameter call, so `PWL(...) r=0` (r after the function) lost the flag;
   the driven-circuit scan in `dcpss.c` now also accepts `VSRCfunctionType != 0`. Dense shoot on 413 states exceeding 8 GB is unchanged.

### What was not done

- No fix to the engine (card scope). No `.pac/.pstb` validation is possible until bug 1 is fixed or refused; the checks in `tests/test_r7.py` for them FAIL on purpose
  (`--no-known-bad` skips them). No MOS1-without-caps copy of R7 as a positive control (a first version collapsed to Vout 0; not debugged): the positive control is
  the cap-less inverter chain in bug 1.
- FRA grid: 8 frequencies + 3 near crossover for the loop (500 Hz .. 200 kHz), 8 for each transfer; pac_perturbed only at 20/80 kHz (loop), 20 kHz (line, ctrl), 80/200 kHz (ctrl).
- The orbit mean offset (+0.36 mV) is attributed to the period restart; the attribution was not isolated (a phase-shifted clock needs a non-PULSE source, which `.pss` refuses).
- Comparator metastability near the PWM edge: only that the steady state has exactly one pwm edge pair per period at default tolerance; no duty-extreme or noise-driven test.
- No saltation question arises (no `sw` events); gate-level events are smooth and resolved by ns steps.
- Slot/time: R7 capmod=0 only (repo cards); capmod 2/3 BSIM4 charge models untested. `pnoise` MC cross-check infeasible (see table).
- Tests: `tests/test_r7.py` 8 pass / 4 FAIL (the four known-bad), ~4 min; `tests/fixtures/r7_*.cir`, `r7p_*.cir`.

## R8 -- transistor-level (BSIM4) free-running hysteretic synchronous buck, card #544, 2026-10-02

`r8_bsim_hysteretic_buck.cir` (run `ngspice -b r8_bsim_hysteretic_buck.cir` from `refs/`; made by `~/Work/data/ngspice_augmentation/544/mkr8.py`): R7's BSIM4 power stage,
drivers, dead-time generator and comparator, with the ramp, error amplifier and loop probe removed; hysteresis by resistive positive feedback from `q` (Rh1 100k / Rh2 1k,
50 mV window at the comparator), the comparator output drives a native `Cq` (2 pF) through `Rq` (1k) so the hysteresis memory is a native state. 5 V -> 1.7055 V, fsw 395.5 kHz
(settled `.tran` period 2.5266 us; `.pss trigger=cx` 2.5284 us), ripple ~95 mV. Fixtures `tests/fixtures/r8_pss.cir`, `r8_pac_ctrl2out.cir`, `r8_pac_line2out.cir`.
**Trigger**: `trigger=cx trigval=1.7055 trigdir=rise` (Cout's positive node). `trigger=q` does NOT work: the comparator output is hidden state that a native-state restart resets
and `q` re-crosses 20 ns after the restart. `.pac solver=td` through the generic model-charge path vs `tools/pac_perturbed.py --freerun`: <= 0.021 dB / 0.052 deg at 20/50/80 kHz
(docs/pss-pac.md, card #544).

## R9 -- digitally controlled synchronous buck from XSPICE primitives, card #551, 2026-10-02

`r9_digital_buck.cir` (GENERATED by `tools/gen_r9.py`; `python3 tools/gen_r9.py --nbits 7` makes the limit-cycling variant, `--inject` adds a summer + `Vinj` in the DPWM
control path, `--lgprobe` sets the include path; run `ngspice -b r9_digital_buck.cir` from `refs/`: 60 ms settled `.tran`, ~90 s). Power stage = R3 (12 V -> 3.3 V, L 10 uH,
C 100 uF, ESR 10 mOhm, Rload 3.3, smooth tanh switches, 500 kHz, `lgprobe` XP in the feedback divider). Controller, all XSPICE: two `adc_bridge` comparators on the feedback
node (fb > vref+bin, fb < vref-bin; bin 3 mV, i.e. +-9.9 mV on Vout: a 3-level quantiser with a zero-error bin) sampled by `d_dff` on the falling clock edge; a slow enable
clock `d_fdiv` (/600 = 1.2 ms); an 11-bit synchronous up/down counter (`d_tff` + `d_and`/`d_or` carry logic, one LSB per enable edge while the sampled error is outside the bin)
= the digital integrator; DPWM = multi-bit `dac_bridge` (code 0..1) -> `d_pwm` (duty = 0.01 + 0.98 code, 500 kHz) -> `dac_bridge` -> switch drive. Initial duty 0.2745 (counter ic).
Why the slow step: the plant is a Q ~ 10 LC (5 kHz, ring time constant 0.66 ms); with a step every 96 us the sign-based integrator rings the output across the bin and hunts
+-1 LSB (measured, 48 clocks); with 1.2 ms between steps every step has settled before the next sample. No-limit-cycle condition (Peterchev/Sanders): DPWM step < ADC zero-error bin:
11 bits = 5.7 mV at Vout (1.7 mV at fb) < 6 mV bin width; `--nbits 7` = 92 mV >> bin: the integrator hunts between two codes (period 2 x 1.2 ms, Vout 3.225 / 3.315 V in alternate
windows) -- the second variant, kept as the refusal/diagnostic case. Fixtures `tests/fixtures/r9_fine.cir`, `r9_lc.cir`, `r9_inj.cir`; suite `tests/test_xspice.py` (section R9).

Numbers (build aaba50ba3f49-551-x, development workstation, OMP 2): settled tran: code 556 from 6 ms on (constant), Vout mean 3.2926 V (in the bin), ripple 3.72 mV. `.pss method=shoot`, T = 1.2 ms
(600 cycles, Points 24000, stabilisation 12 ms): converged in 3 iterations, 152 s, spectral radius 0.029; orbit vs the last 1.2 ms of the settled tran (100 ns grid): max 0.012 mV
(0.33 % of ripple), mean offset +0.005 mV, DPWM code 556 in both; `XSPICE hidden state ... digital node mismatches 0`. With `PSS_XS_NO_STATE=1` (register back at ic 553 every period):
orbit 17.15 mV low (4.6 x ripple). The 7-bit variant: `.pss` at 1.2 ms does not converge (residual 0.55) -- no periodic orbit with that period.
FRA / loop gain: control-to-output (DPWM code -> Vout, `tools/fra.py --mode transfer --source vinj`, `tests/fixtures/r9_inj.cir`, 24 ms settle, 20 periods, `--tmax 1e-7`): 22.81 dB / -4.1 deg at 2 kHz, 35.84 dB / -83.1 deg at 5 kHz,
11.95 / -168.5 at 10 kHz, -1.89 / -168.7 at 20 kHz, -17.84 / -157.9 at 50 kHz, amplitude-INDEPENDENT (injection 2.5 / 10 / 40 mV of the 0..1 code: 22.82 / 22.81 / 22.79 dB, 11.91 / 11.95 / 11.94, -2.00 / -1.89 / -1.99
at 2 / 10 / 20 kHz, spread <= 0.11 dB, 0.6 deg) -- d_pwm is continuous in its control. Against the averaged LC model (K = 0.98 x 12 V, L + Rdcr, Cout + ESR, Rload; no switch Ron, no sampling delay):
-0.04 dB / -0.5 deg at 2 kHz, -0.01 / +1.8 at 10 kHz, +0.07 / +1.7 at 20 kHz, +0.17 / +3.7 at 50 kHz, -1.57 dB at the 5 kHz resonance (Q ~ 10: the switch Ron, 10 mOhm each, is not in the model).
LOOP GAIN (`--mode loop --probe XP --startup uic --settle 24e-3`, fb-node injection, 300 / 1000 / 3000 Hz): NOT a small-signal quantity here. |T| at 300 / 1000 / 3000 Hz in dB:
0.5 mV: -142 / -151 / -21; 2 mV: -15 / -157 / -33; 5 mV: -87 / -118 / -52; 20 mV: -17 / -26 / -30 (phase -90 .. -113 deg: the integrator). The injection is smaller than or comparable to the 6 mV zero-error
bin and the comparator output is sampled every 2 us and applied every 1.2 ms: the response is zero, aliased or amplitude-dependent by 100 dB; only a large injection (20 mV, > bin) shows the sign-integrator's
describing function (1/(A f) slope). So no PAC reference exists for this loop, `.pac/.pstb solver=td` REFUSES the deck (adc_bridge / d_tff / dac_bridge named), and the FRA table above is the diagnostic that
shows why (data `~/Work/data/ngspice_augmentation/551/r9/fra_loop_*`, `c2o/`). `tools/pac_perturbed.py` was started on `r9_inj` (LCM-period backend) and killed after 40 min without a point: not validated.
Not done: loop-gain FRA above 3 kHz (the digital loop has nothing there), the 7-bit variant's FRA, FRA/pac_perturbed uncertainty beyond the amplitude spread.
