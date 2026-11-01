# Adjoint sensitivity in ngspice_augmentation

Card #475 (plan T10.7), FIRST HALF, 2026-09-29. Adjoint sensitivity of
native `.stb`'s loop gain T(f), phase/gain margin and crossover frequency,
about a fixed DC operating point. The SECOND half (adjoint sensitivity on
top of the shooting-PSS/periodic-AC orbit, `pacshoot.c`'s recording pass)
is separate, later work and not covered here.

Private-fork worktree `ngspice-src-sens`, branch `aug-sens` = pin +
`patches/0001`-`0016` (promoted `.stb`/numparam/accel series) +
`patches-wip/sens/0400`-`0403` (this card). Built at
`~/Work/apps/ngspice_augmentation/aaba50ba3f49-sens6/bin/ngspice-aug`.

**Update 2026-10-11 (card #595, patch 0133, s133): the operating point is no
longer held fixed.** `dT/dp` and the margin sensitivities now include the
shift of the DC operating point, and for nonlinear devices every quantity the
device computes in its DC load (gm, gds, charges). Parameters that move the
bias (a MOSFET's `w`, a bias current, a model `vto`) were wrong before, by
40-200 %. See "The bias-point term and nonlinear devices" below; the sections
in between describe the card #475 design they build on.

## What this computes

For native `.stb`'s Tian two-injection loop probe -- one complex matrix
factorisation per swept frequency, two back-solves giving `(v1,i1)` and
`(v2,i2)`, combined into `s = 2*(i1*v2 - v1*i2) + v1 + i2`, `T = s/(1-s)`
(`docs/stb.md`, `src/spicelib/analysis/stban.c`) -- `sens=all|<name>` adds,
per requested parameter `p`:

- `dT_d<name>` (a complex vector vs frequency, same plot family as
  `loopgain`/`t`, output plot title `"STB dT/dp"`)
- `dPM_d<name>`, `dfc_d<name>`, `dGM_d<name>` (real scalars, output plot
  title `"STB dMargins/dp"`) -- the phase margin, crossover frequency and
  gain margin sensitivities, at the SAME worst-margin crossing(s) native
  `.stb`'s own `STB Margins` plot already reports. `NaN` when the
  underlying margin has no crossing (never a fabricated number).

## Syntax

```
.stb dec|oct|lin <n> <fstart> <fstop> probe=<probespec> [sens=all | sens=<name> [sens=<name> ...]]
```

`sens=` is loop-mode (single-ended) only -- refused (`E_PARMVAL`) together
with a device-mode probe or `probe2=`/`mode=dm|cm`. It also needs the
Sparse solver (`.options sparse`, or simply not passing `.option klu`):
the complex transposed solve this needs (`SMPcaSolve`) and the sparse
matrix-vector product (`SMPmultiply`) are Sparse1.3-only in this fork;
refused, not silently wrong, under KLU.

`<name>` is a parameter name in the SAME form
`src/spicelib/analysis/cktsens.c`'s existing DC/AC `.sens` analysis uses
(and produces, if you run `.sens` on the same circuit): the bare instance
name for a device's principal parameter (`R1o` for a resistor's own
resistance, `Cc` for a capacitor's own capacitance), or
`<instance>_<keyword>` for anything else (`M2_w` for a MOSFET's channel
width). `'*'`/`'?'` glob wildcards are accepted (`sens=R*`). Repeat
`sens=` for more than one parameter -- a comma-separated list was
deliberately NOT used (`sens=Vp,Vn`-style syntax on a netlist card or the
interactive `stb` command both risk exactly the tokenizer comma-splitting
hazard `probe=Vp,Vn` already hit, card #457; `probe=<Vp> probe2=<Vn>` was
that card's own fix, and `sens=` follows the same pattern).

`sens=all` matches every `IF_REAL`, `IF_SET|IF_ASK`, AC-capable
instance/model parameter of every device in the circuit (R, C, L,
VCCS/VCVS gains, and MOSFET model parameters such as `w`/`l`/`vto`) --
the same generic walk (`sgen_init`/`sgen_next`, `cktsgen.c`) `.sens`
already uses, reused here rather than a fresh, device-type-specific list.

## Design

For a small-signal reading `y = c^T x` at one frequency, with `Y(p) x = b`
(`b` independent of `p` here -- the probe's own unit injections):

```
dy/dp = lambda^T (db/dp - (dY/dp) x),   Y^T lambda = c
```

**One extra TRANSPOSED solve per OUTPUT reading** (`SMPcaSolve`, the SAME
factorisation `SMPcLUfac` already produced for the forward solves -- no
new factorisation), reused across every requested parameter. Four probe
readings (`v1,i1,v2,i2`) but only TWO distinct `c` vectors --
`v1`/`v2` share `c_v = e_{node_b}`, `i1`/`i2` share `c_i = -e_{branch_probe}`
-- so exactly **2 adjoint solves per frequency**, independent of how many
parameters are swept.

`(dY/dp) x` for one isolated device/parameter is built with
`cktsens.c`'s OWN proven technique: `DEVsetup` into a throwaway,
per-parameter matrix (an isolated single-instance model, via `sgen`'s own
temporary instance-unlinking trick), then `DEVacLoad` at `p` and at
`p+delta` (finite difference of the device's own analytic stamp, not of
the circuit's solution), subtracted. This is not new device-stamp math --
every device type `.sens` already covers is covered here for free, with
no per-device-type code in `stban.c`.

**The cheap part**: given the throwaway `dY/dp` matrix, `(dY/dp) x` is
ONE sparse matrix-vector product (`SMPmultiply`) -- no solve. `dv/dp`,
`di/dp` are then a plain dot product against the precomputed `lambda_v`/
`lambda_i` (`stb_sens_dot()`, an ordinary complex "bilinear" sum, NOT a
conjugated inner product -- the adjoint identity needs `Y^T`, not `Y^H`).
Chain rule through Tian's formula (`s`, `T = s/(1-s)`) gives `dT/dp`.

**Margins**: `dPM/dp` and `dfc/dp` come from implicit differentiation of
the crossing condition `|T(fc)| = 1` through the SAME bracket (the two
swept frequencies either side of the crossing) and log-frequency
interpolation fraction the crossing search itself used:

```
d ln|T| /dp + d ln|T| /dx * dx/dp = 0     (x = ln f)
dx/dp = -(d ln|T|/dp) / (d ln|T|/dx)
dfc/dp = fc * dx/dp
dPM/dp = d(arg T)/dp + d(arg T)/dx * dx/dp
```

`d(.)/dp` uses the exact, per-frequency-point `dT/dp` (interpolated to
the crossing with the same fraction as `T` itself); `d(.)/dx` is a plain
bracket finite difference of the SWEPT `T(f)`/`dT/dp(f)` traces --
**this is the one place this implementation is not exact**: its accuracy
is set by the sweep's own frequency resolution near the crossing, not by
the adjoint formula (which is otherwise accurate to floating-point
noise -- see Validation). `dGM/dp` needs the SAME two-term treatment,
since the `-180deg` gain-margin crossing's own frequency shifts with `p`
too (`d(arg T)/dp + d(arg T)/dx * df_gm/dp = 0` first, then the value's
own chain rule) -- an earlier version of this code computed only the
first ("evaluated at a fixed f") term and was wrong by O(1); central-FD
validation caught it (see "Bugs found", below).

## What is / is not included

- **Superseded by card #595 (below) -- kept as the record of card #475.**
  **The DC operating point is treated as FIXED** (linearised circuit
  only): `dx_op/dp` (how the bias point itself shifts with `p`) is NOT
  computed or folded in. For a parameter that does not affect the DC
  operating point at all (a pure AC-coupling/load capacitor, for example
  -- no DC current path), this is not an approximation, it is EXACT: the
  small-signal matrix `Y` genuinely does not depend on `p` through any
  bias-point term. For a parameter that DOES move the bias point (a
  MOSFET's `w`, for instance -- shifts `Id`, hence `gm`/`gds`/every
  bias-dependent capacitance of every transistor sharing that bias
  branch), this implementation's `dT/dp` is missing a real term and WILL
  disagree with a true (bias-point-included) finite difference -- see
  `tests/test_sens.py`'s `M2_w` case, deliberately kept in the test suite
  (labelled `EXPECT-DIVERGE`, not a FAIL) as an honest, checked
  demonstration of this boundary, not a hidden gap.
- Loop mode only (single-ended `probe=<Vname>` or `probe=<inst>:<terminal>`).
  Device-mode probes (`probe=<MOS1/BSIM4 instance>`) and differential/
  common-mode (`probe2=`/`mode=dm|cm`) are refused together with `sens=`.
- Sparse solver only (KLU refused, see Syntax).
- Element/model PARAMETER sensitivity (R, C, L, gain, W/L/model
  parameters), not device-level small-signal QUANTITY sensitivity
  (treating `gm`/`gds` themselves as independent knobs, which the card's
  own brief floated as a "device parameters first" option). Real settable
  parameters are what `.sens`'s existing, generic, proven mechanism
  already reaches, so that is what this reuses; `gm`/`gds` are
  operating-point-DERIVED, not independent circuit unknowns, and treating
  them as if they were would need new, unvalidated math this pass did not
  attempt.
- Pause/resume (`SPfrontEnd->IFpauseTest()` mid-sweep) is refused
  together with `sens=` -- the per-parameter sensitivity arrays are plain
  local variables in `STBan()`, not the function-static `stbPlot` pattern
  the rest of the analysis uses for resuming, and were not made to
  survive a pause.

## Worked example: which components most affect PM of the R3 loop

`tests/fixtures/sens_r3_averaged_native.cir` (native-probe copy of
`tests/fixtures/compopt_r3_averaged_detuned.cir`, the spec's own R3
DC-op-half reference -- a state-space-averaged linear stand-in for the
real switching R3 buck's type-III-compensated voltage loop, `docs/stb.md`
and that fixture's own header for the averaging derivation). At its
DETUNED starting point (`r2=100k`, `c2=0.3n`), native
`.stb sens=R2f sens=C2f sens=Resr` (`dec 300 3e4 8e4`, bracketing the
`~55.8` kHz PM crossing) gives, at the crossing:

| Parameter | Role | `dPM/dp` | `p0` | `dPM/d(ln p)` (deg per e-fold) | deg per 1% |
|---|---|---|---|---|---|
| `R2f` (compensator mid-band gain resistor) | sets the type-III compensator's flat-gain region | `-2.860e-4 deg/ohm` | 100 kohm | **-28.6** | -0.286 |
| `Resr` (output cap ESR) | sets the power-stage ESR zero | `+1667 deg/ohm` | 10 mohm | **+16.7** | +0.167 |
| `C2f` (compensator high-freq pole cap) | sets the type-III compensator's HF pole | `+6.089e9 deg/F` | 0.3 nF | **+1.83** | +0.018 |

**`R2f` dominates** -- about 16x `Resr`'s sensitivity and about 16x
`C2f`'s, at this operating point. This matches hand intuition for a
type-III compensator: `R2f` sets the mid-band gain that directly
determines how far the loop's 0 dB crossing sits from the compensator's
poles/zeros, so a small change in it moves the crossing frequency (and
therefore the phase there) more than an equally-sized change in either
reactive element. `Resr`'s comparatively large effect (more than `C2f`'s)
reflects that the ESR zero sits close to this design's crossover already
(`compopt_r3_averaged_detuned.cir`'s own header notes the ESR zero and
crossover interact) -- a genuine, not accidental, finding, not merely a
units artefact (the "per 1%" column normalises that out).

Reproduce: `NGSPICE=<prefix>/bin/ngspice-aug python3 tests/test_sens.py`
(this table's numbers are this script's own `R3avg` red-proof/FD-compared
run, re-derivable from its printed `adjoint=...` values times each
parameter's nominal value).

## Validation

`tests/test_sens.py`, against
`~/Work/apps/ngspice_augmentation/aaba50ba3f49-sens6/bin/ngspice-aug`
(pin + `patches/0001`-`0016` + `patches-wip/sens/0400`-`0403`):
**27/27 gating checks pass** (0 FAIL, 0 COULD-NOT-RUN; 2 `skip` -- R1d has
no gain-margin crossing, so `dGM_dCL` has nothing to compare against; 5
`info` -- the `M2_w` op-point-dependence demonstration above, and the
`DELTA` stability sanity check, neither gates pass/fail) plus **2/2
red-proofs pass** (see below; these need two throwaway builds via
`NGSPICE_REDPROOF_YT`/`NGSPICE_REDPROOF_SIGN` env vars, not run by
default in the same invocation as the main suite -- see the script's
`check_red_proofs()`).

**R1b** (`refs/r1b_bilateral.cir`'s topology, native-probe copy
`tests/fixtures/sens_r1b_native.cir`, linear): `R1o`, `Cc`, `Rf1` --
`dT/dp` matches central FD to **~1-2e-6 relative** at every sampled
frequency (well inside 1e-4); `dPM/dp`, `dfc/dp`, `dGM/dp` match to
**3e-5 to 5e-5 relative**.

**R3 averaged linear stand-in** (`tests/fixtures/sens_r3_averaged_native.cir`,
the spec's own DC-op-half reference for R3): `R2f`, `C2f`, `Resr` --
`dT/dp` matches to **~1-2e-6 relative**; `dPM/dp`, `dfc/dp`, `dGM/dp`
match to **1.5e-6 to 3.6e-5 relative**.

**R1d nonlinear** (`refs/r1d_mos_follower.cir`'s topology, native-probe
copy `tests/fixtures/sens_r1d_native.cir`, level-1 MOSFETs): `CL` (the
op-point-INDEPENDENT load capacitor) matches central FD to **~1e-6
relative** on `dT/dp`, `dPM/dp`, `dfc/dp` (no GM crossing on this deck,
so `dGM` is `skip`, never a fabricated pass) -- confirms the adjoint
machinery is correct on a REAL nonlinear small-signal matrix (BSIM/MOS1
`gm`/`gds`/`Cgs`/`Cgd` stamps), not just on linear R/C/L/VCCS elements.
`M2_w` (the op-point-DEPENDENT case) is `EXPECT-DIVERGE`, **~40-100%
relative "error"** -- exactly the missing-op-point-term limitation
documented above, checked to actually be present and large, not silently
inside tolerance.

**`DELTA` choice**: central differences use `DELTA = 1e-3` (0.1%
relative parameter step). Checked this pass: FD gradients at
`DELTA=1e-3` and `DELTA=1e-4` agree to **~1e-6 relative** (nowhere near
floating-point-roundoff-limited yet); `DELTA=1e-2` already differs from
`DELTA=1e-3` by **~1e-4 relative** (early curvature) -- `1e-3` sits
comfortably in the stable middle, confirmed empirically, not merely
assumed.

**Margin-sensitivity sweep density**: `dT/dp` itself is accurate to
floating-point noise (~1e-6 relative) at ANY reasonable sweep density,
since it is read off directly per swept frequency point -- no
interpolation involved. `dPM/dp`/`dfc/dp`/`dGM/dp`, by contrast, need a
bracket finite difference of the swept trace for the implicit-
differentiation frequency-derivative term, so THEIR accuracy scales with
sweep density near the crossing -- measured this pass on R1b's `R1o`:
`dec 25` gives ~1e-3 relative error, `dec 100` ~1e-4, `dec 400` ~1e-5
(and going far past that, e.g. `dec 1000` over a narrow range, can
actually get WORSE again from floating-point cancellation in the
bracket's now-tiny frequency step -- also observed this pass on R3avg).
`tests/test_sens.py` uses `dec 200`-`dec 300` over roughly one decade
bracketing each fixture's known crossing (found first with a coarse, wide
scan), which was dense enough to clear the 1e-4 bound with margin
everywhere it was tested.

### Bugs found and fixed validating this (all four are in
`patches-wip/sens/0400`-`0403`, see that directory's `README.md` for the
full commit messages)

1. **Early `break` out of the per-parameter `sgen` walk once `k>=nsens`**
   corrupted `ngspice`'s own temporarily-unlinked instance chains
   (`sgen_next()`'s isolation trick needs the walk to run to completion to
   restore them) -- silently dropped a device from every SUBSEQUENT
   frequency's real load. Found because `dT/dp` for the FIRST parameters
   requested went exactly to zero from the second swept frequency onward
   while the LAST one stayed correct. Fixed: `continue`, never `break`
   (the same pattern `cktsens.c`'s own per-frequency `sgen` walk already
   uses for its filtered-out parameters).
2. **`DEVsetup` on the throwaway per-parameter matrix silently repoints a
   device's own PERSISTENT cached matrix-element pointers** (some device
   types, e.g. capacitors, cache pointers directly on the instance struct
   at `DEVsetup` time and never call `SMPmakeElt` again) at the throwaway
   matrix, not the real one -- restoring `ckt->CKTmatrix` back to the real
   `Y` afterward does NOT undo this. Found because `T` ITSELF (not just
   `dT/dp`) went wrong from the second swept frequency onward once 2+
   `sens=` parameters with this property were requested (a lone resistor
   alone did not show it). `cktsens.c`'s own DC/AC `.sens` never hits this
   because it calls `CKTunsetup()+CKTsetup()` fresh before every frequency
   point, which happens to re-fix every pointer as a side effect;
   `STBan()`'s per-frequency loop does not do that (by design -- one
   factor, two solves, no full circuit re-setup). Fixed: re-run each
   device's `DEVsetup` against the REAL matrix immediately after its
   sensitivity work, before moving to the next parameter.
3. **`dGM/dp` was missing the gain-crossing frequency's own implicit
   shift with `p`** (only the "evaluated at a fixed f" partial term was
   computed) -- wrong by O(1) (relative errors of 0.5-1.2 against central
   FD, far above `dPM/dp`'s plain bracket-discretisation-level error).
   Fixed: the same two-term implicit-differentiation structure `dfc/dp`
   already had, applied to the `-180deg` phase constraint instead of the
   `|T|=1` one.

### Red-proofs

Both built as throwaway one-line-diff variants of the SAME commit
(`aaba50ba3f49-sens-redproof-yt`, `aaba50ba3f49-sens-redproof-sign` --
not kept; rebuildable from `patches-wip/sens/0400`-`0403` plus the one-line
diff each name describes), run through the SAME `tests/test_sens.py`
comparison machinery via `NGSPICE_REDPROOF_YT`/`NGSPICE_REDPROOF_SIGN`:

- **`Y` instead of `Y^T`** for the two adjoint solves (`SMPcSolve` in
  place of `SMPcaSolve`): `dT/dR1o` on R1b comes out with **100% relative
  error** (correctly fails).
- **Dropped `dY/dp` stamp sign** (`+` instead of `-` in the
  `db/dp - (dY/dp) x` contraction): **200% relative error** (correctly
  fails).

Neither red-proof is silently inside the 1e-4 bound -- both are exactly
the kind of "plausible-looking but structurally wrong" bug the method is
vulnerable to (get the transpose or the sign backwards and you still get
SOME complex number out, not a crash), and both are caught by a factor of
1000-2000x over the tolerance.

## Speed: adjoint vs N finite-difference reruns

Not separately benchmarked with wall-clock numbers this pass (see "what
was not done" below) -- the DESIGN's own claimed complexity is
architectural, not measured: **2 extra transposed solves per frequency,
total, regardless of the parameter count**, versus finite differences'
**2 full extra `.stb` reruns (a fresh operating point, N frequency
points, a factorisation each) PER PARAMETER**. For R1b/R3avg's 3-parameter
validation runs, this is already a re-run-count ratio of 6-to-1 in the
adjoint's favour (2 adjoint solves total vs 6 full `.stb` reruns for
central differences on 3 parameters) that only grows with more
parameters; for `sens=all` on a circuit with dozens of device parameters
the ratio would be dozens-to-one. The per-parameter marginal cost inside
one adjoint run (one sparse matrix build + one matvec + a few dot
products, all `O(nnz)`, no solve) is architecturally cheap relative to a
whole extra circuit solve, but this pass did not instrument or report an
actual wall-clock number for either side -- stated as a real gap, not
implied to be measured.

## The bias-point term and nonlinear devices (card #595, 2026-10-11)

Card #475's derivative of the small-signal matrix was the difference of the
device's own AC stamp at p and p+delta, at fixed node voltages. Two terms
were missing:

1. **The operating-point shift.** `dY/dp|total = dY/dp|x + (dY/dx)(dx_op/dp)`.
   `dx_op/dp` now comes from one DC Newton step from the converged point
   with the parameter at p+delta and at p-delta (central difference; two
   real factorisations per parameter, once, before the sweep). At the
   converged point the derivative of the Newton map is `-J^-1 dF/dp`, which
   is `dx_op/dp`.
2. **Quantities a nonlinear device computes in its DC load.** A MOSFET's gm
   and gds, and every junction or gate charge, are computed in `DEVload`
   and only read by `DEVacLoad`. Re-running the AC stamp at p+delta
   therefore did not see their dependence on p, even at a fixed operating
   point (gm on `w`, a junction capacitance on `cbd`). The old `M2_w` value
   was only the overlap-capacitance part.

For every parameter of a device that is not a linear element (R, C, L, K,
E, F, G, H, V, I), the whole derivative is now taken through the full
linearisation: the circuit is re-linearised at (p + h, x_op + h dx_op/dp)
and at (p - h, x_op - h dx_op/dp) (`CKTload` in MODEINITFLOAT, then
MODEINITSMSIG, then `CKTacLoad`), each matrix is multiplied by the two
forward solutions, and the difference replaces the stamp difference. h
makes the largest node-voltage shift 1e-5 V (`STB_SENS_BIAS_DV`), at most
1e-4 of the parameter; a parameter that does not move the bias uses 1e-5 of
its value. Linear elements keep the exact stamp difference and add the bias
term when `dx_op/dp` is not zero (a bias resistor, a source's DC value).
After every shift the device states are restored from a copy and the
operating-point linearisation is re-run, so T itself is bit-identical with
and without `sens=` (checked).

Two older defects were found on the way and fixed in the same patch:

- **A source's AC magnitude had a non-zero sensitivity** (`Ib_acmag`
  2.8e3 on R1d). `.stb` drops every AC stimulus, so `db/dp` is zero; the
  device's stimulus stamp is no longer used.
- **`spMultiply` writes only the rows present in the scratch matrix's
  structure.** The rest of the output vector kept whatever was there. Card
  #475's code was safe only because nothing else wrote those rows; it is
  now cleared before every product.

**Validation** (`tests/test_sens.py`, 46 pass, 0 fail, 7 skip, 2 info of 55
checks; against full `.stb` reruns with the parameter at +-0.1 %):

| Deck | Parameter | What it exercises | dT/dp worst rel. error |
|---|---|---|---|
| R1d | `M2_w` | instance parameter, moves the bias (was EXPECT-DIVERGE, 103 %) | 2.8e-7 |
| R1dbias | `Ib` | linear source, bias term only | 1.5e-7 |
| R1dbias | `M2_w` | as above | 2.7e-7 |
| R1dbias | `M2:vto` | model parameter (own model), moves the bias | 4.9e-6 |
| R1dbias | `M2:cbd` | junction capacitance: charge only, bias unchanged | 3.3e-9 |

dPM/dp and dfc/dp of the same parameters agree to 1e-7 to 5.2e-5 (the bracket derivative of the sweep sets that, see "Margin-sensitivity sweep density"). Red proof:
`STB_SENS_FIXED_OP=1` (card #475 behaviour) gives 1.03 relative error on
`dT/dM2_w`. Also checked: T bit-identical with and without `sens=`;
`Vin_acmag`, `Ib_acmag` exactly 0; the warning below appears.

**Model parameters of a model with several instances keep card #475's
result, which is incomplete**, and one warning lists them:
`sens=M1:vto` is per instance (the parameter walker isolates one instance),
but the full linearisation needs every instance loaded, and setting a model
parameter then moves all of them. Give the device a model of its own (as
`tests/fixtures/sens_r1d_bias.cir` does for M2) for the full derivative.

**Cost**: `sens=all` on R1dbias (150 parameters, 301 frequencies): 0.25 s
against 0.083 s with the operating point fixed and 0.005 s without `sens=`
(development workstation: AMD Ryzen AI Max+ 395, Linux; 2026-10-11). Each nonlinear-device parameter costs four `CKTload` and
two `CKTacLoad` per frequency.

Not done (card #595): no BSIM3/BSIM4 deck validated (only level-1 MOSFETs);
no shared-model full derivative; the absolute step for a parameter whose
value is 0 is still 1e-6 in its own unit (for a length such as `ld` that is
1 um, which makes the effective channel length negative: ~1,800 warnings
under `sens=all` on R1d, already there on s131); instant mode keeps the
fixed operating point (warned); KLU still refused; the periodic half of
#475 (PSS) not started.

## What was not done

- **No wall-clock speed benchmark** (adjoint vs N-reruns), as above --
  design/architecture only.
- **(Done by card #595, above.)** **DC operating-point sensitivity dx_op/dp is not included** (see "What
  is / is not included") -- the single largest scope boundary of this
  pass, checked and documented (the `M2_w` `EXPECT-DIVERGE` case), not
  merely disclaimed.
- **Model-parameter sensitivity beyond what `sens=all`'s generic walk
  already reaches was not specially validated** -- `M2_w` (MOSFET
  channel width) is the only non-principal/model-style parameter
  exercised, and it is the op-point-DEPENDENT (expected-divergence) case,
  not a clean pass/FAIL data point for THAT class of parameter on its
  own merits.
- **Device-mode probes, differential/common-mode probes, and KLU are all
  refused, not implemented** for `sens=` -- a scoping decision (see
  "What is / is not included"), not an oversight found late.
- **Pause/resume mid-sweep is refused for `sens=` runs** -- the
  sensitivity accumulator arrays are plain locals in `STBan()`, not
  designed to survive the function returning and being re-entered.
- **No independent second-reader verification** -- this pass, its bug
  fixes, and its validation numbers are all from one session.
- **The periodic/PSS half of card #475 (adjoint sensitivity on a
  shooting-PSS orbit, building on `pacshoot.c`'s recording pass) is
  explicitly out of scope for this pass** -- per the card's own split
  into a DC-op/AC-small-signal first half (this document) and a periodic
  second half, later.
  (Pointer, 2026-10-14: the periodic half exists as finite differences, not as an adjoint -- `tools/pstbsens.py`, `docs/pstb-sens.md`.)
- **Sweep density for the margin-sensitivity bracket derivative is a
  per-deck tuning knob, not automatic** -- this pass found (empirically,
  see "Margin-sensitivity sweep density" above) that too coarse a sweep
  under-resolves `dPM/dp`/`dfc/dp`/`dGM/dp`, and too fine a sweep can
  reintroduce floating-point-cancellation noise; native `.stb` does not
  warn about either condition itself, and this pass did not add such a
  warning (out of scope for a first pass; a real usability gap for
  anyone running `sens=` without also having validated their own sweep
  density against a spot-check the way `tests/test_sens.py` did here).
