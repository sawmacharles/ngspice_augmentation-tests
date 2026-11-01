# Bug 1 (card #539, 2026-10-01): `solver=td` .pac/.pstb/.pxf/.pnoise drop the history of device-model charge storage -- silently wrong on every transistor-level converter

Engine: `aaba50ba3f49-s88` (and by source reading every build with `pacshoot.c` as in `ngspice-src-537`).
Severity: silently wrong results (no warning, no refusal). `.pss method=shoot` itself is NOT affected (orbit checked, see
the R7 section of `refs/README.md`).

## Symptom on R7 (BSIM4 buck, `refs/r7_bsim_buck.cir`)

| quantity | `.pac`/`.pstb` solver=td (500 Hz .. 200 kHz) | converged FRA reference |
|---|---|---|
| control-to-output (Vctrl at the comparator input), 20 kHz | 4e-7 V/V (|H| ~ -128 dB) | 1.502 V/V (+3.6 dB, -7.5 deg) |
| line-to-output, 500 Hz | 0.304 V/V (= the open-loop duty, i.e. the loop is absent) | see `fra/curves-r7-*/` (closed-loop, tens of dB lower) |
| `.pstb` loop gain |T| | 2e-10 (no 0 dB crossing, "gain margin 201.8 dB") | crossover in the 10-50 kHz class |

Walking the sideband-0 response down the chain comparator -> inverters -> NOR/delay -> driver (`pac_c2o_n.cir`,
1 kHz, V/V per V of control): `cout` 1.14, `cpb` 0.61, `pwm` 3.1e-2, `pwm_b` 1.7e-4, `hs_on` 2.1e-6, `ghs` 1.3e-8, `sw` 2.4e-7.
Each BSIM4 stage destroys most of the timing information the previous one carried.

## Minimal reproducer: `run_repro.py` (writes only to `./work`)

A chain of N inverters driven by `ramp + Vc` (2 us). Low-frequency control-to-average-output gain = static d<v>/dVc of the
periodic steady state (two `.tran` runs, central difference). Correct answer: |gain| ~ 0.50-0.56 for every N (the chain only
delays an edge). `.pac solver=td` at 1 kHz (`pacpoints=400 / 40000`):

```
models         stages     static     pac(400)   pac(40000)
level1_cg0     1          0.5360       0.5325       0.5355      <- MOS1, no capacitance: correct (grid-converging)
level1_cg0     2          0.4865       0.4677       0.4981
level1_cg0     3          0.4320       0.3300       0.5228
level1_cg1e-8  1          0.5355       0.4991       0.2507      <- MOS1 + overlap caps cgso=cgdo=1e-8 F/m: WRONG, worse on a finer grid
level1_cg1e-8  2          0.5020       0.1102       0.0600
level1_cg1e-8  3          0.5250       0.0035       0.0138
bsim4          1          0.5650       0.5546       0.0214      <- BSIM4 (the d_analog_block card): same
bsim4          2          0.5020       0.3125       0.0201
bsim4          3          0.4990       0.0406       0.0028
```

(`run_repro.out`.) Not BSIM4-specific: any MOS (or diode with Cj/tt, BJT, ...) whose capacitance lives in the model.

## What it is NOT (each ruled out by experiment, files in `../inv/`)

- not step-size/tolerance: `reltol` 1e-3..1e-6 and `trtol` change the numbers by <= factor 2 and never fix it;
- not the stale Jacobian: `PAC_EXACT_JACOBIAN=1` gives identical numbers to 3 digits at pacpoints 400/4000/40000;
- not explicit capacitors: adding 100 fF / 1 pF `C` elements on every stage output (these ARE reduced states) leaves
  the decay unchanged (`cx_*`: 0.554, 0.30, 0.04-0.07, 0.001-0.003 for N=1,2,3,4);
- not the extraction/demodulation: level-1 without capacitance converges to the right value as `pacpoints` grows.

## Diagnosis (source reading, `src/spicelib/analysis/pacshoot.c` header "DESIGN")

The reduced state is "one entry per native capacitor's terminal voltage, one per native inductor's branch current"
(same map as `pssshoot.c`'s `PSSstateMap`; only `cap/inddefs`, `vsrc`, `isrc` are included). The per-step recursion is
`J_k y_k = (companion-history term built from y_{k-1}, same recursion capload.c/indload.c use) + b e^{jwt_k}` with
`J_k = G_k + ag0_k C_k` taken from the real transient's last factorisation. `J_k` therefore contains the device-model
`ag0*C` terms (BSIM4 Cgs/Cgd/Cdb/Csb/Cgb, MOS1 overlap, diode Cj) but `y_{k-1}`'s contribution through those
charges (the trapezoid `i_k = 2/h (q_k - q_{k-1}) - i_{k-1}`) is not in the history term: every model capacitance acts as a
conductance `ag0 C` to a perturbation-free reference. A conductance `2C/h` that grows as the step shrinks pins the node:
the perturbation is damped more the finer the grid (the 40000-point column above), and a signal that has to traverse
several nodes (a gate chain) is damped per stage (0.55 -> 0.31 -> 0.04 -> 0.001). The header's "WHAT THIS DOES NOT DO"
list does not mention it and nothing at run time checks for it.

## Suggested fixes (not done, card scope)

1. Cheap and honest first: at `.pac/.pstb/.pxf/.pnoise ... solver=td` start, scan the device list for models with
   charge storage (any device whose `DEVacLoad` capacitance is non-zero: MOS/BSIM, diode with CJO/TT, BJT, JFET) and
   either REFUSE with a clear message or print `PAC-TD: WARNING device-model charges are not part of the reduced state`.
2. Proper: add the device charges to the state (qd, qg, qb... from `CKTstate0`, restored per step as `PSSsensRecordStep`
   does for C/L) and propagate their history `-i_{k-1}`/`q_{k-1}` perturbation in the recursion; the monodromy
   (`Phi`) then has the same states.
3. Regression: this reproducer as a test (static derivative vs `.pac` at 1 kHz, N=1..3, MOS1 with cgso/cgdo).
