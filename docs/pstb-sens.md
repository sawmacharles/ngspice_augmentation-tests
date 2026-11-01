# Parameter sensitivities of the periodic loop gain (`tools/pstbsens.py`)

Card #475, periodic half, 2026-10-14.
Tool `tools/pstbsens.py`, test `tests/test_pstbsens.py`, fixture `tests/fixtures/pstbsens_r1d_clk.cir`.
No C changes: runs the fork build `~/Work/apps/ngspice_augmentation/aaba50ba3f49-s136/bin/ngspice` (env `NGSPICE` overrides).

## What it computes

For a deck with a `.pstb` card (switching converter, sideband-0 Tian loop gain `T(f)` on the PSS orbit) and a
parameter `p`:

- `dT/dp` at every sweep frequency (complex; `T` is the `t` vector `.pstb` writes, same sign convention);
- `dPM/dp`, `dfc/dp` (phase margin and 0 dB frequency), `dGM/dp` when a -180 deg crossing exists;
- normalised sensitivities `p/PM * dPM/dp` (and `p/fc * dfc/dp` in `summary.json`, key `nfc`);
- with `--orbit NODE`: `d mean/dp` and `d peak-to-peak/dp` of `v(NODE)` over one period of the PSS orbit.
  The orbit is reachable in the control block: after `run` the PSS plot is the only plot holding `v(node)`, so
  `foreach p $plots / setplot $p / wrdata orb.txt v(node) / end` writes it (one period, `time` as scale; the mean uses
  the periodic trapezoid with `v(0) = v(T)`).

Margin definition (identical to `pacshoot.c`, `pac_td_pstb`): phase unwrapped from the first point, `|T|` in dB;
0 dB and -180 deg crossings by linear interpolation in `ln f`; `PM = 180 + phase` at the 0 dB crossing; the lowest of several
crossings is reported; `GM = -|T|_dB` at the lowest -180 deg crossing. The tool recomputes PM and fc from `T` in every
run and refuses if they differ from the `pstb: phase margin X deg at F Hz` line by more than the line's 4 significant digits
(1e-3 relative PM, 1e-5 relative fc). The line itself is not used (4 digits would be a 1e-4 quantisation of a difference).

## Method

Central finite differences over the complete `.pss` + `.pstb` run, one deck per evaluation:

```
D(s) = (f(p+s) - f(p-s)) / (2 s)        s = h p   (relative; default h = 1e-2)   or s = --abs-step when p = 0
D    = (4 D(h/2) - D(h)) / 3            Richardson: cancels the s^2 term
e    = |D(h/2) - D(h)|                  error estimate, per output (per frequency for dT/dp)
```

4 runs per parameter plus one nominal. Runs go through `systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0` with
`OMP_NUM_THREADS=2`, `--jobs N` at a time (default 4, at most 5). The results do not depend on `--jobs` (checked: byte-identical).

**Why finite differences and not an adjoint.** The `.pstb` result rests on the time-domain recording pass with all its
switch-event machinery (ideal-switch jumps, fast nodes, routed jumps, hand-over, charge history; patches 0120-0133).
A finite difference of the whole run inherits every one of those paths with no new code, so it cannot disagree with
what `.pstb` itself reports; an adjoint of the periodic map would need each event path differentiated again. The adjoint
exists for the fixed-bias-point half (`.stb sens=`, `docs/sens.md`), and the two are checked against each other below.
The price is cost (4 full runs per parameter) and noise (below).

## Syntax

```
python3 tools/pstbsens.py DECK --param P [--param Q ...] [--h 1e-2] [--abs-step A] [--jobs 4]
        [--orbit NODE] [--out DIR] [--pss-reltol X] [--margin-interp chord|cubic] [--timeout S] [--no-systemd]
```

`--param` names (mirrors `sens=` in `docs/sens.md`):

| form | meaning | how it is perturbed |
|---|---|---|
| `cl_val` | a `.param` name | its value on the `.param` line (an expression value `{expr}` becomes `{(expr)+(delta)}`) |
| `CL`, `R1`, `Vdd` | element name: principal value (R, C, L; V/I `DC` value) | the value token on the element line |
| `M2:w`, `M2_w` | instance parameter | the `w=...` token on the element line |
| `M2:vto` | parameter of the instance's model | the `vto=...` token on the `.model` line; a model shared by several instances moves all of them (a note is printed, as `sens=` warns) |

The deck is rewritten per run (text only; continuation lines joined, relative `.include` made absolute). The deck's own
`.pstb` card is kept; **every `.control` block of the deck is replaced** by `run; wrdata t` (and the orbit write), and a note
says so. Outputs in `--out` (default `~/Work/data/ngspice_augmentation/pstbsens/<deck>-<time>/`): `dtdp.csv`
(`freq, T, and per parameter dT/dp re, im, err`), `summary.json` (all scalars, with errors), `timing.json`, `runs/<run>/` (decks and logs).
Table on stdout: param, value, dPM/dp, +-err, normalised, dfc/dp, +-err, dGM/dp. Exit 0 = ran, 2 = refused.

Refused, with a message: no `.pstb` card; unknown parameter; a parameter whose `+h` run gives a bit-identical `T` ("had no effect": not found, or not used by
the loop/orbit); PSS not converged in any run (no `PSS-SHOOT: CONVERGED`); PM line missing (no 0 dB crossing in the band); value 0 without `--abs-step`; `--jobs` outside 1..5.

### Two knobs that matter

- `--pss-reltol 1e-8` (sets `PSS_SHOOT_RELTOL` for every run; engine default 1e-5). On a switching deck the engine's `T` is not a smooth function of
  `p`: the adaptive time grid of the PSS changes with `p` and `T` jumps by about 1e-3 relative, PM by 0.02 deg, fc by 25 Hz, for a 1e-9 relative change of a parameter
  (measured on `rd_pstb_r3d.cir`, `cval`). A 1e-2 step produces a PM change of 0.2 deg: the error estimate came out as large as D or larger (e up to 40 x D on `cval`). With 1e-8 the noise is 5e-6 deg (rms residual of a cubic through
  9 points of `cval`) and e is 0.003-0.7 %. It is not the default because it makes the linear fixture fail in the PSS (timestep too small); use it for switching decks.
  **Since s136 (card #625) it is usually not needed.** The noise came from the adaptive grid of the `.pstb` recording passes, not from the PSS. Those passes now run at
  reltol 1e-8 on decks without SW/CSW/XSPICE. On `rd_pstb_r3d.cir` (h = 0.03, no `--pss-reltol`, jobs 2), e/|D| for (dPM, dfc) is: `cval` 0.034 %, 0.063 %;
  `esrval` 0.016 %, 0.39 %; `r2` 0.005 %, 0.06 %. That is equal to or better than with `--pss-reltol 1e-8`, at 6.8 s for 13 runs instead of 18.4 s
  (`~/Work/data/ngspice_augmentation/625/pstbsens_s136/`). Keep the option for decks with ideal switches or code models, which keep the inherited recording tolerance.
- `--margin-interp`. Default `chord` is the `.pstb` definition. The chord PM has a kink every time the crossing moves into the next sweep cell, so
  a step that spans cells gives a wrong `dPM/dp` when `dPM/dp` is small (OTA fixture, `dPM/dIb`: 0.5 % off at h = 0.05 with 200 points/decade; at h = 0.003 it agrees with the adjoint to 5e-5). Keep the step small
  enough that fc moves less than one cell, or use `cubic` (cubic through the four nearest points in `ln f`; smooth in `p`, converges to a value 0.5 % different from the chord adjoint on that case, i.e. the chord's own interpolation error).
  Both are always computed; `summary.json` has `dPM_alt`, `dfc_alt` for the other one.

## Worked example (`tests/fixtures/rd_pstb_r3d.cir`, voltage-mode buck with dead time)

```
PSS_SHOOT_RELTOL set by --pss-reltol 1e-8
python3 tools/pstbsens.py tests/fixtures/rd_pstb_r3d.cir --param cval --param esrval --param r2 --h 0.03 --pss-reltol 1e-8 --orbit vout
param               value      dPM/dp      +-err  p/PM dPM/dp      dfc/dp      +-err      dGM/dp
cval           1.0000e-04  2.4219e+05 9.5009e+01   6.2709e-01 -2.0293e+08 1.4572e+05         n/a
esrval         1.0000e-02  1.1386e+03 3.1410e-01   2.9480e-01  6.9829e+04 4.7913e+02         n/a
r2             2.7510e+04 -7.3561e-04 1.9212e-08  -5.2398e-01  7.1622e-01 3.1422e-04         n/a
(PM = 38.6215 deg at fc = 33447 Hz, nominal; central difference + Richardson)

param            d mean/dp      +-err    d p-p/dp      +-err   (v(vout) over one PSS period)
cval           -6.1645e-04 1.0187e-04 -3.9675e-02 2.7080e-05
esrval          4.8576e-06 5.2830e-07  4.7942e-01 2.0000e-06
r2              2.5231e-12 2.1175e-12  9.4656e-13 3.0272e-13
13 runs, 18.4 s wall, jobs=4 (mean run 4.6 s)
```

Reading: `p/PM dPM/dp = +0.63` for `cval`: a 1 % larger `cval` raises PM by 0.63 % (0.24 deg) and lowers fc by 0.61 %;
a 1 % larger compensator `R2` lowers PM by 0.52 %. ESR raises PM (zero) with a small effect on fc. `r2` does not move the regulated orbit (the d mean, d p-p of `r2` are noise at 1e-12);
the output ripple depends on ESR (0.48 V/ohm) and, weakly, on `cval`. `h = 0.03`, not the default 1e-2: at 1e-2 the PSS does not converge in the `cval` and `r2` runs (refused, see limits).

## Validation (`python3 tests/test_pstbsens.py`: 32 pass, 0 FAIL, 6 info of 38 checks, 2026-10-14)

1. **Equivalence with the `.stb sens=` adjoint** (fixture `pstbsens_r1d_clk.cir` = `sens_r1d_bias.cir` + a clock decoupled from the loop; `.pstb 1meg 0 clkgnd 1000 1 5 0.001 dec 200 5e6 1e8 ... pacpoints=200000`; h = 0.003):
   `.pstb` T equals `.stb` T to 8.4e-7 (worst, at 99.8 MHz; `pstb_linear` gate 1e-6). Worst relative difference of dT/dp over the 261 frequencies (gate 1e-4):

   | parameter | `sens=` | worst dT/dp diff | dPM/dp diff (gate 1e-3) | dfc/dp diff |
   |---|---|---|---|---|
   | `ib_val` | `Ib` | 7.9e-7 | 5.4e-5 | 2.2e-7 |
   | `m2w` | `M2_w` | 8.9e-7 | 3.8e-6 | 2.8e-7 |
   | `M2:vto` | `M2:vto` | 6.0e-6 | 2.0e-5 | 5.2e-6 |
   | `CL` | `CL` | 8.4e-7 | 7.4e-6 | 2.8e-7 |

   The agreement includes the bias-point term of patch 0133. The tool's own e/|D| was at most 6.5e-6.
2. **Switching orbit** (`rd_pstb_r3d.cir`, h = 0.03, `--pss-reltol 1e-8`): e/|D| (dPM, dfc): `cval` 0.039 %, 0.072 %; `esrval` 0.028 %, 0.686 %; `r2` 0.003 %, 0.044 %.
   Independent check, two plain ngspice runs at p(1 +- 0.005) with the test's own T parsing and margin code: difference to the tool 0.001 % (cval dPM), 0.003 % (cval dfc), 0.012 %, 0.286 % (esrval), 0.007 %, 0.012 % (r2); all inside max(3e, 2 %|D|).
3. `--jobs 1` and `--jobs 4`: `dtdp.csv` and `summary.json` byte-identical (`esrval`, `r3`, with `--orbit vout`).
4. Refusals: unknown parameter, parameter without effect, no `.pstb`, PM line missing, PSS not converged (cval at h = 1e-2), `--jobs 6`.
5. **Red proof**: a one-sided forward difference at h = 0.1 without Richardson (hidden `--redproof-forward`) misses check 1's dT/dp gate by a factor 893 (worst 8.9e-2 on `CL`; 2.5e-2 `ib_val`, 4.9e-2 `m2w`, 3.5e-3 `M2:vto`), against 6e-6 for the real method.

## Cost (measured)

| deck | run | one parameter, jobs = 1 | one parameter, amortised at jobs = 4 |
|---|---|---|---|
| `rd_pstb_r3d.cir` (9 states, 14 frequencies, `--pss-reltol 1e-8`) | 4.6 s | about 18 s | about 5 s (13 runs in 18.4 s for 3 parameters) |
| `pstbsens_r1d_clk.cir` (`pacpoints=200000`, 261 frequencies) | 37 s | about 150 s | about 45 s (17 runs in 181 s for 4 parameters) |

The cost of a run is the cost of `.pss + .pstb`; the sweep length enters through the recording pass (0.06 s/point at `pacpoints=200000`).

## Limits

- **Noise on switching decks.** See `--pss-reltol` above. Even with 1e-8 a PSS run at a given `p` may fail to converge (observed: `cval` and `r2` at h = 1e-2, `esrval` at h = 0.2, `rload` at h = 0.2): the tool refuses; change `h` (`0.03` worked for all three here). A parameter whose effect on PM or fc is
  near the noise (rload on fc: 22 vs 43 Hz/ohm at h = 0.02 / 0.01) has a wide `e`; read `e`, not the digits.
- **Step choice.** `h` must be large enough that the change of PM/fc clears the noise and small enough that Richardson's remainder (h^4) and the chord kink stay small (see `--margin-interp`). `e` shows which side you are on.
- **PM kinks**: if the number of 0 dB crossings changes between the runs a warning is printed (margin derivatives are not smooth there).
- Only the sideband-0 `T` that `.pstb` produces; the limits of `.pstb` itself (smooth-switch/ideal-switch scope, probe a plain top-level `DC 0` V-source, no `lgprobe`, free-running orbits) apply (docs/pss-pac.md). The probe trap applies to the fixture: `Vprobe` is a plain source.
- Perturbable values: `.param`, R/C/L value, V/I `DC` value, `name=value` tokens of element and `.model` lines. Not `.subckt`-internal `.param` overrides, not `.options`, not `m=`-style multipliers inside parentheses, not parameters given only by position (MOSFET W/L positional), not `.include`d files' contents (the text of the main deck only).
- Nominal value of a parameter comes from a small expression evaluator (numbers with SPICE suffixes, `.param` names, `+ - * / ^`, `sqrt exp log log10 sin cos tan pow min max abs`); other functions are refused.

## Not done

- `dGM/dp` is implemented and printed, but **no test exercises it**: neither fixture has a -180 deg crossing (the table shows `n/a`).
- `--abs-step` (parameter equal to 0) is implemented, **not tested**.
- `--orbit` is tested only for being present in the output; its numbers were **not checked against an independent calculation** (no plain-run cross-check of mean / peak-to-peak).
- Only `PSS_SHOOT_RELTOL` was tried as a remedy for the noise; `Points`, `pacpoints`, the PSS tolerance field and `.options reltol` were compared on `cval` only (documented in the card notes, not kept as tests). Diagnosed in card #625: the recording-pass grid, not the PSS orbit (see above and docs/pss-pac.md).
- No periodic adjoint (the second half of the original card #475 idea); no sensitivity of the closed-loop transfer functions (`.pac`) or of `.pnoise`; no multi-parameter simultaneous step (cross terms); no Hessian.
- Ran on the development workstation (AMD Ryzen AI Max+ 395, Linux) with s133 only; other builds untested. One session wrote the tool, the test and this page; no independent second reader.
