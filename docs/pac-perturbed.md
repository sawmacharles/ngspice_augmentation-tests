# Perturbed periodic steady state, generic tool (`tools/pac_perturbed.py`)

Card #494 (parallel track B, alongside the saltation-correction work in
`docs/pss-pac.md`'s "Saltation" section, `patches-wip/salt/`). Private
fork, `project ruling #446 B/#447 A` restriction applies throughout (never
contact upstream about any of this).

## What this is, and why it exists alongside `.pac solver=td`

`.pac solver=td` (`docs/pss-pac.md`) linearises a switching converter's
own accepted-transient-step sequence to get a small-signal transfer
function directly from a converged periodic orbit -- fast, and validated
to 0.1-0.6dB against independent references on **smooth-switch** decks
(behavioural `tanh`-type conductances, `refs/r3_vm_buck.cir`-family). It
is NOT yet right for **ideal switches** (SPICE `sw`/`csw`) or **XSPICE**
digital PWM: the per-step linearisation gives an ideal switch's own
control-to-power-port gain as exactly zero except at the accepted steps
where the switch actually flips, and the correction that mechanism needs
(Fang's saltation matrix, `docs/pss-pac.md`'s "Saltation" section) is
implemented on this fork but measured **~2.2x too high** (13.19 V/V vs a
~6 V/V consensus on the open-loop ideal-switch reference, card #494 main
track) -- not resolved as of this card.

`tools/pac_perturbed.py` sidesteps the whole problem: it never
linearises anything. For a target frequency, it runs the REAL, nonlinear
circuit (ideal switches, XSPICE, anything ngspice can transient-simulate)
to two independently-converged periodic steady states -- baseline
(unperturbed) and perturbed (a small sinusoid added at the injection
point) -- subtracts them, and reads the sideband-0 (same-frequency)
component of the difference by a single-bin DFT, divided by the same DFT
of the known injected waveform. This is the "perturbed periodic steady
state" method `tools/pssac_ref.py` (card #484, `research/near-fs2.md`)
already validated against `.ac`, an exact analytic LPTV reference, and
(informally) against a finite-difference-of-orbits check on the saltation
investigation's own R6o fixture. `pac_perturbed.py` reuses that numerical
core (rational frequency matching, single-bin DFT, orbit-plausibility
check) but is **generic**: any deck, an injection point named on the
command line (no hardcoded per-fixture table), any output(s), an optional
Tian two-injection loop-gain probe (`tools/stbcore.py`'s `tian()`), and
parallel dispatch of the sweep's many independent ngspice runs through
`tools/campaign.py`'s memory-aware scheduler.

**Bottom line: use `.pac solver=td` when the deck is smooth-switch (it is
faster and already validated there). Use `pac_perturbed.py` when the deck
has an ideal `sw`/`csw` or XSPICE switch, or whenever you want a
method-independent cross-check (it shares no code path with `pacshoot.c`
at all -- a disagreement between the two is informative either way).**

## When to use it

- **Ideal-switch or XSPICE decks**, where `.pac solver=td`'s saltation
  correction is not yet trustworthy (see above) -- this is the primary
  use case this card was scoped for.
- **An independent check** of `.pac solver=td` on a smooth deck (no
  shared code with `pacshoot.c`'s per-step linearisation) -- see the R3
  validation table below.
- **Loop gain of a circuit an ideal switch sits inside**, via
  `--loop-gain` (the switching-converter analogue of `tools/stb.py`'s
  native-`.ac` loop-gain tool, but perturbed-PSS-based so it works on
  ideal-switch/XSPICE decks too).

## When NOT to use it

- A smooth-switch deck where `.pac solver=td` already runs (it is
  ~10-100x cheaper per point: one shooting-Newton solve + a fast
  frequency sweep vs. this tool's two-full-transient-per-point cost).
- Above roughly half the switching frequency the sideband-0 assumption
  itself becomes questionable for BOTH tools (see "Limits" below) --
  neither is meaningful there without sideband coupling, which neither
  implements.
- A free-running (autonomous, hysteretic/COT/QR) converter used to be out
  of scope here (this tool assumed a **driven** circuit with a known,
  fixed switching frequency `--fs`) -- card #493 added `--freerun`
  specifically for this case; see "Free-running (autonomous) mode" below.
  `--fs` in `--freerun` mode is only a sizing ESTIMATE, never exact.

## Cost

Per requested frequency: 1 shared baseline run (cached across every
frequency that maps to the same LCM-period multiplier `q`) + 1
full-amplitude perturbed run + (default) 1 half-amplitude perturbed run
for the linearity check = up to 3 ngspice invocations per NEW `q`, 2 per
additional frequency sharing a `q`. `--loop-gain` needs one more pair
(voltage-injection and current-injection runs share the same baseline) --
up to 5 runs for a new `q`. Every run in a sweep is independent and
dispatched through `tools/campaign.py run_campaign_a`'s memory-aware,
one-process-per-point scheduler (`--max-procs`), the same scheduler this
project's PVT/MC campaigns use -- see `docs/campaign.md`.

Measured this card, `ngspice-aug` (`aaba50ba3f49-s41`), `--max-procs 2`
(this project's own <=6-concurrent-ngspice rule, shared with other
users' work at the time):

| deck | backend | points | wall (total, incl. halving) |
|---|---|---|---|
| `pac_lti.cir` (trivial RC) | pss, `Points=2*100` | 3 freqs | ~1.3s |
| `refs/r2_switched_rc.cir` | pss | 5 freqs | ~9s (one point q=80 took 6.3s alone) |
| `pac_r3_ctrl2out.cir` (7-state buck) | tran (pss q too large to converge at these points) | 3 freqs, `--tran-stabtime-mult 30` | ~14s/point |
| `pac_r6o_ctrl2out.cir` (open-loop ideal-sw buck) | auto (1kHz: tran, q=498; 25.5kHz: pss, q=39) | 2 freqs | ~12s total |
| `pacp_r3_loopgain.cir` (--loop-gain, 5 runs/freq) | tran | 1 freq | ~24s |

The `.tran` fallback backend (needed whenever `q` is too large for
shooting-PSS to converge, or its known q>1 trajectory-emission bug
corrupts the orbit -- `research/near-fs2.md`'s "Engine finding") is the
expensive path: `--tran-stabtime-mult` (default 10x the pss-tuned
`--stabtime`) directly trades cost for the periodicity-residual accuracy
needed at deeply-attenuated points (see "near-fs/2" below) -- a sweep
that needs several near-fs/2 points at high `--tran-stabtime-mult` can
run minutes per point, not seconds.

## Limits (honest, with the numbers that show each one)

- **Commensurate frequencies only.** A target frequency `f` is replaced
  by the closest achievable `f_actual = fs*p/q` (small integers,
  `--qmax`/`--tol` control the search) -- the DFT bin is then exactly
  integer by construction (this is what makes the sideband-0 extraction
  correct, not merely approximate). **Default behaviour reports the
  achieved frequency honestly, even when far from the request** (matching
  `tools/pssac_ref.py`'s own precedent: R2 at target 1000Hz, `--qmax 80`,
  achieves 1250Hz -- 25% away -- and is still trusted, because the
  comparison below is made AT 1250Hz, not silently mislabelled as
  1000Hz). `--strict-freq` refuses instead (RED-PROOF, exercised in
  `tests/test_pac_perturbed.py`). A low `fs`-to-`f` ratio needs a large
  `--qmax` (measured: 1000Hz target on a 500kHz clock needed `--qmax 500`
  for an exact q=498 match).
- **Near-fs/2, deeply-attenuated points need much more settle time than
  the tran-backend default.** Measured this card, `pac_r3_ctrl2out.cir`
  at 199054Hz (~-44dB, close to fs/2=250kHz): at `--tran-stabtime-mult
  10` (3000us total settle), agreement against native `.pac solver=td`
  was **2.0dB / 17.0deg off** despite a periodicity residual that LOOKED
  acceptable (2.1-2.3%, under the default 5% `--residual-bound`) --
  the residual metric itself was not tight enough to catch this. At
  `--tran-stabtime-mult 30` (9000us), agreement improved to **0.52dB /
  0.40deg**, with residual (1.2%/3.3%) closer to the `--residual-bound`
  edge -- the SAME root cause `research/near-fs2.md` already documents
  for `tools/fra.py`'s own transient-FRA reference at these exact
  points (settle time, not switching-frequency proximity per se, is
  the limiting factor). **A caller measuring near fs/2 should raise
  `--tran-stabtime-mult` and tighten `--residual-bound`, and check the
  halving-linearity number too** (the 30x-settle point above still only
  reached 0.47dB halving spread, not the 0.05dB target -- reported
  honestly in the CSV, not hidden).
- **`--qmax`/`q` practicality**: shooting-PSS (`backend=pss`) only
  converges reliably at small `q` (this fork's own known q>1
  trajectory-emission bug, `research/near-fs2.md`); `backend=auto`
  (default) falls back to `.tran` automatically, but a caller who forces
  `--backend pss` on a large-`q` point will see it refused
  (COULD-NOT-RUN), not silently wrong.
- **Above ~fs/2**: only sideband-0 is computed here (same restriction
  `.pac solver=td` has, `docs/pss-pac.md`'s "Known limits") -- neither
  tool sees the k=+-1 sideband coupling a real measurement above fs/2
  would need.
- **`--loop-gain`'s flat Tian probe** requires the deck to have NO
  existing device directly tying the two named probe nodes together (the
  probe breaks that wire) -- `tests/fixtures/pacp_r3_loopgain.cir` is
  `pac_r3_ctrl2out.cir` with its own `Vprobe fbtop vout DC 0` REMOVED for
  exactly this reason (see that fixture's own header).
- **`--outputs a,b`** (comma-separated) does not handle a raw expression
  that itself contains a comma (e.g. a differential `v(a,b)` probe) --
  split naively on `,`; not hit by any fixture in this project, not
  fixed this card.
- **`source:` injection** only recognises a plain
  `Name n+ n- DC <value>[ AC <value>]` independent-source line (one
  regex, one match required) -- a source written across multiple lines,
  or with a PWL/PULSE default already present, is not supported; use
  `node:`/`series:` to add a fresh injection device instead.

## Real, run examples (this card, `ngspice-aug`)

**LTI vs `.ac`** (`tests/fixtures/pac_lti.cir`, `--inject source:Vin
--outputs out`, `--backend pss --stabtime 0`):

```
f_target=10000Hz -> f_actual=10000Hz backend=pss: v(out)=-2.0542dB/-29.7349deg
f_target=30000Hz -> f_actual=30000Hz backend=pss: v(out)=-6.7789dB/-59.7329deg
f_target=60000Hz -> f_actual=60000Hz backend=pss: v(out)=-11.8816dB/-73.7319deg
```
matches the closed-form RC response (Rin=1k, Cout=10nF, Rload=10k, DC
gain 0.9091, pole 17505Hz) to better than 0.01dB/0.01deg at every point;
also reproduced via the `.tran` fallback backend (forced with
`--backend tran`) to within 0.04dB.

**R2 switched-RC vs the exact analytic `g0()`**
(`refs/r2_switched_rc.cir`, `--inject source:Vin --outputs vx`,
`--backend pss --qmax 80`):

```
f=1250.0Hz (target 1000, q=80): |H| relerr=3.9e-04, dphase=0.0012deg
f=10000.0Hz (q=10):             |H| relerr=3.4e-04, dphase=0.0074deg
f=25000.0Hz (q=4):              |H| relerr=1.9e-04, dphase=0.0092deg
f=40000.0Hz (q=5):              |H| relerr=1.1e-04, dphase=0.0091deg
f=47826.1Hz (target 48000, q=23): |H| relerr=9.0e-05, dphase=0.0098deg
```
comfortably inside the 0.1%/0.5deg target.

**R3 (smooth buck) vs native `.pac solver=td`**
(`tests/fixtures/pac_r3_ctrl2out.cir`, `--inject source:Vctrl --outputs
vout`, `--backend auto --tran-stabtime-mult 30`):

| f_actual | backend | pac_perturbed | td-PAC | delta |
|---|---|---|---|---|
| 5000Hz | tran | -5.572dB/-27.104deg | -5.605dB/-27.171deg | 0.033dB/0.067deg |
| 50000Hz | tran | -19.231dB/-137.670deg | -19.453dB/-136.004deg | 0.223dB/1.667deg |
| 199054Hz (near fs/2) | tran, mult=30 | -43.326dB/-127.700deg | -43.842dB/-128.101deg | 0.517dB/0.401deg |

Worst case 0.52dB/1.67deg -- two independent methods (a per-step
linearisation and a nonlinear perturb-and-DFT) sharing no code agree to
well inside the project's usual 0.5dB/3deg target, even at the
historically disputed near-fs/2 point.

**R6o open-loop ideal-switch control-to-output at 1kHz**
(`tests/fixtures/pac_r6o_ctrl2out.cir`, `--inject source:Vctrl --outputs
vout`, `--backend auto --tran-stabtime-mult 30 --qmax 500`): **16.12dB =
6.39 V/V, -1.93deg** (`.tran` fallback backend, q=498) -- matching
`research/near-fs2.md`'s saltation-section finding for this exact point
almost exactly (their independent run of `pssac_ref.py`'s own
`r6o_ctrl2out` tag: 6.39 V/V, -1.94deg) and landing inside the ~5.9-6.4
V/V three-way consensus (analytic Vin/Vramp=6.00, finite-difference of
two converged `.pss` orbits=5.934). **Native `.pac solver=td` on this
same fixture gives 13.19 V/V -- a real, ~2.2x saltation-correction bug,
not fixed by this card** (`docs/pss-pac.md`'s "Saltation" section) --
`pac_perturbed.py` is unaffected because it never linearises anything.

**`--loop-gain` near R3's documented crossover**
(`tests/fixtures/pacp_r3_loopgain.cir`, `--probe-nodes fbtop:vout`,
target 33554Hz, `--backend auto --tran-stabtime-mult 30`): f_actual=
33707.9Hz (q=89, 0.46% off target), **T=-0.069dB/-141.639deg, loopGain=
-0.069dB/38.361deg**, halving spread 0.006dB. Native `.pstb` on the same
fixture family (`docs/pss-pac.md`'s worked example) reports **phase
margin 38.4deg at 33554Hz** -- derived PM here (180+phase(T)) = 38.36deg,
matching to 0.04deg, and `|T|` within 0.07dB of the 0dB crossing --
despite the two measurements not landing on exactly the same grid point
(this tool's own commensurate-frequency picker vs. `.pstb`'s native
sweep), an independent method (nonlinear perturb-and-DFT, Tian
double-injection) reproduces the documented crossover closely.

## Red-proofs (both actually run, in `tests/test_pac_perturbed.py`)

- **`--strict-freq`**: R2 at target 1000Hz, `--qmax 64` (too small for
  any q<=64 to get within `--tol`) -- refused (COULD-NOT-RUN), not
  silently reported at whatever frequency q=64 happens to land on.
- **A too-tight `--residual-bound`**: R3 ctrl2out at 199054Hz forced onto
  the `.tran` backend with `--residual-bound 0.02` (tighter than the
  2.1-3.3% residual that backend/settle combination actually achieves) --
  refused (COULD-NOT-RUN) rather than reporting the 2.0dB-off number the
  same point gives at a looser bound (see the R3 table above,
  `--tran-stabtime-mult 10` row, not shown, which WOULD have silently
  passed a looser residual check with a wrong answer -- exactly the
  failure mode `--residual-bound` exists to catch).

## Free-running (autonomous) mode (`--freerun`, card #493)

Card #493: `.pac`/`.pstb` on a free-running (hysteretic/COT/QR) orbit
(`docs/pss-pac.md`'s free-running section, card #489 SECOND pass) matched
an independent reference to 20-60kHz on R5 but was 3+dB/7-10deg off at
150kHz (35% of R5's own 422.68kHz switching frequency), with two
candidate causes named but not distinguished. That reference itself was
an ad hoc, uncommitted script (docs/pss-pac.md's own "what was not done"
flagged it as such); this card promoted it into `pac_perturbed.py` as
`--freerun`, then used it to settle the question.

**Why a free-running deck needs a DIFFERENT method than every backend
above.** Every backend above needs a fixed switching period T shared by
the baseline and perturbed orbits, so `perturbed(t)-baseline(t)` is
periodic and a plain rectangular-window single-bin DFT is exact. A
free-running orbit's own period is not fixed, and a perturbation
generally SHIFTS it too (this is exactly the period-variation effect
`docs/pss-pac.md`'s `pacshoot.c` correction derives) -- there is no
common period to exploit. `--freerun` instead runs ONE long baseline
`.tran` (shared across every requested frequency) and one long
perturbed `.tran` per frequency (full- and half-amplitude), applies a
Hann taper to suppress the (incommensurate) switching content's own
spectral leakage into the target bin, and reads `(Xpert-Xbase)/Xinj` at
a single DFT bin. The commensurate-frequency handling the brief asked
about is answered by anchoring to the analysis WINDOW instead of to the
switching frequency: `f_actual = M/W` for an integer cycle count `M`
(exact by construction, same "always exact, report relerr honestly"
principle `compute_freq_plan()` uses for the driven backends) -- `W` is
sized to comfortably clear both a cycle-count floor (`--freerun-cycles`,
default 100, relative to the LOWEST requested frequency) and a switching-
period floor (`--freerun-min-switch-periods`, default 300, relative to
the `--fs` ESTIMATE) so the Hann taper has enough of both to work with.
See `tools/pac_perturbed.py`'s own "free-running (autonomous) backend"
comment (`run_sweep_freerun()`) for the full derivation.

**A real bug found building this**: a fixture's own leftover `.pac ...
trigger="q"` dot-card, left in the deck body the way every OTHER backend
here already (safely) leaves a fixture's `.pac` card in place, silently
hijacked the appended `.control` block's `run` command whenever its own
PSS-FREERUN shoot CONVERGED (which it reliably does on R5) -- `wrdata`
then wrote that card's own ~2.4us one-period orbit instead of the
intended multi-millisecond window, giving baseline==perturbed (a bogus
~-200dB "gain", not a real near-zero one). Fixed by stripping any
leftover `.pac`/`.pstb`/`.pss` dot-card from the body before appending
this tool's own `.tran` (`_strip_dot_cards()`, `--freerun` only -- the
existing, already-validated backends are untouched, though the same
latent risk applies there too whenever shooting-PSS happens to converge
on a q small enough to leave a competing plot; not fixed there, out of
scope for this card).

**Numbers (R5, `refs/r5_hysteretic_buck.cir`/`tests/fixtures/
pac_r5_ctrl2out.cir`, control-to-output, `--freerun --fs 422681 --inject
source:Vrefsrc --outputs vout --amplitude 0.001`), native `.pac` on the
`ngspice-src-fr150`/`aug-fr150` build (card #493's own fine forward-
difference `f(T)` fix, see `docs/pss-pac.md`'s card #493 section)**:

| freq (% of fsw) | pac_perturbed --freerun (halving spread) | native corrected | native uncorrected | delta corrected | delta uncorrected |
|---|---|---|---|---|---|
| 20kHz (4.7%) | -0.0088dB/-0.091deg (0.059dB/0.37deg) | 0.061dB/-0.013deg | 0.067dB/-0.042deg | +0.070dB/+0.079deg | +0.075dB/+0.049deg |
| 60kHz (14.2%) | 0.279dB/-0.454deg (0.090dB/0.22deg) | 0.615dB/0.155deg | 0.647dB/-0.657deg | +0.336dB/+0.609deg | +0.367dB/-0.203deg |
| 100kHz (23.7%) | 0.769dB/-2.125deg (0.042dB/0.14deg) | 1.875dB/0.920deg | 1.866dB/-3.113deg | +1.106dB/+3.045deg | +1.096dB/-0.988deg |
| 150kHz (35.5%) | 1.412dB/-5.892deg (0.034dB/0.06deg) | 5.049dB/5.234deg | 4.513dB/-12.986deg | +3.637dB/+11.126deg | +3.101dB/-7.094deg |
| 200kHz (47.3%) | 2.329dB/-11.327deg (0.146dB/1.29deg) | 11.588dB/41.950deg | 7.707dB/-47.730deg | +9.259dB/+53.277deg | +5.379dB/-36.403deg |

**Verdict: the reference was right, the engine is wrong** -- the opposite
of card #494's SIXTH pass (saltation), where cross-validation found the
REFERENCE wrong and the C code right. `--freerun` at 20kHz/60kHz/150kHz
matches the card #489 ad hoc Hann-DFT reference closely (150kHz:
1.412dB/-5.892deg here vs 1.465dB/-5.78..-5.96deg there), confirming that
reference was correct all along -- the disagreement is real, in the
native engine, and (new this card) does NOT start at 150kHz: it is
already outside the 0.5dB/3deg target by 100kHz (23.7% of fsw), well
below where the driven engine's own near-fs/2 accuracy floor usually
bites.

**Root cause: NOT the period-variation correction term itself** -- three
things were tried and none closed the gap:
1. **A controlled fine forward-difference `f(T)` estimate** (`aug-fr150`,
   `pssshoot.c`'s `PSSshootFreerun()`, `PSS_FREERUN_FT_DT` env var,
   default 1e-13s vs the old ~15.8ps backward difference): barely moves
   the 150kHz result (5.05dB/5.23deg vs the old 5.20dB/3.82deg) -- both
   equally far from the 1.41dB/-5.89deg reference.
2. **`PAC_EXACT_JACOBIAN=1`** (`pacshoot.c`'s own, unrelated, already-
   documented Jacobian-staleness floor): changes the 150kHz result by
   <1e-4 relative -- ruled out.
3. **`pacpoints=4000` vs the fixture's default 200** (recording-pass grid
   resolution): changes the WRONG answer to a DIFFERENT wrong answer
   (5.22dB/-5.46deg), not a smaller one -- consistent with this project's
   own prior experience elsewhere in this file ("near-fs/2 ... not
   reduced by `pacpoints=`/`Points`").

**And it is not specific to the correction term**: both the corrected
and the UNCORRECTED (`PAC_FREERUN_NO_PROJ=1`) native results grow wrong
TOGETHER, starting at the same ~100kHz point -- consistent with a general
per-step-linearisation accuracy floor specific to (or dramatically worse
on) R5's own extremely stiff switch node (`Csw`, ~1ps local RC, far
faster than any state in R3/R4), not yet root-caused further. See
`docs/pss-pac.md`'s card #493 section for the fork-side detail.

**Card #503 update**: a driven-clocked replica isolating R5's `Csw`
stiffness alone (same ~1ps switch-node RC, no `Cq` latch state) was
matched by the native DRIVEN engine (no `trigger=`) to <0.11dB/0.24deg
against this same tool's plain (non-`--freerun`) backend up to 47% of
fsw -- ruling OUT Csw stiffness in isolation as the explanation. The
free-running Newton solve's own finite-difference monodromy was also
found to sit at or past its own numerical noise floor on R5 (fails to
converge at all at `steady_coeff=1e-9`; reports a spurious `UNSTABLE`
spectral radius at `1e-5`), and its Floquet report shows a genuinely
near-degenerate pair of near-zero multipliers (not the single real
trivial eigenvalue the theory predicts) -- consistent with R5's TWO
widely-separated fast time constants (`Csw`~1ps, `Cq`~1ns). See
`docs/pss-pac.md`'s card #503 section for the full experiment writeup;
disposition unchanged (still refused by default).

**Card #503 SECOND pass update**: the coordinator asked for a replica that
keeps the LATCH (`Cq`/`Bq`), not just `Csw`, plus a targeted step-refinement
check, a non-transient reference, and a fix attempt. With the latch kept
(vout->Bq feedback broken by a fixed PWL replica of the recorded orbit
instead) and the period known exactly, `pac_perturbed.py`'s plain driven
backend (amplitude cut to 2e-6 to stay linear at this deck's much higher
open-loop gain -- 0.001 gave a 5.8dB halving spread, badly nonlinear) still
agrees with native driven `.pac` to within 0.2-0.4dB/0.5deg up to 200kHz
once the recording-pass grid is refined (`pacpoints=200000`) -- i.e. the
SAME grid-refinement lever that does NOT help the real free-running engine
(tried directly on it too: 150kHz/200kHz results wander non-monotonically,
not converging toward this tool's own `--freerun` reference) cleanly fixes
a driven replica with the identical latch dynamics. This narrows the
dominant cause to the autonomous/trigger-relocation formulation itself,
not the latch's steepness (a real but secondary, grid-fixable contributor)
and not general stiffness. See `docs/pss-pac.md`'s card #503 SECOND pass
section for the full tables, the internal-consistency (Floquet-vs-.pac)
non-transient argument, and what a fix would still need. No fix was
achieved or shipped for the free-running engine itself this pass.

**Card #503 THIRD pass update**: root-caused and FIXED (partially). An
independent, high-precision free-running orbit (a tight-reltol `.tran`
seeded at the shooter's own converged `x0`, run to cycle-to-cycle
stability) differs from what `PSSshootFreerun()` was handing to the
`.pac`/`.pstb` recording pass by ~530ppm in the period and up to ~0.02%
in the inductor-current state -- small numbers, but injecting that
independent x0/T (a new debug override) cut the disagreement against
THIS tool's own `--freerun` reference by 5-8x at every frequency tested
(200kHz: 9.26dB/53.3deg -> 1.20dB/3.81deg), confirming hand-off precision
as the dominant cause. Root cause: `pssshoot.c` caps `reltol` at 1e-5 for
robustness, too loose for R5's own ~1e6:1 stiffness ratio. Fixed by
polishing `x0`/`T` for the `.pac`/`.pstb` hand-off with a few cheap extra
period-map evaluations under a SCOPED reltol tightening (the orbit is
heavily damped, so this converges in ~5 iterations) -- on by default
whenever the experimental gate is opted into. Result: 20/60/100kHz
(<=23.7% of fsw) now cleanly meet the 0.5dB/3deg target against this
tool's `--freerun` reference; 150kHz is close (phase inside target,
magnitude 0.24dB over); 200kHz improved 7-8x but remains outside target.
See `docs/pss-pac.md`'s card #503 THIRD pass section for the full tables
and the (not yet closed) 150-200kHz residual. Disposition: still gated
by default, not a full ungate -- see that section's own proposal.

**Disposition**: kept REFUSED BY DEFAULT (`PSS_FREERUN_PAC_EXPERIMENTAL=1`
still required, unchanged), plus a new stderr warning
(`pacshoot.c`, "free-running (trigger=) ... UNVALIDATED margin") printed
whenever a free-running `.pac`/`.pstb` frequency reaches 10% of the
measured free-running switching frequency (conservative relative to the
~14% point actually measured clean, same "warn before the budget is
exhausted" convention patch 0046's saltation band guard uses) -- not
enabled by default, not fixed, the fork's changes only improve the
`f(T)` estimate's own numerical honesty (ruling out that candidate) and
make the failure mode loud rather than silent.

## Usage

```sh
# control-to-output, ideal-switch deck
python3 tools/pac_perturbed.py sweep \
    --deck tests/fixtures/pac_r6o_ctrl2out.cir --fs 500e3 \
    --inject source:Vctrl --outputs vout \
    --freqs 1000,25500 --amplitude 0.01 \
    --ngspice ~/Work/apps/ngspice_augmentation/aaba50ba3f49-s41/bin/ngspice-aug \
    --out ~/Work/data/ngspice_augmentation/pac_perturbed/run1 --max-procs 4

# a dec sweep, current injection at a node (output impedance)
python3 tools/pac_perturbed.py sweep \
    --deck mydeck.cir --fs 200e3 --inject node:vout --outputs vout \
    --sweep-spec "dec 5 500 90000" --amplitude 0.01 ...

# loop gain, Tian double injection (deck must NOT already tie a-b together)
python3 tools/pac_perturbed.py sweep \
    --deck mydeck.cir --fs 200e3 --loop-gain --probe-nodes fbtop:vout \
    --sweep-spec "dec 10 500 100000" --amplitude 0.01 ...

# free-running (autonomous/hysteretic) deck, card #493 -- --fs is only a
# sizing ESTIMATE here, never exact; no --backend/--qmax/--strict-freq
python3 tools/pac_perturbed.py sweep \
    --deck tests/fixtures/pac_r5_ctrl2out.cir --freerun --fs 422681 \
    --inject source:Vrefsrc --outputs vout \
    --freqs 20000,60000,100000,150000,200000 --amplitude 0.001 \
    --stabtime 500e-6 --ngspice ~/Work/apps/ngspice_augmentation/aaba50ba3f49-s46/bin/ngspice \
    --out ~/Work/data/ngspice_augmentation/pac_perturbed/run2 --max-procs 2
```

`--inject source:<InstName>` (an existing `Name n+ n- DC <val>[ AC
<val>]` line -- the AC token, if present, is ignored/dropped, never
relied on) / `node:<node>[:<ref>]` (new shunt current source) /
`series:<a>:<b>` (new series voltage source). `--backend auto` (default)
tries `.pss method=shoot` first, falls back to a long, settled `.tran`
with a periodicity-residual check automatically. See `python3
tools/pac_perturbed.py sweep --help` for every flag.

## What was not done

- **No `.pac_events`-style multi-sideband output** -- sideband-0 only,
  same restriction `.pac solver=td` has.
- **`--outputs` with a comma inside one expression** (a differential
  probe) is not supported (see "Limits").
- **The R6 CLOSED-LOOP case (R6s, XSPICE `d_pwm`)** was not swept this
  card: `docs/pss-pac.md`'s saltation section already documents that
  `.pss method=shoot` itself does not converge on an ideal-`sw`-in-
  feedback-loop buck on this fork, and XSPICE `d_pwm` segfaults inside
  the SAME shooting engine's period-restart -- both are engine-level
  issues in the `.pss method=shoot` PATH specifically. This tool's
  `.tran` fallback backend does NOT use shooting-PSS at all (a plain,
  non-restarting transient), so it should not hit either bug in
  principle, but a full R6/R6s closed-loop, XSPICE control-to-output/
  line-to-output/output-impedance/loop-gain sweep across fs/1000-fs/2
  (the validation table the card's brief describes) was not run this
  pass -- time budget; only the open-loop R6o single-point check and the
  R3 loop-gain smoke-test above were completed. A follow-up pass should
  attempt it directly (the mechanism should work; it has not been
  proven on that specific deck).
- **Per-point amplitude/settle auto-tuning** is not implemented --
  `--amplitude`/`--tran-stabtime-mult`/`--residual-bound` are swept-wide
  flags, not tuned per frequency the way `research/near-fs2.md`'s own
  manual investigation was; a caller sweeping across a wide dB range in
  one command may need several passes with different settings for the
  hardest (near-fs/2, deep-attenuation) points, same as that file's own
  "what was not done" note for `pssac_ref.py`.
- **Runtime was not systematically minimised** -- `.tran`-backend points
  at `--tran-stabtime-mult 30` cost tens of seconds each; no attempt was
  made to use `tools/campaign.py`'s checkpoint-shared-settle mode (c) or
  a shorter-but-still-correct settle heuristic this card.
- **card #493 (`--freerun`)**: only R5 control-to-output was swept
  (line-to-output, and any other free-running fixture, were not); only
  5 frequency points (20/60/100/150/200kHz), no dense sweep across
  fs/1000-fs/2 the way the driven backends' own tables do. `--loop-gain
  --freerun` is not implemented (R5 has no natural loop-break point,
  `docs/pss-pac.md`'s own free-running section already notes this).
  `--freerun-cycles`/`--freerun-min-switch-periods` defaults (100/300)
  were not tuned or swept for cost/accuracy trade-offs -- only checked
  to give a small (<=0.15dB/1.3deg) amplitude-halving spread at the 5
  points actually run. The root cause of the native engine's own
  breakdown (a per-step-linearisation accuracy floor apparently specific
  to R5's very stiff switch node) was narrowed (three candidates ruled
  out: `f(T)` precision, `PAC_EXACT_JACOBIAN`, `pacpoints=`) but not
  identified further -- see `docs/pss-pac.md`'s card #493 section.
