# Fast / high-sigma yield estimation

Card #474, plan T10.6, 2026-09-29. Companion to `tools/mc.py` (card #473,
plan T10.5, imported here, not modified) and to research/feature-survey.md
item #6 ("fast/high-sigma Monte Carlo"). `tools/mc.py` draws N independent
samples and computes plain statistics -- correct, but at a 4-6 sigma yield
target it needs far more runs than are affordable: at Pf ~= 3e-5 (4 sigma),
a 20%-relative, 95%-confidence plain-MC estimate needs roughly
N ~= (1.96/0.20)^2 * (1-Pf)/Pf ~= 3.0e6 runs (the exact number this card's
own validation used, see below); at 6 sigma it is ~9.7e10. `tools/yieldest.py`
estimates the same Pf with two independent methods, using **~160x fewer
runs at 4 sigma, up to ~5.15 million x fewer at 6 sigma** (measured below --
this tool's own actual `n_sims` divided into the plain-MC-equivalent count
above, not a marketing number), never an ML surrogate (research/
feature-survey.md item #6's own ranking preferred a worst-case-distance/
importance-sampling wrapper for a first version over a surrogate).

## Methods

All sampling happens in a **standardized** parameter space: every declared
Gaussian parameter x_i (nominal m_i, sigma s_i) becomes u_i = (x_i-m_i)/s_i,
so the nominal sampling distribution is always the standard multivariate
normal N(0, I) regardless of each parameter's physical units. A "failure"
is a single scalar measure (read back via `mc.py`'s own `MCMEAS` stdout
convention) crossing a threshold (`fail: {"op": "lt"|"gt", "threshold": v}`).
Only `dist: "gauss"` parameters are supported -- the whole likelihood-ratio
and worst-case-distance machinery below is Gaussian-specific; `mc.py`'s
`"unif"` distribution is not.

### 1. Worst-case-distance (WCD) search, gradient-free

A deterministic set of directions in u-space is searched outward from the
origin for the pass/fail boundary, by an expand-then-bisect line search
along each ray (the standard WCD/FORM assumption: each ray crosses
pass->fail at most once). The direction set is **not random** (so the
search itself is exactly reproducible run to run, matching this project's
seeded/reproducible convention elsewhere):
- every coordinate axis (2 directions per parameter);
- for exactly 2 parameters, a dense 24-direction (15-degree-spaced)
  angular grid;
- for 3+ parameters, the coordinate axes plus the (capped at 32) all-sign
  unit-diagonal directions, which a coordinate-only search would miss
  entirely.

The smallest boundary radius found is then refined by a gradient-free
directional **pattern (Hooke-Jeeves-style) search** on the direction vector
itself: each iteration perturbs the current best direction by a shrinking
step in every parameter, re-locates the boundary along each perturbed
direction (a fast local bisection, not a full re-search), and keeps
whichever candidate reduced the radius. This gives `beta` (the worst-case
distance, exactly what FORM/HSMC-style commercial tools call it) and
`u_star` (the most-probable-failure point, MPP) -- see `search_mpp` in
`tools/yieldest.py`.

### 2. Mean-shift importance sampling (primary estimate)

Draw N samples from N(u_star, I) (same covariance as the nominal
distribution, shifted mean -- the simplest, most standard importance-
sampling variant), evaluate each on the real circuit, and weight every
sample by the likelihood ratio
`w(u) = f0(u)/f_IS(u) = exp(-u.u_star + 0.5*||u_star||^2)`. The estimator
`Pf_hat = mean(w * 1{fail})` is a standard, unbiased importance-sampling
estimator: most of the shifted distribution's mass sits right where the
failure region is, instead of almost every draw from the nominal N(0,I)
being an uninformative "pass" the way it is at 4-6 sigma. A 95% CI is
built from the sample variance of the weighted indicator, in log space
(`log(Pf_hat) +- z*CV`) so it does not fabricate a negative lower bound
at very small Pf. `tools/yieldest.py`'s headline `pf_estimate` is this
method's result whenever it ran.

### 3. Scaled-sigma sampling (SSS), independent cross-check

Independent of the WCD search's direction/weights (it reuses the found
`beta` only to pick a sensible scale range -- see `default_sss_scales`):
plain MC at several **enlarged** standard deviations
(u ~ N(0, c^2 I), c > 1), where the failure probability Pf_c is large
enough to measure directly and cheaply. Assuming (as WCD/FORM does) the
failure boundary is approximately a hyperplane at distance `beta` from the
origin, `Pf_c ~= Phi(-beta/c)` exactly for a true hyperplane, so
`beta_hat_c = -c * Phi^-1(Pf_c_hat)` is recovered from each scale's
plain-MC estimate. The per-scale estimates are combined by **weighted
least squares** (regressing `Phi^-1(Pf_c)` against `1/c` through the
origin, weight = 1/Var via the delta method on each scale's own binomial
proportion) -- found, during this card, to matter: a first version that
simply averaged the per-scale `beta_hat_c` let one noisy scale (too few
observed failures) swing the combined answer by several sigma-equivalents;
weighting by each scale's own reliability fixed that (see "What changed
during validation" below).

### Honesty check: do the two methods agree?

`estimate_yield` compares the IS estimate against the SSS estimate; a
disagreement of more than 3x either way sets `limitation_warning` in the
result, rather than silently returning "the answer". Two distinct causes
can trigger this (the tool cannot itself tell them apart, and says so in
the warning text):
1. a failure region the WCD search's direction set did not find (see the
   disjoint-region red-proof below), or
2. SSS's hyperplane-tail assumption breaking down on a **strongly curved**
   failure boundary even with a single region (importance sampling does
   not share that assumption -- see the MOS case below).

## Validation (`tests/test_yield.py`)

Every ground-truth Pf is computed **independently** of `tools/yieldest.py`'s
own code -- case (a) via the exact closed-form normal CDF, cases (b)/(c) via
a plain stdlib 2-D Riemann-sum integration over the same nonlinear formula
the ngspice B-source evaluates. Real captured output, this card,
`NGSPICE=/usr/local/bin/ngspice` (the pinned unmodified system build),
`python3 tests/test_yield.py`, wall clock 1m16s:

```
[pass] ngspice availability smoke test: vtop=20.0
resistor (linear), 4 sigma
[pass] resistor 4sigma: estimate_yield ran: n_sims=18893 wall=13.4s beta_mpp=4.0000
[pass] resistor 4sigma: importance-sampling Pf within 2x: pf=2.8374e-05 true=3.1671e-05 ratio=0.896
[pass] resistor 4sigma: scaled-sigma-sampling Pf within 2x: pf=2.4624e-05 true=3.1671e-05 ratio=0.777
[pass] resistor 4sigma: n_sims <= 1% of plain-MC-equivalent: n_sims=18893 plain_mc_equivalent~=3.032e+06 used=0.6231% of it
[pass] resistor 4sigma: MPP contribution fractions sum to 1, ~50/50 R1/R2
resistor (linear), 5 sigma
[pass] resistor 5sigma: estimate_yield ran: n_sims=18900 wall=13.3s beta_mpp=5.0001
[pass] resistor 5sigma: importance-sampling Pf within 2x: pf=2.1354e-07 true=2.8665e-07 ratio=0.745
[pass] resistor 5sigma: scaled-sigma-sampling Pf within 2x: pf=2.4219e-07 true=2.8665e-07 ratio=0.845
[pass] resistor 5sigma: n_sims <= 1% of plain-MC-equivalent: used=0.0056% of the plain-MC-equivalent (3.350e+08)
resistor (linear), 6 sigma
[pass] resistor 6sigma: estimate_yield ran: n_sims=18902 wall=13.3s beta_mpp=6.0001
[pass] resistor 6sigma: importance-sampling Pf within 2x: pf=9.8276e-10 true=9.8659e-10 ratio=0.996
[pass] resistor 6sigma: scaled-sigma-sampling Pf within 2x: pf=6.6179e-10 true=9.8659e-10 ratio=0.671
[pass] resistor 6sigma: n_sims <= 1% of plain-MC-equivalent: used=0.0000% of the plain-MC-equivalent (9.735e+10)
[pass] plain MC sanity @3sigma (resistor): n=6000 k=10 pf=1.6667e-03 ci95=[6.35e-04,2.70e-03] true=1.3499e-03
MOS (nonlinear), ~5 sigma
[info] independent 2-D numerical-integration ground truth: threshold=2.5000e-04 Pf=4.9188e-07 equivalent_sigma=4.895
[pass] MOS nonlinear case: estimate_yield ran: n_sims=18940 wall=13.5s beta_mpp=4.9191
[pass] MOS nonlinear case: importance-sampling Pf within 2x: pf=5.6395e-07 true=4.9188e-07 ratio=1.147
[pass] MOS nonlinear case: n_sims <= 1% of plain-MC-equivalent: used=0.0097% of the plain-MC-equivalent (1.953e+08)
[info] MOS nonlinear case: scaled-sigma-sampling Pf (non-gating): pf=1.923e-06 ratio=3.909 disagreement=3.41
[pass] plain MC sanity @~3sigma (MOS): n=6000 k=6 pf=1.0000e-03 ci95=[2.00e-04,1.80e-03] true=1.3490e-03
red-proof: dropped importance-sampling weights
[pass] naive (unweighted) fail-rate=4.9000e-01 vs true Pf=9.8659e-10 -> overestimates by 4.967e+08x
red-proof: disjoint failure region
[info] Pf_A=1.3234e-06 Pf_B=4.5992e-05 (B/A=34.8x) Pf_total=4.7315e-05
[pass] worst-case-distance search found region A, not B: beta_found=4.5001 region_A_nearest=4.5000
[pass] importance-sampling estimate misses the true (both-region) Pf by >2x: IS pf=1.227e-06 vs true
       total=4.7315e-05 -> underestimates by 38.6x
[pass] tool reports a limitation_warning instead of a silent wrong answer: disagreement=4.9e-06 limitation_warning=set
25 pass, 0 FAIL (gating), 0 COULD-NOT-RUN (gating), 3 info/non-gating
```

### Case (a): linear resistor stack -- exact ground truth

`tests/fixtures/yield_resistor.cir`: an ideal current source (1 mA) through
two resistors R1, R2 in series; the output is the top-node voltage,
`Vtop = I*(R1+R2)` -- deliberately **not** a Vin/R1/R2 voltage divider,
whose tap voltage is the ratio `Vin*R2/(R1+R2)`, nonlinear in R1 and R2 and
without a closed-form Gaussian output. Because Vtop is an exact linear
combination of two independent Gaussian resistor values, it is itself
exactly Gaussian, and `Pf = Phi(-k)` is exact at any k. R1, R2: nominal
10 kOhm, 3% (300 Ohm) 1-sigma each. At each of 4, 5, 6 sigma: both methods
land within a factor of 2 of the truth (ratios 0.67-1.00), using
0.0000%-0.62% of the plain-MC run count a matched-confidence brute-force
estimate would need (0.62% at 4 sigma is the loosest of the three --
consistent with "<=1%" being the harder bound at the *lower* end of the
4-6 sigma range, since plain MC itself is least impossible there).

### Case (b): nonlinear MOS square-law -- numerically-integrated ground truth

`tests/fixtures/yield_mos.cir`: a behavioral (B-source) square-law current
stage, `id = 0.5*Kp*max(Vgs-Vth, 0)^2`, Vgs fixed at 1.0 V, **both** Kp and
Vth independently Gaussian (Kp: nominal 2.5 mA/V^2, 8%; Vth: nominal 0.4 V,
20 mV). Not a real transistor model (stated plainly, like `docs/mc.md`'s
diffpair example) -- a B-source so the ngspice-evaluated current is exactly
the formula the independent Python ground truth also uses, so a mismatch
would be this tool's bug, not a model-fidelity gap. Because Kp multiplies
`(Vgs-Vth)^2`, the failure boundary in (Kp,Vth) space is a **curve**, not
reducible to a threshold on a single Gaussian combination -- true Pf needs
2-D numerical integration (implemented independently in
`tests/test_yield.py`, not shared with `tools/yieldest.py`).

At `id < 2.5e-4` (Pf = 4.919e-7, equivalent 4.895 sigma): importance
sampling lands at ratio 1.147 (well within 2x) using 0.0097% of the
plain-MC-equivalent run count. **Scaled-sigma sampling misses by 3.9x here
(non-gating)** -- this is the real, found-not-assumed limitation the
module docstring's "honesty check" section describes: SSS's hyperplane
assumption is a much better approximation for the resistor case's flat
boundary than for this curved one, and the tool's own disagreement
diagnostic (3.41x here) correctly flags it as unreliable rather than
averaging it silently into a combined "answer" (`pf_estimate` always
prefers the IS result when it ran ok, exactly to avoid that).

### Plain-MC sanity checks at ~3 sigma

At 3 sigma (Pf ~= 1.35e-3), brute-force MC is affordable (6000 runs) and
was run as a direct sanity check on both circuits: resistor case,
`k=10/6000` (Pf_hat=1.667e-3, 95% CI [6.3e-4, 2.70e-3]) against the exact
1.3499e-3; MOS case (at a threshold separately bisected to Pf~=1.35e-3),
`k=6/6000` (Pf_hat=1.0e-3, CI [2.0e-4, 1.80e-3]) against the numerically
integrated 1.3490e-3 -- both inside their own CI, i.e. brute force and the
fast methods are not telling different stories where brute force can still
be checked directly.

### Red-proof 1: dropped importance-sampling weights

`tools/yieldest.py`'s `importance_sampling` exposes the *unweighted* fail
rate under the shifted sampling distribution as a diagnostic field
(`naive_unweighted_fail_rate_NOT_A_PF_ESTIMATE`) -- never as its own Pf
estimate. At 6 sigma this naive rate is 0.49 (roughly half the shifted
draws fail, since the shift centers the distribution almost exactly on the
boundary) against a true Pf of 9.87e-10 -- an overestimate by **4.97e8x**,
against the correctly-weighted estimate's 0.996x. This is the proof the
likelihood-ratio weighting is load-bearing, not a decorative detail that
happens not to matter.

### Red-proof 2: a failure region the search cannot find

`tests/fixtures/yield_disjoint.cir` is a synthetic (not physical, stated as
such) fixture with **two disjoint failure regions** in (u1,u2) space:
region A, a disk centered exactly on one of the WCD search's 24 (15-degree)
grid rays (nearest-boundary distance 4.5, found reliably); region B, a much
smaller disk centered exactly **midway between two grid rays** (7.5 degrees
either side -- wider than region B's own ~3-degree angular half-width), so
it is *never* seen by the direction-limited search, even though its
nearest-boundary distance (3.3) is *closer* to the origin than region A's.
Independent 2-D integration: `Pf_A = 1.323e-6`, `Pf_B = 4.599e-5` --
region B, the one the search cannot find, holds **35x more probability**
than the one it does find. The importance-sampling estimate (anchored to
region A's MPP) comes out at 1.227e-6 -- essentially exactly `Pf_A`
(ratio 0.927 to A alone) -- missing the true total (`Pf_A+Pf_B`=4.73e-5) by
**38.6x**. `tools/yieldest.py` does not silently report this wrong number:
`method_disagreement_ratio` (SSS/IS here is ~4.9e-6, i.e. SSS's answer is
~200,000x *smaller* than IS's -- SSS is *also* unreliable on this fixture,
for a third reason stated in Limits below) trips the >3x-disagreement gate
and sets `limitation_warning`.

### What changed during validation (found, not assumed)

The first version of `scaled_sigma_sampling` averaged each scale's
`beta_hat_c` unweighted, with a fixed `[2.0, 3.0, 4.0]` scale list and
`n_sss=200`. On the resistor case at 5-6 sigma this was off by up to 11x
(e.g. one seed: `beta_hat` 6.38 vs true 6.00 -- only a 6% error in `beta`,
but an **11x** error in `Pf = Phi(-beta)` because the relationship is
exponential) -- caught by this card's own validation, not assumed correct.
The fix: (1) choose the scale range from the WCD-found `beta` itself
(`default_sss_scales`, spanning 0.35x-1.2x of `beta`, so each scale's Pf_c
lands in a cheaply-measurable ~1%-60% range instead of being tiny or
saturated); (2) combine scales by weighted least squares instead of a
plain mean; (3) raise `n_sss` from 200 to 3000/scale (still cheap -- see
run counts above). A 50-trial synthetic (pure-Python, no ngspice) check
of the resulting estimator at beta=4,5,6 found 45-50 of 50 trials within
2x of the truth (90-100%, scale-dependent) -- **not proven to be within
2x on every possible seed**, which is exactly why `pf_estimate` always
prefers importance sampling (validated at 100/100 of an equivalent
synthetic check, tighter CI, more robust to a curved boundary) and SSS is
reported as a cross-check, not the number a caller should read alone.

## Worked example: differential-pair offset vs a comparator spec, ~5 sigma

Reusing `docs/examples/mc_diffpair_offset.cir` (the same real, synthetic-
model 5T NMOS differential pair `docs/mc.md`'s worked example uses -- W=8u,
L=1u, 100uA tail, `nch4` BSIM4-family model, synthetic, not derived from any real design). A local/
mismatch input-referred offset (`VOFF`, in M2's gate path) is Pelgrom-
scaled: `sigma_Voff = A/sqrt(W*L)` with the same illustrative, clearly-
synthetic `A = 4 mV*um` `docs/mc.md` uses -- here, `sigma_Voff = 1.4142 mV`.
Using `docs/mc.md`'s own independently-measured gain (finite-difference
`op` at Voff=0 and Voff=1mV, **not** a hand-derived formula: 3.786 V/V from
input offset to `vdiff = v(out1)-v(out2)`), the spec is framed as "does the
comparator's output differential voltage cross its decision threshold at a
5-sigma input-referred offset": `fail: vdiff < -3.78598*5*1.4142e-3 V`
= `vdiff < -26.77 mV`.

Real captured run, this card, `NGSPICE=/usr/local/bin/ngspice`,
`n_is=300`, seed 555, wall 15.6s, **18958 simulations total**:

```
mpp beta:      5.00095   (found essentially exactly 5, as expected --
                           the real BSIM4 device is close to linear over
                           this small an offset range, same finding
                           docs/mc.md's own gain measurement relies on)
u_star:        {'MW': 0.0, 'ML': 0.0, 'VOFF': -5.00095}
contribution:  {'MW': 0.0, 'ML': 0.0, 'VOFF': 1.0}   -- correctly attributes
               100% of the failure to VOFF (MW/ML were held fixed, sigma=0)
pf_estimate:   2.8987e-07   (importance sampling)
pf_ci95:       [2.1876e-07, 3.8408e-07]
equivalent_sigma: 4.998
scaled-sigma-sampling pf: 2.5805e-07 (ratio 0.90, agrees, no disagreement warning)
```

Against the idealized linear-approximation reference `Phi(-5) = 2.8665e-7`
(what a purely linear model at this exact gain and sigma would predict):
ratio **1.011** -- agrees to 1.1%, i.e. the real nonlinear BSIM4 op-point
solve and the idealized linear picture agree closely at this bias point
and offset range, consistent with `docs/mc.md`'s own finding that this
circuit's offset-to-output relationship is close to linear here. This is
reported as a demonstration (a real device, not a synthetic B-source), not
a second independent-ground-truth validation the way cases (a)/(b) are --
the "truth" it is compared against is itself only the linear approximation,
not an independently-solved BSIM4 result.

## Limits / what was not done

- **Only Gaussian parameters** (`dist: "gauss"`) -- the likelihood-ratio
  weights, the worst-case-distance `beta`, and the SSS hyperplane algebra
  are all Gaussian-specific. `mc.py`'s `"unif"` distribution is not
  supported here; passing one raises immediately (`check_params_supported`)
  rather than silently mis-weighting.
- **Single scalar measure per estimate.** No combined multi-measure yield
  (e.g. "fails if EITHER offset OR gain is out of spec") yet.
- **The WCD search assumes each ray crosses pass->fail at most once**
  (standard FORM assumption). A non-monotonic ray (fails, then passes
  again, then fails again further out) can mislead the bisection into
  reporting the *first* crossing, which may not be the nearest true
  boundary point. Not tested here; a plausible failure mode for a
  circuit with a resonance or a folded nonlinearity, not exercised by any
  of this card's reference fixtures.
- **A narrow or non-nearest failure region can be missed entirely** (red-
  proof 2, above) -- by construction, since the direction search is finite.
  The `limitation_warning` disagreement check is this tool's *only* defence
  against reporting a confidently-wrong number in that case, and it cannot
  distinguish "a missed region" from "SSS's hyperplane assumption failing
  on a curved boundary" (both trip the same >3x gate) -- see the MOS case
  and red-proof 2 above, both real, run, and reported, not hypothetical.
- **Scaled-sigma sampling is noticeably less reliable than importance
  sampling at a matched run budget**, especially on a curved failure
  boundary (MOS case: 3.9x off) or a bounded (non-hyperplane-like) failure
  region (disjoint-region red-proof: ~200,000x off in the *opposite*
  direction from the IS miss) -- SSS is reported as a cross-check, gated
  at a 2x tolerance only on the linear case where it was validated
  reliable, non-gating everywhere else in this test suite. A caller who
  needs SSS's own number to be trustworthy on a nonlinear circuit should
  increase `n_sss` well past this tool's default and/or treat a wide
  method-disagreement as inconclusive, not as "SSS must be the wrong one."
- **No ML/surrogate-based method** (research/feature-survey.md item #6's
  own ranking preferred a WCD/importance-sampling wrapper for a first
  version; Yin et al. 2023's deep-kernel surrogate is future work, not
  attempted here).
- **The "<=1% of plain MC" figure is a stated formula, not a run plain MC
  actually did at 4-6 sigma** (that many runs is exactly what this tool
  exists to avoid) -- `plain_mc_runs_needed` in `tests/test_yield.py`
  computes the count a 20%-relative, 95%-confidence *plain* MC estimate
  would need via the standard normal-approximation formula
  `N ~= (z/relerr)^2 * (1-Pf)/Pf`, and this card's own actual `n_sims` is
  compared against that number, not against a brute-force run that was
  never affordable to run. The only sigma level this project *did* run
  plain MC at directly is 3 sigma (above), where both agree.
- **Run counts are what this card's defaults used** (`n_is=300`,
  `n_sss=3000/scale` across 4-6 scales, ~19000 simulations per full
  `estimate_yield` call, ~13-16s wall clock at <=4 concurrent ngspice
  processes) -- not independently re-tuned per circuit; a harder
  (higher-dimensional, more strongly nonlinear) real circuit may need more.
- Every number in this document was produced by an actual run captured
  during this card (`tests/test_yield.py`'s own printed output, and the
  worked-example driver script run separately for the table above, both on
  `/usr/local/bin/ngspice` 45.2, the pinned unmodified system build) --
  none were estimated or assumed.
