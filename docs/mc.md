# Monte Carlo: global (process) vs local (mismatch) variation

Card #473, plan T10.5, 2026-09-29. Companion tool to `tools/loss.py` and
`tools/stb.py` in spirit (a script layer over an unmodified ngspice), but
for statistical yield/offset work: `tools/mc.py` runs a Monte Carlo
campaign that correctly separates a **global/process** parameter (moves
every device in a run together) from a **local/mismatch** parameter (an
independent draw per device instance, optionally Pelgrom-scaled by
`1/sqrt(W*L)`), on the **unmodified system ngspice-45.2**. No C change was
made or was found necessary for this card -- see "Why a tools-level
solution, not a core patch" below for what was checked and what a future
native-side fix would look like.

## What ngspice already does (and does not) -- with evidence

ngspice does have real Monte Carlo primitives: `agauss`/`gauss`/`aunif`/
`unif`/`limit` usable inside a `.param` expression or directly in a device
value, `setseed`, `.option seed=<n>|random`, and `mc_source` (re-source and
re-parse the circuit from ngspice's own internal copy, re-evaluating every
random function) for use inside a `.control dowhile` loop. This project's
own copy of the upstream examples shows the intended idiom directly:
`ngspice-src/examples/Monte_Carlo/MC_2_control.sp`, `MonteCarlo.sp`,
`ngspice-src/examples/various/agauss_test.cir`,
`ngspice-src/examples/Monte_Carlo/rand_numb_test.cir` (manual ch. 17.8.6 /
18 / 22, `extras.html`).

**Finding 1 (the one this tool is built around): a `.param` containing a
random function is not bound once -- it is re-substituted and re-evaluated
at every point of use.** This was not assumed from the manual; it was
probed directly on `/usr/local/bin/ngspice` (45.2, the pinned system
build), first suspected from research/feature-survey.md item #7's flag
("no distinct local-mismatch primitive separate from a plain
global-parameter MC"), then confirmed to be worse than that flag states.

Single parse, no `.control` loop, no `mc_source` at all:

```
.param proc_off = agauss(0, 10m, 3)
.param a = proc_off
.param b = proc_off
R1 1 0 R='1k*(1+proc_off)'
R2 2 0 R='1k*(1+proc_off)'
R5 5 0 R='1k*(1+a)'
R6 6 0 R='1k*(1+b)'
...
.control
op
print @r1[resistance] @r2[resistance] @r5[resistance] @r6[resistance]
.endc
```
```
@r1[resistance] = 1.010073e+03
@r2[resistance] = 1.000933e+03
@r5[resistance] = 1.003192e+03
@r6[resistance] = 9.984058e+02
```

R1 and R2 both reference the *same* `.param proc_off` by name -- if
`.param` bound a value once at declaration (what every other SPICE dialect
does, and what the manual's phrasing implies), they would read back
identical. They do not. Indirecting through `a`/`b` (plain copies, `agauss`
typed only once, in `proc_off`'s own declaration) does not fix it either --
`a` and `b` still diverge from each other. This means **the natural,
documented idiom for "shared/global process variation" -- one `.param`,
referenced by several devices -- silently gives independent, not
correlated, draws.** A 20-run loop with `setseed $i` / `mc_source` /
`op` per iteration (`/tmp/.../probe1.cir` during this card's own
investigation, not kept in the repo per the no-simulation-output rule)
showed the same divergence run after run -- this is not a one-off parser
glitch, it is numparam's actual substitution model for any `.param` whose
RHS contains a call with side effects (a fresh random draw).

**Finding 2 (independently documented by ngspice's own shipped example,
not found by this card):** `ngspice-src/examples/Monte_Carlo/MC_2_control.sp`
contains this comment about looping `tran` inside one process via
`mc_source`:
> `* without the reset switch there is some strange drift`
> `* towards lower and lower frequencies`

i.e. even ngspice's own example distrusts the in-process reload path
enough to work around it, without fully explaining the mechanism.

**What this means for "is agauss in a `.param` re-drawn per instance or
once per run":** neither, cleanly -- it is re-drawn at **every textual
reference**, whether or not that reference is inside a subcircuit, inside
another `.param`, or repeated verbatim. There is no ngspice mechanism that
gives "one shared draw, referenced twice" for a random `.param`. Per-
instance independent mismatch is achievable today (put the `agauss()` call
directly in each instance's own value/parameter, or inside a `.subckt`
body so each `X` instantiation gets its own local scope) -- that part
already works, is not new, and is not what this card had to fix. What does
not work, and is silently wrong rather than merely absent, is the *global*
half of "global vs local".

## Why a tools-level solution, not a core patch

The card's brief prefers a tools-level solution when the core already has
the needed primitives. Given Finding 1, a tools-level solution is not just
preferable but *necessary* for the global side to be correct at all --
patching numparam to memoize a `.param`'s random draw once at declaration
would be a real, surgical, plausible fix (S-M effort, matches research/
feature-survey.md item #7's own estimate), but it is a **core semantic
change to every existing deck that uses `.param ... agauss(...)`**, is out
of scope for this card once the tools layer is shown sufficient, and was
therefore not attempted -- left as a candidate for its own future card, not
created here per this session's instructions.

`tools/mc.py` instead:
- does **all** sampling in Python (stdlib `random`, no numpy/scipy on this
  box -- confirmed absent), with a deterministic, hash-derived seed per
  `(master_seed, run_index, param_name)` (`mc._seed_int` /
  `mc._rng_for`) -- independent of any other run's or parameter's state,
  so replaying run *k* alone needs nothing but the spec, the seed, and
  *k*;
- substitutes the sampled values into the deck as **plain numeric
  literals** via `@@NAME@@` placeholders (never `.param ... agauss(...)`),
  so no ngspice numparam behaviour -- correct or broken -- is depended on
  for the actual randomness;
- runs each Monte Carlo point as its own ngspice **process**
  (`subprocess.run(["ngspice", "-b", run_k.cir])`), never the in-process
  `.control dowhile`/`mc_source` loop, sidestepping Finding 2 as well as
  Finding 1;
- makes "global" vs "local" purely a naming convention in the spec: the
  same `@@NAME@@` token used in more than one device's expression is
  global (one shared value per run); a distinct token per instance is
  local/independent (optionally Pelgrom-scaled). This is deliberately the
  same convention sky130's public PDK documentation describes for its own
  `mc_pr_switch` (process, shared) / `mc_mm_switch` (mismatch, per-
  instance) idiom -- not copied from any private deck, read from sky130's own public PDK
  documentation, and re-implemented here in a way that does not hit
  Finding 1.

## `tools/mc.py`

```
tools/mc.py run SPEC.json [--out DIR]
tools/mc.py replay MANIFEST.json --run K --spec SPEC.json
```

A spec (JSON file, or a plain dict when used as a library -- see
`tests/test_mc.py` for real examples) has:

- `template` / `template_path`: deck text with `@@NAME@@` placeholders.
  `@@...@@` was chosen because it collides with nothing ngspice itself
  emits or parses (unlike `{...}` or `$...`, both real ngspice syntax).
- `seed` (int), `n` (point count).
- `params`: list of `{name, scope, dist, nominal, sigma}` or
  `{..., pelgrom: {coeff, w, l}}` (sigma = coeff / sqrt(w*l)). `dist` is
  `gauss` (default, `sigma` = 1-sigma) or `unif` (`sigma` = half-range).
  `scope` is `"global"`/`"local"`/anything else -- purely documentation/
  bookkeeping, see above; the actual global/local behaviour comes from
  whether the `name` is reused across the template.
- `measures`: names the deck's own `.control` block must report via a
  line of the exact form `MCMEAS <name>=<value>` on stdout, e.g.:
  ```
  .control
  op
  let x = @r1[resistance]
  set xs = $&x
  echo MCMEAS r1=$xs
  .endc
  ```
  (this satisfies "read back per run via .meas of the device values or
  currents" without fighting `.meas`'s awkward pairing with a bare `.op`,
  and without a rawfile parser).
- `checks` (optional): `{"type": "correlation"|"mean"|"sigma", ...,
  "expect": v}`, evaluated against the 99% CI of the resulting statistic
  -- `pass` if `expect` is inside the CI, else `FAIL`.
- `max_workers` (hard-capped at 4 regardless of what the spec asks for --
  this card's CPU-modest rule), `timeout` (seconds/process), `ngspice`
  (binary path, default `$NGSPICE` or `ngspice`).

Every run's deck, log, sampled parameter values (with enough information
to re-derive them, since they are also purely a function of `seed`/
`run_index`/`param name`) and measured results are written to
`manifest.json` under the output directory (never the repository --
`~/Work/data/ngspice_augmentation/mc/...` by CLI default). `mc.replay_run`
re-*samples* (not re-reads) one run's parameters from the spec alone, reruns
it, and diffs against the manifest's recorded result.

### Statistics (stdlib only)

- mean, sample stdev: `statistics.mean`/`statistics.stdev` (ddof=1).
- 99% CI for a mean: normal approximation (`statistics.NormalDist`,
  z=2.5758). At N~=500 the Student-t correction (df=499) changes z by
  <0.4% -- stated here, not hidden.
- 99% CI for a sigma: Wilson-Hilferty cube-root normal approximation to
  the chi-square distribution -- accurate to a fraction of a percent in
  the N~=500 regime this tool uses it in; not claimed accurate at small N.
- 99% CI for a Pearson correlation: Fisher z-transform, standard normal CI
  in z-space, inverse-transformed back.

## Validation (`tests/test_mc.py`)

Matched **resistor** pair (the spec's own "matched MOS or resistor pair"
option), read back via `@r1[resistance]`/`@r2[resistance]` (device-value
readback, per the spec's own "e.g." wording). All statistical checks at
N=500. Real captured output, this card, `NGSPICE=/usr/local/bin/ngspice`
(the pinned unmodified system build):

```
[pass           ] ngspice availability smoke test: 2/2 runs ok
[pass           ] global-only: correlation r1: r=1.0000 ci99=[1.0, 1.0] expect=1.0 n=500
[pass           ] global-only: mean r1: mean=1000.18 ci99=[997.86, 1002.50] expect=1000.0
[pass           ] global-only: sigma r1: sigma=20.1351 ci99=[18.61, 21.91] expect=20.0
[pass           ] global-only: correlation exactly 1 (identical formula, not just correlated): r=1.0
[pass           ] local-only: correlation r1: r=0.0250 ci99=[-0.090, 0.140] expect=0.0 n=500
[pass           ] local-only: sigma r1: sigma=9.964 ci99=[9.209, 10.843] expect=10.0
[pass           ] local-only: sigma r2: sigma=9.812 ci99=[9.069, 10.677] expect=10.0
[pass           ] both: correlation r1 (analytic expect=0.8000): r=0.8061 ci99=[0.762, 0.843] expect=0.8 n=500
[pass           ] both: sigma r1 (analytic expect=22.3607): sigma=21.977 ci99=[20.313, 23.916] expect=22.361
[pass           ] Pelgrom scaling: sigma_small=21.167 sigma_big=5.292 ratio=4.000 expect=4.000 rel_err=0.0%
[pass           ] reproducibility: 50 runs compared, 0 mismatches
[pass           ] reproducibility: replaying run 17 alone reproduces its recorded result: matched
[pass           ] red-proof: checker state=FAIL (expected FAIL) r=1.0 detail=...expect=0.0 n=500
[pass           ] worked example (non-gating): W=8u/L=1u rel_err=3.8%; W=32u/L=4u rel_err=5.3%
15 pass, 0 FAIL (gating), 0 COULD-NOT-RUN (gating), 0 info/non-gating
```
(`python3 tests/test_mc.py`, ~3 seconds wall clock, ~3000 ngspice processes
total across all checks, 4 concurrent max.)

The **"both" mixture check is the one that most directly answers the spec
criterion**: with global sigma=2% and local sigma=1% (per instance), the
analytic prediction for the correlation between two devices sharing the
global term is `sigma_g^2/(sigma_g^2+sigma_l^2) = 0.8` -- measured 0.806,
inside its 99% CI.

The **red-proof** constructs a deck where a "local" (mismatch) parameter is
deliberately wired to the *same* token for both devices (the exact bug
Finding 1 causes if you rely on ngspice's own `.param` sharing) and checks
that `mc.py`'s own independence assertion (`expect correlation=0.0`)
reports `FAIL`, not that the tool always says pass -- proving the checker
actually distinguishes a genuinely independent local draw from an
accidentally-shared one.

## Worked example: differential-pair offset, MC vs analytic Pelgrom

`docs/examples/mc_diffpair_offset.cir`: a textbook 5T-style NMOS
differential pair (`nch4`, the same synthetic BSIM4-family model as
`bench/gen_analog_block.py` uses -- synthetic, not derived from any real design), resistive loads,
100uA tail current. A local/mismatch input-referred offset
(`@@VOFF@@`, series voltage source in M2's gate path) is Pelgrom-scaled:
`sigma_Voff = A / sqrt(W*L)`, illustrative `A = 4 mV*um` (a made-up,
clearly-labelled-as-synthetic number for this synthetic model -- not a
real process's characterised mismatch coefficient).

**How the gain figure was measured** (not assumed): a deterministic,
non-random finite-difference run (`op` at `Voff=0`, `alter voff = 1m`,
`op` again) gives `Vout_diff(1mV) - Vout_diff(0) = 3.78598 mV`, i.e. a
measured gain of **3.786 V/V** from input-referred offset to output
differential voltage, at this bias point (`W=8u/L=1u`, tail 100uA). This
empirical gain (not a hand-derived BSIM4 formula, deliberately, per this
project's "never report a result you did not run") is what converts the
MC's *output*-referred sigma back to an input-referred number for
comparison against the Pelgrom prediction.

Two device sizes, N=500 each, same campaign machinery as the resistor
tests above:

| W, L | W*L | predicted `sigma_Voff` (Pelgrom) | measured MC `sigma(Vout_diff)` | back-converted (`/3.786`) | agreement |
|---|---|---|---|---|---|
| 8u, 1u | 8 | 1.4142 mV | 5.2950 mV (99% CI [4.894, 5.762] mV) | 1.3986 mV | 1.1% |
| 32u, 4u | 128 | 0.3536 mV | 1.3429 mV (99% CI [1.241, 1.461] mV) | 0.3547 mV | 0.3% |

Ratio of the two back-converted sigmas: 1.3986/0.3547 = 3.94, vs the
Pelgrom prediction `sqrt(128/8) = 4.00` (1.4% off). This is the same
Pelgrom-scaling check as the resistor test, but propagated through a real
nonlinear BSIM4 device op-point solve rather than a linear resistor
formula, and cross-checked against an independently, empirically measured
gain rather than the same number used to inject the mismatch -- i.e. it is
not circular. (`tests/test_mc.py`'s check 7 reruns this, non-gating, at
N=300 with a looser 15% tolerance for everyday regression use; the numbers
above are the N=500 run captured for this document.)

## Limits / what was not done

- **No C/numparam change.** Finding 1 (broken shared-`.param` random
  re-evaluation) is real and reproducible but was not patched -- the
  tools-level workaround is sufficient for this card's spec criterion, and
  a core fix touches every deck using `.param ... agauss(...)` today,
  which is a bigger, separate decision. Left as a documented finding for a
  future card, not created here (out of scope for this session).
- **Statistics use approximations**, all stated above: normal (not
  Student-t) CI for the mean, Wilson-Hilferty chi-square approximation for
  the sigma CI, Fisher z for the correlation CI. All are standard and
  accurate at the N~=500 this tool's validation uses; none were checked at
  small N (say N<30), and the tool does not warn if used there.
- **Pelgrom check is ratio/tolerance-based**, not a formally derived joint
  confidence interval on the ratio of two sigmas -- a 20% (resistor) / 15%
  (MOS worked example, test suite) band, stated as a deliberately generous
  margin, not tuned to the observed numbers.
- **`unif` distribution's `sigma` field means half-range**, `gauss`'s means
  1-sigma -- different units for the same field name, documented here and
  in the module docstring, easy to trip over if not read.
- **No importance-sampling / high-sigma acceleration** -- this is plain
  brute-force MC (research/feature-survey.md item #6, "fast/high-sigma
  MC," is a separate, unbuilt item; this card is item #7 only, "native
  MC/mismatch primitive review").
- **`MCMEAS` is a convention this tool defines**, not an ngspice feature --
  every deck used with `mc.py` must add the `let`/`set`/`echo MCMEAS`
  lines itself; there is no automatic extraction from `.meas` or a
  rawfile.
- Every number in this document was produced by an actual run captured
  during this card (`tests/test_mc.py`'s own printed output, and the
  ad-hoc driver script used for the worked-example table above, both run
  on `/usr/local/bin/ngspice` 45.2, the pinned unmodified system build) --
  none were estimated or assumed.

## Draw-once form (card #482, private fork only)

This card's own "no C/numparam change" limit above was for card #473's
scope (a tools-level workaround was sufficient there). Card #482 picks up
the "future card" it deferred to: a small, private-fork **core** patch
(`patches-wip/numparam/0300`, branch `aug-numparam`) that adds an explicit
draw-once form to numparam itself, so a plain, in-deck `.param` (no
`tools/mc.py`, no `@@NAME@@` template layer) can express "one shared draw
per circuit instantiation" directly. The default per-use re-evaluation
(Finding 1 above) is **unchanged** -- this is a strict addition, opt-in per
`.param`.

### Mechanics found (reading `src/frontend/numparam/*.c` and
`src/frontend/inpcom.c`, confirmed by tracing actual calls to
`agauss()`/`nupa_assignment()`/`fetchnumentry()`/`nupa_define()` against
the Finding-1 reproducer, not inferred from the manual)

Finding 1's root cause is **not** inside numparam's expression evaluator
(`xpressn.c`) at all, and it is not "a `.param` re-reads its own random
call on every symbol lookup" -- `fetchnumentry()` (`xpressn.c`) is a plain
cached-value read (`entry->vl`, written once by `nupa_define()`); a bare
identifier reference is always a cache hit, never a re-draw. The actual
mechanism is a **deliberate, textual, HSPICE-compatibility preprocessing
pass**, `inp_fix_agauss_in_param()` (`src/frontend/inpcom.c`, called once
per statistical function name -- `agauss`, `gauss`, `aunif`, `unif`,
`limit` -- from `inp_readall()`, well before `inp_subcktexpand()` /
numparam ever runs):

1. Any `.param name = <expr>` whose RHS textually contains one of those
   five function names is rewritten to `.func name() <expr>` (a callable,
   zero-argument numparam macro), tracked with its subcircuit-nesting
   depth/count so a subckt-local `.param` only affects references inside
   that subckt.
2. Every *other*, bare textual occurrence of `name` elsewhere in the deck
   (in another `.param`'s RHS, in a `{...}`/`'...'` device value, anywhere)
   is rewritten to `name()`.
3. `inp_expand_macros_in_deck()` (called right after, still before
   `inp_subcktexpand()`) then **textually inlines** each `name()` call
   with the `.func`'s own body -- a real, textual, preprocessor-style
   macro expansion, not a runtime function call. This is why a `.param a =
   proc_off` (RHS = a bare copy, no `agauss` of its own) still ends up
   independently redrawn: by the time numparam evaluates that line its RHS
   already textually reads `(agauss(0,10m,3))`, freshly inlined -- a
   different, brand-new copy of the call text from R1's or R2's own
   inlined copy. Confirmed directly: instrumenting `nupa_eval()` shows
   the deck, by the time numparam's own single evaluation pass runs, has
   already been rewritten line-for-line (e.g. `.param a=proc_off` arrives
   at numparam already reading `.param a={(agauss(0,10m,3))}`), and
   `agauss()` fires once per such inlined copy, each producing a different
   draw -- exactly Finding 1's 1010.07/1000.93/1003.19/998.41 (`R1/R2/R5/R6`)
   divergence, all from ONE textual `.param proc_off = agauss(...)`
   declaration.
4. `reset` / `mc_source` / a fresh `source` all re-run this whole pipeline
   from scratch (`nupa_signal(NUPADECKCOPY)` -> `nupa_init()`, a brand new,
   empty symbol dictionary, once per deck-copy cycle) -- which is exactly
   why a fresh draw per reset "just works" for both the default and the
   draw-once form below, with no extra reset-hook code needed.

### Syntax chosen, and the three options considered

1. **A dot-card qualifier, `.param global p = agauss(...)`** (the card's
   own second candidate). Rejected: `.global` is already ngspice's own
   dot-card for global *node* names, and this project's own `tools/mc.py`
   vocabulary already uses "global" for a *shared* Monte Carlo parameter
   in the spec-file sense (see above) -- reusing the same word for a
   third, syntactically different thing (a `.param` qualifier keyword) was
   judged more likely to confuse a reader than to align with that
   vocabulary. It is also a slightly bigger patch: it needs a new keyword
   scan and strip step ahead of `inp_fix_agauss_in_param()`, rather than
   reusing infrastructure already in that function.
2. **A wrapper function, `once(expr)`, with its own runtime memoization**
   (cache keyed on source line number or symbol name, cleared on
   `nupa_signal(NUPADECKCOPY)`). Considered because it seemed like the
   "natural" place to put a freeze given `once()` reads as an expression-
   level operation. Rejected once the real mechanism (above) was
   understood: the freeze this card actually needs is "don't textually
   inline/redraw at every reference," which happens entirely in
   `inp_fix_agauss_in_param()`'s preprocessing, before any runtime
   evaluation -- a `once()` with its own cache would be solving the wrong
   layer (and would silently fail to help if used at a *reference* site
   rather than the `.param` declaration, since by then the identifier has
   already been macro-inlined away -- see the implementation note below).
3. **`once(expr)` as a pure syntactic marker, no runtime memoization**
   (chosen). `once()` is registered as an ordinary, trivial identity
   numparam function (`XFU_ONCE`, `mathfunction()`, `xpressn.c` --
   `y = x`, nothing else). The entire freeze effect is
   `inp_fix_agauss_in_param()` skipping its `.func`-rewrite (step 1 above)
   for any `.param` line containing the literal token `once(`
   (`search_identifier()`, the same word-boundary-aware helper the
   function already uses to find `agauss` etc.) -- so a draw-once `.param`
   simply never enters the textual-macro pipeline at all, and is
   evaluated exactly the way an ordinary constant `.param k = 5` already
   is: once, by `nupa_assignment()`, cached in the symbol table, read by
   every reference via the ordinary (never-broken) `fetchnumentry()` path.
   Chosen because it is the smallest change (one `continue;` plus one
   identity function), touches nothing at runtime (no new cache, no new
   state to reset), and composes for free with everything numparam
   already gets right for a plain `.param` (subckt scope, `reset`/
   `mc_source`/`alterparam`+`reset`, `.control` `dowhile` loops).

**Usage note (a real caveat, not a hypothetical one):** `once(...)` must
wrap the **declaration**'s RHS (`.param p = once(agauss(...))`), not a
reference site. `R1 ... once(p) ...` where `p` was declared as a *plain*
`.param p = agauss(...)` does nothing useful -- `p` has already been
turned into a `.func`/textually inlined by the time `once()` would run,
so each reference site's `once()` call only freezes *that one site's own*
already-independent draw, not a value shared across sites.

### Implementation

`patches-wip/numparam/0300-numparam-draw-once-form.patch` (branch
`aug-numparam`, on top of `aug`), two files:
- `src/frontend/numparam/xpressn.c`: adds `once` to `fmathS` and a new
  `XFU_ONCE` enum value/`mathfunction()` case (identity).
- `src/frontend/inpcom.c`: `inp_fix_agauss_in_param()` gains a
  `search_identifier(curr_line, "once", curr_line)` check that `continue`s
  (skips the `.func` rewrite) before doing anything else with the line.

### Tests (`tests/test_numparam.py`, fixtures `tests/fixtures/np_*.cir`)

Built in a separate fork worktree (`ngspice-src-numparam`, branch
`aug-numparam`), prefix
`~/Work/apps/ngspice_augmentation/aaba50ba3f49-np1`. A second prefix,
`~/Work/apps/ngspice_augmentation/aaba50ba3f49-np1-nopatch`, built from the
same worktree with the two patched files reverted to their pre-#482
content (verified before AND after each build, via `git status`/`grep`
before compiling and `strings <object>.o | grep once` after -- a first
attempt at this got the revert/rebuild ordering wrong and silently linked
a stale patched object into the "reverted" binary; caught by exactly this
`strings` check, not assumed), used only for the red-proof check.

Actual captured output, this card:

```
pass          shared-once single-run: r1==r2==r5==r6                                          values=[1001.42, 1001.42, 1001.42, 1001.42]
pass          shared-once reset-loop: r1==r2 in every one of 400 iterations                   400 iterations, 0 mismatched pair(s)
pass          shared-once reset-loop: mean/sigma of R match agauss(0,10m,3) within 99% CI     n=400 mean=1000.1147 CI=[999.6912,1000.5382] (expect 1000.0000); sigma=3.2881 CI=[3.0120,3.6157] (expect 3.3333)
pass          regression: plain .param agauss(...) still gives independent per-use draws      values=[1002.735, 1001.145, 999.4908, 997.4097]
pass          regression: sky130-style per-instance mismatch (subckt-local agauss) unchanged  x1=1000.954 x2=995.1721
pass          red-proof: patch reverted -> shared-once guarantee does not hold                reverted build did not produce the expected vectors (...) -- guarantee does not hold, as expected

6/6 pass
```

The red-proof's reverted build does not merely give different values, it
fails to run the deck at all (`Undefined parameter [once]` /
`Formula() error`, exit 1) -- a stronger demonstration that the shared-
value guarantee depends entirely on the patch, not merely a check that
happens to still pass.

Full regression, same `aaba50ba3f49-np1` build: `tests/test_stb.py`
60/60, `tests/test_stb_native.py` 28/29 (1 skip, pre-existing/unrelated),
`tests/test_stb_native2.py` 54/54, `tests/test_stb_device.py` 16/16,
`tests/test_stb_followups.py` 29/29, `tests/test_stb_cli.py` 42/42,
`tests/test_compopt.py` 20/20, `tests/test_mc.py` 15/15 -- all unchanged
from the parent's own baseline numbers, confirming the patch (2 files,
~30 lines) does not disturb anything outside the `once(` opt-in path.

### What was not done

- **Promoted 2026-09-29**: now `patches/0016-numparam-draw-once-form.patch`
  in the series and present in every build from s16 on (the `patches-wip/
  numparam/0300` name below is its original working name).
- **No `.func`-level general fix** -- a `.func` whose body calls a random
  function has the identical per-call-site redraw behaviour as a bare
  `.param`, unaddressed by this card (out of scope: the card's problem
  statement is specifically about `.param`).
- **No attempt to make `once()` useful at a reference site** (see the
  usage note above) -- would need real runtime memoization keyed off the
  *declaration* identity, which the chosen design deliberately avoids
  needing; left as a known limitation, not a bug.
- **Not validated against `alter`/`alterparam` directly** -- only
  `reset`/`mc_source`, per the card's own two named mechanisms tested
  end-to-end; `alterparam ... ; reset` (`ngspice-src/examples/various/
  alterparam.sp`'s own idiom) is not itself a distinct code path from
  `reset` (it only stages a value that `reset` then picks up), so this is
  believed covered by the `reset` testing above, but was not separately
  exercised with an `alterparam` line in the test suite.
- **`tests/test_numparam.py`'s red-proof needs a second, reverted build**
  (`NGSPICE_NUMPARAM_NOPATCH`) to actually run -- without it the check
  reports `COULD-NOT-RUN`, honestly, rather than a silent skip; that
  reverted prefix is a local build artifact (`~/Work/apps/
  ngspice_augmentation/aaba50ba3f49-np1-nopatch`), not committed anywhere,
  and will need rebuilding if this is re-verified later.
