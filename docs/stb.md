# Loop stability (`stb`) in ngspice_augmentation

Written for analog designers who know Cadence Spectre's `stb` but not this
project's code. Card #463 (plan T7.1), 2026-09-28; updated by card #464,
2026-09-28 (R1c positive-feedback fix, the `low_f_positive_feedback`
warning, and falling/rising crossing labelling -- see "Interpreting
margins" and example (c) below for what changed and why); updated again by
card #466, 2026-09-29, fixing an independent reader's findings
(`~/Work/data/ngspice_augmentation/verify-v-g4/`, 2026-09-29): a real code
bug where native `.stb`/`stb` silently defaulted an explicitly-bad
`fstart`/`fstop`/numsteps instead of refusing (see "Frequency and
point-count validation" below), netlist-card `probe=`/`probe2=`/`mode=`/
`localgnd=` now auto-quote themselves (no longer *require* manual quoting
-- see "Netlist-card quoting"), and several corrections to this doc's own
examples and numbers (magnitude narrative, crossing labels, native-vs-
script agreement, `patches/` status, terminal-name reference, the
`.stb`-in-batch-mode-needs-`.print` note). **All example output in this
doc was re-captured on card #466's build**,
`~/Work/apps/ngspice_augmentation/aaba50ba3f49-stb5/bin/ngspice-aug` (pin
+ `patches/0001`-`0011` + `patches-wip/stb5/0012` -- see "How this was
validated" for the full regression table); updated again by card #476,
2026-09-29 (private fork, WORKSPACE `ngspice-src-stb6`/branch `aug-stb6`),
**narrowing** the `low_f_positive_feedback` warning card #464 added: it
false-triggered on healthy integrator-type loops (a type-III compensated
buck, `tests/fixtures/compopt_r3_averaged_detuned.cir`, at `fstart=10 Hz`
-- `Re(T)` genuinely straddles 0 a few degrees past -90 deg on a
near-integrator branch, not a real instability). The warning now fires
only when `|arg T(fmin)| > 150 deg` **and** `|T(fmin)| > 1` (positive
feedback WITH loop gain); a neutral `low_f_integrator` label appears when
`arg T(fmin)` is within -90 +/- 30 deg; and a note (not the warning) when
`|arg T(fmin)| > 150 deg` but `|T(fmin)| <= 1` -- see "Interpreting
margins" below for the corrected explanation and real captured output on
this card's build,
`~/Work/apps/ngspice_augmentation/aaba50ba3f49-stb6/bin/ngspice-aug` (pin +
`patches/0001`-`0012` + `patches-wip/stb6/0013`).

There are **two separate tools**, not one:

1. **`tools/stb.py`** -- a testbench + control script that runs on the
   **unmodified system ngspice** (`/usr/local/bin/ngspice`, 45.2). It
   inserts a probe subcircuit (`stb/lgprobe.lib`, or `stb/diffprobe.lib`
   for differential/common-mode), runs two ordinary `ac` analyses, and
   applies Tian's two-injection formula in Python (`tools/stbcore.py`).
   Nothing about this needs a custom build -- if you can run `ngspice -b`,
   you can run this.
2. **Native `.stb`** -- a real analysis type (`STBinfo`, `.stb` netlist
   card or interactive `stb` command) built into this project's **private
   fork** of ngspice, at `~/Work/apps/ngspice_augmentation/<prefix>/bin/ngspice-aug`.
   It does the same Tian two-injection math inside the simulator (one
   matrix factorization and two back-solves per frequency, in
   `src/spicelib/analysis/stban.c`), plus a **device-based ("gain-nulling")
   mode** that the script-based tool cannot do at all: probing a MOSFET
   directly to see local loops (through Cgd, for example) that no wire
   break can reach.

Both compute the **same T** (Tian, Visvanathan, Hantgan, Kundert, IEEE
Circuits & Devices 17(1), 2001, eq. 30) and report it in the **same
Spectre sign convention** (`loopGain = -T`). Where they overlap, they have
been checked against each other to roughly 1e-9 to 1e-14 relative error
(see "How this was validated" below) -- for a **loop-mode** probe placed
at the same wire, treat them as interchangeable. Device mode and native
DM/CM (`mode=dm|cm`) exist **only** in the native tool.

Background theory (Tian's method, Spectre's `stb`/`pstb`, the math behind
every formula below): `research/cadence-stb-pstb.md`.

## Which tool to use when

| Situation | Use |
|---|---|
| You have `/usr/local/bin/ngspice` only, no private-fork build | `tools/stb.py` |
| Ordinary loop-gain probe on a wire (op-amp feedback, LDO feedback, ...) | Either -- native is faster (see "How this was validated") and needs no probe subcircuit in the netlist |
| Differential / common-mode loop gain | Either -- `tools/stb.py --mode dm\|cm` with `stb/diffprobe.lib`, or native `probe=<Vp> probe2=<Vn> mode=dm\|cm` |
| A local loop **inside** a transistor (e.g. through Cgd) that no wire break can isolate | **Native device mode only** (`probe=<MOS1 or BSIM4 instance>`) -- `tools/stb.py` has no equivalent |
| You want to probe a terminal without hand-editing the netlist | **Native only** (`probe=<inst>:<terminal>` auto-splices a probe at parse time) |
| You need Spectre-style `values`/`valuesfile` sweeps, `.dc`-linked param sweeps, or `pstb` (loop gain of a switching converter's periodic steady state) | **Neither** -- not implemented (see "Known limits and caveats") |

## Quick start

### `tools/stb.py` (system ngspice)

Put an `lgprobe` instance in series with the wire that breaks the loop:

```spice
XP <a> <b> lgprobe
.include <path-to-repo>/stb/lgprobe.lib
```

`a` is **downstream** (the side the loop signal travels to), `b` is
**upstream** (the side it comes from) -- see "Probe placement" below.

```sh
python3 tools/stb.py DECK.cir --probe XP --sweep "dec 50 1 1e9" [--out DIR]
```

Prints phase/gain margins to stdout; writes `stb.csv` (T and Spectre-sign
`loopGain` per frequency) and `summary.json` under
`~/Work/data/ngspice_augmentation/stb/<deck>-<timestamp>/` (never the
repository) unless `--out` is given. Exit status 0 = ran, 2 = could not run.

For differential/common-mode, use `stb/diffprobe.lib` instead and add
`--mode dm|cm`:

```spice
XP <a1> <b1> <a2> <b2> diffprobe
.include <path-to-repo>/stb/diffprobe.lib
```
```sh
python3 tools/stb.py DECK.cir --probe XP --mode dm --sweep "dec 50 1 1e9"
```

### Native `.stb` (private fork build)

Build once (`buildrecipe/README.md` has the full recipe and rollback):

```sh
buildrecipe/build.sh --patches --prefix-name=<name>
```

installs `~/Work/apps/ngspice_augmentation/<name>/bin/ngspice-aug` (never
touches `/usr/local`). Point `NGSPICE` at it, or call it directly:

```sh
~/Work/apps/ngspice_augmentation/<name>/bin/ngspice-aug -b DECK.cir
```

No probe subcircuit needed -- the probe is an ordinary 0 V V-source (or an
existing device terminal, or a MOSFET instance for device mode) already in
the netlist. Interactive form, inside a `.control` block (no quoting
needed):

```spice
Vprobe <a> <b> DC 0
.control
stb dec 50 1 1e9 probe=Vprobe
.endc
```

or as a netlist card (quoting is optional as of card #466 -- `probe=`,
`probe2=`, `mode=`, and `localgnd=` all auto-quote themselves now; see
"Netlist-card quoting" below):

```spice
.stb dec 50 1 1e9 probe=Vprobe
.print stb loopgain
```

**A netlist-card `.stb` in batch mode (`-b`) needs a `.print`/`.plot`
line, or a `.control` block, or it prints "incomplete or empty netlist ...
no simulations run!" even though the analysis itself ran fine** -- the
same rule as `.ac`/`.tran`/every other batch-mode dot-card analysis, not
specific to `.stb`, but easy to trip over the first time (found while
re-running this doc's own examples for card #466). The `.print stb
loopgain` line above is enough; a `.control` block with `run`/`wrdata`
(what this project's own tests use) works too.

Both print the same phase/gain margin lines `STBan()` always prints, and
expose vectors `loopgain`, `t` (per frequency) and a second scalar plot
`phaseMargin`/`phaseMarginFreq`/`gainMargin`/`gainMarginFreq` (`NaN` when
there is no crossing), readable with `wrdata`/`print`/`plot` like any other
analysis's output.

## Syntax reference: native `.stb`

```
.stb dec|oct|lin <n> <fstart> <fstop> probe=<probespec> [probe2=<Vn>] [localgnd=<node>] [mode=dm|cm]
```

Same `dec`/`oct`/`lin`/`n`/`fstart`/`fstop` prefix as `.ac`, but **not**
the same "invalid values fall back to a default" behaviour `.ac` has --
see "Frequency and point-count validation" right below (card #466 fixed a
real bug here: the fork used to silently default a bad value with only a
generic warning, exactly like `.ac`, but the analysis's own E_PARMVAL
checks were supposed to refuse it and could not, because the value was
already overwritten before they ever ran). Only a field that is truly
**omitted** (no token there at all -- e.g. `.stb dec probe=Vprobe`, `n`/
`fstart`/`fstop` all left out) falls back to a default (`n`->10,
`fstart`->1, `fstop`->1000x fstart) with a printed warning ("Warning,
ngspice assumes default parameter(s) for stb simulation"); an *explicit*
bad value refuses instead (see below).

### Frequency and point-count validation

`fstart`, `fstop`, and `n` (numsteps), when given explicitly, are
validated and a bad one **refuses the whole run** (`E_PARMVAL`, no
analysis, no vectors -- verified in both the netlist-card and interactive
forms, card #466):

| Condition | Refused with |
|---|---|
| `fstart <= 0` | `Frequency of <= 0 is invalid for STB start` |
| `fstop <= 0` | `Frequency of <= 0 is invalid for STB stop` |
| `fstart > fstop` (both given) | `STB stop frequency must not be less than STB start frequency` |
| `n <= 0` | `Number of points must be > 0 for STB numsteps` |

A field that is **omitted** entirely (see above) still defaults, with the
generic warning, exactly as before -- this table is only about a value
that *is* present and out of range. `tests/test_stb_cli.py` is the
COULD-NOT-RUN-proof for this table: every row above asserts no `loopgain`
vector is produced, in both forms, not merely a nonzero exit code.

### `probe=` (required), three forms

| Form | Meaning | Mode |
|---|---|---|
| `probe=<Vname>` | an existing 0 V V-source already in the netlist, in series with the loop-break wire | loop |
| `probe=<inst>:<terminal>` | auto-insert: splits the named terminal of an existing instance onto a new internal node, in series with a freshly-created 0 V V-source, **at netlist-parse time** (before the instance's own `.stb` card is even reached -- the instance must already exist in the deck; a forward reference is a hard parse error, not a silent skip) | loop |
| `probe=<MOS1 or BSIM4 instance>` | device-based ("gain-nulling") mode: nulls the device's gm and reports the resulting return ratio -- no wire is broken at all | device |

`probe=<Vp>,<Vn>` (comma form) together with `mode=dm|cm` is accepted but
**not recommended**: the interactive `stb` command is tokenized by
ngspice's own command-line lexer *before* it is reassembled into a `.stb`
line, and that lexer also splits on comma, so `probe=Vp,Vn` can arrive
mangled. Use `probe=<Vp> probe2=<Vn> mode=dm|cm` instead (three ordinary
`key=value` pairs, nothing for any tokenizer to split on) -- this is the
form this project's own tests and this doc's examples use.

For `probe=<inst>:<terminal>`: `a` (downstream, where run 1's voltage is
read) is the new internal node now bound to the terminal; `b` (upstream,
where run 2's current is injected) is the original net the terminal used
to connect straight to. Idempotent (probing the same inst:terminal twice
reuses the already-inserted probe). A nonexistent instance or terminal
name is a hard parse error (never a silent default probe).

`<terminal>` must match one of that device type's own terminal names
(case-insensitive), read directly from its `IFdevice` table
(`stb_probe_splice()`, `src/spicelib/parser/inp2dot.c`, which works
generically for **any** device type that has one -- virtually every
device in the fork, not just the ones below). A few commonly-probed ones
(exact spelling from each device's own `*names[]` array in its
`src/spicelib/devices/<dev>/<dev>.c`):

| Device | Valid `<terminal>` names |
|---|---|
| V-source (`vsrc`) | `V+`, `V-` |
| Resistor (`res`) | `R+`, `R-` |
| Capacitor (`cap`) | `C+`, `C-` |
| Inductor (`ind`) | `L+`, `L-` |
| I-source (`isrc`) | `I+`, `I-` |
| MOS1/MOS2/MOS3 | `Drain`, `Gate`, `Source`, `Bulk` |
| BSIM4 | `Drain`, `Gate`, `Source`, `Bulk`, `Charge` (the NQS charge node) |
| BJT | `collector`, `base`, `emitter`, `substrate` |
| JFET | `Drain`, `Gate`, `Source` |

`tests/test_stb_device.py`'s `Rf:R+` (a resistor) and this doc's own
device-mode examples' `probe=M1` (a whole-instance device-mode probe, not
a terminal splice) both exercise this code path; not every device type in
the table above has a dedicated test, since the splice code itself does
not special-case any of them.

### Other parameters

| Parameter | Type | Meaning | Default |
|---|---|---|---|
| `probe2=<Vn>` | instance | branch-n probe for native DM/CM; requires `mode=` | none |
| `mode=dm\|cm` | keyword | differential (antiphase) or common-mode (in-phase) correlated two-injection with `probe2=`; requires `probe2=` | none (single-ended) |
| `localgnd=<node>` | node | reference the current-injection run to this node instead of true ground | true ground (no effect) |

### Errors and warnings (never a silent number)

- **`fstart <= 0`, `fstop <= 0`, `fstart > fstop`, or `n <= 0`** (card
  #466): refused (`E_PARMVAL`, no analysis, no vectors) -- see "Frequency
  and point-count validation" above for the exact message per case. Only
  a genuinely *omitted* field still defaults with the generic "assumes
  default parameter(s)" warning.
- Probe instance not found, or not a V-source/MOS1/BSIM4: `E_NOTFOUND`.
- `probe2=` given without `mode=`, or `mode=` given without `probe2=`:
  `E_PARMVAL`.
- `probe2=`/`mode=` combined with a **device-mode** probe: refused
  (`E_PARMVAL`) -- the two features cannot be combined.
- `localgnd=` given together with a device-mode probe: warning only, not
  fatal ("localgnd= has no effect in device mode -- ignored"); the run
  proceeds.
- **Either probe terminal is literal ground (node 0), loop mode**:
  warning, not refusal -- `T` reads out as **exactly 0 at every
  frequency**, honestly, because Tian's two-injection method is
  structurally degenerate there (see "Probe placement" below). Do not
  read `T=0` here as "provably stable, huge margin."
- **`|arg T(fmin)| > 150 deg` AND `|T(fmin)| > 1`** (card #464, narrowed by
  card #476 -- see "Interpreting margins" below for why): warning, not
  refusal -- a `WARNING: positive feedback at low frequency: ...` line is
  printed after the margins, and `lowFPositiveFeedback` (a vector on the
  `"STB Margins"` plot, alongside `phaseMargin`/etc.) is set to `1`.
  This is not a math error in the analysis; it means the *loop itself*
  has not settled to the positive-real DC value a negative-feedback loop
  should have, AND still has loop gain there (found on a real circuit bug,
  not a corner case of this tool -- see "Interpreting margins" below). The
  margins are still computed and printed, but must not be read as real
  stability margins when this fires. Two related, non-warning labels can
  also appear at the same point in the output: `low_f_integrator` (a
  neutral note, `arg T(fmin)` within -90 +/- 30 deg -- ordinary
  integrator-like behaviour, e.g. a type-III compensator's near-integrator
  branch) and `low_f_positive_feedback_low_gain` (a note, not a warning:
  `|arg T(fmin)| > 150 deg` but `|T(fmin)| <= 1`, so there is positive
  feedback but not enough loop gain there to call it unstable). Every 0 dB
  crossing is also now labelled falling (the usual stability-margin
  crossing) or rising (`phaseMarginFalling` vector; see "Interpreting
  margins").
- **BSIM4 device mode with NQS on** (`acnqsmod=1`, instance or model):
  refused outright (`E_NOTFOUND`) rather than silently nulling a
  real-only scalar where a complex, frequency-shaped one is needed.
- Probing a passive element or a nonexistent device in device mode:
  refused, never a number.

### Netlist-card quoting

On a **netlist** `.stb` card, ngspice's numparam evaluator runs over every
dot-card line before `.stb` itself is parsed, and tries to evaluate an
unquoted `key=value` RHS as a number -- `probe=Vprobe` would then fail
with `Undefined parameter [vprobe]` (this is not specific to `.stb`, or
even to a recognized card: an unrelated, made-up `.foocard bar=Vprobe`
line fails the identical way, confirmed while investigating this for card
#466). **As of card #466, this is fixed for `.stb`**: `dot_stb`
auto-quotes `probe=`, `probe2=`, `mode=`, and `localgnd=`'s value before
numparam ever sees the line (`inp_stb_autoquote()`,
`src/frontend/inpcom.c`), so all four now work **either quoted or
unquoted** on a netlist card, exactly like the **interactive** `stb`
control-block command (never run through numparam to begin with, so it
always accepted these unquoted). Quoting is still accepted -- harmless,
not required -- so existing decks/docs that already quote these need no
change. `tests/test_stb_cli.py` checks all four options both ways
(unquoted netlist-card form gives an identical `T` to the hand-quoted
form, for `probe=`, `probe=`+`probe2=`+`mode=` together, and `localgnd=`).

## Instant mode: linearising at a transient time point (`instant`)

Card #500, 2026-09-29. `.stb ... instant` (a bare flag, no value -- same
convention as `dec`/`oct`/`lin`) skips the internal DC operating-point
solve (`CKTop()`) that native `.stb`/`.ac` otherwise **always** run before
loading the small-signal matrix, and instead linearises about whatever
operating point is **already loaded in the circuit** -- node voltages
(`CKTrhsOld`) and every device's own persistent internal small-signal
state (transconductance, an incremental capacitance, ...), exactly as the
last thing that touched the circuit left them.

### Why this exists

`docs/stb-tran.md` (`tools/stb_tran.py`, card #497) found and proved,
three independent ways, that **without a fork change** native `.stb`/`.ac`
cannot be made to linearise about an arbitrary non-equilibrium point:
`CKTop()`'s Newton solve always converges to the circuit's TRUE DC
equilibrium given the present source values, discarding any prior
`.tran`/`IC=`/`.nodeset`/`savestate` state regardless of the starting
guess. `instant` is the fork change that closes that gap: it is the
`optime=<t>` / `.options acattran`-style hook card #500's own brief asked
for, implemented as the smallest correct change -- skip the ONE step
(`CKTop()`) that discards the state, not a new time-argument or a new
breakpoint-scheduling mechanism.

### How to use it: `.tran` + `.stb ... instant` as netlist cards, run together

```spice
.tran 50n 100u 0 uic
.stb dec 30 1 1e8 probe=Vprobe instant
.control
run
print stb2.phasemargin stb2.phasemarginfreq
.endc
.end
```

**Both must be top-level netlist (`.` prefixed) cards, run together via a
single `run`** (or an implicit batch-mode run) -- **not** two sequential
interactive `.control` commands (`tran ...` then `stb ... instant`).
This is not a style preference; it is required for `instant` to see the
transient's own state at all, and is worth understanding precisely
(found and verified empirically while validating this card, with a debug
dump of `ckt->CKTrhsOld`/`ckt->CKTtime` at `instant`'s own entry point --
`STB_DEBUG_INSTANT=1` reproduces it):

- **Every interactive analysis command** (`tran`, `ac`, `stb`, ... --
  `runcoms.c`'s `com_*` -> `dosim()` -> `if_run()`) calls
  `doAnalyses(ckt, RESTART, task)` (`spiceif.c`), i.e. `reset=1`.
  `CKTdoJob()` (`src/spicelib/analysis/cktdojob.c`), given `reset=1`,
  unconditionally sets `ckt->CKTtime = 0` and re-runs `CKTunsetup()` +
  `CKTsetup()` -- a full circuit reinitialisation -- **before running any
  job**, regardless of what analysis type is being started or whether
  that analysis itself calls `CKTop()`. So a `tran` command followed by a
  *separate* `stb ... instant` command still loses the transient's ending
  state: the reset happens one step before `STBan()` even starts, not
  inside it, so skipping `CKTop()` does not help.
- **Only `resume`** (`doAnalyses(ckt, RESUME, task)`, `reset=0`) skips
  that reset -- but `resume` continues the SAME paused job, it does not
  let a different analysis type pick up where a finished `tran` left off.
- **Jobs run together in one `doAnalyses(ckt, RESTART, task)` call** (one
  netlist, one `run`/batch invocation, several `.` analysis cards) hit the
  reset exactly ONCE, at the very start, before ANY job runs --
  `CKTdoJob()`'s per-job loop (`for (i = 0; i < ANALmaxnum; i++) for (job
  = task->jobs; ...) if (job->JOBtype == i) ...`) then runs each job in
  turn **without** an intervening reset. Crucially, this loop iterates by
  **analysis TYPE INDEX** (`analInfo[]`'s registration order in
  `src/spicelib/analysis/analysis.c`: `AC`, ..., `TRAN`, ..., `STB` --
  `STB` last, being the newest fork addition), **not netlist declaration
  order** -- so `.tran` always runs before `.stb` in the same task
  regardless of which card appears first in the deck. A `.tran` job
  ending at exactly `t_k` (its own `tstop`) therefore leaves `CKTrhsOld`
  holding the accepted transient solution at `t=t_k` when the following
  `.stb ... instant` job starts, with no reset in between.

Net effect: to linearise at `t_k`, set the `.tran` card's own `tstop` to
`t_k` (the stepper always lands exactly on `tstop`), and run it together
with `.stb ... instant` via one `run`/batch invocation, per t_k, one
process per t_k (`tools/stb_tran.py --mode=instant` does exactly this,
see `docs/stb-tran.md`).

### What "instant" computes, and why the CKTload(MODEINITSMSIG) call is still valid at a non-equilibrium point

`STBan()`'s own structure, unchanged either way: (default) `CKTop()` ->
`CKTload(MODEDCOP|MODEINITSMSIG)` -> per-frequency `CKTacLoad()` + solve.
`instant` removes only the `CKTop()` call; the `CKTload(MODEDCOP|
MODEINITSMSIG)` call right after it still runs, now evaluating the device
load routines at whatever node voltages are already present (the
transient's own `t_k` solution) instead of a freshly-converged DC
equilibrium. This is architecturally sound, not merely convenient,
because:

- **`CKTload` is one evaluation, not an iteration.** Nothing about it
  assumes the input voltages are a converged fixed point; it computes each
  device's stamp (conductances, an incremental capacitance, a
  transconductance) as a pure function of the terminal voltages currently
  in `CKTrhsOld` -- e.g. `asrcload.c`'s B-source handler reads
  `ckt->CKTrhsOld[...]` and differentiates its own expression tree fresh,
  every call, regardless of DC/transient/AC mode; `mos1load.c`/`bsim4...c`
  compute `gm`/`gds`/capacitances the same way, as algebraic functions of
  `Vgs`/`Vds`/`Vbs`. A device's incremental parameters are correct at
  `t_k` for the same reason they are correct at a DC equilibrium: they are
  evaluated AT the voltages present, not derived from HOW those voltages
  were reached.
- **The KCL residual this reload computes is not zero at a genuine
  transient point (a capacitor's real current is nonzero there, but this
  MODEDCOP-style reload treats reactive elements as open, folding them
  into the matrix only via `CKTacLoad`'s separate `jwC`/`jwL` stamp) --
  and this does not matter.** AC/STB analysis is a pure LINEARISATION: one
  complex-matrix factorisation and one or two back-solves for the probe's
  own excitation, never a Newton search that would need a converged
  residual to start from. A nonzero "DC-style" residual at the frozen
  point is simply never used -- `CKTacLoad`'s own RHS and the
  probe-injection RHS overwrite `CKTrhs`/`CKTirhs` completely before any
  solve happens.
- **Reactive elements are still exactly `jwC`/`jwL` in the small-signal
  matrix**, same as a DC-equilibrium `.stb`/`.ac` -- `instant` changes
  nothing about how `CKTacLoad` builds the frequency-dependent part of the
  matrix, only where the algebraic (real, DC-like) part's operating point
  comes from.

### Validation (`tests/test_stb_tran.py`, card #500's own done-when)

Both checked against `tests/fixtures/stbtran_ldo_softstart.cir` (the same
LDO soft-start fixture `docs/stb-tran.md` uses), this build
(`aaba50ba3f49-instant-N`):

1. **Settled point equals DC `.stb`.** At `t=1200us` (the ramp ends at
   `tss=500us`; the fixture's own header already shows `v(vout)` flat from
   `~600us` on), `instant`'s `phasemargin`/`phasemarginfreq` match a plain
   DC-op `.stb` of the unmodified deck (`Vss` at its final `0.45`) to
   **~1e-9 relative** (`19.55988778...` vs `19.55988776...` deg,
   `60744.02708` vs `60744.02725` Hz) -- three orders of magnitude inside
   the card's own 1e-6 bound.
2. **Mid-ramp point differs from quasistatic in the expected direction,
   and matches an independent reference.** At `t=100us` (only 1/5 of the
   way through the 500us ramp): `instant` gives PM=18.36 deg,
   fc=69243 Hz; `tools/stb_tran.py --mode=quasistatic` (frozen `Vss=0.09`,
   `docs/stb-tran.md`'s own table) gives PM=17.29 deg, fc=81379 Hz -- a
   real ~15% shift in crossover frequency, not noise. Mechanism: the
   fixture's nonlinear pass device (`Bpass`, `I =
   gm2*tanh(k2*(V(ea)-vth2))`) has an incremental gm that depends on
   `v(ea)`; quasistatic mode converges `ea` to the equilibrium a
   Vss=0.09-forever loop would settle to (`ea~0.045`, sech^2~0.98, near
   the tanh's linear region), while the REAL transient's `ea` at t=100us
   (measured directly, `meas tran ea find v(ea) at=100u`) is `~0.198`
   (sech^2~0.72) -- more than 4x larger, because the loop is still
   correcting a real tracking error the frozen-forever equilibrium never
   has. Lower incremental gm -> lower loop gain -> lower crossover, which
   is exactly the direction observed. Independently verified: measuring
   `v(ea)` at t=100us, hand-computing Bpass's EXACT incremental
   transconductance there (`gm2*k2*sech^2(k2*(ea-vth2))`, closed-form,
   Python), substituting a linear VCCS of that gain for the nonlinear
   `Bpass`, and running a plain (non-`instant`) DC-op `.stb` on the
   resulting fully-linear circuit gives PM=18.3645032... deg,
   fc=69243.107 Hz -- matching `instant`'s own PM=18.3645074... deg,
   fc=69243.069 Hz to **~1e-6 relative**, far inside the card's 0.5dB/3deg
   bound.

Both are in `tests/test_stb_tran.py`'s `check_instant_settled()`/
`check_instant_mid_ramp()` ("layer 4"), plus CLI red-proofs
(`check_instant_cli_red_proofs()`: `--mode=instant` without `--tstep`,
or combined with `--freeze`/`--check-tran`, exits 2/COULD-NOT-RUN, never
silently misinterpreted).

### Known limits

- **XSPICE event-driven instances are refused outright** (`E_PARMVAL`),
  not silently wrong: `instant` only skips the analog `CKTop()`; it does
  not skip `EVTop()`'s own event-side equilibrium re-solve (a separate,
  unaddressed instance of the same underlying problem), so an event-mixed
  netlist would otherwise silently combine a frozen analog state with a
  freshly re-initialised event state.
- **No literal syntax for "linearise at time t" on the `.stb` card
  itself** -- `instant` relies entirely on the CALLER (a `.tran` card with
  `tstop=t_k`, run together in the same `run`) to have already put the
  circuit at the desired instant; `.stb ... instant` on its own, with no
  preceding `.tran` in the same task, linearises about the circuit's
  power-up state (all-zero `CKTrhsOld`, `CKTtime=0`) -- not an error
  (there is nothing wrong with that state), but almost certainly not what
  is wanted; read the deck's own vectors/`meas` output to confirm the
  intended `t_k` was actually reached before trusting an `instant` result.
- **Only `.stb`, not `.ac`, was given this flag** -- `acan.c` has the
  identical `CKTop()`-then-`CKTload(MODEINITSMSIG)`-then-frequency-loop
  structure, and the same patch would apply there essentially unchanged
  (extending `.ac`'s own netlist-card parser, `dot_ac()`, to accept a
  trailing keyword tail at all is the only nontrivial part, since `.ac`
  currently takes only the bare `dec|oct|lin n fstart fstop` prefix and no
  keyword tail whatsoever) -- not implemented, time budget; the
  independent validation above used a plain (non-`instant`) `.ac`-style
  linear-circuit `.stb` instead, which does not need it.
- **No literal energy-storage-state freeze beyond what `.tran` itself
  already computed** -- this is not a limitation of `instant` specifically,
  it is the definition of the feature: `instant` reports whatever `.tran`
  actually converged to at `t_k` (correct, physical), it does not offer
  any other state to linearise about.
- Not validated on any circuit with XSPICE code models, KLU (`.option
  klu`), or a device family other than the ones the regression suite's
  fixtures already exercise (BSIM3/4, MOS1, diode, R/L/C, B-source) --
  the underlying claim (device loads are pure functions of terminal
  voltages, independent of analysis-type mode bits, for the ALGEBRAIC
  part of their stamp) is a general property of every `DEVload`
  implementation in this fork, not exercised device-by-device here.

## How results map to Spectre `stb`

| | This project | Spectre |
|---|---|---|
| Formula | Tian eq. 30 (both tools) | Tian eq. 30 (Tian, Visvanathan, Hantgan, Kundert 2001) |
| Sign | `loopGain = -T` (T = Bode return ratio, textbook sign) | `loopGain = -T` (same, per Wiedmann 2005 citing Tian) |
| Vector names | `loopgain`, `t` (native); `T`, `loopGain` columns (`tools/stb.py`'s `stb.csv`) | `loopGain`; `T` not separately exposed in ADE |
| Margin names | `phaseMargin`, `phaseMarginFreq`, `gainMargin`, `gainMarginFreq` | same names (ADE `getData`, per Wiedmann; not independently confirmed from a fetched Cadence doc, see `research/cadence-stb-pstb.md` sec 1) |
| Orientation | Tian T is symmetric under swapping the probe's two nodes (validated: "probe reversed = forward" passes on every loop-mode fixture) | same claim (Wiedmann) |
| `probe=` | `<Vname>` \| `<inst>:<terminal>` \| MOS1/BSIM4 instance for device mode | probe instance, or `probe=X:n` (terminal **index**, not name) |
| `localgnd=` | implemented, wired into the current-injection RHS | "Node name of local ground. If not specified, referenced to global ground." (same idea) |
| DM/CM | `mode=dm\|cm` + `probe2=` is a **fork-specific convenience parameter** -- Spectre has no `mode=` parameter; it achieves the same result by building a `diffprobe` subcircuit (two iprobes + a VCVS/CCCS pair) and altering its gain to +-1 before running plain `stb` | `diffprobe`/`diffstbprobe`/`cmdmprobe` subcircuit + gain alter |
| Device mode | explicitly documented: **gm only**, not gm+gmbs, for MOS1 and BSIM4 (quasi-static only) | "the dominant controlled source" -- which terms Spectre actually nulls is undocumented publicly (`research/cadence-stb-pstb.md` sec 2.5, `[DOC-gap]`) |
| Sweep types | `dec`/`oct`/`lin` frequency only | frequency (`dec`/`oct`/`lin`/`log`/`values`/`valuesfile`) plus `dev`/`mod`/`param`/`temp` sweeps at fixed `freq` |
| `pstb` (loop gain of a PSS periodic operating point) | **not implemented** -- native `.stb` here is small-signal about one ordinary DC operating point, same restriction plain Spectre `stb` has | separate analysis, `pstb`, needs `pss` first |

## Probe placement

Same rules as Spectre `stb` (`research/cadence-stb-pstb.md` sec 5):

- **Break every loop.** The probe wire must be a "critical wire" -- one
  whose removal opens every feedback path. If it only interrupts one loop
  among several, you get a necessary-but-not-sufficient per-loop check,
  not the whole circuit's stability.
- **No parallel path around the break.** If signal can get from the
  loop's input side to its output side any other way, Tian's break-point
  invariance stops holding.
- **Never place a probe terminal on literal ground.** Tian's two-injection
  method is structurally degenerate there: the probe's own branch equation
  `V(a) - V(b) = <this run's injection>` either pins the *other* terminal
  to the injection value outright, or the injected current has nowhere to
  go but the probe's own branch. Native `.stb` warns and reads `T = 0`
  (proved three ways during validation: algebraically, empirically with
  `STB_DEBUG=1`, and independently via `tools/stb.py`'s real 3-element
  `lgprobe` subckt, which comes out identically degenerate for the same
  reason).
- **Both sides of the break need real, finite admittance.** An idealized
  control input (a VCCS's control port, for example) can have *zero*
  admittance to the rest of the circuit -- breaking a loop right at one of
  those makes the current-injection run degenerate too. **This is not
  hypothetical**: building this doc's own example (a) originally broke the
  loop right at an op-amp-style VCCS's bare control input and hit exactly
  this (`ZeroDivisionError` in `tools/stbcore.middlebrook`, current
  injection run giving `i2=0`); adding a small input capacitance
  (`Cfb fb 0 2p`, modeling a realistic amplifier input capacitance) fixed
  it. If your loop-gain run blows up with a division-by-zero or a
  suspiciously flat `T`, check for exactly this.
- Loop-mode `T` does **not** depend on which physical node is called `a`
  vs `b` (validated on every fixture: "probe reversed = forward" passes to
  1e-9 to 1e-12 relative). Middlebrook's comparison-only output (available
  from `tools/stb.py`'s `run()`/`run_diff()` results, not printed by
  default) *does* care about orientation -- `a` must be downstream.

## Differential and common-mode

Both tools implement DM/CM as **two correlated single-ended probes**, one
per differential branch, driven in lockstep: antiphase (`acphase=180`) for
DM, in phase for CM. `T` is read from **branch 1 (`p`) only** -- this is
exact, not approximate, for a network that is *exactly* symmetric between
its `p`- and `n`-attached halves (the standard differential/common-mode
half-circuit theorem: a purely antisymmetric or purely symmetric
excitation makes the whole network's response equally antisymmetric or
symmetric, so the ratio `T` is the same read at either symmetric point).
It is **not** exact for an asymmetric network, and it does not couple the
two branches with an explicit mismatch-dependent source the way Spectre's
own `diffprobe` (an E/F pair) does -- see `stb/diffprobe.lib`'s header
comment and `patches/0006-*`'s commit message for the full derivation.

Both branches of the break need the **same downstream/upstream sense**
(matching `a`/`b` orientation on both wires) for "antiphase = DM,
in-phase = CM" to mean what it says.

`tools/stb.py --mode dm|cm` uses `stb/diffprobe.lib` (a 4-terminal
subcircuit: `Xname a1 b1 a2 b2 diffprobe`, deck-prepared into the
`diffprobe_dm`/`diffprobe_cm` variant by `tools.stb.select_diff_mode()`).
Native `.stb` uses two plain 0 V V-sources plus `probe=<Vp> probe2=<Vn>
mode=dm|cm` -- it does not expand the `diffprobe.lib` subcircuit directly
(the probe must resolve to an actual V-source instance, not a
subcircuit). See `docs/examples/c_diffamp_dmcm_diffprobe.cir` /
`c_diffamp_dmcm_native.cir` for the same circuit run both ways.

## Device mode

`probe=<MOS1 or BSIM4 instance>` nulls the device's **gm alone** (not
gm+gmbs) -- a rank-1 update to the already-factored AC matrix, one extra
back-solve, `T = 1/(1 - gm * v^T Y^-1 u) - 1` (research/cadence-stb-pstb.md
sec 2.5's "cheaper alternative" to a second LU factorization). Documented,
deliberate choice, not an oversight: gmbs is a real second controlled
source and **not folding it in is measurably not negligible** -- on a
source-degenerated stage, `tests/test_stb_device.py`'s body-effect check
found `|T(gmbs)|/|T(gm)|` up to (this run, 2026-09-28)
**6.4x** at some frequencies. If your circuit has meaningful body effect
at the probed device, read the gm-only device-mode number as one
component of the true loop gain, not the whole thing.

Node/sign convention (`stb_mos1_probe()`/`stb_bsim4_probe()` in
`stban.c`) was derived by hand against each device model's own AC stamp,
corner by corner, not merely assumed, and handles **reverse mode**
(`Vds < 0`, the model's own notion of drain/source swapped) automatically
-- validated in `tests/test_stb_followups.py`'s `check_reverse_mos`
(matches an independent brute-force reference using the *physical*, not
netlist-declared, drain/source roles; using the netlist-declared roles
literally gives errors of order 1-600, confirming the swap is load-bearing,
not a no-op).

BSIM4 device mode is **quasi-static only** in this fork -- `acnqsmod=1` is
refused outright, not silently computed with the wrong (real-only) scalar.
Validated to <=1e-6 relative up to 100 kHz (`tests/test_stb_followups.py`'s
`check_bsim4_cgd`, sweep `dec 15 1 1e5`); the same script's own comment
notes the residual (non-reciprocal Cdg transcapacitance the brute-force
reference does not fully model) crosses 1e-6 relative around 2e5 Hz, so
**do not extrapolate BSIM4 device-mode accuracy much past 100 kHz**
without re-checking.

`tools/stb.py` has **no device-mode equivalent** -- there is no wire to
put a probe subcircuit in series with; this is native-`.stb`-only. See
`docs/examples/d_device_mode_cs_cgd.cir`.

## Interpreting margins, including multiple crossings

Both tools compute margins the same way (`tools/stbcore.margins()` in
Python, mirrored exactly in `stban.c`'s end-of-sweep post-processing):

- **Phase margin**: at every frequency where `|T|` crosses 0 dB,
  `PM = 180 + phase(T)` (phase **unwrapped** continuously from the first
  swept point, not wrapped to any particular range). If several crossings
  exist, **all** are found; the worst (lowest) margin is what gets
  reported and printed, and native prints a `(N crossings found; the
  lowest margin is reported)` note when `N > 1` (same in `tools/stb.py`'s
  `describe()`).
- **Gain margin**: at every frequency where `phase(T)` crosses `-180 + 360k`
  for any integer `k`, `GM = -20*log10|T|`. Same worst-of-several handling.
- **No crossing is reported as `None` (Python) / `NaN` (native), never as
  a number** -- "phase margin: no 0 dB crossing of T" / "gain margin: no
  -180 deg crossing of T" are legitimate, common outputs, not errors.
  Every example in this doc that doesn't cross -180 deg says so plainly
  rather than printing a fake gain margin.
- **Falling vs rising 0 dB crossings** (card #464): every crossing is
  labelled falling (`|T|` decreasing through it -- the crossing classic
  Bode-plot stability intuition is about) or rising (`|T|` increasing
  through it). The printed phase margin names which one it came from
  (`"(falling 0 dB crossing)"` / `"(rising 0 dB crossing)"`); a rising
  crossing is not the textbook stability margin even when its `PM` happens
  to read under 180 deg, and a rising crossing whose `PM` is over 180 deg
  is called out explicitly in the printed text. See example (d) below for
  a genuine rising-crossing case.
- **`low_f_positive_feedback` warning** (card #464, **narrowed by card
  #476** -- see below for why): let `arg0 = arg(T(fmin))` in degrees and
  `mag0 = |T(fmin)|`, both at the **lowest swept frequency**. If
  `|arg0| > 150 deg` **and** `mag0 > 1` (positive feedback WITH loop
  gain), both tools print a `WARNING: positive feedback at low frequency:
  loop is unstable if |T|>1 there; phase/gain margins are not meaningful`
  line (and `tools/stbcore.margins()` sets `low_f_positive_feedback: True`
  / `note` in its returned dict; native sets a matching
  `lowFPositiveFeedback` vector on the `"STB Margins"` plot) in addition
  to the usual PM/GM output. **This is not a hypothetical** -- it is
  exactly what this doc's own example (c) circuit did before its card
  #464 fix (see below): a genuine positive-feedback bug in the circuit
  produced a plausible-looking `PM = 274 deg`, which this section
  previously (and wrongly) explained away as "far from the Nyquist
  point." A loop with `arg0` past +/-150 deg has NOT settled to the sign
  a negative-feedback loop should have, and still has gain there; whether
  it is actually unstable depends on the full Nyquist plot, not just this
  one point (Nyquist, not Bode, is the right tool once this fires) -- the
  PM/GM numbers printed alongside the warning are still computed (never
  suppressed) but must not be read as real margins.

  **Card #476: why this narrowed, and what replaced the gap.** The
  original card #464 trigger was `Re(T(fmin)) < 0` alone (equivalent to
  `|arg0| > 90 deg`), which false-triggered on healthy integrator-type
  loops: a type-III compensated buck's near-integrator compensator branch
  (`R1in`/`C1f` in `tests/fixtures/compopt_r3_averaged_detuned.cir`) has
  `T` genuinely near -90 deg at low frequency, and a few degrees of extra
  lag from the real pole pushes `Re(T)` slightly negative with no real
  instability -- found sweeping that fixture from `fstart=10 Hz`
  (`docs/compopt.md`'s old "Sweep note", now corrected). DECISION (option A,
  card #476): only warn when there is BOTH positive feedback
  (`|arg0| > 150 deg`, not just `> 90 deg`) AND loop gain there
  (`mag0 > 1`); replace the gap with two informational, non-warning
  labels so nothing that used to be visible is silently dropped:
  - **`low_f_integrator`** (`low_f_integrator_note`, a `note:` line, never
    a warning): `arg0` within -90 +/- 30 deg -- ordinary integrator-like
    low-frequency behaviour, ubiquitous on a type-III/PI compensator's
    near-integrator branch.
  - **`low_f_positive_feedback_low_gain`**
    (`low_f_positive_feedback_low_gain_note`, a `note:` line, not the
    warning): `|arg0| > 150 deg` but `mag0 <= 1` -- genuine positive
    feedback at low frequency, but not enough loop gain there to call it
    unstable by this test alone.

  Real captured output, this card's build (`aaba50ba3f49-stb6`,
  `tools/stb.py`, `"dec 50 10 1e10"`, matching `tests/test_stb.py`):

  `tests/fixtures/compopt_r3_averaged_detuned.cir` (the near-integrator
  fixture that motivated this card, its own DETUNED starting point,
  `fstart=10 Hz` -- the exact point that used to false-trigger):
  ```
  phase margin: -5.391 deg at 55838.8 Hz (falling 0 dB crossing)
  gain margin: -3.314 dB at 46778.3 Hz
  note: 2 -180 deg crossings; the lowest margin is reported
  note: integrator-like low-frequency behaviour
  ```
  No `WARNING:` line -- `arg(T(10 Hz))` is close to -90 deg (a near-
  integrator, not positive feedback with gain), so `low_f_integrator`
  fires as a neutral note instead. (This fixture's starting point is a
  genuinely bad, low-margin compensator -- `PM = -5.4 deg` -- but that is
  a real finding from the ordinary falling-crossing margin, not from the
  low-frequency flag; see `docs/compopt.md`.)

  `tests/fixtures/stb6_r3_tuned.cir` (the same near-integrator family,
  TUNED to a healthy, well-margined point -- `docs/compopt.md`'s worked
  example, `r2=27950.5, c2=2.4425e-9`):
  ```
  phase margin: 37.998 deg at 34000.3 Hz (falling 0 dB crossing)
  gain margin: no -180 deg crossing
  note: integrator-like low-frequency behaviour
  ```
  Same neutral `low_f_integrator` note on a genuinely good loop -- proof
  this is a label about low-frequency *shape*, not a judgment on the
  loop's health.

  `tests/fixtures/stb6_posfb_low_gain.cir` (new card #476 fixture: a
  unilateral amplifier, `refs/r1a_unilateral.cir`'s topology, with its
  VCCS control-node order swapped and `gm` reduced so `T(fmin) = -0.333`
  -- real positive feedback, `arg0 = 180 deg`, but `mag0 = 0.333 < 1`):
  ```
  phase margin: no 0 dB crossing
  gain margin: no -180 deg crossing
  note: positive feedback at low frequency, loop gain < 1
  ```
  The note, not the warning -- `mag0 <= 1` there.

  Native `.stb` (same three fixtures, `probe=Vprobe`, same build,
  `dec 50 10 1e10`) prints the same three labels at matching values (see
  `tests/test_stb_native2.py`'s `check_low_f_classification_native`):
  ```
  stb STB Loop Stability Analysis: phase margin -5.391 deg at 55838.8 Hz (falling 0 dB crossing)
  stb STB Loop Stability Analysis: gain margin -3.314 dB at 46778.3 Hz (2 crossings found; the lowest margin is reported)
  stb STB Loop Stability Analysis: note: integrator-like low-frequency behaviour
  ```
  ```
  stb STB Loop Stability Analysis: phase margin 38 deg at 34000.3 Hz (falling 0 dB crossing)
  stb STB Loop Stability Analysis: gain margin: no -180 deg crossing of T
  stb STB Loop Stability Analysis: note: integrator-like low-frequency behaviour
  ```
  ```
  stb STB Loop Stability Analysis: phase margin: no 0 dB crossing of T
  stb STB Loop Stability Analysis: gain margin: no -180 deg crossing of T
  stb STB Loop Stability Analysis: note: positive feedback at low frequency, loop gain < 1
  ```
  And the archived, still-buggy R1c (`tests/fixtures/stb4_r1c_posfb.cir`,
  DM mode) still warns, unchanged, re-checked at `fstart=10 Hz` too (not
  just this doc's usual 1 Hz sweeps) as this card's own red-proof that the
  narrowing did not also swallow the real bug:
  ```
  stb STB Loop Stability Analysis: phase margin 274.4 deg at 6.34704e+06 Hz (falling 0 dB crossing)
  stb STB Loop Stability Analysis: gain margin: no -180 deg crossing of T
  stb STB Loop Stability Analysis: WARNING: positive feedback at low frequency: loop is unstable if |T|>1 there; phase/gain margins are not meaningful
  ```

**A phase margin above 180 deg is possible and is not automatically a
bug**, but it is also not automatically fine -- check which of the two
cases below it is before trusting it:

- **Genuinely far from the Nyquist point, on a real negative-feedback
  loop, at a RISING crossing.** `docs/examples/d_device_mode_cs_cgd.cir`
  (example (d) below): `T` starts near 0 at DC with about +90 deg of lead
  (not a positive-feedback loop -- `Re(T) > 0` throughout, no
  `low_f_positive_feedback` warning), so its 0 dB crossing is *rising*
  (`|T|` increasing through it) and the resulting `PM = 258 deg` reflects
  genuine phase headroom on an upward crossing, not a falling, textbook
  margin. The tools now label this explicitly: `"(rising 0 dB crossing)"`.
- **A positive-feedback bug being misread as a large margin.** This is
  what example (c) (the DM/CM differential amplifier, below) used to be,
  in this exact doc, before card #464: `T(1 Hz) = -12.9`/`-3333`
  (`Re(T) < 0` at DC -- both VCCS control-node orders were swapped, a
  circuit bug, not a math bug in this analysis), printing `PM = 274`/`186`
  deg from ordinary **falling** crossings. A falling crossing with `PM`
  over 180 deg on a loop that never had `Re(T) > 0` at low frequency is
  exactly what `low_f_positive_feedback` now catches -- example (c) is
  fixed now (see below) and no longer triggers it, but the archived,
  still-broken circuit (`tests/fixtures/stb4_r1c_posfb.cir`) is kept
  specifically so this doc's own red-proof tests can confirm the warning
  fires on it, in both tools.

So: **always look at the phase trace (or at least the crossing count and
direction, and the `low_f_positive_feedback` flag) before trusting a
single scalar margin**, especially an unusually large or small one.

## Known limits and caveats

- **DM/CM is exact only for an exactly symmetric network** (see above);
  not validated, and not expected to be exact, on an asymmetric one.
- **`probe=<Vp>,<Vn>` (comma form) is unreliable on the interactive `stb`
  command** -- use `probe=<Vp> probe2=<Vn> mode=dm|cm`.
- **A probe terminal on literal ground is degenerate** (loop mode): `T`
  reads exactly 0, with a warning, not a refusal -- move the break one
  component further from true ground.
- **Device mode nulls gm only**, never gm+gmbs, for both MOS1 and BSIM4;
  body effect can matter (see "Device mode" above).
- **BSIM4 device mode is quasi-static only** (NQS refused) and validated
  only to 100 kHz.
- **`localgnd=` is implemented and validated** (`tests/test_stb_native2.py`'s
  `check_localgnd`: `localgnd=0` matches the default to 1e-12, a nonzero
  local ground gives a measurably different, independently-referenced `T`)
  -- an earlier patch's commit message called it unvalidated; that was
  fixed by a later card (see "How this was validated").
- **`probe2=` on a netlist card was never tested** here (see "Netlist-card
  quoting" above) -- use the interactive form for DM/CM if in doubt.
- No Spectre-style `values`/`valuesfile`/`log` frequency sweeps, and no
  `dev`/`mod`/`param`/`temp` sweep at fixed frequency.
- **No `pstb`** -- no periodic/PSS-linked stability analysis in this
  project as of this card; native `.stb` is small-signal about one
  ordinary DC operating point, same restriction plain Spectre `stb` has.
- `probe=<inst>:<terminal>` needs the instance to already exist earlier in
  the netlist (forward references are a hard parse error).

## How this was validated

Five independent test scripts, all re-run 2026-09-28 against this card's
(#464) build (`~/Work/apps/ngspice_augmentation/aaba50ba3f49-stb4/bin/ngspice-aug`,
patches 0001-0010 plus `patches-wip/stb4/0011` applied on top of the
`pre-master-48` pin -- `buildrecipe/README.md`, `patches-wip/stb4/README.md`):
`tests/test_stb.py` (the pure Python/system-ngspice testbench against
closed-form, independent nodal, and ngspice Y-parameter references, plus
this card's new `Re(T(lowest f))>0`/`low_f_positive_feedback` sanity and
red-proof checks) **34/34 pass** (was 22/22 before this card);
`tests/test_stb_native.py` (native loop mode vs `tools/stb.py`,
DC-identity, probe-reversal, KLU option, quoted netlist-card form)
**28/29 pass** (1 `skip`, by design: R1c's 4-terminal diffprobe has no
single-ended native equivalent -- unchanged by this card);
`tests/test_stb_native2.py` (`probe=<inst>:<terminal>`, `localgnd=`,
native `mode=dm|cm`, plus this card's new native warning/sign checks)
**44/44 pass** (was 30/30); `tests/test_stb_device.py` (MOS1 device mode
vs an independent brute-force small-signal Bode reference) **16/16 pass**
(unchanged); `tests/test_stb_followups.py` (grounded-probe degeneracy fix,
BSIM4 device mode, MOS1 reverse mode, plus the full regression chain
through all four prior scripts) **29/29 pass** (unchanged). Every script's
checks include red-proofs (a deliberately mutated or bogus comparison,
which must then fail) so a vacuously loose bound cannot pass silently --
this card's own red-proof is `tests/fixtures/stb4_r1c_posfb.cir`, the
archived pre-fix R1c, which must (and does) trip `low_f_positive_feedback`
in both tools, in both DM and CM mode. Typical native-vs-script agreement
on loop mode is 1e-9 to 1e-14 relative (one exception found and explained
this pass: the fixed R1c's CM leg falls to `|T|~3e-9` by the top of a
`1e10` Hz sweep, where native/script agreement is ~1e-7 relative but
~3e-16 absolute -- double-precision noise on an already-negligible signal,
not a real mismatch; see `tests/test_stb_native2.py`'s `check_dm_cm`);
device mode's own bound is 1e-6 (met, generally with margin -- see each
script's printed `max rel err`).

**Card #466, 2026-09-29** (this pass): re-run against
`~/Work/apps/ngspice_augmentation/aaba50ba3f49-stb5/bin/ngspice-aug`,
patches 0001-0011 (now all promoted to `patches/` -- see below) plus
`patches-wip/stb5/0012` (this card's fix) applied on top of the pin.
Same five scripts, unchanged counts: `test_stb.py` **34/34**,
`test_stb_native.py` **28/29** (1 skip, unchanged), `test_stb_native2.py`
**44/44**, `test_stb_device.py` **16/16**, `test_stb_followups.py`
**29/29** -- this card's fix (the `dot_stb` parameter-validation bypass
and the netlist-card auto-quoting) touches only `.stb`'s own frequency/
numsteps parsing and a new, narrowly-scoped preprocessing step, and
regresses none of the existing behaviour these scripts check. Plus a new
sixth script, `tests/test_stb_cli.py` (this card): the COULD-NOT-RUN-proof
for the parameter-validation fix (all six bad-value cases, both forms,
asserting no `loopgain` vector is produced) and the unquoted-vs-quoted
numeric-equivalence checks for the auto-quote fix -- **42/42 pass**.

**Card #476, 2026-09-29** (this pass, private fork worktree
`ngspice-src-stb6`, branch `aug-stb6` = pin + `patches/0001`-`0012`):
re-run against `~/Work/apps/ngspice_augmentation/aaba50ba3f49-stb6/bin/ngspice-aug`,
patches 0001-0012 plus `patches-wip/stb6/0013` (this card's narrowed
`low_f_positive_feedback` trigger, `low_f_integrator` label, and
`low_f_positive_feedback_low_gain` note -- see "Interpreting margins"
above) applied on top of the pin. All six scripts green, with new checks
added: `tests/test_stb.py` **60/60** (was 34/34 -- +26: fstart=10 Hz
re-checks of every R1a/R1b/R1c/R1d reference and the archived-posfb
red-proof, the compopt R3 fixture (detuned and its new TUNED fixture,
`tests/fixtures/stb6_r3_tuned.cir`), a new low-gain-positive-feedback
fixture (`tests/fixtures/stb6_posfb_low_gain.cir`), and six synthetic
`stbcore.margins()` classification checks); `tests/test_stb_native.py`
**28/29** (1 skip, unchanged); `tests/test_stb_native2.py` **54/54** (was
44/44 -- +10: `check_low_f_classification_native`, the same
detuned/tuned/low-gain fixtures run natively and compared to
`tools/stb.py`, plus the fstart=10 Hz archived-R1c red-proof);
`tests/test_stb_device.py` **16/16** (unchanged); `tests/test_stb_followups.py`
**29/29** (unchanged); `tests/test_stb_cli.py` **42/42** (unchanged). Plus
`tests/test_compopt.py` (not part of the original five/six -- validates
`tools/compopt.py`, which imports `stbcore.margins()` unmodified and reads
its `low_f_positive_feedback`/`note` fields to refuse degenerate
candidates) **20/20 pass**, confirming the narrower flag does not regress
that consumer (`compopt_r3_averaged_detuned.cir`'s own optimizer run no
longer needs to route around a false positive -- `docs/compopt.md`'s old
"Sweep note" is corrected). Every new low-frequency-classification check
has a red-proof or an independent cross-check (native vs `tools/stb.py`,
or a closed-form `T(0)` for the new fixtures) -- see
`tools/stbcore.margins()`'s docstring and this section above for the
actual thresholds and real captured output.

## Building and running the fork binary

`buildrecipe/README.md` has the full recipe (clone, dependencies, build,
prefix layout, and the rollback line: delete the prefix directory --
nothing under `/usr/local` is ever touched). Short version:

```sh
buildrecipe/build.sh --patches --prefix-name=<name>
NGSPICE=~/Work/apps/ngspice_augmentation/<name>/bin/ngspice-aug python3 tests/test_stb_native.py
```

This doc's own examples were re-run against the `aaba50ba3f49-stb5`
prefix for card #466 (see the top of this doc and "How this was
validated" above for exactly which patches that build has).

## Examples

All four run and their real printed output pasted in below. Decks are
under `docs/examples/`; any output files they produce land only under
`~/Work/data/ngspice_augmentation/docs-examples/` (never the repository).

### (a) LDO-style loop, `tools/stb.py` + `lgprobe`

`docs/examples/a_ldo_loop_lgprobe.cir` -- a two-stage LDO-shaped loop
(single-pole error amp, transconductance pass device into an output
filter with ESR, resistive feedback divider), probed on the feedback wire.

```sh
python3 tools/stb.py docs/examples/a_ldo_loop_lgprobe.cir --probe XP --sweep "dec 30 10 10meg"
```

```
phase margin: 35.819 deg at 43488.8 Hz
gain margin: no -180 deg crossing
```

(DC loop gain 59.9 dB; the ESR zero at ~67.7 kHz sits close enough to the
43.5 kHz crossover to erode phase margin -- a realistic, moderately-margined
loop, not a cherry-picked textbook 60 deg case.)

### (b) Same LDO loop, native `.stb`

`docs/examples/b_ldo_loop_native.cir` -- the **same** circuit and the
**same** break (feedback wire), with the `lgprobe` subcircuit replaced by
a single `Vprobe fb fbtop DC 0` (native `.stb`'s `probe=` resolves a plain
V-source instance, not a subcircuit).

```sh
~/Work/apps/ngspice_augmentation/aaba50ba3f49-stb5/bin/ngspice-aug -b docs/examples/b_ldo_loop_native.cir
```

```
stb STB Loop Stability Analysis: phase margin 35.82 deg at 43488.8 Hz (falling 0 dB crossing)
stb STB Loop Stability Analysis: gain margin: no -180 deg crossing of T
```

**Agreement with (a):** same break, same printed phase margin frequency
(43488.8 Hz) and value to the 4 significant figures native prints
(35.819 vs 35.82 deg), and now both label it the same way
(`falling 0 dB crossing`). Comparing the full complex `T(f)` traces of
both runs point-by-point over the swept range: this doc previously
claimed a maximum relative error of 3.8e-10, sourced (it turns out) from
`tools/stb.py`'s `stb.csv`, which only writes 9 significant digits --
that bound was measuring CSV-text rounding, not the tools' true
agreement. Recomputed for card #466 comparing full double-precision `T`
values on both sides (`tools/stb.py`'s `run()` return value directly, not
the rounded CSV; native's `wrdata`, full `numdgt=15` precision) gives a
**maximum relative error of 2.3e-13** over all 181 swept points --
comfortably inside this project's stated 1e-9-to-1e-14 native-vs-script
range (see "How this was validated" below). The independent reader whose
findings prompted this pass (verify-v-g4, 2026-09-29) separately measured
3.3e-9 for this same comparison; that number could not be reproduced with
the full-precision methodology above (frequency grids match exactly
point-for-point between the two runs; see the worked example in this
card's report), and its exact derivation is not known -- flagging this
rather than silently dropping it, since this doc must not report a
number it cannot independently stand behind.

### (c) Differential amplifier, DM/CM

A fully differential transconductor with an odd-symmetric DM path and a
separate CMFB loop (so DM and CM are genuinely different loops with
different compensation, not the same number twice) -- same topology
family as `refs/r1c_fully_diff.cir` (independently re-verified here, not
copied: rebuilding this exact circuit and running it gave the identical
printed numbers below).

**Card #464, 2026-09-28: this example was wrong until this pass.** The
original `G1p`/`G1n`/`Gcmfb` control-node order gave POSITIVE feedback in
BOTH the DM and CMFB loops (`T(1 Hz) = -12.9` DM, `-3333` CM), which this
very doc previously (and wrongly) explained as a legitimate "PM > 180 deg"
case -- see "Interpreting margins" below for the corrected explanation and
how the new `low_f_positive_feedback` warning catches exactly this. Fixed
by swapping each VCCS's control-node order (flips the sign of the current
it injects, touches no pole) and lowering `gmcorr` (1m -> 8u) so the CM
loop's margin lands in a genuine, non-degenerate range instead of near
zero at the loop's now much higher effective crossover. Stability of the
fixed circuit is also confirmed in the time domain: a 1 uA transient kick
into `outp`/`outn` (DM: opposite polarity; CM: same polarity) decays
smoothly back to numerical noise (`<1e-20`) within a few microseconds in
both modes, never diverging (the pre-fix circuit's equivalent kick
diverged to `~1e29 V` in 2 us) -- see the card #464 report for the traces.

`docs/examples/c_diffamp_dmcm_diffprobe.cir` (`tools/stb.py`,
`stb/diffprobe.lib`):

```sh
python3 tools/stb.py docs/examples/c_diffamp_dmcm_diffprobe.cir --probe XP --mode dm --sweep "dec 30 1 100meg"
python3 tools/stb.py docs/examples/c_diffamp_dmcm_diffprobe.cir --probe XP --mode cm --sweep "dec 30 1 100meg"
```

```
mode=dm: phase margin: 94.446 deg at 6.34704e+06 Hz (falling 0 dB crossing)
         gain margin: no -180 deg crossing
mode=cm: phase margin: 56.880 deg at 343956 Hz (falling 0 dB crossing)
         gain margin: no -180 deg crossing
```

`docs/examples/c_diffamp_dmcm_native.cir` (native, `probe=Vp probe2=Vn
mode=dm|cm`, same circuit, two plain V-sources instead of `diffprobe`):

```
mode=dm: phase margin 94.45 deg at 6.34704e+06 Hz (falling 0 dB crossing)
         gain margin: no -180 deg crossing of T
mode=cm: phase margin 56.88 deg at 343956 Hz (falling 0 dB crossing)
         gain margin: no -180 deg crossing of T
```

Native and script agree to the precision printed (already validated
project-wide to ~1e-9 for this exact comparison outside CM's last decade,
where `|T_cm|` falls to `~3e-9` and the comparison is dominated by
double-precision noise -- see `tests/test_stb_native2.py`'s `check_dm_cm`
for the ~1e-7-at-~3e-9-magnitude detail). Both margins are now ordinary,
falling-crossing phase margins under 180 deg (DM 94.4 deg -- a single real
pole's crossing is always just above 90 deg, see below; CM 56.9 deg, a
clearly different, lower crossover from the added CMFB pole) and **neither
run raises the `low_f_positive_feedback` warning** -- `Re(T(1 Hz))` is
`+12.9` (DM) and `+26.7` (CM), a genuine negative-feedback loop in both
modes. Compare `tests/fixtures/stb4_r1c_posfb.cir` (the archived, still-
buggy circuit): running the same two commands against it reproduces the
old 274.4/185.6 deg numbers exactly, now WITH the warning printed --
this is this project's own red-proof for the fix, not just a narrative.

### (d) Device mode, common-source with Cgd

`docs/examples/d_device_mode_cs_cgd.cir` -- a common-source stage with a
large discrete gate-drain (Miller) capacitor: the only "loop" is inside
the transistor, through Cgd, invisible to any wire-break probe (Tian
applicability case 2). `tools/stb.py` has no equivalent for this.

```sh
~/Work/apps/ngspice_augmentation/aaba50ba3f49-stb5/bin/ngspice-aug -b docs/examples/d_device_mode_cs_cgd.cir
```

```
stb STB Loop Stability Analysis: phase margin 258.3 deg at 3288.26 Hz (rising 0 dB crossing)
     -- RISING |T| crossing, not the usual falling-crossing stability margin
        (a value over 180 deg here is not a real margin)
stb STB Loop Stability Analysis: gain margin: no -180 deg crossing of T
```

(re-run against card #466's build, `aaba50ba3f49-stb5` -- same
258.3 deg/3288.26 Hz as every earlier run of this example.)

The deck sweeps `dec 20 1 1e8` -- **1 Hz to 100 MHz, not 1 GHz** (this doc
previously said 1 GHz here; corrected). `|T|` (read from the `stb1.loopgain`
vector, `wrdata`) rises from **-70.2 dB at 1 Hz to +13.9 dB, asymptotic
from about 10 MHz to the sweep's actual end at 100 MHz** (this doc
previously said -86.9 dB to -0.4 dB at 1 GHz -- both the endpoint values
and the sweep's actual top frequency were wrong; recomputed directly from
this pass's `wrdata stb1.loopgain` trace: |loopgain(1 Hz)| = 3.106e-4 =
-70.156 dB, |loopgain(1e7 Hz)| = 4.9295 = +13.856 dB, essentially flat
through |loopgain(1e8 Hz)| = 4.9295 = +13.856 dB) -- the local Cgd loop
clearly showing up, exactly the behavior a loop-mode probe could never
see. Cross-checked against an independent brute-force
small-signal reference (same technique as `tests/test_stb_device.py`:
replace the transistor with its own op-point-derived Rds/Cgd/gm
equivalent, disconnect gm's control port, drive it with 1 V AC, read back
`T = -v_returned`): **maximum relative error 6.4e-7** over the same sweep,
comfortably inside this project's 1e-6 device-mode bound.

**This is the legitimate rising-crossing case, not a positive-feedback
bug** (contrast with example (c)'s history, above): `Re(T) > 0` throughout
this sweep (`T` starts near 0 at DC with about +90 deg of lead, never
negative), so `low_f_positive_feedback` does **not** fire here -- the
circuit itself is a genuine, correctly-behaved negative-feedback-ish local
loop. The `(rising 0 dB crossing)` label and the "not the usual
falling-crossing stability margin" note are about how to *read* a PM
computed at an upward `|T|` crossing (it is not the classic Bode-plot
distance-to-instability number), not a claim that anything is wrong with
this circuit or this `T` -- `T` itself is independently verified to 6.4e-7
above. Keep this example; just don't read its `258.3 deg` the way you
would a falling-crossing op-amp margin.

## What was not done

- No independent re-verification by a second reader -- this doc,
  its examples, and their numbers are all from a single session.
- Not every native `.stb` parameter combination was exercised by this
  doc's own examples (in particular `localgnd=` and `probe=<inst>:<terminal>`
  are documented from reading `tests/test_stb_native2.py` and its results,
  not re-run fresh into a dedicated docs example here -- time budget).
- No `pstb`/PSS-linked stability work exists in this project yet, so
  nothing about it could be documented beyond "not implemented."

**Card #464 additions (2026-09-28), what was not done there:**

- The `low_f_positive_feedback`/falling-vs-rising warning was added and
  red-proofed only for **loop mode** (`tools/stb.py`'s `run`/`run_diff`
  and native `.stb` loop mode, including `mode=dm|cm`) -- **not** wired
  into native **device mode** (`probe=<MOS1/BSIM4 instance>`, `stban.c`'s
  separate device-mode branch); device-mode's own margin post-processing
  shares the same code path this pass edited, so the warning likely fires
  there too as a side effect, but this was not specifically tested against
  a positive-feedback device-mode case (none of this project's existing
  device-mode reference decks are positive feedback) and is not claimed.
- No search was made for OTHER pre-existing positive-feedback bugs in this
  project's reference decks beyond R1c -- R1a/R1b/R1d were checked (all
  `Re(T)>0` at the lowest swept frequency, see `tests/test_stb.py`'s new
  checks) because this card's brief named them, not because every deck in
  `refs/` or `docs/examples/` was audited.
- `patches-wip/stb4/0011` was built and tested against `ngspice-src-stb4`
  only at the time this note was written, and had not yet been promoted
  to `patches/`, nor cross-checked by regenerating it from a fresh
  pin+series `git am` -- **both since done** (found stale during card
  #466's pass, 2026-09-29): patch 0011 is now `patches/0011-*` (part of
  the promoted series, `patches/series`), and its reproducibility (pin +
  `patches/0001`-`0011`, `git am`, diff against `aug-stb5` HEAD before
  this card's own commit) was re-checked as part of card #466's own
  reproducibility check (`patches-wip/stb5/README.md`).
- The transient stability kicks (DM and CM, 1 uA, decay confirmed within
  a few us) were run as one-off decks under
  `~/Work/data/ngspice_augmentation/stb4/` (`r1c_kick_dm.cir`,
  `r1c_kick_cm.cir`), not checked into the repo as a test -- the
  frequency-domain half-circuit comparisons in `tests/test_stb.py` and
  `tests/test_stb_native2.py` are the checked-in proof of stability;
  the transient kicks are corroborating, not gating, evidence.

**Card #466 additions (2026-09-29), what was not done there:**

- `probe2=`/`mode=` on a **netlist-card** `.stb` line used to be flagged
  above as "inferred, not verified, likely needs quoting" -- this is now
  resolved (auto-quoting; `tests/test_stb_cli.py`'s
  `check_autoquote_probe2_mode`), so that caveat is removed rather than
  carried forward.
- Did not trace *why* `.options` (and presumably other stock dot-cards
  with unquoted string keywords) escapes the same numparam evaluation
  path that used to trip up an unquoted `.stb probe=` -- confirmed this is
  a general numparam behaviour, not `.stb`-specific (a bogus, unrecognized
  `.foocard bar=<bareword>` line fails identically), but did not trace the
  exact mechanism that exempts `.options`; see `patches-wip/stb5/README.md`
  for the full note. The fix (auto-quoting `.stb`'s own options before
  numparam runs) does not depend on understanding this.
- Did not audit every other dot-card in this fork for the same unquoted-
  string-keyword problem -- only `.stb` was in scope for this card.
- The device-mode magnitude/crossing numbers in example (d) above were
  re-verified (real, re-run output), but the **6.4e-7** device-mode
  reference-agreement bound quoted there was not independently
  recomputed this pass (not one of the reader's flagged numbers) --
  carried forward from the prior pass's measurement, unchanged since
  device mode's own code was not touched by this card's fix.
- This doc's own claim of a 3.3e-9 native-vs-script agreement figure
  (attributed to the independent reader whose findings prompted this
  pass) could not be reproduced with this pass's full-precision
  methodology (2.3e-13 measured instead, see example (b) above); the
  discrepancy's cause was not chased down further -- both figures are
  reported, honestly, rather than one being silently preferred.

**Card #476 additions (2026-09-29), what was not done there:**

- Same device-mode gap as card #464 left open, still open: the narrowed
  `low_f_positive_feedback`/`low_f_integrator`/`low_f_positive_feedback_low_gain`
  classification was added and tested only for **loop mode** (`tools/stb.py`'s
  `run`/`run_diff` and native `.stb` loop mode, including `mode=dm|cm`) --
  not specifically exercised against native **device mode**
  (`probe=<MOS1/BSIM4 instance>`); device mode's margin post-processing
  shares the same `stban.c` code this card edited, so the same
  classification likely applies there too as a side effect, but this was
  not tested (none of this project's device-mode reference decks are
  near an integrator or positive feedback at low frequency).
- The two new thresholds (`150 deg` for positive feedback, `-90 +/- 30 deg`
  for the integrator label) are the DECISION's own numbers (option
  A, card #476), not independently re-derived from a wider survey of
  circuit types -- they were checked to correctly classify every fixture
  this card's brief named (the archived buggy R1c, the healthy R1a/R1b/R1c/R1d
  references, the near-integrator R3 compensator detuned/tuned, and a new
  low-gain-positive-feedback fixture), not swept against an exhaustive
  range of `arg(T(fmin))`/`|T(fmin)|` combinations beyond the synthetic
  boundary checks in `tests/test_stb.py` (arg = -119/-121/-90/0/170 deg).
- `docs/compopt.md` was updated only where this card's fix directly
  changed a claim already in that file (the `fstart=10 Hz` "Sweep note"
  and the FRA-verification-window note, both about the old, now-replaced
  `low_f_positive_feedback` trigger) -- the rest of that doc (the general
  `low_f_positive_feedback` description in "How an evaluation works",
  scoring/cost-function text, etc.) was left as-is, out of this card's
  scope.
- No independent re-verification by a second reader -- this pass's
  additions, like the rest of this doc, are from a single session.
