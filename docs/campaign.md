# Campaign-level scheduling (`tools/campaign.py`)

Card #491, plan T9.6, transient-ACCELERATION round 6. Five prior rounds
(`research/transient-acceleration.md` sections 1-10, `patches-wip/accel/`)
squeezed single-run speedups out of the fork itself. This round's lever is
CAMPAIGN level: a PVT sweep, a Monte Carlo run set, or any other
many-point sweep is many independent ngspice invocations, and the naive
"run them one after another" approach leaves two things on the table: idle
cores (parallel dispatch) and repeated parse/setup cost (a single-process
`alterparam`/`option temp`/`reset` loop, already validated by this repo's
own `bench/gen_pvt.py`). Memory is the ceiling on how far parallel dispatch
can go, not core count alone.

`tools/campaign.py` is a **generic scheduler/runner engine** (see its own
module docstring for the full design). A "campaign" is a small Python
module supplying `POINTS` + `build_point_deck()` (mode a) and/or
`build_loop_prelude()`/`build_loop_entry()`/`build_loop_tail()` (mode b).
Two worked campaigns ship with this round:

- `bench/campaign_pvt_r3d.py` -- a new 36-point PVT campaign (temp x Vin x
  Rload) on `refs/r3d_vm_buck_deadtime.cir` (card #479's dead-time buck).
- `bench/campaign_mc_diffpair.py` -- a 200-run Monte Carlo campaign on the
  differential-pair offset example (`docs/examples/mc_diffpair_offset.cir`,
  card #473), reusing `tools/mc.py`'s own seeded sampling scheme.
- `bench/campaign_pvt_bench_adapter.py` -- wraps the EXISTING `c_pvt` bench
  artifact (`bench/gen_pvt.py`, card #467 T9.1) so it can be run through
  the new scheduler too, unmodified.

## Modes

- **(a) one process per point.** `run_campaign_a()`: probes the first
  `--probe-n` points serially under `/usr/bin/time -v` for peak RSS, then
  runs the rest at `min(--max-procs, cores, floor((mem_budget -
  headroom) / per_point_rss))` concurrency (`pick_concurrency()`).
  `--mem-budget-mb` defaults to the box's *currently free* memory
  (`/proc/meminfo MemAvailable`), not total RAM -- other users' resident
  work is not this campaign's to evict. `--per-point-mem-cap-mb` further
  wraps each point in its own `systemd-run --user --scope -p
  MemoryMax=...` scope (the same idiom `bench/runner.py` already uses) as
  a safety net against one point ballooning past the probe's estimate.
- **(b) single-process loop.** `run_campaign_b()`: one ngspice session,
  `alterparam`/`option temp`/`reset`/`run`/`meas` per point -- the exact
  idiom `bench/gen_pvt.py` already validated, generalised so any campaign
  module can use it. Card #473's own finding (an in-deck `.param =
  agauss(...)` is re-evaluated per *use*, not bound once) does not apply
  here: every value a campaign loop pushes in via `alterparam` is a plain
  Python-sampled literal, never an in-deck random function (`once()`,
  patch 0016, is the fix for the different case of wanting an in-ngspice
  random draw frozen -- not needed by anything in this module).
- **(c) checkpoint-shared corner settling.** `run_campaign_c()`: for the
  narrower case where several DIFFERENT analyses share the same
  (T, Vin, Rload, ...) corner, `savestate` once and `loadstate` per
  downstream view instead of re-running the shared settle prefix each
  time -- the same mechanism `tools/fra.py --checkpoint` already uses in
  production (docs/checkpoint.md: 1.4-3.1x measured there). Not
  separately re-measured this round (the two worked campaigns above are
  PVT/MC point sweeps across DIFFERENT corners, not several views of the
  SAME corner -- `run_campaign_c`'s plumbing is built and documented, its
  own correctness rests entirely on docs/checkpoint.md's already-measured
  checkpoint/restart behaviour, not re-derived here). **What was not
  done**, stated plainly.

## Measured: naive serial vs campaign.py, each mode

Every number below: same binary
(`~/Work/apps/ngspice_augmentation/aaba50ba3f49-s39/bin/ngspice-aug`, the
current series build, pin + patches/0001-0039), development workstation (AMD Ryzen AI Max+ 395, Linux), `--max-procs 6` for
every scheduled run (this project's own concurrency cap is <=12 concurrent
ngspice total, shared with other users' free-running PSS/R3d/R4d FRA work
running concurrently on this box -- 6 leaves headroom). Single runs, not
n=5 medians (this round's time budget; each number is one measured
wall-clock, not extrapolated -- see "what was not done" below for the
honest limitation this implies).

| campaign | naive serial | mode (a) scheduled | mode (b) loop | best speedup |
|---|---|---|---|---|
| `c_pvt` (existing bench, 50 pts) | 3.70s | 0.81s (concurrency=6) | 3.09s (existing `c_pvt_loop.cir`) | **4.6x** (mode a) |
| `campaign_pvt_r3d` (new, 36 pts) | 107.75s | 28.77s (concurrency=6) | 116.83s (**slower** than naive) | **3.75x** (mode a) |
| `campaign_mc_diffpair` (200 MC runs, `.op` only) | 0.852s | 0.180s (concurrency=6) | 0.095s | **8.97x** (mode b) |

**Read the r3d row carefully -- mode (b) is a real regression there, not a
rounding error.** r3d's dead-time network makes per-point transient cost
large relative to parse/setup cost, so amortising setup (mode b's whole
mechanism) barely matters, while `alterparam`+`reset`'s own per-point
overhead (re-elaborating the circuit in place) turns out to cost slightly
MORE than a fresh process's parse for this deck -- 116.83s vs 107.75s
naive. Mode (a)'s win comes entirely from exploiting idle cores, which
works regardless of the parse/compute ratio; mode (b)'s win comes from
amortising a cost that has to actually be large relative to the point's
own compute to matter. `campaign_mc_diffpair`'s `.op`-only points are the
opposite extreme (compute is trivial, process-launch overhead is
everything), which is exactly why mode (b) wins biggest there (8.97x,
beating even 6-way parallel dispatch of the same 200 launches).
**Takeaway, stated as a rule of thumb because it is is exactly what these
three campaigns demonstrate, not asserted in the abstract: mode (a) is the
safe default (idle-core exploitation is close to workload-shape-agnostic,
per the table above); mode (b) is worth trying only when parse/setup is a
real fraction of a point's own cost, and must be validated per-deck (see
below) before being trusted.**

## Identity checks

`accel_same_results`-style discipline: per-point results compared against
the naive run, tolerance stated BEFORE the run (see `tests/test_campaign.py`
module docstring for the full reasoning), not fitted after.

- **`c_pvt`**: mode (a) vs naive -- **bit-identical**, all 50 points x 4
  measures (200/200 exact). mode (b) (the existing `c_pvt_loop.cir`) vs
  naive -- the two forms of this deck use DIFFERENT `.meas` names for the
  same quantity (`vout_avg_final` in the separate-process form,
  `p{NNN}_vout_avg` in the loop form -- a pre-existing `bench/gen_pvt.py`
  artifact quirk, not touched this round); the two measures with matching
  names (`wave_mid`, `wave_end`) agree to **1e-5 relative, all 50 points**
  (0 exact matches -- expected, see below).
- **`campaign_pvt_r3d`**: mode (a) vs naive -- **bit-identical**, all 36
  points x 4 measures (144/144 exact; `pick_concurrency`'s scheduling
  choice cannot change a single-process, single-point result, and doesn't).
  mode (b) vs naive -- **NOT within the originally-stated tolerance at
  full scale.** A first-pass 5-point SAMPLE (what `tests/test_campaign.py`
  actually gates on, by design -- see that file's own docstring for why
  it checks a sample, not the full campaign) passed comfortably within
  1e-3 (smooth measures) / 1e-2 (edge-crossing measures: `t_period`,
  `duty_final`, and `vout_ripple` -- a peak-to-peak measure, grouped with
  the edge-crossing tier because it depends on the stepper landing a
  sample near the true extremum, the same class of timing-sensitivity as
  an edge crossing). The FULL 36-point run does not: `duty_final` and
  `vout_ripple` exceed 1e-2 relative at several points, concentrated at
  the `temp=85` corners (up to ~2% and ~4% relative respectively) --
  plausibly because the dead-time steering diodes (`ddt`, a real,
  temperature-dependent exponential `D` device) are more sensitive to the
  same floating-point-non-associativity effect docs/checkpoint.md already
  characterises for checkpoint/restart on this exact family of decks, at
  higher temperature. **`vout_avg_final` and `t_period` stay within
  tolerance across all 36 points; `duty_final`/`vout_ripple` do not.**
  Stated as a finding, not smoothed over: **a small-sample identity check
  is not sufficient to certify mode (b) for a new deck at production
  scale -- the divergence here is corner-dependent, and the sample that
  passed did not happen to include the worst corners.** Combined with
  mode (b) also being *slower* than naive on this deck (table above), the
  recommendation for `campaign_pvt_r3d`-shaped campaigns is unambiguous:
  **use mode (a).** It is bit-identical AND faster; mode (b) offers
  neither advantage here.

  **Update, card #492: root cause found, and it is not floating-point
  non-associativity.** The paragraph above's "plausibly" was wrong to
  reach for docs/checkpoint.md's non-associativity explanation --
  research/reset-leak.md has the real cause, a state-leak bug in
  `reset` itself (upstream ngspice, reproduced identically on stock
  `/usr/local/bin/ngspice` 45.2 and the unpatched pin, **not** introduced
  by any patch in this fork's series): `option temp = {temp}` (the ONLY
  way this loop varies temperature -- `alterparam` only touches
  `vin_v`/`rload`) is silently discarded by the very next `reset`, so
  EVERY point in the loop actually ran at ngspice's hard-coded default
  (27C), not its intended temperature. `t_period` (set purely by the
  fixed-frequency `Vramp` source) and `vout_avg_final` (inside the
  closed voltage-mode feedback loop, which servos out a wrong duty cycle
  almost entirely) could not show this regardless; `duty_final`/
  `vout_ripple` (downstream of the one genuinely temperature-sensitive
  element in the timing path, the dead-time steering diodes' exponential
  `IS(T)`) show it in exact proportion to how far each point's *actual*
  (bugged) temperature was from its *intended* one -- which is why the
  divergence was "concentrated at temp=85" (58C from the 27C default,
  vs. 47C for temp=-20, an asymmetry the diode's own exponential
  amplifies) rather than uniform.

  Fixed in `patches-wip/reset/0800-...patch` (`com_rset()` now replays
  `ft_curckt->ci_vars` -- the interactive-`option`-command cache --
  across `reset` instead of silently dropping it; see
  research/reset-leak.md for the full trace). Re-measured on the FIXED
  build (`aaba50ba3f49-reset1` = `aaba50ba3f49-s39` + this one commit):
  mode (b) vs mode (a), **bit-identical, the FULL 36 points x 4 measures
  (144/144 exact)** -- not merely "within tolerance": the fix makes this
  an exact match, the same as mode (a) vs naive already was. The
  slowness finding (table above, mode (b) 116.83s vs naive 107.75s) is
  UNCHANGED by this fix (it was never about correctness) -- **use mode
  (a) for this deck** remains the recommendation, now for the
  speed reason alone rather than speed-and-correctness.

  **Re-confirmed, card #496 (transient-ACCELERATION round 7)**, on this
  round's own mandated baseline `aaba50ba3f49-s46` (pin + 46 patches,
  patch 0041 = this section's fix, a later build than `reset1` above):
  mode (a) vs mode (b), full 36-point `campaign_pvt_r3d`, **still
  bit-identical -- 0/144 mismatches at rtol=0/atol=0** (not just within
  round 6's own 1e-3/1e-2 tolerance tiers, checked separately, also 0
  mismatches). This also fills in the number this section's own table
  never stated directly: **mode (b) vs mode (a)** (not vs naive). Under
  this round's own "at most 2 concurrent ngspice processes" rule (a
  tighter cap than the `--max-procs 6` used for the table above, so not
  directly comparable to the 28.77s mode-(a) figure), a fresh measurement
  gave naive 123.67s, mode (a) (`--max-procs 2`) 62.66s, mode (b) 129.32s
  -- **mode (b) is 2.06x SLOWER than mode (a)** on this deck (and would be
  further behind at higher mode-(a) concurrency, since mode (b) is
  single-process by construction and cannot benefit from more cores at
  all, while mode (a)'s advantage only grows with the concurrency cap).
  Research doc: `research/transient-acceleration.md` section 12.4.
- **`campaign_mc_diffpair`**: mode (a) AND mode (b), both vs naive --
  **bit-identical, all 200 points.** A bare `.op` solve has no adaptive-
  step history to diverge across (unlike a multi-period `.tran`), so
  neither scheduling nor `alterparam`+`reset` has anything to perturb --
  consistent with, not contradicting, the r3d finding above (the
  mechanism that causes r3d's mode-(b) divergence, hundreds of
  independent switching-related step decisions, simply does not exist in
  a single DC operating-point solve).

## Memory-aware scheduling: the red-proof

`tests/test_campaign.py`'s `check_memcap_red_proof` (see that file):
`bench/decks/campaign_hungry.cir` (`bench/gen_pex_buck.py --mesh-n 150
--tstop-us 1`, ~141MB probed RSS/point, ~6.6s/point) run 4-way under a
REAL `systemd-run --user --scope -p MemoryMax=220M` cgroup cap (220MB, ~1.5x
one point's own RSS -- 4-way naive concurrency needs 560MB+ and cannot
fit):

| path | concurrency | result |
|---|---|---|
| memory-aware (`campaign.py`'s own picker, `--mem-budget-mb 220`) | **1** (picked from the probe) | all 4 points complete, 0 OOM |
| naive (`--force-concurrency 4`, bypassing the picker, SAME 220MB cap) | 4 | killed under the cap (systemd terminates the scope; no manifest written) |

The picker is doing real, load-bearing work here, not passing a cap that
was loose enough to not matter either way.

## Solver choice (opt-in): `probe_solver()` / `solver-table`

Runs a deck's own first (short) point twice -- SPARSE forced, KLU forced,
built symmetrically from the SAME base deck text (`_force_sparse()`/
`_inject_option()`) so this is meaningful even when a deck already bakes a
solver choice into its own `.options` line (found while building this:
EVERY deck `bench/gen_pex_buck.py` generates, grid OR tree topology,
unconditionally bakes `klu` into its `.options` line -- an earlier,
naive "as-is vs +klu" version of this prober degenerated to a useless
"klu vs klu" comparison on exactly the two decks, `a_pex_buck`/
`e_pex_tree`, this table most needs to say something about).

| deck | chosen | speedup (klu/sparse) | note |
|---|---|---|---|
| `a_pex_buck` | **klu** | n/a | SPARSE hit a 90s probe timeout without finishing -- consistent with round 5's own ~3.7h extrapolation for this deck; KLU: 50.4s |
| `e_pex_tree` | sparse | 0.94x (klu 6% slower) | matches round 5: same power stage, well-ordered tree mesh, "solver-indifferent" but SPARSE default is right |
| `b_long_switch` | sparse | 0.92x | matches round 5 |
| `c_pvt_loop` | sparse | 0.90x | matches round 5 |
| `d_analog_block` | klu | 1.12x | **disagrees with round 5's own n=5-median finding (SPARSE, klu -10.7%)** -- single-sample only, see caveat below |
| `d_analog_block_long` | klu | 1.09x | consistent direction with `d_analog_block` above, same caveat |
| `f_analog_block_bsim3` | klu | 1.14x | same caveat (sub-second deck) |
| `f_analog_block_bsim3_long` | sparse | 0.99x | essentially a wash, matches round 5's own reasoning (neither solver should matter much here) |

**Every deck with a clear (>5%) margin reproduces round 5's own truth
table exactly** (`a_pex_buck` needs KLU, essentially required; everything
else defaults right with SPARSE). **Caveat, stated plainly, not hidden:
`probe_solver()` runs each side ONCE.** On sub-second decks
(`d_analog_block`, `f_analog_block_bsim3`), a single sample falls inside
this project's own already-documented run-to-run noise band (round 1:
`d_analog_block` measured 0.53-0.66s across repeated runs of the SAME
binary/deck) -- the `klu`-favouring picks on those two rows are NOT
trusted over round 5's own n=5-median SPARSE recommendation; they are
recorded here as an honest limitation of this round's `probe_solver`
(single-sample), not silently smoothed into agreement. **Production use
of `--solver-probe` on a sub-second-per-point campaign should not trust a
single-sample probe; this was not fixed this round (see "what was not
done").**

## Usage

```sh
# mode (a), memory-aware scheduling
python3 tools/campaign.py run-a --module bench/campaign_pvt_r3d.py \
    --ngspice ~/Work/apps/ngspice_augmentation/aaba50ba3f49-s39/bin/ngspice-aug \
    --out ~/Work/data/ngspice_augmentation/campaign/r3d_run1 \
    --max-procs 6 --probe-n 3

# mode (b), single-process loop
python3 tools/campaign.py run-b --module bench/campaign_pvt_r3d.py \
    --ngspice ~/Work/apps/ngspice_augmentation/aaba50ba3f49-s39/bin/ngspice-aug \
    --out ~/Work/data/ngspice_augmentation/campaign/r3d_loop1

# compare two result sets (identity check)
python3 tools/campaign.py compare --a .../manifest_a.json --b .../manifest_b.json \
    --measures vout_avg_final vout_ripple --rtol 1e-3

# solver-choice truth table across a list of decks
python3 tools/campaign.py solver-table --ngspice <bin> \
    --decks bench/decks/a_pex_buck.cir bench/decks/e_pex_tree.cir
```

## What was not done (stated, not hidden)

- **Mode (c) (checkpoint-shared corner settling) was built but not
  separately re-measured this round.** Its correctness rests on
  docs/checkpoint.md's already-measured checkpoint/restart behaviour
  (1.4-3.1x on `tools/fra.py --checkpoint`, the same underlying
  mechanism); neither worked campaign this round has the "several views,
  one corner" shape mode (c) targets. A future round wanting to measure
  it needs a campaign with that shape (e.g. a PVT sweep that ALSO wants a
  load-step response at the same settled corner).
- **`--per-point-mem-cap-mb` (the per-point systemd-run safety net) is
  implemented but not separately red-proofed** -- `check_memcap_red_proof`
  exercises the AGGREGATE `--mem-budget-mb` picker, not the per-point cap.
- **Every wall-clock number in this document is a single run, not an n=5
  median** (this round's time budget, and this project's own concurrency-
  cap constraint on a shared box) -- unlike `bench/README.md`'s own
  n=5/n=7 convention for single-run candidates. Directionally consistent
  with expectations (parallel dispatch beats serial; amortised setup only
  helps when parse cost is a real fraction of point cost) but a single
  sample, stated as such, not as a medianed, noise-characterised number.
  The solver-choice table's sub-second-deck disagreement with round 5
  above is a direct, visible consequence of this.
- **No cross-corner parallelism in mode (c)** -- `run_campaign_c` runs one
  corner at a time (by design, stated in its own docstring); combining it
  with mode (a)'s scheduler to parallelise ACROSS corners while still
  settling once PER corner is a natural next step, not built this round.
- **`--self-cgroup-mb` (re-exec the whole campaign.py process under one
  shared cgroup) was considered and deliberately NOT built** -- the
  red-proof test wraps `campaign.py`'s own CLI invocation with
  `systemd-run` from the caller side instead (simpler, no re-exec
  fragility inside the tool itself); documented here so a future round
  doesn't re-derive the same design choice from scratch.
- **`refs/r4d_pcm_buck_deadtime.cir`** -- not used this round (r3d only,
  per the card's own instruction); a peak-current-mode PVT campaign would
  be a natural second worked example.
