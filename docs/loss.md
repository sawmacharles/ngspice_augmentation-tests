# loss.py -- per-element efficiency / loss breakdown for switching converters

Card #470, plan T10.2, 2026-09-29; extended by card #481 (native `D`/`M`
device coverage, same day) -- see "Card #481" below. Companion tool to
`tools/stb.py` (loop gain) and `tools/fra.py` (frequency response): given a
converged switching period of a converter deck, computes input power, output
power, per-element average dissipated power grouped into categories, the
energy-balance residual, and efficiency.

## Usage

```
tools/loss.py DECK --vin Vin --load Rload --period 2e-6 --periods 5 --tstop 1.5e-3 \
    [--dcr Rdcr] [--esr Resr] [--switch Bhs] [--switch Bls] [--gate-drive NAME] \
    [--diode Dname] [--switch-low 0.10] [--switch-high 0.90] [--out DIR]

tools/loss.py DECK --pss --vin Vin --load Rload --period 2e-6 \
    --fguess 500e3 --stabtime 300e-6 --oscnode vout \
    --ngspice-aug ~/Work/apps/ngspice_augmentation/<fork>/bin/ngspice-aug
```

`--vin` names the input voltage source; `--load` the output load resistor.
`--dcr`/`--esr` name resistors to bucket as inductor-DCR / capacitor-ESR loss
(repeatable). `--switch` names B-source, G, or (card #481) native `M`
(MOSFET) elements that are power switches, split into conduction vs.
transition loss (repeatable). `--diode` names a native `D` element that is a
power diode, bucketed into `diode_conduction` (card #481; includes
reverse-recovery energy -- see "Card #481" below) instead of
`resistive_other` (repeatable). Everything else in the deck (feedback
network, compensator op-amp, probes, comparators, ...) is summed into
`resistive_other`. `--gate-drive` names an explicit gate-driver loss element,
if the deck models one; R3/R4 do not, so that category reports 0 with a "not
modelled" note rather than silently omitting it.

**Exit status (card #481): 0 means a COMPLETE, trustworthy breakdown -- every
device in the deck was accounted for.** 2 means either the run could not be
made at all, OR it ran but at least one device could not be measured (an
unhandled kind, an unresolved subcircuit, a MOSFET inside a `.subckt`, ...);
`report.txt`/`summary.json` still show what WAS computed, and list every such
device under `excluded`, but the run is never reported as a pass in that
case. Before card #481, an unsupported/unresolved device was silently
dropped from the sum and the run still exited 0 -- see "Card #481" below.

Two ways to supply the converged window, matching the spec:
- **`--tstop`/`--periods`** (default path, works on any ngspice): runs a
  settled `.tran` to `--tstop` and analyses the last `--periods` periods
  (need not be an integer -- see Limits).
- **`--pss`** (fork build only, `ngspice-aug`): runs `.pss ... method=shoot`
  and analyses the one converged period it returns.

Outputs under `~/Work/data/ngspice_augmentation/loss/<run>/` (never the
repository): `wrapper.cir` (the deck plus this tool's own control block --
the original deck's own `.control`/analysis lines are stripped and never
edited on disk), `ngspice.log`, `trace.txt` (raw `wrdata`), `summary.json`,
`report.txt`. Exit 0 = computed a breakdown, 2 = could not run/compute.

## Per-element power: how each device type is measured

The single biggest finding of card #470 (2026-09-29) was that **`@device[i]`
and `@device[p]`, read back through `wrdata` inside a `tran ... uic` run on
this ngspice-45.2 build, appeared to not update per accepted timestep** --
frozen at whatever value they had at t=0, for every device type tried (R, L,
C, V, B, G). Card #481 (same day) found the actual cause and the fix, needed
to bring native diodes (`D`) and MOSFETs (`M`) into this tool: **those
vectors were never `.save`d**. See "Card #481" below for the corrected
reproducer and the two traps found getting it working.

Given that, every element's power in this tool is built from two kinds of
primitive: ngspice's own per-timestep-tracked circuit vectors -- node
voltages `v(node)`, native MNA branch currents `i(Vxxx)`/`i(Lxxx)` (including
their subcircuit-qualified form `i(<type>.<inst>.<name>)`, e.g. `i(v.xp.vi)`
for `Vi` inside instance `XP` of `.subckt ... lgprobe`), neither of which
ever needed `save` -- and (card #481) `save`d `@device[param]` vectors, used
only where no native branch current exists and the underlying model equation
is impractical to re-derive in Python (diodes, and MOSFETs' own sense-source
currents). Each device type's own defining relationship then turns those
into P(t):

| device type | how P(t) is obtained |
|---|---|
| V (the input source, and any other independent voltage source, top-level or inside a `.subckt`) | native `i(name)`; `P = (V(n1)-V(n2)) * i(name)` |
| L (inductor) | native `i(name)`; same formula |
| R (resistor) | `@name[resistance]` queried ONCE via a static `print` (never needed `save`, since it is not a per-timestep query) -> `P = (V(n1)-V(n2))^2 / R` from the dynamic voltage trace |
| C (capacitor) | `@name[capacitance]` queried once, statically; charge/energy via `dE = V_avg * C*dV` on each trace interval (equivalent to `I=C dV/dt`, but avoids an explicit division so a repeated-timestamp row in the trace contributes exactly 0 instead of `0/0`) |
| B, `I=expr` (smooth switches Bhs/Bls, comparators, the RS-latch driver) | the controlling expression is parsed from the netlist text, `{param}` references resolved via this tool's own `.param` evaluator, and `V(node)`/`V(n1,n2)` reads rewritten to index a per-timestep voltage dict; the resulting Python expression is evaluated at every sample using the (correctly dynamic) `v(node)` trace -- i.e. exactly "compute v*i from terminal voltages and the controlled current" as the task asked, re-deriving the current from the deck's own formula rather than trusting `@dev[i]` |
| B, `V=expr` (comparators like `Bcmp`, current-sense `Bcs`) | current is undetermined from the expression alone (it defines a voltage, not a current); if the device's own two terminals are the *only* thing electrically attached to those nodes (checked from the full device list, not from what the expression merely *reads* -- reading `V(node)` inside an expression does not draw current), its current must be exactly 0 by KCL, and that is what is reported. If some other device also has a terminal there, the device is excluded from the sum and reported under `excluded` (not silently folded into 0) -- this never happens on R3/R4 (`pwm`, `q`, `vcs`, `reset_sig` are each touched only by their own comparator); it DOES happen on `refs/r3d_vm_buck_deadtime.cir`'s `pwm`/`pwm_inv` (see "Card #481" below) |
| G/VCCS (the compensator op-amp `Goa`) | same idea as `B, I=`: `I = gain*(Vc1-Vc2)`, gain resolved via the `.param` evaluator, evaluated per sample from the dynamic control-node voltages |
| X (subcircuit instances, e.g. `XP ... lgprobe`) | expanded one level (recursively) from the `.subckt` definition text, using ngspice's own `<type>.<inst>.<name>` qualified naming for both node references and `i(...)` queries; supports the 2-terminal R/L/C/V/I/B/G/E/D devices this project's subcircuits (`stb/lgprobe.lib`) actually use |
| I (independent current source) | only a constant `DC 0` is supported (matches `lgprobe.lib`'s `Ii`); anything else is excluded and reported, not silently assumed |
| D (native diode, card #481) | `@name[id]`, `save`d before the run (see "Card #481"); `P = (V(anode)-V(cathode)) * id`. Works at any `.subckt` nesting depth, same as the static `@name[resistance]`/`[capacitance]` queries above |
| M (native MOSFET, any compact-model level, card #481) | NOT `@name[id/is/ig/ib/p]` -- see "Card #481" for why (BSIM4 doesn't expose most of those). Instead: 4 inserted 0V sense sources (one per terminal: drain, gate, source, bulk) in a wrapper-only copy of the deck, `P = sum` over the 4 terminals of `V(original node) * i(sense source there)`. Top-level `M` lines only; one inside a `.subckt` is reported `excluded`, not measured |

Everything is validated against ngspice's own uniform sign convention,
confirmed empirically (isolated single-point and transient reproducers, card
#470): `@device[p]` (where it works) and the formulas above agree, `P>0`
means the device absorbs/dissipates power, `P<0` means it delivers it. Input
power is therefore `P_in = -avg(P_Vin)`.

Switch conduction vs. transition (overlap) split uses the criterion the task
suggested: for each `--switch` device, `rail` = that device's own max
`|V(n1)-V(n2)|` over the window; each trace interval is bucketed by the
*midpoint* fraction `|V(n1)-V(n2)|/rail` -- below `--switch-low` (default
0.10) or above `--switch-high` (0.90) is "conduction" (on-state Ron loss, or
off-state leakage, both lumped there since off-state is negligible either
way), the band in between is "transition". **Honest limitation**: on a
no-dead-time smooth-switch deck like R3/R4 (both switches driven by the same
complementary control signal, so they overlap right at the crossover -- see
`refs/README.md`'s own "no dead-time was modelled" note), a large
shoot-through current can still be flowing at an instant where *one*
switch's own terminal drop happens to sit outside the transition band
(because the drop reflects a voltage divider between BOTH switches, not just
this one's own conduction state) -- so some shoot-through energy gets
attributed to "conduction" for whichever switch has the smaller drop at that
instant. This does not affect the energy balance (every element's total
power is still correct and accounted for exactly once) or a switch's own
TOTAL power (conduction+transition), only how that total is split between
the two buckets. Documented here rather than hidden; see the R3/R4 numbers
below and Limits.

## Card #481: native diode and MOSFET device coverage (2026-09-29)

Task: extend `loss.py` to account native `D` (diode) and `M` (MOSFET, level=1
and BSIM4) elements, and never silently skip an element. `refs/README.md`'s
R3d/R4d section (card #479) had already found the gap directly: a trial
low-side body diode built as a native `D` gave a 0.5249% energy-balance
residual, just over the 0.5% bound, because `loss.py` had no handling for
kind `D` at all and silently dropped it.

### 1. Per-timestep device currents: the `save` fix, and its own trap

Card #470's finding ("`@device[i/p]` looks frozen at t=0 in a `tran ...
uic`") was re-investigated per this card's instruction to check whether it
was simply about `.save`. It was. Reproducer, any plain ngspice-45.2:
```
Vin vin 0 12
L1 vin swl 10u IC=0
Rdcr swl vout 10m
Rload vout 0 3.3
.control
set filetype=ascii
save @l1[p] i(l1) v(vout)
tran 5n 5u uic
wrdata out.txt @l1[p] i(l1) v(vout)
.endc
.end
```
With that `save` line (issued BEFORE `tran`), `@l1[p]` tracks correctly
across the whole run -- card #470's own reproducer (no `save` line at all)
still reproduces the freeze, confirming this was the actual cause the whole
time: ngspice's default "no `.save`/`save` at all -> keep every ordinary
circuit vector" behaviour does NOT extend to `@device[param]` references,
which are only accumulated per timestep if explicitly named in a `save`
issued before the analysis starts. Confirmed for a native diode too
(`@d1[id]`/`@d1[vd]`/`@d1[p]`, `save`d, all track correctly; `@d1[p]` was
found to equal `(V(anode)-V(cathode)) * id` -- the external terminal
voltage, including the model's own series resistance drop -- to 5+
significant figures at several sampled points, confirming `id` is the true
terminal current, not some internal-only quantity) and for a MOS1 switch
(`@m1[id/is/ig/ib/p]`, `save`d, `p` matched `id*Vd + ig*Vg + is*Vs + ib*Vb`
computed independently from the same trace, again to 5+ figures).

**Second trap, found the same day, building the R3d-native-diode fixture
below: `save` is NOT additive.** Issuing ANY `save` command REPLACES
ngspice's default "save everything ordinary" behaviour with exactly the
named list -- `save @dbody[id]` alone made even `i(l1)`/`v(vout)` (normally
free) unavailable to the later `wrdata`, failing with `Error: no such
function as i, or i(l1) is not available.` `build_wrapper()` therefore
issues one `save` line containing the FULL signal request (every `i(...)`
and `v(node)` too, not just the `@device[...]` ones) whenever any
`@device[...]` vector is needed at all.

### 2. Diode: `@name[id]`, same formula as V/L

`D` uses `@name[id]` (`save`d) exactly like a native branch current:
`P = (V(anode)-V(cathode)) * id`, bucketed into a new `diode_conduction`
category (one bucket for conduction AND reverse-recovery together -- a
diode's reverse-recovery current is just another few nanoseconds of real
current in the same `id` trace, not a separate criterion the way switch
on/off state is, so it does not need its own transition-style split the way
`--switch` devices do). `--diode NAME` (repeatable) opts a `D` device into
this bucket; an unflagged diode still gets measured (never excluded for that
reason alone) but falls into `resistive_other`.

### 3. MOSFET: why NOT `@name[id/is/ig/ib/p]`, and the sense-source method

MOS1 (level=1) exposes `is`/`ig`/`ib`/`p` and they were validated above --
but **BSIM4 (level=54) does not**. Checked directly against ngspice's own
device source (`ngspice-src/src/spicelib/devices/bsim4/b4.c`'s `BSIM4Table`
parameter list): BSIM4's OP table has `id` (`"Ids"`), `ibd`, `ibs`, `isub`,
`igs`/`igd`/`igb`/`igcs`/`igcd` (gate-leakage components) -- but no `is`, no
`ib`, no `p`. A tool reading only `@name[id]` for a BSIM4 switch would
silently miss its gate current and body/bulk current entirely, worst right
where it matters most (the fast dV/dt of a switching transition, where gate
charge and body-diode conduction are largest).

So MOSFETs use a THIRD primitive instead, one that never depends on what a
given compact model's own OP table happens to expose:
`inject_mosfet_sense()` rewrites every TOP-LEVEL `M` line in a wrapper-only
copy of the deck (never the file on disk -- same pattern `build_wrapper()`
already uses for the control block), routing each of the 4 terminals
(drain, gate, source, bulk) through an inserted 0V sense voltage source:
```
Mname nd ng ns nb model ...
```
becomes
```
Mname __lsnd_mname __lsng_mname __lsns_mname __lsnb_mname model ...
V__lsd_Mname nd __lsnd_mname DC 0
V__lsg_Mname ng __lsng_mname DC 0
V__lss_Mname ns __lsns_mname DC 0
V__lsb_Mname nb __lsnb_mname DC 0
```
and `P(t) = sum` over the 4 terminals of `V(original node) * i(sense source
there)` -- `i(Vxxx)`, the same native MNA branch-current primitive this tool
already trusted for V/L, just applied once per terminal instead of once per
2-terminal device. A 0V source adds zero impedance, so this never changes
the circuit's electrical behaviour; only `wrapper.cir` (never the deck file)
contains the rewritten lines. Only top-level `M` lines are rewritten (one
inside a `.subckt` instance is reported `excluded`, not measured -- out of
scope for this card; the diode method above works at any nesting depth,
MOSFETs do not yet).

**Validation, independent of loss.py's own `P=sum(V*I)` formula**: by
construction, ngspice's own Newton solution must satisfy KCL at each of the
4 sense sources (`id+ig+is+ib = 0`, current into the device). Checked
directly from `trace.txt` at ~200 sampled points on both MOSFET fixtures
below: worst-case `|id+ig+is+ib|` was **~1.4e-13 A** (Mhs) and **~8.2e-14 A**
(Mls) on the level=1 fixture -- machine precision, confirming the
instrumentation is wired correctly (this is now `tests/test_loss.py`'s
`KCL` check, bound at a generous 1e-6 A).

### 4. "Never silently skip an element" (device enumeration)

Before this card, `enumerate_devices()` silently `continue`d past any device
line whose type letter wasn't already in `TYPEMAP` -- `M`, `D`, `J`, an
XSPICE `A` block, etc. were invisible to the tool, not even reported under
`excluded`. It now enumerates EVERY device line (any type letter), and an
unresolved `.subckt` reference (`X` instance whose subcircuit definition
wasn't found among included files) is recorded too instead of silently
skipped. `build_plan()` decides what it can and can't compute; anything it
can't is listed under `excluded` -- and (the actual behaviour change) `main()`
now refuses to call the run a "pass" whenever `excluded` is non-empty: it
still writes `report.txt`/`summary.json` (so what WAS computed, and exactly
what wasn't, are both visible) but exits 2, not 0. See the red-proof below.

One incidental parsing bug found and fixed while building the BSIM4 fixture
(a multi-line `.model` card, `+`-continuation lines): those continuation
lines were ALSO invisible before this card's enumeration fix (their first
token, `+`, wasn't in `TYPEMAP` either, so they were silently dropped the
same way) -- once real device-line enumeration no longer silently drops
unrecognized kinds, a `+`-continuation line would have been wrongly reported
as an "unhandled device kind +". Fixed properly: `merge_continuations()`
joins every `+`-prefixed physical line into its logical predecessor before
any other pass (deck text and every `.subckt`-defining included file) sees
it, matching ordinary SPICE line-continuation semantics.

### 5. A pre-existing, SEPARATE limitation found while building the diode fixture

The task's literal first ask -- take `refs/r3d_vm_buck_deadtime.cir` and
replace its smooth `Bbody` freewheel diode with a native `D` -- was tried
(`tests/fixtures/loss_r3d_native_diode.cir`, since removed) and it ran
correctly, but it also correctly EXCLUDED `Bcmp` and `Bpwminv`. This turned
out to be nothing to do with `D`/`M` support: R3d's OWN dead-time network
(`Rchg_hs`/`Rdis_hs` between `pwm` and `dly_hs`, mirrored for `pwm_inv`)
gives node `pwm` a fanout of 3 (its driving `V=expr` B-source `Bcmp`'s own 2
terminals, PLUS the RC network's 2 more terminal-slots there) -- confirmed
directly against the UNCHANGED `refs/r3d_vm_buck_deadtime.cir` itself (not
introduced by the diode substitution):
```
>>> fan.get('pwm'), fan.get('pwm_inv')
(3, 3)
```
This trips the existing (card #470) "V-defining B-source with external
loading" exclusion -- current cannot be inferred from a `V=expr` alone when
something besides the device's own 2 terminals draws current from the node
it defines. A genuine, separate `loss.py` limitation at the time (B-source current
inference on a resistively-loaded node); fixed 2026-09-29 by reading the
B-source's own branch current `i(name)`.
`tests/fixtures/loss_diode_buck.cir` (below) sidesteps it with a topology
that never creates the situation at all: a plain non-synchronous buck, one
`I=expr` high-side switch (never subject to that check) and a native diode
freewheel path, no dead-time steering network needed.

### 6. New fixtures and results (real output, ngspice-45.2, 2026-09-29)

All three settled-`.tran` path, last 5 periods (10 us) of a 1.5ms run,
`--period 2e-6 --periods 5 --tstop 1.5e-3 --dcr Rdcr --esr Resr`:

| | `loss_diode_buck.cir` (native `D`) | `loss_mos1_buck.cir` (level=1 `M`) | `loss_bsim4_buck.cir` (BSIM4/level=54 `M`) |
|---|---:|---:|---:|
| topology | non-sync buck: 1 `I=expr` switch (`Bhs`) + native `D` freewheel (`Ddiode`) | sync buck: PMOS `Mhs` + NMOS `Mls`, hand-timed dead-time | same as level=1, BSIM4 models |
| P_in (Vin) | 3.768472 W | 3.242987 W | 3.213927 W |
| P_out (Rload) | 3.196314 W | 3.165508 W | 3.106597 W |
| efficiency | 84.82 % | 97.61 % | 96.66 % |
| switch_conduction | 0.004910 W | 0.043289 W | 0.058003 W |
| switch_transition | 0.003718 W | 0.007791 W | 0.040925 W |
| diode_conduction | 0.572572 W | 0 (n/a) | 0 (n/a) |
| inductor_dcr | 0.009863 W | 0.009876 W | 0.009603 W |
| capacitor_esr | 0.000259 W | 0.000188 W | 0.000187 W |
| resistive_other | -0.015154 W | 0.013908 W | -0.001450 W |
| **residual** | **-0.1064 % of P_in** | **+0.0748 % of P_in** | **+0.0019 % of P_in** |
| excluded devices | none | none | none |
| exit code | 0 | 0 | 0 |

All three comfortably inside the 0.5% bound, and all exit 0 (no excluded
devices). `resistive_other`'s small negative values are `L1`/`Cout`'s own
directly-measured average power over the window (not assumed zero -- same
"why `loss.py` measures L/C directly" reasoning as R3/R4's own red-proof 2),
small but not exactly zero since the window is not a perfectly steady-state
integer-period slice of an infinitely-settled waveform.

**Red-proof (`tests/fixtures/loss_unsupported_jfet.cir`)**: a small buck-like
circuit, otherwise fully supported, with one native JFET (`J1`, drain=vout,
gate=source=0 -- a real, non-negligible current load, not a degenerate
no-op) added in parallel with the load. Required and confirmed behaviour:
`J1` is enumerated, listed under `excluded` (`"unhandled device kind J"`),
and the run exits 2 -- even though the naive/untrusted residual computed
from the incomplete sum (+0.1963% of P_in) would, on its own, have looked
like a clean pass well inside the 0.5% bound. This is the point of the
red-proof: a small-looking residual is not sufficient evidence of a
complete sum, which is exactly why `excluded` now gates the pass/fail
decision rather than the residual number alone.

## R3 and R4 results (real output, ngspice-45.2, 2026-09-29)

Settled-`.tran` path, last 5 periods (10 us) of the decks' own established
runs (R3: `tran ... 1.5e-3 0 uic`; R4: `tran ... 6e-3 0 uic`, both settled
per `refs/README.md`), `--switch-low 0.10 --switch-high 0.90` (defaults):

| | R3 (voltage-mode) | R4 (peak-current-mode) |
|---|---:|---:|
| P_in (Vin) | 11.520161 W | 11.395223 W |
| P_out (Rload) | 3.299971 W | 3.300037 W |
| efficiency | 28.65 % | 28.96 % |
| switch_conduction | 2.266604 W | 1.616473 W |
| switch_transition | 5.933799 W | 6.414134 W |
| gate_drive | 0 (not modelled) | 0 (not modelled) |
| inductor_dcr | 0.010197 W | 0.010202 W |
| capacitor_esr | 0.000191 W | 0.000181 W |
| resistive_other | 0.002715 W | 0.054303 W |
| sum of losses | 8.213505 W | 8.095293 W |
| **residual** | **+0.0580 % of P_in** | **-0.0009 % of P_in** |
| wall time | 13.8 s | 66.3 s |

Both comfortably inside the 0.5% target. R4's `resistive_other` is larger
than R3's mostly because of its current-sense resistor `Rsense` (0.05 ohm,
carrying the full load current, ~1A -> ~0.05W, matches).

**The headline number here is efficiency ~29%, not ~90%+.** This is real,
not a bug: R3/R4's smooth (tanh) switches have no dead-time between Bhs and
Bls, so at every transition both switches are significantly conductive at
once, producing a genuine shoot-through current spike (peaking in the
hundreds of amps momentarily, confirmed against the native, ngspice-solved
`i(Vin)` -- not an artifact of this tool's own current computation, see
Validation below). `switch_transition` dominates both decks' loss budget.
`refs/README.md`'s own "what was not done" section already flagged the
missing dead-time as a modelling gap; this tool is the first to quantify its
cost. The DCR/ESR/total-switch-power numbers -- the genuinely small,
well-behaved quantities -- validate cleanly (next section); it is the
conduction/transition SPLIT of that large transition loss, on this specific
no-dead-time topology, that carries the "honest limitation" noted above.

PSS-shoot orbit path (fork build `aaba50ba3f49-pss9`, R3, one converged
period, `Fguess=500kHz StabTime=300us Points=200`; PSS-SHOOT reported
`CONVERGED iter=4 residual=4.1e-09 spectral_radius=0.9803 stable`, matching
`patches-wip/pss/README.md`'s own R3 result): P_in=11.201866 W, P_out=3.300037
W, efficiency=29.46%, switch_conduction=1.702530 W, switch_transition=
6.185882 W, inductor_dcr=0.010213 W, capacitor_esr=0.000191 W,
resistive_other=0.003000 W, **residual = 0.00011% of P_in** -- even tighter
than the tran path, and consistent with it (same shoot-through-dominated
picture). This run is real (this exact wrapper netlist, run directly) but
**the `.pss method=shoot` analysis itself was observed to be non-
deterministic**: the identical wrapper failed at
`PSS-SHOOT: [stage 1/probe] monodromy perturbation run for state C:c2f (+h)
did not complete ... (ngspice error 106)` on other attempts, with no change
to the netlist or command. `patches-wip/pss/README.md` already documents
real Newton-stability sensitivity in this WIP analysis (breakpoint forcing,
`eps`/`reltol` tuning, "Newton stall/fail" cases called out explicitly), so
this is treated as a property of the fork build under test, not of
`loss.py` -- `tests/test_loss.py` attempts it once, best-effort, and treats
a convergence failure as COULD-NOT-RUN (non-gating), never FAIL.

## Validation (tests/test_loss.py)

- **Energy balance <= 0.5% of P_in**, R3 and R4, settled-`.tran` path:
  **pass** (+0.058%, -0.0009%; see table above).
- **DCR loss vs. hand calc** `Rdcr * mean(i_L1^2)`, using the *native*
  `i(L1)` (independent of this tool's own `(V(swl)-V(vout))^2/Rdcr`
  computation for the same element): R3 rel_err 0.0000%, R4 rel_err 0.0000%
  (both pass, well inside 1%).
- **ESR loss vs. hand calc** `Resr * mean(i_C^2)`, `i_C = Cout*dV(cx)/dt`
  from the capacitor's own two terminals (independent of this tool's own
  `(V(vout)-V(cx))^2/Resr`): R3 rel_err 0.0198%, R4 rel_err 0.0031% (pass).
- **Bhs total dissipated power vs. hand calc** `(V(vin)-V(sw)) * (-i(vin))`,
  using ngspice's own Newton-solved native source current -- independent of
  this tool's Python re-evaluation of Bhs's tanh conductance expression, and
  valid because Bhs is the only device at node `vin` in both decks (checked):
  R3 rel_err 0.0455%, R4 rel_err 0.0088% (pass). This is the strongest single
  check in the suite: it confirms the B-expression evaluator (parsing,
  `.param` substitution, `V(...)` rewriting) reproduces what ngspice itself
  solved, to within numerical/windowing noise, not just that the tool is
  internally consistent.
- **Red-proof 1 (omit one element)**: on R3, dropping Bhs (~35% of P_in, not
  Rdcr -- Rdcr is only ~0.1% of P_in on this shoot-through-dominated deck and
  does not by itself move the residual past 0.5%, found while building this
  suite) from the sum: residual jumps from +0.058% to +35.68% -- **pass**.
- **Red-proof 2 (non-integer-period window)**: `loss.py`'s own balance does
  NOT break for a non-integer window (R3, 1.5 periods: residual +0.099%,
  still inside 0.5%) -- **and that is correct**, because this tool measures
  every element (including every L/C) directly from the trace every time,
  rather than assuming a reactive element's long-run average is ~0 (which is
  only true over an exact integer number of periods). The red-proof is
  demonstrated against the failure mode it actually targets instead: with
  every L/C device DROPPED from the sum (simulating a tool that made that
  assumption), the integer-period window still balances fine (+0.053%,
  because the assumption happens to be true there) but the 1.5-period window
  does not (-14.47%, because it is not) -- **pass**. This is why `loss.py`
  measures L/C directly rather than assuming.
- **Sanity**: efficiency in (0,1) for both decks (pass, with the honest note
  above); DCR loss R3 vs R4 within 0.06% (shared power stage, expected very
  close: pass); ESR loss within 5.3% (same power stage, some spread expected
  from the different control schemes' ripple shape: pass against a 15%
  bound); total switch loss (conduction+transition) R3 vs R4 ratio 1.02
  (pass, same switches).
- **PSS-orbit path**: attempted once on R3, best-effort/non-gating (see
  above) -- COULD-NOT-RUN in the recorded run, PASS in an earlier manual run
  with the identical wrapper (residual 0.00011%).
- **Card #481 (native D/M device coverage)**:
  - `loss_diode_buck.cir` (native `D` freewheel diode): energy balance
    **pass** (-0.1064%), `diode_conduction` genuinely nonzero **pass**
    (0.572572 W).
  - `loss_mos1_buck.cir` (level=1 MOSFET switches): energy balance **pass**
    (+0.0748%); independent KCL check at each MOSFET's own 4 sense sources
    (`id+ig+is+ib ~= 0`, nothing to do with loss.py's own power formula):
    Mhs **pass** (max residual 1.475e-13 A), Mls **pass** (8.081e-14 A).
  - `loss_bsim4_buck.cir` (BSIM4/level=54 MOSFET switches): energy balance
    **pass** (+0.0019%); KCL: Mhs **pass** (8.349e-14 A), Mls **pass**
    (9.284e-15 A).
  - Red-proof (`loss_unsupported_jfet.cir`, a native JFET `J1` with no power
    method): exit code is 2, never 0 -- **pass**; `J1` listed under
    `excluded` -- **pass**; recorded as non-gating info that the naive
    residual alone (+0.1963%) would have looked like a pass, which is why
    `excluded` (not the residual) has to gate.

Run: `python3 tests/test_loss.py` (needs `ngspice` on PATH, or `NGSPICE=...`;
the PSS check additionally wants `NGSPICE_AUG=<fork>/bin/ngspice-aug`, default
`~/Work/apps/ngspice_augmentation/aaba50ba3f49-pss9/bin/ngspice-aug`). Exit 0
= every gating check passed (the recorded run, card #481: **30 pass, 0 FAIL,
0 COULD-NOT-RUN gating, 2 non-gating info** -- the original 17 R3/R4 checks
unchanged plus 13 new checks for D/M coverage and the red-proof).

## Limits / what was not done

- **Device coverage**: R, L, C, V, I (constant-0 only), B (`I=`/`V=`), G/E
  (simple linear 2-terminal-output form only, `Gname n1 n2 nc1 nc2 gain`), D
  (native diode, any model level, card #481) and M (native MOSFET, any
  compact-model level, TOP-LEVEL lines only, card #481), and one level of
  `.subckt` expansion (recursive, but only for decks built from R/L/C/V/I/
  B/G/E/D -- enough for `stb/lgprobe.lib`, R3 and R4; a MOSFET inside a
  `.subckt` instance is enumerated and reported `excluded`, never silently
  measured wrong or dropped -- see "Card #481"). Still no BJT compact
  models, no `POLY(...)` B/G/E forms, no time-varying independent current
  sources, and no XSPICE code-model blocks. Card #481's own hard requirement
  is that none of THIS is silent any more: any device kind/situation this
  tool cannot account is enumerated, listed under `excluded`, and blocks a
  "pass" (exit 2) -- see the JFET red-proof above.
- **V-defining B-source with external resistive loading** (fixed 2026-09-29):
  such a `B, V=expr` device is now measured by ngspice's own branch current
  `i(name)` (a V-defining B-source has a branch, same sign convention as a V
  source), so R3d/R4d's `Bcmp`/`Bpwminv` (`Bqbuf`/`Bqinv`) are accounted and
  their energy balance runs. Before this fix they were excluded (card #470/#481).
- **MOSFET fixtures are open-loop** (`loss_mos1_buck.cir`, `loss_bsim4_buck.cir`
  -- fixed duty, hand-timed dead-time PULSE gate drives, no compensator):
  built to validate `loss.py`'s device accounting, not as converter-design
  references the way R3/R4/R3d/R4d are. Only two model levels were tried
  (level=1, BSIM4/level=54, per the task); nothing else in ngspice's MOSFET
  family (VDMOS, BSIM3/level=8 or 9, BSIM6, ...) was checked, though the
  sense-source method does not depend on model level by construction.
- **`diode_conduction` is one bucket**, not split into separate conduction
  and reverse-recovery categories the way `--switch` devices split into
  conduction/transition: a diode's reverse-recovery current is just a few
  more nanoseconds inside the same measured `id` trace, not a distinct
  time-domain criterion like a switch's own on/off state, so it did not seem
  to warrant the same kind of split. If a specific reverse-recovery number
  turns out to be needed later, the switch-style rail-fraction split
  wouldn't directly apply (a diode has no analogous "blocking voltage
  fraction" state machine) and would need its own criterion.
- **The conduction/transition split** on a no-dead-time, shoot-through-heavy
  smooth-switch deck like R3/R4 has the attribution caveat documented above
  (total switch power is exact; the split between the two buckets is not
  perfectly clean for this specific topology). A deck WITH real dead-time
  (e.g. a future R3/R4 variant, or R6's ideal-`sw` construction) would not
  have this issue, since there would be no simultaneous-conduction interval
  to misattribute.
- **PSS-orbit path**: implemented and shown to work (residual 0.00011% on a
  successful R3 run), but not included as a gating/required test because the
  fork's `.pss method=shoot` analysis itself was observed to be
  non-deterministic in this environment (see above) -- consistent with, not
  contradicting, this project's own `patches-wip/pss/README.md`. R4's PSS
  path was not attempted at all (time budget went to nailing down the tran
  path and the `@device[i]/[p]` freeze bug, which is the load-bearing finding
  for this whole tool); nothing suggests it would behave differently from R3.
- **No independent re-verification (V-G4-style)**: every number in this
  report and in `tests/test_loss.py` was produced and read by this same
  card; a second pass by another reviewer has not happened.
- **CPU/memory**: the first working version of this tool read the *entire*
  settled-`.tran` trace into Python memory regardless of the requested
  window (R4's ~6ms/2ns-step run: several GB, became the dominant cost, not
  ngspice) -- found and fixed in this same card by streaming the trace and
  keeping only the requested window (`read_wrdata_window`); the numbers in
  this doc are from the fixed version (R3 14s/360MB peak, R4 66s/1.7GB peak).
- `--tstep` defaults to `period/1000`; finer than the decks' own `10n`
  default step but coarser than an earlier `period/2000` attempt that was
  needlessly slow/memory-heavy for no meaningful accuracy gain (checked: R3's
  residual moved from 0.024% to 0.058% between those two step sizes, both
  comfortably inside the 0.5% bound).
