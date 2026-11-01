# Bug 2 (card #539, robustness, not silent): `.tran` aborts "Timestep too small ... trouble with node vramp#branch" at tight tolerances on R7

`repro.cir` (R7 deck body + `option reltol=1e-6 abstol=1e-9`, `tran 1n 60u 0 uic`, engine s88): aborts at t = 4 us,
`Timestep too small; time = 4e-06, timestep = 1.25e-21: trouble with node "vramp#branch"`. Variants: reltol 1e-6 +
vntol 1e-7 aborts at 12 us; reltol 1e-6 + abstol 1e-9 + vntol 1e-7 at 26 us; reltol 1e-3 + the same abstol/vntol runs; reltol 1e-6
alone (as `.options` card or control `option`) also ran 60 us in other trials, so the failure depends on the combination
and is not monotonic (a LOOSER abstol makes it fail). The unknown is the branch current of the ramp source, which only
charges one BSIM4 gate (fA): the Newton convergence test on it is dominated by model noise at these tolerances.
Consequence: reltol <= 1e-5 settled-tran references need `.options` only with reltol alone (what `tests/fixtures/r7_tran.cir`
uses, 1e-5) and `tools/pac_perturbed.py`'s default `--reltol 1e-5` aborted on R7 at 42 us (use `--reltol 1e-4`, validated
here against FRA to 0.004 dB). Not diagnosed further.

RESOLVED (card #541, 2026-10-02): a deck bug, not an engine one: the ramp PULSE has tr+pw+tf = per+1n; stock ngspice 45.2 and the unpatched np1 build abort identically. See refs/README.md R7 "Engine findings" item 2.
