# Bug 3 (card #539): `.pxf adjoint=1` gives ~0 for a source that sets a MOS threshold crossing; `adjoint=0` and `.pac` agree with each other

`python3 run_repro.py` (`run_repro.out`): inverter chain driven by ramp + DC source Vc, `.pxf ... sources="Vc" adjoint=0|1`
at 1 kHz. Truth (static derivative of the periodic state): ~0.5.

```
level1_cg0     N=1  forward 0.5325  adjoint 0.0194      <- MOS1, NO capacitance anywhere: bug 3 is not bug 1
level1_cg0     N=2  forward 0.4677  adjoint 0.0000
bsim4          N=1  forward 0.5546  adjoint 0.0161
```
On R7 (`../bug1_pac_model_charge/r7_pxf_fwd_vs_adjoint.txt`): Iinj agrees (0.01-0.27 dB), Vctrl differs 0.5 dB (both ~1e-6, meaningless),
Vin differs by a constant 24 dB at 5/20/200 kHz (adjoint lower). `PAC_EXACT_JACOBIAN=1` changes nothing. Not diagnosed; the
adjoint engine (`pac_engine_run_adjoint`) was validated on smooth B-source decks only (docs/pss-pac.md card #504), where no
node sets a timing through a MOS gate. Consequence for users: forward-vs-adjoint "consistency" cannot be used as evidence on
these decks; here it fails, and the forward side is itself wrong on R7 (bug 1).

RESOLVED (card #541, 2026-10-02): see refs/README.md R7 "Engine findings" item 3 and docs/pss-pac.md "Native path: J_k stored". With the fix: level1_cg0 N=1 forward 0.5325 adjoint 0.5325; N=2 0.4677 / 0.4677. Run from a copy outside the repo (the script writes ./work).
