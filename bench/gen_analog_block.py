#!/usr/bin/env python3
"""gen_analog_block.py -- generate bench/decks/d_analog_block.cir
(ngspice_augmentation card #467, plan T9.1, bench task 1d).

A synthetic analog block with many BSIM4 devices (200-500 transistors, per
the task): N independent 5-transistor OTA cells (differential pair + PMOS
current-mirror load + NMOS tail current source -- a standard textbook OTA,
not any real design) sharing global Vdd/bias/ground rails, each driven by
its own small differential step so every cell does real transient work.
Default N=60 -> 300 BSIM4 instances (60 NMOS pairs + 60 PMOS mirrors + 60
NMOS tails = 3*60 = 180... corrected below to 5*N).

This is deliberately a *parallel array* of small independent cells rather
than one giant coupled circuit: it gives many BSIM4 evaluations per Newton
iteration (the property the task asks to stress -- "device evaluation
dominates") without the topology risk of a single large multi-stage
op-amp/bandgap converging. Not derived from any real op-amp/bandgap design
(#454 A: synthetic only).

Output: bench/decks/d_analog_block.cir (checked in; regenerate with this
script). Simulation output goes to
~/Work/data/ngspice_augmentation/accel/bench/d_analog_block/, never into
the repo.

--device {bsim4,bsim3} (card #483 T9.3 task 2): the same 5T-OTA topology,
same tstop/tstep/meas structure, with the .model cards (and only the
.model cards -- topology, sizes, sources are unchanged) swapped for a
different device model, so the same benchmark shape can exercise a
different OpenMP device loader (b3ld.c's #pragma omp parallel for,
candidate 0105) instead of BSIM4's. BSIM3 model card
(level=49 version=3.3.0) uses the same minimal-parameter style as this
fork's own examples/svg/nmos_out_BSIM330_ps.sp (a checked-in, known-good
BSIM3v3.3.0 invocation) -- tox/nch/nsub/vth0/u0/vsat only, not a full
foundry-style parameter deck, same spirit as the BSIM4 model cards below.
"""
import argparse
import pathlib

HERE = pathlib.Path(__file__).resolve().parent

CELL = """\
* --- cell {n}: 5T OTA (M1/M2 diff pair, M3/M4 PMOS mirror load, M5 NMOS tail) ---
Vinp{n} inp{n} 0 PULSE({vcm} {vhi:.6f} 200n 5n 5n 800n 2000n)
M1_{n} outn{n} inp{n} tail{n} 0 {nmod} w=4u l=0.5u
M2_{n} out{n} vcm tail{n} 0 {nmod} w=4u l=0.5u
M3_{n} outn{n} outn{n} vdd vdd {pmod} w=8u l=0.5u
M4_{n} out{n} outn{n} vdd vdd {pmod} w=8u l=0.5u
M5_{n} tail{n} vbias_tail 0 0 {nmod} w=8u l=1u
Cload{n} out{n} 0 200f
"""

MODELS = {
    "bsim4": {
        "nmod": "nch4",
        "pmod": "pch4",
        "cards": """\
.model nch4 nmos level=54 version=4.8.2
+ mobmod=0 capmod=0 igcmod=0 igbmod=0 rdsmod=0 rbodymod=0 rgatemod=0 acnqsmod=0
+ toxe=3e-9 toxp=3e-9 toxm=3e-9
+ vth0=0.5 k1=0.4 k2=0 k3=0 w0=2.5e-6
+ u0=300 vsat=1.1e5
+ cgso=2e-10 cgdo=2e-10 cgbo=0 cgsl=0 cgdl=0
+ cf=0 cjs=1e-4 cjd=1e-4 cjswd=1e-10 cjswgd=1e-10

.model pch4 pmos level=54 version=4.8.2
+ mobmod=0 capmod=0 igcmod=0 igbmod=0 rdsmod=0 rbodymod=0 rgatemod=0 acnqsmod=0
+ toxe=3e-9 toxp=3e-9 toxm=3e-9
+ vth0=-0.5 k1=0.4 k2=0 k3=0 w0=2.5e-6
+ u0=100 vsat=1.0e5
+ cgso=2e-10 cgdo=2e-10 cgbo=0 cgsl=0 cgdl=0
+ cf=0 cjs=1e-4 cjd=1e-4 cjswd=1e-10 cjswgd=1e-10
""",
    },
    # card #483 T9.3 task 2: exercises b3ld.c's OMP load loop (candidate
    # 0105) instead of b4ld.c's. Minimal BSIM3v3.3.0 parameters, same style
    # as examples/svg/nmos_out_BSIM330_ps.sp (checked-in, known-good).
    "bsim3": {
        "nmod": "nch3",
        "pmod": "pch3",
        "cards": """\
.model nch3 nmos level=49 version=3.3.0
+ tox=3.5e-9 nch=2.4e17 nsub=5e16
+ vth0=0.5 u0=300 vsat=1.1e5

.model pch3 pmos level=49 version=3.3.0
+ tox=3.5e-9 nch=2.5e17 nsub=5e16
+ vth0=-0.5 u0=100 vsat=1.0e5
""",
    },
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-cells", type=int, default=60)
    ap.add_argument("--tstop-ns", type=float, default=2000.0)
    ap.add_argument("--tstep", default="0.5n")
    ap.add_argument("--device", choices=sorted(MODELS), default="bsim4")
    ap.add_argument("--out", default=str(HERE / "decks" / "d_analog_block.cir"))
    args = ap.parse_args()
    n = args.n_cells
    ntrans = 5 * n
    deck_name = pathlib.Path(args.out).stem
    dev = MODELS[args.device]
    devlabel = args.device.upper()

    L = [f"""\
* {deck_name} -- {n} independent 5T-OTA cells = {ntrans} {devlabel} instances
* (ngspice_augmentation card #467/#478/#483, plan T9.1/T9.2/T9.3; generated
* by bench/gen_analog_block.py --tstop-ns {args.tstop_ns} --device {args.device}
* -- do not hand-edit, regenerate instead). Synthetic textbook OTA array,
* synthetic, not derived from any real design.
* run: ngspice -b <this file> (writes to
*  ~/Work/data/ngspice_augmentation/accel/bench/{deck_name}/, never
*  into the repo)

.param vdd_v=3.3 vcm=1.65

{dev['cards']}
Vdd vdd 0 DC {{vdd_v}}
Vcm vcm 0 DC {{vcm}}
Vbias_tail vbias_tail 0 DC 1.0
Rbias vbias_tail 0 1Meg

"""]

    for k in range(n):
        vhi = 1.65 + 0.05 + 0.0002 * k
        L.append(CELL.format(n=k, vcm="{vcm}", vhi=vhi, nmod=dev["nmod"], pmod=dev["pmod"]))

    tstop_ns = args.tstop_ns
    t_from = tstop_ns * 0.75
    t_mid = tstop_ns * 0.5
    t_end = tstop_ns * 0.995

    meas = []
    meas.append(f"\n.tran {args.tstep} {tstop_ns}n 0\n.options reltol=1e-3\n\n.control\nset filetype=ascii\nrusage cputime\nrun\nrusage cputime\n")
    for k in (0, n // 2, n - 1):
        meas.append(f"meas tran wave_out{k}_mid find v(out{k}) at={t_mid}n\n")
        meas.append(f"meas tran wave_out{k}_end find v(out{k}) at={t_end}n\n")
    meas.append(f"meas tran wave_out0_settle AVG v(out0) FROM={t_from}n TO={tstop_ns}n\n")
    meas.append("print wave_out0_mid wave_out0_end\n")
    meas.append(f"wrdata $HOME/Work/data/ngspice_augmentation/accel/bench/{deck_name}/wave.txt v(out0) v(out{n//2}) v(out{n-1})\n")
    meas.append("rusage cputime\nquit\n.endc\n.end\n")
    L.append("".join(meas))

    outpath = pathlib.Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text("".join(L))
    print(f"wrote {outpath} ({n} cells, {ntrans} {devlabel} instances, tstop={tstop_ns}ns)")


if __name__ == "__main__":
    main()
