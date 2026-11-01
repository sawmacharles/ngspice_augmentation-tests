#!/usr/bin/env python3
"""gen_trancorpus -- writes the transient benchmark corpus tests/fixtures/tran/*.cir (card #641).

Every deck is self-contained (built-in models only), small (stock run well under ~20 s on one
core), and starts with
    * <name> -- <what it stresses>
    * stock: <outcome of the stock ngspice-559 binary, measured by tranbench>
    * probe: <vectors tranbench compares>
Decks marked FAIL below were verified to make the stock binary fail ('timestep too small' /
no convergence); the note is kept in STOCK and is the only thing that needs updating when the
base binary changes.  Run:  python3 tools/gen_trancorpus.py [outdir]
"""
import os
import sys

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests", "fixtures", "tran")

# measured with tools/tranbench.py on ngspice-src-559/release/src/ngspice, label base_stock (2026-10-15)
STOCK = {   # only the decks the stock binary does NOT complete; every other deck is noted as 'completes'
    "tr_chatter_bsrc": "FAIL: 'Transient op failed, timestep too small' (the operating point sits on the comparator discontinuity)",
    "tr_chatter_bsrc_uic": "FAIL: timestep too small after ~136 points (0.7% of the run): sliding-mode chatter",
    "tr_float_cap_switch": "FAIL: no convergence at the operating point (singular matrix: nodes x, y have no DC path)",
    "tr_srlatch_bsrc": "FAIL: 'Transient op failed, timestep too small' (algebraic loop through two discontinuous B-sources)",
    "tr_sw_nohyst_osc": "FAIL: timestep too small at ~5% of the run (switch controlled by its own node, VH=0)",
    "tr_sw_selfctl": "FAIL: timestep too small at ~2.5% of the run (switch controlled by its own node, VH=0)",
}

MODELS = {
    "d": ".model dsch D(IS=1e-8 N=1.05 RS=0.05 CJO=50p TT=5n BV=40 IBV=1m)",
    "sw": ".model swm SW(RON=10m ROFF=10Meg VT=2.5 VH=0.5)",
    "n1": ".model nm1 NMOS(LEVEL=1 VTO=0.7 KP=120u LAMBDA=0.03 CGSO=0.3n CGDO=0.3n CBD=10f CBS=10f PB=0.8)",
    "p1": ".model pm1 PMOS(LEVEL=1 VTO=-0.7 KP=50u LAMBDA=0.03 CGSO=0.3n CGDO=0.3n CBD=10f CBS=10f PB=0.8)",
    "npn": ".model qn NPN(IS=1e-15 BF=120 BR=3 VAF=80 CJE=1p CJC=0.8p TF=0.3n TR=5n RB=50 RC=5 RE=1)",
    "pnp": ".model qp PNP(IS=1e-15 BF=60 BR=2 VAF=50 CJE=1p CJC=0.8p TF=0.5n TR=10n RB=50)",
    "n3": ".model nm3 NMOS(LEVEL={lvl} VERSION=3.3.0 TOX=7.5n VTH0=0.5 U0=400 VSAT=1.2e5 CGSO=2e-10 CGDO=2e-10 CJ=1e-3 CJSW=2e-10 MJ=0.5 MJSW=0.33 PB=0.9 PBSW=0.9 XJ=0.15u NCH=1.7e17)",
    "p3": ".model pm3 PMOS(LEVEL={lvl} VERSION=3.3.0 TOX=7.5n VTH0=-0.5 U0=130 VSAT=1.0e5 CGSO=2e-10 CGDO=2e-10 CJ=1e-3 CJSW=2e-10 MJ=0.5 MJSW=0.33 PB=0.9 PBSW=0.9 XJ=0.15u NCH=1.7e17)",
}

D = {}   # name -> (stress description, probe, body)


def deck(name, stress, probe, body):
    D[name] = (stress, probe, body.strip("\n") + "\n")


# ---------------------------------------------------------------- linear
deck("tr_rc_step", "RC low-pass driven by a pulse: linear sanity deck, exact answer 1-exp(-t/RC) (tau 1 ms)", "v(out) v(in)", """
V1 in 0 PULSE(0 1 100u 1n 1n 2m 4m)
R1 in out 1k
C1 out 0 1u
.meas tran v_at_1m FIND v(out) AT=1.1m
.meas tran v_at_end FIND v(out) AT=4.9m
.tran 5u 5m
""")

deck("tr_rlc_underdamped", "series RLC step response, Q ~ 3.2, f0 5 kHz: ringing with exact solution; checks damping of the integration method", "v(c) l1#branch", """
V1 in 0 PULSE(0 1 20u 1n 1n 1 2)
R1 in a 10
L1 a c 1m IC=0
C1 c 0 1u IC=0
.meas tran vpk MAX v(c) FROM=0 TO=2m
.tran 2u 2m
""")

deck("tr_lc_highq", "nearly lossless parallel LC, 100 cycles of ringing (trapezoidal damping/phase error); step current kick", "v(n)", """
I1 0 n PULSE(0 1m 1u 1n 1n 10u 1)
L1 n 0 1m
C1 n 0 1u
R1 n 0 1Meg
.meas tran vpk MAX v(n) FROM=18m TO=20m
.tran 1u 20m
""")

deck("tr_rc_sine", "sine steady state through an RC low-pass (fc 159 Hz, 1 kHz drive), 20 cycles: phase and amplitude accuracy at the default tolerance", "v(out)", """
V1 in 0 SIN(0 1 1k)
R1 in out 1k
C1 out 0 1u
.meas tran vmax MAX v(out) FROM=15m TO=20m
.tran 2u 20m
""")

deck("tr_rc_ladder20", "20-section RC ladder (diffusion line) step response: many time constants spread over a decade or two, sparse linear system", "v(n20) v(n10)", ("""
V1 in 0 PULSE(0 1 1u 10n 10n 1 2)
R0 in n1 100
C0 n1 0 100n
""" + "".join(f"R{i} n{i} n{i+1} 100\nC{i} n{i+1} 0 {int(100*(1.0+0.05*i))}n\n" for i in range(1, 20)) + """
.tran 100n 2m
"""))

# ---------------------------------------------------------------- rectifiers / clamps
deck("tr_bridge_rect", "diode bridge rectifier with a big smoothing C (15 V peak, 50 Hz, 1 mF, 100 ohm): stiff charging pulses every half cycle", "v(out) v(a)", MODELS["d"] + """
Vs a b SIN(0 15 50)
Rga a 0 1Meg
Rgb b 0 1Meg
D1 a out dsch
D2 b out dsch
D3 0 a dsch
D4 0 b dsch
Cb out 0 1m
RL out 0 100
.meas tran vout_end FIND v(out) AT=95m
.tran 20u 100m
""")

deck("tr_halfwave_ideal", "half-wave rectifier with an almost ideal diode (N=0.02 -> exp(v/0.5 mV), RS=1 uOhm): extreme exponential, Newton limiting", "v(out) v(in)", """
.model didl D(IS=1e-14 N=0.02 RS=1u)
V1 in 0 SIN(0 5 1k)
D1 in out didl
RL out 0 1k
.meas tran vpk MAX v(out) FROM=1m TO=3m
.tran 1u 3m
""")

deck("tr_zener_clamp", "avalanche/zener clamp: 20 V sine through 1 kOhm into a 6.8 V BV diode in parallel with a signal diode (reverse breakdown knee)", "v(clamp) v(in)", """
.model dz D(IS=1e-14 N=1 RS=2 BV=6.8 IBV=5m NBV=1.5 CJO=30p TT=2n)
V1 in 0 SIN(0 20 1k)
R1 in clamp 1k
Dz 0 clamp dz
Dc clamp 0 dz
.meas tran vmax MAX v(clamp) FROM=0 TO=3m
.meas tran vmin MIN v(clamp) FROM=0 TO=3m
.tran 1u 3m
""")

# ---------------------------------------------------------------- converters
deck("tr_buck_sw", "buck converter, SW-model high-side switch driven by PULSE (100 kHz, D=0.4, 12 V -> 5 V), freewheel Schottky-like diode, 100 cycles", "v(out) l1#branch", MODELS["d"] + "\n" + MODELS["sw"] + """
Vin in 0 DC 12
Vg g 0 PULSE(0 5 0 20n 20n 4u 10u)
S1 in sw g 0 swm
D1 0 sw dsch
L1 sw o 22u IC=0
Resr o out 20m
C1 out 0 47u IC=0
RL out 0 5
.meas tran vavg AVG v(out) FROM=0.9m TO=1m
.tran 50n 1m uic
""")

deck("tr_buck_mos", "buck converter with an NMOS (level 1) high-side switch and floating gate drive, MOS capacitances, diode freewheel", "v(out) l1#branch", MODELS["d"] + "\n" + MODELS["n1"] + """
Vin in 0 DC 12
Vdrv g sw PULSE(0 10 0 10n 10n 4u 10u)
M1 in g sw 0 nm1 W=20000u L=1u
D1 0 sw dsch
L1 sw o 22u IC=0
Resr o out 20m
C1 out 0 47u IC=0
RL out 0 5
.meas tran vavg AVG v(out) FROM=0.9m TO=1m
.tran 50n 1m uic
""")

deck("tr_boost_sw", "boost converter, SW-model low-side switch, 5 V -> ~9 V open loop, 100 kHz, diode rectifier, 100 cycles", "v(out) l1#branch", MODELS["d"] + "\n" + MODELS["sw"] + """
Vin in 0 DC 5
Vg g 0 PULSE(0 5 0 20n 20n 4u 10u)
L1 in sw 22u IC=0
S1 sw 0 g 0 swm
D1 sw out dsch
C1 out 0 47u IC=0
RL out 0 25
.meas tran vavg AVG v(out) FROM=0.9m TO=1m
.tran 50n 1m uic
""")

deck("tr_boost_mos", "boost converter with a grounded NMOS (level 1) switch driven by a 0-10 V gate pulse; drain ringing with the diode capacitance", "v(out) l1#branch", MODELS["d"] + "\n" + MODELS["n1"] + """
Vin in 0 DC 5
Vg g 0 PULSE(0 10 0 10n 10n 4u 10u)
L1 in sw 22u IC=0
M1 sw g 0 0 nm1 W=20000u L=1u
D1 sw out dsch
C1 out 0 47u IC=0
RL out 0 25
.meas tran vavg AVG v(out) FROM=0.9m TO=1m
.tran 50n 1m uic
""")

deck("tr_ind_kick_clamp", "inductive kick: 10 mH relay-like load carrying 1 A is interrupted by a MOS-sized switch; diode + 20 V zener clamp across the load", "v(sw) l1#branch", """
.model dcl D(IS=1e-12 N=1.1 RS=0.1 BV=20 IBV=10m CJO=20p TT=20n)
.model swk SW(RON=50m ROFF=100Meg VT=2.5 VH=0.1)
V1 vcc 0 DC 12
Vg g 0 PULSE(5 0 20u 100n 100n 60u 200u)
L1 vcc a 10m IC=0.24
Rl a sw 50m
Dc a vcc dcl
S1 sw 0 g 0 swk
Csw sw 0 100p
.tran 50n 150u uic
""")

# ---------------------------------------------------------------- oscillators
deck("tr_bjt_astable", "BJT astable multivibrator (2N2222-like, 5 V, Rb 47k, C 1 nF, ~66 us period), 8 periods: hard regenerative edges; .ic breaks the symmetry", "v(c1) v(c2) v(b1)", MODELS["npn"] + """
Vcc vcc 0 DC 5
Rc1 vcc c1 1k
Rc2 vcc c2 1k
Rb1 vcc b1 47k
Rb2 vcc b2 47k
C1 c1 b2 1n
C2 c2 b1 1n
Q1 c1 b1 0 qn
Q2 c2 b2 0 qn
.ic v(c1)=0.2 v(c2)=5 v(b1)=0.75 v(b2)=0
.meas tran c1_max MAX v(c1) FROM=200u TO=500u
.tran 20n 500u
""")

_ring = """
Vdd vdd 0 DC 3.3
{inv}
.ic v(n1)=0 v(n2)=3.3 v(n3)=0 v(n4)=3.3 v(n5)=0
.tran {tstep} {tstop}
"""


def ring(models, nmod, pmod, wn, wp, lch, cl, tstep, tstop):
    inv = ""
    for i in range(1, 6):
        a, b = f"n{i}", f"n{i % 5 + 1}"
        inv += f"Mp{i} {b} {a} vdd vdd {pmod} W={wp} L={lch}\nMn{i} {b} {a} 0 0 {nmod} W={wn} L={lch}\nCl{i} {b} 0 {cl}\n"
    return models + _ring.format(inv=inv, tstep=tstep, tstop=tstop)


deck("tr_ring5_l1", "5-stage CMOS ring oscillator with level-1 MOSFETs (3.3 V, 10 fF loads, period ~1 ns), .ic breaks the symmetry: fast regenerative edges, many cycles", "v(n1) v(n3)",
     ring(MODELS["n1"] + "\n" + MODELS["p1"] + "\n", "nm1", "pm1", "10u", "25u", "1u", "50f", "5p", "30n") + ".meas tran vmax MAX v(n1) FROM=10n TO=30n\n")
deck("tr_ring5_bsim3_l8", "5-stage CMOS ring oscillator with BSIM3v3 (LEVEL=8) models (0.35u class, 3.3 V): heavy device evaluation, capacitance model, charge conservation", "v(n1) v(n3)",
     ring(MODELS["n3"].format(lvl=8) + "\n" + MODELS["p3"].format(lvl=8) + "\n", "nm3", "pm3", "2u", "5u", "0.35u", "20f", "5p", "20n") + ".meas tran vmax MAX v(n1) FROM=8n TO=20n\n")
deck("tr_ring5_bsim3_l49", "same ring oscillator with the BSIM3v3 LEVEL=49 (instance-parameter) alias: identical physics, different code path in the model card parser", "v(n1) v(n3)",
     ring(MODELS["n3"].format(lvl=49) + "\n" + MODELS["p3"].format(lvl=49) + "\n", "nm3", "pm3", "2u", "5u", "0.35u", "20f", "5p", "20n") + ".meas tran vmax MAX v(n1) FROM=8n TO=20n\n")

deck("tr_lc_negres", "LC oscillator start-up: parallel LC (1.6 MHz) with a saturating negative conductance (B-source, tanh): amplitude growth from 1 mV to a limit cycle", "v(n)", """
L1 n 0 10u
C1 n 0 1n
R1 n 0 5k
Bneg n 0 I = -0.4m*tanh(v(n)/0.3)*2.5
.ic v(n)=1m
.meas tran vpk MAX v(n) FROM=40u TO=50u
.tran 5n 50u uic
""")

deck("tr_colpitts", "BJT Colpitts oscillator (~4 MHz), supply ramp as the start-up kick, amplitude limited by the transistor: nonlinear start-up", "v(c) v(e)", MODELS["npn"] + """
Vcc vcc 0 PWL(0 0 2u 9 1)
Rb1 vcc b 47k
Rb2 b 0 15k
Cb b 0 100n
Q1 c b e qn
Re e 0 1k
Lc vcc c 10u
Rlc vcc c 20k
C1 c e 470p
C2 e 0 1n
.tran 5n 40u
""")

deck("tr_relax_sgn", "opamp-less relaxation oscillator with an ideal comparator (B-source sgn, positive feedback divider, RC timing): discontinuous source switching at threshold crossings", "v(c) v(out)", """
Bcmp out 0 V = sgn(v(c) - 0.5*v(out) )*-1
R1 out c 10k
C1 c 0 10n
.ic v(c)=0.1
.meas tran vmax MAX v(c) FROM=0.5m TO=1m
.tran 100n 1m
""")

# ---------------------------------------------------------------- B-sources
deck("tr_bsrc_tanh", "B-source with a huge-gain tanh (gain 1e6): a 1 kHz sine becomes a square wave with ~1 ns transitions that the step control has to find, then RC-filtered", "v(sq) v(f)", """
V1 in 0 SIN(0 1 1k)
Bsq sq 0 V = tanh(1e6*v(in))
R1 sq f 1k
C1 f 0 100n
.meas tran vf FIND v(f) AT=2.75m
.tran 1u 3m
""")

deck("tr_bsrc_abs", "B-source abs(): ideal full-wave rectifier with a kink at every zero crossing, RC smoothing (non-differentiable, not discontinuous)", "v(f) v(r)", """
V1 in 0 SIN(0 2 1k)
Br r 0 V = abs(v(in))
R1 r f 1k
C1 f 0 1u
RL f 0 5k
.meas tran vavg AVG v(f) FROM=3m TO=5m
.tran 2u 5m
""")

deck("tr_bsrc_sgn", "B-source sgn(): ideal hard limiter (discontinuous at each zero crossing of a 1 kHz sine), RC loaded", "v(f) v(s)", """
V1 in 0 SIN(0 1 1k)
Bs s 0 V = sgn(v(in))
R1 s f 1k
C1 f 0 100n
.meas tran vf FIND v(f) AT=2.6m
.tran 1u 3m
""")

deck("tr_bsrc_expdiode", "B-source diode I = 1e-14*(exp(v/0.0258)-1) without any limiting, 10 V sine through 100 ohm: exp() overflow protection in Newton", "v(out) v(in)", """
V1 in 0 SIN(0 10 1k)
R1 in out 100
Bd out 0 I = 1e-14*(exp(v(out)/0.0258)-1)
Cp out 0 10p
.meas tran vpk MAX v(out) FROM=0 TO=3m
.tran 1u 3m
""")

# ---------------------------------------------------------------- sharp edges / stiff
deck("tr_pwl_edges", "PWL source with 1 ps edges inside a 1 ms run feeding two RC branches (tau 1 ns and 100 us): the step control must resolve picosecond edges and then stretch to microseconds", "v(fast) v(slow)", """
V1 in 0 PWL(0 0 100u 0 100.001n 1 300u 1 300.001n 0 500u 0 500.001n 1 520u 1 520.001n 0 700u 0 700.001n 0.5 1m 0.5)
R1 in fast 100
C1 fast 0 10p
R2 in slow 10k
C2 slow 0 10n
.meas tran vslow FIND v(slow) AT=0.95m
.tran 1u 1m
""")

deck("tr_stiff_mix", "stiff network, time constants 1 ps and 1 ms together (R=1 ohm,C=1 pF branch; R=1k,C=1 uF branch) behind a 1 ps edge and a diode", "v(a) v(b) v(c)", """
.model dd D(IS=1e-14 N=1 RS=1 CJO=1p)
V1 in 0 PULSE(0 2 10n 1p 1p 3m 6m)
R1 in a 1
C1 a 0 1p
R2 a b 100
D1 b c dd
R3 c d 1k
C3 d 0 1u
R4 c 0 100k
.meas tran vd FIND v(d) AT=2.5m
.tran 100n 4m
""")

deck("tr_tmax_tiny", "very small tmax relative to the span: tmax = 2 ns on a 400 us run with an RC (200k forced points) -- cost scales with the step ceiling, not the dynamics", "v(out)", """
V1 in 0 PULSE(0 1 10u 100n 100n 100u 300u)
R1 in out 1k
C1 out 0 10n
.tran 2n 400u 0 2n
""")

# ---------------------------------------------------------------- latches / thyristor / switched networks / tline
deck("tr_latch_sr_cmos", "CMOS NOR SR latch (level-1 MOS, 3.3 V) toggled by set/reset pulses every 2 us: bistable regeneration, feedback, floating internal nodes of the series PMOS", "v(q) v(qb)", MODELS["n1"] + "\n" + MODELS["p1"] + """
Vdd vdd 0 DC 3.3
Vs s 0 PULSE(0 3.3 1u 1n 1n 500n 4u)
Vr r 0 PULSE(0 3.3 3u 1n 1n 500n 4u)
* NOR(a,b) = x : series PMOS to vdd, parallel NMOS
Mp1 vdd s m1 vdd pm1 W=20u L=1u
Mp2 m1 qb q vdd pm1 W=20u L=1u
Mn1 q s 0 0 nm1 W=8u L=1u
Mn2 q qb 0 0 nm1 W=8u L=1u
Mp3 vdd r m2 vdd pm1 W=20u L=1u
Mp4 m2 q qb vdd pm1 W=20u L=1u
Mn3 qb r 0 0 nm1 W=8u L=1u
Mn4 qb q 0 0 nm1 W=8u L=1u
Cq q 0 20f
Cqb qb 0 20f
.ic v(q)=0 v(qb)=3.3
.tran 1n 8u
""")

deck("tr_scr_bjt", "thyristor (SCR) made of a PNP and an NPN BJT: a gate pulse fires it, the anode supply dips below the holding current to turn it off, a second gate pulse refires it: regenerative latch-up and drop-out", "v(a) v(k) v(g2)", MODELS["npn"] + "\n" + MODELS["pnp"] + """
Va va 0 PWL(0 20 150u 20 200u 0.2 300u 0.2 350u 20 600u 20)
Ra va a 100
Qp g2 g1 a qp
Qn g1 g2 k qn
Rk k 0 50
Rb1 a g1 100k
Rb2 g2 k 10k
Vgate gt 0 PULSE(0 5 20u 100n 100n 5u 400u)
Rg gt g2 1k
.tran 100n 600u
""")

deck("tr_chargeshare_sw", "switched-capacitor charge sharing: 2-phase non-overlapping clocks, SW models with RON = 1 uOhm (time constant ~1 ps) moving charge between 1 uF, 100 nF and 1 uF: zero-resistance-like switches, charge conservation", "v(c2) v(c3)", """
.model swz SW(RON=1u ROFF=1e12 VT=2.5 VH=0.1)
V1 in 0 DC 5
Vp1 p1 0 PULSE(0 5 0 5n 5n 800n 2u)
Vp2 p2 0 PULSE(0 5 1u 5n 5n 800n 2u)
S1 in c1 p1 0 swz
S2 c1 c2 p2 0 swz
S3 c2 c3 p1 0 swz
C1 c1 0 100n IC=0
C2 c2 0 1u IC=0
C3 c3 0 1u IC=0
.meas tran v3 FIND v(c3) AT=190u
.tran 10n 200u uic
""")

deck("tr_charge_pump", "4-stage Dickson charge pump (diode-capacitor ladder, 1 MHz complementary 5 V clocks, 1 nF stage caps, 100 kOhm load): diode conduction bursts every half cycle, charge pumping to ~15 V", "v(out) v(n2)", MODELS["d"] + """
V1 clk 0 PULSE(0 5 0 10n 10n 490n 1u)
V2 clkb 0 PULSE(5 0 0 10n 10n 490n 1u)
Vin in 0 DC 5
D0 in n1 dsch
C1 n1 clk 1n
D1 n1 n2 dsch
C2 n2 clkb 1n
D2 n2 n3 dsch
C3 n3 clk 1n
D3 n3 out dsch
Co out 0 1n
RL out 0 100k
.meas tran vout_end FIND v(out) AT=95u
.tran 10n 100u
""")

deck("tr_tline_mismatch", "lossless transmission line (T element, Z0 50, TD 2 ns) from a 50 ohm source into a 200 ohm || 2 pF load (reflection coefficient 0.6): delayed-history interpolation, breakpoints at echoes", "v(a) v(b)", """
V1 in 0 PULSE(0 1 1n 100p 100p 20n 60n)
Rs in a 50
T1 a 0 b 0 Z0=50 TD=2n
RL b 0 200
Cl b 0 2p
.meas tran vb FIND v(b) AT=30n
.tran 20p 60n
""")

# ---------------------------------------------------------------- decks aimed at stock failures
deck("tr_chatter_bsrc", "sliding-mode chatter: ideal comparator B-source (v(c)>0.5 ? 0 : 1) closing a loop through an RC; the solution sits ON the discontinuity", "v(c) v(out)", """
Bcmp out 0 V = v(c)>0.5 ? 0 : 1
R1 out c 1k
C1 c 0 1n
.tran 10n 100u
""")

deck("tr_chatter_bsrc_uic", "same sliding-mode comparator loop, started with uic from c=0 (no operating point): chatter at t > 0 instead of in the DC solution", "v(c) v(out)", """
Bcmp out 0 V = v(c)>0.5 ? 0 : 1
R1 out c 1k
C1 c 0 1n
.ic v(c)=0
.tran 10n 100u uic
""")

deck("tr_ind_interrupt", "inductor current interrupted by a switch without a clamp: 1 mH carrying 10 A, SW opens (RON 1 mOhm, ROFF 1e12) into 1 fF: dI/dt energy has nowhere to go (kV-scale spike)", "v(b) l1#branch", """
.model swi SW(RON=1m ROFF=1e12 VT=2.5 VH=0.1)
V1 in 0 DC 10
Vg g 0 PULSE(5 0 10u 10n 10n 1 2)
R1 in a 1m
L1 a b 1m IC=10
S1 b 0 g 0 swi
Csw b 0 1f
.tran 10n 30u uic
""")

deck("tr_float_cap_switch", "floating nodes after switching: a B-source current switch (exactly zero conductance when off) leaves node y connected to the rest only through capacitors; no DC path at the operating point either", "v(x) v(y)", """
V1 in 0 PULSE(0 5 10u 1u 1u 30u 100u)
C1 in x 1n
C2 x y 1n
Vc ctl 0 PULSE(0 5 50u 1u 1u 20u 200u)
Bsw y 0 I = (v(ctl) > 2.5) ? v(y)*1e-3 : 0
.tran 100n 150u
""")

deck("tr_sw_ron_tiny", "ideal-ish switch RON = 1e-12 ohm, ROFF = 1e12 (24 decades of conductance) shorting a charged 1 nF into an LC every 10 us", "v(a) v(m)", """
.model swt SW(RON=1e-12 ROFF=1e12 VT=2.5 VH=0)
V1 in 0 DC 5
Vg g 0 PULSE(0 5 1u 1n 1n 4u 10u)
R1 in a 1k
C1 a 0 1n
S1 a m g 0 swt
L1 m n 1u IC=0
C2 n 0 1n IC=0
.tran 10n 100u
""")

deck("tr_sw_nohyst_osc", "RC oscillator whose SW discharges the capacitor with VH = 0 (no hysteresis): the switch controls itself from the node it switches", "v(c)", """
.model swn SW(RON=100 ROFF=1e9 VT=3 VH=0)
V1 in 0 DC 5
R1 in c 10k
C1 c 0 1n IC=0
S1 c 0 c 0 swn
.tran 10n 200u uic
""")

deck("tr_ring3_hardtanh", "3-inverter ring of B-source inverters V = -tanh(1e9*v) (gain 1e9) with RC delays: oscillates with ~1e9 V/V slopes at every crossing", "v(n1) v(n2)", """
B1 m2 0 V = -tanh(1e9*v(n1))
B2 m3 0 V = -tanh(1e9*v(n2))
B3 m1 0 V = -tanh(1e9*v(n3))
R1 m1 n1 1k
R2 m2 n2 1k
R3 m3 n3 1k
C1 n1 0 100p
C2 n2 0 100p
C3 n3 0 100p
.ic v(n1)=0.1 v(n2)=0 v(n3)=0
.tran 1n 5u
""")


deck("tr_schmitt_bjt", "emitter-coupled BJT Schmitt trigger driven by a very slow triangle (6 V in 10 ms): hysteresis, regenerative snap at each threshold with a nearly static input", "v(out) v(in)", MODELS["npn"] + """
Vcc vcc 0 DC 10
Vin in 0 PWL(0 0 10m 6 20m 0)
Rc1 vcc c1 2k
R1 c1 b2 5k
R2 b2 0 4k
Rc2 vcc out 1k
Re e 0 500
Q1 c1 in e qn
Q2 out b2 e qn
.meas tran vhi MAX v(out) FROM=0 TO=20m
.tran 10u 20m
""")

deck("tr_srlatch_bsrc", "SR latch of two delay-less B-source NOR gates (V = cond ? 0 : 1) cross-coupled directly: algebraic loop through discontinuities, even the operating point is ambiguous", "v(q) v(qb)", """
Vs s 0 PULSE(0 1 1u 1n 1n 1u 4u)
Vr r 0 PULSE(0 1 3u 1n 1n 1u 4u)
Bq q 0 V = (v(r)+v(qb)) > 0.5 ? 0 : 1
Bqb qb 0 V = (v(s)+v(q)) > 0.5 ? 0 : 1
.tran 10n 8u
""")

deck("tr_sw_selfctl", "SW whose control node is the node it switches (VT 1 V, VH 0): a 2 V pulse through 1 kOhm charges a 1 pF node until the switch clamps it: self-referential switch state, chatter at the threshold", "v(a)", """
.model swn SW(RON=1 ROFF=1e9 VT=1 VH=0)
V1 in 0 PULSE(0 2 1u 1n 1n 20u 40u)
R1 in a 1k
S1 a 0 a 0 swn
C1 a 0 1p
.tran 1n 40u
""")


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    for name, (stress, probe, body) in sorted(D.items()):
        head = [f"* {name} -- {stress}"]
        head.append(f"* stock: {STOCK.get(name, 'completes at default settings')}")
        head.append(f"* probe: {probe}")
        with open(os.path.join(outdir, name + ".cir"), "w") as f:
            f.write("\n".join(head) + "\n" + body + ".end\n")
    print(f"{len(D)} decks written to {outdir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else OUT)
