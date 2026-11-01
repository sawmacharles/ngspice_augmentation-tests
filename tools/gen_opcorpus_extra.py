# second batch of corpus decks (exec'd by gen_opcorpus.py; shares its namespace).
# Added after the first baseline showed that Newton alone solved nearly every first-batch deck.

deck("bsrc_exp_50v", "two B-source exponential diodes (Is=1e-15, VT=26 mV) in series, 50 V through 0.1 ohm",
     "no junction limiting exists for B-source exp(): the first step lands at tens of volts and exp(V/VT) overflows; the root has ~0.8 V across each device",
     """V1 in 0 50
R1 in a 0.1
Bd1 a b I=1e-15*(exp((V(a)-V(b))/0.026)-1)
Bd2 b 0 I=1e-15*(exp(V(b)/0.026)-1)
""")

deck("bsrc_exp_chain6", "six exponential B-source diodes in series from 30 V through 0.5 ohm, plus a shunt exp device at each tap through 1 kohm",
     "chain of unlimited exponentials where each tap is only pinned by the next device; the Newton path overshoots at every node simultaneously",
     "V1 in 0 30\nR1 in n0 0.5\n" + "".join(
         f"Bd{i} n{i} {'n'+str(i+1) if i < 5 else '0'} I=1e-15*(exp((V(n{i})-{'V(n'+str(i+1)+')' if i < 5 else '0'})/0.026)-1)\n"
         f"Rt{i} n{i} t{i} 1k\nBt{i} t{i} 0 I=1e-15*(exp(V(t{i})/0.026)-1)\n" for i in range(6)))

s = "Vdd vdd 0 1.8\nVin s0 0 0.78\n"
for i in range(60):
    s += inv(str(i + 1), f"s{i}", f"s{i+1}", "vdd", "nb3", "pb3", "1u", "2.5u", "0.18u")
s += "Rld s60 0 1Meg\n" + NB3 + "\n" + PB3 + "\n"
deck("inv_chain60_bsim3", "60-stage CMOS inverter chain (BSIM3), input parked at the first stage's switching region",
     "60 stages of gain ~10-20 in series (1e60 overall); every stage is in a different region on the path from the all-zero start",
     s)

s = "Vdd vdd 0 1.8\nVin s0 0 0.82\n"
for i in range(100):
    s += inv(str(i + 1), f"s{i}", f"s{i+1}", "vdd", "n1", "p1", "1u", "2.5u", "0.18u")
s += "Rld s100 0 1Meg\n" + N1 + "\n" + P1 + "\n"
deck("inv_chain100_l1", "100-stage CMOS inverter chain (level 1, lambda 0.05), input near the trip point",
     "level-1 inverters have a steep transfer and no leakage; 100 stages give ~1e100 gain and 100 coupled switching points",
     s)

rng2 = random.Random(6370)
R_, C_ = 4, 4
s = "Vdd vdd 0 1.8\nVbl bl 0 1.8\nVblb blb 0 1.8\nVwl wl 0 0\n"
for r_ in range(R_):
    for c_ in range(C_):
        k = f"{r_}_{c_}"
        w = lambda lo, hi: round(rng2.uniform(lo, hi), 2)
        s += (f"Mn1_{k} qb{k} q{k} 0 0 nb4 w={w(0.4,0.7)}u l=0.18u\nMp1_{k} qb{k} q{k} vdd vdd pb4 w={w(0.25,0.45)}u l=0.18u\n"
              f"Mn2_{k} q{k} qb{k} 0 0 nb4 w={w(0.4,0.7)}u l=0.18u\nMp2_{k} q{k} qb{k} vdd vdd pb4 w={w(0.25,0.45)}u l=0.18u\n"
              f"Ma1_{k} bl wl q{k} 0 nb4 w=0.3u l=0.18u\nMa2_{k} blb wl qb{k} 0 nb4 w=0.3u l=0.18u\n")
s += NB4 + "\n" + PB4 + "\n"
deck("sram_array4x4", "4x4 SRAM array (96 BSIM4 transistors), random per-cell device widths (seed 6370), word lines off",
     "16 independent bistable cells, each with its own (asymmetric) metastable point; the global Newton path has to resolve 16 symmetry-breaking decisions at once",
     s)

deck("bgr_cmos_opamp", "Kuijk/Banba-style bandgap: PMOS mirror (3 branches), ideal tanh amplifier equalising the two branch voltages, 1x and 8x PNP (substrate-style) with R, output branch R3+PNP",
     "zero-current root (mirror off) plus the intended root; the amplifier drives the PMOS gate through a saturating tanh so Newton from 0 sits on a flat stretch",
     f"""Vdd vdd 0 3.3
Bop pg 0 V=1.65*(1+tanh(500*(V(b)-V(a))))
Mp1 a pg vdd vdd p3 w=40u l=1u
Mp2 b pg vdd vdd p3 w=40u l=1u
Mp3 ref pg vdd vdd p3 w=40u l=1u
Q1 0 0 a qs
Rb b e2 6k
Q2 0 0 e2 qs area=8
R3 ref r3 2k
Q3 0 0 r3 qs
Rs ref 0 1e9
{P3}
.model qs pnp is=1e-16 bf=50 vaf=50 rb=100 re=2
""")

s = "Vdd vdd 0 1.8\nVin in0 0 0.9\nIb vdd bias 20u\nMb bias bias 0 0 nb3 w=4u l=0.5u\n"
prev = "in0"
for k in range(5):
    s += (f"Mt{k} tail{k} bias 0 0 nb3 w=16u l=0.5u\nM1_{k} d1_{k} {prev} tail{k} 0 nb3 w=20u l=0.5u\n"
          f"M2_{k} d2_{k} o{k} tail{k} 0 nb3 w=20u l=0.5u\nM3_{k} d1_{k} d1_{k} vdd vdd pb3 w=20u l=0.5u\n"
          f"M4_{k} d2_{k} d1_{k} vdd vdd pb3 w=20u l=0.5u\nM6_{k} o{k} d2_{k} vdd vdd pb3 w=60u l=0.5u\n"
          f"M7_{k} o{k} bias 0 0 nb3 w=64u l=0.5u\nCc{k} d2_{k} o{k} 1p\n")
    prev = f"o{k}"
s += NB3 + "\n" + PB3 + "\n"
deck("opamp2_follow_chain5", "five two-stage Miller opamps (BSIM3) cascaded as unity-gain followers, input 0.9 V",
     "five nested high-gain (1e5) feedback loops, each with an output at a rail until its own loop closes; 55 transistors",
     s)

s = "Vcc vcc 0 15\nVin in 0 0.6\nRin in b0 10k\n"
for k in range(8):
    s += f"Rc{k} vcc c{k} {2000 + 300 * k}\nQ{k} c{k} b{k} e{k} qn\nRe{k} e{k} 0 {200 + 50 * k}\n"
    if k < 7:
        s += f"Rco{k} c{k} b{k+1} 1k\n"
s += "Rfb c7 b0 200k\n" + QN + "\n"
deck("bjt_cascade8_fb", "eight direct-coupled common-emitter BJT stages with a global DC feedback resistor from the last collector to the first base",
     "total gain 1e8+ in a global loop; every stage swings from cutoff to saturation along the Newton path, wild Vbe values at the start",
     s)

n = 40
deck("diode_stack_200v", f"{n} diodes in series across 200 V through 0.1 ohm",
     "currents of ~1500 A at about 1.2 V per junction; per-junction step limiting is taxed (5 V of drive per junction at the start)",
     "V1 n0 0 200\nR1 n0 n1 0.1\n" + "".join(f"D{i} n{i} {'n'+str(i+1) if i < n else '0'} dm\n" for i in range(1, n + 1)) + DM + "\n")

s = "Vcc vcc 0 5\n"
for k in range(6):
    s += (f"Rl{k}a vcc a{k} 5k\nRl{k}b vcc b{k} 5k\nQa{k} a{k} b{k} e{k} qn\nQb{k} b{k} a{k} e{k} qn\nIt{k} e{k} 0 DC 200u\n")
    nxt = (k + 1) % 6
    s += f"Rx{k} a{k} b{nxt} 30k\nRy{k} b{k} a{nxt} 47k\n"
s += QN + "\n"
deck("latch_ring_bjt6", "six cross-coupled BJT latches (200 uA tails) coupled in a ring by 30k/47k resistors",
     "2^6 latch states plus metastable mixtures; the coupling is asymmetric so the all-metastable point is not a root and Newton must break symmetry",
     s)

s = "Vdd vdd 0 3.3\nVg g 0 3.3\n"
prev = "vdd"
for k in range(12):
    nx = "0" if k == 11 else f"x{k}"
    s += f"M{k} {prev} g {nx} 0 n3 w=6u l=0.35u\n"
    prev = nx
s += N3 + "\n"
deck("mos_stack12", "twelve series NMOS (level 3, body at ground) from a 3.3 V rail to ground, all gates at 3.3 V",
     "NAND12-like pull-down stack: body effect grows up the stack, the intermediate nodes start at 0 with most devices deep in saturation and cut off",
     s)

rng3 = random.Random(6371)
NN = 14
wts = []
for i in range(NN):
    wts.append(([(j, round(rng3.uniform(-1.5, 1.5), 3)) for j in range(NN) if j != i and rng3.random() < 0.5],
                round(rng3.uniform(-0.2, 0.2), 3)))
for g in (25, 3):
    s = "Vb vb 0 0.1\n"
    for i, (tl, b0) in enumerate(wts):
        terms = "+".join(f"{w}*V(h{j})" for j, w in tl)
        s += f"Bh{i} h{i} 0 V=tanh({g}*({terms or '0'}+{b0}))\nRh{i} h{i} 0 1Meg\n"
    deck("tanh_hopfield14" + ("" if g == 25 else "_g3"),
         f"14-unit Hopfield-style network of tanh B-sources (gain {g}) with random +-1.5 couplings (seed 6371)",
         "dense recurrent tanh network: many fixed points, many saturated flat regions where the Newton derivative vanishes"
         + ("; gain 25 is expected to defeat every fallback (stress control)" if g == 25 else ""),
         s)

deck("bjt_exp_loop", "BJT stage whose base current is a B-source exp() of the input voltage (times a (1+V(b)) factor), 5 V rail",
     "unlimited exp() in the base drive on top of a limited BJT junction; the two exponentials multiply",
     f"""Vcc vcc 0 5
Vin in 0 0.5
Bb vcc b I=1e-12*(exp(V(in)/0.0259)-1)*(1+V(b))
Rb b e1 10
Q1 c1 b e1 qn
Rc vcc c1 100
Re e1 0 10
Ce c1 0 1p
{QN}
""")

s = "Vdd vdd 0 1.8\nVin s0 0 0.79\n"
for i in range(100):
    s += inv(str(i + 1), f"s{i}", f"s{i+1}", "vdd", "nb3", "pb3", "1u", "2.5u", "0.18u")
s += "Rld s100 0 1Meg\n" + NB3 + "\n" + PB3 + "\n"
deck("inv_chain100_bsim3", "100-stage CMOS inverter chain (BSIM3), input at the first stage's switching region",
     "BSIM3 with subthreshold tails: 100 stages (~1e100 gain) whose Newton path crosses every stage's transition region in turn",
     s)

s = "Vdd vdd 0 1.8\n"
for i in range(20):
    s += inv(str(i + 1), f"s{i}", f"s{(i + 1) % 20}", "vdd", "nb3", "pb3", "1u", f"{2.3 + 0.02 * i:.2f}u", "0.18u")
s += NB3 + "\n" + PB3 + "\n"
deck("ring_even20_bsim3", "20-stage (even) BSIM3 inverter ring with slightly graded PMOS widths (breaks the symmetry of the metastable point)",
     "two stable alternating states plus a metastable one; with the asymmetry the all-equal start is not a root and 20 stages of gain must resolve consistently",
     s)

M = 6
s = "Vs top 0 20\nRs top m0_0 1\n"
for i in range(M):
    for j in range(M):
        if j + 1 < M:
            s += f"Rh{i}_{j} m{i}_{j} m{i}_{j+1} 10\n"
        if i + 1 < M:
            s += f"Rv{i}_{j} m{i}_{j} m{i+1}_{j} 10\n"
        s += f"Bd{i}_{j} m{i}_{j} 0 I=1e-15*(exp(V(m{i}_{j})/0.026)-1)\n"
deck("bsrc_exp_mesh6", "6x6 resistive mesh (10 ohm links) with a B-source exponential diode at every node, fed by 20 V through 1 ohm",
     "36 unlimited exponentials all coupled through a resistive mesh: the first step puts ~20 V on exp(V/26 mV) at the feed corner and propagates it",
     s)

s = "V1 vcc 0 40\nVin in 0 12\nR1 in b1 50\n"
for k in range(10):
    s += f"Q{k} vcc b{k+1 if k else 1} e{k} qn\n"
s = "V1 vcc 0 40\nVin in 0 7.5\nRin in b0 50\n"
for k in range(10):
    s += f"Q{k} vcc b{k} b{k+1} qn\n"
s += "Re b10 0 100\n" + QN + "\n"
deck("bjt_darlington10", "ten emitter-follower BJTs in a chain (10-deep Darlington), 40 V supply, base driven to 7.5 V through 50 ohm",
     "ten series Vbe drops each needing ~0.75 V from an initially zero chain; collector currents span 1e-15..1 A and the base resistances run in high injection",
     s)

_t = decks["folded_cascode_ol"]
_t = _t.replace("Vn inn 0 0.9\n", "").replace("M2 f2 inn tail", "M2 f2 out tail")
_t = _t.replace("* op_folded_cascode_ol: folded-cascode OTA (BSIM4) open loop, inputs at mid supply", "* op_folded_cascode_follow: the folded-cascode OTA (BSIM4) as a unity-gain follower, input 0.9 V")
_t = _t.replace("* Hard for DC OP because: four stacked", "* Hard for DC OP because: DC feedback closes through gain ~1e3-1e4 of four stacked cascode branches; four stacked")
decks["folded_cascode_follow"] = _t

_t = decks["casc_selfbias"].replace("op_casc_selfbias:", "op_casc_selfbias_start:").replace("self-biased cascoded beta-multiplier", "self-biased cascoded beta-multiplier with a 50 Mohm startup bleed")
_t = _t.replace(".model n3", "Rst vdd ng 50Meg\n.model n3", 1)
decks["casc_selfbias_start"] = _t
