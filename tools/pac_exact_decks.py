# tools/pac_exact_decks.py (cards #713, #728, 2026-10-19): deck texts copied from tests/test_pac_charge.py (BOOST_579, BOOSTC_580, BB_583) for tools/pac_exact.py
BOOST_579 = ("* boost, low-side ideal switch, output diode\n.param fs=500k tclk={1/fs}\nVin vin 0 DC 5\nRdcr vin vl 10m\nL1 vl sw 10u IC=0\n"
             "S1 sw 0 ctrlnode ramp SWMOD\n.model SWMOD SW(vt=0 vh=0.005 ron=10m roff=1Meg)\nD1 sw vout DMOD\n"
             ".model DMOD D(IS=1e-9 RS=0.05 N=1.05 CJO=%s VJ=0.7 M=0.5 TT=%s)\nCout vout cx 22u IC=0\nResr cx 0 10m\nRload vout 0 10\n"
             "Vramp ramp 0 PULSE(0 2 0 {0.97*tclk} {0.03*tclk} 1n {tclk})\nVctrl ctrlnode 0 DC 0.8 AC 1\n")
BOOSTC_580 = BOOST_579.replace("Cout vout cx 22u IC=0", "Csw sw 0 %s\nCout vout cx 22u IC=0")
BB_583 = ("* inverting buck-boost\n.param fs=500k tclk={1/fs}\nVin vin 0 DC 12\nS1 vin sw ctrlnode ramp SWMOD\n.model SWMOD SW(vt=0 vh=0.005 ron=10m roff=1Meg)\n"
          "L1 sw lg 10u IC=0\nRdcr lg 0 10m\nD1 vout sw DMOD\n.model DMOD D(IS=1e-9 RS=0.05 N=1.05 CJO=%s VJ=0.7 M=0.5 TT=%s)\n%s"
          "Cout vout cx 22u IC=0\nResr cx 0 10m\nRload vout 0 10\nVramp ramp 0 PULSE(0 2 0 {0.97*tclk} {0.03*tclk} 1n {tclk})\nVctrl ctrlnode 0 DC 0.8 AC 1\n")
# name: (topology, cjo, tt, csw, deck text, pinned ref {f: (dB, deg)} from the test)
DECKS = {
    "bst_cl_csw": ("boost", 0.0, 0.0, 1e-9, BOOSTC_580 % ("0", "0", "1n"), {20000: (-1.6851, 169.0429), 50000: (-16.532, 146.3875), 100000: (-25.1428, 130.1333)}),
    "bst_q_csw": ("boost", 100e-12, 2e-9, 1e-9, BOOSTC_580 % ("100p", "2n", "1n"), {20000: (-1.6989, 169.0004), 50000: (-16.5319, 146.3029), 100000: (-25.1456, 130.1234)}),
    "bst_cl_csw100p": ("boost", 0.0, 0.0, 100e-12, BOOSTC_580 % ("0", "0", "100p"), {20000: (-1.6795, 168.8301), 50000: (-16.5465, 146.5302), 100000: (-25.1741, 130.3579)}),
    "bst_q_csw3n": ("boost", 100e-12, 2e-9, 3e-9, BOOSTC_580 % ("100p", "2n", "3n"), {20000: (-1.7107, 169.4522), 50000: (-16.4893, 145.9996), 100000: (-25.0794, 129.4944)}),
    "bb_cl_csw": ("bb", 0.0, 0.0, 1e-9, BB_583 % ("0", "0", "Csw sw 0 1n\n"), {20000: (5.625, 1.2047), 50000: (-10.6306, -11.9392), 100000: (-21.5722, -24.3383)}),
    "bb_q_csw": ("bb", 100e-12, 2e-9, 1e-9, BB_583 % ("100p", "2n", "Csw sw 0 1n\n"), {20000: (5.6084, 1.2172), 50000: (-10.644, -12.002), 100000: (-21.5716, -24.4017)}),
}
