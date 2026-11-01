#!/usr/bin/env python3
"""bench/campaign_pvt_bench_adapter.py -- adapts the EXISTING `c_pvt`
bench artifact (`bench/gen_pvt.py`'s 50-point grid, card #467 T9.1) to
`tools/campaign.py`'s module interface, for card #491's "measure c_pvt
(bench)" task. Does not regenerate or modify `bench/decks/pvt_points/*.cir`
or `bench/decks/c_pvt_loop.cir` -- reads the already-checked-in files
verbatim, so the campaign.py-measured numbers are directly comparable to
`bench/README.md`'s own recorded baseline (50 processes serial = 5.51s,
same 50 points in one session = 3.9s) for the SAME files, not a
regenerated near-copy.
"""
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
POINTS_DIR = HERE / "decks" / "pvt_points"
LOOP_DECK = HERE / "decks" / "c_pvt_loop.cir"

NGSPICE_DEFAULT = None
TIMEOUT = 30

N = 50
POINTS = [{"idx": i} for i in range(N)]
MEASURES = ["vout_avg_final", "vout_ripple", "wave_mid", "wave_end"]


def build_point_deck(point, idx):
    path = POINTS_DIR / f"pt_{idx:03d}.cir"
    return path.read_text()
