#!/usr/bin/env python3
"""16_yield_run.py -- handbook example driver for tools/yieldest.py
(docs/yield.md). Linear resistor stack (16_yield_deck.cir), 4-sigma target,
run against the pinned system ngspice (this tool needs no private-fork
feature). Ground truth (Phi(-4)) is printed alongside for comparison, the
same way tests/test_yield.py checks it -- computed independently here, not
read from tools/yieldest.py's own output.
"""
import math
import os
import sys
from statistics import NormalDist

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import yieldest  # noqa: E402

NGSPICE = os.environ.get("NGSPICE", "ngspice")
TEMPLATE = open(os.path.join(HERE, "16_yield_deck.cir")).read()

R_NOM = 10000.0
R_SIGMA = R_NOM * 0.03
I_SRC = 1e-3
VTOP_NOM = I_SRC * 2 * R_NOM
SIGMA_VTOP = I_SRC * math.sqrt(2) * R_SIGMA
K_SIGMA = 4.0
THRESHOLD = VTOP_NOM - K_SIGMA * SIGMA_VTOP

spec = {
    "template": TEMPLATE, "seed": 1, "ngspice": NGSPICE,
    "params": [
        {"name": "R1", "nominal": R_NOM, "sigma": R_SIGMA},
        {"name": "R2", "nominal": R_NOM, "sigma": R_SIGMA},
    ],
    "measure": "vtop", "fail": {"op": "lt", "threshold": THRESHOLD},
    "n_is": 300, "max_workers": 2,
}

outdir = os.path.expanduser("~/Work/data/ngspice_augmentation/handbook/16run")
os.makedirs(outdir, exist_ok=True)
result = yieldest.estimate_yield(spec, outdir)
nd = NormalDist()
pf_true = nd.cdf(-K_SIGMA)

print(f"n_sims={result.get('n_sims')} beta_mpp={result['mpp_search']['beta']:.4f}")
print(f"importance-sampling pf_estimate={result['pf_estimate']:.4e} "
      f"ci95={result.get('pf_ci95')}")
print(f"true Phi(-{K_SIGMA})={pf_true:.4e} "
      f"ratio={result['pf_estimate']/pf_true:.3f}")
print(f"scaled-sigma-sampling pf={result.get('scaled_sigma_sampling', {}).get('pf')}")
print(f"limitation_warning={result.get('limitation_warning')}")
