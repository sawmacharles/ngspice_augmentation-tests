#!/usr/bin/env python3
"""summarise thread_sweep.py JSON: median/min/max wall per deck x cfg, % vs a chosen cfg"""
import json, sys, statistics, collections
d = json.load(open(sys.argv[1])); ref = sys.argv[2] if len(sys.argv) > 2 else None
g = collections.OrderedDict()
for r in d: g.setdefault(r["deck"], collections.OrderedDict()).setdefault(r["cfg"], []).append(r)
for deck, cs in g.items():
    base = statistics.median(x["wall"] for x in cs[ref]) if ref in cs else None
    print(deck)
    for c, rs in cs.items():
        w = [x["wall"] for x in rs]; m = statistics.median(w)
        pc = "" if base is None else " %+6.1f%%" % (100 * (m / base - 1))
        print("  %-10s med %8.3f min %8.3f max %8.3f n=%d cpu %7.2f load %4.1f-%4.1f dig %s rc %s%s" % (
            c, m, min(w), max(w), len(w), statistics.median(x["cpu"] for x in rs),
            min(x["load1"] for x in rs), max(x["load1"] for x in rs),
            len(set(x["digest"] for x in rs)), sorted(set(x["rc"] for x in rs)), pc))
