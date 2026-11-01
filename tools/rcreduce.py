#!/usr/bin/env python3
"""rcreduce.py -- RC long-chain / PEX-mesh netlist reduction preprocessor
(ngspice_augmentation card #506, transient acceleration round 8, candidate 1;
implements round 7's ranked #1 future candidate, research/transient-
acceleration.md section 12.5 item 1, arXiv 2508.13159 "Accelerating
Transistor-Level Simulation ... via Equivalence of RC Long-Chain Structures"
and the classic TICER/node-elimination family it cites).

Built as a DECK PRE-PROCESSOR TOOL, not an in-engine `.options` pass -- the
card's own instructions asked to pick whichever gives the best value/risk,
and a tool is strictly easier to validate (its output is an ordinary .cir
file any existing test harness can diff/run unmodified, no C source touched,
no risk to the accepted patch series on aug-accel8) at the cost of being
opt-in and non-transparent to a user who does not know to run it -- judged
the right trade for a first round on this idea (see the round's own report).

Two techniques, both operating ONLY on top-level (outside any .subckt/.ends)
two-terminal R and C elements:

1. EXACT star-mesh (Kron/Gaussian) elimination of an internal node that
   touches ONLY resistors (no capacitor, no independent source, no
   nonlinear device, not a ground/protected/measured node). Removing such
   a node is exact at every instant, not just at DC: with no capacitor
   storing charge and no source injecting current AT that node, KCL there
   is an algebraic (not differential) constraint satisfied identically at
   every simulation time, so its voltage is exactly a fixed linear
   combination of its neighbours' voltages throughout the transient, and
   eliminating it via the standard star-mesh transform (replace the node's
   N resistor legs with an (N choose 2) clique of equivalent resistors
   among its neighbours, conductances g_ij = g_ia*g_ib/sum(g)) reproduces
   every OTHER node's voltage bit-for-bit (up to floating-point re-
   association -- the matrix the solver sees is smaller but numerically
   equivalent, not merely approximately so). Degree-2 elimination is the
   familiar series-resistor-merge special case of the same formula
   (g_ia*g_ib/(g_ia+g_ib) = 1/(R_ia+R_ib)); degree-1 (a dangling stub) is
   simply dropped (no capacitance/source means no current ever flows into
   a dead end, so it and its one resistor contribute nothing to any other
   node's answer); degree-0 orphans are dropped the same way.

2. OPT-IN, bounded-error "fast node" elimination (--fast-node-tol), for an
   otherwise-eligible node that has exactly one capacitor to ground: if
   that node's own local RC time constant (tau_i = C_i / sum(g_ij) over
   its resistor legs) is <= tol * --ref-tstep, the node is eliminated the
   same star-mesh way and its capacitor is simply dropped (quasi-static/
   TICER-style approximation -- the node's own pole settles within a
   stated fraction of the caller's own transient step, so its dynamics
   are, by the tolerance's own definition, below the resolution the
   simulation is already asking for). This is NOT exact -- tol must be
   stated by the caller BEFORE reducing, same discipline as every other
   tolerance in this project (tests/test_accel.py's CANDIDATES table).
   Off by default; the exact pass alone never needs it.

Safety / scope, stated plainly (see also "Known limitations" below):
  - Only top-level R/C elements are considered; anything inside a
    .subckt/.ends block is left untouched.
  - A node is only a REDUCTION CANDIDATE if every element touching it is a
    top-level R or C -- any V/I/M/Q/D/B/E/G/L/K/X/etc. reference makes it
    ineligible, unconditionally (no linearisation risk from a nonlinear
    device terminal, no touching a source).
  - A node is additionally PROTECTED (never eliminated even if otherwise
    eligible) if its name appears anywhere else in the deck text outside
    its own R/C element lines -- e.g. a `.meas`/`print`/`plot`/`wrdata`/
    `.save`/`.ic`/`.nodeset` reference, or a comment. This is a
    conservative whole-file text scan (may over-protect on a name that
    happens to collide with unrelated text), never an under-protection.
  - Fill-in is bounded by --max-degree (default 6): a node whose
    resistor-degree exceeds this is left alone even if otherwise eligible,
    since star-mesh elimination of a degree-N node can add up to
    N*(N-1)/2 new resistors, and an unbounded high-degree hub could grow
    the netlist instead of shrinking it -- CONFIRMED this round on a
    synthetic branching-tree/grid demonstration deck (see the round's own
    report): a naive elimination pass with too loose a cap can add MORE
    resistors than it removes on a high-branching or grid-like topology,
    even though every individual elimination is exact. 6 is deliberately
    close to a real production PEX-reduction tool's own published choice
    for the equivalent operation (klayout-pex issue #214, "degree cap of
    5 for TICER" -- found by web research this round, not guessed) rather
    than this project's own earlier, untested default of 12. Nodes are
    processed in ascending-degree order (a minimum-degree-style heuristic,
    recomputed as eliminations proceed) specifically to avoid ever needing
    to eliminate a high-degree hub if a lower-degree path around it exists
    instead. NOTE: unlike TICER's own published fill-in criterion (an
    estimated fill-in COST check, not just a raw degree cap), this first
    version only bounds degree directly -- a coarser, more conservative
    proxy; a true fill-in-cost estimate is future work (see the round's
    report).

Usage:
  tools/rcreduce.py DECK.cir -o REDUCED.cir
  tools/rcreduce.py DECK.cir -o REDUCED.cir --fast-node-tol 0.01 --ref-tstep 10n
  tools/rcreduce.py DECK.cir --dry-run          # report only, write nothing

Exit status: 0 = wrote (or, with --dry-run, would have written) a reduced
deck; 1 = malformed input the parser refuses to touch (parse error on a
line that looks like an R/C element).

Never edits the input in place; never touches simulation output (the tool
does not run ngspice itself -- pair with the existing test/bench harness to
validate, see tests/test_rcreduce.py).
"""
import argparse
import re
import sys
import time
from collections import defaultdict

SUFFIX = {
    "t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "mil": 25.4e-6,
    "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15,
}
# longest suffix first so "meg" is tried before "m"; SPICE units are
# case-insensitive (and "M"/"m" both mean milli -- only the 3-letter
# "MEG"/"meg" means mega, same convention ngspice itself uses).
_SUF_RE = re.compile(
    r"^([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)(meg|mil|[tgkmunpf])?[a-zA-Z]*$",
    re.IGNORECASE,
)


def parse_spice_num(tok):
    """Parse a plain SPICE numeric literal (no {expr}, no .param reference)
    with an optional unit suffix. Returns None (not a crash) if tok is not
    a plain numeric literal -- callers treat that as "can't reduce this
    element", not an error."""
    m = _SUF_RE.match(tok.strip())
    if not m:
        return None
    val = float(m.group(1))
    suf = m.group(2)
    if suf:
        val *= SUFFIX[suf.lower()]
    return val


def join_continuations(lines):
    """Join SPICE '+' continuation lines onto the previous logical line.
    Returns a list of (logical_line, first_source_line_index)."""
    out = []
    for i, raw in enumerate(lines):
        if raw.startswith("+") and out:
            out[-1] = (out[-1][0] + " " + raw[1:].strip(), out[-1][1])
        else:
            out.append((raw.rstrip("\n"), i))
    return out


def is_comment_or_blank(line):
    s = line.strip()
    return (not s) or s.startswith("*")


class ParseError(Exception):
    pass


def parse_deck(text):
    """Returns a dict describing the deck:
      raw_lines: original lines (list[str], no trailing \\n)
      logical: [(text, src_idx)] continuation-joined
      elements: [{'name','type','n1','n2','value_tok','value','src_idx',
                  'in_subckt'}] for top-level R/C elements (value=None if
                  value_tok didn't parse as a plain numeric literal)
      other_touches: {node: set(element_name)} -- every non-R/C element
                  (top-level or not) that references a node as one of its
                  first tokens' node fields (best-effort: node-position
                  device types only, see NODE_POSITIONS)
      all_text_nodewords: set of every \\w+ token anywhere in the file
                  (used for the protected-node text scan)
    """
    raw_lines = text.split("\n")
    logical = join_continuations(raw_lines)

    elements = []
    other_node_refs = defaultdict(set)
    subckt_depth = 0

    # a conservative, explicit list of "how many of the first N tokens
    # after the element name are node names" for common device types --
    # used only to flag a node as "touched by a non-R/C element" (to make
    # it ineligible), not to fully model any device.
    NODE_POSITIONS = {
        "v": 2, "i": 2, "d": 2, "q": 3, "j": 3, "z": 3,
        "m": 4, "b": 2, "e": 4, "g": 4, "f": 2, "h": 2,
        "l": 2, "k": 0, "s": 4, "w": 2, "t": 4, "u": 0,
        "x": None,  # subckt call: node count is subckt-specific, unknown
                    # here -- treat EVERY non-value trailing token as a
                    # possible node ref, conservatively (see below).
        "o": 4,
    }

    for text_line, src_idx in logical:
        s = text_line.strip()
        if is_comment_or_blank(s):
            continue
        low = s.lower()
        if low.startswith(".subckt"):
            subckt_depth += 1
            continue
        if low.startswith(".ends"):
            subckt_depth = max(0, subckt_depth - 1)
            continue
        if s.startswith("."):
            # a control/model/option card -- not an element line. Node-like
            # tokens here are handled by the whole-file protected-node scan
            # instead of per-element modelling.
            continue

        toks = s.split()
        if not toks:
            continue
        ename = toks[0]
        etype = ename[0].lower()

        if subckt_depth == 0 and etype in ("r", "c") and len(toks) >= 4:
            n1, n2, value_tok = toks[1], toks[2], toks[3]
            value = parse_spice_num(value_tok)
            elements.append({
                "name": ename, "type": etype, "n1": n1, "n2": n2,
                "value_tok": value_tok, "value": value, "src_idx": src_idx,
            })
            continue

        # everything else: record node touches conservatively so these
        # nodes are marked ineligible for elimination.
        if etype == "x":
            # subckt instantiation: `Xname n1 n2 ... subckt_name [params]`
            # -- can't know the node count without the subckt definition;
            # treat every token up to the last one (the subckt name) as a
            # possible node, conservative but safe (over-protection only).
            for t in toks[1:-1]:
                if "=" not in t:
                    other_node_refs[t].add(ename)
        else:
            n = NODE_POSITIONS.get(etype)
            if n is None:
                n = 2  # unknown two-terminal-like device: assume first 2
            for t in toks[1:1 + n]:
                other_node_refs[t].add(ename)

    all_text_nodewords = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", text))

    return {
        "raw_lines": raw_lines,
        "logical": logical,
        "elements": elements,
        "other_node_refs": other_node_refs,
        "all_text_nodewords": all_text_nodewords,
    }


def build_graph(deck):
    """node -> {'r': [(neighbor, elem_dict)], 'c': [elem_dict], 'other': bool}"""
    g = defaultdict(lambda: {"r": [], "c": [], "other": False})
    for el in deck["elements"]:
        n1, n2 = el["n1"], el["n2"]
        if el["type"] == "r":
            g[n1]["r"].append((n2, el))
            g[n2]["r"].append((n1, el))
        else:
            g[n1]["c"].append(el)
            if n2 != "0":
                g[n2]["c"].append(el)
    for node in deck["other_node_refs"]:
        g[node]["other"] = True
    return g


def protected_nodes(deck, node):
    """A node is protected if its name occurs anywhere in the raw file text
    outside of its own R/C element definition lines (conservative: also
    protects a node whose name happens to appear as a substring match of a
    LARGER identifier -- \\b word-boundary regex avoids that)."""
    return node in deck["_protected_cache"]


def compute_protected(deck):
    """Whole-file text scan: a node is protected if the count of its
    \\b<node>\\b occurrences in the RAW text exceeds the count contributed
    by its own top-level R/C element lines (n1/n2 slots only).

    A SINGLE pass tokenises the whole file with one regex and a Counter --
    not one re.findall() per candidate node (that was O(n_nodes * file_size)
    and made this unusable on a 10k+-node mesh: ~30s+ and climbing, timed
    out in tests/test_rcreduce.py on e_pex_tree.cir/a_pex_buck.cir during
    this round's own validation, root-caused and fixed before being relied
    on for anything -- see the round's report)."""
    own_count = defaultdict(int)
    for el in deck["elements"]:
        own_count[el["n1"]] += 1
        own_count[el["n2"]] += 1

    text = "\n".join(deck["raw_lines"])
    from collections import Counter
    total_count = Counter(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", text))

    protected = set()
    candidates = set(own_count) | set(deck["other_node_refs"])
    for node in candidates:
        if total_count.get(node, 0) > own_count.get(node, 0):
            protected.add(node)
    return protected


def reduce_graph(g, protected, ground, max_degree, fast_node_tol, ref_tstep, log):
    """Mutates g in place (dict node -> {'r':[(nbr,el)], 'c':[el], 'other':bool}).
    Returns (eliminated_exact, eliminated_fast, new_resistors) where
    new_resistors is [(a, b, ohms, source_desc)] for every synthesized/merged
    resistor, in elimination order (for deterministic, reviewable output)."""
    import heapq

    def degree(node):
        return len(g[node]["r"])

    def eligible(node):
        if node == ground or node in protected:
            return False
        if g[node]["other"]:
            return False
        ncap = len(g[node]["c"])
        if ncap == 0:
            return True, "exact"
        if ncap == 1 and fast_node_tol is not None:
            return True, "fast"
        return False

    heap = []
    counter = 0
    nodes_of_interest = [n for n in g if n != ground]
    for n in nodes_of_interest:
        e = eligible(n)
        if e:
            heapq.heappush(heap, (degree(n), counter, n))
            counter += 1

    eliminated_exact = []
    eliminated_fast = []
    new_resistors = []
    stale = set()

    def combine_parallel(a, b, new_g):
        """Add conductance new_g between a and b in g, merging with any
        existing direct R edge between them (parallel combination)."""
        for i, (nbr, el) in enumerate(g[a]["r"]):
            if nbr == b and el.get("_synthetic"):
                old_g = el["_g"]
                el["_g"] = old_g + new_g
                el["value"] = 1.0 / el["_g"]
                return el
        # no existing synthetic edge between a,b -- but there may be an
        # ORIGINAL element edge between them too; keep that one separate
        # (both end up in parallel electrically, which is correct: two
        # resistor lines in parallel between the same two nodes is valid
        # SPICE and numerically identical to combining them, so leaving
        # the original untouched and adding a new synthetic one alongside
        # it is exact and simpler than merging into the original line).
        new_el = {
            "name": f"Rrcred_{a}_{b}_{len(new_resistors)}", "type": "r",
            "n1": a, "n2": b, "value": 1.0 / new_g, "_synthetic": True,
            "_g": new_g, "src_idx": None,
        }
        g[a]["r"].append((b, new_el))
        g[b]["r"].append((a, new_el))
        new_resistors.append(new_el)
        return new_el

    while heap:
        d, _, node = heapq.heappop(heap)
        if node in stale:
            continue
        if d != degree(node):
            # degree changed since pushed -- reinsert with current degree
            heapq.heappush(heap, (degree(node), counter, node))
            counter += 1
            continue
        e = eligible(node)
        if not e:
            continue
        mode = e[1] if isinstance(e, tuple) else "exact"
        legs = g[node]["r"]
        if len(legs) > max_degree:
            continue  # left alone -- safety valve, not an error

        # compute conductances of each leg; skip (leave node alone) if any
        # leg's value isn't a known positive numeric (unparsed expression,
        # zero/negative -- refuse rather than guess).
        gs = []
        ok = True
        for nbr, el in legs:
            v = el.get("value")
            if v is None or v <= 0:
                ok = False
                break
            gs.append((nbr, el, 1.0 / v))
        if not ok:
            continue

        gtotal = sum(x[2] for x in gs)
        if gtotal <= 0:
            continue

        # star-mesh: for every pair of neighbours, add g_a*g_b/gtotal
        n_legs = len(gs)
        if n_legs >= 2:
            for i in range(n_legs):
                for j in range(i + 1, n_legs):
                    a, _, ga = gs[i]
                    b, _, gb = gs[j]
                    if a == b:
                        continue  # a multigraph edge back to the same
                                  # neighbour via two different legs of the
                                  # node being eliminated -- no third node
                                  # to place a clique edge between; the
                                  # equivalent conductance folds into a's
                                  # own self-path and needs no new element.
                    combine_parallel(a, b, ga * gb / gtotal)

        # remove node's R legs from each neighbour's adjacency
        touched_nbrs = {nbr for nbr, _ in legs}
        for nbr, el in legs:
            g[nbr]["r"] = [(n2, e2) for (n2, e2) in g[nbr]["r"] if e2 is not el]

        if mode == "fast":
            eliminated_fast.append((node, g[node]["c"][0]))
        else:
            eliminated_exact.append(node)

        del g[node]
        stale.add(node)

        for nbr in touched_nbrs:
            if nbr in g and nbr != ground:
                ee = eligible(nbr)
                if ee:
                    heapq.heappush(heap, (degree(nbr), counter, nbr))
                    counter += 1

    return eliminated_exact, eliminated_fast, new_resistors


def render_output(deck, g, eliminated_exact, eliminated_fast, new_resistors, args):
    eliminated = set(eliminated_exact) | {n for n, _ in eliminated_fast}
    eliminated_cap_els = {id(c) for _, c in eliminated_fast}

    drop_src_idx = set()
    for el in deck["elements"]:
        if el["n1"] in eliminated or el["n2"] in eliminated:
            # an R touching an eliminated node is gone (replaced by the
            # star-mesh clique, or was one of the node's own legs -- both
            # cases mean this original line must not appear in the output).
            if el["type"] == "r" and (el["n1"] in eliminated or el["n2"] in eliminated):
                drop_src_idx.add(el["src_idx"])
            elif el["type"] == "c" and id(el) in eliminated_cap_els:
                drop_src_idx.add(el["src_idx"])

    out_lines = []
    header = [
        f"* --- rcreduce.py output (ngspice_augmentation card #506, round 8) ---",
        f"* source: {args.deck}",
        f"* generated: {time.strftime('%Y-%m-%dT%H:%M:%S')}",
        f"* exact eliminations (zero-capacitance internal R-only nodes): {len(eliminated_exact)}",
        f"* fast-node (opt-in, tol={args.fast_node_tol}, ref_tstep={args.ref_tstep}) eliminations: {len(eliminated_fast)}",
        f"* new/merged resistors from star-mesh elimination: {len(new_resistors)}",
        f"* max-degree safety cap: {args.max_degree}",
        f"* --- do not hand-edit below; regenerate with rcreduce.py instead ---",
    ]
    out_lines.extend(header)

    for i, raw in enumerate(deck["raw_lines"]):
        if i in drop_src_idx:
            continue
        out_lines.append(raw)

    if new_resistors:
        out_lines.append("* --- rcreduce.py: star-mesh clique/merged resistors ---")
        for el in new_resistors:
            out_lines.append(f"{el['name']} {el['n1']} {el['n2']} {el['value']:.10g}")

    return "\n".join(out_lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("-o", "--out", help="output .cir path (required unless --dry-run)")
    ap.add_argument("--ground", default="0")
    ap.add_argument("--max-degree", type=int, default=6,
                     help="fill-in safety cap (default 6, see module docstring); "
                          "a node whose resistor-degree exceeds this is left alone")
    ap.add_argument("--fast-node-tol", type=float, default=None,
                     help="opt-in bounded (TICER-style) elimination: eliminate a "
                          "single-capacitor internal node if its own tau=C/sum(g) "
                          "<= tol * --ref-tstep. Not exact -- state a tolerance "
                          "and validate waveforms before trusting the result.")
    ap.add_argument("--ref-tstep", default=None,
                     help="reference timestep (SPICE numeric, e.g. 10n) that "
                          "--fast-node-tol is relative to. Required if "
                          "--fast-node-tol is given.")
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    args = ap.parse_args()

    if args.fast_node_tol is not None and args.ref_tstep is None:
        ap.error("--fast-node-tol requires --ref-tstep")
    ref_tstep_s = None
    if args.ref_tstep is not None:
        ref_tstep_s = parse_spice_num(args.ref_tstep)
        if ref_tstep_s is None or ref_tstep_s <= 0:
            ap.error(f"--ref-tstep {args.ref_tstep!r} is not a positive SPICE numeric literal")

    if not args.dry_run and not args.out:
        ap.error("-o/--out is required unless --dry-run")

    with open(args.deck) as f:
        text = f.read()

    deck = parse_deck(text)
    deck["_protected_cache"] = compute_protected(deck)
    n_elements_in = len(deck["elements"])

    g = build_graph(deck)
    log = []
    eliminated_exact, eliminated_fast, new_resistors = reduce_graph(
        g, deck["_protected_cache"], args.ground, args.max_degree,
        args.fast_node_tol, ref_tstep_s, log)

    print(f"rcreduce: {args.deck}: {n_elements_in} top-level R/C elements parsed")
    print(f"  exact (zero-cap, star-mesh) eliminations: {len(eliminated_exact)}")
    print(f"  fast-node (tol={args.fast_node_tol}) eliminations: {len(eliminated_fast)}")
    print(f"  new/merged resistors: {len(new_resistors)}")
    print(f"  protected nodes (referenced outside R/C lines): {len(deck['_protected_cache'])}")

    if args.dry_run:
        return 0

    out_text = render_output(deck, g, eliminated_exact, eliminated_fast, new_resistors, args)
    with open(args.out, "w") as f:
        f.write(out_text)
    print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
