#!/usr/bin/env python3
"""gen_pex_buck.py -- generate bench/decks/a_pex_buck.cir (ngspice_augmentation
card #467, plan T9.1, bench task 1a).

A "PEX-sized" synchronous buck: same power-stage topology as refs/r3_vm_buck.cir
(12V->3.3V, fs=500kHz, type-III compensator, smooth PWM comparator -- see
that file's header for the design derivation, copied here structurally) but
with two changes to make it look like the real BSIM4-PEX decks the archived
ngspice_accel project profiled (device evaluation dominates the load phase --
carried-over fact, research/ngspice-internals-pss-pac.md):

1. The high-side and low-side switches are each split into HS_FINGERS /
   LS_FINGERS parallel BSIM4 (level=54) instances, each with its own small
   parasitic gate/drain/source resistor to a common bus node -- a plausible
   stand-in for a multi-finger power MOSFET layout extracted from a real PEX
   flow: many device instances, each individually cheap but BSIM4-expensive
   to evaluate. Gate drive is an ideal bootstrap-style E-source (not itself
   BSIM4 -- the driver is not the device-eval target here).
2. A synthetic parasitic decoupling mesh (an N x N grid of R/C unit cells,
   trace resistance + plane capacitance to ground) is coupled onto the vout
   rail through one resistor. This is the "thousands of interconnect
   elements" that make the linear system big (target ~10-20k unknowns) while
   staying sparse and not participating in device evaluation -- matching how
   a real PEX netlist's node count is dominated by extracted parasitics, not
   by the transistor count.

All values here are synthetic / plausible, not derived from or copied out of
any real design (#454 A).

Output: bench/decks/a_pex_buck.cir (checked in -- it is a netlist, not
simulation data; regenerate with this script rather than hand-editing).
Simulation output goes to
~/Work/data/ngspice_augmentation/accel/bench/a_pex_buck/ when run, never
into the repo.
"""
import argparse
import pathlib

HERE = pathlib.Path(__file__).resolve().parent


def mesh_lines(n, rtrace, cplane):
    lines = [f"* {n}x{n} RC mesh: rtrace={rtrace} between neighbours, cplane={cplane} to ground per node\n"]
    for i in range(n):
        for j in range(n):
            node = f"mesh_{i}_{j}"
            lines.append(f"Cp_{i}_{j} {node} 0 {cplane}\n")
            if j + 1 < n:
                lines.append(f"Rx_{i}_{j} {node} mesh_{i}_{j+1} {rtrace}\n")
            if i + 1 < n:
                lines.append(f"Ry_{i}_{j} {node} mesh_{i+1}_{j} {rtrace}\n")
    return lines


def tree_lines(branch, depth, rtrace0, rgrow, cplane0, cshrink,
               local_mesh_stride, coupling_caps, coupling_c, sparse_cap=False):
    """A tree-shaped synthetic parasitic network (card #478 T9.2 task 2),
    standing in for a real PEX-extracted net far better than
    mesh_lines()'s regular grid does (research/transient-acceleration.md
    section 2 / bench/README.md's "known limitation" note: a real
    extraction tool builds per-net RC TREES from the routing, not a
    lattice, and orders far better for a sparse direct solver -- a regular
    grid is close to a worst case for AMD/COLAMD-family fill-in without
    nested-dissection).

    Structure, BFS by level from one root (the coupling point onto vout):
      - `branch`-way fan-out per level, `depth` levels deep (level 0 = root).
      - Each tree EDGE (parent -> child) is a series trace resistor; each
        NODE gets a capacitor to ground (plane/junction capacitance stand-
        in). Both scale with depth: resistance grows (`rgrow`x per level,
        thinner wire deeper into the net) and capacitance shrinks
        (`cshrink`x per level, smaller area deeper into the net) --
        the same qualitative shape as a real extracted trunk-to-leaf RC
        tree (Elmore-delay style), not uniform like the grid.
      - "local meshes": every `local_mesh_stride`'th pair of adjacent
        siblings under the same parent gets an extra resistor directly
        between them (not through the parent) -- a small loop, standing in
        for local routing meshes (e.g. a strapped power grid segment) that
        a pure tree wouldn't have. Adds a bounded, deterministic number of
        non-tree edges (does not change the node count).
      - "coupling caps": a deterministic, evenly-spaced set of small
        capacitors between nodes in DIFFERENT subtrees at the same depth
        (adjacent-track capacitive coupling stand-in), count set by
        `coupling_caps`.

    `sparse_cap` (card #508, round 9, plan T9.9 -- round 8's own
    `research/transient-acceleration.md` section 13.1.4/13.1.5 flagged this
    exact change as "the clearest concrete next step"): when True, a node
    gets its `Cp_*` capacitor to ground ONLY at EVEN depths (0, 2, 4, ...)
    -- odd-depth nodes are resistor-only (no capacitor at all), mimicking a
    real PEX extraction's via/contact layer (no explicit capacitance,
    round 8's own `iic-jku/klayout-pex` finding) alternating with a
    capacitive metal-plane layer. Default False reproduces this function's
    ORIGINAL behaviour byte-for-byte (a capacitor on every node, round 7's
    own `e_pex_tree.cir`) -- additive, not a change to any existing bench
    deck. Round 8 independently found branch=3/depth=7 with this exact
    alternating pattern gives a genuine ~40% net element-count reduction
    under `tools/rcreduce.py`'s exact pass (2459 eliminations, no net
    fill-in increase), while branch=3/depth=8 reproduces the grid's own
    fill-in pathology instead (deeper odd levels are no longer the leaf
    level, so eliminating them creates fill-in edges directly between
    surviving even-level nodes) -- `--tree-depth` stays a caller choice,
    this function does not clamp it, but see `bench/decks/h_pex_tree_
    sparse.cir`'s own header for the round-9-recommended depth=7.
    Fully deterministic (no RNG/seed) so the generator's output is
    reproducible from its arguments alone, same as mesh_lines().
    Returns (lines, root_node_name, deepest_last_node_name, total_nodes).
    """
    lines = [f"* tree-shaped synthetic PEX net: branch={branch} depth={depth}, "
             f"rtrace0={rtrace0} (x{rgrow}/level), cplane0={cplane0} (x{cshrink}/level), "
             f"local-mesh-stride={local_mesh_stride}, coupling-caps={coupling_caps}\n"]

    # BFS build: level 0 is [root]; level k+1 is each level-k node's `branch` children.
    node_name = lambda nid: f"tree_{nid}"
    levels = [[0]]          # levels[d] = list of node ids at depth d
    parent_of = {0: None}
    next_id = 1
    lines.append(f"Cp_0 {node_name(0)} 0 {cplane0:.6g}\n")
    for d in range(depth):
        rtrace = rtrace0 * (rgrow ** d)
        cplane = cplane0 * (cshrink ** (d + 1))
        this_level = []
        for parent in levels[d]:
            for _ in range(branch):
                cid = next_id
                next_id += 1
                parent_of[cid] = parent
                this_level.append(cid)
                lines.append(f"Rtr_{parent}_{cid} {node_name(parent)} {node_name(cid)} {rtrace:.6g}\n")
                if not sparse_cap or (d + 1) % 2 == 0:
                    lines.append(f"Cp_{cid} {node_name(cid)} 0 {cplane:.6g}\n")
        levels.append(this_level)

    # local meshes: within each parent's sibling group, chain every
    # local_mesh_stride'th adjacent pair directly (loop, not through parent).
    if local_mesh_stride > 0:
        for d in range(1, len(levels)):
            siblings_by_parent = {}
            for nid in levels[d]:
                siblings_by_parent.setdefault(parent_of[nid], []).append(nid)
            rtrace = rtrace0 * (rgrow ** (d - 1))
            for parent, sibs in siblings_by_parent.items():
                for i in range(0, len(sibs) - 1, local_mesh_stride):
                    a, b = sibs[i], sibs[i + 1]
                    lines.append(f"Rlm_{a}_{b} {node_name(a)} {node_name(b)} {rtrace * 4:.6g}\n")

    # coupling caps: deterministic evenly spaced pairs of nodes at the same
    # depth, offset into a different parent's subtree.
    if coupling_caps > 0:
        deepest = levels[-1]
        m = len(deepest)
        if m >= 2:
            stride = max(1, m // coupling_caps)
            k = 0
            for i in range(0, m - stride, stride):
                a, b = deepest[i], deepest[i + stride]
                lines.append(f"Ccpl_{a}_{b} {node_name(a)} {node_name(b)} {coupling_c:.6g}\n")
                k += 1
                if k >= coupling_caps:
                    break

    total_nodes = next_id
    return lines, node_name(0), node_name(levels[-1][-1]), total_nodes


def finger_lines(prefix, n, drain_bus, gate_bus, source_bus, model, w_total, l, rfing):
    lines = []
    w_each = w_total / n
    for k in range(n):
        d, g, s = f"{prefix}_d{k}", f"{prefix}_g{k}", f"{prefix}_s{k}"
        lines.append(f"M{prefix}{k} {d} {g} {s} {s} {model} w={w_each:.6g} l={l}\n")
        lines.append(f"Rfd{prefix}{k} {drain_bus} {d} {rfing}\n")
        lines.append(f"Rfg{prefix}{k} {gate_bus} {g} {rfing}\n")
        lines.append(f"Rfs{prefix}{k} {source_bus} {s} {rfing}\n")
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--topology", choices=["grid", "tree"], default="grid",
                     help="grid (default, a_pex_buck.cir's original regular NxN mesh) or "
                          "tree (card #478 T9.2 task 2: branches/local meshes/coupling caps, "
                          "closer to a real PEX-extracted net -- see tree_lines() docstring)")
    ap.add_argument("--mesh-n", type=int, default=110, help="[grid] mesh grid side (n*n nodes)")
    ap.add_argument("--tree-branch", type=int, default=10, help="[tree] fan-out per level")
    ap.add_argument("--tree-depth", type=int, default=4, help="[tree] levels below the root")
    ap.add_argument("--tree-rtrace0", type=float, default=0.2, help="[tree] root-level series R (ohm)")
    ap.add_argument("--tree-rgrow", type=float, default=1.6, help="[tree] R multiplier per level deeper")
    ap.add_argument("--tree-cplane0", type=float, default=20e-15, help="[tree] root node's own C to ground (F)")
    ap.add_argument("--tree-cshrink", type=float, default=0.6, help="[tree] C multiplier per level deeper")
    ap.add_argument("--tree-local-mesh-stride", type=int, default=3,
                     help="[tree] extra sibling-to-sibling R every Nth adjacent pair (0 disables)")
    ap.add_argument("--tree-coupling-caps", type=int, default=40,
                     help="[tree] count of cross-subtree coupling caps at the deepest level")
    ap.add_argument("--tree-coupling-c", type=float, default=2e-15, help="[tree] coupling cap value (F)")
    ap.add_argument("--tree-sparse-cap", action="store_true",
                     help="[tree] card #508 round 9: only even-depth nodes get a Cp_* capacitor "
                          "(odd depths are resistor-only, a via/contact-layer stand-in) -- makes "
                          "tools/rcreduce.py's exact star-mesh pass find eliminable nodes on this "
                          "generator's own output (round 8 found the default, capacitor-on-every-"
                          "node shape has zero eliminable nodes). Default off: unchanged output.")
    ap.add_argument("--hs-fingers", type=int, default=30)
    ap.add_argument("--ls-fingers", type=int, default=30)
    ap.add_argument("--tstop-us", type=float, default=200.0)
    ap.add_argument("--tstep", default="10n")
    ap.add_argument("--out", default=str(HERE / "decks" / "a_pex_buck.cir"))
    args = ap.parse_args()

    n = args.mesh_n
    tstop_us = args.tstop_us
    if args.topology == "tree":
        b, d = args.tree_branch, args.tree_depth
        nodes = (b ** (d + 1) - 1) // (b - 1) if b > 1 else d + 1
        parasitic_desc = (
            f"a tree-shaped synthetic PEX net ({nodes} nodes: branch={b}, depth={d}, "
            f"plus local-mesh loops and cross-subtree coupling caps -- see "
            f"bench/gen_pex_buck.py's tree_lines() docstring) is coupled onto vout"
        )
    else:
        nodes = n * n
        parasitic_desc = f"a {n}x{n}={nodes}-node synthetic RC decoupling mesh is coupled onto vout"

    deck_name = pathlib.Path(args.out).stem
    L = []
    L.append(f"""\
* {deck_name} -- PEX-sized synchronous buck, synthetic
* (ngspice_augmentation card #467/#478, plan T9.1/T9.2; generated by
* bench/gen_pex_buck.py --topology {args.topology} -- do not hand-edit,
* regenerate instead). Not derived from any real deck.
*
* Power stage / compensator topology as refs/r3_vm_buck.cir (12V->3.3V,
* fs=500kHz; see that file for the full compensator design derivation) but:
*  - HS/LS switches are {args.hs_fingers}/{args.ls_fingers}-finger BSIM4 (level=54)
*    arrays instead of smooth B-source ideal switches (device-eval load).
*  - {parasitic_desc}
*    (matrix-size load; ~{nodes} extra unknowns, sparse).
* run: ngspice -b <this file>  (writes to
*  ~/Work/data/ngspice_augmentation/accel/bench/<deck>/, never into the repo)

.param vin_v=12 fs=500k tclk={{1/fs}}
.param lval=10u rdcr=10m cval=100u esrval=10m rload=3.3
.param ra=2300 rb=1000 vref=1
.param r1=10k r2=27.51k r3=101.7 c1=36.35p c2=2.299n c3=6.261n
.param vramp=2
.param gm_oa=1000 ro_oa=100 co_oa=100p
.param kc=30
.param vgate=5

* --- plausible inline BSIM4 model, NOT calibrated to any real process
* (deliberately simple: capmod=0, most extra-term flags off, same
* minimisation style as tests/fixtures/stb3_bsim4_cgd.cir) ---
.model nch4 nmos level=54 version=4.8.2
+ mobmod=0 capmod=0 igcmod=0 igbmod=0 rdsmod=0 rbodymod=0 rgatemod=0 acnqsmod=0
+ toxe=5e-9 toxp=5e-9 toxm=5e-9
+ vth0=0.7 k1=0.4 k2=0 k3=0 w0=2.5e-6
+ u0=300 vsat=1.1e5
+ cgso=2e-10 cgdo=2e-10 cgbo=0 cgsl=0 cgdl=0
+ cf=0 cjs=1e-4 cjd=1e-4 cjswd=1e-10 cjswgd=1e-10

Vin vin 0 DC {{vin_v}}

* --- PWM comparator (as r3): pwm in [0,1], high = HS on ---
Vramp ramp 0 PULSE(0 {{vramp}} 0 {{0.97*tclk}} {{0.03*tclk}} 1n {{tclk}})
Bcmp pwm 0 V = 0.5*(1+tanh({{kc}}*(V(opout)-V(ramp))))
Bpwmls pwm_ls 0 V = 1-V(pwm)

* --- gate drive: ideal bootstrap-style level shifters (not BSIM4 -- the
* driver itself is not the device-eval target here) ---
Ehsgate hsgate_bus sw pwm 0 {{vgate}}
Elsgate lsgate_bus 0 pwm_ls 0 {{vgate}}

* --- HS switch: {args.hs_fingers}-finger BSIM4 array, drain=vin bus, source=sw bus ---
""")
    L += finger_lines("hs", args.hs_fingers, "vin", "hsgate_bus", "sw", "nch4", 0.06, "0.35u", "0.02")
    L.append(f"\n* --- LS switch: {args.ls_fingers}-finger BSIM4 array, drain=sw bus, source=gnd ---\n")
    L += finger_lines("ls", args.ls_fingers, "sw", "lsgate_bus", "0", "nch4", 0.06, "0.35u", "0.02")

    L.append("""
Csw sw 0 100p
L1 sw swl {lval} IC=0
Rdcr swl vout {rdcr}
Resr vout cx {esrval}
Cout cx 0 {cval} IC=0
Rload vout 0 {rload}

* --- feedback divider + type-III error amp (as r3, inline op-amp) ---
Ra vout fbnode {ra}
Rb fbnode 0 {rb}

Vrefsrc vrefnode 0 DC {vref}
Goa opout 0 vrefnode compin {gm_oa}
Roa opout 0 {ro_oa}
Coa opout 0 {co_oa}

R1in fbnode compin {r1}
R3in fbnode n3 {r3}
C3in n3 compin {c3}
C1f compin opout {c1}
R2f compin n2 {r2}
C2f n2 opout {c2}
""")

    tstop = f"{tstop_us}u"
    t_last10_from = f"{tstop_us*0.9}u"
    t_mid_from, t_mid_to = f"{tstop_us*0.45}u", f"{tstop_us*0.5}u"
    t_25, t_50, t_75 = f"{tstop_us*0.25}u", f"{tstop_us*0.5}u", f"{tstop_us*0.75}u"

    if args.topology == "tree":
        L.append(f"\n* --- synthetic PEX net: tree topology, {nodes} nodes, coupled onto vout ---\n")
        tlines, root_node, far_node, actual_nodes = tree_lines(
            args.tree_branch, args.tree_depth,
            args.tree_rtrace0, args.tree_rgrow,
            args.tree_cplane0, args.tree_cshrink,
            args.tree_local_mesh_stride, args.tree_coupling_caps, args.tree_coupling_c,
            sparse_cap=args.tree_sparse_cap)
        assert actual_nodes == nodes, f"node-count mismatch: predicted {nodes}, built {actual_nodes}"
        L.append(f"Rtreecouple vout {root_node} 10\n")
        L += tlines
        far_point_meas = f"meas tran wave_mesh_far find v({far_node}) at={tstop}\n"
    else:
        L.append(f"\n* --- synthetic PEX decoupling mesh: {n}x{n}={nodes} nodes, coupled onto vout ---\n")
        L.append(f"Rmeshcouple vout mesh_0_0 10\n")
        L += mesh_lines(n, rtrace="1", cplane="1p")
        far_point_meas = f"meas tran wave_mesh_far find v(mesh_{n-1}_{n-1}) at={tstop}\n"

    L.append(f"""
.tran {args.tstep} {tstop} 0 uic
.options reltol=1e-3 abstol=1e-9 vntol=1e-6 klu

.control
set filetype=ascii
rusage cputime
run
rusage cputime
meas tran vout_avg_final AVG v(vout) FROM={t_last10_from} TO={tstop}
meas tran vout_max_final MAX v(vout) FROM={t_last10_from} TO={tstop}
meas tran vout_min_final MIN v(vout) FROM={t_last10_from} TO={tstop}
let vout_ripple = vout_max_final - vout_min_final
meas tran duty_meas AVG v(pwm) FROM={t_mid_from} TO={t_mid_to}
meas tran wave_t1 find v(vout) at={t_25}
meas tran wave_t2 find v(vout) at={t_50}
meas tran wave_t3 find v(vout) at={t_75}
meas tran wave_t4 find v(vout) at={tstop}
meas tran wave_sw_t1 find v(sw) at={t_50}
{far_point_meas}print vout_avg_final vout_max_final vout_min_final vout_ripple duty_meas wave_t1 wave_t2 wave_t3 wave_t4 wave_sw_t1 wave_mesh_far
rusage cputime
wrdata $HOME/Work/data/ngspice_augmentation/accel/bench/{deck_name}/wave.txt v(vout) v(sw) i(l1)
quit
.endc
.end
""")

    outpath = pathlib.Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text("".join(L))
    print(f"wrote {outpath} ({nodes} mesh nodes, {args.hs_fingers}+{args.ls_fingers} BSIM4 fingers, tstop={tstop})")


if __name__ == "__main__":
    main()
