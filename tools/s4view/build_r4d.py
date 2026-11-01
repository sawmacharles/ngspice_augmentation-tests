#!/usr/bin/env python3
"""build_r4d.py -- the R4d peak-current-mode buck (refs/r4d_pcm_buck_deadtime.cir) as a schematics v4 project
(card #681), for looking at the design in the s4 UI.

    python3 tools/s4view/build_r4d.py [--root DIR] [--s4core DIR]

Writes the project r4d_pcm under --root (default ~/Work/data/ngspice_augmentation/s4view): a `devices` library
(one cell per symbol; bodies in tools/s4view/r4d_devices.lib) and a `pcm` library with the testbench cell r4d.
Connectivity is by pin_overrides (net per pin); wires and labels are drawn for reading only. Rebuilds from scratch
each run. It never writes outside --root.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = Path(os.path.expanduser("~/Work/data/ngspice_augmentation/s4view"))
DEFAULT_CORE = Path(os.path.expanduser("~/schematicsv4/app/core"))

# ---------------------------------------------------------------- symbols


def two_pin(names, graphics, params):
    return {"pins": [(names[0], (0, 5)), (names[1], (40, 5))], "box": (0, 0, 40, 10), "graphics": graphics,
            "params": params, "label": (-2, 16), "plabel": (2, -10)}


ZIG = [[0, 5], [10, 5], [13, 1], [16, 9], [19, 1], [22, 9], [25, 1], [28, 9], [30, 5], [40, 5]]


def coil():
    pts = [[0, 5], [6, 5]]
    for k in range(4):
        x0 = 6 + 7 * k
        pts += [[x0 + 1, 9], [x0 + 3.5, 11], [x0 + 6, 9], [x0 + 7, 5]]
    pts += [[40, 5]]
    return [("polyline", {"points": pts})]


SYMBOLS = {
    "resistor": two_pin(("p", "n"), [("polyline", {"points": ZIG})], [("r", "1k")]),
    "capacitor": two_pin(("p", "n"), [("line", {"points": [[0, 5], [18, 5]]}), ("line", {"points": [[18, -2], [18, 12]]}),
                                      ("line", {"points": [[22, -2], [22, 12]]}), ("line", {"points": [[22, 5], [40, 5]]})],
                         [("c", "1p"), ("ic", "0")]),
    "inductor": two_pin(("p", "n"), coil(), [("l", "1u"), ("ic", "0")]),
    "diode": two_pin(("a", "k"), [("line", {"points": [[0, 5], [14, 5]]}), ("polygon", {"points": [[14, 0], [14, 10], [24, 5]]}),
                                  ("line", {"points": [[24, 0], [24, 10]]}), ("line", {"points": [[24, 5], [40, 5]]})], []),
    "body_diode": two_pin(("a", "k"), [("line", {"points": [[0, 5], [14, 5]]}), ("polygon", {"points": [[14, 0], [14, 10], [24, 5]]}),
                                       ("line", {"points": [[24, 0], [24, 10]]}), ("line", {"points": [[24, 5], [40, 5]]}),
                                       ("rect", {"x": 10, "y": -4, "w": 18, "h": 18})],
                          [("g", "50"), ("vf", "0.6"), ("kd", "50")]),
    "vsource": two_pin(("n_plus", "n_minus"), [("line", {"points": [[0, 5], [10, 5]]}), ("circle", {"cx": 20, "cy": 5, "r": 10}),
                                               ("line", {"points": [[30, 5], [40, 5]]}), ("line", {"points": [[13, 2], [13, 8]]}),
                                               ("line", {"points": [[10, 5], [16, 5]]}), ("line", {"points": [[24, 5], [28, 5]]})],
                       [("dc", "0")]),
    "vpulse": two_pin(("p", "n"), [("line", {"points": [[0, 5], [10, 5]]}), ("circle", {"cx": 20, "cy": 5, "r": 10}),
                                   ("line", {"points": [[30, 5], [40, 5]]}),
                                   ("polyline", {"points": [[14, 1], [14, 8], [20, 8], [20, 1], [26, 1]]})],
                      [("v1", "0"), ("v2", "1"), ("td", "0"), ("tr", "1n"), ("tf", "1n"), ("pw", "1u"), ("per", "2u")]),
    "sw_smooth": {"pins": [("a", (0, 20)), ("b", (40, 20)), ("ctl", (20, 0))], "box": (0, 0, 40, 30),
                  "graphics": [("line", {"points": [[0, 20], [10, 20]]}), ("line", {"points": [[10, 20], [28, 28]]}),
                               ("line", {"points": [[30, 20], [40, 20]]}), ("circle", {"cx": 10, "cy": 20, "r": 1.5}),
                               ("circle", {"cx": 30, "cy": 20, "r": 1.5}), ("line", {"points": [[20, 0], [20, 22]]})],
                  "params": [("gon", "70"), ("goff", "1e-6"), ("k", "150")], "label": (24, 8), "plabel": (24, 2)},
}


def block(pins_in, pin_out, h, params, shape="box"):
    w = 60
    g = [("rect", {"x": 10, "y": 0, "w": w - 20, "h": h})] if shape == "box" else \
        [("polygon", {"points": [[10, 0], [10, h], [w - 10, h / 2]]})]
    pins = []
    for name, y in pins_in:
        pins.append((name, (0, y)))
        g.append(("line", {"points": [[0, y], [10, y]]}))
    pins.append((pin_out, (w, h / 2)))
    g.append(("line", {"points": [[w - 10, h / 2], [w, h / 2]]}))
    return {"pins": pins, "box": (0, 0, w, h), "graphics": g, "params": params, "label": (12, h + 6),
            "plabel": (12, -10), "cellname": (14, h / 2 - 3) if shape == "box" else None}


SYMBOLS.update({
    "csamp": block([("p", 30), ("n", 10)], "out", 40, [("gain", "8")]),
    "cmp_sum": block([("a", 50), ("b", 30), ("c", 10)], "out", 60, [("k", "100")]),
    "rs_latch": block([("s", 40), ("r", 20)], "q", 60, [("cq", "1p"), ("gset", "1e-4"), ("greset", "1")]),
    "ota_gm": block([("inp", 30), ("inn", 10)], "out", 40, [("gm", "1000")], shape="tri"),
    "vbuf": block([("in", 10)], "out", 20, [], shape="tri"),
    "vinv": block([("in", 10)], "out", 20, [], shape="tri"),
})
# the output pin is listed last in block(); netlist pin order must be out first (matches r4d_devices.lib)
for _n in ("csamp", "cmp_sum", "rs_latch", "ota_gm", "vbuf", "vinv"):
    p = SYMBOLS[_n]["pins"]
    SYMBOLS[_n]["pins"] = [p[-1]] + p[:-1]
SYMBOLS["vinv"]["graphics"].append(("circle", {"cx": 52, "cy": 10, "r": 2}))

# the one parameter each symbol shows on the sheet (all of them are in the instance's properties)
SHOW = {"sw_smooth": ["gon"], "vpulse": ["v2"], "rs_latch": [], "body_diode": ["vf"]}

# ---------------------------------------------------------------- placement
# (name, cell, pins-in-symbol-order -> nets, params, placement)
# placement: ("h", x, y) horizontal, first pin at (x, y), second at (x+40, y); ("hr", x, y) reversed: first pin at
# (x, y), second at (x-40, y); ("v", x, y) vertical, first pin at (x, y), second at (x, y-40); ("vu", x, y) first pin
# at (x, y), second at (x, y+40); ("at", x, y, rot) raw position.

E = [
    # power stage
    ("Vin", "vsource", ["vin", "gnd"], {"dc": "11"}, ("v", 0, 400)),
    ("XHS", "sw_smooth", ["vin", "sw", "dly_hs"], {"gon": "70", "goff": "1e-6", "k": "150"}, ("at", 40, 380, 0)),
    ("XLS", "sw_smooth", ["sw", "gnd", "dly_ls"], {"gon": "70", "goff": "1e-6", "k": "150"}, ("at", 100, 360, 270)),
    ("XBODY", "body_diode", ["gnd", "sw"], {"g": "50", "vf": "0.6", "kd": "50"}, ("vu", 150, 300)),
    ("XCSW", "capacitor", ["sw", "gnd"], {"c": "100p", "ic": "0"}, ("v", 190, 360)),
    ("XRSENSE", "resistor", ["sw", "swcs"], {"r": "0.05"}, ("h", 120, 400)),
    ("XL1", "inductor", ["swcs", "swl"], {"l": "10u", "ic": "0"}, ("h", 180, 400)),
    ("XRDCR", "resistor", ["swl", "vout"], {"r": "10m"}, ("h", 240, 400)),
    ("XRESR", "resistor", ["vout", "cx"], {"r": "10m"}, ("v", 320, 360)),
    ("XCOUT", "capacitor", ["cx", "gnd"], {"c": "100u", "ic": "0"}, ("v", 320, 310)),
    ("XRLOAD", "resistor", ["vout", "gnd"], {"r": "3.3"}, ("v", 370, 360)),
    ("Vprobe", "vsource", ["fbtop", "vout"], {"dc": "0"}, ("hr", 460, 400)),
    ("XRA", "resistor", ["fbtop", "fbnode"], {"r": "2300"}, ("v", 500, 400)),
    ("XRB", "resistor", ["fbnode", "gnd"], {"r": "1000"}, ("v", 500, 350)),
    # voltage loop (type II)
    ("XR1IN", "resistor", ["fbnode", "compin"], {"r": "10k"}, ("h", 500, 190)),
    ("Vrefsrc", "vsource", ["vrefnode", "gnd"], {"dc": "1"}, ("v", 560, 260)),
    ("XOTA", "ota_gm", ["vc", "vrefnode", "compin"], {"gm": "1000"}, ("at", 580, 180, 0)),
    ("XC1F", "capacitor", ["compin", "vc"], {"c": "25.21p", "ic": "0"}, ("h", 590, 140)),
    ("XR2F", "resistor", ["compin", "n2"], {"r": "50.51k"}, ("h", 560, 100)),
    ("XC2F", "capacitor", ["n2", "vc"], {"c": "6.534n", "ic": "0"}, ("h", 620, 100)),
    ("XROA", "resistor", ["vc", "gnd"], {"r": "100"}, ("v", 700, 200)),
    ("XCOA", "capacitor", ["vc", "gnd"], {"c": "100p", "ic": "0"}, ("v", 740, 200)),
    # current loop: sense, slope, comparator, latch
    ("XCS", "csamp", ["vcs", "sw", "swcs"], {"gain": "8"}, ("at", 0, 0, 0)),
    ("XVSLOPE", "vpulse", ["slope", "gnd"], {"v1": "0", "v2": "0", "td": "0", "tr": "1.94u", "tf": "60n", "pw": "1n",
                                              "per": "2u"}, ("v", 110, -30)),
    ("XCMP", "cmp_sum", ["reset_sig", "vcs", "slope", "vc"], {"k": "100"}, ("at", 160, -10, 0)),
    ("XVSETP", "vpulse", ["setp", "gnd"], {"v1": "0", "v2": "1", "td": "0", "tr": "1n", "tf": "1n", "pw": "20n",
                                            "per": "2u"}, ("v", 260, -30)),
    ("XLATCH", "rs_latch", ["q", "setp", "reset_sig"], {"cq": "1p", "gset": "1e-4", "greset": "1"}, ("at", 300, -10, 0)),
    ("XQBUF", "vbuf", ["qbuf", "q"], {}, ("at", 420, 10, 0)),
    ("XQINV", "vinv", ["qbuf_inv", "q"], {}, ("at", 420, -100, 0)),
    # dead-time network, high side
    ("XRCHG_HS", "resistor", ["qbuf", "dly_hs"], {"r": "430"}, ("h", 500, 20)),
    ("XDDIS_HS", "diode", ["dly_hs", "dis_hs"], {}, ("hr", 580, 60)),
    ("XRDIS_HS", "resistor", ["dis_hs", "qbuf"], {"r": "60"}, ("hr", 540, 60)),
    ("XCDLY_HS", "capacitor", ["dly_hs", "gnd"], {"c": "100p", "ic": "0"}, ("v", 580, 20)),
    # dead-time network, low side
    ("XRCHG_LS", "resistor", ["qbuf_inv", "dly_ls"], {"r": "430"}, ("h", 500, -90)),
    ("XDDIS_LS", "diode", ["dly_ls", "dis_ls"], {}, ("hr", 580, -50)),
    ("XRDIS_LS", "resistor", ["dis_ls", "qbuf_inv"], {"r": "60"}, ("hr", 540, -50)),
    ("XCDLY_LS", "capacitor", ["dly_ls", "gnd"], {"c": "100p", "ic": "0"}, ("v", 580, -90)),
]

NOTES = [
    ("R4d peak-current-mode buck: 11 V -> 3.3 V, 500 kHz, Rload 3.3 ohm (refs/r4d_pcm_buck_deadtime.cir)", (0, 470)),
    ("POWER STAGE", (0, 440)),
    ("VOLTAGE LOOP (type II, Vref 1 V)", (500, 280)),
    ("CURRENT LOOP: sense x8, slope, comparator, RS latch", (0, 80)),
    ("DEAD TIME (break-before-make)", (620, 40)),
    (".pstb probe: Vprobe (fbtop - vout). Analyses are not stored in this view; see tools/s4view/build_r4d.py.", (0, -160)),
]


def rot(dx, dy, r):
    if r == 90:
        return -dy, dx
    if r == 180:
        return -dx, -dy
    if r == 270:
        return dy, -dx
    return dx, dy


def place(cell, how):
    """(position, rotation) for a placement spec"""
    pins = SYMBOLS[cell]["pins"]
    kind = how[0]
    if kind == "at":
        return (how[1], how[2]), how[3]
    x, y = how[1], how[2]
    if kind == "h":
        return (x, y - 5), 0
    if kind == "hr":
        return (x, y + 5), 180
    if kind == "v":   # rot 270: (0,5)->(5,0); (40,5)->(5,-40)
        return (x - 5, y), 270
    if kind == "vu":  # rot 90: (0,5)->(-5,0); (40,5)->(-5,40)
        return (x + 5, y), 90
    raise ValueError(kind)


def pin_world(pos, r, local):
    dx, dy = rot(local[0], local[1], r)
    return (pos[0] + dx, pos[1] + dy)


def build(root, core):
    sys.path.insert(0, str(core))
    from s4core import io  # noqa: E402
    from s4core.model import (CellParams, Config, Graphic, Instance, InstanceRef, NetLabel, ParamDef, Pin, Port,  # noqa
                              PropertyText, SchematicView, Sheet, SymbolView, TextNote, ViewChoice, Wire)

    pdir = Path(root) / "r4d_pcm"
    if pdir.exists():
        shutil.rmtree(pdir)
    pdir.mkdir(parents=True)
    project = io.create_project(pdir, name="r4d_pcm")
    io.create_library(pdir, project, "devices", kind="work", priority=0)
    for name, s in SYMBOLS.items():
        cp = CellParams(parameters=[ParamDef(name=k, default=v, unit="") for k, v in s["params"]])
        cdir = io.create_cell(pdir, project, "devices", name, cp)
        props = [PropertyText(bind="system:instance_name", position=s["label"])]
        if s.get("cellname"):
            props.append(PropertyText(bind="system:cell_name", position=s["cellname"]))
        shown = SHOW.get(name, [k for k, _ in s["params"][:1]])
        for k in shown:
            props.append(PropertyText(bind=f"param:{k}", position=s["plabel"]))
        pins = [Pin(name=n, direction="out", position=p) for n, p in s["pins"]]
        io.save_symbol_view(cdir, "symbol", SymbolView(origin=(0, 0), boundary=s["box"], pins=pins,
                                                       graphics=[Graphic(k, a) for k, a in s["graphics"]],
                                                       properties=props))
    io.create_library(pdir, project, "pcm", kind="work", priority=1)
    cdir = io.create_cell(pdir, project, "pcm", "r4d", CellParams())

    instances, pinpos = [], {}
    for name, cell, nets, params, how in E:
        pos, r = place(cell, how)
        pins = SYMBOLS[cell]["pins"]
        assert len(pins) == len(nets), name
        ov = {}
        for (pn, local), net in zip(pins, nets):
            ov[pn] = {"net": net}
            pinpos.setdefault(net, []).append(pin_world(pos, r, local))
        instances.append(Instance(name=name, ref=InstanceRef(cell=cell, symbol_view="symbol", library="devices"),
                                  position=pos, rotation=r, mirror="none", params=params, pin_overrides=ov))
    wires, labels = [], []
    for net, pts in sorted(pinpos.items()):
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        if net != "gnd" and max(xs) - min(xs) <= 100 and max(ys) - min(ys) <= 100 and len(pts) > 1:
            pts = sorted(set(pts))
            for a, b in zip(pts, pts[1:]):
                # vertical leg first, horizontal into the pin: wires never run along a block's input edge
                path = [a, b] if a[0] == b[0] or a[1] == b[1] else [a, (a[0], b[1]), b]
                for u, v in zip(path, path[1:]):
                    wires.append(Wire(net=net, points=[u, v]))
            labels.append(NetLabel(net=net, position=(pts[0][0] + 2, pts[0][1] + 3)))
        else:
            for p in sorted(set(pts)):
                labels.append(NetLabel(net=net, position=(p[0] + 2, p[1] + 3)))
    notes = [TextNote(text=t, position=p) for t, p in NOTES]
    sheet = Sheet(name="1", instances=instances, wires=wires, labels=labels, notes=notes, graphics=[])
    ports = [Port(name="gnd", direction="power", position=(0, -200), power_polarity="negative", is_reference=True)]
    io.save_schematic_view(cdir, "testbench", SchematicView(kind="testbench", ports=ports, sheets=[sheet]))
    io.save_config(pdir, Config(name="default", view_priority=["schematic", "testbench", "veriloga"],
                                cell_overrides={f"devices/{n}": ViewChoice(disposition="opaque") for n in SYMBOLS}))
    return pdir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(DEFAULT_ROOT))
    ap.add_argument("--s4core", default=str(DEFAULT_CORE))
    a = ap.parse_args()
    print(build(a.root, Path(a.s4core)))


if __name__ == "__main__":
    main()
