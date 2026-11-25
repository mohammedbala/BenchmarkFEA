"""Multi-bay multi-storey moment frames (plane model in 3D space, out-of-plane DOFs restrained).

Rigid beam-column joints, fixed / pinned / partially fixed bases, floor loads in two pattern cases
(alternate bays), wind as storey forces.  Optionally second order, one bay with pinned beams (a
gravity bay leaning on the moment frame) and rigid end offsets of the beams (only without second
order).
"""
from __future__ import annotations

from ._common import (BUCKLING_CHANCE, BUCKLING_FLAG, BRACING_PATTERNS, Builder, add_combinations, flag, grid_lines, line_id, service_ids,
                      storey_levels)

NAME = "multistorey_2d"
LEVELS = (2, 3, 4)
FEATURES = {
    "buckling": BUCKLING_FLAG,
    "secondOrder": "strength combinations analysed second order",
    "pinnedBay": "one bay whose beams are pinned at both ends (gravity bay)",
    "baseSprings": "partially fixed bases (rotational springs)",
    "rigidOffsets": "rigid end offsets of the beams at the columns (only without second order)",
}
_SHAPES = {2: ((1, 2), (2, 1), (1, 2)), 3: ((2, 2), (2, 3), (3, 2), (1, 3), (1, 4)), 4: ((3, 4), (4, 3), (3, 5), (4, 4), (2, 5))}


def build(draw, level: int, unit_system: str, vertical: str, features: dict) -> Builder:
    B = Builder(unit_system, vertical, NAME)
    plane = draw.choice(("xv", "yv"))
    nb, ns = draw.choice(_SHAPES[level])
    xs = grid_lines(B, draw, nb, 5.0, 9.0, 0.5, equal=draw.chance(0.6))
    vs = storey_levels(B, draw, ns)
    base = draw.choice(("fixed", "pinned"))
    second = flag(features, "secondOrder", draw.chance({2: 0.0, 3: 0.8, 4: 0.9}[level]))
    pinned_bay = flag(features, "pinnedBay", draw.chance({2: 0.0, 3: 0.3, 4: 0.7}[level])) and nb >= 2
    springs = flag(features, "baseSprings", draw.chance({2: 0.0, 3: 0.3, 4: 0.4}[level]))
    offsets = flag(features, "rigidOffsets", draw.chance({2: 0.2, 3: 0.3, 4: 0.3}[level])) and not second
    n_combos = {2: draw.integer(3, 5), 3: draw.integer(4, 7), 4: draw.integer(6, 8)}[level]
    spacing = draw.uniform(5.0, 8.0)
    qd, ql, qr = draw.uniform(3.5, 5.5), draw.uniform(2.0, 5.0), draw.uniform(0.6, 1.5)
    qw = draw.uniform(0.6, 1.4)
    k_base = draw.uniform(15000, 80000)
    gravity_bay = draw.integer(0, max(nb - 1, 0))
    bracing = draw.choice(BRACING_PATTERNS)

    def xyv(s, v):
        return (s, 0.0, v) if plane == "xv" else (0.0, s, v)

    along, rot = ("x", "ry") if plane == "xv" else ("y", "rx")
    roll = 0.0 if plane == "xv" else 90.0
    for i, x in enumerate(xs):
        for k, v in enumerate(vs):
            B.node(f"{line_id(i)}{k}", *xyv(x, v))
    split = ns >= 3
    lower = B.group("columnLower" if split else "column", "column", 5)
    upper = B.group("columnUpper", "column", 4) if split else lower
    half = (ns + 1) // 2
    gb, gr = B.group("beam", "beam", 8), B.group("roofBeam", "beam", 7)
    columns, beams = {lower: [], upper: []}, {gb: [], gr: []}
    for i in range(nb + 1):
        for k in range(ns):
            g = lower if k < half else upper
            columns[g].append(B.member(f"C{line_id(i)}{k + 1}", f"{line_id(i)}{k}", f"{line_id(i)}{k + 1}", g, roll=roll))
    d = (1.0, 0.0, 0.0) if plane == "xv" else (0.0, 1.0, 0.0)
    for k in range(1, ns + 1):
        g = gr if k == ns else gb
        for i in range(nb):
            pinned = pinned_bay and i == gravity_bay
            mid = B.member(f"B{line_id(i)}{k}", f"{line_id(i)}{k}", f"{line_id(i + 1)}{k}", g,
                           releases={"i": {"rz": True}, "j": {"rz": True}} if pinned else None)
            beams[g].append(mid)
            if offsets:
                col = lower if k - 1 < half else upper
                B.offset_half_depth(mid, "i", col, d)
                B.offset_half_depth(mid, "j", col, tuple(-c for c in d))
    for i in range(nb + 1):
        nid = f"{line_id(i)}0"
        B.support(nid, "pinned")
        if springs:
            B.support(nid, **{rot: {"spring": B.k_rot(k_base)}})
        elif base == "fixed":
            B.support(nid, **{rot: "fixed"})
    B.restrain_out_of_plane(plane)

    B.case("D", "dead load and self weight")
    B.self_weight("D")
    B.case("L1", "live load, bays 1, 3, ...")
    if nb > 1:
        B.case("L2", "live load, bays 2, 4, ...")
    for k in range(1, ns + 1):
        roof = k == ns
        for i in range(nb):
            mid = f"B{line_id(i)}{k}"
            B.line("D", mid, "v", -B.W((qd if not roof else 0.7 * qd) * spacing))
            B.line("L1" if i % 2 == 0 else "L2", mid, "v", -B.W((ql if not roof else qr) * spacing))
    B.case("W", "wind")
    for k in range(1, ns + 1):
        trib = ((vs[k] - vs[k - 1]) + (vs[k + 1] - vs[k] if k < ns else 0.0)) / 2.0 / B.us["per_m"]
        B.nodal("W", f"A{k}", **{"F" + along: B.F(qw * spacing * trib * (0.8 + 0.4 * k / ns))})
    live = [c for c in ("L1", "L2") if c in B.cases]
    full16, full1, half5 = {c: 1.6 for c in live}, {c: 1.0 for c in live}, {c: 0.5 for c in live}
    cands = [("ULS1", {"D": 1.2, **full16}, "strength"),
             ("SLS1", {"D": 1.0, **full1}, "serviceability"),
             ("ULS2", {"D": 1.2, **half5, "W": 1.0}, "strength"),
             ("ULS3", {"D": 0.9, "W": 1.0}, "strength"),
             ("SLS2", {"D": 1.0, "W": 0.6}, "serviceability"),
             ("ULS4", {"D": 1.2, "L1": 1.6}, "strength"),
             ("ULS5", {"D": 1.2, **half5, "W": -1.0}, "strength"),
             ("ULS6", {"D": 1.4}, "strength")]
    add_combinations(B, cands, n_combos, second_order=second)
    buckling = B.request_buckling(flag(features, "buckling", draw.chance(BUCKLING_CHANCE[level])), modes=1)
    sls = service_ids(cands, n_combos)
    kx = 1.0 if second else 1.5
    for g, ids in columns.items():
        if ids:
            B.design_group(g + "s", ids, role="column", Kx=kx)
    for g, ids in beams.items():
        B.design_group(g + "s", ids, role="beam", bracing=dict(bracing[1]),
                       deflectionLimits=[{"spanRatio": 300, "direction": "y", "combinations": [sls[0]]}])
    B.parameters = {
        "buckling": buckling,
        "plane": plane, "bays": nb, "storeys": ns, "gridLines": xs, "levels": vs,
        "bases": "spring" if springs else base, "baseSpring": B.k_rot(k_base) if springs else None,
        "frameSpacing_m": round(spacing, 3), "secondOrder": second,
        "pinnedBay": gravity_bay + 1 if pinned_bay else None, "rigidOffsets": offsets,
        "nCombinations": n_combos, "bracing": bracing[0], "columnKx": kx,
    }
    return B
