"""Single-bay portal frames, flat or pitched (plane model in 3D space, out-of-plane DOFs restrained).

Two columns and a beam or two rafters; pinned, fixed or partially fixed bases; dead load, roof
live load (projected on pitched rafters) and wind on the columns; optionally second-order
analysis, semi-rigid eaves connections (rotational springs) and rigid end offsets of the beam
(half the column depth; only in models without second-order combinations, see families/__init__).
"""
from __future__ import annotations

from ._common import BUCKLING_CHANCE, BUCKLING_FLAG, BRACING_PATTERNS, Builder, add_combinations, flag, service_ids

NAME = "portal_frame"
LEVELS = (1, 2, 3)
FEATURES = {
    "buckling": BUCKLING_FLAG,
    "pitched": "pitched roof with a ridge node (two rafters)",
    "secondOrder": "strength combinations analysed second order",
    "baseSprings": "partially fixed bases (rotational springs)",
    "semiRigid": "semi-rigid eaves connections (rotational springs at the beam / rafter ends)",
    "rigidOffsets": "rigid end offsets of the beam at the columns (only without second order)",
    "craneLoad": "point load with a moment on one column (bracket load)",
}


def build(draw, level: int, unit_system: str, vertical: str, features: dict) -> Builder:
    B = Builder(unit_system, vertical, NAME)
    plane = draw.choice(("xv", "yv"))
    pitched = flag(features, "pitched", draw.chance(0.5))
    span = draw.step(8.0, 24.0, 1.0) if pitched else draw.step(6.0, 12.0, 0.5)
    height = draw.step(4.0, 8.0, 0.5)
    rise = draw.step(0.5, 2.5, 0.25)
    base = draw.choice(("pinned", "fixed"))
    second = flag(features, "secondOrder", draw.chance({1: 0.0, 2: 0.6, 3: 0.7}[level]))
    springs = flag(features, "baseSprings", draw.chance({1: 0.0, 2: 0.3, 3: 0.8}[level]))
    semi = flag(features, "semiRigid", draw.chance({1: 0.0, 2: 0.2, 3: 0.7}[level]))
    offsets = flag(features, "rigidOffsets", draw.chance({1: 0.0, 2: 0.4, 3: 0.9}[level])) and not second
    crane = flag(features, "craneLoad", draw.chance({1: 0.0, 2: 0.4, 3: 0.6}[level]))
    n_combos = {1: 2, 2: draw.integer(3, 5), 3: draw.integer(6, 8)}[level]
    spacing = draw.uniform(5.0, 8.0)
    qd = draw.uniform(0.6, 1.6) if pitched else draw.uniform(3.0, 5.0)        # kPa
    ql = draw.uniform(0.6, 2.0) if pitched else draw.uniform(2.0, 5.0)
    qw = draw.uniform(0.5, 1.2)
    k_base, k_eave = draw.uniform(10000, 60000), draw.uniform(20000, 150000)
    crane_p, crane_e, crane_at = draw.uniform(40, 150), draw.choice((0.3, 0.4, 0.5)), draw.choice((0.6, 0.7, 0.75))
    bracing = draw.choice(BRACING_PATTERNS)

    def xyv(s, v):
        return (s, 0.0, v) if plane == "xv" else (0.0, s, v)

    along, rot = ("x", "ry") if plane == "xv" else ("y", "rx")
    roll = 0.0 if plane == "xv" else 90.0
    Ls, H = B.L(span, 0.5), B.L(height, 0.5)
    B.node("N1", *xyv(0.0, 0.0))
    B.node("N2", *xyv(0.0, H))
    if pitched:
        B.node("N3", *xyv(Ls / 2.0, H + B.L(rise, 0.25)))
    B.node("N4", *xyv(Ls, H))
    B.node("N5", *xyv(Ls, 0.0))
    col, raf = B.group("column", "column", 5), B.group("rafter", "beam", 8)
    B.member("C1", "N1", "N2", col, roll=roll)
    B.member("C2", "N5", "N4", col, roll=roll)
    k = B.k_rot(k_eave)
    rel = {"rz": {"spring": k}}
    if pitched:
        B.member("R1", "N2", "N3", raf, releases={"i": rel} if semi else None)
        B.member("R2", "N3", "N4", raf, releases={"j": rel} if semi else None)
        beams = ["R1", "R2"]
    else:
        B.member("B1", "N2", "N4", raf, releases={"i": rel, "j": rel} if semi else None)
        beams = ["B1"]
    if offsets:
        d = (1.0, 0.0, 0.0) if plane == "xv" else (0.0, 1.0, 0.0)
        B.offset_half_depth(beams[0], "i", col, d)
        B.offset_half_depth(beams[-1], "j", col, tuple(-c for c in d))
    for nid in ("N1", "N5"):
        B.support(nid, "pinned")
        if springs:
            B.support(nid, **{rot: {"spring": B.k_rot(k_base)}})
        elif base == "fixed":
            B.support(nid, **{rot: "fixed"})
    B.restrain_out_of_plane(plane)

    B.case("D", "dead load and self weight")
    B.self_weight("D")
    B.case("L", "roof live load" if pitched else "floor live load")
    for mid in beams:
        B.line("D", mid, "v", -B.W(qd * spacing))
        B.line("L", mid, "v", -B.W(ql * spacing), projected=pitched)
    lateral = n_combos > 2
    if lateral:
        B.case("W", "wind")
        B.line("W", "C1", along, B.W(qw * spacing))
        B.line("W", "C2", along, B.W(0.5 * qw * spacing))
    if crane and lateral:
        B.case("C", "bracket load on column C1")
        B.point("C", "C1", "v", -B.F(crane_p), crane_at)
        # moment of the eccentric bracket load about the out-of-plane axis (vector along the plane normal)
        moment = B.M(crane_p * crane_e)
        normal = "y" if plane == "xv" else "x"
        B.point("C", "C1", normal, moment if plane == "xv" else -moment, crane_at, kind="moment")
    else:
        crane = False
    c = {"C": 1.0} if crane else {}
    c16 = {"C": 1.6} if crane else {}
    cands = [("ULS1", {"D": 1.2, "L": 1.6, **c16}, "strength"),
             ("SLS1", {"D": 1.0, "L": 1.0, **c}, "serviceability"),
             ("ULS2", {"D": 1.2, "L": 0.5, "W": 1.0, **c}, "strength"),
             ("ULS3", {"D": 0.9, "W": 1.0}, "strength"),
             ("ULS4", {"D": 1.4}, "strength"),
             ("SLS2", {"D": 1.0, "W": 0.6}, "serviceability"),
             ("ULS5", {"D": 1.2, "L": 0.5, "W": -1.0, **c}, "strength"),
             ("ULS6", {"D": 0.9, "W": -1.0}, "strength")]
    add_combinations(B, cands, n_combos, second_order=second)
    buckling = B.request_buckling(flag(features, "buckling", draw.chance(BUCKLING_CHANCE[level])), modes=2)
    sls = service_ids(cands, n_combos)
    B.design_group("columns", ["C1", "C2"], role="column", Kx=1.0 if second else (2.0 if base == "pinned" and not springs else 1.5),
                   bracing={"minor": [0.5], "torsion": [0.5]})
    B.design_group("rafters", beams, role="beam", bracing=dict(bracing[1]),
                   deflectionLimits=[{"spanRatio": 240 if pitched else 300, "direction": "y", "combinations": [sls[0]]}])
    B.parameters = {
        "buckling": buckling,
        "plane": plane, "span": Ls, "height": H, "rise": B.L(rise, 0.25) if pitched else None,
        "bases": "spring" if springs else base, "baseSpring": B.k_rot(k_base) if springs else None,
        "frameSpacing_m": round(spacing, 3), "deadLoad": B.W(qd * spacing), "liveLoad": B.W(ql * spacing),
        "windLoad": B.W(qw * spacing) if lateral else None, "secondOrder": second, "semiRigid": k if semi else None,
        "rigidOffsets": offsets, "craneLoad": crane, "nCombinations": n_combos, "bracing": bracing[0],
    }
    return B
