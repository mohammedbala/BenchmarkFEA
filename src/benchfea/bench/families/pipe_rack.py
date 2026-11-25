"""Pipe racks (true 3D): transverse moment bents joined by longitudinal struts, braced lengthwise.

Bents (two columns and one or two tier beams, rigid in the plane of the bent) stand in planes normal to x; pinned
struts along x join them at every tier on both column lines, and one bay of each column line has
vertical bracing (X or single diagonal; on request -- ``features={"tensionOnly": True}``, never by
default, because no foreign solver reproduces it -- tension-only X-bracing; then every combination
contains the friction case ``F``, so that one diagonal of each panel is always stretched).  Loads:
pipes as point loads on the tier beams (empty weight ``D``, contents ``O``), friction of the pipes
along x (``F``, horizontal point loads on the tier beams: weak-axis bending) and transverse wind on
the columns and at the tiers (``W``).
"""
from __future__ import annotations

from ._common import BUCKLING_CHANCE, BUCKLING_FLAG, TENSION_ONLY_CHANCE, Builder, add_combinations, flag, line_id, service_ids, tidy

NAME = "pipe_rack"
LEVELS = (3, 4, 5)
FEATURES = {
    "buckling": BUCKLING_FLAG,
    "tensionOnly": "longitudinal X-bracing with tension-only diagonals (opt-in: not reproduced by the foreign solvers)",
    "secondOrder": "strength combinations analysed second order",
    "baseSprings": "partially fixed bases (rotational springs about both horizontal axes)",
}
_SHAPES = {3: ((2, 1), (3, 1)), 4: ((3, 1), (2, 2), (3, 2), (4, 1)), 5: ((3, 2), (4, 2), (4, 1))}   # (bents, tiers)


def build(draw, level: int, unit_system: str, vertical: str, features: dict) -> Builder:
    B = Builder(unit_system, vertical, NAME)
    nb, nt = draw.choice(_SHAPES[level])
    spacing = B.L(draw.step(5.0, 8.0, 0.5), 0.5)
    width = B.L(draw.step(4.0, 8.0, 0.5), 0.5)
    h1, h2 = B.L(draw.step(4.0, 6.0, 0.5), 0.5), B.L(draw.step(2.0, 3.0, 0.5), 0.5)
    pattern = draw.choice(("X", "diagonal"))
    tension_only = flag(features, "tensionOnly", draw.chance(TENSION_ONLY_CHANCE))
    second = flag(features, "secondOrder", draw.chance({3: 0.2, 4: 0.7, 5: 0.9}[level]))
    springs = flag(features, "baseSprings", draw.chance({3: 0.1, 4: 0.3, 5: 0.7}[level]))
    n_combos = {3: draw.integer(3, 5), 4: draw.integer(4, 7), 5: draw.integer(6, 8)}[level]
    n_pipes = draw.integer(2, 4)
    p_empty, p_content = draw.uniform(4.0, 12.0), draw.uniform(4.0, 14.0)        # kN per pipe per bent
    qw, friction = draw.uniform(0.8, 1.6), draw.choice((0.1, 0.2, 0.3))
    k_base = draw.uniform(15000, 80000)
    braced_bay = draw.integer(0, nb - 2)
    brace_kind = draw.choice(("pipe", "box"))
    if tension_only:                                         # as in braced_frame: X, first order, buckling in the file
        pattern, second = "X", False
    levels = [0.0, h1] + ([tidy(h1 + h2)] if nt == 2 else [])
    per_m = B.us["per_m"]

    def nid(b, side, k):
        return f"{line_id(b)}{side}{k}"

    for b in range(nb):
        for side, y in (("L", 0.0), ("R", width)):
            for k, v in enumerate(levels):
                B.node(nid(b, side, k), tidy(b * spacing), y, v)
    col, beam = B.group("column", "column", 4), B.group("tierBeam", "beam", 6)
    strut, brace = B.group("strut", "box", 5), B.group("brace", brace_kind, 6)
    columns, beams, struts, braces = [], [], [], []
    for b in range(nb):
        for side in "LR":
            for k in range(nt):
                # web in the plane of the bent: strong axis resists transverse sway
                columns.append(B.member(f"C{line_id(b)}{side}{k + 1}", nid(b, side, k), nid(b, side, k + 1), col, roll=90.0))
        for k in range(1, nt + 1):
            # rigid in the plane of the bent; simply supported for the horizontal (friction) loads, so
            # that no weak-axis end moment of the beam twists the columns
            beams.append(B.member(f"T{line_id(b)}{k}", nid(b, "L", k), nid(b, "R", k), beam,
                                  releases={"i": {"ry": True}, "j": {"ry": True}}))
    pin = {"ry": True, "rz": True}
    for b in range(nb - 1):
        for side in "LR":
            for k in range(1, nt + 1):
                # pinned both ways at both ends and free to twist at end j: a strut carries no torsion
                struts.append(B.member(f"S{line_id(b)}{side}{k}", nid(b, side, k), nid(b + 1, side, k), strut,
                                       releases={"i": pin, "j": {"rx": True, **pin}}))
    kind = "tensionOnly" if tension_only else "truss"
    a, c = braced_bay, braced_bay + 1
    for side in "LR":
        for k in range(1, nt + 1):
            if pattern == "X":
                braces.append(B.member(f"X{side}{k}a", nid(a, side, k - 1), nid(c, side, k), brace, type=kind))
                braces.append(B.member(f"X{side}{k}b", nid(c, side, k - 1), nid(a, side, k), brace, type=kind))
            else:
                lo, hi = (a, c) if k % 2 else (c, a)
                braces.append(B.member(f"X{side}{k}", nid(lo, side, k - 1), nid(hi, side, k), brace, type=kind))
    base = draw.choice(("pinned", "fixed"))
    for b in range(nb):
        for side in "LR":
            n = nid(b, side, 0)
            B.support(n, "pinned", rv="fixed")                # the base plate holds the twist of the column
            if springs:
                k_file = B.k_rot(k_base)
                B.support(n, rx={"spring": k_file}, ry={"spring": k_file}, rv="fixed")
            elif base == "fixed":
                B.support(n, "fixed")

    B.case("D", "self weight and empty pipes")
    B.self_weight("D")
    B.case("O", "pipe contents (operating)")
    B.case("F", "pipe friction along the rack")
    B.case("W", "transverse wind")
    positions = [tidy(round((i + 1) / (n_pipes + 1), 3)) for i in range(n_pipes)]
    for b in range(nb):
        end = 0.5 if b in (0, nb - 1) else 1.0                # end bents carry half a bay of pipe
        for k in range(1, nt + 1):
            mid = f"T{line_id(b)}{k}"
            for pos in positions:
                B.point("D", mid, "v", -B.F(p_empty * end), pos)
                B.point("O", mid, "v", -B.F(p_content * end), pos)
                B.point("F", mid, "x", max(B.F(friction * (p_empty + p_content) * end), 0.1), pos)
            trib = spacing / per_m * end
            B.nodal("W", nid(b, "L", k), Fy=B.F(qw * 1.2 * trib * n_pipes * 0.3))
        for k in range(nt):
            B.line("W", f"C{line_id(b)}L{k + 1}", "y", B.W(qw * 0.4))
            B.line("W", f"C{line_id(b)}R{k + 1}", "y", B.W(qw * 0.4))
    cands = [("ULS1", {"D": 1.4, "O": 1.4, "F": 1.4}, "strength"),
             ("SLS1", {"D": 1.0, "O": 1.0}, "serviceability"),
             ("ULS2", {"D": 1.2, "O": 1.2, "W": 1.0}, "strength"),
             ("ULS3", {"D": 0.9, "W": 1.0}, "strength"),
             ("SLS2", {"D": 1.0, "O": 1.0, "W": 0.6}, "serviceability"),
             ("ULS4", {"D": 1.2, "O": 1.2, "F": -1.2}, "strength"),
             ("ULS5", {"D": 1.2, "O": 1.2, "F": 1.2, "W": 1.0}, "strength"),
             ("ULS6", {"D": 0.9, "W": -1.0}, "strength")]
    if tension_only:                                         # a lengthwise load in every combination
        cands = [(cid, factors if "F" in factors else {**factors, "F": 1.0}, purpose) for cid, factors, purpose in cands]
    ids = add_combinations(B, cands, n_combos, second_order=second)
    used = {c for cid in ids for c in B.combinations[cid]["factors"]}
    for c in [c for c in B.cases if c not in used]:
        del B.cases[c]
    buckling = B.request_buckling(flag(features, "buckling", draw.chance(BUCKLING_CHANCE[level])) or tension_only,
                                  modes=1)
    sls = service_ids(cands, n_combos)
    k_major = 1.0 if second else 1.5
    B.design_group("columns", columns, role="column", Kx=k_major, Ky=1.0 if base != "fixed" or second else 1.0)
    B.design_group("tierBeams", beams, role="beam", bracing={"topFlange": positions, "bottomFlange": [0.5]},
                   deflectionLimits=[{"spanRatio": 240, "direction": "y", "combinations": [sls[0]]}])
    B.design_group("struts", struts, role="beam")
    B.design_group("braces", braces, role="brace")
    B.parameters = {
        "buckling": buckling,
        "bents": nb, "tiers": nt, "bentSpacing": spacing, "width": width, "levels": levels,
        "bracing": pattern, "braceType": kind, "braceSection": brace_kind, "bracedBay": braced_bay + 1,
        "bases": "spring" if springs else base, "baseSpring": B.k_rot(k_base) if springs else None,
        "pipesPerTier": n_pipes, "pipePositions": positions, "friction": friction,
        "secondOrder": second, "nCombinations": n_combos, "columnKx": k_major,
    }
    return B
