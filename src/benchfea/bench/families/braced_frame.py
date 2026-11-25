"""Braced multi-storey frames (plane model in 3D space, out-of-plane DOFs restrained).

Continuous columns, beams pinned to the columns, one braced bay with X, chevron (inverted V) or
single-diagonal bracing.  Braces are pin-ended axial members (``truss``).  On request
(``features={"tensionOnly": True}``) the bracing is X-bracing with tension-only diagonals -- a
feature no foreign solver reproduces, so such a model cannot be cross-checked (recorded in the
metadata, ``tierC_expected``) and is never drawn by default.  Every combination contains a lateral load -- wind, or the minimum
lateral load case ``N`` in the gravity combinations -- so that one diagonal of every tension-only
panel is always stretched (with no lateral load at all both diagonals would be slack and the pinned
frame a mechanism).
"""
from __future__ import annotations

from ._common import (BUCKLING_CHANCE, BUCKLING_FLAG, TENSION_ONLY_CHANCE, BRACING_PATTERNS, Builder, add_combinations, flag, grid_lines, line_id, service_ids,
                      storey_levels)

NAME = "braced_frame"
LEVELS = (2, 3, 4)
FEATURES = {
    "buckling": BUCKLING_FLAG,
    "tensionOnly": "X-bracing with tension-only diagonals (opt-in: not reproduced by the foreign solvers)",
    "secondOrder": "strength combinations analysed second order",
    "baseSprings": "partially fixed bases (rotational springs)",
}
_SHAPES = {2: ((1, 1), (1, 2), (2, 1)), 3: ((2, 2), (2, 3), (3, 2), (1, 3)), 4: ((3, 3), (3, 4), (2, 5), (4, 3))}
PATTERNS = ("X", "chevron", "diagonal")


def build(draw, level: int, unit_system: str, vertical: str, features: dict) -> Builder:
    B = Builder(unit_system, vertical, NAME)
    plane = draw.choice(("xv", "yv"))
    nb, ns = draw.choice(_SHAPES[level])
    xs = grid_lines(B, draw, nb, 5.0, 8.0, 0.5, equal=draw.chance(0.7))
    vs = storey_levels(B, draw, ns, first=(3.5, 5.0), other=(3.0, 4.0))
    pattern = draw.choice(PATTERNS)
    tension_only = flag(features, "tensionOnly", draw.chance(TENSION_ONLY_CHANCE))
    second = flag(features, "secondOrder", draw.chance({2: 0.0, 3: 0.6, 4: 0.9}[level]))
    springs = flag(features, "baseSprings", draw.chance({2: 0.0, 3: 0.3, 4: 0.5}[level]))
    n_combos = {2: draw.integer(3, 5), 3: draw.integer(4, 7), 4: draw.integer(6, 8)}[level]
    spacing = draw.uniform(5.0, 8.0)
    qd, ql, qr = draw.uniform(3.5, 5.5), draw.uniform(2.0, 5.0), draw.uniform(0.6, 1.5)
    qw = draw.uniform(0.8, 1.6)
    k_base = draw.uniform(15000, 80000)
    braced_bay = draw.integer(0, nb - 1)
    brace_kind = draw.choice(("pipe", "box"))
    bracing = draw.choice(BRACING_PATTERNS)
    if tension_only:
        # X pattern; first order only, with the stability check by a linear buckling analysis in the file:
        # the engine refuses a second-order analysis of such a frame under gravity with a small lateral
        # load ("beyond the critical load", although the buckling factor on the working diagonals is
        # large) -- see the report of this family in docs/DIFFICULTY.md, section 5
        pattern, second = "X", False

    def xyv(s, v):
        return (s, 0.0, v) if plane == "xv" else (0.0, s, v)

    along, rot = ("x", "ry") if plane == "xv" else ("y", "rx")
    roll = 0.0 if plane == "xv" else 90.0
    for i, x in enumerate(xs):
        for k, v in enumerate(vs):
            B.node(f"{line_id(i)}{k}", *xyv(x, v))
    split = ns >= 3
    half = (ns + 1) // 2
    col = B.group("columnLower" if split else "column", "column", 4)
    col_up = B.group("columnUpper", "column", 3) if split else col
    gb, gr = B.group("beam", "beam", 8), B.group("roofBeam", "beam", 7)
    gx = B.group("braceLower" if split else "brace", brace_kind, 6)
    gx_up = B.group("braceUpper", brace_kind, 5) if split else gx
    columns, beams, braces, chains = {col: [], col_up: []}, {gb: [], gr: []}, {gx: [], gx_up: []}, {}
    for i in range(nb + 1):
        for k in range(ns):
            g = col if k < half else col_up
            columns[g].append(B.member(f"C{line_id(i)}{k + 1}", f"{line_id(i)}{k}", f"{line_id(i)}{k + 1}", g, roll=roll))
    pin = {"rz": True}
    a, b = line_id(braced_bay), line_id(braced_bay + 1)
    for k in range(1, ns + 1):
        g = gr if k == ns else gb
        for i in range(nb):
            mid, ni, nj = f"B{line_id(i)}{k}", f"{line_id(i)}{k}", f"{line_id(i + 1)}{k}"
            rigid = False
            if pattern == "chevron" and i == braced_bay:
                m = B.node(f"M{k}", *xyv((xs[i] + xs[i + 1]) / 2.0, vs[k]))
                B.member(mid + "a", ni, m, g, releases={"i": pin})
                B.member(mid + "b", m, nj, g, releases={"j": pin})
                chains[mid] = [mid + "a", mid + "b"]
            else:
                B.member(mid, ni, nj, g, releases=None if rigid else {"i": pin, "j": pin})
            beams[g].append(mid)
    kind = "tensionOnly" if tension_only else "truss"
    for k in range(1, ns + 1):
        g = gx if k - 1 < half else gx_up
        if pattern == "X":
            braces[g].append(B.member(f"X{k}a", f"{a}{k - 1}", f"{b}{k}", g, type=kind))
            braces[g].append(B.member(f"X{k}b", f"{b}{k - 1}", f"{a}{k}", g, type=kind))
        elif pattern == "chevron":
            braces[g].append(B.member(f"X{k}a", f"{a}{k - 1}", f"M{k}", g, type=kind))
            braces[g].append(B.member(f"X{k}b", f"{b}{k - 1}", f"M{k}", g, type=kind))
        else:
            lo, hi = (a, b) if k % 2 else (b, a)                # zig-zag up the bay
            braces[g].append(B.member(f"X{k}", f"{lo}{k - 1}", f"{hi}{k}", g, type=kind))
    for dm, chain in chains.items():
        B.design_member(dm, chain)
    base = "pinned"
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
    B.case("L", "live load")
    gravity = [0.0] * (ns + 1)                                   # unfactored D + L per level, kN
    for k in range(1, ns + 1):
        roof = k == ns
        wd, wl = (qd if not roof else 0.7 * qd) * spacing, (ql if not roof else qr) * spacing
        for i in range(nb):
            for mid in chains.get(f"B{line_id(i)}{k}", [f"B{line_id(i)}{k}"]):
                B.line("D", mid, "v", -B.W(wd))
                B.line("L", mid, "v", -B.W(wl))
        gravity[k] = (wd + wl) * (xs[-1] - xs[0]) / B.us["per_m"]
    B.case("W", "wind")
    B.case("N", "minimum lateral load (0.2 % of the dead and live load of each level)")
    for k in range(1, ns + 1):
        trib = ((vs[k] - vs[k - 1]) + (vs[k + 1] - vs[k] if k < ns else 0.0)) / 2.0 / B.us["per_m"]
        B.nodal("W", f"A{k}", **{"F" + along: B.F(qw * spacing * trib * (0.8 + 0.4 * k / ns) * 1.5)})
        B.nodal("N", f"A{k}", **{"F" + along: max(B.F(0.002 * gravity[k]), B.F(0.5) or 0.1)})
    cands = [("ULS1", {"D": 1.2, "L": 1.6, "N": 1.6}, "strength"),
             ("SLS1", {"D": 1.0, "L": 1.0, "N": 1.0}, "serviceability"),
             ("ULS2", {"D": 1.2, "L": 0.5, "W": 1.0}, "strength"),
             ("ULS3", {"D": 0.9, "W": 1.0}, "strength"),
             ("SLS2", {"D": 1.0, "W": 0.6}, "serviceability"),
             ("ULS4", {"D": 1.2, "L": 0.5, "W": -1.0}, "strength"),
             ("ULS5", {"D": 0.9, "W": -1.0}, "strength"),
             ("ULS6", {"D": 1.4, "N": 1.4}, "strength")]
    add_combinations(B, cands, n_combos, second_order=second)
    buckling = B.request_buckling(flag(features, "buckling", draw.chance(BUCKLING_CHANCE[level])) or tension_only,
                                  modes=1)
    sls = service_ids(cands, n_combos)
    for g, ids in columns.items():
        B.design_group(g + "s", ids, role="column", Kx=1.0)
    for g, ids in beams.items():
        if not ids:
            continue
        B.design_group(g + "s", ids, role="beam", bracing=dict(bracing[1]),
                       deflectionLimits=[{"spanRatio": 300, "direction": "y", "combinations": [sls[0]]}])
    for g, ids in braces.items():
        B.design_group(g + "s", ids, role="brace")
    B.parameters = {
        "buckling": buckling,
        "plane": plane, "bays": nb, "storeys": ns, "gridLines": xs, "levels": vs, "bracing": pattern,
        "bracedBay": braced_bay + 1, "braceType": kind, "braceSection": brace_kind,
        "bases": "spring" if springs else base, "baseSpring": B.k_rot(k_base) if springs else None,
        "frameSpacing_m": round(spacing, 3), "secondOrder": second, "nCombinations": n_combos,
        "beamBracing": bracing[0],
    }
    return B
