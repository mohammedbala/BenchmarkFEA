"""Three-dimensional building frame on a rectangular grid: shared by ``multistorey_3d`` and
``diaphragm_frame``.

Grid lines A, B, ... along x and 1, 2, ... along y; node ``B2.1`` is on line B / 2 at level 1
(level 0 = base).  Columns ``C<line><n>.<storey>``, beams along x ``BX<line><n>.<level>`` (from
the node named to the next line along x), beams along y ``BY<line><n>.<level>``.

Floors span one way, along y, onto the beams along x; the beams along y on the two end lines carry
a cladding line load.  Columns on lines A, C, ... have their web along x (strong axis resists sway
along x); those on lines B, D, ... are turned by 90 degrees.
"""
from __future__ import annotations

from ._common import (BRACING_PATTERNS, Builder, GenerationRejected, add_combinations, grid_lines, line_id,
                      service_ids, storey_levels)


def nid(i: int, j: int, k: int) -> str:
    return f"{line_id(i)}{j + 1}.{k}"


def _check_master_lines(cx: float, cy: float, points) -> None:
    """Reject a master position from which the straight line to a slave passes through another slave."""
    for ax, ay in points:
        dx, dy = ax - cx, ay - cy
        length2 = dx * dx + dy * dy
        if length2 == 0.0:
            raise GenerationRejected("diaphragm master on a grid node")
        for bx, by in points:
            if (bx, by) == (ax, ay):
                continue
            t = ((bx - cx) * dx + (by - cy) * dy) / length2
            off = abs((bx - cx) * dy - (by - cy) * dx) / length2 ** 0.5
            if 0.0 < t < 1.0 and off < 1e-6 * length2 ** 0.5:
                raise GenerationRejected("a diaphragm slave lies on the line from the master to another slave")


def space_frame(B: Builder, draw, shape, *, second: bool, pinned_secondary: bool, springs: bool,
                diaphragm: bool, braced_side: bool, torsion: bool, n_combos: int, buckling: bool) -> dict:
    """Build the frame into ``B``; returns the parameters drawn.  The draws are always consumed in
    the same order, whatever the options."""
    nx, ny, ns = shape
    xs = grid_lines(B, draw, nx, 5.0, 8.0, 0.5, equal=draw.chance(0.6))
    ys = grid_lines(B, draw, ny, 4.0, 7.0, 0.5, equal=draw.chance(0.6))
    vs = storey_levels(B, draw, ns, first=(3.5, 5.0), other=(3.0, 4.0))
    base = draw.choice(("fixed", "fixed", "pinned"))
    qd, ql, qr = draw.uniform(3.5, 5.5), draw.uniform(2.0, 5.0), draw.uniform(0.6, 1.5)
    qc, qw = draw.uniform(2.0, 5.0), draw.uniform(0.8, 1.6)
    k_base = draw.uniform(15000, 80000)
    ecc = draw.choice((0.05, 0.1, 0.15)), draw.choice((0.05, 0.1, 0.15))
    if ecc[0] == ecc[1]:
        # equal relative eccentricities put the master on a diagonal of the grid: the line from the
        # master to a corner node would pass through other grid nodes (PyNite, which emulates the
        # diaphragm with members from the master to each slave, would connect them)
        ecc = ecc[0], {0.05: 0.1, 0.1: 0.15, 0.15: 0.05}[ecc[1]]
    brace_kind = draw.choice(("pipe", "box"))
    bracing = draw.choice(BRACING_PATTERNS)
    per_m = B.us["per_m"]
    if pinned_secondary or braced_side:
        base = "fixed" if not springs else base

    for k, v in enumerate(vs):
        for i, x in enumerate(xs):
            for j, y in enumerate(ys):
                B.node(nid(i, j, k), x, y, v)
    split = ns >= 3
    half = (ns + 1) // 2
    lower = B.group("columnLower" if split else "column", "column", 5)
    upper = B.group("columnUpper", "column", 4) if split else lower
    gx, gxr = B.group("beamX", "beam", 8), B.group("roofBeamX", "beam", 7)
    gy = B.group("beamY", "beam", 6)
    columns, beams_x, beams_y, braces = {lower: [], upper: []}, {gx: [], gxr: []}, [], []
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(ns):
                g = lower if k < half else upper
                columns[g].append(B.member(f"C{line_id(i)}{j + 1}.{k + 1}", nid(i, j, k), nid(i, j, k + 1), g,
                                           roll=90.0 if i % 2 else 0.0))
    pin = {"rz": True}
    for k in range(1, ns + 1):
        for j in range(ny + 1):
            for i in range(nx):
                g = gxr if k == ns else gx
                beams_x[g].append(B.member(f"BX{line_id(i)}{j + 1}.{k}", nid(i, j, k), nid(i + 1, j, k), g))
        for i in range(nx + 1):
            for j in range(ny):
                interior = 0 < i < nx
                rel = {"i": pin, "j": pin} if (pinned_secondary and interior) else None
                beams_y.append(B.member(f"BY{line_id(i)}{j + 1}.{k}", nid(i, j, k), nid(i, j + 1, k), gy, releases=rel))
    if braced_side:
        gb = B.group("brace", brace_kind, 8)
        for k in range(1, ns + 1):                                   # one diagonal per storey in the face y = 0
            a, b = (0, 1) if k % 2 else (1, 0)
            braces.append(B.member(f"X{k}", nid(a, 0, k - 1), nid(b, 0, k), gb, type="truss"))
    for i in range(nx + 1):
        for j in range(ny + 1):
            n = nid(i, j, 0)
            B.support(n, "pinned")
            if springs:
                k_file = B.k_rot(k_base)
                B.support(n, rx={"spring": k_file}, ry={"spring": k_file}, rv="fixed")
            elif base == "fixed":
                B.support(n, "fixed")

    # ---- gravity
    B.case("D", "dead load and self weight")
    B.self_weight("D")
    B.case("L", "live load")
    for k in range(1, ns + 1):
        roof = k == ns
        for j in range(ny + 1):
            trib = ((ys[j] - ys[j - 1] if j > 0 else 0.0) + (ys[j + 1] - ys[j] if j < ny else 0.0)) / 2.0 / per_m
            for i in range(nx):
                mid = f"BX{line_id(i)}{j + 1}.{k}"
                B.line("D", mid, "v", -B.W((qd if not roof else 0.7 * qd) * trib))
                B.line("L", mid, "v", -B.W((ql if not roof else qr) * trib))
        if not roof:
            for i in (0, nx):
                for j in range(ny):
                    B.line("D", f"BY{line_id(i)}{j + 1}.{k}", "v", -B.W(qc))
    # ---- wind: storey forces
    B.case("WX", "wind along x")
    B.case("WY", "wind along y")
    width_x, width_y = (xs[-1] - xs[0]) / per_m, (ys[-1] - ys[0]) / per_m
    masters = []
    for k in range(1, ns + 1):
        trib = ((vs[k] - vs[k - 1]) + (vs[k + 1] - vs[k] if k < ns else 0.0)) / 2.0 / per_m
        shape_k = 0.8 + 0.4 * k / ns
        fx, fy = qw * 1.3 * width_y * trib * shape_k, qw * 1.3 * width_x * trib * shape_k
        if diaphragm:
            # master at the (eccentric) mass centre, never on a grid line
            cx = (xs[0] + xs[-1]) / 2.0 + ecc[0] * (xs[-1] - xs[0])
            cy = (ys[0] + ys[-1]) / 2.0 + ecc[1] * (ys[-1] - ys[0])
            m = B.node(f"CM{k}", cx, cy, vs[k])
            masters.append(m)
            _check_master_lines(cx, cy, [(x, y) for x in xs for y in ys])
            B.constraint(f"floor{k}", "rigidDiaphragm", m,
                         [nid(i, j, k) for i in range(nx + 1) for j in range(ny + 1)])
            # accidental torsion: the storey force times 5 % of the plan dimension normal to it
            B.nodal("WX", m, Fx=B.F(fx), Mv=B.M(0.05 * width_y * fx) if torsion else 0.0)
            B.nodal("WY", m, Fy=B.F(fy), Mv=-B.M(0.05 * width_x * fy) if torsion else 0.0)
        else:
            for j in range(ny + 1):
                share = (1.0 if 0 < j < ny else 0.5) / ny
                B.nodal("WX", nid(0, j, k), Fx=B.F(fx * share))
            for i in range(nx + 1):
                share = (1.0 if 0 < i < nx else 0.5) / nx
                B.nodal("WY", nid(i, 0, k), Fy=B.F(fy * share))
    cands = [("ULS1", {"D": 1.2, "L": 1.6}, "strength"),
             ("SLS1", {"D": 1.0, "L": 1.0}, "serviceability"),
             ("ULS2", {"D": 1.2, "L": 0.5, "WX": 1.0}, "strength"),
             ("ULS3", {"D": 1.2, "L": 0.5, "WY": 1.0}, "strength"),
             ("SLS2", {"D": 1.0, "WX": 0.6}, "serviceability"),
             ("ULS4", {"D": 0.9, "WX": 1.0}, "strength"),
             ("ULS5", {"D": 0.9, "WY": 1.0}, "strength"),
             ("ULS6", {"D": 1.4}, "strength")]
    ids = add_combinations(B, cands, n_combos, second_order=second)
    used = {c for cid in ids for c in B.combinations[cid]["factors"]}
    for c in [c for c in B.cases if c not in used]:
        del B.cases[c]
    buckling = B.request_buckling(buckling)
    sls = service_ids(cands, n_combos)
    k_col = 1.0 if second else 1.5
    for g, members in columns.items():
        B.design_group(g + "s", members, role="column", Kx=k_col, Ky=k_col)
    for g, members in beams_x.items():
        if members:
            B.design_group(g + "s", members, role="beam", bracing=dict(bracing[1]),
                           deflectionLimits=[{"spanRatio": 300, "direction": "y", "combinations": [sls[0]]}])
    B.design_group("beamsY", beams_y, role="beam", bracing={"topFlange": [0.5], "bottomFlange": [0.5]},
                   deflectionLimits=[{"spanRatio": 300, "direction": "y", "combinations": [sls[0]]}])
    if braces:
        B.design_group("braces", braces, role="brace")
    return {
        "buckling": buckling,
        "baysX": nx, "baysY": ny, "storeys": ns, "gridX": xs, "gridY": ys, "levels": vs,
        "bases": "spring" if springs else base, "baseSpring": B.k_rot(k_base) if springs else None,
        "secondOrder": second, "pinnedSecondaryBeams": bool(pinned_secondary and nx >= 2),
        "rigidDiaphragm": diaphragm, "diaphragmMasters": masters,
        "massCentreEccentricity": list(ecc) if diaphragm else None, "accidentalTorsion": bool(diaphragm and torsion),
        "bracedSide": bool(braced_side), "braceSection": brace_kind if braced_side else None,
        "nCombinations": n_combos, "beamBracing": bracing[0], "columnK": k_col,
    }
