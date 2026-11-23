"""Continuous beams on point supports (plane model in 3D space, out-of-plane DOFs restrained).

Spans of 4-10 m carrying uniform dead and live load (live load in two pattern cases: odd and even
spans), optionally a Gerber hinge inside a span, a rotational spring (partial fixity) at an end,
an elastic support, a settlement case, a cantilever overhang, a trapezoidal load and point loads.
"""
from __future__ import annotations

from ._common import BRACING_PATTERNS, Builder, GenerationRejected, add_combinations, flag, service_ids

NAME = "continuous_beam"
LEVELS = (1, 2, 3)
FEATURES = {
    "hinge": "Gerber hinge inside a span (full release of the bending moment)",
    "partialFixity": "rotational spring at an end support",
    "springSupport": "elastic vertical support at an interior support",
    "settlement": "load case with an imposed settlement of an interior support",
    "overhang": "cantilever overhang beyond the last support",
    "varyingLoad": "trapezoidal (linearly varying) live load on one span",
    "pointLoads": "point loads inside a span",
}


def build(draw, level: int, unit_system: str, vertical: str, features: dict) -> Builder:
    B = Builder(unit_system, vertical, NAME)
    plane = draw.choice(("xv", "yv"))
    n_spans = {1: draw.integer(1, 2), 2: draw.integer(2, 4), 3: draw.integer(6, 7)}[level]
    base_span = draw.step(4.0, 10.0, 0.5)
    spans = [max(3.0, base_span * draw.choice((0.8, 1.0, 1.0, 1.2))) for _ in range(n_spans)]
    ends = [draw.choice(("pinned", "fixed") if level < 3 else ("pinned", "fixed", "fixed")) for _ in range(2)]
    if n_spans == 1 and ends == ["pinned", "pinned"] and level > 1:
        ends[0] = "fixed"
    hinge = flag(features, "hinge", draw.chance({1: 0.0, 2: 0.6, 3: 0.8}[level]))
    partial = flag(features, "partialFixity", draw.chance({1: 0.0, 2: 0.4, 3: 0.5}[level]))
    spring = flag(features, "springSupport", draw.chance({1: 0.0, 2: 0.3, 3: 0.6}[level]))
    settle = flag(features, "settlement", draw.chance({1: 0.0, 2: 0.3, 3: 0.6}[level]))
    overhang = flag(features, "overhang", draw.chance({1: 0.0, 2: 0.2, 3: 0.4}[level]))
    varying = flag(features, "varyingLoad", draw.chance({1: 0.0, 2: 0.4, 3: 0.6}[level]))
    points = flag(features, "pointLoads", draw.chance({1: 0.3, 2: 0.6, 3: 0.8}[level]))
    n_combos = {1: 2, 2: draw.integer(3, 5), 3: draw.integer(6, 8)}[level]
    wd, wl = draw.uniform(8, 25), draw.uniform(6, 30)
    hinge_span, hinge_at = draw.integer(0, n_spans - 1), draw.choice((0.2, 0.25, 0.3))
    spring_at, settle_at = draw.integer(1, max(1, n_spans - 1)), draw.integer(1, max(1, n_spans - 1))
    k_spring, settlement = draw.uniform(15, 60), draw.choice((5, 10, 15, 20))
    k_end = draw.uniform(20000, 120000)
    load_span, load_at, load_p = draw.integer(0, n_spans - 1), draw.choice((0.25, 0.4, 0.5, 0.6)), draw.uniform(20, 90)
    vary_span = draw.integer(0, n_spans - 1)
    over_len, tip_load = draw.step(1.5, 3.0, 0.5), draw.uniform(5, 25)
    bracing = draw.choice(BRACING_PATTERNS)

    if n_spans < 2:
        spring = settle = False
    if hinge and n_spans == 1 and "fixed" not in ends:
        hinge = False
    if spring and settle and spring_at == settle_at:
        settle = False

    def xyv(s):                                              # position along the beam -> abstract coordinates
        return (s, 0.0, 0.0) if plane == "xv" else (0.0, s, 0.0)

    roll = 0.0
    lengths = [B.L(s, 0.5) for s in spans]
    x, stations = 0.0, [0.0]
    for length in lengths:
        x += length
        stations.append(x)
    g = B.group("beam", "beam", 4)
    members, chains = [], {}
    for k in range(n_spans + 1):
        B.node(f"N{k + 1}", *xyv(stations[k]))
    for k in range(n_spans):
        a, b = f"N{k + 1}", f"N{k + 2}"
        if hinge and k == hinge_span:
            h = B.node(f"H{k + 1}", *xyv(stations[k] + B.L(spans[k] * hinge_at, 0.25)))
            B.member(f"B{k + 1}a", a, h, g, roll=roll)
            B.member(f"B{k + 1}b", h, b, g, roll=roll, releases={"i": {"rz": True}})
            chains[f"B{k + 1}"] = [f"B{k + 1}a", f"B{k + 1}b"]
            members.append(f"B{k + 1}")
        else:
            B.member(f"B{k + 1}", a, b, g, roll=roll)
            members.append(f"B{k + 1}")
    for dm, chain in chains.items():
        B.design_member(dm, chain)
    # supports
    for k in range(n_spans + 1):
        nid = f"N{k + 1}"
        B.support(nid, uv="fixed")
    B.support("N1", **{"ux" if plane == "xv" else "uy": "fixed"})
    rot = "ry" if plane == "xv" else "rx"
    for nid, kind, is_first in (("N1", ends[0], True), (f"N{n_spans + 1}", ends[1], False)):
        if overhang and not is_first:
            continue
        if kind == "fixed":
            if partial and is_first:
                B.support(nid, **{rot: {"spring": B.k_rot(k_end)}})
            else:
                B.support(nid, **{rot: "fixed"})
    if spring:
        B.support(f"N{spring_at + 1}", uv={"spring": B.k_lin(k_spring)})
    tip = None
    if overhang:
        tip = B.node("T1", *xyv(stations[-1] + B.L(over_len, 0.5)))
        B.member("BC", f"N{n_spans + 1}", tip, g, roll=roll)
    B.restrain_out_of_plane(plane)

    # loads
    B.case("D", "dead load and self weight")
    B.self_weight("D")
    span_members = {k: chains.get(f"B{k + 1}", [f"B{k + 1}"]) for k in range(n_spans)}
    for k in range(n_spans):
        for mid in span_members[k]:
            B.line("D", mid, "v", -B.W(wd))
            B.line("L1" if k % 2 == 0 else "L2", mid, "v", -B.W(wl))
    if overhang:
        B.line("D", "BC", "v", -B.W(wd))
        B.nodal("L2" if n_spans % 2 == 0 else "L1", tip, Fv=-B.F(tip_load))
    B.case("L1", "live load on spans 1, 3, ...")
    if n_spans > 1 or overhang:
        B.case("L2", "live load on spans 2, 4, ...")
    live = [c for c in ("L1", "L2") if c in B.cases and B.cases[c]["loads"]]
    if points and not (hinge and load_span == hinge_span):
        B.point("L1" if load_span % 2 == 0 else "L2", f"B{load_span + 1}", "v", -B.F(load_p), load_at)
    else:
        points = False
    if varying and not (hinge and vary_span == hinge_span):
        case = "L1" if vary_span % 2 == 0 else "L2"
        B.line(case, f"B{vary_span + 1}", "v", -B.W(0.6 * wl), 0.0, start=0.5, end=1.0)
    else:
        varying = False
    if settle:
        B.case("S", "support settlement")
        B.settlement("S", f"N{settle_at + 1}", uv=-B.small(settlement))

    full = {c: 1.6 for c in live}
    cands = [("ULS1", {"D": 1.2, **full}, "strength"),
             ("SLS1", {"D": 1.0, **{c: 1.0 for c in live}}, "serviceability"),
             ("ULS2", {"D": 1.4}, "strength"),
             ("SLS2", {c: 1.0 for c in live}, "serviceability")]
    if len(live) == 2:
        cands += [("ULS3", {"D": 1.2, "L1": 1.6}, "strength"), ("ULS4", {"D": 1.2, "L2": 1.6}, "strength")]
    if settle:
        cands += [("ULS5", {"D": 1.2, **full, "S": 1.0}, "strength"), ("ULS6", {"D": 0.9, "S": 1.0}, "strength")]
    if len(live) == 2:
        cands += [("ULS7", {"D": 0.9, "L1": 1.6}, "strength"), ("ULS8", {"D": 0.9, "L2": 1.6}, "strength")]
    if settle and n_combos < 7 and level == 3:
        n_combos = 7                                         # so that the settlement case is used
    if len(cands) < n_combos:
        raise GenerationRejected("not enough load cases for the number of combinations")
    ids = add_combinations(B, cands, n_combos)
    used = {c for cid in ids for c in B.combinations[cid]["factors"]}
    for c in [c for c in B.cases if c not in used]:
        del B.cases[c]
    sls = service_ids(cands, n_combos)
    limits = [{"spanRatio": 240, "direction": "y", "combinations": [sls[0]]}]
    if len(sls) > 1:
        limits.append({"spanRatio": 360, "direction": "y", "combinations": [sls[1]]})
    B.design_group("beams", members, role="beam", bracing=dict(bracing[1]), deflectionLimits=limits)
    if overhang:
        B.design_member("BC", ["BC"], role="beam", Cb=1.0, deflectionLimits=[])
    B.parameters = {
        "plane": plane, "spans": lengths, "endSupports": ends, "deadLoad": B.W(wd), "liveLoad": B.W(wl),
        "nCombinations": n_combos, "bracing": bracing[0],
        "hinge": {"span": hinge_span + 1, "at": hinge_at} if hinge else None,
        "partialFixity": B.k_rot(k_end) if (partial and ends[0] == "fixed") else None,
        "springSupport": {"node": f"N{spring_at + 1}", "k": B.k_lin(k_spring)} if spring else None,
        "settlement": {"node": f"N{settle_at + 1}", "value": B.small(settlement)} if settle and "S" in B.cases else None,
        "overhang": B.L(over_len, 0.5) if overhang else None, "varyingLoad": varying, "pointLoads": points,
    }
    return B
