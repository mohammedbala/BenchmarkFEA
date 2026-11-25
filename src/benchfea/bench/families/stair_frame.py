"""Stair frames: inclined stringers and landings.

Levels 1-2: one stringer line in a vertical plane (plane model in 3D space) -- lower flight, flat
landing, upper flight, optionally a post under each end of the landing.  Levels 3-4: a dog-leg
stair in 3D -- two flights of two stringers each, running in opposite directions, and a landing
frame (front beam, two side beams, back beam) on posts.  Floor ends are pinned.  Dead and live
load are area loads on the plan projection (``projected`` line loads on the inclined stringers).
"""
from __future__ import annotations

from ._common import BUCKLING_CHANCE, BUCKLING_FLAG, Builder, add_combinations, flag, service_ids, tidy

NAME = "stair_frame"
LEVELS = (1, 2, 3, 4)
FEATURES = {
    "buckling": BUCKLING_FLAG + " (stairs with posts)",
    "posts": "posts under the landing (always present in the dog-leg stair)",
    "hingedFlights": "stringers pinned (bending released) where they meet the landing",
    "springSupport": "elastic vertical support at the lower floor end (dog-leg) / under the landing (plane stair)",
    "frontPosts": "dog-leg stair: posts under the front corners of the landing as well",
}
_TREADS = {"topFlange": [0.25, 0.5, 0.75], "bottomFlange": [0.25, 0.5, 0.75]}


def build(draw, level: int, unit_system: str, vertical: str, features: dict) -> Builder:
    B = Builder(unit_system, vertical, NAME)
    plane = draw.choice(("xv", "yv"))
    run = B.L(draw.step(2.5, 4.0, 0.25), 0.25)
    rise = B.L(draw.step(1.5, 2.0, 0.25), 0.25)                # per flight
    landing = B.L(draw.step(1.25, 2.0, 0.25), 0.25)
    width_m = draw.step(1.0, 1.5, 0.25)
    posts = flag(features, "posts", draw.chance({1: 0.0, 2: 0.8, 3: 1.0, 4: 1.0}[level]))
    hinged = flag(features, "hingedFlights", draw.chance({1: 0.0, 2: 0.5, 3: 0.3, 4: 0.8}[level]))
    spring = flag(features, "springSupport", draw.chance({1: 0.0, 2: 0.3, 3: 0.2, 4: 0.7}[level]))
    front = flag(features, "frontPosts", draw.chance({1: 0.0, 2: 0.0, 3: 0.3, 4: 0.7}[level]))
    n_combos = {1: 2, 2: draw.integer(3, 5), 3: draw.integer(3, 5), 4: draw.integer(6, 7)}[level]
    qd, ql = draw.uniform(2.5, 4.5), draw.uniform(3.0, 5.0)     # kPa on plan
    k_spring = draw.uniform(5, 30)
    base = draw.choice(("fixed", "pinned"))
    kind = draw.choice(("box", "beam"))
    per_m = B.us["per_m"]
    B.case("D", "dead load and self weight")
    B.self_weight("D")
    B.case("L1", "live load on the flights")
    B.case("L2", "live load on the landing")
    pin = {"rz": True}
    if level <= 2:
        def xyv(s, v):
            return (s, 0.0, v) if plane == "xv" else (0.0, s, v)

        roll = 0.0
        g = B.group("stringer", kind, 4)
        B.node("N1", *xyv(0.0, 0.0))
        B.node("N2", *xyv(run, rise))
        B.node("N3", *xyv(tidy(run + landing), rise))
        B.node("N4", *xyv(tidy(2 * run + landing), tidy(2 * rise)))
        B.member("S1", "N1", "N2", g, roll=roll, releases={"j": pin} if hinged and posts else None)
        B.member("S2", "N2", "N3", g, roll=roll)
        B.member("S3", "N3", "N4", g, roll=roll, releases={"i": pin} if hinged and posts else None)
        along, rot = ("ux", "ry") if plane == "xv" else ("uy", "rx")
        B.support("N1", "pinned")
        B.support("N4", uv="fixed", **({along: "fixed"} if level == 2 else {}))
        post_members = []
        if posts:
            gp = B.group("post", "box", 4)
            for n, pid in (("N2", "P1"), ("N3", "P2")):
                b = B.node("B" + pid[1], B.nodes[n][0], B.nodes[n][1], 0.0)
                post_members.append(B.member(pid, b, n, gp, roll=roll))
                B.support(b, "pinned")
                if base == "fixed":
                    B.support(b, **{rot: "fixed"})
            if spring:
                B.support("B2", uv={"spring": B.k_lin(k_spring)})
        B.restrain_out_of_plane(plane)
        trib = width_m / 2.0
        for mid, case in (("S1", "L1"), ("S2", "L2"), ("S3", "L1")):
            B.line("D", mid, "v", -B.W(qd * trib), projected=mid != "S2")
            B.line(case, mid, "v", -B.W(ql * trib), projected=mid != "S2")
        stringers, landing_beams = ["S1", "S3"], ["S2"]
        hinged = hinged and posts
        spring = spring and posts
        front = False
    else:
        gap = B.L(0.25, 0.25)
        w = B.L(width_m, 0.25)
        ys = [0.0, w, tidy(w + gap), tidy(2 * w + gap)]
        kind = "beam"                                         # open sections: little compatibility torsion
        g = B.group("stringer", kind, 4)
        gl, gp = B.group("landingBeam", "beam", 5), B.group("post", "box", 5)
        for j, y in enumerate(ys):
            B.node(f"P{j}", run, y, rise)
        for j in (0, 1):
            B.node(f"F{j}", 0.0, ys[j], 0.0)
            B.support(f"F{j}", "pinned")
        for j in (2, 3):
            B.node(f"T{j}", 0.0, ys[j], tidy(2 * rise))
            B.support(f"T{j}", "pinned")
        B.node("Q0", tidy(run + landing), ys[0], rise)
        B.node("Q3", tidy(run + landing), ys[3], rise)
        rel = {"j": pin} if hinged else None
        stringers = [B.member(f"S{j}", f"F{j}", f"P{j}", g, releases=rel) for j in (0, 1)]
        stringers += [B.member(f"S{j}", f"T{j}", f"P{j}", g, releases=rel) for j in (2, 3)]
        landing_beams = [B.member(f"LF{j}", f"P{j}", f"P{j + 1}", gl) for j in range(3)]
        landing_beams += [B.member("LS0", "P0", "Q0", gl), B.member("LS3", "P3", "Q3", gl),
                          B.member("LB", "Q0", "Q3", gl)]
        post_members = []
        for n in ["Q0", "Q3"] + (["P0", "P3"] if front else []):
            b = B.node("B" + n, B.nodes[n][0], B.nodes[n][1], 0.0)
            post_members.append(B.member("C" + n, b, n, gp))
            B.support(b, base)
        if spring:
            B.support("F0", uv={"spring": B.k_lin(k_spring)})
            B.support("F1", uv={"spring": B.k_lin(k_spring)})
        trib = width_m / 2.0
        for mid in stringers:
            B.line("D", mid, "v", -B.W(qd * trib), projected=True)
            B.line("L1", mid, "v", -B.W(ql * trib), projected=True)
        half = landing / per_m / 2.0
        for mid in ("LF0", "LF1", "LF2", "LB"):
            B.line("D", mid, "v", -B.W(qd * half))
            B.line("L2", mid, "v", -B.W(ql * half))
        posts = True
    cands = [("ULS1", {"D": 1.2, "L1": 1.6, "L2": 1.6}, "strength"),
             ("SLS1", {"D": 1.0, "L1": 1.0, "L2": 1.0}, "serviceability"),
             ("ULS2", {"D": 1.4}, "strength"),
             ("ULS3", {"D": 1.2, "L1": 1.6}, "strength"),
             ("ULS4", {"D": 1.2, "L2": 1.6}, "strength"),
             ("SLS2", {"L1": 1.0, "L2": 1.0}, "serviceability"),
             ("ULS5", {"D": 0.9, "L1": 1.6}, "strength")]
    add_combinations(B, cands, n_combos)
    sls = service_ids(cands, n_combos)
    buckling = B.request_buckling(flag(features, "buckling", draw.chance(BUCKLING_CHANCE[level])) and bool(post_members))
    limits = [{"spanRatio": 360, "direction": "y", "combinations": [sls[0]]}]
    B.design_group("stringers", stringers, role="beam", bracing=dict(_TREADS), deflectionLimits=limits)
    B.design_group("landingBeams", landing_beams, role="beam", deflectionLimits=limits)
    if post_members:
        B.design_group("posts", post_members, role="column", Kx=1.0 if level <= 2 else 1.2, Ky=1.0 if level <= 2 else 1.2)
    B.parameters = {
        "buckling": buckling,
        "layout": "straight, plane " + plane if level <= 2 else "dog-leg", "run": run, "risePerFlight": rise,
        "landing": landing, "width_m": width_m, "stringerSection": kind, "posts": len(post_members),
        "postBases": base if post_members else None, "hingedFlights": hinged, "springSupport": spring,
        "frontPosts": front, "deadLoad_kPa": round(qd, 3), "liveLoad_kPa": round(ql, 3), "nCombinations": n_combos,
    }
    return B
