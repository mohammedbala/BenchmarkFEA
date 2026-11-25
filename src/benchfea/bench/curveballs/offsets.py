"""Curve balls on rigid end offsets."""
from __future__ import annotations

import numpy as np

from ._util import frame_ends_at, member_axis, section_depth, tidy
from .registry import Candidate, register

CATEGORY = "offsets"
_WHY = ("A rigid end offset says where the flexible length of a member starts: the part inside the joint (half "
        "the depth of the column or beam it frames into) does not bend. It shortens the flexible span and moves "
        "the point where the member end forces act on the joint, which adds the moment V x e to the joint.")


def _vectors(member: dict) -> dict:
    o = member.get("offsets") or {}
    return {e: list(o[e]) for e in ("i", "j") if e in o and any(c != 0 for c in o[e])}


@register(
    name="offset_added", category=CATEGORY, title="rigid end offset added",
    changes="Adds a rigid end offset at one end of a member that frames into another member: half the section "
            "depth of that member, along the member axis into the span (`/members/<id>/offsets`).",
    why=_WHY + " An offset that the structure does not have makes the member stiffer and shorter than it is.",
    moves="Member end moments and shears at the joint, span moments, joint rotations and deflections (smaller), "
          "sway stiffness of the frame.",
    effect_classes=("changes-analysis",),
    precondition="a frame member end without an offset at a joint with a non-parallel frame member whose section "
                 "has a depth (I, box, pipe or rectangle section)")
def offset_added(ctx):
    doc, model = ctx.doc, ctx.model
    out = []
    for mid, m in model.members.items():
        if m.type != "frame":
            continue
        dm = doc["members"][mid]
        have = _vectors(dm)
        axis = member_axis(model, mid)
        for end, node in (("i", m.i), ("j", m.j)):
            if end in have:
                continue
            depth, other = 0.0, None
            for oid, _ in frame_ends_at(model, node):
                if oid == mid or abs(float(member_axis(model, oid) @ axis)) > 0.9:
                    continue
                d = section_depth(doc["sections"][doc["members"][oid]["section"]])
                if d and d > depth:
                    depth, other = float(d), oid
            length = float(model.units.from_internal(m.length, "length"))
            if other is None or depth / 2.0 >= 0.25 * length:
                continue
            sign = 1.0 if end == "i" else -1.0
            axes = (dm.get("offsets") or {}).get("axes", "global") if "offsets" in dm else "local"
            vec = tidy([sign * depth / 2.0, 0, 0]) if axes == "local" else tidy(list(sign * depth / 2.0 * axis))

            def edit(d, mid=mid, end=end, vec=vec, axes=axes):
                mem = d["members"][mid]
                if "offsets" not in mem:
                    mem["offsets"] = {"axes": axes}
                mem["offsets"][end] = list(vec)
            out.append(Candidate(f"member {mid}, end {end} (half the depth of {other})", edit,
                                 {"members": [mid], "nodes": [node]}, {"offsetLength": tidy(depth / 2.0)}))
    return out


@register(
    name="offset_removed", category=CATEGORY, title="rigid end offset removed",
    changes="Deletes the rigid end offsets of a member, at both ends or at one (`/members/<id>/offsets`).",
    why=_WHY + " Without it the member spans centre to centre: it is longer and more flexible, and the joint "
        "loses the eccentricity moment.",
    moves="Span and end moments (larger span), deflections and drift (larger), joint moments.",
    effect_classes=("changes-analysis",),
    precondition="a member with a nonzero rigid end offset")
def offset_removed(ctx):
    out = []
    for mid, dm in ctx.doc["members"].items():
        have = _vectors(dm)
        if not have:
            continue
        m = ctx.model.members[mid]
        nodes = {"i": m.i, "j": m.j}
        out.append(Candidate(f"member {mid}: all offsets", lambda d, mid=mid: d["members"][mid].pop("offsets"),
                             {"members": [mid], "nodes": [nodes[e] for e in have]}))
        if len(have) == 2:
            for end in have:
                out.append(Candidate(f"member {mid}: offset at end {end}",
                                     lambda d, mid=mid, end=end: d["members"][mid]["offsets"].pop(end),
                                     {"members": [mid], "nodes": [nodes[end]]}))
    return out


@register(
    name="offset_sign_reversed", category=CATEGORY, title="rigid end offset with reversed sign",
    changes="Negates the offset vectors of a member (`/members/<id>/offsets/i`, `/j`).",
    why=_WHY + " With the sign reversed the flexible end lies on the wrong side of the node: the flexible "
        "length grows by twice the offset instead of shrinking, and the eccentricity moment changes sign.",
    moves="Span length and with it end moments, shears and deflections; joint moments from the eccentricity.",
    effect_classes=("changes-analysis",),
    precondition="a member with a nonzero rigid end offset")
def offset_sign_reversed(ctx):
    out = []
    for mid, dm in ctx.doc["members"].items():
        have = _vectors(dm)
        if not have:
            continue

        def edit(d, mid=mid, have=have):
            for end, v in have.items():
                d["members"][mid]["offsets"][end] = tidy([-c for c in v])
        m = ctx.model.members[mid]
        out.append(Candidate(f"member {mid}", edit, {"members": [mid], "nodes": [{"i": m.i, "j": m.j}[e] for e in have]}))
    return out


@register(
    name="offset_wrong_axis", category=CATEGORY, title="rigid end offset on the wrong axis",
    changes="Puts the offset of a member on another axis: the vector components are rotated one or two places "
            "(X->Y->Z), or the `axes` flag is switched between global and local with the numbers kept "
            "(`/members/<id>/offsets`).",
    why=_WHY + " On the wrong axis the offset no longer shortens the member along its axis: it shifts the "
        "flexible end sideways or vertically, which introduces an eccentricity (axial force x offset = moment, "
        "shear x offset = torsion) and changes the member's direction slightly.",
    moves="Member end moments and torsion, joint equilibrium moments, deflections; out-of-plane forces in a "
          "plane frame.",
    effect_classes=("changes-analysis",),
    precondition="a member with a nonzero rigid end offset")
def offset_wrong_axis(ctx):
    out = []
    for mid, dm in ctx.doc["members"].items():
        have = _vectors(dm)
        if not have:
            continue
        m = ctx.model.members[mid]
        ent = {"members": [mid], "nodes": [{"i": m.i, "j": m.j}[e] for e in have]}
        for shift in (1, 2):
            new = {e: [v[(k - shift) % 3] for k in range(3)] for e, v in have.items()}
            if new == have:
                continue

            def edit(d, mid=mid, new=new):
                d["members"][mid]["offsets"].update({e: list(v) for e, v in new.items()})
            out.append(Candidate(f"member {mid}: components moved {shift} axis place(s)", edit, ent))
        axes = dm["offsets"].get("axes", "global")
        other = "local" if axes == "global" else "global"
        r = np.asarray(m.R)
        if not all(np.allclose(r @ np.asarray(v, float), np.asarray(v, float)) for v in have.values()):
            out.append(Candidate(f"member {mid}: axes {axes} -> {other}",
                                 lambda d, mid=mid, other=other: d["members"][mid]["offsets"].__setitem__("axes", other), ent))
    return out
