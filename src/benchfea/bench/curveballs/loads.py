"""Curve balls on loads: sign, direction, and the member or node a load is applied to."""
from __future__ import annotations

from ._util import is_number, member_axis, tidy
from .registry import Candidate, register

CATEGORY = "loads"
_NODAL = (("Fx", "Fy", "Fz"), ("Mx", "My", "Mz"))
_DISP = (("ux", "uy", "uz"), ("rx", "ry", "rz"))
_VALUE_KEYS = {"nodal": ("Fx", "Fy", "Fz", "Mx", "My", "Mz"), "memberPoint": ("value",),
               "memberDistributed": ("w1", "w2"), "supportDisplacement": ("ux", "uy", "uz", "rx", "ry", "rz")}


def _loads(doc):
    """(case id, index, load) of every load except self weight."""
    for cid, case in doc["loadCases"].items():
        for k, load in enumerate(case["loads"]):
            if load["type"] in _VALUE_KEYS:
                yield cid, k, load


def _describe(cid, k, load) -> str:
    where = load.get("member") or load.get("node")
    return f"load {k} of case {cid} ({load['type']} on {where})"


def _entities(cid, load, extra=()):
    key = "members" if "member" in load else "nodes"
    return {"loadCases": [cid], key: [load[key[:-1]], *extra]}


@register(
    name="load_sign_flipped", category=CATEGORY, title="load sign flipped",
    changes="Negates every component of one load (`/loadCases/<case>/loads/<k>/...`).",
    why="Loads are signed along the positive axes: gravity on a model with Z up is a negative Z value. A load "
        "with the wrong sign acts upwards or to the other side; in combination with the other loads it relieves "
        "where it should add, so the governing forces, the governing combination and even the sign of the "
        "reactions (uplift) change.",
    moves="Reactions, moments, shears and deflections of every combination that contains the case; which "
          "combination governs; design ratios.",
    effect_classes=("changes-analysis",),
    precondition="a nodal load, member load or imposed support displacement")
def load_sign_flipped(ctx):
    out = []
    for cid, k, load in _loads(ctx.doc):
        keys = [key for key in _VALUE_KEYS[load["type"]] if is_number(load.get(key)) and load[key] != 0]
        if not keys:
            continue

        def edit(d, cid=cid, k=k, keys=tuple(keys)):
            target = d["loadCases"][cid]["loads"][k]
            for key in keys:
                target[key] = tidy(-target[key])
        out.append(Candidate(_describe(cid, k, load), edit, _entities(cid, load)))
    return out


@register(
    name="load_direction_changed", category=CATEGORY, title="load applied in the wrong direction",
    changes="Turns one load to another axis of the same system: the `direction` of a member load (X/Y/Z or "
            "x/y/z), or the components of a nodal load moved one or two places (Fx->Fy->Fz, Mx->My->Mz) "
            "(`/loadCases/<case>/loads/<k>/...`).",
    why="The same magnitude along another axis is a different load case: gravity applied horizontally becomes a "
        "sway load, a load along the member axis produces axial force instead of bending, and a load out of the "
        "plane of a plane frame goes straight into the out-of-plane restraints or bends the members about their "
        "weak axis. With Y-up and Z-up conventions both in use this is an easy slip.",
    moves="Which members bend and about which axis, axial forces, reactions (direction and magnitude), "
          "deflections; design ratios.",
    effect_classes=("changes-analysis",),
    precondition="a nodal load or member load")
def load_direction_changed(ctx):
    out = []
    for cid, k, load in _loads(ctx.doc):
        if load["type"] in ("memberPoint", "memberDistributed"):
            family = "XYZ" if load["direction"] in "XYZ" else "xyz"
            for new in family:
                if new != load["direction"]:
                    out.append(Candidate(
                        f"{_describe(cid, k, load)}: direction {load['direction']} -> {new}",
                        lambda d, cid=cid, k=k, new=new: d["loadCases"][cid]["loads"][k].__setitem__("direction", new),
                        _entities(cid, load)))
        elif load["type"] == "nodal":
            for shift in (1, 2):
                new = {}
                for group in _NODAL:
                    for i, key in enumerate(group):
                        if key in load:
                            new[group[(i + shift) % 3]] = load[key]
                if all(load.get(key) == v for key, v in new.items()) and len(new) == sum(key in load for g in _NODAL for key in g):
                    continue

                def edit(d, cid=cid, k=k, new=new):
                    old = d["loadCases"][cid]["loads"][k]
                    rest = {key: v for key, v in old.items() if key not in _NODAL[0] and key not in _NODAL[1]}
                    comps = {key: new[key] for g in _NODAL for key in g if key in new}
                    d["loadCases"][cid]["loads"][k] = {**rest, **comps}
                names = ", ".join(f"{key}" for g in _NODAL for key in g if key in new)
                out.append(Candidate(f"{_describe(cid, k, load)}: components moved {shift} axis place(s) (now {names})",
                                     edit, _entities(cid, load)))
    return out


@register(
    name="load_moved", category=CATEGORY, title="load moved to the wrong member or node",
    changes="Applies one load to another member (member loads) or another node (nodal loads) with its values "
            "kept (`/loadCases/<case>/loads/<k>/member` or `/node`).",
    why="The total load on the structure is unchanged, so the sum of the reactions still checks -- but the load "
        "path does not: the member that should carry the load is unloaded, another carries twice its share, and "
        "the support reactions shift towards the wrong bay or storey.",
    moves="Moments, shears and deflections of the two members involved, column axial forces and individual "
          "reactions (their sum stays), design ratios of both members.",
    effect_classes=("changes-analysis",),
    precondition="a nodal load or member load, and another node / member to put it on")
def load_moved(ctx):
    model = ctx.model
    member_ids, node_ids = list(model.members), list(model.nodes)
    connected = {n for m in model.members.values() for n in (m.i, m.j)}
    out = []
    for cid, k, load in _loads(ctx.doc):
        if "member" in load:
            a = load["member"]
            start = member_ids.index(a)
            ring = member_ids[start + 1:] + member_ids[:start]
            axis = member_axis(model, a)
            same = [b for b in ring if model.members[b].type == model.members[a].type]
            parallel = [b for b in same if abs(float(member_axis(model, b) @ axis)) > 0.99]
            for b in (parallel or same)[:2]:
                out.append(Candidate(f"{_describe(cid, k, load)}: moved to member {b}",
                                     lambda d, cid=cid, k=k, b=b: d["loadCases"][cid]["loads"][k].__setitem__("member", b),
                                     _entities(cid, load, [b])))
        elif load["type"] == "nodal":
            a = load["node"]
            start = node_ids.index(a)
            ring = [n for n in node_ids[start + 1:] + node_ids[:start] if n in connected]
            for b in ring[:2]:
                out.append(Candidate(f"{_describe(cid, k, load)}: moved to node {b}",
                                     lambda d, cid=cid, k=k, b=b: d["loadCases"][cid]["loads"][k].__setitem__("node", b),
                                     _entities(cid, load, [b])))
    return out
