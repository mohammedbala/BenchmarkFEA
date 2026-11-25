"""Curve balls on rigid links and rigid diaphragms."""
from __future__ import annotations

import numpy as np

from ._util import AXIS_OF, DOFS, level_groups, new_id, support_dofs
from .registry import Candidate, register

CATEGORY = "constraints"
_WHY_DIAPHRAGM = ("A rigid diaphragm or link is a kinematic assumption that decides how lateral load is shared: "
                  "with it the connected nodes move as one rigid body and the load goes to the frames in "
                  "proportion to their stiffness and position; without it each node moves on its own and the "
                  "load follows the members.")


def _drop_constraint(doc, cid):
    # the (possibly empty) collection is kept, so that the recorded path names the constraint
    del doc["constraints"][cid]


def _set_nodes(doc) -> set:
    out = set()
    for c in doc.get("constraints", {}).values():
        out.add(c["master"])
        out.update(c["slaves"])
    return out


@register(
    name="constraint_removed", category=CATEGORY, title="rigid link / diaphragm removed",
    changes="Deletes one rigid constraint (`/constraints/<id>`).",
    why=_WHY_DIAPHRAGM + " Removing it frees the slaves from the master: storey forces applied at the master no "
        "longer reach the frame (a master without members is then a loaded node with no stiffness -- a mechanism), "
        "or the floor racks and twists instead of moving rigidly.",
    moves="Storey drifts and floor rotation, the split of storey shear between frames, column and beam end "
          "moments; or no results at all (mechanism at the former master).",
    effect_classes=("changes-analysis", "creates-mechanism"),
    precondition="the model has at least one rigid constraint (rigidLink or rigidDiaphragm)")
def constraint_removed(ctx):
    out = []
    for cid, c in ctx.doc.get("constraints", {}).items():
        out.append(Candidate(f"constraint {cid}", lambda d, cid=cid: _drop_constraint(d, cid),
                             {"constraints": [cid], "nodes": [c["master"], *c["slaves"]]}))
    return out


@register(
    name="constraint_added", category=CATEGORY, title="rigid link / diaphragm added",
    changes="Adds a rigid constraint that the structure does not have: a rigid diaphragm over the nodes of one "
            "level, or a rigid link between the two end nodes of a member (`/constraints/<new id>`).",
    why="A rigid constraint removes deformation the real structure has. A diaphragm that is not there forces all "
        "nodes of a level to translate and twist together, so flexible frames are propped by stiff ones; a rigid "
        "link across a member makes that member infinitely stiff, so it attracts load and stops deforming.",
    moves="Drifts (smaller), distribution of lateral load between frames, forces in the member bridged by a "
          "link (its deformation vanishes), reactions.",
    effect_classes=("changes-analysis",),
    precondition="a level with at least two nodes that are not supported in the DOFs a diaphragm ties and are "
                 "not part of a constraint, or a frame member whose end node is unsupported and unconstrained "
                 "(plane frames modelled with out-of-plane supports at every node do not qualify)")
def constraint_added(ctx):
    doc, model = ctx.doc, ctx.model
    existing = doc.get("constraints", {})
    taken = _set_nodes(doc)
    v = AXIS_OF[ctx.vertical]
    tied = [d for k, d in enumerate(DOFS[:3]) if k != v] + [DOFS[3 + v]]
    supports = doc.get("supports", {})

    def free_in(nid, dofs):
        s = supports.get(nid)
        return s is None or all(support_dofs(s, ctx.vertical)[d] == "free" for d in dofs)

    out = []
    for n, (_, nodes) in enumerate(level_groups(model).items()):
        ok = [x for x in nodes if x not in taken and free_in(x, tied) and
              any(x in (m.i, m.j) for m in model.members.values())]
        if len(ok) < 2 or len(ok) != len(nodes):
            continue
        cid = new_id(existing, f"diaphragm-L{n}")

        def edit(d, cid=cid, ok=tuple(ok)):
            d.setdefault("constraints", {})[cid] = {"type": "rigidDiaphragm", "master": ok[0], "slaves": list(ok[1:])}
        out.append(Candidate(f"rigid diaphragm over level {n} (master {ok[0]})", edit,
                             {"constraints": [cid], "nodes": list(ok)}, {"kind": "rigidDiaphragm"}))
    for mid, m in model.members.items():
        if m.type != "frame" or m.j in taken or m.i in taken or not free_in(m.j, DOFS):
            continue
        cid = new_id(existing, f"link-{mid}")

        def edit(d, cid=cid, m=m):
            d.setdefault("constraints", {})[cid] = {"type": "rigidLink", "master": m.i, "slaves": [m.j]}
        out.append(Candidate(f"rigid link {m.i} -> {m.j} (across member {mid})", edit,
                             {"constraints": [cid], "nodes": [m.i, m.j], "members": [mid]}, {"kind": "rigidLink"}))
    return out


@register(
    name="constraint_slaves_dropped", category=CATEGORY, title="constraint attached to too few nodes",
    changes="Shortens the slave list of a rigid constraint (`/constraints/<id>/slaves`): one slave left out, half "
            "of them, or all but one.",
    why=_WHY_DIAPHRAGM + " A node left out of the list is no longer carried by the diaphragm: the master's load "
        "enters the structure only through the remaining slaves and has to travel through the members to reach "
        "the frames that were disconnected.",
    moves="Split of storey shear between frames, floor rotation, axial forces in the beams that now act as "
          "collectors, drift of the disconnected nodes.",
    effect_classes=("changes-analysis",),
    precondition="a rigid constraint with at least two slaves")
def constraint_slaves_dropped(ctx):
    out = []
    for cid, c in ctx.doc.get("constraints", {}).items():
        slaves = list(c["slaves"])
        if len(slaves) < 2:
            continue
        variants = [(f"slave {s} left out", [x for x in slaves if x != s]) for s in slaves]
        if len(slaves) >= 3:
            variants += [(f"only slave {s} kept", [s]) for s in slaves]
        if len(slaves) >= 4:
            h = len(slaves) // 2
            variants += [("second half of the slaves left out", slaves[:h]), ("first half of the slaves left out", slaves[h:])]
        for label, keep in variants:
            dropped = [s for s in slaves if s not in keep]
            out.append(Candidate(f"constraint {cid}: {label}",
                                 lambda d, cid=cid, keep=tuple(keep): d["constraints"][cid].__setitem__("slaves", list(keep)),
                                 {"constraints": [cid], "nodes": dropped}))
    return out


def _nearest_outside(model, node: str, inside: set, count: int = 2) -> list:
    x = model.nodes[node].xyz
    others = [(float(np.linalg.norm(n.xyz - x)), k, nid) for k, (nid, n) in enumerate(model.nodes.items())
              if nid not in inside]
    return [nid for _, _, nid in sorted(others)[:count]]


@register(
    name="constraint_wrong_slave", category=CATEGORY, title="constraint attached to a wrong node",
    changes="Replaces one slave of a rigid constraint by a nearby node that does not belong to it "
            "(`/constraints/<id>/slaves/<k>`).",
    why="The constraint now ties a node of another part of the structure (typically the node above or below on "
        "the same column line) to the master and lets the intended node go. Where that node is supported or "
        "already belongs to another constraint the model is contradictory and is refused; otherwise two levels "
        "are locked together and the intended node drifts freely.",
    moves="Either no analysis (validation refuses a slave that is supported in a tied DOF or has two masters), "
          "or storey drifts, column shears and moments between the two levels that were locked together.",
    effect_classes=("changes-analysis", "invalid-model"),
    precondition="the model has a rigid constraint and a node outside it")
def constraint_wrong_slave(ctx):
    out = []
    for cid, c in ctx.doc.get("constraints", {}).items():
        inside = {c["master"], *c["slaves"]}
        for k, s in enumerate(c["slaves"]):
            for n in _nearest_outside(ctx.model, s, inside):
                out.append(Candidate(f"constraint {cid}: slave {s} replaced by {n}",
                                     lambda d, cid=cid, k=k, n=n: d["constraints"][cid]["slaves"].__setitem__(k, n),
                                     {"constraints": [cid], "nodes": [s, n]}))
    return out


@register(
    name="constraint_wrong_master", category=CATEGORY, title="constraint with the wrong master node",
    changes="Changes the master of a rigid constraint (`/constraints/<id>/master`): one of its slaves is made "
            "the master (and leaves the slave list), or a nearby node outside the constraint is.",
    why="Loads and masses of a floor are applied at the master. If another node is declared master, the node "
        "that carries the storey loads is no longer part of the rigid body: a memberless master is left as a "
        "loaded node with no stiffness (a mechanism), and a master taken from outside ties the floor to a node "
        "that should move independently.",
    moves="No results (mechanism at the abandoned master, reported by the diagnosis), or drifts and the load "
          "path of the storey forces.",
    effect_classes=("changes-analysis", "creates-mechanism"),
    precondition="the model has a rigid constraint")
def constraint_wrong_master(ctx):
    out = []
    for cid, c in ctx.doc.get("constraints", {}).items():
        old = c["master"]
        for s in c["slaves"]:
            def edit(d, cid=cid, s=s):
                cc = d["constraints"][cid]
                cc["master"] = s
                cc["slaves"] = [x for x in cc["slaves"] if x != s]
            if len(c["slaves"]) >= 2:
                out.append(Candidate(f"constraint {cid}: slave {s} made master, {old} detached", edit,
                                     {"constraints": [cid], "nodes": [old, s]}))
        for n in _nearest_outside(ctx.model, old, {old, *c["slaves"]}):
            out.append(Candidate(f"constraint {cid}: master {old} replaced by {n}",
                                 lambda d, cid=cid, n=n: d["constraints"][cid].__setitem__("master", n),
                                 {"constraints": [cid], "nodes": [old, n]}))
    return out


@register(
    name="constraint_made_flexible", category=CATEGORY, title="rigid constraint replaced by flexible members",
    changes="Deletes a rigid constraint and adds one frame member of a section of the model from the master to "
            "each slave (`/constraints/<id>` removed, `/members/<id>-flex-<k>` added).",
    why="Members of finite stiffness do not enforce rigid-body motion: the connected nodes can move relative to "
        "each other, so load is shared according to the flexibility of the connecting members instead of the "
        "rigid-body rule. Unlike a diaphragm, the stand-in members also couple the out-of-plane DOFs, and they "
        "have self weight and are design members like any other.",
    moves="Floor rotation and drifts, split of storey shear between frames, forces in the columns under the "
          "floor; the new members carry forces of their own.",
    effect_classes=("changes-analysis",),
    precondition="the model has a rigid constraint whose master does not coincide with a slave")
def constraint_made_flexible(ctx):
    doc, model = ctx.doc, ctx.model
    materials = [m["material"] for m in doc["members"].values()]
    material = max(dict.fromkeys(materials), key=materials.count)
    out = []
    for cid, c in doc.get("constraints", {}).items():
        xm = model.nodes[c["master"]].xyz
        if any(float(np.linalg.norm(model.nodes[s].xyz - xm)) <= 1e-9 for s in c["slaves"]):
            continue
        for sec in doc["sections"]:
            ids = []
            for k, _ in enumerate(c["slaves"], start=1):
                ids.append(new_id(set(doc["members"]) | set(ids), f"{cid}-flex-{k}"))

            def edit(d, cid=cid, c=c, sec=sec, ids=tuple(ids)):
                _drop_constraint(d, cid)
                for mid, s in zip(ids, c["slaves"]):
                    d["members"][mid] = {"i": c["master"], "j": s, "section": sec, "material": material}
            out.append(Candidate(f"constraint {cid} replaced by {len(ids)} members of section {sec}", edit,
                                 {"constraints": [cid], "members": list(ids), "nodes": [c["master"], *c["slaves"]],
                                  "sections": [sec]}, {"standInSection": sec}))
    return out
