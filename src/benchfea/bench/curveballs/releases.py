"""Curve balls on member end releases."""
from __future__ import annotations

from ._util import frame_ends_at, release_state, set_release, vertical_member
from .registry import Candidate, register

CATEGORY = "releases"


def _connected(state, springs: bool) -> bool:
    """Is this DOF still (partly) connected?  ``springs``: a connection spring counts as connected."""
    return state is False or (springs and state is not True)


def _release(doc, ends, dofs, springs: bool = False) -> None:
    """Fully release ``dofs`` at the member ends (``springs``: connection springs are released too)."""
    for mid, end in ends:
        for dof in dofs:
            if _connected(release_state(doc["members"][mid], end, dof), springs):
                set_release(doc["members"][mid], end, dof, True)


def _open(doc, ends, dofs, springs: bool = False) -> list:
    """The (member, end) of ``ends`` where at least one of ``dofs`` is still connected."""
    return [(mid, end) for mid, end in ends
            if any(_connected(release_state(doc["members"][mid], end, d), springs) for d in dofs)]


@register(
    name="release_added", category=CATEGORY, title="moment release added",
    changes="Releases both bending rotations (local ry and rz) at one end of a frame member that was rigidly "
            "connected (`/members/<id>/releases/<end>`).",
    why="A released end cannot transmit bending moment. The moment the connection carried has to go somewhere "
        "else: span moments grow, the joint and the members framing into it lose that moment, and a frame that "
        "resisted sway through that connection becomes more flexible.",
    moves="Moment at the released end (to zero), span moment and deflection of the member, moments of the "
          "neighbouring members, sway and support reactions of an indeterminate structure.",
    effect_classes=("changes-analysis",),
    precondition="a frame member end whose bending rotations are rigidly connected")
def release_added(ctx):
    out = []
    for mid, m in ctx.model.members.items():
        if m.type != "frame":
            continue
        for end, node in (("i", m.i), ("j", m.j)):
            if _open(ctx.doc, [(mid, end)], ("ry", "rz")):
                out.append(Candidate(f"member {mid}, end {end}",
                                     lambda d, mid=mid, end=end: _release(d, [(mid, end)], ("ry", "rz")),
                                     {"members": [mid], "nodes": [node]}))
    return out


@register(
    name="release_removed", category=CATEGORY, title="release removed (connection made rigid)",
    changes="Deletes the releases at one end of a member -- full releases and connection springs alike "
            "(`/members/<id>/releases/<end>`).",
    why="A pinned or semi-rigid connection modelled as rigid attracts moment it cannot carry in reality: the "
        "member becomes continuous with the joint, its span moment drops, and the joint, the column or the "
        "neighbouring span receive a moment they were not designed for.",
    moves="End moment of the member (from zero or a spring value to the rigid value), span moment and deflection "
          "(smaller), moments in the members at the joint, frame sway (smaller).",
    effect_classes=("changes-analysis",),
    precondition="a member with a release or connection spring at an end")
def release_removed(ctx):
    out = []
    for mid, dm in ctx.doc["members"].items():
        m = ctx.model.members[mid]
        for end, node in (("i", m.i), ("j", m.j)):
            rel = (dm.get("releases") or {}).get(end) or {}
            if not any(v is not False for v in rel.values()):
                continue

            def edit(d, mid=mid, end=end):
                # the (possibly empty) `releases` object is kept, so that the recorded path names the end
                del d["members"][mid]["releases"][end]
            out.append(Candidate(f"member {mid}, end {end}", edit, {"members": [mid], "nodes": [node]}))
    return out


@register(
    name="release_mechanism", category=CATEGORY, title="releases that create a mechanism",
    changes="Releases the rotations of EVERY frame member end at one node (local ry and rz, or rx, ry and rz), "
            "leaving the node itself without rotational stiffness (`/members/<id>/releases/<end>` of each "
            "member at the node).",
    why="Releases belong to member ends, the rotation belongs to the node. If every member end at a node is "
        "released, nothing holds the node's rotation: the stiffness matrix is singular. The structure may still "
        "look sensible (each member is pinned), but there is no solution; the correct answer is the diagnosis of "
        "the mechanism, and the fix is to keep one member end connected.",
    moves="There are no results: the reference is the solver's mechanism diagnosis (node, DOFs, cause "
          "`all-member-ends-released`, members involved).",
    effect_classes=("creates-mechanism",),
    precondition="a node with frame member ends whose rotations are still connected and whose own rotations are "
                 "not all held by a support")
def release_mechanism(ctx):
    out = []
    for node in ctx.model.nodes:
        ends = frame_ends_at(ctx.model, node)
        if not ends:
            continue
        s = ctx.model.supports.get(node)
        if s is not None and all(s.fixed[3:]):
            continue
        for dofs, label in ((("ry", "rz"), "bending rotations"), (("rx", "ry", "rz"), "all rotations")):
            todo = _open(ctx.doc, ends, dofs, True)
            if not todo:
                continue
            out.append(Candidate(f"node {node}: {label} released at every member end "
                                 f"({', '.join(f'{m}.{e}' for m, e in ends)})",
                                 lambda d, ends=tuple(ends), dofs=dofs: _release(d, ends, dofs, True),
                                 {"nodes": [node], "members": sorted({m for m, _ in todo})},
                                 {"releasedDofs": list(dofs)}))
    return out


@register(
    name="release_pinned_frame", category=CATEGORY, title="moment frame pinned at every beam-column joint",
    changes="Fully releases the bending rotations (rz, or ry and rz; connection springs included) at every beam "
            "end that frames into a column "
            "(`/members/<beam>/releases/<end>` of all beams).",
    why="A moment frame resists sway by the bending of beams and columns through rigid joints. With every "
        "beam-column joint pinned the beams become simply supported and no longer restrain the columns: on "
        "pinned bases the frame is a sway mechanism (no solution, the answer is the diagnosis); on fixed or "
        "partially fixed bases the columns act as cantilevers from the base and the drift and base moments "
        "grow several times.",
    moves="Mechanism (pinned bases), or: storey drift and column base moments (much larger), beam end moments "
          "(zero), beam span moments (simply supported values).",
    effect_classes=("creates-mechanism", "changes-analysis"),
    precondition="non-vertical frame members (beams) with a rigidly connected end at a node where a vertical "
                 "frame member (column) ends")
def release_pinned_frame(ctx):
    model = ctx.model
    column_nodes = set()
    for mid, m in model.members.items():
        if m.type == "frame" and vertical_member(model, mid):
            column_nodes.update((m.i, m.j))
    ends = []
    for mid, m in model.members.items():
        if m.type != "frame" or vertical_member(model, mid):
            continue
        ends += [(mid, e) for e, n in (("i", m.i), ("j", m.j)) if n in column_nodes]
    out = []
    for dofs, label in ((("ry", "rz"), "both bending rotations"), (("rz",), "major-axis bending rotation")):
        todo = _open(ctx.doc, ends, dofs, True)
        if not todo:
            continue
        nodes = sorted({(model.members[m].i if e == "i" else model.members[m].j) for m, e in todo})
        out.append(Candidate(f"{label} released at {len(todo)} beam end(s) at the columns",
                             lambda d, todo=tuple(todo), dofs=dofs: _release(d, todo, dofs, True),
                             {"members": sorted({m for m, _ in todo}), "nodes": nodes}, {"releasedDofs": list(dofs)}))
    return out
