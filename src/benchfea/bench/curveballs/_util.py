"""Helpers shared by the curve-ball modules (document-level edits and model queries)."""
from __future__ import annotations

import numpy as np

from ...geometry import AXIS, is_vertical
from ..families._common import tidy

DOFS = ("ux", "uy", "uz", "rx", "ry", "rz")
TRANSLATIONS, ROTATIONS = DOFS[:3], DOFS[3:]
AXIS_OF = {"X": 0, "Y": 1, "Z": 2}

__all__ = ["DOFS", "TRANSLATIONS", "ROTATIONS", "tidy", "support_dofs", "support_object", "set_support",
           "plane_dofs", "frame_ends_at", "vertical_member", "member_axis", "new_id", "release_state",
           "set_release", "level_groups", "section_depth", "is_number"]


def is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def support_dofs(support: dict, vertical: str) -> dict:
    """DOF -> "fixed" | "free" | {"spring": k}, shorthand expanded (the rule of the validator)."""
    v = {"X": "ux", "Y": "uy", "Z": "uz"}[vertical]
    base = {"fixed": DOFS, "pinned": DOFS[:3], "roller": (v,)}.get(support.get("type"), ())
    out = {d: ("fixed" if d in base else "free") for d in DOFS}
    out.update(support.get("dofs") or {})
    return out


def support_object(dofs: dict):
    """The support entry for a DOF table (None: nothing restrained): shorthand where one fits."""
    rest = {d: dofs[d] for d in DOFS if dofs.get(d, "free") != "free"}
    if not rest:
        return None
    fixed = {d for d, v in rest.items() if v == "fixed"}
    if fixed == set(DOFS):
        return {"type": "fixed"}
    if set(TRANSLATIONS) <= fixed:
        extra = {d: v for d, v in rest.items() if d not in TRANSLATIONS}
        return {"type": "pinned", **({"dofs": extra} if extra else {})}
    return {"dofs": rest}


def set_support(doc: dict, node: str, dofs: dict) -> None:
    obj = support_object(dofs)
    if obj is None:
        doc["supports"].pop(node, None)
    else:
        doc["supports"][node] = obj


def plane_dofs(doc: dict) -> set:
    """DOFs fixed at EVERY node of the model: the out-of-plane restraint of a plane frame modelled
    in 3D space (empty for a space model).  Support curve balls never touch these."""
    vertical = doc["conventions"]["verticalAxis"]
    supports = doc.get("supports", {})
    common = set(DOFS)
    for nid in doc["nodes"]:
        s = supports.get(nid)
        if s is None:
            return set()
        common &= {d for d, v in support_dofs(s, vertical).items() if v == "fixed"}
    return common


def frame_ends_at(model, node: str) -> list:
    """[(member id, "i" | "j")] of the frame members that end at ``node``."""
    out = []
    for mid, m in model.members.items():
        if m.type != "frame":
            continue
        if m.i == node:
            out.append((mid, "i"))
        if m.j == node:
            out.append((mid, "j"))
    return out


def member_axis(model, mid: str) -> np.ndarray:
    return np.asarray(model.members[mid].R[0], dtype=float)


def vertical_member(model, mid: str) -> bool:
    return is_vertical(member_axis(model, mid), model.vertical_axis)


def new_id(existing, base: str) -> str:
    cand, k = base[:56], 1
    while cand in existing:
        k += 1
        cand = f"{base[:56]}-{k}"
    return cand


def release_state(member: dict, end: str, dof: str):
    """False (rigid) | True (released) | {"spring": k} as written in the document."""
    return (member.get("releases") or {}).get(end, {}).get(dof, False)


def set_release(member: dict, end: str, dof: str, value) -> None:
    member.setdefault("releases", {}).setdefault(end, {})[dof] = value


def level_groups(model, tol: float = 1e-6) -> dict:
    """elevation (mm, rounded) -> node ids at that elevation, lowest first."""
    k = AXIS_OF[model.vertical_axis]
    size = max(float(np.ptp([n.xyz[k] for n in model.nodes.values()])), 1.0)
    groups: dict = {}
    for nid, n in model.nodes.items():
        z = float(n.xyz[k])
        for key in groups:
            if abs(key - z) <= tol * size:
                groups[key].append(nid)
                break
        else:
            groups[z] = [nid]
    return dict(sorted(groups.items()))


def section_depth(section: dict):
    """Overall depth (along local y) of a parametric section definition, file units; None if unknown."""
    return {"I": section.get("d"), "box": section.get("H"), "pipe": section.get("D"),
            "rectangle": section.get("h")}.get(section.get("type"))


def unit(v) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    return v / n if n > 0 else v


_ = AXIS  # re-exported for the modules
