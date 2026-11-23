"""Measured difficulty of a model: ``difficulty(document) -> {properties..., "score", "level"}``.

Every quantity is measured on the model itself (docs/DIFFICULTY.md gives the definitions and the
reasons); nothing is taken from the generator's intentions.  The level is a function of six
measured components:

=============== ================================================================ ============
component       measured property                                                points
=============== ================================================================ ============
size            free DOFs (not fixed, not slaves, not inactive; in-plane DOFs    0..4
                only for a planar model): <=12, <=30, <=72, <=150, more
indeterminacy   degree of static indeterminacy: 0, 1-3, 4-12, 13-48, more        0..4
dimension       three-dimensional (not planar)                                   0..1
features        kinds present among: member end releases; rigid end offsets;    0..3
                rigid constraints; tension-/compression-only members; support
                springs or imposed support displacements (at most 3 counted)
secondOrder     no second-order combination 0; some 1; and the largest          0..2
                amplification (second-order / first-order largest
                translation) >= 1.10: 2
combinations    number of combinations: <=2, 3-5, more                           0..2
=============== ================================================================ ============

``level`` = 1 for a score of 0-2, 2 for 3-5, 3 for 6-8, 4 for 9-11, 5 for 12 and more.
"""
from __future__ import annotations

import copy
import math

import numpy as np

from ..model import build_model

__all__ = ["difficulty", "feature_flags", "static_properties", "amplification", "level_of", "score_components",
           "SIZE_LIMITS", "INDETERMINACY_LIMITS", "LEVEL_LIMITS", "AMPLIFICATION_BAND"]

SIZE_LIMITS = (12, 30, 72, 150)            # free DOFs: points = number of limits exceeded
INDETERMINACY_LIMITS = (0, 3, 12, 48)
LEVEL_LIMITS = (2, 5, 8, 11)               # score: level = 1 + number of limits exceeded
AMPLIFICATION_BAND = 1.10
FEATURE_CAP = 3
_TOL = 1e-9


def _points(value, limits) -> int:
    return sum(1 for lim in limits if value > lim)


def level_of(score: int) -> int:
    return 1 + _points(score, LEVEL_LIMITS)


# ------------------------------------------------------------------------------------- planarity
def _plane(model):
    """Unit normal (a global axis index) of the plane of a planar model, else None.

    Planar means: all nodes lie in one plane normal to a global axis; every member has one
    principal axis along that normal; every load acts in the plane (forces in the plane, moments
    about the normal); rigid end offsets lie in the plane.  A straight line of nodes is taken to lie
    in the plane that also contains the loads if one global-axis plane fits, vertical planes first.
    """
    xyz = np.array([n.xyz for n in model.nodes.values()], float)
    size = float(np.linalg.norm(np.ptp(xyz, axis=0))) or 1.0
    vertical = "XYZ".index(model.vertical_axis)
    candidates = [k for k in range(3) if np.ptp(xyz[:, k]) <= 1e-9 * size]
    candidates.sort(key=lambda k: k == vertical)            # prefer a plane that contains the vertical axis
    for k in candidates:
        if _in_plane(model, k):
            return k
    return None


def _in_plane(model, k: int) -> bool:
    n = np.zeros(3)
    n[k] = 1.0

    def force_ok(v):
        return abs(float(np.dot(v, n))) <= _TOL * (float(np.linalg.norm(v)) or 1.0)

    def moment_ok(v):
        v = np.asarray(v, float)
        return float(np.linalg.norm(v - np.dot(v, n) * n)) <= _TOL * (float(np.linalg.norm(v)) or 1.0)

    for m in model.members.values():
        if not (abs(abs(m.R[1] @ n) - 1.0) < 1e-9 or abs(abs(m.R[2] @ n) - 1.0) < 1e-9):
            return False
        if not (force_ok(m.offset_i) and force_ok(m.offset_j)):
            return False
    for case in model.load_cases.values():
        for load in case.loads:
            kind = type(load).__name__
            if kind == "NodalLoad":
                if not (force_ok(load.values[:3]) and moment_ok(load.values[3:])):
                    return False
            elif kind == "MemberPointLoad":
                g = model.members[load.member].R.T @ load.vector
                if not (force_ok(g) if load.kind == "force" else moment_ok(g)):
                    return False
            elif kind == "MemberDistributedLoad":
                R = model.members[load.member].R
                if not (force_ok(R.T @ load.w1) and force_ok(R.T @ load.w2)):
                    return False
            elif kind == "SelfWeight":
                if not force_ok(load.direction):
                    return False
            elif kind == "SupportDisplacement":
                for dof, value in load.values.items():
                    if value and (dof == k or (dof >= 3 and dof - 3 != k)):
                        return False
    return True


# ------------------------------------------------------------------------------------- counts
def static_properties(model) -> dict:
    """Counts that need no analysis: members, nodes, free DOFs, degree of static indeterminacy,
    planarity.

    Degree of static indeterminacy = unknown forces - equilibrium equations, with

    * unknowns: 6 per frame member (3 in a planar model) minus its fully released end DOFs
      (a spring-connected end is not a release), 1 per truss-type member; one per support DOF that
      is fixed or sprung; one per DOF of a slave node that a rigid constraint eliminates;
    * equations: 6 per node (3 planar), but only 3 (2 planar) at a node joined only by truss-type
      members and by no constraint, and only the DOFs that drive slaves at a memberless master.

    The count equals the dimension of the self-stress space when the model is stable (which the
    engine establishes); it is a property of the model as posed, so a tension-only member counts as
    a member.
    """
    k = _plane(model)
    planar = k is not None
    if planar:
        trans = [d for d in range(3) if d != k]
        dofs = trans + [3 + k]                              # two in-plane translations, rotation about the normal
    else:
        trans, dofs = [0, 1, 2], list(range(6))
    names = ("ux", "uy", "uz", "rx", "ry", "rz")
    frame_nodes, truss_nodes = set(), set()
    unknowns = 0
    for m in model.members.values():
        if m.type == "frame":
            frame_nodes.update((m.i, m.j))
            if planar:
                n = np.zeros(3)
                n[k] = 1.0
                local = [0] + [a for a in (1, 2) if abs(m.R[a] @ n) < 1e-6] + [3 + a for a in (1, 2) if abs(abs(m.R[a] @ n) - 1) < 1e-6]
            else:
                local = list(range(6))
            released = sum(1 for end in (m.k_end_i, m.k_end_j) for a in local if end[a] == 0.0)
            unknowns += len(local) - released
        else:
            truss_nodes.update((m.i, m.j))
            unknowns += 1
    slave_dofs, master_dofs = {}, {}
    for c in model.constraints.values():
        cd = [names.index(d) for d in c.dofs if names.index(d) in dofs]
        for s in c.slaves:
            slave_dofs.setdefault(s, set()).update(cd)
        master_dofs.setdefault(c.master, set()).update(cd)
    equations = free = 0
    active = {}
    for nid in model.nodes:
        if nid in frame_nodes or nid in slave_dofs:
            act = set(dofs)
        elif nid in truss_nodes:
            act = set(trans)
        else:
            act = set()
        act |= master_dofs.get(nid, set())
        active[nid] = act
        equations += len(act)
    reactions = 0
    for nid, act in active.items():
        s = model.supports.get(nid)
        sl = slave_dofs.get(nid, set())
        for d in act:
            if d in sl:
                unknowns += 1                               # constraint force
                continue
            if s is not None and (s.fixed[d] or s.spring[d] != 0.0):
                reactions += 1
            if s is None or not s.fixed[d]:
                free += 1
    degree = unknowns + reactions - equations
    return {"nNodes": len(model.nodes), "nMembers": len(model.members), "planar": planar,
            "plane": None if not planar else "normal " + "XYZ"[k], "freeDofs": free,
            "staticIndeterminacy": degree, "nReactions": reactions}


def feature_flags(document: dict, model=None) -> dict:
    """Which modelling features a document uses (plain booleans; no analysis)."""
    model = model or build_model(document)
    rel = {"full": False, "fullOther": False, "rot": False, "other": False}
    for m in model.members.values():
        if m.type != "frame":
            continue
        for end in (m.k_end_i, m.k_end_j):
            for a, kk in enumerate(end):
                if kk == 0.0:
                    rel["full"] = True
                    if a < 4:
                        rel["fullOther"] = True
                elif not math.isinf(kk):
                    rel["rot" if a >= 4 else "other"] = True
    loads = {"NodalLoad": False, "pointForce": False, "pointMoment": False, "uniform": False, "partial": False,
             "varying": False, "SelfWeight": False, "SupportDisplacement": False}
    for case in model.load_cases.values():
        for load in case.loads:
            kind = type(load).__name__
            if kind == "MemberPointLoad":
                loads["pointForce" if load.kind == "force" else "pointMoment"] = True
            elif kind == "MemberDistributedLoad":
                L = model.members[load.member].length
                full = load.x1 <= 1e-9 * L and load.x2 >= L * (1 - 1e-9)
                if not np.allclose(load.w1, load.w2, rtol=1e-12, atol=0.0):
                    loads["varying"] = True
                elif full:
                    loads["uniform"] = True
                else:
                    loads["partial"] = True
            else:
                loads[kind] = True
    projected = any(l.get("projected") for c in document.get("loadCases", {}).values() for l in c.get("loads", []))
    types = {m.type for m in model.members.values()}
    combos = model.combinations.values()
    return {
        "secondOrder": any(c.analysis_type == "secondOrder" for c in combos),
        "buckling": bool(model.analysis.buckling_combinations),
        "shearDeformation": bool(model.analysis.shear_deformation),
        "trussMembers": "truss" in types,
        "tensionOnly": "tensionOnly" in types,
        "compressionOnly": "compressionOnly" in types,
        "releasesFull": rel["full"],
        "releasesFullNonBending": rel["fullOther"],
        "releasesPartialRotational": rel["rot"],
        "releasesPartialOther": rel["other"],
        "rigidOffsets": any(np.any(m.offset_i != 0) or np.any(m.offset_j != 0) for m in model.members.values()),
        "rigidLink": any(c.type == "rigidLink" for c in model.constraints.values()),
        "rigidDiaphragm": any(c.type == "rigidDiaphragm" for c in model.constraints.values()),
        "supportSprings": any(any(k != 0.0 for k in s.spring) for s in model.supports.values()),
        "prescribedDisplacements": loads["SupportDisplacement"],
        "loadNodal": loads["NodalLoad"],
        "loadMemberPointForce": loads["pointForce"],
        "loadMemberPointMoment": loads["pointMoment"],
        "loadMemberUniform": loads["uniform"],
        "loadMemberPartialUniform": loads["partial"],
        "loadMemberLinearlyVarying": loads["varying"],
        "loadProjected": bool(projected),
        "selfWeight": loads["SelfWeight"],
        "stiffnessReduction": any(c.stiffness_reduction for c in combos),
        "notionalLoads": any(c.notional is not None for c in combos),
        "designBlock": bool(document.get("design")),
    }


def _feature_kinds(flags: dict) -> list:
    kinds = []
    if flags["releasesFull"] or flags["releasesPartialRotational"] or flags["releasesPartialOther"]:
        kinds.append("releases")
    if flags["rigidOffsets"]:
        kinds.append("rigidOffsets")
    if flags["rigidLink"] or flags["rigidDiaphragm"]:
        kinds.append("constraints")
    if flags["tensionOnly"] or flags["compressionOnly"]:
        kinds.append("oneSidedMembers")
    if flags["supportSprings"] or flags["prescribedDisplacements"]:
        kinds.append("springsOrSettlements")
    return kinds


# ------------------------------------------------------------------------------------- amplification
def amplification(document: dict, results=None, first_order=None) -> dict:
    """Second-order amplification of every second-order combination: largest nodal translation
    of the combination as analysed / the same with the combination switched to first order.
    Returns {combination id: factor}; empty when no combination is second order.
    ``results``: results of ``document`` if already computed; ``first_order``: results of the
    same model with every combination first order, if already computed."""
    from ..solve import analyze
    so = [cid for cid, c in document["combinations"].items()
          if c.get("analysisType", document.get("analysis", {}).get("type", "firstOrder")) == "secondOrder"]
    if not so:
        return {}
    if results is None:
        results = analyze(build_model(document))
    lin = first_order
    if lin is None:
        first = copy.deepcopy(document)
        first["combinations"] = {cid: {**document["combinations"][cid], "analysisType": "firstOrder"} for cid in so}
        first.setdefault("analysis", {}).pop("buckling", None)
        first.pop("design", None)
        lin = analyze(build_model(first))
    out = {}
    for cid in so:
        u2 = float(np.max(np.linalg.norm(results.combinations[cid].displacements[:, :3], axis=1)))
        u1 = float(np.max(np.linalg.norm(lin.combinations[cid].displacements[:, :3], axis=1)))
        out[cid] = u2 / u1 if u1 > 0 else 1.0
    return out


# ------------------------------------------------------------------------------------- level
def score_components(static: dict, kinds, n_combinations: int, n_second_order: int, amp: float | None) -> dict:
    second = 0
    if n_second_order:
        second = 2 if (amp is not None and amp >= AMPLIFICATION_BAND) else 1
    return {"size": _points(static["freeDofs"], SIZE_LIMITS),
            "indeterminacy": _points(static["staticIndeterminacy"], INDETERMINACY_LIMITS),
            "dimension": 0 if static["planar"] else 1,
            "features": min(len(kinds), FEATURE_CAP),
            "secondOrder": second,
            "combinations": 0 if n_combinations <= 2 else 1 if n_combinations <= 5 else 2}


def difficulty(document: dict, results=None, analyse: bool = True, first_order=None) -> dict:
    """Measured properties, score and level of a model document.

    ``analyse=False`` skips the two analyses behind the amplification factor; a model with
    second-order combinations then gets ``amplification = None`` and ``level = None`` together with
    ``levelRange`` (the two levels the amplification could give)."""
    model = build_model(document)
    static = static_properties(model)
    flags = feature_flags(document, model)
    kinds = _feature_kinds(flags)
    n_so = sum(1 for c in model.combinations.values() if c.analysis_type == "secondOrder")
    amps = amplification(document, results, first_order) if (analyse and n_so) else {}
    amp = max(amps.values()) if amps else None
    comp = score_components(static, kinds, len(model.combinations), n_so, amp)
    score = sum(comp.values())
    out = {**static, "threeDimensional": not static["planar"], "featureKinds": kinds,
           "nCombinations": len(model.combinations), "nSecondOrderCombinations": n_so,
           "amplification": None if amp is None else round(amp, 4),
           "amplificationByCombination": {k: round(v, 4) for k, v in amps.items()},
           "components": comp, "score": score}
    if n_so and amp is None:
        out["level"] = None
        out["levelRange"] = [level_of(score), level_of(score + 1)]
    else:
        out["level"] = level_of(score)
        out["levelRange"] = [out["level"], out["level"]]
    return out
