"""Tier E: equivalence checks for the modelling features ("curve balls").

Each check rewrites the INPUT DOCUMENT into a model that must give the same answer by a different
route through the engine, solves it, and compares on the original nodes and members:

designVariables   every design variable changed (K factors, unbraced lengths, Cb, ...): every analysis
                  result must be bitwise identical.  Always applied.
memberSplit       frame members cut in two at a result station (loads, releases, offsets carried to the
                  pieces): node results, end forces, and the station forces / displacement at the cut
                  must be reproduced.
rigidOffsets      every rigid end offset replaced by an added node and a stub member whose stiffness is
                  raised by a factor 10 per step; the extrapolated limit must equal the engine's result.
constraints       every rigid link / diaphragm replaced by stiff members (diaphragm links released out of
                  plane); extrapolated limit against the engine.
fullReleases      every full release replaced by a connection spring going to zero.
fixedEnds         on members that carry a release: every rigid bending connection replaced by a spring
                  going to infinity.

A series is adaptive (see ``_series``): the stiffness is changed by a factor 10 per level until the
answer stops changing in the expected way; the error of a stiff (soft) substitute is an expansion
a e + b e^2 + ... in e = 1/s (e = s), so two Richardson steps on the last three levels give the limit
(q^2 r3 - 2 q r2 + r1) / (q - 1)^2.  The observed order log10(step_{k-1} / step_k) must lie in
``POLICY["E"]["observedOrder"]`` on the last step unless the series has converged to round-off.

A check that applies but cannot be constructed or solved is NOT passed (the reason is recorded): the
feature then has no equivalence evidence.
"""
from __future__ import annotations

import copy
import math
import time

import numpy as np

from ..crosscheck.compare import relative_difference
from ..model import RIGID
from .common import (ONE_SIDED, TRUSS_TYPES, SolveFailure, arrays_identical, combine, compare_uniform,
                     first_order_counterpart, new_id, snapshot, solve, worst_of)
from .policy import POLICY
from .record import TierResult

__all__ = ["tier_e", "CHECKS", "split_document", "offsets_document", "constraints_document",
           "releases_document", "fixed_ends_document", "perturb_design_document"]

DOFS = ("ux", "uy", "uz", "rx", "ry", "rz")
AXIS = {"X": (1.0, 0.0, 0.0), "Y": (0.0, 1.0, 0.0), "Z": (0.0, 0.0, 1.0)}
STIFF_START = 1.0e1                 # first level of a stiff series (x10 per level)
SOFT_START = 1.0e-3                 # first level of a soft series (/10 per level)
MAX_LEVELS = 8
STEP_TARGET = 1.0e-6                # a series ends when a step is below this with the order in range
ROUND_OFF = 1.0e-10                 # a series whose steps are below this has converged to round-off
MAX_SPLIT = 40                      # members cut by the split check (largest models: a spread selection)


# ====================================================================================== helpers
def _flen(model, v, dim="length"):
    return float(model.units.from_internal(v, dim))


def _member_scale(model, mid):
    """(E, A, I, L, G, J) of a member in internal units, I = the larger second moment."""
    m = model.members[mid]
    mat, s = model.materials[m.material], model.sections[m.section]
    return float(mat.E), float(s.A), float(max(s.Iy, s.Iz)), float(m.length), float(mat.G), float(s.J)


def _reference_stiffness(model):
    """(E, kA, kI): a Young's modulus and the largest EA/L and EI/L of the members (internal units)."""
    E0, kA, kI = None, 0.0, 0.0
    for mid in model.members:
        E, A, I, L, _, _ = _member_scale(model, mid)
        E0 = E if E0 is None else E0
        kA, kI = max(kA, E * A / L), max(kI, E * I / L)
    return E0, kA, kI


def _weightless(doc, model, E):
    mat_id = "verify.weightless"
    if mat_id not in doc["materials"]:
        doc["materials"][mat_id] = {"E": _flen(model, E, "force/length^2"),
                                    "G": _flen(model, E / 2.6, "force/length^2")}
    return mat_id


def _section(doc, model, tag, A, I):
    sid = new_id(doc["sections"], f"{tag}.stiff")
    doc["sections"][sid] = {"type": "properties", "A": _flen(model, A, "length^2"),
                            "Iy": _flen(model, I, "length^4"), "Iz": _flen(model, I, "length^4"),
                            "J": _flen(model, 2.0 * I, "length^4")}
    return sid


def _stub_parts(doc, model, mid, length, s, tag):
    """Section and material of the stub that replaces a rigid end offset of member ``mid``: isotropic,
    weightless, EI/L and EA/L equal to ``s`` times those of the member."""
    E, A, I, L, _, _ = _member_scale(model, mid)
    return _section(doc, model, tag, s * A * length / L, s * I * length / L), _weightless(doc, model, E)


def _link_parts(doc, model, length, s, tag):
    """Section and material of a stiff member that replaces a rigid constraint: every end stiffness of
    the link (axial EA/L, transverse 3EI/L^3, rotational EI/L) is at least ``s`` times the largest
    axial (EA/L) resp. rotational (EI/L) member stiffness of the model."""
    E, kA, kI = _reference_stiffness(model)
    A = s * kA * length / E
    I = max(s * kA * length ** 3 / (3.0 * E), s * kI * length / E)
    return _section(doc, model, tag, A, I), _weightless(doc, model, E)


def _strip_design(doc):
    doc.pop("design", None)
    return doc


def _stiff_numerics(doc):
    """A substitute with very stiff members needs a lower pivot threshold than the default (which
    refuses a stiffness contrast beyond 1e10 as "almost a mechanism"); the threshold is a numerical
    setting of the derived model, not part of the structure."""
    an = doc.setdefault("analysis", {})
    an["pivotTolerance"] = 1.0e-15
    # round-off of an iteration with very stiff members sits near 1e-9 relative: an iteration tolerance
    # tighter than 1e-8 cannot be met by the substitute (it is 100 times below the tolerance of the check)
    if an.get("convergenceTolerance", 1.0e-6) < 1.0e-8:
        an["convergenceTolerance"] = 1.0e-8
    return doc


def _without_buckling(doc):
    """A stiffness series compares the combination results only; the linear buckling analysis of the
    substitute (not compared, and slow with very stiff members) is not run."""
    an = doc.get("analysis")
    if an and "buckling" in an:
        an = dict(an)
        an.pop("buckling")
        doc = {**doc, "analysis": an}
    return doc


def _same_axes(model, derived, pairs) -> str | None:
    for a, b in pairs:
        if not np.allclose(model.members[a].R, derived.members[b].R, atol=1e-12):
            return f"member {a!r}: the local axes of the derived member differ"
        if abs(model.members[a].length - sum(derived.members[x].length for x in ([b] if isinstance(b, str) else b))) \
                > 1e-9 * model.members[a].length and isinstance(b, str):
            return f"member {a!r}: the flexible length of the derived member differs"
    return None


# ====================================================================================== design variables
def perturb_design_document(raw: dict) -> dict:
    """Change every design variable that has no place in the analysis: K factors, unbraced lengths,
    Cb, Cm, method, bracing, slenderness limits -- in the defaults and in every group and member."""
    d = copy.deepcopy(raw)
    design = d.setdefault("design", {})
    defaults = design.setdefault("defaults", {})
    blocks = [defaults] + list(design.get("groups", {}).values()) + list(design.get("members", {}).values())
    for b in blocks:
        for name, f in (("Kx", 1.9), ("Ky", 0.55), ("Kz", 1.3)):
            b[name] = round(float(b.get(name, 1.0)) * f, 6) if isinstance(b.get(name, 1.0), (int, float)) else f
        for name in ("Lx", "Ly", "Lz", "Lb"):
            if isinstance(b.get(name), (int, float)) and not isinstance(b.get(name), bool):
                b[name] = b[name] * 0.61
        if isinstance(b.get("Cb"), (int, float)):
            b["Cb"] = 1.0 + 0.5 * float(b["Cb"])
        elif b is defaults:
            b["Cb"] = 1.37
        for name in ("Cmx", "Cmy"):
            if isinstance(b.get(name), (int, float)):
                b[name] = 0.85
    if "method" in defaults:
        defaults["method"] = "ASD" if defaults["method"] == "LRFD" else "LRFD"
    return d


def _results_identical(model, a, b) -> tuple[bool, str]:
    if set(a.combinations) != set(b.combinations) or set(a.buckling) != set(b.buckling):
        return False, "different set of combinations / buckling analyses"
    for cid, x in a.combinations.items():
        y = b.combinations[cid]
        for name in ("displacements", "reactions"):
            if not arrays_identical(getattr(x, name), getattr(y, name)):
                return False, f"combination {cid}: {name} differ"
        if (x.constraint_forces is None) != (y.constraint_forces is None) or (
                x.constraint_forces is not None and not arrays_identical(x.constraint_forces, y.constraint_forces)):
            return False, f"combination {cid}: constraint forces differ"
        if x.converged != y.converged or x.iterations != y.iterations or x.tau_b != y.tau_b:
            return False, f"combination {cid}: iteration record differs"
        for mid, mx in x.members.items():
            my = y.members[mid]
            for name in ("end_forces", "x", "forces", "displacements"):
                if not arrays_identical(getattr(mx, name), getattr(my, name)):
                    return False, f"combination {cid}, member {mid}: {name} differ"
            if mx.active != my.active:
                return False, f"combination {cid}, member {mid}: active flag differs"
    for cid, x in a.buckling.items():
        y = b.buckling[cid]
        if not (arrays_identical(x.factors, y.factors) and arrays_identical(x.modes, y.modes)
                and arrays_identical(x.negative_factors, y.negative_factors)):
            return False, f"buckling {cid}: factors or modes differ"
    return True, ""


def _check_design_variables(ctx) -> dict:
    out = {"name": "designVariables", "applies": True, "tolerance": 0.0,
           "what": "K factors, unbraced lengths, Cb, Cm, method changed in the design block; every analysis "
                   "result (displacements, reactions, member end forces, station forces and displacements, "
                   "iteration record, buckling factors and modes) must be bitwise identical"}
    try:
        _, res2 = solve(perturb_design_document(ctx["raw"]))
    except SolveFailure as e:
        return {**out, "passed": False, "reason": f"the perturbed model could not be analysed: {e.reason}"}
    same, why = _results_identical(ctx["model"], ctx["results"], res2)
    return {**out, "passed": same, "value": 0.0 if same else None, "reason": why}


# ====================================================================================== member split
def _load_span(model, load, L):
    """(x1, x2) of a distributed load / x of a point load along the flexible length, in FILE units."""
    rel = load.get("positionType", "relative") == "relative"
    f = L if rel else 1.0
    if load["type"] == "memberPoint":
        return load["position"] * f
    return load.get("start", 0.0) * f, (load["end"] * f if "end" in load else L)


def _split_candidates(raw, model):
    """{member: reason it is not split} and the list of members to split."""
    skipped, loaded = {}, set()
    self_weight = False
    for case in raw.get("loadCases", {}).values():
        for ld in case.get("loads", []):
            if ld["type"] in ("memberPoint", "memberDistributed"):
                loaded.add(ld["member"])
            elif ld["type"] == "selfWeight":
                self_weight = True
                loaded.update(ld.get("members", list(raw["members"])))
    per_member_rules = (any(c.notional is not None for c in model.combinations.values())
                        or (any(c.stiffness_reduction for c in model.combinations.values())
                            and model.analysis.tau_b == "computed"))
    stations = model.analysis.result_stations
    todo = []
    for mid, m in model.members.items():
        if m.type in TRUSS_TYPES:
            skipped[mid] = "truss-type member: an added joint on a pin-ended member is a mechanism"
        elif per_member_rules and mid in loaded:
            skipped[mid] = ("the model uses notional loads or a computed tau_b, which the engine defines per "
                            "analysis member (docs/THEORY.md 14): cutting a LOADED member is a different model")
        elif stations < 3:
            skipped[mid] = "fewer than three result stations: no interior station to cut at"
        else:
            todo.append(mid)
    if len(todo) > MAX_SPLIT:
        keep = {todo[int(round(k))] for k in np.linspace(0, len(todo) - 1, MAX_SPLIT)}
        for mid in todo:
            if mid not in keep:
                skipped[mid] = f"not selected (the split check cuts at most {MAX_SPLIT} members)"
        todo = [m for m in todo if m in keep]
    return todo, skipped


def split_document(raw: dict, model, members=None):
    """Cut the given frame members in two at an interior result station.

    Returns ``(document, cuts, skipped)`` with ``cuts[member] = {"pieces": (id1, id2), "node": id,
    "x": position in internal units, "station": index of the original result station}``."""
    doc = _strip_design(copy.deepcopy(raw))
    todo, skipped = _split_candidates(raw, model)
    if members is not None:
        todo = [m for m in todo if m in members]
    n_st = model.analysis.result_stations
    point_x = {}
    for case in raw.get("loadCases", {}).values():
        for ld in case.get("loads", []):
            if ld["type"] == "memberPoint" and ld["member"] in model.members:
                L = _flen(model, model.members[ld["member"]].length)
                point_x.setdefault(ld["member"], []).append(_load_span(model, ld, L))
    cuts = {}
    for mid in todo:
        m = model.members[mid]
        L = _flen(model, m.length)
        order = sorted(range(1, n_st - 1), key=lambda k: abs(k / (n_st - 1) - 0.37))
        k = next((k for k in order
                  if all(abs(L * k / (n_st - 1) - x) > 1e-9 * L for x in point_x.get(mid, []))), None)
        if k is None:
            skipped[mid] = "every interior result station carries a point load"
            continue
        frac = k / (n_st - 1)
        xs = L * frac
        nid = new_id(doc["nodes"], f"{mid}.cut")
        p = (np.asarray(m.start) + frac * (np.asarray(m.end) - np.asarray(m.start)))
        doc["nodes"][nid] = {"x": _flen(model, p[0]), "y": _flen(model, p[1]), "z": _flen(model, p[2])}
        src = doc["members"].pop(mid)
        a, b = new_id(doc["members"], f"{mid}.s1"), new_id(list(doc["members"]) + [f"{mid}.s1"], f"{mid}.s2")
        pa = {k2: copy.deepcopy(v) for k2, v in src.items() if k2 not in ("releases", "offsets")}
        pb = copy.deepcopy(pa)
        pa["j"], pb["i"] = nid, nid
        rel, off = src.get("releases", {}), src.get("offsets", {})
        if "i" in rel:
            pa["releases"] = {"i": rel["i"]}
        if "j" in rel:
            pb["releases"] = {"j": rel["j"]}
        # offsets as global vectors: "local" offsets refer to the axes of the member they belong to
        if np.any(m.offset_i):
            pa["offsets"] = {"axes": "global", "i": [_flen(model, v) for v in m.offset_i]}
        if np.any(m.offset_j):
            pb["offsets"] = {"axes": "global", "j": [_flen(model, v) for v in m.offset_j]}
        doc["members"][a], doc["members"][b] = pa, pb
        cuts[mid] = {"pieces": (a, b), "node": nid, "x": float(m.length) * frac, "xFile": xs, "L": L, "station": k}
    # loads
    for case in doc.get("loadCases", {}).values():
        new_loads = []
        for ld in case.get("loads", []):
            mid = ld.get("member")
            if ld["type"] == "selfWeight" and "members" in ld:
                ld = dict(ld)
                ld["members"] = [x for mm in ld["members"] for x in (cuts[mm]["pieces"] if mm in cuts else (mm,))]
                new_loads.append(ld)
            elif ld["type"] == "memberPoint" and mid in cuts:
                c = cuts[mid]
                x = _load_span(model, ld, c["L"])
                xs = c["xFile"]
                piece, xp = (c["pieces"][0], x / xs) if x < xs else (c["pieces"][1], (x - xs) / (c["L"] - xs))
                new_loads.append({**ld, "member": piece, "position": min(max(xp, 0.0), 1.0),
                                  "positionType": "relative"})
            elif ld["type"] == "memberDistributed" and mid in cuts:
                c = cuts[mid]
                x1, x2 = _load_span(model, ld, c["L"])
                w1, w2 = ld["w1"], ld.get("w2", ld["w1"])
                w = lambda x: w1 + (w2 - w1) * ((x - x1) / (x2 - x1) if x2 > x1 else 0.0)
                xs = c["xFile"]
                if x1 < xs and x2 > x1:
                    e = min(x2, xs)
                    new_loads.append({**ld, "member": c["pieces"][0], "start": x1 / xs, "end": min(e / xs, 1.0),
                                      "w1": w1, "w2": w(e), "positionType": "relative"})
                if x2 > xs:
                    s = max(x1, xs)
                    rest = c["L"] - xs
                    new_loads.append({**ld, "member": c["pieces"][1], "start": max((s - xs) / rest, 0.0),
                                      "end": min((x2 - xs) / rest, 1.0), "w1": w(s), "w2": w2,
                                      "positionType": "relative"})
            else:
                new_loads.append(ld)
        case["loads"] = new_loads
    sr = doc.get("analysis", {}).get("stiffnessReduction", {})
    if "members" in sr:
        sr["members"] = [x for mm in sr["members"] for x in (cuts[mm]["pieces"] if mm in cuts else (mm,))]
    return doc, cuts, skipped


def _check_split(ctx) -> dict:
    raw, model, results = ctx["raw"], ctx["model"], ctx["results"]
    out = {"name": "memberSplit", "what": "frame members cut in two at a result station"}
    doc, cuts, skipped = split_document(raw, model)
    if not cuts:
        return {**out, "applies": False, "passed": None, "reason": "no member can be cut",
                "skipped": skipped}
    nonlinear = any(c.analysis_type != "firstOrder" or c.stiffness_reduction or c.notional is not None
                    for c in model.combinations.values()) or any(m.type in ONE_SIDED for m in model.members.values())
    tol = POLICY["E"]["splitNonlinear" if nonlinear else "splitFirstOrder"]
    out.update(applies=True, tolerance=tol, membersCut=len(cuts), skipped=skipped)
    try:
        m2, r2 = solve(doc)
    except SolveFailure as e:
        return {**out, "passed": False, "reason": f"the split model could not be analysed: {e.reason}"}
    for mid, c in cuts.items():
        a, b = c["pieces"]
        if not (np.allclose(m2.members[a].R, model.members[mid].R, atol=1e-12)
                and np.allclose(m2.members[b].R, model.members[mid].R, atol=1e-12)):
            return {**out, "passed": False, "reason": f"member {mid!r}: the local axes of the pieces differ"}
    n = len(model.nodes)
    worst, per_combo, ok = None, {}, True
    for cid, cres in results.combinations.items():
        if cid not in r2.combinations:
            return {**out, "passed": False, "reason": f"combination {cid!r} has no result in the split model"}
        s2 = snapshot(m2, r2, cid, n)
        for mid, c in cuts.items():
            a, b = c["pieces"]
            s2.end_forces[mid] = np.r_[s2.end_forces[a][:6], s2.end_forces[b][6:]]
        s1 = snapshot(model, results, cid)
        cmp = compare_uniform(s1, s2, model, tol)
        w = worst_of(cmp)
        # at the cut: station forces and displacement of the original member against the pieces
        fscale = max(cmp.scales["endForce"], 1e-300)
        mscale = max(cmp.scales["endMoment"], 1e-300)
        # translations: the largest of the node translations and the member-axis translations (a beam
        # between two supports has no node translation at all)
        dscale = max(cmp.scales["displacement"],
                     max((float(np.abs(mr.displacements).max(initial=0.0)) for mr in cres.members.values()),
                         default=0.0), 1e-300)
        cut_worst = {"value": 0.0}
        for mid, c in cuts.items():
            mr = cres.members[mid]
            k = int(np.argmin(np.abs(mr.x - c["x"])))
            a, b = c["pieces"]
            f_cut = np.asarray(r2.combinations[cid].members[a].end_forces[6:], float)
            diff = np.abs(mr.forces[k] - f_cut) / np.array([fscale] * 3 + [mscale] * 3)
            u_node = np.asarray(r2.combinations[cid].displacements[m2.nodes[c["node"]].index][:3], float)
            du = np.abs(np.asarray(model.members[mid].R).T @ mr.displacements[k] - u_node) / dscale
            for val, what in ((float(diff.max()), f"member {mid} station forces at the cut"),
                              (float(du.max()), f"member {mid} displacement at the cut")):
                if val > cut_worst["value"]:
                    cut_worst = {"value": val, "where": what, "quantity": "station"}
        if cut_worst["value"] > w["value"]:
            w = cut_worst
        per_combo[cid] = w
        ok = ok and cmp.passed and cut_worst["value"] <= tol
        if worst is None or w["value"] > worst["value"]:
            worst = {**w, "combination": cid}
    for cid, b in results.buckling.items():
        f2 = r2.buckling.get(cid)
        if f2 is None or len(f2.factors) < len(b.factors):
            return {**out, "passed": False, "reason": f"buckling of {cid!r} has no result in the split model"}
        nmodes = len(b.factors)
        rel = float(np.max(np.abs(np.asarray(f2.factors[:nmodes]) - b.factors) / np.abs(b.factors))) if nmodes else 0.0
        per_combo[f"buckling:{cid}"] = {"value": rel, "where": "buckling factors", "quantity": "bucklingFactor"}
        ok = ok and rel <= tol
        if worst is None or rel > worst["value"]:
            worst = {"value": rel, "where": "buckling factors", "quantity": "bucklingFactor", "combination": cid}
    return {**out, "passed": bool(ok), "value": worst["value"] if worst else 0.0, "worst": worst,
            "perCombination": per_combo}


# ====================================================================================== stiffness series
def _series(ctx, name, build, start, direction, base=None):
    """Adaptive stiffness series: solve ``build(s)`` for s = start, start x q, start x q^2, ... (q = 10
    for ``direction = +1``, 1/10 for -1) until the answer stops changing in the expected way, then
    compare the extrapolated limit with the engine.  The stopping rule does not look at the engine's
    result:

    * step_k = largest relative difference between the results of levels k-1 and k (all combinations);
      order_k = log_q(step_{k-1} / step_k), accepted range ``POLICY["E"]["observedOrder"]``;
    * a level whose order leaves the range after a level inside it (or after a step already below
      ``STEP_TARGET``) marks the round-off floor: it is discarded and the series ends; so does a level
      that the engine refuses (an invariant fails by round-off) once two levels exist;
    * the series also ends when a step falls below ``STEP_TARGET`` with the order in range, when a step
      is below ``ROUND_OFF``, or after ``MAX_LEVELS`` levels.

    The limit uses the last three levels (two Richardson steps).  ``base`` = (model, results) to compare
    against when the check runs on a counterpart of the model."""
    model, results = base or (ctx["model"], ctx["results"])
    tol = POLICY["E"]["limit"]
    ratio = POLICY["E"]["seriesRatio"]
    lo, hi = POLICY["E"]["observedOrder"]
    out = {"applies": True, "tolerance": tol}
    n = len(model.nodes)
    cids = list(results.combinations)
    runs, levels, steps, orders, stop = [], [], [None], [None], None
    in_range = lambda o: o is not None and lo <= o <= hi
    s = start
    for _ in range(MAX_LEVELS):
        try:
            m2, r2 = solve(_without_buckling(build(s)))
            missing = [c for c in cids if c not in r2.combinations]
            if missing:
                raise SolveFailure(f"combination {missing[0]!r} has no result in the substitute")
        except SolveFailure as e:
            if len(runs) >= 2:
                stop = f"level {s:g} could not be analysed ({e.reason[:160]}); the series ends at {levels[-1]:g}"
                break
            return {**out, "levels": levels + [s], "passed": False,
                    "reason": f"the substitute model at level {s:g} could not be analysed: {e.reason}"}
        snap = {cid: snapshot(m2, r2, cid, n) for cid in cids}
        if runs:
            step = max(relative_difference(runs[-1][cid], snap[cid], model) for cid in cids)
            order = math.log(steps[-1] / step, ratio) if (steps[-1] and step > 0.0) else None
            if order is not None and not in_range(order) and (in_range(orders[-1]) or steps[-1] <= STEP_TARGET):
                stop = (f"level {s:g} left the asymptotic range (observed order {order:.2f}): round-off floor; "
                        f"discarded, the series ends at {levels[-1]:g}")
                break
            steps.append(step), orders.append(order)
        runs.append(snap), levels.append(s)
        if len(runs) > 1 and (steps[-1] <= ROUND_OFF or (steps[-1] <= STEP_TARGET and (in_range(orders[-1])
                                                                                   or len(runs) == 2))):
            stop = (f"step {steps[-1]:.1e} below {STEP_TARGET:g}"
                    + (" with the observed order in range" if in_range(orders[-1]) else ""))
            break
        s = s * ratio if direction > 0 else s / ratio
    else:
        stop = f"{MAX_LEVELS} levels"
    out.update(levels=levels, stepDifferences=steps[1:], observedOrders=orders[2:], stopped=stop)
    if len(runs) < 2:
        return {**out, "passed": False, "reason": "fewer than two levels could be analysed"}
    # the limit: two Richardson steps on the last three levels when the last step converges at the
    # expected order; the last level itself when the series has stopped changing (last step <= STEP_TARGET)
    richardson = len(runs) >= 3 and in_range(orders[-1])
    order_ok = richardson or steps[-1] <= STEP_TARGET
    out["limit"] = ("Richardson extrapolation of the last three levels" if richardson
                    else "the last level (the series has stopped changing)")
    q = ratio - 1.0
    ok, worst, per_combo = True, None, {}
    for cid in cids:
        eng = snapshot(model, results, cid)
        if richardson:
            r1, r2, r3 = (r[cid] for r in runs[-3:])
            # error = a e + b e^2 + ...: two Richardson steps remove both terms
            limit = combine([(ratio * ratio / (q * q), r3), (-2.0 * ratio / (q * q), r2), (1.0 / (q * q), r1)],
                            "limit")
        else:
            limit = runs[-1][cid]
        cmp = compare_uniform(eng, limit, model, tol)
        w = worst_of(cmp)
        per_combo[cid] = {**w, "differenceToEngineByLevel": [relative_difference(eng, r[cid], model) for r in runs]}
        ok = ok and cmp.passed
        if worst is None or w["value"] > worst["value"]:
            worst = {**w, "combination": cid}
    res = {**out, "passed": bool(ok and order_ok), "value": worst["value"] if worst else 0.0, "worst": worst,
           "perCombination": per_combo}
    if not order_ok:
        res["reason"] = (f"the series neither converges at the expected order (observed {orders[2:]}, accepted "
                         f"[{lo}, {hi}]) nor stops changing (last step {steps[-1]:.1e}); no limit can be stated")
    return res


# ---------------------------------------------------------------------------------- rigid offsets
def offsets_document(raw: dict, model, s: float) -> dict:
    """Every rigid end offset replaced by an added node at the flexible end and a stub member from the
    joint to it, with EI/L and EA/L equal to ``s`` times those of the member."""
    doc = _strip_design(copy.deepcopy(raw))
    for mid, m in model.members.items():
        if not m.has_offsets:
            continue
        d = doc["members"][mid]
        for end, off, node, pos in (("i", m.offset_i, m.i, m.start), ("j", m.offset_j, m.j, m.end)):
            if not np.any(off):
                continue
            nid = new_id(doc["nodes"], f"{mid}.o{end}")
            doc["nodes"][nid] = {"x": _flen(model, pos[0]), "y": _flen(model, pos[1]), "z": _flen(model, pos[2])}
            sid, mat = _stub_parts(doc, model, mid, float(np.linalg.norm(off)), s, f"{mid}.o{end}")
            stub = new_id(doc["members"], f"{mid}.o{end}")
            doc["members"][stub] = {"i": node, "j": nid, "section": sid, "material": mat}
            d[end] = nid
        d.pop("offsets", None)
    return _stiff_numerics(doc)


def _check_offsets(ctx) -> dict:
    model = ctx["model"]
    out = {"name": "rigidOffsets", "what": "rigid end offsets replaced by an added node and a stub member, "
                                           "stiffness x10 per level; extrapolated limit against the engine"}
    if not any(m.has_offsets for m in model.members.values()):
        return {**out, "applies": False, "passed": None, "reason": "no rigid offsets"}
    res = _series(ctx, "rigidOffsets", lambda s: offsets_document(ctx["raw"], model, s), STIFF_START, +1)
    return {**out, **res, "members": [m for m, v in model.members.items() if v.has_offsets]}


# ---------------------------------------------------------------------------------- constraints
def constraints_document(raw: dict, model, s: float):
    """Every rigid link replaced by a stiff frame member master -> slave; every diaphragm link by a stiff
    member with local y along the normal whose slave end releases the out-of-plane shear, the torsion and
    the bending about the in-plane axis (so that only the two in-plane forces and the moment about the
    normal are transmitted).  A master without members gets its out-of-plane DOFs fixed (the engine
    reports them as zero).  Raises ``SolveFailure`` when a link cannot be built."""
    doc = _strip_design(copy.deepcopy(raw))
    attached = {nid: 0 for nid in model.nodes}
    for m in model.members.values():
        attached[m.i] += 1
        attached[m.j] += 1
    for cid, c in model.constraints.items():
        xm = np.asarray(model.nodes[c.master].xyz, float)
        for sl in c.slaves:
            r = np.asarray(model.nodes[sl].xyz, float) - xm
            length = float(np.linalg.norm(r))
            if length <= 1e-9 * max(1.0, float(np.abs(xm).max())):
                raise SolveFailure(f"constraint {cid!r}: slave {sl!r} is at the master's position; a stiff member "
                                   "of zero length does not exist")
            sid, mat = _link_parts(doc, model, length, s, f"{cid}.{sl}")
            link = {"i": c.master, "j": sl, "section": sid, "material": mat}
            if c.type == "rigidDiaphragm":
                nrm = np.array(AXIS[c.normal])
                if abs(float(r @ nrm)) > 1e-9 * length:
                    raise SolveFailure(f"constraint {cid!r}: slave {sl!r} is not in the master's plane; the "
                                       "diaphragm rule is then not a rigid-body relation (docs/THEORY.md 4)")
                link["orientation"] = {"referenceVector": list(AXIS[c.normal])}
                link["releases"] = {"j": {"uy": True, "rx": True, "rz": True}}
            doc["members"][new_id(doc["members"], f"{cid}.{sl}")] = link
        if c.type == "rigidDiaphragm" and attached[c.master] == 0:
            k = "XYZ".index(c.normal)
            free = [DOFS[k]] + [DOFS[3 + a] for a in range(3) if a != k]
            sup = doc.setdefault("supports", {}).setdefault(c.master, {})
            dofs = sup.setdefault("dofs", {})
            fixed_by_type = {"fixed": DOFS, "pinned": DOFS[:3]}.get(sup.get("type"), ())
            for name in free:
                if name not in dofs and name not in fixed_by_type:
                    dofs[name] = "fixed"
    doc.pop("constraints", None)
    return _stiff_numerics(doc)


def _check_constraints(ctx) -> dict:
    model, raw = ctx["model"], ctx["raw"]
    out = {"name": "constraints", "what": "rigid links / diaphragms replaced by stiff members, stiffness x10 "
                                          "per level; extrapolated limit against the engine"}
    if not model.constraints:
        return {**out, "applies": False, "passed": None, "reason": "no constraints"}
    base, note = None, None
    second = any(c.analysis_type == "secondOrder" for c in model.combinations.values()) \
        or bool(model.analysis.buckling_combinations)
    if second:
        # the engine's constraints carry no geometric stiffness (docs/THEORY.md 16) while a stiff member
        # under axial force does, and the released diaphragm link is refused in a second-order analysis:
        # the kinematic relation is checked on the first-order counterpart of the same combinations
        raw = first_order_counterpart(raw)
        note = ("checked on the first-order counterpart of the combinations: the engine's rigid constraints have "
                "no geometric stiffness (docs/THEORY.md 16), a stiff member has; the second-order results of this "
                "model are therefore covered by tiers C and D only")
        try:
            base = solve(raw)
        except SolveFailure as e:
            return {**out, "applies": True, "passed": False, "tolerance": POLICY["E"]["limit"],
                    "reason": f"the first-order counterpart could not be analysed: {e.reason}"}
    cmodel = base[0] if base else model

    def build(s):
        return constraints_document(raw, cmodel, s)
    res = _series(ctx, "constraints", build, STIFF_START, +1, base)
    if note:
        res["note"] = note
    return {**out, **res, "constraints": list(model.constraints)}


# ---------------------------------------------------------------------------------- releases
def _spring_value(model, mid, dof, s):
    """Connection stiffness of member end DOF ``dof`` = ``s`` x the member's own stiffness of that kind,
    in FILE units."""
    E, A, I, L, G, J = _member_scale(model, mid)
    if dof == "ux":
        return _flen(model, s * E * A / L, "force/length")
    if dof in ("uy", "uz"):
        return _flen(model, s * E * I / L ** 3, "force/length")
    if dof == "rx":
        return _flen(model, s * G * J / L, "moment/rotation")
    return _flen(model, s * E * I / L, "moment/rotation")


def releases_document(raw: dict, model, s: float) -> dict:
    """Every full release (``true``) of a frame member replaced by a connection spring of ``s`` times
    the member's own stiffness."""
    doc = copy.deepcopy(raw)
    for mid, d in doc["members"].items():
        if model.members[mid].type in TRUSS_TYPES:
            continue
        for end in ("i", "j"):
            for dof, v in d.get("releases", {}).get(end, {}).items():
                if v is True or (isinstance(v, dict) and v.get("spring") == 0):
                    d["releases"][end][dof] = {"spring": _spring_value(model, mid, dof, s)}
    return doc


def fixed_ends_document(raw: dict, model, s: float) -> dict:
    """On every frame member that carries a release: each rigid bending connection (ry, rz at both
    ends) replaced by a connection spring of ``s`` times the member's own stiffness."""
    doc = copy.deepcopy(raw)
    for mid, d in doc["members"].items():
        m = model.members[mid]
        if m.type in TRUSS_TYPES or not m.has_releases:
            continue
        rel = d.setdefault("releases", {})
        for end, k_end in (("i", m.k_end_i), ("j", m.k_end_j)):
            for a in (4, 5):
                if k_end[a] == RIGID:
                    rel.setdefault(end, {})[DOFS[a]] = {"spring": _spring_value(model, mid, DOFS[a], s)}
    return doc


def _check_full_releases(ctx) -> dict:
    model = ctx["model"]
    out = {"name": "fullReleases", "what": "full releases replaced by connection springs going to zero; "
                                           "extrapolated limit against the engine"}
    members = [mid for mid, m in model.members.items()
               if m.type not in TRUSS_TYPES and any(k == 0.0 for k in m.k_end_i + m.k_end_j)]
    if not members:
        return {**out, "applies": False, "passed": None, "reason": "no full release on a frame member"}
    res = _series(ctx, "fullReleases", lambda s: releases_document(ctx["raw"], model, s), SOFT_START, -1)
    return {**out, **res, "members": members}


def _check_fixed_ends(ctx) -> dict:
    model = ctx["model"]
    out = {"name": "fixedEnds", "what": "rigid bending connections of released members replaced by connection "
                                        "springs going to infinity; extrapolated limit against the engine"}
    members = [mid for mid, m in model.members.items()
               if m.type not in TRUSS_TYPES and m.has_releases
               and any(k[a] == RIGID for k in (m.k_end_i, m.k_end_j) for a in (4, 5))]
    if not members:
        return {**out, "applies": False, "passed": None,
                "reason": "no released frame member with a rigid bending connection"}
    res = _series(ctx, "fixedEnds", lambda s: fixed_ends_document(ctx["raw"], model, s), STIFF_START, +1)
    return {**out, **res, "members": members}


CHECKS = {"designVariables": _check_design_variables, "memberSplit": _check_split,
          "rigidOffsets": _check_offsets, "constraints": _check_constraints,
          "fullReleases": _check_full_releases, "fixedEnds": _check_fixed_ends}


def tier_e(raw: dict, model, results, deadline: float | None = None) -> TierResult:
    t0 = time.perf_counter()
    ctx = {"raw": raw, "model": model, "results": results}
    checks, largest = {}, []
    for name, fn in CHECKS.items():
        if deadline is not None and time.perf_counter() > deadline:
            checks[name] = {"name": name, "applies": True, "passed": False,
                            "reason": "not run: the time budget of verify_model was used up"}
            continue
        t1 = time.perf_counter()
        try:
            checks[name] = fn(ctx)
        except SolveFailure as e:
            checks[name] = {"name": name, "applies": True, "passed": False,
                            "reason": f"the equivalent model cannot be constructed: {e.reason}"}
        checks[name]["seconds"] = round(time.perf_counter() - t1, 4)
    applied = [c for c in checks.values() if c.get("applies")]
    failed = [c["name"] for c in applied if c.get("passed") is not True]
    for c in applied:
        if c.get("value") is not None:
            largest.append({"what": c["name"], "value": c["value"], "tolerance": c.get("tolerance"),
                            **({"where": c["worst"].get("where"), "combination": c["worst"].get("combination")}
                               if c.get("worst") else {})})
    largest.sort(key=lambda d: -(d["value"] or 0.0))
    passed = not failed
    summary = (f"{len(applied)} check(s) applied ({', '.join(c['name'] for c in applied)}); "
               + ("all passed" if passed else "NOT passed: " + "; ".join(
                   f"{n} ({checks[n].get('reason') or 'difference above tolerance'})" for n in failed)))
    details = {"checks": {k: {kk: vv for kk, vv in v.items() if kk != "seconds"} for k, v in checks.items()}}
    tr = TierResult("E", True, passed, summary, largest, dict(POLICY["E"]), details, time.perf_counter() - t0)
    tr.timings = {k: v.get("seconds") for k, v in checks.items()}
    return tr
