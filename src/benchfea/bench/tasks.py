"""Benchmark tasks: verified cases and the four task types built on them (docs/BENCHMARK.md).

A CASE is one curve-balled model with its verified reference:

    case = build_case("portal_frame", "support_fixed_to_pinned", seed=3, difficulty=2)
    case["status"]            # "verified" | "unavailable"  (with case["skips"]: why attempts were dropped)

A TASK is a case seen through a task type and a mode; it comes in two parts:

    task, answers = make_task(case, "predict", "tool-allowed")
    task       # PUBLIC: statement, the model file, the answer format, quantities / criteria
    answers    # HIDDEN: change record, reference values, verification record -- never shown before scoring

Task types (``TASK_TYPES``): predict | diagnose | repair | design.  Modes (``MODES``): tool-allowed
(the agent may run the analysis) | no-solver (the run-analysis tool is disabled).

Only VERIFIED references are served.  ``build_case`` puts the curve-balled model through
``benchfea.verify.verify_model``; an attempt whose reference is not ``verified`` (or whose model
cannot be generated, or to which the curve ball does not apply) is skipped, the skip is logged, and
the next seed of the derivation chain (``derived_seed``) is tried, at most ``MAX_ATTEMPTS`` times.
"verified" means only: passed the tiers listed in the record, at the tolerances stated in it.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path

from .. import __version__
from ..model import build_model
from ..problems import ModelInvalid
from . import curveballs
from .families._common import CATALOGUE, UNIT_SYSTEMS, Draw, section_definition, stable_seed
from .generate import GENERATOR_VERSION, GenerationFailed, generate
from .splits import split_of

__all__ = ["TASK_TYPES", "MODES", "TASK_VERSION", "MAX_ATTEMPTS", "TaskUnavailable", "derived_seed", "build_case",
           "make_task", "task_types_for", "describe", "extract", "select_quantities", "change_direction",
           "reference_design", "design_document", "convert_value", "quantity_dimension", "model_weight", "design_feasible", "public_verification",
           "verification_summary", "case_key", "QUANTITY_KINDS", "DIRECTIONS", "CATEGORY_TEXT", "MECHANISM_KINDS",
           "MECHANISM_MOTIONS", "REPAIR_TOLERANCE", "catalogue_for", "task_id"]

TASK_TYPES = ("predict", "diagnose", "repair", "design")
MODES = ("tool-allowed", "no-solver")
TASK_VERSION = "1"             # bump when a change alters any task or answer for the same case
MAX_ATTEMPTS = 3               # seeds of the derivation chain tried before a case is unavailable
QUANTITY_KINDS = ("displacement", "reaction", "memberEndForce", "memberMaxMoment", "bucklingFactor", "designRatio")
DIRECTIONS = ("added", "removed", "increased", "decreased", "reversed", "changed")
MECHANISM_KINDS = ("node", "rigid-body", "distributed")
MECHANISM_MOTIONS = ("translation", "rotation", "translation+rotation")
CATEGORY_TEXT = dict(curveballs.CATEGORIES)
SIGNIFICANT = 0.05             # a quantity is asked for only if it is at least this fraction of the largest of its kind
MOVED = 0.05                   # ... and preferably one the curve ball moved by at least this (relative)
REPAIR_TOLERANCE = {"full": 0.02, "zero": 0.10}
VERIFY_BUDGET = 300.0          # seconds per verification run of a case (what is not run is not passed)
DESIGN_MAX_PASSES = 4

_DISP = ("ux", "uy", "uz", "rx", "ry", "rz")
_FORCE = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")


class TaskUnavailable(Exception):
    """No verified task for the request; ``reason`` is one line, ``log`` the skipped attempts."""

    def __init__(self, reason: str, log: list | None = None):
        super().__init__(reason)
        self.reason, self.log = reason, list(log or [])


# ============================================================================================ cases
def derived_seed(family: str, seed: int, difficulty: int, attempt: int) -> int:
    """Generator / curve-ball seed of attempt ``attempt`` of the instance ``(family, seed, difficulty)``.
    Attempt 0 uses ``seed`` itself; attempt k > 0 uses the first 31 bits of
    SHA-256("benchfea.bench|task-chain|<TASK_VERSION>|family|seed|difficulty|k")."""
    if attempt == 0:
        return int(seed)
    return stable_seed("task-chain", TASK_VERSION, family, int(seed), int(difficulty), int(attempt)) >> 97


def case_key(family: str, curve_ball: str, seed: int, difficulty: int) -> dict:
    return {"family": family, "curveBall": curve_ball, "seed": int(seed), "difficulty": int(difficulty),
            "versions": {"task": TASK_VERSION, "generator": GENERATOR_VERSION,
                         "curveBall": curveballs.CURVEBALL_VERSION, "engine": __version__}}


def _key_hash(key: dict) -> str:
    return hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()


def task_id(key: dict, task_type: str, mode: str) -> str:
    """Opaque id of a task (it must not name the curve ball)."""
    return "t" + hashlib.sha256(json.dumps([key, task_type, mode], sort_keys=True).encode()).hexdigest()[:14]


def _read_cache(path: Path):
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError, EOFError):
        return None


def _write_cache(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    with gzip.open(tmp, "wt", encoding="utf-8") as fh:
        json.dump(data, fh)
    os.replace(tmp, path)


def _base_model(family: str, gseed: int, difficulty: int, cache_dir, foreign_cache) -> dict:
    """Generated base model of one attempt with its verification -> {"status": "ok" | "failed", ...}.
    The model is renamed (the generator's name states the seed) BEFORE anything is verified."""
    path = None
    if cache_dir:
        path = Path(cache_dir) / f"base-{family}-L{difficulty}-g{gseed}-{GENERATOR_VERSION}-{TASK_VERSION}.json.gz"
        hit = _read_cache(path)
        if hit is not None:
            return hit
    try:
        g = generate(family, gseed, difficulty)
    except GenerationFailed as e:
        out = {"status": "failed", "category": "generation-failed",
               "reason": f"no acceptable model: {e.reasons[-1] if e.reasons else e}"}
    except Exception as e:                                   # a defect of the generator: logged, next seed
        out = {"status": "failed", "category": "generator-error",
               "reason": f"the generator raised {type(e).__name__}: {str(e)[:300]}"}
    else:
        doc = curveballs.record.plain_copy(g.document)
        doc["name"] = f"{family.replace('_', ' ')} model"
        v, dv = _verify(doc, foreign_cache)
        out = {"status": "ok", "document": doc, "units": g.metadata["units"], "verification": v,
               "designVerification": dv}
    if path:
        _write_cache(path, out)
    return out


def _verify(doc: dict, foreign_cache) -> tuple:
    """(analysis record, design record or None) of a model, both as dicts.

    The analysis record is ``verify_model(doc)``: the reference is the results document (or, for a
    mechanism, the diagnosis).  When it is ``verified`` and of kind "results", the model is verified
    once more WITH its design results (``design=True``: tier F applies); that second record says
    whether the design results may be served as a reference too."""
    from ..verify import verify_model
    v = verify_model(doc, budget=VERIFY_BUDGET, cache=foreign_cache).to_dict(volatile=False)
    dv = None
    if v["status"] == "verified" and v["kind"] == "results":
        dv = verify_model(doc, design=True, budget=VERIFY_BUDGET, cache=foreign_cache).to_dict(volatile=False)
        dv.pop("reference", None)                                # the same results document as in ``v``
    return v, dv


def _verified(v: dict | None, kind: str = "results") -> bool:
    return bool(v) and v["status"] == "verified" and v["kind"] == kind


def _quarantine_text(v: dict) -> str:
    q = v.get("quarantine") or {}
    if v["status"] == "refused":
        return "refused: " + "; ".join(p.get("message", "") for p in (v.get("refusal") or {}).get("problems", [])[:2])
    return f"{v['status']} ({q.get('category')}): {str(q.get('diagnosis'))[:400]}"


def build_case(family: str, curve_ball: str, seed: int, difficulty: int, *, cache_dir=None,
               foreign_cache=None) -> dict:
    """One curve-balled model with its verified reference (see the module docstring).

    Returns a plain dict (JSON-serialisable): ``status`` "verified" with ``base`` (document, verification),
    ``perturbed`` (document), ``record`` (the change record), ``effectClass``, ``verification`` (the full
    record of the perturbed model, reference included), or ``status`` "unavailable"; ``skips`` lists every
    attempt that was dropped: ``{"attempt", "generatorSeed", "stage", "category", "reason"}``.
    Deterministic; cached on disk when ``cache_dir`` is given.  Never raises for a model problem."""
    curveballs.get(curve_ball)                                   # KeyError for an unknown name
    key = case_key(family, curve_ball, seed, difficulty)
    path = None
    if cache_dir:
        path = Path(cache_dir) / f"case-{family}-L{difficulty}-s{int(seed)}-{curve_ball}-{_key_hash(key)[:10]}.json.gz"
        hit = _read_cache(path)
        if hit is not None:
            return hit
    skips = []
    case = {"key": key, "status": "unavailable", "skips": skips}

    def skip(k, gseed, stage, category, reason):
        skips.append({"attempt": k, "generatorSeed": gseed, "stage": stage, "category": category,
                      "reason": str(reason)[:600]})

    for k in range(MAX_ATTEMPTS):
        gseed = derived_seed(family, seed, difficulty, k)
        base = _base_model(family, gseed, difficulty, cache_dir, foreign_cache)
        if base["status"] != "ok":
            skip(k, gseed, "generate", base.get("category", "generation-failed"), base["reason"])
            continue
        try:
            pert = curveballs.apply(base["document"], curve_ball, gseed)
        except curveballs.NotApplicable as e:
            skip(k, gseed, "curve-ball", "not-applicable", e.reason)
            continue
        cls = pert.record.effect_class
        if cls == "invalid-model":
            skip(k, gseed, "verify", "refused", "the curve-balled model is refused by validation: there is no "
                                                "analysis and nothing the verification pipeline can verify")
            continue
        v, dv = _verify(pert.document, foreign_cache)
        if cls == "creates-mechanism":
            bad = None if _verified(v, "diagnosis") else v
        elif cls == "design-only":                               # the design results ARE the reference
            bad = v if not _verified(v) else (None if _verified(dv) else dv)
        else:
            bad = None if _verified(v) else v
        if bad is not None:
            skip(k, gseed, "verify", (bad.get("quarantine") or {}).get("category") or bad["status"],
                 ("with the design results (tier F): " if bad is dv else "") + _quarantine_text(bad))
            continue
        case.update(status="verified", attempt=k, generatorSeed=gseed, effectClass=cls, units=base["units"],
                    base={"document": base["document"], "verification": base["verification"],
                          "designVerification": base["designVerification"]},
                    perturbed={"document": pert.document}, record=pert.record.to_dict(), verification=v,
                    designVerification=dv)
        break
    if path:
        _write_cache(path, case)
    return case


# ===================================================================================== verification
def verification_summary(v: dict) -> dict:
    """The record of a verification run without the reference documents and the per-check details:
    status, tiers passed / not passed, and per tier the summary, the tolerances and the largest
    differences found."""
    return {"format": v["format"], "kind": v["kind"], "status": v["status"], "meaning": v["meaning"],
            "modelHash": v["modelHash"], "tiersPassed": v["tiersPassed"], "tiersNotPassed": v["tiersNotPassed"],
            "tiersNotApplicable": v["tiersNotApplicable"], "quarantine": v.get("quarantine"),
            "tiers": {k: {"applies": t["applies"], "passed": t["passed"], "summary": t["summary"],
                          "tolerances": t["tolerances"], "largestDifferences": t["largestDifferences"][:3]}
                      for k, t in v["tiers"].items()},
            "versions": v["versions"], "notes": v.get("notes", [])}


def public_verification(v: dict, scope: str) -> dict:
    """What a task says about its reference BEFORE scoring: that it is verified and what the word means --
    not the kind of reference or the tiers (they would tell a mechanism from a stable model)."""
    return {"status": "verified", "scope": scope, "meaning": v["meaning"],
            "policy": "docs/VERIFICATION.md; the tiers passed and the tolerances are in the answers file "
                      "and in the info returned when the episode ends"}


# ======================================================================================= quantities
def extract(q: dict, results: dict | None, design: dict | None = None) -> float:
    """Value of a named quantity (file units) from a results document (``benchfea-results``) or, for
    ``designRatio``, a design document.  Raises ``KeyError`` when the entity does not exist."""
    kind = q["kind"]
    if kind == "designRatio":
        g = design["members"][q["designMember"]]["governing"]
        if g is None or g.get("ratio") is None:
            raise KeyError("no governing ratio")
        return float(g["ratio"])
    if kind == "bucklingFactor":
        return float(results["buckling"][q["combination"]]["factors"][int(q.get("mode", 1)) - 1])
    c = results["combinations"][q["combination"]]
    if kind == "displacement":
        return float(c["displacements"][q["node"]][q["component"]])
    if kind == "reaction":
        return float(c["reactions"][q["node"]][q["component"]])
    if kind == "memberEndForce":
        return float(c["members"][q["member"]]["endForces"][q["end"]][q["component"]])
    if kind == "memberMaxMoment":
        return max(abs(float(x)) for x in c["members"][q["member"]]["stations"][q["component"]])
    raise KeyError(f"unknown quantity kind {kind!r}")


def convert_value(value: float, dimension: str, from_units: dict, to_units: dict) -> float:
    """A number of ``dimension`` (length | rotation | force | moment | dimensionless) given in the unit
    system ``from_units`` of one file, expressed in the unit system ``to_units`` of another."""
    if dimension == "dimensionless" or from_units == to_units:
        return value
    from ..units import UnitSystem
    a, b = UnitSystem.from_dict(from_units), UnitSystem.from_dict(to_units)
    return float(b.from_internal(a.to_internal(value, dimension), dimension))


def quantity_dimension(q: dict) -> str:
    if q["kind"] in ("bucklingFactor", "designRatio"):
        return "dimensionless"
    return "moment" if q["kind"] == "memberMaxMoment" else _dim(q["component"])


def _dim(component: str) -> str:
    return {"u": "length", "r": "rotation", "F": "force", "M": "moment"}[component[0]]


def _candidates(results: dict) -> list:
    out = []
    for cid, c in results["combinations"].items():
        for nid, row in c["displacements"].items():
            out += [{"kind": "displacement", "combination": cid, "node": nid, "component": k} for k in _DISP]
        for nid in c["reactions"]:
            out += [{"kind": "reaction", "combination": cid, "node": nid, "component": k} for k in _FORCE]
        for mid, m in c["members"].items():
            if not m.get("active", True):
                continue
            for end in ("i", "j"):
                out += [{"kind": "memberEndForce", "combination": cid, "member": mid, "end": end, "component": k}
                        for k in ("Fx", "My", "Mz")]
            out += [{"kind": "memberMaxMoment", "combination": cid, "member": mid, "component": k}
                    for k in ("My", "Mz")]
    return out


def _entity(q: dict):
    return q["kind"], q.get("node") or q.get("member") or q.get("designMember") or q.get("combination")


def select_quantities(reference: dict, other: dict | None, draw: Draw, plan=None, total: int = 6,
                      moved_by: float = MOVED, allow_unmoved: bool = True) -> list:
    """Named quantities of ``reference`` (a results document).

    Candidates: nodal displacements, reactions, member end forces (Fx, My, Mz) and the largest bending
    moment along a member, of every combination, that are SIGNIFICANT -- at least ``SIGNIFICANT`` of
    the largest value of the same dimension in their combination.  With ``other`` (the results of the
    model before / after the curve ball) only entities both models have are used, and a candidate has
    MOVED when the two values differ physically (after unit conversion) by at least ``moved_by``
    (relative to the larger of the two).
    Choice (seeded, one quantity per entity): per kind of ``plan`` among the moved candidates; then,
    up to ``total``, among the moved candidates of any kind; only if fewer than three quantities were
    found that way -- and ``allow_unmoved`` -- among the candidates that did not move.
    Each entry: the quantity, ``value``, ``unit``, ``dimension`` and ``floor`` (the absolute error
    floor: 0.001 of the scale)."""
    plan = plan or (("displacement", 2), ("reaction", 2), ("memberEndForce", 1), ("memberMaxMoment", 1))
    labels = reference["unitLabels"]
    cands, scale = [], {}
    for q in _candidates(reference):
        try:
            v = extract(q, reference)
        except (KeyError, IndexError, TypeError, ValueError):
            continue
        dim = quantity_dimension(q)
        key = (q["combination"], dim)
        scale[key] = max(scale.get(key, 0.0), abs(v))
        cands.append((q, v, dim))
    pool, moved = [], []
    for q, v, dim in cands:
        s = scale[(q["combination"], dim)]
        if s <= 0.0 or abs(v) < SIGNIFICANT * s:
            continue
        entry = dict(q, value=v, unit=labels[dim], dimension=dim, floor=1.0e-3 * s)
        if other is None:
            pool.append(entry)
            continue
        try:
            w = convert_value(extract(q, other), dim, other["units"], reference["units"])
        except (KeyError, IndexError, TypeError, ValueError):
            continue                                             # not an entity of both models
        (moved if abs(v - w) >= moved_by * max(abs(v), abs(w)) else pool).append(entry)

    def shuffled(group):
        group = list(group)
        for i in range(len(group) - 1, 0, -1):                   # Fisher-Yates on the stable draw
            j = draw.integer(0, i)
            group[i], group[j] = group[j], group[i]
        return group

    chosen, seen = [], set()

    def take(group, count):
        for e in group:
            if count <= 0 or len(chosen) >= total:
                break
            if _entity(e) in seen:
                continue
            seen.add(_entity(e))
            chosen.append(e)
            count -= 1

    first = moved if other is not None else pool
    for kind, count in plan:
        take(shuffled([e for e in first if e["kind"] == kind]), count)
    take(shuffled(first), total)
    if len(chosen) < 3 and allow_unmoved:
        take(shuffled(pool), total)
    if reference.get("buckling"):
        cid = next(iter(reference["buckling"]))
        if reference["buckling"][cid]["factors"]:
            v = float(reference["buckling"][cid]["factors"][0])
            chosen.append({"kind": "bucklingFactor", "combination": cid, "mode": 1, "value": v, "unit": "-",
                           "dimension": "dimensionless", "floor": 1.0e-3 * abs(v)})
    for n, e in enumerate(chosen, 1):
        e["id"] = f"q{n}"
    return chosen


def _has(results: dict, q: dict) -> bool:
    try:
        extract(q, results)
        return True
    except (KeyError, IndexError, TypeError, ValueError):
        return False


def _design_quantity(design: dict, other: dict | None, n: int) -> dict | None:
    """The design member whose governing ratio the curve ball moved most (the largest ratio otherwise)."""
    best = None
    for mid, m in design["members"].items():
        g = m.get("governing")
        if not g or g.get("ratio") is None:
            continue
        r = float(g["ratio"])
        move = 0.0
        if other is not None:
            g0 = (other["members"].get(mid) or {}).get("governing")
            if g0 and g0.get("ratio") is not None:
                move = abs(r - float(g0["ratio"])) / max(abs(r), abs(float(g0["ratio"])), 1e-12)
        if best is None or (move, r) > best[0]:
            best = ((move, r), mid, r)
    if best is None:
        return None
    return {"kind": "designRatio", "designMember": best[1], "value": best[2], "unit": "-",
            "dimension": "dimensionless", "floor": 1.0e-3 * max(abs(best[2]), 1e-12), "id": f"q{n}"}


def _quantity_text(q: dict) -> str:
    k = q["kind"]
    if k == "displacement":
        what = "translation" if q["component"][0] == "u" else "rotation"
        return f"{what} {q['component']} of node {q['node']} in combination {q['combination']} [{q['unit']}], global axes"
    if k == "reaction":
        what = "reaction force" if q["component"][0] == "F" else "reaction moment"
        return (f"{what} {q['component']} at node {q['node']} in combination {q['combination']} [{q['unit']}], global "
                "axes, acting on the structure")
    if k == "memberEndForce":
        return (f"end force {q['component']} at end {q['end']} of member {q['member']} in combination "
                f"{q['combination']} [{q['unit']}], member local axes, acting on the member (docs/CONVENTIONS.md)")
    if k == "memberMaxMoment":
        return (f"largest absolute bending moment {q['component']} along member {q['member']} in combination "
                f"{q['combination']} [{q['unit']}] (a non-negative number; local axes)")
    if k == "bucklingFactor":
        return f"lowest positive elastic buckling factor of combination {q['combination']} [-]"
    return f"governing design ratio (demand / capacity, AISC 360-22) of design member {q['designMember']} [-]"


def _public_quantity(q: dict) -> dict:
    return {k: v for k, v in q.items() if k not in ("value", "floor", "dimension")} | {"text": _quantity_text(q)}


# ====================================================================================== description
_LABELS = {"i": "start node", "j": "end node", "E": "E", "G": "G", "Fy": "Fy", "Fu": "Fu", "w1": "intensity at start",
           "w2": "intensity at end", "unitWeight": "unit weight", "dofs": "degrees of freedom", "Kx": "Kx", "Ky": "Ky",
           "Kz": "Kz", "Lb": "Lb", "Cb": "Cb", "analysisType": "analysis type", "notionalLoad": "notional load",
           "stiffnessReduction": "stiffness reduction", "shearDeformation": "shear deformation",
           "deflectionLimits": "deflection limits", "spanRatio": "span / ", "rollAngle": "roll angle",
           "topFlange": "top flange braced at", "bottomFlange": "bottom flange braced at",
           "minor": "minor-axis buckling braced at", "major": "major-axis buckling braced at",
           "torsion": "torsional buckling braced at", "tensionOnly": "tension-only", "compressionOnly": "compression-only"}


def _num(v) -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, (int, float)):
        return json.dumps(v)
    return str(v)


def _words(key: str) -> str:
    if key in _LABELS:
        return _LABELS[key]
    out = ""
    for ch in key:
        out += (" " + ch.lower()) if (ch.isupper() and out and not out[-1].isupper() and len(key) > 2) else ch
    return out


def _say(v, top: bool = True) -> str:
    if isinstance(v, dict):
        if not v:
            return "none"
        text = ", ".join(f"{_words(k)} {_say(x, False)}" for k, x in v.items())
        return text if top else f"({text})"
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_say(x, False) for x in v) + "]"
    return _num(v)


def describe(document: dict) -> str:
    """The model as an engineer would describe it in words: every value of the document, section by
    section, in the units of the document.  Used for the INTENDED design of a diagnose / repair task
    (built from the model before the curve ball; the file itself is not given)."""
    d = document
    u = d.get("units", {})
    L = []
    L.append("Units: " + _say(u) + " (every number below is in these units: moments in force x length, line loads in "
             "force / length, moduli and stresses in force / length^2, unit weights in force / length^3, rotational "
             "springs in force x length / rad).")
    for k, v in d.items():
        if k in ("units", "nodes", "supports", "materials", "sections", "members", "constraints", "loadCases",
                 "combinations", "analysis", "design", "name"):
            continue
        L.append(f"{_words(k).capitalize()}: {_say(v)}.")
    L.append(f"Nodes ({len(d.get('nodes', {}))}), coordinates:")
    for nid, n in d.get("nodes", {}).items():
        rest = {k: v for k, v in n.items() if k not in ("x", "y", "z")}
        L.append(f"  {nid} at ({_num(n.get('x', 0))}, {_num(n.get('y', 0))}, {_num(n.get('z', 0))})"
                 + (f"; {_say(rest)}" if rest else ""))
    L.append("Supports (a degree of freedom not listed is free):")
    for nid, s in d.get("supports", {}).items():
        L.append(f"  node {nid}: {_say(s)}")
    L.append("Materials:")
    for mid, m in d.get("materials", {}).items():
        L.append(f"  {mid}: {_say(m)}")
    L.append("Cross-sections:")
    for sid, s in d.get("sections", {}).items():
        L.append(f"  {sid}: {_say(s)}")
    L.append(f"Members ({len(d.get('members', {}))}); a member end without a stated release is rigidly connected, a "
             "member without stated offsets runs from node to node:")
    for mid, m in d.get("members", {}).items():
        rest = {k: v for k, v in m.items() if k not in ("i", "j", "section", "material")}
        L.append(f"  {mid}: from node {m.get('i')} to node {m.get('j')}, section {m.get('section')}, material "
                 f"{m.get('material')}" + (f"; {_say(rest)}" if rest else ""))
    if d.get("constraints"):
        L.append("Rigid links and diaphragms:")
        for cid, c in d["constraints"].items():
            L.append(f"  {cid}: {_say(c)}")
    else:
        L.append("Rigid links and diaphragms: none.")
    L.append("Load cases:")
    for cid, c in d.get("loadCases", {}).items():
        rest = {k: v for k, v in c.items() if k not in ("loads", "description")}
        L.append(f"  {cid}" + (f" ({c['description']})" if c.get("description") else "") + (f", {_say(rest)}" if rest else "")
                 + ":")
        for n, load in enumerate(c.get("loads", []), 1):
            rest = {k: v for k, v in load.items() if k != "type"}
            L.append(f"    load {n}: {_words(str(load.get('type')))}" + (f", {_say(rest)}" if rest else ""))
    L.append("Load combinations:")
    for cid, c in d.get("combinations", {}).items():
        rest = {k: v for k, v in c.items() if k != "factors"}
        terms = " + ".join(f"{_num(f)} x {lc}" for lc, f in c.get("factors", {}).items())
        L.append(f"  {cid} = {terms}" + (f"; {_say(rest)}" if rest else ""))
    L.append("Analysis settings: " + (_say(d["analysis"]) if d.get("analysis") else "defaults") + ".")
    if d.get("design"):
        L.append("Design data (AISC 360-22 member checks):")
        for k, v in d["design"].items():
            if isinstance(v, dict) and v and all(isinstance(x, dict) for x in v.values()):
                L.append(f"  {_words(k)}:")
                for kk, vv in v.items():
                    L.append(f"    {kk}: {_say(vv)}")
            else:
                L.append(f"  {_words(k)}: {_say(v)}")
    return "\n".join(L)


# ========================================================================================= diagnose
def change_direction(changes: list) -> str:
    """Direction of a change record, derived mechanically from its edits (docs/BENCHMARK.md):
    ``added`` every edit adds something; ``removed`` every edit removes something; ``reversed`` every
    edit replaces a number by its negative; ``increased`` / ``decreased`` every edit replaces a number
    by one of larger / smaller magnitude; ``changed`` anything else (a substituted name, a mixture)."""
    kinds = set()
    for c in changes:
        op = c["op"]
        if op == "add":
            kinds.add("added")
        elif op == "remove":
            kinds.add("removed")
        else:
            b, a = c.get("before"), c.get("after")
            num = all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in (b, a))
            if not num:
                kinds.add("changed")
            elif a == -b and a != 0:
                kinds.add("reversed")
            elif abs(a) > abs(b):
                kinds.add("increased")
            elif abs(a) < abs(b):
                kinds.add("decreased")
            else:
                kinds.add("changed")
    return kinds.pop() if len(kinds) == 1 else "changed"


# =========================================================================================== design
def catalogue_for(document: dict, unit_system: str) -> dict:
    """Sizing groups and the sections each may take.

    A group is the set of members that share a section in the given model; a group whose section is
    an I-shape may take any I-shape of the generic catalogue (``bench.CATALOGUE``: beam and column
    series), a box any box, a pipe any pipe; any other group is fixed.  Returns ``{"groups": {section
    id: [members]}, "options": {section id: [names, lightest first]}, "sections": {name: definition}}``."""
    from ..sections import from_definition
    metric = "metric" if UNIT_SYSTEMS[unit_system]["metric"] else "imperial"
    groups, options, sections = {}, {}, {}
    for mid, m in document["members"].items():
        groups.setdefault(m["section"], []).append(mid)
    for sid in groups:
        kind = document["sections"][sid].get("type")
        series = {"I": ("beam", "column"), "box": ("box",), "pipe": ("pipe",)}.get(kind, ())
        entries = {}
        for s in series:
            for i in range(len(CATALOGUE[metric][s])):
                name, definition, _ = section_definition(s, i, unit_system)
                entries[name] = definition
        if sid in entries or not entries:
            entries.setdefault(sid, document["sections"][sid])
        order = sorted(entries, key=lambda n: (from_definition(dict(entries[n])).A, n))
        options[sid] = order
        sections.update(entries)
    return {"groups": groups, "options": options, "sections": sections}


def design_document(document: dict, catalogue: dict, assignment: dict) -> dict:
    """The model with the section of each member replaced by ``assignment[member]`` (members not
    named keep their section).  Raises ``ValueError`` for an unknown member or a section the member's
    group may not take."""
    doc = curveballs.record.plain_copy(document)
    group_of = {mid: sid for sid, members in catalogue["groups"].items() for mid in members}
    for mid, name in assignment.items():
        if mid not in group_of:
            raise ValueError(f"unknown member {mid!r}")
        if not isinstance(name, str) or name not in catalogue["options"][group_of[mid]]:
            raise ValueError(f"member {mid}: section {name!r} is not one of the sections allowed for it")
        doc["sections"].setdefault(name, curveballs.record.plain_copy(catalogue["sections"][name]))
        doc["members"][mid]["section"] = name
    used = {m["section"] for m in doc["members"].values()}
    doc["sections"] = {k: v for k, v in doc["sections"].items() if k in used}
    return doc


def model_weight(document: dict) -> dict:
    """Steel weight of a model: sum of unit weight x area x flexible length over the members, in the
    force unit of the file (the volume, in length^3, when no material has a unit weight)."""
    model = build_model(curveballs.record.plain_copy(document))
    w = sum(model.materials[m.material].unit_weight * model.sections[m.section].A * m.length
            for m in model.members.values())
    if w > 0.0:
        return {"value": float(model.units.from_internal(w, "force")), "unit": model.units.label("force"),
                "what": "weight"}
    vol = sum(model.sections[m.section].A * m.length for m in model.members.values())
    return {"value": float(vol), "unit": "mm^3", "what": "volume"}


def design_feasible(document: dict) -> tuple:
    """-> (feasible, reason): the model is analysed (strict) and every AISC 360-22 member check passes
    (design status ``pass`` or ``conditional``)."""
    from .. import design as design_module
    from ..solve import analyze
    try:
        model = build_model(curveballs.record.plain_copy(document))
        results = analyze(model, strict=True)
    except ModelInvalid as e:
        return False, "analysis refused: " + ", ".join(sorted({p.code for p in e.problems}))
    dres = design_module.check(model, results)
    if dres.status not in ("pass", "conditional"):
        failing = sorted(k for k, m in dres.members.items() if m.status not in ("pass", "conditional"))
        return False, f"design status {dres.status} (members {', '.join(failing[:6])})"
    return True, ""


def reference_design(document: dict, catalogue: dict) -> dict | None:
    """The REFERENCE design of a design task: best found by this deterministic search, NOT a proven
    optimum.

    1. every group takes the heaviest section it may take; if that is not feasible there is no
       reference (returns None);
    2. the groups are visited in the order of the file; for a group the sections lighter than its
       current one are tried from the lightest upwards with the other groups unchanged, and the first
       feasible one is taken;
    3. step 2 is repeated until a pass changes nothing (at most ``DESIGN_MAX_PASSES`` passes).

    Feasible: ``design_feasible``.  Returns ``{"sections": {member: name}, "groupSections", "weight",
    "evaluations", "search"}``."""
    opts = catalogue["options"]
    current = {sid: len(o) - 1 for sid, o in opts.items()}
    n_eval = 0

    def assignment(idx):
        return {mid: opts[sid][idx[sid]] for sid, members in catalogue["groups"].items() for mid in members}

    def feasible(idx):
        nonlocal n_eval
        n_eval += 1
        return design_feasible(design_document(document, catalogue, assignment(idx)))[0]

    if not feasible(current):
        return None
    for _ in range(DESIGN_MAX_PASSES):
        changed = False
        for sid in opts:
            for i in range(current[sid]):
                trial = dict(current, **{sid: i})
                if feasible(trial):
                    current, changed = trial, True
                    break
        if not changed:
            break
    final = assignment(current)
    return {"sections": final, "groupSections": {sid: opts[sid][i] for sid, i in current.items()},
            "weight": model_weight(design_document(document, catalogue, final)), "evaluations": n_eval,
            "search": "heaviest sections, then group by group the lightest feasible section with the other groups "
                      f"fixed, repeated until nothing changes (at most {DESIGN_MAX_PASSES} passes); best found by this "
                      "search, not a proven optimum"}


# ============================================================================================ tasks
def task_types_for(case: dict) -> dict:
    """{task type: None if it can be built on the case, else the reason it cannot}."""
    cls = case["effectClass"]
    bv, bdv, dv = case["base"]["verification"], case["base"]["designVerification"], case.get("designVerification")
    out = {}
    out["predict"] = ("the curve-balled model is a mechanism: it has no results to predict"
                      if cls == "creates-mechanism" else None)
    out["diagnose"] = None
    if not _verified(bv):
        out["repair"] = "the reference of the intended model is not verified: " + _quarantine_text(bv)
    elif cls == "design-only" and not _verified(bdv):
        out["repair"] = ("the repair of a design-only change is judged on design results, and the design results "
                         "of the intended model are not verified: " + _quarantine_text(bdv))
    else:
        out["repair"] = None
    if cls == "creates-mechanism":
        out["design"] = "the curve-balled model is a mechanism: it cannot be designed as given"
    elif not _verified(dv):
        out["design"] = ("the design results of the given model are not verified (with the design results, tier F): "
                         + _quarantine_text(dv))
    else:
        out["design"] = None
    return out


_TOOLS = {
    "tool-allowed": "Mode tool-allowed: you may run the analysis tool (run_analysis) and read its results.",
    "no-solver": "Mode no-solver: the analysis tool is disabled; inspect the file and reason from it.",
}
_CONVENTIONS = ("Sign conventions, axes and units are those of docs/CONVENTIONS.md; every number you give is in the "
                "units of the model file.")


def _statement(task_type: str, mode: str, family: str, body: list, answer_example: dict) -> str:
    head = {"predict": "TASK (predict): give the named results of the model in the file below, exactly as it is given.",
            "diagnose": "TASK (diagnose): the model file below was derived from the intended design described here, "
                        "but something in it was changed in a way that matters. Say what is wrong and where.",
            "repair": "TASK (repair): the model file below was derived from the intended design described here, but "
                      "something in it was changed in a way that matters. Repair the model.",
            "design": "TASK (design): choose the member sizes of the model below so that every design check passes "
                      "with the least steel weight."}[task_type]
    return "\n".join([head, f"Structure: {family.replace('_', ' ')}. {_TOOLS[mode]}", _CONVENTIONS, ""] + body
                     + ["", "Submit your answer as JSON of this form:", json.dumps(answer_example, indent=1)])


def make_task(case: dict, task_type: str, mode: str, *, foreign_cache=None) -> tuple:
    """(task, answers) of a verified case; raises :class:`TaskUnavailable` when the task type cannot
    be served on it with a verified reference (the reason says why)."""
    if task_type not in TASK_TYPES:
        raise ValueError(f"task type must be one of {TASK_TYPES}, not {task_type!r}")
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, not {mode!r}")
    if case.get("status") != "verified":
        raise TaskUnavailable("no verified case", case.get("skips"))
    why = task_types_for(case)[task_type]
    if why:
        raise TaskUnavailable(why)
    key, cls, record = case["key"], case["effectClass"], case["record"]
    family = key["family"]
    v, bv = case["verification"], case["base"]["verification"]
    dv, bdv = case.get("designVerification"), case["base"]["designVerification"]
    pdoc, bdoc = case["perturbed"]["document"], case["base"]["document"]
    draw = Draw(stable_seed("task", TASK_VERSION, json.dumps(key, sort_keys=True), task_type))
    tid = task_id(key, task_type, mode)
    task = {"format": "benchfea-task", "formatVersion": "1.0.0", "id": tid, "taskType": task_type, "mode": mode,
            "family": family, "difficulty": key["difficulty"],
            "split": split_of(family, key["curveBall"], key["seed"]), "model": pdoc}
    answers = {"format": "benchfea-task-answers", "formatVersion": "1.0.0", "id": tid, "taskType": task_type,
               "mode": mode, "key": key, "generatorSeed": case["generatorSeed"], "skips": case["skips"],
               "effectClass": cls, "changeRecord": record, "verification": verification_summary(v),
               "designVerification": None if dv is None else verification_summary(dv)}
    scope = "analysis results of the given model"

    if task_type == "predict":
        design_ok = _verified(dv)
        if cls == "design-only":
            qs = []
        else:
            other = bv.get("reference") if bv["kind"] == "results" else None
            qs = select_quantities(v["reference"], other, draw)
        if design_ok:
            dq = _design_quantity(dv["designReference"], (bdv or {}).get("designReference"), len(qs) + 1)
            if dq:
                qs.append(dq)
                scope = "analysis and design results of the given model"
        if not qs:
            raise TaskUnavailable("no significant quantity with a verified reference to ask for")
        task["quantities"] = [_public_quantity(q) for q in qs]
        task["answerFormat"] = {"values": {q["id"]: "number" for q in qs}}
        body = ["Quantities:"] + [f"  {q['id']}: {_quantity_text(q)}" for q in qs]
        task["statement"] = _statement(task_type, mode, family, body, {"values": {q["id"]: 0.0 for q in qs}})
        answers["reference"] = {"values": {q["id"]: {"value": q["value"], "floor": q["floor"], "unit": q["unit"],
                                                      "kind": q["kind"], "dimension": q["dimension"]} for q in qs},
                                "resultsHash": _doc_hash(v["reference"])}

    elif task_type == "diagnose":
        cats = "; ".join(f"{k} ({t})" for k, t in CATEGORY_TEXT.items())
        body = ["INTENDED DESIGN:", describe(bdoc), "",
                "Answer with: the category of the change (one of: " + cats + "); the ids of the entities that were "
                "changed (members, nodes, supports, load cases, constraints, sections, design members ... as they are "
                "named in the file); the direction of the change (one of: " + ", ".join(DIRECTIONS) + " -- 'added' / "
                "'removed': the file has something the intended design does not / lacks something it has; 'reversed': "
                "numbers with the opposite sign; 'increased' / 'decreased': numbers of larger / smaller magnitude; "
                "'changed': anything else). If the model as given is a mechanism (unstable), also describe the "
                "mechanism: its kind (" + ", ".join(MECHANISM_KINDS) + "), its motion (" + ", ".join(MECHANISM_MOTIONS)
                + ") and the nodes that move; otherwise give null for 'mechanism'."]
        example = {"category": "supports", "entities": {"nodes": ["N1"]}, "direction": "changed", "mechanism": None}
        task["answerFormat"] = {"category": list(CATEGORY_TEXT), "entities": "{kind: [ids]}",
                                "direction": list(DIRECTIONS),
                                "mechanism": {"kind": list(MECHANISM_KINDS), "motion": list(MECHANISM_MOTIONS),
                                              "nodes": "[ids]", "or": None}}
        task["intendedDesign"] = describe(bdoc)
        task["statement"] = _statement(task_type, mode, family, body, example)
        ref = {"category": record["category"], "entities": record["entities"],
               "direction": change_direction(record["changes"]), "target": record["target"]}
        ref["mechanism"] = v["reference"]["diagnosis"] if cls == "creates-mechanism" else None
        # one wording for every diagnose task: the scope must not tell a mechanism from a stable model
        scope = ("the reference of the given model: its analysis results (with the design results where the "
                 "change is one of design data only) or, if it is unstable, its mechanism diagnosis")
        answers["reference"] = ref

    elif task_type == "repair":
        other = v.get("reference") if v["kind"] == "results" else None
        # only quantities the curve ball moved beyond the zero-credit error are asked for: the model as given
        # earns nothing on them, so a repair can be told from no repair
        qs = [] if cls == "design-only" else select_quantities(
            bv["reference"], other, draw, moved_by=1.5 * REPAIR_TOLERANCE["zero"], allow_unmoved=other is None)
        if not qs:
            raise TaskUnavailable("the curve ball moves no significant result of the intended design beyond the "
                                  "zero-credit error of a repair: a repair could not be told from no repair")
        design_ok = _verified(bdv)
        if design_ok:
            dq = _design_quantity(bdv["designReference"], (dv or {}).get("designReference"), len(qs) + 1)
            if dq:
                qs.append(dq)
        tol = REPAIR_TOLERANCE
        criteria = {"stable": "the submitted model is accepted by the validator and analysed by the engine",
                    "verified": "the reference of the submitted model passes the verification pipeline "
                                "(benchfea verify: every applicable tier)",
                    "quantities": f"each named quantity of the submitted model is within {tol['full']:.0%} of its "
                                  f"value in the intended design (partial credit up to {tol['zero']:.0%})",
                    "design": ("every design check passes (design status pass or conditional)" if design_ok else
                               "not a criterion of this task: the design results of the intended model are not "
                               "verified (tier F), so nothing is scored on them")}
        task["quantities"] = [_public_quantity(q) for q in qs]
        task["criteria"] = criteria
        task["intendedDesign"] = describe(bdoc)
        task["answerFormat"] = {"model": "the complete repaired model document (in the environment: omit it to "
                                         "submit the working copy)"}
        body = ["INTENDED DESIGN:", describe(bdoc), "", "Criteria:"] + [f"  {k}: {t}" for k, t in criteria.items()] + [
            "Named quantities (their values in the intended design are not given; they are compared physically, "
            "whatever units your model declares):"] + [
            f"  {q['id']}: {_quantity_text(q)}" for q in qs]
        task["statement"] = _statement(task_type, mode, family, body, {"model": {"formatVersion": "...", "...": "..."}})
        answers["reference"] = {"values": {q["id"]: {"value": q["value"], "floor": q["floor"], "unit": q["unit"],
                                                      "kind": q["kind"], "dimension": q["dimension"]} for q in qs},
                                "quantities": [_public_quantity(q) for q in qs], "designCriterion": design_ok,
                                "units": bdoc["units"], "intendedModel": bdoc, "resultsHash": _doc_hash(bv["reference"]),
                                "tolerance": dict(tol)}
        answers["baseVerification"] = verification_summary(bv)
        answers["baseDesignVerification"] = None if bdv is None else verification_summary(bdv)
        scope = "analysis results of the intended design" + (" and its design results" if design_ok else "")

    else:                                                         # design
        from ..verify import verify_model
        catalogue = catalogue_for(pdoc, case["units"])
        ref = reference_design(pdoc, catalogue)
        if ref is None:
            raise TaskUnavailable("no feasible design with the heaviest catalogue sections: no reference design")
        ref_doc = design_document(pdoc, catalogue, ref["sections"])
        rv = verify_model(ref_doc, design=True, budget=VERIFY_BUDGET, cache=foreign_cache).to_dict(volatile=False)
        if rv["status"] != "verified":
            raise TaskUnavailable("the reference design is not verified: " + _quarantine_text(rv))
        given = model_weight(pdoc)
        given["feasible"] = design_feasible(pdoc)[0]
        task["catalogue"] = catalogue
        task["criteria"] = {"feasible": "the model with your sections is analysed by the engine, every AISC 360-22 "
                                        "member check passes (design status pass or conditional) and its reference "
                                        "passes the verification pipeline with the design results (benchfea verify "
                                        "--design)",
                            "weight": "score = (reference weight) / (your weight), at most 1; the reference weight is "
                                      "the best found by the documented deterministic search, not a proven optimum"}
        task["answerFormat"] = {"sections": "{member id: section name from the catalogue of its group}"}
        body = ["Only the sections of the members may change. Members that share a section in the file form a group; "
                "each member may take any section listed for its group (they need not stay equal):"]
        body += [f"  group {sid} (members {', '.join(members)}): {', '.join(catalogue['options'][sid])}"
                 for sid, members in catalogue["groups"].items()]
        body += [f"The definitions of these sections are in the task under catalogue/sections. Weight of the model as "
                 f"given: {given['value']:.6g} {given['unit']}; with the sizes as given "
                 f"{'every design check passes' if given['feasible'] else 'at least one design check fails'}.",
                 "Criteria:"] + [
            f"  {k}: {t}" for k, t in task["criteria"].items()]
        first = next(iter(pdoc["members"]))
        task["statement"] = _statement(task_type, mode, family, body, {"sections": {first: pdoc["members"][first]["section"]}})
        answers["reference"] = {"sections": ref["sections"], "weight": ref["weight"], "search": ref["search"],
                                "evaluations": ref["evaluations"], "givenWeight": given}
        answers["referenceDesignVerification"] = verification_summary(rv)
        scope = "analysis and design results of the reference design"
    task["verification"] = public_verification(v, scope)
    return task, answers


def _doc_hash(doc: dict) -> str:
    return hashlib.sha256(json.dumps(doc, separators=(",", ":"), sort_keys=True).encode()).hexdigest()
