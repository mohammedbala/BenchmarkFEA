"""Registry of the curve balls and the engine that applies one.

A curve ball is a named, seeded, deterministic transformation of a base model document:

    perturbed = apply(document, "release_added", seed=3)
    perturbed.document          # the new input document
    perturbed.record            # ChangeRecord: exact changes, entities, why it matters, measured effect

How one is applied (the same for every curve ball):

1. the base model is analysed and designed; a base model that is not analysable is refused;
2. the curve ball lists its CANDIDATE targets in the base model (members, ends, supports, loads,
   design blocks, ...) -- none: :class:`NotApplicable` with the precondition that is not met;
3. the candidates are put in a seeded order (SHA-256 of name and seed -> PCG64, as the generators);
   candidates that keep the model inside what the cross-check can follow (``tier_c_expected``) come
   before those that take it outside;
4. candidate after candidate (at most ``max_candidates``) the perturbed document is built, analysed
   and designed, and the effect is MEASURED against the base model (``measure.compare``).  The first
   candidate whose effect class is one the curve ball declares and whose effect is material
   (``measure.THRESHOLDS``) is taken; the others are listed in the record under ``rejected``;
5. no candidate left: :class:`NotApplicable` saying what was tried.

The change list of the record is the structural difference of the two documents (``record.diff``).
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Callable, NamedTuple

import numpy as np

from ...model import build_model
from ...problems import ModelInvalid
from ...verify.common import document_hash, raw_document
from ..families._common import Draw, stable_seed
from . import measure
from .record import ChangeRecord, diff, plain_copy

__all__ = ["CurveBall", "Candidate", "Context", "NotApplicable", "Perturbed", "register", "apply", "available",
           "get", "applicability", "CURVEBALL_VERSION", "MAX_CANDIDATES"]

CURVEBALL_VERSION = "1"        # bump when a change alters the output of any (document, name, seed)
MAX_CANDIDATES = 12            # candidates measured before giving up
_REGISTRY: "dict[str, CurveBall]" = {}


class NotApplicable(Exception):
    """The curve ball cannot be applied to this model; ``reason`` says why (missing feature, or no
    target with a material effect); ``tried`` lists the candidates measured."""

    def __init__(self, name: str, reason: str, tried: list | None = None):
        super().__init__(f"curve ball {name!r} is not applicable: {reason}")
        self.name, self.reason, self.tried = name, reason, list(tried or [])


class Candidate(NamedTuple):
    target: str                 # one line naming the target, e.g. "member B1, end i"
    edit: Callable              # edit(document) -> None, mutates a deep copy of the base document
    entities: dict              # {"members": [...], "nodes": [...], ...}
    info: dict = {}             # extra facts for the record (e.g. the factor applied)


@dataclass(frozen=True)
class CurveBall:
    name: str
    category: str               # group, = module name
    title: str                  # short phrase
    changes: str                # what it changes in the file
    why: str                    # why it matters physically
    moves: str                  # which result quantities it moves
    effect_classes: tuple       # allowed effect classes, e.g. ("changes-analysis",)
    precondition: str           # what the base model must have
    candidates: Callable = field(repr=False, compare=False, default=None)   # candidates(ctx) -> list[Candidate]
    #: materiality criteria that count (measure.CRITERIA); None = the default of the effect class WITHOUT
    #: "designRatioChange" for changes-analysis: a design ratio alone makes material only a curve ball
    #: that changes the resistance of a member (it then lists the criterion explicitly)
    criteria: tuple | None = None

    def criteria_for(self, effect_class: str) -> tuple:
        if self.criteria is not None:
            return tuple(c for c in measure.CRITERIA.get(effect_class, ()) if c in self.criteria)
        if effect_class == "changes-analysis":
            return ("analysisChange", "designPassFailChanged")
        return tuple(measure.CRITERIA.get(effect_class, ()))


class Context:
    """What a curve ball sees of the base model."""

    def __init__(self, document: dict, base: measure.Outcome, draw: Draw):
        self.doc, self.base, self.draw = document, base, draw
        self.model, self.results, self.design = base.model, base.results, base.design
        self.vertical = document["conventions"]["verticalAxis"]

    def per_metre(self) -> float:
        """Length of one metre in file units."""
        return float(self.model.units.from_internal(1000.0, "length"))


class Perturbed(NamedTuple):
    document: dict
    record: ChangeRecord

    def json(self) -> str:
        from ...io import dumps
        return dumps(self.document, "json")


def register(**spec):
    """Decorator: ``@register(name=..., category=..., ...)`` on the candidates function."""
    def wrap(fn):
        cb = CurveBall(candidates=fn, **spec)
        if cb.name in _REGISTRY:
            raise KeyError(f"duplicate curve ball {cb.name}")
        unknown = set(cb.effect_classes) - set(measure.EFFECT_CLASSES)
        if unknown or not cb.effect_classes:
            raise ValueError(f"{cb.name}: effect classes {sorted(unknown)}")
        _REGISTRY[cb.name] = cb
        return fn
    return wrap


def available() -> list:
    """Names of the curve balls, in registration order (grouped by category)."""
    return list(_REGISTRY)


def get(name: str) -> CurveBall:
    if name not in _REGISTRY:
        raise KeyError(f"unknown curve ball {name!r}; available: {', '.join(_REGISTRY)}")
    return _REGISTRY[name]


def _shuffle(items: list, draw: Draw) -> list:
    items = list(items)
    for i in range(len(items) - 1, 0, -1):               # Fisher-Yates on the stable integer draw
        j = draw.integer(0, i)
        items[i], items[j] = items[j], items[i]
    return items


def _node_on_member(model) -> list:
    """Nodes that lie on the axis of a member (its flexible length or a rigid end offset) without being
    an end of that piece.  The PyNite translator of ``benchfea.crosscheck`` refuses such a model
    (PyNite would connect the node to the member), so only one foreign solver is left."""
    size = max(float(np.ptp(np.array([n.xyz for n in model.nodes.values()]), axis=0).max()), 1.0)
    tol = 1e-9 * size
    hits = []
    for mid, m in model.members.items():
        xi, xj = model.nodes[m.i].xyz, model.nodes[m.j].xyz
        pieces = [(m.start, m.end, {m.i} if not np.any(m.offset_i) else set(), {m.j} if not np.any(m.offset_j) else set())]
        if np.any(m.offset_i):
            pieces.append((xi, m.start, {m.i}, set()))
        if np.any(m.offset_j):
            pieces.append((m.end, xj, set(), {m.j}))
        for a, b, ends_a, ends_b in pieces:
            ab = b - a
            L2 = float(ab @ ab)
            if L2 <= 0.0:
                continue
            for nid, n in model.nodes.items():
                if nid in ends_a or nid in ends_b:
                    continue
                t = float((n.xyz - a) @ ab) / L2
                if 0.0 < t < 1.0 and float(np.linalg.norm(n.xyz - (a + t * ab))) <= tol:
                    hits.append((nid, mid))
    return hits


def _cross_checkable(document: dict):
    """``bench.tier_c_expected`` of a document, plus the topology the PyNite translator refuses (a node
    on a member's axis).  An expectation derived from the model's features; None if it cannot be made."""
    from ..generate import tier_c_expected
    try:
        out = dict(tier_c_expected(document))
        hits = _node_on_member(build_model(plain_copy(document)))
    except ModelInvalid:
        return None
    except Exception:                                    # an expectation only; never blocks a curve ball
        return None
    if hits:
        nid, mid = hits[0]
        out["crossCheckable"] = False
        out["pynite"] = "unsupported"
        out["notes"] = list(out.get("notes", [])) + [
            f"node {nid} lies on the axis of member {mid} without being one of its ends (and {len(hits) - 1} more "
            "such pair(s)): the PyNite translator refuses this topology, which leaves one foreign solver"]
    return out


def _verification(effect_class: str, base_c, pert_c, diagnosis: dict | None = None) -> dict:
    if effect_class == "creates-mechanism":
        out = {"reference": "diagnosis", "route": "diagnosis tiers N, K, C, E of benchfea.verify "
               "(docs/VERIFICATION.md section 4): null space of an independently assembled stiffness, engine "
               "matrix, two foreign solvers under a probing load, invariance of the diagnosis"}
        if diagnosis and any(m.get("cause") == "no-stiffness" for m in diagnosis.get("modes", [])):
            out["expectedVerifiable"] = False
            out["route"] += ("; NOT expected to pass: the diagnosis reports a loaded node without any member "
                             "(cause no-stiffness), and the diagnosis tiers have been measured not to confirm this "
                             "kind of mechanism (docs/CURVEBALLS.md section 4) -- the pipeline quarantines it")
        return out
    if effect_class == "invalid-model":
        return {"reference": "validation-problems", "route": "none: the model is refused by validation, there "
                "is no analysis to verify (benchfea verify reports 'refused'); the reference answer is the "
                "validator's list of problems, reproducible with 'benchfea validate'"}
    out = {"reference": "results" if effect_class == "changes-analysis" else "results + design"}
    if pert_c is None:
        out["route"] = "benchfea verify (tiers as applicable); cross-check expectation not available"
        return out
    out["crossCheckable"] = {"base": None if base_c is None else bool(base_c["crossCheckable"]),
                             "perturbed": bool(pert_c["crossCheckable"])}
    if pert_c["crossCheckable"]:
        out["route"] = ("benchfea verify: every applicable tier; tier C (two independent solvers) is expected to "
                        "be able to follow the perturbed model")
    else:
        out["route"] = ("benchfea verify: tier C is NOT expected to pass -- fewer than two independent solvers "
                        "are expected to represent the perturbed model; the pipeline then quarantines the reference "
                        "(category not-representable), and a quarantined reference is not publishable as verified")
        out["notes"] = list(pert_c.get("notes", []))
    if effect_class == "design-only":
        out["route"] += "; design results: tier F (--design)"
    return out


def _verify_run(doc: dict, mode: str, budget) -> dict:
    """Run ``benchfea.verify.verify_model`` on a perturbed document -> plain summary."""
    from ...verify import verify_model
    rec = verify_model(doc, design=(mode == "design"), budget=budget)
    q = rec.quarantine or {}
    return {"mode": mode, "status": rec.status, "kind": rec.kind, "tiersPassed": list(rec.tiers_passed),
            "tiersNotPassed": [k for k, t in rec.tiers.items() if t.applies and t.passed is not True],
            "category": q.get("category"), "diagnosis": q.get("diagnosis")}


def apply(document, name: str, seed: int = 0, *, effect_class: str | None = None,
          max_candidates: int = MAX_CANDIDATES, verify=False, verify_budget: float | None = None) -> Perturbed:
    """Apply curve ball ``name`` with ``seed`` to a model document (dict, path or ``io.Document``).

    Deterministic: the same ``(document, name, seed)`` gives the same perturbed document and record.
    ``effect_class`` restricts the accepted outcome to one of the classes the curve ball declares.
    ``verify``: False (default) -- the record only carries the EXPECTED verification route;
    ``"analysis"`` / ``"design"`` / True -- the reference solution of every candidate that would be
    accepted is put through ``benchfea.verify.verify_model`` (analysis results or mechanism diagnosis;
    ``"design"``: with the design results, tier F; True = ``"design"`` for a design-only effect, else
    ``"analysis"``) and a candidate whose reference is not ``verified`` is dropped like an immaterial
    one.  The outcome is stored in ``record.verification["run"]``.  An invalid model has nothing to
    verify.  Verification costs far more than the measurement (foreign solvers, equivalence series).
    Raises ``KeyError`` for an unknown name and :class:`NotApplicable` (with ``.reason``) when the
    base model lacks the feature or no target has a material effect (and a verified reference)."""
    if verify not in (False, True, "analysis", "design"):
        raise ValueError("verify must be False, True, 'analysis' or 'design'")
    cb = get(name)
    seed = int(seed)
    allowed = cb.effect_classes
    if effect_class is not None:
        if effect_class not in allowed:
            raise NotApplicable(name, f"it never gives the effect class {effect_class!r} (it gives: "
                                      f"{', '.join(allowed)})")
        allowed = (effect_class,)
    base_doc = plain_copy(raw_document(document))          # no shared sub-objects: an edit touches one place
    base = measure.evaluate(base_doc)
    if not base.ok:
        codes = ", ".join(sorted({p["code"] for p in base.problems})) or base.status
        raise NotApplicable(name, f"the base model is not analysable ({base.status}: {codes})")
    if base.design_error:
        raise NotApplicable(name, f"the design checks of the base model failed ({base.design_error})")
    ctx = Context(base_doc, base, Draw(stable_seed("curveball", CURVEBALL_VERSION, name, seed)))
    candidates = list(cb.candidates(ctx) or [])
    if not candidates:
        raise NotApplicable(name, "precondition not met: " + cb.precondition)
    candidates = _shuffle(candidates, ctx.draw)
    base_c = _cross_checkable(base_doc)

    # build the perturbed documents of the first candidates; those that stay cross-checkable first
    pool = []
    for cand in candidates[:3 * max_candidates]:
        doc = plain_copy(base_doc)
        cand.edit(doc)
        changes = diff(base_doc, doc)
        if not changes:
            continue
        pert_c = _cross_checkable(doc)
        degrades = bool(base_c and base_c["crossCheckable"] and pert_c is not None and not pert_c["crossCheckable"])
        pool.append((degrades, cand, doc, changes, pert_c))
    pool.sort(key=lambda t: t[0])                         # stable: the seeded order is kept within each part
    tried = []
    for degrades, cand, doc, changes, pert_c in pool[:max_candidates]:
        pert = measure.evaluate(doc, cache=False)
        if pert.status == "refused":
            codes = ", ".join(sorted({p["code"] for p in pert.problems}))
            tried.append({"target": cand.target, "reason": f"the perturbed model is refused by the engine ({codes}): "
                                                           "no reference solution"})
            continue
        if pert.ok and pert.design_error:
            tried.append({"target": cand.target, "reason": f"design checks failed ({pert.design_error})"})
            continue
        effect = measure.compare(base, pert) if pert.ok else None
        cls = measure.classify(base, pert, effect)
        if cls not in allowed:
            why = f"effect class {cls!r} is not {' / '.join(allowed)}"
            if pert.status in ("invalid", "unstable"):
                why += " (" + ", ".join(sorted({p["code"] for p in pert.problems})) + ")"
            tried.append({"target": cand.target, "reason": why})
            continue
        material, criterion = measure.is_material(cls, effect, cb.criteria_for(cls))
        if not material:
            tried.append({"target": cand.target, "reason": criterion,
                          "analysisChange": round(effect["analysisChange"], 6),
                          "designRatioChange": None if effect.get("designRatioChange") is None
                          else round(effect["designRatioChange"], 6)})
            continue
        measured = {"class": cls, "material": True, "criterion": criterion,
                    "criteria": list(cb.criteria_for(cls)) or [cls], "thresholds": dict(measure.THRESHOLDS),
                    "base": base.summary()}
        if pert.ok:
            measured.update(effect)
            measured["perturbed"] = pert.summary()
            reference = {"kind": "results", "note": "the reference solution is the analysis (and design) of the "
                                                    "perturbed model, once verified"}
        elif cls == "creates-mechanism":
            measured["perturbed"] = {"status": "unstable", "nMechanisms": pert.diagnosis["nMechanisms"]}
            reference = {"kind": "diagnosis", "note": "the perturbed model is a mechanism: the reference solution is "
                         "the solver's structured mechanism diagnosis, not a set of results",
                         "diagnosis": pert.diagnosis, "problems": pert.problems}
        else:
            measured["perturbed"] = {"status": "invalid"}
            reference = {"kind": "validation-problems", "note": "the perturbed model is refused by validation: the "
                         "reference answer is the list of problems", "problems": pert.problems}
        measured.update(cand.info)
        verification = _verification(cls, base_c, pert_c, pert.diagnosis)
        if verify and cls != "invalid-model":
            mode = verify if verify in ("analysis", "design") else ("design" if cls == "design-only" else "analysis")
            ran = _verify_run(doc, mode, verify_budget)
            if ran["status"] != "verified":
                tried.append({"target": cand.target, "reason": f"the reference of the perturbed model is not verified "
                              f"({ran['status']}, category {ran['category']}; tiers not passed: "
                              f"{', '.join(ran['tiersNotPassed']) or '-'})", "verification": ran})
                continue
            verification["run"] = ran
        record = ChangeRecord(
            name=cb.name, seed=seed, category=cb.category, changes=changes,
            entities={k: list(v) for k, v in cand.entities.items()}, target=cand.target, why=cb.why,
            effect_class=cls, precondition=cb.precondition, effect=measured, reference=reference,
            verification=verification, rejected=tried,
            base_hash=document_hash(base_doc), perturbed_hash=document_hash(doc), title=cb.title,
            version=CURVEBALL_VERSION)
        return Perturbed(doc, record)
    n = len(pool)
    more = "" if n <= max_candidates else f" (the first {max_candidates} of {n} candidates were measured)"
    if not tried:
        raise NotApplicable(name, "no candidate changes the document" + more)
    raise NotApplicable(name, f"none of the {len(tried)} target(s) tried has a material effect of class "
                              f"{' / '.join(allowed)}{' with a verified reference' if verify else ''}{more}: "
                        + "; ".join(f"{t['target']}: {t['reason']}" for t in tried[:4])
                        + (" ..." if len(tried) > 4 else ""), tried)


def applicability(document, names=None, seed: int = 0) -> dict:
    """{curve ball: {"applicable": bool, "effectClass" | "reason": ...}} for one model."""
    out = {}
    for name in (names or available()):
        try:
            p = apply(document, name, seed)
        except NotApplicable as e:
            out[name] = {"applicable": False, "reason": e.reason}
        else:
            out[name] = {"applicable": True, "effectClass": p.record.effect_class, "target": p.record.target}
    return out
