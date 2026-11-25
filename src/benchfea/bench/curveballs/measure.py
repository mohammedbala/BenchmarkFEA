"""Measuring what a curve ball does: analyse (and design) the base and the perturbed model and
compare them.  A curve ball is accepted only if the measured effect is MATERIAL.

Outcome of a document (:func:`evaluate`):

    "ok"        analysed, every invariant passed, every iteration converged (strict analysis)
    "unstable"  refused as a mechanism (``stability.mechanism``); carries the structured diagnosis
    "invalid"   refused by validation (``build_model``); carries the validator's problems
    "refused"   anything else the engine refuses (beyond the elastic critical load, not converged,
                invariant failed, slack-member mechanism, ill-conditioned, ...): there is no reference
                solution, so a curve ball never accepts such a target

Effect of base -> perturbed (:func:`compare`), all in physical (internal) units N, mm:

* ``analysisChange``: for every combination, the largest difference of
  - the nodal translations and the peak displacement magnitude of every member axis, relative to
    the largest base displacement of the combination,
  - the reaction forces and the member forces (end forces as global vectors, peak |N| and peak
    resultant shear of every member), relative to the largest base force of the combination,
  - the reaction moments and the member moments (end moments as global vectors, peak |T| and peak
    resultant bending moment of every member), relative to the largest base moment of the combination,
  over the nodes and members both models have; the largest of these over the combinations.  The
  difference is taken entity by entity, so a load path that moves is measured even where the peak
  values stay the same; member quantities are taken in global axes or as resultants, so turning a
  member's local axes alone is not a change.  (If a base scale is zero the perturbed one is used.)
* ``analysisIdentical``: every result array of every combination is bitwise identical.
* ``designRatioChange``: largest relative change of the governing design ratio of a design member.
* ``designPassFailChanged``: design members whose status moved between {pass, conditional} and fail.

Materiality (``THRESHOLDS``), per effect class:

    changes-analysis    the analysis is not identical AND (analysisChange >= 5 % OR a design member
                        changes pass / fail OR -- only for curve balls that change the RESISTANCE of
                        a member (its section or its orientation) -- designRatioChange >= 5 %).
                        For a curve ball that changes only the demands, a design ratio moves with
                        the member forces; a ratio that jumps while no governing force moves by 5 %
                        would be a discontinuity of the design checks, not an effect of the change.
    design-only         the analysis is bitwise identical AND (designRatioChange >= 5 % OR a design
                        member changes pass / fail)
    creates-mechanism   the base model is analysed and the perturbed one is refused as a mechanism
                        (the reference solution is the solver's diagnosis)
    invalid-model       the base model is valid and the perturbed one is refused by validation
                        (the reference answer is the list of validation problems)
"""
from __future__ import annotations

import copy
from collections import OrderedDict
from dataclasses import dataclass, field

import numpy as np

from ...model import build_model
from ...problems import ERROR, ModelInvalid, ModelUnstable
from ...verify.common import document_hash

__all__ = ["Outcome", "evaluate", "compare", "classify", "is_material", "THRESHOLDS", "CRITERIA", "EFFECT_CLASSES",
           "clear_cache"]

EFFECT_CLASSES = ("changes-analysis", "design-only", "creates-mechanism", "invalid-model")
#: materiality thresholds (relative changes); see the module docstring for the definitions
THRESHOLDS = {"analysisChange": 0.05, "designRatioChange": 0.05}
_PASSING = ("pass", "conditional")
_CACHE: "OrderedDict[tuple, Outcome]" = OrderedDict()
_CACHE_SIZE = 6


@dataclass
class Outcome:
    status: str                                   # "ok" | "unstable" | "invalid" | "refused"
    model: object = None
    results: object = None
    design: object = None
    problems: list = field(default_factory=list)  # [{"code", "path", "message"}] of a refusal
    diagnosis: dict | None = None                 # MechanismDiagnosis.to_dict() of an unstable model
    design_error: str = ""

    @property
    def ok(self) -> bool:
        return self.status == "ok"

    def summary(self) -> dict:
        """Plain numbers of an analysed model, in FILE units: per combination the largest nodal
        translation component, the largest displacement of a member axis (magnitude, stations
        included), the largest reaction force component and the largest member bending moment; the
        design status and largest ratio."""
        if not self.ok:
            return {"status": self.status, "problems": self.problems[:5]}
        u = self.model.units
        out = {"status": "ok", "units": {"force": u.force, "length": u.length}, "combinations": {}}
        for cid, c in self.results.combinations.items():
            mom = max((float(np.abs(m.forces[:, 4:]).max()) for m in c.members.values() if len(m.forces)),
                      default=0.0)
            defl = max((float(np.linalg.norm(m.displacements, axis=1).max()) for m in c.members.values()
                        if len(m.displacements)), default=0.0)
            out["combinations"][cid] = {
                "analysisType": c.analysis_type,
                "maxTranslation": float(u.from_internal(float(np.abs(c.displacements[:, :3]).max()), "length")),
                "maxMemberDisplacement": float(u.from_internal(defl, "length")),
                "maxReactionForce": float(u.from_internal(float(np.abs(c.reactions[:, :3]).max()), "force")),
                "maxMemberMoment": float(u.from_internal(mom, "moment"))}
        if self.design is not None:
            ratios = {k: float(m.governing.ratio) for k, m in self.design.members.items()
                      if m.governing is not None and m.governing.ratio is not None}
            worst = max(ratios, key=ratios.get) if ratios else None
            out["design"] = {"status": self.design.status, "maxRatio": ratios.get(worst), "maxRatioMember": worst,
                             "failing": sorted(k for k, m in self.design.members.items() if m.status == "fail")}
        return out


def clear_cache() -> None:
    _CACHE.clear()


def _problem_list(e: ModelInvalid) -> list:
    return [{"code": p.code, "path": p.path, "message": p.message} for p in e.problems if p.severity == ERROR]


def evaluate(document: dict, design: bool = True, cache: bool = True) -> Outcome:
    """Analyse (strict) and design one document.  Outcomes are cached on the document hash (a small
    least-recently-used cache: the base model of several curve balls is analysed once)."""
    from ... import design as design_module
    from ...solve import analyze

    key = (document_hash(document), bool(design))
    if cache and key in _CACHE:
        _CACHE.move_to_end(key)
        return _CACHE[key]
    try:
        model = build_model(copy.deepcopy(document))
    except ModelInvalid as e:
        out = Outcome("invalid", problems=_problem_list(e))
    else:
        try:
            results = analyze(model)
        except ModelUnstable as e:
            diagnosis = getattr(e, "diagnosis", None)
            codes = {p.code for p in e.problems}
            pure = (diagnosis is not None and codes == {"stability.mechanism"}
                    and diagnosis.classification == "mechanism")
            out = Outcome("unstable" if pure else "refused", model=model, problems=_problem_list(e),
                          diagnosis=diagnosis.to_dict() if diagnosis is not None else None)
        except ModelInvalid as e:
            out = Outcome("refused", model=model, problems=_problem_list(e))
        else:
            out = Outcome("ok", model=model, results=results)
            if design and document.get("design") is not None:
                try:
                    out.design = design_module.check(model, results)
                except Exception as exc:                  # a design refusal is reported, never hidden
                    out.design_error = f"{type(exc).__name__}: {exc}"
    if cache:
        _CACHE[key] = out
        while len(_CACHE) > _CACHE_SIZE:
            _CACHE.popitem(last=False)
    return out


# --------------------------------------------------------------------------------------- compare
def _member_arrays(model, c):
    """member id -> (forces (8,), moments (8,), displacement (1,)), independent of the member's roll
    angle: the end forces at i and j as GLOBAL vectors (signed), the peak |N| and peak resultant
    shear along the member; the end moments as global vectors, the peak |T| and peak resultant
    bending moment; the peak magnitude of the displacement of the member axis."""
    out = {}
    for mid, m in c.members.items():
        rt = np.asarray(model.members[mid].R, dtype=float).T
        e = np.asarray(m.end_forces, dtype=float)
        if len(m.forces):
            f = np.asarray(m.forces, dtype=float)
            peaks = (float(np.abs(f[:, 0]).max()), float(np.hypot(f[:, 1], f[:, 2]).max()),
                     float(np.abs(f[:, 3]).max()), float(np.hypot(f[:, 4], f[:, 5]).max()))
            u = float(np.linalg.norm(np.asarray(m.displacements, dtype=float), axis=1).max())
        else:
            peaks, u = (0.0, 0.0, 0.0, 0.0), 0.0
        out[mid] = (np.concatenate([rt @ e[0:3], rt @ e[6:9], peaks[:2]]),
                    np.concatenate([rt @ e[3:6], rt @ e[9:12], peaks[2:]]), np.array([u]))
    return out


def _identical(a, b) -> bool:
    ra, rb = a.results, b.results
    if list(ra.combinations) != list(rb.combinations) or list(ra.buckling) != list(rb.buckling):
        return False
    for cid, ca in ra.combinations.items():
        cb = rb.combinations[cid]
        if list(ca.members) != list(cb.members):
            return False
        for x, y in ((ca.displacements, cb.displacements), (ca.reactions, cb.reactions)):
            if not np.array_equal(x, y):
                return False
        for mid, ma in ca.members.items():
            mb = cb.members[mid]
            for x, y in ((ma.end_forces, mb.end_forces), (ma.x, mb.x), (ma.forces, mb.forces),
                         (ma.displacements, mb.displacements)):
                if not np.array_equal(x, y):
                    return False
    for cid, ba in ra.buckling.items():
        if not np.array_equal(ba.factors, rb.buckling[cid].factors):
            return False
    return True


def compare(base: Outcome, pert: Outcome) -> dict:
    """The measured effect base -> perturbed (both analysed).  Plain JSON data."""
    if not (base.ok and pert.ok):
        raise ValueError("compare needs two analysed models")
    nodes_b, nodes_p = list(base.model.nodes), list(pert.model.nodes)
    common_nodes = [n for n in nodes_b if n in pert.model.nodes]
    ib = [nodes_b.index(n) for n in common_nodes]
    ip = [nodes_p.index(n) for n in common_nodes]
    worst = {"value": 0.0}
    by_combination = {}
    for cid, cb in base.results.combinations.items():
        cp = pert.results.combinations.get(cid)
        if cp is None:
            continue
        mb, mp = _member_arrays(base.model, cb), _member_arrays(pert.model, cp)
        common = [m for m in mb if m in mp]
        fields = []          # (quantity, scale kind, ids, base array, perturbed array)
        fields.append(("nodal translation", "u", common_nodes, cb.displacements[ib, :3], cp.displacements[ip, :3]))
        fields.append(("reaction force", "f", common_nodes, cb.reactions[ib, :3], cp.reactions[ip, :3]))
        fields.append(("reaction moment", "m", common_nodes, cb.reactions[ib, 3:], cp.reactions[ip, 3:]))
        if common:
            fields.append(("member force", "f", common, np.array([mb[m][0] for m in common]),
                           np.array([mp[m][0] for m in common])))
            fields.append(("member moment", "m", common, np.array([mb[m][1] for m in common]),
                           np.array([mp[m][1] for m in common])))
            fields.append(("member displacement", "u", common, np.array([mb[m][2] for m in common]),
                           np.array([mp[m][2] for m in common])))
        scale_b, scale_p = {"u": 0.0, "f": 0.0, "m": 0.0}, {"u": 0.0, "f": 0.0, "m": 0.0}
        for _, kind, _, a, b in fields:
            if a.size:
                scale_b[kind] = max(scale_b[kind], float(np.abs(a).max()))
                scale_p[kind] = max(scale_p[kind], float(np.abs(b).max()))
        entry = {"value": 0.0}
        for quantity, kind, ids, a, b in fields:
            scale = scale_b[kind] if scale_b[kind] > 0.0 else scale_p[kind]
            if not a.size or scale <= 0.0:
                continue
            d = np.abs(b - a).max(axis=1)
            k = int(d.argmax())
            rel = float(d[k] / scale)
            if rel > entry["value"]:
                entry = {"value": rel, "quantity": quantity, "where": ids[k]}
        by_combination[cid] = round(entry["value"], 6)
        if entry["value"] > worst["value"]:
            worst = {**entry, "combination": cid}
    out = {"analysisIdentical": _identical(base, pert),
           "analysisChange": float(worst["value"]),
           "analysisChangeAt": {k: v for k, v in worst.items() if k != "value"},
           "analysisChangeByCombination": by_combination}
    # ---- design
    ratio_change, where = None, None
    flips, status_changes = [], []
    if base.design is not None and pert.design is not None:
        ratio_change = 0.0
        for mid, a in base.design.members.items():
            b = pert.design.members.get(mid)
            if b is None:
                continue
            if a.status != b.status:
                status_changes.append({"member": mid, "base": a.status, "perturbed": b.status})
                if (a.status in _PASSING and b.status == "fail") or (a.status == "fail" and b.status in _PASSING):
                    flips.append(mid)
            ra = None if a.governing is None or a.governing.ratio is None else float(a.governing.ratio)
            rb = None if b.governing is None or b.governing.ratio is None else float(b.governing.ratio)
            if ra is not None and rb is not None and ra > 0.0:
                rel = abs(rb - ra) / ra
                if rel > ratio_change:
                    ratio_change = rel
                    where = {"member": mid, "base": round(ra, 6), "perturbed": round(rb, 6),
                             "baseLimitState": a.governing.name, "perturbedLimitState": b.governing.name}
        out["designStatus"] = {"base": base.design.status, "perturbed": pert.design.status}
    out["designRatioChange"] = ratio_change
    out["designRatioChangeAt"] = where
    out["designPassFailChanged"] = flips
    out["designStatusChanges"] = status_changes
    return out


def classify(base: Outcome, pert: Outcome, effect: dict | None) -> str | None:
    """The effect class the pair belongs to (None: no class -- e.g. the perturbed model is refused)."""
    if not base.ok:
        return None
    if pert.status == "unstable":
        return "creates-mechanism"
    if pert.status == "invalid":
        return "invalid-model"
    if not pert.ok:
        return None
    return "design-only" if effect["analysisIdentical"] else "changes-analysis"


#: criteria that can make an effect material, per effect class (a curve ball may narrow them)
CRITERIA = {"changes-analysis": ("analysisChange", "designPassFailChanged", "designRatioChange"),
            "design-only": ("designPassFailChanged", "designRatioChange")}


def is_material(effect_class: str, effect: dict | None, criteria=None) -> tuple:
    """-> (material?, criterion met or the reason why not).  ``criteria``: the criteria that count for
    this curve ball (default: all of ``CRITERIA[effect_class]``)."""
    if effect_class in ("creates-mechanism", "invalid-model"):
        return True, effect_class
    allowed = [c for c in CRITERIA[effect_class] if criteria is None or c in criteria]
    t = THRESHOLDS
    rc = effect.get("designRatioChange")
    for c in allowed:
        if c == "analysisChange" and effect["analysisChange"] >= t["analysisChange"]:
            return True, c
        if c == "designPassFailChanged" and effect.get("designPassFailChanged"):
            return True, c
        if c == "designRatioChange" and rc is not None and rc >= t["designRatioChange"]:
            return True, c
    parts = [f"analysis change {effect['analysisChange']:.3%}"]
    if rc is not None:
        parts.append(f"design ratio change {rc:.3%}" + ("" if "designRatioChange" in allowed
                                                         else " (does not count for this curve ball)"))
    return False, "below the materiality threshold: " + ", ".join(parts)
