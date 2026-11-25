"""Scoring of an answer against the verified reference of a task (docs/BENCHMARK.md section 4).

    from benchfea.bench.scoring import score_answer
    out = score_answer(task, answers, answer)
    out["score"]        # float in [0, 1]
    out["valid"]        # False: the answer could not be read as an answer of this task (score 0)
    out["reasons"]      # machine-readable: per-quantity errors and credits, which criteria failed, ...
    out["summary"]      # one line

``score_answer`` never raises for a bad answer: a malformed answer scores 0 with the reason.

Tolerances (``TOLERANCES``): a quantity is judged on its relative error

    e = |submitted - reference| / max(|reference|, floor)

where ``floor`` is 0.001 of the largest value of the same dimension in the same combination (stored
with the reference).  Credit: 1 for ``e <= full``; 0 for ``e >= zero``; in between it falls linearly
in log(e): ``1 - log(e / full) / log(zero / full)``.

    predict, tool-allowed    full 1 %    zero  5 %     (every kind of quantity)
    predict, no-solver       full 5 %    zero 25 %
    repair (either mode)     full 2 %    zero 10 %     (against the intended design)

The full-credit tolerance is never tighter than twice the largest tolerance at which a reference can
have been verified (tier C admits a foreign solver up to an estimated error of 0.5 %).
"""
from __future__ import annotations

import json
import math

from .tasks import (convert_value, DIRECTIONS, MECHANISM_KINDS, MECHANISM_MOTIONS, REPAIR_TOLERANCE, VERIFY_BUDGET, design_document,
                    design_feasible, extract, model_weight)

__all__ = ["score_answer", "credit", "relative_error", "TOLERANCES", "DIAGNOSE_WEIGHTS", "MECHANISM_WEIGHTS",
           "FALSE_MECHANISM_FACTOR", "tolerance_for", "entity_match", "SCORING_VERSION", "reference_answer"]

SCORING_VERSION = "1"
_ALL = ("displacement", "reaction", "memberEndForce", "memberMaxMoment", "bucklingFactor", "designRatio")
#: (task type, mode) -> quantity kind -> {"full", "zero"} relative errors
TOLERANCES = {
    ("predict", "tool-allowed"): {k: {"full": 0.01, "zero": 0.05} for k in _ALL},
    ("predict", "no-solver"): {k: {"full": 0.05, "zero": 0.25} for k in _ALL},
    ("repair", "tool-allowed"): {k: dict(REPAIR_TOLERANCE) for k in _ALL},
    ("repair", "no-solver"): {k: dict(REPAIR_TOLERANCE) for k in _ALL},
}
#: diagnosis of a change: entities named, category, direction (direction counts only with an entity found)
DIAGNOSE_WEIGHTS = {"entities": 0.5, "category": 0.3, "direction": 0.2}
#: diagnosis of a mechanism: nodes that move, kind and motion of the best-matching mode of the verified diagnosis
MECHANISM_WEIGHTS = {"nodes": 0.6, "kind": 0.25, "motion": 0.15}
FALSE_MECHANISM_FACTOR = 0.5   # a mechanism claimed for a model that is stable


def tolerance_for(task_type: str, mode: str, kind: str) -> dict:
    return TOLERANCES[(task_type, mode)][kind]


def credit(error: float, full: float, zero: float) -> float:
    """Partial-credit curve (module docstring)."""
    if not math.isfinite(error):
        return 0.0
    if error <= full:
        return 1.0
    if error >= zero:
        return 0.0
    return 1.0 - math.log(error / full) / math.log(zero / full)


def relative_error(submitted: float, reference: float, floor: float) -> float:
    return abs(submitted - reference) / max(abs(reference), floor, 1e-300)


def _number(x):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        return None
    return float(x)


def _result(task_type, score, valid, reasons, summary) -> dict:
    score = float(min(1.0, max(0.0, score))) if math.isfinite(score) else 0.0
    return {"format": "benchfea-score", "formatVersion": "1.0.0", "scoringVersion": SCORING_VERSION,
            "taskType": task_type, "score": score, "valid": bool(valid), "reasons": reasons, "summary": summary}


def _invalid(task_type, why: str) -> dict:
    return _result(task_type, 0.0, False, {"error": why}, "score 0: " + why)


# ------------------------------------------------------------------------------------------ predict
def _score_quantities(values: dict, submitted_of, tol_of) -> tuple:
    rows, total = {}, 0.0
    for qid, ref in values.items():
        tol = tol_of(ref["kind"])
        row = {"reference": ref["value"], "unit": ref.get("unit"), "tolerance": dict(tol), "credit": 0.0}
        sub, problem = submitted_of(qid)
        if problem:
            row.update(submitted=None, error=problem)
        else:
            e = relative_error(sub, ref["value"], ref["floor"])
            row.update(submitted=sub, relativeError=e, credit=credit(e, tol["full"], tol["zero"]))
        total += row["credit"]
        rows[qid] = row
    return rows, (total / len(values) if values else 0.0)


def _predict(task, answers, answer) -> dict:
    values = answer.get("values") if isinstance(answer.get("values"), dict) else None
    if values is None:
        return _invalid("predict", "the answer has no object 'values' {quantity id: number}")
    mode = task["mode"]

    def submitted(qid):
        if qid not in values:
            return None, "missing"
        x = _number(values[qid])
        return (x, None) if x is not None else (None, "not a finite number")

    rows, score = _score_quantities(answers["reference"]["values"], submitted,
                                    lambda kind: tolerance_for("predict", mode, kind))
    unknown = sorted(str(k) for k in values if k not in rows)
    n_full = sum(1 for r in rows.values() if r["credit"] == 1.0)
    return _result("predict", score, True, {"quantities": rows, "unknownIds": unknown},
                   f"{n_full} of {len(rows)} quantities within tolerance; score {score:.3f}")


# ----------------------------------------------------------------------------------------- diagnose
def _ids(entities) -> set | None:
    """Ids named in an answer: {kind: [ids]} (kinds are not compared), a list of ids, or one id."""
    if isinstance(entities, str):
        return {entities}
    if isinstance(entities, (list, tuple)):
        items = entities
    elif isinstance(entities, dict):
        items = []
        for v in entities.values():
            items += list(v) if isinstance(v, (list, tuple)) else [v]
    else:
        return None
    out = set()
    for x in items:
        if isinstance(x, bool) or not isinstance(x, (str, int)):
            return None
        out.add(str(x))
    return out


def entity_match(named: set, record_entities: dict) -> dict:
    """Entities of an answer against the entities of a change record ``{kind: [ids]}``.

    precision = named ids that are entities of the record (of any kind) / named ids;
    recall    = the best, over the kinds of the record, of (ids of that kind named) / (ids of that kind)
                -- naming every member touched, or every node touched, is a complete answer;
    score     = F1 of the two."""
    every = {str(i) for ids in record_entities.values() for i in ids}
    if not named or not every:
        return {"precision": 0.0, "recall": 0.0, "score": 0.0}
    precision = len(named & every) / len(named)
    recall = max((len(named & {str(i) for i in ids}) / len(ids)) for ids in record_entities.values() if ids)
    f1 = 0.0 if precision + recall == 0.0 else 2.0 * precision * recall / (precision + recall)
    return {"precision": precision, "recall": recall, "score": f1}


def _f1(named: set, truth: set) -> float:
    if not named or not truth:
        return 0.0
    hit = len(named & truth)
    if not hit:
        return 0.0
    p, r = hit / len(named), hit / len(truth)
    return 2.0 * p * r / (p + r)


def _diagnose(task, answers, answer) -> dict:
    ref = answers["reference"]
    mech = answer.get("mechanism")
    if mech is not None and not isinstance(mech, dict):
        return _invalid("diagnose", "'mechanism' must be null or an object {kind, motion, nodes}")
    if ref.get("mechanism") is not None:
        modes = ref["mechanism"]["modes"]
        reasons = {"reference": "mechanism", "nMechanisms": ref["mechanism"]["nMechanisms"],
                   "weights": dict(MECHANISM_WEIGHTS)}
        if mech is None:
            reasons["error"] = "the model as given is a mechanism; the answer names none"
            return _result("diagnose", 0.0, True, reasons, "score 0: the mechanism was not named")
        named = _ids(mech.get("nodes", []))
        if named is None:
            return _invalid("diagnose", "'mechanism.nodes' must be a list of node ids")
        best, best_row = 0.0, None
        for n, m in enumerate(modes):
            loc = _f1(named, {str(x["node"]) for x in m.get("nodes", [])})
            kind = float(loc > 0.0 and mech.get("kind") == m.get("kind"))
            motion = float(loc > 0.0 and mech.get("motion") == m.get("motion"))
            s = MECHANISM_WEIGHTS["nodes"] * loc + MECHANISM_WEIGHTS["kind"] * kind + MECHANISM_WEIGHTS["motion"] * motion
            if best_row is None or s > best:
                best, best_row = s, {"mode": n, "nodes": loc, "kind": kind, "motion": motion,
                                     "description": m.get("description")}
        reasons["bestMode"] = best_row
        reasons["failed"] = [k for k in ("nodes", "kind", "motion") if best_row[k] < 1.0]
        return _result("diagnose", best, True, reasons, f"mechanism: nodes {best_row['nodes']:.2f}, kind "
                       f"{best_row['kind']:.0f}, motion {best_row['motion']:.0f}; score {best:.3f}")
    named = _ids(answer.get("entities", {}))
    if named is None:
        return _invalid("diagnose", "'entities' must be {kind: [ids]} or a list of ids")
    ent = entity_match(named, ref["entities"])
    category = float(answer.get("category") == ref["category"])
    direction = float(ent["score"] > 0.0 and answer.get("direction") == ref["direction"])
    w = DIAGNOSE_WEIGHTS
    score = w["entities"] * ent["score"] + w["category"] * category + w["direction"] * direction
    reasons = {"reference": "change", "weights": dict(w), "entities": ent, "category": category,
               "direction": direction, "expected": {"category": ref["category"], "entities": ref["entities"],
                                                     "direction": ref["direction"]},
               "failed": [k for k, x in (("entities", ent["score"]), ("category", category), ("direction", direction))
                          if x < 1.0]}
    if answer.get("direction") not in DIRECTIONS + (None,):
        reasons["note"] = f"direction is not one of {list(DIRECTIONS)}"
    if mech is not None:
        score *= FALSE_MECHANISM_FACTOR
        reasons["falseMechanism"] = (f"a mechanism was claimed but the model as given is stable: the score is "
                                     f"multiplied by {FALSE_MECHANISM_FACTOR}")
        reasons["failed"].append("mechanism")
    return _result("diagnose", score, True, reasons, f"entities {ent['score']:.2f}, category {category:.0f}, "
                   f"direction {direction:.0f}; score {score:.3f}")


# ------------------------------------------------------------------------------------------- repair
def _pipeline(document, design: bool, foreign_cache, budget) -> dict:
    from ..verify import verify_model
    return verify_model(document, design=design, budget=budget, cache=foreign_cache).to_dict(volatile=False)


def _gate_reason(v: dict) -> dict:
    out = {"status": v["status"], "tiersPassed": v["tiersPassed"], "tiersNotPassed": v["tiersNotPassed"]}
    if v["status"] == "refused":
        out["problems"] = [{"code": p.get("code"), "path": p.get("path"), "message": p.get("message")}
                           for p in (v.get("refusal") or {}).get("problems", [])[:8]]
    elif v.get("quarantine"):
        out["category"] = v["quarantine"].get("category")
        out["diagnosis"] = str(v["quarantine"].get("diagnosis"))[:600]
    if v["kind"] == "diagnosis":
        out["mechanism"] = [m.get("description") for m in v["reference"]["diagnosis"]["modes"][:4]]
    return out


def _repair(task, answers, answer, foreign_cache, budget) -> dict:
    model = answer.get("model")
    if not isinstance(model, dict):
        return _invalid("repair", "the answer has no object 'model' (the complete repaired model document)")
    ref = answers["reference"]
    want_design = bool(ref.get("designCriterion"))
    v = _pipeline(model, want_design, foreign_cache, budget)
    criteria = {"stable": v["kind"] == "results" and v["status"] != "refused",
                "verified": v["kind"] == "results" and v["status"] == "verified"}
    reasons = {"criteria": criteria, "pipeline": _gate_reason(v), "tolerance": ref["tolerance"]}
    if not criteria["verified"]:
        reasons["failed"] = [k for k, ok in criteria.items() if not ok]
        why = ("the submitted model is refused or unstable" if not criteria["stable"] else
               "the reference of the submitted model does not pass the verification pipeline")
        return _result("repair", 0.0, True, reasons, "score 0: " + why)
    results, design = v["reference"], v.get("designReference")
    by_id = {q["id"]: q for q in ref["quantities"]}

    def submitted(qid):
        try:
            x = float(extract(by_id[qid], results, design))
            return convert_value(x, ref["values"][qid]["dimension"], results["units"], ref["units"]), None
        except (KeyError, IndexError, TypeError, ValueError):
            return None, "the submitted model does not have this node / member / combination"

    tol = ref["tolerance"]
    rows, score = _score_quantities(ref["values"], submitted, lambda kind: tol)
    reasons["quantities"] = rows
    failed = [qid for qid, r in rows.items() if r["credit"] < 1.0]
    if want_design:
        criteria["design"] = bool(design) and design.get("status") in ("pass", "conditional")
        if not criteria["design"]:
            reasons["failed"] = ["design"] + failed
            return _result("repair", 0.0, True, reasons, "score 0: a design check of the submitted model fails")
    reasons["failed"] = failed
    return _result("repair", score, True, reasons,
                   f"{len(rows) - len(failed)} of {len(rows)} quantities within {tol['full']:.0%} of the intended "
                   f"design; score {score:.3f}")


# ------------------------------------------------------------------------------------------- design
def _design(task, answers, answer, foreign_cache, budget) -> dict:
    sections = answer.get("sections")
    if not isinstance(sections, dict):
        return _invalid("design", "the answer has no object 'sections' {member id: section name}")
    try:
        doc = design_document(task["model"], task["catalogue"], sections)
    except (ValueError, KeyError, TypeError) as e:
        return _invalid("design", str(e))
    ref = answers["reference"]
    ok, why = design_feasible(doc)
    criteria = {"feasible": ok, "verified": False}
    reasons = {"criteria": criteria, "referenceWeight": ref["weight"], "search": ref["search"]}
    if not ok:
        reasons.update(failed=["feasible"], infeasible=why)
        return _result("design", 0.0, True, reasons, "score 0: not feasible: " + why)
    v = _pipeline(doc, True, foreign_cache, budget)
    criteria["verified"] = v["kind"] == "results" and v["status"] == "verified"
    reasons["pipeline"] = _gate_reason(v)
    if not criteria["verified"]:
        reasons["failed"] = ["verified"]
        return _result("design", 0.0, True, reasons, "score 0: the reference of the submitted design does not pass "
                                                      "the verification pipeline with the design results")
    w = model_weight(doc)
    ratio = ref["weight"]["value"] / w["value"] if w["value"] > 0.0 else 0.0
    reasons.update(weight=w, weightRatio=ratio, failed=[] if ratio >= 1.0 else ["weight"])
    return _result("design", min(1.0, ratio), True, reasons,
                   f"feasible; weight {w['value']:.6g} {w['unit']} against the reference {ref['weight']['value']:.6g} "
                   f"(best found by the stated search); score {min(1.0, ratio):.3f}")


# --------------------------------------------------------------------------------------------- main
def score_answer(task: dict, answers: dict, answer, *, foreign_cache=None, budget: float = VERIFY_BUDGET) -> dict:
    """Score ``answer`` (a dict, or JSON text) for ``task`` against ``answers`` (the hidden part of the
    task).  Returns the dict described in the module docstring; never raises for a bad answer.

    Repairs and designs are re-solved through ``benchfea.verify.verify_model`` (``budget`` seconds;
    what the pipeline could not run in that time is not passed)."""
    task_type = task.get("taskType") if isinstance(task, dict) else None
    if isinstance(answer, (str, bytes)):
        try:
            answer = json.loads(answer)
        except (ValueError, RecursionError):
            return _invalid(task_type, "the answer is not valid JSON")
    if not isinstance(answer, dict):
        return _invalid(task_type, "the answer must be a JSON object")
    try:
        if task_type == "predict":
            return _predict(task, answers, answer)
        if task_type == "diagnose":
            return _diagnose(task, answers, answer)
        if task_type == "repair":
            return _repair(task, answers, answer, foreign_cache, budget)
        if task_type == "design":
            return _design(task, answers, answer, foreign_cache, budget)
    except Exception as e:                                   # a malformed answer must never crash the scorer
        return _invalid(task_type, f"the answer could not be scored ({type(e).__name__}: {str(e)[:300]})")
    return _invalid(task_type, f"unknown task type {task_type!r}")


def reference_answer(task: dict, answers: dict) -> dict:
    """An answer built from the hidden reference that scores 1.0 (used by the scripted "oracle" agent and
    the tests to show that every served task is solvable; it is not available to an agent)."""
    ref = answers["reference"]
    tt = task["taskType"]
    if tt == "predict":
        return {"values": {qid: r["value"] for qid, r in ref["values"].items()}}
    if tt == "diagnose":
        out = {"category": ref["category"], "entities": ref["entities"], "direction": ref["direction"],
               "mechanism": None}
        if ref.get("mechanism") is not None:
            m = ref["mechanism"]["modes"][0]
            out["mechanism"] = {"kind": m.get("kind"), "motion": m.get("motion"),
                                "nodes": [n["node"] for n in m.get("nodes", [])]}
        return out
    if tt == "repair":
        return {"model": ref["intendedModel"]}
    return {"sections": ref["sections"]}


_ = (MECHANISM_KINDS, MECHANISM_MOTIONS)
