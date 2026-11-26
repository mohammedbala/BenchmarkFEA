"""``verify_model``: run every tier that applies to a model and decide whether its reference solution
may be published as verified.  Policy and wording: docs/VERIFICATION.md.

A reference solution is "verified" only if EVERY tier that applies passed.  There is no vote: one
applicable tier that is not passed -- because two computations differ, because too few independent
solvers can represent the model, because a check could not be run -- makes the record "quarantined",
with the reason.  The record always lists which tiers passed.
"""
from __future__ import annotations

import copy
import platform
import sys
import time
from pathlib import Path

from ..model import build_model
from ..problems import ModelInvalid, ModelUnstable
from ..solve import analyze
from .common import SolveFailure, document_hash, raw_document
from .record import (DIAGNOSIS_TITLES, QUARANTINED, REFUSED, TIERS, VERIFIED, TierResult, VerificationRecord)

__all__ = ["verify_model", "software_versions", "tableau_evidence"]

FLAGGED_CODES = ("solve.invariant-failed", "solve.not-converged")
_CATEGORY_TEXT = {
    "disagreement": "two computations of the same quantity differ",
    "not-representable": "fewer than two independent solvers can represent the model, so it is not "
                         "cross-checked; such a reference solution is not publishable as verified",
    "not-covered": "a design limit state used by the reference solution lacks tier-F coverage",
    "not-run": "a tier was not run (excluded by the caller or stopped by the time budget)",
    "error": "a check could not be carried out",
    "not-verifiable": "this kind of reference solution cannot be verified by the pipeline",
}


def software_versions() -> dict:
    import numpy
    import scipy

    from .. import __version__
    from .tier_c import foreign_versions
    out = {"benchfea": __version__, "python": platform.python_version(), "numpy": numpy.__version__,
           "scipy": scipy.__version__}
    out.update({f"foreign:{k}": v for k, v in foreign_versions().items()})
    return out


def _tier_set(tiers) -> tuple:
    if tiers in (None, "all"):
        return TIERS
    if isinstance(tiers, str):
        tiers = [t for t in tiers.replace(",", "").replace(" ", "").upper()]
    bad = [t for t in tiers if t not in TIERS]
    if bad:
        raise ValueError(f"unknown tier(s) {bad}; tiers are {', '.join(TIERS)}")
    return tuple(t for t in TIERS if t in tiers)


def _problems(e: ModelInvalid) -> list:
    return [p.to_dict() for p in e.problems]


def _skipped(tier: str, applies: bool = True) -> TierResult:
    tr = TierResult(tier, applies, None, "not run (excluded by the caller)", skipped=True)
    tr.category = "not-run"
    return tr


def tableau_evidence(raw: dict, model, results_document: dict) -> dict | None:
    """EXTRA EVIDENCE ONLY (never one of the two solvers tier C requires): the independent "tableau"
    reference solver of ``tests/independent`` (a different method written for the review of the engine:
    no stiffness matrix; see its docstring), when it can be found: next to the package in a repository
    checkout (editable install), or under the current directory when the command is run from the root of
    a checkout with the package installed normally.  Returns
    ``{combination: {"largestRelativeDifference", "quantity"}}`` or None when it is not available."""
    import importlib
    candidates = (Path(__file__).resolve().parents[3] / "tests" / "independent",
                  Path.cwd() / "tests" / "independent")
    root = next((c for c in candidates if (c / "harness.py").exists() and (c / "tableau.py").exists()), None)
    if root is None:
        return None
    added = str(root) not in sys.path
    if added:
        sys.path.insert(0, str(root))
    try:
        harness = importlib.import_module("harness")
        out = {}
        for cid, combo in model.combinations.items():
            if cid not in results_document["combinations"]:
                continue
            try:
                d = harness.compare(copy.deepcopy(raw), cid, second_order=combo.analysis_type == "secondOrder",
                                    doc=results_document)
                worst = max(((k, v) for k, v in d.items() if not k.startswith("_")), key=lambda kv: kv[1])
                out[cid] = {"largestRelativeDifference": float(worst[1]), "quantity": worst[0]}
            except Exception as e:                         # extra evidence must never break the pipeline
                out[cid] = {"error": f"{type(e).__name__}: {str(e)[:200]}"}
        return out
    except Exception as e:
        return {"error": f"{type(e).__name__}: {str(e)[:200]}"}
    finally:
        if added and str(root) in sys.path:
            sys.path.remove(str(root))


def _decide(record: VerificationRecord) -> None:
    """status from the tiers: verified only if every applicable tier passed."""
    bad = [(k, t) for k, t in record.tiers.items() if t.applies and t.passed is not True]
    if not bad:
        record.status = VERIFIED
        return
    record.status = QUARANTINED
    cats = []
    for k, t in bad:
        if not t.category:
            t.category = "not-run" if t.skipped else "disagreement"
        if t.category not in cats:
            cats.append(t.category)
    # the most serious category names the quarantine
    order = ["disagreement", "error", "not-verifiable", "not-covered", "not-representable", "not-run"]
    category = min(cats, key=order.index)
    record.quarantine = {
        "category": category, "meaning": _CATEGORY_TEXT[category], "tiers": [k for k, _ in bad],
        "diagnosis": "; ".join(f"tier {k}: {t.summary}" for k, t in bad)}


def verify_model(source, *, tiers="all", design: bool = False, budget: float | None = None,
                 foreign_timeout: float = 120.0, cache=None, targets=None, extra_evidence: bool = False,
                 name: str | None = None) -> VerificationRecord:
    """Verify the reference solution of one model.

    source           path, dict, ``io.Document`` or a ``Model`` built from a file
    tiers            "all" or a selection such as "ACD" / ["C", "D"]; a tier that is left out is recorded
                     as not run and the record can then not be "verified"
    design           the reference solution includes the design results (tier F applies)
    budget           seconds for the whole call; tiers C and E stop starting new work when it is used up
                     (what was not run is not passed)
    foreign_timeout  seconds per foreign solve
    cache            ``ForeignCache`` (default: a process-wide in-memory cache)
    targets          a ``problems.ProblemSpec`` whose targets are evaluated as tier A or B (used by
                     ``verify_problem``)
    extra_evidence   also run the independent tableau solver of tests/independent, when present, and
                     record its largest difference under tier C (never counted as one of the two solvers)

    Returns a :class:`VerificationRecord`; never raises for a model problem.
    """
    from .tier_a import tier_a
    from .tier_c import tier_c
    from .tier_d import tier_d
    from .tier_e import tier_e
    t_start = time.perf_counter()
    deadline = None if budget is None else t_start + float(budget)
    run = _tier_set(tiers)
    versions = software_versions()
    timings = {}
    try:
        raw = raw_document(source)
    except ModelInvalid as e:
        return VerificationRecord("results", REFUSED, {}, name or str(source), "", refusal={"problems": _problems(e)},
                                  versions=versions, timings={"total": time.perf_counter() - t_start})
    label = name or raw.get("name") or (str(source) if isinstance(source, (str, Path)) else "model")
    h = document_hash(raw)

    def finish(rec: VerificationRecord) -> VerificationRecord:
        timings["total"] = time.perf_counter() - t_start
        rec.timings, rec.versions, rec.model, rec.model_hash = timings, versions, label, h
        return rec

    # ---------------------------------------------------------------- the engine's reference solution
    t0 = time.perf_counter()
    try:
        model = build_model(copy.deepcopy(raw))
    except ModelInvalid as e:
        return finish(VerificationRecord("results", REFUSED, {}, refusal={"stage": "validation", "problems": _problems(e)}))
    flagged = None
    try:
        results = analyze(model, strict=True)
    except ModelInvalid as e:
        timings["engine"] = time.perf_counter() - t0
        diagnosis = getattr(e, "diagnosis", None)
        partial = getattr(e, "refused", None)
        codes = [p.code for p in e.problems]
        if isinstance(e, ModelUnstable) and diagnosis is not None:
            return finish(_verify_diagnosis(raw, model, e, run, cache, foreign_timeout, deadline, timings))
        if partial is None and getattr(e, "results", None) is not None and all(c in FLAGGED_CODES for c in codes):
            flagged, results = e, e.results          # results exist but the engine itself flags them
        else:
            return finish(VerificationRecord("results", REFUSED, {},
                                             refusal={"stage": "analysis", "problems": _problems(e)}))
    timings["engine"] = time.perf_counter() - t0
    reference = results.to_document(model)
    record = VerificationRecord("results", QUARANTINED, {}, reference=reference)
    tr = {}
    # ---------------------------------------------------------------- D first (cheap, and decisive)
    tr["D"] = tier_d(model, results) if "D" in run else _skipped("D")
    if flagged is not None:
        tr["D"].passed = False
        tr["D"].category = "disagreement"
        tr["D"].summary = ("the engine flags its own results: " + "; ".join(p.message for p in flagged.problems[:3]))
        for k in "ACE":
            t = TierResult(k, True, None, "not run: the engine's results are flagged (tier D)", skipped=True)
            t.category = "not-run"
            tr[k] = t
    else:
        tr["A"] = tier_a(model, results) if "A" in run else _skipped("A")
        if "C" in run:
            tr["C"] = tier_c(raw, model, results, cache=cache, timeout=foreign_timeout, deadline=deadline)
            if extra_evidence:
                t1 = time.perf_counter()
                tr["C"].details["extraEvidence"] = {
                    "tableau": tableau_evidence(raw, model, reference),
                    "note": "independent tableau solver of tests/independent: extra evidence only, not one of "
                            "the two open-source solvers the tier requires"}
                tr["C"].timings["extraEvidence"] = round(time.perf_counter() - t1, 4)
        else:
            tr["C"] = _skipped("C")
        tr["E"] = tier_e(raw, model, results, deadline=deadline) if "E" in run else _skipped("E")
        if tr["E"].applies and tr["E"].passed is False:
            reasons = [c.get("reason") or "" for c in tr["E"].details["checks"].values()
                       if c.get("applies") and c.get("passed") is not True]
            tr["E"].category = ("not-run" if reasons and all("time budget" in r for r in reasons) else
                                "error" if reasons and all(r for r in reasons) else "disagreement")
    # ---------------------------------------------------------------- targets of a registered problem
    if targets is not None:
        from .problems import evaluate_targets
        design_doc = None
        if targets.needs_design:
            from .. import design as design_module
            design_doc = design_module.check(model, results).to_document(model)
        rows = evaluate_targets(targets, reference, design_doc)
        _merge_targets(tr, targets, rows)
    else:
        tr["B"] = TierResult("B", False, None, "not a registered published problem (see verify_problem)")
    if "A" not in tr:
        tr["A"] = _skipped("A")
    # ---------------------------------------------------------------- F
    if design:
        if "F" in run:
            from .. import design as design_module
            from .tier_f import tier_f
            t1 = time.perf_counter()
            dres = design_module.check(model, results)
            record.design_reference = dres.to_document(model)
            tr["F"] = tier_f(dres)
            tr["F"].seconds = time.perf_counter() - t1
        else:
            tr["F"] = _skipped("F")
    else:
        tr["F"] = TierResult("F", False, None, "design results are not part of this reference solution")
    record.tiers = {k: tr[k] for k in TIERS}
    for k, t in record.tiers.items():
        timings[f"tier{k}"] = t.seconds
    _decide(record)
    if record.status == VERIFIED:
        record.notes.append("verified: passed tier(s) " + ", ".join(record.tiers_passed)
                            + " at the tolerances stated in the record")
    return finish(record)


def _merge_targets(tr: dict, spec, rows: list) -> None:
    """Targets of a registered problem are tier A (derived here) or tier B (published) evidence."""
    missed = [r for r in rows if not r.passed]
    largest = sorted(({"what": r.label, "value": r.relative if r.relative is not None else r.difference,
                       "tolerance": r.tolerance, "expected": r.expected, "actual": r.actual} for r in rows
                      if r.difference is not None),
                     key=lambda d: -(d["value"] or 0.0))[:6]
    summary = (f"{len(rows)} target(s) of problem {spec.id!r} [{spec.basis}]: "
               + ("all reproduced" if not missed else f"{len(missed)} MISSED: " + "; ".join(
                   f"{r.label}: expected {r.expected:.8g}, engine {r.actual if r.actual is None else format(r.actual, '.8g')}"
                   + (f" ({r.error})" if r.error else "") for r in missed[:3])))
    details = {"problem": spec.id, "basis": spec.basis, "source": spec.source,
               "targets": [r.to_dict() for r in rows]}
    tier = spec.tier
    t = TierResult(tier, True, not missed, summary, largest, {"rule": "tolerance declared per target (docs/VERIFICATION.md)"},
                   details)
    t.category = "" if not missed else "disagreement"
    if tier == "A" and tr.get("A") is not None and tr["A"].applies and not tr["A"].skipped:
        auto = tr["A"]                                   # the automatic oracle comparison AND the targets
        t.passed = bool(auto.passed and not missed)
        t.summary = auto.summary + "; " + summary
        t.largest = (auto.largest + largest)[:8]
        t.tolerances = {**auto.tolerances, **t.tolerances}
        t.details = {**auto.details, **details}
        t.seconds = auto.seconds
        t.category = "" if t.passed else "disagreement"
    tr[tier] = t
    other = "B" if tier == "A" else "A"
    if other == "B":
        tr["B"] = TierResult("B", False, None, "not a published problem (targets derived here: tier A)")
    elif other not in tr:
        tr[other] = _skipped(other)


# ====================================================================================== unstable models
def _verify_diagnosis(raw, model, e, run, cache, timeout, deadline, timings) -> VerificationRecord:
    from .diagnosis import verify_diagnosis_tiers
    diag = e.diagnosis.to_dict()
    partial = getattr(e, "refused", None)
    reference = {"diagnosis": diag, "problems": _problems(e)}
    slack = None
    if partial:
        infos = [r.info for r in partial.values() if getattr(r, "diagnosis", None) is not None]
        slack = sorted({m for i in infos for m in list(i.get("slack", [])) + list(i.get("idle", []))})
        reference["refusedCombinations"] = {what: {"problems": [p.to_dict() for p in r.problems], "info": r.info}
                                            for what, r in partial.items()}

    def analyse(doc):
        try:
            m2 = build_model(copy.deepcopy(doc))
            analyze(m2, strict=True)
        except ModelInvalid as e2:
            d2 = getattr(e2, "diagnosis", None)
            return None if d2 is None else d2.to_dict()
        return None

    keys = tuple(k for k in ("N", "K", "C", "E") if (k in run or k in ("N", "K")))
    t0 = time.perf_counter()
    try:
        tiers = verify_diagnosis_tiers(raw, model, diag, analyse, exclude=tuple(slack or ()), cache=cache,
                                       timeout=timeout, deadline=deadline, run=keys)
    except SolveFailure as err:                              # pragma: no cover  (defensive)
        tiers = {"N": TierResult("N", True, False, f"could not be carried out: {err.reason}")}
    timings["diagnosis"] = time.perf_counter() - t0
    for k, t in tiers.items():
        t.title = DIAGNOSIS_TITLES.get(k, k)
        if t.skipped and not t.category:
            t.category = "not-run"
    rec = VerificationRecord("diagnosis", QUARANTINED, tiers, reference=reference)
    _decide(rec)
    if partial:
        # a mechanism that exists only with slack one-sided members: never "verified"
        rec.status = QUARANTINED
        rec.quarantine = {
            "category": "not-verifiable", "meaning": _CATEGORY_TEXT["not-verifiable"],
            "tiers": [k for k, t in tiers.items() if t.passed is not True],
            "diagnosis": ("the mechanism appears only with the tension-/compression-only member(s) "
                          + ", ".join(slack or []) + " inactive.  Checked: the structure without them is a "
                          "mechanism (tiers N, K: " + ("confirmed" if all(tiers[k].passed for k in ("N", "K")) else
                                                        "NOT confirmed") + ").  Not checked by any independent "
                          "computation: that these members are slack under this loading.")}
    elif rec.status == VERIFIED:
        rec.notes.append("verified diagnosis: passed tier(s) " + ", ".join(rec.tiers_passed)
                         + " (docs/VERIFICATION.md, 'Unstable models')")
    return rec
