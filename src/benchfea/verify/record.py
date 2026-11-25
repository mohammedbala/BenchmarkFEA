"""Records produced by the verification pipeline.  Plain data; ``to_dict()`` is JSON-serialisable.

Wording rule (docs/VERIFICATION.md): ``status == "verified"`` means only that the reference solution
passed every tier that applies, at the tolerances stated in the record.  It is a statement about
agreement between computations; nobody has reviewed, approved or certified anything.
"""
from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["TierResult", "VerificationRecord", "ProblemRecord", "TargetResult", "LibraryReport",
           "VERIFIED", "QUARANTINED", "REFUSED", "TIERS", "MEANING", "TIER_TITLES", "DIAGNOSIS_TITLES"]

VERIFIED, QUARANTINED, REFUSED = "verified", "quarantined", "refused"
TIERS = ("A", "B", "C", "D", "E", "F")
MEANING = ("'verified' means only: the reference solution passed every tier listed as applicable in this record, "
           "at the tolerances stated in it. It does not mean reviewed, approved or certified by anyone.")
TIER_TITLES = {"A": "closed form / ODE oracle", "B": "published problem", "C": "independent solvers",
               "D": "invariants", "E": "equivalence checks", "F": "design-check coverage"}
DIAGNOSIS_TITLES = {"N": "null space of an independently assembled stiffness", "K": "engine matrix: K phi = 0",
                    "C": "foreign-solver evidence", "E": "invariance of the diagnosis"}


def _clean(o):
    """numpy scalars/arrays -> plain JSON data; keys sorted by insertion (deterministic)."""
    import numpy as np
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, float) and o != o:
        return None
    return o


@dataclass
class TierResult:
    """Outcome of one tier.

    applies   does the tier apply to this reference solution?
    passed    True | False when it applies (None when it does not apply, or when it was not run because
              the caller excluded it -- then the record cannot be "verified")
    summary   one sentence
    largest   the largest differences found: [{"what", "value", "tolerance", ...}]
    tolerances the tolerances used
    details   everything else (per combination / per check)
    """
    tier: str
    applies: bool
    passed: bool | None = None
    summary: str = ""
    largest: list = field(default_factory=list)
    tolerances: dict = field(default_factory=dict)
    details: dict = field(default_factory=dict)
    seconds: float = 0.0
    skipped: bool = False             # excluded by the caller (tiers=...) or by the time budget
    timings: dict = field(default_factory=dict)       # seconds per part (volatile)
    title: str = ""
    # why an applicable tier is not passed: "disagreement" (two computations differ), "not-representable"
    # (too few independent solvers can take the model), "not-covered" (tier F), "not-run", "error"
    category: str = ""

    def to_dict(self, volatile: bool = True) -> dict:
        d = {"tier": self.tier, "title": self.title or TIER_TITLES.get(self.tier, self.tier), "applies": self.applies,
             "category": self.category,
             "passed": self.passed, "skipped": self.skipped, "summary": self.summary,
             "largestDifferences": _clean(self.largest), "tolerances": _clean(self.tolerances),
             "details": _clean(self.details)}
        if volatile:
            d["seconds"] = round(self.seconds, 4)
            if self.timings:
                d["timings"] = _clean(self.timings)
        return d


@dataclass
class VerificationRecord:
    """Result of ``verify_model``.

    kind        "results"   the reference solution is the results document of a stable model
                "diagnosis" the model is unstable; the reference solution is the mechanism diagnosis
    status      "verified"    every applicable tier passed
                "quarantined" an applicable tier failed, could not be completed, or was not run; the
                              reference solution must not be published; ``quarantine`` says why
                "refused"     the engine produced no reference solution (invalid or unsupported model)
    """
    kind: str
    status: str
    tiers: dict                       # "A".."F" (results) or the diagnosis tiers -> TierResult
    model: str = ""
    model_hash: str = ""
    quarantine: dict | None = None    # {"category", "diagnosis", "tiers": [...]}
    refusal: dict | None = None       # {"problems": [...]} when refused
    reference: dict | None = None     # results document ("benchfea-results") or {"diagnosis": ...}
    design_reference: dict | None = None
    timings: dict = field(default_factory=dict)
    versions: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    @property
    def publishable(self) -> bool:
        return self.status == VERIFIED

    @property
    def tiers_passed(self) -> list:
        return [k for k, t in self.tiers.items() if t.applies and t.passed is True]

    def to_dict(self, volatile: bool = True, reference: bool = True) -> dict:
        """``volatile=False`` leaves out timings (the rest is deterministic for a given model and
        software versions); ``reference=False`` leaves out the reference documents."""
        d = {"format": "benchfea-verification", "formatVersion": "1.0.0", "kind": self.kind,
             "status": self.status, "publishable": self.publishable, "meaning": MEANING,
             "model": self.model, "modelHash": self.model_hash,
             "tiersPassed": self.tiers_passed,
             "tiersNotPassed": [k for k, t in self.tiers.items() if t.applies and t.passed is not True],
             "tiersNotApplicable": [k for k, t in self.tiers.items() if not t.applies],
             "tiers": {k: t.to_dict(volatile) for k, t in self.tiers.items()},
             "quarantine": _clean(self.quarantine), "refusal": _clean(self.refusal),
             "notes": list(self.notes), "versions": _clean(self.versions)}
        if volatile:
            d["timings"] = {k: round(v, 4) for k, v in self.timings.items()}
        if reference:
            d["reference"] = _clean(self.reference)
            if self.design_reference is not None:
                d["designReference"] = _clean(self.design_reference)
        return d

    def summary(self) -> str:
        lines = [f"{self.model or 'model'}: {self.status.upper()} ({self.kind})"]
        for k, t in self.tiers.items():
            state = "n/a" if not t.applies else ("not run" if t.skipped else ("passed" if t.passed else "NOT PASSED"))
            lines.append(f"  tier {k} [{t.title or TIER_TITLES.get(k, k)}]: {state} -- {t.summary}")
        if self.quarantine:
            lines.append(f"  quarantine ({self.quarantine.get('category')}): {self.quarantine.get('diagnosis')}")
        if self.refusal:
            for p in self.refusal.get("problems", [])[:5]:
                lines.append(f"  refused: {p.get('message')}")
        return "\n".join(lines)


@dataclass
class TargetResult:
    index: int
    label: str                        # the path or extractor as text
    expected: float
    actual: float | None
    difference: float | None          # |actual - expected| (file units)
    relative: float | None            # difference / |expected| (None if expected == 0)
    tolerance: dict
    passed: bool
    unit: str = ""
    error: str = ""

    def to_dict(self) -> dict:
        return _clean({"index": self.index, "target": self.label, "expected": self.expected, "actual": self.actual,
                       "difference": self.difference, "relative": self.relative, "tolerance": self.tolerance,
                       "passed": self.passed, "unit": self.unit, "error": self.error})


@dataclass
class ProblemRecord:
    """Result of ``verify_problem``.

    outcome   "passed"       active problem, every target reproduced and no tier that ran on its model failed
              "failed"       active problem, at least one target missed or a tier failed -> must be fixed
                             or quarantined by its author (the CLI exits nonzero)
              "quarantined"  the problem file says quarantined (targets are still evaluated and reported)
              "invalid"      the problem directory could not be loaded
    """
    id: str
    directory: str
    outcome: str
    declared_status: str = ""
    tier: str = ""
    title: str = ""
    source: dict = field(default_factory=dict)
    basis: str = ""
    targets: list = field(default_factory=list)       # TargetResult
    model_record: VerificationRecord | None = None
    failures: list = field(default_factory=list)      # text
    quarantine: dict | None = None
    errors: list = field(default_factory=list)
    seconds: float = 0.0
    notes: list = field(default_factory=list)          # remarks on the tolerance policy found while loading
    limits_note: str = ""

    @property
    def cross_checked(self) -> bool | None:
        """Tier C of the model passed (None: not run / no results)."""
        t = None if self.model_record is None else self.model_record.tiers.get("C")
        return None if t is None or t.skipped else bool(t.passed)

    @property
    def largest_relative(self) -> float | None:
        v = [t.relative for t in self.targets if t.relative is not None]
        return max(v) if v else None

    def to_dict(self, volatile: bool = True, reference: bool = False) -> dict:
        d = {"id": self.id, "directory": self.directory, "outcome": self.outcome,
             "declaredStatus": self.declared_status, "tier": self.tier, "title": self.title,
             "source": _clean(self.source), "basis": self.basis,
             "targets": [t.to_dict() for t in self.targets], "largestRelativeDifference": self.largest_relative,
             "failures": list(self.failures), "quarantine": _clean(self.quarantine), "errors": list(self.errors),
             "notes": list(self.notes), "limitsNote": self.limits_note, "crossChecked": self.cross_checked,
             "modelRecord": None if self.model_record is None else self.model_record.to_dict(volatile, reference)}
        if volatile:
            d["seconds"] = round(self.seconds, 4)
        return d


@dataclass
class LibraryReport:
    root: str
    problems: list                    # ProblemRecord
    versions: dict = field(default_factory=dict)
    seconds: float = 0.0

    @property
    def ok(self) -> bool:
        """Nothing expected-active fails and nothing is unreadable."""
        return all(p.outcome in ("passed", "quarantined") for p in self.problems)

    def counts(self) -> dict:
        out = {"passed": 0, "failed": 0, "quarantined": 0, "invalid": 0}
        for p in self.problems:
            out[p.outcome] = out.get(p.outcome, 0) + 1
        return out

    def to_dict(self, volatile: bool = True) -> dict:
        d = {"format": "benchfea-verification-library", "formatVersion": "1.0.0", "root": self.root, "ok": self.ok,
             "meaning": MEANING, "counts": self.counts(), "versions": _clean(self.versions),
             "problems": [p.to_dict(volatile) for p in self.problems]}
        if volatile:
            d["seconds"] = round(self.seconds, 4)
        return d
