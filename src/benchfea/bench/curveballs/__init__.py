"""Curve balls: named, seeded, deterministic changes of a model that matter.

    from benchfea.bench import curveballs
    curveballs.available()                                  # names, grouped by category
    p = curveballs.apply(document, "release_added", seed=3) # Perturbed(document, record)
    p.record.changes                                        # JSON-pointer paths with before / after
    p.record.effect                                         # the MEASURED effect on analysis and design
    p.record.apply(document) == p.document                  # the record describes the difference exactly
    p.record.revert(p.document) == document

A curve ball refuses a model it does not fit (:class:`NotApplicable` with ``.reason``): the base model
lacks the feature, or no target has a material effect (thresholds and definitions: ``measure``).
Effect classes: ``changes-analysis`` | ``design-only`` | ``creates-mechanism`` | ``invalid-model``.
Each curve ball, with measured examples: docs/CURVEBALLS.md.
"""
from . import (constraints, offsets, releases, design_variables, supports, members, units, loads,  # noqa: F401
               second_order)
from .measure import CRITERIA, EFFECT_CLASSES, THRESHOLDS, Outcome, compare, evaluate
from .record import (RECORD_FORMAT, RECORD_VERSION, Change, ChangeRecord, RecordMismatch, apply_changes, diff,
                     invert_changes)
from .registry import (CURVEBALL_VERSION, MAX_CANDIDATES, Candidate, CurveBall, NotApplicable, Perturbed,
                       applicability, apply, available, get)

#: category -> one-line description, in the order of the task list
CATEGORIES = {
    "constraints": "rigid links and diaphragms",
    "offsets": "rigid end offsets",
    "releases": "member end releases",
    "design": "design variables (design-only)",
    "supports": "supports",
    "members": "member orientation and sections",
    "units": "units",
    "loads": "loads",
    "second_order": "second-order effects",
}

__all__ = ["apply", "available", "get", "applicability", "NotApplicable", "Perturbed", "CurveBall", "Candidate",
           "ChangeRecord", "Change", "RecordMismatch", "diff", "apply_changes", "invert_changes", "evaluate",
           "compare", "Outcome", "THRESHOLDS", "CRITERIA", "EFFECT_CLASSES", "CATEGORIES", "CURVEBALL_VERSION",
           "MAX_CANDIDATES", "RECORD_FORMAT", "RECORD_VERSION"]
