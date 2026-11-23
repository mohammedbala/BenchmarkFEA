"""AISC 360-22 Chapter L -- serviceability (pp. 16.1-174, 16.1-175): deflections.

Chapter L gives NO numeric limit, load combination or reference length: L1
asks for limits "chosen with due regard to the intended function of the
structure" and L2 that deflections do not impair serviceability.  Every limit
checked here is therefore the USER's (design variable ``deflectionLimits``);
without one, no deflection check exists (reported as a default).  A deflection
limit is an enforced check (a member that exceeds it is not "passed"); the
slenderness recommendations of the User Notes of D1 and E2 are advisory checks
of Chapters D and E (limits = design variables ``maxSlendernessTension`` /
``maxSlendernessCompression``; the input format has no switch to enforce them).

Definition of the checked quantity (docs/CONVENTIONS.md section 4; this is the
definition of the input format, not of AISC):

* the station displacements of the design member (member axis, local axes,
  rigid-body motion included) are taken along the chain of its analysis
  members, x from the design member start;
* the CHORD is the straight line through the displaced positions of the two
  design-member ends; the deflection at x is the displacement minus the chord:

      delta(x) = u(x) - [ u(0) + (u(L) - u(0)) x / L ]

  for the local y and z components (the axial component is not a deflection);
* ``direction`` "y" -> |delta_y|, "z" -> |delta_z|, "resultant" ->
  sqrt(delta_y^2 + delta_z^2); the largest value over the result stations of
  each listed combination is compared with the limit;
* limit = L / spanRatio and / or the absolute value; when both are given the
  smaller one governs.  L = design member length.

Limits of this check:

* the maximum is taken over the result stations (``analysis.resultStations``
  per analysis member), not between them; a design member with no station
  between its ends (``resultStations`` = 2) has no measurable deflection
  relative to its chord -> ``notChecked`` (never a vacuous pass);
* a design member with a FREE end (cantilever, overhang): its deflection
  relative to the chord through its two ends does not measure the tip
  deflection, and Chapter L defines no reference for it -> ``notChecked``;
* a listed combination that is not in the results -> ``notChecked``;
* the combination results are used as analysed.  The User Note of L1 says the
  reduced stiffness of the direct analysis method is not intended for
  serviceability: a combination analysed with the stiffness reduction is
  checked but carries a note.
"""
from __future__ import annotations

import numpy as np

from ...results import DesignCheck
from ..base import LENGTH, NONE, OK, TEXT, DesignContext, not_checked

CLAUSE = "AISC 360-22 L1, L2 (the Specification sets no numeric limit: the limit is the user's)"
DIRECTIONS = ("y", "z", "resultant")


def name_of(index: int, direction: str) -> str:
    return f"deflection limit {index + 1} ({direction if direction == 'resultant' else 'local ' + direction})"


def member_displacements(model, combination_result, design_member_id: str):
    """Station displacements of a design member -> (x (n,), u (n, 3) local ux, uy, uz), x from its start."""
    p = model.design_members[design_member_id]
    first = model.members[p.members[0]]
    xs, us = [], []
    for mid in p.members:
        m, r = model.members[mid], combination_result.members[mid]
        xs.append(float(np.linalg.norm(m.start - first.start)) + np.asarray(r.x, float))
        us.append(np.asarray(r.displacements, float))
    return np.concatenate(xs), np.vstack(us)


def chord_deflection(x, u):
    """Displacement relative to the chord through the first and the last station -> (n, 3)."""
    x, u = np.asarray(x, float), np.asarray(u, float)
    span = x[-1] - x[0]
    s = ((x - x[0]) / span)[:, None]
    return u - (u[0] + (u[-1] - u[0]) * s)


def deflection_value(x, u, direction: str):
    """Largest deflection relative to the chord -> (value, x at which it occurs)."""
    d = chord_deflection(x, u)
    if direction == "y":
        v = np.abs(d[:, 1])
    elif direction == "z":
        v = np.abs(d[:, 2])
    elif direction == "resultant":
        v = np.hypot(d[:, 1], d[:, 2])
    else:
        raise ValueError(f"direction must be one of {DIRECTIONS}")
    k = int(np.argmax(v))
    return float(v[k]), float(x[k] - x[0])


def deflection_checks(model, results, design_member_id: str, ctx: DesignContext | None = None) -> list:
    """One check per entry of the design member's ``deflectionLimits`` (governing combination)."""
    p = model.design_members[design_member_id]
    limits = p.deflectionLimits or ()
    if not limits:
        return []
    ctx = ctx or DesignContext.from_model(model, design_member_id)
    source = p.source("deflectionLimits")
    L = p.length
    out = []
    for k, lim in enumerate(limits):
        name = name_of(k, lim.direction)
        candidates = [v for v in (L / lim.span_ratio if lim.span_ratio else None, lim.absolute) if v is not None]
        inputs = {"L": (L, LENGTH, "derived: design member length"),
                  "direction": (lim.direction, TEXT, source),
                  "combinations": (", ".join(lim.combinations), TEXT, source)}
        if lim.span_ratio:
            inputs["spanRatio"] = (lim.span_ratio, NONE, source)
        if lim.absolute is not None:
            inputs["absolute"] = (lim.absolute, LENGTH, source)
        kw = dict(dimension=LENGTH, inputs=inputs)
        if not candidates:
            out.append(not_checked(name, CLAUSE, "the limit gives neither spanRatio nor absolute", ctx, **kw))
            continue
        limit = min(candidates)
        if any(ctx.end_free):
            out.append(not_checked(
                name, CLAUSE, "the design member has a free end (cantilever / overhang): the deflection relative "
                "to the chord through its two ends does not measure the tip deflection, and neither the input "
                "format nor AISC 360-22 Chapter L defines another reference", ctx, **kw))
            continue
        if not lim.combinations:
            out.append(not_checked(name, CLAUSE, "no combination to check: the limit lists none and the model has "
                                   "no combination with purpose = serviceability", ctx, **kw))
            continue
        absent = [c for c in lim.combinations if c not in results.combinations]
        if absent:
            out.append(not_checked(name, CLAUSE, f"combination(s) {', '.join(absent)} are not in the analysis "
                                   "results", ctx, **kw))
            continue
        rows, interior = [], True
        for cid in lim.combinations:
            xs, us = member_displacements(model, results.combinations[cid], design_member_id)
            interior = interior and bool(np.any((xs > xs[0]) & (xs < xs[-1])))
            value, x = deflection_value(xs, us, lim.direction)
            rows.append((value, x, cid))
        if not interior:
            # only the two end stations: the deflection relative to the chord through them is zero by
            # definition, which would be a vacuous "pass"
            out.append(not_checked(name, CLAUSE, "the design member has no result station between its ends "
                                   "(analysis.resultStations = 2): the deflection relative to the chord through "
                                   "the ends is measured at the stations only and cannot be evaluated", ctx, **kw))
            continue
        value, x, cid = max(rows, key=lambda r: r[0])
        notes = ["deflection relative to the chord through the design member ends; maximum over the result "
                 "stations", "the limit is the user's: AISC 360-22 Chapter L sets no numeric limit"]
        reduced = [c for c in lim.combinations if c in model.combinations and model.combinations[c].stiffness_reduction]
        if reduced:
            notes.append(f"combination(s) {', '.join(reduced)} were analysed with the stiffness reduction of C2.3, "
                         "which is not intended for serviceability checks (L1 User Note)")
        unconverged = [c for c in lim.combinations if not results.combinations[c].converged]
        if unconverged:
            out.append(not_checked(name, CLAUSE, f"combination(s) {', '.join(unconverged)} did not converge", ctx,
                                   **kw))
            continue
        details = {"limit": (limit, LENGTH), "L / deflection": (L / value if value > 0.0 else float("inf"), NONE),
                   "per combination": ("; ".join(f"{c}: {v:.6g} at x = {xx:.6g}" for v, xx, c in rows), TEXT)}
        if details["L / deflection"][0] == float("inf"):
            details["L / deflection"] = ("no deflection", TEXT)
        ratio = value / limit
        out.append(DesignCheck(name=name, clause=CLAUSE, demand=value, capacity=limit, ratio=ratio,
                               passed=bool(ratio <= 1.0), combination=cid, x=x, dimension=LENGTH, advisory=False,
                               details=details, notes=notes, status=OK, method=ctx.method, inputs=inputs))
    return out


__all__ = ["deflection_checks", "deflection_value", "chord_deflection", "member_displacements", "name_of", "CLAUSE"]
