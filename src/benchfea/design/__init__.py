"""Steel member design checks.

Public entry::

    from benchfea import design
    design_results = design.check(model, results)      # -> benchfea.results.DesignResults

Lower-level entry (demands given directly, no analysis needed)::

    ctx = design.DesignContext(member="C1", section=design.DesignSection.i_shape("rolled", ...),
                               method="LRFD", E=..., G=..., Fy=..., Fu=..., L=..., Lx=..., Ly=..., Lz=...)
    member_result = design.check_member(ctx, [design.Demand(combination="LC1", Pc=...,
                                                            extra={"Mrx": ..., "Mry": ...})])

Structure: ``registry`` (codes / editions; unknown code -> refusal), ``base``
(context, demand, check builders), ``section`` (design view of a section),
``demands`` (required strengths from analysis results) and one package per
edition (``aisc360_22``) with one module per chapter.  See
docs/ARCHITECTURE.md, "Contract for benchfea.design".

Orchestration rules (this module):

* DESIGN NEVER CHANGES THE ANALYSIS.  ``check`` only reads ``model`` and
  ``results``; unbraced lengths, K factors, Cb, bracing and every other design
  variable change available strengths only.
* Every design member is checked for every combination with purpose
  "strength" that is in the results: the demands come from the result stations
  (``member_demands``: maxima for Chapters D, E, G; the stations themselves in
  ``Demand.extra`` for Chapters F and H), every limit state of the edition is
  evaluated, every check is kept per combination
  (``DesignMemberResult.by_combination``) and the governing check per limit
  state forms the envelope (``DesignMemberResult.checks``).  Deflection limits
  are checked for their own (serviceability) combinations.
* Status of a member and of the whole design, the first that applies:
  ``"fail"``        a non-advisory check has ratio > 1;
  ``"incomplete"``  an applicable limit state has no result (notChecked /
                    refused with a demand), a strength combination is missing
                    or did not converge, or the demands carry no moments;
  ``"conditional"`` every check passes but the required strengths of at least
                    one strength combination come from a FIRST-ORDER analysis,
                    while AISC 360-22 C1 / C2.1(b) require second-order
                    effects (no B1/B2 amplification is applied);
  ``"pass"``        otherwise.
  The status is never "pass" when anything applicable was not checked.
* Analysis method and K factors.  K = 1 (C3) presumes the direct analysis
  method.  The engine records per combination whether it was analysed
  second-order, with the stiffness reduction of C2.3 and with notional loads
  (C2.2b), and WARNS when members in compression are checked with results that
  do not have all three (or with K < 1); it cannot verify the conditions of
  Appendix 7 and does not change K.
"""
from __future__ import annotations

from dataclasses import replace

from ..results import DesignCheck, DesignMemberResult, DesignResults
from . import demands as _demands
from . import registry
from . import aisc360_22                      # noqa: F401  (registers the edition)
from .base import ASD, LRFD, NOT_CHECKED, OK, REFUSED, DesignContext, Demand
from .demands import demand_from_stations
from .registry import Edition, UnknownDesignCode
from .section import DesignSection

__all__ = ["check", "check_member", "DesignContext", "DesignSection", "Demand", "Edition", "UnknownDesignCode",
           "registry", "member_demands", "demand_from_stations", "analysis_basis", "LRFD", "ASD", "OK",
           "NOT_CHECKED", "REFUSED", "PASS", "CONDITIONAL", "INCOMPLETE", "FAIL"]

PASS, CONDITIONAL, INCOMPLETE, FAIL = "pass", "conditional", "incomplete", "fail"
_RANK = {PASS: 0, CONDITIONAL: 1, INCOMPLETE: 2, FAIL: 3}
_MOMENT_KEYS = ("Mrx", "Mry", "Tr")


def worst_status(statuses) -> str:
    """The least favourable of several statuses ("" entries are ignored; none -> "incomplete")."""
    statuses = [s for s in statuses if s]
    return max(statuses, key=_RANK.__getitem__) if statuses else INCOMPLETE


def member_demands(model, results, design_member_id: str, purpose: str | None = "strength") -> dict:
    """{combination id: Demand} of a design member for the combinations of ``purpose`` (``None``: all)
    that are in ``results``: the maxima of :func:`demand_from_stations` plus, in ``Demand.extra``, the
    stations themselves (``"x"`` from the design member start, ``"forces"`` (n, 6) = N, Vy, Vz, T, My, Mz)
    and the ``"analysisType"`` of the combination -- the data of the flexure and interaction checks."""
    out = {}
    for cid, combo in model.combinations.items():
        if (purpose is None or combo.purpose == purpose) and cid in results.combinations:
            r = results.combinations[cid]
            x, f = _demands.member_stations(model, r, design_member_id)
            out[cid] = replace(demand_from_stations(x, f, cid),
                               extra={"x": x, "forces": f, "analysisType": r.analysis_type})
    return out


def _has_moments(demand: Demand) -> bool:
    e = demand.extra or {}
    return ("x" in e and "forces" in e) or any(k in e for k in _MOMENT_KEYS)


def _governing(candidates: list) -> DesignCheck:
    """Of the checks of one limit state over the combinations: a check without a result for
    a nonzero demand comes first (never hidden by another combination), else the largest ratio."""
    blocked = [c for c in candidates if c.status != OK and (c.demand is None or c.demand > 0.0)]
    if blocked:
        return max(blocked, key=lambda c: -1.0 if c.demand is None else c.demand)
    done = [c for c in candidates if c.status == OK]
    if done:
        return max(done, key=lambda c: c.ratio)
    return candidates[0]


def check_member(ctx: DesignContext, demands, edition: Edition | None = None, extra_checks=(),
                 conditions=(), incomplete=()) -> DesignMemberResult:
    """Check one design member for the demands of several combinations.

    ``demands``: iterable of :class:`Demand` (or a dict id -> Demand).  Every limit state of
    the edition is evaluated for every combination; per limit state the governing combination
    is kept (``checks``), and every check of every combination is kept in ``by_combination``.
    With no demand at all the strengths are still reported (zero demand).

    ``extra_checks``: checks made elsewhere (serviceability) to include in the member result.
    ``conditions`` / ``incomplete``: reasons (text) that make the status at best "conditional" /
    "incomplete" (set by :func:`check` from the analysis basis); they are added to the warnings.
    """
    warnings: list[str] = []
    try:
        edition = edition or registry.get(ctx.code)
    except UnknownDesignCode as exc:
        reason = str(exc)
        refusal = DesignCheck(name="design code", clause="", demand=None, capacity=None, ratio=None, passed=None,
                              status=REFUSED, reason=reason, method=ctx.method,
                              inputs={"code": (ctx.code, "n/a", ctx.source("code"))})
        return DesignMemberResult(ctx.member, [refusal], False, None, [f"design code: {reason}"], [reason],
                                  status=INCOMPLETE)
    demands = list(demands.values() if isinstance(demands, dict) else demands)
    incomplete = list(incomplete)
    if not demands:
        demands = [Demand()]
        warnings.append("no strength combination: strengths are reported for zero demand")
    else:
        bare = [d.combination or "(unnamed)" for d in demands if not _has_moments(d)]
        if bare:
            incomplete.append(f"the demands of combination(s) {', '.join(bare)} carry no moments (Demand.extra): "
                              "flexure (Chapter F), combined forces and torsion (Chapter H) were not evaluated")
    by_name: dict[str, list] = {}
    by_combination: dict[str, list] = {}
    for demand in demands:
        found = edition.check(ctx, demand)
        by_combination.setdefault(demand.combination, []).extend(found)
        for c in found:
            by_name.setdefault(c.name, []).append(c)
    checks = [_governing(cands) for cands in by_name.values()] + list(extra_checks)
    strength = [c for c in checks if not c.advisory]
    failed = [c for c in strength if c.status == OK and not c.passed]
    blocked = [c for c in strength if c.status != OK and (c.demand is None or c.demand > 0.0)]
    done = [c for c in strength if c.status == OK]
    governing = max(done, key=lambda c: c.ratio) if done else None
    not_checked = [f"{c.name}: {c.status}: {c.reason}" for c in checks if c.status != OK]
    for note in ctx.section.notes:
        if note not in warnings:
            warnings.append(note)
    if blocked:
        warnings.append("not all limit states with a nonzero demand could be checked: the member is reported "
                        "as not passed")
    warnings += [w for w in (*incomplete, *conditions) if w not in warnings]
    status = FAIL if failed else INCOMPLETE if (blocked or incomplete) else CONDITIONAL if conditions else PASS
    return DesignMemberResult(ctx.member, checks, passed=not failed and not blocked, governing=governing,
                              not_checked=not_checked, warnings=warnings, status=status,
                              by_combination=by_combination)


# ============================================================================== analysis basis
def analysis_basis(model, results):
    """How each combination was analysed, and what AISC 360-22 Chapter C says about it.

    -> (basis, warnings, conditions, incomplete):
    ``basis[cid]`` = {purpose, inResults, analysisType, secondOrder, converged, stiffnessReduction,
    notionalLoad, alpha, directAnalysisSettings}; ``warnings`` = texts for the design document;
    ``conditions`` / ``incomplete`` = reasons that limit the status of every member (first-order strength
    combination / strength combination missing or not converged)."""
    basis, warnings, conditions, incomplete = {}, [], [], []
    methods = {p.method for p in model.design_members.values()}
    for cid, combo in model.combinations.items():
        r = results.combinations.get(cid)
        second = r is not None and r.analysis_type == "secondOrder"
        row = {"purpose": combo.purpose, "inResults": r is not None,
               "analysisType": r.analysis_type if r is not None else None, "secondOrder": second,
               "converged": bool(r.converged) if r is not None else None,
               "stiffnessReduction": bool(combo.stiffness_reduction), "notionalLoad": combo.notional is not None,
               "alpha": float(r.alpha) if r is not None else None}
        row["directAnalysisSettings"] = bool(second and row["stiffnessReduction"] and row["notionalLoad"])
        basis[cid] = row
        if combo.purpose != "strength":
            continue
        if r is None:
            incomplete.append(f"strength combination {cid!r} is not in the analysis results: not checked")
            continue
        if not r.converged:
            incomplete.append(f"combination {cid!r} did not converge: its demands are not reliable")
        if not second:
            conditions.append(
                f"combination {cid!r} was analysed first-order: AISC 360-22 C1 and C2.1(b) require second-order "
                "effects (P-Delta and P-delta) in the required strengths; the checks use the first-order results "
                "as given (the B1/B2 amplification of Appendix 8 is not applied)")
            continue
        if not row["stiffnessReduction"]:
            warnings.append(
                f"combination {cid!r} was analysed second-order WITHOUT the stiffness reduction of AISC 360-22 "
                "C2.3: this is not a direct analysis; K = 1 is then justified only under Appendix 7.2.3 (braced "
                "systems, or second-order / first-order drift <= 1.1), which the engine does not verify")
        if not row["notionalLoad"]:
            warnings.append(
                f"combination {cid!r} has no notional load (AISC 360-22 C2.2b): initial imperfections must then "
                "be modelled directly in the geometry (C2.2a), which the engine cannot verify")
        if ASD in methods and abs(r.alpha - 1.6) > 1e-9:
            warnings.append(
                f"combination {cid!r}: design method ASD with a second-order analysis at alpha = {r.alpha:g}; "
                "AISC 360-22 C2.1(d) requires the analysis at 1.6 times the ASD loads, results divided by 1.6")
    return basis, warnings, conditions, incomplete


def _k_warnings(ctx: DesignContext, demands: dict, basis: dict) -> list:
    """Effective length factors vs the analysis method, for a member with axial compression."""
    used = [cid for cid, d in demands.items() if d.Pc > 0.0]
    if not used:
        return []
    out = []
    ks = f"Kx = {ctx.Kx:g}, Ky = {ctx.Ky:g}, Kz = {ctx.Kz:g}"
    other = [cid for cid in used if not basis[cid]["directAnalysisSettings"]]
    if other:
        out.append(f"effective length factors {ks}: K = 1 (AISC 360-22 C3) presumes the direct analysis method, "
                   f"but combination(s) {', '.join(other)} with compression in this member were not analysed "
                   "second-order with stiffness reduction (C2.3) and notional loads (C2.2b); with any other "
                   "method the K factors must come from Appendix 7.2.3 (sidesway buckling analysis for columns "
                   "of moment frames) -- not verified by the engine")
    if min(ctx.Kx, ctx.Ky, ctx.Kz) < 1.0:
        out.append(f"effective length factors {ks}: a K below 1 must be justified by rational analysis "
                   "(AISC 360-22 C3, Appendix 7.2.3(a)) -- not verified by the engine")
    return out


def check(model, results) -> DesignResults:
    """Design checks of every design member of ``model`` for the analysis ``results``: all strength
    combinations (demands from the result stations), the deflection limits, the envelope per limit state
    and the status (module docstring).  Reads ``model`` and ``results`` only -- it never changes them.
    Never raises for a member it cannot check: such limit states are ``notChecked`` / ``refused`` with
    the reason."""
    basis, warnings, conditions, incomplete = analysis_basis(model, results)
    warnings = [*incomplete, *conditions, *warnings]
    members: dict[str, DesignMemberResult] = {}
    codes = []
    for dmid in model.design_members:
        ctx = DesignContext.from_model(model, dmid)
        if ctx.code not in codes:
            codes.append(ctx.code)
        demands = member_demands(model, results, dmid)
        service = []
        try:
            for fn in registry.get(ctx.code).service:
                service += fn(model, results, dmid, ctx)
        except UnknownDesignCode:
            pass                                    # check_member reports the refusal
        res = check_member(ctx, demands, extra_checks=service, conditions=conditions, incomplete=incomplete)
        res.warnings += [w for w in _k_warnings(ctx, demands, basis) if w not in res.warnings]
        members[dmid] = res
    if not members:
        warnings.append("the model has no design member: nothing was checked")
    return DesignResults(code=", ".join(codes), members=members, defaults_used=list(model.defaults_used),
                         warnings=warnings, status=worst_status(m.status for m in members.values()),
                         combinations=basis)
