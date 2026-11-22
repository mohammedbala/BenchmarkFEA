"""Steel member design checks.

Public entry::

    from benchfea import design
    design_results = design.check(model, results)      # -> benchfea.results.DesignResults

Lower-level entry (demands given directly, no analysis needed)::

    ctx = design.DesignContext(member="C1", section=design.DesignSection.i_shape("rolled", ...),
                               method="LRFD", E=..., G=..., Fy=..., Fu=..., L=..., Lx=..., Ly=..., Lz=...)
    member_result = design.check_member(ctx, [design.Demand(combination="LC1", Pc=...)])

Structure: ``registry`` (codes / editions; unknown code -> refusal), ``base``
(context, demand, check builders), ``section`` (design view of a section),
``demands`` (required strengths from analysis results) and one package per
edition (``aisc360_22``) with one module per chapter.  See
docs/ARCHITECTURE.md, "Design checks".
"""
from __future__ import annotations

from ..results import DesignCheck, DesignMemberResult, DesignResults
from . import registry
from . import aisc360_22                      # noqa: F401  (registers the edition)
from .base import ASD, LRFD, NOT_CHECKED, OK, REFUSED, DesignContext, Demand
from .demands import demand_from_stations, member_demands
from .registry import Edition, UnknownDesignCode
from .section import DesignSection

__all__ = ["check", "check_member", "DesignContext", "DesignSection", "Demand", "Edition", "UnknownDesignCode",
           "registry", "member_demands", "demand_from_stations", "LRFD", "ASD", "OK", "NOT_CHECKED", "REFUSED"]


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


def check_member(ctx: DesignContext, demands, edition: Edition | None = None) -> DesignMemberResult:
    """Check one design member for the demands of several combinations.

    ``demands``: iterable of :class:`Demand` (or a dict id -> Demand).  Every limit state of
    the edition is evaluated for every combination; per limit state the governing combination
    is kept.  With no demand at all the strengths are still reported (zero demand).
    """
    warnings: list[str] = []
    try:
        edition = edition or registry.get(ctx.code)
    except UnknownDesignCode as exc:
        reason = str(exc)
        refusal = DesignCheck(name="design code", clause="", demand=None, capacity=None, ratio=None, passed=None,
                              status=REFUSED, reason=reason, method=ctx.method,
                              inputs={"code": (ctx.code, "n/a", ctx.source("code"))})
        return DesignMemberResult(ctx.member, [refusal], False, None, [f"design code: {reason}"], [reason])
    demands = list(demands.values() if isinstance(demands, dict) else demands)
    if not demands:
        demands = [Demand()]
        warnings.append("no strength combination: strengths are reported for zero demand")
    by_name: dict[str, list] = {}
    for demand in demands:
        for c in edition.check(ctx, demand):
            by_name.setdefault(c.name, []).append(c)
    checks = [_governing(cands) for cands in by_name.values()]
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
    return DesignMemberResult(ctx.member, checks, passed=not failed and not blocked, governing=governing,
                              not_checked=not_checked, warnings=warnings)


def check(model, results) -> DesignResults:
    """Design checks of every design member of ``model`` for the strength combinations in
    ``results`` (demands = maxima over the result stations).  Never raises for a member it
    cannot check: such limit states are reported as ``notChecked`` / ``refused`` with the reason."""
    members: dict[str, DesignMemberResult] = {}
    warnings: list[str] = []
    codes = []
    for cid, combo in model.combinations.items():
        if combo.purpose == "strength" and cid in results.combinations and not results.combinations[cid].converged:
            warnings.append(f"combination {cid!r} did not converge: its demands are not reliable")
    for dmid in model.design_members:
        ctx = DesignContext.from_model(model, dmid)
        if ctx.code not in codes:
            codes.append(ctx.code)
        members[dmid] = check_member(ctx, member_demands(model, results, dmid))
    return DesignResults(code=", ".join(codes), members=members, defaults_used=list(model.defaults_used),
                         warnings=warnings)
