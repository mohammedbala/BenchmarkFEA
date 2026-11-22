"""Loads of a load case / combination in solver form.

A :class:`LoadSet` is linear in the loads, so a combination is ``sum(factor * case)``.
Member loads stay member loads (local axes, flexible length) until the element turns them into
fixed-end forces; they are never lumped.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..model import MemberDistributedLoad, MemberPointLoad, NodalLoad, SelfWeight, SupportDisplacement
from ..problems import pointer
from .beam import SpanLoads
from .element import TRUSS_TYPES
from .errors import refuse

__all__ = ["LoadSet", "case_loads", "combine", "combination_loads", "applied_resultant", "member_load_positions"]


@dataclass
class LoadSet:
    nodal: np.ndarray                                   # (n_dof,) loads given at nodes, global
    prescribed: dict[int, float] = field(default_factory=dict)      # global DOF -> imposed displacement
    member: dict[str, SpanLoads] = field(default_factory=dict)      # member id -> span loads (local)


def _span(load, model):
    if isinstance(load, MemberPointLoad):
        f = np.zeros(6)
        if load.kind == "force":
            f[:3] = load.vector
        else:
            f[3:] = load.vector
        return load.member, SpanLoads(points=((float(load.x), f),))
    return load.member, SpanLoads(dists=((float(load.x1), float(load.x2), np.asarray(load.w1, float),
                                          np.asarray(load.w2, float)),))


def case_loads(model, case) -> LoadSet:
    """Loads of one load case (unfactored)."""
    out = LoadSet(np.zeros(model.n_dof))
    terms: dict[str, list] = {}
    for k, ld in enumerate(case.loads):
        if isinstance(ld, NodalLoad):
            i = 6 * model.nodes[ld.node].index
            out.nodal[i:i + 6] += ld.values
        elif isinstance(ld, SupportDisplacement):
            for dof, v in ld.values.items():
                g = model.dof(ld.node, dof)
                out.prescribed[g] = out.prescribed.get(g, 0.0) + float(v)
        elif isinstance(ld, SelfWeight):
            for mid in ld.members:
                m = model.members[mid]
                w = ld.factor * model.materials[m.material].unit_weight * model.sections[m.section].A
                if w != 0.0:
                    wl = w * (np.asarray(m.R, float) @ np.asarray(ld.direction, float))
                    terms.setdefault(mid, []).append(SpanLoads(dists=((0.0, float(m.length), wl, wl),)))
        elif isinstance(ld, (MemberPointLoad, MemberDistributedLoad)):
            mid, span = _span(ld, model)
            if model.members[mid].type in TRUSS_TYPES and isinstance(ld, MemberPointLoad) \
                    and ld.kind == "moment" and abs(ld.vector[0]) > 0.0:
                raise refuse("load.torque-on-truss", f"{pointer('loadCases', case.id)}/loads/{k}",
                             f"a moment about the axis of the {model.members[mid].type} member {mid!r} cannot be "
                             "carried: a truss-type member has no torsional end connection",
                             "Apply the torque to a frame member or to a node.")
            terms.setdefault(mid, []).append(span)
        else:                                            # pragma: no cover - guarded by the model layer
            raise TypeError(f"unknown load {ld!r}")
    out.member = {mid: SpanLoads.combine((1.0, s) for s in spans) for mid, spans in terms.items()}
    return out


def combine(model, terms) -> LoadSet:
    """``sum(factor * LoadSet)`` over ``terms`` = iterable of (factor, LoadSet)."""
    terms = list(terms)
    out = LoadSet(np.zeros(model.n_dof))
    spans: dict[str, list] = {}
    for f, ls in terms:
        out.nodal += f * ls.nodal
        for g, v in ls.prescribed.items():
            out.prescribed[g] = out.prescribed.get(g, 0.0) + f * v
        for mid, s in ls.member.items():
            spans.setdefault(mid, []).append((f, s))
    out.member = {mid: SpanLoads.combine(t) for mid, t in spans.items()}
    return out


def combination_loads(model, combination, cases: dict[str, LoadSet] | None = None) -> LoadSet:
    """Factored loads of a combination (``cases`` = cache of :func:`case_loads` by case id)."""
    if cases is None:
        cases = {}
    terms = []
    for cid, f in combination.factors.items():
        if cid not in cases:
            cases[cid] = case_loads(model, model.load_cases[cid])
        terms.append((float(f), cases[cid]))
    return combine(model, terms)


def member_load_positions(model) -> dict[str, tuple[set, set]]:
    """For every member: (positions of point loads, ends of line loads) over ALL load cases, so that the
    result stations of a member are the same in every combination."""
    out = {mid: (set(), set()) for mid in model.members}
    for case in model.load_cases.values():
        for ld in case.loads:
            if isinstance(ld, MemberPointLoad):
                out[ld.member][0].add(float(ld.x))
            elif isinstance(ld, MemberDistributedLoad):
                out[ld.member][1].update((float(ld.x1), float(ld.x2)))
    return out


def applied_resultant(model, loads: LoadSet):
    """Resultant of the applied loads computed directly from their definition (not from equivalent
    nodal loads): ``(R, A)`` with R = (Fx, Fy, Fz, Mx, My, Mz) about the global origin and A the sums
    of the absolute values of the contributions (scales for relative errors)."""
    R, A = np.zeros(6), np.zeros(6)

    def add(point, force, moment):
        m = np.cross(point, force)
        R[:3] += force
        R[3:] += moment + m
        A[:3] += np.abs(force)
        A[3:] += np.abs(moment) + np.abs(m)

    for node in model.nodes.values():
        v = loads.nodal[6 * node.index:6 * node.index + 6]
        if v.any():
            add(node.xyz, v[:3], v[3:])
    for mid, span in loads.member.items():
        m = model.members[mid]
        Rt = np.asarray(m.R, float).T
        ex = Rt[:, 0]
        for x, f in span.points:
            add(m.start + x * ex, Rt @ f[:3], Rt @ f[3:])
        for x1, x2, w1, w2 in span.dists:
            h = x2 - x1
            for w, xc in ((w1, x1 + h / 3.0), (w2, x1 + 2.0 * h / 3.0)):
                add(m.start + xc * ex, Rt @ (w * (h / 2.0)), np.zeros(3))
    return R, A
