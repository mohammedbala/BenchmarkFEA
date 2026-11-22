"""Solver-neutral description of ONE combination of a benchfea model: what every translator starts from.

``prepare(model, combination_id)`` resolves, in internal units (N, mm, rad) and benchfea conventions:

* factored loads: nodal loads (global), imposed support displacements, member point loads and line
  loads in member local axes (self weight turned into uniform line loads
  ``factor * unit_weight * A * (R @ direction)``), all multiplied by ``alpha`` for a second-order
  combination (the mapped results are divided by ``alpha`` again, AISC 360-22 C2.1(d) as stated in
  docs/ARCHITECTURE.md);
* member stiffnesses after the stiffness reduction of the combination (EA, EI, GJ, GAs times
  ``reduction_factor`` for the reduction members, EI additionally times tau_b);
* the degrees of freedom that benchfea leaves out of the solution and reports as zero (rotations of
  nodes joined only by truss-type members, free DOFs of nodes without members);
* features no translator may guess (``issues``): notional loads, a "computed" tau_b, tension-/
  compression-only members.

Nothing here imports a foreign solver or ``benchfea.solve``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..model import (DOF_INDEX, MemberDistributedLoad, MemberPointLoad, Model, NodalLoad, SelfWeight,
                     SupportDisplacement)
from .fidelity import EXACT, UNSUPPORTED

TRUSS_TYPES = ("truss", "tensionOnly", "compressionOnly")

__all__ = ["Prepared", "PMember", "prepare", "TRUSS_TYPES"]


@dataclass
class PMember:
    id: str
    type: str
    i: str
    j: str
    E: float
    G: float
    A: float
    Iy: float
    Iz: float
    J: float
    GAy: float | None                 # shear stiffness for shear along local y (None = no shear deformation)
    GAz: float | None
    EIy: float                        # flexural stiffnesses incl. tau_b
    EIz: float
    R: np.ndarray
    start: np.ndarray
    end: np.ndarray
    L: float
    offset_i: np.ndarray
    offset_j: np.ndarray
    k_i: tuple
    k_j: tuple
    points: list = field(default_factory=list)     # (kind, vector(3) local, x)
    dists: list = field(default_factory=list)      # (w1(3), w2(3), x1, x2) local

    @property
    def truss(self) -> bool:
        return self.type in TRUSS_TYPES

    @property
    def loaded(self) -> bool:
        return bool(self.points or self.dists)


@dataclass
class Prepared:
    model: Model
    combination: str
    second_order: bool
    alpha: float
    shear: bool
    node_ids: list
    xyz: np.ndarray                   # (n, 3)
    fixed: np.ndarray                 # (n, 6) bool
    spring: np.ndarray                # (n, 6)
    prescribed: np.ndarray            # (n, 6) imposed displacement of fixed DOFs (times alpha)
    load: np.ndarray                  # (n, 6) nodal loads (times alpha)
    inactive: dict                    # node id -> DOF indices without stiffness that benchfea reports as 0
                                      # (rotations of truss-only nodes, free DOFs of nodes without members)
    members: dict
    constraints: list
    issues: list                      # (feature, level, note) valid for every solver

    def index(self, node_id: str) -> int:
        return self.model.nodes[node_id].index


def prepare(model: Model, combination_id: str, tau_b: dict | None = None,
            extra_nodal_loads: dict | None = None) -> Prepared:
    """``tau_b``: member id -> tau_b to use when the model asks for a "computed" tau_b (it depends on the
    axial forces of the analysis being checked, so it has to be supplied, e.g. from
    ``CombinationResult.tau_b``).  ``extra_nodal_loads``: node id -> 6 global load components added
    (unfactored by alpha) to the combination, e.g. the notional loads of the analysis being checked."""
    combo = model.combinations[combination_id]
    an = model.analysis
    second = combo.analysis_type == "secondOrder"
    alpha = float(an.alpha) if second else 1.0
    n = len(model.nodes)
    ids = list(model.nodes)
    issues = []

    fixed = np.zeros((n, 6), bool)
    spring = np.zeros((n, 6))
    for nid, s in model.supports.items():
        k = model.nodes[nid].index
        fixed[k] = s.fixed
        spring[k] = s.spring
    load = np.zeros((n, 6))
    prescribed = np.zeros((n, 6))

    # ---- members, stiffness reduction
    members: dict[str, PMember] = {}
    reduce = set(model.members if an.reduction_members is None else an.reduction_members) \
        if combo.stiffness_reduction else set()
    for mid, m in model.members.items():
        mat, sec = model.materials[m.material], model.sections[m.section]
        f = float(an.reduction_factor) if mid in reduce else 1.0
        tb = 1.0
        if mid in reduce:
            if an.tau_b == "computed":
                if tau_b is not None and mid in tau_b:
                    tb = float(tau_b[mid])
                    issues.append(("stiffnessReduction", EXACT,
                                   "tau_b values supplied by the caller (not recomputed by the foreign solver)"))
                else:
                    issues.append(("stiffnessReduction", UNSUPPORTED,
                                   "tau_b = 'computed' depends on the axial forces of the analysis (AISC 360-22 "
                                   "Eq. C2-2a/b); pass the values used as tau_b={member: value}"))
            elif an.tau_b != "none":
                tb = float(an.tau_b)
        E, G = mat.E * f, mat.G * f
        shear = an.shear_deformation and not m.type in TRUSS_TYPES
        members[mid] = PMember(
            mid, m.type, m.i, m.j, E, G, float(sec.A), float(sec.Iy), float(sec.Iz), float(sec.J),
            G * sec.Asy if shear and sec.Asy else None, G * sec.Asz if shear and sec.Asz else None,
            E * sec.Iy * tb, E * sec.Iz * tb, np.asarray(m.R, float), np.asarray(m.start, float),
            np.asarray(m.end, float), float(m.length), np.asarray(m.offset_i, float),
            np.asarray(m.offset_j, float), tuple(m.k_end_i), tuple(m.k_end_j))
        if m.type in ("tensionOnly", "compressionOnly"):
            issues.append(("memberTypes." + m.type, UNSUPPORTED,
                           "one-way members need an activation loop that the translators do not implement"))
    if reduce:
        issues.append(("stiffnessReduction", EXACT, "stiffnesses scaled in the translator"))

    # ---- loads
    for cid, factor in combo.factors.items():
        fac = float(factor) * alpha
        for ld in model.load_cases[cid].loads:
            if isinstance(ld, NodalLoad):
                load[model.nodes[ld.node].index] += fac * ld.values
            elif isinstance(ld, SupportDisplacement):
                for dof, v in ld.values.items():
                    prescribed[model.nodes[ld.node].index, dof] += fac * float(v)
            elif isinstance(ld, SelfWeight):
                for mid in ld.members:
                    m = members[mid]
                    w = ld.factor * model.materials[model.members[mid].material].unit_weight * m.A
                    if w != 0.0:
                        wl = fac * w * (m.R @ np.asarray(ld.direction, float))
                        m.dists.append((wl, wl.copy(), 0.0, m.L))
            elif isinstance(ld, MemberPointLoad):
                members[ld.member].points.append((ld.kind, fac * np.asarray(ld.vector, float), float(ld.x)))
            elif isinstance(ld, MemberDistributedLoad):
                if ld.x2 > ld.x1:
                    members[ld.member].dists.append((fac * np.asarray(ld.w1, float), fac * np.asarray(ld.w2, float),
                                                     float(ld.x1), float(ld.x2)))
            else:                                   # pragma: no cover
                raise TypeError(f"unknown load {ld!r}")
    for m in members.values():
        if m.truss and any(kind == "moment" and abs(v[0]) > 0.0 for kind, v, _ in m.points):
            issues.append(("loads.memberPoint", UNSUPPORTED,
                           f"torque on the truss-type member {m.id!r} (benchfea refuses it: load.torque-on-truss)"))
    if combo.notional is not None:
        if extra_nodal_loads is None:
            issues.append(("notionalLoad", UNSUPPORTED,
                           "the notional loads are derived quantities of the analysis being checked (sign for "
                           "uplift and the lumping of member loads are engine decisions); pass them as "
                           "extra_nodal_loads={node: [Fx, Fy, Fz, Mx, My, Mz]}"))
        else:
            issues.append(("notionalLoad", EXACT, "notional loads supplied by the caller as nodal loads"))
    for nid, v in (extra_nodal_loads or {}).items():
        load[model.nodes[nid].index] += alpha * np.asarray(v, float)

    # ---- rotations without stiffness (benchfea reports them as 0): nodes joined only by truss-type
    #      members whose rigid offset at that node is zero, that are not constraint masters
    attached: dict[str, list] = {nid: [] for nid in ids}
    for m in members.values():
        attached[m.i].append((m, m.offset_i))
        attached[m.j].append((m, m.offset_j))
    masters = {c.master for c in model.constraints.values()}
    slave_dofs: dict[str, set] = {}
    for c in model.constraints.values():
        for sl in c.slaves:
            slave_dofs.setdefault(sl, set()).update(DOF_INDEX[d] for d in c.dofs)
    drive_dofs: dict[str, set] = {}
    for c in model.constraints.values():
        drive_dofs.setdefault(c.master, set()).update(DOF_INDEX[d] for d in c.dofs)
    inactive: dict[str, list] = {}
    for nid in ids:
        k = model.nodes[nid].index
        free = [d for d in range(6) if not fixed[k, d] and spring[k, d] == 0.0 and d not in slave_dofs.get(nid, ())]
        if not attached[nid]:
            # a node without members: everything that does not drive the slaves of a constraint
            dofs = [d for d in free if d not in drive_dofs.get(nid, ())]
        elif nid not in masters and all(m.truss and not off.any() for m, off in attached[nid]):
            dofs = [d for d in free if d >= 3]
        else:
            dofs = []
        if dofs:
            inactive[nid] = dofs
            if np.any(load[k, dofs] != 0.0):
                issues.append(("loads.nodal", UNSUPPORTED,
                               f"load on a degree of freedom of node {nid!r} that has no stiffness"))

    return Prepared(model, combination_id, second, alpha, bool(an.shear_deformation), ids,
                    np.array([nd.xyz for nd in model.nodes.values()], float), fixed, spring, prescribed, load,
                    inactive, members, list(model.constraints.values()), issues)
