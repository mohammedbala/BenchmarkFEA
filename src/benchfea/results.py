"""Result data structures shared by the solver, the design checks, the
verification suite and the benchmark environment, and the writer of the
output document.

Everything in these dataclasses is in INTERNAL units (N, mm, rad).  Only
``Results.to_document`` / ``DesignResults.to_document`` convert to the units
of the input file.

Conventions (docs/CONVENTIONS.md is authoritative):

* nodal arrays have shape (n_nodes, 6), rows in ``Model.nodes`` order, columns
  ``ux, uy, uz, rx, ry, rz`` (or ``Fx, Fy, Fz, Mx, My, Mz``), GLOBAL axes;
* ``reactions`` are forces exerted BY the supports ON the structure (spring
  supports included), zero rows for unsupported nodes;
* ``MemberResult.end_forces`` (12,) = ``Fx, Fy, Fz, Mx, My, Mz`` at end i then
  end j: forces exerted by the joints ON the member at the ends of its
  flexible length, member LOCAL axes;
* ``MemberResult.forces`` (n, 6) = ``N, Vy, Vz, T, My, Mz`` at the stations:
  stress resultants on the cross-section face with outward normal +x.  Hence
  ``forces[-1] == end_forces[6:]`` and ``forces[0] == -end_forces[:6]``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

OUTPUT_FORMAT = "benchfea-results"
OUTPUT_VERSION = "1.0.0"
DESIGN_OUTPUT_FORMAT = "benchfea-design"
DESIGN_OUTPUT_VERSION = "1.0.0"

NODE_DISPLACEMENTS = ("ux", "uy", "uz", "rx", "ry", "rz")
NODE_FORCES = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")
STATION_FORCES = ("N", "Vy", "Vz", "T", "My", "Mz")
_DISP_DIMS = ("length",) * 3 + ("rotation",) * 3
_FORCE_DIMS = ("force",) * 3 + ("moment",) * 3
_STATION_DIMS = ("force", "force", "force", "moment", "moment", "moment")


@dataclass
class MemberResult:
    end_forces: np.ndarray            # (12,) local, on the member, end i then end j
    x: np.ndarray                     # (n,) station positions from the flexible start, mm
    forces: np.ndarray                # (n, 6) N, Vy, Vz, T, My, Mz
    displacements: np.ndarray         # (n, 3) total displacement of the member axis, LOCAL ux, uy, uz
    active: bool = True               # False: tension/compression-only member switched off


@dataclass
class InvariantCheck:
    name: str                         # e.g. "global force equilibrium"
    value: float                      # measured (relative) error
    tolerance: float
    passed: bool
    detail: str = ""


@dataclass
class InvariantsReport:
    """Self-checks every analysis must report (the solver fills these in):
    global force and moment equilibrium (applied loads + reactions), joint
    equilibrium residual, and for linear analyses external work = strain energy."""
    checks: list[InvariantCheck] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks)


@dataclass
class CombinationResult:
    id: str
    analysis_type: str                # "firstOrder" | "secondOrder"
    displacements: np.ndarray         # (n_nodes, 6) global
    reactions: np.ndarray             # (n_nodes, 6) global
    members: dict[str, MemberResult]
    converged: bool = True
    iterations: int = 1
    alpha: float = 1.0                # loads were multiplied and results divided by this (already divided here)
    tau_b: dict[str, float] = field(default_factory=dict)     # member id -> tau_b used (stiffness reduction)
    invariants: InvariantsReport = field(default_factory=InvariantsReport)
    warnings: list[str] = field(default_factory=list)


@dataclass
class BucklingResult:
    combination: str
    factors: np.ndarray               # (n_modes,) load multipliers, ascending, positive
    modes: np.ndarray                 # (n_modes, n_nodes, 6) mode shapes, max |translation| = 1


@dataclass
class Envelope:
    """Max/min of the station forces of one member over all combinations."""
    x: np.ndarray                     # (n,)
    max: np.ndarray                   # (n, 6)
    min: np.ndarray                   # (n, 6)
    max_combination: list[list[str]]  # (n, 6) governing combination ids
    min_combination: list[list[str]]


@dataclass
class Results:
    model_name: str
    combinations: dict[str, CombinationResult]
    buckling: dict[str, BucklingResult] = field(default_factory=dict)
    envelopes: dict[str, Envelope] = field(default_factory=dict)
    defaults_used: list = field(default_factory=list)          # list[model.DefaultUsed]
    warnings: list[str] = field(default_factory=list)
    engine_version: str = ""

    def compute_envelopes(self) -> None:
        """Fill ``envelopes`` from ``combinations``."""
        self.envelopes = {}
        combos = list(self.combinations.values())
        if not combos:
            return
        for mid in combos[0].members:
            stack = np.array([c.members[mid].forces for c in combos])         # (n_combo, n, 6)
            imax, imin = stack.argmax(axis=0), stack.argmin(axis=0)
            ids = [c.id for c in combos]
            self.envelopes[mid] = Envelope(
                combos[0].members[mid].x.copy(), stack.max(axis=0), stack.min(axis=0),
                [[ids[k] for k in row] for row in imax], [[ids[k] for k in row] for row in imin])

    def to_document(self, model) -> dict:
        """The output document (format ``benchfea-results`` 1.0.0) in the units
        of the input file. Plain JSON-serialisable data."""
        u = model.units
        nodes = list(model.nodes)

        def table(arr, names, dims):
            return {nid: {n: float(u.from_internal(arr[k, c], d)) for c, (n, d) in enumerate(zip(names, dims))}
                    for k, nid in enumerate(nodes)}

        def cols(arr, dims):
            return [[float(u.from_internal(v, d)) for v, d in zip(row, dims)] for row in arr]

        out: dict[str, Any] = {
            "format": OUTPUT_FORMAT, "formatVersion": OUTPUT_VERSION, "engineVersion": self.engine_version,
            "model": self.model_name, "units": u.to_dict(),
            "unitLabels": {d: u.label(d) for d in ("length", "rotation", "force", "moment")},
            "conventions": {"verticalAxis": model.vertical_axis},
            "combinations": {}, "buckling": {}, "envelopes": {},
            "defaultsUsed": [x.to_dict() for x in self.defaults_used],
            "warnings": list(self.warnings),
        }
        for cid, c in self.combinations.items():
            supported = [k for k, nid in enumerate(nodes) if nid in model.supports]
            reactions = table(c.reactions, NODE_FORCES, _FORCE_DIMS)
            out["combinations"][cid] = {
                "analysisType": c.analysis_type, "converged": c.converged, "iterations": c.iterations,
                "alpha": c.alpha, "tauB": dict(c.tau_b),
                "displacements": table(c.displacements, NODE_DISPLACEMENTS, _DISP_DIMS),
                "reactions": {nodes[k]: reactions[nodes[k]] for k in supported},
                "members": {
                    mid: {
                        "active": m.active,
                        "endForces": {
                            "i": dict(zip(NODE_FORCES, cols([m.end_forces[:6]], _FORCE_DIMS)[0])),
                            "j": dict(zip(NODE_FORCES, cols([m.end_forces[6:]], _FORCE_DIMS)[0]))},
                        "stations": {
                            "x": [float(u.from_internal(v, "length")) for v in m.x],
                            **{n: [r[k] for r in cols(m.forces, _STATION_DIMS)] for k, n in enumerate(STATION_FORCES)},
                            **{n: [r[k] for r in cols(m.displacements, ("length",) * 3)]
                               for k, n in enumerate(("ux", "uy", "uz"))}},
                    } for mid, m in c.members.items()},
                "invariants": {"passed": c.invariants.passed,
                               "checks": [vars(ch) for ch in c.invariants.checks]},
                "warnings": list(c.warnings),
            }
        for cid, b in self.buckling.items():
            out["buckling"][cid] = {
                "factors": [float(f) for f in b.factors],
                "modes": [table(mode, NODE_DISPLACEMENTS, ("dimensionless",) * 6) for mode in b.modes]}
        for mid, e in self.envelopes.items():
            out["envelopes"][mid] = {
                "x": [float(u.from_internal(v, "length")) for v in e.x],
                "max": {n: [r[k] for r in cols(e.max, _STATION_DIMS)] for k, n in enumerate(STATION_FORCES)},
                "min": {n: [r[k] for r in cols(e.min, _STATION_DIMS)] for k, n in enumerate(STATION_FORCES)},
                "maxCombination": {n: [r[k] for r in e.max_combination] for k, n in enumerate(STATION_FORCES)},
                "minCombination": {n: [r[k] for r in e.min_combination] for k, n in enumerate(STATION_FORCES)}}
        return out


# ------------------------------------------------------------------------------ design
@dataclass
class DesignCheck:
    """One limit-state check of one design member."""
    name: str                         # e.g. "flexure, major axis: lateral-torsional buckling"
    clause: str                       # citation, e.g. "AISC 360-22 F2.2, Eq. F2-2"
    demand: float                     # required strength (internal units) or measured quantity
    capacity: float                   # available strength (phi*Rn or Rn/Omega) or limit
    ratio: float                      # demand / capacity
    passed: bool
    combination: str = ""             # governing combination id
    x: float | None = None            # governing location from the design member start, mm
    dimension: str = "dimensionless"  # dimension of demand/capacity for unit conversion
    advisory: bool = False            # True: does not affect `passed` of the member (e.g. slenderness user notes)
    details: dict = field(default_factory=dict)   # intermediate values {name: (value, dimension)}
    notes: list[str] = field(default_factory=list)


@dataclass
class DesignMemberResult:
    id: str
    checks: list[DesignCheck]
    passed: bool
    governing: DesignCheck | None = None
    not_checked: list[str] = field(default_factory=list)   # limit states outside the implemented scope, with reason
    warnings: list[str] = field(default_factory=list)


@dataclass
class DesignResults:
    code: str
    members: dict[str, DesignMemberResult]
    defaults_used: list = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_document(self, model) -> dict:
        """The design output document (format ``benchfea-design`` 1.0.0) in file units."""
        u = model.units

        def conv(v, dim):
            return None if v is None else float(u.from_internal(v, dim))

        def check(c: DesignCheck) -> dict:
            return {"name": c.name, "clause": c.clause, "demand": conv(c.demand, c.dimension),
                    "capacity": conv(c.capacity, c.dimension), "unit": u.label(c.dimension), "ratio": float(c.ratio),
                    "passed": c.passed, "advisory": c.advisory, "combination": c.combination,
                    "x": conv(c.x, "length"),
                    "details": {k: {"value": conv(v, d) if isinstance(v, (int, float)) else v, "unit": u.label(d)}
                                for k, (v, d) in c.details.items()},
                    "notes": list(c.notes)}

        return {
            "format": DESIGN_OUTPUT_FORMAT, "formatVersion": DESIGN_OUTPUT_VERSION, "code": self.code,
            "model": model.name, "units": u.to_dict(),
            "members": {
                mid: {"passed": m.passed, "governing": check(m.governing) if m.governing else None,
                      "checks": [check(c) for c in m.checks], "notChecked": list(m.not_checked),
                      "parameters": {n: {"value": _plain(r.value), "source": r.source, "path": r.path}
                                     for n, r in model.design_members[mid].values.items()},
                      "warnings": list(m.warnings)}
                for mid, m in self.members.items()},
            "defaultsUsed": [x.to_dict() for x in self.defaults_used],
            "warnings": list(self.warnings),
        }


def _plain(v):
    if isinstance(v, tuple):
        return [_plain(x) for x in v]
    if hasattr(v, "__dataclass_fields__"):
        return {k: _plain(getattr(v, k)) for k in v.__dataclass_fields__}
    return v
