"""``analyze(model) -> Results``: the driver.

Structure of the code (and where a nonlinear driver plugs in)
------------------------------------------------------------
``LinearAnalysis``  builds the elements (``build_element``), the ``Structure`` (constraints, supports,
                    assembly), factorises the elastic stiffness once and solves any ``LoadSet``:
                    ``solve_loads(loads) -> (Solution, {member: MemberRecovery})``.
``ANALYSIS_TYPES``  maps ``Combination.analysis_type`` to a function
                    ``(LinearAnalysis, Combination) -> CombinationResult``.  Only ``"firstOrder"`` is
                    registered here; every other type is refused with ``analysis.not-implemented``.
                    A second-order driver registers itself and reuses ``Structure.assemble(matrices)``
                    (geometric stiffness = one more set of per-member 12x12 matrices),
                    ``Structure.factorize`` and ``Structure.solve`` inside its iteration, and
                    ``build_element(model, member, factors, active)`` for reduced stiffness and for
                    tension/compression-only members that are switched off.
``UNSUPPORTED``     the list of model features this linear core refuses (one function each).

This module never reads the design-only data of the model.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from .. import __version__
from ..problems import Problem, pointer
from ..results import CombinationResult, InvariantCheck, InvariantsReport, MemberResult, Results
from . import invariants as inv
from .element import MemberElement, build_element
from .errors import NOT_IMPLEMENTED, ModelRefused
from .loads import LoadSet, applied_resultant, case_loads, combine, member_load_positions
from .system import Structure

__all__ = ["analyze", "LinearAnalysis", "ANALYSIS_TYPES", "station_layout", "unsupported_features"]


def station_layout(length: float, n: int, point_positions=(), line_ends=()) -> list[tuple[float, int]]:
    """Result stations of a member as (x, side): ``n`` equally spaced stations, plus the ends of line
    loads, plus TWO stations at every interior point load (side -1 = just before, +1 = just after), so
    that the jump of the diagram is reported.  x = 0 is taken just before the start (side -1) and
    x = L just after the end (side +1): the values there are minus the end forces at i and the end
    forces at j."""
    L = float(length)
    base = [L * k / (n - 1) for k in range(n)]
    base[-1] = L
    pts = {float(x) for x in point_positions if 0.0 < x < L}
    xs = sorted(set(base) | pts | {float(x) for x in line_ends if 0.0 < x < L})
    out = []
    for x in xs:
        if x in pts:
            out += [(x, -1), (x, +1)]
        else:
            out.append((x, +1 if x == L else -1))
    return out


def unsupported_features(model) -> list[Problem]:
    """Model features that this (linear) solver core does not analyse.  Each is refused rather than
    ignored."""
    out = []

    def add(path, what, hint="This build analyses first-order (linear) combinations of frame and truss members."):
        out.append(Problem(path, what + " is not implemented in this build", hint, code=NOT_IMPLEMENTED))

    for mid, m in model.members.items():
        if m.type in ("tensionOnly", "compressionOnly"):
            add(pointer("members", mid, "type"), f"member type {m.type!r} (member {mid!r})")
    for cid, c in model.combinations.items():
        if c.analysis_type not in ANALYSIS_TYPES:
            add(pointer("combinations", cid, "analysisType"),
                f"analysis type {c.analysis_type!r} (combination {cid!r})")
        if c.stiffness_reduction:
            add(pointer("combinations", cid, "stiffnessReduction"), f"stiffness reduction (combination {cid!r})")
        if c.notional is not None:
            add(pointer("combinations", cid, "notionalLoad"), f"notional loads (combination {cid!r})")
    if model.analysis.buckling_combinations:
        add("/analysis/buckling", "linear buckling analysis")
    return out


class LinearAnalysis:
    """Elements, structure and the factorised elastic stiffness of a model; solves load sets."""

    def __init__(self, model, elements: dict[str, MemberElement] | None = None):
        self.model = model
        self.elements = elements if elements is not None else {mid: build_element(model, m)
                                                               for mid, m in model.members.items()}
        self.structure = Structure(model, self.elements)
        self.K_members = self.structure.assemble()
        self.K = self.K_members.copy()
        self.K[np.diag_indices(self.structure.n)] += self.structure.spring
        self._factor = None
        self._cases: dict[str, LoadSet] = {}
        pos = member_load_positions(model)
        n = int(model.analysis.result_stations)
        self.stations = {mid: station_layout(m.length, n, *pos[mid]) for mid, m in model.members.items()}
        self._model_checks = None

    @property
    def factor(self):
        """Factorisation of the elastic stiffness (raises ModelUnstable with a diagnosis)."""
        if self._factor is None:
            self._factor = self.structure.factorize(self.K)
        return self._factor

    @property
    def model_checks(self) -> list[InvariantCheck]:
        if self._model_checks is None:
            self._model_checks = inv.model_checks(self.structure, self.K_members, self.factor)
        return self._model_checks

    def case(self, case_id: str) -> LoadSet:
        if case_id not in self._cases:
            self._cases[case_id] = case_loads(self.model, self.model.load_cases[case_id])
        return self._cases[case_id]

    def combination_loads(self, combination) -> LoadSet:
        return combine(self.model, [(float(f), self.case(cid)) for cid, f in combination.factors.items()])

    def load_vector(self, loads: LoadSet) -> np.ndarray:
        """Nodal loads + equivalent nodal loads of the member loads (n,)."""
        F = loads.nodal.copy()
        for mid, span in loads.member.items():
            e = self.elements[mid]
            F[e.dofs] += e.equivalent_nodal_loads(span)
        return F

    def recover(self, u: np.ndarray, loads: LoadSet) -> dict:
        from .beam import SpanLoads
        empty = SpanLoads()
        return {mid: e.recover(u[e.dofs], loads.member.get(mid, empty), self.stations[mid])
                for mid, e in self.elements.items()}

    def solve_loads(self, loads: LoadSet):
        """(Solution, {member id: MemberRecovery}) for one load set (linear)."""
        self.structure.check_inactive_loads(loads.nodal)
        sol = self.structure.solve(self.factor, self.load_vector(loads), loads.prescribed)
        return sol, self.recover(sol.u, loads)

    def solve_combination(self, combination):
        return self.solve_loads(self.combination_loads(combination))

    def checks(self, loads: LoadSet, sol, rec) -> InvariantsReport:
        return InvariantsReport(list(self.model_checks) + inv.solution_checks(
            self.structure, self.K_members, loads, sol, rec, applied_resultant(self.model, loads)))

    def result(self, combination, loads: LoadSet, sol, rec, **kw) -> CombinationResult:
        """Package one solution in the shared result structure."""
        members = {mid: MemberResult(r.end_forces, r.x, r.forces, r.displacements[:, :3].copy(),
                                     active=self.elements[mid].active) for mid, r in rec.items()}
        return CombinationResult(combination.id, combination.analysis_type, sol.u.reshape(-1, 6).copy(),
                                 sol.reactions.reshape(-1, 6).copy(), members,
                                 invariants=self.checks(loads, sol, rec),
                                 constraint_forces=sol.constraint_forces.reshape(-1, 6).copy(), **kw)


def _first_order(la: LinearAnalysis, combination) -> CombinationResult:
    """Linear elastic analysis of one combination: the factored loads are solved directly on the
    (once) factorised stiffness, which is the same as superposing the cases (checked by
    ``invariants.superposition_check``)."""
    loads = la.combination_loads(combination)
    sol, rec = la.solve_loads(loads)
    return la.result(combination, loads, sol, rec)


#: analysis type -> solver of one combination.  Other modules add their types here.
ANALYSIS_TYPES: dict[str, Callable[[LinearAnalysis, object], CombinationResult]] = {"firstOrder": _first_order}


def analyze(model, strict: bool = True) -> Results:
    """Analyse every combination of ``model``.

    Raises
    ------
    ModelRefused    unsupported feature (code ``analysis.not-implemented``), member-level mechanism,
                    or -- with ``strict`` (default) -- a failed invariant (code
                    ``solve.invariant-failed``; the flagged results are attached as ``.results``).
    ModelUnstable   mechanism / singular stiffness, with ``.diagnosis`` (MechanismDiagnosis).

    With ``strict=False`` results whose invariants failed are returned, flagged
    (``CombinationResult.invariants.passed`` is False and a warning is added); they are never
    returned unflagged.
    """
    problems = unsupported_features(model)
    if problems:
        raise ModelRefused(problems, model.source)
    la = LinearAnalysis(model)
    combos = {}
    for cid, combination in model.combinations.items():
        combos[cid] = ANALYSIS_TYPES[combination.analysis_type](la, combination)
    warnings = [f"[{p.code}] {p.path}: {p.message}" for p in model.warnings] + list(la.structure.warnings)
    res = Results(model.name, combos, defaults_used=list(model.defaults_used), warnings=warnings,
                  engine_version=__version__)
    res.compute_envelopes()
    failed = []
    for cid, c in combos.items():
        for ch in c.invariants.checks:
            if not ch.passed:
                msg = (f"invariant '{ch.name}' failed in combination {cid!r}: value {ch.value:.3e}, "
                       f"tolerance {ch.tolerance:.3e} ({ch.detail})")
                c.warnings.append(msg)
                failed.append(Problem(pointer("combinations", cid), msg, _hint(model, ch),
                                      code="solve.invariant-failed"))
    if failed:
        res.warnings.append("results are FLAGGED: at least one invariant failed (see the combinations)")
        if strict:
            e = ModelRefused(failed, model.source)
            e.results = res
            raise e
    return res


def _hint(model, check: InvariantCheck) -> str:
    if check.name == "global moment equilibrium" and any(c.type == "rigidDiaphragm" for c in model.constraints.values()):
        return ("A rigid diaphragm whose slave nodes are not in the plane of the master transfers in-plane forces "
                "between levels without the matching moment; put the slaves in the master's plane. Otherwise the "
                "system is ill-conditioned.")
    return "The system is ill-conditioned (very different stiffnesses, or close to a mechanism); check the model."
