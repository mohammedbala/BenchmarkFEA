"""``analyze(model) -> Results``: the driver.

Structure of the code (and where a nonlinear driver plugs in)
------------------------------------------------------------
``LinearAnalysis``  builds the elements (``build_element``), the ``Structure`` (constraints, supports,
                    assembly), factorises the elastic stiffness once and solves any ``LoadSet``:
                    ``solve_loads(loads) -> (Solution, {member: MemberRecovery})``.
``ANALYSIS_TYPES``  maps ``Combination.analysis_type`` to a function
                    ``(LinearAnalysis, Combination) -> CombinationResult``: ``"firstOrder"`` (here)
                    and ``"secondOrder"`` (``nonlinear.py``); a type without an entry is refused with
                    ``analysis.not-implemented``.  First-order combinations with tension/compression-
                    only members, stiffness reduction or notional loads are also solved by
                    ``nonlinear.NonlinearSolver`` (active set, tau_b iteration); linear buckling by
                    ``buckling.py``.

This module never reads the design-only data of the model.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from .. import __version__
from ..problems import ModelUnstable, Problem, pointer
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
    """Model features that no registered solver analyses.  Each is refused rather than ignored.
    (Second-order analysis, tension/compression-only members, stiffness reduction, notional loads
    and linear buckling are implemented by ``nonlinear.py`` / ``buckling.py``.)"""
    out = []
    for cid, c in model.combinations.items():
        if c.analysis_type not in ANALYSIS_TYPES:
            out.append(Problem(pointer("combinations", cid, "analysisType"),
                               f"analysis type {c.analysis_type!r} (combination {cid!r}) is not implemented in this build",
                               "Registered analysis types: " + ", ".join(sorted(ANALYSIS_TYPES)) + ".",
                               code=NOT_IMPLEMENTED))
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


def _nonlinear(la: LinearAnalysis):
    """The nonlinear solver of this model (one per LinearAnalysis)."""
    if getattr(la, "_nonlinear", None) is None:
        from .nonlinear import NonlinearSolver
        la._nonlinear = NonlinearSolver(la)
    return la._nonlinear


def _first_order(la: LinearAnalysis, combination) -> CombinationResult:
    """Linear elastic analysis of one combination: the factored loads are solved directly on the
    (once) factorised stiffness, which is the same as superposing the cases (checked by
    ``invariants.superposition_check``).  With tension/compression-only members, stiffness reduction
    or notional loads the combination goes through ``nonlinear.NonlinearSolver`` instead."""
    from .nonlinear import needs_nonlinear
    if needs_nonlinear(la.model, combination):
        return _nonlinear(la).solve(combination)
    loads = la.combination_loads(combination)
    sol, rec = la.solve_loads(loads)
    return la.result(combination, loads, sol, rec)


def _second_order(la: LinearAnalysis, combination) -> CombinationResult:
    """Second-order analysis of one combination (``nonlinear.NonlinearSolver``)."""
    return _nonlinear(la).solve(combination)


#: analysis type -> solver of one combination.  Other modules add their types here.
ANALYSIS_TYPES: dict[str, Callable[[LinearAnalysis, object], CombinationResult]] = {
    "firstOrder": _first_order, "secondOrder": _second_order}


def analyze(model, strict: bool = True) -> Results:
    """Analyse every combination of ``model`` (each one on its own), then the linear buckling
    analyses listed in ``analysis.buckling``.

    Raises
    ------
    ModelRefused    unsupported feature (code ``analysis.not-implemented``), member-level mechanism,
                    a combination that cannot be analysed (``stability.beyond-critical-load``,
                    ``solve.active-set-cycling``, ``analysis.release-second-order``,
                    ``analysis.tau-b-zero``), or -- with ``strict`` (default) -- a failed invariant
                    (``solve.invariant-failed``) or an iteration that did not converge
                    (``solve.not-converged``).  The results of the other combinations (and the flagged
                    ones) are attached as ``.results``.
    ModelUnstable   mechanism / singular stiffness, with ``.diagnosis`` (MechanismDiagnosis); also when
                    slack tension/compression-only members leave a mechanism in a combination
                    (``stability.slack-mechanism``; ``.results`` then holds the other combinations).

    With ``strict=False`` flagged results (failed invariant, not converged) are returned with a
    warning, and a refused combination is left out of the results with a warning; they are never
    returned unflagged.
    """
    from .nonlinear import CombinationRefused
    problems = unsupported_features(model)
    if problems:
        raise ModelRefused(problems, model.source)
    la = LinearAnalysis(model)
    combos, refused = {}, {}
    for cid, combination in model.combinations.items():
        try:
            combos[cid] = ANALYSIS_TYPES[combination.analysis_type](la, combination)
        except CombinationRefused as r:
            refused[f"combination {cid!r}"] = r
    warnings = [f"[{p.code}] {p.path}: {p.message}" for p in model.warnings] + list(la.structure.warnings)
    res = Results(model.name, combos, defaults_used=list(model.defaults_used), warnings=warnings,
                  engine_version=__version__)
    res.compute_envelopes()
    if model.analysis.buckling_combinations:
        from .buckling import buckling_analysis
        for cid in model.analysis.buckling_combinations:
            try:
                res.buckling[cid] = buckling_analysis(_nonlinear(la), model.combinations[cid])
            except CombinationRefused as r:
                refused[f"buckling analysis of combination {cid!r}"] = r
    failed = []
    for cid, c in combos.items():
        for ch in c.invariants.checks:
            if not ch.passed:
                msg = (f"invariant '{ch.name}' failed in combination {cid!r}: value {ch.value:.3e}, "
                       f"tolerance {ch.tolerance:.3e} ({ch.detail})")
                c.warnings.append(msg)
                failed.append(Problem(pointer("combinations", cid), msg, _hint(model, ch),
                                      code="solve.invariant-failed"))
        if not c.converged:
            failed.append(Problem(pointer("combinations", cid),
                                  f"combination {cid!r} did not converge: " + "; ".join(c.warnings[:1]),
                                  "Increase analysis.maxIterations or check whether the structure is close to "
                                  "its elastic critical load.", code="solve.not-converged"))
    if failed:
        res.warnings.append("results are FLAGGED: at least one invariant failed or an iteration did not converge "
                            "(see the combinations)")
    if refused:
        for what, r in refused.items():
            res.warnings.append(f"NO RESULTS for {what}: " + "; ".join(p.message for p in r.problems))
        diagnosis = next((r.diagnosis for r in refused.values() if r.diagnosis is not None), None)
        allp = [p for r in refused.values() for p in r.problems] + failed
        e = (ModelUnstable if diagnosis is not None else ModelRefused)(allp, model.source)
        if diagnosis is not None:
            e.diagnosis = diagnosis
        e.results = res
        e.refused = {what: r for what, r in refused.items()}
        if strict:
            raise e
        return res
    if failed and strict:
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
