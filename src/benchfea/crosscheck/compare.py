"""Comparison of analysis results in benchfea conventions, and the automatic tier-C decision.

Quantities compared and their scales (differences are relative to a model-wide scale of the same
physical kind, so that a small reaction is not compared against itself):

* ``displacement`` (nodal translations):   scale = largest |translation| of the two results;
* ``rotation`` (nodal rotations):          scale = max(largest |rotation|, translation scale / model size);
* ``reaction`` forces and ``endForce`` forces:   scale = largest |force| among reactions and member end forces;
* ``reaction`` moments and ``endForce`` moments: scale = max(largest |moment|, 1e-3 x force scale x model size).

``compare`` returns, per quantity, the largest relative difference with its location and both values,
and the list of the largest differences overall.  ``cross_check`` runs the foreign solvers and decides
whether the model counts as cross-checked: at least two foreign solvers must represent it (fidelity
"exact", or "approximate" with an estimated error not above ``max_approximation``) and agree with each
other -- and with the benchfea result, if one is given -- within the tolerance plus the estimated
errors of the two results being compared.

``equilibrium`` is the hand check that needs no solver: applied loads (from their definition) plus
reactions must vanish.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .fidelity import EXACT, UNSUPPORTED, Fidelity
from .foreign import ForeignResult
from .prepare import prepare

__all__ = ["Tolerances", "Difference", "Comparison", "CrossCheckReport", "compare", "relative_difference",
           "from_benchfea", "cross_check", "equilibrium", "applied_resultant"]

DOF_NAMES = ("ux", "uy", "uz", "rx", "ry", "rz")
FORCE_NAMES = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")
QUANTITIES = ("displacement", "rotation", "reactionForce", "reactionMoment", "endForce", "endMoment")
MOMENT_FLOOR = 1e-3


@dataclass
class Tolerances:
    """Relative tolerances per quantity (see module docstring for the scales)."""
    displacement: float = 1e-6
    rotation: float = 1e-6
    reactionForce: float = 1e-6
    reactionMoment: float = 1e-6
    endForce: float = 1e-6
    endMoment: float = 1e-6

    def widened(self, extra: float) -> "Tolerances":
        return Tolerances(*(getattr(self, q) + extra for q in QUANTITIES))


@dataclass
class Difference:
    quantity: str
    where: str                    # "node N3 uz" | "member B1 end j Mz"
    a: float
    b: float
    relative: float               # |a - b| / scale

    def to_dict(self) -> dict:
        return {"quantity": self.quantity, "where": self.where, "a": self.a, "b": self.b, "relative": self.relative}


@dataclass
class Comparison:
    a: str
    b: str
    combination: str
    worst: dict                   # quantity -> Difference (largest of that quantity)
    largest: list                 # the largest differences overall, descending
    tolerances: Tolerances
    scales: dict
    passed: bool
    failures: list = field(default_factory=list)      # quantities above tolerance

    @property
    def max_relative(self) -> float:
        return max((d.relative for d in self.worst.values()), default=0.0)

    def to_dict(self) -> dict:
        return {"a": self.a, "b": self.b, "combination": self.combination, "passed": self.passed,
                "maxRelative": self.max_relative, "failures": self.failures, "scales": self.scales,
                "tolerances": {q: getattr(self.tolerances, q) for q in QUANTITIES},
                "worst": {q: d.to_dict() for q, d in self.worst.items()},
                "largest": [d.to_dict() for d in self.largest]}

    def summary(self) -> str:
        lines = [f"{self.a} vs {self.b} [{self.combination}]: {'agree' if self.passed else 'DIFFER'} "
                 f"(max relative difference {self.max_relative:.3e})"]
        for q, d in self.worst.items():
            mark = "" if q not in self.failures else "  <-- above tolerance"
            lines.append(f"  {q:15s} {d.relative:.3e} at {d.where}: {d.a:.9g} vs {d.b:.9g} "
                         f"(tol {getattr(self.tolerances, q):.1e}){mark}")
        return "\n".join(lines)


def _size(model) -> float:
    pts = np.array([n.xyz for n in model.nodes.values()], float)
    return float(np.linalg.norm(np.ptp(pts, axis=0))) or 1.0


def _differences(a: ForeignResult, b: ForeignResult, model):
    ids = list(model.nodes)
    size = _size(model)
    ua, ub = np.asarray(a.displacements, float), np.asarray(b.displacements, float)
    ra, rb = np.asarray(a.reactions, float), np.asarray(b.reactions, float)
    common = [m for m in model.members if m in a.end_forces and m in b.end_forces]
    ea = np.array([a.end_forces[m] for m in common], float).reshape(len(common), 12)
    eb = np.array([b.end_forces[m] for m in common], float).reshape(len(common), 12)
    fcol, mcol = [0, 1, 2, 6, 7, 8], [3, 4, 5, 9, 10, 11]
    amax = lambda *arrs: max((float(np.abs(x).max()) for x in arrs if x.size), default=0.0)
    s_t = amax(ua[:, :3], ub[:, :3])
    s_r = max(amax(ua[:, 3:], ub[:, 3:]), s_t / size)
    s_f = amax(ra[:, :3], rb[:, :3], ea[:, fcol], eb[:, fcol])
    s_m = max(amax(ra[:, 3:], rb[:, 3:], ea[:, mcol], eb[:, mcol]), MOMENT_FLOOR * s_f * size)
    scales = {"displacement": s_t, "rotation": s_r, "reactionForce": s_f, "reactionMoment": s_m,
              "endForce": s_f, "endMoment": s_m}
    out = []

    def nodal(quantity, x, y, cols, names):
        s = scales[quantity]
        if s == 0.0:
            return
        d = np.abs(x[:, cols] - y[:, cols]) / s
        for i, j in zip(*np.nonzero(d > 0.0)):
            out.append(Difference(quantity, f"node {ids[i]} {names[cols[j]]}", float(x[i, cols[j]]),
                                  float(y[i, cols[j]]), float(d[i, j])))

    nodal("displacement", ua, ub, [0, 1, 2], DOF_NAMES)
    nodal("rotation", ua, ub, [3, 4, 5], DOF_NAMES)
    nodal("reactionForce", ra, rb, [0, 1, 2], FORCE_NAMES)
    nodal("reactionMoment", ra, rb, [3, 4, 5], FORCE_NAMES)
    for quantity, cols in (("endForce", fcol), ("endMoment", mcol)):
        s = scales[quantity]
        if s == 0.0 or not common:
            continue
        d = np.abs(ea[:, cols] - eb[:, cols]) / s
        for i, j in zip(*np.nonzero(d > 0.0)):
            c = cols[j]
            out.append(Difference(quantity, f"member {common[i]} end {'ij'[c // 6]} {FORCE_NAMES[c % 6]}",
                                  float(ea[i, c]), float(eb[i, c]), float(d[i, j])))
    return out, scales


def compare(a: ForeignResult, b: ForeignResult, model, tolerances: Tolerances | None = None,
            n_largest: int = 10) -> Comparison:
    """Compare two results of the same combination (both must be ``ok``)."""
    if not (a.ok and b.ok):
        raise ValueError("cannot compare a failed result: " + (a.error or b.error))
    tol = tolerances or Tolerances()
    diffs, scales = _differences(a, b, model)
    worst = {}
    for q in QUANTITIES:
        of_q = [d for d in diffs if d.quantity == q]
        worst[q] = max(of_q, key=lambda d: d.relative) if of_q else Difference(q, "-", 0.0, 0.0, 0.0)
    failures = [q for q in QUANTITIES if not worst[q].relative <= getattr(tol, q)]
    largest = sorted(diffs, key=lambda d: -d.relative)[:n_largest]
    return Comparison(a.solver, b.solver, a.combination, worst, largest, tol, scales, not failures, failures)


def relative_difference(a: ForeignResult, b: ForeignResult, model) -> float:
    """Largest relative difference over all quantities (scales of the module docstring)."""
    diffs, _ = _differences(a, b, model)
    return max((d.relative for d in diffs), default=0.0)


def from_benchfea(model, results, combination: str) -> ForeignResult:
    """A ``benchfea.results.Results`` (or one ``CombinationResult``) as a comparable record."""
    c = results.combinations[combination] if hasattr(results, "combinations") else results
    fid = Fidelity("benchfea")
    fid.add("engine", EXACT, "the analysis being checked")
    return ForeignResult("benchfea", combination, True, fid, np.asarray(c.displacements, float),
                         np.asarray(c.reactions, float),
                         {mid: np.asarray(m.end_forces, float) for mid, m in c.members.items()},
                         info={"converged": c.converged, "iterations": c.iterations})


# ====================================================================================== hand equilibrium
def applied_resultant(model, combination: str, extra_nodal_loads: dict | None = None):
    """Resultant (Fx, Fy, Fz, Mx, My, Mz about the global origin) of the factored loads of a
    combination, straight from the load definitions (no element, no equivalent nodal load), and the
    sum of the absolute contributions (scale).  For a second-order combination the loads are NOT
    multiplied by alpha here (results are reported divided by alpha)."""
    p = prepare(model, combination, extra_nodal_loads=extra_nodal_loads)
    R, S = np.zeros(6), np.zeros(6)

    def add(point, force, moment):
        mo = np.cross(point, force)
        R[:3] += force
        R[3:] += moment + mo
        S[:3] += np.abs(force)
        S[3:] += np.abs(moment) + np.abs(mo)

    for k in range(len(p.node_ids)):
        if p.load[k].any():
            add(p.xyz[k], p.load[k, :3], p.load[k, 3:])
    for m in p.members.values():
        Rt = m.R.T
        for kind, vec, x in m.points:
            g = Rt @ vec
            add(m.start + x * Rt[:, 0], g if kind == "force" else np.zeros(3), g if kind == "moment" else np.zeros(3))
        for w1, w2, x1, x2 in m.dists:
            h = x2 - x1
            # a linear load = two triangles with resultants w1 h/2 at x1 + h/3 and w2 h/2 at x1 + 2h/3
            add(m.start + (x1 + h / 3.0) * Rt[:, 0], Rt @ w1 * (h / 2.0), np.zeros(3))
            add(m.start + (x1 + 2.0 * h / 3.0) * Rt[:, 0], Rt @ w2 * (h / 2.0), np.zeros(3))
    return R / p.alpha, S / p.alpha


def equilibrium(model, result: ForeignResult, extra_nodal_loads: dict | None = None) -> dict:
    """Global equilibrium of applied loads + reactions: relative residuals ``force`` and ``moment``
    (about the origin).  In a second-order analysis the moment balance on the undeformed geometry does
    not hold (P-Delta couples), so ``moment`` is reported as ``None`` there."""
    R, S = applied_resultant(model, result.combination, extra_nodal_loads)
    xyz = np.array([n.xyz for n in model.nodes.values()], float)
    rf, rm = result.reactions[:, :3], result.reactions[:, 3:]
    total = R + np.r_[rf.sum(axis=0), rm.sum(axis=0) + np.cross(xyz, rf).sum(axis=0)]
    fs = max(float(S[:3].max()), float(np.abs(rf).max()), 1e-300)
    ms = max(float(S[3:].max()), float(np.abs(rm).max()), float(np.abs(np.cross(xyz, rf)).max()), 1e-300)
    second = model.combinations[result.combination].analysis_type == "secondOrder"
    return {"force": float(np.abs(total[:3]).max()) / fs,
            "moment": None if second else float(np.abs(total[3:]).max()) / ms}


# ====================================================================================== tier C
@dataclass
class CrossCheckReport:
    combination: str
    results: dict                 # solver -> ForeignResult (benchfea included when given)
    comparisons: list             # Comparison of every pair that could be compared
    eligible: list                # foreign solvers that represent the model (exactly / converged approximation)
    cross_checked: bool
    reasons: list                 # why (not)
    equilibrium: dict             # solver -> {"force", "moment"}

    def to_dict(self) -> dict:
        return {"combination": self.combination, "crossChecked": self.cross_checked, "reasons": self.reasons,
                "eligible": self.eligible, "equilibrium": self.equilibrium,
                "results": {k: {"ok": r.ok, "error": r.error, "fidelity": r.fidelity.to_dict(), "info": r.info}
                            for k, r in self.results.items()},
                "comparisons": [c.to_dict() for c in self.comparisons]}

    def summary(self) -> str:
        lines = [f"combination {self.combination}: {'CROSS-CHECKED' if self.cross_checked else 'not cross-checked'}"]
        lines += ["  - " + r for r in self.reasons]
        for c in self.comparisons:
            lines.append(c.summary())
        return "\n".join(lines)


def cross_check(model, combination: str, benchfea_result=None, solvers=("opensees", "pynite"),
                tolerances: Tolerances | None = None, max_approximation: float = 5e-3, options: dict | None = None,
                timeout: float | None = None, foreign: dict | None = None) -> CrossCheckReport:
    """Run the foreign solvers on one combination and decide whether it is cross-checked.

    ``benchfea_result``: ``Results`` / ``CombinationResult`` of the analysis being checked (optional).
    ``options``: ``{solver: {translate options}}``.  ``foreign``: already computed
    ``{solver: ForeignResult}`` (then the solver is not run again).
    """
    from .runner import DEFAULT_TIMEOUT, run

    tol = tolerances or Tolerances()
    results: dict[str, ForeignResult] = {}
    reasons, eligible = [], []
    for s in solvers:
        if foreign and s in foreign:
            results[s] = foreign[s]
        else:
            results[s] = run(s, model, combination, timeout=timeout or DEFAULT_TIMEOUT, **(options or {}).get(s, {}))
        r = results[s]
        if not r.ok:
            reasons.append(f"{s}: failed: {r.error}")
            continue
        bad = r.fidelity.unsupported
        if bad:
            reasons.append(f"{s}: cannot represent the model: " +
                           "; ".join(f"{i.feature} ({i.note})" for i in bad))
            continue
        err = r.fidelity.estimated_error
        if err is None:
            reasons.append(f"{s}: approximation without an error estimate: " +
                           "; ".join(i.feature for i in r.fidelity.approximate if i.error is None))
            continue
        if err > max_approximation:
            reasons.append(f"{s}: estimated approximation error {err:.2e} exceeds {max_approximation:.1e}")
            continue
        eligible.append(s)
    extra_loads = ((options or {}).get(solvers[0], {}) or {}).get("extra_nodal_loads") if solvers else None
    eq = {s: equilibrium(model, r, extra_loads) for s, r in results.items() if r.ok}
    comparisons = []
    names = list(eligible)
    if benchfea_result is not None:
        results["benchfea"] = from_benchfea(model, benchfea_result, combination)
        eq["benchfea"] = equilibrium(model, results["benchfea"], extra_loads)
        names = ["benchfea"] + names
    agree = True
    for i, x in enumerate(names):
        for y in names[i + 1:]:
            widen = sum(results[z].fidelity.estimated_error or 0.0 for z in (x, y))
            c = compare(results[x], results[y], model, tol.widened(widen))
            comparisons.append(c)
            if not c.passed:
                agree = False
                d = max(c.worst.values(), key=lambda v: v.relative)
                reasons.append(f"{x} and {y} differ: {d.quantity} {d.relative:.2e} at {d.where}")
    if len(eligible) < 2:
        reasons.append(f"only {len(eligible)} foreign solver(s) represent this model; two are required")
    ok = len(eligible) >= 2 and agree
    if ok:
        reasons.append(f"{' and '.join(eligible)} represent the model and agree" +
                       (" with each other and with benchfea" if benchfea_result is not None else ""))
    return CrossCheckReport(combination, results, comparisons, eligible, ok, reasons, eq)
