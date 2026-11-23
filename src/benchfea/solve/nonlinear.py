"""Nonlinear analysis of one combination: tension/compression-only members (active set), second-order
analysis (geometric stiffness, iterated axial forces), the stiffness reduction of AISC 360-22 C2.3 and
the notional loads of C2.2b.  Theory and rules: docs/THEORY.md sections 10-14.

Every combination is solved on its own from its factored loads (no superposition).

Loop structure
--------------
active set (outer)      which tension/compression-only members are slack.  Start: none.  A member is
                        slack when the axial force it carries (or would carry) has the wrong sign; the
                        test quantity is its mean axial force E A (elongation) / L.  A slack member has
                        zero stiffness; its member loads still go to its end nodes.
state (inner)           first order: one linear solve per value of tau_b.
                        second order: K_T(N) u = F(N) with the axial forces N (and tau_b) of the
                        previous iterate, until the displacement and axial-force changes are below
                        ``analysis.convergenceTolerance``; members whose axial force varies along
                        their length are subdivided and the subdivision is doubled until the results
                        change by less than the same tolerance.

This module never reads the design-only data of the model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..problems import ModelUnstable, Problem, pointer
from ..results import CombinationResult, InvariantCheck, InvariantsReport, MemberResult
from . import invariants as inv
from .beam import SpanLoads
from .beamcolumn import BeamColumn, axial_state
from .element import TRUSS_TYPES, StiffnessFactors, build_element, section_flex
from .loads import LoadSet, applied_resultant, combine
from .refine import piece_stations, refine, split_loads
from .system import Structure

__all__ = ["NonlinearSolver", "CombinationRefused", "notional_forces", "tau_b", "needs_nonlinear",
           "ONE_SIDED", "SUBDIVISION_START", "SUBDIVISION_MAX", "KH_TARGET", "KH_LIMIT", "SLACK_TOL"]

ONE_SIDED = ("tensionOnly", "compressionOnly")
SUBDIVISION_START = 2      # pieces of a member whose axial force varies along its length (then doubled)
SUBDIVISION_MAX = 32
KH_TARGET = 1.5            # a-priori piece length: h sqrt(|N| / EI) <= KH_TARGET (with 2 x the first-order N)
KH_LIMIT = 3.0             # verified after convergence; otherwise the subdivision is doubled
SLACK_TOL = 1.0e-9         # a member keeps its state while |strain| <= SLACK_TOL x largest member strain
BEYOND = "stability.beyond-critical-load"
SLACK_MECHANISM = "stability.slack-mechanism"
NOT_CONVERGED = "solve.not-converged"


class CombinationRefused(Exception):
    """One combination cannot be analysed.  ``problems`` say why; ``diagnosis`` is the mechanism
    diagnosis when slack members leave a mechanism; ``info`` holds numbers (e.g. the buckling factor)."""

    def __init__(self, problems, diagnosis=None, info=None):
        super().__init__("; ".join(p.message for p in problems))
        self.problems = list(problems)
        self.diagnosis = diagnosis
        self.info = info or {}


class _Critical(Exception):
    def __init__(self, what, mesh=None, rloads=None, state=None, iteration=0, cause=None):
        super().__init__(what)
        self.what, self.mesh, self.rloads, self.state, self.iteration, self.cause = what, mesh, rloads, state, iteration, cause


def needs_nonlinear(model, combination) -> bool:
    """True when a FIRST-ORDER combination cannot be solved by the plain linear path."""
    return (combination.stiffness_reduction or combination.notional is not None
            or any(m.type in ONE_SIDED for m in model.members.values()))


def tau_b(ratio: float) -> float:
    """AISC 360-22 C2.3(b): Eq. C2-2a (ratio = alpha Pr / Pns <= 0.5: 1.0) and Eq. C2-2b
    (ratio > 0.5: 4 ratio (1 - ratio))."""
    return 1.0 if ratio <= 0.5 else 4.0 * ratio * (1.0 - ratio)


def notional_forces(model, loads: LoadSet, notional) -> np.ndarray:
    """(n_nodes, 3) notional loads, AISC 360-22 C2.2b Eq. C2-1 as the model layer defines it:
    at every node ``coefficient x Y`` along ``direction``, Y = the downward (minus vertical) component
    of the loads applied at the node, member loads lumped to the two end nodes of their member by the
    lever rule (statically).  Y is signed: a net upward load gives a notional load against
    ``direction``.  Moments have no vertical component and imposed displacements are not loads."""
    up = model.vertical
    Y = np.zeros(len(model.nodes))
    for node in model.nodes.values():
        Y[node.index] -= float(loads.nodal[6 * node.index:6 * node.index + 3] @ up)
    for mid, span in loads.member.items():
        m = model.members[mid]
        ez = np.asarray(m.R, float) @ up                     # vertical unit vector in local components
        i, j, L = model.nodes[m.i].index, model.nodes[m.j].index, float(m.length)

        def put(x, f):
            y = -float(np.asarray(f, float)[:3] @ ez)
            Y[i] += y * (1.0 - x / L)
            Y[j] += y * (x / L)
        for x, f in span.points:
            put(min(max(float(x), 0.0), L), f)
        for x1, x2, w1, w2 in span.dists:
            h = x2 - x1                                      # a linear load = two triangles
            put(x1 + h / 3.0, np.asarray(w1, float) * (h / 2.0))
            put(x1 + 2.0 * h / 3.0, np.asarray(w2, float) * (h / 2.0))
    return notional.coefficient * Y[:, None] * np.asarray(notional.direction, float)[None, :]


@dataclass
class _State:
    """A solved state of one combination for one active set."""
    second_order: bool
    u: np.ndarray                         # (n_dof,) original nodes
    reactions: np.ndarray
    constraint_forces: np.ndarray
    members: dict                         # member id -> MemberRecovery-like (end_forces, x, forces, displacements)
    tau: dict
    checks: list
    converged: bool = True
    iterations: int = 1
    history: list = field(default_factory=list)
    subdivisions: dict = field(default_factory=dict)
    info: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    linear: object = None                 # LinearAnalysis (first-order states)
    sol: object = None
    rec: object = None
    n_extreme: dict = field(default_factory=dict)      # member id -> (n_min, n_max)


class NonlinearSolver:
    def __init__(self, la):
        self.la = la
        self.model = la.model
        self._linear_cache: dict = {}
        m = self.model
        red = m.analysis.reduction_members
        self.reduced = set(m.members) if red is None else set(red)
        self.size = la.structure.size if la.structure.size > 0.0 else 1.0

    # ------------------------------------------------------------ stiffness reduction
    def _factors(self, mid, combination, tau) -> StiffnessFactors | None:
        """AISC 360-22 C2.3: (a) 0.80 on all stiffnesses of the listed members, (b) tau_b additionally on
        their flexural stiffness (the User Note of C2.3(b): 0.8 tau_b EI, 0.8 x the other stiffnesses)."""
        if not combination.stiffness_reduction or mid not in self.reduced:
            return None
        f = float(self.model.analysis.reduction_factor)
        return StiffnessFactors(axial=f, flexural=f * tau.get(mid, 1.0), torsional=f, shear=f)

    def _tau_start(self, combination) -> dict:
        if not combination.stiffness_reduction:
            return {}
        t = self.model.analysis.tau_b
        return {mid: (float(t) if isinstance(t, (int, float)) else 1.0) for mid in self.model.members if mid in self.reduced}

    def _tau_update(self, combination, n_extreme: dict, tau: dict) -> dict:
        """tau_b of every reduced member from its largest compression (alpha Pr = the compression in
        this analysis, which is run at alpha x the loads; Pns = Fy A)."""
        if not combination.stiffness_reduction or self.model.analysis.tau_b != "computed":
            return dict(tau)
        out = {}
        for mid in tau:
            m = self.model.members[mid]
            # Pns = Fy A with the yield stress of the member's MATERIAL (analysis data; a design-block
            # override of Fy never reaches the solver -- docs/ARCHITECTURE.md "Analysis vs design data")
            pns = getattr(self.model.materials[m.material], "Fy") * self.model.sections[m.section].A
            ratio = max(0.0, -n_extreme.get(mid, (0.0, 0.0))[0]) / pns
            if ratio >= 1.0:
                raise CombinationRefused([Problem(
                    pointer("combinations", combination.id, "stiffnessReduction"),
                    f"member {mid!r}: the axial compression is {ratio:.3f} times its cross-section compressive "
                    f"strength Pns = Fy A, so tau_b (AISC 360-22 Eq. C2-2b) is zero or negative and the member "
                    "has no flexural stiffness left",
                    "The member is overloaded in compression; increase the section or reduce the load.",
                    code="analysis.tau-b-zero")])
            out[mid] = tau_b(ratio)
        return out

    # ----------------------------------------------------------------------- loads
    def loads(self, combination, second_order: bool):
        """(LoadSet at the analysis level, alpha, notional (n_nodes, 3) or None)."""
        model = self.model
        loads = self.la.combination_loads(combination)
        alpha = float(model.analysis.alpha) if second_order else 1.0
        if alpha != 1.0:
            loads = combine(model, [(alpha, loads)])
        notional = None
        if combination.notional is not None:
            notional = notional_forces(model, loads, combination.notional)
            nodal = loads.nodal.copy()
            nodal.reshape(-1, 6)[:, :3] += notional
            loads = LoadSet(nodal, loads.prescribed, loads.member)
        return loads, alpha, notional

    # ------------------------------------------------------------------ first order
    def _linear(self, combination, slack, tau):
        from .analysis import LinearAnalysis
        if not slack and not combination.stiffness_reduction:
            return self.la
        key = (frozenset(slack), tuple(sorted(tau.items())) if combination.stiffness_reduction else None)
        if key not in self._linear_cache:
            model = self.model
            self._linear_cache[key] = LinearAnalysis(model, {
                mid: build_element(model, m, self._factors(mid, combination, tau), active=mid not in slack)
                for mid, m in model.members.items()})
        return self._linear_cache[key]

    def first_order(self, combination, loads, slack) -> _State:
        """Linear solve for a given active set; iterates tau_b when it is computed."""
        tol, itmax = float(self.model.analysis.convergence_tolerance), int(self.model.analysis.max_iterations)
        tau = self._tau_start(combination)
        history, converged = [], True
        for it in range(1, itmax + 1):
            la = self._linear(combination, slack, tau)
            sol, rec = la.solve_loads(loads)
            ext = {mid: (float(r.forces[:, 0].min()), float(r.forces[:, 0].max())) for mid, r in rec.items()}
            new = self._tau_update(combination, ext, tau)
            change = max([abs(new[k] - tau[k]) for k in tau], default=0.0)
            history.append({"iteration": it, "tauBChange": change})
            if change <= tol:
                break
            tau = new
        else:
            converged = False
        return _State(False, sol.u, sol.reactions, sol.constraint_forces, rec, tau, la.checks(loads, sol, rec).checks,
                      converged=converged, iterations=len(history), history=history, linear=la, sol=sol, rec=rec,
                      n_extreme=ext)

    # ----------------------------------------------------------------- second order
    def second_order(self, combination, loads, slack) -> _State:
        tol = float(self.model.analysis.convergence_tolerance)
        first = self.first_order(combination, loads, slack)         # ModelUnstable if the elastic structure is a mechanism
        varying = {mid for mid, span in loads.member.items() if mid not in slack
                   and any(w1[0] != 0.0 or w2[0] != 0.0 for _, _, w1, w2 in span.dists)}
        level, prev, check = SUBDIVISION_START, None, None
        while True:
            run = self._iterate(combination, loads, slack, first, varying, level)
            if not varying:
                break
            if prev is not None:
                change = _difference(prev, run, self.size)
                check = {"pieces": [level // 2, level], "change": change, "tolerance": tol}
                if change <= tol or not run.converged:
                    break
                if level >= SUBDIVISION_MAX:
                    run.converged = False
                    run.warnings.append(
                        f"subdivision of the members with a varying axial force did not converge: results still "
                        f"change by {change:.2e} between {level // 2} and {level} pieces per member "
                        f"(tolerance {tol:.1e})")
                    break
            prev, level = run, level * 2
        if check:
            run.info["subdivisionCheck"] = check
        return run

    def _divisions(self, combination, loads, slack, first, varying, level, boost):
        model, out = self.model, {}
        for mid, m in model.members.items():
            if mid in slack:
                continue
            flex = section_flex(model, m, False, self._factors(mid, combination, first.tau))
            nmax = max(abs(v) for v in first.n_extreme[mid])
            kh = float(m.length) * math.sqrt(2.0 * boost * nmax * max(flex.bend_y, flex.bend_z))
            out[mid] = max(int(math.ceil(kh / KH_TARGET)), 1) * (level if mid in varying else 1)
        return out

    def _iterate(self, combination, loads, slack, first, varying, level, boost=1.0) -> _State:
        model = self.model
        tol, itmax = float(model.analysis.convergence_tolerance), int(model.analysis.max_iterations)
        mesh = refine(model, self._divisions(combination, loads, slack, first, varying, level, boost), unsplit=slack)
        rmodel = mesh.model
        rloads = split_loads(mesh, loads)
        owner = {pid: mid for mid, ids in mesh.pieces.items() for pid in ids}
        nb = mesh.n_dof_base
        empty = SpanLoads()
        # axial force just inside the start of every piece, from the first-order solution
        Ns = {}
        for mid, ids in mesh.pieces.items():
            if mid in slack:
                continue
            span = loads.member.get(mid, empty)
            n0 = -first.rec[mid].end_forces[0]
            for k, pid in enumerate(ids):
                Ns[pid] = _axial_after(span, n0, float(mesh.bounds[mid][k]))
        tau = dict(first.tau)
        u_prev = first.u
        scale = np.where(np.arange(nb) % 6 < 3, 1.0, self.size)
        history, converged = [], False
        elems = st = Km = fac = sol = None
        for it in range(1, itmax + 1):
            elems = {}
            for pid, pm in rmodel.members.items():
                mid = owner[pid]
                fct = self._factors(mid, combination, tau)
                if mid in slack:
                    elems[pid] = build_element(rmodel, pm, fct, active=False)
                else:
                    elems[pid] = BeamColumn.build(rmodel, pm, rloads.member.get(pid), Ns[pid], fct)
            state = (mesh, rloads, dict(Ns), dict(tau))
            bad = sorted({owner[pid] for pid, e in elems.items() if isinstance(e, BeamColumn) and not e.stable})
            if bad:
                raise _Critical("member(s) " + ", ".join(bad) + " exceed their own elastic critical load between "
                                "their end connections", *state[:2], state, it)
            st = Structure(rmodel, elems)
            Km = st.assemble()
            K = Km.copy()
            K[np.diag_indices(st.n)] += st.spring
            try:
                fac = st.factorize(K)
            except ModelUnstable as e:
                raise _Critical("the tangent stiffness K + Kg of the structure is not positive definite",
                                mesh, rloads, state, it, e) from None
            F = rloads.nodal.copy()
            for pid, e in elems.items():
                F[e.dofs] += e.equivalent_nodal_loads(rloads.member.get(pid, empty))
            st.check_inactive_loads(rloads.nodal)
            sol = st.solve(fac, F, rloads.prescribed)
            new_ns, ext = {}, {}
            for pid, e in elems.items():
                if isinstance(e, BeamColumn):
                    r = e.recover(sol.u[e.dofs])
                    new_ns[pid] = r.n_start
                    lo, hi = ext.get(owner[pid], (math.inf, -math.inf))
                    ext[owner[pid]] = (min(lo, r.n_min), max(hi, r.n_max))
            new_tau = self._tau_update(combination, ext, tau)
            du = _rel(np.linalg.norm((sol.u[:nb] - u_prev) * scale), np.linalg.norm(sol.u[:nb] * scale))
            nref = max((abs(v) for v in new_ns.values()), default=0.0)
            dn = _rel(max((abs(new_ns[p] - Ns[p]) for p in Ns), default=0.0), nref)
            dt = max([abs(new_tau[k] - tau[k]) for k in tau], default=0.0)
            history.append({"iteration": it, "displacementChange": du, "axialForceChange": dn, "tauBChange": dt})
            if not all(math.isfinite(x) for x in (du, dn, dt)):
                break
            if du <= tol and dn <= tol and dt <= tol:
                converged = True
                break
            Ns, tau, u_prev = new_ns, new_tau, sol.u[:nb].copy()
        # ---- a-posteriori check of the piece lengths
        kh = 0.0
        for pid, e in elems.items():
            if isinstance(e, BeamColumn):
                nmax = max(abs(e.axial[0]), abs(e.axial[0] - e.axial[1] * e.length - 0.5 * e.axial[2] * e.length ** 2))
                kh = max(kh, e.length * math.sqrt(nmax * max(e.flex.bend_y, e.flex.bend_z)))
        if kh > KH_LIMIT and boost < 60.0:
            return self._iterate(combination, loads, slack, first, varying, level, boost * 4.0)
        # ---- final recovery at the stations, mapped back to the original members
        recs, members = {}, {}
        for mid, ids in mesh.pieces.items():
            stations = self.la.stations[mid]
            if mid in slack:
                e = elems[mid]
                recs[mid] = members[mid] = e.recover(sol.u[e.dofs], rloads.member.get(mid, empty), stations)
                continue
            where = piece_stations(mesh.bounds[mid], stations)
            forces, disp = np.zeros((len(stations), 6)), np.zeros((len(stations), 6))
            for k, pid in enumerate(ids):
                e = elems[pid]
                loc = [(x, side) for _, x, side in where[k]]
                r = recs[pid] = e.recover(sol.u[e.dofs], loc, checks=True)
                for row, (n, _, _) in enumerate(where[k]):
                    forces[n], disp[n] = r.forces[row], r.displacements[row]
            ends = np.concatenate([recs[ids[0]].end_forces[:6], recs[ids[-1]].end_forces[6:]])
            members[mid] = _Recovery(ends, np.array([x for x, _ in stations]), forces, disp)
        checks = second_order_checks(st, Km, fac, rloads, sol, elems, recs)
        out = _State(True, sol.u[:nb].copy(), sol.reactions[:nb].copy(), sol.constraint_forces[:nb].copy(), members,
                     tau, checks, converged=converged, iterations=len(history), history=history,
                     subdivisions={mid: len(ids) for mid, ids in mesh.pieces.items()},
                     info={"largestPieceKh": kh}, n_extreme=ext)
        if not converged:
            last = history[-1]
            out.warnings.append(
                f"second-order iteration did not converge in {len(history)} iteration(s): last displacement change "
                f"{last['displacementChange']:.2e}, axial-force change {last['axialForceChange']:.2e} "
                f"(tolerance {tol:.1e})")
        return out

    # -------------------------------------------------------- beyond the critical load
    def _critical_factor(self, combination, crit: _Critical, slack):
        """Largest factor lam <= 1 on the axial forces of the last iterate for which the tangent stiffness
        is still positive definite (bisection; exact pieces, so no linearisation in lam)."""
        mesh, rloads, Ns, tau = crit.state
        rmodel = mesh.model
        owner = {pid: mid for mid, ids in mesh.pieces.items() for pid in ids}
        tolp = float(self.model.analysis.pivot_tolerance)

        def stable(lam: float) -> bool:
            elems = {}
            for pid, pm in rmodel.members.items():
                mid = owner[pid]
                fct = self._factors(mid, combination, tau)
                if mid in slack:
                    elems[pid] = build_element(rmodel, pm, fct, active=False)
                    continue
                ax = axial_state(rloads.member.get(pid) or SpanLoads(), Ns[pid], float(pm.length))
                elems[pid] = BeamColumn.build(rmodel, pm, None, 0.0, fct, axial=tuple(lam * a for a in ax))
                if not elems[pid].stable:
                    return False
            st = Structure(rmodel, elems)
            K = st.assemble()
            K[np.diag_indices(st.n)] += st.spring
            Kff = st.reduce(K)[np.ix_(st.free_r, st.free_r)]
            d = np.diag(Kff)
            if not (d > 0.0).all():
                return False
            s = np.sqrt(d)
            try:
                c = np.linalg.cholesky(Kff / np.outer(s, s))
            except np.linalg.LinAlgError:
                return False
            return float(np.min(np.diag(c)) ** 2) >= tolp

        if not stable(0.0):
            return None
        lo, hi = 0.0, 1.0
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if stable(mid):
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def _refuse_critical(self, combination, crit: _Critical, slack):
        lam = self._critical_factor(combination, crit, slack)
        if lam is None and crit.cause is not None:
            raise crit.cause                                 # not a buckling problem: a mechanism
        factor = "unknown" if lam is None else f"{lam:.4g}"
        basis = ("the first-order axial forces" if crit.iteration <= 1
                 else f"the axial forces of second-order iteration {crit.iteration - 1}")
        msg = (f"combination {combination.id!r}: the load exceeds the elastic critical load ({crit.what}); "
               f"buckling factor ~ {factor} (the largest multiple of {basis} at which the stiffness is still "
               "positive definite; a value below 1 means the loads are beyond the critical level)")
        raise CombinationRefused([Problem(
            pointer("combinations", combination.id), msg,
            "Reduce the loads, stiffen or brace the structure; run a buckling analysis (analysis.buckling) of "
            "this combination to see the mode.", code=BEYOND)],
            info={"bucklingFactor": lam, "iteration": crit.iteration})

    # ------------------------------------------------------------------ active set
    def solve(self, combination) -> CombinationResult:
        """Analyse one combination (any analysis type) with all nonlinear features."""
        second = combination.analysis_type == "secondOrder"
        if second:
            self._check_second_order_scope(combination)
        loads, alpha, notional = self.loads(combination, second)
        state, slack, active_info = self.active_set(combination, loads, second)
        return self._package(combination, state, slack, alpha, notional, active_info, loads)

    def _check_second_order_scope(self, combination, what="second-order analysis"):
        problems = []
        for mid, m in self.model.members.items():
            if m.type in TRUSS_TYPES:
                continue
            kk = list(m.k_end_i) + list(m.k_end_j)
            bad = [("i" if a < 6 else "j") + "." + ("ux", "uy", "uz", "rx", "ry", "rz")[a % 6]
                   for a in (1, 2, 7, 8) if kk[a] != math.inf]
            if bad:
                problems.append(Problem(
                    pointer("members", mid, "releases"),
                    f"{what} of combination {combination.id!r}: member {mid!r} has a released or spring-connected "
                    f"transverse translation ({', '.join(bad)}); the moment that the axial force produces across "
                    "such a sliding connection depends on how the connection is built, which the model does not say",
                    "Use rotational or axial releases, or model the sliding connection explicitly.",
                    code="analysis.release-second-order"))
        if problems:
            raise CombinationRefused(problems)

    def inner(self, combination, loads, slack, second: bool) -> _State:
        if not second:
            return self.first_order(combination, loads, slack)
        try:
            return self.second_order(combination, loads, slack)
        except _Critical as crit:
            self._refuse_critical(combination, crit, slack)
            raise                                            # pragma: no cover

    def violations(self, state: _State, slack) -> dict:
        """One-sided members whose state contradicts their mean axial force: {member: |force|}."""
        model = self.model
        strain = {mid: _strain(m, model, state.u) for mid, m in model.members.items()}
        ref = max((abs(e) for e in strain.values()), default=0.0)
        out = {}
        for mid, m in model.members.items():
            if m.type not in ONE_SIDED:
                continue
            e = strain[mid] if m.type == "tensionOnly" else -strain[mid]     # e > 0: the member works
            if abs(e) <= SLACK_TOL * ref:
                continue
            if (mid in slack) == (e > 0.0):
                out[mid] = abs(e) * model.materials[m.material].E * model.sections[m.section].A
        return out

    def active_set(self, combination, loads, second: bool):
        model = self.model
        one_sided = [mid for mid, m in model.members.items() if m.type in ONE_SIDED]
        itmax = int(model.analysis.max_iterations)
        cid = combination.id
        mode = "all"
        slack: set = set()
        seen = [frozenset()]
        trail = []
        it = 0
        while True:
            it += 1
            try:
                state = self.inner(combination, loads, slack, second)
            except ModelUnstable as e:
                if not slack:
                    raise
                if mode == "all":
                    mode, slack, seen = "single", set(), [frozenset()]
                    trail.append({"restart": "one member at a time (a mechanism appeared)"})
                    continue
                names = ", ".join(sorted(slack))
                diag = getattr(e, "diagnosis", None)
                raise CombinationRefused(
                    [Problem(pointer("combinations", cid),
                             f"combination {cid!r}: with the tension-only / compression-only member(s) {names} slack "
                             "(their axial force would have the wrong sign), the structure is a mechanism: "
                             + (diag.summary if diag is not None else str(e)),
                             "The structure relies on members that do not work under this loading; add members "
                             "that can take the opposite force or change the member type.", code=SLACK_MECHANISM)]
                    + list(e.problems), diagnosis=diag, info={"slack": sorted(slack)}) from None
            viol = self.violations(state, slack) if one_sided else {}
            trail.append({"iteration": it, "slack": sorted(slack), "violations": sorted(viol)})
            if not viol:
                break
            flip = set(viol) if mode == "all" else {max(sorted(viol), key=lambda k: viol[k])}
            new = set(slack) ^ flip
            if frozenset(new) in seen or it >= itmax:
                if mode == "all" and it < itmax:
                    mode, slack, seen = "single", set(), [frozenset()]
                    trail.append({"restart": "one member at a time (the active set cycled)"})
                    continue
                why = "cycles" if frozenset(new) in seen else f"did not settle in {itmax} iterations"
                raise CombinationRefused([Problem(
                    pointer("combinations", cid),
                    f"combination {cid!r}: the set of slack tension-only / compression-only members {why}; "
                    f"members changing state: {', '.join(sorted(flip))}; last slack set: "
                    f"{', '.join(sorted(slack)) or 'none'}",
                    "The loading leaves these members at the limit between working and slack (or the structure is "
                    "close to a mechanism); check the bracing layout.", code="solve.active-set-cycling")],
                    info={"history": trail})
            slack = new
            seen.append(frozenset(new))
        if one_sided:
            self._check_unique(combination, state, slack)
        info = {"iterations": it, "slack": sorted(slack), "history": trail, "strategy": mode} if one_sided else None
        return state, slack, info

    def _check_unique(self, combination, state: _State, slack):
        """A converged active set can still be a point of a whole range of equilibrium states: when
        one-sided members sit exactly at zero force and the structure without them (and without the
        slack ones) is a mechanism that can move without stretching any of them.  Then the result is
        not determined and the combination is refused (THEORY 12)."""
        from scipy.linalg import eigh
        model = self.model
        strain = {mid: _strain(m, model, state.u) for mid, m in model.members.items()}
        ref = max((abs(e) for e in strain.values()), default=0.0)
        idle = [mid for mid, m in model.members.items()
                if m.type in ONE_SIDED and mid not in slack and abs(strain[mid]) <= SLACK_TOL * ref]
        if not idle:
            return
        la = self._linear(combination, set(slack) | set(idle), state.tau)
        try:
            la.factor
            return                                           # stable without the idle members: unique
        except ModelUnstable as e:
            unstable = e
        st = la.structure
        f = st.free_r
        Kff = st.reduce(la.K)[np.ix_(f, f)]
        d = np.sqrt(np.where(np.diag(Kff) > 0.0, np.diag(Kff), 1.0))
        lam, vec = eigh(Kff / np.outer(d, d))
        null = [k for k in range(len(lam)) if lam[k] <= float(model.analysis.pivot_tolerance) * max(lam[-1], 1.0)]
        if len(null) == 1:
            ur = np.zeros(len(st.retained))
            ur[f] = vec[:, null[0]] / d
            v = st.T @ ur
            engage = []
            for mid in idle:
                g = _strain(model.members[mid], model, v) * (1.0 if model.members[mid].type == "tensionOnly" else -1.0)
                engage.append(g)
            top = max((abs(g) for g in engage), default=0.0)
            if any(g > 1.0e-9 * top for g in engage) and any(g < -1.0e-9 * top for g in engage):
                return                                       # either sense of the motion loads one of them: unique
        diag = getattr(unstable, "diagnosis", None)
        cid = combination.id
        raise CombinationRefused(
            [Problem(pointer("combinations", cid),
                     f"combination {cid!r}: the tension-only / compression-only member(s) "
                     f"{', '.join(sorted(slack)) or 'none'} are slack and {', '.join(sorted(idle))} carry no force; "
                     "without them the structure is a mechanism that can move without loading them, so the "
                     "displaced state is not determined: " + (diag.summary if diag is not None else str(unstable)),
                     "The structure relies on members that do not work under this loading; add members that can "
                     "take the opposite force or change the member type.", code=SLACK_MECHANISM)]
            + list(unstable.problems), diagnosis=diag, info={"slack": sorted(slack), "idle": sorted(idle)})

    # --------------------------------------------------------------------- result
    def _package(self, combination, state: _State, slack, alpha, notional, active_info, loads) -> CombinationResult:
        model = self.model
        a = float(alpha)
        members = {mid: MemberResult(r.end_forces / a, np.asarray(r.x, float).copy(), r.forces / a,
                                     r.displacements[:, :3] / a, active=mid not in slack)
                   for mid, r in state.members.items()}
        solver: dict = {}
        if combination.analysis_type == "secondOrder":
            solver["secondOrder"] = {
                "iterations": state.iterations, "converged": state.converged,
                "tolerance": float(model.analysis.convergence_tolerance),
                "maxIterations": int(model.analysis.max_iterations), "history": state.history,
                "subdivisions": dict(state.subdivisions), **state.info}
        elif len(state.history) > 1:
            solver["tauBIterations"] = {"iterations": state.iterations, "converged": state.converged,
                                        "history": state.history}
        if active_info:
            solver["activeSet"] = active_info
        if combination.stiffness_reduction:
            solver["stiffnessReduction"] = {
                "clause": "AISC 360-22 C2.3(a), C2.3(b) Eq. C2-2a / C2-2b",
                "factor": float(model.analysis.reduction_factor), "tauB": model.analysis.tau_b,
                "members": sorted(state.tau)}
        warnings = list(state.warnings)
        if notional is not None:
            u = model.units
            solver["notionalLoads"] = {
                "clause": "AISC 360-22 C2.2b, Eq. C2-1", "coefficient": combination.notional.coefficient,
                "direction": [float(x) for x in combination.notional.direction], "unit": u.label("force"),
                "level": "analysis loads (already multiplied by alpha)",
                "nodes": {nid: [float(u.from_internal(v, "force")) for v in notional[n.index]]
                          for nid, n in model.nodes.items() if notional[n.index].any()}}
        if combination.analysis_type == "secondOrder" and model.constraints:
            warnings.append("rigid links / diaphragms are linear kinematic constraints: the second-order effect of "
                            "the forces they transmit is not included")
        return CombinationResult(
            combination.id, combination.analysis_type, state.u.reshape(-1, 6) / a, state.reactions.reshape(-1, 6) / a,
            members, converged=state.converged, iterations=state.iterations, alpha=a, tau_b=dict(state.tau),
            invariants=InvariantsReport(list(state.checks)), warnings=warnings,
            constraint_forces=state.constraint_forces.reshape(-1, 6) / a, solver=solver)


@dataclass
class _Recovery:
    end_forces: np.ndarray
    x: np.ndarray
    forces: np.ndarray
    displacements: np.ndarray


def _rel(err: float, scale: float) -> float:
    return float(err / scale) if scale > 0.0 else (0.0 if err == 0.0 else math.inf)


def _strain(m, model, u: np.ndarray) -> float:
    """Elongation / length of a member from the node displacements (first-order kinematics, rigid
    offsets included) = its mean axial force / EA when it is active."""
    ends = []
    for nid, off in ((m.i, m.offset_i), (m.j, m.offset_j)):
        k = 6 * model.nodes[nid].index
        ends.append(u[k:k + 3] + np.cross(u[k + 3:k + 6], off))
    return float((ends[1] - ends[0]) @ np.asarray(m.R, float)[0]) / float(m.length)


def _axial_after(span: SpanLoads, n0: float, x: float) -> float:
    """Axial force just after x in a member whose axial force just before x = 0 is ``n0`` (statics)."""
    n = float(n0)
    for xp, f in span.points:
        if xp <= x:
            n -= float(f[0])
    for x1, x2, w1, w2 in span.dists:
        if x > x1 and x2 > x1:
            c = min(x, x2) - x1
            n -= float(w1[0]) * c + 0.5 * float(w2[0] - w1[0]) / (x2 - x1) * c * c
    return n


def _difference(a: _State, b: _State, size: float) -> float:
    """Largest relative change between two solutions of the same combination (node displacements and
    station forces / moments of every member, each group relative to its own largest value)."""
    out = 0.0
    ua, ub = a.u.reshape(-1, 6), b.u.reshape(-1, 6)
    for sl in (slice(0, 3), slice(3, 6)):
        out = max(out, _rel(float(np.abs(ua[:, sl] - ub[:, sl]).max()), float(np.abs(ub[:, sl]).max())))
    fa = np.concatenate([r.forces for r in a.members.values()])
    fb = np.concatenate([r.forces for r in b.members.values()])
    for sl in (slice(0, 3), slice(3, 6)):
        out = max(out, _rel(float(np.abs(fa[:, sl] - fb[:, sl]).max()), float(np.abs(fb[:, sl]).max())))
    return out


def second_order_checks(structure, K_members, factor, loads: LoadSet, solution, elems, recs) -> list[InvariantCheck]:
    """Invariants of a second-order solution, evaluated with the terms of the theory (THEORY 15):
    forces are components in the undeformed axes; moments include, for every piece, the couple of its
    axial force acting through the transverse displacements, and the rigid-arm terms."""
    model = structure.model
    tol = float(model.analysis.equilibrium_tolerance)
    n = structure.n
    Lc = structure.size if structure.size > 0.0 else 1.0
    out = []
    kmax = float(np.abs(K_members).max()) if K_members.size else 0.0
    sym = inv._ratio(float(np.abs(K_members - K_members.T).max()) if K_members.size else 0.0, kmax)
    out.append(InvariantCheck("stiffness symmetry", sym, inv.SYMMETRY_TOL, sym <= inv.SYMMETRY_TOL,
                              "max |K - K^T| / max |K| of the assembled tangent stiffness (elastic + geometric)"))
    B = inv.rigid_body_modes(structure.xyz)[:, :3]
    trans = (np.arange(n) % 6) < 3
    num, den = np.abs(K_members @ B), np.abs(K_members) @ np.abs(B)
    rb = 0.0
    for rows in (trans, ~trans):
        if rows.any():
            rb = max(rb, inv._ratio(float(num[rows].max()), float(den[rows].max())))
    out.append(InvariantCheck("rigid-body modes", rb, inv.RIGID_BODY_TOL, rb <= inv.RIGID_BODY_TOL,
                              "max |K B| / max (|K| |B|) for the 3 rigid-body TRANSLATIONS (a rigid rotation of a "
                              "member that carries an axial force does produce forces in second-order theory)"))
    ptol = float(model.analysis.pivot_tolerance)
    out.append(InvariantCheck(
        "positive definiteness", factor.lambda_min, ptol, factor.lambda_min >= ptol and factor.min_pivot >= ptol,
        "smallest eigenvalue of the free-free TANGENT stiffness (elastic + geometric) scaled to unit diagonal "
        f"({'exact' if factor.exact else 'estimated'}), must be >= tolerance; smallest pivot {factor.min_pivot:.3e}; "
        f"{len(factor.free)} free DOFs (internal nodes of subdivided members included)"))
    out.append(InvariantCheck("linear solve residual", solution.residual, inv.RESIDUAL_TOL,
                              solution.residual <= inv.RESIDUAL_TOL,
                              "backward error of the last solve: max |K_ff u_f - rhs| / max (|K_ff| |u_f| + |rhs|)"))
    # ---- global equilibrium
    R, A = applied_resultant(model, loads)
    react = solution.reactions.reshape(-1, 6)
    Rr, Ar = np.zeros(6), np.zeros(6)
    for k in range(len(react)):
        f, mo = react[k, :3], react[k, 3:]
        arm = np.cross(structure.xyz[k], f)
        Rr[:3] += f
        Rr[3:] += mo + arm
        Ar[:3] += np.abs(f)
        Ar[3:] += np.abs(mo) + np.abs(arm)
    geo, geo_abs = np.zeros(3), np.zeros(3)
    u = solution.u
    g = loads.nodal + solution.reactions + solution.constraint_forces
    absum = np.abs(loads.nodal) + np.abs(solution.reactions) + np.abs(solution.constraint_forces)
    fc = mc = dc = rc = fs = ms = ds = rs = 0.0
    dres = dscale = 0.0
    for pid, e in elems.items():
        rec = recs[pid]
        if isinstance(e, BeamColumn):
            fn = rec.node_forces
            arms = e.arm[0] @ u[e.dofs[3:6]] + e.arm[1] @ u[e.dofs[9:12]]
            geo += rec.couple + arms
            geo_abs += np.abs(rec.couple) + np.abs(arms)
            dres, dscale = max(dres, rec.deformed_residual), max(dscale, rec.moment_scale)
            ue = e.G @ u[e.dofs]
            ds, rs = max(ds, float(np.abs(ue[[0, 1, 2, 6, 7, 8]]).max())), max(rs, float(np.abs(ue[[3, 4, 5, 9, 10, 11]]).max()))
        else:
            fn = e.G.T @ rec.end_forces
            ds, rs = max(ds, float(np.abs(rec.displacements[:, :3]).max())), max(rs, float(np.abs(rec.displacements[:, 3:]).max()))
        g[e.dofs] -= fn
        absum[e.dofs] += np.abs(fn)
        fc, mc = max(fc, float(np.abs(rec.force_closure[:3]).max())), max(mc, float(np.abs(rec.force_closure[3:]).max()))
        dc, rc = max(dc, float(np.abs(rec.disp_closure[:3]).max())), max(rc, float(np.abs(rec.disp_closure[3:]).max()))
        fs, ms = max(fs, float(np.abs(rec.end_forces[[0, 1, 2, 6, 7, 8]]).max())), max(ms, float(np.abs(rec.end_forces[[3, 4, 5, 9, 10, 11]]).max()))
    tot = R + Rr
    tot[3:] -= geo
    fscale = max(float(A[:3].max()), float(Ar[:3].max()))
    mscale = max(float(A[3:].max()), float(Ar[3:].max()), float(geo_abs.max()))
    fscale, mscale = max(fscale, mscale / Lc), max(mscale, fscale * Lc)
    v = inv._ratio(float(np.abs(tot[:3]).max()), fscale)
    out.append(InvariantCheck("global force equilibrium", v, tol, v <= tol,
                              "sum of applied loads + reactions (components in the undeformed axes)"))
    v = inv._ratio(float(np.abs(tot[3:]).max()), mscale)
    out.append(InvariantCheck(
        "global moment equilibrium", v, tol, v <= tol,
        "second-order: applied loads + reactions about the global origin at the undeformed positions, minus the "
        "couples of the member axial forces acting through the transverse displacements of the member axis "
        "(integral of N x slope over every piece) and of the rigid end offsets"))
    fmax, mmax = float(absum[trans].max(initial=0.0)), float(absum[~trans].max(initial=0.0))
    fmax, mmax = max(fmax, mmax / Lc), max(mmax, fmax * Lc)
    jf = inv._ratio(float(np.abs(g[trans]).max(initial=0.0)), fmax)
    jm = inv._ratio(float(np.abs(g[~trans]).max(initial=0.0)), mmax)
    worst = int(np.argmax(np.abs(g) / np.where(trans, max(fmax, 1e-300), max(mmax, 1e-300)))) if n else 0
    v = max(jf, jm)
    out.append(InvariantCheck("joint equilibrium", v, tol, v <= tol,
                              "nodal loads + reactions + constraint forces - second-order member end forces "
                              "(rigid-arm terms included) at every node, internal nodes included; largest residual "
                              f"at node {structure.node_ids[worst // 6] if n else '-'}"))
    fs, ms = max(fs, ms / Lc), max(ms, fs * Lc)
    ds, rs = max(ds, rs * Lc), max(rs, ds / Lc)
    v = max(inv._ratio(fc, fs), inv._ratio(mc, ms))
    out.append(InvariantCheck("member equilibrium", v, tol, v <= tol,
                              "resultants carried from end i to end j of every piece by the beam-column transfer "
                              "minus the end forces at j"))
    v = max(inv._ratio(dc, ds), inv._ratio(rc, rs))
    out.append(InvariantCheck("member compatibility", v, tol, v <= tol,
                              "displacements carried from end i to end j of every piece minus the displacement of end j"))
    v = inv._ratio(dres, max(dscale, ms * 1.0e-300))
    out.append(InvariantCheck("member equilibrium on the deformed axis", v, tol, v <= tol,
                              "moment equilibrium of every piece with its axial force acting on the deflected axis: "
                              "end moments + end shears x length + load moments - integral of N x slope, the "
                              "integral evaluated from the recovered deflection"))
    U = float(u @ (K_members @ u))
    terms = (solution.F + solution.reactions + solution.constraint_forces) * u
    W = float(terms.sum())
    v = inv._ratio(abs(U - W), max(abs(U), float(np.abs(terms).sum())))
    out.append(InvariantCheck("work = strain energy", v, tol, v <= tol,
                              "u^T K_T u versus (equivalent nodal loads + reactions) . u, K_T = tangent stiffness: "
                              "strain energy plus the second-order work of the axial forces"))
    return out
