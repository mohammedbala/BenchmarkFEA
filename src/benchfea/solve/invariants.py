"""Self-checks computed on every solve, and linearity checks callable by tests and ``verify``.

Every check is an :class:`benchfea.results.InvariantCheck` (name, value, tolerance, passed, detail).

Model level (once per stiffness matrix)
    stiffness symmetry              max |K - K^T| / max |K|                              <= 1e-10
    rigid-body modes                max |K_free B| / max (|K_free| |B|) over the 6 rigid-body
                                    motions B of the unsupported, unconstrained structure  <= 1e-9
    positive definiteness           smallest eigenvalue of the free-free matrix scaled to unit
                                    diagonal                                               >= analysis.pivotTolerance
Per load vector (tolerance = analysis.equilibriumTolerance unless stated).  Force residuals are
relative to the largest force-type quantity in the check (at least moment scale / model size), moment
residuals to the largest moment-type quantity (at least force scale x model size), so that quantities
that are zero by symmetry do not produce a 0/0 ratio; "model size" = diagonal of the bounding box of the nodes.
    linear solve residual           max |K_ff u_f - rhs| / max (|K_ff| |u_f| + |rhs|)      <= 1e-9
    global force equilibrium        applied loads (from their definition) + reactions
    global moment equilibrium       the same, moments about the global origin
    joint equilibrium               at every node: nodal loads + reactions + constraint forces
                                    - member end forces (recovered element by element)
    member equilibrium              resultants integrated along each member from end i reproduce
                                    the end forces at j
    member compatibility            displacements integrated along each member from end i
                                    (through the end springs) reproduce the displacement of end j
    work = strain energy            u^T K_members u = (loads + reactions) . u
"""
from __future__ import annotations

import numpy as np

from ..results import InvariantCheck

__all__ = ["model_checks", "solution_checks", "superposition_check", "reciprocity_check",
           "SYMMETRY_TOL", "RIGID_BODY_TOL", "RESIDUAL_TOL"]

SYMMETRY_TOL = 1.0e-10
RIGID_BODY_TOL = 1.0e-9
RESIDUAL_TOL = 1.0e-9


def _ratio(err: float, scale: float) -> float:
    return float(err / scale) if scale > 0.0 else (0.0 if err == 0.0 else float("inf"))


def rigid_body_modes(xyz: np.ndarray) -> np.ndarray:
    """(6 n, 6): unit translations along and unit rotations about the global axes through the origin."""
    n = len(xyz)
    B = np.zeros((6 * n, 6))
    for k in range(n):
        for a in range(3):
            e = np.eye(3)[a]
            B[6 * k + a, a] = 1.0
            B[6 * k:6 * k + 3, 3 + a] = np.cross(e, xyz[k])
            B[6 * k + 3 + a, 3 + a] = 1.0
    return B


def model_checks(structure, K_members: np.ndarray, factor) -> list[InvariantCheck]:
    out = []
    kmax = float(np.abs(K_members).max()) if K_members.size else 0.0
    sym = _ratio(float(np.abs(K_members - K_members.T).max()) if K_members.size else 0.0, kmax)
    out.append(InvariantCheck("stiffness symmetry", sym, SYMMETRY_TOL, sym <= SYMMETRY_TOL,
                              "max |K - K^T| / max |K| of the assembled member stiffness"))
    B = rigid_body_modes(structure.xyz)
    num = np.abs(K_members @ B)
    den = np.abs(K_members) @ np.abs(B)
    trans = (np.arange(len(B)) % 6) < 3
    rb = 0.0
    for rows in (trans, ~trans):                 # force rows and moment rows have different units
        if rows.any():
            rb = max(rb, _ratio(float(num[rows].max()), float(den[rows].max())))
    out.append(InvariantCheck("rigid-body modes", rb, RIGID_BODY_TOL, rb <= RIGID_BODY_TOL,
                              "max |K B| / max (|K| |B|), force rows and moment rows separately, for the 6 "
                              "rigid-body motions B of the unsupported, unconstrained structure (K B = 0 means "
                              "zero forces and zero strain energy)"))
    tol = float(structure.model.analysis.pivot_tolerance)
    out.append(InvariantCheck(
        "positive definiteness", factor.lambda_min, tol, factor.lambda_min >= tol and factor.min_pivot >= tol,
        f"Cholesky factorisation succeeded; value = smallest eigenvalue of the free-free stiffness scaled to unit "
        f"diagonal ({'exact' if factor.exact else 'estimated'}), must be >= tolerance; smallest pivot "
        f"{factor.min_pivot:.3e}; condition estimate {factor.condition:.3e}; {len(factor.free)} free DOFs"))
    return out


def solution_checks(structure, K_members, loads, solution, recoveries, applied) -> list[InvariantCheck]:
    """Checks of one solved load vector.  ``recoveries``: member id -> MemberRecovery; ``applied`` =
    ``loads.applied_resultant(model, loads)``."""
    model = structure.model
    tol = float(model.analysis.equilibrium_tolerance)
    n = structure.n
    out = [InvariantCheck("linear solve residual", solution.residual, RESIDUAL_TOL, solution.residual <= RESIDUAL_TOL,
                          "backward error of the solve: max |K_ff u_f - rhs| / max (|K_ff| |u_f| + |rhs|)")]
    # ---- global equilibrium: applied loads (from their definitions) + reactions
    R, A = applied
    react = solution.reactions.reshape(-1, 6)
    Rr, Ar = np.zeros(6), np.zeros(6)
    for k in range(len(react)):
        f, mo = react[k, :3], react[k, 3:]
        arm = np.cross(structure.xyz[k], f)
        Rr[:3] += f
        Rr[3:] += mo + arm
        Ar[:3] += np.abs(f)
        Ar[3:] += np.abs(mo) + np.abs(arm)
    tot = R + Rr
    Lc = structure.size if structure.size > 0.0 else 1.0
    fscale = max(float(A[:3].max()), float(Ar[:3].max()))
    mscale = max(float(A[3:].max()), float(Ar[3:].max()))
    fscale, mscale = max(fscale, mscale / Lc), max(mscale, fscale * Lc)
    for name, sl in (("global force equilibrium", slice(0, 3)), ("global moment equilibrium", slice(3, 6))):
        scale = mscale if sl.start else fscale
        v = _ratio(float(np.abs(tot[sl]).max()), scale)
        out.append(InvariantCheck(name, v, tol, v <= tol,
                                  "sum of applied loads + reactions" + (" about the global origin" if sl.start else "")
                                  + ", relative to the largest sum of absolute contributions"
                                  + (" (at least force scale x model size)" if sl.start
                                     else " (at least moment scale / model size)")))
    # ---- joint equilibrium from the recovered member end forces
    g = loads.nodal + solution.reactions + solution.constraint_forces
    absum = np.abs(loads.nodal) + np.abs(solution.reactions) + np.abs(solution.constraint_forces)
    fmax = mmax = 0.0
    for mid, rec in recoveries.items():
        e = structure.elements[mid]
        fn = e.G.T @ rec.end_forces                          # forces on the member at its nodes, global
        g[e.dofs] -= fn
        absum[e.dofs] += np.abs(fn)
    trans = (np.arange(n) % 6) < 3
    fmax, mmax = float(absum[trans].max(initial=0.0)), float(absum[~trans].max(initial=0.0))
    fmax, mmax = max(fmax, mmax / Lc), max(mmax, fmax * Lc)
    jf, jm = _ratio(float(np.abs(g[trans]).max(initial=0.0)), fmax), _ratio(float(np.abs(g[~trans]).max(initial=0.0)), mmax)
    worst = int(np.argmax(np.abs(g) / np.where(trans, max(fmax, 1e-300), max(mmax, 1e-300)))) if n else 0
    v = max(jf, jm)
    out.append(InvariantCheck("joint equilibrium", v, tol, v <= tol,
                              "nodal loads + reactions + constraint forces - member end forces at every node; "
                              f"largest residual at node {structure.node_ids[worst // 6] if n else '-'}"))
    # ---- member closures
    fc = mc = dc = rc = 0.0
    fs = ms = ds = rs = 0.0
    for mid, rec in recoveries.items():
        fc, mc = max(fc, float(np.abs(rec.force_closure[:3]).max())), max(mc, float(np.abs(rec.force_closure[3:]).max()))
        dc, rc = max(dc, float(np.abs(rec.disp_closure[:3]).max())), max(rc, float(np.abs(rec.disp_closure[3:]).max()))
        fs, ms = max(fs, float(np.abs(rec.forces[:, :3]).max())), max(ms, float(np.abs(rec.forces[:, 3:]).max()))
        ds, rs = max(ds, float(np.abs(rec.displacements[:, :3]).max())), max(rs, float(np.abs(rec.displacements[:, 3:]).max()))
    fs, ms = max(fs, ms / Lc), max(ms, fs * Lc)
    ds, rs = max(ds, rs * Lc), max(rs, ds / Lc)
    v = max(_ratio(fc, fs), _ratio(mc, ms))
    out.append(InvariantCheck("member equilibrium", v, tol, v <= tol,
                              "resultants integrated from end i to end j minus the end forces at j"))
    v = max(_ratio(dc, ds), _ratio(rc, rs))
    out.append(InvariantCheck("member compatibility", v, tol, v <= tol,
                              "displacements integrated from end i (through the end springs) minus the "
                              "displacement of end j"))
    # ---- work = strain energy
    u = solution.u
    U = float(u @ (K_members @ u))
    terms = (solution.F + solution.reactions + solution.constraint_forces) * u
    W = float(terms.sum())
    v = _ratio(abs(U - W), max(abs(U), float(np.abs(terms).sum())))
    out.append(InvariantCheck("work = strain energy", v, tol, v <= tol,
                              "u^T K_members u versus (equivalent nodal loads + reactions) . u"))
    return out


def superposition_check(model, tolerance: float = 1.0e-9) -> list[InvariantCheck]:
    """Linearity: for every first-order combination, the factored sum of the separately solved load
    cases equals the directly solved combination (displacements, reactions, station forces and station
    displacements).  One check per combination."""
    from .analysis import LinearAnalysis
    la = LinearAnalysis(model)
    out = []
    for cid, combo in model.combinations.items():
        direct = la.solve_combination(combo)
        acc = None
        for case_id, f in combo.factors.items():
            sol, rec = la.solve_loads(la.case(case_id))
            parts = [sol.u, sol.reactions] + [x for r in rec.values() for x in (r.forces.ravel(), r.displacements.ravel())]
            acc = [f * p for p in parts] if acc is None else [a + f * p for a, p in zip(acc, parts)]
        sol, rec = direct
        ref = [sol.u, sol.reactions] + [x for r in rec.values() for x in (r.forces.ravel(), r.displacements.ravel())]
        v = 0.0
        for a, b in zip(acc or [], ref):
            v = max(v, _ratio(float(np.abs(a - b).max(initial=0.0)), float(np.abs(b).max(initial=0.0))))
        out.append(InvariantCheck(f"superposition: {cid}", v, tolerance, v <= tolerance,
                                  "sum of factor x case results versus the combination solved directly"))
    return out


def reciprocity_check(model, pairs=None, tolerance: float = 1.0e-9, max_pairs: int = 12) -> InvariantCheck:
    """Maxwell-Betti: the displacement at DOF b due to a unit load at DOF a equals the displacement at
    a due to a unit load at b.  ``pairs`` = iterable of ((node, dof), (node, dof)); by default a
    deterministic selection of free DOFs spread over the model.  The value is the largest
    |d_ab - d_ba| / sqrt(d_aa d_bb)."""
    from .analysis import LinearAnalysis
    la = LinearAnalysis(model)
    st = la.structure
    if pairs is None:
        free = st.retained[st.free_r]
        if len(free) < 2:
            return InvariantCheck("reciprocity (Maxwell-Betti)", 0.0, tolerance, True, "fewer than two free DOFs")
        pick = free[np.unique(np.linspace(0, len(free) - 1, min(len(free), 2 * max_pairs)).astype(int))]
        gp = [(int(pick[k]), int(pick[-1 - k])) for k in range(len(pick) // 2)]
    else:
        gp = [(model.dof(*a), model.dof(*b)) for a, b in pairs]
    cols = {}
    for g in sorted({g for pair in gp for g in pair}):
        F = np.zeros(st.n)
        F[g] = 1.0
        cols[g] = st.solve(la.factor, F).u
    v, worst = 0.0, None
    for a, b in gp:
        dab, dba = cols[a][b], cols[b][a]
        r = _ratio(abs(dab - dba), float(np.sqrt(abs(cols[a][a] * cols[b][b]))))
        if r >= v:
            v, worst = r, (a, b)
    names = "" if worst is None else " and ".join(f"{st.node_ids[g // 6]}.{'ux uy uz rx ry rz'.split()[g % 6]}"
                                                    for g in worst)
    return InvariantCheck("reciprocity (Maxwell-Betti)", v, tolerance, v <= tolerance,
                          f"{len(gp)} DOF pairs; largest |d_ab - d_ba| / sqrt(d_aa d_bb) for {names}")
