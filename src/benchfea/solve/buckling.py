"""Linear (eigenvalue) buckling: (K + lambda Kg) phi = 0 for a reference combination.

Theory: docs/THEORY.md section 11.  Procedure:

1. Reference state = FIRST-ORDER solution of the combination's factored loads (not multiplied by
   ``analysis.alpha``; notional loads, stiffness reduction and the slack set of tension/compression-
   only members as the combination defines them).  It gives the axial force N(x) of every member.
2. Every member is cut into n pieces (``refine``; load break points are piece boundaries).  K is the elastic stiffness
   of the pieces (exact), Kg their consistent geometric stiffness for the reference axial forces
   (``beamcolumn.consistent_geometric_stiffness``), which is linear in the load factor.
3. The generalised symmetric eigenproblem is solved on the free DOFs (constraints and supports
   eliminated exactly as in the static analysis).  lambda > 0: buckling under lambda x the loads;
   lambda < 0: buckling under the reversed loads (reported separately).
4. Each eigenpair is then taken to the eigenpair of K_T(lambda) phi = 0 with the EXACT beam-column
   pieces (``_polish``): this is the same problem without the discretisation error of the consistent
   geometric stiffness (which decreases with the fourth power of the piece length).
5. n = 4, 8, ... until the positive factors change by less than ``BUCKLING_TOL`` between two
   subdivisions (at least two subdivisions are always solved); the last change is reported as the
   estimated relative error.

Modes are normalised to a largest translation of 1 (over nodes and member axes) and the sign is
fixed by making the first translation component that reaches the largest magnitude positive.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh, lu_factor, lu_solve

from ..results import BucklingResult
from .beam import SpanLoads
from .beamcolumn import BeamColumn, axial_state, consistent_geometric_stiffness
from .element import build_element
from .nonlinear import _axial_after
from .refine import piece_stations, refine, split_loads
from .system import Structure

__all__ = ["buckling_analysis", "BUCKLING_TOL", "PIECES_START", "PIECES_MAX", "DENSE_LIMIT"]

BUCKLING_TOL = 1.0e-7      # relative change of the reported factors between the last two subdivisions
EXTRA_MODES = 1            # eigenpairs computed beyond the requested number (the ordering can change when polished)
POLISH_ITERATIONS = 12
POLISH_TOL = 1.0e-10       # relative change of the factor at which the polishing stops
POLISH_ACCEPT = 1.0e-9     # ... and up to which a polished factor is accepted after POLISH_ITERATIONS
POLISH_SHIFT = 2.0e-3
PIECES_START = 4
PIECES_MAX = 32
DENSE_LIMIT = 1500         # free DOFs up to which all eigenvalues are computed
MU_TOL = 1.0e-10           # eigenvalues 1/lambda below this fraction of the largest are "no buckling"
SIGN_TOL = 1.0e-9


def _level(nl, combination, loads, state, slack, n: int, n_modes: int, previous=None):
    model = nl.model
    mesh = refine(model, n, unsplit=slack)
    rmodel = mesh.model
    rloads = split_loads(mesh, loads)
    empty = SpanLoads()
    owner = {pid: mid for mid, ids in mesh.pieces.items() for pid in ids}
    elems, kg, axial, factors = {}, {}, {}, {}
    for mid, ids in mesh.pieces.items():
        fct = nl._factors(mid, combination, state.tau)
        if mid in slack:
            elems[mid] = build_element(rmodel, rmodel.members[mid], fct, active=False)
            continue
        n0 = -state.rec[mid].end_forces[0]
        span = loads.member.get(mid, empty)
        for k, pid in enumerate(ids):
            pm = rmodel.members[pid]
            elems[pid] = BeamColumn.build(rmodel, pm, None, 0.0, fct)
            ns = _axial_after(span, n0, float(mesh.bounds[mid][k]))
            axial[pid] = axial_state(rloads.member.get(pid) or empty, ns, float(pm.length))
            factors[pid] = fct
            kg[pid] = consistent_geometric_stiffness(elems[pid], axial[pid])
    st = Structure(rmodel, elems)
    KE = st.assemble()
    KE[np.diag_indices(st.n)] += st.spring
    st.factorize(KE)                                         # ModelUnstable if the elastic structure is a mechanism
    f = st.free_r
    Ke = st.reduce(KE)[np.ix_(f, f)]
    Kg = st.reduce(st.assemble(kg))[np.ix_(f, f)]
    d = np.sqrt(np.diag(Ke))
    A, B = -Kg / np.outer(d, d), Ke / np.outer(d, d)
    nf = len(f)
    want = n_modes + EXTRA_MODES
    if nf <= DENSE_LIMIT:
        mu, V = eigh(A, B)
    else:
        m = min(want, nf)
        mu_lo, V_lo = eigh(A, B, subset_by_index=[0, m - 1])
        mu_hi, V_hi = eigh(A, B, subset_by_index=[nf - m, nf - 1])
        mu, V = np.concatenate([mu_lo, mu_hi]), np.hstack([V_lo, V_hi])
    top = float(np.abs(mu).max()) if len(mu) else 0.0
    pos = [k for k in np.argsort(-mu) if mu[k] > MU_TOL * top][:want]
    neg = [k for k in np.argsort(mu) if mu[k] < -MU_TOL * top][:want]

    def tangent(lam: float):
        """Free-free tangent stiffness (scaled to the unit diagonal of K) for lam x the reference axial
        forces, with the EXACT beam-column pieces; None if a piece is beyond its own critical load."""
        mats, cache = {}, {}
        for pid, ax in axial.items():
            pm = rmodel.members[pid]
            # equal pieces of one member (same axial force, length, connections) share their matrix
            key = (owner[pid], ax, float(pm.length), pm.k_end_i, pm.k_end_j, pm.offset_i.any(), pm.offset_j.any())
            if key not in cache or key[-1] or key[-2]:
                piece = BeamColumn.build(rmodel, pm, None, 0.0, factors[pid], axial=tuple(lam * a for a in ax))
                if not piece.stable:
                    return None
                cache[key] = piece.k_global
            mats[pid] = cache[key]
        K = st.assemble(mats)
        K[np.diag_indices(st.n)] += st.spring
        return st.reduce(K)[np.ix_(f, f)] / np.outer(d, d)

    def solve(idx, known):
        """(factor, mode on all DOFs, polished?) for the eigenpairs ``idx`` of the linearised problem;
        ``known`` = polished factors of the previous subdivision (starting values)."""
        out = []
        for q, k in enumerate(idx):
            lam, phi, exact = 1.0 / mu[k], V[:, k] * 1.0, False
            start = known[q] if q < len(known) and abs(known[q] / lam - 1.0) < 0.05 else None
            better = _polish(tangent, B, lam, phi, start)
            if better is not None:
                lam, phi, exact = better[0], better[1], True
            ur = np.zeros(len(st.retained))
            ur[f] = phi / d
            out.append((lam, st.T @ ur, exact))
        return out
    pos_s = sorted(solve(pos, previous["pos"] if previous else []), key=lambda t: t[0])[:n_modes]
    neg_s = sorted(solve(neg[:n_modes], previous["neg"] if previous else []), key=lambda t: -t[0])
    return {"mesh": mesh, "elems": elems, "st": st, "nf": nf,
            "pos": [t[0] for t in pos_s], "neg": [t[0] for t in neg_s],
            "pos_modes": [t[1] for t in pos_s], "neg_modes": [t[1] for t in neg_s],
            "exact": all(t[2] for t in pos_s + neg_s)}


def _polish(tangent, B: np.ndarray, lam0: float, phi0: np.ndarray, start: float | None = None):
    """Take an eigenpair of the linearised problem (K + lam Kg) phi = 0 to the eigenpair of the exact
    problem K_T(lam) phi = 0 (the limit of the linearised problem for vanishing piece length; THEORY
    11.2) by residual inverse iteration: lam from the Rayleigh functional p(lam) = phi^T K_T(lam) phi = 0
    (Newton / secant; first slope -phi^T K phi / lam0 from the linearised problem),
    phi <- phi - K_T(sigma)^-1 K_T(lam) phi.  Returns (lam, phi) or None if it does not converge."""
    # factorise slightly below the eigenvalue: K_T(lam0) itself is nearly singular on a fine mesh
    K0 = tangent(lam0 * (1.0 - POLISH_SHIFT))
    if K0 is None:
        return None
    try:
        lu = lu_factor(K0)
    except (ValueError, np.linalg.LinAlgError):
        return None
    phi = phi0 / np.linalg.norm(phi0)
    lam = lam0 if start is None else start
    slope = -float(phi @ B @ phi) / lam0
    change = np.inf
    for _ in range(POLISH_ITERATIONS):
        a, Ka = lam, tangent(lam)
        if Ka is None or slope == 0.0:
            return None
        fa = float(phi @ Ka @ phi)
        b, Kb = a, Ka
        for _ in range(8):
            b = a - fa / slope
            if not 0.5 < b / lam0 < 2.0:                     # left the neighbourhood of this eigenvalue
                return None
            if abs(b - a) <= 1.0e-12 * abs(b):
                Kb = Ka
                break
            Kb = tangent(b)
            if Kb is None:
                return None
            fb = float(phi @ Kb @ phi)
            if fb != fa:
                slope = (fb - fa) / (b - a)
            a, fa, Ka = b, fb, Kb
        x = lu_solve(lu, Kb @ phi)
        if not np.all(np.isfinite(x)):
            return None
        new_phi = phi - x
        new_phi /= np.linalg.norm(new_phi)
        change = abs(b - lam) / abs(b)
        lam, phi = b, new_phi
        if change <= POLISH_TOL:
            return lam, phi
    return (lam, phi) if change <= POLISH_ACCEPT else None


def _normalise(u: np.ndarray) -> np.ndarray:
    t = u.reshape(-1, 6)[:, :3].ravel()
    top = float(np.abs(t).max())
    if top == 0.0:
        return u
    first = int(np.flatnonzero(np.abs(t) >= (1.0 - SIGN_TOL) * top)[0])
    return u / (top if t[first] > 0.0 else -top)


def buckling_analysis(nl, combination) -> BucklingResult:
    """Lowest ``analysis.buckling.modes`` positive buckling factors of ``combination`` and their modes."""
    model = nl.model
    la = nl.la
    n_modes = int(model.analysis.buckling_modes)
    nl._check_second_order_scope(combination, what="buckling analysis")
    loads, _, _ = nl.loads(combination, second_order=False)
    state, slack, active = nl.active_set(combination, loads, second=False)
    prev, level, n, est, est_neg = None, None, PIECES_START, None, []
    while True:
        level = _level(nl, combination, loads, state, slack, n, n_modes, prev)
        if prev is not None:
            est = [abs(a - b) / abs(a) for a, b in zip(level["pos"], prev["pos"])]
            est_neg = [abs(a - b) / abs(a) for a, b in zip(level["neg"], prev["neg"])]
            if max(est, default=0.0) <= BUCKLING_TOL or n >= PIECES_MAX:
                break
        prev, n = level, n * 2
    mesh, elems = level["mesh"], level["elems"]
    nb = mesh.n_dof_base
    pos = [_normalise(u) for u in level["pos_modes"]]
    neg = [_normalise(u) for u in level["neg_modes"]]
    empty = SpanLoads()
    member_shapes = {}
    for mid, ids in mesh.pieces.items():
        stations = la.stations[mid]
        arr = np.zeros((len(pos), len(stations), 3))
        for q, u in enumerate(pos):
            if mid in slack:
                e = elems[mid]
                arr[q] = e.recover(u[e.dofs], empty, stations).displacements[:, :3]
                continue
            where = piece_stations(mesh.bounds[mid], stations)
            for k, pid in enumerate(ids):
                e = elems[pid]
                r = e.recover(u[e.dofs], [(x, side) for _, x, side in where[k]])
                for row, (s, _, _) in enumerate(where[k]):
                    arr[q, s] = r.displacements[row, :3]
        member_shapes[mid] = arr
    warnings = []
    if not pos:
        warnings.append("no positive buckling factor: no member is in compression under this combination")
    if est is not None and max(est, default=0.0) > BUCKLING_TOL:
        warnings.append(f"buckling factors not converged to {BUCKLING_TOL:.0e} with {n} pieces per member "
                        f"(they still change by up to {max(est):.1e} between two subdivisions)")
    if model.constraints:
        warnings.append("rigid links / diaphragms are linear kinematic constraints: the geometric stiffness of the "
                        "forces they transmit is not included")
    info = {
        "basis": "first-order reference state of the combination's factored loads (not multiplied by alpha); "
                 "consistent geometric stiffness of the member axial forces; flexural buckling only",
        "piecesPerMember": n, "subdivisions": {mid: len(ids) for mid, ids in mesh.pieces.items()},
        "freeDofs": int(level["nf"]),
        "estimatedRelativeError": [float(e) for e in (est or [])],
        "estimatedRelativeErrorNegative": [float(e) for e in est_neg],
        "tolerance": BUCKLING_TOL, "exactMemberStiffness": bool(level["exact"]), "slack": sorted(slack), "tauB": dict(state.tau),
        "stiffnessReduction": bool(combination.stiffness_reduction), "warnings": warnings}
    if active:
        info["activeSet"] = active
    n_nodes = nb // 6

    def stack(modes):
        return np.array([u[:nb].reshape(-1, 6) for u in modes], float).reshape(len(modes), n_nodes, 6)
    return BucklingResult(combination.id, np.array(level["pos"], float), stack(pos),
                          negative_factors=np.array(level["neg"], float), negative_modes=stack(neg),
                          member_shapes=member_shapes, info=info)
