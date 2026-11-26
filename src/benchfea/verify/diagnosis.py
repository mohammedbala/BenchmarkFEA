"""Verification of a MECHANISM DIAGNOSIS (the reference solution of an unstable model).

For a model the engine refuses as unstable, the reference answer is the diagnosis
(``MechanismDiagnosis.to_dict()``): how many independent zero-stiffness motions there are and what they
are.  "Verifying" it means, here (docs/VERIFICATION.md, "Unstable models"):

tier N  null space of an INDEPENDENTLY ASSEMBLED stiffness matrix.  ``independent_stiffness`` builds the
        stiffness of the structure from the model data by a route that shares no code with
        ``benchfea.solve``: closed-form 12x12 member matrices (the engine inverts a cantilever
        flexibility), end connections by static condensation of series springs (the engine adds
        connection flexibilities), rigid offsets and constraints by its own kinematic relations.
        Checked: (1) the number of zero eigenvalues of the scaled matrix equals ``nMechanisms``;
        (2) every mode the diagnosis describes by a direction / an axis (rigid-body and single-node
        modes) satisfies K phi = 0 on that matrix and violates no support or constraint; (3) these
        modes are linearly independent and, together with the count of the "distributed" modes, add up
        to ``nMechanisms``; (4) every node the diagnosis lists takes part in the independent null space.
        A node that is LOADED in a degree of freedom without any stiffness (no member, or only truss
        members for a rotation) is such a mode too: these DOFs are kept in the independent matrix (the
        unloaded ones are not part of the solution, by the engine's documented rule), found there as
        zero eigenvalues, and each DOF the diagnosis lists must be loaded according to the load
        definitions.
tier K  the same vectors on the ENGINE's own matrix (K phi = 0), and the engine's eigenvalue count
        recomputed with numpy from the engine's assembled matrix.
tier C  evidence from the foreign solvers (OpenSeesPy, PyNite) under a PROBING LOAD built from the null
        vectors of the independent matrix (a load that does work on every mechanism, so that no
        equilibrium state exists): a solver that fails, reports a singular matrix or returns
        displacements beyond any physical size agrees with the diagnosis; one that returns a finite,
        equilibrated solution contradicts it.  (The model's own loading is no test: a mechanism it does
        not excite is solved without complaint -- measured with both solvers.)
tier E  the diagnosis is a property of the structure: it must be identical after every design variable
        is changed, and identical when the model is analysed again (determinism).

Not verifiable here (the record is then "quarantined", never "verified"):
* a mechanism that appears only because tension-/compression-only members go slack
  (``stability.slack-mechanism``): that these members are slack is a statement about a nonlinear
  problem without a unique solution; only "the structure without them is a mechanism" is checked;
* an "ill-conditioned" classification whose smallest eigenvalue sits at the tolerance.
"""
from __future__ import annotations

import copy
import time

import numpy as np

from ..model import RIGID
from .common import TRUSS_TYPES
from .policy import POLICY
from .record import TierResult

__all__ = ["independent_stiffness", "IndependentSystem", "verify_diagnosis_tiers", "mode_vector"]

DOFS = ("ux", "uy", "uz", "rx", "ry", "rz")
AX = {"X": 0, "Y": 1, "Z": 2}


# ====================================================================================== member matrices
def _bending_block(EI, L, GAs):
    """4x4 stiffness of a prismatic beam for (v1, t1, v2, t2) with t = dv/dx positive, shear deformation
    through phi = 12 EI / (GAs L^2) (0 without).  Closed form."""
    phi = 12.0 * EI / (GAs * L * L) if GAs else 0.0
    c = EI / (L ** 3 * (1.0 + phi))
    return c * np.array([[12.0, 6.0 * L, -12.0, 6.0 * L],
                         [6.0 * L, (4.0 + phi) * L * L, -6.0 * L, (2.0 - phi) * L * L],
                         [-12.0, -6.0 * L, 12.0, -6.0 * L],
                         [6.0 * L, (2.0 - phi) * L * L, -6.0 * L, (4.0 + phi) * L * L]])


def _local_frame_matrix(E, G, A, Iy, Iz, J, L, GAy, GAz):
    """12x12 local stiffness, DOF order (ux, uy, uz, rx, ry, rz) at i then j.
    Plane x-y: deflection uy, slope rz = +duy/dx, inertia Iz, shear area along y.
    Plane x-z: deflection uz, slope ry = -duz/dx, inertia Iy, shear area along z."""
    k = np.zeros((12, 12))
    for a, b, s in ((0, 6, E * A / L), (3, 9, G * J / L)):
        k[a, a] += s
        k[b, b] += s
        k[a, b] -= s
        k[b, a] -= s
    idx = [1, 5, 7, 11]
    k[np.ix_(idx, idx)] += _bending_block(E * Iz, L, GAy)
    idx = [2, 4, 8, 10]
    S = np.diag([1.0, -1.0, 1.0, -1.0])                     # ry = -dw/dx
    k[np.ix_(idx, idx)] += S @ _bending_block(E * Iy, L, GAz) @ S
    return k


def _condense_connections(k, springs):
    """Member with a spring of stiffness ``springs[a]`` between its end DOF a and the joint (inf = rigid,
    0 = released): stiffness seen from the joints, by static condensation of the member-side DOFs."""
    s = np.asarray(springs, float)
    soft = np.flatnonzero(np.isfinite(s))
    if not len(soft):
        return k
    rigid = np.flatnonzero(~np.isfinite(s))
    ks = np.diag(s[soft])
    # unknowns: joint side d (12) and member side q at the soft DOFs; energy = 1/2 q_all' k q_all +
    # 1/2 (q - d)' ks (q - d) with q_all = d at the rigid DOFs
    A = k[np.ix_(soft, soft)] + ks                          # member-side block
    B = np.zeros((len(soft), 12))                           # coupling member side <-> joint side
    B[:, rigid] = k[np.ix_(soft, rigid)]
    B[np.arange(len(soft)), soft] -= s[soft]
    D = np.zeros((12, 12))
    D[np.ix_(rigid, rigid)] = k[np.ix_(rigid, rigid)]
    D[soft, soft] += s[soft]
    return D - B.T @ np.linalg.pinv(A, rcond=1e-13, hermitian=True) @ B


def _skew(r):
    return np.array([[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]])


def _end_transfer(offset):
    """(u, theta) at the flexible end from (u, theta) of the joint: u_end = u + theta x offset."""
    T = np.eye(6)
    T[:3, 3:] = -_skew(np.asarray(offset, float))
    return T


class IndependentSystem:
    """Stiffness of the structure assembled independently of ``benchfea.solve``.

    K      (n, n) nodal matrix, support springs included
    T      (n, n_retained) slaves expressed by masters
    free   retained DOFs that are neither fixed nor "inactive" (the engine's documented rule for DOFs
           that are not part of the solution: free DOFs of memberless nodes that drive no slave, and
           rotations of nodes joined only by truss-type members -- docs/ARCHITECTURE.md)
    """

    def __init__(self, model, exclude=()):
        self.model = model
        n = 6 * len(model.nodes)
        self.n = n
        xyz = np.array([nd.xyz for nd in model.nodes.values()], float)
        self.xyz = xyz
        self.ids = list(model.nodes)
        idx = {nid: k for k, nid in enumerate(self.ids)}
        self.idx = idx
        self.size = float(np.linalg.norm(xyz.max(axis=0) - xyz.min(axis=0))) or 1.0
        K = np.zeros((n, n))
        shear = bool(model.analysis.shear_deformation)
        attached = {nid: [] for nid in self.ids}
        self.member_nodes = {}
        for mid, m in model.members.items():
            if mid in exclude:
                continue
            mat, sec = model.materials[m.material], model.sections[m.section]
            L = float(m.length)
            R = np.asarray(m.R, float)
            truss = m.type in TRUSS_TYPES
            if truss:
                k = np.zeros((12, 12))
                s = mat.E * sec.A / L
                k[0, 0] = k[6, 6] = s
                k[0, 6] = k[6, 0] = -s
            else:
                GAy = mat.G * sec.Asy if (shear and sec.Asy) else None
                GAz = mat.G * sec.Asz if (shear and sec.Asz) else None
                k = _local_frame_matrix(mat.E, mat.G, sec.A, sec.Iy, sec.Iz, sec.J, L, GAy, GAz)
                k = _condense_connections(k, tuple(m.k_end_i) + tuple(m.k_end_j))
            R6 = np.zeros((6, 6))
            R6[:3, :3] = R6[3:, 3:] = R
            G12 = np.zeros((12, 12))
            G12[:6, :6] = R6 @ _end_transfer(m.offset_i)
            G12[6:, 6:] = R6 @ _end_transfer(m.offset_j)
            dofs = np.r_[6 * idx[m.i] + np.arange(6), 6 * idx[m.j] + np.arange(6)]
            K[np.ix_(dofs, dofs)] += G12.T @ k @ G12
            attached[m.i].append((mid, truss and not np.any(m.offset_i)))
            attached[m.j].append((mid, truss and not np.any(m.offset_j)))
            self.member_nodes[mid] = (m.i, m.j)
        fixed = np.zeros(n, bool)
        spring = np.zeros(n)
        for nid, sup in model.supports.items():
            fixed[6 * idx[nid]:6 * idx[nid] + 6] = sup.fixed
            spring[6 * idx[nid]:6 * idx[nid] + 6] = sup.spring
        K[np.arange(n), np.arange(n)] += spring
        self.K, self.fixed, self.spring, self.attached = 0.5 * (K + K.T), fixed, spring, attached
        # ---- constraints: slave DOF = row of master DOFs
        slave = np.zeros(n, bool)
        driven = np.zeros(n, bool)
        T = np.eye(n)
        masters = set()
        for c in model.constraints.values():
            masters.add(c.master)
            m0 = 6 * idx[c.master]
            cd = [DOFS.index(d) for d in c.dofs]
            for sid in c.slaves:
                r = xyz[idx[sid]] - xyz[idx[c.master]]
                if c.type == "rigidDiaphragm":
                    nrm = np.eye(3)[AX[c.normal]]
                    r = r - (r @ nrm) * nrm
                C = _end_transfer(r)                         # u_s = u_m + theta_m x r, theta_s = theta_m
                cols = cd if c.type == "rigidDiaphragm" else range(6)
                for a in cd:
                    g = 6 * idx[sid] + a
                    T[g, :] = 0.0
                    for b in cols:
                        T[g, m0 + b] = C[a, b]
                        if C[a, b]:
                            driven[m0 + b] = True
                    slave[g] = True
        self.slave, self.masters = slave, masters
        self.retained = np.flatnonzero(~slave)
        self.T = T[:, self.retained]
        # ---- DOFs that are not part of the solution (the engine's rule, re-implemented)
        Kc = self.T.T @ self.K @ self.T
        diag0 = np.zeros(n, bool)
        diag0[self.retained] = np.abs(Kc).sum(axis=1) == 0.0
        inactive = np.zeros(n, bool)
        for nid in self.ids:
            i = 6 * idx[nid]
            ok = ~fixed[i:i + 6] & (spring[i:i + 6] == 0.0) & ~slave[i:i + 6] & diag0[i:i + 6]
            if not attached[nid]:
                inactive[i:i + 6] = ok & ~driven[i:i + 6]
            elif all(t for _, t in attached[nid]) and nid not in masters:
                ok[:3] = False
                inactive[i:i + 6] = ok
        # A LOADED degree of freedom without stiffness is not dropped: the engine refuses such a model
        # ("no-stiffness" / "only-truss-members" mode of kind "node"), and the motion must then be found
        # here as a zero eigenvalue.  Which DOFs are loaded comes from the load definitions, not from the
        # diagnosis: the first combination (model order) that loads a DOF without stiffness.
        loaded = _loaded_dofs(model, inactive)
        self.loaded_inactive = inactive & loaded
        inactive = inactive & ~loaded
        self.inactive = inactive
        keep = ~fixed[self.retained] & ~inactive[self.retained]
        self.free = np.flatnonzero(keep)
        self.Kff = Kc[np.ix_(self.free, self.free)]

    # ---------------------------------------------------------------- eigenvalues of the scaled matrix
    def scaled_eigen(self):
        d = np.diag(self.Kff).copy()
        d = np.sqrt(np.where(d > 0.0, d, 1.0))
        Ks = self.Kff / np.outer(d, d)
        lam, vec = np.linalg.eigh(0.5 * (Ks + Ks.T))
        return lam, vec, d

    def null_space(self, tol):
        """Nodal null vectors (n, m) of the reduced matrix: eigenvalues of the scaled matrix below
        ``tol``; and the eigenvalues."""
        lam, vec, d = self.scaled_eigen()
        sel = lam < tol
        ur = np.zeros((len(self.retained), int(sel.sum())))
        ur[self.free] = vec[:, sel] / d[:, None]
        return self.T @ ur, lam

    def weights(self):
        return np.tile([1.0, 1.0, 1.0, self.size, self.size, self.size], self.n // 6)

    def admissible(self, v):
        """Largest violation of supports and constraints by a nodal vector (relative to its size)."""
        w = self.weights()
        vw = w * v
        full = max(float(np.abs(vw).max()), 1e-300)
        viol = float(np.abs(vw[self.fixed]).max(initial=0.0))
        img = self.T @ v[self.retained]
        viol = max(viol, float(np.abs(w * (img - v)).max(initial=0.0)))
        return viol / full

    def residual(self, v):
        """|K v| on the free DOFs relative to |K| |v| (0 = zero-energy motion)."""
        vr = v[self.retained]
        f = self.T.T @ (self.K @ (self.T @ vr))
        num = float(np.abs(f[self.free]).max(initial=0.0))
        den = float((np.abs(self.T.T @ np.abs(self.K) @ np.abs(self.T)) @ np.abs(vr))[self.free].max(initial=0.0))
        return num / den if den > 0.0 else 0.0

    def components(self):
        """Connected parts (members and constraints), as lists of node indices."""
        parent = list(range(len(self.ids)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        for i, j in self.member_nodes.values():
            parent[find(self.idx[i])] = find(self.idx[j])
        for c in self.model.constraints.values():
            for s in c.slaves:
                parent[find(self.idx[s])] = find(self.idx[c.master])
        out = {}
        for k in range(len(self.ids)):
            out.setdefault(find(k), []).append(k)
        return list(out.values())


def _loaded_dofs(model, inactive) -> np.ndarray:
    """Nodal DOFs (n,) that carry a nodal load in the first combination (model order) that loads one of
    the DOFs in ``inactive``; all False when no combination does.  From the load definitions
    (``crosscheck.prepare``: no solver code)."""
    from ..crosscheck.prepare import prepare
    for cid in model.combinations:
        try:
            mask = (prepare(model, cid, extra_nodal_loads={}).load != 0.0).ravel()
        except Exception:                                    # a combination that cannot be prepared loads nothing here
            continue
        if (mask & inactive).any():
            return mask
    return np.zeros(len(inactive), bool)


def independent_stiffness(model, exclude=()) -> IndependentSystem:
    return IndependentSystem(model, exclude)


def loaded_dof_vectors(system: IndependentSystem, mode: dict):
    """Unit vectors of a mode of kind "node" that the diagnosis describes only by the loaded DOFs of a
    node without stiffness there (no direction, no axis): one vector per listed DOF, with the global
    index.  Returns ``[(dof name, global index, vector)]`` or ``[]`` when the mode is not of that form."""
    if mode.get("kind") != "node" or mode.get("direction") or mode.get("axis"):
        return []
    if mode.get("cause") not in ("no-stiffness", "only-truss-members"):
        return []
    nodes = mode.get("nodes", [])
    if len(nodes) != 1 or not nodes[0].get("dofs"):
        return []
    i = 6 * system.idx[nodes[0]["node"]]
    out = []
    for name in nodes[0]["dofs"]:
        if name not in DOFS:
            return []
        v = np.zeros(system.n)
        v[i + DOFS.index(name)] = 1.0
        out.append((name, i + DOFS.index(name), v))
    return out


# ====================================================================================== modes -> vectors
def mode_vector(system: IndependentSystem, mode: dict):
    """Nodal vector (n,) of a mode of the diagnosis that is described geometrically, or None:
    rigid-body translation (direction) / rotation (axis, through) of the connected part that contains
    the mode's nodes, single-node translation (direction) / rotation (axis)."""
    n = system.n
    v = np.zeros(n)
    nodes = [p["node"] for p in mode.get("nodes", [])]
    if not nodes:
        return None
    motion = mode.get("motion")
    if mode["kind"] == "node":
        i = 6 * system.idx[nodes[0]]
        if motion == "translation" and mode.get("direction"):
            v[i:i + 3] = mode["direction"]
        elif motion == "rotation" and mode.get("axis"):
            v[i + 3:i + 6] = mode["axis"]
        else:
            return None
        return v
    if mode["kind"] != "rigid-body":
        return None
    comp = next(c for c in system.components() if system.idx[nodes[0]] in c)
    if motion == "translation" and mode.get("direction"):
        for k in comp:
            v[6 * k:6 * k + 3] = mode["direction"]
    elif motion == "rotation" and mode.get("axis") and mode.get("through"):
        axis, through = np.array(mode["axis"], float), np.array(mode["through"], float)
        for k in comp:
            v[6 * k:6 * k + 3] = np.cross(axis, system.xyz[k] - through)
            v[6 * k + 3:6 * k + 6] = axis
    else:
        return None                                          # screw motion: the pitch is not in the record
    v[system.inactive] = 0.0
    return v


# ====================================================================================== tiers
def _engine_system(model, exclude=()):
    """The engine's own reduced matrix (for tier K)."""
    from ..solve.element import build_element
    from ..solve.system import Structure
    elements = {mid: build_element(model, m) for mid, m in model.members.items() if mid not in exclude}
    st = Structure(model, elements)
    K = st.stiffness()
    Kr = st.reduce(K)
    return st, K, Kr


def tier_null_space(model, diagnosis: dict, exclude=()) -> TierResult:
    t0 = time.perf_counter()
    pol = POLICY["diagnosis"]
    tol_engine = float(model.analysis.pivot_tolerance)
    sysm = independent_stiffness(model, exclude)
    lam, _, _ = sysm.scaled_eigen()
    true_mechanism = diagnosis.get("classification") == "mechanism"
    thr = pol["singular"] if true_mechanism else tol_engine
    count = int((lam < thr).sum())
    failed, largest = [], []
    n_mech = int(diagnosis["nMechanisms"])
    at_threshold = [float(x) for x in lam if 0.5 * thr <= x <= 2.0 * thr] if not true_mechanism else []
    if at_threshold:
        failed.append(f"an eigenvalue of the independent matrix ({at_threshold[0]:.2e}) lies at the tolerance "
                      f"{thr:.1e}: the count of ill-conditioned modes is not determined")
    if count != n_mech:
        failed.append(f"the independent matrix has {count} eigenvalue(s) below {thr:.1e}; the diagnosis says "
                      f"{n_mech} mechanism(s)")
    N, _ = sysm.null_space(thr)
    w = sysm.weights()
    Q = np.linalg.qr(w[:, None] * N)[0] if N.shape[1] else np.zeros((sysm.n, 0))
    vectors, rows, distributed = [], [], 0
    for k, mode in enumerate(diagnosis["modes"]):
        row = {"mode": k, "kind": mode["kind"], "description": mode["description"]}
        unit = loaded_dof_vectors(sysm, mode)
        if mode["kind"] == "distributed":
            distributed += int(mode.get("count", 1))
            row["check"] = "count and participating nodes only (no single vector is defined)"
        elif unit:
            # a node loaded where it has no stiffness: every listed DOF must be, by the independent rule,
            # a loaded DOF without stiffness, and its unit motion a zero-energy motion
            row["check"] = "each loaded DOF: no stiffness in the independent matrix, loaded by the load definitions"
            worst = 0.0
            for name, g, v in unit:
                if not sysm.loaded_inactive[g]:
                    failed.append(f"mode {k}: by the independent rule {mode['nodes'][0]['node']} {name} is not a "
                                  "loaded degree of freedom without stiffness")
                    continue
                worst = max(worst, sysm.residual(v), sysm.admissible(v))
                vectors.append(w * v / max(np.linalg.norm(w * v), 1e-300))
            row["residual"] = worst
            largest.append({"what": f"mode {k}: |K phi| / (|K| |phi|) of the loaded DOFs on the independent matrix",
                            "value": worst, "tolerance": pol["nullResidual"]})
            if worst > pol["nullResidual"]:
                failed.append(f"mode {k}: a loaded DOF has stiffness or is restrained in the independent matrix "
                              f"({worst:.2e})")
        else:
            v = mode_vector(sysm, mode)
            if v is None:
                row["check"] = "not reconstructable from the record (combined motion): counted, not tested"
                distributed += 1
            else:
                res, adm = sysm.residual(v), sysm.admissible(v)
                row.update(check="K phi = 0", residual=res, constraintViolation=adm)
                # an "ill-conditioned" mode has a small but nonzero stiffness: its residual is bounded by
                # the engine's tolerance, not by round-off
                lim = pol["nullResidual"] if true_mechanism else max(pol["nullResidual"], 10.0 * tol_engine)
                largest.append({"what": f"mode {k}: |K phi| / (|K| |phi|) on the independent matrix",
                                "value": res, "tolerance": lim})
                if res > lim:
                    failed.append(f"mode {k} ({mode['description'][:60]}): K phi = {res:.2e} on the independent matrix")
                if adm > pol["nullResidual"]:
                    failed.append(f"mode {k}: the described motion violates a support or constraint ({adm:.2e})")
                vectors.append(w * v / max(np.linalg.norm(w * v), 1e-300))
        # every listed node takes part in the independent null space
        missing = []
        for p in mode.get("nodes", []):
            i = 6 * sysm.idx[p["node"]]
            if Q.shape[1] == 0 or float(np.abs(Q[i:i + 6]).max()) <= 1e-7:
                missing.append(p["node"])
        if missing:
            row["nodesNotInNullSpace"] = missing
            failed.append(f"mode {k}: node(s) {', '.join(missing[:5])} do not move in the independent null space")
        rows.append(row)
    rank = int(np.linalg.matrix_rank(np.column_stack(vectors), tol=1e-6)) if vectors else 0
    if rank != len(vectors):
        failed.append("the modes described by the diagnosis are not linearly independent")
    if rank + distributed != n_mech:
        failed.append(f"the modes of the diagnosis add up to {rank + distributed}, nMechanisms is {n_mech}")
    passed = not failed
    summary = (f"independent stiffness: {count} zero-stiffness motion(s) (diagnosis: {n_mech}); "
               f"{len(vectors)} described mode(s) tested for K phi = 0: "
               + ("agree" if passed else "NOT confirmed: " + "; ".join(failed[:3])))
    return TierResult("N", True, passed, summary, largest[:6],
                      {"singular": thr, "nullResidual": pol["nullResidual"]},
                      {"eigenvaluesBelowThreshold": [float(x) for x in lam[lam < thr][:20]],
                       "smallestEigenvalueAbove": float(lam[lam >= thr][0]) if (lam >= thr).any() else None,
                       "modes": rows, "excludedMembers": sorted(exclude), "failed": failed},
                      time.perf_counter() - t0)


def tier_engine_matrix(model, diagnosis: dict, exclude=()) -> TierResult:
    t0 = time.perf_counter()
    pol = POLICY["diagnosis"]
    st, K, Kr = _engine_system(model, exclude)
    f = st.free_r
    Kff = Kr[np.ix_(f, f)]
    d = np.diag(Kff).copy()
    d = np.sqrt(np.where(d > 0.0, d, 1.0))
    lam = np.linalg.eigvalsh(Kff / np.outer(d, d))
    true_mechanism = diagnosis.get("classification") == "mechanism"
    thr = pol["singular"] if true_mechanism else float(model.analysis.pivot_tolerance)
    count = int((lam < thr).sum())
    failed, largest = [], []
    sysm = independent_stiffness(model, exclude)             # geometry of the mode vectors, loaded DOFs
    # the engine leaves DOFs without stiffness out of its reduced matrix; the loaded ones among them are the
    # mechanisms of a "no-stiffness" diagnosis: counted from the engine's own list of such DOFs and checked
    # on its unreduced matrix (the row and column must be zero)
    loaded_out = np.asarray(st.inactive, bool) & sysm.loaded_inactive
    n_loaded = int(loaded_out.sum())
    kmax = float(np.abs(K).max()) or 1.0
    for g in np.flatnonzero(loaded_out):
        res = float(np.abs(K[:, g]).max()) / kmax
        largest.append({"what": f"loaded DOF {sysm.ids[g // 6]} {DOFS[g % 6]}: column of the engine's matrix / |K|",
                        "value": res, "tolerance": pol["nullResidual"]})
        if res > pol["nullResidual"]:
            failed.append(f"{sysm.ids[g // 6]} {DOFS[g % 6]} has stiffness in the engine's matrix")
    if count + n_loaded != int(diagnosis["nMechanisms"]):
        failed.append(f"the engine's matrix has {count} eigenvalue(s) below {thr:.1e}"
                      + (f" and {n_loaded} loaded DOF(s) without stiffness" if n_loaded else "")
                      + f"; the diagnosis says {diagnosis['nMechanisms']}")
    absK = np.abs(Kr)
    for k, mode in enumerate(diagnosis["modes"]):
        v = mode_vector(sysm, mode)
        if v is None:
            continue
        v = v.copy()
        v[st.inactive] = 0.0
        vr = v[st.retained]
        num = float(np.abs((Kr @ vr)[f]).max(initial=0.0))
        den = float((absK @ np.abs(vr))[f].max(initial=0.0))
        res = num / den if den > 0.0 else 0.0
        lim = pol["nullResidual"] if true_mechanism else max(pol["nullResidual"], 10.0 * thr)
        largest.append({"what": f"mode {k}: |K phi| / (|K| |phi|) on the engine's matrix", "value": res,
                        "tolerance": lim})
        if res > lim:
            failed.append(f"mode {k}: K phi = {res:.2e} on the engine's matrix")
    passed = not failed
    return TierResult("K", True, passed,
                      f"engine matrix: {count} eigenvalue(s) below {thr:.1e}"
                      + (f", {n_loaded} loaded DOF(s) without stiffness" if n_loaded else "") + "; "
                      + ("the described modes satisfy K phi = 0" if passed else "; ".join(failed[:3])),
                      largest[:6], {"singular": thr, "nullResidual": pol["nullResidual"]},
                      {"eigenvaluesBelowThreshold": [float(x) for x in lam[lam < thr][:20]], "failed": failed},
                      time.perf_counter() - t0)


def probe_document(raw: dict, system: IndependentSystem, threshold: float):
    """The model with ONE load case that does work on every zero-stiffness motion: nodal loads
    proportional to sum_k (k + 1) n_k over an orthonormal basis n_k of the independent null space.  Under
    such a load no equilibrium state exists, so a solver that returns a finite, equilibrated solution
    contradicts the diagnosis.  Returns ``(document, combination id)`` or ``(None, None)`` when the
    independent matrix has no null space."""
    N, _ = system.null_space(threshold)
    if N.shape[1] == 0:
        return None, None
    w = system.weights()
    Q = np.linalg.qr(w[:, None] * N)[0]
    F = (Q @ np.arange(1.0, Q.shape[1] + 1.0)) / w           # work F . n_k = (k + 1) != 0 for every k
    F[system.fixed | system.inactive] = 0.0
    F = F / max(float(np.abs(F).max()), 1e-300)
    loads = []
    for k, nid in enumerate(system.ids):
        f = F[6 * k:6 * k + 6]
        if np.abs(f).max() > 1e-12:
            loads.append({"type": "nodal", "node": nid,
                          **{n: float(v) for n, v in zip(("Fx", "Fy", "Fz", "Mx", "My", "Mz"), f) if abs(v) > 1e-12}})
    doc = copy.deepcopy(raw)
    doc.pop("design", None)
    an = doc.get("analysis")
    if an:
        an.pop("type", None)
        an.pop("buckling", None)
        an.pop("stiffnessReduction", None)
    doc["loadCases"] = {"verify.probe": {"loads": loads}}
    doc["combinations"] = {"verify.probe": {"factors": {"verify.probe": 1.0}}}
    return doc, "verify.probe"


def tier_foreign_evidence(raw, model, diagnosis: dict, *, cache=None, timeout=120.0, deadline=None,
                          exclude=()) -> TierResult:
    """Foreign solvers on the unstable model under a PROBING LOAD that does work on every zero-stiffness
    motion (``probe_document``).  The original loading is no test: a mechanism that the loads do not
    excite is solved without complaint by a solver that does not test for singularity (measured with
    both).  Under the probing load a solver agrees with the diagnosis when it fails, reports a singular
    system or returns displacements beyond ``saneDisplacement`` x model size; it contradicts a
    "mechanism" classification when it returns a finite solution in global equilibrium."""
    from ..crosscheck.compare import equilibrium
    from ..model import build_model
    from ..problems import ModelInvalid
    from .common import document_hash
    from .tier_c import DEFAULT_CACHE, SOLVERS, _foreign, foreign_versions
    t0 = time.perf_counter()
    cache = cache if cache is not None else DEFAULT_CACHE
    versions = foreign_versions()
    have = [s for s in SOLVERS if versions.get(s)]
    true_mechanism = diagnosis.get("classification") == "mechanism"

    def skipped(why, category="not-run", applies=True):
        t = TierResult("C", applies, None, why, [], {}, {}, time.perf_counter() - t0, skipped=applies)
        t.category = category
        return t
    if exclude:
        return skipped("not applicable: the mechanism exists only with one-sided members inactive, which the "
                       "foreign solvers do not model", applies=False)
    if not true_mechanism:
        return skipped("not applicable: an ill-conditioned (not singular) model can be solved by any solver",
                       applies=False)
    if not have:
        return skipped("not run: no foreign solver is installed", "not-representable")
    if deadline is not None and time.perf_counter() > deadline:
        return skipped("not run (time budget)")
    sysm = independent_stiffness(model)
    doc, cid = probe_document(raw, sysm, POLICY["diagnosis"]["singular"])
    if doc is None:
        t = TierResult("C", True, False, "no probing load: the independent matrix has no null space (tier N)")
        t.category = "disagreement"
        return t
    try:
        pmodel = build_model(copy.deepcopy(doc))
    except ModelInvalid as e:
        t = TierResult("C", True, False, "the probing model is not valid: " + "; ".join(p.message for p in e.problems[:2]))
        t.category = "error"
        return t
    size = sysm.size
    sane = POLICY["diagnosis"]["saneDisplacement"]
    h = document_hash(doc)
    evidence, contradicts = {}, []
    for s in have:
        res, _ = _foreign(h, pmodel, cid, s, {}, cache, timeout, versions)
        if not res.ok or res.displacements is None:
            verdict, note = "consistent", f"no solution: {(res.error or 'failed')[:200]}"
        else:
            umax = float(np.abs(res.displacements[:, :3]).max())
            rmax = float(np.abs(res.displacements[:, 3:]).max()) * size
            big = max(umax, rmax)
            if not np.isfinite(big) or big > sane * size:
                verdict, note = "consistent", f"displacements of {big:.3e} mm (model size {size:.3e} mm)"
            elif res.fidelity.level == "unsupported":
                verdict, note = "no evidence", "the solver cannot represent the model"
            else:
                eq = equilibrium(pmodel, res)
                r = max(eq.get("force") or 0.0, eq.get("moment") or 0.0)
                if r > 1e-6:
                    verdict, note = "consistent", (f"the returned numbers are not in global equilibrium with the "
                                                   f"probing load (residual {r:.2e})")
                else:
                    verdict, note = "contradicts", (f"returned an equilibrated solution with displacements up to "
                                                    f"{big:.3e} mm under the probing load")
        evidence[s] = {"verdict": verdict, "note": note}
        if verdict == "contradicts":
            contradicts.append(f"{s}: {note}")
    n_cons = sum(1 for e in evidence.values() if e["verdict"] == "consistent")
    need = POLICY["C"]["requiredSolvers"]
    passed = not contradicts and n_cons >= need
    if passed:
        summary = (f"{n_cons} foreign solver(s) find no equilibrium state under a load that does work on the "
                   "mechanism(s)")
    elif contradicts:
        summary = "a foreign solver found an equilibrium state under the probing load: " + "; ".join(contradicts[:2])
    else:
        summary = f"only {n_cons} foreign solver(s) gave evidence; {need} are required"
    t = TierResult("C", True, passed, summary, [], {"saneDisplacement": f"{sane:g} x model size"},
                   {"evidence": evidence, "foreignVersions": versions,
                    "probe": "nodal loads proportional to a combination of the null vectors of the independent "
                             "stiffness with nonzero work on each"}, time.perf_counter() - t0)
    t.category = "" if passed else ("disagreement" if contradicts else "not-representable")
    return t


def tier_invariance(raw, diagnosis: dict, analyse) -> TierResult:
    """The diagnosis after every design variable is changed, and on a second run."""
    from .tier_e import perturb_design_document
    t0 = time.perf_counter()
    failed = []
    for what, doc in (("second run", copy.deepcopy(raw)), ("design variables changed", perturb_design_document(raw))):
        d2 = analyse(doc)
        if d2 != diagnosis:
            failed.append(f"{what}: the diagnosis differs")
    passed = not failed
    return TierResult("E", True, passed, "diagnosis identical on a second run and after changing every design "
                      "variable" if passed else "; ".join(failed), [], {"designVariables": 0.0}, {"failed": failed},
                      time.perf_counter() - t0)


def verify_diagnosis_tiers(raw, model, diagnosis: dict, analyse, *, exclude=(), cache=None, timeout=120.0,
                           deadline=None, run=("N", "K", "C", "E")) -> dict:
    tiers = {}
    for key, fn in (("N", lambda: tier_null_space(model, diagnosis, exclude)),
                    ("K", lambda: tier_engine_matrix(model, diagnosis, exclude)),
                    ("C", lambda: tier_foreign_evidence(raw, model, diagnosis, cache=cache, timeout=timeout,
                                                        deadline=deadline, exclude=exclude)),
                    ("E", lambda: tier_invariance(raw, diagnosis, analyse))):
        if key in run:
            tiers[key] = fn()
        else:
            tiers[key] = TierResult(key, True, None, "not run (excluded by the caller)", skipped=True)
    return tiers
