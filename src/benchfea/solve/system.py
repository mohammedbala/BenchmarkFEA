"""Assembly, exact kinematic constraints, supports and the linear solve.

DOF spaces
----------
nodal      all ``6 * n_nodes`` DOFs, ``Model.dof`` numbering.
retained   nodal DOFs that are not slaves of a constraint; ``u = T u_r`` (master-slave transformation,
           THEORY.md section 4).  T has an identity row for every retained DOF.
free       retained DOFs that are neither fixed by a support nor inactive.

``K_r = T^T K T`` and ``F_r = T^T F`` are partitioned into free (f) and prescribed (p) DOFs:
``K_ff u_f = F_f - K_fp u_p``.  Nothing is penalised; fixed DOFs are eliminated and slaves are
eliminated.

The pieces are separate on purpose so that a nonlinear driver can reuse them: ``Structure.assemble``
adds per-member 12x12 global matrices (elastic by default, any other matrix set -- e.g. geometric
stiffness -- by passing ``matrices``), ``Structure.factorize`` reduces + factorises (and diagnoses a
singular matrix), ``Structure.solve`` back-substitutes one load vector.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import cho_solve, cholesky, eigh, eigvalsh, solve_triangular

from ..model import DOFS
from .element import MemberElement, rigid_transfer
from .errors import MechanismDiagnosis, MechanismMode, NodeParticipation, unstable

__all__ = ["Structure", "Factor", "Solution", "EXACT_EIGEN_LIMIT"]

EXACT_EIGEN_LIMIT = 1200       # free DOFs up to which the smallest eigenvalue is computed exactly (dense)
AXES = ("X", "Y", "Z")


@dataclass
class Factor:
    K: np.ndarray                 # (n, n) nodal stiffness that was factorised (members + support springs)
    Kr: np.ndarray                # (nr, nr) retained
    free: np.ndarray              # indices into the retained numbering
    pres: np.ndarray
    scale: np.ndarray             # sqrt(diag(K_ff))
    chol: np.ndarray              # lower Cholesky factor of the diagonally scaled K_ff
    lambda_min: float             # smallest / largest eigenvalue of the scaled K_ff
    lambda_max: float
    min_pivot: float              # smallest Cholesky pivot of the scaled K_ff (its diagonal is 1)
    exact: bool                   # eigenvalues exact (dense) or estimated (inverse / power iteration)

    @property
    def condition(self) -> float:
        return self.lambda_max / self.lambda_min

    def solve_free(self, rhs: np.ndarray) -> np.ndarray:
        if not len(self.free):
            return np.zeros(0)
        return cho_solve((self.chol, True), rhs / self.scale) / self.scale


@dataclass
class Solution:
    u: np.ndarray                 # (n,) nodal displacements
    reactions: np.ndarray         # (n,) forces by the supports on the structure (fixed DOFs and springs)
    constraint_forces: np.ndarray  # (n,) forces by the rigid constraints on the nodes
    F: np.ndarray                 # (n,) nodal load vector used (nodal loads + equivalent member loads)
    residual: float = 0.0         # max |K_ff u_f - rhs| / max (|K_ff| |u_f| + |rhs|) of the linear solve
    extra: dict = field(default_factory=dict)


class Structure:
    def __init__(self, model, elements: dict[str, MemberElement]):
        self.model = model
        self.elements = elements
        n = self.n = model.n_dof
        self.node_ids = list(model.nodes)
        self.xyz = np.array([nd.xyz for nd in model.nodes.values()], float).reshape(-1, 3)
        self.size = float(np.linalg.norm(np.ptp(self.xyz, axis=0))) if len(self.xyz) else 0.0
        self.fixed = np.zeros(n, bool)
        self.spring = np.zeros(n)
        for nid, s in model.supports.items():
            i = 6 * model.nodes[nid].index
            self.fixed[i:i + 6] = s.fixed
            self.spring[i:i + 6] = s.spring
        self.warnings: list[str] = []
        self._constraints()
        self._inactive()
        r = self.retained
        self.fixed_r = np.flatnonzero(self.fixed[r])
        self.free_r = np.flatnonzero(~self.fixed[r] & ~self.inactive[r])

    # ---------------------------------------------------------------- constraints
    def _constraints(self):
        """u = T u_r.  Rigid link: u_s = u_m + theta_m x r, theta_s = theta_m (r = x_s - x_m).
        Diaphragm with normal n: the same relation for the two in-plane translations and the rotation
        about n only, with the in-plane part of r, driven by the master's in-plane translations and its
        rotation about n."""
        model, n = self.model, self.n
        slave = np.zeros(n, bool)
        rows = {}
        involved = np.zeros(n, bool)
        self.master_of: dict[str, str] = {}
        for c in model.constraints.values():
            m = model.nodes[c.master]
            cd = [DOFS.index(d) for d in c.dofs]
            for sid in c.slaves:
                s = model.nodes[sid]
                r = s.xyz - m.xyz
                if c.type == "rigidDiaphragm":
                    nrm = np.eye(3)[AXES.index(c.normal)]
                    r = r - (r @ nrm) * nrm
                C = rigid_transfer(r)
                self.master_of[sid] = c.master
                for a in cd:
                    row = np.zeros(6)
                    cols = cd if c.type == "rigidDiaphragm" else range(6)
                    for b in cols:
                        row[b] = C[a, b]
                    slave[6 * s.index + a] = True
                    rows[6 * s.index + a] = (6 * m.index, row)
                    involved[6 * s.index + a] = True
                    involved[6 * m.index + np.flatnonzero(row)] = True
        self.slave = slave
        self.constraint_dofs = involved
        self.retained = np.flatnonzero(~slave)
        pos = -np.ones(n, int)
        pos[self.retained] = np.arange(len(self.retained))
        T = np.zeros((n, len(self.retained)))
        T[self.retained, np.arange(len(self.retained))] = 1.0
        for g, (m0, row) in rows.items():
            T[g, pos[m0:m0 + 6]] = row
        self.T = T
        self.pos = pos

    def _inactive(self):
        """Free DOFs that no member, spring or constraint gives any stiffness and that are meaningless
        rather than a mechanism of the structure:
        (a) every free DOF of a node without members, unless it drives the slaves of a constraint;
        (b) the rotations of a node to which only truss-type members (without rigid offsets) connect.
        They are taken out of the system and reported as zero (with a warning).  A load on one of them
        is refused."""
        model, n = self.model, self.n
        attached: dict[str, list] = {nid: [] for nid in model.nodes}
        for e in self.elements.values():
            for nid in e.nodes:
                attached[nid].append(e)
        self.attached = attached
        diag_zero = np.ones(n, bool)
        for e in self.elements.values():
            diag_zero[e.dofs[np.abs(e.k_global).sum(axis=1) > 0.0]] = False
        inactive = np.zeros(n, bool)
        truss_nodes, loose_nodes = [], []
        for nid, nd in model.nodes.items():
            i = 6 * nd.index
            ok = ~self.fixed[i:i + 6] & (self.spring[i:i + 6] == 0.0) & ~self.slave[i:i + 6] & diag_zero[i:i + 6]
            if not attached[nid]:
                ok &= ~self.constraint_dofs[i:i + 6]         # a master DOF that drives slaves stays
                inactive[i:i + 6] = ok
                if ok.any():
                    loose_nodes.append(nid)
            elif attached[nid] and all(e.is_truss for e in attached[nid]) and nid not in self._masters():
                ok[:3] = False
                inactive[i:i + 6] = ok
                if ok.any():
                    truss_nodes.append(nid)
        self.inactive = inactive
        self.inactive_reason = {nid: "only-truss-members" for nid in truss_nodes}
        self.inactive_reason.update({nid: "no-stiffness" for nid in loose_nodes})
        if truss_nodes:
            self.warnings.append("rotations with no stiffness (only truss-type members connect) are not part of "
                                 "the solution and are reported as 0 at node(s): " + ", ".join(truss_nodes))
        if loose_nodes:
            self.warnings.append("node(s) not connected to any member are not part of the solution; their "
                                 "displacements are reported as 0: " + ", ".join(loose_nodes))

    def _masters(self):
        return {c.master for c in self.model.constraints.values()}

    # ------------------------------------------------------------------- assembly
    def assemble(self, matrices: dict[str, np.ndarray] | None = None) -> np.ndarray:
        """Sum of the per-member 12x12 GLOBAL matrices into an (n, n) nodal matrix.  Default: the
        elastic stiffness of every element.  Pass ``matrices`` (member id -> 12x12 in the element's
        nodal global DOFs, e.g. ``G.T @ kg_local @ G``) to assemble anything else."""
        K = np.zeros((self.n, self.n))
        for mid, e in self.elements.items():
            k = e.k_global if matrices is None else matrices.get(mid)
            if k is not None:
                K[np.ix_(e.dofs, e.dofs)] += k
        return K

    def stiffness(self) -> np.ndarray:
        """Elastic nodal stiffness: members + support springs."""
        K = self.assemble()
        K[np.diag_indices(self.n)] += self.spring
        return K

    def reduce(self, K: np.ndarray) -> np.ndarray:
        Kr = self.T.T @ K @ self.T
        return 0.5 * (Kr + Kr.T)

    # ---------------------------------------------------------------------- solve
    def factorize(self, K: np.ndarray) -> Factor:
        """Reduce and factorise ``K`` (nodal, springs included).  Raises ``ModelUnstable`` with a
        :class:`MechanismDiagnosis` when the free-free matrix is not positive definite to within
        ``analysis.pivot_tolerance`` (smallest eigenvalue of the matrix scaled to unit diagonal)."""
        tol = float(self.model.analysis.pivot_tolerance)
        Kr = self.reduce(K)
        f = self.free_r
        Kff = Kr[np.ix_(f, f)]
        nf = len(f)
        if nf == 0:
            return Factor(K, Kr, f, self.fixed_r, np.ones(0), np.zeros((0, 0)), 1.0, 1.0, 1.0, True)
        diag = np.diag(Kff).copy()
        bad = ~(diag > 0.0)
        d = np.sqrt(np.where(bad, 1.0, diag))
        Ks = Kff / np.outer(d, d)
        chol = None
        if not bad.any():
            try:
                chol = cholesky(Ks, lower=True)
            except np.linalg.LinAlgError:
                chol = None
        lam_min = lam_max = None
        exact = True
        if chol is not None:
            pivot = float(np.min(np.diag(chol)) ** 2)
            if nf <= EXACT_EIGEN_LIMIT:
                lam = eigvalsh(Ks)
                lam_min, lam_max = float(lam[0]), float(lam[-1])
            else:
                exact = False
                lam_min, lam_max = _extreme_eigenvalues(Ks, chol)
            if lam_min >= tol and pivot >= tol:
                return Factor(K, Kr, f, self.fixed_r, d, chol, lam_min, lam_max, pivot, exact)
        from .stability import diagnose
        lam, vec = eigh(Ks)
        raise unstable(diagnose(self, lam, vec, d, tol), self.model.source)

    def solve(self, factor: Factor, F: np.ndarray, prescribed: dict[int, float] | None = None) -> Solution:
        """Solve for the nodal load vector ``F`` (n,) and imposed displacements at fixed DOFs."""
        r, f, p = self.retained, factor.free, factor.pres
        Fr = self.T.T @ F
        ur = np.zeros(len(r))
        for g, v in (prescribed or {}).items():
            ur[self.pos[g]] = v
        rhs = Fr[f] - factor.Kr[np.ix_(f, p)] @ ur[p]
        Kff = factor.Kr[np.ix_(f, f)]
        uf = factor.solve_free(rhs)
        for _ in range(2):                                   # iterative refinement (fixed count: deterministic)
            uf = uf + factor.solve_free(rhs - Kff @ uf)
        ur[f] = uf
        res = Kff @ uf - rhs
        # normwise backward error: residual relative to |K| |u| + |rhs|
        scale = max(float(np.max(np.abs(Kff) @ np.abs(uf) + np.abs(rhs))) if len(rhs) else 0.0, np.finfo(float).tiny)
        u = self.T @ ur
        full = factor.K @ u - F                              # force each DOF still needs from outside
        Rp = factor.Kr[p] @ ur - Fr[p]                       # at the fixed retained DOFs (masters: slaves included)
        fixed_part = np.zeros(self.n)
        fixed_part[r[p]] = Rp
        lam = np.where(self.constraint_dofs, full - fixed_part, 0.0)
        reactions = fixed_part - self.spring * u
        return Solution(u, reactions, lam, F.copy(), float(np.max(np.abs(res))) / scale if len(res) else 0.0)

    # -------------------------------------------------------------- inactive DOFs
    def check_inactive_loads(self, nodal: np.ndarray):
        """A load on a DOF without any stiffness cannot be carried: refuse with a diagnosis."""
        hit = np.flatnonzero(self.inactive & (nodal != 0.0))
        if not len(hit):
            return
        modes = []
        for nidx in sorted({int(g // 6) for g in hit}):
            nid = self.node_ids[nidx]
            comps = [int(g % 6) for g in hit if g // 6 == nidx]
            cause = self.inactive_reason.get(nid, "no-stiffness")
            why = ("only truss-type members connect to it" if cause == "only-truss-members"
                   else "no member connects to it")
            names = ", ".join(DOFS[c] for c in comps)
            share = {DOFS[c]: 1.0 / len(comps) for c in comps}
            modes.append(MechanismMode(
                "node", f"node {nid} is loaded in {names}, where it has no stiffness: {why}",
                "rotation" if all(c >= 3 for c in comps) else ("translation" if all(c < 3 for c in comps)
                                                               else "translation+rotation"),
                cause, (NodeParticipation(nid, 1.0, share),), tuple(e.id for e in self.attached[nid])))
        raise unstable(MechanismDiagnosis(len(hit), tuple(modes), float(self.model.analysis.pivot_tolerance),
                                          len(self.free_r)), self.model.source)


def _extreme_eigenvalues(Ks: np.ndarray, chol: np.ndarray, iterations: int = 60):
    """Estimates of the smallest and largest eigenvalue of an SPD matrix by inverse and direct power
    iteration (deterministic start vector)."""
    n = len(Ks)
    v = np.cos(1.0 + 0.7 * np.arange(n))            # deterministic, no special symmetry
    v /= np.linalg.norm(v)
    w = v.copy()
    lam_min = lam_max = 1.0
    for _ in range(iterations):
        y = solve_triangular(chol, v, lower=True)
        y = solve_triangular(chol.T, y, lower=False)
        lam_min = 1.0 / float(v @ y)
        v = y / np.linalg.norm(y)
        z = Ks @ w
        lam_max = float(w @ z)
        w = z / np.linalg.norm(z)
    lam_min = float(v @ (Ks @ v))
    return lam_min, lam_max
