"""The member element: 12-DOF frame member with optional shear deformation, exact full/partial end
releases and rigid end offsets; truss-type members as the axial-only special case.

Formulation (derived in docs/THEORY.md sections 1-3; nothing is taken from a reference):

* The flexible length is treated as a cantilever clamped at end i.  Its end-j flexibility ``Fb`` (6x6)
  follows from the principle of complementary virtual work, ``Fb = int B(x)^T C B(x) dx`` with the
  section resultants ``s(x) = B(x) Q_j`` given by statics and ``C = diag(1/EA, 1/GAsy, 1/GAsz, 1/GJ,
  1/EIy, 1/EIz)``.
* End connections are springs in series with the beam.  A stiff spring (k >= the beam's own
  stiffness in that end DOF; k = inf is the rigid connection) of flexibility c = 1/k adds
  ``c s s^T`` to the cantilever flexibility (s = the column that carries the spring deformation to
  end j).  A soft spring (k below that; k = 0 is the full release) is kept as a spring between the
  node-side DOF and the beam-side DOF, and the beam-side DOFs are eliminated by static condensation.
  Both are exact; choosing by magnitude means neither 1/k nor k is ever subtracted from something
  much smaller, for any k in [0, inf].  The condensed block ``K_bb + diag(k)`` is singular exactly
  when the released member is a mechanism on its own; that is refused.
* The 12x12 stiffness before condensation is ``Gamma^T K0 Gamma`` with ``Gamma = [-H, I]``
  (H = rigid-body transfer from end i to end j), which is equilibrium + contragredience.
* Fixed-end forces: the same cantilever carries the span loads (integrated exactly by
  :mod:`benchfea.solve.beam`); the end-j displacement it produces is removed by ``K0``, then the
  soft beam-side DOFs are released by the same condensation.
* Rigid offsets: ``u_end = u_node + theta_node x offset`` (exact small-displacement kinematics).

Element DOF order: (ux, uy, uz, rx, ry, rz) at end i then end j, member local axes.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigh

from ..model import DOFS, RIGID
from ..problems import pointer
from .beam import SectionFlex, SpanLoads, integrate
from .errors import refuse

__all__ = ["MemberElement", "MemberRecovery", "build_element", "StiffnessFactors", "cantilever_flexibility",
           "rigid_transfer", "skew", "release_stiffness", "TRUSS_TYPES", "MECHANISM_TOL"]

TRUSS_TYPES = ("truss", "tensionOnly", "compressionOnly")
MECHANISM_TOL = 1.0e-12       # relative eigenvalue below which the end releases are a member-level mechanism
# resultants (N, Vy, Vz, T, My, Mz) and end forces / displacement states (x, y, z, rx, ry, rz) share
# one index order for the work-conjugate pairs:
#   N-u (0), Vy-v (1), Vz-w (2), T-theta_x (3), My-theta_y (4), Mz-theta_z (5)


@dataclass(frozen=True)
class StiffnessFactors:
    """Multipliers on EA, EI (both axes), GJ and G*As of one member (1.0 = unreduced).
    Hook for the stiffness reduction of AISC 360-22 C2.3 (applied by the second-order module)."""
    axial: float = 1.0
    flexural: float = 1.0
    torsional: float = 1.0
    shear: float = 1.0


def skew(v) -> np.ndarray:
    """[v]x with [v]x a = v x a."""
    x, y, z = (float(c) for c in v)
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def rigid_transfer(r) -> np.ndarray:
    """6x6 matrix H with (u_b, theta_b) = H (u_a, theta_a) for a rigid arm r from a to b:
    u_b = u_a + theta_a x r, theta_b = theta_a.  Its transpose carries forces from b to a:
    (F, M)_a = H^T (F, M)_b, i.e. M_a = M_b + r x F."""
    H = np.eye(6)
    H[:3, 3:] = -skew(r)
    return H


def cantilever_flexibility(L: float, c: SectionFlex) -> np.ndarray:
    """End-j flexibility of the cantilever clamped at i: d_j = Fb Q_j, DOF order (x, y, z, rx, ry, rz).

    From complementary virtual work with N = Fx, Vy = Fy, Vz = Fz, T = Mx,
    My(x) = My_j - (L - x) Fz, Mz(x) = Mz_j + (L - x) Fy   (THEORY.md section 1.3)."""
    F = np.zeros((6, 6))
    F[0, 0] = c.axial * L
    F[3, 3] = c.torsion * L
    F[1, 1] = c.bend_z * L ** 3 / 3.0 + c.shear_y * L
    F[1, 5] = F[5, 1] = c.bend_z * L ** 2 / 2.0
    F[5, 5] = c.bend_z * L
    F[2, 2] = c.bend_y * L ** 3 / 3.0 + c.shear_z * L
    F[2, 4] = F[4, 2] = -c.bend_y * L ** 2 / 2.0
    F[4, 4] = c.bend_y * L
    return F


def _spd_inverse(M: np.ndarray) -> np.ndarray:
    """Inverse of a symmetric positive definite matrix with mixed units (diagonal scaling + Cholesky)."""
    d = np.sqrt(np.diag(M))
    inv = cho_solve(cho_factor(M / np.outer(d, d), lower=True), np.eye(len(M)))
    inv = inv / np.outer(d, d)
    return 0.5 * (inv + inv.T)


def release_stiffness(Fb: np.ndarray, H: np.ndarray, k_i, k_j):
    """12x12 local stiffness of the member with end springs, and the pieces needed for loads/recovery.

    Stiff connections (k >= the beam's own diagonal stiffness for that end DOF; k = inf included) are
    added in series as flexibilities to the cantilever flexibility; soft connections (k = 0 included)
    are springs between the node-side DOF and the beam-side DOF, and the beam-side DOFs are condensed
    out.  Both steps are exact; the split only decides which of the two equivalent algebraic forms is
    used, so that no large numbers are subtracted for any k in [0, inf].

    Raises ``ValueError(names)`` with the released DOFs that form a member-level mechanism.
    """
    Gamma = np.hstack([-H, np.eye(6)])
    Kbeam12 = Gamma.T @ _spd_inverse(Fb) @ Gamma
    kk = np.array(list(k_i) + list(k_j), float)
    own = np.diag(Kbeam12)
    stiff = [a for a in range(12) if kk[a] != RIGID and kk[a] >= own[a]]
    soft = [a for a in range(12) if kk[a] != RIGID and kk[a] < own[a]]
    F0 = Fb.copy()
    for a in stiff:
        col = H[:, a] if a < 6 else np.eye(6)[a - 6]         # carries the spring deformation to end j
        F0 += np.outer(col, col) / kk[a]
    K0 = _spd_inverse(F0)                                    # basic stiffness of beam + stiff springs
    K0_12 = Gamma.T @ K0 @ Gamma
    rest = [a for a in range(12) if a not in soft]
    ks = kk[soft]
    K = K0_12.copy()
    Kbb_inv = np.zeros((0, 0))
    if soft:
        Kbb = K0_12[np.ix_(soft, soft)] + np.diag(ks)
        scale = np.array([1.0 if a % 6 < 3 else abs(H[1, 5]) for a in soft])     # rotations * L
        d = np.sqrt(np.diag(Kbb))
        if not (d > 0.0).all():
            lam, vec, dd = np.zeros(1), np.eye(len(soft))[:, [int(np.argmin(d))]], np.ones(len(soft))
        else:
            lam, vec = eigh(Kbb / np.outer(d, d))
            dd = d
        if lam[0] < MECHANISM_TOL * max(lam[-1], 1.0):
            v = np.abs(vec[:, 0] / dd) * scale
            raise ValueError([("i" if a < 6 else "j", a % 6) for n, a in enumerate(soft) if v[n] > 1.0e-6 * v.max()])
        Kbb_inv = _spd_inverse(Kbb)
        Kbr = K0_12[np.ix_(soft, rest)]
        K[np.ix_(rest, rest)] = K0_12[np.ix_(rest, rest)] - Kbr.T @ Kbb_inv @ Kbr
        K[np.ix_(soft, rest)] = (ks[:, None] * Kbb_inv) @ Kbr
        K[np.ix_(rest, soft)] = K[np.ix_(soft, rest)].T
        K[np.ix_(soft, soft)] = np.diag(ks) - (ks[:, None] * Kbb_inv) * ks[None, :]
    parts = {"k": kk, "stiff": stiff, "soft": soft, "rest": rest, "K0": K0, "K0_12": K0_12, "Kbb_inv": Kbb_inv}
    return K, parts


@dataclass
class MemberRecovery:
    end_forces: np.ndarray        # (12,) local, on the member
    x: np.ndarray                 # (n,)
    forces: np.ndarray            # (n, 6) N, Vy, Vz, T, My, Mz
    displacements: np.ndarray     # (n, 6) u, v, w, theta_x, theta_y, theta_z of the member axis, local
    force_closure: np.ndarray     # (6,) resultants integrated from end i to x = L, minus the end forces at j
    disp_closure: np.ndarray      # (6,) displacements integrated from end i to x = L (+ spring j), minus end j


@dataclass
class MemberElement:
    id: str
    type: str
    nodes: tuple[str, str]
    dofs: np.ndarray              # (12,) global DOF numbers of the two NODES
    length: float
    flex: SectionFlex
    k_end_i: tuple
    k_end_j: tuple
    H: np.ndarray                 # (6,6) rigid transfer end i -> end j in local axes
    k_local: np.ndarray           # (12,12) at the (node side of the) flexible ends, local axes
    G: np.ndarray                 # (12,12) u_ends,local = G @ u_nodes,global   (rotation + rigid offsets)
    k_global: np.ndarray          # (12,12) = G^T k_local G
    parts: dict = field(repr=False, default_factory=dict)
    active: bool = True

    @property
    def is_truss(self) -> bool:
        return self.type in TRUSS_TYPES

    # ------------------------------------------------------------------ loads
    def _fixed_end(self, loads: SpanLoads):
        """(f_fix, f0, ub0): fixed-end forces at the node-side DOFs; fixed-end forces of beam + stiff
        springs with the soft beam-side DOFs held; displacement of those DOFs when released."""
        p = self.parts
        soft = p["soft"]
        if not loads:
            return np.zeros(12), np.zeros(12), np.zeros(len(soft))
        P0 = loads.resultant()                               # = resultants at x = 0- of the cantilever
        cant = integrate(self.length, self.flex, loads, P0, np.zeros(6))
        u0 = cant.y_end.copy()                               # end-j displacement of the clamped-free member
        for a in p["stiff"]:
            if a < 6:                                        # spring at i carries P0[a] in the cantilever state
                u0 += self.H[:, a] * (P0[a] / p["k"][a])
        Qj = -p["K0"] @ u0
        f0 = np.concatenate([-self.H.T @ Qj - P0, Qj])
        f = f0.copy()
        ub0 = np.zeros(0)
        if soft:
            ub0 = -p["Kbb_inv"] @ f0[soft]
            f = f0 + p["K0_12"][:, soft] @ ub0
            f[soft] = -p["k"][soft] * ub0                    # force through a soft spring (exactly 0 for k = 0)
        return f, f0, ub0

    def fixed_end(self, loads: SpanLoads) -> np.ndarray:
        """Fixed-end forces (12,) of the span loads: forces on the member at its ends, local axes, with
        both nodes held."""
        return self._fixed_end(loads)[0]

    def equivalent_nodal_loads(self, loads: SpanLoads) -> np.ndarray:
        """(12,) global forces on the two nodes that replace the span loads."""
        return -self.G.T @ self.fixed_end(loads)

    # --------------------------------------------------------------- recovery
    def recover(self, u_nodes: np.ndarray, loads: SpanLoads, stations) -> MemberRecovery:
        """End forces and exact station values for the node displacements ``u_nodes`` (12, global).
        ``stations`` = sequence of (x, side)."""
        p = self.parts
        soft, k = p["soft"], p["k"]
        ue = self.G @ u_nodes                                # node side of the end connections
        f_fix, f0, _ = self._fixed_end(loads)
        ub = ue.copy()                                       # beam side (differs at soft connections)
        if not self.active:
            f = f_fix.copy()
        else:
            if soft:
                rest = p["rest"]
                ub[soft] = p["Kbb_inv"] @ (k[soft] * ue[soft] - p["K0_12"][np.ix_(soft, rest)] @ ue[rest] - f0[soft])
            f = p["K0_12"] @ ub + f0
            if soft:
                f[soft] = k[soft] * (ue[soft] - ub[soft])
            if self.is_truss:                                # exactly axial (same simplification as k_local)
                f = self.k_local @ ue + f_fix
        for a in p["stiff"]:
            # spring deformation = (resultant through it) / k.  At i the resultant is -f[a] and the beam
            # side is ahead of the node side; at j it is +f[a] and the node side is ahead: same formula.
            ub[a] = ue[a] - f[a] / k[a]
        run = integrate(self.length, self.flex, loads, -f[:6], ub[:6])
        n = len(stations)
        forces, disp = np.zeros((n, 6)), np.zeros((n, 6))
        for m, (x, side) in enumerate(stations):
            forces[m], disp[m] = run.at(x, side)
        closure = (run.y_end - ub[6:]) if self.active else np.zeros(6)
        return MemberRecovery(f, np.array([x for x, _ in stations]), forces, disp, run.s_end - f[6:], closure)


def section_flex(model, member, shear: bool, factors: StiffnessFactors | None = None) -> SectionFlex:
    f = factors or StiffnessFactors()
    mat, sec = model.materials[member.material], model.sections[member.section]
    stiff = {"E*A": mat.E * sec.A * f.axial, "G*J": mat.G * sec.J * f.torsional,
             "E*Iy": mat.E * sec.Iy * f.flexural, "E*Iz": mat.E * sec.Iz * f.flexural}
    bad = [k for k, v in stiff.items() if not (v > 0.0 and math.isfinite(v))]
    if bad and member.type not in TRUSS_TYPES or "E*A" in bad:
        raise refuse("section.zero-stiffness", pointer("members", member.id),
                     f"member {member.id!r} has no {', '.join(bad)} stiffness",
                     "Give the section positive A, Iy, Iz and J (and the material positive E and G).")
    def sh(As):
        return 1.0 / (mat.G * As * f.shear) if (shear and As is not None and As > 0.0) else 0.0
    inv = {k: (1.0 / v if k not in bad else 0.0) for k, v in stiff.items()}
    return SectionFlex(inv["E*A"], inv["G*J"], inv["E*Iy"], inv["E*Iz"], sh(sec.Asy), sh(sec.Asz))


def build_element(model, member, factors: StiffnessFactors | None = None, active: bool = True) -> MemberElement:
    """Element matrices of one member.  ``factors`` scale the section stiffnesses; ``active = False``
    gives a member without stiffness (hook for tension/compression-only members that are switched off:
    its span loads still go to the nodes)."""
    L = float(member.length)
    truss = member.type in TRUSS_TYPES
    flex = section_flex(model, member, model.analysis.shear_deformation, factors)
    H = rigid_transfer([L, 0.0, 0.0])
    if truss:
        # Axial-only: the member is a simply supported beam between its ends (both bending rotations
        # released at both ends); the twist is held at end i only, so no torque is transmitted.
        if flex.bend_y == 0.0 or flex.bend_z == 0.0 or flex.torsion == 0.0:
            # a truss section without bending/torsion properties: any positive value serves the
            # statically determinate span; use the axial flexibility scaled to the length
            flex = SectionFlex(flex.axial, flex.torsion or flex.axial / L ** 2, flex.bend_y or flex.axial / L ** 2,
                               flex.bend_z or flex.axial / L ** 2, flex.shear_y, flex.shear_z)
        k_i = (RIGID, RIGID, RIGID, RIGID, 0.0, 0.0)
        k_j = (RIGID, RIGID, RIGID, 0.0, 0.0, 0.0)
    else:
        k_i, k_j = tuple(member.k_end_i), tuple(member.k_end_j)
    Fb = cantilever_flexibility(L, flex)
    try:
        k_local, parts = release_stiffness(Fb, H, k_i, k_j)
    except ValueError as e:
        names = ", ".join(f"{end}.{DOFS[a]}" for end, a in e.args[0])
        raise refuse("member.release-mechanism", pointer("members", member.id, "releases"),
                     f"the end releases of member {member.id!r} leave it free to move without deforming "
                     f"(member-level mechanism); releases involved: {names}",
                     "Keep enough end connections to hold the member: e.g. do not release the same "
                     "translation or the twist at both ends, or both end moments together with a shear.") from None
    if truss:
        k_local = np.zeros((12, 12))                # exactly axial: EA / L
        k_local[np.ix_([0, 6], [0, 6])] = np.array([[1.0, -1.0], [-1.0, 1.0]]) / (flex.axial * L)
    if not active:
        k_local = np.zeros((12, 12))
    R = np.asarray(member.R, float)
    G = np.zeros((12, 12))
    for blk, off in ((0, member.offset_i), (6, member.offset_j)):
        E = rigid_transfer(off)                    # node -> flexible end, global components
        rot = np.zeros((6, 6))
        rot[:3, :3] = R
        rot[3:, 3:] = R
        G[blk:blk + 6, blk:blk + 6] = rot @ E
    k_global = G.T @ k_local @ G        # not symmetrised: symmetry is an invariant that is checked
    dofs = np.concatenate([6 * model.nodes[member.i].index + np.arange(6),
                           6 * model.nodes[member.j].index + np.arange(6)])
    return MemberElement(member.id, member.type, (member.i, member.j), dofs, L, flex, k_i, k_j, H, k_local,
                         G, k_global, parts, active)
