"""The second-order member piece: a straight prismatic beam-column with a given axial force.

Derivation: docs/THEORY.md section 10 (nothing is taken from a reference).  Summary:

* Small-displacement second-order theory.  The axial force N(x) (tension positive, known from the
  previous iterate) keeps the direction of the undeformed member axis; equilibrium of a slice is
  written in its displaced position.  With V the force component along the undeformed local y (or z)
  on the +x face, in one bending plane

      v'     = (theta + beta V) / (1 + beta N)        beta = 1 / (G As)   (0: no shear deformation)
      theta' = M / EI
      M'     = -V + N v' = (N theta - V) / (1 + beta N)
      V'     = -q(x)

  (the shear strain is produced by the force normal to the deformed AXIS, V - N v': this is what the
  second-order work of the axial force, -N * integral(v'^2 / 2), gives; THEORY 10.2).  The xz plane obeys
  the same equations in the variables (w, -theta_y, -My, Vz).
* The state z = (v, theta, M, V) is carried along the piece by the matrix exponential of the
  first-order system above, augmented by two states (1, x) for a linearly varying load.  For a
  constant N this is the exact solution operator (the "stability functions" in numerical form, valid
  for compression, tension and N = 0 alike).  For N varying along the piece (axial member loads)
  the two-term Magnus series is used (fourth-order accurate; the subdivision of the member is then
  refined until the results stop changing, see ``nonlinear.py``).
* End stiffness and fixed-end forces follow from the transfer relation by solving for the start
  resultants; end connections (releases, partial springs) are springs between the node side and the
  beam side of an end DOF and are condensed exactly; rigid end offsets add the exact linear
  kinematics plus the second-order work of the axial end force on the rotating rigid arm.
* Axial and torsional behaviour stay first order (uniform bar, uniform St. Venant torsion).

Not included (limits): the effect of bending moments, shears and torque on the stiffness (so
lateral-torsional and torsional / flexural-torsional buckling are not captured), bowing (shortening
due to bending), large rotations.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import eigvalsh, expm

from ..model import RIGID
from .beam import SpanLoads
from .element import TRUSS_TYPES, MemberRecovery, rigid_transfer, section_flex

__all__ = ["BeamColumn", "PieceRecovery", "axial_state", "consistent_geometric_stiffness", "STABLE_TOL"]

SQ3 = math.sqrt(3.0)
_T4, _W4 = np.polynomial.legendre.leggauss(4)
GAUSS4 = tuple(zip((0.5 * (1.0 + _T4)).tolist(), (0.5 * _W4).tolist()))    # points and weights on [0, 1]
STABLE_TOL = 1.0e-10      # smallest scaled eigenvalue of the condensed block below which the piece is unstable


class _Plane:
    """One bending plane of a piece: transfer of z = (v, theta, M, V) along x under the axial force
    N(x) = n0 - na x - ns x^2 / 2 and the load q(x) = qa + qs x."""

    def __init__(self, h, EI, beta, axial, qa=0.0, qs=0.0):
        self.h, self.EI, self.beta = float(h), float(EI), float(beta)
        self.n0, self.na, self.ns = (float(a) for a in axial[:3])
        self.om = self.beta * self.EI / self.h ** 2
        self.qa, self.qs = float(qa), float(qs)
        self.S = np.array([self.h, 1.0, self.EI / self.h, self.EI / self.h ** 2])    # z = S * scaled state
        self._phi: dict[float, np.ndarray] = {}
        self.constant = self.na == 0.0 and self.ns == 0.0
        # 1 + beta N must stay positive (N > -G As): below that the member is unstable in shear
        self.bad = not all(1.0 + self.beta * self.N(x) > 0.0 for x in (0.0, 0.5, 1.0))

    def N(self, xi: float) -> float:
        x = xi * self.h
        return self.n0 - self.na * x - 0.5 * self.ns * x * x

    def _A(self, xi: float) -> np.ndarray:
        """Coefficient matrix of d/dxi (v/h, theta, M h/EI, V h^2/EI, 1, xi)."""
        nu = self.N(xi) * self.h ** 2 / self.EI
        c = 1.0 + self.om * nu
        A = np.zeros((6, 6))
        A[0, 1] = 1.0 / c
        A[0, 3] = self.om / c
        A[1, 2] = 1.0
        A[2, 1] = nu / c
        A[2, 3] = -1.0 / c
        A[3, 4] = -self.qa * self.h ** 3 / self.EI
        A[3, 5] = -self.qs * self.h ** 4 / self.EI
        A[5, 4] = 1.0
        return A

    def phi(self, xi: float) -> np.ndarray:
        """Transfer matrix of the scaled, augmented state from 0 to xi (0 <= xi <= 1)."""
        xi = float(xi)
        if xi not in self._phi:
            if self.constant:
                om = xi * self._A(0.0)
            else:
                # two-term Magnus series with the two Gauss points of [0, xi] (THEORY 10.4)
                a1, a2 = self._A(xi * (0.5 - SQ3 / 6.0)), self._A(xi * (0.5 + SQ3 / 6.0))
                om = 0.5 * xi * (a1 + a2) - (SQ3 * xi * xi / 12.0) * (a1 @ a2 - a2 @ a1)
            self._phi[xi] = expm(om)
        return self._phi[xi]

    def state(self, z0: np.ndarray, xi: float) -> np.ndarray:
        """(v, theta, M, V) at xi for the start state z0 (just inside x = 0)."""
        y = self.phi(xi) @ np.concatenate([np.asarray(z0, float) / self.S, [1.0, 0.0]])
        return y[:4] * self.S

    def slope(self, z: np.ndarray, xi: float) -> float:
        """v' at xi for the state z."""
        return (z[1] + self.beta * z[3]) / (1.0 + self.beta * self.N(xi))

    def stiffness(self):
        """(k, f): forces (F0v, F0theta, F1v, F1theta) on the piece = k (v0, theta0, v1, theta1) + f."""
        P = self.phi(1.0)
        d, s = [0, 1], [3, 2]                              # s = (V, M): work-conjugate to (v, theta)
        A, B, C, D = P[np.ix_(d, d)], P[np.ix_(d, s)], P[np.ix_(s, d)], P[np.ix_(s, s)]
        pd, ps = P[d, 4], P[s, 4]
        Bi = np.linalg.inv(B)
        k = np.block([[Bi @ A, -Bi], [C - D @ Bi @ A, D @ Bi]])
        f = np.concatenate([Bi @ pd, ps - D @ Bi @ pd])
        sd = np.array([self.h, 1.0, self.h, 1.0])
        sf = np.array([self.EI / self.h ** 2, self.EI / self.h] * 2)
        return (sf[:, None] * k) / sd[None, :], sf * f


def axial_state(loads: SpanLoads, n_start: float, length: float):
    """(n0, na, ns, n_i, n_j): axial force N(x) = n0 - na x - ns x^2/2 inside a piece whose axial force
    just inside the start is ``n_start``; n_i / n_j = the axial force just outside the start / the end
    (they differ from the inside values by the axial point loads at the ends)."""
    dist, p_i, p_j = _piece_loads(loads, length)
    n_end = n_start - dist[0][0] * length - 0.5 * dist[1][0] * length ** 2
    return float(n_start), float(dist[0][0]), float(dist[1][0]), float(n_start + p_i[0]), float(n_end - p_j[0])


def _piece_loads(loads: SpanLoads | None, length: float):
    """((wa, ws), P_i, P_j): a piece carries line loads over its whole length (intensity wa + ws x) and
    point loads at its two ends only (the subdivision puts a piece boundary at every load break)."""
    wa, ws, p_i, p_j = np.zeros(3), np.zeros(3), np.zeros(6), np.zeros(6)
    if loads:
        for x, f in loads.points:
            if x <= 0.5 * length:
                p_i += np.asarray(f, float)
            else:
                p_j += np.asarray(f, float)
        for x1, x2, w1, w2 in loads.dists:
            if x2 <= x1:
                continue
            slope = (np.asarray(w2, float) - np.asarray(w1, float)) / (x2 - x1)
            wa += np.asarray(w1, float) - slope * x1
            ws += slope
    return (wa, ws), p_i, p_j


@dataclass
class PieceRecovery(MemberRecovery):
    n_start: float = 0.0              # axial force just inside the start (from the end forces of THIS solution)
    n_min: float = 0.0                # smallest / largest axial force along the piece
    n_max: float = 0.0
    node_forces: np.ndarray | None = None       # (12,) forces on the piece at its two NODES, global
    couple: np.ndarray | None = None            # (3,) global: moment about the start of all forces on the piece
    #                                             taken at the undeformed positions (= integral of N x slope)
    deformed_residual: float = 0.0    # moment equilibrium of the piece on its deformed axis (should be 0)
    moment_scale: float = 0.0


@dataclass
class BeamColumn:
    """One piece (a member or a part of a subdivided member).  Same public attributes as
    ``element.MemberElement`` where the assembly needs them (``dofs, k_global, G, nodes, is_truss``)."""
    id: str
    type: str
    nodes: tuple
    dofs: np.ndarray
    length: float
    flex: object
    kk: np.ndarray                    # (12,) connection stiffnesses (inf = rigid)
    axial: tuple                      # (n0, na, ns, n_i, n_j) used in the bending equations
    loads: tuple                      # ((wa, ws), P_i, P_j)
    R: np.ndarray
    G: np.ndarray
    k_local: np.ndarray               # (12,12) node side, local axes
    f_fixed: np.ndarray               # (12,) fixed-end forces, node side, local axes
    k_global: np.ndarray              # (12,12) = G^T k_local G + rigid-arm terms
    arm: tuple                        # (3x3, 3x3): geometric stiffness of the rigid arms on the node rotations
    stable: bool                      # False: the piece itself has buckled between its end connections
    parts: dict = field(repr=False, default_factory=dict)
    active: bool = True

    @property
    def is_truss(self) -> bool:
        return self.type in TRUSS_TYPES

    # ------------------------------------------------------------------ build
    @classmethod
    def build(cls, model, member, loads: SpanLoads | None = None, n_start: float = 0.0, factors=None,
              axial: tuple | None = None) -> "BeamColumn":
        """``n_start`` = axial force just inside the start (the force along the piece then follows from
        the axial loads on it); ``axial`` overrides the whole axial state (n0, na, ns, n_i, n_j)."""
        h = float(member.length)
        truss = member.type in TRUSS_TYPES
        flex = section_flex(model, member, model.analysis.shear_deformation, factors)
        pl = _piece_loads(loads, h)
        if axial is None:
            axial = axial_state(loads, n_start, h) if loads else (float(n_start), 0.0, 0.0, float(n_start), float(n_start))
        kk = np.array(list(member.k_end_i) + list(member.k_end_j), float)
        if truss:
            kk[[3, 9]] = RIGID                              # no torsional stiffness at all (set to zero below)
        (wa, ws), p_i, p_j = pl
        Kb, f0 = np.zeros((12, 12)), np.zeros(12)
        ea = 1.0 / (flex.axial * h)
        Kb[np.ix_([0, 6], [0, 6])] = ea * np.array([[1.0, -1.0], [-1.0, 1.0]])
        if not truss:
            Kb[np.ix_([3, 9], [3, 9])] = np.array([[1.0, -1.0], [-1.0, 1.0]]) / (flex.torsion * h)
        # uniform bar with both ends held under the line load wa + ws x (THEORY 10.5)
        f0[0] = -(wa[0] * h / 2.0 + ws[0] * h * h / 6.0)
        f0[6] = -(wa[0] * h / 2.0 + ws[0] * h * h / 3.0)
        planes = {}
        stable = True
        for name, idx, sg, EI, beta, c in (("xy", [1, 5, 7, 11], [1.0, 1.0, 1.0, 1.0], 1.0 / flex.bend_z, flex.shear_y, 1),
                                           ("xz", [2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0], 1.0 / flex.bend_y, flex.shear_z, 2)):
            p = _Plane(h, EI, beta, axial, wa[c], ws[c])
            planes[name] = p
            if p.bad:
                stable = False
                continue
            k4, f4 = p.stiffness()
            sg = np.array(sg)
            Kb[np.ix_(idx, idx)] = sg[:, None] * k4 * sg[None, :]
            f0[idx] = sg * f4
        f0[:6] -= p_i
        f0[6:] -= p_j
        # ---- end connections: springs between the node side and the beam side, condensed exactly
        c = [a for a in range(12) if kk[a] != RIGID]
        r = [a for a in range(12) if kk[a] == RIGID]
        K, f, S = Kb.copy(), f0.copy(), np.zeros((0, 0))
        if c:
            D = kk[c]
            M = Kb[np.ix_(c, c)] + np.diag(D)
            dg = np.diag(M)
            if not (dg > 0.0).all():
                stable = False
                S = np.linalg.pinv(M)
            else:
                sc = np.sqrt(dg)
                Ms = M / np.outer(sc, sc)
                if eigvalsh(0.5 * (Ms + Ms.T))[0] < STABLE_TOL:
                    stable = False
                    S = np.linalg.pinv(M)
                else:
                    S = np.linalg.inv(Ms) / np.outer(sc, sc)
            Kcr = Kb[np.ix_(c, r)]
            K[np.ix_(r, r)] = Kb[np.ix_(r, r)] - Kcr.T @ S @ Kcr
            K[np.ix_(c, r)] = (D[:, None] * S) @ Kcr
            K[np.ix_(r, c)] = Kcr.T @ (S * D[None, :])
            K[np.ix_(c, c)] = (D[:, None] * S) @ Kb[np.ix_(c, c)]
            f[r] = f0[r] - Kcr.T @ S @ f0[c]
            f[c] = (D[:, None] * S) @ f0[c]
        Rm = np.asarray(member.R, float)
        G = np.zeros((12, 12))
        for blk, off in ((0, member.offset_i), (6, member.offset_j)):
            rot = np.zeros((6, 6))
            rot[:3, :3] = Rm
            rot[3:, 3:] = Rm
            G[blk:blk + 6, blk:blk + 6] = rot @ rigid_transfer(off)
        kg = G.T @ K @ G
        arm = (_arm(member.offset_i, -axial[3] * Rm[0]), _arm(member.offset_j, axial[4] * Rm[0]))
        kg[3:6, 3:6] += arm[0]
        kg[9:12, 9:12] += arm[1]
        dofs = np.concatenate([6 * model.nodes[member.i].index + np.arange(6),
                               6 * model.nodes[member.j].index + np.arange(6)])
        parts = {"Kb": Kb, "f0": f0, "c": c, "r": r, "S": S, "planes": planes,
                 "offsets": (np.asarray(member.offset_i, float), np.asarray(member.offset_j, float))}
        return cls(member.id, member.type, (member.i, member.j), dofs, h, flex, kk, tuple(axial), pl, Rm, G, K, f,
                   kg, arm, stable, parts)

    # ------------------------------------------------------------------ loads
    def equivalent_nodal_loads(self, loads=None) -> np.ndarray:
        """(12,) global forces on the two nodes that replace the loads of the piece."""
        return -self.G.T @ self.f_fixed

    # --------------------------------------------------------------- recovery
    def _ends(self, u_nodes: np.ndarray):
        """(ue, ub, F): node-side and beam-side end displacements (local), forces on the piece at its ends."""
        p = self.parts
        c, r, Kb, f0 = p["c"], p["r"], p["Kb"], p["f0"]
        ue = self.G @ u_nodes
        ub = ue.copy()
        if c:
            D = self.kk[c]
            ub[c] = p["S"] @ (D * ue[c] - Kb[np.ix_(c, r)] @ ue[r] - f0[c])
        F = Kb @ ub + f0
        if c:
            F[c] = self.kk[c] * (ue[c] - ub[c])            # force through the connection (exactly 0 for a release)
        return ue, ub, F

    def _starts(self, ub, F):
        (_, _), p_i, _ = self.loads
        zy = np.array([ub[1], ub[5], -F[5] - p_i[5], -F[1] - p_i[1]])
        zz = np.array([ub[2], -ub[4], F[4] + p_i[4], -F[2] - p_i[2]])
        return zy, zz

    def recover(self, u_nodes: np.ndarray, stations=(), checks: bool = False) -> PieceRecovery:
        """End forces and station values for the node displacements ``u_nodes`` (12, global);
        ``stations`` = sequence of (x, side) in the piece's own x."""
        (wa, ws), p_i, p_j = self.loads
        h = self.length
        py, pz = self.parts["planes"]["xy"], self.parts["planes"]["xz"]
        ue, ub, F = self._ends(u_nodes)
        zy0, zz0 = self._starts(ub, F)
        n0 = -F[0] - p_i[0]                                  # axial force just inside the start, from THIS solution
        T = -F[3] - p_i[3]
        gj = 0.0 if self.is_truss else self.flex.torsion

        def axial(x):
            return n0 - wa[0] * x - ws[0] * x * x / 2.0

        def inside(x):
            xi = min(max(x / h, 0.0), 1.0)
            zy, zz = py.state(zy0, xi), pz.state(zz0, xi)
            s = np.array([axial(x), zy[3], zz[3], T, -zz[2], zy[2]])
            y = np.array([ub[0] + (n0 * x - wa[0] * x * x / 2.0 - ws[0] * x ** 3 / 6.0) * self.flex.axial,
                          zy[0], zz[0], ub[3] + T * x * gj, -zz[1], zy[1]])
            return s, y

        s1, y1 = inside(h)
        n = len(stations)
        forces, disp = np.zeros((n, 6)), np.zeros((n, 6))
        for m, (x, side) in enumerate(stations):
            if x <= 0.0:
                if side < 0:
                    forces[m], disp[m] = -F[:6], ub[:6]
                else:
                    forces[m], disp[m] = inside(0.0)
            elif x >= h:
                if side > 0:
                    forces[m], disp[m] = F[6:], ub[6:]
                else:
                    forces[m], disp[m] = s1, y1
            else:
                forces[m], disp[m] = inside(x)
        ends = (axial(0.0), axial(h), axial(0.0) + p_i[0], axial(h) - p_j[0])
        if ws[0] != 0.0 and 0.0 < -wa[0] / ws[0] < h:        # extremum of a quadratic N inside the piece
            ends += (axial(-wa[0] / ws[0]),)
        node_forces = self.G.T @ F
        node_forces[3:6] += self.arm[0] @ u_nodes[3:6]
        node_forces[9:12] += self.arm[1] @ u_nodes[9:12]
        rec = PieceRecovery(F, np.array([x for x, _ in stations]), forces, disp, (s1 - p_j) - F[6:], y1 - ub[6:],
                            n_start=float(n0), n_min=float(min(ends)), n_max=float(max(ends)), node_forces=node_forces)
        if checks:
            # moment equilibrium of the piece written on its deformed axis, with the integral of N v'
            # evaluated by quadrature of the recovered field (independent of the stiffness algebra)
            if py.constant:                                  # constant N: the integral of N v' is N (v_j - v_i)
                iy, iz = py.n0 * (y1[1] - zy0[0]), pz.n0 * (y1[2] - zz0[0])
            else:
                iy = iz = 0.0
                for xi, w in GAUSS4:
                    iy += w * h * py.N(xi) * py.slope(py.state(zy0, xi), xi)
                    iz += w * h * pz.N(xi) * pz.slope(pz.state(zz0, xi), xi)
            res_y = (F[5] + F[11] + p_i[5] + p_j[5] - (F[1] + p_i[1]) * h
                     - (wa[1] * h * h / 2.0 + ws[1] * h ** 3 / 6.0) - iy)
            res_z = (-F[4] - F[10] - p_i[4] - p_j[4] - (F[2] + p_i[2]) * h
                     - (wa[2] * h * h / 2.0 + ws[2] * h ** 3 / 6.0) - iz)
            rec.couple = self.R.T @ np.array([0.0, -iz, iy])
            rec.deformed_residual = float(max(abs(res_y), abs(res_z)))
            rec.moment_scale = float(max(np.abs(F[[4, 5, 10, 11]]).max(), np.abs(F[[1, 2, 7, 8]]).max() * h,
                                         abs(iy), abs(iz)))
        return rec

    def slopes(self, xi: float) -> tuple[np.ndarray, np.ndarray]:
        """(sv, sw): v'(xi) and w'(xi) of the piece without its loads, as linear forms (12,) of the
        node-side local end displacements (the displacement field behind the consistent geometric
        stiffness)."""
        p = self.parts
        c, r, Kb = p["c"], p["r"], p["Kb"]
        Tb = np.eye(12)                                      # beam-side = Tb @ node-side
        if c:
            Tb[c, :] = 0.0
            Tb[np.ix_(c, c)] = p["S"] * self.kk[c][None, :]
            Tb[np.ix_(c, r)] = -p["S"] @ Kb[np.ix_(c, r)]
        F = Kb @ Tb
        out = []
        for plane, Z in ((p["planes"]["xy"], np.array([Tb[1], Tb[5], -F[5], -F[1]])),
                         (p["planes"]["xz"], np.array([Tb[2], -Tb[4], F[4], -F[2]]))):
            z = plane.S[:, None] * (plane.phi(xi)[:4, :4] @ (Z / plane.S[:, None]))
            out.append((z[1] + plane.beta * z[3]) / (1.0 + plane.beta * plane.N(xi)))
        return out[0], out[1]


def _arm(offset, f) -> np.ndarray:
    """Geometric stiffness (3x3, node rotations, global) of a rigid arm ``offset`` (node -> flexible
    end) that transmits the force ``f`` to the member end: second-order work of f on the arm end,
    u_end = u + theta x o + theta x (theta x o) / 2   (THEORY 10.7)."""
    o, f = np.asarray(offset, float), np.asarray(f, float)
    if not o.any() or not f.any():
        return np.zeros((3, 3))
    return 0.5 * (np.outer(o, f) + np.outer(f, o)) - float(o @ f) * np.eye(3)


def consistent_geometric_stiffness(piece: BeamColumn, axial: tuple) -> np.ndarray:
    """Geometric stiffness (12x12, nodal global DOFs) of an UNLOADED, N = 0 piece for the axial force
    ``axial`` = (n0, na, ns, n_i, n_j): the second-order work of N on the first-order displacement
    field of the piece, kg[a, b] = integral N (v_a' v_b' + w_a' w_b') dx, plus the rigid-arm terms.
    It equals dK/dN of the exact piece at N = 0 (THEORY 11.1); it is LINEAR in the axial force."""
    h = piece.length
    n0, na, ns, n_i, n_j = axial
    kg = np.zeros((12, 12))
    for xi, w in GAUSS4:
        x = xi * h
        sv, sw = piece.slopes(xi)
        kg += w * h * (n0 - na * x - 0.5 * ns * x * x) * (np.outer(sv, sv) + np.outer(sw, sw))
    out = piece.G.T @ kg @ piece.G
    ex = piece.R[0]
    oi, oj = piece.parts["offsets"]
    out[3:6, 3:6] += _arm(oi, -n_i * ex)
    out[9:12, 9:12] += _arm(oj, n_j * ex)
    return out
