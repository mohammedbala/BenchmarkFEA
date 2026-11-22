"""Direct integration of the member equations (multi-point BVP, ``scipy.integrate.solve_bvp``).

Equations (derived here; signs are those of docs/CONVENTIONS.md section 4)
==========================================================================
Bending in the local x-y plane.  Unknowns along x: deflection ``v`` (along +y), rotation of the
cross-section ``th`` (about +z), bending moment ``M`` (= Mz, about +z on the +x face) and ``Q``, the
force along the UNDEFORMED +y direction on the +x face.  ``N`` is the axial force (tension positive),
``w(x)`` the line load along +y.

*Equilibrium of a slice dx in its displaced position.*  On the +x face at (x+dx, v+dv) act the force
(N, Q) and the moment M+dM; on the -x face the opposites of (N, Q) and M.  Forces along y:
``dQ + w dx = 0``.  Moments about the -x face centre: ``dM + Q dx - N dv = 0`` (the lever arms of Q
and N are dx and dv).  Hence

    Q' = -w                         M' = -Q + N v'

With N = 0, or in a first-order analysis (equilibrium written on the undisplaced slice, lever arm of
N equal to zero), the second equation is ``M' = -Q`` and Q is the shear force Vy.

*Material law.*  Plane sections: axial strain = -y th', so ``M = -int(y sigma) dA = E I th'``:

    th' = M / EI

*Shear deformation.*  The shear angle is ``gamma = v' - th`` and ``V = GAs gamma`` with V the shear
force carried by the section, so ``v' = th + V/GAs``.  Without shear deformation (``GAs = None``)
``v' = th``.  In a first-order analysis V = Q.  In a second-order analysis the transverse force on the
+x face has to be resolved in the displaced configuration and two definitions of "the shear force" are
in use; they differ only when shear deformation AND axial force are both present:

    "engesser":  V = Q - N v'     (force normal to the displaced member axis)
                 ->  v' = (th + Q/GAs) / (1 + N/GAs)
    "haringx":   V = Q - N th     (force in the plane of the rotated cross-section)
                 ->  v' = th + (Q - N th)/GAs

The mechanics does not single one out, so a second-order problem with shear deformation and a nonzero
N must name the model explicitly (``shear_model=``); otherwise :class:`OracleError` is raised.

Bending in the local x-z plane obeys  ``w' = -th_y + Vz/GAs, th_y' = My/EIy, My' = +Vz, Vz' = -w_z``
(CONVENTIONS section 4).  The substitution th = -th_y, M = -My turns it into the x-y system, which is
how ``plane="xz"`` is implemented: rotations, moments, rotational stiffnesses' signs are mapped on
input and output, nothing else changes.

*Joints (stations).*  At a station with applied force P (along +y), applied moment Mo (about +z),
translational spring k and rotational spring kr to ground, equilibrium of the joint under the forces
of the left part (-Q-, -M-), the right part (+Q+, +M+), the load and the springs gives

    Q+ = Q- - P + k v              M+ = M- - Mo + kr th

A fixed translation (rotation) replaces the first (second) relation by v = v0 (th = th0), and the
reaction follows from the same equilibrium: ``R = Q- - Q+ - P``, ``Mr = M- - M+ - Mo``.  An internal
hinge with rotational stiffness kh (0 = free hinge) replaces continuity of th by
``M- = M+ = kh (th+ - th-)``.  At the two ends the missing side contributes Q = M = 0.

Axial / torsional bar (:class:`Bar`):  ``u' = N/EA, N' = -p`` with the joint relation
``N+ = N- - P + k u`` (the same derivation with one force); torsion is the same system with GJ.

Numerical method
================
The member is cut at every station, load end and property change.  Each segment is mapped to [0, 1]
and the states of all segments are stacked into one two-point BVP whose boundary conditions are the
joint relations above.  States are scaled (v/L, th, M L/EI, Q L^2/EI) so that the tolerance of
``solve_bvp`` acts on quantities of comparable size.  The result reports the maximum relative
collocation residual (``rms_residual``) achieved; ``solve`` raises if the tolerance is not met.
``tol`` bounds that RESIDUAL, not the error of the solution: measured differences between solutions
at ``tol`` and at a much tighter tolerance are up to about 20 x tol, so treat the results as accurate
to about 100 x tol (1e-8 relative with the default ``tol = 1e-10``).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_bvp

__all__ = ["Beam", "Bar", "BeamSolution", "BarSolution", "OracleError"]

FIXED = "fixed"


class OracleError(ValueError):
    """The problem is ill-posed, ambiguous, or the integration did not reach the tolerance."""


@dataclass
class _Station:
    x: float
    v: object = None          # None (free) | "fixed" | spring stiffness
    r: object = None
    v0: float = 0.0           # imposed translation / rotation of a fixed DOF
    r0: float = 0.0
    P: float = 0.0
    Mo: float = 0.0
    hinge: object = None      # None | rotational stiffness (0.0 = free hinge)


def _piecewise(spec, L, name):
    """scalar | [(x1, x2, value), ...] -> list of (x1, x2, value) covering [0, L]."""
    if spec is None:
        return None
    if np.isscalar(spec):
        return [(0.0, L, float(spec))]
    out = sorted((float(a), float(b), float(v)) for a, b, v in spec)
    x = 0.0
    for a, b, _ in out:
        if abs(a - x) > 1e-12 * L or b <= a:
            raise OracleError(f"{name}: the pieces must cover [0, L] without gaps or overlaps")
        x = b
    if abs(x - L) > 1e-12 * L:
        raise OracleError(f"{name}: the pieces must end at L")
    return out


def _value(pieces, xm):
    for a, b, v in pieces:
        if a <= xm <= b:
            return v
    raise OracleError("position outside the member")


class _Line:
    """Shared bookkeeping: stations, line loads, segmentation."""

    def __init__(self, length):
        if not length > 0.0:
            raise OracleError("length must be positive")
        self.L = float(length)
        self._st: dict[float, _Station] = {}
        self._dist: list[tuple[float, float, float, float]] = []

    def _station(self, x) -> _Station:
        x = float(x)
        if not -1e-12 * self.L <= x <= self.L * (1 + 1e-12):
            raise OracleError(f"position {x} outside [0, {self.L}]")
        x = min(max(x, 0.0), self.L)
        for key in self._st:
            if abs(key - x) <= 1e-12 * self.L:
                return self._st[key]
        self._st[x] = _Station(x)
        return self._st[x]

    def dist(self, x1, x2, w1, w2=None):
        """Line load varying linearly from ``w1`` at x1 to ``w2`` at x2 (default uniform)."""
        x1, x2 = float(x1), float(x2)
        if not 0.0 <= x1 < x2 <= self.L * (1 + 1e-12):
            raise OracleError("line load needs 0 <= x1 < x2 <= L")
        self._dist.append((x1, min(x2, self.L), float(w1), float(w1 if w2 is None else w2)))
        return self

    def _breaks(self, extra=()):
        xs = {0.0, self.L} | set(self._st) | {a for a, *_ in self._dist} | {b for _, b, *_ in self._dist} | set(extra)
        xs = sorted(xs)
        out = [xs[0]]
        for x in xs[1:]:
            if x - out[-1] > 1e-12 * self.L:
                out.append(x)
        return np.array(out)

    def _load(self, x, k, breaks):
        """Line load intensity at positions x (array) inside segment k."""
        a, b = breaks[k], breaks[k + 1]
        xm = 0.5 * (a + b)
        w = np.zeros_like(x)
        for x1, x2, w1, w2 in self._dist:
            if x1 <= xm <= x2:
                w = w + w1 + (w2 - w1) * (x - x1) / (x2 - x1)
        return w

    def _at(self, x):
        for key, st in self._st.items():
            if abs(key - x) <= 1e-12 * self.L:
                return st
        return _Station(x)


# ====================================================================================== bending
class Beam(_Line):
    """Straight beam / beam-column in one bending plane.

    ``EI``, ``GAs`` (``None`` = no shear deformation) and ``N`` (axial force, tension positive) are
    scalars or lists ``[(x1, x2, value), ...]`` (piecewise constant).  ``second_order=True`` includes
    the ``N v'`` term (see module docstring); with shear deformation and N != 0 the ``shear_model``
    ("engesser" | "haringx") must then be given.  ``plane`` is "xy" (v, th_z, Mz, Vy) or "xz"
    (w, th_y, My, Vz): moments, rotations and their prescribed values are given and returned with the
    signs of that plane (CONVENTIONS section 4).
    """

    def __init__(self, length, EI, GAs=None, N=0.0, second_order=False, shear_model=None, plane="xy"):
        super().__init__(length)
        if plane not in ("xy", "xz"):
            raise OracleError("plane must be 'xy' or 'xz'")
        self.plane = plane
        self._s = 1.0 if plane == "xy" else -1.0       # sign of rotation/moment quantities
        self.EI = _piecewise(EI, self.L, "EI")
        self.GAs = _piecewise(GAs, self.L, "GAs")
        self.N = _piecewise(N, self.L, "N")
        self.second_order = bool(second_order)
        if shear_model not in (None, "engesser", "haringx"):
            raise OracleError("shear_model must be 'engesser' or 'haringx'")
        if self.second_order and self.GAs is not None and any(v != 0.0 for *_, v in self.N) and shear_model is None:
            raise OracleError("second-order analysis with shear deformation and axial force: the shear force "
                              "definition is ambiguous; pass shear_model='engesser' or 'haringx'")
        self.shear_model = shear_model

    # ------------------------------------------------------------------ problem definition
    def support(self, x, v=None, r=None, settlement=0.0, rotation=0.0):
        """Support at x: ``v`` / ``r`` = None (free), "fixed" or a spring stiffness (to ground)."""
        st = self._station(x)
        if v is not None:
            st.v, st.v0 = v, float(settlement)
        if r is not None:
            st.r = r if r == FIXED else self._abs_spring(r)
            st.r0 = self._s * float(rotation)
        return self

    @staticmethod
    def _abs_spring(k):
        if k < 0.0:
            raise OracleError("negative spring stiffness")
        return float(k)

    def point(self, x, P=0.0, M=0.0):
        """Concentrated force P (along the transverse axis) and moment M (about the plane normal axis:
        +z for "xy", +y for "xz") at x."""
        st = self._station(x)
        st.P += float(P)
        st.Mo += self._s * float(M)
        return self

    def hinge(self, x, k=0.0):
        """Internal hinge at 0 < x < L with rotational stiffness k (0 = free)."""
        if not 0.0 < x < self.L:
            raise OracleError("an internal hinge must be strictly inside the member")
        self._station(x).hinge = self._abs_spring(k)
        return self

    # ------------------------------------------------------------------ solution
    def solve(self, tol=1e-10, max_nodes=200000) -> "BeamSolution":
        L = self.L
        prop_breaks = [a for p in (self.EI, self.GAs, self.N) if p for a, _, _ in p]
        br = self._breaks(prop_breaks)
        m = len(br) - 1
        h = np.diff(br)
        mid = 0.5 * (br[:-1] + br[1:])
        EI = np.array([_value(self.EI, x) for x in mid])
        c = np.array([1.0 / _value(self.GAs, x) if self.GAs else 0.0 for x in mid])
        N = np.array([_value(self.N, x) if self.second_order else 0.0 for x in mid])
        for st in self._st.values():
            if st.hinge is not None and (st.Mo != 0.0 or st.r is not None):
                raise OracleError("a moment or a rotational support exactly at an internal hinge is ambiguous "
                                  "(which side?): place it at a small distance")
        EIr = float(EI.max())
        sc = np.array([L, 1.0, EIr / L, EIr / L ** 2])            # physical = sc * scaled
        haringx = self.shear_model == "haringx"

        def phys_rhs(k, x, y):
            v, th, M, Q = y
            if haringx:
                dv = th + c[k] * (Q - N[k] * th)
            else:
                dv = (th + c[k] * Q) / (1.0 + c[k] * N[k])
            return np.array([dv, M / EI[k], -Q + N[k] * dv, -self._load(x, k, br)])

        def rhs(xi, Y):
            out = np.empty_like(Y)
            for k in range(m):
                y = Y[4 * k:4 * k + 4] * sc[:, None]
                out[4 * k:4 * k + 4] = h[k] * phys_rhs(k, br[k] + xi * h[k], y) / sc[:, None]
            return out

        stations = [self._at(x) for x in br]

        def joint(st, left, right):
            """Residuals (scaled) at a station; ``left``/``right`` = physical states or None at the ends."""
            res = []
            Qm, Mm = (left[3], left[2]) if left is not None else (0.0, 0.0)
            Qp, Mp = (right[3], right[2]) if right is not None else (0.0, 0.0)
            one = right if left is None else left
            v, th = one[0], one[1]
            inner = left is not None and right is not None
            if inner:
                res.append((right[0] - left[0]) / sc[0])
            if st.v == FIXED:
                res.append((v - st.v0) / sc[0])
            else:
                res.append((Qp - Qm + st.P - (st.v or 0.0) * v) / sc[3])
            if st.hinge is not None:
                res.append((Mm - st.hinge * (right[1] - left[1])) / sc[2])
                res.append((Mp - Mm) / sc[2])
                return res
            if inner:
                res.append(right[1] - left[1])
            if st.r == FIXED:
                res.append(th - st.r0)
            else:
                res.append((Mp - Mm + st.Mo - (st.r or 0.0) * th) / sc[2])
            return res

        def bc(Ya, Yb):
            ya = [Ya[4 * k:4 * k + 4] * sc for k in range(m)]
            yb = [Yb[4 * k:4 * k + 4] * sc for k in range(m)]
            res = joint(stations[0], None, ya[0])
            for k in range(1, m):
                res += joint(stations[k], yb[k - 1], ya[k])
            res += joint(stations[m], yb[m - 1], None)
            return np.array(res)

        xi = np.linspace(0.0, 1.0, 9)
        sol = solve_bvp(rhs, bc, xi, np.zeros((4 * m, xi.size)), tol=tol, max_nodes=max_nodes)
        if not sol.success or float(np.max(sol.rms_residuals)) > 10.0 * tol:
            raise OracleError(f"integration did not converge: {sol.message} (a mechanism, or a load at/above "
                              "the buckling load, makes the boundary value problem singular)")
        return BeamSolution(self, br, sol, sc, c, N, stations, float(np.max(sol.rms_residuals)))


@dataclass
class BeamSolution:
    beam: Beam
    breaks: np.ndarray
    _sol: object
    _sc: np.ndarray
    _c: np.ndarray
    _N: np.ndarray
    _stations: list
    rms_residual: float
    reactions: dict = field(init=False)

    def __post_init__(self):
        # reaction (force, moment) BY the support ON the beam at every station with a support
        self.reactions = {}
        s = self.beam._s
        for x, st in zip(self.breaks, self._stations):
            if st.v is None and st.r is None:
                continue
            left = self._state(x, -1) if x > 0.0 else None
            right = self._state(x, +1) if x < self.beam.L else None
            Qm, Mm = (left[3], left[2]) if left is not None else (0.0, 0.0)
            Qp, Mp = (right[3], right[2]) if right is not None else (0.0, 0.0)
            one = right if left is None else left
            R = 0.0 if st.v is None else (Qm - Qp - st.P if st.v == FIXED else -st.v * one[0])
            Mr = 0.0 if st.r is None else (Mm - Mp - st.Mo if st.r == FIXED else -st.r * one[1])
            self.reactions[float(x)] = (float(R), float(s * Mr))

    def _segment(self, x, side):
        br = self.breaks
        L = self.beam.L
        if not -1e-12 * L <= x <= L * (1 + 1e-12):
            raise OracleError("position outside the member")
        idx = np.where(np.abs(br - x) <= 1e-12 * L)[0]
        if idx.size:
            k = int(idx[0]) - 1 if side < 0 else int(idx[0])
            k = min(max(k, 0), len(br) - 2)
        else:
            k = int(np.searchsorted(br, x)) - 1
        return k

    def _state(self, x, side=-1):
        k = self._segment(x, side)
        xi = min(max((x - self.breaks[k]) / (self.breaks[k + 1] - self.breaks[k]), 0.0), 1.0)
        return self._sol.sol(xi)[4 * k:4 * k + 4] * self._sc

    # physical results with the signs of the chosen plane -------------------------------------
    def v(self, x, side=-1):
        """Deflection along the transverse axis."""
        return float(self._state(x, side)[0])

    def theta(self, x, side=-1):
        """Rotation of the cross-section (about +z for "xy", about +y for "xz")."""
        return float(self.beam._s * self._state(x, side)[1])

    def M(self, x, side=-1):
        """Bending moment (Mz or My), on the +x face; ``side`` = -1 just left of x, +1 just right."""
        return float(self.beam._s * self._state(x, side)[2])

    def Q(self, x, side=-1):
        """Transverse force along the UNDEFORMED transverse axis on the +x face (= shear force in a
        first-order analysis)."""
        return float(self._state(x, side)[3])

    def V(self, x, side=-1):
        """Shear force carried by the section (enters the shear strain): Q in a first-order analysis,
        Q - N v' ("engesser", also used when there is no shear deformation) or Q - N th ("haringx")."""
        k = self._segment(x, side)
        v, th, M, Q = self._state(x, side)
        if self.beam.shear_model == "haringx":
            return float(Q - self._N[k] * th)
        dv = (th + self._c[k] * Q) / (1.0 + self._c[k] * self._N[k])
        return float(Q - self._N[k] * dv)


# ====================================================================================== axial / torsion
class Bar(_Line):
    """Axial bar (``EA``) or shaft in uniform torsion (pass ``GJ`` as ``EA``):
    ``u' = N/EA, N' = -p``.  Results: ``u(x)``, ``N(x, side)`` (tension / torque on the +x face),
    ``reactions[x]`` = force by the support on the bar."""

    def __init__(self, length, EA):
        super().__init__(length)
        self.EA = _piecewise(EA, self.L, "EA")

    def support(self, x, u=FIXED, settlement=0.0):
        st = self._station(x)
        st.v, st.v0 = u, float(settlement)
        return self

    def point(self, x, P):
        self._station(x).P += float(P)
        return self

    def solve(self, tol=1e-10, max_nodes=200000) -> "BarSolution":
        L = self.L
        br = self._breaks([a for a, _, _ in self.EA])
        m = len(br) - 1
        h = np.diff(br)
        EA = np.array([_value(self.EA, x) for x in 0.5 * (br[:-1] + br[1:])])
        sc = np.array([L, float(EA.max())])
        stations = [self._at(x) for x in br]

        def rhs(xi, Y):
            out = np.empty_like(Y)
            for k in range(m):
                out[2 * k] = h[k] * Y[2 * k + 1] * sc[1] / EA[k] / sc[0]
                out[2 * k + 1] = -h[k] * self._load(br[k] + xi * h[k], k, br) / sc[1]
            return out

        def joint(st, left, right):
            res = []
            Nm = left[1] if left is not None else 0.0
            Np = right[1] if right is not None else 0.0
            u = (right if left is None else left)[0]
            if left is not None and right is not None:
                res.append((right[0] - left[0]) / sc[0])
            if st.v == FIXED:
                res.append((u - st.v0) / sc[0])
            else:
                res.append((Np - Nm + st.P - (st.v or 0.0) * u) / sc[1])
            return res

        def bc(Ya, Yb):
            ya = [Ya[2 * k:2 * k + 2] * sc for k in range(m)]
            yb = [Yb[2 * k:2 * k + 2] * sc for k in range(m)]
            res = joint(stations[0], None, ya[0])
            for k in range(1, m):
                res += joint(stations[k], yb[k - 1], ya[k])
            res += joint(stations[m], yb[m - 1], None)
            return np.array(res)

        xi = np.linspace(0.0, 1.0, 5)
        sol = solve_bvp(rhs, bc, xi, np.zeros((2 * m, xi.size)), tol=tol, max_nodes=max_nodes)
        if not sol.success or float(np.max(sol.rms_residuals)) > 10.0 * tol:
            raise OracleError(f"integration did not converge: {sol.message}")
        return BarSolution(self, br, sol, sc, stations, float(np.max(sol.rms_residuals)))


@dataclass
class BarSolution:
    bar: Bar
    breaks: np.ndarray
    _sol: object
    _sc: np.ndarray
    _stations: list
    rms_residual: float
    reactions: dict = field(init=False)

    def __post_init__(self):
        self.reactions = {}
        for x, st in zip(self.breaks, self._stations):
            if st.v is None:
                continue
            Nm = self.N(x, -1) if x > 0.0 else 0.0
            Np = self.N(x, +1) if x < self.bar.L else 0.0
            self.reactions[float(x)] = float(Nm - Np - st.P if st.v == FIXED else -st.v * self.u(x))

    def _state(self, x, side):
        br = self.breaks
        idx = np.where(np.abs(br - x) <= 1e-12 * self.bar.L)[0]
        if idx.size:
            k = min(max(int(idx[0]) - 1 if side < 0 else int(idx[0]), 0), len(br) - 2)
        else:
            k = int(np.searchsorted(br, x)) - 1
        xi = min(max((x - br[k]) / (br[k + 1] - br[k]), 0.0), 1.0)
        return self._sol.sol(xi)[2 * k:2 * k + 2] * self._sc

    def u(self, x, side=-1):
        return float(self._state(x, side)[0])

    def N(self, x, side=-1):
        return float(self._state(x, side)[1])
