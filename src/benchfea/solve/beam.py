"""Exact integration of the equations of a straight prismatic member in its local axes.

Sign conventions are those of docs/CONVENTIONS.md section 4; the derivation is in docs/THEORY.md
section 1.  With the section resultants s = (N, Vy, Vz, T, My, Mz) on the +x face and the axis
displacements y = (u, v, w, theta_x, theta_y, theta_z):

    equilibrium       N' = -w_x      Vy' = -w_y     Vz' = -w_z     T' = 0
                      Mz' = -Vy      My' = +Vz
    kinematics +      u' = N / EA              theta_x' = T / GJ
    elasticity        theta_z' = Mz / EIz      v' = theta_z + Vy / (G Asy)
                      theta_y' = My / EIy      w' = -theta_y + Vz / (G Asz)

Loads are point forces/moments and piecewise linear line loads, so between two load break points the
resultants are polynomials of degree <= 3 and the displacements polynomials of degree <= 5: they are
integrated in closed form, segment by segment (no quadrature, no shape-function approximation).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["SectionFlex", "SpanLoads", "integrate", "Integration"]


@dataclass(frozen=True)
class SectionFlex:
    """Flexibilities per unit length: 1/EA, 1/GJ, 1/EIy, 1/EIz, 1/(G Asy), 1/(G Asz) (0 = no shear deformation)."""
    axial: float
    torsion: float
    bend_y: float
    bend_z: float
    shear_y: float = 0.0
    shear_z: float = 0.0


@dataclass(frozen=True)
class SpanLoads:
    """Loads on the flexible length in member LOCAL axes.

    points  tuple of (x, f) with f = (Fx, Fy, Fz, Mx, My, Mz) applied at distance x from the start
    dists   tuple of (x1, x2, w1, w2): line load varying linearly from w1(3) at x1 to w2(3) at x2
    """
    points: tuple = ()
    dists: tuple = ()

    def __bool__(self) -> bool:
        return bool(self.points or self.dists)

    @staticmethod
    def combine(terms) -> "SpanLoads":
        """Sum of ``factor * loads`` over ``terms`` = iterable of (factor, SpanLoads)."""
        pts, dst = [], []
        for f, ld in terms:
            if f == 0.0:
                continue
            pts += [(x, f * v) for x, v in ld.points]
            dst += [(x1, x2, f * w1, f * w2) for x1, x2, w1, w2 in ld.dists]
        return SpanLoads(tuple(pts), tuple(dst))

    def resultant(self) -> np.ndarray:
        """Force and moment (6,) of all loads, reduced to the start of the flexible length (x = 0)."""
        r = np.zeros(6)
        for x, f in self.points:
            r[:3] += f[:3]
            r[3:] += f[3:] + np.cross([x, 0.0, 0.0], f[:3])
        for x1, x2, w1, w2 in self.dists:
            h = x2 - x1
            # a linearly varying load = two triangles: w1 h/2 at x1 + h/3 and w2 h/2 at x1 + 2h/3
            for w, xc in ((w1, x1 + h / 3.0), (w2, x1 + 2.0 * h / 3.0)):
                p = w * (h / 2.0)
                r[:3] += p
                r[3:] += np.cross([xc, 0.0, 0.0], p)
        return r

    def breakpoints(self) -> list[float]:
        return [x for x, _ in self.points] + [x for d in self.dists for x in d[:2]]


@dataclass
class Integration:
    breaks: np.ndarray        # (m+1,) segment boundaries
    force_poly: np.ndarray    # (m, 6, 4) ascending coefficients in xi = x - breaks[k]
    disp_poly: np.ndarray     # (m, 6, 6)
    s_start: np.ndarray       # resultants at x = 0-   (before point loads at x = 0)
    y_start: np.ndarray
    s_end: np.ndarray         # resultants at x = L+   (after point loads at x = L)
    y_end: np.ndarray

    def at(self, x: float, side: int = -1):
        """(resultants, displacements) at x; ``side`` -1 = limit from the left, +1 = from the right
        (they differ only at a point load)."""
        b = self.breaks
        k = int(np.searchsorted(b, x, side="left"))
        if k < len(b) and b[k] == x:                       # on a break point
            if side < 0:
                if k == 0:
                    return self.s_start.copy(), self.y_start.copy()
                k, xi = k - 1, b[k] - b[k - 1]
            else:
                if k == len(b) - 1:
                    return self.s_end.copy(), self.y_end.copy()
                xi = 0.0
        else:
            k = min(max(k - 1, 0), len(b) - 2)
            xi = x - b[k]
        return _polyval(self.force_poly[k], xi), _polyval(self.disp_poly[k], xi)


def _polyval(c: np.ndarray, xi: float) -> np.ndarray:
    out = c[:, -1].copy()
    for k in range(c.shape[1] - 2, -1, -1):
        out = out * xi + c[:, k]
    return out


def _integral(c: np.ndarray) -> np.ndarray:
    """Coefficients of the integral from 0 of the polynomial with ascending coefficients c."""
    return np.concatenate([[0.0], c / np.arange(1, len(c) + 1)])


def integrate(length: float, flex: SectionFlex, loads: SpanLoads, s_start, y_start) -> Integration:
    """Integrate the member equations from x = 0 to x = L.

    ``s_start`` = resultants (N, Vy, Vz, T, My, Mz) just before x = 0 (i.e. minus the forces acting on
    the member at its start), ``y_start`` = (u, v, w, theta_x, theta_y, theta_z) at x = 0.
    """
    L = float(length)
    pts = sorted({0.0, L} | {min(max(float(x), 0.0), L) for x in loads.breakpoints()})
    breaks = np.array(pts)
    jumps = {}
    for x, f in loads.points:
        x = min(max(float(x), 0.0), L)
        jumps[x] = jumps.get(x, 0.0) + np.asarray(f, float)
    s = np.array(s_start, float)
    y = np.array(y_start, float)
    s0, y0 = s.copy(), y.copy()
    m = len(breaks) - 1
    FP, DP = np.zeros((m, 6, 4)), np.zeros((m, 6, 6))
    if 0.0 in jumps:
        s = s - jumps[0.0]                 # a load on the member at x: the resultant drops by the load
    for k in range(m):
        a, b = breaks[k], breaks[k + 1]
        wa, ws = np.zeros(3), np.zeros(3)
        for x1, x2, w1, w2 in loads.dists:
            if x1 <= a and b <= x2:
                slope = (w2 - w1) / (x2 - x1)
                wa += w1 + slope * (a - x1)
                ws += slope
        F = FP[k]
        for c in range(3):                                   # N' = -wx, Vy' = -wy, Vz' = -wz
            F[c, :3] = (s[c], -wa[c], -ws[c] / 2.0)
        F[3, 0] = s[3]                                       # T' = 0
        F[4] = (s[4], s[2], -wa[2] / 2.0, -ws[2] / 6.0)      # My' = +Vz
        F[5] = (s[5], -s[1], wa[1] / 2.0, ws[1] / 6.0)       # Mz' = -Vy
        D = DP[k]
        D[0, :5] = flex.axial * _integral(F[0])              # u' = N / EA
        D[3, :5] = flex.torsion * _integral(F[3])            # theta_x' = T / GJ
        D[4, :5] = flex.bend_y * _integral(F[4])             # theta_y' = My / EIy
        D[5, :5] = flex.bend_z * _integral(F[5])             # theta_z' = Mz / EIz
        D[[0, 3, 4, 5], 0] = y[[0, 3, 4, 5]]
        D[1] = _integral(D[5, :5] + np.append(flex.shear_y * F[1], 0.0))     # v' = theta_z + Vy / GAsy
        D[2] = _integral(-D[4, :5] + np.append(flex.shear_z * F[2], 0.0))    # w' = -theta_y + Vz / GAsz
        D[1, 0], D[2, 0] = y[1], y[2]
        h = b - a
        s, y = _polyval(F, h), _polyval(D, h)
        if b in jumps:
            s = s - jumps[b]
    return Integration(breaks, FP, DP, s0, y0, s, y)
