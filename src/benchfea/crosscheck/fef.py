"""Fixed-end forces of a prismatic member under span loads, by the force (flexibility) method.

Why this exists: OpenSees 3.8 has no member point moments and ignores the end intensities of a
linearly varying load in 3D (measured, docs/CROSSCHECK.md).  The OpenSees translator therefore turns
such loads into EXACT equivalent nodal loads: with the element ends held, the loads produce the end
forces ``f`` computed here; applying ``-f`` at the element nodes and adding ``f`` to the element end
forces of that nodal-load analysis reproduces the loaded member exactly (superposition of the
"ends held" state and the "ends moving" state, both linear).

Derivation (x-y plane; v, th_z, Mz, Vy with the signs of docs/CONVENTIONS.md section 4)
----------------------------------------------------------------------------------------
Primary structure: cantilever clamped at end i (x = 0), free at end j (x = L).

1. Statics of the cantilever.  The stress resultants on the +x face at x equal the resultant of
   everything applied at xi > x:  ``V0(x) = sum F(xi)``,  ``M0(x) = sum (xi - x) F(xi) + sum Mo(xi)``
   (forces F along +y, moments Mo about +z).
2. Displacements of the free end by the unit-load theorem (virtual work with the internal forces of a
   unit load at j: unit force -> ``m = L - x, q = 1``; unit moment -> ``m = 1, q = 0``), with the
   curvature ``M/EI`` and the shear angle ``V/GAs``:
        v0  = int M0 (L - x)/EI dx + int V0/GAs dx          th0 = int M0/EI dx
   and the flexibility coefficients from the same integrals of the unit states:
        f_vv = L^3/(3EI) + L/GAs,   f_vt = L^2/(2EI),   f_tt = L/EI.
   A rotational connection spring k_i between the member end and joint i is a concentrated rotation
   ``M(0+)/k_i`` at x = 0 (the same sign rule as the curvature M/EI), a spring k_j one at x = L; in the
   unit-load integrals they add  ``v0 += M0(0+) L/k_i``, ``th0 += M0(0+)/k_i``,
   ``f_vv += L^2/k_i``, ``f_vt += L/k_i``, ``f_tt += 1/k_i + 1/k_j``.
3. Compatibility for the unknown end forces (Xv, Xt) that the joint exerts on the member at j:
        j connected:      v0 + f_vv Xv + f_vt Xt = 0,   th0 + f_vt Xv + f_tt Xt = 0
        j moment-free:    Xt = 0,                        v0 + f_vv Xv = 0
4. End i from equilibrium of the whole member:
        F_i = -(Xv + sum F),      M_i = -(Xt + L Xv + sum (xi F + Mo)).
   i moment-free, j clamped: the mirror image (x -> L - x; forces keep their sign, moments about z
   change sign).  Both ends moment-free: statics alone (``Xt = 0, M_i = 0``).

The integrands are piecewise polynomials (degree <= 4 for a linearly varying load), integrated exactly
by 4-point Gauss-Legendre quadrature between consecutive load positions.

The x-z plane (w, th_y, My, Vz) is the same system after th = -th_y, M = -My (CONVENTIONS section 4),
so point moments about y enter with the opposite sign and the resulting end moments are negated.

Axial and torsional loads (bar clamped at both ends, redundant X at j): ``N0(x) = sum P(xi > x)``,
``u0 = int N0/EA dx``, ``X = -u0 EA/L = -(1/L) int N0 dx``, ``P_i = -(X + sum P)``.

Verified in tests/test_crosscheck_fef.py against direct integration of the ODE (benchfea.oracle).
"""
from __future__ import annotations

import math

import numpy as np

__all__ = ["fixed_end_forces", "clip_loads", "plane_fef", "bar_fef", "bending_split"]

_GX, _GW = np.polynomial.legendre.leggauss(4)


def _gauss(f, a, b):
    """Integral of f over [a, b] (exact for polynomials of degree <= 7)."""
    h = 0.5 * (b - a)
    return h * sum(w * f(0.5 * (a + b) + h * x) for x, w in zip(_GX, _GW))


def _right_resultant(points, dists, x):
    """(sum F, sum (xi - x) F + sum Mo) of everything applied at xi > x.
    points: [(xi, F, Mo)], dists: [(x1, x2, w1, w2)]."""
    V = M = 0.0
    for xi, F, Mo in points:
        if xi > x:
            V += F
            M += (xi - x) * F + Mo
    for x1, x2, w1, w2 in dists:
        a = max(x, x1)
        if a < x2:
            w = lambda s: w1 + (w2 - w1) * (s - x1) / (x2 - x1)
            V += _gauss(w, a, x2)
            M += _gauss(lambda s: (s - x) * w(s), a, x2)
    return V, M


def _total(points, dists):
    """(sum F, sum (xi F + Mo)) about x = 0 of all loads."""
    V = sum(F for _, F, _ in points)
    M = sum(xi * F + Mo for xi, F, Mo in points)
    for x1, x2, w1, w2 in dists:
        h = x2 - x1
        V += 0.5 * (w1 + w2) * h
        M += 0.5 * w1 * h * (x1 + h / 3.0) + 0.5 * w2 * h * (x1 + 2.0 * h / 3.0)
    return V, M


def plane_fef(L, EI, GAs, points, dists, k_i=math.inf, k_j=math.inf):
    """End forces (F_i, M_i, F_j, M_j) ON the member (x-y plane signs) with the joints at both ends
    held.  ``k_i``, ``k_j``: rotational stiffness of the connection between the member end and the
    joint (``inf`` rigid, ``0`` hinge, else a spring).  ``GAs = None``: no shear deformation.
    Loads exactly at an end act on the joint."""
    if k_i == 0.0 and k_j != 0.0:
        # mirror image: x' = L - x, z' = -z
        pm = [(L - xi, F, -Mo) for xi, F, Mo in points]
        dm = [(L - x2, L - x1, w2, w1) for x1, x2, w1, w2 in dists]
        Fi, Mi, Fj, Mj = plane_fef(L, EI, GAs, pm, dm, k_j, 0.0)
        return Fj, -Mj, Fi, -Mi
    sumF, sumM = _total(points, dists)
    if k_i == 0.0 and k_j == 0.0:
        Xv = -sumM / L
        return -(Xv + sumF), 0.0, Xv, 0.0
    c = 0.0 if GAs is None else 1.0 / GAs
    xs = sorted({0.0, L} | {xi for xi, _, _ in points} | {x for x1, x2, _, _ in dists for x in (x1, x2)})
    xs = [x for x in xs if 0.0 <= x <= L]
    v0 = th0 = 0.0
    for a, b in zip(xs[:-1], xs[1:]):
        if b <= a:
            continue
        # Gauss points lie strictly inside (a, b): they never coincide with a load position
        v0 += _gauss(lambda x: _right_resultant(points, dists, x)[1] * (L - x) / EI
                     + _right_resultant(points, dists, x)[0] * c, a, b)
        th0 += _gauss(lambda x: _right_resultant(points, dists, x)[1] / EI, a, b)
    fi = 0.0 if math.isinf(k_i) else 1.0 / k_i
    fj = 0.0 if math.isinf(k_j) or k_j == 0.0 else 1.0 / k_j
    # connection spring at i = concentrated rotation M0(0+)/k_i at x = 0 (lever arm L to the free end)
    M00 = sum(xi * F + Mo for xi, F, Mo in points if xi > 0.0) + _total([], dists)[1]
    v0 += fi * M00 * L
    th0 += fi * M00
    fvv, fvt, ftt = L ** 3 / (3.0 * EI) + L * c + fi * L * L, L ** 2 / (2.0 * EI) + fi * L, L / EI + fi + fj
    if k_j == 0.0:
        Xv, Xt = -v0 / fvv, 0.0
    else:
        Xv, Xt = np.linalg.solve([[fvv, fvt], [fvt, ftt]], [-v0, -th0])
    return -(Xv + sumF), -(Xt + L * Xv + sumM), Xv, Xt


def bending_split(L, EI, GAs, k_i, k_j):
    """Coefficients (a, b, c) with which the bending stiffness of a member with rotational end
    connections ``k_i``, ``k_j`` in one plane equals the sum of three members in parallel between the
    same joints:  ``a`` x (rigidly connected member with EI and GAs both multiplied by a)
    + ``b`` x (member hinged at j) + ``c`` x (member hinged at i), the hinged members having the
    bending stiffness EI' = EI / (1 + 3 EI / (GAs L^2)) and no shear deformation.

    Derivation.  With joint i held, the member is described by the 2x2 stiffness K = F^-1 between
    (v_j, th_j) and the end forces at j; the rest of its stiffness follows from equilibrium.  The
    flexibility is that of the cantilever (unit-load theorem, see module docstring) plus the two
    connection rotations:  F = F0 + (1/k_i) [[L^2, L], [L, 1]] + (1/k_j) [[0, 0], [0, 1]],
    F0 = [[L^3/3EI + L/GAs, L^2/2EI], [L^2/2EI, L/EI]].  The three candidates have
    K_a = F0^-1,  K_b = (1/f) [[1, 0], [0, 0]],  K_c = (1/f) [[1, -L], [-L, L^2]]  with
    f = L^3/3EI + L/GAs (a member hinged at one end is a cantilever of tip flexibility f seen from
    the other end; an Euler-Bernoulli member with EI' has exactly this f).  K is symmetric, so
    K = a K_a + b K_b + c K_c is three linear equations for (a, b, c).  For a single spring on an
    Euler-Bernoulli member the solution is a = 1/(1 + 4EI/(kL)), b = 1 - a.
    """
    if math.isinf(k_i) and math.isinf(k_j):
        return 1.0, 0.0, 0.0
    c = 0.0 if GAs is None else 1.0 / GAs
    f = L ** 3 / (3.0 * EI) + L * c
    if k_i == 0.0 and k_j == 0.0:
        return 0.0, 0.0, 0.0
    if k_i == 0.0:
        return 0.0, 0.0, (1.0 if math.isinf(k_j) else f / (f + L * L / k_j))
    if k_j == 0.0:
        return 0.0, (1.0 if math.isinf(k_i) else f / (f + L * L / k_i)), 0.0
    F0 = np.array([[f, L ** 2 / (2.0 * EI)], [L ** 2 / (2.0 * EI), L / EI]])
    fi = 0.0 if math.isinf(k_i) else 1.0 / k_i
    fj = 0.0 if math.isinf(k_j) else 1.0 / k_j
    K = np.linalg.inv(F0 + fi * np.array([[L * L, L], [L, 1.0]]) + fj * np.array([[0.0, 0.0], [0.0, 1.0]]))
    Ka = np.linalg.inv(F0)
    A = np.array([[Ka[0, 0], 1.0 / f, 1.0 / f], [Ka[0, 1] / L, 0.0, -1.0 / f], [Ka[1, 1] / L ** 2, 0.0, 1.0 / f]])
    a, b, cc = np.linalg.solve(A, [K[0, 0], K[0, 1] / L, K[1, 1] / L ** 2])
    return float(a), float(b), float(cc)


def bar_fef(L, points, dists):
    """End forces (P_i, P_j) ON a bar clamped at both ends; points: [(xi, P)], dists as above."""
    pts = [(xi, P, 0.0) for xi, P in points]
    xs = sorted({0.0, L} | {xi for xi, _ in points} | {x for x1, x2, _, _ in dists for x in (x1, x2)})
    integral = sum(_gauss(lambda x: _right_resultant(pts, dists, x)[0], a, b)
                   for a, b in zip(xs[:-1], xs[1:]) if b > a)
    X = -integral / L
    return -(X + _total(pts, dists)[0]), X


def clip_loads(points, dists, xa, xb, last, tol=0.0):
    """Loads of the sub-span [xa, xb] with positions measured from xa.  A point load at a sub-span
    boundary belongs to the sub-span that starts there (to the last sub-span if at the member end).
    ``points`` = [(kind, vector, x)], ``dists`` = [(w1, w2, x1, x2)] (prepare.PMember layout)."""
    pts = []
    for kind, vec, x in points:
        if (xa - tol <= x < xb - tol) or (last and x >= xb - tol):
            pts.append((kind, vec, min(max(x - xa, 0.0), xb - xa)))
    ds = []
    for w1, w2, x1, x2 in dists:
        a, b = max(x1, xa), min(x2, xb)
        if b - a > tol:
            wa = w1 + (w2 - w1) * (a - x1) / (x2 - x1)
            wb = w1 + (w2 - w1) * (b - x1) / (x2 - x1)
            ds.append((wa, wb, a - xa, b - xa))
    return pts, ds


def fixed_end_forces(L, EIy, EIz, GAy, GAz, points, dists, ky=(math.inf, math.inf), kz=(math.inf, math.inf)):
    """12 end forces ON the member (local axes, ``Fx, Fy, Fz, Mx, My, Mz`` at i then j) with both joints
    held.  ``kz = (k_i, k_j)``: rotational connection stiffness about local z at the two ends (bending in
    the x-y plane; ``inf`` rigid, ``0`` hinge), ``ky``: about local y (x-z plane).  ``points`` =
    [(kind, vector(3), x)], ``dists`` = [(w1(3), w2(3), x1, x2)] in local axes."""
    f = np.zeros(12)
    force = [(x, v) for kind, v, x in points if kind == "force"]
    moment = [(x, v) for kind, v, x in points if kind == "moment"]
    comp = lambda k: [(x1, x2, w1[k], w2[k]) for w1, w2, x1, x2 in dists if w1[k] != 0.0 or w2[k] != 0.0]
    # axial
    px = [(x, v[0]) for x, v in force if v[0] != 0.0]
    if px or comp(0):
        f[0], f[6] = bar_fef(L, px, comp(0))
    # torsion (no distributed torque in the model)
    tx = [(x, v[0]) for x, v in moment if v[0] != 0.0]
    if tx:
        f[3], f[9] = bar_fef(L, tx, [])
    # x-y plane: forces along y, moments about z
    py = [(x, v[1], 0.0) for x, v in force if v[1] != 0.0] + [(x, 0.0, v[2]) for x, v in moment if v[2] != 0.0]
    if py or comp(1):
        f[1], f[5], f[7], f[11] = plane_fef(L, EIz, GAy, py, comp(1), *kz)
    # x-z plane: forces along z, moments about y enter and leave with the opposite sign
    pz = [(x, v[2], 0.0) for x, v in force if v[2] != 0.0] + [(x, 0.0, -v[1]) for x, v in moment if v[1] != 0.0]
    if pz or comp(2):
        Fi, Mi, Fj, Mj = plane_fef(L, EIy, GAz, pz, comp(2), *ky)
        f[2], f[4], f[8], f[10] = Fi, -Mi, Fj, -Mj
    return f
