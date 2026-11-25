#!/usr/bin/env python3
"""Reference solutions for the vendor verification problems (verification/problems/staad, .../ansys),
"recomputed from the published problem statement".

The vendor pages cite textbooks that were not opened for this project, so no textbook number is used:
every reference value is computed HERE from the problem statement on the vendor page, by means that
share no code with ``benchfea.solve``:

* single members and continuous beams, first and second order: the independent ODE oracle
  ``benchfea.oracle`` (direct integration of the member differential equations, scipy ``solve_bvp``);
* Euler load of a pin-ended column: smallest root of the characteristic equation sin(kL) = 0 of
  EI v'''' + P v'' = 0 with v = v'' = 0 at both ends (general solution v = A sin kx + B cos kx + C x + D;
  the end conditions give B = D = 0, C = 0, A sin kL = 0), found numerically;
* plane frames and trusses: ``PlaneFrame`` below, a direct stiffness program written for this purpose
  (closed-form 6 x 6 Euler-Bernoulli member matrix with axial deformation, no shear deformation; uniform
  member loads through fixed-end forces; moment releases by static condensation; prescribed support
  displacements);
* space frames and space trusses under joint loads: ``SpaceFrame`` below (closed-form 12 x 12 member
  matrix; pin-ended bars carry axial force only).

``python3 verification/derive/vendor_reference.py`` prints every reference value used by the builders
(``reference(name)`` returns them as a dict).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from benchfea.oracle import Bar, Beam  # noqa: E402


# =========================================================================================== plane
class PlaneFrame:
    """Plane frame / truss, DOFs (ux, uy, rz) per node.  Member end forces are returned as the forces by
    the joints ON the member in member axes (x from i to j, y = z x X, i.e. x rotated by +90 degrees):
    [Fx_i, Fy_i, M_i, Fx_j, Fy_j, M_j]."""

    def __init__(self):
        self.nodes, self.members, self.fixed, self.F, self.q = {}, {}, {}, {}, {}

    def node(self, name, x, y):
        self.nodes[name] = (float(x), float(y))
        return self

    def member(self, name, i, j, EA, EI, hinge_i=False, hinge_j=False, bar=False):
        self.members[name] = (i, j, float(EA), float(EI), hinge_i or bar, hinge_j or bar)
        return self

    def support(self, node, dofs, values=None):
        for k, d in enumerate(dofs):
            self.fixed[(node, "xyr".index(d))] = 0.0 if values is None else float(values[k])
        return self

    def load(self, node, fx=0.0, fy=0.0, m=0.0):
        cur = self.F.get(node, (0.0, 0.0, 0.0))
        self.F[node] = (cur[0] + fx, cur[1] + fy, cur[2] + m)
        return self

    def uniform(self, member, gx=0.0, gy=0.0):
        """Uniform load per unit length of the member, global components."""
        self.q[member] = (float(gx), float(gy))
        return self

    def _element(self, name):
        i, j, EA, EI, hi, hj = self.members[name]
        (xi, yi), (xj, yj) = self.nodes[i], self.nodes[j]
        L = math.hypot(xj - xi, yj - yi)
        c, s = (xj - xi) / L, (yj - yi) / L
        a, b = EA / L, EI / L ** 3
        k = np.array([[a, 0, 0, -a, 0, 0],
                      [0, 12 * b, 6 * b * L, 0, -12 * b, 6 * b * L],
                      [0, 6 * b * L, 4 * b * L * L, 0, -6 * b * L, 2 * b * L * L],
                      [-a, 0, 0, a, 0, 0],
                      [0, -12 * b, -6 * b * L, 0, 12 * b, -6 * b * L],
                      [0, 6 * b * L, 2 * b * L * L, 0, -6 * b * L, 4 * b * L * L]])
        gx, gy = self.q.get(name, (0.0, 0.0))
        qx, qy = c * gx + s * gy, -s * gx + c * gy
        fef = np.array([-qx * L / 2, -qy * L / 2, -qy * L * L / 12, -qx * L / 2, -qy * L / 2, qy * L * L / 12])
        rel = [d for d, h in ((2, hi), (5, hj)) if h]
        if rel:                                             # static condensation of the released rotations
            keep = [d for d in range(6) if d not in rel]
            krr = k[np.ix_(rel, rel)]
            kc = np.zeros((6, 6))
            fc = np.zeros(6)
            kc[np.ix_(keep, keep)] = k[np.ix_(keep, keep)] - k[np.ix_(keep, rel)] @ np.linalg.solve(krr, k[np.ix_(rel, keep)])
            fc[keep] = fef[keep] - k[np.ix_(keep, rel)] @ np.linalg.solve(krr, fef[rel])
            k, fef = kc, fc
        T = np.zeros((6, 6))
        T[:3, :3] = T[3:, 3:] = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
        return k, fef, T, L

    def solve(self):
        names = list(self.nodes)
        idx = {n: 3 * k for k, n in enumerate(names)}
        n = 3 * len(names)
        K, F = np.zeros((n, n)), np.zeros(n)
        parts = {}
        for name, (i, j, *_rest) in self.members.items():
            k, fef, T, L = self._element(name)
            d = [idx[i], idx[i] + 1, idx[i] + 2, idx[j], idx[j] + 1, idx[j] + 2]
            K[np.ix_(d, d)] += T.T @ k @ T
            F[d] -= T.T @ fef
            parts[name] = (k, fef, T, d)
        for node, f in self.F.items():
            F[idx[node]:idx[node] + 3] += f
        u = np.zeros(n)
        fixed = {idx[node] + d: v for (node, d), v in self.fixed.items()}
        unused = [d for d in range(n) if K[d, d] == 0.0 and d not in fixed]      # rotations of pure bar joints
        free = [d for d in range(n) if d not in fixed and d not in unused]
        fd = list(fixed)
        u[fd] = [fixed[d] for d in fd]
        u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free] - K[np.ix_(free, fd)] @ u[fd])
        R = K @ u - F
        out = {"u": {nm: u[idx[nm]:idx[nm] + 3].copy() for nm in names},
               "R": {nm: R[idx[nm]:idx[nm] + 3].copy() for nm in names if any((nm, d) in self.fixed for d in range(3))},
               "end": {name: k @ (T @ u[d]) + fef for name, (k, fef, T, d) in parts.items()}}
        return out


# =========================================================================================== space
class SpaceFrame:
    """Space frame / truss under joint loads, DOFs (ux, uy, uz, rx, ry, rz).  ``member`` = frame member with
    EA, GJ, EIy, EIz (local y, z: any orthonormal pair -- only used with EIy = EIz); ``bar`` = pin-ended bar."""

    def __init__(self):
        self.nodes, self.members, self.fixed, self.F = {}, {}, set(), {}

    def node(self, name, x, y, z):
        self.nodes[name] = np.array([x, y, z], float)
        return self

    def member(self, name, i, j, EA, GJ, EIy, EIz):
        self.members[name] = (i, j, EA, GJ, EIy, EIz)
        return self

    def bar(self, name, i, j, EA):
        self.members[name] = (i, j, EA, 0.0, 0.0, 0.0)
        return self

    def support(self, node, dofs):
        for d in dofs:
            self.fixed.add((node, d))
        return self

    def load(self, node, f):
        self.F[node] = np.array(f, float)
        return self

    def solve(self):
        names = list(self.nodes)
        idx = {n: 6 * k for k, n in enumerate(names)}
        n = 6 * len(names)
        K, F = np.zeros((n, n)), np.zeros(n)
        parts = {}
        for name, (i, j, EA, GJ, EIy, EIz) in self.members.items():
            dx = self.nodes[j] - self.nodes[i]
            L = float(np.linalg.norm(dx))
            ex = dx / L
            ref = np.array([0.0, 0.0, 1.0]) if abs(ex[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
            ey = np.cross(ref, ex)
            ey /= np.linalg.norm(ey)
            ez = np.cross(ex, ey)
            R = np.array([ex, ey, ez])
            k = np.zeros((12, 12))
            for a, b, v in ((0, 6, EA / L), (3, 9, GJ / L)):
                k[a, a] = k[b, b] = v
                k[a, b] = k[b, a] = -v
            for (v1, t1, v2, t2, EI, sg) in ((1, 5, 7, 11, EIz, 1.0), (2, 4, 8, 10, EIy, -1.0)):
                c = EI / L ** 3
                blk = c * np.array([[12, sg * 6 * L, -12, sg * 6 * L], [sg * 6 * L, 4 * L * L, -sg * 6 * L, 2 * L * L],
                                    [-12, -sg * 6 * L, 12, -sg * 6 * L], [sg * 6 * L, 2 * L * L, -sg * 6 * L, 4 * L * L]])
                k[np.ix_([v1, t1, v2, t2], [v1, t1, v2, t2])] += blk
            T = np.zeros((12, 12))
            for q in range(4):
                T[3 * q:3 * q + 3, 3 * q:3 * q + 3] = R
            d = list(range(idx[i], idx[i] + 6)) + list(range(idx[j], idx[j] + 6))
            K[np.ix_(d, d)] += T.T @ k @ T
            parts[name] = (k, T, d)
        for node, f in self.F.items():
            F[idx[node]:idx[node] + len(f)] += f
        fixed = {idx[node] + "xyzabc".index(d) for node, d in self.fixed}
        unused = [d for d in range(n) if K[d, d] == 0.0 and d not in fixed]
        free = [d for d in range(n) if d not in fixed and d not in unused]
        u = np.zeros(n)
        u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        R = K @ u - F
        return {"u": {nm: u[idx[nm]:idx[nm] + 6].copy() for nm in names},
                "R": {nm: R[idx[nm]:idx[nm] + 6].copy() for nm in names if any((nm, d) in self.fixed for d in "xyzabc")},
                # axial force, tension positive
                "N": {name: float((k @ (T @ u[d]))[6]) for name, (k, T, d) in parts.items()},
                "end": {name: k @ (T @ u[d]) for name, (k, T, d) in parts.items()}}


def euler_pinned(EI, L):
    """Smallest root of sin(kL) = 0 (k > 0) -> P = k^2 EI."""
    k = brentq(lambda kk: math.sin(kk * L), 2.0 / L, 4.0 / L, xtol=1e-15 / L, rtol=1e-15)
    return k * k * EI


# =========================================================================================== cases
REFERENCES = {}


def case(fn):
    REFERENCES[fn.__name__] = fn
    return fn


def reference(name):
    return REFERENCES[name]()


# ------------------------------------------------------------------------------------------ STAAD
@case
def staad_axially_loaded_column():
    """Bar of 10 in held axially at both ends; -1000 lb at 7 in, -500 lb at 4 in (x upward). EA = 30e6 x 1."""
    s = Bar(10.0, 30.0e6).support(0.0).support(10.0).point(7.0, -1000.0).point(4.0, -500.0).solve()
    return {"R_bottom": -s.N(0.0, +1), "R_top": s.N(10.0, -1), "u4": s.u(4.0, -1), "u7": s.u(7.0, -1)}


@case
def staad_beam_trapezoid():
    """7.5 m, pinned at 0, fixed at 7.5; load 4 -> 7 kN/m (down) from 3 m to 7.5 m. EI = 2e8 x 5e-5 kN m^2."""
    L = 7.5
    s = (Beam(L, 2.0e8 * 5.0e-5).support(0.0, v="fixed").support(L, v="fixed", r="fixed")
         .dist(3.0, L, -4.0, -7.0).solve())
    return {"thetaA": s.theta(0.0, +1), "RA": s.reactions[0.0][0], "RB": s.reactions[L][0],
            "MB_reaction": s.reactions[L][1], "v3": s.v(3.0), "M_B_internal": s.M(L, -1)}


@case
def staad_propped_cantilever_selfweight():
    """100 in, pinned at 0, fixed at 100; w = 0.282072 lb/in^3 x 4 in^2 (down; the density of the page's input
    file, which its text rounds to 0.2821); EI = 30e6 x 1.3333."""
    L, w = 100.0, 0.282072 * 4.0
    s = Beam(L, 30.0e6 * 1.3333).support(0.0, v="fixed").support(L, v="fixed", r="fixed").dist(0.0, L, -w).solve()
    return {"R1": s.reactions[0.0][0], "R3": s.reactions[L][0], "M3_reaction": s.reactions[L][1],
            "theta1": s.theta(0.0, +1), "v_4215": s.v(42.15), "M_4215": s.M(42.15), "M3_internal": s.M(L, -1)}


@case
def staad_stepped_shaft():
    """Three 500 mm segments, solid circular d = 80, 60, 40 mm, G = 80 kN/mm^2; torques 3000, 2000, 800
    kN mm at the ends of the segments (units of the page's input file).  Segment torques by statics."""
    G = 80.0
    twist, out = 0.0, {}
    for k, (d, t) in enumerate(((80.0, 5800.0), (60.0, 2800.0), (40.0, 800.0))):
        J = math.pi * d ** 4 / 32.0
        s = Bar(500.0, G * J).support(0.0).point(500.0, t).solve()
        twist += s.u(500.0, -1)
        out[f"twist{k + 2}"] = twist
        out[f"J{k + 1}"] = J
    return out


def _portal1():
    f = PlaneFrame()
    for n, (x, y) in {"n1": (0, 0), "n2": (0, 100), "n3": (60, 100), "n4": (120, 100), "n5": (120, 0)}.items():
        f.node(n, x, y)
    E = 30.0e6
    f.member("c1", "n1", "n2", E * 4.0, E * 2.0 ** 4 / 12.0).member("c2", "n5", "n4", E * 4.0, E * 2.0 ** 4 / 12.0)
    f.member("b1", "n2", "n3", E * 8.0, E * 2.0 * 4.0 ** 3 / 12.0).member("b2", "n3", "n4", E * 8.0, E * 2.0 * 4.0 ** 3 / 12.0)
    f.support("n1", "xy").support("n5", "xy").load("n3", fy=-1000.0)
    return f.solve()


@case
def staad_portal_reactions_1():
    return _portal1()


@case
def staad_portal_reactions_2():
    """20 ft x 20 ft portal, fixed bases, A = 1 ft^2, I = 0.05 ft^4, E = 4.132e6 ksf, 2 kip/ft on the beam."""
    f = PlaneFrame()
    for n, (x, y) in {"n1": (0, 0), "n2": (0, 20), "m": (10, 20), "n3": (20, 20), "n4": (20, 0)}.items():
        f.node(n, x, y)
    E = 4.132e6
    f.member("c1", "n1", "n2", E, E * 0.05).member("c2", "n4", "n3", E, E * 0.05)
    f.member("b1", "n2", "m", E, E * 0.05).member("b2", "m", "n3", E, E * 0.05)
    f.uniform("b1", gy=-2.0).uniform("b2", gy=-2.0).support("n1", "xyr").support("n4", "xyr")
    return f.solve()


@case
def staad_frame_1x2_lateral():
    f = PlaneFrame()
    for n, (x, y) in {"n1": (0, 0), "n2": (20, 0), "n3": (0, 20), "n4": (20, 20), "n5": (0, 40), "n6": (20, 40)}.items():
        f.node(n, x, y)
    E = 4.176e6
    for name, i, j in (("m1", "n1", "n3"), ("m2", "n2", "n4"), ("m3", "n3", "n5"), ("m4", "n4", "n6"),
                       ("m5", "n3", "n4"), ("m6", "n5", "n6")):
        f.member(name, i, j, E * 0.2, E * 0.1)
    f.support("n1", "xyr").support("n2", "xyr").load("n3", fx=20.0).load("n5", fx=10.0)
    return f.solve()


def _two_bay(case_no):
    f = PlaneFrame()
    for n, (x, y) in {"n1": (0, 0), "n2": (216, 0), "n3": (432, 0), "n4": (0, 120), "n5": (216, 120),
                      "n6": (432, 120), "ma": (108, 120), "mb": (324, 120)}.items():
        f.node(n, x, y)
    E, A = 30000.0, 1.0e7
    for name, i, j in (("c1", "n1", "n4"), ("c2", "n2", "n5"), ("c3", "n3", "n6")):
        f.member(name, i, j, E * A, E * 13824.0)
    f.member("a1", "n4", "ma", E * A, E * 27000.0, hinge_i=True).member("a2", "ma", "n5", E * A, E * 27000.0)
    f.member("b1", "n5", "mb", E * A, E * 27000.0).member("b2", "mb", "n6", E * A, E * 27000.0, hinge_j=True)
    f.support("n1", "xyr").support("n3", "xyr").support("n2", "xy")
    if case_no == 1:
        f.load("n4", fy=-50.0).load("n6", fy=-50.0).load("n5", fy=-100.0).load("ma", fy=-100.0).load("mb", fy=-100.0)
    else:
        for name in ("a1", "a2", "b1", "b2"):
            f.uniform(name, gy=-10.0 / 12.0)
    return f.solve()


@case
def staad_two_bay_frame():
    return {"LC1": _two_bay(1), "LC2": _two_bay(2)}


def _sidesway(case_no):
    f = PlaneFrame()
    for n, (x, y) in {"n1": (0, 0), "n2": (0, 750), "n3": (600, 750), "n4": (600, 250)}.items():
        f.node(n, x, y)
    E = 2171.84
    f.member("m1", "n1", "n2", E * 2400.0, E * 720000.0).member("m2", "n2", "n3", E * 4800.0, E * 1.44e6)
    f.member("m3", "n3", "n4", E * 2400.0, E * 720000.0)
    f.support("n1", "xy", (0.0, -1.0 if case_no == 1 else 0.0)).support("n4", "y", (-1.0 if case_no == 2 else 0.0,))
    return f.solve()


@case
def staad_settlement_frame():
    return {"LC1": _sidesway(1), "LC2": _sidesway(2)}


@case
def staad_space_frame():
    f = SpaceFrame()
    f.node("n1", 0, 120, 0).node("n2", 240, 120, 0).node("n3", 0, 0, 0).node("n4", 360, 0, 120)
    E, G = 30000.0, 30000.0 / (2.0 * 1.25)
    for name, i, j in (("m1", "n1", "n2"), ("m2", "n3", "n1"), ("m3", "n2", "n4")):
        f.member(name, i, j, E * 11.0, G * 83.0, E * 56.0, E * 56.0)
    f.support("n3", "xyzabc").support("n4", "xyzabc")
    f.load("n1", [2.0, 0, 0, 0, 0, 0]).load("n2", [0, -1.0, 0, 0, 0, -120.0])
    return f.solve()


@case
def staad_simple_frame():
    f = PlaneFrame()
    for n, (x, y) in {"n1": (0, 0), "n2": (0, 10), "n3": (0, 20), "n4": (10, 20)}.items():
        f.node(n, x, y)
    E = 3000.0
    for name, i, j in (("m1", "n1", "n2"), ("m2", "n2", "n3"), ("m3", "n3", "n4")):
        f.member(name, i, j, E * 10.0, E * 100.0)
    f.support("n1", "xyr").support("n2", "x").load("n4", fy=-1.0)
    return f.solve()


def _cantilever_column(L, EI, P, H):
    s = Beam(L, EI, N=-P, second_order=True).support(0.0, v="fixed", r="fixed").point(L, P=H).solve()
    s1 = Beam(L, EI).support(0.0, v="fixed", r="fixed").point(L, P=H).solve()
    return {"M_base": s.M(0.0, +1), "tip": s.v(L), "M_base_first": s1.M(0.0, +1), "tip_first": s1.v(L)}


@case
def staad_pdelta_column():
    """12 ft cantilever, EI = 29,000 x 30.80, P = 3 kip, H = 5 kip."""
    return _cantilever_column(144.0, 29000.0 * 30.80, 3.0, 5.0)


@case
def staad_direct_analysis_column():
    """25 ft cantilever, EI = 29,000 x 61.9, P = 30 kip, H = 0.09 kip."""
    return _cantilever_column(300.0, 29000.0 * 61.9, 30.0, 0.09)


@case
def staad_column_buckling():
    return {"Pcr": euler_pinned(29730.0 * 248.0, 360.0), "factor": euler_pinned(29730.0 * 248.0, 360.0) / 160.0}


@case
def staad_direct_analysis_beam():
    """40 ft simply supported, w = 0.1 kip/ft, P = 100 kip compression, EI = 29,000 x 341."""
    L, w = 480.0, 0.1 / 12.0
    s = Beam(L, 29000.0 * 341.0, N=-100.0, second_order=True).support(0.0, v="fixed").support(L, v="fixed").dist(0.0, L, -w).solve()
    return {"M_mid": s.M(L / 2), "v_mid": s.v(L / 2)}


def _two_bar():
    a, b = 180.0 * math.cos(math.radians(30.0)), 180.0 * math.sin(math.radians(30.0))
    f = PlaneFrame()
    f.node("n1", 0, 0).node("n2", a, -b).node("n3", 2 * a, 0)
    f.member("m1", "n1", "n2", 30.0e6 * 0.5, 1.0, bar=True).member("m2", "n2", "n3", 30.0e6 * 0.5, 1.0, bar=True)
    f.support("n1", "xy").support("n3", "xy").load("n2", fy=-5000.0)
    return f.solve()


@case
def two_bar_truss():
    return _two_bar()


@case
def staad_truss_reactions_3():
    f = PlaneFrame()
    for n, (x, y) in {"n1": (0, 0), "n2": (150, 100), "n3": (150, 50), "n4": (300, 0)}.items():
        f.node(n, x, y)
    E = 30000.0
    for name, i, j, A in (("m1", "n1", "n2", 5.0), ("m2", "n1", "n3", 3.0), ("m3", "n2", "n3", 2.0),
                          ("m4", "n2", "n4", 5.0), ("m5", "n3", "n4", 3.0)):
        f.member(name, i, j, E * A, 1.0, bar=True)
    f.support("n1", "xy").support("n4", "xy").load("n2", fy=-10.0)
    return f.solve()


TRUSS20_NODES = {"n1": (0, 0), "n2": (180, 0), "n3": (360, 0), "n4": (540, 0), "n5": (720, 0), "n6": (180, 90),
                 "n7": (360, 180), "n8": (540, 90)}
TRUSS20_MEMBERS = (("m1", "n2", "n6", 1.0), ("m2", "n3", "n4", 1.0), ("m3", "n4", "n8", 1.0), ("m4", "n4", "n5", 1.0),
                   ("m5", "n1", "n2", 2.0), ("m6", "n2", "n3", 2.0), ("m7", "n3", "n6", 1.5), ("m8", "n3", "n8", 1.5),
                   ("m9", "n3", "n7", 3.0), ("m10", "n1", "n6", 3.0), ("m11", "n5", "n8", 3.0), ("m12", "n6", "n7", 4.0),
                   ("m13", "n7", "n8", 4.0))


@case
def staad_truss_deflection():
    f = PlaneFrame()
    for n, (x, y) in TRUSS20_NODES.items():
        f.node(n, x, y)
    for name, i, j, A in TRUSS20_MEMBERS:
        f.member(name, i, j, 30000.0 * A, 1.0, bar=True)
    f.support("n1", "xy").support("n3", "y")
    for n in ("n2", "n4", "n5"):
        f.load(n, fy=-20.0)
    return f.solve()


@case
def staad_space_truss():
    f = SpaceFrame()
    f.node("n1", 0, 0, 0).node("n2", 96, 0, 0).node("n3", 0, 72, 0).node("n4", 48, 24, 72)
    for name, i in (("m1", "n1"), ("m2", "n2"), ("m3", "n3")):
        f.bar(name, i, "n4", 30.0e6 * 1.0)
        f.support(i, "xyz")
    f.load("n4", [0, 0, 50.0])
    return f.solve()


# ------------------------------------------------------------------------------------------ ANSYS
@case
def ansys_vm2():
    """Beam of 480 in on supports at 120 and 360 in; 10000/12 lb/in (down) on both overhangs; EI = 30e6 x 7892."""
    w = 10000.0 / 12.0
    s = (Beam(480.0, 30.0e6 * 7892.0).support(120.0, v="fixed").support(360.0, v="fixed")
         .dist(0.0, 120.0, -w).dist(360.0, 480.0, -w).solve())
    return {"v_mid": s.v(240.0), "M_mid": s.M(240.0), "v_end": s.v(0.0, +1), "R": s.reactions[120.0][0]}


@case
def ansys_vm21():
    """Tie rod 200 in, 2.5 in square, pin-ended, tension 21,972.6 lb, lateral load 1.79253 lb/in (down)."""
    L, EI, p = 200.0, 30.0e6 * 2.5 ** 4 / 12.0, 1.79253
    out = {}
    for key, kw in (("stiffened", {"N": 21972.6, "second_order": True}), ("unstiffened", {})):
        s = Beam(L, EI, **kw).support(0.0, v="fixed").support(L, v="fixed").dist(0.0, L, -p).solve()
        out[key] = {"v_mid": s.v(L / 2), "theta0": s.theta(0.0, +1), "M_mid": s.M(L / 2)}
    return out


@case
def ansys_vm41():
    """Beam of 10 in (EI = 30e6 x 1000), held at A against translation, rotational spring 10,000 in lb/rad at A,
    10 lb (down) at the free end B."""
    s = Beam(10.0, 30.0e6 * 1000.0).support(0.0, v="fixed", r=10000.0).point(10.0, P=-10.0).solve()
    return {"v_B": s.v(10.0), "theta_B": s.theta(10.0, -1), "theta_A": s.theta(0.0, +1)}


@case
def ansys_vm127():
    return {"Pcr": euler_pinned(30.0e6 * 0.5 ** 4 / 12.0, 200.0)}


@case
def ansys_vm217():
    """Portal, columns 400 in (I = 20,300 in^4), span 800 in (I = 5 x 20,300), fixed bases, 500 lb/in (down) on the
    span, E = 30e6 psi; axial deformation suppressed (EA = 1e7 in^2 x E), as in the page's formulas."""
    f = PlaneFrame()
    for n, (x, y) in {"n1": (0, 0), "n2": (0, 400), "m": (400, 400), "n3": (800, 400), "n4": (800, 0)}.items():
        f.node(n, x, y)
    E, A = 30.0e6, 1.0e7
    f.member("c1", "n1", "n2", E * A, E * 20300.0).member("c2", "n4", "n3", E * A, E * 20300.0)
    f.member("b1", "n2", "m", E * A, E * 101500.0).member("b2", "m", "n3", E * A, E * 101500.0)
    f.uniform("b1", gy=-500.0).uniform("b2", gy=-500.0).support("n1", "xyr").support("n4", "xyr")
    return f.solve()


def main():
    np.set_printoptions(precision=9, suppress=False, linewidth=160)

    def show(v, indent="   "):
        if isinstance(v, dict):
            for k, x in v.items():
                if isinstance(x, dict):
                    print(f"{indent}{k}:")
                    show(x, indent + "   ")
                else:
                    print(f"{indent}{k}: {np.asarray(x) if not isinstance(x, float) else repr(x)}")
    for name, fn in REFERENCES.items():
        print(name)
        show(fn())


if __name__ == "__main__":
    main()
