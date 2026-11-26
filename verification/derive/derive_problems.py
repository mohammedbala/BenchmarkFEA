#!/usr/bin/env python3
"""Regenerate the tier-A problems of ``verification/problems/derived/`` ("derived here").

    python3 verification/derive/derive_problems.py           # (re)write the problem directories
    python3 verification/derive/derive_problems.py --check   # exit 1 if the stored targets differ

Every target is computed HERE, not taken from any book:

* beams and beam-columns: direct integration of the member differential equations by the independent
  ODE oracle ``benchfea.oracle`` (``scipy.integrate.solve_bvp``; equations derived in its docstring).
  The oracle problem of each case is written down by hand below from the same numbers as the model
  file -- it does not read the model file and shares no code with ``benchfea.solve``;
* Euler columns: the characteristic equation of each end condition is formed from the general solution
  v = A sin kx + B cos kx + C x + D of EI v'''' + P v'' = 0 (k^2 = P / EI) and the four end conditions,
  and its smallest root is found numerically (4 x 4 determinant, sign scan + Brent) -- no tabulated
  effective-length factor is used.

Geometry of every beam model: member along global +X, Z up, default orientation, so that local y = +Z
and local z = -Y (docs/CONVENTIONS.md 2).  Hence, for bending in the vertical plane:
deflection v = uz;  member rotation about local z = -(global ry);  a moment about local z = -(global My).

Tolerance of every target: relative 1e-6 (``benchfea.verify.POLICY["targets"]["derivedHere"]``; the
oracle is accurate to about 1e-8), absolute for a target whose value is zero.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from benchfea.oracle import Bar, Beam  # noqa: E402

OUT = ROOT / "verification" / "problems" / "derived"
REL = 1.0e-6

E = 2.0e8            # kN/m^2
NU = 0.3
G = E / (2.0 * (1.0 + NU))
A, IZ, IY, J = 0.01, 2.0e-4, 1.0e-4, 1.5e-4      # m^2, m^4
EI = E * IZ


def base(nodes, supports, members, loads, *, factor=1.0, analysis=None, section=None, material=None, combo=None):
    sec = {"type": "properties", "A": A, "Iy": IY, "Iz": IZ, "J": J}
    sec.update(section or {})
    mat = {"E": E, "nu": NU}
    mat.update(material or {})
    an = {"shearDeformation": False, "resultStations": 11}
    an.update(analysis or {})
    c = {"factors": {"L": factor}}
    c.update(combo or {})
    return {"formatVersion": "1.0.0", "units": {"force": "kN", "length": "m"}, "conventions": {"verticalAxis": "Z"},
            "nodes": nodes, "supports": supports, "materials": {"M": mat}, "sections": {"S": sec},
            "members": {k: {"section": "S", "material": "M", **v} for k, v in members.items()},
            "loadCases": {"L": {"loads": loads}}, "combinations": {"C": c}, "analysis": an}


def xn(*xs):
    return {name: {"x": float(x), "y": 0.0, "z": 0.0} for name, x in zip("abcdefgh", xs)}


PIN = {"dofs": {"ux": "fixed", "uy": "fixed", "uz": "fixed", "rx": "fixed"}}
ROLLER = {"dofs": {"uy": "fixed", "uz": "fixed"}}
FIXED = {"type": "fixed"}


def T(value, unit, **how):
    t = dict(how)
    t["value"] = float(value)
    t["unit"] = unit
    t["tolerance"] = {"rel": REL} if value != 0.0 else {"abs": 1.0e-9}
    return t


def disp(node, dof, value, unit="m"):
    return T(value, unit, path=f"/combinations/C/displacements/{node}/{dof}")


def reac(node, comp, value, unit="kN"):
    return T(value, unit, path=f"/combinations/C/reactions/{node}/{comp}")


def station(member, quantity, rel, value, unit, side="left"):
    return T(value, unit, extract="memberStation", combination="C", member=member, quantity=quantity,
             relative=rel, side=side)


PROBLEMS = {}


def problem(pid, title, analysis_text):
    def deco(fn):
        PROBLEMS[pid] = (title, analysis_text, fn)
        return fn
    return deco


# ============================================================================================ beams
@problem("beam-point-load", "Simply supported beam, point load at 0.3 L (factored combination)",
         "First order, no shear deformation. Span 6 m, EI = 40 000 kN m^2, P = 25 kN down at x = 1.8 m, "
         "combination factor 1.5.")
def beam_point_load():
    L, P, f = 6.0, 25.0, 1.5
    model = base(xn(0, L), {"a": PIN, "b": ROLLER}, {"m": {"i": "a", "j": "b"}},
                 [{"type": "memberPoint", "member": "m", "direction": "Z", "value": -P, "position": 0.3}], factor=f)
    s = Beam(L, EI).support(0, v="fixed").support(L, v="fixed").point(0.3 * L, P=-P * f).solve()
    return model, [
        reac("a", "Fz", s.reactions[0.0][0]), reac("b", "Fz", s.reactions[L][0]),
        disp("a", "ry", -s.theta(0.0, +1), "rad"), disp("b", "ry", -s.theta(L, -1), "rad"),
        station("m", "uy", 0.5, s.v(0.5 * L), "m"), station("m", "uy", 0.3, s.v(0.3 * L), "m"),
        station("m", "Mz", 0.3, s.M(0.3 * L, -1), "kN m"),
        station("m", "Vy", 0.3, s.Q(0.3 * L, -1), "kN", side="left"),
        station("m", "Vy", 0.3, s.Q(0.3 * L, +1), "kN", side="right"),
    ]


@problem("beam-trapezoid-partial", "Propped cantilever, partial linearly varying load",
         "First order, no shear deformation. Span 5 m, fixed at a, roller at b; line load from 12 kN/m at "
         "x = 1 m to 4 kN/m at x = 4 m (down).")
def beam_trapezoid():
    L, w1, w2 = 5.0, 12.0, 4.0
    model = base(xn(0, L), {"a": FIXED, "b": ROLLER}, {"m": {"i": "a", "j": "b"}},
                 [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -w1, "w2": -w2,
                   "start": 0.2, "end": 0.8}])
    s = (Beam(L, EI).support(0, v="fixed", r="fixed").support(L, v="fixed").dist(1.0, 4.0, -w1, -w2).solve())
    return model, [
        reac("a", "Fz", s.reactions[0.0][0]), reac("a", "My", -s.reactions[0.0][1], "kN m"),
        reac("b", "Fz", s.reactions[L][0]), disp("b", "ry", -s.theta(L, -1), "rad"),
        station("m", "uy", 0.5, s.v(2.5), "m"), station("m", "Mz", 0.5, s.M(2.5), "kN m"),
        station("m", "Mz", 0.0, s.M(0.0, +1), "kN m"), station("m", "Vy", 0.8, s.Q(4.0), "kN"),
    ]


@problem("beam-point-moment", "Simply supported beam, concentrated moment at 0.4 L",
         "First order, no shear deformation. Span 6 m; moment of 30 kN m about the member's local z axis at "
         "x = 2.4 m.")
def beam_point_moment():
    L, M0 = 6.0, 30.0
    model = base(xn(0, L), {"a": PIN, "b": ROLLER}, {"m": {"i": "a", "j": "b"}},
                 [{"type": "memberPoint", "member": "m", "kind": "moment", "direction": "z", "value": M0,
                   "position": 0.4}])
    s = Beam(L, EI).support(0, v="fixed").support(L, v="fixed").point(0.4 * L, M=M0).solve()
    return model, [
        reac("a", "Fz", s.reactions[0.0][0]), reac("b", "Fz", s.reactions[L][0]),
        station("m", "Mz", 0.4, s.M(2.4, -1), "kN m", side="left"),
        station("m", "Mz", 0.4, s.M(2.4, +1), "kN m", side="right"),
        station("m", "uy", 0.4, s.v(2.4), "m"), station("m", "uy", 0.7, s.v(4.2), "m"),
        disp("a", "ry", -s.theta(0.0, +1), "rad"),
    ]


@problem("cantilever-nodal-loads", "Cantilever with tip force, tip moment, axial force and torque",
         "First order, no shear deformation. Length 3 m; at the tip: Fz = -10 kN, My = +5 kN m (global), "
         "Fx = 20 kN, Mx = 2 kN m. Bending, axial and torsional responses are uncoupled.")
def cantilever_nodal():
    L = 3.0
    model = base(xn(0, L), {"a": FIXED}, {"m": {"i": "a", "j": "b"}},
                 [{"type": "nodal", "node": "b", "Fz": -10.0, "My": 5.0, "Fx": 20.0, "Mx": 2.0}])
    # global My = +5 is a moment of -5 about local z (local z = -Y)
    s = Beam(L, EI).support(0, v="fixed", r="fixed").point(L, P=-10.0, M=-5.0).solve()
    ax = Bar(L, E * A).support(0).point(L, 20.0).solve()
    tw = Bar(L, G * J).support(0).point(L, 2.0).solve()
    return model, [
        disp("b", "uz", s.v(L)), disp("b", "ry", -s.theta(L, -1), "rad"),
        disp("b", "ux", ax.u(L, -1)), disp("b", "rx", tw.u(L, -1), "rad"),
        reac("a", "Fz", s.reactions[0.0][0]), reac("a", "My", -s.reactions[0.0][1], "kN m"),
        reac("a", "Fx", -20.0), reac("a", "Mx", -2.0, "kN m"),
        station("m", "Mz", 0.5, s.M(1.5), "kN m"), station("m", "N", 0.5, ax.N(1.5, -1), "kN"),
        station("m", "T", 0.5, tw.N(1.5, -1), "kN m"),
    ]


@problem("beam-self-weight", "Simply supported beam under its own weight (factor 1.35)",
         "First order, no shear deformation. Span 8 m, A = 0.01 m^2, unit weight 78.5 kN/m^3: line load "
         "0.785 kN/m along -Z, combination factor 1.35.")
def beam_self_weight():
    L, gamma, f = 8.0, 78.5, 1.35
    model = base(xn(0, L), {"a": PIN, "b": ROLLER}, {"m": {"i": "a", "j": "b"}}, [{"type": "selfWeight"}],
                 factor=f, material={"unitWeight": gamma})
    w = gamma * A * f
    s = Beam(L, EI).support(0, v="fixed").support(L, v="fixed").dist(0.0, L, -w).solve()
    return model, [
        reac("a", "Fz", s.reactions[0.0][0]), station("m", "uy", 0.5, s.v(4.0), "m"),
        station("m", "Mz", 0.5, s.M(4.0), "kN m"), station("m", "Vy", 0.2, s.Q(1.6), "kN"),
        disp("a", "ry", -s.theta(0.0, +1), "rad"),
    ]


@problem("continuous-beam-settlement", "Two-span continuous beam: support settlement and uniform load",
         "First order, no shear deformation. Spans 5 m + 5 m on three supports; the middle support settles "
         "10 mm; uniform load 6 kN/m on the first span.")
def continuous_settlement():
    q, d = 6.0, -0.010
    model = base(xn(0, 5, 10), {"a": PIN, "b": ROLLER, "c": ROLLER},
                 {"m1": {"i": "a", "j": "b"}, "m2": {"i": "b", "j": "c"}},
                 [{"type": "supportDisplacement", "node": "b", "uz": d},
                  {"type": "memberDistributed", "member": "m1", "direction": "Z", "w1": -q}])
    s = (Beam(10.0, EI).support(0, v="fixed").support(5, v="fixed", settlement=d).support(10, v="fixed")
         .dist(0.0, 5.0, -q).solve())
    return model, [
        reac("a", "Fz", s.reactions[0.0][0]), reac("b", "Fz", s.reactions[5.0][0]),
        reac("c", "Fz", s.reactions[10.0][0]), disp("b", "ry", -s.theta(5.0, -1), "rad"),
        station("m1", "Mz", 1.0, s.M(5.0, -1), "kN m"), station("m1", "uy", 0.5, s.v(2.5), "m"),
        station("m2", "uy", 0.5, s.v(7.5), "m"),
    ]


@problem("timoshenko-fixed-fixed", "Deep fixed-fixed beam with shear deformation: uniform and point load",
         "First order WITH shear deformation (analysis.shearDeformation = true). Span 2 m, EI = 40 000 kN m^2, "
         "shear area 0.008 m^2 (G As = 615 385 kN): 12 EI / (G As L^2) = 0.195. Uniform load 40 kN/m and a "
         "point load of 60 kN at 0.3 L.")
def timoshenko():
    L, q, P, As = 2.0, 40.0, 60.0, 0.008
    model = base(xn(0, L), {"a": FIXED, "b": FIXED}, {"m": {"i": "a", "j": "b"}},
                 [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q},
                  {"type": "memberPoint", "member": "m", "direction": "Z", "value": -P, "position": 0.3}],
                 analysis={"shearDeformation": True}, section={"Asy": As})
    s = (Beam(L, EI, GAs=G * As).support(0, v="fixed", r="fixed").support(L, v="fixed", r="fixed")
         .dist(0.0, L, -q).point(0.3 * L, P=-P).solve())
    return model, [
        reac("a", "Fz", s.reactions[0.0][0]), reac("a", "My", -s.reactions[0.0][1], "kN m"),
        reac("b", "Fz", s.reactions[L][0]), reac("b", "My", -s.reactions[L][1], "kN m"),
        station("m", "uy", 0.5, s.v(1.0), "m"), station("m", "uy", 0.3, s.v(0.6), "m"),
        station("m", "Mz", 0.5, s.M(1.0), "kN m"),
    ]


@problem("beam-on-spring-supports", "Beam on translational and rotational spring supports",
         "First order, no shear deformation. Nodes at 0, 4 and 9 m; a: vertical spring 5000 kN/m and "
         "rotational spring 2000 kN m/rad, b: vertical spring 8000 kN/m, c: fixed vertically. Uniform load "
         "10 kN/m on both members and 30 kN at 0.4 of the first.")
def spring_supports():
    q, P = 10.0, 30.0
    ka, kra, kb = 5000.0, 2000.0, 8000.0
    sup = {"a": {"dofs": {"ux": "fixed", "uy": "fixed", "rx": "fixed", "uz": {"spring": ka}, "ry": {"spring": kra}}},
           "b": {"dofs": {"uy": "fixed", "uz": {"spring": kb}}}, "c": ROLLER}
    model = base(xn(0, 4, 9), sup, {"m1": {"i": "a", "j": "b"}, "m2": {"i": "b", "j": "c"}},
                 [{"type": "memberDistributed", "member": "m1", "direction": "Z", "w1": -q},
                  {"type": "memberDistributed", "member": "m2", "direction": "Z", "w1": -q},
                  {"type": "memberPoint", "member": "m1", "direction": "Z", "value": -P, "position": 0.4}])
    s = (Beam(9.0, EI).support(0, v=ka, r=kra).support(4, v=kb).support(9, v="fixed")
         .dist(0.0, 9.0, -q).point(1.6, P=-P).solve())
    return model, [
        disp("a", "uz", s.v(0.0, +1)), disp("b", "uz", s.v(4.0)), disp("a", "ry", -s.theta(0.0, +1), "rad"),
        reac("a", "Fz", s.reactions[0.0][0]), reac("a", "My", -s.reactions[0.0][1], "kN m"),
        reac("b", "Fz", s.reactions[4.0][0]), reac("c", "Fz", s.reactions[9.0][0]),
        station("m2", "Mz", 0.0, s.M(4.0, +1), "kN m"),
    ]


@problem("propped-cantilever-partial-release", "Propped cantilever with a rotational spring connection at the built-in end",
         "First order, no shear deformation. Span 6 m; the member end at the fixed support a is connected "
         "through a rotational spring of 30 000 kN m/rad (partial release of rz); roller at b; uniform load "
         "15 kN/m. For the beam this is a rotational spring support.")
def partial_release():
    L, q, k = 6.0, 15.0, 3.0e4
    model = base(xn(0, L), {"a": FIXED, "b": ROLLER},
                 {"m": {"i": "a", "j": "b", "releases": {"i": {"rz": {"spring": k}}}}},
                 [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q}])
    s = Beam(L, EI).support(0, v="fixed", r=k).support(L, v="fixed").dist(0.0, L, -q).solve()
    return model, [
        reac("a", "Fz", s.reactions[0.0][0]), reac("a", "My", -s.reactions[0.0][1], "kN m"),
        reac("b", "Fz", s.reactions[L][0]), station("m", "Mz", 0.0, s.M(0.0, +1), "kN m"),
        station("m", "uy", 0.5, s.v(3.0), "m"), station("m", "Mz", 0.5, s.M(3.0), "kN m"),
        disp("b", "ry", -s.theta(L, -1), "rad"),
    ]


# ============================================================================================ beam-columns
def _beam_column(ratio):
    """Simply supported beam-column, equal inertias (so that the other plane has the same critical load),
    uniform load 8 kN/m, axial force = ratio x pi^2 EI / L^2 (compression for ratio > 0)."""
    L, q = 6.0, 8.0
    Pe = math.pi ** 2 * EI / L ** 2
    P = ratio * Pe
    model = base(xn(0, L), {"a": PIN, "b": ROLLER}, {"m": {"i": "a", "j": "b"}},
                 [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q},
                  {"type": "nodal", "node": "b", "Fx": -P}],
                 section={"Iy": IZ}, analysis={"convergenceTolerance": 1.0e-10},
                 combo={"analysisType": "secondOrder"})
    s = (Beam(L, EI, N=-P, second_order=True).support(0, v="fixed").support(L, v="fixed")
         .dist(0.0, L, -q).solve())
    first = Beam(L, EI).support(0, v="fixed").support(L, v="fixed").dist(0.0, L, -q).solve()
    assert (abs(s.v(3.0)) > abs(first.v(3.0))) == (ratio > 0)        # compression amplifies, tension stiffens
    return model, [
        station("m", "uy", 0.5, s.v(3.0), "m"), station("m", "Mz", 0.5, s.M(3.0), "kN m"),
        station("m", "uy", 0.2, s.v(1.2), "m"), station("m", "Mz", 0.2, s.M(1.2), "kN m"),
        disp("a", "ry", -s.theta(0.0, +1), "rad"), reac("a", "Fz", s.reactions[0.0][0]),
        station("m", "N", 0.5, -P, "kN"),
    ]


for _r, _name in ((0.2, "02"), (0.5, "05"), (0.8, "08"), (-0.5, "tension-05")):
    problem(f"beam-column-{_name}",
            f"Simply supported beam-column, uniform load, axial {'compression' if _r > 0 else 'tension'} "
            f"{abs(_r):g} P_e",
            "Second order (member equations with the axial force acting on the deflected axis), no shear "
            "deformation. Span 6 m, EI = 40 000 kN m^2 in both planes, uniform load 8 kN/m, axial force "
            f"{abs(_r):g} x pi^2 EI / L^2 in {'compression' if _r > 0 else 'tension'}.")(lambda r=_r: _beam_column(r))


# ============================================================================================ Euler columns
def _rows(kind, xi, kappa):
    """Rows of the end conditions at xi (0 or 1) for v = A sin(kappa xi) + B cos(kappa xi) + C xi + D,
    lengths in units of L.  v = 0;  slope v' = 0;  moment v'' = 0;  free end: transverse force along the
    undeformed axis v''' + kappa^2 v' = 0 (from M' = -Q + N v' with N = -P and Q = 0)."""
    s, c = math.sin(kappa * xi), math.cos(kappa * xi)
    v = [s, c, xi, 1.0]
    d1 = [kappa * c, -kappa * s, 1.0, 0.0]
    d2 = [-kappa ** 2 * s, -kappa ** 2 * c, 0.0, 0.0]
    d3 = [-kappa ** 3 * c, kappa ** 3 * s, 0.0, 0.0]
    force = [a + kappa ** 2 * b for a, b in zip(d3, d1)]
    return {"pinned": [v, d2], "fixed": [v, d1], "free": [d2, force]}[kind]


def critical_kappa(bottom, top):
    """Smallest positive root kappa = L sqrt(P / EI) of the characteristic determinant."""
    def det(kappa):
        M = np.array(_rows(bottom, 0.0, kappa) + _rows(top, 1.0, kappa))
        M = M / np.abs(M).max(axis=1, keepdims=True)
        return float(np.linalg.det(M))
    grid = np.linspace(0.2, 12.0, 2361)
    vals = [det(k) for k in grid]
    for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if fa * fb < 0.0:
            return brentq(det, a, b, xtol=1e-14, rtol=1e-14)
    raise RuntimeError("no root found")


def _column(bottom, top):
    L, Pref = 4.0, 100.0
    kappa = critical_kappa(bottom, top)
    Pcr = kappa ** 2 * EI / L ** 2
    nodes = {"a": {"x": 0.0, "y": 0.0, "z": 0.0}, "b": {"x": 0.0, "y": 0.0, "z": L}}
    sup = {"a": {"fixed": {"type": "fixed"},
                 "pinned": {"dofs": {"ux": "fixed", "uy": "fixed", "uz": "fixed", "rz": "fixed"}}}[bottom]}
    t = {"pinned": {"dofs": {"ux": "fixed", "uy": "fixed"}},
         "fixed": {"dofs": {"ux": "fixed", "uy": "fixed", "rx": "fixed", "ry": "fixed"}}, "free": None}[top]
    if t:
        sup["b"] = t
    model = base(nodes, sup, {"m": {"i": "a", "j": "b"}}, [{"type": "nodal", "node": "b", "Fz": -Pref}],
                 section={"Iy": IZ}, analysis={"buckling": {"combinations": ["C"], "modes": 1}})
    t = T(Pcr / Pref, "-", extract="bucklingFactor", combination="C", mode=1)
    t["note"] = f"kappa = L sqrt(Pcr / EI) = {kappa:.12f}; Pcr = {Pcr:.9g} kN; reference load {Pref:g} kN"
    return model, [t]


for _b, _t in (("pinned", "pinned"), ("fixed", "free"), ("fixed", "pinned"), ("fixed", "fixed")):
    problem(f"euler-column-{_b}-{_t}", f"Elastic critical load of a column, {_b} base, {_t} top",
            "Linear buckling of a prismatic column, no shear deformation. Length 4 m, EI = 40 000 kN m^2 in "
            "both planes, reference load 100 kN. Target = smallest root of the characteristic determinant of "
            "EI v'''' + P v'' = 0 with the four end conditions, solved numerically.")(
        lambda b=_b, t=_t: _column(b, t))


# ============================================================================================ quarantine
# Problems whose TARGETS are reproduced (tier A) but whose model fails another tier.  They stay in the
# library, quarantined with the diagnosis, and are never counted as passed.
QUARANTINE = {}


# ============================================================================================ output
def build(pid):
    title, text, fn = PROBLEMS[pid]
    model, targets = fn()
    model = {"name": f"derived/{pid}", **model}
    prob = {
        "id": pid, "title": title, "tier": "A",
        "source": {"work": "derived here (no publication): benchfea.oracle direct integration of the member "
                           "equations / numerically solved characteristic equation",
                   "edition": "", "location": f"verification/derive/derive_problems.py, problem '{pid}'",
                   "pages": "", "notes": "regenerate with: python3 verification/derive/derive_problems.py"},
        "basis": "derived-here", "analysis": text, "targets": targets, "status": "active",
    }
    if pid in QUARANTINE:
        prob["status"] = "quarantined"
        prob["quarantine"] = {"diagnosis": QUARANTINE[pid]}
    return model, prob


def close(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(close(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(close(x, y) for x, y in zip(a, b))
    if isinstance(a, float) and isinstance(b, (int, float)) and not isinstance(b, bool):
        return abs(a - b) <= 1e-9 * max(abs(a), abs(b), 1e-300) or abs(a - b) <= 1e-15
    return a == b


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    check = "--check" in argv
    bad = []
    for pid in PROBLEMS:
        model, prob = build(pid)
        d = OUT / pid
        if check:
            try:
                old_m = json.loads((d / "model.json").read_text(encoding="utf-8"))
                old_p = json.loads((d / "problem.json").read_text(encoding="utf-8"))
            except OSError:
                bad.append(f"{pid}: missing")
                continue
            if not (close(model, old_m) and close(prob, old_p)):
                bad.append(f"{pid}: stored files differ from the regenerated ones")
            continue
        d.mkdir(parents=True, exist_ok=True)
        (d / "model.json").write_text(json.dumps(model, indent=1) + "\n", encoding="utf-8")
        (d / "problem.json").write_text(json.dumps(prob, indent=1) + "\n", encoding="utf-8")
    if check:
        print("\n".join(bad) if bad else f"{len(PROBLEMS)} derived problems reproduce")
        return 1 if bad else 0
    print(f"wrote {len(PROBLEMS)} problems to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
