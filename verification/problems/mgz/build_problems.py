#!/usr/bin/env python3
"""Writes the problem directories of verification/problems/mgz/ (model.json + problem.json).

Source of every problem: W. McGuire, R. H. Gallagher, R. D. Ziemian, "Matrix Structural Analysis",
2nd edition.  Nothing of the book is reproduced here: each entry gives the input data that define the
problem and the printed target numbers with example number and page, in this project's own words.

Rules followed by this file (docs/VERIFICATION.md, section 5):

* a published target is entered AS PRINTED (``P("2.41", ...)``); its tolerance is half a unit in the
  last printed digit, computed by ``benchfea.verify.printed_precision_tolerance`` -- never by hand;
* a target "recomputed from a published statement" (a closed-form expression printed in the book,
  evaluated here for a numeric instance chosen here) carries the policy floor 1e-6;
* ``status`` is set by hand from the outcome of ``check_local.py``: a problem that does not reproduce is
  "quarantined" with a diagnosis; no tolerance, model or target is changed to make it pass.

Run:  python3 verification/problems/mgz/build_problems.py     (rewrites the directories)
"""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

from benchfea.verify.policy import printed_precision_tolerance

HERE = Path(__file__).resolve().parent
WORK = "McGuire, W., Gallagher, R. H., Ziemian, R. D., Matrix Structural Analysis"
EDITION = "2nd edition (copyright 2014 R. D. Ziemian)"
FLOOR = 1.0e-6

PROBLEMS: dict = {}


# ----------------------------------------------------------------------------------------- helpers
def P(printed: str, where: dict, unit: str, *, sign: float = 1.0, factor: float = 1.0, note: str = "",
      signnote: str = "same sign convention as the source (global axes, right-hand rule)") -> dict:
    """A published target: value = sign * factor * (the number as printed)."""
    tol = printed_precision_tolerance(printed)
    if "abs" in tol:
        tol = {"abs": tol["abs"] * abs(factor)}
    pn = f'printed "{printed}"'
    if factor != 1.0:
        pn += f" x {factor:g} (unit conversion / common factor of the source)"
    pn += "; tolerance = half a unit in the last printed digit"
    if note:
        pn += "; " + note
    t = dict(where)
    t.update({"value": sign * factor * float(printed.replace(",", "")), "unit": unit, "tolerance": tol,
              "precisionNote": pn, "sign": signnote})
    return t


def R(value: float, where: dict, unit: str, formula: str, signnote: str = "engine convention") -> dict:
    """A target recomputed here from a closed-form statement printed in the source."""
    t = dict(where)
    t.update({"value": float(value), "unit": unit, "tolerance": {"rel": FLOOR},
              "precisionNote": "evaluated here in double precision from the expression printed in the source ("
                               + formula + "); tolerance = policy floor 1e-6", "sign": signnote})
    return t


def disp(node, dof, comb=None):
    d = {"extract": "nodeDisplacement", "node": node, "dof": dof}
    if comb:
        d["combination"] = comb
    return d


def reac(node, comp, comb=None):
    d = {"extract": "reaction", "node": node, "component": comp}
    if comb:
        d["combination"] = comb
    return d


def station(member, q, rel=None, x=None, side=None, comb=None):
    d = {"extract": "memberStation", "member": member, "quantity": q}
    if rel is not None:
        d["relative"] = rel
    else:
        d["x"] = x
    if side:
        d["side"] = side
    if comb:
        d["combination"] = comb
    return d


def endf(member, end, comp, comb=None):
    d = {"extract": "memberEndForce", "member": member, "end": end, "component": comp}
    if comb:
        d["combination"] = comb
    return d


def buck(mode=1, comb="C"):
    return {"extract": "bucklingFactor", "combination": comb, "mode": mode}


def model(name, units, vertical, nodes, supports, materials, sections, members, load_cases, combinations,
          analysis=None, description=""):
    m = {"name": "mgz/" + name, "formatVersion": "1.0.0"}
    if description:
        m["description"] = description
    m.update({"units": {"force": units[0], "length": units[1]}, "conventions": {"verticalAxis": vertical},
              "nodes": {k: {"x": float(v[0]), "y": float(v[1]), "z": float(v[2]) if len(v) > 2 else 0.0}
                        for k, v in nodes.items()},
              "supports": supports, "materials": materials, "sections": sections, "members": members,
              "loadCases": load_cases, "combinations": combinations,
              "analysis": analysis if analysis is not None else {"shearDeformation": False}})
    return m


def sup(*dofs):
    return {"dofs": {d: "fixed" for d in dofs}}


FIXED = {"type": "fixed"}
PIN2D = sup("ux", "uy", "uz")            # plane truss joint held in x and y (and out of plane)
OUT = ("uz", "rx", "ry")                 # out-of-plane DOFs of a plane (x-y) model


def bar(a):
    """Section of an axial-force member: only A is used by a truss member."""
    return {"type": "properties", "description": "axial-force member: Iy, Iz, J are placeholders, not used",
            "A": a, "Iy": 1.0, "Iz": 1.0, "J": 1.0}


def sec(a, iz, j=None, iy=None, **kw):
    s = {"type": "properties", "A": a, "Iy": iy if iy is not None else iz, "Iz": iz, "J": j if j is not None else iz}
    s.update(kw)
    return s


def truss(i, j, s, mat="M"):
    return {"type": "truss", "i": i, "j": j, "section": s, "material": mat}


def frame(i, j, s, mat="M", **kw):
    m = {"i": i, "j": j, "section": s, "material": mat}
    m.update(kw)
    return m


def one_combo(loads, kind=None):
    c = {"factors": {"L": 1.0}}
    if kind:
        c["analysisType"] = kind
    return {"L": {"loads": loads}}, {"C": c}


def nodal(node, **f):
    d = {"type": "nodal", "node": node}
    d.update(f)
    return d


def settle(node, **u):
    d = {"type": "supportDisplacement", "node": node}
    d.update(u)
    return d


STEEL_SI = {"M": {"E": 200.0, "nu": 0.3}}            # kN/mm^2 (the book's E = 200,000 MPa)
STEEL_US = {"M": {"E": 29000.0, "nu": 0.3}}          # ksi

PLANE_NOTE = ("Plane (x-y) problem of the source modelled in 3D with the global Y axis vertical, so that the "
              "source's x, y and theta_z are the model's X, Y and rz. In-plane bending is about local z for every "
              "member (default orientation), hence the source's I is the section's Iz. ")
TRUSS_NOTE = ("Pin-jointed plane truss in the x-y plane (global Y vertical): members of type 'truss'; the "
              "out-of-plane translation of every joint is held by a support (it has no stiffness in a plane truss). ")
UNUSED = ("Section properties the source does not give because they do not enter its plane model (Iy, J; A where "
          "stated) are placeholders that do not influence the targets. ")
BAR_SIGN = "bar force: tension positive in the source and in the engine (N)"
A_RIGID = 1.0e8          # mm^2: member area used where the source neglects axial deformation
RIGID_NOTE = ("The engine always includes axial deformation; the source's assumption is represented by a member area "
              "of 1e8 mm^2 for every member (four orders of magnitude above the areas the source uses for the same "
              "members elsewhere, 4e3 and 6e3 mm^2), one value for all problems of this kind. The remaining axial "
              "flexibility changes the targets by less than 1e-5 of their value (measured by repeating the analysis "
              "with 1e9 mm^2; see README.md), far below their printed precision.")


def add(pid, title, location, pages, analysis, mdl, targets, *, basis="published", status="active",
        diagnosis=None, notes="", limits=""):
    src = {"work": WORK, "edition": EDITION, "location": location, "pages": pages,
           "notes": notes or "input data and targets read on the cited pages; nothing of the source is reproduced"}
    if pid in QUARANTINE:
        status, diagnosis = "quarantined", QUARANTINE[pid]
    p = {"id": pid, "title": title, "tier": "B", "source": src, "basis": basis, "analysis": analysis,
         "targets": targets, "status": status}
    if status == "quarantined":
        p["quarantine"] = {"diagnosis": diagnosis}
    elif diagnosis:
        raise ValueError(pid + ": a diagnosis belongs to a quarantined problem")
    if limits:
        p["limitsNote"] = limits
    PROBLEMS[pid] = (p, mdl)


from quarantine import QUARANTINE            # noqa: E402
from problems_ch2_5 import define as _d1     # noqa: E402
from problems_ch7_13 import define as _d2    # noqa: E402


def main():
    import sys
    ns = sys.modules[__name__]
    _d1(ns)
    _d2(ns)
    unknown = set(QUARANTINE) - set(PROBLEMS)
    if unknown:
        raise SystemExit("quarantine.py names unknown problems: " + ", ".join(sorted(unknown)))
    for d in HERE.iterdir():                               # remove directories of problems no longer defined
        if d.is_dir() and (d / "problem.json").is_file() and d.name not in PROBLEMS:
            shutil.rmtree(d)
    for pid, (p, mdl) in PROBLEMS.items():
        d = HERE / pid
        d.mkdir(exist_ok=True)
        (d / "problem.json").write_text(json.dumps(p, indent=1) + "\n", encoding="utf-8")
        (d / "model.json").write_text(json.dumps(mdl, indent=1) + "\n", encoding="utf-8")
    print(f"{len(PROBLEMS)} problems written to {HERE}")


if __name__ == "__main__":
    main()
