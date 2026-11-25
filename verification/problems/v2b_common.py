"""Helpers shared by the problem builders of aisc-commentary/, ziemian/, aisc-de/, staad/ and ansys/.

Tolerance of a target (docs/VERIFICATION.md section 5; nothing is set by hand per case):

* ``P(printed, ...)``  a published number, entered AS PRINTED: half a unit in its last printed digit
  (``benchfea.verify.policy.printed_precision_tolerance``; floor 1e-6);
* ``S(printed, rel, statement, ...)``  a published number for which the SOURCE ITSELF states an
  acceptance tolerance (AISC 360-22 Commentary C2.1: moments 3 %, deflections 5 %): that tolerance;
* ``R(value, ...)``  a number computed here (ODE oracle / own derivation) from a published problem
  statement: the policy floor 1e-6.

``status`` of a problem is "active" unless the builder passes a quarantine diagnosis; a diagnosis is
written after ``check_local_v2b.py`` showed that the problem does not reproduce (or cannot be used).
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from benchfea.verify.policy import printed_precision_tolerance

FLOOR = 1.0e-6


def P(printed: str, where: dict, unit: str, *, sign: float = 1.0, factor: float = 1.0, note: str = "",
      signnote: str = "same sign as the source") -> dict:
    """A published target: value = sign * factor * (the number as printed)."""
    tol = printed_precision_tolerance(printed)
    if "abs" in tol:
        tol = {"abs": tol["abs"] * abs(factor)}
    pn = f'printed "{printed}"'
    if factor != 1.0:
        pn += f" x {factor:.12g} (unit conversion / common factor of the source)"
    pn += "; tolerance = half a unit in the last printed digit"
    if note:
        pn += "; " + note
    t = dict(where)
    t.update({"value": sign * factor * float(printed.replace(",", "")), "unit": unit, "tolerance": tol,
              "precisionNote": pn, "sign": signnote})
    return t


def S(printed: str, rel: float, statement: str, where: dict, unit: str, *, sign: float = 1.0, factor: float = 1.0,
      note: str = "", signnote: str = "same sign as the source") -> dict:
    """A published target with an acceptance tolerance stated by the source itself."""
    pn = f'stated tolerance {rel:g} on printed "{printed}"'
    if factor != 1.0:
        pn += f" x {factor:.12g} (unit conversion / common factor of the source)"
    pn += "; " + statement
    if note:
        pn += "; " + note
    t = dict(where)
    t.update({"value": sign * factor * float(printed.replace(",", "")), "unit": unit, "tolerance": {"rel": rel},
              "precisionNote": pn, "sign": signnote})
    return t


def R(value: float, where: dict, unit: str, how: str, signnote: str = "engine convention (docs/CONVENTIONS.md)") -> dict:
    """A target computed here from the published problem statement."""
    t = dict(where)
    t.update({"value": float(value), "unit": unit,
              "tolerance": {"rel": FLOOR} if value != 0.0 else {"abs": 1.0e-9},
              "precisionNote": "evaluated here in double precision (" + how + "); tolerance = policy floor 1e-6",
              "sign": signnote})
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


def extreme(member, q, kind="absmax", comb=None):
    d = {"extract": "memberExtreme", "member": member, "quantity": q, "kind": kind}
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


def sup(*dofs):
    return {"dofs": {d: "fixed" for d in dofs}}


OUT = ("uz", "rx", "ry")                 # out-of-plane DOFs of a plane (x-y) model
FIXED = {"type": "fixed"}


def model(name, units, vertical, nodes, supports, materials, sections, members, load_cases, combinations,
          analysis=None, description=""):
    m = {"name": name, "formatVersion": "1.0.0"}
    if description:
        m["description"] = description
    m.update({"units": {"force": units[0], "length": units[1]}, "conventions": {"verticalAxis": vertical},
              "nodes": {k: {"x": float(v[0]), "y": float(v[1]), "z": float(v[2]) if len(v) > 2 else 0.0}
                        for k, v in nodes.items()},
              "supports": supports, "materials": materials, "sections": sections, "members": members,
              "loadCases": load_cases, "combinations": combinations,
              "analysis": analysis if analysis is not None else {"shearDeformation": False}})
    return m


class Library:
    """Collects the problems of one source directory and writes them."""

    def __init__(self, here: Path, slug: str):
        self.here, self.slug, self.problems = Path(here), slug, {}

    def add(self, pid, title, source, analysis, mdl, targets, *, basis="published", tier="B", diagnosis=None,
            limits=""):
        if pid in self.problems:
            raise ValueError("duplicate problem id " + pid)
        p = {"id": pid, "title": title, "tier": tier, "source": source, "basis": basis, "analysis": analysis,
             "targets": targets, "status": "quarantined" if diagnosis else "active"}
        if diagnosis:
            p["quarantine"] = {"diagnosis": diagnosis}
        if limits:
            p["limitsNote"] = limits
        mdl = dict(mdl)
        mdl["name"] = f"{self.slug}/{pid}"
        self.problems[pid] = (p, mdl)

    def write(self):
        for d in self.here.iterdir():                      # remove directories of problems no longer defined
            if d.is_dir() and (d / "problem.json").is_file() and d.name not in self.problems:
                shutil.rmtree(d)
        for pid, (p, mdl) in self.problems.items():
            d = self.here / pid
            d.mkdir(exist_ok=True)
            (d / "problem.json").write_text(json.dumps(p, indent=1) + "\n", encoding="utf-8")
            (d / "model.json").write_text(json.dumps(mdl, indent=1) + "\n", encoding="utf-8")
        print(f"{len(self.problems)} problems written to {self.here}")
