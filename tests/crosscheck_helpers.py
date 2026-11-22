"""Shared helpers of the cross-check tests (model documents in N, mm; solver availability)."""
import importlib.util

import numpy as np
import pytest

E, G = 200000.0, 77000.0
SEC = {"A": 5000.0, "Iy": 8.0e6, "Iz": 4.0e7, "J": 3.0e5, "Asy": 2500.0, "Asz": 3000.0}
EA, EIy, EIz, GJ = E * SEC["A"], E * SEC["Iy"], E * SEC["Iz"], G * SEC["J"]
GAy, GAz = G * SEC["Asy"], G * SEC["Asz"]

HAVE = {"opensees": importlib.util.find_spec("openseespy") is not None,
        "pynite": importlib.util.find_spec("Pynite") is not None}
SOLVER_PARAMS = [pytest.param(s, marks=pytest.mark.skipif(not HAVE[s], reason=f"{s} is not installed"))
                 for s in ("opensees", "pynite")]
needs_opensees = pytest.mark.skipif(not HAVE["opensees"], reason="openseespy is not installed")
needs_pynite = pytest.mark.skipif(not HAVE["pynite"], reason="PyNiteFEA is not installed")
needs_both = pytest.mark.skipif(not (HAVE["opensees"] and HAVE["pynite"]), reason="needs openseespy and PyNiteFEA")


def document(nodes, members, supports, cases, combos=None, analysis=None, constraints=None, sections=None,
             materials=None, vertical="Z", units=None):
    """Input document.  nodes {id: (x, y, z)}, members {id: (i, j) | dict}, supports {id: "fixed" | dict},
    cases {id: [load, ...]}; one combination per case unless ``combos`` is given; shear deformation is
    OFF unless ``analysis`` says otherwise."""
    doc = {
        "formatVersion": "1.0.0", "units": units or {"force": "N", "length": "mm"},
        "conventions": {"verticalAxis": vertical},
        "nodes": {k: {"x": float(v[0]), "y": float(v[1]), "z": float(v[2])} for k, v in nodes.items()},
        "supports": {k: ({"type": v} if isinstance(v, str) else v) for k, v in supports.items()},
        "materials": materials or {"M": {"E": E, "G": G, "unitWeight": 7.7e-5}},
        "sections": sections or {"S": {"type": "properties", **SEC}},
        "members": {}, "loadCases": {k: {"loads": list(v)} for k, v in cases.items()},
        "combinations": combos or {k: {"factors": {k: 1.0}} for k in cases},
        "analysis": {"shearDeformation": False, **(analysis or {})},
    }
    for k, v in members.items():
        m = {"i": v[0], "j": v[1]} if isinstance(v, (tuple, list)) else dict(v)
        m.setdefault("section", next(iter(doc["sections"])))
        m.setdefault("material", next(iter(doc["materials"])))
        doc["members"][k] = m
    if constraints:
        doc["constraints"] = constraints
    return doc


def solve(solver, model, combination, **options):
    """Run a foreign solver: PyNite in-process (pure Python), OpenSees in the subprocess."""
    from benchfea.crosscheck import run

    return run(solver, model, combination, in_process=(solver == "pynite"), **options)


def local_node(model, result, member, node):
    """(u_local(3), theta_local(3)) of a node in the axes of a member."""
    R = np.asarray(model.members[member].R, float)
    u = result.displacements[model.nodes[node].index]
    return R @ u[:3], R @ u[3:]


def rel(a, b, scale=None):
    return abs(a - b) / (scale if scale else max(abs(b), 1e-300))
