"""Helpers shared by staad/build_problems.py and ansys/build_problems.py.

A vendor verification problem becomes

* ``<id>``: basis "recomputed-from-published-statement".  The problem STATEMENT is taken from the vendor
  page; the targets are (a) reference values computed here (verification/derive/vendor_reference.py;
  tolerance: policy floor 1e-6) and (b) the vendor's printed "theory" / "target" values that reproduce at
  their printed precision (half a unit in the last printed digit);
* ``<id>-vendor-printed`` (only where needed): the vendor's printed values that do NOT reproduce at their
  printed precision, quarantined with a diagnosis, so that the recomputed targets stay active.

The textbooks the vendor pages cite were not opened for this project; nothing is taken from them.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "derive"))
from v2b_common import (FIXED, OUT, Library, P, R, buck, disp, endf, extreme, model, reac, station,  # noqa: E402,F401
                        sup)
from vendor_reference import reference  # noqa: E402,F401

DERIVED = "verification/derive/vendor_reference.py"


def plane(nodes, supports=None):
    """Supports of a plane (x-y) model: the out-of-plane DOFs of every node are held."""
    out = {}
    for n in nodes:
        dofs = {d: "fixed" for d in OUT}
        for d, v in (supports or {}).get(n, {}).items():
            dofs[d] = v
        out[n] = {"dofs": dofs}
    return out


def fix(*dofs):
    return {d: "fixed" for d in dofs}


def truss_supports(nodes, supports):
    """Plane truss: uz of every joint held; ``supports``: node -> tuple of held in-plane translations."""
    return {n: {"dofs": {d: "fixed" for d in tuple(supports.get(n, ())) + ("uz",)}} for n in nodes}


def sec(a, iz, iy=None, j=None, **kw):
    s = {"type": "properties", "A": float(a), "Iz": float(iz), "Iy": float(iy if iy is not None else iz),
         "J": float(j if j is not None else iz)}
    s.update(kw)
    return s


def bar(a):
    return {"type": "properties", "description": "axial-force member: Iy, Iz, J are placeholders, not used",
            "A": float(a), "Iy": 1.0, "Iz": 1.0, "J": 1.0}


def frame(i, j, s="S", **kw):
    m = {"i": i, "j": j, "section": s, "material": "M"}
    m.update(kw)
    return m


def truss(i, j, s):
    return {"type": "truss", "i": i, "j": j, "section": s, "material": "M"}


def one(loads, kind=None, extra=None):
    c = {"factors": {"L": 1.0}}
    if kind:
        c["analysisType"] = kind
    return {"L": {"loads": loads}}, {"C": c}


def add_vendor(lib, source_fn, pid, title, page, cites, analysis, mdl, recomputed, printed_ok=(), printed_bad=(),
               diagnosis=None, how="", limits="", whole_diagnosis=None):
    """Registers ``<pid>`` (active, or quarantined as a whole when ``whole_diagnosis`` is given: the engine
    cannot solve it or another tier fails) and, when ``printed_bad`` is given, ``<pid>-vendor-printed``
    (quarantined)."""
    if bool(printed_bad) != bool(diagnosis):
        raise ValueError(pid + ": printed_bad targets and their diagnosis go together")
    text = (analysis + " Reference values: computed here from this statement (" + DERIVED
            + (", " + how if how else "") + "); the textbook cited by the vendor page was not opened and no number "
              "is taken from it.")
    lib.add(pid, title, source_fn(page, cites), text, mdl, list(recomputed) + list(printed_ok),
            basis="recomputed-from-published-statement", limits=limits, diagnosis=whole_diagnosis)
    if printed_bad:
        lib.add(pid + "-vendor-printed", title + " - values printed by the vendor that do not reproduce",
                source_fn(page, cites),
                analysis + " Same model as problem '" + pid + "' (whose recomputed targets are active); the targets "
                "here are the vendor's printed values only.", mdl, list(printed_bad), basis="published",
                diagnosis=diagnosis, limits=limits)
