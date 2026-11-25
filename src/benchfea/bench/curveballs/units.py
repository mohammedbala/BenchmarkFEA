"""Curve ball: units relabelled without converting the numbers."""
from __future__ import annotations

from ...units import FORCE_UNITS, LENGTH_UNITS
from .registry import Candidate, register

CATEGORY = "units"


@register(
    name="units_relabelled", category=CATEGORY, title="units relabelled without conversion",
    changes="Changes the declared force unit or length unit and keeps every number (`/units/force` or "
            "`/units/length`).",
    why="Every dimensional number of the file is read in the declared units, so relabelling them describes a "
        "different physical structure: lengths, sections, loads, stiffnesses and strengths all scale by powers "
        "of the unit factor. Because the unit system is consistent the model stays self-similar: the NUMBERS of "
        "the results in file units and all design ratios are unchanged -- nothing inside the file looks wrong. "
        "The error shows only against physical data stated outside the file (a span given in metres, a load in "
        "kN, a deflection limit in millimetres): every physical result is off by the unit factor.",
    moves="Every physical quantity: lengths and displacements by the length factor, forces and reactions by the "
          "force factor, moments by their product. Design ratios and numbers in file units do not move.",
    effect_classes=("changes-analysis",),
    precondition="any model (the declared units are part of every file)")
def units_relabelled(ctx):
    units = ctx.doc["units"]
    out = []
    for key, table in (("length", LENGTH_UNITS), ("force", FORCE_UNITS)):
        for new in table:
            if new == units[key]:
                continue
            factor = float(table[new] / table[units[key]])
            out.append(Candidate(f"{key} unit {units[key]} -> {new} (numbers kept)",
                                 lambda d, key=key, new=new: d["units"].__setitem__(key, new), {"units": [key]},
                                 {"relabelled": key, "from": units[key], "to": new, "physicalFactor": factor}))
    return out
