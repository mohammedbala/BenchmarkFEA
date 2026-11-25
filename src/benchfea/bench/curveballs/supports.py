"""Curve balls on supports."""
from __future__ import annotations

from ._util import ROTATIONS, TRANSLATIONS, plane_dofs, set_support, support_dofs, tidy
from .registry import Candidate, register

CATEGORY = "supports"
_PLANE = (" (DOFs that are fixed at every node of the model -- the out-of-plane restraint of a plane frame "
          "modelled in 3D space -- are never touched.)")


def _tables(ctx):
    """[(node, dof table, plane DOFs)] of the supports that restrain something besides the plane DOFs."""
    plane = plane_dofs(ctx.doc)
    out = []
    for node, s in ctx.doc.get("supports", {}).items():
        t = support_dofs(s, ctx.vertical)
        if any(t[d] != "free" for d in t if d not in plane):
            out.append((node, t, plane))
    return out


def _candidate(label, node, table, info=None):
    return Candidate(label, lambda d, node=node, table=dict(table): set_support(d, node, table),
                     {"nodes": [node], "supports": [node]}, info or {})


@register(
    name="support_fixed_to_pinned", category=CATEGORY, title="fixed (or partially fixed) support made pinned",
    changes="Frees the rotational restraints of one support -- fixed rotations and rotational springs; all of "
            "them, or all but the rotation about the vertical axis (`/supports/<node>`)." + _PLANE,
    why="A fixed base takes moment and restrains the column end; a pinned base does not. Pinning it removes the "
        "base moment, roughly doubles the effective length of a sway column, and moves the moment to the other "
        "end of the member and into the beams; the frame becomes markedly more flexible.",
    moves="Reaction moment at the support (to zero), column end moments, drift and deflections (larger), "
          "moments in the connected beams; second-order amplification.",
    effect_classes=("changes-analysis",),
    precondition="a support with a fixed or sprung rotation (other than the plane restraint of a plane frame)")
def support_fixed_to_pinned(ctx):
    out = []
    twist = {"X": "rx", "Y": "ry", "Z": "rz"}[ctx.vertical]
    for node, t, plane in _tables(ctx):
        rot = [d for d in ROTATIONS if d not in plane and t[d] != "free"]
        if rot:
            out.append(_candidate(f"support {node}: rotation(s) {', '.join(rot)} freed", node,
                                  {**t, **{d: "free" for d in rot}}))
        # variant: the rotation about the vertical axis (twist of a column) stays held
        bending = [d for d in rot if d != twist]
        if bending and bending != rot:
            out.append(_candidate(f"support {node}: rotation(s) {', '.join(bending)} freed ({twist} kept)", node,
                                  {**t, **{d: "free" for d in bending}}))
    return out


@register(
    name="support_pinned_to_fixed", category=CATEGORY, title="pinned (or partially fixed) support made fixed",
    changes="Fixes the free or sprung rotations of one support whose three translations are restrained "
            "(`/supports/<node>`)." + _PLANE,
    why="Full fixity is rarely there unless the foundation is designed for it. A base modelled as fixed takes a "
        "moment that relieves the rest of the frame: drift and the moments at the top of the columns are under-"
        "estimated, and the foundation receives a moment it may not be able to take.",
    moves="Reaction moment at the support (appears), column end moments, drift and deflections (smaller), beam "
          "moments.",
    effect_classes=("changes-analysis",),
    precondition="a support with the three translations restrained and a free or sprung rotation")
def support_pinned_to_fixed(ctx):
    out = []
    for node, t, plane in _tables(ctx):
        if not all(t[d] != "free" for d in TRANSLATIONS):
            continue
        rot = [d for d in ROTATIONS if t[d] != "fixed"]
        if rot:
            out.append(_candidate(f"support {node}: rotation(s) {', '.join(rot)} fixed", node,
                                  {**t, **{d: "fixed" for d in rot}}))
    return out


@register(
    name="support_roller_direction", category=CATEGORY, title="roller turned to another direction",
    changes="Moves the restraint of a roller support from the translation it holds to one it leaves free "
            "(`/supports/<node>`)." + _PLANE,
    why="A roller carries load in one direction and slides in the others. Turned by 90 degrees it no longer "
        "carries the load it was there for (the span loses a support) and instead restrains a movement that "
        "should be free (thermal or axial movement, horizontal thrust). The structure may become a mechanism.",
    moves="Reactions (the turned support carries none of the original load), span moments and deflections of the "
          "members that lose the support, axial forces from the new restraint; or a mechanism.",
    effect_classes=("changes-analysis", "creates-mechanism"),
    precondition="a support that restrains some but not all translations (a roller), other than by the plane "
                 "restraint")
def support_roller_direction(ctx):
    out = []
    for node, t, plane in _tables(ctx):
        held = [d for d in TRANSLATIONS if d not in plane and t[d] != "free"]
        free = [d for d in TRANSLATIONS if d not in plane and t[d] == "free"]
        for a in held:
            for b in free:
                out.append(_candidate(f"support {node}: restraint moved from {a} to {b}", node,
                                      {**t, a: "free", b: t[a]}))
    return out


@register(
    name="support_spring_stiffness", category=CATEGORY, title="support spring stiffness changed",
    changes="Multiplies one support spring stiffness by 0.01, 0.1, 10 or 100 (`/supports/<node>/dofs/<dof>/spring`).",
    why="A spring support shares load with the structure in proportion to stiffness. A spring ten times too "
        "stiff behaves almost like a rigid support, one ten times too soft almost like no support; a slip of "
        "the decimal point or of the unit (kN/m against kN/mm) therefore redistributes reactions and moments.",
    moves="Reaction at the spring, displacement of the supported node, moments and shears of the members that "
          "span over or frame into it, drift for rotational base springs.",
    effect_classes=("changes-analysis",),
    precondition="a support with a spring")
def support_spring_stiffness(ctx):
    out = []
    for node, t, _ in _tables(ctx):
        for dof, v in t.items():
            if isinstance(v, dict):
                for factor in (0.01, 0.1, 10.0, 100.0):
                    new = tidy(v["spring"] * factor)
                    out.append(_candidate(f"support {node}: spring {dof} {v['spring']:g} -> {new:g}", node,
                                          {**t, dof: {"spring": new}},
                                          {"factor": factor, "baseValue": v["spring"], "newValue": new}))
    return out


@register(
    name="support_removed", category=CATEGORY, title="support removed",
    changes="Deletes one support (`/supports/<node>`)." + _PLANE,
    why="Every support is part of the load path. Without it the load it carried must reach the remaining "
        "supports through the members -- spans double, cantilevers appear, moments and deflections multiply -- "
        "or, if the structure was not redundant enough, there is no equilibrium at all (a mechanism).",
    moves="All reactions, moments and deflections of the members that lose the support (often several times "
          "larger); or no results: mechanism diagnosis (rigid-body mode).",
    effect_classes=("changes-analysis", "creates-mechanism"),
    precondition="a support (other than the plane restraint of a plane frame)")
def support_removed(ctx):
    out = []
    for node, t, plane in _tables(ctx):
        out.append(_candidate(f"support {node}", node, {d: (t[d] if d in plane else "free") for d in t}))
    return out
