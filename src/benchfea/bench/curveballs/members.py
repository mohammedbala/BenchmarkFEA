"""Curve balls on member orientation and sections."""
from __future__ import annotations

from ..families._common import CATALOGUE, UNIT_SYSTEMS, section_definition
from ._util import tidy
from .registry import Candidate, register

CATEGORY = "members"


@register(
    name="member_rolled_90", category=CATEGORY, title="member rolled 90 degrees (weak axis)",
    changes="Adds 90 degrees to the roll angle of one frame member whose section has different stiffness about "
            "its two axes (`/members/<id>/orientation/roll`).",
    why="The roll angle decides which principal axis of the section resists the bending. Rolled by 90 degrees an "
        "I-shape bends about its weak axis: for typical beam shapes the stiffness drops by a factor of 10 to 30 "
        "and the flexural strength by 3 to 6; a column turned the wrong way buckles and sways about the wrong "
        "axis. Loads given in member local axes turn with the member.",
    moves="Deflections and drift (much larger), distribution of moments in an indeterminate frame, the axis "
          "about which the member's moments act (Mz becomes My), flexure and compression design ratios.",
    effect_classes=("changes-analysis",),
    criteria=("analysisChange", "designPassFailChanged", "designRatioChange"),      # the resistance changes
    precondition="a frame member whose section has Iy different from Iz by more than 5 %")
def member_rolled_90(ctx):
    out = []
    for mid, m in ctx.model.members.items():
        p = ctx.model.sections[m.section].props
        if m.type != "frame" or abs(p.Iy - p.Iz) <= 0.05 * max(p.Iy, p.Iz):
            continue
        roll = (ctx.doc["members"][mid].get("orientation") or {}).get("roll", 0)
        new = tidy((roll + 90.0) % 180.0)

        def edit(d, mid=mid, new=new):
            d["members"][mid].setdefault("orientation", {})["roll"] = new
        out.append(Candidate(f"member {mid}: roll {roll:g} -> {new:g} degrees", edit, {"members": [mid]},
                             {"baseRoll": roll, "newRoll": new}))
    return out


def _catalogue_position(doc, sid):
    """(unit system, kind, index) of a section that is an entry of the generators' generic catalogue."""
    for us_name, us in UNIT_SYSTEMS.items():
        if us["units"] != doc["units"]:
            continue
        table = CATALOGUE["metric" if us["metric"] else "imperial"]
        for kind, entries in table.items():
            for i in range(len(entries)):
                name, definition, _ = section_definition(kind, i, us_name)
                if name == sid and definition == doc["sections"][sid]:
                    return us_name, kind, i
    return None


def _section_swaps(ctx, direction: int):
    """Candidates replacing the section of one member by a lighter (-1) or heavier (+1) one of the
    same kind: a section of the model, or the neighbours of a generic catalogue section."""
    doc, model = ctx.doc, ctx.model
    out = []
    cache = {}
    for mid, m in model.members.items():
        sid = m.section
        a0, kind = model.sections[sid].props.A, model.sections[sid].kind
        options = []                                     # (new id, definition or None)
        for other, s in model.sections.items():
            if other != sid and s.kind == kind and (s.props.A - a0) * direction > 1e-9 * a0:
                options.append((other, None))
        if sid not in cache:
            cache[sid] = _catalogue_position(doc, sid)
        pos = cache[sid]
        if pos is not None:
            us_name, ckind, i = pos
            n = len(CATALOGUE["metric" if UNIT_SYSTEMS[us_name]["metric"] else "imperial"][ckind])
            for step in (1, 2):
                j = i + direction * step
                if 0 <= j < n:
                    name, definition, _ = section_definition(ckind, j, us_name)
                    if name not in doc["sections"]:
                        options.append((name, definition))
                    elif (name, None) not in options and doc["sections"][name] == definition:
                        options.append((name, None))
        for new, definition in dict(options).items():
            def edit(d, mid=mid, new=new, definition=definition):
                if definition is not None:
                    d["sections"][new] = dict(definition)
                d["members"][mid]["section"] = new
            out.append(Candidate(f"member {mid}: section {sid} -> {new}", edit,
                                 {"members": [mid], "sections": [sid, new]}, {"baseSection": sid, "newSection": new}))
    return out


@register(
    name="section_lighter", category=CATEGORY, title="section swapped for a lighter one",
    changes="Gives one member a lighter section of the same kind -- another section of the model or the next "
            "smaller sizes of the generators' generic catalogue (`/members/<id>/section`, and `/sections/<new>` "
            "if the section is new to the file).",
    why="A lighter section has less area, inertia and strength. The member's own design ratios rise (capacity "
        "falls faster than the demand), and in an indeterminate structure the softer member sheds force to its "
        "neighbours, so their demands change too; deflections grow.",
    moves="Design ratios and pass / fail of the member, deflections and drift, distribution of moments between "
          "members, self weight.",
    effect_classes=("changes-analysis",),
    criteria=("analysisChange", "designPassFailChanged", "designRatioChange"),      # the resistance changes
    precondition="a member for which a lighter section of the same kind exists (in the model or in the generic "
                 "catalogue of the generators)")
def section_lighter(ctx):
    return _section_swaps(ctx, -1)


@register(
    name="section_heavier", category=CATEGORY, title="section swapped for a heavier one",
    changes="Gives one member a heavier section of the same kind -- another section of the model or the next "
            "larger sizes of the generators' generic catalogue (`/members/<id>/section`, and `/sections/<new>` "
            "if the section is new to the file).",
    why="A heavier section looks safe, but stiffness attracts force: in an indeterminate structure the stiffer "
        "member takes a larger share of the moment, so the demands on the joints and on the neighbouring members "
        "change, and the member's own ratios drop -- a model that reports a comfortable ratio for a section that "
        "is not the one built.",
    moves="Design ratios of the member (down) and of its neighbours, distribution of moments, deflections and "
          "drift (smaller), self weight.",
    effect_classes=("changes-analysis",),
    criteria=("analysisChange", "designPassFailChanged", "designRatioChange"),      # the resistance changes
    precondition="a member for which a heavier section of the same kind exists (in the model or in the generic "
                 "catalogue of the generators)")
def section_heavier(ctx):
    return _section_swaps(ctx, +1)


@register(
    name="sections_swapped", category=CATEGORY, title="sections swapped between two members",
    changes="Exchanges the sections of two members that have different sections (`/members/<a>/section`, "
            "`/members/<b>/section`).",
    why="Each member was sized for its own forces. With the sections exchanged (a beam section on a column, the "
        "upper-storey column size at the base) one member is under-sized and the other over-sized; the stiffness "
        "distribution of the frame changes with them, and so do the forces every member attracts.",
    moves="Design ratios and pass / fail of both members, moment distribution, deflections and drift.",
    effect_classes=("changes-analysis",),
    criteria=("analysisChange", "designPassFailChanged", "designRatioChange"),      # the resistance changes
    precondition="two members with different sections")
def sections_swapped(ctx):
    ids = list(ctx.model.members)
    out = []
    for k, a in enumerate(ids):
        sa = ctx.doc["members"][a]["section"]
        for b in ids[k + 1:] + ids[:k]:
            sb = ctx.doc["members"][b]["section"]
            if sb != sa:
                def edit(d, a=a, b=b, sa=sa, sb=sb):
                    d["members"][a]["section"], d["members"][b]["section"] = sb, sa
                out.append(Candidate(f"members {a} ({sa}) and {b} ({sb})", edit,
                                     {"members": [a, b], "sections": [sa, sb]}))
                break
    return out
