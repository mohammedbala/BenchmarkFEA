"""Curve balls on design variables (design-only: the analysis must stay bitwise identical)."""
from __future__ import annotations

from ._util import is_number, tidy
from .record import split_pointer
from .registry import Candidate, register

CATEGORY = "design"
_ONLY = (" Nothing in the `design` block can change the analysis: member forces, displacements and reactions "
         "are identical; only the available strengths, and with them the design ratios and pass / fail, move.")
_BRACE = {"Lx": "braceMajor", "Ly": "braceMinor", "Lz": "braceTorsion"}
_FLANGES = ("topFlange", "bottomFlange")


def _set(doc: dict, path: str, value) -> None:
    parts = split_pointer(path)
    node = doc
    for key in parts[:-1]:
        node = node.setdefault(key, {})
    node[parts[-1]] = value


def _block(doc: dict, path: str):
    node = doc
    for key in split_pointer(path):
        node = node.get(key) if isinstance(node, dict) else None
        if node is None:
            return None
    return node


def _round(x: float):
    return tidy(float(f"{x:.4g}"))


def _longest_segment(points, length: float) -> float:
    xs = [0.0] + sorted(points) + [length]
    return max(b - a for a, b in zip(xs, xs[1:]))


def _absolute(p, name: str) -> list:
    """Brace positions of a resolved list as distances from the design member start (internal units).
    The model stores them as fractions of the length."""
    return [x * p.length for x in p.values[name].value]


def _users(model, name: str, path: str) -> list:
    """Design members whose variable ``name`` is taken from ``path``."""
    return [k for k, p in model.design_members.items() if p.values[name].path == path]


def _steel_members(model):
    """Design members that the checks can say something about (they have a governing ratio)."""
    return model.design_members.items()


@register(
    name="design_unbraced_length_changed", category=CATEGORY, title="unbraced length changed",
    changes="Gives one design member an explicit unbraced length -- Lx or Ly (flexural buckling) or Lb (lateral-"
            "torsional buckling) -- that is half or twice the length the bracing really gives "
            "(`/design/members/<id>/L*`).",
    why="The available compression strength depends on (K L / r) squared in the elastic range and the lateral-"
        "torsional strength on Lb: a wrong unbraced length changes the capacity, not the demand. Half the length "
        "is unconservative (the member is credited with bracing it does not have), twice is conservative." + _ONLY,
    moves="Available compression strength (E3) or flexural strength (F2 lateral-torsional buckling), the H1 "
          "interaction ratio, pass / fail of the member. No analysis quantity.",
    effect_classes=("design-only",),
    precondition="a `design` block and a design member")
def design_unbraced_length_changed(ctx):
    model, u = ctx.model, ctx.model.units
    out = []
    for dm, p in _steel_members(model):
        member_block = _block(ctx.doc, f"/design/members/{dm}") or {}
        for var in ("Lx", "Ly", "Lb"):
            r = p.values[var]
            if r.source == "member" or var in member_block:
                current = r.value
            elif r.value is not None:
                current = r.value
            else:                                          # defined by brace points: the longest segment
                names = ("braceTopFlange", "braceBottomFlange") if var == "Lb" else (_BRACE[var],)
                current = max(_longest_segment(_absolute(p, n), p.length) for n in names)
            lists = (member_block.get("bracing") or {})
            if (var == "Lb" and any(f in lists for f in _FLANGES)) or \
                    (var != "Lb" and {"Lx": "major", "Ly": "minor"}[var] in lists):
                continue                                   # same level as its bracing list: an input error
            for factor in (0.5, 2.0):
                new = _round(float(u.from_internal(current * factor, "length")))
                out.append(Candidate(
                    f"design member {dm}: {var} = {new} ({factor:g} x the longest unbraced segment of the base model)",
                    lambda d, dm=dm, var=var, new=new: _set(d, f"/design/members/{dm}/{var}", new),
                    {"designMembers": [dm], "members": list(p.members)},
                    {"variable": var, "factor": factor,
                     "baseValue": _round(float(u.from_internal(current, "length"))), "newValue": new}))
    return out


@register(
    name="design_k_changed", category=CATEGORY, title="effective length factor changed",
    changes="Halves or doubles an effective length factor Kx or Ky where it is set (group, defaults) or, if it "
            "is a built-in default, sets it on one design member (`/design/.../Kx`, `/Ky`).",
    why="K turns the member length into the buckling length. K = 1 for a column of a sway frame analysed first "
        "order ignores that the frame can sway (K > 1); K = 2 where the analysis already includes second-order "
        "effects with K = 1 penalises the column twice. The compression strength changes roughly with 1 / K "
        "squared for slender members." + _ONLY,
    moves="Available compression strength (E3), the H1 interaction ratio of members with axial compression, "
          "pass / fail. No analysis quantity.",
    effect_classes=("design-only",),
    precondition="a `design` block and a design member")
def design_k_changed(ctx):
    model = ctx.model
    out, seen = [], set()
    for dm, p in _steel_members(model):
        for var in ("Kx", "Ky"):
            r = p.values[var]
            if r.path in seen or not is_number(r.value):
                continue
            seen.add(r.path)
            users = _users(model, var, r.path)
            for factor in (0.5, 2.0):
                new = tidy(r.value * factor)
                where = r.source if r.source != "default" else f"design member {dm} (built-in default {r.value:g})"
                out.append(Candidate(f"{var} {r.value:g} -> {new:g} at {where}",
                                     lambda d, path=r.path, new=new: _set(d, path, new),
                                     {"designMembers": users, "members": sorted({m for k in users for m in model.design_members[k].members})},
                                     {"variable": var, "factor": factor, "baseValue": r.value, "newValue": new}))
    return out


def _brace_lists(ctx):
    """{path: (name, tuple of points, users)} of the brace lists written in the file."""
    found = {}
    for dm, p in ctx.model.design_members.items():
        for name in ("braceMajor", "braceMinor", "braceTorsion", "braceTopFlange", "braceBottomFlange"):
            r = p.values[name]
            if r.source in ("default", "bracing") or not r.value:
                continue
            if _block(ctx.doc, r.path) is None:
                continue
            found.setdefault(r.path, (name, list(_block(ctx.doc, r.path)), []))[2].append(dm)
    return found


@register(
    name="design_bracing_points_changed", category=CATEGORY, title="bracing points changed",
    changes="Edits a list of brace points: all points but the one nearest midspan removed, one point removed, or "
            "a single point moved to half its distance (`/design/.../bracing/<list>`).",
    why="Brace points divide the member into unbraced segments; the longest segment governs flexural buckling "
        "(minor / major lists) and each segment has its own lateral-torsional strength (flange lists). Removing "
        "or moving a point lengthens a segment and lowers the capacity; the member forces are untouched." + _ONLY,
    moves="Available compression strength about the braced axis, lateral-torsional buckling strength and Cb of "
          "the segments, interaction ratio, pass / fail. No analysis quantity.",
    effect_classes=("design-only",),
    precondition="a brace-point list in the `design` block (bracing.major / minor / torsion / topFlange / "
                 "bottomFlange)")
def design_bracing_points_changed(ctx):
    out = []
    for path, (name, points, users) in _brace_lists(ctx).items():
        variants = []
        if len(points) >= 2:
            mid = (min(points) + max(points)) / 2.0
            keep = min(points, key=lambda x: (abs(x - mid), x))
            variants.append(("only the point nearest the middle kept", [keep]))
            for k, x in enumerate(points):
                variants.append((f"point {x:g} removed", points[:k] + points[k + 1:]))
        else:
            variants.append((f"point {points[0]:g} moved to {tidy(points[0] / 2.0):g}", [tidy(points[0] / 2.0)]))
        members = sorted({m for k in users for m in ctx.model.design_members[k].members})
        for label, new in variants:
            out.append(Candidate(f"{path}: {label}", lambda d, path=path, new=tuple(new): _set(d, path, list(new)),
                                 {"designMembers": list(users), "members": members}, {"list": name}))
    return out


@register(
    name="design_lb_full_where_braced", category=CATEGORY, title="Lb left at full length where the flange is braced",
    changes="Deletes the flange brace lists of a design block (`bracing.topFlange`, `bracing.bottomFlange`), so "
            "that Lb falls back to the full member length.",
    why="The lateral-torsional buckling strength of a beam falls quickly with the unbraced length Lb. A beam "
        "whose compression flange is held by the floor or by purlins but is checked with Lb = full span is "
        "heavily under-rated: the check fails or a much larger section is chosen. The moments do not change." + _ONLY,
    moves="Flexural strength for lateral-torsional buckling (F2-2 / F2-3 instead of the plastic moment), the "
          "flexure and interaction ratios (up), pass / fail. No analysis quantity.",
    effect_classes=("design-only",),
    precondition="flange brace points (bracing.topFlange / bottomFlange) in the `design` block")
def design_lb_full_where_braced(ctx):
    blocks = {}
    for path, (name, _, users) in _brace_lists(ctx).items():
        if name in ("braceTopFlange", "braceBottomFlange"):
            blocks.setdefault(path.rsplit("/", 1)[0], set()).update(users)
    out = []
    for bracing_path, users in blocks.items():
        def edit(d, bracing_path=bracing_path):
            b = _block(d, bracing_path)
            for f in _FLANGES:                             # the `bracing` object itself is kept
                b.pop(f, None)
        users = [k for k in ctx.model.design_members if k in users]
        out.append(Candidate(f"{bracing_path}: flange brace points deleted", edit,
                             {"designMembers": users,
                              "members": sorted({m for k in users for m in ctx.model.design_members[k].members})}))
    return out


@register(
    name="design_lb_short_where_unbraced", category=CATEGORY, title="Lb short where the flange is not braced",
    changes="Adds flange brace points (quarter points, third points or midspan, both flanges) to a design block "
            "whose members have no flange bracing, so that Lb becomes a fraction of the member length "
            "(`/design/.../bracing/topFlange`, `/bottomFlange`).",
    why="The reverse error: a member whose compression flange is free over its whole length is credited with "
        "bracing that does not exist. The lateral-torsional buckling strength is over-estimated -- an "
        "unconservative design that the analysis cannot reveal, because the forces are the same." + _ONLY,
    moves="Flexural strength for lateral-torsional buckling (up), flexure and interaction ratios (down), "
          "pass / fail. No analysis quantity.",
    effect_classes=("design-only",),
    precondition="a design member whose Lb is the full member length (no flange brace points and no Lb given)")
def design_lb_short_where_unbraced(ctx):
    model, doc = ctx.model, ctx.doc
    groups = (doc.get("design") or {}).get("groups") or {}
    blocks = {}
    for dm, p in model.design_members.items():
        if p.values["Lb"].source != "default" or p.values["braceTopFlange"].value or p.values["braceBottomFlange"].value:
            continue
        home = [g for g, b in groups.items() if dm in b["members"]]
        if home:
            g = max(home, key=lambda x: groups[x].get("priority", 0))
            path = f"/design/groups/{g}"
            if "Lb" in groups[g] or (groups[g].get("bracing") or {}).get("positionType") == "absolute":
                continue
        else:
            path = f"/design/members/{dm}"
        blocks.setdefault(path, []).append(dm)
    out = []
    patterns = (("quarter points", [0.25, 0.5, 0.75]), ("third points", [0.3333, 0.6667]), ("midspan", [0.5]))
    for path, users in blocks.items():
        if path.startswith("/design/groups/"):
            g = path.rsplit("/", 1)[1]
            if any(model.design_members[k].values["Lb"].source != "default" for k in groups[g]["members"]
                   if k in model.design_members):
                continue                                   # a member of the group has its own Lb: leave it
        for label, pts in patterns:
            def edit(d, path=path, pts=tuple(pts)):
                b = d
                for key in split_pointer(path):
                    b = b.setdefault(key, {})
                br = b.setdefault("bracing", {})
                for f in _FLANGES:
                    br[f] = list(pts)
            out.append(Candidate(f"{path}: both flanges braced at the {label}", edit,
                                 {"designMembers": list(users),
                                  "members": sorted({m for k in users for m in model.design_members[k].members})},
                                 {"pattern": label}))
    return out
