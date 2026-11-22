"""Resolution of design variables: design member > group(s) > design.defaults >
built-in default, with the provenance of every value.

DESIGN ONLY.  Nothing computed here is read by the analysis: the solver works
from ``Model.members`` etc. and never from ``Model.design_members``.

The same function serves validation (raw file units, collecting problems) and
the model builder (document already converted to internal units).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from . import defaults as D
from .geometry import is_vertical
from .problems import ERROR, WARNING, pointer

#: unbraced length -> bracing lists that define the same thing
PAIRS = {"Lx": ("major",), "Ly": ("minor",), "Lz": ("torsion",), "Lb": ("topFlange", "bottomFlange")}
BRACE_KEYS = {"major": "braceMajor", "minor": "braceMinor", "torsion": "braceTorsion",
              "topFlange": "braceTopFlange", "bottomFlange": "braceBottomFlange"}
SIMPLE = ("code", "method", "role", "Kx", "Ky", "Kz", "Cb", "Cmx", "Cmy", "sideswayInhibited", "An", "U",
          "maxSlendernessCompression", "maxSlendernessTension", "deflectionLimits")
#: every resolved variable name, in report order
NAMES = SIMPLE[:3] + ("Lx", "Ly", "Lz", "Lb") + SIMPLE[3:6] + tuple(BRACE_KEYS.values()) + SIMPLE[6:]
#: dimension of variables that are not dimensionless (for the defaults-used report)
DIMENSIONS = {"Lx": "length", "Ly": "length", "Lz": "length", "Lb": "length", "An": "length^2"}

COLLINEAR_SINE = 1.0e-6


@dataclass(frozen=True)
class Resolved:
    """A resolved design variable.

    value   the value (None = to be derived by the design module, see `source`)
    source  "member" | "group:<id>" | "global" | "default" | "bracing"
            ("bracing": an unbraced length that is defined by brace points)
    path    JSON pointer where the value is set (or would be set, for defaults)
    reason  for defaults: why this default applies
    """

    value: Any
    source: str
    path: str
    reason: str = ""


def design_member_table(data: dict, geoms: dict, add: Callable) -> dict[str, dict]:
    """Design member id -> {"members": [...], "length": L | None, "explicit": bool}."""
    design = data.get("design") or {}
    members = data.get("members") or {}
    explicit = design.get("members") or {}
    table: dict[str, dict] = {}
    used: dict[str, str] = {}
    for dmid, dm in explicit.items():
        base = pointer("design", "members", dmid)
        given = isinstance(dm, dict) and "members" in dm
        mids = list(dm["members"]) if given else [dmid]
        ok = True
        for k, mid in enumerate(mids):
            if mid not in members:
                ok = False
                if given:
                    add("ref.unknown-member", f"{base}/members/{k}", f"design member {dmid!r} refers to unknown member {mid!r}",
                        "Use an id from /members.")
                else:
                    add("ref.unknown-member", base, f"design member {dmid!r} has no `members` list and no analysis "
                        f"member has the id {dmid!r}", "Add `members: [...]` or use the id of an analysis member.")
            elif mid in used:
                ok = False
                add("design.member-reused", f"{base}/members/{k}" if given else base,
                    f"member {mid!r} already belongs to design member {used[mid]!r}",
                    "An analysis member can be part of only one design member.")
            else:
                used[mid] = dmid
        if dmid in members and mids != [dmid]:
            ok = False
            add("design.id-collision", base, f"design member id {dmid!r} is also the id of an analysis member "
                "that is not its only member", "Rename the design member.")
        length = None
        if ok and all(m in geoms for m in mids):
            first, last = geoms[mids[0]], geoms[mids[-1]]
            for k in range(1, len(mids)):
                a, b = members[mids[k - 1]], members[mids[k]]
                ga, gb = geoms[mids[k - 1]], geoms[mids[k]]
                p = f"{base}/members/{k}"
                if a["j"] != b["i"]:
                    ok = False
                    add("design.not-connected", p, f"member {mids[k]!r} does not start at the end node of "
                        f"{mids[k - 1]!r}", "List the members in order, all running the same way (end j of one = "
                        "end i of the next).")
                    continue
                if np.linalg.norm(np.cross(ga[2][0], gb[2][0])) > COLLINEAR_SINE or ga[2][0] @ gb[2][0] < 0:
                    ok = False
                    add("design.not-collinear", p, f"member {mids[k]!r} is not collinear with {mids[k - 1]!r}",
                        "A design member must be straight; split it into separate design members.")
                elif np.linalg.norm(ga[2][1] - gb[2][1]) > COLLINEAR_SINE:
                    ok = False
                    add("design.orientation-differs", p, f"member {mids[k]!r} has a different local y axis than "
                        f"{mids[k - 1]!r}", "Give all members of a design member the same orientation.")
                if any(a.get(key, dflt) != b.get(key, dflt)
                       for key, dflt in (("section", None), ("material", None), ("type", "frame"))):
                    ok = False
                    add("design.not-prismatic", p, f"member {mids[k]!r} differs from {mids[k - 1]!r} in section, "
                        "material or type", "A design member must be prismatic and of one material and type.")
            if ok:
                length = float(np.linalg.norm(last[1] - first[0]))
        table[dmid] = {"members": mids, "length": length, "explicit": True, "valid": ok}
    if design.get("autoMembers", True):
        for mid in members:
            if mid not in used and mid not in table:
                g = geoms.get(mid)
                table[mid] = {"members": [mid], "length": g[3] if g else None, "explicit": False,
                              "valid": g is not None}
    return table


def resolve(data: dict, geoms: dict, add: Callable) -> dict[str, dict]:
    """Resolve all design variables.

    ``geoms``: member id -> (start, end, R, length, ...) from
    :func:`benchfea.geometry.member_geometry`.
    ``add(code, path, message, hint, severity=ERROR)`` receives problems.

    Returns design member id -> {"members", "length", "values": {name: Resolved}}.
    ``An`` (default: gross area) is left as ``Resolved(None, "default", ...)``
    for the model builder, which knows the section.
    """
    design = data.get("design") or {}
    members = data.get("members") or {}
    combos = data.get("combinations") or {}
    vertical = (data.get("conventions") or {}).get("verticalAxis", "Z")
    if vertical not in ("Y", "Z"):
        vertical = "Z"
    table = design_member_table(data, geoms, add)
    glob = design.get("defaults") or {}
    groups = design.get("groups") or {}
    explicit = design.get("members") or {}

    # references made inside any level
    def scan_limits(level: dict, base: str):
        for k, lim in enumerate(level.get("deflectionLimits") or ()):
            for q, cid in enumerate(lim.get("combinations") or ()):
                if cid not in combos:
                    add("ref.unknown-combination", f"{base}/deflectionLimits/{k}/combinations/{q}",
                        f"deflection limit refers to unknown combination {cid!r}", "Use an id from /combinations.")

    scan_limits(glob, "/design/defaults")
    membership: dict[str, list[str]] = {}
    for gid, grp in groups.items():
        gbase = pointer("design", "groups", gid)
        scan_limits(grp, gbase)
        for k, dmid in enumerate(grp.get("members") or ()):
            if dmid not in table:
                owner = next((o for o, t in table.items() if dmid in t["members"]), None)
                if owner:
                    add("ref.unknown-design-member", f"{gbase}/members/{k}",
                        f"design group {gid!r} refers to {dmid!r}, which is not a design member: it is part of "
                        f"design member {owner!r}", f"Refer to {owner!r} instead.")
                    continue
                add("ref.unknown-design-member", f"{gbase}/members/{k}",
                    f"design group {gid!r} refers to unknown design member {dmid!r}",
                    "Use the id of a design member (/design/members) or of an analysis member"
                    + ("" if design.get("autoMembers", True) else " (autoMembers is false, so analysis members "
                                                                    "are not design members automatically)") + ".")
            else:
                membership.setdefault(dmid, []).append(gid)
    for dmid, dm in explicit.items():
        scan_limits(dm, pointer("design", "members", dmid))

    def get(level: dict, key: str):
        if "/" in key:
            a, b = key.split("/")
            sub = level.get(a) or {}
            return (True, sub[b]) if b in sub else (False, None)
        return (True, level[key]) if key in level else (False, None)

    out: dict[str, dict] = {}
    for dmid, info in table.items():
        own = explicit.get(dmid) or {}
        base = pointer("design", "members", dmid)
        L = info["length"]

        def find(key: str):
            """-> (value, source, path, rank, level dict) or None"""
            has, val = get(own, key)
            if has:
                return val, "member", f"{base}/{key}", (0, 0), own
            cands = []
            for gid in membership.get(dmid, ()):
                has, val = get(groups[gid], key)
                if has:
                    cands.append(((1, -int(groups[gid].get("priority", 0))), gid, val))
            if cands:
                cands.sort(key=lambda c: c[0])
                rank, gid, val = cands[0]
                for r2, g2, v2 in cands[1:]:
                    if r2 == rank and v2 != val:
                        add("design.group-conflict", f"{pointer('design', 'groups', g2)}/{key}",
                            f"groups {gid!r} and {g2!r} both set {key.replace('/', '.')} for design member "
                            f"{dmid!r} with different values and equal priority",
                            "Give one group a higher `priority`, or set the value on the design member.")
                return val, f"group:{gid}", f"{pointer('design', 'groups', gid)}/{key}", rank, groups[gid]
            has, val = get(glob, key)
            if has:
                return val, "global", f"/design/defaults/{key}", (2, 0), glob
            return None

        values: dict[str, Resolved] = {}
        for name in SIMPLE:
            hit = find(name)
            if hit:
                values[name] = Resolved(hit[0], hit[1], hit[2])
            elif name in D.DESIGN:
                val, reason = D.DESIGN[name]
                values[name] = Resolved(val, "default", f"{base}/{name}", reason)
            elif name == "role":
                mtypes = {members[m].get("type", "frame") for m in info["members"] if m in members}
                g = geoms.get(info["members"][0])
                if mtypes - {"frame"}:
                    role, why = "brace", "truss-type member"
                elif g is not None and is_vertical(g[2][0], vertical):
                    role, why = "column", "frame member parallel to the vertical axis"
                else:
                    role, why = "beam", "frame member not parallel to the vertical axis"
                values[name] = Resolved(role, "default", f"{base}/role", f"role from geometry: {why}")
            elif name == "An":
                values[name] = Resolved(None, "default", f"{base}/An", "net area = gross area (no holes)")

        for lname, subs in PAIRS.items():
            lhit = find(lname)
            bhits = {s: find(f"bracing/{s}") for s in subs}
            branks = [h[3] for h in bhits.values() if h]
            if lhit and branks:
                if lhit[3] == min(branks):
                    add("design.length-and-bracing", lhit[2],
                        f"{lname} and bracing.{'/'.join(subs)} are both given at the same level for design "
                        f"member {dmid!r}", f"Give either the unbraced length {lname} or the brace locations.")
                if lhit[3] <= min(branks):
                    bhits = {s: None for s in subs}       # explicit length wins over lower-level bracing
                else:
                    lhit = None
            for s, h in bhits.items():
                key = BRACE_KEYS[s]
                if not h:
                    values[key] = Resolved((), "default", f"{base}/bracing/{s}",
                                           "no intermediate brace points (ends only)")
                    continue
                absolute = (h[4].get("bracing") or {}).get("positionType", "relative") == "absolute"
                fr = []
                for k, p in enumerate(h[0]):
                    lim = L if absolute else 1.0
                    if lim is not None and not 0.0 < p < lim:
                        add("design.bracing-range", f"{h[2]}/{k}",
                            f"brace location {p:g} is not strictly inside the design member "
                            f"(0 < position < {lim:g})", "Brace locations are intermediate points; the member "
                            "ends are braced by definition. Check positionType (relative/absolute).")
                    elif lim is not None:
                        fr.append(p / lim)
                values[key] = Resolved(tuple(sorted(fr)), h[1], h[2])
            if lhit:
                values[lname] = Resolved(lhit[0], lhit[1], lhit[2])
            elif any(bhits.values()):
                values[lname] = Resolved(None, "bracing", f"{base}/{lname}",
                                         "unbraced lengths are the segments between brace points")
            else:
                values[lname] = Resolved(L, "default", f"{base}/{lname}",
                                         "unbraced length = design member length (no brace points given)")
        out[dmid] = {"members": info["members"], "length": L, "valid": info["valid"],
                     "values": {n: values[n] for n in NAMES}}
    return out
