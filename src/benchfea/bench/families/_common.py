"""Shared tools of the model families: the stable random draw, unit systems with engineering
rounding, the generic section catalogue and the document builder.

Families describe a structure in ABSTRACT axes ``x, y`` (horizontal) and ``v`` (up).  The builder
maps them to the file's axes for the chosen vertical-axis convention by a proper rotation:

    verticalAxis "Z":  X = x,  Y = y,  Z = v
    verticalAxis "Y":  X = x,  Y = v,  Z = -y

so the same family produces the same physical structure in both conventions.

Nothing here is reference data: the section catalogue is a list of GENERIC round dimensions defined
in this file (not AISC shape data), the material is nominal structural steel, and loads are drawn
from ranges chosen for sensible proportions.
"""
from __future__ import annotations

import hashlib
import math

import numpy as np

__all__ = ["Draw", "stable_seed", "UNIT_SYSTEMS", "CATALOGUE", "section_definition", "Builder", "tidy",
           "GenerationRejected"]


class GenerationRejected(Exception):
    """A drawn parameter set cannot be turned into an acceptable model (the caller redraws)."""


# --------------------------------------------------------------------------------- random draws
def stable_seed(*parts) -> int:
    """128-bit seed from the text of ``parts`` (SHA-256; never Python's salted ``hash``)."""
    text = "benchfea.bench|" + "|".join(str(p) for p in parts)
    return int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest()[:16], "big")


class Draw:
    """Random draws from ``numpy.random.Generator(PCG64(seed))``.

    Only the raw 64-bit output of the bit generator is used (``random_raw``), whose stream numpy
    guarantees to be stable; every distribution below is plain integer / IEEE arithmetic on it, so
    the draws do not depend on the numpy version or the platform.
    """

    def __init__(self, seed: int):
        self.generator = np.random.Generator(np.random.PCG64(seed))
        self._raw = self.generator.bit_generator.random_raw

    def u(self) -> float:
        """Uniform in [0, 1) with 53 random bits."""
        return (int(self._raw()) >> 11) / 9007199254740992.0

    def integer(self, lo: int, hi: int) -> int:
        """Uniform integer in [lo, hi] (both included)."""
        if hi < lo:
            raise ValueError("empty range")
        return lo + min(int(self.u() * (hi - lo + 1)), hi - lo)

    def choice(self, options):
        options = list(options)
        return options[self.integer(0, len(options) - 1)]

    def uniform(self, lo: float, hi: float) -> float:
        return lo + (hi - lo) * self.u()

    def chance(self, p: float) -> bool:
        return self.u() < p

    def step(self, lo: float, hi: float, step: float) -> float:
        """A value of the grid lo, lo + step, ... <= hi."""
        n = int(math.floor((hi - lo) / step + 1e-9))
        return tidy(lo + step * self.integer(0, n))


def tidy(x):
    """Round away binary noise (12 significant digits) and write whole numbers as integers."""
    if isinstance(x, (list, tuple)):
        return [tidy(v) for v in x]
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return x
    if x == 0:
        return 0
    y = float(f"{float(x):.12g}")
    return int(y) if y == int(y) and abs(y) < 1e15 else y


def _grid(value: float, step: float) -> float:
    return tidy(round(value / step) * step)


def _sig(value: float, digits: int = 2) -> float:
    if value == 0:
        return 0
    return tidy(float(f"{value:.{digits}g}"))


# --------------------------------------------------------------------------------- unit systems
#: name -> units block, length of 1 m, force of 1 kN in file units, rounding grids, material (nominal
#: structural steel; E and Fy/Fu are model data of the generated file, not reference results).
UNIT_SYSTEMS = {
    "kN-m": {"units": {"force": "kN", "length": "m"}, "metric": True, "per_m": 1.0,
             "material": {"E": 200000000, "G": 77000000, "Fy": 345000, "Fu": 450000, "unitWeight": 77}},
    "kN-mm": {"units": {"force": "kN", "length": "mm"}, "metric": True, "per_m": 1000.0,
              "material": {"E": 200, "G": 77, "Fy": 0.345, "Fu": 0.45, "unitWeight": 7.7e-8}},
    "kip-in": {"units": {"force": "kip", "length": "in"}, "metric": False, "per_m": 1.0 / 0.0254,
               "material": {"E": 29000, "G": 11200, "Fy": 50, "Fu": 65, "unitWeight": 0.000284}},
}
KIP = 4.4482216152605          # kN per kip (exact factor used by benchfea.units)


# --------------------------------------------------------------------------------- section catalogue
# GENERIC shapes with round dimensions, defined here.  They are NOT rolled shapes of any standard and
# no dimension is taken from a shapes table.  Each list is ordered by increasing strength.
# I-shapes: (d, bf, tf, tw, r); boxes: (H, B, t); pipes: (D, t).  Metric in mm, imperial in in.
_METRIC = {
    "beam": [(100, 55, 6, 4, 5), (120, 64, 6, 4.5, 6), (150, 75, 7, 5, 7), (180, 90, 8, 5, 8), (200, 100, 8, 5, 8), (250, 125, 9, 6, 10), (300, 150, 10, 7, 12), (350, 175, 11, 7, 14),
             (400, 180, 13, 8, 16), (450, 200, 14, 9, 18), (500, 200, 16, 10, 20), (550, 210, 18, 11, 20),
             (600, 220, 20, 12, 22), (700, 250, 22, 13, 24), (800, 300, 25, 14, 26), (900, 300, 30, 16, 28),
             (1000, 350, 35, 18, 30)],
    "column": [(100, 100, 8, 6, 6), (120, 120, 9, 6, 8), (150, 150, 10, 7, 8), (200, 200, 12, 8, 12), (250, 250, 14, 9, 14), (300, 300, 16, 10, 16),
               (350, 350, 19, 12, 18), (400, 400, 22, 13, 20), (400, 400, 30, 18, 20), (450, 400, 38, 22, 22),
               (500, 450, 45, 26, 24), (550, 500, 55, 32, 26)],
    "box": [(40, 40, 3), (50, 50, 3), (60, 60, 4), (80, 80, 4), (100, 100, 5), (120, 120, 6), (150, 150, 6), (180, 180, 8), (200, 200, 10),
            (250, 250, 10), (300, 300, 12), (350, 350, 14), (400, 400, 16)],
    "pipe": [(21, 2.5), (27, 3), (34, 3), (42, 3), (48, 4), (60, 4), (76, 4), (89, 5), (114, 5), (140, 6), (168, 6), (194, 8), (219, 8), (273, 10),
             (324, 10), (356, 12), (406, 12)],
}
_IMPERIAL = {
    "beam": [(4, 2.2, 0.25, 0.15, 0.2), (5, 2.5, 0.25, 0.18, 0.25), (6, 3, 0.28, 0.2, 0.3), (7, 3.5, 0.3, 0.2, 0.3),
             (8, 4, 0.3, 0.2, 0.3), (10, 5, 0.35, 0.25, 0.4), (12, 6, 0.4, 0.25, 0.5), (14, 7, 0.45, 0.3, 0.5),
             (16, 7, 0.5, 0.3, 0.6), (18, 8, 0.55, 0.35, 0.7), (20, 8, 0.65, 0.4, 0.8), (22, 8.5, 0.7, 0.45, 0.8),
             (24, 9, 0.8, 0.5, 0.9), (28, 10, 0.9, 0.5, 1), (32, 12, 1, 0.55, 1), (36, 12, 1.2, 0.65, 1.1),
             (40, 14, 1.4, 0.7, 1.2)],
    "column": [(4, 4, 0.3, 0.2, 0.25), (5, 5, 0.35, 0.25, 0.3), (6, 6, 0.4, 0.25, 0.3), (8, 8, 0.45, 0.3, 0.5), (10, 10, 0.55, 0.35, 0.5), (12, 12, 0.65, 0.4, 0.6),
               (14, 14, 0.75, 0.45, 0.7), (16, 16, 0.85, 0.5, 0.8), (16, 16, 1.2, 0.7, 0.8), (18, 16, 1.5, 0.9, 0.9),
               (20, 18, 1.8, 1, 1), (22, 20, 2.2, 1.25, 1)],
    "box": [(1.5, 1.5, 0.12), (2, 2, 0.12), (2.5, 2.5, 0.15), (3, 3, 0.15), (4, 4, 0.2), (5, 5, 0.25), (6, 6, 0.25), (7, 7, 0.3), (8, 8, 0.4), (10, 10, 0.4),
            (12, 12, 0.5), (14, 14, 0.55), (16, 16, 0.65)],
    "pipe": [(0.8, 0.1), (1.05, 0.12), (1.3, 0.12), (1.7, 0.12), (1.9, 0.15), (2.4, 0.15), (3, 0.15), (3.5, 0.2), (4.5, 0.2), (5.5, 0.25), (6.6, 0.25), (7.6, 0.3), (8.6, 0.3),
             (10.8, 0.4), (12.8, 0.4), (14, 0.5), (16, 0.5)],
}
CATALOGUE = {"metric": _METRIC, "imperial": _IMPERIAL}


def _fmt(v) -> str:
    return f"{tidy(v)}"


def section_definition(kind: str, index: int, unit_system: str):
    """-> (generic name, section definition in file units, depth in file units) of catalogue entry."""
    us = UNIT_SYSTEMS[unit_system]
    dims = CATALOGUE["metric" if us["metric"] else "imperial"][kind][index]
    scale = 0.001 if unit_system == "kN-m" else 1.0
    val = [tidy(d * scale) for d in dims]
    label = "x".join(_fmt(d) for d in (dims[:-1] if kind in ("beam", "column") else dims))
    if kind in ("beam", "column"):
        return f"I-{label}", {"type": "I", "d": val[0], "bf": val[1], "tf": val[2], "tw": val[3], "r": val[4]}, val[0]
    if kind == "box":
        return f"Box-{label}", {"type": "box", "H": val[0], "B": val[1], "t": val[2]}, val[0]
    if kind == "pipe":
        return f"Pipe-{label}", {"type": "pipe", "D": val[0], "t": val[1]}, val[0]
    raise KeyError(kind)


def catalogue_size(kind: str, unit_system: str) -> int:
    return len(CATALOGUE["metric" if UNIT_SYSTEMS[unit_system]["metric"] else "imperial"][kind])


# --------------------------------------------------------------------------------- builder
_DOF = {"Z": {"ux": ("ux", 1), "uy": ("uy", 1), "uv": ("uz", 1), "rx": ("rx", 1), "ry": ("ry", 1), "rv": ("rz", 1)},
        "Y": {"ux": ("ux", 1), "uy": ("uz", -1), "uv": ("uy", 1), "rx": ("rx", 1), "ry": ("rz", -1), "rv": ("ry", 1)}}
_LOAD = {"Fx": "ux", "Fy": "uy", "Fv": "uv", "Mx": "rx", "My": "ry", "Mv": "rv"}
_FILE_LOAD = {"ux": "Fx", "uy": "Fy", "uz": "Fz", "rx": "Mx", "ry": "My", "rz": "Mz"}
_DIR = {"Z": {"x": ("X", 1), "y": ("Y", 1), "v": ("Z", 1)}, "Y": {"x": ("X", 1), "y": ("Z", -1), "v": ("Y", 1)}}
_ORDER = ("ux", "uy", "uz", "rx", "ry", "rz")
TRANSLATIONS = ("ux", "uy", "uv")
ALL_DOFS = ("ux", "uy", "uv", "rx", "ry", "rv")


class Builder:
    """Collects a model in abstract axes and renders the input document for given section sizes."""

    def __init__(self, unit_system: str, vertical: str, name: str, description: str = ""):
        if unit_system not in UNIT_SYSTEMS:
            raise ValueError(f"unknown unit system {unit_system!r}; use one of {sorted(UNIT_SYSTEMS)}")
        if vertical not in ("Y", "Z"):
            raise ValueError("vertical axis must be 'Y' or 'Z'")
        self.unit_system, self.vertical = unit_system, vertical
        self.us = UNIT_SYSTEMS[unit_system]
        self.name, self.description = name, description
        self.nodes: dict[str, tuple] = {}
        self.supports: dict[str, dict] = {}
        self.members: dict[str, dict] = {}
        self.member_group: dict[str, str] = {}
        self.constraints: dict[str, dict] = {}
        self.cases: dict[str, dict] = {}
        self.combinations: dict[str, dict] = {}
        self.analysis: dict = {"shearDeformation": False}
        self.groups: dict[str, dict] = {}           # sizing group -> {"kind", "index", "fixed"}
        self.design_defaults: dict = {"method": "LRFD"}
        self.design_groups: dict[str, dict] = {}
        self.design_members: dict[str, dict] = {}
        self.auto_members = True
        #: combination whose elastic stability the generator checks (None: nothing to buckle)
        self.stability_combination: str | None = "ULS1"
        self._late = []
        self.parameters: dict = {}

    # ---- quantities in file units, rounded to engineering values -----------------------------
    def L(self, metres: float, step: float = 0.25) -> float:
        """Length given in metres -> file units, on a grid of ``step`` m (metric) or 0.5 ft (imperial;
        1 ft for step >= 0.5)."""
        if self.us["metric"]:
            return tidy(_grid(metres, step) * self.us["per_m"])
        ft = 1.0 if step >= 0.5 else 0.5
        return tidy(_grid(metres / 0.3048, ft) * 12)

    def small(self, mm: float) -> float:
        """Small length given in mm (settlement, eccentricity) -> file units: 1 mm / 0.05 in grid."""
        if self.us["metric"]:
            return tidy(_grid(mm, 1.0) * self.us["per_m"] / 1000.0)
        return _grid(mm / 25.4, 0.05)

    def F(self, kn: float) -> float:
        """Force given in kN -> file units: 0.5 kN / 0.1 kip grid."""
        return _grid(kn, 0.5) if self.us["metric"] else _grid(kn / KIP, 0.1)

    def W(self, kn_per_m: float) -> float:
        """Line load given in kN/m -> file units: 0.5 kN/m / 0.005 kip/in grid."""
        if self.us["metric"]:
            return tidy(_grid(kn_per_m, 0.5) / self.us["per_m"])
        return _grid(kn_per_m / KIP * 0.0254, 0.005)

    def M(self, knm: float) -> float:
        """Moment given in kN m -> file units: 1 kN m / 10 kip in grid."""
        if self.us["metric"]:
            return tidy(_grid(knm, 1.0) * self.us["per_m"])
        return _grid(knm / KIP / 0.0254, 10.0)

    def k_rot(self, knm_per_rad: float) -> float:
        """Rotational stiffness given in kN m/rad -> file units, two significant digits."""
        if self.us["metric"]:
            return tidy(_sig(knm_per_rad) * self.us["per_m"])
        return _sig(knm_per_rad / KIP / 0.0254)

    def k_lin(self, kn_per_mm: float) -> float:
        """Translational stiffness given in kN/mm -> file units, two significant digits."""
        if self.us["metric"]:
            return tidy(_sig(kn_per_mm) * 1000.0 / self.us["per_m"])
        return _sig(kn_per_mm / KIP * 25.4)

    # ---- axes ------------------------------------------------------------------------------
    def vec(self, x: float, y: float, v: float) -> list:
        return tidy([x, y, v] if self.vertical == "Z" else [x, v, -y])

    def _components(self, values: dict, table: dict | None = None) -> dict:
        """Abstract DOF / load components -> file components, in file DOF order."""
        out = {}
        for key, value in values.items():
            name, sign = _DOF[self.vertical][table[key] if table else key]
            out[name] = value if not isinstance(value, (int, float)) else tidy(sign * value)
        return {k: out[k] for k in _ORDER if k in out}

    # ---- geometry --------------------------------------------------------------------------
    def node(self, nid: str, x: float, y: float, v: float) -> str:
        if nid in self.nodes:
            raise KeyError(f"duplicate node {nid}")
        self.nodes[nid] = (tidy(x), tidy(y), tidy(v))
        return nid

    def support(self, nid: str, kind: str | None = None, **dofs) -> None:
        """``kind`` "fixed" | "pinned"; ``dofs``: abstract DOF -> "fixed" | "free" | {"spring": k}."""
        s = self.supports.setdefault(nid, {})
        if kind == "fixed":
            s.update({d: "fixed" for d in ALL_DOFS})
        elif kind == "pinned":
            s.update({d: "fixed" for d in TRANSLATIONS})
        elif kind is not None:
            raise ValueError(kind)
        s.update(dofs)

    def restrain_out_of_plane(self, plane: str, nodes=None) -> None:
        """Hold every node in the plane "xv" or "yv" (out-of-plane translation and the two
        out-of-plane rotations fixed): a plane frame in 3D space."""
        dofs = ("uy", "rx", "rv") if plane == "xv" else ("ux", "ry", "rv")
        for nid in (nodes if nodes is not None else self.nodes):
            s = self.supports.setdefault(nid, {})
            for d in dofs:
                s.setdefault(d, "fixed")

    def group(self, name: str, kind: str, index: int, fixed: bool = False) -> str:
        """A sizing group: all its members share one catalogue section, chosen by the sizing loop."""
        n = catalogue_size(kind, self.unit_system)
        self.groups[name] = {"kind": kind, "index": max(0, min(index, n - 1)), "fixed": fixed}
        return name

    def member(self, mid: str, i: str, j: str, group: str, type: str = "frame", releases: dict | None = None,
               roll: float = 0.0, offsets: dict | None = None) -> str:
        """``releases``: {"i": {"rz": True | {"spring": k}}, ...} in member local DOFs."""
        if mid in self.members:
            raise KeyError(f"duplicate member {mid}")
        m = {"i": i, "j": j, "section": "@" + group, "material": "steel"}
        if type != "frame":
            m = {"type": type, **m}
        if roll:
            m["orientation"] = {"roll": tidy(roll)}
        if releases:
            m["releases"] = releases
        if offsets:
            m["offsets"] = offsets
        self.members[mid] = m
        self.member_group[mid] = group
        return mid

    def offset_half_depth(self, mid: str, end: str, group: str, direction) -> None:
        """Rigid end offset of member ``mid`` at ``end`` ("i" | "j") = half the depth of the section
        of sizing group ``group`` along the abstract unit vector ``direction`` (resolved when the
        sections are known)."""
        def fix(doc, depth):
            o = doc["members"][mid].setdefault("offsets", {"axes": "global"})
            o[end] = self.vec(*(c * depth[group] / 2.0 for c in direction))
        self._late.append(fix)

    def constraint(self, cid: str, type: str, master: str, slaves) -> None:
        self.constraints[cid] = {"type": type, "master": master, "slaves": list(slaves)}

    # ---- loads -----------------------------------------------------------------------------
    def case(self, cid: str, description: str = "") -> str:
        c = self.cases.setdefault(cid, {"loads": []})
        if description and "description" not in c:
            self.cases[cid] = {"description": description, **c}
        return cid

    def _add(self, case: str, load: dict) -> None:
        self.cases[self.case(case)]["loads"].append(load)

    def nodal(self, case: str, node: str, **components) -> None:
        """Components Fx, Fy, Fv, Mx, My, Mv in abstract axes."""
        comp = self._components({k: v for k, v in components.items() if v}, _LOAD)
        if comp:
            self._add(case, {"type": "nodal", "node": node, **{_FILE_LOAD[k]: v for k, v in comp.items()}})

    def settlement(self, case: str, node: str, **components) -> None:
        """Imposed support displacement, abstract DOF names (ux, uy, uv, rx, ry, rv)."""
        self._add(case, {"type": "supportDisplacement", "node": node, **self._components(components)})

    def _direction(self, direction: str):
        if direction in ("lx", "ly", "lz"):
            return direction[1], 1
        return _DIR[self.vertical][direction]

    def line(self, case: str, member: str, direction: str, w1: float, w2: float | None = None,
             start: float | None = None, end: float | None = None, projected: bool = False) -> None:
        """Line load; ``direction`` "x" | "y" | "v" (abstract global) or "lx" | "ly" | "lz" (local)."""
        d, sign = self._direction(direction)
        load = {"type": "memberDistributed", "member": member, "direction": d, "w1": tidy(sign * w1)}
        if w2 is not None and w2 != w1:
            load["w2"] = tidy(sign * w2)
        if start is not None:
            load["start"] = tidy(start)
        if end is not None:
            load["end"] = tidy(end)
        if projected:
            load["projected"] = True
        self._add(case, load)

    def point(self, case: str, member: str, direction: str, value: float, position: float,
              kind: str = "force") -> None:
        d, sign = self._direction(direction)
        load = {"type": "memberPoint", "member": member, "direction": d}
        if kind != "force":
            load["kind"] = kind
        load.update({"value": tidy(sign * value), "position": tidy(position)})
        self._add(case, load)

    def self_weight(self, case: str) -> None:
        self._add(case, {"type": "selfWeight"})

    def combination(self, cid: str, factors: dict, purpose: str = "strength", second_order: bool = False,
                    description: str = "") -> None:
        c = {}
        if description:
            c["description"] = description
        c["factors"] = {k: tidy(v) for k, v in factors.items() if v}
        if second_order:
            c["analysisType"] = "secondOrder"
        if purpose != "strength":
            c["purpose"] = purpose
        self.combinations[cid] = c

    # ---- design ----------------------------------------------------------------------------
    def design_group(self, gid: str, members, **variables) -> None:
        self.design_groups[gid] = {"members": list(members), **variables}

    def design_member(self, dmid: str, members, **variables) -> None:
        self.design_members[dmid] = {"members": list(members), **variables}

    # ---- output ----------------------------------------------------------------------------
    @property
    def nonlinear(self) -> bool:
        """True when the complete document differs from the fast one."""
        return "buckling" in self.analysis or any("analysisType" in c for c in self.combinations.values())

    def section_names(self) -> dict:
        return {g: section_definition(v["kind"], v["index"], self.unit_system)[0] for g, v in self.groups.items()}

    def request_buckling(self, wanted: bool, modes: int = 1) -> bool:
        """Ask for a linear buckling analysis of the stability combination in the file."""
        if wanted and self.stability_combination:
            self.analysis["buckling"] = {"combinations": [self.stability_combination], "modes": modes}
        return "buckling" in self.analysis

    def render(self, fast: bool = False) -> dict:
        """The input document for the current catalogue indices of the sizing groups.
        ``fast=True``: every combination first order and no buckling analysis -- the model the
        sizing iterations use before the complete document is analysed."""
        names, sections, depth = {}, {}, {}
        for g, v in self.groups.items():
            name, definition, d = section_definition(v["kind"], v["index"], self.unit_system)
            names[g], depth[g] = name, d
            sections.setdefault(name, definition)
        doc = {"formatVersion": "1.0.0", "name": self.name}
        if self.description:
            doc["description"] = self.description
        doc["units"] = dict(self.us["units"])
        doc["conventions"] = {"verticalAxis": self.vertical}
        doc["nodes"] = {}
        for nid, (x, y, v) in self.nodes.items():
            c = self.vec(x, y, v)
            doc["nodes"][nid] = {"x": c[0], "y": c[1], "z": c[2]}
        doc["supports"] = {nid: self._support(s) for nid, s in self.supports.items() if s}
        doc["materials"] = {"steel": {"description": "structural steel (nominal properties)", **self.us["material"]}}
        doc["sections"] = sections
        doc["members"] = {}
        for mid, m in self.members.items():
            m = {k: (dict(v) if isinstance(v, dict) else v) for k, v in m.items()}
            m["section"] = names[m["section"][1:]]
            doc["members"][mid] = m
        if self.constraints:
            doc["constraints"] = {k: dict(v, slaves=list(v["slaves"])) for k, v in self.constraints.items()}
        doc["loadCases"] = {cid: {**{k: v for k, v in c.items() if k != "loads"},
                                  "loads": [dict(x) for x in c["loads"]]} for cid, c in self.cases.items()}
        doc["combinations"] = {cid: {**c, "factors": dict(c["factors"])} for cid, c in self.combinations.items()}
        doc["analysis"] = _copy(self.analysis)
        if fast:
            doc["analysis"].pop("buckling", None)
            for c in doc["combinations"].values():
                c.pop("analysisType", None)
        design = {}
        if not self.auto_members:
            design["autoMembers"] = False
        design["defaults"] = _copy(self.design_defaults)
        if self.design_groups:
            design["groups"] = _copy(self.design_groups)
        if self.design_members:
            design["members"] = _copy(self.design_members)
        doc["design"] = design
        for fix in self._late:
            fix(doc, depth)
        return doc

    def _support(self, s: dict) -> dict:
        fixed = {d for d, v in s.items() if v == "fixed"}
        if fixed == set(ALL_DOFS):
            return {"type": "fixed"}
        rest = {d: v for d, v in s.items() if v != "free"}
        if set(TRANSLATIONS) <= fixed:
            extra = self._components({d: v for d, v in rest.items() if d not in TRANSLATIONS})
            return {"type": "pinned", **({"dofs": extra} if extra else {})}
        return {"dofs": self._components(rest)}


def _copy(x):
    if isinstance(x, dict):
        return {k: _copy(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_copy(v) for v in x]
    return x


# --------------------------------------------------------------------------------- family helpers
def flag(features: dict, name: str, drawn: bool) -> bool:
    """A feature flag: the forced value if the caller gave one, else the drawn one (the draw is
    always consumed first by the caller, so forcing a flag does not shift the random stream)."""
    return bool(features[name]) if name in features else bool(drawn)


def add_combinations(builder: "Builder", candidates, count: int, second_order: bool = False) -> list:
    """Add the first ``count`` of ``candidates`` [(id, factors, purpose)]; strength combinations are
    second order when asked.  Returns the ids added."""
    ids = []
    for cid, factors, purpose in list(candidates)[:count]:
        builder.combination(cid, factors, purpose, second_order=second_order and purpose == "strength")
        ids.append(cid)
    return ids


def service_ids(candidates, count: int) -> list:
    return [cid for cid, _, purpose in list(candidates)[:count] if purpose == "serviceability"]


#: Probability that a model of a level asks for a linear buckling analysis in the file: zero.  No
#: foreign solver reproduces linear buckling, so a file that asks for it cannot pass the cross-check
#: (verification tier C); the feature is therefore OPT-IN (``features={"buckling": True}``).  The draw
#: is still consumed, so that the flag does not change the rest of the model.
BUCKLING_CHANCE = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}
BUCKLING_FLAG = ("linear buckling analysis of combination ULS1 requested in the file (opt-in: not reproduced "
                 "by the foreign solvers)")
#: Same for tension-only bracing: opt-in (``features={"tensionOnly": True}``).
TENSION_ONLY_CHANCE = 0.0

BRACING_PATTERNS = (
    ("both flanges braced at the quarter points", {"topFlange": [0.25, 0.5, 0.75], "bottomFlange": [0.25, 0.5, 0.75]}),
    ("both flanges braced at the third points", {"topFlange": [0.3333, 0.6667], "bottomFlange": [0.3333, 0.6667]}),
    ("both flanges braced at midspan", {"topFlange": [0.5], "bottomFlange": [0.5]}),
    ("top flange braced at the fifth points, bottom flange at midspan",
     {"topFlange": [0.2, 0.4, 0.6, 0.8], "bottomFlange": [0.5]}),
)


def grid_lines(builder: "Builder", draw: "Draw", n: int, lo: float, hi: float, step: float = 0.5,
               equal: bool = True) -> list:
    """Coordinates of n + 1 grid lines with bay widths drawn in [lo, hi] m (file units)."""
    first = draw.step(lo, hi, step)
    widths = [first] * n
    if not equal:
        widths = [draw.step(lo, hi, step) for _ in range(n)]
    out, x = [0.0], 0.0
    for w in widths:
        x += builder.L(w, step)
        out.append(tidy(x))
    return out


def storey_levels(builder: "Builder", draw: "Draw", n: int, first=(4.0, 5.5), other=(3.0, 4.0)) -> list:
    """Elevations of the base and n floors (file units): a taller first storey, equal upper storeys."""
    h1, h = draw.step(first[0], first[1], 0.5), draw.step(other[0], other[1], 0.5)
    out, v = [0.0], 0.0
    for k in range(n):
        v += builder.L(h1 if k == 0 else h, 0.5)
        out.append(tidy(v))
    return out


def line_id(i: int) -> str:
    """Grid line letter: A, B, C, ..."""
    return chr(ord("A") + i)
