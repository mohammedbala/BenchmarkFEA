"""The resolved internal model: plain dataclasses in INTERNAL units (N, mm, rad).

``build_model(source)`` = load -> validate (all problems; refuse on error) ->
convert every dimensional number to internal units (schema-driven, exact
factors) -> resolve defaults, geometry, sections and design variables.

Separation of analysis and design data
--------------------------------------
The solver reads ``nodes, supports, materials, sections, members, constraints,
load_cases, combinations, analysis``.  Design-only variables live exclusively
in ``design_members`` and must never be read by ``benchfea.solve``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from . import defaults as D
from . import designparams, sections, shapes
from .designparams import Resolved
from .geometry import AXIS, is_vertical, local_triad, member_geometry
from .io import Document, load
from .problems import ERROR, ModelInvalid, Problem, pointer
from .sections import SectionProps
from .unitmap import convert_document
from .units import INTERNAL, UnitSystem
from .validate import DIAPHRAGM_DOFS, validate

DOFS = D.DOFS                         # ("ux", "uy", "uz", "rx", "ry", "rz") -- order of the 6 DOFs of a node
DOF_INDEX = {name: k for k, name in enumerate(DOFS)}
RIGID = math.inf                      # end-connection stiffness of an unreleased DOF

__all__ = ["Model", "Node", "Support", "Material", "Section", "Member", "Constraint", "NodalLoad",
           "MemberPointLoad", "MemberDistributedLoad", "SelfWeight", "SupportDisplacement", "LoadCase",
           "Combination", "NotionalLoad", "AnalysisOptions", "DesignParams", "DeflectionLimit", "DefaultUsed",
           "Resolved", "SectionProps", "build_model", "local_triad", "DOFS", "DOF_INDEX", "RIGID"]


@dataclass(frozen=True)
class DefaultUsed:
    """One default the engine applied because the file was silent.
    ``path`` = JSON pointer where the value could be set; ``value`` is in FILE units."""
    path: str
    value: Any
    reason: str

    def to_dict(self) -> dict:
        v = self.value
        return {"path": self.path, "value": list(v) if isinstance(v, tuple) else v, "reason": self.reason}


@dataclass(frozen=True)
class Node:
    id: str
    index: int                        # position in Model.nodes; global DOF k of this node = 6*index + k
    xyz: np.ndarray                   # (3,) mm


@dataclass(frozen=True)
class Support:
    node: str
    fixed: tuple[bool, ...]           # 6 flags, DOFS order: displacement prescribed (0 unless a load imposes one)
    spring: tuple[float, ...]         # 6 spring stiffnesses to ground (N/mm, N*mm/rad); 0 = none


@dataclass(frozen=True)
class Material:
    id: str
    E: float                          # N/mm^2
    G: float                          # N/mm^2
    nu: float | None
    unit_weight: float                # N/mm^3
    Fy: float | None                  # N/mm^2
    Fu: float | None
    description: str = ""


@dataclass(frozen=True)
class Section:
    """Section in internal units. Local y = depth direction, local z = width
    direction; Iz is the AISC major-axis (x-x) inertia of an I-shape."""
    id: str
    kind: str                         # "standard" | "properties" | "I" | "box" | "pipe" | "rectangle"
    props: SectionProps
    shape: str | None = None          # database designation for kind == "standard"

    def __getattr__(self, name):      # A, Iy, Iz, J, Asy, Asz, Cw, Sy, Sz, Zy, Zz, ry, rz, family, dims
        if name == "props":
            raise AttributeError(name)
        return getattr(self.props, name)


@dataclass(frozen=True)
class Member:
    id: str
    type: str                         # "frame" | "truss" | "tensionOnly" | "compressionOnly"
    i: str
    j: str
    section: str
    material: str
    start: np.ndarray                 # (3,) flexible start = node i + offset_i, mm
    end: np.ndarray                   # (3,) flexible end   = node j + offset_j
    length: float                     # flexible length, mm
    R: np.ndarray                     # (3,3) rows = local x, y, z in global components (v_local = R @ v_global)
    offset_i: np.ndarray              # (3,) rigid offset node i -> flexible start, GLOBAL components
    offset_j: np.ndarray
    k_end_i: tuple[float, ...]        # 6 connection stiffnesses at end i, local DOFS order:
    k_end_j: tuple[float, ...]        #   RIGID (inf) = not released, 0.0 = fully released, k = partial spring

    @property
    def has_releases(self) -> bool:
        return any(k != RIGID for k in self.k_end_i + self.k_end_j)

    @property
    def has_offsets(self) -> bool:
        return bool(self.offset_i.any() or self.offset_j.any())


@dataclass(frozen=True)
class Constraint:
    """u_slave = u_master + theta_master x (x_slave - x_master), theta_slave = theta_master,
    restricted to ``dofs`` (all six for a rigid link; the two in-plane translations and the
    rotation about the normal for a diaphragm, with only the in-plane part of the offset)."""
    id: str
    type: str                         # "rigidLink" | "rigidDiaphragm"
    master: str
    slaves: tuple[str, ...]
    normal: str | None                # "X" | "Y" | "Z" for a diaphragm
    dofs: tuple[str, ...]             # constrained DOF names of each slave


@dataclass(frozen=True)
class NodalLoad:
    node: str
    values: np.ndarray                # (6,) Fx, Fy, Fz [N], Mx, My, Mz [N*mm], global


@dataclass(frozen=True)
class MemberPointLoad:
    member: str
    kind: str                         # "force" | "moment"
    vector: np.ndarray                # (3,) force [N] or moment [N*mm] in member LOCAL axes
    x: float                          # distance from the flexible start, mm


@dataclass(frozen=True)
class MemberDistributedLoad:
    member: str
    w1: np.ndarray                    # (3,) intensity at x1 in LOCAL axes, N per mm of member length
    w2: np.ndarray                    # (3,) intensity at x2 (linear in between); projection already applied
    x1: float                         # mm from the flexible start
    x2: float


@dataclass(frozen=True)
class SelfWeight:
    factor: float
    direction: np.ndarray             # (3,) unit vector, global
    members: tuple[str, ...]          # resolved member ids
    # line load on each member = factor * material.unit_weight * section.A along `direction`


@dataclass(frozen=True)
class SupportDisplacement:
    node: str
    values: dict[int, float]          # DOF index (0..5) -> imposed displacement (mm or rad), global


@dataclass(frozen=True)
class LoadCase:
    id: str
    loads: tuple
    description: str = ""


@dataclass(frozen=True)
class NotionalLoad:
    coefficient: float
    direction: np.ndarray             # (3,) unit vector, global, horizontal


@dataclass(frozen=True)
class Combination:
    id: str
    factors: dict[str, float]
    analysis_type: str                # "firstOrder" | "secondOrder"
    purpose: str                      # "strength" | "serviceability"
    stiffness_reduction: bool
    notional: NotionalLoad | None
    description: str = ""


@dataclass(frozen=True)
class AnalysisOptions:
    type: str
    shear_deformation: bool
    result_stations: int
    convergence_tolerance: float
    max_iterations: int
    alpha: float
    equilibrium_tolerance: float
    pivot_tolerance: float
    geometry_tolerance: float
    buckling_combinations: tuple[str, ...]
    buckling_modes: int
    reduction_factor: float
    tau_b: Any                        # "computed" | "none" | float
    reduction_members: tuple[str, ...] | None     # None = all members


@dataclass(frozen=True)
class DeflectionLimit:
    span_ratio: float | None          # limit = design member length / span_ratio
    absolute: float | None            # mm
    direction: str                    # "y" | "z" | "resultant"
    combinations: tuple[str, ...]


@dataclass(frozen=True)
class DesignParams:
    """Design variables of one design member, each with its provenance.

    ``values[name]`` is a :class:`Resolved` (value, source, path, reason); names are
    ``benchfea.designparams.NAMES``.  ``p.Kx`` is shorthand for ``p.values["Kx"].value``.
    Lengths in mm, areas in mm^2, E/G/Fy/Fu in N/mm^2; brace locations are fractions of
    ``length``.  ``E, G, Fy, Fu`` default to the member's material (source ``"material:<id>"``).
    AISC axis names: x = major = local z, y = minor = local y, z = torsion.
    """
    id: str
    members: tuple[str, ...]          # analysis member ids, start to end
    length: float                     # mm, flexible start of the first to flexible end of the last member
    section: str
    material: str
    values: dict[str, Resolved]

    def __getattr__(self, name):
        if name == "values":
            raise AttributeError(name)
        try:
            return self.values[name].value
        except KeyError:
            raise AttributeError(name) from None

    def source(self, name: str) -> str:
        return self.values[name].source


@dataclass
class Model:
    name: str
    format_version: str
    units: UnitSystem                 # units of the FILE (results are reported back in these)
    vertical_axis: str                # "Y" | "Z"
    nodes: dict[str, Node]
    supports: dict[str, Support]
    materials: dict[str, Material]
    sections: dict[str, Section]
    members: dict[str, Member]
    constraints: dict[str, Constraint]
    load_cases: dict[str, LoadCase]
    combinations: dict[str, Combination]
    analysis: AnalysisOptions
    design_members: dict[str, DesignParams]        # DESIGN ONLY -- never read by the solver
    defaults_used: list[DefaultUsed] = field(default_factory=list)
    warnings: list[Problem] = field(default_factory=list)
    source: str | None = None

    @property
    def n_dof(self) -> int:
        return 6 * len(self.nodes)

    def dof(self, node_id: str, dof: "str | int") -> int:
        """Global DOF number: 6 * node index + (ux, uy, uz, rx, ry, rz)."""
        k = DOF_INDEX[dof] if isinstance(dof, str) else dof
        return 6 * self.nodes[node_id].index + k

    @property
    def vertical(self) -> np.ndarray:
        return AXIS[self.vertical_axis].copy()


# ============================================================================== builder
def build_model(source) -> Model:
    """Load, validate and resolve a model. Raises :class:`ModelInvalid` carrying
    every problem if the document has any error."""
    doc: Document = load(source)
    problems = validate(doc)
    if any(p.severity == ERROR for p in problems):
        raise ModelInvalid(problems, doc.source)
    file_units = UnitSystem.from_dict(doc.data["units"])
    raw = doc.data
    d = convert_document(raw, INTERNAL)
    used: list[DefaultUsed] = []

    def default(path, value, reason, dim=None):
        if dim and value is not None:
            value = file_units.from_internal(value, dim)
        if isinstance(value, tuple):
            value = list(value)             # plain JSON data in the report
        used.append(DefaultUsed(path, value, reason))

    vertical = d["conventions"]["verticalAxis"]
    for key, val in D.CONVENTIONS.items():
        if key not in d["conventions"]:
            default(f"/conventions/{key}", val, "only supported convention")

    # ---- analysis options
    an = d.get("analysis") or {}
    for key, val in D.ANALYSIS.items():
        if key not in an:
            default(f"/analysis/{key}", val, D.ANALYSIS_REASONS[key])
    opt = {**D.ANALYSIS, **{k: v for k, v in an.items() if k in D.ANALYSIS}}
    buck = an.get("buckling") or {}
    if buck and "modes" not in buck:
        default("/analysis/buckling/modes", D.BUCKLING["modes"], "number of buckling modes")
    sr = an.get("stiffnessReduction") or {}
    combos_raw = d["combinations"]
    if any(c.get("stiffnessReduction") for c in combos_raw.values()):
        for key, val in D.STIFFNESS_REDUCTION.items():
            if key not in sr:
                default(f"/analysis/stiffnessReduction/{key}", val, D.STIFFNESS_REDUCTION_REASONS[key])
        if "members" not in sr:
            default("/analysis/stiffnessReduction/members", list(d["members"]),
                    "stiffness reduction applied to all members (AISC 360-22 C2.3(a) permits it for all "
                    "stiffnesses, p. 16.1-30)")
    analysis = AnalysisOptions(
        opt["type"], opt["shearDeformation"], opt["resultStations"], opt["convergenceTolerance"],
        opt["maxIterations"], opt["alpha"], opt["equilibriumTolerance"], opt["pivotTolerance"],
        opt["geometryTolerance"], tuple(buck.get("combinations", ())), buck.get("modes", D.BUCKLING["modes"]),
        sr.get("factor", D.STIFFNESS_REDUCTION["factor"]), sr.get("tauB", D.STIFFNESS_REDUCTION["tauB"]),
        tuple(sr["members"]) if "members" in sr else None)

    # ---- nodes, supports
    nodes = {nid: Node(nid, k, np.array([n["x"], n["y"], n.get("z", 0.0)], float))
             for k, (nid, n) in enumerate(d["nodes"].items())}
    for nid, n in d["nodes"].items():
        if "z" not in n:
            default(pointer("nodes", nid, "z"), 0.0, "no z coordinate given: the node lies in the plane z = 0")
    vdof = {"X": "ux", "Y": "uy", "Z": "uz"}[vertical]
    supports = {}
    for nid, s in (d.get("supports") or {}).items():
        base = {"fixed": DOFS, "pinned": DOFS[:3], "roller": (vdof,)}.get(s.get("type"), ())
        spec = {dof: ("fixed" if dof in base else "free") for dof in DOFS}
        spec.update(s.get("dofs") or {})
        supports[nid] = Support(nid, tuple(spec[k] == "fixed" for k in DOFS),
                                tuple(float(spec[k]["spring"]) if isinstance(spec[k], dict) else 0.0 for k in DOFS))

    # ---- materials
    any_self_weight = any(ld["type"] == "selfWeight" for c in d["loadCases"].values() for ld in c.get("loads", ()))
    materials = {}
    for mid, m in d["materials"].items():
        if "G" in m:
            G = m["G"]
        else:
            G = m["E"] / (2.0 * (1.0 + m["nu"]))
            default(pointer("materials", mid, "G"), G, "G = E / (2 (1 + nu))", "force/length^2")
        if "unitWeight" not in m and any_self_weight:
            default(pointer("materials", mid, "unitWeight"), 0.0, "no unit weight given: self weight is zero")
        materials[mid] = Material(mid, float(m["E"]), float(G), m.get("nu"), float(m.get("unitWeight", 0.0)),
                                  m.get("Fy"), m.get("Fu"), m.get("description", ""))

    # ---- sections
    secs = {}
    for sid, s in d["sections"].items():
        base = pointer("sections", sid)
        props = sections.from_definition(s, lambda name: shapes.lookup(name, INTERNAL.length))
        if s["type"] == "properties":
            for key in ("Asy", "Asz"):
                if key not in s:
                    default(f"{base}/{key}", None, "no shear area given: no shear deformation in this direction")
            if "Cw" not in s:
                default(f"{base}/Cw", 0.0, "warping constant taken as zero")
        if s["type"] == "I" and "r" not in s:
            default(f"{base}/r", 0.0, "no fillets")
        secs[sid] = Section(sid, s["type"], props, s.get("shape"))

    # ---- members
    size = float(np.linalg.norm(np.ptp(np.array([n.xyz for n in nodes.values()]), axis=0)))
    tol = analysis.geometry_tolerance * size
    members, geoms = {}, {}
    for mid, m in d["members"].items():
        base = pointer("members", mid)
        if "type" not in m:
            default(f"{base}/type", "frame", "member type frame")
        mtype = m.get("type", "frame")
        g = member_geometry(m, d["nodes"], vertical, tol)
        geoms[mid] = g
        start, end, R, length, o_i, o_j = g
        if "referenceVector" not in (m.get("orientation") or {}):
            if is_vertical(R[0], vertical):
                default(f"{base}/orientation/referenceVector", [1.0, 0.0, 0.0],
                        "default orientation of a vertical member: local y = global +X (before roll)")
            else:
                default(f"{base}/orientation/referenceVector", [float(x) for x in AXIS[vertical]],
                        "default orientation: local y in the plane of the member and the vertical axis, "
                        "on the upward side (before roll)")
        if "offsets" in m and "axes" not in m["offsets"]:
            default(f"{base}/offsets/axes", "global", "rigid end offsets are global vectors")
        rel = m.get("releases") or {}

        def stiff(end_name):
            if mtype != "frame":
                # truss-type: only the axial force is transmitted (rotations and shears free)
                return (RIGID, RIGID, RIGID, 0.0, 0.0, 0.0)
            out = []
            for dof in DOFS:
                v = (rel.get(end_name) or {}).get(dof, False)
                out.append(float(v["spring"]) if isinstance(v, dict) else (0.0 if v else RIGID))
            return tuple(out)

        members[mid] = Member(mid, mtype, m["i"], m["j"], m["section"], m["material"], start, end, length, R,
                              o_i, o_j, stiff("i"), stiff("j"))

    # ---- constraints
    constraints = {}
    for cid, c in (d.get("constraints") or {}).items():
        if c["type"] == "rigidDiaphragm":
            if "normal" not in c:
                default(pointer("constraints", cid, "normal"), vertical, "diaphragm normal = vertical axis")
            normal = c.get("normal", vertical)
            dofs = DIAPHRAGM_DOFS[normal]
        else:
            normal, dofs = None, DOFS
        constraints[cid] = Constraint(cid, c["type"], c["master"], tuple(c["slaves"]), normal, tuple(dofs))

    # ---- load cases
    cases = {}
    for cid, case in d["loadCases"].items():
        loads = []
        for k, ld in enumerate(case.get("loads", ())):
            base = f"{pointer('loadCases', cid)}/loads/{k}"
            t = ld["type"]
            if t == "nodal":
                loads.append(NodalLoad(ld["node"], np.array([ld.get(c, 0.0) for c in
                                                             ("Fx", "Fy", "Fz", "Mx", "My", "Mz")], float)))
            elif t == "supportDisplacement":
                loads.append(SupportDisplacement(ld["node"], {DOF_INDEX[k2]: float(ld[k2]) for k2 in DOFS if k2 in ld}))
            elif t == "selfWeight":
                if "direction" in ld:
                    g = np.asarray(ld["direction"], float)
                    g = g / np.linalg.norm(g)
                else:
                    g = -AXIS[vertical]
                    default(f"{base}/direction", [float(x) for x in g], "gravity acts along the negative vertical axis")
                if "factor" not in ld:
                    default(f"{base}/factor", 1.0, "self weight = 1.0 x unit weight x area")
                if "members" not in ld:
                    default(f"{base}/members", list(members), "self weight acts on all members")
                loads.append(SelfWeight(float(ld.get("factor", 1.0)), g, tuple(ld.get("members", members))))
            else:
                mem = members[ld["member"]]
                if "positionType" not in ld:
                    default(f"{base}/positionType", "relative", "positions are fractions of the flexible length")
                scale = 1.0 if ld.get("positionType", "relative") == "absolute" else mem.length
                direction = ld["direction"]
                if direction in AXIS:
                    unit_global = AXIS[direction]
                    unit_local = mem.R @ unit_global
                else:
                    unit_local = np.eye(3)["xyz".index(direction)]
                    unit_global = None
                if t == "memberPoint":
                    if "kind" not in ld:
                        default(f"{base}/kind", "force", "a concentrated force (not a moment)")
                    loads.append(MemberPointLoad(mem.id, ld.get("kind", "force"), ld["value"] * unit_local,
                                                 ld["position"] * scale))
                else:
                    f = 1.0
                    if ld.get("projected", False):
                        # per unit projected length -> per unit member length
                        f = float(np.linalg.norm(np.cross(mem.R[0], unit_global)))
                    x2 = ld["end"] * scale if "end" in ld else mem.length
                    loads.append(MemberDistributedLoad(mem.id, ld["w1"] * f * unit_local,
                                                       ld.get("w2", ld["w1"]) * f * unit_local,
                                                       ld.get("start", 0.0) * scale, min(x2, mem.length)))
        cases[cid] = LoadCase(cid, tuple(loads), case.get("description", ""))

    # ---- combinations
    combos = {}
    for cid, c in combos_raw.items():
        base = pointer("combinations", cid)
        if "analysisType" not in c:
            default(f"{base}/analysisType", analysis.type, "analysis.type")
        for key, val in D.COMBINATION.items():
            if key not in c:
                default(f"{base}/{key}", val, {"purpose": "combinations are strength combinations unless stated",
                                               "stiffnessReduction": "no stiffness reduction unless switched on"}[key])
        nl = c.get("notionalLoad")
        notional = None
        if nl:
            if "coefficient" not in nl:
                default(f"{base}/notionalLoad/coefficient", D.NOTIONAL["coefficient"], D.NOTIONAL_REASON)
            dv = np.asarray(nl["direction"], float)
            notional = NotionalLoad(float(nl.get("coefficient", D.NOTIONAL["coefficient"])), dv / np.linalg.norm(dv))
        combos[cid] = Combination(cid, {k: float(v) for k, v in c["factors"].items()},
                                  c.get("analysisType", analysis.type), c.get("purpose", D.COMBINATION["purpose"]),
                                  bool(c.get("stiffnessReduction", D.COMBINATION["stiffnessReduction"])), notional,
                                  c.get("description", ""))

    # ---- design (design-only data)
    def no_problem(*args, **kw):          # all design problems were reported by validate()
        pass

    design_raw = d.get("design") or {}
    if "autoMembers" not in design_raw:
        default("/design/autoMembers", True, "every analysis member not listed in an explicit design member "
                                             "is a design member by itself")
    levels = [("/design/defaults", design_raw.get("defaults") or {})]
    levels += [(pointer("design", "groups", k), v) for k, v in (design_raw.get("groups") or {}).items()]
    levels += [(pointer("design", "members", k), v) for k, v in (design_raw.get("members") or {}).items()]
    for lbase, level in levels:
        if "bracing" in level and "positionType" not in level["bracing"]:
            default(f"{lbase}/bracing/positionType", "relative",
                    "brace locations are fractions of the design member length")
    service = tuple(cid for cid, c in combos.items() if c.purpose == "serviceability")
    design_members = {}
    for dmid, info in designparams.resolve(d, geoms, no_problem).items():
        first = members[info["members"][0]]
        values = dict(info["values"])
        if values["An"].value is None:
            r = values["An"]
            values["An"] = Resolved(secs[first.section].A, r.source, r.path, r.reason)
        lims = []
        for k, lim in enumerate(values["deflectionLimits"].value):
            r = values["deflectionLimits"]
            if r.source != "default":
                if "direction" not in lim:
                    default(f"{r.path}/{k}/direction", "resultant", "resultant of local y and z deflection")
                if "combinations" not in lim:
                    default(f"{r.path}/{k}/combinations", list(service), "all serviceability combinations")
            lims.append(DeflectionLimit(lim.get("spanRatio"), lim.get("absolute"), lim.get("direction", "resultant"),
                                        tuple(lim.get("combinations", service))))
        r = values["deflectionLimits"]
        values["deflectionLimits"] = Resolved(tuple(lims), r.source, r.path, r.reason)
        for name, r in values.items():
            if r.source == "default":
                default(r.path, r.value, r.reason, designparams.DIMENSIONS.get(name))
        design_members[dmid] = DesignParams(dmid, tuple(info["members"]), info["length"], first.section,
                                            first.material, values)

    return Model(
        name=d.get("name", ""), format_version=d["formatVersion"], units=file_units, vertical_axis=vertical,
        nodes=nodes, supports=supports, materials=materials, sections=secs, members=members,
        constraints=constraints, load_cases=cases, combinations=combos, analysis=analysis,
        design_members=design_members, defaults_used=used,
        warnings=[p for p in problems if p.severity != ERROR], source=doc.source)
