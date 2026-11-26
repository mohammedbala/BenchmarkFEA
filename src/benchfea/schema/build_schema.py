"""Generator of the published JSON Schema (draft 2020-12).

Run ``python -m benchfea.schema.build_schema`` to rewrite
``benchfea-input-1.0.schema.json`` next to this file.  A test checks that the
committed JSON equals the generator output.

Annotation keywords used in addition to the standard vocabulary:

``x-unit``  dimension of the value: a product/quotient of ``force``, ``length``,
            ``moment`` (force*length), ``rotation`` (radian), ``degree``,
            ``dimensionless``; ``position`` (fraction of the member length or a
            length, depending on the sibling ``positionType``); or ``n/a`` for
            identifiers, enumerations, flags and containers.
``x-hint``  text shown as the hint when the value fails an anyOf/oneOf.
"""
from __future__ import annotations

import json
from pathlib import Path

from .. import defaults as D

SCHEMA_FILE = Path(__file__).with_name("benchfea-input-1.0.schema.json")
SCHEMA_ID = "urn:benchfea:schema:input:1.0"
ID_PATTERN = r"^[A-Za-z0-9_][A-Za-z0-9_.\-]*$"


def num(desc, unit, default=None, **kw):
    s = {"type": "number", "description": desc, "x-unit": unit}
    s.update(kw)
    if default is not None:
        s["default"] = default
    return s


def integer(desc, default=None, **kw):
    s = {"type": "integer", "description": desc, "x-unit": "dimensionless"}
    s.update(kw)
    if default is not None:
        s["default"] = default
    return s


def text(desc, **kw):
    return {"type": "string", "description": desc, "x-unit": "n/a", **kw}


def enum(desc, values, default=None):
    s = {"enum": list(values), "description": desc, "x-unit": "n/a"}
    if default is not None:
        s["default"] = default
    return s


def flag(desc, default=None):
    s = {"type": "boolean", "description": desc, "x-unit": "n/a"}
    if default is not None:
        s["default"] = default
    return s


def ref(name, desc):
    return {"$ref": f"#/$defs/{name}", "description": desc, "x-unit": "n/a"}


def idref(desc):
    return {"type": "string", "pattern": ID_PATTERN, "description": desc, "x-unit": "n/a"}


def idlist(desc, **kw):
    return {"type": "array", "description": desc, "x-unit": "n/a", "uniqueItems": True,
            "items": {"type": "string", "pattern": ID_PATTERN, "description": "An id.", "x-unit": "n/a"}, **kw}


def idmap(desc, item, **kw):
    return {"type": "object", "description": desc, "x-unit": "n/a",
            "propertyNames": {"pattern": ID_PATTERN, "maxLength": 64},
            "additionalProperties": item, **kw}


def obj(desc, props, required=(), **kw):
    s = {"type": "object", "description": desc, "x-unit": "n/a", "properties": props}
    if required:
        s["required"] = list(required)
    s.setdefault("additionalProperties", False)
    s.update(kw)
    return s


def vec3(desc, unit, **kw):
    return {"type": "array", "description": desc, "x-unit": unit, "minItems": 3, "maxItems": 3,
            "items": {"type": "number", "description": "Vector component (X, Y, Z order).", "x-unit": unit}, **kw}


def variant(type_name, props, required=()):
    then = {"properties": props}
    if required:
        then["required"] = list(required)
    return {"if": {"properties": {"type": {"const": type_name}}, "required": ["type"]}, "then": then}


def build() -> dict:
    A, SR = D.ANALYSIS, D.STIFFNESS_REDUCTION
    defs: dict = {}

    # ------------------------------------------------------------------ supports
    def support_dof(kind, unit):
        return {
            "description": f"Restraint of this {kind} DOF: \"fixed\", \"free\", or {{\"spring\": k}} "
                           "(linear elastic spring to ground). Absent: free, unless the shorthand `type` "
                           "fixes this DOF.",
            "x-unit": "n/a", "examples": ["fixed", {"spring": 5000}],
            "x-hint": "Use \"fixed\", \"free\" or {\"spring\": <stiffness > 0>}.",
            "anyOf": [
                {"enum": ["fixed", "free"]},
                obj("Elastic spring support.", {"spring": num("Spring stiffness to ground.", unit, exclusiveMinimum=0,
                                                        examples=[5000])},
                    required=["spring"]),
            ],
        }
    tsup = support_dof("translational", "force/length")
    rsup = support_dof("rotational", "moment/rotation")
    defs["support"] = obj(
        "Support at a node (the key is the node id). Give `type` (shorthand), `dofs` (per-DOF), or both; "
        "entries in `dofs` override the shorthand. DOFs not mentioned are free. Imposed support "
        "displacements are loads (`supportDisplacement`) acting on fixed DOFs.",
        {
            "type": enum("Shorthand: fixed = all six DOFs fixed; pinned = ux, uy, uz fixed; roller = only the "
                         "translation along the vertical axis fixed. Absent: only the entries of `dofs` "
                         "apply.", ["fixed", "pinned", "roller"]),
            "dofs": obj("Per-DOF restraint in global axes. Absent: only the shorthand `type` applies.", {
                "ux": tsup, "uy": tsup, "uz": tsup, "rx": rsup, "ry": rsup, "rz": rsup}),
        },
        anyOf=[{"required": ["type"]}, {"required": ["dofs"]}],
        **{"x-hint": "A support needs `type` and/or `dofs`."})

    # ----------------------------------------------------------------- materials
    defs["material"] = obj(
        "Linear elastic isotropic material. Give G or nu (or both; G is used for stiffness). The design "
        "checks use the same E, G, Fy, Fu unless a design variable of the same name overrides them.",
        {
            "description": text("Free text, e.g. the steel grade designation. Absent: no text."),
            "E": num("Modulus of elasticity.", "force/length^2", exclusiveMinimum=0),
            "G": num("Shear modulus. Default: E / (2 (1 + nu)).", "force/length^2", exclusiveMinimum=0),
            "nu": num("Poisson's ratio, used only to derive G when G is absent. Absent: G must be given.",
                      "dimensionless", minimum=0, exclusiveMaximum=0.5, examples=[0.3]),
            "unitWeight": num("Weight per unit volume, used by selfWeight loads.", "force/length^3", 0.0, minimum=0),
            "Fy": num("Specified minimum yield stress. Required for design checks and for computed tau_b. "
                      "Absent: no value.", "force/length^2", exclusiveMinimum=0),
            "Fu": num("Specified minimum tensile strength. Required for design checks. Absent: no value.",
                      "force/length^2", exclusiveMinimum=0),
        },
        required=["E"], anyOf=[{"required": ["G"]}, {"required": ["nu"]}],
        **{"x-hint": "A material needs `G` or `nu` in addition to `E`."})

    # ------------------------------------------------------------------ sections
    keep = "Absent: the computed or looked-up value is kept."
    overrides = obj(
        "Explicit values that replace the computed/looked-up section properties (e.g. a tabulated J or a "
        "shear area from a refined theory). Absent: every property is the computed or looked-up one.",
        {
            "A": num("Cross-sectional area. " + keep, "length^2", exclusiveMinimum=0, examples=[0.0125]),
            "Iy": num("Second moment of area about local y. " + keep, "length^4", exclusiveMinimum=0,
                      examples=[2.1e-5]),
            "Iz": num("Second moment of area about local z. " + keep, "length^4", exclusiveMinimum=0,
                      examples=[4.2e-4]),
            "J": num("Saint-Venant torsion constant. " + keep, "length^4", exclusiveMinimum=0,
                     examples=[6.5e-7]),
            "Asy": num("Shear area for shear force along local y. " + keep, "length^2", exclusiveMinimum=0,
                       examples=[0.004]),
            "Asz": num("Shear area for shear force along local z. " + keep, "length^2", exclusiveMinimum=0,
                       examples=[0.0065]),
            "Cw": num("Warping constant. " + keep, "length^6", minimum=0, examples=[1.9e-6]),
        }, examples=[{"J": 6.5e-7}])
    L = "length"
    pos = {"exclusiveMinimum": 0}
    defs["section"] = {
        "type": "object", "x-unit": "n/a",
        "description": "Prismatic cross-section. `type` selects how it is defined: by standard shape name "
                       "(user-supplied AISC Shapes Database), by properties, or by parametric dimensions from "
                       "which the engine computes the properties. Section axes: y = depth direction (web of an "
                       "I-shape), z = width direction (parallel to the flanges); AISC x-x (major) axis = "
                       "local z, AISC y-y (minor) axis = local y.",
        "required": ["type"],
        "properties": {
            "type": enum("Kind of definition.", ["standard", "properties", "I", "box", "pipe", "rectangle"]),
            "description": text("Free text. Absent: no text.", examples=["floor beam"]),
        },
        "allOf": [
            variant("standard", {
                "shape": text("Shape designation exactly as in the AISC Shapes Database v16.0 "
                              "(e.g. \"W14X90\"). Needs the user-supplied database file.", minLength=1,
                              examples=["W14X90"]),
                "overrides": overrides}, ["shape"]),
            variant("properties", {
                "A": num("Cross-sectional area.", "length^2", examples=[0.0125], **pos),
                "Iy": num("Second moment of area about local y (AISC minor axis for an I-shape).", "length^4",
                          examples=[2.1e-5], **pos),
                "Iz": num("Second moment of area about local z (AISC major axis for an I-shape).", "length^4",
                          examples=[4.2e-4], **pos),
                "J": num("Saint-Venant torsion constant.", "length^4", examples=[6.5e-7], **pos),
                "Asy": num("Shear area for shear along local y. Absent: no shear deformation in that direction.",
                           "length^2", examples=[0.004], **pos),
                "Asz": num("Shear area for shear along local z. Absent: no shear deformation in that direction.",
                           "length^2", examples=[0.0065], **pos),
                "Cw": num("Warping constant (design only; warping torsion is not analysed).", "length^6", 0.0,
                          minimum=0),
                "Sy": num("Elastic section modulus about local y (design only). Absent: no value.", "length^3",
                          examples=[2.1e-4], **pos),
                "Sz": num("Elastic section modulus about local z (design only). Absent: no value.", "length^3",
                          examples=[1.87e-3], **pos),
                "Zy": num("Plastic section modulus about local y (design only). Absent: no value.", "length^3",
                          examples=[3.3e-4], **pos),
                "Zz": num("Plastic section modulus about local z (design only). Absent: no value.", "length^3",
                          examples=[2.1e-3], **pos),
            }, ["A", "Iy", "Iz", "J"]),
            variant("I", {
                "d": num("Overall depth (along local y).", L, **pos),
                "bf": num("Flange width (along local z); both flanges equal.", L, **pos),
                "tf": num("Flange thickness.", L, **pos),
                "tw": num("Web thickness.", L, **pos),
                "r": num("Web-to-flange fillet radius.", L, 0.0, minimum=0),
                "overrides": overrides}, ["d", "bf", "tf", "tw"]),
            variant("box", {
                "H": num("Overall depth (along local y).", L, examples=[0.3], **pos),
                "B": num("Overall width (along local z).", L, examples=[0.2], **pos),
                "t": num("Wall thickness (uniform, sharp corners).", L, examples=[0.01], **pos),
                "overrides": overrides}, ["H", "B", "t"]),
            variant("pipe", {
                "D": num("Outside diameter.", L, **pos),
                "t": num("Wall thickness.", L, **pos),
                "overrides": overrides}, ["D", "t"]),
            variant("rectangle", {
                "h": num("Depth (along local y).", L, examples=[0.4], **pos),
                "b": num("Width (along local z).", L, examples=[0.25], **pos),
                "overrides": overrides}, ["h", "b"]),
        ],
        "unevaluatedProperties": False,
    }

    # ------------------------------------------------------------------- members
    def release(kind, unit):
        return {
            "description": f"Release of this {kind} DOF between the member end and the node (member local "
                           "axes): false = rigidly connected, true = fully released, {\"spring\": k} = "
                           "partial release through a linear spring.",
            "x-unit": "n/a", "default": False,
            "x-hint": "Use true, false or {\"spring\": <stiffness >= 0>}.",
            "anyOf": [
                {"type": "boolean"},
                obj("Partial release.", {"spring": num("Connection spring stiffness (0 = fully released).", unit,
                                                       minimum=0, examples=[20000])}, required=["spring"]),
            ],
        }
    trel, rrel = release("translational", "force/length"), release("rotational", "moment/rotation")
    defs["endRelease"] = obj("Releases at one member end. DOFs not mentioned are rigidly connected.",
                             {"ux": trel, "uy": trel, "uz": trel, "rx": rrel, "ry": rrel, "rz": rrel})
    defs["member"] = obj(
        "Straight prismatic member between two nodes.",
        {
            "type": enum("frame = axial, shear, bending, torsion; truss = axial force only (pin-ended); "
                         "tensionOnly / compressionOnly = truss that is removed when it would carry the other "
                         "sign (iterative analysis).", ["frame", "truss", "tensionOnly", "compressionOnly"], "frame"),
            "i": idref("Start node id (local x runs from i to j)."),
            "j": idref("End node id."),
            "section": idref("Section id."),
            "material": idref("Material id."),
            "orientation": obj(
                "Orientation of the local y and z axes about local x. Default rule: local y lies in the "
                "plane of local x and the global vertical axis and points to the +vertical side; for a "
                "vertical member local y = global +X. See docs/CONVENTIONS.md. Absent: the default rule, "
                "no roll.",
                {
                    "referenceVector": vec3("Global vector replacing the vertical axis in the default rule: "
                                            "local y lies in the plane of local x and this vector, on the "
                                            "vector's side. Must not be parallel to the member. Absent: the "
                                            "vertical axis is used.", "dimensionless", examples=[[0, 1, 0]]),
                    "roll": num("Rotation of local y and z about +local x (right-hand rule), applied after "
                                "the default rule or the reference vector.", "degree", 0.0),
                }, examples=[{"roll": 90}]),
            "releases": obj("End releases in member local axes. Absent: both ends rigidly connected.",
                            {"i": ref("endRelease", "Releases at end i. Absent: end i rigidly connected."),
                             "j": ref("endRelease", "Releases at end j. Absent: end j rigidly connected.")}),
            "offsets": obj(
                "Rigid end offsets: vectors from each node to the corresponding end of the flexible length. "
                "Absent: no offsets.",
                {
                    "axes": enum("Axes of the offset vectors: global, or local = the triad of the "
                                 "node-to-node line (same orientation rule).", ["global", "local"], "global"),
                    "i": vec3("Offset from node i to the flexible start.", L, default=[0, 0, 0]),
                    "j": vec3("Offset from node j to the flexible end.", L, default=[0, 0, 0]),
                }),
        },
        required=["i", "j", "section", "material"])

    # --------------------------------------------------------------- constraints
    defs["constraint"] = {
        "type": "object", "x-unit": "n/a",
        "description": "Exact kinematic constraint (small-displacement rigid-body relation between a master "
                       "node and slave nodes).",
        "required": ["type", "master", "slaves"],
        "properties": {
            "type": enum("rigidLink = all six DOFs of every slave follow the master as a rigid body; "
                         "rigidDiaphragm = the two in-plane translations and the rotation about the normal "
                         "follow the master, other DOFs stay independent.", ["rigidLink", "rigidDiaphragm"]),
            "master": idref("Master node id."),
            "slaves": idlist("Slave node ids.", minItems=1),
        },
        "allOf": [variant("rigidDiaphragm", {
            "normal": enum("Global axis normal to the diaphragm plane. Default: the vertical axis.",
                           ["X", "Y", "Z"])})],
        "unevaluatedProperties": False,
    }

    # --------------------------------------------------------------------- loads
    direction = enum("Load direction: X, Y, Z = global axes; x, y, z = member local axes. Positive values act "
                     "along the positive axis.", ["X", "Y", "Z", "x", "y", "z"])
    ptype = enum("How positions along the member are given: relative = fraction of the flexible length "
                 "(0 at end i, 1 at end j); absolute = distance from the flexible start.",
                 ["relative", "absolute"], "relative")
    defs["load"] = {
        "type": "object", "x-unit": "n/a",
        "description": "One load. `type` selects the kind.",
        "required": ["type"],
        "properties": {"type": enum("Kind of load.", ["nodal", "memberPoint", "memberDistributed", "selfWeight",
                                                       "supportDisplacement"])},
        "allOf": [
            variant("nodal", {
                "node": idref("Loaded node id."),
                **{k: num(f"Force along global {k[1]}.", "force", 0.0) for k in ("Fx", "Fy", "Fz")},
                **{k: num(f"Moment about global {k[1]} (right-hand rule).", "moment", 0.0) for k in ("Mx", "My", "Mz")},
            }, ["node"]),
            variant("memberPoint", {
                "member": idref("Loaded member id."),
                "direction": direction,
                "kind": enum("force = concentrated force along `direction`; moment = concentrated moment about "
                             "`direction`.", ["force", "moment"], "force"),
                "value": {"type": "number", "x-unit": "force",
                          "description": "Magnitude: a force, or a moment (force*length) when kind = moment."},
                "position": num("Location of the load along the member.", "position", minimum=0),
                "positionType": ptype,
            }, ["member", "direction", "value", "position"]),
            variant("memberDistributed", {
                "member": idref("Loaded member id."),
                "direction": direction,
                "w1": num("Intensity at `start` (force per unit length of member, or per unit projected "
                          "length when projected = true).", "force/length"),
                "w2": num("Intensity at `end`; linear variation between. Default: w1 (uniform).", "force/length",
                          examples=[-6]),
                "start": num("Start of the loaded length.", "position", 0.0, minimum=0),
                "end": num("End of the loaded length. Default: the member end.", "position", minimum=0,
                           examples=[0.5]),
                "positionType": ptype,
                "projected": flag("Global directions only: intensities are per unit length of the member's "
                                  "projection on the plane normal to the load direction (e.g. snow on a "
                                  "sloping rafter).", False),
            }, ["member", "direction", "w1"]),
            variant("selfWeight", {
                "factor": num("Multiplier on the weight (unitWeight * A per unit length of flexible member).",
                              "dimensionless", 1.0),
                "direction": vec3("Global direction of gravity (normalised by the engine). Default: the "
                                  "negative vertical axis.", "dimensionless", examples=[[0, 0, -1]]),
                "members": idlist("Members that receive self weight. Default: all members.",
                                  examples=[["B1", "B2"]]),
            }),
            variant("supportDisplacement", {
                "node": idref("Supported node id; every DOF given here must be \"fixed\" at that support."),
                **{k: num(f"Imposed translation along global {k[1].upper()}. Absent: this DOF is not imposed.",
                          "length", examples=[-0.01]) for k in ("ux", "uy", "uz")},
                **{k: num(f"Imposed rotation about global {k[1].upper()}. Absent: this DOF is not imposed.",
                          "rotation", examples=[0.002]) for k in ("rx", "ry", "rz")},
            }, ["node"]),
        ],
        "unevaluatedProperties": False,
    }
    # a point moment has dimension force*length: switch the unit of `value`
    defs["load"]["allOf"].append({
        "if": {"properties": {"type": {"const": "memberPoint"}, "kind": {"const": "moment"}},
               "required": ["type", "kind"]},
        "then": {"properties": {"value": {"x-unit": "moment"}}}})

    defs["loadCase"] = obj("A named set of loads that act together.", {
        "description": text("Free text. Absent: no text."),
        "loads": {"type": "array", "description": "The loads of this case.", "x-unit": "n/a",
                  "items": {"$ref": "#/$defs/load"}, "default": []},
    })
    defs["combination"] = obj(
        "Factored sum of load cases, analysed as one loading.",
        {
            "description": text("Free text. Absent: no text."),
            "factors": {"type": "object", "x-unit": "n/a", "minProperties": 1,
                        "description": "Load factor per load case id.",
                        "propertyNames": {"pattern": ID_PATTERN},
                        "additionalProperties": num("Load factor.", "dimensionless")},
            "analysisType": enum("Analysis type for this combination. Default: analysis.type.",
                                 ["firstOrder", "secondOrder"]),
            "purpose": enum("strength = used for strength design checks; serviceability = used for "
                            "deflection checks.", ["strength", "serviceability"], D.COMBINATION["purpose"]),
            "stiffnessReduction": flag("Apply analysis.stiffnessReduction (AISC 360-22 C2.3) in this "
                                       "combination. AFFECTS ANALYSIS.", D.COMBINATION["stiffnessReduction"]),
            "notionalLoad": obj(
                "Notional lateral loads added in this combination (AISC 360-22 C2.2b): at every node, "
                "coefficient * (factored gravity load applied at that node, including member loads lumped "
                "statically to the member ends), acting along `direction`. AFFECTS ANALYSIS. Absent: none.",
                {
                    "coefficient": num("Notional load coefficient.", "dimensionless", D.NOTIONAL["coefficient"],
                                       exclusiveMinimum=0),
                    "direction": vec3("Global horizontal direction of the notional loads (normalised by the "
                                      "engine; must be perpendicular to the vertical axis).", "dimensionless",
                                      examples=[[1, 0, 0]]),
                }, required=["direction"], examples=[{"direction": [1, 0, 0]}]),
        },
        required=["factors"])

    # -------------------------------------------------------------------- design
    def brace_list(what):
        return {"type": "array", "x-unit": "position", "uniqueItems": True,
                "description": f"Locations of intermediate points braced against {what}, measured from the "
                               "start of the design member (its ends are always braced points). Absent: no "
                               "intermediate brace point of this kind at this level.",
                "examples": [[0.5]],
                "items": {"type": "number", "minimum": 0, "x-unit": "position",
                          "description": "Brace location.", "examples": [0.5]}}

    def computed_or_number(desc, name):
        return {"description": desc, "x-unit": "dimensionless", "default": D.DESIGN[name][0],
                "x-hint": "Give a positive number or the string \"computed\".",
                "anyOf": [{"const": "computed"}, {"type": "number", "exclusiveMinimum": 0}]}

    dd = D.DESIGN
    defs["designVariables"] = {
        "type": "object", "x-unit": "n/a",
        "description": "Design variables (DESIGN ONLY: none of these can influence the analysis). Axis names "
                       "follow AISC: x = major axis (= member local z), y = minor axis (= member local y), "
                       "z = longitudinal axis (torsion). Each variable is resolved in the order: design "
                       "member > group(s) > design.defaults > built-in default.",
        "properties": {
            "code": enum("Design code and edition.", ["AISC 360-22"], dd["code"][0]),
            "method": enum("Design method.", ["LRFD", "ASD"], dd["method"][0]),
            "role": enum("Member role. Default: truss-type members = brace; frame members parallel to the "
                         "vertical axis = column; other frame members = beam.", ["beam", "column", "brace"]),
            "E": num("Modulus of elasticity used by the design checks. Default: E of the member's material. "
                     "An override is DESIGN ONLY: the analysis always uses the material's E.", "force/length^2",
                     examples=[200000000], **pos),
            "G": num("Shear modulus used by the design checks. Default: G of the member's material. An "
                     "override is DESIGN ONLY: the analysis always uses the material's G.", "force/length^2",
                     examples=[77000000], **pos),
            "Fy": num("Specified minimum yield stress used by the design checks. Default: Fy of the member's "
                      "material. DESIGN ONLY: a computed tau_b (analysis.stiffnessReduction) always uses the "
                      "material's Fy.", "force/length^2", examples=[345000], **pos),
            "Fu": num("Specified minimum tensile strength used by the design checks. Default: Fu of the "
                      "member's material.", "force/length^2", examples=[450000], **pos),
            "Lx": num("Unbraced length for flexural buckling about the major (x) axis. Default: from "
                      "bracing.major, else the design member length.", L, examples=[4.0], **pos),
            "Ly": num("Unbraced length for flexural buckling about the minor (y) axis. Default: from "
                      "bracing.minor, else the design member length.", L, examples=[2.0], **pos),
            "Lz": num("Unbraced length for torsional buckling. Default: from bracing.torsion, else the "
                      "design member length.", L, examples=[4.0], **pos),
            "Lb": num("Unbraced length for lateral-torsional buckling. Default: from bracing.topFlange / "
                      "bracing.bottomFlange, else the design member length.", L, examples=[2.0], **pos),
            "Kx": num("Effective length factor, major axis flexural buckling.", "dimensionless", dd["Kx"][0], **pos),
            "Ky": num("Effective length factor, minor axis flexural buckling.", "dimensionless", dd["Ky"][0], **pos),
            "Kz": num("Effective length factor, torsional buckling.", "dimensionless", dd["Kz"][0], **pos),
            "bracing": obj(
                "Intermediate bracing points along the design member. Within one level (member, group or "
                "defaults) an unbraced length and the matching bracing list must not both be given. Absent: no "
                "intermediate brace points at this level.",
                {
                    "positionType": enum("relative = fractions of the design member length; absolute = "
                                         "distances from its start.", ["relative", "absolute"], "relative"),
                    "major": brace_list("major-axis (x) flexural buckling"),
                    "minor": brace_list("minor-axis (y) flexural buckling"),
                    "torsion": brace_list("twist (torsional buckling)"),
                    "topFlange": brace_list("lateral movement of the top flange (the flange on the +local y "
                                            "side), for lateral-torsional buckling"),
                    "bottomFlange": brace_list("lateral movement of the bottom flange (the flange on the "
                                               "-local y side), for lateral-torsional buckling"),
                }),
            "Cb": computed_or_number("Lateral-torsional buckling modification factor: a number, or "
                                     "\"computed\" from the moment diagram of each unbraced segment "
                                     "(AISC 360-22 Eq. F1-1).", "Cb"),
            "Cmx": computed_or_number("Equivalent uniform moment factor Cm for major-axis bending "
                                      "(AISC 360-22 App. 8.1.2), a number or \"computed\".", "Cmx"),
            "Cmy": computed_or_number("Equivalent uniform moment factor Cm for minor-axis bending "
                                      "(AISC 360-22 App. 8.1.2), a number or \"computed\".", "Cmy"),
            "sideswayInhibited": flag("Sidesway flag: true when relative lateral translation of the member "
                                      "ends is prevented (braced frame).", dd["sideswayInhibited"][0]),
            "An": num("Net area for tension rupture. Default: the gross area (no holes).", "length^2",
                     examples=[0.0105], **pos),
            "U": num("Shear lag factor for tension rupture.", "dimensionless", dd["U"][0], exclusiveMinimum=0,
                     maximum=1),
            "maxSlendernessCompression": num("Limit on Lc/r for members in compression (advisory check).",
                                             "dimensionless", dd["maxSlendernessCompression"][0], **pos),
            "maxSlendernessTension": num("Limit on L/r for members in tension (advisory check).",
                                         "dimensionless", dd["maxSlendernessTension"][0], **pos),
            "deflectionLimits": {
                "type": "array", "x-unit": "n/a", "default": [],
                "description": "Deflection limits to check. Default: none.",
                "items": obj(
                    "One deflection limit; give spanRatio and/or absolute.",
                    {
                        "spanRatio": num("n in the limit L/n, L = design member length. Absent: only "
                                         "`absolute` limits.", "dimensionless", **pos),
                        "absolute": num("Limit as a length. Absent: only `spanRatio` limits.", L,
                                        examples=[0.025], **pos),
                        "direction": enum("Deflection component checked, in member local axes, measured "
                                          "relative to the chord through the design member ends.",
                                          ["y", "z", "resultant"], "resultant"),
                        "combinations": idlist("Combination ids the limit applies to. Default: every "
                                               "combination with purpose = serviceability."),
                    },
                    anyOf=[{"required": ["spanRatio"]}, {"required": ["absolute"]}],
                    **{"x-hint": "A deflection limit needs `spanRatio` and/or `absolute`."}),
            },
        },
    }
    dv = {"$ref": "#/$defs/designVariables"}
    design = {
        "type": "object", "x-unit": "n/a",
        "description": "Design input. Design-only data: nothing in this block can change analysis results. "
                       "Absent: every analysis member is a design member with the built-in defaults.",
        "additionalProperties": False,
        "properties": {
            "autoMembers": flag("true: every analysis member that is not part of an explicit design member "
                                "is a design member by itself, with the same id.", True),
            "defaults": {**dv, "description": "Design variables applied to every design member unless "
                                              "overridden by a group or the member. Absent: the built-in "
                                              "defaults.", "x-unit": "n/a",
                         "unevaluatedProperties": False},
            "groups": idmap("Design groups: variables shared by several design members. Absent: no groups.", {
                "type": "object", "x-unit": "n/a", "description": "A design group.",
                "allOf": [dv], "required": ["members"],
                "properties": {
                    "members": idlist("Design member ids in this group.", minItems=1),
                    "priority": integer("When a design member is in several groups that set the same "
                                        "variable, the group with the highest priority wins; equal priority "
                                        "with different values is an error.", 0),
                },
                "unevaluatedProperties": False}),
            "members": idmap("Design members. A design member is one analysis member or a chain of "
                             "collinear analysis members (a physical member split at nodes). Absent: only "
                             "the automatic design members of `autoMembers`.", {
                "type": "object", "x-unit": "n/a", "description": "A design member with its own variables.",
                "allOf": [dv],
                "properties": {
                    "members": idlist("Analysis member ids forming this design member, in order from its "
                                      "start to its end. Default: the single analysis member whose id equals "
                                      "this design member's id.", minItems=1),
                },
                "unevaluatedProperties": False}),
        },
    }

    # ------------------------------------------------------------------ analysis
    analysis = obj("Analysis options. Everything here AFFECTS ANALYSIS. Absent: every option takes its "
                   "default.", {
        "type": enum("Default analysis type of combinations: firstOrder = linear elastic; secondOrder = "
                     "geometrically nonlinear (geometric stiffness updated iteratively, P-Delta and "
                     "P-delta).", ["firstOrder", "secondOrder"], A["type"]),
        "shearDeformation": flag("Include shear deformation (Timoshenko beam) using the section shear areas.",
                                 A["shearDeformation"]),
        "resultStations": integer("Number of equally spaced result stations per member, ends included.",
                                  A["resultStations"], minimum=2, maximum=1001),
        "convergenceTolerance": num("Iterative analyses stop when the norm of the displacement change is "
                                    "below this fraction of the displacement norm.", "dimensionless",
                                    A["convergenceTolerance"], exclusiveMinimum=0),
        "maxIterations": integer("Maximum iterations (second-order, tension/compression-only).",
                                 A["maxIterations"], minimum=1),
        "alpha": num("Second-order analyses are run at alpha times the combination loads and the results "
                     "divided by alpha (AISC 360-22 C2.1(d): 1.0 for LRFD, 1.6 for ASD).", "dimensionless",
                     A["alpha"], exclusiveMinimum=0),
        "equilibriumTolerance": num("Largest relative out-of-balance force accepted in the invariants report.",
                                    "dimensionless", A["equilibriumTolerance"], exclusiveMinimum=0),
        "pivotTolerance": num("Relative pivot size below which the stiffness matrix is declared singular "
                              "(mechanism).", "dimensionless", A["pivotTolerance"], exclusiveMinimum=0),
        "geometryTolerance": num("Lengths smaller than this fraction of the model size count as zero.",
                                 "dimensionless", A["geometryTolerance"], exclusiveMinimum=0),
        "buckling": obj("Linear (eigenvalue) buckling analysis. Absent: not performed.", {
            "combinations": idlist("Reference combinations: buckling factors multiply these loadings.", minItems=1),
            "modes": integer("Number of buckling modes per reference combination.", D.BUCKLING["modes"], minimum=1),
        }, required=["combinations"]),
        "stiffnessReduction": obj(
            "Stiffness reduction for stability design (AISC 360-22 C2.3). Used only in combinations with "
            "stiffnessReduction = true. Absent: `factor` and `tauB` take their defaults and every member is "
            "reduced.", {
                "factor": num("Factor applied to all member stiffnesses (EA, EI, GJ, GAs).", "dimensionless",
                              SR["factor"], exclusiveMinimum=0, maximum=1),
                "tauB": {"description": "Additional factor tau_b on flexural stiffness: \"computed\" (Eq. C2-2a/b "
                                        "from the member axial force, Pns = Fy*A), \"none\" (1.0), or a number.",
                         "x-unit": "dimensionless", "default": SR["tauB"],
                         "x-hint": "Use \"computed\", \"none\" or a number in (0, 1].",
                         "anyOf": [{"enum": ["computed", "none"]},
                                   {"type": "number", "exclusiveMinimum": 0, "maximum": 1}]},
                "members": idlist("Members whose stiffness is reduced. Default: all members.",
                                  examples=[["C1", "C2"]]),
            }, examples=[{"tauB": "computed"}]),
    })

    conv = D.CONVENTIONS
    conventions = obj(
        "Sign and axis conventions of this file. Only the listed values are supported; anything else is "
        "refused. Explained in docs/CONVENTIONS.md.",
        {
            "verticalAxis": enum("Global axis that points up. Sets the default member orientation, the "
                                 "default gravity direction (negative vertical) and the default diaphragm "
                                 "normal.", ["Y", "Z"]),
            "handedness": enum("Global and local axes are right-handed.", [conv["handedness"]], conv["handedness"]),
            "rotations": enum("Rotations and moments are positive by the right-hand rule about the positive "
                              "axis.", [conv["rotations"]], conv["rotations"]),
            "memberLocalX": enum("Member local x runs from node i to node j.", [conv["memberLocalX"]],
                                 conv["memberLocalX"]),
            "defaultOrientation": enum("Default rule for local y and z (see member.orientation).",
                                       [conv["defaultOrientation"]], conv["defaultOrientation"]),
            "memberEndForces": enum("Member end forces are the forces and moments the joints exert on the "
                                    "member ends, in member local axes.", [conv["memberEndForces"]],
                                    conv["memberEndForces"]),
            "internalForces": enum("Internal forces are the stress resultants acting on the cross-section "
                                   "face whose outward normal is +local x, as components along the local "
                                   "axes (axial tension positive).", [conv["internalForces"]],
                                   conv["internalForces"]),
            "reactions": enum("Reactions are the forces and moments the supports exert on the structure, "
                              "in global axes.", [conv["reactions"]], conv["reactions"]),
            "loads": enum("Load values are positive along (forces) or about (moments) the positive axis "
                          "of their direction.", [conv["loads"]], conv["loads"]),
        },
        required=["verticalAxis"])

    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": SCHEMA_ID,
        "title": "BenchFEA input file, format 1.0",
        "description": "Everything needed to analyse and design a 3D frame structure. Collections are objects "
                       "keyed by id, so JSON pointers such as /members/B1/section stay valid when items are "
                       "added or removed. All dimensional numbers are in the consistent units built from "
                       "units.force and units.length.",
        "type": "object", "x-unit": "n/a",
        "additionalProperties": False,
        "required": ["formatVersion", "units", "conventions", "nodes", "materials", "sections", "members",
                     "loadCases", "combinations"],
        "properties": {
            "$schema": text("Optional pointer to this schema for editors. Absent: nothing changes (the engine "
                            "does not read it).", examples=[SCHEMA_ID]),
            "formatVersion": text("Version of the file format (semantic versioning, see "
                                  "docs/FORMAT_VERSIONING.md).", pattern=r"^\d+\.\d+\.\d+$"),
            "name": text("Short model name. Absent: no name."),
            "description": text("Free text. Absent: no text."),
            "units": obj("Units of every dimensional number in the file.", {
                "force": enum("Force unit.", ["N", "kN", "lbf", "kip"]),
                "length": enum("Length unit.", ["mm", "m", "in", "ft"]),
            }, required=["force", "length"]),
            "conventions": conventions,
            "nodes": idmap("Nodes, keyed by node id.", obj("Node position in global axes.", {
                "x": num("Global X coordinate.", L), "y": num("Global Y coordinate.", L),
                "z": num("Global Z coordinate.", L, 0.0)}, required=["x", "y"]), minProperties=1),
            "supports": idmap("Supports, keyed by the id of the supported node.", {"$ref": "#/$defs/support"},
                              default={}),
            "materials": idmap("Materials, keyed by id.", {"$ref": "#/$defs/material"}),
            "sections": idmap("Sections, keyed by id.", {"$ref": "#/$defs/section"}),
            "members": idmap("Members, keyed by id.", {"$ref": "#/$defs/member"}, minProperties=1),
            "constraints": idmap("Rigid links and rigid diaphragms, keyed by id.", {"$ref": "#/$defs/constraint"},
                                 default={}),
            "loadCases": idmap("Load cases, keyed by id.", {"$ref": "#/$defs/loadCase"}),
            "combinations": idmap("Load combinations, keyed by id. Results are produced per combination.",
                                  {"$ref": "#/$defs/combination"}, minProperties=1),
            "analysis": analysis,
            "design": design,
        },
        "$defs": defs,
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Write the JSON Schema of the input format from its generator.")
    ap.add_argument("--check", action="store_true", help="do not write; exit 1 if the file on disk differs")
    args = ap.parse_args(argv)
    text_ = json.dumps(build(), indent=2) + "\n"
    if args.check:
        if not SCHEMA_FILE.exists() or SCHEMA_FILE.read_text(encoding="utf-8") != text_:
            print(f"{SCHEMA_FILE.name} is not up to date: run python -m benchfea.schema.build_schema")
            return 1
        print(f"{SCHEMA_FILE.name} is up to date")
        return 0
    SCHEMA_FILE.write_text(text_, encoding="utf-8")
    print(f"wrote {SCHEMA_FILE}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
