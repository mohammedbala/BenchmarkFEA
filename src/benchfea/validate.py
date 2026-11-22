"""Validation of input documents: collect EVERY problem, with its location and a hint.

Two layers:

1. schema   -- the document against the published JSON Schema (draft 2020-12);
2. semantic -- cross-references, geometry, releases, constraints, loads,
               combinations, design groups (things a schema cannot express).

Both layers always run; semantic checks skip only the items too malformed to
inspect.  Stability of the assembled structure (global mechanisms) is the
solver's job: it raises :class:`benchfea.problems.ModelUnstable` with problems
of code ``stability.mechanism`` (see docs/ARCHITECTURE.md).

Public API::

    problems = validate(doc)          # list[Problem], errors and warnings
    warnings = check(doc)             # raises ModelInvalid if any error
"""
from __future__ import annotations

import difflib
import math
import re
from contextlib import contextmanager

import numpy as np

from . import defaults as D
from . import designparams, sections, shapes
from .geometry import AXIS, GeometryError, is_vertical, member_geometry, model_size
from .io import Document, load
from .problems import ERROR, WARNING, ModelInvalid, ModelUnstable, Problem, pointer
from .schema import load_schema, validator

__all__ = ["Problem", "ModelInvalid", "ModelUnstable", "validate", "check", "validate_schema",
           "validate_semantics", "ERROR", "WARNING"]

DIAPHRAGM_DOFS = {"X": ("uy", "uz", "rx"), "Y": ("ux", "uz", "ry"), "Z": ("ux", "uy", "rz")}
TRUSS_TYPES = ("truss", "tensionOnly", "compressionOnly")


# =========================================================================== schema layer
def _jtype(v) -> str:
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "boolean"
    if isinstance(v, (int, float)):
        return "number"
    if isinstance(v, str):
        return "string"
    return "array" if isinstance(v, list) else "object"


def _deref(schema):
    root = load_schema()
    while isinstance(schema, dict) and "$ref" in schema:
        target = root
        for part in schema["$ref"].lstrip("#/").split("/"):
            target = target[part]
        schema = target
    return schema


def _known_properties(schema, depth=0) -> dict:
    """All property names (with their schemas) a (sub)schema can accept."""
    schema = _deref(schema)
    out = {}
    if not isinstance(schema, dict) or depth > 6:
        return out
    out.update(schema.get("properties", {}))
    for key in ("allOf", "anyOf", "oneOf"):
        for sub in schema.get(key, ()):
            out.update(_known_properties(sub, depth + 1))
    for key in ("then", "else"):
        if key in schema:
            out.update(_known_properties(schema[key], depth + 1))
    return out


def _applicable_properties(schema, inst, depth=0) -> dict:
    """Property names that apply to ``inst``: own properties plus those of the
    allOf members and of the if/then variants whose condition ``inst`` meets."""
    schema = _deref(schema)
    out = {}
    if not isinstance(schema, dict) or depth > 6:
        return out
    out.update(schema.get("properties", {}))
    for sub in schema.get("allOf", ()):
        out.update(_applicable_properties(sub, inst, depth + 1))
    if "if" in schema:
        ok = validator().evolve(schema={**schema["if"], "$defs": load_schema()["$defs"]}).is_valid(inst)
        branch = schema.get("then" if ok else "else")
        if branch:
            out.update(_applicable_properties(branch, inst, depth + 1))
    return out


def _describe(schema) -> str:
    schema = _deref(schema) if isinstance(schema, dict) else {}
    desc = schema.get("description", "")
    unit = schema.get("x-unit")
    if unit and unit != "n/a":
        desc += f" [unit: {unit}]"
    return desc.strip()


def _schema_error_problems(err, out: list, suppress: set):
    path = pointer(*err.absolute_path)
    kw, inst, schema = err.validator, err.instance, err.schema

    def add(code, p, message, hint=""):
        out.append(Problem(p, message, hint, ERROR, "schema." + code))

    if kw == "required":
        name = re.match(r"'(.*)' is a required property", err.message).group(1)
        hint = f"Add `{name}`"
        desc = _describe(schema.get("properties", {}).get(name) or _known_properties(schema).get(name) or {})
        add("required", path, f"missing required property '{name}'", hint + (f": {desc}" if desc else "."))
        if name == "type":
            suppress.add(path)
    elif kw in ("additionalProperties", "unevaluatedProperties"):
        # Work out the applicable property names ourselves: jsonschema reports every
        # property as unevaluated as soon as the matching variant has any other error.
        known = _applicable_properties(schema, inst)
        for name in [k for k in inst if k not in known]:
            near = difflib.get_close_matches(name, list(known), n=1)
            hint = f"Did you mean '{near[0]}'? " if near else ""
            if known:
                hint += "Allowed here: " + ", ".join(sorted(known)) + "."
            out.append(Problem(path + pointer(name), f"unknown property '{name}'", hint, ERROR,
                               "schema.unknown-property"))
            if kw == "unevaluatedProperties":
                suppress.add(id(out[-1]))      # may be a mere consequence of a bad `type` discriminator
    elif "propertyNames" in err.absolute_schema_path:
        out.append(Problem(path + pointer(inst), f"invalid id {inst!r}", "Ids are 1-64 characters: letters, "
                           "digits, '_', '.', '-' (not starting with '.' or '-').", ERROR, "id.invalid"))
    elif kw == "pattern":
        if path == "/formatVersion":
            add("pattern", path, f"formatVersion {inst!r} is not of the form MAJOR.MINOR.PATCH",
                "Use a semantic version such as \"1.0.0\" (docs/FORMAT_VERSIONING.md).")
        else:
            add("pattern", path, f"{inst!r} is not a valid id", "Ids are 1-64 characters: letters, digits, "
                "'_', '.', '-' (not starting with '.' or '-').")
    elif kw in ("enum", "const"):
        allowed = schema["enum"] if kw == "enum" else [schema["const"]]
        near = difflib.get_close_matches(str(inst), [str(a) for a in allowed], n=1)
        hint = "Allowed: " + ", ".join(repr(a) for a in allowed) + "."
        if len(allowed) == 1 and path.startswith("/conventions"):
            hint = f"The engine supports only {allowed[0]!r} for this convention; convert the model to it."
        elif near:
            hint = f"Did you mean {near[0]!r}? " + hint
        add("enum", path, f"value {inst!r} is not allowed here", hint)
        if path.endswith("/type"):
            suppress.add(path[:-5])
    elif kw == "type":
        hint = _describe(schema)
        if schema["type"] == "string" and _jtype(inst) in ("number", "boolean", "null"):
            hint = (hint + " " if hint else "") + "Write the value as a string; in YAML put it in quotes " \
                "(unquoted 1.0, yes, no, on, off, null are not strings)."
        add("type", path, f"expected {schema['type']}, got {_jtype(inst)}", hint)
    elif kw in ("minimum", "exclusiveMinimum", "maximum", "exclusiveMaximum"):
        op = {"minimum": ">=", "exclusiveMinimum": ">", "maximum": "<=", "exclusiveMaximum": "<"}[kw]
        add("range", path, f"value {inst!r} must be {op} {schema[kw]}", _describe(schema))
    elif kw in ("anyOf", "oneOf"):
        branches = schema[kw]
        for k, sub in enumerate(branches):
            sub = _deref(sub)
            if isinstance(inst, dict) and sub.get("type") == "object" and len(
                    [b for b in branches if _deref(b).get("type") == "object"]) == 1:
                for e in err.context:
                    if e.relative_schema_path[0] == k:
                        _schema_error_problems(e, out, suppress)
                return
        if all(set(_deref(b)) <= {"required"} for b in branches):
            names = " or ".join("`" + "`+`".join(b["required"]) + "`" for b in branches)
            add("required", path, f"missing property: one of {names} is required", schema.get("x-hint", ""))
        else:
            add("invalid", path, f"value {inst!r} is not valid here" if not isinstance(inst, (dict, list))
                else "value is not valid here", schema.get("x-hint", _describe(schema)))
    elif kw in ("minItems", "maxItems"):
        lo, hi = schema.get("minItems"), schema.get("maxItems")
        want = f"exactly {lo}" if lo == hi else (f"at least {lo}" if kw == "minItems" else f"at most {hi}")
        add(kw, path, f"list has {len(inst)} item(s), needs {want}", _describe(schema))
    elif kw == "uniqueItems":
        dup = sorted({repr(v) for k, v in enumerate(inst) if v in inst[:k]})
        add(kw, path, f"list has repeated item(s): {', '.join(dup)}", "Remove the repetition. " + _describe(schema))
    elif kw == "minProperties":
        add(kw, path, f"needs at least {schema[kw]} entr{'y' if schema[kw] == 1 else 'ies'}, has {len(inst)}",
            _describe(schema))
    else:
        add(kw, path, err.message, _describe(schema))


def validate_schema(data) -> list[Problem]:
    """Layer 1. Returns unlocated problems (``validate`` adds line/column)."""
    out: list[Problem] = []
    suppress: set[str] = set()
    for err in validator().iter_errors(data):
        _schema_error_problems(err, out, suppress)
    # an object whose `type` discriminator is wrong/missing would also list every
    # other property as unknown: keep only the root cause
    return [p for p in out if not (id(p) in suppress and p.path.rsplit("/", 1)[0] in suppress)]


# ========================================================================= semantic layer
class _Semantics:
    def __init__(self, data: dict, strict: bool):
        self.d = data
        self.strict = strict
        self.out: list[Problem] = []
        self.geoms: dict[str, tuple] = {}

    def add(self, code, path, message, hint="", severity=ERROR):
        self.out.append(Problem(path, message, hint, severity, code))

    @contextmanager
    def guard(self):
        """Skip an item that is too malformed to inspect (the schema layer has
        already reported it); in strict mode (schema-valid input) re-raise."""
        try:
            yield
        except (KeyError, TypeError, AttributeError, ValueError, IndexError):
            if self.strict:
                raise

    def table(self, name) -> dict:
        t = self.d.get(name) if isinstance(self.d, dict) else None
        return t if isinstance(t, dict) else {}

    def items(self, name):
        for key, val in self.table(name).items():
            if isinstance(val, dict):
                yield key, val, pointer(name, key)

    # ------------------------------------------------------------------ run
    def run(self) -> list[Problem]:
        with self.guard():
            self.prepare()
        for step in (self.nodes, self.supports, self.materials, self.sections, self.members, self.constraints,
                     self.loads, self.combinations, self.analysis, self.design):
            with self.guard():
                step()
        return self.out

    def prepare(self):
        self.vertical = "Z"
        with self.guard():
            v = self.d["conventions"]["verticalAxis"]
            if v in ("Y", "Z"):
                self.vertical = v
        self.node_ok = {k: n for k, n in self.table("nodes").items()
                        if isinstance(n, dict) and all(isinstance(n.get(c, 0.0), (int, float))
                                                       and not isinstance(n.get(c, 0.0), bool)
                                                       and math.isfinite(n.get(c, 0.0)) for c in "xyz")
                        and "x" in n and "y" in n}
        self.size = model_size(self.node_ok) if self.node_ok else 0.0
        gt = D.ANALYSIS["geometryTolerance"]
        with self.guard():
            g = self.d.get("analysis", {}).get("geometryTolerance", gt)
            if isinstance(g, (int, float)) and g > 0:
                gt = g
        self.tol = gt * self.size

    def xyz(self, nid):
        n = self.node_ok[nid]
        return np.array([n["x"], n["y"], n.get("z", 0.0)], float)

    # ---------------------------------------------------------------- nodes
    def nodes(self):
        seen: dict[tuple, str] = {}
        q = self.tol if self.tol > 0 else 1.0
        for nid in self.node_ok:
            key = tuple(np.round(self.xyz(nid) / q).astype(np.int64)) if self.tol > 0 else tuple(self.xyz(nid))
            if key in seen:
                self.add("node.coincident", pointer("nodes", nid),
                         f"node {nid!r} has the same position as node {seen[key]!r}",
                         "Coincident nodes are not connected to each other; merge them or connect them with a "
                         "member or constraint if that is intended.", WARNING)
            else:
                seen[key] = nid
        used = set()
        for _, m, _ in self.items("members"):
            used.update(v for v in (m.get("i"), m.get("j")) if isinstance(v, str))
        for _, c, _ in self.items("constraints"):
            if isinstance(c.get("master"), str):
                used.add(c["master"])
            used.update(s for s in c.get("slaves") or () if isinstance(s, str))
        for nid in self.table("nodes"):
            if nid not in used:
                self.add("node.unused", pointer("nodes", nid), f"node {nid!r} is not used by any member or constraint",
                         "Remove it, or connect it; an unconnected, unrestrained node makes the structure singular.",
                         WARNING)

    # ------------------------------------------------------------- supports
    def supports(self):
        nodes = self.table("nodes")
        restrained = False
        for nid, s, p in self.items("supports"):
            if nid not in nodes:
                self.add("ref.unknown-node", p, f"support refers to unknown node {nid!r}",
                         "Supports are keyed by the id of the supported node; use an id from /nodes.")
            with self.guard():
                if all(v == "free" for v in self.support_dofs(s).values()):
                    self.add("support.no-restraint", p, f"support at node {nid!r} restrains nothing",
                             "Give `type` or set at least one DOF to \"fixed\" or a spring.", WARNING)
                else:
                    restrained = True
        raw = self.d.get("supports", {})
        if isinstance(raw, dict) and not restrained and all(isinstance(v, dict) for v in raw.values()):
            self.add("support.none", "/supports" if "supports" in self.d else "",
                     "the structure has no support: no DOF of any node is fixed or on a spring",
                     "Add `supports` (keyed by node id), e.g. {\"N1\": {\"type\": \"fixed\"}}; an unsupported "
                     "structure is a mechanism.")

    def support_dofs(self, s: dict) -> dict:
        """DOF name -> "fixed" | "free" | {"spring": k} with the shorthand expanded."""
        v = {"X": "ux", "Y": "uy", "Z": "uz"}[self.vertical]
        base = {"fixed": D.DOFS, "pinned": D.DOFS[:3], "roller": (v,)}.get(s.get("type"), ())
        out = {dof: ("fixed" if dof in base else "free") for dof in D.DOFS}
        out.update(s.get("dofs") or {})
        return out

    # ------------------------------------------------------------ materials
    def materials(self):
        for mid, m, p in self.items("materials"):
            with self.guard():
                if "G" in m and "nu" in m:
                    g = m["E"] / (2.0 * (1.0 + m["nu"]))
                    if abs(g - m["G"]) > 0.01 * g:
                        self.add("material.g-nu-inconsistent", p + "/G",
                                 f"G = {m['G']:g} differs from E/(2(1+nu)) = {g:g} by more than 1%",
                                 "G is used for stiffness and nu is ignored; remove one of them if unintended.",
                                 WARNING)

    # ------------------------------------------------------------- sections
    def sections(self):
        for sid, s, p in self.items("sections"):
            with self.guard():
                for prop, msg in sections.check_dimensions(s):
                    self.add("section.geometry", f"{p}/{prop}", f"section {sid!r}: {msg}",
                             "Check the dimensions (d = overall depth, bf = flange width, H/B/D = outside "
                             "dimensions).")
            if s.get("type") == "standard" and isinstance(s.get("shape"), str):
                try:
                    shapes.lookup_raw(s["shape"])
                except shapes.ShapesDatabaseError as e:
                    self.add("section.shapes-db", p + "/shape", f"section {sid!r}: {e}", e.hint)

    # -------------------------------------------------------------- members
    def members(self):
        nodes, secs, mats = self.table("nodes"), self.table("sections"), self.table("materials")
        for mid, m, p in self.items("members"):
            ok = True
            for end in ("i", "j"):
                if isinstance(m.get(end), str) and m[end] not in nodes:
                    ok = False
                    self.add("ref.unknown-node", f"{p}/{end}", f"member {mid!r} refers to unknown node {m[end]!r}",
                             "Use an id from /nodes.")
            if isinstance(m.get("section"), str) and m["section"] not in secs:
                self.add("ref.unknown-section", p + "/section",
                         f"member {mid!r} refers to unknown section {m['section']!r}", "Use an id from /sections.")
            if isinstance(m.get("material"), str) and m["material"] not in mats:
                self.add("ref.unknown-material", p + "/material",
                         f"member {mid!r} refers to unknown material {m['material']!r}", "Use an id from /materials.")
            mtype = m.get("type", "frame")
            with self.guard():
                rel = m.get("releases") or {}
                if mtype in TRUSS_TYPES and any(rel.get(e) for e in ("i", "j")):
                    self.add("member.release-on-truss", p + "/releases",
                             f"member {mid!r} is of type {mtype} and cannot have releases",
                             "Truss-type members are pin-ended by definition; remove `releases` or use type frame.")
                elif mtype == "frame":
                    self.release_mechanisms(mid, rel, p)
            if not ok or m.get("i") not in self.node_ok or m.get("j") not in self.node_ok:
                continue
            with self.guard():
                if m["i"] == m["j"]:
                    self.add("member.same-node", p + "/j", f"member {mid!r} starts and ends at node {m['i']!r}",
                             "A member needs two different end nodes.")
                    continue
                try:
                    g = member_geometry(m, self.node_ok, self.vertical, self.tol)
                except GeometryError as e:
                    where = {"member.zero-length": p, "member.offset-zero-length": p + "/offsets",
                             "member.orientation-collinear": p + "/orientation/referenceVector"}[e.code]
                    hint = {"member.zero-length": "The end nodes coincide; move a node or remove the member.",
                            "member.offset-zero-length": "Shorten the rigid end offsets.",
                            "member.orientation-collinear": "Choose a reference vector that is not parallel to "
                                                            "the member (it fixes the plane of local x and y)."}[e.code]
                    self.add(e.code, where, f"member {mid!r}: {e}", hint)
                    continue
                self.geoms[mid] = g
                x = g[2][0]
                sine = float(np.linalg.norm(np.cross(x, AXIS[self.vertical])))
                if "referenceVector" not in (m.get("orientation") or {}) and not is_vertical(x, self.vertical) \
                        and sine < 1.0e-3:
                    self.add("member.near-vertical", p, f"member {mid!r} is almost but not exactly vertical; its "
                             "default local y direction is then decided by a very small horizontal offset",
                             "Give orientation.referenceVector explicitly, or make the member exactly vertical.",
                             WARNING)

    def release_mechanisms(self, mid, rel, p):
        def retained(end, dof):
            v = (rel.get(end) or {}).get(dof, False)
            if isinstance(v, dict):
                return v.get("spring", 0) > 0
            return v is False
        # rigid-body modes of the member relative to its end nodes; rows = retained connections
        blocks = {
            "axial translation (ux released at both ends)": ([("ux", [1.0])], 1),
            "rotation about its own axis (rx released at both ends)": ([("rx", [1.0])], 1),
            "rigid-body motion in the local x-y plane (uy / rz releases)": ([("uy", None), ("rz", [0.0, 1.0])], 2),
            "rigid-body motion in the local x-z plane (uz / ry releases)": ([("uz", None), ("ry", [0.0, 1.0])], 2),
        }
        for what, (spec, nmodes) in blocks.items():
            rows = []
            for end, x in (("i", 0.0), ("j", 1.0)):
                for dof, row in spec:
                    if retained(end, dof):
                        rows.append(row if row is not None else [1.0, x])
            rank = np.linalg.matrix_rank(np.array(rows, float).reshape(-1, nmodes)) if rows else 0
            if rank < nmodes:
                self.add("member.release-mechanism", p + "/releases",
                         f"the releases of member {mid!r} make it a mechanism: {what}",
                         "Keep enough end connections to hold the member: e.g. do not release the same "
                         "translation at both ends, or both end moments together with a shear.")

    # ---------------------------------------------------------- constraints
    def constraints(self):
        nodes = self.table("nodes")
        slave_of: dict[str, str] = {}
        masters = {c.get("master"): cid for cid, c, _ in self.items("constraints") if isinstance(c.get("master"), str)}
        for cid, c, p in self.items("constraints"):
            master = c.get("master")
            if isinstance(master, str) and master not in nodes:
                self.add("ref.unknown-node", p + "/master", f"constraint {cid!r} refers to unknown node {master!r}",
                         "Use an id from /nodes.")
            normal = c.get("normal", self.vertical)
            slaved = D.DOFS if c.get("type") == "rigidLink" else DIAPHRAGM_DOFS.get(normal, ())
            for k, s in enumerate(c.get("slaves") or ()):
                sp = f"{p}/slaves/{k}"
                if not isinstance(s, str) or s in c["slaves"][:k]:
                    continue                    # malformed or repeated (reported by the schema layer)
                if s not in nodes:
                    self.add("ref.unknown-node", sp, f"constraint {cid!r} refers to unknown node {s!r}",
                             "Use an id from /nodes.")
                    continue
                if s == master:
                    self.add("constraint.master-is-slave", sp, f"node {s!r} is both master and slave of "
                             f"constraint {cid!r}", "Remove the master from `slaves`.")
                    continue
                if s in slave_of:
                    self.add("constraint.slave-two-masters", sp, f"node {s!r} is already a slave in constraint "
                             f"{slave_of[s]!r}", "A node can be the slave of only one constraint; merge the "
                             "constraints under one master.")
                else:
                    slave_of[s] = cid
                if s in masters and masters[s] != cid:
                    self.add("constraint.chained", sp, f"node {s!r} is a slave here but the master of constraint "
                             f"{masters[s]!r}", "Chained constraints are not supported; use one master for all.")
                sup = self.table("supports").get(s)
                if isinstance(sup, dict):
                    with self.guard():
                        hit = [d for d, v in self.support_dofs(sup).items() if d in slaved and v != "free"]
                        if hit:
                            self.add("constraint.slave-restrained", sp,
                                     f"slave node {s!r} is also supported in constrained DOF(s) {', '.join(hit)}",
                                     "A constrained DOF follows the master and cannot be restrained separately; "
                                     "support the master instead, or free these DOFs.")
                if c.get("type") == "rigidDiaphragm" and normal in AXIS and s in self.node_ok \
                        and master in self.node_ok:
                    if abs(float((self.xyz(s) - self.xyz(master)) @ AXIS[normal])) > max(self.tol, 0.0):
                        self.add("constraint.diaphragm-out-of-plane", sp,
                                 f"slave node {s!r} is not in the diaphragm plane of master {master!r}",
                                 "Allowed (only the in-plane offset is used), but check the node levels.", WARNING)

    # ---------------------------------------------------------------- loads
    def loads(self):
        nodes, members = self.table("nodes"), self.table("members")
        used_cases = set()
        for _, c, _ in self.items("combinations"):
            if isinstance(c.get("factors"), dict):
                used_cases.update(c["factors"])
        for cid, case, cp in self.items("loadCases"):
            if cid not in used_cases:
                self.add("case.unused", cp, f"load case {cid!r} is not used by any combination",
                         "Results are produced per combination; add a combination that uses it.", WARNING)
            loads = case.get("loads")
            for k, ld in enumerate(loads if isinstance(loads, list) else ()):
                if isinstance(ld, dict):
                    with self.guard():
                        self.load(ld, f"{cp}/loads/{k}", nodes, members)

    def load(self, ld, p, nodes, members):
        t = ld.get("type")
        if t in ("nodal", "supportDisplacement"):
            nid = ld.get("node")
            if isinstance(nid, str) and nid not in nodes:
                self.add("ref.unknown-node", p + "/node", f"load refers to unknown node {nid!r}",
                         "Use an id from /nodes.")
                return
        if t == "supportDisplacement":
            sup = self.table("supports").get(ld.get("node"))
            dofs = self.support_dofs(sup) if isinstance(sup, dict) else {}
            given = [d for d in D.DOFS if d in ld]
            if not given:
                self.add("load.empty", p, "supportDisplacement load gives no displacement",
                         "Give at least one of ux, uy, uz, rx, ry, rz.", WARNING)
            for dof in given:
                if dofs.get(dof) != "fixed":
                    self.add("load.displacement-dof-not-fixed", f"{p}/{dof}",
                             f"imposed displacement {dof} at node {ld.get('node')!r} needs that DOF to be "
                             "\"fixed\" at a support", "Add a support with this DOF fixed (springs and free DOFs "
                             "cannot take an imposed displacement).")
        if t in ("memberPoint", "memberDistributed"):
            mid = ld.get("member")
            if isinstance(mid, str) and mid not in members:
                self.add("ref.unknown-member", p + "/member", f"load refers to unknown member {mid!r}",
                         "Use an id from /members.")
                return
            L = self.geoms[mid][3] if mid in self.geoms else None
            absolute = ld.get("positionType", "relative") == "absolute"
            lim = L if absolute else 1.0
            what = "the flexible member length" if absolute else "1 (positions are relative)"
            names = ("position",) if t == "memberPoint" else ("start", "end")
            if lim is not None:
                for name in names:
                    if name in ld and ld[name] > lim * (1 + 1e-12):
                        self.add("load.position-range", f"{p}/{name}",
                                 f"{name} = {ld[name]:g} is beyond {what}" + (f" = {lim:g}" if absolute else ""),
                                 "Positions run from 0 at end i to the member length (absolute) or 1 (relative); "
                                 "check positionType.")
            if t == "memberDistributed":
                end = ld.get("end", lim)
                if end is not None and not ld.get("start", 0.0) < end:
                    self.add("load.position-range", p + ("/end" if "end" in ld else "/start"),
                             f"loaded length is empty: start = {ld.get('start', 0.0):g} is not before end = {end:g}",
                             "Give start < end.")
                if ld.get("projected") and ld.get("direction") in ("x", "y", "z"):
                    self.add("load.projected-local", p + "/projected",
                             "projected = true is only meaningful for global directions (X, Y, Z)",
                             "Use a global direction, or remove `projected`.")
        if t == "selfWeight":
            if "direction" in ld and not np.linalg.norm(np.asarray(ld["direction"], float)) > 0:
                self.add("load.zero-direction", p + "/direction", "self-weight direction is the zero vector",
                         "Give the direction of gravity, e.g. [0, 0, -1].")
            targets = ld.get("members")
            for k, mid in enumerate(targets or ()):
                if mid not in members:
                    self.add("ref.unknown-member", f"{p}/members/{k}",
                             f"selfWeight refers to unknown member {mid!r}", "Use an id from /members.")
            mats = self.table("materials")
            zero = sorted({m.get("material") for mid, m in members.items()
                           if isinstance(m, dict) and (targets is None or mid in targets)
                           and isinstance(mats.get(m.get("material")), dict)
                           and not mats[m["material"]].get("unitWeight")})
            if zero:
                self.add("load.self-weight-zero", p, "selfWeight acts on members whose material has no unitWeight: "
                         + ", ".join(map(repr, zero)), "Set `unitWeight` (force per volume) on these materials.",
                         WARNING)

    # --------------------------------------------------------- combinations
    def combinations(self):
        cases = self.table("loadCases")
        an = self.table("analysis")
        for cid, c, p in self.items("combinations"):
            for case in (c.get("factors") if isinstance(c.get("factors"), dict) else {}):
                if case not in cases:
                    self.add("ref.unknown-case", f"{p}/factors{pointer(case)}",
                             f"combination {cid!r} refers to unknown load case {case!r}", "Use an id from /loadCases.")
            with self.guard():
                nl = c.get("notionalLoad")
                if nl:
                    d = np.asarray(nl["direction"], float)
                    n = float(np.linalg.norm(d))
                    if not n > 0:
                        self.add("combination.notional-direction", p + "/notionalLoad/direction",
                                 "notional load direction is the zero vector", "Give a horizontal direction.")
                    elif abs(float(d @ AXIS[self.vertical])) > 1e-9 * n:
                        self.add("combination.notional-direction", p + "/notionalLoad/direction",
                                 f"notional load direction is not horizontal (vertical axis is {self.vertical})",
                                 "Notional loads are lateral loads: the direction must be perpendicular to the "
                                 "vertical axis.")
            with self.guard():
                if c.get("stiffnessReduction") and c.get("analysisType", an.get("type", "firstOrder")) != "secondOrder":
                    self.add("combination.reduction-first-order", p + "/stiffnessReduction",
                             f"combination {cid!r} uses stiffness reduction with a first-order analysis",
                             "AISC 360-22 C2 pairs the reduced stiffness with a second-order analysis; set "
                             "analysisType: secondOrder if that is intended.", WARNING)

    # ------------------------------------------------------------- analysis
    def analysis(self):
        an = self.table("analysis")
        combos, members, mats = self.table("combinations"), self.table("members"), self.table("materials")
        with self.guard():
            for k, cid in enumerate((an.get("buckling") or {}).get("combinations") or ()):
                if cid not in combos:
                    self.add("ref.unknown-combination", f"/analysis/buckling/combinations/{k}",
                             f"buckling analysis refers to unknown combination {cid!r}", "Use an id from /combinations.")
        with self.guard():
            sr = an.get("stiffnessReduction") or {}
            targets = sr.get("members")
            for k, mid in enumerate(targets or ()):
                if mid not in members:
                    self.add("ref.unknown-member", f"/analysis/stiffnessReduction/members/{k}",
                             f"stiffness reduction refers to unknown member {mid!r}", "Use an id from /members.")
            active = any(isinstance(c, dict) and c.get("stiffnessReduction") for c in combos.values())
            if active and sr.get("tauB", D.STIFFNESS_REDUCTION["tauB"]) == "computed":
                need = {m.get("material") for mid, m in members.items() if isinstance(m, dict)
                        and m.get("type", "frame") == "frame" and (targets is None or mid in targets)}
                for mat in sorted(x for x in need if isinstance(mats.get(x), dict) and "Fy" not in mats[x]):
                    self.add("analysis.tau-b-needs-fy", pointer("materials", mat),
                             f"computed tau_b needs Fy of material {mat!r} (Pns = Fy*A)",
                             "Add `Fy` to the material, or set analysis.stiffnessReduction.tauB to \"none\" or a number.")

    # --------------------------------------------------------------- design
    def design(self):
        if not isinstance(self.d.get("design", {}), dict):
            return
        res = designparams.resolve(self.d, self.geoms, self.add, tolerant=not self.strict)
        if "design" not in self.d:
            return                      # analysis-only file: no design warnings
        mats = self.table("materials")
        missing: dict[str, list] = {}
        for dmid, info in res.items():
            vals = info["values"]
            lack = [k for k in ("Fy", "Fu") if vals[k].value is None]
            mat = self.table("members").get(info["members"][0], {}).get("material")
            if lack and isinstance(mats.get(mat), dict):
                missing.setdefault(mat, [])
                missing[mat] += [k for k in lack if k not in missing[mat]]
            for k in ("E", "G"):
                r = vals[k]
                ref = designparams.material_value(self.d, info["members"][0], k)
                if not r.source.startswith("material:") and ref is not None and isinstance(r.value, (int, float)) \
                        and abs(r.value - ref) > 1e-9 * abs(ref):
                    self.add("design.stiffness-differs", r.path,
                             f"design {k} = {r.value:g} for design member {dmid!r} differs from {k} = {ref:g} of "
                             "its material, which the analysis uses",
                             "Design variables never change the analysis; change the material if the analysis "
                             "should use this value too.", WARNING)
        for mat, lack in missing.items():
            self.add("design.missing-strength", pointer("materials", mat),
                     f"material {mat!r} has no {' and '.join(sorted(lack, reverse=True))}; members of this "
                     "material cannot be designed",
                     "Add Fy and Fu to the material (or set the design variables Fy and Fu).", WARNING)
        an = self.table("analysis")
        second = any(isinstance(c, dict) and c.get("analysisType", an.get("type", "firstOrder")) == "secondOrder"
                     for c in self.table("combinations").values())
        methods = {info["values"]["method"].value for info in res.values()}
        if second and "ASD" in methods and an.get("alpha", 1.0) == 1.0:
            self.add("design.asd-alpha", "/analysis/alpha" if "alpha" in an else "/analysis",
                     "design method ASD with second-order analysis but analysis.alpha = 1.0",
                     "AISC 360-22 C2.1(d): for ASD run the second-order analysis at 1.6 times the loads; set "
                     "analysis.alpha to 1.6. (Design variables never change the analysis by themselves.)", WARNING)


def validate_semantics(data, strict: bool = True) -> list[Problem]:
    """Layer 2. ``strict=False`` when the schema layer found errors."""
    if not isinstance(data, dict):
        return []
    return _Semantics(data, strict).run()


# ================================================================================ driver
def _version_problems(data) -> list[Problem]:
    v = data.get("formatVersion") if isinstance(data, dict) else None
    if not isinstance(v, str) or not re.fullmatch(r"\d+\.\d+\.\d+", v):
        return []
    major, minor, _ = (int(x) for x in v.split("."))
    if major != D.SUPPORTED_MAJOR or minor > D.SUPPORTED_MINOR:
        return [Problem("/formatVersion", f"format version {v} is not supported by this engine "
                        f"(it reads {D.SUPPORTED_MAJOR}.0.0 to {D.SUPPORTED_MAJOR}.{D.SUPPORTED_MINOR}.x)",
                        "Convert the file to a supported version (docs/FORMAT_VERSIONING.md) or upgrade benchfea.",
                        ERROR, "version.unsupported")]
    return []


def validate(source) -> list[Problem]:
    """Validate a document (path, dict or Document). Returns ALL problems,
    errors and warnings, each with path, line/column (when loaded from text),
    message and hint.  Raises ModelInvalid only for unreadable/unparsable input.
    """
    doc: Document = load(source)
    problems = list(doc.parse_problems)
    version = _version_problems(doc.data)
    problems += version
    if not version:
        schema = validate_schema(doc.data)
        problems += schema
        problems += validate_semantics(doc.data, strict=not schema and not doc.parse_problems)
    located = []
    seen = set()
    for p in problems:
        if p.line is None:
            line, col = doc.locate(p.path)
            p = Problem(p.path, p.message, p.hint, p.severity, p.code, line, col)
        key = (p.path, p.code, p.message)
        if key not in seen:
            seen.add(key)
            located.append(p)
    located.sort(key=lambda p: (p.severity != ERROR, p.line if p.line is not None else 10**9, p.col or 0, p.path))
    return located


def check(source) -> list[Problem]:
    """Validate and raise :class:`ModelInvalid` (carrying every problem) if
    there is any error. Returns the warnings otherwise."""
    doc = load(source)
    problems = validate(doc)
    if any(p.severity == ERROR for p in problems):
        raise ModelInvalid(problems, doc.source)
    return problems
