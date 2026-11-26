#!/usr/bin/env python3
"""Build docs/INPUT_REFERENCE.md from the JSON Schema of the input format.

    python3 tools/build_input_reference.py            # writes docs/INPUT_REFERENCE.md
    python3 tools/build_input_reference.py --check    # exit 1 if the file on disk differs or a field is incomplete

Everything in the document comes from ``src/benchfea/schema/benchfea-input-1.0.schema.json`` (written
by ``python -m benchfea.schema.build_schema``), from ``benchfea.defaults`` (the reasons of the built-in
defaults) and from the files ``examples/*.json`` (the example of a field is its first occurrence
there).  Nothing about a field is typed in this tool.  A field is *incomplete*, and the tool fails,
when the schema gives it

* no ``description``,
* no ``x-unit`` (the dimension, or ``n/a`` for names, switches and containers),
* no statement of what happens when it is left out: it is neither required, nor has a ``default``,
  nor does its description hold a sentence that begins with ``Default:`` or ``Absent:``,
* no example: none in ``examples/*.json``, no ``examples`` keyword in the schema, no default and no
  list of allowed values.

The remedy is to complete the schema source (``src/benchfea/schema/build_schema.py``), never to edit
the generated file.  Every ``examples`` value of the schema is validated against its own sub-schema.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from jsonschema import Draft202012Validator  # noqa: E402

from benchfea import defaults as D  # noqa: E402
from benchfea.schema import load_schema  # noqa: E402
from benchfea.units import POSITION, UnitSystem  # noqa: E402

OUT = ROOT / "docs" / "INPUT_REFERENCE.md"
EXAMPLES = sorted((ROOT / "examples").glob("[0-9]*.json"))
SI, US = UnitSystem("kN", "m"), UnitSystem("kip", "in")
ID, IDX = "<id>", "<n>"
DESIGN_LEVELS = "/design/{defaults | groups/<id> | members/<id>}"
SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


# ------------------------------------------------------------------------------------- schema walk
class Field:
    def __init__(self, path, schema, required, conditions, note=""):
        self.path, self.schema, self.required = path, schema, required
        self.conditions = conditions            # ((path of the object, key, value), ...): variant of that object
        self.also = []                          # further variants with the identical field
        self.note = note

    @property
    def shown_path(self):
        return "/" + "/".join(self.path)


def walk(schema, root):
    """Every field of the document described by ``schema`` -> list of Field, in schema order."""
    out = []

    def deref(s):
        seen = dict(s)
        while "$ref" in seen:
            target = root["$defs"][seen.pop("$ref").rsplit("/", 1)[1]]
            seen = {**target, **seen}            # keywords next to the $ref (description) win
        return seen

    def children(s, path, cond):
        s = deref(s)
        req = set(s.get("required", ()))
        for name, sub in s.get("properties", {}).items():
            field(sub, path + (name,), name in req, cond)
        extra = s.get("additionalProperties")
        if isinstance(extra, dict):
            field(extra, path + (ID,), False, cond)
        if isinstance(s.get("items"), dict):
            item = deref(s["items"])
            if "examples" in s and "examples" not in item and isinstance(s["examples"][0], list):
                item = {**item, "examples": s["examples"][0][:1]}      # the first item of the list's example
            field(item, path + (IDX,), False, cond)
        for part in s.get("allOf", ()):
            if "if" in part:
                when = tuple((path, k, v["const"]) for k, v in part["if"]["properties"].items())
                then = part["then"]
                treq = set(then.get("required", ()))
                for name, sub in then.get("properties", {}).items():
                    if "description" not in sub:           # an annotation on a field of another variant
                        for f in out:
                            if f.path == path + (name,):
                                f.note += (f" Unit `{sub['x-unit']}` when "
                                           + " and ".join(f"`{k}` = `{v}`" for _, k, v in when) + ".")
                        continue
                    field(sub, path + (name,), name in treq, cond + when)
            else:
                children(part, path, cond)
        for alt in s.get("anyOf", ()):
            if alt.get("type") == "object":
                children(alt, path, cond)

    def field(s, path, required, cond):
        s = deref(s)
        for f in out:                              # the identical field of another variant: one row
            if f.path == path and f.schema == s and f.required == required and f.conditions[:-1] == cond[:-1] \
                    and f.conditions and cond and f.conditions[-1][:2] == cond[-1][:2]:
                f.also.append(cond[-1][2])
                return
        out.append(Field(path, s, required, cond))
        children(s, path, cond)

    children(schema, (), ())
    return out


# ---------------------------------------------------------------------------------- field columns
def js(v):
    return json.dumps(v, ensure_ascii=False)


def code(textv):
    return "`" + textv.replace("|", "\\|") + "`"


def bounds(s):
    b = []
    for key, sym in (("exclusiveMinimum", ">"), ("minimum", "≥"), ("exclusiveMaximum", "<"), ("maximum", "≤")):
        if key in s:
            b.append(f"{sym} {s[key]:g}")
    return ", ".join(b)


def kind(s):
    """Type and allowed values, in words."""
    if "enum" in s:
        if len(s["enum"]) == 1:
            return "always " + code(js(s["enum"][0]))
        return "one of " + ", ".join(code(js(v)) for v in s["enum"])
    if "const" in s:
        return code(js(s["const"]))
    if "anyOf" in s and not all(set(a) <= {"required"} for a in s["anyOf"]):
        return " or ".join(kind(a) for a in s["anyOf"])
    t = s.get("type")
    if t in ("number", "integer"):
        b = bounds(s)
        return t + (f", {b}" if b else "")
    if t == "string":
        extra = []
        if "pattern" in s:
            extra.append("an id" if s["pattern"].startswith("^[A-Za-z0-9_]") else f"pattern {code(s['pattern'])}")
        if s.get("minLength"):
            extra.append("not empty")
        return "text" + (" (" + ", ".join(extra) + ")" if extra else "")
    if t == "boolean":
        return "`true` or `false`"
    if t == "array":
        n0, n1 = s.get("minItems"), s.get("maxItems")
        size = (f"exactly {n0}" if n0 is not None and n0 == n1 else f"at least {n0}" if n0 else "")
        item = s.get("items", {})
        what = {"number": "numbers", "string": "ids" if "pattern" in item else "texts"}.get(item.get("type"), "items")
        return "list of " + (size + " " if size else "") + what + (", no repeats" if s.get("uniqueItems") else "")
    if t == "object":
        extra = ""
        if isinstance(s.get("additionalProperties"), dict):
            extra = " keyed by id" + (f", at least {s['minProperties']}" if s.get("minProperties") else "")
        elif list(s.get("properties", {})) and "required" in s and len(s["properties"]) == 1:
            return "object " + code("{" + js(s["required"][0]) + ": …}")
        return "object" + extra
    raise ValueError(f"no type: {s}")


def unit(s):
    u = s.get("x-unit")
    if u is None:
        return None
    if u == "n/a":
        return "–"
    if u == "dimensionless":
        return "number without unit"
    if u == "degree":
        return "degrees (never converted)"
    if u == "rotation":
        return "rad (never converted)"
    if u == POSITION:
        return "fraction of the length, or a length, as `positionType` says"
    pretty = re.sub(r"\^(\d+)", lambda m: m.group(1).translate(SUP), u)
    return f"{pretty}: {lab(SI, u)} · {lab(US, u)}"


def lab(system, dim):
    return re.sub(r"\^(\d+)", lambda m: m.group(1).translate(SUP), system.label(dim)).replace("*", "·")


_STATEMENT = re.compile(r"\b(Default|Absent): (.*?)(?:\.(?=\s|$)|$)")


def default(f):
    s = f.schema
    if f.path[-1] in (ID, IDX):
        return "–"                               # an item of a collection, not a field that can be left out
    if f.required:
        return "required"
    if "default" in s:
        return code(js(s["default"]))
    m = _STATEMENT.search(s.get("description", ""))
    if m:
        textv = m.group(2).strip()
        return ("if absent: " if m.group(1) == "Absent" else "") + textv
    return None


def short(v, limit=70):
    """Compact JSON of ``v``; nested containers are abbreviated when the text gets long."""
    textv = js(v)
    if len(textv) <= limit:
        return textv

    def abbreviate(x, depth):
        if isinstance(x, dict):
            return "{…}" if depth > 0 else "{" + ", ".join(f"{js(k)}: {abbreviate(w, depth + 1)}" for k, w in x.items()) + "}"
        if isinstance(x, list):
            return "[…]" if depth > 0 else "[" + ", ".join(abbreviate(w, depth + 1) for w in x) + "]"
        return js(x)

    textv = abbreviate(v, 0)
    return textv if len(textv) <= limit else textv[:limit - 2] + " …"


def instances(doc, f):
    """The values of field ``f`` in a document (variant conditions respected)."""
    found = []

    def go(node, depth, trail):
        for cpath, key, value in f.conditions:
            if len(cpath) == depth and (not isinstance(node, dict) or node.get(key, _variant_default(key)) != value):
                return
        if depth == len(f.path):
            found.append(node)
            return
        seg = f.path[depth]
        if seg == ID and isinstance(node, dict):
            for w in node.values():
                go(w, depth + 1, trail + (seg,))
        elif seg == IDX and isinstance(node, list):
            for w in node:
                go(w, depth + 1, trail + (seg,))
        elif isinstance(node, dict) and seg in node:
            go(node[seg], depth + 1, trail + (seg,))

    go(doc, 0, ())
    if not found and f.also:
        for value in f.also:
            g = Field(f.path, f.schema, f.required, f.conditions[:-1] + (f.conditions[-1][:2] + (value,),))
            found.extend(instances(doc, g))
    return found


def _variant_default(key):
    return {"kind": "force"}.get(key)


def example(f, docs):
    s = f.schema
    name = f.path[-1]
    lead = "" if name in (ID, IDX) else js(name) + ": "
    if "examples" in s:
        return code(lead + short(s["examples"][0])), "schema"
    alternatives = [f]
    if f.path[:2] == ("design", "defaults"):      # a design variable: the same field at the two other levels
        alternatives += [Field(("design", level, ID) + f.path[2:], s, f.required, ()) for level in ("members", "groups")]
    for fname, doc in docs:
        vals = [v for g in alternatives for v in instances(doc, g)]
        if vals:
            return code(lead + short(min(vals, key=lambda v: len(js(v))))), fname
    if "default" in s:
        return code(lead + short(s["default"])), "default"
    if "enum" in s:
        return code(lead + js(s["enum"][0])), "allowed value"
    return None, None


def condition(f):
    if not f.conditions:
        return ""
    parts = [f"`{k}` is `{v}`" for _, k, v in f.conditions]
    if f.also:
        parts[-1] += "".join(f" or `{v}`" for v in f.also)
    return "only when " + " and ".join(parts)


def cell(textv):
    return textv.replace("|", "\\|").replace("\n", " ")


# ----------------------------------------------------------------------------------------- document
INTRO = """# Input reference (format {version})

Generated by `tools/build_input_reference.py` from the JSON Schema of the input format
(`src/benchfea/schema/benchfea-input-1.0.schema.json`); do not edit by hand. `--check` tells whether
this file is up to date, and `tests/test_input_reference.py` runs it.

An input file is one JSON object (a YAML file with the same content is accepted). This page lists
every field: where it is in the file, what kind of value it takes, its unit, what happens when it is
left out, and an example taken from the files in [`examples/`](../examples/README.md). Sign and
axis conventions are in [CONVENTIONS.md](CONVENTIONS.md); what the engine does not do is in
[LIMITS.md](../LIMITS.md).

## How to read the tables

* **Field** is the place of the value in the file, written as a path: `/members/<id>/section` is the
  entry `section` of a member inside `members`. `<id>` stands for a name you choose (a node id, a
  member id, …) and `<n>` for the position of an item in a list. An id may hold letters, digits,
  `_`, `.` and `-`, must not begin with `.` or `-`, and has at most 64 characters.
* **Value** is the kind of value and the values allowed.
* **Unit** gives the dimension and, after the colon, the unit that results in a file with
  `"units": {{"force": "kN", "length": "m"}}` and in one with `"units": {{"force": "kip", "length": "in"}}`.
  Every dimensional number of a file is in the consistent units built from the file's force unit
  ({forces}) and length unit ({lengths}): a modulus of elasticity in a kN, m file is in kN/m², not in
  MPa. Rotations are always in radians and the member roll angle in degrees; these are never
  converted. Results are written in the units of the file.
* **If left out** says `required`, or the value used when the field is absent. A value in code
  type is the built-in default.
* **Example** is the first occurrence of the field in the example files (the file is named in
  brackets), or the default where no example file uses the field.
* Unknown fields are refused: a misspelt name is reported with its line, not ignored.

A complete small file to read next to this page: [`examples/01_simple_beam.json`](../examples/01_simple_beam.json),
explained in [`examples/01_simple_beam.md`](../examples/01_simple_beam.md).
"""

SECTION_TITLES = {
    "": "File header",
    "units": "Units",
    "conventions": "Conventions",
    "nodes": "Nodes",
    "supports": "Supports",
    "materials": "Materials",
    "sections": "Sections",
    "members": "Members",
    "constraints": "Constraints (rigid links and rigid diaphragms)",
    "loadCases": "Load cases and loads",
    "combinations": "Combinations",
    "analysis": "Analysis options",
    "design": "Design input",
}
HEADER_FIELDS = ("$schema", "formatVersion", "name", "description")


def build_document():
    schema = load_schema()
    fields = walk(schema, schema)
    docs = [(p.name, json.loads(p.read_text(encoding="utf-8"))) for p in EXAMPLES]
    problems = []

    # the design variables are the same object at three places: listed once
    dv_names = list(schema["$defs"]["designVariables"]["properties"])
    listed, dv_rows = [], []
    for f in fields:
        p = f.path
        if p[:1] == ("design",) and len(p) >= 2:
            if p[1] == "defaults" and len(p) >= 3:
                dv_rows.append(f)
                continue
            if p[1] in ("groups", "members") and len(p) >= 4 and p[3] in dv_names:
                continue
        listed.append(f)

    def row(f, shown=None):
        s = f.schema
        where = shown or f.shown_path
        desc = s.get("description")
        u, dflt = unit(s), default(f)
        ex, src = example(f, docs)
        for label, value in (("description", desc), ("x-unit", u), ("default statement", dflt), ("example", ex)):
            if value is None:
                problems.append(f"{where}: no {label}")
        if "examples" in s:
            validator = Draft202012Validator({**s, "$defs": schema["$defs"]})
            for v in s["examples"]:
                for e in validator.iter_errors(v):
                    problems.append(f"{where}: schema example {js(v)} is not valid: {e.message}")
        desc = (desc or "") + f.note
        cond = condition(f)
        if cond:
            desc = f"*{cond[0].upper() + cond[1:]}.* " + desc
        if ex and src not in ("default", "allowed value", "schema"):
            ex += f" ({src[:2]})"
        elif ex and src == "default":
            ex += " (default)"
        return (f"| {code(where)} | {cell(kind(s))} | {cell(u or '?')} | {cell(dflt or '?')} | {cell(desc)} | "
                f"{ex or '?'} |")

    head = ("| Field | Value | Unit | If left out | Description | Example |\n"
            "|---|---|---|---|---|---|")
    forces = ", ".join(schema["properties"]["units"]["properties"]["force"]["enum"])
    lengths = ", ".join(schema["properties"]["units"]["properties"]["length"]["enum"])
    out = [INTRO.format(version=D.FORMAT_VERSION, forces=forces, lengths=lengths)]
    out.append("Example files, by the two-digit number shown in brackets: "
               + ", ".join(f"({p.name[:2]}) [`{p.name}`](../examples/{p.name})" for p in EXAMPLES) + ".\n")

    first = docs[0][1]
    for key, title in SECTION_TITLES.items():
        rows = [f for f in listed if (f.path[0] in HEADER_FIELDS and len(f.path) == 1 if key == "" else f.path[0] == key)]
        out.append(f"## {title}\n")
        if key:
            top = next(f for f in rows if len(f.path) == 1)
            out.append(top.schema["description"] + "\n")
            if key in first:
                snippet = json.dumps({key: first[key]}, indent=2, ensure_ascii=False)
                lines = [ln[2:] for ln in snippet.splitlines()[1:-1]]      # without the outer braces
                if len(lines) > 26:
                    lines = lines[:25] + ["  …"]
                out.append(f"From `{docs[0][0]}`:\n\n```json\n" + "\n".join(lines) + "\n```\n")
        out.append(head)
        out.extend(row(f) for f in rows)
        out.append("")
        if key == "design":
            out.append(design_variables(schema, dv_rows, row, head))

    out.append(defaults_section())
    return "\n".join(out).rstrip() + "\n", problems


def design_variables(schema, dv_rows, row, head):
    dv = schema["$defs"]["designVariables"]
    out = ["### Design variables\n", dv["description"] + "\n",
           "The same variables may be given at three levels. For each variable of a design member the value is "
           "taken from the first level that gives it:\n\n"
           "1. the design member itself, `/design/members/<id>/…`;\n"
           "2. the group or groups that list the design member, `/design/groups/<id>/…` (when several groups give "
           "the same variable, the one with the highest `priority` is used; equal priority with different values "
           "is an error);\n"
           "3. `/design/defaults/…`;\n"
           "4. the built-in default of the table below.\n\n"
           "In short: member > group > defaults. The design results name, for every variable, the level its "
           "value came from.\n",
           f"The paths below are written `{DESIGN_LEVELS}/…`: the field may stand at any of the three places.\n",
           head]
    for f in dv_rows:
        out.append(row(f, DESIGN_LEVELS + "/" + "/".join(f.path[2:])))
    out.append("")
    out.append("#### Built-in defaults of the design variables\n")
    out.append("| Variable | Built-in default | Why |\n|---|---|---|")
    for name, sub in dv["properties"].items():
        if name in D.DESIGN:
            value, reason = D.DESIGN[name]
            out.append(f"| `{name}` | {code(js(value))} | {cell(reason)} |")
        else:
            m = _STATEMENT.search(sub.get("description", ""))
            if m and m.group(1) == "Default":
                out.append(f"| `{name}` | from the model: {cell(m.group(2).strip())} | derived, not a fixed value |")
            elif m:
                out.append(f"| `{name}` | none: {cell(m.group(2).strip())} | – |")
    out.append("")
    return "\n".join(out)


def defaults_section():
    out = ["## Built-in defaults of the analysis options, and their reasons\n",
           "Every default the engine used for a file is listed in the results under `defaultsUsed`, and "
           "`benchfea validate FILE --defaults` prints them.\n",
           "| Field | Built-in default | Why |\n|---|---|---|"]
    for name, value in D.ANALYSIS.items():
        out.append(f"| `/analysis/{name}` | {code(js(value))} | {cell(D.ANALYSIS_REASONS[name])} |")
    for name, value in D.STIFFNESS_REDUCTION.items():
        out.append(f"| `/analysis/stiffnessReduction/{name}` | {code(js(value))} | "
                   f"{cell(D.STIFFNESS_REDUCTION_REASONS[name])} |")
    out.append(f"| `/analysis/buckling/modes` | {code(js(D.BUCKLING['modes']))} | engine default |")
    out.append(f"| `/combinations/<id>/notionalLoad/coefficient` | {code(js(D.NOTIONAL['coefficient']))} | "
               f"{cell(D.NOTIONAL_REASON)} |")
    for name, value in D.COMBINATION.items():
        out.append(f"| `/combinations/<id>/{name}` | {code(js(value))} | engine default |")
    out.append("")
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if docs/INPUT_REFERENCE.md differs or a field is incomplete")
    args = ap.parse_args(argv)
    textv, problems = build_document()
    if problems:
        print("the schema is incomplete (fix src/benchfea/schema/build_schema.py, not the generated file):")
        for p in problems:
            print("  " + p)
        return 1
    if args.check:
        if not OUT.exists() or OUT.read_text(encoding="utf-8") != textv:
            print(f"{OUT.relative_to(ROOT)} is not up to date: run python3 tools/build_input_reference.py")
            return 1
        print(f"{OUT.relative_to(ROOT)} is up to date")
        return 0
    OUT.write_text(textv, encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
