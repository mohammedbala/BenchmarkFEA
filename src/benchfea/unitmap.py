"""Schema-driven unit handling: find the dimension of every number in a
document from the ``x-unit`` annotations of the JSON Schema, and convert whole
documents between unit systems.

The schema is the single source of truth for units; no number is ever scaled
by hand-written per-field code.
"""
from __future__ import annotations

import copy
from typing import Any

from .problems import pointer
from .schema import load_schema, validator
from .units import POSITION, UnitSystem, convert

_NO_SCALE = {"n/a", "dimensionless", "degree", "rotation"}


def unit_map(data: Any) -> dict[str, str]:
    """Map JSON pointer -> dimension string for every number in ``data`` that
    carries a dimension (numbers that are dimensionless are omitted)."""
    root = load_schema()
    base = validator()
    out: dict[str, str] = {}

    def deref(s):
        while isinstance(s, dict) and "$ref" in s:
            target = root
            for part in s["$ref"].lstrip("#/").split("/"):
                target = target[part]
            s = {**target, **{k: v for k, v in s.items() if k != "$ref"}}
        return s

    def walk(inst, schema, ptr, owner):
        schema = deref(schema)
        if not isinstance(schema, dict):
            return
        unit = schema.get("x-unit")
        if unit and isinstance(inst, (int, float)) and not isinstance(inst, bool):
            if unit == POSITION:
                unit = "length" if isinstance(owner, dict) and owner.get("positionType") == "absolute" \
                    else "dimensionless"
            if unit in _NO_SCALE:
                out.pop(ptr, None)
            else:
                out[ptr] = unit
        if isinstance(inst, dict):
            props = schema.get("properties", {})
            addl = schema.get("additionalProperties")
            for k, v in inst.items():
                if k in props:
                    walk(v, props[k], ptr + pointer(k), inst)
                elif isinstance(addl, dict):
                    walk(v, addl, ptr + pointer(k), inst)
        elif isinstance(inst, list) and isinstance(schema.get("items"), dict):
            for k, v in enumerate(inst):
                walk(v, schema["items"], f"{ptr}/{k}", owner)
        for sub in schema.get("allOf", ()):
            walk(inst, sub, ptr, owner)
        for key in ("anyOf", "oneOf"):
            for sub in schema.get(key, ()):
                if base.evolve(schema={**sub, "$defs": root["$defs"]}).is_valid(inst):
                    walk(inst, sub, ptr, owner)
        if "if" in schema:
            ok = base.evolve(schema={**schema["if"], "$defs": root["$defs"]}).is_valid(inst)
            branch = schema.get("then" if ok else "else")
            if branch:
                walk(inst, branch, ptr, owner)

    walk(data, root, "", None)
    return out


def _set(data, ptr: str, fn):
    parts = [p.replace("~1", "/").replace("~0", "~") for p in ptr.split("/")[1:]]
    cur = data
    for p in parts[:-1]:
        cur = cur[int(p)] if isinstance(cur, list) else cur[p]
    last = int(parts[-1]) if isinstance(cur, list) else parts[-1]
    cur[last] = fn(cur[last])


def convert_document(data: dict, target: UnitSystem) -> dict:
    """Return a copy of a (schema-valid) document with every dimensional number
    converted to ``target`` and the ``units`` block updated."""
    src = UnitSystem.from_dict(data["units"])
    out = copy.deepcopy(data)
    if src == target:
        return out
    for ptr, dim in unit_map(data).items():
        _set(out, ptr, lambda v, dim=dim: convert(v, dim, src, target))
    out["units"] = target.to_dict()
    return out
