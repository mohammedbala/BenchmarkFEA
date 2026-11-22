"""The published JSON Schema: well-formed, complete, and in sync with its generator."""
import json

import pytest
from jsonschema import Draft202012Validator

from benchfea import defaults as D
from benchfea.schema import SCHEMA_PATH, load_schema
from benchfea.schema.build_schema import build
from benchfea.units import POSITION, parse_dimension


def test_schema_is_valid_draft_2020_12():
    schema = load_schema()
    assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    Draft202012Validator.check_schema(schema)


def test_committed_schema_matches_generator():
    assert json.loads(SCHEMA_PATH.read_text()) == build(), \
        "run `python -m benchfea.schema.build_schema` and commit the result"


def _property_schemas(schema, path="#"):
    """Yield (path, subschema) for every schema that describes a value."""
    if isinstance(schema, dict):
        for key in ("properties", "$defs"):
            for name, sub in schema.get(key, {}).items():
                yield f"{path}/{key}/{name}", sub
                yield from _property_schemas(sub, f"{path}/{key}/{name}")
        for key in ("additionalProperties", "items", "then", "if"):
            if isinstance(schema.get(key), dict):
                if key in ("additionalProperties", "items"):
                    yield f"{path}/{key}", schema[key]
                yield from _property_schemas(schema[key], f"{path}/{key}")
        for key in ("allOf", "anyOf", "oneOf"):
            for k, sub in enumerate(schema.get(key, ())):
                yield from _property_schemas(sub, f"{path}/{key}/{k}")


def test_every_field_has_description_and_unit():
    missing = []
    count = 0
    for path, sub in _property_schemas(load_schema()):
        if "/if/" in path or set(sub) <= {"const", "x-unit"}:
            continue                      # discriminator tests and unit overrides, not fields
        count += 1
        if "$ref" in sub and "description" not in sub:
            continue                      # described at the target
        if not sub.get("description"):
            missing.append(path + " (description)")
        if "x-unit" not in sub:
            missing.append(path + " (x-unit)")
    assert count > 150
    assert not missing, missing


def test_every_unit_annotation_is_parseable():
    for path, sub in _property_schemas(load_schema()):
        unit = sub.get("x-unit")
        if unit in (None, "n/a", POSITION):
            continue
        parse_dimension(unit)             # raises on an unknown dimension


def test_numbers_have_a_dimension_not_na():
    bad = [p for p, s in _property_schemas(load_schema())
           if s.get("type") in ("number", "integer") and s.get("x-unit") in (None, "n/a")]
    assert not bad, bad


def test_schema_defaults_come_from_defaults_module():
    s = load_schema()
    an = s["properties"]["analysis"]["properties"]
    for key, val in D.ANALYSIS.items():
        assert an[key]["default"] == val, key
    dv = s["$defs"]["designVariables"]["properties"]
    for key, (val, _reason) in D.DESIGN.items():
        assert dv[key]["default"] == val, key
    conv = s["properties"]["conventions"]["properties"]
    for key, val in D.CONVENTIONS.items():
        assert conv[key]["enum"] == [val] and conv[key]["default"] == val
    assert an["stiffnessReduction"]["properties"]["factor"]["default"] == 0.8
    assert s["$defs"]["combination"]["properties"]["notionalLoad"]["properties"]["coefficient"]["default"] == 0.002


def test_collections_are_id_keyed_objects():
    s = load_schema()
    for name in ("nodes", "supports", "materials", "sections", "members", "constraints", "loadCases",
                 "combinations"):
        coll = s["properties"][name]
        assert coll["type"] == "object" and "propertyNames" in coll, name


@pytest.mark.parametrize("ident,ok", [("N1", True), ("b-2.x_3", True), ("1", True), ("a b", False),
                                       ("a/b", False), ("-a", False), ("", False), ("a~b", False)])
def test_id_pattern(ident, ok):
    import re
    pattern = load_schema()["properties"]["nodes"]["propertyNames"]["pattern"]
    assert bool(re.search(pattern, ident)) is ok
