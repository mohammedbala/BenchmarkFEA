"""Loading documents: JSON/YAML equivalence, positions, canonical dump."""
import json

import pytest

from benchfea import io
from benchfea.model import build_model
from benchfea.problems import ModelInvalid
from benchfea.validate import validate

from helpers import assert_same


def test_json_and_yaml_are_equivalent(fixtures):
    j = io.load(fixtures / "valid" / "beam.json")
    y = io.load(fixtures / "valid" / "beam.yaml")
    assert j.format == "json" and y.format == "yaml"
    assert j.data == y.data
    assert_same(build_model(j), build_model(y), ignore=("source", "warnings"))


def test_yaml_reads_json_style_exponent_numbers():
    doc = io.loads("a: 1e5\nb: 2.5E-3\nc: -1e+2\nd: 1.0e5\ne: 12\nf: abc1e5\n", "yaml")
    assert doc.data == {"a": 1e5, "b": 2.5e-3, "c": -100.0, "d": 1e5, "e": 12, "f": "abc1e5"}
    assert all(isinstance(doc.data[k], float) for k in "abcd") and isinstance(doc.data["e"], int)


def test_dump_round_trip_json_and_yaml(fixtures, tmp_path):
    data = io.load(fixtures / "valid" / "frame_kip_in.json").data
    for ext in ("json", "yaml"):
        path = tmp_path / f"out.{ext}"
        io.dump(data, path)
        assert io.load(path).data == data
    text = io.dumps(data)
    assert text.endswith("\n") and json.loads(text) == data
    assert list(json.loads(text)) == list(data)          # key order preserved


def test_json_positions_are_one_based_line_and_column():
    text = '{\n  "nodes": {\n    "A": {"x": 0,\n          "y": 1.5}\n  },\n  "list": [10,\n   20]\n}\n'
    doc = io.loads(text, "json")
    assert doc.positions[""] == (1, 1)
    assert doc.positions["/nodes"] == (2, 12)
    assert doc.positions["/nodes/A"] == (3, 10)
    assert doc.positions["/nodes/A/x"] == (3, 16)
    assert doc.positions["/nodes/A/y"] == (4, 16)
    assert doc.positions["/list/0"] == (6, 12)
    assert doc.positions["/list/1"] == (7, 4)
    assert doc.locate("/nodes/A/missing/deeper") == (3, 10)   # nearest ancestor


def test_yaml_positions():
    doc = io.loads("nodes:\n  A: {x: 0, y: 1}\nlist:\n  - 5\n  - 6\n", "yaml")
    assert doc.positions["/nodes/A"] == (2, 6)
    assert doc.positions["/nodes/A/y"] == (2, 16)
    assert doc.positions["/list/1"] == (5, 5)


def test_dict_source_has_no_positions_and_is_copied():
    data = {"a": [1, 2]}
    doc = io.load(data)
    assert doc.format == "dict" and doc.locate("/a/0") == (None, None)
    doc.data["a"].append(3)
    assert data == {"a": [1, 2]}


@pytest.mark.parametrize("text,line,col", [
    ('{"a": 1,}', 1, 9),
    ('{"a": NaN}', 1, 7),
    ('{"a": 1} x', 1, 10),
    ('{\n "a": [1, 2\n}', 3, 1),
    ("{'a': 1}", 1, 2),
    ('{"a": 01}', 1, 8),
])
def test_json_syntax_errors_are_located(text, line, col):
    with pytest.raises(ModelInvalid) as e:
        io.loads(text, "json")
    (p,) = e.value.problems
    assert (p.code, p.line, p.col) == ("parse.json", line, col) and p.hint


def test_yaml_syntax_error_is_located():
    with pytest.raises(ModelInvalid) as e:
        io.loads("a: [1, 2\nb: 3\n", "yaml")
    (p,) = e.value.problems
    assert p.code == "parse.yaml" and p.line is not None


def test_json_strings_and_escapes():
    doc = io.loads(r'{"a": "x\"yé\n", "b": [true, false, null, -0.5e-2, 7]}', "json")
    assert doc.data == {"a": 'x"yé\n', "b": [True, False, None, -0.005, 7]}


def test_unsupported_extension_and_missing_file(tmp_path):
    with pytest.raises(ModelInvalid) as e:
        io.load(tmp_path / "model.txt")
    assert e.value.problems[0].code == "parse.extension"
    with pytest.raises(ModelInvalid) as e:
        io.load(tmp_path / "nope.json")
    assert e.value.problems[0].code == "parse.io"


def test_yaml_non_json_values_are_refused():
    doc = io.loads("a: 2024-01-01\nb: .nan\nc: {<<: {x: 1}}\n", "yaml")
    codes = sorted((p.code, p.path) for p in doc.parse_problems)
    assert codes == [("parse.value", "/a"), ("parse.value", "/b"), ("parse.yaml", "/c")]


def test_valid_fixtures_have_no_errors(fixtures):
    for path in sorted((fixtures / "valid").iterdir()):
        errors = [p for p in validate(path) if p.severity == "error"]
        assert not errors, (path.name, errors)
