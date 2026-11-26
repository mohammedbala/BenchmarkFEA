"""docs/INPUT_REFERENCE.md is what tools/build_input_reference.py gives from the schema, and the schema
states for every field a description, a unit, what happens when the field is left out, and an example.

A failure means "complete src/benchfea/schema/build_schema.py and run the tool again", never "edit the
generated file".
"""
import copy
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _tool():
    sys.path.insert(0, str(ROOT / "tools"))
    try:
        import build_input_reference as B
    finally:
        sys.path.pop(0)
    return B


def test_committed_input_reference_is_up_to_date():
    p = subprocess.run([sys.executable, str(ROOT / "tools" / "build_input_reference.py"), "--check"],
                       capture_output=True, text=True, cwd=ROOT)
    assert p.returncode == 0, p.stdout + p.stderr


def test_every_field_of_the_schema_is_complete_and_listed():
    B = _tool()
    text, problems = B.build_document()
    assert problems == []
    schema = B.load_schema()
    fields = B.walk(schema, schema)
    assert len(fields) > 300
    for f in fields:
        p = f.path
        if p[:1] == ("design",) and len(p) >= 4 and p[1] in ("groups", "members") \
                and p[3] in schema["$defs"]["designVariables"]["properties"]:
            continue                                   # listed once, under the three-level path
        shown = (B.DESIGN_LEVELS + "/" + "/".join(p[2:]) if p[:2] == ("design", "defaults") and len(p) > 2
                 else f.shown_path)
        assert "| `" + shown.replace("|", "\\|") + "` |" in text, shown


def test_every_design_variable_is_in_the_defaults_table_or_derived():
    B = _tool()
    text, _ = B.build_document()
    table = text.split("#### Built-in defaults of the design variables", 1)[1].split("\n## ", 1)[0]
    for name, (value, reason) in B.D.DESIGN.items():
        assert f"| `{name}` | `{B.js(value)}` |" in table, name
        assert reason.replace("|", "\\|") in table, name
    assert "member > group > defaults" in text


def test_an_incomplete_field_is_reported(monkeypatch):
    B = _tool()
    schema = copy.deepcopy(B.load_schema())       # load_schema() is cached: never change its result
    del schema["properties"]["nodes"]["additionalProperties"]["properties"]["x"]["x-unit"]
    schema["properties"]["name"]["description"] = "Short model name."
    monkeypatch.setattr(B, "load_schema", lambda: schema)
    _, problems = B.build_document()
    assert "/nodes/<id>/x: no x-unit" in problems
    assert "/name: no default statement" in problems
