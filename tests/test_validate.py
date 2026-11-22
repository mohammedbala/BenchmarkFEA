"""Validation: every invalid fixture yields ALL expected problems, each located."""
import copy

import pytest

from benchfea import io
from benchfea.model import build_model
from benchfea.problems import ModelInvalid, Problem
from benchfea.validate import check, validate, validate_schema, validate_semantics

# fixture -> every expected ERROR as (code, JSON pointer, line in the file)
EXPECTED = {
    "schema_errors.json": [
        ("schema.pattern", "/formatVersion", 2),
        ("schema.required", "/units", 3),                      # force
        ("schema.enum", "/units/length", 3),
        ("schema.enum", "/conventions/internalForces", 4),     # unsupported convention
        ("id.invalid", "/nodes/n 1", 6),
        ("schema.type", "/nodes/n2/x", 7),
        ("schema.required", "/nodes/n3", 8),                   # y
        ("schema.invalid", "/supports/n2/dofs/ux", 11),
        ("schema.range", "/supports/n2/dofs/rz/spring", 11),
        ("schema.required", "/supports/n3", 12),               # type or dofs
        ("schema.required", "/materials/m", 15),               # G or nu
        ("schema.range", "/materials/m/E", 15),
        ("schema.unknown-property", "/materials/m2/density", 16),
        ("schema.required", "/sections/s1", 19),               # tw
        ("schema.enum", "/sections/s2/type", 20),
        ("schema.unknown-property", "/sections/s3/bf", 21),
        ("schema.required", "/sections/s4", 22),               # type
        ("schema.required", "/members/b1", 25),                # section
        ("schema.unknown-property", "/members/b1/sction", 25),
        ("schema.enum", "/members/b1/type", 25),
        ("schema.invalid", "/members/b2/releases/i/rz", 26),
        ("schema.type", "/loadCases/D/loads/0/Fz", 30),
        ("schema.required", "/loadCases/D/loads/1", 31),       # position
        ("schema.enum", "/loadCases/D/loads/1/direction", 31),
        ("schema.enum", "/loadCases/D/loads/2/type", 32),
        ("schema.minProperties", "/combinations", 35),
        ("schema.range", "/analysis/resultStations", 36),
        ("schema.enum", "/analysis/type", 36),
        ("schema.unknown-property", "/analysis/pDelta", 36),
        ("schema.range", "/design/defaults/Kx", 37),
        ("schema.invalid", "/design/defaults/Cb", 37),
        ("schema.required", "/design/defaults/deflectionLimits/0", 37),
    ],
    "references.json": [
        ("ref.unknown-node", "/supports/Q", 11),
        ("ref.unknown-node", "/members/M2/j", 17),
        ("ref.unknown-section", "/members/M2/section", 17),
        ("ref.unknown-material", "/members/M2/material", 17),
        ("ref.unknown-node", "/constraints/link/master", 20),
        ("ref.unknown-node", "/constraints/link/slaves/1", 20),
        ("ref.unknown-node", "/loadCases/D/loads/0/node", 24),
        ("ref.unknown-member", "/loadCases/D/loads/1/member", 25),
        ("ref.unknown-member", "/loadCases/D/loads/2/members/1", 26),
        ("ref.unknown-node", "/loadCases/D/loads/3/node", 27),
        ("ref.unknown-case", "/combinations/C1/factors/LL", 31),
        ("ref.unknown-combination", "/analysis/buckling/combinations/1", 34),
        ("ref.unknown-member", "/analysis/stiffnessReduction/members/1", 35),
        ("ref.unknown-combination", "/design/defaults/deflectionLimits/0/combinations/0", 38),
        ("ref.unknown-design-member", "/design/groups/g/members/0", 39),   # M1 is part of DM1
        ("ref.unknown-design-member", "/design/groups/g/members/1", 39),
        ("ref.unknown-member", "/design/members/DM1/members/1", 41),
        ("ref.unknown-member", "/design/members/DM2", 42),
    ],
    "geometry.json": [
        ("section.geometry", "/sections/I1/tf", 14),
        ("section.geometry", "/sections/I1/tw", 14),
        ("section.geometry", "/sections/I1/r", 14),
        ("section.geometry", "/sections/P1/t", 15),
        ("section.geometry", "/sections/X1/t", 16),
        ("member.zero-length", "/members/zero", 20),
        ("member.same-node", "/members/same/j", 21),
        ("member.offset-zero-length", "/members/offs/offsets", 23),
        ("member.orientation-collinear", "/members/refpar/orientation/referenceVector", 25),
        ("member.orientation-collinear", "/members/refzero/orientation/referenceVector", 27),
    ],
    "releases.json": [
        ("member.release-mechanism", "/members/axial/releases", 14),
        ("member.release-mechanism", "/members/spin/releases", 16),
        ("member.release-mechanism", "/members/shear/releases", 18),
        ("member.release-mechanism", "/members/swing/releases", 20),
        ("member.release-on-truss", "/members/truss/releases", 22),
    ],
    "constraints.json": [
        ("constraint.slave-restrained", "/constraints/d1/slaves/1", 26),
        ("constraint.master-is-slave", "/constraints/d1/slaves/3", 26),
        ("constraint.slave-two-masters", "/constraints/d2/slaves/0", 27),
        ("constraint.chained", "/constraints/d3/slaves/0", 28),
    ],
    "loads.json": [
        ("load.position-range", "/loadCases/D/loads/0/position", 20),
        ("load.position-range", "/loadCases/D/loads/1/position", 21),
        ("load.position-range", "/loadCases/D/loads/2/end", 23),
        ("load.projected-local", "/loadCases/D/loads/3/projected", 24),
        ("load.zero-direction", "/loadCases/D/loads/4/direction", 25),
        ("load.displacement-dof-not-fixed", "/loadCases/D/loads/5/ry", 26),
        ("load.displacement-dof-not-fixed", "/loadCases/D/loads/6/uz", 27),
        ("combination.notional-direction", "/combinations/C1/notionalLoad/direction", 33),
        ("combination.notional-direction", "/combinations/C2/notionalLoad/direction", 34),
    ],
    "design.json": [
        ("design.length-and-bracing", "/design/groups/g1/Lb", 37),
        ("design.group-conflict", "/design/groups/g2/Kx", 38),
        ("design.not-collinear", "/design/members/bent/members/1", 42),
        ("design.not-prismatic", "/design/members/stepped/members/1", 43),
        ("design.not-connected", "/design/members/reversed/members/1", 44),
        ("design.id-collision", "/design/members/M7", 45),
        ("design.length-and-bracing", "/design/members/M1/Ly", 46),
        ("design.bracing-range", "/design/members/M1/bracing/torsion/0", 46),
        ("design.bracing-range", "/design/members/M1/bracing/torsion/2", 46),
    ],
    "duplicate_keys.json": [
        ("id.duplicate", "/nodes/A", 8),
        ("id.duplicate", "/members/M1", 15),
    ],
    "yaml_keys.yaml": [
        ("parse.key-type", "/nodes", 5),
        ("parse.key-type", "/nodes", 6),
        ("id.duplicate", "/nodes/top", 8),
        ("ref.unknown-node", "/supports/1", 10),               # consequences of the refused keys
        ("ref.unknown-node", "/members/M1/i", 16),
        ("ref.unknown-node", "/members/M1/j", 16),
    ],
    "version.json": [
        ("version.unsupported", "/formatVersion", 2),
    ],
    "standard_shape.json": [
        ("section.shapes-db", "/sections/W/shape", 12),
    ],
    "many_errors.json": [        # one file with schema, reference, geometry, release, constraint, load, analysis and design errors
        ("schema.enum", "/conventions/handedness", 4),
        ("schema.type", "/nodes/E/x", 10),
        ("ref.unknown-node", "/supports/Z", 14),
        ("schema.range", "/supports/C/dofs/rx/spring", 15),
        ("schema.range", "/materials/t/nu", 19),
        ("section.geometry", "/sections/i/tf", 22),
        ("section.geometry", "/sections/i/r", 22),
        ("schema.range", "/sections/i/overrides/A", 22),
        ("schema.unknown-property", "/sections/i/overrides/Q", 22),
        ("section.geometry", "/sections/p/t", 23),
        ("schema.unknown-property", "/sections/p/H", 23),
        ("schema.required", "/sections/q", 24),
        ("member.release-mechanism", "/members/m1/releases", 27),
        ("ref.unknown-section", "/members/m2/section", 28),
        ("ref.unknown-material", "/members/m2/material", 28),
        ("member.release-on-truss", "/members/m2/releases", 28),
        ("schema.minItems", "/members/m2/offsets/i", 29),
        ("schema.enum", "/members/m2/offsets/axes", 29),
        ("member.zero-length", "/members/m3", 30),
        ("schema.type", "/members/m3/orientation/roll", 30),
        ("member.offset-zero-length", "/members/m5/offsets", 32),
        ("schema.uniqueItems", "/constraints/c1/slaves", 35),
        ("constraint.master-is-slave", "/constraints/c1/slaves/0", 35),
        ("constraint.chained", "/constraints/c1/slaves/1", 35),
        ("constraint.slave-restrained", "/constraints/c1/slaves/1", 35),
        ("schema.unknown-property", "/constraints/c1/normal", 35),
        ("schema.enum", "/constraints/c2/normal", 36),
        ("schema.unknown-property", "/loadCases/D/loads/0/Fw", 40),
        ("load.position-range", "/loadCases/D/loads/1/end", 41),
        ("load.position-range", "/loadCases/D/loads/2/end", 42),
        ("load.projected-local", "/loadCases/D/loads/2/projected", 42),
        ("schema.range", "/loadCases/D/loads/3/position", 43),
        ("ref.unknown-member", "/loadCases/D/loads/4/member", 44),
        ("schema.enum", "/loadCases/D/loads/4/kind", 44),
        ("load.displacement-dof-not-fixed", "/loadCases/D/loads/5/ux", 45),
        ("schema.minItems", "/loadCases/D/loads/6/direction", 46),
        ("load.zero-direction", "/loadCases/D/loads/6/direction", 46),
        ("schema.uniqueItems", "/loadCases/D/loads/6/members", 46),
        ("ref.unknown-member", "/loadCases/D/loads/6/members/2", 46),
        ("ref.unknown-case", "/combinations/c/factors/Q", 51),
        ("schema.enum", "/combinations/c/analysisType", 51),
        ("schema.required", "/combinations/c/notionalLoad", 51),
        ("schema.range", "/combinations/c/notionalLoad/coefficient", 51),
        ("schema.minProperties", "/combinations/d/factors", 52),
        ("schema.range", "/analysis/maxIterations", 54),
        ("ref.unknown-combination", "/analysis/buckling/combinations/0", 54),
        ("schema.range", "/analysis/buckling/modes", 54),
        ("schema.range", "/analysis/stiffnessReduction/factor", 55),
        ("schema.invalid", "/analysis/stiffnessReduction/tauB", 55),
        ("ref.unknown-member", "/analysis/stiffnessReduction/members/0", 55),
        ("schema.enum", "/design/defaults/method", 57),
        ("schema.invalid", "/design/defaults/Cb", 57),
        ("schema.range", "/design/defaults/U", 57),
        ("schema.required", "/design/defaults/deflectionLimits/0", 57),
        ("schema.enum", "/design/defaults/deflectionLimits/0/direction", 57),
        ("ref.unknown-design-member", "/design/groups/g/members/0", 58),
        ("ref.unknown-design-member", "/design/groups/g/members/1", 58),
        ("schema.range", "/design/groups/g/Lx", 58),
        ("design.length-and-bracing", "/design/groups/g/Lx", 58),
        ("schema.uniqueItems", "/design/groups/g/bracing/major", 58),
        ("design.bracing-range", "/design/groups/g/bracing/major/2", 58),
        ("schema.range", "/design/members/m1/Kx", 59),
        ("schema.minItems", "/design/members/m1/members", 59),
    ],
    "many_errors.yaml": [        # YAML: block and flow style, YAML-only pitfalls (octal, hex, sexagesimal, yes/no, dates, aliases, merge keys)
        ("schema.type", "/formatVersion", 1),
        ("schema.enum", "/conventions/verticalAxis", 3),
        ("parse.value", "/nodes/n2/x", 6),
        ("parse.value", "/nodes/n2/y", 6),
        ("parse.value", "/nodes/n2/z", 6),
        ("parse.key-type", "/nodes", 7),
        ("parse.value", "/nodes/n3/x", 8),
        ("parse.value", "/nodes/n4/y", 9),
        ("parse.key-type", "/nodes", 12),
        ("parse.value", "/nodes/n7/y", 13),
        ("schema.type", "/nodes/n7/y", 13),
        ("parse.key-type", "/nodes", 14),
        ("schema.type", "/supports/n1", 17),
        ("schema.invalid", "/supports/n2/dofs/uy", 18),
        ("id.duplicate", "/materials/s", 21),
        ("schema.enum", "/members/m1/type", 25),
        ("schema.type", "/members/m2/j", 26),
        ("member.release-mechanism", "/members/m3/releases", 33),
        ("ref.unknown-node", "/constraints/c1/slaves/1", 46),
        ("schema.enum", "/constraints/c1/normal", 47),
        ("load.position-range", "/loadCases/D/loads/1/end", 56),
        ("parse.yaml", "/loadCases/D/loads/2", 62),
        ("load.position-range", "/loadCases/D/loads/3/position", 67),
        ("ref.unknown-case", "/combinations/c1/factors/L", 70),
        ("schema.range", "/design/defaults/E", 73),
        ("ref.unknown-design-member", "/design/groups/g1/members/0", 76),
        ("ref.unknown-design-member", "/design/groups/g1/members/1", 76),
        ("design.length-and-bracing", "/design/groups/g1/Lb", 77),
        ("design.bracing-range", "/design/groups/g1/bracing/minor/1", 82),
        ("design.group-conflict", "/design/groups/g3/Kx", 88),
        ("schema.type", "/design/members/m3", 90),
        ("design.not-connected", "/design/members/chain/members/1", 94),
    ],
}


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_invalid_fixture_reports_every_problem_with_location(fixtures, name):
    problems = validate(fixtures / "invalid" / name)
    errors = sorted((p.code, p.path, p.line) for p in problems if p.severity == "error")
    assert errors == sorted(EXPECTED[name])
    for p in problems:
        assert p.line is not None and p.col is not None and p.col >= 1, p
        assert p.message and p.hint, p                      # every problem explains how to fix it
        assert p.path == "" or p.path.startswith("/")


def test_every_invalid_fixture_is_covered(fixtures):
    names = {p.name for p in (fixtures / "invalid").iterdir()}
    assert names == set(EXPECTED) | {"syntax.json"}


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_build_model_refuses_with_all_problems(fixtures, name):
    with pytest.raises(ModelInvalid) as e:
        build_model(fixtures / "invalid" / name)
    got = sorted((p.code, p.path, p.line) for p in e.value.errors)
    assert got == sorted(EXPECTED[name])
    assert name in str(e.value) and "hint:" in str(e.value)


def test_columns_point_at_the_value(fixtures):
    path = fixtures / "invalid" / "references.json"
    lines = path.read_text().splitlines()
    by_path = {p.path: p for p in validate(path)}
    p = by_path["/members/M2/section"]
    assert lines[p.line - 1][p.col - 1:].startswith('"R2"')
    p = by_path["/combinations/C1/factors/LL"]
    assert lines[p.line - 1][p.col - 1:].startswith("1.6")
    p = by_path["/constraints/link/slaves/1"]
    assert lines[p.line - 1][p.col - 1:].startswith('"X2"')


def test_syntax_error_is_one_located_problem(fixtures):
    with pytest.raises(ModelInvalid) as e:
        validate(fixtures / "invalid" / "syntax.json")
    (p,) = e.value.problems
    assert (p.code, p.line, p.col) == ("parse.json", 4, 34)


def test_hints_are_specific(fixtures):
    by_path = {(p.code, p.path): p for p in validate(fixtures / "invalid" / "schema_errors.json")}
    assert "Did you mean 'section'" in by_path[("schema.unknown-property", "/members/b1/sction")].hint
    assert "'rectangle'" in by_path[("schema.enum", "/sections/s2/type")].hint
    assert "positive-face" in by_path[("schema.enum", "/conventions/internalForces")].hint
    assert "force" in by_path[("schema.required", "/units")].message
    assert "tw" in by_path[("schema.required", "/sections/s1")].message
    assert "`G` or `nu`" in by_path[("schema.required", "/materials/m")].message
    by_path = {(p.code, p.path): p for p in validate(fixtures / "invalid" / "standard_shape.json")}
    hint = by_path[("section.shapes-db", "/sections/W/shape")].hint
    assert "BENCHFEA_SHAPES_DB" in hint and "data/user" in hint and "not redistributable" in hint


def test_warnings_do_not_refuse(fixtures):
    problems = validate(fixtures / "valid" / "frame_kip_in.json")
    assert {p.severity for p in problems} == {"warning"}
    assert sorted((p.code, p.path) for p in problems) == [
        ("design.missing-strength", "/materials/other"),
        ("load.self-weight-zero", "/loadCases/D/loads/0"),
    ]
    assert check(fixtures / "valid" / "frame_kip_in.json") == problems
    model = build_model(fixtures / "valid" / "frame_kip_in.json")
    assert [p.code for p in model.warnings] == ["design.missing-strength", "load.self-weight-zero"]


def test_expected_warnings(fixtures):
    def warnings(name):
        return sorted((p.code, p.path) for p in validate(fixtures / "invalid" / name) if p.severity == "warning")

    assert ("node.coincident", "/nodes/A2") in warnings("geometry.json")
    assert ("constraint.diaphragm-out-of-plane", "/constraints/d3/slaves/0") in warnings("constraints.json")
    assert ("node.unused", "/nodes/top") in warnings("yaml_keys.yaml")
    assert ("support.no-restraint", "/supports/n3") in warnings("schema_errors.json")


def base():
    return {
        "formatVersion": "1.0.0", "units": {"force": "kN", "length": "m"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0}, "B": {"x": 4, "y": 0}, "C": {"x": 4, "y": 0, "z": 3}},
        "supports": {"A": {"type": "fixed"}},
        "materials": {"S": {"E": 2e8, "nu": 0.3, "Fy": 355e3, "Fu": 470e3}},
        "sections": {"R": {"type": "rectangle", "h": 0.3, "b": 0.2}},
        "members": {"M1": {"i": "A", "j": "B", "section": "R", "material": "S"},
                    "M2": {"i": "B", "j": "C", "section": "R", "material": "S"}},
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fz": -1}]}},
        "combinations": {"C1": {"factors": {"D": 1.0}}},
    }


def test_base_document_is_valid_and_dict_problems_have_no_line():
    assert validate(base()) == []
    d = base()
    d["members"]["M1"]["j"] = "nowhere"
    (p,) = validate(d)
    assert (p.code, p.path, p.line, p.col) == ("ref.unknown-node", "/members/M1/j", None, None)


def test_schema_and_semantic_problems_are_reported_together():
    d = base()
    d["materials"]["S"]["E"] = "stiff"                      # schema
    d["members"]["M2"]["section"] = "nope"                  # semantic
    d["combinations"]["C1"]["factors"]["X"] = 1             # semantic
    d["nodes"]["B"]["x"] = None                             # schema; M1/M2 geometry cannot be inspected
    got = sorted((p.code, p.path) for p in validate(d) if p.severity == "error")
    assert got == [("ref.unknown-case", "/combinations/C1/factors/X"),
                   ("ref.unknown-section", "/members/M2/section"),
                   ("schema.type", "/materials/S/E"), ("schema.type", "/nodes/B/x")]


def test_semantic_layer_survives_garbage():
    for garbage in ({}, {"nodes": 3, "members": [1, 2]}, {"members": {"a": None}, "design": 5},
                    {"nodes": {"a": {"x": "q"}}, "members": {"m": {"i": "a", "j": 7}},
                     "loadCases": {"c": {"loads": [3, {"type": "nodal"}, {"type": "memberPoint", "member": 1}]}}}):
        validate(garbage)                                    # must not raise
        validate_semantics(garbage, strict=False)


def test_tau_b_needs_fy():
    d = base()
    del d["materials"]["S"]["Fy"]
    d["combinations"]["C1"].update(stiffnessReduction=True, analysisType="secondOrder")
    got = [(p.code, p.path) for p in validate(d) if p.severity == "error"]
    assert got == [("analysis.tau-b-needs-fy", "/materials/S")]
    d["analysis"] = {"stiffnessReduction": {"tauB": "none"}}
    assert not [p for p in validate(d) if p.severity == "error"]


def test_warnings_for_questionable_input():
    d = base()
    d["combinations"]["C1"]["stiffnessReduction"] = True                # with first-order analysis
    d["nodes"]["D"] = {"x": 4.0003, "y": 0, "z": -3}                   # almost vertical member
    d["members"]["M3"] = {"i": "D", "j": "B", "section": "R", "material": "S"}
    d["materials"]["S"]["G"] = 5e7                                      # inconsistent with nu
    d["loadCases"]["unused"] = {"loads": []}
    d["design"] = {"defaults": {"method": "ASD"}}
    d["analysis"] = {"type": "secondOrder"}
    d["combinations"]["C1"]["analysisType"] = "firstOrder"
    d["combinations"]["C2"] = {"factors": {"D": 1}}
    got = sorted((p.code, p.path, p.severity) for p in validate(d))
    assert got == [
        ("case.unused", "/loadCases/unused", "warning"),
        ("combination.reduction-first-order", "/combinations/C1/stiffnessReduction", "warning"),
        ("design.asd-alpha", "/analysis", "warning"),
        ("material.g-nu-inconsistent", "/materials/S/G", "warning"),
        ("member.near-vertical", "/members/M3", "warning"),
    ]


def test_release_combinations_that_are_fine():
    d = base()
    d["members"]["M1"]["releases"] = {"i": {"rz": True, "ry": True}, "j": {"rz": True, "ry": True}}
    d["members"]["M2"]["releases"] = {"i": {"uy": True, "rz": True}, "j": {"ux": {"spring": 5}}}
    assert validate(d) == []
    d["members"]["M2"]["releases"] = {"i": {"uy": True, "rz": True}, "j": {"rz": True}}
    assert [(p.code, p.path) for p in validate(d)] == [("member.release-mechanism", "/members/M2/releases")]
    d["members"]["M2"]["releases"] = {"i": {"rz": {"spring": 0}}, "j": {"rz": True, "uy": True}}
    assert [(p.code, p.path) for p in validate(d)] == [("member.release-mechanism", "/members/M2/releases")]


def test_offsets_in_local_axes_need_distinct_nodes():
    d = base()
    d["nodes"]["B2"] = {"x": 4, "y": 0}
    d["members"]["M3"] = {"i": "B", "j": "B2", "section": "R", "material": "S",
                          "offsets": {"axes": "local", "j": [1, 0, 0]}}
    assert ("member.zero-length", "/members/M3") in [(p.code, p.path) for p in validate(d)]
    d["members"]["M3"]["offsets"]["axes"] = "global"           # global offset separates the ends: accepted
    assert not [p for p in validate(d) if p.severity == "error"]


def test_newer_minor_version_is_refused_and_patch_is_accepted():
    d = base()
    d["formatVersion"] = "1.0.7"
    assert validate(d) == []
    d["formatVersion"] = "1.1.0"
    assert [p.code for p in validate(d)] == ["version.unsupported"]


def test_problem_formatting():
    p = Problem("/a/b", "bad", "fix it", "error", "x.y", 3, 4)
    assert p.format("f.json") == "f.json:3:4: error [x.y] at /a/b: bad\n    hint: fix it"
    assert p.to_dict()["path"] == "/a/b"


def test_validate_schema_only_layer():
    d = base()
    d["extra"] = 1
    (p,) = validate_schema(d)
    assert (p.code, p.path) == ("schema.unknown-property", "/extra")


# ------------------------------------------------------------------------------------ CLI
def test_cli_validate(fixtures, capsys):
    import json
    from benchfea.cli import main

    assert main(["validate", str(fixtures / "valid" / "beam.yaml")]) == 0
    assert "valid (0 error(s), 0 warning(s))" in capsys.readouterr().out
    assert main(["validate", str(fixtures / "invalid" / "references.json")]) == 1
    out = capsys.readouterr().out
    assert "references.json:17:27: error [ref.unknown-node] at /members/M2/j" in out and "hint:" in out
    assert "INVALID (18 error(s), 1 warning(s))" in out
    assert main(["validate", "--json", str(fixtures / "invalid" / "syntax.json")]) == 1
    rep = json.loads(capsys.readouterr().out)
    assert rep["valid"] is False and rep["problems"][0]["code"] == "parse.json" and rep["problems"][0]["line"] == 4
    assert main(["validate", "--json", "--defaults", str(fixtures / "valid" / "beam.json")]) == 0
    text = capsys.readouterr().out
    assert '"valid": true' in text and '"/design/members/m1/Kx"' in text


def test_cli_help_and_unimplemented_commands(capsys):
    from benchfea.cli import build_parser, main

    for cmd in ("validate", "run", "design", "perturb", "verify", "score"):
        with pytest.raises(SystemExit) as e:
            build_parser().parse_args([cmd, "-h"])
        assert e.value.code == 0
        assert "example:" in capsys.readouterr().out
    assert main(["run", "model.json"]) == 2
