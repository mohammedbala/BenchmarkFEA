"""Independent checks of the input-format foundation (review pass).

Each block re-derives a result by a method different from the one the engine
uses: Gram-Schmidt + Rodrigues for the local axes, strip integration for the
section properties, and hand-written tricky documents for the loaders and the
validator.
"""
import copy
import math

import numpy as np
import pytest
from scipy.integrate import quad

from benchfea import io, sections
from benchfea.geometry import local_triad
from benchfea.model import build_model
from benchfea.problems import ModelInvalid
from benchfea.validate import validate

from helpers import assert_same

AXES = {"X": np.array([1.0, 0, 0]), "Y": np.array([0, 1.0, 0]), "Z": np.array([0, 0, 1.0])}


def base():
    return {
        "formatVersion": "1.0.0", "units": {"force": "kip", "length": "in"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0}, "B": {"x": 100, "y": 0}, "C": {"x": 100, "y": 0, "z": 120}},
        "supports": {"A": {"type": "fixed"}},
        "materials": {"S": {"E": 29000, "nu": 0.3, "Fy": 50, "Fu": 65}, "T": {"E": 10000, "G": 4000}},
        "sections": {"R": {"type": "rectangle", "h": 12, "b": 8}},
        "members": {"M1": {"i": "A", "j": "B", "section": "R", "material": "S"},
                    "M2": {"i": "B", "j": "C", "section": "R", "material": "T"}},
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fz": -1}]}},
        "combinations": {"C1": {"factors": {"D": 1.0}}},
    }


# ------------------------------------------------------------------------------ local axes
def independent_triad(p_i, p_j, vertical, ref=None, roll=0.0):
    """Local y = the part of the reference vector perpendicular to x (Gram-Schmidt),
    z = x cross y, then a Rodrigues rotation about x by the roll angle."""
    x = np.asarray(p_j, float) - np.asarray(p_i, float)
    x /= np.linalg.norm(x)
    if ref is None:
        ref = AXES[vertical] if np.linalg.norm(np.cross(x, AXES[vertical])) > 1e-6 else AXES["X"]
    v = np.asarray(ref, float)
    y = v - (v @ x) * x
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    th = math.radians(roll)
    K = np.array([[0, -x[2], x[1]], [x[2], 0, -x[0]], [-x[1], x[0], 0]])
    rot = np.eye(3) + math.sin(th) * K + (1 - math.cos(th)) * K @ K
    return np.vstack([x, rot @ y, rot @ z])


def test_triads_of_skewed_members_with_roll_and_reference_vectors_against_independent_construction():
    rng = np.random.default_rng(7)
    for _ in range(2000):
        p_i = rng.normal(size=3) * 10
        p_j = p_i + rng.normal(size=3) * rng.choice([1e-3, 1.0, 1e3])
        vertical = str(rng.choice(["Y", "Z"]))
        roll = float(rng.choice([0, 30, 90, -90, 180, 270, 360, rng.uniform(-720, 720)]))
        ref = None if rng.random() < 0.5 else rng.normal(size=3)
        R = local_triad(p_i, p_j, vertical, ref, roll)
        np.testing.assert_allclose(R, independent_triad(p_i, p_j, vertical, ref, roll), atol=1e-12)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-13)
        assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-13)             # right-handed
        np.testing.assert_allclose(np.cross(R[0], R[1]), R[2], atol=1e-13)   # x cross y = z


@pytest.mark.parametrize("vertical", ["Y", "Z"])
@pytest.mark.parametrize("sign", [1, -1])
@pytest.mark.parametrize("roll", [0.0, 90.0, 37.0, -90.0, 180.0])
def test_triads_of_vertical_members(vertical, sign, roll):
    p_i = np.array([1.0, 2.0, 3.0])
    p_j = p_i + sign * 7.0 * AXES[vertical]
    R = local_triad(p_i, p_j, vertical, None, roll)
    np.testing.assert_allclose(R, independent_triad(p_i, p_j, vertical, None, roll), atol=1e-14)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-14)
    np.testing.assert_allclose(R[0], sign * AXES[vertical], atol=1e-15)
    if roll == 0.0:
        np.testing.assert_allclose(R[1], AXES["X"], atol=1e-15)             # local y = global +X
    if roll == 90.0:                                                        # y' = z, z' = -y
        np.testing.assert_allclose(R[1], np.cross(sign * AXES[vertical], AXES["X"]), atol=1e-15)
        np.testing.assert_allclose(R[2], -AXES["X"], atol=1e-15)


def test_vertical_detection_threshold_and_near_vertical_warning():
    # sine of the angle to the vertical <= 1e-6: treated as vertical (y = +X, also when leaning towards +X)
    for eps in (0.0, 5e-7, 1e-6):
        np.testing.assert_allclose(local_triad((0, 0, 0), (eps, 0, 1), "Z")[1], (1, 0, -eps), atol=1e-12)
    # beyond the threshold the general rule applies: y points up-ish, i.e. to -X for a member leaning to +X
    np.testing.assert_allclose(local_triad((0, 0, 0), (1e-4, 0, 1), "Z")[1], (-1, 0, 1e-4), atol=1e-8)
    d = base()
    d["nodes"]["C"]["x"] = 100.01                           # column 0.01 in 120 out of plumb
    codes = [(p.code, p.path, p.severity) for p in validate(d)]
    assert ("member.near-vertical", "/members/M2", "warning") in codes


# ----------------------------------------------------------- section properties by strips
def strip_properties(outer, inner, half_depth, breaks):
    """A, Iz, Iy, Zz, Zy of a doubly symmetric section whose strip at height y is a
    centred segment of width outer(y) with a centred hole of width inner(y)."""
    pts = sorted({-half_depth, half_depth, 0.0, *breaks, *(-b for b in breaks)})

    def integral(f):
        return sum(quad(f, a, b, epsabs=0, epsrel=1e-13, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))

    return {
        "A": integral(lambda y: outer(y) - inner(y)),
        "Iz": integral(lambda y: y * y * (outer(y) - inner(y))),
        "Iy": integral(lambda y: (outer(y)**3 - inner(y)**3) / 12.0),
        "Zz": integral(lambda y: abs(y) * (outer(y) - inner(y))),
        "Zy": integral(lambda y: (outer(y)**2 - inner(y)**2) / 4.0),
    }


def compare(props, ref, c_y, c_z):
    ref = dict(ref, Sz=ref["Iz"] / c_y, Sy=ref["Iy"] / c_z)
    for name, value in ref.items():
        assert getattr(props, name) == pytest.approx(value, rel=1e-10), name
    assert props.rz == pytest.approx(math.sqrt(ref["Iz"] / ref["A"]), rel=1e-10)
    assert props.ry == pytest.approx(math.sqrt(ref["Iy"] / ref["A"]), rel=1e-10)


@pytest.mark.parametrize("d,bf,tf,tw,r", [(300, 150, 12, 8, 15), (14, 10, 0.7, 0.44, 0.5), (18, 7.5, 0.57, 0.355, 0),
                                          (200, 200, 15, 9, 18), (100, 60, 40, 10, 10)])
def test_i_shape_by_strip_integration(d, bf, tf, tw, r):
    hw = d - 2 * tf

    def outer(y):
        a = abs(y)
        if a > hw / 2:
            return bf
        s = hw / 2 - a                                      # distance from the inner face of the flange
        return tw + 2 * (r - math.sqrt(max(r * r - (r - s)**2, 0.0))) if s < r else tw

    ref = strip_properties(outer, lambda y: 0.0, d / 2, [hw / 2, hw / 2 - r])
    compare(sections.i_shape(d, bf, tf, tw, r), ref, d / 2, bf / 2)


@pytest.mark.parametrize("H,B,t", [(8, 6, 0.375), (200, 100, 10), (100, 300, 20)])
def test_box_by_strip_integration(H, B, t):
    ref = strip_properties(lambda y: B, lambda y: B - 2 * t if abs(y) < H / 2 - t else 0.0, H / 2, [H / 2 - t])
    compare(sections.box(H, B, t), ref, H / 2, B / 2)


@pytest.mark.parametrize("D,t", [(6, 0.28), (200, 10), (100, 45)])
def test_pipe_by_strip_integration(D, t):
    ro, ri = D / 2, D / 2 - t
    ref = strip_properties(lambda y: 2 * math.sqrt(max(ro * ro - y * y, 0.0)),
                           lambda y: 2 * math.sqrt(max(ri * ri - y * y, 0.0)), ro, [ri])
    p = sections.pipe(D, t)
    compare(p, ref, ro, ro)
    assert p.J == pytest.approx(ref["Iy"] + ref["Iz"], rel=1e-10)            # polar moment (exact for a tube)


@pytest.mark.parametrize("h,b", [(2, 1), (1, 3), (5, 5)])
def test_rectangle_by_strip_integration(h, b):
    compare(sections.rectangle(h, b), strip_properties(lambda y: b, lambda y: 0.0, h / 2, []), h / 2, b / 2)


@pytest.mark.parametrize("h,b", [(1, 1), (2, 1), (1, 4), (10, 1)])
def test_rectangle_torsion_constant_against_double_sine_series(h, b):
    """J = 2 * integral of the Prandtl stress function, expanded in the double
    sine series of the Poisson problem (a different series than the engine's)."""
    m = np.arange(1, 800, 2)[:, None]
    n = np.arange(1, 800, 2)[None, :]
    J = 256 / math.pi**6 * h * b * np.sum(1 / (m * m * n * n * ((m / h)**2 + (n / b)**2)))
    assert sections.rectangle_torsion_constant(h, b) == pytest.approx(J, rel=2e-6)   # truncated series


# ------------------------------------------------------------------------------- loaders
YAML_HEAD = "formatVersion: '1.0.0'\nunits: {force: kN, length: m}\nconventions: {verticalAxis: Z}\n"


@pytest.mark.parametrize("literal,read_as", [("010", 8), ("0x10", 16), ("1_000", 1000), ("1:30", 90), ("0o17", None)])
def test_yaml_numbers_that_are_not_plain_decimals_are_refused(literal, read_as):
    doc = io.loads(YAML_HEAD + f"nodes:\n  a: {{x: {literal}, y: 0}}\n", "yaml")
    if read_as is None:                                     # PyYAML reads this one as a string
        assert doc.data["nodes"]["a"]["x"] == literal and not doc.parse_problems
        return
    (p,) = doc.parse_problems
    assert (p.code, p.path, p.line, p.col) == ("parse.value", "/nodes/a/x", 5, 10)
    assert literal in p.message and repr(read_as) in p.message and "JSON" in p.hint


def test_yaml_accepts_the_json_number_grammar_and_quoted_strings():
    doc = io.loads("a: [0, -0, 10, -3.5, 1e3, 1E-3, 2.5e+2, +5, .5, 5., 0.0]\nb: '010'\nc: \"1_000\"\n", "yaml")
    assert not doc.parse_problems
    assert doc.data == {"a": [0, 0, 10, -3.5, 1000.0, 0.001, 250.0, 5, 0.5, 5.0, 0.0], "b": "010", "c": "1_000"}


def test_yaml_alias_is_located_at_its_own_key():
    doc = io.loads("nodes:\n  a: &p {x: 1, y: 2}\n  b: *p\n", "yaml")
    assert doc.data["nodes"]["b"] == {"x": 1, "y": 2} and not doc.parse_problems
    assert doc.locate("/nodes/a") == (2, 6) and doc.locate("/nodes/b") == (3, 3)


def test_empty_documents_and_wrong_top_level_type_are_located():
    for text, fmt in (("", "yaml"), ("# only a comment\n", "yaml"), ("[]", "json"), ("null", "json"), ("3", "json")):
        (p,) = validate(io.loads(text, fmt))
        assert (p.code, p.path, p.line, p.col) == ("schema.type", "", 1, 1), text


def test_byte_order_mark_bad_encoding_and_deep_nesting(tmp_path):
    f = tmp_path / "bom.json"
    f.write_bytes(b"\xef\xbb\xbf" + io.dumps(base()).encode())
    assert validate(f) == []
    f = tmp_path / "bom.yaml"
    f.write_bytes(b"\xef\xbb\xbf" + io.dumps(base(), "yaml").encode())
    assert validate(f) == []
    f = tmp_path / "latin1.json"
    f.write_bytes('{"name": "caf\xe9"}'.encode("latin-1"))
    with pytest.raises(ModelInvalid) as e:
        validate(f)
    assert e.value.problems[0].code == "parse.encoding" and "UTF-8" in e.value.problems[0].hint
    for text, fmt in (("[" * 100000, "json"), ("{\"a\":" * 100000, "json"), ("[" * 100000, "yaml")):
        with pytest.raises(ModelInvalid) as e:
            io.loads(text, fmt)
        assert e.value.problems[0].code in ("parse.json", "parse.yaml")


def test_dict_input_is_held_to_the_same_rules_as_text():
    d = base()
    d["nodes"]["A"]["x"] = float("nan")
    d["nodes"]["B"]["y"] = float("inf")
    d["nodes"][7] = {"x": 0, "y": 0}
    d["name"] = {1, 2}
    got = sorted((p.code, p.path) for p in validate(d) if p.severity == "error")
    assert got == [("parse.key-type", "/nodes"), ("parse.value", "/name"), ("parse.value", "/nodes/A/x"),
                   ("parse.value", "/nodes/B/y")]              # and no bogus geometry errors from the NaN
    with pytest.raises(ModelInvalid):
        build_model(d)
    # tuples and numpy numbers are plain data
    d = base()
    d["members"]["M1"]["offsets"] = {"i": (1, 0, 0)}
    d["nodes"]["C"]["z"] = np.float64(120.0)
    d["combinations"]["C1"]["factors"]["D"] = np.int64(2)
    m = build_model(d)
    assert m.members["M1"].length == pytest.approx(99 * 25.4) and m.combinations["C1"].factors["D"] == 2.0


# ----------------------------------------------------------------------------- validation
def test_structure_without_any_support_is_refused():
    for supports in (None, {}, {"A": {"dofs": {"ux": "free"}}}):
        d = base()
        if supports is None:
            del d["supports"]
        else:
            d["supports"] = supports
        errors = [(p.code, p.path) for p in validate(d) if p.severity == "error"]
        assert errors == [("support.none", "" if supports is None else "/supports")]
    d = base()
    d["supports"] = {"A": {"dofs": {"uz": {"spring": 10}}}}      # a spring is a support (stability: solver)
    assert not [p for p in validate(d) if p.severity == "error"]


def test_schema_messages_say_what_is_wrong_in_plain_words():
    d = base()
    d["members"]["M1"]["offsets"] = {"i": [1, 2]}
    d["constraints"] = {"c": {"type": "rigidLink", "master": "A", "slaves": ["B", "B"]}}
    d["combinations"]["C2"] = {"factors": {}}
    d["members"]["M2"]["i"] = False
    by = {p.path: p for p in validate(d)}
    assert by["/members/M1/offsets/i"].message == "list has 2 item(s), needs exactly 3"
    assert by["/constraints/c/slaves"].message == "list has repeated item(s): 'B'"
    assert by["/combinations/C2/factors"].message == "needs at least 1 entry, has 0"
    assert "quotes" in by["/members/M2/i"].hint
    # the repeated slave is reported once, not again by every semantic check
    assert [p.path for p in validate(d) if p.path.startswith("/constraints")] == ["/constraints/c/slaves"]


def test_one_malformed_design_entry_does_not_hide_the_other_design_problems():
    d = base()
    d["nodes"]["D"] = {"x": 200, "y": 0}
    d["members"]["M3"] = {"i": "B", "j": "D", "section": "R", "material": "S"}
    d["design"] = {
        "members": {"M1": {"members": []},                           # schema error (empty list)
                    "bad": 5,                                        # schema error (not an object)
                    "girder": {"members": ["M1", "M2"]},             # semantic: not collinear, not prismatic
                    "M3": {"Lx": "long", "bracing": {"major": ["x", 2.0]}}},
        "groups": {"g": {"members": ["girder", 3], "priority": "high", "Ly": 5, "bracing": {"minor": [0.5]}}},
        "defaults": {"deflectionLimits": [{"spanRatio": 300, "combinations": ["nope"]}, 7]},
    }
    got = {(p.code, p.path) for p in validate(d) if p.severity == "error"}
    assert {("design.not-collinear", "/design/members/girder/members/1"),
            ("design.not-prismatic", "/design/members/girder/members/1"),
            ("design.bracing-range", "/design/members/M3/bracing/major/1"),
            ("design.length-and-bracing", "/design/members/M3/Lx"),
            ("design.length-and-bracing", "/design/groups/g/Ly"),
            ("ref.unknown-combination", "/design/defaults/deflectionLimits/0/combinations/0")} <= got
    assert {c for c, _ in got if c.startswith("schema.")} >= {"schema.minItems", "schema.type"}


def test_bracing_checks_do_not_need_a_design_member():
    d = base()
    d["design"] = {"groups": {"unused": {"members": ["nobody"], "Lz": 10, "bracing": {"torsion": [0.0, 0.5, 1.2]}}}}
    got = sorted((p.code, p.path) for p in validate(d) if p.severity == "error")
    assert got == [("design.bracing-range", "/design/groups/unused/bracing/torsion/0"),
                   ("design.bracing-range", "/design/groups/unused/bracing/torsion/2"),
                   ("design.length-and-bracing", "/design/groups/unused/Lz"),
                   ("ref.unknown-design-member", "/design/groups/unused/members/0")]


def test_length_and_bracing_from_two_groups_of_equal_priority():
    d = base()
    d["design"] = {"groups": {"g1": {"members": ["M1"], "Lb": 50}, "g2": {"members": ["M1"], "bracing": {"topFlange": [0.5]}}}}
    assert [(p.code, p.path) for p in validate(d) if p.severity == "error"] == [
        ("design.length-and-bracing", "/design/groups/g1/Lb")]
    d["design"]["groups"]["g2"]["priority"] = 1
    m = build_model(d)
    assert m.design_members["M1"].source("Lb") == "bracing" and m.design_members["M1"].braceTopFlange == (0.5,)


# ------------------------------------------------- material properties as design variables
KSI = 4448.2216152605 / 25.4**2


def test_design_material_properties_default_to_the_material_with_provenance():
    m = build_model(base())
    m1, m2 = m.design_members["M1"], m.design_members["M2"]
    assert (m1.E, m1.Fy, m1.Fu) == (pytest.approx(29000 * KSI), pytest.approx(50 * KSI), pytest.approx(65 * KSI))
    assert m1.G == pytest.approx(29000 * KSI / 2.6) and m1.G == m.materials["S"].G
    for name in ("E", "G", "Fy", "Fu"):
        assert m1.source(name) == "material:S" and m1.values[name].path == f"/materials/S/{name}"
        assert m2.source(name) == "material:T"
    assert (m2.E, m2.G, m2.Fy, m2.Fu) == (pytest.approx(10000 * KSI), pytest.approx(4000 * KSI), None, None)
    assert "no Fy" in m2.values["Fy"].reason
    # values taken from the material are not defaults
    assert not [x for x in m.defaults_used if x.path.endswith(("/E", "/Fy", "/Fu")) or "/design/" in x.path
                and x.path.endswith("/G")]


def test_design_material_overrides_resolve_by_level_convert_units_and_never_touch_the_analysis():
    plain = base()
    plain["design"] = {}
    over = copy.deepcopy(plain)
    over["design"] = {"defaults": {"Fy": 36, "Fu": 58},
                      "groups": {"g": {"members": ["M1", "M2"], "Fy": 42, "E": 29000}},
                      "members": {"M2": {"G": 3000, "Fy": 46}}}
    a, b = build_model(plain), build_model(over)
    m1, m2 = b.design_members["M1"], b.design_members["M2"]
    assert (m1.Fy, m1.source("Fy")) == (pytest.approx(42 * KSI), "group:g")
    assert (m2.Fy, m2.source("Fy")) == (pytest.approx(46 * KSI), "member")
    assert (m1.Fu, m1.source("Fu")) == (pytest.approx(58 * KSI), "global")
    assert (m2.Fu, m2.source("Fu")) == (pytest.approx(58 * KSI), "global")      # the material had none
    assert (m2.E, m2.source("E")) == (pytest.approx(29000 * KSI), "group:g")
    assert (m2.G, m2.source("G"), m2.values["G"].path) == (pytest.approx(3000 * KSI), "member", "/design/members/M2/G")
    assert (m1.G, m1.source("G")) == (pytest.approx(29000 * KSI / 2.6), "material:S")
    # the analysis objects are identical: the overrides are design-only
    assert_same(a, b, ignore=("design_members", "defaults_used", "warnings"))
    assert b.materials["T"].E == pytest.approx(10000 * KSI) and b.materials["T"].G == pytest.approx(4000 * KSI)
    assert b.materials["S"].Fy == pytest.approx(50 * KSI) and b.materials["T"].Fy is None
    # a design stiffness that differs from the analysis stiffness is flagged (the equal one is not)
    warned = sorted((p.code, p.path) for p in b.warnings)
    assert warned == [("design.stiffness-differs", "/design/groups/g/E"),
                      ("design.stiffness-differs", "/design/members/M2/G")]
    # without the overrides the missing strength of material T is reported instead
    assert [(p.code, p.path) for p in a.warnings] == [("design.missing-strength", "/materials/T")]


def test_design_block_is_refused_for_bad_material_overrides():
    d = base()
    d["design"] = {"defaults": {"Fy": 0, "E": "steel", "Fz": 1}}
    got = sorted((p.code, p.path) for p in validate(d) if p.severity == "error")
    assert got == [("schema.range", "/design/defaults/Fy"), ("schema.type", "/design/defaults/E"),
                   ("schema.unknown-property", "/design/defaults/Fz")]


# -------------------------------------------------------------------- defaults-used report
def test_defaults_report_lists_geometry_orientation_and_load_defaults():
    d = base()
    d["materials"]["S"]["unitWeight"] = 2.84e-4
    d["members"]["M1"]["offsets"] = {"i": [2, 0, 0]}
    d["loadCases"]["D"]["loads"] += [{"type": "selfWeight"},
                                     {"type": "memberPoint", "member": "M1", "direction": "Z", "value": -1,
                                      "position": 0.5}]
    d["combinations"]["C1"].update(analysisType="secondOrder", stiffnessReduction=True)
    d["analysis"] = {"stiffnessReduction": {"tauB": "none"}}
    d["design"] = {"groups": {"g": {"members": ["M1"], "bracing": {"minor": [0.5]}}}}
    rep = {x.path: x for x in build_model(d).defaults_used}
    assert rep["/nodes/A/z"].value == 0.0 and rep["/nodes/B/z"].value == 0.0 and "/nodes/C/z" not in rep
    assert rep["/members/M1/orientation/referenceVector"].value == [0.0, 0.0, 1.0]
    assert rep["/members/M2/orientation/referenceVector"].value == [1.0, 0.0, 0.0]      # vertical member
    assert "vertical member" in rep["/members/M2/orientation/referenceVector"].reason
    assert rep["/members/M1/offsets/axes"].value == "global"
    assert rep["/loadCases/D/loads/1/factor"].value == 1.0
    assert rep["/loadCases/D/loads/1/members"].value == ["M1", "M2"]
    assert rep["/loadCases/D/loads/1/direction"].value == [0.0, 0.0, -1.0]
    assert rep["/loadCases/D/loads/2/kind"].value == "force"
    assert rep["/analysis/stiffnessReduction/members"].value == ["M1", "M2"]
    assert rep["/design/autoMembers"].value is True
    assert rep["/design/groups/g/bracing/positionType"].value == "relative"
    # nothing that the file states is reported
    d["members"]["M1"]["orientation"] = {"referenceVector": [0, 0, 1]}
    d["members"]["M1"]["offsets"]["axes"] = "local"
    d["design"]["autoMembers"] = True
    d["design"]["groups"]["g"]["bracing"]["positionType"] = "relative"
    rep = {x.path for x in build_model(d).defaults_used}
    assert not rep & {"/members/M1/orientation/referenceVector", "/members/M1/offsets/axes", "/design/autoMembers",
                      "/design/groups/g/bracing/positionType"}


def test_design_document_reports_the_design_variables_in_file_units():
    from benchfea.results import DesignMemberResult, DesignResults
    d = base()                                                       # kip, in
    d["design"] = {"members": {"M1": {"Lx": 60, "An": 90, "Fy": 36, "bracing": {"minor": [0.25]},
                                      "deflectionLimits": [{"absolute": 0.5}, {"spanRatio": 240}]}}}
    m = build_model(d)
    assert m.design_members["M1"].Lx == pytest.approx(60 * 25.4)      # internal: mm
    doc = DesignResults("AISC 360-22", {"M1": DesignMemberResult("M1", [], True)}).to_document(m)
    par = doc["members"]["M1"]["parameters"]
    assert doc["units"] == {"force": "kip", "length": "in"}
    assert (par["Lx"]["value"], par["Lx"]["unit"], par["Lx"]["source"]) == (pytest.approx(60.0), "in", "member")
    assert (par["An"]["value"], par["An"]["unit"]) == (pytest.approx(90.0), "in^2")
    assert (par["Fy"]["value"], par["Fy"]["unit"]) == (pytest.approx(36.0), "kip/in^2")
    assert (par["E"]["value"], par["E"]["source"]) == (pytest.approx(29000.0), "material:S")
    assert par["Ly"]["value"] is None and par["Ly"]["source"] == "bracing"
    assert par["Lz"]["value"] == pytest.approx(100.0) and par["braceMinor"]["value"] == [0.25]
    assert par["Kx"] == {"value": 1.0, "unit": "", "source": "default", "path": "/design/members/M1/Kx"}
    lims = par["deflectionLimits"]["value"]
    assert lims[0]["absolute"] == pytest.approx(0.5) and lims[1]["absolute"] is None and lims[1]["spanRatio"] == 240
    assert lims[0] == {"spanRatio": None, "absolute": pytest.approx(0.5), "direction": "resultant", "combinations": []}
    import json
    json.dumps(doc)
