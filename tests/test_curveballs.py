"""Curve balls (``benchfea.bench.curveballs``): the change record, determinism, measured materiality,
applicability per generator family, design-only invariance, the CLI.

The base models are generated here (``benchfea.bench.generate``), the cheapest model of each family
that has the feature a curve ball needs.  ``EXPECTED`` names, for every curve ball, one generated
model on which it must be applicable and the effect class it must give there; the same pairs are
rows of the applicability table that ``tools/build_curveball_docs.py`` writes into
docs/CURVEBALLS.md.
"""
import copy
import functools
import json
import re
from pathlib import Path

import numpy as np
import pytest

from benchfea import design as design_module
from benchfea.bench import curveballs as cb
from benchfea.bench import families, generate, tier_c_expected
from benchfea.bench.curveballs import measure
from benchfea.bench.curveballs.record import Change, ChangeRecord, RecordMismatch, apply_changes, diff, invert_changes
from benchfea.cli import main
from benchfea.model import build_model
from benchfea.problems import ModelInvalid, ModelUnstable
from benchfea.solve import analyze
from benchfea.validate import validate

ROOT = Path(__file__).resolve().parents[1]

# key -> generate() arguments (family, seed, difficulty, forced features)
MODELS = {
    "beam1": ("continuous_beam", 1, 1, {}),
    "beam2": ("continuous_beam", 1, 2, {}),
    "beamOverhang": ("continuous_beam", 1, 2, {"overhang": True}),
    "portal1": ("portal_frame", 1, 1, {}),
    "portalOffsets": ("portal_frame", 2, 2, {"rigidOffsets": True}),
    "portalSecondOrder": ("portal_frame", 1, 3, {}),
    "frame2d": ("multistorey_2d", 1, 2, {}),
    "frame3d": ("multistorey_3d", 1, 3, {}),
    "braced": ("braced_frame", 1, 2, {}),
    "rack": ("pipe_rack", 1, 3, {}),
    "stair": ("stair_frame", 1, 1, {}),
    "diaphragm": ("diaphragm_frame", 1, 3, {}),
}

# curve ball -> (model key, effect class it must give there with seed 0)
EXPECTED = {
    "constraint_removed": ("diaphragm", "creates-mechanism"),
    "constraint_added": ("frame3d", "changes-analysis"),
    "constraint_slaves_dropped": ("diaphragm", "changes-analysis"),
    "constraint_wrong_slave": ("diaphragm", "invalid-model"),
    "constraint_wrong_master": ("diaphragm", "creates-mechanism"),
    "constraint_made_flexible": ("diaphragm", "changes-analysis"),
    "offset_added": ("portal1", "changes-analysis"),
    "offset_removed": ("portalOffsets", "changes-analysis"),
    "offset_sign_reversed": ("portalOffsets", "changes-analysis"),
    "offset_wrong_axis": ("portalOffsets", "changes-analysis"),
    "release_added": ("portal1", "changes-analysis"),
    "release_removed": ("beam2", "changes-analysis"),
    "release_mechanism": ("portal1", "creates-mechanism"),
    "release_pinned_frame": ("frame2d", "creates-mechanism"),
    "design_unbraced_length_changed": ("portal1", "design-only"),
    "design_k_changed": ("portal1", "design-only"),
    "design_bracing_points_changed": ("beam2", "design-only"),
    "design_lb_full_where_braced": ("beam2", "design-only"),
    "design_lb_short_where_unbraced": ("beamOverhang", "design-only"),
    "support_fixed_to_pinned": ("portal1", "changes-analysis"),
    "support_pinned_to_fixed": ("beam1", "changes-analysis"),
    "support_roller_direction": ("beam2", "changes-analysis"),
    "support_spring_stiffness": ("beam2", "changes-analysis"),
    "support_removed": ("beam2", "changes-analysis"),
    "member_rolled_90": ("portal1", "changes-analysis"),
    "section_lighter": ("portal1", "changes-analysis"),
    "section_heavier": ("portal1", "changes-analysis"),
    "sections_swapped": ("portal1", "changes-analysis"),
    "units_relabelled": ("portal1", "changes-analysis"),
    "load_sign_flipped": ("portal1", "changes-analysis"),
    "load_direction_changed": ("portal1", "changes-analysis"),
    "load_moved": ("portal1", "changes-analysis"),
    "second_order_off": ("portalSecondOrder", "changes-analysis"),
}
# one cheap curve ball per family: every family takes curve balls
FAMILY_MODEL = {"continuous_beam": "beam1", "portal_frame": "portal1", "multistorey_2d": "frame2d",
                "multistorey_3d": "frame3d", "braced_frame": "braced", "pipe_rack": "rack", "stair_frame": "stair",
                "diaphragm_frame": "diaphragm"}

TASK_LIST = {
    "constraints": ["constraint_removed", "constraint_added", "constraint_slaves_dropped", "constraint_wrong_slave",
                    "constraint_wrong_master", "constraint_made_flexible"],
    "offsets": ["offset_added", "offset_removed", "offset_sign_reversed", "offset_wrong_axis"],
    "releases": ["release_added", "release_removed", "release_mechanism", "release_pinned_frame"],
    "design": ["design_unbraced_length_changed", "design_k_changed", "design_bracing_points_changed",
               "design_lb_full_where_braced", "design_lb_short_where_unbraced"],
    "supports": ["support_fixed_to_pinned", "support_pinned_to_fixed", "support_roller_direction",
                 "support_spring_stiffness", "support_removed"],
    "members": ["member_rolled_90", "section_lighter", "section_heavier", "sections_swapped"],
    "units": ["units_relabelled"],
    "loads": ["load_sign_flipped", "load_direction_changed", "load_moved"],
    "second_order": ["second_order_off"],
}


@functools.lru_cache(maxsize=None)
def _model_json(key: str) -> str:
    family, seed, level, features = MODELS[key]
    return generate(family, seed, level, features=dict(features)).json()


def model(key: str) -> dict:
    return json.loads(_model_json(key))


@functools.lru_cache(maxsize=None)
def _applied(name: str, key: str, seed: int = 0):
    return cb.apply(model(key), name, seed)


def text(doc) -> str:
    """Canonical text: equal texts <=> equal documents, key order included."""
    return json.dumps(doc)


def run(doc):
    m = build_model(copy.deepcopy(doc))
    return m, analyze(m)


# =========================================================================================== record
def test_diff_apply_invert_roundtrip_with_key_order():
    base = {"a": 1, "b": {"x": [1, 2, 3], "y": "s", "z": {"k": 1}}, "c": [{"p": 1}, {"p": 2}], "d": 4.0}
    new = {"a": 1, "n0": True, "b": {"x": [1, 2], "z": {"k": 2, "m": None}, "w": 0}, "c": [{"p": 1}, {"p": 3, "q": 1}],
           "n1": {"deep": [1]}}
    changes = diff(base, new)
    ops = {(c.op, c.path) for c in changes}
    assert ("remove", "/d") in ops and ("add", "/n0") in ops and ("add", "/n1") in ops
    assert ("replace", "/b/x") in ops                       # array of another length: replaced as a whole
    assert ("remove", "/b/y") in ops and ("add", "/b/w") in ops and ("add", "/b/z/m") in ops
    assert ("replace", "/b/z/k") in ops and ("replace", "/c/1/p") in ops and ("add", "/c/1/q") in ops
    assert ("replace", "/a") not in ops
    assert text(apply_changes(base, changes)) == text(new)                       # forward, key order included
    assert text(apply_changes(new, invert_changes(changes))) == text(base)       # and back
    assert base["b"]["x"] == [1, 2, 3]                                           # inputs untouched
    # through JSON
    again = [Change.from_dict(json.loads(json.dumps(c.to_dict()))) for c in changes]
    assert text(apply_changes(base, again)) == text(new)


def test_diff_keeps_insertion_position_of_object_members():
    base = {"n1": 1, "n2": 2, "n3": 3}
    new = {"n0": 0, "n1": 1, "nx": 9, "n3": 3}
    changes = diff(base, new)
    assert [(c.op, c.path, c.index) for c in changes] == [("remove", "/n2", 1), ("add", "/n0", 0), ("add", "/nx", 2)]
    assert list(apply_changes(base, changes)) == ["n0", "n1", "nx", "n3"]
    assert list(apply_changes(new, invert_changes(changes))) == ["n1", "n2", "n3"]


def test_diff_escapes_pointer_tokens_and_distinguishes_number_spelling():
    base = {"a/b": {"~k": 1}, "v": 1}
    new = {"a/b": {"~k": 2}, "v": 1.0}
    changes = diff(base, new)
    assert {c.path for c in changes} == {"/a~1b/~0k", "/v"}
    assert text(apply_changes(base, changes)) == text(new)
    assert diff(base, copy.deepcopy(base)) == []


def test_apply_changes_checks_before_values_and_paths():
    base = {"a": {"x": 1}, "l": [1, 2]}
    with pytest.raises(RecordMismatch, match="holds"):
        apply_changes(base, [Change("replace", "/a/x", 5, 6)])
    with pytest.raises(RecordMismatch, match="does not exist"):
        apply_changes(base, [Change("remove", "/a/y", before=1)])
    with pytest.raises(RecordMismatch, match="already exists"):
        apply_changes(base, [Change("add", "/a/x", after=1)])
    with pytest.raises(RecordMismatch, match="never added or removed"):
        apply_changes(base, [Change("remove", "/l/0", before=1)])
    with pytest.raises(RecordMismatch):
        apply_changes(base, [Change("replace", "/nope/x", 1, 2)])
    assert apply_changes(base, [Change("replace", "/l/1", 2, 7)]) == {"a": {"x": 1}, "l": [1, 7]}


def test_shared_sub_objects_of_a_python_document_are_not_edited_together():
    """A document built in Python may use ONE dict in two places (the generators do, for the two end
    releases of a member).  An edit of one place must not leak into the other."""
    shared = {"rz": {"spring": 5.0}}
    base = {"members": {"B": {"releases": {"i": shared, "j": shared}}}}
    out = apply_changes(base, [Change("add", "/members/B/releases/i/ry", after=True, index=1)])
    assert out["members"]["B"]["releases"] == {"i": {"rz": {"spring": 5.0}, "ry": True}, "j": {"rz": {"spring": 5.0}}}
    assert shared == {"rz": {"spring": 5.0}}
    # the generated document itself (not its JSON text) goes through a curve ball unharmed
    family, seed, level, features = MODELS["portalSecondOrder"]
    g = generate(family, seed, level, features=dict(features))
    rel = g.document["members"]["B1"]["releases"]
    p = cb.apply(g.document, "release_mechanism", 0)
    assert text(p.record.apply(g.document)) == text(p.document)
    assert text(p.record.revert(p.document)) == text(json.loads(g.json()))
    assert g.document["members"]["B1"]["releases"] == rel == json.loads(g.json())["members"]["B1"]["releases"]
    assert p.record.to_dict() == cb.apply(json.loads(g.json()), "release_mechanism", 0).record.to_dict()


def test_reordered_object_is_replaced_as_a_whole():
    base, new = {"o": {"a": 1, "b": 2}}, {"o": {"b": 2, "a": 1}}
    changes = diff(base, new)
    assert [(c.op, c.path) for c in changes] == [("replace", "/o")]
    assert text(apply_changes(base, changes)) == text(new)


# ========================================================================================= registry
def test_every_curve_ball_of_the_task_list_is_registered_and_documented():
    names = cb.available()
    assert len(names) == len(set(names))
    for category, expected in TASK_LIST.items():
        assert [n for n in names if cb.get(n).category == category] == expected
    assert set(names) == {n for group in TASK_LIST.values() for n in group} == set(EXPECTED)
    assert set(cb.CATEGORIES) == set(TASK_LIST)
    for n in names:
        c = cb.get(n)
        assert c.effect_classes and set(c.effect_classes) <= set(cb.EFFECT_CLASSES)
        for field in (c.title, c.changes, c.why, c.moves, c.precondition):
            assert isinstance(field, str) and len(field) > 10
        assert len(c.why) > 120, f"{n}: say why it matters physically"
    # design variables are design-only and nothing else is
    for n in names:
        assert (cb.get(n).effect_classes == ("design-only",)) == (cb.get(n).category == "design")
    assert cb.get("release_mechanism").effect_classes == ("creates-mechanism",)
    assert cb.THRESHOLDS == {"analysisChange": 0.05, "designRatioChange": 0.05}


def test_unknown_curve_ball_and_impossible_effect_class():
    with pytest.raises(KeyError, match="unknown curve ball"):
        cb.apply(model("beam1"), "no_such_curve_ball", 0)
    with pytest.raises(cb.NotApplicable, match="never gives the effect class"):
        cb.apply(model("beam1"), "load_sign_flipped", 0, effect_class="design-only")


# ================================================================================ every curve ball
@pytest.mark.parametrize("name", list(EXPECTED))
def test_curve_ball_applies_record_describes_the_difference_and_effect_is_material(name):
    key, expected_class = EXPECTED[name]
    base = model(key)
    untouched = text(base)
    p = _applied(name, key)
    rec = p.record
    assert text(model(key)) == untouched                              # the base document is never modified
    assert rec.name == name and rec.seed == 0 and rec.category == cb.get(name).category
    assert rec.effect_class == expected_class and expected_class in cb.get(name).effect_classes
    assert rec.why == cb.get(name).why and rec.precondition == cb.get(name).precondition and rec.target

    # ---- the record IS the difference: forward and inverse, key order included
    assert rec.changes, "a curve ball changes something"
    assert text(rec.apply(base)) == text(p.document)
    assert text(rec.revert(p.document)) == untouched
    assert text(p.document) != untouched
    assert [c.to_dict() for c in diff(base, p.document)] == [c.to_dict() for c in rec.changes]
    for c in rec.changes:
        assert c.path.startswith("/") and c.op in ("replace", "add", "remove")
    # ... and it survives JSON
    d = json.loads(json.dumps(rec.to_dict()))
    assert d["format"] == cb.RECORD_FORMAT and d["curveBall"] == name and d["effectClass"] == expected_class
    again = ChangeRecord.from_dict(d)
    assert again.to_dict() == rec.to_dict()
    assert text(again.apply(base)) == text(p.document)

    # ---- the entities named exist in the base or in the perturbed model
    assert rec.entities
    table = {"members": "members", "nodes": "nodes", "supports": "supports", "constraints": "constraints",
             "sections": "sections", "loadCases": "loadCases", "combinations": "combinations"}
    for kind, ids in rec.entities.items():
        assert ids, kind
        if kind in table:
            for i in ids:
                assert i in base.get(table[kind], {}) or i in p.document.get(table[kind], {}), (kind, i)
    # every change path is under an entity kind the record names, or a section / design block it implies
    tops = {c.path.split("/")[1] for c in rec.changes}
    assert tops <= {"members", "supports", "constraints", "sections", "loadCases", "combinations", "analysis",
                    "design", "units"}

    # ---- the effect, measured again here from the two documents
    b = measure.evaluate(base)
    q = measure.evaluate(p.document, cache=False)
    assert b.ok
    e = rec.effect
    assert e["material"] is True and e["class"] == expected_class and e["thresholds"] == cb.THRESHOLDS
    if expected_class == "creates-mechanism":
        assert q.status == "unstable" and rec.reference["kind"] == "diagnosis"
        with pytest.raises(ModelUnstable) as info:
            run(p.document)
        assert {pr.code for pr in info.value.problems} == {"stability.mechanism"}
        assert info.value.diagnosis.to_dict() == rec.reference["diagnosis"]      # the reference IS the diagnosis
        assert rec.reference["diagnosis"]["nMechanisms"] >= 1 and rec.reference["diagnosis"]["modes"]
        assert "diagnosis" in rec.verification["route"]
    elif expected_class == "invalid-model":
        assert q.status == "invalid" and rec.reference["kind"] == "validation-problems"
        errors = [pr for pr in validate(p.document) if pr.severity == "error"]
        assert errors and {pr.code for pr in errors} == {pr["code"] for pr in rec.reference["problems"]}
        assert not [pr for pr in validate(base) if pr.severity == "error"]
        with pytest.raises(ModelInvalid):
            build_model(copy.deepcopy(p.document))
    else:
        assert q.ok and rec.reference["kind"] == "results"
        fresh = measure.compare(b, q)
        for k in ("analysisIdentical", "analysisChange", "designRatioChange", "designPassFailChanged"):
            assert fresh[k] == e[k], k
        if expected_class == "design-only":
            assert e["analysisIdentical"] is True and e["analysisChange"] == 0.0
            assert e["criterion"] in ("designRatioChange", "designPassFailChanged")
            assert e["designPassFailChanged"] or e["designRatioChange"] >= 0.05
        else:
            assert e["analysisIdentical"] is False
            assert e["criterion"] in cb.get(name).criteria_for("changes-analysis") and e["criterion"] in e["criteria"]
            if e["criterion"] == "analysisChange":
                assert e["analysisChange"] >= 0.05
            elif e["criterion"] == "designPassFailChanged":
                assert e["designPassFailChanged"]
            else:
                assert cb.get(name).category == "members" and e["designRatioChange"] >= 0.05
        # every perturbed model that has results passed the engine's own invariants (strict analysis)
        assert all(c.invariants.passed and c.converged for c in q.results.combinations.values())
    for dropped in rec.rejected:
        assert dropped["target"] and dropped["reason"]
    assert rec.base_hash != rec.perturbed_hash and len(rec.base_hash) == 64
    assert rec.verification["route"]


@pytest.mark.parametrize("name", list(EXPECTED))
def test_curve_ball_is_deterministic(name):
    key, _ = EXPECTED[name]
    first = _applied(name, key)
    second = cb.apply(model(key), name, 0)                # a fresh call on a fresh copy of the document
    assert text(second.document) == text(first.document)
    assert json.dumps(second.record.to_dict()) == json.dumps(first.record.to_dict())
    assert second.json() == first.json()


def test_seed_selects_the_target_and_nothing_else_decides():
    base = model("portal1")
    targets = {}
    for seed in range(6):
        p = cb.apply(base, "load_sign_flipped", seed)
        targets[seed] = p.record.target
        assert p.record.seed == seed
    assert len(set(targets.values())) >= 2                # different seeds reach different targets
    # the same seed again, after other seeds were used in between: no hidden state
    assert cb.apply(base, "load_sign_flipped", 3).record.target == targets[3]
    assert cb.apply(model("portal1"), "load_sign_flipped", 3).record.to_dict() == \
        cb.apply(base, "load_sign_flipped", 3).record.to_dict()


@pytest.mark.parametrize("family", list(FAMILY_MODEL))
def test_every_generator_family_takes_curve_balls(family):
    assert set(FAMILY_MODEL) == set(families())
    key = FAMILY_MODEL[family]
    assert MODELS[key][0] == family
    p = cb.apply(model(key), "load_sign_flipped", 0)
    assert p.record.effect_class == "changes-analysis" and p.record.effect["analysisChange"] >= 0.05
    assert text(p.record.revert(p.document)) == text(model(key))


def test_every_curve_ball_is_applicable_to_a_generator_family():
    covered = {MODELS[key][0] for key, _ in EXPECTED.values()}
    assert covered <= set(families())
    for name, (key, _) in EXPECTED.items():
        assert _applied(name, key).record.name == name


# ============================================================================================ refusals
@pytest.mark.parametrize("name, key, fragment", [
    ("constraint_removed", "portal1", "precondition not met: the model has at least one rigid constraint"),
    ("constraint_made_flexible", "beam1", "precondition not met"),
    ("offset_removed", "portal1", "precondition not met: a member with a nonzero rigid end offset"),
    ("offset_sign_reversed", "beam1", "precondition not met"),
    ("release_removed", "portal1", "precondition not met: a member with a release"),
    ("support_spring_stiffness", "portal1", "precondition not met: a support with a spring"),
    ("second_order_off", "portal1", "precondition not met: the model has second-order combinations"),
    ("sections_swapped", "beam1", "precondition not met: two members with different sections"),
    ("release_pinned_frame", "beam1", "precondition not met"),
])
def test_refuses_a_model_that_lacks_the_feature(name, key, fragment):
    with pytest.raises(cb.NotApplicable) as info:
        cb.apply(model(key), name, 0)
    assert fragment in info.value.reason and info.value.name == name
    assert info.value.tried == []


def test_refuses_when_no_target_matters_and_reports_what_was_measured():
    """A single-span beam carries no axial force: no effective length factor matters."""
    with pytest.raises(cb.NotApplicable) as info:
        cb.apply(model("beam1"), "design_k_changed", 0)
    e = info.value
    assert "material effect" in e.reason and "below the materiality threshold" in e.reason
    assert e.tried and all(t["analysisChange"] == 0.0 for t in e.tried)
    assert all(t["designRatioChange"] < 0.05 for t in e.tried)


def test_second_order_off_only_where_the_amplification_is_material():
    base = model("portalSecondOrder")
    p = _applied("second_order_off", "portalSecondOrder")
    e = p.record.effect
    assert e["criterion"] == "analysisChange" and e["analysisChange"] >= 0.05
    amps = e["baseAmplification"]
    assert amps and max(amps.values()) > 1.0
    assert all("analysisType" not in c for c in p.document["combinations"].values())
    assert any(c.get("analysisType") == "secondOrder" for c in base["combinations"].values())
    assert set(p.record.entities["combinations"]) == {cid for cid, c in base["combinations"].items()
                                                      if c.get("analysisType") == "secondOrder"}
    # the same frame under a twentieth of the loads: second-order effects are there but do not govern
    light = copy.deepcopy(base)
    for c in light["combinations"].values():
        c["factors"] = {k: v * 0.05 for k, v in c["factors"].items()}
    with pytest.raises(cb.NotApplicable) as info:
        cb.apply(light, "second_order_off", 0)
    assert "below the materiality threshold" in info.value.reason
    assert 0.0 < info.value.tried[0]["analysisChange"] < 0.05


def test_base_model_that_cannot_be_analysed_is_refused():
    broken = model("beam1")
    broken["supports"] = {k: v for k, v in list(broken["supports"].items())[:0]}
    with pytest.raises(cb.NotApplicable, match="base model is not analysable"):
        cb.apply(broken, "load_sign_flipped", 0)


# ================================================================================= design-only
@pytest.mark.parametrize("name", TASK_LIST["design"])
def test_design_only_curve_ball_leaves_the_analysis_bitwise_identical(name):
    key, _ = EXPECTED[name]
    base = model(key)
    p = _applied(name, key)
    assert p.record.effect_class == "design-only"
    assert {c.path.split("/")[1] for c in p.record.changes} == {"design"}       # nothing outside the design block
    (m0, r0), (m1, r1) = run(base), run(p.document)
    assert list(r0.combinations) == list(r1.combinations)
    for cid, c0 in r0.combinations.items():
        c1 = r1.combinations[cid]
        assert np.array_equal(c0.displacements, c1.displacements) and np.array_equal(c0.reactions, c1.reactions)
        for mid, a in c0.members.items():
            b = c1.members[mid]
            assert np.array_equal(a.end_forces, b.end_forces) and np.array_equal(a.forces, b.forces)
            assert np.array_equal(a.displacements, b.displacements) and np.array_equal(a.x, b.x)
    assert json.dumps(r0.to_document(m0)["combinations"]) == json.dumps(r1.to_document(m1)["combinations"])
    # ... and the design moved materially
    d0, d1 = design_module.check(m0, r0), design_module.check(m1, r1)
    changes = [abs(d1.members[k].governing.ratio - v.governing.ratio) / v.governing.ratio
               for k, v in d0.members.items() if v.governing is not None and d1.members[k].governing is not None]
    flipped = [k for k, v in d0.members.items() if (v.status == "fail") != (d1.members[k].status == "fail")]
    assert max(changes) >= 0.05 or flipped
    assert max(changes) == pytest.approx(p.record.effect["designRatioChange"], rel=1e-12)
    assert sorted(flipped) == sorted(p.record.effect["designPassFailChanged"])


# ====================================================================================== measure
def _cantilever(load=-10.0, design=False):
    doc = {"formatVersion": "1.0.0", "units": {"force": "kN", "length": "m"}, "conventions": {"verticalAxis": "Z"},
           "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": 3, "y": 0, "z": 0}},
           "supports": {"A": {"type": "fixed"}},
           "materials": {"s": {"E": 2.0e8, "G": 7.7e7, "Fy": 345000, "Fu": 450000}},
           "sections": {"sec": {"type": "I", "d": 0.3, "bf": 0.15, "tf": 0.01, "tw": 0.007}},
           "members": {"M": {"i": "A", "j": "B", "section": "sec", "material": "s"}},
           "loadCases": {"P": {"loads": [{"type": "nodal", "node": "B", "Fz": load}]}},
           "combinations": {"C1": {"factors": {"P": 1.0}}}}
    if design:
        doc["design"] = {"defaults": {"method": "LRFD"}}
    return doc


def test_measure_compare_on_a_hand_computable_case():
    """Tip load on a cantilever: doubling the load doubles every result (change = 1.0 of the base
    peak); flipping it reverses every result (change = 2.0)."""
    base = measure.evaluate(_cantilever(-10.0), cache=False)
    double = measure.evaluate(_cantilever(-20.0), cache=False)
    flipped = measure.evaluate(_cantilever(10.0), cache=False)
    same = measure.evaluate(_cantilever(-10.0), cache=False)
    e = measure.compare(base, double)
    assert e["analysisChange"] == pytest.approx(1.0, rel=1e-9) and not e["analysisIdentical"]
    assert e["analysisChangeByCombination"] == {"C1": pytest.approx(1.0, rel=1e-6)}
    assert e["designRatioChange"] is None
    assert measure.compare(base, flipped)["analysisChange"] == pytest.approx(2.0, rel=1e-9)
    z = measure.compare(base, same)
    assert z["analysisIdentical"] is True and z["analysisChange"] == 0.0
    # a 4 % change is not material, a 6 % change is
    for load, material in ((-10.4, False), (-10.6, True)):
        e = measure.compare(base, measure.evaluate(_cantilever(load), cache=False))
        assert e["analysisChange"] == pytest.approx(abs(load + 10.0) / 10.0, rel=1e-6)
        assert measure.is_material("changes-analysis", e)[0] is material


def test_measure_metric_does_not_depend_on_the_roll_of_a_member():
    """A square box rolled by 90 degrees is the same structure: local-axis results swap places, the
    metric (global vectors and resultants) must not move."""
    a = _cantilever()
    a["sections"]["sec"] = {"type": "box", "H": 0.2, "B": 0.2, "t": 0.01}
    b = copy.deepcopy(a)
    b["members"]["M"]["orientation"] = {"roll": 90}
    ea, eb = measure.evaluate(a, cache=False), measure.evaluate(b, cache=False)
    ma, mb = ea.results.combinations["C1"].members["M"], eb.results.combinations["C1"].members["M"]
    assert abs(ma.forces[0, 5]) > 1.0 and abs(mb.forces[0, 5]) < 1e-6 * abs(ma.forces[0, 5])   # Mz became My
    assert measure.compare(ea, eb)["analysisChange"] < 1e-9
    with pytest.raises(cb.NotApplicable, match="precondition not met"):
        cb.apply(a, "member_rolled_90", 0)                 # Iy = Iz: rolling it cannot matter


def test_measure_classifies_outcomes():
    ok = measure.evaluate(_cantilever(), cache=False)
    free = _cantilever()
    free["supports"] = {"A": {"type": "pinned"}}
    unstable = measure.evaluate(free, cache=False)
    assert unstable.status == "unstable" and unstable.diagnosis["nMechanisms"] >= 1
    bad = _cantilever()
    bad["members"]["M"]["section"] = "missing"
    invalid = measure.evaluate(bad, cache=False)
    assert invalid.status == "invalid" and invalid.problems[0]["code"] == "ref.unknown-section"
    assert measure.classify(ok, unstable, None) == "creates-mechanism"
    assert measure.classify(ok, invalid, None) == "invalid-model"
    assert measure.classify(unstable, ok, None) is None
    assert measure.is_material("creates-mechanism", None) == (True, "creates-mechanism")


def test_changes_analysis_through_the_design_ratio_alone():
    """A lighter section on a statically determinate cantilever: the forces do not change (well
    below 5 %), the design ratio does -- material through the design ratio."""
    base = _cantilever(design=True)
    lighter = copy.deepcopy(base)
    lighter["sections"]["sec"]["tf"] = 0.008
    e = measure.compare(measure.evaluate(base, cache=False), measure.evaluate(lighter, cache=False))
    assert not e["analysisIdentical"]
    forces_only = {k: v for k, v in e.items()}
    assert e["designRatioChange"] > 0.05
    assert measure.is_material("changes-analysis", forces_only)[0]
    small = {**e, "analysisChange": 0.0}
    assert measure.is_material("changes-analysis", small)[1] == "designRatioChange"
    assert measure.is_material("changes-analysis", small, ("analysisChange",))[0] is False
    ok, why = measure.is_material("changes-analysis", small, ("analysisChange", "designPassFailChanged"))
    assert ok is False and "does not count for this curve ball" in why
    assert measure.is_material("changes-analysis", {**small, "designPassFailChanged": ["M"]},
                               ("analysisChange", "designPassFailChanged")) == (True, "designPassFailChanged")
    # which criteria count: a design ratio alone only where the member's resistance changes
    full = ("analysisChange", "designPassFailChanged", "designRatioChange")
    for name in cb.available():
        c = cb.get(name)
        if "changes-analysis" in c.effect_classes:
            expected = full if c.category == "members" else (("analysisChange",) if name == "second_order_off"
                                                             else full[:2])
            assert c.criteria_for("changes-analysis") == expected, name
        if c.category == "design":
            assert c.criteria_for("design-only") == ("designPassFailChanged", "designRatioChange")


# ============================================================================ specific curve balls
def test_units_relabelled_keeps_the_numbers_and_changes_the_physics():
    base = model("portal1")
    p = _applied("units_relabelled", "portal1")
    e = p.record.effect
    assert [c.path for c in p.record.changes] in (["/units/force"], ["/units/length"])
    factor = e["physicalFactor"]
    assert factor != 1.0 and e["relabelled"] in ("force", "length")
    (m0, r0), (m1, r1) = run(base), run(p.document)
    d0, d1 = r0.to_document(m0), r1.to_document(m1)
    for cid in d0["combinations"]:                         # numbers in FILE units: the same
        for node, row in d0["combinations"][cid]["displacements"].items():
            assert row == pytest.approx(d1["combinations"][cid]["displacements"][node], rel=1e-6, abs=1e-12)
        for node, row in d0["combinations"][cid]["reactions"].items():
            assert row == pytest.approx(d1["combinations"][cid]["reactions"][node], rel=1e-6, abs=1e-9)
    c0, c1 = r0.combinations["ULS1"], r1.combinations["ULS1"]
    if e["relabelled"] == "force":                         # physical forces scale, physical lengths do not
        assert np.abs(c1.reactions[:, :3]).max() == pytest.approx(factor * np.abs(c0.reactions[:, :3]).max(), rel=1e-6)
        assert np.abs(c1.displacements[:, :3]).max() == pytest.approx(np.abs(c0.displacements[:, :3]).max(), rel=1e-6)
    else:
        assert np.abs(c1.displacements[:, :3]).max() == pytest.approx(factor * np.abs(c0.displacements[:, :3]).max(), rel=1e-6)
    assert e["designRatioChange"] < 1e-6                   # self-similar: no design ratio moves
    assert e["criterion"] == "analysisChange"


def test_mechanism_reference_is_the_structured_diagnosis():
    p = _applied("release_mechanism", "portal1")
    d = p.record.reference["diagnosis"]
    assert d["stable"] is False and d["classification"] == "mechanism"
    node = p.record.entities["nodes"][0]
    assert any(n["node"] == node for mode in d["modes"] for n in mode["nodes"])
    assert any(mode["cause"] == "all-member-ends-released" for mode in d["modes"])
    assert "not a set of results" in p.record.reference["note"]
    released = [c for c in p.record.changes if "/releases/" in c.path or c.path.endswith("/releases")]
    assert released and len(released) == len(p.record.changes)


def test_pinned_frame_is_a_mechanism_on_pinned_bases_and_a_cantilever_frame_on_fixed_ones():
    p = _applied("release_pinned_frame", "frame2d")
    assert p.record.effect_class == "creates-mechanism"
    assert set(cb.get("release_pinned_frame").effect_classes) == {"creates-mechanism", "changes-analysis"}
    q = cb.apply(model("portal1"), "release_pinned_frame", 0)          # this portal has fixed bases
    assert q.record.effect_class == "changes-analysis" and q.record.effect["analysisChange"] > 0.3
    with pytest.raises(cb.NotApplicable):
        cb.apply(model("portal1"), "release_pinned_frame", 0, effect_class="creates-mechanism")


def test_support_curve_balls_leave_the_plane_restraint_alone():
    from benchfea.bench.curveballs._util import plane_dofs, support_dofs
    base = model("portal1")
    plane = plane_dofs(base)
    assert len(plane) == 3                                             # a plane frame in 3D space
    for name in ("support_fixed_to_pinned", "support_removed"):
        p = cb.apply(base, name, 0)
        node = p.record.entities["supports"][0]
        after = support_dofs(p.document["supports"][node], base["conventions"]["verticalAxis"])
        assert all(after[d] == "fixed" for d in plane)
    assert plane_dofs(model("frame3d")) == set()


def test_verification_route_is_recorded_and_degrading_targets_come_last():
    base = model("portalSecondOrder")
    assert tier_c_expected(base)["crossCheckable"] is True
    # rigid offsets with second-order combinations: no two foreign solvers are expected to follow
    p = cb.apply(base, "offset_added", 0)
    v = p.record.verification
    assert v["crossCheckable"] == {"base": True, "perturbed": False}
    assert "not-representable" in v["route"] and v["notes"]
    assert tier_c_expected(p.document)["crossCheckable"] is False
    # a curve ball that stays inside says so
    q = cb.apply(base, "load_sign_flipped", 0)
    assert q.record.verification["crossCheckable"] == {"base": True, "perturbed": True}
    assert "expected to be able to follow" in q.record.verification["route"]
    assert "diagnosis tiers" in _applied("release_mechanism", "portal1").record.verification["route"]
    assert "benchfea validate" in _applied("constraint_wrong_slave", "diaphragm").record.verification["route"]


def test_targets_that_keep_the_model_cross_checkable_are_tried_first(monkeypatch):
    """Two candidates, both material: one adds a rigid offset to a second-order model (no two foreign
    solvers are expected to follow), one scales a load.  Whatever the seed, the second is taken."""
    from benchfea.bench.curveballs import registry

    def candidates(ctx):
        def offset(d):
            d["members"]["B1"]["offsets"] = {"axes": "local", "i": [0.3, 0, 0]}

        def load(d):
            d["loadCases"]["D"]["loads"][1]["w1"] *= 3
        return [registry.Candidate("offset", offset, {"members": ["B1"]}),
                registry.Candidate("load", load, {"members": ["B1"]})]

    spec = registry.CurveBall(name="_test_order", category="loads", title="t", changes="c", why="w", moves="m",
                              effect_classes=("changes-analysis",), precondition="p", candidates=candidates)
    monkeypatch.setitem(registry._REGISTRY, "_test_order", spec)
    base = model("portalSecondOrder")
    seen = set()
    for seed in range(6):
        p = cb.apply(base, "_test_order", seed)
        seen.add(p.record.target)
        assert p.record.rejected == []
    assert seen == {"load"}
    # the offset candidate alone is accepted too (it is material), and is recorded as not cross-checkable
    alone = registry.CurveBall(name="_test_offset", category="offsets", title="t", changes="c", why="w", moves="m",
                               effect_classes=("changes-analysis",), precondition="p",
                               candidates=lambda ctx: candidates(ctx)[:1])
    monkeypatch.setitem(registry._REGISTRY, "_test_offset", alone)
    q = cb.apply(base, "_test_offset", 0)
    assert q.record.target == "offset" and q.record.verification["crossCheckable"] == {"base": True, "perturbed": False}
    # the order of the two candidates in the seeded shuffle does vary with the seed
    orders = set()
    for seed in range(12):
        draw = registry.Draw(registry.stable_seed("curveball", registry.CURVEBALL_VERSION, "_test_order", seed))
        orders.add(tuple(registry._shuffle(["offset", "load"], draw)))
    assert len(orders) == 2


def test_apply_with_verification_accepts_only_a_verified_reference():
    """``verify=True`` runs the verification pipeline on the perturbed model.  Needs both foreign
    solvers (tier C)."""
    from benchfea.crosscheck.runner import available
    if not (available("opensees") and available("pynite")):
        pytest.skip("OpenSeesPy and PyNite are needed for tier C")
    base = model("portalOffsets")
    p = cb.apply(base, "offset_removed", 0, verify=True)
    ran = p.record.verification["run"]
    assert ran["status"] == "verified" and ran["mode"] == "analysis" and ran["tiersNotPassed"] == []
    assert {"C", "D"} <= set(ran["tiersPassed"])
    assert text(p.record.apply(base)) == text(p.document)
    # a reversed offset puts a node on the member axis: one foreign solver left, not verifiable -- refused
    assert _applied("offset_sign_reversed", "portalOffsets").record.verification["crossCheckable"]["perturbed"] is False
    with pytest.raises(cb.NotApplicable, match="with a verified reference") as info:
        cb.apply(base, "offset_sign_reversed", 0, verify=True)
    assert info.value.tried[0]["verification"]["status"] == "quarantined"
    assert info.value.tried[0]["verification"]["category"] == "not-representable"
    # a mechanism: the diagnosis is the reference and goes through the diagnosis tiers
    q = cb.apply(model("portal1"), "release_mechanism", 0, verify=True)
    assert q.record.verification["run"]["kind"] == "diagnosis" and q.record.verification["run"]["status"] == "verified"
    with pytest.raises(ValueError):
        cb.apply(base, "offset_removed", 0, verify="yes")


def test_loaded_memberless_master_is_recorded_as_not_expected_to_verify():
    """QUARANTINED CASE (docs/CURVEBALLS.md section 4): a diaphragm removed leaves its loaded master
    without members; the engine's diagnosis (cause no-stiffness) is not confirmed by the diagnosis
    tiers.  The record must say so instead of promising a verifiable reference."""
    p = _applied("constraint_removed", "diaphragm")
    modes = p.record.reference["diagnosis"]["modes"]
    assert any(m["cause"] == "no-stiffness" and m["kind"] == "node" and not m["members"] for m in modes)
    assert p.record.verification["expectedVerifiable"] is False
    assert "NOT expected to pass" in p.record.verification["route"]
    # a mechanism made of released member ends carries no such warning
    assert "expectedVerifiable" not in _applied("release_mechanism", "portal1").record.verification


def test_constraint_made_flexible_replaces_the_constraint_by_members():
    base = model("diaphragm")
    p = _applied("constraint_made_flexible", "diaphragm")
    cid = p.record.entities["constraints"][0]
    c = base["constraints"][cid]
    assert cid not in p.document.get("constraints", {})
    new = [m for m in p.document["members"] if m not in base["members"]]
    assert sorted(new) == sorted(p.record.entities["members"]) and len(new) == len(c["slaves"])
    assert {(p.document["members"][m]["i"], p.document["members"][m]["j"]) for m in new} == \
        {(c["master"], s) for s in c["slaves"]}
    assert all(p.document["members"][m]["section"] in base["sections"] for m in new)


# =============================================================================================== docs
def test_docs_table_lists_every_curve_ball_and_family():
    path = ROOT / "docs" / "CURVEBALLS.md"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(line for line in lines if line.startswith("| curve ball |"))
    fams = re.findall(r"`([a-z0-9_]+)`", header)
    assert fams == families()
    rows = {}
    for line in lines:
        m = re.match(r"\| `([a-z0-9_]+)` \|(.*)\|$", line)
        if m and m.group(1) in EXPECTED:
            rows[m.group(1)] = [c.strip() for c in m.group(2).split("|")]
    assert set(rows) == set(cb.available())
    for name, cells in rows.items():
        assert len(cells) == len(fams)
        assert any(c != "-" for c in cells), f"{name}: applicable to no family in the table"
        key, cls = EXPECTED[name]
        letter = {"changes-analysis": "A", "design-only": "D", "creates-mechanism": "M", "invalid-model": "I"}[cls]
        assert letter in cells[fams.index(MODELS[key][0])].split(), (name, key)
    body = "\n".join(lines)
    for name in cb.available():
        assert f"#### `{name}`" in body
        assert cb.get(name).why in body
    # nothing in the document claims that anybody checked, approved or certified a result
    for line in lines:
        if re.search(r"certif|approved|reviewed|signed off|checked by", line, re.I):
            assert "Nothing here says that a result has been reviewed or approved by anyone" in line, line
    assert "passed the tiers listed, at the tolerances stated" in body


# ================================================================================================ CLI
def test_cli_list_and_help(capsys):
    assert main(["perturb", "--list"]) == 0
    out = capsys.readouterr().out
    for name in cb.available():
        assert name in out
    assert "creates-mechanism" in out and "materiality" in out
    assert main(["perturb", "--list", "--json"]) == 0
    listed = json.loads(capsys.readouterr().out)
    assert [x["name"] for x in listed] == cb.available() and all(x["whyItMatters"] for x in listed)
    with pytest.raises(SystemExit) as info:
        main(["perturb", "-h"])
    assert info.value.code == 0
    text_ = capsys.readouterr().out
    assert "example:" in text_ and "benchfea perturb" in text_ and "--show-record" in text_ and "--seed" in text_


def test_cli_applies_and_writes_model_and_record(tmp_path, capsys):
    src = tmp_path / "base.json"
    src.write_text(_model_json("portal1"))
    out, rec = tmp_path / "variant.json", tmp_path / "variant.record.json"
    code = main(["perturb", str(src), "support_fixed_to_pinned", "--seed", "0", "-o", str(out), "--record", str(rec)])
    captured = capsys.readouterr()
    assert code == 0 and captured.out == ""
    assert "support_fixed_to_pinned seed 0" in captured.err and "effect class changes-analysis" in captured.err
    record = ChangeRecord.from_dict(json.loads(rec.read_text()))
    variant = json.loads(out.read_text())
    assert text(record.apply(model("portal1"))) == text(variant)
    assert out.read_text() == _applied("support_fixed_to_pinned", "portal1").json()       # same as the API
    assert main(["validate", str(out)]) == 0
    capsys.readouterr()
    # --show-record prints the record; without -o the model is not printed with it
    assert main(["perturb", str(src), "support_fixed_to_pinned", "--show-record"]) == 0
    shown = json.loads(capsys.readouterr().out)
    assert shown == json.loads(rec.read_text())
    # no -o and no --show-record: the perturbed model on standard output
    assert main(["perturb", str(src), "support_fixed_to_pinned"]) == 0
    assert capsys.readouterr().out == out.read_text()


def test_cli_exit_codes(tmp_path, capsys):
    src = tmp_path / "base.json"
    src.write_text(_model_json("portal1"))
    assert main(["perturb", str(src), "offset_removed"]) == 7
    assert "not applicable" in capsys.readouterr().err
    assert main(["perturb", str(src), "offset_removed", "--json"]) == 7
    refusal = json.loads(capsys.readouterr().out)
    assert refusal["ok"] is False and "precondition" in refusal["notApplicable"]
    assert main(["perturb", str(src), "no_such_curve_ball"]) == 2
    assert "unknown curve ball" in capsys.readouterr().err
    assert main(["perturb"]) == 2
    assert main(["perturb", str(src)]) == 2
    assert main(["perturb", str(tmp_path / "missing.json"), "load_moved"]) == 2
    capsys.readouterr()
