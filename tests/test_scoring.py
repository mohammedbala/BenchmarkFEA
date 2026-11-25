"""Scoring (benchfea.bench.scoring): tolerances and the partial-credit curve, matching of a diagnosis to
the change record / the verified mechanism diagnosis, repairs re-solved through the verification
pipeline, and robustness to malformed answers.  docs/BENCHMARK.md section 4."""
import copy
import json
import math

import pytest
from bench_helpers import CELLS, case, full_coverage_registry, task

from benchfea.bench import scoring as S
from benchfea.bench import tasks as T


# ------------------------------------------------------------------------------------ credit curve
def test_credit_curve():
    assert S.credit(0.0, 0.01, 0.05) == 1.0 and S.credit(0.01, 0.01, 0.05) == 1.0
    assert S.credit(0.05, 0.01, 0.05) == 0.0 and S.credit(0.3, 0.01, 0.05) == 0.0
    # linear in log(e) between the two: the geometric mean of full and zero earns one half
    assert S.credit(math.sqrt(0.01 * 0.05), 0.01, 0.05) == pytest.approx(0.5)
    assert S.credit(0.02, 0.01, 0.05) == pytest.approx(1 - math.log(2) / math.log(5))
    assert S.credit(float("nan"), 0.01, 0.05) == 0.0 and S.credit(float("inf"), 0.01, 0.05) == 0.0
    xs = [0.001 * 1.2 ** k for k in range(40)]
    ys = [S.credit(x, 0.01, 0.05) for x in xs]
    assert all(a >= b for a, b in zip(ys, ys[1:]))               # monotone


def test_relative_error_uses_the_floor_for_small_references():
    assert S.relative_error(1.02, 1.0, 1e-3) == pytest.approx(0.02)
    assert S.relative_error(-1.0, 1.0, 1e-3) == pytest.approx(2.0)            # the sign counts
    assert S.relative_error(2e-4, 1e-4, 1e-2) == pytest.approx(1e-2)          # |ref| below the floor


def test_tolerance_table():
    for kind in T.QUANTITY_KINDS:
        assert S.tolerance_for("predict", "tool-allowed", kind) == {"full": 0.01, "zero": 0.05}
        assert S.tolerance_for("predict", "no-solver", kind) == {"full": 0.05, "zero": 0.25}
        assert S.tolerance_for("repair", "no-solver", kind) == {"full": 0.02, "zero": 0.10}
        # never tighter than twice the largest tolerance of the verification (tier C: 0.5 %)
        from benchfea.verify import POLICY
        assert S.tolerance_for("predict", "tool-allowed", kind)["full"] >= 2 * POLICY["C"]["maxApproximation"]


# ----------------------------------------------------------------------------------------- predict
def test_predict_scores():
    public, answers = task("change", "predict")
    perfect = S.reference_answer(public, answers)
    out = S.score_answer(public, answers, perfect)
    assert out["score"] == 1.0 and out["valid"] and out["taskType"] == "predict"
    assert all(r["credit"] == 1.0 and r["relativeError"] == 0.0 for r in out["reasons"]["quantities"].values())
    ids = list(perfect["values"])
    n = len(ids)
    # one value 2 % off: partial credit on it, per the curve; the rest full
    a = copy.deepcopy(perfect)
    a["values"][ids[0]] *= 1.02
    out = S.score_answer(public, answers, a)
    part = 1 - math.log(2) / math.log(5)
    assert out["reasons"]["quantities"][ids[0]]["credit"] == pytest.approx(part, abs=1e-9)
    assert out["reasons"]["quantities"][ids[0]]["relativeError"] == pytest.approx(0.02)
    assert out["score"] == pytest.approx((n - 1 + part) / n)
    # wrong sign, a missing value, a text, NaN: no credit for these, reasons given, no crash
    a = copy.deepcopy(perfect)
    a["values"][ids[0]] *= -1
    del a["values"][ids[1]]
    a["values"][ids[2]] = "12 kN"
    out = S.score_answer(public, answers, a)
    q = out["reasons"]["quantities"]
    assert q[ids[0]]["credit"] == 0.0 and q[ids[1]]["error"] == "missing" and q[ids[2]]["error"] == "not a finite number"
    assert out["score"] == pytest.approx((n - 3) / n) and out["valid"]
    a = copy.deepcopy(perfect)
    a["values"][ids[0]] = float("nan")
    a["values"]["q99"] = 1.0
    out = S.score_answer(public, answers, a)
    assert out["reasons"]["quantities"][ids[0]]["credit"] == 0.0 and out["reasons"]["unknownIds"] == ["q99"]
    json.dumps(out)                                              # the record is plain JSON
    # the values of the model BEFORE the curve ball do not pass: the curve ball matters
    base = case("change")["base"]["verification"]["reference"]
    stale = {"values": {x["id"]: T.extract(x, base) for x in public["quantities"]}}
    assert S.score_answer(public, answers, stale)["score"] < 0.2


def test_predict_no_solver_mode_is_judged_more_loosely():
    tool, tool_answers = task("change", "predict", "tool-allowed")
    hand, hand_answers = task("change", "predict", "no-solver")
    a = S.reference_answer(tool, tool_answers)
    a["values"] = {k: v * 1.04 for k, v in a["values"].items()}
    assert S.score_answer(hand, hand_answers, a)["score"] == 1.0
    assert 0.0 < S.score_answer(tool, tool_answers, a)["score"] < 0.2


# ---------------------------------------------------------------------------------------- diagnose
def test_diagnose_matching_to_the_change_record():
    public, answers = task("change", "diagnose")
    ref = answers["reference"]
    perfect = S.reference_answer(public, answers)
    assert S.score_answer(public, answers, perfect)["score"] == 1.0
    w = S.DIAGNOSE_WEIGHTS
    assert sum(w.values()) == pytest.approx(1.0)
    wrong_cat = dict(perfect, category="loads")
    out = S.score_answer(public, answers, wrong_cat)
    assert out["score"] == pytest.approx(w["entities"] + w["direction"]) and out["reasons"]["failed"] == ["category"]
    wrong_dir = dict(perfect, direction="removed")
    assert S.score_answer(public, answers, wrong_dir)["score"] == pytest.approx(w["entities"] + w["category"])
    # a wrong entity: the category alone counts; the direction is not credited without an entity found
    wrong_ent = dict(perfect, entities={"members": ["no-such-member"]})
    out = S.score_answer(public, answers, wrong_ent)
    assert out["score"] == pytest.approx(w["category"]) and set(out["reasons"]["failed"]) == {"entities", "direction"}
    assert out["reasons"]["expected"]["category"] == ref["category"]
    # entities may be given as a list or a single id; the kind is not compared
    ids = [i for v in ref["entities"].values() for i in v]
    assert S.score_answer(public, answers, dict(perfect, entities=ids))["score"] == 1.0
    assert S.score_answer(public, answers, dict(perfect, entities={"anything": ids}))["score"] == 1.0
    # naming extra entities costs precision
    extra = dict(perfect, entities=ids + ["X1", "X2", "X3"])
    assert w["category"] + w["direction"] < S.score_answer(public, answers, extra)["score"] < 1.0
    # a mechanism claimed for a stable model halves the score
    claimed = dict(perfect, mechanism={"kind": "node", "motion": "rotation", "nodes": ["N1"]})
    out = S.score_answer(public, answers, claimed)
    assert out["score"] == pytest.approx(S.FALSE_MECHANISM_FACTOR) and "falseMechanism" in out["reasons"]


def test_entity_match():
    rec = {"nodes": ["N1"], "supports": ["N1"], "members": ["B1", "B2"]}
    assert S.entity_match({"N1"}, rec)["score"] == 1.0                      # every node touched
    assert S.entity_match({"B1", "B2"}, rec)["score"] == 1.0                # every member touched
    half = S.entity_match({"B1"}, rec)
    assert half["precision"] == 1.0 and half["recall"] == 0.5 and half["score"] == pytest.approx(2 / 3)
    assert S.entity_match({"B1", "Z9"}, rec)["precision"] == 0.5
    assert S.entity_match(set(), rec)["score"] == 0.0 and S.entity_match({"Q"}, rec)["score"] == 0.0


def test_diagnose_mechanism_is_matched_to_the_verified_diagnosis():
    public, answers = task("mechanism", "diagnose")
    perfect = S.reference_answer(public, answers)
    assert perfect["mechanism"]["nodes"] and perfect["mechanism"]["kind"] in T.MECHANISM_KINDS
    assert S.score_answer(public, answers, perfect)["score"] == 1.0
    w = S.MECHANISM_WEIGHTS
    assert sum(w.values()) == pytest.approx(1.0)
    # the change named correctly but no mechanism: nothing (the answer of a mechanism case is the mechanism)
    out = S.score_answer(public, answers, dict(perfect, mechanism=None))
    assert out["score"] == 0.0 and out["valid"] and "mechanism" in out["reasons"]["error"]
    m = perfect["mechanism"]
    other_kind = next(k for k in T.MECHANISM_KINDS if k != m["kind"])
    out = S.score_answer(public, answers, dict(perfect, mechanism=dict(m, kind=other_kind)))
    assert out["score"] == pytest.approx(w["nodes"] + w["motion"]) and out["reasons"]["failed"] == ["kind"]
    out = S.score_answer(public, answers, dict(perfect, mechanism=dict(m, nodes=["nowhere"])))
    assert out["score"] == 0.0                                    # kind and motion alone earn nothing
    out = S.score_answer(public, answers, dict(perfect, mechanism=dict(m, nodes=m["nodes"] + ["extra"])))
    assert w["kind"] + w["motion"] < out["score"] < 1.0


# ------------------------------------------------------------------------------------------ repair
def test_repair_is_resolved_through_the_pipeline_and_compared_with_the_intended_design():
    public, answers = task("change", "repair")
    out = S.score_answer(public, answers, S.reference_answer(public, answers))
    assert out["score"] == 1.0 and out["reasons"]["criteria"] == {"stable": True, "verified": True}
    assert out["reasons"]["pipeline"]["status"] == "verified" and out["reasons"]["failed"] == []
    # the model as given (no repair): stable and verified, but every named quantity is off
    out = S.score_answer(public, answers, {"model": public["model"]})
    assert out["score"] == 0.0 and out["reasons"]["criteria"]["verified"]
    assert set(out["reasons"]["failed"]) == set(answers["reference"]["values"])
    assert all(r["relativeError"] > 0.10 for r in out["reasons"]["quantities"].values())


def test_repair_compares_physical_values_whatever_the_units_of_the_submission():
    """A repaired model written in other units is the same structure: its results are converted to the
    units of the intended design before they are compared."""
    from benchfea.unitmap import convert_document
    from benchfea.units import UnitSystem
    public, answers = task("change", "repair")
    base = answers["reference"]["intendedModel"]
    target = {"force": "N", "length": "mm"} if base["units"] != {"force": "N", "length": "mm"} else {"force": "kN", "length": "m"}
    converted = convert_document(copy.deepcopy(base), UnitSystem.from_dict(target))
    assert converted["units"] == target and converted["nodes"] != base["nodes"]
    out = S.score_answer(public, answers, {"model": converted})
    assert out["reasons"]["criteria"]["verified"], out["reasons"]["pipeline"]
    assert out["score"] == pytest.approx(1.0)
    assert max(r["relativeError"] for r in out["reasons"]["quantities"].values()) < 1e-9


def test_repair_of_a_mechanism():
    public, answers = task("mechanism", "repair")
    assert S.score_answer(public, answers, S.reference_answer(public, answers))["score"] == 1.0
    out = S.score_answer(public, answers, {"model": public["model"]})      # still a mechanism
    assert out["score"] == 0.0 and out["reasons"]["criteria"] == {"stable": False, "verified": False}
    assert out["reasons"]["pipeline"]["mechanism"] and out["reasons"]["failed"] == ["stable", "verified"]


def test_repair_with_a_broken_model_scores_zero_with_the_validator_problems():
    public, answers = task("change", "repair")
    broken = copy.deepcopy(answers["reference"]["intendedModel"])
    first = next(iter(broken["members"]))
    broken["members"][first]["section"] = "no-such-section"
    out = S.score_answer(public, answers, {"model": broken})
    assert out["score"] == 0.0 and out["valid"] and not out["reasons"]["criteria"]["stable"]
    assert out["reasons"]["pipeline"]["status"] == "refused"
    assert out["reasons"]["pipeline"]["problems"][0]["code"] == "ref.unknown-section"
    # a repaired model that lost a node / member that is asked for: no credit for those quantities
    public_m, answers_m = task("mechanism", "repair")
    out = S.score_answer(public_m, answers_m, {"model": {"formatVersion": "1.0.0"}})
    assert out["score"] == 0.0 and out["reasons"]["pipeline"]["status"] == "refused"


# ------------------------------------------------------------------------------------------ design
def test_design_scoring_with_the_coverage_test_device(monkeypatch):
    """Scoring of a design task (feasibility gate x weight ratio), exercised with the tier-F test device
    of bench_helpers (with the real registry this cell has no verified design results)."""
    import benchfea.verify.tier_f as tier_f
    monkeypatch.setattr(tier_f, "coverage_registry", full_coverage_registry)
    c = T.build_case(*CELLS["change"])
    public, answers = T.make_task(c, "design", "tool-allowed")
    ref = answers["reference"]
    out = S.score_answer(public, answers, S.reference_answer(public, answers))
    assert out["score"] == 1.0 and out["reasons"]["criteria"] == {"feasible": True, "verified": True}
    cat = public["catalogue"]
    sid, members = next(iter(cat["groups"].items()))
    # the heaviest sections: feasible, score = reference weight / weight < 1
    heavy = {m: cat["options"][sid][-1] for m in members}
    out = S.score_answer(public, answers, {"sections": heavy})
    w = T.model_weight(T.design_document(public["model"], cat, heavy))["value"]
    assert out["score"] == pytest.approx(ref["weight"]["value"] / w) and 0 < out["score"] < 1
    assert out["reasons"]["failed"] == ["weight"] and "not a proven optimum" in out["reasons"]["search"]
    # the lightest sections: not feasible, gate 0
    light = {m: cat["options"][sid][0] for m in members}
    out = S.score_answer(public, answers, {"sections": light})
    assert out["score"] == 0.0 and out["reasons"]["failed"] == ["feasible"] and out["valid"]
    # sections outside the catalogue, unknown members, wrong types: invalid, score 0
    for bad in ({"sections": {members[0]: "W99x999"}}, {"sections": {"nope": cat["options"][sid][0]}},
                {"sections": {members[0]: 7}}, {"sections": [1, 2]}, {}):
        out = S.score_answer(public, answers, bad)
        assert out["score"] == 0.0 and not out["valid"] and out["reasons"]["error"]
    # without the device the same submission does not pass the gate: its design results are not verified
    monkeypatch.undo()
    out = S.score_answer(public, answers, S.reference_answer(public, answers))
    assert out["score"] == 0.0 and out["reasons"]["failed"] == ["verified"]
    assert out["reasons"]["pipeline"]["category"] == "not-covered"


# --------------------------------------------------------------------------------------- malformed
MALFORMED = [None, 0, 1.5, True, "", "not json", "{", "[1, 2]", "null", [], [{"values": {}}], {}, {"values": None},
             {"values": []}, {"values": "q1"}, {"answer": {"values": {}}}, {"values": {"q1": {"deep": [1]}}},
             {"values": {"q1": [1.0]}, "entities": 3, "category": ["x"], "direction": {}, "mechanism": "yes"},
             {"entities": {"members": [["B1"]]}}, {"entities": {"members": {"a": 1}}}, {"entities": 5},
             {"mechanism": {"nodes": 5}}, {"mechanism": {"nodes": [1.5, None]}}, {"mechanism": []},
             {"model": None}, {"model": "model.json"}, {"model": []}, {"model": {}},
             {"model": {"nodes": 1, "members": "x"}}, {"model": {"formatVersion": "1.0.0", "units": []}},
             {"sections": None}, {"sections": {"": ""}}, b"\xff\xfe", "\x00", {"values": {"q1": 1e400}},
             {"category": "x" * 10000, "entities": {"k": ["v"] * 1000}}]


@pytest.mark.parametrize("name,task_type", [("change", "predict"), ("change", "diagnose"), ("mechanism", "diagnose"),
                                            ("change", "repair"), ("mechanism", "repair")])
def test_malformed_answers_score_zero_with_a_reason_and_never_raise(name, task_type):
    public, answers = task(name, task_type)
    for bad in MALFORMED:
        out = S.score_answer(public, answers, bad)
        assert out["score"] == 0.0, bad
        assert 0.0 <= out["score"] <= 1.0 and isinstance(out["reasons"], dict) and out["summary"]
        assert out["taskType"] == task_type
        json.dumps(out)


def test_deeply_nested_answer_does_not_crash():
    public, answers = task("change", "repair")
    deep = {}
    node = deep
    for _ in range(3000):
        node["a"] = {}
        node = node["a"]
    out = S.score_answer(public, answers, {"model": deep})
    assert out["score"] == 0.0
    text = "[" * 100000
    assert S.score_answer(public, answers, text)["score"] == 0.0


def test_every_served_task_is_solvable_by_the_reference_answer():
    for name, types in (("change", ("predict", "diagnose", "repair")), ("mechanism", ("diagnose", "repair"))):
        for tt in types:
            for mode in T.MODES:
                public, answers = task(name, tt, mode)
                assert S.score_answer(public, answers, S.reference_answer(public, answers))["score"] == 1.0
