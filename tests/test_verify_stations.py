"""Regression tests of the audit of the verification layer.

* Tier D recomputes the values at the result stations from the finished reference solution
  (``benchfea.verify.stations``): a station value that is wrong in the results -- consistently, in every
  solve of the engine, so that neither the bitwise re-solve nor the member split sees it, and with
  every node value right, so that the foreign solvers cannot see it -- must fail the tier.
* The tier-F registry credits a hand calculation (a) only for checks the hand test's own code read.
* The report says what the passed problems compare against.
* What the limits notes of the AISC Commentary problems with shear deformation state is true.
"""
import copy
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from benchfea.crosscheck import available
from benchfea.model import build_model
from benchfea.solve import analyze
from benchfea.verify import POLICY, QUARANTINED, VERIFIED, ForeignCache, verify_model
from benchfea.verify import common, pipeline
from benchfea.verify.problems import evaluate_targets, load_problem
from benchfea.verify.stations import station_deflections, station_statics
from benchfea.verify.tier_d import tier_d

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent / "independent"))
from models import LINEAR, ONE_SIDED, SECOND          # noqa: E402

foreign = pytest.mark.skipif(not (available("opensees") and available("pynite")),
                             reason="OpenSeesPy and PyNite are needed for tier C")

PORTAL = {
    "formatVersion": "1.0.0", "units": {"force": "kN", "length": "m"}, "conventions": {"verticalAxis": "Z"},
    "nodes": {"a": {"x": 0, "y": 0, "z": 0}, "b": {"x": 0, "y": 0, "z": 4}, "c": {"x": 6, "y": 0, "z": 4},
              "d": {"x": 6, "y": 0, "z": 0}},
    "supports": {"a": {"type": "fixed"}, "d": {"type": "pinned"}},
    "materials": {"S": {"E": 2.1e8, "nu": 0.3}},
    "sections": {"R": {"type": "rectangle", "h": 0.3, "b": 0.15}},
    "members": {"c1": {"i": "a", "j": "b", "section": "R", "material": "S"},
                "bm": {"i": "b", "j": "c", "section": "R", "material": "S"},
                "c2": {"i": "d", "j": "c", "section": "R", "material": "S"}},
    "loadCases": {"G": {"loads": [{"type": "memberDistributed", "member": "bm", "direction": "Z", "w1": -12.0},
                                  {"type": "nodal", "node": "b", "Fx": 10.0}]}},
    "combinations": {"U": {"factors": {"G": 1.35}}},
    "analysis": {"shearDeformation": False, "resultStations": 5},
}


def _solve(raw):
    raw = raw[0] if isinstance(raw, tuple) else raw
    model = build_model(copy.deepcopy(raw))
    return model, analyze(model, strict=True)


def _bank(prefix, source):
    return next(v for k, v in source.items() if k.startswith(prefix))()


# ============================================================================== the checks pass
@pytest.mark.parametrize("prefix, source", [
    ("L01", LINEAR), ("L04", LINEAR), ("L05", LINEAR), ("L07", LINEAR), ("L13", LINEAR), ("L15", LINEAR),
    ("L16", LINEAR), ("L17", LINEAR), ("L18", LINEAR), ("S01", SECOND), ("S02", SECOND), ("S04", SECOND),
    ("S05", SECOND), ("S06", SECOND), ("S07", SECOND), ("T01", ONE_SIDED), ("T03", ONE_SIDED)])
def test_station_values_of_the_model_bank_are_reproduced(prefix, source):
    """Statics at every station to round-off, integrated member equations to 1e-9 -- first and second order,
    shear deformation, releases, offsets, constraints, point moments, trapezoid loads, one-sided members."""
    model, results = _solve(_bank(prefix, source))
    st = station_statics(model, results, POLICY["D"]["stationStatics"], POLICY["D"]["stationStaticsIterated"])
    assert st["passed"], st["failed"][:3] + st["notEvaluated"][:3]
    assert st["largest"].get("value", 0.0) < 1e-10 and st["stations"] > 0
    sd = station_deflections(model, results, POLICY["D"]["stationIntegration"])
    assert sd["passed"], sd["failed"][:3] + sd["notEvaluated"][:3]
    assert sd["largest"].get("value", 0.0) < 1e-9


def test_axial_line_load_in_second_order_is_left_to_the_integration_and_passes():
    """Columns under their own weight, second order: the statics check cannot take the lever of the axial
    line load from station values alone, says so (``deferred``), and the integration covers the moments."""
    model, results = _solve(_bank("S07", SECOND))
    st = station_statics(model, results, 1e-9, 1e-6)
    assert st["deferred"] and "axial line load" in st["deferred"][0] and st["passed"]
    assert tier_d(model, results).passed is True


def _three_stations_at_one_point(results):
    return [(cid, mid, k) for cid, c in results.combinations.items() for mid, mr in c.members.items()
            for k in range(len(mr.x) - 2) if mr.x[k + 2] - mr.x[k] < 1e-9 * mr.x[-1]]


def test_a_regular_station_on_a_point_load_left_left_right():
    """Found by the audit: a regular station one ulp BELOW a point load gives three stations at one point
    holding (left, left, right) limits.  Ziemian frame 9, case 'right', has such members."""
    p = ROOT / "verification" / "problems" / "ziemian" / "frame-9-first-order-displacements-right" / "model.json"
    if not p.exists():
        pytest.skip("verification library not present")
    model, results = _solve(json.loads(p.read_text()))
    cid, mid, k = _three_stations_at_one_point(results)[0]
    v = results.combinations[cid].members[mid].forces[k:k + 3, 1]
    assert abs(v[0] - v[1]) < 1e-9 * abs(v[0]) and abs(v[1] - v[2]) > 1e-3 * abs(v[0])
    assert tier_d(model, results).passed is True


def test_a_regular_station_on_a_point_load_left_right_right_and_a_missing_limit():
    """... and one ulp ABOVE it gives (left, right, right): which limit a station at a point load holds is
    decided by round-off of its position, so the check accepts either -- and requires that the stations
    of the point show both.  (This generated beam was quarantined by the first version of the check.)"""
    from benchfea.bench.curveballs.record import plain_copy
    from benchfea.bench.generate import generate
    model, results = _solve(plain_copy(generate("continuous_beam", 435514608, 1).document))
    cid, mid, k = _three_stations_at_one_point(results)[0]
    v = results.combinations[cid].members[mid].forces[k:k + 3, 1]
    assert abs(v[1] - v[2]) < 1e-9 * abs(v[1]) and abs(v[0] - v[1]) > 1e-3 * abs(v[0])
    assert tier_d(model, results).passed is True
    for c in results.combinations.values():                   # flatten the jump: no limit from the right
        c.members[mid].forces[k + 1:k + 3] = c.members[mid].forces[k].copy()
    t = tier_d(model, results)
    assert t.passed is False and "do not show both the limit from the left and the limit from the right" in " ".join(t.details["failed"])


# ============================================================================== the checks fail
def _corrupt(member, array, station, column, factor=1.01):
    def apply(model, results):
        for c in results.combinations.values():
            if member in c.members:
                a = getattr(c.members[member], array)
                assert abs(a[station, column]) > 1e-9
                a[station, column] *= factor
    return apply


def _everywhere(monkeypatch, corrupt):
    """The same error in EVERY solve of the engine (reference, bitwise re-solve, derived models)."""
    real = pipeline.analyze

    def bad(model, strict=True):
        res = real(model, strict=strict)
        corrupt(model, res)
        return res
    monkeypatch.setattr(pipeline, "analyze", bad)
    monkeypatch.setattr(common, "analyze", bad)


@pytest.mark.parametrize("array, station, column, what", [
    ("forces", 2, 5, "Mz"),              # the span moment at midspan
    ("forces", 3, 5, "Mz"),              # ... at a quarter point (not the station the member split cuts at)
    ("forces", 3, 1, "Vy"),
    ("displacements", 3, 1, "deflection uy"),
    ("displacements", 2, 1, "deflection uy")])
def test_a_wrong_station_value_fails_tier_d(monkeypatch, array, station, column, what):
    model, results = _solve(PORTAL)
    assert tier_d(model, results).passed is True
    _corrupt("bm", array, station, column)(model, results)
    t = tier_d(model, results)
    assert t.passed is False and "station" in t.summary and what in " ".join(t.details["failed"])


@foreign
@pytest.mark.parametrize("array, station, column", [("forces", 3, 5), ("displacements", 3, 1)])
def test_a_consistently_wrong_station_value_is_never_verified(monkeypatch, array, station, column):
    """Before the station checks existed this model came out VERIFIED: tier C compares node values and end
    forces (right), the bitwise re-solve repeats the error, the member split cuts at another station."""
    assert verify_model(copy.deepcopy(PORTAL), cache=ForeignCache()).status == VERIFIED
    _everywhere(monkeypatch, _corrupt("bm", array, station, column))
    rec = verify_model(copy.deepcopy(PORTAL), cache=ForeignCache())
    assert rec.status == QUARANTINED and "D" in rec.quarantine["tiers"]
    assert rec.tiers["C"].passed is True                      # the foreign solvers cannot see a station value
    assert rec.quarantine["category"] == "disagreement"


def test_station_ends_must_equal_the_end_forces_and_the_node_displacements():
    model, results = _solve(PORTAL)
    c = next(iter(results.combinations.values()))
    c.members["bm"].forces[0, 5] *= 1.001                     # station at x = 0 no longer -(end force at i)
    st = station_statics(model, results, 1e-9, 1e-6)
    assert not st["passed"] and "end forces" in st["failed"][0]
    model, results = _solve(PORTAL)
    c = next(iter(results.combinations.values()))
    c.members["bm"].displacements[-1, 1] *= 1.001             # member end no longer where node c is
    sd = station_deflections(model, results, 1e-6)
    assert not sd["passed"] and any("node c" in f for f in sd["failed"])


def test_a_wrong_station_value_in_a_second_order_solution_fails():
    raw = _bank("S03", SECOND)
    model, results = _solve(raw)
    assert tier_d(model, results).passed is True
    c = next(iter(results.combinations.values()))
    # the largest bending moment at an interior station of the model
    mid, k, col = max(((m, k, col) for m, mr in c.members.items() for k in range(1, mr.forces.shape[0] - 1)
                       for col in (4, 5)), key=lambda t: abs(c.members[t[0]].forces[t[1], t[2]]))
    c.members[mid].forces[k, col] *= 1.0 + 1e-4
    t = tier_d(model, results)
    assert t.passed is False and "station" in t.summary


# ============================================================================== tier F registry
def test_hand_calculation_credit_needs_the_test_to_look_at_the_check():
    """(a) used to be credited to every check the engine produced inside a hand test.  Eq. G2-1 (web shear)
    was credited to a deflection test, Eq. D2-1 (tension yielding) to a compression slenderness test."""
    reg = json.loads((ROOT / "src" / "benchfea" / "verify" / "design_coverage.json").read_text())
    items = reg["items"]
    assert not any("test_design_service" in t for t in items["Eq. G2-1"]["a"])
    assert not any("test_hand_e2_slenderness_is_advisory" in t for t in items["Eq. D2-1"]["a"])
    assert any("test_hand_tensile_yielding_d2a" in t for t in items["Eq. D2-1"]["a"])
    assert any("test_design_shear" in t for t in items["Eq. G2-1"]["a"])
    # function-level hand calculation, listed explicitly
    assert reg["functionLevelA"] == [{"item": "Eq. F1-1",
                                      "test": "tests/test_design_flexure.py::test_hand_f1_cb_values_of_the_user_note",
                                      "function": "benchfea.design.aisc360_22.f_flexure.cb_f1_1"}]
    assert items["Eq. F1-1"]["a"] == ["tests/test_design_flexure.py::test_hand_f1_cb_values_of_the_user_note"]
    # a hand test that reads only an intermediate value of a check is not (a) coverage
    assert items["Eq. F7-5"]["a"] == [] and items["Eq. F7-5"]["aDetail"]
    for e in items.values():
        assert not set(e["a"]) & set(e["aDetail"])
        assert all("::test_hand_" in t for t in e["a"] + e["aDetail"])


# ============================================================================== report and library
def _report_tool():
    spec = importlib.util.spec_from_file_location("build_verification_report", ROOT / "tools" / "build_verification_report.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_report_says_what_the_passed_problems_compare_against():
    from benchfea.verify import verify_library
    lib = ROOT / "verification" / "problems"
    if not lib.is_dir():
        pytest.skip("verification library not present")
    tool = _report_tool()
    rep = verify_library(only=["staad/portal-fixed-uniform-load", "aisc-commentary/case-2-us-without-shear-p0",
                               "derived/beam-point-load", "mgz/ex-5-9"], tiers="D")
    by_id = {r.id: r for r in rep.problems}
    assert tool._target_kinds(by_id["portal-fixed-uniform-load"]) == {"printed": 0, "byConstruction": 0,
                                                                      "recomputed": 5, "derived": 0}
    assert tool._target_kinds(by_id["case-2-us-without-shear-p0"]) == {"printed": 1, "byConstruction": 1,
                                                                       "recomputed": 0, "derived": 0}
    assert tool._target_kinds(by_id["ex-5-9-portal-axially-rigid"])["printed"] == 6
    text = "\n".join(tool.evidence_section(rep))
    assert "not always a reproduction of a published number" in text
    assert "| staad | recomputed-from-published-statement | 1 | 1 | 0 | 0 | 0 | 0 | 5 | 0 |" in text
    assert "1 passed problem(s) with a published statement have no printed target at all" in text


def test_aisc_commentary_shear_problems_say_that_they_do_not_test_shear_deformation_and_it_is_true():
    """The source's tolerance (3 % / 5 %) is wider than the whole shear effect of its tables: every
    with-shear target is also met by the same model WITHOUT shear deformation.  The limits note says so."""
    lib = ROOT / "verification" / "problems" / "aisc-commentary"
    if not lib.is_dir():
        pytest.skip("verification library not present")
    dirs = sorted(d for d in lib.iterdir() if d.is_dir() and "-with-shear-" in d.name)
    assert len(dirs) == 16
    for d in dirs:
        spec = load_problem(d)
        assert "does NOT test shear deformation" in spec.limits_note
        raw = json.loads(spec.model_path.read_text())
        raw["analysis"]["shearDeformation"] = False
        model, results = _solve(raw)
        rows = evaluate_targets(spec, results.to_document(model))
        assert all(r.passed for r in rows), d.name


# ============================================================================== dataset audit tool
def test_dataset_audit_tool_finds_nothing_on_the_first_models_of_the_dataset(capsys):
    """tools/audit_dataset.py: split rule, hashes, engine re-solve, independent tableau solver, mechanism
    count, scoring of right and wrong answers, leak probe -- here on the first models only."""
    if not (ROOT / "data" / "tasks" / "manifest.json").exists():
        pytest.skip("dataset not present")
    spec = importlib.util.spec_from_file_location("audit_dataset", ROOT / "tools" / "audit_dataset.py")
    tool = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tool)
    code = tool.main([str(ROOT / "data"), "--limit", "3"])
    out = json.loads(capsys.readouterr().out)
    assert code == 0 and out["findings"] == []
    assert out["stableModels"]["models"] == 3 and out["stableModels"]["largestDifferenceEngineTableau"] < 1e-9
    assert out["stableModels"]["largestDifferenceReferenceFreshSolve"] == 0.0
    assert out["mechanisms"]["countConfirmedByTheTableau"] == out["mechanisms"]["models"] == 3
    assert out["split"]["curveBalledModelsInBothSplits"] == 0
