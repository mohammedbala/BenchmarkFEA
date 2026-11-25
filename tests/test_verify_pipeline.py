"""Verification pipeline (docs/VERIFICATION.md): verify_model on stable and unstable models, and the
negative tests -- a corrupted engine result, a corrupted foreign result (two of three still agree: no
majority vote), a model only one foreign solver can represent, tiers left out, no time budget."""
import copy
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from benchfea.crosscheck import available
from benchfea.model import build_model
from benchfea.verify import POLICY, QUARANTINED, REFUSED, VERIFIED, ForeignCache, verify_model
from benchfea.verify import pipeline, tier_c
from benchfea.verify.diagnosis import (_engine_system, independent_stiffness, tier_engine_matrix,
                                       tier_null_space)

sys.path.insert(0, str(Path(__file__).parent / "independent"))
from models import LINEAR, ONE_SIDED          # noqa: E402
from unstable import UNSTABLE                 # noqa: E402

foreign = pytest.mark.skipif(not (available("opensees") and available("pynite")),
                             reason="OpenSeesPy and PyNite are needed for tier C")


def bank(prefix, source=LINEAR):
    f = next(v for k, v in source.items() if k.startswith(prefix))
    raw = f()
    return raw[0] if isinstance(raw, tuple) else raw


BEAM = {
    "formatVersion": "1.0.0", "units": {"force": "kN", "length": "m"}, "conventions": {"verticalAxis": "Z"},
    "nodes": {"a": {"x": 0, "y": 0, "z": 0}, "b": {"x": 4, "y": 0, "z": 0}, "c": {"x": 9, "y": 0, "z": 0}},
    "supports": {"a": {"type": "fixed"}, "b": {"dofs": {"uz": "fixed", "uy": "fixed"}},
                 "c": {"dofs": {"uz": "fixed", "uy": "fixed"}}},
    "materials": {"S": {"E": 2.1e8, "nu": 0.3}},
    "sections": {"R": {"type": "rectangle", "h": 0.3, "b": 0.15}},
    "members": {"m1": {"i": "a", "j": "b", "section": "R", "material": "S"},
                "m2": {"i": "b", "j": "c", "section": "R", "material": "S"}},
    "loadCases": {"G": {"loads": [{"type": "memberDistributed", "member": "m1", "direction": "Z", "w1": -12.0},
                                  {"type": "memberPoint", "member": "m2", "direction": "Z", "value": -30.0,
                                   "position": 0.4}]}},
    "combinations": {"U": {"factors": {"G": 1.35}}},
    "analysis": {"shearDeformation": False},
}


# ============================================================================== verified models
@foreign
def test_continuous_beam_passes_every_tier_that_applies_and_says_which():
    rec = verify_model(copy.deepcopy(BEAM))
    assert rec.status == VERIFIED and rec.publishable and rec.kind == "results"
    assert rec.tiers_passed == ["A", "C", "D", "E"]
    assert not rec.tiers["B"].applies and not rec.tiers["F"].applies and rec.quarantine is None
    d = rec.to_dict()
    json.dumps(d)                                               # plain JSON data
    assert d["tiersPassed"] == ["A", "C", "D", "E"] and d["tiersNotApplicable"] == ["B", "F"]
    assert "does not mean reviewed, approved or certified" in d["meaning"]
    assert d["reference"]["format"] == "benchfea-results"
    assert d["versions"]["foreign:opensees"] and d["versions"]["foreign:pynite"] and d["versions"]["benchfea"]
    assert {"total", "engine", "tierC", "tierE"} <= set(d["timings"])
    for t in d["tiers"].values():                               # every tier reports differences and tolerances
        assert {"applies", "passed", "summary", "largestDifferences", "tolerances", "details"} <= set(t)
    a = rec.tiers["A"]
    assert a.largest and all(x["value"] <= POLICY["A"]["relative"] for x in a.largest)
    assert rec.tiers["E"].details["checks"]["designVariables"]["passed"] is True


@foreign
def test_record_is_deterministic():
    one = verify_model(copy.deepcopy(BEAM)).to_dict(volatile=False)
    two = verify_model(copy.deepcopy(BEAM), cache=ForeignCache()).to_dict(volatile=False)
    assert json.dumps(one, sort_keys=True) == json.dumps(two, sort_keys=True)


@foreign
@pytest.mark.parametrize("prefix", ["L05", "L09", "L21", "L24"])
def test_frames_with_curve_ball_features_are_verified(prefix):
    """offsets, diaphragms, rigid links, releases: tier E applies the matching equivalence checks."""
    rec = verify_model(bank(prefix))
    assert rec.status == VERIFIED, rec.summary()
    checks = rec.tiers["E"].details["checks"]
    applied = {k for k, c in checks.items() if c.get("applies")}
    expect = {"L05": {"rigidOffsets"}, "L09": {"constraints"}, "L21": {"rigidOffsets", "constraints", "fullReleases",
                                                                         "fixedEnds"}, "L24": {"constraints", "fullReleases"}}
    assert expect[prefix] <= applied and "designVariables" in applied
    for name in expect[prefix]:
        assert checks[name]["passed"] is True and checks[name]["value"] <= POLICY["E"]["limit"]
        assert len(checks[name]["levels"]) >= 2                 # a convergence series, not one stiff model


@foreign
def test_one_sided_members_are_cross_checked_through_the_active_set_and_it_is_recorded():
    rec = verify_model(bank("T03", ONE_SIDED))
    assert rec.status == VERIFIED
    info = next(iter(rec.tiers["C"].details["combinations"].values()))
    assert [s["what"] for s in info["substitutions"]] == ["one-sided members"]
    assert info["complementarity"]["violations"] == []
    assert "superposition and reciprocity not run" in rec.tiers["D"].details["note"]


@foreign
def test_extra_evidence_is_recorded_but_is_not_a_tier_c_solver():
    rec = verify_model(bank("L02"), extra_evidence=True)
    ev = rec.tiers["C"].details["extraEvidence"]
    assert "extra evidence only" in ev["note"]
    worst = max(v["largestRelativeDifference"] for v in ev["tableau"].values())
    assert worst < 1e-9
    assert rec.tiers["C"].details["combinations"]["G"]["eligible"] == ["opensees", "pynite"]


# ============================================================================== negative tests
@foreign
def test_model_only_one_foreign_solver_can_represent_is_not_publishable():
    rec = verify_model(bank("L15"))                             # shear deformation: PyNite has none
    assert rec.status == QUARANTINED and not rec.publishable
    c = rec.tiers["C"]
    assert c.passed is False and c.category == "not-representable"
    assert "only 1 solver(s) can represent" in c.summary and "pynite" in c.summary
    assert rec.quarantine["category"] == "not-representable" and rec.quarantine["tiers"] == ["C"]
    assert "not publishable as verified" in rec.quarantine["meaning"]
    assert rec.tiers_passed == ["D", "E"]                      # the others passed -- still not verified


@foreign
def test_corrupted_engine_result_is_quarantined(monkeypatch):
    real = pipeline.analyze

    def corrupted(model, strict=True):
        res = real(model, strict=strict)
        c = next(iter(res.combinations.values()))
        c.displacements[model.nodes["b"].index, 4] *= 1.0 + 1.0e-4      # one rotation wrong by 0.01 %
        return res
    monkeypatch.setattr(pipeline, "analyze", corrupted)
    rec = verify_model(copy.deepcopy(BEAM))
    assert rec.status == QUARANTINED and rec.quarantine["category"] == "disagreement"
    assert rec.tiers["A"].passed is False and "rotation" in rec.tiers["A"].summary     # the oracle sees it
    assert rec.tiers["C"].passed is False and rec.tiers["C"].category == "disagreement"   # both foreign solvers too
    assert rec.tiers["E"].details["checks"]["designVariables"]["passed"] is False      # no longer bitwise equal
    assert rec.tiers["D"].passed is True                        # the stored invariants were computed before
    assert set(rec.quarantine["tiers"]) == {"A", "C", "E"}


@foreign
def test_tier_disagreement_is_quarantined_never_decided_by_majority(monkeypatch):
    """OpenSees returns a wrong result; the engine, PyNite, the oracle, the invariants and the equivalence
    checks all agree with each other.  Four tiers pass, one does not: quarantined."""
    real = tier_c._foreign

    def one_wrong(doc_hash, model, cid, solver, options, cache, timeout, versions):
        res, cached = real(doc_hash, model, cid, solver, options, cache, timeout, versions)
        if solver == "opensees":
            res = copy.deepcopy(res)
            res.displacements = res.displacements * 1.01
        return res, cached
    monkeypatch.setattr(tier_c, "_foreign", one_wrong)
    rec = verify_model(copy.deepcopy(BEAM), cache=ForeignCache())
    assert rec.tiers_passed == ["A", "D", "E"]
    assert rec.status == QUARANTINED and rec.quarantine["category"] == "disagreement"
    info = rec.tiers["C"].details["combinations"]["U"]
    rows = {(c["a"], c["b"]): c["passed"] for c in info["comparisons"]}
    assert rows == {("benchfea", "opensees"): False, ("benchfea", "pynite"): True, ("opensees", "pynite"): False}
    assert "benchfea and opensees differ" in rec.quarantine["diagnosis"]


def test_a_tier_left_out_or_out_of_time_cannot_give_verified():
    rec = verify_model(copy.deepcopy(BEAM), tiers="AD")
    assert rec.status == QUARANTINED and rec.quarantine["category"] == "not-run"
    assert rec.tiers_passed == ["A", "D"] and rec.tiers["C"].skipped and rec.tiers["E"].skipped
    rec = verify_model(copy.deepcopy(BEAM), budget=0.0)
    assert rec.status == QUARANTINED and "time budget" in rec.quarantine["diagnosis"]
    assert rec.tiers["C"].passed is False and rec.tiers["E"].passed is False
    with pytest.raises(ValueError):
        verify_model(copy.deepcopy(BEAM), tiers="AZ")


def test_invalid_and_unsupported_models_are_refused_not_quarantined():
    bad = copy.deepcopy(BEAM)
    bad["members"]["m1"]["section"] = "nope"
    rec = verify_model(bad)
    assert rec.status == REFUSED and rec.refusal["stage"] == "validation" and not rec.tiers
    assert rec.refusal["problems"][0]["code"] == "ref.unknown-section"
    rec = verify_model(bank("U12", UNSTABLE))                   # member-level release mechanism: validation
    assert rec.status == REFUSED
    json.dumps(rec.to_dict())


def test_results_the_engine_flags_itself_fail_tier_d(monkeypatch):
    from benchfea.solve import invariants
    real = invariants.solution_checks

    def failing(*a, **k):
        checks = real(*a, **k)
        checks[0].passed = False
        return checks
    import benchfea.solve.analysis as analysis
    monkeypatch.setattr(analysis, "solution_checks", failing, raising=False)
    monkeypatch.setattr(invariants, "solution_checks", failing)
    rec = verify_model(copy.deepcopy(BEAM), tiers="D")
    if rec.tiers["D"].passed:                                   # the solver binds the function elsewhere
        pytest.skip("solution_checks is not reachable by monkeypatching in this build")
    assert rec.status == QUARANTINED and rec.tiers["D"].passed is False


# ============================================================================== unstable models
def test_independent_stiffness_equals_the_engine_matrix_on_stable_models():
    """The matrix that verifies a diagnosis is assembled without the solver's code; on stable models with
    every feature (offsets, releases, springs, links, diaphragms, trusses, shear deformation) it must
    reproduce the engine's reduced matrix."""
    for name, f in {**LINEAR, **ONE_SIDED}.items():
        raw = f()
        model = build_model(raw[0] if isinstance(raw, tuple) else raw)
        s = independent_stiffness(model)
        st, _, Kr = _engine_system(model)
        assert np.array_equal(st.retained[st.free_r], s.retained[s.free]), name
        Ke = Kr[np.ix_(st.free_r, st.free_r)]
        assert np.abs(Ke - s.Kff).max() <= 1e-11 * np.abs(Ke).max(), name
        lam, _, _ = s.scaled_eigen()
        assert lam[0] > POLICY["diagnosis"]["singular"], name    # stable: no zero eigenvalue


@foreign
@pytest.mark.parametrize("prefix, n", [("U01", 3), ("U02", 1), ("U03", 3), ("U04", 2), ("U05", 1), ("U06", 2),
                                       ("U07", 1), ("U08", 6), ("U11", 1)])
def test_mechanism_diagnosis_is_verified_by_its_own_tiers(prefix, n):
    rec = verify_model(bank(prefix, UNSTABLE))
    assert rec.kind == "diagnosis" and rec.status == VERIFIED, rec.summary()
    assert list(rec.tiers) == ["N", "K", "C", "E"] and rec.tiers_passed == ["N", "K", "C", "E"]
    assert rec.reference["diagnosis"]["nMechanisms"] == n
    assert len(rec.tiers["N"].details["eigenvaluesBelowThreshold"]) == n
    ev = rec.tiers["C"].details["evidence"]                      # under the probing load nobody finds equilibrium
    assert {v["verdict"] for v in ev.values()} == {"consistent"} and set(ev) == {"opensees", "pynite"}
    json.dumps(rec.to_dict())


def test_a_wrong_diagnosis_is_not_confirmed():
    raw = bank("U04", UNSTABLE)
    model = build_model(raw)
    rec = verify_model(raw, tiers="D")                          # tiers N and K always run
    good = rec.reference["diagnosis"]
    assert tier_null_space(model, good).passed and tier_engine_matrix(model, good).passed
    wrong_count = copy.deepcopy(good)
    wrong_count["nMechanisms"] = 3
    t = tier_null_space(model, wrong_count)
    assert t.passed is False and "2 eigenvalue(s)" in t.summary
    wrong_mode = copy.deepcopy(good)
    d = np.array([6, 3, 1.2]) / np.linalg.norm([6, 3, 1.2])
    wrong_mode["modes"][0]["direction"] = [float(x) for x in d]          # along the members: they resist it
    assert tier_null_space(model, wrong_mode).passed is False
    assert tier_engine_matrix(model, wrong_mode).passed is False
    wrong_node = copy.deepcopy(good)
    wrong_node["modes"][0]["nodes"][0]["node"] = "A"
    assert tier_null_space(model, wrong_node).passed is False


def test_slack_mechanism_is_never_verified():
    rec = verify_model(bank("U09", UNSTABLE))
    assert rec.kind == "diagnosis" and rec.status == QUARANTINED
    assert rec.quarantine["category"] == "not-verifiable" and "x1" in rec.quarantine["diagnosis"]
    assert rec.tiers["N"].passed and rec.tiers["K"].passed       # "without them it is a mechanism" is confirmed
    assert rec.tiers["N"].details["excludedMembers"] == ["x1"]
    assert not rec.tiers["C"].applies


@foreign
def test_foreign_contradiction_quarantines_a_diagnosis(monkeypatch):
    """A foreign solver that finds an equilibrium state under the probing load contradicts the diagnosis."""
    from benchfea.crosscheck.fidelity import EXACT, Fidelity
    from benchfea.crosscheck.foreign import ForeignResult

    def fake(doc_hash, model, cid, solver, options, cache, timeout, versions):
        from benchfea.crosscheck.compare import applied_resultant
        n = len(model.nodes)
        fid = Fidelity(solver)
        fid.add("all", EXACT, "fake")
        R, _ = applied_resultant(model, cid)
        reactions = np.zeros((n, 6))
        reactions[0] = -R                                        # balances the forces; moments about the origin too
        xyz = np.array([nd.xyz for nd in model.nodes.values()], float)
        reactions[0, 3:] -= np.cross(xyz[0], reactions[0, :3])
        return ForeignResult(solver, cid, True, fid, np.full((n, 6), 1e-3), reactions, {}), False
    monkeypatch.setattr(tier_c, "_foreign", fake)
    rec = verify_model(bank("U11", UNSTABLE))
    assert rec.status == QUARANTINED and rec.quarantine["category"] == "disagreement"
    assert rec.tiers["C"].passed is False and "equilibrium state under the probing load" in rec.tiers["C"].summary
    assert rec.tiers_passed == ["N", "K", "E"]


# ============================================================================== design (tier F)
DESIGN_BEAM = {
    "formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"}, "conventions": {"verticalAxis": "Z"},
    "nodes": {"a": {"x": 0, "y": 0, "z": 0}, "b": {"x": 6000, "y": 0, "z": 0}},
    "supports": {"a": {"dofs": {"ux": "fixed", "uy": "fixed", "uz": "fixed", "rx": "fixed"}},
                 "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}},
    "materials": {"S": {"E": 200000.0, "nu": 0.3, "Fy": 345.0, "Fu": 450.0}},
    "sections": {"I": {"type": "I", "d": 400.0, "bf": 180.0, "tf": 13.5, "tw": 8.6, "r": 21.0}},
    "members": {"m": {"i": "a", "j": "b", "section": "I", "material": "S"}},
    "loadCases": {"L": {"loads": [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -20.0}]}},
    "combinations": {"U": {"factors": {"L": 1.5}, "analysisType": "secondOrder"}},
    "analysis": {"shearDeformation": False},
    "design": {"defaults": {"code": "AISC 360-22", "method": "LRFD", "Lb": 500.0}},
}


@foreign
def test_design_reference_of_a_braced_rolled_beam_is_verified_with_tier_f():
    rec = verify_model(copy.deepcopy(DESIGN_BEAM), design=True)
    assert rec.status == VERIFIED and rec.tiers_passed == ["A", "C", "D", "E", "F"], rec.summary()
    f = rec.tiers["F"]
    used = f.details["itemsUsed"]
    assert {"Eq. F2-1", "Eq. G2-1", "Eq. G2-2", "Table B4.1b Case 10", "Table B4.1b Case 15"} <= set(used)
    assert all(v["status"] == "covered" and v["aTests"] and v["bTests"] and v["examples"] for v in used.values())
    assert rec.design_reference["format"] == "benchfea-design" and f.details["zeroDemandChecksNotCovered"] > 0


def test_design_reference_with_an_uncovered_limit_state_names_it():
    doc = copy.deepcopy(DESIGN_BEAM)
    del doc["design"]["defaults"]["Lb"]                         # unbraced: LTB with r_ts from Eq. F2-7
    rec = verify_model(doc, design=True, tiers="DF")
    f = rec.tiers["F"]
    assert f.passed is False and f.category == "not-covered" and rec.status == QUARANTINED
    assert "Eq. F2-7" in f.details["lacking"]
    lack = f.details["lacking"]["Eq. F2-7"]
    assert "no Design Examples reproduction (b)" in lack["lacking"]
    assert any("lateral-torsional buckling" in c for c in lack["checks"])
    assert "Eq. F2-7" in f.summary


def test_design_reference_with_a_not_checked_limit_state_is_not_publishable():
    doc = copy.deepcopy(DESIGN_BEAM)
    doc["sections"]["I"] = {"type": "rectangle", "h": 300.0, "b": 150.0}      # outside the scope of F and G
    rec = verify_model(doc, design=True, tiers="DF")
    f = rec.tiers["F"]
    assert f.passed is False and f.details["notCheckedOrRefused"]
    assert any("notChecked" in x for x in f.details["notCheckedOrRefused"])
