"""Unstable structures are refused with a structured, deterministic diagnosis (never solved)."""
import json

import numpy as np
import pytest

from benchfea.model import build_model
from benchfea.problems import ModelInvalid, ModelUnstable
from benchfea.solve import MechanismDiagnosis, ModelRefused, analyze

from solve_reference import SEC, make

L = 4000.0
LOAD = {"c": [{"type": "nodal", "node": "b", "Fz": -1000.0}]}


def diagnose(doc) -> MechanismDiagnosis:
    with pytest.raises(ModelUnstable) as e:
        analyze(build_model(doc))
    assert isinstance(e.value, ModelInvalid)
    assert all(p.code == "stability.mechanism" for p in e.value.problems)
    assert len(e.value.problems) == len(e.value.diagnosis.modes)
    json.dumps(e.value.diagnosis.to_dict(), allow_nan=False)
    return e.value.diagnosis


def beam(sa, sb, **member):
    return make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": {"i": "a", "j": "b", **member}}, {"a": sa, "b": sb}, LOAD)


def test_missing_axial_restraint_is_a_rigid_body_translation_in_x():
    d = diagnose(beam({"dofs": {"uy": "fixed", "uz": "fixed", "rx": "fixed"}}, {"dofs": {"uy": "fixed", "uz": "fixed"}}))
    assert d.classification == "mechanism" and d.n_mechanisms == 1 and len(d.modes) == 1
    m = d.modes[0]
    assert (m.kind, m.motion, m.cause) == ("rigid-body", "translation", "no-restraint")
    assert m.description == "rigid-body translation in X: no restraint"
    assert m.direction == (1.0, 0.0, 0.0)
    assert [(n.node, round(n.share, 6)) for n in m.nodes] == [("a", 0.5), ("b", 0.5)]
    assert all(n.dofs == {"ux": 1.0, "uy": 0.0, "uz": 0.0, "rx": 0.0, "ry": 0.0, "rz": 0.0} for n in m.nodes)


def test_only_vertical_supports_leave_three_rigid_body_motions_in_canonical_form():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (0, L, 0)}, {"m": ("a", "b"), "n": ("a", "c"), "o": ("b", "c")},
               {"a": "roller", "b": "roller", "c": "roller"}, LOAD)
    d = diagnose(doc)
    assert d.n_mechanisms == 3
    assert [m.description for m in d.modes] == [
        "rigid-body translation in X: no restraint", "rigid-body translation in Y: no restraint",
        "rigid-body rotation about an axis parallel to Z through node a: no restraint"]
    assert d.modes[2].axis == (0.0, 0.0, 1.0) and d.modes[2].through == (0.0, 0.0, 0.0)


def test_member_free_to_spin_about_its_own_axis():
    doc = make({"a": (0, 0, 0), "b": (0, 0, L)}, {"m": ("a", "b")}, {"a": "pinned", "b": {"dofs": {"ux": "fixed", "uy": "fixed"}}}, LOAD)
    d = diagnose(doc)
    assert d.n_mechanisms == 1
    assert d.modes[0].description == "rigid-body rotation about an axis parallel to Z through node a: no restraint"


def test_rotation_about_a_skew_axis_through_two_pinned_nodes():
    doc = make({"a": (0, 0, 0), "b": (3000, 4000, 0), "c": (3000, 4000, 2000)}, {"m": ("a", "b"), "n": ("b", "c")},
               {"a": "pinned", "b": "pinned"}, {"c": [{"type": "nodal", "node": "c", "Fz": -1.0}]})
    d = diagnose(doc)
    assert d.n_mechanisms == 1 and d.modes[0].motion == "rotation"
    np.testing.assert_allclose(d.modes[0].axis, (0.6, 0.8, 0.0), atol=1e-9)
    assert d.modes[0].description == "rigid-body rotation about the axis (0.6, 0.8, 0) through node a: no restraint"
    assert d.modes[0].nodes[0].node == "c"                      # the node off the axis moves most


def test_node_where_all_member_ends_are_released_spins():
    """Release 'ry' is about member LOCAL y, which is global +Z for a beam along X."""
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (2 * L, 0, 0)},
               {"m": {"i": "a", "j": "b", "releases": {"j": {"ry": True}}}, "n": {"i": "b", "j": "c", "releases": {"i": {"ry": True}}}},
               {"a": "fixed", "c": "fixed"}, LOAD)
    d = diagnose(doc)
    assert d.n_mechanisms == 1
    m = d.modes[0]
    assert m.description == "node b rotation about Z unrestrained: all connected member ends released"
    assert (m.kind, m.motion, m.cause, m.axis, m.members) == ("node", "rotation", "all-member-ends-released", (0.0, 0.0, 1.0), ("m", "n"))
    assert m.nodes[0].to_dict() == {"node": "b", "share": 1.0, "dofs": {"rz": 1.0}}


def test_all_rotations_released_at_a_node_gives_three_node_mechanisms():
    full = {"rx": True, "ry": True, "rz": True}
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (L, L, 0), "d": (L, 0, L)},
               {"m": {"i": "a", "j": "b", "releases": {"j": full}}, "n": {"i": "c", "j": "b", "releases": {"j": full}},
                "o": {"i": "d", "j": "b", "releases": {"j": full}}},
               {"a": "fixed", "c": "fixed", "d": "fixed"}, LOAD)
    d = diagnose(doc)
    assert [m.description for m in d.modes] == [
        f"node b rotation about {ax} unrestrained: all connected member ends released" for ax in "XYZ"]


def test_three_hinges_in_line_are_a_distributed_mechanism_naming_the_nodes_and_members():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (2 * L, 0, 0)},
               {"m": {"i": "a", "j": "b", "releases": {"j": {"rz": True}}}, "n": ("b", "c")},
               {"a": {"type": "pinned", "dofs": {"rx": "fixed"}}, "c": "pinned"}, LOAD)
    d = diagnose(doc)
    assert d.n_mechanisms == 1 and d.modes[0].kind == "distributed" and d.modes[0].count == 1
    m = d.modes[0]
    assert [n.node for n in m.nodes] == ["b", "a", "c"] and m.members == ("m", "n")
    assert set(k for k, v in m.nodes[0].dofs.items() if v > 1e-9) == {"uz", "ry"}     # b drops and rotates
    assert set(k for k, v in m.nodes[1].dofs.items() if v > 1e-9) == {"ry"}           # the ends only rotate
    assert sum(n.share for n in m.nodes) == pytest.approx(1.0)
    assert m.description == "1 mechanism(s) involving several nodes; largest participation: b (uz/ry), a (ry), c (ry)"


def test_portal_with_four_hinges_sways():
    doc = make({"a": (0, 0, 0), "b": (0, 0, 3000), "c": (4000, 0, 3000), "d": (4000, 0, 0)},
               {"c1": {"i": "a", "j": "b", "releases": {"j": {"rz": True}}}, "bm": ("b", "c"),
                "c2": {"i": "d", "j": "c", "releases": {"j": {"rz": True}}}},
               {"a": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}}, "d": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}}},
               LOAD)
    d = diagnose(doc)
    assert d.n_mechanisms == 1 and d.modes[0].kind == "distributed"
    assert {n.node for n in d.modes[0].nodes} == {"a", "b", "c", "d"}
    by = {n.node: n for n in d.modes[0].nodes}
    assert by["b"].dofs["ux"] == pytest.approx(1.0) and by["a"].dofs["ry"] == pytest.approx(1.0)
    assert by["b"].share == pytest.approx(by["c"].share)


def test_collinear_truss_bars_leave_the_middle_node_free_sideways():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (2 * L, 0, 0)},
               {"m": {"i": "a", "j": "b", "type": "truss"}, "n": {"i": "b", "j": "c", "type": "truss"}},
               {"a": "pinned", "c": "pinned"}, {"c": [{"type": "nodal", "node": "b", "Fx": 1.0}]})
    d = diagnose(doc)
    assert [m.description for m in d.modes] == [
        "node b translation in Y unrestrained: the connected truss-type members do not restrain it",
        "node b translation in Z unrestrained: the connected truss-type members do not restrain it"]
    assert all(m.cause == "only-truss-members" and m.kind == "node" for m in d.modes)


def test_disconnected_part_has_six_rigid_body_motions_and_is_named():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (0, L, 0), "d": (L, L, 0)}, {"m": ("a", "b"), "n": ("c", "d")},
               {"a": "fixed"}, LOAD)
    d = diagnose(doc)
    assert d.n_mechanisms == 6 and len(d.modes) == 6
    assert all(m.kind == "rigid-body" and "of the part containing node c" in m.description for m in d.modes)
    assert all({n.node for n in m.nodes} <= {"c", "d"} for m in d.modes)
    assert [m.motion for m in d.modes] == ["translation"] * 3 + ["rotation"] * 3


def test_rigid_body_and_local_mechanisms_are_reported_together():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (2 * L, 0, 0)},
               {"m": {"i": "a", "j": "b", "releases": {"j": {"ry": True}}}, "n": {"i": "b", "j": "c", "releases": {"i": {"ry": True}}}},
               {"a": {"dofs": {"uy": "fixed", "uz": "fixed", "rx": "fixed", "ry": "fixed", "rz": "fixed"}},
                "c": {"dofs": {"uy": "fixed", "uz": "fixed"}}}, LOAD)
    d = diagnose(doc)
    assert d.n_mechanisms == 2
    assert [m.description for m in d.modes] == [
        "rigid-body translation in X: no restraint",
        "node b rotation about Z unrestrained: all connected member ends released"]


def test_diagnosis_is_deterministic_and_structured():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (0, L, 0), "d": (L, L, 0), "e": (L, L, L)},
               {"m": ("a", "b"), "n": ("a", "c"), "o": ("b", "c"), "p": {"i": "c", "j": "d", "releases": {"j": {"ry": True, "rz": True}}},
                "q": {"i": "d", "j": "e", "releases": {"i": {"rx": True, "ry": True, "rz": True}}}},
               {"a": "roller", "b": "roller", "c": "roller"}, LOAD)
    runs = [diagnose(doc).to_dict() for _ in range(3)]
    assert runs[0] == runs[1] == runs[2]
    d = runs[0]
    assert set(d) == {"stable", "classification", "nMechanisms", "nFreeDofs", "tolerance", "modes"} and d["stable"] is False
    assert set(d["modes"][0]) == {"kind", "description", "motion", "cause", "count", "direction", "axis", "through", "nodes", "members"}
    # the same structure with the nodes listed in another order gives the same modes for the same node ids
    other = make({"e": (L, L, L), "d": (L, L, 0), "c": (0, L, 0), "b": (L, 0, 0), "a": (0, 0, 0)},
                 {"m": ("a", "b"), "n": ("a", "c"), "o": ("b", "c"), "p": {"i": "c", "j": "d", "releases": {"j": {"ry": True, "rz": True}}},
                  "q": {"i": "d", "j": "e", "releases": {"i": {"rx": True, "ry": True, "rz": True}}}},
                 {"a": "roller", "b": "roller", "c": "roller"}, LOAD)
    d2 = diagnose(other).to_dict()
    assert d2["nMechanisms"] == d["nMechanisms"]
    assert sorted(m["kind"] for m in d2["modes"]) == sorted(m["kind"] for m in d["modes"])


def test_extreme_stiffness_contrast_is_refused_as_ill_conditioned_not_called_a_mechanism():
    """A tie 1e8 times stiffer than the columns: the smallest eigenvalue of the unit-diagonal-scaled matrix
    is about 7e-11, between 1e-12 (zero to round-off = mechanism) and the default pivotTolerance 1e-10.
    (A contrast that brings it below 1e-12 cannot be told from a mechanism in double precision and is
    reported as one.)"""
    secs = {"S": {"type": "properties", **SEC}, "huge": {"type": "properties", **{k: v * 1e8 for k, v in SEC.items()}}}
    doc = make({"a": (0, 0, 0), "b": (0, 0, 3000), "c": (2000, 0, 3000), "d": (2000, 0, 0)},
               {"c1": ("a", "b"), "tie": {"i": "b", "j": "c", "section": "huge"}, "c2": ("d", "c")},
               {"a": "fixed", "d": "fixed"}, LOAD, sections=secs)
    with pytest.raises(ModelUnstable) as e:
        analyze(build_model(doc))
    d = e.value.diagnosis
    assert d.classification == "ill-conditioned"
    assert "ill-conditioned" in e.value.problems[0].message and "rigid links" in e.value.problems[0].hint
    assert {n.node for n in d.modes[0].nodes} == {"b", "c"}


def test_stable_structures_with_soft_parts_are_not_refused():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed", "b": {"dofs": {"uz": {"spring": 1e-9}}}}, LOAD)
    res = analyze(build_model(doc))
    assert res.combinations["c"].invariants.passed


def test_pivot_tolerance_is_the_threshold_and_is_reported():
    secs = {"S": {"type": "properties", **SEC}, "stiff": {"type": "properties", **{k: v * 1e7 for k, v in SEC.items()}}}
    def doc(tol):
        return make({"a": (0, 0, 0), "b": (0, 0, 3000), "c": (2000, 0, 3000), "d": (2000, 0, 0)},
                    {"c1": ("a", "b"), "tie": {"i": "b", "j": "c", "section": "stiff"}, "c2": ("d", "c")},
                    {"a": "fixed", "d": "fixed"}, LOAD, sections=secs, analysis={"pivotTolerance": tol})
    res = analyze(build_model(doc(1e-12)))
    pd = next(ch for ch in res.combinations["c"].invariants.checks if ch.name == "positive definiteness")
    assert pd.passed and 1e-12 < pd.value < 1e-6 and pd.tolerance == 1e-12 and "condition estimate" in pd.detail
    with pytest.raises(ModelUnstable) as e:
        analyze(build_model(doc(1e-6)))
    assert e.value.diagnosis.tolerance == 1e-6 and e.value.diagnosis.classification == "ill-conditioned"


def test_unsupported_features_are_refused_listing_every_one(fixtures, monkeypatch):
    """Tension-only members, second order, stiffness reduction, notional loads and buckling are analysed
    now (nonlinear.py, buckling.py): what remains refused is an analysis type without a registered solver."""
    from benchfea.solve import ANALYSIS_TYPES
    model = build_model(fixtures / "valid" / "frame_kip_in.json")
    res = analyze(model)
    assert set(res.combinations) == {"LRFD1", "LRFD2", "LRFD3", "SLS"} and set(res.buckling) == {"LRFD1"}
    assert all(c.invariants.passed and c.converged for c in res.combinations.values())
    monkeypatch.delitem(ANALYSIS_TYPES, "secondOrder")
    with pytest.raises(ModelRefused) as e:
        analyze(model)
    assert set(e.value.codes) == {"analysis.not-implemented"}
    assert [p.path for p in e.value.problems] == ["/combinations/LRFD2/analysisType"]
    assert not isinstance(e.value, ModelUnstable)


def test_large_system_uses_the_estimated_eigenvalue_path_and_still_detects_a_mechanism(monkeypatch):
    """Above EXACT_EIGEN_LIMIT free DOFs the smallest eigenvalue is estimated by inverse iteration."""
    import benchfea.solve.system as system
    monkeypatch.setattr(system, "EXACT_EIGEN_LIMIT", 10)
    n = 12
    nodes = {f"n{k}": (k * 500.0, 0, 0) for k in range(n + 1)}
    members = {f"m{k}": (f"n{k}", f"n{k + 1}") for k in range(n)}
    load = {"c": [{"type": "nodal", "node": f"n{n}", "Fz": -1000.0}]}
    res = analyze(build_model(make(nodes, members, {"n0": "fixed"}, load)))
    pd = next(ch for ch in res.combinations["c"].invariants.checks if ch.name == "positive definiteness")
    assert "estimated" in pd.detail and pd.passed
    monkeypatch.setattr(system, "EXACT_EIGEN_LIMIT", 100000)
    exact = analyze(build_model(make(nodes, members, {"n0": "fixed"}, load)))
    pe = next(ch for ch in exact.combinations["c"].invariants.checks if ch.name == "positive definiteness")
    assert "exact" in pe.detail and pd.value == pytest.approx(pe.value, rel=1e-6)
    np.testing.assert_array_equal(res.combinations["c"].displacements, exact.combinations["c"].displacements)
    monkeypatch.setattr(system, "EXACT_EIGEN_LIMIT", 10)
    members["m5"] = {"i": "n5", "j": "n6", "releases": {"j": {"ry": True}}}
    members["m6"] = {"i": "n6", "j": "n7", "releases": {"i": {"ry": True}}}
    with pytest.raises(ModelUnstable) as e:
        analyze(build_model(make(nodes, members, {"n0": "fixed", f"n{n}": "fixed"}, load)))
    assert e.value.diagnosis.modes[0].description == "node n6 rotation about Z unrestrained: all connected member ends released"
