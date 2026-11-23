"""Independent review: unstable models must be refused with a diagnosis that names what moves and why;
refusals must carry reasons; the invariants must fire when an imbalance is injected."""
import copy
from unittest import mock

import numpy as np
import pytest

import benchfea.solve.analysis as analysis_module
import benchfea.solve.beamcolumn as beamcolumn
import benchfea.solve.element as element
import benchfea.solve.system as system
from benchfea.problems import ModelInvalid, ModelUnstable
from benchfea.solve import ModelRefused
from harness import case, compare, member, nodes, run_benchfea, worst
from models import LINEAR, SECOND, sec, start
from unstable import UNSTABLE


def refused(raw):
    with pytest.raises(ModelInvalid) as e:
        run_benchfea(raw)
    for p in e.value.problems:                       # every refusal says what, where and what to do
        assert p.code and p.path and p.message and p.hint, p
    return e.value


def modes(raw):
    e = refused(raw)
    assert isinstance(e, ModelUnstable)
    return e, e.diagnosis.to_dict()


def test_missing_supports_are_named_as_rigid_body_motions():
    e, d = modes(UNSTABLE["U01_supports_leave_a_rigid_translation"]())
    assert d["classification"] == "mechanism" and d["nMechanisms"] == 3
    assert [p.message.split(": ", 1)[1] for p in e.problems] == [
        "rigid-body translation in X: no restraint", "rigid-body translation in Y: no restraint",
        "rigid-body rotation about an axis parallel to Z through node A: no restraint"]
    assert all(p.path == "/supports" for p in e.problems)
    e, d = modes(UNSTABLE["U11_beam_free_to_slide_along_its_axis"]())
    assert [m["description"] for m in d["modes"]] == ["rigid-body translation in X: no restraint"]


def test_plane_frame_on_two_pins_rotates_about_the_line_of_its_supports():
    e, d = modes(UNSTABLE["U02_plane_frame_on_two_pins_rotates_about_the_support_line"]())
    (m,) = d["modes"]
    assert (m["kind"], m["motion"], m["cause"]) == ("rigid-body", "rotation", "no-restraint")
    assert np.allclose(np.abs(m["axis"]), [1, 0, 0]) and "through node A" in m["description"]


def test_joint_with_every_member_end_released():
    e, d = modes(UNSTABLE["U03_joint_with_all_member_ends_released"]())
    assert len(d["modes"]) == 3
    for m, axis in zip(d["modes"], "XYZ"):
        assert (m["kind"], m["cause"]) == ("node", "all-member-ends-released")
        assert m["nodes"][0]["node"] == "B" and f"rotation about {axis}" in m["description"]
        assert set(m["members"]) == {"m1", "m2", "m3"}
    assert all(p.path == "/nodes/B" for p in e.problems)


def test_collinear_pin_ended_members_leave_the_two_transverse_translations():
    e, d = modes(UNSTABLE["U04_collinear_pin_ended_members"]())
    axis = np.array([6, 3, 1.2]) / np.linalg.norm([6, 3, 1.2])
    assert len(d["modes"]) == 2
    for m in d["modes"]:
        assert (m["kind"], m["motion"], m["cause"]) == ("node", "translation", "only-truss-members")
        assert m["nodes"][0]["node"] == "B" and abs(np.dot(m["direction"], axis)) < 1e-5


def test_member_spinning_about_its_own_axis_is_named():
    """Regression (review): this used to read 'no member, support or constraint resists it'."""
    e, d = modes(UNSTABLE["U05_member_spins_about_its_own_axis"]())
    (m,) = d["modes"]
    axis = np.array([3, 2, 4]) / np.linalg.norm([3, 2, 4])
    assert (m["kind"], m["motion"], m["cause"]) == ("node", "rotation", "member-spin")
    assert abs(abs(np.dot(m["axis"], axis)) - 1.0) < 1e-6 and m["members"] == ["strut"]
    assert "member 'strut' is free to spin about its own axis" in e.problems[0].message


def test_storey_held_only_by_a_diaphragm_sways():
    e, d = modes(UNSTABLE["U06_storey_held_by_the_diaphragm_only"]())
    (m,) = d["modes"]
    assert m["kind"] == "distributed" and m["count"] == 2 and d["nMechanisms"] == 2      # sway in X and in Y
    tops = {n["node"]: n["dofs"] for n in m["nodes"] if n["node"].startswith("u")}
    assert set(tops) == {"u0", "u1", "u2", "u3"} and all(set(v) == {"ux", "uy"} for v in tops.values())
    assert {"c0", "c1", "c2", "c3"} <= set(m["members"])


def test_pin_jointed_portal_is_a_sway_mechanism():
    e, d = modes(UNSTABLE["U07_pin_jointed_portal_sways"]())
    (m,) = d["modes"]
    assert m["kind"] == "distributed" and m["count"] == 1
    assert {n["node"] for n in m["nodes"]} == {"A", "B", "C", "D"} and set(m["members"]) == {"c1", "g", "c2"}


def test_unsupported_part_is_identified_by_one_of_its_nodes():
    e, d = modes(UNSTABLE["U08_second_part_not_connected_to_the_supported_part"]())
    assert d["nMechanisms"] == 6 and all(m["kind"] == "rigid-body" for m in d["modes"])
    assert all("of the part containing node P" in m["description"] for m in d["modes"])
    assert all({n["node"] for n in m["nodes"]} <= {"P", "Q", "R"} for m in d["modes"])


def test_slack_tension_only_bracing_is_refused_naming_the_slack_members():
    e, d = modes(UNSTABLE["U09_tension_only_bracing_goes_slack"]())
    assert e.problems[0].code == "stability.slack-mechanism" and "member(s) x1 slack" in e.problems[0].message
    assert {n["node"] for n in d["modes"][0]["nodes"]} == {"B", "C"}
    e, d = modes(UNSTABLE["U10_gravity_only_on_tension_only_x_bracing_is_not_determined"]())
    assert e.problems[0].code == "stability.slack-mechanism" and "not determined" in e.problems[0].message


def test_member_level_release_mechanism_is_refused_before_the_solve():
    e = refused(UNSTABLE["U12_member_release_mechanism"]())
    assert [p.code for p in e.problems] == ["member.release-mechanism"] and e.problems[0].path == "/members/m/releases"


def test_rigid_body_motion_held_only_by_negligible_springs_is_named():
    """Regression (review): six springs of 1e-14 kN/m under a frame used to be reported as
    '6 mechanism(s) involving several nodes'."""
    raw, _ = LINEAR["L05_offsets_off_axis_skew"]()
    raw["supports"] = {"A": {"dofs": {k: {"spring": 1e-14} for k in ("ux", "uy", "uz", "rx", "ry", "rz")}}}
    e, d = modes(raw)
    assert d["nMechanisms"] == 6 and len(d["modes"]) == 6
    assert all((m["kind"], m["cause"]) == ("rigid-body", "negligible-springs") for m in d["modes"])
    assert sorted(m["motion"] for m in d["modes"]) == ["rotation"] * 3 + ["translation"] * 3
    assert all("restrained only by support spring(s) that are negligible" in p.message for p in e.problems)


def test_springs_that_do_restrain_are_not_reported_as_negligible():
    """Vertical rollers + one finite spring in X: only Y and the spin about Z through the spring are
    free; the mode held by the spring must not appear."""
    raw = UNSTABLE["U01_supports_leave_a_rigid_translation"]()
    raw["supports"]["A"] = {"type": "roller", "dofs": {"ux": {"spring": 500.0}}}
    e, d = modes(raw)
    assert d["nMechanisms"] == 2 and all(m["cause"] == "no-restraint" for m in d["modes"])


# ----------------------------------------------------------------------------- refusals with reasons
def _l11():
    return LINEAR["L11_support_springs_stiff_and_soft"]()[0]


@pytest.mark.parametrize("build, code", [
    (lambda: (lambda d: (d["loadCases"]["Q"]["loads"].append({"type": "supportDisplacement", "node": "B", "uz": -5.0}), d)[1])(_l11()),
     "load.displacement-dof-not-fixed"),                                  # imposed displacement on a spring DOF
    (lambda: (lambda d: (d["loadCases"]["Q"]["loads"].append({"type": "supportDisplacement", "node": "C", "uz": -5.0}), d)[1])(_l11()),
     "load.displacement-dof-not-fixed"),                                  # ... on a free DOF
    (lambda: (lambda d: (d["supports"].update({"B": {"dofs": {"ux": {"spring": 40.0}}}}), d)[1])(
        LINEAR["L08_diaphragm_z_up_unequal_columns"]()[0]), "constraint.slave-restrained"),
    (lambda: (lambda d: (d["constraints"].update({"l2": {"type": "rigidLink", "master": "S1", "slaves": ["T"]}}), d)[1])(
        LINEAR["L07_rigid_links_with_lever_arms"]()[0]), "constraint.chained"),
    (lambda: (lambda d: (d["constraints"].update({"l2": {"type": "rigidLink", "master": "C", "slaves": ["S2"]}}), d)[1])(
        LINEAR["L07_rigid_links_with_lever_arms"]()[0]), "constraint.slave-two-masters"),
    (lambda: (lambda d: (d["loadCases"]["W"]["loads"].append({"type": "nodal", "node": "T", "Mx": 5.0}), d)[1])(
        LINEAR["L16_space_truss_with_member_loads"]()[0]), "stability.mechanism"),          # moment on a truss-only node
    (lambda: (lambda d: (d["nodes"]["M"].update({"z": 5.0}), d)[1])(
        LINEAR["L08_diaphragm_z_up_unequal_columns"]()[0]), "solve.invariant-failed"),      # diaphragm master off its plane
    (lambda: (lambda d: (d["members"]["g1"].update({"releases": {"i": {"uy": {"spring": 500.0}}}}), d)[1])(
        SECOND["S03_portal_with_leaning_column"]()[0]), "analysis.release-second-order"),
    (lambda: (lambda d: (d["analysis"].update({"maxIterations": 2}), d)[1])(
        SECOND["S03_portal_with_leaning_column"]()[0]), "solve.not-converged"),
])
def test_refusals_carry_a_code_a_location_and_a_reason(build, code):
    e = refused(build())
    assert code in [p.code for p in e.problems], [p.code for p in e.problems]


def test_diaphragm_master_off_its_plane_is_explained():
    raw, _ = LINEAR["L08_diaphragm_z_up_unequal_columns"]()
    raw["nodes"]["M"]["z"] = 5.0
    e = refused(raw)
    assert "global moment equilibrium" in e.problems[0].message and "plane of the master" in e.problems[0].hint


# ----------------------------------------------------------------------------- invariants fire
def _failed(raw):
    with pytest.raises(ModelRefused) as e:
        run_benchfea(raw)
    assert all(p.code == "solve.invariant-failed" for p in e.value.problems), [p.code for p in e.value.problems]
    return {p.message.split("'")[1] for p in e.value.problems}


EPS = 1.0e-4            # relative size of the injected error (the equilibrium tolerance is 1e-6)


@pytest.mark.parametrize("where", ["linear", "constraints", "second order"])
def test_injected_errors_are_caught_by_the_invariants(where):
    raw = {"linear": LINEAR["L05_offsets_off_axis_skew"], "constraints": LINEAR["L08_diaphragm_z_up_unequal_columns"],
           "second order": SECOND["S03_portal_with_leaning_column"]}[where]()[0]
    run_benchfea(raw)                                             # clean: passes
    solve = system.Structure.solve

    def wrong_reactions(self, *a, **k):
        s = solve(self, *a, **k)
        s.reactions *= 1.0 + EPS
        return s

    def wrong_displacements(self, *a, **k):
        s = solve(self, *a, **k)
        s.u *= 1.0 + EPS
        return s
    with mock.patch.object(system.Structure, "solve", wrong_reactions):
        assert {"global force equilibrium", "joint equilibrium"} <= _failed(raw)
    with mock.patch.object(system.Structure, "solve", wrong_displacements):
        assert {"joint equilibrium", "work = strain energy"} <= _failed(raw)
    if where == "second order":                                   # member loads: fixed-end forces of the pieces
        loads = beamcolumn.BeamColumn.equivalent_nodal_loads
        with mock.patch.object(beamcolumn.BeamColumn, "equivalent_nodal_loads",
                               lambda self, *a, **k: loads(self, *a, **k) * (1.0 - EPS)):
            assert {"global force equilibrium", "joint equilibrium"} <= _failed(raw)
    else:
        loads = element.MemberElement.equivalent_nodal_loads
        with mock.patch.object(element.MemberElement, "equivalent_nodal_loads",
                               lambda self, ld: loads(self, ld) * (1.0 - EPS)):
            assert {"global force equilibrium", "joint equilibrium"} <= _failed(raw)
        recover = element.MemberElement.recover

        def wrong_end_forces(self, *a, **k):
            r = recover(self, *a, **k)
            r.end_forces *= 1.0 + EPS
            return r
        with mock.patch.object(element.MemberElement, "recover", wrong_end_forces):
            assert "joint equilibrium" in _failed(raw)


def test_unsymmetric_member_stiffness_is_caught():
    raw, _ = LINEAR["L05_offsets_off_axis_skew"]()
    build = element.build_element

    def unsymmetric(*a, **k):
        e = build(*a, **k)
        e.k_global[0, 7] *= 1.0 + EPS
        return e
    with mock.patch.object(analysis_module, "build_element", unsymmetric):
        assert "stiffness symmetry" in _failed(raw)


def test_a_wrong_station_value_is_not_seen_by_the_invariants_but_by_the_reference():
    """Limit found by the review: the invariants close the member integration at x = L, they do not
    re-derive the values reported at interior stations.  Those are covered by the independent
    reference (every station of every member in tests/independent)."""
    raw, cid = LINEAR["L05_offsets_off_axis_skew"]()
    recover = element.MemberElement.recover

    def wrong_station(self, *a, **k):
        r = recover(self, *a, **k)
        r.forces[len(r.forces) // 2] *= 1.0 + 1.0e-3
        return r
    with mock.patch.object(element.MemberElement, "recover", wrong_station):
        _, res, doc = run_benchfea(raw)                           # not refused
    assert res.combinations[cid].invariants.passed
    assert worst(compare(raw, cid, doc=doc)) > 1.0e-4             # ... but the reference sees it
