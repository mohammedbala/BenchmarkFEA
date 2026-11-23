"""Tension-only / compression-only members: active-set iteration per combination.

Frame: pin-jointed rectangle a-b-c-d in the X-Z plane (columns a-b and d-c, beam b-c, all truss
members) with X-bracing a-c and d-b.  Rule under test (docs/THEORY.md section 12): a one-sided member
is slack -- zero stiffness -- when the axial force it would carry has the wrong sign.
"""
import numpy as np
import pytest

from benchfea.model import build_model
from benchfea.problems import ModelUnstable
from benchfea.solve import ModelRefused, analyze

from solve_reference import E, SEC, make

H_FRAME, W_FRAME = 3000.0, 4000.0
DIAG = float(np.hypot(H_FRAME, W_FRAME))
EA = E * SEC["A"]
NODES = {"a": (0, 0, 0), "b": (0, 0, H_FRAME), "c": (W_FRAME, 0, H_FRAME), "d": (W_FRAME, 0, 0)}
SUPPORTS = {"a": "pinned", "d": "pinned", "b": {"dofs": {"uy": "fixed"}}, "c": {"dofs": {"uy": "fixed"}}}


def frame(cases, combos=None, brace="tensionOnly", braces=("ac", "db"), analysis=None, **kw):
    members = {"col1": {"i": "a", "j": "b", "type": "truss"}, "col2": {"i": "d", "j": "c", "type": "truss"},
               "beam": {"i": "b", "j": "c", "type": "truss"}}
    if "ac" in braces:
        members["ac"] = {"i": "a", "j": "c", "type": brace}
    if "db" in braces:
        members["db"] = {"i": "d", "j": "b", "type": brace}
    return make(NODES, members, SUPPORTS, cases, combos, analysis, **kw)


WIND = {"W": [{"type": "nodal", "node": "b", "Fx": 10000.0}]}
REVERSING = {"plus": {"factors": {"W": 1.0}}, "minus": {"factors": {"W": -1.0}}}


def test_x_bracing_under_reversing_load_only_the_stretched_brace_works():
    res = analyze(build_model(frame(WIND, REVERSING)))
    plus, minus = res.combinations["plus"], res.combinations["minus"]
    assert plus.members["ac"].active and not plus.members["db"].active
    assert minus.members["db"].active and not minus.members["ac"].active
    # statics of the pin-jointed panel: the working diagonal carries the whole shear in tension
    assert plus.members["ac"].forces[:, 0] == pytest.approx(10000.0 * DIAG / W_FRAME, rel=1e-12)
    assert minus.members["db"].forces[:, 0] == pytest.approx(10000.0 * DIAG / W_FRAME, rel=1e-12)
    assert np.abs(plus.members["db"].forces).max() == 0.0 and np.abs(plus.members["db"].end_forces).max() == 0.0
    # the same answer as the frame that only has the working diagonal, as an ordinary truss member
    one = analyze(build_model(frame(WIND, {"plus": {"factors": {"W": 1.0}}}, brace="truss", braces=("ac",))))
    assert plus.displacements == pytest.approx(one.combinations["plus"].displacements, rel=1e-12, abs=1e-15)
    assert plus.reactions == pytest.approx(one.combinations["plus"].reactions, rel=1e-12, abs=1e-9)
    other = analyze(build_model(frame(WIND, {"minus": {"factors": {"W": -1.0}}}, brace="truss", braces=("db",))))
    assert minus.displacements == pytest.approx(other.combinations["minus"].displacements, rel=1e-12, abs=1e-15)
    assert plus.invariants.passed and minus.invariants.passed
    rec = plus.solver["activeSet"]
    assert rec["slack"] == ["db"] and rec["iterations"] == 2
    assert rec["history"][0] == {"iteration": 1, "slack": [], "violations": ["db"]}
    assert rec["history"][1] == {"iteration": 2, "slack": ["db"], "violations": []}


def test_with_both_diagonals_as_ordinary_members_the_answer_is_different():
    both = analyze(build_model(frame(WIND, REVERSING, brace="truss"))).combinations["plus"]
    one = analyze(build_model(frame(WIND, REVERSING))).combinations["plus"]
    assert both.members["db"].forces[0, 0] < 0.0                               # the second diagonal would be compressed
    assert one.displacements[1, 0] > 1.5 * both.displacements[1, 0]


def test_compression_only_members_mirror_the_rule():
    res = analyze(build_model(frame(WIND, REVERSING, brace="compressionOnly")))
    plus = res.combinations["plus"]
    assert plus.members["db"].active and not plus.members["ac"].active
    assert plus.members["db"].forces[:, 0] == pytest.approx(-10000.0 * DIAG / W_FRAME, rel=1e-12)


def test_no_superposition_each_combination_is_solved_on_its_own():
    """Two wind cases of opposite sign: their sum is zero load, but the sum of the separately solved
    results is not zero -- so combinations must be (and are) solved directly."""
    cases = {"A": [{"type": "nodal", "node": "b", "Fx": 10000.0}], "B": [{"type": "nodal", "node": "b", "Fx": -4000.0}]}
    combos = {"A": {"factors": {"A": 1.0}}, "B": {"factors": {"B": 1.0}}, "AB": {"factors": {"A": 1.0, "B": 1.0}},
              "none": {"factors": {"A": 1.0, "B": 2.5}}}
    res = analyze(build_model(frame(cases, combos)))
    a, b, ab, none = (res.combinations[k] for k in ("A", "B", "AB", "none"))
    assert ab.members["ac"].forces[0, 0] == pytest.approx(6000.0 * DIAG / W_FRAME, rel=1e-12)
    assert a.members["db"].forces[0, 0] + b.members["db"].forces[0, 0] == pytest.approx(4000.0 * DIAG / W_FRAME, rel=1e-12)
    assert not ab.members["db"].active                                        # superposition would give a force here
    # zero load: nothing is stretched or shortened, both diagonals keep their state, all results are zero
    assert none.members["ac"].active and none.members["db"].active and np.abs(none.displacements).max() == 0.0


def test_gravity_only_leaves_both_diagonals_slack_and_the_frame_a_mechanism():
    cases = {**WIND, "G": [{"type": "nodal", "node": "b", "Fz": -50000.0}, {"type": "nodal", "node": "c", "Fz": -50000.0}]}
    combos = {"wind": {"factors": {"W": 1.0}}, "gravity": {"factors": {"G": 1.0}}}
    model = build_model(frame(cases, combos))
    with pytest.raises(ModelUnstable) as e:
        analyze(model)
    codes = [p.code for p in e.value.problems]
    assert codes[0] == "stability.slack-mechanism" and "stability.mechanism" in codes
    msg = e.value.problems[0].message
    # both diagonals shorten when the columns shorten; the frame can then sway freely between the two
    # positions in which one of them becomes taut: no unique answer
    assert "'gravity'" in msg and "ac" in msg and "db" in msg and "mechanism" in msg and "not determined" in msg
    assert e.value.problems[0].path == "/combinations/gravity"
    assert e.value.diagnosis.n_mechanisms == 1                                # the sway of the unbraced panel
    assert set(e.value.results.combinations) == {"wind"}                      # the other combination was solved
    partial = analyze(model, strict=False)
    assert set(partial.combinations) == {"wind"} and any("NO RESULTS for combination 'gravity'" in w for w in partial.warnings)


def test_gravity_plus_enough_wind_is_stable_and_gravity_reaches_the_working_brace():
    cases = {**WIND, "G": [{"type": "nodal", "node": "b", "Fz": -50000.0}, {"type": "nodal", "node": "c", "Fz": -50000.0}]}
    c = analyze(build_model(frame(cases, {"c": {"factors": {"W": 1.0, "G": 1.0}}}))).combinations["c"]
    assert c.members["ac"].active and not c.members["db"].active
    assert c.members["ac"].forces[0, 0] == pytest.approx(10000.0 * DIAG / W_FRAME, rel=1e-12)
    assert c.reactions[:, 2].sum() == pytest.approx(100000.0, rel=1e-12) and c.invariants.passed


def test_member_loads_of_a_slack_member_still_reach_its_nodes():
    mats = {"M": {"E": E, "G": 77000.0, "unitWeight": 7.7e-5}}
    cases = {**WIND, "SW": [{"type": "selfWeight", "members": ["ac", "db"]}]}
    c = analyze(build_model(frame(cases, {"c": {"factors": {"W": 1.0, "SW": 1.0}}}, materials=mats))).combinations["c"]
    weight = 7.7e-5 * SEC["A"] * DIAG
    assert not c.members["db"].active
    assert c.reactions[:, 2].sum() == pytest.approx(2 * weight, rel=1e-12)     # both diagonals weigh on the supports
    slack = c.members["db"]
    assert slack.forces[0, 0] != 0.0 and abs(slack.forces[len(slack.x) // 2, 5]) > 0.0   # its own weight: simple beam
    assert c.invariants.passed


def test_second_order_combination_with_tension_only_bracing():
    """Gravity on the columns + wind, second order: the same result as the frame with the working
    diagonal as an ordinary truss member (and more sway than first order)."""
    cases = {**WIND, "G": [{"type": "nodal", "node": "b", "Fz": -1.0e6}, {"type": "nodal", "node": "c", "Fz": -1.0e6}]}
    so = {"c": {"factors": {"W": 1.0, "G": 1.0}, "analysisType": "secondOrder"}}
    fo = {"c": {"factors": {"W": 1.0, "G": 1.0}}}
    c = analyze(build_model(frame(cases, so))).combinations["c"]
    ref = analyze(build_model(frame(cases, so, brace="truss", braces=("ac",)))).combinations["c"]
    first = analyze(build_model(frame(cases, fo))).combinations["c"]
    assert not c.members["db"].active and c.solver["activeSet"]["slack"] == ["db"]
    assert c.displacements == pytest.approx(ref.displacements, rel=1e-10, abs=1e-14)
    # panel statics with P-Delta: the diagonal resists H + (sum of gravity) x drift / height
    assert c.displacements[1, 0] > 1.005 * first.displacements[1, 0]
    T = c.members["ac"].forces[0, 0]
    C1, C2 = -c.members["col1"].forces[0, 0], -c.members["col2"].forces[0, 0]
    # horizontal equilibrium of the two top joints on the displaced frame: the compressed columns push
    # in the direction of the sway (C x drift / height); the stretched diagonal, whose chord has turned
    # by t / DIAG (t = displacement of c normal to the diagonal), pulls back along its new direction
    sway = (C1 * c.displacements[1, 0] + C2 * c.displacements[2, 0]) / H_FRAME
    t = (-H_FRAME * c.displacements[2, 0] + W_FRAME * c.displacements[2, 2]) / DIAG
    assert sway > 0.05 * 10000.0 and t < 0.0
    assert T * W_FRAME / DIAG - T * (t / DIAG) * (H_FRAME / DIAG) == pytest.approx(10000.0 + sway, rel=1e-9)
    assert c.reactions[:, 0].sum() == pytest.approx(-10000.0, rel=1e-9)
    assert c.converged and c.invariants.passed and c.solver["secondOrder"]["iterations"] >= 2


def test_cycling_active_set_is_detected_and_refused(monkeypatch):
    """Forced: a violation test that always asks to flip the same brace back and forth."""
    from benchfea.solve.nonlinear import NonlinearSolver
    monkeypatch.setattr(NonlinearSolver, "violations", lambda self, state, slack: {"ac": 1.0})
    with pytest.raises((ModelRefused, ModelUnstable)) as e:
        analyze(build_model(frame(WIND, {"plus": {"factors": {"W": 1.0}}})))
    p = e.value.problems[0]
    assert p.code in ("solve.active-set-cycling", "stability.slack-mechanism")
    assert "'plus'" in p.message and "ac" in p.message


def test_max_iterations_bounds_the_active_set_iteration(monkeypatch):
    from benchfea.solve.nonlinear import NonlinearSolver
    calls = []

    def never_settles(self, state, slack):
        calls.append(1)
        return {"db": 1.0} if len(calls) % 2 else {"ac": 1.0}
    monkeypatch.setattr(NonlinearSolver, "violations", never_settles)
    with pytest.raises((ModelRefused, ModelUnstable)) as e:
        analyze(build_model(frame(WIND, {"plus": {"factors": {"W": 1.0}}}, analysis={"maxIterations": 6})))
    assert len(calls) <= 2 * 6 + 2 and e.value.problems[0].path == "/combinations/plus"


def test_document_marks_slack_members_and_records_the_iteration():
    model = build_model(frame(WIND, REVERSING))
    out = analyze(model).to_document(model)["combinations"]["plus"]
    assert out["members"]["db"]["active"] is False and out["members"]["ac"]["active"] is True
    assert out["solver"]["activeSet"]["slack"] == ["db"]
