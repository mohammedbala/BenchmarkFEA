"""The error bound of the foreign solvers in second order (docs/CROSSCHECK.md 1.3 and 3.4) and the
gaps of the verification pipeline closed with it.

* ``estimate_from_levels``: the rule itself, on sequences with a known limit;
* acceptance: on models with an answer that depends neither on the engine nor on a foreign solver (ODE
  oracle, tableau solver of tests/independent) the estimated error of each foreign solver is not below
  its true error, and not vacuous;
* a corrupted engine result (1e-3 on a displacement) does not pass tier C;
* the measured facts the translators rely on (OpenSees: axial force of the PDelta transformation with
  an axial element load; PyNite: the terms of its geometric stiffness);
* rigid constraints in second order, a truss member with an axial member load, a slack member with
  member loads, a loaded node without members, a member split next to a point load.
"""
import copy
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

import crosscheck_acceptance as acceptance
from crosscheck_helpers import E, EA, EIz, G, HAVE, SEC, document, needs_both, needs_opensees, needs_pynite
from benchfea.crosscheck import relative_difference, run
from benchfea.crosscheck.runner import ROUNDOFF_FLOOR, estimate_from_levels
from benchfea.model import build_model
from benchfea.solve import analyze
from benchfea.verify import verify_model
from benchfea.verify.policy import POLICY
from benchfea.verify.tier_c import ForeignCache, slack_member_loads, tier_c

ROOT = Path(__file__).resolve().parents[1]


def build(doc):
    return build_model(copy.deepcopy(doc))


# ====================================================================================== the rule
def _levels(error, n=8):
    """(d1, d2, true error of the finest) of a scalar sequence with the given error function of n."""
    e = [error(n), error(2 * n), error(4 * n)]
    return abs(e[0] - e[1]), abs(e[1] - e[2]), abs(e[2])


@pytest.mark.parametrize("p, factor", [(1.0, 1.0), (2.0, 3.0)])
def test_bound_on_a_pure_power_law_is_the_last_difference(p, factor):
    d1, d2, true = _levels(lambda n: 0.3 * n ** -p)
    est = estimate_from_levels(d1, d2, 2.0, (1.0, 2.0))
    assert est["status"] == "asymptotic" and est["observedOrder"] == pytest.approx(p)
    assert est["estimate"] == pytest.approx(d2) and est["estimate"] / true == pytest.approx(factor)


def test_bound_holds_for_mixed_first_and_second_order_errors():
    """error = a / n^2 + b / n (a, b >= 0): the contraction of the differences looks better than 1/2 on
    the levels run and tends to 1/2; the old Richardson estimate with the order 2 was too small."""
    rng = np.random.default_rng(1)
    worst_ratio = 0.0
    for _ in range(500):
        a, b = rng.uniform(0.0, 1.0), rng.uniform(0.0, 0.05)
        d1, d2, true = _levels(lambda n: a / n ** 2 + b / n)
        est = estimate_from_levels(d1, d2, 2.0, (1.0, 2.0))["estimate"]
        assert est >= true * (1.0 - 1e-12)
        worst_ratio = max(worst_ratio, est / true)
        assert d2 / 3.0 <= true * (1.0 + 1e-12)                   # the estimate this rule replaces: never above
    assert worst_ratio <= 3.0 + 1e-9                              # never more than the pure second-order factor


def test_bound_does_not_trust_a_contraction_beyond_the_highest_order():
    # errors +4, -1, -1.01 (a sign change between the levels): the last difference is an accident
    d1, d2 = 5.0, 0.01
    est = estimate_from_levels(d1, d2, 2.0, (1.0, 2.0))
    assert est["estimate"] == pytest.approx(d1 / 4.0) and est["estimate"] >= 1.01
    assert estimate_from_levels(d1, d2, 2.0, (1.0, 4.0))["estimate"] == pytest.approx(d1 / 16.0)


def test_bound_with_an_order_below_one_uses_the_observed_contraction():
    d1, d2, true = _levels(lambda n: n ** -0.5)
    est = estimate_from_levels(d1, d2, 2.0, (1.0, 2.0))
    rho = d2 / d1
    assert est["estimate"] == pytest.approx(d2 * rho / (1.0 - rho)) and est["estimate"] >= true * (1 - 1e-12)


def test_no_estimate_when_the_levels_do_not_contract_and_round_off():
    assert estimate_from_levels(1e-4, 1e-4)["estimate"] is None
    out = estimate_from_levels(1e-4, 3e-4)
    assert out["estimate"] is None and out["status"] == "not-asymptotic"
    assert estimate_from_levels(0.0, 1e-5)["estimate"] is None
    r = estimate_from_levels(3e-13, 7e-13)
    assert r["status"] == "round-off" and r["estimate"] == 7e-13 and 7e-13 < ROUNDOFF_FLOOR
    assert estimate_from_levels(0.0, 0.0)["estimate"] == 0.0
    two = estimate_from_levels(2e-4, None, 2.0, (1.0, 2.0))       # two levels only: lowest order assumed
    assert two["status"] == "two-levels" and two["estimate"] == pytest.approx(2e-4)


def test_refinement_levels_of_the_translators():
    from benchfea.crosscheck.foreign import refinement_levels
    assert refinement_levels(32, (1, 2), (8, 16, 32))["levels"] == [8, 16, 32]
    assert refinement_levels(64, (1, 2), (8, 16, 32))["levels"] == [16, 32, 64]
    assert refinement_levels(6, (1, 2))["levels"] == [3, 6]
    assert refinement_levels(5, (1, 2))["levels"] == [5]


def test_one_refinement_parameter_is_counted_once_in_the_estimated_error():
    from benchfea.crosscheck.fidelity import APPROXIMATE, Fidelity
    f = Fidelity("x")
    f.add("secondOrder", APPROXIMATE, control="subdivisions")
    f.add("releases.full", APPROXIMATE, control="subdivisions")
    f.add("other", APPROXIMATE, error=1e-5)
    assert f.estimated_error is None
    f.set_error("subdivisions", 2e-4)
    assert f.estimated_error == pytest.approx(2.1e-4)


# ====================================================================================== acceptance
@pytest.fixture(scope="module")
def table():
    if not (HAVE["opensees"] and HAVE["pynite"]):
        pytest.skip("needs openseespy and PyNiteFEA")
    return acceptance.table(quick=True)


def test_estimated_error_is_not_below_the_true_error(table):
    """The criterion of the estimator: estimate >= true error (foreign solver against the ODE oracle /
    the tableau solver), on every row where the solver represents the model and has an estimate."""
    assert len(table) >= 10 and not [r["model"] for r in table if "error" in r]
    checked = 0
    for r in table:
        for s, d in r["solvers"].items():
            assert d["ok"], (r["model"], s, d)
            if d["unsupported"] or d["estimate"] is None:
                continue
            checked += 1
            assert d["estimate"] + r["referenceUncertainty"] >= d["true"], (r["model"], s, d)
    assert checked >= 18


def test_estimate_is_not_vacuous(table):
    """est / true: 3 for OpenSees (second-order convergence against a bound that assumes first order
    from the last level on); at most 2^4 for PyNite on these models (fourth-order convergence of the
    cubic element); and every estimate stays far below the eligibility limit on these models."""
    for r in table:
        o, p = r["solvers"]["opensees"], r["solvers"]["pynite"]
        assert 2.5 <= o["ratio"] <= 3.5 and 1.8 <= o["observedOrder"] <= 2.1, (r["model"], o)
        assert o["estimate"] < 2.5e-3 and o["eligible"]
        if p["ratio"] is not None:
            assert 1.0 <= p["ratio"] <= 17.0, (r["model"], p)
        assert p["estimate"] < 2.5e-3 and p["eligible"]
        assert p["true"] < 2e-3 and o["true"] < 1e-3


@needs_both
@pytest.mark.parametrize("label", ["derived/beam-column-02", "bank/S03", "bank/S04"])
def test_corrupted_engine_result_does_not_pass_tier_c(label):
    """1e-3 on the largest displacement of a second-order model: caught, because the smaller of the two
    foreign estimates is far below 1e-3; the unchanged result passes."""
    _, raw, cid, _, _ = next(c for c in acceptance.cases(quick=True) if c[0] == label)
    model = build_model(copy.deepcopy(raw))
    results = analyze(model)
    cache = ForeignCache()
    good = tier_c(raw, model, results, cache=cache)
    info = good.details["combinations"][cid]
    assert good.passed, good.summary
    est = [v["estimatedError"] for v in info["solvers"].values()]
    assert min(est) < 2e-4 and max(est) < 1e-3
    assert all(c["value"] <= c["tolerance"] for c in info["comparisons"])
    bad = copy.deepcopy(results)
    u = bad.combinations[cid].displacements
    k = np.unravel_index(np.argmax(np.abs(u[:, :3])), u[:, :3].shape)
    u[k] *= 1.0 + 1.0e-3
    tr = tier_c(raw, model, bad, cache=cache)
    assert tr.passed is False and tr.category == "disagreement"
    assert any("benchfea" in d for d in tr.details["disagreements"])


@needs_both
def test_solver_without_contraction_is_not_eligible(monkeypatch):
    """Three levels that do not contract give no estimate: the solver is reported as not eligible with
    the reason, never as passing."""
    from benchfea.crosscheck import compare as cmp_mod, cross_check, runner
    _, raw, cid, _, _ = next(c for c in acceptance.cases(quick=True) if c[0] == "derived/beam-column-02")
    model = build_model(copy.deepcopy(raw))
    monkeypatch.setattr(runner, "estimate_from_levels",
                        lambda *a, **k: {"estimate": None, "observedOrder": None, "contraction": 1.2,
                                         "status": "not-asymptotic"})
    rep = cross_check(model, cid, benchfea_result=analyze(model))
    assert rep.eligible == [] and not rep.cross_checked
    assert sum("approximation without an error estimate" in r and "not-asymptotic" in r for r in rep.reasons) == 2


def test_policy_numbers_of_tier_c_are_unchanged():
    assert POLICY["C"] == {"relative": 1.0e-6, "maxApproximation": 5.0e-3, "requiredSolvers": 2,
                           "estimateSafety": 1.0, "scale": "model-wide scales of benchfea.crosscheck.compare"}


# ====================================================================================== measured facts
def _column_document(loads, analysis=None, member=None, nodes=None, supports=None):
    L = 4000.0
    return document(nodes or {"a": (0, 0, 0), "b": (0, 0, L)}, {"c": member or ("a", "b")},
                    supports or {"a": "fixed"}, {"L": loads},
                    combos={"U": {"factors": {"L": 1.0}, "analysisType": "secondOrder"}},
                    analysis={"alpha": 1.0, **(analysis or {})})


@needs_opensees
def test_opensees_pdelta_takes_the_axial_force_of_end_j_with_an_axial_element_load():
    """MEASURED (OpenSees 3.8): one cantilever element with the PDelta transformation.  The lateral tip
    stiffness is 3EI/L^3 - N/L with N the axial force the transformation uses.  With a tip load P it is
    P; with the same total as an axial ``eleLoad`` -beamUniform it is 0 -- the force of end j, not the
    mean w L / 2.  This is why the translator applies axial member loads at the nodes in second order."""
    import subprocess
    code = r'''
import json, sys
import openseespy.opensees as ops
E, A, I, L, H = 2.0e5, 5000.0, 4.0e7, 4000.0, 1000.0
def used(P, w):
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0, 0, 0); ops.node(2, 0, 0, L); ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("PDelta", 1, 0.0, 1.0, 0.0)
    ops.element("elasticBeamColumn", 1, 1, 2, A, E, 7.7e4, 3.0e5, I, I, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(2, H, 0, -P, 0, 0, 0)
    if w:
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", 0.0, 0.0, -w)
    ops.system("UmfPack"); ops.numberer("RCM"); ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-12, 50); ops.algorithm("Newton"); ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return (3 * E * I / L ** 3 - H / ops.nodeDisp(2, 1)) * L
print(json.dumps([used(2.0e5, 0.0), used(0.0, 50.0)]))
'''
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    with_tip_load, with_element_load = json.loads(out.stdout.strip().splitlines()[-1])
    assert with_tip_load == pytest.approx(2.0e5, rel=1e-9)
    assert abs(with_element_load) < 1e-6 * 2.0e5                  # w L = 2e5 N on the element, none of it used


@needs_opensees
def test_opensees_second_order_translation_has_no_axial_element_load_and_converges_at_second_order():
    """A column under its own weight (axial line load) with a lateral tip load: the job carries no
    ``eleLoad`` with an axial component, and the observed order of the three levels is 2 (it was 1 with
    the axial component as an element load)."""
    from benchfea.crosscheck import opensees as ops_tr
    model = build(_column_document([{"type": "nodal", "node": "b", "Fx": 2000.0},
                                    {"type": "memberDistributed", "member": "c", "direction": "Z", "w1": -120.0},
                                    {"type": "memberPoint", "member": "c", "direction": "Z", "value": -5.0e4,
                                     "position": 0.45},
                                    {"type": "memberDistributed", "member": "c", "direction": "X", "w1": 0.4}]))
    t = ops_tr.translate(model, "U")
    native = [c for c in t.job["pattern"] if c[0] == "eleLoad"]
    assert {c[4] for c in native} == {"-beamUniform"}             # the transverse line load stays native
    for c in native:                                              # ..., kind, wy, wz, wx[, a, b] / Py, Pz, x, Px
        assert (c[7] if c[4] == "-beamUniform" else c[8]) == 0.0, c
    r = run("opensees", model, "U")
    conv = r.info["convergence"]["subdivisions"]
    assert r.ok and conv["status"] == "asymptotic" and 1.8 < conv["observedOrder"] < 2.3
    d1, d2 = conv["differences"]
    assert conv["estimatedError"] == pytest.approx(max(d2, d1 / 4.0))


@needs_pynite
def test_pynite_geometric_stiffness_has_an_axial_and_a_torsional_term():
    """MEASURED (PyNiteFEA 3.2.1): ``Member3D.kg(P)`` has P/L on the axial DOFs and P (Iy + Iz) / (A L)
    on the twist -- two terms that benchfea's theory does not have and that the worker removes through
    the section properties."""
    from Pynite import FEModel3D
    m = FEModel3D()
    m.add_node("a", 0, 0, 0)
    m.add_node("b", 1000.0, 0, 0)
    m.add_material("M", E, G, 0.3, 0.0)
    m.add_section("S", SEC["A"], SEC["Iy"], SEC["Iz"], SEC["J"])
    m.add_member("c", "a", "b", "M", "S")
    kg = np.asarray(m.members["c"].kg(7.0e4), float)
    assert kg[0, 0] == pytest.approx(7.0e4 / 1000.0) and kg[0, 6] == pytest.approx(-7.0e4 / 1000.0)
    assert kg[3, 3] == pytest.approx(7.0e4 * (SEC["Iy"] + SEC["Iz"]) / SEC["A"] / 1000.0)
    assert kg[1, 1] == pytest.approx(1.2 * 7.0e4 / 1000.0)


@needs_pynite
def test_pynite_route_iterates_the_axial_forces_and_removes_its_extra_terms():
    """Beam-column with a strong second-order effect: the PyNite route (axial forces iterated with
    PyNite's own matrices, axial and torsional term of its geometric stiffness removed) reproduces the
    ODE oracle within its estimate, which is now the discretisation bound: no N/EA term of 1e-4."""
    from benchfea.oracle import Beam
    Lm, P, H = 4000.0, 0.5 * math.pi ** 2 * EIz / (2 * 4000.0) ** 2, 1000.0
    model = build(_column_document([{"type": "nodal", "node": "b", "Fz": -P, "Fx": H}]))
    ref = Beam(Lm, EIz, N=-P, second_order=True).support(0.0, v="fixed", r="fixed").point(Lm, P=H).solve()
    r = run("pynite", model, "U", in_process=True)
    so = r.info["secondOrder"]
    assert r.ok and so["iterations"] >= 1 and so["iterationChange"] <= 1e-10
    assert so["formulationResidual"] < 1e-9 and so["postProcessing"] < 1e-4
    est = r.fidelity.estimated_error
    assert P / EA > 5e-5 and est < 5e-5                           # the removed axial term alone was N / EA
    assert abs(r.displacements[1, 0] / ref.v(Lm) - 1.0) <= est
    assert abs(r.displacements[1, 2] + P * Lm / EA) <= 1e-6 * P * Lm / EA       # axial stiffness is EA / L
    items = {i.feature for i in r.fidelity.approximate}
    assert {"secondOrder", "secondOrder.axialForceIteration", "secondOrder.geometricStiffnessTerms",
            "secondOrder.postProcessing"} <= items


# ====================================================================================== constraints
def _diaphragm_storey():
    sys.path.insert(0, str(ROOT / "tests" / "independent"))
    try:
        import models as bank
    finally:
        sys.path.remove(str(ROOT / "tests" / "independent"))
    return bank.SECOND["S08_diaphragm_storey_with_leaning_columns"]()


@needs_both
def test_rigid_diaphragm_in_second_order_pynite_is_unsupported_and_the_omitted_term_is_real():
    """Bank model S08.  PyNite's stiff members carry geometric stiffness from the forces the diaphragm
    transmits; benchfea and OpenSees use a linear kinematic constraint without it.  The record says
    "unsupported" for PyNite (not "approximate"), so the model is not cross-checked.

    The term itself, from the mechanics (docs/CROSSCHECK.md 3.5): a rigid body that turns by theta about
    the normal moves a slave at r by theta x r - theta^2 r / 2; the in-plane force f_s it applies to the
    slave does the second-order work -theta^2 (f_s . r_s) / 2, i.e. the rotational stiffness
    -sum f_s . r_s at the master.  Added to the engine's model as a rotational spring (iterated to its
    fixed point), it changes the largest displacement of this model by several per cent: the linear
    constraint and the rigid body are different problems in second order."""
    raw, cid = _diaphragm_storey()
    model = build_model(copy.deepcopy(raw))
    res = analyze(model)
    c = res.combinations[cid]
    assert any("linear kinematic constraints" in w for w in c.warnings)
    ops, pyn = run("opensees", model, cid), run("pynite", model, cid, in_process=True)
    assert ops.ok and ops.fidelity.level == "approximate" and not ops.fidelity.unsupported
    assert [i.feature for i in pyn.fidelity.unsupported] == ["constraints.rigidDiaphragm.secondOrder"]
    rec = tier_c(raw, model, res, cache=ForeignCache())
    assert rec.passed is False and rec.category == "not-representable"
    assert rec.details["combinations"][cid]["eligible"] == ["opensees"]
    # the engine agrees with OpenSees (same statement of the problem) within the OpenSees bound
    from benchfea.crosscheck import from_benchfea
    assert relative_difference(from_benchfea(model, res, cid), ops, model) <= ops.fidelity.estimated_error

    ids = list(model.nodes)
    xm = np.array(model.nodes["u4"].xyz)

    def solve(doc):
        m = build_model(copy.deepcopy(doc))
        cc = analyze(m).combinations[cid]
        cf, u = np.asarray(cc.constraint_forces), np.asarray(cc.displacements)
        k = -sum(cf[ids.index(s), :2] @ (np.array(m.nodes[s].xyz) - xm)[:2] for s in ("u0", "u1", "u2", "u3"))
        return u[ids.index("u4")], k                              # N mm / rad

    u_lin, k = solve(raw)
    for _ in range(8):
        doc = copy.deepcopy(raw)
        doc["supports"]["u4"] = {"dofs": {"rz": {"spring": k / 1.0e6}}}       # kN m / rad
        u_body, k = solve(doc)
    assert k > 0.0                                                # the links are in tension: stabilising
    change = abs(u_body[0] / u_lin[0] - 1.0)
    assert 0.05 < change < 0.10                                   # 7.6 % measured
    # PyNite's stiff members contain the term (their axial forces are iterated like those of the members):
    # its result is the rigid-body statement to within 1 %, and 8 % away from the linear constraint
    u_p = pyn.displacements[ids.index("u4")]
    assert abs(u_p[0] / u_body[0] - 1.0) < 0.02 and abs(u_p[0] / u_lin[0] - 1.0) > 0.05


# ====================================================================================== truss with axial load
@needs_both
def test_truss_member_with_an_axial_member_load_is_a_beam_column_in_the_foreign_models():
    """Bank model S06 (cause of its former disagreement).  A truss-type member that carries a member load
    with a component along its axis is a pin-ended beam-column in benchfea: it bows, and the axial load
    acts on the bowed axis.  A single foreign element without bending cannot show this (difference
    2.5e-4, with the tableau solver converging to the ENGINE as its pieces are refined).  The
    translators now subdivide such a member with its own bending stiffness and hinges at the ends."""
    sys.path.insert(0, str(ROOT / "tests" / "independent"))
    try:
        import models as bank
    finally:
        sys.path.remove(str(ROOT / "tests" / "independent"))
    raw, cid = bank.SECOND["S06_truss_members_with_transverse_load_and_partial_releases"]()
    model = build_model(copy.deepcopy(raw))
    from benchfea.crosscheck import from_benchfea, opensees as ops_tr, pynite as pyn_tr
    from benchfea.crosscheck.prepare import prepare
    p = prepare(model, cid)
    assert p.members["strut"].bows(True) and p.members["strut2"].bows(True) and not p.members["c1"].bows(True)
    assert not p.members["strut"].bows(False)
    t = ops_tr.translate(model, cid)
    n_nodes = sum(1 for c in t.job["model"] if c[0] == "node")
    assert n_nodes >= len(model.nodes) + 5 * 31                   # five members of 32 elements each
    assert pyn_tr.translate(model, cid).controls["subdivisions"]["levels"] == [2, 4, 8]
    eng = from_benchfea(model, analyze(model), cid)
    ref = acceptance.tableau_reference(raw, model, cid, 64)
    assert relative_difference(eng, ref, model) < 1e-6            # engine = tableau (64 pieces)
    for s in ("opensees", "pynite"):
        r = run(s, model, cid, in_process=(s == "pynite"))
        assert r.ok and r.fidelity.estimated_error >= relative_difference(r, ref, model)
        assert relative_difference(r, eng, model) <= r.fidelity.estimated_error + 1e-6
    # a truss member whose section has no bending stiffness cannot be modelled as a beam-column
    raw2 = copy.deepcopy(raw)
    raw2["sections"]["T"] = {"type": "properties", "A": raw["sections"]["S"]["A"]}
    raw2["members"]["strut"]["section"] = "T"
    try:
        m2 = build_model(raw2)
    except Exception:
        return                                                    # the model layer refuses such a section
    assert [i.feature for i in ops_tr.translate(m2, cid).fidelity.unsupported] == ["members.truss"]


# ====================================================================================== pipeline gaps
def _braced_bay(self_weight=True, second_tie_load=False):
    d = document({"a": (0, 0, 0), "b": (4000, 0, 0), "A": (0, 0, 3000), "B": (4000, 0, 3000)},
                 {"c1": ("a", "A"), "c2": ("b", "B"), "g": ("A", "B"),
                  "x1": {"i": "a", "j": "B", "type": "tensionOnly"}, "x2": {"i": "b", "j": "A", "type": "tensionOnly"}},
                 {"a": "fixed", "b": "fixed"},
                 {"W": [{"type": "nodal", "node": "A", "Fx": 30000.0}] +
                       ([{"type": "selfWeight"}] if self_weight else []) +
                       ([{"type": "memberPoint", "member": "x2", "direction": "Z", "value": -900.0, "position": 0.3},
                         {"type": "memberDistributed", "member": "x2", "direction": "x", "w1": 0.2, "w2": 0.5}]
                        if second_tie_load else [])})
    return d


@needs_both
def test_slack_one_sided_member_with_member_loads_is_cross_checked():
    """examples/05 in small: the slack diagonal carries its self weight, a point load and an axial line
    load.  Tier C gives its loads to the foreign solvers as the end reactions of a simply supported beam
    (docs/THEORY.md 13), computed from the load definitions; the engine agrees to round-off."""
    raw = _braced_bay(second_tie_load=True)
    model = build_model(copy.deepcopy(raw))
    res = analyze(model)
    c = res.combinations["W"]
    assert c.members["x1"].active and not c.members["x2"].active
    loads = slack_member_loads(model, "W", ["x2"])
    total = np.sum([v[:3] for v in loads.values()], axis=0)
    m = model.members["x2"]
    L = float(m.length)
    w = 7.7e-5 * SEC["A"]
    ex = np.asarray(m.R)[0]
    expect = np.array([0.0, 0.0, -w * L - 900.0]) + ex * 0.35 * L           # self weight + point + axial line load
    assert set(loads) == {"b", "A"} and total == pytest.approx(expect, rel=1e-12, abs=1e-9)
    tr = tier_c(raw, model, res, cache=ForeignCache())
    info = tr.details["combinations"]["W"]
    assert tr.passed, tr.summary
    assert info["substitutions"][0]["slackMembersWithLoads"] == ["x2"]
    assert max(cmp["value"] for cmp in info["comparisons"]) < 1e-9
    assert info["complementarity"]["violations"] == []
    # tier D: the slack member's end forces are the simply supported reactions, its stations follow by statics
    rec = verify_model(raw, tiers="D")
    assert rec.tiers["D"].passed, rec.tiers["D"].summary
    bad = copy.deepcopy(res)
    bad.combinations["W"].members["x2"].end_forces[1] *= 1.01
    from benchfea.verify.stations import station_statics
    out = station_statics(model, bad, POLICY["D"]["stationStatics"], POLICY["D"]["stationStaticsIterated"])
    assert not out["passed"] and any("simply supported" in f for f in out["failed"])


@needs_both
def test_loaded_node_without_members_is_a_verifiable_diagnosis():
    """A load on a node that no member holds (and a moment on a node joined only by truss members): the
    engine refuses with a "no-stiffness" diagnosis.  Tier N keeps such loaded degrees of freedom in the
    independent stiffness, so the diagnosis is confirmed (it could not be before: they were dropped)."""
    d = document({"a": (0, 0, 0), "b": (0, 0, 4000), "x": (2000, 0, 4000), "t": (3000, 0, 0)},
                 {"c": ("a", "b"), "tr": {"i": "t", "j": "b", "type": "truss"}},
                 {"a": "fixed", "t": "pinned"},
                 {"P": [{"type": "nodal", "node": "b", "Fx": 12000.0},
                        {"type": "nodal", "node": "x", "Fx": 1000.0, "Mz": 2.0e6}]})
    rec = verify_model(d)
    assert rec.kind == "diagnosis" and rec.status == "verified" and "".join(rec.tiers_passed) == "NKCE"
    assert rec.reference["diagnosis"]["nMechanisms"] == 2
    assert "2 zero-stiffness motion(s)" in rec.tiers["N"].summary
    assert "2 loaded DOF(s) without stiffness" in rec.tiers["K"].summary
    ev = rec.tiers["C"].details["evidence"]
    assert {v["verdict"] for v in ev.values()} == {"consistent"}
    # an unloaded node without members is still not part of the solution (no mechanism)
    d2 = copy.deepcopy(d)
    d2["loadCases"]["P"]["loads"].pop()
    assert verify_model(d2, tiers="D").kind == "results"
    # the independent rule must find the load itself: a diagnosis that names an unloaded DOF is not confirmed
    from benchfea.verify.diagnosis import tier_null_space
    diag = copy.deepcopy(rec.reference["diagnosis"])
    diag["modes"][0]["nodes"][0]["dofs"] = {"ux": 0.5, "uy": 0.5}
    out = tier_null_space(build_model(copy.deepcopy(d)), diag)
    assert out.passed is False and "not a loaded degree of freedom" in out.summary + " ".join(out.details["failed"])


def test_member_split_keeps_clear_of_load_positions():
    """A point load 1e-5 L beside a result station (Ziemian frame 2: loads at 0.25 +- 1e-5 of inclined
    rafters, five stations): the cut is not made at that station, where it would leave a sliver between
    the new joint and the load.  The split tolerance is unchanged."""
    from benchfea.verify.tier_e import MIN_PIECE, split_document
    loads = [{"type": "memberPoint", "member": "g", "direction": "Z", "value": -20000.0, "position": 0.25 + 1e-5},
             {"type": "memberPoint", "member": "g", "direction": "Z", "value": -20000.0, "position": 0.75 - 1e-5},
             {"type": "nodal", "node": "A", "Fx": 5000.0}]
    d = document({"a": (0, 0, 0), "A": (0, 0, 3000), "B": (6000, 0, 4000), "b": (6000, 0, 0)},
                 {"c1": ("a", "A"), "g": ("A", "B"), "c2": ("b", "B")}, {"a": "fixed", "b": "fixed"}, {"L": loads},
                 combos={"U": {"factors": {"L": 1.0}, "analysisType": "secondOrder"}},
                 analysis={"alpha": 1.0, "resultStations": 5})
    for n in ("A", "B"):
        d["supports"][n] = {"dofs": {"uy": "fixed"}}
    model = build_model(copy.deepcopy(d))
    doc, cuts, skipped = split_document(d, model)
    assert cuts["g"]["station"] == 2 and "g" not in skipped         # 0.5 L, not 0.25 L
    assert MIN_PIECE == 1e-3 and POLICY["E"]["splitNonlinear"] == 1.0e-6
    rec = verify_model(d, tiers="E")
    split = rec.tiers["E"].details["checks"]["memberSplit"]
    assert split["passed"] is True, split
    # every interior station next to a load: the member is skipped with the reason, not cut
    d3 = copy.deepcopy(d)
    d3["analysis"]["resultStations"] = 3
    d3["loadCases"]["L"]["loads"][0]["position"] = 0.5 + 1e-5
    _, cuts3, skipped3 = split_document(d3, build_model(copy.deepcopy(d3)))
    assert "g" not in cuts3 and "closer than" in skipped3["g"]
    # the end of a line load exactly at the station is no obstacle
    d4 = copy.deepcopy(d)
    d4["loadCases"]["L"]["loads"] = [{"type": "memberDistributed", "member": "g", "direction": "Z", "w1": -5.0,
                                      "start": 0.25, "end": 0.75}, loads[2]]
    _, cuts4, _ = split_document(d4, build_model(copy.deepcopy(d4)))
    assert cuts4["g"]["station"] == 1
