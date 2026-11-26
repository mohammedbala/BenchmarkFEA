"""Regression tests for defects of the engine's self-checks and of the active-set iteration.

* invariants of a solution whose forces vanish (rigid-body motion under a support settlement): the
  residuals are measured against the round-off floor ROUND_OFF x max(|K| |u|), not against forces that are
  rounding errors of zero -- and a loaded case is measured exactly as before (the floor is far below);
* superposition check: a member without force in a combination no longer gives a ratio of two rounding
  errors, and a corrupted result is still caught at the same size;
* second order with tension-only bracing: a trial active set in which tension-only members would carry
  compression is not final, so its instability is not the structure's;
* truss-type members are pin-ended struts with the I of their section: the buckling result says so when
  that governs.
"""
import copy
import json
import pathlib
import sys

import numpy as np
import pytest

from benchfea.model import build_model
from benchfea.solve import LinearAnalysis, ModelRefused, analyze, superposition_check
from benchfea.solve import invariants as inv

ROOT = pathlib.Path(__file__).resolve().parents[1]
SETTLEMENT = ROOT / "verification" / "problems" / "staad" / "frame-support-settlement" / "model.json"
EQUILIBRIUM = ("global force equilibrium", "global moment equilibrium", "joint equilibrium",
               "member equilibrium", "work = strain energy")


def _doc(nodes, members, supports, cases, combos, **kw):
    d = {"formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"}, "conventions": {"verticalAxis": "Z"},
         "nodes": {k: dict(zip("xyz", map(float, v))) for k, v in nodes.items()},
         "supports": supports, "materials": {"S": {"E": 200000.0, "nu": 0.3}},
         "sections": {"sec": {"type": "properties", "A": 8000.0, "Iy": 2.0e7, "Iz": 1.5e8, "J": 5.0e5}},
         "members": {k: {"section": "sec", "material": "S", **v} for k, v in members.items()},
         "loadCases": cases, "combinations": combos}
    d.update(kw)
    return d


# ====================================================================================== invariants
def test_rigid_body_settlement_of_a_determinate_frame_is_accepted():
    """Statically determinate frame, support settlement only: every force is zero in exact arithmetic.
    Before the fix the equilibrium checks returned 0.3 ... 1.0 (round-off / round-off) and the engine
    refused the solution."""
    model = build_model(json.loads(SETTLEMENT.read_text()))
    res = analyze(model)                                       # used to raise ModelRefused (solve.invariant-failed)
    for c in res.combinations.values():
        assert float(np.abs(c.displacements).max()) > 1.0      # the frame moves (mm) ...
        assert float(np.abs(c.reactions).max()) < 1.0e-6       # ... without forces (N)
        checks = {ch.name: ch for ch in c.invariants.checks}
        for name in EQUILIBRIUM:
            assert checks[name].passed and checks[name].value <= 1.0e-7, (name, checks[name].value)


def _portal(settle=False):
    cases = {"G": {"loads": [{"type": "memberDistributed", "member": "beam", "direction": "Z", "w1": -20.0},
                             {"type": "nodal", "node": "b", "Fx": 30.0e3}]}}
    if settle:
        cases = {"S": {"loads": [{"type": "supportDisplacement", "node": "d", "uz": -12.5}]}}
    sup = {"a": {"type": "pinned"}, "d": {"dofs": {"uz": "fixed", "uy": "fixed", "rx": "fixed"}}}
    return _doc({"a": (0, 0, 0), "b": (0, 0, 4000), "c": (6000, 0, 4000), "d": (6000, 0, 0)},
                {"col1": {"i": "a", "j": "b"}, "beam": {"i": "b", "j": "c"}, "col2": {"i": "c", "j": "d"}},
                sup, cases, {"C": {"factors": {next(iter(cases)): 1.0}}})


def _checks(doc, corrupt=None):
    model = build_model(doc)
    la = LinearAnalysis(model)
    loads = la.combination_loads(model.combinations["C"])
    sol, rec = la.solve_loads(loads)
    if corrupt is not None:
        corrupt(la, sol)
    return la, sol, {c.name: c for c in la.checks(loads, sol, rec).checks}


def test_the_floor_does_not_touch_a_loaded_case_and_a_wrong_reaction_is_caught_at_the_same_size():
    la, sol, clean = _checks(_portal())
    Sf, Sm, Se = inv.stiffness_displacement_scale(la.K_members, sol.u, la.structure.size)
    tol = float(la.model.analysis.equilibrium_tolerance)
    fscale = float(np.abs(sol.reactions.reshape(-1, 6)[:, :3]).sum(axis=0).max())
    assert inv.ROUND_OFF * Sf / tol < 1.0e-3 * fscale          # the floor is orders below the forces of the check
    for name in EQUILIBRIUM:
        assert clean[name].passed and clean[name].value < 1.0e-12
    k = int(np.argmax(np.abs(sol.reactions)))                  # a vertical reaction
    assert k % 6 == 2

    def wrong(rel):
        def f(la_, s):
            s.reactions[k] *= 1.0 + rel
        return f
    # a relative error of 3e-6 in one reaction is caught, 3e-7 is not: tolerance 1e-6 of the force scale
    # (sum of |reactions| and |loads| in that direction), exactly the old measure err / fscale
    r = float(abs(sol.reactions[k]))
    _, _, bad = _checks(_portal(), wrong(3.0e-6))
    assert not bad["global force equilibrium"].passed
    assert bad["global force equilibrium"].value == pytest.approx(3.0e-6 * r / fscale, rel=1e-5)
    _, _, ok = _checks(_portal(), wrong(3.0e-7))
    assert ok["global force equilibrium"].passed


def test_a_spurious_force_in_a_force_free_solution_is_still_caught():
    """Rigid-body settlement: forces are zero.  The floor accepts ROUND_OFF x S (S = stiffness x
    displacement scale); a spurious reaction of 1e-9 S -- far below anything a force-based scale could see
    before the fix, when the scale itself was round-off -- fails the check."""
    la, sol, clean = _checks(_portal(settle=True))
    Sf, _, _ = inv.stiffness_displacement_scale(la.K_members, sol.u, la.structure.size)
    assert float(np.abs(sol.reactions).max()) < 1.0e-10 * Sf and all(clean[n].passed for n in EQUILIBRIUM)

    def spurious(la_, s):
        s.reactions[2] += 1.0e-9 * Sf
    _, _, bad = _checks(_portal(settle=True), spurious)
    assert not bad["global force equilibrium"].passed
    tol = float(la.model.analysis.equilibrium_tolerance)
    assert bad["global force equilibrium"].value == pytest.approx(1.0e-9 * tol / inv.ROUND_OFF, rel=1e-3)


def test_ratio_floor_never_reduces_a_scale():
    assert inv._ratio(2.0, 4.0) == 0.5
    assert inv._ratio(2.0, 4.0, floor=1.0, tol=1.0e-6) == 0.5                       # floor 3.6e-9 < 4
    assert inv._ratio(2.0, 0.0, floor=1.0, tol=1.0e-6) == pytest.approx(2.0e-6 / inv.ROUND_OFF)
    assert inv._ratio(0.0, 0.0) == 0.0 and inv._ratio(1.0, 0.0) == float("inf")
    assert inv.ROUND_OFF == 16.0 * np.finfo(float).eps


# ====================================================================================== superposition
def _two_case_frame():
    """Portal with a cantilever stub: the stub carries load in case Q only, so the combination without Q
    leaves it with zero forces while its nodes move."""
    nodes = {"a": (0, 0, 0), "b": (0, 0, 4000), "c": (6000, 0, 4000), "d": (6000, 0, 0), "e": (8000, 0, 4000)}
    members = {"col1": {"i": "a", "j": "b"}, "beam": {"i": "b", "j": "c"}, "col2": {"i": "c", "j": "d"},
               "stub": {"i": "c", "j": "e"}}
    sup = {"a": {"type": "fixed"}, "d": {"type": "fixed"}}
    cases = {"G": {"loads": [{"type": "memberDistributed", "member": "beam", "direction": "Z", "w1": -20.0}]},
             "W": {"loads": [{"type": "nodal", "node": "b", "Fx": 30.0e3}]},
             "Q": {"loads": [{"type": "nodal", "node": "e", "Fz": -5.0e3}]}}
    combos = {"U1": {"factors": {"G": 1.2, "W": 1.6}}, "U2": {"factors": {"G": 1.2, "W": 1.0, "Q": 1.6}}}
    return _doc(nodes, members, sup, cases, combos)


def _old_superposition(model):
    """The check as it was before the fix (array by array, own scale, forces and moments in one array)."""
    la = LinearAnalysis(model)
    out = {}
    for cid, combo in model.combinations.items():
        sol, rec = la.solve_combination(combo)
        ref = [sol.u, sol.reactions] + [x for r in rec.values() for x in (r.forces.ravel(), r.displacements.ravel())]
        acc = None
        for case_id, f in combo.factors.items():
            s, r_ = la.solve_loads(la.case(case_id))
            parts = [s.u, s.reactions] + [x for r in r_.values() for x in (r.forces.ravel(), r.displacements.ravel())]
            acc = [f * p for p in parts] if acc is None else [a + f * p for a, p in zip(acc, parts)]
        out[cid] = max(inv._ratio(float(np.abs(a - b).max(initial=0.0)), float(np.abs(b).max(initial=0.0)))
                       for a, b in zip(acc, ref))
    return out


def test_superposition_with_a_member_that_carries_no_force_in_a_combination():
    model = build_model(_two_case_frame())
    rows = {name: (err, own, floor) for name, err, own, floor in inv.superposition_arrays(model)["U1"]}
    err, own, floor = rows["member stub: station moments"]
    assert floor > 1.0e6 and own < 1.0e-9 * floor              # "forces" of the stub: round-off of k u
    assert err <= inv.ROUND_OFF * floor                        # and the difference is round-off of the same sums
    for ch in superposition_check(model):
        assert ch.passed and ch.value < 1.0e-11, (ch.name, ch.value, ch.detail)
    # the old measure: ratio of two rounding errors, of order one whenever both are not exactly zero
    old = _old_superposition(model)
    assert old["U2"] < 1.0e-11
    assert old["U1"] < 1.0e-11 or old["U1"] > 1.0e-3


@pytest.mark.parametrize("rel, caught", [(3.0e-9, True), (3.0e-10, False)])
def test_superposition_still_catches_a_wrong_force_of_a_loaded_member_at_the_same_size(monkeypatch, rel, caught):
    """Corrupt the moments of the directly solved combination in one loaded member by a relative amount:
    caught above the tolerance 1e-9 of the member's own largest moment, exactly as before the fix."""
    model = build_model(_two_case_frame())
    orig = LinearAnalysis.solve_combination

    def corrupted(self, combination):
        sol, rec = orig(self, combination)
        if combination.id == "U2":
            rec["beam"].forces[:, 5] *= 1.0 + rel
        return sol, rec
    monkeypatch.setattr(LinearAnalysis, "solve_combination", corrupted)
    new = {c.name: c for c in superposition_check(model)}
    old = _old_superposition(model)
    assert new["superposition: U1"].passed
    assert new["superposition: U2"].passed is (not caught)
    assert (old["U2"] > 1.0e-9) is caught
    if caught:
        assert "member beam: station moments" in new["superposition: U2"].detail
        assert new["superposition: U2"].value == pytest.approx(rel, rel=0.05)
        assert old["U2"] == pytest.approx(rel, rel=0.05)       # the same number as the old measure


def test_superposition_catches_a_spurious_force_in_a_force_free_member(monkeypatch):
    """The stub is force free in U1.  A spurious moment of 1e-12 of its stiffness x displacement scale is
    300 times the round-off bound and is caught."""
    model = build_model(_two_case_frame())
    floor = {n: f for n, _, _, f in inv.superposition_arrays(model)["U1"]}["member stub: station moments"]
    orig = LinearAnalysis.solve_combination

    def corrupted(self, combination):
        sol, rec = orig(self, combination)
        if combination.id == "U1":
            rec["stub"].forces[:, 5] += 1.0e-12 * floor
        return sol, rec
    monkeypatch.setattr(LinearAnalysis, "solve_combination", corrupted)
    ch = {c.name: c for c in superposition_check(model)}["superposition: U1"]
    assert not ch.passed and "member stub: station moments" in ch.detail


# ====================================================================================== tension-only, second order
def _braced_example(self_weight):
    raw = json.loads((ROOT / "examples" / "05_braced_frame.json").read_text())
    raw["combinations"] = {"ULS1": raw["combinations"]["ULS1"]}
    raw["combinations"]["ULS1"]["analysisType"] = "secondOrder"
    raw["analysis"]["convergenceTolerance"] = 1.0e-11
    raw["analysis"].pop("buckling", None)
    raw.pop("design", None)
    frame = [m for m, v in raw["members"].items() if v.get("type") != "tensionOnly"]
    loads = raw["loadCases"]["D"]["loads"]
    if self_weight:
        for ld in loads:
            if ld["type"] == "selfWeight":
                ld["members"] = frame                         # (slack members with own loads are outside tier C)
    else:
        raw["loadCases"]["D"]["loads"] = [ld for ld in loads if ld["type"] != "selfWeight"]
    return raw


def test_second_order_with_tension_only_bracing_under_gravity_is_solved_not_refused():
    """examples/05_braced_frame.json, gravity combination ULS1 in second order.  The first trial set has all
    six diagonals active, all in compression; the slender pipes then exceed their own Euler load and the
    tangent stiffness of THAT set is not positive definite.  The engine used to refuse the combination as
    beyond its critical load (factor 0.38) although the structure, with the diagonals that can work, has a
    buckling factor of 7.65."""
    model = build_model(_braced_example(self_weight=True))
    res = analyze(model)
    c = res.combinations["ULS1"]
    assert c.converged and c.invariants.passed
    slack = sorted(m for m, v in c.members.items() if not v.active)
    assert slack == ["X1b", "X2b", "X3b"]
    for mid in ("X1a", "X2a", "X3a"):
        assert c.members[mid].forces[0, 0] > 0.0              # the working diagonals are in tension
    hist = c.solver["activeSet"]["history"]
    assert any("no second-order solution of this trial set" in h.get("note", "") for h in hist)
    # the linear buckling analysis of the same combination agrees that the load is far below critical
    raw = _braced_example(self_weight=True)
    raw["combinations"]["ULS1"]["analysisType"] = "firstOrder"
    raw["analysis"]["buckling"] = {"combinations": ["ULS1"], "modes": 1}
    lam = analyze(build_model(raw)).buckling["ULS1"].factors[0]
    assert lam > 5.0


def test_an_admissible_active_set_that_is_unstable_is_still_refused():
    """The refusal stays for a set without wrong-signed one-sided members: here the gravity load is
    multiplied beyond the buckling factor of the braced frame."""
    raw = _braced_example(self_weight=True)
    raw["combinations"]["ULS1"]["factors"] = {k: 12.0 * v for k, v in raw["combinations"]["ULS1"]["factors"].items()}
    with pytest.raises(ModelRefused) as err:
        analyze(build_model(raw))
    assert any(p.code == "stability.beyond-critical-load" for p in err.value.problems)


def test_second_order_tension_only_result_against_the_tableau_reference_with_the_same_active_set():
    """Independent reference (tests/independent/tableau.py, no stiffness matrix): the model of the engine's
    final active set (working diagonals as truss members, slack ones removed), second order, constant axial
    forces (no self weight) so that the reference is exact.  The sign conditions of all six diagonals are
    evaluated from the REFERENCE's own displacements and forces."""
    sys.path.insert(0, str(ROOT / "tests" / "independent"))
    try:
        from harness import compare, run_benchfea, worst
    finally:
        sys.path.remove(str(ROOT / "tests" / "independent"))
    from benchfea.verify.tier_c import active_set_document
    raw = _braced_example(self_weight=False)
    model, res, doc = run_benchfea(raw)
    c = res.combinations["ULS1"]
    twin, slack, problem = active_set_document(raw, model, c)
    assert problem is None and sorted(slack) == ["X1b", "X2b", "X3b"]
    u = model.units
    for sid, s in model.sections.items():                      # the reference reads explicit properties only
        twin["sections"][sid] = {"type": "properties", "A": float(u.from_internal(s.A, "length^2")),
                                 **{k: float(u.from_internal(getattr(s, k), "length^4")) for k in ("Iy", "Iz", "J")}}
    doc2 = copy.deepcopy(doc)
    for mid in slack:
        del doc2["combinations"]["ULS1"]["members"][mid]
    out = compare(twin, "ULS1", second_order=True, doc=doc2)
    assert worst(out) <= 1.0e-9, {k: v for k, v in out.items() if not k.startswith("_")}
    ref, ids = out["_ref"], list(raw["nodes"])
    for mid, m in raw["members"].items():
        if m.get("type") != "tensionOnly":
            continue
        xi, xj = (np.array([raw["nodes"][n][k] for k in "xyz"], float) for n in (m["i"], m["j"]))
        e = (xj - xi) / np.linalg.norm(xj - xi)
        elongation = float(e @ (ref["displacements"][ids.index(m["j"])][:3] - ref["displacements"][ids.index(m["i"])][:3]))
        if mid in slack:
            assert elongation < 0.0, mid                       # it would be shortened: rightly slack
        else:
            assert ref["members"][mid]["forces"][0, 0] > 0.0 and elongation > 0.0, mid


def test_second_order_tension_only_result_against_opensees_and_pynite_with_the_same_active_set():
    """Tier C of the verification pipeline on the same combination (with self weight on the frame members):
    OpenSeesPy and PyNite get the linear model of the engine's final active set; the sign conditions are
    evaluated from their results; tolerances are the pipeline's (benchfea.verify.policy), unchanged."""
    from benchfea import crosscheck
    if not (crosscheck.available("opensees") and crosscheck.available("pynite")):
        pytest.skip("OpenSeesPy / PyNite are not installed")
    from benchfea.verify.tier_c import tier_c
    raw = _braced_example(self_weight=True)
    model = build_model(copy.deepcopy(raw))
    res = analyze(model)
    r = tier_c(raw, model, res)
    info = r.details["combinations"]["ULS1"]
    assert r.passed and info["crossChecked"] and info["complementarity"]["violations"] == []
    assert info["substitutions"][0]["slack"] == ["X1b", "X2b", "X3b"]
    by = {(c["a"], c["b"]): c for c in info["comparisons"]}
    assert by[("benchfea", "opensees")]["passed"] and by[("benchfea", "opensees")]["value"] < 1.0e-5


# ====================================================================================== truss-type members as struts
def _two_bar(i_bar):
    doc = _doc({"a": (0, 0, 0), "b": (3000, 0, 0), "c": (3000, 0, -4000)},
               {"tie": {"i": "a", "j": "b", "type": "truss"}, "strut": {"i": "c", "j": "b", "type": "truss"}},
               {"a": {"type": "pinned"}, "c": {"type": "pinned"},
                "b": {"dofs": {"uy": "fixed"}}},
               {"P": {"loads": [{"type": "nodal", "node": "b", "Fz": -1000.0}]}}, {"C": {"factors": {"P": 1.0}}},
               analysis={"shearDeformation": False, "buckling": {"combinations": ["C"], "modes": 1}})
    doc["sections"]["sec"].update(Iy=i_bar, Iz=i_bar)
    return doc


def test_buckling_result_says_when_the_lowest_mode_is_the_own_buckling_of_a_truss_member():
    """A truss-type member is a pin-ended strut with the Iy, Iz of its section (THEORY 10.6).  With a
    placeholder I the lowest buckling factor is the strut's own Euler factor pi^2 E I / (L^2 N); the result
    is unchanged (intended behaviour) and now carries a warning that names the member."""
    b = analyze(build_model(_two_bar(1.0))).buckling["C"]
    euler = np.pi ** 2 * 200000.0 * 1.0 / 4000.0 ** 2 / 1000.0          # strut: N = -1000 N, L = 4000
    assert b.factors[0] == pytest.approx(euler, rel=1e-6)
    w = b.info["warnings"]
    assert len(w) == 1 and "truss-type member 'strut'" in w[0] and "pin-ended strut" in w[0]
    # second order: the refusal carries the same explanation
    doc = _two_bar(1.0)
    doc["combinations"]["C"]["analysisType"] = "secondOrder"
    with pytest.raises(ModelRefused) as err:
        analyze(build_model(doc))
    assert any("truss-type member 'strut'" in p.message for p in err.value.problems)
    # the factor follows the section's I (here the strut's own mode is the only buckling mode of the truss)
    b = analyze(build_model(_two_bar(1.0e10))).buckling["C"]
    assert b.factors[0] == pytest.approx(1.0e10 * euler, rel=1e-6)
    # a truss whose lowest mode is a mode of the STRUCTURE (McGuire et al. Example 9.1 as modelled in
    # verification/problems: bar bc with I = 1e10 mm^4) carries no such warning
    p = ROOT / "verification" / "problems" / "mgz" / "ex-9-1-truss-critical-load" / "model.json"
    b = analyze(build_model(json.loads(p.read_text()))).buckling["C"]
    assert not any("truss-type" in x for x in b.info["warnings"])
