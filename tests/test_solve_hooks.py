"""The extension points a nonlinear driver builds on: analysis-type dispatch, per-member stiffness
factors, inactive members, assembly of arbitrary per-member matrices, per-load-set solve."""
import numpy as np
import pytest

from benchfea.model import build_model
from benchfea.solve import (ANALYSIS_TYPES, LinearAnalysis, ModelRefused, StiffnessFactors, analyze, build_element,
                            combination_loads, station_layout)

from solve_reference import EA, EIz, GAy, beam_bvp, clamped, make

L = 3000.0


def doc(**kw):
    return make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed"},
                {"c": [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -4.0},
                       {"type": "nodal", "node": "b", "Fx": 7000.0}]}, **kw)


def test_other_analysis_types_are_refused_until_a_solver_registers_for_them(monkeypatch):
    d = doc()
    d["combinations"]["c"]["analysisType"] = "secondOrder"
    model = build_model(d)
    with pytest.raises(ModelRefused) as e:
        analyze(model)
    assert e.value.codes == ["analysis.not-implemented"] and e.value.problems[0].path == "/combinations/c/analysisType"
    calls = []

    def fake(la, combination):
        calls.append(combination.id)
        loads = la.combination_loads(combination)
        sol, rec = la.solve_loads(loads)
        return la.result(combination, loads, sol, rec, iterations=7)
    monkeypatch.setitem(ANALYSIS_TYPES, "secondOrder", fake)
    res = analyze(model)
    assert calls == ["c"] and res.combinations["c"].iterations == 7 and res.combinations["c"].analysis_type == "secondOrder"


def test_stiffness_factors_scale_each_section_stiffness():
    model = build_model(doc())
    m = model.members["m"]
    el = build_element(model, m, StiffnessFactors(axial=0.5, flexural=0.8, torsional=0.25, shear=0.9))
    la = LinearAnalysis(model, {"m": el})
    sol, rec = la.solve_loads(la.combination_loads(model.combinations["c"]))
    ref = beam_bvp(L, 0.8 * EIz, 0.9 * GAy, clamped, lambda s: [s[2], s[3]], [(0, L, -4.0, -4.0)])
    assert sol.u[6 + 2] == pytest.approx(ref(L)[0], rel=1e-9)
    assert sol.u[6 + 0] == pytest.approx(7000.0 * L / (0.5 * EA), rel=1e-12)
    base = build_element(model, m)
    assert el.k_local[3, 3] == pytest.approx(0.25 * base.k_local[3, 3], rel=1e-12)


def test_inactive_member_has_no_stiffness_but_still_hands_its_loads_to_the_nodes():
    d = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b"), "tie": {"i": "a", "j": "b", "type": "truss"}}, {"a": "fixed"},
             {"c": [{"type": "nodal", "node": "b", "Fx": 7000.0},
                    {"type": "memberPoint", "member": "tie", "direction": "Z", "value": -900.0, "position": 0.5}]})
    model = build_model(d)
    both = LinearAnalysis(model)
    off = LinearAnalysis(model, {"m": build_element(model, model.members["m"]),
                                 "tie": build_element(model, model.members["tie"], active=False)})
    assert not off.elements["tie"].k_global.any() and not off.elements["tie"].active
    loads = combination_loads(model, model.combinations["c"])
    s_both, r_both = both.solve_loads(loads)
    s_off, r_off = off.solve_loads(loads)
    assert s_both.u[6] == pytest.approx(7000.0 * L / (2 * EA), rel=1e-12)       # two bars in parallel
    assert s_off.u[6] == pytest.approx(7000.0 * L / EA, rel=1e-12)              # the tie carries nothing
    assert abs(r_off["tie"].end_forces[0]) <= 1e-9 and r_both["tie"].end_forces[6] == pytest.approx(3500.0, rel=1e-12)
    assert s_off.reactions[2] == pytest.approx(900.0, rel=1e-12)                # its span load still reaches the nodes
    assert off.checks(loads, s_off, r_off).passed
    assert off.result(model.combinations["c"], loads, s_off, r_off).members["tie"].active is False


def test_assemble_accepts_any_set_of_per_member_matrices():
    model = build_model(doc())
    la = LinearAnalysis(model)
    st = la.structure
    np.testing.assert_array_equal(st.assemble(), st.assemble({mid: e.k_global for mid, e in la.elements.items()}))
    extra = {"m": 0.1 * la.elements["m"].k_global}                  # e.g. a geometric stiffness contribution
    K = st.assemble() + st.assemble(extra)
    K[np.diag_indices(st.n)] += st.spring
    f = st.factorize(K)
    F = np.zeros(st.n)
    F[6 + 2] = -1000.0
    assert st.solve(f, F).u[8] == pytest.approx(st.solve(la.factor, F).u[8] / 1.1, rel=1e-12)
    assert not st.assemble({}).any()


def test_station_layout():
    assert station_layout(10.0, 3) == [(0.0, -1), (5.0, -1), (10.0, 1)]
    lay = station_layout(10.0, 3, point_positions=[5.0, 2.0, 0.0, 10.0], line_ends=[2.0, 7.5])
    assert lay == [(0.0, -1), (2.0, -1), (2.0, 1), (5.0, -1), (5.0, 1), (7.5, -1), (10.0, 1)]
