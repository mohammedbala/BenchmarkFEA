"""Three-way comparison benchfea.solve vs OpenSees vs PyNite on models built through the public input
format, and the automatic tier-C decision (``cross_check``).

Skips cleanly when a foreign package or ``benchfea.solve.analyze`` is missing; the second-order cases
activate by themselves once ``benchfea.solve`` registers the "secondOrder" analysis type.
"""
import numpy as np
import pytest

from benchfea.crosscheck import compare, cross_check, from_benchfea
from benchfea.crosscheck.fidelity import EXACT
from benchfea.model import build_model

from crosscheck_helpers import SEC, document, needs_both, solve

solve_pkg = pytest.importorskip("benchfea.solve")
pytestmark = needs_both

NODES = {"A": (0, 0, 0), "B": (0, 0, 3000), "C": (4000, 1000, 3500), "D": (4000, 1000, 0)}
MEMBERS = {"c1": ("A", "B"), "b": ("B", "C"), "c2": ("D", "C")}
SUPPORTS = {"A": "fixed", "D": "pinned"}
LOADS = [
    {"type": "nodal", "node": "B", "Fx": 2e4, "My": 3e6},
    {"type": "memberDistributed", "member": "b", "direction": "Z", "w1": -5, "w2": -12, "start": 0.1, "end": 0.8},
    {"type": "memberDistributed", "member": "c1", "direction": "y", "w1": 2.0},
    {"type": "memberDistributed", "member": "c2", "direction": "X", "w1": 1.0, "start": 0.3, "end": 0.6},
    {"type": "memberDistributed", "member": "b", "direction": "Z", "w1": -1.5, "projected": True},
    {"type": "memberPoint", "member": "b", "direction": "y", "value": -8e3, "position": 0.4},
    {"type": "memberPoint", "member": "c1", "direction": "x", "value": -8e3, "position": 0.7},
    {"type": "memberPoint", "member": "c2", "direction": "Y", "kind": "moment", "value": 4e6, "position": 0.5},
    {"type": "memberPoint", "member": "b", "direction": "x", "kind": "moment", "value": 4e6, "position": 0.25},
    {"type": "selfWeight"},
    {"type": "supportDisplacement", "node": "A", "uz": -2.0, "ry": 0.001},
]


def _frame(members=None, supports=None, loads=None, **kw):
    return document(kw.pop("nodes", NODES), {**MEMBERS, **(members or {})}, supports or SUPPORTS,
                    {"L": LOADS if loads is None else loads}, **kw)


def m01_all_load_types():
    return _frame(), "both"


def m02_rolled_member_full_releases():
    return _frame({"b": {"i": "B", "j": "C", "orientation": {"roll": 25},
                         "releases": {"j": {"rz": True}, "i": {"ry": True}}}}), "both"


def m03_partial_rotational_springs():
    return _frame({"b": {"i": "B", "j": "C", "orientation": {"roll": -40},
                         "releases": {"j": {"rz": {"spring": 5e9}, "ry": {"spring": 1e9}},
                                      "i": {"rz": {"spring": 2e10}, "ry": True}}}}), "both"


def m04_offsets_and_releases():
    return _frame({"c2": {"i": "D", "j": "C", "offsets": {"j": [0, 0, -200], "i": [0, 0, 150]},
                          "releases": {"j": {"rz": {"spring": 5e9}, "ry": True}}},
                   "b": {"i": "B", "j": "C", "offsets": {"axes": "local", "i": [180, 0, 0]}}},
                  supports={"A": "fixed", "D": "fixed"}), "both"


def m05_axial_spring_and_torsion_release():
    return _frame({"c1": {"i": "A", "j": "B", "releases": {"j": {"ux": {"spring": 5e4}, "rx": True,
                                                                 "uy": {"spring": 2e4}}}}}), "opensees"


def m06_support_springs():
    return _frame(supports={"A": "fixed", "D": {"type": "pinned", "dofs": {"rz": {"spring": 1e9},
                                                                            "ux": {"spring": 5e3}}}}), "both"


def m07_truss_and_frame():
    loads = [ld for ld in LOADS if ld.get("member") != "c2"] + [
        {"type": "memberDistributed", "member": "t", "direction": "Z", "w1": -1, "w2": -3}]
    return _frame({"t": {"i": "A", "j": "C", "type": "truss"}, "c2": {"i": "D", "j": "C", "type": "truss"}},
                  loads=loads), "both"


def m08_rigid_diaphragm():
    nodes = {**NODES, "B": (0, 0, 3500), "M": (1500, 2500, 3500), "E": (2000, 3000, 3900)}
    return _frame({"e": ("C", "E")}, {**SUPPORTS, "E": "pinned"}, nodes=nodes,
                  constraints={"d": {"type": "rigidDiaphragm", "master": "M", "slaves": ["B", "C"]}}), "both"


def m09_rigid_link():
    nodes = {**NODES, "B": (0, 0, 3500), "M": (1500, 2500, 3500), "E": (2000, 3000, 3900)}
    doc = document(nodes, {"c1": ("A", "B"), "c2": ("D", "C"), "e": ("M", "E")}, {**SUPPORTS, "E": "pinned"},
                   {"L": [ld for ld in LOADS if ld.get("member") != "b"]},
                   constraints={"d": {"type": "rigidLink", "master": "B", "slaves": ["C", "M"]}})
    return doc, "both"


def m10_shear_deformation_with_releases_and_offsets():
    doc, _ = m04_offsets_and_releases()
    doc["analysis"]["shearDeformation"] = True
    return doc, "opensees"


def m11_y_up_kip_inch():
    doc = {
        "formatVersion": "1.0.0", "units": {"force": "kip", "length": "in"}, "conventions": {"verticalAxis": "Y"},
        "nodes": {"N1": {"x": 0, "y": 0, "z": 0}, "N2": {"x": 0, "y": 144, "z": 0}, "N3": {"x": 240, "y": 156, "z": 30},
                  "N4": {"x": 240, "y": 0, "z": 30}},
        "supports": {"N1": {"type": "fixed"}, "N4": {"type": "pinned", "dofs": {"ry": "fixed"}}},
        "materials": {"steel": {"E": 29000, "G": 11200, "unitWeight": 0.000284}},
        "sections": {"col": {"type": "I", "d": 14, "bf": 10, "tf": 0.7, "tw": 0.44},
                     "beam": {"type": "I", "d": 18, "bf": 7.5, "tf": 0.57, "tw": 0.355}},
        "members": {"C1": {"i": "N1", "j": "N2", "section": "col", "material": "steel"},
                    "C2": {"i": "N4", "j": "N3", "section": "col", "material": "steel", "orientation": {"roll": 90}},
                    "B1": {"i": "N2", "j": "N3", "section": "beam", "material": "steel",
                           "releases": {"j": {"rz": {"spring": 50000}}}}},
        "loadCases": {"D": {"loads": [{"type": "selfWeight"},
                                      {"type": "memberDistributed", "member": "B1", "direction": "Y", "w1": -0.1,
                                       "w2": -0.05, "start": 0.2, "end": 0.9}]},
                      "W": {"loads": [{"type": "nodal", "node": "N2", "Fx": 8, "Fz": -3},
                                      {"type": "memberDistributed", "member": "C1", "direction": "y", "w1": 0.03}]}},
        "combinations": {"L": {"factors": {"D": 1.2, "W": 1.6}}},
        "analysis": {"shearDeformation": False},
    }
    return doc, "both"


def m12_continuous_beam_with_combination_factors():
    nodes = {f"n{k}": (3000.0 * k, 0, 0) for k in range(4)}
    members = {f"s{k}": (f"n{k}", f"n{k + 1}") for k in range(3)}
    supports = {"n0": {"type": "pinned", "dofs": {"rx": "fixed"}}, "n1": "roller", "n2": "roller",
                "n3": {"type": "roller", "dofs": {"uy": "fixed"}}}
    cases = {"D": [{"type": "memberDistributed", "member": m, "direction": "Z", "w1": -6.0} for m in members],
             "L1": [{"type": "memberDistributed", "member": "s0", "direction": "Z", "w1": -9.0, "w2": -2.0},
                    {"type": "memberPoint", "member": "s2", "direction": "Z", "value": -3e4, "position": 0.35}],
             "T": [{"type": "nodal", "node": "n1", "Fy": 4e3, "Mx": 2e5}]}
    doc = document(nodes, members, supports, cases, combos={"L": {"factors": {"D": 1.35, "L1": 1.5, "T": -0.7}}})
    return doc, "both"


def m13_grid_out_of_plane():
    nodes = {"a": (0, 0, 0), "b": (4000, 0, 0), "c": (4000, 3000, 0), "d": (0, 3000, 0), "e": (2000, 1500, 0)}
    members = {"ab": ("a", "b"), "bc": {"i": "b", "j": "c", "orientation": {"roll": 90}}, "dc": ("d", "c"),
               "ad": ("a", "d"), "ae": {"i": "a", "j": "e", "orientation": {"roll": 30}},
               "ec": {"i": "e", "j": "c", "releases": {"i": {"rx": True}}}}
    loads = [{"type": "nodal", "node": "e", "Fz": -5e4, "Mx": 3e6},
             {"type": "memberDistributed", "member": "bc", "direction": "Z", "w1": -4.0, "w2": -9.0},
             {"type": "memberPoint", "member": "ab", "direction": "x", "kind": "moment", "value": 2e6, "position": 0.6}]
    return document(nodes, members, {"a": "fixed", "c": "pinned", "d": "roller"}, {"L": loads}), "both"


def m14_two_storey_space_frame():
    nodes, members, supports = {}, {}, {}
    plan = {"1": (0, 0), "2": (6000, 0), "3": (6000, 5000), "4": (0, 5000)}
    for lvl, z in enumerate((0, 3500, 7000)):
        for k, (x, y) in plan.items():
            nodes[f"{k}{lvl}"] = (x, y, z)
    sections = {"S": {"type": "properties", **SEC},
                "G": {"type": "properties", "A": 8000.0, "Iy": 2.0e7, "Iz": 1.2e8, "J": 6.0e5}}
    for k in plan:
        supports[f"{k}0"] = "fixed" if k in "13" else "pinned"
        for lvl in (0, 1):
            members[f"c{k}{lvl}"] = {"i": f"{k}{lvl}", "j": f"{k}{lvl + 1}",
                                     "orientation": {"roll": 90} if k in "24" else {}}
    for lvl in (1, 2):
        for a, b in (("1", "2"), ("2", "3"), ("4", "3"), ("1", "4")):
            members[f"g{a}{b}{lvl}"] = {"i": f"{a}{lvl}", "j": f"{b}{lvl}", "section": "G"}
    members["brace"] = {"i": "10", "j": "21", "type": "truss"}
    members["g121"]["releases"] = {"i": {"rz": True}, "j": {"rz": {"spring": 8e9}}}
    cases = {"D": [{"type": "selfWeight"}] + [{"type": "memberDistributed", "member": m, "direction": "Z", "w1": -12.0}
                                             for m in members if m.startswith("g")],
             "W": [{"type": "nodal", "node": "12", "Fx": 3e4, "Fy": -1e4}, {"type": "nodal", "node": "11", "Fx": 2e4},
                   {"type": "memberDistributed", "member": "c10", "direction": "X", "w1": 2.0, "w2": 5.0}]}
    doc = document(nodes, members, supports, cases, sections=sections,
                   combos={"L": {"factors": {"D": 1.2, "W": 1.0}}, "S": {"factors": {"D": 1.0, "W": 0.6},
                                                                         "purpose": "serviceability"}})
    return doc, "both"


MODELS = [m01_all_load_types, m02_rolled_member_full_releases, m03_partial_rotational_springs,
          m04_offsets_and_releases, m05_axial_spring_and_torsion_release, m06_support_springs,
          m07_truss_and_frame, m08_rigid_diaphragm, m09_rigid_link,
          m10_shear_deformation_with_releases_and_offsets, m11_y_up_kip_inch,
          m12_continuous_beam_with_combination_factors, m13_grid_out_of_plane, m14_two_storey_space_frame]


@pytest.mark.parametrize("builder", MODELS, ids=lambda f: f.__name__)
def test_three_way(builder):
    doc, expect = builder()
    model = build_model(doc)
    results = solve_pkg.analyze(model)
    foreign = {"pynite": solve("pynite", model, "L")}
    rep = cross_check(model, "L", benchfea_result=results, foreign=foreign)
    for s in ("opensees", "pynite"):
        assert rep.results[s].ok, rep.results[s].error
    assert rep.eligible == (["opensees", "pynite"] if expect == "both" else ["opensees"]), rep.summary()
    assert rep.cross_checked == (expect == "both"), rep.summary()
    # benchfea agrees with every solver that represents the model, the solvers agree with each other
    assert rep.comparisons and all(c.passed for c in rep.comparisons), rep.summary()
    for c in rep.comparisons:
        exact = all(rep.results[x].fidelity.level == EXACT for x in (c.a, c.b))
        if exact:
            assert c.max_relative < 1e-8, c.summary()        # same mathematical problem: round-off only
    # hand check that needs no solver: applied loads + reactions = 0
    for s, eq in rep.equilibrium.items():
        if rep.results[s].fidelity.level == EXACT:
            assert eq["force"] < 1e-9 and eq["moment"] < 1e-9, (s, eq)
    if expect != "both":
        # the solver that cannot represent the model really differs (the refusal is not a formality)
        bad = compare(from_benchfea(model, results, "L"), rep.results["pynite"], model)
        assert bad.max_relative > 1e-5


def test_report_structure_and_other_combination():
    doc, _ = m14_two_storey_space_frame()
    model = build_model(doc)
    results = solve_pkg.analyze(model)
    rep = cross_check(model, "S", benchfea_result=results, foreign={"pynite": solve("pynite", model, "S")})
    assert rep.cross_checked
    d = rep.to_dict()
    assert d["crossChecked"] and set(d["results"]) == {"opensees", "pynite", "benchfea"}
    assert {(c["a"], c["b"]) for c in d["comparisons"]} == {("benchfea", "opensees"), ("benchfea", "pynite"),
                                                           ("opensees", "pynite")}
    for c in d["comparisons"]:
        assert set(c["worst"]) == {"displacement", "rotation", "reactionForce", "reactionMoment", "endForce",
                                   "endMoment"}
    assert "CROSS-CHECKED" in rep.summary()
    import json
    json.dumps(d)                                             # plain JSON data


def test_a_wrong_engine_result_is_detected():
    """Perturb one benchfea member end moment by 1e-4 of the moment scale: not cross-checked, and the
    report names the member."""
    doc, _ = m01_all_load_types()
    model = build_model(doc)
    results = solve_pkg.analyze(model)
    comb = results.combinations["L"]
    scale = max(np.abs(m.end_forces[[3, 4, 5, 9, 10, 11]]).max() for m in comb.members.values())
    comb.members["b"].end_forces[11] += 1e-4 * scale
    rep = cross_check(model, "L", benchfea_result=results, foreign={"pynite": solve("pynite", model, "L")})
    assert not rep.cross_checked
    assert any("benchfea" in r and "member b end j Mz" in r for r in rep.reasons), rep.reasons


def test_second_order_three_way_when_the_engine_supports_it():
    if "secondOrder" not in solve_pkg.ANALYSIS_TYPES:
        pytest.skip("benchfea.solve has no second-order analysis yet (activates automatically)")
    h, b = 4000.0, 6000.0
    doc = document({"a": (0, 0, 0), "A": (0, 0, h), "B": (b, 0, h), "b": (b, 0, 0)},
                   {"c1": ("a", "A"), "g": ("A", "B"), "c2": ("b", "B")},
                   {"a": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}},
                    "b": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}},
                    # The frame is braced out of plane at the column tops.  Without these two
                    # restraints the columns are out-of-plane cantilevers and the 300 kN column loads
                    # exceed the elastic critical load (buckling factor 0.709): the engine then refuses
                    # the combination, correctly.  See tests/independent/test_indep_second_order.py
                    # (test_refusal_beyond_the_critical_load_reports_the_buckling_factor_of_the_reference).
                    "A": {"dofs": {"uy": "fixed"}}, "B": {"dofs": {"uy": "fixed"}}},
                   {"L": [{"type": "nodal", "node": "A", "Fz": -3e5, "Fx": 2e4}, {"type": "nodal", "node": "B", "Fz": -3e5},
                          {"type": "memberDistributed", "member": "g", "direction": "Z", "w1": -10.0}]},
                   combos={"L": {"factors": {"L": 1.0}, "analysisType": "secondOrder"}}, analysis={"alpha": 1.0})
    model = build_model(doc)
    results = solve_pkg.analyze(model)
    rep = cross_check(model, "L", benchfea_result=results, foreign={"pynite": solve("pynite", model, "L")})
    # end shears of a second-order analysis are compared along the UNDEFORMED member axes; if the
    # engine reports them otherwise this fails and docs/CROSSCHECK.md ("open questions") applies
    assert rep.cross_checked, rep.summary()
