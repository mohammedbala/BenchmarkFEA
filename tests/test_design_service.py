"""Serviceability: deflection limits of the user, relative to the chord between the
design-member ends (AISC 360-22 Chapter L sets no numeric limit; pp. 16.1-174, -175).

Tier F: (a) hand calculation with the beam deflection derived in the test from the
beam equation; (b) AISC Design Examples v16.0 F.8A (deflection of a simply
supported beam against L/240; tests/design_examples/part2, the ONE tolerance
policy), solved with the engine's own solver.
"""
import math

import numpy as np
import pytest

from benchfea import design
from benchfea.design.aisc360_22 import serviceability as sv
from benchfea.model import build_model
from benchfea.results import CombinationResult, MemberResult, Results
from design_examples.part2 import FT, assert_printed, example

E, I = 200000.0, 2.0e8
L = 8000.0


def _data(limits, members=None, combos=None, nodes=None, supports=None, design_extra=None):
    nodes = nodes or {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": L, "y": 0, "z": 0}}
    members = members or {"beam": {"i": "A", "j": "B", "section": "I400", "material": "S"}}
    dm = {} if limits is None else {"deflectionLimits": limits}
    dm.update(design_extra or {})
    return {
        "formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"}, "conventions": {"verticalAxis": "Z"},
        "nodes": nodes,
        "supports": supports or {"A": {"type": "pinned", "dofs": {"rx": "fixed"}},
                                 "B": {"type": "roller", "dofs": {"uy": "fixed"}}},
        "materials": {"S": {"E": E, "nu": 0.3, "Fy": 345, "Fu": 450}},
        "sections": {"I400": {"type": "I", "d": 400, "bf": 180, "tf": 14, "tw": 9, "r": 20}},
        "members": members,
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fx": 1}]},
                      "Lv": {"loads": [{"type": "nodal", "node": "B", "Fx": 1}]}},
        "combinations": combos or {"U1": {"factors": {"D": 1.4}},
                                   "S1": {"factors": {"D": 1.0}, "purpose": "serviceability"},
                                   "S2": {"factors": {"D": 1.0, "Lv": 1.0}, "purpose": "serviceability"}},
        "design": {"members": {"beam": dm}},
    }


def _udl_deflection(x, w, span=L):
    """Simply supported beam, uniform load w (derived from E I v'''' = w with v = v'' = 0 at both ends):
    v(x) = w x (span^3 - 2 span x^2 + x^3) / (24 E I); at midspan 5 w span^4 / (384 E I)."""
    return w * x * (span ** 3 - 2.0 * span * x ** 2 + x ** 3) / (24.0 * E * I)


def _results(model, shapes, converged=True):
    """Hand-made results: ``shapes[cid][member id]`` = (x, uy, uz) station displacements (local axes)."""
    combos = {}
    for cid, per_member in shapes.items():
        mr = {}
        for mid, (x, uy, uz) in per_member.items():
            x = np.asarray(x, float)
            u = np.column_stack([np.zeros(len(x)), np.asarray(uy, float), np.asarray(uz, float)])
            mr[mid] = MemberResult(np.zeros(12), x, np.zeros((len(x), 6)), u)
        n = len(model.nodes)
        combos[cid] = CombinationResult(cid, "firstOrder", np.zeros((n, 6)), np.zeros((n, 6)), mr, converged=converged)
    return Results("hand", combos)


def _by_name(checks):
    return {c.name: c for c in checks}


X = np.linspace(0.0, L, 11)
ZERO = np.zeros(11)


# ============================================================================== definition
def test_chord_deflection_removes_rigid_body_motion():
    x = np.linspace(0.0, 4000.0, 9)
    bend = np.column_stack([np.zeros(9), -3.0 * np.sin(math.pi * x / 4000.0), 2.0 * np.sin(math.pi * x / 4000.0)])
    # translation and rotation of the chord (support settlement, sway): u = a + b x in every component
    rigid = np.column_stack([5.0 + 0.0 * x, 7.0 - 0.004 * x, -2.0 + 0.003 * x])
    d = sv.chord_deflection(x, bend + rigid)
    assert d[:, 1:] == pytest.approx(bend[:, 1:], abs=1e-12) and np.all(d[0] == 0.0) and d[-1] == pytest.approx([0, 0, 0], abs=1e-12)
    assert sv.deflection_value(x, bend + rigid, "y") == (pytest.approx(3.0), 2000.0)
    assert sv.deflection_value(x, bend + rigid, "z")[0] == pytest.approx(2.0)
    assert sv.deflection_value(x, bend + rigid, "resultant")[0] == pytest.approx(math.hypot(3.0, 2.0))
    assert sv.deflection_value(x, rigid, "resultant")[0] == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError):
        sv.deflection_value(x, bend, "x")


# ============================================================================== tier F (a): hand calculation
def test_hand_deflection_against_span_ratio_and_absolute_limit():
    # w = 20 N/mm on 8 m, E I = 4e13: midspan deflection 5 w L^4/(384 E I) = 5*20*8000^4/(384*4e13) = 26.667 mm
    delta = 5.0 * 20.0 * 8000.0 ** 4 / (384.0 * 200000.0 * 2.0e8)
    assert delta == pytest.approx(26.6667, abs=1e-4) and _udl_deflection(4000.0, 20.0) == pytest.approx(delta, rel=1e-14)
    m = build_model(_data([{"spanRatio": 250, "direction": "y"}, {"spanRatio": 360}, {"absolute": 20.0, "direction": "z"},
                           {"spanRatio": 250, "absolute": 30.0, "direction": "y", "combinations": ["S2"]}]))
    # S1: gravity deflection (local y, downwards) on top of a support settlement of 4 mm at B (rigid-body part);
    # S2: 0.8 x that, plus 15 mm sideways (local z)
    uy = -_udl_deflection(X, 20.0) - 4.0 * X / L
    res = _results(m, {"S1": {"beam": (X, uy, ZERO)},
                       "S2": {"beam": (X, 0.8 * uy, 15.0 * np.sin(math.pi * X / L))},
                       "U1": {"beam": (X, 10.0 * uy, ZERO)}})
    c1, c2, c3, c4 = sv.deflection_checks(m, res, "beam")
    # limit 1: L/250 = 32 mm, direction y: governing combination S1 (26.667 > 0.8 * 26.667): ratio 0.8333
    assert c1.name == "deflection limit 1 (local y)" and (c1.status, c1.advisory) == ("ok", False)
    assert c1.capacity == 32.0 and c1.demand == pytest.approx(delta, rel=1e-13) and c1.ratio == pytest.approx(delta / 32.0)
    assert (c1.combination, c1.x, c1.passed, c1.dimension) == ("S1", 4000.0, True, "length")
    assert c1.clause.startswith("AISC 360-22 L1, L2") and "no numeric limit" in c1.clause
    assert c1.inputs["spanRatio"] == (250, "dimensionless", "member") and c1.inputs["L"][0] == 8000.0
    assert c1.inputs["combinations"][0] == "S1, S2"                         # default: every serviceability combination
    assert c1.details["L / deflection"][0] == pytest.approx(8000.0 / delta)
    # limit 2: L/360 = 22.22 mm, resultant: S2 governs: sqrt((0.8*26.667)^2 + 15^2) = 26.08 < S1 = 26.667 -> S1
    assert c2.name == "deflection limit 2 (resultant)" and c2.capacity == pytest.approx(8000.0 / 360.0)
    assert math.hypot(0.8 * delta, 15.0) < delta
    assert c2.combination == "S1" and c2.ratio == pytest.approx(delta / (8000.0 / 360.0)) and c2.passed is False
    # limit 3: 20 mm absolute, direction z: only S2 moves sideways: 15/20
    assert (c3.combination, c3.capacity) == ("S2", 20.0) and c3.ratio == pytest.approx(0.75) and "spanRatio" not in c3.inputs
    # limit 4: both given -> the smaller governs: min(32, 30) = 30; listed combination S2 only: 0.8 * 26.667 / 30
    assert c4.capacity == 30.0 and c4.combination == "S2" and c4.ratio == pytest.approx(0.8 * delta / 30.0)
    assert c4.inputs["combinations"][0] == "S2" and c4.inputs["absolute"][0] == 30.0
    # the strength combination U1 (10 x the deflection) is never used for a deflection limit
    assert all(c.combination != "U1" for c in (c1, c2, c3, c4))
    # through design.check: an exceeded limit fails the member (it is not advisory)
    out = design.check(m, res)
    names = _by_name(out.members["beam"].checks)
    assert names["deflection limit 2 (resultant)"].passed is False
    assert out.members["beam"].passed is False and out.members["beam"].status == "fail" and out.status == "fail"
    doc = out.to_document(m)
    d = {c["name"]: c for c in doc["members"]["beam"]["checks"]}["deflection limit 1 (local y)"]
    assert d["unit"] == "mm" and d["capacity"] == 32.0 and d["combination"] == "S1"


def test_chord_is_between_the_design_member_ends_not_the_analysis_nodes():
    # design member of two analysis members (A-M, M-B): the deflection at the intermediate node counts
    nodes = {"A": {"x": 0, "y": 0, "z": 0}, "M": {"x": 5000, "y": 0, "z": 0}, "B": {"x": L, "y": 0, "z": 0}}
    members = {"m1": {"i": "A", "j": "M", "section": "I400", "material": "S"},
               "m2": {"i": "M", "j": "B", "section": "I400", "material": "S"}}
    data = _data([{"spanRatio": 300, "direction": "y"}], members=members, nodes=nodes)
    data["design"] = {"members": {"girder": {"members": ["m1", "m2"], "deflectionLimits": [{"spanRatio": 300, "direction": "y"}]}}}
    m = build_model(data)
    x1, x2 = np.linspace(0.0, 5000.0, 6), np.linspace(0.0, 3000.0, 4)
    res = _results(m, {cid: {"m1": (x1, -_udl_deflection(x1, 20.0), np.zeros(6)),
                             "m2": (x2, -_udl_deflection(5000.0 + x2, 20.0), np.zeros(4))} for cid in ("S1", "S2", "U1")})
    c, = sv.deflection_checks(m, res, "girder")
    assert c.status == "ok" and c.inputs["L"][0] == 8000.0 and c.capacity == pytest.approx(8000.0 / 300.0)
    # stations at 0, 1000 ... 8000: the maximum is at x = 4000 (inside m1)
    assert c.demand == pytest.approx(_udl_deflection(4000.0, 20.0), rel=1e-13) and c.x == 4000.0
    x, u = sv.member_displacements(m, res.combinations["S1"], "girder")
    assert list(x) == [0, 1000, 2000, 3000, 4000, 5000, 5000, 6000, 7000, 8000]


def test_cases_without_a_result_are_not_checked_with_the_reason():
    uy = -_udl_deflection(X, 20.0)
    # no deflection limit: no check at all (reported as a default of the model)
    m = build_model(_data(None))
    assert sv.deflection_checks(m, _results(m, {"S1": {"beam": (X, uy, ZERO)}}), "beam") == []
    assert any(d.path.endswith("/deflectionLimits") for d in m.defaults_used)
    # a listed combination is missing from the results
    m = build_model(_data([{"spanRatio": 300}]))
    c, = sv.deflection_checks(m, _results(m, {"S1": {"beam": (X, uy, ZERO)}}), "beam")
    assert c.status == "notChecked" and "S2" in c.reason and "not in the analysis results" in c.reason
    # a combination that did not converge
    res = _results(m, {"S1": {"beam": (X, uy, ZERO)}, "S2": {"beam": (X, uy, ZERO)}}, converged=False)
    c, = sv.deflection_checks(m, res, "beam")
    assert c.status == "notChecked" and "did not converge" in c.reason
    # no serviceability combination in the model and none listed
    m = build_model(_data([{"spanRatio": 300}], combos={"U1": {"factors": {"D": 1.4}}}))
    c, = sv.deflection_checks(m, _results(m, {"U1": {"beam": (X, uy, ZERO)}}), "beam")
    assert c.status == "notChecked" and "no combination" in c.reason
    # a member with a notChecked deflection limit is not "pass"
    out = design.check(m, _results(m, {"U1": {"beam": (X, uy, ZERO)}}))
    assert out.members["beam"].status == "incomplete" and any("deflection limit 1" in s for s in out.members["beam"].not_checked)


def test_cantilever_deflection_is_not_checked():
    # free end at B: the chord through both ends does not measure the tip deflection -> notChecked
    m = build_model(_data([{"spanRatio": 180}], supports={"A": {"type": "fixed"}}))
    tip = 50.0 * (X / L) ** 2                      # tip deflection 50 mm; relative to the chord it would be 12.5 mm
    res = _results(m, {"S1": {"beam": (X, -tip, ZERO)}, "S2": {"beam": (X, -tip, ZERO)}})
    assert sv.deflection_value(X, np.column_stack([ZERO, -tip, ZERO]), "y")[0] == pytest.approx(12.5)
    c, = sv.deflection_checks(m, res, "beam")
    assert c.status == "notChecked" and "free end" in c.reason and "tip deflection" in c.reason


def test_reduced_stiffness_combination_is_flagged():
    combos = {"U1": {"factors": {"D": 1.4}},
              "S1": {"factors": {"D": 1.0}, "purpose": "serviceability", "analysisType": "secondOrder",
                     "stiffnessReduction": True}}
    m = build_model(_data([{"spanRatio": 300}], combos=combos))
    c, = sv.deflection_checks(m, _results(m, {"S1": {"beam": (X, -_udl_deflection(X, 20.0), ZERO)}}), "beam")
    assert c.status == "ok" and any("L1 User Note" in n and "S1" in n for n in c.notes)


# ============================================================================== tier F (b): Design Examples v16.0
def test_design_example_f8a_live_load_deflection_with_the_solver():
    """Example F.8A (F-37 to F-38, PDF 199-200): simply supported 21 ft beam, I = 54.4 in.^4, live load
    0.375 kip/ft, limit L/240: the beam is analysed by the engine (no shear deformation, as in the example)
    and the deflection check reproduces the printed deflection and limit."""
    solve = pytest.importorskip("benchfea.solve")
    ex = example("F.8A")
    i, p = ex["inputs"], ex["printed"]
    span = i["span_ft"] * FT
    data = {
        "formatVersion": "1.0.0", "units": {"force": "kip", "length": "in"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": span, "y": 0, "z": 0}},
        "supports": {"A": {"type": "pinned", "dofs": {"rx": "fixed"}}, "B": {"type": "roller", "dofs": {"uy": "fixed"}}},
        "materials": {"S": {"E": i["E"], "nu": 0.3, "Fy": 50}},
        # only the listed moment of inertia matters for the deflection; the other values are placeholders
        "sections": {"sec": {"type": "properties", "A": 5.0, "Iy": i["I"], "Iz": i["I"], "J": 80.0}},
        "members": {"beam": {"i": "A", "j": "B", "section": "sec", "material": "S"}},
        "loadCases": {"L": {"loads": [{"type": "memberDistributed", "member": "beam", "direction": "Z",
                                       "w1": -i["wL_kip_per_ft"] / FT}]}},
        "combinations": {"live": {"factors": {"L": 1.0}, "purpose": "serviceability"}},
        "analysis": {"shearDeformation": False},
        "design": {"members": {"beam": {"deflectionLimits": [{"spanRatio": i["spanRatio"], "direction": "y"}]}}},
    }
    m = build_model(data)
    res = solve.analyze(m)
    c, = sv.deflection_checks(m, res, "beam")
    kip_in = m.units
    assert c.status == "ok" and c.combination == "live"
    assert_printed(float(kip_in.from_internal(c.capacity, "length")), p["delta_max"], "L/240")
    assert_printed(float(kip_in.from_internal(c.demand, "length")), p["delta"], "live load deflection")
    assert c.passed and c.x == pytest.approx(m.design_members["beam"].length / 2.0)
    # and against the closed form with the listed numbers
    closed = 5.0 * (i["wL_kip_per_ft"] / FT) * span ** 4 / (384.0 * i["E"] * i["I"])
    assert float(kip_in.from_internal(c.demand, "length")) == pytest.approx(closed, rel=1e-9)
