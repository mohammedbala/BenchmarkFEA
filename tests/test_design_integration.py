"""Design orchestration end to end: model -> analysis (the engine's solver) -> design.check ->
design document and CLI.

* a simple beam and a beam-column analysed by ``benchfea.solve`` and checked against
  independent arithmetic (AISC 360-22 F2, E3, E4, H1 written out in the test from the
  model's section properties);
* the rule "unbraced lengths, K factors and every other design variable affect design
  capacities only and never an analysis result";
* ``design.check`` never mutates the model or the results;
* status rules ("pass" is never reported when something applicable was not checked;
  first-order required strengths give "conditional"), analysis-method / K-factor warnings;
* the ``benchfea design`` command.

The solver-based tests need ``benchfea.solve``; the second-order ones activate by
themselves once the solver registers the "secondOrder" analysis type.
"""
import copy
import json
import math

import numpy as np
import pytest

from benchfea import design
from benchfea.cli import main
from benchfea.design.aisc360_22 import e_compression, f_flexure as ff, h_combined as hc
from benchfea.model import build_model
from benchfea.results import CombinationResult, MemberResult, Results
from helpers import assert_same

solve = pytest.importorskip("benchfea.solve")
SECOND_ORDER = "secondOrder" in getattr(solve.analysis, "ANALYSIS_TYPES", {})
needs_second_order = pytest.mark.skipif(not SECOND_ORDER, reason="the solver has not registered second-order analysis")

E, NU, FY, FU = 200000.0, 0.3, 345.0, 450.0
G = E / (2.0 * (1.0 + NU))
PI = math.pi
SECTION = {"type": "I", "d": 400, "bf": 180, "tf": 14, "tw": 9, "r": 20}       # made-up rolled I (fillet r > 0)


def _beam(design_block=None, combos=None, w=-20.0, extra=None):
    """Simply supported 6 m beam along X (y local = +Z), uniform load w along Z (N/mm)."""
    data = {
        "formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": 6000, "y": 0, "z": 0}},
        "supports": {"A": {"type": "pinned", "dofs": {"rx": "fixed"}}, "B": {"type": "roller", "dofs": {"uy": "fixed"}}},
        "materials": {"S": {"E": E, "nu": NU, "Fy": FY, "Fu": FU}},
        "sections": {"I400": dict(SECTION)},
        "members": {"beam": {"i": "A", "j": "B", "section": "I400", "material": "S"}},
        "loadCases": {"D": {"loads": [{"type": "memberDistributed", "member": "beam", "direction": "Z", "w1": w}]}},
        "combinations": combos or {"ULS": {"factors": {"D": 1.5}},
                                   "SLS": {"factors": {"D": 1.0}, "purpose": "serviceability"}},
    }
    if design_block is not None:
        data["design"] = design_block
    data.update(extra or {})
    return data


def _column(design_block=None, P=-400.0e3, w=15.0, combos=None):
    """Pin-ended 4 m column along Z (y local = +X), axial load P at the top and a uniform load w along X."""
    data = {
        "formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": 0, "y": 0, "z": 4000}},
        "supports": {"A": {"type": "pinned", "dofs": {"rz": "fixed"}},
                     "B": {"dofs": {"ux": "fixed", "uy": "fixed"}}},
        "materials": {"S": {"E": E, "nu": NU, "Fy": FY, "Fu": FU}},
        "sections": {"I400": dict(SECTION)},
        "members": {"col": {"i": "A", "j": "B", "section": "I400", "material": "S"}},
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fz": P},
                                      {"type": "memberDistributed", "member": "col", "direction": "X", "w1": w}]}},
        "combinations": combos or {"U1": {"factors": {"D": 1.0}}},
    }
    if design_block is not None:
        data["design"] = design_block
    return data


def _by_name(checks):
    return {c.name: c for c in checks}


def _arrays(results):
    """Every numeric array of the analysis results, by path."""
    out = {}
    for cid, c in results.combinations.items():
        out[f"{cid}/displacements"], out[f"{cid}/reactions"] = c.displacements, c.reactions
        for mid, m in c.members.items():
            out[f"{cid}/{mid}/end_forces"], out[f"{cid}/{mid}/x"] = m.end_forces, m.x
            out[f"{cid}/{mid}/forces"], out[f"{cid}/{mid}/displacements"] = m.forces, m.displacements
    return out


def _bitwise_equal(a, b):
    A, B = _arrays(a), _arrays(b)
    assert list(A) == list(B)
    for k in A:
        assert A[k].shape == B[k].shape and A[k].tobytes() == B[k].tobytes(), k


# ------------------------------------------------------------------------------ independent AISC arithmetic
def _f2_strength(sec, Lb, Cb):
    """phi Mn of F2 (yielding and LTB) from the model's section properties (x = local z)."""
    Zx, Sx, Iy, Cw, Jx, A = sec.Zz, sec.Sz, sec.Iy, sec.Cw, sec.J, sec.A
    Mp = FY * Zx                                                            # Eq. F2-1
    ry = math.sqrt(Iy / A)
    Lp = 1.76 * ry * math.sqrt(E / FY)                                      # Eq. F2-5
    rts = math.sqrt(math.sqrt(Iy * Cw) / Sx)                                # Eq. F2-7
    ho = SECTION["d"] - SECTION["tf"]
    X = Jx / (Sx * ho)                                                      # c = 1 (Eq. F2-8a)
    Lr = 1.95 * rts * E / (0.7 * FY) * math.sqrt(X + math.sqrt(X * X + 6.76 * (0.7 * FY / E) ** 2))     # Eq. F2-6
    if Lb <= Lp:
        Mn = Mp
    elif Lb <= Lr:
        Mn = min(Mp, Cb * (Mp - (Mp - 0.7 * FY * Sx) * (Lb - Lp) / (Lr - Lp)))                          # Eq. F2-2
    else:
        Fcr = Cb * PI ** 2 * E / (Lb / rts) ** 2 * math.sqrt(1.0 + 0.078 * X * (Lb / rts) ** 2)         # Eq. F2-4
        Mn = min(Mp, Fcr * Sx)                                                                          # Eq. F2-3
    return 0.9 * Mn, 0.9 * Mp, Lp, Lr


def _e_strength(sec, Lcx, Lcy, Lcz):
    """phi Pn of E3 (both axes) and E4 (torsional) for a nonslender section -> the smallest."""
    A, Ix, Iy = sec.A, sec.Iz, sec.Iy
    out = []
    for Lc, I in ((Lcx, Ix), (Lcy, Iy)):
        slender = Lc / math.sqrt(I / A)
        Fe = PI ** 2 * E / slender ** 2                                     # Eq. E3-4
        Fn = 0.658 ** (FY / Fe) * FY if slender <= 4.71 * math.sqrt(E / FY) else 0.877 * Fe      # E3-2 / E3-3
        out.append(0.9 * Fn * A)
    Fe = (PI ** 2 * E * sec.Cw / Lcz ** 2 + G * sec.J) / (Ix + Iy)          # Eq. E4-2
    Fn = 0.658 ** (FY / Fe) * FY if FY / Fe <= 2.25 else 0.877 * Fe
    out.append(0.9 * Fn * A)
    return min(out), out


# ============================================================================== simple beam
def test_simple_beam_analysed_and_checked():
    m = build_model(_beam({"members": {"beam": {"bracing": {"topFlange": [0.5]},
                                               "deflectionLimits": [{"spanRatio": 360, "direction": "y"}]}}}))
    res = solve.analyze(m)
    out = design.check(m, res)
    beam = out.members["beam"]
    sec = m.sections["I400"]
    checks = _by_name(beam.checks)
    # demands from the stations: Mu = 1.5 * 20 * 6000^2 / 8 = 135e6 N mm at midspan; Vu = 1.5 * 20 * 3000 = 90e3
    y = checks[ff.NAME_Y["x"]]
    assert y.demand == pytest.approx(135.0e6, rel=1e-9) and y.x == pytest.approx(3000.0) and y.combination == "ULS"
    assert checks["shear along local y"].demand == pytest.approx(90.0e3, rel=1e-9)
    # F2: gravity load puts the top (+y) flange in compression; braced at midspan: Lb = 3000, Cb by Eq. F1-1 on
    # half a parabola: 12.5 / (2.5 + 3*0.4375 + 4*0.75 + 3*0.9375) = 1.2987
    Cb = 12.5 / (2.5 + 3 * 0.4375 + 4 * 0.75 + 3 * 0.9375)
    cap, cap_y, Lp, Lr = _f2_strength(sec, 3000.0, Cb)
    assert Lp < 3000.0 < Lr
    ltb = checks[ff.NAME_LTB]
    assert ltb.inputs["Cb"][0] == pytest.approx(Cb, rel=1e-9) and ltb.inputs["Lb"][0] == pytest.approx(3000.0)
    assert ltb.capacity == pytest.approx(cap, rel=1e-9) and y.capacity == pytest.approx(cap_y, rel=1e-12)
    assert ltb.ratio == pytest.approx(135.0e6 / cap, rel=1e-9)
    # without the brace the whole span is unbraced: Cb = 12.5/11, Lb = 6000 > Lr -> elastic LTB governs the member
    m2 = build_model(_beam({"members": {"beam": {"deflectionLimits": [{"spanRatio": 360, "direction": "y"}]}}}))
    out2 = design.check(m2, solve.analyze(m2))
    cap2, _, _, Lr2 = _f2_strength(sec, 6000.0, 12.5 / 11.0)
    ltb2 = _by_name(out2.members["beam"].checks)[ff.NAME_LTB]
    assert 6000.0 > Lr2 and ltb2.capacity == pytest.approx(cap2, rel=1e-9) and cap2 < cap
    assert out2.members["beam"].governing is ltb2 and ltb2.inputs["Lb"][2] == "default"
    # deflection: service load 20 N/mm; 5 w L^4/(384 E I) plus the shear part w L^2/(8 G As) (shear deformation
    # is on by default: integral of V/(G As) from the support to midspan); limit 6000/360 = 16.667 mm
    c = checks["deflection limit 1 (local y)"]
    delta = 5.0 * 20.0 * 6000.0 ** 4 / (384.0 * E * sec.Iz)
    if m.analysis.shear_deformation and sec.Asy:
        delta += 20.0 * 6000.0 ** 2 / (8.0 * G * sec.Asy)
    assert c.demand == pytest.approx(delta, rel=1e-9) and c.capacity == pytest.approx(6000.0 / 360.0) and c.combination == "SLS"
    # no axial force, no minor-axis moment, no torsion: no interaction and no torsion record
    assert hc.NAME_H1["c"] not in checks and hc.NAME_TORSION not in checks
    # status: every check passes, but the strength combination is first-order -> "conditional", never "pass"
    assert beam.passed is True and all(c.passed for c in beam.checks if c.status == "ok" and not c.advisory)
    assert beam.status == "conditional" and out.status == "conditional" and not beam.not_checked
    assert any("first-order" in w and "C1" in w for w in beam.warnings) and any("first-order" in w for w in out.warnings)
    assert out.combinations["ULS"] == {"purpose": "strength", "inResults": True, "analysisType": "firstOrder",
                                       "secondOrder": False, "converged": True, "stiffnessReduction": False,
                                       "notionalLoad": False, "alpha": 1.0, "directAnalysisSettings": False}
    # every check of every strength combination is in the result and in the document
    assert set(beam.by_combination) == {"ULS"} and len(beam.by_combination["ULS"]) >= 12
    doc = out.to_document(m)
    text = json.dumps(doc)
    assert "NaN" not in text and "Infinity" not in text
    assert doc["status"] == "conditional" and doc["members"]["beam"]["status"] == "conditional"
    assert doc["combinations"]["ULS"]["secondOrder"] is False
    names = [c["name"] for c in doc["members"]["beam"]["byCombination"]["ULS"]]
    assert ff.NAME_LTB in names and "shear along local y" in names
    d = {c["name"]: c for c in doc["members"]["beam"]["checks"]}[ff.NAME_LTB]
    assert d["unit"] == "N*mm"
    assert d["inputs"]["Cb"]["source"] == "derived: Eq. F1-1" and "F1-1" in d["clause"]
    assert doc["members"]["beam"]["governing"]["name"] == beam.governing.name
    for c in doc["members"]["beam"]["checks"]:
        assert c["clause"].startswith("AISC 360-22") and all(v["source"] for v in c["inputs"].values())
    # nothing in the document claims an approval
    assert not any(word in text.lower() for word in ("approved", "certified", "checked by"))


def test_design_member_with_defaults_reports_them():
    m = build_model(_beam())                    # no design block at all: every design variable is a default
    out = design.check(m, solve.analyze(m))
    ltb = _by_name(out.members["beam"].checks)[ff.NAME_LTB]
    assert ltb.inputs["Lb"] == (pytest.approx(6000.0), "length", "default") and ltb.inputs["method"][2] == "default"
    paths = {d.path for d in out.defaults_used}
    assert {"/design/members/beam/Lb", "/design/members/beam/Cb", "/design/members/beam/method"} <= paths
    doc = out.to_document(m)
    assert doc["members"]["beam"]["parameters"]["Lb"]["source"] == "default"
    assert any(d["path"] == "/design/members/beam/Kx" for d in doc["defaultsUsed"])


# ============================================================================== beam-column
def test_beam_column_interaction_from_the_solver():
    m = build_model(_column())
    res = solve.analyze(m)
    out = design.check(m, res)
    col = out.members["col"]
    sec = m.sections["I400"]
    checks = _by_name(col.checks)
    # first-order demands: Pr = 400e3 (constant), Mrx = w L^2/8 = 15 * 4000^2/8 = 30e6 at mid-height
    Pc, parts = _e_strength(sec, 4000.0, 4000.0, 4000.0)
    assert Pc == parts[1]                                                   # minor-axis flexural buckling governs
    assert checks[e_compression.NAME_FB["y"]].capacity == pytest.approx(Pc, rel=1e-9)
    assert checks[e_compression.NAME_FB["y"]].demand == pytest.approx(400.0e3, rel=1e-9)
    # the load along +X = +y local bends the column with the -y flange in compression (Mz < 0): one unbraced
    # segment, Cb = 12.5/11 (parabola), Lb = 4000
    Mcx, _, Lp, Lr = _f2_strength(sec, 4000.0, 12.5 / 11.0)
    assert checks[ff.NAME_LTB].capacity == pytest.approx(Mcx, rel=1e-9)
    assert checks[ff.NAME_LTB].demand == pytest.approx(30.0e6, rel=1e-9) and checks[ff.NAME_LTB].x == pytest.approx(2000.0)
    # H1.1 at mid-height: Pr/Pc >= 0.2 -> Eq. H1-1a
    assert 400.0e3 / Pc >= 0.2
    value = 400.0e3 / Pc + 8.0 / 9.0 * 30.0e6 / Mcx
    h = checks[hc.NAME_H1["c"]]
    assert h.ratio == pytest.approx(value, rel=1e-9) and h.x == pytest.approx(2000.0) and "Eq. H1-1a" in h.clause
    assert h.inputs["Pc"][2] == "derived: " + e_compression.NAME_FB["y"]
    assert h.inputs["Mcx"][2] == "derived: " + ff.NAME_LTB
    assert col.governing is h and col.passed == (value <= 1.0)
    # first-order analysis of a member in compression: said in the check, in the member and in the design
    assert any("FIRST-ORDER" in n for n in h.notes)
    assert col.status == ("conditional" if value <= 1.0 else "fail") and out.combinations["U1"]["secondOrder"] is False
    assert any("first-order" in w for w in col.warnings)
    # K = 1 is used although the analysis is not a direct analysis: warned, never silent
    assert any("K = 1" in w and "direct analysis" in w and "U1" in w for w in col.warnings)
    # the engine did not amplify: B1 of Appendix 8 for this member would be > 1
    B1 = hc.b1_a_8_3(1.0, 1.0, 400.0e3, hc.pe1_a_8_5(E * sec.Iz, 4000.0))
    assert B1 > 1.0 and h.details["Mrx"][0] == pytest.approx(30.0e6, rel=1e-9)


def test_k_factor_below_one_and_tension_only_member_warnings():
    m = build_model(_column({"members": {"col": {"Ky": 0.7}}}))
    col = design.check(m, solve.analyze(m)).members["col"]
    assert any("below 1" in w and "rational analysis" in w for w in col.warnings)
    # a member without compression gets no K-factor warning
    m = build_model(_column(P=+400.0e3))
    col = design.check(m, solve.analyze(m)).members["col"]
    assert not any("effective length factors" in w for w in col.warnings)
    assert hc.NAME_H1["t"] in _by_name(col.checks)


# ============================================================================== design never changes the analysis
def test_design_variables_change_capacities_but_never_analysis_results():
    """Unbraced lengths, K factors, bracing, Cb, method, material overrides: the analysis results are
    bitwise identical; the available strengths change."""
    base_model = build_model(_column())
    base = solve.analyze(base_model)
    base_checks = _by_name(design.check(base_model, base).members["col"].checks)
    variants = {
        "Ky": ({"members": {"col": {"Ky": 2.0}}}, e_compression.NAME_FB["y"]),
        "Kx": ({"members": {"col": {"Kx": 3.0}}}, e_compression.NAME_FB["x"]),
        "Kz": ({"defaults": {"Kz": 2.5}}, e_compression.NAME_TB),
        "Ly": ({"members": {"col": {"Ly": 2000}}}, e_compression.NAME_FB["y"]),
        "Lx": ({"groups": {"g": {"members": ["col"], "Lx": 9000}}}, e_compression.NAME_FB["x"]),
        "Lb": ({"members": {"col": {"Lb": 9000}}}, ff.NAME_LTB),
        "bracing": ({"members": {"col": {"bracing": {"bottomFlange": [0.5], "minor": [0.5]}}}}, ff.NAME_LTB),
        "Cb": ({"members": {"col": {"Cb": 1.0}}}, ff.NAME_LTB),
        "method": ({"defaults": {"method": "ASD"}}, ff.NAME_Y["x"]),
        "Fy": ({"members": {"col": {"Fy": 250}}}, ff.NAME_Y["x"]),
        "E": ({"members": {"col": {"E": 150000}}}, e_compression.NAME_FB["y"]),
    }
    for name, (block, changed) in variants.items():
        m = build_model(_column(block))
        r = solve.analyze(m)
        _bitwise_equal(base, r)                                             # the analysis did not notice
        checks = _by_name(design.check(m, r).members["col"].checks)
        assert checks[changed].capacity != base_checks[changed].capacity, name
        # ... and the demands of the checks are the same numbers
        assert checks[changed].demand == base_checks[changed].demand, name
    # the direction of the change, for two of them: K = 2 about the minor axis quarters Fe (Eq. E3-4)
    m = build_model(_column({"members": {"col": {"Ky": 2.0}}}))
    c = _by_name(design.check(m, solve.analyze(m)).members["col"].checks)[e_compression.NAME_FB["y"]]
    assert c.details["Fe"][0] == pytest.approx(base_checks[e_compression.NAME_FB["y"]].details["Fe"][0] / 4.0, rel=1e-12)
    assert c.capacity < base_checks[e_compression.NAME_FB["y"]].capacity
    # same for the beam: bracing and deflection limits
    b0 = solve.analyze(build_model(_beam()))
    b1 = solve.analyze(build_model(_beam({"defaults": {"Lb": 1500, "Kx": 2.0, "deflectionLimits": [{"absolute": 5}]},
                                          "members": {"beam": {"bracing": {"minor": [0.25, 0.75]}}}})))
    _bitwise_equal(b0, b1)


def test_check_never_mutates_model_or_results():
    for data in (_beam({"members": {"beam": {"bracing": {"topFlange": [0.5]}, "Ky": 0.9,
                                             "deflectionLimits": [{"spanRatio": 360}]}}}), _column()):
        m = build_model(data)
        r = solve.analyze(m)
        m_before, r_before = copy.deepcopy(m), copy.deepcopy(r)
        first = design.check(m, r)
        assert_same(m_before, m)                                            # the model is untouched ...
        _bitwise_equal(r_before, r)                                         # ... and so is every result array
        assert_same(r_before.combinations, r.combinations)
        assert [str(w) for w in r_before.warnings] == [str(w) for w in r.warnings]
        # checking again gives the same design; mutating the returned design does not reach model / results
        second = design.check(m, r)
        assert [(c.name, c.ratio, c.capacity) for mm in first.members.values() for c in mm.checks] == \
               [(c.name, c.ratio, c.capacity) for mm in second.members.values() for c in mm.checks]
        first.defaults_used.clear()
        first.warnings.append("x")
        for mm in first.members.values():
            for c in mm.checks:
                c.inputs.clear()
                c.details.clear()
        assert_same(m_before, m)
        _bitwise_equal(r_before, r)
        # the demands hold copies of the station arrays, not views into the results
        d = design.member_demands(m, r, next(iter(m.design_members)))
        cid = next(iter(d))
        d[cid].extra["forces"][:] = 0.0
        d[cid].extra["x"][:] = -1.0
        _bitwise_equal(r_before, r)


# ============================================================================== status rules on hand-built results
def _hand_results(model, analysis_type="firstOrder", n=-100.0e3, mz=20.0e6, t=0.0, converged=True, combos=("U1",)):
    x = np.linspace(0.0, 4000.0, 5)
    f = np.zeros((5, 6))
    f[:, 0], f[:, 3], f[:, 5] = n, t, mz * np.sin(PI * x / 4000.0)
    f[:, 1] = -mz * PI / 4000.0 * np.cos(PI * x / 4000.0)                   # dMz/dx = -Vy
    out = {}
    for cid in combos:
        mr = MemberResult(np.concatenate([-f[0], f[-1]]), x.copy(), f.copy(), np.zeros((5, 3)))
        out[cid] = CombinationResult(cid, analysis_type, np.zeros((2, 6)), np.zeros((2, 6)), {"col": mr},
                                     converged=converged, alpha=1.0)
    return Results("hand", out)


def test_status_pass_only_with_second_order_results_and_everything_checked():
    second = {"U1": {"factors": {"D": 1.0}, "analysisType": "secondOrder", "stiffnessReduction": True,
                     "notionalLoad": {"direction": [1, 0, 0]}}}
    m = build_model(_column(combos=second))
    # second-order, stiffness reduction, notional load, everything checked and passing -> "pass", no K warning
    out = design.check(m, _hand_results(m, "secondOrder"))
    col = out.members["col"]
    assert col.status == "pass" and out.status == "pass" and col.passed
    assert out.combinations["U1"]["directAnalysisSettings"] is True and out.warnings == []
    assert not any("effective length" in w or "first-order" in w for w in col.warnings)
    # the same results flagged first-order -> "conditional"
    out = design.check(m, _hand_results(m, "firstOrder"))
    assert out.status == "conditional" and out.members["col"].status == "conditional"
    # not converged -> "incomplete"
    out = design.check(m, _hand_results(m, "secondOrder", converged=False))
    assert out.status == "incomplete" and any("did not converge" in w for w in out.warnings)
    # a torsional moment that is not negligible (H3 not implemented) -> "incomplete", never "pass"
    out = design.check(m, _hand_results(m, "secondOrder", t=5.0e6))
    col = out.members["col"]
    assert col.status == "incomplete" and col.passed is False and any(s.startswith("torsion:") for s in col.not_checked)
    assert any("not all limit states" in w for w in col.warnings)
    # a failing check dominates
    out = design.check(m, _hand_results(m, "secondOrder", mz=900.0e6, t=5.0e6))
    assert out.status == "fail"
    # a strength combination missing from the results -> "incomplete"
    two = dict(second, U2={"factors": {"D": 1.2}})
    m2 = build_model(_column(combos=two))
    out = design.check(m2, _hand_results(m2, "secondOrder"))
    assert out.status == "incomplete" and any("'U2'" in w and "not in the analysis results" in w for w in out.warnings)
    assert out.combinations["U2"]["inResults"] is False


def test_analysis_method_warnings():
    # second-order without stiffness reduction and without notional loads: not a direct analysis
    m = build_model(_column(combos={"U1": {"factors": {"D": 1.0}, "analysisType": "secondOrder"}}))
    out = design.check(m, _hand_results(m, "secondOrder"))
    assert out.status == "pass"                                             # second-order required strengths
    assert any("WITHOUT the stiffness reduction" in w and "C2.3" in w for w in out.warnings)
    assert any("no notional load" in w and "C2.2" in w for w in out.warnings)
    assert out.combinations["U1"]["directAnalysisSettings"] is False
    assert any("K = 1" in w and "direct analysis" in w for w in out.members["col"].warnings)
    # ASD with a second-order analysis at alpha = 1.0: C2.1(d)
    m = build_model(_column({"defaults": {"method": "ASD"}},
                            combos={"U1": {"factors": {"D": 1.0}, "analysisType": "secondOrder",
                                           "stiffnessReduction": True, "notionalLoad": {"direction": [1, 0, 0]}}}))
    out = design.check(m, _hand_results(m, "secondOrder"))
    assert any("ASD" in w and "1.6" in w and "C2.1(d)" in w for w in out.warnings)
    basis, warnings, conditions, incomplete = design.analysis_basis(m, _hand_results(m, "firstOrder"))
    assert len(conditions) == 1 and "first-order" in conditions[0] and not incomplete and warnings == []


def test_envelope_over_combinations_keeps_the_governing_one_per_limit_state():
    combos = {c: {"factors": {"D": 1.0}, "analysisType": "secondOrder", "stiffnessReduction": True,
                  "notionalLoad": {"direction": [1, 0, 0]}} for c in ("U1", "U2", "U3")}
    m = build_model(_column(combos=combos))
    r = _hand_results(m, "secondOrder", combos=("U1", "U2", "U3"))
    r.combinations["U2"].members["col"].forces[:, 5] *= 3.0                 # largest moment in U2
    r.combinations["U2"].members["col"].forces[:, 1] *= 3.0
    r.combinations["U3"].members["col"].forces[:, 0] *= 4.0                 # largest compression in U3
    col = design.check(m, r).members["col"]
    checks = _by_name(col.checks)
    assert checks[ff.NAME_Y["x"]].combination == "U2" and checks[e_compression.NAME_FB["y"]].combination == "U3"
    per = {cid: _by_name(cs) for cid, cs in col.by_combination.items()}
    assert set(per) == {"U1", "U2", "U3"}
    for name, c in checks.items():
        if c.status == "ok" and not c.advisory and name in per["U1"]:
            assert c.ratio == max(per[cid][name].ratio for cid in per if name in per[cid])
    h = checks[hc.NAME_H1["c"]]
    assert h.ratio == max(per[cid][hc.NAME_H1["c"]].ratio for cid in per) and col.governing.ratio >= h.ratio


# ============================================================================== second-order (activates with the solver)
@needs_second_order
def test_second_order_analysis_gives_larger_moment_and_no_first_order_condition():
    first = build_model(_column())
    second = build_model(_column(combos={"U1": {"factors": {"D": 1.0}, "analysisType": "secondOrder"}}))
    o1 = design.check(first, solve.analyze(first))
    o2 = design.check(second, solve.analyze(second))
    h1 = _by_name(o1.members["col"].checks)[hc.NAME_H1["c"]]
    h2 = _by_name(o2.members["col"].checks)[hc.NAME_H1["c"]]
    assert h2.details["Mrx"][0] > h1.details["Mrx"][0] and h2.ratio > h1.ratio      # P-delta amplifies the moment
    assert o2.combinations["U1"]["secondOrder"] is True and o2.status != "conditional"
    assert not any("FIRST-ORDER" in n for n in h2.notes)
    # the P-delta amplification of the mid-height moment is close to B1 with Cm = 1 (App. 8, Eq. A-8-3)
    sec = second.sections["I400"]
    B1 = hc.b1_a_8_3(1.0, 1.0, 400.0e3, hc.pe1_a_8_5(E * sec.Iz, 4000.0))
    assert h2.details["Mrx"][0] / h1.details["Mrx"][0] == pytest.approx(B1, rel=0.02)


# ============================================================================== CLI
def test_cli_design_writes_the_document_and_prints_the_table(tmp_path, capsys):
    path = tmp_path / "beam.json"
    path.write_text(json.dumps(_beam({"members": {"beam": {"bracing": {"topFlange": [0.5]},
                                                          "deflectionLimits": [{"spanRatio": 360}]}}})))
    out_path, res_path = tmp_path / "beam.design.json", tmp_path / "beam.results.json"
    code = main(["design", str(path), "--results", str(res_path)])
    text = capsys.readouterr().out
    assert code == 0                                            # "conditional" (first-order) is not a failure
    assert "beam" in text and "conditional" in text and "AISC 360-22" in text and "CONDITIONAL" in text
    assert "governing limit state [clause]" in text and "warning:" in text and "first-order" in text
    doc = json.loads(out_path.read_text())
    assert doc["format"] == "benchfea-design" and doc["status"] == "conditional" and doc["code"] == "AISC 360-22"
    assert doc["members"]["beam"]["governing"]["name"] in text
    assert json.loads(res_path.read_text())["format"] == "benchfea-results"
    # --json: machine-readable summary; a failing design gives exit code 4
    heavy = tmp_path / "heavy.json"
    heavy.write_text(json.dumps(_beam(w=-200.0)))
    code = main(["design", str(heavy), "-o", str(tmp_path / "h.json"), "--json"])
    summary = json.loads(capsys.readouterr().out)
    assert code == 4 and summary["status"] == "fail" and summary["ok"] is False
    row, = summary["members"]
    assert row["member"] == "beam" and row["ratio"] > 1.0 and row["clause"].startswith("AISC 360-22")
    assert summary["combinations"]["ULS"]["analysisType"] == "firstOrder"
    # a section outside the scope: "incomplete", exit code 4, the reasons are printed
    odd = _beam()
    odd["sections"]["I400"] = {"type": "rectangle", "h": 300, "b": 150}
    p3 = tmp_path / "odd.json"
    p3.write_text(json.dumps(odd))
    code = main(["design", str(p3), "-o", str(tmp_path / "o.json")])
    text = capsys.readouterr().out
    assert code == 4 and "incomplete" in text and "not checked --" in text and "outside the implemented scope" in text
    # an invalid model is refused: exit code 1, nothing written
    bad = _beam()
    del bad["supports"]
    p4 = tmp_path / "bad.json"
    p4.write_text(json.dumps(bad))
    assert main(["design", str(p4), "-o", str(tmp_path / "none.json")]) == 1
    assert not (tmp_path / "none.json").exists()
    capsys.readouterr()
    assert main(["design", str(tmp_path / "absent.json")]) == 2


# ============================================================================== free ends (F1 cantilever rule)
def test_cantilever_and_overhang_from_the_model():
    # --- cantilever: built-in at A, free tip B with a point load: Mz = -P (L - x) (bottom flange in compression)
    data = _beam({"members": {"beam": {"deflectionLimits": [{"spanRatio": 180}]}}})
    data["supports"] = {"A": {"type": "fixed"}}
    data["loadCases"] = {"D": {"loads": [{"type": "nodal", "node": "B", "Fz": -10.0e3}]}}
    m = build_model(data)
    ctx = design.DesignContext.from_model(m, "beam")
    assert ctx.end_free == (False, True) and ctx.end_built_in == (True, False)
    out = design.check(m, solve.analyze(m))
    checks = _by_name(out.members["beam"].checks)
    ltb = checks[ff.NAME_LTB]
    # F1: cantilever with the free end unbraced, warping prevented at the support: Cb = 1.0 (not Eq. F1-1,
    # which would give 1.67 for this triangular diagram); Lb = the cantilever length
    cap, _, _, _ = _f2_strength(m.sections["I400"], 6000.0, 1.0)
    assert ltb.status == "ok" and ltb.inputs["Cb"][0] == 1.0 and ltb.capacity == pytest.approx(cap, rel=1e-9)
    assert ltb.demand == pytest.approx(1.5 * 10.0e3 * 6000.0, rel=1e-9) and ltb.x == pytest.approx(0.0)
    assert any("cantilever" in n and "ASSUMED" in n for n in ltb.notes)
    # the deflection limit of a cantilever is not checked (chord through both ends) -> never "pass"
    d = checks["deflection limit 1 (resultant)"]
    assert d.status == "notChecked" and "free end" in d.reason and out.members["beam"].status == "incomplete"
    # --- overhang: A --- B --- C with supports at A and B, load at the free end C
    data = _beam()
    data["nodes"] = {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": 4000, "y": 0, "z": 0}, "C": {"x": 6000, "y": 0, "z": 0}}
    data["members"] = {"span": {"i": "A", "j": "B", "section": "I400", "material": "S"},
                       "tip": {"i": "B", "j": "C", "section": "I400", "material": "S"}}
    data["loadCases"] = {"D": {"loads": [{"type": "nodal", "node": "C", "Fz": -10.0e3}]}}
    m = build_model(data)
    assert design.DesignContext.from_model(m, "tip").end_free == (False, True)
    assert design.DesignContext.from_model(m, "tip").end_built_in == (False, False)
    assert design.DesignContext.from_model(m, "span").end_free == (False, False)
    out = design.check(m, solve.analyze(m))
    tip = _by_name(out.members["tip"].checks)[ff.NAME_LTB]
    # the overhang is not a cantilever with warping prevented at the support: F1 gives no Cb -> notChecked
    assert tip.status == "notChecked" and "free end" in tip.reason and tip.demand == pytest.approx(1.5 * 10.0e3 * 2000.0, rel=1e-9)
    assert out.members["tip"].status == "incomplete" and out.members["tip"].passed is False
    assert _by_name(out.members["span"].checks)[ff.NAME_LTB].status == "ok"
    # ... unless the engineer gives Cb (and Lb)
    data["design"] = {"members": {"tip": {"Cb": 1.0, "Lb": 2000}}}
    m = build_model(data)
    tip = _by_name(design.check(m, solve.analyze(m)).members["tip"].checks)[ff.NAME_LTB]
    assert tip.status == "ok" and tip.inputs["Cb"] == (1.0, "dimensionless", "member")
