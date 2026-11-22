"""AISC 360-22 Chapter D (tension): D1 slenderness (advisory), D2 yielding and
rupture with Ae = An U -- and the path analysis results -> demands -> checks ->
design output document.

Tier F: (a) hand calculations from the Specification (pp. 16.1-32, 16.1-33) as
explicit arithmetic; (b) AISC Design Examples v16.0 D.1, D.4, D.5.
"""
import json
import math

import numpy as np
import pytest

from benchfea import design
from benchfea.design import DesignContext, DesignSection, Demand
from benchfea.design.aisc360_22 import d_tension
from benchfea.model import build_model
from benchfea.results import CombinationResult, MemberResult, Results
from design_examples import FT, assert_printed, example


def _ctx(method, sec=None, **kw):
    sec = sec or DesignSection.i_shape("rolled", A=5000.0, rx=80.0, ry=25.0, flange_ratio=8.0, web_ratio=30.0)
    base = dict(member="T", section=sec, method=method, E=200000.0, G=77000.0, Fy=345.0, Fu=450.0, L=6000.0,
                Lx=6000.0, Ly=6000.0, Lz=6000.0)
    base.update(kw)
    return DesignContext(**base)


def _by_name(checks):
    return {c.name: c for c in checks}


# ============================================================================== tier F (a): hand calculations
def test_hand_tensile_yielding_d2a():
    # Eq. D2-1: Pn = Fy Ag = 345 * 5000 = 1 725 000 N
    Pn = 345.0 * 5000.0
    assert Pn == 1725000.0
    demand = Demand("U1", Pt=900000.0, x_Pt=1500.0)
    c = d_tension.yielding(_ctx("LRFD"), demand)[0]
    # phi_t = 0.90: phi Pn = 1 552 500 N; ratio = 900000 / 1552500 = 0.57971
    assert c.nominal == pytest.approx(Pn, rel=1e-15) and (c.factor, c.factor_name) == (0.90, "phi")
    assert c.capacity == pytest.approx(0.90 * 1725000.0, rel=1e-15)
    assert c.ratio == pytest.approx(900000.0 / 1552500.0, rel=1e-14) and c.passed is True
    assert c.clause == "AISC 360-22 D2(a), Eq. D2-1; B3.1 Eq. B3-1"
    assert (c.combination, c.x, c.dimension) == ("U1", 1500.0, "force")
    assert c.inputs["Fy"] == (345.0, "force/length^2", "given") and c.inputs["Ag"][0] == 5000.0
    # Omega_t = 1.67: Pn / Omega = 1 032 934.13 N
    c = d_tension.yielding(_ctx("ASD"), demand)[0]
    assert c.capacity == pytest.approx(1725000.0 / 1.67, rel=1e-15) and (c.factor, c.factor_name) == (1.67, "Omega")
    assert 1032934.0 < c.capacity < 1032935.0
    assert c.clause.endswith("B3.2 Eq. B3-2")


def test_hand_tensile_rupture_d2b_with_an_and_u():
    # Eq. D3-1: Ae = An U = 4200 * 0.85 = 3570 mm^2; Eq. D2-2: Pn = Fu Ae = 450 * 3570 = 1 606 500 N
    Ae = 4200.0 * 0.85
    Pn = 450.0 * Ae
    assert (Ae, Pn) == (3570.0, 1606500.0)
    kw = dict(An=4200.0, U=0.85, sources={"An": "member", "U": "group:g1", "Fu": "material:S"})
    c = d_tension.rupture(_ctx("LRFD", **kw), Demand("U1", Pt=1.3e6))[0]
    # phi_t = 0.75: 1 204 875 N; ratio = 1.3e6 / 1 204 875 = 1.07895 -> fails
    assert c.nominal == pytest.approx(Pn, rel=1e-15) and c.capacity == pytest.approx(0.75 * 1606500.0, rel=1e-15)
    assert c.ratio == pytest.approx(1.3e6 / 1204875.0, rel=1e-14) and c.passed is False and c.factor == 0.75
    assert c.clause == "AISC 360-22 D2(b), Eq. D2-2; D3, Eq. D3-1; B3.1 Eq. B3-1"
    assert c.details["Ae = An U (Eq. D3-1)"] == (pytest.approx(3570.0), "length^2")
    assert c.inputs["An"] == (4200.0, "length^2", "member") and c.inputs["U"] == (0.85, "dimensionless", "group:g1")
    assert c.inputs["Fu"][2] == "material:S" and not c.notes
    # Omega_t = 2.00: 803 250 N
    c = d_tension.rupture(_ctx("ASD", **kw), Demand("U1", Pt=8.0e5))[0]
    assert c.capacity == pytest.approx(1606500.0 / 2.00, rel=1e-15) and c.factor == 2.00 and c.passed is True


def test_hand_governing_tension_limit_state_differs_between_lrfd_and_asd():
    """D2: the lower AVAILABLE strength governs.  With An U / Ag = 0.92:
    LRFD yield 0.9 Fy Ag = 1 552 500 vs rupture 0.75 * 450 * 4600 = 1 552 500 (tie);
    ASD  yield 1725000/1.67 = 1 032 934 vs rupture 2 070 000/2.00 = 1 035 000 -> yielding governs."""
    kw = dict(An=4600.0, U=1.0)
    lr = _by_name(d_tension.check(_ctx("LRFD", **kw), Demand("U", Pt=1.0e6)))
    asd = _by_name(d_tension.check(_ctx("ASD", **kw), Demand("U", Pt=1.0e6)))
    y, r = d_tension.NAME_YIELD, d_tension.NAME_RUPTURE
    assert lr[y].capacity == pytest.approx(lr[r].capacity, rel=1e-14)
    assert asd[y].capacity < asd[r].capacity and asd[y].ratio > asd[r].ratio
    res = design.check_member(_ctx("ASD", **kw), [Demand("U", Pt=1.0e6)])
    assert res.governing.name == y


def test_hand_slenderness_d1_is_advisory():
    # D1 User Note: L / r_min = 6000 / 25 = 240 <= 300; reported ratio 240/300 = 0.8
    c = d_tension.slenderness(_ctx("LRFD"), Demand())[0]
    assert (c.demand, c.capacity, c.ratio, c.passed, c.advisory) == (240.0, 300.0, 0.8, True, True)
    assert c.clause == "AISC 360-22 D1 (User Note)" and c.dimension == "dimensionless" and c.nominal is None
    assert c.inputs["r_min"][0] == 25.0 and c.inputs["L"][0] == 6000.0
    # L = 9000: 360 > 300 -> advisory check fails, the member does not
    ctx = _ctx("LRFD", L=9000.0)
    c = d_tension.slenderness(ctx, Demand())[0]
    assert c.demand == 360.0 and c.passed is False
    res = design.check_member(ctx.with_(section=DesignSection.i_shape(
        "rolled", A=5000.0, rx=80.0, ry=25.0, J=1.0e5, Cw=1.0e10, flange_ratio=8.0, web_ratio=30.0)),
        [Demand("U", Pt=1.0e5)])
    assert res.passed is True


def test_defaults_an_equals_ag_and_u_one_are_reported():
    ctx = _ctx("LRFD", sources={"U": "default"})          # An None -> gross area; U default 1.0
    c = d_tension.rupture(ctx, Demand("U", Pt=1.0))[0]
    assert c.nominal == pytest.approx(450.0 * 5000.0, rel=1e-15)
    assert c.inputs["An"][0] == 5000.0 and c.inputs["An"][2].startswith("default")
    assert c.inputs["U"] == (1.0, "dimensionless", "default")
    assert any("defaults" in n and "end connection" in n for n in c.notes)


def test_refusals_and_not_checked():
    d = Demand("U", Pt=1.0e5)
    c = d_tension.rupture(_ctx("LRFD", An=5200.0), d)[0]              # An > Ag
    assert (c.status, c.capacity, c.ratio, c.passed) == ("refused", None, None, None) and "5200" in c.reason
    assert c.demand == 1.0e5
    for U in (0.0, -0.2, 1.01):
        assert d_tension.rupture(_ctx("LRFD", U=U), d)[0].status == "refused"
    c = d_tension.rupture(_ctx("LRFD", Fu=None), d)[0]
    assert c.status == "notChecked" and "Fu" in c.reason
    c = d_tension.yielding(_ctx("LRFD", Fy=None), d)[0]
    assert c.status == "notChecked" and "Fy" in c.reason
    # Chapter D needs only the area: a shape outside the scope is still checked in tension
    other = DesignSection.unsupported("angle section", "L1", A=2000.0, Ix=3.0e6, Iy=1.0e6)
    checks = _by_name(d_tension.check(_ctx("LRFD", sec=other), d))
    assert checks[d_tension.NAME_YIELD].capacity == pytest.approx(0.9 * 345.0 * 2000.0, rel=1e-15)
    assert checks[d_tension.NAME_SLENDER].demand == pytest.approx(6000.0 / math.sqrt(1.0e6 / 2000.0), rel=1e-14)
    nothing = DesignSection.unsupported("no data", "X")
    assert all(c.status == "notChecked" for c in d_tension.check(_ctx("LRFD", sec=nothing), d))


# ============================================================================== tier F (b): Design Examples
def _de_tension(eid, sec):
    x = example(eid)
    i, p = x["inputs"], x["printed"]
    out = {}
    for method, demand in (("LRFD", i["Pu"]), ("ASD", i["Pa"])):
        ctx = DesignContext(member=eid, section=sec, method=method, E=29000.0, G=11200.0, Fy=i["Fy"], Fu=i["Fu"],
                            L=i["L_ft"] * FT, Lx=None, Ly=None, Lz=None, An=i["An"], U=i["U"])
        out[method] = _by_name(d_tension.check(ctx, Demand(method, Pt=demand)))
    y, r = d_tension.NAME_YIELD, d_tension.NAME_RUPTURE
    assert_printed(out["LRFD"][y].capacity, p["yield_phiPn"], f"{eid} yielding phi Pn")
    assert_printed(out["ASD"][y].capacity, p["yield_Pn_over_Omega"], f"{eid} yielding Pn/Omega")
    assert_printed(out["LRFD"][r].details["Ae = An U (Eq. D3-1)"][0], p["Ae"], f"{eid} Ae")
    assert_printed(out["LRFD"][r].nominal, p["rupture_Pn"], f"{eid} rupture Pn")
    assert_printed(out["LRFD"][r].capacity, p["rupture_phiPn"], f"{eid} rupture phi Pn")
    assert_printed(out["ASD"][r].capacity, p["rupture_Pn_over_Omega"], f"{eid} rupture Pn/Omega")
    # the examples state that rupture governs and that the member is adequate
    for method in ("LRFD", "ASD"):
        assert out[method][r].ratio > out[method][y].ratio and out[method][r].passed
    return out, p


def test_de_d1_w_shape_tension_member():
    """Example D.1, W8x21 (D-2..D-4, PDF 38-40): yielding 277/184, rupture Pn 281 -> 211/141 kips.
    Listed: Ag, bf, tf, d, ry.  An = 4.76 in.^2 and U = 0.908 are the example's connection values (inputs)."""
    i = example("D.1")["inputs"]
    sec = DesignSection.i_shape("rolled", A=i["Ag"], ry=i["ry"], d=i["d"], bf=i["bf"], tf=i["tf"])
    out, p = _de_tension("D.1", sec)
    # the example lists only ry (the least radius): without rx the engine does not assume which is least
    assert out["LRFD"][d_tension.NAME_SLENDER].status == "notChecked"
    assert_printed(i["L_ft"] * FT / i["ry"], p["L_over_r"], "D.1 L/r (plain arithmetic)")


def test_de_d4_rectangular_hss_tension_member():
    """Example D.4, HSS6x4x3/8 (D-11..D-13, PDF 47-49): yielding 278/185, rupture Pn 330 -> 248/165 kips."""
    i = example("D.4")["inputs"]
    sec = DesignSection.rectangular("rectHSS", A=i["Ag"], ry=i["ry"], H=6.0, B=4.0, t=i["t"])
    _de_tension("D.4", sec)


def test_de_d5_round_hss_tension_member_and_slenderness():
    """Example D.5, HSS6.000x0.500 (D-14..D-16, PDF 50-52): yielding 364/242, rupture Pn 465 -> 349/233 kips;
    L/r = 184 < 300 (Section D1)."""
    i = example("D.5")["inputs"]
    sec = DesignSection.round_hss(A=i["Ag"], r=i["r"], D=i["D"], t=i["t"])
    out, p = _de_tension("D.5", sec)
    s = out["LRFD"][d_tension.NAME_SLENDER]
    assert_printed(s.demand, p["L_over_r"], "D.5 L/r")
    assert s.capacity == 300.0 and s.passed and s.advisory


# ============================================================================== results -> demands -> checks
def _chain_model(units=None, design=None):
    kip_in = units == {"force": "kip", "length": "in"}
    L = 120.0 if kip_in else 3000.0
    data = {
        "formatVersion": "1.0.0", "units": units or {"force": "N", "length": "mm"},
        "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": L, "y": 0, "z": 0}, "C": {"x": 2 * L, "y": 0, "z": 0}},
        "supports": {"A": {"type": "pinned", "dofs": {"rx": "fixed"}}, "C": {"type": "roller", "dofs": {"uy": "fixed"}}},
        "materials": {"S": ({"E": 29000, "nu": 0.3, "Fy": 50, "Fu": 65} if kip_in
                            else {"E": 200000, "nu": 0.3, "Fy": 345, "Fu": 450})},
        "sections": {"chs": ({"type": "pipe", "D": 8.0, "t": 0.25} if kip_in
                             else {"type": "pipe", "D": 200, "t": 6})},
        "members": {"m1": {"i": "A", "j": "B", "section": "chs", "material": "S"},
                    "m2": {"i": "B", "j": "C", "section": "chs", "material": "S"}},
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fz": -1}]}},
        "combinations": {"U1": {"factors": {"D": 1.4}}, "U2": {"factors": {"D": 1.2}},
                         "S1": {"factors": {"D": 1.0}, "purpose": "serviceability"}},
        "design": {"members": {"tie": {"members": ["m1", "m2"], **(design or {})}}},
    }
    return build_model(data)


def _results(model, fields):
    """Results with station forces given by ``fields[combination](X) -> (N, Vy, Vz)`` in terms of the
    distance X from node A (members m1, m2 run along +X), 11 stations per member."""
    combos = {}
    for cid, fn in fields.items():
        members = {}
        for mid, m in model.members.items():
            x = np.linspace(0.0, m.length, 11)
            f = np.zeros((11, 6))
            for k, xk in enumerate(x):
                f[k, :3] = fn(m.start[0] + xk)
            members[mid] = MemberResult(np.concatenate([-f[0], f[-1]]), x, f, np.zeros((11, 3)))
        n = len(model.nodes)
        combos[cid] = CombinationResult(cid, "firstOrder", np.zeros((n, 6)), np.zeros((n, 6)), members)
    return Results(model.name, combos)


def test_demands_of_a_chain_of_members_per_combination():
    m = _chain_model()
    L = 6000.0
    res = _results(m, {
        # tension rising linearly along the chain, shear of a simply supported beam under uniform load
        "U1": lambda X: (1.0e5 + 20.0 * X, 3.0e4 * (1.0 - 2.0 * X / L), 0.0),
        # compression with its maximum at midspan (the junction of the two members), constant shear in z
        "U2": lambda X: (-5.0e4 - 30.0 * min(X, L - X), 0.0, -8.0e3),
        "S1": lambda X: (9.9e9, 9.9e9, 9.9e9),        # serviceability: never a strength demand
    })
    d = design.member_demands(m, res, "tie")
    assert list(d) == ["U1", "U2"]
    u1, u2 = d["U1"], d["U2"]
    assert (u1.Pt, u1.x_Pt, u1.Pc, u1.x_Pc) == (pytest.approx(2.2e5), pytest.approx(L), 0.0, None)
    assert (u1.Vy, u1.x_Vy) == (pytest.approx(3.0e4), 0.0)
    assert u1.Lv_y == pytest.approx(L / 2.0) and u1.Lv == pytest.approx(L / 2.0)      # zero shear at midspan
    assert (u1.Vz, u1.x_Vz, u1.Lv_z) == (0.0, None, None)
    assert (u2.Pc, u2.x_Pc, u2.Pt) == (pytest.approx(5.0e4 + 30.0 * 3000.0), pytest.approx(3000.0), 0.0)
    assert (u2.Vz, u2.Lv_z, u2.Lv) == (pytest.approx(8.0e3), None, None)            # shear never vanishes
    assert u2.V == pytest.approx(8.0e3)
    assert list(design.member_demands(m, res, "tie", purpose=None)) == ["U1", "U2", "S1"]


def test_demand_resultant_shear_and_distance_to_zero():
    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    f = np.zeros((5, 6))
    f[:, 1] = [30.0, 15.0, 0.0, -15.0, -30.0]          # Vy
    f[:, 2] = [40.0, 20.0, 0.0, -20.0, -40.0]          # Vz
    d = design.demand_from_stations(x, f, "c")
    assert d.V == pytest.approx(50.0) and d.x_V == 0.0 and d.Lv == pytest.approx(2.0)
    # a jump through zero between two stations at the same position (concentrated load)
    x = np.array([0.0, 1.0, 1.0, 3.0])
    f = np.zeros((4, 6))
    f[:, 1] = [20.0, 20.0, -10.0, -10.0]
    d = design.demand_from_stations(x, f)
    assert (d.Vy, d.x_Vy, d.Lv_y) == (20.0, 0.0, pytest.approx(1.0))
    # sign change inside an interval: linear interpolation, V = 10 at x = 0 and -30 at x = 4 -> zero at x = 1
    f = np.zeros((2, 6))
    f[:, 1] = [10.0, -30.0]
    d = design.demand_from_stations(np.array([0.0, 4.0]), f)
    assert (d.Vy, d.x_Vy, d.Lv_y) == (30.0, 4.0, pytest.approx(3.0))


def test_check_end_to_end_and_design_document_in_file_units():
    """kip-in model: demands from (hand-made) results, checks in internal units, document in kip, in."""
    m = _chain_model({"force": "kip", "length": "in"}, {"An": 5.0, "U": 0.9, "method": "ASD"})
    kip = 4448.2216152605
    res = _results(m, {"U1": lambda X: (100.0 * kip, 0.0, 0.0), "U2": lambda X: (60.0 * kip, 0.0, 0.0),
                       "S1": lambda X: (500.0 * kip, 0.0, 0.0)})
    out = design.check(m, res)
    assert out.code == "AISC 360-22" and list(out.members) == ["tie"]
    tie = out.members["tie"]
    A = math.pi * (8.0 ** 2 - 7.5 ** 2) / 4.0                       # pipe 8 x 0.25: 6.0868 in.^2
    checks = _by_name(tie.checks)
    y, r = checks[d_tension.NAME_YIELD], checks[d_tension.NAME_RUPTURE]
    # ASD: yielding 50 A / 1.67; rupture 65 * (5.0 * 0.9) / 2.00 = 146.25 kips; demand 100 kips (U1)
    assert y.capacity / kip == pytest.approx(50.0 * A / 1.67, rel=1e-12)
    assert r.capacity / kip == pytest.approx(65.0 * 5.0 * 0.9 / 2.0, rel=1e-12)
    assert (r.combination, r.method) == ("U1", "ASD") and r.ratio == pytest.approx(100.0 / 146.25, rel=1e-12)
    assert tie.governing is r and tie.passed
    doc = out.to_document(m)
    json.dumps(doc)                                                 # plain JSON data
    d = {c["name"]: c for c in doc["members"]["tie"]["checks"]}[d_tension.NAME_RUPTURE]
    assert d["capacity"] == pytest.approx(146.25, rel=1e-12) and d["unit"] == "kip"
    assert d["demand"] == pytest.approx(100.0, rel=1e-12) and d["nominal"] == pytest.approx(292.5, rel=1e-12)
    assert (d["status"], d["method"], d["factor"], d["factorName"]) == ("ok", "ASD", 2.0, "Omega")
    assert d["inputs"]["An"] == {"value": pytest.approx(5.0), "unit": "in^2", "source": "member"}
    assert d["inputs"]["Fu"] == {"value": pytest.approx(65.0), "unit": "kip/in^2", "source": "material:S"}
    assert d["details"]["Ae = An U (Eq. D3-1)"]["value"] == pytest.approx(4.5)
    assert d["clause"] == "AISC 360-22 D2(b), Eq. D2-2; D3, Eq. D3-1; B3.2 Eq. B3-2"
    s = {c["name"]: c for c in doc["members"]["tie"]["checks"]}[d_tension.NAME_SLENDER]
    assert s["advisory"] and s["inputs"]["L"]["value"] == pytest.approx(240.0) and s["inputs"]["L"]["unit"] == "in"
    assert doc["members"]["tie"]["governing"]["name"] == d_tension.NAME_RUPTURE


def test_check_reports_not_checked_limit_states_in_the_document():
    data_units = {"force": "N", "length": "mm"}
    m = _chain_model(data_units)
    # replace the design view: a material without Fu makes rupture impossible to check
    m.design_members["tie"].values["Fu"] = type(m.design_members["tie"].values["Fu"])(None, "material:S", "/x")
    res = _results(m, {"U1": lambda X: (1.0e4, 0.0, 0.0), "U2": lambda X: (0.0, 0.0, 0.0)})
    out = design.check(m, res)
    tie = out.members["tie"]
    assert tie.passed is False and any("rupture" in s for s in tie.not_checked)
    doc = out.to_document(m)
    json.dumps(doc)
    c = {c["name"]: c for c in doc["members"]["tie"]["checks"]}[d_tension.NAME_RUPTURE]
    assert (c["status"], c["capacity"], c["ratio"], c["passed"]) == ("notChecked", None, None, None)
    assert c["demand"] == pytest.approx(1.0e4) and "Fu" in c["reason"]
    assert doc["members"]["tie"]["notChecked"] and doc["members"]["tie"]["passed"] is False
    # a non-converged strength combination is flagged
    res.combinations["U1"].converged = False
    assert any("U1" in w and "converge" in w for w in design.check(m, res).warnings)
