"""AISC 360-22 Chapter G (shear): G2.1 webs of I-shapes (incl. the phi = 1.00 /
Omega = 1.50 case and the Cv1 branches), G4 rectangular HSS / box, G5 round HSS,
G6 minor-axis shear of I-shapes.  Tension field action is not used.

Tier F: (a) hand calculations from the Specification (pp. 16.1-75 to 16.1-81) as
explicit arithmetic; (b) AISC Design Examples v16.0 G.1B, G.4, G.5, G.6.
"""
import math

import numpy as np
import pytest

from benchfea import design
from benchfea.design import DesignContext, DesignSection, Demand
from benchfea.design.aisc360_22 import g_shear as gs
from benchfea.model import build_model
from benchfea.results import CombinationResult, MemberResult, Results
from design_examples import assert_printed, example

E, FY = 200000.0, 345.0                 # MPa
ROOT = math.sqrt(E / FY)                # 24.0772
VY, VZ, V = gs.NAME_VY, gs.NAME_VZ, gs.NAME_V


def _ctx(sec, method="LRFD", **kw):
    base = dict(member="B", section=sec, method=method, E=E, G=77200.0, Fy=FY, Fu=450.0, L=6000.0, Lx=6000.0,
                Ly=6000.0, Lz=6000.0)
    base.update(kw)
    return DesignContext(**base)


def _checks(sec, demand, method="LRFD", **kw):
    return {c.name: c for c in gs.check(_ctx(sec, method, **kw), demand)}


def _i(fabrication, d, bf, tf, tw, **kw):
    return DesignSection.i_shape(fabrication, A=2 * bf * tf + (d - 2 * tf) * tw, rx=d / 2.4, ry=bf / 4.5, d=d, bf=bf,
                                 tf=tf, tw=tw, **kw)


# ============================================================================== tier F (a): hand calculations
def test_hand_g21a_rolled_web_phi_1_00_omega_1_50():
    # rolled I 400 x 200, tf 12, tw 8, r 18: h = 400 - 24 - 36 = 340, h/tw = 42.5 <= 2.24 sqrt(E/Fy) = 53.93
    assert 53.9 < 2.24 * ROOT < 54.0
    # G2.1(a): Cv1 = 1.0 (Eq. G2-2); Aw = d tw = 400 * 8 = 3200; Eq. G2-1: Vn = 0.6 * 345 * 3200 * 1.0 = 662 400 N
    Vn = 0.6 * 345.0 * 3200.0 * 1.0
    assert Vn == pytest.approx(662400.0, rel=1e-15)
    sec = _i("rolled", 400.0, 200.0, 12.0, 8.0, r=18.0)
    c = _checks(sec, Demand("U1", Vy=5.0e5, x_Vy=0.0))[VY]
    assert c.nominal == pytest.approx(Vn, rel=1e-14) and (c.factor, c.factor_name) == (1.00, "phi")
    assert c.capacity == pytest.approx(1.00 * Vn, rel=1e-14) and c.ratio == pytest.approx(5.0e5 / 662400.0, rel=1e-14)
    assert c.clause == "AISC 360-22 G2.1(a), Eq. G2-1, G2-2; B3.1 Eq. B3-1" and c.passed
    assert c.inputs["section fabrication"][0] == "rolled" and c.inputs["h/tw"][0] == pytest.approx(42.5)
    a = _checks(sec, Demand("U1", Vy=4.5e5), "ASD")[VY]
    # Omega_v = 1.50: 662 400 / 1.50 = 441 600 N; 450 000 > 441 600 -> fails
    assert a.capacity == pytest.approx(441600.0, rel=1e-14) and (a.factor, a.factor_name) == (1.50, "Omega")
    assert a.passed is False


def test_hand_g21b_built_up_web_cv1_one_phi_0_90():
    # the same dimensions welded (no fillet): h = 400 - 24 = 376, h/tw = 47; G2.1(a) is for ROLLED shapes only
    # G2.1(b): kv = 5.34 (no transverse stiffeners); 1.10 sqrt(kv E/Fy) = 1.10 * 55.638 = 61.20 >= 47
    limit = 1.10 * math.sqrt(5.34 * E / FY)
    assert 61.1 < limit < 61.3
    # Cv1 = 1.0 (Eq. G2-3); Vn = 662 400; phi_v = 0.90 -> 596 160; Omega_v = 1.67 -> 396 646.7
    sec = _i("builtUp", 400.0, 200.0, 12.0, 8.0)
    c = _checks(sec, Demand("U1", Vy=1.0e5))[VY]
    assert c.nominal == pytest.approx(662400.0, rel=1e-14) and c.capacity == pytest.approx(596160.0, rel=1e-14)
    assert c.clause == "AISC 360-22 G2.1(b), Eq. G2-1, G2-3; G1; B3.1 Eq. B3-1" and c.factor == 0.90
    assert c.details["kv (web without transverse stiffeners)"][0] == 5.34
    a = _checks(sec, Demand("U1", Vy=1.0e5), "ASD")[VY]
    assert a.capacity == pytest.approx(662400.0 / 1.67, rel=1e-14) and a.factor == 1.67


def test_hand_g21b_cv1_below_one_built_up_and_rolled():
    limit = 1.10 * math.sqrt(5.34 * E / FY)       # 61.20
    # welded I 600 x 300, tf 10, tw 6: h = 580, h/tw = 96.667 > 61.20 -> Eq. G2-4: Cv1 = 61.20 / 96.667 = 0.6331
    cv1 = limit / (580.0 / 6.0)
    assert 0.633 < cv1 < 0.634
    Vn = 0.6 * 345.0 * (600.0 * 6.0) * cv1        # Eq. G2-1: 745 200 * 0.6331 = 471 800 N
    assert 471000.0 < Vn < 472500.0
    c = _checks(_i("builtUp", 600.0, 300.0, 10.0, 6.0), Demand("U1", Vy=1.0e5))[VY]
    assert c.nominal == pytest.approx(Vn, rel=1e-14) and c.capacity == pytest.approx(0.9 * Vn, rel=1e-14)
    assert "Eq. G2-1, G2-4" in c.clause and gs.cv1(580.0 / 6.0, 5.34, E, FY) == (pytest.approx(cv1), "G2-4")
    # rolled I 500 deep, tf 10, tw 6, r 10: h = 460, h/tw = 76.67 > 2.24 sqrt(E/Fy) = 53.93: NOT G2.1(a);
    # G2.1(b): 76.67 > 61.20 -> Cv1 = 61.20 / 76.667 = 0.7983, phi_v = 0.90
    cv1 = limit / (460.0 / 6.0)
    Vn = 0.6 * 345.0 * (500.0 * 6.0) * cv1
    c = _checks(_i("rolled", 500.0, 200.0, 10.0, 6.0, r=10.0), Demand("U1", Vy=1.0e5))[VY]
    assert c.nominal == pytest.approx(Vn, rel=1e-14) and c.factor == 0.90 and "G2.1(b)" in c.clause
    assert any("outside G2.1(a)" in n for n in c.notes)
    # rolled, between the two limits (53.93 < h/tw = 58 <= 61.20): G2.1(b) with Cv1 = 1.0 but phi_v = 0.90
    sec = DesignSection.i_shape("rolled", A=9000.0, rx=160.0, ry=40.0, d=400.0, bf=200.0, tf=12.0, tw=8.0, web_ratio=58.0)
    c = _checks(sec, Demand("U1", Vy=1.0e5))[VY]
    assert c.nominal == pytest.approx(662400.0, rel=1e-14) and c.factor == 0.90 and "G2-3" in c.clause
    # exactly at the limit of G2.1(a): "<=" -> (a) applies
    sec = DesignSection.i_shape("rolled", A=9000.0, rx=160.0, ry=40.0, d=400.0, bf=200.0, tf=12.0, tw=8.0,
                                web_ratio=2.24 * ROOT)
    assert _checks(sec, Demand("U1", Vy=1.0))[VY].factor == 1.00


def test_hand_cv2_three_branches_g2_9_to_g2_11():
    # kv = 5: 1.10 sqrt(kv E/Fy) = 1.10 * 53.838 = 59.22; 1.37 sqrt(kv E/Fy) = 73.76
    root = math.sqrt(5.0 * E / FY)
    assert 59.2 < 1.10 * root < 59.3 and 73.7 < 1.37 * root < 73.8
    assert gs.cv2(59.0, 5.0, E, FY) == (1.0, "G2-9")
    c, eq = gs.cv2(70.0, 5.0, E, FY)              # Eq. G2-10: 59.22 / 70 = 0.8460
    assert eq == "G2-10" and c == pytest.approx(1.10 * root / 70.0, rel=1e-15) and 0.845 < c < 0.847
    c, eq = gs.cv2(100.0, 5.0, E, FY)             # Eq. G2-11: 1.51 * 5 * 200000 / (100^2 * 345) = 0.4377
    assert eq == "G2-11" and c == pytest.approx(1.51 * 5.0 * E / (100.0 ** 2 * FY), rel=1e-15) and 0.437 < c < 0.438


def test_hand_g4_rectangular_hss_both_directions_and_box():
    # HSS 200 x 100 x 5, corner radius not known: h = 200 - 3*5 = 185 (shear along y), b = 100 - 15 = 85 (along z)
    hss = DesignSection.rectangular("rectHSS", A=2800.0, rx=72.0, ry=41.0, H=200.0, B=100.0, t=5.0)
    # h/t = 37 <= 59.22 -> Cv2 = 1.0 (Eq. G2-9); Aw = 2 h t = 2 * 185 * 5 = 1850
    # Eq. G4-1: Vn = 0.6 * 345 * 1850 * 1.0 = 382 950 N; phi_v = 0.90 -> 344 655 N
    c = _checks(hss, Demand("U1", Vy=2.0e5, Vz=1.0e5))
    assert c[VY].nominal == pytest.approx(382950.0, rel=1e-14) and c[VY].capacity == pytest.approx(344655.0, rel=1e-14)
    assert c[VY].clause == "AISC 360-22 G4, Eq. G4-1; G2.2, Eq. G2-9; G1; B3.1 Eq. B3-1"
    assert c[VY].details["kv"][0] == 5.0 and c[VY].details["Aw = 2 h t"][0] == pytest.approx(1850.0)
    # shear along z: the two 100-wide walls: Aw = 2 * 85 * 5 = 850; Vn = 0.6 * 345 * 850 = 175 950 N
    assert c[VZ].nominal == pytest.approx(175950.0, rel=1e-14) and c[VZ].demand == 1.0e5
    a = _checks(hss, Demand("U1", Vy=2.0e5), "ASD")[VY]
    assert a.capacity == pytest.approx(382950.0 / 1.67, rel=1e-14)
    # box section: h = clear distance between flanges = 200 - 2*5 = 190: Aw = 1900, Vn = 393 300 N
    box = DesignSection.rectangular("box", A=2900.0, rx=72.0, ry=41.0, H=200.0, B=100.0, t=5.0)
    assert _checks(box, Demand("U1", Vy=1.0))[VY].nominal == pytest.approx(0.6 * 345.0 * 1900.0, rel=1e-14)
    # slender walls: h/t = 100 (h = 200, t = 2): Cv2 = 0.4377 (Eq. G2-11); Vn = 0.6 * 345 * (2*200*2) * 0.4377
    thin = DesignSection.rectangular("rectHSS", A=1200.0, rx=80.0, ry=40.0, h=200.0, b=80.0, t=2.0)
    cv2 = 1.51 * 5.0 * E / (100.0 ** 2 * FY)
    c = _checks(thin, Demand("U1", Vy=1.0))[VY]
    assert c.nominal == pytest.approx(0.6 * 345.0 * 800.0 * cv2, rel=1e-14) and "Eq. G2-11" in c.clause


def test_hand_g5_round_hss_yielding_and_buckling():
    # D = 300, t = 5: D/t = 60; Ag = pi (300^2 - 290^2) / 4 = 4633.85 mm^2
    Ag = math.pi * (300.0 ** 2 - 290.0 ** 2) / 4.0
    # Eq. G5-2b: 0.78 E / (D/t)^1.5 = 156 000 / 464.76 = 335.66 > 0.6 Fy = 207 -> Fcr = 207 whatever Lv is
    fb = 0.78 * E / 60.0 ** 1.5
    assert 335.0 < fb < 336.0 and 0.6 * FY == pytest.approx(207.0)
    Vn = 207.0 * Ag / 2.0                         # Eq. G5-1: 479 603 N
    sec = DesignSection.round_hss(A=Ag, r=104.0, D=300.0, t=5.0)
    c = _checks(sec, Demand("U1", V=3.0e5, x_V=0.0, Lv=None))[V]
    assert c.nominal == pytest.approx(Vn, rel=1e-14) and c.capacity == pytest.approx(0.90 * Vn, rel=1e-14)
    assert c.clause == "AISC 360-22 G5, Eq. G5-1, G5-2a, G5-2b; G1; B3.1 Eq. B3-1"
    assert any("whatever Lv" in n for n in c.notes)
    assert _checks(sec, Demand("U1", V=3.0e5), "ASD")[V].capacity == pytest.approx(Vn / 1.67, rel=1e-14)
    # thin tube D = 500, t = 2: D/t = 250. Eq. G5-2b: 156 000 / 250^1.5 = 39.47 MPa (< 207)
    fb = 0.78 * E / 250.0 ** 1.5
    assert 39.4 < fb < 39.5
    # Lv = 3000 -> Lv/D = 6: Eq. G5-2a: 1.60 E / (sqrt(6) * 250^1.25) = 320 000 / (2.4495 * 994.1) = 131.4 MPa
    fa = 1.60 * E / (math.sqrt(3000.0 / 500.0) * 250.0 ** 1.25)
    assert 131.0 < fa < 132.0                     # larger of the two, below 0.6 Fy -> Fcr = 131.4
    A2 = math.pi * (500.0 ** 2 - 496.0 ** 2) / 4.0
    thin = DesignSection.round_hss(A=A2, r=176.0, D=500.0, t=2.0)
    c = _checks(thin, Demand("U1", V=1.0e5, Lv=3000.0))[V]
    assert c.nominal == pytest.approx(fa * A2 / 2.0, rel=1e-14) and c.details["Fcr (Eq. G5-2a)"][0] == pytest.approx(fa)
    # shear that never vanishes along the member: Lv undefined -> Eq. G5-2a is not credited (lower bound)
    c = _checks(thin, Demand("U1", V=1.0e5, Lv=None))[V]
    assert c.status == "ok" and c.nominal == pytest.approx(fb * A2 / 2.0, rel=1e-14)
    assert any("G5-2a not credited" in n and "Lv is undefined" in n for n in c.notes)
    # a very long shear span: Eq. G5-2a falls below G5-2b, which then governs (larger of the two)
    c = _checks(thin, Demand("U1", V=1.0e5, Lv=1.0e6))[V]
    assert c.nominal == pytest.approx(fb * A2 / 2.0, rel=1e-14)
    assert gs.round_fcr(E, FY, 250.0, 6.0) == (pytest.approx(fa), pytest.approx(fa), pytest.approx(fb))


def test_hand_g6_minor_axis_shear_of_i_shape():
    # bf = 200, tf = 12: G6: h/tw = bf/(2 tf) = 8.333, kv = 1.2: 1.10 sqrt(1.2 E/Fy) = 1.10 * 26.375 = 29.01
    root = math.sqrt(1.2 * E / FY)
    assert 29.0 < 1.10 * root < 29.02
    # Cv2 = 1.0 (Eq. G2-9); Eq. G6-1 per flange: 0.6 * 345 * 200 * 12 = 496 800 N; two flanges: 993 600 N
    sec = _i("rolled", 400.0, 200.0, 12.0, 8.0, r=18.0)
    c = _checks(sec, Demand("U1", Vz=3.0e5, x_Vz=6000.0))[VZ]
    assert c.nominal == pytest.approx(993600.0, rel=1e-14) and c.capacity == pytest.approx(0.90 * 993600.0, rel=1e-14)
    assert c.clause == "AISC 360-22 G6, Eq. G6-1; G2.2, Eq. G2-9; G1; B3.1 Eq. B3-1"
    assert c.details["Vn per flange = 0.6 Fy bf tf Cv2 (Eq. G6-1)"][0] == pytest.approx(496800.0)
    assert (c.factor, c.x, c.demand) == (0.90, 6000.0, 3.0e5)
    # the phi = 1.00 of G2.1(a) is not used for the flanges
    assert _checks(sec, Demand("U1", Vz=3.0e5), "ASD")[VZ].capacity == pytest.approx(993600.0 / 1.67, rel=1e-14)
    # wide thin flanges, bf = 700, tf = 10: ratio 35 in (29.01, 1.37 * 26.375 = 36.13] -> Eq. G2-10: 29.01/35
    cv2 = 1.10 * root / 35.0
    c = _checks(_i("builtUp", 400.0, 700.0, 10.0, 8.0), Demand("U1", Vz=1.0))[VZ]
    assert c.nominal == pytest.approx(2 * 0.6 * 345.0 * 700.0 * 10.0 * cv2, rel=1e-14) and "Eq. G2-10" in c.clause
    # bf = 800, tf = 10: ratio 40 > 36.13 -> Eq. G2-11: 1.51 * 1.2 * 200000 / (40^2 * 345) = 0.6565
    cv2 = 1.51 * 1.2 * E / (40.0 ** 2 * FY)
    assert 0.656 < cv2 < 0.657
    c = _checks(_i("builtUp", 400.0, 800.0, 10.0, 8.0), Demand("U1", Vz=1.0))[VZ]
    assert c.nominal == pytest.approx(2 * 0.6 * 345.0 * 800.0 * 10.0 * cv2, rel=1e-14) and "Eq. G2-11" in c.clause


def test_shape_dispatch_and_not_checked():
    d = Demand("U", Vy=1.0e4, Vz=2.0e3, V=1.1e4)
    assert list(_checks(_i("rolled", 400.0, 200.0, 12.0, 8.0, r=18.0), d)) == [VY, VZ]
    assert list(_checks(DesignSection.round_hss(A=4000.0, r=100.0, D=300.0, t=5.0), d)) == [V]
    # shape outside the scope: both directions notChecked with the demand kept
    other = DesignSection.unsupported("tee: outside the scope", "T1", A=1000.0)
    c = _checks(other, d)
    assert [c[VY].status, c[VZ].status] == ["notChecked", "notChecked"] and c[VY].demand == 1.0e4
    assert "tee" in c[VZ].reason and c[VY].capacity is None
    res = design.check_member(_ctx(other), [d])
    assert res.passed is False and len([s for s in res.not_checked if s.startswith("shear")]) == 2
    # rolled shape without the fillet dimension: h/tw unknown -> web shear notChecked, flanges still checked
    sec = _i("rolled", 400.0, 200.0, 12.0, 8.0)
    c = _checks(sec, d)
    assert c[VY].status == "notChecked" and "h/tw" in c[VY].reason and c[VZ].status == "ok"
    # ratios only (no d, tw / bf, tf / t): areas cannot be formed
    sec = DesignSection.i_shape("rolled", A=1.0, rx=1.0, ry=1.0, flange_ratio=8.0, web_ratio=30.0)
    c = _checks(sec, d)
    assert c[VY].status == "notChecked" and c[VZ].status == "notChecked"
    sec = DesignSection.rectangular("rectHSS", A=1.0, rx=1.0, ry=1.0, h_ratio=30.0, b_ratio=20.0)
    assert _checks(sec, d)[VY].status == "notChecked"
    # no material
    assert _checks(DesignSection.round_hss(A=4000.0, r=100.0, D=300.0, t=5.0), d, Fy=None)[V].status == "notChecked"
    # very slender web: note about the F13.2 proportioning limit, still a Chapter G result
    c = _checks(_i("builtUp", 2000.0, 300.0, 20.0, 7.0), d)[VY]
    assert c.status == "ok" and any("260" in n for n in c.notes)


# ============================================================================== tier F (b): Design Examples
def _de(eid, sec, name, demand_key, **demand_kw):
    x = example(eid)
    i, p = x["inputs"], x["printed"]
    out = {}
    for method, key in (("LRFD", "Vu"), ("ASD", "Va")):
        ctx = DesignContext(member=eid, section=sec, method=method, E=i["E"], G=11200.0, Fy=i["Fy"], Fu=None,
                            L=1.0, Lx=1.0, Ly=1.0, Lz=1.0)
        out[method] = {c.name: c for c in gs.check(ctx, Demand(method, **{demand_key: i[key]}, **demand_kw))}[name]
        assert out[method].status == "ok" and out[method].passed            # every example ends "o.k."
    assert_printed(out["LRFD"].nominal, p["Vn"], f"{eid} Vn")
    assert_printed(out["LRFD"].capacity, p["phiVn"], f"{eid} phi Vn")
    assert_printed(out["ASD"].capacity, p["Vn_over_Omega"], f"{eid} Vn/Omega")
    return out, i, p


def test_de_g1b_w24x62_major_axis_shear_g21a():
    """Example G.1B (G-4, PDF 240): d = 23.7, tw = 0.430; Cv1 = 1.0, Vn = 306 kips, phi_v = 1.00 -> 306,
    Omega_v = 1.50 -> 204 kips.  The example does not print h/tw (it states that G2.1(a) applies): the test
    supplies h/tw at the limit 2.24 sqrt(E/Fy), the largest value for which the statement is true."""
    i = example("G.1B")["inputs"]
    at_limit = 2.24 * math.sqrt(i["E"] / i["Fy"])
    sec = DesignSection.i_shape("rolled", A=1.0, rx=1.0, ry=1.0, d=i["d"], tw=i["tw"], web_ratio=at_limit,
                                sources={"web_ratio": "given: G2.1(a) limit (not printed in the example)"})
    out, i, p = _de("G.1B", sec, VY, "Vy")
    assert (out["LRFD"].factor, out["ASD"].factor) == (p["phi"], p["Omega"])
    assert "G2.1(a), Eq. G2-1, G2-2" in out["LRFD"].clause
    assert out["LRFD"].inputs["h/tw"][2].startswith("given: G2.1(a) limit")


def test_de_g4_rectangular_hss_shear():
    """Example G.4, HSS6x4x3/8 (G-9..G-10, PDF 245-246): h = 4.95, h/t = 14.2, Aw = 3.46, Vn = 104 ->
    93.6 / 62.3 kips (rounded chain; full precision 103.7 / 93.3 / 62.1, inside the policy)."""
    i = example("G.4")["inputs"]
    sec = DesignSection.rectangular("rectHSS", A=1.0, rx=1.0, ry=1.0, H=i["H"], B=i["B"], t=i["t"])
    out, i, p = _de("G.4", sec, VY, "Vy")
    d = out["LRFD"].details
    assert_printed(d["h"][0], p["h"], "G.4 h")
    assert_printed(d["h/t"][0], p["h_t"], "G.4 h/t")
    assert_printed(d["Aw = 2 h t"][0], p["Aw"], "G.4 Aw")
    assert d["Cv2 (Eq. G2-9)"][0] == 1.0 and d["kv"][0] == 5.0


def test_de_g5_round_hss_shear_with_documented_asd_value():
    """Example G.5, HSS16.000x0.375 (G-11..G-12, PDF 247-248): Fcr 112 (G5-2a), 73.0 (G5-2b), 0.6 Fy = 30.0
    controls; Vn = 258 -> 232 / 155 kips.  Documented printed inconsistency: 258 / 1.67 = 154.49 (-> 154),
    printed 155; inside the policy (0.3 %)."""
    i = example("G.5")["inputs"]
    sec = DesignSection.round_hss(A=i["Ag"], r=1.0, D=i["D"], ratio=i["D_t"])
    out, i, p = _de("G.5", sec, V, "V", Lv=i["Lv"])
    d = out["LRFD"].details
    assert_printed(d["Fcr (Eq. G5-2a)"][0], p["Fcr_a"], "G.5 Fcr (G5-2a)")
    assert_printed(d["Fcr (Eq. G5-2b)"][0], p["Fcr_b"], "G.5 Fcr (G5-2b)")
    assert_printed(d["Fcr"][0], p["Fcr"], "G.5 Fcr")
    # the engine value for ASD is the arithmetic one
    assert float(f"{out['ASD'].capacity:.4g}") == 154.5
    assert float(f"{p['Vn'] / 1.67:.4g}") == 154.5 and p["Vn_over_Omega"] == 155


def test_de_g6_w21x48_minor_axis_shear():
    """Example G.6 (G-13..G-14, PDF 249-250): bf = 8.14, tf = 0.430; bf/2tf = 9.47; Cv2 = 1.0;
    Vn = 210 (two flanges) -> 189 / 126 kips."""
    i = example("G.6")["inputs"]
    sec = DesignSection.i_shape("rolled", A=1.0, rx=1.0, ry=1.0, bf=i["bf"], tf=i["tf"])
    out, i, p = _de("G.6", sec, VZ, "Vz")
    d = out["LRFD"].details
    assert_printed(d["bf/(2 tf)"][0], p["bf_2tf"], "G.6 bf/2tf")
    assert d["Cv2 (Eq. G2-9)"][0] == 1.0 and d["kv"][0] == 1.2
    assert (out["LRFD"].factor, out["ASD"].factor) == (0.90, 1.67)


def test_de_g8b_second_panel_cv2_by_eq_g2_11():
    """Example G.8B (G-23..G-24, PDF 259-260), second panel of the built-up girder: kv = 5 + 5/(a/h)^2 = 5.67
    (Eq. G2-5, a = 90.0, h = 33.0), 1.37 sqrt(kv E/Fy) = 78.6 < h/tw = 106, so Cv2 = 1.51 kv E/((h/tw)^2 Fy)
    = 0.442 (Eq. G2-11).  The engine's function for Eq. G2-9 to G2-11 is called with the numbers the example
    carries; the example's strength (tension field, Eq. G2-7) is outside the engine and is not reproduced."""
    x = example("G.8B")
    i, p = x["inputs"], x["printed"]
    assert_printed(5.0 + 5.0 / (i["a"] / i["h"]) ** 2, p["kv"], "G.8B kv (Eq. G2-5)")
    assert_printed(1.37 * math.sqrt(i["kv"] * i["E"] / i["Fy"]), p["Cv2_limit_1.37"], "G.8B 1.37 sqrt(kv E/Fy)")
    value, eq = gs.cv2(i["h_tw"], i["kv"], i["E"], i["Fy"])
    assert eq == "G2-11"
    assert_printed(value, p["Cv2"], "G.8B Cv2 (Eq. G2-11)")


# ============================================================================== results -> Lv -> G5
def test_round_member_shear_from_results_uses_resultant_and_lv():
    data = {
        "formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": 6000, "y": 0, "z": 0}},
        "supports": {"A": {"type": "fixed"}, "B": {"type": "pinned"}},
        "materials": {"S": {"E": 200000, "nu": 0.3, "Fy": 345, "Fu": 450}},
        "sections": {"thin": {"type": "pipe", "D": 500, "t": 2}},
        "members": {"beam": {"i": "A", "j": "B", "section": "thin", "material": "S"}},
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fz": -1}]}},
        "combinations": {"U1": {"factors": {"D": 1.4}}},
    }
    m = build_model(data)
    x = np.linspace(0.0, 6000.0, 11)
    f = np.zeros((11, 6))
    f[:, 1] = 3.0e4 * (1.0 - x / 3000.0)          # Vy: simply supported, uniform load: zero at midspan
    f[:, 2] = 4.0e4 * (1.0 - x / 3000.0)          # Vz likewise -> resultant 50 kN at the ends, Lv = 3000
    mr = MemberResult(np.concatenate([-f[0], f[-1]]), x, f, np.zeros((11, 3)))
    res = Results("m", {"U1": CombinationResult("U1", "firstOrder", np.zeros((2, 6)), np.zeros((2, 6)), {"beam": mr})})
    out = design.check(m, res)
    c = {c.name: c for c in out.members["beam"].checks}[V]
    A = m.sections["thin"].A
    fa = 1.60 * E / (math.sqrt(3000.0 / 500.0) * 250.0 ** 1.25)     # Eq. G5-2a (131.4 MPa, see the hand test)
    assert c.demand == pytest.approx(5.0e4) and c.inputs["Lv"][0] == pytest.approx(3000.0)
    assert c.nominal == pytest.approx(fa * A / 2.0, rel=1e-12) and c.combination == "U1"
    assert c.ratio == pytest.approx(5.0e4 / (0.9 * fa * A / 2.0), rel=1e-12)
    assert any("0.93" in w for w in out.members["beam"].warnings)   # B4.2 statement travels with the member
