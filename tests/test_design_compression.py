"""AISC 360-22 Chapter E (compression): E2 slenderness (advisory), E3 flexural
buckling about both axes, E4 torsional buckling (Eq. E4-2), E7 slender elements.

Tier F: (a) hand calculations from the Specification (pp. 16.1-38 to 16.1-49) as
explicit arithmetic; (b) AISC Design Examples v16.0 E.1C, E.1D, E.1E, E.2, E.3,
E.9, E.10, E.11 (tests/design_examples; ONE tolerance policy).
"""
import math

import pytest

from benchfea import design
from benchfea.design import DesignContext, DesignSection, Demand
from benchfea.design.aisc360_22 import e_compression as ec
from design_examples import FT, assert_printed, example

E, G, FY = 200000.0, 77200.0, 345.0     # MPa
ROOT = math.sqrt(E / FY)
PI2 = math.pi ** 2
FBX, FBY, TB = ec.NAME_FB["x"], ec.NAME_FB["y"], ec.NAME_TB


def _ctx(sec, method="LRFD", **kw):
    base = dict(member="C", section=sec, method=method, E=E, G=G, Fy=FY, Fu=450.0, L=6000.0, Lx=6000.0, Ly=6000.0,
                Lz=6000.0)
    base.update(kw)
    return DesignContext(**base)


def _stocky(**kw):
    """Rolled I with nonslender elements (flange 8 <= 13.48, web 30 <= 35.87): A = 8000, rx = 150, ry = 50."""
    base = dict(A=8000.0, rx=150.0, ry=50.0, flange_ratio=8.0, web_ratio=30.0)
    base.update(kw)
    return DesignSection.i_shape("rolled", **base)


def _by_name(checks):
    return {c.name: c for c in checks}


def _detail(c, start):
    return next(v[0] for k, v in c.details.items() if k.startswith(start))


# ============================================================================== tier F (a): hand calculations
def test_hand_e3_inelastic_flexural_buckling_major_axis():
    # E2: Lcx = Kx Lx = 1.0 * 6000; Lcx/rx = 6000/150 = 40 <= 4.71 sqrt(E/Fy) = 113.40 -> Eq. E3-2
    limit = 4.71 * ROOT
    assert 113.3 < limit < 113.5
    Fe = PI2 * E / 40.0 ** 2                      # Eq. E3-4: 1233.70 MPa
    assert 1233.6 < Fe < 1233.8
    Fn = 0.658 ** (FY / Fe) * FY                  # Eq. E3-2: 0.658^0.27965 * 345 = 306.90 MPa
    assert 306.8 < Fn < 307.0
    Pn = Fn * 8000.0                              # Eq. E3-1
    c = ec.flexural_buckling_major(_ctx(_stocky()), Demand("U1", Pc=1.5e6, x_Pc=0.0))[0]
    assert c.nominal == pytest.approx(Pn, rel=1e-13) and c.capacity == pytest.approx(0.90 * Pn, rel=1e-13)
    assert c.ratio == pytest.approx(1.5e6 / (0.90 * Pn), rel=1e-13) and c.passed
    assert c.clause == "AISC 360-22 E3, Eq. E3-2, E3-4; E2; B4.1, Table B4.1a; B3.1 Eq. B3-1"
    assert c.details["Fe"][0] == pytest.approx(Fe, rel=1e-13) and c.details["Fn"][0] == pytest.approx(Fn, rel=1e-13)
    assert c.details["Lcx/rx"] == (40.0, "dimensionless") and c.inputs["Kx"][0] == 1.0
    assert (c.factor, c.combination, c.x) == (0.90, "U1", 0.0)
    # ASD: Omega_c = 1.67
    a = ec.flexural_buckling_major(_ctx(_stocky(), "ASD"), Demand("U1", Pc=1.0e6))[0]
    assert a.capacity == pytest.approx(Pn / 1.67, rel=1e-13) and (a.factor, a.factor_name) == (1.67, "Omega")


def test_hand_e3_elastic_flexural_buckling_minor_axis_and_k_factor():
    # Lcy/ry = 6000/50 = 120 > 113.40 -> Eq. E3-3: Fn = 0.877 Fe, Fe = pi^2 E / 120^2 = 137.08 MPa
    Fe = PI2 * E / 120.0 ** 2
    assert 137.0 < Fe < 137.1
    Fn = 0.877 * Fe                               # 120.22 MPa
    c = ec.flexural_buckling_minor(_ctx(_stocky()), Demand("U1", Pc=5.0e5))[0]
    assert c.nominal == pytest.approx(Fn * 8000.0, rel=1e-13)
    assert c.clause.startswith("AISC 360-22 E3, Eq. E3-3, E3-4")
    # Ky = 0.8: Lcy = 4800, Lcy/ry = 96 <= 113.40 -> Eq. E3-2 with Fe = pi^2 E / 96^2 = 214.18
    Fe = PI2 * E / 96.0 ** 2
    Fn = 0.658 ** (FY / Fe) * FY
    c = ec.flexural_buckling_minor(_ctx(_stocky(), Ky=0.8, sources={"Ky": "member"}), Demand("U1", Pc=5.0e5))[0]
    assert c.nominal == pytest.approx(Fn * 8000.0, rel=1e-13)
    assert c.details["Lcy = Ky Ly (E2)"][0] == pytest.approx(4800.0) and c.inputs["Ky"] == (0.8, "dimensionless", "member")


def test_hand_commentary_limits_of_e3():
    """4.71 sqrt(E/Fy) and Fe at that slenderness for E = 29,000 ksi, as tabulated to 3 figures in the
    Commentary of E3 (Table C-E3.1, p. 16.1-357): Fy 36 -> 134, 16.0; 50 -> 113, 22.2; 65 -> 99.5, 28.9;
    70 -> 95.9, 31.1."""
    for Fy, slender, Fe in ((36, 134, 16.0), (50, 113, 22.2), (65, 99.5, 28.9), (70, 95.9, 31.1)):
        lim = 4.71 * math.sqrt(29000.0 / Fy)
        assert float(f"{lim:.3g}") == slender
        assert float(f"{ec.elastic_flexural_stress(29000.0, lim):.3g}") == Fe
        assert ec.branch_by_slenderness(29000.0, Fy, lim) and not ec.branch_by_slenderness(29000.0, Fy, lim * 1.0001)


def test_e3_branch_forms_disagree_only_in_a_sliver_and_it_is_noted():
    # Lc/r = 4.711 sqrt(E/Fy): above 4.71 sqrt(E/Fy) (Eq. E3-3) while Fy/Fe = 4.711^2/pi^2 = 2.2487 <= 2.25
    r = 50.0
    L = 4.711 * ROOT * r
    c = ec.flexural_buckling_minor(_ctx(_stocky(), Ly=L), Demand("U", Pc=1.0))[0]
    Fe = PI2 * E / (L / r) ** 2
    assert FY / Fe < 2.25 and c.nominal == pytest.approx(0.877 * Fe * 8000.0, rel=1e-13)
    assert any("disagree" in n for n in c.notes)
    # the two formulas differ by less than 0.05 % there
    assert 0.877 * Fe == pytest.approx(0.658 ** (FY / Fe) * FY, rel=5e-4)
    assert not any("disagree" in n for n in ec.flexural_buckling_minor(_ctx(_stocky()), Demand())[0].notes)


def test_hand_e4_torsional_buckling_inelastic_and_elastic():
    # Ix = A rx^2 = 1.8e8, Iy = A ry^2 = 2.0e7 mm^4; J = 5.0e5 mm^4, Cw = 4.0e11 mm^6, Lcz = 1.0 * 6000
    sec = _stocky(J=5.0e5, Cw=4.0e11)
    # Eq. E4-2: Fe = [pi^2 E Cw / Lcz^2 + G J] / (Ix + Iy)
    #              = [pi^2 * 2e5 * 4e11 / 3.6e7 + 77200 * 5e5] / 2.0e8 = [2.1932e10 + 3.86e10] / 2e8 = 302.66 MPa
    Fe = (PI2 * E * 4.0e11 / 6000.0 ** 2 + G * 5.0e5) / (1.8e8 + 2.0e7)
    assert 302.6 < Fe < 302.7
    assert FY / Fe <= 2.25                        # 1.1399 -> Eq. E3-2
    Fn = 0.658 ** (FY / Fe) * FY                  # 214.1 MPa
    assert 214.0 < Fn < 214.2
    c = ec.torsional_buckling(_ctx(sec), Demand("U1", Pc=1.0e6))[0]
    assert c.nominal == pytest.approx(Fn * 8000.0, rel=1e-13)        # Eq. E4-1
    assert c.clause == "AISC 360-22 E4, Eq. E4-1, E4-2; E3, Eq. E3-2; B4.1, Table B4.1a; B3.1 Eq. B3-1"
    assert c.details["Fe"][0] == pytest.approx(Fe, rel=1e-13) and c.inputs["G"][0] == G
    assert c.inputs["Cw"][:2] == (4.0e11, "length^6") and c.inputs["J"][:2] == (5.0e5, "length^4")
    assert any("does not exceed" in n for n in c.notes)
    # elastic branch: J = 5.0e4, Kz Lz = 1.0 * 12000:
    # Fe = [pi^2 * 2e5 * 4e11 / 1.44e8 + 77200 * 5e4] / 2e8 = [5.483e9 + 3.86e9] / 2e8 = 46.72 MPa
    Fe = (PI2 * E * 4.0e11 / 12000.0 ** 2 + G * 5.0e4) / 2.0e8
    assert 46.7 < Fe < 46.8 and FY / Fe > 2.25    # 7.38 -> Eq. E3-3
    c = ec.torsional_buckling(_ctx(_stocky(J=5.0e4, Cw=4.0e11), Lz=12000.0), Demand("U1", Pc=1.0e5))[0]
    assert c.nominal == pytest.approx(0.877 * Fe * 8000.0, rel=1e-13) and "E3-3" in c.clause
    assert any("Lz exceeds" in n for n in c.notes)
    # Kz
    c2 = ec.torsional_buckling(_ctx(_stocky(J=5.0e4, Cw=4.0e11), Lz=24000.0, Kz=0.5), Demand("U1", Pc=1.0e5))[0]
    assert c2.nominal == pytest.approx(c.nominal, rel=1e-13)
    # closed section: Cw = 0 -> Fe = G J / (Ix + Iy); round tube J = 2 I -> Fe = G
    tube = DesignSection.round_hss(A=3000.0, r=70.0, ratio=30.0, J=2.0 * 3000.0 * 70.0 ** 2)
    c = ec.torsional_buckling(_ctx(tube), Demand("U1", Pc=1.0))[0]
    assert c.details["Fe"][0] == pytest.approx(G, rel=1e-13) and any("Cw = 0" in n for n in c.notes)


def _slender_web_section():
    # rolled I 400 x 200, tf 12, tw 8, r 18: h = 340, h/tw = 42.5 > 1.49 sqrt(E/Fy) = 35.875 (slender web);
    # flange 100/12 = 8.33 <= 13.48 (nonslender)
    return DesignSection.i_shape("rolled", A=7700.0, rx=160.0, ry=45.0, d=400.0, bf=200.0, tf=12.0, tw=8.0, r=18.0)


def test_hand_e7_slender_web_effective_width():
    sec = _slender_web_section()
    # minor axis, Lcy = 1800: Lcy/ry = 40 -> Fe = 1233.70, Fn = 306.90 MPa (Eq. E3-4, E3-2; gross section)
    Fn = 0.658 ** (FY / (PI2 * E / 40.0 ** 2)) * FY
    lam, lam_r = 340.0 / 8.0, 1.49 * ROOT
    # E7.1(b): lam = 42.5 > lam_r sqrt(Fy/Fn) = 35.875 * 1.0603 = 38.04
    gate = lam_r * math.sqrt(FY / Fn)
    assert 38.0 < gate < 38.1 and lam > gate
    # Table E7.1 case (a) (stiffened element): c1 = 0.18, c2 = 1.31
    Fel = (1.31 * lam_r / lam) ** 2 * FY          # Eq. E7-5: (1.1058)^2 * 345 = 421.9 MPa
    assert 421.5 < Fel < 422.3
    q = math.sqrt(Fel / Fn)
    he = 340.0 * (1.0 - 0.18 * q) * q             # Eq. E7-3: 314.5 mm
    assert 314.0 < he < 315.0
    Ae = 7700.0 - (340.0 - he) * 8.0              # E7 User Note: Ag - (b - be) t = 7496 mm^2
    Pn = Fn * Ae                                  # Eq. E7-1
    c = ec.flexural_buckling_minor(_ctx(sec, Ly=1800.0), Demand("U1", Pc=1.0e6))[0]
    assert c.nominal == pytest.approx(Pn, rel=1e-13) and c.capacity == pytest.approx(0.9 * Pn, rel=1e-13)
    assert c.clause == ("AISC 360-22 E7, Eq. E7-1; E7.1, Eq. E7-2, E7-3, E7-5, Table E7.1; E3, Eq. E3-2, E3-4; E2; "
                        "B4.1, Table B4.1a; B3.1 Eq. B3-1")
    assert _detail(c, "be web") == pytest.approx(he, rel=1e-13) and _detail(c, "Fel web") == pytest.approx(Fel)
    assert _detail(c, "Ae") == pytest.approx(Ae, rel=1e-13) and c.details["class web"][0] == "slender"
    # long column, Lcy/ry = 150 > 113.4: Fn = 0.877 * pi^2 E / 150^2 = 76.94; gate = 35.875 sqrt(345/76.94) = 75.97
    # > 42.5 -> Eq. E7-2, be = b, Ae = Ag: the slender web is fully effective at this stress
    Fn = 0.877 * PI2 * E / 150.0 ** 2
    assert lam_r * math.sqrt(FY / Fn) > lam
    c = ec.flexural_buckling_minor(_ctx(sec, Ly=6750.0), Demand("U1", Pc=1.0e5))[0]
    assert c.nominal == pytest.approx(Fn * 7700.0, rel=1e-13) and "Eq. E7-1; E7.1, Eq. E7-2;" in c.clause
    # each limit state uses its own Fn for Ae: major axis (Lcx/rx = 37.5) has a higher Fn, hence a smaller Ae
    major = ec.flexural_buckling_major(_ctx(sec, Ly=1800.0), Demand("U1", Pc=1.0e6))[0]
    assert major.details["Fn"][0] > 306.9 and _detail(major, "Ae") < Ae


def test_hand_e7_slender_flanges_and_web_of_a_built_up_shape():
    # welded I 600 x 300, tf 10, tw 6: flange b = 150, lam = 15 > 0.64 sqrt(kc E/Fy) = 9.829 (kc = 0.40684);
    # web h = 580, lam = 96.667 > 35.875.  Ag = 2*300*10 + 580*6 = 9480 mm^2
    sec = DesignSection.i_shape("builtUp", A=9480.0, rx=250.0, ry=68.0, d=600.0, bf=300.0, tf=10.0, tw=6.0)
    # minor axis, Lcy = 3400 -> Lcy/ry = 50: Fe = pi^2 E / 2500 = 789.57; Fn = 0.658^(345/789.57) * 345 = 287.3
    Fn = 0.658 ** (FY / (PI2 * E / 50.0 ** 2)) * FY
    assert 287.0 < Fn < 287.6
    kc = 4.0 / math.sqrt(580.0 / 6.0)
    lr_f, lr_w = 0.64 * math.sqrt(kc * E / FY), 1.49 * ROOT
    lam_f, lam_w = 15.0, 580.0 / 6.0
    assert lam_f > lr_f * math.sqrt(FY / Fn) and lam_w > lr_w * math.sqrt(FY / Fn)
    # flange: Table E7.1 case (c): c1 = 0.22, c2 = 1.49
    Fel_f = (1.49 * lr_f / lam_f) ** 2 * FY       # (0.9763)^2 * 345 = 328.9
    qf = math.sqrt(Fel_f / Fn)
    be = 150.0 * (1.0 - 0.22 * qf) * qf           # 122.6 mm
    assert 328.0 < Fel_f < 330.0 and 122.0 < be < 123.5
    # web: case (a): c1 = 0.18, c2 = 1.31
    Fel_w = (1.31 * lr_w / lam_w) ** 2 * FY       # (0.48616)^2 * 345 = 81.54
    qw = math.sqrt(Fel_w / Fn)
    he = 580.0 * (1.0 - 0.18 * qw) * qw           # 279.4 mm
    assert 81.0 < Fel_w < 82.0 and 278.0 < he < 281.0
    Ae = 9480.0 - 4 * (150.0 - be) * 10.0 - (580.0 - he) * 6.0      # four flange outstands, one web
    c = ec.flexural_buckling_minor(_ctx(sec, Ly=3400.0), Demand("U1", Pc=1.0e6))[0]
    assert c.nominal == pytest.approx(Fn * Ae, rel=1e-13)
    assert _detail(c, "be flange") == pytest.approx(be, rel=1e-13) and _detail(c, "be web") == pytest.approx(he, rel=1e-13)
    assert c.details["c1, c2 flange (Table E7.1 case (c))"][0] == "0.22, 1.49"
    assert c.inputs["section fabrication"][0] == "builtUp"


def test_hand_e7_rectangular_hss_walls():
    # HSS 200 x 100 x 5: h = 185 (37.0 > 1.40 sqrt(E/Fy) = 33.708 slender), b = 85 (17, nonslender); Ag = 2800
    sec = DesignSection.rectangular("rectHSS", A=2800.0, rx=72.0, ry=41.0, H=200.0, B=100.0, t=5.0)
    # minor axis, Lcy = 2050 -> Lcy/ry = 50 -> Fn = 287.3 (as above)
    Fn = 0.658 ** (FY / (PI2 * E / 50.0 ** 2)) * FY
    lam, lam_r = 37.0, 1.40 * ROOT
    assert lam > lam_r * math.sqrt(FY / Fn)       # 37 > 33.708 * 1.0958 = 36.94
    # Table E7.1 case (b) (walls of rectangular HSS): c1 = 0.20, c2 = 1.38
    Fel = (1.38 * lam_r / lam) ** 2 * FY          # (1.2572)^2 * 345 = 545.3
    q = math.sqrt(Fel / Fn)
    he = 185.0 * (1.0 - 0.20 * q) * q             # 184.6 mm
    assert 545.0 < Fel < 546.0 and 184.0 < he < 185.0
    Ae = 2800.0 - 2 * (185.0 - he) * 5.0          # two depth walls
    c = ec.flexural_buckling_minor(_ctx(sec, Ly=2050.0), Demand("U1", Pc=1.0e5))[0]
    assert c.nominal == pytest.approx(Fn * Ae, rel=1e-13)
    assert c.details["c1, c2 wall (depth) (Table E7.1 case (b))"][0] == "0.2, 1.38"
    # the same outside dimensions as a built-up BOX: clear widths 190 / 90, Case 8, Table E7.1 case (a)
    box = DesignSection.rectangular("box", A=2900.0, rx=72.0, ry=41.0, H=200.0, B=100.0, t=5.0)
    lam, lam_r = 38.0, 1.49 * ROOT
    # at Fn = 287.3: lam_r sqrt(Fy/Fn) = 35.875 * 1.0958 = 39.31 > 38 -> Eq. E7-2, fully effective, Ae = Ag
    assert lam <= lam_r * math.sqrt(FY / Fn)
    c = ec.flexural_buckling_minor(_ctx(box, Ly=2050.0), Demand("U1", Pc=1.0e5))[0]
    assert c.nominal == pytest.approx(Fn * 2900.0, rel=1e-13) and c.details["class wall (depth)"][0] == "slender"
    # shorter column, Lcy/ry = 20: Fe = pi^2 E / 400 = 4934.8, Fn = 0.658^(0.06991) * 345 = 335.05;
    # lam_r sqrt(Fy/Fn) = 35.875 * 1.01474 = 36.40 < 38 -> Eq. E7-3 with c1 = 0.18, c2 = 1.31
    Fn = 0.658 ** (FY / (PI2 * E / 20.0 ** 2)) * FY
    assert 335.0 < Fn < 335.1 and lam > lam_r * math.sqrt(FY / Fn)
    Fel = (1.31 * lam_r / lam) ** 2 * FY          # (1.23674)^2 * 345 = 527.7
    q = math.sqrt(Fel / Fn)
    he = 190.0 * (1.0 - 0.18 * q) * q             # 184.5 mm
    assert 527.0 < Fel < 528.5 and 184.0 < he < 185.0
    c = ec.flexural_buckling_minor(_ctx(box, Ly=820.0), Demand("U1", Pc=1.0e5))[0]
    assert c.nominal == pytest.approx(Fn * (2900.0 - 2 * (190.0 - he) * 5.0), rel=1e-13)
    assert c.details["c1, c2 wall (depth) (Table E7.1 case (a))"][0] == "0.18, 1.31"


def test_e7_effective_width_is_never_more_than_b_and_continuous_enough():
    """With the rounded c2 of Table E7.1, Eq. E7-3 gives slightly more than b just above the limit of
    Eq. E7-2 ((1 - 0.18 * 1.31) * 1.31 = 1.0011): be is taken as b there."""
    lam_r, Fy, Fn = 35.0, 345.0, 300.0
    at = lam_r * math.sqrt(Fy / Fn)
    assert ec.effective_width(100.0, at, lam_r, Fy, Fn, 0.18, 1.31) == (100.0, None, "E7-2")
    be, Fel, eq = ec.effective_width(100.0, at * 1.0001, lam_r, Fy, Fn, 0.18, 1.31)
    assert eq == "E7-3" and be == 100.0 and (1.0 - 0.18 * 1.31) * 1.31 > 1.0
    be, _, _ = ec.effective_width(100.0, at * 1.2, lam_r, Fy, Fn, 0.18, 1.31)
    assert be < 100.0


def test_hand_e7_round_hss():
    # D/t = 80: 0.11 E/Fy = 63.77 < 80 < 0.45 E/Fy = 260.87 -> Eq. E7-7:
    # Ae = [0.038 E / (Fy (D/t)) + 2/3] Ag = [0.27536 + 0.66667] Ag = 0.94203 Ag
    factor = 0.038 * E / (FY * 80.0) + 2.0 / 3.0
    assert 0.9420 < factor < 0.9421
    sec = DesignSection.round_hss(A=6200.0, r=140.0, D=400.0, t=5.0)
    Fn = 0.658 ** (FY / (PI2 * E / (6000.0 / 140.0) ** 2)) * FY
    c = ec.flexural_buckling_major(_ctx(sec), Demand("U1", Pc=1.0e6))[0]
    assert c.nominal == pytest.approx(Fn * factor * 6200.0, rel=1e-13) and "E7.2, Eq. E7-7" in c.clause
    assert c.details["class wall"][0] == "slender"
    # D/t = 64 (just above 63.77): Eq. E7-7 gives 0.038 * 579.71/64 + 2/3 = 1.0109 Ag -> Ae is not taken above Ag
    assert 0.038 * E / (FY * 64.0) + 2.0 / 3.0 > 1.0
    assert ec.round_effective_area(1000.0, 64.0, E, FY) == (1000.0, "E7-7")
    sec = DesignSection.round_hss(A=6200.0, r=140.0, ratio=64.0)
    c = ec.flexural_buckling_major(_ctx(sec), Demand("U1", Pc=1.0e6))[0]
    assert c.nominal == pytest.approx(Fn * 6200.0, rel=1e-13) and any("Ae taken as Ag" in n for n in c.notes)
    # D/t <= 0.11 E/Fy: nonslender, Eq. E3-1
    assert ec.round_effective_area(1000.0, 63.0, E, FY) == (1000.0, "E7-6")
    # D/t >= 0.45 E/Fy = 260.87: no provision -> notChecked (never a number)
    sec = DesignSection.round_hss(A=6200.0, r=140.0, ratio=261.0)
    c = ec.flexural_buckling_major(_ctx(sec), Demand("U1", Pc=1.0))[0]
    assert (c.status, c.capacity, c.passed) == ("notChecked", None, None) and "0.45 E/Fy" in c.reason
    assert ec.round_effective_area(1000.0, 0.45 * E / FY, E, FY)[0] is None      # the limit itself is excluded


def test_hand_e2_slenderness_is_advisory():
    # Lcx/rx = 6000/150 = 40, Lcy/ry = 6000/50 = 120 -> 120 <= 200 (E2 User Note); ratio 0.6
    c = ec.slenderness(_ctx(_stocky()), Demand())[0]
    assert (c.demand, c.capacity, c.ratio, c.passed, c.advisory) == (120.0, 200.0, 0.6, True, True)
    assert c.clause == "AISC 360-22 E2 (User Note)" and c.details["governing axis"][0] == "y"
    c = ec.slenderness(_ctx(_stocky(), Ky=2.0), Demand())[0]
    assert c.demand == 240.0 and c.passed is False
    res = design.check_member(_ctx(_stocky(J=5.0e5, Cw=4.0e11), Ky=2.0), [Demand("U", Pc=1.0e4)])
    assert res.passed is True                     # the advisory check does not fail the member


def test_not_checked_cases_are_explicit():
    d = Demand("U", Pc=1.0e5)
    # no J / Cw: torsional buckling cannot be evaluated -> notChecked, and the member is not "passed"
    res = design.check_member(_ctx(_stocky()), [d])
    tb = _by_name(res.checks)[TB]
    assert tb.status == "notChecked" and "J" in tb.reason and "Cw" in tb.reason and res.passed is False
    assert _by_name(res.checks)[FBX].status == "ok"
    # slender element known only by its ratio and without a thickness: the area reduction cannot be computed
    sec = DesignSection.i_shape("rolled", A=7700.0, rx=160.0, ry=45.0, flange_ratio=8.0, web_ratio=42.5)
    c = ec.flexural_buckling_minor(_ctx(sec, Ly=1800.0), d)[0]
    assert c.status == "notChecked" and "thickness" in c.reason
    # ... but when it is fully effective at the column stress no dimension is needed
    assert ec.flexural_buckling_minor(_ctx(sec, Ly=6750.0), d)[0].status == "ok"
    # missing material data, missing G, out-of-scope shape
    assert ec.flexural_buckling_major(_ctx(_stocky(), E=None), d)[0].status == "notChecked"
    c = ec.torsional_buckling(_ctx(_stocky(J=1.0, Cw=1.0), G=None), d)[0]
    assert c.status == "notChecked" and "G" in c.reason
    other = DesignSection.unsupported("channel: outside the scope", "C1", A=1000.0, Ix=1.0e6, Iy=1.0e5)
    checks = ec.check(_ctx(other), d)
    assert [c.status for c in checks[:3]] == ["notChecked"] * 3 and all("channel" in c.reason for c in checks[:3])
    assert all(c.demand == 1.0e5 for c in checks[:3])
    # rolled shape without the fillet dimension: classification impossible -> notChecked
    sec = DesignSection.i_shape("rolled", A=7700.0, rx=160.0, ry=45.0, d=400.0, bf=200.0, tf=12.0, tw=8.0)
    assert "web depth h not available" in ec.flexural_buckling_major(_ctx(sec), d)[0].reason


# ============================================================================== tier F (b): Design Examples
def _de_ctx(eid, sec, method, Lx, Ly, Lz=None, G=11200.0):
    i = example(eid)["inputs"]
    return DesignContext(member=eid, section=sec, method=method, E=i["E"], G=i.get("G", G), Fy=i["Fy"], Fu=None,
                         L=max(Lx, Ly), Lx=Lx, Ly=Ly, Lz=Lz if Lz is not None else Ly)


def _de_strengths(eid, sec, name, printed, Lx, Ly, Lz=None, demands=(0.0, 0.0)):
    """LRFD and ASD check ``name``; compares Pn, phi Pn, Pn/Omega with the printed values (policy)."""
    lr = _by_name(ec.check(_de_ctx(eid, sec, "LRFD", Lx, Ly, Lz), Demand("LRFD", Pc=demands[0])))[name]
    asd = _by_name(ec.check(_de_ctx(eid, sec, "ASD", Lx, Ly, Lz), Demand("ASD", Pc=demands[1])))[name]
    assert lr.status == "ok", lr.reason
    if "Pn" in printed:
        assert_printed(lr.nominal, printed["Pn"], f"{eid} Pn")
    assert_printed(lr.capacity, printed["phiPn"], f"{eid} phi Pn")
    assert_printed(asd.capacity, printed["Pn_over_Omega"], f"{eid} Pn/Omega")
    assert (lr.factor, asd.factor) == (0.90, 1.67)
    return lr, asd


def test_de_e1d_w14x90_flexural_buckling_e3():
    """Example E.1D (E-10..E-12, PDF 74-76): Lcx = 30 ft, Lcy = 15 ft; x-x governs; Fe 83.3, Fn 38.9 ksi,
    Pn 1,030 -> 927 / 617 kips.  Listed: Ag, rx, ry, bf/2tf, h/tw."""
    x = example("E.1D")
    i, p = x["inputs"], x["printed"]
    sec = DesignSection.i_shape("rolled", A=i["Ag"], rx=i["rx"], ry=i["ry"], flange_ratio=i["bf_2tf"],
                                web_ratio=i["h_tw"])
    Lx, Ly, Lz = i["Lcx_ft"] * FT, i["Lcy_ft"] * FT, i["Lcz_ft"] * FT
    lr, asd = _de_strengths("E.1D", sec, FBX, p, Lx, Ly, Lz, (i["Pu"], i["Pa"]))
    assert_printed(lr.details["Lcx/rx"][0], p["Lcx_over_rx"], "E.1D Lcx/rx")
    assert_printed(lr.details["Fe"][0], p["Fe"], "E.1D Fe")
    assert_printed(lr.details["Fn"][0], p["Fn"], "E.1D Fn")
    assert lr.details["Fn equation"][0] == "E3-2" and lr.passed and asd.passed
    minor = _by_name(ec.check(_de_ctx("E.1D", sec, "LRFD", Lx, Ly, Lz), Demand("LRFD", Pc=i["Pu"])))[FBY]
    assert_printed(minor.details["Lcy/ry"][0], p["Lcy_over_ry"], "E.1D Lcy/ry")
    assert minor.ratio < lr.ratio                 # "x-x buckling controls"
    # the example does not evaluate torsional buckling and lists no J, Cw: the engine reports it notChecked
    res = design.check_member(_de_ctx("E.1D", sec, "LRFD", Lx, Ly, Lz), [Demand("LRFD", Pc=i["Pu"])])
    assert _by_name(res.checks)[TB].status == "notChecked" and res.governing.name == FBX


def test_de_e1c_w14x132_flexural_buckling_e3():
    """Example E.1C (E-8..E-9, PDF 72-73): Lc = 30 ft both axes; 892 / 598 kips.  The example takes the
    available stresses from Manual Table 4-14; documented case: Eq. E3-2 gives Pn/Omega = 594 kips (0.6 %
    below the printed 598, inside the policy)."""
    x = example("E.1C")
    i, p = x["inputs"], x["printed"]
    sec = DesignSection.i_shape("rolled", A=i["Ag"], rx=i["rx"], ry=i["ry"], flange_ratio=i["bf_2tf"],
                                web_ratio=i["h_tw"])
    lr, asd = _de_strengths("E.1C", sec, FBY, p, i["Lcx_ft"] * FT, i["Lcy_ft"] * FT, demands=(i["Pu"], i["Pa"]))
    assert_printed(lr.details["Lcy/ry"][0], p["Lcy_over_ry"], "E.1C Lcy/ry")
    assert float(f"{asd.capacity:.3g}") == 594.0  # the equation value, as quoted from Table 4-1a in Example E.1A


@pytest.mark.parametrize("case", ["a", "b", "c"])
def test_de_e1e_w16x31_slender_web_e7(case):
    """Example E.1E (E-13..E-19, PDF 77-83): W16x31, Lc = 5, 10, 15 ft; printed Pn 349 / 210 / 96.8 kips.
    Listed: Ag, bf, d, tf, tw, rx, ry, kdes, bf/2tf, h/tw.  Case a reduces the web (Eq. E7-3), b and c do not."""
    x = example("E.1E")
    i, p = x["inputs"], x["printed"][case]
    sec = DesignSection.i_shape("rolled", A=i["Ag"], rx=i["rx"], ry=i["ry"], d=i["d"], bf=i["bf"], tf=i["tf"],
                                tw=i["tw"], k=i["kdes"], flange_ratio=i["bf_2tf"], web_ratio=i["h_tw"])
    L = p["Lc_ft"] * FT
    lr, _ = _de_strengths("E.1E", sec, FBY, p, L, L)
    assert_printed(lr.details["Lcy/ry"][0], p["Lc_over_r"], "E.1E Lc/r")
    assert_printed(lr.details["Fe"][0], p["Fe"], "E.1E Fe")
    assert_printed(lr.details["Fn"][0], p["Fn"], "E.1E Fn")
    assert_printed(_detail(lr, "lambda_r sqrt(Fy/Fn) web"), p["gate"], "E.1E lambda_r sqrt(Fy/Fn)")
    assert_printed(_detail(lr, "Ae"), p["Ae"], "E.1E Ae")
    assert lr.details["class web"][0] == "slender" and "E7, Eq. E7-1" in lr.clause
    assert lr.details["Fn equation"][0] == ("E3-3" if case == "c" else "E3-2")
    if case == "a":
        assert_printed(_detail(lr, "Fel web"), p["Fel"], "E.1E Fel")
        assert_printed(_detail(lr, "be web"), p["he"], "E.1E he")
        assert lr.details["c1, c2 web (Table E7.1 case (a))"][0] == "0.18, 1.31"
    else:
        assert not any(k.startswith("be web") for k in lr.details)       # Eq. E7-2: fully effective


def _built_up(eid):
    i = example(eid)["inputs"]
    return i, DesignSection.i_shape("builtUp", A=i["Ag"], Ix=i["Ix"], Iy=i["Iy"], ry=i["ry"], J=i["J"], Cw=i["Cw"],
                                    d=i["d"], bf=i["bf"], tf=i["tf"], tw=i["tw"], h=i["h"])


def test_de_e2_built_up_column_slender_web_e7_and_e4_stress():
    """Example E.2 (E-20..E-24, PDF 84-88): built-up I, slender web; Fe flexural 38.3, torsional (Eq. E4-2)
    91.9 ksi; Fel 30.7, he 12.6, Ae 19.2; Pn 555 -> 500 / 332 kips."""
    i, sec = _built_up("E.2")
    p = example("E.2")["printed"]
    L = i["Lc_ft"] * FT
    lr, asd = _de_strengths("E.2", sec, FBY, p, L, L, L, (i["Pu"], i["Pa"]))
    assert_printed(lr.details["Lcy/ry"][0], p["Lcy_over_ry"], "E.2 Lcy/ry")
    assert_printed(lr.details["Fe"][0], p["Fe_flexural"], "E.2 Fe (E3-4)")
    assert_printed(lr.details["Fn"][0], p["Fn"], "E.2 Fn")
    assert_printed(_detail(lr, "lambda_r sqrt(Fy/Fn) web"), p["gate"], "E.2 gate")
    assert_printed(_detail(lr, "Fel web"), p["Fel"], "E.2 Fel")
    assert_printed(_detail(lr, "be web"), p["he"], "E.2 he")
    assert_printed(_detail(lr, "Ae"), p["Ae"], "E.2 Ae")
    assert lr.passed and asd.passed
    # torsional buckling: the example prints the elastic stress of Eq. E4-2 and finds that flexure controls
    tb = _by_name(ec.check(_de_ctx("E.2", sec, "LRFD", L, L, L), Demand("LRFD", Pc=i["Pu"])))[TB]
    assert tb.status == "ok" and tb.details["Fe equation"][0] == "E4-2"
    assert_printed(tb.details["Fe"][0], p["Fe_torsional"], "E.2 Fe (E4-2)")
    assert tb.nominal > lr.nominal
    res = design.check_member(_de_ctx("E.2", sec, "LRFD", L, L, L), [Demand("LRFD", Pc=i["Pu"])])
    assert res.passed and res.governing.name == FBY
    assert all("rupture" in s for s in res.not_checked)      # (no Fu in a compression example; zero tension)


def test_de_e3_built_up_column_slender_flanges_fully_effective():
    """Example E.3 (E-25..E-29, PDF 89-93): slender flanges (14.0 > 13.3) but be = b at Fn = 36.4 ksi
    (Eq. E7-2); Fe flexural 65.9, torsional (Eq. E4-2) 71.2 ksi; Pn 353 -> 318 / 211 kips."""
    i, sec = _built_up("E.3")
    p = example("E.3")["printed"]
    L = i["Lc_ft"] * FT
    lr, asd = _de_strengths("E.3", sec, FBY, p, L, L, L, (i["Pu"], i["Pa"]))
    assert_printed(lr.details["Lcy/ry"][0], p["Lcy_over_ry"], "E.3 Lcy/ry")
    assert_printed(lr.details["Fe"][0], p["Fe_flexural"], "E.3 Fe (E3-4)")
    assert_printed(lr.details["Fn"][0], p["Fn"], "E.3 Fn")
    assert_printed(_detail(lr, "lambda_r sqrt(Fy/Fn) flange"), p["gate"], "E.3 gate")
    assert lr.details["class flange"][0] == "slender" and _detail(lr, "Ae") == i["Ag"]
    assert "E7, Eq. E7-1; E7.1, Eq. E7-2;" in lr.clause
    tb = _by_name(ec.check(_de_ctx("E.3", sec, "LRFD", L, L, L), Demand("LRFD", Pc=i["Pu"])))[TB]
    assert_printed(tb.details["Fe"][0], p["Fe_torsional"], "E.3 Fe (E4-2)")
    assert tb.nominal > lr.nominal                # "use Fe = 65.9 ksi": flexural buckling controls


def test_de_e9_rectangular_hss_without_slender_elements():
    """Example E.9, HSS12x10x3/8 (E-59..E-61, PDF 123-125): Lc = 16 ft; Fe 125, Fn 42.3 ksi; Pn 618 ->
    556 / 370 kips.  Listed: Ag, rx, ry, b/t, h/t."""
    x = example("E.9")
    i, p = x["inputs"], x["printed"]
    sec = DesignSection.rectangular("rectHSS", A=i["Ag"], rx=i["rx"], ry=i["ry"], h_ratio=i["h_t"], b_ratio=i["b_t"])
    L = i["Lc_ft"] * FT
    lr, asd = _de_strengths("E.9", sec, FBY, p, L, L, demands=(i["Pu"], i["Pa"]))
    assert_printed(lr.details["Lcy/ry"][0], p["Lcy_over_ry"], "E.9 Lcy/ry")
    assert_printed(lr.details["Fe"][0], p["Fe"], "E.9 Fe")
    assert_printed(lr.details["Fn"][0], p["Fn"], "E.9 Fn")
    assert "E3, Eq. E3-2" in lr.clause and "E7" not in lr.clause and lr.passed and asd.passed


@pytest.mark.parametrize("case", ["18", "24", "40"])
def test_de_e10_rectangular_hss_with_slender_walls_e7(case):
    """Example E.10, HSS12x8x3/16 (E-62..E-66, PDF 126-130; Table E.10-1): Lc = 18, 24, 40 ft; printed Pn
    197 / 168 / 83.1 kips.  Listed: Ag, t, rx, ry, b/t, h/t (flat widths = ratio x t, as in the example)."""
    x = example("E.10")
    i, p = x["inputs"], x["printed"][case]
    sec = DesignSection.rectangular("rectHSS", A=i["Ag"], rx=i["rx"], ry=i["ry"], t=i["t"], h_ratio=i["h_t"],
                                    b_ratio=i["b_t"])
    L = p["Lc_ft"] * FT
    lr, _ = _de_strengths("E.10", sec, FBY, p, L, L)
    assert_printed(lr.details["Fn"][0], p["Fn"], "E.10 Fn")
    assert_printed(_detail(lr, "lambda_r sqrt(Fy/Fn) wall (depth)"), p["gate"], "E.10 gate")
    assert_printed(_detail(lr, "Ae"), p["Ae"], "E.10 Ae")
    assert lr.details["class wall (depth)"][0] == lr.details["class wall (width)"][0] == "slender"
    if "he" in p:
        assert_printed(_detail(lr, "Fel wall (depth)"), p["Fel_wide"], "E.10 Fel wide wall")
        assert_printed(_detail(lr, "be wall (depth)"), p["he"], "E.10 he")
        assert lr.details["c1, c2 wall (depth) (Table E7.1 case (b))"][0] == "0.2, 1.38"
    if "be" in p:
        assert_printed(_detail(lr, "Fel wall (width)"), p["Fel_narrow"], "E.10 Fel narrow wall")
        assert_printed(_detail(lr, "be wall (width)"), p["be"], "E.10 be")
    assert any(k.startswith("be wall (depth)") for k in lr.details) == ("he" in p)
    assert any(k.startswith("be wall (width)") for k in lr.details) == ("be" in p)
    if case == "24":
        assert_printed(lr.details["Lcy/ry"][0], p["Lc_over_r"], "E.10 Lc/r")
        assert_printed(lr.details["Fe"][0], p["Fe"], "E.10 Fe")


def test_de_e11_pipe_compression_member():
    """Example E.11, Pipe 10 Std. (E-67..E-69, PDF 131-133): Fy = 35 ksi, Lcx = 30 ft, Lcy = 15 ft; Fe 29.9,
    Fn 21.4 ksi; Pn 246 -> 221 / 147 kips.  Listed: Ag, r, D/t."""
    x = example("E.11")
    i, p = x["inputs"], x["printed"]
    sec = DesignSection.round_hss(A=i["Ag"], r=i["r"], ratio=i["D_t"])
    lr, asd = _de_strengths("E.11", sec, FBX, p, i["Lcx_ft"] * FT, i["Lcy_ft"] * FT, demands=(i["Pu"], i["Pa"]))
    assert_printed(lr.details["Lcx/rx"][0], p["Lc_over_r"], "E.11 Lc/r")
    assert_printed(lr.details["Fe"][0], p["Fe"], "E.11 Fe")
    assert_printed(lr.details["Fn"][0], p["Fn"], "E.11 Fn")
    assert lr.details["class wall"][0] == "nonslender" and lr.passed and asd.passed
    s = _by_name(ec.check(_de_ctx("E.11", sec, "LRFD", i["Lcx_ft"] * FT, i["Lcy_ft"] * FT), Demand()))[ec.NAME_SLENDER]
    assert s.details["governing axis"][0] == "x" and s.passed     # "x-x controls"
