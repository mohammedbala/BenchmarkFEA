"""Independent review of the AISC 360-22 member checks (benchfea.design).

(a) Spot-check: worked cases of the AISC Design Examples v16.0 that the other design tests do NOT use,
    run through the public API (tests/design_examples/part3; the ONE tolerance policy of
    tests/design_examples/__init__.py, never adjusted per case):
    D.2, D.3, D.6 (tension); E.4A, E.4B (compression); F.14 (built-up girder: F4 -> notChecked, Cb, shear);
    G.8A (built-up web shear); H.5C (rectangular HSS: G4 and F7); A-6.5 (Cb, reverse curvature);
    Part III moment-frame column W14x90 by three analysis methods (E3 + Eq. H1-1a).
(b) Regression tests of the defects found by the review:
    - the interaction of H1 (and the resultant moment of round sections) was evaluated at the stations only,
      so a moment that peaks between two stations was missed (no H1 check at all with two stations);
    - a rectangular section without elastic modulus was notChecked for LTB although Lb <= Lp (Eq. F7-10
      needs no elastic modulus).

Units: the Design Examples cases run in kip, in. (every equation used has dimensionless coefficients);
the regression tests in N, mm.
"""
import math

import numpy as np
import pytest

from benchfea import design
from benchfea.design import DesignContext, DesignSection, Demand
from benchfea.design.aisc360_22 import d_tension as dt, e_compression as ec, f_flexure as ff, g_shear as gs, \
    h_combined as hc
from design_examples.part3 import FT, assert_printed, example

E_KSI, G_KSI = 29000.0, 11200.0
METHODS = ("LRFD", "ASD")


def _ctx(sec, method, L, Fy=50.0, Fu=65.0, **kw):
    base = dict(member="M", section=sec, method=method, E=E_KSI, G=G_KSI, Fy=Fy, Fu=Fu, L=L, Lx=L, Ly=L, Lz=L)
    base.update(kw)
    return DesignContext(**base)


def _by_name(checks):
    return {c.name: c for c in checks}


# ============================================================================== (a) Chapter D
@pytest.mark.parametrize("eid", ["D.2", "D.3", "D.6"])
def test_de_tension_members_outside_the_shape_scope(eid):
    """Examples D.2 (single angle), D.3 (WT), D.6 (double angle): Chapter D needs only Ag (and r for the
    D1 recommendation), so it applies to sections outside the shape scope.  D2-1, D2-2, D3-1."""
    x = example(eid)
    i, p = x["inputs"], x["printed"]
    r = i["r_least"]
    # section outside the shape scope: only A and the radii of gyration (the least one as listed)
    sec = DesignSection.unsupported("shape outside the implemented scope", eid, A=i["Ag"], Ix=i["Ag"] * r * r,
                                    Iy=i["Ag"] * (2.0 * r) ** 2)
    L = (i["L_ft"] if "L_ft" in i else p["Lmax_ft_for_L_over_r_300"]) * FT
    for method, y_key, r_key in (("LRFD", "yield_phiPn", "rupture_phiPn"),
                                 ("ASD", "yield_Pn_over_Omega", "rupture_Pn_over_Omega")):
        c = _by_name(dt.check(_ctx(sec, method, L, Fy=i["Fy"], Fu=i["Fu"], An=i["An"], U=i["U"]),
                              Demand(method, Pt=1.0)))
        y, ru, sl = c[dt.NAME_YIELD], c[dt.NAME_RUPTURE], c[dt.NAME_SLENDER]
        assert (y.status, ru.status, sl.status) == ("ok", "ok", "ok")
        if "yield_Pn" in p:
            assert_printed(y.nominal, p["yield_Pn"], f"{eid} Pn yielding")
        assert_printed(y.capacity, p[y_key], f"{eid} {method} yielding")
        assert_printed(ru.details["Ae = An U (Eq. D3-1)"][0], p["Ae"], f"{eid} Ae")
        assert_printed(ru.nominal, p["rupture_Pn"], f"{eid} Pn rupture")
        assert_printed(ru.capacity, p[r_key], f"{eid} {method} rupture")
        assert "D2-1" in y.clause and "D2-2" in ru.clause and "D3-1" in ru.clause
        assert (y.factor, ru.factor) == ((0.90, 0.75) if method == "LRFD" else (1.67, 2.00))
        assert ru.capacity < y.capacity                       # "governed by the tensile rupture limit state"
        assert_printed(sl.demand, p.get("L_over_r", 300.0), f"{eid} L/r")
        assert sl.advisory


# ============================================================================== (a) Chapter E
def _w14x82(i):
    """W14x82 of Examples E.4A / E.4B: listed Ag, Ix (and rx/ry in E.4A).  The examples list no
    width-to-thickness ratio and read Manual Table 4-1a (no slender-element reduction): the ratios below are
    PLACEHOLDERS for a nonslender section -- for such a section they enter neither Fn nor Pn (Eq. E3-1)."""
    rx = math.sqrt(i["Ix"] / i["Ag"])
    ry = rx / i["rx_over_ry"] if "rx_over_ry" in i else None
    return DesignSection.i_shape("rolled", A=i["Ag"], rx=rx, ry=ry, flange_ratio=1.0, web_ratio=1.0), rx


def test_de_e4b_moment_frame_column_effective_length():
    """Example E.4B (E-34, E-35; PDF 98-99): Lcx = K L = 2.0 x 14 ft = 28 ft, y-axis continuously supported;
    861 / 573 kips.  E2 (Lc = K L), E3-2, E3-4."""
    x = example("E.4B")
    i, p = x["inputs"], x["printed"]
    sec, _ = _w14x82(i)
    for method, key in (("LRFD", "phiPn"), ("ASD", "Pn_over_Omega")):
        c = _by_name(ec.check(_ctx(sec, method, i["L_ft"] * FT, Kx=i["K"], Ly=0.0, Lz=0.0),
                              Demand(method, Pc=1.0)))[ec.NAME_FB["x"]]
        assert c.status == "ok" and c.details["Fn equation"][0] == "E3-2"
        assert_printed(c.details["Lcx = Kx Lx (E2)"][0], p["Lcx_ft"] * FT, "E.4B Lcx")
        assert_printed(c.capacity, p[key], f"E.4B {method}")


def test_de_e4a_moment_frame_column_table_entry():
    """Example E.4A (E-30..E-33; PDF 94-97).  The printed 940 / 626 kips are the Manual Table 4-1a entry at
    Lc = 9 ft about the y-axis (the example rounds the equivalent lengths 8.61 ft and 8.03 ft up), so the
    engine is compared at THAT length: ry = rx / 2.44, Lcy = 9 ft.  At the actual Lcx = 1.5 x 14 = 21 ft the
    equation gives more (the example calls its table reading conservative); this is the example's rounding,
    not a discrepancy of the engine, and it is asserted here as such (no tolerance involved)."""
    x = example("E.4A")
    i, p = x["inputs"], x["printed"]
    sec, rx = _w14x82(i)
    assert_printed(i["K_upper"] * i["L_ft"] / i["rx_over_ry"], p["Lc_equiv_upper_ft"], "E.4A equivalent Lc")
    for method, key in (("LRFD", "phiPn"), ("ASD", "Pn_over_Omega")):
        at_table = _by_name(ec.check(_ctx(sec, method, i["L_ft"] * FT, Lx=0.0, Ly=i["Lc_table_ft"] * FT, Lz=0.0),
                                     Demand(method, Pc=1.0)))[ec.NAME_FB["y"]]
        assert_printed(at_table.capacity, p[key], f"E.4A {method} at the table entry Lc = 9 ft")
        actual = _by_name(ec.check(_ctx(sec, method, i["L_ft"] * FT, Kx=i["K_upper"], Ly=0.0, Lz=0.0),
                                   Demand(method, Pc=1.0)))[ec.NAME_FB["x"]]
        assert_printed(actual.details["Lcx = Kx Lx (E2)"][0], p["Lcx_upper_ft"] * FT, "E.4A Lcx")
        assert actual.capacity > p[key] and actual.capacity > at_table.capacity


# ============================================================================== (a) Part III: E3 + H1-1a
@pytest.mark.parametrize("eid", ["III.MF.DA", "III.MF.ELM", "III.MF.B2"])
@pytest.mark.parametrize("method", METHODS)
def test_de_part3_moment_frame_column_e3_and_h1_1a(eid, method):
    """Part III, interior column W14x90 of the moment frame, three analysis methods (III-84..86, III-90..91,
    III-93..94; PDF 1024-1034): Pc by E3 through the engine, then Eq. H1-1a with the engine's Pc and the
    Mcx the example reads from Manual Table 6-1 (the document lists no flexural property of the W14x90, so
    Chapter F cannot be run for it)."""
    x = example(eid)
    i, p = x["inputs"], x["printed"][method]
    r = i[method]
    sec = DesignSection.i_shape("rolled", A=i["Ag"], rx=i["rx"], ry=i["ry"], flange_ratio=i["bf_2tf"],
                                web_ratio=i["h_tw"])
    c = _by_name(ec.check(_ctx(sec, method, r["Lc_ft"] * FT, Lx=0.0, Lz=0.0), Demand(method, Pc=r["Pr"])))
    minor = c[ec.NAME_FB["y"]]
    assert minor.status == "ok" and minor.details["Fn equation"][0] == "E3-2"
    assert_printed(minor.capacity, p["Pc"], f"{eid} {method} Pc")
    assert_printed(r["Pr"] / minor.capacity, p["Pr_over_Pc"], f"{eid} {method} Pr/Pc")
    value, eq = hc.h1_interaction(r["Pr"], minor.capacity, r["Mrx_kipft"], r["Mcx_kipft"])
    assert eq == "H1-1a"
    assert_printed(value, p["interaction"], f"{eid} {method} Eq. H1-1a")


# ============================================================================== (a) Chapter G
@pytest.mark.parametrize("method,key,demand", [("LRFD", "phiVn", "Vu"), ("ASD", "Vn_over_Omega", "Va")])
def test_de_g8a_built_up_girder_web_without_stiffeners(method, key, demand):
    """Example G.8A, first step (G-17, G-18; PDF 253-254): built-up girder, web 33 x 5/16, no stiffeners:
    kv = 5.34, Cv1 by Eq. G2-4, Vn = 196 kips, 176 / 117 kips (the example concludes n.g.)."""
    x = example("G.8A")
    i, p = x["inputs"], x["printed"]
    A = 2.0 * i["bf"] * i["tf"] + i["h"] * i["tw"]
    sec = DesignSection.i_shape("builtUp", A=A, d=i["d"], bf=i["bf"], tf=i["tf"], tw=i["tw"])
    c = _by_name(gs.check(_ctx(sec, method, 56.0 * FT, Fy=i["Fy"]), Demand(method, Vy=i[demand])))[gs.NAME_VY]
    assert c.status == "ok" and "G2.1(b)" in c.clause and "G2-4" in c.clause
    assert sec.elements["web"].width == pytest.approx(i["h"])          # B4.1b(b): clear distance between flanges
    assert_printed(c.details["Aw = d tw"][0], p["Aw"], "G.8A Aw")
    assert_printed(c.details["h/tw"][0], p["h_tw"], "G.8A h/tw")
    assert c.details["kv (web without transverse stiffeners)"][0] == 5.34
    assert_printed(c.details["1.10 sqrt(kv E/Fy)"][0], p["shear_limit"], "G.8A 1.10 sqrt(kv E/Fy)")
    assert_printed(c.details["Cv1 (Eq. G2-4)"][0], p["Cv1"], "G.8A Cv1")
    assert_printed(c.nominal, p["Vn"], "G.8A Vn")
    assert_printed(c.capacity, p[key], f"G.8A {method}")
    assert c.factor == (0.90 if method == "LRFD" else 1.67) and c.passed is False      # "n.g."


# ============================================================================== (a) F.14: plate girder
def _f14_section(i):
    A = 2.0 * i["bf"] * i["tf"] + i["h"] * i["tw"]
    return DesignSection.i_shape("builtUp", A=A, d=i["d"], bf=i["bf"], tf=i["tf"], tw=i["tw"], Ix=i["Ix"],
                                 Iy=i["Iy"], Sx=i["Sx"], Zx=i["Zx"], J=i["J"], ho=i["ho"])


def _f14_stations(i, method):
    """Stations of the simply supported girder (kip, in.): point load P at midspan, uniform load w.
    Statics: R = P/2 + w L/2; M(x) = R x - w x^2/2 - P <x - L/2>; station forces with dMz/dx = -Vy."""
    r = i[method]
    L, P, w = i["L_ft"] * FT, r["P"], r["w_kip_per_ft"] / FT
    R = P / 2.0 + w * L / 2.0
    xs = np.array([0.0, L / 4.0, L / 2.0, L / 2.0, 3.0 * L / 4.0, L])
    f = np.zeros((len(xs), 6))
    for k, xv in enumerate(xs):
        right = k >= 3                                         # right of the point load
        shear = R - w * xv - (P if right else 0.0)             # dM/dx
        f[k, 1] = -shear
        f[k, 5] = R * xv - w * xv * xv / 2.0 - (P * (xv - L / 2.0) if right else 0.0)
    return xs, f, R


@pytest.mark.parametrize("method,key", [("LRFD", "phiVn"), ("ASD", "Vn_over_Omega")])
def test_de_f14_plate_girder_shear_and_flexure_not_checked(method, key):
    """Example F.14 (F-65..F-73; PDF 227-235).  Shear: G2.1(b), Eq. G2-4, 489 -> 440 / 293 kips.  Flexure:
    the web is noncompact (124 between 90.6 and 137) and the example uses Section F4, which the engine does
    not implement -> the flexure must be notChecked and the member must NOT be reported as passing."""
    x = example("F.14")
    i, p = x["inputs"], x["printed"]
    sec = _f14_section(i)
    xs, f, R = _f14_stations(i, method)
    assert_printed(R, i[method]["Vr"], "F.14 end shear")
    ctx = _ctx(sec, method, i["L_ft"] * FT, Fy=i["Fy"], braceTop=tuple(i["brace_fractions"]),
               braceBottom=tuple(i["brace_fractions"]))
    demand = design.demand_from_stations(xs, f, method)
    res = design.check_member(ctx, [Demand(method, Vy=demand.Vy, x_Vy=demand.x_Vy, extra={"x": xs, "forces": f})])
    c = _by_name(res.checks)
    # --- shear
    v = c[gs.NAME_VY]
    assert v.status == "ok" and "G2-4" in v.clause
    assert_printed(v.details["h/tw"][0], p["h_tw"], "F.14 h/tw")
    assert_printed(v.details["1.10 sqrt(kv E/Fy)"][0], p["shear_limit"], "F.14 1.10 sqrt(kv E/Fy)")
    assert_printed(v.details["Cv1 (Eq. G2-4)"][0], p["Cv1"], "F.14 Cv1")
    assert_printed(v.nominal, p["Vn"], "F.14 Vn")
    assert_printed(v.capacity, p[key], f"F.14 {method} shear")
    assert v.passed
    # --- flexure: classification as printed, Section F4 not implemented
    y = c[ff.NAME_Y["x"]]
    assert y.status == "notChecked" and "F4" in y.clause and "noncompact" in y.reason
    assert y.capacity is None and y.ratio is None and y.passed is None
    assert_printed(y.demand, i[method]["Mr_kipft"] * FT, f"F.14 {method} required flexural strength")
    assert_printed(y.details["lambda_p web (Table B4.1b Case 15)"][0], p["lambda_pw"], "F.14 lambda_pw")
    assert_printed(y.details["lambda_r web (Table B4.1b Case 15)"][0], p["lambda_rw"], "F.14 lambda_rw")
    assert_printed(y.details["lambda flange"][0], p["flange_ratio"], "F.14 bf/2tf")
    assert_printed(y.details["lambda_p flange (Table B4.1b Case 11)"][0], p["lambda_pf"], "F.14 lambda_pf")
    assert y.details["class flange"][0] == "compact" and y.details["class web"][0] == "noncompact"
    assert res.status == "incomplete" and res.passed is False
    assert any("F4" in text for text in res.not_checked)


def test_de_f14_cb_of_the_middle_unbraced_segment():
    """Example F.14 (F-71; PDF 233): Cb = 1.25 for the segment between the quarter point and midspan, with
    MA, MB, MC at its quarter points -- Eq. F1-1 evaluated by the engine from the station forces."""
    x = example("F.14")
    i, p = x["inputs"], x["printed"]
    xs, f, _ = _f14_stations(i, "LRFD")
    mo = ff.moments(Demand("LRFD", extra={"x": xs, "forces": f}))
    L = i["L_ft"] * FT
    Cb, Mmax, MA, MB, MC, x_max = ff.segment_cb(mo.Mx, L / 4.0, L / 2.0)
    assert_printed(Mmax, i["LRFD"]["Mr_kipft"] * FT, "F.14 Mmax")
    assert_printed(MA, p["MA_kipft"] * FT, "F.14 MA")
    assert_printed(MB, p["MB_kipft"] * FT, "F.14 MB")
    assert_printed(MC, p["MC_kipft"] * FT, "F.14 MC")
    assert_printed(Cb, p["Cb"], "F.14 Cb")
    assert x_max == pytest.approx(L / 2.0)


def test_de_a65_cb_reverse_curvature():
    """Example A-6.5 (A6-19, A6-21; PDF 529, 531): Cb = 2.46 for a segment in reverse curvature.  The moments
    are printed as multiples of the uniform load w; with w = 1 the diagram is M(x) = M0 + V0 x - x^2/2 and V0
    follows from the two end moments.  The engine gets the two END stations only (values and slopes) and
    must reproduce the quarter-point moments and Cb (Eq. F1-1, absolute values)."""
    x = example("A-6.5")
    i, p = x["inputs"], x["printed"]
    Ls, M0, M1 = i["segment_ft"], i["M_start_over_w"], i["M_end_over_w"]
    V0 = (M1 - M0 + Ls * Ls / 2.0) / Ls                      # dM/dx at the start (w = 1)
    D = ff.Diagram([0.0, Ls], [M0, M1], [V0, V0 - Ls])
    Cb, Mmax, MA, MB, MC, _ = ff.segment_cb(D, 0.0, Ls)
    assert Mmax == pytest.approx(abs(M0))
    assert_printed(MA, abs(i["M_quarter_over_w"]), "A-6.5 MA")
    assert_printed(MB, abs(i["M_mid_over_w"]), "A-6.5 MB")
    assert_printed(MC, abs(i["M_three_quarter_over_w"]), "A-6.5 MC")
    assert_printed(Cb, p["Cb"], "A-6.5 Cb")
    assert Cb == pytest.approx(ff.cb_f1_1(M0, MA, MB, MC), rel=1e-14)


# ============================================================================== (a) H.5C: rectangular HSS
def _h5c_section(i):
    """HSS6x4x1/4 with the listed properties.  The example lists neither Ix nor Sx; the engine needs to know
    which axis is the stronger one (F7.4), so Ix and Iy of the sharp-cornered tube with the listed H, B, t
    are given (derived from the listed dimensions; they enter no compared value)."""
    H, B, t = i["H"], i["B"], i["t"]
    Ix = (B * H ** 3 - (B - 2.0 * t) * (H - 2.0 * t) ** 3) / 12.0
    Iy = i["Ag"] * i["ry"] ** 2
    assert Ix > Iy
    return DesignSection.rectangular("rectHSS", A=i["Ag"], H=H, B=B, t=t, h_ratio=i["h_t"], b_ratio=i["b_t"],
                                     Ix=Ix, ry=i["ry"], Zx=i["Zx"], J=i["J"])


@pytest.mark.parametrize("method,key", [("LRFD", "phiVn"), ("ASD", "Vn_over_Omega")])
def test_de_h5c_rectangular_hss_shear_g4(method, key):
    """Example H.5C (H-21, H-22; PDF 283-284): h = 6.00 - 3 (0.233), Aw = 2 h t, kv = 5, Cv2 = 1.0
    (Eq. G2-9), Vn = 74.1 kips (Eq. G4-1), 66.7 / 44.4 kips."""
    x = example("H.5C")
    i, p = x["inputs"], x["printed"]
    c = _by_name(gs.check(_ctx(_h5c_section(i), method, i["Lb_ft"] * FT, Fy=i["Fy"]),
                          Demand(method, Vy=i[method]["Vr"])))[gs.NAME_VY]
    assert c.status == "ok" and "G4-1" in c.clause and "G2-9" in c.clause
    assert_printed(c.details["h"][0], p["h"], "H.5C h")
    assert_printed(c.details["Aw = 2 h t"][0], p["Aw"], "H.5C Aw")
    assert c.details["kv"][0] == 5.0
    assert_printed(c.details["1.10 sqrt(kv E/Fy)"][0], p["shear_limit"], "H.5C 1.10 sqrt(kv E/Fy)")
    assert_printed(c.nominal, p["Vn"], "H.5C Vn")
    assert_printed(c.capacity, p[key], f"H.5C {method}")


@pytest.mark.parametrize("method,key", [("LRFD", "phiMn_kipin"), ("ASD", "Mn_over_Omega_kipin")])
def test_de_h5c_rectangular_hss_flexure_f7(method, key):
    """Example H.5C (H-22, H-23; PDF 284-285): Mp = 427 kip-in. (Eq. F7-1), compact flange and web (Cases 17,
    19), Lp = 143 in. (Eq. F7-10; the example's citation reads F7-12) > Lb = 8 ft: no LTB; 384 / 256 kip-in.
    The example lists no Sx: before the review fix the engine answered notChecked here."""
    x = example("H.5C")
    i, p = x["inputs"], x["printed"]
    ctx = _ctx(_h5c_section(i), method, i["Lb_ft"] * FT, Fy=i["Fy"], Lb=i["Lb_ft"] * FT)
    res = ff.axis_strength(ctx, Demand(method, extra={"Mrx": 84.6}), "x")
    assert [c.name for c in res.checks] == [ff.NAME_Y["x"]]              # yielding only: no FLB, WLB, LTB
    y = res.checks[0]
    assert y.status == "ok" and "F7-1" in y.clause
    assert_printed(y.nominal, p["Mp_kipin"], "H.5C Mp")
    assert_printed(y.details["lambda_p flange (width wall) (Table B4.1b Case 17)"][0], p["lambda_p_flange"],
                   "H.5C lambda_p flange")
    assert_printed(y.details["lambda_p web (depth wall) (Table B4.1b Case 19)"][0], p["lambda_p_web"],
                   "H.5C lambda_p web")
    assert_printed(y.details["Lp (Eq. F7-10)"][0], p["Lp_in"], "H.5C Lp")
    assert "Lb <= Lp" in y.details["LTB"][0]
    assert_printed(res.Mc, p[key], f"H.5C {method} Mc")
    assert res.governing == ff.NAME_Y["x"]


# ============================================================================== (b) regression: fixes
E, G, FY, FU = 200000.0, 77200.0, 345.0, 450.0
HC, HT = hc.NAME_H1["c"], hc.NAME_H1["t"]


def _i_section():
    """Rolled I (made-up numbers), nonslender in compression and compact in flexure."""
    return DesignSection.i_shape("rolled", A=9000.0, rx=160.0, ry=40.0, Zx=1.2e6, Sx=1.05e6, Zy=2.4e5, Sy=1.4e5,
                                 J=5.0e5, Cw=7.0e11, rts=48.0, ho=380.0, flange_ratio=7.0, web_ratio=30.0, tf=14.0,
                                 tw=9.0, d=394.0, bf=196.0)


def _round_section():
    return DesignSection.round_hss(A=11000.0, D=406.0, t=8.9, I=2.2e8, J=4.4e8, S=1.08e6, Z=1.40e6)


def _si_ctx(sec, L=6000.0, **kw):
    base = dict(member="M", section=sec, method="LRFD", E=E, G=G, Fy=FY, Fu=FU, L=L, Lx=L, Ly=L, Lz=L, Lb=0.0)
    base.update(kw)
    return DesignContext(**base)


def _uniform_load_stations(n, L, w, N):
    """Simply supported member under a uniform load w along -local y and a constant axial force N:
    Mz(x) = w x (L - x)/2, dMz/dx = -Vy."""
    x = np.linspace(0.0, L, n)
    f = np.zeros((n, 6))
    f[:, 0] = N
    f[:, 5] = w * x * (L - x) / 2.0
    f[:, 1] = -(w * (L / 2.0 - x))
    return x, f


@pytest.mark.parametrize("n", [2, 3, 4, 6, 11])
def test_h1_interaction_does_not_depend_on_the_number_of_stations(n):
    """Beam-column, uniform load, constant compression: the largest moment w L^2/8 is at midspan.  Hand value
    (Eq. H1-1a): Pr/Pc + 8/9 (w L^2/8)/Mcx.  With two stations (member ends, both M = 0) the engine used to
    produce NO interaction check and report "pass"; with an even number of stations it used the largest
    STATION moment."""
    L, w, N = 6000.0, 40.0, -600.0e3
    ctx = _si_ctx(_i_section())
    x, f = _uniform_load_stations(n, L, w, N)
    res = design.check_member(ctx, [Demand("U", Pc=-N, extra={"x": x, "forces": f})])
    c = _by_name(res.checks)
    Pc = min(c[k].capacity for k in (ec.NAME_FB["x"], ec.NAME_FB["y"], ec.NAME_TB))
    Mcx = 0.9 * FY * 1.2e6                                    # Lb = 0: yielding, Eq. F2-1
    Mmax = w * L * L / 8.0                                    # 180e6 N mm
    expected = 600.0e3 / Pc + 8.0 / 9.0 * Mmax / Mcx
    h = c[HC]
    assert 600.0e3 / Pc >= 0.2 and h.details["equation"][0] == "H1-1a"
    assert h.ratio == pytest.approx(expected, rel=1e-12) and expected > 1.0
    assert h.x == pytest.approx(L / 2.0) and h.details["Mrx"][0] == pytest.approx(Mmax, rel=1e-12)
    assert res.status == "fail" and res.passed is False and res.governing.name == HC


def test_round_section_resultant_moment_between_stations():
    """Round section, uniform loads in both planes, two stations (ends): the resultant moment at midspan is
    sqrt(wy^2 + wz^2) L^2/8.  It used to be taken at the stations only (here zero)."""
    L, wy, wz = 6000.0, 30.0, 40.0
    x = np.array([0.0, L])
    f = np.zeros((2, 6))
    f[:, 1] = -(wy * (L / 2.0 - x))                           # dMz/dx = -Vy
    f[:, 2] = wz * (L / 2.0 - x)                              # dMy/dx = +Vz
    mo = ff.moments(Demand("U", extra={"x": x, "forces": f}))
    assert mo.Mr == pytest.approx(50.0 * L * L / 8.0, rel=1e-12) and mo.x_Mr == pytest.approx(L / 2.0)
    res = design.check_member(_si_ctx(_round_section()), [Demand("U", extra={"x": x, "forces": f})])
    y = _by_name(res.checks)[ff.NAME_Y_ROUND]
    assert y.demand == pytest.approx(50.0 * L * L / 8.0, rel=1e-12) and y.ratio > 0.0


def _brute_force(Nf, Mzf, Myf, L, round_, kind, Pc, Mcx, Mcy):
    xx = np.linspace(0.0, L, 100001)
    Nn = Nf(xx)
    mask = (Nn <= 0.0) if kind == "c" else (Nn >= 0.0)
    if not mask.any():
        return None
    mx, my = (np.hypot(Mzf(xx), Myf(xx)), 0.0 * xx) if round_ else (np.abs(Mzf(xx)), np.abs(Myf(xx)))
    axial = np.abs(Nn) / Pc
    flex = mx / Mcx + my / Mcy
    return float(np.where(axial >= 0.2, axial + 8.0 / 9.0 * flex, axial / 2.0 + flex)[mask].max())


@pytest.mark.parametrize("seed", range(40))
def test_h1_interaction_between_stations_against_dense_sampling(seed):
    """Random uniform loads in both planes, end moments and a linearly varying axial force (which may change
    sign), 2 to 4 stations: the engine's interaction is never below the largest value found by evaluating
    Eq. H1-1a / H1-1b with the same strengths at 100 001 points (I-shape: two flexural terms; round section:
    resultant moment)."""
    rng = np.random.default_rng(1000 + seed)
    round_ = bool(seed % 2)
    L = 6000.0
    ctx = _si_ctx(_round_section() if round_ else _i_section())
    n = int(rng.integers(2, 5))
    x = np.linspace(0.0, L, n)
    wy, wz = rng.uniform(-60.0, 60.0), rng.uniform(-25.0, 25.0)
    m0, m1 = rng.uniform(-150.0e6, 150.0e6, 2)
    q0, q1 = rng.uniform(-45.0e6, 45.0e6, 2)
    n0 = rng.uniform(-900.0e3, 300.0e3)
    n1 = n0 if seed % 4 >= 2 else rng.uniform(-900.0e3, 300.0e3)

    def Mzf(s):
        return m0 + (m1 - m0) * s / L + wy * s * (L - s) / 2.0

    def Myf(s):
        return q0 + (q1 - q0) * s / L + wz * s * (L - s) / 2.0

    def Nf(s):
        return n0 + (n1 - n0) * s / L

    f = np.zeros((n, 6))
    f[:, 0], f[:, 5], f[:, 4] = Nf(x), Mzf(x), Myf(x)
    f[:, 1] = -((m1 - m0) / L + wy * (L / 2.0 - x))           # dMz/dx = -Vy
    f[:, 2] = (q1 - q0) / L + wz * (L / 2.0 - x)              # dMy/dx = +Vz
    demand = Demand("U", Pc=max(0.0, -min(n0, n1)), Pt=max(0.0, max(n0, n1)), extra={"x": x, "forces": f})
    checks = hc.interaction(ctx, demand)
    assert checks and all(c.status == "ok" for c in checks)
    for c in checks:
        kind = "c" if c.name == HC else "t"
        Pc = c.inputs["Pc"][0] or math.inf
        Mcx = c.inputs["Mcx"][0] or math.inf
        Mcy = (c.inputs["Mcy"][0] if "Mcy" in c.inputs else None) or math.inf
        brute = _brute_force(Nf, Mzf, Myf, L, round_, kind, Pc, Mcx, Mcy)
        assert brute is not None
        assert c.ratio >= brute * (1.0 - 1.0e-5), (c.name, c.ratio, brute)
        # and the engine's value is a value of the interaction at a point of the member (not an over-estimate):
        s = c.x
        mx = math.hypot(Mzf(s), Myf(s)) if round_ else abs(Mzf(s))
        my = 0.0 if round_ else abs(Myf(s))
        assert c.details["Mrx"][0] == pytest.approx(mx, rel=1e-9, abs=1e-3)
        assert c.details["Mry"][0] == pytest.approx(my, rel=1e-9, abs=1e-3)
        assert c.details["Pr"][0] == pytest.approx(abs(Nf(s)), rel=1e-9, abs=1e-3)


def _tube(**kw):
    """Rectangular HSS (made-up numbers), compact in flexure, stronger about local z."""
    base = dict(A=5000.0, H=250.0, B=150.0, t=8.0, Ix=4.2e7, Iy=1.9e7, J=4.4e7, Zx=4.1e5)
    base.update(kw)
    return DesignSection.rectangular("rectHSS", **base)


def test_rectangular_ltb_needs_the_elastic_modulus_only_beyond_lp():
    """F7.4(a): no LTB when Lb <= Lp, and Lp (Eq. F7-10) = 0.13 E ry sqrt(J Ag)/Mp needs no elastic section
    modulus.  A section without Sx is therefore checked when Lb <= Lp and notChecked (with the reason) only
    when a segment is longer than Lp."""
    sec = _tube()
    Mp = FY * 4.1e5
    ry = math.sqrt(1.9e7 / 5000.0)
    # Eq. F7-10: ry = sqrt(3800) = 61.644; sqrt(J Ag) = 469 042; Mp = 141.45e6 -> Lp = 5314.6 mm
    Lp = 0.13 * E * ry * math.sqrt(4.4e7 * 5000.0) / Mp
    assert 5314.0 < Lp < 5315.0
    demand = Demand("U", extra={"Mrx": 50.0e6})
    short = ff.axis_strength(_si_ctx(sec, L=5000.0, Lb=5000.0), demand, "x")
    assert [c.status for c in short.checks] == ["ok"] and short.Mc == pytest.approx(0.9 * Mp, rel=1e-13)
    assert short.checks[0].details["Lp (Eq. F7-10)"][0] == pytest.approx(Lp, rel=1e-13)
    long_ = ff.axis_strength(_si_ctx(sec, L=14000.0, Lb=14000.0), demand, "x")
    ltb = _by_name(long_.checks)[ff.NAME_LTB]
    assert ltb.status == "notChecked" and "Sx" in ltb.reason and "F7-8" in ltb.reason and long_.Mc is None
    # with Sx the same member is checked by Eq. F7-8
    full = ff.axis_strength(_si_ctx(_tube(Sx=3.4e5), L=14000.0, Lb=14000.0), demand, "x")
    ltb = _by_name(full.checks)[ff.NAME_LTB]
    Lr = 2.0 * E * ry * math.sqrt(4.4e7 * 5000.0) / (0.7 * FY * 3.4e5)      # Eq. F7-11
    Mn = min(Mp, 1.0 * (Mp - (Mp - 0.7 * FY * 3.4e5) * (14000.0 - Lp) / (Lr - Lp)))
    assert ltb.status == "ok" and "F7-8" in ltb.clause and ltb.nominal == pytest.approx(Mn, rel=1e-13)


def test_deflection_limit_with_end_stations_only_is_not_checked():
    """With analysis.resultStations = 2 a design member has only its two end stations; the deflection relative
    to the chord through the ends is then zero by definition.  This used to be reported as a satisfied limit
    (ratio 0); it is notChecked, and the limit is checked as soon as there is a station in between."""
    from benchfea.design.aisc360_22 import serviceability as sv
    from benchfea.model import build_model
    from benchfea.results import CombinationResult, MemberResult, Results

    span = 8000.0
    data = {
        "formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": span, "y": 0, "z": 0}},
        "supports": {"A": {"type": "pinned", "dofs": {"rx": "fixed"}}, "B": {"type": "roller", "dofs": {"uy": "fixed"}}},
        "materials": {"S": {"E": 200000.0, "nu": 0.3, "Fy": 345, "Fu": 450}},
        "sections": {"I400": {"type": "I", "d": 400, "bf": 180, "tf": 14, "tw": 9, "r": 20}},
        "members": {"beam": {"i": "A", "j": "B", "section": "I400", "material": "S"}},
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fx": 1}]}},
        "combinations": {"U1": {"factors": {"D": 1.4}}, "S1": {"factors": {"D": 1.0}, "purpose": "serviceability"}},
        "design": {"members": {"beam": {"deflectionLimits": [{"direction": "z", "spanRatio": 360}]}}},
    }
    model = build_model(data)

    def results(x, uz):
        x = np.asarray(x, float)
        u = np.column_stack([np.zeros(len(x)), np.zeros(len(x)), np.asarray(uz, float)])
        mr = {"beam": MemberResult(np.zeros(12), x, np.zeros((len(x), 6)), u)}
        n = len(model.nodes)
        return Results("hand", {"S1": CombinationResult("S1", "firstOrder", np.zeros((n, 6)), np.zeros((n, 6)), mr)})

    two = sv.deflection_checks(model, results([0.0, span], [0.0, 0.0]), "beam")[0]
    assert two.status == "notChecked" and two.passed is None and "no result station between its ends" in two.reason
    three = sv.deflection_checks(model, results([0.0, span / 2.0, span], [0.0, -40.0, 0.0]), "beam")[0]
    assert three.status == "ok" and three.demand == pytest.approx(40.0) and three.passed is False     # 40 > 8000/360
