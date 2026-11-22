"""AISC 360-22 Chapter H: H1.1 / H1.2 interaction (Eq. H1-1a, H1-1b), the torsion
screening rule of the engine (H3 is not implemented), and the formulas of
Appendix 8, Section 8.1.2 (Eq. A-8-3 to A-8-5).

Tier F: (a) hand calculations from the Specification (pp. 16.1-82, 16.1-83,
16.1-284, 16.1-285) as explicit arithmetic; (b) AISC Design Examples v16.0 H.1B,
H.3, H.4 (tests/design_examples/part2; the ONE tolerance policy).
"""
import math

import numpy as np
import pytest

from benchfea import design
from benchfea.design import DesignContext, DesignSection, Demand
from benchfea.design.aisc360_22 import d_tension, e_compression, f_flexure as ff, h_combined as hc
from design_examples.part2 import FT, assert_printed, example

E, G, FY, FU = 200000.0, 77200.0, 345.0, 450.0
PI = math.pi
HC, HT = hc.NAME_H1["c"], hc.NAME_H1["t"]


def _sec(**kw):
    """Rolled I (made-up numbers): nonslender in compression (flange 7 <= 13.48, web 30 <= 35.87, Table B4.1a)
    and compact in flexure; A = 9000, rx = 160, ry = 40 -> Ix = 2.304e8, Iy = 1.44e7."""
    base = dict(A=9000.0, rx=160.0, ry=40.0, Zx=1.2e6, Sx=1.05e6, Zy=2.4e5, Sy=1.4e5, J=5.0e5, Cw=7.0e11, rts=48.0,
                ho=380.0, flange_ratio=7.0, web_ratio=30.0, tf=14.0, tw=9.0)
    base.update(kw)
    return DesignSection.i_shape("rolled", **base)


def _ctx(sec=None, method="LRFD", **kw):
    base = dict(member="BC", section=sec if sec is not None else _sec(), method=method, E=E, G=G, Fy=FY, Fu=FU,
                L=6000.0, Lx=6000.0, Ly=6000.0, Lz=6000.0, Lb=0.0)
    base.update(kw)
    return DesignContext(**base)


def _d(Pc=0.0, Pt=0.0, Mrx=0.0, Mry=0.0, Tr=0.0, **extra):
    return Demand("U1", Pc=Pc, Pt=Pt, extra={"Mrx": Mrx, "Mry": Mry, "Tr": Tr, **extra})


def _stations(x, n, mz=None, my=None, t=None, analysis=None):
    x = np.asarray(x, float)
    f = np.zeros((len(x), 6))
    for col, v in ((0, n), (3, t), (4, my), (5, mz)):
        if v is not None:
            f[:, col] = v
    d = design.demand_from_stations(x, f, "U1")
    extra = {"x": x, "forces": f}
    if analysis:
        extra["analysisType"] = analysis
    return Demand("U1", Pt=d.Pt, x_Pt=d.x_Pt, Pc=d.Pc, x_Pc=d.x_Pc, extra=extra)


def _by_name(checks):
    return {c.name: c for c in checks}


# hand values of the available strengths of _sec() / _ctx() (LRFD)
def _hand_strengths():
    # E3 minor axis: Lc/r = 6000/40 = 150 > 4.71 sqrt(E/Fy) = 113.4 -> Fe = pi^2 E/150^2 = 87.73 (E3-4),
    # Fn = 0.877 Fe = 76.94 (E3-3), Pn = Fn Ag = 692.4e3 (E3-1); phi_c = 0.90
    Fe = PI ** 2 * 200000.0 / 150.0 ** 2
    Pc = 0.90 * 0.877 * Fe * 9000.0
    assert 87.7 < Fe < 87.8 and 623.1e3 < Pc < 623.3e3
    # D2: yielding 0.90 * 345 * 9000 = 2794.5e3 (D2-1); rupture 0.75 * 450 * 9000 = 3037.5e3 (D2-2) -> 2794.5e3
    Pt = min(0.90 * 345.0 * 9000.0, 0.75 * 450.0 * 9000.0)
    # F2-1 (Lb = 0): 0.90 * 345 * 1.2e6 = 372.6e6; F6-1: 0.90 * min(345 * 2.4e5, 1.6 * 345 * 1.4e5) = 69.552e6
    Mcx = 0.90 * 345.0 * 1.2e6
    Mcy = 0.90 * min(345.0 * 2.4e5, 1.6 * 345.0 * 1.4e5)
    assert (Pt, Mcx) == (2794.5e3, 372.6e6) and Mcy == pytest.approx(69.552e6)
    return Pc, Pt, Mcx, Mcy


# ============================================================================== tier F (a): H1
def test_hand_h1_interaction_equations():
    # Eq. H1-1a (Pr/Pc >= 0.2): 0.5 + 8/9 (0.3 + 0.1) = 0.85556
    v, eq = hc.h1_interaction(50.0, 100.0, 30.0, 100.0, 10.0, 100.0)
    assert eq == "H1-1a" and v == pytest.approx(0.5 + 8.0 / 9.0 * 0.4, rel=1e-15)
    # Eq. H1-1b (Pr/Pc < 0.2): 0.1/2 + (0.3 + 0.1) = 0.45
    v, eq = hc.h1_interaction(10.0, 100.0, 30.0, 100.0, 10.0, 100.0)
    assert eq == "H1-1b" and v == pytest.approx(0.05 + 0.4, rel=1e-15)
    # exactly 0.2 belongs to H1-1a ("when Pr/Pc >= 0.2"): 0.2 + 8/9 * 0.4 = 0.5556
    v, eq = hc.h1_interaction(20.0, 100.0, 40.0, 100.0)
    assert eq == "H1-1a" and v == pytest.approx(0.2 + 8.0 / 9.0 * 0.4)
    # the two equations meet on the limit: at Pr/Pc = 0.2 and Mr/Mc = 0.9 both give exactly 1.0
    assert hc.h1_interaction(20.0, 100.0, 90.0, 100.0)[0] == pytest.approx(1.0, rel=1e-15) == 0.2 / 2.0 + 0.9
    # all terms positive (User Note)
    assert hc.h1_interaction(-50.0, 100.0, -30.0, 100.0, -10.0, 100.0)[0] == pytest.approx(0.5 + 8.0 / 9.0 * 0.4)


def test_hand_h1_1_flexure_and_compression():
    Pc, Pt, Mcx, Mcy = _hand_strengths()
    # Pr = 200e3: Pr/Pc = 0.3209 >= 0.2 -> Eq. H1-1a with Mrx = 150e6, Mry = 10e6:
    # 0.3209 + 8/9 (0.40258 + 0.14378) = 0.8066
    value = 200.0e3 / Pc + 8.0 / 9.0 * (150.0e6 / Mcx + 10.0e6 / Mcy)
    assert 0.806 < value < 0.807
    c, = hc.interaction(_ctx(), _d(Pc=200.0e3, Mrx=150.0e6, Mry=10.0e6))
    assert c.name == HC and c.status == "ok" and c.ratio == pytest.approx(value, rel=1e-13) and c.passed
    assert c.clause == "AISC 360-22 H1.1, Eq. H1-1a; Chapter E; Chapter F"
    assert (c.capacity, c.dimension, c.demand) == (1.0, "dimensionless", c.ratio)
    assert c.inputs["Pc"][0] == pytest.approx(Pc, rel=1e-13) and c.inputs["Pc"][2] == "derived: " + e_compression.NAME_FB["y"]
    assert c.inputs["Mcx"] == (pytest.approx(Mcx), "moment", "derived: " + ff.NAME_Y["x"])
    assert c.inputs["Mcy"] == (pytest.approx(Mcy), "moment", "derived: " + ff.NAME_Y["y"])
    assert c.details["Pr/Pc"][0] == pytest.approx(200.0e3 / Pc) and c.details["equation"][0] == "H1-1a"
    # Pr = 50e3: Pr/Pc = 0.0802 < 0.2 -> Eq. H1-1b: 0.0401 + 0.5464 = 0.5865
    value = 50.0e3 / (2.0 * Pc) + (150.0e6 / Mcx + 10.0e6 / Mcy)
    c, = hc.interaction(_ctx(), _d(Pc=50.0e3, Mrx=150.0e6, Mry=10.0e6))
    assert 0.586 < value < 0.587 and c.ratio == pytest.approx(value, rel=1e-13) and "Eq. H1-1b" in c.clause
    # ASD: the available strengths carry Omega (1.67)
    c, = hc.interaction(_ctx(method="ASD"), _d(Pc=100.0e3, Mrx=100.0e6))
    PcA, McxA = 0.877 * PI ** 2 * E / 150.0 ** 2 * 9000.0 / 1.67, 345.0 * 1.2e6 / 1.67
    assert c.ratio == pytest.approx(100.0e3 / PcA + 8.0 / 9.0 * (100.0e6 / McxA), rel=1e-13) and c.method == "ASD"
    # Mcx is the LOWEST Chapter F strength: with Lb = 6000 > Lr = 5276 the elastic LTB strength governs
    lam2 = (6000.0 / 48.0) ** 2
    Fcr = PI ** 2 * E / lam2 * math.sqrt(1.0 + 0.078 * 5.0e5 / (1.05e6 * 380.0) * lam2)       # Eq. F2-4, Cb = 1.0
    c, = hc.interaction(_ctx(Lb=6000.0), _d(Pc=200.0e3, Mrx=100.0e6))
    assert c.inputs["Mcx"][2] == "derived: " + ff.NAME_LTB
    assert c.ratio == pytest.approx(200.0e3 / Pc + 8.0 / 9.0 * 100.0e6 / (0.9 * Fcr * 1.05e6), rel=1e-12)


def test_hand_h1_2_flexure_and_tension():
    Pc, Pt, Mcx, Mcy = _hand_strengths()
    # Pr = 1000e3 tension: Pr/Pc = 0.3578 -> H1-1a: 0.3578 + 8/9 (100e6/372.6e6) = 0.5964
    value = 1000.0e3 / Pt + 8.0 / 9.0 * (100.0e6 / Mcx)
    c, = hc.interaction(_ctx(), _d(Pt=1000.0e3, Mrx=100.0e6))
    assert c.name == HT and c.ratio == pytest.approx(value, rel=1e-13) and 0.596 < value < 0.597
    assert c.clause == "AISC 360-22 H1.2, Eq. H1-1a; Chapter D; Chapter F"
    assert c.inputs["Pc"] == (pytest.approx(Pt), "force", "derived: " + d_tension.NAME_YIELD)
    # rupture governs Pc when An U is small: An = 5000 -> 0.75 * 450 * 5000 = 1687.5e3
    c, = hc.interaction(_ctx(An=5000.0), _d(Pt=1000.0e3, Mrx=100.0e6))
    assert c.inputs["Pc"] == (pytest.approx(1687.5e3), "force", "derived: " + d_tension.NAME_RUPTURE)
    # the optional Cb increase of H1.2 (Eq. H1-2) is not applied: tension does not raise Mcx
    a, = hc.interaction(_ctx(Lb=6000.0), _d(Pt=1000.0e3, Mrx=100.0e6))
    b, = hc.interaction(_ctx(Lb=6000.0), _d(Pt=1.0, Mrx=100.0e6))
    assert a.inputs["Mcx"][0] == b.inputs["Mcx"][0]


def test_interaction_at_every_station_with_concurrent_demands():
    Pc, Pt, Mcx, Mcy = _hand_strengths()
    # column 6 m: compression 300 kN over the lower half, 100 kN above; moment linear 0 ... 120e6 ... 0 (peak at
    # mid-height, where the two axial values meet at a double station)
    x = [0.0, 1500.0, 3000.0, 3000.0, 4500.0, 6000.0]
    n = [-300.0e3, -300.0e3, -300.0e3, -100.0e3, -100.0e3, -100.0e3]
    mz = [0.0, 60.0e6, 120.0e6, 120.0e6, 60.0e6, 0.0]
    c, = hc.interaction(_ctx(), _stations(x, n, mz=mz))
    # governing: x = 3000 with Pr = 300e3 (H1-1a): 0.4814 + 8/9 * 0.32206 = 0.7676
    value = 300.0e3 / Pc + 8.0 / 9.0 * 120.0e6 / Mcx
    assert c.ratio == pytest.approx(value, rel=1e-13) and c.x == 3000.0 and c.details["Pr"][0] == 300.0e3
    assert any("every station" in s for s in c.notes)
    # the maxima need not be concurrent: maximum compression where the moment is zero
    mz2 = [0.0, 0.0, 0.0, 0.0, 60.0e6, 120.0e6]
    # x = 6000: Pr = 100e3 with M = 120e6 -> H1-1b: 0.0802 + 0.3221 = 0.4023;
    # x = 0 ... 3000: Pr = 300e3 without moment -> H1-1a: 0.4814 (the axial ratio itself) -> governs
    c, = hc.interaction(_ctx(), _stations(x, n, mz=mz2))
    assert 100.0e3 / (2.0 * Pc) + 120.0e6 / Mcx < 300.0e3 / Pc
    assert c.x == 0.0 and c.ratio == pytest.approx(300.0e3 / Pc, rel=1e-13)
    assert c.ratio < 300.0e3 / Pc + 8.0 / 9.0 * 120.0e6 / Mcx              # (what non-concurrent maxima would give)
    # axial force changing sign along the member: compression stations -> H1.1, tension stations -> H1.2
    n3 = [-200.0e3, -200.0e3, -200.0e3, 400.0e3, 400.0e3, 400.0e3]
    checks = _by_name(hc.interaction(_ctx(), _stations(x, n3, mz=mz)))
    assert set(checks) == {HC, HT}
    assert checks[HC].ratio == pytest.approx(200.0e3 / Pc + 8.0 / 9.0 * 120.0e6 / Mcx, rel=1e-13)
    assert checks[HT].ratio == pytest.approx(400.0e3 / (2.0 * Pt) + 120.0e6 / Mcx, rel=1e-13)
    # biaxial flexure without axial force: Eq. H1-1b with Pr = 0
    c, = hc.interaction(_ctx(), _d(Mrx=100.0e6, Mry=20.0e6))
    assert c.ratio == pytest.approx(100.0e6 / Mcx + 20.0e6 / Mcy, rel=1e-13) and c.inputs["Pc"][0] is None
    # nothing to combine: axial force only, or one moment only
    assert hc.interaction(_ctx(), _d(Pc=1.0e5)) == [] and hc.interaction(_ctx(), _d(Mrx=1.0e6)) == []
    assert hc.interaction(_ctx(), Demand("U1", Pc=1.0e5)) == []            # no moment data at all
    # first-order results: said in the check
    c, = hc.interaction(_ctx(), _stations(x, n, mz=mz, analysis="firstOrder"))
    assert any("FIRST-ORDER" in s and "C1" in s for s in c.notes)
    c, = hc.interaction(_ctx(), _stations(x, n, mz=mz, analysis="secondOrder"))
    assert not any("first-order" in s.lower() for s in c.notes)


def test_round_section_uses_the_resultant_moment_once():
    sec = DesignSection.round_hss(A=5000.0, r=70.0, ratio=30.0, Z=3.0e5, S=2.3e5, J=4.9e7)
    ctx = _ctx(sec)
    Mc = 0.9 * 345.0 * 3.0e5
    Pc = e_compression.flexural_buckling_major(ctx, Demand("U1", Pc=1.0))[0].capacity
    c, = hc.interaction(ctx, _d(Pc=0.5 * Pc, Mrx=3.0e7, Mry=4.0e7))
    assert c.ratio == pytest.approx(0.5 + 8.0 / 9.0 * 5.0e7 / Mc, rel=1e-12) and "Mcy" not in c.inputs
    assert any("resultant" in s for s in c.notes)
    # one moment component only: still an interaction with the axial force
    c, = hc.interaction(ctx, _d(Pc=0.5 * Pc, Mry=4.0e7))
    assert c.ratio == pytest.approx(0.5 + 8.0 / 9.0 * 4.0e7 / Mc, rel=1e-12)


def test_interaction_is_not_checked_when_a_strength_is_missing():
    # torsional buckling cannot be evaluated without Cw -> Pc is not available -> interaction notChecked
    c, = hc.interaction(_ctx(_sec(Cw=None)), _d(Pc=1.0e5, Mrx=1.0e8))
    assert c.status == "notChecked" and "torsional buckling" in c.reason and c.demand is None
    # flexure not available (noncompact web -> F4 not implemented)
    c, = hc.interaction(_ctx(_sec(web_ratio=100.0)), _d(Pt=1.0e5, Mrx=1.0e8))
    assert c.status == "notChecked" and "F4" in c.reason
    # a missing strength of an axis without moment does not block
    c, = hc.interaction(_ctx(_sec(Zy=None)), _d(Pt=1.0e5, Mrx=1.0e8))
    assert c.status == "ok" and c.inputs["Mcy"][2].startswith("not needed")
    # the member can then not be reported as passing
    res = design.check_member(_ctx(_sec(Cw=None)), [_d(Pc=1.0e5, Mrx=1.0e8)])
    assert res.passed is False and res.status == "incomplete" and any(HC in s for s in res.not_checked)


def test_interaction_through_check_member_governs_and_is_enveloped():
    Pc, Pt, Mcx, Mcy = _hand_strengths()
    d1 = Demand("U1", Pc=200.0e3, extra={"Mrx": 150.0e6, "Mry": 10.0e6})
    d2 = Demand("U2", Pc=400.0e3, extra={"Mrx": 20.0e6, "Mry": 0.0})
    res = design.check_member(_ctx(), [d1, d2])
    c = _by_name(res.checks)[HC]
    v1 = 200.0e3 / Pc + 8.0 / 9.0 * (150.0e6 / Mcx + 10.0e6 / Mcy)
    v2 = 400.0e3 / Pc + 8.0 / 9.0 * (20.0e6 / Mcx)
    assert v1 > v2 and (c.combination, c.ratio) == ("U1", pytest.approx(v1, rel=1e-13))
    assert res.governing is c and res.passed and res.status == "pass"
    # every check of every combination is kept
    assert set(res.by_combination) == {"U1", "U2"}
    assert _by_name(res.by_combination["U2"])[HC].ratio == pytest.approx(v2, rel=1e-13)
    # an interaction above 1 fails the member although every single strength check passes
    res = design.check_member(_ctx(), [Demand("U1", Pc=400.0e3, extra={"Mrx": 250.0e6, "Mry": 0.0})])
    assert all(c.passed for c in res.checks if c.name != HC and not c.advisory and c.status == "ok")
    assert _by_name(res.checks)[HC].passed is False and res.status == "fail" and res.governing.name == HC


# ============================================================================== torsion: engine screening rule
def test_torsional_shear_stress_formulas_by_numerical_integration():
    """The two formulas of the screening rule, checked against direct integration (module docstring of
    h_combined): open thin-walled tau = T t/J; closed single cell tau = T/(2 A0 t)."""
    # narrow rectangle b x t, twist rate th: phi = G th (t^2/4 - n^2); T = 2 * integral(phi) dA
    b, t, Gth = 200.0, 12.0, 3.0
    n = np.linspace(-t / 2.0, t / 2.0, 20001)
    phi = Gth * (t ** 2 / 4.0 - n ** 2)
    T = 2.0 * b * float(np.sum((phi[1:] + phi[:-1]) / 2.0 * np.diff(n)))
    tau_max = Gth * t                                  # |d phi/dn| = 2 G th |n| at n = t/2
    Jrect = b * t ** 3 / 3.0
    assert T == pytest.approx(Gth * Jrect, rel=1e-8) and tau_max == pytest.approx(T * t / Jrect, rel=1e-8)
    # I-shape of three rectangles with a common twist rate: tau in the thickest element = T t_max / J
    sec = _sec(J=(2 * 180.0 * 14.0 ** 3 + 360.0 * 9.0 ** 3) / 3.0)
    tau, how = hc.torsional_shear_stress(sec, 1.0e6)
    assert tau == pytest.approx(1.0e6 * 14.0 / sec.J, rel=1e-15) and "max(tf, tw)" in how
    # closed cell: constant shear flow q along the wall mid-line; T = integral q (x dy - y dx) = 2 A0 q
    H, B, tw, q = 300.0, 200.0, 8.0, 5.0
    a, c = (B - tw) / 2.0, (H - tw) / 2.0
    corners = np.array([[a, -c], [a, c], [-a, c], [-a, -c], [a, -c]])
    pts = np.vstack([np.linspace(corners[k], corners[k + 1], 400, endpoint=False) for k in range(4)] + [corners[:1]])
    T = q * float(np.sum(pts[:-1, 0] * np.diff(pts[:, 1]) - pts[:-1, 1] * np.diff(pts[:, 0])))
    assert T == pytest.approx(2.0 * (H - tw) * (B - tw) * q, rel=1e-12)
    box = DesignSection.rectangular("rectHSS", A=7744.0, H=H, B=B, t=tw)
    assert hc.torsional_shear_stress(box, T)[0] == pytest.approx(q / tw, rel=1e-12)       # tau = q / t
    # round tube: A0 = pi (D - t)^2/4
    pipe = DesignSection.round_hss(A=5000.0, D=200.0, t=8.0)
    assert hc.torsional_shear_stress(pipe, 1.0e7)[0] == pytest.approx(1.0e7 / (2 * PI * 192.0 ** 2 / 4 * 8.0), rel=1e-15)
    # no dimensions: cannot be screened
    assert hc.torsional_shear_stress(DesignSection.round_hss(A=5000.0, ratio=25.0), 1.0)[0] is None


def test_torsion_is_not_checked_unless_negligible_by_the_engine_rule():
    sec = _sec()
    limit = 0.01 * 0.6 * 345.0                         # engine rule: 1 % of 0.6 Fy = 2.07 MPa
    T_limit = limit * sec.J / 14.0                     # tau = T tf / J
    # negligible: advisory record, nothing blocks
    c, = hc.torsion(_ctx(), _d(Tr=0.5 * T_limit))
    assert c.name == hc.NAME_TORSION_SCREEN and c.advisory and c.status == "ok" and c.passed
    assert c.ratio == pytest.approx(0.5) and c.capacity == pytest.approx(limit) and c.dimension == "force/length^2"
    assert "engine" in c.clause and any("not an AISC rule" in s for s in c.notes) and c.combination == "U1"
    res = design.check_member(_ctx(), [_d(Pt=1.0e5, Mrx=1.0e7, Tr=0.5 * T_limit)])
    assert res.passed and res.status == "pass"
    # not negligible: H3 is not implemented -> notChecked with the torsional moment as demand
    c, = hc.torsion(_ctx(), _d(Tr=2.0 * T_limit))
    assert c.name == hc.NAME_TORSION and c.status == "notChecked" and c.demand == pytest.approx(2.0 * T_limit)
    assert "H3" in c.reason and "engine rule" in c.reason and c.dimension == "moment" and not c.advisory
    res = design.check_member(_ctx(), [_d(Pt=1.0e5, Mrx=1.0e7, Tr=2.0 * T_limit)])
    assert res.passed is False and res.status == "incomplete" and any("torsion" in s for s in res.not_checked)
    # no torsional moment: no torsion record at all
    assert hc.torsion(_ctx(), _d(Mrx=1.0e7)) == [] and hc.torsion(_ctx(), Demand("U1", Pc=1.0)) == []
    # torque that cannot be screened (no dimensions, or no Fy) is never neglected
    assert hc.torsion(_ctx(_sec(tf=None, tw=None)), _d(Tr=1.0))[0].status == "notChecked"
    assert hc.torsion(_ctx(Fy=None), _d(Tr=1.0))[0].status == "notChecked"
    # from stations: the largest |T| and where it occurs
    c, = hc.torsion(_ctx(), _stations([0.0, 1000.0, 2000.0], [0.0, 0.0, 0.0], t=[0.0, -3.0 * T_limit, T_limit]))
    assert c.demand == pytest.approx(3.0 * T_limit) and c.x == 1000.0


# ============================================================================== Appendix 8, Section 8.1.2
def test_hand_appendix_8_b1_formulas():
    # Eq. A-8-4: Cm = 0.6 - 0.4 (M1/M2), M1/M2 positive for reverse curvature, negative for single curvature.
    # With signed station moments (same sign convention at both ends) equal signs = single curvature.
    assert hc.cm_a_8_4(100.0, 100.0) == pytest.approx(0.6 + 0.4)            # uniform moment: M1/M2 = -1 -> 1.0
    assert hc.cm_a_8_4(100.0, -100.0) == pytest.approx(0.6 - 0.4)           # full reverse curvature: +1 -> 0.2
    assert hc.cm_a_8_4(0.0, 80.0) == pytest.approx(0.6)                     # one end moment zero
    assert hc.cm_a_8_4(-50.0, 100.0) == pytest.approx(0.6 - 0.4 * 0.5)      # reverse, M1/M2 = +0.5 -> 0.4
    assert hc.cm_a_8_4(100.0, -50.0) == pytest.approx(0.4)                  # order of the ends does not matter
    with pytest.raises(ValueError):
        hc.cm_a_8_4(0.0, 0.0)
    # Eq. A-8-5: Pe1 = pi^2 EI*/Lc1^2: E = 200000, I = 1e8, Lc1 = 5000 -> 7.8957e6 N
    Pe1 = hc.pe1_a_8_5(200000.0 * 1.0e8, 5000.0)
    assert Pe1 == pytest.approx(PI ** 2 * 2.0e13 / 25.0e6) and 7.89e6 < Pe1 < 7.90e6
    # Eq. A-8-3: B1 = Cm/(1 - alpha Pr/Pe1) >= 1: Cm = 1, Pr = Pe1/5 -> 1.25; ASD alpha = 1.6 -> 1/(1 - 0.32)
    assert hc.b1_a_8_3(1.0, 1.0, Pe1 / 5.0, Pe1) == pytest.approx(1.25)
    assert hc.b1_a_8_3(1.0, 1.6, Pe1 / 5.0, Pe1) == pytest.approx(1.0 / 0.68)
    assert hc.b1_a_8_3(0.4, 1.0, Pe1 / 5.0, Pe1) == 1.0                     # 0.5 -> not less than 1
    assert hc.b1_a_8_3(1.0, 1.0, Pe1, Pe1) == math.inf                      # alpha Pr >= Pe1: not covered
    # the engine does NOT amplify: a first-order demand is checked as given
    Pc, Pt, Mcx, Mcy = _hand_strengths()
    c, = hc.interaction(_ctx(), _d(Pc=200.0e3, Mrx=150.0e6, analysisType="firstOrder"))
    assert c.ratio == pytest.approx(200.0e3 / Pc + 8.0 / 9.0 * 150.0e6 / Mcx, rel=1e-13)


# ============================================================================== tier F (b): Design Examples v16.0
@pytest.mark.parametrize("method", ["LRFD", "ASD"])
def test_design_example_h1b_interaction_with_tabulated_strengths(method):
    """Example H.1B (H-4 to H-5, PDF 266-267): W14x99, Eq. H1-1a with the available strengths the example
    takes from Manual Table 6-1 (no section property is listed, so only Chapter H is reproduced)."""
    ex = example("H.1B")
    i, p = ex["inputs"][method], ex["printed"][method]
    assert_printed(i["Pr"] / i["Pc"], p["Pr_over_Pc"], "Pr/Pc")
    value, eq = hc.h1_interaction(i["Pr"], i["Pc"], i["Mrx_kipft"], i["Mcx_kipft"], i["Mry_kipft"], i["Mcy_kipft"])
    assert eq == "H1-1a"
    assert_printed(value, p["interaction"], "Eq. H1-1a")


@pytest.mark.parametrize("method", ["LRFD", "ASD"])
def test_design_example_h3_tension_and_biaxial_flexure(method):
    """Example H.3 (H-8 to H-11, PDF 270-273): W14x82, tension + biaxial flexure through the engine's checks
    (D2-1, F2-2 with the listed Lp / Lr, F6-1, H1-1b).  Cb = 1.41 is the example's value (it includes the
    optional H1.2 increase, which the engine does not apply by itself) and is given as the design variable."""
    ex = example("H.3")
    i, p = ex["inputs"], ex["printed"]
    r, q = i[method], p[method]
    # --- the example's Eq. H1-2 arithmetic (independent of the engine): Pey and the Cb factors
    Pey = PI ** 2 * i["E"] * i["Iy"] / (i["Lb_ft"] * FT) ** 2
    assert_printed(Pey, p["Pey"], "Pey (Eq. H1-2)")
    assert_printed(math.sqrt(1.0 + 1.0 * i["LRFD"]["Pr"] / Pey), p["Cb_factor_LRFD"], "Cb factor LRFD")
    assert_printed(math.sqrt(1.0 + 1.6 * i["ASD"]["Pr"] / Pey), p["Cb_factor_ASD"], "Cb factor ASD")
    assert_printed(p["Cb_factor_LRFD"] * i["Cb_table"], i["Cb_used"], "Cb used")
    # --- strengths
    factor = 0.90 if method == "LRFD" else 1.0 / 1.67
    Pn = i["Fy"] * i["Ag"]                                                  # Eq. D2-1
    assert_printed(Pn, p["Pn"], "Pn")
    Mp = i["Fy"] * i["Zx"]
    assert_printed(Mp, p["Mp_kipin"], "Mp")
    Mnx = ff.inelastic_ltb(i["Cb_used"], Mp, i["Fy"], i["Sx"], i["Lb_ft"], i["Lp_ft"], i["Lr_ft"])      # Eq. F2-2
    assert_printed(Mnx, p["Mnx_kipin"], "Mnx")
    assert_printed(Mnx / FT, p["Mnx_kipft"], "Mnx")
    Mny = ff.mp_minor(i["Fy"], i["Zy"], i["Sy"])                            # Eq. F6-1
    assert_printed(Mny, p["Mny_kipin"], "Mny")
    assert_printed(Mny / FT, p["Mny_kipft"], "Mny")
    assert_printed(factor * Pn, q["Pc"], "Pc")
    assert_printed(factor * Mnx / FT, q["Mcx_kipft"], "Mcx")
    assert_printed(factor * Mny / FT, q["Mcy_kipft"], "Mcy")
    assert_printed(r["Pr"] / (factor * Pn), q["Pr_over_Pc"], "Pr/Pc")
    value, eq = hc.h1_interaction(r["Pr"], factor * Pn, r["Mrx_kipft"] * FT, factor * Mnx, r["Mry_kipft"] * FT,
                                  factor * Mny)
    assert eq == "H1-1b"
    assert_printed(value, q["interaction"], "Eq. H1-1b")
    # --- the same through the limit-state functions of the engine: the D2 and F6 checks and the H1.2 check
    # (Lb = 0 in the context so that Mcx is supplied by Eq. F2-2 above through a section whose Zx gives Mnx)
    root = math.sqrt(i["E"] / i["Fy"])
    sec = DesignSection.i_shape("rolled", A=i["Ag"], Zx=Mnx / i["Fy"], Zy=i["Zy"], Sy=i["Sy"],
                                flange_ratio=0.5 * 0.38 * root, web_ratio=0.5 * 3.76 * root)
    ctx = DesignContext(member="H.3", section=sec, method=method, E=i["E"], G=11200.0, Fy=i["Fy"], Fu=65.0,
                        L=i["Lb_ft"] * FT, Lx=None, Ly=None, Lz=None, Lb=0.0)
    c, = hc.interaction(ctx, Demand(method, Pt=r["Pr"], extra={"Mrx": r["Mrx_kipft"] * FT, "Mry": r["Mry_kipft"] * FT}))
    assert c.name == HT and "Eq. H1-1b" in c.clause
    assert_printed(c.ratio, q["interaction"], "H1.2 check")
    assert c.ratio == pytest.approx(value, rel=1e-12) and c.passed


@pytest.mark.parametrize("method", ["LRFD", "ASD"])
def test_design_example_h4_compression_biaxial_flexure_with_b1(method):
    """Example H.4 (H-12 to H-15, PDF 274-277): W10x33; first-order moments amplified by B1 (Eq. A-8-3, A-8-5,
    Cm = 1.0), F2-2 with the listed Lp / Lr and Cb, F6-1, Eq. H1-1b with Pc from Manual Table 6-1."""
    ex = example("H.4")
    i, p = ex["inputs"], ex["printed"]
    r, q = i[method], p[method]
    Lc1 = i["Lc1_ft"] * FT
    Pe1x, Pe1y = hc.pe1_a_8_5(i["E"] * i["Ix"], Lc1), hc.pe1_a_8_5(i["E"] * i["Iy"], Lc1)      # Eq. A-8-5, EI* = EI
    assert_printed(Pe1x, p["Pe1x"], "Pe1x")
    assert_printed(Pe1y, p["Pe1y"], "Pe1y")
    B1x = hc.b1_a_8_3(i["Cm"], r["alpha"], r["Pr"], Pe1x)                   # Eq. A-8-3
    B1y = hc.b1_a_8_3(i["Cm"], r["alpha"], r["Pr"], Pe1y)
    assert_printed(B1x, q["B1x"], "B1x")
    assert_printed(B1y, q["B1y"], "B1y")
    Mrx, Mry = B1x * r["Mx1_kipft"], B1y * r["My1_kipft"]
    assert_printed(Mrx, q["Mrx_kipft"], "B1x Mx")
    assert_printed(Mry, q["Mry_kipft"], "B1y My")
    factor = 0.90 if method == "LRFD" else 1.0 / 1.67
    Mp = i["Fy"] * i["Zx"]
    assert_printed(Mp, p["Mp_kipin"], "Mp")
    Mnx = ff.inelastic_ltb(i["Cb"], Mp, i["Fy"], i["Sx"], i["Lb_ft"], i["Lp_ft"], i["Lr_ft"])            # Eq. F2-2
    assert_printed(Mnx, p["Mnx_kipin"], "Mnx")
    assert_printed(Mnx / FT, p["Mnx_kipft"], "Mnx")
    Mny = ff.mp_minor(i["Fy"], i["Zy"], i["Sy"])                            # Eq. F6-1
    assert_printed(Mny, p["Mny_kipin"], "Mny")
    assert_printed(Mny / FT, p["Mny_kipft"], "Mny")
    assert_printed(factor * Mnx / FT, q["Mcx_kipft"], "Mcx")
    assert_printed(factor * Mny / FT, q["Mcy_kipft"], "Mcy")
    assert_printed(r["Pr"] / r["Pc"], q["Pr_over_Pc"], "Pr/Pc")
    value, eq = hc.h1_interaction(r["Pr"], r["Pc"], Mrx * FT, factor * Mnx, Mry * FT, factor * Mny)
    assert eq == "H1-1b"
    assert_printed(value, q["interaction"], "Eq. H1-1b")
