"""AISC 360-22 Chapter F (flexure): F1 (Cb), F2, F3, F6, F7, F8.

Tier F: (a) hand calculations from the Specification (pp. 16.1-50 to 16.1-65) as
explicit arithmetic in the test (independent of the engine's functions), with the
clause of every step; (b) AISC Design Examples v16.0 F.1-1B, F.1-2B, F.1-3B, F.3B,
F.5, F.7B, F.8B, F.9B (tests/design_examples/part2; the ONE tolerance policy of
tests/design_examples).

Units of the hand tests: N, mm (MPa); of the Design Examples: kip, in.
"""
import math

import numpy as np
import pytest

from benchfea import design
from benchfea.design import DesignContext, DesignSection, Demand
from benchfea.design.aisc360_22 import f_flexure as ff
from design_examples.part2 import FT, assert_printed, example

E, G, FY = 200000.0, 77200.0, 345.0     # MPa
ROOT = math.sqrt(E / FY)                # 24.0772
PI = math.pi
Y, FLB, WLB, LTB = ff.NAME_Y, ff.NAME_FLB, ff.NAME_WLB, ff.NAME_LTB

# hand section (made-up numbers, not a catalogue shape): compact flange (7 <= 0.38 ROOT = 9.149,
# Table B4.1b Case 10) and compact web (40 <= 3.76 ROOT = 90.53, Case 15)
ZX, SX, RY, RTS, HO, J = 1.2e6, 1.05e6, 40.0, 48.0, 380.0, 5.0e5
ZY, SY = 2.4e5, 1.4e5


def _i(**kw):
    base = dict(A=9000.0, rx=160.0, ry=RY, Zx=ZX, Sx=SX, Zy=ZY, Sy=SY, J=J, rts=RTS, ho=HO, flange_ratio=7.0,
                web_ratio=40.0)
    base.update(kw)
    return DesignSection.i_shape("rolled", **base)


def _ctx(sec=None, method="LRFD", L=6000.0, **kw):
    base = dict(member="B", section=sec if sec is not None else _i(), method=method, E=E, G=G, Fy=FY, Fu=450.0,
                L=L, Lx=L, Ly=L, Lz=L)
    base.update(kw)
    return DesignContext(**base)


def _m(Mrx=0.0, Mry=0.0, c="U1"):
    return Demand(c, extra={"Mrx": Mrx, "Mry": Mry})


def _stations(x, mz=None, vy=None, my=None, vz=None, n=None, t=None):
    """Demand with stations; forces columns N, Vy, Vz, T, My, Mz (dMz/dx = -Vy, dMy/dx = +Vz)."""
    x = np.asarray(x, float)
    f = np.zeros((len(x), 6))
    for col, v in ((0, n), (1, vy), (2, vz), (3, t), (4, my), (5, mz)):
        if v is not None:
            f[:, col] = v
    return Demand("U1", extra={"x": x, "forces": f})


def _uniform(L, w, n=11):
    """Simply supported beam, uniform load w: Mz = w x (L - x)/2, Vy = -dMz/dx = -w (L - 2x)/2."""
    x = np.linspace(0.0, L, n)
    return _stations(x, mz=w * x * (L - x) / 2.0, vy=-w * (L - 2.0 * x) / 2.0)


def _by_name(checks):
    return {c.name: c for c in checks}


def _detail(c, start):
    return next(v[0] for k, v in c.details.items() if k.startswith(start))


# ============================================================================== moment diagram (derivation check)
def test_moment_interpolation_is_exact_for_linearly_varying_load():
    """Simply supported beam, load increasing linearly from 0 to w0: equilibrium gives R_A = w0 L/6 and
    M(x) = w0 L x/6 - w0 x^3/(6 L) (a cubic); V = dM/dx.  With stations only at the tenth points the cubic
    through station moments and slopes must return M anywhere, and find the maximum at x = L/sqrt(3)."""
    L, w0 = 7000.0, 12.0

    def M(x):
        return w0 * L * x / 6.0 - w0 * x ** 3 / (6.0 * L)

    def dM(x):
        return w0 * L / 6.0 - w0 * x ** 2 / (2.0 * L)

    x = np.linspace(0.0, L, 11)
    d = ff.Diagram(x, M(x), dM(x))
    for xq in (123.4, 0.25 * L, L / 3.0, 5555.5, 0.999 * L):
        assert d.at(xq) == pytest.approx(M(xq), rel=1e-13)
    Mmax, x_max, lo, hi = d.extremes(0.0, L)
    assert x_max == pytest.approx(L / math.sqrt(3.0), rel=1e-12)            # dM/dx = 0
    assert Mmax == pytest.approx(w0 * L ** 2 / (9.0 * math.sqrt(3.0)), rel=1e-13) and lo == 0.0 and hi == Mmax
    assert Mmax > max(M(x))                                                 # the peak lies between two stations
    # through Demand.extra: the station shear gives the slope (dMz/dx = -Vy); minor axis dMy/dx = +Vz
    mo = ff.moments(_stations(x, mz=M(x), vy=-dM(x), my=-2.0 * M(x), vz=-2.0 * dM(x)))
    assert mo.Mrx == pytest.approx(Mmax, rel=1e-13) and mo.x_Mrx == pytest.approx(L / math.sqrt(3.0), rel=1e-12)
    assert mo.Mry == pytest.approx(2.0 * Mmax, rel=1e-13)
    # a jump (two stations at the same x, e.g. an applied moment): the larger magnitude is reported there
    dj = ff.Diagram([0.0, 1.0, 1.0, 2.0], [0.0, 10.0, -30.0, 0.0], [10.0, 10.0, 30.0, 30.0])
    assert dj.at(1.0) == -30.0 and dj.extremes(0.0, 2.0)[0] == 30.0
    with pytest.raises(ValueError):
        d.at(L + 1.0)


def test_demand_without_moment_data_gives_no_flexure_check_and_is_reported():
    assert ff.moments(Demand("U1", Pc=1.0)) is None and ff.check(_ctx(), Demand("U1", Pc=1.0)) == []
    res = design.check_member(_ctx(), [Demand("U1", Pc=1.0)])
    assert res.status == "incomplete" and any("carry no moments" in w for w in res.warnings)
    assert not any(c.name.startswith("flexure") for c in res.checks)
    # numbers given: flexure is evaluated (zero moment -> ratio 0)
    c = _by_name(ff.check(_ctx(Lb=0.0), _m()))[Y["x"]]
    assert c.status == "ok" and c.ratio == 0.0 and c.x is None


# ============================================================================== F1: Cb
def test_hand_f1_cb_values_of_the_user_note():
    # Eq. F1-1: Cb = 12.5 Mmax / (2.5 Mmax + 3 MA + 4 MB + 3 MC); User Note of F1(c):
    assert ff.cb_f1_1(1.0, 1.0, 1.0, 1.0) == pytest.approx(12.5 / 12.5)                 # uniform moment: 1.0
    # equal end moments, reverse curvature: +M ... -M linear: MA = MC = M/2, MB = 0 -> 12.5/5.5 = 2.27
    assert ff.cb_f1_1(1.0, 0.5, 0.0, 0.5) == pytest.approx(12.5 / (2.5 + 1.5 + 1.5)) == pytest.approx(2.2727, abs=5e-5)
    # one end moment zero: M ... 0 linear: MA = 0.75, MB = 0.5, MC = 0.25 -> 12.5/7.5 = 1.67
    assert ff.cb_f1_1(1.0, 0.75, 0.5, 0.25) == pytest.approx(12.5 / 7.5)
    assert ff.cb_f1_1(-1.0, -0.75, 0.5, -0.25) == pytest.approx(12.5 / 7.5)              # absolute values
    # from a diagram: linear +M to -M over the segment [0, 4000]
    x = np.array([0.0, 4000.0])
    d = ff.Diagram(x, [8.0e7, -8.0e7], [-4.0e4, -4.0e4])
    Cb, Mmax, MA, MB, MC, x_max = ff.segment_cb(d, 0.0, 4000.0)
    assert (Mmax, MA, MB, MC) == pytest.approx((8.0e7, 4.0e7, 0.0, 4.0e7), abs=1e-6) and Cb == pytest.approx(12.5 / 5.5)
    # simply supported beam, uniform load, unbraced over the span: M at the quarter points = 0.75 Mmax
    # -> Cb = 12.5 / (2.5 + 2.25 + 4 + 2.25) = 12.5/11 = 1.136
    mo = ff.moments(_uniform(6000.0, 10.0))
    seg, = ff.ltb_segments(_ctx(), mo)
    assert seg.Cb == pytest.approx(12.5 / 11.0, rel=1e-13) and (seg.a, seg.b, seg.Lb) == (0.0, 6000.0, 6000.0)
    assert seg.Mr == pytest.approx(10.0 * 6000.0 ** 2 / 8.0) and seg.x == pytest.approx(3000.0)
    # no moment in the segment: Eq. F1-1 is 0/0 -> the segment is not checked
    assert ff.ltb_segments(_ctx(), ff.moments(_uniform(6000.0, 0.0))) == []
    assert ff.segment_cb(ff.Diagram(x, [0.0, 0.0], [0.0, 0.0]), 0.0, 4000.0)[0] is None


def test_unbraced_segments_from_brace_points_and_cb_per_segment():
    # uniform load, top flange (in compression: Mz > 0) braced at the third points
    L, w = 9000.0, 10.0
    d = _uniform(L, w)
    Mmax = w * L ** 2 / 8.0
    segs = ff.ltb_segments(_ctx(L=L, braceTop=(1 / 3, 2 / 3)), ff.moments(d))
    assert [(round(s.a), round(s.b)) for s in segs] == [(0, 3000), (3000, 6000), (6000, 9000)]

    def M(x):
        return w * x * (L - x) / 2.0

    # centre segment: Mmax at midspan, MA = MC = M(3750), MB = Mmax
    cb_c = 12.5 * Mmax / (2.5 * Mmax + 3 * M(3750.0) + 4 * Mmax + 3 * M(3750.0))
    # end segment: Mmax = M(3000), MA = M(750), MB = M(1500), MC = M(2250)
    cb_e = 12.5 * M(3000.0) / (2.5 * M(3000.0) + 3 * M(750.0) + 4 * M(1500.0) + 3 * M(2250.0))
    assert 1.01 < cb_c < 1.02 and 1.45 < cb_e < 1.47
    assert [s.Cb for s in segs] == pytest.approx([cb_e, cb_c, cb_e], rel=1e-12)
    assert segs[0].Mr == pytest.approx(M(3000.0)) and segs[1].Mr == pytest.approx(Mmax)
    # the bottom flange is never in compression: bracing only the bottom flange leaves one segment (the member)
    segs = ff.ltb_segments(_ctx(L=L, braceBottom=(1 / 3, 2 / 3)), ff.moments(d))
    assert [(s.a, s.b) for s in segs] == [(0.0, L)]
    # hogging (Mz < 0): now the bottom flange is the compression flange
    hog = _stations(d.extra["x"], mz=-d.extra["forces"][:, 5], vy=-d.extra["forces"][:, 1])
    assert len(ff.ltb_segments(_ctx(L=L, braceBottom=(1 / 3, 2 / 3)), ff.moments(hog))) == 3
    assert len(ff.ltb_segments(_ctx(L=L, braceTop=(1 / 3, 2 / 3)), ff.moments(hog))) == 1


def test_reverse_curvature_counts_only_points_where_both_flanges_are_braced():
    # fixed-fixed beam, uniform load: Mz = w (6 L x - L^2 - 6 x^2)/12, hogging at the ends, sagging inside,
    # zero at x = 0.2113 L and 0.7887 L
    L, w = 10000.0, 12.0
    x = np.linspace(0.0, L, 21)
    d = _stations(x, mz=w * (6 * L * x - L ** 2 - 6 * x ** 2) / 12.0, vy=-w * (6 * L - 12 * x) / 12.0)
    mo = ff.moments(d)
    # top flange braced at 0.5 L only: in [0, 0.5 L] the moment changes sign -> the top-flange brace does not
    # brace the compression flange everywhere -> the segment between points where BOTH flanges are braced
    # (the member ends) is used for both flanges: one segment, the whole member
    segs = ff.ltb_segments(_ctx(L=L, braceTop=(0.5,)), mo)
    assert [(s.a, s.b) for s in segs] == [(0.0, L)] and "both flanges are braced" in segs[0].note
    assert segs[0].Mr == pytest.approx(w * L ** 2 / 12.0)                   # the end moment governs
    # Eq. F1-1 on the whole member: MA = MC = |M(L/4)| = w L^2/96, MB = w L^2/24
    cb = 12.5 / 12.0 / (2.5 / 12.0 + 3 / 96.0 + 4 / 24.0 + 3 / 96.0)
    assert segs[0].Cb == pytest.approx(cb, rel=1e-12) and 2.3 < cb < 2.4
    # both flanges braced at 0.5 L: two segments, each in reverse curvature, Eq. F1-1 applies (F1(c))
    segs = ff.ltb_segments(_ctx(L=L, braceTop=(0.5,), braceBottom=(0.5,)), mo)
    assert [(s.a, s.b) for s in segs] == [(0.0, 5000.0), (5000.0, L)] and all(not s.note for s in segs)
    # top flange braced at 0.3 and 0.7 (sagging throughout between them): that segment is checked on its own
    # with the top flange; the end regions fall back to the whole member
    segs = ff.ltb_segments(_ctx(L=L, braceTop=(0.3, 0.7)), mo)
    assert sorted((round(s.a), round(s.b)) for s in segs) == [(0, 10000), (3000, 7000)]


def test_given_lb_and_cb_and_free_end_rules():
    L, w = 6000.0, 10.0
    mo = ff.moments(_uniform(L, w))
    Mmax = w * L ** 2 / 8.0
    # Lb given, different from the member length: brace positions unknown -> Cb = 1.0, largest moment
    seg, = ff.ltb_segments(_ctx(Lb=2000.0), mo)
    assert (seg.Lb, seg.Cb, seg.a) == (2000.0, 1.0, None) and seg.Mr == pytest.approx(Mmax) and "not known" in seg.cb_text
    # Cb given: used as given in every segment
    seg, = ff.ltb_segments(_ctx(Lb=2000.0, Cb=1.3), mo)
    assert seg.Cb == 1.3 and seg.cb_text == "Cb given"
    assert [s.Cb for s in ff.ltb_segments(_ctx(braceTop=(0.5,), Cb=1.3), mo)] == [1.3, 1.3]
    # numbers only (no diagram): Cb = 1.0 unless given; Lb = the longest segment between brace points
    mo_n = ff.moments(_m(Mrx=5.0e7))
    # (which flange is in compression is not known without a diagram: the longer of the two flanges' segments)
    seg, = ff.ltb_segments(_ctx(braceTop=(0.25,), braceBottom=(0.25,)), mo_n)
    assert (seg.Lb, seg.Cb) == (4500.0, 1.0) and "no moment diagram" in seg.cb_text
    assert ff.ltb_segments(_ctx(braceTop=(0.25,)), mo_n)[0].Lb == 6000.0
    assert ff.ltb_segments(_ctx(Cb=1.14), mo_n)[0].Cb == 1.14
    # cantilever: tip free; Mz = -P (L - x) (hogging), Vy = +P ... built-in root -> F1: Cb = 1.0, Lb = L
    x = np.linspace(0.0, L, 11)
    cant = ff.moments(_stations(x, mz=-1.0e4 * (L - x), vy=np.full(11, -1.0e4)))
    seg, = ff.ltb_segments(_ctx(end_free=(False, True), end_built_in=(True, False)), cant)
    assert (seg.Cb, seg.Lb, seg.blocked) == (1.0, L, "") and "cantilever" in seg.cb_text and "ASSUMED" in seg.cb_text
    # the same diagram between braced ends would give Cb = 1.67 by Eq. F1-1 -- not used for a free end
    assert ff.ltb_segments(_ctx(), cant)[0].Cb == pytest.approx(12.5 / 7.5)
    # free end but the other end is not built in (an overhang): no rule in F1 -> blocked
    seg, = ff.ltb_segments(_ctx(end_free=(False, True)), cant)
    assert seg.Cb is None and "free end" in seg.blocked
    # free end beyond a brace point: the end segment is blocked, the inner one is a normal segment
    segs = ff.ltb_segments(_ctx(end_free=(False, True), end_built_in=(True, False), braceBottom=(0.5,)), cant)
    assert [bool(s.blocked) for s in segs] == [False, True]
    # ... unless Cb is given
    assert all(not s.blocked for s in ff.ltb_segments(_ctx(end_free=(False, True), Cb=1.0), cant))
    # the blocked segment makes LTB notChecked (never a silent pass)
    c = _by_name(ff.check(_ctx(end_free=(False, True)), _stations(x, mz=-1.0e4 * (L - x), vy=np.full(11, -1.0e4))))[LTB]
    assert c.status == "notChecked" and "free end" in c.reason and c.demand == pytest.approx(6.0e7)


# ============================================================================== tier F (a): F2
def test_hand_f2_yielding():
    # F2.1 Eq. F2-1: Mn = Mp = Fy Zx = 345 * 1.2e6 = 414.0e6 N mm; F1(a): phi_b = 0.90, Omega_b = 1.67
    Mp = 345.0 * 1.2e6
    assert Mp == 414.0e6
    c = _by_name(ff.check(_ctx(Lb=0.0), _m(Mrx=2.0e8)))[Y["x"]]
    assert c.nominal == Mp and c.capacity == pytest.approx(0.90 * Mp, rel=1e-15)
    assert c.ratio == pytest.approx(2.0e8 / (0.90 * 414.0e6), rel=1e-14) and c.passed and c.dimension == "moment"
    assert c.clause == "AISC 360-22 F2.1, Eq. F2-1; B4.1, Table B4.1b; F1(a); B3.1 Eq. B3-1"
    assert c.details["class flange"][0] == "compact" and c.details["class web"][0] == "compact"
    assert "does not apply" in c.details["LTB"][0]                          # Lb = 0 <= Lp: F2.2(a)
    assert LTB not in _by_name(ff.check(_ctx(Lb=0.0), _m(Mrx=2.0e8)))
    a = _by_name(ff.check(_ctx(method="ASD", Lb=0.0), _m(Mrx=2.0e8)))[Y["x"]]
    assert a.capacity == pytest.approx(Mp / 1.67, rel=1e-15) and (a.factor, a.factor_name) == (1.67, "Omega")
    assert a.clause.endswith("F1(a); B3.2 Eq. B3-2")


def _lp_lr():
    Lp = 1.76 * 40.0 * math.sqrt(200000.0 / 345.0)                          # Eq. F2-5: 1695.03 mm
    X = 5.0e5 * 1.0 / (1.05e6 * 380.0)                                      # Jc/(Sx ho) = 1.25313e-3
    Lr = 1.95 * 48.0 * 200000.0 / (0.7 * 345.0) * math.sqrt(X + math.sqrt(X ** 2 + 6.76 * (0.7 * 345.0 / 200000.0) ** 2))
    return Lp, X, Lr                                                        # Eq. F2-6: 5276.46 mm


def test_hand_f2_inelastic_ltb():
    Lp, X, Lr = _lp_lr()
    assert 1695.0 < Lp < 1695.1 and 1.2531e-3 < X < 1.2532e-3 and 5276.4 < Lr < 5276.5
    Mp = 414.0e6
    # F2.2(b): Lp < Lb = 4000 <= Lr -> Eq. F2-2 with Cb = 1.0:
    # Mn = Mp - (Mp - 0.7 Fy Sx)(Lb - Lp)/(Lr - Lp), 0.7 Fy Sx = 253.575e6
    Mn = Mp - (Mp - 0.7 * 345.0 * 1.05e6) * (4000.0 - Lp) / (Lr - Lp)
    assert 310.7e6 < Mn < 310.8e6
    c = _by_name(ff.check(_ctx(Lb=4000.0, Cb=1.0), _m(Mrx=2.0e8)))[LTB]
    assert c.nominal == pytest.approx(Mn, rel=1e-13) and c.capacity == pytest.approx(0.9 * Mn, rel=1e-13)
    assert c.clause == ("AISC 360-22 F2.2, Eq. F2-2, F2-5, F2-6, F2-8a; F2.1, Eq. F2-1; B4.1, Table B4.1b; F1; "
                        "F1(a); B3.1 Eq. B3-1")
    assert _detail(c, "Lp") == pytest.approx(Lp, rel=1e-14) and _detail(c, "Lr") == pytest.approx(Lr, rel=1e-14)
    assert c.inputs["Lb"][0] == 4000.0 and c.inputs["Cb"] == (1.0, "dimensionless", "given")
    assert c.inputs["rts"][0] == 48.0 and c.inputs["ho"][0] == 380.0
    # Cb = 1.2: Cb x the bracket, 372.9e6 < Mp
    c = _by_name(ff.check(_ctx(Lb=4000.0, Cb=1.2), _m(Mrx=2.0e8)))[LTB]
    assert c.nominal == pytest.approx(1.2 * Mn, rel=1e-13)
    # Cb = 1.5: 466.1e6 > Mp -> capped at Mp (the "<= Mp" of Eq. F2-2)
    c = _by_name(ff.check(_ctx(Lb=4000.0, Cb=1.5), _m(Mrx=2.0e8)))[LTB]
    assert c.nominal == Mp
    # at Lb = Lr with Cb = 1 Eq. F2-2 gives 0.7 Fy Sx; just beyond, Eq. F2-3 continues from there.  Eq. F2-6 is
    # Eq. F2-4 solved for Fcr = 0.7 Fy (algebra done here: Lr = pi^2 sqrt(0.039) rts E/(0.7 Fy) sqrt(X +
    # sqrt(X^2 + (0.7 Fy/E)^2/(0.039^2 pi^4)))), with pi^2 sqrt(0.039) = 1.9491 printed as 1.95 and
    # 1/(0.039^2 pi^4) = 6.7495 printed as 6.76 -- so the two branches meet within about 0.1 %, not exactly
    assert PI ** 2 * math.sqrt(0.039) == pytest.approx(1.95, abs=1e-3)
    assert 1.0 / (0.039 ** 2 * PI ** 4) == pytest.approx(6.76, abs=0.011)
    Mn_r = ff.ltb_i_shape(E, FY, Mp, SX, RY, RTS, J, 1.0, HO, Lr, 1.0)
    assert Mn_r[1] == "F2-2" and Mn_r[0] == pytest.approx(0.7 * FY * SX, rel=1e-13)
    Mn_e = ff.ltb_i_shape(E, FY, Mp, SX, RY, RTS, J, 1.0, HO, Lr * (1 + 1e-9), 1.0)
    assert Mn_e[1] == "F2-3, F2-4" and Mn_e[0] == pytest.approx(0.7 * FY * SX, rel=2e-3)
    # Lb <= Lp: LTB does not apply
    assert ff.ltb_i_shape(E, FY, Mp, SX, RY, RTS, J, 1.0, HO, Lp, 1.0)[:2] == (None, "F2.2(a)")


def test_hand_f2_elastic_ltb_and_the_classical_formula():
    Lp, X, Lr = _lp_lr()
    # F2.2(c): Lb = 9000 > Lr -> Eq. F2-4: Fcr = Cb pi^2 E/(Lb/rts)^2 sqrt(1 + 0.078 X (Lb/rts)^2), (Lb/rts) = 187.5
    lam2 = (9000.0 / 48.0) ** 2
    Fcr = PI ** 2 * 200000.0 / lam2 * math.sqrt(1.0 + 0.078 * X * lam2)
    assert 118.2 < Fcr < 118.3
    Mn = Fcr * 1.05e6                                                       # Eq. F2-3: 124.17e6 <= Mp
    assert 124.1e6 < Mn < 124.2e6
    c = _by_name(ff.check(_ctx(L=9000.0, Lb=9000.0, Cb=1.0), _m(Mrx=1.0e8)))[LTB]
    assert c.nominal == pytest.approx(Mn, rel=1e-13) and "Eq. F2-3, F2-4" in c.clause
    assert _detail(c, "Fcr") == pytest.approx(Fcr, rel=1e-13) and c.passed is True
    assert c.ratio == pytest.approx(1.0e8 / (0.9 * Mn), rel=1e-13)
    # User Note of F2.2: Eq. F2-3 / F2-4 give the same as  Mcr = Cb (pi/Lb) sqrt(E Iy G J + (pi E/Lb)^2 Iy Cw).
    # Algebra (done here): with rts^2 = sqrt(Iy Cw)/Sx (Eq. F2-7) and Cw = Iy ho^2/4 the two agree exactly when
    # G = 0.039 pi^2 E (0.078 = 2 x 0.039), i.e. G/E = 0.3849 -- the value 11,200/29,000 = 0.3862 rounded.
    Iy, ho, Sx, Jx, Lb = 2.0e7, 380.0, 1.05e6, 5.0e5, 9000.0
    Cw = Iy * ho ** 2 / 4.0
    rts = math.sqrt(math.sqrt(Iy * Cw) / Sx)                                # Eq. F2-7
    assert ff.rts_f2_7(Iy, Cw, Sx) == pytest.approx(rts, rel=1e-15)
    Gx = 0.039 * PI ** 2 * E
    Mcr = PI / Lb * math.sqrt(E * Iy * Gx * Jx + (PI * E / Lb) ** 2 * Iy * Cw)
    assert ff.fcr_f2_4(1.0, E, Lb, rts, Jx, 1.0, Sx, ho) * Sx == pytest.approx(Mcr, rel=1e-13)
    Mcr_g = PI / Lb * math.sqrt(E * Iy * (11200.0 / 29000.0 * E) * Jx + (PI * E / Lb) ** 2 * Iy * Cw)
    assert ff.fcr_f2_4(1.0, E, Lb, rts, Jx, 1.0, Sx, ho) * Sx == pytest.approx(Mcr_g, rel=2e-3)
    # the engine derives rts by Eq. F2-7 and ho = d - tf when they are not given
    sec = _i(rts=None, ho=None, Iy=Iy, ry=None, Cw=Cw, d=392.0, tf=12.0, bf=168.0, tw=7.6, flange_ratio=None,
             web_ratio=None, r=14.0)
    c = _by_name(ff.check(_ctx(sec, L=Lb, Lb=Lb, Cb=1.0), _m(Mrx=1.0e8)))[LTB]
    assert _detail(c, "rts (Eq. F2-7)") == pytest.approx(rts) and _detail(c, "ho = d - tf") == 380.0
    assert c.nominal == pytest.approx(Mcr, rel=1e-13) and "F2-7" in c.clause


def test_ltb_from_stations_governing_segment_and_summary():
    # uniform load, top flange braced at midspan: two equal segments, Mmax at the brace; Cb by Eq. F1-1
    L, w = 8000.0, 30.0
    d = _uniform(L, w)
    Mmax = w * L ** 2 / 8.0                                                 # 240e6

    def M(x):
        return w * x * (L - x) / 2.0

    Cb = 12.5 * Mmax / (2.5 * Mmax + 3 * M(1000.0) + 4 * M(2000.0) + 3 * M(3000.0))     # 1.2987
    assert Cb == pytest.approx(12.5 / (2.5 + 3 * 0.4375 + 4 * 0.75 + 3 * 0.9375))
    Lp, X, Lr = _lp_lr()
    Mn = min(414.0e6, Cb * (414.0e6 - (414.0e6 - 0.7 * 345.0 * 1.05e6) * (4000.0 - Lp) / (Lr - Lp)))
    c = _by_name(ff.check(_ctx(L=L, braceTop=(0.5,)), d))[LTB]
    assert c.nominal == pytest.approx(Mn, rel=1e-12) and c.demand == pytest.approx(Mmax) and c.x == pytest.approx(4000.0)
    assert c.inputs["Cb"] == (pytest.approx(Cb), "dimensionless", "derived: Eq. F1-1")
    assert c.inputs["Lb"] == (pytest.approx(4000.0), "length", "bracing: segment between brace points")
    assert "Eq. F1-1" in c.clause and "from 0 to 4000" in c.details["all segments"][0]
    assert "from 4000 to 8000" in c.details["all segments"][0] and any("Eq. F1-1" in n for n in c.notes)
    # yielding is checked with the largest moment of the member
    y = _by_name(ff.check(_ctx(L=L, braceTop=(0.5,)), d))[Y["x"]]
    assert y.demand == pytest.approx(Mmax) and y.x == pytest.approx(4000.0)
    # data for Lr missing while a segment is longer than Lp: notChecked, not a silent pass
    c = _by_name(ff.check(_ctx(_i(rts=None, J=None), L=L), d))[LTB]
    assert c.status == "notChecked" and "rts" in c.reason and "J" in c.reason and c.demand == pytest.approx(Mmax)
    # ... but not needed when every segment is within Lp
    names = _by_name(ff.check(_ctx(_i(rts=None, J=None), L=L, braceTop=tuple(k / 8 for k in range(1, 8))), d))
    assert LTB not in names and names[Y["x"]].status == "ok"


# ============================================================================== tier F (a): F3
def test_hand_f3_flange_local_buckling():
    # flange lambda = 12: 0.38 ROOT = 9.149 < 12 <= 1.0 ROOT = 24.077 -> noncompact (Table B4.1b Case 10)
    lam_p, lam_r = 0.38 * math.sqrt(200000.0 / 345.0), math.sqrt(200000.0 / 345.0)
    Mp = 414.0e6
    # F3.2(a) Eq. F3-1: Mn = Mp - (Mp - 0.7 Fy Sx)(lam - lam_pf)/(lam_rf - lam_pf) = 383.36e6
    Mn = Mp - (Mp - 0.7 * 345.0 * 1.05e6) * (12.0 - lam_p) / (lam_r - lam_p)
    assert 383.3e6 < Mn < 383.4e6
    checks = _by_name(ff.check(_ctx(_i(flange_ratio=12.0), Lb=0.0), _m(Mrx=3.0e8)))
    c = checks[FLB["x"]]
    assert c.nominal == pytest.approx(Mn, rel=1e-13) and Y["x"] not in checks      # F3: LTB and FLB only
    assert c.clause == "AISC 360-22 F3.2, Eq. F3-1; F2.1, Eq. F2-1; B4.1, Table B4.1b; F1(a); B3.1 Eq. B3-1"
    assert c.details["class flange"][0] == "noncompact"
    # LTB by F2.2 also applies to F3 members (F3.1)
    c = _by_name(ff.check(_ctx(_i(flange_ratio=12.0), Lb=4000.0, Cb=1.0), _m(Mrx=3.0e8)))[LTB]
    assert c.clause.startswith("AISC 360-22 F3.1; F2.2, Eq. F2-2") and 310.7e6 < c.nominal < 310.8e6
    # slender flange, lambda = 30 > 24.077: F3.2(b) Eq. F3-2: Mn = 0.9 E kc Sx / lam^2,
    # kc = 4/sqrt(h/tw) = 4/sqrt(40) = 0.6325 (between 0.35 and 0.76)
    kc = 4.0 / math.sqrt(40.0)
    Mn = 0.9 * 200000.0 * kc * 1.05e6 / 30.0 ** 2
    assert 0.632 < kc < 0.633 and 132.8e6 < Mn < 132.9e6
    c = _by_name(ff.check(_ctx(_i(flange_ratio=30.0), Lb=0.0), _m(Mrx=1.0e8)))[FLB["x"]]
    assert c.nominal == pytest.approx(Mn, rel=1e-13) and "Eq. F3-2" in c.clause
    assert _detail(c, "kc") == pytest.approx(kc) and c.details["class flange"][0] == "slender"
    # kc limits: h/tw = 20 -> 4/sqrt(20) = 0.894 -> 0.76 (web stays compact)
    c = _by_name(ff.check(_ctx(_i(flange_ratio=30.0, web_ratio=20.0), Lb=0.0), _m(Mrx=1.0e8)))[FLB["x"]]
    assert c.nominal == pytest.approx(0.9 * E * 0.76 * SX / 900.0, rel=1e-13)
    # at lambda = lambda_rf Eq. F3-1 gives 0.7 Fy Sx
    assert ff.flb_i_major(E, FY, Mp, SX, lam_r, lam_p, lam_r, kc) == (pytest.approx(0.7 * FY * SX), "F3-1")


# ============================================================================== tier F (a): F6
def test_hand_f6_minor_axis():
    # F6.1 Eq. F6-1: Mp = Fy Zy <= 1.6 Fy Sy: 345 * 2.4e5 = 82.8e6 > 1.6 * 345 * 1.4e5 = 77.28e6 -> 77.28e6
    assert 345.0 * 2.4e5 == pytest.approx(82.8e6) and 1.6 * 345.0 * 1.4e5 == pytest.approx(77.28e6)
    Mp = 77.28e6
    checks = _by_name(ff.check(_ctx(), _m(Mry=5.0e7)))
    c = checks[Y["y"]]
    assert c.nominal == pytest.approx(Mp, rel=1e-15) and FLB["y"] not in checks        # compact flange: no FLB
    assert c.clause == "AISC 360-22 F6.1, Eq. F6-1; B4.1, Table B4.1b; F1(a); B3.1 Eq. B3-1"
    assert c.ratio == pytest.approx(5.0e7 / (0.9 * Mp), rel=1e-14)
    # the cap not governing: Zy = 2.0e5 -> 69.0e6 < 77.28e6
    assert _by_name(ff.check(_ctx(_i(Zy=2.0e5)), _m(Mry=1.0)))[Y["y"]].nominal == pytest.approx(69.0e6)
    # noncompact flange (lambda = 12, Table B4.1b Case 13: 9.149 / 24.077):
    # F6.2(b) Eq. F6-2: Mn = Mp - (Mp - 0.70 Fy Sy)(lam - lam_pf)/(lam_rf - lam_pf) = 68.98e6
    lam_p, lam_r = 0.38 * ROOT, ROOT
    Mn = Mp - (Mp - 0.70 * 345.0 * 1.4e5) * (12.0 - lam_p) / (lam_r - lam_p)
    assert 68.9e6 < Mn < 69.0e6
    checks = _by_name(ff.check(_ctx(_i(flange_ratio=12.0), Lb=0.0), _m(Mry=5.0e7)))
    assert checks[FLB["y"]].nominal == pytest.approx(Mn, rel=1e-13) and "Eq. F6-2" in checks[FLB["y"]].clause
    assert checks[Y["y"]].nominal == pytest.approx(Mp)                      # F6 lists yielding and FLB
    # slender flange (lambda = 30): F6.2(c) Eq. F6-3, F6-4: Mn = 0.70 E/(b/tf)^2 Sy = 21.78e6
    Mn = 0.70 * 200000.0 / 30.0 ** 2 * 1.4e5
    assert 21.7e6 < Mn < 21.8e6
    c = _by_name(ff.check(_ctx(_i(flange_ratio=30.0), Lb=0.0), _m(Mry=1.0e7)))[FLB["y"]]
    assert c.nominal == pytest.approx(Mn, rel=1e-13) and "Eq. F6-3, F6-4" in c.clause
    # no LTB and no web limit state for minor-axis bending: an unbraced length changes nothing
    assert _by_name(ff.check(_ctx(Lb=20000.0, L=20000.0), _m(Mry=5.0e7)))[Y["y"]].nominal == pytest.approx(Mp)


# ============================================================================== tier F (a): F7
def _box(kind="rectHSS", **kw):
    """300 x 200 x 8 (H x B x t): A, I, S, Z of the sharp-cornered tube computed here."""
    H, B, t = 300.0, 200.0, 8.0
    hi, bi = H - 2 * t, B - 2 * t
    base = dict(A=H * B - hi * bi, H=H, B=B, t=t, Ix=(B * H ** 3 - bi * hi ** 3) / 12.0,
                Iy=(H * B ** 3 - hi * bi ** 3) / 12.0, Zx=(B * H ** 2 - bi * hi ** 2) / 4.0,
                Zy=(H * B ** 2 - hi * bi ** 2) / 4.0, J=2.0 * t * ((H - t) * (B - t)) ** 2 / ((H - t) + (B - t)))
    base["Sx"], base["Sy"] = base["Ix"] / (H / 2.0), base["Iy"] / (B / 2.0)
    base.update(kw)
    return DesignSection.rectangular(kind, **base)


def test_hand_f7_yielding_and_local_buckling():
    sec = _box()                                    # HSS: flat widths = outside - 3 t: b/t = 22, h/t = 34.5
    assert sec.elements["wall_z"].ratio == 22.0 and sec.elements["wall_y"].ratio == 34.5
    # Table B4.1b Case 17 flange: 1.12 ROOT = 26.97 / 1.40 ROOT = 33.71; Case 19 web: 2.42 ROOT = 58.27 / 137.2
    # -> major axis: flange b/t = 22 and web h/t = 34.5 compact: only yielding, Eq. F7-1: Mn = Fy Z
    # -> minor axis: the flanges are the depth walls, h/t = 34.5 > 33.71: slender (flange local buckling as well)
    checks = _by_name(ff.check(_ctx(sec, Lb=0.0), _m(Mrx=1.0e8, Mry=5.0e7)))
    assert set(checks) == {Y["x"], Y["y"], FLB["y"]}
    assert checks[FLB["y"]].details["class flange (depth wall)"][0] == "slender"
    assert checks[Y["x"]].nominal == pytest.approx(345.0 * sec.Zx, rel=1e-15)
    assert checks[Y["y"]].nominal == pytest.approx(345.0 * sec.Zy, rel=1e-15)
    assert checks[Y["x"]].clause == "AISC 360-22 F7.1, Eq. F7-1; B4.1, Table B4.1b; F1(a); B3.1 Eq. B3-1"
    # noncompact flange: b/t = 30 (given): F7.2(b) Eq. F7-2: Mn = Mp - (Mp - Fy S)(lam - lam_p)/(lam_r - lam_p)
    sec = _box(b_ratio=30.0)
    Mp, lam_p, lam_r = 345.0 * sec.Zx, 1.12 * ROOT, 1.40 * ROOT
    Mn = Mp - (Mp - 345.0 * sec.Sx) * (30.0 - lam_p) / (lam_r - lam_p)
    assert 0.8 * Mp < Mn < Mp
    c = _by_name(ff.check(_ctx(sec, Lb=0.0), _m(Mrx=1.0e8)))[FLB["x"]]
    assert c.nominal == pytest.approx(Mn, rel=1e-13) and "F7.2(b), Eq. F7-2" in c.clause
    # box section: the flange limit lambda_r is 1.49 ROOT = 35.87 (Case 21): b/t = 35 is still noncompact
    c = _by_name(ff.check(_ctx(_box("box", b_ratio=35.0), Lb=0.0), _m(Mrx=1.0e8)))[FLB["x"]]
    assert c.nominal == pytest.approx(Mp - (Mp - 345.0 * sec.Sx) * (35.0 - lam_p) / (1.49 * ROOT - lam_p), rel=1e-13)
    # bending about the minor axis: the flanges are the depth walls (h/t), the webs the width walls
    c = _by_name(ff.check(_ctx(_box(h_ratio=30.0), Lb=0.0), _m(Mry=1.0e7)))
    Mpy = 345.0 * sec.Zy
    assert c[FLB["y"]].nominal == pytest.approx(Mpy - (Mpy - 345.0 * sec.Sy) * (30.0 - lam_p) / (lam_r - lam_p), rel=1e-13)
    # noncompact web: h/t = 80 (58.27 < 80 <= 137.2): F7.3(b) Eq. F7-6 with the web limits
    lam_pw, lam_rw = 2.42 * ROOT, 5.70 * ROOT
    c = _by_name(ff.check(_ctx(_box(h_ratio=80.0), Lb=0.0), _m(Mrx=1.0e8)))[WLB["x"]]
    assert c.nominal == pytest.approx(Mp - (Mp - 345.0 * sec.Sx) * (80.0 - lam_pw) / (lam_rw - lam_pw), rel=1e-13)
    assert "F7.3(b), Eq. F7-6" in c.clause
    # slender web: Eq. F7-7 is not implemented -> notChecked with the reason
    c = _by_name(ff.check(_ctx(_box("box", h_ratio=150.0), Lb=0.0), _m(Mrx=1.0e8)))[WLB["x"]]
    assert c.status == "notChecked" and "F7-7" in c.reason and c.demand == 1.0e8


def test_hand_f7_slender_flange_effective_width():
    # HSS 300 x 200 x 4: b = 200 - 3*4 = 188, b/t = 47 > 1.40 ROOT = 33.71 -> slender flange (major axis)
    H, B, t = 300.0, 200.0, 4.0
    sec = _box(H=H, B=B, t=t, A=3936.0, Ix=5.2e7, Zx=4.0e5, Sx=5.2e7 / 150.0, Iy=2.8e7, Zy=3.0e5, Sy=2.8e5)
    assert sec.elements["wall_z"].ratio == 47.0
    root = math.sqrt(200000.0 / 345.0)
    # F7.2(c) Eq. F7-4: be = 1.92 t sqrt(E/Fy) [1 - 0.38/(b/t) sqrt(E/Fy)] = 184.91 * 0.80533 = 148.9 <= b
    be = 1.92 * 4.0 * root * (1.0 - 0.38 / 47.0 * root)
    assert 148.8 < be < 149.0
    # effective section modulus, ineffective width removed from both flanges (engine method):
    # Ieff = I - 2 [(b - be) t^3/12 + (b - be) t ((H - t)/2)^2], Se = Ieff/(H/2); Eq. F7-3: Mn = Fy Se
    lost = 188.0 - be
    Ieff = 5.2e7 - 2.0 * (lost * 4.0 ** 3 / 12.0 + lost * 4.0 * ((300.0 - 4.0) / 2.0) ** 2)
    Mn = 345.0 * Ieff / 150.0
    c = _by_name(ff.check(_ctx(sec, Lb=0.0), _m(Mrx=5.0e7)))[FLB["x"]]
    assert c.nominal == pytest.approx(Mn, rel=1e-13) and "F7.2(c), Eq. F7-3, F7-4" in c.clause
    assert _detail(c, "be") == pytest.approx(be, rel=1e-14) and any("both flanges" in n for n in c.notes)
    # box section: Eq. F7-5 (0.34) and clear width b = B - 2 t = 192, b/t = 48 > 1.49 ROOT
    be_box = 1.92 * 4.0 * root * (1.0 - 0.34 / 48.0 * root)
    box = _box("box", H=H, B=B, t=t, A=3936.0, Ix=5.2e7, Zx=4.0e5, Sx=5.2e7 / 150.0, Iy=2.8e7, Zy=3.0e5, Sy=2.8e5)
    c = _by_name(ff.check(_ctx(box, Lb=0.0), _m(Mrx=5.0e7)))[FLB["x"]]
    assert _detail(c, "be") == pytest.approx(be_box, rel=1e-14) and "F7-5" in c.clause
    # at b/t = lambda_r = 1.40 sqrt(E/Fy) Eq. F7-4 returns be = 0.999 b (the flange is just fully effective)
    assert ff.effective_width_f7(1.40 * root * 4.0, 4.0, E, FY, True) == pytest.approx(1.40 * root * 4.0, rel=2e-3)
    # ratio only (no dimensions): the effective modulus cannot be formed -> notChecked
    c = _by_name(ff.check(_ctx(DesignSection.rectangular("rectHSS", A=3936.0, Zx=4.0e5, Sx=3.4e5, b_ratio=47.0,
                                                         h_ratio=60.0), Lb=0.0), _m(Mrx=5.0e7)))[FLB["x"]]
    assert c.status == "notChecked" and "effective section modulus" in c.reason


def test_hand_f7_lateral_torsional_buckling():
    sec = _box()
    Mp = 345.0 * sec.Zx
    # F7.4: Lp = 0.13 E ry sqrt(J Ag)/Mp (Eq. F7-10); Lr = 2 E ry sqrt(J Ag)/(0.7 Fy Sx) (Eq. F7-11)
    ry = math.sqrt(sec.Iy / sec.A)
    rootJA = math.sqrt(sec.J * sec.A)
    Lp = 0.13 * 200000.0 * ry * rootJA / Mp
    Lr = 2.0 * 200000.0 * ry * rootJA / (0.7 * 345.0 * sec.Sx)
    assert 4000.0 < Lp < 20000.0 and Lr > 20 * Lp
    # Lb <= Lp: no LTB check (F7.4(a))
    assert LTB not in _by_name(ff.check(_ctx(sec, L=Lp, Lb=Lp, Cb=1.0), _m(Mrx=1.0e8)))
    # Lp < Lb <= Lr: Eq. F7-8 = Cb [Mp - (Mp - 0.7 Fy Sx)(Lb - Lp)/(Lr - Lp)] <= Mp
    Lb = 5.0 * Lp
    Mn = Mp - (Mp - 0.7 * 345.0 * sec.Sx) * (Lb - Lp) / (Lr - Lp)
    c = _by_name(ff.check(_ctx(sec, L=Lb, Lb=Lb, Cb=1.0), _m(Mrx=1.0e8)))[LTB]
    assert c.nominal == pytest.approx(Mn, rel=1e-13) and Mn < Mp
    assert c.clause == "AISC 360-22 F7.4, Eq. F7-8, F7-10, F7-11; F7.1, Eq. F7-1; F1; F1(a); B3.1 Eq. B3-1"
    assert _detail(c, "Lp") == pytest.approx(Lp, rel=1e-14) and _detail(c, "Lr") == pytest.approx(Lr, rel=1e-14)
    # Lb > Lr: Eq. F7-9: Mn = 2 E Cb sqrt(J Ag)/(Lb/ry) <= Mp; at Lb = Lr this equals 0.7 Fy Sx (continuity)
    assert ff.ltb_rectangular(E, FY, Mp, sec.Sx, ry, sec.J, sec.A, Lr * (1 + 1e-12), 1.0)[0] == \
        pytest.approx(0.7 * FY * sec.Sx, rel=1e-9)
    Lb = 2.0 * Lr
    Mn, eq, _, _ = ff.ltb_rectangular(E, FY, Mp, sec.Sx, ry, sec.J, sec.A, Lb, 1.0)
    assert eq == "F7-9" and Mn == pytest.approx(2.0 * 200000.0 * rootJA / (Lb / ry), rel=1e-14)
    # minor-axis bending of this section (Iy < Ix): no LTB (F7.4 User Note), whatever Lb
    checks = _by_name(ff.check(_ctx(sec, L=50000.0, Lb=50000.0), _m(Mry=1.0e7)))
    assert set(checks) == {Y["x"], Y["y"], FLB["y"], LTB} and checks[LTB].demand == 0.0     # LTB: major axis only
    assert "smaller or equal inertia" in checks[Y["y"]].details["LTB"][0]
    # a square section has no LTB about either axis
    sq = _box(H=200.0, B=200.0, Ix=1.0e7, Iy=1.0e7, Zx=1.0e5, Zy=1.0e5, Sx=9.0e4, Sy=9.0e4)
    assert set(_by_name(ff.check(_ctx(sq, L=50000.0, Lb=50000.0), _m(Mrx=1.0e7, Mry=1.0e7)))) == {Y["x"], Y["y"]}
    # section turned (stronger about local y) and bent about local y: LTB cannot be checked with the +y/-y bracing
    turned = _box(H=200.0, B=300.0, Ix=sec.Iy, Iy=sec.Ix, Zx=sec.Zy, Zy=sec.Zx, Sx=sec.Sy, Sy=sec.Sx)
    c = _by_name(ff.check(_ctx(turned, Lb=0.0), _m(Mry=1.0e7)))[ff.NAME_LTB_MINOR]
    assert c.status == "notChecked" and "stronger about local y" in c.reason


# ============================================================================== tier F (a): F8
def test_hand_f8_round_hss():
    # limits: 0.07 E/Fy = 40.58 (compact), 0.31 E/Fy = 179.7 (noncompact), F8 applies below 0.45 E/Fy = 260.9
    assert 0.07 * 200000.0 / 345.0 == pytest.approx(40.5797, abs=1e-4)
    Z, S = 3.0e5, 2.3e5

    def sec(ratio):
        return DesignSection.round_hss(A=5000.0, r=70.0, ratio=ratio, Z=Z, S=S)

    # compact: F8.1 Eq. F8-1: Mn = Fy Z = 103.5e6; demand = resultant of Mx = 3e7, My = 4e7 -> 5e7
    checks = _by_name(ff.check(_ctx(sec(30.0)), _m(Mrx=3.0e7, Mry=4.0e7)))
    assert set(checks) == {ff.NAME_Y_ROUND}
    c = checks[ff.NAME_Y_ROUND]
    assert c.nominal == pytest.approx(345.0 * 3.0e5) and c.demand == pytest.approx(5.0e7)
    assert c.clause == "AISC 360-22 F8.1, Eq. F8-1; B4.1, Table B4.1b; F1(a); B3.1 Eq. B3-1"
    # noncompact, D/t = 100: F8.2(b) Eq. F8-2: Mn = (0.021 E/(D/t) + Fy) S = (42 + 345) * 2.3e5 = 89.01e6
    Mn = (0.021 * 200000.0 / 100.0 + 345.0) * 2.3e5
    assert Mn == pytest.approx(89.01e6)
    checks = _by_name(ff.check(_ctx(sec(100.0)), _m(Mrx=5.0e7)))
    assert checks[ff.NAME_LB_ROUND].nominal == pytest.approx(Mn, rel=1e-14) and "Eq. F8-2" in checks[ff.NAME_LB_ROUND].clause
    assert checks[ff.NAME_Y_ROUND].status == "ok"                           # lower of yielding and local buckling
    # slender walls, D/t = 200: F8.2(c) Eq. F8-3, F8-4: Mn = 0.33 E/(D/t) S = 330 * 2.3e5 = 75.9e6
    c = _by_name(ff.check(_ctx(sec(200.0)), _m(Mrx=5.0e7)))[ff.NAME_LB_ROUND]
    assert c.nominal == pytest.approx(0.33 * 200000.0 / 200.0 * 2.3e5, rel=1e-14) and "Eq. F8-3, F8-4" in c.clause
    # D/t >= 0.45 E/Fy: F8 does not apply -> notChecked
    c = ff.check(_ctx(sec(261.0)), _m(Mrx=5.0e7))
    assert len(c) == 1 and c[0].status == "notChecked" and "0.45 E/Fy" in c[0].reason
    # resultant from stations: the largest sqrt(Mz^2 + My^2) over the stations
    x = np.array([0.0, 1000.0, 2000.0])
    d = _stations(x, mz=[0.0, 3.0e7, 1.0e7], my=[1.0e7, -4.0e7, 0.0])
    c = ff.check(_ctx(sec(30.0)), d)[0]
    assert c.demand == pytest.approx(5.0e7) and c.x == 1000.0
    # the elastic modulus is needed only when the section is not compact
    assert ff.check(_ctx(DesignSection.round_hss(A=5000.0, r=70.0, ratio=30.0, Z=Z)), _m(Mrx=1.0))[0].status == "ok"
    lb = _by_name(ff.check(_ctx(DesignSection.round_hss(A=5000.0, r=70.0, ratio=100.0, Z=Z)), _m(Mrx=1.0)))
    assert lb[ff.NAME_LB_ROUND].status == "notChecked"


# ============================================================================== scope: notChecked, never a guess
def test_out_of_scope_sections_are_not_checked():
    d = _m(Mrx=1.0e8, Mry=1.0e7)
    # noncompact web (h/tw = 100: 90.53 < 100 <= 137.2) -> F4 not implemented; slender web -> F5
    c = _by_name(ff.check(_ctx(_i(web_ratio=100.0), Lb=0.0), d))[Y["x"]]
    assert c.status == "notChecked" and "F4" in c.reason and "noncompact" in c.reason and c.demand == 1.0e8
    c = _by_name(ff.check(_ctx(_i(web_ratio=140.0), Lb=0.0), d))[Y["x"]]
    assert c.status == "notChecked" and "F5" in c.reason
    # the minor axis of those sections is still checked (F6 has no web limit state)
    assert _by_name(ff.check(_ctx(_i(web_ratio=140.0), Lb=0.0), d))[Y["y"]].status == "ok"
    # no section moduli
    checks = _by_name(ff.check(_ctx(_i(Zx=None, Zy=None), Lb=0.0), d))
    assert checks[Y["x"]].status == checks[Y["y"]].status == "notChecked" and "Zx" in checks[Y["x"]].reason
    # an I-shape whose local z axis is not its major axis
    c = _by_name(ff.check(_ctx(_i(Ix=1.0e7, Iy=2.0e7, rx=None, ry=None), Lb=0.0), d))[Y["x"]]
    assert c.status == "notChecked" and "presume x = major axis" in c.reason
    # other shapes: both axes notChecked with the demand kept
    other = DesignSection.unsupported("section 'T' (a tee) is outside the implemented scope", A=5000.0)
    checks = ff.check(_ctx(other), d)
    assert [c.status for c in checks] == ["notChecked", "notChecked"] and [c.demand for c in checks] == [1.0e8, 1.0e7]
    # no Fy
    assert ff.check(_ctx(Fy=None), d)[0].status == "notChecked"
    # a member with a flexural demand that cannot be checked is never "pass"
    res = design.check_member(_ctx(_i(web_ratio=100.0), Lb=0.0), [d])
    assert res.passed is False and res.status == "incomplete"


def test_axis_strength_is_the_lowest_available_strength_including_all_segments():
    L, w = 8000.0, 30.0
    d = _uniform(L, w)
    res = ff.axis_strength(_ctx(L=L, braceTop=(0.5,)), d, "x")
    ltb = _by_name(res.checks)[LTB]
    assert res.Mc == pytest.approx(ltb.capacity) and res.governing == LTB and res.blocked == ""
    # continuously braced: yielding governs
    res = ff.axis_strength(_ctx(L=L, Lb=0.0), d, "x")
    assert res.Mc == pytest.approx(0.9 * 414.0e6) and res.governing == Y["x"]
    # a blocked limit state blocks the axis strength
    res = ff.axis_strength(_ctx(_i(web_ratio=100.0), L=L), d, "x")
    assert res.Mc is None and "F4" in res.blocked
    assert ff.axis_strength(_ctx(), Demand("U1"), "x") is None


# ============================================================================== tier F (b): Design Examples v16.0
def _de_ctx(sec, i, method, L, **kw):
    return DesignContext(member="DE", section=sec, method=method, E=i["E"], G=11200.0, Fy=i["Fy"], Fu=None, L=L,
                         Lx=L, Ly=L, Lz=L, **kw)


def _compact(i):
    """Ratios for an example that ASSERTS compactness without listing bf/2tf and h/tw (User Notes of F2 / F6):
    half the compact limits of Table B4.1b; they only select the F2 / F6 branch."""
    root = math.sqrt(i["E"] / i["Fy"])
    return dict(flange_ratio=0.5 * 0.38 * root, web_ratio=0.5 * 3.76 * root)


def _both(cap_fn, printed, what):
    """LRFD and ASD available strengths (kip-in from the engine) against the printed kip-ft values."""
    assert_printed(cap_fn("LRFD") / FT, printed["phiMn_kipft"], f"{what}: phi Mn")
    assert_printed(cap_fn("ASD") / FT, printed["Mn_over_Omega_kipft"], f"{what}: Mn/Omega")


def test_design_example_f1_1b_w18x50_yielding():
    """Example F.1-1B (F-8, PDF 170): W18x50 continuously braced, Mn = Mp = Fy Zx (Eq. F2-1)."""
    ex = example("F.1-1B")
    i, p = ex["inputs"], ex["printed"]
    sec = DesignSection.i_shape("rolled", A=1.0, Zx=i["Zx"], **_compact(i))        # A: not listed, not used

    def cap(method):
        c = _by_name(ff.check(_de_ctx(sec, i, method, 35 * FT, Lb=0.0), Demand(method, extra={"Mrx": i["Mu_kipft"] * FT})))
        return c[Y["x"]]

    assert_printed(cap("LRFD").nominal, p["Mn_kipin"], "Mn")
    assert_printed(cap("LRFD").nominal / FT, p["Mn_kipft"], "Mn")
    _both(lambda m: cap(m).capacity, p, "F.1-1B")
    assert cap("LRFD").passed and cap("LRFD").ratio == pytest.approx(266.0 / (0.9 * 5050.0 / 12.0), rel=1e-12)


def _w18x50(i):
    return DesignSection.i_shape("rolled", A=1.0, ry=i["ry"], Zx=i["Zx"], Sx=i["Sx"], J=i["J"], rts=i["rts"],
                                 ho=i["ho"], **_compact(i))


def test_design_example_f1_2b_w18x50_inelastic_ltb_with_cb_from_the_diagram():
    """Example F.1-2B (F-10 to F-11, PDF 172-173): braced at the third points; the engine computes Cb of every
    segment from the moment diagram of the uniformly loaded 35 ft span (stations at the tenth points only)."""
    ex = example("F.1-2B")
    i, p = ex["inputs"], ex["printed"]
    L = i["span_ft"] * FT
    sec = _w18x50(i)

    def run(method, w_kip_ft):
        d = _uniform(L, w_kip_ft / FT)
        d = Demand(method, extra=d.extra)
        return _by_name(ff.check(_de_ctx(sec, i, method, L, braceTop=(1 / 3, 2 / 3)), d))

    lrfd, asd = run("LRFD", i["wu_kip_per_ft"]), run("ASD", i["wa_kip_per_ft"])
    c = lrfd[LTB]
    assert_printed(c.demand / FT, i["Mu_kipft"], "Mu from the load")
    assert_printed(c.inputs["Cb"][0], p["Cb_center"], "Cb, centre segment")
    assert_printed(_detail(c, "Lp"), p["Lp_in"], "Lp")
    assert_printed(_detail(c, "Lr"), p["Lr_in"], "Lr")
    assert_printed(c.inputs["Lb"][0] / FT, p["Lb_ft"], "Lb")
    assert c.details["governing segment"][0] == "from 140 to 280" and "Eq. F2-2" in c.clause
    # end segments: Cb = 1.46
    segs = ff.ltb_segments(_de_ctx(sec, i, "LRFD", L, braceTop=(1 / 3, 2 / 3)), ff.moments(_uniform(L, 1.0)))
    assert_printed(segs[0].Cb, p["Cb_end"], "Cb, end segment")
    assert_printed(segs[2].Cb, p["Cb_end"], "Cb, end segment")
    assert_printed(c.nominal, p["Mn_kipin"], "Mn")
    assert_printed(c.nominal / FT, p["Mn_kipft"], "Mn")
    assert_printed(c.capacity / FT, p["phiMn_kipft"], "phi Mn")
    assert_printed(asd[LTB].capacity / FT, p["Mn_over_Omega_kipft"], "Mn/Omega")
    assert_printed(asd[LTB].demand / FT, i["Ma_kipft"], "Ma from the load")
    assert c.passed and asd[LTB].passed
    # yielding: Mp = 5,050 kip-in. (the example refers to F.1-1B)
    assert lrfd[Y["x"]].nominal == pytest.approx(5050.0)


def test_design_example_f1_3b_w18x50_elastic_ltb():
    """Example F.1-3B (F-14 to F-15, PDF 176-177): braced at midspan, Lb = 17.5 ft > Lr: Eq. F2-3 / F2-4."""
    ex = example("F.1-3B")
    i, p = ex["inputs"], ex["printed"]
    L = i["span_ft"] * FT
    sec = _w18x50(i)
    out = {}
    for method, w in (("LRFD", i["wu_kip_per_ft"]), ("ASD", i["wa_kip_per_ft"])):
        d = Demand(method, extra=_uniform(L, w / FT).extra)
        out[method] = _by_name(ff.check(_de_ctx(sec, i, method, L, braceTop=(0.5,)), d))[LTB]
    c = out["LRFD"]
    assert_printed(c.inputs["Cb"][0], p["Cb"], "Cb")
    assert_printed(c.inputs["Lb"][0] / FT, p["Lb_ft"], "Lb")
    assert_printed(_detail(c, "Lp") / FT, p["Lp_ft"], "Lp")
    assert_printed(_detail(c, "Lr") / FT, p["Lr_ft"], "Lr")
    assert_printed(_detail(c, "Fcr"), p["Fcr"], "Fcr")
    assert_printed(c.nominal, p["Mn_kipin"], "Mn")
    assert_printed(c.nominal / FT, p["Mn_kipft"], "Mn")
    assert "Eq. F2-3, F2-4" in c.clause
    _both(lambda m: out[m].capacity, p, "F.1-3B")


def test_design_example_f3b_w21x48_noncompact_flange():
    """Example F.3B (F-24 to F-25, PDF 186-187): continuously braced, flange local buckling Eq. F3-1."""
    ex = example("F.3B")
    i, p = ex["inputs"], ex["printed"]
    web = 0.5 * 3.76 * math.sqrt(i["E"] / i["Fy"])          # the web is not classified in the example (compact)
    sec = DesignSection.i_shape("rolled", A=1.0, Zx=i["Zx"], Sx=i["Sx"], flange_ratio=i["bf_2tf"], web_ratio=web)

    def cap(method):
        d = Demand(method, extra={"Mrx": i["Mu_kipft" if method == "LRFD" else "Ma_kipft"] * FT})
        checks = ff.flexure_major(_de_ctx(sec, i, method, 40 * FT, Lb=0.0), d)
        assert [c.name for c in checks] == [FLB["x"]]         # continuously braced: flange local buckling only
        return checks[0]

    c = cap("LRFD")
    assert_printed(_detail(c, "lambda_p flange"), p["lambda_pf"], "lambda_pf")
    assert_printed(_detail(c, "lambda_r flange"), p["lambda_rf"], "lambda_rf")
    assert_printed(_detail(c, "Mp"), p["Mp_kipin"], "Mp")
    assert_printed(c.nominal, p["Mn_kipin"], "Mn")
    assert_printed(c.nominal / FT, p["Mn_kipft"], "Mn")
    assert "F3.2, Eq. F3-1" in c.clause and c.details["class flange"][0] == "noncompact"
    _both(lambda m: cap(m).capacity, p, "F.3B")
    assert cap("LRFD").passed and cap("ASD").passed


def test_design_example_f5_w12x58_minor_axis():
    """Example F.5 (F-28 to F-29, PDF 190-191): minor-axis yielding, Eq. F6-1 with the 1.6 Fy Sy limit."""
    ex = example("F.5")
    i, p = ex["inputs"], ex["printed"]
    sec = DesignSection.i_shape("rolled", A=1.0, Zy=i["Zy"], Sy=i["Sy"], **_compact(i))

    def cap(method):
        d = Demand(method, extra={"Mry": i["Mu_kipft" if method == "LRFD" else "Ma_kipft"] * FT})
        return _by_name(ff.check(_de_ctx(sec, i, method, 15 * FT), d))[Y["y"]]

    c = cap("LRFD")
    assert_printed(c.nominal, p["Mn_kipin"], "Mn")
    assert_printed(_detail(c, "1.6 Fy Sy"), p["cap_kipin"], "1.6 Fy Sy")
    assert_printed(c.nominal / FT, p["Mn_kipft"], "Mn")
    _both(lambda m: cap(m).capacity, p, "F.5")
    assert c.passed and cap("ASD").passed


def test_design_example_f7b_hss10x6_noncompact_flange_and_ltb():
    """Example F.7B (F-34 to F-36, PDF 196-198): HSS10x6x3/16 braced at the ends, uniform load: flange local
    buckling (Eq. F7-2) governs; LTB (Eq. F7-8, cited as F7-10 in the example) exceeds Mp and is capped."""
    ex = example("F.7B")
    i, p = ex["inputs"], ex["printed"]
    L = i["span_ft"] * FT
    sec = DesignSection.rectangular("rectHSS", A=i["Ag"], Zx=i["Zx"], Sx=i["Sx"], Ix=i["Ix"], ry=i["ry"], J=i["J"],
                                    b_ratio=i["b_t"], h_ratio=i["h_t"])
    assert sec.Ix > sec.Iy                                   # Iy = Ag ry^2 = 34.1 in.^4: strong axis = local z
    out = {}
    for method in ("LRFD", "ASD"):
        d = Demand(method, extra=_uniform(L, i["wu_kip_per_ft"] / FT).extra)
        out[method] = ff.axis_strength(_de_ctx(sec, i, method, L), d, "x")
    checks = _by_name(out["LRFD"].checks)
    assert set(checks) == {Y["x"], FLB["x"], LTB}
    assert_printed(checks[Y["x"]].nominal, p["Mp_kipin"], "Mp")
    assert_printed(checks[FLB["x"]].nominal, p["Mn_flb_kipin"], "Mn (FLB)")
    assert_printed(checks[FLB["x"]].nominal / FT, p["Mn_kipft"], "Mn (FLB)")
    ltb = checks[LTB]
    assert_printed(_detail(ltb, "Lp"), p["Lp_in"], "Lp")
    assert_printed(_detail(ltb, "Lr"), p["Lr_in"], "Lr")
    assert_printed(ltb.inputs["Cb"][0], p["Cb"], "Cb")
    assert_printed(ltb.nominal, p["Mn_ltb_kipin"], "Mn (LTB, capped at Mp)")
    # the uncapped value of the bracket (printed 1,020 kip-in.), by the engine's own formula pieces
    Lp, Lr, Cb, Mp = _detail(ltb, "Lp"), _detail(ltb, "Lr"), ltb.inputs["Cb"][0], 900.0
    assert_printed(Cb * (Mp - (Mp - 0.7 * i["Fy"] * i["Sx"]) * (L - Lp) / (Lr - Lp)), p["Mn_ltb_uncapped_kipin"], "LTB")
    assert out["LRFD"].governing == FLB["x"]                 # flange local buckling controls
    _both(lambda m: out[m].Mc, p, "F.7B")


def test_design_example_f8b_hss8x8_slender_flange():
    """Example F.8B (F-39 to F-41, PDF 201-203): slender flange, be by Eq. F7-4, Se with the ineffective width
    removed from both flanges, Mn = Fy Se (Eq. F7-3)."""
    ex = example("F.8B")
    i, p = ex["inputs"], ex["printed"]
    sec = DesignSection.rectangular("rectHSS", A=1.0, H=i["H"], B=i["B"], t=i["t"], Ix=i["I"], Iy=i["I"],
                                    Zx=i["Z"], Sx=i["S"], Zy=i["Z"], Sy=i["S"], b_ratio=i["b_t"], h_ratio=i["h_t"])

    def cap(method):
        d = Demand(method, extra={"Mrx": i["Mu_kipft" if method == "LRFD" else "Ma_kipft"] * FT})
        return ff.axis_strength(_de_ctx(sec, i, method, 21 * FT), d, "x")

    res = cap("LRFD")
    c = _by_name(res.checks)[FLB["x"]]
    assert_printed(_detail(c, "b"), p["b"], "b = B - 3 t")
    assert_printed(_detail(c, "be"), p["be"], "be")
    assert_printed(_detail(c, "Se"), p["Se"], "Se")
    assert_printed(_detail(c, "Se") * i["H"] / 2.0, p["Ieff"], "Ieff")
    assert_printed(c.nominal, p["Mn_kipin"], "Mn")
    assert_printed(c.nominal / FT, p["Mn_kipft"], "Mn")
    assert "F7.2(c), Eq. F7-3, F7-4" in c.clause and LTB not in _by_name(res.checks)      # square: no LTB
    assert res.governing == FLB["x"]
    _both(lambda m: cap(m).Mc, p, "F.8B")
    assert c.passed


def test_design_example_f9b_pipe8_yielding():
    """Example F.9B (F-43 to F-44, PDF 205-206): Pipe 8 x-Strong, compact, Mn = Fy Z (Eq. F8-1)."""
    ex = example("F.9B")
    i, p = ex["inputs"], ex["printed"]
    sec = DesignSection.round_hss(A=1.0, ratio=i["D_t"], Z=i["Z"])

    def cap(method):
        d = Demand(method, extra={"Mrx": i["Mu_kipft" if method == "LRFD" else "Ma_kipft"] * FT})
        checks = ff.check(_de_ctx(sec, i, method, 16 * FT), d)
        assert [c.name for c in checks] == [ff.NAME_Y_ROUND]
        return checks[0]

    c = cap("LRFD")
    assert_printed(_detail(c, "0.45 E/Fy"), p["limit_F8"], "0.45 E/Fy")
    assert_printed(c.nominal, p["Mn_kipin"], "Mn")
    assert_printed(c.nominal / FT, p["Mn_kipft"], "Mn")
    _both(lambda m: cap(m).capacity, p, "F.9B")
    assert c.passed and cap("ASD").passed
