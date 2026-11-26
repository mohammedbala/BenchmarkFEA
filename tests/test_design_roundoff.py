"""Design decisions must not depend on rounding errors of the analysis (design/demands.py: ZERO_FORCE),
and the scope of Section E4 (torsional buckling) follows its first paragraph (AISC 360-22 p. 16.1-41).

Regression for: LTB checked on the wrong unbraced segment because a moment of -2e-7 N mm at a released
beam end (next to 2e8 N mm in the span) counted as a change of sign.
"""
import numpy as np
import pytest

from benchfea import design
from benchfea.design import DesignContext, DesignSection, Demand
from benchfea.design import demands as dm
from benchfea.design.aisc360_22 import e_compression as ec
from benchfea.design.aisc360_22 import f_flexure as ff
from benchfea.design.aisc360_22 import h_combined as hc
from benchfea.model import build_model
from benchfea.solve import analyze

E, G, FY = 200000.0, 77200.0, 345.0


def _i(fabrication="rolled", **kw):
    base = dict(A=9000.0, rx=160.0, ry=40.0, Zx=1.2e6, Sx=1.05e6, Zy=2.4e5, Sy=1.4e5, J=5.0e5, Cw=4.0e11, rts=48.0,
                ho=380.0, flange_ratio=7.0, web_ratio=30.0)
    base.update(kw)
    return DesignSection.i_shape(fabrication, **base)


def _ctx(sec=None, L=9000.0, **kw):
    base = dict(member="B", section=sec if sec is not None else _i(), method="LRFD", E=E, G=G, Fy=FY, Fu=450.0,
                L=L, Lx=L, Ly=L, Lz=L)
    base.update(kw)
    return DesignContext(**base)


def _beam(L=9000.0, w=20.0, n=13, end_moment=0.0, clean=False):
    """Simply supported beam under a uniform load: Mz = w x (L - x)/2 >= 0 (top flange in compression),
    Vy = -dMz/dx; ``end_moment`` replaces the exact zero at the far end (a released end as the analysis
    returns it)."""
    x = np.linspace(0.0, L, n)
    f = np.zeros((n, 6))
    f[:, 5] = w * x * (L - x) / 2.0
    f[:, 1] = -w * (L - 2.0 * x) / 2.0
    f[-1, 5] = end_moment
    if clean:
        f = dm.clean_stations(x, f)
    return Demand("U1", extra={"x": x, "forces": f})


# ====================================================================================== LTB segment selection
@pytest.mark.parametrize("end_moment", [0.0, -2.0e-7, +2.0e-7, -1.0e-3])
def test_ltb_segments_do_not_flip_on_a_rounding_error_at_a_released_end(end_moment):
    """Top flange braced at the third points, bottom flange unbraced, positive moment throughout: the three
    top-flange segments (Lb = 3000) are checked.  With -2e-7 N mm at the end (Mmax = 2.0e8) the end segment
    used to count as "moment changes sign", and the whole member (Lb = 9000, the segment between points
    where both flanges are braced) was checked instead."""
    ctx = _ctx(braceTop=(1.0 / 3.0, 2.0 / 3.0))
    segs = ff.ltb_segments(ctx, ff.moments(_beam(end_moment=end_moment)))
    assert [round(s.Lb) for s in segs] == [3000, 3000, 3000]
    assert all(s.note == "" for s in segs)
    exact = ff.ltb_segments(ctx, ff.moments(_beam()))
    for a, b in zip(segs, exact):
        assert a.Cb == pytest.approx(b.Cb, rel=1e-9) and a.Mr == pytest.approx(b.Mr, rel=1e-12)
    # the governing LTB ratio is that of the exact diagram, not the three-times-longer segment's
    ltb = {c.name: c for c in ff.check(ctx, _beam(end_moment=end_moment))}[ff.NAME_LTB]
    ref = {c.name: c for c in ff.check(ctx, _beam())}[ff.NAME_LTB]
    assert ltb.ratio == pytest.approx(ref.ratio, rel=1e-9)
    whole = {c.name: c for c in ff.check(_ctx(), _beam())}[ff.NAME_LTB]         # Lb = 9000: what the defect checked
    assert whole.ratio > 1.5 * ref.ratio


def test_a_real_reversal_of_the_moment_still_selects_the_segment_braced_on_both_flanges():
    """The scale does not hide a real change of sign: an end moment of -1e-6 Mmax (1000 times the
    threshold) is a reversal, and the engine rule (segment between points where both flanges are braced)
    applies as before."""
    ctx = _ctx(braceTop=(1.0 / 3.0, 2.0 / 3.0))
    mmax = 20.0 * 9000.0 ** 2 / 8.0
    segs = ff.ltb_segments(ctx, ff.moments(_beam(end_moment=-1.0e-6 * mmax)))
    assert sorted(round(s.Lb) for s in segs) == [3000, 3000, 9000]
    assert any("changes sign" in s.note for s in segs)
    assert dm.ZERO_FORCE == 1.0e-9


def test_a_segment_whose_moment_is_a_rounding_error_is_a_segment_without_moment():
    """Cantilever-like diagram: moment only in the first half; the second half carries 1e-8 N mm of noise.
    No LTB segment (and no blocked free-end segment) is produced for the noise."""
    x = np.array([0.0, 1000.0, 2000.0, 3000.0, 4000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0])   # jump of V at 4000
    f = np.zeros((10, 6))
    f[:5, 5] = 1.0e8 * (1.0 - x[:5] / 4000.0)
    f[:5, 1] = 1.0e8 / 4000.0
    f[5:, 5] = 1.0e-8 * np.array([1.0, -1.0, 1.0, -1.0, 1.0])
    ctx = _ctx(L=8000.0, braceTop=(0.5,), braceBottom=(0.5,), end_free=(False, True))
    segs = ff.ltb_segments(ctx, ff.moments(Demand("U1", extra={"x": x, "forces": f})))
    assert [(s.a, s.b) for s in segs] == [(0.0, 4000.0)] and not segs[0].blocked


# ====================================================================================== zero of a computed force
def test_clean_stations_sets_rounding_errors_to_zero_and_nothing_else():
    x = np.linspace(0.0, 6000.0, 4)
    f = np.array([[-1.0e-9, 5.0e4, 3.0e-12, 0.0, 2.0e-7, -2.0e-7],
                  [4.0e-10, 1.0e4, -1.0e-11, 1.0e-9, -1.0e-8, 9.0e7],
                  [-2.0e-10, -1.0e4, 0.0, 0.0, 6.0e-9, 9.0e7],
                  [1.0e-9, -5.0e4, 2.0e-12, 0.0, 0.2, 0.5]])
    fref, mref = dm.force_scales(x, f)
    assert fref == 5.0e4 and mref == pytest.approx(5.0e4 * 6000.0)       # the largest force x L exceeds Mmax
    g = dm.clean_stations(x, f)
    assert np.all(g[:, 0] == 0.0) and np.all(g[:, 2] == 0.0) and np.all(g[:, 3] == 0.0)
    assert np.array_equal(g[:, 1], f[:, 1])                              # real shears untouched
    assert np.all(g[:, 4] == 0.0)                                        # 0.2 <= 1e-9 x 3e8 = 0.3 N mm
    assert g[0, 5] == 0.0 and g[1, 5] == 9.0e7 and g[3, 5] == 0.5        # 0.5 > 0.3: kept
    assert f[0, 5] == -2.0e-7                                            # the input is not modified
    # a member with moments only still has a force scale (Mmax / L), one with forces only a moment scale
    m_only = np.zeros((2, 6))
    m_only[:, 5] = 6.0e7
    m_only[0, 0] = 1.0e-6                                                # 1e-6 N <= 1e-9 x 6e7 / 6000 = 1e-5 N
    assert dm.clean_stations([0.0, 6000.0], m_only)[0, 0] == 0.0
    assert np.array_equal(dm.clean_stations(x, np.zeros((4, 6))), np.zeros((4, 6)))


def test_demands_from_an_analysis_carry_no_rounding_noise():
    """Simply supported beam with released ends in a plane frame: the analysis returns the end moments, the
    axial force and the out-of-plane moments as rounding errors; the design sees exact zeros, so there is no
    interaction check (no axial force, one moment) and no compression demand."""
    data = {
        "formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": 9000, "y": 0, "z": 0}, "C": {"x": 9000, "y": 0, "z": -4000},
                  "D": {"x": 0, "y": 0, "z": -4000}},
        "supports": {"C": {"type": "fixed"}, "D": {"type": "fixed"}},
        "materials": {"S": {"E": E, "nu": 0.3, "Fy": FY, "Fu": 450.0}},
        "sections": {"I": {"type": "I", "d": 400, "bf": 200, "tf": 14, "tw": 9, "r": 0}},
        "members": {"beam": {"i": "A", "j": "B", "section": "I", "material": "S",
                             "releases": {"i": {"rz": True}, "j": {"rz": True}}},
                    "c1": {"i": "D", "j": "A", "section": "I", "material": "S"},
                    "c2": {"i": "C", "j": "B", "section": "I", "material": "S"}},
        "loadCases": {"D": {"loads": [{"type": "memberDistributed", "member": "beam", "direction": "Z", "w1": -20.0}]}},
        "combinations": {"U1": {"factors": {"D": 1.0}}},
        "design": {"members": {"beam": {"bracing": {"topFlange": [1.0 / 3.0, 2.0 / 3.0]}}}},
    }
    m = build_model(data)
    r = analyze(m)
    raw = np.asarray(r.combinations["U1"].members["beam"].forces)
    d = design.member_demands(m, r, "beam")["U1"]
    f = d.extra["forces"]
    assert np.abs(f[:, 5]).max() == pytest.approx(20.0 * 9000.0 ** 2 / 8.0, rel=1e-9)
    assert f[0, 5] == 0.0 and f[-1, 5] == 0.0 and np.all(f[:, 4] == 0.0) and np.all(f[:, 3] == 0.0)
    assert np.all(np.abs(raw - f) <= 1.0e-9 * np.array([1.0, 1.0, 1.0, 9000.0, 9000.0, 9000.0]) * 9.0e4)
    assert d.Pc == 0.0 or d.Pc > 1.0e-3                       # an axial force is either real or exactly zero
    out = design.check(m, r).members["beam"]
    segs = [c for c in out.checks if c.name == ff.NAME_LTB]
    assert segs and segs[0].inputs["Lb"][0] == pytest.approx(3000.0)


# ====================================================================================== E4 scope
def test_e4_is_evaluated_only_where_its_first_paragraph_or_the_e3_user_note_call_for_it():
    d = Demand("U1", Pc=1.0e5)
    L = 6000.0
    # rolled I, Lz = Lx = Ly: the torsional unbraced length does not exceed the lateral one -> no E4 check,
    # and the flexural buckling checks say why
    assert ec.e4_scope(_ctx(L=L)) [0] is False
    checks = {c.name: c for c in ec.check(_ctx(L=L), d)}
    assert ec.NAME_TB not in checks
    for axis in "xy":
        text = checks[ec.NAME_FB[axis]].details["torsional buckling (E4)"][0]
        assert text.startswith("not evaluated") and "16.1-41" in text
    # Lz exceeds both lateral lengths: required
    on, why = ec.e4_scope(_ctx(L=L, Lz=7000.0))
    assert on and "Lx and Ly" in why
    c = {c.name: c for c in ec.check(_ctx(L=L, Lz=7000.0), d)}[ec.NAME_TB]
    assert c.status == "ok" and "E4-2" in c.clause and any("scope of E4" in n for n in c.notes)
    assert "torsional buckling (E4)" not in {c.name: c for c in ec.check(_ctx(L=L, Lz=7000.0), d)}[ec.NAME_FB["y"]].details
    # Lz exceeds only the smaller lateral length: the text does not say which length is meant -> evaluated
    on, why = ec.e4_scope(_ctx(L=L, Ly=3000.0, Lz=4000.0))
    assert on and "does not say which" in why
    # Lz below both lateral lengths
    assert ec.e4_scope(_ctx(L=L, Ly=3000.0, Lz=3000.0))[0] is False
    assert ec.e4_scope(_ctx(L=L, Ly=3000.0, Lz=2000.0))[0] is False
    # equal lengths that differ in the last digit (brace fractions) are equal
    assert ec.e4_scope(_ctx(L=L, Lz=L * (1.0 + 1.0e-12)))[0] is False
    # effective lengths (User Note of E3): Kz Lz > Ky Ly -> evaluated
    on, why = ec.e4_scope(_ctx(L=L, Kz=1.2))
    assert on and "EFFECTIVE" in why
    assert ec.e4_scope(_ctx(L=L, Kz=0.8))[0] is False
    assert ec.e4_scope(_ctx(L=L, Kx=2.0, Ky=2.0, Kz=1.5))[0] is False
    # built-up I-shape and box: "certain doubly symmetric members, such as ... built-up members" -> evaluated
    assert ec.e4_scope(_ctx(_i("builtUp"), L=L))[0] is True
    box = DesignSection.rectangular("box", H=300.0, B=200.0, t=10.0, A=9600.0, Ix=1.2e8, Iy=6.4e7, J=1.3e8)
    assert ec.e4_scope(_ctx(box, L=L))[0] is True
    hss = DesignSection.rectangular("rectHSS", H=300.0, B=200.0, t=10.0, A=9600.0, Ix=1.2e8, Iy=6.4e7, J=1.3e8)
    assert ec.e4_scope(_ctx(hss, L=L))[0] is False and ec.e4_scope(_ctx(hss, L=L, Lz=6500.0))[0] is True
    # a missing length: the scope cannot be decided -> the check is produced (and reports what is missing)
    assert ec.e4_scope(_ctx(L=L, Lz=None))[0] is True
    c = ec.torsional_buckling(_ctx(L=L, Lz=None), d)[0]
    assert c.status == "notChecked"


def test_e4_out_of_scope_does_not_enter_the_axial_strength_of_the_interaction():
    """H1: Pc is the lowest of the limit states that apply.  With Lz = Lx = Ly torsional buckling does not
    apply, so a section without Cw no longer blocks the interaction."""
    sec = _i(Cw=None)
    d = Demand("U1", Pc=1.0e5, extra={"Mrx": 1.0e8, "Mry": 0.0})
    c, = hc.interaction(_ctx(sec, L=6000.0, Lb=0.0), d)
    assert c.status == "ok" and c.inputs["Pc"][2] == "derived: " + ec.NAME_FB["y"]
    c, = hc.interaction(_ctx(sec, L=6000.0, Lz=9000.0, Lb=0.0), d)
    assert c.status == "notChecked" and "torsional buckling" in c.reason
