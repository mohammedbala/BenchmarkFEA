"""Single members through the whole pipeline (file -> model -> analyze) against the independent ODE
solution: every support arrangement x every member load type x both bending planes, shear deformation
on/off; axial and torsional loads; global, local and projected load directions; self weight; a skew
member with roll; exactness of the stations including the jumps at point loads."""
import math

import numpy as np
import pytest

from benchfea.model import build_model
from benchfea.solve import analyze

from solve_reference import (EA, EIy, EIz, GAy, GAz, GJ, SEC, bar_bvp, beam_bvp, clamped, free, make, pinned)

L = 5000.0
REL = 1e-9          # the reference (Chebyshev collocation) and the engine are both exact to round-off

SUPPORTS = {
    # name: (support a, support b, bc at 0, bc at L)
    "cantilever": ("fixed", None, clamped, free),
    "simply supported": ({"type": "pinned", "dofs": {"rx": "fixed"}}, {"dofs": {"uy": "fixed", "uz": "fixed"}}, pinned, pinned),
    "fixed-fixed": ("fixed", "fixed", clamped, clamped),
    "propped": ("fixed", {"dofs": {"uy": "fixed", "uz": "fixed"}}, clamped, pinned),
}

# name: (file loads for local direction d, reference line loads, reference point loads)
def LOADS(d, mdir):
    return {
        "uniform full": ([dist(d, -4.0)], [(0, L, -4.0, -4.0)], []),
        "uniform partial": ([dist(d, 6.0, start=0.25, end=0.6)], [(0.25 * L, 0.6 * L, 6.0, 6.0)], []),
        "trapezoid partial": ([dist(d, -2.0, w2=-9.0, start=0.1, end=0.85)], [(0.1 * L, 0.85 * L, -2.0, -9.0)], []),
        "triangle full": ([dist(d, 0.0, w2=-8.0)], [(0, L, 0.0, -8.0)], []),
        "sign-changing": ([dist(d, 5.0, w2=-5.0, start=0.2, end=0.8)], [(0.2 * L, 0.8 * L, 5.0, -5.0)], []),
        "point force": ([point(d, -12000.0, 0.35)], [], [(0.35 * L, -12000.0, 0.0)]),
        "two point forces": ([point(d, -7000.0, 0.2), point(d, 3000.0, 0.7)], [], [(0.2 * L, -7000.0, 0.0), (0.7 * L, 3000.0, 0.0)]),
        "point moment": ([point(mdir, 8.0e6, 0.45, kind="moment")], [], [(0.45 * L, 0.0, 8.0e6)]),
        "mixed": ([dist(d, -3.0, w2=-1.0, start=0.0, end=0.5), point(d, -5000.0, 0.5), point(mdir, -2.0e6, 0.8, kind="moment"),
                   dist(d, 2.0, start=0.5, end=1.0)],
                  [(0, 0.5 * L, -3.0, -1.0), (0.5 * L, L, 2.0, 2.0)], [(0.5 * L, -5000.0, 0.0), (0.8 * L, 0.0, -2.0e6)]),
    }


def dist(direction, w1, w2=None, start=None, end=None, **kw):
    ld = {"type": "memberDistributed", "member": "m", "direction": direction, "w1": w1, **kw}
    if w2 is not None:
        ld["w2"] = w2
    if start is not None:
        ld["start"] = start
    if end is not None:
        ld["end"] = end
    return ld


def point(direction, value, position, kind="force"):
    return {"type": "memberPoint", "member": "m", "direction": direction, "value": value, "position": position,
            "kind": kind}


def solve(doc):
    model = build_model(doc)
    res = analyze(model)
    return model, res


def sides(x):
    """Side of each station: the first of two equal x is the left limit, the second the right limit;
    x = 0 is taken after the start (the reference has no state before it) and x = L before the end."""
    out = []
    for k, v in enumerate(x):
        if k + 1 < len(x) and x[k + 1] == v:
            out.append(-1)
        elif k > 0 and x[k - 1] == v:
            out.append(+1)
        else:
            out.append(+1 if v == 0.0 else -1)
    return out


def beam_doc(support, loads, shear=True, stations=11):
    sa, sb = support[0], support[1]
    sup = {"a": sa}
    if sb is not None:
        sup["b"] = sb
    return make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, sup, {"c": loads},
                analysis={"shearDeformation": shear, "resultStations": stations})


@pytest.mark.parametrize("support", list(SUPPORTS))
@pytest.mark.parametrize("load", list(LOADS("y", "z")))
@pytest.mark.parametrize("plane,shear", [("xy", True), ("xy", False), ("xz", True)])
def test_beam_against_ode_solution(support, load, plane, shear):
    """Beam along +X, Z up: local y = +Z, local z = -Y.  Loads in local axes."""
    d, mdir = ("y", "z") if plane == "xy" else ("z", "y")
    file_loads, ref_dist, ref_points = LOADS(d, mdir)[load]
    sup = SUPPORTS[support]
    model, res = solve(beam_doc(sup, file_loads, shear))
    EI, GA = (EIz, GAy) if plane == "xy" else (EIy, GAz)
    ref = beam_bvp(L, EI, GA if shear else None, sup[2], sup[3], ref_dist, ref_points, plane=plane)
    c = res.combinations["c"]
    assert c.invariants.passed
    mr = c.members["m"]
    iV, iM, iD = (1, 5, 1) if plane == "xy" else (2, 4, 2)
    xs = np.linspace(0, L, 201)
    vals = np.array([ref(x) for x in xs])
    sd, sth, sM, sV = (max(np.abs(vals[:, k]).max(), 1e-300) for k in range(4))
    sV = max(sV, sM / L)
    for k, (x, side) in enumerate(zip(mr.x, sides(mr.x))):
        dref, th, M, V = ref(x, side)
        assert abs(mr.forces[k, iV] - V) <= REL * sV, ("V", x, side)
        assert abs(mr.forces[k, iM] - M) <= REL * sM, ("M", x, side)
        assert abs(mr.displacements[k, iD] - dref) <= REL * sd, ("d", x, side)
    # the other plane, axial force and torsion stay exactly zero
    other = [0, 3] + ([2, 4] if plane == "xy" else [1, 5])
    assert np.abs(mr.forces[:, other]).max() <= 1e-9 * sV
    # end forces: at i minus the resultants just inside (no end loads here), at j the resultants
    a, b = ref(0.0, 1), ref(L)
    assert abs(mr.end_forces[iV] + a[3]) <= REL * sV and abs(mr.end_forces[iM] + a[2]) <= REL * sM
    assert abs(mr.end_forces[6 + iV] - b[3]) <= REL * sV and abs(mr.end_forces[6 + iM] - b[2]) <= REL * sM
    np.testing.assert_allclose(mr.forces[-1], mr.end_forces[6:], rtol=1e-9, atol=1e-9 * sM)
    np.testing.assert_allclose(mr.forces[0], -mr.end_forces[:6], rtol=0, atol=0)
    # node displacement and rotation at b (global): uz = +v, ry = -theta_z   |   uy = -w, rz = +theta_y (local y = +Z)
    ub = c.displacements[1]
    if plane == "xy":
        assert abs(ub[2] - b[0]) <= REL * sd and abs(-ub[4] - b[1]) <= REL * sth
    else:
        assert abs(-ub[1] - b[0]) <= REL * sd and abs(ub[5] - b[1]) <= REL * sth
    # reactions: forces by the supports on the structure = minus (end forces on the member)... i.e. the
    # joint gives the member its end force and receives the reaction
    R = c.reactions
    if plane == "xy":
        assert abs(R[0, 2] - mr.end_forces[1]) <= 1e-9 * sV
    else:
        assert abs(-R[0, 1] - mr.end_forces[2]) <= 1e-9 * sV


def test_stations_report_both_sides_of_a_point_load_and_the_kinks_of_line_loads():
    loads = [point("y", -6000.0, 0.3), point("y", -6000.0, 0.5), point("z", 4.0e6, 0.5, kind="moment"),
             dist("y", -3.0, start=0.12, end=0.66)]
    model, res = solve(beam_doc(SUPPORTS["simply supported"], loads, stations=11))
    mr = res.combinations["c"].members["m"]
    x = list(mr.x)
    assert x.count(0.3 * L) == 2 and x.count(0.5 * L) == 2          # both sides (0.5 L is also a tenth point)
    assert x.count(0.12 * L) == 1 and x.count(0.66 * L) == 1        # ends of the line load
    assert x == sorted(x) and x[0] == 0.0 and x[-1] == L and len(x) == 11 + 1 + 1 + 2 + 0
    k = x.index(0.3 * L)
    assert mr.forces[k + 1, 1] - mr.forces[k, 1] == pytest.approx(6000.0, rel=1e-12)       # Vy' = -w: drop by P
    assert mr.forces[k + 1, 5] == pytest.approx(mr.forces[k, 5], rel=1e-12)
    k = x.index(0.5 * L)
    assert mr.forces[k + 1, 1] - mr.forces[k, 1] == pytest.approx(6000.0, rel=1e-12)
    assert mr.forces[k + 1, 5] - mr.forces[k, 5] == pytest.approx(-4.0e6, rel=1e-12)
    np.testing.assert_allclose(mr.displacements[k + 1], mr.displacements[k], rtol=1e-13)   # no jump in deflection


def test_station_layout_is_the_same_in_every_combination():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed"},
               {"P": [point("y", -1.0, 0.37)], "Q": [dist("y", -1.0, start=0.21, end=0.5)]},
               combos={"onlyP": {"factors": {"P": 1.0}}, "onlyQ": {"factors": {"Q": 2.0}}, "both": {"factors": {"P": 1, "Q": 1}}})
    model, res = solve(doc)
    xs = [res.combinations[c].members["m"].x for c in ("onlyP", "onlyQ", "both")]
    np.testing.assert_array_equal(xs[0], xs[1])
    np.testing.assert_array_equal(xs[0], xs[2])
    assert list(xs[0]).count(0.37 * L) == 2
    mq = res.combinations["onlyQ"].members["m"]
    k = list(mq.x).index(0.37 * L)
    np.testing.assert_array_equal(mq.forces[k], mq.forces[k + 1])   # no jump where this combination has no load


def test_point_loads_exactly_at_the_member_ends_go_to_the_joints():
    """A force at x = 0 of a cantilever goes straight into the support; a force at x = L is the tip load."""
    doc = beam_doc(SUPPORTS["cantilever"], [point("y", -500.0, 0.0), point("y", -900.0, 1.0)])
    model, res = solve(doc)
    c = res.combinations["c"]
    mr = c.members["m"]
    ref = beam_bvp(L, EIz, GAy, clamped, lambda s: [s[2], s[3] + 900.0])       # V(L-) = -900 on the +x face... = load
    assert c.displacements[1, 2] == pytest.approx(ref(L)[0], rel=REL)
    assert c.reactions[0, 2] == pytest.approx(1400.0, rel=1e-12)
    assert mr.end_forces[1] == pytest.approx(1400.0, rel=1e-12)             # the joint holds both loads
    assert mr.forces[0, 1] == pytest.approx(-1400.0, rel=1e-12)             # station at 0 = minus the end force
    assert mr.forces[1, 1] == pytest.approx(-900.0, rel=1e-12)
    assert abs(mr.forces[-1, 1]) <= 1e-9 * 900 and abs(mr.end_forces[7]) <= 1e-9 * 900    # free end


@pytest.mark.parametrize("support", ["fixed-free", "fixed-fixed"])
def test_axial_and_torsional_member_loads_against_ode_solution(support):
    loads = [dist("x", 3.0, w2=1.0, start=0.1, end=0.7), point("x", -4000.0, 0.4), point("x", 6.0e5, 0.55, kind="moment")]
    sup = ("fixed", None) if support == "fixed-free" else ("fixed", "fixed")
    model, res = solve(beam_doc(sup, loads))
    mr = res.combinations["c"].members["m"]
    end = (lambda s: [s[1]]) if support == "fixed-free" else (lambda s: [s[0]])
    ax = bar_bvp(L, EA, lambda s: [s[0]], end, [(0.1 * L, 0.7 * L, 3.0, 1.0)], [(0.4 * L, -4000.0)])
    to = bar_bvp(L, GJ, lambda s: [s[0]], end, (), [(0.55 * L, 6.0e5)])
    sN = max(abs(ax(x, s)[1]) for x in np.linspace(0, L, 21) for s in (-1, 1))
    su = max(abs(ax(x)[0]) for x in np.linspace(0, L, 21))
    for k, (x, side) in enumerate(zip(mr.x, sides(mr.x))):
        assert abs(mr.forces[k, 0] - ax(x, side)[1]) <= REL * sN
        assert abs(mr.displacements[k, 0] - ax(x, side)[0]) <= REL * su
        assert abs(mr.forces[k, 3] - to(x, side)[1]) <= REL * 6.0e5
    if support == "fixed-free":
        c = res.combinations["c"]
        assert c.displacements[1, 0] == pytest.approx(ax(L)[0], rel=REL)
        assert c.displacements[1, 3] == pytest.approx(to(L)[0], rel=REL)


# ----------------------------------------------------- inclined member: global / projected / local
ALPHA = math.atan2(3.0, 4.0)
LI = 5000.0


def inclined(loads, supports=None):
    return make({"a": (0, 0, 0), "b": (4000, 0, 3000)}, {"m": ("a", "b")},
                supports or {"a": {"type": "pinned", "dofs": {"rx": "fixed"}}, "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}},
                {"c": loads})


@pytest.mark.parametrize("projected", [False, True])
def test_global_and_projected_gravity_load_on_an_inclined_member(projected):
    """Member rises 3 in 4 (cos = 0.8, sin = 0.6).  A vertical line load q per unit MEMBER length has local
    components (-q sin, -q cos) along (x, y); per unit HORIZONTAL length it is multiplied by cos."""
    q = 5.0
    ld = {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q}
    if projected:
        ld["projected"] = True
    # roller at b frees ux: support there holds uz only -> vertical reaction; the member carries axial force
    model, res = solve(inclined([ld]))
    f = 0.8 if projected else 1.0
    wx, wy = -q * f * 0.6, -q * f * 0.8
    c = res.combinations["c"]
    mr = c.members["m"]
    total = q * f * LI
    assert c.reactions[:, 2].sum() == pytest.approx(total, rel=1e-12)
    assert c.reactions[0, 2] == pytest.approx(total / 2, rel=1e-10) and abs(c.reactions[0, 0]) <= 1e-9 * total
    # local problem: end b can slide horizontally: reaction at b is vertical = (Rb sin, Rb cos) in (x, y)
    Rb = total / 2
    ref = beam_bvp(LI, EIz, GAy, pinned, pinned, [(0, LI, wy, wy)])
    axial = bar_bvp(LI, EA, lambda s: [s[0]], lambda s: [s[1] - Rb * 0.6], [(0, LI, wx, wx)])
    sM = abs(ref(LI / 2)[2])
    for k, x in enumerate(mr.x):
        assert abs(mr.forces[k, 5] - ref(x)[2]) <= REL * sM
        assert abs(mr.forces[k, 0] - axial(x)[1]) <= REL * total
    # chord-relative deflection equals the pinned-pinned ODE deflection
    chord = mr.displacements[0, 1] + (mr.displacements[-1, 1] - mr.displacements[0, 1]) * mr.x / LI
    for k, x in enumerate(mr.x):
        assert abs(mr.displacements[k, 1] - chord[k] - ref(x)[0]) <= REL * abs(ref(LI / 2)[0])


def test_global_direction_point_load_and_local_direction_agree_with_hand_resolution():
    P = 10000.0
    g = solve(inclined([{"type": "memberPoint", "member": "m", "direction": "X", "value": P, "position": 0.4}]))[1]
    l = solve(inclined([{"type": "memberPoint", "member": "m", "direction": "x", "value": P * 0.8, "position": 0.4},
                        {"type": "memberPoint", "member": "m", "direction": "y", "value": -P * 0.6, "position": 0.4}]))[1]
    a, b = g.combinations["c"], l.combinations["c"]
    np.testing.assert_allclose(a.displacements, b.displacements, rtol=1e-11, atol=1e-14)
    np.testing.assert_allclose(a.members["m"].forces, b.members["m"].forces, rtol=1e-11, atol=1e-7)
    assert a.reactions[:, 0].sum() == pytest.approx(-P, rel=1e-12)


def test_self_weight_is_the_line_load_gamma_A_along_gravity():
    gamma, A = 7.7e-5, SEC["A"]
    sw = solve(inclined([{"type": "selfWeight"}]))[1].combinations["c"]
    eq = solve(inclined([{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -gamma * A}]))[1].combinations["c"]
    np.testing.assert_allclose(sw.displacements, eq.displacements, rtol=1e-12, atol=1e-18)
    np.testing.assert_allclose(sw.members["m"].forces, eq.members["m"].forces, rtol=1e-12, atol=1e-9)
    assert sw.reactions[:, 2].sum() == pytest.approx(gamma * A * LI, rel=1e-12)
    # factor and direction
    s2 = solve(inclined([{"type": "selfWeight", "factor": 1.5, "direction": [1, 0, 0]}]))[1].combinations["c"]
    assert s2.reactions[:, 0].sum() == pytest.approx(-1.5 * gamma * A * LI, rel=1e-12)


def test_self_weight_acts_on_the_flexible_length_only():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": {"i": "a", "j": "b", "offsets": {"i": [400, 0, 0], "j": [-600, 0, 0]}}},
               {"a": "fixed"}, {"c": [{"type": "selfWeight"}]})
    c = solve(doc)[1].combinations["c"]
    assert c.reactions[0, 2] == pytest.approx(7.7e-5 * SEC["A"] * (L - 1000.0), rel=1e-12)


# ------------------------------------------------------------------ skew cantilever with roll
@pytest.mark.parametrize("roll", [0.0, 37.0, -110.0])
@pytest.mark.parametrize("shear", [True, False])
def test_skew_cantilever_with_roll_tip_loads_and_member_loads(roll, shear):
    """Arbitrary direction in space.  The load is resolved in local axes with the member triad (tested in
    test_model), each plane is solved by the ODE reference, and the tip motion is rotated back."""
    pj = np.array([2200.0, -1400.0, 3100.0])
    Lm = float(np.linalg.norm(pj))
    F = np.array([900.0, -1500.0, 700.0])
    M = np.array([2.0e5, 3.0e5, -1.0e5])
    qg = np.array([0.4, 0.2, -1.1])                    # global line load, two global components via two loads + Z
    doc = make({"a": (0, 0, 0), "b": tuple(pj)}, {"m": {"i": "a", "j": "b", "orientation": {"roll": roll}}}, {"a": "fixed"},
               {"c": [{"type": "nodal", "node": "b", "Fx": F[0], "Fy": F[1], "Fz": F[2], "Mx": M[0], "My": M[1], "Mz": M[2]},
                      {"type": "memberDistributed", "member": "m", "direction": "X", "w1": qg[0]},
                      {"type": "memberDistributed", "member": "m", "direction": "Y", "w1": qg[1]},
                      {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": qg[2], "start": 0.3, "end": 0.9}]},
               analysis={"shearDeformation": shear})
    model, res = solve(doc)
    R = model.members["m"].R
    Fl, Ml, ql = R @ F, R @ M, R @ np.diag(qg)          # columns of ql: local vectors of the three global loads
    GAY, GAZ = (GAy, GAz) if shear else (None, None)
    full = lambda c: [(0.0, Lm, ql[c, 0] + ql[c, 1], ql[c, 0] + ql[c, 1]), (0.3 * Lm, 0.9 * Lm, ql[c, 2], ql[c, 2])]
    xy = beam_bvp(Lm, EIz, GAY, clamped, lambda s: [s[2] - Ml[2], s[3] - Fl[1]], full(1))
    xz = beam_bvp(Lm, EIy, GAZ, clamped, lambda s: [s[2] - Ml[1], s[3] - Fl[2]], full(2), plane="xz")
    ax = bar_bvp(Lm, EA, lambda s: [s[0]], lambda s: [s[1] - Fl[0]], full(0))
    to = bar_bvp(Lm, GJ, lambda s: [s[0]], lambda s: [s[1] - Ml[0]])
    u_local = np.array([ax(Lm)[0], xy(Lm)[0], xz(Lm)[0]])
    th_local = np.array([to(Lm)[0], xz(Lm)[1], xy(Lm)[1]])
    c = res.combinations["c"]
    np.testing.assert_allclose(c.displacements[1, :3], R.T @ u_local, rtol=0, atol=REL * np.abs(u_local).max())
    np.testing.assert_allclose(c.displacements[1, 3:], R.T @ th_local, rtol=0, atol=REL * np.abs(th_local).max())
    mr = c.members["m"]
    for k, x in enumerate(mr.x):
        ref = np.array([ax(x)[1], xy(x)[3], xz(x)[3], to(x)[1], xz(x)[2], xy(x)[2]])
        scale = np.array([np.abs(Fl).max() + np.abs(qg).max() * Lm] * 3 + [abs(xy(0.0, 1)[2]) + abs(xz(0.0, 1)[2])] * 3)
        assert (np.abs(mr.forces[k] - ref) <= REL * scale).all(), x
        assert abs(mr.displacements[k, 1] - xy(x)[0]) <= REL * np.abs(u_local).max()
        assert abs(mr.displacements[k, 2] - xz(x)[0]) <= REL * np.abs(u_local).max()
    assert c.invariants.passed


def test_vertical_axis_y_gives_the_same_physics_as_z():
    """The same cantilever described in a Y-up and a Z-up file."""
    ld = lambda d: [{"type": "memberDistributed", "member": "m", "direction": d, "w1": -3.0}, {"type": "selfWeight"}]
    z = solve(make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed"}, {"c": ld("Z")}))[1].combinations["c"]
    y = solve(make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed"}, {"c": ld("Y")},
                   vertical="Y"))[1].combinations["c"]
    assert z.displacements[1, 2] == pytest.approx(y.displacements[1, 1], rel=1e-12)
    np.testing.assert_allclose(z.members["m"].forces, y.members["m"].forces, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(z.members["m"].displacements, y.members["m"].displacements, rtol=1e-12, atol=1e-15)
