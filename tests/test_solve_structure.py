"""Structures: supports, springs, prescribed displacements, releases in a frame, rigid offsets, rigid
links, diaphragms, trusses, combinations and envelopes -- and the equivalences that tie them together
(split members, release limits, offset vs. stiff stub, constraint vs. stiff members, rigid rotation of
the whole model)."""
import numpy as np
import pytest

from benchfea.model import build_model
from benchfea.problems import ModelUnstable
from benchfea.solve import ModelRefused, analyze

from solve_reference import (E, EA, EIy, EIz, G, GAy, GAz, GJ, SEC, bar_bvp, beam_bvp, clamped, make, pinned)

L = 4000.0
REL = 1e-9


def run(doc, **kw):
    model = build_model(doc)
    return model, analyze(model, **kw)


def nodal(node, **kw):
    return {"type": "nodal", "node": node, **kw}


def dist(member, direction, w1, **kw):
    return {"type": "memberDistributed", "member": member, "direction": direction, "w1": w1, **kw}


def point(member, direction, value, position, **kw):
    return {"type": "memberPoint", "member": member, "direction": direction, "value": value, "position": position, **kw}


# ------------------------------------------------------------------------------- supports
@pytest.mark.parametrize("kfac", [1e-6, 1e-2, 1.0, 1e3, 1e9])
def test_translational_spring_support_against_ode_solution(kfac):
    k = kfac * EIz / L ** 3
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed", "b": {"dofs": {"uz": {"spring": k}}}},
               {"c": [dist("m", "Z", -5.0)]})
    model, res = run(doc)
    c = res.combinations["c"]
    # the spring pushes back on the member end: end force on the member at j, Fy = V(L) = -k v(L)
    ref = beam_bvp(L, EIz, GAy, clamped, lambda s: [s[2], s[3] + k * s[0]], [(0, L, -5.0, -5.0)])
    assert c.displacements[1, 2] == pytest.approx(ref(L)[0], rel=REL)
    assert c.reactions[1, 2] == pytest.approx(-k * c.displacements[1, 2], rel=1e-12)      # reported as -k u
    assert c.reactions[1, 2] == pytest.approx(ref(L)[3], rel=REL)          # the joint passes the reaction to the member end
    assert c.reactions[:, 2].sum() == pytest.approx(5.0 * L, rel=1e-10)


@pytest.mark.parametrize("kfac", [1e-4, 0.5, 20.0, 1e6])
def test_rotational_spring_support_against_ode_solution(kfac):
    k = kfac * EIz / L
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")},
               {"a": {"type": "pinned", "dofs": {"rx": "fixed", "ry": {"spring": k}}}, "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}},
               {"c": [dist("m", "Z", -5.0)]})
    c = run(doc)[1].combinations["c"]
    # at x = 0 the spring gives the member the end moment Mz_i = -k theta, and Mz_i = -M(0)
    ref = beam_bvp(L, EIz, GAy, lambda s: [s[0], s[2] - k * s[1]], pinned, [(0, L, -5.0, -5.0)])
    theta = -c.displacements[0, 4]                       # local theta_z = -global ry for a member along +X
    assert theta == pytest.approx(ref(0.0, 1)[1], rel=REL)
    assert c.reactions[0, 4] == pytest.approx(-k * c.displacements[0, 4], rel=1e-12)


def test_spring_support_limits_reach_free_and_fixed():
    def tip(support_b):
        sup = {"a": "fixed"}
        if support_b:
            sup["b"] = support_b
        return run(make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, sup, {"c": [dist("m", "Z", -5.0)]}))[1].combinations["c"]
    free_tip, fixed_tip = tip(None), tip({"dofs": {"uz": "fixed"}})
    soft, hard = tip({"dofs": {"uz": {"spring": 1e-9}}}), tip({"dofs": {"uz": {"spring": 1e15}}})
    assert soft.displacements[1, 2] == pytest.approx(free_tip.displacements[1, 2], rel=1e-9)
    assert hard.reactions[1, 2] == pytest.approx(fixed_tip.reactions[1, 2], rel=1e-7)
    assert abs(hard.displacements[1, 2]) < 1e-9


def test_prescribed_support_displacement_against_ode_solution():
    """Settlement and rotation imposed at the right end of a fixed-fixed beam, combined with a load."""
    dz, ry = -12.0, 0.004
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed", "b": "fixed"},
               {"S": [{"type": "supportDisplacement", "node": "b", "uz": dz, "ry": ry}], "Q": [dist("m", "Z", -5.0)]},
               combos={"c": {"factors": {"S": 1.0, "Q": 1.0}}, "half": {"factors": {"S": 0.5}}})
    model, res = run(doc)
    c = res.combinations["c"]
    np.testing.assert_allclose(c.displacements[1], [0, 0, dz, 0, ry, 0], atol=1e-15)
    ref = beam_bvp(L, EIz, GAy, clamped, lambda s: [s[0] - dz, s[1] + ry], [(0, L, -5.0, -5.0)])    # theta_z = -ry
    mr = c.members["m"]
    sM = max(abs(ref(x)[2]) for x in np.linspace(0, L, 21))
    for k, x in enumerate(mr.x):
        assert abs(mr.forces[k, 5] - ref(x)[2]) <= REL * sM
        assert abs(mr.displacements[k, 1] - ref(x)[0]) <= REL * abs(dz)
    assert c.reactions[:, 2].sum() == pytest.approx(5.0 * L, rel=1e-9)
    h = res.combinations["half"]
    np.testing.assert_allclose(h.displacements[1], [0, 0, dz / 2, 0, ry / 2, 0], atol=1e-15)
    assert abs(h.reactions[:, 2].sum()) <= 1e-9 * abs(h.reactions[0, 2])          # self-equilibrated
    assert c.invariants.passed and h.invariants.passed


def test_support_types_fixed_pinned_roller_and_per_dof():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (2 * L, 0, 0)}, {"m1": ("a", "b"), "m2": ("b", "c")},
               {"a": {"type": "pinned", "dofs": {"rx": "fixed"}}, "b": "roller", "c": {"type": "roller", "dofs": {"uy": "fixed"}}},
               {"c": [dist("m1", "Z", -4.0), dist("m2", "Z", -4.0), nodal("b", Fy=100.0)]})
    c = run(doc)[1].combinations["c"]
    R = c.reactions
    assert R[1, 0] == R[1, 1] == 0.0 and R[1, 3:].tolist() == [0, 0, 0]        # roller: vertical only
    assert R[:, 2].sum() == pytest.approx(8.0 * L, rel=1e-10)
    assert R[1, 2] > R[0, 2] and R[0, 2] == pytest.approx(R[2, 2], rel=1e-9)   # symmetric two-span beam
    assert R[:, 1].sum() == pytest.approx(-100.0, rel=1e-9)
    assert c.displacements[1, 2] == 0.0 and c.displacements[1, 1] != 0.0


# ------------------------------------------------------------------------- frame releases
def portal(releases=None, **kw):
    mem = {"c1": ("a", "b"), "bm": {"i": "b", "j": "c", **({"releases": releases} if releases else {})}, "c2": ("d", "c")}
    return make({"a": (0, 0, 0), "b": (0, 0, 3000), "c": (5000, 0, 3000), "d": (5000, 0, 0)}, mem,
                {"a": "fixed", "d": "fixed"},
                {"c": [dist("bm", "Z", -6.0), nodal("b", Fx=9000.0), point("bm", "Z", -4000.0, 0.3)]}, **kw)


def test_full_moment_releases_give_zero_end_moments_and_keep_equilibrium():
    model, res = run(portal({"i": {"rz": True}, "j": {"rz": True}}))
    c = res.combinations["c"]
    ef = c.members["bm"].end_forces
    assert abs(ef[5]) <= 1e-12 * 1e7 and abs(ef[11]) <= 1e-12 * 1e7
    assert c.invariants.passed
    assert c.reactions[:, 2].sum() == pytest.approx(6.0 * 5000 + 4000, rel=1e-10)
    assert c.reactions[:, 0].sum() == pytest.approx(-9000.0, rel=1e-10)
    # the released beam is simply supported between its (moving) ends: chord-relative deflection = ODE solution
    mr = c.members["bm"]
    ref = beam_bvp(5000.0, EIz, GAy, pinned, pinned, [(0, 5000.0, -6.0, -6.0)], [(1500.0, -4000.0, 0.0)])
    chord = mr.displacements[0, 1] + (mr.displacements[-1, 1] - mr.displacements[0, 1]) * mr.x / 5000.0
    side = [-1 if (k + 1 < len(mr.x) and mr.x[k + 1] == x) else 1 for k, x in enumerate(mr.x)]
    for k, x in enumerate(mr.x):
        assert abs(mr.displacements[k, 1] - chord[k] - ref(x)[0]) <= REL * abs(ref(2500.0)[0])
        assert abs(mr.forces[k, 5] - ref(x, side[k] if 0 < x < 5000 else (1 if x == 0 else -1))[2]) <= REL * abs(ref(2500.0)[2])


def test_partial_release_tends_to_the_rigid_and_to_the_pinned_solution_and_is_monotonic():
    rigid = run(portal())[1].combinations["c"]
    hinged = run(portal({"i": {"rz": True}, "j": {"rz": True}}))[1].combinations["c"]
    kb = EIz / 5000.0
    prev = None
    for kfac in [1e-12, 1e-6, 1e-3, 0.1, 1.0, 10.0, 1e3, 1e6, 1e12]:
        k = kfac * kb
        c = run(portal({"i": {"rz": {"spring": k}}, "j": {"rz": {"spring": k}}}))[1].combinations["c"]
        assert c.invariants.passed
        m_end = abs(c.members["bm"].end_forces[11])
        if prev is not None:
            assert m_end >= prev                         # end moment grows with the connection stiffness
        prev = m_end
        if kfac == 1e-12:
            np.testing.assert_allclose(c.displacements, hinged.displacements, rtol=1e-9, atol=1e-12)
            np.testing.assert_allclose(c.members["bm"].forces, hinged.members["bm"].forces, rtol=1e-9, atol=1e-3)
        if kfac == 1e12:
            np.testing.assert_allclose(c.displacements, rigid.displacements, rtol=1e-9, atol=1e-12)
            np.testing.assert_allclose(c.members["bm"].forces, rigid.members["bm"].forces, rtol=1e-9, atol=1e-3)


@pytest.mark.parametrize("kfac", [1e-3, 1.0, 1e3])
def test_partial_release_equals_an_explicit_rotational_spring_between_two_coincident_nodes(kfac):
    """A semi-rigid connection of stiffness k at mid-span of a fixed-fixed beam, modelled (a) as a release
    spring of the left member's end j and (b) solved by the ODE with a rotation jump M / k at mid-span."""
    k = kfac * EIz / L
    doc = make({"a": (0, 0, 0), "b": (L / 2, 0, 0), "c": (L, 0, 0)},
               {"m1": {"i": "a", "j": "b", "releases": {"j": {"rz": {"spring": k}}}}, "m2": ("b", "c")},
               {"a": "fixed", "c": "fixed"}, {"c": [point("m1", "Z", -9000.0, 0.5), dist("m2", "Z", -2.0)]})
    c = run(doc)[1].combinations["c"]
    # ODE reference by superposition (the problem is linear): the continuous beam under the loads, plus
    # phi times the beam with a unit rotation jump at mid-span and no load; phi from M(mid) = k phi.
    h = L / 2
    base = beam_bvp(L, EIz, GAy, clamped, clamped, [(h, L, -2.0, -2.0)], [(h / 2, -9000.0, 0.0)])
    # unit rotation jump: left half clamped at 0 and right half clamped at L, each loaded at mid-span by the
    # unknown interface resultants (M0, V0); these follow from  v_L = v_R  and  theta_R - theta_L = 1
    def left(M0, V0):
        s = beam_bvp(h, EIz, GAy, clamped, lambda q: [q[2] - M0, q[3] - V0])(h)
        return s[0], s[1]

    def right(M0, V0):
        s = beam_bvp(h, EIz, GAy, lambda q: [q[2] - M0, q[3] - V0], clamped)(0.0, 1)
        return s[0], s[1]
    # linear in (M0, V0): build the 2x2 system  v_L = v_R,  th_R - th_L = 1
    cols = []
    for M0, V0 in ((1.0, 0.0), (0.0, 1.0)):
        (vl, tl), (vr, tr) = left(M0, V0), right(M0, V0)
        cols.append([vl - vr, tr - tl])
    A = np.array(cols).T
    M0, V0 = np.linalg.solve(A, [0.0, 1.0])
    # with rotation jump phi: M(mid) = base M(mid) + phi * M0 = k phi
    Mb = base(h)[2]
    phi = Mb / (k - M0)
    M_mid = k * phi
    mr = c.members["m1"]
    assert mr.end_forces[11] == pytest.approx(M_mid, rel=1e-9)
    v_mid = base(h)[0] + phi * left(M0, V0)[0]
    assert c.displacements[1, 2] == pytest.approx(v_mid, rel=1e-9)


def test_translational_release_transmits_no_force_in_that_direction():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (2 * L, 0, 0)},
               {"m1": {"i": "a", "j": "b", "releases": {"j": {"ux": True}}}, "m2": ("b", "c")},
               {"a": "fixed", "c": "fixed"}, {"c": [nodal("b", Fx=5000.0, Fz=-3000.0), dist("m1", "X", 2.0)]})
    c = run(doc)[1].combinations["c"]
    assert abs(c.members["m1"].end_forces[6]) <= 1e-9 * 5000
    assert c.members["m2"].end_forces[0] == pytest.approx(5000.0, rel=1e-10)        # all of Fx goes right
    assert c.reactions[0, 0] == pytest.approx(-2.0 * L, rel=1e-10)                  # m1's own axial load goes left
    assert c.invariants.passed


def test_releases_in_a_skew_rolled_member_are_in_member_local_axes():
    """Releasing local rz at the tip-side end: the end moment about LOCAL z vanishes, about local y not."""
    doc = make({"a": (0, 0, 0), "b": (2000, 1500, 2500), "c": (5000, 1500, 2500)},
               {"m1": {"i": "a", "j": "b", "orientation": {"roll": 33}, "releases": {"j": {"rz": True}}}, "m2": ("b", "c")},
               {"a": "fixed", "c": "fixed"}, {"c": [nodal("b", Fx=3000.0, Fy=-2000.0, Fz=-5000.0, Mx=1e6)]})
    c = run(doc)[1].combinations["c"]
    ef = c.members["m1"].end_forces
    assert abs(ef[11]) <= 1e-10 * abs(ef[10]) and abs(ef[10]) > 1e3
    assert c.invariants.passed


# ------------------------------------------------------------------------------- splitting
@pytest.mark.parametrize("shear", [True, False])
def test_member_split_into_pieces_gives_identical_results(shear):
    """A loaded skew beam with an overhanging tip, as 1 member and as 4 collinear pieces."""
    a, b = np.array([0.0, 0.0, 0.0]), np.array([4000.0, 2000.0, 1500.0])
    Lm = float(np.linalg.norm(b - a))
    ana = {"shearDeformation": shear, "resultStations": 9}
    ori = {"orientation": {"roll": 20}}
    loads1 = [dist("m", "Z", -3.0, w2=-7.0, start=0.1, end=0.9), point("m", "y", 5000.0, 0.3), point("m", "x", 2000.0, 0.6),
              point("m", "x", 3.0e5, 0.55, kind="moment"), point("m", "y", -4.0e5, 0.8, kind="moment")]
    one = make({"a": a, "b": b}, {"m": {"i": "a", "j": "b", **ori}}, {"a": "fixed", "b": {"dofs": {"uz": {"spring": 800.0}}}},
               {"c": loads1}, analysis=ana)
    cuts = [0.0, 0.25, 0.5, 0.75, 1.0]
    nodes = {f"n{k}": a + t * (b - a) for k, t in enumerate(cuts)}
    members = {f"p{k}": {"i": f"n{k}", "j": f"n{k + 1}", **ori} for k in range(4)}

    def piece(t):
        k = min(int(t / 0.25), 3)
        return f"p{k}", (t - cuts[k]) / 0.25
    w = lambda t: -3.0 + (-4.0) * (t - 0.1) / 0.8
    loads4 = []
    for k in range(4):
        t1, t2 = max(cuts[k], 0.1), min(cuts[k + 1], 0.9)
        loads4.append(dist(f"p{k}", "Z", w(t1), w2=w(t2), start=(t1 - cuts[k]) / 0.25, end=(t2 - cuts[k]) / 0.25))
    for t, d, v, kind in ((0.3, "y", 5000.0, "force"), (0.6, "x", 2000.0, "force"), (0.55, "x", 3.0e5, "moment"),
                          (0.8, "y", -4.0e5, "moment")):
        m, pos = piece(t)
        loads4.append(point(m, d, v, pos, kind=kind))
    four = make(nodes, members, {"n0": "fixed", "n4": {"dofs": {"uz": {"spring": 800.0}}}}, {"c": loads4}, analysis=ana)
    c1, c4 = run(one)[1].combinations["c"], run(four)[1].combinations["c"]
    scale_u = np.abs(c1.displacements).max()
    np.testing.assert_allclose(c4.displacements[-1], c1.displacements[-1], rtol=1e-9, atol=1e-12 * scale_u)
    np.testing.assert_allclose(c4.reactions[0], c1.reactions[0], rtol=1e-9, atol=1e-6)
    m1 = c1.members["m"]
    fscale = np.abs(m1.forces).max(axis=0)
    fscale = np.maximum(fscale, [fscale[:3].max()] * 3 + [fscale[3:].max()] * 3)
    for k in range(4):
        mp = c4.members[f"p{k}"]
        for s, xs in enumerate(mp.x):
            x = cuts[k] * Lm + xs
            # the same station (same side) of the single member, by position
            idx = [i for i, v in enumerate(m1.x) if abs(v - x) <= 1e-9 * Lm]
            if not idx:
                continue
            # at a point load the single member has two stations (left and right limit): take the matching side
            best = min(idx, key=lambda i: np.abs((m1.forces[i] - mp.forces[s]) / fscale).max())
            assert np.abs((m1.forces[best] - mp.forces[s]) / fscale).max() <= 1e-9
            np.testing.assert_allclose(mp.displacements[s], m1.displacements[best], rtol=1e-9, atol=1e-10 * scale_u)
    # internal node displacements = station displacements of the single member (local -> global)
    R = build_model(one).members["m"].R
    for k in (1, 2, 3):
        i = int(np.argmin(np.abs(m1.x - cuts[k] * Lm)))
        np.testing.assert_allclose(R @ c4.displacements[k, :3], m1.displacements[i], rtol=1e-9, atol=1e-10 * scale_u)


# --------------------------------------------------------------------------- rigid offsets
def test_offset_along_the_axis_makes_the_ends_rigid_ode_check():
    """Cantilever of node-to-node length L with a rigid first part of length a: tip deflection = that of a
    cantilever of length L - a (ODE), because the rigid part neither deforms nor rotates."""
    a = 700.0
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": {"i": "a", "j": "b", "offsets": {"i": [a, 0, 0]}}}, {"a": "fixed"},
               {"c": [nodal("b", Fz=-8000.0), dist("m", "Z", -2.0)]})
    model, res = run(doc)
    assert model.members["m"].length == pytest.approx(L - a)
    c = res.combinations["c"]
    ref = beam_bvp(L - a, EIz, GAy, clamped, lambda s: [s[2], s[3] + 8000.0], [(0, L - a, -2.0, -2.0)])
    assert c.displacements[1, 2] == pytest.approx(ref(L - a)[0], rel=REL)
    assert c.reactions[0, 2] == pytest.approx(8000.0 + 2.0 * (L - a), rel=1e-10)
    # reaction moment about global Y at the node, by statics about the node (the rigid part is a lever arm):
    # tip load -8000 at x = L and the line load -2 on [a, L] with its centroid at (L + a) / 2
    My_loads = 8000.0 * L + 2.0 * (L - a) * (L + a) / 2.0
    assert c.reactions[0, 4] == pytest.approx(-My_loads, rel=1e-10)
    # the member end moment at the flexible start is the ODE value (no lever arm)
    assert c.members["m"].end_forces[5] == pytest.approx(-ref(0.0, 1)[2], rel=REL)
    assert c.invariants.passed


def test_offset_not_along_the_axis_eccentric_beam_statics_and_kinematics():
    """Beam whose flexible axis is shifted sideways and upwards from the node line by equal offsets at both
    ends (an eccentric connection).  Node a is fixed; the tip load at node b acts with lever arms."""
    off = np.array([0.0, 150.0, -200.0])
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)},
               {"m": {"i": "a", "j": "b", "offsets": {"i": list(off), "j": list(off)}}}, {"a": "fixed"},
               {"c": [nodal("b", Fx=4000.0, Fy=1500.0, Fz=-6000.0)]})
    model, res = run(doc)
    c = res.combinations["c"]
    F = np.array([4000.0, 1500.0, -6000.0])
    # statics: reactions at a balance F at b
    np.testing.assert_allclose(c.reactions[0, :3], -F, rtol=1e-10)
    np.testing.assert_allclose(c.reactions[0, 3:], -np.cross([L, 0, 0], F), rtol=1e-10, atol=1e-6)
    # member end forces at the flexible end j: the node force moved rigidly to the flexible end (arm -off)
    m = model.members["m"]
    Mj = np.cross(-off, F)
    np.testing.assert_allclose(c.members["m"].end_forces[6:9], m.R @ F, rtol=1e-10, atol=1e-7)
    np.testing.assert_allclose(c.members["m"].end_forces[9:12], m.R @ Mj, rtol=1e-10, atol=1e-4)
    # kinematics: node b = flexible end j carried back through the rigid arm; flexible end from the ODE
    Fl, Ml = m.R @ F, m.R @ Mj
    xy = beam_bvp(L, EIz, GAy, clamped, lambda s: [s[2] - Ml[2], s[3] - Fl[1]])
    xz = beam_bvp(L, EIy, GAz, clamped, lambda s: [s[2] - Ml[1], s[3] - Fl[2]], plane="xz")
    ax = bar_bvp(L, EA, lambda s: [s[0]], lambda s: [s[1] - Fl[0]])
    to = bar_bvp(L, GJ, lambda s: [s[0]], lambda s: [s[1] - Ml[0]])
    u_end = m.R.T @ np.array([ax(L)[0], xy(L)[0], xz(L)[0]])
    th = m.R.T @ np.array([to(L)[0], xz(L)[1], xy(L)[1]])
    u_node = u_end - np.cross(th, off)
    np.testing.assert_allclose(c.displacements[1, 3:], th, rtol=0, atol=REL * np.abs(th).max())
    np.testing.assert_allclose(c.displacements[1, :3], u_node, rtol=0, atol=REL * np.abs(u_node).max())


def test_offset_converges_to_a_very_stiff_stub_member_with_an_added_node():
    """Rigid offset (exact) vs. a stub member of increasing stiffness between the node and the flexible end."""
    off_i, off_j = np.array([300.0, 100.0, -150.0]), np.array([-250.0, -80.0, 120.0])
    a, b = np.array([0.0, 0.0, 0.0]), np.array([4000.0, 1000.0, 500.0])
    loads = [nodal("b", Fx=2000.0, Fy=-3000.0, Fz=-5000.0, Mx=4e5), dist("m", "Z", -3.0)]
    ori = {"orientation": {"referenceVector": [0, 0, 1]}}
    exact = run(make({"a": a, "b": b}, {"m": {"i": "a", "j": "b", "offsets": {"i": list(off_i), "j": list(off_j)}, **ori}},
                     {"a": "fixed", "b": {"dofs": {"uz": {"spring": 500.0}, "uy": {"spring": 300.0}}}}, {"c": loads}))[1].combinations["c"]
    errs = []
    for factor in (1e2, 1e4, 1e6):
        secs = {"S": {"type": "properties", **SEC},
                "stub": {"type": "properties", **{k: v * factor for k, v in SEC.items()}}}
        doc = make({"a": a, "b": b, "ea": a + off_i, "eb": b + off_j},
                   {"m": {"i": "ea", "j": "eb", "section": "S", **ori}, "si": {"i": "a", "j": "ea", "section": "stub"},
                    "sj": {"i": "eb", "j": "b", "section": "stub"}},
                   {"a": "fixed", "b": {"dofs": {"uz": {"spring": 500.0}, "uy": {"spring": 300.0}}}}, {"c": loads}, sections=secs)
        c = run(doc)[1].combinations["c"]
        e_u = np.abs(c.displacements[1] - exact.displacements[1]).max() / np.abs(exact.displacements[1]).max()
        e_f = np.abs(c.members["m"].forces - exact.members["m"].forces).max() / np.abs(exact.members["m"].forces).max()
        errs.append(max(e_u, e_f))
    assert errs[0] > errs[1] > errs[2] and errs[2] < 1e-5
    assert errs[1] / errs[0] < 0.02 and errs[2] / errs[1] < 0.02          # error proportional to 1 / stiffness


# ----------------------------------------------------------------------------- rigid links
def link_model(kind="link", factor=None):
    """Two cantilever columns whose tops are tied; the right top (slave) is offset in all three directions."""
    nodes = {"a": (0, 0, 0), "m": (0, 0, 3000), "d": (2500, 800, 0), "s": (2500, 800, 3400)}
    members = {"c1": ("a", "m"), "c2": ("d", "s")}
    loads = [nodal("s", Fx=5000.0, Fy=-2000.0, Fz=-9000.0, Mz=3e6), nodal("m", Fy=1000.0)]
    kw = {}
    secs = None
    if kind == "link":
        kw["constraints"] = {"L": {"type": "rigidLink", "master": "m", "slaves": ["s"]}}
    else:
        secs = {"S": {"type": "properties", **SEC}, "stiff": {"type": "properties", **{k: v * factor for k, v in SEC.items()}}}
        members["tie"] = {"i": "m", "j": "s", "section": "stiff"}
    return make(nodes, members, {"a": "fixed", "d": "fixed"}, {"c": loads}, sections=secs, **kw)


def test_rigid_link_is_the_exact_rigid_body_relation_and_reports_constraint_forces():
    model, res = run(link_model())
    c = res.combinations["c"]
    um, us = c.displacements[1], c.displacements[3]
    r = model.nodes["s"].xyz - model.nodes["m"].xyz
    np.testing.assert_allclose(us[3:], um[3:], rtol=0, atol=0)
    np.testing.assert_allclose(us[:3], um[:3] + np.cross(um[3:], r), rtol=1e-14, atol=1e-18)
    assert c.invariants.passed
    lam = c.constraint_forces
    assert not lam[[0, 2]].any()
    # constraint forces are self-equilibrated: force and moment (about the master) of the pair vanish
    np.testing.assert_allclose(lam[1, :3] + lam[3, :3], 0, atol=1e-8 * np.abs(lam).max())
    np.testing.assert_allclose(lam[1, 3:] + lam[3, 3:] + np.cross(r, lam[3, :3]), 0, atol=1e-8 * np.abs(lam).max())
    # slave joint: applied load + constraint force = what column c2 takes
    ef = c.members["c2"].end_forces
    R2 = model.members["c2"].R
    np.testing.assert_allclose(R2.T @ ef[6:9], np.array([5000.0, -2000.0, -9000.0]) + lam[3, :3], rtol=1e-9, atol=1e-6)
    doc = res.to_document(model)
    assert set(doc["combinations"]["c"]["constraintForces"]) == {"m", "s"}


def test_rigid_link_is_the_limit_of_a_very_stiff_member():
    exact = run(link_model())[1].combinations["c"]
    errs = []
    for factor in (1e2, 1e4, 1e6):
        c = run(link_model("member", factor))[1].combinations["c"]
        errs.append(max(np.abs(c.displacements - exact.displacements).max() / np.abs(exact.displacements).max(),
                        np.abs(c.reactions - exact.reactions).max() / np.abs(exact.reactions).max()))
    # a stiff member is a penalty: its error is proportional to 1 / stiffness (factor 100 per step)
    assert errs[0] > errs[1] > errs[2] and errs[2] < 1e-4 and errs[1] / errs[0] < 0.02 and errs[2] / errs[1] < 0.02


def test_rigid_link_with_coincident_master_and_slave_ties_two_nodes_together():
    """Zero-length link = merged nodes: same result as one node."""
    tied = make({"a": (0, 0, 0), "b": (L, 0, 0), "b2": (L, 0, 0), "c": (2 * L, 0, 0)}, {"m1": ("a", "b"), "m2": ("b2", "c")},
                {"a": "fixed", "c": "pinned"}, {"c": [nodal("b2", Fz=-7000.0, My=2e6)]},
                constraints={"t": {"type": "rigidLink", "master": "b", "slaves": ["b2"]}})
    one = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (2 * L, 0, 0)}, {"m1": ("a", "b"), "m2": ("b", "c")},
               {"a": "fixed", "c": "pinned"}, {"c": [nodal("b", Fz=-7000.0, My=2e6)]})
    ct, co = run(tied)[1].combinations["c"], run(one)[1].combinations["c"]
    np.testing.assert_allclose(ct.displacements[1], co.displacements[1], rtol=1e-11, atol=1e-16)
    np.testing.assert_allclose(ct.displacements[2], co.displacements[1], rtol=1e-11, atol=1e-16)
    np.testing.assert_allclose(ct.members["m2"].forces, co.members["m2"].forces, rtol=1e-10, atol=1e-6)


def test_master_with_a_support_carries_the_slave_loads_into_its_reaction():
    doc = make({"m": (0, 0, 0), "s": (1000, 500, 200), "t": (1000, 500, 3000)}, {"col": ("s", "t")},
               {"m": "fixed"}, {"c": [nodal("t", Fx=4000.0, Fz=-1000.0)]},
               constraints={"L": {"type": "rigidLink", "master": "m", "slaves": ["s"]}})
    model, res = run(doc)
    c = res.combinations["c"]
    F = np.array([4000.0, 0, -1000.0])
    np.testing.assert_allclose(c.reactions[0, :3], -F, rtol=1e-10, atol=1e-9)
    np.testing.assert_allclose(c.reactions[0, 3:], -np.cross([1000, 500, 3000], F), rtol=1e-10, atol=1e-6)
    assert not c.displacements[1].any() and c.invariants.passed


# ------------------------------------------------------------------------------ diaphragm
def building(diaphragm=True, tie_factor=None, shift=0.0):
    """One-storey 3D frame: four columns of different orientation, perimeter beams, master at the centre."""
    H, a, b = 3200.0, 6000.0, 4000.0
    nodes = {"b1": (0, 0, 0), "b2": (a, 0, 0), "b3": (a, b, 0), "b4": (0, b, 0),
             "t1": (0, 0, H), "t2": (a, 0, H), "t3": (a, b, H), "t4": (0, b, H + shift), "cm": (a / 2, b / 2, H)}
    members = {f"c{k}": {"i": f"b{k}", "j": f"t{k}", "orientation": {"roll": 90.0 if k % 2 else 0.0}} for k in (1, 2, 3, 4)}
    members.update({"g1": ("t1", "t2"), "g2": ("t2", "t3"), "g3": ("t4", "t3"), "g4": ("t1", "t4")})
    secs = {"S": {"type": "properties", **SEC}}
    kw = {}
    if diaphragm:
        kw["constraints"] = {"floor": {"type": "rigidDiaphragm", "master": "cm", "slaves": ["t1", "t2", "t3", "t4"]}}
    elif tie_factor:
        # in-plane rigid, out-of-plane soft "slab": spokes that are stiff axially and in in-plane bending only
        secs["spoke"] = {"type": "properties", "A": SEC["A"] * tie_factor, "Iy": SEC["Iz"] * tie_factor,
                         "Iz": 1e-3, "J": 1e-3}
        for k in (1, 2, 3, 4):
            members[f"s{k}"] = {"i": "cm", "j": f"t{k}", "section": "spoke",
                                "releases": {"j": {"rx": True, "rz": True}}}
    loads = [nodal("cm", Fx=40000.0, Fy=15000.0, Mz=3.0e7), nodal("t2", Fz=-20000.0),
             {"type": "memberDistributed", "member": "g1", "direction": "Z", "w1": -5.0}]
    sup = {f"b{k}": "fixed" for k in (1, 2, 3, 4)}
    if not diaphragm:
        sup["cm"] = {"dofs": {"uz": "fixed", "rx": "fixed", "ry": "fixed"}}
    return make(nodes, members, sup, {"c": loads}, sections=secs, **kw)


def test_rigid_diaphragm_kinematics_and_free_out_of_plane_dofs():
    model, res = run(building())
    c = res.combinations["c"]
    assert c.invariants.passed
    um = c.displacements[model.nodes["cm"].index]
    for nid in ("t1", "t2", "t3", "t4"):
        us = c.displacements[model.nodes[nid].index]
        r = model.nodes[nid].xyz - model.nodes["cm"].xyz
        assert us[0] == pytest.approx(um[0] - um[5] * r[1], rel=1e-13, abs=1e-18)
        assert us[1] == pytest.approx(um[1] + um[5] * r[0], rel=1e-13, abs=1e-18)
        assert us[5] == um[5]
    # out-of-plane DOFs are independent: the loaded corner goes down more than the others, slabs rotate freely
    uz = {n: c.displacements[model.nodes[n].index, 2] for n in ("t1", "t2", "t3", "t4")}
    assert uz["t2"] < uz["t4"] and len({round(v, 12) for v in uz.values()}) == 4
    # the master has no member: its uz, rx, ry are not part of the solution (reported 0, with a warning)
    assert um[2] == um[3] == um[4] == 0.0
    assert any("not connected to any member" in w for w in res.warnings)
    # storey shear: base reactions balance the diaphragm loads
    np.testing.assert_allclose(c.reactions[:, 0].sum(), -40000.0, rtol=1e-9)
    np.testing.assert_allclose(c.reactions[:, 1].sum(), -15000.0, rtol=1e-9)
    # constraint forces: in-plane only, self-equilibrated in force and in moment about the normal
    lam = c.constraint_forces
    assert not lam[:, [2, 3, 4]].any()
    assert abs(lam[:, 0].sum()) <= 1e-8 * np.abs(lam).max() and abs(lam[:, 1].sum()) <= 1e-8 * np.abs(lam).max()
    mz = sum(lam[n.index, 5] + n.xyz[0] * lam[n.index, 1] - n.xyz[1] * lam[n.index, 0] for n in model.nodes.values())
    assert abs(mz) <= 1e-9 * 3.0e7


def test_rigid_diaphragm_is_the_limit_of_stiff_in_plane_ties():
    exact_model, exact_res = run(building())
    exact = exact_res.combinations["c"]
    idx = [exact_model.nodes[n].index for n in ("t1", "t2", "t3", "t4")]
    errs = []
    for factor in (1e2, 1e4, 1e6):
        model, res = run(building(False, factor))
        c = res.combinations["c"]
        jdx = [model.nodes[n].index for n in ("t1", "t2", "t3", "t4")]
        errs.append(np.abs(c.displacements[jdx][:, [0, 1, 5]] - exact.displacements[idx][:, [0, 1, 5]]).max()
                    / np.abs(exact.displacements[idx][:, [0, 1, 5]]).max())
    assert errs[0] > errs[1] > errs[2] and errs[2] < 1e-4


def test_diaphragm_with_out_of_plane_slave_and_unbalanced_moment_is_refused_not_published():
    """A slave above the master's plane receives in-plane constraint forces at its own level: the overall
    moment balance is violated by force x level difference.  That result is flagged, never returned as good."""
    doc = building(shift=400.0)
    model = build_model(doc)
    assert any(p.code == "constraint.diaphragm-out-of-plane" for p in model.warnings)
    with pytest.raises(ModelRefused) as e:
        analyze(model)
    assert e.value.codes == ["solve.invariant-failed"] and "global moment equilibrium" in e.value.problems[0].message
    assert "diaphragm" in e.value.problems[0].hint
    flagged = analyze(model, strict=False)
    c = flagged.combinations["c"]
    assert not c.invariants.passed and any("FLAGGED" in w for w in flagged.warnings)
    failed = [ch.name for ch in c.invariants.checks if not ch.passed]
    assert failed == ["global moment equilibrium"]


def test_diaphragm_slave_keeps_its_own_vertical_support_and_load():
    doc = make({"m": (0, 0, 0), "s": (2000, 0, 0), "t": (2000, 0, 3000)}, {"col": ("s", "t"), "bm": ("m", "s")},
               {"m": {"type": "pinned", "dofs": {"rz": "fixed", "rx": "fixed"}}, "s": {"dofs": {"uz": "fixed"}}},
               {"c": [nodal("t", Fz=-5000.0, Fy=800.0)]},
               constraints={"d": {"type": "rigidDiaphragm", "master": "m", "slaves": ["s"]}})
    c = run(doc)[1].combinations["c"]
    assert c.invariants.passed
    assert c.reactions[0, 2] + c.reactions[1, 2] == pytest.approx(5000.0, rel=1e-9)       # vertical equilibrium
    assert c.reactions[1, 2] != 0.0 and not c.reactions[1, [0, 1, 3, 4, 5]].any()         # the slave's own support: uz only
    assert c.displacements[1, 0] == 0.0 and c.displacements[1, 1] == 0.0 and c.displacements[1, 2] == 0.0


# ---------------------------------------------------------------------------------- truss
def truss_doc(extra_loads=()):
    nodes = {"a": (0, 0, 0), "b": (4000, 0, 0), "c": (2000, 0, 3000), "d": (2000, 2500, 1000)}
    members = {k: {"i": i, "j": j, "type": "truss"} for k, (i, j) in
               {"ab": ("a", "b"), "ac": ("a", "c"), "bc": ("b", "c"), "ad": ("a", "d"), "bd": ("b", "d"), "cd": ("c", "d")}.items()}
    return make(nodes, members, {"a": "pinned", "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}, "d": {"dofs": {"uy": "fixed"}}},
                {"c": [nodal("c", Fx=10000.0, Fz=-20000.0), *extra_loads]})


def test_space_truss_bar_forces_by_the_method_of_joints():
    """Statically determinate: bar forces from joint equilibrium solved independently with numpy."""
    model, res = run(truss_doc())
    c = res.combinations["c"]
    assert c.invariants.passed
    # independent statics: unknowns = 6 bar forces + 6 reaction components; 12 joint equations
    X = {k: n.xyz for k, n in model.nodes.items()}
    bars = [(m.id, m.i, m.j) for m in model.members.values()]
    react = [("a", 0), ("a", 1), ("a", 2), ("b", 1), ("b", 2), ("d", 1)]
    names = list(X)
    A = np.zeros((12, 12))
    rhs = np.zeros(12)
    for col, (_, i, j) in enumerate(bars):
        e = (X[j] - X[i]) / np.linalg.norm(X[j] - X[i])
        A[3 * names.index(i):3 * names.index(i) + 3, col] += e          # tension pulls node i towards j
        A[3 * names.index(j):3 * names.index(j) + 3, col] -= e
    for col, (n, a) in enumerate(react):
        A[3 * names.index(n) + a, 6 + col] = 1.0
    rhs[3 * names.index("c") + 0] = -10000.0
    rhs[3 * names.index("c") + 2] = 20000.0
    sol = np.linalg.solve(A, rhs)
    for col, (mid, _, _) in enumerate(bars):
        mr = c.members[mid]
        assert mr.forces[0, 0] == pytest.approx(sol[col], rel=1e-9, abs=1e-6)
        assert np.abs(mr.forces[:, 1:]).max() <= 1e-9 * 20000.0          # axial only
    for col, (n, a) in enumerate(react):
        assert c.reactions[model.nodes[n].index, a] == pytest.approx(sol[6 + col], rel=1e-9, abs=1e-6)
    # rotations of truss-only nodes are not part of the solution
    assert not c.displacements[:, 3:].any()
    assert any("only truss-type members" in w for w in res.warnings)


def test_truss_member_load_is_carried_to_the_nodes_as_simply_supported_beam_reactions():
    extra = [{"type": "memberPoint", "member": "ab", "direction": "Z", "value": -6000.0, "position": 0.25},
             {"type": "selfWeight", "members": ["ab"]}]
    model, res = run(truss_doc(extra))
    base = run(truss_doc())[1].combinations["c"]
    w = 7.7e-5 * SEC["A"]
    equiv = run(truss_doc([nodal("a", Fz=-(4500.0 + w * 2000.0)), nodal("b", Fz=-(1500.0 + w * 2000.0))]))[1].combinations["c"]
    c = res.combinations["c"]
    np.testing.assert_allclose(c.displacements, equiv.displacements, rtol=1e-10, atol=1e-13)
    for mid in ("ac", "bc", "ad", "bd", "cd"):
        np.testing.assert_allclose(c.members[mid].forces, equiv.members[mid].forces, rtol=1e-9, atol=1e-6)
    mr = c.members["ab"]
    ref = beam_bvp(4000.0, EIz, GAy, pinned, pinned, [(0, 4000.0, -w, -w)], [(1000.0, -6000.0, 0.0)])
    side = [-1 if (k + 1 < len(mr.x) and mr.x[k + 1] == x) else 1 for k, x in enumerate(mr.x)]
    for k, x in enumerate(mr.x):
        s = side[k] if 0 < x < 4000 else (1 if x == 0 else -1)
        assert abs(mr.forces[k, 5] - ref(x, s)[2]) <= REL * abs(ref(1000.0)[2])
        assert abs(mr.forces[k, 1] - ref(x, s)[3]) <= REL * 6000.0
    assert mr.end_forces[5] == pytest.approx(0.0, abs=1e-6) and mr.end_forces[11] == pytest.approx(0.0, abs=1e-6)
    assert c.invariants.passed and base.invariants.passed


def test_torque_on_a_truss_member_is_refused():
    with pytest.raises(ModelRefused) as e:
        run(truss_doc([{"type": "memberPoint", "member": "ab", "direction": "x", "kind": "moment", "value": 1e5, "position": 0.5}]))
    assert e.value.codes == ["load.torque-on-truss"]


def test_moment_on_a_truss_only_node_is_refused_as_unstable():
    with pytest.raises(ModelUnstable) as e:
        run(truss_doc([nodal("c", My=1e5)]))
    d = e.value.diagnosis
    assert d.modes[0].kind == "node" and d.modes[0].cause == "only-truss-members" and d.modes[0].nodes[0].node == "c"
    assert e.value.problems[0].code == "stability.mechanism" and e.value.problems[0].path == "/nodes/c"


def test_truss_bracing_in_a_frame():
    def frame(brace):
        members = {"c1": ("a", "b"), "bm": ("b", "c"), "c2": ("d", "c")}
        if brace:
            members["br"] = {"i": "a", "j": "c", "type": "truss"}
        return make({"a": (0, 0, 0), "b": (0, 0, 3000), "c": (4000, 0, 3000), "d": (4000, 0, 0)}, members,
                    {"a": "fixed", "d": "fixed"}, {"c": [nodal("b", Fx=30000.0)]})
    c = run(frame(True))[1].combinations["c"]
    br = c.members["br"]
    assert br.forces[0, 0] > 0 and np.abs(br.forces[:, 1:]).max() <= 1e-9 * 30000        # tension, axial only
    np.testing.assert_allclose(br.forces[:, 0], br.forces[0, 0], rtol=1e-12)
    unbraced = run(frame(False))[1].combinations["c"]
    assert abs(c.displacements[1, 0]) < 0.2 * abs(unbraced.displacements[1, 0])
    # the brace force follows from its elongation: N = EA / L * (u_c - u_a) . e
    e = np.array([4000.0, 0.0, 3000.0]) / 5000.0
    assert br.forces[0, 0] == pytest.approx(EA / 5000.0 * (c.displacements[2, :3] @ e), rel=1e-9)
    assert c.invariants.passed


def test_plane_frame_on_pinned_bases_is_unstable_out_of_plane_and_says_so():
    doc = make({"a": (0, 0, 0), "b": (0, 0, 3000), "c": (4000, 0, 3000), "d": (4000, 0, 0)},
               {"c1": ("a", "b"), "bm": ("b", "c"), "c2": ("d", "c")}, {"a": "pinned", "d": "pinned"},
               {"c": [nodal("b", Fx=30000.0)]})
    with pytest.raises(ModelUnstable) as e:
        run(doc)
    d = e.value.diagnosis
    assert d.n_mechanisms == 1 and d.modes[0].kind == "rigid-body" and d.modes[0].motion == "rotation"
    assert "rotation about an axis parallel to X through node a" in d.modes[0].description


# --------------------------------------------------------------- combinations, envelopes
def test_combinations_are_the_factored_sum_of_cases_and_envelopes_name_the_governing_combination():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (2 * L, 0, 0)}, {"m1": ("a", "b"), "m2": ("b", "c")},
               {"a": "pinned", "b": {"dofs": {"uz": "fixed", "uy": "fixed"}}, "c": {"dofs": {"uz": "fixed", "uy": "fixed", "rx": "fixed"}}},
               {"D": [dist("m1", "Z", -4.0), dist("m2", "Z", -4.0)], "L1": [dist("m1", "Z", -10.0)],
                "L2": [dist("m2", "Z", -10.0)], "W": [nodal("b", Fx=1000.0), point("m2", "Z", 5000.0, 0.5)]},
               combos={"D": {"factors": {"D": 1.0}}, "L1": {"factors": {"L1": 1.0}}, "L2": {"factors": {"L2": 1.0}},
                       "W": {"factors": {"W": 1.0}},
                       "U1": {"factors": {"D": 1.2, "L1": 1.6}}, "U2": {"factors": {"D": 1.2, "L2": 1.6}},
                       "U3": {"factors": {"D": 0.9, "W": -1.0, "L1": 0.5}}})
    model, res = run(doc)
    C = res.combinations
    for cid, fac in (("U1", {"D": 1.2, "L1": 1.6}), ("U2", {"D": 1.2, "L2": 1.6}), ("U3", {"D": 0.9, "W": -1.0, "L1": 0.5})):
        for attr in ("displacements", "reactions"):
            expect = sum(f * getattr(C[k], attr) for k, f in fac.items())
            np.testing.assert_allclose(getattr(C[cid], attr), expect, rtol=1e-10, atol=1e-9 * np.abs(expect).max())
        for mid in ("m1", "m2"):
            expect = sum(f * C[k].members[mid].forces for k, f in fac.items())
            np.testing.assert_allclose(C[cid].members[mid].forces, expect, rtol=1e-10, atol=1e-9 * np.abs(expect).max())
    env = res.envelopes["m1"]
    ids = list(C)
    stack = np.array([C[k].members["m1"].forces for k in ids])
    np.testing.assert_array_equal(env.max, stack.max(axis=0))
    np.testing.assert_array_equal(env.min, stack.min(axis=0))
    k_mid = int(np.argmin(np.abs(env.x - 0.4 * L)))
    assert env.max_combination[k_mid][5] == "U1" and env.min_combination[k_mid][5] in ("L2", "U2")
    for s in range(len(env.x)):
        for q in range(6):
            assert C[env.max_combination[s][q]].members["m1"].forces[s, q] == env.max[s, q]
            assert C[env.min_combination[s][q]].members["m1"].forces[s, q] == env.min[s, q]
    ne = res.nodal_envelopes
    b = model.nodes["b"].index
    assert ne["reactions"].max[b, 2] == max(C[k].reactions[b, 2] for k in ids)
    assert C[ne["reactions"].max_combination[b][2]].reactions[b, 2] == ne["reactions"].max[b, 2]
    assert C[ne["displacements"].min_combination[b][4]].displacements[b, 4] == ne["displacements"].min[b, 4]
    doc_out = res.to_document(model)
    assert doc_out["envelopes"]["m1"]["maxCombination"]["Mz"][k_mid] == "U1"
    assert doc_out["nodalEnvelopes"]["reactions"]["b"]["maxCombination"]["Fz"] == ne["reactions"].max_combination[b][2]


def test_zero_load_combination_gives_zero_results_and_passes_invariants():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed"}, {"c": [nodal("b", Fz=-1.0)]},
               combos={"zero": {"factors": {"c": 0.0}}})
    c = run(doc)[1].combinations["zero"]
    assert not c.displacements.any() and not c.reactions.any() and not c.members["m"].forces.any()
    assert c.invariants.passed


# --------------------------------------------------------- rigid rotation of the whole model
def test_results_are_objective_under_a_rigid_rotation_of_the_whole_model():
    """The same space frame described in a rotated global frame: displacements and reactions rotate, member
    forces (local) are unchanged.  Orientation is kept by rotating the reference vectors too."""
    rng = np.random.default_rng(11)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    P = {"a": (0, 0, 0), "b": (500, 300, 3200), "c": (4500, 800, 3000), "d": (5000, 0, 0), "e": (4800, 3500, 2800), "f": (5200, 4000, 0)}
    refs = {"c1": (1, 0, 0), "bm": (0, 0, 1), "c2": (1, 0.2, 0), "g": (0, 0, 1), "c3": (0, 1, 0)}
    conn = {"c1": ("a", "b"), "bm": ("b", "c"), "c2": ("d", "c"), "g": ("c", "e"), "c3": ("f", "e")}
    F = {"b": np.array([8000.0, -3000.0, -12000.0, 2e5, -4e5, 6e5]), "e": np.array([-2000.0, 5000.0, -7000.0, 0, 0, 0])}
    rel = {"bm": {"j": {"rz": {"spring": 5.0e8}}}, "g": {"i": {"ry": True}}}
    off = {"bm": {"i": [120.0, 0.0, -80.0]}}

    def build(Rm):
        nodes = {k: Rm @ np.array(v, float) for k, v in P.items()}
        members = {}
        for k, (i, j) in conn.items():
            m = {"i": i, "j": j, "orientation": {"referenceVector": list(Rm @ np.array(refs[k], float)), "roll": 15.0}}
            if k in rel:
                m["releases"] = rel[k]
            if k in off:
                m["offsets"] = {kk: list(Rm @ np.array(v)) for kk, v in off[k].items()}
            members[k] = m
        loads = [nodal(n, **dict(zip(("Fx", "Fy", "Fz", "Mx", "My", "Mz"), np.concatenate([Rm @ v[:3], Rm @ v[3:]])))) for n, v in F.items()]
        loads += [dist("bm", "y", -4.0, w2=-1.0, start=0.2), point("g", "z", 3000.0, 0.4)]
        sup = {"a": "fixed", "d": "fixed", "f": "fixed"}
        return make(nodes, members, sup, {"c": loads})

    m0, r0 = run(build(np.eye(3)))
    m1, r1 = run(build(Q))
    c0, c1 = r0.combinations["c"], r1.combinations["c"]
    su, sr = np.abs(c0.displacements[:, :3]).max(), np.abs(c0.displacements[:, 3:]).max()
    np.testing.assert_allclose(c1.displacements[:, :3], c0.displacements[:, :3] @ Q.T, rtol=0, atol=1e-9 * su)
    np.testing.assert_allclose(c1.displacements[:, 3:], c0.displacements[:, 3:] @ Q.T, rtol=0, atol=1e-9 * sr)
    np.testing.assert_allclose(c1.reactions[:, :3], c0.reactions[:, :3] @ Q.T, rtol=0, atol=1e-9 * np.abs(c0.reactions[:, :3]).max())
    for mid in conn:
        s = np.abs(c0.members[mid].forces).max(axis=0)
        s = np.maximum(s, [s[:3].max()] * 3 + [s[3:].max()] * 3)
        assert (np.abs(c1.members[mid].forces - c0.members[mid].forces) <= 1e-9 * s).all()
        np.testing.assert_allclose(c1.members[mid].displacements, c0.members[mid].displacements, rtol=0, atol=1e-9 * su)


def test_results_document_is_in_file_units_and_carries_defaults_and_warnings(tmp_path):
    """The same beam in (kN, m) and in (N, mm): identical physics, numbers scaled by the unit factors."""
    def doc(force, length):
        f, l = {"N": 1.0, "kN": 1e-3}[force], {"mm": 1.0, "m": 1e-3}[length]
        return {"formatVersion": "1.0.0", "units": {"force": force, "length": length}, "conventions": {"verticalAxis": "Z"},
                "nodes": {"a": {"x": 0, "y": 0, "z": 0}, "b": {"x": L * l, "y": 0, "z": 0}, "lonely": {"x": 0, "y": 1000 * l, "z": 0}},
                "supports": {"a": {"type": "fixed"}},
                "materials": {"M": {"E": E * f / l ** 2, "G": G * f / l ** 2}},
                "sections": {"S": {"type": "properties", **{k: v * l ** (4 if k[0] in "IJ" else 2) for k, v in SEC.items()}}},
                "members": {"m": {"i": "a", "j": "b", "section": "S", "material": "M"}},
                "loadCases": {"c": {"loads": [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -3.0 * f / l},
                                              {"type": "nodal", "node": "b", "My": 2.0e6 * f * l}]}},
                "combinations": {"c": {"factors": {"c": 1.0}}}}
    mn, rn = run(doc("N", "mm"))
    mk, rk = run(doc("kN", "m"))
    dn, dk = rn.to_document(mn), rk.to_document(mk)
    assert dk["units"] == {"force": "kN", "length": "m"}
    cn, ck = dn["combinations"]["c"], dk["combinations"]["c"]
    assert ck["displacements"]["b"]["uz"] == pytest.approx(cn["displacements"]["b"]["uz"] * 1e-3, rel=1e-10)
    assert ck["displacements"]["b"]["ry"] == pytest.approx(cn["displacements"]["b"]["ry"], rel=1e-10)
    assert ck["reactions"]["a"]["Fz"] == pytest.approx(cn["reactions"]["a"]["Fz"] * 1e-3, rel=1e-10)
    assert ck["reactions"]["a"]["My"] == pytest.approx(cn["reactions"]["a"]["My"] * 1e-6, rel=1e-10)
    st_n, st_k = cn["members"]["m"]["stations"], ck["members"]["m"]["stations"]
    np.testing.assert_allclose(st_k["x"], np.array(st_n["x"]) * 1e-3, rtol=1e-12)
    np.testing.assert_allclose(st_k["Mz"], np.array(st_n["Mz"]) * 1e-6, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(st_k["uy"], np.array(st_n["uy"]) * 1e-3, rtol=1e-10, atol=1e-15)
    assert list(ck["reactions"]) == ["a"]
    assert ck["invariants"]["passed"] is True and len(ck["invariants"]["checks"]) >= 9
    assert all({"name", "value", "tolerance", "passed", "detail"} <= set(ch) for ch in ck["invariants"]["checks"])
    assert len(dk["defaultsUsed"]) == len(mk.defaults_used) > 0
    assert any("node.unused" in w for w in dk["warnings"]) and any("not connected to any member" in w for w in dk["warnings"])
    import json
    json.dumps(dk, allow_nan=False)
