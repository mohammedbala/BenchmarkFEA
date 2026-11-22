"""Element level: exact integrator, cantilever flexibility, stiffness, releases, fixed-end forces.

Reference values come from an independent collocation solution of the governing first-order system
(``solve_reference.beam_bvp`` / ``bar_bvp``), never from quoted formulas.
"""
import dataclasses

import numpy as np
import pytest

from benchfea.model import RIGID, build_model
from benchfea.solve import ModelRefused, SectionFlex, SpanLoads, build_element, cantilever_flexibility, integrate
from benchfea.solve.element import release_stiffness, rigid_transfer

from solve_reference import (EA, EIy, EIz, GAy, GAz, GJ, bar_bvp, beam_bvp, clamped, free, make)

L = 3200.0
FLEX = SectionFlex(1 / EA, 1 / GJ, 1 / EIy, 1 / EIz, 1 / GAy, 1 / GAz)
FLEX_EB = SectionFlex(1 / EA, 1 / GJ, 1 / EIy, 1 / EIz)
INF = RIGID
RR = (INF,) * 6


def close(a, b, rel=1e-9, scale=None):
    a, b = np.asarray(a, float), np.asarray(b, float)
    s = scale if scale is not None else max(np.abs(b).max(), 1e-300)
    assert np.abs(a - b).max() <= rel * s, (a, b)


def f6(**kw):
    v = np.zeros(6)
    for k, x in kw.items():
        v["Fx Fy Fz Mx My Mz".split().index(k)] = x
    return v


def w3(x=0.0, y=0.0, z=0.0):
    return np.array([x, y, z], float)


# ----------------------------------------------------------------------------- integrator
@pytest.mark.parametrize("shear", [True, False])
def test_integrator_reproduces_the_ode_solution_of_a_cantilever_in_the_xy_plane(shear):
    loads = SpanLoads(points=((0.35 * L, f6(Fy=-4000.0)), (0.7 * L, f6(Mz=2.5e6))),
                      dists=((0.1 * L, 0.9 * L, w3(y=-2.0), w3(y=-6.5)), (0.0, L, w3(y=1.0), w3(y=1.0))))
    flex = FLEX if shear else FLEX_EB
    run = integrate(L, flex, loads, loads.resultant(), np.zeros(6))
    ref = beam_bvp(L, EIz, GAy if shear else None, clamped, free,
                   [(0.1 * L, 0.9 * L, -2.0, -6.5), (0.0, L, 1.0, 1.0)], [(0.35 * L, -4000.0, 0.0), (0.7 * L, 0.0, 2.5e6)])
    close(run.s_end, np.zeros(6), scale=np.abs(run.s_start).max(), rel=1e-12)        # free end: no resultants
    scale_d, scale_m, scale_v = abs(ref(L)[0]), abs(ref(0.0, +1)[2]), abs(ref(0.0, +1)[3])
    for x in np.linspace(0, L, 17):
        for side in (-1, 1):
            s, y = run.at(float(x), side)
            d, th, M, V = ref(x, side)
            close([y[1]], [d], scale=scale_d)
            close([y[5]], [th], scale=abs(ref(L)[1]))
            close([s[5]], [M], scale=scale_m)
            close([s[1]], [V], scale=scale_v)
    for x, side in ((0.35 * L, -1), (0.35 * L, 1), (0.7 * L, -1), (0.7 * L, 1)):         # the jumps
        s, _ = run.at(x, side)
        d, th, M, V = ref(x, side)
        close([s[1]], [V], scale=scale_v)
        close([s[5]], [M], scale=scale_m)
    s_before, s_after = run.at(0.35 * L, -1)[0], run.at(0.35 * L, 1)[0]
    assert s_after[1] - s_before[1] == pytest.approx(4000.0, rel=1e-12)           # V drops by the load (-4000)
    assert run.at(0.7 * L, 1)[0][5] - run.at(0.7 * L, -1)[0][5] == pytest.approx(-2.5e6, rel=1e-9)


@pytest.mark.parametrize("shear", [True, False])
def test_integrator_reproduces_the_ode_solution_of_a_cantilever_in_the_xz_plane(shear):
    loads = SpanLoads(points=((0.5 * L, f6(Fz=3000.0)), (0.25 * L, f6(My=-1.5e6))),
                      dists=((0.2 * L, 0.6 * L, w3(z=4.0), w3(z=-1.0)),))
    run = integrate(L, FLEX if shear else FLEX_EB, loads, loads.resultant(), np.zeros(6))
    ref = beam_bvp(L, EIy, GAz if shear else None, clamped, free, [(0.2 * L, 0.6 * L, 4.0, -1.0)],
                   [(0.5 * L, 3000.0, 0.0), (0.25 * L, 0.0, -1.5e6)], plane="xz")
    sd, sm, sv = abs(ref(L)[0]), abs(ref(0.0, 1)[2]), 3000.0
    for x in np.linspace(0, L, 13):
        for side in (-1, 1):
            s, y = run.at(float(x), side)
            d, th, M, V = ref(x, side)
            close([y[2]], [d], scale=sd)
            close([y[4]], [th], scale=abs(ref(L)[1]))
            close([s[4]], [M], scale=sm)
            close([s[2]], [V], scale=sv)


def test_integrator_axial_and_torsion():
    loads = SpanLoads(points=((0.4 * L, f6(Fx=900.0, Mx=4.0e5)),), dists=((0.0, 0.5 * L, w3(x=1.5), w3(x=0.5)),))
    run = integrate(L, FLEX, loads, loads.resultant(), np.zeros(6))
    ax = bar_bvp(L, EA, lambda y: [y[0]], lambda y: [y[1]], [(0.0, 0.5 * L, 1.5, 0.5)], [(0.4 * L, 900.0)])
    to = bar_bvp(L, GJ, lambda y: [y[0]], lambda y: [y[1]], (), [(0.4 * L, 4.0e5)])
    for x in np.linspace(0, L, 11):
        for side in (-1, 1):
            s, y = run.at(float(x), side)
            close([y[0]], [ax(x, side)[0]], scale=abs(ax(L)[0]))
            close([s[0]], [ax(x, side)[1]], scale=abs(ax(0.0, 1)[1]))
            close([y[3]], [to(x, side)[0]], scale=abs(to(L)[0]))
            close([s[3]], [to(x, side)[1]], scale=4.0e5)


def test_resultant_of_span_loads_equals_numerical_integration():
    from scipy.integrate import quad
    x1, x2, a, b = 0.15 * L, 0.85 * L, w3(1.0, -2.0, 3.0), w3(-0.5, 4.0, 1.0)
    r = SpanLoads(dists=((x1, x2, a, b),)).resultant()
    w = lambda x: a + (b - a) * (x - x1) / (x2 - x1)
    for c in range(3):
        assert r[c] == pytest.approx(quad(lambda x: w(x)[c], x1, x2)[0], rel=1e-12)
    assert r[4] == pytest.approx(quad(lambda x: -x * w(x)[2], x1, x2)[0], rel=1e-12)      # (x e_x) x (w_z e_z)
    assert r[5] == pytest.approx(quad(lambda x: x * w(x)[1], x1, x2)[0], rel=1e-12)
    assert r[3] == 0.0


# ---------------------------------------------------------------- flexibility and stiffness
@pytest.mark.parametrize("flex", [FLEX, FLEX_EB])
def test_closed_form_cantilever_flexibility_equals_exact_integration_of_unit_end_loads(flex):
    Fb = cantilever_flexibility(L, flex)
    for a in range(6):
        q = np.zeros(6)
        q[a] = 1.0
        loads = SpanLoads(points=((L, q),))
        run = integrate(L, flex, loads, loads.resultant(), np.zeros(6))
        np.testing.assert_allclose(run.y_end, Fb[:, a], rtol=1e-13, atol=1e-30)
    np.testing.assert_array_equal(Fb, Fb.T)


@pytest.mark.parametrize("shear", [True, False])
def test_cantilever_flexibility_against_the_independent_ode_solution(shear):
    Fb = cantilever_flexibility(L, FLEX if shear else FLEX_EB)
    # unit end force / end moment enter the BVP as boundary conditions (resultants at x = L- equal the end loads)
    y = beam_bvp(L, EIz, GAy if shear else None, clamped, lambda s: [s[2], s[3] - 1.0])(L)
    close([Fb[1, 1], Fb[5, 1]], y[:2])
    y = beam_bvp(L, EIz, GAy if shear else None, clamped, lambda s: [s[2] - 1.0, s[3]])(L)
    close([Fb[1, 5], Fb[5, 5]], y[:2])
    y = beam_bvp(L, EIy, GAz if shear else None, clamped, lambda s: [s[2], s[3] - 1.0], plane="xz")(L)
    close([Fb[2, 2], Fb[4, 2]], y[:2])
    y = beam_bvp(L, EIy, GAz if shear else None, clamped, lambda s: [s[2] - 1.0, s[3]], plane="xz")(L)
    close([Fb[2, 4], Fb[4, 4]], y[:2])
    assert Fb[0, 0] == pytest.approx(bar_bvp(L, EA, lambda s: [s[0]], lambda s: [s[1] - 1.0])(L)[0], rel=1e-9)
    assert Fb[3, 3] == pytest.approx(bar_bvp(L, GJ, lambda s: [s[0]], lambda s: [s[1] - 1.0])(L)[0], rel=1e-9)


def local_k(flex=FLEX, k_i=RR, k_j=RR):
    return release_stiffness(cantilever_flexibility(L, flex), rigid_transfer([L, 0, 0]), k_i, k_j)[0]


def scaled(K):
    d = np.array([1, 1, 1, L, L, L] * 2, float)
    return K / np.outer(d, d)


@pytest.mark.parametrize("shear", [True, False])
def test_stiffness_columns_are_the_end_forces_of_the_ode_solution_with_unit_end_displacements(shear):
    K = local_k(FLEX if shear else FLEX_EB)
    GA_y, GA_z = (GAy, GAz) if shear else (None, None)
    # end forces on the member: at i  Fy = -V(0), Mz = -M(0);  at j  Fy = V(L), Mz = M(L)  (CONVENTIONS section 4)
    for plane, EI, GA, (iv, ith, jv, jth) in (("xy", EIz, GA_y, (1, 5, 7, 11)), ("xz", EIy, GA_z, (2, 4, 8, 10))):
        for col, bc0, bcL in ((iv, lambda s: [s[0] - 1.0, s[1]], clamped), (ith, lambda s: [s[0], s[1] - 1.0], clamped),
                              (jv, clamped, lambda s: [s[0] - 1.0, s[1]]), (jth, clamped, lambda s: [s[0], s[1] - 1.0])):
            ref = beam_bvp(L, EI, GA, bc0, bcL, plane=plane)
            a, b = ref(0.0, 1), ref(L)
            expect = np.array([-a[3], -a[2], b[3], b[2]])
            close(K[[iv, ith, jv, jth], col], expect, rel=1e-9)
    for col, stiff, (i, j) in ((0, EA, (0, 6)), (6, EA, (0, 6)), (3, GJ, (3, 9)), (9, GJ, (3, 9))):
        bc0 = (lambda s: [s[0] - 1.0]) if col == i else (lambda s: [s[0]])
        bcL = (lambda s: [s[0] - 1.0]) if col == j else (lambda s: [s[0]])
        ref = bar_bvp(L, stiff, bc0, bcL)
        close(K[[i, j], col], [-ref(0.0, 1)[1], ref(L)[1]])


@pytest.mark.parametrize("k_i,k_j", [
    (RR, RR),
    ((INF, INF, INF, INF, 0.0, 0.0), RR),
    ((INF, INF, INF, INF, INF, 3.0e9), (INF, INF, INF, INF, 2.0e8, 0.0)),
    ((INF, 0.0, INF, INF, INF, INF), (5.0e4, INF, 700.0, 0.0, INF, INF)),
    ((1.0e5, 2.0e3, 3.0e3, 4.0e8, 5.0e9, 6.0e9), (7.0e5, 8.0e3, 9.0e3, 1.0e8, 2.0e9, 3.0e9)),
])
def test_stiffness_is_symmetric_positive_semidefinite_with_exactly_six_rigid_body_modes(k_i, k_j):
    """Zero-energy motions of the 12 end DOFs: the 6 rigid-body motions, plus one per fully released
    connection (the node side of a k = 0 connection moves freely) -- and nothing else."""
    K = local_k(FLEX, k_i, k_j)
    n_zero = 6 + sum(k == 0.0 for k in k_i + k_j)
    Ks = scaled(K)
    assert np.abs(Ks - Ks.T).max() <= 1e-12 * np.abs(Ks).max()
    lam = np.linalg.eigvalsh(0.5 * (Ks + Ks.T))
    assert np.sum(np.abs(lam) <= 1e-10 * lam.max()) == n_zero and lam.min() > -1e-10 * lam.max()
    B = np.zeros((12, 6))                                   # rigid-body motions of the two ends
    for a in range(3):
        e = np.eye(3)[a]
        B[a, a] = B[6 + a, a] = 1.0
        B[3 + a, 3 + a] = B[9 + a, 3 + a] = 1.0
        B[6:9, 3 + a] = np.cross(e, [L, 0, 0])
    assert np.abs(K @ B).max() <= 1e-9 * (np.abs(K) @ np.abs(B)).max()


# -------------------------------------------------------------------------------- releases
def test_release_limits_k_to_zero_and_k_to_infinity_are_reached_continuously_and_without_ill_conditioning():
    rigid = local_k(FLEX, RR, RR)
    pinned_j = local_k(FLEX, RR, (INF, INF, INF, INF, INF, 0.0))
    EIL = EIz / L
    prev = None
    for k in [0.0, 1e-300, 1e-30, 1e-12, 1e-6, 1e-3, 0.1, 1.0, 10.0, 1e3, 1e6, 1e12, 1e30, 1e300]:
        K = local_k(FLEX, RR, (INF, INF, INF, INF, INF, k * EIL))
        assert np.isfinite(K).all()
        Ks = scaled(K)
        assert np.abs(Ks - Ks.T).max() <= 1e-12 * np.abs(Ks).max()
        if k <= 1e-12:
            assert np.abs(scaled(K - pinned_j)).max() <= 1e-11 * np.abs(scaled(rigid)).max()
        if k >= 1e12:
            assert np.abs(scaled(K - rigid)).max() <= 1e-11 * np.abs(scaled(rigid)).max()
        m_end = K[11, 11]                                   # moment at j per unit rotation of j: grows with k
        if prev is not None:
            assert m_end >= prev - 1e-9 * rigid[11, 11]
        prev = m_end
    assert pinned_j[11, 11] == pytest.approx(0.0, abs=1e-9 * rigid[11, 11])


@pytest.mark.parametrize("kfac", [1e-9, 1e-4, 1e-2, 0.3, 1.0, 4.0, 1e2, 1e5, 1e10])
@pytest.mark.parametrize("end", ["i", "j"])
def test_partial_rotational_release_against_the_ode_solution_with_a_spring_boundary_condition(kfac, end):
    """Unit rotation of the node at the sprung end, the other end clamped.  The spring sits between node
    and beam: moment through the spring (on the +x face) = k * (rotation right of it - rotation left of it)."""
    k = kfac * EIz / L
    if end == "j":
        K = local_k(FLEX, RR, (INF, INF, INF, INF, INF, k))
        ref = beam_bvp(L, EIz, GAy, clamped, lambda s: [s[0], s[2] - k * (1.0 - s[1])])
        col = 11
    else:
        K = local_k(FLEX, (INF, INF, INF, INF, INF, k), RR)
        ref = beam_bvp(L, EIz, GAy, lambda s: [s[0], s[2] - k * (s[1] - 1.0)], clamped)
        col = 5
    a, b = ref(0.0, 1), ref(L)
    close(K[[1, 5, 7, 11], col], [-a[3], -a[2], b[3], b[2]], rel=1e-9)


@pytest.mark.parametrize("kfac", [1e-8, 1e-2, 1.0, 50.0, 1e8])
def test_partial_translational_release_against_the_ode_solution(kfac):
    """Shear spring at end j (uy): shear through the spring = k * (v_node - v_beam)."""
    k = kfac * EIz / L ** 3
    K = local_k(FLEX, RR, (INF, k, INF, INF, INF, INF))
    ref = beam_bvp(L, EIz, GAy, clamped, lambda s: [s[3] - k * (1.0 - s[0]), s[1]])
    a, b = ref(0.0, 1), ref(L)
    close(K[[1, 5, 7, 11], 7], [-a[3], -a[2], b[3], b[2]], rel=1e-9)
    ka = kfac * EA / L                                       # axial spring at end i
    K = local_k(FLEX, (ka, INF, INF, INF, INF, INF), RR)
    ref = bar_bvp(L, EA, lambda s: [s[1] - ka * (s[0] - 1.0)], lambda s: [s[0]])
    close(K[[0, 6], 0], [-ref(0.0, 1)[1], ref(L)[1]], rel=1e-9)


def beam_member(**member):
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": {"i": "a", "j": "b", **member}}, {"a": "fixed", "b": "fixed"},
               {"c": [{"type": "nodal", "node": "a", "Fx": 1.0}]})
    model = build_model(doc)
    return model, model.members["m"]


@pytest.mark.parametrize("k_i,k_j,expected", [
    ((INF, INF, INF, 0.0, INF, INF), (INF, INF, INF, 0.0, INF, INF), {"i.rx", "j.rx"}),
    ((0.0, INF, INF, INF, INF, INF), (0.0, INF, INF, INF, INF, INF), {"i.ux", "j.ux"}),
    ((INF, 0.0, INF, INF, INF, INF), (INF, 0.0, INF, INF, INF, INF), {"i.uy", "j.uy"}),
    ((INF, INF, INF, INF, INF, 0.0), (INF, 0.0, INF, INF, INF, 0.0), {"i.rz", "j.uy", "j.rz"}),
    ((INF, INF, 0.0, INF, 0.0, INF), (INF, INF, INF, INF, 0.0, INF), {"i.uz", "i.ry", "j.ry"}),
])
def test_member_level_mechanism_is_refused_with_the_releases_involved(k_i, k_j, expected):
    model, m = beam_member()
    bad = dataclasses.replace(m, k_end_i=k_i, k_end_j=k_j)
    with pytest.raises(ModelRefused) as e:
        build_element(model, bad)
    p = e.value.problems[0]
    assert p.code == "member.release-mechanism" and p.path == "/members/m/releases"
    named = set(p.message.split("releases involved: ")[1].split(", "))
    assert named == expected


def test_a_nearly_free_spring_pair_is_refused_as_ill_conditioned_not_solved():
    model, m = beam_member()
    tiny = 1e-20
    bad = dataclasses.replace(m, k_end_i=(INF, INF, INF, tiny, INF, INF), k_end_j=(INF, INF, INF, tiny, INF, INF))
    with pytest.raises(ModelRefused):
        build_element(model, bad)
    ok = dataclasses.replace(m, k_end_i=(INF, INF, INF, tiny, INF, INF))          # one tiny spring: fine
    el = build_element(model, ok)
    assert el.k_local[3, 3] == pytest.approx(tiny, rel=1e-9)


# ----------------------------------------------------------------------- fixed-end forces
LOADS_XY = {
    "uniform": (SpanLoads(dists=((0.0, L, w3(y=-3.0), w3(y=-3.0)),)), [(0.0, L, -3.0, -3.0)], []),
    "partial trapezoid": (SpanLoads(dists=((0.2 * L, 0.75 * L, w3(y=-1.0), w3(y=-5.0)),)), [(0.2 * L, 0.75 * L, -1.0, -5.0)], []),
    "point force": (SpanLoads(points=((0.3 * L, f6(Fy=-8000.0)),)), [], [(0.3 * L, -8000.0, 0.0)]),
    "point moment": (SpanLoads(points=((0.6 * L, f6(Mz=3.0e6)),)), [], [(0.6 * L, 0.0, 3.0e6)]),
}


@pytest.mark.parametrize("name", list(LOADS_XY))
@pytest.mark.parametrize("config", ["fixed-fixed", "hinge at i", "hinge at j", "spring at j", "shear release at j",
                                    "no shear deformation"])
def test_fixed_end_forces_against_the_ode_solution(name, config):
    loads, dist, points = LOADS_XY[name]
    model, m = beam_member()
    k_i, k_j, GA, bc0, bcL = RR, RR, GAy, clamped, clamped
    ks = 0.8 * EIz / L
    if config == "hinge at i":
        k_i, bc0 = (INF, INF, INF, INF, INF, 0.0), (lambda s: [s[0], s[2]])
    elif config == "hinge at j":
        k_j, bcL = (INF, INF, INF, INF, INF, 0.0), (lambda s: [s[0], s[2]])
    elif config == "spring at j":
        k_j, bcL = (INF, INF, INF, INF, INF, ks), (lambda s: [s[0], s[2] - ks * (0.0 - s[1])])
    elif config == "shear release at j":
        k_j, bcL = (INF, 0.0, INF, INF, INF, INF), (lambda s: [s[1], s[3]])
    elif config == "no shear deformation":
        GA = None
        model.analysis = dataclasses.replace(model.analysis, shear_deformation=False)
    el = build_element(model, dataclasses.replace(m, k_end_i=k_i, k_end_j=k_j))
    f = el.fixed_end(loads)
    ref = beam_bvp(L, EIz, GA, bc0, bcL, dist, points)
    a, b = ref(0.0, 1), ref(L)
    close(f[[1, 5, 7, 11]], [-a[3], -a[2], b[3], b[2]], rel=1e-9)
    assert np.abs(f[[0, 2, 3, 4, 6, 8, 9, 10]]).max() == 0.0
    # recovery with zero node displacements: stations follow the ODE solution, closures vanish
    rec = el.recover(np.zeros(12), loads, [(x, -1) for x in np.linspace(0, L, 9)])
    for k, x in enumerate(rec.x):
        d, th, M, V = ref(x, -1) if x > 0 else ref(x, 1)
        close([rec.forces[k, 1]], [V], scale=max(abs(a[3]), abs(b[3]), abs(a[2]) / L), rel=1e-9)  # floor: V = 0 everywhere for a moment load on a shear-released member
        close([rec.forces[k, 5]], [M], scale=max(abs(a[2]), abs(b[2]), abs(ref(L / 2)[2])), rel=1e-9)
        close([rec.displacements[k, 1]], [d], scale=max(abs(ref(xx)[0]) for xx in np.linspace(0, L, 9)), rel=1e-9)
    assert np.abs(rec.force_closure).max() <= 1e-9 * np.abs(rec.forces).max()
    assert np.abs(rec.disp_closure[:3]).max() <= 1e-9 * np.abs(rec.displacements[:, :3]).max()


def test_fixed_end_forces_axial_with_partial_axial_release():
    model, m = beam_member()
    ka = 0.6 * EA / L
    el = build_element(model, dataclasses.replace(m, k_end_j=(ka, INF, INF, INF, INF, INF)))
    loads = SpanLoads(points=((0.25 * L, f6(Fx=5000.0)),), dists=((0.5 * L, L, w3(x=2.0), w3(x=4.0)),))
    f = el.fixed_end(loads)
    ref = bar_bvp(L, EA, lambda s: [s[0]], lambda s: [s[1] - ka * (0.0 - s[0])], [(0.5 * L, L, 2.0, 4.0)], [(0.25 * L, 5000.0)])
    close(f[[0, 6]], [-ref(0.0, 1)[1], ref(L)[1]], rel=1e-9)


def test_equivalent_nodal_loads_balance_the_span_loads_for_any_release_pattern():
    model, m = beam_member()
    loads = SpanLoads(points=((0.3 * L, f6(Fx=10.0, Fy=-20.0, Fz=30.0, Mx=4000.0, My=5000.0, Mz=-6000.0)),),
                      dists=((0.1 * L, 0.8 * L, w3(1.0, 2.0, -3.0), w3(-2.0, 0.5, 1.0)),))
    P0 = loads.resultant()
    for k_i, k_j in [(RR, RR), ((INF, INF, INF, INF, 0.0, 0.0), (INF, INF, INF, 0.0, 0.0, 0.0)),
                     ((INF, 2e3, INF, 5e7, INF, 1e9), (3e4, INF, 0.0, INF, 2e8, INF))]:
        el = build_element(model, dataclasses.replace(m, k_end_i=k_i, k_end_j=k_j))
        f = el.fixed_end(loads)
        total = f[:6] + rigid_transfer([L, 0, 0]).T @ f[6:] + P0          # everything reduced to end i
        assert np.abs(total[:3]).max() <= 1e-9 * np.abs(P0[:3]).max()
        assert np.abs(total[3:]).max() <= 1e-9 * np.abs(P0[3:]).max()
        for a in range(6):                                               # a fully released DOF carries nothing
            if k_i[a] == 0.0:
                assert abs(f[a]) <= 1e-9 * np.abs(P0).max()
            if k_j[a] == 0.0:
                assert abs(f[6 + a]) <= 1e-9 * np.abs(P0).max()


def test_truss_element_is_exactly_axial_and_carries_span_loads_as_a_simply_supported_beam():
    model, m = beam_member(type="truss")
    el = build_element(model, m)
    K = el.k_local
    expect = np.zeros((12, 12))
    expect[np.ix_([0, 6], [0, 6])] = EA / L * np.array([[1, -1], [-1, 1]])
    np.testing.assert_allclose(K, expect, rtol=1e-14, atol=0)
    assert not K[[1, 2, 3, 4, 5, 7, 8, 9, 10, 11]].any()                  # exactly zero rows
    loads = SpanLoads(points=((0.3 * L, f6(Fy=-900.0)),))
    f = el.fixed_end(loads)
    ref = beam_bvp(L, EIz, GAy, lambda s: [s[0], s[2]], lambda s: [s[0], s[2]], (), [(0.3 * L, -900.0, 0.0)])
    close(f[[1, 7]], [-ref(0.0, 1)[3], ref(L)[3]], rel=1e-9)
    assert np.abs(f[[5, 11]]).max() <= 1e-9 * 900.0 * L


def test_rigid_offsets_are_the_exact_small_displacement_kinematics():
    doc = make({"a": (0, 0, 0), "b": (5000, 0, 0)},
               {"m": {"i": "a", "j": "b", "offsets": {"i": [300, 40, -70], "j": [-200, -30, 120]},
                      "orientation": {"roll": 25}}},
               {"a": "fixed", "b": "fixed"}, {"c": [{"type": "nodal", "node": "a", "Fx": 1.0}]})
    model = build_model(doc)
    m = model.members["m"]
    el = build_element(model, m)
    rng = np.random.default_rng(7)
    u = rng.normal(size=12) * np.array([1, 1, 1, 1e-3, 1e-3, 1e-3] * 2)
    ue = el.G @ u
    for blk, off in ((0, m.offset_i), (6, m.offset_j)):
        disp = u[blk:blk + 3] + np.cross(u[blk + 3:blk + 6], off)
        np.testing.assert_allclose(ue[blk:blk + 3], m.R @ disp, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(ue[blk + 3:blk + 6], m.R @ u[blk + 3:blk + 6], rtol=1e-13, atol=1e-16)
    # contragredience: nodal forces do the same work as the end forces
    f_end = rng.normal(size=12)
    assert (el.G.T @ f_end) @ u == pytest.approx(f_end @ ue, rel=1e-12)
    # rigid-body motion of the two nodes: no forces
    B = np.zeros((12, 6))
    xi, xj = model.nodes["a"].xyz, model.nodes["b"].xyz
    for a in range(3):
        e = np.eye(3)[a]
        for blk, x in ((0, xi), (6, xj)):
            B[blk + a, a] = 1.0
            B[blk:blk + 3, 3 + a] = np.cross(e, x)
            B[blk + 3 + a, 3 + a] = 1.0
    assert np.abs(el.k_global @ B).max() <= 1e-9 * (np.abs(el.k_global) @ np.abs(B)).max()


def test_zero_stiffness_section_is_refused():
    model, m = beam_member()
    sec = model.sections["S"]
    model.sections["S"] = dataclasses.replace(sec, props=dataclasses.replace(sec.props, J=0.0))
    with pytest.raises(ModelRefused) as e:
        build_element(model, m)
    assert e.value.codes == ["section.zero-stiffness"]


def test_the_two_independent_ode_references_agree():
    """The Chebyshev-collocation reference used throughout against scipy's adaptive collocation solver."""
    from solve_reference import beam_bvp_scipy
    args = (L, EIz, GAy, clamped, lambda s: [s[0], s[2] + 3.0e8 * s[1]])
    kw = dict(dist=[(0.2 * L, 0.8 * L, -3.0, -7.0)], points=[(0.6 * L, -5000.0, 2.0e6)])
    a, b = beam_bvp(*args, **kw), beam_bvp_scipy(*args, **kw)
    ref = np.array([max(abs(a(x)[k]) for x in np.linspace(0, L, 9)) for k in range(4)])
    for x in np.linspace(0, L, 9):
        for side in (-1, 1):
            assert (np.abs(a(x, side) - b(x, side)) <= 1e-9 * ref).all()
