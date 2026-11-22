"""Fixed-end forces and the parallel-member decomposition of benchfea.crosscheck.fef, verified against
direct integration of the governing ODE (benchfea.oracle) -- a different numerical path."""
import math

import numpy as np
import pytest

from benchfea.crosscheck.fef import bar_fef, bending_split, clip_loads, fixed_end_forces, plane_fef
from benchfea.oracle import Bar, Beam

L, EI, GA = 4000.0, 8.0e12, 3.0e8
INF = math.inf
POINTS = [(1200.0, -5e3, 2e6), (3100.0, 7e3, 0.0)]                  # (x, force, moment)
DISTS = [(500.0, 2900.0, -3.0, -8.0), (0.0, L, 1.5, 1.5)]           # (x1, x2, w1, w2)
CONNECTIONS = [(INF, INF), (INF, 0.0), (0.0, INF), (0.0, 0.0), (3e9, INF), (INF, 2e10), (3e9, 2e10),
               (0.0, 2e10), (3e9, 0.0)]


def _support(k):
    return "fixed" if math.isinf(k) else (None if k == 0.0 else k)


def _oracle(ki, kj, ga, points=POINTS, dists=DISTS):
    """A rotational connection spring between the member end and a held joint is a rotational spring
    to ground."""
    b = Beam(L, EI, ga).support(0.0, v="fixed", r=_support(ki)).support(L, v="fixed", r=_support(kj))
    for x, F, M in points:
        b.point(x, P=F, M=M)
    for d in dists:
        b.dist(*d)
    s = b.solve(tol=1e-11)
    return np.r_[s.reactions[0.0], s.reactions[L]]


@pytest.mark.parametrize("ki, kj", CONNECTIONS)
@pytest.mark.parametrize("ga", [None, GA])
def test_plane_fixed_end_forces_against_ode(ki, kj, ga):
    f = np.array(plane_fef(L, EI, ga, POINTS, DISTS, ki, kj))
    o = _oracle(ki, kj, ga)
    scale = np.array([np.abs(o[[0, 2]]).max(), np.abs(o[[1, 3]]).max() or np.abs(o[[0, 2]]).max() * L] * 2)
    assert np.abs(f - o).max() / scale.max() < 1e-9 and np.all(np.abs(f - o) / scale < 1e-8)
    # statics: forces and moments about x = 0 balance the loads
    sumF = sum(F for _, F, _ in POINTS) + sum(0.5 * (w1 + w2) * (x2 - x1) for x1, x2, w1, w2 in DISTS)
    assert abs(f[0] + f[2] + sumF) < 1e-9 * abs(sumF)


def test_known_values_from_statics_and_symmetry():
    """Uniform load q on a clamped-clamped span: by symmetry both end forces are -qL/2 (on the member,
    opposing the load) and the end moments are equal and opposite.  Simply supported: no moments."""
    q = 2.0
    Fi, Mi, Fj, Mj = plane_fef(L, EI, None, [], [(0.0, L, q, q)])
    assert Fi == pytest.approx(-q * L / 2, rel=1e-12) and Fj == pytest.approx(-q * L / 2, rel=1e-12)
    assert Mi == pytest.approx(-Mj, rel=1e-12) and Mi < 0.0
    Fi, Mi, Fj, Mj = plane_fef(L, EI, GA, [(1000.0, 10.0, 0.0)], [], 0.0, 0.0)
    assert (Fi, Mi, Fj, Mj) == pytest.approx((-7.5, 0.0, -2.5, 0.0), abs=1e-12)


def test_bar_fixed_end_forces_against_ode():
    pts, ds = [(900.0, 4000.0), (3000.0, -1500.0)], [(200.0, 3500.0, 1.0, 4.0)]
    Pi, Pj = bar_fef(L, pts, ds)
    b = Bar(L, 1.0e9).support(0.0).support(L)
    for x, P in pts:
        b.point(x, P)
    for d in ds:
        b.dist(*d)
    s = b.solve(tol=1e-11)
    assert Pi == pytest.approx(s.reactions[0.0], rel=1e-9) and Pj == pytest.approx(s.reactions[L], rel=1e-9)


def test_3d_assembly_and_sign_of_the_xz_plane():
    """fixed_end_forces assembles the two planes, the axial and the torsional part.  The x-z plane
    is checked against the oracle in its own sign convention (plane="xz")."""
    EIy, EIz = 3.0e12, EI
    points = [("force", np.array([2000.0, -5e3, 3e3]), 1200.0), ("moment", np.array([4e5, 2e6, -1e6]), 2500.0)]
    dists = [(np.array([0.5, -3.0, 1.0]), np.array([1.5, -8.0, 4.0]), 500.0, 2900.0)]
    f = fixed_end_forces(L, EIy, EIz, GA, None, points, dists, ky=(INF, 0.0), kz=(2e10, INF))
    # x-y plane: EIz, shear stiffness GA, spring at i
    b = Beam(L, EIz, GA).support(0.0, v="fixed", r=2e10).support(L, v="fixed", r="fixed")
    b.point(1200.0, P=-5e3).point(2500.0, M=-1e6).dist(500.0, 2900.0, -3.0, -8.0)
    s = b.solve(tol=1e-11)
    assert f[[1, 5, 7, 11]] == pytest.approx(np.r_[s.reactions[0.0], s.reactions[L]], rel=1e-8)
    # x-z plane: EIy, no shear deformation, hinge at j
    b = Beam(L, EIy, plane="xz").support(0.0, v="fixed", r="fixed").support(L, v="fixed")
    b.point(1200.0, P=3e3).point(2500.0, M=2e6).dist(500.0, 2900.0, 1.0, 4.0)
    s = b.solve(tol=1e-11)
    assert f[[2, 4, 8]] == pytest.approx(np.r_[s.reactions[0.0], s.reactions[L][0]], rel=1e-8)
    assert f[10] == 0.0
    # axial and torsion
    assert f[0] + f[6] == pytest.approx(-(2000.0 + 0.5 * (0.5 + 1.5) * 2400.0), rel=1e-12)
    assert f[3] + f[9] == pytest.approx(-4e5, rel=1e-12) and f[9] == pytest.approx(-4e5 * 2500.0 / L, rel=1e-12)


def _stiffness(ki, kj, EI_, ga):
    """2x2 stiffness between (v_j, th_j) of joint j and the forces the joint exerts on the member at j,
    joint i held, by the ODE: impose v_j = 1, then th_j = 1 (a spring to a rotated joint = spring to
    ground plus the moment k * th_joint)."""
    K = np.zeros((2, 2))
    for col, (dv, dth) in enumerate(((1.0, 0.0), (0.0, 1.0))):
        b = Beam(L, EI_, ga).support(0.0, v="fixed", r=_support(ki))
        if math.isinf(kj):
            b.support(L, v="fixed", r="fixed", settlement=dv, rotation=dth)
        else:
            b.support(L, v="fixed", r=_support(kj), settlement=dv)
            if kj != 0.0 and dth:
                b.point(L, M=kj * dth)
        s = b.solve(tol=1e-11)
        R, Mr = s.reactions[L]
        if not math.isinf(kj) and kj != 0.0:
            Mr = kj * (dth - s.theta(L))                    # moment the joint puts on the member
        K[:, col] = (R, Mr)
    return K


@pytest.mark.parametrize("ki, kj", CONNECTIONS)
@pytest.mark.parametrize("ga", [None, GA])
def test_parallel_member_decomposition_against_ode(ki, kj, ga):
    """K(member with connection springs) = a K(rigid ends, EI a, GA a) + b K(hinge at j, EI') +
    c K(hinge at i, EI'), with EI' = EI/(1 + 3EI/(GA L^2)) and no shear deformation in the hinged ones."""
    a, b, c = bending_split(L, EI, ga, ki, kj)
    assert min(a, b, c) > -1e-12
    EIh = EI / (1.0 + 3.0 * EI / (ga * L * L)) if ga else EI
    target = _stiffness(ki, kj, EI, ga)
    total = np.zeros((2, 2))
    if a > 0.0:
        total += _stiffness(INF, INF, EI * a, ga * a if ga else None)
    if b > 0.0:
        total += _stiffness(INF, 0.0, EIh * b, None)
    if c > 0.0:
        total += _stiffness(0.0, INF, EIh * c, None)
    scale = np.array([[1.0, L], [L, L * L]]) * (12.0 * EI / L ** 3)
    assert np.abs(target - total).max() < 1e-12 or np.all(np.abs(target - total) / scale < 1e-8)


def test_single_spring_closed_form_and_positivity():
    """Euler-Bernoulli, one spring k: a = 1/(1 + 4EI/(kL)), the partner is hinged at the spring end."""
    k = 5.0e9
    a, b, c = bending_split(L, EI, None, INF, k)
    assert a == pytest.approx(1.0 / (1.0 + 4.0 * EI / (k * L)), rel=1e-12)
    assert b == pytest.approx(1.0 - a, rel=1e-12) and abs(c) < 1e-12
    rng = np.random.default_rng(7)
    for _ in range(3000):
        ki, kj = 10.0 ** rng.uniform(5, 14, 2)
        ga = None if rng.random() < 0.3 else 10.0 ** rng.uniform(6, 11)
        assert min(bending_split(L, EI, ga, ki, kj)) > -1e-9


def test_clip_loads():
    pts = [("force", np.array([0.0, 1.0, 0.0]), 1000.0), ("force", np.array([0.0, 2.0, 0.0]), 2000.0),
           ("moment", np.array([0.0, 0.0, 3.0]), 4000.0)]
    ds = [(np.array([0.0, 2.0, 0.0]), np.array([0.0, 6.0, 0.0]), 1000.0, 3000.0)]
    p1, d1 = clip_loads(pts, ds, 0.0, 2000.0, last=False)
    p2, d2 = clip_loads(pts, ds, 2000.0, 4000.0, last=True)
    assert [x for *_, x in p1] == [1000.0] and [x for *_, x in p2] == [0.0, 2000.0]   # boundary load -> next span
    (wa, wb, xa, xb), = d1
    assert (wa[1], wb[1], xa, xb) == (2.0, 4.0, 1000.0, 2000.0)
    (wa, wb, xa, xb), = d2
    assert (wa[1], wb[1], xa, xb) == (4.0, 6.0, 0.0, 1000.0)
    # resultants are preserved by the split
    total = sum(0.5 * (a[1] + b[1]) * (x2 - x1) for a, b, x1, x2 in d1 + d2)
    assert total == pytest.approx(0.5 * (2.0 + 6.0) * 2000.0)
