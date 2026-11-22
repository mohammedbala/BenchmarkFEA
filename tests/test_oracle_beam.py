"""The ODE oracle (benchfea.oracle) against closed forms.

Every closed form used here is derived in the docstring of its test from statics and the beam
equation ``E I v'''' = w`` (or the equations stated in benchfea/oracle/beam.py); none is quoted from
a reference.
"""
import math

import numpy as np
import pytest

from benchfea.oracle import Bar, Beam, OracleError

L, EI, GA = 3000.0, 8.0e12, 2.0e8
TOL = 2e-9


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def test_cantilever_tip_load_with_shear_deformation():
    """Cantilever clamped at 0, force P at L.  M(x) = P (L - x), V = P.
    th' = M/EI -> th(L) = P L^2 / 2EI;  v' = th + V/GA -> v(L) = P L^3 / 3EI + P L / GA."""
    P = -1000.0
    s = Beam(L, EI, GA).support(0.0, v="fixed", r="fixed").point(L, P=P).solve()
    assert rel(s.v(L), P * L ** 3 / (3 * EI) + P * L / GA) < TOL
    assert rel(s.theta(L), P * L ** 2 / (2 * EI)) < TOL
    assert rel(s.M(0.0, +1), P * L) < TOL and rel(s.Q(0.5 * L), P) < TOL
    R, Mr = s.reactions[0.0]
    assert rel(R, -P) < TOL and rel(Mr, -P * L) < TOL        # support force -P, moment -(L x P)
    assert s.rms_residual < 1e-9


def test_simply_supported_uniform_load():
    """w = -q.  Reactions qL/2; M(x) = q x (L - x)/2 (CONVENTIONS section 4 example);
    integrating EI v'' = M twice with v(0) = v(L) = 0: v(L/2) = -5 q L^4 / (384 EI)."""
    q = 2.5
    s = Beam(L, EI).support(0.0, v="fixed").support(L, v="fixed").dist(0.0, L, -q).solve()
    assert rel(s.v(L / 2), -5 * q * L ** 4 / (384 * EI)) < TOL
    assert rel(s.M(L / 2), q * L ** 2 / 8) < TOL
    assert rel(s.Q(0.0, +1), -q * L / 2) < TOL
    assert rel(s.reactions[0.0][0], q * L / 2) < TOL and rel(s.reactions[L][0], q * L / 2) < TOL


def test_triangular_load_statics_and_equilibrium():
    """Linear load 0 -> w2 on a simply supported span: resultant w2 L/2 at 2L/3, so
    R(0) = w2 L / 6 and R(L) = w2 L / 3 (moments about the supports)."""
    w2 = -4.0
    s = Beam(L, EI, GA).support(0.0, v="fixed").support(L, v="fixed").dist(0.0, L, 0.0, w2).solve()
    assert rel(s.reactions[0.0][0], -w2 * L / 6) < TOL
    assert rel(s.reactions[L][0], -w2 * L / 3) < TOL


def test_two_span_continuous_beam():
    """Two equal spans L, uniform q, three rigid supports.  Remove the middle support: midspan
    deflection of the span 2L is 5 q (2L)^4 / 384EI; a force R at midspan lifts it by R (2L)^3 / 48EI
    (both from EI v'''' = w).  Compatibility: R = 5 q L / 4, the end reactions are 3 q L / 8 and the
    moment over the middle support is 3qL/8 * L - q L^2 / 2 = -q L^2 / 8."""
    q = 2.0
    b = Beam(2 * L, EI).dist(0.0, 2 * L, -q)
    for x in (0.0, L, 2 * L):
        b.support(x, v="fixed")
    s = b.solve()
    assert rel(s.reactions[L][0], 1.25 * q * L) < TOL
    assert rel(s.reactions[0.0][0], 0.375 * q * L) < TOL
    assert rel(s.M(L), -q * L ** 2 / 8) < TOL
    assert abs(s.v(L)) < 1e-9


def test_spring_support_and_settlement():
    """Cantilever with a spring k under the tip and a tip force P: the spring force is k v, so
    v = P / (k + 3EI/L^3).  Imposed settlement d of a rigid tip support on an unloaded cantilever:
    the support force is 3 EI d / L^3."""
    P, k = -5000.0, 400.0
    s = Beam(L, EI).support(0.0, v="fixed", r="fixed").support(L, v=k).point(L, P=P).solve()
    v = P / (k + 3 * EI / L ** 3)
    assert rel(s.v(L), v) < TOL and rel(s.reactions[L][0], -k * v) < TOL
    d = -3.0
    s = Beam(L, EI).support(0.0, v="fixed", r="fixed").support(L, v="fixed", settlement=d).solve()
    assert rel(s.v(L), d) < TOL and rel(s.reactions[L][0], 3 * EI * d / L ** 3) < TOL


def test_rotational_spring_and_imposed_rotation():
    """Beam simply supported at both ends with a rotational spring kr at 0 and a moment Mo applied at L.
    With th0 the rotation at 0: the spring moment on the beam is -kr th0, the beam is a simply
    supported span under end moments.  From EI v'' = M with M linear between the end moments Ma (at 0)
    and Mb (at L):  th(0) = -(2 Ma + Mb) L / 6EI.  Joint relations: Ma = kr th0 (M+ = kr th at the left
    end) and Mb = Mo ... solved for th0:  th0 = -Mo L / (6EI + 2 kr L)."""
    kr, Mo = 4.0e9, 6.0e6
    s = Beam(L, EI).support(0.0, v="fixed", r=kr).support(L, v="fixed").point(L, M=Mo).solve()
    th0 = -Mo * L / (6 * EI + 2 * kr * L)
    assert rel(s.theta(0.0, +1), th0) < TOL
    assert rel(s.M(0.0, +1), kr * th0) < TOL and rel(s.M(L), Mo) < TOL
    # imposed rotation of a clamped end, far end clamped: end moments 4EI th/L and 2EI th/L (from the
    # cubic v with v(0)=v(L)=0, v'(0)=th, v'(L)=0: v'' (0) = -4 th/L, v''(L) = 2 th/L)
    th = 0.002
    s = Beam(L, EI).support(0.0, v="fixed", r="fixed", rotation=th).support(L, v="fixed", r="fixed").solve()
    assert rel(s.M(0.0, +1), -4 * EI * th / L) < TOL and rel(s.M(L), 2 * EI * th / L) < TOL


def test_internal_hinge_free_and_sprung():
    """Clamped at both ends, hinge at midspan, force P at the hinge: by symmetry each half is a
    cantilever carrying P/2, so v = (P/2) (L/2)^3 / 3EI.  With a rotational spring kh in the hinge and
    a NON-symmetric load the rotation jump must equal M/kh."""
    P = -3000.0
    s = Beam(L, EI).support(0.0, v="fixed", r="fixed").support(L, v="fixed", r="fixed") \
        .hinge(L / 2).point(L / 2, P=P).solve()
    assert rel(s.v(L / 2), 0.5 * P * (L / 2) ** 3 / (3 * EI)) < TOL
    assert abs(s.M(L / 2)) < 1e-6 * abs(P * L)
    kh = 3.0e9
    s = Beam(L, EI).support(0.0, v="fixed", r="fixed").support(L, v="fixed", r="fixed") \
        .hinge(0.4 * L, k=kh).point(0.7 * L, P=P).solve()
    jump = s.theta(0.4 * L, +1) - s.theta(0.4 * L, -1)
    assert abs(jump) > 1e-6 and rel(s.M(0.4 * L, -1), kh * jump) < 1e-8
    assert rel(s.M(0.4 * L, +1), s.M(0.4 * L, -1)) < 1e-9
    # kh -> very stiff reproduces the continuous beam
    stiff = Beam(L, EI).support(0.0, v="fixed", r="fixed").support(L, v="fixed", r="fixed") \
        .hinge(0.4 * L, k=1e20).point(0.7 * L, P=P).solve()
    cont = Beam(L, EI).support(0.0, v="fixed", r="fixed").support(L, v="fixed", r="fixed").point(0.7 * L, P=P).solve()
    assert rel(stiff.v(0.7 * L), cont.v(0.7 * L)) < 1e-7


@pytest.mark.parametrize("ratio", [0.2, 0.5, 0.8])
def test_beam_column_cantilever_closed_form(ratio):
    """Cantilever, axial compression P (N = -P), tip force H, no shear deformation.
    M' = -Q + N v' with Q = H:  EI v'' = M = H (L - x) + P (d - v), d = v(L).  With k^2 = P/EI the
    solution satisfying v(0) = v'(0) = 0 gives  d = (H / P k) (tan kL - kL), i.e.
    d = H L^3/3EI * 3 (tan u - u)/u^3 with u = kL, and M(0) = H L + P d = (H/k) tan u."""
    Pe = math.pi ** 2 * EI / (4 * L ** 2)
    P, H = ratio * Pe, 1000.0
    s = Beam(L, EI, N=-P, second_order=True).support(0.0, v="fixed", r="fixed").point(L, P=H).solve()
    u = L * math.sqrt(P / EI)
    assert rel(s.v(L), H * L ** 3 / (3 * EI) * 3 * (math.tan(u) - u) / u ** 3) < 1e-8
    assert rel(s.M(0.0, +1), H * math.tan(u) / math.sqrt(P / EI)) < 1e-8
    assert rel(s.reactions[0.0][0], -H) < 1e-9             # Q is the force along the undeformed axis
    # tension stiffens: same derivation with hyperbolic functions, d = (H / P k)(kL - tanh kL)
    s = Beam(L, EI, N=+P, second_order=True).support(0.0, v="fixed", r="fixed").point(L, P=H).solve()
    k = math.sqrt(P / EI)
    assert rel(s.v(L), H / (P * k) * (k * L - math.tanh(k * L))) < 1e-8


def test_second_order_flag_off_ignores_axial_force():
    s1 = Beam(L, EI, N=-1e6).support(0.0, v="fixed", r="fixed").point(L, P=1.0).solve()
    assert rel(s1.v(L), L ** 3 / (3 * EI)) < TOL


def test_blows_up_at_the_buckling_load():
    """At P = Pe = pi^2 EI / 4L^2 the homogeneous cantilever problem has a nontrivial solution; just
    below it the amplification 3 (tan u - u)/u^3 is huge."""
    Pe = math.pi ** 2 * EI / (4 * L ** 2)
    s = Beam(L, EI, N=-0.999 * Pe, second_order=True).support(0.0, v="fixed", r="fixed").point(L, P=1.0).solve()
    assert s.v(L) / (L ** 3 / (3 * EI)) > 500.0


def test_shear_and_axial_force_needs_an_explicit_shear_model():
    with pytest.raises(OracleError, match="ambiguous"):
        Beam(L, EI, GA, N=-1e5, second_order=True)
    Beam(L, EI, GA, N=0.0, second_order=True)             # no axial force: no ambiguity
    Beam(L, EI, GA, N=-1e5)                               # first order: no ambiguity


def test_engesser_and_haringx_columns():
    """Pin-ended column with shear flexibility under a sinusoidal-free check: both models must reduce
    to the Euler-Bernoulli solution for GA -> infinity, and differ from each other for a soft GA.
    Engesser model, simply supported, uniform w, compression P: substituting v' = (th + Q/GA)/(1 - P/GA)
    into M' = -Q - P v' and th' = M/EI gives  EI v'''' (1 - P/GA) ... = w (1 + ...); instead of that
    algebra the test uses the buckling limit, which follows from the same substitution with w = 0 and
    v = sin(pi x/L):  P_cr,Engesser = Pe/(1 + Pe/GA),  Pe = pi^2 EI/L^2.  Just below P_cr the midspan
    deflection under a small load grows without bound; just above 0.9 P_cr it is still moderate."""
    Pe = math.pi ** 2 * EI / L ** 2
    soft = 5.0 * Pe
    pcr = Pe / (1 + Pe / soft)

    def mid(P, model, ga):
        return Beam(L, EI, ga, N=-P, second_order=True, shear_model=model).support(0.0, v="fixed") \
            .support(L, v="fixed").dist(0.0, L, -1.0).solve().v(L / 2)

    v0 = Beam(L, EI, soft).support(0.0, v="fixed").support(L, v="fixed").dist(0.0, L, -1.0).solve().v(L / 2)
    amp_e = mid(0.99 * pcr, "engesser", soft) / v0
    assert amp_e > 80.0                                    # close to its own critical load
    amp_h = mid(0.99 * pcr, "haringx", soft) / v0
    assert 1.0 < amp_h < 0.5 * amp_e                       # Haringx critical load is higher
    # stiff in shear: both equal the Euler-Bernoulli beam-column
    eb = Beam(L, EI, N=-0.5 * Pe, second_order=True).support(0.0, v="fixed").support(L, v="fixed") \
        .dist(0.0, L, -1.0).solve().v(L / 2)
    for model in ("engesser", "haringx"):
        assert rel(mid(0.5 * Pe, model, 1e9 * Pe), eb) < 1e-7


def test_xz_plane_signs():
    """Cantilever in the x-z plane (CONVENTIONS section 4: My' = +Vz, My = -EIy w'').  Tip force Pz and
    tip moment My: My(0) = My_tip - Pz L (moment about +y of a force Pz at x = L is -L Pz);
    w(L) = Pz L^3/3EI - My L^2/2EI;  th_y(L) = -Pz L^2/2EI + My L/EI."""
    Pz, My = -1000.0, 5.0e5
    s = Beam(L, EI, plane="xz").support(0.0, v="fixed", r="fixed").point(L, P=Pz, M=My).solve()
    assert rel(s.M(0.0, +1), My - Pz * L) < TOL
    assert rel(s.v(L), Pz * L ** 3 / (3 * EI) - My * L ** 2 / (2 * EI)) < TOL
    assert rel(s.theta(L), -Pz * L ** 2 / (2 * EI) + My * L / EI) < TOL
    assert rel(s.reactions[0.0][1], -(My - Pz * L)) < TOL


def test_piecewise_properties():
    """Stepped cantilever (EI1 on [0, a], EI2 beyond), tip force P: by the unit-load theorem
    v(L) = int M m / EI = P [ (L^3 - (L-a)^3) / 3EI1 + (L-a)^3 / 3EI2 ]."""
    a, EI2, P = 1000.0, 2.0e12, 700.0
    s = Beam(L, [(0.0, a, EI), (a, L, EI2)]).support(0.0, v="fixed", r="fixed").point(L, P=P).solve()
    assert rel(s.v(L), P * ((L ** 3 - (L - a) ** 3) / (3 * EI) + (L - a) ** 3 / (3 * EI2))) < TOL


def test_mesh_independence_and_tolerance():
    b = Beam(L, EI, GA).support(0.0, v="fixed", r="fixed").support(L, v="fixed").dist(200.0, 2500.0, -1.0, -6.0)
    b.point(1700.0, P=-2e3, M=3e5)
    # ``tol`` bounds the collocation RESIDUAL, not the error of the solution: the two solutions may
    # differ by a multiple of the looser tolerance (documented in docs/CROSSCHECK.md; bound used: 100 tol)
    a, c = b.solve(tol=1e-8), b.solve(tol=1e-11)
    assert rel(a.v(1500.0), c.v(1500.0)) < 1e-6 and rel(a.reactions[L][0], c.reactions[L][0]) < 1e-6
    d = b.solve(tol=1e-10)
    assert rel(d.v(1500.0), c.v(1500.0)) < 1e-8 and rel(d.reactions[L][0], c.reactions[L][0]) < 1e-8
    total = sum(r[0] for r in c.reactions.values())
    assert rel(total, 0.5 * (1.0 + 6.0) * 2300.0 + 2e3) < 1e-9


def test_ill_posed_problems_are_refused():
    with pytest.raises(OracleError):
        Beam(L, EI).support(0.0, v="fixed").point(L, P=1.0).solve()          # mechanism
    with pytest.raises(OracleError):
        Beam(L, EI).hinge(0.0)
    with pytest.raises(OracleError):
        Beam(L, EI).support(0.0, v="fixed", r="fixed").hinge(1000.0).point(1000.0, M=1.0).solve()
    with pytest.raises(OracleError):
        Beam(L, EI, plane="yz")
    with pytest.raises(OracleError):
        Beam(L, [(0.0, 1000.0, EI)])


def test_bar_axial_and_torsion():
    """Bar fixed at 0, free at L, uniform p: N(x) = p (L - x), u(L) = p L^2 / 2EA.  Fixed-fixed bar with
    a force P at a: N jumps by P; compatibility (N1 a + N2 (L - a) = 0, N1 - N2 = P) gives
    N1 = P (L - a)/L left of the load and N2 = -P a/L right of it; reactions -N1... on the bar:
    R0 = -P (L - a)/L, RL = -P a/L."""
    EA, p = 1.0e9, 3.0
    s = Bar(L, EA).support(0.0).dist(0.0, L, p).solve()
    assert rel(s.u(L), p * L ** 2 / (2 * EA)) < TOL and rel(s.N(0.0, +1), p * L) < TOL
    assert rel(s.reactions[0.0], -p * L) < TOL
    P, a = 5000.0, 1000.0
    s = Bar(L, EA).support(0.0).support(L).point(a, P).solve()
    assert rel(s.N(a, -1), P * (L - a) / L) < TOL and rel(s.N(a, +1), -P * a / L) < TOL
    assert rel(s.reactions[0.0], -P * (L - a) / L) < TOL and rel(s.reactions[L], -P * a / L) < TOL
    # spring end: u = P / (k + EA/L)
    k = 2.0e5
    s = Bar(L, EA).support(0.0).support(L, u=k).point(L, P).solve()
    assert rel(s.u(L), P / (k + EA / L)) < TOL
