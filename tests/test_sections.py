"""Parametric section properties against hand-derived closed forms and
independent numerical integration."""
import math

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from benchfea import sections
from benchfea.sections import box, i_shape, pipe, rectangle, rectangle_torsion_constant


# --------------------------------------------------------------- independent polygon integrals
def polygon_integrals(pts):
    """Area, first and second moments of a simple polygon (Green's theorem);
    pts = [(y, z), ...] counter-clockwise.  Returns A, int y, int z, int y^2, int z^2."""
    y, z = np.array(pts, float).T
    y1, z1 = np.roll(y, -1), np.roll(z, -1)
    c = y * z1 - y1 * z
    return (c.sum() / 2, ((y + y1) * c).sum() / 6, ((z + z1) * c).sum() / 6,
            ((y * y + y * y1 + y1 * y1) * c).sum() / 12, ((z * z + z * z1 + z1 * z1) * c).sum() / 12)


def i_quarter_outline(d, bf, tf, tw, r, n=20000):
    """Outline of the quarter of an I-shape in y >= 0, z >= 0 (fillet as a fine polyline)."""
    hw2 = d / 2 - tf
    pts = [(0, 0), (d / 2, 0), (d / 2, bf / 2), (hw2, bf / 2)]       # up the web axis, across the flange
    if r > 0:
        yc, zc = hw2 - r, tw / 2 + r                                 # fillet arc centre
        pts.append((hw2, zc))
        for k in range(1, n):
            a = math.pi / 2 * k / n                                  # from the flange down to the web face
            pts.append((yc + r * math.cos(a), zc - r * math.sin(a)))
        pts.append((yc, tw / 2))
    else:
        pts.append((hw2, tw / 2))
    pts.append((0, tw / 2))
    # orientation: make it counter-clockwise in (y, z)
    if polygon_integrals(pts)[0] < 0:
        pts = pts[::-1]
    return pts


def test_rectangle_closed_forms():
    p = rectangle(300.0, 200.0)                      # h along y, b along z
    assert p.A == 60000.0
    assert p.Iz == pytest.approx(200 * 300**3 / 12) == pytest.approx(4.5e8)
    assert p.Iy == pytest.approx(2.0e8)
    assert p.Sz == pytest.approx(3.0e6) and p.Sy == pytest.approx(2.0e6)
    assert p.Zz == pytest.approx(4.5e6) and p.Zy == pytest.approx(3.0e6)
    assert p.Asy == pytest.approx(50000.0) and p.Asz == pytest.approx(50000.0)
    assert p.Cw == 0.0
    assert p.rz == pytest.approx(300 / math.sqrt(12)) and p.ry == pytest.approx(200 / math.sqrt(12))


def test_pipe_closed_forms():
    p = pipe(100.0, 10.0)
    assert p.A == pytest.approx(900 * math.pi)
    i = math.pi * (100**4 - 80**4) / 64
    assert p.Iy == pytest.approx(i) and p.Iz == pytest.approx(i)
    assert p.J == pytest.approx(2 * i)                       # polar second moment (exact for a circle)
    assert p.Zy == pytest.approx((100**3 - 80**3) / 6) == pytest.approx(81333.3333333)
    assert p.Sz == pytest.approx(i / 50)
    assert p.Asy == pytest.approx(450 * math.pi)


def test_pipe_plastic_modulus_and_inertia_by_numerical_integration():
    D, t = 100.0, 10.0
    n = 400000
    theta = (np.arange(n) + 0.5) * 2 * math.pi / n
    ro, ri = D / 2, D / 2 - t
    # int |y| dA and int y^2 dA over the annulus in polar coordinates (radial integrals in closed form)
    z = (np.abs(np.sin(theta)) * (ro**3 - ri**3) / 3).sum() * 2 * math.pi / n
    i = (np.sin(theta)**2 * (ro**4 - ri**4) / 4).sum() * 2 * math.pi / n
    p = pipe(D, t)
    assert p.Zz == pytest.approx(z, rel=1e-9)
    assert p.Iz == pytest.approx(i, rel=1e-9)


def test_box_closed_forms():
    p = box(200.0, 100.0, 10.0)
    assert p.A == pytest.approx(5600.0)
    assert p.Iz == pytest.approx((100 * 200**3 - 80 * 180**3) / 12) == pytest.approx(27786666.6667)
    assert p.Iy == pytest.approx((200 * 100**3 - 180 * 80**3) / 12) == pytest.approx(8986666.6667)
    assert p.Zz == pytest.approx(352000.0)
    assert p.Zy == pytest.approx((200 * 100**2 - 180 * 80**2) / 4) == pytest.approx(212000.0)
    assert p.Sz == pytest.approx(27786666.6667 / 100) and p.Sy == pytest.approx(8986666.6667 / 50)
    # Bredt: J = 4 Am^2 / (perimeter / t), Am = 190 * 90, perimeter = 2 (190 + 90)
    assert p.J == pytest.approx(4 * (190 * 90)**2 / (2 * (190 + 90) / 10)) == pytest.approx(20886428.5714)
    assert p.Asy == 4000.0 and p.Asz == 2000.0 and p.Cw == 0.0
    sq = box(150.0, 150.0, 8.0)
    assert sq.Iy == pytest.approx(sq.Iz) and sq.Zy == pytest.approx(sq.Zz)


def test_i_shape_without_fillets_closed_forms():
    p = i_shape(300.0, 150.0, 10.0, 6.0)
    assert p.A == pytest.approx(4680.0)
    assert p.Iz == pytest.approx(74076000.0)                 # (150*300^3 - 144*280^3)/12
    assert p.Iy == pytest.approx(5630040.0)                  # (2*10*150^3 + 280*6^3)/12
    assert p.Zz == pytest.approx(552600.0)                   # 150*10*290 + 6*280^2/4
    assert p.Zy == pytest.approx(115020.0)                   # 10*150^2/2 + 280*6^2/4
    assert p.Sz == pytest.approx(74076000.0 / 150) and p.Sy == pytest.approx(5630040.0 / 75)
    assert p.J == pytest.approx(120880.0)                    # (2*150*10^3 + 290*6^3)/3
    assert p.Cw == pytest.approx(1.18265625e11)              # 10*150^3*290^2/24
    assert p.Asy == pytest.approx(1800.0) and p.Asz == pytest.approx(2500.0)
    assert p.dims == {"d": 300.0, "bf": 150.0, "tf": 10.0, "tw": 6.0, "r": 0.0}


def test_i_shape_warping_constant_is_flange_inertia_times_h0_squared_over_4():
    d, bf, tf, tw = 14.0, 10.0, 0.7, 0.44
    p = i_shape(d, bf, tf, tw)
    iy_flanges = 2 * tf * bf**3 / 12
    assert p.Cw == pytest.approx(iy_flanges * (d - tf)**2 / 4)


@pytest.mark.parametrize("d,bf,tf,tw,r", [(300.0, 150.0, 12.0, 8.0, 15.0), (14.0, 10.0, 0.7, 0.44, 0.5),
                                          (200.0, 200.0, 15.0, 9.0, 18.0), (300.0, 150.0, 12.0, 8.0, 0.0)])
def test_i_shape_with_fillets_against_polygon_integration(d, bf, tf, tw, r):
    a, qy, qz, iyy, izz = polygon_integrals(i_quarter_outline(d, bf, tf, tw, r))
    p = i_shape(d, bf, tf, tw, r)
    assert p.A == pytest.approx(4 * a, rel=1e-8)
    assert p.Iz == pytest.approx(4 * iyy, rel=1e-8)          # int y^2 dA
    assert p.Iy == pytest.approx(4 * izz, rel=1e-8)          # int z^2 dA
    assert p.Zz == pytest.approx(4 * qy, rel=1e-8)           # int |y| dA
    assert p.Zy == pytest.approx(4 * qz, rel=1e-8)           # int |z| dA
    if r > 0:
        plain = i_shape(d, bf, tf, tw)
        assert p.A - plain.A == pytest.approx(4 * r * r * (1 - math.pi / 4))
        assert p.Iz > plain.Iz and p.Iy > plain.Iy
        assert (p.J, p.Cw, p.Asy) == (plain.J, plain.Cw, plain.Asy)     # documented: fillets ignored


def prandtl_torsion_constant(a, b, n):
    """J of an a x b rectangle from a finite-difference solution of
    laplace(phi) = -2, phi = 0 on the boundary; J = 2 * integral(phi)."""
    nx, ny = n, max(2, round(n * b / a))
    hx, hy = a / nx, b / ny
    mx, my = nx - 1, ny - 1

    def lap(m, h):
        return sp.diags([np.full(m - 1, 1.0), np.full(m, -2.0), np.full(m - 1, 1.0)], [-1, 0, 1]) / h**2

    op = sp.kron(lap(mx, hx), sp.identity(my)) + sp.kron(sp.identity(mx), lap(my, hy))
    phi = spla.spsolve(op.tocsc(), np.full(mx * my, -2.0))
    return 2.0 * phi.sum() * hx * hy


@pytest.mark.parametrize("a,b", [(1.0, 1.0), (2.0, 1.0), (4.0, 1.0)])
def test_rectangle_torsion_constant_against_finite_differences(a, b):
    coarse, fine = prandtl_torsion_constant(a, b, 60), prandtl_torsion_constant(a, b, 120)
    richardson = fine + (fine - coarse) / 3.0                # second-order scheme
    assert rectangle_torsion_constant(a, b) == pytest.approx(richardson, rel=2e-4)
    assert rectangle_torsion_constant(b, a) == pytest.approx(rectangle_torsion_constant(a, b), rel=1e-14)


def test_rectangle_torsion_limits():
    # thin strip: J -> h b^3 / 3 from below; and J never exceeds the polar second moment
    assert 0.99 < rectangle_torsion_constant(1000.0, 1.0) / (1000.0 / 3) < 1.0
    for h, b in ((1.0, 1.0), (3.0, 1.0), (1.0, 7.0)):
        assert rectangle_torsion_constant(h, b) < (h * b**3 + b * h**3) / 12


def test_rectangle_shear_area_is_the_energy_value_for_parabolic_shear_stress():
    # tau(y) = V Q(y) / (I b); strain energy  int tau^2 / (2G) dA  ==  V^2 / (2 G As)
    h, b, V = 0.3, 0.2, 1.0
    n = 200000
    y = (np.arange(n) + 0.5) / n * h - h / 2
    q = b * (h * h / 4 - y * y) / 2
    tau = V * q / ((b * h**3 / 12) * b)
    As = V * V / ((tau**2).sum() * (h / n) * b)
    assert rectangle(h, b).Asy == pytest.approx(As, rel=1e-8)


def test_pipe_shear_area_is_the_thin_walled_energy_value():
    # thin ring: shear flow q = V sin(theta) / (pi R); energy int q^2/(2 G t) ds == V^2 / (2 G As)
    R, t, V = 100.0, 0.01, 1.0
    n = 100000
    th = (np.arange(n) + 0.5) * 2 * math.pi / n
    q = V * np.sin(th) / (math.pi * R)
    As = V * V / ((q**2 / t).sum() * R * 2 * math.pi / n)
    assert pipe(2 * R + t, t).Asy == pytest.approx(As, rel=1e-3)      # thin-wall limit


def test_overrides_and_definition_dispatch():
    p = sections.from_definition({"type": "I", "d": 18, "bf": 7.5, "tf": 0.57, "tw": 0.355, "overrides": {"J": 1.25}})
    assert p.J == 1.25 and p.overridden == ("J",) and p.family == "I"
    assert p.Iz == i_shape(18, 7.5, 0.57, 0.355).Iz
    g = sections.from_definition({"type": "properties", "A": 10, "Iy": 50, "Iz": 200, "J": 1.5, "Asy": 4})
    assert (g.family, g.A, g.Iy, g.Iz, g.J, g.Asy, g.Asz, g.Cw, g.Zz) == ("general", 10, 50, 200, 1.5, 4, None, 0.0, None)
    assert sections.from_definition({"type": "box", "H": 8, "B": 6, "t": 0.375}).family == "box"
    assert sections.from_definition({"type": "pipe", "D": 6, "t": 0.28}).family == "pipe"
    assert sections.from_definition({"type": "rectangle", "h": 2, "b": 1}).family == "rectangle"
    with pytest.raises(ValueError):
        sections.from_definition({"type": "standard", "shape": "X"})


def test_every_parametric_shape_lists_its_approximations():
    for p in (i_shape(300, 150, 10, 6), box(200, 100, 10), pipe(100, 10), rectangle(300, 200)):
        assert p.approximations


def test_dimension_checks():
    assert sections.check_dimensions({"type": "I", "d": 300, "bf": 150, "tf": 10, "tw": 6, "r": 12}) == []
    bad = sections.check_dimensions({"type": "I", "d": 300, "bf": 150, "tf": 150, "tw": 200, "r": 100})
    assert [k for k, _ in bad] == ["tf", "tw", "r"]
    assert [k for k, _ in sections.check_dimensions({"type": "box", "H": 100, "B": 50, "t": 25})] == ["t"]
    assert [k for k, _ in sections.check_dimensions({"type": "pipe", "D": 100, "t": 50})] == ["t"]
    assert sections.check_dimensions({"type": "rectangle", "h": 1, "b": 2}) == []
