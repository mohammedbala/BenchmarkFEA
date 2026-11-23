"""Linear buckling against INDEPENDENT solutions of the column equation.

References (``second_order_reference``): the characteristic determinant of EI v'''' + P v'' = 0 built
from the general solution and the boundary conditions of each case and solved numerically for its
roots (no effective-length factor is quoted), and for shear-deformable columns the eigenvalue of the
beam-column equations found by ``scipy.integrate.solve_bvp``.

Boundary conditions used below (N = -P, k^2 = P / EI; V = force along the undeformed transverse axis
on the +x face = -EI v''' - P v', M = EI v''):
    held translation   v = 0            free translation   v''' + k^2 v' = 0   (V = 0)
    held rotation      v' = 0           free rotation      v'' = 0             (M = 0)

Column a-b along +Z, length L; local y = X.  The section is ten times stiffer about local y, so the
lowest modes are in the X-Z plane (local xy plane, EI = E Iz).  A unit load is applied, so the
buckling factor is the critical load.
"""
import math

import numpy as np
import pytest

from benchfea.model import build_model
from benchfea.solve import analyze
from benchfea.solve.buckling import BUCKLING_TOL

from second_order_reference import column_shape, critical_k, critical_load_bvp
from solve_reference import E, G, make

SEC = {"A": 5000.0, "Iy": 4.0e8, "Iz": 4.0e7, "J": 3.0e5, "Asy": 2500.0, "Asz": 3000.0}
EI, GA, EA = E * SEC["Iz"], G * SEC["Asy"], E * SEC["A"]
L = 4000.0
ACC = 3.0 * BUCKLING_TOL          # stated accuracy: the factors change by <= BUCKLING_TOL between subdivisions

PINNED = {"type": "pinned", "dofs": {"rz": "fixed"}}
ROLLER = {"dofs": {"ux": "fixed", "uy": "fixed"}}                       # lateral roller, free to shorten
CLAMPED_ROLLER = {"dofs": {"ux": "fixed", "uy": "fixed", "rx": "fixed", "ry": "fixed"}}
IN_PLANE = {"dofs": {"uy": "fixed", "rx": "fixed", "rz": "fixed"}}      # keeps a node in the X-Z plane

# rows of the characteristic matrix for each end condition (see the module docstring)
HELD = lambda k, x: column_shape(k, x)[0]                                # noqa: E731
NO_ROTATION = lambda k, x: column_shape(k, x)[1]                         # noqa: E731
NO_MOMENT = lambda k, x: column_shape(k, x)[2]                           # noqa: E731
NO_SHEAR = lambda k, x: column_shape(k, x)[3] + k * k * column_shape(k, x)[1]     # noqa: E731


def model(supports, loads=None, modes=1, shear=False, members=None, nodes=None, sections=None, analysis=None, **kw):
    doc = make(nodes or {"a": (0, 0, 0), "b": (0, 0, L)}, members or {"m": ("a", "b")}, supports,
               {"c": loads or [{"type": "nodal", "node": "b", "Fz": -1.0}]},
               sections=sections or {"S": {"type": "properties", **SEC}},
               analysis={"shearDeformation": shear, "buckling": {"combinations": ["c"], "modes": modes},
                         **(analysis or {})}, **kw)
    return doc


def buckle(doc):
    return analyze(build_model(doc)).buckling["c"]


CASES = {
    "pinned-pinned": ({"a": PINNED, "b": ROLLER}, [HELD, NO_MOMENT], [HELD, NO_MOMENT]),
    "fixed-free": ({"a": "fixed"}, [HELD, NO_ROTATION], [NO_MOMENT, NO_SHEAR]),
    "fixed-pinned": ({"a": "fixed", "b": ROLLER}, [HELD, NO_ROTATION], [HELD, NO_MOMENT]),
    "fixed-fixed": ({"a": "fixed", "b": CLAMPED_ROLLER}, [HELD, NO_ROTATION], [HELD, NO_ROTATION]),
}


def characteristic(case, count=1, length=L):
    _, start, end = CASES[case]
    return critical_k(lambda k: [f(k, 0.0) for f in start] + [f(k, length) for f in end], 12.0 / length, count=count)


@pytest.mark.parametrize("case", list(CASES))
def test_euler_column_four_classic_end_conditions(case):
    k = characteristic(case)
    b = buckle(model(CASES[case][0]))
    assert b.factors[0] == pytest.approx(k * k * EI, rel=ACC)
    assert b.info["estimatedRelativeError"][0] <= BUCKLING_TOL
    # the four cases are genuinely different problems (ratios to the pinned-pinned load: 1, 1/4, ~2.05, 4)
    assert b.factors[0] / (characteristic("pinned-pinned") ** 2 * EI) == pytest.approx(
        {"pinned-pinned": 1.0, "fixed-free": 0.25, "fixed-pinned": 2.0457, "fixed-fixed": 4.0}[case], rel=1e-4)


@pytest.mark.parametrize("case", ["pinned-pinned", "fixed-pinned"])
def test_higher_modes_and_stated_accuracy(case):
    ks = characteristic(case, count=3)
    b = buckle(model(CASES[case][0], modes=3))
    assert b.factors == pytest.approx([k * k * EI for k in ks], rel=ACC)
    assert max(b.info["estimatedRelativeError"]) <= BUCKLING_TOL
    assert list(b.factors) == sorted(b.factors) and len(b.negative_factors) == 0


def test_out_of_plane_mode_uses_the_other_inertia():
    """Hold the weak plane: the first mode is then in the local xz plane with EI = E Iy."""
    b = buckle(model({"a": "fixed", "b": {"dofs": {"ux": "fixed", "ry": "fixed"}}}))
    in_plane = characteristic("fixed-fixed") ** 2 * EI
    out_of_plane = characteristic("fixed-free") ** 2 * E * SEC["Iy"]
    assert out_of_plane < in_plane
    assert b.factors[0] == pytest.approx(out_of_plane, rel=ACC)
    assert abs(b.modes[0][1, 1]) == pytest.approx(1.0) and abs(b.modes[0][1, 0]) < 1e-9      # tip moves along Y


@pytest.mark.parametrize("case,guess,norm", [("pinned-pinned", lambda s: np.sin(np.pi * s), ("start", 1)),
                                             ("fixed-free", lambda s: 1.0 - np.cos(np.pi * s / 2), ("end", 0))])
def test_shear_deformable_column(case, guess, norm):
    """With shear deformation the critical load is lower; reference = eigenvalue of the beam-column
    equations (solve_bvp).  A stocky column is used so that the effect is large (about 20 %)."""
    Ls = 1200.0
    k = characteristic(case, length=Ls)
    euler = k * k * EI
    bcs = {"pinned-pinned": lambda s, e: [s[0], s[2], e[0], e[2]], "fixed-free": lambda s, e: [s[0], s[1], e[2], e[3]]}
    ref = critical_load_bvp(Ls, EI, GA, bcs[case], euler, guess, norm)
    b = buckle(model(CASES[case][0], shear=True, nodes={"a": (0, 0, 0), "b": (0, 0, Ls)}))
    assert b.factors[0] == pytest.approx(ref, rel=ACC)
    assert ref < 0.99 * euler
    no_shear = buckle(model(CASES[case][0], shear=False, nodes={"a": (0, 0, 0), "b": (0, 0, Ls)}))
    assert no_shear.factors[0] == pytest.approx(euler, rel=ACC)


def portal_reference(h, b, EIc, EIb, EAc, base):
    """Sway buckling load per column of a symmetric portal (equal loads on the two column tops).

    In the sway mode both columns have the same shape, so the two column-top shears are equal and,
    summing to zero, vanish: V(h) = 0.  Both joints rotate by theta and the beam bends antisymmetrically.
    Beam (EI v'''' = 0) with end slopes phi relative to its chord: v = phi x (b - x)(b - 2x) / b^2, hence
    end moment 6 EIb phi / b and end shear 12 EIb phi / b^2 (computed below from the polynomial).  The
    shears stretch one column and shorten the other by Q h / EAc each, which rotates the chord of the
    beam by psi = 2 Q h / (EAc b) in the sense of theta: phi = theta - psi.  The beam therefore restrains
    each column top with the rotational stiffness kr = (6 EIb / b) / (1 + kappa), kappa = 24 EIb h / (EAc b^3):
    EIc v''(h) + kr v'(h) = 0."""
    shape = np.poly1d([2.0, -3.0 * b, b * b, 0.0]) / b ** 2            # x (b - x)(b - 2x) / b^2
    assert shape.deriv()(0.0) == pytest.approx(1.0) and shape.deriv()(b) == pytest.approx(1.0)
    m_end = -EIb * shape.deriv(2)(0.0)                                   # per unit phi
    q_end = EIb * shape.deriv(3)(0.0)
    assert m_end == pytest.approx(6 * EIb / b) and q_end == pytest.approx(12 * EIb / b ** 2)
    kappa = 2.0 * q_end * h / (EAc * b)
    kr = m_end / (1.0 + kappa)
    foot = [HELD, NO_ROTATION] if base == "fixed" else [HELD, NO_MOMENT]

    def rows(k):
        return [f(k, 0.0) for f in foot] + [NO_SHEAR(k, h), column_shape(k, h)[2] + kr / EIc * column_shape(k, h)[1]]
    k = critical_k(rows, 6.0 / h)
    return k * k * EIc


@pytest.mark.parametrize("base", ["fixed", "pinned"])
def test_portal_frame_sway_buckling(base):
    h, b = 4000.0, 6000.0
    beam_sec = {"A": 8000.0, "Iy": 9.0e8, "Iz": 9.0e7, "J": 5.0e5}
    doc = model({"a": "fixed" if base == "fixed" else {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}},
                 "d": "fixed" if base == "fixed" else {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}},
                 "b": IN_PLANE, "c": IN_PLANE},
                loads=[{"type": "nodal", "node": "b", "Fz": -1.0}, {"type": "nodal", "node": "c", "Fz": -1.0}],
                nodes={"a": (0, 0, 0), "b": (0, 0, h), "c": (b, 0, h), "d": (b, 0, 0)},
                members={"c1": ("a", "b"), "g": {"i": "b", "j": "c", "section": "B"}, "c2": ("d", "c")},
                sections={"S": {"type": "properties", **SEC}, "B": {"type": "properties", **beam_sec}}, modes=2)
    res = buckle(doc)
    ref = portal_reference(h, b, EI, E * beam_sec["Iz"], EA, base)
    assert res.factors[0] == pytest.approx(ref, rel=ACC)
    mode = res.modes[0]
    assert mode[1, 0] == pytest.approx(1.0) and mode[2, 0] == pytest.approx(1.0, rel=1e-6)   # both tops sway together
    assert mode[1, 4] == pytest.approx(mode[2, 4], rel=1e-6)                                 # equal joint rotations
    assert res.factors[1] > 2.0 * res.factors[0]                                             # the non-sway mode is far above
    # the restraint of the beam matters: the load lies between the two limiting cases of the column alone
    lo = critical_k(lambda k: [HELD(k, 0.0), (NO_ROTATION if base == "fixed" else NO_MOMENT)(k, 0.0),
                               NO_SHEAR(k, h), NO_MOMENT(k, h)], 6.0 / h) ** 2 * EI if base == "fixed" else 0.0
    assert lo < res.factors[0] < characteristic("fixed-fixed", length=h) ** 2 * EI


def test_translational_spring_support():
    """Pinned base, top held by a lateral spring ks: V(L) = -ks v(L), i.e. v''' + k^2 v' - (ks / EI) v = 0."""
    euler = characteristic("pinned-pinned") ** 2 * EI
    for ks, expect_rigid in ((0.5 * euler / L, True), (4.0 * euler / L, False)):
        def rows(k):
            sh = column_shape(k, L)
            return [HELD(k, 0.0), NO_MOMENT(k, 0.0), NO_MOMENT(k, L), sh[3] + k * k * sh[1] - ks / EI * sh[0]]
        ref = critical_k(rows, 12.0 / L) ** 2 * EI
        b = buckle(model({"a": PINNED, "b": {"dofs": {"ux": {"spring": ks}, "uy": "fixed"}}}))
        assert b.factors[0] == pytest.approx(ref, rel=ACC)
        # soft spring: the column tilts as a rigid bar (P = ks L); stiff spring: it bows between its ends
        assert (ref == pytest.approx(ks * L, rel=1e-9)) == expect_rigid
        assert (ref == pytest.approx(euler, rel=1e-9)) != expect_rigid


@pytest.mark.parametrize("how", ["support spring", "partial release", "full release"])
def test_rotational_spring_at_the_base_as_support_and_as_member_release(how):
    """Free top; the base rotation is held by a spring k_theta: EI v''(0) = k_theta v'(0).  The same
    problem is posed with a spring support, with a spring-connected member end at a fixed node, and (k = 0
    with a lateral roller at the top) with a full release."""
    kt = 0.0 if how == "full release" else 3.0 * EI / L
    top = [HELD, NO_MOMENT] if how == "full release" else [NO_MOMENT, NO_SHEAR]

    def rows(k):
        sh = column_shape(k, 0.0)
        return [HELD(k, 0.0), EI * sh[2] - kt * sh[1]] + [f(k, L) for f in top]
    ref = critical_k(rows, 12.0 / L) ** 2 * EI
    if how == "support spring":            # rotation about global Y = local z of the column
        doc = model({"a": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed", "ry": {"spring": kt}}}})
    elif how == "partial release":
        doc = model({"a": "fixed"}, members={"m": {"i": "a", "j": "b", "releases": {"i": {"rz": {"spring": kt}}}}})
    else:
        doc = model({"a": "fixed", "b": ROLLER}, members={"m": {"i": "a", "j": "b", "releases": {"i": {"rz": True}}}})
    b = buckle(doc)
    assert b.factors[0] == pytest.approx(ref, rel=ACC)
    if how == "full release":
        assert ref == pytest.approx(characteristic("pinned-pinned") ** 2 * EI, rel=1e-9)


def test_leaning_column_through_a_rigid_diaphragm():
    """Column 1 is a cantilever; column 2 is pin-ended and leans on it through a rigid diaphragm.  The
    leaner (load P2, straight in the mode) pushes the top of column 1 with P2 v / L:
    V(L) = P2 v(L) / L, i.e. v''' + k^2 v' + (P2 / (EI L)) v = 0 with P2 = 2 P1."""
    def rows(k):
        sh = column_shape(k, L)
        return [HELD(k, 0.0), NO_ROTATION(k, 0.0), NO_MOMENT(k, L), sh[3] + k * k * sh[1] + 2.0 * k * k / L * sh[0]]
    ref = critical_k(rows, 6.0 / L) ** 2 * EI
    doc = model({"a": "fixed", "c": {"type": "pinned", "dofs": {"rz": "fixed"}}, "b": {"dofs": {"uy": "fixed", "rz": "fixed"}}},        # no twist of the floor about the column
                loads=[{"type": "nodal", "node": "b", "Fz": -1.0}, {"type": "nodal", "node": "d", "Fz": -2.0}],
                nodes={"a": (0, 0, 0), "b": (0, 0, L), "c": (3000, 0, 0), "d": (3000, 0, L)},
                members={"m": ("a", "b"), "lean": {"i": "c", "j": "d", "type": "truss"}},
                constraints={"floor": {"type": "rigidDiaphragm", "master": "b", "slaves": ["d"]}})
    b = buckle(doc)
    assert b.factors[0] == pytest.approx(ref, rel=ACC)
    alone = characteristic("fixed-free") ** 2 * EI
    assert alone / 3.0 < ref < alone / 2.0      # worse than the column alone, better than all load on the column
    assert b.modes[0][3, 0] == pytest.approx(b.modes[0][1, 0])                  # the diaphragm ties the tops
    assert any("rigid links / diaphragms" in w for w in b.info["warnings"])


def test_rigid_end_offset_carrying_the_load():
    """The top a = 0.2 L of the column is a rigid offset and the load sits on its tip.  For the flexible
    length l = L - a:  V(l) = 0 and M(l) = P a v'(l), i.e. v'' - k^2 a v' = 0 (the rigid arm tilts with the
    member end and the load moves sideways by a v')."""
    a = 0.2 * L
    ell = L - a

    def rows(k):
        sh = column_shape(k, ell)
        return [HELD(k, 0.0), NO_ROTATION(k, 0.0), NO_SHEAR(k, ell), sh[2] - k * k * a * sh[1]]
    ref = critical_k(rows, 6.0 / ell) ** 2 * EI
    b = buckle(model({"a": "fixed"}, members={"m": {"i": "a", "j": "b", "offsets": {"j": [0, 0, -a]}}}))
    assert b.factors[0] == pytest.approx(ref, rel=ACC)
    short = characteristic("fixed-free", length=ell) ** 2 * EI          # load at the end of the flexible length
    full = characteristic("fixed-free") ** 2 * EI                        # no rigid part
    assert full < ref < short


def test_negative_factors_are_reported_separately():
    """A tie in tension cannot buckle; reversing the load would: negative factor."""
    euler = characteristic("fixed-free") ** 2 * EI
    b = buckle(model({"a": "fixed"}, loads=[{"type": "nodal", "node": "b", "Fz": +1.0}], modes=2))
    assert len(b.factors) == 0 and b.modes.shape == (0, 2, 6)
    assert b.negative_factors[0] == pytest.approx(-euler, rel=ACC)
    assert abs(b.negative_factors[1]) > abs(b.negative_factors[0])
    assert any("no positive buckling factor" in w for w in b.info["warnings"])
    # one member in compression and one in tension: both signs
    doc = model({"a": "fixed", "c": "fixed"},
                loads=[{"type": "nodal", "node": "b", "Fz": -1.0}, {"type": "nodal", "node": "d", "Fz": 2.0}],
                nodes={"a": (0, 0, 0), "b": (0, 0, L), "c": (3000, 0, 0), "d": (3000, 0, L)},
                members={"m": ("a", "b"), "t": ("c", "d")})
    b = buckle(doc)
    assert b.factors[0] == pytest.approx(euler, rel=ACC) and b.negative_factors[0] == pytest.approx(-euler / 2, rel=ACC)


def test_modes_are_normalised_with_a_deterministic_sign():
    doc = model({"a": "fixed"}, modes=3)
    b1, b2 = buckle(doc), buckle(doc)
    assert np.array_equal(b1.modes, b2.modes) and np.array_equal(b1.factors, b2.factors)      # bitwise repeatable
    for mode in b1.modes:
        t = mode[:, :3].ravel()
        k = int(np.flatnonzero(np.abs(t) >= (1 - 1e-9) * np.abs(t).max())[0])
        assert t[k] > 0.0
    assert np.abs(b1.modes[0][:, :3]).max() == pytest.approx(1.0)        # cantilever: the tip is the largest point
    # mirror the structure: the mode keeps the positive sign on its largest component
    assert b1.modes[0][1, 0] == pytest.approx(1.0)
    # a mode with no node translation (pin-ended column): normalised on the member axis instead
    b = buckle(model({"a": PINNED, "b": ROLLER}))
    assert np.abs(b.modes[0][:, :3]).max() < 1e-9 and np.abs(b.modes[0][:, 3:]).max() > 0.0
    shape = b.member_shapes["m"][0]
    assert shape[:, 1].max() == pytest.approx(1.0, rel=1e-3) and shape[5, 1] == pytest.approx(1.0, rel=1e-3)
    assert shape[0, 1] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("pieces", [2, 3])
def test_member_split_by_the_user_gives_the_same_factors(pieces):
    whole = buckle(model({"a": "fixed"}, modes=2))
    nodes = {f"n{k}": (0, 0, L * k / pieces) for k in range(pieces + 1)}
    members = {f"m{k}": (f"n{k}", f"n{k + 1}") for k in range(pieces)}
    split = buckle(model({"n0": "fixed"}, loads=[{"type": "nodal", "node": f"n{pieces}", "Fz": -1.0}], modes=2,
                         nodes=nodes, members=members))
    assert split.factors == pytest.approx(whole.factors, rel=ACC)


def test_axial_force_varying_along_the_member_self_weight_buckling():
    """Column under its own weight only (N varies linearly): reference = eigenvalue of the equations with
    N(x) = -lambda w (L - x), found with solve_bvp for the factor lambda."""
    from scipy.integrate import solve_bvp
    gamma = 7.7e-5
    w = gamma * SEC["A"]

    def rhs(s, y, p):                                       # y = (v / L, theta, M L / EI, V L^2 / EI); no shear
        n = -p[0] * (1.0 - s)                               # N L^2 / EI
        return np.vstack([y[1], y[2], -y[3] + n * y[1], np.zeros_like(s)])

    def bc(ya, yb, p):
        return np.array([ya[0], ya[1], yb[2], yb[3], yb[0] - 1.0])
    s = np.linspace(0, 1, 81)
    v = 1 - np.cos(np.pi * s / 2)
    sol = solve_bvp(rhs, bc, s, np.vstack([v, np.gradient(v, s), np.gradient(np.gradient(v, s), s), 0 * s]), p=[8.0],
                    tol=1e-10, max_nodes=100000)
    assert sol.success
    ref = sol.p[0] * EI / L ** 2 / (w * L)                  # factor on the self weight (critical w L^3 / EI = 7.84)
    assert sol.p[0] == pytest.approx(7.84, rel=2e-3)
    b = buckle(model({"a": "fixed"}, loads=[{"type": "selfWeight"}]))
    assert b.factors[0] == pytest.approx(ref, rel=ACC)


def test_buckling_reference_state_is_first_order_and_ignores_alpha():
    base = buckle(model({"a": "fixed"}))
    doc = model({"a": "fixed"}, analysis={"alpha": 1.6})
    doc["combinations"]["c"]["analysisType"] = "secondOrder"
    doc["loadCases"]["c"]["loads"][0]["Fz"] = -1000.0
    b = buckle(doc)
    assert b.factors[0] * 1000.0 == pytest.approx(base.factors[0], rel=1e-9)
    assert "first-order reference state" in b.info["basis"] and b.info["subdivisions"]["m"] >= 4


def test_document_carries_factors_modes_and_accuracy():
    m = build_model(model({"a": "fixed"}, modes=2))
    out = analyze(m).to_document(m)["buckling"]["c"]
    assert out["factors"][0] == pytest.approx(characteristic("fixed-free") ** 2 * EI, rel=ACC)
    assert set(out["modes"][0]["b"]) == {"ux", "uy", "uz", "rx", "ry", "rz"} and out["modes"][0]["b"]["ux"] == pytest.approx(1.0)
    assert out["info"]["tolerance"] == BUCKLING_TOL and len(out["info"]["estimatedRelativeError"]) == 2
    assert "memberShapes" in out and "negativeFactors" not in out
    import json
    json.dumps(out)
