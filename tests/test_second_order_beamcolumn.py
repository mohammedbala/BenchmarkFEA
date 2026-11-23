"""Second-order analysis of single members against INDEPENDENT integrations of the beam-column
equations (``second_order_reference``: scipy.integrate.solve_bvp; nothing shared with the engine).

Members: a column a-b along +Z (local x = Z, local y = X, local z = Y) or a beam a-b along +X (local
y = Z, local z = -Y).  The section is stiffer about local y than about local z, so that the in-plane
(local xy) behaviour under study is not pre-empted by out-of-plane buckling.
"""
import math

import numpy as np
import pytest

from benchfea.model import build_model
from benchfea.solve import analyze

from second_order_reference import beam_column_bvp, euler_bernoulli_bvp
from solve_reference import E, G, make

SEC = {"A": 5000.0, "Iy": 4.0e8, "Iz": 4.0e7, "J": 3.0e5, "Asy": 2500.0, "Asz": 3000.0}
EI, GA, EA = E * SEC["Iz"], G * SEC["Asy"], E * SEC["A"]
L = 4000.0
PE_CANT = math.pi ** 2 * EI / (2.0 * L) ** 2          # only used to choose load levels (ratio P / Pe)
PE_PIN = math.pi ** 2 * EI / L ** 2
TOL = 2.0e-8                                            # reference accuracy (solve_bvp tolerance 1e-10)


def column(loads, shear=False, supports=None, members=None, nodes=None, analysis=None, **kw):
    doc = make(nodes or {"a": (0, 0, 0), "b": (0, 0, L)}, members or {"m": ("a", "b")},
               supports or {"a": "fixed"}, {"c": loads}, sections={"S": {"type": "properties", **SEC}},
               analysis={"shearDeformation": shear, **(analysis or {})}, **kw)
    doc["combinations"]["c"]["analysisType"] = "secondOrder"
    return doc


def beam(loads, shear=False, **kw):
    return column(loads, shear, nodes={"a": (0, 0, 0), "b": (L, 0, 0)},
                  supports={"a": {"type": "pinned", "dofs": {"rx": "fixed"}}, "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}},
                  **kw)


def solve(doc, cid="c"):
    return analyze(build_model(doc)).combinations[cid]


def compare(member, ref, scale_v, scale_m, tol=TOL, right=()):
    """Deflection and moment at every station against the reference."""
    worst = 0.0
    for k, x in enumerate(member.x):
        r = ref(x, right=(k in right))
        worst = max(worst, abs(member.displacements[k, 1] - r[0]) / scale_v, abs(member.forces[k, 5] - r[2]) / scale_m)
    assert worst <= tol, worst
    return worst


def test_the_two_references_agree_with_each_other():
    """First-order-system form versus the classical fourth-order form EI v'''' - N v'' = q."""
    P, q = 0.6 * PE_PIN, 3.0
    a = beam_column_bvp(L, EI, None, -P, -q, lambda s, e: [s[0, 0], s[0, 2], e[0, 0], e[0, 2]])
    b = euler_bernoulli_bvp(L, EI, -P, -q, lambda s, e: [s[0], s[2], e[0], e[2]])
    for x in np.linspace(0.0, L, 9):
        assert a(x)[0] == pytest.approx(b(x)[0], rel=1e-8, abs=1e-9 * abs(b(L / 2)[0]))
        assert a(x)[2] == pytest.approx(EI * b(x)[2], rel=1e-8, abs=1e-9 * abs(EI * b(L / 2)[2]))


@pytest.mark.parametrize("shear", [False, True])
@pytest.mark.parametrize("ratio", [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
def test_cantilever_beam_column_across_the_load_range(ratio, shear):
    """Fixed-base cantilever, axial load P = ratio x Pe and a lateral tip load: P-Delta and P-delta."""
    P, H = ratio * PE_CANT, 1000.0
    c = solve(column([{"type": "nodal", "node": "b", "Fz": -P, "Fx": H}], shear))
    ref = beam_column_bvp(L, EI, GA if shear else None, -P, 0.0,
                          lambda s, e: [s[0, 0], s[0, 1], e[0, 2], e[0, 3] - H])
    m = c.members["m"]
    compare(m, ref, abs(ref(L)[0]), abs(ref(0.0)[2]))
    assert c.displacements[1, 0] == pytest.approx(ref(L)[0], rel=TOL)
    assert c.reactions[0, 4] == pytest.approx(-ref(0.0)[2], rel=TOL)         # base moment about +Y: -(H L + P Delta)
    first = H * L ** 3 / (3 * EI) + (H * L / GA if shear else 0.0)
    assert c.displacements[1, 0] / first > 1.0 / (1.0 - ratio) * 0.98                  # the amplification is real
    assert c.converged and c.invariants.passed and c.analysis_type == "secondOrder"


@pytest.mark.parametrize("shear", [False, True])
@pytest.mark.parametrize("ratio", [0.2, 0.5, 0.8])
def test_simply_supported_beam_column_station_moments(ratio, shear):
    """No joint translation at all: only the member curvature effect (P-delta).  Every station moment
    and deflection is second-order accurate with a single element."""
    P, q = ratio * PE_PIN, 3.0
    c = solve(beam([{"type": "nodal", "node": "b", "Fx": -P},
                    {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q}], shear,
                   analysis={"resultStations": 21}))
    ref = beam_column_bvp(L, EI, GA if shear else None, -P, -q, lambda s, e: [s[0, 0], s[0, 2], e[0, 0], e[0, 2]])
    m = c.members["m"]
    compare(m, ref, abs(ref(L / 2)[0]), abs(ref(L / 2)[2]))
    assert c.solver["secondOrder"]["subdivisions"]["m"] <= 3          # a-priori rule k h <= 1.5 only
    mid = m.forces[10, 5]
    assert mid / (q * L ** 2 / 8.0) > 1.0 / (1.0 - ratio)                     # more than the first-order moment
    # end shears are components along the UNDEFORMED axis: they still carry half the load each
    assert m.end_forces[1] == pytest.approx(q * L / 2, rel=1e-9) and m.end_forces[7] == pytest.approx(q * L / 2, rel=1e-9)


@pytest.mark.parametrize("shear", [False, True])
def test_tension_stiffened_beam(shear):
    """Axial TENSION reduces deflections and moments (hyperbolic solution branch)."""
    T, q = 2.5 * PE_PIN, 3.0
    c = solve(beam([{"type": "nodal", "node": "b", "Fx": T},
                    {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q}], shear))
    ref = beam_column_bvp(L, EI, GA if shear else None, T, -q, lambda s, e: [s[0, 0], s[0, 2], e[0, 0], e[0, 2]])
    m = c.members["m"]
    compare(m, ref, abs(ref(L / 2)[0]), abs(ref(L / 2)[2]))
    assert m.forces[5, 5] < 0.4 * q * L ** 2 / 8.0 and m.forces[5, 0] == pytest.approx(T, rel=1e-12)
    assert c.invariants.passed


def test_very_large_tension_is_handled_by_subdividing():
    """k L = 40: the member is cut into pieces (a-priori rule k h <= 1.5) and the solution stays accurate."""
    T, q = 1600.0 * EI / L ** 2, 3.0
    c = solve(beam([{"type": "nodal", "node": "b", "Fx": T},
                    {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q}]))
    ref = beam_column_bvp(L, EI, None, T, -q, lambda s, e: [s[0, 0], s[0, 2], e[0, 0], e[0, 2]])
    m = c.members["m"]
    compare(m, ref, abs(ref(L / 2)[0]), q * EI / T, tol=1e-7)
    assert c.solver["secondOrder"]["subdivisions"]["m"] >= 27 and c.solver["secondOrder"]["largestPieceKh"] <= 3.0
    assert c.invariants.passed


@pytest.mark.parametrize("shear", [False, True])
def test_member_loads_inside_a_propped_beam_column(shear):
    """Fixed at a, roller at b, compression, a point load, a point moment and a partial trapezoidal load:
    fixed-end forces and station values with the second-order effect of the loads inside the member."""
    P = 0.9 * PE_PIN                                   # 0.44 of the fixed-pinned critical load
    xp, F, xm, Mo, x1, x2, w1, w2 = 0.3 * L, -9000.0, 0.55 * L, 4.0e6, 0.4 * L, 0.85 * L, -2.0, -6.0
    # beam along +X: local y = +Z, local z = -Y; all loads act in the local xy plane
    doc = column([{"type": "nodal", "node": "b", "Fx": -P},
                  {"type": "memberPoint", "member": "m", "direction": "Z", "value": F, "position": 0.3},
                  {"type": "memberPoint", "member": "m", "direction": "z", "kind": "moment", "value": Mo, "position": 0.55},
                  {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": w1, "w2": w2, "start": 0.4, "end": 0.85}],
                 shear, nodes={"a": (0, 0, 0), "b": (L, 0, 0)},
                 supports={"a": "fixed", "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}})
    c = solve(doc)

    def trapezoid(x):
        return w1 + (w2 - w1) * (x - x1) / (x2 - x1)
    q = [0.0, 0.0, trapezoid, trapezoid, 0.0]            # segments: 0-xp, xp-x1, x1-xm, xm-x2, x2-L

    cuts = [xp, x1, xm, x2]
    lengths = np.diff([0.0] + cuts + [L])
    jumps = {0: (0.0, F), 2: (Mo, 0.0)}                 # after segment k: (moment drop, shear drop)

    def bc(s, e):
        r = [s[0, 0], s[0, 1], e[-1, 0], e[-1, 2]]
        for k in range(len(lengths) - 1):
            dm, dv = jumps.get(k, (0.0, 0.0))
            r += [s[k + 1, 0] - e[k, 0], s[k + 1, 1] - e[k, 1], s[k + 1, 2] - (e[k, 2] - dm), s[k + 1, 3] - (e[k, 3] - dv)]
        return r
    ref = beam_column_bvp(lengths, EI, GA if shear else None, -P, q, bc)
    m = c.members["m"]
    scale_v = max(abs(ref(x)[0]) for x in np.linspace(0, L, 41))
    scale_m = max(abs(ref(x)[2]) for x in np.linspace(0, L, 41))
    right = {k for k in range(1, len(m.x)) if m.x[k] == m.x[k - 1]}          # second station at a point load
    compare(m, ref, scale_v, scale_m, right=right)
    for k in right:                                                          # the jumps are reported
        if m.x[k] == xp:
            assert m.forces[k, 1] - m.forces[k - 1, 1] == pytest.approx(-F, rel=1e-9)
        if m.x[k] == xm:
            assert m.forces[k, 5] - m.forces[k - 1, 5] == pytest.approx(-Mo, rel=1e-9)
    assert c.reactions[0, 4] == pytest.approx(ref(0.0)[2], rel=1e-7)          # fixing moment (+Y = -local z)
    assert c.invariants.passed


@pytest.mark.parametrize("shear", [False, True])
def test_axial_force_varying_along_the_member(shear):
    """A heavy column: its own weight makes the axial force vary linearly.  The reference integrates the
    equations with N(x); the engine subdivides until its results change by less than the tolerance."""
    gamma = 7.7e-5 * 400.0
    w = gamma * SEC["A"]                                 # N per mm; total weight = 0.5 Pe
    P, H, qx = 0.4 * PE_CANT, 1000.0, 0.5
    doc = column([{"type": "nodal", "node": "b", "Fz": -P, "Fx": H}, {"type": "selfWeight"},
                  {"type": "memberDistributed", "member": "m", "direction": "X", "w1": qx}], shear,
                 materials={"M": {"E": E, "G": G, "unitWeight": gamma}})
    c = solve(doc)
    ref = beam_column_bvp(L, EI, GA if shear else None, lambda x: -(P + w * (L - x)), qx,
                          lambda s, e: [s[0, 0], s[0, 1], e[0, 2], e[0, 3] - H])
    m = c.members["m"]
    worst = compare(m, ref, abs(ref(L)[0]), abs(ref(0.0)[2]), tol=1e-6)
    assert worst < 1e-8                                   # actual accuracy, far inside the stated tolerance
    so = c.solver["secondOrder"]
    assert so["subdivisionCheck"]["change"] <= so["tolerance"] and so["subdivisions"]["m"] >= 4
    assert m.forces[0, 0] == pytest.approx(-(P + w * L), rel=1e-12) and m.forces[-1, 0] == pytest.approx(-P, rel=1e-12)
    first = analyze(build_model({**doc, "combinations": {"c": {"factors": {"c": 1.0}}}})).combinations["c"]
    assert c.displacements[1, 0] > 2.0 * first.displacements[1, 0]            # the comparison is not trivial
    assert c.invariants.passed


def test_axial_load_in_the_span_changes_the_axial_force_piecewise():
    """Stepped axial force: a point load along the member axis at mid-height."""
    P1, P2, H = 0.3 * PE_CANT, 0.5 * PE_CANT, 1000.0
    c = solve(column([{"type": "nodal", "node": "b", "Fz": -P1, "Fx": H},
                      {"type": "memberPoint", "member": "m", "direction": "Z", "value": -P2, "position": 0.5}]))

    def bc(s, e):
        return [s[0, 0], s[0, 1], e[1, 2], e[1, 3] - H,
                s[1, 0] - e[0, 0], s[1, 1] - e[0, 1], s[1, 2] - e[0, 2], s[1, 3] - e[0, 3]]
    ref = beam_column_bvp([L / 2, L / 2], EI, None, [-(P1 + P2), -P1], 0.0, bc)
    m = c.members["m"]
    right = {k for k in range(1, len(m.x)) if m.x[k] == m.x[k - 1]}
    compare(m, ref, abs(ref(L)[0]), abs(ref(0.0)[2]), right=right)
    assert c.solver["secondOrder"]["subdivisions"]["m"] == 2


@pytest.mark.parametrize("pieces", [2, 3, 5])
def test_a_member_split_by_the_user_gives_the_same_answer(pieces):
    """Constant axial force: one element is already exact, so splitting changes nothing (round-off)."""
    P, H, q = 0.7 * PE_CANT, 1000.0, 0.8
    whole = solve(column([{"type": "nodal", "node": "b", "Fz": -P, "Fx": H},
                          {"type": "memberDistributed", "member": "m", "direction": "X", "w1": q}], True))
    nodes = {f"n{k}": (0, 0, L * k / pieces) for k in range(pieces + 1)}
    members = {f"m{k}": (f"n{k}", f"n{k + 1}") for k in range(pieces)}
    loads = [{"type": "nodal", "node": f"n{pieces}", "Fz": -P, "Fx": H}]
    loads += [{"type": "memberDistributed", "member": mid, "direction": "X", "w1": q} for mid in members]
    split = solve(column(loads, True, nodes=nodes, members=members, supports={"n0": "fixed"}))
    assert split.displacements[-1] == pytest.approx(whole.displacements[1], rel=1e-9, abs=1e-12)
    assert split.reactions[0] == pytest.approx(whole.reactions[0], rel=1e-9, abs=1e-6)
    assert split.members["m0"].forces[0] == pytest.approx(whole.members["m"].forces[0], rel=1e-9, abs=1e-6)


def test_split_equivalence_with_a_varying_axial_force_holds_within_the_tolerance():
    gamma = 7.7e-5 * 400.0
    P, H = 0.4 * PE_CANT, 1000.0
    mats = {"M": {"E": E, "G": G, "unitWeight": gamma}}
    whole = solve(column([{"type": "nodal", "node": "b", "Fz": -P, "Fx": H}, {"type": "selfWeight"}], materials=mats))
    nodes = {f"n{k}": (0, 0, L * k / 3) for k in range(4)}
    members = {f"m{k}": (f"n{k}", f"n{k + 1}") for k in range(3)}
    split = solve(column([{"type": "nodal", "node": "n3", "Fz": -P, "Fx": H}, {"type": "selfWeight"}],
                         nodes=nodes, members=members, supports={"n0": "fixed"}, materials=mats))
    tol = 1e-6                                           # analysis.convergenceTolerance (default)
    assert split.displacements[-1, 0] == pytest.approx(whole.displacements[1, 0], rel=tol)
    assert split.reactions[0, 4] == pytest.approx(whole.reactions[0, 4], rel=tol)


def test_second_order_in_the_other_bending_plane_and_for_a_rolled_member():
    """The same cantilever rotated by 90 degrees about its axis bends in its local xz plane: the result
    must be the mirror image (this exercises the sign mapping of the second plane)."""
    P, H = 0.6 * PE_CANT, 1000.0
    ref = solve(column([{"type": "nodal", "node": "b", "Fz": -P, "Fx": H}], True))
    sec = {**SEC, "Iy": SEC["Iz"], "Iz": SEC["Iy"], "Asy": SEC["Asz"], "Asz": SEC["Asy"]}
    doc = column([{"type": "nodal", "node": "b", "Fz": -P, "Fx": H}], True,
                 members={"m": {"i": "a", "j": "b", "orientation": {"roll": 90}}})
    doc["sections"]["S"] = {"type": "properties", **sec}
    rolled = solve(doc)
    assert rolled.displacements[1] == pytest.approx(ref.displacements[1], rel=1e-10, abs=1e-14)
    # local z of the rolled member = -X: the tip load is -z, the base moment is about local y
    assert rolled.members["m"].forces[0, 4] == pytest.approx(ref.members["m"].forces[0, 5], rel=1e-10)
    assert rolled.members["m"].displacements[-1, 2] == pytest.approx(-ref.members["m"].displacements[-1, 1], rel=1e-10)


def test_axial_shortening_is_first_order_and_reported():
    P, H = 0.5 * PE_CANT, 1000.0
    c = solve(column([{"type": "nodal", "node": "b", "Fz": -P, "Fx": H}]))
    assert c.displacements[1, 2] == pytest.approx(-P * L / EA, rel=1e-12)     # no bowing term (stated limit)
