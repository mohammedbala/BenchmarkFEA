"""The OpenSees and PyNite translators against (i) the ODE oracle, (ii) hand statics / kinematics and
(iii) each other.  benchfea.solve is not used here.  Tests skip when a foreign package is absent.
"""
import math

import numpy as np
import pytest

from benchfea.crosscheck import compare, cross_check, equilibrium
from benchfea.crosscheck.fidelity import APPROXIMATE, EXACT, UNSUPPORTED
from benchfea.model import build_model
from benchfea.oracle import Bar, Beam

from crosscheck_helpers import (EA, EIy, EIz, GAy, GAz, GJ, SEC, SOLVER_PARAMS, document, local_node,
                                needs_both, needs_opensees, rel, solve)

# direction of the member axis, orientation -- aligned, skew, skew + roll, vertical with reference vector
ORIENTATIONS = {
    "along-X": ((1.0, 0.0, 0.0), None),
    "skew": ((2.0, 1.0, 0.5), None),
    "skew-rolled": ((2.0, -1.0, 1.5), {"roll": 33.0}),
    "vertical-ref": ((0.0, 0.0, 1.0), {"referenceVector": [1.0, 2.0, 0.0], "roll": -20.0}),
}


def _member(i, j, orientation, **extra):
    m = {"i": i, "j": j, **extra}
    if orientation:
        m["orientation"] = orientation
    return m


# ====================================================================================== local axes
AXES_CASES = [((1, 0, 0), None, "Z"), ((-1, 0, 0), None, "Z"), ((0, 1, 0), None, "Z"), ((0, 0, 1), None, "Z"),
              ((0, 0, -1), None, "Z"), ((1, 2, 3), None, "Z"), ((1, 2, 3), {"roll": 40.0}, "Z"),
              ((3, -1, 0.5), {"referenceVector": [0, 1, 0], "roll": -75.0}, "Z"),
              ((0, 0, 1), {"roll": 90.0}, "Z"), ((1e-3, 0, 1), {"referenceVector": [0, 1, 1]}, "Z"),
              ((1, 0, 0), None, "Y"), ((0, 1, 0), None, "Y"), ((0, -1, 0), {"roll": 30.0}, "Y"),
              ((0, 0, 1), None, "Y"), ((2, 3, -1), {"roll": 200.0}, "Y"), ((1, 1e-4, 0), None, "Y")]


@pytest.mark.parametrize("solver", SOLVER_PARAMS)
@pytest.mark.parametrize("vertical", ["Z", "Y"])
def test_member_axes_reproduced_exactly(solver, vertical):
    """Independent cantilevers (Iy != Iz), each loaded at the tip by a force along its LOCAL y and a
    torque about its local x.  Any error in the triad given to the foreign solver changes the result:
    in benchfea's axes the tip must move by P L^3/3EIz along y only, rotate by P L^2/2EIz about z and
    T L/GJ about x; the end forces must be Fy_i = -P, Mz_i = -P L, Mx_i = -T (statics)."""
    Lm, P, T = 2000.0, 1500.0, 2.0e5
    nodes, members, supports, loads = {}, {}, {}, []
    k = 0
    for d, ori, vert in AXES_CASES:
        if vert != vertical:
            continue
        d = np.array(d, float) / np.linalg.norm(d)
        base = np.array([5000.0 * k, 0.0, 0.0])
        nodes[f"a{k}"], nodes[f"b{k}"] = base, base + Lm * d
        members[f"m{k}"] = _member(f"a{k}", f"b{k}", ori)
        supports[f"a{k}"] = "fixed"
        loads += [{"type": "memberPoint", "member": f"m{k}", "direction": "y", "value": P, "position": 1.0},
                  {"type": "memberPoint", "member": f"m{k}", "direction": "x", "kind": "moment", "value": T,
                   "position": 1.0}]
        k += 1
    model = build_model(document(nodes, members, supports, {"L": loads}, vertical=vertical))
    r = solve(solver, model, "L")
    assert r.ok, r.error
    assert r.fidelity.level == EXACT
    for j in range(k):
        u, th = local_node(model, r, f"m{j}", f"b{j}")
        v = P * Lm ** 3 / (3 * EIz)
        assert u == pytest.approx([0.0, v, 0.0], abs=1e-9 * v)
        assert th == pytest.approx([T * Lm / GJ, 0.0, P * Lm ** 2 / (2 * EIz)], abs=1e-9 * T * Lm / GJ)
        f = r.end_forces[f"m{j}"]
        assert f[:6] == pytest.approx([0.0, -P, 0.0, -T, 0.0, -P * Lm], abs=1e-8 * P * Lm)
    if solver == "pynite":
        assert r.info["axesError"] < 1e-9


# ====================================================================================== against the oracle
def _two_member_model(direction, orientation, shear=False, releases=True):
    """Clamped - B - clamped line of two members of length L along ``direction`` with a rotational
    spring about local z and a hinge about local y at B, loaded in both bending planes, axially and in
    torsion."""
    Lm = 2500.0
    d = np.array(direction, float) / np.linalg.norm(direction)
    nodes = {"A": np.zeros(3), "B": Lm * d, "C": 2 * Lm * d}
    rel1 = {"j": {"ry": True}} if releases else {}
    rel2 = {"i": {"rz": {"spring": 6.0e9}}} if releases else {}
    members = {"m1": _member("A", "B", orientation, releases=rel1), "m2": _member("B", "C", orientation, releases=rel2)}
    loads = [
        {"type": "memberDistributed", "member": "m1", "direction": "y", "w1": -2.0, "w2": -7.0, "start": 0.2, "end": 0.9},
        {"type": "memberPoint", "member": "m2", "direction": "y", "value": -6e3, "position": 0.6},
        {"type": "memberPoint", "member": "m1", "direction": "z", "kind": "moment", "value": 3e6, "position": 0.5},
        {"type": "memberDistributed", "member": "m2", "direction": "z", "w1": 3.0},
        {"type": "memberPoint", "member": "m1", "direction": "z", "value": 4e3, "position": 0.3},
        {"type": "memberPoint", "member": "m2", "direction": "y", "kind": "moment", "value": -2e6, "position": 0.4},
        {"type": "memberDistributed", "member": "m1", "direction": "x", "w1": 1.0, "w2": 4.0},
        {"type": "memberPoint", "member": "m2", "direction": "x", "value": 9e3, "position": 0.5},
        {"type": "memberPoint", "member": "m1", "direction": "x", "kind": "moment", "value": 5e5, "position": 0.7},
    ]
    doc = document(nodes, members, {"A": "fixed", "C": "fixed"}, {"L": loads}, analysis={"shearDeformation": shear})
    return build_model(doc), Lm


def _two_member_oracle(Lm, shear, releases=True):
    xy = Beam(2 * Lm, EIz, GAy if shear else None).support(0.0, v="fixed", r="fixed").support(2 * Lm, v="fixed", r="fixed")
    xy.dist(0.2 * Lm, 0.9 * Lm, -2.0, -7.0).point(Lm + 0.6 * Lm, P=-6e3).point(0.5 * Lm, M=3e6)
    xz = Beam(2 * Lm, EIy, GAz if shear else None, plane="xz").support(0.0, v="fixed", r="fixed") \
        .support(2 * Lm, v="fixed", r="fixed")
    xz.dist(Lm, 2 * Lm, 3.0).point(0.3 * Lm, P=4e3).point(Lm + 0.4 * Lm, M=-2e6)
    if releases:
        xy.hinge(Lm, k=6.0e9)
        xz.hinge(Lm)
    ax = Bar(2 * Lm, EA).support(0.0).support(2 * Lm).dist(0.0, Lm, 1.0, 4.0).point(1.5 * Lm, 9e3)
    to = Bar(2 * Lm, GJ).support(0.0).support(2 * Lm).point(0.7 * Lm, 5e5)
    return xy.solve(tol=1e-11), xz.solve(tol=1e-11), ax.solve(tol=1e-11), to.solve(tol=1e-11)


@pytest.mark.parametrize("solver", SOLVER_PARAMS)
@pytest.mark.parametrize("name", list(ORIENTATIONS))
def test_two_members_with_spring_and_hinge_against_oracle(solver, name):
    direction, orientation = ORIENTATIONS[name]
    model, Lm = _two_member_model(direction, orientation)
    xy, xz, ax, to = _two_member_oracle(Lm, False)
    r = solve(solver, model, "L")
    assert r.ok, r.error
    assert r.fidelity.level == EXACT, r.fidelity.to_dict()
    u, th = local_node(model, r, "m1", "B")
    assert rel(u[1], xy.v(Lm)) < 1e-8 and rel(u[2], xz.v(Lm)) < 1e-8 and rel(u[0], ax.u(Lm)) < 1e-8
    assert rel(th[0], to.u(Lm)) < 1e-8
    # the joint carries the rotation of the side that is rigidly attached to it:
    # about z the spring is at the start of m2 -> joint follows m1 (left of the hinge);
    # about y the hinge is at the end of m1 -> joint follows m2 (right of the hinge)
    assert rel(th[2], xy.theta(Lm, -1)) < 1e-8 and rel(th[1], xz.theta(Lm, +1)) < 1e-8
    f1, f2 = r.end_forces["m1"], r.end_forces["m2"]
    F, M = 2.0e4, 1.5e7                                      # force / moment scales of the problem
    expect1_i = [-ax.N(0.0, +1), -xy.Q(0.0, +1), -xz.Q(0.0, +1), -to.N(0.0, +1), -xz.M(0.0, +1), -xy.M(0.0, +1)]
    expect1_j = [ax.N(Lm, -1), xy.Q(Lm, -1), xz.Q(Lm, -1), to.N(Lm, -1), xz.M(Lm, -1), xy.M(Lm, -1)]
    expect2_j = [ax.N(2 * Lm), xy.Q(2 * Lm), xz.Q(2 * Lm), to.N(2 * Lm), xz.M(2 * Lm), xy.M(2 * Lm)]
    scale = np.array([F] * 3 + [M] * 3)
    assert np.all(np.abs(f1[:6] - expect1_i) / scale < 1e-8)
    assert np.all(np.abs(f1[6:] - expect1_j) / scale < 1e-8)
    assert np.all(np.abs(f2[6:] - expect2_j) / scale < 1e-8)
    assert abs(f1[10]) < 1e-9 * M                            # hinge about y at the end of m1
    eq = equilibrium(model, r)
    assert eq["force"] < 1e-10 and eq["moment"] < 1e-10


@needs_opensees
@pytest.mark.parametrize("name", ["along-X", "skew-rolled"])
def test_opensees_shear_deformation_with_releases_against_oracle(name):
    """Timoshenko members with a rotational spring and a hinge: exact through parallel elements."""
    direction, orientation = ORIENTATIONS[name]
    model, Lm = _two_member_model(direction, orientation, shear=True)
    xy, xz, ax, to = _two_member_oracle(Lm, True)
    r = solve("opensees", model, "L")
    assert r.ok and r.fidelity.level == EXACT
    u, th = local_node(model, r, "m1", "B")
    assert rel(u[1], xy.v(Lm)) < 1e-8 and rel(u[2], xz.v(Lm)) < 1e-8
    f1 = r.end_forces["m1"]
    assert rel(f1[5], -xy.M(0.0, +1)) < 1e-8 and rel(f1[4], -xz.M(0.0, +1)) < 1e-8
    assert rel(f1[1], -xy.Q(0.0, +1)) < 1e-8 and rel(f1[2], -xz.Q(0.0, +1)) < 1e-8
    # the shear deformation matters here (otherwise the test would prove nothing)
    eb = _two_member_oracle(Lm, False)[0]
    assert rel(eb.v(Lm), xy.v(Lm)) > 1e-3


@pytest.mark.skipif(not __import__("crosscheck_helpers").HAVE["pynite"], reason="PyNiteFEA is not installed")
def test_pynite_refuses_shear_deformation():
    model, _ = _two_member_model((1.0, 0.0, 0.0), None, shear=True)
    r = solve("pynite", model, "L")
    assert r.ok and r.fidelity.level == UNSUPPORTED
    assert "shearDeformation" in [i.feature for i in r.fidelity.unsupported]


@needs_opensees
@pytest.mark.parametrize("native", [True, False])
def test_opensees_native_and_equivalent_loads_agree(native):
    """Both load routes (eleLoad where OpenSees is exact / equivalent nodal loads everywhere) must give
    the oracle's answer; the unreleased model exercises eleLoad."""
    model, Lm = _two_member_model((2.0, 1.0, 0.5), {"roll": 15.0}, releases=False)
    xy, xz, ax, to = _two_member_oracle(Lm, False, releases=False)
    r = solve("opensees", model, "L", native_loads=native)
    assert r.ok and r.fidelity.level == EXACT
    notes = " ".join(i.note for i in r.fidelity.items.values())
    assert ("eleLoad" in notes) == native
    u, th = local_node(model, r, "m1", "B")
    assert rel(u[1], xy.v(Lm)) < 1e-8 and rel(u[2], xz.v(Lm)) < 1e-8 and rel(th[2], xy.theta(Lm)) < 1e-8
    assert rel(r.end_forces["m1"][5], -xy.M(0.0, +1)) < 1e-8
    assert rel(r.end_forces["m2"][10], xz.M(2 * Lm)) < 1e-8


# ====================================================================================== offsets
@pytest.mark.parametrize("solver", SOLVER_PARAMS)
def test_rigid_offsets_hand_solution(solver):
    """Cantilever along X from A (fixed).  Rigid offsets: e_i along the axis at A, and at the far end
    the node B sits e_z above the flexible end (offset_j = (0, 0, -e_z) from node to member end).
    A force Fx at B then bends the member: at the flexible end it is Fx plus the moment of the
    lever arm,  M_end = (r x F) with r = -offset_j = (0, 0, e_z):  My = e_z Fx.
    Flexible length Lf: member local axes y = +Z, z = -Y, so My(global) = -Mz(local):
    Mz_local = -e_z Fx constant along the member;  th_z,local(end) = Mz Lf / EIz,
    v(end) = Mz Lf^2 / 2EIz, u(end) = Fx Lf / EA.
    Node B:  u_B = u_end - theta x offset_j  (CONVENTIONS section 2)."""
    Ltot, ei, ez, Fx = 3000.0, 250.0, 180.0, 4.0e4
    model = build_model(document({"A": (0, 0, 0), "B": (Ltot, 0, ez)},
                                 {"m": {"i": "A", "j": "B", "offsets": {"i": [ei, 0, 0], "j": [0, 0, -ez]}}},
                                 {"A": "fixed"}, {"L": [{"type": "nodal", "node": "B", "Fx": Fx}]}))
    Lf = Ltot - ei
    assert model.members["m"].length == pytest.approx(Lf)
    r = solve(solver, model, "L")
    assert r.ok, r.error
    Mz = -ez * Fx
    th_local = np.array([0.0, 0.0, Mz * Lf / EIz])
    u_local = np.array([Fx * Lf / EA, Mz * Lf ** 2 / (2 * EIz), 0.0])
    R = model.members["m"].R
    th = R.T @ th_local
    uB = R.T @ u_local - np.cross(th, np.array([0.0, 0.0, -ez]))
    tol = 1e-9 if solver == "opensees" else 5e-4            # PyNite: stiff-stub emulation
    assert np.abs(r.displacements[1, :3] - uB).max() < tol * np.abs(uB).max()
    assert np.abs(r.displacements[1, 3:] - th).max() < tol * np.abs(th).max()
    f = r.end_forces["m"]                                    # at the ends of the FLEXIBLE length
    assert f[6:] == pytest.approx([Fx, 0.0, 0.0, 0.0, 0.0, Mz], abs=tol * abs(Mz))
    assert r.reactions[0] == pytest.approx([-Fx, 0, 0, 0, -ez * Fx, 0], abs=tol * abs(Mz))
    item = r.fidelity.items["offsets"]
    if solver == "opensees":
        assert item.level == EXACT
    else:
        assert item.level == APPROXIMATE and item.error is not None
        assert np.abs(r.displacements[1, :3] - uB).max() / np.abs(uB).max() <= item.error


# ====================================================================================== supports
@pytest.mark.parametrize("solver", SOLVER_PARAMS)
def test_support_springs_and_imposed_displacements_against_oracle(solver):
    """Beam along X (local y = +Z, local z = -Y): A clamped with an imposed settlement uz and rotation
    ry, B on a translational spring (uz) and a rotational spring (ry), uniform load in -Z.
    Rotations about local z are -ry."""
    Lm, k, kr, d, phi, q = 4000.0, 800.0, 3.0e9, -4.0, 0.0015, 3.0
    model = build_model(document(
        {"A": (0, 0, 0), "B": (Lm, 0, 0)}, {"m": ("A", "B")},
        {"A": "fixed", "B": {"dofs": {"uz": {"spring": k}, "ry": {"spring": kr}, "uy": "fixed", "rx": "fixed"}}},
        {"L": [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q},
               {"type": "supportDisplacement", "node": "A", "uz": d, "ry": phi}]}))
    s = Beam(Lm, EIz).support(0.0, v="fixed", r="fixed", settlement=d, rotation=-phi).support(Lm, v=k, r=kr) \
        .dist(0.0, Lm, -q).solve(tol=1e-11)
    r = solve(solver, model, "L")
    assert r.ok and r.fidelity.level == EXACT
    assert r.displacements[0, 2] == pytest.approx(d, rel=1e-12) and r.displacements[0, 4] == pytest.approx(phi, rel=1e-12)
    assert rel(r.displacements[1, 2], s.v(Lm)) < 1e-8 and rel(r.displacements[1, 4], -s.theta(Lm)) < 1e-8
    RA, MA = s.reactions[0.0]
    RB, MB = s.reactions[Lm]
    assert rel(r.reactions[0, 2], RA) < 1e-8 and rel(r.reactions[0, 4], -MA) < 1e-8
    assert rel(r.reactions[1, 2], RB) < 1e-8 and rel(r.reactions[1, 4], -MB) < 1e-8
    assert rel(r.reactions[1, 2], -k * r.displacements[1, 2]) < 1e-9         # spring reaction = -k u


# ====================================================================================== trusses
@pytest.mark.parametrize("solver", SOLVER_PARAMS)
def test_truss_hand_solution(solver):
    """Two bars A-C and B-C, pinned supports, load (H, -P) at C(a, h).  Joint equilibrium with the unit
    vectors e1 = (a, h)/L, e2 = (-a, h)/L from the supports to C (tension N pulls C towards the
    support): -N1 e1 - N2 e2 + (H, -P) = 0  ->  N1 = (H/a - P/h) L/2,  N2 = (-H/a - P/h) L/2.
    Displacement by the unit-load theorem: u = sum N n L / EA with n from a unit load."""
    a, h, H, P = 2000.0, 1500.0, 3.0e4, 8.0e4
    Lb = math.hypot(a, h)
    model = build_model(document({"A": (0, 0, 0), "B": (2 * a, 0, 0), "C": (a, 0, h)},
                                 {"t1": {"i": "A", "j": "C", "type": "truss"}, "t2": {"i": "B", "j": "C", "type": "truss"}},
                                 {"A": "pinned", "B": "pinned", "C": {"dofs": {"uy": "fixed"}}},
                                 {"L": [{"type": "nodal", "node": "C", "Fx": H, "Fz": -P}]}))
    N1, N2 = (H / a - P / h) * Lb / 2, (-H / a - P / h) * Lb / 2
    r = solve(solver, model, "L")
    assert r.ok and r.fidelity.level == EXACT
    assert r.end_forces["t1"][6] == pytest.approx(N1, rel=1e-10)          # Fx at end j = tension
    assert r.end_forces["t2"][6] == pytest.approx(N2, rel=1e-10)
    assert np.abs(np.delete(r.end_forces["t1"], [0, 6])).max() < 1e-9 * abs(N1)
    n1x, n2x = Lb / (2 * a), -Lb / (2 * a)                                 # unit load along +x at C
    n1z, n2z = Lb / (2 * h), Lb / (2 * h)                                  # unit load along +z at C
    ux = (N1 * n1x + N2 * n2x) * Lb / EA
    uz = (N1 * n1z + N2 * n2z) * Lb / EA
    assert r.displacements[2, 0] == pytest.approx(ux, rel=1e-10) and r.displacements[2, 2] == pytest.approx(uz, rel=1e-10)
    assert np.all(r.displacements[2, 3:] == 0.0)             # rotations of a truss-only node are reported as 0
    assert r.reactions[0, :3] == pytest.approx(-N1 * np.array([a, 0, h]) / Lb, rel=1e-10)


@pytest.mark.parametrize("solver", SOLVER_PARAMS)
def test_truss_member_loads_go_to_the_joints_as_simple_beam_reactions(solver):
    """A transverse load on a truss-type member is carried by a simply supported span between its ends
    (ARCHITECTURE): end shears = simple-beam reactions, no end moments, no axial force."""
    Lb, q = 3000.0, 2.0
    model = build_model(document({"A": (0, 0, 0), "B": (Lb, 0, 0)}, {"t": {"i": "A", "j": "B", "type": "truss"}},
                                 {"A": "pinned", "B": "pinned"},
                                 {"L": [{"type": "memberDistributed", "member": "t", "direction": "Z", "w1": 0.0, "w2": -q}]}))
    r = solve(solver, model, "L")
    assert r.ok and r.fidelity.level == EXACT
    assert r.reactions[0, 2] == pytest.approx(q * Lb / 6, rel=1e-10)
    assert r.reactions[1, 2] == pytest.approx(q * Lb / 3, rel=1e-10)
    f = r.end_forces["t"]
    assert f[1] == pytest.approx(q * Lb / 6, rel=1e-10) and f[7] == pytest.approx(q * Lb / 3, rel=1e-10)
    assert abs(f[5]) < 1e-9 * q * Lb * Lb and abs(f[11]) < 1e-9 * q * Lb * Lb


# ====================================================================================== constraints
@pytest.mark.parametrize("solver", SOLVER_PARAMS)
def test_rigid_link_hand_solution(solver):
    """Column A-B (fixed base, height h) with a slave S rigidly linked to B at r = (e, d, 0), loaded by
    Fz = -P.  Statics: at B the force -P and the moment r x F = (-d P, e P, 0).  Column axes: x = +Z,
    y = +X, z = +Y, so bending about global Y uses EIz and about global X uses EIy:
    ry = e P h / EIz,  ux = e P h^2 / 2EIz,  rx = -d P h / EIy,  uy = -rx h/2 ... = d P h^2 / 2EIy,
    uz = -P h / EA.  Slave: u_S = u_B + theta x r, theta_S = theta_B."""
    h, e, d, P = 3000.0, 700.0, 400.0, 5.0e4
    model = build_model(document({"A": (0, 0, 0), "B": (0, 0, h), "S": (e, d, h)}, {"c": ("A", "B")}, {"A": "fixed"},
                                 {"L": [{"type": "nodal", "node": "S", "Fz": -P}]},
                                 constraints={"link": {"type": "rigidLink", "master": "B", "slaves": ["S"]}}))
    ry, rx = e * P * h / EIz, -d * P * h / EIy
    uB = np.array([ry * h / 2, -rx * h / 2, -P * h / EA])
    th = np.array([rx, ry, 0.0])
    uS = uB + np.cross(th, [e, d, 0.0])
    r = solve(solver, model, "L")
    assert r.ok, r.error
    tol = 1e-10 if solver == "opensees" else 1e-3
    assert np.abs(r.displacements[1, :3] - uB).max() < tol * np.abs(uS).max()
    assert np.abs(r.displacements[2, :3] - uS).max() < tol * np.abs(uS).max()
    assert np.abs(r.displacements[2, 3:] - th).max() < tol * np.abs(th).max()
    assert r.reactions[0] == pytest.approx([0, 0, P, d * P, -e * P, 0], abs=tol * e * P)
    item = r.fidelity.items["constraints.rigidLink"]
    assert item.level == (EXACT if solver == "opensees" else APPROXIMATE)
    if solver == "pynite":
        assert np.abs(r.displacements[2, :3] - uS).max() / np.abs(uS).max() <= max(item.error, 1e-12)


@needs_both
def test_rigid_diaphragm_kinematics_and_agreement():
    """One-storey frame with unequal columns and a diaphragm (normal Z).  OpenSees must satisfy the
    diaphragm relation exactly (in-plane translations and rz follow the master; uz, rx, ry of the
    slaves stay independent); PyNite's stiff-member emulation must agree within its own estimate."""
    h = 3000.0
    nodes = {"a": (0, 0, 0), "b": (6000, 0, 0), "c": (6000, 4000, 0), "d": (0, 4000, 0),
             "A": (0, 0, h), "B": (6000, 0, h), "C": (6000, 4000, h), "D": (0, 4000, h), "M": (2500, 1500, h)}
    sections = {"S": {"type": "properties", **SEC},
                "T": {"type": "properties", "A": 9000.0, "Iy": 3.0e7, "Iz": 9.0e7, "J": 8.0e5}}
    members = {"ca": ("a", "A"), "cb": {"i": "b", "j": "B", "section": "T"},
               "cc": {"i": "c", "j": "C", "orientation": {"roll": 90.0}}, "cd": ("d", "D"),
               "g1": ("A", "B"), "g2": {"i": "D", "j": "C", "section": "T"}}
    model = build_model(document(nodes, members, {k: "fixed" for k in "abcd"},
                                 {"L": [{"type": "nodal", "node": "M", "Fx": 9e4, "Fy": 4e4, "Mz": 6e7},
                                        {"type": "memberDistributed", "member": "g1", "direction": "Z", "w1": -8.0}]},
                                 sections=sections,
                                 constraints={"floor": {"type": "rigidDiaphragm", "master": "M",
                                                        "slaves": ["A", "B", "C", "D"]}}))
    ro = solve("opensees", model, "L")
    assert ro.ok and ro.fidelity.level == EXACT and ro.info["constraintViolation"] < 1e-12
    um = ro.displacements[model.nodes["M"].index]
    for s in "ABCD":
        us = ro.displacements[model.nodes[s].index]
        arm = model.nodes[s].xyz - model.nodes["M"].xyz
        expect = um[:2] + um[5] * np.array([-arm[1], arm[0]])           # theta_z e_z x arm
        assert us[:2] == pytest.approx(expect, abs=1e-11 * np.abs(um[:2]).max())
        assert us[5] == pytest.approx(um[5], rel=1e-11)
    assert abs(ro.displacements[model.nodes["A"].index, 4]) > 1e-6       # out-of-plane rotations are free
    assert np.all(um[[2, 3, 4]] == 0.0)                      # master DOFs without stiffness are reported as 0
    eq = equilibrium(model, ro)
    assert eq["force"] < 1e-11 and eq["moment"] < 1e-11
    rp = solve("pynite", model, "L")
    assert rp.ok and rp.fidelity.level == APPROXIMATE
    est = rp.fidelity.estimated_error
    c = compare(ro, rp, model)
    assert 0.0 < c.max_relative <= est < 5e-3


@pytest.mark.parametrize("solver", SOLVER_PARAMS)
def test_diaphragm_slave_out_of_plane_is_unsupported(solver):
    nodes = {"a": (0, 0, 0), "A": (0, 0, 3000), "b": (4000, 0, 0), "B": (4000, 0, 3300), "M": (2000, 500, 3000)}
    model = build_model(document(nodes, {"ca": ("a", "A"), "cb": ("b", "B")}, {"a": "fixed", "b": "fixed"},
                                 {"L": [{"type": "nodal", "node": "M", "Fx": 1e4}]},
                                 constraints={"f": {"type": "rigidDiaphragm", "master": "M", "slaves": ["A", "B"]}}))
    r = solve(solver, model, "L", estimate_errors=False)
    assert "constraints.rigidDiaphragm" in [i.feature for i in r.fidelity.unsupported]


# ====================================================================================== releases needing a node chain
@needs_opensees
def test_opensees_axial_spring_and_torsion_release_exact_when_axes_are_global():
    """Column along Z with an axial connection spring k and a torsion release at the top.
    Axial load P: shortening = P h/EA + P/k.  Torque at the top node: carried by nothing but ...
    here a second, unreleased column shares the node, so the torque goes entirely to it."""
    h, P, T, k = 3000.0, 6.0e4, 4.0e5, 9.0e4
    model = build_model(document(
        {"A": (0, 0, 0), "B": (0, 0, h), "C": (0, 0, 2 * h)},
        {"lower": {"i": "A", "j": "B", "releases": {"j": {"ux": {"spring": k}, "rx": True}}}, "upper": ("C", "B")},
        {"A": "fixed", "C": {"dofs": {"rx": "fixed", "ry": "fixed", "rz": "fixed", "ux": "fixed", "uy": "fixed"}}},
        {"L": [{"type": "nodal", "node": "B", "Fz": -P, "Mz": T}]}))
    r = solve("opensees", model, "L")
    assert r.ok and r.fidelity.level == EXACT
    assert r.displacements[1, 2] == pytest.approx(-(P * h / EA + P / k), rel=1e-10)
    assert r.displacements[1, 5] == pytest.approx(T * h / GJ, rel=1e-10)        # twist of the upper column only
    assert r.reactions[0] == pytest.approx([0, 0, P, 0, 0, 0], abs=1e-8 * P)
    assert r.reactions[2, 5] == pytest.approx(-T, rel=1e-10)


@needs_opensees
def test_opensees_skew_axial_spring_is_flagged_approximate_with_a_valid_estimate():
    """The same connection on a skew member needs penalty springs: flagged, and the estimate bounds the
    difference to the exact answer (shortening of bar + spring in series under an axial tip load)."""
    d = np.array([2.0, 1.0, 2.0]) / 3.0
    Lm, P, k = 3000.0, 5.0e4, 7.0e4
    model = build_model(document({"A": (0, 0, 0), "B": tuple(Lm * d)},
                                 {"m": {"i": "A", "j": "B", "releases": {"i": {"ux": {"spring": k}}}}}, {"A": "fixed"},
                                 {"L": [{"type": "memberPoint", "member": "m", "direction": "x", "value": P,
                                         "position": 1.0}]}))
    r = solve("opensees", model, "L")
    assert r.ok
    item = r.fidelity.items["releases.partial"]
    assert item.level == APPROXIMATE and item.control == "penalty" and item.error is not None
    u, _ = local_node(model, r, "m", "B")
    exact = P * Lm / EA + P / k
    assert abs(u[0] - exact) / exact <= max(item.error, 1e-9)
    assert abs(u[0] - exact) / exact < 1e-5
    assert "convergence" in r.info and "penalty" in r.info["convergence"]


# ====================================================================================== second order
def _beam_column(ratio=0.4):
    Lm = 3000.0
    Pe = math.pi ** 2 * EIz / (4 * Lm ** 2)
    P, H = ratio * Pe, 5.0e3
    model = build_model(document(
        {"A": (0, 0, 0), "B": (0, 0, Lm)}, {"c": {"i": "A", "j": "B"}}, {"A": "fixed"},
        {"L": [{"type": "nodal", "node": "B", "Fz": -P, "Fx": H},
               {"type": "memberDistributed", "member": "c", "direction": "X", "w1": 1.0, "w2": 3.0}]},
        combos={"L": {"factors": {"L": 1.0}, "analysisType": "secondOrder"}}, analysis={"alpha": 1.0}))
    # column axes: y = +X, bending about z -> EIz
    s = Beam(Lm, EIz, N=-P, second_order=True).support(0.0, v="fixed", r="fixed").point(Lm, P=H) \
        .dist(0.0, Lm, 1.0, 3.0).solve(tol=1e-11)
    return model, s, Lm, P, H


@pytest.mark.parametrize("solver", SOLVER_PARAMS)
def test_second_order_beam_column_against_oracle(solver):
    """Cantilever beam-column at 40 % of its buckling load: each solver must reproduce the ODE solution
    within its own estimated error; the amplification is far above that error."""
    model, s, Lm, P, H = _beam_column()
    r = solve(solver, model, "L")
    assert r.ok, r.error
    est = r.fidelity.estimated_error
    assert r.fidelity.level == APPROXIMATE and est is not None and est < 2e-3
    tip, M0 = s.v(Lm), s.M(0.0, +1)
    first = Beam(Lm, EIz).support(0.0, v="fixed", r="fixed").point(Lm, P=H).dist(0.0, Lm, 1.0, 3.0).solve()
    assert tip / first.v(Lm) > 1.5                           # strong second-order effect
    assert rel(r.displacements[1, 0], tip) <= est + 1e-7
    f = r.end_forces["c"]
    assert rel(f[5], -M0) <= est + 1e-7                      # base moment H L + ... + P delta
    assert f[0] == pytest.approx(P, rel=1e-5)                # axial force on the member at i (compression)
    # end forces along the UNDEFORMED axes: base shear = total horizontal load
    W = H + 0.5 * (1.0 + 3.0) * Lm
    assert f[1] == pytest.approx(-W, rel=1e-3 if solver == "pynite" else 1e-9)
    assert r.reactions[0, 0] == pytest.approx(-W, rel=1e-3 if solver == "pynite" else 1e-9)
    assert "subdivisions" in r.info["convergence"]


@needs_opensees
def test_opensees_subdivision_convergence_is_second_order():
    from benchfea.crosscheck import convergence_study

    model, s, Lm, P, H = _beam_column()
    study = convergence_study("opensees", model, "L", "subdivisions", [2, 4, 8, 64])
    d = [row["difference"] for row in study]
    assert d[-1] == 0.0 and d[0] > d[1] > d[2] > 0.0
    assert 3.0 < d[0] / d[1] < 5.0 and 3.0 < d[1] / d[2] < 5.5     # error ~ 1/n^2


@pytest.mark.skipif(not __import__("crosscheck_helpers").HAVE["pynite"], reason="PyNiteFEA is not installed")
def test_pynite_subdivision_converges_faster_than_second_order():
    """Consistent geometric stiffness: the differences to a fine solution fall by more than 8 per
    halving of the sub-member length, so the 1/n^2 law used for the estimate is conservative."""
    from benchfea.crosscheck import convergence_study

    model, s, Lm, P, H = _beam_column()
    study = convergence_study("pynite", model, "L", "subdivisions", [1, 2, 4, 16], in_process=True)
    d = [row["difference"] for row in study]
    assert d[0] / d[1] > 8.0 and d[1] / d[2] > 8.0 and d[2] > 0.0


@needs_both
def test_second_order_portal_cross_check_between_solvers():
    """Pinned-base portal with gravity on the columns and a lateral load: OpenSees (chord-rotation Kg,
    iterated N, 16 elements) and PyNite (consistent Kg, first-order N, 4 sub-members) are different
    formulations; they must agree within the sum of their estimated errors."""
    h, b = 4000.0, 6000.0
    model = build_model(document(
        {"a": (0, 0, 0), "A": (0, 0, h), "B": (b, 0, h), "b": (b, 0, 0)},
        {"c1": ("a", "A"), "g": ("A", "B"), "c2": ("b", "B")},
        {"a": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}},
         "b": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}}},
        {"L": [{"type": "nodal", "node": "A", "Fz": -3e5, "Fx": 2e4}, {"type": "nodal", "node": "B", "Fz": -3e5},
               {"type": "memberDistributed", "member": "g", "direction": "Z", "w1": -10.0}]},
        combos={"L": {"factors": {"L": 1.0}, "analysisType": "secondOrder"},
                "F": {"factors": {"L": 1.0}, "analysisType": "firstOrder"}}, analysis={"alpha": 1.0}))
    rep = cross_check(model, "L", options={"pynite": {}}, foreign={"pynite": solve("pynite", model, "L")})
    assert rep.cross_checked, rep.summary()
    ro = rep.results["opensees"]
    assert rep.comparisons[0].max_relative < 2e-3
    first = solve("opensees", model, "F")
    sway2, sway1 = ro.displacements[1, 0], first.displacements[1, 0]
    assert sway2 / sway1 > 1.2                               # the P-Delta effect is what is being compared
    assert rep.equilibrium["opensees"]["force"] < 1e-9
    assert rep.equilibrium["pynite"]["force"] < 5e-3         # PyNite's P-Delta reactions are not equilibrated


@pytest.mark.parametrize("solver", SOLVER_PARAMS)
def test_alpha_scales_loads_and_results(solver):
    """alpha = 1.6: the loads are multiplied by 1.6 for the analysis and the results divided by 1.6;
    the result equals the alpha = 1 analysis of 1.6 times the loads, divided by 1.6."""
    def model(alpha, factor):
        Lm = 3000.0
        P = 0.25 * math.pi ** 2 * EIz / (4 * Lm ** 2)
        return build_model(document(
            {"A": (0, 0, 0), "B": (0, 0, Lm)}, {"c": ("A", "B")}, {"A": "fixed"},
            {"L": [{"type": "nodal", "node": "B", "Fz": -P, "Fx": 4e3}]},
            combos={"L": {"factors": {"L": factor}, "analysisType": "secondOrder"}}, analysis={"alpha": alpha}))
    a = solve(solver, model(1.6, 1.0), "L", estimate_errors=False)
    b = solve(solver, model(1.0, 1.6), "L", estimate_errors=False)
    assert a.ok and b.ok and a.info["alpha"] == 1.6
    assert a.displacements == pytest.approx(b.displacements / 1.6, rel=1e-9, abs=1e-12)
    assert a.reactions == pytest.approx(b.reactions / 1.6, rel=1e-9, abs=1e-6)


@needs_opensees
def test_opensees_second_order_with_shear_deformation_is_unsupported():
    model = build_model(document(
        {"A": (0, 0, 0), "B": (0, 0, 3000)}, {"c": ("A", "B")}, {"A": "fixed"},
        {"L": [{"type": "nodal", "node": "B", "Fz": -1e5, "Fx": 4e3}]},
        combos={"L": {"factors": {"L": 1.0}, "analysisType": "secondOrder"}},
        analysis={"alpha": 1.0, "shearDeformation": True}))
    r = solve("opensees", model, "L", estimate_errors=False)
    assert "shearDeformation" in [i.feature for i in r.fidelity.unsupported]
    rep = cross_check(model, "L", solvers=("opensees",), foreign={"opensees": r})
    assert not rep.cross_checked and any("cannot represent" in x for x in rep.reasons)
