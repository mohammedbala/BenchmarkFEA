"""The beam-column piece (``benchfea.solve.beamcolumn``) on its own: limits, symmetry, the consistent
geometric stiffness as the derivative of the exact stiffness, the transfer with a varying axial force
against an independent integration, the subdivision and its load / station mapping."""
import math

import numpy as np
import pytest
import scipy.linalg

from benchfea.model import build_model
from benchfea.solve import build_element
from benchfea.solve.beam import SpanLoads
from benchfea.solve.beamcolumn import BeamColumn, _Plane, consistent_geometric_stiffness, expm
from benchfea.solve.loads import combination_loads
from benchfea.solve.refine import piece_stations, refine, split_loads

from second_order_reference import beam_column_bvp
from solve_reference import EIy, EIz, GAy, make

L = 3000.0
RELEASES = {"i": {"rz": True, "ry": {"spring": 2.0e9}}, "j": {"rz": {"spring": 5.0e9}, "ux": {"spring": 4.0e5}}}


def member(shear=True, **m):
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0)}, {"m": {"i": "a", "j": "b", **m}}, {"a": "fixed"},
               {"c": [{"type": "nodal", "node": "b", "Fx": 1.0}]}, analysis={"shearDeformation": shear})
    model = build_model(doc)
    return model, model.members["m"]


@pytest.mark.parametrize("shear", [False, True])
@pytest.mark.parametrize("m", [{}, {"releases": RELEASES}, {"offsets": {"i": [100, 0, 50], "j": [-80, 30, 0]}},
                               {"type": "truss"}], ids=["plain", "releases", "offsets", "truss"])
def test_without_axial_force_the_piece_is_the_first_order_element(m, shear):
    """Stiffness, fixed-end forces and station values of the exact first-order element (THEORY 1-5),
    reached by a completely different route (transfer matrix instead of cantilever flexibility)."""
    model, mem = member(shear, **m)
    first, piece = build_element(model, mem), BeamColumn.build(model, mem)
    scale = np.abs(first.k_global).max()
    assert np.abs(piece.k_global - first.k_global).max() <= 1e-12 * scale
    assert np.abs(piece.k_global - piece.k_global.T).max() <= 1e-13 * scale
    h = mem.length
    span = SpanLoads(points=((0.0, np.array([0.0, 5.0, 3.0, 0.0, 2000.0, 1000.0])),
                             (h, np.array([1.0, 2.0, 3.0, 0.0, 5.0, 6.0]))),
                     dists=((0.0, h, np.array([0.0, 1.0, 2.0]), np.array([0.0, -1.0, 4.0])),))
    loaded = BeamColumn.build(model, mem, span)
    f_ref = first.fixed_end(span)
    assert loaded.f_fixed == pytest.approx(f_ref, rel=1e-10, abs=1e-10 * np.abs(f_ref).max())
    u = np.random.default_rng(3).normal(size=12) * 1e-3
    stations = [(0.0, -1), (0.0, 1), (0.3 * h, -1), (h, -1), (h, 1)]
    a, b = first.recover(u, span, stations), loaded.recover(u, stations, checks=True)
    assert b.forces == pytest.approx(a.forces, rel=1e-9, abs=1e-9 * np.abs(a.forces).max())
    skip = [0, 1, 2, 4, 5] if m.get("type") == "truss" else slice(None)       # a truss has no twist of its own
    assert b.displacements[:, skip] == pytest.approx(a.displacements[:, skip], rel=1e-9, abs=1e-12)
    assert b.deformed_residual <= 1e-10 * b.moment_scale


@pytest.mark.parametrize("axial", [-2.0e5, 3.0e6])
@pytest.mark.parametrize("m", [{}, {"releases": RELEASES}, {"offsets": {"i": [100, 0, 50], "j": [-80, 30, 0]}}],
                         ids=["plain", "releases", "offsets"])
def test_stiffness_with_axial_force_is_symmetric_and_keeps_translations_free(m, axial):
    model, mem = member(True, **m)
    piece = BeamColumn.build(model, mem, n_start=axial)
    k = piece.k_global
    assert piece.stable and np.abs(k - k.T).max() <= 1e-12 * np.abs(k).max()
    for a in range(3):                                   # a rigid translation produces no force
        t = np.zeros(12)
        t[[a, 6 + a]] = 1.0
        assert np.abs(k @ t).max() <= 1e-9 * np.abs(k).max()
    spin = np.zeros(12)                                  # ... a rigid rotation about a transverse axis does
    spin[[4, 10]] = 1.0                                  # rotation about global Y through node a
    spin[8] = -mem.length - (m.get("offsets", {}).get("j", [0])[0] - m.get("offsets", {}).get("i", [0])[0])
    if not m:
        force = k @ spin
        assert force[8] == pytest.approx(axial * spin[8] / mem.length, rel=1e-9)   # N x chord rotation


@pytest.mark.parametrize("m", [{}, {"releases": RELEASES}, {"offsets": {"i": [100, 0, 50], "j": [-80, 30, 0]}},
                               {"type": "truss"}], ids=["plain", "releases", "offsets", "truss"])
def test_consistent_geometric_stiffness_is_the_derivative_of_the_exact_stiffness(m):
    """kg = integral N (v_a' v_b' + w_a' w_b') on the first-order field  ==  dK_T/dN at N = 0."""
    model, mem = member(True, **m)
    kg = consistent_geometric_stiffness(BeamColumn.build(model, mem), (1.0, 0.0, 0.0, 1.0, 1.0))
    dn = 2.0e3
    fd = (BeamColumn.build(model, mem, n_start=dn).k_global - BeamColumn.build(model, mem, n_start=-dn).k_global) / (2 * dn)
    assert kg == pytest.approx(fd, rel=1e-6, abs=1e-7 * np.abs(kg).max())
    assert np.abs(kg - kg.T).max() <= 1e-12 * np.abs(kg).max()


def test_piece_beyond_its_own_critical_load_is_flagged():
    model, mem = member(False, releases={"i": {"ry": True, "rz": True}, "j": {"ry": True, "rz": True}})
    euler = math.pi ** 2 * EIy / mem.length ** 2          # pin-ended, weak plane (value used as a load level only)
    assert BeamColumn.build(model, mem, n_start=-0.95 * euler).stable
    assert not BeamColumn.build(model, mem, n_start=-1.05 * euler).stable


def test_matrix_exponential_against_scipy():
    rng = np.random.default_rng(0)
    for scale in (0.05, 1.0, 12.0, 40.0):
        A = rng.normal(size=(6, 6)) * scale
        A[4:] = 0.0
        ref = scipy.linalg.expm(A)
        assert np.abs(expm(A) - ref).max() <= 1e-11 * np.abs(ref).max()


@pytest.mark.parametrize("beta", [0.0, 1.0 / GAy])
def test_transfer_with_a_varying_axial_force_is_fourth_order_accurate(beta):
    """Initial value problem over one piece with N(x) = n0 - na x - ns x^2 / 2 and a trapezoidal load,
    against solve_bvp; halving the piece divides the error by about 2^5 per piece (local error)."""
    n0, na, ns = -1.5e5, -60.0, 0.02
    z0 = np.array([0.3, 1.0e-3, 2.0e6, 500.0])
    errors = []
    for h in (2000.0, 1000.0):
        plane = _Plane(h, EIz, beta, (n0, na, ns), 0.7, 3.0e-4)
        ref = beam_column_bvp(h, EIz, 1.0 / beta if beta else None, lambda x: n0 - na * x - 0.5 * ns * x * x,
                              lambda x: 0.7 + 3.0e-4 * x,
                              lambda s, e: [s[0, 0] - z0[0], s[0, 1] - z0[1], s[0, 2] - z0[2], s[0, 3] - z0[3]], tol=1e-12)
        z, r = plane.state(z0, 1.0), ref(h)
        errors.append(np.abs((z - r) / r).max())
        zi, ri = plane.state(z0, 0.37), ref(0.37 * h)     # inside the piece
        assert np.abs((zi - ri) / ri).max() <= 1e-7
    assert errors[0] <= 1e-6 and errors[1] <= errors[0] / 12.0


def test_subdivision_keeps_nodes_loads_and_stations_consistent():
    doc = make({"a": (0, 0, 0), "b": (L, 0, 0), "c": (L, 0, 2000)},
               {"m": {"i": "a", "j": "b", "releases": {"j": {"rz": True}}, "offsets": {"i": [100, 0, 0]}},
                "t": {"i": "b", "j": "c", "type": "tensionOnly"}}, {"a": "fixed", "c": "pinned"},
               {"c": [{"type": "memberPoint", "member": "m", "direction": "Z", "value": -900.0, "position": 0.25},
                      {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -1.0, "w2": -3.0,
                       "start": 0.25, "end": 0.75},
                      {"type": "memberPoint", "member": "m", "direction": "Z", "value": 50.0, "position": 0.0}]})
    model = build_model(doc)
    mesh = refine(model, {"m": 4, "t": 3}, unsplit={"t"})
    m = model.members["m"]
    h = m.length
    assert mesh.pieces["t"] == ["t"] and list(mesh.model.nodes)[:3] == ["a", "b", "c"]          # original nodes first
    assert mesh.bounds["m"] == pytest.approx([0, 0.25 * h, 0.5 * h, 0.75 * h, h])              # breaks are boundaries
    pieces = [mesh.model.members[p] for p in mesh.pieces["m"]]
    assert pieces[0].offset_i.any() and not pieces[1].offset_i.any() and not pieces[0].offset_j.any()
    assert pieces[-1].k_end_j[5] == 0.0 and pieces[0].k_end_j[5] == math.inf and pieces[-1].k_end_i[5] == math.inf
    assert sum(p.length for p in pieces) == pytest.approx(h)
    loads = combination_loads(model, model.combinations["c"])
    split = split_loads(mesh, loads)
    total = sum((split.member[p].resultant()[:3] for p in mesh.pieces["m"] if p in split.member), np.zeros(3))
    assert total == pytest.approx(loads.member["m"].resultant()[:3], rel=1e-12)
    first = split.member[mesh.pieces["m"][0]]
    assert sorted(x for x, _ in first.points) == pytest.approx([0.0, 0.25 * h])               # start load, load at the break
    where = piece_stations(mesh.bounds["m"], [(0.0, -1), (0.25 * h, -1), (0.25 * h, 1), (0.6 * h, -1), (h, 1)])
    assert [[(n, round(x, 6), s) for n, x, s in w] for w in where] == [
        [(0, 0.0, -1), (1, round(0.25 * h, 6), -1)], [(2, 0.0, 1)], [(3, round(0.1 * h, 6), -1)], [(4, round(0.25 * h, 6), 1)]]
