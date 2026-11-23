"""Independent review: sign and axis conventions of docs/CONVENTIONS.md checked by direct equilibrium,
with no solver at all on the reference side.

A statically determinate tree (one fixed support, skew members with rolls) is cut at every station:
the stress resultants on the +x face are the resultant of all loads beyond the cut, expressed in the
member's local axes (triad from ``tableau.triad``, the reviewer's reading of CONVENTIONS section 2).
"""
import numpy as np
import pytest

from benchfea.model import build_model
from harness import case, member, nodes, run_benchfea
from models import start
from tableau import triad


def _tree(vertical):
    d = start(("kN", "m"), vertical)
    nodes(d, A=(0, 0, 0), B=(1.0, 2.0, 3.0), C=(4.0, 3.0, 2.5), D=(0.5, 5.0, 4.0))
    d["supports"] = {"A": {"type": "fixed"}}
    member(d, "ab", "A", "B", orientation={"roll": 20})
    member(d, "bc", "B", "C", orientation={"roll": -35})
    member(d, "bd", "B", "D", orientation={"referenceVector": [1, 0, 0]})
    case(d, "L", {"type": "nodal", "node": "C", "Fx": 3, "Fy": -4, "Fz": -6, "Mx": 1, "My": 2, "Mz": -1.5},
         {"type": "nodal", "node": "D", "Fx": -2, "Fz": 5, "My": 0.7},
         {"type": "memberDistributed", "member": "bc", "direction": "Z", "w1": -2.0},
         {"type": "memberDistributed", "member": "bd", "direction": "X", "w1": 1.5},
         {"type": "memberPoint", "member": "ab", "direction": "Y", "value": 4.0, "position": 0.5})
    return d


def _loads_beyond(d, mid, x):
    """[(point, force, moment)] of everything outboard of the cut at distance x on member ``mid``."""
    X = {n: np.array([v["x"], v["y"], v["z"]], float) for n, v in d["nodes"].items()}
    ends = {"ab": ("A", "B"), "bc": ("B", "C"), "bd": ("B", "D")}
    out = []

    def nodal(n):
        for ld in d["loadCases"]["L"]["loads"]:
            if ld["type"] == "nodal" and ld["node"] == n:
                out.append((X[n], np.array([ld.get(k, 0.0) for k in ("Fx", "Fy", "Fz")]),
                            np.array([ld.get(k, 0.0) for k in ("Mx", "My", "Mz")])))

    def span(m, a):                                    # loads on member m beyond distance a from its start
        i, j = ends[m]
        L = np.linalg.norm(X[j] - X[i])
        ex = (X[j] - X[i]) / L
        for ld in d["loadCases"]["L"]["loads"]:
            if ld.get("member") != m:
                continue
            e = np.eye(3)["XYZ".index(ld["direction"])]
            if ld["type"] == "memberDistributed" and a < L:
                out.append((X[i] + ex * (a + L) / 2, e * ld["w1"] * (L - a), np.zeros(3)))
            elif ld["type"] == "memberPoint" and ld["position"] * L > a:
                out.append((X[i] + ex * ld["position"] * L, e * ld["value"], np.zeros(3)))
    span(mid, x)
    if mid == "ab":
        for m in ("bc", "bd"):
            span(m, 0.0)
        for n in ("B", "C", "D"):
            nodal(n)
    else:
        nodal(ends[mid][1])
    return out, X, ends


@pytest.mark.parametrize("vertical", ["Z", "Y"])
def test_station_forces_are_the_resultant_of_the_loads_beyond_the_cut(vertical):
    d = _tree(vertical)
    _, res, doc = run_benchfea(d)
    c = doc["combinations"]["L"]
    for mid, m in c["members"].items():
        md = d["members"][mid]["orientation"]
        _, X, ends = _loads_beyond(d, mid, 0.0)
        i, j = ends[mid]
        R = triad(X[i], X[j], vertical, md.get("referenceVector"), md.get("roll", 0.0))
        L = np.linalg.norm(X[j] - X[i])
        xs = m["stations"]["x"]
        for k, x in enumerate(xs):
            if 0 < k < len(xs) - 1 and (xs[k - 1] == x or xs[k + 1] == x):
                continue                               # the two sides of a point load are checked by neighbours
            loads, _, _ = _loads_beyond(d, mid, min(x, L))
            cut = X[i] + R[0] * x
            F = sum(f for _, f, _ in loads)
            M = sum(mo + np.cross(p - cut, f) for p, f, mo in loads)
            expected = np.concatenate([R @ F, R @ M])                 # N, Vy, Vz, T, My, Mz on the +x face
            got = np.array([m["stations"][q][k] for q in ("N", "Vy", "Vz", "T", "My", "Mz")])
            assert np.abs(got - expected).max() <= 1e-9 * 40.0, (mid, x, got, expected)
    # reaction = minus the resultant of all loads about the support (forces ON the structure)
    loads, X, _ = _loads_beyond(d, "ab", 0.0)
    F = sum(f for _, f, _ in loads)
    M = sum(mo + np.cross(p - X["A"], f) for p, f, mo in loads)
    got = np.array([c["reactions"]["A"][q] for q in ("Fx", "Fy", "Fz", "Mx", "My", "Mz")])
    assert np.abs(got + np.concatenate([F, M])).max() <= 1e-9 * 40.0
    # end forces: on the member; at i they are minus the station at 0, at j the station at L
    for mid, m in c["members"].items():
        for end, k, sign in (("i", 0, -1.0), ("j", -1, 1.0)):
            for q, f in zip(("N", "Vy", "Vz", "T", "My", "Mz"), ("Fx", "Fy", "Fz", "Mx", "My", "Mz")):
                assert abs(m["endForces"][end][f] - sign * m["stations"][q][k]) <= 1e-9 * 40.0


def test_the_worked_example_of_conventions_section_4():
    """Simply supported beam, y up, uniform load w_y = -q: end forces Fy = +qL/2 at both ends,
    Vy(x) = -qL/2 + q x, Mz(x) = q x (L - x) / 2 >= 0, positive vertical reactions."""
    q, L = 7.0, 6.0
    d = start()
    nodes(d, A=(0, 0, 0), B=(L, 0, 0))
    d["supports"] = {"A": {"type": "pinned", "dofs": {"rx": "fixed"}}, "B": {"type": "roller", "dofs": {"uy": "fixed"}}}
    member(d, "b", "A", "B")
    case(d, "L", {"type": "memberDistributed", "member": "b", "direction": "y", "w1": -q})
    _, res, doc = run_benchfea(d)
    m = doc["combinations"]["L"]["members"]["b"]
    assert abs(m["endForces"]["i"]["Fy"] - q * L / 2) < 1e-9 and abs(m["endForces"]["j"]["Fy"] - q * L / 2) < 1e-9
    x = np.array(m["stations"]["x"])
    assert np.allclose(m["stations"]["Vy"], -q * L / 2 + q * x, atol=1e-9)
    assert np.allclose(m["stations"]["Mz"], q * x * (L - x) / 2, atol=1e-9)
    assert abs(doc["combinations"]["L"]["reactions"]["A"]["Fz"] - q * L / 2) < 1e-9
    # the beam deflects along -y (down); local z = -Y carries nothing
    assert min(m["stations"]["uy"]) < 0.0 and np.allclose(m["stations"]["uz"], 0.0)


def test_default_member_axes_of_conventions_section_2():
    d = start()
    nodes(d, O=(0, 0, 0), X=(3, 0, 0), Y=(0, 3, 0), U=(0, 0, 3), D=(0, 0, -3))
    d["supports"] = {"O": {"type": "fixed"}}
    member(d, "x", "O", "X")
    member(d, "y", "O", "Y")
    member(d, "up", "O", "U")
    member(d, "down", "O", "D")
    member(d, "rolled", "O", "X", orientation={"roll": 90})
    case(d, "L", {"type": "nodal", "node": "X", "Fz": -1})
    model = build_model(d)
    R = {mid: np.asarray(m.R) for mid, m in model.members.items()}
    assert np.allclose(R["x"], [[1, 0, 0], [0, 0, 1], [0, -1, 0]])           # y = +Z, z = -Y
    assert np.allclose(R["y"], [[0, 1, 0], [0, 0, 1], [1, 0, 0]])            # y = +Z, z = +X
    assert np.allclose(R["up"][1], [1, 0, 0]) and np.allclose(R["down"][1], [1, 0, 0])
    assert np.allclose(R["rolled"], [[1, 0, 0], [0, -1, 0], [0, 0, -1]])     # y' = z, z' = -y
    for mid, m in model.members.items():
        o = d["members"][mid].get("orientation", {})
        assert np.allclose(R[mid], triad(m.start, m.end, "Z", o.get("referenceVector"), o.get("roll", 0.0)))
