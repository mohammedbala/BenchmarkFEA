"""Independent review, second order.

References
* tableau (reviewer's own): the member ODE with the axial force acting on the displaced slice, solved
  through exp(A L) of the full 12-state system -- the same *theory* as docs/THEORY.md 10 (P-Delta and
  P-delta, both planes, shear deformation in the Engesser form, rigid-arm term), a different method.
  With a constant axial force per member the reference is exact: 1e-9 is demanded (measured 5e-13).
  With an axial force that varies along a member the reference uses a piecewise constant N (midpoint
  rule, error O(h^2)); its convergence towards the engine at that order is what is asserted.
* OpenSees ``PDelta`` with 32 elements per member: chord (string) geometric stiffness only, no P-delta
  inside an element, error O(h^2) -- bound 1e-3 (measured <= 2.1e-4).  End shears are not compared (OpenSees
  reports them normal to the element chord).  ``-jntOffset`` carries no P-Delta of the rigid arm
  (measured: 25 % low on S05), so offsets are checked with stiff stub members instead.
* PyNite ``analyze_PDelta`` with 8 sub-members: consistent cubic Kg, ONE pass with first-order axial
  forces, an extra axial Kg term -- bound 1e-3 (measured <= 1.7e-4).
"""
import copy

import numpy as np
import pytest

import oracles
from benchfea.solve import ModelRefused
from harness import compare, compare_foreign, doc_arrays, nodes, member, rel, run_benchfea, worst
from models import SECOND, sec
from tableau import Tableau

CONSTANT_N = ["S01_cantilever_beam_column_both_planes", "S02_cantilever_beam_column_with_shear_deformation_and_tension",
              "S03_portal_with_leaning_column", "S05_offsets_along_the_axis_carrying_compression",
              "S08_diaphragm_storey_with_leaning_columns"]


def _numbers(out):
    return {k: v for k, v in out.items() if not k.startswith("_")}


@pytest.mark.parametrize("name", CONSTANT_N)
def test_second_order_against_the_tableau_reference(name):
    raw, cid = SECOND[name]()
    _, res, doc = run_benchfea(raw)
    c = doc["combinations"][cid]
    assert c["converged"] and c["invariants"]["passed"]
    out = compare(raw, cid, second_order=True, doc=doc)
    assert worst(out) <= 1.0e-9, _numbers(out)
    # the second-order effect is large in these models: the test is not vacuous
    first = copy.deepcopy(raw)
    first["combinations"][cid]["analysisType"] = "firstOrder"
    _, _, d1 = run_benchfea(first)
    U2, U1 = doc_arrays(raw, c)[0], doc_arrays(raw, d1["combinations"][cid])[0]
    amplification = np.abs(U2[:, :3]).max() / np.abs(U1[:, :3]).max()
    assert abs(amplification - 1.0) > 0.1, amplification


@pytest.mark.parametrize("name, n, bound", [("S04_skew_space_frame_y_up_kip_in", 8, 3.0e-6),
                                            ("S06_truss_members_with_transverse_load_and_partial_releases", 16, 2.0e-5),
                                            ("S07_columns_under_self_weight_varying_axial_force", 16, 1.0e-3)])
def test_second_order_with_varying_axial_force_reference_converges_to_the_engine(name, n, bound):
    """Members with axial load components (inclined beams under gravity, a loaded truss member, columns
    under self weight at an amplification of 3).  The reference (piecewise constant N) must approach the
    engine's result at (at least) its own order: halving its step divides the difference by 4 or more
    (measured ratios 4.0 for the self-weight columns, 7 and more where the load breaks make it better)."""
    raw, cid = SECOND[name]()
    _, res, doc = run_benchfea(raw)
    coarse = worst(compare(raw, cid, second_order=True, doc=doc, nsub=n))
    fine = worst(compare(raw, cid, second_order=True, doc=doc, nsub=2 * n))
    assert fine <= bound, (coarse, fine)
    assert coarse / fine >= 3.0, (coarse, fine)


@pytest.mark.parametrize("solver, name, nsub", [
    ("opensees", "S03_portal_with_leaning_column", 32), ("opensees", "S04_skew_space_frame_y_up_kip_in", 32),
    ("opensees", "S08_diaphragm_storey_with_leaning_columns", 32),
    ("pynite", "S01_cantilever_beam_column_both_planes", 8), ("pynite", "S03_portal_with_leaning_column", 8),
    ("pynite", "S04_skew_space_frame_y_up_kip_in", 8)])
def test_second_order_against_opensees_and_pynite(solver, name, nsub):
    if not oracles.have(solver):
        pytest.skip(f"{solver} is not installed")
    raw, cid = SECOND[name]()
    _, res, doc = run_benchfea(raw)
    f = (oracles.run_opensees if solver == "opensees" else oracles.run_pynite)(raw, cid, second_order=True, nsub=nsub)
    assert f is not None
    for q, v in compare_foreign(raw, doc["combinations"][cid], f, shears=False).items():
        assert v <= 1.0e-3, f"{name} vs {solver}: {q} differs by {v:.2e}"


def test_rigid_offsets_carrying_compression_equal_stiff_stub_members():
    """S05 (amplification 3): the rigid end offsets replaced by members 100 times larger (EA x 1e2,
    EI x 1e4).  Engine with offsets = engine with stubs = OpenSees PDelta with stubs to 1e-3; this is
    the second-order term of the rotating rigid arm (THEORY 10.7)."""
    raw, cid = SECOND["S05_offsets_along_the_axis_carrying_compression"]()
    _, _, doc = run_benchfea(raw)
    U = doc_arrays(raw, doc["combinations"][cid])[0]
    stub = copy.deepcopy(raw)
    stub["analysis"]["convergenceTolerance"] = 1.0e-7
    stub["sections"]["R"] = sec(("kN", "m"), scale=100.0)
    nodes(stub, a1=(0, 0, 0.6), A1=(0, 0, 4.6), B1=(6, 0, 4.5), A2=(0.3, 0, 5), B2=(5.7, 0, 5))
    stub["members"] = {}
    member(stub, "s1", "a", "a1", sec="R")
    member(stub, "c1", "a1", "A1")
    member(stub, "s2", "A1", "A", sec="R")
    member(stub, "c2", "b", "B1", orientation={"roll": 90})
    member(stub, "s3", "B1", "B", sec="R", orientation={"roll": 90})
    member(stub, "s4", "A", "A2", sec="R")
    member(stub, "g", "A2", "B2", sec="B")
    member(stub, "s5", "B2", "B", sec="R")
    _, _, doc2 = run_benchfea(stub)
    ids, ids2 = list(raw["nodes"]), list(stub["nodes"])
    U2 = doc_arrays(stub, doc2["combinations"][cid])[0][[ids2.index(n) for n in ids]]
    assert rel(U[:, :3], U2[:, :3]) <= 1.0e-3 and rel(U[:, 3:], U2[:, 3:]) <= 1.0e-3
    if oracles.have("opensees"):
        f = oracles.run_opensees(stub, cid, second_order=True, nsub=32)
        Uo = f["displacements"][[ids2.index(n) for n in ids]]
        assert rel(U[:, :3], Uo[:, :3]) <= 1.0e-3 and rel(U[:, 3:], Uo[:, 3:]) <= 1.0e-3


def test_off_axis_offsets_follow_the_documented_arm_rule():
    """Offsets not along the member (THEORY 10.7, 16: only the axial end force acts on the rotating
    arm; symmetric form 1/2 (o f^T + f o^T) - (o.f) I).  The reference with the same rule agrees; without
    the term it does not, so the term is really there."""
    from models import LINEAR, _so
    raw, cid = LINEAR["L05_offsets_off_axis_skew"]()
    loads = raw["loadCases"]["W"]["loads"]
    raw["loadCases"]["W"]["loads"] = [ld for ld in loads if ld["type"] == "nodal"] + \
        [{"type": "nodal", "node": n, "Fz": -900.0} for n in "BCF"]
    _so(raw)
    _, _, doc = run_benchfea(raw)
    with_arm = worst(compare(raw, cid, second_order=True, doc=doc))
    without = worst(compare(raw, cid, second_order=True, doc=doc, arm=False))
    assert with_arm <= 1.0e-9, with_arm
    assert without >= 1.0e-3, without


def _reference(raw, cid, doc, scale=1.0, second=True):
    """Largest relative difference of displacements, reactions and end forces between the engine's
    output ``doc`` and the tableau solution of ``raw`` multiplied by ``scale``."""
    c = doc["combinations"][cid]
    ref = Tableau(raw).solve(cid, second_order=second)
    U, R, E = doc_arrays(raw, c)
    Er = np.array([ref["members"][m]["end_forces"] for m in E])
    Eb = np.array([E[m] for m in E])
    return max(rel(U, ref["displacements"] * scale), rel(R, ref["reactions"] * scale), rel(Eb, Er * scale))


def test_alpha_multiplies_the_loads_and_divides_the_results():
    raw, cid = SECOND["S03_portal_with_leaning_column"]()
    raw["analysis"]["alpha"] = 1.6
    _, _, doc = run_benchfea(raw)
    ref = copy.deepcopy(raw)
    ref["combinations"][cid]["factors"] = {"U": 1.6}
    assert _reference(ref, cid, doc, 1 / 1.6) <= 1.0e-9
    assert doc["combinations"][cid]["alpha"] == 1.6


def test_stiffness_reduction_is_the_analysis_of_the_reduced_structure():
    """AISC 360-22 C2.3: 0.8 on all stiffnesses and tau_b on the flexural stiffness -- here as a given
    number, as a member subset, and computed (Eq. C2-2a/b with P_ns = F_y A)."""
    raw, cid = SECOND["S03_portal_with_leaning_column"]()
    raw["combinations"][cid]["stiffnessReduction"] = True
    raw["analysis"]["stiffnessReduction"] = {"factor": 0.8, "tauB": 0.9}
    _, _, doc = run_benchfea(raw)
    ref = copy.deepcopy(raw)
    ref["materials"]["M"]["E"] *= 0.8
    ref["materials"]["M"]["G"] *= 0.8
    for s in ref["sections"].values():
        s["Iy"] *= 0.9
        s["Iz"] *= 0.9
    assert _reference(ref, cid, doc) <= 1.0e-9

    raw["analysis"]["stiffnessReduction"] = {"factor": 0.7, "tauB": "none", "members": ["c1", "c2", "lean"]}
    _, _, doc = run_benchfea(raw)
    ref = copy.deepcopy(raw)
    ref["materials"]["R"] = {"E": 2.0e8 * 0.7, "G": 7.7e7 * 0.7}
    for m in ("c1", "c2", "lean"):
        ref["members"][m]["material"] = "R"
    assert _reference(ref, cid, doc) <= 1.0e-9

    raw["analysis"].pop("stiffnessReduction")
    raw["materials"]["M"]["Fy"] = 120000.0          # P_ns = 1440 kN; the leaning column carries about 890 kN
    _, _, doc = run_benchfea(raw)
    c = doc["combinations"][cid]
    tau = c["tauB"]
    assert tau["c1"] == tau["c2"] == tau["g1"] == 1.0 and 0.9 < tau["lean"] < 0.97
    for m, t in tau.items():
        ratio = -min(c["members"][m]["stations"]["N"]) / (120000.0 * raw["sections"][raw["members"][m]["section"]]["A"])
        expected = 1.0 if ratio <= 0.5 else 4.0 * ratio * (1.0 - ratio)       # AISC 360-22 Eq. C2-2a / C2-2b
        assert abs(t - expected) <= 1.0e-9
    ref = copy.deepcopy(raw)
    ref["materials"]["M"]["E"] *= 0.8
    ref["materials"]["M"]["G"] *= 0.8
    for m, t in tau.items():
        s = copy.deepcopy(raw["sections"][raw["members"][m]["section"]])
        s["Iy"] *= t
        s["Iz"] *= t
        ref["sections"]["S_" + m] = s
        ref["members"][m]["section"] = "S_" + m
    assert _reference(ref, cid, doc) <= 1.0e-9


@pytest.mark.parametrize("second", [True, False])
def test_notional_loads_are_the_stated_nodal_loads(second):
    """AISC 360-22 C2.2b Eq. C2-1 as the engine defines it: coefficient x gravity load at each node,
    member loads lumped by the lever rule (here by hand: the girder's 160 kN go half to A, half to B)."""
    raw, cid = SECOND["S03_portal_with_leaning_column"]()
    raw["combinations"][cid]["notionalLoad"] = {"coefficient": 0.003, "direction": [1, 0, 0]}
    if not second:
        raw["combinations"][cid]["analysisType"] = "firstOrder"
    _, _, doc = run_benchfea(raw)
    ref = copy.deepcopy(raw)
    ref["loadCases"]["U"]["loads"] += [{"type": "nodal", "node": n, "Fx": 0.003 * y}
                                       for n, y in (("A", 380.0), ("B", 480.0), ("C", 900.0))]
    assert _reference(ref, cid, doc, second=second) <= 1.0e-9


def test_refusal_beyond_the_critical_load_reports_the_buckling_factor_of_the_reference():
    """A pinned-base portal whose columns can buckle out of plane as cantilevers (Pe = 247 kN) under
    300 kN each: the engine refuses the second-order combination; its factor is the first root of the
    reference's determinant (0.709).  (This is the model of
    tests/test_crosscheck_threeway.py::test_second_order_three_way_when_the_engine_supports_it.)"""
    from models import start
    from harness import case
    d = start(("N", "mm"))
    d["sections"]["S"] = {"type": "properties", "A": 5000.0, "Iy": 8.0e6, "Iz": 4.0e7, "J": 3.0e5}
    nodes(d, a=(0, 0, 0), A=(0, 0, 4000.0), B=(6000.0, 0, 4000.0), b=(6000.0, 0, 0))
    for n in "ab":
        d["supports"][n] = {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}}
    member(d, "c1", "a", "A")
    member(d, "g", "A", "B")
    member(d, "c2", "b", "B")
    case(d, "L", {"type": "nodal", "node": "A", "Fz": -3e5, "Fx": 2e4}, {"type": "nodal", "node": "B", "Fz": -3e5},
         {"type": "memberDistributed", "member": "g", "direction": "Z", "w1": -10.0}, analysisType="secondOrder")
    with pytest.raises(ModelRefused) as e:
        run_benchfea(d)
    assert [p.code for p in e.value.problems] == ["stability.beyond-critical-load"]
    (refused,) = e.value.refused.values()
    root = Tableau(d).buckling("L", n_modes=1, lam_max=2.0, steps=40)[0]
    assert abs(root - 0.70905) < 1e-4
    # the engine's factor scales the axial forces of its last iterate, the reference the first-order ones
    assert abs(refused.info["bucklingFactor"] - root) <= 2.0e-3 * root
