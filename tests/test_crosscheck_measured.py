"""Measured behaviour of the installed foreign solvers that the translators are built around.

Each test reproduces one behaviour recorded in docs/CROSSCHECK.md (OpenSeesPy 3.8.0, PyNiteFEA 3.2.1)
with a minimal model.  They are not tests of benchfea: if one FAILS, the installed solver behaves
differently from the version the translators were written for, and the corresponding work-around in
``benchfea.crosscheck`` (and its documentation) must be reviewed -- the translators stay correct as long
as they avoid the facility in question.
"""
import numpy as np
import pytest

from benchfea.crosscheck import run_job

from crosscheck_helpers import E, G, SEC, needs_opensees, needs_pynite

A, Iy, Iz, J, Avy, Avz = (SEC[k] for k in ("A", "Iy", "Iz", "J", "Asy", "Asz"))
L = 3000.0
EB = ["element", "elasticBeamColumn", 1, 1, 2, A, E, G, J, Iy, Iz, 1]
TIMO = ["element", "ElasticTimoshenkoBeam", 1, 1, 2, E, G, A, J, Iy, Iz, Avy, Avz, 1]


def _job(model, pattern, second=False, **extra):
    """Beam 1-2 along X; vecxz = -Y gives local y = +Z, local z = -Y."""
    job = {"solver": "opensees", "secondOrder": second, "tolerance": 1e-12, "maxIterations": 50, "system": "UmfPack",
           "model": [["node", 1, 0.0, 0.0, 0.0], ["node", 2, L, 0.0, 0.0], ["fix", 1, 1, 1, 1, 1, 1, 1]] + model,
           "pattern": pattern, "nodes": [1, 2], "elements": [1], "beams": [1]}
    job.update(extra)
    out = run_job(job)
    assert out["ok"], out.get("error")
    return out


LINEAR = ["geomTransf", "Linear", 1, 0.0, -1.0, 0.0]


@needs_opensees
def test_opensees_ignores_end_intensity_of_a_varying_load_in_3d():
    """eleLoad -beamUniform Wy Wz Wx a b Wyb Wzb Wxb: the end values are ignored without a message; the
    response equals that of the uniform load Wy on [a, b] (so a triangular load starting at 0 does
    nothing)."""
    tip = lambda args: _job([LINEAR, EB], [["eleLoad", "-ele", 1, "-type", "-beamUniform", *args]])["disp"]["2"][2]
    uniform = tip([-2.0, 0.0, 0.0, 0.2, 0.8])
    varying = tip([-2.0, 0.0, 0.0, 0.2, 0.8, -7.0, 0.0, 0.0])
    assert varying == pytest.approx(uniform, rel=1e-12)
    assert tip([0.0, 0.0, 0.0, 0.0, 1.0, -7.0, 0.0, 0.0]) == 0.0


@needs_opensees
def test_opensees_release_flags_condense_only_full_span_uniform_loads():
    """Both ends held, -releasez at j.  A full-span uniform load gives no moment at the released end
    (propped cantilever), but a point load or a partial uniform load leaves the clamped-clamped
    fixed-end moment there: the translator never puts eleLoad on a released element."""
    hold = [["fix", 2, 1, 1, 1, 0, 1, 1]]

    def mz_j(load):
        out = _job([LINEAR, EB + ["-releasez", 2]] + hold, [["eleLoad", "-ele", 1, "-type", *load]])
        return out["localForce"]["1"][11]

    q = 2.0
    assert abs(mz_j(["-beamUniform", q, 0.0, 0.0])) < 1e-6 * q * L * L
    a = 0.4 * L
    assert mz_j(["-beamPoint", 1000.0, 0.0, 0.4]) == pytest.approx(1000.0 * a * a * (L - a) / L ** 2, rel=1e-9)
    assert abs(mz_j(["-beamUniform", q, 0.0, 0.0, 0.2, 0.7])) > 1e-2 * q * L * L


@needs_opensees
def test_opensees_timoshenko_element_mishandles_joint_offsets_and_ignores_pdelta():
    """Cantilever, tip force P, rigid offsets 300 at i and 200 at j (flexible length Lf = 2500).
    Correct (hand): at the flexible end act P and the moment 200 P, so with the shear term
    v_node = P Lf^3/3EI + 200 P Lf^2/2EI + P Lf/GA + 200 (P Lf^2/2EI + 200 P Lf/EI), and the moment at
    node 1 is P L.  elasticBeamColumn gives exactly that (without the shear term).
    ElasticTimoshenkoBeam with -jntOffset uses the shortened length but NOT the rigid arms: it returns
    P Lf^3/3EI + P Lf/GA and a node moment P Lf (moment equilibrium violated).  The translator therefore
    gives Timoshenko elements explicit nodes at the flexible ends, tied by rigidLink('beam')."""
    P, Lf, e = 1000.0, 2500.0, 200.0
    EI, GA = E * Iz, G * Avy
    load = [["load", 2, 0.0, 0.0, P, 0.0, 0.0, 0.0]]
    plain = _job([LINEAR, TIMO], load)["disp"]["2"][2]
    assert plain == pytest.approx(P * L ** 3 / (3 * EI) + P * L / GA, rel=1e-12)
    off = LINEAR + ["-jntOffset", 300.0, 0.0, 0.0, -e, 0.0, 0.0]
    bending = P * Lf ** 3 / (3 * EI) + e * P * Lf ** 2 / (2 * EI) + e * (P * Lf ** 2 / (2 * EI) + e * P * Lf / EI)
    eb = _job([off, EB], load)
    assert eb["disp"]["2"][2] == pytest.approx(bending, rel=1e-12)
    assert eb["eleForce"]["1"][4] == pytest.approx(P * L, rel=1e-12)       # moment about Y at node 1
    timo = _job([off, TIMO], load)
    assert timo["disp"]["2"][2] == pytest.approx(P * Lf ** 3 / (3 * EI) + P * Lf / GA, rel=1e-12)
    assert timo["eleForce"]["1"][4] == pytest.approx(P * Lf, rel=1e-12)    # not P L: equilibrium violated
    eb_plain = _job([LINEAR, EB], load)["disp"]["2"][2]
    # geometric nonlinearity: compression of half the buckling load, PDelta transformation
    N = 0.5 * np.pi ** 2 * E * Iz / (4 * L ** 2)
    pd = ["geomTransf", "PDelta", 1, 0.0, -1.0, 0.0]
    load2 = [["load", 2, -N, 0.0, P, 0.0, 0.0, 0.0]]
    timo = _job([pd, TIMO], load2, second=True)["disp"]["2"][2]
    assert timo == pytest.approx(plain, rel=1e-9)                         # first-order answer
    eb = _job([pd, EB], load2, second=True)["disp"]["2"][2]
    assert eb / eb_plain > 1.5                                            # elasticBeamColumn amplifies


@needs_opensees
def test_opensees_single_element_pdelta_is_the_chord_rotation_stiffness():
    """One elasticBeamColumn with PDelta: tip deflection H / (3EI/L^3 - P/L) exactly, i.e. only the
    N/L 'string' term; the P-delta within the element is missing (hence the subdivision)."""
    N = 0.5 * np.pi ** 2 * E * Iz / (4 * L ** 2)
    H = 1000.0
    out = _job([["geomTransf", "PDelta", 1, 0.0, -1.0, 0.0], EB], [["load", 2, -N, 0.0, H, 0.0, 0.0, 0.0]], second=True)
    assert out["disp"]["2"][2] == pytest.approx(H / (3 * E * Iz / L ** 3 - N / L), rel=1e-9)


@needs_opensees
def test_opensees_reactions_and_chained_constraints_with_transformation_handler():
    """(a) nodeReaction is zero at a fixed node that is the retained node of a rigid link.
    (b) A constrained node used as the retained node of a second constraint is not handled reliably.
    The translator therefore computes reactions from joint equilibrium and never chains constraints
    (and verifies the constraints from the displacements it gets back)."""
    e = 300.0
    model = [["node", 10, e, 0.0, 0.0], LINEAR, ["element", "elasticBeamColumn", 1, 10, 2, A, E, G, J, Iy, Iz, 1],
             ["rigidLink", "beam", 1, 10]]
    P = 1000.0
    out = _job(model, [["load", 2, 0.0, 0.0, P, 0.0, 0.0, 0.0]], nodeReaction=[1], nodes=[1, 2, 10])
    assert out["disp"]["2"][2] == pytest.approx(P * (L - e) ** 3 / (3 * E * Iz), rel=1e-10)   # the answer is right
    assert np.abs(out["nodeReaction"]["1"]).max() == 0.0                                   # ... the reaction is not
    assert out["eleForce"]["1"][2] == pytest.approx(-P, rel=1e-10)                         # element forces are
    # (b) free node 3 --rigidLink--> 10 --equalDOF--> 11, element from 11: the second constraint is violated
    model = [["node", 3, 0.5 * e, 0.0, 0.0], ["node", 10, e, 0.0, 0.0], ["node", 11, e, 0.0, 0.0], LINEAR,
             ["element", "elasticBeamColumn", 1, 1, 3, A, E, G, J, Iy, Iz, 1],
             ["geomTransf", "Linear", 2, 0.0, -1.0, 0.0],
             ["element", "elasticBeamColumn", 2, 11, 2, A, E, G, J, Iy, Iz, 2],
             ["rigidLink", "beam", 3, 10], ["equalDOF", 10, 11, 1, 2, 3, 4, 5, 6]]
    job_nodes = [1, 2, 3, 10, 11]
    out = _job(model, [["load", 2, 0.0, 0.0, P, 0.0, 0.0, 0.0]], nodes=job_nodes, elements=[1, 2], beams=[1, 2])
    u10, u11 = np.array(out["disp"]["10"]), np.array(out["disp"]["11"])
    assert np.abs(u10 - u11).max() > 1e-6 * np.abs(u11).max()


@needs_pynite
def test_pynite_varying_load_on_a_member_split_at_intermediate_nodes():
    """A linearly varying line load on a member that PyNite splits into three or more pieces at nodes
    lying on it violates statics (uniform loads and two pieces are fine).  The translator therefore
    never lets PyNite split a member: sub-members are explicit and loads are clipped by the translator."""
    from Pynite import FEModel3D

    def reaction(n, w1, w2):
        m = FEModel3D()
        m.add_material("m", E, G, 0.3, 0.0)
        m.add_section("s", A, Iy, Iz, J)
        m.add_node("a", 0, 0, 0)
        m.add_node("b", L, 0, 0)
        for k in range(1, n):
            m.add_node(f"x{k}", L * k / n, 0, 0)
        m.add_member("M", "a", "b", "m", "s")
        m.def_support("a", True, True, True, True, True, True)
        m.add_member_dist_load("M", "Fy", w1, w2)
        m.analyze_linear()
        return m.nodes["a"].RxnFY["Combo 1"]

    total = -0.5 * (1.0 + 3.0) * L
    assert reaction(1, 1.0, 3.0) == pytest.approx(total, rel=1e-12)
    assert reaction(2, 1.0, 3.0) == pytest.approx(total, rel=1e-12)
    assert abs(reaction(4, 1.0, 3.0) - total) > 0.01 * abs(total)          # wrong total load
    assert reaction(4, 2.0, 2.0) == pytest.approx(-2.0 * L, rel=1e-12)    # uniform: fine


@needs_pynite
def test_pynite_pdelta_is_one_pass_with_first_order_axial_forces():
    """analyze_PDelta on a cantilever whose axial force is statically determinate converges to the
    beam-column solution (consistent geometric stiffness: one element is already within 0.5 %), and its
    reactions are post-processed with axial forces from the final displacements, so the horizontal
    reaction differs slightly from the applied load."""
    from Pynite import FEModel3D

    N = 0.5 * np.pi ** 2 * E * Iz / (4 * L ** 2)
    H = 1000.0
    m = FEModel3D()
    m.add_material("m", E, G, 0.3, 0.0)
    m.add_section("s", A, Iy, Iz, J)
    m.add_node("a", 0, 0, 0)
    m.add_node("b", L, 0, 0)
    m.add_member("M", "a", "b", "m", "s")
    m.def_support("a", True, True, True, True, True, True)
    m.add_node_load("b", "FX", -N)
    m.add_node_load("b", "FY", H)
    m.analyze_PDelta()
    u = L * np.sqrt(N / (E * Iz))
    exact = H * L ** 3 / (3 * E * Iz) * 3 * (np.tan(u) - u) / u ** 3      # derived in tests/test_oracle_beam.py
    assert m.nodes["b"].DY["Combo 1"] == pytest.approx(exact, rel=5e-3)
    rxn = m.nodes["a"].RxnFY["Combo 1"]
    assert rxn == pytest.approx(-H, rel=2e-3) and abs(rxn + H) > 1e-7 * H


@needs_opensees
def test_opensees_timoshenko_element_takes_no_release_flags_and_no_point_loads():
    """ElasticTimoshenkoBeam rejects ``-releasez`` (argument error) and has no ``-beamPoint`` element
    load ("load type unknown"); both end without results.  The translator uses neither."""
    base = {"solver": "opensees", "secondOrder": False, "tolerance": 1e-12, "maxIterations": 10, "system": "UmfPack",
            "nodes": [1, 2], "elements": [1], "beams": [1]}
    nodes = [["node", 1, 0.0, 0.0, 0.0], ["node", 2, L, 0.0, 0.0], ["fix", 1, 1, 1, 1, 1, 1, 1]]
    out = run_job(dict(base, model=nodes + [LINEAR, TIMO + ["-releasez", 2]],
                       pattern=[["load", 2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]]))
    assert out["ok"] is False
    out = run_job(dict(base, model=nodes + [LINEAR, TIMO],
                       pattern=[["eleLoad", "-ele", 1, "-type", "-beamPoint", 1000.0, 0.0, 0.3]]))
    assert out["ok"] is False or abs(out["disp"]["2"][2]) == 0.0
