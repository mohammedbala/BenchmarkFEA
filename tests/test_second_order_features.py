"""Second-order analysis: assembled members, connections, refusals, AISC 360-22 C2.2b / C2.3 options,
the record of the iteration, invariants, independence from the design block, CLI.

References are independent integrations of the beam-column equations (``second_order_reference``) or
statics written out in the test.
"""
import copy
import json
import math

import numpy as np
import pytest

from benchfea.cli import main
from benchfea.model import build_model
from benchfea.solve import ModelRefused, analyze

from second_order_reference import HELD_ROWS, beam_column_bvp, critical_k
from solve_reference import E, G, make

SEC = {"A": 5000.0, "Iy": 4.0e8, "Iz": 4.0e7, "J": 3.0e5, "Asy": 2500.0, "Asz": 3000.0}
EI, GA, EA = E * SEC["Iz"], G * SEC["Asy"], E * SEC["A"]
L = 4000.0
PE = math.pi ** 2 * EI / (2.0 * L) ** 2                  # load levels only
TOL = 2.0e-8


def column(loads, shear=False, supports=None, members=None, nodes=None, analysis=None, combos=None, **kw):
    doc = make(nodes or {"a": (0, 0, 0), "b": (0, 0, L)}, members or {"m": ("a", "b")},
               supports or {"a": "fixed"}, {"c": loads}, combos, sections=kw.pop("sections", {"S": {"type": "properties", **SEC}}),
               analysis={"shearDeformation": shear, **(analysis or {})}, **kw)
    if combos is None:
        doc["combinations"]["c"]["analysisType"] = "secondOrder"
    return doc


def solve(doc, cid="c", **kw):
    return analyze(build_model(doc), **kw).combinations[cid]


TIP = [{"type": "nodal", "node": "b", "Fz": -0.5 * PE, "Fx": 1000.0}]


# ------------------------------------------------------------------ assembled members and connections
def test_two_span_continuous_beam_column():
    """Two members, three supports, compression and a line load on one span: the joint rotation couples
    two members that both have a second-order stiffness."""
    L1, L2 = 3000.0, 4500.0
    P, q = 0.5 * math.pi ** 2 * EI / L2 ** 2, 4.0
    doc = column([{"type": "nodal", "node": "c", "Fx": -P},
                  {"type": "memberDistributed", "member": "m1", "direction": "Z", "w1": -q}], True,
                 nodes={"a": (0, 0, 0), "b": (L1, 0, 0), "c": (L1 + L2, 0, 0)},
                 members={"m1": ("a", "b"), "m2": ("b", "c")},
                 supports={"a": {"type": "pinned", "dofs": {"rx": "fixed"}}, "b": {"dofs": {"uy": "fixed", "uz": "fixed"}},
                           "c": {"dofs": {"uy": "fixed", "uz": "fixed"}}})
    c = solve(doc)

    def bc(s, e):
        return [s[0, 0], s[0, 2], e[1, 0], e[1, 2],                         # simple supports at the two ends
                e[0, 0], s[1, 0], s[1, 1] - e[0, 1], s[1, 2] - e[0, 2]]     # support at b: v = 0, slope and moment continuous
    ref = beam_column_bvp([L1, L2], EI, GA, -P, [-q, 0.0], bc)
    sv, sm = abs(ref(L1 / 2)[0]), abs(ref(L1)[2])
    for mid, x0 in (("m1", 0.0), ("m2", L1)):
        m = c.members[mid]
        for k, x in enumerate(m.x):
            r = ref(x0 + x, right=(mid == "m2" and k == 0))
            assert abs(m.displacements[k, 1] - r[0]) <= TOL * sv and abs(m.forces[k, 5] - r[2]) <= TOL * sm
    assert c.members["m2"].forces[5, 5] != 0.0 and c.invariants.passed


@pytest.mark.parametrize("shear", [False, True])
def test_rigid_end_offset_carrying_axial_load(shear):
    """The top 0.2 L of the column is a rigid offset with the loads on its tip.  For the flexible length l:
    V(l) = H and M(l) = (H + P theta(l)) a  -- the rigid arm turns with the member end, so the axial load
    acts with the lever arm a theta."""
    a = 0.2 * L
    ell = L - a
    P, H = 0.5 * PE, 1000.0
    c = solve(column(TIP, shear, members={"m": {"i": "a", "j": "b", "offsets": {"j": [0, 0, -a]}}}))
    ref = beam_column_bvp(ell, EI, GA if shear else None, -P, 0.0,
                          lambda s, e: [s[0, 0], s[0, 1], e[0, 3] - H, e[0, 2] - (H + P * e[0, 1]) * a])
    m = c.members["m"]
    assert m.x[-1] == pytest.approx(ell)
    for k, x in enumerate(m.x):
        assert m.displacements[k, 1] == pytest.approx(ref(x)[0], rel=TOL, abs=TOL * abs(ref(ell)[0]))
        assert m.forces[k, 5] == pytest.approx(ref(x)[2], rel=TOL)
    theta = ref(ell)[1]
    assert c.displacements[1, 0] == pytest.approx(ref(ell)[0] + a * theta, rel=TOL)          # node = tip of the arm
    assert c.displacements[1, 4] == pytest.approx(theta, rel=TOL)
    assert c.invariants.passed


def test_rigid_offset_equals_a_very_stiff_stub_member():
    a = 0.2 * L
    off = solve(column(TIP, members={"m": {"i": "a", "j": "b", "offsets": {"j": [0, 0, -a]}}}))
    stiff = {"type": "properties", **{k: v * (1.0e5 if k in ("A", "Iy", "Iz", "J") else 1.0) for k, v in SEC.items()}}
    doc = column([{"type": "nodal", "node": "t", "Fz": -0.5 * PE, "Fx": 1000.0}],
                 nodes={"a": (0, 0, 0), "b": (0, 0, L - a), "t": (0, 0, L)},
                 members={"m": ("a", "b"), "stub": {"i": "b", "j": "t", "section": "X"}},
                 sections={"S": {"type": "properties", **SEC}, "X": stiff})
    stub = solve(doc)
    assert stub.displacements[2, 0] == pytest.approx(off.displacements[1, 0], rel=2e-5)
    assert stub.reactions[0, 4] == pytest.approx(off.reactions[0, 4], rel=2e-5)


@pytest.mark.parametrize("kind", ["release", "support"])
def test_rotational_spring_connection(kind):
    """Base rotation held by a spring k: M(0) = k theta(0), posed as a spring-connected member end at a
    fixed node and as a spring support."""
    kt = 4.0 * EI / L
    P, H = 0.3 * PE, 1000.0
    loads = [{"type": "nodal", "node": "b", "Fz": -P, "Fx": H}]
    if kind == "release":
        doc = column(loads, True, members={"m": {"i": "a", "j": "b", "releases": {"i": {"rz": {"spring": kt}}}}})
    else:
        doc = column(loads, True, supports={"a": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed", "ry": {"spring": kt}}}})
    c = solve(doc)
    ref = beam_column_bvp(L, EI, GA, -P, 0.0, lambda s, e: [s[0, 0], s[0, 2] - kt * s[0, 1], e[0, 2], e[0, 3] - H])
    assert c.displacements[1, 0] == pytest.approx(ref(L)[0], rel=TOL)
    assert c.members["m"].forces[0, 5] == pytest.approx(ref(0.0)[2], rel=TOL)
    fixed = solve(column(loads, True))
    assert c.displacements[1, 0] > 1.5 * fixed.displacements[1, 0]


def test_truss_member_with_a_transverse_load_is_a_pin_ended_beam_column():
    """A compressed truss-type member carrying a load between its ends: its own bending moment is
    amplified (P-delta inside an axial-only member)."""
    P, q = 0.5 * math.pi ** 2 * EI / L ** 2, 2.0
    doc = column([{"type": "nodal", "node": "b", "Fx": -P},
                  {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -q}],
                 nodes={"a": (0, 0, 0), "b": (L, 0, 0)}, members={"m": {"i": "a", "j": "b", "type": "truss"}},
                 supports={"a": "pinned", "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}})
    c = solve(doc)
    ref = beam_column_bvp(L, EI, None, -P, -q, lambda s, e: [s[0, 0], s[0, 2], e[0, 0], e[0, 2]])
    m = c.members["m"]
    for k, x in enumerate(m.x):
        assert m.forces[k, 5] == pytest.approx(ref(x)[2], rel=TOL, abs=TOL * abs(ref(L / 2)[2]))
    assert m.forces[5, 5] > 1.9 * q * L ** 2 / 8 and c.invariants.passed


# ----------------------------------------------------------------------------------- refusals
def weak_axis_load():
    """Critical load of the cantilever in its weak plane (local xz... here the stiff plane is xz, so the
    weak plane is xy with EI = E Iz), from the characteristic determinant."""
    k = critical_k(lambda k: [HELD_ROWS["held"](k, 0.0), HELD_ROWS["no rotation"](k, 0.0),
                              HELD_ROWS["no moment"](k, L), HELD_ROWS["no shear"](k, L)], 6.0 / L)
    return k * k * EI


def test_load_beyond_the_elastic_critical_load_is_refused_with_the_buckling_factor():
    pcr = weak_axis_load()
    loads = {"ok": [{"type": "nodal", "node": "b", "Fz": -0.5 * pcr, "Fx": 1000.0}],
             "over": [{"type": "nodal", "node": "b", "Fz": -1.25 * pcr, "Fx": 1000.0}]}
    doc = make({"a": (0, 0, 0), "b": (0, 0, L)}, {"m": ("a", "b")}, {"a": "fixed"}, loads,
               {k: {"factors": {k: 1.0}, "analysisType": "secondOrder"} for k in loads},
               sections={"S": {"type": "properties", **SEC}}, analysis={"shearDeformation": False})
    model = build_model(doc)
    with pytest.raises(ModelRefused) as e:
        analyze(model)
    p = e.value.problems[0]
    assert p.code == "stability.beyond-critical-load" and p.path == "/combinations/over"
    assert "exceeds the elastic critical load" in p.message and "buckling factor ~ 0.8" in p.message
    assert e.value.refused["combination 'over'"].info["bucklingFactor"] == pytest.approx(1.0 / 1.25, rel=1e-6)
    assert set(e.value.results.combinations) == {"ok"}                   # the other combination is delivered
    partial = analyze(model, strict=False)
    assert set(partial.combinations) == {"ok"} and any("NO RESULTS for combination 'over'" in w for w in partial.warnings)


def test_just_below_the_critical_load_the_analysis_still_converges():
    pcr = weak_axis_load()
    c = solve(column([{"type": "nodal", "node": "b", "Fz": -0.98 * pcr, "Fx": 1000.0}]))
    first = 1000.0 * L ** 3 / (3 * EI)
    assert c.converged and c.displacements[1, 0] / first == pytest.approx(1.0 / (1.0 - 0.98), rel=0.03)


def test_member_beyond_its_own_critical_load_between_released_ends_is_refused():
    """A truss-type member is pin-ended: beyond its own Euler load it has buckled even though the joints
    around it are held."""
    euler = critical_k(lambda k: [HELD_ROWS["held"](k, 0.0), HELD_ROWS["no moment"](k, 0.0),
                                  HELD_ROWS["held"](k, L), HELD_ROWS["no moment"](k, L)], 6.0 / L) ** 2 * E * SEC["Iz"]
    for factor, ok in ((0.9, True), (1.1, False)):
        doc = column([{"type": "nodal", "node": "b", "Fz": -factor * euler}],
                     members={"m": {"i": "a", "j": "b", "type": "truss"}},
                     supports={"a": "pinned", "b": {"dofs": {"ux": "fixed", "uy": "fixed"}}})
        if ok:
            assert solve(doc).converged
        else:
            with pytest.raises(ModelRefused) as e:
                analyze(build_model(doc))
            assert e.value.problems[0].code == "stability.beyond-critical-load"
            assert e.value.refused["combination 'c'"].info["bucklingFactor"] == pytest.approx(1.0 / factor, rel=1e-6)


def test_released_transverse_translation_is_refused_for_second_order_and_buckling():
    doc = column(TIP, members={"m": {"i": "a", "j": "b", "releases": {"j": {"uy": {"spring": 500.0}}}}},
                 supports={"a": "fixed", "b": {"dofs": {"ux": "fixed"}}})
    with pytest.raises(ModelRefused) as e:
        analyze(build_model(doc))
    assert e.value.problems[0].code == "analysis.release-second-order" and e.value.problems[0].path == "/members/m/releases"
    assert "j.uy" in e.value.problems[0].message
    doc["combinations"]["c"].pop("analysisType")
    doc["analysis"]["buckling"] = {"combinations": ["c"]}
    with pytest.raises(ModelRefused) as e:
        analyze(build_model(doc))
    assert e.value.problems[0].code == "analysis.release-second-order" and "buckling analysis" in e.value.problems[0].message
    doc["analysis"].pop("buckling")
    assert analyze(build_model(doc)).combinations["c"].invariants.passed      # first order: fine


# -------------------------------------------------------------- record of the iteration, convergence
def portal(analysis=None, factor=1.0):
    """Fixed-base portal, gravity on the columns and a lateral load: the column axial forces depend on
    the displacements (overturning), so the iteration needs more than two passes."""
    h, b = 4000.0, 6000.0
    keep = {"dofs": {"uy": "fixed", "rx": "fixed", "rz": "fixed"}}
    doc = make({"a": (0, 0, 0), "b": (0, 0, h), "c": (b, 0, h), "d": (b, 0, 0)},
               {"c1": ("a", "b"), "g": ("b", "c"), "c2": ("d", "c")}, {"a": "fixed", "d": "fixed", "b": keep, "c": keep},
               {"c": [{"type": "nodal", "node": "b", "Fz": -1.5e6 * factor, "Fx": 5.0e4 * factor},
                      {"type": "nodal", "node": "c", "Fz": -1.5e6 * factor},
                      {"type": "memberDistributed", "member": "g", "direction": "Z", "w1": -20.0 * factor}]},
               sections={"S": {"type": "properties", **SEC}}, analysis={"shearDeformation": True, **(analysis or {})})
    doc["combinations"]["c"]["analysisType"] = "secondOrder"
    return doc


def test_iterations_history_and_subdivisions_are_recorded():
    c = solve(portal())
    so = c.solver["secondOrder"]
    assert c.converged and c.iterations == so["iterations"] == len(so["history"]) >= 3
    changes = [h["displacementChange"] for h in so["history"]]
    assert changes == sorted(changes, reverse=True) and changes[0] > 1e-2          # contracting from a real change
    assert changes[-1] <= so["tolerance"] == 1e-6 and so["history"][-1]["axialForceChange"] <= 1e-6
    assert so["history"][0]["iteration"] == 1 and set(so["subdivisions"]) == {"c1", "g", "c2"}
    assert so["maxIterations"] == 50 and so["largestPieceKh"] <= 3.0
    tight = solve(portal({"convergenceTolerance": 1e-11}))
    assert tight.iterations > c.iterations
    assert tight.displacements[1, 0] == pytest.approx(c.displacements[1, 0], rel=1e-6)   # the stated tolerance holds
    first = solve({**portal(), "combinations": {"c": {"factors": {"c": 1.0}}}})
    assert c.displacements[1, 0] > 1.1 * first.displacements[1, 0]
    assert c.members["c1"].forces[0, 0] != pytest.approx(first.members["c1"].forces[0, 0], rel=1e-4)   # N changed


def test_iteration_cap_gives_a_flagged_not_converged_result():
    model = build_model(portal({"maxIterations": 2}))
    with pytest.raises(ModelRefused) as e:
        analyze(model)
    assert e.value.problems[0].code == "solve.not-converged" and e.value.problems[0].path == "/combinations/c"
    c = analyze(model, strict=False).combinations["c"]
    assert not c.converged and c.iterations == 2 and any("did not converge in 2 iteration" in w for w in c.warnings)
    assert c.solver["secondOrder"]["converged"] is False


def test_invariants_of_a_second_order_solution_are_evaluated_with_the_second_order_terms():
    c = solve(portal())
    names = [ch.name for ch in c.invariants.checks]
    assert names == ["stiffness symmetry", "rigid-body modes", "positive definiteness", "linear solve residual",
                     "global force equilibrium", "global moment equilibrium", "joint equilibrium",
                     "member equilibrium", "member compatibility", "member equilibrium on the deformed axis",
                     "work = strain energy"]
    assert c.invariants.passed
    # reactions balance the loads as forces; as moments about the origin on the UNDEFORMED geometry they do
    # not (the difference is the P-Delta couple the check accounts for)
    assert c.reactions[:, :3].sum(axis=0) == pytest.approx([-5.0e4, 0.0, 3.0e6 + 20.0 * 6000.0], rel=1e-9, abs=1e-6)
    xyz = np.array([[0, 0, 0], [0, 0, 4000.0], [6000.0, 0, 4000.0], [6000.0, 0, 0]])
    m_react = (np.cross(xyz, c.reactions[:, :3]) + c.reactions[:, 3:]).sum(axis=0)[1]
    m_loads = 4000.0 * 5.0e4 + 6000.0 * 1.5e6 + 20.0 * 6000.0 * 3000.0          # about +Y (r x F)
    # ... they miss exactly the couples of the member axial forces N acting through the transverse
    # displacement of the member ends: columns N x (sway of the top), beam N x (difference of the vertical
    # end displacements); about +Y, with N < 0 in compression
    n1, n2, ng = (c.members[k].forces[0, 0] for k in ("c1", "c2", "g"))
    couple = n1 * c.displacements[1, 0] + n2 * c.displacements[2, 0] - ng * (c.displacements[2, 2] - c.displacements[1, 2])
    assert m_react + m_loads == pytest.approx(couple, rel=1e-8)
    assert abs(couple) > 1e-3 * m_loads and couple < 0.0                        # and that is not a small number


# ------------------------------------------------------------------- AISC 360-22 C2.1(d): alpha
def test_alpha_multiplies_the_loads_and_divides_the_results():
    """C2.1(d): ASD combinations are analysed at 1.6 x the loads and the results divided by 1.6."""
    asd = solve(column(TIP, True, analysis={"alpha": 1.6}))
    loads = [{"type": "nodal", "node": "b", "Fz": -0.5 * PE * 1.6, "Fx": 1600.0}]
    scaled = solve(column(loads, True))
    plain = solve(column(TIP, True))
    assert asd.alpha == 1.6 and plain.alpha == 1.0
    assert asd.displacements == pytest.approx(scaled.displacements / 1.6, rel=1e-12, abs=1e-18)
    assert asd.members["m"].forces == pytest.approx(scaled.members["m"].forces / 1.6, rel=1e-12, abs=1e-9)
    assert asd.reactions[0, 2] == pytest.approx(0.5 * PE, rel=1e-12)             # the required strength level
    assert asd.displacements[1, 0] > 1.5 * plain.displacements[1, 0]             # more amplification at 1.6 P
    first_order = solve({**column(TIP, True, analysis={"alpha": 1.6}), "combinations": {"c": {"factors": {"c": 1.0}}}})
    assert first_order.alpha == 1.0                                              # alpha is for second order only


# ------------------------------------------------------------ AISC 360-22 C2.3: stiffness reduction
LS = 2000.0                                               # a stockier column, so that P can exceed 0.5 Py
FY = 345.0
PY = FY * SEC["A"]


def reduced(P, tau, shear=True, members=None, **combo):
    doc = column([{"type": "nodal", "node": "b", "Fz": -P, "Fx": 1000.0}], shear,
                 nodes={"a": (0, 0, 0), "b": (0, 0, LS)}, materials={"M": {"E": E, "G": G, "Fy": FY}},
                 analysis={"stiffnessReduction": {"tauB": tau, **({"members": members} if members is not None else {})}})
    doc["combinations"]["c"].update({"stiffnessReduction": True, **combo})
    return doc


@pytest.mark.parametrize("tau,ratio,expected", [("none", 0.3, 1.0), (0.7, 0.3, 0.7), ("computed", 0.3, 1.0),
                                                ("computed", 0.5, 1.0), ("computed", 0.7, 0.84), ("computed", 0.8, 0.64)])
def test_stiffness_reduction_c2_3(tau, ratio, expected):
    """C2.3(a): 0.80 on all stiffnesses; C2.3(b): tau_b on the flexural stiffness, Eq. C2-2a (alpha Pr / Pns
    <= 0.5: 1.0) / Eq. C2-2b (4 x (1 - x)), Pns = Fy A.  Reference: the cantilever with 0.8 tau_b EI and
    0.8 G As (tau_b does not apply to the shear, axial or torsional stiffness)."""
    P = ratio * PY
    c = solve(reduced(P, tau))
    assert c.tau_b == pytest.approx({"m": expected}, rel=1e-12)
    ref = beam_column_bvp(LS, 0.8 * expected * EI, 0.8 * GA, -P, 0.0,
                          lambda s, e: [s[0, 0], s[0, 1], e[0, 2], e[0, 3] - 1000.0])
    assert c.displacements[1, 0] == pytest.approx(ref(LS)[0], rel=TOL)
    assert c.members["m"].forces[0, 5] == pytest.approx(ref(0.0)[2], rel=TOL)
    assert c.displacements[1, 2] == pytest.approx(-P * LS / (0.8 * EA), rel=1e-12)       # 0.8 EA, no tau_b
    rec = c.solver["stiffnessReduction"]
    assert rec["factor"] == 0.8 and rec["tauB"] == tau and rec["members"] == ["m"] and "C2.3" in rec["clause"]
    assert c.invariants.passed and c.converged


def test_tau_b_uses_the_compression_of_the_analysis_at_alpha_times_the_loads():
    """alpha Pr / Pns with alpha = 1.6: Pr = 0.4 Py gives 0.64 -> tau_b = 4 x 0.64 x 0.36."""
    doc = reduced(0.4 * PY, "computed")
    doc["analysis"]["alpha"] = 1.6
    c = solve(doc)
    assert c.tau_b["m"] == pytest.approx(4 * 0.64 * 0.36, rel=1e-12) and c.alpha == 1.6


def test_tau_b_is_one_in_tension_and_the_reduction_can_be_limited_to_listed_members():
    c = solve(reduced(-0.9 * PY, "computed"))
    assert c.tau_b == {"m": 1.0}
    nodes = {"a": (0, 0, 0), "b": (0, 0, LS / 2), "c": (0, 0, LS)}
    doc = column([{"type": "nodal", "node": "c", "Fz": -0.7 * PY, "Fx": 1000.0}], True, nodes=nodes,
                 members={"lower": ("a", "b"), "upper": ("b", "c")}, materials={"M": {"E": E, "G": G, "Fy": FY}},
                 analysis={"stiffnessReduction": {"members": ["lower"]}})
    doc["combinations"]["c"]["stiffnessReduction"] = True
    c = solve(doc)
    assert c.tau_b == pytest.approx({"lower": 0.84}) and c.solver["stiffnessReduction"]["members"] == ["lower"]
    # shortening: lower part with 0.8 EA, upper part unreduced
    assert c.displacements[2, 2] == pytest.approx(-0.7 * PY * (LS / 2) / EA * (1 / 0.8 + 1.0), rel=1e-12)


def test_compression_at_or_above_the_cross_section_strength_is_refused_with_computed_tau_b():
    with pytest.raises(ModelRefused) as e:
        analyze(build_model(reduced(1.02 * PY, "computed")))
    p = e.value.problems[0]
    assert p.code == "analysis.tau-b-zero" and "Pns = Fy A" in p.message and "C2-2b" in p.message


# ------------------------------------------------------------- AISC 360-22 C2.2b: notional loads
def notional_frame(second=True, coefficient=None, direction=(1, 0, 0), alpha=1.0):
    h, b = 4000.0, 6000.0
    keep = {"dofs": {"uy": "fixed", "rx": "fixed", "rz": "fixed"}}
    loads = [{"type": "nodal", "node": "b", "Fz": -3.0e5}, {"type": "nodal", "node": "c", "Fz": -1.0e5, "Fx": 7000.0},
             {"type": "memberDistributed", "member": "g", "direction": "Z", "w1": -10.0},
             {"type": "memberPoint", "member": "g", "direction": "Z", "value": -8.0e4, "position": 0.25},
             {"type": "memberPoint", "member": "c1", "direction": "X", "value": 9000.0, "position": 0.5}]
    doc = make({"a": (0, 0, 0), "b": (0, 0, h), "c": (b, 0, h), "d": (b, 0, 0)},
               {"c1": ("a", "b"), "g": ("b", "c"), "c2": ("d", "c")}, {"a": "fixed", "d": "fixed", "b": keep, "c": keep},
               {"c": loads}, sections={"S": {"type": "properties", **SEC}}, analysis={"alpha": alpha})
    nl = {"direction": list(direction)}
    if coefficient is not None:
        nl["coefficient"] = coefficient
    doc["combinations"]["c"].update({"notionalLoad": nl, **({"analysisType": "secondOrder"} if second else {})})
    return doc


@pytest.mark.parametrize("second", [True, False])
def test_notional_loads_c2_2b(second):
    """Eq. C2-1: Ni = 0.002 alpha Yi at every node, Yi = gravity load at the node with member loads lumped
    statically to the member ends (lever rule): b: 3e5 + 10 x 6000 / 2 + 0.75 x 8e4; c: 1e5 + 3e4 + 0.25 x 8e4."""
    doc = notional_frame(second)
    c = solve(doc)
    rec = c.solver["notionalLoads"]
    Yb, Yc = 3.0e5 + 3.0e4 + 6.0e4, 1.0e5 + 3.0e4 + 2.0e4
    assert rec["coefficient"] == 0.002 and "C2-1" in rec["clause"] and rec["unit"] == "N"
    assert set(rec["nodes"]) == {"b", "c"}                                       # horizontal loads carry no Y
    assert rec["nodes"]["b"] == pytest.approx([0.002 * Yb, 0.0, 0.0]) and rec["nodes"]["c"] == pytest.approx([0.002 * Yc, 0.0, 0.0])
    # they are loads: the base shear contains them
    assert c.reactions[:, 0].sum() == pytest.approx(-(7000.0 + 9000.0 + 0.002 * (Yb + Yc)), rel=1e-9)
    plain = copy.deepcopy(doc)
    plain["combinations"]["c"].pop("notionalLoad")
    explicit = copy.deepcopy(plain)
    explicit["loadCases"]["c"]["loads"] += [{"type": "nodal", "node": "b", "Fx": 0.002 * Yb},
                                            {"type": "nodal", "node": "c", "Fx": 0.002 * Yc}]
    assert c.displacements == pytest.approx(solve(explicit).displacements, rel=1e-10, abs=1e-14)
    assert c.displacements[1, 0] > solve(plain).displacements[1, 0]
    assert c.invariants.passed


def test_notional_load_coefficient_direction_and_alpha():
    c = solve(notional_frame(coefficient=0.003, direction=(0, -2, 0)))           # direction is normalised
    assert c.solver["notionalLoads"]["nodes"]["b"] == pytest.approx([0.0, -0.003 * 3.9e5, 0.0])
    assert c.solver["notionalLoads"]["direction"] == [0.0, -1.0, 0.0]
    c = solve(notional_frame(alpha=1.6))                                         # Ni = 0.002 alpha Yi at the analysis level
    assert c.solver["notionalLoads"]["nodes"]["b"] == pytest.approx([0.002 * 1.6 * 3.9e5, 0.0, 0.0])
    assert c.reactions[:, 0].sum() == pytest.approx(-(16000.0 + 0.002 * (3.9e5 + 1.5e5)), rel=1e-9)   # results / alpha


def test_direct_analysis_settings_together():
    """Second order + notional loads + 0.8 tau_b: the usual set of AISC 360-22 Chapter C, on one frame."""
    doc = notional_frame()
    doc["materials"]["M"]["Fy"] = FY
    doc["combinations"]["c"]["stiffnessReduction"] = True
    c = solve(doc)
    assert c.converged and c.invariants.passed and set(c.tau_b) == {"c1", "g", "c2"} and set(c.tau_b.values()) == {1.0}
    stiff = solve(notional_frame())
    assert c.displacements[1, 0] > 1.25 * stiff.displacements[1, 0]              # 1 / 0.8 and more amplification


# -------------------------------------------------------------- design data never changes the analysis
HEAVY_DESIGN = {
    "defaults": {"method": "ASD", "Kx": 2.1, "Ky": 0.5, "Kz": 3.0, "Cb": 1.7, "Fy": 250.0, "E": 150000.0},
    "members": {"c1": {"Lx": 9000.0, "Ly": 1000.0, "Lb": 500.0, "Kx": 3.3},
                "c2": {"bracing": {"minor": [0.25, 0.5], "torsion": [0.5]}, "Ky": 0.65}},
}


def test_unbraced_lengths_and_k_factors_do_not_change_second_order_or_buckling_results_bitwise():
    doc = portal({"buckling": {"combinations": ["c"], "modes": 3}}, factor=0.5)
    doc["materials"]["M"]["Fy"] = FY
    doc["materials"]["M"]["Fu"] = 450.0
    doc["combinations"]["c"]["stiffnessReduction"] = True
    doc["combinations"]["c"]["notionalLoad"] = {"direction": [1, 0, 0]}
    with_design = copy.deepcopy(doc)
    with_design["design"] = HEAVY_DESIGN
    a, b = analyze(build_model(doc)), analyze(build_model(with_design))
    ca, cb = a.combinations["c"], b.combinations["c"]
    assert np.array_equal(ca.displacements, cb.displacements) and np.array_equal(ca.reactions, cb.reactions)
    for mid in ca.members:
        assert np.array_equal(ca.members[mid].forces, cb.members[mid].forces)
        assert np.array_equal(ca.members[mid].displacements, cb.members[mid].displacements)
    assert ca.iterations == cb.iterations and ca.tau_b == cb.tau_b and ca.solver == cb.solver
    assert np.array_equal(a.buckling["c"].factors, b.buckling["c"].factors)
    assert np.array_equal(a.buckling["c"].modes, b.buckling["c"].modes)
    assert len(a.buckling["c"].factors) == 3


# ------------------------------------------------------------------------------------------- CLI
def test_cli_run_reports_iterations_and_buckling_factors(tmp_path, capsys):
    doc = portal({"buckling": {"combinations": ["c"], "modes": 2}})
    path = tmp_path / "portal.json"
    path.write_text(json.dumps(doc))
    assert main(["run", str(path)]) == 0
    out = capsys.readouterr().out
    assert "c [secondOrder]" in out and "second-order iteration converged in" in out and "iteration(s)" in out
    assert "buckling of combination c: factors" in out
    assert main(["run", str(path), "--json"]) == 0
    summary = json.loads(capsys.readouterr().out)
    so = summary["combinations"]["c"]["secondOrder"]
    assert so["converged"] is True and so["iterations"] >= 3 and so["lastDisplacementChange"] <= 1e-6
    res = analyze(build_model(doc))
    assert summary["buckling"]["c"]["factors"] == pytest.approx(list(res.buckling["c"].factors), rel=1e-12)
    written = json.loads((tmp_path / "portal.results.json").read_text())
    assert written["combinations"]["c"]["solver"]["secondOrder"]["iterations"] == so["iterations"]
    assert written["buckling"]["c"]["factors"] == pytest.approx(summary["buckling"]["c"]["factors"])


def test_cli_run_refuses_a_combination_beyond_its_critical_load(tmp_path, capsys):
    doc = column([{"type": "nodal", "node": "b", "Fz": -1.25 * weak_axis_load(), "Fx": 1000.0}])
    path = tmp_path / "over.json"
    path.write_text(json.dumps(doc))
    assert main(["run", str(path)]) == 1
    err = capsys.readouterr().err
    assert "stability.beyond-critical-load" in err and "buckling factor ~ 0.8" in err and "REFUSED" in err
    assert not (tmp_path / "over.results.json").exists()
