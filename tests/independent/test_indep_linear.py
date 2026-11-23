"""Independent review, first order: 27 irregular models against three references that share no code
with ``benchfea.solve`` -- the reviewer's tableau solver (every model, every output incl. stations),
OpenSeesPy and PyNite (the models each can represent exactly).

Tolerances are set from the fidelity of the references, not from the measured differences: all three
solve the same linear problem exactly, so 1e-9 (relative to the largest value of each quantity kind)
is demanded; the ill-conditioned model (20 mm stubs next to 30 m spans) gets 1e-8.
Measured: <= 4e-12 everywhere, <= 2e-10 on the ill-conditioned model.
"""
from functools import lru_cache

import numpy as np
import pytest

import oracles
from harness import compare, compare_foreign, run_benchfea
from models import LINEAR

TOL = {"L19_very_short_and_very_long_members": 1.0e-8}
QUANTITIES = ("translations", "rotations", "reaction forces", "reaction moments", "end forces", "end moments",
              "station forces", "station moments", "station displacements", "constraint forces", "constraint moments")


@lru_cache(maxsize=None)
def solved(name):
    raw, cid = LINEAR[name]()
    _, res, doc = run_benchfea(raw)
    return raw, cid, res, doc


@pytest.mark.parametrize("name", sorted(LINEAR))
def test_first_order_against_the_tableau_reference(name):
    raw, cid, res, doc = solved(name)
    out = compare(raw, cid, doc=doc)
    assert out["_tableau"].residual < 1e-10, "the reference itself did not solve its equations"
    for q in QUANTITIES:
        assert out[q] <= TOL.get(name, 1.0e-9), f"{name}: {q} differs by {out[q]:.2e}"
    assert doc["combinations"][cid]["invariants"]["passed"]


# models each foreign solver can represent EXACTLY (the translators decline the others: trapezoids and
# member point moments in OpenSees, shear deformation / offsets / constraints / partial releases in PyNite, ...)
COVERED = {"opensees": ["L02", "L03", "L04", "L05", "L07", "L08", "L09", "L10", "L12", "L14", "L18", "L19", "L20",
                        "L21", "L23", "L24", "L26", "L27"],
           "pynite": ["L02", "L04", "L11", "L12", "L13", "L16", "L17", "L18", "L19"]}
_BY_NUMBER = {name[:3]: name for name in LINEAR}
_RUN = {"opensees": oracles.run_opensees, "pynite": oracles.run_pynite}


@pytest.mark.parametrize("solver, number", [(s, n) for s in COVERED for n in COVERED[s]])
def test_first_order_against_opensees_and_pynite(solver, number):
    pytest.importorskip({"opensees": "openseespy.opensees", "pynite": "Pynite"}[solver])
    name = _BY_NUMBER[number]
    raw, cid, res, doc = solved(name)
    f = _RUN[solver](raw, cid)
    assert f is not None, f"{solver} translator declined {name}"
    for q, v in compare_foreign(raw, doc["combinations"][cid], f).items():
        assert v <= TOL.get(name, 1.0e-9), f"{name} vs {solver}: {q} differs by {v:.2e}"


@pytest.mark.parametrize("solver", ["opensees", "pynite"])
def test_models_outside_the_exact_scope_of_a_foreign_solver_are_declined_not_approximated(solver):
    pytest.importorskip({"opensees": "openseespy.opensees", "pynite": "Pynite"}[solver])
    for number, name in _BY_NUMBER.items():
        if number not in COVERED[solver]:
            assert _RUN[solver](*LINEAR[name]()) is None, name


def test_zero_loads_give_exact_zeros():
    raw, cid, res, doc = solved("L20_zero_loads")
    c = res.combinations[cid]
    assert not c.displacements.any() and not c.reactions.any()
    assert all(not m.end_forces.any() and not m.forces.any() for m in c.members.values())
    assert c.invariants.passed


def test_results_do_not_depend_on_the_unit_system_of_the_file():
    """The same structure written in kN-m and in kip-in gives the same physical results (the reference
    works in the file's units, the engine in N-mm)."""
    raw, cid, res, doc = solved("L05_offsets_off_axis_skew")
    kn, m = 4.4482216152605, 0.0254                      # 1 kip = 4.448... kN, 1 in = 0.0254 m
    other, _ = LINEAR["L05_offsets_off_axis_skew"]()
    other["units"] = {"force": "kip", "length": "in"}
    for n in other["nodes"].values():
        for k in "xyz":
            n[k] /= m
    for mat in other["materials"].values():
        mat["E"] *= m * m / kn
        mat["G"] *= m * m / kn
        mat["unitWeight"] *= m ** 3 / kn
    for s in other["sections"].values():
        s["A"] /= m ** 2
        for k in ("Iy", "Iz", "J"):
            s[k] /= m ** 4
    for mem in other["members"].values():
        if "offsets" in mem:
            for e in ("i", "j"):
                if e in mem["offsets"]:
                    mem["offsets"][e] = [v / m for v in mem["offsets"][e]]
    for ld in other["loadCases"]["W"]["loads"]:
        if ld["type"] == "memberDistributed":
            ld["w1"] *= m / kn
        elif ld["type"] == "memberPoint":
            ld["value"] /= kn
        elif ld["type"] == "nodal":
            for k in ("Fx", "Fy", "Fz"):
                if k in ld:
                    ld[k] /= kn
            for k in ("Mx", "My", "Mz"):
                if k in ld:
                    ld[k] /= kn * m
    _, res2, _ = run_benchfea(other)
    a, b = res.combinations[cid], res2.combinations[cid]
    for x, y in ((a.displacements, b.displacements), (a.reactions, b.reactions)):
        assert np.abs(x - y).max() <= 1e-9 * np.abs(x).max()
    for mid in a.members:
        x, y = a.members[mid].forces, b.members[mid].forces
        assert np.abs(x - y).max() <= 1e-9 * max(np.abs(x).max(), 1.0)
