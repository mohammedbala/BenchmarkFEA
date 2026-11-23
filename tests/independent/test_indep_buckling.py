"""Independent review, linear buckling: the engine's factors against the roots of the determinant of
the reviewer's tableau with N = factor x (first-order axial forces) -- the transcendental problem,
found by a sign scan and bisection; no geometric stiffness matrix, no subdivision.

For a constant axial force per member both are exact for the same theory: 1e-6 is demanded
(the engine states convergence to 1e-7); measured <= 7e-11.
"""
import pytest

from harness import run_benchfea
from models import ONE_SIDED, SECOND
from tableau import Tableau


def _buckling(raw, cid, modes):
    raw["combinations"][cid]["analysisType"] = "firstOrder"
    raw.setdefault("analysis", {})["buckling"] = {"combinations": [cid], "modes": modes}
    raw["analysis"].pop("convergenceTolerance", None)
    _, res, doc = run_benchfea(raw)
    return doc["buckling"][cid]


@pytest.mark.parametrize("name, modes", [
    ("S01_cantilever_beam_column_both_planes", 3),                        # both planes of a rolled column
    ("S02_cantilever_beam_column_with_shear_deformation_and_tension", 2),  # shear-flexible (Engesser), a tie in tension
    ("S03_portal_with_leaning_column", 2),                                 # sway frame + leaning column, releases
    ("S05_offsets_along_the_axis_carrying_compression", 2),                # rigid arms carrying the load
    ("S08_diaphragm_storey_with_leaning_columns", 2)])                     # torsional sway of a diaphragm storey
def test_buckling_factors_are_the_roots_of_the_reference_determinant(name, modes):
    raw, cid = SECOND[name]()
    b = _buckling(raw, cid, modes)
    factors = b["factors"]
    assert len(factors) == modes
    roots = Tableau(raw).buckling(cid, n_modes=modes, lam_max=factors[-1] * 1.02, steps=30 * modes)
    assert len(roots) == modes, (factors, roots)                 # no root below the engine's that it missed
    for f, r in zip(factors, roots):
        assert abs(f - r) <= 1.0e-6 * r, (factors, roots)


def test_first_buckling_factor_of_the_cantilever_is_the_weak_axis_euler_load():
    raw, cid = SECOND["S01_cantilever_beam_column_both_planes"]()
    import math
    b = _buckling(raw, cid, 1)
    euler = math.pi ** 2 * 2.0e8 * 4.0e-5 / (2 * 4.0) ** 2             # kN, effective length 2 L derived from v'' + k^2 v = 0
    assert abs(b["factors"][0] - euler / 450.0) <= 1e-6 * euler / 450.0


def test_buckling_with_partial_releases_truss_members_and_reversed_loads():
    """Partial rotational releases, pin-ended truss members that buckle on their own, and the factor
    under reversed loads (``negativeFactors``)."""
    raw, cid = SECOND["S06_truss_members_with_transverse_load_and_partial_releases"]()
    raw["loadCases"]["U"]["loads"] = [ld for ld in raw["loadCases"]["U"]["loads"] if ld["type"] == "nodal"]
    b = _buckling(raw, cid, 2)
    t = Tableau(raw)
    roots = t.buckling(cid, n_modes=2, lam_max=b["factors"][-1] * 1.02, steps=60)
    assert len(roots) == 2
    for f, r in zip(b["factors"], roots):
        assert abs(f - r) <= 1.0e-6 * r
    neg = t.buckling(cid, n_modes=1, lam_max=abs(b["negativeFactors"][0]) * 1.02, steps=40, negative=True)
    assert abs(b["negativeFactors"][0] - neg[0]) <= 1.0e-6 * abs(neg[0])


def test_buckling_with_an_axial_force_that_varies_along_loaded_members():
    """Same frame with its member loads (axial components along the truss members): the reference's
    piecewise constant N (8 steps) is O(h^2) accurate; 1e-5 is demanded (measured 1.2e-7)."""
    raw, cid = SECOND["S06_truss_members_with_transverse_load_and_partial_releases"]()
    b = _buckling(raw, cid, 1)
    t = Tableau(raw)
    t.nsub = 8
    (root,) = t.buckling(cid, n_modes=1, lam_max=b["factors"][0] * 1.02, steps=30)
    assert abs(b["factors"][0] - root) <= 1.0e-5 * root


def test_buckling_keeps_the_slack_set_of_the_reference_state():
    raw, cid = ONE_SIDED["T01_two_storey_x_braced_frame"]()
    b = _buckling(raw, cid, 1)
    assert sorted(b["info"]["slack"]) == ["x1b", "x2a"]
    (root,) = Tableau(raw).buckling(cid, n_modes=1, lam_max=b["factors"][0] * 1.02, steps=20)
    assert abs(b["factors"][0] - root) <= 1.0e-6 * root
