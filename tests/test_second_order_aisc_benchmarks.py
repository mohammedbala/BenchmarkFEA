"""AISC 360-22 Commentary C2.1 benchmark problems for second-order analysis: Case 1 (Fig. C-C2.2,
simply supported beam-column, P-delta only) and Case 2 (Fig. C-C2.3, fixed-base cantilever, P-Delta and
P-delta).  Every tabulated moment and deflection, with and without shear deformation, US and SI tables,
against the Commentary's acceptance criterion (p. 16.1-337): moments within 3 %, deflections within 5 %.

Given in the figures: W14x48 (W360x72), major-axis bending, E = 29,000 ksi (200 000 MPa), length
28.0 ft (8.53 m); Case 1: 0.200 kip/ft (2.92 kN/m) transverse load; Case 2: 1.0 kip (4.45 kN) at the top.

NOT given: the section properties (I_x, A, shear area) and G.  No property is taken from memory here:

* I_x is CALIBRATED on the Commentary's own first-order (P = 0) value without shear deformation of
  Case 2, I = H L^3 / (3 E Delta_tip).  The P = 0 columns therefore test nothing new for I; the
  P > 0 columns test what the benchmark is for, the second-order amplification.
* The shear flexibility 1 / (G A_s) is CALIBRATED per case on the difference between the P = 0
  deflections with and without shear deformation (printed to one significant digit only), so the
  "with shear" rows are a consistency test of the amplification, not a test of a shear area.
* A does not enter these results in this theory (asserted below).

The US and the SI tables are independently rounded and are tested as two separate problems.
"""
import pytest

from benchfea.model import build_model
from benchfea.solve import analyze

MOMENT_TOL, DEFLECTION_TOL = 0.03, 0.05           # AISC 360-22 Commentary C2.1, p. 16.1-337

# (P, value with shear deformation, [value without shear deformation])
US = {"units": {"force": "kip", "length": "in"}, "E": 29000.0, "L": 28.0 * 12.0, "w": 0.200 / 12.0, "H": 1.0,
      "moment": 1.0, "force": 1.0,
      1: {"P": [0.0, 150.0, 300.0, 450.0],
          "M": [(235.0, 235.0), (270.0, 269.0), (316.0, 313.0), (380.0, 375.0)],
          "D": [(0.202, 0.197), (0.230, 0.224), (0.269, 0.261), (0.322, 0.311)]},
      2: {"P": [0.0, 100.0, 150.0, 200.0],
          "M": [(336.0, 336.0), (470.0, 469.0), (601.0, 598.0), (856.0, 848.0)],
          "D": [(0.907, 0.901), (1.34, 1.33), (1.77, 1.75), (2.60, 2.56)]}}
SI = {"units": {"force": "N", "length": "mm"}, "E": 200000.0, "L": 8530.0, "w": 2.92, "H": 4450.0,
      "moment": 1.0e6, "force": 1.0e3,
      1: {"P": [0.0, 667.0, 1330.0, 2000.0],
          "M": [(26.6, 26.6), (30.5, 30.4), (35.7, 35.4), (42.9, 42.4)],
          "D": [(5.13, 5.00), (5.84, 5.69), (6.83, 6.63), (8.18, 7.90)]},
      2: {"P": [0.0, 445.0, 667.0, 890.0],
          "M": [(38.0, 38.0), (53.1, 53.0), (67.9, 67.6), (96.7, 95.8)],
          "D": [(23.0, 22.9), (34.0, 33.8), (45.0, 44.5), (66.0, 65.0)]}}


def calibrated(table, case):
    """(I, G A_s) from the P = 0 columns of the Commentary (see the module docstring)."""
    L, E = table["L"], table["E"]
    inertia = table["H"] * L ** 3 / (3.0 * E * table[2]["D"][0][1])
    with_shear, without = table[case]["D"][0]
    if case == 1:
        gas = table["w"] * L ** 2 / (8.0 * (with_shear - without))       # midspan shear deflection w L^2 / (8 G As)
    else:
        gas = table["H"] * L / (with_shear - without)                    # tip shear deflection H L / (G As)
    return inertia, gas


def document(table, case, P, shear, area=None):
    L, E = table["L"], table["E"]
    inertia, gas = calibrated(table, case)
    G = E / 2.6                                          # only the product G As matters; As is set from it
    sec = {"type": "properties", "A": area or inertia / (L / 50.0) ** 2, "Iy": inertia, "Iz": inertia, "J": inertia / 100.0}
    if shear:
        sec["Asy"] = gas / G
    doc = {"formatVersion": "1.0.0", "units": table["units"], "conventions": {"verticalAxis": "Z"},
           "nodes": {"a": {"x": 0, "y": 0, "z": 0}, "b": {"x": 0, "y": 0, "z": L}},
           "materials": {"M": {"E": E, "G": G}}, "sections": {"S": sec},
           "members": {"m": {"i": "a", "j": "b", "section": "S", "material": "M"}},     # local y = +X: major axis = Iz
           "analysis": {"shearDeformation": shear, "resultStations": 3},
           "combinations": {"c": {"factors": {"c": 1.0}, "analysisType": "secondOrder"}}}
    if case == 1:                                        # pin at the bottom, lateral roller at the loaded top
        doc["supports"] = {"a": {"type": "pinned", "dofs": {"rz": "fixed"}}, "b": {"dofs": {"ux": "fixed", "uy": "fixed"}}}
        doc["loadCases"] = {"c": {"loads": [{"type": "nodal", "node": "b", "Fz": -P},
                                            {"type": "memberDistributed", "member": "m", "direction": "X", "w1": table["w"]}]}}
    else:                                                # fixed base, free top (held out of plane only)
        doc["supports"] = {"a": {"type": "fixed"}, "b": {"dofs": {"uy": "fixed", "rx": "fixed"}}}
        doc["loadCases"] = {"c": {"loads": [{"type": "nodal", "node": "b", "Fz": -P, "Fx": table["H"]}]}}
    return doc


def results(table, case, P, shear, area=None):
    """(moment, deflection) at the Commentary's locations, in the file units."""
    model = build_model(document(table, case, P * table["force"], shear, area))
    res = analyze(model)
    comb = res.to_document(model)["combinations"]["c"]
    assert comb["converged"] and comb["invariants"]["passed"]
    st = comb["members"]["m"]["stations"]
    if case == 1:
        return abs(st["Mz"][1]) / table["moment"], abs(st["uy"][1])        # midspan
    return abs(st["Mz"][0]) / table["moment"], abs(comb["displacements"]["b"]["ux"])   # base moment, tip deflection


@pytest.mark.parametrize("table", [US, SI], ids=["US", "SI"])
@pytest.mark.parametrize("case", [1, 2])
@pytest.mark.parametrize("shear", [False, True], ids=["without-shear", "with-shear"])
def test_commentary_benchmark(table, case, shear, record_property):
    worst_m = worst_d = 0.0
    for k, P in enumerate(table[case]["P"]):
        M, D = results(table, case, P, shear)
        M_ref, D_ref = table[case]["M"][k][0 if shear else 1], table[case]["D"][k][0 if shear else 1]
        assert M == pytest.approx(M_ref, rel=MOMENT_TOL), (P, M, M_ref)
        assert D == pytest.approx(D_ref, rel=DEFLECTION_TOL), (P, D, D_ref)
        worst_m, worst_d = max(worst_m, abs(M / M_ref - 1.0)), max(worst_d, abs(D / D_ref - 1.0))
    record_property("largest moment difference", worst_m)
    record_property("largest deflection difference", worst_d)
    # the actual differences are of the size of the rounding of the printed values (3 significant
    # figures; the SI shear calibration rests on a single digit): far inside the acceptance criterion
    assert worst_m <= (0.004 if table is US else 0.012), worst_m
    assert worst_d <= (0.009 if table is US else 0.03), worst_d


def test_the_benchmarks_do_not_depend_on_the_area():
    for case, P in ((1, 300.0), (2, 150.0)):
        a = results(US, case, P, True, area=5.0)
        b = results(US, case, P, True, area=50.0)
        assert a == pytest.approx(b, rel=1e-11)


def test_second_order_effect_is_what_is_being_tested():
    """A first-order analysis fails the benchmark by far more than the tolerance at the highest loads."""
    for case in (1, 2):
        M0, D0 = results(US, case, 0.0, False)
        M_ref, D_ref = US[case]["M"][-1][1], US[case]["D"][-1][1]
        assert abs(M0 / M_ref - 1.0) > 10 * MOMENT_TOL and abs(D0 / D_ref - 1.0) > 5 * DEFLECTION_TOL
