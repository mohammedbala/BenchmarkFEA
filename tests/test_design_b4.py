"""AISC 360-22 B3 (design basis), B4.1 (classification) and the design framework
(registry, section design view, check record).

Tier F tests:
(a) hand calculations written out as explicit arithmetic from the Specification
    (Tables B4.1a / B4.1b, pp. 16.1-21 to 16.1-23), independent of the code path;
(b) reproductions of printed limits of the AISC Design Examples v16.0
    (tests/design_examples; ONE tolerance policy, see design_examples.POLICY).
"""
import math

import pytest

from benchfea import design
from benchfea.design import DesignContext, DesignSection, Demand, registry
from benchfea.design.aisc360_22 import b4_classification as b4
from benchfea.design.base import refused, strength_check
from benchfea.model import build_model
from design_examples import assert_printed, example

E_SI, FY_SI = 200000.0, 345.0           # MPa
ROOT = math.sqrt(E_SI / FY_SI)          # sqrt(E/Fy) = 24.0772...


def _elements(cls):
    assert cls.status == "ok", cls.reason
    return {e.label: e for e in cls.elements}


# ============================================================================== tier F (a): hand calculations
def test_hand_rolled_i_compression_table_b41a_cases_1_and_5():
    # Rolled I 400 x 200, tf = 12, tw = 8, fillet r = 18 (mm), Fy = 345 MPa.
    sec = DesignSection.i_shape("rolled", A=7700.0, Ix=2.0e8, Iy=1.6e7, d=400.0, bf=200.0, tf=12.0, tw=8.0, r=18.0)
    el = _elements(b4.classify_compression(sec, E_SI, FY_SI))
    # B4.1a(a): b = bf/2 = 100; lambda = 100/12
    assert el["flange"].ratio == pytest.approx(100.0 / 12.0, rel=1e-14)
    # Table B4.1a Case 1: lambda_r = 0.56 sqrt(E/Fy) = 0.56 * 24.0772 = 13.483
    assert el["flange"].lambda_r == pytest.approx(0.56 * ROOT, rel=1e-14)
    assert 13.48 < el["flange"].lambda_r < 13.49
    assert (el["flange"].case, el["flange"].cls) == (1, "nonslender")          # 8.33 <= 13.48
    # B4.1b(a): h = clear distance between flanges less the fillet at each flange = 400 - 2*12 - 2*18 = 340
    assert el["web"].element.width == pytest.approx(340.0, rel=1e-14)
    assert el["web"].ratio == pytest.approx(340.0 / 8.0, rel=1e-14)            # 42.5
    # Table B4.1a Case 5: lambda_r = 1.49 sqrt(E/Fy) = 35.875
    assert el["web"].lambda_r == pytest.approx(1.49 * ROOT, rel=1e-14)
    assert (el["web"].case, el["web"].cls) == (5, "slender")                   # 42.5 > 35.875
    # Table E7.1 cases that go with the elements: unstiffened flange (c), stiffened web (a)
    assert el["flange"].e7 == ("(c)", 0.22, 1.49) and el["web"].e7 == ("(a)", 0.18, 1.31)


def test_hand_built_up_i_compression_table_b41a_case_2_with_kc():
    # Welded I 600 x 300, tf = 10, tw = 6: h = clear distance between flanges = 600 - 20 = 580 (B4.1b(b))
    sec = DesignSection.i_shape("builtUp", A=9480.0, Ix=5.0e8, Iy=4.5e7, d=600.0, bf=300.0, tf=10.0, tw=6.0)
    cls = b4.classify_compression(sec, E_SI, FY_SI)
    el = _elements(cls)
    h_tw = 580.0 / 6.0                                                          # 96.667
    assert el["web"].ratio == pytest.approx(h_tw, rel=1e-14)
    # footnote [a]: kc = 4 / sqrt(h/tw) = 4 / 9.8319 = 0.40684 (inside 0.35 ... 0.76)
    kc = 4.0 / math.sqrt(h_tw)
    assert 0.4068 < kc < 0.4069 and el["flange"].kc == pytest.approx(kc, rel=1e-14)
    # Case 2: lambda_r = 0.64 sqrt(kc E / Fy) = 0.64 * sqrt(235.85) = 9.829
    lam_r = 0.64 * math.sqrt(kc * E_SI / FY_SI)
    assert 9.82 < lam_r < 9.84 and el["flange"].lambda_r == pytest.approx(lam_r, rel=1e-14)
    assert (el["flange"].case, el["flange"].ratio, el["flange"].cls) == (2, 15.0, "slender")    # 150/10 = 15 > 9.83
    assert el["web"].cls == "slender" and cls.slender                           # 96.7 > 35.875


def test_hand_kc_limits_footnote_a():
    assert b4.kc_factor(20.0) == 0.76             # 4/sqrt(20) = 0.894 -> not greater than 0.76
    assert b4.kc_factor(200.0) == 0.35            # 4/sqrt(200) = 0.283 -> not less than 0.35
    assert b4.kc_factor(64.0) == pytest.approx(0.5, rel=1e-15)                  # 4/8


def test_hand_rectangular_hss_and_box_compression_cases_6_and_8():
    # HSS 200 x 100 x 5 (t = design wall thickness, an input): corner radius not known ->
    # B4.1b(d): h = 200 - 3*5 = 185, b = 100 - 3*5 = 85; ratios 37 and 17
    hss = DesignSection.rectangular("rectHSS", A=2800.0, Ix=1.4e7, Iy=4.6e6, H=200.0, B=100.0, t=5.0)
    el = _elements(b4.classify_compression(hss, E_SI, FY_SI))
    assert el["wall (depth)"].ratio == pytest.approx(185.0 / 5.0) and el["wall (width)"].ratio == pytest.approx(17.0)
    # Case 6: lambda_r = 1.40 sqrt(E/Fy) = 33.708
    assert el["wall (depth)"].lambda_r == pytest.approx(1.40 * ROOT, rel=1e-14)
    assert [e.cls for e in el.values()] == ["slender", "nonslender"] and el["wall (depth)"].case == 6
    assert el["wall (depth)"].e7 == ("(b)", 0.20, 1.38)
    # Box section of the same outside size: B4.1b(e) clear distance: h = 200 - 2*5 = 190, b = 90; ratios 38, 18
    box = DesignSection.rectangular("box", A=2900.0, Ix=1.5e7, Iy=5.0e6, H=200.0, B=100.0, t=5.0)
    el = _elements(b4.classify_compression(box, E_SI, FY_SI))
    assert el["wall (depth)"].ratio == pytest.approx(38.0) and el["wall (width)"].ratio == pytest.approx(18.0)
    # Case 8: lambda_r = 1.49 sqrt(E/Fy) = 35.875 -> 38 slender, 18 nonslender
    assert el["wall (depth)"].lambda_r == pytest.approx(1.49 * ROOT, rel=1e-14)
    assert [e.cls for e in el.values()] == ["slender", "nonslender"] and el["wall (depth)"].case == 8
    assert el["wall (depth)"].e7 == ("(a)", 0.18, 1.31)


def test_hand_round_hss_compression_case_9():
    # Case 9: lambda_r = 0.11 E/Fy = 0.11 * 579.71 = 63.77 (no square root)
    lam_r = 0.11 * E_SI / FY_SI
    assert 63.7 < lam_r < 63.8
    a = _elements(b4.classify_compression(DesignSection.round_hss(A=4600.0, I=5.0e7, D=300.0, t=5.0), E_SI, FY_SI))
    assert (a["wall"].ratio, a["wall"].cls, a["wall"].case) == (60.0, "nonslender", 9)
    assert a["wall"].lambda_r == pytest.approx(lam_r, rel=1e-14)
    b = _elements(b4.classify_compression(DesignSection.round_hss(A=6200.0, I=1.2e8, D=400.0, t=5.0), E_SI, FY_SI))
    assert (b["wall"].ratio, b["wall"].cls) == (80.0, "slender")


def test_hand_flexure_table_b41b_i_shapes_cases_10_11_13_15():
    rolled = DesignSection.i_shape("rolled", A=7700.0, Ix=2.0e8, Iy=1.6e7, d=400.0, bf=200.0, tf=12.0, tw=8.0, r=18.0)
    el = _elements(b4.classify_flexure(rolled, E_SI, FY_SI, "major"))
    # Case 10: lambda_p = 0.38 sqrt(E/Fy) = 9.149, lambda_r = 1.0 sqrt(E/Fy) = 24.077; flange 8.33 -> compact
    assert el["flange"].lambda_p == pytest.approx(0.38 * ROOT, rel=1e-14)
    assert el["flange"].lambda_r == pytest.approx(1.0 * ROOT, rel=1e-14)
    assert (el["flange"].case, el["flange"].cls) == (10, "compact")
    # Case 15: lambda_p = 3.76 sqrt(E/Fy) = 90.53, lambda_r = 5.70 sqrt(E/Fy) = 137.24; web 42.5 -> compact
    assert el["web"].lambda_p == pytest.approx(3.76 * ROOT, rel=1e-14)
    assert el["web"].lambda_r == pytest.approx(5.70 * ROOT, rel=1e-14)
    assert (el["web"].case, el["web"].cls) == (15, "compact")
    # minor axis: Case 13 for the flanges only (Table B4.1b has no case for the web on the neutral axis)
    minor = b4.classify_flexure(rolled, E_SI, FY_SI, "minor")
    assert [(e.label, e.case, e.cls) for e in minor.elements] == [("flange", 13, "compact")]

    built = DesignSection.i_shape("builtUp", A=9480.0, Ix=5.0e8, Iy=4.5e7, d=600.0, bf=300.0, tf=10.0, tw=6.0)
    cls = b4.classify_flexure(built, E_SI, FY_SI, "major")
    el = _elements(cls)
    # Case 11: lambda_r = 0.95 sqrt(kc E / FL), kc = 4/sqrt(96.667) = 0.40684, FL = 0.7 Fy = 241.5 (footnote [b],
    # doubly symmetric: Sxt/Sxc = 1 >= 0.7): 0.95 * sqrt(0.40684 * 200000 / 241.5) = 0.95 * 18.356 = 17.44
    kc = 4.0 / math.sqrt(580.0 / 6.0)
    lam_r = 0.95 * math.sqrt(kc * E_SI / (0.7 * FY_SI))
    assert 17.43 < lam_r < 17.45 and el["flange"].lambda_r == pytest.approx(lam_r, rel=1e-14)
    assert (el["flange"].case, el["flange"].cls) == (11, "noncompact")          # 9.149 < 15 <= 17.44
    assert el["web"].cls == "noncompact" and cls.cls == "noncompact"            # 90.53 < 96.67 <= 137.24


def test_hand_flexure_table_b41b_hss_box_round_cases_17_19_20_21():
    hss = DesignSection.rectangular("rectHSS", A=2800.0, Ix=1.4e7, Iy=4.6e6, H=200.0, B=100.0, t=5.0)
    major = _elements(b4.classify_flexure(hss, E_SI, FY_SI, "major"))
    # major axis: flange = width wall (b/t = 17), web = depth wall (h/t = 37)
    f, w = major["flange (width wall)"], major["web (depth wall)"]
    assert (f.case, f.ratio, f.cls) == (17, 17.0, "compact")                    # lambda_p = 1.12 * 24.077 = 26.97
    assert f.lambda_p == pytest.approx(1.12 * ROOT, rel=1e-14) and f.lambda_r == pytest.approx(1.40 * ROOT, rel=1e-14)
    assert (w.case, w.ratio, w.cls) == (19, 37.0, "compact")                    # lambda_p = 2.42 * 24.077 = 58.27
    assert w.lambda_p == pytest.approx(2.42 * ROOT, rel=1e-14) and w.lambda_r == pytest.approx(5.70 * ROOT, rel=1e-14)
    minor = _elements(b4.classify_flexure(hss, E_SI, FY_SI, "minor"))
    # minor axis: the depth walls are the flanges: 37 > lambda_r = 1.40 * 24.077 = 33.71 -> slender
    assert minor["flange (depth wall)"].cls == "slender" and minor["web (width wall)"].cls == "compact"

    box = DesignSection.rectangular("box", A=2900.0, Ix=1.5e7, Iy=5.0e6, H=200.0, B=100.0, t=5.0)
    f = _elements(b4.classify_flexure(box, E_SI, FY_SI, "minor"))["flange (depth wall)"]
    # Case 21: lambda_r = 1.49 sqrt(E/Fy) = 35.875; clear width 190/5 = 38 -> slender
    assert (f.case, f.ratio, f.cls) == (21, 38.0, "slender") and f.lambda_r == pytest.approx(1.49 * ROOT, rel=1e-14)

    # Case 20: lambda_p = 0.07 E/Fy = 40.58, lambda_r = 0.31 E/Fy = 179.71; D/t = 60 -> noncompact
    r = _elements(b4.classify_flexure(DesignSection.round_hss(A=4600.0, I=5.0e7, D=300.0, t=5.0), E_SI, FY_SI))["wall"]
    assert r.lambda_p == pytest.approx(0.07 * E_SI / FY_SI, rel=1e-14)
    assert r.lambda_r == pytest.approx(0.31 * E_SI / FY_SI, rel=1e-14)
    assert (r.case, r.cls) == (20, "noncompact")


def test_boundary_equal_to_limit_is_the_favourable_class():
    """B4.1: slender only if the ratio EXCEEDS lambda_r; compact if it does not exceed lambda_p."""
    at = 0.56 * ROOT
    sec = DesignSection.i_shape("rolled", A=1.0, rx=1.0, ry=1.0, flange_ratio=at, web_ratio=10.0)
    assert _elements(b4.classify_compression(sec, E_SI, FY_SI))["flange"].cls == "nonslender"
    above = DesignSection.i_shape("rolled", A=1.0, rx=1.0, ry=1.0, flange_ratio=math.nextafter(at, math.inf),
                                  web_ratio=10.0)
    assert _elements(b4.classify_compression(above, E_SI, FY_SI))["flange"].cls == "slender"
    lam_p = 0.38 * ROOT
    sec = DesignSection.i_shape("rolled", A=1.0, rx=1.0, ry=1.0, flange_ratio=lam_p, web_ratio=10.0)
    assert _elements(b4.classify_flexure(sec, E_SI, FY_SI))["flange"].cls == "compact"


# ============================================================================== tier F (b): Design Examples
def test_de_e1d_and_e1e_w_shape_limits_table_b41a():
    """Examples E.1D (E-10, PDF 74) and E.1E (E-13..E-14, PDF 77-78): Cases 1 and 5."""
    for eid, web_class in (("E.1D", "nonslender"), ("E.1E", "slender"), ("E.1C", "nonslender")):
        x = example(eid)
        i, p = x["inputs"], x["printed"]
        sec = DesignSection.i_shape("rolled", A=i["Ag"], rx=i["rx"], ry=i["ry"], flange_ratio=i["bf_2tf"],
                                    web_ratio=i["h_tw"])
        el = _elements(b4.classify_compression(sec, i["E"], i["Fy"]))
        assert_printed(el["flange"].lambda_r, p["flange_lambda_r"], f"{eid} flange lambda_r")
        assert_printed(el["web"].lambda_r, p["web_lambda_r"], f"{eid} web lambda_r")
        assert (el["flange"].cls, el["web"].cls) == ("nonslender", web_class)


def test_de_e2_and_e3_built_up_limits_with_kc():
    """Examples E.2 (E-22..E-23, PDF 86-87) and E.3 (E-26..E-27, PDF 90-91): Case 2 with kc, Case 5."""
    for eid, classes in (("E.2", ("nonslender", "slender")), ("E.3", ("slender", "nonslender"))):
        x = example(eid)
        i, p = x["inputs"], x["printed"]
        sec = DesignSection.i_shape("builtUp", A=i["Ag"], Ix=i["Ix"], Iy=i["Iy"], d=i["d"], bf=i["bf"], tf=i["tf"],
                                    tw=i["tw"], h=i["h"])
        el = _elements(b4.classify_compression(sec, i["E"], i["Fy"]))
        assert_printed(el["flange"].ratio, p["flange_lambda"], f"{eid} flange lambda")
        assert_printed(el["web"].ratio, p["web_lambda"], f"{eid} web lambda")
        assert_printed(el["flange"].kc, p["kc"], f"{eid} kc")
        assert_printed(el["flange"].lambda_r, p["flange_lambda_r"], f"{eid} flange lambda_r")
        assert_printed(el["web"].lambda_r, p["web_lambda_r"], f"{eid} web lambda_r")
        assert (el["flange"].cls, el["web"].cls) == classes and el["flange"].case == 2


def test_de_e9_e10_e11_hss_and_pipe_limits():
    """Examples E.9 (E-60, PDF 124), E.10 (E-62, PDF 126): Case 6; E.11 (E-68, PDF 132): Case 9."""
    for eid, classes in (("E.9", ["nonslender", "nonslender"]), ("E.10", ["slender", "slender"])):
        x = example(eid)
        i = x["inputs"]
        sec = DesignSection.rectangular("rectHSS", A=i["Ag"], rx=i["rx"], ry=i["ry"], h_ratio=i["h_t"],
                                        b_ratio=i["b_t"])
        cls = b4.classify_compression(sec, i["E"], i["Fy"])
        assert [e.cls for e in cls.elements] == classes
        for e in cls.elements:
            assert_printed(e.lambda_r, x["printed"]["lambda_r"], f"{eid} lambda_r")
    x = example("E.11")
    i = x["inputs"]
    cls = b4.classify_compression(DesignSection.round_hss(A=i["Ag"], r=i["r"], ratio=i["D_t"]), i["E"], i["Fy"])
    assert_printed(cls.elements[0].lambda_r, x["printed"]["lambda_r"], "E.11 lambda_r")
    assert cls.elements[0].cls == "nonslender"


def test_de_f3b_f7b_f9b_flexure_limits_table_b41b():
    """Examples F.3B (F-24, PDF 186): Case 10; F.7B (F-34..F-35, PDF 196-197): Cases 17, 19;
    F.9B (F-43, PDF 205): Case 20."""
    x = example("F.3B")
    i, p = x["inputs"], x["printed"]
    sec = DesignSection.i_shape("rolled", A=1.0, rx=1.0, ry=1.0, flange_ratio=i["bf_2tf"])
    cls = b4.classify_flexure(sec, i["E"], i["Fy"], "minor")        # flange only: the example prints no h/tw
    # (Case 13, minor axis, has the same limits as Case 10; the major-axis call needs the web ratio)
    assert cls.elements[0].cls == p["flange_class"]
    sec = DesignSection.i_shape("rolled", A=1.0, rx=1.0, ry=1.0, flange_ratio=i["bf_2tf"], web_ratio=1.0)
    f = _elements(b4.classify_flexure(sec, i["E"], i["Fy"], "major"))["flange"]
    assert f.case == 10 and f.cls == p["flange_class"]
    assert_printed(f.lambda_p, p["lambda_pf"], "F.3B lambda_pf")
    assert_printed(f.lambda_r, p["lambda_rf"], "F.3B lambda_rf")

    x = example("F.7B")
    i, p = x["inputs"], x["printed"]
    sec = DesignSection.rectangular("rectHSS", A=1.0, rx=1.0, ry=1.0, h_ratio=i["h_t"], b_ratio=i["b_t"])
    el = _elements(b4.classify_flexure(sec, i["E"], i["Fy"], "major"))
    f, w = el["flange (width wall)"], el["web (depth wall)"]
    assert_printed(f.lambda_p, p["flange_lambda_p"], "F.7B flange lambda_p")
    assert_printed(f.lambda_r, p["flange_lambda_r"], "F.7B flange lambda_r")
    assert_printed(w.lambda_p, p["web_lambda_p"], "F.7B web lambda_p")
    assert (f.cls, w.cls) == (p["flange_class"], p["web_class"])

    x = example("F.9B")
    i, p = x["inputs"], x["printed"]
    r = b4.classify_flexure(DesignSection.round_hss(A=1.0, r=1.0, ratio=i["D_t"]), i["E"], i["Fy"]).elements[0]
    assert_printed(r.lambda_p, p["lambda_p"], "F.9B lambda_p")
    assert r.cls == p["class"]


# ============================================================================== not classifiable -> notChecked
def test_missing_dimension_is_not_checked_never_a_silent_class():
    # rolled shape without fillet dimension, h or h/tw: web ratio unknown
    sec = DesignSection.i_shape("rolled", A=7700.0, Ix=2.0e8, Iy=1.6e7, d=400.0, bf=200.0, tf=12.0, tw=8.0)
    cls = b4.classify_compression(sec, E_SI, FY_SI)
    assert cls.status == "notChecked" and "web" in cls.reason and cls.cls is None
    assert b4.classify_flexure(sec, E_SI, FY_SI, "major").status == "notChecked"
    assert b4.classify_flexure(sec, E_SI, FY_SI, "minor").status == "ok"       # flanges only
    # no Fy
    assert b4.classify_compression(DesignSection.round_hss(A=1.0, r=1.0, ratio=30.0), E_SI, None).status == "notChecked"
    # shape outside the scope
    out = DesignSection.unsupported("angle", "L1")
    assert b4.classify_compression(out, E_SI, FY_SI).reason == "angle"
    with pytest.raises(ValueError):
        b4.classify_flexure(sec, E_SI, FY_SI, "x")


# ============================================================================== section design view from the model
def _model(sections, design=None, units=None):
    data = {
        "formatVersion": "1.0.0", "units": units or {"force": "N", "length": "mm"},
        "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0, "z": 0}, "B": {"x": 0, "y": 0, "z": 4000}},
        "supports": {"A": {"type": "fixed"}},
        "materials": {"S": {"E": 200000, "nu": 0.3, "Fy": 345, "Fu": 450}},
        "sections": sections,
        "members": {name: {"i": "A", "j": "B", "section": name, "material": "S"} for name in sections},
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fz": -1000}]}},
        "combinations": {"U1": {"factors": {"D": 1.4}}},
    }
    if design:
        data["design"] = design
    return build_model(data)


def test_design_view_of_model_sections_and_its_stated_mapping():
    m = _model({"rolled": {"type": "I", "d": 400, "bf": 200, "tf": 12, "tw": 8, "r": 18},
                "welded": {"type": "I", "d": 600, "bf": 300, "tf": 10, "tw": 6},
                "tube": {"type": "box", "H": 200, "B": 100, "t": 5},
                "chs": {"type": "pipe", "D": 300, "t": 5},
                "bar": {"type": "rectangle", "h": 100, "b": 50}})
    rolled = DesignSection.from_model(m.sections["rolled"])
    assert (rolled.shape, rolled.fabrication) == ("I", "rolled") and "fillet radius" in rolled.source("fabrication")
    assert rolled.elements["web"].width == pytest.approx(340.0) and rolled.source("A") == "section:rolled"
    assert (rolled.Ix, rolled.Iy) == (m.sections["rolled"].Iz, m.sections["rolled"].Iy)     # x = major = local z
    welded = DesignSection.from_model(m.sections["welded"])
    assert welded.fabrication == "builtUp" and welded.elements["web"].width == pytest.approx(580.0)
    assert any("built-up" in n for n in welded.notes)
    tube = DesignSection.from_model(m.sections["tube"])
    assert tube.shape == "box" and tube.elements["wall_y"].width == pytest.approx(190.0)
    assert any("not as an HSS" in n for n in tube.notes)
    chs = DesignSection.from_model(m.sections["chs"])
    assert chs.shape == "roundHSS" and chs.elements["wall"].ratio == pytest.approx(60.0)
    assert any("0.93" in n and "NOT applied" in n for n in chs.notes)          # B4.2: t is an input
    bar = DesignSection.from_model(m.sections["bar"])
    assert bar.shape is None and "outside the implemented scope" in bar.reason and bar.A == pytest.approx(5000.0)


def test_design_view_of_database_shapes(tmp_path, monkeypatch):
    """A user database row: W -> rolled with h = d - 2 kdes; HSS -> rectangular HSS with tdes; without the
    kdes column the web depth is unknown (notChecked), never assumed."""
    rows = ["AISC_Manual_Label,Type,A,d,bf,tw,tf,kdes,Ht,B,tdes,OD,Ix,Zx,Sx,Iy,Zy,Sy,J,Cw",
            "WTEST,W,10,16,6,0.3,0.5,0.9,,,,,400,60,50,20,8,6,0.5,1200",
            "HSSTEST,HSS,6,,,,,,6,4,0.35,,28,11,9,15,8,7,30,0"]
    db = tmp_path / "db.csv"
    db.write_text("\n".join(rows))
    monkeypatch.setenv("BENCHFEA_SHAPES_DB", str(db))
    m = _model({"w": {"type": "standard", "shape": "WTEST"}, "h": {"type": "standard", "shape": "HSSTEST"}},
               units={"force": "kip", "length": "in"})
    w = DesignSection.from_model(m.sections["w"])
    assert (w.shape, w.fabrication) == ("I", "rolled") and "database" in w.source("fabrication")
    assert w.elements["web"].width == pytest.approx((16 - 2 * 0.9) * 25.4)    # h = d - 2 kdes, internal mm
    h = DesignSection.from_model(m.sections["h"])
    assert h.shape == "rectHSS" and h.elements["wall_y"].width == pytest.approx((6 - 3 * 0.35) * 25.4)

    db.write_text("\n".join(",".join(c for k, c in enumerate(r.split(",")) if k != 7) for r in rows))
    m = _model({"w": {"type": "standard", "shape": "WTEST"}}, units={"force": "kip", "length": "in"})
    w = DesignSection.from_model(m.sections["w"])
    assert w.elements["web"].ratio is None and any("kdes" in n for n in w.notes)
    assert b4.classify_compression(w, 200000.0, 345.0).status == "notChecked"


# ============================================================================== B3 design basis and the record
def _ctx(method="LRFD", **kw):
    sec = DesignSection.round_hss(A=1000.0, r=30.0, ratio=20.0, J=1.8e6)       # J = 2 I = 2 A r^2
    base = dict(member="M", section=sec, method=method, E=E_SI, G=77000.0, Fy=FY_SI, Fu=450.0, L=3000.0,
                Lx=3000.0, Ly=3000.0, Lz=3000.0)
    base.update(kw)
    return DesignContext(**base)


def test_b3_lrfd_and_asd_available_strength():
    # B3.1 Eq. B3-1: Ru <= phi Rn -> available = 0.90 * 1000 = 900; ratio = 450/900 = 0.5
    c = strength_check("x", "AISC 360-22 X", _ctx("LRFD"), nominal=1000.0, phi=0.90, omega=1.67, demand=450.0)
    assert (c.capacity, c.ratio, c.passed, c.factor, c.factor_name, c.method) == (900.0, 0.5, True, 0.90, "phi", "LRFD")
    assert c.clause.endswith("B3.1 Eq. B3-1") and c.status == "ok" and c.nominal == 1000.0
    # B3.2 Eq. B3-2: Ra <= Rn / Omega -> available = 1000/1.67 = 598.80; ratio = 600/598.80 > 1 -> fails
    c = strength_check("x", "AISC 360-22 X", _ctx("ASD"), nominal=1000.0, phi=0.90, omega=1.67, demand=600.0)
    assert c.capacity == pytest.approx(1000.0 / 1.67, rel=1e-15) and c.factor_name == "Omega" and c.factor == 1.67
    assert c.ratio == pytest.approx(600.0 * 1.67 / 1000.0, rel=1e-15) and c.passed is False
    assert c.clause.endswith("B3.2 Eq. B3-2")
    # demand exactly equal to the available strength passes (B3-1 is "<=")
    assert strength_check("x", "c", _ctx(), nominal=1000.0, phi=0.9, omega=1.67, demand=900.0).passed is True
    # anything but LRFD / ASD is refused, not guessed
    c = strength_check("x", "c", _ctx("LSD"), nominal=1000.0, phi=0.9, omega=1.67, demand=1.0)
    assert c.status == "refused" and c.capacity is None and c.passed is None and "LSD" in c.reason


def test_registry_has_one_edition_and_refuses_any_other():
    assert registry.registered() == ("AISC 360-22",)
    ed = registry.get("AISC 360-22")
    assert ed.code == "AISC 360-22" and len(ed.limit_states) >= 8
    for other in ("AISC 360-16", "AISC 360", "Eurocode 3", "", None):
        with pytest.raises(design.UnknownDesignCode) as err:
            registry.get(other)
        assert "AISC 360-22" in str(err.value)
    with pytest.raises(ValueError):
        registry.register(registry.Edition("AISC 360-22", "duplicate", ()))
    res = design.check_member(_ctx(code="AISC 360-16"), [Demand("U1", Pt=1.0)])
    assert res.passed is False and res.governing is None
    assert [c.status for c in res.checks] == ["refused"] and "AISC 360-16" in res.checks[0].reason
    assert res.not_checked and "not available" in res.not_checked[0]


def test_every_check_carries_the_full_record():
    # Lz > Lx = Ly so that torsional buckling (E4) is one of the checks
    res = design.check_member(_ctx(Lz=4000.0), [Demand("U1", Pt=1.0e5, Pc=1.0e5, V=3.0e4)])
    assert len(res.checks) == 8
    assert len(design.check_member(_ctx(), [Demand("U1", Pt=1.0e5, Pc=1.0e5, V=3.0e4)]).checks) == 7
    for c in res.checks:
        assert c.name and c.clause.startswith("AISC 360-22") and c.status == "ok" and c.method == "LRFD"
        assert c.inputs and all(len(v) == 3 and isinstance(v[2], str) and v[2] for v in c.inputs.values())
        assert isinstance(c.passed, bool) and c.ratio == pytest.approx(c.demand / c.capacity, rel=1e-14)
        if not c.advisory:
            assert c.combination == "U1" and c.factor_name == "phi" and c.details
            assert c.capacity == pytest.approx(c.factor * c.nominal, rel=1e-14)
            assert "Eq." in c.clause and "B3.1 Eq. B3-1" in c.clause
    assert res.passed and res.governing is max((c for c in res.checks if not c.advisory), key=lambda c: c.ratio)


def test_governing_combination_and_blocked_limit_states():
    demands = [Demand("U1", Pt=1.0e5), Demand("U2", Pt=2.5e5), Demand("U3", Pt=5.0e4)]
    res = design.check_member(_ctx(), demands)
    y = next(c for c in res.checks if c.name == "tension: yielding in the gross section")
    assert (y.combination, y.demand) == ("U2", 2.5e5)
    # a failing advisory check does not fail the member; a failing strength check does
    res = design.check_member(_ctx(maxSlendernessTension=50.0), demands)
    adv = next(c for c in res.checks if c.name.startswith("tension: slenderness"))
    assert adv.advisory and adv.passed is False and res.passed is True
    assert design.check_member(_ctx(), [Demand("U1", Pt=1.0e6)]).passed is False
    # a limit state that cannot be checked blocks the member only when it has a demand
    no_fu = _ctx(Fu=None)
    res = design.check_member(no_fu, [Demand("U1", Pc=1.0e4)])
    assert res.passed is True and any("rupture" in s and "notChecked" in s for s in res.not_checked)
    res = design.check_member(no_fu, [Demand("U1", Pc=1.0e4), Demand("U2", Pt=1.0)])
    rup = next(c for c in res.checks if "rupture" in c.name)
    assert res.passed is False and (rup.status, rup.combination, rup.passed) == ("notChecked", "U2", None)
    # no combination at all: strengths are still reported, with a warning
    res = design.check_member(_ctx(), [])
    assert all(c.ratio == 0.0 for c in res.checks if not c.advisory) and any("zero demand" in w for w in res.warnings)
    assert refused("n", "c", "why").status == "refused"


def test_context_from_model_carries_provenance_and_bracing():
    m = _model({"col": {"type": "I", "d": 400, "bf": 200, "tf": 12, "tw": 8, "r": 18}},
               design={"defaults": {"method": "ASD"},
                       "members": {"col": {"Ky": 0.8, "bracing": {"minor": [0.25]}, "An": 7000, "U": 0.9}}})
    ctx = DesignContext.from_model(m, "col")
    assert (ctx.method, ctx.source("method")) == ("ASD", "global")
    assert (ctx.Ky, ctx.source("Ky"), ctx.Kx, ctx.source("Kx")) == (0.8, "member", 1.0, "default")
    assert ctx.source("Fy") == "material:S" and ctx.Fy == 345.0
    assert ctx.Lx == pytest.approx(4000.0) and ctx.source("Lx") == "default"
    # brace at 0.25 L about the minor axis: longest segment = 0.75 * 4000
    assert ctx.Ly == pytest.approx(3000.0) and ctx.source("Ly").startswith("bracing")
    assert (ctx.An, ctx.U, ctx.source("U")) == (7000.0, 0.9, "member")
    assert ctx.inputs("Ky", "Fy") == {"Ky": (0.8, "dimensionless", "member"),
                                      "Fy": (345.0, "force/length^2", "material:S")}


def test_design_document_of_every_shape_family_is_plain_json(fixtures):
    """All checks of all section kinds of a kip-in fixture convert to the file's units and serialise."""
    import json

    from benchfea.results import DesignResults

    m = build_model(fixtures / "valid" / "frame_kip_in.json")
    kip = 4448.2216152605
    demand = Demand("U1", Pt=10 * kip, Pc=20 * kip, Vy=5 * kip, Vz=1 * kip, V=5.1 * kip, Lv=1000.0)
    members = {dm: design.check_member(DesignContext.from_model(m, dm), [demand]) for dm in m.design_members}
    shapes = {DesignContext.from_model(m, dm).section.shape for dm in m.design_members}
    assert shapes == {"I", "box", "roundHSS", None}
    doc = DesignResults("AISC 360-22", members).to_document(m)
    text = json.dumps(doc)
    assert "NaN" not in text and "Infinity" not in text
    for dm, res in doc["members"].items():
        for c in res["checks"]:
            assert c["status"] in ("ok", "notChecked", "refused") and c["clause"].startswith("AISC 360-22")
            if c["status"] == "ok":
                assert c["capacity"] > 0 and c["ratio"] == pytest.approx(c["demand"] / c["capacity"], rel=1e-12)
                assert all(set(v) == {"value", "unit", "source"} for v in c["inputs"].values())
            else:
                assert c["reason"] and c["capacity"] is None and c["ratio"] is None and c["passed"] is None
    # stresses are reported in the file's units (kip/in^2): E of the steel is 29000 there
    some = next(c for c in doc["members"]["C1"]["checks"] if "flexural buckling" in c["name"])
    assert some["inputs"]["E"]["unit"] == "kip/in^2" and some["inputs"]["E"]["value"] == pytest.approx(29000.0)
    assert some["details"]["Fe"]["unit"] == "kip/in^2"
