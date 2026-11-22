"""The resolved internal model: local axes, normalised loads, design-variable
resolution with provenance, defaults-used report, shapes database reader."""
import copy
import math

import numpy as np
import pytest

from benchfea import io
from benchfea import shapes
from benchfea.geometry import GeometryError, local_triad, member_geometry
from benchfea.model import DOFS, RIGID, build_model
from benchfea.problems import ModelInvalid

from helpers import assert_same

KIP, IN = 4448.2216152605, 25.4


# ------------------------------------------------------------------------- local axes
@pytest.mark.parametrize("vertical,pj,x,y,z", [
    ("Z", (5, 0, 0), (1, 0, 0), (0, 0, 1), (0, -1, 0)),
    ("Z", (0, 5, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0)),
    ("Z", (-5, 0, 0), (-1, 0, 0), (0, 0, 1), (0, 1, 0)),
    ("Z", (0, 0, 5), (0, 0, 1), (1, 0, 0), (0, 1, 0)),          # vertical up: local y = global +X
    ("Z", (0, 0, -5), (0, 0, -1), (1, 0, 0), (0, -1, 0)),       # vertical down: local y = global +X
    ("Z", (3, 0, 4), (0.6, 0, 0.8), (-0.8, 0, 0.6), (0, -1, 0)),
    ("Y", (5, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ("Y", (0, 0, 5), (0, 0, 1), (0, 1, 0), (-1, 0, 0)),
    ("Y", (0, 5, 0), (0, 1, 0), (1, 0, 0), (0, 0, -1)),         # vertical (Y up): local y = global +X
    ("Y", (3, 4, 0), (0.6, 0.8, 0), (-0.8, 0.6, 0), (0, 0, 1)),
])
def test_default_local_axes(vertical, pj, x, y, z):
    R = local_triad((0, 0, 0), pj, vertical)
    np.testing.assert_allclose(R, [x, y, z], atol=1e-15)


def test_roll_rotates_y_and_z_about_x_by_the_right_hand_rule():
    R = local_triad((0, 0, 0), (5, 0, 0), "Z", roll_deg=90)
    np.testing.assert_allclose(R, [(1, 0, 0), (0, -1, 0), (0, 0, -1)], atol=1e-15)   # y' = z, z' = -y
    R30 = local_triad((0, 0, 0), (5, 0, 0), "Z", roll_deg=30)
    c, s = math.cos(math.pi / 6), math.sin(math.pi / 6)
    np.testing.assert_allclose(R30[1], (0, -s, c), atol=1e-15)
    np.testing.assert_allclose(R30[2], (0, -c, -s), atol=1e-15)


def test_reference_vector_defines_the_xy_plane():
    R = local_triad((0, 0, 0), (5, 0, 0), "Z", reference_vector=(0, 2, 0))
    np.testing.assert_allclose(R, [(1, 0, 0), (0, 1, 0), (0, 0, 1)], atol=1e-15)
    R = local_triad((0, 0, 0), (5, 0, 0), "Z", reference_vector=(7, 0, 3))    # not perpendicular: projected
    np.testing.assert_allclose(R, [(1, 0, 0), (0, 0, 1), (0, -1, 0)], atol=1e-15)
    R = local_triad((0, 0, 0), (0, 0, 5), "Z", reference_vector=(0, 1, 0))    # column, web along Y
    np.testing.assert_allclose(R, [(0, 0, 1), (0, 1, 0), (-1, 0, 0)], atol=1e-15)


def test_triads_are_right_handed_and_orthonormal():
    rng = np.random.default_rng(1)
    for _ in range(200):
        a, b = rng.normal(size=3), rng.normal(size=3)
        R = local_triad(a, b, "Z", roll_deg=float(rng.uniform(-180, 180)))
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(R) == pytest.approx(1.0)
        np.testing.assert_allclose(R[0], (b - a) / np.linalg.norm(b - a), atol=1e-12)
        R0 = local_triad(a, b, "Z")
        assert R0[1][2] > 0                                    # default: local y points to the +vertical side
        assert R0[2][2] == pytest.approx(0.0, abs=1e-12)       # and local z is horizontal


def test_triad_errors():
    with pytest.raises(GeometryError) as e:
        local_triad((1, 2, 3), (1, 2, 3))
    assert e.value.code == "member.zero-length"
    for ref in ((2, 0, 0), (-1, 0, 0), (0, 0, 0)):
        with pytest.raises(GeometryError) as e:
            local_triad((0, 0, 0), (5, 0, 0), reference_vector=ref)
        assert e.value.code == "member.orientation-collinear"


def test_offsets_local_and_global_are_equivalent():
    nodes = {"a": {"x": 0, "y": 0, "z": 0}, "b": {"x": 3, "y": 0, "z": 4}}
    loc = {"i": "a", "j": "b", "offsets": {"axes": "local", "i": [0.5, 0, 0.1], "j": [-1, 0.2, 0]}}
    # local axes of the node line: x = (.6,0,.8), y = (-.8,0,.6), z = (0,-1,0)
    glo = {"i": "a", "j": "b", "offsets": {"axes": "global", "i": [0.3, -0.1, 0.4], "j": [-0.76, 0, -0.68]}}
    ga, gb = member_geometry(loc, nodes, "Z"), member_geometry(glo, nodes, "Z")
    for p, q in zip(ga, gb):
        np.testing.assert_allclose(p, q, atol=1e-14)
    start, end, R, length = ga[:4]
    np.testing.assert_allclose(start, (0.3, -0.1, 0.4), atol=1e-15)
    assert length == pytest.approx(np.linalg.norm(np.array([2.24, 0, 3.32]) - np.array([0.3, -0.1, 0.4])))
    np.testing.assert_allclose(R[0], (end - start) / length)   # local x follows the flexible length


# ------------------------------------------------------------------- the fixture model
@pytest.fixture(scope="module")
def frame():
    from conftest import FIXTURES
    return build_model(FIXTURES / "valid" / "frame_kip_in.json")


def test_nodes_and_dof_numbering(frame):
    assert DOFS == ("ux", "uy", "uz", "rx", "ry", "rz")
    assert list(frame.nodes) == ["N1", "N2", "N3", "N4", "N5", "N6", "N7"]
    assert frame.nodes["N3"].index == 2
    np.testing.assert_allclose(frame.nodes["N3"].xyz, (240 * IN, 0, 144 * IN))
    assert frame.n_dof == 42
    assert frame.dof("N1", "ux") == 0 and frame.dof("N3", "rz") == 17 and frame.dof("N2", 4) == 10
    assert frame.vertical_axis == "Z" and frame.units.to_dict() == {"force": "kip", "length": "in"}


def test_supports(frame):
    s = frame.supports
    assert s["N1"].fixed == (True,) * 6 and s["N1"].spring == (0.0,) * 6
    assert s["N4"].fixed == (False, True, True, False, False, False)         # ux overridden by a spring
    assert s["N4"].spring[0] == pytest.approx(250 * KIP / IN)                # kip/in -> N/mm
    assert s["N4"].spring[4] == pytest.approx(60000 * KIP * IN)              # kip*in/rad -> N*mm/rad
    assert s["N6"].fixed == (True, True, True, False, False, False)
    assert s["N7"].fixed == (True, True, True, False, False, False)
    assert set(s) == {"N1", "N4", "N6", "N7"}


def test_roller_restrains_the_vertical_translation(fixtures):
    m = build_model(fixtures / "valid" / "beam.json")                        # vertical axis Y
    assert m.supports["b"].fixed == (False, True, True, False, False, False)  # uy (roller) + uz (given)
    assert m.supports["a"].fixed == (True, True, True, True, False, False)


def test_materials_and_sections(frame):
    st, ot = frame.materials["steel"], frame.materials["other"]
    assert st.E == pytest.approx(29000 * KIP / IN**2) and st.G == pytest.approx(11200 * KIP / IN**2)
    assert st.Fy == pytest.approx(50 * KIP / IN**2) and st.Fu == pytest.approx(65 * KIP / IN**2)
    assert st.unit_weight == pytest.approx(0.000284 * KIP / IN**3)
    assert ot.G == pytest.approx(ot.E / 2.5) and ot.Fy is None and ot.unit_weight == 0.0
    col = frame.sections["col"]
    assert col.kind == "I" and col.family == "I"
    assert col.dims["d"] == pytest.approx(14 * IN) and col.dims["r"] == pytest.approx(0.5 * IN)
    a_in2 = 2 * 10 * 0.7 + (14 - 1.4) * 0.44 + 4 * 0.25 * (1 - math.pi / 4)
    assert col.A == pytest.approx(a_in2 * IN**2)
    assert frame.sections["beam"].J == pytest.approx(1.25 * IN**4)           # override
    gen = frame.sections["general"]
    assert (gen.A, gen.Iz, gen.Asy) == pytest.approx((10 * IN**2, 200 * IN**4, 4 * IN**2))
    assert gen.Asz is None and gen.Zz == pytest.approx(30 * IN**3) and gen.Zy is None
    assert frame.sections["round"].J == pytest.approx(2 * frame.sections["round"].Iy)


def test_members(frame):
    m = frame.members
    assert m["C1"].type == "frame" and not m["C1"].has_releases and not m["C1"].has_offsets
    np.testing.assert_allclose(m["C1"].R, [(0, 0, 1), (1, 0, 0), (0, 1, 0)], atol=1e-15)
    np.testing.assert_allclose(m["C2"].R, [(0, 0, 1), (0, 1, 0), (-1, 0, 0)], atol=1e-15)   # roll 90
    b = m["B1a"]
    np.testing.assert_allclose(b.start, (7 * IN, 0, 144 * IN))               # local offset along x
    np.testing.assert_allclose(b.offset_i, (7 * IN, 0, 0), atol=1e-12)
    assert b.length == pytest.approx(113 * IN)
    assert b.k_end_i == (RIGID, RIGID, RIGID, RIGID, RIGID, 0.0) and b.k_end_j == (RIGID,) * 6
    c = m["B1b"]
    assert c.length == pytest.approx(113 * IN) and c.has_offsets and c.has_releases
    assert c.k_end_j[0] == pytest.approx(900 * KIP / IN) and c.k_end_j[5] == pytest.approx(50000 * KIP * IN)
    assert c.k_end_j[1:5] == (RIGID,) * 4
    assert m["BR1"].type == "tensionOnly" and m["BR2"].type == "truss"
    assert m["BR1"].k_end_i == (RIGID, RIGID, RIGID, 0.0, 0.0, 0.0)
    np.testing.assert_allclose(m["BR2"].R[1] @ m["BR2"].R[0], 0, atol=1e-15)
    assert m["BR2"].R[1][0] > 0                                              # y on the reference vector's side
    assert m["BR1"].length == pytest.approx(math.hypot(120, 144) * IN)


def test_constraints(frame):
    c = frame.constraints["floor"]
    assert (c.type, c.master, c.slaves, c.normal, c.dofs) == ("rigidDiaphragm", "N5", ("N2", "N3"), "Z",
                                                              ("ux", "uy", "rz"))


def test_loads_are_normalised_to_local_axes_and_absolute_positions(frame):
    d, l, w, s = (frame.load_cases[k].loads for k in "DLWS")
    sw = d[0]
    np.testing.assert_allclose(sw.direction, (0, 0, -1))
    assert sw.factor == 1.0 and sw.members == tuple(frame.members)
    np.testing.assert_allclose(d[1].values, (0, 0, -10 * KIP, 0, 20 * KIP * IN, 0))
    # global Z load on a horizontal beam whose local y is +Z
    np.testing.assert_allclose(d[2].w1, (0, -0.1 * KIP / IN, 0), atol=1e-12)
    assert (d[2].x1, d[2].x2) == (0.0, pytest.approx(113 * IN))
    np.testing.assert_allclose(d[3].w2, (0, -0.05 * KIP / IN, 0), atol=1e-12)
    assert (d[3].x1, d[3].x2) == (pytest.approx(0.2 * 113 * IN), pytest.approx(0.9 * 113 * IN))
    np.testing.assert_allclose(l[0].vector, (0, -5 * KIP, 0))
    assert l[0].x == pytest.approx(30 * IN) and l[0].kind == "force"
    # moment about global Y on a member with local z = -Y
    np.testing.assert_allclose(l[1].vector, (0, 0, -100 * KIP * IN), atol=1e-6)
    assert l[1].kind == "moment" and l[1].x == pytest.approx(0.5 * 113 * IN)
    # projected load on the inclined brace: member (0,-120,144)/187.45; projection factor = 120/187.45
    length = math.hypot(120, 144)
    wz = -0.02 * (120 / length) * KIP / IN                                   # per unit member length, global Z
    np.testing.assert_allclose(l[2].w1, (wz * 144 / length, wz * 120 / length, 0), atol=1e-12)
    assert (l[2].x1, l[2].x2) == (pytest.approx(12 * IN), pytest.approx(100 * IN))
    np.testing.assert_allclose(w[0].values, (8 * KIP, 0, 0, 0, 0, 0))
    np.testing.assert_allclose(w[1].w1, (0, 0.03 * KIP / IN, 0))             # local direction: unchanged
    assert s[0].node == "N1" and s[0].values == {2: pytest.approx(-0.5 * IN), 4: 0.001}


def test_combinations_and_analysis_options(frame):
    c = frame.combinations
    assert c["LRFD1"].analysis_type == "firstOrder" and c["LRFD1"].purpose == "strength"
    assert not c["LRFD1"].stiffness_reduction and c["LRFD1"].notional is None
    assert c["LRFD2"].analysis_type == "secondOrder" and c["LRFD2"].stiffness_reduction
    assert c["LRFD2"].notional.coefficient == 0.002
    np.testing.assert_allclose(c["LRFD2"].notional.direction, (1, 0, 0))
    assert c["LRFD3"].factors == {"D": 1.2, "W": 1.0, "L": 0.5, "S": 1.0}
    assert c["SLS"].purpose == "serviceability"
    a = frame.analysis
    assert (a.type, a.shear_deformation, a.result_stations, a.max_iterations) == ("firstOrder", True, 21, 50)
    assert (a.buckling_combinations, a.buckling_modes) == (("LRFD1",), 2)
    assert (a.reduction_factor, a.tau_b, a.reduction_members) == (0.8, "computed", ("C1", "B1a", "B1b"))
    assert (a.alpha, a.convergence_tolerance) == (1.0, 1e-6)


# -------------------------------------------------------------------- design variables
def test_design_members_and_provenance(frame):
    dm = frame.design_members
    assert list(dm) == ["B1", "BR1", "C1", "C2", "BR2"]            # explicit first, then automatic
    b1 = dm["B1"]
    assert b1.members == ("B1a", "B1b") and b1.length == pytest.approx(226 * IN)
    assert (b1.section, b1.material) == ("beam", "steel")
    assert (b1.Lx, b1.source("Lx")) == (pytest.approx(226 * IN), "member")
    assert (b1.Cb, b1.source("Cb")) == (1.0, "group:beams")
    assert (b1.method, b1.source("method")) == ("LRFD", "global")
    assert (b1.Kx, b1.source("Kx")) == (1.0, "default")
    assert (b1.role, b1.source("role")) == ("beam", "default")
    # bracing given in absolute positions by the group -> fractions of the design member length
    assert b1.braceTopFlange == pytest.approx((116.5 / 226,)) and b1.source("braceTopFlange") == "group:beams"
    assert b1.braceMinor == pytest.approx((58.25 / 226, 116.5 / 226))
    assert b1.braceBottomFlange == () and b1.source("braceBottomFlange") == "default"
    assert (b1.Ly, b1.source("Ly")) == (None, "bracing")           # defined by the brace points
    assert (b1.Lb, b1.source("Lb")) == (None, "bracing")
    assert (b1.Lz, b1.source("Lz")) == (pytest.approx(226 * IN), "default")
    assert b1.values["Cb"].path == "/design/groups/beams/Cb"
    assert b1.values["Kx"].path == "/design/members/B1/Kx"
    lim = b1.deflectionLimits
    assert [(x.span_ratio, x.absolute, x.direction, x.combinations) for x in lim] == [
        (360, None, "resultant", ("SLS",)), (None, pytest.approx(0.75 * IN), "y", ("SLS",))]

    c1 = dm["C1"]
    assert (c1.Kx, c1.source("Kx")) == (2.0, "group:columns")
    assert (c1.role, c1.source("role")) == ("column", "group:columns")
    assert (c1.Ly, c1.source("Ly")) == (pytest.approx(72 * IN), "member")
    assert c1.braceTorsion == (0.5,) and (c1.Lz, c1.source("Lz")) == (None, "bracing")
    assert c1.An == pytest.approx(frame.sections["col"].A) and c1.source("An") == "default"

    br1 = dm["BR1"]
    assert (br1.An, br1.U, br1.source("U")) == (pytest.approx(4.5 * IN**2), 0.85, "member")
    assert (br1.role, br1.source("role")) == ("brace", "default")
    assert br1.maxSlendernessTension == 300.0 and br1.maxSlendernessCompression == 200.0
    assert dm["C2"].role == "column" and dm["BR2"].role == "brace"
    with pytest.raises(AttributeError):
        b1.nonexistent


def base():
    return {
        "formatVersion": "1.0.0", "units": {"force": "kN", "length": "m"}, "conventions": {"verticalAxis": "Z"},
        "nodes": {"A": {"x": 0, "y": 0}, "B": {"x": 4, "y": 0}, "C": {"x": 10, "y": 0}},
        "supports": {"A": {"type": "fixed"}},
        "materials": {"S": {"E": 2e8, "G": 8e7, "Fy": 355e3, "Fu": 470e3}},
        "sections": {"R": {"type": "rectangle", "h": 0.3, "b": 0.2}},
        "members": {"M1": {"i": "A", "j": "B", "section": "R", "material": "S"},
                    "M2": {"i": "B", "j": "C", "section": "R", "material": "S"}},
        "loadCases": {"D": {"loads": [{"type": "nodal", "node": "B", "Fz": -1}]}},
        "combinations": {"C1": {"factors": {"D": 1.0}}},
    }


def test_resolution_order_member_group_global_default():
    d = base()
    d["design"] = {
        "defaults": {"Kx": 1.1, "Ky": 1.1, "Kz": 1.1},
        "groups": {"lo": {"members": ["M1", "M2"], "Kx": 1.2, "Ky": 1.2, "priority": 1},
                   "hi": {"members": ["M1"], "Kx": 1.3, "priority": 5},
                   "same": {"members": ["M1"], "Ky": 1.2, "priority": 1}},      # equal value: no conflict
        "members": {"M1": {"Kx": 1.4}},
    }
    m = build_model(d)
    m1, m2 = m.design_members["M1"], m.design_members["M2"]
    assert (m1.Kx, m1.source("Kx")) == (1.4, "member")
    assert (m2.Kx, m2.source("Kx")) == (1.2, "group:lo")
    assert (m1.Ky, m1.source("Ky")) == (1.2, "group:lo")
    assert (m1.Kz, m1.source("Kz")) == (1.1, "global")
    assert (m1.U, m1.source("U")) == (1.0, "default")
    del d["design"]["members"]
    assert build_model(d).design_members["M1"].source("Kx") == "group:hi"        # higher priority wins


def test_explicit_length_beats_lower_level_bracing_and_vice_versa():
    d = base()
    d["design"] = {"defaults": {"bracing": {"major": [0.5], "topFlange": [0.25, 0.75]}, "Lz": 2.0},
                   "members": {"M1": {"Lx": 3.0, "bracing": {"torsion": [0.5]}}}}
    m = build_model(d)
    m1, m2 = m.design_members["M1"], m.design_members["M2"]
    assert (m1.Lx, m1.source("Lx"), m1.braceMajor) == (3000.0, "member", ())     # member Lx beats global bracing
    assert (m2.Lx, m2.source("Lx"), m2.braceMajor) == (None, "bracing", (0.5,))
    assert (m1.Lz, m1.source("Lz"), m1.braceTorsion) == (None, "bracing", (0.5,))  # member bracing beats global Lz
    assert (m2.Lz, m2.source("Lz")) == (2000.0, "global")
    assert m1.braceTopFlange == (0.25, 0.75) and m1.Lb is None
    assert (m1.Ly, m1.source("Ly")) == (4000.0, "default") and m2.Ly == 6000.0


def test_auto_members_can_be_switched_off():
    d = base()
    d["design"] = {"autoMembers": False, "members": {"girder": {"members": ["M1", "M2"]}}}
    m = build_model(d)
    assert list(m.design_members) == ["girder"] and m.design_members["girder"].length == pytest.approx(10000.0)
    d["design"] = {"autoMembers": False}
    assert build_model(d).design_members == {}


def test_design_variables_cannot_influence_the_analysis_model():
    plain = base()
    loaded = copy.deepcopy(plain)
    loaded["design"] = {
        "defaults": {"method": "ASD", "Kx": 3.0, "Cb": 2.0, "Lb": 1.0, "An": 0.01, "U": 0.5, "role": "brace",
                     "sideswayInhibited": True, "Cmx": 0.6, "maxSlendernessCompression": 50,
                     "deflectionLimits": [{"spanRatio": 100}]},
        "groups": {"g": {"members": ["M1"], "Ly": 0.5, "bracing": {"torsion": [0.3]}}},
        "members": {"girder": {"members": ["M1", "M2"], "Kz": 0.5}}, "autoMembers": False,
    }
    loaded["design"]["groups"]["g"]["members"] = ["girder"]
    a, b = build_model(plain), build_model(loaded)
    assert_same(a, b, ignore=("design_members", "defaults_used", "warnings"))
    assert a.design_members.keys() != b.design_members.keys()


# ------------------------------------------------------------------- defaults-used report
def test_defaults_used_report(fixtures):
    m = build_model(fixtures / "valid" / "beam.json")
    rep = {x.path: x for x in m.defaults_used}
    assert len(rep) == len(m.defaults_used)                        # one entry per path
    for x in m.defaults_used:
        assert x.path.startswith("/") and x.reason and set(x.to_dict()) == {"path", "value", "reason"}
    # conventions, analysis options
    assert rep["/conventions/internalForces"].value == "positive-face"
    assert rep["/analysis/type"].value == "firstOrder" and rep["/analysis/resultStations"].value == 11
    assert rep["/analysis/shearDeformation"].value is True
    # derived material value, in FILE units (kN/m^2)
    assert rep["/materials/S355/G"].value == pytest.approx(2.1e8 / 2.6)
    assert rep["/members/m1/type"].value == "frame"
    assert rep["/combinations/ULS/analysisType"].value == "firstOrder"
    assert rep["/combinations/ULS/purpose"].value == "strength"
    assert rep["/loadCases/G/loads/0/positionType"].value == "relative"
    # every design variable of the automatic design member, lengths in FILE units (m)
    assert rep["/design/members/m1/Lb"].value == pytest.approx(6.0)
    assert rep["/design/members/m1/An"].value == pytest.approx(0.045)
    assert rep["/design/members/m1/Kx"].value == 1.0 and "C3" in rep["/design/members/m1/Kx"].reason
    assert rep["/design/members/m1/Cb"].value == "computed" and "F1-1" in rep["/design/members/m1/Cb"].reason
    assert rep["/design/members/m1/maxSlendernessCompression"].value == 200.0
    assert rep["/design/members/m1/method"].value == "LRFD"
    assert rep["/design/members/m1/deflectionLimits"].value == []
    # values that ARE given are not reported
    assert "/materials/S355/E" not in rep and "/conventions/verticalAxis" not in rep


def test_given_values_are_not_reported_as_defaults(frame):
    paths = {x.path for x in frame.defaults_used}
    assert "/analysis/resultStations" not in paths and "/analysis/buckling/modes" not in paths
    assert "/materials/steel/G" not in paths and "/materials/other/G" in paths
    assert "/design/members/B1/Lx" not in paths and "/design/members/B1/Cb" not in paths
    assert "/design/members/B1/Kx" in paths and "/design/members/C1/Kx" not in paths
    assert "/analysis/stiffnessReduction/factor" in paths and "/analysis/stiffnessReduction/tauB" in paths
    assert "/combinations/LRFD2/notionalLoad/coefficient" in paths
    assert "/constraints/floor/normal" in paths and "/loadCases/D/loads/0/direction" in paths
    assert "/sections/general/Asz" in paths and "/sections/general/Asy" not in paths
    assert "/sections/beam/r" in paths and "/sections/col/r" not in paths


# ------------------------------------------------------------------------ shapes database
CSV = ("Type,AISC_Manual_Label,A,d,bf,tw,tf,Ix,Zx,Sx,Iy,Zy,Sy,J,Cw,Ht,B,tdes,OD\n"
       "W,TESTW10,10.0,10.0,8.0,0.3,0.5,170.0,38.0,34.0,43.0,16.0,10.75,0.6,1000.0,,,,\n"
       "HSS,TESTHSS6X4,4.0,,,,,20.0,8.0,6.7,10.0,6.0,5.0,22.0,,6.0,4.0,0.25,\n"
       "PIPE,TESTPIPE,3.0,,,,,7.0,3.0,2.5,7.0,3.0,2.5,14.0,,,,0.2,5.0\n"
       "L,TESTL,2.0,,,,,3.0,1.5,1.0,3.0,1.5,1.0,0.05,,,,,\n")       # made-up numbers, not AISC data


def test_standard_shape_without_database_is_refused_with_hint(fixtures):
    with pytest.raises(ModelInvalid) as e:
        build_model(fixtures / "invalid" / "standard_shape.json")
    (p,) = e.value.errors
    assert p.code == "section.shapes-db" and "BENCHFEA_SHAPES_DB" in p.hint
    with pytest.raises(shapes.ShapesDatabaseError):
        shapes.lookup("W14X90", "mm")


@pytest.mark.parametrize("kind", ["csv", "xlsx"])
def test_user_supplied_database_is_read_and_converted(fixtures, tmp_path, monkeypatch, kind):
    path = tmp_path / f"shapes.{kind}"
    if kind == "csv":
        path.write_text(CSV)
    else:
        openpyxl = pytest.importorskip("openpyxl")
        wb = openpyxl.Workbook()
        wb.active.title = "Readme"
        wb.active.append(["some notes"])
        ws = wb.create_sheet("Database")
        for line in CSV.strip().splitlines():
            ws.append([(float(c) if c.replace(".", "").isdigit() else (c or None)) for c in line.split(",")])
        wb.save(path)
    monkeypatch.setenv("BENCHFEA_SHAPES_DB", str(path))
    m = build_model(fixtures / "invalid" / "standard_shape.json")      # valid once the database exists
    s = m.sections["W"]
    assert (s.kind, s.shape, s.family) == ("standard", "TESTW10", "I")
    assert s.A == pytest.approx(10 * IN**2)
    assert s.Iz == pytest.approx(170 * IN**4) and s.Iy == pytest.approx(43 * IN**4)   # AISC x-x -> local z
    assert s.Zz == pytest.approx(38 * IN**3) and s.Sy == pytest.approx(10.75 * IN**3)
    assert s.J == pytest.approx(0.6 * IN**4) and s.Cw == pytest.approx(1000 * IN**6)
    assert s.dims == pytest.approx({"d": 254.0, "bf": 203.2, "tf": 12.7, "tw": 7.62})
    assert s.Asy == pytest.approx(10 * 0.3 * IN**2) and s.Asz == pytest.approx(5 / 6 * 2 * 8 * 0.5 * IN**2)
    hss = shapes.lookup("testhss6x4", "in")                           # case-insensitive, native units
    assert (hss.family, hss.dims, hss.Asy, hss.Asz) == ("box", {"H": 6.0, "B": 4.0, "t": 0.25}, 3.0, 2.0)
    pipe = shapes.lookup("TESTPIPE", "in")
    assert (pipe.family, pipe.dims, pipe.Asy) == ("pipe", {"D": 5.0, "t": 0.2}, 1.5)
    angle = shapes.lookup("TESTL", "in")
    assert (angle.family, angle.Asy, angle.Asz, angle.dims) == ("L", None, None, {})
    with pytest.raises(shapes.ShapesDatabaseError) as e:
        shapes.lookup("W99X999", "in")
    assert "not in the shapes database" in str(e.value)


def test_database_with_unexpected_headers_is_refused(tmp_path, monkeypatch):
    path = tmp_path / "shapes.csv"
    path.write_text("name,area\nX,1\n")
    monkeypatch.setenv("BENCHFEA_SHAPES_DB", str(path))
    with pytest.raises(shapes.ShapesDatabaseError) as e:
        shapes.lookup("X", "in")
    assert "BENCHFEA_SHAPES_COLUMNS" in e.value.hint


# --------------------------------------------------------------------------- results stubs
def test_results_document_is_written_in_file_units(fixtures):
    from benchfea.results import CombinationResult, MemberResult, Results
    m = build_model(fixtures / "valid" / "beam.json")                 # kN, m
    n = 3
    forces = np.zeros((n, 6))
    forces[:, 0] = 5000.0                                             # N
    forces[:, 5] = [0.0, 2.0e6, 0.0]                                  # N*mm
    end = np.concatenate([-forces[0], forces[-1]])
    disp = np.zeros((2, 6))
    disp[1, 1] = -12.0                                                # mm
    disp[1, 5] = 0.01                                                 # rad
    react = np.zeros((2, 6))
    react[0, 1] = 7500.0
    mr = MemberResult(end, np.array([0.0, 3000.0, 6000.0]), forces, np.zeros((n, 3)))
    res = Results(m.name, {"ULS": CombinationResult("ULS", "firstOrder", disp, react, {"m1": mr})},
                  defaults_used=m.defaults_used)
    res.compute_envelopes()
    doc = res.to_document(m)
    assert (doc["format"], doc["formatVersion"], doc["units"]) == ("benchfea-results", "1.0.0",
                                                                   {"force": "kN", "length": "m"})
    c = doc["combinations"]["ULS"]
    assert c["displacements"]["b"]["uy"] == pytest.approx(-0.012) and c["displacements"]["b"]["rz"] == 0.01
    assert c["reactions"] == {"a": {"Fx": 0, "Fy": 7.5, "Fz": 0, "Mx": 0, "My": 0, "Mz": 0},
                              "b": dict.fromkeys(("Fx", "Fy", "Fz", "Mx", "My", "Mz"), 0.0)}
    st = c["members"]["m1"]["stations"]
    assert st["x"] == [0.0, 3.0, 6.0] and st["N"] == [5.0, 5.0, 5.0] and st["Mz"] == pytest.approx([0, 2.0, 0])
    assert c["members"]["m1"]["endForces"]["i"]["Fx"] == -5.0 and c["members"]["m1"]["endForces"]["j"]["Fx"] == 5.0
    assert doc["envelopes"]["m1"]["max"]["Mz"] == pytest.approx([0, 2.0, 0])
    assert doc["envelopes"]["m1"]["maxCombination"]["Mz"] == ["ULS"] * 3
    assert len(doc["defaultsUsed"]) == len(m.defaults_used)
    import json
    json.dumps(doc)                                                   # plain JSON data
