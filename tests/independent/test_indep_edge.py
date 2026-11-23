"""Independent review: edge cases, each compared with the tableau reference (1e-9 unless stated)."""
import pytest

from harness import case, compare, member, nodes, run_benchfea, worst
from models import LINEAR, SECOND, _so, start


def agree(raw, cid, tol=1.0e-9, **kw):
    _, res, doc = run_benchfea(raw)
    out = compare(raw, cid, doc=doc, **kw)
    assert worst(out) <= tol, {k: v for k, v in out.items() if not k.startswith("_")}
    return res


@pytest.mark.parametrize("k", [1e-12, 1e-3, 1e3, 1e9, 1e15, 1e25])
def test_member_end_springs_from_almost_free_to_almost_rigid(k):
    raw, cid = LINEAR["L03_partial_releases_axis_aligned_kn_mm"]()
    raw["members"]["b1"]["releases"] = {"i": {"rz": {"spring": k}, "ux": {"spring": k}}, "j": {"uy": {"spring": k}}}
    agree(raw, cid)


def test_support_springs_over_33_decades_are_in_the_bank():
    raw, cid = LINEAR["L11_support_springs_stiff_and_soft"]()
    ks = [v["spring"] for s in raw["supports"].values() for v in s.get("dofs", {}).values() if isinstance(v, dict)]
    assert min(ks) <= 1e-9 and max(ks) >= 1e18
    res = agree(raw, cid)
    c = res.combinations[cid]
    # a spring reaction is -k u (CONVENTIONS 1): the 1e15 N/mm spring under node A carries load, the 1e-9 one none
    assert abs(c.reactions[0, 2]) > 1.0 and abs(c.reactions[0, 4]) < 1e-6


def test_supports_and_imposed_displacements_on_the_free_dofs_of_diaphragm_nodes():
    raw, cid = LINEAR["L08_diaphragm_z_up_unequal_columns"]()
    raw["supports"]["C"] = {"dofs": {"uz": "fixed", "rx": {"spring": 2000.0}}}        # slave: out-of-plane DOFs are its own
    raw["loadCases"]["H"]["loads"] += [{"type": "supportDisplacement", "node": "M", "uz": -0.01, "rx": 0.002},
                                       {"type": "supportDisplacement", "node": "C", "uz": 0.004}]
    agree(raw, cid)


def test_point_loads_at_the_member_ends_and_abutting_line_loads():
    raw, cid = LINEAR["L02_space_frame_y_up_kip_in"]()
    raw["loadCases"]["G"]["loads"] += [
        {"type": "memberPoint", "member": "bm1", "direction": "Y", "value": -5, "position": 0.0},
        {"type": "memberPoint", "member": "bm1", "direction": "z", "value": 3, "position": 1.0},
        {"type": "memberPoint", "member": "bm2", "direction": "x", "kind": "moment", "value": 30, "position": 1.0},
        {"type": "memberDistributed", "member": "bm2", "direction": "Y", "w1": -0.1, "w2": 0.0, "start": 0.5, "end": 1.0},
        {"type": "memberDistributed", "member": "bm2", "direction": "Y", "w1": -0.1, "w2": -0.3, "start": 0.0, "end": 0.5},
        {"type": "memberPoint", "member": "bm2", "direction": "Y", "value": -4, "position": 0.5}]
    res = agree(raw, cid)
    m = res.combinations[cid].members["bm2"]
    # CONVENTIONS 4: the station at 0 is minus the end forces at i, the station at L the end forces at j
    scale = abs(m.end_forces).max()
    assert (m.forces[0] == -m.end_forces[:6]).all()
    assert abs(m.forces[-1] - m.end_forces[6:]).max() <= 1e-12 * scale      # closure of the integration: round-off
    assert (m.x[:-1] == m.x[1:]).sum() == 1                                  # two stations at the interior point load


def test_column_pointing_down_reference_vector_with_roll_and_nearly_vertical_member():
    raw, cid = LINEAR["L02_space_frame_y_up_kip_in"]()
    raw["members"]["col0"].update({"i": "A", "j": "a"})
    raw["members"]["col2"]["orientation"] = {"referenceVector": [1, 0.2, 0.3], "roll": 77}
    agree(raw, cid)
    d = start()
    nodes(d, A=(0, 0, 0), B=(1e-5, 0, 4), C=(3, 0, 4))
    d["supports"] = {"A": {"type": "fixed"}}
    member(d, "c", "A", "B")
    member(d, "b", "B", "C")
    case(d, "L", {"type": "nodal", "node": "C", "Fy": 5, "Fz": -3},
         {"type": "memberDistributed", "member": "c", "direction": "y", "w1": 2})
    res = agree(d, "L")
    assert any("near-vertical" in w for w in res.warnings)


def test_offsets_that_overlap_reverse_the_flexible_length_consistently():
    d = start()
    nodes(d, A=(0, 0, 0), B=(4, 0, 0))
    d["supports"] = {"A": {"type": "fixed"}}
    member(d, "b", "A", "B", offsets={"i": [3, 0, 0], "j": [-2, 0, 0]})
    case(d, "L", {"type": "nodal", "node": "B", "Fz": -3, "Fy": 2, "Mx": 1},
         {"type": "memberDistributed", "member": "b", "direction": "Z", "w1": -4})
    agree(d, "L")


def test_unconnected_node_is_reported_and_left_at_zero():
    raw, cid = LINEAR["L05_offsets_off_axis_skew"]()
    nodes(raw, LOOSE=(1, 1, 1))
    res = agree(raw, cid)
    assert any("not connected to any member" in w for w in res.warnings)
    assert not res.combinations[cid].displacements[-1].any()


def test_second_order_with_shear_deformation_and_partial_releases():
    raw, cid = SECOND["S06_truss_members_with_transverse_load_and_partial_releases"]()
    raw["analysis"]["shearDeformation"] = True
    for s in raw["sections"].values():
        s["Asy"], s["Asz"] = 0.08 * s["A"], 0.1 * s["A"]
    raw["loadCases"]["U"]["loads"] = [ld for ld in raw["loadCases"]["U"]["loads"] if ld["type"] == "nodal"]
    agree(raw, cid, second_order=True)


def test_second_order_with_axial_and_torsional_releases():
    raw, cid = SECOND["S03_portal_with_leaning_column"]()
    raw["members"]["g1"]["releases"] = {"i": {"ux": {"spring": 5000.0}, "rx": True}}
    agree(raw, cid, second_order=True)


def test_second_order_with_settlements_a_support_spring_and_an_axial_point_load():
    raw, cid = LINEAR["L12_prescribed_support_displacements"]()
    raw["supports"]["D"] = {"type": "pinned", "dofs": {"rx": {"spring": 5000.0}}}
    raw["loadCases"]["S"]["loads"] += [{"type": "nodal", "node": n, "Fy": -120.0} for n in "BCF"] + \
        [{"type": "memberPoint", "member": "c1", "direction": "x", "value": -80, "position": 0.5}]
    _so(raw)
    # the axial force of c1 jumps at the point load and that of the inclined girders varies: the
    # reference steps N piecewise (16 steps, O(h^2)); measured 2e-8
    agree(raw, cid, tol=1.0e-6, second_order=True, nsub=16)
