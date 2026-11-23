"""Independent review: tension-only / compression-only members against the reviewer's active-set loop
(tableau) and against PyNite's own tension/compression-only iteration."""
import pytest

import oracles
from harness import compare, compare_foreign, run_benchfea, worst
from models import ONE_SIDED

SLACK = {"T01_two_storey_x_braced_frame": {"x1b", "x2a"},
         "T02_space_tower_with_cross_bracing_y_up_kip_in": {"xb10", "xb11", "xa12", "xa13", "xb20", "xa21", "xa22", "xa23"},
         "T03_compression_only_props_and_tension_only_ties_with_settlement": {"prop1", "prop2", "tie2"}}


@pytest.mark.parametrize("second_order", [False, True])
@pytest.mark.parametrize("name", sorted(ONE_SIDED))
def test_one_sided_members_against_the_tableau_reference(name, second_order):
    raw, cid = ONE_SIDED[name]()
    if second_order:
        raw["combinations"][cid]["analysisType"] = "secondOrder"
        raw.setdefault("analysis", {})["convergenceTolerance"] = 1.0e-11
    _, res, doc = run_benchfea(raw)
    c = doc["combinations"][cid]
    assert {m for m, v in c["members"].items() if not v["active"]} == SLACK[name]
    out = compare(raw, cid, second_order=second_order, doc=doc)     # also asserts the same active set
    assert worst(out) <= 1.0e-9, {k: v for k, v in out.items() if not k.startswith("_")}


@pytest.mark.parametrize("name", sorted(ONE_SIDED))
def test_one_sided_members_against_pynite(name):
    if not oracles.have("pynite"):
        pytest.skip("PyNite is not installed")
    raw, cid = ONE_SIDED[name]()
    _, res, doc = run_benchfea(raw)
    f = oracles.run_pynite(raw, cid)
    assert f["slack"] == SLACK[name]
    for q, v in compare_foreign(raw, doc["combinations"][cid], f).items():
        assert v <= 1.0e-9, f"{name} vs PyNite: {q} differs by {v:.2e}"


def test_a_slack_member_still_hands_its_own_loads_to_its_ends():
    """Self weight, a point load and a line load on members that go slack: the engine carries them to
    the end nodes as a simply supported beam (THEORY 13); the reference keeps 1e-13 of the axial
    stiffness of a slack member, which gives the same statics."""
    raw, cid = ONE_SIDED["T01_two_storey_x_braced_frame"]()
    raw["loadCases"]["W"]["loads"] += [
        {"type": "selfWeight", "factor": 30.0},
        {"type": "memberPoint", "member": "x1b", "direction": "Z", "value": -8, "position": 0.3},
        {"type": "memberDistributed", "member": "x2a", "direction": "X", "w1": 2.0, "w2": 5.0}]
    _, res, doc = run_benchfea(raw)
    assert {m for m, v in doc["combinations"][cid]["members"].items() if not v["active"]} == {"x1b", "x2a"}
    out = compare(raw, cid, doc=doc)
    assert worst(out) <= 1.0e-9, {k: v for k, v in out.items() if not k.startswith("_")}
    slack = res.combinations[cid].members["x1b"]
    assert abs(slack.end_forces[0] + slack.end_forces[6]) > 0.0       # it carries its own axial load component
