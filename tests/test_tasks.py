"""Tasks of the benchmark (benchfea.bench.tasks): cases with a verified reference, the four task types,
what is public and what is hidden.  docs/BENCHMARK.md."""
import copy
import json

import pytest
from bench_helpers import CELLS, case, full_coverage_registry, leaves, task

from benchfea.bench import curveballs, generate
from benchfea.bench import tasks as T
from benchfea.bench.families._common import Draw
from benchfea.model import build_model
from benchfea.solve import analyze


# ------------------------------------------------------------------------------------------- cases
def test_derived_seed_chain_is_deterministic_and_starts_at_the_seed():
    assert T.derived_seed("portal_frame", 7, 2, 0) == 7
    chain = [T.derived_seed("portal_frame", 7, 2, k) for k in range(1, 4)]
    assert chain == [T.derived_seed("portal_frame", 7, 2, k) for k in range(1, 4)]
    assert len(set(chain)) == 3 and all(0 <= s < 2 ** 31 for s in chain)
    assert T.derived_seed("portal_frame", 8, 2, 1) != chain[0] and T.derived_seed("stair_frame", 7, 2, 1) != chain[0]


def test_case_is_verified_and_the_record_is_the_difference():
    c = case("change")
    assert c["status"] == "verified" and c["effectClass"] == "changes-analysis"
    v = c["verification"]
    assert v["status"] == "verified" and v["kind"] == "results" and set("CDE") <= set(v["tiersPassed"])
    assert not v["tiersNotPassed"]
    record = curveballs.ChangeRecord.from_dict(c["record"])
    assert record.apply(c["base"]["document"]) == c["perturbed"]["document"]
    assert record.revert(c["perturbed"]["document"]) == c["base"]["document"]
    # the served model does not carry the generator's name (it states the seed)
    assert "seed" not in c["perturbed"]["document"]["name"]
    # the reference is the results document the pipeline verified: the engine reproduces it
    model = build_model(copy.deepcopy(c["perturbed"]["document"]))
    assert analyze(model).to_document(model)["combinations"] == v["reference"]["combinations"]


def test_case_is_deterministic():
    a = T.build_case(*CELLS["mechanism"])                       # no cache: built again from nothing
    b = case("mechanism")
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_mechanism_case_reference_is_the_verified_diagnosis():
    c = case("mechanism")
    v = c["verification"]
    assert c["effectClass"] == "creates-mechanism" and v["kind"] == "diagnosis" and v["status"] == "verified"
    assert v["tiersPassed"] == ["N", "K", "C", "E"]
    assert v["reference"]["diagnosis"]["nMechanisms"] >= 1
    types = T.task_types_for(c)
    assert types["predict"] and types["design"] and types["diagnose"] is None and types["repair"] is None
    with pytest.raises(T.TaskUnavailable):
        T.make_task(c, "predict", "tool-allowed")


def test_design_only_case_is_quarantined_not_served():
    """The reference of a design-only curve ball is the design results; with the present tier-F coverage
    they are not verified, so every attempt of the chain is skipped and the skip is logged."""
    c = case("design_only")
    assert c["status"] == "unavailable" and "perturbed" not in c
    assert len(c["skips"]) == T.MAX_ATTEMPTS
    assert [s["attempt"] for s in c["skips"]] == list(range(T.MAX_ATTEMPTS))
    assert c["skips"][0]["generatorSeed"] == CELLS["design_only"][2]
    assert c["skips"][1]["generatorSeed"] == T.derived_seed("continuous_beam", 0, 2, 1)
    covered = [s for s in c["skips"] if s["stage"] == "verify"]
    assert covered and all(s["category"] in ("not-covered", "disagreement", "not-representable") for s in covered)
    assert any(s["category"] == "not-covered" and "tier F" in s["reason"] for s in covered)
    with pytest.raises(T.TaskUnavailable):
        T.make_task(c, "diagnose", "tool-allowed")


def test_unknown_names_are_errors():
    with pytest.raises(KeyError):
        T.build_case("continuous_beam", "no_such_curve_ball", 0, 1)
    with pytest.raises(ValueError):
        T.make_task(case("change"), "guess", "tool-allowed")
    with pytest.raises(ValueError):
        T.make_task(case("change"), "predict", "with-solver")


# ----------------------------------------------------------------------------------------- predict
HIDDEN_KEYS = ("changeRecord", "curveBall", "generatorSeed", "reference", "intendedModel", "skips", "tiersPassed")


def _keys(value):
    if isinstance(value, dict):
        for k, v in value.items():
            yield k
            yield from _keys(v)
    elif isinstance(value, list):
        for v in value:
            yield from _keys(v)


@pytest.mark.parametrize("name,task_type", [("change", "predict"), ("change", "diagnose"), ("change", "repair"),
                                            ("mechanism", "diagnose"), ("mechanism", "repair")])
def test_public_task_holds_nothing_of_the_hidden_part(name, task_type):
    public, answers = task(name, task_type)
    keys = set(_keys(public))
    assert not keys & set(HIDDEN_KEYS)
    text = json.dumps(public)
    cb = CELLS[name][1]
    assert cb not in text and answers["changeRecord"]["target"] not in text
    assert public["id"] == answers["id"] and cb not in public["id"]
    assert public["verification"]["status"] == "verified" and "tiers" not in json.dumps(public["verification"]["status"])
    assert "kind" not in public["verification"]                   # would tell a mechanism from a stable model
    assert answers["verification"]["status"] == "verified" and answers["verification"]["tiersPassed"]
    assert all("tolerances" in t for t in answers["verification"]["tiers"].values())
    assert public["model"] == case(name)["perturbed"]["document"]


def test_predict_quantities_are_those_of_the_given_model():
    public, answers = task("change", "predict")
    ref = answers["reference"]["values"]
    assert 3 <= len(public["quantities"]) <= 6 and [q["id"] for q in public["quantities"]] == list(ref)
    model = build_model(copy.deepcopy(public["model"]))
    results = analyze(model).to_document(model)
    base = case("change")["base"]["verification"]["reference"]
    moved = 0
    for q in public["quantities"]:
        assert "value" not in q and "floor" not in q and q["id"] in public["statement"]
        value = T.extract(q, results)
        assert value == pytest.approx(ref[q["id"]]["value"], rel=1e-12)
        assert ref[q["id"]]["floor"] > 0 and abs(value) >= 50 * ref[q["id"]]["floor"] * 0.999   # significant: >= 5 % of scale
        moved += abs(value - T.extract(q, base)) >= 0.05 * max(abs(value), abs(T.extract(q, base)))
    assert moved == len(public["quantities"])                    # the curve ball matters for what is asked
    for mode in T.MODES:
        assert T.make_task(case("change"), "predict", mode)[0]["quantities"] == public["quantities"]
    assert T.make_task(case("change"), "predict", "no-solver")[0]["id"] != public["id"]


def test_select_quantities_is_seeded_and_one_per_entity():
    ref = case("change")["verification"]["reference"]
    a = T.select_quantities(ref, None, Draw(1))
    assert a == T.select_quantities(ref, None, Draw(1))
    entities = [(q["kind"], q.get("node") or q.get("member")) for q in a]
    assert len(set(entities)) == len(entities)
    assert {q["kind"] for q in a} <= set(T.QUANTITY_KINDS)


def test_extract_and_units():
    ref = case("change")["verification"]["reference"]
    cid = next(iter(ref["combinations"]))
    mid = next(iter(ref["combinations"][cid]["members"]))
    q = {"kind": "memberMaxMoment", "combination": cid, "member": mid, "component": "My"}
    assert T.extract(q, ref) == max(abs(x) for x in ref["combinations"][cid]["members"][mid]["stations"]["My"])
    with pytest.raises(KeyError):
        T.extract({"kind": "displacement", "combination": cid, "node": "nope", "component": "ux"}, ref)
    kn_m, kip_in = {"force": "kN", "length": "m"}, {"force": "kip", "length": "in"}
    assert T.convert_value(1.0, "length", kn_m, {"force": "kN", "length": "mm"}) == pytest.approx(1000.0)
    assert T.convert_value(1.0, "force", kip_in, kn_m) == pytest.approx(4.4482216152605)
    assert T.convert_value(1.0, "moment", kip_in, kn_m) == pytest.approx(4.4482216152605 * 0.0254)
    assert T.convert_value(0.37, "dimensionless", kip_in, kn_m) == 0.37


# ---------------------------------------------------------------------------------------- describe
@pytest.mark.parametrize("family,seed,level", [("continuous_beam", 1, 2), ("portal_frame", 1, 2),
                                               ("stair_frame", 1, 1), ("braced_frame", 1, 2)])
def test_description_states_every_value_of_the_model(family, seed, level):
    doc = generate(family, seed, level).document
    text = T.describe(doc)
    name = doc.pop("name")
    # (a load type such as "memberDistributed" is written in words: "member distributed")
    missing = [v for v in leaves(doc) if not isinstance(v, bool) and T._num(v) not in text
               and T._words(str(v)) not in text]
    assert not missing
    for section in ("nodes", "members", "sections", "materials", "loadCases", "combinations", "supports"):
        assert all(str(k) in text for k in doc[section])
    assert name not in text


def test_diagnose_task_gives_the_intended_design_in_words_not_the_file():
    public, answers = task("change", "diagnose")
    c = case("change")
    assert public["intendedDesign"] == T.describe(c["base"]["document"]) and public["intendedDesign"] in public["statement"]
    assert json.dumps(c["base"]["document"]) not in json.dumps(public)
    ref = answers["reference"]
    assert ref["category"] == "releases" and ref["direction"] == "added" and ref["mechanism"] is None
    assert ref["entities"] == c["record"]["entities"] and ref["entities"]["members"]
    for word in list(T.CATEGORY_TEXT) + list(T.DIRECTIONS):
        assert word in public["statement"]
    _, m = task("mechanism", "diagnose")
    assert m["reference"]["mechanism"] == case("mechanism")["verification"]["reference"]["diagnosis"]


def test_change_direction():
    d = T.change_direction
    assert d([{"op": "add", "path": "/a", "after": 1}]) == "added"
    assert d([{"op": "remove", "path": "/a", "before": 1}, {"op": "remove", "path": "/b", "before": 2}]) == "removed"
    assert d([{"op": "replace", "path": "/a", "before": 2.0, "after": -2.0}]) == "reversed"
    assert d([{"op": "replace", "path": "/a", "before": 2.0, "after": 5}]) == "increased"
    assert d([{"op": "replace", "path": "/a", "before": -2.0, "after": 1}]) == "decreased"
    assert d([{"op": "replace", "path": "/a", "before": "B1", "after": "B2"}]) == "changed"
    assert d([{"op": "replace", "path": "/a", "before": True, "after": False}]) == "changed"
    assert d([{"op": "add", "path": "/a", "after": 1}, {"op": "replace", "path": "/b", "before": 1, "after": 2}]) == "changed"


# ------------------------------------------------------------------------------------------ repair
def test_repair_task_is_judged_on_the_verified_intended_design():
    public, answers = task("change", "repair")
    c = case("change")
    assert c["base"]["verification"]["status"] == "verified"
    ref = answers["reference"]
    assert ref["intendedModel"] == c["base"]["document"] and answers["baseVerification"]["status"] == "verified"
    assert ref["designCriterion"] is False and "not a criterion" in public["criteria"]["design"]
    base, given = c["base"]["verification"]["reference"], c["verification"]["reference"]
    for q in ref["quantities"]:
        v = T.extract(q, base)
        assert v == ref["values"][q["id"]]["value"]
        # the model as given is beyond the zero-credit error on every quantity asked
        assert abs(T.extract(q, given) - v) > T.REPAIR_TOLERANCE["zero"] * abs(v)
    assert json.dumps(ref["values"]) not in json.dumps(public)
    assert all(str(r["value"]) not in public["statement"] for r in ref["values"].values())


# ------------------------------------------------------------------------------------------ design
def test_design_task_is_not_served_without_verified_design_results():
    c = case("change")
    assert c["designVerification"]["status"] == "quarantined"
    assert c["designVerification"]["quarantine"]["category"] == "not-covered"
    with pytest.raises(T.TaskUnavailable) as e:
        T.make_task(c, "design", "tool-allowed")
    assert "tier F" in e.value.reason


def test_design_task_and_reference_search_with_the_coverage_test_device(monkeypatch):
    """The machinery of design tasks, exercised with a tier-F registry that counts everything as covered
    (a test device; see bench_helpers.FullCoverage).  Checks the reference search, not any design value."""
    import benchfea.verify.tier_f as tier_f
    monkeypatch.setattr(tier_f, "coverage_registry", full_coverage_registry)
    c = T.build_case(*CELLS["change"])
    assert c["designVerification"]["status"] == "verified"
    public, answers = T.make_task(c, "design", "tool-allowed")
    cat, ref = public["catalogue"], answers["reference"]
    assert set(cat["groups"]) == {m["section"] for m in public["model"]["members"].values()}
    assert all(len(o) > 5 for o in cat["options"].values())
    assert "not a proven optimum" in ref["search"] and "not a proven optimum" in public["criteria"]["weight"]
    doc = T.design_document(public["model"], cat, ref["sections"])
    assert T.design_feasible(doc)[0] and T.model_weight(doc) == ref["weight"]
    given_ok = T.design_feasible(public["model"])[0]
    assert ref["givenWeight"]["feasible"] is given_ok
    assert ("every design check passes" if given_ok else "at least one design check fails") in public["statement"]
    # one size lighter in any single group is not feasible (the search stops only there), unless a
    # lighter section further down was infeasible too -- so check the first lighter one of each group
    for sid, members in cat["groups"].items():
        current = ref["sections"][members[0]]
        i = cat["options"][sid].index(current)
        for lighter in cat["options"][sid][:i]:
            trial = dict(ref["sections"], **{m: lighter for m in members})
            assert not T.design_feasible(T.design_document(public["model"], cat, trial))[0]
    # the search is deterministic
    assert T.reference_design(public["model"], cat)["sections"] == ref["sections"]
    with pytest.raises(ValueError):
        T.design_document(public["model"], cat, {"nope": "x"})
    with pytest.raises(ValueError):
        T.design_document(public["model"], cat, {next(iter(public["model"]["members"])): "I-not-in-the-catalogue"})
    assert answers["referenceDesignVerification"]["status"] == "verified"
    # a predict task then also asks for a design ratio
    p, pa = T.make_task(c, "predict", "tool-allowed")
    assert p["quantities"][-1]["kind"] == "designRatio" and 0 < pa["reference"]["values"][p["quantities"][-1]["id"]]["value"]
