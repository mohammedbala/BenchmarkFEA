"""Model generators (``benchfea.bench``): determinism, validity, stability, proportions, measured
difficulty, feature flags, the tier-C expectation.

Generating a model means sizing it with the engine (analysis + member checks, several times), so a
model costs from a fraction of a second (a level-1 beam) to minutes (the largest 3D frames).  The
default run uses the cases of ``CASES`` -- every family, the cheapest seeds; the whole matrix (every
family x every level x ``FULL_SEEDS``) runs with ``BENCHFEA_GENERATOR_FULL=1``.

``GOLDEN`` pins the SHA-256 of the JSON text of some cases: generation must give the same bytes on
every run and platform.  A deliberate change of a family means bumping
``benchfea.bench.GENERATOR_VERSION`` and updating the hashes.
"""
import functools
import hashlib
import json
import os

import numpy as np
import pytest

from benchfea import bench
from benchfea import design as design_module
from benchfea.bench import (BUCKLING_MIN, FAMILIES, GenerationFailed, Draw, difficulty, feature_flags, generate,
                            level_of, levels, section_definition, stable_seed, static_properties, tier_c_expected)
from benchfea.bench.families._common import CATALOGUE, UNIT_SYSTEMS, catalogue_size
from benchfea.crosscheck import cross_check
from benchfea.crosscheck.runner import available
from benchfea.io import loads
from benchfea.model import build_model
from benchfea.solve import analyze
from benchfea.validate import validate

from solve_reference import make

ALL_FAMILIES = ("continuous_beam", "portal_frame", "multistorey_2d", "multistorey_3d", "braced_frame", "pipe_rack",
                "stair_frame", "diaphragm_frame")

# (family, level, seed): the default sample -- see the module docstring
CASES = [
    ("continuous_beam", 1, 0), ("continuous_beam", 1, 1), ("continuous_beam", 2, 0), ("continuous_beam", 2, 1),
    ("continuous_beam", 3, 0), ("continuous_beam", 3, 1),
    ("portal_frame", 1, 0), ("portal_frame", 2, 0), ("portal_frame", 3, 1),
    ("stair_frame", 1, 0), ("stair_frame", 1, 1), ("stair_frame", 2, 0), ("stair_frame", 3, 1), ("stair_frame", 4, 0),
    ("braced_frame", 2, 2), ("braced_frame", 3, 1),
    ("multistorey_2d", 2, 0), ("multistorey_2d", 3, 0),
    ("multistorey_3d", 3, 0), ("multistorey_3d", 4, 0),
    ("pipe_rack", 3, 0),
    ("diaphragm_frame", 3, 2), ("diaphragm_frame", 4, 2), ("diaphragm_frame", 5, 0),
]
FULL_SEEDS = (0, 1)
if os.environ.get("BENCHFEA_GENERATOR_FULL"):
    CASES = sorted(set(CASES) | {(f, lv, s) for f in ALL_FAMILIES for lv in levels(f) for s in FULL_SEEDS})

GOLDEN = {
    ('continuous_beam', 1, 0): "bc10d0ef896b8117a0c8c2209d30d7e1ccc093f967d9a7c6968a017de31b142a",
    ('continuous_beam', 2, 1): "5fad5ccc03631bf41d7f03310599b996e1f1c4df9e3a923dabbcf06796ec9d81",
    ('continuous_beam', 3, 0): "eed0f1ed5cd1ecef95cf6dcecc86c0d82e13dda0f2fafd5acd7e8dcf66c501f8",
    ('portal_frame', 1, 0): "59208e968c7583dc3cfa9c1a9f4a93f0f88bc01f34b0b72b9aeebbafda1ee75c",
    ('portal_frame', 3, 1): "1833d0ac6aaffe90e988c004a4a4e201cedcc04f9e58c6f62ddeab7516fe70fd",
    ('stair_frame', 1, 1): "2e745f9eec35be5d4b011166b4627276f388914599a24ac50fbca3862a01fe4d",
    ('stair_frame', 3, 1): "02741dba3b64ec7119eb283a701f1adfafe55ca5b4c607347064c96c0602ac25",
    ('braced_frame', 3, 1): "4e69fd42c23fc45880c032f27031fe6a825fe72aa5237d7580d3e30e39d05824",
    ('multistorey_2d', 3, 0): "441de670612a62ed9f13b6e515b9a941f683d71217da5dc34bd7183e3d753d75",
    ('multistorey_3d', 3, 0): "92c52b2d06664a5f00f5096320e60272f8ce85682306c0c428177a3d4b542292",
    ('pipe_rack', 3, 0): "e69e5a2fef1e6a2365746a0239961bdd0e89173d9670f4f71a88c7a548283302",
    ('diaphragm_frame', 3, 2): "af308084d6a77d3544c4d988a39ad863d97fbe57e501c7e22e3c5b45fb39bd79",
    ('diaphragm_frame', 5, 0): "3e40c74c0ab5e76a8bc88cf0ec05d64308fb3e8046613f8af40db89d2933d46f",
}


@functools.lru_cache(maxsize=None)
def gen(family, level, seed):
    return generate(family, seed, level)


@functools.lru_cache(maxsize=None)
def solved(family, level, seed):
    g = gen(family, level, seed)
    model = build_model(g.document)
    results = analyze(model)
    return model, results, design_module.check(model, results)


def case_id(case):
    return "{}-L{}-s{}".format(*case)


# ------------------------------------------------------------------------------------ random draws
def test_stable_seed_is_a_fixed_function_of_its_text():
    # SHA-256 of "benchfea.bench|1|portal_frame|7|2|0", first 16 bytes, big endian (not Python's hash())
    text = "benchfea.bench|1|portal_frame|7|2|0"
    assert stable_seed("1", "portal_frame", 7, 2, 0) == int.from_bytes(hashlib.sha256(text.encode()).digest()[:16], "big")
    assert stable_seed("1", "portal_frame", 7, 2, 0) == 185008142103257040451120527542217913913
    assert stable_seed("1", "portal_frame", 7, 2, 1) != stable_seed("1", "portal_frame", 7, 2, 0)


def test_draws_are_pinned():
    """The stream of PCG64 and the arithmetic on it: pinned values (a change here changes every model)."""
    d = Draw(stable_seed("x", 1))
    assert isinstance(d.generator, np.random.Generator) and isinstance(d.generator.bit_generator, np.random.PCG64)
    assert [d.integer(0, 99) for _ in range(6)] == [8, 12, 82, 57, 5, 84]
    assert d.u() == 0.724715475196404
    assert d.step(4, 10, 0.5) == 8.5
    assert d.choice("abc") == "a"


def test_draw_ranges():
    d = Draw(3)
    xs = [d.integer(2, 5) for _ in range(400)]
    assert set(xs) == {2, 3, 4, 5}
    assert all(0.0 <= d.u() < 1.0 for _ in range(100))
    assert all(v in (1.0, 1.5, 2.0) for v in (d.step(1.0, 2.0, 0.5) for _ in range(50)))
    with pytest.raises(ValueError):
        d.integer(3, 2)


# ------------------------------------------------------------------------------------ catalogue
@pytest.mark.parametrize("unit_system", sorted(UNIT_SYSTEMS))
def test_catalogue_is_generic_valid_and_ordered(unit_system):
    """Every catalogue entry is a valid parametric section with a generic name; within a kind the
    plastic modulus / area grows with the index (so 'one size up' is stronger)."""
    system = "metric" if UNIT_SYSTEMS[unit_system]["metric"] else "imperial"
    for kind, entries in CATALOGUE[system].items():
        sections, order = {}, []
        for index in range(len(entries)):
            name, definition, depth = section_definition(kind, index, unit_system)
            assert name.split("-")[0] in ("I", "Box", "Pipe") and depth > 0
            assert not any(name.upper().startswith(p) for p in ("W", "HSS", "HP", "IPE", "HE", "UB", "UC"))
            sections[name] = definition
            order.append(name)
        assert len(sections) == len(entries) == catalogue_size(kind, unit_system)       # no duplicate names
        doc = make({"a": (0, 0, 0), "b": (1000, 0, 0)}, {"m": {"i": "a", "j": "b"}}, {"a": "fixed"},
                   {"P": [{"type": "nodal", "node": "b", "Fz": -1.0}]})
        doc["units"] = dict(UNIT_SYSTEMS[unit_system]["units"])
        doc["sections"] = sections
        doc["materials"] = {"steel": dict(UNIT_SYSTEMS[unit_system]["material"])}
        doc["members"]["m"].update(section=order[0], material="steel")
        assert [p for p in validate(doc) if p.severity == "error"] == []
        model = build_model(doc)
        props = [model.sections[n].props for n in order]
        strength = [p.A if kind == "pipe" else p.Zz for p in props]
        assert all(b > a for a, b in zip(strength, strength[1:])), (kind, strength)


def test_material_is_the_same_steel_in_every_unit_system():
    ref = None
    for name, us in UNIT_SYSTEMS.items():
        doc = make({"a": (0, 0, 0), "b": (1, 0, 0)}, {"m": {"i": "a", "j": "b"}}, {"a": "fixed"},
                   {"P": [{"type": "nodal", "node": "b", "Fz": -1.0}]})
        doc["units"] = dict(us["units"])
        doc["materials"] = {"steel": dict(us["material"])}
        doc["sections"] = {"s": section_definition("box", 5, name)[1]}
        doc["members"]["m"].update(section="s", material="steel")
        m = build_model(doc).materials["steel"]
        values = np.array([m.E, m.G, m.Fy, m.Fu, m.unit_weight])
        if name == "kip-in":                       # nominal imperial steel: within 1 % of the metric one
            assert values == pytest.approx(ref, rel=0.012)
        elif ref is None:
            ref = values
        else:
            assert values == pytest.approx(ref, rel=1e-12)


# ------------------------------------------------------------------------------------ API
def test_families_and_levels():
    assert tuple(bench.families()) == ALL_FAMILIES
    assert set(lv for f in ALL_FAMILIES for lv in levels(f)) == {1, 2, 3, 4, 5}
    for f in ALL_FAMILIES:
        assert FAMILIES[f].NAME == f and FAMILIES[f].FEATURES and (FAMILIES[f].__doc__ or "").strip()


def test_bad_arguments_are_refused():
    with pytest.raises(ValueError, match="unknown family"):
        generate("bridge", 0, 1)
    with pytest.raises(ValueError, match="difficulty levels"):
        generate("continuous_beam", 0, 5)
    with pytest.raises(ValueError, match="no feature flag"):
        generate("continuous_beam", 0, 1, features={"secondOrder": True})
    with pytest.raises(ValueError, match="unit system"):
        generate("continuous_beam", 0, 1, units="N-m")
    with pytest.raises(ValueError, match="vertical_axis"):
        generate("continuous_beam", 0, 1, vertical_axis="X")


def test_a_level_that_a_forced_feature_makes_unreachable_fails_with_the_reasons():
    # a portal frame with second order, base springs and semi-rigid eaves scores above level 1 in every draw
    with pytest.raises(GenerationFailed) as e:
        generate("portal_frame", 0, 1, features={"secondOrder": True, "baseSprings": True, "semiRigid": True},
                 max_attempts=6)
    assert len(e.value.reasons) == 6 and all("draw gives level" in r for r in e.value.reasons)


# ------------------------------------------------------------------------------------ determinism
@pytest.mark.parametrize("case", [c for c in CASES if c[0] in ("continuous_beam", "stair_frame") and c[1] <= 2],
                         ids=case_id)
def test_same_seed_gives_byte_identical_json(case):
    family, level, seed = case
    first = gen(*case)
    again = generate(family, seed, level)                    # a second, independent generation
    assert again.json() == first.json()
    assert json.dumps(again.metadata, sort_keys=True) == json.dumps(first.metadata, sort_keys=True)
    assert again.json().endswith("\n") and loads(again.json(), "json").data == first.document


@pytest.mark.parametrize("case", sorted(GOLDEN), ids=case_id)
def test_output_is_pinned(case):
    """Same bytes as when the hash was recorded (another run, possibly another platform)."""
    if case not in CASES:
        pytest.skip("not in the default sample")
    assert hashlib.sha256(gen(*case).json().encode("utf-8")).hexdigest() == GOLDEN[case]


def test_different_seeds_and_levels_give_different_models():
    texts = {gen(*c).json() for c in CASES if c[0] == "continuous_beam"}
    assert len(texts) == len([c for c in CASES if c[0] == "continuous_beam"])


def test_coordinates_and_loads_are_round_numbers():
    """Readable files: no binary noise -- every number of the document has at most 6 significant digits."""
    def numbers(x):
        if isinstance(x, dict):
            for v in x.values():
                yield from numbers(v)
        elif isinstance(x, list):
            for v in x:
                yield from numbers(v)
        elif isinstance(x, (int, float)) and not isinstance(x, bool):
            yield x
    for case in CASES[:8]:
        for v in numbers(gen(*case).document):
            assert float(f"{v:.6g}") == v, (case, v)


# ------------------------------------------------------------------------------------ validity
@pytest.mark.parametrize("case", CASES, ids=case_id)
def test_generated_model_is_valid_stable_and_sized(case):
    family, level, seed = case
    g = gen(*case)
    doc, meta = g.document, g.metadata
    json.dumps(meta)                                           # plain data
    assert [p for p in validate(doc)] == [], "schema + semantic validation, warnings included"
    model, results, design = solved(*case)                     # an analysis independent of the generator's
    assert set(results.combinations) == set(doc["combinations"])
    assert all(c.invariants.passed and c.converged for c in results.combinations.values())
    # design: no failure, nothing above 1.0, most members between 0.3 and 1.1
    assert design.status in ("pass", "conditional"), design.status
    ratios = [m.governing.ratio for m in design.members.values() if m.governing is not None]
    assert len(ratios) == len(design.members) == meta["design"]["nDesignMembers"]
    assert max(ratios) <= 1.0 and max(ratios) == pytest.approx(meta["design"]["ratioMax"], abs=1e-3)
    inside = sum(1 for r in ratios if 0.3 <= r <= 1.1)
    assert inside >= 0.5 * len(ratios)
    assert max(ratios) >= 0.5, "the governing member is not oversized"
    # stability of the gravity combination
    st = meta["stability"]
    if family == "continuous_beam":
        assert st is None
    else:
        assert st["combination"] == "ULS1" and st["minimum"] == BUCKLING_MIN
        if "buckling" in doc["analysis"]:
            assert float(results.buckling["ULS1"].factors[0]) >= BUCKLING_MIN
            assert st["bucklingFactor"] == pytest.approx(float(results.buckling["ULS1"].factors[0]), abs=1e-3)
        else:
            assert st["bucklingFactor"] is None and "second-order analysis at 2 x" in st["basis"]
    # the design block gives curve balls something to act on
    assert doc["design"]["groups"] and all("role" in grp for grp in doc["design"]["groups"].values())
    assert any("deflectionLimits" in grp for grp in doc["design"]["groups"].values())
    assert any("bracing" in grp for grp in doc["design"]["groups"].values())
    # envelope of the cross-check: no shear deformation, no notional loads, no stiffness reduction
    assert doc["analysis"]["shearDeformation"] is False
    assert not any("notionalLoad" in c or c.get("stiffnessReduction") for c in doc["combinations"].values())
    assert meta["units"] in UNIT_SYSTEMS and doc["conventions"]["verticalAxis"] == meta["verticalAxis"]
    for name in doc["sections"]:
        assert name.split("-")[0] in ("I", "Box", "Pipe")


@pytest.mark.parametrize("case", CASES, ids=case_id)
def test_measured_difficulty_is_the_requested_level(case):
    family, level, seed = case
    g = gen(*case)
    _, results, _ = solved(*case)
    measured = difficulty(g.document, results)                 # measured again, from the document alone
    assert measured["level"] == level == g.metadata["difficulty"]
    assert measured == g.metadata["measured"]
    assert measured["score"] == sum(measured["components"].values()) and level_of(measured["score"]) == level
    assert g.metadata["features"] == feature_flags(g.document)


def test_stability_probe_agrees_with_a_buckling_analysis():
    """A model accepted on the second-order probe (no buckling analysis in the file) has a buckling
    factor of ULS1 above the minimum when the buckling analysis is run after all."""
    case = next(c for c in CASES if c[0] == "portal_frame" and "buckling" not in gen(*c).document["analysis"])
    doc = json.loads(gen(*case).json())
    doc["analysis"]["buckling"] = {"combinations": ["ULS1"], "modes": 1}
    doc.pop("design")
    assert float(analyze(build_model(doc)).buckling["ULS1"].factors[0]) >= BUCKLING_MIN


def test_units_and_vertical_axes_all_appear_and_describe_the_same_structure():
    """Forcing the unit system or the vertical axis does not change the structure drawn: the same
    seed gives the same physical model, so the largest deflection (as a length) and the design
    ratios agree -- exactly between the two axis conventions and between kN-m and kN-mm, and within
    the rounding of the imperial grid sizes otherwise (not compared)."""
    seen = set()
    out = {}
    for units in sorted(UNIT_SYSTEMS):
        for axis in ("Y", "Z"):
            g = generate("continuous_beam", 4, 2, units=units, vertical_axis=axis)
            assert g.document["units"] == UNIT_SYSTEMS[units]["units"] and g.document["conventions"]["verticalAxis"] == axis
            seen.add((g.metadata["units"], g.metadata["verticalAxis"]))
            model = build_model(g.document)
            res = analyze(model)
            sls = next(c for c in res.combinations if c.startswith("SLS"))
            out[units, axis] = (float(np.abs(res.combinations[sls].displacements[:, :3]).max()),     # mm
                                g.metadata["design"]["ratioMax"], g.metadata["sizing"]["sections"])
    assert len(seen) == 6
    for units in sorted(UNIT_SYSTEMS):
        assert out[units, "Y"][0] == pytest.approx(out[units, "Z"][0], rel=1e-9)
        assert out[units, "Y"][1:] == out[units, "Z"][1:]
    assert out["kN-m", "Z"][0] == pytest.approx(out["kN-mm", "Z"][0], rel=1e-9)
    drawn = {(gen(*c).metadata["units"], gen(*c).metadata["verticalAxis"]) for c in CASES}
    assert {u for u, _ in drawn} == set(UNIT_SYSTEMS) and {a for _, a in drawn} == {"Y", "Z"}


def test_three_dimensional_families_are_measured_as_3d_and_plane_ones_as_planar():
    for case in CASES:
        m = gen(*case).metadata["measured"]
        if case[0] in ("multistorey_3d", "pipe_rack", "diaphragm_frame") or (case[0] == "stair_frame" and case[1] >= 3):
            assert m["threeDimensional"]
        else:
            assert m["planar"], case


# ------------------------------------------------------------------------------------ feature flags
def test_forced_features_appear_in_the_document_and_the_metadata():
    g = generate("continuous_beam", 2, 3, features={"hinge": True, "settlement": True})      # 6-7 spans
    f = g.metadata["features"]
    assert f["releasesFull"] and f["prescribedDisplacements"] and g.metadata["parameters"]["hinge"]
    assert g.metadata["requestedFeatures"] == {"hinge": True, "settlement": True}
    for seed in range(3):
        g = generate("continuous_beam", seed, 2, features={"hinge": False, "springSupport": False, "partialFixity": False})
        f = g.metadata["features"]
        assert not f["releasesFull"] and not f["supportSprings"]
    g = generate("portal_frame", 1, 2, features={"secondOrder": True})
    assert g.metadata["features"]["secondOrder"] and g.metadata["measured"]["amplification"] >= 1.0
    assert all(c.get("analysisType") == "secondOrder" for c in g.document["combinations"].values()
               if c.get("purpose", "strength") == "strength")


def test_family_specific_features():
    for case in CASES:
        g = gen(*case)
        f, doc = g.metadata["features"], g.document
        types = {m.get("type", "frame") for m in doc["members"].values()}
        if case[0] == "diaphragm_frame":
            assert f["rigidDiaphragm"] and all(c["type"] == "rigidDiaphragm" for c in doc["constraints"].values())
            masters = {c["master"] for c in doc["constraints"].values()}
            assert not any(m["i"] in masters or m["j"] in masters for m in doc["members"].values())
            lateral = [l for c in doc["loadCases"].values() for l in c["loads"] if l["type"] == "nodal"]
            assert lateral and all(l["node"] in masters for l in lateral)
        if case[0] in ("braced_frame", "pipe_rack"):
            assert "truss" in types and f["trussMembers"]
        if case[0] == "stair_frame":
            model = solved(*case)[0]
            inclined = [m for m in model.members.values()
                        if 0.05 < abs(m.R[0]["XYZ".index(model.vertical_axis)]) < 0.95]
            assert inclined and f["loadProjected"]
        if f["secondOrder"]:
            assert not f["rigidOffsets"], "second order with rigid offsets is outside the cross-check envelope"


def test_opt_in_features_outside_the_cross_check_envelope():
    """Tension-only bracing and a buckling analysis in the file are produced only on request, and the
    metadata then says the model cannot be cross-checked."""
    g = generate("braced_frame", 2, 2, features={"tensionOnly": True})
    doc, meta = g.document, g.metadata
    assert meta["features"]["tensionOnly"] and meta["parameters"]["bracing"] == "X"
    assert meta["tierC_expected"]["crossCheckable"] is False
    assert meta["tierC_expected"]["opensees"] == meta["tierC_expected"]["pynite"] == "unsupported"
    assert validate(doc) == []
    results = analyze(build_model(doc))
    assert all(c.invariants.passed for c in results.combinations.values())
    # every combination carries a lateral load, so one diagonal of every panel works and one is slack
    n_braces = sum(1 for m in doc["members"].values() if m.get("type") == "tensionOnly")
    for cid, c in results.combinations.items():
        slack = [m for m, r in c.members.items() if not r.active]
        assert len(slack) == n_braces // 2, cid
    assert difficulty(doc, results)["level"] == 2
    # a tension-only model is first order and carries its stability check as a buckling analysis in the file
    assert not meta["features"]["secondOrder"] and doc["analysis"]["buckling"]["combinations"] == ["ULS1"]
    assert meta["stability"]["bucklingFactor"] >= BUCKLING_MIN
    assert meta["design"]["fractionBetween0.3And1.1"] >= 0.5

    g = generate("portal_frame", 0, 1, features={"buckling": True})
    assert g.document["analysis"]["buckling"] == {"combinations": ["ULS1"], "modes": 2}
    assert g.metadata["tierC_expected"]["crossCheckable"] is False and g.metadata["tierC_expected"]["buckling"] == "unsupported"
    assert g.metadata["stability"]["basis"] == "linear buckling analysis of the file"
    assert g.metadata["stability"]["bucklingFactor"] >= BUCKLING_MIN
    assert g.metadata["bucklingFactors"]["ULS1"] == g.metadata["stability"]["bucklingFactor"]
    # the flag changes nothing else: same structure as without it
    plain = generate("portal_frame", 0, 1)
    a, b = json.loads(g.json()), json.loads(plain.json())
    a["analysis"].pop("buckling")
    assert a == b


# ------------------------------------------------------------------------------------ difficulty()
BEAM = {"a": (0, 0, 0), "b": (4000, 0, 0)}
LOAD = {"P": [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -1.0}]}


@pytest.mark.parametrize("supports, degree, free", [
    ({"a": {"dofs": {"ux": "fixed", "uy": "fixed", "uz": "fixed", "rx": "fixed"}}, "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}}, 0, 3),
    ({"a": "fixed", "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}}, 1, 2),
    ({"a": "fixed", "b": "fixed"}, 3, 0),
    ({"a": "fixed"}, 0, 3),
])
def test_static_indeterminacy_of_single_beams(supports, degree, free):
    """Plane beams (counted in their plane): simply supported 0, propped cantilever 1, fixed-fixed 3,
    cantilever 0; free DOFs = the in-plane DOFs no support fixes."""
    p = static_properties(build_model(make(BEAM, {"m": {"i": "a", "j": "b"}}, supports, LOAD)))
    assert p["planar"] and p["staticIndeterminacy"] == degree and p["freeDofs"] == free


def test_static_indeterminacy_in_3d_and_with_a_hinge():
    load_y = {"P": [{"type": "memberDistributed", "member": "m", "direction": "Y", "w1": 1.0},
                    {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -1.0}]}
    p = static_properties(build_model(make(BEAM, {"m": {"i": "a", "j": "b"}}, {"a": "fixed", "b": "fixed"}, load_y)))
    assert not p["planar"] and p["staticIndeterminacy"] == 6 and p["freeDofs"] == 0
    nodes = {"a": (0, 0, 0), "h": (2000, 0, 0), "b": (4000, 0, 0)}
    members = {"m": {"i": "a", "j": "h"}, "n": {"i": "h", "j": "b", "releases": {"i": {"rz": True}}}}
    p = static_properties(build_model(make(nodes, members, {"a": "fixed", "b": "fixed"}, LOAD)))
    assert p["planar"] and p["staticIndeterminacy"] == 2          # fixed-fixed (3) minus one hinge


def test_difficulty_components_and_level_limits():
    assert [level_of(s) for s in (0, 2, 3, 5, 6, 8, 9, 11, 12, 16)] == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    doc = make(BEAM, {"m": {"i": "a", "j": "b"}}, {"a": "fixed", "b": {"dofs": {"uy": "fixed", "uz": "fixed"}}}, LOAD)
    d = difficulty(doc)
    assert d["components"] == {"size": 0, "indeterminacy": 1, "dimension": 0, "features": 0, "secondOrder": 0,
                               "combinations": 0}
    assert d["score"] == 1 and d["level"] == 1 and d["amplification"] is None and d["featureKinds"] == []


def test_amplification_decides_the_second_order_points():
    """A cantilever column with a lateral load H and an axial load P (derived here from the
    beam-column equation EI v'' + P v = H (L - x)): tip deflection = H (tan kL - kL) / (P k),
    k = sqrt(P / EI); first order H L^3 / (3 EI).  The measured amplification is their ratio."""
    from solve_reference import E, SEC
    L, EI = 3000.0, E * SEC["Iz"]              # a column has its local y along global X: bending about local z
    for fraction, points in ((0.02, 1), (0.30, 2)):               # P as a fraction of pi^2 EI / (4 L^2)
        P = fraction * np.pi ** 2 * EI / (4 * L ** 2)
        # (the top is held out of plane so that the weak axis does not buckle first)
        doc = make({"a": (0, 0, 0), "b": (0, 0, L)}, {"m": {"i": "a", "j": "b"}}, {"a": "fixed", "b": {"dofs": {"uy": "fixed"}}},
                   {"P": [{"type": "nodal", "node": "b", "Fz": -P, "Fx": 1000.0}]},
                   {"c": {"factors": {"P": 1.0}, "analysisType": "secondOrder"}}, {"shearDeformation": False})
        d = difficulty(doc)
        k = np.sqrt(P / EI)
        exact = (np.tan(k * L) - k * L) / (P * k) / (L ** 3 / (3 * EI))
        assert d["amplification"] == pytest.approx(exact, rel=1e-3)   # rounded to 4 decimals in the output
        assert d["components"]["secondOrder"] == points
        pre = difficulty(doc, analyse=False)
        assert pre["level"] is None and pre["amplification"] is None and pre["levelRange"] == [1, 1]


def test_worked_examples_have_the_documented_levels():
    """docs/DIFFICULTY.md quotes these."""
    from pathlib import Path
    from benchfea.io import load
    ex = Path(__file__).resolve().parents[1] / "examples"
    got = {}
    for name in ("01_simple_beam", "02_continuous_beam", "03_portal_frame", "04_diaphragm_frame", "05_braced_frame"):
        d = difficulty(load(ex / f"{name}.json").data)
        got[name] = (d["score"], d["level"])
    assert got == {"01_simple_beam": (1, 1), "02_continuous_beam": (4, 2), "03_portal_frame": (4, 2),
                   "04_diaphragm_frame": (7, 3), "05_braced_frame": (7, 3)}


# ------------------------------------------------------------------------------------ tier C
def test_tier_c_expectation_follows_the_features():
    plain = gen("continuous_beam", 1, 0)
    assert plain.metadata["tierC_expected"]["crossCheckable"] is True
    assert plain.metadata["tierC_expected"]["opensees"] == plain.metadata["tierC_expected"]["pynite"] == "exact"
    doc = json.loads(gen("braced_frame", 2, 2).json())
    for m in doc["members"].values():
        if m.get("type") == "truss":
            m["type"] = "tensionOnly"
    t = tier_c_expected(doc)
    assert t["crossCheckable"] is False and t["opensees"] == t["pynite"] == "unsupported"
    doc = json.loads(plain.json())
    doc["analysis"]["shearDeformation"] = True
    t = tier_c_expected(doc)
    assert t["pynite"] == "unsupported" and t["opensees"] == "exact" and t["crossCheckable"] is False
    for case in CASES:
        g = gen(*case)
        t = g.metadata["tierC_expected"]
        # default draws stay inside the envelope of the cross-check: the two opt-in features are absent
        assert not g.metadata["features"]["tensionOnly"] and not g.metadata["features"]["buckling"], case
        assert t["crossCheckable"] is True and t["buckling"] == "none", case
        so = [c for c, v in g.document["combinations"].items() if v.get("analysisType") == "secondOrder"]
        assert all(t["byCombination"][c]["opensees"] == "approximate" for c in so)


@pytest.mark.skipif(not (available("opensees") and available("pynite")), reason="needs OpenSeesPy and PyNite")
@pytest.mark.parametrize("case", [("continuous_beam", 2, 0), ("stair_frame", 1, 0), ("portal_frame", 2, 0)], ids=case_id)
def test_expected_cross_checkable_models_are_cross_checked(case):
    """The expectation is borne out on small models: both foreign solvers are eligible and agree with
    the engine, in every combination, at the default tolerances of ``cross_check``."""
    g = gen(*case)
    assert g.metadata["tierC_expected"]["crossCheckable"]
    model, results, _ = solved(*case)
    for cid in model.combinations:
        report = cross_check(model, cid, benchfea_result=results)
        assert report.cross_checked, (cid, report.reasons)
