"""Problem library of the verification pipeline: loader, targets, extractors, verify_problem,
verify_library, the CLI, and the negative test of a deliberately corrupted target."""
import copy
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from benchfea import cli
from benchfea.crosscheck import available
from benchfea.model import build_model
from benchfea.solve import analyze
from benchfea.verify import (POLICY, ProblemError, find_problems, library_root, load_problem,
                             printed_precision_tolerance, verify_library, verify_problem)
from benchfea.verify import problems as P

ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "verification" / "problems" / "derived"
foreign = pytest.mark.skipif(not (available("opensees") and available("pynite")),
                             reason="OpenSeesPy and PyNite are needed for tier C")


def clone(tmp_path, pid, edit=None, model_edit=None):
    d = tmp_path / "lib" / "derived" / pid
    shutil.copytree(DERIVED / pid, d)
    if edit:
        data = json.loads((d / "problem.json").read_text())
        edit(data)
        (d / "problem.json").write_text(json.dumps(data))
    if model_edit:
        data = json.loads((d / "model.json").read_text())
        model_edit(data)
        (d / "model.json").write_text(json.dumps(data))
    return d


# ============================================================================== tolerance policy
def test_printed_precision_gives_half_a_unit_in_the_last_digit():
    assert printed_precision_tolerance("0.3125") == {"rel": pytest.approx(0.00005 / 0.3125)}
    assert printed_precision_tolerance("1250") == {"rel": pytest.approx(0.5 / 1250)}
    assert printed_precision_tolerance("1,250") == {"rel": pytest.approx(0.5 / 1250)}
    assert printed_precision_tolerance("-1.20e3") == {"rel": pytest.approx(5.0 / 1200)}
    assert printed_precision_tolerance("123456789.25") == {"rel": POLICY["B"]["relFloor"]}       # floor 1e-6
    assert printed_precision_tolerance("0.00") == {"abs": 0.005}


# ============================================================================== loader
def test_library_is_found_and_every_derived_problem_loads():
    assert library_root(ROOT) == ROOT / "verification" / "problems"
    dirs = [d for d in find_problems() if d.parent.name == "derived"]
    assert len(dirs) >= 17
    for d in dirs:
        spec = load_problem(d)
        assert spec.tier == "A" and spec.basis == "derived-here" and spec.id == d.name and spec.targets
        assert all(t["tolerance"] in ({"rel": POLICY["targets"]["derivedHere"]}, {"abs": 1e-9}) for t in spec.targets)
        assert not spec.notes and "derived here" in spec.source["work"]
        build_model(spec.model_path)                             # a valid input file


@pytest.mark.parametrize("edit, message", [
    (lambda p: p.pop("source"), "'source' is required"),
    (lambda p: p.update(tier="C"), "tier must be 'A' or 'B'"),
    (lambda p: p.update(basis="memory"), "basis must be one of"),
    (lambda p: p.update(id="other"), "must equal the directory name"),
    (lambda p: p.update(status="quarantined"), "needs quarantine.diagnosis"),
    (lambda p: p["targets"][0].update(tolerance={"rel": 1e-9}), "below the policy floor"),
    (lambda p: p["targets"][0].update(tolerance={"rel": 0.01, "abs": 1.0}), 'must be {"rel": x} or {"abs": y}'),
    (lambda p: p["targets"][0].update(tolerance={"rel": -1.0}), "positive finite number"),
    (lambda p: p["targets"][0].update(value=0.0), "can never be met"),
    (lambda p: p["targets"][0].update(value="12"), "'value' must be a finite number"),
    (lambda p: p["targets"][0].update(extract="magic"), "exactly one of 'path' and 'extract'"),
    (lambda p: p["targets"].__setitem__(0, {"extract": "magic", "value": 1.0, "tolerance": {"rel": 1e-6}}),
     "unknown extractor"),
    (lambda p: p.update(targets=[]), "non-empty list"),
    (lambda p: p.update(tier="B"), "basis 'derived-here' is tier A"),
])
def test_loader_refuses_what_does_not_follow_the_format(tmp_path, edit, message):
    d = clone(tmp_path, "beam-point-load", edit)
    with pytest.raises(ProblemError) as e:
        load_problem(d)
    assert message in str(e.value)
    assert verify_problem(d).outcome == "invalid"


def test_loader_notes_departures_from_the_tolerance_policy(tmp_path):
    d = clone(tmp_path, "beam-point-load", lambda p: p["targets"][1].update(tolerance={"rel": 1e-3}))
    assert "looser than the policy" in load_problem(d).notes[0]

    def published(p):
        p.update(tier="B", basis="published")
        p["source"]["location"] = "Example 1"
    d2 = clone(tmp_path / "x", "beam-point-load", published)
    notes = load_problem(d2).notes
    assert len(notes) == 9 and all("no precisionNote" in n for n in notes)


# ============================================================================== targets and extractors
@pytest.fixture(scope="module")
def beam_doc():
    model = build_model(DERIVED / "beam-point-load" / "model.json")
    return analyze(model).to_document(model)


def test_paths_and_extractors(beam_doc):
    c = beam_doc["combinations"]["C"]
    assert P.resolve_path(beam_doc, "/combinations/C/reactions/a/Fz") == c["reactions"]["a"]["Fz"]
    assert P.resolve_path(beam_doc, "/combinations/C/members/m/stations/Mz/3") == c["members"]["m"]["stations"]["Mz"][3]
    for bad, msg in (("combinations/C", "must start with"), ("/combinations/X", "no key 'X'"),
                     ("/combinations/C/members/m/stations/Mz/99", "not an index"), ("/combinations/C", "not a number")):
        with pytest.raises(P.TargetError, match=msg):
            P.resolve_path(beam_doc, bad)
    X = P.EXTRACTORS
    st = c["members"]["m"]["stations"]
    k = st["x"].index(1.8) if 1.8 in st["x"] else min(range(len(st["x"])), key=lambda i: abs(st["x"][i] - 1.8))
    t = {"extract": "memberStation", "member": "m", "quantity": "Vy", "relative": 0.3}
    left, right = X["memberStation"](beam_doc, None, {**t, "side": "left"}), X["memberStation"](beam_doc, None, {**t, "side": "right"})
    assert (left, right) == (st["Vy"][k], st["Vy"][k + 1]) and abs(left - right) == pytest.approx(37.5)   # 1.5 x 25 kN
    by_x = {"extract": "memberStation", "member": "m", "quantity": "Vy", "x": st["x"][k]}
    assert X["memberStation"](beam_doc, None, by_x) == left
    with pytest.raises(P.TargetError, match="not a result station"):
        X["memberStation"](beam_doc, None, {**t, "relative": 0.33})
    with pytest.raises(P.TargetError, match="exactly one of"):
        X["memberStation"](beam_doc, None, {"extract": "memberStation", "member": "m", "quantity": "Vy"})
    assert X["memberExtreme"](beam_doc, None, {"extract": "memberExtreme", "member": "m", "quantity": "Mz",
                                               "kind": "max"}) == max(st["Mz"])
    assert X["memberExtreme"](beam_doc, None, {"extract": "x", "member": "m", "quantity": "uy"}) == max(map(abs, st["uy"]))
    assert X["nodeDisplacement"](beam_doc, None, {"extract": "x", "node": "a", "dof": "ry"}) == c["displacements"]["a"]["ry"]
    assert X["reaction"](beam_doc, None, {"extract": "x", "node": "b", "component": "Fz"}) == c["reactions"]["b"]["Fz"]
    assert X["sumReactions"](beam_doc, None, {"extract": "x", "component": "Fz"}) == pytest.approx(37.5)
    assert X["memberEndForce"](beam_doc, None, {"extract": "x", "member": "m", "end": "i", "component": "Fy"}) \
        == c["members"]["m"]["endForces"]["i"]["Fy"]
    with pytest.raises(P.TargetError, match="no buckling analysis"):
        X["bucklingFactor"](beam_doc, None, {"extract": "bucklingFactor", "combination": "C"})
    with pytest.raises(P.TargetError, match="no design results"):
        X["designCheck"](beam_doc, None, {"extract": "designCheck", "member": "m", "check": "shear"})


def test_buckling_and_design_extractors():
    model = build_model(DERIVED / "euler-column-fixed-pinned" / "model.json")
    doc = analyze(model).to_document(model)
    f = P.EXTRACTORS["bucklingFactor"](doc, None, {"extract": "bucklingFactor", "mode": 1})
    assert f == doc["buckling"]["C"]["factors"][0]
    with pytest.raises(P.TargetError, match="mode 5"):
        P.EXTRACTORS["bucklingFactor"](doc, None, {"extract": "bucklingFactor", "mode": 5})
    from benchfea import design
    from test_verify_pipeline import DESIGN_BEAM
    m = build_model(copy.deepcopy(DESIGN_BEAM))
    r = analyze(m)
    ddoc = design.check(m, r).to_document(m)
    t = {"extract": "designCheck", "member": "m", "check": "flexure, major axis: yielding", "field": "nominal"}
    Mp = P.EXTRACTORS["designCheck"](r.to_document(m), ddoc, t)
    assert Mp > 0 and P.EXTRACTORS["designCheck"](None, ddoc, {**t, "field": "ratio"}) == pytest.approx(
        1.5 * 20.0 * 6000.0 ** 2 / 8.0 / (0.9 * Mp), rel=1e-3)      # second order without axial force = first order
    with pytest.raises(P.TargetError, match="match"):
        P.EXTRACTORS["designCheck"](None, ddoc, {**t, "check": "flexure"})


# ============================================================================== verify_problem
@foreign
def test_derived_problem_passes_with_its_model_verified():
    r = verify_problem(DERIVED / "beam-point-load")
    assert r.outcome == "passed" and not r.failures and len(r.targets) == 9 and all(t.passed for t in r.targets)
    assert r.largest_relative < 1e-9 and r.cross_checked is True
    rec = r.model_record
    assert rec.status == "verified" and rec.tiers_passed == ["A", "C", "D", "E"]
    a = rec.tiers["A"]                                          # the automatic oracle comparison AND the targets
    assert "against the ODE oracle" in a.summary and "9 target(s)" in a.summary
    d = r.to_dict()
    json.dumps(d)
    assert d["outcome"] == "passed" and d["modelRecord"]["status"] == "verified" and d["crossChecked"] is True


def test_corrupted_target_fails_the_problem(tmp_path):
    def corrupt(p):
        p["targets"][4]["value"] *= 1.0 + 1.0e-4               # midspan deflection off by 0.01 %
    r = verify_problem(clone(tmp_path, "beam-point-load", corrupt), tiers="AD")
    assert r.outcome == "failed" and [t.passed for t in r.targets].count(False) == 1
    bad = r.targets[4]
    assert not bad.passed and bad.relative == pytest.approx(1.0e-4, rel=1e-3) and bad.tolerance == {"rel": 1e-6}
    assert len(r.failures) == 1 and "memberStation" in r.failures[0] and "expected" in r.failures[0]
    assert r.model_record.tiers["A"].passed is False and r.model_record.status == "quarantined"
    assert "1 MISSED" in r.model_record.tiers["A"].summary


def test_target_that_cannot_be_evaluated_fails_with_the_reason(tmp_path):
    def edit(p):
        p["targets"][0]["path"] = "/combinations/C/reactions/zz/Fz"
    r = verify_problem(clone(tmp_path, "beam-point-load", edit), tiers="AD")
    assert r.outcome == "failed" and r.targets[0].actual is None and "no key 'zz'" in r.targets[0].error


def test_corrupted_model_fails_the_targets(tmp_path):
    def stiffer(m):
        m["sections"]["S"]["Iz"] *= 1.001
    r = verify_problem(clone(tmp_path, "beam-point-load", model_edit=stiffer), tiers="AD")
    assert r.outcome == "failed" and any(not t.passed for t in r.targets)
    assert all(t.passed for t in r.targets[:2])                 # the reactions of this beam do not depend on EI


def test_quarantined_problem_is_reported_and_never_counted_as_passed(tmp_path):
    def quarantine(p):
        p.update(status="quarantined", quarantine={"diagnosis": "test: suspected misprint"})
    r = verify_problem(clone(tmp_path, "beam-point-load", quarantine), tiers="AD")
    assert r.outcome == "quarantined" and all(t.passed for t in r.targets)      # still evaluated
    assert r.quarantine["diagnosis"] == "test: suspected misprint"
    lib = verify_library(tmp_path / "lib", tiers="AD")
    assert lib.ok and lib.counts() == {"passed": 0, "failed": 0, "quarantined": 1, "invalid": 0}


def test_unstable_or_refused_model_fails_the_problem(tmp_path):
    def unsupported(m):
        m["supports"] = {"a": {"dofs": {"uz": "fixed"}}}
    r = verify_problem(clone(tmp_path, "beam-point-load", model_edit=unsupported), tiers="AD")
    assert r.outcome == "failed" and "no results" in r.failures[0]


@foreign
def test_problem_whose_model_only_one_solver_represents_passes_but_is_not_cross_checked():
    r = verify_problem(DERIVED / "timoshenko-fixed-fixed")
    assert r.outcome == "passed" and r.cross_checked is False
    assert r.model_record.status == "quarantined" and r.model_record.quarantine["category"] == "not-representable"
    r = verify_problem(DERIVED / "euler-column-fixed-free")
    assert r.outcome == "passed" and r.model_record.quarantine["category"] == "not-representable"
    assert "linear buckling" in r.model_record.tiers["C"].summary


@foreign
def test_library_of_derived_problems_and_its_quarantined_cases():
    lib = verify_library(DERIVED.parent, only=["derived/"])
    by = {r.id: r for r in lib.problems}
    assert lib.ok and len(by) >= 17
    assert {k for k, r in by.items() if r.outcome == "quarantined"} == {"beam-column-02", "beam-column-05"}
    for k in ("beam-column-02", "beam-column-05"):               # targets reproduced, tier C not: documented
        assert all(t.passed for t in by[k].targets) and by[k].model_record.tiers["A"].passed
        assert by[k].model_record.tiers["C"].passed is False and "LEADING-ORDER" in by[k].quarantine["diagnosis"]
    assert all(r.outcome == "passed" for k, r in by.items() if not k.startswith("beam-column-0"))
    assert max(r.largest_relative for r in by.values()) < 1e-7
    d = lib.to_dict()
    json.dumps(d)
    assert d["counts"]["failed"] == 0 and "certified" in d["meaning"]


def test_derived_targets_are_reproducible_from_the_script():
    out = subprocess.run([sys.executable, str(ROOT / "verification" / "derive" / "derive_problems.py"), "--check"],
                         capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stdout + out.stderr
    assert "derived problems reproduce" in out.stdout


# ============================================================================== command line
def test_cli_help_has_an_example(capsys):
    with pytest.raises(SystemExit):
        cli.main(["verify", "-h"])
    out = capsys.readouterr().out
    assert "benchfea verify --all" in out and "--design" in out and "--budget" in out


@foreign
def test_cli_model_file_problem_directory_and_exit_codes(tmp_path, capsys):
    model = DERIVED / "beam-point-load" / "model.json"
    assert cli.main(["verify", str(model)]) == 0
    out = capsys.readouterr().out
    assert "VERIFIED" in out and "tier C [independent solvers]: passed" in out and "certified" in out
    rec_file = tmp_path / "rec.json"
    assert cli.main(["verify", str(model), "--json", "-o", str(rec_file)]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["status"] == "verified" and "reference" not in data
    assert json.loads(rec_file.read_text())["reference"]["format"] == "benchfea-results"
    assert cli.main(["verify", str(model), "--tiers", "AD"]) == 5           # a tier left out: not verified
    capsys.readouterr()
    assert cli.main(["verify", str(DERIVED / "beam-point-load")]) == 0
    assert "PASSED" in capsys.readouterr().out
    bad = clone(tmp_path, "beam-point-load", lambda p: p["targets"][0].update(value=p["targets"][0]["value"] * 1.01))
    assert cli.main(["verify", str(bad), "--tiers", "AD"]) == 5
    assert "MISS" in capsys.readouterr().out
    assert cli.main(["verify", str(tmp_path / "lib"), "--tiers", "AD", "--json"]) == 5      # a library directory
    assert json.loads(capsys.readouterr().out)["counts"]["failed"] == 1
    assert cli.main(["verify"]) == 2
    invalid = tmp_path / "bad.json"
    invalid.write_text("{}")
    assert cli.main(["verify", str(invalid)]) == 1                          # refused: no reference solution
