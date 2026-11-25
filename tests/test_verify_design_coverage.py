"""Tier F: the machine-readable coverage registry of the design limit states and the clause items."""
import ast
import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchfea.results import DesignCheck
from benchfea.verify import clause_items, coverage_registry, item_coverage
from benchfea.verify.tier_f import COVERAGE_FILE, tier_f
from benchfea.results import DesignMemberResult, DesignResults

ROOT = Path(__file__).resolve().parents[1]


def test_clause_items_identify_the_equations_tables_and_cases_used():
    c = {"clause": "AISC 360-22 E7, Eq. E7-1; E7.1, Eq. E7-2, E7-3, E7-5, Table E7.1; E3, Eq. E3-2, E3-4; E2; "
                   "B4.1, Table B4.1a; B3.1 Eq. B3-1",
         "details": {"lambda_r web (Table B4.1a Case 5)": 1, "lambda_r flange (Table B4.1a Case 2)": 2}}
    assert sorted(clause_items(c)) == sorted([
        "Eq. E7-1", "Eq. E7-2", "Eq. E7-3", "Eq. E7-5", "Table E7.1", "Eq. E3-2", "Eq. E3-4", "Sec. E2",
        "Table B4.1a", "Eq. B3-1", "Table B4.1a Case 5", "Table B4.1a Case 2"])
    assert clause_items({"clause": "AISC 360-22 D1 (User Note)"}) == ["Sec. D1"]
    assert clause_items({"clause": "AISC 360-22 H1.1, Eq. H1-1b; Chapter E; Chapter F"}) == ["Eq. H1-1b"]
    assert clause_items({"clause": "AISC 360-22 G5, Eq. G5-1, G5-2a, G5-2b; G1; B3.2 Eq. B3-2"}) == [
        "Eq. G5-1", "Eq. G5-2a", "Eq. G5-2b", "Sec. G1", "Eq. B3-2"]
    assert clause_items({"clause": "AISC 360-22 H3 (not implemented); engine screening rule"}) == [
        "Sec. H3", "engine rule: screening rule"]
    assert clause_items({"clause": "AISC 360-22 L1, L2 (the Specification sets no numeric limit: the limit is "
                                   "the user's)"}) == ["Sec. L1", "Sec. L2"]
    assert clause_items({"clause": "AISC 360-22 F2.2, Eq. F2-2, F2-5, F2-6, F2-7, F2-8a; F2.1, Eq. F2-1; B4.1, "
                                   "Table B4.1b; F1, Eq. F1-1; F1(a); B3.1 Eq. B3-1"})[:5] == [
        "Eq. F2-2", "Eq. F2-5", "Eq. F2-6", "Eq. F2-7", "Eq. F2-8a"]


def test_registry_refers_only_to_tests_that_exist_and_cites_the_examples():
    reg = coverage_registry()
    assert reg["format"] == "benchfea-design-coverage" and len(reg["items"]) >= 70
    names = {}
    for f in (ROOT / "tests").glob("test_design_*.py"):
        tree = ast.parse(f.read_text(encoding="utf-8"))
        names[f"tests/{f.name}"] = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    examples = {}
    for f in (ROOT / "tests" / "design_examples").glob("*.json"):
        examples.update(json.loads(f.read_text(encoding="utf-8"))["examples"])
    for item, e in reg["items"].items():
        for test in e["a"]:
            file, func = test.split("::")
            assert func.split("[")[0] in names[file] and func.startswith("test_hand_"), (item, test)
        for b in e["b"]:
            file, func = b["test"].split("::")
            assert func.split("[")[0] in names[file], (item, b["test"])
            assert b["level"] in ("strength", "detail", "function") and b["examples"], (item, b)
            for x in b["examples"]:
                assert x["example"] in examples and x["pages"], (item, x)


def test_coverage_of_known_items():
    full = item_coverage("Eq. E3-2")
    assert full["status"] == "covered" and full["a"] and full["b"] == "strength" and "E.1D" in full["examples"]
    assert item_coverage("Eq. F2-2")["status"] == "covered" and "F.1-2B" in item_coverage("Eq. F2-2")["examples"]
    assert item_coverage("Eq. H1-1a")["b"] == "function"         # reproduced by the equation function (H.1B)
    tb = item_coverage("Eq. E4-1")                               # the examples print only the elastic stress
    assert tb["status"] == "partial" and tb["b"] == "detail" and set(tb["examples"]) == {"E.2", "E.3"}
    assert item_coverage("Eq. G2-3")["status"] == "partial" and item_coverage("Eq. G2-3")["b"] is None
    none = item_coverage("Eq. Z9-9")
    assert none["status"] == "none" and not none["a"] and none["b"] is None


def _check(name, clause, ratio=0.5, status="ok", demand=10.0, **kw):
    return DesignCheck(name=name, clause=clause, demand=demand, capacity=None if status != "ok" else demand / ratio,
                       ratio=None if status != "ok" else ratio, passed=None if status != "ok" else True,
                       combination="U", status=status, reason="out of scope" if status != "ok" else "", **kw)


def _design(checks):
    m = DesignMemberResult("m", checks, True, by_combination={"U": checks}, status="pass")
    return DesignResults("AISC 360-22", {"m": m}, status="pass")


def test_tier_f_rules():
    ok = _check("compression: flexural buckling", "AISC 360-22 E3, Eq. E3-2, E3-4; B3.1 Eq. B3-1")
    assert tier_f(_design([ok])).passed is True
    # a limit state with partial (b) coverage, used: not publishable, and the record says which
    tors = _check("compression: torsional buckling", "AISC 360-22 E4, Eq. E4-1, E4-2; B3.1 Eq. B3-1")
    t = tier_f(_design([ok, tors]))
    assert t.passed is False and t.category == "not-covered" and set(t.details["lacking"]) == {"Eq. E4-1", "Eq. E4-2"}
    assert "only an intermediate value" in t.details["lacking"]["Eq. E4-1"]["lacking"]
    assert "torsional buckling" in t.details["lacking"]["Eq. E4-1"]["checks"][0]
    # the same limit state with a zero demand is not used: it does not block, and it is counted
    idle = _check("compression: torsional buckling", "AISC 360-22 E4, Eq. E4-1, E4-2; B3.1 Eq. B3-1", ratio=1e-16,
                  demand=1e-13)
    t = tier_f(_design([ok, idle]))
    assert t.passed is True and t.details["zeroDemandChecksNotCovered"] == 1
    # notChecked with a demand blocks; without a demand it does not
    nc = _check("shear along local y", "AISC 360-22 Chapter G", status="notChecked", demand=5.0)
    t = tier_f(_design([ok, nc]))
    assert t.passed is False and "notChecked (out of scope)" in t.details["notCheckedOrRefused"][0]
    nc0 = _check("shear along local y", "AISC 360-22 Chapter G", status="notChecked", demand=0.0)
    assert tier_f(_design([ok, nc0])).passed is True
    # an item nobody has tested
    new = _check("something new", "AISC 360-22 J4, Eq. J4-1")
    t = tier_f(_design([new]))
    assert t.passed is False and "no hand calculation (a); no Design Examples reproduction (b)" in \
        t.details["lacking"]["Eq. J4-1"]["lacking"]
    # a missing registry never passes
    assert tier_f(_design([ok]), {"items": {}, "missing": True}).passed is False


def test_registry_is_up_to_date_with_the_design_tests():
    assert COVERAGE_FILE.is_file()
    out = subprocess.run([sys.executable, str(ROOT / "tools" / "build_design_coverage.py"), "--check"],
                         capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stdout[-2000:] + out.stderr[-2000:]
