#!/usr/bin/env python3
"""Build the tier-F coverage registry ``src/benchfea/verify/design_coverage.json`` FROM THE DESIGN TESTS.

    python3 tools/build_design_coverage.py            # writes the registry
    python3 tools/build_design_coverage.py --check    # exit 1 if the stored registry is out of date

The registry maps every element of a design check's clause (an equation, a table case, a section --
``benchfea.verify.tier_f.clause_items``) to

(a) the hand-calculation tests (``tests/test_design_*.py::test_hand_*``: arithmetic written out from
    AISC 360-22 in the test body) in which a check with that item was produced with a result AND its nominal
    or available strength (measured value of a limit check) was read by the test code -- a check that the
    engine merely produced on the way (the other limit states of the same member) earns no credit --, and
(b) the tests that reproduce a printed value of the AISC Design Examples v16.0: a check carrying the item
    whose NOMINAL or AVAILABLE STRENGTH (or measured value / limit of a limit check) was passed to
    ``design_examples.assert_printed`` (the one tolerance policy of those tests) -- level "strength" --
    or of which only an intermediate value was compared -- level "detail".  The example ids are the ones
    the test read with ``design_examples.example`` and their page references come from the data files.

It is derived by running the design tests with three hooks (nothing in the tests is changed):
``DesignCheck.__init__`` records every check a test creates, ``assert_printed`` records every value a
test compares with a printed number, ``example`` records the examples a test reads.  A value is
attributed to a check when it is the same float (relative 1e-12).  Only tests that PASS are counted.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "src" / "benchfea" / "verify" / "design_coverage.json"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))


class Trace:
    def __init__(self):
        self.current = None
        self.tests = {}

    def start(self, nodeid):
        self.current = {"checks": [], "asserted": [], "examples": [], "read": set(), "readDetail": set()}
        self.tests[nodeid] = self.current

    def stop(self, nodeid, passed):
        if not passed:
            self.tests.pop(nodeid, None)
        self.current = None


TRACE = Trace()


def install_hooks():
    import design_examples
    from benchfea import results

    orig_assert = design_examples.assert_printed

    def assert_printed(value, printed, what=""):
        if TRACE.current is not None:
            TRACE.current["asserted"].append((float(value), what))
        return orig_assert(value, printed, what)

    design_examples.assert_printed = assert_printed
    from design_examples import part2, part3
    sources = {}
    for mod in (design_examples, part2, part3):
        mod.assert_printed = assert_printed
        data = mod._DATA
        for eid, ex in data["examples"].items():
            sources[(mod.__name__, eid)] = {"example": eid, "pages": ex.get("pages", ""), "title": ex.get("title", "")}

        def example(example_id, _orig=mod.example, _name=mod.__name__):
            if TRACE.current is not None and (_name, example_id) not in TRACE.current["examples"]:
                TRACE.current["examples"].append((_name, example_id))
            return _orig(example_id)
        mod.example = example

    orig_init = results.DesignCheck.__init__

    def __init__(self, *a, **k):
        orig_init(self, *a, **k)
        if TRACE.current is not None:
            TRACE.current["checks"].append(self)

    results.DesignCheck.__init__ = __init__

    # (a): which checks does the TEST CODE of a hand-calculation test look at?  A read of ``nominal``,
    # ``capacity`` or ``ratio`` (for a limit check without a nominal value also ``demand``) from a frame of
    # a file under tests/ is recorded ("read"); a read of ``details`` / ``inputs`` only is recorded
    # separately ("readDetail": an intermediate value, not full (a) coverage).  Reads by the engine itself
    # (e.g. the interaction check reading a capacity) are not recorded.
    tests_dir = str(ROOT / "tests") + os.sep
    orig_get = results.DesignCheck.__getattribute__

    def __getattribute__(self, name):
        if name in ("nominal", "capacity", "ratio", "demand", "details", "inputs") and TRACE.current is not None:
            if sys._getframe(1).f_code.co_filename.startswith(tests_dir):
                if name in ("details", "inputs"):
                    TRACE.current["readDetail"].add(id(self))
                elif name != "demand" or orig_get(self, "nominal") is None:
                    TRACE.current["read"].add(id(self))
        return orig_get(self, name)

    results.DesignCheck.__getattribute__ = __getattribute__
    return sources, design_examples.SOURCE.get("document", "AISC Design Examples v16.0")


class Plugin:
    def pytest_runtest_call(self, item):
        TRACE.start(item.nodeid)

    def pytest_runtest_logreport(self, report):
        if report.when == "call":
            TRACE.stop(report.nodeid, report.passed)
        elif report.when == "setup" and not report.passed:
            TRACE.stop(report.nodeid, False)


# unit factors between the value a test compares and the value stored in the check: the Design Examples
# tests work in kip, in and convert ft (x 12); tests that run the solver get N, mm back
_KIP, _IN, _FT = 4448.2216152605, 25.4, 304.8
FACTORS = (1.0, 12.0, _IN, _FT, _KIP, _KIP * _IN, _KIP * _FT, _KIP / _IN ** 2)
FACTORS = FACTORS + tuple(1.0 / f for f in FACTORS[1:])

# Equations that a Design Examples test reproduces by calling the equation FUNCTION of the engine (the one
# the limit-state check calls) instead of producing a check: (item, test, label passed to assert_printed,
# function).  An entry is used only if that test passed and compared a value under that label.
FUNCTION_LEVEL = (
    ("Eq. H1-1a", "tests/test_design_combined.py::test_design_example_h1b_interaction_with_tabulated_strengths",
     "Eq. H1-1a", "benchfea.design.aisc360_22.h_combined.h1_interaction"),
)


# Equations whose hand calculation calls the equation FUNCTION of the engine with hand-computed numbers
# (no check object is read): (item, hand test, function).  Used only if that test passed.
FUNCTION_LEVEL_A = (
    ("Eq. F1-1", "tests/test_design_flexure.py::test_hand_f1_cb_values_of_the_user_note",
     "benchfea.design.aisc360_22.f_flexure.cb_f1_1"),
)


def same(a, b):
    if a is None or b is None:
        return False
    return any(abs(a - b * f) <= 1e-12 * max(abs(a), abs(b * f)) for f in FACTORS)


def build():
    import pytest
    sources, document = install_hooks()
    from benchfea.verify.tier_f import clause_items

    files = sorted(str(p) for p in (ROOT / "tests").glob("test_design_*.py"))
    code = pytest.main(["-q", "-p", "no:cacheprovider", "--no-header", "-x", *files], plugins=[Plugin()])
    if code != 0:
        raise SystemExit(f"the design tests did not pass (pytest exit code {code}); no registry written")
    items = {}

    def entry(item):
        return items.setdefault(item, {"a": [], "aDetail": [], "b": []})

    unmatched = 0
    for nodeid, t in sorted(TRACE.tests.items()):
        test = nodeid.split("/")[-1] if "/" in nodeid else nodeid
        test = "tests/" + test.split("tests/")[-1]
        name = test.split("::")[-1]
        ok = [c for c in t["checks"] if c.status == "ok"]
        if name.startswith("test_hand_"):
            for c in ok:
                # a check the engine merely produced in this test and the test never looked at earns nothing
                key = "a" if id(c) in t["read"] else ("aDetail" if id(c) in t["readDetail"] else None)
                if key is None:
                    continue
                for it in clause_items(c):
                    if test not in entry(it)[key]:
                        entry(it)[key].append(test)
        if not t["asserted"]:
            continue
        examples = [sources[k] for k in t["examples"]]
        credit = {}                                   # item -> ("strength" | "detail", what)
        for value, what in t["asserted"]:
            hit = False
            for c in ok:
                strength = [v for v in (c.nominal, c.capacity) if isinstance(v, (int, float))]
                if c.nominal is None and isinstance(c.demand, (int, float)):
                    strength.append(c.demand)               # limit check: the measured value against a limit
                detail = [v[0] for v in c.details.values()
                          if isinstance(v, (tuple, list)) and v and isinstance(v[0], (int, float))
                          and not isinstance(v[0], bool)]
                level = ("strength" if any(same(value, v) for v in strength)
                         else "detail" if any(same(value, v) for v in detail) else None)
                if level is None:
                    continue
                hit = True
                for it in clause_items(c):
                    if level == "strength" or it not in credit:
                        if credit.get(it, ("", ""))[0] != "strength":
                            credit[it] = (level, what)
            unmatched += not hit
        for it, (level, what) in credit.items():
            entry(it)["b"].append({"test": test, "level": level, "printed": what, "examples": examples})
    for item, test, label, function in FUNCTION_LEVEL:
        for nodeid, t in sorted(TRACE.tests.items()):
            if nodeid.split("[")[0].endswith(test.split("tests/")[-1]) and any(w == label for _, w in t["asserted"]):
                entry(item)["b"].append({"test": "tests/" + nodeid.split("tests/")[-1], "level": "function",
                                         "printed": label, "function": function,
                                         "examples": [sources[k] for k in t["examples"]]})
    passed_tests = {"tests/" + n.split("tests/")[-1] for n in TRACE.tests}
    function_a = []
    for item, test, function in FUNCTION_LEVEL_A:
        if test in passed_tests:
            if test not in entry(item)["a"]:
                entry(item)["a"].append(test)
            function_a.append({"item": item, "test": test, "function": function})
    for e in items.values():
        e["aDetail"] = [t for t in e["aDetail"] if t not in e["a"]]
    registry = {
        "format": "benchfea-design-coverage", "formatVersion": "1.0.0",
        "generatedBy": "tools/build_design_coverage.py (from tests/test_design_*.py; do not edit by hand)",
        "designExamples": document,
        "ruleA": "(a) lists the test_hand_* tests whose own code read the nominal / available strength or the "
                 "utilisation (measured value of a limit check) of a check citing the item, plus the "
                 "function-level entries of 'functionLevelA'; 'aDetail' lists hand tests that read only an "
                 "intermediate value of such a check (not counted as (a)); a check the engine produced in a "
                 "hand test without the test looking at it earns no credit",
        "functionLevelA": function_a,
        "rule": "an item is covered when it has at least one (a) test and at least one (b) entry of level "
                "'strength' (a printed nominal / available strength, measured value or limit of a check that "
                "cites the item) or 'function' (the engine's function for that equation reproduces the printed "
                "value); (b) entries of level 'detail' only (an intermediate value) = partial",
        "comparedValuesNotAttributedToACheck": unmatched,
        "items": {k: items[k] for k in sorted(items)},
    }
    return registry


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    registry = build()
    text = json.dumps(registry, indent=1, ensure_ascii=False) + "\n"
    if "--check" in argv:
        if not OUT.exists() or OUT.read_text(encoding="utf-8") != text:
            print(f"{OUT} is out of date: run tools/build_design_coverage.py")
            return 1
        print("design coverage registry is up to date")
        return 0
    OUT.write_text(text, encoding="utf-8")
    n = len(registry["items"])
    full = sum(1 for v in registry["items"].values()
               if v["a"] and any(b["level"] in ("strength", "function") for b in v["b"]))
    print(f"wrote {OUT}: {n} clause items, {full} with (a) and (b, strength) coverage; "
          f"{registry['comparedValuesNotAttributedToACheck']} printed comparisons not attributed to a check")
    return 0


if __name__ == "__main__":
    sys.exit(main())
