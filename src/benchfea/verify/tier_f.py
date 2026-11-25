"""Tier F: coverage of the design checks that are part of a reference solution.

A design result is publishable only if every limit state it USES has tier-F coverage in the test
suite: (a) a hand calculation worked from AISC 360-22 and (b) a reproduction of a printed value of the
AISC Design Examples v16.0.  The coverage is machine readable (``design_coverage.json``, generated from
the design tests by ``tools/build_design_coverage.py``) and is keyed by the ITEMS of a check's clause:

* every equation and table the check cites (``"Eq. E3-2"``, ``"Table E7.1"``),
* a section cited without an equation (``"Sec. E2"``, ``"Sec. D1"``),
* the cases of Tables B4.1a / B4.1b used to classify the section (``"Table B4.1a Case 5"``),
* ``"engine rule: ..."`` for a rule that is the engine's and not the Specification's.

The clause of a check lists the branch actually taken (Eq. E3-2 or E3-3, F2-2 or F2-3, ...), so the
items identify the limit state at the level of the equations evaluated.

"Used" (the rule of ``DesignMemberResult.passed``, docs/ARCHITECTURE.md: a limit state with a nonzero
demand): a check whose utilisation exceeds ``POLICY["F"]["negligibleUtilisation"]`` (a zero force is a
round-off number in a results document).  A check that is not used contributes no utilisation; its
capacity is counted as not covered by this tier and is not part of what "verified" covers.

The tier fails, naming the checks, when a used check (1) is ``notChecked`` or ``refused``, or (2) cites
an item without (a) coverage or without (b) coverage at strength level (an item of which only an
intermediate value was reproduced counts as PARTIAL and fails as well).
"""
from __future__ import annotations

import json
import re
import time
from functools import lru_cache
from pathlib import Path

from .policy import POLICY
from .record import TierResult

__all__ = ["tier_f", "clause_items", "coverage_registry", "item_coverage", "COVERAGE_FILE"]

COVERAGE_FILE = Path(__file__).with_name("design_coverage.json")
_CASE = re.compile(r"Table B4\.1[ab] Case \d+")
_EDITION = re.compile(r"^\s*AISC 360-\d+\s+")


def clause_items(check) -> list:
    """Items of a design check (see the module docstring); ``check`` is a ``DesignCheck`` or the dict
    of a design document (``clause`` and ``details``)."""
    clause = check["clause"] if isinstance(check, dict) else check.clause
    details = check.get("details", {}) if isinstance(check, dict) else check.details
    items = []

    def add(x):
        if x and x not in items:
            items.append(x)
    for part in _EDITION.sub("", clause).split(";"):
        part = part.strip()
        if not part:
            continue
        if part.lower().startswith("engine"):
            add("engine rule: " + part.split(None, 1)[-1] if " " in part else "engine rule")
            continue
        if part.startswith("Chapter "):
            continue                                   # refers to strengths that have their own checks
        note = re.sub(r"\(.*?\)", "", part).strip()
        head, _, rest = note.partition("Eq.")
        section = head.strip().rstrip(",").strip()
        found = False
        tables = re.findall(r"Table\s+[A-Z0-9][\w.\-]*", note)
        for t in tables:
            add(re.sub(r"\s+", " ", t))
            found = True
        if rest:
            for eq in re.findall(r"[A-Z]+\d*-\d+[a-z]?(?:-\d+[a-z]?)?", re.sub(r"Table\s+\S+", "", rest)):
                add("Eq. " + eq)
                found = True
        if not found:
            sec = re.sub(r"Table\s+\S+", "", section).strip().rstrip(",").strip()
            for sname in [x.strip() for x in sec.split(",") if x.strip()]:
                add("Sec. " + sname)
    for key in details:
        for m in _CASE.findall(str(key)):
            add(m)
    return items


@lru_cache(maxsize=1)
def coverage_registry() -> dict:
    """The registry as stored (``{"items": {item: {"a": [test ids], "b": [{"test", "level", "printed",
    "examples": [{"example", "pages", "title"}]}]}}}``); an empty registry if the file is missing."""
    try:
        return json.loads(COVERAGE_FILE.read_text(encoding="utf-8"))
    except OSError:
        return {"items": {}, "missing": True}


def item_coverage(item: str, registry: dict | None = None) -> dict:
    """``{"a": bool, "b": "strength" | "function" | "detail" | None, "status": "covered" | "partial" | "none", ...}``."""
    reg = registry if registry is not None else coverage_registry()
    e = reg.get("items", {}).get(item, {"a": [], "b": []})
    levels = {b["level"] for b in e["b"]}
    b = "strength" if "strength" in levels else ("function" if "function" in levels else ("detail" if levels else None))
    status = "covered" if (e["a"] and b in ("strength", "function")) else ("partial" if (e["a"] or b) else "none")
    examples = sorted({x["example"] for bb in e["b"] if bb["level"] == (b or "") for x in bb["examples"]})
    return {"a": bool(e["a"]), "b": b, "status": status, "aTests": list(e["a"]),
            "bTests": sorted({bb["test"] for bb in e["b"] if bb["level"] == (b or "")}), "examples": examples}


def _lacking(cov: dict) -> str:
    miss = []
    if not cov["a"]:
        miss.append("no hand calculation (a)")
    if cov["b"] is None:
        miss.append("no Design Examples reproduction (b)")
    elif cov["b"] == "detail":
        miss.append("Design Examples reproduce only an intermediate value (b partial: "
                    + ", ".join(cov["examples"][:4]) + ")")
    return "; ".join(miss)


def tier_f(design_results, registry: dict | None = None) -> TierResult:
    """``design_results``: ``benchfea.results.DesignResults``."""
    t0 = time.perf_counter()
    reg = registry if registry is not None else coverage_registry()
    used, unused, blocked, lacking = {}, [], [], {}
    n_checks = 0
    eps = POLICY["F"]["negligibleUtilisation"]
    scale = {}                                         # largest demand per dimension (for checks without a ratio)
    for m in design_results.members.values():
        for c in [c for cs in m.by_combination.values() for c in cs] + list(m.checks):
            if isinstance(c.demand, (int, float)):
                scale[c.dimension] = max(scale.get(c.dimension, 0.0), abs(c.demand))
    for mid, m in design_results.members.items():
        seen = set()
        rows = [(cid, c) for cid, cs in m.by_combination.items() for c in cs] + [(c.combination, c) for c in m.checks]
        for cid, c in rows:
            key = (c.name, cid)
            if key in seen:
                continue
            seen.add(key)
            n_checks += 1
            where = f"member {mid}, combination {cid or '-'}: {c.name}"
            if c.status == "ok" and isinstance(c.ratio, (int, float)):
                negligible = abs(c.ratio) <= eps
            else:
                negligible = (isinstance(c.demand, (int, float))
                              and abs(c.demand) <= eps * scale.get(c.dimension, 0.0))
            if negligible:
                unused.append(where)
                continue
            if c.status != "ok":
                blocked.append(f"{where}: {c.status} ({c.reason})")
                continue
            for it in clause_items(c):
                used.setdefault(it, []).append(where)
    for it, wheres in used.items():
        cov = item_coverage(it, reg)
        if cov["status"] != "covered":
            lacking[it] = {"lacking": _lacking(cov), "coverage": cov, "checks": sorted(set(wheres))[:8],
                           "nChecks": len(set(wheres))}
    failed = []
    if reg.get("missing"):
        failed.append("the coverage registry design_coverage.json is missing (run tools/build_design_coverage.py)")
    failed += [f"{it}: {v['lacking']} -- used by {v['nChecks']} check(s), e.g. {v['checks'][0]}" for it, v in lacking.items()]
    failed += blocked
    passed = not failed
    if passed:
        summary = (f"{n_checks} design check(s); every limit state used ({len(used)} clause item(s)) has (a) hand and "
                   "(b) Design Examples coverage")
    else:
        summary = (f"NOT publishable: {len(lacking)} clause item(s) without full coverage, {len(blocked)} used "
                   "check(s) notChecked / refused: " + "; ".join(failed[:3]))
    details = {"itemsUsed": {it: item_coverage(it, reg) for it in sorted(used)}, "lacking": lacking,
               "notCheckedOrRefused": blocked, "zeroDemandChecksNotCovered": len(unused), "failed": failed,
               "designStatus": design_results.status,
               "registry": {"file": COVERAGE_FILE.name, "designExamples": reg.get("designExamples")}}
    tr = TierResult("F", True, passed, summary, [], dict(POLICY["F"]), details, time.perf_counter() - t0)
    tr.category = "" if passed else "not-covered"
    return tr
