"""docs/OUTPUT_REFERENCE.md names every key of the results and design documents of the five worked
examples.  The test runs the examples and walks the two documents; it does not check the descriptions,
only that no key is written without being named in the reference.

Keys that are identifiers of the input file (nodes, members, combinations, design members, groups) and
the open-ended maps (names of inputs and intermediate values of a check, names of design variables, the
per-iteration history of a solver) are not expected in the document.
"""
import json
import re
from pathlib import Path

import pytest

from benchfea import design as design_module
from benchfea.model import build_model
from benchfea.solve import analyze

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = (ROOT / "docs" / "OUTPUT_REFERENCE.md").read_text(encoding="utf-8")
NAMED = set(re.findall(r"[A-Za-z][A-Za-z0-9]*", " ".join(re.findall(r"`([^`]+)`", REFERENCE))))
NAMES = ["01_simple_beam", "02_continuous_beam", "03_portal_frame", "04_diaphragm_frame", "05_braced_frame"]
# children of these keys are names chosen by the input file or by a check, not keys of the format
OPEN_MAPS = {"inputs", "details", "parameters", "tauB", "subdivisions", "history", "subdivisionCheck",
             "stiffnessReduction", "notionalLoads"}
# layout of the entries of the open maps that the reference does state
ENTRY_KEYS = {"inputs": {"value", "unit", "source"}, "details": {"value", "unit"},
              "parameters": {"value", "unit", "source", "path"}}


def _keys(node, ids, found, parent=None):
    if isinstance(node, dict):
        for key, value in node.items():
            if parent in OPEN_MAPS:
                if parent in ENTRY_KEYS and isinstance(value, dict):
                    assert set(value) == ENTRY_KEYS[parent], (parent, key, sorted(value))
                continue
            if key not in ids:
                found.add(key)
            _keys(value, ids, found, key)
    elif isinstance(node, list):
        for value in node:
            _keys(value, ids, found, parent)


@pytest.mark.parametrize("name", NAMES)
def test_every_key_of_the_output_documents_is_named_in_the_reference(name):
    raw = json.loads((ROOT / "examples" / f"{name}.json").read_text(encoding="utf-8"))
    ids = set()
    for part in ("nodes", "members", "combinations", "loadCases", "sections", "materials", "constraints"):
        ids |= set(raw.get(part, {}))
    ids |= set(raw.get("design", {}).get("members", {})) | set(raw.get("design", {}).get("groups", {}))
    model = build_model(ROOT / "examples" / f"{name}.json")
    results = analyze(model)
    documents = {"results": results.to_document(model),
                 "design": design_module.check(model, results).to_document(model)}
    for kind, document in documents.items():
        found = set()
        _keys(document, ids, found)
        missing = sorted(found - NAMED)
        assert not missing, f"{kind} document of {name}: keys not named in docs/OUTPUT_REFERENCE.md: {missing}"


def test_reference_states_formats_and_claims_no_review():
    assert "`benchfea-results`" in REFERENCE and "`benchfea-design`" in REFERENCE
    assert "Nothing in either document says that a result has been reviewed, approved or certified" in REFERENCE
