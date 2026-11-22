"""Data and the single tolerance policy for the reproductions of the AISC
Design Examples (Version 16.0) in the design tests.

``example(id)`` returns the record of one example (inputs the example lists,
values it prints, page reference).  ``assert_printed(value, printed)`` applies
THE tolerance policy -- it is defined once, in ``tolerancePolicy`` of the JSON
file, and never adjusted per case.
"""
import json
from pathlib import Path

import pytest

_DATA = json.loads((Path(__file__).parent / "aisc_design_examples_v16.json").read_text(encoding="utf-8"))

POLICY = _DATA["tolerancePolicy"]
REL_TOL = POLICY["relative"]
SOURCE = _DATA["source"]
FT = 12.0                               # in. per ft


def example(example_id: str) -> dict:
    return _DATA["examples"][example_id]


def assert_printed(value, printed, what=""):
    """``value`` (engine, full precision) agrees with the printed value within the policy."""
    assert value == pytest.approx(printed, rel=REL_TOL), \
        f"{what}: engine {value!r} vs printed {printed!r} (policy: relative {REL_TOL})"
