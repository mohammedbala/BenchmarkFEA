"""Design Examples data of the independent review spot-check (tests/test_design_review.py).

Same document, same SINGLE tolerance policy as ``tests/design_examples/__init__.py``
(``assert_printed`` is re-exported unchanged); the records -- worked cases that the other
design tests do not use -- live in ``aisc_design_examples_v16_part3.json``.
"""
import json
from pathlib import Path

from . import FT, POLICY, REL_TOL, assert_printed          # noqa: F401  (the one policy)

_DATA = json.loads((Path(__file__).parent / "aisc_design_examples_v16_part3.json").read_text(encoding="utf-8"))
SOURCE = _DATA["source"]


def example(example_id: str) -> dict:
    return _DATA["examples"][example_id]
