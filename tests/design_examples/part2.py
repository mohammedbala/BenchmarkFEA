"""Design Examples data of the flexure / combined-forces tests (Chapters F and H).

Same document, same SINGLE tolerance policy as ``tests/design_examples/__init__.py``
(``assert_printed`` is re-exported unchanged); the records live in their own file,
``aisc_design_examples_v16_part2.json``.
"""
import json
from pathlib import Path

from . import FT, POLICY, REL_TOL, assert_printed          # noqa: F401  (the one policy)

_DATA = json.loads((Path(__file__).parent / "aisc_design_examples_v16_part2.json").read_text(encoding="utf-8"))
SOURCE = _DATA["source"]


def example(example_id: str) -> dict:
    return _DATA["examples"][example_id]
