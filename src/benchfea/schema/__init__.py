"""The published JSON Schema of the input format and helpers that read it."""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

SCHEMA_PATH = Path(__file__).with_name("benchfea-input-1.0.schema.json")


@lru_cache(maxsize=1)
def load_schema() -> dict:
    """Return the input-file JSON Schema (draft 2020-12) as a dict."""
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def validator():
    from jsonschema import Draft202012Validator

    schema = load_schema()
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)
