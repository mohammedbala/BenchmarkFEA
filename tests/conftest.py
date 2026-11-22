import os
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def _no_user_shapes_db(monkeypatch, tmp_path):
    """Tests never see a real user-supplied shapes database."""
    monkeypatch.setenv("BENCHFEA_SHAPES_DB", str(tmp_path / "absent-shapes-db.csv"))
    monkeypatch.delenv("BENCHFEA_SHAPES_COLUMNS", raising=False)
    monkeypatch.delenv("BENCHFEA_SHAPES_LENGTH_UNIT", raising=False)


@pytest.fixture
def fixtures():
    return FIXTURES
