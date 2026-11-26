"""Generated files are what their builders give: the dataset numbers in data/DATA_CARD.md and
docs/BENCHMARK.md (and ``servedCombinations`` of splits.json), and the generated parts of examples/.

The builders are run with ``--check`` in a subprocess; a failure means "run the tool again", never
"edit the file".
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run(tool, *args):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    return subprocess.run([sys.executable, str(ROOT / "tools" / tool), *args], capture_output=True, text=True,
                          cwd=ROOT, env=env)


def test_dataset_numbers_in_the_docs_are_those_of_the_manifest():
    p = run("build_dataset_docs.py", "--check")
    assert p.returncode == 0, p.stdout + p.stderr


def test_examples_are_what_the_builder_gives():
    p = run("build_examples.py", "--check")
    assert p.returncode == 0, p.stdout + p.stderr


def test_the_docs_hold_no_dataset_count_outside_the_generated_blocks():
    """The headline numbers of the manifest do not occur in the prose of the two files."""
    sys.path.insert(0, str(ROOT / "tools"))
    try:
        import build_dataset_docs as B
    finally:
        sys.path.pop(0)
    counts = json.loads((ROOT / "data" / "tasks" / "manifest.json").read_text())["counts"]
    for path in B.DOCS:
        prose = B.MARK.sub("", path.read_text(encoding="utf-8"))
        assert "<!-- generated:" not in prose and "<!-- /generated:" not in prose, path.name
        for key in ("tasks", "cellsVerified"):
            assert f"{counts[key]} " not in prose, (path.name, key, counts[key])


def test_checks_file_belongs_to_the_manifest():
    """data/CHECKS.json: both checks were run on the manifest on disk and found nothing."""
    import hashlib
    sha = hashlib.sha256((ROOT / "data" / "tasks" / "manifest.json").read_bytes()).hexdigest()
    checks = json.loads((ROOT / "data" / "CHECKS.json").read_text())
    for key in ("scores", "audit"):
        assert checks[key]["manifestSha256"] == sha, key
        assert checks[key]["findings"] == [], key
    s = checks["scores"]
    assert s["referenceAnswerScores1"] == s["tasks"] > 0
    assert s["repair"]["givenModelScoresBelowReference"] == s["repair"]["tasks"]
    assert checks["audit"]["split"]["heldOutTasksOnABaseModelOfTheTrainSplit"] == 0
    assert checks["audit"]["split"]["curveBalledModelsInBothSplits"] == 0


@pytest.mark.parametrize("document, order", [
    ({"combinations": {"a": {}}}, "first-order"),
    ({"combinations": {"a": {}, "b": {"analysisType": "secondOrder"}}}, "second-order"),
    ({"analysis": {"type": "secondOrder"}, "combinations": {"a": {}}}, "second-order"),
    ({"analysis": {"type": "secondOrder"}, "combinations": {"a": {"analysisType": "firstOrder"}}}, "first-order"),
])
def test_analysis_order_of_a_document(document, order):
    from benchfea.bench.dataset import analysis_order
    assert analysis_order(document) == order


def test_split_file_is_the_rule_applied_to_the_applicability_table():
    p = run("apply_split_rule.py", "--check")
    assert p.returncode == 0, p.stdout + p.stderr
