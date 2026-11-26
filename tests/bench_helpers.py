"""Shared pieces of the benchmark tests (tasks, scoring, environment, dataset): three small cells of the
continuous-beam family, built once per test session into a temporary cache."""
import atexit
import functools
import os
import shutil
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from benchfea.bench.tasks import build_case, make_task  # noqa: E402

CACHE = tempfile.mkdtemp(prefix="benchfea-bench-tests-")
atexit.register(shutil.rmtree, CACHE, True)

#: name -> (family, curve ball, seed, difficulty).  The curve balls are pairs of the TRAIN split for this
#: family (the held-out combinations of benchfea/bench/splits.json follow a rule on names; when the rule or
#: the applicability table changes, pick other train pairs here -- test_helper_cells_are_in_the_expected_split).
CELLS = {
    "change": ("continuous_beam", "support_fixed_to_pinned", 0, 2),          # changes-analysis, train
    "mechanism": ("continuous_beam", "support_removed", 0, 1),               # creates-mechanism, train
    "design_only": ("continuous_beam", "design_bracing_points_changed", 0, 2),   # design-only: tier F not passed
    "held_out_seed": ("continuous_beam", "support_fixed_to_pinned", 4, 2),   # seed % 5 == 4
}


@functools.lru_cache(maxsize=None)
def case(name: str) -> dict:
    return build_case(*CELLS[name], cache_dir=CACHE)


@functools.lru_cache(maxsize=None)
def task(name: str, task_type: str, mode: str = "tool-allowed") -> tuple:
    return make_task(case(name), task_type, mode)


class FullCoverage(dict):
    """TEST DEVICE ONLY: a tier-F registry in which every clause item counts as covered.  With the real
    registry few generated models have verified design results (docs/BENCHMARK.md section 9) and the
    small cells of these tests do not, so the code paths of design references -- design tasks,
    design-ratio quantities -- are exercised here by pretending the coverage exists.  Nothing built
    with it is a verified reference."""

    def get(self, key, default=None):
        return {"a": ["test-device"], "b": [{"test": "test-device", "level": "strength", "printed": "",
                                             "examples": []}]}


def full_coverage_registry() -> dict:
    return {"items": FullCoverage(), "designExamples": "TEST DEVICE: everything counted as covered"}


def leaves(value):
    """Every scalar of a JSON value."""
    if isinstance(value, dict):
        for v in value.values():
            yield from leaves(v)
    elif isinstance(value, (list, tuple)):
        for v in value:
            yield from leaves(v)
    else:
        yield value
