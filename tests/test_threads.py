"""One BLAS thread: the matrices of this package are small and dense, and a multi-threaded BLAS was
measured to make the suite hundreds of times slower on a many-core machine.  The variables must be set
before numpy is loaded, so they are set by ``benchfea/__init__.py`` (first thing any entry point of the
package imports: the CLI, the dataset workers started with ``spawn``, ``python -m
benchfea.crosscheck.worker``) and by ``tests/conftest.py``."""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def _run(code: str, **extra) -> str:
    env = {k: v for k, v in os.environ.items() if k not in VARIABLES}
    env["PYTHONPATH"] = str(SRC) + os.pathsep + env.get("PYTHONPATH", "")
    env.update(extra)
    out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    return out.stdout.strip()


def test_the_test_session_runs_with_one_blas_thread():
    assert all(os.environ.get(v) for v in VARIABLES)              # set by conftest.py, or by the caller


@pytest.mark.parametrize("entry", ["benchfea", "benchfea.cli", "benchfea.bench.dataset", "benchfea.crosscheck.worker",
                                   "benchfea.crosscheck.runner"])
def test_package_sets_one_thread_before_numpy_is_imported(entry):
    code = ("import builtins, os, sys\n"
            "real = builtins.__import__\n"
            "state = {}\n"
            "def spy(name, *a, **k):\n"
            "    if name.split('.')[0] in ('numpy', 'scipy') and 'at' not in state:\n"
            f"        state['at'] = [os.environ.get(v) for v in {VARIABLES!r}]\n"
            "    return real(name, *a, **k)\n"
            "builtins.__import__ = spy\n"
            f"import {entry}\n"
            "import json\n"
            "print(json.dumps([state.get('at'), [os.environ.get(v) for v in " + repr(VARIABLES) + "]]))\n")
    at_numpy_import, afterwards = json.loads(_run(code))
    assert afterwards == ["1", "1", "1"]
    assert at_numpy_import in (None, ["1", "1", "1"])             # None: this entry point loads no numpy


def test_a_value_set_by_the_caller_is_kept():
    out = _run("import os, benchfea; print(os.environ['OMP_NUM_THREADS'], os.environ['OPENBLAS_NUM_THREADS'])",
               OMP_NUM_THREADS="3")
    assert out == "3 1"
