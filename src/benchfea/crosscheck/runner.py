"""Run a foreign solver on a translated model, in a SUBPROCESS with a timeout.

OpenSees terminates the interpreter on some input errors (measured: a ``vecxz`` parallel to the member
axis, a system without free DOFs), so a foreign solver never runs in the caller's process unless
``in_process=True`` is asked for.  The job goes to the worker as a JSON file and the raw results come
back as a JSON file (stdout is not used: the libraries print to it).

    result = run("opensees", model, "LC1")          # -> ForeignResult
    python -m benchfea.crosscheck.worker job.json out.json     # the worker process

For every refinement parameter of an approximate feature (``Translation.controls``) the model is run
a second time at the coarse level.  Where the error law is known (``order`` p: error proportional to
``value**-p``, e.g. p = 2 for the number of elements per member) the estimated error of the fine run
is ``|fine - coarse| / (ratio**p - 1)``; where round-off competes with the approximation (penalty
stiffnesses, ``order = None``) twice the whole difference ``|fine - coarse|`` is taken.  The estimate is
written into the fidelity record.
"""
from __future__ import annotations

import importlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .foreign import ForeignResult

__all__ = ["run", "run_job", "available", "convergence_study", "SOLVERS", "DEFAULT_TIMEOUT"]

SOLVERS = {"opensees": ("benchfea.crosscheck.opensees", "openseespy"),
           "pynite": ("benchfea.crosscheck.pynite", "Pynite")}
DEFAULT_TIMEOUT = 300.0
PENALTY_SAFETY = 2.0


def available(solver: str) -> bool:
    """Is the optional package of this solver installed? (Does not import it.)"""
    try:
        return importlib.util.find_spec(SOLVERS[solver][1]) is not None
    except (ImportError, ValueError):
        return False


def _translator(solver: str):
    if solver not in SOLVERS:
        raise ValueError(f"unknown solver {solver!r}; known: {', '.join(SOLVERS)}")
    return importlib.import_module(SOLVERS[solver][0])


def run_job(job: dict, timeout: float = DEFAULT_TIMEOUT, in_process: bool = False) -> dict:
    """Execute a job; never raises for a solver failure: returns ``{"ok": False, "error": ...}``."""
    if not available(job["solver"]):
        return {"ok": False, "error": f"package {SOLVERS[job['solver']][1]!r} is not installed", "missing": True}
    if in_process:
        try:
            return _translator(job["solver"]).execute(job)
        except Exception as e:                               # noqa: BLE001 - any library failure is a result
            return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    with tempfile.TemporaryDirectory(prefix="benchfea-crosscheck-") as tmp:
        src, dst = Path(tmp) / "job.json", Path(tmp) / "out.json"
        src.write_text(json.dumps(job))
        env = dict(os.environ)
        root = str(Path(__file__).resolve().parents[2])      # directory that contains the benchfea package
        env["PYTHONPATH"] = root + os.pathsep + env.get("PYTHONPATH", "")
        try:
            proc = subprocess.run([sys.executable, "-m", "benchfea.crosscheck.worker", str(src), str(dst)],
                                  capture_output=True, text=True, timeout=timeout, env=env, cwd=tmp)
        except subprocess.TimeoutExpired:
            return {"ok": False, "error": f"{job['solver']} did not finish within {timeout:g} s (killed)",
                    "timeout": True}
        tail = (proc.stderr or "").strip()[-2000:]
        if not dst.exists():
            return {"ok": False, "error": f"{job['solver']} worker ended with exit code {proc.returncode} without "
                                          f"results; stderr: {tail or '(empty)'}", "returncode": proc.returncode}
        out = json.loads(dst.read_text())
        if not out.get("ok") and tail:
            out["error"] = f"{out.get('error', '')} stderr: {tail}"
        return out


def run(solver: str, model, combination: str, *, timeout: float = DEFAULT_TIMEOUT, estimate_errors: bool = True,
        in_process: bool = False, **options) -> ForeignResult:
    """Translate, run and map one combination.  ``options`` go to the solver's ``translate``."""
    from .compare import relative_difference

    mod = _translator(solver)
    tr = mod.translate(model, combination, **options)
    res = tr.results(run_job(tr.job, timeout, in_process))
    if not res.ok or not estimate_errors:
        return res
    convergence = {}
    for name, ctl in tr.controls.items():
        coarse = mod.translate(model, combination, **{**options, name: ctl["coarse"]})
        res2 = coarse.results(run_job(coarse.job, timeout, in_process))
        if not res2.ok:
            convergence[name] = {"coarse": ctl["coarse"], "fine": ctl["fine"], "error": None,
                                 "failure": res2.error}
            continue
        diff = relative_difference(res, res2, model)
        if ctl.get("order"):
            ratio = float(ctl["fine"]) / float(ctl["coarse"])
            err = diff / (ratio ** ctl["order"] - 1.0)
        else:                       # no reliable error law (round-off competes): twice the difference
            err = PENALTY_SAFETY * diff
        res.fidelity.set_error(name, err)
        convergence[name] = {"coarse": ctl["coarse"], "fine": ctl["fine"], "order": ctl.get("order"),
                             "difference": diff, "estimatedError": err}
    extra = getattr(mod, "finalise", None)
    if extra is not None:
        extra(res)
    res.info["convergence"] = convergence
    return res


def convergence_study(solver: str, model, combination: str, control: str, values, *, timeout: float = DEFAULT_TIMEOUT,
                      in_process: bool = False, **options) -> list:
    """Solve with each value of one refinement parameter (``"subdivisions"``, ``"penalty"``,
    ``"stiffness"``) and report the largest relative difference of every run to the LAST one:
    ``[{"value": v, "difference": d, "ok": bool}, ...]`` (the last entry has difference 0)."""
    from .compare import relative_difference

    runs = [run(solver, model, combination, timeout=timeout, estimate_errors=False, in_process=in_process,
                **{**options, control: v}) for v in values]
    ref = runs[-1]
    return [{"value": v, "ok": r.ok,
             "difference": relative_difference(r, ref, model) if (r.ok and ref.ok) else None}
            for v, r in zip(values, runs)]
