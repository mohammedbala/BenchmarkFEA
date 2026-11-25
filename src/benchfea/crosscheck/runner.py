"""Run a foreign solver on a translated model, in a SUBPROCESS with a timeout.

OpenSees terminates the interpreter on some input errors (measured: a ``vecxz`` parallel to the member
axis, a system without free DOFs), so a foreign solver never runs in the caller's process unless
``in_process=True`` is asked for.  The job goes to the worker as a JSON file and the raw results come
back as a JSON file (stdout is not used: the libraries print to it).

    result = run("opensees", model, "LC1")          # -> ForeignResult
    python -m benchfea.crosscheck.worker job.json out.json     # the worker process

Error of an approximate feature (``Translation.controls``), from runs of the FOREIGN solver only:

* discretisation parameter (elements per member; ``levels`` = n/4, n/2, n): the model is solved at
  the three levels and the error of the finest run is bounded by :func:`estimate_from_levels` from
  the two level-to-level differences.  No order of convergence is assumed beyond "at least
  ``orders[0]``" for the levels that were not run; a sequence that does not contract gives no
  estimate, and a solver without an estimate is not eligible for a cross-check.
* penalty stiffness / stiff emulation member (``order = None``): round-off competes with the
  approximation, no error law is used: twice the difference ``|fine - coarse|`` of two levels.

The estimate is written into the fidelity record; the levels, differences and the observed order go
to ``info["convergence"]``.
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

__all__ = ["run", "run_job", "available", "convergence_study", "estimate_from_levels", "SOLVERS",
           "DEFAULT_TIMEOUT", "ROUNDOFF_FLOOR"]

SOLVERS = {"opensees": ("benchfea.crosscheck.opensees", "openseespy"),
           "pynite": ("benchfea.crosscheck.pynite", "Pynite")}
DEFAULT_TIMEOUT = 300.0
PENALTY_SAFETY = 2.0
ROUNDOFF_FLOOR = 1.0e-9       # level-to-level differences below this are round-off of the foreign solves (three
                              # decades under the comparison tolerance 1e-6): no order can be observed from them


def estimate_from_levels(d1: float, d2: float | None, ratio: float = 2.0, orders=(1.0, 2.0),
                         floor: float = ROUNDOFF_FLOOR) -> dict:
    """Bound of the discretisation error of the finest of three runs at n, ratio*n, ratio^2*n.

    ``d1`` = |r(n) - r(ratio n)|, ``d2`` = |r(ratio n) - r(ratio^2 n)| (largest relative difference over
    all compared quantities).  With r(inf) the limit,

        |r(fine) - r(inf)|  <=  sum of all further level-to-level differences
                            <=  d2 (rho + rho^2 + ...) = d2 rho / (1 - rho)

    if every further difference contracts at least by rho.  Used here:

    * ``rho = max(d2/d1, ratio**-p_min)``: the observed contraction, but never better than the lowest
      order ``p_min = orders[0]`` the method is known to have -- an error a/n^2 + b/n contracts by
      nearly 1/4 on coarse levels and by 1/2 in the limit, so the observed contraction must not be
      extrapolated when it is better than the guaranteed one;
    * ``d2`` is replaced by ``max(d2, d1 * ratio**-p_max)``: a contraction faster than the highest
      order ``p_max = orders[1]`` the method can have is a cancellation between two levels, not
      convergence;
    * ``d2 >= d1`` (no contraction): the levels are not in the asymptotic range -> ``estimate`` is
      ``None`` (the solver is then not eligible; never a pass);
    * both differences below ``floor``: round-off, the estimate is the larger of the two;
    * ``d2 is None`` (only two levels were run): no order can be observed; ``d1`` is the last
      difference and the bound uses ``p_min``.

    For ``p_min = 1`` and ratio 2 the bound is the last difference itself whenever the observed order
    is between 1 and ``p_max``; it exceeds a pure second-order error by the factor 3 and a pure
    first-order error by the factor 1.  Returns ``{"estimate", "observedOrder", "contraction",
    "status"}`` with status ``"asymptotic"`` | ``"round-off"`` | ``"two-levels"`` | ``"not-asymptotic"``.
    """
    import math

    p_min, p_max = float(orders[0]), float(orders[1])
    rho_hi, rho_lo = ratio ** -p_min, ratio ** -p_max
    if d2 is None:
        return {"estimate": d1 * rho_hi / (1.0 - rho_hi), "observedOrder": None, "contraction": None,
                "status": "two-levels"}
    if d1 <= floor and d2 <= floor:
        return {"estimate": max(d1, d2), "observedOrder": None, "contraction": None, "status": "round-off"}
    if not d2 < d1:
        return {"estimate": None, "observedOrder": None, "contraction": (d2 / d1 if d1 > 0.0 else None),
                "status": "not-asymptotic"}
    rho = d2 / d1
    order = math.log(1.0 / rho) / math.log(ratio) if rho > 0.0 else None
    used = max(rho, rho_hi)
    return {"estimate": max(d2, d1 * rho_lo) * used / (1.0 - used), "observedOrder": order, "contraction": rho,
            "status": "asymptotic"}


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

    def at(name, value):
        other = mod.translate(model, combination, **{**options, name: value})
        return other.results(run_job(other.job, timeout, in_process))

    for name, ctl in tr.controls.items():
        if "levels" not in ctl:     # penalty / emulation stiffness: no error law (round-off competes)
            res2 = at(name, ctl["coarse"])
            if not res2.ok:
                convergence[name] = {"coarse": ctl["coarse"], "fine": ctl["fine"], "error": None,
                                     "failure": res2.error}
                continue
            diff = relative_difference(res, res2, model)
            err = PENALTY_SAFETY * diff
            res.fidelity.set_error(name, err)
            convergence[name] = {"coarse": ctl["coarse"], "fine": ctl["fine"], "order": None,
                                 "difference": diff, "estimatedError": err}
            continue
        levels = list(ctl["levels"])
        rec = {"levels": levels, "fine": ctl["fine"], "orders": list(ctl["orders"]), "differences": [],
               "estimatedError": None}
        convergence[name] = rec
        if len(levels) < 2:
            rec["status"] = "one-level"
            continue
        runs = [at(name, v) for v in levels[:-1]] + [res]
        bad = [r.error for r in runs if not r.ok]
        if bad:
            rec.update(status="failed", failure=bad[0])
            continue
        diffs = [relative_difference(runs[k + 1], runs[k], model) for k in range(len(runs) - 1)]
        est = estimate_from_levels(diffs[-2] if len(diffs) > 1 else diffs[-1], diffs[-1] if len(diffs) > 1 else None,
                                   float(levels[-1]) / float(levels[-2]), ctl["orders"])
        rec.update(differences=diffs, estimatedError=est["estimate"], observedOrder=est["observedOrder"],
                   contraction=est["contraction"], status=est["status"])
        if est["estimate"] is not None:
            res.fidelity.set_error(name, est["estimate"])
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
