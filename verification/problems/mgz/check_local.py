#!/usr/bin/env python3
"""Self-check of the problems in this directory, without the verification pipeline.

For every problem directory: load model.json with benchfea, run ``benchfea.solve.analyze``, write the
results document in memory, evaluate every target of problem.json and compare with its tolerance.  Also
checks that the tolerance of every published target is the one the written policy gives for the number as
printed (docs/VERIFICATION.md section 5: half a unit in the last printed digit, floor 1e-6) -- the printed
form is taken from the target's ``precisionNote`` -- and that a recomputed target carries the floor.

This runs the engine only (tier B targets).  The other tiers (cross-check with foreign solvers,
invariants, convergence) are run by ``benchfea verify verification/problems/mgz/<id>`` or ``--all``.

Usage:  python3 verification/problems/mgz/check_local.py [-v] [--table] [id-substring ...]
        (--table prints the index table of README.md from the current files and results)
Exit code 0: every active problem reproduces all its targets and every tolerance follows the policy.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from benchfea.model import build_model
from benchfea.problems import ModelInvalid
from benchfea.solve import analyze
from benchfea.verify.policy import printed_precision_tolerance
from benchfea.verify.problems import ProblemError, evaluate_targets, load_problem

HERE = Path(__file__).resolve().parent


def policy_findings(spec) -> list:
    out = []
    for k, t in enumerate(spec.targets):
        note = t.get("precisionNote", "")
        m = re.match(r'printed "([^"]+)"(?: x ([0-9.eE+-]+))?', note)
        if m:
            tol = printed_precision_tolerance(m.group(1))
            factor = float(m.group(2)) if m.group(2) else 1.0
            if "abs" in tol:
                tol = {"abs": tol["abs"] * abs(factor)}
            key = next(iter(tol))
            if key not in t["tolerance"] or abs(t["tolerance"][key] - tol[key]) > 1e-12 * tol[key]:
                out.append(f"target {k}: tolerance {t['tolerance']} is not the policy value {tol} for {m.group(1)!r}")
            if abs(abs(t["value"]) - abs(factor) * abs(float(m.group(1).replace(",", "")))) > 1e-9 * abs(t["value"]):
                out.append(f"target {k}: value {t['value']} is not the printed number {m.group(1)!r} x {factor:g}")
        elif note.startswith("evaluated here") or note.startswith("exact rational number"):
            if t["tolerance"] != {"rel": 1e-6}:
                out.append(f"target {k}: a recomputed target must carry the policy floor 1e-6")
        else:
            out.append(f"target {k}: precisionNote does not state the printed form of the number")
    return out


def check(directory: Path, verbose: bool) -> tuple:
    """Returns (status, n_targets, n_missed, worst relative difference / tolerance, lines)."""
    lines = []
    try:
        spec = load_problem(directory)
    except ProblemError as e:
        return "invalid", 0, 0, None, ["  " + x for x in e.errors]
    for f in policy_findings(spec):
        lines.append("  POLICY: " + f)
    try:
        model = build_model(str(spec.model_path))
        results = analyze(model)
    except ModelInvalid as e:
        return "no-results", len(spec.targets), len(spec.targets), None, ["  " + p.format(directory.name)
                                                                           for p in e.problems[:5]]
    doc = results.to_document(model)
    flagged = [cid for cid, c in results.combinations.items() if not (c.invariants.passed and c.converged)]
    if flagged:
        lines.append("  FLAGGED combinations (invariants / convergence): " + ", ".join(flagged))
    rows = evaluate_targets(spec, doc)
    missed, worst = 0, 0.0
    for r in rows:
        tol = r.tolerance
        if r.actual is None:
            missed += 1
            lines.append(f"  MISS  {r.label}: {r.error}")
            continue
        ratio = r.difference / tol["abs"] if "abs" in tol else r.relative / tol["rel"]
        worst = max(worst, ratio)
        if not r.passed:
            missed += 1
        if verbose or not r.passed:
            rel = f"{r.relative:.2e}" if r.relative is not None else "-"
            lines.append(f"  {'ok  ' if r.passed else 'MISS'}  {r.label}: expected {r.expected:.8g}, engine "
                         f"{r.actual:.8g}, diff {r.difference:.2e} (rel {rel}), tolerance {tol}")
    bad_policy = any(x.startswith("  POLICY") for x in lines)
    ok = missed == 0 and not flagged and not bad_policy
    if spec.status == "quarantined":
        status = "quarantined" + (" (all targets met)" if ok else "")
    else:
        status = "passed" if ok else "FAILED"
    return status, len(rows), missed, worst, lines


def table(dirs) -> None:
    """Markdown index (README.md): one row per problem with the measured differences."""
    print("| problem | source location | pages | what it verifies | targets | missed | largest rel. difference "
          "| largest difference / tolerance | status |")
    print("|---|---|---|---|---|---|---|---|---|")
    for d in dirs:
        spec = load_problem(d)
        model = build_model(str(spec.model_path))
        rows = evaluate_targets(spec, analyze(model).to_document(model))
        rel = max((r.relative for r in rows if r.relative is not None), default=0.0)
        ratio = max((r.difference / r.tolerance["abs"] if "abs" in r.tolerance else r.relative / r.tolerance["rel"])
                    for r in rows)
        missed = sum(not r.passed for r in rows)
        print(f"| `{spec.id}` | {spec.source['location']} | {spec.source['pages']} | {spec.title} | {len(rows)} | "
              f"{missed} | {rel:.1e} | {ratio:.2f} | {spec.status} |")


def main(argv) -> int:
    verbose = "-v" in argv
    only = [a for a in argv if not a.startswith("-")]
    dirs = sorted(p.parent for p in HERE.glob("*/problem.json"))
    if only:
        dirs = [d for d in dirs if any(o in d.name for o in only)]
    if "--table" in argv:
        table(dirs)
        return 0
    counts, rc = {}, 0
    for d in dirs:
        status, n, missed, worst, lines = check(d, verbose)
        counts[status.split(" ")[0]] = counts.get(status.split(" ")[0], 0) + 1
        w = f"worst diff/tol {worst:8.3f}" if worst is not None else ""
        print(f"{d.name:44s} {status:28s} targets {n:2d}  missed {missed:2d}  {w}")
        for x in lines:
            print(x)
        if status in ("FAILED", "invalid", "no-results") or any(x.startswith("  POLICY") for x in lines):
            rc = 1
    print("\n" + ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())) + f"  (of {len(dirs)} problems)")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
