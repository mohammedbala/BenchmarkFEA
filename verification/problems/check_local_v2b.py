#!/usr/bin/env python3
"""Self-check of the problems under aisc-commentary/, ziemian/, aisc-de/, staad/ and ansys/, without the
verification pipeline.

For every problem directory: load model.json with benchfea, run ``benchfea.solve.analyze``, evaluate every
target of problem.json on the results document and compare with its tolerance.  Also checks that every
tolerance is the one the written policy gives (docs/VERIFICATION.md section 5), from the target's
``precisionNote``:

* ``printed "<number>" [x factor]``           half a unit in the last printed digit (floor 1e-6);
* ``stated tolerance <rel> on printed "..."``  an acceptance tolerance stated by the source itself; only
                                               the values of AISC 360-22 Commentary C2.1 are accepted
                                               (0.03 on a moment, 0.05 on a deflection);
* ``evaluated here ...``                       the policy floor 1e-6 (absolute 1e-9 for a zero target).

This runs the engine only (the targets).  The other tiers (cross-check with foreign solvers, invariants,
convergence) are run by ``benchfea verify verification/problems/<source>/<id>`` or ``--all``.

Usage:  python3 verification/problems/check_local_v2b.py [-v] [--table] [substring of "<source>/<id>" ...]
        (--table prints the index table of the README.md files from the current files and results)
Exit code 0: every ACTIVE problem reproduces all its targets and every tolerance follows the policy.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from benchfea.model import build_model
from benchfea.problems import ModelInvalid
from benchfea.solve import analyze
from benchfea.verify.policy import printed_precision_tolerance
from benchfea.verify.problems import ProblemError, evaluate_targets, load_problem

HERE = Path(__file__).resolve().parent
SOURCES = ("aisc-commentary", "ziemian", "aisc-de", "staad", "ansys")
STATED = (0.03, 0.05)                    # AISC 360-22 Commentary C2.1, p. 16.1-337


def policy_findings(spec) -> list:
    out = []
    for k, t in enumerate(spec.targets):
        note, tol = t.get("precisionNote", ""), t["tolerance"]
        m = re.match(r'printed "([^"]+)"(?: x ([0-9.eE+-]+))?', note)
        s = re.match(r'stated tolerance ([0-9.eE+-]+) on printed "([^"]+)"(?: x ([0-9.eE+-]+))?', note)
        if m:
            want = printed_precision_tolerance(m.group(1))
            factor = float(m.group(2)) if m.group(2) else 1.0
            if "abs" in want:
                want = {"abs": want["abs"] * abs(factor)}
            key = next(iter(want))
            if key not in tol or abs(tol[key] - want[key]) > 1e-12 * want[key]:
                out.append(f"target {k}: tolerance {tol} is not the policy value {want} for {m.group(1)!r}")
            printed = abs(factor) * abs(float(m.group(1).replace(",", "")))
            if abs(abs(t["value"]) - printed) > 1e-9 * max(abs(t["value"]), 1e-300):
                out.append(f"target {k}: value {t['value']} is not the printed number {m.group(1)!r} x {factor:g}")
        elif s:
            rel = float(s.group(1))
            if spec.directory.parent.name != "aisc-commentary" or rel not in STATED or tol != {"rel": rel}:
                out.append(f"target {k}: a stated tolerance is accepted only for AISC 360-22 Commentary C2.1 "
                           f"(0.03 / 0.05), found {tol}")
            factor = float(s.group(3)) if s.group(3) else 1.0
            printed = abs(factor) * abs(float(s.group(2).replace(",", "")))
            if abs(abs(t["value"]) - printed) > 1e-9 * max(abs(t["value"]), 1e-300):
                out.append(f"target {k}: value {t['value']} is not the printed number {s.group(2)!r} x {factor:g}")
        elif note.startswith("evaluated here"):
            if tol != ({"rel": 1e-6} if t["value"] != 0.0 else {"abs": 1e-9}):
                out.append(f"target {k}: a recomputed target must carry the policy floor 1e-6")
        else:
            out.append(f"target {k}: precisionNote does not state the printed form of the number")
    return out


def run(directory: Path):
    spec = load_problem(directory)
    model = build_model(str(spec.model_path))
    results = analyze(model)
    return spec, results, evaluate_targets(spec, results.to_document(model))


def ratio(r):
    if r.actual is None:
        return float("inf")
    return r.difference / r.tolerance["abs"] if "abs" in r.tolerance else r.relative / r.tolerance["rel"]


def check(directory: Path, verbose: bool) -> tuple:
    lines = []
    try:
        spec = load_problem(directory)
    except ProblemError as e:
        return "invalid", 0, 0, None, ["  " + x for x in e.errors]
    lines += ["  POLICY: " + f for f in policy_findings(spec)]
    try:
        spec, results, rows = run(directory)
    except ModelInvalid as e:
        return ("quarantined (no results)" if spec.status == "quarantined" else "no-results"), len(spec.targets), \
            len(spec.targets), None, lines + ["  " + p.format(directory.name) for p in e.problems[:3]]
    flagged = [cid for cid, c in results.combinations.items() if not (c.invariants.passed and c.converged)]
    if flagged:
        lines.append("  FLAGGED combinations (invariants / convergence): " + ", ".join(flagged))
    missed, worst = 0, 0.0
    for r in rows:
        if r.actual is None:
            missed += 1
            lines.append(f"  MISS  {r.label}: {r.error}")
            continue
        worst = max(worst, ratio(r))
        missed += not r.passed
        if verbose or not r.passed:
            rel = f"{r.relative:.2e}" if r.relative is not None else "-"
            lines.append(f"  {'ok  ' if r.passed else 'MISS'}  {r.label}: expected {r.expected:.8g}, engine "
                         f"{r.actual:.8g}, diff {r.difference:.2e} (rel {rel}), tolerance {r.tolerance}")
    ok = missed == 0 and not flagged and not any(x.startswith("  POLICY") for x in lines)
    if spec.status == "quarantined":
        status = "quarantined" + (" (all targets met)" if ok else "")
    else:
        status = "passed" if ok else "FAILED"
    return status, len(rows), missed, worst, lines


def table(dirs) -> None:
    print("| problem | source location | what it verifies | basis | targets | missed | largest rel. difference "
          "| largest difference / tolerance | status |")
    print("|---|---|---|---|---|---|---|---|---|")
    for d in dirs:
        spec = load_problem(d)
        where = spec.source["location"].split(" (http")[0]           # the URL is in problem.json
        try:
            _, _, rows = run(d)
        except ModelInvalid:
            print(f"| `{spec.id}` | {where} | {spec.title} | {spec.basis} | {len(spec.targets)} | "
                  f"- | - | - | {spec.status} (the engine refuses the model) |")
            continue
        rel = max((r.relative for r in rows if r.relative is not None), default=0.0)
        print(f"| `{spec.id}` | {where} | {spec.title} | {spec.basis} | {len(rows)} | "
              f"{sum(not r.passed for r in rows)} | {rel:.1e} | {max(ratio(r) for r in rows):.2f} | {spec.status} |")


def main(argv) -> int:
    verbose = "-v" in argv
    only = [a for a in argv if not a.startswith("-")]
    dirs = sorted(p.parent for s in SOURCES for p in (HERE / s).glob("*/problem.json"))
    if only:
        dirs = [d for d in dirs if any(o in f"{d.parent.name}/{d.name}" for o in only)]
    if "--table" in argv:
        table(dirs)
        return 0
    counts, rc = {}, 0
    for d in dirs:
        status, n, missed, worst, lines = check(d, verbose)
        counts[status.split(" ")[0]] = counts.get(status.split(" ")[0], 0) + 1
        w = f"worst diff/tol {worst:9.3f}" if worst is not None else ""
        print(f"{d.parent.name + '/' + d.name:62s} {status:28s} targets {n:2d}  missed {missed:2d}  {w}")
        for x in lines:
            print(x)
        if status in ("FAILED", "invalid", "no-results") or any(x.startswith("  POLICY") for x in lines):
            rc = 1
    print("\n" + ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())) + f"  (of {len(dirs)} problems)")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
