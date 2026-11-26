#!/usr/bin/env python3
"""Rebuild a sample of the cells of the committed dataset from nothing and compare the files byte for byte.

    python3 tools/compare_dataset_sample.py                       # the default sample (about 2 minutes)
    python3 tools/compare_dataset_sample.py --families continuous_beam,portal_frame \\
            --curve-balls section_lighter,support_removed --seeds 0,4 --workers 2

The cells (family x curve ball x seed) of the sample are built into a temporary directory with
``benchfea bench make-dataset`` -- generation, curve ball and the whole verification pipeline run
again, without a cache -- and then

* every cell of the sample must have the status and the task ids it has in ``data/tasks/manifest.json``;
* every task, model and answers file written must equal the committed file byte for byte.

The answers files record the versions of python, numpy, scipy and the two foreign solvers they were
built with.  When this run uses other versions, a file that differs ONLY in those recorded versions
is counted separately ("equal apart from the recorded versions") and is not a difference; with
``--strict`` it is.  Any other difference -- a number in its last digit included -- means that the
committed dataset is not what the code of the repository builds in this environment, and is listed.
Exit code 0 when there is no difference, 1 otherwise.  This is a comparison of two runs of the same code; it adds no
independent evidence about the reference solutions (that is the job of ``tools/audit_dataset.py``).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = {"families": "continuous_beam", "curve_balls": "section_lighter,support_removed,load_sign_flipped",
           "seeds": "0,4"}


def without_versions(x):
    if isinstance(x, dict):
        return {k: without_versions(v) for k, v in x.items() if k != "versions"}
    if isinstance(x, list):
        return [without_versions(v) for v in x]
    return x


def cell_key(c):
    return (c["family"], c["curveBall"], c["seed"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", default=str(ROOT / "data"), help="the committed dataset (default: data)")
    ap.add_argument("--families", default=DEFAULT["families"], help="families of the sample, comma separated")
    ap.add_argument("--curve-balls", default=DEFAULT["curve_balls"], help="curve balls of the sample, comma separated")
    ap.add_argument("--seeds", default=DEFAULT["seeds"], help="instance seeds of the sample, comma separated")
    ap.add_argument("--workers", type=int, default=1, help="parallel processes")
    ap.add_argument("--strict", action="store_true",
                    help="a file that differs only in the recorded software versions is a difference too")
    args = ap.parse_args(argv)

    data = Path(args.data)
    committed = json.loads((data / "tasks" / "manifest.json").read_text(encoding="utf-8"))
    env = dict(os.environ)
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        env.setdefault(var, "1")
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    problems = []
    with tempfile.TemporaryDirectory(prefix="benchfea-sample-") as tmp:
        cmd = [sys.executable, "-m", "benchfea.cli", "bench", "make-dataset", "--out", tmp, "--seeds", args.seeds,
               "--families", args.families, "--curve-balls", args.curve_balls, "--workers", str(args.workers),
               "--quiet"]
        p = subprocess.run(cmd, env=env, cwd=ROOT, capture_output=True, text=True)
        if p.returncode != 0:
            print(p.stdout + p.stderr)
            print("the sample could not be built")
            return 1
        sample = json.loads((Path(tmp) / "tasks" / "manifest.json").read_text(encoding="utf-8"))
        if sample["versions"] != committed["versions"]:
            print("versions differ -- committed:", committed["versions"])
            print("                   this run: ", sample["versions"])
        there = {cell_key(c): c for c in committed["cells"]}
        for c in sample["cells"]:
            ref = there.get(cell_key(c))
            if ref is None:
                problems.append(f"cell {cell_key(c)}: not in the committed dataset")
            elif (c["status"], c["tasks"]) != (ref["status"], ref["tasks"]):
                problems.append(f"cell {cell_key(c)}: rebuilt {c['status']} {c['tasks']}, committed "
                                f"{ref['status']} {ref['tasks']}")
        files = equal = versions_only = 0
        for path in sorted(Path(tmp).rglob("*.json")):
            rel = path.relative_to(tmp)
            if rel.name == "manifest.json":
                continue
            files += 1
            other = data / rel
            if not other.is_file():
                problems.append(f"{rel}: not in the committed dataset")
            elif other.read_bytes() == path.read_bytes():
                equal += 1
            elif not args.strict and without_versions(json.loads(other.read_text(encoding="utf-8"))) == \
                    without_versions(json.loads(path.read_text(encoding="utf-8"))):
                versions_only += 1
            else:
                problems.append(f"{rel}: differs from the committed file")
        tasks = sum(len(c["tasks"]) for c in sample["cells"])
    for line in problems:
        print(line)
    print(f"{len(sample['cells'])} cell(s) rebuilt, {tasks} task(s), {files} file(s): {equal} equal to the committed "
          f"file byte for byte, {versions_only} equal apart from the recorded versions, {len(problems)} "
          f"difference(s)")
    if not files:
        print("the sample holds no task: nothing was compared")
        return 1
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
