#!/usr/bin/env python3
"""After a dataset build: put every task through the normal scoring path.

    python3 tools/check_dataset_scores.py data --workers 4 [--cache DIR] [--limit N]

For every task of the dataset (``benchfea.bench.scoring.score_answer``; repair and design tasks run
the verification pipeline, so this takes a while):

* the answer built from the hidden reference (``scoring.reference_answer``) must score 1.0;
* for a repair task, the model as given, handed in as the "repair", must score less than the
  reference answer does.

The summary is printed and stored under the key ``scores`` of ``<root>/CHECKS.json`` together with
the SHA-256 of the manifest it belongs to (``tools/build_dataset_docs.py`` quotes it only when that
hash is the one of the manifest on disk).  With ``--limit`` nothing is stored.  Exit code 0 when
nothing is found, 1 otherwise.  This is a computation; it does not say that anybody reviewed,
approved or certified a task.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_FOREIGN = None


def one(args) -> dict:
    """Worker: the scores of one task file."""
    global _FOREIGN
    path, cache = args
    from benchfea.bench import dataset, scoring
    if _FOREIGN is None:
        from benchfea.verify import ForeignCache
        _FOREIGN = ForeignCache(Path(cache) / "foreign" if cache else None)
    task, answers = dataset.load_task(path), dataset.load_answers(path)
    t0 = time.perf_counter()
    out = {"id": task["id"], "taskType": task["taskType"], "given": None}
    try:
        out["reference"] = scoring.score_answer(task, answers, scoring.reference_answer(task, answers),
                                                foreign_cache=_FOREIGN)["score"]
        if task["taskType"] == "repair":
            out["given"] = scoring.score_answer(task, answers, {"model": task["model"]},
                                                foreign_cache=_FOREIGN)["score"]
    except Exception as e:                                    # reported as a finding, never swallowed
        out["error"] = f"{type(e).__name__}: {str(e)[:300]}"
    out["seconds"] = round(time.perf_counter() - t0, 1)
    return out


def manifest_hash(root: Path) -> str:
    return hashlib.sha256((root / "tasks" / "manifest.json").read_bytes()).hexdigest()


def update_checks(root: Path, key: str, summary: dict) -> None:
    """Store ``summary`` under ``key`` of ``<root>/CHECKS.json`` with the hash of the manifest."""
    path = root / "CHECKS.json"
    data = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    data[key] = {"manifestSha256": manifest_hash(root), **summary}
    path.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--cache", default=None, help="directory of the foreign-solve cache of the build")
    ap.add_argument("--limit", type=int, default=None, help="only the first N tasks (nothing is stored)")
    args = ap.parse_args(argv)
    root = Path(args.root)
    manifest = json.loads((root / "tasks" / "manifest.json").read_text(encoding="utf-8"))
    files = [str(root / "tasks" / t["file"]) for t in manifest["tasks"]][: args.limit]
    jobs = [(f, args.cache) for f in files]
    rows = []
    if args.workers > 1:
        import multiprocessing
        with multiprocessing.get_context("spawn").Pool(args.workers) as pool:
            for k, r in enumerate(pool.imap_unordered(one, jobs, chunksize=2)):
                rows.append(r)
                if (k + 1) % 50 == 0:
                    print(f"  {k + 1} / {len(jobs)}", file=sys.stderr, flush=True)
    else:
        rows = [one(j) for j in jobs]
    rows.sort(key=lambda r: r["id"])
    findings = []
    for r in rows:
        if "error" in r:
            findings.append(f"{r['id']} ({r['taskType']}): {r['error']}")
            continue
        if r["reference"] != 1.0:
            findings.append(f"{r['id']} ({r['taskType']}): the reference answer scores {r['reference']}")
        if r["given"] is not None and not r["given"] < r["reference"]:
            findings.append(f"{r['id']}: the model as given scores {r['given']}, the reference answer {r['reference']}")
    by_type = {}
    for r in rows:
        o = by_type.setdefault(r["taskType"], {"tasks": 0, "referenceScores1": 0})
        o["tasks"] += 1
        o["referenceScores1"] += r.get("reference") == 1.0
    given = [r["given"] for r in rows if r.get("given") is not None]
    summary = {
        "tasks": len(rows), "referenceAnswerScores1": sum(r.get("reference") == 1.0 for r in rows),
        "byTaskType": dict(sorted(by_type.items())),
        "repair": {"tasks": len(given), "givenModelScoresBelowReference": sum(g < 1.0 for g in given),
                   "givenModelScoresZero": sum(g == 0.0 for g in given),
                   "largestScoreOfAGivenModel": max(given) if given else None},
        "findings": findings,
    }
    print(json.dumps(summary, indent=1))
    if args.limit is None:
        update_checks(root, "scores", summary)
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
