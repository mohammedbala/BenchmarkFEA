"""Dataset builder: verified tasks written to disk (docs/BENCHMARK.md section 8, data/DATA_CARD.md).

    from benchfea.bench import dataset
    cells = dataset.plan(seeds=(0, 1, 4))                       # (family, curve ball, seed, difficulty) cells
    manifest = dataset.build("data", cells, workers=4, cache_dir="/tmp/benchfea-cache")

Layout under ``root``:

    tasks/manifest.json                       counts, every task, the quarantine log, versions
    tasks/{train,heldout}/<id>.task.json      PUBLIC: statement, answer format, quantities / criteria
    tasks/{train,heldout}/<id>.model.json     PUBLIC: the input file of the task
    answers/{train,heldout}/<id>.answers.json HIDDEN: change record, reference solution, verification record

A cell is one curve-balled model (``tasks.build_case``); every task type that can be served on it with
a verified reference becomes one task, in the mode ``mode_for`` gives.  A cell without a verified case
and a task type that cannot be served are written to the quarantine log with the reason; nothing
unverified is written as a task.  The build is deterministic, resumable (cases are cached in
``cache_dir``; a second run re-reads them) and parallel (``workers`` processes).
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

from . import curveballs
from .families import FAMILIES
from .families._common import stable_seed
from .splits import split_of
from .tasks import MODES, TASK_TYPES, TASK_VERSION, TaskUnavailable, build_case, make_task

__all__ = ["plan", "build", "mode_for", "level_for", "load_task", "load_answers", "iter_tasks", "SPLIT_DIRS",
           "DATASET_FORMAT"]

SPLIT_DIRS = {"train": "train", "held-out": "heldout"}
DATASET_FORMAT = "benchfea-dataset"


def level_for(family: str, curve_ball: str, seed: int) -> int:
    """Difficulty level of a cell: the levels of the family in rotation over the seed and the position
    of the curve ball in the registry, so that every level a family has is covered."""
    levels = list(FAMILIES[family].LEVELS)
    return levels[(int(seed) + curveballs.available().index(curve_ball)) % len(levels)]


def mode_for(key: dict, task_type: str) -> str:
    """Mode of the task of a cell: a fixed pseudo-random choice (SHA-256 of the cell and task type)."""
    return MODES[stable_seed("dataset-mode", TASK_VERSION, json.dumps(key, sort_keys=True), task_type) % 2]


def plan(seeds=(0, 1, 4), families=None, curve_balls=None) -> list:
    """Every (family, curve ball, seed, difficulty) cell for the given seeds, in a fixed order."""
    out = []
    for family in (families or list(FAMILIES)):
        if family not in FAMILIES:
            raise ValueError(f"unknown family {family!r}")
        for cb in (curve_balls or curveballs.available()):
            curveballs.get(cb)
            for seed in seeds:
                out.append((family, cb, int(seed), level_for(family, cb, seed)))
    return out


def _write(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    os.replace(tmp, path)


_FOREIGN = None


def _cell(args) -> dict:
    """Worker: build one cell and write its tasks; returns the manifest row."""
    global _FOREIGN
    (family, cb, seed, difficulty), root, cache_dir, task_types = args
    t0 = time.perf_counter()
    if _FOREIGN is None:
        from ..verify import ForeignCache
        _FOREIGN = ForeignCache(Path(cache_dir) / "foreign" if cache_dir else None)
    row = {"family": family, "curveBall": cb, "seed": seed, "difficulty": difficulty,
           "split": split_of(family, cb, seed), "status": "unavailable", "tasks": [], "notServed": {}}
    try:
        case = build_case(family, cb, seed, difficulty, cache_dir=cache_dir, foreign_cache=_FOREIGN)
    except Exception as e:                                   # a failure of one cell must not stop the build
        row["skips"] = [{"attempt": 0, "stage": "build", "category": "error",
                         "reason": f"{type(e).__name__}: {str(e)[:400]}"}]
        row["seconds"] = round(time.perf_counter() - t0, 1)
        return row
    row["status"], row["skips"] = case["status"], case["skips"]
    if case["status"] == "verified":
        row["effectClass"] = case["effectClass"]
        row["tiersPassed"] = case["verification"]["tiersPassed"]
        sub = SPLIT_DIRS[row["split"]]
        for tt in task_types:
            mode = mode_for(case["key"], tt)
            try:
                task, answers = make_task(case, tt, mode, foreign_cache=_FOREIGN)
            except TaskUnavailable as e:
                row["notServed"][tt] = e.reason[:600]
                continue
            except Exception as e:
                row["notServed"][tt] = f"error: {type(e).__name__}: {str(e)[:400]}"
                continue
            tid = task["id"]
            public = dict(task)
            model = public.pop("model")
            public["modelFile"] = f"{tid}.model.json"
            _write(Path(root) / "tasks" / sub / f"{tid}.task.json", public)
            _write(Path(root) / "tasks" / sub / f"{tid}.model.json", model)
            _write(Path(root) / "answers" / sub / f"{tid}.answers.json", answers)
            row["tasks"].append({"id": tid, "taskType": tt, "mode": mode})
    row["seconds"] = round(time.perf_counter() - t0, 1)
    return row


def _count(rows, key) -> dict:
    out = {}
    for r in rows:
        k = key(r)
        out[k] = out.get(k, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: str(kv[0])))


def build(root, cells, *, workers: int = 1, cache_dir=None, task_types=TASK_TYPES, clean: bool = False,
          progress=None) -> dict:
    """Build the tasks of ``cells`` under ``root`` and write ``tasks/manifest.json``; returns the manifest.

    ``clean``: remove task / model / answers files under ``root`` that this run did not write."""
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):      # small dense matrices: threads only cost time
        os.environ.setdefault(var, "1")
    root = Path(root)
    jobs = [(tuple(c), str(root), None if cache_dir is None else str(cache_dir), tuple(task_types)) for c in cells]
    rows = []
    if workers > 1:
        import multiprocessing
        with multiprocessing.get_context("spawn").Pool(workers) as pool:
            for row in pool.imap_unordered(_cell, jobs, chunksize=1):
                rows.append(row)
                if progress:
                    progress(row, len(rows), len(jobs))
    else:
        for job in jobs:
            rows.append(_cell(job))
            if progress:
                progress(rows[-1], len(rows), len(jobs))
    order = {tuple(c): i for i, c in enumerate(cells)}
    rows.sort(key=lambda r: order[(r["family"], r["curveBall"], r["seed"], r["difficulty"])])
    manifest = _manifest(rows)
    if clean:
        for split, sub in SPLIT_DIRS.items():
            keep = {t["id"] for t in manifest["tasks"] if t["split"] == split}
            for d, suffixes in ((root / "tasks" / sub, (".task.json", ".model.json")),
                                (root / "answers" / sub, (".answers.json",))):
                for f in (sorted(d.glob("*.json")) if d.is_dir() else []):
                    if f.name.endswith(suffixes) and f.name.split(".")[0] not in keep:
                        f.unlink()
    _write(root / "tasks" / "manifest.json", manifest)
    return manifest


def _manifest(rows: list) -> dict:
    from ..verify import MEANING, software_versions
    tasks, quarantine = [], []
    for r in rows:
        for t in r["tasks"]:
            tasks.append({"id": t["id"], "split": r["split"], "taskType": t["taskType"], "mode": t["mode"],
                          "family": r["family"], "difficulty": r["difficulty"], "effectClass": r["effectClass"],
                          "tiersPassed": r["tiersPassed"],
                          "file": f"{SPLIT_DIRS[r['split']]}/{t['id']}.task.json"})
        cell = {k: r[k] for k in ("family", "curveBall", "seed", "difficulty", "split")}
        for s in r.get("skips", []):
            quarantine.append({**cell, "what": "attempt", **s,
                               "outcome": "replaced by the next seed of the chain" if r["status"] == "verified"
                               else "cell not served"})
        for tt, why in r["notServed"].items():
            quarantine.append({**cell, "what": "task type", "taskType": tt, "reason": why})
    verified = [r for r in rows if r["status"] == "verified"]
    pairs = {}
    for r in rows:
        p = pairs.setdefault(f"{r['family']}|{r['curveBall']}", {"cells": 0, "verified": 0, "tasks": 0, "reasons": {}})
        p["cells"] += 1
        p["verified"] += r["status"] == "verified"
        p["tasks"] += len(r["tasks"])
        if r["status"] != "verified":
            last = (r.get("skips") or [{}])[-1]
            cat = f"{last.get('stage')}:{last.get('category')}"
            p["reasons"][cat] = p["reasons"].get(cat, 0) + 1

    def by(rows_, key):
        out = {}
        for r in rows_:
            o = out.setdefault(key(r), {"cells": 0, "verified": 0, "tasks": 0, "notVerified": {}})
            o["cells"] += 1
            o["verified"] += r["status"] == "verified"
            o["tasks"] += len(r["tasks"])
            if r["status"] != "verified":
                last = (r.get("skips") or [{}])[-1]
                cat = f"{last.get('stage')}:{last.get('category')}"
                o["notVerified"][cat] = o["notVerified"].get(cat, 0) + 1
        return dict(sorted(out.items()))

    return {
        "format": DATASET_FORMAT, "formatVersion": "1.0.0", "meaning": MEANING,
        "versions": {**software_versions(), "task": TASK_VERSION},
        "counts": {"cells": len(rows), "cellsVerified": len(verified), "tasks": len(tasks),
                   "bySplit": _count(tasks, lambda t: t["split"]),
                   "byTaskType": _count(tasks, lambda t: t["taskType"]),
                   "byMode": _count(tasks, lambda t: t["mode"]),
                   "byFamily": _count(tasks, lambda t: t["family"]),
                   "byDifficulty": _count(tasks, lambda t: t["difficulty"]),
                   "byEffectClass": _count(tasks, lambda t: t["effectClass"]),
                   "bySplitAndTaskType": _count(tasks, lambda t: f"{t['split']}|{t['taskType']}")},
        "cellsByFamily": by(rows, lambda r: r["family"]),
        "cellsByCurveBall": by(rows, lambda r: r["curveBall"]),
        "cellsByFamilyAndCurveBall": dict(sorted(pairs.items())),
        "taskTypesNotServed": _count([(tt, why.split(":")[0][:90]) for r in rows for tt, why in r["notServed"].items()],
                                     lambda x: f"{x[0]}: {x[1]}"),
        "tasks": tasks,
        "cells": [{k: r.get(k) for k in ("family", "curveBall", "seed", "difficulty", "split", "status",
                                         "effectClass")} | {"tasks": [t["id"] for t in r["tasks"]]} for r in rows],
        "quarantineLog": quarantine,
    }


# ------------------------------------------------------------------------------------------ reading
def load_task(path) -> dict:
    """A task file with its model embedded (``task["model"]``)."""
    path = Path(path)
    task = json.loads(path.read_text(encoding="utf-8"))
    if "model" not in task:
        task["model"] = json.loads((path.parent / task["modelFile"]).read_text(encoding="utf-8"))
    return task


def load_answers(task_path) -> dict:
    """The hidden answers of a task file of the dataset layout (``answers/<split>/<id>.answers.json``)."""
    p = Path(task_path)
    name = p.name.replace(".task.json", ".answers.json")
    return json.loads((p.parent.parent.parent / "answers" / p.parent.name / name).read_text(encoding="utf-8"))


def iter_tasks(root, split: str = "train"):
    """Paths of the task files of a split, sorted."""
    return sorted((Path(root) / "tasks" / SPLIT_DIRS[split]).glob("*.task.json"))
