#!/usr/bin/env python3
"""Audit a dataset built by ``benchfea bench make-dataset`` with means that are NOT the verification
pipeline that produced it.

    python3 tools/audit_dataset.py data                 # every check, every task
    python3 tools/audit_dataset.py data --limit 20      # the first 20 distinct models only
    python3 tools/audit_dataset.py data --json out.json

What it checks, and with what:

split        every task lies in the split that ``src/benchfea/bench/splits.json`` assigns to its
             (family, curve ball, seed) -- the rule is re-evaluated here from the file, not through
             ``benchfea.bench.splits``; no curve-balled model and no base model (the model before
             the curve ball, by the hash in the change record) occurs in both splits.
change       the change record applied to nothing: the model file hashes to ``perturbedHash``.
references   every distinct STABLE model is solved again by the engine and every reference value of
             its predict tasks is recomputed; the same model is solved by the independent tableau
             solver of ``tests/independent`` (a different method: no stiffness matrix; section
             properties from ``benchfea.sections``) and compared with the engine on node displacements,
             reactions, end forces and station values, every combination.  A difference in a
             first-order combination is a finding.  For a second-order combination the tableau takes
             the axial force piecewise constant along a member (its documented limit), so the
             difference is reported separately and is not a finding.
mechanisms   every distinct UNSTABLE model: the number of independent mechanisms of the reference
             diagnosis against the tableau (null space of its matrix + the degrees of freedom it finds
             without any stiffness at nodes that carry frame members).
scoring      the answer built from the hidden reference scores 1.0, a wrong one less than 1.0
             (predict and diagnose tasks; repair and design tasks run the pipeline and are left to
             ``tests/test_dataset.py`` and the build's own check).
leaks        the first observation and info of the environment for a sample of tasks contain neither
             the name of the curve ball, nor the target of the change, nor a hash of the record.

Without ``--limit`` the summary is also stored under the key ``audit`` of ``<root>/CHECKS.json`` with
the SHA-256 of the manifest it belongs to.  Exit code 0 when nothing is found, 1 otherwise.  This is an audit by computation; it does not say that
anybody reviewed, approved or certified a task.
"""
from __future__ import annotations

import collections
import copy
import glob
import hashlib
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests" / "independent"))
sys.path.insert(0, str(ROOT / "tools"))

TABLEAU_TOLERANCE = 1.0e-8          # engine against the independent tableau solver (first order; round-off)
SINGULAR = 1.0e-9                   # singular value of the equilibrated tableau counted as zero


def explicit_sections(raw: dict):
    """The document with every section given by its properties (the tableau solver reads properties
    only), and the model."""
    from benchfea.model import build_model
    model = build_model(copy.deepcopy(raw))
    u = model.units
    doc = copy.deepcopy(raw)
    for sid, s in model.sections.items():
        new = {"type": "properties", "A": u.from_internal(s.A, "length^2"), "Iy": u.from_internal(s.Iy, "length^4"),
               "Iz": u.from_internal(s.Iz, "length^4"), "J": u.from_internal(s.J, "length^4")}
        for k in ("Asy", "Asz"):
            if getattr(s, k, None):
                new[k] = u.from_internal(getattr(s, k), "length^2")
        doc["sections"][sid] = new
    doc.pop("design", None)
    return doc, model


def audit_stable(model_path: str, tasks: list) -> dict:
    """Engine re-solve + tableau comparison of one stable model; ``tasks``: [(task, answers)]."""
    import harness
    from benchfea import design as design_module
    from benchfea.bench.tasks import extract
    from benchfea.solve import analyze
    out = {"model": model_path, "findings": [], "combinations": 0, "tableauWorst": 0.0, "secondOrder": 0,
           "referenceWorst": 0.0, "quantities": 0, "tableauWorstSecondOrder": 0.0}
    raw = json.loads(Path(model_path).read_text(encoding="utf-8"))
    doc, model = explicit_sections(raw)
    results = analyze(model, strict=True)
    rdoc = results.to_document(model)
    ddoc = None
    for task, ans in tasks:
        if task["taskType"] != "predict":
            continue
        qs = {q["id"]: q for q in task["quantities"]}
        for qid, rv in ans["reference"]["values"].items():
            q = qs[qid]
            if q["kind"] == "designRatio" and ddoc is None:
                ddoc = design_module.check(model, results).to_document(model)
            value = extract(q, rdoc, ddoc)
            d = abs(value - rv["value"]) / max(abs(rv["value"]), rv["floor"])
            out["referenceWorst"] = max(out["referenceWorst"], d)
            out["quantities"] += 1
            if d > 1e-9:
                out["findings"].append(f"{task['id']} {qid}: reference {rv['value']!r}, fresh engine solve {value!r}")
    for cid, c in model.combinations.items():
        second = c.analysis_type == "secondOrder"
        out["secondOrder"] += second
        try:
            diff = harness.compare(copy.deepcopy(doc), cid, second_order=second, doc=rdoc)
        except Exception as e:                                   # the tableau cannot take the model
            out["findings"].append(f"combination {cid}: tableau solver failed ({type(e).__name__}: {str(e)[:120]})")
            continue
        worst = max(v for k, v in diff.items() if not k.startswith("_"))
        out["combinations"] += 1
        if second:
            out["tableauWorstSecondOrder"] = max(out["tableauWorstSecondOrder"], float(worst))
            continue
        out["tableauWorst"] = max(out["tableauWorst"], float(worst))
        if worst > TABLEAU_TOLERANCE:
            out["findings"].append(f"combination {cid}: engine and tableau differ by {worst:.2e}")
    return out


def audit_mechanism(model_path: str, n_reference: int) -> dict:
    """Number of independent mechanisms of an unstable model by the tableau solver."""
    import numpy as np
    from tableau import Tableau
    raw = json.loads(Path(model_path).read_text(encoding="utf-8"))
    doc, _ = explicit_sections(raw)
    for c in doc["combinations"].values():
        c.pop("analysisType", None)
    doc.get("analysis", {}).pop("type", None)
    t = Tableau(doc)
    t.set_combination(next(iter(doc["combinations"])))
    M, _ = t.assemble()
    A, _, _ = t._equilibrate(M)
    s = np.linalg.svd(A, compute_uv=False)
    null = int((s < SINGULAR * s[0]).sum())
    frame_nodes = {n for m in doc["members"].values() if m.get("type", "frame") == "frame" for n in (m["i"], m["j"])}
    n_u = 6 * len(t.node_ids)
    # DOFs the tableau took out because nothing gives them stiffness, at nodes that carry frame members
    idle = sum(1 for g in range(n_u) if M[g, g] == 1.0 and np.count_nonzero(M[g, :]) == 1
               and np.count_nonzero(M[:, g]) == 1 and not t.fixed[g] and t.node_ids[g // 6] in frame_nodes)
    out = {"model": model_path, "reference": n_reference, "nullSpace": null, "withoutStiffness": idle, "findings": []}
    if null + idle != n_reference:
        out["findings"].append(f"reference diagnosis has {n_reference} mechanism(s), the tableau finds "
                               f"{null} + {idle}")
    return out


def _held_out(rule: dict, family: str, curve_ball: str, seed: int) -> bool:
    return [family, curve_ball] in rule["heldOutCombinations"] or \
        seed % rule["seedRule"]["modulus"] in rule["seedRule"]["heldOutResidues"]


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0].startswith("-"):
        print(__doc__)
        return 2
    root = Path(argv[0])
    limit = int(argv[argv.index("--limit") + 1]) if "--limit" in argv else None
    from benchfea.bench import dataset
    from benchfea.bench.scoring import reference_answer, score_answer
    from benchfea.verify.common import document_hash
    rule = json.loads((ROOT / "src" / "benchfea" / "bench" / "splits.json").read_text(encoding="utf-8"))
    findings, summary = [], collections.OrderedDict()
    tasks = []
    for split, sub in (("train", "train"), ("held-out", "heldout")):
        for p in sorted(glob.glob(str(root / "tasks" / sub / "*.task.json"))):
            task = dataset.load_task(p)
            ans = dataset.load_answers(p)
            tasks.append((split, p, task, ans))
    summary["tasks"] = len(tasks)
    # ---------------------------------------------------------------- split, hashes
    by_pert, by_base = {}, {}
    for split, p, task, ans in tasks:
        key = ans["key"]
        want = "held-out" if _held_out(rule, key["family"], key["curveBall"], key["seed"]) else "train"
        if not (want == split == task["split"]):
            findings.append(f"{task['id']}: split {task['split']} in directory {split}, rule says {want}")
        rec = ans["changeRecord"]
        if document_hash(task["model"]) != rec["perturbedHash"]:
            findings.append(f"{task['id']}: the model file does not hash to the record's perturbedHash")
        by_pert.setdefault(rec["perturbedHash"], set()).add(split)
        by_base.setdefault(rec["baseHash"], {}).setdefault(split, []).append(task)
    both = [h for h, v in by_pert.items() if len(v) == 2]
    if both:
        findings.append(f"{len(both)} curve-balled model(s) occur in both splits")
    shared = [t for v in by_base.values() if len(v) == 2 for t in v["held-out"]]
    if shared:
        findings.append(f"{len(shared)} held-out task(s) are built on a base model that train tasks use as well")
    summary["split"] = {"curveBalledModelsInBothSplits": len(both),
                        "heldOutTasksOnABaseModelOfTheTrainSplit": len(shared),
                        "ofThoseByTaskType": dict(collections.Counter(t["taskType"] for t in shared)),
                        "baseModels": len(by_base)}
    # ---------------------------------------------------------------- references of the distinct models
    stable, unstable = collections.OrderedDict(), collections.OrderedDict()
    for split, p, task, ans in tasks:
        mp = str(Path(p).with_name(task["modelFile"]))
        h = ans["changeRecord"]["perturbedHash"]
        if ans["verification"]["kind"] == "diagnosis":
            if task["taskType"] == "diagnose":
                unstable.setdefault(h, (mp, ans["reference"]["mechanism"]["nMechanisms"]))
        else:
            stable.setdefault(h, (mp, []))[1].append((task, ans))
    items = list(stable.values())[:limit]
    worst_tab = worst_ref = worst_second = 0.0
    n_comb = n_q = n_second = 0
    for k, (mp, tl) in enumerate(items):
        r = audit_stable(mp, tl)
        findings += [f"{Path(mp).name}: {f}" for f in r["findings"]]
        worst_tab, worst_ref = max(worst_tab, r["tableauWorst"]), max(worst_ref, r["referenceWorst"])
        worst_second = max(worst_second, r["tableauWorstSecondOrder"])
        n_comb, n_q, n_second = n_comb + r["combinations"], n_q + r["quantities"], n_second + r["secondOrder"]
        if (k + 1) % 20 == 0:
            print(f"  {k + 1} / {len(items)} stable models", flush=True)
    summary["stableModels"] = {"models": len(items), "combinationsComparedWithTheTableau": n_comb,
                               "ofThoseSecondOrder": n_second, "largestDifferenceEngineTableau": worst_tab,
                               "largestDifferenceEngineTableauSecondOrder": worst_second,
                               "predictQuantitiesRecomputed": n_q, "largestDifferenceReferenceFreshSolve": worst_ref}
    mech = list(unstable.values())[:limit]
    agree = 0
    for mp, n_ref in mech:
        r = audit_mechanism(mp, n_ref)
        findings += [f"{Path(mp).name}: {f}" for f in r["findings"]]
        agree += not r["findings"]
    summary["mechanisms"] = {"models": len(mech), "countConfirmedByTheTableau": agree}
    # ---------------------------------------------------------------- scoring (predict, diagnose)
    n_scored = 0
    for split, p, task, ans in tasks[: None if limit is None else 40 * limit]:
        if task["taskType"] not in ("predict", "diagnose"):
            continue
        ref = reference_answer(task, ans)
        if score_answer(task, ans, ref)["score"] != 1.0:
            findings.append(f"{task['id']}: the reference answer does not score 1.0")
        wrong = copy.deepcopy(ref)
        if task["taskType"] == "predict":
            k0 = sorted(wrong["values"])[0]
            wrong["values"][k0] = wrong["values"][k0] * 1.5 if wrong["values"][k0] else 1.0     # 50 % off
        elif wrong.get("mechanism") is not None:
            wrong["mechanism"] = None
        else:
            wrong["category"] = next(c for c in task["answerFormat"]["category"] if c != ref["category"])
        if not score_answer(task, ans, wrong)["score"] < 1.0:
            findings.append(f"{task['id']}: a wrong answer scores 1.0")
        n_scored += 1
    summary["scoring"] = {"predictAndDiagnoseTasks": n_scored}
    # ---------------------------------------------------------------- leaks (a spread sample)
    from benchfea.bench.env import BenchFEAEnv
    sample = tasks[:: max(1, len(tasks) // 25)][:25]
    envs = {s: BenchFEAEnv(split=s, dataset_dir=str(root)) for s in ("train", "held-out")}
    n_leak = 0
    for split, p, task, ans in sample:
        key, rec = ans["key"], ans["changeRecord"]
        opts = {"family": key["family"], "curve_ball": key["curveBall"], "instance_seed": key["seed"],
                "difficulty": key["difficulty"], "task_type": task["taskType"], "mode": task["mode"]}
        obs, info = envs[split].reset(seed=0, options=opts)
        blob = json.dumps(obs) + json.dumps(info)
        for needle in (key["curveBall"], rec["title"], rec["target"], rec["baseHash"], rec["perturbedHash"],
                       "changeRecord", "whyItMatters", "generatorSeed"):
            if needle and needle in blob:
                findings.append(f"{task['id']}: the first observation / info contains {needle[:40]!r}")
                n_leak += 1
        other = "train" if split == "held-out" else "held-out"
        try:
            envs[other].reset(seed=0, options=opts)
            findings.append(f"{task['id']}: served by the environment of the {other} split")
        except Exception:
            pass
    summary["leaks"] = {"tasksProbed": len(sample), "found": n_leak}
    summary["findings"] = findings
    print(json.dumps(summary, indent=1))
    if "--json" in argv:
        Path(argv[argv.index("--json") + 1]).write_text(json.dumps(summary, indent=1), encoding="utf-8")
    if limit is None:
        from check_dataset_scores import update_checks
        update_checks(root, "audit", dict(summary))
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
