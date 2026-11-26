#!/usr/bin/env python3
"""Rewrite everything that is a NUMBER OF THE DATASET in the documentation, from the dataset on disk.

    python3 tools/build_dataset_docs.py            # rewrite
    python3 tools/build_dataset_docs.py --check    # exit 1 if a file would change

Sources: ``data/tasks/manifest.json`` (counts, cells, quarantine log, versions), the task files it
lists (only to count the tasks that use design results), ``data/CHECKS.json`` (the summaries stored by
``tools/check_dataset_scores.py`` and ``tools/audit_dataset.py``; quoted only when they carry the
SHA-256 of the manifest on disk) and ``benchfea.verify.POLICY`` (the tolerances).

Rewritten:

* in ``data/DATA_CARD.md`` and ``docs/BENCHMARK.md``: every block between the markers
  ``<!-- generated:NAME -->`` and ``<!-- /generated:NAME -->`` (the blocks are defined in ``blocks``);
  text outside the markers is prose and is left alone;
* ``servedCombinations`` of ``src/benchfea/bench/splits.json``: the (family, curve ball) pairs with a
  verified cell in the manifest.

Nothing here says that anybody reviewed, approved or certified a task: the blocks report counts and the
outcomes of computations.
"""
from __future__ import annotations

import collections
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
DOCS = (DATA / "DATA_CARD.md", ROOT / "docs" / "BENCHMARK.md")
SPLITS_FILE = ROOT / "src" / "benchfea" / "bench" / "splits.json"
MARK = re.compile(r"<!-- generated:([a-z-]+) -->\n.*?<!-- /generated:\1 -->", re.S)
TASK_TYPES = ("predict", "diagnose", "repair", "design")
FEW = 3                      # a curve ball with at most this many tasks is listed as a gap

# what a (stage, category) of the quarantine log means; the texts follow benchfea.bench.tasks.build_case
# and benchfea.verify.pipeline
CATEGORY = {
    "curve-ball:not-applicable": "the curve ball does not fit the model (precondition not met, or no measured "
                                 "material effect); not a quarantine",
    "verify:disagreement": "two computations of the same quantity differ by more than the tolerance of the tier",
    "verify:not-covered": "tier F: a design limit state used by the reference lacks coverage",
    "verify:not-representable": "fewer than two foreign solvers can represent the model or give an eligible "
                                "error bound, so it is not cross-checked",
    "verify:not-verifiable": "this kind of reference cannot be verified by the pipeline",
    "verify:not-run": "a tier was not run within the time budget",
    "verify:error": "a check could not be carried out",
    "verify:refused": "the curve-balled model is refused by validation (no analysis to verify)",
    "generate:generation-failed": "the generator found no model of the requested level",
    "generate:generator-error": "the generator raised an error",
    "build:error": "the cell raised an error",
}


def table(head, rows) -> str:
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def ranked(counter) -> str:
    return ", ".join(f"{k} {n}" for k, n in sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))) or "-"


def seeds_text(seeds) -> str:
    s = [str(x) for x in seeds]
    return s[0] if len(s) == 1 else ", ".join(s[:-1]) + " and " + s[-1]


def failing_tier(reason: str) -> str:
    m = re.search(r"tier ([A-FNK]):", reason or "")
    return f"tier {m.group(1)}" if m else "no tier named"


def load():
    manifest_path = DATA / "tasks" / "manifest.json"
    raw = manifest_path.read_bytes()
    m = json.loads(raw)
    checks_path = DATA / "CHECKS.json"
    checks = json.loads(checks_path.read_text(encoding="utf-8")) if checks_path.exists() else {}
    sha = hashlib.sha256(raw).hexdigest()
    checks = {k: v for k, v in checks.items() if v.get("manifestSha256") == sha}
    return m, checks


def design_use(m) -> dict:
    """Tasks that use design results, counted from the task files."""
    out = collections.Counter()
    for t in m["tasks"]:
        if t["taskType"] == "design":
            out["design"] += 1
            continue
        if t["taskType"] not in ("predict", "repair"):
            continue
        task = json.loads((DATA / "tasks" / t["file"]).read_text(encoding="utf-8"))
        if any(q.get("kind") == "designRatio" for q in task.get("quantities", [])):
            out[t["taskType"]] += 1
    return out


def blocks(m, checks) -> dict:
    from benchfea.verify import POLICY

    c = m["counts"]
    cells = m["cells"]
    by_id = {tid: cell for cell in cells for tid in cell["tasks"]}
    seeds = sorted({cell["seed"] for cell in cells})
    families = list(m["cellsByFamily"])
    curve_balls = list(m["cellsByCurveBall"])
    attempts = [e for e in m["quarantineLog"] if e["what"] == "attempt"]
    not_served = [e for e in m["quarantineLog"] if e["what"] == "task type"]
    splits = json.loads(SPLITS_FILE.read_text(encoding="utf-8"))
    held = [tuple(p) for p in splits["heldOutCombinations"]]
    pairs = m["cellsByFamilyAndCurveBall"]
    B = {}

    # ------------------------------------------------------------------ headline, build
    B["headline"] = (f"**{c['tasks']} verified tasks** from {c['cellsVerified']} verified cells of {c['cells']} "
                     f"(every family × every curve ball × seeds {seeds_text(seeds)}): "
                     + ", ".join(f"{c['byTaskType'].get(t, 0)} {t}" for t in TASK_TYPES) + "; "
                     + ", ".join(f"{c['bySplit'].get(s, 0)} {s}" for s in ("train", "held-out")) + ".")
    command = (f"benchfea bench make-dataset --out data --seeds {','.join(str(s) for s in seeds)} --workers 4 "
               f"--cache <dir> --clean")
    v = m["versions"]
    B["build"] = "\n".join([
        "```", command, "python3 tools/check_dataset_scores.py data --workers 4 --cache <dir>",
        "python3 tools/audit_dataset.py data", "python3 tools/build_dataset_docs.py", "```", "",
        "Versions of this build: " + ", ".join(f"{k} {x}" for k, x in v.items()) + ".", "",
        f"Cells planned: {len(families)} families × {len(curve_balls)} curve balls × {len(seeds)} seeds "
        f"({seeds_text(seeds)}) = {c['cells']}; with a verified instance: {c['cellsVerified']}. "
        f"Entries in the quarantine log of the manifest: {len(m['quarantineLog'])}."])

    # ------------------------------------------------------------------ composition
    def count(key):
        out = collections.Counter()
        for t in m["tasks"]:
            out[key(t)] += 1
        return out

    L = [table(["split", "tasks"] + list(TASK_TYPES),
               [[s, c["bySplit"].get(s, 0)] + [c["bySplitAndTaskType"].get(f"{s}|{t}", 0) for t in TASK_TYPES]
                for s in ("train", "held-out")]
               + [["all", c["tasks"]] + [c["byTaskType"].get(t, 0) for t in TASK_TYPES]]), ""]
    mode = count(lambda t: (t["taskType"], t["mode"]))
    L += [table(["mode", "tasks"] + list(TASK_TYPES),
                [[f"`{k}`", n] + [mode.get((t, k), 0) for t in TASK_TYPES] for k, n in c["byMode"].items()]), ""]
    fam = count(lambda t: (t["family"], t["split"]))
    L += [table(["family", "tasks", "train", "held-out", "cells", "cells verified"],
                [[f"`{f}`", o["tasks"], fam.get((f, "train"), 0), fam.get((f, "held-out"), 0), o["cells"], o["verified"]]
                 for f, o in m["cellsByFamily"].items()]), ""]
    lvl = count(lambda t: (t["difficulty"], t["split"]))
    L += [table(["difficulty level", "tasks", "train", "held-out"],
                [[k, n, lvl.get((int(k), "train"), 0), lvl.get((int(k), "held-out"), 0)]
                 for k, n in c["byDifficulty"].items()]), ""]
    eff = count(lambda t: (t["effectClass"], t["taskType"]))
    L += [table(["effect class", "tasks"] + list(TASK_TYPES),
                [[f"`{k}`", n] + [eff.get((k, t), 0) for t in TASK_TYPES] for k, n in c["byEffectClass"].items()]), ""]
    order = count(lambda t: (t["analysisOrder"], t["split"]))
    order_family = count(lambda t: (t["analysisOrder"], t["family"]))
    L += [table(["analysis order of the task's model", "tasks", "train", "held-out", "by family"],
                [[k, n, order.get((k, "train"), 0), order.get((k, "held-out"), 0),
                  ", ".join(f"`{f}` {order_family[(k, f)]}" for f in families if order_family.get((k, f))) or "-"]
                 for k, n in c["byAnalysisOrder"].items()]),
          "", "A model counts as second-order when at least one of its combinations is analysed second order."]
    B["composition"] = "\n".join(L)

    cb_split = collections.Counter()
    for t in m["tasks"]:
        cb_split[(by_id[t["id"]]["curveBall"], t["split"])] += 1
    B["curve-balls"] = table(
        ["curve ball", "cells", "verified", "tasks", "train", "held-out", "cells not verified: stage:category of the last attempt"],
        [[f"`{k}`", o["cells"], o["verified"], o["tasks"], cb_split.get((k, "train"), 0), cb_split.get((k, "held-out"), 0),
          ranked(o["notVerified"])] for k, o in m["cellsByCurveBall"].items()])
    B["families"] = table(
        ["family", "cells", "verified", "tasks", "cells not verified: stage:category of the last attempt"],
        [[f"`{k}`", o["cells"], o["verified"], o["tasks"], ranked(o["notVerified"])]
         for k, o in m["cellsByFamily"].items()])

    # ------------------------------------------------------------------ verified
    tiers = count(lambda t: "".join(t["tiersPassed"]))
    P = POLICY
    L = [table(["tiers passed by the reference", "tasks"], [[f"`{k}`", n] for k, n in sorted(tiers.items())]), "",
         "Tolerances (`benchfea.verify.POLICY`; relative to the scale stated in docs/VERIFICATION.md):", "",
         f"* tier A (ODE oracle, straight member lines): {P['A']['relative']:g};",
         f"* tier C (engine against OpenSeesPy and PyNite): {P['C']['relative']:g} plus the error bounds of the two "
         f"results compared; a foreign result with a bound above {P['C']['maxApproximation']:g} is not eligible; "
         f"{P['C']['requiredSolvers']} eligible foreign solvers are required;",
         f"* tier D (invariants): superposition and reciprocity {P['D']['superposition']:g}, station values from "
         f"statics {P['D']['stationStatics']:g} ({P['D']['stationStaticsIterated']:g} in an iterated analysis), "
         f"station values from an integration of the member equations {P['D']['stationIntegration']:g}; equilibrium "
         f"and energy at the tolerances the solver stores with each invariant;",
         f"* tier E (equivalent models): member split {P['E']['splitFirstOrder']:g} "
         f"({P['E']['splitNonlinear']:g} in an iterated analysis), limit of a stiffness series {P['E']['limit']:g}, "
         f"design variables bitwise;",
         "* tier F (design results): every clause item of every design check used has a hand-calculation test and a "
         "reproduced strength of the AISC Design Examples v16.0 in the test suite (a coverage rule, no tolerance of "
         "its own);",
         f"* mechanism diagnoses (tiers N, K, C, E of the diagnosis route): null-space residual "
         f"{P['diagnosis']['nullResidual']:g}, singular threshold {P['diagnosis']['singular']:g}."]
    B["verified"] = "\n".join(L)

    # ------------------------------------------------------------------ gaps
    du = design_use(m)
    zero = [(k, o) for k, o in m["cellsByCurveBall"].items() if o["tasks"] == 0]
    few = [(k, o) for k, o in m["cellsByCurveBall"].items() if 0 < o["tasks"] <= FEW]

    def why(name):
        reasons = collections.Counter()
        for e in attempts:
            if e["curveBall"] != name or e["outcome"] != "cell not served":
                continue
            key = f"{e['stage']}:{e['category']}"
            if e["stage"] == "verify" and e["category"] not in ("refused",):
                key += f" ({failing_tier(e['reason'])})"
            reasons[key] += 1
        return ranked(reasons)

    L = [f"**Curve balls without any task ({len(zero)} of {len(curve_balls)})** — attempts of their cells, by stage, "
         "category and the tier named in the log:", ""]
    L += [table(["curve ball", "cells", "why (attempts of cells not served)"],
                [[f"`{k}`", o["cells"], why(k)] for k, o in zero]) if zero else "None.", ""]
    L += [f"**Curve balls with 1 to {FEW} tasks ({len(few)})**:", ""]
    L += [table(["curve ball", "tasks", "cells verified", "why (attempts of cells not served)"],
                [[f"`{k}`", o["tasks"], f"{o['verified']} of {o['cells']}", why(k)] for k, o in few]) if few else "None.", ""]
    design_cb = [k for k in curve_balls if k.startswith("design_")]
    design_tasks = sum(m["cellsByCurveBall"][k]["tasks"] for k in design_cb)
    L += [f"**Design results.** Design tasks: {du.get('design', 0)}. Predict tasks that ask for a design ratio: "
          f"{du.get('predict', 0)}. Repair tasks with a design-ratio quantity: {du.get('repair', 0)}. Tasks of the "
          f"{len(design_cb)} design-only curve balls (`design_*`): {design_tasks}. Task types not served on verified "
          f"cells because design results did not verify are counted in the quarantine summary.", ""]
    so = c["byAnalysisOrder"].get("second-order", 0)
    L += [f"**Second order.** {so} of {c['tasks']} tasks are on a model with a second-order combination "
          f"({100.0 * so / max(c['tasks'], 1):.1f} %).", ""]
    served_share = []
    for f, o in m["cellsByFamily"].items():
        applicable = o["cells"] - sum(n for k, n in o["notVerified"].items() if k == "curve-ball:not-applicable")
        served_share.append([f"`{f}`", applicable, o["verified"],
                             f"{100.0 * o['verified'] / applicable:.0f} %" if applicable else "-"])
    L += ["**Share of the applicable cells that are served, by family** (applicable = the last attempt was not "
          "\"curve ball not applicable\"):", "",
          table(["family", "applicable cells", "verified", "share"], served_share), ""]
    held_rows = []
    for f, cb in held:
        o = pairs.get(f"{f}|{cb}", {"cells": 0, "verified": 0, "tasks": 0, "reasons": {}})
        held_rows.append([f"`{f}`", f"`{cb}`", o["verified"], o["tasks"], ranked(o["reasons"])])
    unserved = sum(1 for r in held_rows if r[2] == 0)
    L += [f"**Held-out combinations** ({len(held)}; {len(held) - unserved} with a verified instance, {unserved} "
          "without — those are not served and were not replaced):", "",
          table(["family", "curve ball", "cells verified", "tasks", "cells not verified"], held_rows)]
    B["gaps"] = "\n".join(L)

    # ------------------------------------------------------------------ quarantine
    att = collections.Counter(f"{e['stage']}:{e['category']}" for e in attempts)
    att_tier = collections.Counter()
    for e in attempts:
        if e["stage"] == "verify" and e["category"] != "refused":
            att_tier[(f"{e['stage']}:{e['category']}", failing_tier(e["reason"]))] += 1
    replaced = collections.Counter(f"{e['stage']}:{e['category']}" for e in attempts
                                   if e["outcome"] != "cell not served")
    rows = []
    for k, n in sorted(att.items(), key=lambda kv: (-kv[1], kv[0])):
        tiers_ = ", ".join(f"{t} {x}" for (kk, t), x in sorted(att_tier.items(), key=lambda kv: -kv[1]) if kk == k)
        rows.append([f"`{k}`", n, replaced.get(k, 0), tiers_ or "-", CATEGORY.get(k, "")])
    L = [f"{len(m['quarantineLog'])} entries: {len(attempts)} attempts that were skipped and {len(not_served)} task "
         "types not served on a verified cell. Every entry carries its reason text.", "",
         table(["attempts skipped: stage:category", "attempts", "of those replaced by a later attempt",
                "tier named in the reason", "meaning"], rows), "",
         table(["task type not served on a verified cell: reason", "cells"],
               [[k, n] for k, n in m["taskTypesNotServed"].items()]) if m["taskTypesNotServed"] else ""]
    B["quarantine"] = "\n".join(L).rstrip()

    # ------------------------------------------------------------------ split
    seed_rule = splits["seedRule"]
    B["split"] = "\n".join([
        table(["split", "tasks", "cells verified", "seeds"],
              [[s, c["bySplit"].get(s, 0), sum(1 for x in cells if x["split"] == s and x["status"] == "verified"),
                seeds_text(sorted({x["seed"] for x in cells if x["split"] == s})) or "-"]
               for s in ("train", "held-out")]), "",
        f"Held-out: {len(held)} (family, curve ball) combinations with every seed, and of every other combination "
        f"the seeds with `seed % {seed_rule['modulus']}` in {seed_rule['heldOutResidues']}."])

    # ------------------------------------------------------------------ checks
    L = []
    s = checks.get("scores")
    if s:
        r = s["repair"]
        L += [f"* `tools/check_dataset_scores.py` (normal scoring path, repairs and designs re-solved through the "
              f"pipeline): the answer built from the hidden reference scores 1.0 on {s['referenceAnswerScores1']} of "
              f"{s['tasks']} tasks; the model as given, handed in as the repair, scores below the reference answer "
              f"on {r['givenModelScoresBelowReference']} of {r['tasks']} repair tasks (0 on "
              f"{r['givenModelScoresZero']}; largest score {r['largestScoreOfAGivenModel']}). "
              f"Findings: {len(s['findings'])}."]
    else:
        L += ["* `tools/check_dataset_scores.py`: **not run on this build**."]
    a = checks.get("audit")
    if a:
        sm, sp, mech = a["stableModels"], a["split"], a["mechanisms"]
        first = sm["combinationsComparedWithTheTableau"] - sm["ofThoseSecondOrder"]
        L += [f"* `tools/audit_dataset.py` (not the pipeline that produced the references): {a['tasks']} tasks on "
              f"{sp['baseModels']} base models; curve-balled models in both splits: "
              f"{sp['curveBalledModelsInBothSplits']}; held-out tasks on a base model that train tasks also use: "
              f"{sp['heldOutTasksOnABaseModelOfTheTrainSplit']}. {sm['models']} stable models solved again by the "
              f"engine ({sm['predictQuantitiesRecomputed']} reference values recomputed, largest relative "
              f"difference {sm['largestDifferenceReferenceFreshSolve']:.1e}) and by the independent tableau solver "
              f"of `tests/independent`: {first} first-order combinations, largest relative difference "
              f"{sm['largestDifferenceEngineTableau']:.1e}"
              + (f"; {sm['ofThoseSecondOrder']} second-order combinations with the tableau refined (it takes the "
                 f"axial force constant over a piece of a member; "
                 + ", ".join(str(n) for n in sm["secondOrderPiecesPerMember"]) + " pieces per member, raised while the "
                 f"difference to the engine is above 1e-3): {sm['secondOrderCombinationsRefinedBeyondTheFirstPair']} "
                 f"needed more than 16 pieces, the largest difference at the last count run is "
                 f"{sm['largestDifferenceEngineTableauSecondOrder']:.1e}, and between the last two counts run the "
                 f"difference falls by a factor of at least {sm['smallestReductionOfTheSecondOrderDifference']:.1f} "
                 f"(a finding below 4)"
                 if sm["ofThoseSecondOrder"] else "")
              + f". Number of mechanisms of {mech['models']} unstable models from the null space of the tableau: "
              f"{mech['countConfirmedByTheTableau']} agree. Right and wrong answers through the scorer: "
              f"{a['scoring']['predictAndDiagnoseTasks']} predict and diagnose tasks. First observation of "
              f"{a['leaks']['tasksProbed']} tasks probed for parts of the change record: {a['leaks']['found']} "
              f"found. Findings: {len(a['findings'])}."]
    else:
        L += ["* `tools/audit_dataset.py`: **not run on this build**."]
    L += ["", "These are computations by this project's code. Nobody has reviewed, approved or certified a task, "
          "a model or a reference."]
    B["checks"] = "\n".join(L)
    return B


def served_combinations(m) -> list:
    return [k.split("|") for k, v in m["cellsByFamilyAndCurveBall"].items() if v["verified"]]


def outputs() -> dict:
    m, checks = load()
    B = blocks(m, checks)
    out = {}
    for path in DOCS:
        text = path.read_text(encoding="utf-8")
        used = set()

        def sub(match):
            name = match.group(1)
            if name not in B:
                raise SystemExit(f"{path.name}: unknown generated block {name!r}")
            used.add(name)
            return f"<!-- generated:{name} -->\n{B[name]}\n<!-- /generated:{name} -->"

        new = MARK.sub(sub, text)
        if not used:
            raise SystemExit(f"{path.name}: no generated block found")
        if new.count("<!-- generated:") != len(MARK.findall(new)):
            raise SystemExit(f"{path.name}: a generated block is not closed")
        out[path] = new
    d = json.loads(SPLITS_FILE.read_text(encoding="utf-8"))
    d["servedCombinations"] = served_combinations(m)
    out[SPLITS_FILE] = json.dumps(d, indent=1) + "\n"
    return out


def main(argv) -> int:
    check = "--check" in argv
    changed = []
    for path, text in outputs().items():
        if path.read_text(encoding="utf-8") != text:
            changed.append(path.relative_to(ROOT).as_posix())
            if not check:
                path.write_text(text, encoding="utf-8")
    if check and changed:
        print("out of date (python3 tools/build_dataset_docs.py):", ", ".join(changed))
        return 1
    print("up to date" if not changed else "rewritten: " + ", ".join(changed))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
