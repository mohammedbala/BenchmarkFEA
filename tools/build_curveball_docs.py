#!/usr/bin/env python3
"""Build docs/CURVEBALLS.md from the curve-ball registry and from ACTUAL runs.

    python3 tools/build_curveball_docs.py                 # writes docs/CURVEBALLS.md  (several minutes)
    python3 tools/build_curveball_docs.py --cache DIR     # keep the generated base models and the runs in DIR
    python3 tools/build_curveball_docs.py --check         # exit 1 if the file on disk differs

Everything numeric in the document is produced here by running the code: the base models are
generated (``benchfea.bench.generate`` with the arguments in ``BANK``), every curve ball is applied
to every base model with seed 0, and the applicability table and the before / after numbers of the
examples are taken from those runs.  The prose of each curve ball (what it changes, why it matters,
what it moves, precondition) is the text of its registration in ``benchfea.bench.curveballs``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from benchfea import bench  # noqa: E402
from benchfea.bench import curveballs as cb  # noqa: E402

#: base models: (family, seed, difficulty, forced features).  At least one per family; a second one
#: where a family's optional features (offsets, second order, springs) open further curve balls.
BANK = [
    ("continuous_beam", 1, 1, {}),
    ("continuous_beam", 1, 2, {}),
    ("continuous_beam", 1, 2, {"overhang": True}),
    ("portal_frame", 1, 1, {}),
    ("portal_frame", 2, 2, {"rigidOffsets": True}),
    ("portal_frame", 1, 3, {}),
    ("multistorey_2d", 1, 2, {}),
    ("multistorey_2d", 1, 3, {}),
    ("multistorey_3d", 1, 3, {}),
    ("braced_frame", 1, 2, {}),
    ("pipe_rack", 1, 3, {}),
    ("stair_frame", 1, 1, {}),
    ("stair_frame", 1, 2, {}),
    ("diaphragm_frame", 1, 3, {}),
]
SEED = 0
ABBREVIATION = {"changes-analysis": "A", "design-only": "D", "creates-mechanism": "M", "invalid-model": "I"}


def model_key(family, seed, level, features) -> str:
    flags = "".join(f"+{k}" if v else f"-{k}" for k, v in sorted(features.items()))
    return f"{family}/seed{seed}/L{level}{flags}"


def model_call(family, seed, level, features) -> str:
    extra = f", features={features!r}" if features else ""
    return f'generate("{family}", seed={seed}, difficulty={level}{extra})'


def base_model(entry, cache: Path | None) -> dict:
    name = model_key(*entry).replace("/", "_") + ".json"
    if cache is not None and (cache / name).is_file():
        return json.loads((cache / name).read_text())
    g = bench.generate(entry[0], entry[1], entry[2], features=entry[3])
    if cache is not None:
        cache.mkdir(parents=True, exist_ok=True)
        (cache / name).write_text(g.json())
    return json.loads(g.json())


def run_all(cache: Path | None, progress=True) -> dict:
    """{model key: {curve ball: {"record": dict} | {"reason": text}}} (cached per model in ``cache``)."""
    out = {}
    for entry in BANK:
        key = model_key(*entry)
        run_file = None if cache is None else cache / (key.replace("/", "_") + f".runs-v{cb.CURVEBALL_VERSION}.json")
        if run_file is not None and run_file.is_file():
            out[key] = json.loads(run_file.read_text())
            continue
        doc = base_model(entry, cache)
        rows = {}
        for name in cb.available():
            try:
                p = cb.apply(doc, name, SEED)
            except cb.NotApplicable as e:
                rows[name] = {"reason": e.reason}
            else:
                # the record must describe the difference exactly -- checked on every run of this tool
                assert json.dumps(p.record.apply(doc)) == json.dumps(p.document), (key, name)
                assert json.dumps(p.record.revert(p.document)) == json.dumps(doc), (key, name)
                rows[name] = {"record": p.record.to_dict()}
            if progress:
                state = rows[name].get("record", {}).get("effectClass", "not applicable")
                print(f"  {key:44s} {name:32s} {state}", file=sys.stderr, flush=True)
        out[key] = rows
        if run_file is not None:
            run_file.write_text(json.dumps(rows))
    return out


def examples_of(runs: dict, name: str) -> list:
    """[(BANK entry, record)] shown for a curve ball: the first model of the bank for each effect class."""
    shown, out = set(), []
    for entry in BANK:
        row = runs[model_key(*entry)][name]
        if "record" in row and row["record"]["effectClass"] not in shown:
            shown.add(row["record"]["effectClass"])
            out.append((entry, row["record"]))
    return out


def verify_examples(runs: dict, cache: Path | None, budget: float, progress=True) -> dict:
    """Run the verification pipeline (``benchfea.verify.verify_model``) on the perturbed model of every
    example: {(model key, curve ball): {"analysis": {...}, "design": {...} | None}}.  ``analysis`` is
    the record without the design results; ``design`` (design-only curve balls) with tier F."""
    from benchfea import verify as V

    def brief(rec) -> dict:
        q = rec.quarantine or {}
        return {"status": rec.status, "kind": rec.kind, "tiersPassed": list(rec.tiers_passed),
                "tiersNotPassed": [k for k, t in rec.tiers.items() if t.applies and t.passed is not True],
                "category": q.get("category"), "diagnosis": q.get("diagnosis")}

    file = None if cache is None else cache / f"verification-v{cb.CURVEBALL_VERSION}.json"
    stored = json.loads(file.read_text()) if file is not None and file.is_file() else {}
    out = {}
    for name in cb.available():
        for entry, record in examples_of(runs, name):
            key = model_key(*entry)
            tag = f"{key}|{name}|{record['perturbedHash']}"
            if tag not in stored or (record["effectClass"] == "design-only" and "baseDesign" not in stored[tag]):
                doc = cb.ChangeRecord.from_dict(record).apply(base_model(entry, cache))
                row = {"analysis": brief(V.verify_model(doc, budget=budget)), "design": None}
                if record["effectClass"] == "design-only":
                    row["design"] = brief(V.verify_model(doc, design=True, budget=budget))
                    # the base model with its design results, for comparison (tier F is about the checks used)
                    row["baseDesign"] = brief(V.verify_model(base_model(entry, cache), design=True, budget=budget))
                stored[tag] = row
                if file is not None:
                    file.write_text(json.dumps(stored))
                if progress:
                    print(f"  verify {key:40s} {name:32s} {row['analysis']['status']}", file=sys.stderr, flush=True)
            out[(key, name)] = stored[tag]
    return out


def verification_line(v: dict | None) -> list:
    if v is None:
        return []
    a = v["analysis"]
    tiers = ", ".join(a["tiersPassed"]) or "none"
    if a["status"] == "verified":
        what = "mechanism diagnosis" if a["kind"] == "diagnosis" else "analysis results"
        text = (f"* Verification run on this perturbed model (`benchfea verify`, {what}): `verified` -- passed tiers "
                f"{tiers}, which are all the tiers that apply (meaning and tolerances: docs/VERIFICATION.md).")
    elif a["status"] == "refused":
        text = ("* Verification run on this perturbed model (`benchfea verify`): `refused` -- the engine produces no "
                "reference solution for it, as expected for an invalid model.")
    else:
        reason = (a["diagnosis"] or "")
        reason = reason if len(reason) < 300 else reason[:296] + " ..."
        what = "mechanism diagnosis" if a["kind"] == "diagnosis" else "analysis results"
        text = (f"* Verification run on this perturbed model (`benchfea verify`, {what}): `quarantined` "
                f"(category `{a['category']}`; tiers passed: {tiers}; not passed: {', '.join(a['tiersNotPassed'])}) -- "
                f"{reason} **This reference must not be published.**")
    lines = [text]
    d = v.get("design")
    if d is not None:
        if d["status"] == "verified":
            lines.append(f"* The same with the design results (`--design`, tier F): `verified` -- passed tiers "
                         f"{', '.join(d['tiersPassed'])}.")
        else:
            reason = (d["diagnosis"] or "")
            reason = reason if len(reason) < 300 else reason[:296] + " ..."
            lines.append(f"* The same with the design results (`--design`, tier F): `{d['status']}` (category "
                         f"`{d['category']}`; not passed: {', '.join(d['tiersNotPassed'])}) -- {reason}")
    return lines


# ------------------------------------------------------------------------------------- formatting
def num(x, digits=4) -> str:
    if x is None:
        return "n/a"
    return f"{float(x):.{digits}g}"


def pct(x) -> str:
    return "n/a" if x is None else (f"{100 * x:.1f} %" if x < 10 else f"{x:.3g} x the base value")


def short(value, limit=70) -> str:
    text = json.dumps(value)
    return text if len(text) <= limit else text[:limit - 4] + " ..."


def changes_lines(record, limit=6) -> list:
    lines = []
    for c in record["changes"][:limit]:
        if c["op"] == "replace":
            lines.append(f"  * `{c['path']}`: `{short(c['before'])}` -> `{short(c['after'])}`")
        elif c["op"] == "add":
            lines.append(f"  * `{c['path']}` added: `{short(c['after'])}`")
        else:
            lines.append(f"  * `{c['path']}` removed (was `{short(c['before'])}`)")
    if len(record["changes"]) > limit:
        lines.append(f"  * ... and {len(record['changes']) - limit} more change(s) of the same kind")
    return lines


def example_block(key, entry, record, verification=None) -> list:
    e = record["effect"]
    lines = [f"Example (actual run, seed {SEED}) on `{model_call(*entry)}`; target: {record['target']}.", "",
             "* Changes recorded:"]
    lines += changes_lines(record)
    cls = record["effectClass"]
    if cls == "creates-mechanism":
        d = record["reference"]["diagnosis"]
        lines.append(f"* Effect class **creates-mechanism**: the base model is analysed, the perturbed one is refused "
                     f"as a mechanism ({d['nMechanisms']} mechanism(s)). Reference = the diagnosis; first mode: "
                     f"*{d['modes'][0]['description']}* (kind `{d['modes'][0]['kind']}`, cause `{d['modes'][0]['cause']}`).")
    elif cls == "invalid-model":
        probs = record["reference"]["problems"]
        lines.append("* Effect class **invalid-model**: validation refuses the perturbed model. Reference = the "
                     "problems: " + "; ".join(f"`{p['code']}` at `{p['path']}`" for p in probs[:3]) + ".")
    else:
        at = e["analysisChangeAt"]
        if e["analysisIdentical"]:
            lines.append("* Analysis: every result array of every combination is **bitwise identical** to the base "
                         "model.")
        else:
            lines.append(f"* Analysis change (largest, relative to the base peak of its kind): **{pct(e['analysisChange'])}** "
                         f"-- {at.get('quantity')} at `{at.get('where')}`, combination `{at.get('combination')}`.")
            cid = at.get("combination")
            b, p = e["base"]["combinations"].get(cid), e["perturbed"]["combinations"].get(cid)
            if b and p:
                ub, up = e["base"]["units"], e["perturbed"]["units"]
                lines.append(
                    f"* Combination `{cid}`, before -> after: largest displacement of a member axis "
                    f"{num(b['maxMemberDisplacement'])} {ub['length']} -> {num(p['maxMemberDisplacement'])} "
                    f"{up['length']}; largest reaction force "
                    f"{num(b['maxReactionForce'])} {ub['force']} -> {num(p['maxReactionForce'])} {up['force']}; largest "
                    f"member bending moment {num(b['maxMemberMoment'])} {ub['force']}·{ub['length']} -> "
                    f"{num(p['maxMemberMoment'])} {up['force']}·{up['length']}.")
        if e.get("designRatioChange") is not None:
            w = e.get("designRatioChangeAt")
            db, dp = e["base"].get("design", {}), e["perturbed"].get("design", {})
            text = (f"* Design: status {db.get('status')} -> {dp.get('status')}; largest ratio {num(db.get('maxRatio'), 3)} "
                    f"(`{db.get('maxRatioMember')}`) -> {num(dp.get('maxRatio'), 3)} (`{dp.get('maxRatioMember')}`)")
            if w:
                text += (f"; largest change of a governing ratio: design member `{w['member']}` {num(w['base'], 3)} -> "
                         f"{num(w['perturbed'], 3)} ({pct(e['designRatioChange'])})")
            if e.get("designPassFailChanged"):
                text += "; pass / fail changed for: " + ", ".join(f"`{m}`" for m in e["designPassFailChanged"])
            lines.append(text + ".")
        extra = {k: e[k] for k in ("baseAmplification", "physicalFactor", "relabelled", "factor") if k in e}
        if "baseAmplification" in extra:
            lines.append("* Measured second-order amplification of the base model (largest horizontal translation, "
                         "second order / first order): "
                         + ", ".join(f"`{k}` {v:.3f}" for k, v in extra["baseAmplification"].items()) + ".")
        if "physicalFactor" in extra:
            lines.append(f"* The {extra['relabelled']} unit factor is {num(extra['physicalFactor'], 6)}: the numbers "
                         "of the results in file units are unchanged, their physical meaning is not.")
        lines.append(f"* Materiality criterion met: `{e['criterion']}`.")
    v = record["verification"]
    lines.append(f"* Expected verification route: {v['route']}.")
    lines += verification_line(verification)
    if record["rejected"]:
        lines.append(f"* {len(record['rejected'])} target(s) were tried and dropped before this one, e.g. "
                     f"*{record['rejected'][0]['target']}* -- {record['rejected'][0]['reason']}.")
    return lines


INTRO = """# Curve balls

A **curve ball** is a named, seeded, deterministic change of a model that *matters*. The benchmark
hands an agent a model in which something has been changed; the agent's answer is scored against
the reference solution of the changed model and, for a diagnosis task, against the **change
record**. Code: `benchfea.bench.curveballs`; command line: `benchfea perturb`.

> This file is generated by `tools/build_curveball_docs.py`. Every number in it comes from running
> the code on generated models; the prose of each curve ball is the text of its registration.
> Nothing here says that a result has been reviewed or approved by anyone: "verified" in this
> project means only *passed the tiers listed, at the tolerances stated* (docs/VERIFICATION.md), and
> the verification route given for each example is an **expectation**, not a result.

```python
from benchfea.bench import generate, curveballs

base = generate("portal_frame", seed=1, difficulty=1).document
p = curveballs.apply(base, "support_fixed_to_pinned", seed=0)   # Perturbed(document, record)
p.document                     # the changed input document
p.record.changes               # [Change(op, path, before, after, index)] -- JSON pointers
p.record.entities              # {"nodes": [...], "supports": [...]} ids touched
p.record.effect_class          # "changes-analysis"
p.record.effect                # measured: analysisChange, designRatioChange, where, before / after
p.record.apply(base)  == p.document      # the record IS the difference
p.record.revert(p.document) == base
curveballs.apply(base, "offset_removed", seed=0)                # raises NotApplicable(reason)
```

```
benchfea perturb --list
benchfea perturb model.json release_added --seed 3 -o variant.json --record variant.record.json
benchfea perturb model.json support_removed --show-record
benchfea perturb model.json --applicable           # which curve balls fit this model
```

Exit codes of `benchfea perturb`: 0 applied; 7 not applicable (the reason is printed); 2 usage,
unknown curve ball or unreadable model.

## 1. How a curve ball is applied

1. The base model is analysed (strict: all invariants passed, every iteration converged) and
   designed. A base model that is not analysable is refused.
2. The curve ball lists its **candidate targets** in the base model (member ends, supports, loads,
   design blocks, ...). No candidate: `NotApplicable`, naming the precondition that is not met.
3. The candidates are put in a seeded order: SHA-256 of
   `"benchfea.bench|curveball|<version>|<name>|<seed>"` -> PCG64, Fisher-Yates on the raw integer
   stream (the same construction as the generators; never Python's `hash()`). Candidates that keep
   the model inside what the cross-check is expected to follow (`bench.tier_c_expected`) are tried
   before candidates that take a cross-checkable model outside it.
4. Candidate by candidate (at most `MAX_CANDIDATES` = {max_candidates}) the perturbed model is
   analysed and designed and **the effect is measured** against the base model. The first candidate
   whose effect class is one the curve ball declares and whose effect is material is taken. Every
   candidate dropped on the way is listed in the record (`rejected`, with the measured values).
5. No candidate passes: `NotApplicable`, listing what was tried and what was measured.

A perturbed model that the engine refuses for any reason other than a plain mechanism or a
validation error (beyond the elastic critical load, not converged, an invariant failed, a
mechanism that exists only through slack tension-only members, an ill-conditioned matrix) has no
reference solution; such a candidate is never accepted.

## 2. Effect classes and materiality

| effect class | what happened | material when | reference solution of the perturbed model |
|---|---|---|---|
| `changes-analysis` | the analysis results changed | **analysis change ≥ {a:.0%}**, or a design member changes between pass / conditional and fail; for the curve balls that change the *resistance* of a member (its section or orientation) also: the governing **design ratio of a design member changes by ≥ {d:.0%}** | results (and design) of the perturbed model |
| `design-only` | the analysis is **bitwise identical** (every array of every combination), the design changed | governing design ratio of a design member changes by ≥ {d:.0%}, or a design member changes pass / fail | the same results; the design of the perturbed model |
| `creates-mechanism` | the base model is analysed, the perturbed one is refused with `stability.mechanism` | always | the solver's structured **mechanism diagnosis** (`reference.diagnosis` in the record) |
| `invalid-model` | the base model is valid, the perturbed one is refused by validation | always | the validator's list of problems (`reference.problems`) |

**Analysis change** (`effect.analysisChange`), measured in physical units (N, mm) on the nodes and
members both models have, for every combination: the largest difference, entity by entity, of

* the nodal translations and the peak displacement of each member axis, divided by the largest
  such displacement of the base model in that combination;
* the reaction forces and member forces (end forces as global vectors, peak axial force and peak
  resultant shear along each member), divided by the largest such force of the base model in that
  combination;
* the reaction moments and member moments (end moments as global vectors, peak torque and peak
  resultant bending moment along each member), divided by the largest such moment of the base model
  in that combination;

and the largest of these over the combinations. The difference is taken where it occurs, so a load
path that moves is measured even if the peak values stay the same; it is normalised by the governing
(peak) value of its kind, so a large relative change of a negligible quantity does not count; member
quantities are global vectors or resultants, so turning a member's local axes is not by itself a
change.
`effect.analysisChangeAt` names the quantity, the node or member and the combination.

**Design ratio change** (`effect.designRatioChange`): largest relative change of the governing
ratio of a design member (`effect.designRatioChangeAt` names it and its limit states).

Which criteria count for a curve ball is listed with it ("Material through"). A curve ball that
changes only the demands (loads, supports, releases, offsets, constraints, units) is not made
material by a design ratio alone: with the resistance unchanged a ratio moves with the member
forces, and a ratio that jumps while no governing force moves by {a:.0%} would be a discontinuity of the
design checks, not an effect of the change. `second_order_off` counts the analysis change only: it
is applicable only where switching the second-order analysis off moves the analysis results by at
least {a:.0%}.

The thresholds are `curveballs.THRESHOLDS`; they are the same for every curve ball and model, and a
target that misses them is dropped -- never forced.

## 3. The change record

`record.to_dict()` (format `{fmt}` {ver}; `benchfea perturb --record FILE`):

| field | content |
|---|---|
| `curveBall`, `seed`, `category`, `title`, `curveBallVersion` | which curve ball, applied with which seed |
| `effectClass` | the measured class (one of those the curve ball declares) |
| `target` | one line naming what was chosen |
| `changes` | the exact difference of the two documents: `{{"op": "replace" \\| "add" \\| "remove", "path": JSON pointer, "before", "after", "index"}}`. `index` is the position of an added / removed key in its object, so that key order is reproduced. Arrays whose length changes are replaced as a whole |
| `entities` | ids touched, by kind (`members`, `nodes`, `supports`, `constraints`, `sections`, `loadCases`, `combinations`, `designMembers`, `units`) |
| `whyItMatters` | the physical reason, in plain language |
| `precondition` | what the base model must have |
| `effect` | the measured effect: `class`, `criterion`, `thresholds`, `analysisIdentical`, `analysisChange`, `analysisChangeAt`, `analysisChangeByCombination`, `designRatioChange`, `designRatioChangeAt`, `designPassFailChanged`, `designStatusChanges`, `base` / `perturbed` summaries in file units, and facts of the curve ball (factor applied, amplification, ...) |
| `reference` | `kind`: `results` \\| `diagnosis` \\| `validation-problems`; for a mechanism the full diagnosis, for an invalid model the problems |
| `verification` | expected verification route of the perturbed model |
| `rejected` | candidates measured and dropped before the accepted one |
| `baseHash`, `perturbedHash` | SHA-256 of the two documents |

The change list is not written by the curve balls: it is the structural difference of the base and
the perturbed document (`curveballs.diff`). `record.apply(base)` and `record.revert(perturbed)`
check every `before` value while applying; the tests (and this script, on every run) require both
to reproduce the other document exactly, key order included. A diagnosis answer can therefore be
matched against `changes[*].path` and `entities`.

## 4. Verification route

Curve balls prefer targets that keep the perturbed model verifiable. Where that is not possible the
record says so (`verification.route`, `verification.notes`):

* `changes-analysis` / `design-only`: `benchfea verify` on the perturbed model. If the perturbed
  model uses a feature that fewer than two independent solvers can represent (for example rigid end
  offsets together with second-order combinations), tier C cannot pass and the pipeline quarantines
  the reference; a quarantined reference is not published.
* `creates-mechanism`: the diagnosis tiers of docs/VERIFICATION.md section 4.
* `invalid-model`: there is nothing to analyse; the reference answer is the validator's problem
  list, reproducible with `benchfea validate`.

The route in a record is derived from the features of the perturbed model; it is not the outcome of
a verification run. Two cases in which the expectation is negative, both measured on the examples
below and reported there as quarantined:

* **A loaded node without members.** Removing a diaphragm (or detaching its master) leaves the
  master node, which carries the storey loads, connected to nothing. The engine refuses the model
  as a mechanism (mode kind `node`, cause `no-stiffness`: the node is loaded in DOFs where it has
  no stiffness). The diagnosis tiers do not confirm it: the independently assembled stiffness of
  tier N leaves memberless nodes out, finds no zero-stiffness motion, and so disagrees with the
  number of mechanisms; tiers K and C follow. The reference of such a perturbed model is
  quarantined (category `disagreement`). The record says `expectedVerifiable: false`.
* **A node on the axis of a member.** A reversed end offset puts the flexible end beyond its node,
  so the member passes over its own end node. The PyNite translator refuses such a topology (PyNite
  would connect the node to the member), which leaves one foreign solver: tier C is
  `not-representable`. The record's expectation includes this check.

To accept only targets whose reference actually passes, apply with verification:

```python
curveballs.apply(base, "offset_removed", seed=0, verify=True)   # record.verification["run"]
```

```
benchfea perturb model.json offset_removed --require-verified
```

Every candidate that is material is then put through `benchfea.verify.verify_model` (analysis
results or mechanism diagnosis; for a design-only effect, or with `verify="design"`, also the design
results and tier F) and dropped unless the status is `verified`; if none is left the curve ball is
`NotApplicable`. `verified` means only: passed the tiers listed in `record.verification["run"]`, at
the tolerances of docs/VERIFICATION.md. This costs far more than the measurement.
"""


def build(cache: Path | None, progress=True, verify: bool = True, budget: float = 300.0) -> str:
    runs = run_all(cache, progress)
    verified = verify_examples(runs, cache, budget, progress) if verify else {}
    families = list(dict.fromkeys(e[0] for e in BANK))
    names = cb.available()
    intro = INTRO
    for token, value in (("{max_candidates}", str(cb.MAX_CANDIDATES)),
                         ("{a:.0%}", f"{cb.THRESHOLDS['analysisChange']:.0%}"),
                         ("{d:.0%}", f"{cb.THRESHOLDS['designRatioChange']:.0%}"),
                         ("{fmt}", cb.RECORD_FORMAT), ("{ver}", cb.RECORD_VERSION), ("{{", "{"), ("}}", "}")):
        intro = intro.replace(token, value)
    out = [intro]

    # ---- applicability table
    out.append("## 5. Applicability: family x curve ball\n")
    out.append(f"Measured by applying every curve ball with seed {SEED} to these generated models "
               "(`benchfea.bench.generate`):\n")
    for entry in BANK:
        out.append(f"* `{model_call(*entry)}`")
    out.append("")
    out.append("A cell lists the effect classes obtained on the models of the family: **A** = changes-analysis, "
               "**D** = design-only, **M** = creates-mechanism, **I** = invalid-model; `-` = not applicable to any "
               "of the models above (precondition not met, or no target with a material effect). Another seed or "
               "another model of the family (other feature flags, other level) can give a different cell: the "
               "table states what was measured, not what is possible.\n")
    out.append("| curve ball | " + " | ".join(f"`{f}`" for f in families) + " |")
    out.append("|---|" + "---|" * len(families))
    never = []
    for name in names:
        cells, any_ok = [], False
        for fam in families:
            classes = []
            for entry in BANK:
                if entry[0] != fam:
                    continue
                row = runs[model_key(*entry)][name]
                if "record" in row:
                    a = ABBREVIATION[row["record"]["effectClass"]]
                    if a not in classes:
                        classes.append(a)
            cells.append(" ".join(classes) if classes else "-")
            any_ok = any_ok or bool(classes)
        if not any_ok:
            never.append(name)
        out.append(f"| `{name}` | " + " | ".join(cells) + " |")
    out.append("")
    total = sum(1 for k in runs for n in names if "record" in runs[k][n])
    out.append(f"{total} of {len(runs) * len(names)} (model, curve ball) pairs were applicable. "
               + ("Every curve ball is applicable to at least one family." if not never else
                  "NOT applicable to any model above: " + ", ".join(f"`{n}`" for n in never) + "."))
    out.append("")

    # ---- the curve balls
    if verified:
        counts = {}
        for v in verified.values():
            counts[v["analysis"]["status"]] = counts.get(v["analysis"]["status"], 0) + 1
        design_rows = [v["design"] for v in verified.values() if v["design"] is not None]
        out.append("### Verification runs of the examples\n")
        out.append(f"The perturbed model of each of the {len(verified)} examples of section 6 was put through "
                   "`benchfea.verify.verify_model` (all tiers, analysis results or mechanism diagnosis): "
                   + ", ".join(f"{n} `{k}`" for k, n in sorted(counts.items()))
                   + ". The outcome is stated with each example; `verified` means only that the reference passed "
                     "the tiers listed, at the tolerances of docs/VERIFICATION.md. A `quarantined` example shows a "
                     "perturbed model whose reference is not publishable as it stands.")
        if design_rows:
            dc, bc = {}, {}
            for d in design_rows:
                dc[d["status"]] = dc.get(d["status"], 0) + 1
            base_rows = [v["baseDesign"] for v in verified.values() if v.get("baseDesign") is not None]
            for b in base_rows:
                bc[b["status"]] = bc.get(b["status"], 0) + 1
            out.append("")
            out.append(f"For the {len(design_rows)} design-only examples the run was repeated with the design results "
                       "(tier F): " + ", ".join(f"{n} `{k}`" for k, n in sorted(dc.items()))
                       + ". Tier F is about the design checks a model uses, not about the curve ball: the BASE models "
                         "of these examples, verified with their design results in the same way, gave "
                       + ", ".join(f"{n} `{k}`" for k, n in sorted(bc.items()))
                       + " (categories: " + ", ".join(sorted({str(b['category']) for b in base_rows})) + "). As long "
                         "as a limit state used by a model lacks full tier-F coverage, neither the base nor the "
                         "perturbed design results of that model are publishable as verified; the analysis results "
                         "of the design-only examples (identical to those of the base model) are covered by the "
                         "analysis run above.")
        out.append("")
    out.append("## 6. The curve balls\n")
    category = None
    for name in names:
        c = cb.get(name)
        if c.category != category:
            category = c.category
            out.append(f"### {cb.CATEGORIES[category].capitalize()}\n")
        out.append(f"#### `{name}` -- {c.title}\n")
        out.append(f"* **Effect class**: {' or '.join(f'`{x}`' for x in c.effect_classes)}")
        counted = [x for cls in c.effect_classes for x in c.criteria_for(cls)]
        if counted:
            out.append("* **Material through**: " + ", ".join(f"`{x}`" for x in dict.fromkeys(counted)))
        out.append(f"* **What it changes**: {c.changes}")
        out.append(f"* **Why it matters**: {c.why}")
        out.append(f"* **Result quantities it moves**: {c.moves}")
        out.append(f"* **Applicable when**: {c.precondition}.")
        out.append("")
        shown = examples_of(runs, name)
        for entry, record in shown:
            out += example_block(model_key(*entry), entry, record, verified.get((model_key(*entry), name)))
            out.append("")
        if not shown:
            out.append("No example: not applicable to any of the models of section 5.\n")
        refusal = next((runs[model_key(*e)][name]["reason"] for e in BANK if "reason" in runs[model_key(*e)][name]
                        and not runs[model_key(*e)][name]["reason"].startswith("precondition")), None)
        if refusal:
            entry = next(e for e in BANK if runs[model_key(*e)][name].get("reason") == refusal)
            text = refusal if len(refusal) < 420 else refusal[:416] + " ..."
            out.append(f"Refused where it does not matter, e.g. on `{model_call(*entry)}`: *{text}*\n")
    out.append("## 7. Limits\n")
    out.append(LIMITS)
    return "\n".join(out).rstrip() + "\n"


LIMITS = """* A curve ball changes ONE thing (one target, or one family of ends for `release_pinned_frame`,
  `release_mechanism` and `second_order_off`); combinations of curve balls are obtained by applying
  them one after the other (each record refers to the document it was applied to).
* At most `MAX_CANDIDATES` targets are measured per call; a model in which only a later candidate
  would matter is reported as not applicable for that seed.
* Materiality is measured on the combinations of the file and on the result quantities listed in
  section 2 (nodal rotations and buckling factors are not part of the metric).
* `units_relabelled` is measured in physical units. The file is self-similar under a relabelling, so
  results in file units and design ratios do not change; a task built on it must state physical data
  outside the file.
* `section_lighter` / `section_heavier` take a section of the model, or a neighbour in the
  generators' generic catalogue when the member's section is a catalogue entry; they do not invent
  sections for hand-written models whose sections are neither.
* `offset_added` uses half the overall depth (`d`, `H`, `D` or `h`) of the deepest non-parallel frame
  member at the joint, whatever its orientation.
* `constraint_made_flexible` uses a section and the most common material of the model for the
  stand-in members; they carry self weight and are design members.
* Support curve balls leave alone the DOFs that are fixed at every node (the out-of-plane restraint
  of a plane frame in 3D space).
* Design ratios, pass / fail and the design status are those of `benchfea.design`; whatever the
  design checks do at a discontinuity (a change of the governing limit state, of a section class, of
  the lateral-torsional segment where a moment changes sign) carries over into `designRatioChange`
  and `designPassFailChanged`.
* The record's verification route is an expectation derived from the model's features. Whether a
  perturbed model's reference solution is verified is decided only by running `benchfea verify` on it.
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", help="directory for the generated base models and the runs (reused if present)")
    ap.add_argument("--check", action="store_true", help="do not write; exit 1 if docs/CURVEBALLS.md differs")
    ap.add_argument("--no-verify", action="store_true",
                    help="do not run the verification pipeline on the examples (the document then says nothing "
                         "about their verification status)")
    ap.add_argument("--budget", type=float, default=300.0, help="time budget per verification run, seconds")
    ap.add_argument("-o", "--output", default=str(ROOT / "docs" / "CURVEBALLS.md"))
    args = ap.parse_args(argv)
    text = build(Path(args.cache) if args.cache else None, verify=not args.no_verify, budget=args.budget)
    target = Path(args.output)
    if args.check:
        same = target.is_file() and target.read_text(encoding="utf-8") == text
        print("docs/CURVEBALLS.md is up to date" if same else "docs/CURVEBALLS.md DIFFERS from a fresh build")
        return 0 if same else 1
    target.write_text(text, encoding="utf-8")
    print(f"written {target} ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
