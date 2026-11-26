#!/usr/bin/env python3
"""Generate docs/VERIFICATION_REPORT.md from ACTUAL RUNS of the verification pipeline.

    python3 tools/build_verification_report.py                 # library + model bank, writes the report
    python3 tools/build_verification_report.py --no-bank       # the problem library only
    python3 tools/build_verification_report.py --json out.json # also write the machine-readable records

Nothing in the report is typed by hand: every number comes from ``benchfea.verify.verify_library`` (the
problems of verification/problems), ``verify_model`` on the model bank of tests/independent (the models
written for the independent review of the solver, stable and unstable) and the tier-F coverage registry.
Exit code 1 if an active problem fails or a problem directory is invalid.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "docs" / "VERIFICATION_REPORT.md"

from benchfea import verify as V  # noqa: E402
from benchfea.verify.record import DIAGNOSIS_TITLES, TIER_TITLES  # noqa: E402


def esc(text) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def num(x, fmt=".2e") -> str:
    return "n/a" if x is None else format(x, fmt)


def tol_text(t) -> str:
    return ", ".join(f"{k} {v:g}" for k, v in t.items()) if isinstance(t, dict) else str(t)


def source_text(src: dict) -> str:
    parts = [src.get("work", "")]
    for k in ("edition", "location", "pages"):
        if src.get(k):
            parts.append(f"{k}: {src[k]}")
    return "; ".join(p for p in parts if p)


def tier_state(t) -> str:
    if not t.applies:
        return "n/a"
    if t.skipped:
        return "not run"
    return "passed" if t.passed else f"NOT passed ({t.category or 'disagreement'})"


def largest_text(t, n=2) -> str:
    out = []
    for d in t.largest[:n]:
        v, tol = d.get("value"), d.get("tolerance")
        tol_s = tol_text(tol) if isinstance(tol, dict) else (num(tol) if isinstance(tol, (int, float)) else "")
        out.append(f"{esc(d.get('what', ''))}: {num(v) if isinstance(v, (int, float)) else v}"
                   + (f" (tol. {tol_s})" if tol_s else ""))
    return "; ".join(out)


def _target_kinds(r) -> dict:
    """Targets of a problem by kind, read from its problem.json: ``printed`` (a number of the source at its
    printed precision or at a tolerance the source states), ``byConstruction`` (a printed number the model
    was calibrated on: reproduced by construction, no evidence), ``recomputed`` (a number computed in this
    repository from the published statement), ``derived`` (derived here, tier A)."""
    out = {"printed": 0, "byConstruction": 0, "recomputed": 0, "derived": 0}
    try:
        data = json.loads((Path(r.directory) / "problem.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return out
    for t in data.get("targets", []):
        note = str(t.get("precisionNote", ""))
        if data.get("basis") == "derived-here":
            out["derived"] += 1
        elif "reproduced by construction" in note:
            out["byConstruction"] += 1
        elif note.startswith("evaluated here"):
            out["recomputed"] += 1
        elif data.get("basis") == "recomputed-from-published-statement" and not note.startswith(("printed", "stated")):
            out["recomputed"] += 1
        else:
            out["printed"] += 1
    return out


def evidence_section(report) -> list:
    """What the passed problems are evidence of: by source and by basis, with the kinds of targets."""
    L = ["### 1.1 What the problems compare against", "",
         "A problem that *passed* is not always a reproduction of a published number. Bases: `published` = every "
         "target is a number printed in the source; `recomputed-from-published-statement` = the problem STATEMENT "
         "is published, the targets were computed in this repository by means that share no code with the engine "
         "(ODE oracle, a separate small frame program, closed-form expressions evaluated here) -- evidence of the "
         "same kind as tier A, two computations made by this project -- plus those printed values of the source "
         "that reproduce; `derived-here` = tier A. Printed values of such a source that do NOT reproduce are kept "
         "as separate problems `<id>-vendor-printed`, quarantined. Targets *by construction* are printed numbers "
         "the model was calibrated on (no evidence).", "",
         "| source | basis | problems | passed | quarantined | failed | targets of the passed problems: printed in "
         "the source | by construction | recomputed here | derived here |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    rows, total = {}, {"printed": 0, "byConstruction": 0, "recomputed": 0, "derived": 0}
    only_recomputed = 0
    for r in report.problems:
        key = (Path(r.directory).parent.name, r.basis or "-")
        row = rows.setdefault(key, {"n": 0, "passed": 0, "quarantined": 0, "failed": 0, "invalid": 0,
                                    "printed": 0, "byConstruction": 0, "recomputed": 0, "derived": 0})
        row["n"] += 1
        row[r.outcome] = row.get(r.outcome, 0) + 1
        if r.outcome == "passed":
            kinds = _target_kinds(r)
            for k, v in kinds.items():
                row[k] += v
                total[k] += v
            if r.basis != "derived-here" and kinds["printed"] == 0:
                only_recomputed += 1
    for (slug, basis), row in sorted(rows.items()):
        L.append(f"| {slug} | {basis} | {row['n']} | {row['passed']} | {row['quarantined']} | "
                 f"{row['failed'] + row['invalid']} | {row['printed']} | {row['byConstruction']} | {row['recomputed']} | "
                 f"{row['derived']} |")
    n_pass = sum(1 for r in report.problems if r.outcome == "passed")
    L += ["", f"Targets of the {n_pass} passed problem(s): {total['printed']} printed in a source, "
          f"{total['byConstruction']} reproduced by construction, {total['recomputed']} recomputed here from a "
          f"published statement, {total['derived']} derived here. {only_recomputed} passed problem(s) with a "
          "published statement have no printed target at all (their evidence is the recomputation only).", ""]
    return L


def library_section(report) -> list:
    L = ["## 2. Problem library (`verification/problems`)", ""]
    if not report.problems:
        return L + ["No problem directory found.", ""]
    c = report.counts()
    L += [f"{len(report.problems)} problem(s): **{c['passed']} passed**, **{c['failed']} failed**, "
          f"**{c['quarantined']} quarantined**, **{c['invalid']} invalid**.", "",
          "Columns: *outcome* of the problem (targets + tiers that show a disagreement); *largest* = largest "
          "relative difference between a target and the engine; *model record* = status of `verify_model` on "
          "the problem's model (a model is `verified` only if every applicable tier passed; a model that fewer "
          "than two foreign solvers can represent is `quarantined (not-representable)`: it is not cross-checked "
          "and not publishable as verified, although its targets are reproduced).", "",
          "| problem | tier | basis | source | targets | largest | outcome | tiers passed | model record |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in report.problems:
        mr = r.model_record
        slug = Path(r.directory).parent.name
        status = "-" if mr is None else (mr.status + (f" ({mr.quarantine['category']})" if mr.quarantine else ""))
        L.append(f"| `{slug}/{r.id}` | {r.tier} | {r.basis} | {esc(source_text(r.source))} | {len(r.targets)} | "
                 f"{num(r.largest_relative)} | **{r.outcome}** | {''.join(mr.tiers_passed) if mr else '-'} | {status} |")
    L.append("")
    L += ["### 2.1 Targets and tiers, problem by problem", ""]
    for r in report.problems:
        slug = Path(r.directory).parent.name
        L += [f"#### `{slug}/{r.id}` — {esc(r.title)}", "",
              f"Outcome: **{r.outcome}** (declared `{r.declared_status}`). Source: {esc(source_text(r.source))}. "
              f"Basis: {r.basis}.", ""]
        for e in r.errors:
            L.append(f"* INVALID: {esc(e)}")
        if r.targets:
            L += ["| target | expected | engine | difference | relative | tolerance | |", "|---|---|---|---|---|---|---|"]
            for t in r.targets:
                L.append(f"| `{esc(t.label)}` | {t.expected:.9g} {esc(t.unit)} | "
                         f"{'no value' if t.actual is None else format(t.actual, '.9g')} | {num(t.difference)} | "
                         f"{num(t.relative)} | {tol_text(t.tolerance)} | {'ok' if t.passed else 'MISS ' + esc(t.error)} |")
            L.append("")
        if r.model_record is not None and r.model_record.tiers:
            for k, t in r.model_record.tiers.items():
                if t.applies:
                    lt = largest_text(t)
                    L.append(f"* tier {k} ({t.title or TIER_TITLES.get(k, DIAGNOSIS_TITLES.get(k, k))}): "
                             f"{tier_state(t)} — {esc(t.summary)}" + (f" Largest: {lt}." if lt else ""))
            L.append("")
        for f in r.failures:
            L.append(f"* {'recorded' if r.outcome == 'quarantined' else 'FAILURE'}: {esc(f)}")
        for n in r.notes:
            L.append(f"* note (tolerance policy): {esc(n)}")
        if r.limits_note:
            L.append(f"* limits: {esc(r.limits_note)}")
        if r.quarantine:
            L.append(f"* **quarantine diagnosis**: {esc(r.quarantine.get('diagnosis', ''))}")
        L.append("")
    return L


def quarantine_section(report, bank) -> list:
    L = ["## 3. Quarantined, failed and invalid cases", ""]
    rows = [r for r in report.problems if r.outcome != "passed"]
    if not rows:
        L += ["Problem library: none.", ""]
    for r in rows:
        slug = Path(r.directory).parent.name
        L.append(f"* **`{slug}/{r.id}`** — {r.outcome}. "
                 + (esc(r.quarantine.get("diagnosis", "")) if r.quarantine else "")
                 + (" Findings: " + esc("; ".join(r.failures + r.errors)) if (r.failures or r.errors) else ""))
    L.append("")
    not_cross = [r for r in report.problems if r.outcome == "passed" and r.model_record is not None
                 and r.model_record.status != "verified"]
    if not_cross:
        L += ["Problems whose targets are reproduced but whose MODEL is not publishable as verified "
              "(the record says why):", ""]
        for r in not_cross:
            q = r.model_record.quarantine or {}
            L.append(f"* `{Path(r.directory).parent.name}/{r.id}` — {q.get('category')}: {esc(q.get('diagnosis', ''))}")
        L.append("")
    if bank:
        L += ["Model bank: see section 4 (every model that is not `verified` is listed there with its diagnosis).", ""]
    return L


def bank_section(bank) -> list:
    L = ["## 4. Model bank of the independent review (`tests/independent`) through `verify_model`", "",
         "These models were written to review the solver (every modelling feature, alone and combined); they "
         "have no published targets. \"Independent\" means separate code written within this project; nobody outside "
         "the project has reviewed the solver. The table shows what the pipeline decides for each: `verified` = every "
         "applicable tier passed. `tableau` = largest relative difference to the independent tableau solver of "
         "`tests/independent` (EXTRA EVIDENCE ONLY: it is not one of the two foreign solvers tier C requires and "
         "does not influence the status).", "",
         "| model | kind | status | tiers passed | not passed | tableau | time [s] |", "|---|---|---|---|---|---|---|"]
    details = []
    for name, rec in bank:
        not_passed = ", ".join(f"{k} ({t.category or '-'})" for k, t in rec.tiers.items()
                               if t.applies and t.passed is not True)
        tab = ""
        if rec.kind == "results" and rec.tiers:
            ev = (rec.tiers["C"].details.get("extraEvidence") or {}).get("tableau")
            if isinstance(ev, dict):
                vals = [v.get("largestRelativeDifference") for v in ev.values() if isinstance(v, dict)]
                vals = [v for v in vals if v is not None]
                errs = [v.get("error") for v in ev.values() if isinstance(v, dict) and v.get("error")]
                tab = num(max(vals)) if vals else ("error" if errs or ev.get("error") else "")
        L.append(f"| `{name}` | {rec.kind} | **{rec.status}** | {''.join(rec.tiers_passed) or '-'} | {not_passed or '-'} | "
                 f"{tab or '-'} | {rec.timings.get('total', 0.0):.1f} |")
        if rec.status != "verified":
            why = (rec.quarantine or {}).get("diagnosis") or "; ".join(
                p.get("message", "") for p in (rec.refusal or {}).get("problems", [])[:2])
            details.append(f"* `{name}` — {rec.status}"
                           + (f" ({rec.quarantine['category']})" if rec.quarantine else "") + f": {esc(why)}")
    refined = {name: r for name, r in REFINED.items() if r}
    if refined:
        L += ["", "The tableau solver takes the axial force constant over a piece of a member, and the column above is "
              "its run with one piece per member. In a second-order model whose axial force varies along a member "
              "(self weight of columns, loads along a member) that number is the tableau's own discretisation error. "
              f"Every second-order model with a difference above {REFINE_ABOVE:g} was therefore run again with more "
              "pieces per member (extra evidence only; it changes no status):", "",
              "| model | 1 piece | " + " | ".join(f"{n} pieces" for n in REFINED_PIECES) + " |",
              "|---|---|" + "---|" * len(REFINED_PIECES)]
        for name, r in refined.items():
            L.append(f"| `{name}` | " + (r["error"] if "error" in r else
                                         " | ".join(num(r.get(n)) for n in (1,) + REFINED_PIECES)) + " |")
    n = len(bank)
    counts = {}
    for _, rec in bank:
        counts[rec.status] = counts.get(rec.status, 0) + 1
    L += ["", f"{n} model(s): " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) + ".", "",
          "### 4.1 Models that are not verified, with the diagnosis", ""] + (details or ["None."]) + [""]
    return L


def coverage_section() -> list:
    reg = V.coverage_registry()
    L = ["## 5. Tier F: coverage of the design limit states", "",
         f"Registry `src/benchfea/verify/design_coverage.json`, generated from the design tests by "
         f"`tools/build_design_coverage.py`. (b) = {esc(reg.get('designExamples', 'AISC Design Examples'))}. "
         "An item is *covered* when it has (a) a hand-calculation test and (b) a test in which a printed strength "
         "(or measured value / limit) of a check citing the item is reproduced; *partial* = one of the two is "
         "missing or only an intermediate value was reproduced. A design reference that uses an item that is not "
         "covered is not publishable as verified.", "",
         "| clause item | status | (a) tests | (b) level | (b) examples |", "|---|---|---|---|---|"]
    not_cov = []
    for item in reg.get("items", {}):
        c = V.item_coverage(item, reg)
        L.append(f"| {esc(item)} | {'covered' if c['status'] == 'covered' else '**' + c['status'] + '**'} | "
                 f"{len(c['aTests'])} | {c['b'] or '-'} | {', '.join(c['examples']) or '-'} |")
        if c["status"] != "covered":
            not_cov.append(item)
    L += ["", f"{len(reg.get('items', {}))} item(s), {len(not_cov)} without full coverage: "
          + (", ".join(not_cov) or "none") + ".", ""]
    return L


REFINED: dict = {}          # second-order bank model -> {pieces per member: largest difference tableau - engine}
REFINED_PIECES = (8, 32)
REFINE_ABOVE = 1e-6         # a difference above this, with one piece per member, is shown refined


def tableau_refined(builder) -> dict:
    """The tableau solver takes the axial force constant over a piece of a member.  Where the axial force
    varies along a member the one-piece result of the table is therefore the tableau's own discretisation
    error; this reruns the comparison with more pieces per member (extra evidence only)."""
    import harness
    raw, cid = builder()
    out = {}
    try:
        for n in (1,) + REFINED_PIECES:
            d = harness.compare(raw, cid, second_order=True, nsub=n)
            out[n] = max(v for k, v in d.items() if not k.startswith("_"))
            if n == 1 and out[n] <= REFINE_ABOVE:
                return {}
    except Exception as e:                                 # extra evidence must never break the report
        return {"error": f"{type(e).__name__}: {str(e)[:120]}"}
    return out


def run_bank(budget):
    sys.path.insert(0, str(ROOT / "tests" / "independent"))
    try:
        from models import LINEAR, ONE_SIDED, SECOND
        from unstable import UNSTABLE
    except ImportError:
        return []
    out = []
    for name, f in {**LINEAR, **SECOND, **ONE_SIDED, **UNSTABLE}.items():
        raw = f()
        raw = raw[0] if isinstance(raw, tuple) else raw
        rec = V.verify_model(raw, budget=budget, extra_evidence=True, name=name)
        if name in SECOND:
            REFINED[name] = tableau_refined(f)
        out.append((name, rec))
        print(f"  bank {name}: {rec.status} [{''.join(rec.tiers_passed)}] {rec.timings.get('total', 0):.1f} s", flush=True)
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    budget = 600.0
    t0 = time.perf_counter()
    report = V.verify_library(budget=budget, progress=lambda r: print(f"  {r.outcome:<11} {r.id} ({r.seconds:.1f} s)",
                                                                       flush=True))
    bank = [] if "--no-bank" in argv else run_bank(budget)
    versions = report.versions
    c = report.counts()
    L = ["# Verification report", "",
         "Generated by `tools/build_verification_report.py` from actual runs; do not edit by hand. "
         "Policy, tiers and tolerances: [VERIFICATION.md](VERIFICATION.md).", "",
         f"> {V.MEANING}", "",
         "## 1. Summary", "",
         "Software: " + ", ".join(f"{k} {v or 'not installed'}" for k, v in versions.items()) + ".", "",
         f"* Problem library: {len(report.problems)} problem(s) — {c['passed']} passed, {c['failed']} failed, "
         f"{c['quarantined']} quarantined, {c['invalid']} invalid.",
         f"* Of the {c['passed']} passed problem(s), "
         f"{sum(1 for r in report.problems if r.outcome == 'passed' and r.model_record and r.model_record.status == 'verified')} "
         "have a model whose record is `verified` (every applicable tier passed, including the cross-check by two "
         "foreign solvers)."]
    if bank:
        n_v = sum(1 for _, r in bank if r.status == "verified")
        L.append(f"* Model bank: {len(bank)} model(s) — {n_v} verified, "
                 f"{sum(1 for _, r in bank if r.status == 'quarantined')} quarantined, "
                 f"{sum(1 for _, r in bank if r.status == 'refused')} refused by the engine.")
    L += [f"* Run time of this report: {time.perf_counter() - t0:.0f} s (single BLAS thread).", ""]
    L += evidence_section(report)
    L += library_section(report) + quarantine_section(report, bank)
    if bank:
        L += bank_section(bank)
    else:
        L += ["## 4. Model bank", "", "Not run (`--no-bank`).", ""]
    L += coverage_section()
    L += ["## 6. Tolerances in force", "", "```json", json.dumps(V.POLICY, indent=1, default=list), "```", ""]
    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT}")
    if "--json" in argv:
        path = Path(argv[argv.index("--json") + 1])
        path.write_text(json.dumps({"library": report.to_dict(),
                                    "bank": {n: r.to_dict(reference=False) for n, r in bank}}, indent=1), encoding="utf-8")
        print(f"wrote {path}")
    return 0 if report.ok else 1


if __name__ == "__main__":
    sys.exit(main())
