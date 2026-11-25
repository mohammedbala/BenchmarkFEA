"""Command line interface: ``benchfea <command> ...``."""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__

_COMMANDS = {
    "validate": ("Check an input file against the schema and the semantic rules; list every problem.",
                 "benchfea validate examples/03_portal_frame.json\n  benchfea validate model.yaml --json"),
    "run": ("Analyse a model and write the results document.",
            "benchfea run examples/03_portal_frame.json -o portal.results.json"),
    "design": ("Analyse a model and run the AISC 360-22 member checks.",
               "benchfea design examples/03_portal_frame.json -o portal.design.json"),
    "perturb": ("Apply a named curve ball to a model: a seeded, deterministic change that matters (measured on "
                "analysis and design), with a machine-checkable change record (docs/CURVEBALLS.md).",
                "benchfea perturb --list\n"
                "  benchfea perturb examples/03_portal_frame.json release_added --seed 1 -o variant.json "
                "--record variant.record.json\n"
                "  benchfea perturb examples/03_portal_frame.json support_pinned_to_fixed --show-record\n"
                "  benchfea perturb examples/03_portal_frame.json --applicable"),
    "verify": ("Verify a reference solution: run every verification tier that applies to a model file, to one "
               "problem of the verification library, or to the whole library (--all). 'verified' means only: "
               "passed the tiers listed, at the tolerances stated.",
               "benchfea verify examples/03_portal_frame.json\n  benchfea verify examples/03_portal_frame.json --design --json\n"
               "  benchfea verify verification/problems/derived/beam-point-load\n  benchfea verify --all"),
    "score": ("Score an answer file against a task of the benchmark and its verified reference (the hidden "
              "answers file). Prints the score in [0, 1] and the machine-readable reasons (docs/BENCHMARK.md).",
              "benchfea score answer.json --task data/tasks/train/<id>.task.json\n"
              "  benchfea score answer.json --task my.task.json --answers my.answers.json --json"),
    "bench": ("Benchmark helpers: build the dataset of verified tasks, run an agent in the environment.",
              "benchfea bench make-dataset --out data --seeds 0,1,4 --workers 4 --cache /tmp/benchfea-cache\n"
              "  benchfea bench run-agent examples/agents/oracle_agent.py --dataset data --episodes 3 --seed 1\n"
              "  benchfea bench run-agent examples/agents/random_agent.py --task-type predict --mode tool-allowed"),
    "generate": ("Generate a benchmark model of a family for a seed and a difficulty level (deterministic; "
                 "docs/DIFFICULTY.md).",
                 "benchfea generate portal_frame --seed 7 --difficulty 2 -o portal.json --metadata portal.meta.json\n"
                 "  benchfea generate --list"),
}


def _validate(args) -> int:
    from .io import load
    from .problems import ERROR, ModelInvalid
    from .validate import validate

    try:
        doc = load(args.file)
        problems = validate(doc)
    except ModelInvalid as e:
        problems = e.problems
    errors = [p for p in problems if p.severity == ERROR]
    if args.json:
        json.dump({"valid": not errors, "problems": [p.to_dict() for p in problems]}, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        for p in problems:
            print(p.format(args.file))
        print(f"{args.file}: " + ("valid" if not errors else "INVALID")
              + f" ({len(errors)} error(s), {len(problems) - len(errors)} warning(s))")
    if not errors and args.defaults:
        from .model import build_model
        json.dump([x.to_dict() for x in build_model(args.file).defaults_used], sys.stdout, indent=2)
        sys.stdout.write("\n")
    return 1 if errors else 0


def _run(args) -> int:
    """Analyse a model, write the results document, print a short summary."""
    from pathlib import Path

    from .io import dump
    from .model import build_model
    from .problems import ModelInvalid, ModelUnstable
    from .solve import analyze

    out_path = Path(args.output) if args.output else Path(args.file).with_suffix(".results.json")
    try:
        model = build_model(args.file)
        results = analyze(model, strict=not args.allow_flagged)
    except ModelInvalid as e:
        diagnosis = getattr(e, "diagnosis", None)
        if args.json:
            payload = {"ok": False, "refusal": "unstable" if isinstance(e, ModelUnstable) else "invalid",
                       "problems": [p.to_dict() for p in e.problems]}
            if diagnosis is not None:
                payload["diagnosis"] = diagnosis.to_dict()
            json.dump(payload, sys.stdout, indent=2)
            sys.stdout.write("\n")
        else:
            for p in e.problems:
                print(p.format(args.file), file=sys.stderr)
            what = "UNSTABLE (mechanism)" if isinstance(e, ModelUnstable) else "REFUSED"
            print(f"{args.file}: {what}, no results written ({len(e.errors)} error(s))", file=sys.stderr)
        # exit code 2: the input could not be read at all; 1: the model was read and refused
        return 2 if e.problems and all(p.code == "parse.io" for p in e.problems) else 1
    doc = results.to_document(model)
    dump(doc, out_path)
    u = model.units
    # flagged: an invariant failed or an iteration did not converge; with --allow-flagged a combination
    # that could not be analysed at all (e.g. beyond its elastic critical load) is missing from the results
    flagged = [cid for cid, c in results.combinations.items() if not (c.invariants.passed and c.converged)]
    flagged += [cid for cid in model.combinations if cid not in results.combinations]
    flagged += [f"buckling:{cid}" for cid in model.analysis.buckling_combinations if cid not in results.buckling]
    summary = {"ok": not flagged, "output": str(out_path), "model": model.name, "nodes": len(model.nodes),
               "members": len(model.members), "combinations": {}, "warnings": len(results.warnings),
               "defaultsUsed": len(results.defaults_used), "flagged": flagged}
    for cid, c in results.combinations.items():
        t = abs(c.displacements[:, :3])
        k, a = divmod(int(t.argmax()), 3)
        worst = max(c.invariants.checks, key=lambda ch: (not ch.passed, ch.value / ch.tolerance
                                                         if ch.name != "positive definiteness" else 0.0))
        mid, mx, mv = max(((m, float(r.x[int(abs(r.displacements).max(axis=1).argmax())]),
                            float(abs(r.displacements).max())) for m, r in c.members.items()),
                          key=lambda t: t[2], default=("", 0.0, 0.0))
        summary["combinations"][cid] = {
            "analysisType": c.analysis_type,
            "maxDisplacement": {"node": list(model.nodes)[k], "dof": ("ux", "uy", "uz")[a],
                                "value": float(u.from_internal(c.displacements[k, a], "length"))},
            "maxStationDisplacement": {"member": mid, "x": float(u.from_internal(mx, "length")),
                                       "value": float(u.from_internal(mv, "length"))},
            "sumReactions": [float(u.from_internal(v, "force")) for v in c.reactions[:, :3].sum(axis=0)],
            "invariantsPassed": c.invariants.passed,
            "worstInvariant": {"name": worst.name, "value": float(worst.value), "tolerance": float(worst.tolerance)}}
        so = c.solver.get("secondOrder")
        if so:                                   # record of the second-order iteration (docs/THEORY.md section 10)
            last = so["history"][-1] if so["history"] else {}
            summary["combinations"][cid]["secondOrder"] = {
                "converged": bool(so["converged"]), "iterations": int(so["iterations"]),
                "tolerance": so["tolerance"], "lastDisplacementChange": last.get("displacementChange"),
                "lastAxialForceChange": last.get("axialForceChange"), "alpha": float(c.alpha),
                "piecesPerMember": {"min": min(so["subdivisions"].values()), "max": max(so["subdivisions"].values())}}
        if "activeSet" in c.solver:
            summary["combinations"][cid]["slackMembers"] = list(c.solver["activeSet"]["slack"])
        if c.tau_b:
            summary["combinations"][cid]["tauB"] = {"min": min(c.tau_b.values()), "max": max(c.tau_b.values())}
    if results.buckling:
        summary["buckling"] = {
            cid: {"factors": [float(f) for f in b.factors], "negativeFactors": [float(f) for f in b.negative_factors],
                  "estimatedRelativeError": list(b.info.get("estimatedRelativeError", [])),
                  "piecesPerMember": b.info.get("piecesPerMember"), "warnings": list(b.info.get("warnings", []))}
            for cid, b in results.buckling.items()}
    if args.json:
        json.dump(summary, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        print(f"{args.file}: {len(model.nodes)} nodes, {len(model.members)} members, "
              f"{len(results.combinations)} combination(s); units {u.force}, {u.length}")
        for cid, c in summary["combinations"].items():
            d, r, m = c["maxDisplacement"], c["sumReactions"], c["maxStationDisplacement"]
            print(f"  {cid} [{c['analysisType']}]: largest node displacement {d['value']:.6g} {u.label('length')} "
                  f"({d['dof']} at {d['node']}); largest member-axis displacement component {m['value']:.6g} "
                  f"{u.label('length')} (member {m['member']}, x = {m['x']:.6g}); sum of reactions ({r[0]:.6g}, {r[1]:.6g}, {r[2]:.6g}) "
                  f"{u.label('force')}; invariants {'passed' if c['invariantsPassed'] else 'FAILED'} "
                  f"(largest: {c['worstInvariant']['name']} {c['worstInvariant']['value']:.2e})")
            so = c.get("secondOrder")
            if so:
                p = so["piecesPerMember"]
                print(f"      second-order iteration {'converged' if so['converged'] else 'DID NOT CONVERGE'} in "
                      f"{so['iterations']} iteration(s): last displacement change {so['lastDisplacementChange']:.2e}, "
                      f"axial-force change {so['lastAxialForceChange']:.2e} (tolerance {so['tolerance']:.1e}); "
                      f"alpha {so['alpha']:g}; pieces per member {p['min']}..{p['max']}")
            if c.get("slackMembers"):
                print(f"      slack tension/compression-only members: {', '.join(c['slackMembers'])}")
            if "tauB" in c:
                print(f"      stiffness reduction: tau_b {c['tauB']['min']:.4g}..{c['tauB']['max']:.4g}")
        for cid, b in summary.get("buckling", {}).items():
            pos = ", ".join(f"{f:.6g}" for f in b["factors"]) or "none (no member in compression)"
            line = f"  buckling of combination {cid}: factors {pos}"
            if b["negativeFactors"]:
                line += "; under reversed loads " + ", ".join(f"{f:.6g}" for f in b["negativeFactors"])
            if b["estimatedRelativeError"]:
                line += f" (estimated relative error {max(b['estimatedRelativeError']):.1e})"
            print(line)
            for w in b["warnings"]:
                print(f"      note: {w}")
        for w in results.warnings:
            print(f"  warning: {w}")
        print(f"  {len(results.defaults_used)} default(s) used (listed in the output under defaultsUsed)")
        print(f"results written to {out_path}"
              + (" -- FLAGGED (invariant failed, not converged or no results): " + ", ".join(flagged) if flagged else ""))
    return 3 if flagged else 0


def _design(args) -> int:
    """Analyse a model, run the member design checks, write the design document, print one row per
    design member: status, governing limit state, clause, combination, ratio.

    Exit code: 0 design status "pass" or "conditional"; 4 "fail" or "incomplete"; 1 model refused;
    2 input unreadable / usage; 3 analysis flagged (invariants failed; the design is still written)."""
    from pathlib import Path

    from . import design
    from .io import dump
    from .model import build_model
    from .problems import ModelInvalid, ModelUnstable
    from .solve import analyze

    out_path = Path(args.output) if args.output else Path(args.file).with_suffix(".design.json")
    try:
        model = build_model(args.file)
        results = analyze(model, strict=not args.allow_flagged)
    except ModelInvalid as e:
        if args.json:
            json.dump({"ok": False, "refusal": "unstable" if isinstance(e, ModelUnstable) else "invalid",
                       "problems": [p.to_dict() for p in e.problems]}, sys.stdout, indent=2)
            sys.stdout.write("\n")
        else:
            for p in e.problems:
                print(p.format(args.file), file=sys.stderr)
            print(f"{args.file}: REFUSED, no design written ({len(e.errors)} error(s))", file=sys.stderr)
        return 2 if e.problems and all(p.code == "parse.io" for p in e.problems) else 1
    out = design.check(model, results)
    dump(out.to_document(model), out_path)
    if args.results:
        dump(results.to_document(model), Path(args.results))
    flagged = [cid for cid, c in results.combinations.items() if not c.invariants.passed]
    rows = []
    for mid, m in out.members.items():
        g = m.governing
        rows.append({"member": mid, "status": m.status, "governing": g.name if g else None,
                     "clause": g.clause if g else None, "combination": g.combination if g else None,
                     "ratio": float(g.ratio) if g else None, "notChecked": list(m.not_checked),
                     "warnings": list(m.warnings)})
    summary = {"ok": out.status in (design.PASS, design.CONDITIONAL) and not flagged, "status": out.status,
               "code": out.code, "output": str(out_path), "model": model.name, "members": rows,
               "combinations": out.combinations, "warnings": list(out.warnings), "flagged": flagged}
    if args.json:
        json.dump(summary, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        print(f"{args.file}: {len(out.members)} design member(s), code {out.code or 'n/a'}")
        width = max([6] + [len(r["member"]) for r in rows])
        print(f"  {'member':<{width}}  {'status':<11}  {'ratio':>7}  {'comb.':<8}  governing limit state [clause]")
        for r in rows:
            ratio = "   n/a " if r["ratio"] is None else f"{r['ratio'] + 0.0:7.3f}"
            print(f"  {r['member']:<{width}}  {r['status']:<11}  {ratio}  {(r['combination'] or '-'):<8}  "
                  f"{r['governing'] or 'no limit state could be checked'}"
                  + (f" [{r['clause']}]" if r["clause"] else ""))
            for text in r["notChecked"]:
                if ": notChecked: " in text or ": refused: " in text:
                    print(f"  {'':<{width}}    not checked -- {text}")
        for w in out.warnings:
            print(f"  warning: {w}")
        n_member_warnings = sum(len(r["warnings"]) for r in rows)
        if n_member_warnings:
            print(f"  {n_member_warnings} member warning(s) are listed in the output under members/<id>/warnings")
        print(f"design status: {out.status.upper()} -- written to {out_path}"
              + (" -- analysis FLAGGED: invariants failed" if flagged else ""))
    if flagged:
        return 3
    return 0 if out.status in (design.PASS, design.CONDITIONAL) else 4


def _verify(args) -> int:
    """Run the verification pipeline (docs/VERIFICATION.md).

    Exit code: 0 the reference solution is verified / every active problem passed; 5 a model is
    quarantined (an applicable tier is not passed) or an active problem fails; 1 the model was refused
    (invalid or unsupported: there is no reference solution); 2 usage / nothing to verify."""
    import os
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):      # small dense matrices: threads only cost time
        os.environ.setdefault(var, "1")
    from pathlib import Path

    from . import verify as V

    def emit(payload):
        json.dump(payload, sys.stdout, indent=2)
        sys.stdout.write("\n")

    def problem_line(r):
        tiers = "" if r.model_record is None else "".join(r.model_record.tiers_passed)
        rel = r.largest_relative
        cross = {True: "yes", False: "NO", None: "-"}[r.cross_checked]
        return (f"  {r.outcome.upper():<11} {Path(r.directory).parent.name}/{r.id}: {len(r.targets)} target(s), largest "
                f"relative difference {'n/a' if rel is None else format(rel, '.2e')}; tiers passed: {tiers or '-'}; "
                f"cross-checked by two solvers: {cross}")

    if args.all or (args.target and Path(args.target).is_dir() and not (Path(args.target) / "problem.json").is_file()):
        root = args.target if args.target else None
        if root is None and V.library_root() is None:
            print("benchfea verify --all: no verification/problems directory found from here", file=sys.stderr)
            return 2
        report = V.verify_library(root, tiers=args.tiers, budget=args.budget,
                                  progress=None if args.json else (lambda r: print(problem_line(r), flush=True)))
        if args.json:
            emit(report.to_dict())
        else:
            for r in report.problems:
                for text in r.failures + r.errors:
                    print(f"  {r.id}: {text}")
            c = report.counts()
            print(f"{len(report.problems)} problem(s): {c['passed']} passed, {c['failed']} FAILED, "
                  f"{c['quarantined']} quarantined, {c['invalid']} invalid ({report.seconds:.1f} s)")
            print(V.MEANING)
        if not report.problems:
            return 2
        return 0 if report.ok else 5
    if not args.target:
        print("benchfea verify: give a model file, a problem directory, or --all", file=sys.stderr)
        return 2
    target = Path(args.target)
    if target.is_dir():
        r = V.verify_problem(target, tiers=args.tiers, budget=args.budget)
        if args.json:
            emit(r.to_dict(reference=args.reference))
        else:
            print(problem_line(r))
            for t in r.targets:
                print(f"    {'ok  ' if t.passed else 'MISS'} {t.label}: expected {t.expected:.9g}, engine "
                      + ("no value" if t.actual is None else f"{t.actual:.9g}")
                      + (f" ({t.error})" if t.error else f" (difference {t.difference:.2e}, tolerance {t.tolerance})"))
            for text in r.failures + r.errors + r.notes:
                print(f"    {text}")
            if r.quarantine:
                print(f"    quarantined: {r.quarantine.get('diagnosis')}")
            if r.model_record is not None:
                print("    " + r.model_record.summary().replace("\n", "\n    "))
        return 0 if r.outcome in ("passed", "quarantined") else (2 if r.outcome == "invalid" else 5)
    rec = V.verify_model(str(target), tiers=args.tiers, design=args.design, budget=args.budget,
                         extra_evidence=args.extra_evidence)
    if args.output:
        from .io import dump
        dump(rec.to_dict(), Path(args.output))
    if args.json:
        emit(rec.to_dict(reference=args.reference))
    else:
        print(rec.summary())
        print(f"  time: {rec.timings.get('total', 0.0):.1f} s; " + V.MEANING)
    return {V.VERIFIED: 0, V.QUARANTINED: 5, V.REFUSED: 1}[rec.status]


def _generate(args) -> int:
    """Write one generated model (canonical JSON) and optionally its metadata.  Exit code 0 ok, 2 usage,
    6 no acceptable model for these arguments."""
    import os
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):      # small dense matrices: threads only cost time
        os.environ.setdefault(var, "1")
    from . import bench

    if args.list:
        for name in bench.families():
            fam = bench.FAMILIES[name]
            print(f"{name}: levels {', '.join(str(lv) for lv in fam.LEVELS)}; features: {', '.join(fam.FEATURES)}")
        return 0
    if not args.family:
        print("benchfea generate: give a family (see --list)", file=sys.stderr)
        return 2
    features = {}
    for item in args.feature or []:
        name, _, value = item.partition("=")
        if value.lower() not in ("true", "false", "on", "off", ""):
            print(f"benchfea generate: --feature {item!r}: use NAME=true or NAME=false", file=sys.stderr)
            return 2
        features[name] = value.lower() in ("true", "on", "")
    try:
        g = bench.generate(args.family, args.seed, args.difficulty, units=args.units, vertical_axis=args.vertical_axis,
                           features=features)
    except ValueError as e:
        print(f"benchfea generate: {e}", file=sys.stderr)
        return 2
    except bench.GenerationFailed as e:
        print(f"benchfea generate: {e}", file=sys.stderr)
        return 6
    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(g.json())
    else:
        sys.stdout.write(g.json())
    if args.metadata:
        with open(args.metadata, "w", encoding="utf-8") as fh:
            json.dump(g.metadata, fh, indent=2)
            fh.write("\n")
    m = g.metadata
    print(f"{args.family} seed {args.seed}: level {m['measured']['level']} (score {m['measured']['score']}), "
          f"{m['measured']['nMembers']} member(s), units {m['units']}, vertical axis {m['verticalAxis']}, design ratios "
          f"{m['design']['ratioMin']} .. {m['design']['ratioMax']}, cross-checkable (expected): "
          f"{m['tierC_expected']['crossCheckable']}", file=sys.stderr)
    return 0


def _perturb(args) -> int:
    """Apply a curve ball.  Exit code 0 ok; 2 usage, unknown curve ball or unreadable input; 7 the curve
    ball is not applicable to this model (the reason is printed)."""
    import os
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):      # small dense matrices: threads only cost time
        os.environ.setdefault(var, "1")
    from .bench import curveballs as cb
    from .problems import ModelInvalid

    if args.list:
        if args.json:
            json.dump([{"name": n, "category": cb.get(n).category, "title": cb.get(n).title,
                        "effectClasses": list(cb.get(n).effect_classes), "changes": cb.get(n).changes,
                        "whyItMatters": cb.get(n).why, "moves": cb.get(n).moves,
                        "precondition": cb.get(n).precondition} for n in cb.available()], sys.stdout, indent=2)
            sys.stdout.write("\n")
            return 0
        category = None
        for n in cb.available():
            c = cb.get(n)
            if c.category != category:
                category = c.category
                print(f"{category} -- {cb.CATEGORIES.get(category, category)}")
            print(f"  {n:<32} {c.title} [{' | '.join(c.effect_classes)}]")
        print("effect classes: " + ", ".join(cb.EFFECT_CLASSES) + "; materiality: analysis change or design ratio "
              f"change >= {cb.THRESHOLDS['analysisChange']:.0%}, or a design member changes pass / fail "
              "(docs/CURVEBALLS.md)")
        return 0
    if not args.file:
        print("benchfea perturb: give a model file and a curve ball name (see --list)", file=sys.stderr)
        return 2
    try:
        if args.applicable:
            table = cb.applicability(args.file, seed=args.seed)
            if args.json:
                json.dump(table, sys.stdout, indent=2)
                sys.stdout.write("\n")
            else:
                for n, row in table.items():
                    print(f"  {n:<32} " + (f"applicable [{row['effectClass']}]: {row['target']}" if row["applicable"]
                                           else f"not applicable: {row['reason']}"))
            return 0
        if not args.name:
            print("benchfea perturb: give the name of a curve ball (see --list) or --applicable", file=sys.stderr)
            return 2
        try:
            result = cb.apply(args.file, args.name, args.seed, effect_class=args.effect_class,
                              verify={None: False, "auto": True}.get(args.require_verified, args.require_verified))
        except KeyError as e:
            print(f"benchfea perturb: {e.args[0]}", file=sys.stderr)
            return 2
        except cb.NotApplicable as e:
            if args.json:
                json.dump({"ok": False, "curveBall": e.name, "notApplicable": e.reason, "tried": e.tried},
                          sys.stdout, indent=2)
                sys.stdout.write("\n")
            else:
                print(f"benchfea perturb: {e}", file=sys.stderr)
            return 7
    except ModelInvalid as e:
        for p in e.problems:
            print(p.format(args.file), file=sys.stderr)
        print(f"{args.file}: cannot be read as a model", file=sys.stderr)
        return 2
    rec = result.record
    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(result.json())
    if args.record:
        with open(args.record, "w", encoding="utf-8") as fh:
            fh.write(rec.json())
    if args.show_record or args.json:
        sys.stdout.write(rec.json())
    elif not args.output:
        sys.stdout.write(result.json())
    e = rec.effect
    line = f"{rec.name} seed {rec.seed}: {rec.target}; effect class {rec.effect_class}"
    if "analysisChange" in e:
        line += f"; analysis change {e['analysisChange']:.1%}"
        if e.get("designRatioChange") is not None:
            line += f", design ratio change {e['designRatioChange']:.1%}"
    line += f"; {len(rec.changes)} change(s): " + ", ".join(rec.paths[:4]) + (" ..." if len(rec.changes) > 4 else "")
    print(line, file=sys.stderr)
    ran = rec.verification.get("run")
    print(f"reference of the perturbed model: {rec.reference['kind']}; "
          + (f"verification run ({ran['mode']}): {ran['status']}, tiers passed {', '.join(ran['tiersPassed'])}" if ran
             else f"expected verification route (not run): {rec.verification['route']}"), file=sys.stderr)
    return 0


def _score(args) -> int:
    """Score an answer file.  Exit code 0: scored (whatever the score; a malformed answer scores 0 with the
    reason); 2: the task, the answers file or the answer file cannot be read."""
    import os
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(var, "1")
    from pathlib import Path

    from .bench import dataset
    from .bench.scoring import score_answer

    try:
        task = dataset.load_task(args.task)
        if args.answers:
            answers = json.loads(Path(args.answers).read_text(encoding="utf-8"))
        else:
            answers = dataset.load_answers(args.task)
    except (OSError, ValueError, KeyError) as e:
        print(f"benchfea score: cannot read the task or its answers file ({e}); give --answers if the task is not "
              "in the dataset layout", file=sys.stderr)
        return 2
    try:
        text = Path(args.answer).read_text(encoding="utf-8")
    except OSError as e:
        print(f"benchfea score: cannot read {args.answer}: {e}", file=sys.stderr)
        return 2
    if task.get("id") != answers.get("id"):
        print("benchfea score: the answers file belongs to another task", file=sys.stderr)
        return 2
    out = score_answer(task, answers, text, budget=args.budget)
    if args.output:
        Path(args.output).write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    if args.json:
        json.dump(out, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        print(f"task {task.get('id')} [{task.get('taskType')}, {task.get('mode')}]: score {out['score']:.4f}"
              + ("" if out["valid"] else " (answer not valid)"))
        print(f"  {out['summary']}")
        for qid, row in out["reasons"].get("quantities", {}).items():
            err = row.get("error") or f"relative error {row['relativeError']:.3e}"
            print(f"  {qid}: credit {row['credit']:.3f} ({err}; full credit up to {row['tolerance']['full']:.0%})")
        failed = out["reasons"].get("failed")
        if failed:
            print("  not met: " + ", ".join(str(x) for x in failed))
    return 0


def _bench(args) -> int:
    """``benchfea bench make-dataset | run-agent``.  Exit code 0 ok; 2 usage; 8 no verified task."""
    import os
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(var, "1")
    if args.bench_command == "make-dataset":
        from .bench import dataset
        try:
            cells = dataset.plan(seeds=tuple(int(x) for x in args.seeds.split(",")),
                                 families=args.families.split(",") if args.families else None,
                                 curve_balls=args.curve_balls.split(",") if args.curve_balls else None)
        except (ValueError, KeyError) as e:
            print(f"benchfea bench make-dataset: {e.args[0]}", file=sys.stderr)
            return 2

        def progress(row, i, n):
            print(f"[{i}/{n}] {row['family']} / {row['curveBall']} seed {row['seed']} level {row['difficulty']}: "
                  f"{row['status']}, {len(row['tasks'])} task(s), {row['seconds']} s", file=sys.stderr, flush=True)

        m = dataset.build(args.out, cells, workers=args.workers, cache_dir=args.cache, clean=args.clean,
                          progress=None if args.quiet else progress)
        c = m["counts"]
        print(f"{c['tasks']} verified task(s) from {c['cellsVerified']} of {c['cells']} cell(s); by split "
              f"{c['bySplit']}; by task type {c['byTaskType']}; {len(m['quarantineLog'])} entr(ies) in the quarantine "
              f"log; manifest: {args.out}/tasks/manifest.json")
        return 0
    if args.bench_command == "run-agent":
        import importlib.util

        try:
            import gymnasium

            from .bench import env as env_module
        except ImportError:
            print("benchfea bench run-agent: gymnasium is not installed (pip install benchfea[env])", file=sys.stderr)
            return 2
        from .bench.tasks import TaskUnavailable
        path, _, name = args.agent.partition(":")
        try:
            spec = importlib.util.spec_from_file_location("benchfea_agent", path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        except (OSError, AttributeError, ImportError, SyntaxError) as e:
            print(f"benchfea bench run-agent: cannot load {path}: {e}", file=sys.stderr)
            return 2
        factory = getattr(module, name or "make_agent", None)
        if factory is None:
            print(f"benchfea bench run-agent: {path} has no function {name or 'make_agent'}(env)", file=sys.stderr)
            return 2
        env = gymnasium.make(env_module.ENV_ID, split=args.split, max_steps=args.max_steps, cache_dir=args.cache,
                             dataset_dir=args.dataset)
        options = {k: v for k, v in (("family", args.family), ("curve_ball", args.curve_ball),
                                     ("difficulty", args.difficulty), ("task_type", args.task_type),
                                     ("mode", args.mode)) if v is not None}
        rows = []
        for k in range(args.episodes):
            try:
                ep = env_module.run_episode(env, factory(env), seed=args.seed + k, options=options)
            except (TaskUnavailable, ValueError) as e:
                print(f"benchfea bench run-agent: episode {k}: {e}", file=sys.stderr)
                return 8
            info = ep["info"]
            rows.append({"episode": k, "seed": args.seed + k, "task": info["task_id"], "taskType": info["task_type"],
                         "mode": info["mode"], "family": info["family"], "reward": ep["reward"], "steps": ep["steps"],
                         "terminated": ep["terminated"], "truncated": ep["truncated"],
                         "summary": (info.get("score") or {}).get("summary", "no answer submitted")})
            if not args.json:
                r = rows[-1]
                print(f"episode {k} (seed {r['seed']}): {r['family']} {r['taskType']} [{r['mode']}] -> reward "
                      f"{r['reward']:.3f} in {r['steps']} step(s){' (truncated)' if r['truncated'] else ''}; {r['summary']}")
        mean = sum(r["reward"] for r in rows) / len(rows) if rows else 0.0
        if args.json:
            json.dump({"episodes": rows, "meanReward": mean}, sys.stdout, indent=2)
            sys.stdout.write("\n")
        else:
            print(f"mean reward over {len(rows)} episode(s): {mean:.3f}")
        return 0
    print("benchfea bench: give a sub-command: make-dataset or run-agent", file=sys.stderr)
    return 2


def _not_implemented(args) -> int:
    print(f"benchfea {args.command}: not available in this build (module not implemented yet)", file=sys.stderr)
    return 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="benchfea", description="3D frame analysis engine and structural-modelling benchmark.",
        epilog="Run 'benchfea <command> -h' for the options and an example of each command.")
    parser.add_argument("--version", action="version", version=f"benchfea {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="command")
    for name, (text, example) in _COMMANDS.items():
        p = sub.add_parser(name, help=text, description=text, epilog="example:\n  " + example,
                           formatter_class=argparse.RawDescriptionHelpFormatter)
        if name == "validate":
            p.add_argument("file", help="input file (.json, .yaml, .yml)")
            p.add_argument("--json", action="store_true", help="print the problems as JSON")
            p.add_argument("--defaults", action="store_true",
                           help="if valid, also print the defaults-used report (JSON)")
            p.set_defaults(func=_validate)
        elif name == "run":
            p.add_argument("file", help="input file (.json, .yaml, .yml)")
            p.add_argument("-o", "--output", help="results file (.json or .yaml); default: <input>.results.json")
            p.add_argument("--json", action="store_true", help="print the summary (or the refusal) as JSON")
            p.add_argument("--allow-flagged", action="store_true",
                           help="write results even if an invariant fails (they are flagged; exit code 3)")
            p.set_defaults(func=_run)
        elif name == "design":
            p.add_argument("file", help="input file (.json, .yaml, .yml)")
            p.add_argument("-o", "--output", help="design file (.json or .yaml); default: <input>.design.json")
            p.add_argument("--results", help="also write the analysis results document to this file")
            p.add_argument("--json", action="store_true", help="print the summary (or the refusal) as JSON")
            p.add_argument("--allow-flagged", action="store_true",
                           help="design even if an analysis invariant fails (exit code 3)")
            p.set_defaults(func=_design)
        elif name == "verify":
            p.add_argument("target", nargs="?", help="model file (.json, .yaml), a problem directory of the "
                                                     "verification library, or a library directory")
            p.add_argument("--all", action="store_true", help="verify every problem of verification/problems")
            p.add_argument("--design", action="store_true",
                           help="the reference solution includes the design results (tier F applies)")
            p.add_argument("--tiers", default="all",
                           help="tiers to run, e.g. ACD (default: all; a tier left out makes 'verified' impossible)")
            p.add_argument("--budget", type=float, default=None,
                           help="time budget in seconds (per model / per problem); checks not started in time "
                                "count as not passed")
            p.add_argument("--extra-evidence", action="store_true",
                           help="also run the independent tableau solver of tests/independent (extra evidence only)")
            p.add_argument("--json", action="store_true", help="print the record as JSON")
            p.add_argument("--reference", action="store_true",
                           help="with --json: include the reference results document in the record")
            p.add_argument("-o", "--output", help="write the full verification record (with the reference "
                                                  "solution) to this file")
            p.set_defaults(func=_verify)
        elif name == "generate":
            p.add_argument("family", nargs="?", help="model family (see --list)")
            p.add_argument("--list", action="store_true", help="list the families, their levels and feature flags")
            p.add_argument("--seed", type=int, default=0, help="seed (default 0)")
            p.add_argument("--difficulty", type=int, default=1, help="difficulty level 1-5 (default 1)")
            p.add_argument("--units", choices=["kN-m", "kN-mm", "kip-in"], help="unit system (default: from the seed)")
            p.add_argument("--vertical-axis", choices=["Y", "Z"], help="vertical axis (default: from the seed)")
            p.add_argument("--feature", action="append", metavar="NAME=true|false",
                           help="force a feature flag of the family on or off (repeatable)")
            p.add_argument("-o", "--output", help="model file (.json); default: standard output")
            p.add_argument("--metadata", help="also write the metadata (JSON) to this file")
            p.set_defaults(func=_generate)
        elif name == "perturb":
            p.add_argument("file", nargs="?", help="base model file (.json, .yaml, .yml)")
            p.add_argument("name", nargs="?", help="curve ball to apply (see --list)")
            p.add_argument("--list", action="store_true", help="list the curve balls by category with their "
                                                               "effect classes")
            p.add_argument("--applicable", action="store_true",
                           help="try every curve ball on the model and say which apply (analyses the model many "
                                "times)")
            p.add_argument("--seed", type=int, default=0, help="seed choosing the target (default 0); the same "
                                                               "file, name and seed give the same result")
            p.add_argument("--effect-class", choices=["changes-analysis", "design-only", "creates-mechanism",
                                                      "invalid-model"],
                           help="accept only this effect class (for curve balls that can give more than one)")
            p.add_argument("--require-verified", nargs="?", const="auto", choices=["auto", "analysis", "design"],
                           help="accept only a target whose reference solution passes 'benchfea verify' (analysis "
                                "results or mechanism diagnosis; 'design': with the design results, tier F; default "
                                "'auto': design for a design-only effect). Much slower")
            p.add_argument("-o", "--output", help="perturbed model file (.json); default: standard output")
            p.add_argument("--record", help="write the change record (JSON) to this file")
            p.add_argument("--show-record", action="store_true",
                           help="print the change record to standard output (the model then only goes to -o)")
            p.add_argument("--json", action="store_true",
                           help="print JSON: the record, the list (--list), the table (--applicable) or the refusal")
            p.set_defaults(func=_perturb)
        elif name == "score":
            p.add_argument("answer", help="answer file (JSON; the format is stated in the task)")
            p.add_argument("--task", required=True, help="task file (<id>.task.json)")
            p.add_argument("--answers", help="hidden answers file of the task; default: the dataset layout "
                                             "(answers/<split>/<id>.answers.json next to tasks/)")
            p.add_argument("--budget", type=float, default=300.0,
                           help="seconds for the verification pipeline when a repair / design is re-solved")
            p.add_argument("--json", action="store_true", help="print the scoring record as JSON")
            p.add_argument("-o", "--output", help="write the scoring record (JSON) to this file")
            p.set_defaults(func=_score)
        elif name == "bench":
            p.set_defaults(func=_bench, bench_command=None)
            bsub = p.add_subparsers(dest="bench_command", metavar="sub-command")
            b = bsub.add_parser("make-dataset", help="build verified tasks on disk (docs/BENCHMARK.md section 8)",
                                description="Build the tasks of every (family, curve ball) pair for the given seeds: "
                                            "only tasks with a verified reference are written; everything else "
                                            "goes to the quarantine log of the manifest.")
            b.add_argument("--out", default="data", help="root directory (tasks/ and answers/ are written under it)")
            b.add_argument("--seeds", default="0,1,4", help="instance seeds, comma separated (seed %% 5 == 4: held out)")
            b.add_argument("--families", help="only these families, comma separated")
            b.add_argument("--curve-balls", help="only these curve balls, comma separated")
            b.add_argument("--workers", type=int, default=1, help="parallel processes")
            b.add_argument("--cache", help="directory for the cache of generated and verified cases (makes the "
                                           "build resumable)")
            b.add_argument("--clean", action="store_true", help="remove task files this run did not write")
            b.add_argument("--quiet", action="store_true", help="no progress lines")
            b = bsub.add_parser("run-agent", help="run an agent in the Gymnasium environment",
                                description="Run episodes of BenchFEA-v0 with the agent built by make_agent(env) of "
                                            "a Python file (see examples/agents).")
            b.add_argument("agent", help="Python file, optionally FILE.py:factory (default factory: make_agent)")
            b.add_argument("--episodes", type=int, default=1)
            b.add_argument("--seed", type=int, default=0, help="seed of the first episode (then seed + 1, ...)")
            b.add_argument("--split", default="train", choices=["train", "held-out"])
            b.add_argument("--family")
            b.add_argument("--curve-ball")
            b.add_argument("--difficulty", type=int)
            b.add_argument("--task-type", choices=["predict", "diagnose", "repair", "design"])
            b.add_argument("--mode", choices=["tool-allowed", "no-solver"])
            b.add_argument("--max-steps", type=int, default=30)
            b.add_argument("--cache", help="directory for the cache of generated and verified cases")
            b.add_argument("--dataset", help="serve the stored tasks of this dataset directory (e.g. data) instead of "
                                             "generating and verifying instances at each reset")
            b.add_argument("--json", action="store_true")
        else:
            p.add_argument("args", nargs=argparse.REMAINDER, help="command arguments")
            p.set_defaults(func=_not_implemented)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
