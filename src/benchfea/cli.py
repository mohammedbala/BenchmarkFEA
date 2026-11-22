"""Command line interface: ``benchfea <command> ...``."""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__

_COMMANDS = {
    "validate": ("Check an input file against the schema and the semantic rules; list every problem.",
                 "benchfea validate examples/portal.json\n  benchfea validate model.yaml --json"),
    "run": ("Analyse a model and write the results document.",
            "benchfea run examples/portal.json -o portal.results.json"),
    "design": ("Analyse a model and run the AISC 360-22 member checks.",
               "benchfea design examples/portal.json -o portal.design.json"),
    "perturb": ("Create a perturbed variant of a model for the benchmark.",
                "benchfea perturb examples/portal.json --seed 1 -o variant.json"),
    "verify": ("Run the verification suite against reference solutions.",
               "benchfea verify"),
    "score": ("Score an agent's submission against the reference results.",
              "benchfea score submission.json --task tasks/portal"),
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
    flagged = [cid for cid, c in results.combinations.items() if not c.invariants.passed]
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
        for w in results.warnings:
            print(f"  warning: {w}")
        print(f"  {len(results.defaults_used)} default(s) used (listed in the output under defaultsUsed)")
        print(f"results written to {out_path}" + (" -- FLAGGED: invariants failed" if flagged else ""))
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
        else:
            p.add_argument("args", nargs=argparse.REMAINDER, help="command arguments")
            p.set_defaults(func=_not_implemented)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
