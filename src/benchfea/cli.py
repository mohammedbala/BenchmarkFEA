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
        else:
            p.add_argument("args", nargs=argparse.REMAINDER, help="command arguments")
            p.set_defaults(func=_not_implemented)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
