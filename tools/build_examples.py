#!/usr/bin/env python3
"""Rewrite the generated parts of ``examples/`` from what the engine computes for the example files.

    python3 tools/build_examples.py            # rewrite
    python3 tools/build_examples.py --check    # exit 1 if a file would change

Generated, and never edited by hand:

* ``examples/expected.json`` -- the key results of every example (12 significant digits);
* in every ``examples/<name>.md``: the table under ``## Key results`` and the table under
  ``## Design summary``;
* ``examples/01_simple_beam.yaml`` -- the YAML twin of the first example.

The quantities and their labels are defined once, in ``tests/test_examples.py`` (``QUANTITIES``); that
test recomputes every row and compares it with ``expected.json``.  The numbers are what the engine
gives for these files; they are not reference values of any publication.  The prose of the ``.md``
files is not generated: when a number changes, read the text around the tables again.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
EXAMPLES = ROOT / "examples"
YAML_HEAD = "# Example 01 in YAML: exactly the same model as 01_simple_beam.json (JSON is the canonical form).\n"


def replace_table(text: str, heading: str, table: str, where: str) -> str:
    """``text`` with the first markdown table after ``heading`` replaced by ``table``."""
    lines = text.split("\n")
    try:
        start = lines.index(heading)
    except ValueError:
        raise SystemExit(f"{where}: heading {heading!r} not found")
    first = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("|")), None)
    nxt = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
    if first is None or first > nxt:
        raise SystemExit(f"{where}: no table under {heading!r}")
    last = first
    while last + 1 < len(lines) and lines[last + 1].startswith("|"):
        last += 1
    return "\n".join(lines[:first] + table.split("\n") + lines[last + 1:])


def outputs() -> dict:
    """path -> text of every generated file."""
    import test_examples as T
    from benchfea.io import dumps, load

    out, expected = {}, {}
    for name in T.NAMES:
        run = T.Run(name)
        quantities = T.QUANTITIES[name](run)
        expected[name] = {k: {"value": float(f"{v:.12g}"), "unit": u} for k, (v, u) in quantities.items()}
        key = ["| quantity | value | unit |", "|---|---|---|"]
        key += [f"| {k} | {T.format_value(expected[name][k]['value'])} | {u} |" for k, (_, u) in quantities.items()]
        design = ["| design member | status | ratio | combination | governing limit state |", "|---|---|---|---|---|"]
        for k, m in run.design.members.items():
            g = m.governing
            design.append(f"| {k} | {m.status} | {g.ratio:.3f} | {g.combination} | {g.name} |")
        path = EXAMPLES / f"{name}.md"
        text = path.read_text(encoding="utf-8")
        text = replace_table(text, "## Key results", "\n".join(key), path.name)
        text = replace_table(text, "## Design summary", "\n".join(design), path.name)
        out[path] = text
    out[EXAMPLES / "expected.json"] = json.dumps(expected, indent=2) + "\n"
    out[EXAMPLES / "01_simple_beam.yaml"] = YAML_HEAD + dumps(load(EXAMPLES / "01_simple_beam.json").data, "yaml")
    return out


def main(argv) -> int:
    check = "--check" in argv
    changed = []
    for path, text in outputs().items():
        if not path.exists() or path.read_text(encoding="utf-8") != text:
            changed.append(path.relative_to(ROOT).as_posix())
            if not check:
                path.write_text(text, encoding="utf-8")
    if check and changed:
        print("out of date (python3 tools/build_examples.py):", ", ".join(changed))
        return 1
    print("up to date" if not changed else "rewritten: " + ", ".join(changed))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
