#!/usr/bin/env python3
"""Apply the held-out selection rule of docs/BENCHMARK.md section 7 to the applicability table.

    python3 tools/apply_split_rule.py            # rewrite src/benchfea/bench/splits.json
    python3 tools/apply_split_rule.py --check    # exit 1 if the file would change

Reads the table of docs/CURVEBALLS.md section 5 (written by ``tools/build_curveball_docs.py`` from
actual runs), takes as candidates of a family the curve balls whose cell is neither ``-`` nor ``I``
alone, and writes ``heldOutSelection.applicable`` and ``heldOutCombinations``
(``benchfea.bench.splits.select_held_out``) into ``splits.json``.  Run it after the curve-ball
document was rebuilt.  It prints whether the held-out combinations changed: if they did, the split of
existing instances changed and the dataset has to be rebuilt.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def applicable() -> dict:
    text = (ROOT / "docs" / "CURVEBALLS.md").read_text(encoding="utf-8")
    section = text.split("## 5. Applicability", 1)[1].split("\n## ", 1)[0]
    rows = [[c.strip(" `") for c in line.strip().strip("|").split("|")]
            for line in section.splitlines() if line.startswith("|")]
    families = rows[0][1:]
    out = {f: [] for f in families}
    for row in rows[2:]:
        for f, cell in zip(families, row[1:]):
            if cell != "-" and set(cell.split()) - {"I"}:
                out[f].append(row[0])
    return out


def main(argv) -> int:
    from benchfea.bench.splits import SPLITS_FILE, select_held_out
    old_text = SPLITS_FILE.read_text(encoding="utf-8")
    d = json.loads(old_text)
    before = d["heldOutCombinations"]
    d["heldOutSelection"]["applicable"] = applicable()
    d["heldOutCombinations"] = select_held_out(d["heldOutSelection"]["applicable"], d["heldOutSelection"]["perFamily"])
    new_text = json.dumps(d, indent=1) + "\n"
    changed = d["heldOutCombinations"] != before
    print("held-out combinations " + ("CHANGED: rebuild the dataset" if changed else "unchanged"))
    if new_text == old_text:
        print("up to date")
        return 0
    if "--check" in argv:
        print("out of date (python3 tools/apply_split_rule.py): src/benchfea/bench/splits.json")
        return 1
    SPLITS_FILE.write_text(new_text, encoding="utf-8")
    print("rewritten: src/benchfea/bench/splits.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
