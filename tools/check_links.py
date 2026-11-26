#!/usr/bin/env python3
"""Check that the relative links of the Markdown documents resolve.

    python3 tools/check_links.py            # README.md, CONTRIBUTING.md, LIMITS.md, docs/, data/, examples/,
                                            # verification/problems/*/README.md
    python3 tools/check_links.py FILE.md …  # only these files

For every Markdown link or image ``[text](target)`` outside code blocks and code spans:

* a relative target must be an existing file or directory of the repository;
* a ``#fragment`` on a Markdown target (or on the file itself) must be the anchor of a heading of that
  file, formed the way the repository host forms heading anchors (lower case, punctuation dropped,
  spaces to hyphens, ``-1``, ``-2`` … for repeated headings).

Targets with a scheme (``https:``, ``mailto:`` …) are not fetched: this tool makes no network request.
Exit code 1 when a link does not resolve; every one is listed with its file and line.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ["README.md", "CONTRIBUTING.md", "LIMITS.md", "docs/*.md", "data/*.md", "data/user/*.md",
           "examples/*.md", "examples/agents/*.md", "verification/problems/*/README.md"]
LINK = re.compile(r"!?\[(?:[^\]\[]|\[[^\]]*\])*\]\(\s*<?([^)\s>]+)>?(?:\s+\"[^\"]*\")?\s*\)")
FENCE = re.compile(r"^\s*(```|~~~)")
SPAN = re.compile(r"(`+)(?:(?!\1).)+\1")
HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")


def strip_code(lines):
    """(line number, text without code spans) of every line outside fenced code blocks."""
    fenced = False
    for n, line in enumerate(lines, 1):
        if FENCE.match(line):
            fenced = not fenced
            continue
        if not fenced:
            yield n, line


def slug(heading: str) -> str:
    text = re.sub(r"!?\[((?:[^\]\[]|\[[^\]]*\])*)\]\([^)]*\)", r"\1", heading)     # links: their text
    text = text.replace("`", "").replace("*", "").strip().lower()
    text = re.sub(r"<[^>]+>", "", text)
    text = "".join(c for c in text if c.isalnum() or c in " -_")
    return text.replace(" ", "-")


_anchors: dict = {}


def anchors(path: Path) -> set:
    if path not in _anchors:
        seen, out = {}, set()
        for _, line in strip_code(path.read_text(encoding="utf-8").splitlines()):
            m = HEADING.match(line)
            if m:
                s = slug(m.group(2))
                k = seen.get(s, 0)
                seen[s] = k + 1
                out.add(s if k == 0 else f"{s}-{k}")
        _anchors[path] = out
    return _anchors[path]


def check(path: Path) -> list:
    problems = []
    for n, line in strip_code(path.read_text(encoding="utf-8").splitlines()):
        for m in LINK.finditer(SPAN.sub("", line)):
            target = m.group(1)
            if re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:", target):
                continue                                           # a scheme: not fetched
            ref, _, fragment = target.partition("#")
            dest = path if not ref else (path.parent / ref)
            try:
                dest = dest.resolve()
                dest.relative_to(ROOT)
            except ValueError:
                problems.append((n, target, "points outside the repository"))
                continue
            if not dest.exists():
                problems.append((n, target, "no such file"))
                continue
            if fragment and dest.suffix.lower() == ".md" and fragment.lower() not in anchors(dest):
                problems.append((n, target, f"no heading with the anchor #{fragment}"))
    return problems


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    files = [Path(a).resolve() for a in args] if args else sorted({p for g in DEFAULT for p in ROOT.glob(g)})
    bad, links = 0, 0
    for path in files:
        problems = check(path)
        for n, target, why in problems:
            print(f"{path.relative_to(ROOT)}:{n}: {target}: {why}")
        bad += len(problems)
    print(f"{len(files)} file(s) checked: {bad} link(s) do not resolve")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
