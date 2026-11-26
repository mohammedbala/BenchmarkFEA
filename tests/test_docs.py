"""README.md, CONTRIBUTING.md and the index of the documents: links resolve, the numbers quoted from
generated files are the ones in those files, and the commands shown are commands the CLI accepts.

The commands of README.md were also run by hand in a fresh virtual environment; this test keeps the
cheap part of that true (the command line parses, the files it names exist).
"""
import re
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

from benchfea.cli import build_parser

ROOT = Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text(encoding="utf-8")


def _block(text, name):
    m = re.search(rf"<!-- from:{name} -->\n(.*?)\n<!-- /from:{name} -->", text, re.S)
    assert m, name
    return m.group(1)


def test_links_of_the_documents_resolve():
    p = subprocess.run([sys.executable, str(ROOT / "tools" / "check_links.py")], capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr


def test_readme_sections_are_in_the_stated_order():
    headings = re.findall(r"^## (.+)$", README, re.M)
    assert headings == [
        "1. What it is", "2. What it does not do", "3. Install", "4. Run one example and read the result",
        "5. Apply one curve ball and see what changed", "6. Score one agent answer",
        "7. How verification works, and what \"verified\" means",
        "8. Adding a model, a curve ball or a verification case", "What you must supply yourself",
        "The Gymnasium environment", "The dataset", "Licence", "How to cite", "Acknowledgement"]


def test_readme_quotes_the_verification_report_and_the_data_card_verbatim():
    report = (ROOT / "docs" / "VERIFICATION_REPORT.md").read_text(encoding="utf-8")
    for line in _block(README, "VERIFICATION_REPORT").splitlines():
        assert line.startswith("* ") and line in report.splitlines(), line
    card = (ROOT / "data" / "DATA_CARD.md").read_text(encoding="utf-8")
    headline = re.search(r"<!-- generated:headline -->\n(.*?)\n<!-- /generated:headline -->", card, re.S).group(1)
    assert _block(README, "DATA_CARD") == headline


def test_readme_says_what_verified_means_and_claims_nothing_more():
    assert "passed the tiers listed in its record, at the\ntolerances stated there" in README
    assert "Nobody has checked, approved or certified any result" in README
    assert "![" not in README.split("## 1. What it is")[0]            # no badges


def _commands(text):
    for body in re.findall(r"```bash\n(.*?)```", text, re.S):
        for line in body.splitlines():
            if line.startswith("benchfea "):
                yield line


@pytest.mark.parametrize("name", ["README.md", "CONTRIBUTING.md"])
def test_commands_shown_are_accepted_by_the_command_line(name):
    text = (ROOT / name).read_text(encoding="utf-8")
    commands = list(_commands(text))
    assert commands
    parser = build_parser()
    for line in commands:
        argv = shlex.split(line)[1:]
        if argv == ["--version"]:
            continue
        args = parser.parse_args(argv)                     # exits on a wrong option
        for value in vars(args).values():
            if isinstance(value, str) and value.startswith(("examples/", "data/", "verification/")):
                assert (ROOT / value).exists(), (line, value)


def test_every_document_is_in_the_index():
    index = (ROOT / "docs" / "README.md").read_text(encoding="utf-8")
    for path in sorted((ROOT / "docs").glob("*.md")):
        if path.name != "README.md":
            assert f"]({path.name})" in index, path.name
    for rel in ("../README.md", "../LIMITS.md", "../CONTRIBUTING.md", "../data/DATA_CARD.md",
                "../examples/README.md"):
        assert f"]({rel})" in index, rel


def test_every_command_prints_help_with_an_example_that_names_existing_files():
    parser = build_parser()
    sub = next(a for a in parser._actions if hasattr(a, "choices") and a.choices)
    parsers = dict(sub.choices)
    bench = next(a for a in parsers["bench"]._actions if hasattr(a, "choices") and a.choices)
    parsers.update({f"bench {k}": v for k, v in bench.choices.items()})
    assert {"validate", "run", "design", "perturb", "verify", "score", "generate", "bench", "bench make-dataset",
            "bench run-agent"} <= set(parsers)
    for name, p in parsers.items():
        text = p.format_help()
        assert "example:" in text and f"benchfea {name}" in text, name
        for token in re.findall(r"(?:examples|verification|data)/[\w./<>-]+", text.split("example:", 1)[1]):
            if "<" not in token:
                assert (ROOT / token).exists(), (name, token)
