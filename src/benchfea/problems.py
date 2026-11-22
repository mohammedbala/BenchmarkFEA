"""Problem records and the exceptions that carry them.

Every refusal of an input file is expressed as a list of :class:`Problem`
objects: *all* problems found, each with a location and a hint.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any

ERROR = "error"
WARNING = "warning"


@dataclass(frozen=True)
class Problem:
    """One problem found in an input document.

    path      JSON pointer (RFC 6901) to the offending value ("" = document root).
    line,col  1-based position in the source text, or None when the document
              was supplied as a dict (or the position is unknown).
    message   What is wrong.
    hint      How to fix it.
    severity  "error" (model refused) or "warning" (model accepted, reported).
    code      Stable machine-readable identifier, e.g. "ref.unknown-node".
    """

    path: str
    message: str
    hint: str = ""
    severity: str = ERROR
    code: str = ""
    line: int | None = None
    col: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def format(self, source: str | None = None) -> str:
        loc = source or "<document>"
        if self.line is not None:
            loc += f":{self.line}:{self.col}"
        text = f"{loc}: {self.severity} [{self.code}] at {self.path or '/'}: {self.message}"
        if self.hint:
            text += f"\n    hint: {self.hint}"
        return text


class ModelInvalid(Exception):
    """Raised when an input document is refused. Carries every problem found."""

    def __init__(self, problems: list[Problem], source: str | None = None):
        self.problems = list(problems)
        self.source = source
        errors = [p for p in self.problems if p.severity == ERROR]
        head = f"{len(errors)} error(s) in {source or 'model'}"
        super().__init__(head + "\n" + "\n".join(p.format(source) for p in self.problems))

    @property
    def errors(self) -> list[Problem]:
        return [p for p in self.problems if p.severity == ERROR]


class ModelUnstable(ModelInvalid):
    """Hook for the solver: raised after assembly when the structure (or part of
    it) is a mechanism / the stiffness matrix is singular. The solver must build
    Problems with code "stability.mechanism" whose path points at the node
    (``/nodes/<id>``) and whose message names the unrestrained DOF(s)."""


def pointer(*parts: Any) -> str:
    """Build a JSON pointer from path parts (escaping ~ and /)."""
    return "".join("/" + str(p).replace("~", "~0").replace("/", "~1") for p in parts)
