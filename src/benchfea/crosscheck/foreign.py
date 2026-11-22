"""Results of a foreign solver mapped into benchfea conventions, and the base class of a translation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .fidelity import Fidelity

__all__ = ["ForeignResult", "Translation", "TranslationError"]


class TranslationError(ValueError):
    """The model cannot be turned into a foreign model at all (as opposed to a feature that is
    recorded as ``unsupported`` in the fidelity record)."""


@dataclass
class ForeignResult:
    """One combination solved by a foreign solver, in benchfea conventions and internal units
    (docs/CONVENTIONS.md): ``displacements`` and ``reactions`` (n_nodes, 6) global, node order of the
    model, reactions = forces BY the supports ON the structure; ``end_forces[member]`` (12,) local, on
    the member, at the ends of its flexible length.  Results of a second-order combination are already
    divided by ``alpha``.  ``ok = False`` means the solver failed; ``error`` says why."""
    solver: str
    combination: str
    ok: bool
    fidelity: Fidelity
    displacements: np.ndarray | None = None
    reactions: np.ndarray | None = None
    end_forces: dict = field(default_factory=dict)
    error: str = ""
    info: dict = field(default_factory=dict)       # solver diagnostics (iterations, residuals, parameters, ...)

    def to_dict(self) -> dict:
        return {"solver": self.solver, "combination": self.combination, "ok": self.ok, "error": self.error,
                "fidelity": self.fidelity.to_dict(), "info": self.info,
                "displacements": None if self.displacements is None else self.displacements.tolist(),
                "reactions": None if self.reactions is None else self.reactions.tolist(),
                "endForces": {k: v.tolist() for k, v in self.end_forces.items()}}


@dataclass
class Translation:
    """A foreign model ready to run: ``job`` is plain JSON data (the list of library calls and the
    outputs to read) executed by ``benchfea.crosscheck.runner`` in a subprocess; ``controls`` names the
    refinement parameters of approximate features: ``{name: {"coarse": value, "fine": value,
    "order": p}}`` (error assumed proportional to value**-p); ``results(raw)`` maps the raw output of
    the worker to a :class:`ForeignResult`."""
    solver: str
    combination: str
    job: dict
    fidelity: Fidelity
    controls: dict = field(default_factory=dict)
    options: dict = field(default_factory=dict)
    mapper: object = None

    def results(self, raw: dict) -> ForeignResult:
        if not raw.get("ok"):
            return ForeignResult(self.solver, self.combination, False, self.fidelity,
                                 error=raw.get("error", "unknown failure"), info=raw.get("info", {}))
        return self.mapper(raw)
