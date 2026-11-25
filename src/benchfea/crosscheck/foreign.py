"""Results of a foreign solver mapped into benchfea conventions, and the base class of a translation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .fidelity import Fidelity

__all__ = ["ForeignResult", "Translation", "TranslationError", "refinement_levels"]


def refinement_levels(fine: int, orders, default=None) -> dict:
    """The ``controls`` entry of a refinement parameter (elements per member): the levels the runner
    solves, coarsest first, the last being the delivered run.  Three levels with the constant ratio 2
    (``default``, or fine/4, fine/2, fine when ``fine`` is a multiple of 4); two levels when ``fine``
    is only even (no order can be observed: the bound then assumes the lowest order of ``orders``);
    one level otherwise (no estimate: the solver is not eligible)."""
    fine = int(fine)
    if default is not None and int(default[-1]) == fine:
        levels = [int(v) for v in default]
    elif fine % 4 == 0:
        levels = [fine // 4, fine // 2, fine]
    elif fine % 2 == 0:
        levels = [fine // 2, fine]
    else:
        levels = [fine]
    return {"levels": levels, "fine": fine, "orders": [float(orders[0]), float(orders[1])]}


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
    refinement parameters of approximate features.  A discretisation parameter:
    ``{name: {"levels": [n/4, n/2, n], "fine": n, "orders": [p_min, p_max]}}`` (:func:`refinement_levels`;
    the error bound of ``runner.estimate_from_levels``).  A penalty / emulation stiffness:
    ``{name: {"coarse": value, "fine": value, "order": None}}`` (round-off competes: twice the
    difference).  ``results(raw)`` maps the raw output of the worker to a :class:`ForeignResult`."""
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
