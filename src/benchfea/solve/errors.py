"""Refusals of the solver and the structured mechanism diagnosis.

``ModelRefused``   a model the solver will not analyse (unsupported analysis
                   type, member-level mechanism, failed invariant, ...).
``ModelUnstable``  (from :mod:`benchfea.problems`) raised by :func:`unstable`
                   with the attribute ``diagnosis`` (:class:`MechanismDiagnosis`).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..problems import ModelInvalid, ModelUnstable, Problem, pointer

__all__ = ["ModelRefused", "MechanismDiagnosis", "MechanismMode", "NodeParticipation", "refuse", "unstable"]

NOT_IMPLEMENTED = "analysis.not-implemented"


class ModelRefused(ModelInvalid):
    """The solver refuses the model. ``problems`` say why (codes ``analysis.not-implemented``,
    ``member.release-mechanism``, ``section.zero-stiffness``, ``load.torque-on-truss``,
    ``solve.invariant-failed``)."""

    @property
    def codes(self) -> list[str]:
        return [p.code for p in self.problems]


def refuse(code: str, path: str, message: str, hint: str = "") -> ModelRefused:
    return ModelRefused([Problem(path, message, hint, code=code)])


@dataclass(frozen=True)
class NodeParticipation:
    """Share of one node in a mechanism (shares of all nodes of a mode sum to 1) and its split over the
    six global DOFs.  Translations and rotations are made comparable by multiplying rotations by the
    model size (diagonal of the bounding box of the nodes)."""
    node: str
    share: float
    dofs: dict[str, float]

    def to_dict(self) -> dict:
        return {"node": self.node, "share": round(self.share, 6),
                "dofs": {k: round(v, 6) for k, v in self.dofs.items() if round(v, 6) > 0.0}}


@dataclass(frozen=True)
class MechanismMode:
    """One zero-stiffness motion (or one group of them, for ``kind == "distributed"``).

    kind         "rigid-body"  the whole structure (or a disconnected part) moves as a rigid body;
                 "node"        a single node moves or spins on its own;
                 "distributed" the remaining mechanisms, which involve several nodes (``count`` of them).
    motion       "translation" | "rotation" | "translation+rotation" | "mixed"
    direction    unit vector of the translation (global) or None
    axis         unit vector of the rotation axis (global) or None
    through      a point on the rotation axis (rigid-body rotations), mm, or None
    cause        short machine-readable reason: "no-restraint", "negligible-springs" (a rigid-body motion
                 that only support springs oppose, negligible against the structure's stiffness),
                 "all-member-ends-released", "only-truss-members", "member-spin" (a member free to spin
                 about its own axis), "no-stiffness"
    nodes        participating nodes, largest share first
    members      ids of the members connected to the participating nodes
    """
    kind: str
    description: str
    motion: str
    cause: str
    nodes: tuple[NodeParticipation, ...]
    members: tuple[str, ...] = ()
    direction: tuple[float, ...] | None = None
    axis: tuple[float, ...] | None = None
    through: tuple[float, ...] | None = None
    count: int = 1

    def to_dict(self) -> dict:
        def vec(v, nd=6):
            return None if v is None else [round(float(x), nd) + 0.0 for x in v]
        return {"kind": self.kind, "description": self.description, "motion": self.motion, "cause": self.cause,
                "count": self.count, "direction": vec(self.direction), "axis": vec(self.axis),
                "through": vec(self.through, 3), "nodes": [n.to_dict() for n in self.nodes],
                "members": list(self.members)}


@dataclass(frozen=True)
class MechanismDiagnosis:
    """Why the stiffness matrix is singular.  ``to_dict()`` is deterministic for a given model (it
    contains no round-off quantities) and is meant to be compared against as a reference answer;
    ``numerics`` holds the eigenvalue information, which is round-off dependent."""
    n_mechanisms: int                     # dimension of the null space of the reduced stiffness matrix
    modes: tuple[MechanismMode, ...]
    tolerance: float                      # threshold on the eigenvalues of the diagonally scaled matrix
    n_free_dofs: int
    numerics: dict = field(default_factory=dict)
    classification: str = "mechanism"     # "mechanism": zero stiffness to round-off; "ill-conditioned": the
    #                                       smallest scaled eigenvalue is below the tolerance but not zero
    #                                       (stiffness contrast of about 1 / eigenvalue)

    @property
    def summary(self) -> str:
        return "; ".join(m.description for m in self.modes)

    def to_dict(self) -> dict:
        return {"stable": False, "classification": self.classification, "nMechanisms": self.n_mechanisms,
                "nFreeDofs": self.n_free_dofs,
                "tolerance": self.tolerance, "modes": [m.to_dict() for m in self.modes]}


def unstable(diagnosis: MechanismDiagnosis, source: str | None = None) -> ModelUnstable:
    """Build the exception: one Problem per mode, ``.diagnosis`` attached."""
    problems = []
    for m in diagnosis.modes:
        if m.kind == "rigid-body":
            path, hint = "/supports", "Add supports (or members) that restrain this motion."
        else:
            path = pointer("nodes", m.nodes[0].node) if m.nodes else ""
            hint = ("Keep at least one member end connected in this direction, restrain the node, or connect it "
                    "to the rest of the structure.")
        if diagnosis.classification == "mechanism":
            head = "the structure is a mechanism: "
        else:
            head = ("the stiffness matrix is ill-conditioned (stiffness contrast beyond 1 / pivotTolerance), "
                    "almost a mechanism: ")
            hint = ("Avoid members or springs that are many orders of magnitude stiffer than the rest (use rigid "
                    "offsets, rigid links or diaphragms for rigid parts) and check for nearly unrestrained motions. "
                    + hint)
        problems.append(Problem(path, head + m.description, hint, code="stability.mechanism"))
    e = ModelUnstable(problems, source)
    e.diagnosis = diagnosis
    return e
