"""Shared helpers of the verification tiers: documents, hashing, solving, snapshots, comparison."""
from __future__ import annotations

import copy
import hashlib
import json

import numpy as np

from ..crosscheck.compare import Tolerances, compare, from_benchfea
from ..crosscheck.fidelity import EXACT, Fidelity
from ..crosscheck.foreign import ForeignResult
from ..io import Document, load
from ..model import Model, build_model
from ..problems import ModelInvalid
from ..solve import analyze

__all__ = ["raw_document", "document_hash", "solve", "snapshot", "combine", "compare_uniform", "worst_of",
           "new_id", "first_order_counterpart", "is_linear_model", "TRUSS_TYPES", "ONE_SIDED", "SolveFailure"]

TRUSS_TYPES = ("truss", "tensionOnly", "compressionOnly")
ONE_SIDED = ("tensionOnly", "compressionOnly")


class SolveFailure(Exception):
    """A derived model could not be built or analysed; ``reason`` is one line."""

    def __init__(self, reason: str, error: Exception | None = None):
        super().__init__(reason)
        self.reason, self.error = reason, error


def raw_document(source) -> dict:
    """The input document as a plain dict (deep copy).  Accepts a path, a dict, an ``io.Document`` or a
    ``Model`` that was built from a file (``model.source``)."""
    if isinstance(source, Model):
        if not source.source:
            raise TypeError("verify_model needs the input document: this Model was not built from a file. "
                            "Pass the dict or the path instead (the equivalence checks rewrite the document).")
        source = source.source
    if isinstance(source, Document):
        return copy.deepcopy(source.data)
    if isinstance(source, dict):
        return copy.deepcopy(source)
    return copy.deepcopy(load(source).data)


def document_hash(raw: dict) -> str:
    """SHA-256 of the canonical JSON of the document (key order of objects is part of the model: node
    order = file order, so keys are NOT sorted)."""
    return hashlib.sha256(json.dumps(raw, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


def solve(raw: dict, strict: bool = True):
    """Build and analyse a derived document; raise :class:`SolveFailure` with a one-line reason."""
    try:
        model = build_model(copy.deepcopy(raw))
        return model, analyze(model, strict=strict)
    except ModelInvalid as e:
        msgs = "; ".join(f"[{p.code}] {p.path}: {p.message}" for p in e.problems[:3])
        raise SolveFailure(msgs or str(e), e) from e


def snapshot(model, results, cid: str, n_nodes: int | None = None) -> ForeignResult:
    """A combination result as a comparable record, optionally restricted to the first ``n_nodes``
    nodes (derived models keep the original nodes first)."""
    s = from_benchfea(model, results, cid)
    if n_nodes is not None:
        s.displacements = s.displacements[:n_nodes]
        s.reactions = s.reactions[:n_nodes]
    return s


def combine(parts, name: str = "combined") -> ForeignResult:
    """Linear combination ``sum(c * snapshot)`` of comparable records: ``parts = [(c, snapshot), ...]``."""
    first = parts[0][1]
    fid = Fidelity(name)
    fid.add("derived", EXACT, "linear combination of analyses")
    common = set(first.end_forces)
    for _, s in parts[1:]:
        common &= set(s.end_forces)
    return ForeignResult(name, first.combination, True, fid,
                         sum(c * s.displacements for c, s in parts), sum(c * s.reactions for c, s in parts),
                         {m: sum(c * s.end_forces[m] for c, s in parts) for m in first.end_forces if m in common})


def compare_uniform(a: ForeignResult, b: ForeignResult, model, tol: float):
    return compare(a, b, model, Tolerances(*([tol] * 6)))


def worst_of(comparison) -> dict:
    d = max(comparison.worst.values(), key=lambda v: v.relative)
    return {"quantity": d.quantity, "where": d.where, "a": d.a, "b": d.b, "value": d.relative}


def new_id(existing, base: str) -> str:
    """An id not in ``existing`` (ids: letters, digits, '_', '.', '-'; at most 64 characters)."""
    base = base[:56]
    cand, k = base, 1
    while cand in existing:
        k += 1
        cand = f"{base}-{k}"
    return cand


def first_order_counterpart(raw: dict) -> dict:
    """The same document with every combination analysed first order and without buckling analyses."""
    d = copy.deepcopy(raw)
    an = d.get("analysis")
    if an:
        an.pop("type", None)
        an.pop("buckling", None)
    for c in d.get("combinations", {}).values():
        c.pop("analysisType", None)
    return d


def is_linear_model(model) -> bool:
    """Every combination goes through the plain linear path (superposition holds)."""
    return (not any(m.type in ONE_SIDED for m in model.members.values())
            and all(c.analysis_type == "firstOrder" and not c.stiffness_reduction and c.notional is None
                    for c in model.combinations.values()))


def arrays_identical(a, b) -> bool:
    return np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True)
