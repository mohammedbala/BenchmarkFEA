"""The change record of a curve ball and the machinery that makes it machine-checkable.

A :class:`Change` is one edit of the input document addressed by a JSON pointer (RFC 6901):

    op "replace"  the value at ``path`` was ``before`` and is ``after``
    op "add"      ``path`` did not exist and now holds ``after``
    op "remove"   ``path`` held ``before`` and no longer exists

``index`` (add / remove of an object member) is the position of the key among the keys of its
object, so that applying a record and inverting it reproduce the key ORDER too (node order = file
order in this format).  Arrays are never edited element-wise when their length changes: the whole
array is replaced.

    changes = diff(base, perturbed)                     # exact, minimal structural difference
    apply_changes(base, changes)      == perturbed      # forward
    apply_changes(perturbed, invert_changes(changes)) == base

The record is not written by hand by the curve balls: :func:`diff` derives it from the two
documents, so it describes the difference by construction; the tests check both directions.
"""
from __future__ import annotations

import copy
import json
from dataclasses import dataclass, field

__all__ = ["Change", "ChangeRecord", "diff", "apply_changes", "invert_changes", "pointer", "split_pointer",
           "resolve", "plain_copy", "RecordMismatch", "RECORD_FORMAT", "RECORD_VERSION"]

RECORD_FORMAT = "benchfea-curveball-record"
RECORD_VERSION = "1.0.0"
_ABSENT = object()


class RecordMismatch(ValueError):
    """A change does not fit the document it is applied to (path missing, or ``before`` differs)."""


def pointer(*parts) -> str:
    """JSON pointer from path parts (escaping ``~`` and ``/``)."""
    return "".join("/" + str(p).replace("~", "~0").replace("/", "~1") for p in parts)


def split_pointer(path: str) -> list:
    if path == "":
        return []
    if not path.startswith("/"):
        raise ValueError(f"not a JSON pointer: {path!r}")
    return [p.replace("~1", "/").replace("~0", "~") for p in path[1:].split("/")]


def _step(node, key: str, path: str):
    if isinstance(node, dict):
        if key not in node:
            raise RecordMismatch(f"{path}: no member {key!r}")
        return node[key]
    if isinstance(node, list):
        try:
            return node[int(key)]
        except (ValueError, IndexError):
            raise RecordMismatch(f"{path}: no array element {key!r}") from None
    raise RecordMismatch(f"{path}: {key!r} addresses into a value that is neither object nor array")


def resolve(document, path: str):
    """The value at a JSON pointer (raises :class:`RecordMismatch` if it does not exist)."""
    node = document
    for key in split_pointer(path):
        node = _step(node, key, path)
    return node


def plain_copy(value):
    """A copy in which no two places share an object.  ``copy.deepcopy`` keeps aliasing: a document
    built in Python may use one dict for two members (``releases = {"i": r, "j": r}``), and an edit of
    one would then silently edit the other.  Tuples become lists (as in JSON)."""
    if isinstance(value, dict):
        return {k: plain_copy(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain_copy(v) for v in value]
    return value


def _same(a, b) -> bool:
    """Equality of JSON values INCLUDING the key order of objects and the int / float distinction of
    whole numbers as written (``1`` and ``1.0`` are different texts)."""
    return json.dumps(a, sort_keys=False) == json.dumps(b, sort_keys=False)


@dataclass(frozen=True)
class Change:
    op: str                    # "replace" | "add" | "remove"
    path: str                  # JSON pointer into the input document
    before: object = _ABSENT   # absent for "add"
    after: object = _ABSENT    # absent for "remove"
    index: int | None = None   # add / remove of an object member: position of the key in its object

    def to_dict(self) -> dict:
        d = {"op": self.op, "path": self.path}
        if self.before is not _ABSENT:
            d["before"] = copy.deepcopy(self.before)
        if self.after is not _ABSENT:
            d["after"] = copy.deepcopy(self.after)
        if self.index is not None:
            d["index"] = self.index
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "Change":
        return cls(d["op"], d["path"], d.get("before", _ABSENT), d.get("after", _ABSENT), d.get("index"))

    def inverse(self) -> "Change":
        if self.op == "replace":
            return Change("replace", self.path, self.after, self.before)
        if self.op == "add":
            return Change("remove", self.path, before=self.after, index=self.index)
        if self.op == "remove":
            return Change("add", self.path, after=self.before, index=self.index)
        raise ValueError(self.op)


def _as_change(c) -> Change:
    return c if isinstance(c, Change) else Change.from_dict(c)


# ------------------------------------------------------------------------------------------ diff
def diff(before, after, path: str = "") -> list:
    """Structural difference of two JSON values as a list of :class:`Change`, in an order in which
    they can be applied one after the other (removals of a level, then additions, then the changes
    inside the members both documents have)."""
    if isinstance(before, dict) and isinstance(after, dict):
        common = [k for k in before if k in after]
        if common != [k for k in after if k in before]:
            # the shared keys were reordered: not expressible by add / remove, replace the object
            return [Change("replace", path, plain_copy(before), plain_copy(after))]
        out = []
        b_keys, a_keys = list(before), list(after)
        for i in range(len(b_keys) - 1, -1, -1):                     # removals, last key first
            k = b_keys[i]
            if k not in after:
                out.append(Change("remove", path + pointer(k), before=plain_copy(before[k]), index=i))
        for i, k in enumerate(a_keys):                               # additions, first key first
            if k not in before:
                out.append(Change("add", path + pointer(k), after=plain_copy(after[k]), index=i))
        for k in common:
            out.extend(diff(before[k], after[k], path + pointer(k)))
        return out
    if isinstance(before, list) and isinstance(after, list) and len(before) == len(after):
        out = []
        for i, (b, a) in enumerate(zip(before, after)):
            out.extend(diff(b, a, path + pointer(i)))
        return out
    if _same(before, after):
        return []
    return [Change("replace", path, plain_copy(before), plain_copy(after))]


# ----------------------------------------------------------------------------------------- apply
def _insert(obj: dict, key: str, value, index: int | None) -> None:
    if index is None or index >= len(obj):
        obj[key] = value
        return
    items = list(obj.items())
    items.insert(index, (key, value))
    obj.clear()
    obj.update(items)


def apply_changes(document: dict, changes, check: bool = True) -> dict:
    """A new document = ``document`` with ``changes`` applied in order (the input is not modified).

    ``check=True`` (default): every change must fit -- the path of a "replace" / "remove" exists and
    holds exactly ``before``, the path of an "add" does not exist, and a recorded ``index`` is the
    position the key really has; otherwise :class:`RecordMismatch`."""
    doc = plain_copy(document)
    for c in map(_as_change, changes):
        parts = split_pointer(c.path)
        if not parts:
            if c.op != "replace":
                raise RecordMismatch("the document root can only be replaced")
            if check and not _same(doc, c.before):
                raise RecordMismatch("root: 'before' is not the document")
            doc = plain_copy(c.after)
            continue
        parent = doc
        for key in parts[:-1]:
            parent = _step(parent, key, c.path)
        key = parts[-1]
        if isinstance(parent, list):
            if c.op != "replace":
                raise RecordMismatch(f"{c.path}: array elements are replaced, never added or removed")
            current = _step(parent, key, c.path)
            if check and not _same(current, c.before):
                raise RecordMismatch(f"{c.path}: holds {current!r}, the record says {c.before!r}")
            parent[int(key)] = plain_copy(c.after)
        elif isinstance(parent, dict):
            if c.op == "add":
                if key in parent:
                    raise RecordMismatch(f"{c.path}: already exists")
                _insert(parent, key, plain_copy(c.after), c.index)
                if check and c.index is not None and list(parent).index(key) != c.index:
                    raise RecordMismatch(f"{c.path}: cannot be placed at position {c.index}")
            elif c.op in ("replace", "remove"):
                if key not in parent:
                    raise RecordMismatch(f"{c.path}: does not exist")
                if check and not _same(parent[key], c.before):
                    raise RecordMismatch(f"{c.path}: holds {parent[key]!r}, the record says {c.before!r}")
                if c.op == "replace":
                    parent[key] = plain_copy(c.after)
                else:
                    if check and c.index is not None and list(parent).index(key) != c.index:
                        raise RecordMismatch(f"{c.path}: is at position {list(parent).index(key)}, the record "
                                             f"says {c.index}")
                    del parent[key]
            else:
                raise RecordMismatch(f"unknown op {c.op!r}")
        else:
            raise RecordMismatch(f"{c.path}: parent is neither object nor array")
    return doc


def invert_changes(changes) -> list:
    """The changes that undo ``changes`` (reverse order, each one inverted)."""
    return [_as_change(c).inverse() for c in reversed(list(changes))]


# ---------------------------------------------------------------------------------------- record
@dataclass
class ChangeRecord:
    """What a curve ball did to a model -- the reference a diagnosis answer is matched against.

    name, seed          the curve ball and the seed it was applied with
    category            group of the curve ball (constraints, offsets, releases, design, supports, ...)
    changes             exact list of :class:`Change` (JSON-pointer paths with before / after values)
    entities            ids touched, by kind: {"members": [...], "nodes": [...], "constraints": [...], ...}
    target              one line naming what was chosen (e.g. "member B1, end i")
    why                 plain-language reason why the change matters physically
    effect_class        "changes-analysis" | "design-only" | "creates-mechanism" | "invalid-model"
    precondition        what the base model must have for the curve ball to apply
    effect              the MEASURED effect: thresholds, metric values, where the largest change is
    reference           what the reference solution of the perturbed model is: {"kind": "results" |
                        "diagnosis" | "validation-problems", ...} (for a mechanism: the solver's
                        structured diagnosis; for an invalid model: the validator's problems)
    verification        expected verification route of the perturbed model (an expectation, not a result)
    rejected            targets tried before this one and why they were dropped (immaterial, wrong class)
    base_hash, perturbed_hash   SHA-256 of the two documents (``benchfea.verify.common.document_hash``)
    """
    name: str
    seed: int
    category: str
    changes: list
    entities: dict
    target: str
    why: str
    effect_class: str
    precondition: str
    effect: dict = field(default_factory=dict)
    reference: dict = field(default_factory=dict)
    verification: dict = field(default_factory=dict)
    rejected: list = field(default_factory=list)
    base_hash: str = ""
    perturbed_hash: str = ""
    title: str = ""
    version: str = ""

    @property
    def paths(self) -> list:
        return [c.path for c in self.changes]

    def to_dict(self) -> dict:
        return {"format": RECORD_FORMAT, "formatVersion": RECORD_VERSION, "curveBall": self.name,
                "title": self.title, "curveBallVersion": self.version, "seed": self.seed,
                "category": self.category, "effectClass": self.effect_class, "target": self.target,
                "changes": [c.to_dict() for c in self.changes], "entities": copy.deepcopy(self.entities),
                "whyItMatters": self.why, "precondition": self.precondition,
                "effect": copy.deepcopy(self.effect), "reference": copy.deepcopy(self.reference),
                "verification": copy.deepcopy(self.verification), "rejected": copy.deepcopy(self.rejected),
                "baseHash": self.base_hash, "perturbedHash": self.perturbed_hash}

    @classmethod
    def from_dict(cls, d: dict) -> "ChangeRecord":
        return cls(name=d["curveBall"], seed=d["seed"], category=d["category"],
                   changes=[Change.from_dict(c) for c in d["changes"]], entities=d.get("entities", {}),
                   target=d.get("target", ""), why=d.get("whyItMatters", ""), effect_class=d["effectClass"],
                   precondition=d.get("precondition", ""), effect=d.get("effect", {}),
                   reference=d.get("reference", {}), verification=d.get("verification", {}),
                   rejected=d.get("rejected", []), base_hash=d.get("baseHash", ""),
                   perturbed_hash=d.get("perturbedHash", ""), title=d.get("title", ""),
                   version=d.get("curveBallVersion", ""))

    def apply(self, base: dict) -> dict:
        """``base`` with the recorded changes applied (checked against ``before``)."""
        return apply_changes(base, self.changes)

    def revert(self, perturbed: dict) -> dict:
        """The base document recovered from the perturbed one."""
        return apply_changes(perturbed, invert_changes(self.changes))

    def json(self) -> str:
        return json.dumps(self.to_dict(), indent=2) + "\n"
