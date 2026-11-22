"""Fidelity record of a translation: how faithfully a foreign solver represents each feature of a
benchfea model.

Levels
------
``exact``        the foreign model is mathematically the same problem (differences = round-off);
``approximate``  the foreign model converges to the same problem; ``error`` is the ESTIMATED relative
                 error of the delivered results (``None`` until the runner has measured it from two
                 parameter levels, or when no estimate exists);
``unsupported``  the foreign solver cannot represent the feature, or the translator refuses to guess;
                 results (if any were still produced) must not be used as a cross-check.
"""
from __future__ import annotations

from dataclasses import dataclass, field

EXACT, APPROXIMATE, UNSUPPORTED = "exact", "approximate", "unsupported"
_ORDER = {EXACT: 0, APPROXIMATE: 1, UNSUPPORTED: 2}

__all__ = ["Fidelity", "FidelityItem", "EXACT", "APPROXIMATE", "UNSUPPORTED"]


@dataclass
class FidelityItem:
    feature: str
    level: str
    note: str = ""
    error: float | None = None        # estimated relative error (approximate items)
    control: str | None = None        # name of the refinement parameter the estimate comes from

    def to_dict(self) -> dict:
        return {"feature": self.feature, "level": self.level, "note": self.note, "error": self.error,
                "control": self.control}


@dataclass
class Fidelity:
    solver: str
    items: dict[str, FidelityItem] = field(default_factory=dict)

    def add(self, feature: str, level: str, note: str = "", error: float | None = None,
            control: str | None = None) -> None:
        """Record a feature; a repeated feature keeps the worse level and collects the notes."""
        old = self.items.get(feature)
        if old is None:
            self.items[feature] = FidelityItem(feature, level, note, error, control)
            return
        if _ORDER[level] > _ORDER[old.level]:
            old.level, old.control = level, control or old.control
        if note and note not in old.note:
            old.note = (old.note + "; " + note) if old.note else note
        if error is not None:
            old.error = max(error, old.error or 0.0)

    @property
    def level(self) -> str:
        return max((i.level for i in self.items.values()), key=_ORDER.get, default=EXACT)

    @property
    def unsupported(self) -> list[FidelityItem]:
        return [i for i in self.items.values() if i.level == UNSUPPORTED]

    @property
    def approximate(self) -> list[FidelityItem]:
        return [i for i in self.items.values() if i.level == APPROXIMATE]

    @property
    def estimated_error(self) -> float | None:
        """Sum of the estimated errors of the approximate items; ``None`` if one has no estimate."""
        total = 0.0
        for i in self.approximate:
            if i.error is None:
                return None
            total += i.error
        return total

    def set_error(self, control: str, error: float) -> None:
        for i in self.items.values():
            if i.level == APPROXIMATE and i.control == control:
                i.error = error

    def to_dict(self) -> dict:
        return {"solver": self.solver, "level": self.level, "estimatedError": self.estimated_error,
                "items": [i.to_dict() for i in self.items.values()]}

    @classmethod
    def from_dict(cls, d: dict) -> "Fidelity":
        f = cls(d["solver"])
        for i in d["items"]:
            f.items[i["feature"]] = FidelityItem(i["feature"], i["level"], i["note"], i["error"], i["control"])
        return f
