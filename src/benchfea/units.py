"""Unit systems and exact unit conversion.

Internal system (everything inside the engine, the Model and the Results):

    force   newton (N)
    length  millimetre (mm)
    angle   radian
    => stress N/mm^2 (MPa), moment N*mm, line load N/mm, unit weight N/mm^3

A file declares ``units: {force, length}``; every dimensional number in the file
is in the *consistent* unit built from those two (e.g. kip & in -> stress in
kip/in^2 = ksi, moment in kip*in, unit weight in kip/in^3).  Numbers are
converted, never relabelled.

Exact conversion factors (kept as rational numbers):

    1 in  = 25.4 mm   exactly      (NIST SP 811, Appendix B.8, "inch", exact)
    1 ft  = 304.8 mm  exactly      (NIST SP 811, Appendix B.8, "foot", exact)
    1 lbf = 0.45359237 kg * 9.80665 m/s^2 = 4.4482216152605 N exactly
            (NIST SP 811, Appendix B.8 footnotes 22 and 23)
    1 kip = 1000 lbf

Dimension strings (the ``x-unit`` annotation of the schema) are products and
quotients of ``force``, ``length``, ``moment`` (= force*length) and
``rotation`` / ``dimensionless`` / ``degree`` (no scaling), with integer
powers: e.g. ``"force/length^2"``, ``"length^4"``, ``"moment/rotation"``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction

INTERNAL_FORCE = "N"
INTERNAL_LENGTH = "mm"

_LBF = Fraction(45359237, 100000000) * Fraction(980665, 100000)  # N, exact

FORCE_UNITS: dict[str, Fraction] = {
    "N": Fraction(1),
    "kN": Fraction(1000),
    "lbf": _LBF,
    "kip": 1000 * _LBF,
}
LENGTH_UNITS: dict[str, Fraction] = {
    "mm": Fraction(1),
    "m": Fraction(1000),
    "in": Fraction(254, 10),
    "ft": Fraction(3048, 10),
}

#: x-unit value used for positions along a member: a fraction of the length
#: when the sibling ``positionType`` is "relative" (default), a length when it
#: is "absolute".
POSITION = "position"

_TERM = re.compile(r"^(force|length|moment|rotation|dimensionless|degree)(?:\^(\d+))?$")


def parse_dimension(dim: str) -> tuple[int, int]:
    """Return (force exponent, length exponent) of a dimension string."""
    f = l = 0
    sign = 1
    for tok in re.split(r"([*/])", dim.replace(" ", "")):
        if tok == "*":
            sign = 1
            continue
        if tok == "/":
            sign = -1
            continue
        m = _TERM.match(tok)
        if not m:
            raise ValueError(f"unknown dimension term {tok!r} in {dim!r}")
        name, power = m.group(1), int(m.group(2) or 1)
        if name == "force":
            f += sign * power
        elif name == "length":
            l += sign * power
        elif name == "moment":
            f += sign * power
            l += sign * power
        # rotation, dimensionless, degree: no scaling
    return f, l


@dataclass(frozen=True)
class UnitSystem:
    """A consistent unit system defined by a force unit and a length unit."""

    force: str
    length: str

    def __post_init__(self):
        if self.force not in FORCE_UNITS:
            raise ValueError(f"unsupported force unit {self.force!r}; use one of {sorted(FORCE_UNITS)}")
        if self.length not in LENGTH_UNITS:
            raise ValueError(f"unsupported length unit {self.length!r}; use one of {sorted(LENGTH_UNITS)}")

    @classmethod
    def from_dict(cls, d: dict) -> "UnitSystem":
        return cls(d["force"], d["length"])

    def to_dict(self) -> dict:
        return {"force": self.force, "length": self.length}

    def factor(self, dim: str) -> Fraction:
        """Exact factor: internal value = file value * factor."""
        f, l = parse_dimension(dim)
        return FORCE_UNITS[self.force] ** f * LENGTH_UNITS[self.length] ** l

    def to_internal(self, value, dim: str):
        """Convert a number (or numpy array) from this system to N, mm."""
        fac = self.factor(dim)
        return value if fac == 1 else value * float(fac)

    def from_internal(self, value, dim: str):
        """Convert a number (or numpy array) from N, mm to this system."""
        fac = self.factor(dim)
        return value if fac == 1 else value / float(fac)

    def label(self, dim: str) -> str:
        """Human-readable unit label for a dimension, e.g. 'kip/in^2'."""
        if dim == "degree":
            return "deg"
        if dim == "n/a":
            return ""
        if dim == POSITION:
            return f"- or {self.length}"
        f, l = parse_dimension(dim)
        rot = -1 if re.search(r"/\s*rotation", dim) else (1 if "rotation" in dim else 0)

        def term(sym, p):
            return sym if abs(p) == 1 else f"{sym}^{abs(p)}"

        num = [term(s, p) for s, p in ((self.force, f), (self.length, l), ("rad", rot)) if p > 0]
        den = [term(s, p) for s, p in ((self.force, f), (self.length, l), ("rad", rot)) if p < 0]
        if not num and not den:
            return "-"
        return ("*".join(num) or "1") + ("/" + "/".join(den) if den else "")


INTERNAL = UnitSystem(INTERNAL_FORCE, INTERNAL_LENGTH)


def convert(value, dim: str, src: UnitSystem, dst: UnitSystem):
    """Convert ``value`` of dimension ``dim`` from ``src`` to ``dst`` (exact
    rational factor, applied in one floating-point multiplication)."""
    fac = src.factor(dim) / dst.factor(dim)
    return value if fac == 1 else value * float(fac)
