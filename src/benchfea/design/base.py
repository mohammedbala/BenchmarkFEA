"""Records shared by every design code/edition: the context of one design
member, the demand of one load combination, and the builders of
:class:`benchfea.results.DesignCheck`.

Contract of a limit-state function (every chapter module of an edition)::

    def limit_state(ctx: DesignContext, demand: Demand) -> list[DesignCheck]

* it never raises for a section or input it cannot handle: it returns a check
  with ``status="notChecked"`` (outside the implemented scope, or data missing)
  or ``status="refused"`` (input that contradicts itself), with the reason;
* it never returns an empty list for a limit state it owns (no silent pass);
* everything is in one consistent unit system (the engine calls it with N, mm;
  every equation used is dimensionless-coefficient, so tests may call it with
  kip, in.).

Dimensions use the strings of :mod:`benchfea.units` (``"force"``,
``"length"``, ``"length^2"``, ``"force/length^2"``, ``"dimensionless"``; ``"n/a"``
for text).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any

from ..results import DesignCheck
from .section import DesignSection

LRFD, ASD = "LRFD", "ASD"
OK, NOT_CHECKED, REFUSED = "ok", "notChecked", "refused"

STRESS = "force/length^2"
FORCE = "force"
LENGTH = "length"
AREA = "length^2"
NONE = "dimensionless"
TEXT = "n/a"

#: dimension of the context variables reported as inputs
_DIMS = {"E": STRESS, "G": STRESS, "Fy": STRESS, "Fu": STRESS, "Lx": LENGTH, "Ly": LENGTH, "Lz": LENGTH,
         "L": LENGTH, "An": AREA, "Kx": NONE, "Ky": NONE, "Kz": NONE, "U": NONE, "method": TEXT,
         "maxSlendernessCompression": NONE, "maxSlendernessTension": NONE}


@dataclass(frozen=True)
class Demand:
    """Required strengths of ONE load combination for one design member
    (maxima over the stations of all its analysis members).  Forces are
    magnitudes (>= 0); ``x_*`` is the distance from the design member start.

    ``Lv_*`` = distance from the section of maximum shear to the nearest
    section of zero shear (AISC 360-22 G5), ``None`` when the shear force does
    not vanish anywhere along the member.
    """
    combination: str = ""
    Pt: float = 0.0                   # maximum axial tension
    x_Pt: float | None = None
    Pc: float = 0.0                   # maximum axial compression (positive number)
    x_Pc: float | None = None
    Vy: float = 0.0                   # maximum |shear along local y| (web of an I-shape)
    x_Vy: float | None = None
    Lv_y: float | None = None
    Vz: float = 0.0                   # maximum |shear along local z| (flanges of an I-shape)
    x_Vz: float | None = None
    Lv_z: float | None = None
    V: float = 0.0                    # maximum resultant shear sqrt(Vy^2 + Vz^2) (round sections)
    x_V: float | None = None
    Lv: float | None = None
    extra: dict = field(default_factory=dict)     # room for later limit states (moments, ...)


@dataclass(frozen=True)
class DesignContext:
    """Everything a limit-state function may use for one design member.

    AISC axis letters: x = major axis = bending about local z; y = minor axis =
    bending about local y; z = twist (docs/CONVENTIONS.md section 3).
    ``sources[name]`` is the provenance of each value (reported with every check).
    """
    member: str
    section: DesignSection
    method: str                       # "LRFD" | "ASD"
    E: float | None
    G: float | None
    Fy: float | None
    Fu: float | None
    L: float                          # design member length
    Lx: float | None                  # unbraced lengths (largest segment between brace points)
    Ly: float | None
    Lz: float | None
    Kx: float = 1.0
    Ky: float = 1.0
    Kz: float = 1.0
    An: float | None = None           # net area (D3 / B4.3); None = gross area
    U: float = 1.0                    # shear lag factor (Table D3.1)
    maxSlendernessCompression: float = 200.0
    maxSlendernessTension: float = 300.0
    sources: dict = field(default_factory=dict)
    code: str = "AISC 360-22"

    def source(self, name: str) -> str:
        return self.sources.get(name, "given")

    def inputs(self, *names: str) -> dict:
        """{name: (value, dimension, source)} for the named context variables."""
        return {n: (getattr(self, n), _DIMS.get(n, NONE), self.source(n)) for n in names}

    def with_(self, **changes) -> "DesignContext":
        return replace(self, **changes)

    # ------------------------------------------------------------------ from the model
    @classmethod
    def from_model(cls, model, design_member_id: str) -> "DesignContext":
        p = model.design_members[design_member_id]
        sec = model.sections[p.section]
        src = {n: p.source(n) for n in p.values}
        src["L"] = "derived: design member length"

        def unbraced(name, brace):
            """Unbraced length; with brace points: the longest segment between them
            (the axial force of a combination is checked with the longest segment)."""
            val = getattr(p, name)
            if val is not None:
                return val
            cuts = [0.0, *getattr(p, brace), 1.0]
            seg = max(b - a for a, b in zip(cuts, cuts[1:])) * p.length
            src[name] = f"bracing: longest segment between brace points ({p.values[brace].source})"
            return seg

        return cls(member=design_member_id, section=DesignSection.from_model(sec), method=p.method,
                   E=p.E, G=p.G, Fy=p.Fy, Fu=p.Fu, L=p.length,
                   Lx=unbraced("Lx", "braceMajor"), Ly=unbraced("Ly", "braceMinor"),
                   Lz=unbraced("Lz", "braceTorsion"), Kx=p.Kx, Ky=p.Ky, Kz=p.Kz, An=p.An, U=p.U,
                   maxSlendernessCompression=p.maxSlendernessCompression,
                   maxSlendernessTension=p.maxSlendernessTension, sources=src, code=p.code)


# ============================================================================== check builders
def strength_check(name: str, clause: str, ctx: DesignContext, *, nominal: float, phi: float, omega: float,
                   demand: float, dimension: str = FORCE, combination: str = "", x: float | None = None,
                   inputs: dict | None = None, details: dict | None = None,
                   notes: list[str] | None = None) -> DesignCheck:
    """Available strength and utilisation by AISC 360-22 B3.1 (LRFD, Eq. B3-1:
    Ru <= phi Rn) or B3.2 (ASD, Eq. B3-2: Ra <= Rn / Omega).  Other codes with the
    same format can reuse it; the citation of the design basis is added to ``clause``."""
    if ctx.method == LRFD:
        available, factor, fname, basis = phi * nominal, phi, "phi", "B3.1 Eq. B3-1"
    elif ctx.method == ASD:
        available, factor, fname, basis = nominal / omega, omega, "Omega", "B3.2 Eq. B3-2"
    else:
        return refused(name, clause, f"unknown design method {ctx.method!r} (LRFD or ASD, AISC 360-22 B3)",
                       ctx, demand=demand, dimension=dimension, combination=combination, x=x, inputs=inputs)
    ratio = demand / available if available > 0.0 else math.inf
    inp = {**ctx.inputs("method"), **(inputs or {})}
    return DesignCheck(name=name, clause=f"{clause}; {basis}", demand=demand, capacity=available, ratio=ratio,
                       passed=bool(ratio <= 1.0), combination=combination, x=x, dimension=dimension,
                       details=dict(details or {}), notes=list(notes or []), status=OK, method=ctx.method,
                       nominal=nominal, factor=factor, factor_name=fname, inputs=inp)


def limit_check(name: str, clause: str, ctx: DesignContext, *, value: float, limit: float,
                dimension: str = NONE, advisory: bool = True, inputs: dict | None = None,
                details: dict | None = None, notes: list[str] | None = None) -> DesignCheck:
    """A measured quantity against a limit (e.g. the slenderness recommendations)."""
    ratio = value / limit if limit > 0.0 else math.inf
    return DesignCheck(name=name, clause=clause, demand=value, capacity=limit, ratio=ratio,
                       passed=bool(ratio <= 1.0), dimension=dimension, advisory=advisory,
                       details=dict(details or {}), notes=list(notes or []), status=OK, method=ctx.method,
                       inputs=dict(inputs or {}))


def _no_result(status, name, clause, reason, ctx, demand, dimension, combination, x, inputs, details, notes,
               advisory):
    return DesignCheck(name=name, clause=clause, demand=demand, capacity=None, ratio=None, passed=None,
                       combination=combination, x=x, dimension=dimension, advisory=advisory,
                       details=dict(details or {}), notes=list(notes or []), status=status, reason=reason,
                       method=ctx.method if ctx is not None else "", inputs=dict(inputs or {}))


def not_checked(name: str, clause: str, reason: str, ctx: DesignContext | None = None, *,
                demand: float | None = None, dimension: str = FORCE, combination: str = "",
                x: float | None = None, inputs: dict | None = None, details: dict | None = None,
                notes: list[str] | None = None, advisory: bool = False) -> DesignCheck:
    """The limit state is outside the implemented scope or data is missing."""
    return _no_result(NOT_CHECKED, name, clause, reason, ctx, demand, dimension, combination, x, inputs,
                      details, notes, advisory)


def refused(name: str, clause: str, reason: str, ctx: DesignContext | None = None, *,
            demand: float | None = None, dimension: str = FORCE, combination: str = "",
            x: float | None = None, inputs: dict | None = None, details: dict | None = None,
            notes: list[str] | None = None, advisory: bool = False) -> DesignCheck:
    """The input is contradictory or impossible: no result is given."""
    return _no_result(REFUSED, name, clause, reason, ctx, demand, dimension, combination, x, inputs,
                      details, notes, advisory)


def missing(ctx: DesignContext, *names: str) -> list[str]:
    """Names of context values that are None."""
    return [n for n in names if getattr(ctx, n) is None]


def positive(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value > 0
