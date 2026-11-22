"""AISC 360-22 Chapter D -- members in axial tension (pp. 16.1-32 to 16.1-35).

* D1: no maximum slenderness; the User Note recommends L/r <= 300 with L the
  length of the member as fabricated and r the least radius of gyration.
  Reported as an ADVISORY check (limit = design variable
  ``maxSlendernessTension``); L is taken as the design member length.
* D2(a) tensile yielding in the gross section:  Pn = Fy Ag  (Eq. D2-1),
  phi_t = 0.90, Omega_t = 1.67.
* D2(b) tensile rupture in the net section:  Pn = Fu Ae  (Eq. D2-2),
  phi_t = 0.75, Omega_t = 2.00, with Ae = An U (Eq. D3-1).

``An`` (B4.3) and ``U`` (Table D3.1) depend on the end connection, which a
frame model does not describe: they are design variables.  Their defaults
(An = Ag, U = 1.0) are reported as defaults in the inputs of the check.  The
lower bound on U of Section D3 and Table D3.1 itself are not evaluated.

The two limit states are reported as two checks; the available tensile
strength of D2 is the lower of the two AVAILABLE strengths, i.e. the check
with the larger ratio governs.

Chapter D applies to any cross-section: only A (and the radii of gyration for
D1) are needed.  Not covered: D4 (built-up members), D5 (pin-connected
members), D6 (eyebars), block shear (J4.3), eccentric connections.
"""
from __future__ import annotations

from ..base import AREA, FORCE, NONE, STRESS, DesignContext, Demand, limit_check, not_checked, positive, \
    refused, strength_check

PHI_YIELD, OMEGA_YIELD = 0.90, 1.67         # D2(a)
PHI_RUPTURE, OMEGA_RUPTURE = 0.75, 2.00     # D2(b)

NAME_YIELD = "tension: yielding in the gross section"
NAME_RUPTURE = "tension: rupture in the net section"
NAME_SLENDER = "tension: slenderness L/r (advisory)"
CLAUSE_YIELD = "AISC 360-22 D2(a), Eq. D2-1"
CLAUSE_RUPTURE = "AISC 360-22 D2(b), Eq. D2-2; D3, Eq. D3-1"
CLAUSE_SLENDER = "AISC 360-22 D1 (User Note)"


def _area(ctx: DesignContext):
    A = ctx.section.A
    return A if positive(A) else None


def yielding(ctx: DesignContext, demand: Demand) -> list:
    """D2(a): Pn = Fy Ag (Eq. D2-1); phi_t = 0.90 (LRFD), Omega_t = 1.67 (ASD)."""
    kw = dict(demand=demand.Pt, combination=demand.combination, x=demand.x_Pt)
    Ag = _area(ctx)
    if Ag is None:
        return [not_checked(NAME_YIELD, CLAUSE_YIELD, ctx.section.reason or "the section has no gross area",
                            ctx, **kw)]
    if not positive(ctx.Fy):
        return [not_checked(NAME_YIELD, CLAUSE_YIELD, "no yield stress Fy (material or design variable)",
                            ctx, **kw)]
    Pn = ctx.Fy * Ag                                               # Eq. D2-1
    return [strength_check(NAME_YIELD, CLAUSE_YIELD, ctx, nominal=Pn, phi=PHI_YIELD, omega=OMEGA_YIELD,
                           inputs={**ctx.inputs("Fy"), "Ag": (Ag, AREA, ctx.section.source("A"))},
                           details={"Pn = Fy Ag (Eq. D2-1)": (Pn, FORCE)}, **kw)]


def rupture(ctx: DesignContext, demand: Demand) -> list:
    """D2(b): Pn = Fu Ae (Eq. D2-2), Ae = An U (Eq. D3-1); phi_t = 0.75, Omega_t = 2.00."""
    kw = dict(demand=demand.Pt, combination=demand.combination, x=demand.x_Pt)
    Ag = _area(ctx)
    if Ag is None:
        return [not_checked(NAME_RUPTURE, CLAUSE_RUPTURE, ctx.section.reason or "the section has no gross area",
                            ctx, **kw)]
    if not positive(ctx.Fu):
        return [not_checked(NAME_RUPTURE, CLAUSE_RUPTURE,
                            "no tensile strength Fu (material or design variable)", ctx, **kw)]
    An, an_src = (Ag, "default: net area = gross area (no holes, AISC 360-22 B4.3b)") if ctx.An is None \
        else (ctx.An, ctx.source("An"))
    inputs = {**ctx.inputs("Fu", "U"), "An": (An, AREA, an_src), "Ag": (Ag, AREA, ctx.section.source("A"))}
    if not positive(An) or An > Ag * (1.0 + 1e-9):
        return [refused(NAME_RUPTURE, CLAUSE_RUPTURE, f"net area An = {An:g} must be positive and not larger "
                        f"than the gross area Ag = {Ag:g} (AISC 360-22 B4.3)", ctx, inputs=inputs, **kw)]
    if not positive(ctx.U) or ctx.U > 1.0:
        return [refused(NAME_RUPTURE, CLAUSE_RUPTURE, f"shear lag factor U = {ctx.U:g} must be in (0, 1] "
                        "(AISC 360-22 Table D3.1)", ctx, inputs=inputs, **kw)]
    Ae = An * ctx.U                                                # Eq. D3-1
    Pn = ctx.Fu * Ae                                               # Eq. D2-2
    notes = []
    if ctx.source("U") == "default" or ctx.An is None or ctx.source("An") == "default":
        notes.append("An and/or U are defaults (An = Ag, U = 1.0): they depend on the end connection "
                     "(B4.3b, Table D3.1), which the model does not describe")
    return [strength_check(NAME_RUPTURE, CLAUSE_RUPTURE, ctx, nominal=Pn, phi=PHI_RUPTURE, omega=OMEGA_RUPTURE,
                           inputs=inputs, details={"Ae = An U (Eq. D3-1)": (Ae, AREA),
                                                   "Pn = Fu Ae (Eq. D2-2)": (Pn, FORCE)}, notes=notes, **kw)]


def slenderness(ctx: DesignContext, demand: Demand) -> list:
    """D1 User Note: L / r_min preferably not above 300 (advisory; no mandatory limit)."""
    sec = ctx.section
    if sec.rx is None or sec.ry is None:
        return [not_checked(NAME_SLENDER, CLAUSE_SLENDER,
                            sec.reason or "radii of gyration of the section are not available", ctx,
                            dimension=NONE, advisory=True)]
    r = min(sec.rx, sec.ry)
    return [limit_check(NAME_SLENDER, CLAUSE_SLENDER, ctx, value=ctx.L / r, limit=ctx.maxSlendernessTension,
                        inputs={**ctx.inputs("L", "maxSlendernessTension"),
                                "r_min": (r, "length", sec.source("A"))},
                        details={"L/r": (ctx.L / r, NONE)},
                        notes=["advisory: AISC 360-22 D1 sets no maximum slenderness for tension members; "
                               "L = design member length; does not apply to rods"])]


LIMIT_STATES = (yielding, rupture, slenderness)


def check(ctx: DesignContext, demand: Demand) -> list:
    """All Chapter D checks of one combination."""
    return [c for fn in LIMIT_STATES for c in fn(ctx, demand)]


__all__ = ["yielding", "rupture", "slenderness", "check", "LIMIT_STATES", "STRESS"]
