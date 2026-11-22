"""AISC 360-22 Chapter E -- members in axial compression (pp. 16.1-38 to 16.1-49).

Scope: doubly symmetric I-shapes, rectangular HSS / box sections, round HSS.

* E1: phi_c = 0.90 (LRFD), Omega_c = 1.67 (ASD) for every limit state; Pn is the
  lowest value of the applicable limit states -> one check per limit state,
  the largest ratio governs.
* E2: Lc = K L per axis (Lcx = Kx Lx, Lcy = Ky Ly, Lcz = Kz Lz).  The User Note
  recommendation Lc/r <= 200 is an ADVISORY check (limit = design variable
  ``maxSlendernessCompression``).
* E3 flexural buckling, about the major (x) and the minor (y) axis separately:
  Fe = pi^2 E / (Lc/r)^2 (Eq. E3-4); Fn = 0.658^(Fy/Fe) Fy when
  Lc/r <= 4.71 sqrt(E/Fy) (Eq. E3-2), else Fn = 0.877 Fe (Eq. E3-3);
  Pn = Fn Ag (Eq. E3-1).
* E4 torsional buckling of doubly symmetric members twisting about the shear
  centre: Fe = [pi^2 E Cw / Lcz^2 + G J] / (Ix + Iy) (Eq. E4-2), Fn from
  Eq. E3-2 / E3-3 with that Fe (branch on Fy/Fe <= 2.25), Pn = Fn Ag (Eq. E4-1).
* E7 members with slender elements (Table B4.1a): Pn = Fn Ae (Eq. E7-1) with Fn
  of E3 / E4 (gross section) and, per limit state, Ae computed with THAT Fn:
  E7.1 effective widths (Eq. E7-2, E7-3, E7-5, Table E7.1) with
  Ae = Ag - sum (b - be) t (User Note of E7); E7.2 round HSS (Eq. E7-6, E7-7).

Choices where the Specification leaves room (all reported in the check):

* E3 branch test.  The Specification gives ``Lc/r <= 4.71 sqrt(E/Fy)`` and, in
  parentheses, ``Fy/Fe <= 2.25`` and says both give the same result; they differ
  in the sliver 4.71 < (Lc/r)/sqrt(E/Fy) <= 1.5 pi (the two formulas differ by
  0.05 % there).  Flexural buckling uses the Lc/r form (the primary one);
  torsional buckling, which has no Lc/r, uses Fy/Fe.  A note is added when the
  two tests disagree.
* E4 scope.  The Specification applies E4 to doubly symmetric members "when the
  torsional unbraced length exceeds the lateral unbraced length" without
  saying which lateral length.  The torsional mode is evaluated ALWAYS (it is a
  buckling mode of the member whatever the lengths are, so this can only
  lower Pn); the check states whether Lz exceeds Lx or Ly.  Closed sections
  (HSS, box) use Cw = 0 in Eq. E4-2 (lower bound of Fe).
* E7.1 c2: the values of Table E7.1 (1.31, 1.38, 1.49) are used, not the
  unrounded Eq. E7-4.  With the rounded c2, Eq. E7-3 can return slightly more
  than b just above the limit of Eq. E7-2; be is then taken as b.
* E7.1 width b of an element given only by its ratio: b = (b/t) t.
* E7.2: Eq. E7-7 exceeds Ag for D/t just above 0.11 E/Fy; Ae is not taken
  larger than Ag.  D/t >= 0.45 E/Fy: no provision -> notChecked.

Not covered: singly symmetric and unsymmetric members (E4-3, E4-4), bracing
offset from the shear centre (E4-10 to E4-12), single angles (E5), built-up
members (E6), Fe from a buckling analysis.
"""
from __future__ import annotations

import math

from ..base import AREA, FORCE, LENGTH, NONE, STRESS, TEXT, DesignContext, Demand, limit_check, missing, \
    not_checked, positive, strength_check
from ..section import I_SHAPE, ROUND_HSS
from . import b4_classification as b4

PHI_C, OMEGA_C = 0.90, 1.67                 # E1

NAME_FB = {"x": "compression: flexural buckling about the major axis (x)",
           "y": "compression: flexural buckling about the minor axis (y)"}
NAME_TB = "compression: torsional buckling"
NAME_SLENDER = "compression: slenderness Lc/r (advisory)"
CLAUSE_SLENDER = "AISC 360-22 E2 (User Note)"


# ============================================================================== pure functions
def elastic_flexural_stress(E: float, slenderness: float) -> float:
    """Eq. E3-4: Fe = pi^2 E / (Lc/r)^2."""
    return math.inf if slenderness == 0.0 else math.pi ** 2 * E / slenderness ** 2


def nominal_stress(Fy: float, Fe: float, inelastic: bool) -> float:
    """Eq. E3-2 (``inelastic``): Fn = 0.658^(Fy/Fe) Fy; Eq. E3-3: Fn = 0.877 Fe."""
    return 0.658 ** (Fy / Fe) * Fy if inelastic else 0.877 * Fe


def branch_by_slenderness(E: float, Fy: float, slenderness: float) -> bool:
    """E3(a) applies (Eq. E3-2) when Lc/r <= 4.71 sqrt(E/Fy)."""
    return slenderness <= 4.71 * math.sqrt(E / Fy)


def branch_by_stress(Fy: float, Fe: float) -> bool:
    """E3(a) applies (Eq. E3-2) when Fy/Fe <= 2.25."""
    return Fy / Fe <= 2.25


def torsional_stress(E: float, G: float, Cw: float, J: float, Ix: float, Iy: float, Lcz: float) -> float:
    """Eq. E4-2: Fe = [pi^2 E Cw / Lcz^2 + G J] / (Ix + Iy)."""
    warping = math.pi ** 2 * E * Cw / Lcz ** 2 if Cw > 0.0 else 0.0
    return (warping + G * J) / (Ix + Iy)


def effective_width(b: float, lam: float, lam_r: float, Fy: float, Fn: float, c1: float, c2: float):
    """E7.1: effective width of one element -> (be, Fel or None, equation).

    (a) lam <= lam_r sqrt(Fy/Fn): be = b (Eq. E7-2);
    (b) otherwise Fel = (c2 lam_r / lam)^2 Fy (Eq. E7-5) and
        be = b (1 - c1 sqrt(Fel/Fn)) sqrt(Fel/Fn) (Eq. E7-3), not taken larger than b.
    """
    if lam <= lam_r * math.sqrt(Fy / Fn):
        return b, None, "E7-2"
    Fel = (c2 * lam_r / lam) ** 2 * Fy                              # Eq. E7-5
    q = math.sqrt(Fel / Fn)
    return min(b, b * (1.0 - c1 * q) * q), Fel, "E7-3"              # Eq. E7-3


def round_effective_area(Ag: float, D_over_t: float, E: float, Fy: float):
    """E7.2 -> (Ae, equation) or (None, reason) when D/t >= 0.45 E/Fy (no provision)."""
    if D_over_t <= 0.11 * E / Fy:
        return Ag, "E7-6"                                           # Eq. E7-6
    if D_over_t < 0.45 * E / Fy:
        return min(Ag, (0.038 * E / (Fy * D_over_t) + 2.0 / 3.0) * Ag), "E7-7"      # Eq. E7-7
    return None, (f"D/t = {D_over_t:g} is not below 0.45 E/Fy = {0.45 * E / Fy:g}: AISC 360-22 E7.2 gives no "
                  "effective area for such round HSS")


# ============================================================================== helpers
def _prerequisites(ctx: DesignContext):
    """-> (classification, None) or (None, reason)."""
    sec = ctx.section
    if sec.shape is None:
        return None, sec.reason
    miss = missing(ctx, "E", "Fy")
    if miss:
        return None, f"no {' and '.join(miss)} (material or design variable)"
    if not positive(sec.A):
        return None, f"section {sec.id!r} has no gross area"
    cls = b4.classify_compression(sec, ctx.E, ctx.Fy)
    if cls.status != "ok":
        return None, cls.reason
    return cls, None


def _effective_area(ctx: DesignContext, cls: b4.Classification, Fn: float):
    """E7: effective area for the nominal stress Fn -> (Ae, details, notes, clause) or (None, reason, ...)."""
    sec, Ag = ctx.section, ctx.section.A
    details, notes = {}, []
    if sec.shape == ROUND_HSS:
        Ae, eq = round_effective_area(Ag, sec.elements["wall"].ratio, ctx.E, ctx.Fy)
        if Ae is None:
            return None, eq, notes, ""
        details[f"Ae (Eq. {eq})"] = (Ae, AREA)
        if eq == "E7-7" and Ae == Ag:
            notes.append("Eq. E7-7 gives more than Ag for this D/t: Ae taken as Ag")
        return Ae, details, notes, f"E7.2, Eq. {eq}"
    Ae, eqs = Ag, set()
    for e in cls.elements:
        _, c1, c2 = e.e7
        el = e.element
        width = el.width
        gate = e.lambda_r * math.sqrt(ctx.Fy / Fn)
        details[f"lambda_r sqrt(Fy/Fn) {e.label}"] = (gate, NONE)
        if e.ratio <= gate:                                         # Eq. E7-2: be = b
            eqs.add("E7-2")
            continue
        if el.thickness is None:
            return None, (f"the {e.label} of section {sec.id!r} is not fully effective (lambda = {e.ratio:g} > "
                          f"{gate:g}) but its thickness is not given, so the area reduction (b - be) t cannot "
                          "be computed"), notes, ""
        if width is None:
            width = e.ratio * el.thickness
            notes.append(f"width of the {e.label} taken as (width-to-thickness ratio) x t = {width:g}")
        be, Fel, eq = effective_width(width, e.ratio, e.lambda_r, ctx.Fy, Fn, c1, c2)
        eqs.add(eq)
        Ae -= el.count * (width - be) * el.thickness                # E7 User Note: Ae = Ag - sum (b - be) t
        details[f"c1, c2 {e.label} (Table E7.1 case {e.e7[0]})"] = (f"{c1:g}, {c2:g}", TEXT)
        details[f"Fel {e.label} (Eq. E7-5)"] = (Fel, STRESS)
        details[f"b {e.label}"] = (width, LENGTH)
        details[f"be {e.label} (Eq. {eq})"] = (be, LENGTH)
    details["Ae = Ag - sum (b - be) t (E7 User Note)"] = (Ae, AREA)
    return Ae, details, notes, "E7.1, Eq. " + ", ".join(sorted(eqs)) + (", E7-5, Table E7.1" if "E7-3" in eqs else "")


def _strength(ctx, cls, name, base_clause, Fn, Fe, demand: Demand, inputs, details, notes):
    """Pn = Fn Ag (nonslender, Eq. E3-1 / E4-1) or Fn Ae (slender, Eq. E7-1)."""
    sec = ctx.section
    kw = dict(demand=demand.Pc, combination=demand.combination, x=demand.x_Pc)
    inputs = {**inputs, **ctx.inputs("E", "Fy"), "Ag": (sec.A, AREA, sec.source("A")), **cls.inputs(),
              **sec.inputs("shape")}
    if sec.shape == I_SHAPE:
        inputs.update(sec.inputs("fabrication"))
    details = {**details, **cls.details(), "Fe": (Fe, STRESS), "Fn": (Fn, STRESS)}
    notes = [*notes, *sec.notes]
    if cls.slender:
        Ae, more, extra, e7_clause = _effective_area(ctx, cls, Fn)
        if Ae is None:
            return not_checked(name, f"AISC 360-22 E7; {base_clause}", more, ctx, inputs=inputs, details=details,
                               notes=notes, **kw)
        Pn = Fn * Ae                                                # Eq. E7-1
        details.update(more)
        details["Pn = Fn Ae (Eq. E7-1)"] = (Pn, FORCE)
        clause = f"AISC 360-22 E7, Eq. E7-1; {e7_clause}; {base_clause}; B4.1, Table B4.1a"
        notes = [*notes, *extra, "slender-element section (Table B4.1a): Section E7"]
    else:
        Pn = Fn * sec.A                                             # Eq. E3-1 / E4-1
        eq = "E4-1" if "E4" in base_clause else "E3-1"
        details[f"Pn = Fn Ag (Eq. {eq})"] = (Pn, FORCE)
        clause = f"AISC 360-22 {base_clause}; B4.1, Table B4.1a"
    return strength_check(name, clause, ctx, nominal=Pn, phi=PHI_C, omega=OMEGA_C, inputs=inputs,
                          details=details, notes=notes, **kw)


# ============================================================================== limit states
def _flexural(ctx: DesignContext, demand: Demand, axis: str) -> list:
    """E3 (E7 with slender elements): flexural buckling about AISC axis ``axis`` ("x" or "y")."""
    name = NAME_FB[axis]
    kw = dict(demand=demand.Pc, combination=demand.combination, x=demand.x_Pc)
    cls, reason = _prerequisites(ctx)
    if reason:
        return [not_checked(name, "AISC 360-22 E3", reason, ctx, **kw)]
    sec = ctx.section
    L, K, r = getattr(ctx, f"L{axis}"), getattr(ctx, f"K{axis}"), getattr(sec, f"r{axis}")
    if L is None or K is None or not positive(r):
        return [not_checked(name, "AISC 360-22 E3", f"unbraced length L{axis}, K{axis} or the radius of gyration "
                            f"r{axis} is not available", ctx, **kw)]
    Lc = K * L                                                      # E2: Lc = K L
    slender = Lc / r
    Fe = elastic_flexural_stress(ctx.E, slender)                    # Eq. E3-4
    inelastic = branch_by_slenderness(ctx.E, ctx.Fy, slender)       # E3(a) / E3(b)
    Fn = ctx.Fy if math.isinf(Fe) else nominal_stress(ctx.Fy, Fe, inelastic)
    eq = "E3-2" if inelastic else "E3-3"
    notes = []
    if not math.isinf(Fe) and inelastic != branch_by_stress(ctx.Fy, Fe):
        notes.append("the two forms of the E3 branch test disagree here (Lc/r <= 4.71 sqrt(E/Fy) vs "
                     "Fy/Fe <= 2.25); the Lc/r form is used (difference in Fn below 0.05 %)")
    inputs = {**ctx.inputs(f"L{axis}", f"K{axis}"), f"r{axis}": (r, LENGTH, sec.source("A"))}
    details = {f"Lc{axis} = K{axis} L{axis} (E2)": (Lc, LENGTH), f"Lc{axis}/r{axis}": (slender, NONE),
               "4.71 sqrt(E/Fy)": (4.71 * math.sqrt(ctx.E / ctx.Fy), NONE),
               "Fe equation": ("E3-4", TEXT), "Fn equation": (eq, TEXT)}
    return [_strength(ctx, cls, name, f"E3, Eq. {eq}, E3-4; E2", Fn, Fe, demand, inputs, details, notes)]


def flexural_buckling_major(ctx: DesignContext, demand: Demand) -> list:
    return _flexural(ctx, demand, "x")


def flexural_buckling_minor(ctx: DesignContext, demand: Demand) -> list:
    return _flexural(ctx, demand, "y")


def torsional_buckling(ctx: DesignContext, demand: Demand) -> list:
    """E4(a) (E7 with slender elements): torsional buckling of a doubly symmetric member, Eq. E4-2."""
    name, base = NAME_TB, "AISC 360-22 E4"
    kw = dict(demand=demand.Pc, combination=demand.combination, x=demand.x_Pc)
    cls, reason = _prerequisites(ctx)
    if reason:
        return [not_checked(name, base, reason, ctx, **kw)]
    sec = ctx.section
    need = [n for n in ("J", "Ix", "Iy") if not positive(getattr(sec, n))]
    if sec.shape == I_SHAPE and not positive(sec.Cw):
        need.append("Cw")
    if need:
        return [not_checked(name, base, f"section {sec.id!r} gives no {', '.join(need)} (needed by Eq. E4-2)",
                            ctx, **kw)]
    if not positive(ctx.G):
        return [not_checked(name, base, "no shear modulus G (material or design variable)", ctx, **kw)]
    if ctx.Lz is None or ctx.Kz is None or not ctx.Kz * ctx.Lz > 0.0:
        return [not_checked(name, base, "torsional unbraced length Lz / Kz not available", ctx, **kw)]
    Lcz = ctx.Kz * ctx.Lz                                           # Lcz = Kz Lz
    Cw = sec.Cw or 0.0
    Fe = torsional_stress(ctx.E, ctx.G, Cw, sec.J, sec.Ix, sec.Iy, Lcz)      # Eq. E4-2
    inelastic = branch_by_stress(ctx.Fy, Fe)
    Fn = nominal_stress(ctx.Fy, Fe, inelastic)                      # Eq. E3-2 / E3-3 with Fe of E4
    eq = "E3-2" if inelastic else "E3-3"
    laterals = [v for v in (ctx.Lx, ctx.Ly) if v is not None]
    exceeds = bool(laterals) and ctx.Lz > min(laterals)
    notes = ["torsional buckling is evaluated whatever the unbraced lengths are; AISC 360-22 E4 requires it "
             "for doubly symmetric members when the torsional unbraced length exceeds the lateral unbraced "
             f"length (here Lz {'exceeds' if exceeds else 'does not exceed'} the smaller of Lx, Ly)"]
    if sec.shape != I_SHAPE:
        notes.append("closed section: Cw = 0 in Eq. E4-2 (lower bound of Fe)")
    inputs = {**ctx.inputs("Lz", "Kz", "G"), **sec.inputs("J", "Cw", "Ix", "Iy")}
    details = {"Lcz = Kz Lz": (Lcz, LENGTH), "Fy/Fe": (ctx.Fy / Fe, NONE),
               "Fe equation": ("E4-2", TEXT), "Fn equation": (eq, TEXT)}
    return [_strength(ctx, cls, name, f"E4, Eq. E4-1, E4-2; E3, Eq. {eq}", Fn, Fe, demand, inputs, details,
                      notes)]


def slenderness(ctx: DesignContext, demand: Demand) -> list:
    """E2 User Note: Lc/r preferably not above 200 (advisory)."""
    sec = ctx.section
    ratios = {}
    for axis in "xy":
        L, K, r = getattr(ctx, f"L{axis}"), getattr(ctx, f"K{axis}"), getattr(sec, f"r{axis}")
        if L is not None and K is not None and positive(r):
            ratios[axis] = K * L / r
    if len(ratios) < 2:
        return [not_checked(NAME_SLENDER, CLAUSE_SLENDER, sec.reason or "unbraced lengths or radii of gyration "
                            "are not available", ctx, dimension=NONE, advisory=True)]
    worst = max(ratios, key=ratios.get)
    return [limit_check(NAME_SLENDER, CLAUSE_SLENDER, ctx, value=ratios[worst],
                        limit=ctx.maxSlendernessCompression,
                        inputs={**ctx.inputs("Lx", "Ly", "Kx", "Ky", "maxSlendernessCompression"),
                                "rx": (sec.rx, LENGTH, sec.source("A")), "ry": (sec.ry, LENGTH, sec.source("A"))},
                        details={"Lcx/rx": (ratios["x"], NONE), "Lcy/ry": (ratios["y"], NONE),
                                 "governing axis": (worst, TEXT)},
                        notes=["advisory: AISC 360-22 E2 sets no maximum slenderness for compression members"])]


LIMIT_STATES = (flexural_buckling_major, flexural_buckling_minor, torsional_buckling, slenderness)


def check(ctx: DesignContext, demand: Demand) -> list:
    """All Chapter E checks of one combination."""
    return [c for fn in LIMIT_STATES for c in fn(ctx, demand)]
