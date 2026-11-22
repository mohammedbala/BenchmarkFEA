"""AISC 360-22 Chapter G -- members in shear (pp. 16.1-75 to 16.1-81).

G1: phi_v = 0.90 (LRFD), Omega_v = 1.67 (ASD), except G2.1(a).

* G2.1 webs of I-shapes, shear along local y:  Vn = 0.6 Fy Aw Cv1 (Eq. G2-1),
  Aw = d tw.
  (a) rolled I-shape with h/tw <= 2.24 sqrt(E/Fy): Cv1 = 1.0 (Eq. G2-2) with
      phi_v = 1.00, Omega_v = 1.50;
  (b) all other I-shapes: Cv1 = 1.0 when h/tw <= 1.10 sqrt(kv E/Fy) (Eq. G2-3),
      else Cv1 = 1.10 sqrt(kv E/Fy) / (h/tw) (Eq. G2-4); webs WITHOUT
      transverse stiffeners: kv = 5.34 (the model has no stiffeners).
  h is the web depth of the section's design view (rolled: clear distance
  between flanges less the fillets; built-up welded: clear distance between
  flanges).  G2.1(b)(1) words h only for built-up sections; for a rolled shape
  that does not meet (a) the same h as in (a) is used.
* G4 rectangular HSS and box sections, shear along local y or z:
  Vn = 0.6 Fy Aw Cv2 (Eq. G4-1), Aw = 2 h t with h = the flat width of the two
  walls PARALLEL to the shear force ("width resisting the shear force": HSS:
  outside dimension - 3 t when the corner radius is not known; box: clear
  distance between flanges), Cv2 of G2.2 (Eq. G2-9 to G2-11) with kv = 5.
* G5 round HSS: Vn = Fcr Ag / 2 (Eq. G5-1), Fcr = larger of Eq. G5-2a and G5-2b
  but not more than 0.6 Fy.  A round section has no preferred direction, so
  the check uses the RESULTANT shear force sqrt(Vy^2 + Vz^2).  Lv = distance
  from maximum to zero shear force comes from the shear diagram of the
  combination; when the shear does not vanish along the member Lv is
  undefined and Eq. G5-2a is not credited (Fcr = Eq. G5-2b <= 0.6 Fy, a lower
  bound of "the larger of").  When Eq. G5-2b alone reaches 0.6 Fy, Lv is not
  needed.
* G6 minor-axis shear of doubly symmetric I-shapes, shear along local z:
  Vn = 0.6 Fy bf tf Cv2 (Eq. G6-1) for EACH flange, Cv2 of G2.2 with
  h/tw = bf / (2 tf) and kv = 1.2; the member strength is the sum of the two
  equal flanges, 2 x Eq. G6-1.

Each direction is checked on its own (Chapter G has no interaction between
shear in two directions for I-shapes and rectangular sections).  Not used:
tension field action (G2.2, G2.3), transverse stiffeners (G2.4), G3, G7 (web
openings), torsion (H3).
"""
from __future__ import annotations

import math

from ..base import AREA, FORCE, LENGTH, NONE, STRESS, TEXT, DesignContext, Demand, missing, not_checked, \
    positive, strength_check
from ..section import BOX, I_SHAPE, RECT_HSS, ROLLED, ROUND_HSS

PHI_V, OMEGA_V = 0.90, 1.67                 # G1(a)
PHI_V_ROLLED, OMEGA_V_ROLLED = 1.00, 1.50   # G2.1(a)
KV_UNSTIFFENED = 5.34                       # G2.1(b)(2)(i)
KV_HSS = 5.0                                # G4
KV_FLANGE = 1.2                             # G6

NAME_VY = "shear along local y"
NAME_VZ = "shear along local z"
NAME_V = "shear (resultant, round section)"
CLAUSE = "AISC 360-22 Chapter G"


# ============================================================================== pure functions
def cv1(h_tw: float, kv: float, E: float, Fy: float):
    """G2.1(b)(1): web shear strength coefficient -> (Cv1, equation)."""
    limit = 1.10 * math.sqrt(kv * E / Fy)
    if h_tw <= limit:
        return 1.0, "G2-3"                                          # Eq. G2-3
    return limit / h_tw, "G2-4"                                     # Eq. G2-4


def cv2(h_tw: float, kv: float, E: float, Fy: float):
    """G2.2: web shear buckling strength coefficient -> (Cv2, equation)."""
    root = math.sqrt(kv * E / Fy)
    if h_tw <= 1.10 * root:
        return 1.0, "G2-9"                                          # Eq. G2-9
    if h_tw <= 1.37 * root:
        return 1.10 * root / h_tw, "G2-10"                          # Eq. G2-10
    return 1.51 * kv * E / (h_tw ** 2 * Fy), "G2-11"                # Eq. G2-11


def round_fcr(E: float, Fy: float, D_over_t: float, Lv_over_D: float | None):
    """G5: Fcr = larger of Eq. G5-2a and G5-2b, not more than 0.6 Fy -> (Fcr, Fcr_a or None, Fcr_b).
    ``Lv_over_D = None``: Eq. G5-2a not evaluated."""
    fb = 0.78 * E / D_over_t ** 1.5                                 # Eq. G5-2b
    fa = None
    if Lv_over_D is not None and Lv_over_D > 0.0:
        fa = 1.60 * E / (math.sqrt(Lv_over_D) * D_over_t ** 1.25)   # Eq. G5-2a
    return min(0.6 * Fy, max(fb, fa) if fa is not None else fb), fa, fb


# ============================================================================== limit states
def _basic(ctx: DesignContext):
    if ctx.section.shape is None:
        return ctx.section.reason
    miss = missing(ctx, "E", "Fy")
    if miss:
        return f"no {' and '.join(miss)} (material or design variable)"
    return None


def _web_g21(ctx: DesignContext, demand: Demand):
    """G2.1: web of an I-shape, shear along local y."""
    sec = ctx.section
    kw = dict(demand=demand.Vy, combination=demand.combination, x=demand.x_Vy)
    web, d, tw = sec.elements["web"], sec.dims.get("d"), sec.dims.get("tw")
    if d is None or tw is None:
        return not_checked(NAME_VY, "AISC 360-22 G2.1", f"section {sec.id!r} gives no d and tw (Aw = d tw)",
                           ctx, **kw)
    if web.ratio is None:
        return not_checked(NAME_VY, "AISC 360-22 G2.1", f"h/tw of section {sec.id!r} is not available: "
                           f"{web.definition}", ctx, **kw)
    E, Fy, lam = ctx.E, ctx.Fy, web.ratio
    Aw = d * tw                                                     # G2.1: Aw = d tw
    limit_a = 2.24 * math.sqrt(E / Fy)
    notes = [*sec.notes]
    details = {"Aw = d tw": (Aw, AREA), "h/tw": (lam, NONE), "width definition web": (web.definition, TEXT),
               "2.24 sqrt(E/Fy)": (limit_a, NONE)}
    if sec.fabrication == ROLLED and lam <= limit_a:                # G2.1(a)
        c, eq, phi, omega = 1.0, "G2-2", PHI_V_ROLLED, OMEGA_V_ROLLED
        clause = "AISC 360-22 G2.1(a), Eq. G2-1, G2-2"
    else:                                                           # G2.1(b), no transverse stiffeners
        kv = KV_UNSTIFFENED
        c, eq = cv1(lam, kv, E, Fy)
        phi, omega = PHI_V, OMEGA_V
        clause = f"AISC 360-22 G2.1(b), Eq. G2-1, {eq}; G1"
        details["kv (web without transverse stiffeners)"] = (kv, NONE)
        details["1.10 sqrt(kv E/Fy)"] = (1.10 * math.sqrt(kv * E / Fy), NONE)
        notes.append("web without transverse stiffeners: kv = 5.34; tension field action (G2.2, G2.3) not used")
        if sec.fabrication == ROLLED:
            notes.append("rolled shape outside G2.1(a): h of G2.1(a) is also used in G2.1(b)")
    if lam > 260.0:
        notes.append("h/tw exceeds 260, the proportioning limit of AISC 360-22 F13.2 for unstiffened girders "
                     "(not a Chapter G check)")
    Vn = 0.6 * Fy * Aw * c                                          # Eq. G2-1
    details.update({f"Cv1 (Eq. {eq})": (c, NONE), "Vn = 0.6 Fy Aw Cv1 (Eq. G2-1)": (Vn, FORCE)})
    inputs = {**ctx.inputs("E", "Fy"), **sec.inputs("d", "tw", "fabrication"),
              "h/tw": (lam, NONE, web.source)}
    return strength_check(NAME_VY, clause, ctx, nominal=Vn, phi=phi, omega=omega, inputs=inputs,
                          details=details, notes=notes, **kw)


def _flanges_g6(ctx: DesignContext, demand: Demand):
    """G6: flanges of a doubly symmetric I-shape, shear along local z."""
    sec = ctx.section
    kw = dict(demand=demand.Vz, combination=demand.combination, x=demand.x_Vz)
    bf, tf = sec.dims.get("bf"), sec.dims.get("tf")
    if bf is None or tf is None:
        return not_checked(NAME_VZ, "AISC 360-22 G6", f"section {sec.id!r} gives no bf and tf", ctx, **kw)
    lam = bf / (2.0 * tf)                                           # G6: h/tw = bf / (2 tf) for I-shapes
    c, eq = cv2(lam, KV_FLANGE, ctx.E, ctx.Fy)
    Vn1 = 0.6 * ctx.Fy * bf * tf * c                                # Eq. G6-1, one flange
    Vn = 2.0 * Vn1
    details = {"bf/(2 tf)": (lam, NONE), "kv": (KV_FLANGE, NONE),
               "1.10 sqrt(kv E/Fy)": (1.10 * math.sqrt(KV_FLANGE * ctx.E / ctx.Fy), NONE),
               f"Cv2 (Eq. {eq})": (c, NONE), "Vn per flange = 0.6 Fy bf tf Cv2 (Eq. G6-1)": (Vn1, FORCE),
               "Vn = 2 flanges": (Vn, FORCE)}
    return strength_check(NAME_VZ, f"AISC 360-22 G6, Eq. G6-1; G2.2, Eq. {eq}; G1", ctx, nominal=Vn, phi=PHI_V,
                          omega=OMEGA_V, inputs={**ctx.inputs("E", "Fy"), **sec.inputs("bf", "tf")},
                          details=details, notes=["Eq. G6-1 is per shear-resisting element; the two equal "
                                                  "flanges are added", *sec.notes], **kw)


def _walls_g4(ctx: DesignContext, demand: Demand, direction: str):
    """G4: rectangular HSS / box, shear along local ``direction`` ("y" or "z")."""
    sec = ctx.section
    name = NAME_VY if direction == "y" else NAME_VZ
    V, x = (demand.Vy, demand.x_Vy) if direction == "y" else (demand.Vz, demand.x_Vz)
    kw = dict(demand=V, combination=demand.combination, x=x)
    wall = sec.elements[f"wall_{direction}"]                        # the two walls parallel to the shear force
    t = wall.thickness
    if wall.ratio is None or t is None:
        return not_checked(name, "AISC 360-22 G4", f"section {sec.id!r}: flat width and thickness of the walls "
                           f"parallel to local {direction} are not available ({wall.definition})", ctx, **kw)
    notes = [*sec.notes]
    h = wall.width
    if h is None:
        h = wall.ratio * t
        notes.append(f"h taken as (h/t) x t = {h:g}")
    c, eq = cv2(wall.ratio, KV_HSS, ctx.E, ctx.Fy)                  # G4: Cv2 of G2.2 with h/tw = h/t, kv = 5
    Aw = 2.0 * h * t                                                # G4: Aw = 2 h t
    Vn = 0.6 * ctx.Fy * Aw * c                                      # Eq. G4-1
    details = {"h": (h, LENGTH), "width definition": (wall.definition, TEXT), "h/t": (wall.ratio, NONE),
               "kv": (KV_HSS, NONE), "1.10 sqrt(kv E/Fy)": (1.10 * math.sqrt(KV_HSS * ctx.E / ctx.Fy), NONE),
               f"Cv2 (Eq. {eq})": (c, NONE), "Aw = 2 h t": (Aw, AREA), "Vn = 0.6 Fy Aw Cv2 (Eq. G4-1)": (Vn, FORCE)}
    inputs = {**ctx.inputs("E", "Fy"), **sec.inputs("shape"), "t": (t, LENGTH, sec.source("t")),
              "h/t": (wall.ratio, NONE, wall.source)}
    return strength_check(name, f"AISC 360-22 G4, Eq. G4-1; G2.2, Eq. {eq}; G1", ctx, nominal=Vn, phi=PHI_V,
                          omega=OMEGA_V, inputs=inputs, details=details, notes=notes, **kw)


def _round_g5(ctx: DesignContext, demand: Demand):
    """G5: round HSS, resultant shear."""
    sec = ctx.section
    kw = dict(demand=demand.V, combination=demand.combination, x=demand.x_V)
    wall = sec.elements["wall"]
    if wall.ratio is None or not positive(sec.A):
        return not_checked(NAME_V, "AISC 360-22 G5", f"section {sec.id!r}: D/t or the gross area is not "
                           "available", ctx, **kw)
    D = wall.width
    notes = [*sec.notes, "resultant shear sqrt(Vy^2 + Vz^2) (a round section has no preferred direction)"]
    ratio = None
    if demand.Lv is not None and D is not None:
        ratio = demand.Lv / D
    Fcr, fa, fb = round_fcr(ctx.E, ctx.Fy, wall.ratio, ratio)
    cap = 0.6 * ctx.Fy
    details = {"D/t": (wall.ratio, NONE), "Fcr (Eq. G5-2b)": (fb, STRESS), "0.6 Fy": (cap, STRESS)}
    if fa is not None:
        details["Lv"] = (demand.Lv, LENGTH)
        details["Fcr (Eq. G5-2a)"] = (fa, STRESS)
    elif fb < cap:
        why = "the shear force does not vanish along the member, so Lv is undefined" if demand.Lv is None \
            else "the outside diameter D is not given"
        notes.append(f"Eq. G5-2a not credited ({why}): Fcr = Eq. G5-2b, a lower bound of the larger of "
                     "Eq. G5-2a and G5-2b")
    else:
        notes.append("Eq. G5-2b alone reaches 0.6 Fy: Fcr = 0.6 Fy whatever Lv is")
    Vn = Fcr * sec.A / 2.0                                          # Eq. G5-1
    details.update({"Fcr": (Fcr, STRESS), "Vn = Fcr Ag / 2 (Eq. G5-1)": (Vn, FORCE)})
    inputs = {**ctx.inputs("E", "Fy"), "Ag": (sec.A, AREA, sec.source("A")),
              "D/t": (wall.ratio, NONE, wall.source), "Lv": (demand.Lv, LENGTH, "derived: shear diagram")}
    return strength_check(NAME_V, "AISC 360-22 G5, Eq. G5-1, G5-2a, G5-2b; G1", ctx, nominal=Vn, phi=PHI_V,
                          omega=OMEGA_V, inputs=inputs, details=details, notes=notes, **kw)


def shear(ctx: DesignContext, demand: Demand) -> list:
    """The Chapter G checks of the section's shape: I-shape -> G2.1 (local y) and G6 (local z);
    rectangular HSS / box -> G4 in both directions; round HSS -> G5 on the resultant."""
    reason = _basic(ctx)
    shape = ctx.section.shape
    if reason:
        c = demand.combination
        if shape == ROUND_HSS:
            return [not_checked(NAME_V, CLAUSE, reason, ctx, demand=demand.V, combination=c, x=demand.x_V)]
        return [not_checked(NAME_VY, CLAUSE, reason, ctx, demand=demand.Vy, combination=c, x=demand.x_Vy),
                not_checked(NAME_VZ, CLAUSE, reason, ctx, demand=demand.Vz, combination=c, x=demand.x_Vz)]
    if shape == I_SHAPE:
        return [_web_g21(ctx, demand), _flanges_g6(ctx, demand)]
    if shape in (RECT_HSS, BOX):
        return [_walls_g4(ctx, demand, "y"), _walls_g4(ctx, demand, "z")]
    return [_round_g5(ctx, demand)]


LIMIT_STATES = (shear,)


def check(ctx: DesignContext, demand: Demand) -> list:
    """All Chapter G checks of one combination."""
    return [c for fn in LIMIT_STATES for c in fn(ctx, demand)]
