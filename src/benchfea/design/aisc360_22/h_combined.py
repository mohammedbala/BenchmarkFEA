"""AISC 360-22 Chapter H -- combined forces and torsion (pp. 16.1-82 to 16.1-87),
with the utilities of Appendix 8, Section 8.1.2 (pp. 16.1-284, 16.1-285).

Implemented
-----------
* H1.1 flexure and compression, H1.2 flexure and tension, for doubly symmetric
  members (I-shapes, rectangular HSS / box, round HSS):

      Pr/Pc >= 0.2:  Pr/Pc + 8/9 (Mrx/Mcx + Mry/Mcy) <= 1.0        (Eq. H1-1a)
      Pr/Pc <  0.2:  Pr/(2 Pc) + (Mrx/Mcx + Mry/Mcy) <= 1.0        (Eq. H1-1b)

  All terms positive (User Note of H1.1).  Pc = available compressive strength
  of Chapter E (H1.1) or available tensile strength of Chapter D (H1.2);
  Mc = available flexural strength of Chapter F.

  Policy (the Specification does not say at which section the required
  strengths are taken, nor which Mc pairs with which section):

  - the interaction is evaluated AT EVERY STATION with the concurrent demands
    Pr(x), Mrx(x), Mry(x) of the combination, and the largest value governs;
  - against the member's GOVERNING available strengths of that combination:
    Pc = the smallest available strength of all Chapter E limit states
    (flexural buckling about both axes and torsional buckling; Commentary H1:
    Pc is governed by the maximum slenderness whatever the plane of bending)
    or of both Chapter D limit states; Mcx / Mcy = the smallest available
    strength of all Chapter F limit states of the axis, including the lowest
    LTB strength of all unbraced segments checked in the combination;
  - a station in compression uses H1.1, a station in tension H1.2 (each is
    reported as its own check); stations without axial force are evaluated
    with Pr = 0 in the check that exists (Eq. H1-1b; Commentary H1: H1-1 may be
    used for biaxial flexure without axial force);
  - with numbers only (no stations) the given maxima are taken as concurrent;
  - round sections: one flexural term, resultant moment / Mc (see f_flexure);
  - when a strength that the interaction needs is not available, the
    interaction is ``notChecked`` (never evaluated with a missing term).
  This is conservative with respect to pairing each station with the strength
  of its own unbraced segment.

  The check is produced when the combination has axial force and moment, or
  moments about both axes.  The optional H1.3 (Eq. H1-3) and the optional Cb
  increase of H1.2 (Eq. H1-2) are not used.  H2 and H4 are not implemented.

* Second-order effects.  Pr and Mr are "determined in accordance with Chapter
  C", i.e. they must contain second-order effects.  The checks use the analysis
  results AS GIVEN.  When the combination was analysed first-order the check
  says so in its notes (and the orchestration marks the design "conditional").
  The B1 / B2 amplification of Appendix 8 is NOT applied automatically: B2
  needs story data the model does not have, B1 applies to the no-translation
  moments only (Eq. A-8-1) which a single analysis does not separate, Cm with
  transverse load is "by analysis or 1.0", and EI* depends on the analysis
  method.  The formulas of 8.1.2 are provided as plain functions
  (:func:`cm_a_8_4`, :func:`pe1_a_8_5`, :func:`b1_a_8_3`) for callers that have
  that information.

* Torsion.  H3 is NOT implemented.  A torsional moment is neglected only when
  it is negligible by an ENGINE RULE (not an AISC rule): the St. Venant shear
  stress estimate tau_t must not exceed 1 % of 0.6 Fy (0.6 Fy being the shear
  yield stress of Eq. H3-8).  Then an advisory record of the screening is
  given.  Otherwise the member gets a ``notChecked`` check "torsion" with the
  torsional moment as demand, so it cannot be reported as passing.

  tau_t (derived here, uniform torsion of thin-walled sections):
  - open section (I-shape): for a narrow rectangle b x t under twist rate
    theta', Prandtl's stress function phi = G theta' (t^2/4 - n^2) (n across
    the thickness) satisfies laplace(phi) = -2 G theta' and phi = 0 on the long
    edges; tau = |d phi/dn| = 2 G theta' |n|, largest G theta' t at the faces;
    T = 2 integral(phi dA) = G theta' b t^3 / 3.  For several rectangles with a
    common theta': T = G J theta', J = sum(b t^3)/3, so in each rectangle
    tau_max = T t / J, largest in the thickest element:
        tau_t = T max(tf, tw) / J;
  - closed single cell (rectangular, round): the shear flow q = tau t is
    constant along the wall (longitudinal equilibrium of a wall element) and
    T = closed integral(q r ds) = 2 A0 q, A0 = area enclosed by the wall
    mid-line (Bredt): tau_t = T / (2 A0 t) with A0 = (H - t)(B - t) or
    pi (D - t)^2 / 4.
  Tests check both formulas by numerical integration.
"""
from __future__ import annotations

import math

import numpy as np

from ...results import DesignCheck
from ..base import LENGTH, NONE, OK, STRESS, TEXT, DesignContext, Demand, limit_check, not_checked, positive
from ..section import BOX, I_SHAPE, RECT_HSS, ROUND_HSS
from . import d_tension, e_compression, f_flexure

MOMENT = "moment"
FORCE = "force"
NAME_H1 = {"c": "combined flexure and compression", "t": "combined flexure and tension"}
CLAUSE_H1 = {"c": "AISC 360-22 H1.1", "t": "AISC 360-22 H1.2"}
NAME_TORSION = "torsion"
NAME_TORSION_SCREEN = "torsion: screening of the St. Venant shear stress (engine rule)"
CLAUSE_TORSION = "AISC 360-22 H3 (not implemented)"
TORSION_NEGLIGIBLE = 0.01               # engine rule: tau_t <= 1 % of 0.6 Fy


# ============================================================================== pure functions
def h1_interaction(Pr: float, Pc: float, Mrx: float, Mcx: float, Mry: float = 0.0, Mcy: float = math.inf):
    """Eq. H1-1a (Pr/Pc >= 0.2) or Eq. H1-1b -> (value, equation).  All terms positive."""
    axial = abs(Pr) / Pc
    flexure = abs(Mrx) / Mcx + abs(Mry) / Mcy
    if axial >= 0.2:
        return axial + 8.0 / 9.0 * flexure, "H1-1a"                 # Eq. H1-1a
    return axial / 2.0 + flexure, "H1-1b"                           # Eq. H1-1b


def cm_a_8_4(M1: float, M2: float) -> float:
    """Eq. A-8-4: Cm = 0.6 - 0.4 (M1/M2) for a beam-column WITHOUT transverse loading between supports in
    the plane of bending.  ``M1``, ``M2``: the end moments as SIGNED station moments (same sign convention at
    both ends); |M2| >= |M1| is arranged here.  With that convention equal signs mean single curvature, for
    which the Specification's ratio M1/M2 is negative."""
    if abs(M1) > abs(M2):
        M1, M2 = M2, M1
    if M2 == 0.0:
        raise ValueError("Eq. A-8-4 is undefined when both end moments are zero")
    return 0.6 - 0.4 * (-(M1 / M2))


def pe1_a_8_5(EI_star: float, Lc1: float) -> float:
    """Eq. A-8-5: Pe1 = pi^2 EI* / Lc1^2 (EI* = 0.8 tau_b EI with the direct analysis method, else EI)."""
    return math.pi ** 2 * EI_star / Lc1 ** 2


def b1_a_8_3(Cm: float, alpha: float, Pr: float, Pe1: float) -> float:
    """Eq. A-8-3: B1 = Cm / (1 - alpha Pr / Pe1) >= 1; ``math.inf`` when alpha Pr >= Pe1 (not covered by the
    equation: the member is unstable)."""
    den = 1.0 - alpha * Pr / Pe1
    return math.inf if den <= 0.0 else max(1.0, Cm / den)


def torsional_shear_stress(section, T: float):
    """St. Venant shear stress estimate (module docstring) -> (tau, formula) or (None, reason)."""
    d = section.dims
    if section.shape == I_SHAPE:
        if positive(section.J) and "tf" in d and "tw" in d:
            return abs(T) * max(d["tf"], d["tw"]) / section.J, "tau = T max(tf, tw) / J"
        return None, "J, tf and tw are needed"
    if section.shape in (RECT_HSS, BOX):
        if all(k in d for k in ("H", "B", "t")):
            return abs(T) / (2.0 * (d["H"] - d["t"]) * (d["B"] - d["t"]) * d["t"]), "tau = T / (2 (H - t)(B - t) t)"
        return None, "H, B and t are needed"
    if section.shape == ROUND_HSS:
        if all(k in d for k in ("D", "t")):
            return abs(T) / (2.0 * math.pi * (d["D"] - d["t"]) ** 2 / 4.0 * d["t"]), "tau = T / (2 pi (D - t)^2/4 t)"
        return None, "D and t are needed"
    return None, section.reason or "section outside the implemented scope"


# ============================================================================== H1
def _axial_strength(ctx: DesignContext, demand: Demand, kind: str):
    """-> (Pc, governing limit state, blocked reason)."""
    if kind == "c":
        checks = [c for fn in (e_compression.flexural_buckling_major, e_compression.flexural_buckling_minor,
                               e_compression.torsional_buckling) for c in fn(ctx, demand)]
    else:
        checks = [c for fn in (d_tension.yielding, d_tension.rupture) for c in fn(ctx, demand)]
    bad = [c for c in checks if c.status != OK]
    if bad:
        return None, "", f"{bad[0].name}: {bad[0].reason}"
    c = min(checks, key=lambda c: c.capacity)
    return c.capacity, c.name, ""


def interaction(ctx: DesignContext, demand: Demand) -> list:
    """H1.1 / H1.2: Eq. H1-1a / H1-1b at every station (policy in the module docstring)."""
    mo = f_flexure.moments(demand)
    if mo is None or ctx.section.shape is None:
        return []
    round_ = ctx.section.shape == ROUND_HSS
    if mo.x is not None:
        f = np.asarray(demand.extra["forces"], float)
        x0 = float(mo.x[0])
        xs, N = mo.x - x0, f[:, 0]
        Mx, My = (np.hypot(f[:, 5], f[:, 4]), np.zeros(len(xs))) if round_ else (np.abs(f[:, 5]), np.abs(f[:, 4]))
    else:                                                           # numbers: the maxima taken as concurrent
        xs = np.array([np.nan, np.nan])
        N = np.array([demand.Pt, -demand.Pc])
        Mx, My = (np.full(2, mo.Mr), np.zeros(2)) if round_ else (np.full(2, mo.Mrx), np.full(2, mo.Mry))
    has_m = bool(np.any(Mx > 0.0) or np.any(My > 0.0))
    has_n = bool(np.any(N != 0.0))
    if not has_m or not (has_n or (np.any(Mx > 0.0) and np.any(My > 0.0))):
        return []
    # --- available flexural strengths (Chapter F), governing over limit states and unbraced segments
    need = {}
    if round_:
        need["x"] = f_flexure.axis_strength(ctx, demand, "r")
    else:
        if np.any(Mx > 0.0):
            need["x"] = f_flexure.axis_strength(ctx, demand, "x")
        if np.any(My > 0.0):
            need["y"] = f_flexure.axis_strength(ctx, demand, "y")
    Mcx = need["x"].Mc if "x" in need else math.inf
    Mcy = need["y"].Mc if "y" in need else math.inf
    flex_blocked = [f"{_axis_name(a, round_)}: {r.blocked}" for a, r in need.items() if r.Mc is None]
    first_order = (demand.extra or {}).get("analysisType") == "firstOrder"
    out = []
    groups = [("c", N < 0.0), ("t", N > 0.0)]
    zero = N == 0.0
    if not np.any(N < 0.0) and not np.any(N > 0.0):
        groups = [("c", zero)]                                      # biaxial flexure without axial force
    else:
        k0 = 0 if np.any(N < 0.0) else 1
        groups[k0] = (groups[k0][0], groups[k0][1] | zero)
    for kind, mask in groups:
        if not np.any(mask):
            continue
        name, base = NAME_H1[kind], CLAUSE_H1[kind]
        idx = np.flatnonzero(mask)
        Pr_all = np.abs(N[idx])
        if not (np.any(Mx[idx] > 0.0) or np.any(My[idx] > 0.0)):
            continue                                                # axial force only in this group: Chapter D / E
        Pc, p_name, p_blocked = (math.inf, "", "") if not np.any(Pr_all > 0.0) else _axial_strength(ctx, demand, kind)
        blocked = flex_blocked + ([f"available axial strength: {p_blocked}"] if p_blocked else [])
        if blocked:
            out.append(not_checked(name, base, "a strength needed by the interaction is not available -- "
                                   + "; ".join(blocked), ctx, demand=None, dimension=NONE,
                                   combination=demand.combination))
            continue
        best = None
        for k in idx:
            value, eq = h1_interaction(abs(N[k]), Pc, Mx[k], Mcx, My[k], Mcy)
            if best is None or value > best[0]:
                best = (value, eq, int(k))
        value, eq, k = best
        Pr = float(abs(N[k]))
        notes = ["evaluated at every station with the concurrent Pr, Mrx, Mry against the member's governing "
                 "available strengths (engine policy; the Specification does not prescribe the location)"]
        if mo.x is None:
            notes[0] = "the given maxima of Pr, Mrx, Mry are taken as concurrent"
        if kind == "c" and first_order and Pr > 0.0:
            notes.append("required strengths from a FIRST-ORDER analysis: AISC 360-22 C1 and H1 require "
                         "second-order effects (P-Delta and P-delta); the B1/B2 amplification of Appendix 8 is "
                         "not applied")
        elif first_order:
            notes.append("required strengths from a first-order analysis (AISC 360-22 C1 requires second-order "
                         "effects in the analysis of the structure)")
        if round_:
            notes.append("round section: Mrx is the resultant moment, one flexural term")
        inputs = {"Pc": (None if math.isinf(Pc) else Pc, FORCE, f"derived: {p_name}" if p_name else "not needed "
                         "(no axial force)"),
                  "Mcx": (None if math.isinf(Mcx) else Mcx, MOMENT,
                          f"derived: {need['x'].governing}" if "x" in need else "not needed (no moment)")}
        if not round_:
            inputs["Mcy"] = (None if math.isinf(Mcy) else Mcy, MOMENT,
                             f"derived: {need['y'].governing}" if "y" in need else "not needed (no moment)")
        details = {"Pr": (Pr, FORCE), "Mrx": (float(Mx[k]), MOMENT), "Mry": (float(My[k]), MOMENT),
                   "Pr/Pc": (0.0 if math.isinf(Pc) else Pr / Pc, NONE),
                   "Mrx/Mcx": (float(Mx[k]) / Mcx, NONE), "Mry/Mcy": (float(My[k]) / Mcy, NONE),
                   "equation": (eq, TEXT)}
        chapter = "Chapter E" if kind == "c" else "Chapter D"
        out.append(DesignCheck(
            name=name, clause=f"{base}, Eq. {eq}; {chapter}; Chapter F", demand=value, capacity=1.0, ratio=value,
            passed=bool(value <= 1.0), combination=demand.combination,
            x=None if math.isnan(xs[k]) else float(xs[k]), dimension=NONE, details=details, notes=notes,
            status=OK, method=ctx.method, inputs={**ctx.inputs("method"), **inputs}))
    return out


def _axis_name(axis: str, round_: bool) -> str:
    return "flexure (resultant moment)" if round_ else f"flexure, {'major' if axis == 'x' else 'minor'} axis"


# ============================================================================== torsion
def torsion(ctx: DesignContext, demand: Demand) -> list:
    """H3 is not implemented: a torsional moment that is not negligible by the engine rule (module
    docstring) makes the member ``notChecked``; a negligible one is recorded as an advisory screening."""
    mo = f_flexure.moments(demand)
    if mo is None or mo.Tr == 0.0:
        return []
    kw = dict(demand=mo.Tr, dimension=MOMENT, combination=demand.combination, x=mo.x_Tr)
    rule = (f"engine rule (not an AISC rule): torsion is neglected only when the St. Venant shear stress estimate "
            f"does not exceed {TORSION_NEGLIGIBLE:g} x 0.6 Fy")
    if not positive(ctx.Fy):
        return [not_checked(NAME_TORSION, CLAUSE_TORSION, "the member carries a torsional moment and AISC 360-22 "
                            "H3 is not implemented; no Fy to screen it (" + rule + ")", ctx, **kw)]
    tau, how = torsional_shear_stress(ctx.section, mo.Tr)
    if tau is None:
        return [not_checked(NAME_TORSION, CLAUSE_TORSION, "the member carries a torsional moment and AISC 360-22 "
                            f"H3 is not implemented; the torsional shear stress cannot be screened ({how}); "
                            + rule, ctx, **kw)]
    limit = TORSION_NEGLIGIBLE * 0.6 * ctx.Fy
    inputs = {**ctx.inputs("Fy"), "T": (mo.Tr, MOMENT, "derived: analysis results")}
    details = {"tau_t": (tau, STRESS), "formula": (how, TEXT), "limit = 0.01 x 0.6 Fy": (limit, STRESS)}
    if tau <= limit:
        c = limit_check(NAME_TORSION_SCREEN, CLAUSE_TORSION + "; engine screening rule", ctx, value=tau,
                        limit=limit, dimension=STRESS, advisory=True, inputs=inputs, details=details,
                        notes=["the torsional moment is neglected in the member checks: " + rule,
                               "AISC 360-22 H3 (torsion) is not implemented"])
        c.combination, c.x = demand.combination, mo.x_Tr
        return [c]
    return [not_checked(NAME_TORSION, CLAUSE_TORSION, f"the member carries a torsional moment that is not "
                        f"negligible (tau_t = {tau:.4g} > {limit:.4g}; {rule}) and AISC 360-22 H3 (torsion and "
                        "combined torsion, flexure, shear and axial force) is not implemented", ctx,
                        inputs=inputs, details=details, **kw)]


LIMIT_STATES = (interaction, torsion)


def check(ctx: DesignContext, demand: Demand) -> list:
    """All Chapter H checks of one combination."""
    return [c for fn in LIMIT_STATES for c in fn(ctx, demand)]


__all__ = ["h1_interaction", "cm_a_8_4", "pe1_a_8_5", "b1_a_8_3", "torsional_shear_stress", "interaction",
           "torsion", "check", "LIMIT_STATES", "LENGTH"]
