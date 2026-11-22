"""AISC 360-22 Section B4.1 -- classification of sections for local buckling.

Implemented from the Specification (pp. 16.1-19 to 16.1-24):

* Table B4.1a (axial compression): lambda_r for Cases 1, 2 (with k_c, footnote
  [a]), 5, 6, 8 and 9;
* Table B4.1b (flexure): lambda_p and lambda_r for Cases 10, 11 (k_c, F_L,
  footnotes [a], [b]), 13, 15, 17, 19, 20 and 21;
* the width definitions of B4.1a / B4.1b are applied in
  :mod:`benchfea.design.section` (which states, per element, how the width was
  obtained); the thickness of HSS is the DESIGN wall thickness of B4.2, which
  is an input of the section.

Rules of the comparison (B4.1): an element whose ratio does not exceed the
limit is in the more favourable class, i.e. ``lambda == lambda_r`` is
nonslender / noncompact and ``lambda == lambda_p`` is compact.  No numerical
tolerance is applied at the boundaries.

Shape -> case mapping (compression / flexure):

============== ======================= =========================================
shape          Table B4.1a             Table B4.1b
============== ======================= =========================================
rolled I       flange 1, web 5         flange 10 (major) / 13 (minor), web 15
built-up I     flange 2, web 5         flange 11 (major) / 13 (minor), web 15
rect. HSS      walls 6                 flanges 17, webs 19
box section    walls 8                 flanges 21, webs 19
round HSS      9                       20
============== ======================= =========================================

For minor-axis flexure of an I-shape the web lies on the neutral axis and
Table B4.1b has no case for it: only the flanges are classified.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from ..section import BOX, BUILT_UP, I_SHAPE, RECT_HSS, ROLLED, ROUND_HSS, DesignSection, PlateElement

NONSLENDER, SLENDER = "nonslender", "slender"
COMPACT, NONCOMPACT = "compact", "noncompact"
MAJOR, MINOR = "major", "minor"

# Table E7.1 (p. 16.1-49): effective width imperfection adjustment factors (c1, c2), by element type.
E7_STIFFENED = ("(a)", 0.18, 1.31)      # stiffened elements except walls of square and rectangular HSS
E7_HSS_WALL = ("(b)", 0.20, 1.38)       # walls of square and rectangular HSS
E7_OTHER = ("(c)", 0.22, 1.49)          # all other elements


@dataclass(frozen=True)
class ElementClass:
    """Classification of one plate element."""
    key: str                          # element key of the DesignSection
    label: str                        # e.g. "flange", "web", "wall (depth)"
    case: int                         # case number of Table B4.1a / B4.1b
    table: str                        # "B4.1a" | "B4.1b"
    ratio: float
    lambda_p: float | None            # flexure only
    lambda_r: float
    cls: str                          # "nonslender" | "slender" | "compact" | "noncompact"
    element: PlateElement
    e7: tuple | None = None           # (Table E7.1 case, c1, c2) for compression elements
    kc: float | None = None           # Table B4.1a/b footnote [a]
    limit_text: str = ""              # the limit as a formula

    @property
    def citation(self) -> str:
        return f"Table {self.table} Case {self.case}"


@dataclass(frozen=True)
class Classification:
    """Result for a section: ``status`` "ok" or "notChecked" (with ``reason``)."""
    status: str
    elements: tuple = ()
    reason: str = ""

    @property
    def cls(self) -> str | None:
        """Class of the section = the least favourable element class."""
        if self.status != "ok":
            return None
        order = [NONSLENDER, COMPACT, NONCOMPACT, SLENDER]
        return max((e.cls for e in self.elements), key=order.index)

    @property
    def slender(self) -> bool:
        return self.cls == SLENDER

    def details(self) -> dict:
        """Intermediate values for a DesignCheck."""
        out = {}
        for e in self.elements:
            out[f"lambda {e.label}"] = (e.ratio, "dimensionless")
            if e.lambda_p is not None:
                out[f"lambda_p {e.label} ({e.citation})"] = (e.lambda_p, "dimensionless")
            out[f"lambda_r {e.label} ({e.citation})"] = (e.lambda_r, "dimensionless")
            if e.kc is not None:
                out[f"kc {e.label} (Table {e.table} footnote [a])"] = (e.kc, "dimensionless")
            out[f"class {e.label}"] = (e.cls, "n/a")
            out[f"width definition {e.label}"] = (e.element.definition, "n/a")
        return out

    def inputs(self) -> dict:
        return {f"lambda {e.label}": (e.ratio, "dimensionless", e.element.source) for e in self.elements}


def kc_factor(h_over_tw: float) -> float:
    """Tables B4.1a / B4.1b footnote [a]: kc = 4 / sqrt(h/tw), not less than 0.35 nor greater than 0.76."""
    return min(0.76, max(0.35, 4.0 / math.sqrt(h_over_tw)))


def _unavailable(sec: DesignSection, E, Fy, keys) -> str | None:
    if sec.shape is None:
        return sec.reason
    if E is None or Fy is None:
        return "E and Fy are needed for the width-to-thickness limits (AISC 360-22 Tables B4.1a/B4.1b)"
    for key in keys:
        el = sec.elements[key]
        if el.ratio is None:
            return (f"width-to-thickness ratio of the {key.replace('_', ' ')} of section {sec.id!r} is not "
                    f"available: {el.definition}")
        if not el.ratio > 0.0:
            return f"width-to-thickness ratio of the {key.replace('_', ' ')} of section {sec.id!r} is not positive"
    return None


# ============================================================================== Table B4.1a
def classify_compression(sec: DesignSection, E: float | None, Fy: float | None) -> Classification:
    """Table B4.1a: nonslender / slender elements of a member in axial compression."""
    keys = {I_SHAPE: ("flange", "web"), RECT_HSS: ("wall_y", "wall_z"), BOX: ("wall_y", "wall_z"),
            ROUND_HSS: ("wall",)}.get(sec.shape, ())
    reason = _unavailable(sec, E, Fy, keys)
    if reason:
        return Classification("notChecked", (), reason)
    root = math.sqrt(E / Fy)
    out = []

    def add(key, label, case, lam_r, e7, text, kc=None):
        el = sec.elements[key]
        cls = NONSLENDER if el.ratio <= lam_r else SLENDER          # B4.1: slender if the ratio EXCEEDS lambda_r
        out.append(ElementClass(key, label, case, "B4.1a", el.ratio, None, lam_r, cls, el, e7, kc, text))

    if sec.shape == I_SHAPE:
        web = sec.elements["web"]
        if sec.fabrication == ROLLED:
            # Case 1: flanges of rolled I-shaped sections, b/t, lambda_r = 0.56 sqrt(E/Fy)
            add("flange", "flange", 1, 0.56 * root, E7_OTHER, "0.56 sqrt(E/Fy)")
        else:
            # Case 2: flanges of built-up I-shaped sections, lambda_r = 0.64 sqrt(kc E/Fy), footnote [a]
            kc = kc_factor(web.ratio)
            add("flange", "flange", 2, 0.64 * math.sqrt(kc * E / Fy), E7_OTHER, "0.64 sqrt(kc E/Fy)", kc)
        # Case 5: webs of doubly symmetric rolled and built-up I-shaped sections, h/tw, 1.49 sqrt(E/Fy)
        add("web", "web", 5, 1.49 * root, E7_STIFFENED, "1.49 sqrt(E/Fy)")
    elif sec.shape == RECT_HSS:
        # Case 6: walls of rectangular HSS, b/t, lambda_r = 1.40 sqrt(E/Fy)
        add("wall_y", "wall (depth)", 6, 1.40 * root, E7_HSS_WALL, "1.40 sqrt(E/Fy)")
        add("wall_z", "wall (width)", 6, 1.40 * root, E7_HSS_WALL, "1.40 sqrt(E/Fy)")
    elif sec.shape == BOX:
        # Case 8: all other stiffened elements, b/t, lambda_r = 1.49 sqrt(E/Fy)
        add("wall_y", "wall (depth)", 8, 1.49 * root, E7_STIFFENED, "1.49 sqrt(E/Fy)")
        add("wall_z", "wall (width)", 8, 1.49 * root, E7_STIFFENED, "1.49 sqrt(E/Fy)")
    else:
        # Case 9: round HSS, D/t, lambda_r = 0.11 E/Fy (no square root)
        add("wall", "wall", 9, 0.11 * E / Fy, None, "0.11 E/Fy")
    return Classification("ok", tuple(out))


# ============================================================================== Table B4.1b
def classify_flexure(sec: DesignSection, E: float | None, Fy: float | None, axis: str = MAJOR) -> Classification:
    """Table B4.1b: compact / noncompact / slender elements of a member in flexure.

    ``axis``: ``"major"`` = bending about the AISC x-axis (local z); ``"minor"`` = about
    the AISC y-axis (local y).  For a rectangular section the "flanges" are the walls
    perpendicular to the plane of bending and the "webs" the walls parallel to it:
    major-axis bending -> flanges = width walls (``wall_z``), webs = depth walls (``wall_y``).
    """
    if axis not in (MAJOR, MINOR):
        raise ValueError("axis must be 'major' or 'minor'")
    if sec.shape == I_SHAPE:
        # the web ratio is needed for the web (major axis) and for kc of a built-up flange (Case 11)
        keys = ("flange", "web") if axis == MAJOR else ("flange",)
    else:
        keys = {RECT_HSS: ("wall_y", "wall_z"), BOX: ("wall_y", "wall_z"), ROUND_HSS: ("wall",)}.get(sec.shape, ())
    reason = _unavailable(sec, E, Fy, keys)
    if reason:
        return Classification("notChecked", (), reason)
    root = math.sqrt(E / Fy)
    out = []

    def add(key, label, case, lam_p, lam_r, text, kc=None):
        el = sec.elements[key]
        # B4.1: compact if no element exceeds lambda_p; noncompact if one exceeds lambda_p but none lambda_r
        cls = COMPACT if el.ratio <= lam_p else (NONCOMPACT if el.ratio <= lam_r else SLENDER)
        out.append(ElementClass(key, label, case, "B4.1b", el.ratio, lam_p, lam_r, cls, el, None, kc, text))

    if sec.shape == I_SHAPE:
        if axis == MINOR:
            # Case 13: flanges of all I-shaped sections in flexure about the minor axis
            add("flange", "flange", 13, 0.38 * root, 1.0 * root, "0.38 / 1.0 sqrt(E/Fy)")
        elif sec.fabrication == ROLLED:
            # Case 10: flanges of rolled I-shaped sections
            add("flange", "flange", 10, 0.38 * root, 1.0 * root, "0.38 / 1.0 sqrt(E/Fy)")
        else:
            # Case 11: flanges of doubly symmetric built-up I-shapes: lambda_r = 0.95 sqrt(kc E/FL).
            # Footnote [b]: FL = 0.7 Fy for slender-web members and for compact/noncompact-web members with
            # Sxt/Sxc >= 0.7; a doubly symmetric section has Sxt/Sxc = 1, so FL = 0.7 Fy in every case.
            kc = kc_factor(sec.elements["web"].ratio)
            add("flange", "flange", 11, 0.38 * root, 0.95 * math.sqrt(kc * E / (0.7 * Fy)),
                "0.38 sqrt(E/Fy) / 0.95 sqrt(kc E/FL), FL = 0.7 Fy", kc)
        if axis == MAJOR:
            # Case 15: webs of doubly symmetric I-shaped sections, h/tw
            add("web", "web", 15, 3.76 * root, 5.70 * root, "3.76 / 5.70 sqrt(E/Fy)")
    elif sec.shape in (RECT_HSS, BOX):
        flange, web = ("wall_z", "wall_y") if axis == MAJOR else ("wall_y", "wall_z")
        names = {"wall_y": "depth wall", "wall_z": "width wall"}
        if sec.shape == RECT_HSS:
            # Case 17: flanges of rectangular HSS
            add(flange, f"flange ({names[flange]})", 17, 1.12 * root, 1.40 * root, "1.12 / 1.40 sqrt(E/Fy)")
        else:
            # Case 21: flanges of box sections
            add(flange, f"flange ({names[flange]})", 21, 1.12 * root, 1.49 * root, "1.12 / 1.49 sqrt(E/Fy)")
        # Case 19: webs of rectangular HSS and box sections, h/t
        add(web, f"web ({names[web]})", 19, 2.42 * root, 5.70 * root, "2.42 / 5.70 sqrt(E/Fy)")
    else:
        # Case 20: round HSS, D/t: lambda_p = 0.07 E/Fy, lambda_r = 0.31 E/Fy
        add("wall", "wall", 20, 0.07 * E / Fy, 0.31 * E / Fy, "0.07 / 0.31 E/Fy")
    return Classification("ok", tuple(out))
