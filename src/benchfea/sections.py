"""Section properties of parametric shapes, computed from mechanics.

Axes: y = depth direction, z = width direction, origin at the centroid (all
shapes here are doubly symmetric).  ``Iz = int y^2 dA`` (bending about z, the
major axis of an I-shape, AISC "x-x"); ``Iy = int z^2 dA`` (AISC "y-y").

What is exact and what is an approximation (also listed for LIMITS.md):

=========== ============================== =======================================
quantity    exact (area integrals)         approximation / definition
=========== ============================== =======================================
A, I, S, Z  all shapes, fillets included
J           pipe (2I), rectangle (series)  I: sum(b t^3)/3 on the centre-line model
                                           (flanges bf x tf, web (d - tf) x tw,
                                           fillets and junction effects ignored);
                                           box: Bredt single-cell formula on the
                                           wall centre line
Cw          --                             I: tf bf^3 (d - tf)^2 / 24 (thin-walled,
                                           web and fillets ignored; identical to
                                           Iy ho^2/4 with Iy of the two flanges,
                                           cf. AISC 360-22 F2 User Note p. 16.1-54);
                                           box, pipe, rectangle: taken as 0
Asy, Asz    rectangle 5A/6 (energy,        I: Asy = d tw, Asz = (5/6) 2 bf tf;
            parabolic shear stress)        box: Asy = 2 H t, Asz = 2 B t;
                                           pipe: A/2 (thin-walled energy value)
=========== ============================== =======================================

Shear areas are *simple definitions*, not Cowper (1966) coefficients; give
``overrides.Asy/Asz`` to use a value from a refined theory.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace


@dataclass(frozen=True)
class SectionProps:
    """Section properties in whatever consistent length unit the inputs used."""

    family: str                      # "I" | "box" | "pipe" | "rectangle" | "general" | DB type
    A: float
    Iy: float
    Iz: float
    J: float
    Asy: float | None = None         # None: no shear deformation along y
    Asz: float | None = None
    Cw: float = 0.0
    Sy: float | None = None
    Sz: float | None = None
    Zy: float | None = None
    Zz: float | None = None
    dims: dict = field(default_factory=dict)       # design dimensions (d, bf, tf, tw, r / H, B, t / D, t / h, b)
    approximations: tuple[str, ...] = ()           # which values are approximate
    overridden: tuple[str, ...] = ()               # names replaced by `overrides`

    @property
    def ry(self) -> float:
        return math.sqrt(self.Iy / self.A)

    @property
    def rz(self) -> float:
        return math.sqrt(self.Iz / self.A)

    def with_overrides(self, ov: dict | None) -> "SectionProps":
        if not ov:
            return self
        return replace(self, **ov, overridden=tuple(ov))


def i_shape(d: float, bf: float, tf: float, tw: float, r: float = 0.0) -> SectionProps:
    """Doubly symmetric I-shape with rectangular flanges and optional fillets."""
    hw = d - 2.0 * tf                      # clear web height
    # --- rectangles: two flanges + web
    A = 2.0 * bf * tf + hw * tw
    Iz = (bf * d**3 - (bf - tw) * hw**3) / 12.0
    Iy = (2.0 * tf * bf**3 + hw * tw**3) / 12.0
    Zz = bf * tf * (d - tf) + tw * hw**2 / 4.0
    Zy = tf * bf**2 / 2.0 + hw * tw**2 / 4.0
    # --- four fillets: each is a square r x r minus a quarter disc of radius r.
    # Square: y in [hw/2 - r, hw/2], z in [tw/2, tw/2 + r]; disc centre (yc, zc).
    if r > 0.0:
        aq = math.pi * r * r / 4.0          # quarter-disc area
        e = 4.0 * r / (3.0 * math.pi)       # its centroid distance from the centre
        iq = math.pi * r**4 / 16.0          # its second moment about a bounding radius line
        yc, zc = hw / 2.0 - r, tw / 2.0 + r
        ys, zs = hw / 2.0 - r / 2.0, tw / 2.0 + r / 2.0       # square centroid
        A += 4.0 * (r * r - aq)
        # int y^2 dA: square (own + transfer) minus quarter disc (y = yc + eta, eta in [0, r])
        Iz += 4.0 * ((r**4 / 12.0 + r * r * ys**2) - (yc**2 * aq + 2.0 * yc * r**3 / 3.0 + iq))
        # int z^2 dA: quarter disc has z = zc - zeta, zeta in [0, r]
        Iy += 4.0 * ((r**4 / 12.0 + r * r * zs**2) - (zc**2 * aq - 2.0 * zc * r**3 / 3.0 + iq))
        Zz += 4.0 * (r * r * ys - aq * (yc + e))
        Zy += 4.0 * (r * r * zs - aq * (zc - e))
    h0 = d - tf                            # distance between flange centroids
    J = (2.0 * bf * tf**3 + h0 * tw**3) / 3.0
    Cw = tf * bf**3 * h0**2 / 24.0
    return SectionProps(
        "I", A, Iy, Iz, J, Asy=d * tw, Asz=(5.0 / 6.0) * 2.0 * bf * tf, Cw=Cw,
        Sy=Iy / (bf / 2.0), Sz=Iz / (d / 2.0), Zy=Zy, Zz=Zz,
        dims={"d": d, "bf": bf, "tf": tf, "tw": tw, "r": r},
        approximations=("J: sum(b t^3)/3, centre-line model, fillets ignored",
                        "Cw: tf bf^3 (d-tf)^2/24, thin-walled, web and fillets ignored",
                        "Asy = d tw; Asz = (5/6) 2 bf tf (simple definitions)"))


def box(H: float, B: float, t: float) -> SectionProps:
    """Rectangular hollow section, uniform wall, sharp corners."""
    hi, bi = H - 2.0 * t, B - 2.0 * t
    A = H * B - hi * bi
    Iz = (B * H**3 - bi * hi**3) / 12.0
    Iy = (H * B**3 - hi * bi**3) / 12.0
    hm, bm = H - t, B - t                  # wall centre-line dimensions
    J = 2.0 * t * hm**2 * bm**2 / (hm + bm)   # Bredt: 4 Am^2 / (perimeter / t)
    return SectionProps(
        "box", A, Iy, Iz, J, Asy=2.0 * H * t, Asz=2.0 * B * t, Cw=0.0,
        Sy=Iy / (B / 2.0), Sz=Iz / (H / 2.0),
        Zy=(H * B**2 - hi * bi**2) / 4.0, Zz=(B * H**2 - bi * hi**2) / 4.0,
        dims={"H": H, "B": B, "t": t},
        approximations=("J: Bredt single-cell formula on the wall centre line",
                        "Cw taken as 0", "Asy = 2 H t; Asz = 2 B t (simple definitions)",
                        "sharp corners (no corner radii)"))


def pipe(D: float, t: float) -> SectionProps:
    """Circular hollow section."""
    Di = D - 2.0 * t
    A = math.pi * (D**2 - Di**2) / 4.0
    I = math.pi * (D**4 - Di**4) / 64.0
    Z = (D**3 - Di**3) / 6.0
    return SectionProps(
        "pipe", A, I, I, 2.0 * I, Asy=A / 2.0, Asz=A / 2.0, Cw=0.0,
        Sy=I / (D / 2.0), Sz=I / (D / 2.0), Zy=Z, Zz=Z, dims={"D": D, "t": t},
        approximations=("Asy = Asz = A/2 (thin-walled energy value)",))


def rectangle_torsion_constant(h: float, b: float, terms: int = 200) -> float:
    """Saint-Venant torsion constant of a solid rectangle (series solution of
    the Prandtl stress function; verified against a finite-difference solution
    of the same Poisson problem in tests/test_sections.py)."""
    a, c = max(h, b), min(h, b)            # a = long side, c = short side
    s = 0.0
    for n in range(1, 2 * terms, 2):
        s += math.tanh(n * math.pi * a / (2.0 * c)) / n**5
    return a * c**3 / 3.0 * (1.0 - 192.0 * c / (math.pi**5 * a) * s)


def rectangle(h: float, b: float) -> SectionProps:
    """Solid rectangle, depth h along y, width b along z."""
    A = h * b
    Iz, Iy = b * h**3 / 12.0, h * b**3 / 12.0
    return SectionProps(
        "rectangle", A, Iy, Iz, rectangle_torsion_constant(h, b), Asy=5.0 * A / 6.0, Asz=5.0 * A / 6.0,
        Cw=0.0, Sy=h * b**2 / 6.0, Sz=b * h**2 / 6.0, Zy=h * b**2 / 4.0, Zz=b * h**2 / 4.0,
        dims={"h": h, "b": b},
        approximations=("Cw taken as 0", "Asy = Asz = 5A/6 (energy value for parabolic shear stress)"))


def check_dimensions(sec: dict) -> list[tuple[str, str]]:
    """Geometric consistency of a parametric section; returns (property, message)."""
    t = sec.get("type")
    out = []
    if t == "I":
        d, bf, tf, tw, r = sec["d"], sec["bf"], sec["tf"], sec["tw"], sec.get("r", 0.0)
        if not 2.0 * tf < d:
            out.append(("tf", f"two flanges (2*tf = {2 * tf:g}) do not fit in the depth d = {d:g}"))
        if not tw <= bf:
            out.append(("tw", f"web thickness tw = {tw:g} exceeds the flange width bf = {bf:g}"))
        if r > 0 and (r > (bf - tw) / 2.0 or r > (d - 2.0 * tf) / 2.0):
            out.append(("r", f"fillet radius r = {r:g} does not fit between web and flange tips"))
    elif t == "box":
        if not 2.0 * sec["t"] < min(sec["H"], sec["B"]):
            out.append(("t", f"wall thickness t = {sec['t']:g} leaves no hollow (need 2t < min(H, B))"))
    elif t == "pipe":
        if not 2.0 * sec["t"] < sec["D"]:
            out.append(("t", f"wall thickness t = {sec['t']:g} leaves no hollow (need 2t < D)"))
    return out


def from_definition(sec: dict, shapes_lookup=None) -> SectionProps:
    """Build SectionProps from a section definition of the input format (any
    consistent unit).  ``shapes_lookup(name)`` must return a SectionProps in the
    same unit and is needed only for ``type: standard``."""
    t = sec["type"]
    if t == "properties":
        p = SectionProps("general", sec["A"], sec["Iy"], sec["Iz"], sec["J"], sec.get("Asy"), sec.get("Asz"),
                         sec.get("Cw", 0.0), sec.get("Sy"), sec.get("Sz"), sec.get("Zy"), sec.get("Zz"))
        return p
    if t == "I":
        p = i_shape(sec["d"], sec["bf"], sec["tf"], sec["tw"], sec.get("r", 0.0))
    elif t == "box":
        p = box(sec["H"], sec["B"], sec["t"])
    elif t == "pipe":
        p = pipe(sec["D"], sec["t"])
    elif t == "rectangle":
        p = rectangle(sec["h"], sec["b"])
    elif t == "standard":
        if shapes_lookup is None:
            raise ValueError("a shapes lookup is required for standard sections")
        p = shapes_lookup(sec["shape"])
    else:
        raise ValueError(f"unknown section type {t!r}")
    return p.with_overrides(sec.get("overrides"))
