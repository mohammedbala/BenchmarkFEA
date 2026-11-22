"""Design view of a cross-section: what the member checks need beyond the
analysis properties -- the shape family, how it is made (rolled / built-up,
HSS / box) and the plate elements with their widths, thicknesses and
width-to-thickness ratios, each with its provenance.

Code-neutral data only; the limiting ratios live in the edition modules.
The width DEFINITIONS used to derive a width from dimensions are those of
AISC 360-22 Section B4.1 (cited per element); a ratio or width that is given
directly (e.g. a tabulated ``h/tw``) is used as given and reported as such.

AISC axis letters: x = major axis = bending about local z (``Ix = Iz``),
y = minor axis = bending about local y (docs/CONVENTIONS.md section 3; the
mapping is fixed and does not depend on which inertia is larger).

Scope: doubly symmetric I-shapes, rectangular HSS, rectangular box sections,
round HSS / pipe.  Anything else has ``shape = None`` and ``reason`` set; every
check then reports ``notChecked``.

How a section of the input file maps (``from_model``) -- the input format has
no field that says how a shape is made, so the mapping is by what the section
data itself says, and it is reported with every check:

================ ============================================================
section          design view
================ ============================================================
``standard``     database type W, M, S, HP -> rolled I-shape; web depth
                 h = d - 2 k with k = the ``kdes`` column when the user's
                 database has it (otherwise h is unknown and web-dependent
                 checks are notChecked).  HSS with Ht, B -> rectangular HSS
                 (t = ``tdes``); HSS/PIPE with OD -> round HSS (t = ``tdes``).
``I``            fillet radius r > 0 -> rolled I-shape, h = d - 2 tf - 2 r;
                 r = 0 (or absent) -> built-up welded I-shape, h = d - 2 tf.
``box``          box section (built-up, sharp corners): clear widths
                 H - 2 t, B - 2 t.  NOT treated as an HSS.
``pipe``         round HSS (D, t); t is used as the design wall thickness.
other            out of scope (``rectangle``, ``properties``, other database
                 types).
================ ============================================================
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

I_SHAPE, RECT_HSS, BOX, ROUND_HSS = "I", "rectHSS", "box", "roundHSS"
ROLLED, BUILT_UP = "rolled", "builtUp"


@dataclass(frozen=True)
class PlateElement:
    """One kind of plate element of the cross-section."""
    key: str                          # "flange" | "web" | "wall_y" | "wall_z" | "wall"
    ratio: float | None               # width-to-thickness ratio (lambda)
    width: float | None               # b, h (D for a round section); None when only the ratio is given
    thickness: float | None
    count: int                        # number of identical elements in the cross-section
    definition: str                   # how the width / ratio was obtained
    source: str                       # provenance of the ratio


@dataclass(frozen=True)
class DesignSection:
    shape: str | None                 # "I" | "rectHSS" | "box" | "roundHSS" | None (out of scope)
    id: str = ""
    reason: str = ""                  # why the section is out of scope / incomplete
    A: float | None = None
    Ix: float | None = None           # major axis (= Iz of the engine)
    Iy: float | None = None           # minor axis
    rx: float | None = None
    ry: float | None = None
    J: float | None = None
    Cw: float | None = None
    fabrication: str | None = None    # I-shapes: "rolled" | "builtUp"
    dims: dict = field(default_factory=dict)          # d, bf, tf, tw | H, B, t | D, t
    elements: dict = field(default_factory=dict)      # key -> PlateElement
    sources: dict = field(default_factory=dict)       # name -> provenance
    notes: tuple = ()
    # --- flexure (Chapter F): elastic / plastic section moduli; None = not available
    Sx: float | None = None           # major axis (= Sz of the engine)
    Zx: float | None = None
    Sy: float | None = None           # minor axis
    Zy: float | None = None

    # ------------------------------------------------------------------ reporting
    def source(self, name: str) -> str:
        return self.sources.get(name, self.sources.get("*", "given"))

    def inputs(self, *names: str) -> dict:
        """{name: (value, dimension, source)} for section values / dimensions."""
        dim = {"A": "length^2", "Ix": "length^4", "Iy": "length^4", "J": "length^4", "Cw": "length^6",
               "fabrication": "n/a", "shape": "n/a", "Sx": "length^3", "Zx": "length^3", "Sy": "length^3",
               "Zy": "length^3"}
        out = {}
        for n in names:
            v = self.dims[n] if n in self.dims else getattr(self, n)
            out[n if n not in ("shape", "fabrication") else f"section {n}"] = (v, dim.get(n, "length"),
                                                                              self.source(n))
        return out

    # ------------------------------------------------------------------ constructors
    @staticmethod
    def unsupported(reason: str, id: str = "", A: float | None = None, Ix: float | None = None,
                    Iy: float | None = None, source: str = "given") -> "DesignSection":
        """A section outside the shape scope.  Area and radii of gyration are kept because
        the tension checks (which need nothing else) apply to any cross-section."""
        return DesignSection(shape=None, id=id, reason=reason, A=A, Ix=Ix, Iy=Iy,
                             rx=_radius(A, Ix, None) if A else None, ry=_radius(A, Iy, None) if A else None,
                             sources={"*": source})

    @staticmethod
    def i_shape(fabrication: str, *, A: float, Ix: float | None = None, Iy: float | None = None,
                rx: float | None = None, ry: float | None = None, J: float | None = None,
                Cw: float | None = None, d: float | None = None, bf: float | None = None,
                tf: float | None = None, tw: float | None = None, h: float | None = None,
                k: float | None = None, r: float | None = None, flange_ratio: float | None = None,
                web_ratio: float | None = None, source: str = "given", id: str = "",
                sources: dict | None = None, notes: tuple = (), Sx: float | None = None,
                Zx: float | None = None, Sy: float | None = None, Zy: float | None = None,
                rts: float | None = None, ho: float | None = None) -> "DesignSection":
        """Doubly symmetric I-shape.

        ``Sx, Zx, Sy, Zy``: section moduli for flexure.  ``rts`` (AISC 360-22 Eq. F2-7) and
        ``ho`` (distance between flange centroids) may be given (tabulated values); they are
        kept in ``dims`` and take precedence over the values derived in Chapter F.

        ``fabrication``: ``"rolled"`` or ``"builtUp"`` (welded; fillet welds ignored).
        Web depth ``h`` (AISC 360-22 B4.1b): given ``h``; else for a rolled shape
        ``d - 2 k`` (k = outer flange face to web toe of fillet) or ``d - 2 tf - 2 r``
        (clear distance between flanges less the fillet radius at each flange,
        B4.1b(a)); for a built-up welded shape the clear distance between flanges
        ``d - 2 tf`` (B4.1b(b)).  ``flange_ratio`` (= bf / 2tf) and ``web_ratio``
        (= h / tw) may be given directly (tabulated values) and then take precedence.
        """
        if fabrication not in (ROLLED, BUILT_UP):
            raise ValueError("fabrication must be 'rolled' or 'builtUp'")
        src = {"*": source, **(sources or {})}
        Ix, Iy, rx, ry = _inertia(A, Ix, rx), _inertia(A, Iy, ry), _radius(A, Ix, rx), _radius(A, Iy, ry)
        # --- flange: B4.1a(a) width b = bf / 2
        b = bf / 2.0 if bf is not None else None
        if flange_ratio is not None:
            f = PlateElement("flange", flange_ratio, b, tf, 4, "bf/(2 tf) as given", src.get("flange_ratio", source))
        elif b is not None and tf is not None:
            f = PlateElement("flange", b / tf, b, tf, 4, "b = bf/2 (B4.1a(a)), t = tf", "derived: bf, tf")
        else:
            f = PlateElement("flange", None, b, tf, 4, "bf and tf (or bf/2tf) not given", "missing")
        # --- web: B4.1b(a) rolled, B4.1b(b) built-up welded
        hdef = ""
        if h is not None:
            hdef = "h as given"
        elif fabrication == ROLLED:
            if k is not None and d is not None:
                h, hdef = d - 2.0 * k, "h = d - 2 k (clear distance between flanges less the fillet, B4.1b(a))"
            elif r is not None and d is not None and tf is not None:
                h, hdef = d - 2.0 * tf - 2.0 * r, ("h = d - 2 tf - 2 r (clear distance between flanges less "
                                                  "the fillet radius at each flange, B4.1b(a))")
        elif d is not None and tf is not None:
            h, hdef = d - 2.0 * tf, "h = d - 2 tf (clear distance between flanges, welded, B4.1b(b))"
        if web_ratio is not None:
            w = PlateElement("web", web_ratio, h, tw, 1, "h/tw as given", src.get("web_ratio", source))
        elif h is not None and tw is not None:
            w = PlateElement("web", h / tw, h, tw, 1, hdef, "derived: " + ("h, tw" if hdef == "h as given"
                                                                          else "section dimensions"))
        else:
            w = PlateElement("web", None, h, tw, 1,
                             "web depth h not available: a rolled shape needs the fillet dimension (k or r), "
                             "h, or h/tw", "missing")
        dims = {n: v for n, v in (("d", d), ("bf", bf), ("tf", tf), ("tw", tw), ("h", h), ("rts", rts),
                                  ("ho", ho)) if v is not None}
        src.setdefault("fabrication", source)
        return DesignSection(I_SHAPE, id, "", A, Ix, Iy, rx, ry, J, Cw, fabrication, dims,
                             {"flange": f, "web": w}, src, tuple(notes), Sx, Zx, Sy, Zy)

    @staticmethod
    def rectangular(kind: str, *, A: float, H: float | None = None, B: float | None = None,
                    t: float | None = None, Ix: float | None = None, Iy: float | None = None,
                    rx: float | None = None, ry: float | None = None, J: float | None = None,
                    h: float | None = None, b: float | None = None, h_ratio: float | None = None,
                    b_ratio: float | None = None, source: str = "given", id: str = "",
                    sources: dict | None = None, notes: tuple = (), Sx: float | None = None,
                    Zx: float | None = None, Sy: float | None = None,
                    Zy: float | None = None) -> "DesignSection":
        """Rectangular HSS (``kind="rectHSS"``) or built-up box section (``kind="box"``).

        ``H`` = outside depth (along local y), ``B`` = outside width (along local z),
        ``t`` = DESIGN wall thickness (AISC 360-22 B4.2 -- an input; no factor is applied).
        Flat widths: given ``h``/``b``; else HSS: outside dimension - 3 t (corner radius not
        known, B4.1b(d)); box: clear distance, outside dimension - 2 t (B4.1b(e)).
        ``h_ratio`` (= h/t, walls along the depth) and ``b_ratio`` (= b/t) may be given.
        """
        if kind not in (RECT_HSS, BOX):
            raise ValueError("kind must be 'rectHSS' or 'box'")
        src = {"*": source, **(sources or {})}
        Ix, Iy, rx, ry = _inertia(A, Ix, rx), _inertia(A, Iy, ry), _radius(A, Ix, rx), _radius(A, Iy, ry)
        n, rule = (3.0, "outside dimension - 3 t (corner radius not known, B4.1b(d))") if kind == RECT_HSS \
            else (2.0, "clear distance between walls = outside dimension - 2 t (B4.1b(e))")

        def wall(key, width, outside, ratio, name):
            definition = f"{name} as given"
            if width is None and outside is not None and t is not None:
                width, definition = outside - n * t, f"{name} = {rule}"
            if ratio is not None:
                return PlateElement(key, ratio, width, t, 2, f"{name}/t as given", src.get(f"{name}_ratio", source))
            if width is not None and t is not None:
                return PlateElement(key, width / t, width, t, 2, definition, "derived: section dimensions")
            return PlateElement(key, None, width, t, 2, f"{name} and t (or {name}/t) not given", "missing")

        dims = {k: v for k, v in (("H", H), ("B", B), ("t", t)) if v is not None}
        return DesignSection(kind, id, "", A, Ix, Iy, rx, ry, J, 0.0, None, dims,
                             {"wall_y": wall("wall_y", h, H, h_ratio, "h"),
                              "wall_z": wall("wall_z", b, B, b_ratio, "b")}, src, tuple(notes), Sx, Zx, Sy, Zy)

    @staticmethod
    def round_hss(*, A: float, D: float | None = None, t: float | None = None, I: float | None = None,
                  r: float | None = None, J: float | None = None, ratio: float | None = None,
                  source: str = "given", id: str = "", sources: dict | None = None,
                  notes: tuple = (), S: float | None = None, Z: float | None = None) -> "DesignSection":
        """Round HSS / pipe. ``t`` = DESIGN wall thickness (B4.2, an input); ``ratio`` = D/t as given;
        ``S``, ``Z`` = elastic / plastic section modulus (any axis)."""
        src = {"*": source, **(sources or {})}
        I, r = _inertia(A, I, r), _radius(A, I, r)
        if ratio is not None:
            w = PlateElement("wall", ratio, D, t, 1, "D/t as given", src.get("ratio", source))
        elif D is not None and t is not None:
            w = PlateElement("wall", D / t, D, t, 1, "D = outside diameter, t = design wall thickness (B4.1b(g))",
                             "derived: section dimensions")
        else:
            w = PlateElement("wall", None, D, t, 1, "D and t (or D/t) not given", "missing")
        dims = {k: v for k, v in (("D", D), ("t", t)) if v is not None}
        return DesignSection(ROUND_HSS, id, "", A, I, I, r, r, J, 0.0, None, dims, {"wall": w}, src, tuple(notes),
                             S, Z, S, Z)

    # ------------------------------------------------------------------ from the model
    @staticmethod
    def from_model(section) -> "DesignSection":
        """Design view of a :class:`benchfea.model.Section` (see the table in the module docstring)."""
        sid, kind, fam, dm = section.id, section.kind, section.props.family, dict(section.props.dims)
        src = f"section:{sid}"
        common = dict(A=section.A, Ix=section.Iz, Iy=section.Iy, J=section.J, source=src, id=sid)
        moduli = dict(Sx=section.Sz, Zx=section.Zz, Sy=section.Sy, Zy=section.Zy)      # x = major = local z
        if kind == "I" or (kind == "standard" and fam == "I"):
            if kind == "standard":
                fab, why = ROLLED, f"derived: shapes database type of {section.shape} (rolled shape)"
                note = ()
                if "k" not in dm:
                    note = (f"section {sid!r}: the shapes database row gives no fillet dimension (column kdes), "
                            "so the web depth h (AISC 360-22 B4.1b(a)) is unknown",)
                extra = {"k": dm.get("k")}
            elif dm.get("r", 0.0) > 0.0:
                fab, why = ROLLED, "derived: the section has a fillet radius r > 0 (rolled shape)"
                note, extra = (), {"r": dm["r"]}
            else:
                fab, why = BUILT_UP, "derived: the section has no fillet radius (welded plates)"
                note = (f"section {sid!r} has no fillet radius: checked as a built-up (welded) I-shape "
                        "(AISC 360-22 Table B4.1a Case 2, G2.1(b)); give r > 0 for a rolled shape",)
                extra = {}
            return DesignSection.i_shape(fab, Cw=section.Cw, d=dm["d"], bf=dm["bf"], tf=dm["tf"], tw=dm["tw"],
                                         sources={"fabrication": why}, notes=note, **extra, **common, **moduli)
        if kind == "box":
            note = (f"section {sid!r} (type box, sharp corners) is checked as a built-up box section "
                    "(AISC 360-22 B4.1b(e), Table B4.1a Case 8, Table E7.1 case (a), G4 box), not as an HSS",)
            return DesignSection.rectangular(BOX, H=dm["H"], B=dm["B"], t=dm["t"], notes=note, **common, **moduli)
        if kind == "standard" and fam == "box":
            note = ("t is the design wall thickness of the shapes database (AISC 360-22 B4.2)",)
            return DesignSection.rectangular(RECT_HSS, H=dm["H"], B=dm["B"], t=dm["t"], notes=note, **common,
                                             **moduli)
        if kind == "pipe" or (kind == "standard" and fam == "pipe"):
            note = ("t is used as the design wall thickness (AISC 360-22 B4.2): the 0.93 factor for HSS made "
                    "to standards other than ASTM A1065 / A1085 is NOT applied by the engine",)
            common.pop("Ix"), common.pop("Iy")
            return DesignSection.round_hss(D=dm["D"], t=dm["t"], I=section.Iz, notes=note, S=section.Sz,
                                           Z=section.Zz, **common)
        what = f"database shape type {fam!r}" if kind == "standard" else f"section type {kind!r}"
        return DesignSection.unsupported(
            f"section {sid!r} ({what}) is outside the implemented scope (doubly symmetric I-shapes, "
            "rectangular HSS / box sections, round HSS / pipe)", sid, section.A, section.Iz, section.Iy, src)


def _inertia(A, I, r):
    return I if I is not None else (A * r * r if r is not None else None)


def _radius(A, I, r):
    return r if r is not None else (math.sqrt(I / A) if I is not None else None)
