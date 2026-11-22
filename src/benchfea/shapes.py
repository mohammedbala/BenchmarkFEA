"""Reader for a USER-SUPPLIED AISC Shapes Database (v16.0), used only for
sections given by standard shape name.

The database is published by AISC and is not redistributable: this repository
contains none of its data.  Supply your own copy:

* ``$BENCHFEA_SHAPES_DB`` = path to the .xlsx or .csv file, or
* ``data/user/aisc-shapes-database-v16.0.xlsx`` (or ``.csv``) under the current
  working directory or the repository root.

Reading .xlsx needs ``openpyxl`` (``pip install benchfea[shapes]``).

ASSUMPTIONS (the AISC download page could not be opened while this reader was
written, so they are unverified and overridable):

* the first row of the first worksheet (or of the CSV) holds the column
  headers named in ``COLUMNS``; when a header occurs more than once the FIRST
  occurrence is used and is assumed to be the US customary value;
* lengths are inches (``DB_LENGTH_UNIT``).

Override the header names with ``$BENCHFEA_SHAPES_COLUMNS`` = path to a JSON
object ``{"label": [...], "Ix": "...", ...}`` and the length unit with
``$BENCHFEA_SHAPES_LENGTH_UNIT`` (``in`` or ``mm``).
"""
from __future__ import annotations

import csv
import json
import os
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from .sections import SectionProps
from .units import LENGTH_UNITS

ENV_DB = "BENCHFEA_SHAPES_DB"
ENV_COLUMNS = "BENCHFEA_SHAPES_COLUMNS"
ENV_UNIT = "BENCHFEA_SHAPES_LENGTH_UNIT"
DEFAULT_NAMES = ("aisc-shapes-database-v16.0.xlsx", "aisc-shapes-database-v16.0.csv")
DB_LENGTH_UNIT = "in"

#: logical name -> header(s) in the user's file
COLUMNS: dict[str, object] = {
    "label": ["AISC_Manual_Label", "EDI_Std_Nomenclature"],
    "type": "Type",
    "A": "A", "d": "d", "bf": "bf", "tw": "tw", "tf": "tf",
    "Ht": "Ht", "B": "B", "tdes": "tdes", "OD": "OD",
    "Ix": "Ix", "Zx": "Zx", "Sx": "Sx", "Iy": "Iy", "Zy": "Zy", "Sy": "Sy", "J": "J", "Cw": "Cw",
}
I_TYPES = ("W", "M", "S", "HP")

HINT = ("Shape-name lookup needs your own copy of the AISC Shapes Database v16.0 (not redistributable, so "
        f"not shipped). Put it at data/user/{DEFAULT_NAMES[0]} (or .csv) or set {ENV_DB} to its path; "
        "or define the section by `type: properties` or by parametric dimensions (I, box, pipe, rectangle).")


class ShapesDatabaseError(Exception):
    """Database absent, unreadable, or shape not found. ``hint`` says what to do."""

    def __init__(self, message: str, hint: str = HINT):
        super().__init__(message)
        self.hint = hint


def find_database() -> Path | None:
    env = os.environ.get(ENV_DB)
    if env:
        return Path(env)
    roots = [Path.cwd(), Path(__file__).resolve().parents[2]]
    for root in roots:
        for name in DEFAULT_NAMES:
            p = root / "data" / "user" / name
            if p.is_file():
                return p
    return None


def _rows(path: Path):
    if path.suffix.lower() == ".csv":
        with path.open(newline="", encoding="utf-8-sig") as f:
            yield from csv.reader(f)
        return
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        try:
            import openpyxl
        except ImportError:
            raise ShapesDatabaseError(
                "reading an .xlsx shapes database needs the optional package openpyxl",
                "pip install openpyxl (or benchfea[shapes]), or save the database sheet as .csv.") from None
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            labels = _columns()["label"]
            sheets = wb.worksheets
            chosen = sheets[0]
            for ws in sheets:      # the database sheet is the one whose header row has a label column
                first = next(ws.iter_rows(values_only=True), None)
                if first and any(h in first for h in labels):
                    chosen = ws
                    break
            yield from chosen.iter_rows(values_only=True)
        finally:
            wb.close()
        return
    raise ShapesDatabaseError(f"unsupported shapes database file type: {path.suffix}",
                              "Use the .xlsx file as downloaded from AISC, or a .csv export of its database sheet.")


def _columns() -> dict[str, list[str]]:
    cols = dict(COLUMNS)
    env = os.environ.get(ENV_COLUMNS)
    if env:
        cols.update(json.loads(Path(env).read_text(encoding="utf-8")))
    return {k: (v if isinstance(v, list) else [v]) for k, v in cols.items()}


@lru_cache(maxsize=4)
def _load(path_str: str, mtime: float, columns_key: str) -> dict[str, dict]:
    path = Path(path_str)
    it = _rows(path)
    header = next(it, None)
    if header is None:
        raise ShapesDatabaseError(f"shapes database {path} is empty")
    header = [str(h).strip() if h is not None else "" for h in header]
    cols = _columns()
    index: dict[str, int] = {}
    label_idx: list[int] = []
    for logical, names in cols.items():
        for name in names:
            if name in header:
                if logical == "label":
                    label_idx.append(header.index(name))
                else:
                    index.setdefault(logical, header.index(name))
    if not label_idx:
        raise ShapesDatabaseError(
            f"shapes database {path} has no label column (looked for {cols['label']})",
            f"Check that the first row holds the column headers, or map your headers with {ENV_COLUMNS}.")
    table: dict[str, dict] = {}
    for row in it:
        rec = {}
        for logical, k in index.items():
            v = row[k] if k < len(row) else None
            if logical != "type":
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    v = None
            rec[logical] = v
        for k in label_idx:
            name = row[k] if k < len(row) else None
            if name not in (None, ""):
                table.setdefault(str(name).strip().upper(), rec)
    return table


def lookup_raw(name: str) -> dict:
    """Raw record of a shape (database units). Raises ShapesDatabaseError."""
    path = find_database()
    if path is None:
        raise ShapesDatabaseError(f"no AISC shapes database found, cannot look up shape {name!r}")
    if not path.is_file():
        raise ShapesDatabaseError(f"shapes database file {path} does not exist")
    table = _load(str(path), path.stat().st_mtime, os.environ.get(ENV_COLUMNS, ""))
    rec = table.get(name.strip().upper())
    if rec is None:
        raise ShapesDatabaseError(
            f"shape {name!r} is not in the shapes database {path}",
            "Use the designation exactly as in the database label column (e.g. W14X90, HSS6X4X1/4).")
    return rec


def lookup(name: str, length_unit: str) -> SectionProps:
    """Section properties of a standard shape, converted (exactly) from the
    database length unit to ``length_unit``.

    Mapping to engine axes: AISC x-x -> local z (Iz, Sz, Zz); AISC y-y -> local y.
    Shear areas use the same simple definitions as the parametric shapes; for
    shape types without one they are None (no shear deformation).
    """
    rec = lookup_raw(name)
    db_unit = os.environ.get(ENV_UNIT, DB_LENGTH_UNIT)
    f = float(Fraction(LENGTH_UNITS[db_unit]) / Fraction(LENGTH_UNITS[length_unit]))

    def g(key, power):
        v = rec.get(key)
        return None if v is None else v * f**power

    missing = [k for k in ("A", "Ix", "Iy", "J") if rec.get(k) is None]
    if missing:
        raise ShapesDatabaseError(
            f"shape {name!r}: the database row has no numeric value for {missing}",
            f"Check the column mapping ({ENV_COLUMNS}) or define the section by properties.")
    typ = str(rec.get("type") or "").strip().upper()
    dims, asy, asz, family = {}, None, None, typ or "general"
    d, bf, tf, tw = g("d", 1), g("bf", 1), g("tf", 1), g("tw", 1)
    if typ in I_TYPES and None not in (d, bf, tf, tw):
        family, dims = "I", {"d": d, "bf": bf, "tf": tf, "tw": tw}
        asy, asz = d * tw, (5.0 / 6.0) * 2.0 * bf * tf
    elif typ == "HSS" and None not in (g("Ht", 1), g("B", 1), g("tdes", 1)):
        H, B, t = g("Ht", 1), g("B", 1), g("tdes", 1)
        family, dims = "box", {"H": H, "B": B, "t": t}
        asy, asz = 2.0 * H * t, 2.0 * B * t
    elif typ in ("HSS", "PIPE") and None not in (g("OD", 1), g("tdes", 1)):
        family, dims = "pipe", {"D": g("OD", 1), "t": g("tdes", 1)}
        asy = asz = g("A", 2) / 2.0
    return SectionProps(
        family, g("A", 2), g("Iy", 4), g("Ix", 4), g("J", 4), asy, asz, g("Cw", 6) or 0.0,
        Sy=g("Sy", 3), Sz=g("Sx", 3), Zy=g("Zy", 3), Zz=g("Zx", 3), dims=dims,
        approximations=("properties as tabulated in the user-supplied AISC Shapes Database",
                        "shear areas: simple definitions (see benchfea.sections)"))
