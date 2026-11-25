"""Reader of the two spreadsheets of the Ziemian & Ziemian benchmark-frame data set (Mendeley Data,
DOI 10.17632/39sjhchwtx.1, CC BY 4.0): "FEA Model Details_All frames.xlsx" (one tab per frame, the
text listing of the finite element model) and "Simulation Results_All frames.xlsx".

The spreadsheets themselves are not part of this repository; ``build_from_dataset.py`` reads them from
a directory given on the command line.  Nothing here interprets the data: the tabs are parsed into
plain dictionaries.
"""
from __future__ import annotations

import openpyxl


def _clean(v):
    return v.strip() if isinstance(v, str) else v


def _rows(ws):
    out = []
    for r in ws.iter_rows(values_only=True):
        r = [_clean(v) for v in r]
        while r and r[-1] in (None, ""):
            r.pop()
        out.append(r)
    return out


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def read_model(path, frame: int) -> dict:
    """The model tab "Frame <n>" as a dict: units, nodes (perfect / imperfect), fixities, elements,
    connections, sections, materials, loads."""
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    rows = _rows(wb[f"Frame {frame}"])
    m = {"frame": frame, "units": None, "perfect": {}, "imperfect": {}, "fixities": {}, "elements": {},
         "connections": {}, "sections": {}, "materials": {}, "nodalForces": {}, "nodalMoments": {},
         "uniformLoads": {}, "prescribed": {}, "other": []}
    block = None
    for r in rows:
        if not r:
            continue
        head = r[0] if isinstance(r[0], str) else ""
        if head.startswith("UNITS"):
            m["units"] = head.split(":", 1)[1].strip()
        key = None
        for token, name in (("Coordinates (PERFECT", "coords"), ("Fixities", "fix"), ("Connectivity & Attributes", "conn"),
                            ("Connections (cont.)", "cap"), ("Connections", "joints"), ("Part 1:", "s1"),
                            ("Part 2:", "s2"), ("Part 3:", "s3"), ("Part 4:", "s4"), ("(iv)  Material", "mat"),
                            ("Prescribed Displacements", "pd"), ("Prescribed Rotations", "pr"),
                            ("Nodal Forces", "nf"), ("Nodal Moments", "nm"), ("Uniform Loads", "ul"),
                            ("Temperature Effects", "te")):
            if head.startswith(token):
                key = name
                break
        if key:
            block = key
            continue
        if head.startswith("***") and "No " not in head:
            m["other"].append((block, head))
        if not _is_num(r[0]):
            continue
        k = int(r[0])
        if block == "coords":
            m["perfect"][k] = (float(r[1]), float(r[2]), float(r[3]))
            m["imperfect"][int(r[5])] = (float(r[6]), float(r[7]), float(r[8]))
        elif block == "fix":
            m["fixities"][k] = [str(v).upper() for v in r[1:7]]
        elif block == "conn":
            m["elements"][k] = {"i": int(r[1]), "j": int(r[2]), "beta": float(r[3]), "section": int(r[4]),
                                "material": int(r[5])}
        elif block == "joints":
            m["connections"][k] = [v for v in r[1:5]]            # kz_i, ky_i, kz_j, ky_j
        elif block == "s1":
            m["sections"].setdefault(k, {}).update(dict(zip(("A", "Izz", "Iyy", "J", "Cw"), r[1:6])))
        elif block == "s3":
            m["sections"].setdefault(k, {}).update(dict(zip(("Zzz", "Zyy", "Ayy", "Azz", "name"), r[1:6])))
        elif block == "mat":
            m["materials"][k] = dict(zip(("E", "nu", "Fy", "Wt", "name"), r[1:6]))
        elif block == "pd":
            m["prescribed"].setdefault(k, {})["u"] = r[1:4]
        elif block == "pr":
            m["prescribed"].setdefault(k, {})["r"] = r[1:4]
        elif block == "nf":
            m["nodalForces"][k] = tuple(float(v) for v in r[1:4])
        elif block == "nm":
            m["nodalMoments"][k] = tuple(float(v) for v in r[1:4])
        elif block == "ul":
            m["uniformLoads"][k] = tuple(float(v) for v in r[1:4])
        elif block in ("te", "cap", "s2", "s4"):
            if block == "te" or (block == "cap" and len(r) > 1):
                m["other"].append((block, r))
    return m


def read_results(path, frame: int) -> dict:
    """Results tab "Frame <n>": {"displacementUnit", "momentUnit", "joints": {node: {"WC": (right, left),
    "SIPC": ..., "LA": ...}}, "members": {number: {... , "Mp": Z Fy}}} with the stored (unrounded) numbers."""
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    rows = _rows(wb[f"Frame {frame}"])
    out = {"header": [rows[0], rows[1], rows[2]], "joints": {}, "members": {}}
    for r in rows[3:]:
        r = r + [None] * (16 - len(r))
        if _is_num(r[0]) and all(_is_num(v) for v in r[1:7]):
            out["joints"][int(r[0])] = {"WC": (r[1], r[2]), "SIPC": (r[3], r[4]), "LA": (r[5], r[6])}
        if _is_num(r[8]) and all(_is_num(v) for v in r[9:15]):
            out["members"][int(r[8])] = {"WC": (r[9], r[10]), "SIPC": (r[11], r[12]), "LA": (r[13], r[14]), "Mp": r[15]}
    return out


def read_alpha_cr(path) -> dict:
    """SUMMARY tab: {frame number: elastic critical load ratio} (single-valued rows only)."""
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    out = {}
    for r in _rows(wb["SUMMARY"]):
        r = r + [None] * (12 - len(r))
        for a, b in ((3, 6), (7, 10)):
            if _is_num(r[a]) and _is_num(r[b]):
                out[int(r[a])] = float(r[b])
    return out
