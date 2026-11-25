#!/usr/bin/env python3
"""Writes the problem directories of verification/problems/aisc-de/.

Source: AISC, Companion to the AISC Steel Construction Manual, Volume 1: Design Examples, Version 16.0,
Examples C.1A, C.1B and C.1C (pp. C-2 to C-17): a single-storey moment frame (two W12x65 columns with
pinned bases, 20 ft high, and a W18x40 beam of 30 ft) that stabilises a "leaning" column, under the
maximum gravity combination with a notional (C.1A, C.1B) or additional lateral (C.1C) load.  The
examples give frame-level results of an analysis by unnamed software in their figures: first-order and
second-order drift, column top moments, base shears and vertical reactions, for LRFD and ASD.

Nothing of the source is reproduced; input data and the printed results are given with their location.

What the source gives and what it does not (the reason for the status of these problems, see README.md):

given      geometry, loads, E (29,000 ksi is used in the examples' own evaluations), W12x65: Ag = 19.1 in^2
           and Ix = 533 in^4, W18x40: Ix = 612 in^4, the stiffness factor (0.8 in C.1A, 1.0 in C.1B/C.1C),
           the statement that the software accounts for P-Delta and P-delta;
NOT given  the area of the W18x40, the properties of the leaning column and of its link, whether the
           software includes axial and shear deformation, the shear areas.

Model here: axial deformation of the columns with the given Ag; beam, link and leaning column axially
rigid (area 1e5 in^2: the beam's axial force is about 6 kips, its shortening with this area 1e-6 in); no
shear deformation.  With an axially rigid beam the drifts at grid lines B and C coincide, so the
source's "average of the drifts at B and C" is the drift at B.  No property is taken from a shape table.

ASD: the source analyses at 1.6 times the ASD loads and prints reactions and moments divided by 1.6 and
the drift undivided; the targets are 1.6 x the printed reactions and moments (C.1A, C.1B).  In C.1C the
ASD loads are not multiplied.

Run:  python3 verification/problems/aisc-de/build_problems.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from v2b_common import Library, P, disp, model, reac, station, sup  # noqa: E402
from quarantine import QUARANTINE  # noqa: E402

WORK = ("American Institute of Steel Construction, Companion to the AISC Steel Construction Manual, Volume 1: "
        "Design Examples")
NOTES = ("input data and targets read on the cited pages and figures; nothing of the source is reproduced. The "
         "source's results are output of unnamed analysis software; the area of the W18x40, the leaning column and "
         "link properties and the treatment of axial and shear deformation are not stated (see analysis)")
H, SPAN, E = 240.0, 360.0, 29000.0
A_COL, I_COL, I_BEAM, A_RIGID = 19.1, 533.0, 612.0, 1.0e5

# example -> case -> data.  w kip/ft, P kips (column loads), lean kips, N kip (lateral);
# results: drift in., (moment kip-ft, base shear kips, vertical kips) for column B and column C
DATA = {
    "c1a": {"title": "Example C.1A (direct analysis method)", "pages": "pp. C-2 to C-6", "k": 0.8, "lean": True,
            "lrfd": {"w": "2.40", "P": "36.0", "lean": "144", "N": "0.576", "div": 1.0,
                     "first": ("0.181", ("113", "5.64", "71.6"), ("124", "6.21", "72.4")),
                     "second": ("0.290", ("109", "5.52", "71.4"), ("127", "6.26", "72.6"))},
            "asd": {"w": "2.56", "P": "38.4", "lean": "154", "N": "0.614", "div": 1.6,
                    "first": ("0.193", ("75.2", "3.76", "47.7"), ("82.8", "4.14", "48.3")),
                    "second": ("0.321", ("72.2", "3.68", "47.6"), ("84.8", "4.18", "48.4"))},
            "figs": {"first": "Fig. C.1A-3", "second": "Fig. C.1A-4"}},
    "c1b": {"title": "Example C.1B (effective length method)", "pages": "pp. C-7 to C-12", "k": 1.0, "lean": True,
            "lrfd": {"w": "2.40", "P": "36.0", "lean": "144", "N": "0.576", "div": 1.0,
                     "first": ("0.145", ("113", "5.64", "71.6"), ("124", "6.21", "72.4")),
                     "second": ("0.204", ("110", "5.56", "71.5"), ("126", "6.24", "72.5"))},
            "asd": {"w": "2.56", "P": "38.4", "lean": "154", "N": "0.614", "div": 1.6,
                    "first": ("0.155", ("75.2", "3.76", "47.7"), ("82.8", "4.14", "48.3")),
                    "second": ("0.223", ("73.1", "3.70", "47.6"), ("84.1", "4.16", "48.4"))},
            "figs": {"first": "Fig. C.1B-3", "second": "Fig. C-1B-4 (as captioned)"}},
    "c1c": {"title": "Example C.1C (first-order analysis method)", "pages": "pp. C-13 to C-17", "k": 1.0,
            "lean": False,
            "lrfd": {"w": "2.40", "P": "36.0", "N": "1.21", "div": 1.0,
                     "first": ("0.304", ("106", "5.32", "71.2"), ("131", "6.53", "72.8"))},
            "asd": {"w": "1.60", "P": "24.0", "N": "0.806", "div": 1.0,
                    "first": ("0.203", ("70.9", "3.55", "47.5"), ("87.1", "4.35", "48.5"))},
            "figs": {"first": "Fig. C.1C-3"}},
}


def frame_model(ex, case, order):
    d = DATA[ex][case]
    k = DATA[ex]["k"]
    nodes = {"B0": (0, 0), "B1": (0, H), "C0": (SPAN, 0), "C1": (SPAN, H)}
    plane = ("uz", "rx", "ry")
    supports = {"B0": sup("ux", "uy", *plane), "C0": sup("ux", "uy", *plane), "B1": sup(*plane), "C1": sup(*plane)}
    sections = {"col": {"type": "properties", "description": "W12x65: Ag and Ix as given by the source; Iy = 1000 Ix "
                                                             "(plane problem), J placeholder",
                        "A": A_COL, "Iz": I_COL, "Iy": 1000.0 * I_COL, "J": 1.0},
                "beam": {"type": "properties", "description": "W18x40: Ix as given by the source; area not given: "
                                                              "axially rigid; Iy = 1000 Ix, J placeholder",
                         "A": A_RIGID, "Iz": I_BEAM, "Iy": 1000.0 * I_BEAM, "J": 1.0}}
    members = {"colB": {"i": "B0", "j": "B1", "section": "col", "material": "M"},
               "colC": {"i": "C0", "j": "C1", "section": "col", "material": "M"},
               "beam": {"i": "B1", "j": "C1", "section": "beam", "material": "M"}}
    loads = [{"type": "memberDistributed", "member": "beam", "direction": "Y", "w1": -float(d["w"]) / 12.0},
             {"type": "nodal", "node": "B1", "Fy": -float(d["P"]), "Fx": float(d["N"])},
             {"type": "nodal", "node": "C1", "Fy": -float(d["P"])}]
    if DATA[ex]["lean"]:
        nodes.update({"L0": (2 * SPAN, 0), "L1": (2 * SPAN, H)})
        supports.update({"L0": sup("ux", "uy", "uz"), "L1": sup("uz")})
        sections["rigid"] = {"type": "properties", "description": "leaning column and link: properties not given by "
                                                                  "the source; axially rigid pin-ended bars",
                             "A": A_RIGID, "Iz": 1.0e5, "Iy": 1.0e5, "J": 1.0}
        members["lean"] = {"type": "truss", "i": "L0", "j": "L1", "section": "rigid", "material": "M"}
        members["link"] = {"type": "truss", "i": "C1", "j": "L1", "section": "rigid", "material": "M"}
        loads.append({"type": "nodal", "node": "L1", "Fy": -float(d["lean"])})
    return model("", ("kip", "in"), "Y", nodes, supports,
                 {"M": {"E": k * E, "nu": 0.3,
                        "description": f"{k:g} x 29,000 ksi" + (" (AISC 360 C2.3 stiffness reduction, tau_b = 1)"
                                                               if k != 1.0 else "")}},
                 sections, members, {"L": {"loads": loads}},
                 {"C": {"factors": {"L": 1.0}, "analysisType": "firstOrder" if order == "first" else "secondOrder"}},
                 analysis={"shearDeformation": False, "resultStations": 3})


def targets(ex, case, order):
    d = DATA[ex][case]
    drift, colb, colc = d[order]
    f = d["div"]
    asd = " (ASD: the source prints the value divided by 1.6; the analysis is at 1.6 x the ASD loads)" if f != 1 else ""
    out = [P(drift, disp("B1", "ux"), "in", note="drift = average of the drifts at grid lines B and C in the source; "
                                                  "equal at B and C here (axially rigid beam)"
             + (" (ASD: printed 'prior to dividing by 1.6')" if f != 1 else ""),
             signnote="the lateral load acts along +X; drift +X")]
    for name, node, (mom, shear, vert), ssign in (("colB", "B0", colb, 1.0), ("colC", "C0", colc, -1.0)):
        out.append(P(mom, station(name, "Mz", rel=1.0), "kip in", factor=12.0 * f, sign=ssign,
                     note="column moment at the beam, printed in kip-ft" + asd,
                     signnote="magnitude in the source's figure; in the engine's convention (column local x upward, "
                              "local y = global +X, local z = -Z) the top moment Mz is positive in column B and "
                              "negative in column C"))
        out.append(P(shear, reac(node, "Fx"), "kip", factor=f, sign=ssign, note="horizontal base reaction" + asd,
                     signnote="the source's arrows: toward +X at column B, toward -X at column C (force by the "
                              "support on the structure)"))
        out.append(P(vert, reac(node, "Fy"), "kip", factor=f, note="vertical base reaction" + asd,
                     signnote="upward"))
    return out


def build(lib):
    for ex, X in DATA.items():
        for case in ("lrfd", "asd"):
            d = X[case]
            for order in ("first", "second"):
                if order not in d:
                    continue
                pid = f"{ex}-{order}-order-{case}"
                loadtxt = (f"beam load {d['w']} kip/ft, {d['P']} kips on each frame column"
                           + (f", {d['lean']} kips on the leaning column" if X["lean"] else "")
                           + f", lateral load {d['N']} kip at the top of column B toward column C")
                analysis = (
                    f"{X['title']}, {case.upper()}: columns B and C W12x65 (Ag = 19.1 in^2, Ix = 533 in^4, given in "
                    "Examples C.1A-C), 20 ft, pinned bases; beam W18x40 (Ix = 612 in^4, given in C.1C), 30 ft, rigidly "
                    "connected; " + ("a pin-ended leaning column linked to column C by a pin-ended bar at roof level "
                                     "(30 ft to the right; the distance does not matter for rigid bars); "
                                     if X["lean"] else "no leaning column in the model of this example; ")
                    + f"E = 29,000 ksi x {X['k']:g}. Loads: {loadtxt}. "
                    + ("First-order analysis. " if order == "first" else
                       "Second-order analysis (the engine's members carry the exact stability functions: P-Delta and "
                       "P-delta). ")
                    + "Plane frame in the global X-Y plane (Y vertical), out-of-plane DOFs held, Iy = 1000 Ix so that "
                      "no out-of-plane buckling enters (the source's analysis is plane). NOT given by the source and "
                      "chosen here: beam, link and leaning column axially rigid (A = 1e5 in^2); no shear deformation. "
                      "The source's results are output of unnamed software.")
                lib.add(pid, f"AISC Design Example {ex[:2].upper()}.{ex[2:].upper()}: {order}-order frame results, "
                             f"{case.upper()}",
                        {"work": WORK, "edition": "Version 16.0 (2023)",
                         "location": f"{X['title'].split(' (')[0]}, {X['figs'][order]}, {case.upper()} column",
                         "pages": X["pages"], "notes": NOTES},
                        analysis, frame_model(ex, case, order), targets(ex, case, order),
                        diagnosis=QUARANTINE.get(pid),
                        limits="The source does not state the beam area, the leaning column / link properties or the "
                               "treatment of axial and shear deformation in its software; see README.md for the "
                               "measured influence of the choices made here.")
    unknown = set(QUARANTINE) - set(lib.problems)
    if unknown:
        raise SystemExit("quarantine.py names unknown problems: " + ", ".join(sorted(unknown)))


def main():
    lib = Library(HERE, "aisc-de")
    build(lib)
    lib.write()


if __name__ == "__main__":
    main()
