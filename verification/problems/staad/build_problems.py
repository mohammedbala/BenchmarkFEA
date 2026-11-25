#!/usr/bin/env python3
"""Writes the problem directories of verification/problems/staad/.

Source of the problem STATEMENTS: Bentley Systems, STAAD.Pro Help, "Verification Examples" (public
documentation site docs.bentley.com, pages of the STAAD.Pro 2023 samples; read 2026-11-22).  Each page
states a problem, cites a textbook, prints a "theory" / "reference" value and the program's value, and
lists the program input.  Nothing of a page is reproduced here: the data that define the problem are
given in this project's own words, with the page title and identifier.

The cited textbooks were not opened, so the basis of every problem is
"recomputed-from-published-statement" (see ../vendor_common.py): the reference values are computed in
verification/derive/vendor_reference.py; the page's printed theory values are recorded as additional
targets with the printed-precision tolerance, and those that do not reproduce are quarantined
separately with a diagnosis (quarantine.py).  STAAD.Pro's own results are not used as targets.

Where a page's text and its input listing differ, the input listing defines the problem (noted in the
problem).  Sections taken by a page from a shape table without printing the properties make the problem
unusable unless the page prints the property that matters (README.md lists what was excluded).

Run:  python3 verification/problems/staad/build_problems.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from vendor_common import (FIXED, Library, P, R, add_vendor, bar, buck, disp, extreme, fix, frame, model, one,  # noqa: E402
                           plane, reac, reference, sec, station, truss, truss_supports)
from quarantine import QUARANTINE  # noqa: E402

BASE = "https://docs.bentley.com/LiveContent/web/STAAD.Pro%20Help-v21/en/GUID-"
PLANE = ("Plane problem in the global X-Y plane (Y vertical); the out-of-plane DOFs (uz, rx, ry) of every node are "
         "held; in-plane bending about local z. ")
MAG = "the page prints a magnitude; sign from the reference solution, engine convention (docs/CONVENTIONS.md)"


def source(page, cites):
    title, guid = page
    return {"work": "Bentley Systems, STAAD.Pro Help, Verification Examples (public documentation, "
                    "docs.bentley.com)",
            "edition": "STAAD.Pro Help v21 build (pages refer to the STAAD.Pro 2023 sample files), read 2026-11-22",
            "location": f"'V. {title}' ({BASE}{guid}.html)",
            "pages": "web page",
            "notes": f"problem statement taken from the vendor page (text and input listing). The page cites: {cites} "
                     "- NOT opened for this project; no number is taken from it. The page's 'theory' values are "
                     "recorded as printed; the program's own results are not used."}


def build(lib):
    def add(pid, *a, printed_bad=(), **kw):
        add_vendor(lib, source, pid, *a, printed_bad=printed_bad, diagnosis=QUARANTINE.get(pid + "-vendor-printed"),
                   whole_diagnosis=QUARANTINE.get(pid), **kw)

    # ------------------------------------------------------------------ beams
    r = reference("staad_axially_loaded_column")
    n = {"n1": (0, 0), "n2": (0, 4), "n3": (0, 7), "n4": (0, 10)}
    lc, cb = one([{"type": "nodal", "node": "n3", "Fy": -1000.0}, {"type": "nodal", "node": "n2", "Fy": -500.0}])
    add("axially-loaded-column", "Bar held at both ends with two axial loads: support reactions",
        ("Axially Loaded Column", "300512D3-2533-4502-8921-235D1482B4EE"),
        "Timoshenko, Strength of Materials, Part I, 3rd ed. 1956, p. 26, problem 10",
        PLANE + "Vertical bar of 10 in held against axial and lateral translation at both ends (nodes at 0, 4, 7, "
        "10 in), 1000 lb downward at 7 in and 500 lb downward at 4 in; A = 1 in^2, Iz = 0.0001 in^4 (the page's "
        "input), E = 30e6 psi. First-order analysis.",
        model("", ("lbf", "in"), "Y", n, plane(n, {"n1": fix("ux", "uy"), "n4": fix("ux", "uy")}),
              {"M": {"E": 30.0e6, "nu": 0.3}}, {"S": sec(1.0, 0.0001)},
              {"m1": frame("n1", "n2"), "m2": frame("n2", "n3"), "m3": frame("n3", "n4")}, lc, cb),
        [R(r["R_bottom"], reac("n1", "Fy"), "lbf", "ODE oracle, bar"), R(r["R_top"], reac("n4", "Fy"), "lbf", "ODE oracle, bar"),
         R(r["u4"], disp("n2", "uy"), "in", "ODE oracle, bar"), R(r["u7"], disp("n3", "uy"), "in", "ODE oracle, bar")],
        [P("600", reac("n1", "Fy"), "lbf", signnote="upward"), P("900", reac("n4", "Fy"), "lbf", signnote="upward")])

    r = reference("staad_beam_trapezoid")
    n = {"a": (0, 0), "b": (7.5, 0)}
    lc, cb = one([{"type": "memberDistributed", "member": "m", "direction": "Y", "w1": -4.0, "w2": -7.0,
                   "start": 0.4, "end": 1.0}])
    add("beam-partial-trapezoidal-load", "Pinned-fixed beam with a linearly varying load on part of the span",
        ("Deflection and Reactions in a Beam", "42CF0C95-4031-4EB9-8BEF-3E3BA17D9CBC"),
        "Roark's Formulas for Stress and Strain (W. C. Young), 6th ed., Table 3, case 2c, p. 103",
        PLANE + "Beam of 7.5 m, pinned at A (x = 0), built in at B; line load from 4 kN/m at x = 3 m to 7 kN/m at "
        "B (downward); Iz = 5000 cm^4, A = 50 cm^2, E = 200 kN/mm^2. Units kN, m. First-order analysis, no shear "
        "deformation.",
        model("", ("kN", "m"), "Y", n, plane(n, {"a": fix("ux", "uy"), "b": fix("ux", "uy", "rz")}),
              {"M": {"E": 2.0e8, "nu": 0.3}}, {"S": sec(50.0e-4, 5000.0e-8)}, {"m": frame("a", "b")}, lc, cb,
              analysis={"shearDeformation": False, "resultStations": 6}),
        [R(r["thetaA"], disp("a", "rz"), "rad", "ODE oracle"), R(r["RA"], reac("a", "Fy"), "kN", "ODE oracle"),
         R(r["RB"], reac("b", "Fy"), "kN", "ODE oracle"), R(r["MB_reaction"], reac("b", "Mz"), "kN m", "ODE oracle"),
         R(r["v3"], station("m", "uy", rel=0.4), "m", "ODE oracle")],
        [P("0.0020", disp("a", "rz"), "rad", sign=-1.0, signnote=MAG), P("3.2886", reac("a", "Fy"), "kN"),
         P("21.461", reac("b", "Fy"), "kN"), P("25.9605", reac("b", "Mz"), "kN m", sign=-1.0, signnote=MAG)])

    r = reference("staad_propped_cantilever_selfweight")
    n = {"a": (0, 0), "c": (42.15, 0), "b": (100, 0)}
    lc, cb = one([{"type": "selfWeight"}])
    add("propped-cantilever-self-weight", "Propped cantilever under its own weight",
        ("Forces on a Propped Cantilever 2", "A8D232D3-FBB3-4D80-BFD5-B5D25123CD7B"),
        "Roark and Young, Formulas for Stress and Strain, 5th ed., p. 109, no. 23",
        PLANE + "Beam of 100 in, simply supported at node 1 (x = 0), built in at x = 100 in; A = 4 in^2, Iz = "
        "1.3333 in^4, depth 2 in, E = 30e6 psi; load = self weight with the density of the page's input listing, "
        "0.282072 lb/in^3 (the page's text rounds it to 0.2821; the printed results correspond to the listing). A "
        "node at x = 42.15 in. First-order analysis, no shear deformation.",
        model("", ("lbf", "in"), "Y", n, plane(n, {"a": fix("ux", "uy"), "b": fix("ux", "uy", "rz")}),
              {"M": {"E": 30.0e6, "nu": 0.3, "unitWeight": 0.282072}}, {"S": sec(4.0, 1.3333)},
              {"m1": frame("a", "c"), "m2": frame("c", "b")}, lc, cb),
        [R(r["R1"], reac("a", "Fy"), "lbf", "ODE oracle"), R(r["R3"], reac("b", "Fy"), "lbf", "ODE oracle"),
         R(r["M3_reaction"], reac("b", "Mz"), "lbf in", "ODE oracle"), R(r["theta1"], disp("a", "rz"), "rad", "ODE oracle"),
         R(r["v_4215"], disp("c", "uy"), "in", "ODE oracle"), R(r["M_4215"], station("m1", "Mz", rel=1.0), "lbf in", "ODE oracle")],
        [P("42.31", reac("a", "Fy"), "lbf"), P("70.52", reac("b", "Fy"), "lbf"),
         P("1410.4", reac("b", "Mz"), "lbf in", sign=-1.0, signnote=MAG),
         P("585.9", station("m1", "Mz", rel=1.0), "lbf in", factor=1.3333,
           note="the page prints the bending stress at the node at 42.15 in (psi); moment = stress x Iz / (depth/2) = "
                "stress x 1.3333 in^3", signnote="sagging: positive Mz (local y up)"),
         P("0.000588", disp("a", "rz"), "rad", sign=-1.0, signnote="printed negative; same sense (rotation about +Z)"),
         P("0.01528", disp("c", "uy"), "in", sign=-1.0, signnote="printed negative; downward")])

    r = reference("staad_stepped_shaft")
    n = {"n1": (0, 0, 0), "n2": (500, 0, 0), "n3": (1000, 0, 0), "n4": (1500, 0, 0)}
    secs = {f"d{d}": sec(math.pi * d * d / 4.0, math.pi * d ** 4 / 64.0, j=math.pi * d ** 4 / 32.0,
                         description=f"solid circular section, diameter {d} mm") for d in (80, 60, 40)}
    lc, cb = one([{"type": "nodal", "node": "n2", "Mx": 3000.0}, {"type": "nodal", "node": "n3", "Mx": 2000.0},
                  {"type": "nodal", "node": "n4", "Mx": 800.0}])
    add("stepped-shaft-torsion", "Stepped cantilever shaft under three torques: angle of twist",
        ("Torsion on a Stepped Cantilever", "B19AED4F-3D2E-40A4-AA61-15630220B3F5"),
        "Gere and Timoshenko, Mechanics of Materials, 2nd ed., p. 171, problem 3.3-1",
        "Shaft along global X, built in at x = 0, three segments of 500 mm with solid circular sections of 80, 60 and "
        "40 mm diameter (J = pi d^4 / 32), G = 80 kN/mm^2 (E = 200 kN/mm^2, nu = 0.25). Torques about the axis at the "
        "ends of the segments: 3000, 2000 and 800 kN mm - the units of the page's input listing and of its own "
        "arithmetic (0.0090 + 0.0138 + 0.0199 rad); the page's text says N mm, which is inconsistent with its "
        "result. Units kN, mm. First-order analysis (uniform torsion).",
        model("", ("kN", "mm"), "Z", n, {"n1": FIXED}, {"M": {"E": 200.0, "G": 80.0}}, secs,
              {"m1": frame("n1", "n2", "d80"), "m2": frame("n2", "n3", "d60"), "m3": frame("n3", "n4", "d40")}, lc, cb),
        [R(r["twist2"], disp("n2", "rx"), "rad", "ODE oracle, bar in torsion per segment, segment torques by statics"),
         R(r["twist3"], disp("n3", "rx"), "rad", "ODE oracle, bar in torsion per segment"),
         R(r["twist4"], disp("n4", "rx"), "rad", "ODE oracle, bar in torsion per segment")],
        [P("0.0427", disp("n4", "rx"), "rad", signnote="twist in the sense of the torques (about +X)")])

    # ------------------------------------------------------------------ frames
    r = reference("staad_portal_reactions_1")
    n = {"n1": (0, 0), "n2": (0, 100), "n3": (60, 100), "n4": (120, 100), "n5": (120, 0)}
    lc, cb = one([{"type": "nodal", "node": "n3", "Fy": -1000.0}])
    add("portal-pinned-midspan-load", "Pinned-base portal frame with a load at midspan of the beam",
        ("2D Portal Reactions 1", "1EAB2FE4-E21F-462E-A4BC-DDFA86663722"),
        "Timoshenko, Strength of Materials, Part I, 2nd ed. 1940, pp. 188-191",
        PLANE + "Portal: columns 100 in high, 2 in x 2 in (A = 4 in^2, I = 1.3333 in^4), beam of 120 in, 2 in wide "
        "and 4 in deep (A = 8 in^2, I = 10.667 in^4), E = 30e6 psi, pinned bases; 1000 lb at midspan (taken "
        "downward; the page's listing applies it upward, its results are magnitudes). First-order analysis with "
        "axial deformation, without shear deformation.",
        model("", ("lbf", "in"), "Y", n, plane(n, {"n1": fix("ux", "uy"), "n5": fix("ux", "uy")}),
              {"M": {"E": 30.0e6, "nu": 0.3}},
              {"col": sec(4.0, 2.0 ** 4 / 12.0), "beam": sec(8.0, 2.0 * 4.0 ** 3 / 12.0)},
              {"c1": frame("n1", "n2", "col"), "b1": frame("n2", "n3", "beam"), "b2": frame("n3", "n4", "beam"),
               "c2": frame("n5", "n4", "col")}, lc, cb),
        [R(r["R"]["n1"][0], reac("n1", "Fx"), "lbf", "plane frame program"),
         R(r["R"]["n1"][1], reac("n1", "Fy"), "lbf", "plane frame program"),
         R(r["R"]["n5"][0], reac("n5", "Fx"), "lbf", "plane frame program"),
         R(r["u"]["n3"][1], disp("n3", "uy"), "in", "plane frame program"),
         R(-r["end"]["b1"][2], station("b1", "Mz", rel=0.0), "lbf in", "plane frame program")],
        [P("500", reac("n1", "Fy"), "lbf"), P("27.55", reac("n1", "Fx"), "lbf", signnote=MAG)],
        printed_bad=[P("2754.97", station("b1", "Mz", rel=0.0), "lbf in", sign=-1.0, signnote=MAG,
                       note="moment at the ends of the beam; value of the page's closed-form expression")])

    r = reference("staad_portal_reactions_2")
    n = {"n1": (0, 0), "n2": (0, 20), "n3": (20, 20), "n4": (20, 0)}
    lc, cb = one([{"type": "memberDistributed", "member": "b", "direction": "Y", "w1": -2.0}])
    add("portal-fixed-uniform-load", "Fixed-base square portal frame with a uniform load on the beam",
        ("2D Portal Reactions 2", "A208EC6A-C93D-44C7-A0A7-DFF5ED10E8F9"),
        "McCormac, Structural Analysis, 3rd ed. 1975, p. 383, example 22-5",
        PLANE + "Portal 20 ft x 20 ft, fixed bases, all members A = 1 ft^2, I = 0.05 ft^4, E = 4.132e6 kip/ft^2 "
        "(the page's listing); 2 kip/ft downward on the beam. Units kip, ft. First-order analysis with axial "
        "deformation. The page's 'maximum moment' (44.40 kip ft) is the largest member END moment, at the beam-column "
        "joints; the midspan moment of the beam is larger (55.6 kip ft, also a target here).",
        model("", ("kip", "ft"), "Y", n, plane(n, {"n1": fix("ux", "uy", "rz"), "n4": fix("ux", "uy", "rz")}),
              {"M": {"E": 4.132e6, "nu": 0.3}}, {"S": sec(1.0, 0.05)},
              {"c1": frame("n1", "n2"), "b": frame("n2", "n3"), "c2": frame("n4", "n3")}, lc, cb,
              analysis={"shearDeformation": False, "resultStations": 3}),
        [R(r["R"]["n1"][0], reac("n1", "Fx"), "kip", "plane frame program"),
         R(r["R"]["n1"][2], reac("n1", "Mz"), "kip ft", "plane frame program"),
         R(r["u"]["n2"][2], disp("n2", "rz"), "rad", "plane frame program"),
         R(-r["end"]["b1"][2], station("b", "Mz", rel=0.0), "kip ft", "plane frame program"),
         R(r["end"]["b1"][5], station("b", "Mz", rel=0.5), "kip ft", "plane frame program")],
        printed_bad=[P("44.40", station("b", "Mz", rel=0.0), "kip ft", sign=-1.0, signnote=MAG,
                       note="the page's 'maximum moment': member end moment at the beam-column joint")])

    r = reference("staad_frame_1x2_lateral")
    n = {"n1": (0, 0), "n2": (20, 0), "n3": (0, 20), "n4": (20, 20), "n5": (0, 40), "n6": (20, 40)}
    lc, cb = one([{"type": "nodal", "node": "n3", "Fx": 20.0}, {"type": "nodal", "node": "n5", "Fx": 10.0}])
    mm = {k: frame(i, j) for k, i, j in (("m1", "n1", "n3"), ("m2", "n2", "n4"), ("m3", "n3", "n5"),
                                         ("m4", "n4", "n6"), ("m5", "n3", "n4"), ("m6", "n5", "n6"))}
    e1 = r["end"]["m1"]
    add("two-storey-frame-lateral-loads", "One-bay two-storey frame under lateral joint loads",
        ("1x2 Plane Frame Lateral Load", "F24F1D5A-40F6-410E-91EA-A9F842327F1F"),
        "McCormac, Structural Analysis, 3rd ed. 1975, p. 388, example 22-7",
        PLANE + "Bay 20 ft, two storeys of 20 ft, fixed bases, all members A = 0.2 ft^2, I = 0.1 ft^4, E = 4.176e6 "
        "kip/ft^2; 20 kip at the first level and 10 kip at the roof, on the left column line, toward +X. Units kip, "
        "ft. First-order analysis with axial deformation. Target of the page: the largest moment, in the lower left "
        "column (member 1 of the listing).",
        model("", ("kip", "ft"), "Y", n, plane(n, {"n1": fix("ux", "uy", "rz"), "n2": fix("ux", "uy", "rz")}),
              {"M": {"E": 4.176e6, "nu": 0.3}}, {"S": sec(0.2, 0.1)}, mm, lc, cb),
        [R(r["u"]["n5"][0], disp("n5", "ux"), "ft", "plane frame program"),
         R(r["u"]["n3"][0], disp("n3", "ux"), "ft", "plane frame program"),
         R(r["R"]["n1"][0], reac("n1", "Fx"), "kip", "plane frame program"),
         R(r["R"]["n1"][1], reac("n1", "Fy"), "kip", "plane frame program"),
         R(r["R"]["n2"][2], reac("n2", "Mz"), "kip ft", "plane frame program"),
         R(max(abs(e1[2]), abs(e1[5])), extreme("m1", "Mz"), "kip ft", "plane frame program; largest end moment of the "
           "lower left column", "magnitude")],
        printed_bad=[P("176.40", extreme("m1", "Mz"), "kip ft", signnote="magnitude",
                       note="the page's 'theory' value of the largest moment (lower left column)")])

    r = reference("staad_two_bay_frame")
    n = {"n1": (0, 0), "n2": (216, 0), "n3": (432, 0), "n4": (0, 120), "n5": (216, 120), "n6": (432, 120)}
    lcs = {"P": {"loads": [{"type": "nodal", "node": "n4", "Fy": -50.0}, {"type": "nodal", "node": "n6", "Fy": -50.0},
                           {"type": "nodal", "node": "n5", "Fy": -100.0},
                           {"type": "memberPoint", "member": "a", "direction": "Y", "value": -100.0, "position": 0.5},
                           {"type": "memberPoint", "member": "b", "direction": "Y", "value": -100.0, "position": 0.5}]},
           "W": {"loads": [{"type": "memberDistributed", "member": "a", "direction": "Y", "w1": -10.0 / 12.0},
                           {"type": "memberDistributed", "member": "b", "direction": "Y", "w1": -10.0 / 12.0}]}}
    cbs = {"LC1": {"factors": {"P": 1.0}}, "LC2": {"factors": {"W": 1.0}}}
    r1, r2 = r["LC1"], r["LC2"]
    add("two-bay-frame-beam-forces", "Two-bay frame with beams pinned to the outer columns: beam moment and shear",
        ("2 Bay Frame Moments and Shear", "62DF48AF-F8CE-45D2-8907-5D1598FA9969"),
        "AISC, Manual of Steel Construction - Allowable Stress Design, 9th ed. 1989 (no page given)",
        PLANE + "Two bays of 18 ft, columns 10 ft; columns 12 in x 24 in (I = 13,824 in^4), beams 12 in x 30 in (I = "
        "27,000 in^4), every member A = 1e7 in^2 (the page's listing: axial deformation suppressed), E = 30,000 ksi. "
        "Outer column bases fixed, centre column base pinned. The two beams are continuous over the centre column and "
        "moment-released at the outer columns (the listing's member releases). Load case 1: 50 kip on the outer "
        "column tops, 100 kip on the centre column top, 100 kip at midspan of each beam; load case 2: 10 kip/ft on "
        "both beams. Targets at the middle of the first beam; in load case 1 the shear is the value on the side of "
        "the centre column (the page's section at 0.501 of the span). Units kip, in. First-order analysis.",
        model("", ("kip", "in"), "Y", n, plane(n, {"n1": fix("ux", "uy", "rz"), "n3": fix("ux", "uy", "rz"),
                                                    "n2": fix("ux", "uy")}),
              {"M": {"E": 30000.0, "nu": 0.3}}, {"col": sec(1.0e7, 13824.0), "beam": sec(1.0e7, 27000.0)},
              {"c1": frame("n1", "n4", "col"), "c2": frame("n2", "n5", "col"), "c3": frame("n3", "n6", "col"),
               "a": frame("n4", "n5", "beam", releases={"i": {"rz": True}}),
               "b": frame("n5", "n6", "beam", releases={"j": {"rz": True}})}, lcs, cbs,
              analysis={"shearDeformation": False, "resultStations": 3}),
        [R(r1["end"]["a1"][5], station("a", "Mz", rel=0.5, comb="LC1"), "kip in", "plane frame program"),
         R(-r1["end"]["a2"][1], station("a", "Vy", rel=0.5, side="right", comb="LC1"), "kip", "plane frame program"),
         R(r1["end"]["a1"][4], station("a", "Vy", rel=0.5, side="left", comb="LC1"), "kip", "plane frame program"),
         R(r1["R"]["n2"][1], reac("n2", "Fy", "LC1"), "kip", "plane frame program"),
         R(r2["end"]["a1"][5], station("a", "Mz", rel=0.5, comb="LC2"), "kip in", "plane frame program"),
         R(r2["end"]["a1"][4], station("a", "Vy", rel=0.5, comb="LC2"), "kip", "plane frame program"),
         R(r2["R"]["n2"][1], reac("n2", "Fy", "LC2"), "kip", "plane frame program")],
        [P("3375", station("a", "Mz", rel=0.5, comb="LC1"), "kip in", signnote="sagging: positive Mz (local y up)"),
         P("68.75", station("a", "Vy", rel=0.5, side="right", comb="LC1"), "kip", signnote=MAG),
         P("2430", station("a", "Mz", rel=0.5, comb="LC2"), "kip in", signnote="sagging: positive Mz (local y up)"),
         P("22.50", station("a", "Vy", rel=0.5, comb="LC2"), "kip", signnote=MAG)])

    r = reference("staad_settlement_frame")                  # reference computed in kN, cm; model in kN, mm
    n = {"n1": (0, 0), "n2": (0, 7500), "n3": (6000, 7500), "n4": (6000, 2500)}
    lcs = {"A": {"loads": [{"type": "supportDisplacement", "node": "n1", "uy": -10.0}]},
           "D": {"loads": [{"type": "supportDisplacement", "node": "n4", "uy": -10.0}]}}
    cbs = {"LC1": {"factors": {"A": 1.0}}, "LC2": {"factors": {"D": 1.0}}}
    tr, tp = [], []
    for cid in ("LC1", "LC2"):
        u = r[cid]["u"]
        s = 1.0 if cid == "LC2" else -1.0
        tr += [R(10.0 * u["n2"][0], disp("n2", "ux", cid), "mm", "plane frame program"),
               R(u["n2"][2], disp("n2", "rz", cid), "rad", "plane frame program"),
               R(10.0 * u["n4"][0], disp("n4", "ux", cid), "mm", "plane frame program")]
        tp += [P("1.25", disp("n2", "ux", cid), "mm", sign=s, factor=10.0, note="printed in cm", signnote=MAG),
               P("0.001667", disp("n2", "rz", cid), "rad", sign=-s, signnote=MAG),
               P("0.4167", disp("n4", "ux", cid), "mm", sign=s, factor=10.0, note="printed in cm", signnote=MAG)]
    add("frame-support-settlement", "Statically determinate frame: displacements caused by support settlements",
        ("2D Portal Reactions Sidesway 1", "8835257B-5327-43C0-B270-BA561A93D215"),
        "C. K. Wang, Intermediate Structural Analysis, 1983, section 2.11, p. 47",
        PLANE + "Frame with nodes A (0, 0), B (0, 7.5 m), C (6 m, 7.5 m), D (6 m, 2.5 m); A pinned, D on a roller that "
        "holds the vertical displacement only. Members AB and CD: A = 2400 cm^2, I = 720,000 cm^4; BC: A = 4800 cm^2, "
        "I = 1.44e6 cm^4; E = 2171.84 kN/cm^2. Load case 1: support A settles 1 cm; load case 2: support D settles "
        "1 cm. The frame is statically determinate: the displacements are a rigid-body motion. The page works in kN "
        "and cm; the model file is in kN and mm (the input format has no cm).",
        model("", ("kN", "mm"), "Y", n, plane(n, {"n1": fix("ux", "uy"), "n4": fix("uy")}),
              {"M": {"E": 21.7184, "nu": 0.17}}, {"col": sec(2400.0e2, 720000.0e4), "beam": sec(4800.0e2, 1.44e10)},
              {"m1": frame("n1", "n2", "col"), "m2": frame("n2", "n3", "beam"), "m3": frame("n3", "n4", "col")},
              lcs, cbs), tr, tp)

    r = reference("staad_space_frame")
    n = {"n1": (0, 120, 0), "n2": (240, 120, 0), "n3": (0, 0, 0), "n4": (360, 0, 120)}
    lc, cb = one([{"type": "nodal", "node": "n1", "Fx": 2.0}, {"type": "nodal", "node": "n2", "Fy": -1.0, "Mz": -120.0}])
    u, rr = r["u"], r["R"]
    comp = {"ux": 0, "uy": 1, "uz": 2, "rx": 3, "ry": 4, "rz": 5}
    tr = [R(u[node][comp[d]], disp(node, d), "in" if d[0] == "u" else "rad", "space frame program")
          for node, d in (("n1", "ux"), ("n1", "uy"), ("n1", "uz"), ("n1", "rz"), ("n2", "ux"), ("n2", "uy"),
                          ("n2", "uz"), ("n2", "rx"), ("n2", "rz"))]
    tr += [R(rr[node][k], reac(node, c), "kip" if c[0] == "F" else "kip in", "space frame program")
           for node, k, c in (("n3", 0, "Fx"), ("n3", 1, "Fy"), ("n3", 2, "Fz"), ("n3", 5, "Mz"), ("n4", 1, "Fy"),
                              ("n4", 2, "Fz"), ("n4", 3, "Mx"))]
    tr.append(R(r["N"]["m3"], station("m3", "N", rel=0.5), "kip", "space frame program", "tension positive"))
    add("space-frame-three-members", "Space frame of three members under joint loads",
        ("3D Frame Max Forces", "A250C718-03E4-4C7F-AAA6-A446F07B306A"),
        "Weaver, Computer Programs for Structural Analysis, p. 146, problem 8",
        "Space frame, global Y vertical: nodes 1 (0, 120, 0), 2 (240, 120, 0), 3 (0, 0, 0), 4 (360, 0, 120) in; members "
        "1-2, 3-1, 2-4, all A = 11 in^2, torsion constant 83 in^4, Iy = Iz = 56 in^4; E = 30,000 ksi, nu = 0.25 (G = "
        "12,000 ksi); nodes 3 and 4 built in. Loads: 2 kip along +X at node 1; 1 kip downward and a moment of 120 kip in "
        "about -Z at node 2. First-order analysis, no shear deformation. Of the page's three printed values only the "
        "largest axial force (member 2-4) is used: its largest bending moments 'about local y' and 'about local z' "
        "depend on the vendor's local-axis convention for an inclined member, which the page does not state.",
        model("", ("kip", "in"), "Y", n, {"n3": FIXED, "n4": FIXED}, {"M": {"E": 30000.0, "G": 12000.0}},
              {"S": sec(11.0, 56.0, 56.0, 83.0)},
              {"m1": frame("n1", "n2"), "m2": frame("n3", "n1"), "m3": frame("n2", "n4")}, lc, cb,
              analysis={"shearDeformation": False, "resultStations": 3}),
        tr, [P("1.47", station("m3", "N", rel=0.5), "kip", sign=math.copysign(1.0, r["N"]["m3"]), signnote=MAG)])

    r = reference("staad_simple_frame")
    n = {"n1": (0, 0), "n2": (0, 10), "n3": (0, 20), "n4": (10, 20)}
    lc, cb = one([{"type": "nodal", "node": "n4", "Fy": -1.0}])
    add("bent-cantilever-with-prop", "Bent cantilever with an intermediate lateral support: support reactions",
        ("Support Reactions for a Simple Frame", "1A25B902-67CC-4B2C-A469-37D370E85C34"),
        "Timoshenko, Strength of Materials, Part I, 3rd ed. 1956, p. 346, problem 2",
        PLANE + "Vertical leg from (0, 0) to (0, 20 in) with a horizontal arm to (10, 20 in); built in at the base, "
        "held horizontally at midheight (0, 10 in); 1 kip downward at the end of the arm; A = 10 in^2, I = 100 in^4, "
        "E = 3000 ksi. First-order analysis. The page's target is the horizontal reaction (1.5 kip).",
        model("", ("kip", "in"), "Y", n, plane(n, {"n1": fix("ux", "uy", "rz"), "n2": fix("ux")}),
              {"M": {"E": 3000.0, "nu": 0.17}}, {"S": sec(10.0, 100.0)},
              {"m1": frame("n1", "n2"), "m2": frame("n2", "n3"), "m3": frame("n3", "n4")}, lc, cb),
        [R(r["R"]["n2"][0], reac("n2", "Fx"), "kip", "plane frame program"),
         R(r["R"]["n1"][0], reac("n1", "Fx"), "kip", "plane frame program"),
         R(r["R"]["n1"][2], reac("n1", "Mz"), "kip in", "plane frame program"),
         R(r["u"]["n4"][1], disp("n4", "uy"), "in", "plane frame program")],
        [P("1.5", reac("n2", "Fx"), "kip", sign=math.copysign(1.0, r["R"]["n2"][0]), signnote=MAG)])

    # ------------------------------------------------------------------ second order, buckling
    def column(pid, title, page, cites, L, inertia, pload, hload, shape, printed_ok, printed_bad, extra=""):
        rr = reference(pid)
        nn = {"a": (0, 0), "b": (0, L)}
        lcs = {"H": {"loads": [{"type": "nodal", "node": "b", "Fx": hload}]},
               "P": {"loads": [{"type": "nodal", "node": "b", "Fy": -pload}]}}
        cbs = {"second": {"factors": {"H": 1.0, "P": 1.0}, "analysisType": "secondOrder"},
               "first": {"factors": {"H": 1.0}, "analysisType": "firstOrder"}}
        how = "ODE oracle, second-order beam-column"
        add(title[0], title[1], page, cites,
            PLANE + f"Cantilever column of {L:g} in, built in at the base, free at the top; {shape}, E = 29,000 ksi; "
            f"axial load {pload:g} kip and lateral load {hload:g} kip at the top. Combination 'second': second-order "
            "analysis of both loads (the engine's member carries the exact stability functions); combination 'first': "
            "first-order analysis of the lateral load alone. The area does not enter the targets (placeholder 10 "
            "in^2); Iy = 1000 Iz so that no out-of-plane buckling enters. " + extra,
            model("", ("kip", "in"), "Y", nn, plane(nn, {"a": fix("ux", "uy", "rz")}), {"M": {"E": 29000.0, "nu": 0.3}},
                  {"S": sec(10.0, inertia, 1000.0 * inertia, 1.0)}, {"m": frame("a", "b")}, lcs, cbs,
                  analysis={"shearDeformation": False, "resultStations": 3}),
            [R(rr["M_base"], station("m", "Mz", rel=0.0, comb="second"), "kip in", how),
             R(rr["tip"], disp("b", "ux", "second"), "in", how),
             R(rr["M_base_first"], station("m", "Mz", rel=0.0, comb="first"), "kip in", "ODE oracle, first order"),
             R(rr["tip_first"], disp("b", "ux", "first"), "in", "ODE oracle, first order")],
            printed_ok, printed_bad=printed_bad)
        return rr

    rr = reference("staad_pdelta_column")
    column("staad_pdelta_column", ("cantilever-column-p-delta", "Cantilever column with axial and lateral load: "
                                                                "second-order base moment"),
           ("Single Column P-Delta Analysis", "A271BC6F-0313-4D0A-8388-114D259A493F"),
           "no textbook; the page validates with M = H L + P x (first-order tip deflection)",
           144.0, 30.80, 3.0, 5.0, "W8x10 with the moment of inertia the page uses, I = 30.80 in^4",
           [P("5.57", disp("b", "ux", "first"), "in", note="first-order tip deflection H L^3 / (3 E I) of the page")],
           [P("736.72", station("m", "Mz", rel=0.0, comb="second"), "kip in",
              sign=math.copysign(1.0, rr["M_base"]), signnote=MAG,
              note="the page's reference base moment (its table header says ft kip; the value is in kip in)")])
    rr = reference("staad_direct_analysis_column")
    column("staad_direct_analysis_column", ("cantilever-column-notional-load", "Cantilever column with axial load and "
                                            "a small lateral load: second-order moment and deflection"),
           ("Direct Analysis of a Column", "4E631D44-D7D9-497E-85C5-ABFA4E11CF76"),
           "AISC, Specification for Structural Steel Buildings, 2005, p. 16.1-435, Case 2 (the 2005 edition was not "
           "opened; the same benchmark of the 2022 edition is in ../aisc-commentary)",
           300.0, 61.9, 30.0, 0.09, "W8x18 with the moment of inertia the page gives, Iz = 61.9 in^4", [],
           [P("61.37", station("m", "Mz", rel=0.0, comb="second"), "kip in", sign=math.copysign(1.0, rr["M_base"]),
              signnote=MAG, note="the page's reference maximum moment"),
            P("1.145", disp("b", "ux", "second"), "in", note="the page's reference maximum deflection")],
           extra="The lateral load is 0.003 x 30 kip = 0.09 kip (the page's text has '0.03 x 30'; its arithmetic uses "
                 "0.09). The page's closed form uses the full E (no stiffness reduction). ")

    r = reference("staad_column_buckling")
    nn = {"a": (0, 0), "b": (0, 360)}
    lc, cb = one([{"type": "nodal", "node": "b", "Fy": -160.0}])
    add("pinned-column-buckling-factor", "Pin-ended column: elastic buckling factor",
        ("Column Buckling Factor", "B476903C-918F-4E89-B6FD-093CBBF80A7D"),
        "no textbook; the page validates with Pcr = pi^2 E I / L^2",
        PLANE + "Column of 30 ft = 360 in, pinned at both ends (the top free to move axially), W10x45 with the moment "
        "of inertia the page uses, I = 248 in^4, E = 29,730 ksi; axial load 160 kip. Linear buckling analysis; target: "
        "the buckling factor Pcr / 160. The engine's member carries the exact stability functions. The area does not "
        "enter (placeholder 10 in^2); Iy = 1000 Iz so that the in-plane mode is the lowest.",
        model("", ("kip", "in"), "Y", nn, plane(nn, {"a": fix("ux", "uy"), "b": fix("ux")}),
              {"M": {"E": 29730.0, "nu": 0.3}}, {"S": sec(10.0, 248.0, 248000.0, 1.0)}, {"m": frame("a", "b")}, lc, cb,
              analysis={"shearDeformation": False, "buckling": {"combinations": ["C"], "modes": 1}}),
        [R(r["factor"], buck(1, "C"), "-", "smallest root of the characteristic equation sin(kL) = 0, P = k^2 EI")],
        [P("3.509", buck(1, "C"), "-", signnote="positive factor")])

    r = reference("staad_direct_analysis_beam")
    nn = {"a": (0, 0), "b": (480, 0)}
    lc, cb = one([{"type": "memberDistributed", "member": "m", "direction": "Y", "w1": -0.1 / 12.0},
                  {"type": "nodal", "node": "b", "Fx": -100.0}], "secondOrder")
    add("beam-column-uniform-load", "Simply supported beam-column with uniform load: second-order moment and "
        "deflection",
        ("Direct Analysis of a Beam", "98B68FD9-783E-4E16-88E9-FEB8CF028425"),
        "AISC, Specification for Structural Steel Buildings, 2005, p. 16.1-435, Case 1 (the 2005 edition was not "
        "opened; the same benchmark of the 2022 edition is in ../aisc-commentary)",
        PLANE + "Beam of 40 ft = 480 in, simply supported (the far end free to move axially), W10x60 with the moment "
        "of inertia the page gives, Iz = 341 in^4, E = 29,000 ksi (the page's closed form uses the full E); uniform "
        "load 0.1 kip/ft downward and axial compression 100 kip. Second-order analysis (exact stability functions). "
        "The area does not enter the targets (placeholder 10 in^2); Iy = 1000 Iz so that no out-of-plane buckling "
        "enters.",
        model("", ("kip", "in"), "Y", nn, plane(nn, {"a": fix("ux", "uy"), "b": fix("uy")}),
              {"M": {"E": 29000.0, "nu": 0.3}}, {"S": sec(10.0, 341.0, 341000.0, 1.0)}, {"m": frame("a", "b")}, lc, cb,
              analysis={"shearDeformation": False, "resultStations": 3}),
        [R(r["M_mid"], station("m", "Mz", rel=0.5), "kip in", "ODE oracle, second-order beam-column"),
         R(r["v_mid"], station("m", "uy", rel=0.5), "in", "ODE oracle, second-order beam-column")],
        [P("26.36", station("m", "Mz", rel=0.5), "kip in", factor=12.0, note="printed in ft kip",
           signnote="sagging: positive Mz (local y up)"),
         P("0.763", station("m", "uy", rel=0.5), "in", sign=-1.0, signnote="downward")])

    # ------------------------------------------------------------------ trusses
    r = reference("two_bar_truss")
    a, b = 180.0 * math.cos(math.radians(30.0)), 180.0 * math.sin(math.radians(30.0))
    n = {"n1": (0, 0), "n2": (a, -b), "n3": (2 * a, 0)}
    lc, cb = one([{"type": "nodal", "node": "n2", "Fy": -5000.0}])
    add("two-bar-truss", "Two-bar truss: bar force and joint deflection",
        ("Stress in a 2D Truss Model", "08C53183-5FAC-4C64-B3DC-984BE8F8DF38"),
        "Timoshenko, Strength of Materials, Part I, 3rd ed. 1956, p. 10, problem 2",
        "Plane truss (X-Y plane, Y vertical; uz of every joint held): two bars of 180 in, each inclined 30 degrees to "
        "the horizontal, from two pinned supports down to a common joint; A = 0.5 in^2, E = 30e6 psi; 5000 lb downward "
        "at the joint. The page prints the bar stress; the target here is the bar force = stress x 0.5 in^2.",
        model("", ("lbf", "in"), "Y", n, truss_supports(n, {"n1": ("ux", "uy"), "n3": ("ux", "uy")}),
              {"M": {"E": 30.0e6, "nu": 0.3}}, {"S": bar(0.5)}, {"m1": truss("n1", "n2", "S"), "m2": truss("n2", "n3", "S")},
              lc, cb),
        [R(r["end"]["m1"][3], station("m1", "N", rel=0.5), "lbf", "plane truss program", "tension positive"),
         R(r["u"]["n2"][1], disp("n2", "uy"), "in", "plane truss program")],
        [P("10,000.", station("m1", "N", rel=0.5), "lbf", factor=0.5, note="printed bar stress in psi",
           signnote="tension positive"),
         P("0.12", disp("n2", "uy"), "in", sign=-1.0, signnote="downward")])

    r = reference("staad_truss_reactions_3")
    n = {"n1": (0, 0), "n2": (150, 100), "n3": (150, 50), "n4": (300, 0)}
    lc, cb = one([{"type": "nodal", "node": "n2", "Fy": -10.0}])
    add("truss-five-bars-thrust", "Statically indeterminate five-bar truss: horizontal reaction",
        ("Reactions in a 2D Truss Model 3", "A59F1504-07A9-4DE5-B09E-0FD8E1A37AAB"),
        "Timoshenko, Strength of Materials, Part I, 3rd ed. 1956, p. 346, problem 3",
        "Plane truss (X-Y plane, Y vertical; uz of every joint held): joints 1 (0, 0), 2 (150, 100), 3 (150, 50), "
        "4 (300, 0) in; bars 1-2 and 2-4 with A = 5 in^2, 1-3 and 3-4 with A = 3 in^2, 2-3 with A = 2 in^2; E = 30,000 "
        "ksi; joints 1 and 4 pinned; 10 kip downward at joint 2. Target of the page: the horizontal reaction.",
        model("", ("kip", "in"), "Y", n, truss_supports(n, {"n1": ("ux", "uy"), "n4": ("ux", "uy")}),
              {"M": {"E": 30000.0, "nu": 0.3}}, {"a5": bar(5.0), "a3": bar(3.0), "a2": bar(2.0)},
              {"m1": truss("n1", "n2", "a5"), "m2": truss("n1", "n3", "a3"), "m3": truss("n2", "n3", "a2"),
               "m4": truss("n2", "n4", "a5"), "m5": truss("n3", "n4", "a3")}, lc, cb),
        [R(r["R"]["n4"][0], reac("n4", "Fx"), "kip", "plane truss program"),
         R(r["R"]["n1"][0], reac("n1", "Fx"), "kip", "plane truss program"),
         R(r["u"]["n2"][1], disp("n2", "uy"), "in", "plane truss program"),
         R(r["end"]["m3"][3], station("m3", "N", rel=0.5), "kip", "plane truss program", "tension positive")],
        [P("8.77", reac("n4", "Fx"), "kip", sign=math.copysign(1.0, r["R"]["n4"][0]), signnote=MAG)])

    r = reference("staad_truss_deflection")
    from vendor_reference import TRUSS20_MEMBERS, TRUSS20_NODES
    n = dict(TRUSS20_NODES)
    areas = sorted({a for *_x, a in TRUSS20_MEMBERS})
    lc, cb = one([{"type": "nodal", "node": k, "Fy": -20.0} for k in ("n2", "n4", "n5")])
    add("truss-thirteen-bars-deflection", "Thirteen-bar truss with an overhang: joint deflection",
        ("Deflections in a 2D Truss Model", "8E17DF60-937D-4637-B0EA-005517160C74"),
        "McCormac, Structural Analysis, 3rd ed. 1975, p. 271, example 18-2",
        "Plane truss (X-Y plane, Y vertical; uz of every joint held), joints and bars as in the page's input listing: "
        "bottom chord joints 1 to 5 at 15 ft spacing, upper joints 6 (15, 7.5 ft), 7 (30, 15 ft), 8 (45, 7.5 ft); bars "
        "2-6, 3-4, 4-8, 4-5 with A = 1 in^2; 1-2, 2-3 with 2 in^2; 3-6, 3-8 with 1.5 in^2; 3-7, 1-6, 5-8 with 3 in^2; "
        "6-7, 7-8 with 4 in^2; E = 30,000 ksi; joint 1 pinned, joint 3 on a roller (vertical support); 20 kip downward "
        "at joints 2, 4 and 5. Target of the page: the vertical deflection of joint 5. Units kip, in.",
        model("", ("kip", "in"), "Y", n, truss_supports(n, {"n1": ("ux", "uy"), "n3": ("uy",)}),
              {"M": {"E": 30000.0, "nu": 0.3}}, {f"a{a:g}": bar(a) for a in areas},
              {name: truss(i, j, f"a{a:g}") for name, i, j, a in TRUSS20_MEMBERS}, lc, cb),
        [R(r["u"]["n5"][1], disp("n5", "uy"), "in", "plane truss program"),
         R(r["u"]["n5"][0], disp("n5", "ux"), "in", "plane truss program"),
         R(r["u"]["n7"][1], disp("n7", "uy"), "in", "plane truss program"),
         R(r["R"]["n3"][1], reac("n3", "Fy"), "kip", "plane truss program"),
         R(r["end"]["m9"][3], station("m9", "N", rel=0.5), "kip", "plane truss program", "tension positive")],
        [P("2.63", disp("n5", "uy"), "in", sign=-1.0, signnote="downward")])

    r = reference("staad_space_truss")
    n = {"n1": (0, 0, 0), "n2": (96, 0, 0), "n3": (0, 72, 0), "n4": (48, 24, 72)}
    lc, cb = one([{"type": "nodal", "node": "n4", "Fz": 50.0}])
    add("space-truss-three-bars", "Three-bar space truss: bar forces",
        ("Axial Stress on a Truss Model", "ABE6D61C-8A85-4C02-A566-E2D6D9BD1942"),
        "Beer and Johnston, Vector Mechanics for Engineers, Statics and Dynamics, 1962, p. 47, problem 2.70",
        "Space truss: three bars from the pinned joints 1 (0, 0, 0), 2 (8, 0, 0), 3 (0, 6, 0) ft to the joint 4 "
        "(4, 2, 6) ft; A = 1 in^2, E = 30e6 psi; 50 lb along +Z at joint 4 (away from the plane of the supports). The "
        "page prints the bar stresses; with A = 1 in^2 they equal the bar forces in lb. Statically determinate. "
        "Units lb, in.",
        model("", ("lbf", "in"), "Z", n, {k: {"dofs": fix("ux", "uy", "uz")} for k in ("n1", "n2", "n3")},
              {"M": {"E": 30.0e6, "nu": 0.3}}, {"S": bar(1.0)},
              {"m1": truss("n1", "n4", "S"), "m2": truss("n2", "n4", "S"), "m3": truss("n3", "n4", "S")}, lc, cb),
        [R(r["N"][k], station(k, "N", rel=0.5), "lbf", "space truss program", "tension positive") for k in ("m1", "m2", "m3")]
        + [R(r["u"]["n4"][2], disp("n4", "uz"), "in", "space truss program")],
        [P(v, station(k, "N", rel=0.5), "lbf", note="printed bar stress in psi; A = 1 in^2", signnote="tension positive")
         for k, v in (("m1", "10.4"), ("m2", "31.2"), ("m3", "22.9"))])

    unknown = set(QUARANTINE) - set(lib.problems)
    if unknown:
        raise SystemExit("quarantine.py names unknown problems: " + ", ".join(sorted(unknown)))


def main():
    lib = Library(HERE, "staad")
    build(lib)
    lib.write()


if __name__ == "__main__":
    main()
