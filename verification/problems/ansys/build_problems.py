#!/usr/bin/env python3
"""Writes the problem directories of verification/problems/ansys/.

Source of the problem STATEMENTS: ANSYS, Inc., Mechanical APDL Verification Manual, Release 2025 R2
(public help site ansyshelp.ansys.com, pages Hlp_V_VM<n>.html; read 2026-11-22).  Each page states a
test case, cites a reference, and prints a "Target" value next to the program's value.  Nothing of a
page is reproduced here: the data that define the problem are given in this project's own words.

The cited textbooks were not opened, so the basis of every problem is
"recomputed-from-published-statement" (see ../vendor_common.py): reference values are computed in
verification/derive/vendor_reference.py; the page's printed Target values are recorded as additional
targets with the printed-precision tolerance, and those that do not reproduce are quarantined separately
with a diagnosis (quarantine.py).  The program's own results are not used as targets.

The problem sketches of the pages were not readable; a problem is used only where the page's TEXT (test
case, modelling notes, formulas) determines supports and load positions.  README.md lists what was
excluded and why.

Run:  python3 verification/problems/ansys/build_problems.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from vendor_common import (Library, P, R, add_vendor, bar, buck, disp, fix, frame, model, one, plane, reac,  # noqa: E402
                           reference, sec, station, truss, truss_supports)
from quarantine import QUARANTINE  # noqa: E402

BASE = "https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_vm/Hlp_V_"
PLANE = ("Plane problem in the global X-Y plane (Y vertical); the out-of-plane DOFs (uz, rx, ry) of every node are "
         "held; in-plane bending about local z. ")
MAG = "sign from the reference solution, engine convention (docs/CONVENTIONS.md); the page's sign convention is not stated"


def source(page, cites):
    vm, title = page
    return {"work": "ANSYS, Inc., Mechanical APDL Verification Manual (public help site ansyshelp.ansys.com)",
            "edition": "Release 2025 R2, read 2026-11-22",
            "location": f"{vm}, '{title}' ({BASE}{vm}.html)",
            "pages": "web page",
            "notes": f"problem statement taken from the vendor page (text; the sketch was not readable). The page "
                     f"cites: {cites} - NOT opened for this project; no number is taken from it. The page's 'Target' "
                     "values are recorded as printed; the program's own results are not used."}


def build(lib):
    def add(pid, *a, printed_bad=(), **kw):
        add_vendor(lib, source, pid, *a, printed_bad=printed_bad, diagnosis=QUARANTINE.get(pid + "-vendor-printed"),
                   whole_diagnosis=QUARANTINE.get(pid), **kw)

    # ------------------------------------------------------------------ VM2
    r = reference("ansys_vm2")
    w = 10000.0 / 12.0
    n = {"a": (0, 0), "b": (120, 0), "c": (360, 0), "d": (480, 0)}
    lc, cb = one([{"type": "memberDistributed", "member": k, "direction": "Y", "w1": -w} for k in ("o1", "o2")])
    add("vm2-beam-with-overhangs", "Beam with two loaded overhangs: moment in the middle portion and midspan deflection",
        ("VM2", "Beam Stresses and Deflections"),
        "Timoshenko, Strength of Materials, Part I, 3rd ed. 1955, p. 98, problem 4",
        PLANE + "Beam of 480 in on two supports 240 in apart, with two overhangs of 120 in; 10,000 lb/ft downward on "
        "both overhangs only; A = 50.65 in^2, Iz = 7892 in^4, depth 30 in, E = 30e6 psi. First-order analysis without "
        "shear deformation. The page prints the bending stress in the middle portion (compression at the bottom); "
        "the target here is the bending moment = stress x Iz / (depth / 2) = stress x 526.1333 in^3.",
        model("", ("lbf", "in"), "Y", n, plane(n, {"b": fix("ux", "uy"), "c": fix("uy")}), {"M": {"E": 30.0e6, "nu": 0.3}},
              {"S": sec(50.65, 7892.0)}, {"o1": frame("a", "b"), "s": frame("b", "c"), "o2": frame("c", "d")}, lc, cb,
              analysis={"shearDeformation": False, "resultStations": 3}),
        [R(r["M_mid"], station("s", "Mz", rel=0.5), "lbf in", "ODE oracle"),
         R(r["v_mid"], station("s", "uy", rel=0.5), "in", "ODE oracle"),
         R(r["v_end"], disp("a", "uy"), "in", "ODE oracle"), R(r["R"], reac("b", "Fy"), "lbf", "ODE oracle")],
        [P("0.182", station("s", "uy", rel=0.5), "in", signnote="upward (the middle portion bends upward)")],
        printed_bad=[P("11400.000", station("s", "Mz", rel=0.5), "lbf in", sign=-1.0, factor=7892.0 / 15.0,
                       note="printed bending stress in psi (printed negative, at the bottom flange)",
                       signnote="hogging in the middle portion: negative Mz (local y up)")])

    # ------------------------------------------------------------------ VM4
    r = reference("two_bar_truss")
    a, b = 180.0 * math.cos(math.radians(30.0)), 180.0 * math.sin(math.radians(30.0))
    n = {"n1": (0, 0), "n2": (a, -b), "n3": (2 * a, 0)}
    lc, cb = one([{"type": "nodal", "node": "n2", "Fy": -5000.0}])
    add("vm4-two-bar-truss", "Two hinged bars under a load at their common joint: bar force and deflection",
        ("VM4", "Deflection of a Hinged Support"),
        "Timoshenko, Strength of Materials, Part I, 3rd ed. 1955, p. 10, problem 2",
        "Plane truss (X-Y plane, Y vertical; uz of every joint held): two equal bars of 15 ft = 180 in with hinged "
        "ends, each inclined 30 degrees to the horizontal (the page's modelling notes: half span = l cos 30, rise = "
        "l sin 30), A = 0.5 in^2, E = 30e6 psi; 5000 lb at the common joint, in the direction that puts the bars in "
        "tension (taken downward with the joint below the supports). The page prints the bar stress; the target here "
        "is the bar force = stress x 0.5 in^2. The same textbook problem is in ../staad/two-bar-truss.",
        model("", ("lbf", "in"), "Y", n, truss_supports(n, {"n1": ("ux", "uy"), "n3": ("ux", "uy")}),
              {"M": {"E": 30.0e6, "nu": 0.3}}, {"S": bar(0.5)}, {"m1": truss("n1", "n2", "S"), "m2": truss("n2", "n3", "S")},
              lc, cb),
        [R(r["end"]["m1"][3], station("m1", "N", rel=0.5), "lbf", "plane truss program", "tension positive"),
         R(r["u"]["n2"][1], disp("n2", "uy"), "in", "plane truss program")],
        [P("10,000.", station("m1", "N", rel=0.5), "lbf", factor=0.5, note="printed bar stress in psi",
           signnote="tension positive"),
         P("0.120", disp("n2", "uy"), "in", sign=-1.0, signnote="printed negative; in the direction of the load")])

    # ------------------------------------------------------------------ VM21
    r = reference("ansys_vm21")
    n = {"a": (0, 0), "b": (200, 0)}
    inertia = 2.5 ** 4 / 12.0
    lcs = {"p": {"loads": [{"type": "memberDistributed", "member": "m", "direction": "Y", "w1": -1.79253}]},
           "F": {"loads": [{"type": "nodal", "node": "b", "Fx": 21972.6}]}}
    cbs = {"stiffened": {"factors": {"p": 1.0, "F": 1.0}, "analysisType": "secondOrder"},
           "unstiffened": {"factors": {"p": 1.0}, "analysisType": "firstOrder"}}
    tr, tp, tb = [], [], []
    printed = {"stiffened": ("0.19945", "0.0032352", "4580.1"), "unstiffened": ("0.38241", "0.0061185", "8962.7")}
    for cid in ("stiffened", "unstiffened"):
        how = "ODE oracle, second-order beam-column in tension" if cid == "stiffened" else "ODE oracle, first order"
        x = r[cid]
        tr += [R(x["v_mid"], station("m", "uy", rel=0.5, comb=cid), "in", how),
               R(x["theta0"], disp("a", "rz", cid), "rad", how),
               R(x["M_mid"], station("m", "Mz", rel=0.5, comb=cid), "lbf in", how)]
        z, th, mm = printed[cid]
        tp += [P(z, station("m", "uy", rel=0.5, comb=cid), "in", sign=-1.0,
                 signnote="printed negative; in the direction of the lateral load"),
               P(th, disp("a", "rz", cid), "rad", sign=-1.0, signnote="printed positive; " + MAG)]
        tb += [P(mm, station("m", "Mz", rel=0.5, comb=cid), "lbf in", signnote="printed negative; " + MAG)]
    add("vm21-tie-rod-lateral-load", "Tie rod with uniform lateral load: deflection, end slope and moment, with and "
        "without the tensile force",
        ("VM21", "Tie Rod with Lateral Loading"),
        "Timoshenko, Strength of Materials, Part II, 3rd ed. 1956, p. 42, article 6",
        PLANE + "Rod of 200 in, square section 2.5 in x 2.5 in (A = 6.25 in^2, I = 3.2552 in^4), E = 30e6 psi, uniform "
        "lateral load 1.79253 lb/in, tensile force 21,972.6 lb. Supports: both ends held laterally and free to rotate, "
        "one end free to move axially - the page's text asks for the slope at the end and uses symmetry about "
        "midspan (the sketch was not readable). Combination 'stiffened': second-order analysis with the tensile force "
        "(exact stability functions for tension); 'unstiffened': first-order analysis of the lateral load alone.",
        model("", ("lbf", "in"), "Y", n, plane(n, {"a": fix("ux", "uy"), "b": fix("uy")}), {"M": {"E": 30.0e6, "nu": 0.3}},
              {"S": sec(6.25, inertia)}, {"m": frame("a", "b")}, lcs, cbs,
              analysis={"shearDeformation": False, "resultStations": 3}), tr, tp, printed_bad=tb)

    # ------------------------------------------------------------------ VM41
    r = reference("ansys_vm41")
    n = {"a": (0, 0), "b": (10, 0)}
    lc, cb = one([{"type": "nodal", "node": "b", "Fy": -10.0}])
    add("vm41-stiff-beam-on-rotational-spring", "Very stiff beam pinned on a weak rotational spring: small rotation",
        ("VM41", "Small Deflection of a Rigid Beam"),
        "'any basic statics and strength of material book' (the page names no specific reference)",
        PLANE + "Beam AB of 10 in; A is held against translation and restrained against rotation by a spring of "
        "10,000 in lb/rad; lateral load of 10 lb at B. 'Thick beam' data of the page's modelling notes: A = 100 in^2, "
        "I = 1000 in^4, E = 30e6 psi (the beam is very stiff, not rigid: its bending adds 1.1e-7 in at B). First-order "
        "analysis.",
        model("", ("lbf", "in"), "Y", n, plane(n, {"a": {"ux": "fixed", "uy": "fixed", "rz": {"spring": 10000.0}}}),
              {"M": {"E": 30.0e6, "nu": 0.3}}, {"S": sec(100.0, 1000.0)}, {"m": frame("a", "b")}, lc, cb),
        [R(r["v_B"], disp("b", "uy"), "in", "ODE oracle, beam on a rotational spring support"),
         R(r["theta_B"], disp("b", "rz"), "rad", "ODE oracle"), R(r["theta_A"], disp("a", "rz"), "rad", "ODE oracle")],
        [P("0.1", disp("b", "uy"), "in", sign=-1.0, signnote="printed negative; in the direction of the load"),
         P("0.01", disp("b", "rz"), "rad", sign=-1.0, signnote="printed negative; rotation in the sense of the load"),
         P("0.0", disp("b", "ux"), "in", signnote="-")])

    # ------------------------------------------------------------------ VM127
    r = reference("ansys_vm127")
    n = {"a": (0, 0), "b": (0, 200)}
    lc, cb = one([{"type": "nodal", "node": "b", "Fy": -1.0}])
    i127 = 0.5 ** 4 / 12.0
    add("vm127-hinged-bar-buckling", "Slender bar with hinged ends: critical load",
        ("VM127", "Buckling of a Bar with Hinged Ends (Line Elements)"),
        "Timoshenko, Strength of Materials, Part II, 3rd ed. 1956, p. 148, article 29",
        PLANE + "Bar of 200 in with hinged ends (the far end free to move axially), square section 0.5 in x 0.5 in "
        "(A = 0.25 in^2, I = 0.0052083 in^4), E = 30e6 psi; reference load 1 lb, so the buckling factor is the critical "
        "load in lb. Linear buckling analysis of the whole bar (the page models half of it by symmetry); the engine's "
        "member carries the exact stability functions. Iy = 1000 Iz so that the in-plane mode is the lowest (the "
        "square section has no preferred plane; the page's problem is plane).",
        model("", ("lbf", "in"), "Y", n, plane(n, {"a": fix("ux", "uy"), "b": fix("ux")}), {"M": {"E": 30.0e6, "nu": 0.3}},
              {"S": sec(0.25, i127, 1000.0 * i127, 1.0)}, {"m": frame("a", "b")}, lc, cb,
              analysis={"shearDeformation": False, "buckling": {"combinations": ["C"], "modes": 1}}),
        [R(r["Pcr"], buck(1, "C"), "lbf", "smallest root of the characteristic equation sin(kL) = 0, P = k^2 EI")],
        [P("38.553", buck(1, "C"), "lbf", signnote="positive factor")])

    # ------------------------------------------------------------------ VM217
    r = reference("ansys_vm217")
    n = {"n1": (0, 0), "n2": (0, 400), "n3": (800, 400), "n4": (800, 0)}
    lc, cb = one([{"type": "memberDistributed", "member": "b", "direction": "Y", "w1": -500.0}])
    add("vm217-portal-frame-uniform-load", "Fixed-base portal frame with a stiff span under uniform load: joint "
        "rotation and largest moment",
        ("VM217", "Portal Frame Under Symmetric Loading"),
        "N. J. Hoff, The Analysis of Structures, 1956, pp. 115-119",
        PLANE + "Portal: columns of 400 in with I = 20,300 in^4, span of 800 in with I = 5 x 20,300 in^4, E = 30e6 psi, "
        "500 lb/in downward on the span. Column bases built in: the sketch was not readable, but the two formulas the "
        "page gives for its targets (rotation = w a^3 / (27 E I_col), largest moment = (19/54) w a^2) hold for built-in "
        "bases only (pinned bases give (3/8) w a^2 at midspan). The formulas neglect axial and shear deformation; the "
        "model follows them: every member A = 1e7 in^2, no shear deformation (the page's program model with finite "
        "area and shear deformation differs from its own targets by 9 % and 2 %). Targets: rotation of the "
        "beam-column joint, moment at midspan of the span.",
        model("", ("lbf", "in"), "Y", n, plane(n, {"n1": fix("ux", "uy", "rz"), "n4": fix("ux", "uy", "rz")}),
              {"M": {"E": 30.0e6, "nu": 0.3}}, {"col": sec(1.0e7, 20300.0), "span": sec(1.0e7, 101500.0)},
              {"c1": frame("n1", "n2", "col"), "b": frame("n2", "n3", "span"), "c2": frame("n4", "n3", "col")}, lc, cb,
              analysis={"shearDeformation": False, "resultStations": 3}),
        [R(r["u"]["n2"][2], disp("n2", "rz"), "rad", "plane frame program"),
         R(r["end"]["b1"][5], station("b", "Mz", rel=0.5), "lbf in", "plane frame program"),
         R(-r["end"]["b1"][2], station("b", "Mz", rel=0.0), "lbf in", "plane frame program"),
         R(r["R"]["n1"][0], reac("n1", "Fx"), "lbf", "plane frame program"),
         R(r["R"]["n1"][2], reac("n1", "Mz"), "lbf in", "plane frame program")],
        [P("0.195E-2", disp("n2", "rz"), "rad", sign=math.copysign(1.0, r["u"]["n2"][2]), signnote=MAG),
         P("0.281E8", station("b", "Mz", rel=0.5), "lbf in", signnote="sagging at midspan: positive Mz (local y up)")])

    unknown = set(QUARANTINE) - set(lib.problems)
    if unknown:
        raise SystemExit("quarantine.py names unknown problems: " + ", ".join(sorted(unknown)))


def main():
    lib = Library(HERE, "ansys")
    build(lib)
    lib.write()


if __name__ == "__main__":
    main()
