"""Problem definitions, chapters 4 (Ex. 4.14), 7, 8, 9, 12 and 13: member releases, shear deformation,
elastic critical loads, second-order analysis, flexible connections, rigid offsets, substructured frames,
symmetry examples, reanalysis.  Used by build_problems.py; see the notes at the top of that file.
"""
import math

import numpy as np


def define(B):                                               # B = the build_problems module
    P, R, add, model = B.P, B.R, B.add, B.model
    disp, reac, station, endf, buck = B.disp, B.reac, B.station, B.endf, B.buck
    bar, sec, truss, frame, sup, nodal, settle = B.bar, B.sec, B.truss, B.frame, B.sup, B.nodal, B.settle
    one_combo, SI, US, FIXED, PIN2D, OUT = B.one_combo, B.STEEL_SI, B.STEEL_US, B.FIXED, B.PIN2D, B.OUT
    PLANE, TRUSS, UNUSED = B.PLANE_NOTE, B.TRUSS_NOTE, B.UNUSED
    MOM = "moment / rotation about z, counter-clockwise positive in the source = rz / Mz of the engine"
    NOSHEAR = {"shearDeformation": False}
    SECOND = {"shearDeformation": False, "convergenceTolerance": 1.0e-10}

    def buckling(modes=1):
        return {"shearDeformation": False, "buckling": {"combinations": ["C"], "modes": modes}}

    STAB = ("Plane stability problem: the source's model has no out-of-plane degrees of freedom. Here the "
            "out-of-plane DOFs (uz, rx, ry) of every node are held and Iy = 1000 Iz, so that no out-of-plane buckling "
            "mode precedes the in-plane modes compared; J is a placeholder. ")
    EXACT = ("The engine's members carry the exact stability functions (no discretisation error for a prismatic "
             "member with constant axial force); the source's computed values use a stated number of elements with a "
             "polynomial (approximate) geometric stiffness. ")

    def psec(a, iz):
        return sec(a, iz, iy=1000.0 * iz)

    def outs(*dofs):
        return sup(*(dofs + OUT))

    def stiff_targets(cases, unit_dof):
        """cases: [(combination, node, component, value, formula)] -> recomputed targets."""
        return [R(v, reac(n, c, comb=cb), "kN" if c in ("Fx", "Fy") else "kN mm", f,
                  "stiffness coefficient x imposed displacement = reaction at the node (force on the structure by the "
                  "support), same axes and signs as the source's element forces") for cb, n, c, v, f in cases]

    # ------------------------------------------------------------------------------------- Ex. 4.14
    E, I, a, b = 200.0, 200e6, 4000.0, 2000.0
    K = 3.0 * E * I / (a ** 3 + b ** 3)
    th = 1.0e-3
    f414 = "Example 4.14 result: K = 3EI/(a^3 + b^3) [1, a, -1, b; a, a^2, -a, ab; -1, -a, 1, -b; b, ab, -b, b^2]"
    lcs = {"V": {"loads": [settle("n1", uy=1.0)]}, "R": {"loads": [settle("n1", rz=th)]}}
    cbs = {"V": {"factors": {"V": 1.0}}, "R": {"factors": {"R": 1.0}}}
    add("ex-4-14-beam-internal-hinge", "Two beam segments joined by a hinge: condensed stiffness coefficients",
        "Example 4.14 (stiffness matrix of the combined member)", "pp. 86-87",
        PLANE + "The source derives symbolically the stiffness matrix of a member 1-0-2 with a moment hinge at 0 "
        "(segments a and b, equal EI). Numeric instance chosen here: a = 4 m, b = 2 m, E = 200 kN/mm^2, I = 200e6 mm^4. "
        "Model: members n1-h and h-n2, the hinge as a full rz release at end j of the first member; both end nodes "
        "fully held. Combination V imposes v1 = 1 mm, combination R imposes theta1 = 0.001 rad; the reactions at n1 "
        "and n2 are the first and second column of the printed matrix times the imposed displacement. First order, no "
        "shear deformation. " + UNUSED,
        model("ex-4-14-beam-internal-hinge", ("kN", "mm"), "Y", {"n1": (0, 0), "h": (a, 0), "n2": (a + b, 0)},
              {"n1": FIXED, "n2": FIXED}, SI, {"s": sec(1.0e4, I)},
              {"m1": frame("n1", "h", "s", releases={"j": {"rz": True}}), "m2": frame("h", "n2", "s")}, lcs, cbs),
        stiff_targets([("V", "n1", "Fy", K, f414), ("V", "n1", "Mz", K * a, f414), ("V", "n2", "Fy", -K, f414),
                       ("V", "n2", "Mz", K * b, f414), ("R", "n1", "Fy", K * a * th, f414),
                       ("R", "n1", "Mz", K * a * a * th, f414), ("R", "n2", "Fy", -K * a * th, f414),
                       ("R", "n2", "Mz", K * a * b * th, f414)], "n1"),
        basis="recomputed-from-published-statement")

    # ------------------------------------------------------------------------------------- Ex. 7.10
    bb, dd, Pt = 100.0, 200.0, 10.0
    Iz = bb * dd ** 3 / 12.0
    s710 = sec(bb * dd, Iz, Asy=2.0 * bb * dd / 3.0)
    lc, cb = one_combo([nodal("b1", Fy=Pt), nodal("b2", Fy=Pt)])
    add("ex-7-10-cantilever-shear-deformation", "End-loaded cantilever of rectangular section with shear deformation",
        "Example 7.10", "p. 206",
        PLANE + "The source gives the tip deflection of an end-loaded cantilever including shear deformation as "
        "gamma P L^3 / (3 E I), with I = b d^3 / 12, shear area (2/3) b d and nu = 0.3, and prints gamma = 1.039 for "
        "d/L = 1/5 and gamma = 1.010 for d/L = 1/10. Numeric instance chosen here (the example is non-dimensional): "
        "b = 100 mm, d = 200 mm, E = 200 kN/mm^2, nu = 0.3, P = 10 kN upward at the tip; member 1 with L = 1000 mm "
        "(d/L = 1/5), member 2 with L = 2000 mm (d/L = 1/10). Section given by properties with Asy = (2/3) b d (the "
        "source's shear area); analysis with shear deformation. The target is the printed gamma times P L^3/(3 E I) "
        "of the instance. First order.",
        model("ex-7-10-cantilever-shear-deformation", ("kN", "mm"), "Y",
              {"a1": (0, 0), "b1": (1000, 0), "a2": (0, 5000), "b2": (2000, 5000)}, {"a1": FIXED, "a2": FIXED},
              SI, {"s": s710}, {"m1": frame("a1", "b1", "s"), "m2": frame("a2", "b2", "s")}, lc, cb,
              analysis={"shearDeformation": True}),
        [P("1.039", disp("b1", "uy"), "mm", factor=Pt * 1000.0 ** 3 / (3.0 * 200.0 * Iz),
           note="gamma for d/L = 1/5; the factor is P L^3/(3 E I) of the instance",
           signnote="load and deflection along +y"),
         P("1.010", disp("b2", "uy"), "mm", factor=Pt * 2000.0 ** 3 / (3.0 * 200.0 * Iz),
           note="gamma for d/L = 1/10; the factor is P L^3/(3 E I) of the instance",
           signnote="load and deflection along +y")],
        limits="The model has shear deformation: the foreign solvers of tier C cannot represent it "
               "(docs/VERIFICATION.md, Limits).")

    # ------------------------------------------------------------------------- Ex. 8.3 / 8.4 instances
    A93, I93 = 1.27e4, 3.66e7
    EI = 200.0 * I93
    L, al, Pc = 8000.0, 1.25e-4, 600.0
    kk = math.sqrt(Pc / EI)
    delta = al * L * (1.0 / math.cos(kk * L / 2.0) - 1.0)
    lc, cb = one_combo([nodal("b", Fy=-Pc, Mz=al * Pc * L), nodal("a", Mz=-al * Pc * L)], "secondOrder")
    add("ex-8-3-beam-column-end-moments", "Pin-ended beam-column with equal end moments: midheight deflection",
        "Example 8.3, Equation f (the member of Example 9.3)", "pp. 223-224 (data: p. 249)",
        PLANE + "The source gives in closed form (small-displacement theory) the midheight deflection of a pin-ended "
        "bar under an axial load P and end moments M0 = alpha P L (single curvature): delta = alpha L [sec(sqrt(P/EI) "
        "L/2) - 1]. Numeric instance: the member of Example 9.3 (L = 8 m, A = 1.27e4 mm^2, I = 3.66e7 mm^4, "
        "E = 200,000 MPa) with alpha = 1.25e-4 as in the statement of Example 9.3 and P = 600 kN chosen here (the source "
        "shows the response of this member only as a curve). Column along +Y, base a pinned, top b held laterally; "
        "moments -M0 (about Z) at a and +M0 at b bend the member toward +X. Second-order analysis, one member, no shear "
        "deformation. " + STAB,
        model("ex-8-3-beam-column-end-moments", ("kN", "mm"), "Y", {"a": (0, 0), "b": (0, L)},
              {"a": outs("ux", "uy"), "b": outs("ux")}, SI, {"s": psec(A93, I93)}, {"ab": frame("a", "b", "s")},
              lc, cb, analysis=SECOND),
        [R(delta, station("ab", "uy", rel=0.5), "mm", "Example 8.3, Eq. f",
           "deflection toward +X = +local y of the vertical member (default orientation)")],
        basis="recomputed-from-published-statement")

    L4, al4, Pc4, Pt4 = 4000.0, 0.005, 600.0, 5000.0
    kc, kt = math.sqrt(Pc4 / EI) * L4, math.sqrt(Pt4 / EI) * L4
    dc = al4 * L4 * (math.tan(kc) / kc - 1.0)
    dt = al4 * L4 * (1.0 - math.tanh(kt) / kt)
    lcs = {"LC": {"loads": [nodal("b", Fy=-Pc4, Fx=al4 * Pc4)]}, "LT": {"loads": [nodal("b", Fy=Pt4, Fx=al4 * Pt4)]}}
    cbs = {"C": {"factors": {"LC": 1.0}, "analysisType": "secondOrder"},
           "T": {"factors": {"LT": 1.0}, "analysisType": "secondOrder"}}
    add("ex-8-4-cantilever-beam-column", "Cantilever with axial and lateral end load: compression and tension",
        "Example 8.4, Equations e and g (the member of Example 9.4)", "pp. 225-226 (data: p. 250)",
        PLANE + "The source gives in closed form the tip deflection of a cantilever under a vertical load P and a "
        "horizontal load alpha P at the tip: Delta = alpha L [tan(kL)/(kL) - 1] for P downward (Eq. e) and "
        "Delta = alpha L [1 - tanh(kL)/(kL)] for P upward (Eq. g), k = sqrt(P/EI). Numeric instance: the member of "
        "Example 9.4 (L = 4 m, A = 1.27e4 mm^2, I = 3.66e7 mm^4, E = 200,000 MPa, alpha = 0.005) with P = 600 kN "
        "(combination C, compression) and P = 5000 kN (combination T, tension) chosen here; the source shows these "
        "responses only as curves. Second-order analysis, one member, no shear deformation. " + STAB,
        model("ex-8-4-cantilever-beam-column", ("kN", "mm"), "Y", {"a": (0, 0), "b": (0, L4)},
              {"a": FIXED, "b": outs()}, SI, {"s": psec(A93, I93)}, {"ab": frame("a", "b", "s")}, lcs, cbs,
              analysis=SECOND),
        [R(dc, disp("b", "ux", "C"), "mm", "Example 8.4, Eq. e", "deflection in the direction of the lateral load (+X)"),
         R(dt, disp("b", "ux", "T"), "mm", "Example 8.4, Eq. g", "deflection in the direction of the lateral load (+X)")],
        basis="recomputed-from-published-statement")

    # ------------------------------------------------------------------------------------- Ex. 9.1
    lc, cb = one_combo([nodal("b", Fy=-1.0)])
    add("ex-9-1-truss-critical-load", "Two-bar truss: elastic critical load of a bar restrained by a weak tie",
        "Example 9.1, case 1 (alpha = 0), Equation d", "p. 247",
        TRUSS + "a (0, 4) m pinned, b (4, 4) m, c (4, 0) m pinned; axial-force members ab (horizontal, A = 2 mm^2) and "
        "bc (vertical, A = 5e3 mm^2), E = 200,000 MPa. Reference load 1 kN downward at b, so the buckling factor is the "
        "critical load in kN. Linear buckling analysis; the mode is the rotation of bar bc about c, restrained by the "
        "axial stiffness of ab (P = k L = 0.1 kN/mm x 4000 mm). The source's axial-force member has no flexural "
        "buckling of its own; the engine treats a compressed truss member as a pin-ended beam-column with the section's "
        "I, so bar bc is given I = 1e10 mm^4 (member Euler load 1.2e6 kN; with 1e8 or 1e12 mm^4 the system factor is "
        "the same to 1e-7). The second case of the example (alpha = 0.05, "
        "large-displacement response with a limit point) is not applicable to the engine.",
        model("ex-9-1-truss-critical-load", ("kN", "mm"), "Y", {"a": (0, 4000), "b": (4000, 4000), "c": (4000, 0)},
              {"a": PIN2D, "c": PIN2D, "b": sup("uz")}, SI,
              {"ab": bar(2.0), "bc": {"type": "properties", "description": "axial-force member; I = 1e10 mm^4 keeps "
                                      "its own Euler load (1.2e6 kN) out of the comparison", "A": 5e3, "Iy": 1e10,
                                      "Iz": 1e10, "J": 1.0}},
              {"ab": truss("a", "b", "ab"), "bc": truss("c", "b", "bc")}, lc, cb, analysis=buckling()),
        [P("400", buck(), "kN (factor on 1 kN)", note="P_cre = lambda P_ref = 400 kN, an exact closed-form value "
           "(the source also derives it as k L = 2 x 200)", signnote="positive factor on a downward (compressive) load")])

    # ------------------------------------------------------------------------------- Ex. 9.3 / 9.4
    lc, cb = one_combo([nodal("b", Fy=-1.0)])
    add("ex-9-3-column-pinned-critical-load", "Pin-ended column: elastic critical load",
        "Example 9.3, case 1 (alpha = 0): value 'from Example 8.3'", "pp. 249-250",
        PLANE + "Column of 8 m, A = 1.27e4 mm^2, I = 3.66e7 mm^4, E = 200,000 MPa, base pinned, top held laterally. "
        "Reference load 1 kN downward at the top; the buckling factor is the critical load in kN. Linear buckling "
        "analysis, one member. Target: the theoretical value pi^2 EI/L^2 that the source prints as 1,129 kN. " + EXACT
        + "The source's two-element value, 1,137 kN, is its own approximation (0.7 % above its theoretical value) and is "
        "not a target. " + STAB,
        model("ex-9-3-column-pinned-critical-load", ("kN", "mm"), "Y", {"a": (0, 0), "b": (0, 8000)},
              {"a": outs("ux", "uy"), "b": outs("ux")}, SI, {"s": psec(A93, I93)}, {"ab": frame("a", "b", "s")},
              lc, cb, analysis=buckling()),
        [P("1,129", buck(), "kN (factor on 1 kN)", note="theoretical value pi^2 EI/L^2 as printed",
           signnote="positive factor on a compressive load")])
    add("ex-9-4-column-cantilever-critical-load", "Cantilever column: elastic critical load",
        "Example 9.4, case 1 (alpha = 0)", "p. 250",
        PLANE + "Column of 4 m, A = 1.27e4 mm^2, I = 3.66e7 mm^4, E = 200,000 MPa, base fixed, top free. Reference "
        "load 1 kN downward at the top. Linear buckling analysis, one member. Target: 1,129 kN, printed both as the "
        "theoretical value pi^2 EI/(4 L^2) and as the source's two-element result. " + EXACT + STAB,
        model("ex-9-4-column-cantilever-critical-load", ("kN", "mm"), "Y", {"a": (0, 0), "b": (0, 4000)},
              {"a": FIXED, "b": outs()}, SI, {"s": psec(A93, I93)}, {"ab": frame("a", "b", "s")}, lc, cb,
              analysis=buckling()),
        [P("1,129", buck(), "kN (factor on 1 kN)",
           note="theoretical value pi^2 EI/(4 L^2); the source's two-element result is printed with the same digits",
           signnote="positive factor on a compressive load")])

    # ------------------------------------------------------------------------------------- Ex. 9.5
    n95 = {"a": (0, 0), "b": (0, 6000), "c": (0, 10000), "d": (4000, 6000), "e": (4000, 10000)}
    s95 = {"a": outs("ux", "uy"), "b": outs(), "c": outs(), "d": outs("uy"), "e": outs("uy")}
    sec95 = {"col": psec(2.5e4, 6.36e8), "beam": psec(1.76e4, 8.61e8)}
    m95 = {"ab": frame("a", "b", "col"), "bc": frame("b", "c", "col"), "bd": frame("b", "d", "beam"),
           "ce": frame("c", "e", "beam")}
    geo95 = ("a (0, 0), b (0, 6), c (0, 10), d (4, 6), e (4, 10) m. Column ab-bc: A = 2.5e4 mm^2, I = 6.36e8 mm^4; "
             "beams bd, ce (rigidly connected to the column): A = 1.76e4 mm^2, I = 8.61e8 mm^4; E = 200,000 MPa. a is "
             "pinned; d and e are on rollers that hold the vertical displacement only (as drawn in the source). ")
    lc, cb = one_combo([nodal("c", Fy=-1.0)])
    add("ex-9-5-frame-critical-load", "Two-storey column with beams: elastic critical load",
        "Example 9.5, case 1 (alpha = 0)", "p. 251",
        PLANE + geo95 + "Reference load 1 kN downward at c; the buckling factor is the critical load in kN. Linear "
        "buckling analysis. " + EXACT + "The source models each member as a single element. " + STAB,
        model("ex-9-5-frame-critical-load", ("kN", "mm"), "Y", n95, s95, SI, sec95, m95, lc, cb, analysis=buckling()),
        [P("6630", buck(), "kN (factor on 1 kN)",
           note="computed by the source with one element per member and the polynomial geometric stiffness of its "
                "Eq. 9.14: an approximation by construction, no error bound stated",
           signnote="positive factor on a compressive load")])
    Pv = 3978.0
    ld95 = [nodal("c", Fy=-Pv, Fx=0.01 * Pv), nodal("b", Fx=0.01 * Pv)]
    SG95 = ("the source's moment diagram gives magnitudes; sign by statics in the engine's convention: the lateral "
            "loads act along +X, the rollers at d and e react upward, so the beams sag at the column (positive Mz at "
            "their end i), and the column moment at b below the joint is -(base shear) x 6 m (negative Mz at end j of "
            "ab, local y = +X)")
    lc, cb = one_combo(ld95)
    add("ex-9-5-frame-first-order", "Two-storey column with beams: first-order moments under vertical and lateral load",
        "Example 9.5, case 2 (alpha = 0.01, P = 3978 kN), moment diagram for Eq. 8.1", "p. 251",
        PLANE + geo95 + "Loads: P = 3978 kN downward at c, 0.01 P along +x at c and at b. First-order analysis. "
        "Targets: the first-order values of the source's moment diagram (figure labels, kN m): 477 in the column below "
        "b, 461 in beam bd at b, 175 at joint c (beam ce). The assignment of the labels to member ends follows from "
        "statics: base shear 2 x 39.78 kN x 6 m = 477 kN m; 461 + 175 = overturning moment / 1. " + UNUSED,
        model("ex-9-5-frame-first-order", ("kN", "mm"), "Y", n95, s95, SI, sec95, m95, lc, cb, analysis=NOSHEAR),
        [P("477", station("ab", "Mz", rel=1.0), "kN mm", sign=-1, factor=1000.0, note="figure label (kNm)", signnote=SG95),
         P("461", station("bd", "Mz", rel=0.0), "kN mm", factor=1000.0, note="figure label (kNm)", signnote=SG95),
         P("175", station("ce", "Mz", rel=0.0), "kN mm", factor=1000.0, note="figure label (kNm)", signnote=SG95)])
    lc, cb = one_combo(ld95, "secondOrder")
    add("ex-9-5-frame-second-order", "Two-storey column with beams: second-order moments under vertical and lateral load",
        "Example 9.5, case 2 (alpha = 0.01, P = 3978 kN), moment diagram for Eq. 8.3", "p. 251",
        PLANE + geo95 + "Loads: P = 3978 kN (0.6 of the source's critical load) downward at c, 0.01 P along +x at c "
        "and at b. Second-order analysis. Targets: the second-order values of the source's moment diagram (figure "
        "labels, kN m): 932 in the column below b, 869 in beam bd at b, 272 at joint c. " + EXACT + "The source models "
        "each member as a single element and states that column bowing is not included in the displayed results. " + STAB,
        model("ex-9-5-frame-second-order", ("kN", "mm"), "Y", n95, s95, SI, sec95, m95, lc, cb, analysis=SECOND),
        [P("932", station("ab", "Mz", rel=1.0), "kN mm", sign=-1, factor=1000.0,
           note="figure label (kNm); source: one element per member, incremental analysis", signnote=SG95),
         P("869", station("bd", "Mz", rel=0.0), "kN mm", factor=1000.0,
           note="figure label (kNm); source: one element per member, incremental analysis", signnote=SG95),
         P("272", station("ce", "Mz", rel=0.0), "kN mm", factor=1000.0,
           note="figure label (kNm); source: one element per member, incremental analysis", signnote=SG95)])

    # ------------------------------------------------------------------------------------- Ex. 9.6
    n96 = {"a": (0, 0), "b": (0, 288), "c": (240, 288), "d": (720, 288), "e": (720, 0)}
    s96 = {"a": outs("ux", "uy"), "e": outs("ux", "uy"), "b": outs(), "c": outs(), "d": outs()}
    sec96 = {"col": psec(13.3, 248.0), "gir": psec(24.8, 2850.0)}
    m96 = {"ab": frame("a", "b", "col"), "bc": frame("b", "c", "gir"), "cd": frame("c", "d", "gir"),
           "ed": frame("e", "d", "col")}
    geo96 = ("Units kip, in. a (0, 0), b (0, 24 ft), c (20 ft, 24 ft), d (60 ft, 24 ft), e (60 ft, 0); a and e pinned. "
             "Columns ab, ed: A = 13.3 in^2, I = 248 in^4; girder bc-cd: A = 24.8 in^2, I = 2850 in^4; E = 29,000 ksi. "
             "Loads: 6 kips along +x at b, 60 kips downward at c, 120 kips downward at d. ")
    ld96 = [nodal("b", Fx=6.0), nodal("c", Fy=-60.0), nodal("d", Fy=-120.0)]
    SG96 = ("the source's moment diagram gives magnitudes; sign by statics in the engine's convention (column ab from "
            "a to b, column ed from e to d, local y = +X; girder y up): girder moment at c sagging = positive Mz; "
            "column-top moment = -(base reaction along X on the structure... i.e. Mz(j) = -Vy L with Vy(0) = -R_x) ")
    lc, cb = one_combo(ld96)
    add("ex-9-6-portal-critical-load", "Unsymmetrically loaded pin-based portal frame: critical load factor",
        "Example 9.6, part 1", "p. 252",
        PLANE + geo96 + "Linear buckling analysis with the applied loads as reference loading; target: the source's "
        "lambda = 2.20. " + EXACT + "The source models each column with two elements and each girder segment with one. "
        + STAB,
        model("ex-9-6-portal-critical-load", ("kip", "in"), "Y", n96, s96, US, sec96, m96, lc, cb,
              analysis=buckling()),
        [P("2.20", buck(), "-", note="computed by the source with two elements per column and the polynomial "
           "geometric stiffness of its Eq. 9.14: an approximation by construction, no error bound stated",
           signnote="positive factor on the applied loading")])
    add("ex-9-6-portal-first-order", "Unsymmetrically loaded pin-based portal frame: first-order moments",
        "Example 9.6, part 2, moment diagram for Eq. 8.1", "p. 252",
        PLANE + geo96 + "First-order analysis. Targets: the first-order values of the source's moment diagram (figure "
        "labels, in-kips): 316 at b, 8708 at c, 2044 at d. Which label of each pair is the first-order one follows from "
        "statics (the column-top moments divided by the 288 in column height are the base shears, and "
        "2044/288 - 316/288 = 6.00 kips = the applied lateral load). " + UNUSED,
        model("ex-9-6-portal-first-order", ("kip", "in"), "Y", n96, s96, US, sec96, m96, lc, cb, analysis=NOSHEAR),
        [P("316", station("ab", "Mz", rel=1.0), "kip in", note="figure label (in. kips)",
           signnote="magnitude in the source; statics: the base reaction at a acts along +X (gravity thrust 4.1 kips "
                    "inward minus the share of the lateral load), giving Mz = +R_x L at the top of ab (local y = +X)"),
         P("8708", station("bc", "Mz", rel=1.0), "kip in", note="figure label (in. kips)",
           signnote="magnitude in the source; sagging under the 60 kip load = positive Mz (y up)"),
         P("2044", station("ed", "Mz", rel=1.0), "kip in", sign=-1, note="figure label (in. kips)",
           signnote="magnitude in the source; statics: the base reaction at e acts along -X, giving a negative Mz at "
                    "the top of ed (local y = +X)")])
    lc, cb = one_combo(ld96, "secondOrder")
    add("ex-9-6-portal-second-order", "Unsymmetrically loaded pin-based portal frame: second-order moments",
        "Example 9.6, part 2, moment diagram for Eq. 8.3", "p. 252",
        PLANE + geo96 + "Second-order analysis at the applied loads (load factor 1). Targets: the second-order values of "
        "the source's moment diagram (figure labels, in-kips): 519 at b (the source states that the moment at b changes "
        "sign), 9078 at c, 2661 at d. " + EXACT + "The source models each column with two elements and each girder "
        "segment with one. " + STAB,
        model("ex-9-6-portal-second-order", ("kip", "in"), "Y", n96, s96, US, sec96, m96, lc, cb, analysis=SECOND),
        [P("519", station("ab", "Mz", rel=1.0), "kip in", sign=-1,
           note="figure label (in. kips); source: two elements per column, incremental analysis",
           signnote="magnitude in the source, which states that the moment at b changes sign relative to the "
                    "first-order analysis (+316 in the engine's convention)"),
         P("9078", station("bc", "Mz", rel=1.0), "kip in",
           note="figure label (in. kips); source: two elements per column, incremental analysis",
           signnote="magnitude in the source; sagging = positive Mz"),
         P("2661", station("ed", "Mz", rel=1.0), "kip in", sign=-1,
           note="figure label (in. kips); source: two elements per column, incremental analysis",
           signnote="magnitude in the source; same sign as the first-order moment at the top of ed")])

    # ------------------------------------------------------------------------------------- Ex. 9.8
    n98 = {"a": (0, 0), "b": (0, 156), "c": (0, 312)}
    s98 = {"a": outs("ux", "uy"), "b": outs(), "c": outs("ux")}
    m98 = {"ab": frame("a", "b", "col"), "bc": frame("b", "c", "col")}
    geo98 = ("Units kip, in. Column ac of 2 x 13 ft, A = 7.08 in^2, I = 18.2 in^4, E = 29,000 ksi, base a pinned, top c "
             "held laterally, node b at midheight. Reference load 1 kip downward at c; the buckling factor is the "
             "critical load in kips. Linear buckling analysis. ")
    lc, cb = one_combo([nodal("c", Fy=-1.0)])
    add("ex-9-8-column-unbraced-critical-load", "Column with a midheight strut of zero area: critical load",
        "Example 9.8, part 1, sketch (a): value at A_br = 0", "p. 254",
        PLANE + geo98 + "Brace area zero (no strut in the model). Target: the value 53.5 kips labelled on the source's "
        "plot of critical load against brace area (computed by the source with four elements). " + EXACT + STAB,
        model("ex-9-8-column-unbraced-critical-load", ("kip", "in"), "Y", n98, s98, US, {"col": psec(7.08, 18.2)},
              m98, lc, cb, analysis=buckling()),
        [P("53.5", buck(), "kip (factor on 1 kip)", note="axis label of the source's plot; four-element model",
           signnote="positive factor on a compressive load")])
    n98b = dict(n98)
    n98b["d"] = (120, 156)
    s98b = dict(s98)
    s98b["d"] = PIN2D
    m98b = dict(m98)
    m98b["bd"] = truss("b", "d", "strut")
    add("ex-9-8-column-braced-critical-load", "Column fully braced at midheight by an elastic strut: critical load",
        "Example 9.8, part 1, sketch (a): plateau value; part 2 brace area", "p. 254",
        PLANE + geo98 + "Horizontal strut bd at midheight, 10 ft long, A_br = 1.27 in^2 (the brace area of part 2, far "
        "above the 0.01136 in^2 the source derives as the minimum for buckling in two half waves), far end d held. "
        "Target: the plateau value 215.7 kips labelled on the source's plot (four-element model). " + EXACT + STAB,
        model("ex-9-8-column-braced-critical-load", ("kip", "in"), "Y", n98b, s98b, US,
              {"col": psec(7.08, 18.2), "strut": bar(1.27)}, m98b, lc, cb, analysis=buckling()),
        [P("215.7", buck(), "kip (factor on 1 kip)",
           note="axis label of the source's plot; computed by the source with four elements and the polynomial "
                "geometric stiffness of its Eq. 9.14: an approximation by construction, no error bound stated",
           signnote="positive factor on a compressive load")])

    # ------------------------------------------------------------------------------------- Ex. 12.7
    lc, cb = one_combo([nodal("b", Fy=-1.0)])
    add("ex-12-7-cantilever-column-critical-loads", "Cantilever column: first and second elastic critical load",
        "Examples 12.6, 12.7 and 12.9: theoretical values", "pp. 362-363, 367",
        PLANE + "Column of 10 m, A = 7.2e3 mm^2, I = 150e6 mm^4, E = 200,000 MPa, base fixed, top free. Reference load "
        "1 kN downward at the top. Linear buckling analysis, one member, two modes. Targets: the theoretical values the "
        "source prints for comparison, pi^2 EI/(4 L^2) = 740 kN and 9 pi^2 EI/(4 L^2) = 6662 kN. " + EXACT
        + "The source's one-element results (746 kN, stated to be within 1.5 %, and 9651 kN) and its two-element second "
        "mode (6884 kN, within 3.3 %) are its own approximations and are not targets. " + STAB,
        model("ex-12-7-cantilever-column-critical-loads", ("kN", "mm"), "Y", {"a": (0, 0), "b": (0, 10000)},
              {"a": FIXED, "b": outs()}, SI, {"s": psec(7.2e3, 150e6)}, {"ab": frame("a", "b", "s")}, lc, cb,
              analysis=buckling(2)),
        [P("740", buck(1), "kN (factor on 1 kN)", note="theoretical value as printed in Example 12.7",
           signnote="positive factor on a compressive load"),
         P("6662", buck(2), "kN (factor on 1 kN)", note="theoretical value as printed in Example 12.9",
           signnote="positive factor on a compressive load")])

    # ------------------------------------------------------------------------------------- Ex. 13.2
    lc, cb = one_combo([nodal("n1", Fy=-50.0), nodal("n2", Fy=-30.0), nodal("n3", Fy=-10.0)])
    add("ex-13-2-fixed-beam-three-loads", "Fixed-ended beam with three point loads",
        "Example 13.2", "pp. 380-381",
        PLANE + "Beam of 8 m, both ends fixed, I = 200e6 mm^4, E = 200,000 MPa, nodes every 2 m. Loads 50, 30 and "
        "10 kN downward at x = 2, 4 and 6 m. First order, no shear deformation. " + UNUSED,
        model("ex-13-2-fixed-beam-three-loads", ("kN", "mm"), "Y",
              {"l": (0, 0), "n1": (2000, 0), "n2": (4000, 0), "n3": (6000, 0), "r": (8000, 0)},
              {"l": FIXED, "r": FIXED}, SI, {"s": sec(1.0e4, 200e6)},
              {"m1": frame("l", "n1", "s"), "m2": frame("n1", "n2", "s"), "m3": frame("n2", "n3", "s"),
               "m4": frame("n3", "r", "s")}, lc, cb),
        [P("2.54", disp("n1", "uy"), "mm", sign=-1, note="printed -2.54"),
         P("1.56e-3", disp("n1", "rz"), "rad", sign=-1, note="printed -1.56 x 10^-3", signnote=MOM),
         P("4.0", disp("n2", "uy"), "mm", sign=-1, note="printed -4.0"),
         P("0.25e-3", disp("n2", "rz"), "rad", note="printed 0.25 x 10^-3", signnote=MOM),
         P("1.96", disp("n3", "uy"), "mm", sign=-1, note="printed -1.96"),
         P("1.48e-3", disp("n3", "rz"), "rad", note="printed 1.48 x 10^-3", signnote=MOM)])

    # ------------------------------------------------------------------------------------- Ex. 13.3
    AR = B.A_RIGID
    n133 = {"n7": (0, 0), "n8": (8000, 0), "n5": (0, 5000), "n6": (8000, 5000), "n3": (0, 10000),
            "n4": (8000, 10000), "n1": (0, 15000), "n2": (8000, 15000)}
    m133 = {"c75": frame("n7", "n5", "col"), "c53": frame("n5", "n3", "col"), "c31": frame("n3", "n1", "col"),
            "c86": frame("n8", "n6", "col"), "c64": frame("n6", "n4", "col"), "c42": frame("n4", "n2", "col"),
            "b12": frame("n1", "n2", "beam"), "b34": frame("n3", "n4", "beam"), "b56": frame("n5", "n6", "beam")}
    lc, cb = one_combo([nodal(n, Fx=3.0, Mz=2.0e3) for n in ("n1", "n3", "n5")])
    add("ex-13-3-three-storey-frame", "Three-storey single-bay rigid frame, axial deformation neglected",
        "Example 13.3", "pp. 383-385",
        PLANE + "Bay 8 m, three storeys of 5 m: bases 7 (0, 0) and 8 (8, 0) fixed; joints 5, 6 at 5 m, 3, 4 at 10 m, "
        "1, 2 at 15 m. Beams I = 200e6 mm^4, columns I = 50e6 mm^4, E = 200,000 MPa. Loads at each of joints 1, 3, 5: "
        "3 kN along +x and a moment of +2 kN m (the source's load vector entries are 3 and 2 x 10^3, positive). The "
        "source neglects axial deformation. " + B.RIGID_NOTE + " First order, no shear deformation. " + UNUSED,
        model("ex-13-3-three-storey-frame", ("kN", "mm"), "Y", n133, {"n7": FIXED, "n8": FIXED}, SI,
              {"col": sec(AR, 50e6), "beam": sec(AR, 200e6)}, m133, lc, cb),
        [P("11.36", disp("n3", "ux"), "mm"),
         P("3.141e-4", disp("n3", "rz"), "rad", sign=-1, note="printed -3.141 x 10^-4", signnote=MOM),
         P("3.652e-4", disp("n4", "rz"), "rad", sign=-1, note="printed -3.652 x 10^-4", signnote=MOM),
         P("14.05", disp("n1", "ux"), "mm"),
         P("0.573e-4", disp("n1", "rz"), "rad", sign=-1, note="printed -0.573 x 10^-4", signnote=MOM),
         P("1.570e-4", disp("n2", "rz"), "rad", sign=-1, note="printed -1.570 x 10^-4", signnote=MOM),
         P("6.04", disp("n5", "ux"), "mm"),
         P("5.070e-4", disp("n5", "rz"), "rad", sign=-1, note="printed -5.070 x 10^-4", signnote=MOM),
         P("5.759e-4", disp("n6", "rz"), "rad", sign=-1, note="printed -5.759 x 10^-4", signnote=MOM)])

    # ------------------------------------------------------------------------------------- Ex. 13.5
    c30, s30 = math.cos(math.radians(30.0)), math.sin(math.radians(30.0))

    def rot(x, y):                                           # components in the joint axes x'', y'' of b
        return (x * c30 + y * s30, -x * s30 + y * c30)

    xa = 4000.0 / math.tan(math.radians(30.0))
    lc, cb = one_combo([nodal("a", Fx=rot(383.0, 321.4)[0], Fy=rot(383.0, 321.4)[1])])
    add("ex-13-5-truss-inclined-roller", "Three-bar plane truss on an inclined roller",
        "Example 13.5 (truss of Examples 3.1 and 3.2)", "pp. 390-391",
        TRUSS + "Truss of Example 3.2 (c (0, 0), a (4/tan 30 deg, 4), b (4/tan 30 deg + 4, 0) m; A_ac = 15e3, "
        "A_ab = 20e3, A_bc = 18e3 mm^2; E = 200,000 MPa; load at a 383.0 kN in x and 321.4 kN in y), c pinned, b on a "
        "roller whose plane is inclined 30 deg to the x axis. The input format has supports along global axes only. The "
        "model is therefore written in the source's joint axes x'', y'' of b: every coordinate and the load are rotated "
        "by -30 deg, so that the roller plane is the model's X axis and the roller holds uy at b. The rotation of the "
        "data is exact; the structure is the same. Only the displacement the source gives in these axes (u_b'') can be "
        "compared without transforming printed values; the source's u_a = 0.928 mm and v_a = 1.143 mm (global axes) are "
        "not targets (the engine gives 0.9282 and 1.1425 mm after rotating its result back, see README.md). First order.",
        model("ex-13-5-truss-inclined-roller", ("kN", "mm"), "Y",
              {"c": (0, 0), "b": rot(xa + 4000.0, 0.0), "a": rot(xa, 4000.0)},
              {"c": PIN2D, "b": sup("uy", "uz"), "a": sup("uz")}, SI,
              {"ac": bar(15e3), "ab": bar(20e3), "bc": bar(18e3)},
              {"ab": truss("a", "b", "ab"), "bc": truss("b", "c", "bc"), "ac": truss("a", "c", "ac")}, lc, cb),
        [P("0.094", disp("b", "ux"), "mm", sign=-1, note="printed -0.094 (u_b'' along the inclined plane)",
           signnote="the model's X axis is the source's joint axis x'' at b")])

    # ------------------------------------------------------------------------------------- Ex. 13.10
    lc, cb = one_combo([{"type": "memberPoint", "member": "ab", "direction": "Y", "value": -20.0, "position": 2000.0,
                         "positionType": "absolute"},
                        {"type": "memberPoint", "member": "bo", "direction": "Y", "value": -20.0, "position": 2000.0,
                         "positionType": "absolute"}])
    add("ex-13-10-three-span-beam", "Three-span beam with fixed ends and two point loads",
        "Example 13.10", "pp. 406-407",
        PLANE + "Beam of 15 m in three 5 m spans: a (0) and d (15 m) fixed, b (5 m) and c (10 m) on simple supports; "
        "I = 50e6 mm^4, E = 200,000 MPa. Loads: 20 kN downward at x = 2 m and at x = 7 m. A node o is placed at midlength "
        "(7.5 m), where the source reports the deflection. The full beam is analysed (the source combines a symmetric "
        "and an antisymmetric half-model). First order, no shear deformation. " + UNUSED,
        model("ex-13-10-three-span-beam", ("kN", "mm"), "Y",
              {"a": (0, 0), "b": (5000, 0), "o": (7500, 0), "c": (10000, 0), "d": (15000, 0)},
              {"a": FIXED, "d": FIXED, "b": sup("uy"), "c": sup("uy")}, SI, {"s": sec(1.0e4, 50e6)},
              {"ab": frame("a", "b", "s"), "bo": frame("b", "o", "s"), "oc": frame("o", "c", "s"),
               "cd": frame("c", "d", "s")}, lc, cb),
        [P("0.480e-3", disp("b", "rz"), "rad", sign=-1, note="printed -0.480 x 10^-3", signnote=MOM),
         P("0.720e-3", disp("c", "rz"), "rad", note="printed +0.720 x 10^-3", signnote=MOM),
         P("1.92", disp("o", "uy"), "mm", sign=-1, note="printed -1.92"),
         P("0.140e-3", disp("o", "rz"), "rad", note="printed 0.140 x 10^-3 (antisymmetric part and figure label)",
           signnote=MOM)])

    # ------------------------------------------------------------------------------------- Ex. 13.11
    run = 1000.0 * math.sqrt(55.0)
    lc, cb = one_combo([nodal("b", Fx=5.0), nodal("d", Fx=5.0), nodal("c", Fy=-10.0)])
    add("ex-13-11-gable-frame", "Symmetric gable frame under lateral and ridge loads",
        "Example 13.11", "pp. 407-408",
        PLANE + "a (0, 0), b (0, 5), c (sqrt 55, 8), d (2 sqrt 55, 5), e (2 sqrt 55, 0) m (rafters 8 m long with a 3 m "
        "rise; printed run 7.416 m); a and e fixed. Columns ab, ed: A = 4e3 mm^2, I = 50e6 mm^4; rafters bc, cd: "
        "A = 6e3 mm^2, I = 200e6 mm^4; E = 200,000 MPa. Loads: 5 kN along +x at b and at d, 10 kN downward at c. The "
        "full frame is analysed (the source combines a symmetric and an antisymmetric half-model). First order, axial "
        "deformation included. " + UNUSED,
        model("ex-13-11-gable-frame", ("kN", "mm"), "Y",
              {"a": (0, 0), "b": (0, 5000), "c": (run, 8000), "d": (2 * run, 5000), "e": (2 * run, 0)},
              {"a": FIXED, "e": FIXED}, SI, {"col": sec(4e3, 50e6), "raf": sec(6e3, 200e6)},
              {"ab": frame("a", "b", "col"), "bc": frame("b", "c", "raf"), "cd": frame("c", "d", "raf"),
               "ed": frame("e", "d", "col")}, lc, cb),
        [P("4.88", disp("b", "ux"), "mm"), P("0.02", disp("b", "uy"), "mm", sign=-1, note="printed -0.02"),
         P("1.351e-3", disp("b", "rz"), "rad", sign=-1, note="printed -1.351 x 10^-3", signnote=MOM),
         P("7.09", disp("c", "ux"), "mm"), P("5.58", disp("c", "uy"), "mm", sign=-1, note="printed -5.58"),
         P("0.368e-3", disp("c", "rz"), "rad", note="printed 0.368 x 10^-3", signnote=MOM),
         P("9.28", disp("d", "ux"), "mm"), P("0.04", disp("d", "uy"), "mm", sign=-1, note="printed -0.04"),
         P("0.129e-3", disp("d", "rz"), "rad", sign=-1, note="printed -0.129 x 10^-3", signnote=MOM)])

    # ------------------------------------------------------------------------------------- Ex. 13.12
    def chain(pid, title, loc, ks, den, mat, text):
        nodes = {f"n{k}": (1000.0 * (k - 1), 0) for k in range(1, 6)}
        sups = {"n1": PIN2D, "n5": PIN2D, "n2": sup("uy", "uz"), "n3": sup("uy", "uz"), "n4": sup("uy", "uz")}
        secs = {f"s{k}": bar(5.0 * kv) for k, kv in enumerate(ks, 1)}
        mem = {f"m{k}": truss(f"n{k}", f"n{k + 1}", f"s{k}") for k in range(1, 5)}
        lcs_ = {f"P{k}": {"loads": [nodal(f"n{k}", Fx=1.0)]} for k in (2, 3, 4)}
        cbs_ = {f"P{k}": {"factors": {f"P{k}": 1.0}} for k in (2, 3, 4)}
        tg = []
        for r_, row in zip((2, 3, 4), mat):
            for c_, num in zip((2, 3, 4), row):
                t = disp(f"n{r_}", "ux", f"P{c_}")
                t.update({"value": num / den, "unit": "mm", "tolerance": {"rel": 1.0e-6},
                          "precisionNote": f"exact rational number printed in the source ({num}/{den}); tolerance = "
                                           "policy floor 1e-6",
                          "sign": "displacement and load along +x"})
                tg.append(t)
        add(pid, title, loc, "pp. 411-412",
            TRUSS + "Chain of four axial springs along x between fixed ends 1 and 5, free nodes 2, 3, 4. The source's "
            "stiffnesses are dimensionless; here kN/mm: " + text + " Each spring is a truss member of 1000 mm with "
            "E = 200 kN/mm^2 and A = 5 k mm^2 (EA/L = k). The transverse translations of the collinear joints are held. "
            "Three combinations with a unit load (1 kN) at node 2, 3, 4: the displacements are the flexibility "
            "coefficients the source prints as exact fractions. First order.",
            model(pid, ("kN", "mm"), "Y", nodes, sups, SI, secs, mem, lcs_, cbs_), tg)

    chain("ex-13-12-spring-chain-original", "Chain of axial springs: flexibility matrix of the original structure",
          "Example 13.12, Equation 13.1 of the example (flexibility of the original structure)",
          (1.0, 2.0, 3.0, 2.0), 28, [[16, 10, 6], [10, 15, 9], [6, 9, 11]], "k = 1, 2, 3, 2 (members 1-2 ... 4-5).")
    chain("ex-13-12-spring-chain-changed", "Chain of axial springs: flexibility matrix of the changed structure",
          "Example 13.12, summary result (changed structure)",
          (2.0, 2.0, 3.0, 2.0), 44, [[16, 10, 6], [10, 20, 12], [6, 12, 16]],
          "k = 2, 2, 3, 2 (member 1-2 changed from 1 to 2).")

    # ------------------------------------------------------------------------- Eq. 13.31 / 13.32 / 13.33
    E, I, L = 200.0, 200e6, 6000.0
    a1, a2 = 3.0, 6.0
    eil = E * I / L
    al = a1 * a2 / (a1 * a2 + 4 * a1 + 4 * a2 + 12)
    k11 = al * eil * (12 / L ** 2) * (1 + (a1 + a2) / (a1 * a2))
    k12 = al * eil * (6 / L) * (1 + 2 / a2)
    k14 = al * eil * (6 / L) * (1 + 2 / a1)
    k22 = al * eil * 4 * (1 + 3 / a2)
    k24 = al * eil * 2
    th = 1.0e-3
    f31, f32, f33 = "Equation 13.31", "Equation 13.32", "Equation 13.33"
    nodes = {"n1": (0, 0), "n2": (L, 0), "n3": (0, 3000), "n4": (L, 3000), "n5": (0, 6000), "n6": (L, 6000)}
    lcs = {"V": {"loads": [settle(n, uy=1.0) for n in ("n1", "n3", "n5")]},
           "R": {"loads": [settle(n, rz=th) for n in ("n1", "n3", "n5")]}}
    cbs = {"V": {"factors": {"V": 1.0}}, "R": {"factors": {"R": 1.0}}}
    add("eq-13-31-flexible-connections", "Beam with flexible (semi-rigid) and pinned end connections: stiffness coefficients",
        "Equations 13.31, 13.32, 13.33 (Examples 13.6 and 13.7)", "pp. 394-398",
        PLANE + "The source gives the stiffness matrix of a beam with rotational springs k1 = alpha1 EI/L and "
        "k2 = alpha2 EI/L at its ends (Eq. 13.31) and the special cases of a pin at end 1 (Eq. 13.32) and at end 2 "
        "(Eq. 13.33). Numeric instance chosen here: L = 6 m, E = 200 kN/mm^2, I = 200e6 mm^4, alpha1 = 3, alpha2 = 6. "
        "Three separate members between fully held nodes: 'semi' (n1-n2) with partial rz releases (springs k1 at end i, "
        "k2 at end j), 'pin1' (n3-n4) with a full rz release at end i, 'pin2' (n5-n6) with a full rz release at end j. "
        "Combination V imposes v = 1 mm at the first node of each member, combination R a rotation of 0.001 rad there; "
        "the reactions are the first and second columns of the printed matrices times the imposed displacement (zero "
        "coefficients are not listed as targets). First order, no shear deformation. " + UNUSED,
        model("eq-13-31-flexible-connections", ("kN", "mm"), "Y", nodes, {k: FIXED for k in nodes}, SI,
              {"s": sec(1.0e4, I)},
              {"semi": frame("n1", "n2", "s", releases={"i": {"rz": {"spring": a1 * eil}},
                                                           "j": {"rz": {"spring": a2 * eil}}}),
               "pin1": frame("n3", "n4", "s", releases={"i": {"rz": True}}),
               "pin2": frame("n5", "n6", "s", releases={"j": {"rz": True}})}, lcs, cbs),
        stiff_targets([("V", "n1", "Fy", k11, f31), ("V", "n1", "Mz", k12, f31), ("V", "n2", "Fy", -k11, f31),
                       ("V", "n2", "Mz", k14, f31), ("R", "n1", "Fy", k12 * th, f31), ("R", "n1", "Mz", k22 * th, f31),
                       ("R", "n2", "Fy", -k12 * th, f31), ("R", "n2", "Mz", k24 * th, f31),
                       ("V", "n3", "Fy", 3 * eil / L ** 2, f32), ("V", "n4", "Fy", -3 * eil / L ** 2, f32),
                       ("V", "n4", "Mz", 3 * eil / L, f32),
                       ("V", "n5", "Fy", 3 * eil / L ** 2, f33), ("V", "n5", "Mz", 3 * eil / L, f33),
                       ("V", "n6", "Fy", -3 * eil / L ** 2, f33), ("R", "n5", "Fy", 3 * eil / L * th, f33),
                       ("R", "n5", "Mz", 3 * eil * th, f33), ("R", "n6", "Fy", -3 * eil / L * th, f33)], "n1"),
        basis="recomputed-from-published-statement")

    # ------------------------------------------------------------------------------------- Eq. 13.38-13.40
    E, A, I = 200.0, 6.0e3, 200e6
    p1, p2 = np.array([0.0, 0.0]), np.array([5000.0, 3000.0])
    ei, ej = np.array([400.0, 250.0]), np.array([300.0, 350.0])        # source: vector 1 -> i, vector j -> 2
    pi_, pj_ = p1 + ei, p2 - ej
    Le = float(np.linalg.norm(pj_ - pi_))
    cs, sn = (pj_ - pi_) / Le
    kl = np.zeros((6, 6))                                               # Eq. 4.34, plane part [u v th | u v th]
    ea, b12, b6, b4, b2 = E * A / Le, 12 * E * I / Le ** 3, 6 * E * I / Le ** 2, 4 * E * I / Le, 2 * E * I / Le
    kl[0, 0] = kl[3, 3] = ea
    kl[0, 3] = kl[3, 0] = -ea
    for (r_, c_, v) in [(1, 1, b12), (1, 2, b6), (1, 4, -b12), (1, 5, b6), (2, 2, b4), (2, 4, -b6), (2, 5, b2),
                        (4, 4, b12), (4, 5, -b6), (5, 5, b4)]:
        kl[r_, c_] = kl[c_, r_] = v
    g = np.array([[cs, sn, 0.0], [-sn, cs, 0.0], [0.0, 0.0, 1.0]])     # Example 5.2
    G = np.zeros((6, 6))
    G[:3, :3] = G[3:, 3:] = g
    T = np.eye(6)                                                       # Eq. 13.38
    T[2, 0], T[2, 1] = -ei[1], ei[0]
    T[5, 3], T[5, 4] = ej[1], -ej[0]
    K12 = T @ G.T @ kl @ G @ T.T                                        # Eq. 13.40
    f38 = "Equation 13.40 with T of Equation 13.38, the plane part of Equation 4.34 and the rotation matrix of Example 5.2"
    unit = {"U": (3, "ux", 1.0), "V": (4, "uy", 1.0), "R": (5, "rz", 1.0e-3)}
    lcs = {k: {"loads": [settle("n2", **{d: v})]} for k, (_, d, v) in unit.items()}
    cbs = {k: {"factors": {k: 1.0}} for k in unit}
    cases = []
    for k, (col, _, v) in unit.items():
        for row, (n, comp) in enumerate([("n1", "Fx"), ("n1", "Fy"), ("n1", "Mz"), ("n2", "Fx"), ("n2", "Fy"),
                                         ("n2", "Mz")]):
            cases.append((k, n, comp, float(K12[row, col] * v), f38))
    add("eq-13-38-rigid-end-offsets", "Member with rigid end offsets (finite joint size): stiffness coefficients",
        "Equations 13.37 to 13.40 (Example 13.8)", "pp. 398-401",
        PLANE + "The source gives the transformation of the stiffness of a plane member whose ends i, j are joined to "
        "the nodes 1, 2 by rigid links: F12 = T Gamma^T k' Gamma T^T Delta12, with offsets e_i = vector from node 1 to "
        "end i and e_j = vector from end j to node 2. Numeric instance chosen here: nodes 1 (0, 0) and 2 (5, 3) m; "
        "e_i = (400, 250) mm, e_j = (300, 350) mm, i.e. member offsets (400, 250, 0) at end i and (-300, -350, 0) at end "
        "j in the engine's convention (vector from the node to the flexible end); A = 6e3 mm^2, I = 200e6 mm^4, "
        "E = 200 kN/mm^2. Both nodes fully held; combinations U, V, R impose u2 = 1 mm, v2 = 1 mm, theta2 = 0.001 rad; "
        "the reactions at both nodes are columns 4, 5, 6 of the 6 x 6 matrix of Eq. 13.40 times the imposed "
        "displacement. First order, no shear deformation. " + UNUSED,
        model("eq-13-38-rigid-end-offsets", ("kN", "mm"), "Y", {"n1": tuple(p1), "n2": tuple(p2)},
              {"n1": FIXED, "n2": FIXED}, SI, {"s": sec(A, I)},
              {"m": frame("n1", "n2", "s", offsets={"axes": "global", "i": [400.0, 250.0, 0.0],
                                                     "j": [-300.0, -350.0, 0.0]})}, lcs, cbs),
        stiff_targets(cases, "n2"), basis="recomputed-from-published-statement")
