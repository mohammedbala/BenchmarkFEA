"""Problem definitions, chapters 2 to 5 (trusses, beams, plane frames, space truss, grid, member loads,
support settlement, lack of fit).  Used by build_problems.py; see the notes at the top of that file.

Every number below was read on the cited pages of the source.  Coordinates are derived from the printed
dimensions and angles (stated per problem); lengths in mm, forces in kN (the source's units; its
E = 200,000 MPa = 200 kN/mm^2).
"""
import math

S2 = math.sqrt(2.0)
S3 = math.sqrt(3.0)


def define(B):                                               # B = the build_problems module
    P, add, model = B.P, B.add, B.model
    disp, reac, station, endf = B.disp, B.reac, B.station, B.endf
    bar, sec, truss, frame, sup, nodal, settle = B.bar, B.sec, B.truss, B.frame, B.sup, B.nodal, B.settle
    one_combo, SI, FIXED, PIN2D = B.one_combo, B.STEEL_SI, B.FIXED, B.PIN2D
    PLANE, TRUSS, UNUSED, BAR_SIGN = B.PLANE_NOTE, B.TRUSS_NOTE, B.UNUSED, B.BAR_SIGN
    FREE2D = sup("uz")                                        # free joint of a plane truss

    # ------------------------------------------------------------------------------------- Ex. 2.1
    lc, cb = one_combo([nodal("a", Fx=500.0)])
    add("ex-2-1-truss-two-bar", "Two-bar plane truss: joint displacement under a horizontal load",
        "Example 2.1", "pp. 21-22",
        TRUSS + "Joints (m): b (0, 0), a (6, 4), c (9, 0) from the printed dimensions 6 m, 3 m, 4 m; b and c "
        "pinned. Bars ab A = 6000 mm^2, ac A = 8000 mm^2, E = 200,000 MPa. Load 500 kN along +x at a. First order.",
        model("ex-2-1-truss-two-bar", ("kN", "mm"), "Y", {"b": (0, 0), "a": (6000, 4000), "c": (9000, 0)},
              {"b": PIN2D, "c": PIN2D, "a": FREE2D}, SI, {"ab": bar(6000.0), "ac": bar(8000.0)},
              {"ab": truss("a", "b", "ab"), "ac": truss("a", "c", "ac")}, lc, cb),
        [P("2.41", disp("a", "ux"), "mm"), P("0.72", disp("a", "uy"), "mm"),
         P("400.6", station("ab", "N", rel=0.5), "kN", signnote=BAR_SIGN,
           note="bar force given with the problem figure (statics)"),
         P("277.8", station("ac", "N", rel=0.5), "kN", sign=-1, signnote=BAR_SIGN,
           note="bar force given with the problem figure (statics), printed -277.8")])

    # ------------------------------------------------------------------------------------- Ex. 2.2
    lc, cb = one_combo([nodal("a", Fx=1000.0)])
    add("ex-2-2-truss-four-bar", "Statically indeterminate four-bar plane truss",
        "Example 2.2", "pp. 22-23",
        TRUSS + "Joint a at the origin; bars to pinned supports b (+45 deg), c (-45 deg), d (210 deg), e (150 deg); "
        "the supports are 3 m above / below a, so b (3, 3), c (3, -3), d (-3 sqrt 3, -3), e (-3 sqrt 3, 3) m. "
        "A = 3000 mm^2 (ab, ac), 5000 mm^2 (ad, ae); E = 200,000 MPa. Load 1000 kN along +x at a. First order.",
        model("ex-2-2-truss-four-bar", ("kN", "mm"), "Y",
              {"a": (0, 0), "b": (3000, 3000), "c": (3000, -3000), "d": (-3000 * S3, -3000), "e": (-3000 * S3, 3000)},
              {"b": PIN2D, "c": PIN2D, "d": PIN2D, "e": PIN2D, "a": FREE2D}, SI,
              {"s3": bar(3000.0), "s5": bar(5000.0)},
              {"ab": truss("a", "b", "s3"), "ac": truss("a", "c", "s3"), "ad": truss("a", "d", "s5"),
               "ae": truss("a", "e", "s5")}, lc, cb),
        [P("2.55", disp("a", "ux"), "mm"),
         P("368.7", station("ad", "N", rel=0.5), "kN", signnote=BAR_SIGN),
         P("368.7", station("ae", "N", rel=0.5), "kN", signnote=BAR_SIGN),
         P("255.4", station("ab", "N", rel=0.5), "kN", sign=-1, signnote=BAR_SIGN, note="printed -255.4"),
         P("255.4", station("ac", "N", rel=0.5), "kN", sign=-1, signnote=BAR_SIGN, note="printed -255.4")])

    # ------------------------------------------------------------------------------------- Ex. 2.3
    lc, cb = one_combo([settle("a", uy=-5.0)])
    add("ex-2-3-truss-prescribed-displacement", "Force needed for a prescribed joint displacement of a two-bar truss",
        "Example 2.3", "pp. 23-24",
        TRUSS + "Joint a at the origin, pinned supports b (2, 2) m and c (-4, 2) m (printed dimensions 4 m, 2 m, 2 m). "
        "A_ab = 1000 mm^2, A_ac = 4000 mm^2, E = 200,000 MPa. Joint a is moved 5 mm vertically downward with no "
        "horizontal displacement: modelled as a support at a (ux, uy) with the imposed displacement uy = -5 mm; the "
        "force the source asks for is the reaction of that support. First order.",
        model("ex-2-3-truss-prescribed-displacement", ("kN", "mm"), "Y",
              {"a": (0, 0), "b": (2000, 2000), "c": (-4000, 2000)},
              {"a": PIN2D, "b": PIN2D, "c": PIN2D}, SI, {"ab": bar(1000.0), "ac": bar(4000.0)},
              {"ab": truss("a", "b", "ab"), "ac": truss("a", "c", "ac")}, lc, cb),
        [P("181.0", reac("a", "Fx"), "kN",
           signnote="the source's required force P on joint a = the engine's reaction at a (force on the structure)"),
         P("355.7", reac("a", "Fy"), "kN", sign=-1, note="printed -355.7",
           signnote="the source's required force P on joint a = the engine's reaction at a (force on the structure)")])

    # --------------------------------------------------------------------------- Ex. 3.2 / 3.3 / 5.10
    xa = 4000.0 / math.tan(math.radians(30.0))
    n3 = {"c": (0, 0), "b": (xa + 4000.0, 0), "a": (xa, 4000.0)}
    geo3 = ("Joints from the printed geometry (height 4 m, bar ac at 30 deg, bar ab at 45 deg to the horizontal; the "
            "source prints the resulting base length as 10.928 m): c (0, 0), a (4/tan 30 deg, 4), b (4/tan 30 deg + 4, 0) "
            "m. A_ac = 15e3, A_ab = 20e3, A_bc = 18e3 mm^2, E = 200,000 MPa. ")
    s3 = {"ac": bar(15e3), "ab": bar(20e3), "bc": bar(18e3)}
    m3 = {"ab": truss("a", "b", "ab"), "bc": truss("b", "c", "bc"), "ac": truss("a", "c", "ac")}
    lc, cb = one_combo([nodal("a", Fx=383.0, Fy=321.4)])
    add("ex-3-2-truss-three-bar", "Three-bar plane truss: displacements, reactions, bar forces",
        "Example 3.2 (truss of Example 3.1)", "pp. 35-38",
        TRUSS + geo3 + "c pinned, b on a horizontal roller (vertical restraint). Load at a: the source's components "
        "383.0 kN (+x) and 321.4 kN (+y) of 500 kN at 40 deg. First order.",
        model("ex-3-2-truss-three-bar", ("kN", "mm"), "Y", n3,
              {"c": PIN2D, "b": sup("uy", "uz"), "a": FREE2D}, SI, s3, m3, lc, cb),
        [P("0.871", disp("a", "ux"), "mm"), P("1.244", disp("a", "uy"), "mm"),
         P("0.193", disp("b", "ux"), "mm", sign=-1, note="printed -0.193"),
         P("63.6", reac("b", "Fy"), "kN", sign=-1, note="printed -63.6"),
         P("383.4", reac("c", "Fx"), "kN", sign=-1, note="printed -383.4"),
         P("258.1", reac("c", "Fy"), "kN", sign=-1, note="printed -258.1"),
         P("90.0", station("ab", "N", rel=0.5), "kN", signnote=BAR_SIGN),
         P("63.6", station("bc", "N", rel=0.5), "kN", sign=-1, signnote=BAR_SIGN, note="printed -63.6"),
         P("516.2", station("ac", "N", rel=0.5), "kN", signnote=BAR_SIGN)])

    n33 = dict(n3)
    n33["d"] = (0, 4000.0)
    s33 = dict(s3)
    s33["ad"] = bar(20e3)
    m33 = dict(m3)
    m33["ad"] = truss("d", "a", "ad")
    add("ex-3-3-truss-with-tie", "Three-bar plane truss with an added horizontal tie",
        "Example 3.3", "pp. 38-39",
        TRUSS + geo3 + "Added: horizontal tie ad, A = 20e3 mm^2 (the figure prints the unit as 'mm'), length "
        "6.928 m, d held by a wall. Supports and load as Example 3.2. First order.",
        model("ex-3-3-truss-with-tie", ("kN", "mm"), "Y", n33,
              {"c": PIN2D, "b": sup("uy", "uz"), "d": PIN2D, "a": FREE2D}, SI, s33, m33, lc, cb),
        [P("0.383", disp("a", "ux"), "mm"), P("1.228", disp("a", "uy"), "mm"),
         P("0.437", disp("b", "ux"), "mm", sign=-1, note="printed -0.437"),
         P("144.3", reac("b", "Fy"), "kN", sign=-1, note="printed -144.3"),
         P("163.2", reac("c", "Fx"), "kN", sign=-1, note="printed -163.2"),
         P("177.4", reac("c", "Fy"), "kN", sign=-1, note="printed -177.4"),
         P("221.1", reac("d", "Fx"), "kN", sign=-1, note="printed -221.1")])

    lc10, cb10 = one_combo([settle("d", ux=-3.0)])
    add("ex-5-10-truss-lack-of-fit", "Plane truss with a bar fabricated 3 mm too short",
        "Example 5.10 (truss of Example 3.3)", "pp. 123-124",
        TRUSS + geo3 + "Tie ad (A = 20e3 mm^2) is fabricated 3 mm too short; no external load. The input format has "
        "no initial-strain load. The lack of fit is entered as an exactly equivalent loading: the full-length bar with "
        "its wall support d moved 3 mm away from a (imposed ux = -3 mm at d). The force in bar ad is "
        "(EA/L)(u_a - u_d + 3 mm) with u_d = 0 in the source and (EA/L)(u_a - (-3 mm)) here: the same equation, hence "
        "the same joint displacements and reactions. First order.",
        model("ex-5-10-truss-lack-of-fit", ("kN", "mm"), "Y", n33,
              {"c": PIN2D, "b": sup("uy", "uz"), "d": PIN2D, "a": FREE2D}, SI, s33, m33, lc10, cb10),
        [P("1.683", disp("a", "ux"), "mm", sign=-1, note="printed -1.683"),
         P("0.0518", disp("a", "uy"), "mm", sign=-1, note="printed -0.0518"),
         P("0.8448", disp("b", "ux"), "mm", sign=-1, note="printed -0.8448"),
         P("278.1", reac("b", "Fy"), "kN", sign=-1, note="printed -278.1"),
         P("760.0", reac("c", "Fx"), "kN"),
         P("278.1", reac("c", "Fy"), "kN"),
         P("760.2", reac("d", "Fx"), "kN", sign=-1, note="printed -760.2")])

    # ------------------------------------------------------------------------------- Ex. 3.5 / 3.6
    h = 5000.0 * S3
    n5 = {"a": (0, 0), "b": (5000, 0), "c": (5000, h), "d": (10000, h)}
    geo5 = ("Joints from the printed geometry (5 m base, inclined bars at 60 deg): a (0, 0), b (5, 0), c (5, 5 sqrt 3), "
            "d (10, 5 sqrt 3) m. A = 10e3 mm^2 (ab, cd), 15e3 mm^2 (ac, bd, bc); E = 200,000 MPa. Load at d: the "
            "source's components 282.84 kN (+x) and -282.84 kN (y) of 400 kN at 45 deg. ")
    s5 = {"s10": bar(10e3), "s15": bar(15e3)}
    m5 = {"ab": truss("a", "b", "s10"), "cd": truss("c", "d", "s10"), "ac": truss("a", "c", "s15"),
          "bd": truss("b", "d", "s15"), "bc": truss("b", "c", "s15")}
    lc, cb = one_combo([nodal("d", Fx=282.84, Fy=-282.84)])
    add("ex-3-5-truss-five-bar", "Five-bar plane truss: displacements, reactions, bar forces",
        "Example 3.5 (truss of Example 3.4)", "pp. 42-46",
        TRUSS + geo5 + "a pinned, b on a horizontal roller. First order.",
        model("ex-3-5-truss-five-bar", ("kN", "mm"), "Y", n5,
              {"a": PIN2D, "b": sup("uy", "uz"), "c": FREE2D, "d": FREE2D}, SI, s5, m5, lc, cb),
        [P("0.407", disp("b", "ux"), "mm", sign=-1, note="printed -0.407"),
         P("9.809", disp("c", "ux"), "mm"), P("2.232", disp("c", "uy"), "mm", sign=-1, note="printed -2.232"),
         P("10.926", disp("d", "ux"), "mm"), P("7.801", disp("d", "uy"), "mm", sign=-1, note="printed -7.801"),
         P("282.9", reac("a", "Fx"), "kN", sign=-1, note="printed -282.9"),
         P("772.0", reac("a", "Fy"), "kN", sign=-1, note="printed -772.0"),
         P("1056.2", reac("b", "Fy"), "kN"),
         P("162.8", station("ab", "N", rel=0.5), "kN", sign=-1, signnote=BAR_SIGN, note="printed -162.8"),
         P("446.8", station("cd", "N", rel=0.5), "kN", signnote=BAR_SIGN),
         P("891.4", station("ac", "N", rel=0.5), "kN", signnote=BAR_SIGN),
         P("326.8", station("bd", "N", rel=0.5), "kN", sign=-1, signnote=BAR_SIGN, note="printed -326.8"),
         P("773.2", station("bc", "N", rel=0.5), "kN", sign=-1, signnote=BAR_SIGN, note="printed -773.2")])
    add("ex-3-6-truss-added-restraints", "Five-bar plane truss with added horizontal restraints",
        "Example 3.6", "p. 46",
        TRUSS + geo5 + "As Example 3.5 with horizontal restraints added at b and c (a and b pinned, c held "
        "horizontally). First order.",
        model("ex-3-6-truss-added-restraints", ("kN", "mm"), "Y", n5,
              {"a": PIN2D, "b": PIN2D, "c": sup("ux", "uz"), "d": FREE2D}, SI, s5, m5, lc, cb),
        [P("1.114", disp("d", "ux"), "mm"), P("1.901", disp("d", "uy"), "mm", sign=-1, note="printed -1.901")])

    # ---------------------------------------------------------------------- beam of Ex. 4.8 (4.9-4.12)
    nb = {"a": (0, 0), "b": (8000, 0), "c": (13000, 0)}
    sb = {"ab": sec(6e3, 200e6, 300e3), "bc": sec(4e3, 50e6, 100e3)}
    mb = {"ab": frame("a", "b", "ab"), "bc": frame("b", "c", "bc")}
    beam48 = ("Beam of Example 4.8 along x: a (0), b (8 m), c (13 m). ab: A = 6e3 mm^2, Iz = 200e6 mm^4, "
              "J = 300e3 mm^4; bc: A = 4e3 mm^2, Iz = 50e6 mm^4, J = 100e3 mm^4; E = 200,000 MPa, nu = 0.3. The "
              "source omits bending normal to the plane: Iy is a placeholder (= Iz) that does not influence the targets. ")
    MOM = "moment / rotation about z, counter-clockwise positive in the source = rz / Mz of the engine"
    lc, cb = one_combo([nodal("c", Fx=5.0 / S2, Fy=-5.0 / S2)])
    add("ex-4-9-beam-overhang", "Propped cantilever with overhang, inclined end load",
        "Example 4.9", "pp. 80-81",
        PLANE + beam48 + "a fixed, b on a roller (vertical restraint), c free. Load at c: 5 kN at 45 deg, components "
        "+5/sqrt 2 (x) and -5/sqrt 2 (y) kN. First order, no shear deformation.",
        model("ex-4-9-beam-overhang", ("kN", "mm"), "Y", nb, {"a": FIXED, "b": sup("uy")}, SI, sb, mb, lc, cb),
        [P("0.024", disp("b", "ux"), "mm"), P("0.046", disp("c", "ux"), "mm"),
         P("19.15", disp("c", "uy"), "mm", sign=-1, note="printed -19.15"),
         P("0.00088", disp("b", "rz"), "rad", sign=-1, note="printed -0.00088", signnote=MOM),
         P("0.00530", disp("c", "rz"), "rad", sign=-1, note="printed -0.00530", signnote=MOM),
         P("3.60", reac("a", "Fx"), "kN", sign=-1, note="printed -3.60"),
         P("3.30", reac("a", "Fy"), "kN", sign=-1, note="printed -3.30"),
         P("8.80", reac("a", "Mz"), "kN mm", sign=-1, factor=1000.0, note="printed -8.80 kNm", signnote=MOM),
         P("6.85", reac("b", "Fy"), "kN"),
         P("17.68", station("ab", "Mz", rel=1.0), "kN mm", sign=-1, factor=1000.0,
           note="moment-diagram label 17.68 kNm at b (magnitude)",
           signnote="hogging moment at b (downward load on the overhang): negative Mz in the engine (y up)")])

    lc, cb = one_combo([nodal("b", Mz=50.0e3)])
    add("ex-4-10-beam-joint-moment", "Two-span beam with fixed ends, moment applied at the interior support",
        "Example 4.10", "pp. 81-82",
        PLANE + beam48 + "a and c fixed, b on a roller. Load: moment 50 kN m (counter-clockwise) at b. First order.",
        model("ex-4-10-beam-joint-moment", ("kN", "mm"), "Y", nb, {"a": FIXED, "b": sup("uy"), "c": FIXED}, SI, sb,
              mb, lc, cb),
        [P("0.001786", disp("b", "rz"), "rad", signnote=MOM),
         P("6.70", reac("a", "Fy"), "kN"),
         P("17.86", reac("a", "Mz"), "kN mm", factor=1000.0, signnote=MOM),
         P("2.41", reac("b", "Fy"), "kN", sign=-1, note="printed -2.41"),
         P("4.29", reac("c", "Fy"), "kN", sign=-1, note="printed -4.29"),
         P("7.14", reac("c", "Mz"), "kN mm", factor=1000.0, signnote=MOM),
         P("35.72", endf("ab", "j", "Mz"), "kN mm", factor=1000.0,
           signnote="member end moment on member ab at b, element convention of the source = endForces of the engine "
                    "(local z = global Z for these members)"),
         P("14.28", endf("bc", "i", "Mz"), "kN mm", factor=1000.0,
           signnote="member end moment on member bc at b, element convention of the source = endForces of the engine")])

    lc, cb = one_combo([settle("b", uy=-20.0)])
    add("ex-4-11-beam-settlement", "Two-span beam with fixed ends, settlement of the interior support",
        "Example 4.11", "pp. 82-83",
        PLANE + beam48 + "a and c fixed, b on a roller that settles 20 mm (imposed uy = -20 mm). No other load. "
        "First order.",
        model("ex-4-11-beam-settlement", ("kN", "mm"), "Y", nb, {"a": FIXED, "b": sup("uy"), "c": FIXED}, SI, sb,
              mb, lc, cb),
        [P("0.0009643", disp("b", "rz"), "rad", sign=-1, note="printed -0.0009643", signnote=MOM),
         P("15.14", reac("a", "Fy"), "kN"),
         P("65.36", reac("a", "Mz"), "kN mm", factor=1000.0, signnote=MOM),
         P("36.66", reac("b", "Fy"), "kN", sign=-1, note="printed -36.66"),
         P("21.51", reac("c", "Fy"), "kN"),
         P("51.86", reac("c", "Mz"), "kN mm", sign=-1, factor=1000.0, note="printed -51.86", signnote=MOM),
         P("55.71", endf("ab", "j", "Mz"), "kN mm", factor=1000.0,
           signnote="member end moment on member ab at b, element convention = endForces of the engine"),
         P("55.71", endf("bc", "i", "Mz"), "kN mm", sign=-1, factor=1000.0, note="printed -55.71",
           signnote="member end moment on member bc at b, element convention = endForces of the engine")])

    lc, cb = one_combo([nodal("b", Fy=-1.0, Mx=40.0)])
    add("ex-4-12-beam-bending-torsion", "Fixed-ended beam with an eccentric load: bending and torsion",
        "Example 4.12", "pp. 83-84",
        PLANE + beam48 + "a and c fixed, no support at b. A 1 kN downward force acts on a rigid bracket projecting "
        "40 mm from the beam axis at b in the +z direction: entered as the source does, Fy = -1 kN and a torque "
        "Mx = +40 kN mm at b. First order; Saint-Venant torsion.",
        model("ex-4-12-beam-bending-torsion", ("kN", "mm"), "Y", nb, {"a": FIXED, "c": FIXED}, SI, sb, mb, lc, cb),
        [P("0.545", disp("b", "uy"), "mm", sign=-1, note="printed -0.545"),
         P("0.009044", disp("b", "rx"), "rad"),
         P("0.0000263", disp("b", "rz"), "rad", sign=-1, note="printed -0.0000263", signnote=MOM),
         P("0.423", reac("a", "Fy"), "kN"),
         P("26.1", reac("a", "Mx"), "kN mm", sign=-1, note="printed -26.1 kNmm"),
         P("1.781", reac("a", "Mz"), "kN mm", factor=1000.0, signnote=MOM),
         P("0.586", reac("c", "Fy"), "kN"),
         P("13.9", reac("c", "Mx"), "kN mm", sign=-1, note="printed -13.9 kNmm"),
         P("1.413", reac("c", "Mz"), "kN mm", sign=-1, factor=1000.0, note="printed -1.413 kNm", signnote=MOM)])

    # ------------------------------------------------------------------------------------- Ex. 4.13
    lc, cb = one_combo([nodal("b", Fx=100.0 / S2, Fy=-100.0 / S2, Mz=50.0e3)])
    add("ex-4-13-frame-two-member", "Two-member rigid plane frame with a loaded joint",
        "Example 4.13", "pp. 84-85",
        PLANE + "Beam ab horizontal, 8 m (A = 6e3 mm^2, I = 200e6 mm^4); column bc vertical, 5 m, c below b "
        "(A = 4e3 mm^2, I = 50e6 mm^4): a (0, 0), b (8, 0), c (8, -5) m. a and c fixed. E = 200,000 MPa. Loads at b: "
        "100 kN at 45 deg (+100/sqrt 2 in x, -100/sqrt 2 in y) and a moment of 50 kN m, counter-clockwise (the "
        "solution uses 50 x 10^3 kN mm; the problem figure prints the unit as 'kN'). First order, axial deformation "
        "included. " + UNUSED,
        model("ex-4-13-frame-two-member", ("kN", "mm"), "Y", {"a": (0, 0), "b": (8000, 0), "c": (8000, -5000)},
              {"a": FIXED, "c": FIXED}, SI, {"ab": sec(6e3, 200e6), "bc": sec(4e3, 50e6)},
              {"ab": frame("a", "b", "ab"), "bc": frame("b", "c", "bc")}, lc, cb),
        [P("0.4414", disp("b", "ux"), "mm"), P("0.3998", disp("b", "uy"), "mm", sign=-1, note="printed -0.3998"),
         P("0.00169", disp("b", "rz"), "rad", signnote=MOM)])

    # ------------------------------------------------------------------------------------- Ex. 4.15
    US = B.STEEL_US
    nn = {f"n{k}": (15.0 * k, 0.0) for k in range(21)}
    ss = {"n0": sup("ux", "uy", "uz", "rx"), "n20": sup("uy", "uz")}
    mm = {f"s{k}": frame(f"n{k - 1}", f"n{k}", "beam") for k in range(1, 21)}
    for k in range(1, 20):
        nn[f"g{k}"] = (15.0 * k, -15.0)
        ss[f"g{k}"] = PIN2D
        mm[f"f{k}"] = truss(f"g{k}", f"n{k}", "spring")
    lc, cb = one_combo([nodal("n10", Fy=-40.0)])
    add("ex-4-15-beam-on-elastic-supports", "Beam on an elastic foundation modelled by 19 axial members",
        "Example 4.15, part 2 (stiffness method)", "pp. 87-88",
        PLANE + "Units kip, in. Beam ab of 300 in in 20 segments of 15 in, I = 128.5 in^4, E = 29,000 ksi, pin at a, "
        "roller at b. The foundation (modulus 1.5 kip/in^2) is represented as in the source by 19 vertical axial "
        "members, one under each interior node, L = 15 in, 22.5 kip/in each, pinned to the ground; their area is entered "
        "as the source's defining expression 22.5 x 15 / 29,000 in^2 (printed rounded as 0.0116 in^2). "
        "Load 40 kips downward at midspan c. The beam area is not given by the source and carries no force (A = 10 in^2 "
        "placeholder). First order. The targets are the source's values for this discrete model (not its closed-form "
        "continuous-foundation values 0.237 in and 554 in-kips, which the engine cannot represent).",
        model("ex-4-15-beam-on-elastic-supports", ("kip", "in"), "Y", nn, ss, US,
              {"beam": sec(10.0, 128.5), "spring": bar(22.5 * 15.0 / 29000.0)}, mm, lc, cb),
        [P("0.238", disp("n10", "uy"), "in", sign=-1, note="label on the deflected-shape figure (magnitude)",
           signnote="downward deflection: negative uy"),
         P("547", station("s10", "Mz", rel=1.0), "kip in", note="label on the moment-diagram figure (magnitude)",
           signnote="sagging moment under the load: positive Mz in the engine (y up)")])

    # ------------------------------------------------------------------------------------- Ex. 5.3
    run = 1000.0 * math.sqrt(55.0)                           # 8 m rafter with a 3 m rise
    ng = {"a": (0, 0), "b": (0, 5000), "c": (run, 8000), "d": (run + 4000.0, 5000), "e": (run + 4000.0, 0)}
    sg = {"col": sec(4e3, 50e6), "raf": sec(6e3, 200e6)}
    mg = {"ab": frame("a", "b", "col"), "bc": frame("b", "c", "raf"), "cd": frame("c", "d", "col"),
          "ed": frame("e", "d", "col")}
    unit = {"U": ("ux", 1.0), "V": ("uy", 1.0), "R": ("rz", 1.0e-3)}
    lcs, cbs = {}, {}
    for nd in "bcd":
        for k, (dof, val) in unit.items():
            lcs[k + nd] = {"loads": [settle(nd, **{dof: val})]}
            cbs[k + nd] = {"factors": {k + nd: 1.0}}
    printed53 = {"b": [("U", "Fx", "0.6500", 1), ("U", "Fy", "0.2591", 1), ("U", "Mz", "4.969", 1),
                       ("V", "Fy", "0.9095", 1), ("V", "Mz", "17.381", 1), ("R", "Mz", "1.4e5", 1)],
                 "c": [("U", "Fx", "1.1589", 1), ("U", "Fy", "0.1226", -1), ("U", "Mz", "14.231", 1),
                       ("V", "Fy", "0.4006", 1), ("V", "Mz", "7.781", -1), ("R", "Mz", "1.4e5", 1)],
                 "d": [("U", "Fx", "0.5185", 1), ("U", "Fy", "0.3817", -1), ("U", "Mz", "4.80", 1),
                       ("V", "Fy", "1.0911", 1), ("V", "Mz", "9.60", -1), ("R", "Mz", "0.8e5", 1)]}
    t53 = []
    for nd, rows in printed53.items():
        for k, comp, pr, sg_ in rows:
            fac = 200.0 * unit[k][1]
            t53.append(P(pr, reac(nd, comp, comb=k + nd), "kN" if comp != "Mz" else "kN mm", sign=sg_, factor=fac,
                         note=f"assembled global stiffness coefficient K({comp} at {nd}, {unit[k][0]} at {nd}), the "
                              "source prints the matrix with the common factor 200"
                              + ("; imposed rotation 0.001 rad" if k == "R" else "; imposed displacement 1 mm"),
                         signnote="stiffness coefficient = reaction at the node where the unit displacement is imposed, "
                                  "all other joint displacements held at zero; same axes as the source"))
    add("ex-5-3-gable-frame-stiffness", "Gable frame: assembled global stiffness coefficients (rotated members)",
        "Example 5.3 (assembled global stiffness equation)", "pp. 102-104",
        PLANE + "Gable frame, joints (m): a (0, 0), b (0, 5), c (sqrt 55, 8), d (sqrt 55 + 4, 5), e (sqrt 55 + 4, 0) "
        "- the source prints 7.416 m and 11.416 m for the 8 m rafter with a 3 m rise and the 5 m rafter with a 3 m "
        "drop. ab, cd, ed: A = 4e3 mm^2, I = 50e6 mm^4; bc: A = 6e3 mm^2, I = 200e6 mm^4; E = 200,000 MPa. The source "
        "assembles the stiffness matrix and solves nothing. A stiffness coefficient is obtained here as a reaction: "
        "every joint is fully held, one joint displacement (1 mm or 0.001 rad) is imposed per combination, and the "
        "reactions at that joint are the coefficients of the diagonal 3 x 3 block of that joint. First order, axial "
        "deformation included. " + UNUSED,
        model("ex-5-3-gable-frame-stiffness", ("kN", "mm"), "Y", ng, {k: FIXED for k in ng}, SI, sg, mg, lcs, cbs),
        t53)

    # ------------------------------------------------------------------------------------- Ex. 5.4
    lc, cb = one_combo([nodal("a", Fx=200.0, Fy=600.0, Fz=-800.0)])
    t54 = [P("0.1783", disp("a", "ux"), "mm"), P("2.722", disp("a", "uy"), "mm"),
           P("0.4863", disp("a", "uz"), "mm", sign=-1, note="printed -0.4863")]
    for nd, vals in {"b": ("-76.4", "-152.8", "-305.6"), "c": ("170.8", "-113.8", "-227.7"),
                     "d": ("-470.7", "-156.9", "627.8"), "e": ("176.3", "-176.3", "705.5")}.items():
        for comp, pr in zip(("Fx", "Fy", "Fz"), vals):
            neg = pr.startswith("-")
            t54.append(P(pr.lstrip("-"), reac(nd, comp), "kN", sign=-1 if neg else 1,
                         note=("printed " + pr) if neg else ""))
    add("ex-5-4-space-truss", "Pin-jointed space truss: four bars to one loaded joint",
        "Example 5.4", "pp. 104-106",
        "Space truss, global z vertical (model: verticalAxis Z, the source's axes are the model's). Joint a (2, 4, 8) m; "
        "pinned supports b (0, 0, 0), c (8, 0, 0), d (8, 6, 0), e (0, 6, 0) m - from the source's table of member "
        "projections (ab: -2, -4, -8; ac: 6, -4, -8; ad: 6, 2, -8; ae: -2, 2, -8). A = 20e3 (ab), 30e3 (ac), 40e3 (ad), "
        "30e3 (ae) mm^2; E = 200,000 MPa. Loads at a: 200 kN (x), 600 kN (y), -800 kN (z). Members of type 'truss'. "
        "First order.",
        model("ex-5-4-space-truss", ("kN", "mm"), "Z",
              {"a": (2000, 4000, 8000), "b": (0, 0, 0), "c": (8000, 0, 0), "d": (8000, 6000, 0), "e": (0, 6000, 0)},
              {k: PIN2D for k in "bcde"}, SI,
              {"ab": bar(20e3), "ac": bar(30e3), "ad": bar(40e3), "ae": bar(30e3)},
              {"ab": truss("a", "b", "ab"), "ac": truss("a", "c", "ac"), "ad": truss("a", "d", "ad"),
               "ae": truss("a", "e", "ae")}, lc, cb),
        t54)

    # ------------------------------------------------------------------------------- Ex. 5.6 / 5.11
    sup56 = {"a": sup("ux", "uy"), "b": sup("uy"), "c": FIXED}
    lc, cb = one_combo([{"type": "memberDistributed", "member": "ab", "direction": "Y", "w1": -0.002},
                        {"type": "memberPoint", "member": "bc", "direction": "Y", "value": -20.0, "position": 2000.0,
                         "positionType": "absolute"}])
    add("ex-5-6-beam-member-loads", "Two-span beam with a distributed load and a point load between nodes",
        "Example 5.6", "pp. 113-114",
        PLANE + beam48 + "a: pin (the source has no axial degree of freedom; ux is held at a), b: roller, c: fixed. "
        "Loads: 2 kN/m downward on ab; 20 kN downward on bc at 2 m from b. First order, no shear deformation.",
        model("ex-5-6-beam-member-loads", ("kN", "mm"), "Y", nb, sup56, SI, sb, mb, lc, cb),
        [P("5.684e-4", disp("a", "rz"), "rad", sign=-1, note="printed -5.684 x 10^-4", signnote=MOM),
         P("0.698e-4", disp("b", "rz"), "rad", note="printed 0.698 x 10^-4", signnote=MOM),
         P("6.13", reac("a", "Fy"), "kN"), P("23.00", reac("b", "Fy"), "kN"), P("6.87", reac("c", "Fy"), "kN"),
         P("9.32e3", reac("c", "Mz"), "kN mm", sign=-1, note="printed -9.32 x 10^3 kNmm", signnote=MOM),
         P("14.96e3", endf("ab", "j", "Mz"), "kN mm", sign=-1, note="printed -14.96 x 10^3 kNmm",
           signnote="member end moment on ab at b, element convention = endForces of the engine"),
         P("14.96e3", endf("bc", "i", "Mz"), "kN mm", note="printed 14.96 x 10^3 kNmm",
           signnote="member end moment on bc at b, element convention = endForces of the engine"),
         P("11.29", station("bc", "Mz", x=2000.0), "kN mm", factor=1000.0,
           note="moment-diagram label 11.29 (kNm) under the 20 kN load (magnitude)",
           signnote="sagging under the load: positive Mz in the engine (y up)")])

    lc, cb = one_combo([settle("b", uy=-15.0)])
    add("ex-5-11-beam-settlement", "Two-span beam: settlement of the interior support",
        "Example 5.11", "pp. 124-125",
        PLANE + beam48 + "Supports as Example 5.6 (a pin, b roller, c fixed). Support b settles 15 mm (imposed "
        "uy = -15 mm); no other load. First order.",
        model("ex-5-11-beam-settlement", ("kN", "mm"), "Y", nb, sup56, SI, sb, mb, lc, cb),
        [P("0.002984", disp("a", "rz"), "rad", sign=-1, note="printed -0.002984", signnote=MOM),
         P("0.000342", disp("b", "rz"), "rad", signnote=MOM),
         P("4.16", reac("a", "Fy"), "kN"), P("17.74", reac("b", "Fy"), "kN", sign=-1, note="printed -17.74"),
         P("13.58", reac("c", "Fy"), "kN"),
         P("34.63e3", reac("c", "Mz"), "kN mm", sign=-1, note="printed -34.63 x 10^3 kNmm", signnote=MOM)])

    # ------------------------------------------------------------------------------------- Ex. 5.7
    lc, cb = one_combo([nodal("b", Fx=18.75, Fy=-46.35),
                        {"type": "memberDistributed", "member": "bc", "direction": "Y", "w1": -0.004}])
    add("ex-5-7-frame-inclined-member", "Two-member rigid frame with an inclined member, joint load and distributed load",
        "Example 5.7", "pp. 115-116",
        PLANE + "a (0, 0), b (sqrt 55, 3), c (sqrt 55 + 8, 3) m: ab is 8 m long with a 3 m rise (printed run 7.416 m), "
        "bc horizontal, 8 m. Both members A = 6000 mm^2, I = 200e6 mm^4; E = 200,000 MPa. a and c fixed. Loads: the "
        "source's components of the 50 kN joint load at b, 18.75 kN (+x) and -46.35 kN (y); 4 kN/m downward on bc. "
        "First order, flexural and axial deformation. " + UNUSED,
        model("ex-5-7-frame-inclined-member", ("kN", "mm"), "Y",
              {"a": (0, 0), "b": (run, 3000), "c": (run + 8000.0, 3000)}, {"a": FIXED, "c": FIXED}, SI,
              {"s": sec(6000.0, 200e6)}, {"ab": frame("a", "b", "s"), "bc": frame("b", "c", "s")}, lc, cb),
        [P("0.9982", disp("b", "ux"), "mm",
           note="the value used in the source's reaction calculation (p. 116); the solution line on p. 115 prints 0.09982"),
         P("4.996", disp("b", "uy"), "mm", sign=-1, note="printed -4.996"),
         P("0.000534", disp("b", "rz"), "rad", sign=-1, note="printed -0.000534", signnote=MOM),
         P("131.0", reac("a", "Fx"), "kN"), P("55.4", reac("a", "Fy"), "kN"),
         P("13.43e3", reac("a", "Mz"), "kN mm", note="printed 13.43 x 10^3 kNmm", signnote=MOM),
         P("149.7", reac("c", "Fx"), "kN", sign=-1, note="printed -149.7"),
         P("22.7", reac("c", "Fy"), "kN"),
         P("45.41e3", reac("c", "Mz"), "kN mm", sign=-1, note="printed -45.41 x 10^3 kNmm", signnote=MOM)])

    # ------------------------------------------------------------------------------------- Ex. 5.8
    lc, cb = one_combo([{"type": "memberDistributed", "member": "ab", "direction": "Z", "w1": -0.003},
                        nodal("b", Fz=-5.0)])
    add("ex-5-8-grid", "Grid (plane frame loaded normal to its plane): bending and torsion",
        "Example 5.8", "pp. 116-119",
        "Grid in the horizontal x-y plane, z vertical (model: verticalAxis Z; the source's axes are the model's). "
        "a (0, 0, 0), b (5, 0, 0), c (5, 8, 0) m, right angle at b; a and c fixed. Member ab (5 m) has the properties of "
        "member bc of Example 4.8 (I = 50e6 mm^4, J = 100e3 mm^4, A = 4e3 mm^2), member bc (8 m) those of member ab of "
        "Example 4.8 (I = 200e6 mm^4, J = 300e3 mm^4, A = 6e3 mm^2); E = 200,000 MPa, nu = 0.3. Default orientation: "
        "local y vertical for both members, so vertical loads bend them about local z and the source's I is Iz (Iy is "
        "a placeholder). Loads: 3 kN/m downward on ab, 5 kN downward at b. First order, no shear deformation.",
        model("ex-5-8-grid", ("kN", "mm"), "Z", {"a": (0, 0, 0), "b": (5000, 0, 0), "c": (5000, 8000, 0)},
              {"a": FIXED, "c": FIXED}, SI, {"ab": sec(4e3, 50e6, 100e3), "bc": sec(6e3, 200e6, 300e3)},
              {"ab": frame("a", "b", "ab"), "bc": frame("b", "c", "bc")}, lc, cb),
        [P("22.37", disp("b", "uz"), "mm", sign=-1, note="printed -22.37"),
         P("4.195e-3", disp("b", "rx"), "rad", note="printed 4.195 x 10^-3"),
         P("5.931e-3", disp("b", "ry"), "rad", note="printed 5.931 x 10^-3"),
         P("14.74", reac("a", "Fz"), "kN"),
         P("6.45", reac("a", "Mx"), "kN mm", sign=-1, note="printed -6.45 kNmm (torque in ab)"),
         P("36.21e3", reac("a", "My"), "kN mm", sign=-1, note="printed -36.21 x 10^3 kNmm"),
         P("5.25", reac("c", "Fz"), "kN"),
         P("41.94e3", reac("c", "Mx"), "kN mm", sign=-1, note="printed -41.94 x 10^3 kNmm"),
         P("17.11", reac("c", "My"), "kN mm", sign=-1, note="printed -17.11 kNmm (torque in bc)")])

    # ------------------------------------------------------------------------------------- Ex. 5.9
    AR = B.A_RIGID
    lcs = {"W": {"loads": [{"type": "memberDistributed", "member": "bc", "direction": "Y", "w1": -0.002}]},
           "H": {"loads": [nodal("b", Fx=2.5), nodal("c", Fx=2.5)]}}
    cbs = {"W": {"factors": {"W": 1.0}}, "H": {"factors": {"H": 1.0}}}
    add("ex-5-9-portal-axially-rigid", "Portal frame, axial deformation neglected: vertical and horizontal load",
        "Example 5.9", "pp. 119-120",
        PLANE + "a (0, 0), b (0, 5), c (8, 5), d (8, 0) m; a and d fixed. Columns I = 50e6 mm^4, beam I = 200e6 mm^4, "
        "E = 200,000 MPa. The source neglects axial deformation. " + B.RIGID_NOTE + " Combination W: 2 kN/m downward on "
        "the beam; combination H: 2.5 kN (+x) at b and at c. First order, no shear deformation. " + UNUSED,
        model("ex-5-9-portal-axially-rigid", ("kN", "mm"), "Y",
              {"a": (0, 0), "b": (0, 5000), "c": (8000, 5000), "d": (8000, 0)}, {"a": FIXED, "d": FIXED}, SI,
              {"col": sec(AR, 50e6), "beam": sec(AR, 200e6)},
              {"ab": frame("a", "b", "col"), "bc": frame("b", "c", "beam"), "dc": frame("d", "c", "col")}, lcs, cbs),
        [P("0.000593", disp("b", "rz", "W"), "rad", sign=-1, note="printed -0.000593", signnote=MOM),
         P("0.000593", disp("c", "rz", "W"), "rad",
           note="the source states theta_c = -theta_b (symmetry)", signnote=MOM),
         P("3.09", disp("b", "ux", "H"), "mm"),
         P("3.09", disp("c", "ux", "H"), "mm", note="the source states u_c = u_b"),
         P("0.000195", disp("b", "rz", "H"), "rad", sign=-1, note="printed -0.000195", signnote=MOM),
         P("0.000195", disp("c", "rz", "H"), "rad", sign=-1,
           note="printed -0.000195; the source states theta_c = theta_b (antisymmetry)", signnote=MOM)])
