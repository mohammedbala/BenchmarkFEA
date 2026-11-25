#!/usr/bin/env python3
"""Writes the problem directories of verification/problems/aisc-commentary/.

Source: ANSI/AISC 360-22, Specification for Structural Steel Buildings, Commentary to Section C2.1,
benchmark problems Case 1 (Fig. C-C2.2) and Case 2 (Fig. C-C2.3), p. 16.1-337.  Nothing of the source is
reproduced: the input data that define the two problems and the tabulated target numbers are given with
their location, in this project's own words.

Section properties.  The figures name the shape (W14x48 / W360x72, major-axis bending) and E, but give no
I_x, A, shear area or G.  No property is taken from memory or from a source that was not opened; exactly
as in tests/test_second_order_aisc_benchmarks.py:

* I_x is CALIBRATED on the Commentary's own first-order (P = 0) deflection without shear deformation of
  Case 2:  I = H L^3 / (3 E Delta_tip)  (US table: 483.9 in^4; SI table: 4.021e8 mm^4, each from its own
  table).  The P = 0, no-shear deflection of Case 2 is therefore reproduced by construction and tests
  nothing; every other number does.
* the shear flexibility 1/(G A_s) is CALIBRATED per case and per table on the difference between the
  P = 0 deflections with and without shear deformation, which is printed to one significant digit only;
  the rows "with shear deformation" are therefore a consistency test of the amplification, not a test of
  a shear area.  G = E/2.6 is nominal: only the product G A_s enters.
* A does not enter these results in this theory (tests/test_second_order_aisc_benchmarks.py asserts it);
  the model uses A = I / (L/50)^2.

Tolerance: the Commentary's own acceptance criterion on the same page: moments within 3 %, deflections
within 5 %.  The US and the SI tables are independently rounded and are separate problems; so is every
axial load level (one problem per level: 2 cases x 2 tables x with / without shear x 4 levels = 32),
so that a level that cannot be published (quarantine.py) does not take the others with it.

Run:  python3 verification/problems/aisc-commentary/build_problems.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from v2b_common import Library, S, disp, model, station  # noqa: E402
from quarantine import QUARANTINE  # noqa: E402

SOURCE = {"work": "ANSI/AISC 360-22, Specification for Structural Steel Buildings (American Institute of Steel "
                  "Construction), Commentary",
          "edition": "August 1, 2022",
          "pages": "p. 16.1-337 (text p. 16.1-336 to 16.1-337)",
          "notes": "input data and targets read on the cited page; nothing of the source is reproduced. Section "
                   "properties of the W14x48 are not given by the source and are not taken from any table: I_x is "
                   "calibrated on the source's own P = 0 deflection without shear deformation of Case 2, the shear "
                   "flexibility on the source's P = 0 deflections with and without shear deformation (see "
                   "build_problems.py and README.md)"}
CRITERION = ("acceptance criterion stated by the source on the same page (AISC 360-22 Commentary C2.1, p. 16.1-337): "
             "moments within 3 %, deflections within 5 %; value printed to 3 significant figures")

# axial load levels and (value with shear deformation, value without) as printed in the two figures
TABLES = {
    "us": {"units": ("kip", "in"), "E": 29000.0, "L": 28.0 * 12.0, "w": 0.200 / 12.0, "H": 1.0, "kM": 1.0, "kP": 1.0,
           "uM": "kip in", "uD": "in", "text": "US customary table: E = 29,000 ksi, L = 28.0 ft = 336 in, ",
           "loads": {1: "0.200 kip/ft", 2: "1.0 kip"},
           1: {"P": ["0", "150", "300", "450"],
               "M": [("235", "235"), ("270", "269"), ("316", "313"), ("380", "375")],
               "D": [("0.202", "0.197"), ("0.230", "0.224"), ("0.269", "0.261"), ("0.322", "0.311")]},
           2: {"P": ["0", "100", "150", "200"],
               "M": [("336", "336"), ("470", "469"), ("601", "598"), ("856", "848")],
               "D": [("0.907", "0.901"), ("1.34", "1.33"), ("1.77", "1.75"), ("2.60", "2.56")]}},
    "si": {"units": ("N", "mm"), "E": 200000.0, "L": 8530.0, "w": 2.92, "H": 4450.0, "kM": 1.0e6, "kP": 1.0e3,
           "uM": "N mm", "uD": "mm", "text": "SI table: E = 200 000 MPa, L = 8.53 m, ",
           "loads": {1: "2.92 kN/m", 2: "4.45 kN"},
           1: {"P": ["0", "667", "1330", "2000"],
               "M": [("26.6", "26.6"), ("30.5", "30.4"), ("35.7", "35.4"), ("42.9", "42.4")],
               "D": [("5.13", "5.00"), ("5.84", "5.69"), ("6.83", "6.63"), ("8.18", "7.90")]},
           2: {"P": ["0", "445", "667", "890"],
               "M": [("38.0", "38.0"), ("53.1", "53.0"), ("67.9", "67.6"), ("96.7", "95.8")],
               "D": [("23.0", "22.9"), ("34.0", "33.8"), ("45.0", "44.5"), ("66.0", "65.0")]}},
}


def calibrated(t, case):
    """(I, G A_s) from the P = 0 columns of the Commentary (module docstring)."""
    L, E = t["L"], t["E"]
    inertia = t["H"] * L ** 3 / (3.0 * E * float(t[2]["D"][0][1]))
    with_shear, without = (float(x) for x in t[case]["D"][0])
    gas = t["w"] * L ** 2 / (8.0 * (with_shear - without)) if case == 1 else t["H"] * L / (with_shear - without)
    return inertia, gas


def build(lib):
    for key, t in TABLES.items():
        for case in (1, 2):
            for shear in (False, True):
                L, E = t["L"], t["E"]
                inertia, gas = calibrated(t, case)
                G = E / 2.6
                sec = {"type": "properties", "A": inertia / (L / 50.0) ** 2, "Iy": inertia, "Iz": inertia,
                       "J": inertia / 100.0,
                       "description": "W14x48 of the source; I calibrated on the source's own P = 0 deflection, see "
                                      "problem.json (source.notes); A, Iy, J do not enter the targets"}
                if shear:
                    sec["Asy"] = gas / G
                if case == 1:
                    supports = {"a": {"type": "pinned", "dofs": {"rz": "fixed"}},
                                "b": {"dofs": {"ux": "fixed", "uy": "fixed"}}}
                    lcs = {"T": {"loads": [{"type": "memberDistributed", "member": "m", "direction": "X",
                                            "w1": t["w"]}]}}
                else:
                    supports = {"a": {"type": "fixed"}, "b": {"dofs": {"uy": "fixed", "rx": "fixed"}}}
                    lcs = {"T": {"loads": [{"type": "nodal", "node": "b", "Fx": t["H"]}]}}
                lcs["P"] = {"loads": [{"type": "nodal", "node": "b", "Fz": -t["kP"]}]}
                for k, p in enumerate(t[case]["P"]):
                    f = {"T": 1.0}
                    if float(p) != 0.0:
                        f["P"] = float(p)
                    combos = {"C": {"factors": f, "analysisType": "secondOrder"}}
                    m_ref, d_ref = t[case]["M"][k][0 if shear else 1], t[case]["D"][k][0 if shear else 1]
                    msign = ("magnitude in the source; in the engine's convention (member local x upward, local y = "
                             "global +X, loads along +X) the midspan moment Mz of Case 1 is negative and the base "
                             "moment of Case 2 positive")
                    note = ""
                    if case == 2 and not shear and float(p) == 0.0:
                        note = "I_x is calibrated on this number: reproduced by construction, not evidence"
                    if shear and float(p) == 0.0:
                        note = ("the shear flexibility is calibrated on this number minus the value without shear "
                                "deformation: reproduced by construction, not evidence")
                    if case == 1:
                        targets = [S(m_ref, 0.03, CRITERION, station("m", "Mz", rel=0.5), t["uM"], sign=-1.0,
                                     factor=t["kM"], signnote=msign),
                                   S(d_ref, 0.05, CRITERION, station("m", "uy", rel=0.5), t["uD"], note=note,
                                     signnote="magnitude in the source; the load acts along +X, deflection +X "
                                              "(local y)")]
                    else:
                        targets = [S(m_ref, 0.03, CRITERION, station("m", "Mz", rel=0.0), t["uM"], factor=t["kM"],
                                     signnote=msign),
                                   S(d_ref, 0.05, CRITERION, disp("b", "ux"), t["uD"], note=note,
                                     signnote="magnitude in the source; the load acts along +X, deflection +X")]
                    an = {"shearDeformation": shear, "resultStations": 3}
                    mdl = model("", t["units"], "Z", {"a": (0, 0, 0), "b": (0, 0, L)}, supports,
                                {"M": {"E": E, "G": G}}, {"S": sec},
                                {"m": {"i": "a", "j": "b", "section": "S", "material": "M"}}, lcs, combos, analysis=an)
                    what = ("simply supported beam-column (pin at the bottom, lateral roller at the loaded top) under "
                            f"a uniform transverse load of {t['loads'][1]} and the axial compression P; targets: "
                            "moment and deflection at midspan" if case == 1 else
                            f"fixed-base cantilever under a lateral load of {t['loads'][2]} and the axial "
                            "compression P at the top; targets: moment at the base and lateral deflection of the top")
                    analysis = (f"Case {case} of the source: {what}. Axial load level P = {p} x {t['kP']:g} "
                                f"{t['units'][0]} (load case P holds the unit axial load, the combination factor is "
                                "the level). " + t["text"] +
                                f"I = {inertia:.6g} {t['units'][1]}^4 (calibrated, see source.notes), "
                                + (f"G A_s = {gas:.6g} {t['units'][0]} (calibrated on the P = 0 column of this case; "
                                   "shear deformation included: the source's unbracketed values). " if shear else
                                   "shear deformation excluded (the source's values in brackets). ") +
                                "Member vertical along global Z, bending in the X-Z plane about local z (the "
                                "section's Iz = Iy = I). Second-order analysis; the engine's member carries the exact "
                                "stability functions, so one member is the exact solution of the second-order theory. "
                                "The support idealisation is read from the source's sketch.")
                    pid = f"case-{case}-{key}-{'with' if shear else 'without'}-shear-p{p}"
                    lib.add(pid,
                            f"AISC benchmark Case {case} ({'P-delta only' if case == 1 else 'P-Delta and P-delta'}), "
                            f"{key.upper()} table, {'with' if shear else 'without'} shear deformation, P = {p}",
                            dict(SOURCE, location=f"Commentary Section C2.1, Fig. C-C2.{case + 1} (benchmark problem "
                                                  f"Case {case}), {'US customary' if key == 'us' else 'SI'} table, "
                                                  f"column P = {p}, "
                                                  + ("unbracketed values" if shear else "values in brackets")),
                            analysis, mdl, targets, diagnosis=QUARANTINE.get(pid),
                            limits="Section properties are calibrated on the source's own P = 0 values (the source "
                                   "gives none): the problem verifies the second-order amplification, not section "
                                   "properties."
                                   + (" It does NOT test shear deformation: the source's acceptance tolerance (3 % on "
                                      "moments, 5 % on deflections) is wider than the whole effect of shear "
                                      "deformation in its tables (at most 1.3 % on a moment and 3.5 % on a "
                                      "deflection between the bracketed and unbracketed values), so these targets "
                                      "would also be met without shear deformation; and the shear stiffness is "
                                      "calibrated on a difference the source prints to one significant digit."
                                      if shear else ""))
    unknown = set(QUARANTINE) - set(lib.problems)
    if unknown:
        raise SystemExit("quarantine.py names unknown problems: " + ", ".join(sorted(unknown)))


def main():
    lib = Library(HERE, "aisc-commentary")
    build(lib)
    lib.write()


if __name__ == "__main__":
    main()
