#!/usr/bin/env python3
"""Writes the problem directories of verification/problems/ziemian/ from the published data set

    C. W. Ziemian, R. D. Ziemian, "Steel benchmark frames for structural analysis and validation studies:
    Finite element models and numerical simulation data", Data in Brief 39 (2021) 107564,
    doi:10.1016/j.dib.2021.107564; data: Mendeley Data, doi:10.17632/39sjhchwtx.1, licence CC BY 4.0.

The two spreadsheets of the data set are NOT in this repository.  Download them and run

    python3 verification/problems/ziemian/build_from_dataset.py <directory with model.xlsx and results.xlsx>
    python3 verification/problems/ziemian/build_from_dataset.py <directory> --infer-members   (see MEMBERS)

(model.xlsx = "FEA Model Details_All frames.xlsx", results.xlsx = "Simulation Results_All frames.xlsx").

What is taken from the data set, per frame: node coordinates (perfect geometry), fixities, element
connectivity, end connections, section properties A, Izz and J (the listed Iyy is replaced by 1000 Izz:
the data set's analyses are plane, the engine's member would otherwise buckle out of plane between the
joints), nodal loads; from the results file: the ELASTIC results only -- first-order (LA) and
second-order (WC) lateral joint displacements and member design moments, and the elastic critical load
ratio (LBA).  The approximate SIPC results and everything inelastic are not used.

Analysis conditions of the data set (article, Tables 2 and 3) and how they are reproduced here:

* stiffness 0.8 E for every elastic analysis: the model file's element table points to its material 1
  (E = 29,000 ksi), but the published results are obtained with 0.8 E = 23,200 ksi.  Measured for EVERY
  frame used here: the first-order displacements are reproduced to 1e-5 or better with 23,200 ksi and are
  20 % smaller with 29,000 ksi (README.md, "conditions verified per frame");
* no shear deformation (the file's shear areas are "Inf");
* LA and WC on the IMPERFECT geometry: every node shifted horizontally by (its height) / 500; LBA on the
  PERFECT geometry;
* the two published cases.  "Right": nodes at x + y/500 and the lateral loads acting toward +x; "Left":
  nodes at x - y/500 and the lateral loads acting toward -x (the mirrored direction of imperfection AND
  lateral load; gravity loads unchanged).  The model file lists the imperfect coordinates with the
  rightward lean for every frame, and the lateral loads in the direction of the description file -- for
  frames 5 and 7 that direction is leftward, so the file's own combination (rightward lean, leftward
  loads) is neither published case.  The definition used here was confirmed with the first-order results
  of every frame used (README.md);
* the model file lists coordinates and loads with five significant figures (e.g. the quarter points of
  the straight rafters of frame 2 at y = 339.88, 379.75, 419.62); interior nodes that lie on the straight
  line between the member ends within that rounding are taken as on the line;
* each member of the data set is 4 line elements with the distributed loads lumped at the quarter
  points.  Here a straight line of elements of one section with rigid internal joints is ONE engine member
  and the loads at its interior nodes are point loads on the member at the same places: the same
  structure, since the engine's member is exact for point loads (first order) and carries the exact
  stability functions (second order, buckling).  Result stations at the quarter points in the problems
  that have member design moments as targets.

Printed precision of a target: the spreadsheet stores doubles and DISPLAYS them with a fixed number of
decimals (cell number format): displacements "0.000", moments "0.00", critical load ratios "0.000".  The
target is the number as displayed, the tolerance half a unit of its last displayed digit; the stored
value is recorded in the precisionNote.

Member design moments: the data set numbers its members without defining the numbering.  ``MEMBERS``
below is INFERRED (``--infer-members``): for every published member, the straight member line (a) whose
plastic moment Zzz x Fy equals the tab's "plastic section modulus x Fy" and (b) whose largest first-order
moment magnitude is closest to the published LA values of both cases; accepted only if the match is
unique and closer than 1e-4.  Frames whose published member count differs from the number of straight
member lines (3, 5, 7) get no moment targets.  The first-order moment targets are therefore not evidence
independent of the identification; the second-order ones use an identification fixed beforehand.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import dataset  # noqa: E402
from v2b_common import Library, P, buck, disp, extreme, model  # noqa: E402
from quarantine import QUARANTINE  # noqa: E402

FRAMES = {1: "single-storey portal frame, fixed bases",
          2: "single-storey gabled frame, pinned bases",
          3: "two-storey one-bay frame with pin-ended beams and diagonals",
          5: "single-storey two-bay frame with unequal column heights",
          7: "single-storey two-bay frame with a diagonal brace and a lean-on column",
          9: "two-storey two-bay frame, gravity loads only"}
E_NOMINAL, E_USED = 29000.0, 23200.0

#: published member number -> (first element, last element) of the straight member line; inferred, see docstring
MEMBERS = {1: {1: (1, 4), 2: (9, 12), 3: (5, 8)},
           2: {1: (1, 4), 2: (5, 8), 3: (9, 12), 4: (13, 16)},
           9: {1: (1, 4), 2: (5, 8), 3: (9, 12), 4: (21, 24), 5: (25, 28), 6: (29, 32), 7: (13, 16), 8: (17, 20),
               9: (33, 36), 10: (37, 40)}}

ATTRIBUTION = ("C. W. Ziemian and R. D. Ziemian, 'Steel benchmark frames for structural analysis and validation "
               "studies: Finite element models and numerical simulation data', Data in Brief 39 (2021) 107564, "
               "doi:10.1016/j.dib.2021.107564; data set: Mendeley Data, doi:10.17632/39sjhchwtx.1; licence CC BY 4.0 "
               "(https://creativecommons.org/licenses/by/4.0/). The model data of this frame are reproduced from the "
               "data set under that licence, converted to this project's input format (changes: 4 elements per member "
               "merged into one member, plane model embedded in 3D, 0.8 E entered as the modulus).")


# ------------------------------------------------------------------------------------- geometry
def rounding(v):
    """Half a unit in the fifth significant figure of a listed coordinate (the listing's precision)."""
    return 0.0 if v == 0 else 0.5 * 10.0 ** (math.floor(math.log10(abs(v))) - 4)


def lines(m):
    """Straight member lines: lists of element numbers joined end to end through nodes that connect exactly
    two elements, carry no support, have rigid connections and lie on the straight line through the outer
    nodes of the two elements (perfect geometry) within the rounding of the listed coordinates: the model
    file lists coordinates with five significant figures (frame 2: rafter nodes at y = 339.88, 379.75,
    419.62 for the quarter points 339.875, 379.75, 419.625 of a straight rafter)."""
    els, X = m["elements"], m["perfect"]
    at = {}
    for k, e in els.items():
        at.setdefault(e["i"], []).append(k)
        at.setdefault(e["j"], []).append(k)
    supported = {k for k, f in m["fixities"].items() if "FIXED" in f}

    def joins(a, b):
        ea, eb = els[a], els[b]
        n = ea["j"]
        if n != eb["i"] or len(at[n]) != 2 or n in supported or ea["section"] != eb["section"]:
            return False
        if m["connections"][a][2] != "Rigid" or m["connections"][b][0] != "Rigid":
            return False
        p, q, r = X[ea["i"]], X[n], X[eb["j"]]
        cross = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
        return abs(cross) / math.dist(p, r) <= max(rounding(q[0]), rounding(q[1]))

    nxt = {a: b for a in els for b in els if a != b and joins(a, b)}
    starts = sorted(set(els) - set(nxt.values()))
    out = []
    for s in starts:
        line = [s]
        while line[-1] in nxt:
            line.append(nxt[line[-1]])
        out.append(line)
    assert sorted(k for ln in out for k in ln) == sorted(els)
    return out


def member_id(line):
    return f"m{line[0]}" if len(line) == 1 else f"m{line[0]}-{line[-1]}"


def lateral_sign(m):
    total = sum(v[0] for v in m["nodalForces"].values())
    return 1.0 if total >= 0.0 else -1.0


def document(m, case, kind, buckling=False, stations=False):
    """benchfea input document.  case: "right" | "left" (imperfect geometry) | "perfect".  ``stations``: result
    stations at the quarter points of every member (needed only where member design moments are targets)."""
    lean = {"right": 1.0, "left": -1.0, "perfect": 0.0}[case]
    # lateral loads: toward +x in the "right" case, toward -x in the "left" case; as listed for the LBA
    lsign = lateral_sign(m) * lean if lean else 1.0
    X = {k: (x + lean * y / 500.0, y, 0.0) for k, (x, y, z) in m["perfect"].items()}
    lns = lines(m)
    interior = {}
    for ln in lns:
        a, b = m["elements"][ln[0]]["i"], m["elements"][ln[-1]]["j"]
        length = math.dist(m["perfect"][a], m["perfect"][b])
        for k in ln[:-1]:
            n = m["elements"][k]["j"]
            interior[n] = (member_id(ln), math.dist(m["perfect"][a], m["perfect"][n]) / length)
    nodes = {f"n{k}": X[k] for k in sorted(X) if k not in interior}
    supports = {}
    for k in sorted(X):
        if k in interior:
            continue
        f = m["fixities"][k]
        dofs = {"uz": "fixed", "rx": "fixed", "ry": "fixed"}
        for name, flag in (("ux", f[0]), ("uy", f[1]), ("rz", f[5])):
            if flag == "FIXED":
                dofs[name] = "fixed"
        supports[f"n{k}"] = {"dofs": dofs}
    sections = {f"s{k}": {"type": "properties",
                          "description": str(s["name"]).strip() + ": A, Izz and J as listed in the data set; Iy = "
                                         "1000 Izz (plane problem, see problem.json)",
                          "A": float(s["A"]), "Iz": float(s["Izz"]), "Iy": 1000.0 * float(s["Izz"]), "J": float(s["J"])}
                for k, s in sorted(m["sections"].items())}
    members = {}
    for ln in lns:
        first, last = m["elements"][ln[0]], m["elements"][ln[-1]]
        mm = {"i": f"n{first['i']}", "j": f"n{last['j']}", "section": f"s{first['section']}", "material": "M"}
        rel = {}
        if m["connections"][ln[0]][0] == "Pinned":
            rel["i"] = {"rz": True}
        if m["connections"][ln[-1]][2] == "Pinned":
            rel["j"] = {"rz": True}
        if rel:
            mm["releases"] = rel
        members[member_id(ln)] = mm
    loads = []
    for k, (fx, fy, fz) in sorted(m["nodalForces"].items()):
        fx = fx * lsign
        if k in interior:
            mid, pos = interior[k]
            if fy:
                loads.append({"type": "memberPoint", "member": mid, "direction": "Y", "value": fy, "position": pos})
            if fx:
                loads.append({"type": "memberPoint", "member": mid, "direction": "X", "value": fx, "position": pos})
        else:
            ld = {"type": "nodal", "node": f"n{k}"}
            if fx:
                ld["Fx"] = fx
            if fy:
                ld["Fy"] = fy
            loads.append(ld)
    assert not m["nodalMoments"] and not m["uniformLoads"]
    an = {"shearDeformation": False}
    if stations:
        an["resultStations"] = 5
    if buckling:
        an["buckling"] = {"combinations": ["C"], "modes": 1}
    return model("", ("kip", "in"), "Y", nodes, supports,
                 {"M": {"E": E_USED, "nu": 0.3, "description": "0.8 x 29,000 ksi (analysis condition of the data set)"}},
                 sections, members, {"L": {"loads": loads}},
                 {"C": {"factors": {"L": 1.0}, "analysisType": kind}}, analysis=an)


# ------------------------------------------------------------------------------------- targets
def shown(v, decimals):
    return f"{v:.{decimals}f}"


def displacement_targets(res, key, side):
    out = []
    for j, v in sorted(res["joints"].items()):
        x = v[key][side]
        if x == 0:
            continue                       # supports, and entries the data set filtered out (stored as 0)
        out.append(P(shown(x, 3), disp(f"n{j}", "ux"), "in",
                     note=f"number as displayed by the data set's spreadsheet (cell format 0.000); stored value {x!r}",
                     signnote="same as the data set: displacement along +x (to the right) positive"))
    return out


def moment_targets(frame, res, key, side):
    out = []
    for k, v in sorted(res["members"].items()):
        if k not in MEMBERS.get(frame, {}):
            continue
        x = v[key][side]
        a, b = MEMBERS[frame][k]
        out.append(P(shown(x, 2), extreme(f"m{a}-{b}", "Mz"), "kip in",
                     note=f"number as displayed by the data set's spreadsheet (cell format 0.00); stored value {x!r}; "
                          f"member {k} of the data set identified as the member line of elements {a} to {b} "
                          "(identification inferred from the first-order values, see README.md)",
                     signnote="magnitude: 'design moment' = largest moment magnitude in the member (here: over the "
                              "result stations at the quarter points, the element ends of the data set)"))
    return out


def source(frame, where):
    return {"work": "C. W. Ziemian, R. D. Ziemian, Steel benchmark frames for structural analysis and validation "
                    "studies: Finite element models and numerical simulation data, Data in Brief 39 (2021) 107564, "
                    "doi:10.1016/j.dib.2021.107564, with the data set on Mendeley Data, doi:10.17632/39sjhchwtx.1",
            "edition": "data set version 1 (2021), licence CC BY 4.0",
            "location": f"Frame {frame}: 'FEA Model Details_All frames.xlsx', tab 'Frame {frame}'; "
                        f"'Simulation Results_All frames.xlsx', {where}",
            "pages": "article Tables 2 and 3 (modelling and analysis conditions); 'Description Of Frames.pdf'",
            "notes": ATTRIBUTION}


def conditions(m, frame):
    return (f"Frame {frame} of the data set ({FRAMES[frame]}), units kip and inch. Plane frame in the global X-Y plane "
            "(Y vertical); the out-of-plane DOFs (uz, rx, ry) of every node are held (the data set's frames are "
            "'fully braced out-of-plane'); in-plane bending about local z, so the data set's Izz is the section's Iz. "
            f"Modulus 0.8 E = {E_USED:g} ksi (the data set's condition for all elastic analyses; its element table "
            f"points to E = {E_NOMINAL:g} ksi, the published results correspond to 0.8 E - verified for this frame "
            "with the first-order displacements). No shear deformation. A straight line of the data set's elements "
            "is one member here, the loads at its interior nodes are point loads on the member. The data set's model "
            "is plane (no out-of-plane buckling of a member between joints); the engine's member is spatial, so Iy = "
            "1000 Izz is entered for every section instead of the listed Iyy: no out-of-plane mode can precede the "
            "in-plane behaviour compared, and Iy does not enter the in-plane results. ")


CASE_TEXT = {"right": "Case 'Right' of the data set: imperfect geometry with every node at x + y/500, lateral loads "
                      "acting toward +x. ",
             "left": "Case 'Left' of the data set: imperfect geometry with every node at x - y/500, lateral loads "
                     "acting toward -x (gravity loads unchanged). "}
LA_TEXT = "First-order (linear) analysis: the data set's 'LA'. "
WC_TEXT = ("Second-order elastic analysis, compared with the data set's 'WC' results: MASTAN2, work-control "
           "incremental analysis with step size 0.001, updated Lagrangian formulation with geometric stiffness "
           "matrices, 4 elements per member (article, Tables 2 and 3). The engine solves the small-displacement "
           "second-order theory with exact member stability functions (docs/THEORY.md 10, 12). ")
LBA_TEXT = ("Linear buckling analysis on the PERFECT geometry with the loads as listed in the model file; target: the "
            "data set's elastic critical load ratio (MASTAN2 eigenvalue analysis, 4 elements per member with a "
            "polynomial geometric stiffness; 'lowest sway buckling mode'). The engine's members carry the exact "
            "stability functions. ")


def build(lib, folder):
    mx, rx = str(Path(folder) / "model.xlsx"), str(Path(folder) / "results.xlsx")
    acr = dataset.read_alpha_cr(rx)
    for frame in FRAMES:
        m, res = dataset.read_model(mx, frame), dataset.read_results(rx, frame)
        assert m["units"] == "kip, inch" and all(e["beta"] == 0 for e in m["elements"].values())
        cond = conditions(m, frame)
        has_m = bool(MEMBERS.get(frame))
        for side, case in enumerate(("right", "left")):
            tab = f"tab 'Frame {frame}', columns LA / {case.capitalize()}"
            pid = f"frame-{frame}-first-order-displacements-{case}"
            lib.add(pid, f"Ziemian frame {frame} ({FRAMES[frame]}): first-order joint displacements, case '{case}'",
                    source(frame, tab), cond + CASE_TEXT[case] + LA_TEXT, document(m, case, "firstOrder"),
                    displacement_targets(res, "LA", side), diagnosis=QUARANTINE.get(pid))
            if has_m:
                pid = f"frame-{frame}-first-order-moments-{case}"
                lib.add(pid, f"Ziemian frame {frame} ({FRAMES[frame]}): first-order member design moments, case "
                        f"'{case}'", source(frame, tab), cond + CASE_TEXT[case] + LA_TEXT,
                        document(m, case, "firstOrder", stations=True), moment_targets(frame, res, "LA", side),
                        diagnosis=QUARANTINE.get(pid),
                        limits="The member numbering of the data set is not defined in its files; the identification "
                               "of the members was inferred from the first-order values (README.md), so these "
                               "targets are not independent evidence of the identification.")
            pid = f"frame-{frame}-second-order-{case}"
            targets = displacement_targets(res, "WC", side) + moment_targets(frame, res, "WC", side)
            lib.add(pid, f"Ziemian frame {frame} ({FRAMES[frame]}): second-order joint displacements"
                    + (" and member design moments" if has_m else "") + f", case '{case}'",
                    source(frame, f"tab 'Frame {frame}', columns WC / {case.capitalize()}"),
                    cond + CASE_TEXT[case] + WC_TEXT, document(m, case, "secondOrder", stations=has_m), targets,
                    diagnosis=QUARANTINE.get(pid))
        pid = f"frame-{frame}-critical-load-ratio"
        lib.add(pid, f"Ziemian frame {frame} ({FRAMES[frame]}): elastic critical load ratio",
                source(frame, "tab 'SUMMARY', column acr"), cond + LBA_TEXT,
                document(m, "perfect", "firstOrder", buckling=True),
                [P(shown(acr[frame], 3), buck(1, "C"), "-",
                   note="number as displayed by the data set's spreadsheet (cell format 0.000); stored value "
                        f"{acr[frame]!r}", signnote="positive factor on the listed loads")],
                diagnosis=QUARANTINE.get(pid))
    unknown = set(QUARANTINE) - set(lib.problems)
    if unknown:
        raise SystemExit("quarantine.py names unknown problems: " + ", ".join(sorted(unknown)))


def infer_members(folder):
    """Prints, per frame, the identification of the data set's member numbers (module docstring)."""
    from benchfea.model import build_model
    from benchfea.solve import analyze
    mx, rx = str(Path(folder) / "model.xlsx"), str(Path(folder) / "results.xlsx")
    for frame in FRAMES:
        m, res = dataset.read_model(mx, frame), dataset.read_results(rx, frame)
        lns = lines(m)
        fy = m["materials"][1]["Fy"]
        got = {}
        for case in ("right", "left"):
            mdl = build_model(document(m, case, "firstOrder", stations=True))
            c = analyze(mdl).to_document(mdl)["combinations"]["C"]
            got[case] = {member_id(ln): max(abs(x) for x in c["members"][member_id(ln)]["stations"]["Mz"])
                         for ln in lns}
        print(f"frame {frame}: {len(res['members'])} published members, {len(lns)} straight member lines")
        if len(res["members"]) != len(lns):
            print("   member count differs: no moment targets")
            continue
        for k, v in sorted(res["members"].items()):
            cands = []
            for ln in lns:
                mp = m["sections"][m["elements"][ln[0]]["section"]]["Zzz"] * fy
                if abs(mp - v["Mp"]) > 1e-6 * v["Mp"]:
                    continue
                mid = member_id(ln)
                err = max(abs(got[c][mid] - v["LA"][s]) / max(abs(v["LA"][s]), 1e-9)
                          for s, c in enumerate(("right", "left")))
                cands.append((err, ln[0], ln[-1]))
            cands.sort()
            ok = cands and cands[0][0] < 1e-4 and (len(cands) == 1 or cands[1][0] > 1e-2)
            print(f"   member {k}: " + (f"elements {cands[0][1]} to {cands[0][2]} (largest LA difference "
                                        f"{cands[0][0]:.1e})" if ok else f"NOT identified: {cands[:3]}"))


def main(argv):
    args = [a for a in argv if not a.startswith("-")]
    if len(args) != 1:
        raise SystemExit(__doc__)
    if "--infer-members" in argv:
        infer_members(args[0])
        return
    lib = Library(HERE, "ziemian")
    build(lib, args[0])
    lib.write()


if __name__ == "__main__":
    main(sys.argv[1:])
