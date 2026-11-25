#!/usr/bin/env python3
"""Evidence for the quarantine diagnoses of verification/problems/ziemian/ (not a verification tier).

    python3 verification/derive/ziemian_diagnosis.py <directory with model.xlsx and results.xlsx> [frame ...]

Per frame it prints

1. LBA.  The critical load ratio of the data set's OWN discretisation, recomputed here with a small plane
   frame program written for this purpose: every listed element is one cubic element with the standard
   polynomial (consistent) geometric stiffness, axial forces from a linear analysis, perfect geometry,
   0.8 E.  It is compared with the published ratio and with the engine (exact stability functions).
2. Second order, case "Right".  The joint displacements by OpenSeesPy with (a) the linearised "PDelta"
   transformation and (b) the "Corotational" (large displacement) transformation, each listed element
   subdivided into 8, load applied in 20 steps with Newton iteration; compared with the engine and with
   the data set's WC values.  (a) is the theory the engine solves, (b) the theory of the data set.
3. The largest compressive force in a pin-ended single-element member relative to its Euler load
   (frame 7: the diagonal brace), first-order analysis of the case "Right".
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "verification" / "problems" / "ziemian"))
sys.path.insert(0, str(ROOT / "verification" / "problems"))
import dataset  # noqa: E402
import build_from_dataset as bfd  # noqa: E402

E = bfd.E_USED


# ------------------------------------------------------------------ own plane frame (cubic elements)
def plane_frame(m, lean, lsign):
    """Returns (K, Kg(N from linear analysis), u, axial forces) on the free DOFs."""
    X = {k: (x + lean * y / 500.0, y) for k, (x, y, z) in m["perfect"].items()}
    dof = {}
    for k in sorted(X):
        dof[(k, "x")], dof[(k, "y")], dof[(k, "r")] = len(dof), len(dof) + 1, len(dof) + 2
    edofs = {}
    for k, e in m["elements"].items():
        c = m["connections"][k]
        ri, rj = dof[(e["i"], "r")], dof[(e["j"], "r")]
        if c[0] == "Pinned":
            ri = dof[(k, "ri")] = len(dof)
        if c[2] == "Pinned":
            rj = dof[(k, "rj")] = len(dof)
        edofs[k] = [dof[(e["i"], "x")], dof[(e["i"], "y")], ri, dof[(e["j"], "x")], dof[(e["j"], "y")], rj]
    n = len(dof)
    fixed = set()
    for k, f in m["fixities"].items():
        for name, flag in (("x", f[0]), ("y", f[1]), ("r", f[5])):
            if flag == "FIXED":
                fixed.add(dof[(k, name)])
    used = {d for v in edofs.values() for d in v}
    free = [d for d in range(n) if d not in fixed and d in used]
    K = np.zeros((n, n))
    parts = {}
    for k, e in m["elements"].items():
        s = m["sections"][e["section"]]
        A, I = float(s["A"]), float(s["Izz"])
        (xi, yi), (xj, yj) = X[e["i"]], X[e["j"]]
        L = math.hypot(xj - xi, yj - yi)
        c, sn = (xj - xi) / L, (yj - yi) / L
        ke = np.array([[E * A / L, 0, 0, -E * A / L, 0, 0],
                       [0, 12 * E * I / L ** 3, 6 * E * I / L ** 2, 0, -12 * E * I / L ** 3, 6 * E * I / L ** 2],
                       [0, 6 * E * I / L ** 2, 4 * E * I / L, 0, -6 * E * I / L ** 2, 2 * E * I / L],
                       [-E * A / L, 0, 0, E * A / L, 0, 0],
                       [0, -12 * E * I / L ** 3, -6 * E * I / L ** 2, 0, 12 * E * I / L ** 3, -6 * E * I / L ** 2],
                       [0, 6 * E * I / L ** 2, 2 * E * I / L, 0, -6 * E * I / L ** 2, 4 * E * I / L]])
        kg = np.array([[1 / L, 0, 0, -1 / L, 0, 0],
                       [0, 6 / (5 * L), 1 / 10, 0, -6 / (5 * L), 1 / 10],
                       [0, 1 / 10, 2 * L / 15, 0, -1 / 10, -L / 30],
                       [-1 / L, 0, 0, 1 / L, 0, 0],
                       [0, -6 / (5 * L), -1 / 10, 0, 6 / (5 * L), -1 / 10],
                       [0, 1 / 10, -L / 30, 0, -1 / 10, 2 * L / 15]])
        R = np.array([[c, sn, 0], [-sn, c, 0], [0, 0, 1]])
        T = np.zeros((6, 6))
        T[:3, :3] = T[3:, 3:] = R
        parts[k] = (ke, kg, T, edofs[k], E * A / L, E * I, L)
        K[np.ix_(edofs[k], edofs[k])] += T.T @ ke @ T
    F = np.zeros(n)
    for k, (fx, fy, fz) in m["nodalForces"].items():
        F[dof[(k, "x")]] += lsign * fx
        F[dof[(k, "y")]] += fy
    u = np.zeros(n)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    Kg = np.zeros((n, n))
    axial = {}
    for k, (ke, kg, T, d, ea_l, ei, L) in parts.items():
        ul = T @ u[d]
        N = ea_l * (ul[3] - ul[0])                         # tension positive
        axial[k] = (N, math.pi ** 2 * ei / L ** 2)
        Kg[np.ix_(d, d)] += T.T @ (N * kg) @ T
    return K[np.ix_(free, free)], Kg[np.ix_(free, free)], axial


def lba_cubic(m):
    K, Kg, _ = plane_frame(m, 0.0, 1.0)
    w = eigh(-Kg, K, eigvals_only=True)                    # -Kg x = (1/lambda) K x
    return 1.0 / w[-1]


# ------------------------------------------------------------------ OpenSeesPy
def opensees(m, transf, lean, lsign, nsub=8, steps=20):
    import openseespy.opensees as ops
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    X = {k: (x + lean * y / 500.0, y) for k, (x, y, z) in m["perfect"].items()}
    for k, (x, y) in X.items():
        ops.node(k, x, y)
    for k, f in m["fixities"].items():
        if "FIXED" in f:
            ops.fix(k, int(f[0] == "FIXED"), int(f[1] == "FIXED"), int(f[5] == "FIXED"))
    ops.geomTransf(transf, 1)
    nid, eid = 100000, 0
    for k, e in m["elements"].items():
        s, c = m["sections"][e["section"]], m["connections"][k]
        a, b = e["i"], e["j"]
        if c[0] == "Pinned":
            nid += 1
            ops.node(nid, *X[a])
            ops.equalDOF(a, nid, 1, 2)
            a = nid
        if c[2] == "Pinned":
            nid += 1
            ops.node(nid, *X[b])
            ops.equalDOF(b, nid, 1, 2)
            b = nid
        prev = a
        for q in range(1, nsub + 1):
            if q < nsub:
                nid += 1
                t = q / nsub
                ops.node(nid, X[e["i"]][0] + t * (X[e["j"]][0] - X[e["i"]][0]),
                         X[e["i"]][1] + t * (X[e["j"]][1] - X[e["i"]][1]))
                nxt = nid
            else:
                nxt = b
            eid += 1
            ops.element("elasticBeamColumn", eid, prev, nxt, float(s["A"]), E, float(s["Izz"]), 1)
            prev = nxt
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for k, (fx, fy, fz) in m["nodalForces"].items():
        ops.load(k, lsign * fx, fy, 0.0)
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    if ops.analyze(steps) != 0:
        return None
    return {k: ops.nodeDisp(k, 1) for k in X}


def engine(m, case, kind):
    from benchfea.model import build_model
    from benchfea.problems import ModelInvalid
    from benchfea.solve import analyze
    mdl = build_model(bfd.document(m, case, kind, buckling=(case == "perfect")))
    try:
        return analyze(mdl).to_document(mdl)
    except ModelInvalid as e:
        return str(e.problems[0].message)[:160]


def main(argv):
    folder = Path(argv[0])
    frames = [int(a) for a in argv[1:]] or list(bfd.FRAMES)
    acr = dataset.read_alpha_cr(str(folder / "results.xlsx"))
    for f in frames:
        m = dataset.read_model(str(folder / "model.xlsx"), f)
        res = dataset.read_results(str(folder / "results.xlsx"), f)
        print(f"===== frame {f}")
        own = lba_cubic(m)
        eng = engine(m, "perfect", "firstOrder")["buckling"]["C"]["factors"][0]
        print(f"LBA: data set {acr[f]:.6f}; own program, the data set's elements with cubic geometric stiffness "
              f"{own:.6f} (rel. diff. to the data set {own / acr[f] - 1:+.1e}); engine {eng:.6f} "
              f"({eng / acr[f] - 1:+.1e})")
        lsign = bfd.lateral_sign(m)
        doc = engine(m, "right", "secondOrder")
        pd = opensees(m, "PDelta", 1.0, lsign)
        co = opensees(m, "Corotational", 1.0, lsign)
        co1 = opensees(m, "Corotational", 1.0, lsign, nsub=1)
        print("second order, case Right: joint, data set WC, engine, OpenSees PDelta (8 per element), OpenSees "
              "Corotational (8 per element), OpenSees Corotational (1 per element)")
        worst = {"engine-pdelta": 0.0, "engine-wc": 0.0, "corot-wc": 0.0, "corot1-wc": 0.0}
        for j, v in sorted(res["joints"].items()):
            wc = v["WC"][0]
            if wc == 0:
                continue
            en = doc["combinations"]["C"]["displacements"][f"n{j}"]["ux"] if isinstance(doc, dict) else float("nan")
            row = [wc, en] + [(x[j] if x else float("nan")) for x in (pd, co, co1)]
            print(f"   {j:3d} " + " ".join(f"{x:11.6f}" for x in row))
            worst["engine-pdelta"] = max(worst["engine-pdelta"], abs(en / row[2] - 1))
            worst["engine-wc"] = max(worst["engine-wc"], abs(en / wc - 1))
            worst["corot-wc"] = max(worst["corot-wc"], abs(row[3] / wc - 1))
            worst["corot1-wc"] = max(worst["corot1-wc"], abs(row[4] / wc - 1))
        if not isinstance(doc, dict):
            print("   engine: " + doc)
        print("   largest relative differences: " + ", ".join(f"{k} {v:.1e}" for k, v in worst.items()))
        _, _, axial = plane_frame(m, 1.0, lsign)
        for k, (N, pe) in axial.items():
            c = m["connections"][k]
            if c[0] == "Pinned" and c[2] == "Pinned":
                print(f"   pin-ended element {k}: first-order axial force {N:+.2f} kip (tension +), Euler load "
                      f"{pe:.2f} kip, ratio {-N / pe:.3f}")


if __name__ == "__main__":
    main(sys.argv[1:])
