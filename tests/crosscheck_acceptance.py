"""Acceptance of the foreign-error bound of ``benchfea.crosscheck`` (docs/CROSSCHECK.md 1.3).

For every model below an answer exists that depends neither on the engine nor on a foreign solver:

* ``oracle``  -- a straight member line: the ODE oracle (``benchfea.oracle`` through the line
  representation of verification tier A), accurate to about 1e-8;
* ``tableau`` -- a frame: the tableau solver of tests/independent (transfer matrices of the member
  equations, no stiffness matrix).  It is exact for a member with constant axial force; with a varying
  axial force it takes the force at the middle of each of ``nsub`` pieces, so it is run at ``nsub`` and
  ``2 nsub`` and the difference of the two runs is carried as its own uncertainty.

``true error`` of a foreign solver = largest relative difference to that answer in the norm of
``crosscheck.compare`` (nodal displacements and rotations, reactions, member end forces);
``estimate`` = what the solver's fidelity record says.  The criterion is ``estimate >= true error``
(minus the uncertainty of the reference, plus 1e-8 for the oracle).  The engine's answer is not used.

    python3 tests/crosscheck_acceptance.py            # prints the whole table (a few minutes)
    python3 tests/crosscheck_acceptance.py --quick    # the subset the test suite runs
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "independent"))

from benchfea.crosscheck import relative_difference, run                     # noqa: E402
from benchfea.crosscheck.fidelity import Fidelity                            # noqa: E402
from benchfea.crosscheck.foreign import ForeignResult                        # noqa: E402
from benchfea.model import build_model                                       # noqa: E402

ORACLE_ACCURACY = 1.0e-8
SOLVERS = ("opensees", "pynite")
LIBRARY = ROOT / "verification" / "problems"


# ====================================================================================== references
def oracle_reference(model, cid) -> ForeignResult:
    """Nodal displacements, reactions and member end forces of a member line from the ODE oracle."""
    from benchfea.verify.tier_a import line_representation, solve_line

    line = line_representation(model)
    sol = solve_line(model, line, cid)
    s, sign, gdof, R = line["s"], line["sign"], line["gdof"], np.asarray(line["R"], float)
    prob = sol["problem"]
    n = len(model.nodes)
    U, Rx = np.zeros((n, 6)), np.zeros((n, 6))
    names = ("u", "v", "w", "tx", "ty", "tz")
    for k, nid in enumerate(line["nodes"]):
        side = -1 if k == len(line["nodes"]) - 1 else +1
        for a, name in enumerate(names):
            if a in (4, 5) and nid in prob["end_release"][a]:
                raise ValueError("released end rotation: not a node quantity of the oracle")
            sd = prob["follow"][a].get(nid, side) if a in (4, 5) else side
            U[model.nodes[nid].index, gdof[a]] = sign[a] * sol[name](s[nid], sd)
        if nid in sol["reactions"]:
            for a in range(6):
                Rx[model.nodes[nid].index, gdof[a]] = sign[a] * float(sol["reactions"][nid][a])
    ends = {}
    for mid in line["members"]:
        m = model.members[mid]
        xi, xj = s[m.i], s[m.j]
        fi = np.array([sol[q](xi, +1) for q in ("N", "Vy", "Vz", "T", "My", "Mz")], float)
        fj = np.array([sol[q](xj, -1) for q in ("N", "Vy", "Vz", "T", "My", "Mz")], float)
        ends[mid] = np.r_[-fi, fj]
    return ForeignResult("oracle", cid, True, Fidelity("oracle"), U, Rx, ends)


def tableau_reference(raw, model, cid, nsub=1) -> ForeignResult:
    from tableau import Tableau

    t = Tableau(copy.deepcopy(raw))
    t.nsub, t.arm = nsub, True
    ref = t.solve(cid, second_order=model.combinations[cid].analysis_type == "secondOrder")
    fl, ff = float(model.units.to_internal(1.0, "length")), float(model.units.to_internal(1.0, "force"))
    six = np.r_[[ff] * 3, [ff * fl] * 3]
    U = np.array(ref["displacements"], float) * np.r_[[fl] * 3, [1.0] * 3]
    Rx = np.array(ref["reactions"], float) * six
    ends = {mid: np.array(m["end_forces"], float) * np.r_[six, six] for mid, m in ref["members"].items()}
    return ForeignResult("tableau", cid, True, Fidelity("tableau"), U, Rx, ends)


# ====================================================================================== models
def _library(slug, pid):
    return json.loads((LIBRARY / slug / pid / "model.json").read_text())


def _plain_sections(raw):
    """The tableau solver reads sections of type ``properties`` and materials with E and G only: the
    catalogue / parametric sections of a generated model are replaced by the property values of the
    built model (what the foreign translators receive as well), the design block is dropped."""
    model = build_model(copy.deepcopy(raw))
    fl = float(model.units.to_internal(1.0, "length"))
    fp = float(model.units.to_internal(1.0, "force")) / fl ** 2
    raw.pop("design", None)
    for sid, sec in model.sections.items():
        d = {"type": "properties", "A": sec.A / fl ** 2, "Iy": sec.Iy / fl ** 4, "Iz": sec.Iz / fl ** 4,
             "J": sec.J / fl ** 4}
        if model.analysis.shear_deformation and sec.Asy and sec.Asz:
            d["Asy"], d["Asz"] = sec.Asy / fl ** 2, sec.Asz / fl ** 2
        raw["sections"][sid] = d
    for mid, mat in model.materials.items():
        old = raw["materials"][mid]
        new = {"E": mat.E / fp, "G": mat.G / fp}
        for k in ("unitWeight", "Fy", "Fu"):
            if k in old:
                new[k] = old[k]
        if "unitWeight" not in new and mat.unit_weight:
            new["unitWeight"] = mat.unit_weight / (float(model.units.to_internal(1.0, "force")) / fl ** 3)
        raw["materials"][mid] = new
    return raw


def _second_order_combination(raw):
    return next(c for c, v in raw["combinations"].items() if v.get("analysisType") == "secondOrder")


def cases(quick=False) -> list:
    """``(label, raw, combination, reference kind, tableau nsub)``"""
    import models as bank

    out = []
    derived = ["beam-column-02", "beam-column-05", "beam-column-08", "beam-column-tension-05"]
    for pid in derived[:2] if quick else derived:
        raw = _library("derived", pid)
        out.append((f"derived/{pid}", raw, _second_order_combination(raw), "oracle", None))
    aisc = sorted(p.name for p in (LIBRARY / "aisc-commentary").iterdir()
                  if p.is_dir() and "without-shear" in p.name and not p.name.endswith("-p0"))
    for pid in aisc[:2] + aisc[6:8] if quick else aisc:
        raw = _library("aisc-commentary", pid)
        out.append((f"aisc-commentary/{pid}", raw, _second_order_combination(raw), "oracle", None))
    # frames: tableau.  nsub = 1 where every member has a constant axial force (exact), otherwise 32 / 64
    frames = [("S01", 1), ("S03", 1), ("S04", 1), ("S06", 32), ("S07", 32), ("S08", 1)]
    for key, nsub in frames[:3] + frames[3:4] if quick else frames:
        name = next(k for k in bank.SECOND if k.startswith(key))
        raw, cid = bank.SECOND[name]()
        out.append((f"bank/{name[:3]}", raw, cid, "tableau", nsub))
    if not quick:
        from benchfea.bench.generate import generate
        for family, seed, level in (("portal_frame", 1, 3), ("portal_frame", 2, 3), ("multistorey_2d", 1, 3),
                                    ("multistorey_2d", 2, 4), ("pipe_rack", 1, 3), ("braced_frame", 1, 3),
                                    ("multistorey_3d", 1, 4), ("diaphragm_frame", 1, 4)):
            try:
                g = generate(family, seed, level, features={"secondOrder": True})
            except Exception as e:                                   # a family without that flag / level
                out.append((f"generated/{family}-{level}-{seed}", None, f"{type(e).__name__}: {e}", "skip", None))
                continue
            raw = _plain_sections(json.loads(json.dumps(g.document)))
            for cid, v in raw["combinations"].items():
                if v.get("analysisType") == "secondOrder":
                    out.append((f"generated/{family}-{level}-{seed}/{cid}", raw, cid, "tableau", 16))
                    break
    return out


# ====================================================================================== the table
def row(label, raw, cid, kind, nsub, solvers=SOLVERS) -> dict:
    model = build_model(copy.deepcopy(raw))
    out = {"model": label, "combination": cid, "reference": kind, "solvers": {}}
    try:
        if kind == "oracle":
            ref, slack = oracle_reference(model, cid), ORACLE_ACCURACY
        else:
            ref = tableau_reference(raw, model, cid, nsub)
            slack = 0.0
            if nsub > 1:
                slack = relative_difference(ref, tableau_reference(raw, model, cid, 2 * nsub), model)
                ref = tableau_reference(raw, model, cid, 2 * nsub)
                slack /= 3.0                 # mid-point rule: the finer run is 4x closer; (coarse - fine) / 3
    except Exception as e:                                           # the reference itself cannot take the model
        out["error"] = f"{type(e).__name__}: {str(e)[:160]}"
        return out
    out["referenceUncertainty"] = slack
    for s in solvers:
        r = run(s, model, cid, in_process=(s == "pynite"))
        rec = {"ok": r.ok}
        if r.ok:
            est = r.fidelity.estimated_error
            true = relative_difference(r, ref, model)
            conv = (r.info.get("convergence") or {}).get("subdivisions", {})
            rec.update(estimate=est, true=true, unsupported=[i.feature for i in r.fidelity.unsupported],
                       eligible=bool(est is not None and est <= 5e-3 and not r.fidelity.unsupported),
                       observedOrder=conv.get("observedOrder"), status=conv.get("status"),
                       bounded=None if est is None else bool(est + slack >= true),
                       ratio=None if (est is None or true <= 10 * slack or true == 0.0) else est / true)
        else:
            rec["error"] = r.error[:160]
        out["solvers"][s] = rec
    return out


def table(quick=False, progress=None) -> list:
    rows = []
    for label, raw, cid, kind, nsub in cases(quick):
        if kind == "skip":
            rows.append({"model": label, "error": cid, "solvers": {}})
        else:
            rows.append(row(label, raw, cid, kind, nsub))
        if progress:
            progress(rows[-1])
    return rows


def format_row(r) -> str:
    cells = [r["model"], r.get("reference", "-")]
    if "error" in r:
        return "| " + " | ".join(cells + ["reference not available: " + r["error"]] + [""] * 7) + " |"
    for s in SOLVERS:
        d = r["solvers"][s]
        if not d["ok"]:
            cells += ["failed", "", "", ""]
            continue
        if d["unsupported"]:
            cells += ["unsupported", f"{d['true']:.2e}", "", ""]
            continue
        est = "none" if d["estimate"] is None else f"{d['estimate']:.2e}"
        order = "-" if d["observedOrder"] is None else f"{d['observedOrder']:.2f}"
        ratio = "-" if d["ratio"] is None else f"{d['ratio']:.1f}"
        flag = "" if d["bounded"] in (True, None) else " **NOT A BOUND**"
        cells += [est + ("" if d["eligible"] else " (not eligible)") + flag, f"{d['true']:.2e}", ratio, order]
    return "| " + " | ".join(cells) + f" | {r['referenceUncertainty']:.0e} |"


HEADER = ("| model | reference | OpenSees estimate | true error | est/true | observed order | PyNite estimate | "
          "true error | est/true | observed order | reference uncertainty |\n" + "|---" * 11 + "|")

if __name__ == "__main__":
    quick = "--quick" in sys.argv
    print(HEADER, flush=True)
    rows = table(quick, progress=lambda r: print(format_row(r), flush=True))
    if "--json" in sys.argv:
        Path(sys.argv[sys.argv.index("--json") + 1]).write_text(json.dumps(rows, indent=1, default=float))
    bad = [(r["model"], s) for r in rows for s, d in r["solvers"].items()
           if d.get("ok") and not d["unsupported"] and d["bounded"] is False]
    print("\nestimate below the true error:", bad or "none")
