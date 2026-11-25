"""Tier D: invariants on every solve.

The solver stores, per combination: stiffness symmetry, zero energy under rigid-body motion ("rigid-body
modes"), positive definiteness after supports, linear solve residual, global force and moment
equilibrium, joint equilibrium, member equilibrium, member compatibility, work = strain energy (and the
deformed-axis member equilibrium of a second-order combination) -- definitions in docs/THEORY.md 8, 15.
This tier collects them and adds two runs of its own on the first-order operator of the model:

superposition   factored sum of the separately solved load cases = the combination solved directly;
reciprocity     Maxwell-Betti on a deterministic selection of DOF pairs.

Both are properties of a LINEAR system; for a model with one-sided members they do not hold and are
not run (recorded).  For second-order / reduced-stiffness combinations they are run on the first-order
operator (they test the stiffness matrix, the load vectors and the recovery, not the iteration).

station statics  the internal forces at EVERY result station recomputed here from statics (the station
                forces at the member start, the load definitions, and in second order the lever of the
                axial force on the stored deflections), and the station forces at the member ends
                against the end forces (``benchfea.verify.stations``; no solver code).
station         the displacements of the member axis, the bending moments and the shears at every
integration     interior station from an integration of the member differential equations between
                the first and the last station (initial value problems, scipy ``solve_ivp``).
                The stored invariants are computed by the solver during the solve; these two parts
                are computed from the finished reference solution, so that a station value that is
                wrong in the results (and that tier C, which compares node displacements, reactions
                and end forces, cannot see) fails the tier.
"""
from __future__ import annotations

import time

from ..solve.invariants import reciprocity_check, superposition_check
from .common import ONE_SIDED
from .policy import POLICY
from .record import TierResult
from .stations import station_deflections, station_statics

__all__ = ["tier_d", "REQUIRED_INVARIANTS"]

#: invariants that must be present in every combination result (a missing one fails the tier)
REQUIRED_INVARIANTS = ("stiffness symmetry", "rigid-body modes", "positive definiteness",
                       "global force equilibrium", "global moment equilibrium", "joint equilibrium")


def tier_d(model, results) -> TierResult:
    t0 = time.perf_counter()
    failed, largest, per_combo, missing = [], [], {}, []
    worst_by_name = {}
    for cid, c in results.combinations.items():
        names = {ch.name for ch in c.invariants.checks}
        for req in REQUIRED_INVARIANTS:
            if req not in names:
                missing.append(f"{cid}: {req}")
        rows = []
        for ch in c.invariants.checks:
            rows.append({"name": ch.name, "value": float(ch.value), "tolerance": float(ch.tolerance),
                         "passed": bool(ch.passed)})
            if not ch.passed:
                failed.append(f"{cid}: {ch.name} = {ch.value:.3e} (tolerance {ch.tolerance:.1e})")
            if ch.name != "positive definiteness":       # its value is an eigenvalue (larger = better)
                ratio = ch.value / ch.tolerance if ch.tolerance > 0 else float(ch.value > 0)
                if ch.name not in worst_by_name or ratio > worst_by_name[ch.name]["ratio"]:
                    worst_by_name[ch.name] = {"what": ch.name, "combination": cid, "value": float(ch.value),
                                              "tolerance": float(ch.tolerance), "ratio": float(ratio)}
        if not c.converged:
            failed.append(f"{cid}: the iteration did not converge")
        per_combo[cid] = {"converged": bool(c.converged), "checks": rows}
    extra = {}
    if any(m.type in ONE_SIDED for m in model.members.values()):
        extra["note"] = ("superposition and reciprocity not run: the model has tension-/compression-only members "
                         "(the system is not linear)")
    else:
        sup = superposition_check(model, POLICY["D"]["superposition"])
        rec = reciprocity_check(model, tolerance=POLICY["D"]["reciprocity"])
        for ch in list(sup) + [rec]:
            row = {"name": ch.name, "value": float(ch.value), "tolerance": float(ch.tolerance),
                   "passed": bool(ch.passed), "detail": ch.detail}
            extra.setdefault("linearity", []).append(row)
            key = ch.name.split(":")[0]
            ratio = ch.value / ch.tolerance
            if key not in worst_by_name or ratio > worst_by_name[key]["ratio"]:
                worst_by_name[key] = {"what": ch.name, "value": float(ch.value), "tolerance": float(ch.tolerance),
                                      "ratio": float(ratio)}
            if not ch.passed:
                failed.append(f"{ch.name} = {ch.value:.3e} (tolerance {ch.tolerance:.1e})")
        if not all(c.analysis_type == "firstOrder" for c in model.combinations.values()):
            extra["note"] = "superposition and reciprocity run on the first-order operator of the model"
    failed += [f"missing invariant {m}" for m in missing]
    # ---- station forces recomputed from statics on the finished results (no solver code)
    st = station_statics(model, results, POLICY["D"]["stationStatics"], POLICY["D"]["stationStaticsIterated"])
    extra["stationStatics"] = st
    # ---- station deflections, moments and shears by integrating the member equations (no solver code)
    sd = station_deflections(model, results, POLICY["D"]["stationIntegration"])
    extra["stationIntegration"] = sd
    only_not_evaluated = (not failed and not st["failed"] and not sd["failed"]
                          and bool(st["notEvaluated"] or sd["notEvaluated"]))
    failed += ["station statics: " + f for f in st["failed"]]
    failed += ["station statics, not evaluated: " + f for f in st["notEvaluated"]]
    failed += ["station integration: " + f for f in sd["failed"]]
    failed += ["station integration, not evaluated: " + f for f in sd["notEvaluated"]]
    if sd["largest"].get("value", 0.0) > 0.0:
        w = sd["largest"]
        worst_by_name["station integration"] = {
            "what": f"station integration ({w.get('quantity')}, {w.get('where')})", "combination": w.get("combination"),
            "value": float(w["value"]), "tolerance": float(w["tolerance"]), "ratio": float(w["value"] / w["tolerance"])}
    if st["largest"].get("value", 0.0) > 0.0:
        w = st["largest"]
        worst_by_name["station statics"] = {
            "what": f"station statics ({w.get('quantity')}, {w.get('where')})", "combination": w.get("combination"),
            "value": float(w["value"]), "tolerance": float(w["tolerance"]), "ratio": float(w["value"] / w["tolerance"])}
    largest = sorted(worst_by_name.values(), key=lambda d: -d["ratio"])[:6]
    passed = not failed
    n = sum(len(v["checks"]) for v in per_combo.values()) + len(extra.get("linearity", []))
    summary = (f"{n} invariant check(s) over {len(per_combo)} combination(s) and the statics of "
               f"{st['stations']} result station(s): "
               + ("all passed" if passed else "FAILED: " + "; ".join(failed[:4])))
    tr = TierResult("D", True, passed, summary, largest,
                    {**POLICY["D"], "stored": "tolerance stored with each invariant (docs/THEORY.md 8, 15)"},
                    {"combinations": per_combo, **extra, "failed": failed}, time.perf_counter() - t0)
    if only_not_evaluated:
        tr.category = "error"                    # a check that applies could not be carried out
    return tr
