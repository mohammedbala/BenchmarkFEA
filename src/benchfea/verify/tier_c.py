"""Tier C: the same model solved by independent solvers (OpenSeesPy, PyNite) through
``benchfea.crosscheck``.

A combination counts as cross-checked only if AT LEAST TWO foreign solvers represent it -- fidelity
"exact", or "approximate" with an error bound (from runs of the foreign solver alone) not above ``maxApproximation`` -- and agree
with the engine and with each other within the tolerance plus the estimated errors of the results
compared (``benchfea.crosscheck.cross_check``).  Otherwise the tier is not passed and the reason says
how many solvers can represent the model and why the others cannot.  A model that is not cross-checked
is not publishable as verified.

Quantities that the engine DERIVES from its own solution and that a foreign solver cannot recompute are
handled as follows (docs/VERIFICATION.md, "Substitutions in tier C"); each use is written in the record:

notional loads   recomputed here from the load definitions (``crosscheck.prepare`` + the lever rule of
                 docs/THEORY.md 14, no solver code), compared with the loads the engine recorded, then
                 given to the foreign solvers as nodal loads.
computed tau_b   AISC 360-22 Eq. C2-2a/b re-evaluated here from the engine's station axial forces and
                 compared with the tau_b the engine used; the foreign solvers then run with these
                 values.  Agreement of their member forces with the engine's closes the loop: the
                 state satisfies equilibrium at the reduced stiffness (foreign solvers) and the
                 stiffness satisfies Eq. C2-2 at that state (here).
one-sided        tension-/compression-only members: the foreign solvers get the LINEAR model of the
members          engine's final active set (active members as truss members, slack members removed).
                 The solution of the one-sided problem is characterised by: equilibrium of that linear
                 model (foreign solvers), every active member carrying a force of its permitted sign,
                 and every slack member being strained in its forbidden sense (both evaluated here
                 from the FOREIGN results).  This is exact, not an approximation.  Member loads of a
                 slack member (e.g. its self weight) go to its two end nodes as the reactions of the
                 member held at both ends with its bending released there (docs/THEORY.md 13: "the
                 reactions of a simply supported beam"; the axial component by the lever rule of a
                 uniform bar).  They are computed here from the load definitions with the force-method
                 routine of ``crosscheck.fef`` (no solver code) and given to the foreign solvers as
                 nodal loads (``slack_member_loads``).

Not reproduced: linear buckling (neither translator computes buckling factors).  A model whose
reference solution contains buckling results therefore cannot pass tier C.
"""
from __future__ import annotations

import copy
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from ..crosscheck import cross_check
from ..crosscheck.compare import Tolerances
from ..crosscheck.fidelity import Fidelity
from ..crosscheck.foreign import ForeignResult
from ..crosscheck.prepare import prepare
from ..crosscheck.runner import available, run
from ..model import build_model
from ..problems import ModelInvalid
from .common import ONE_SIDED, document_hash
from .policy import POLICY
from .record import TierResult

__all__ = ["tier_c", "ForeignCache", "foreign_versions", "SOLVERS", "independent_notional_loads",
           "independent_tau_b", "active_set_document", "slack_member_loads"]

SOLVERS = ("opensees", "pynite")
_DIST = {"opensees": "openseespy", "pynite": "PyNiteFEA"}
RUN_OPTIONS = {"opensees": {}, "pynite": {"in_process": True}}      # PyNite is pure Python


def foreign_versions() -> dict:
    import importlib.metadata as md
    out = {}
    for s in SOLVERS:
        try:
            out[s] = md.version(_DIST[s]) if available(s) else None
        except md.PackageNotFoundError:
            out[s] = None
    return out


# ====================================================================================== cache
class ForeignCache:
    """Foreign solves keyed by (document hash, combination, solver, options, solver version).

    In memory, and on disk when ``directory`` is given (one JSON file per solve).  A cached result is
    exactly what the solver returned for that document; nothing about the engine is cached."""

    def __init__(self, directory=None, max_items: int = 4096):
        self.directory = Path(directory) if directory else None
        self.max_items = max_items
        self._mem: dict = {}
        self.hits = self.misses = 0
        if self.directory:
            self.directory.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def key(doc_hash: str, combination: str, solver: str, options: dict, version) -> str:
        blob = json.dumps([doc_hash, combination, solver, options, version], sort_keys=True, default=_jsonable)
        return hashlib.sha256(blob.encode()).hexdigest()

    def get(self, key: str):
        if key in self._mem:
            self.hits += 1
            return _from_dict(self._mem[key])
        if self.directory:
            f = self.directory / f"{key}.json"
            if f.exists():
                try:
                    d = json.loads(f.read_text())
                except (OSError, ValueError):
                    return None
                self._mem[key] = d
                self.hits += 1
                return _from_dict(d)
        return None

    def put(self, key: str, result: ForeignResult) -> None:
        d = json.loads(json.dumps(result.to_dict(), default=_jsonable))
        if len(self._mem) >= self.max_items:
            self._mem.pop(next(iter(self._mem)))
        self._mem[key] = d
        self.misses += 1
        if self.directory:
            try:
                (self.directory / f"{key}.json").write_text(json.dumps(d))
            except OSError:
                pass


def _jsonable(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        return o.item()
    return str(o)


def _from_dict(d: dict) -> ForeignResult:
    return ForeignResult(d["solver"], d["combination"], d["ok"], Fidelity.from_dict(d["fidelity"]),
                         None if d["displacements"] is None else np.array(d["displacements"], float),
                         None if d["reactions"] is None else np.array(d["reactions"], float),
                         {k: np.array(v, float) for k, v in d["endForces"].items()}, d.get("error", ""),
                         copy.deepcopy(d.get("info", {})))


DEFAULT_CACHE = ForeignCache()


# ====================================================================================== derived quantities
def independent_notional_loads(model, cid: str) -> dict:
    """Notional loads of a combination from the load definitions, per node, internal units and at the
    level of the unfactored-by-alpha loads: ``coefficient x Y x direction`` with Y = the downward
    component of the loads at the node, member loads carried to the two member ends by the lever rule
    (docs/THEORY.md 14; AISC 360-22 C2.2b Eq. C2-1 with the engine's lumping rule)."""
    combo = model.combinations[cid]
    p = prepare(model, cid, extra_nodal_loads={})
    up = np.asarray(model.vertical, float)
    Y = np.zeros(len(p.node_ids))
    for k in range(len(p.node_ids)):
        Y[k] -= float(p.load[k, :3] @ up)
    for m in p.members.values():
        i, j = p.index(m.i), p.index(m.j)
        up_local = m.R @ up

        def put(x, force_local):
            y = -float(np.asarray(force_local, float) @ up_local)
            Y[i] += y * (1.0 - x / m.L)
            Y[j] += y * (x / m.L)
        for kind, vec, x in m.points:
            if kind == "force":
                put(min(max(x, 0.0), m.L), vec)
        for w1, w2, x1, x2 in m.dists:
            h = x2 - x1
            put(x1 + h / 3.0, np.asarray(w1) * h / 2.0)
            put(x1 + 2.0 * h / 3.0, np.asarray(w2) * h / 2.0)
    d = np.asarray(combo.notional.direction, float)
    return {nid: (combo.notional.coefficient * Y[k] / p.alpha) * d for k, nid in enumerate(p.node_ids) if Y[k] != 0.0}


def independent_tau_b(model, cid: str, cres) -> dict:
    """AISC 360-22 Eq. C2-2a / C2-2b evaluated from the station axial forces of a combination result:
    ratio = alpha P_r / P_ns with alpha P_r = the largest compression along the member in the analysis
    (results are stored divided by alpha) and P_ns = F_y A (engine rule, docs/THEORY.md 14)."""
    an = model.analysis
    alpha = float(an.alpha) if model.combinations[cid].analysis_type == "secondOrder" else 1.0
    members = model.members if an.reduction_members is None else an.reduction_members
    out = {}
    for mid in members:
        m = model.members[mid]
        fy = model.materials[m.material].Fy
        comp = max(0.0, -float(np.min(cres.members[mid].forces[:, 0]))) * alpha
        ratio = comp / (fy * model.sections[m.section].A)
        out[mid] = 1.0 if ratio <= 0.5 else 4.0 * ratio * (1.0 - ratio)
    return out


def slack_member_loads(model, cid: str, slack) -> dict:
    """Nodal loads that stand for the member loads of the slack members of a combination:
    ``{node: [Fx, Fy, Fz, Mx, My, Mz]}`` in internal units, at the level of the loads NOT multiplied by
    alpha (the form of ``extra_nodal_loads``).

    Rule (docs/THEORY.md 13): a slack member has no stiffness and hands its loads to its end nodes as
    the reactions of a simply supported beam -- the member held at both ends of its flexible length with
    both bending rotations released; a force along the axis is shared by the lever rule (uniform bar), a
    point moment becomes a couple of transverse end forces.  The end reactions come from
    ``crosscheck.fef.fixed_end_forces`` (force method; shares nothing with ``benchfea.solve``), the load
    intensities from ``crosscheck.prepare``; a rigid end offset carries the end force to the node with
    its moment ``offset x F``."""
    from ..crosscheck.fef import fixed_end_forces
    p = prepare(model, cid, extra_nodal_loads={})
    out = {}
    for mid in slack:
        m = p.members[mid]
        if not m.loaded:
            continue
        ei_y, ei_z = m.EIy or m.E * m.A ** 2 / 12.0, m.EIz or m.E * m.A ** 2 / 12.0   # statically determinate: any value
        f = fixed_end_forces(m.L, ei_y, ei_z, None, None, m.points, m.dists, ky=(0.0, 0.0), kz=(0.0, 0.0))
        for blk, nid, off in ((0, m.i, m.offset_i), (6, m.j, m.offset_j)):
            F = -m.R.T @ f[blk:blk + 3]
            M = -m.R.T @ f[blk + 3:blk + 6] + np.cross(off, F)
            v = out.setdefault(nid, np.zeros(6))
            v[:3] += F / p.alpha
            v[3:] += M / p.alpha
    return {nid: [float(x) for x in v] for nid, v in out.items() if np.any(v != 0.0)}


def active_set_document(raw: dict, model, cres):
    """The linear model of the engine's final active set of one combination: active one-sided members
    become truss members, slack ones are removed together with their member loads (which
    :func:`slack_member_loads` turns into nodal loads).  Returns ``(document, slack ids, ids of the slack
    members that carry member loads or None)``."""
    slack = [mid for mid, m in model.members.items() if m.type in ONE_SIDED and not cres.members[mid].active]
    doc = copy.deepcopy(raw)
    doc.pop("design", None)
    loaded = set()
    for case in doc.get("loadCases", {}).values():
        keep = []
        for ld in case.get("loads", []):
            if ld["type"] in ("memberPoint", "memberDistributed") and ld["member"] in slack:
                loaded.add(ld["member"])
                continue
            if ld["type"] == "selfWeight":
                ld = dict(ld)
                listed = ld.get("members", list(raw["members"]))
                for mid in listed:
                    if mid in slack and model.materials[model.members[mid].material].unit_weight != 0.0:
                        loaded.add(mid)
                ld["members"] = [mid for mid in listed if mid not in slack]
                if not ld["members"]:
                    continue
            keep.append(ld)
        case["loads"] = keep
    for mid, m in model.members.items():
        if m.type in ONE_SIDED:
            if mid in slack:
                del doc["members"][mid]
            else:
                doc["members"][mid]["type"] = "truss"
    sr = doc.get("analysis", {}).get("stiffnessReduction", {})
    if "members" in sr:
        sr["members"] = [x for x in sr["members"] if x not in slack]
    for case in doc.get("loadCases", {}).values():
        if not case["loads"]:
            # a load case must keep a load: a zero nodal load keeps the document valid
            case["loads"] = [{"type": "nodal", "node": next(iter(doc["nodes"])), "Fx": 0.0}]
    return doc, slack, (sorted(loaded) or None)


def _complementarity(model, slack, result: ForeignResult) -> list:
    """Sign conditions of the one-sided members evaluated from a foreign result of the active-set
    model; returns the violations as text."""
    out = []
    u = result.displacements
    scale = max(float(np.abs(u[:, :3]).max()), 1e-300)
    fscale = max((float(np.abs(v[[0, 6]]).max()) for v in result.end_forces.values()), default=0.0)
    for mid, m in model.members.items():
        if m.type not in ONE_SIDED:
            continue
        sign = 1.0 if m.type == "tensionOnly" else -1.0
        if mid in slack:
            ui, uj = u[model.nodes[m.i].index], u[model.nodes[m.j].index]
            ei = ui[:3] + np.cross(ui[3:], m.offset_i)
            ej = uj[:3] + np.cross(uj[3:], m.offset_j)
            elong = float((ej - ei) @ np.asarray(m.R)[0])
            if sign * elong > 1e-7 * scale:
                out.append(f"{result.solver}: slack {m.type} member {mid} would be "
                           f"{'stretched' if sign > 0 else 'shortened'} by {abs(elong):.3e} mm (it should be active)")
        else:
            n_force = float(result.end_forces[mid][6])
            if sign * n_force < -1e-7 * max(fscale, 1e-300):
                out.append(f"{result.solver}: active {m.type} member {mid} carries N = {n_force:.6g} N of the "
                           "forbidden sign (it should be slack)")
    return out


# ====================================================================================== the tier
TRANSLATOR_REVISION = 2      # part of the cache key: results of an earlier translator / estimator are not reused


def _foreign(doc_hash, model, cid, solver, options, cache, timeout, versions):
    key = cache.key(doc_hash, cid, solver, options, [versions.get(solver), TRANSLATOR_REVISION])
    hit = cache.get(key)
    if hit is not None:
        return hit, True
    res = run(solver, model, cid, timeout=timeout, **RUN_OPTIONS.get(solver, {}), **options)
    if res.ok or "did not finish" not in (res.error or ""):
        cache.put(key, res)                      # a timeout is not cached (it depends on the machine)
    return res, False


def tier_c(raw: dict, model, results, *, cache: ForeignCache | None = None, timeout: float = 120.0,
           deadline: float | None = None, solvers=SOLVERS) -> TierResult:
    t0 = time.perf_counter()
    cache = cache if cache is not None else DEFAULT_CACHE
    pol = POLICY["C"]
    versions = foreign_versions()
    tol = Tolerances(*([pol["relative"]] * 6))
    doc_hash = document_hash(raw)
    per_combo, failed, largest, timings = {}, [], [], {}
    disagreements, not_run = [], []
    one_sided = any(m.type in ONE_SIDED for m in model.members.values())
    for cid, cres in results.combinations.items():
        t1 = time.perf_counter()
        info = {"substitutions": []}
        per_combo[cid] = info
        if deadline is not None and time.perf_counter() > deadline:
            info.update(crossChecked=False, reasons=["not run: the time budget of verify_model was used up"])
            failed.append(f"{cid}: not run (time budget)")
            not_run.append(cid)
            continue
        combo = model.combinations[cid]
        options, problems = {}, []
        c_model, c_hash = model, doc_hash
        # ---- derived quantities
        if combo.notional is not None:
            mine = independent_notional_loads(model, cid)
            rec = (cres.solver or {}).get("notionalLoads", {}).get("nodes", {})
            alpha = float(model.analysis.alpha) if combo.analysis_type == "secondOrder" else 1.0
            fscale = max((float(np.abs(v).max()) for v in mine.values()), default=0.0)
            worst = 0.0
            for nid in set(mine) | set(rec):
                a = np.asarray(mine.get(nid, np.zeros(3)), float) * alpha
                b = np.array([float(model.units.to_internal(v, "force")) for v in rec.get(nid, [0.0, 0.0, 0.0])])
                worst = max(worst, float(np.abs(a - b).max()) / fscale if fscale else float(np.abs(b).max()))
            info["substitutions"].append({
                "what": "notional loads", "how": "recomputed from the load definitions and given to the foreign "
                                                 "solvers as nodal loads",
                "differenceToEngine": worst, "tolerance": 1e-9})
            if worst > 1e-9:
                problems.append(f"the notional loads recomputed from the load definitions differ from the engine's "
                                f"by {worst:.2e}")
            options["extra_nodal_loads"] = {nid: [float(x) for x in v] + [0.0, 0.0, 0.0] for nid, v in mine.items()}
        if combo.stiffness_reduction and model.analysis.tau_b == "computed":
            mine = independent_tau_b(model, cid, cres)
            worst = max((abs(mine[mid] - cres.tau_b.get(mid, 1.0)) for mid in mine), default=0.0)
            tb_tol = 10.0 * float(model.analysis.convergence_tolerance)
            info["substitutions"].append({
                "what": "tau_b", "how": "AISC 360-22 Eq. C2-2a/b re-evaluated from the engine's station axial "
                                        "forces; the foreign solvers run with the engine's values",
                "differenceToEngine": worst, "tolerance": tb_tol})
            if worst > tb_tol:
                problems.append(f"tau_b re-evaluated from the station axial forces differs from the engine's by "
                                f"{worst:.2e}")
            options["tau_b"] = {mid: float(v) for mid, v in cres.tau_b.items()}
        slack, blocked = [], False
        if one_sided:
            doc, slack, loaded_slack = active_set_document(raw, model, cres)
            try:
                c_model, c_hash = build_model(copy.deepcopy(doc)), document_hash(doc)
            except ModelInvalid as e:
                problems.append("the active-set substitute is not a valid model: "
                                + "; ".join(p.message for p in e.problems[:2]))
                blocked = True
            sub = {"what": "one-sided members", "slack": slack,
                   "how": "foreign solvers get the linear model of the engine's final active set; the sign conditions "
                          "of every one-sided member are then evaluated from the foreign results"}
            if loaded_slack and not blocked:
                # only the loads of this combination's slack members that act in this combination
                extra = slack_member_loads(model, cid, loaded_slack)
                merged = {nid: list(v) for nid, v in options.get("extra_nodal_loads", {}).items()}
                for nid, v in extra.items():
                    old = merged.get(nid, [0.0] * 6)
                    merged[nid] = [float(a + b) for a, b in zip(old, v)]
                if merged:
                    options["extra_nodal_loads"] = merged
                sub.update(slackMembersWithLoads=loaded_slack,
                           slackMemberLoads="member loads of the slack members given to the foreign solvers as the "
                                            "end reactions of a simply supported beam (docs/THEORY.md 13), computed "
                                            "from the load definitions")
            info["substitutions"].append(sub)
        if blocked:
            disagreements += [f"{cid}: {x}" for x in problems if "differ" in x]
            info.update(crossChecked=False, reasons=problems, eligible=[])
            failed.append(f"{cid}: " + problems[0])
            continue
        # ---- foreign solves (cached)
        foreign, cached = {}, {}
        opt_key = {k: options[k] for k in sorted(options)}
        for s in solvers:
            foreign[s], cached[s] = _foreign(c_hash, c_model, cid, s, opt_key, cache, timeout, versions)
        report = cross_check(c_model, cid, benchfea_result=cres, solvers=tuple(solvers), tolerances=tol,
                             max_approximation=pol["maxApproximation"], foreign=foreign,
                             options={s: dict(options) for s in solvers})
        # decision: crosscheck's eligibility rule; agreement within the base tolerance plus the estimated
        # errors of the two results compared (``estimateSafety`` = 1.0: exactly crosscheck's rule).  The
        # estimates are bounds computed from three refinement levels of the foreign solver itself
        # (crosscheck.runner.estimate_from_levels; docs/CROSSCHECK.md 1.3) -- never from the engine's result
        base, safety = pol["relative"], pol["estimateSafety"]
        reasons = [r for r in report.reasons if " differ: " not in r and not r.endswith("agree")
                   and "represent the model and agree" not in r] + problems
        comps, agree = [], True
        for c in report.comparisons:
            row = {"a": c.a, "b": c.b, "passed": True, "maxRelative": c.max_relative}
            ratio_worst = -1.0
            for qn, d in c.worst.items():
                allowed = base + safety * (getattr(c.tolerances, qn) - base)
                if d.relative / allowed > ratio_worst:
                    ratio_worst = d.relative / allowed
                    row.update(where=d.where, quantity=qn, value=d.relative, tolerance=allowed)
                if not d.relative <= allowed:
                    row["passed"] = False
            if not row["passed"]:
                agree = False
                disagreements.append(f"{cid}: {c.a} vs {c.b}")
                reasons.append(f"{c.a} and {c.b} differ: {row['quantity']} {row['value']:.2e} at {row['where']} "
                               f"(allowed {row['tolerance']:.2e})")
            comps.append(row)
            largest.append({"what": f"{c.a} vs {c.b}", "combination": cid, "value": row["value"],
                            "tolerance": row["tolerance"], "where": row["where"]})
        ok = len(report.eligible) >= pol["requiredSolvers"] and agree and not problems
        disagreements += [f"{cid}: {x}" for x in problems]
        if ok:
            reasons.append(" and ".join(report.eligible) + " represent the model and agree with each other and "
                           "with the engine")
        if ok and one_sided:
            viol = [v for s in report.eligible for v in _complementarity(model, slack, foreign[s])]
            if viol:
                ok, reasons = False, reasons + viol
                disagreements += [f"{cid}: {v}" for v in viol]
            info["complementarity"] = {"violations": viol}
        info.update(
            crossChecked=bool(ok), eligible=list(report.eligible), reasons=reasons, comparisons=comps,
            equilibrium=report.equilibrium,
            solvers={s: {"ok": r.ok, "error": r.error, "level": r.fidelity.level,
                         "estimatedError": r.fidelity.estimated_error,
                         "unsupported": [f"{i.feature}: {i.note}" for i in r.fidelity.unsupported],
                         "approximate": [{"feature": i.feature, "error": i.error, "note": i.note}
                                         for i in r.fidelity.approximate]}
                     for s, r in foreign.items()})
        if not ok:
            n = len(report.eligible)
            if n < pol["requiredSolvers"]:
                why = [r for r in reasons if not r.startswith("only ")]
                failed.append(f"{cid}: only {n} solver(s) can represent: " + "; ".join(why)[:600])
            else:
                failed.append(f"{cid}: " + "; ".join(r for r in reasons if not r.startswith("only "))[:600])
        timings[cid] = {"seconds": round(time.perf_counter() - t1, 4), "cached": cached}
    buckling = {}
    for cid in results.buckling:
        buckling[cid] = {"crossChecked": False,
                         "reason": "linear buckling: only 0 solver(s) can represent (OpenSees translator: buckling "
                                   "not implemented; PyNite translator: unsupported)"}
        failed.append(f"buckling of {cid}: only 0 solver(s) can represent linear buckling")
    largest.sort(key=lambda d: -d["value"] / d["tolerance"])
    passed = not failed
    n = len(per_combo)
    summary = (f"{n} combination(s) cross-checked by OpenSeesPy and PyNite" if passed
               else "not passed (" + "; ".join(failed[:3]) + ")")
    tr = TierResult("C", True, passed, summary, largest[:6],
                    {"relative": pol["relative"], "maxApproximation": pol["maxApproximation"],
                     "requiredSolvers": pol["requiredSolvers"], "estimateSafety": pol["estimateSafety"],
                     "rule": "relative + (estimated errors of the two results compared), the rule of "
                             "benchfea.crosscheck.cross_check"},
                    {"combinations": per_combo, "buckling": buckling, "failed": failed,
                     "foreignVersions": versions}, time.perf_counter() - t0)
    tr.timings = timings
    # why the tier is not passed: two computations differ (the reference solution is in doubt), or too few
    # independent solvers can take the model (it is merely not cross-checked)
    tr.category = "" if passed else ("disagreement" if disagreements else ("not-run" if not_run and len(not_run) == len(
        [f for f in failed]) else "not-representable"))
    tr.details["disagreements"] = disagreements
    return tr
