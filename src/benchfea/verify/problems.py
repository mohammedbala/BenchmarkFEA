"""The verification problem library: loader, targets, ``verify_problem``, ``verify_library``.

Layout (docs/VERIFICATION.md, "Problem library")::

    verification/problems/<source-slug>/<problem-id>/
        model.json      a normal benchfea input file
        problem.json    {"id", "title", "tier": "A" | "B",
                         "source": {"work", "edition", "location", "pages", "notes"},
                         "basis": "published" | "derived-here" | "recomputed-from-published-statement",
                         "analysis": text,
                         "targets": [{"path": "/combinations/L/displacements/B/uz"  |  "extract": name, ...args,
                                      "value": number in the MODEL FILE's units, "unit": text,
                                      "tolerance": {"rel": x} | {"abs": y},
                                      "precisionNote": text, "sign": text}],
                         "status": "active" | "quarantined", "quarantine": {"diagnosis": text},
                         "limitsNote": text}

A target is either a path into the RESULTS document (the one ``benchfea run`` writes: keys separated by
"/", list items by index; "~1" stands for "/" and "~0" for "~" inside a key) or a named extractor
(``EXTRACTORS``).  ``"document": "design"`` makes a path refer to the design document.

Outcome of a problem (``ProblemRecord.outcome``):

passed       declared active; every target reproduced within its tolerance, and no tier that ran on the
             model shows a disagreement (tier D, tier E, the automatic part of tier A, tier C between the
             solvers that can represent the model).  Whether the model is also cross-checked by two
             foreign solvers is reported in the model record (a model only one solver can represent is
             "not cross-checked": that is a statement about the foreign solvers, not a failure of the
             problem, and it is printed in the report).
failed       declared active but a target is missed or a tier shows a disagreement.  The author must
             fix the problem or quarantine it with a diagnosis; ``benchfea verify`` exits nonzero.
quarantined  declared quarantined (targets are still evaluated and reported, never counted as passed).
invalid      the directory could not be loaded.
"""
from __future__ import annotations

import json
import math
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

from .policy import POLICY, check_target_tolerance
from .record import LibraryReport, ProblemRecord, TargetResult

__all__ = ["ProblemSpec", "ProblemError", "load_problem", "find_problems", "library_root", "evaluate_targets",
           "verify_problem", "verify_library", "EXTRACTORS", "resolve_path"]

BASES = ("published", "derived-here", "recomputed-from-published-statement")
STATION_QUANTITIES = ("N", "Vy", "Vz", "T", "My", "Mz", "ux", "uy", "uz")


class ProblemError(Exception):
    """A problem directory that does not follow the format; ``errors`` lists every finding."""

    def __init__(self, errors):
        super().__init__("; ".join(errors))
        self.errors = list(errors)


class TargetError(Exception):
    """A target that cannot be evaluated on the results (the message says why)."""


@dataclass
class ProblemSpec:
    id: str
    directory: Path
    title: str
    tier: str
    source: dict
    basis: str
    analysis: str
    targets: list
    status: str
    quarantine: dict | None = None
    limits_note: str = ""
    notes: list = field(default_factory=list)          # policy remarks found while loading

    @property
    def model_path(self) -> Path:
        return self.directory / "model.json"

    @property
    def needs_design(self) -> bool:
        return any(t.get("document") == "design" or _norm(t.get("extract", "")) in ("designcheck",)
                   for t in self.targets)


# ====================================================================================== loading
def library_root(start=None) -> Path | None:
    """``verification/problems`` of the repository this package (or ``start`` / the working directory)
    belongs to, or None."""
    for base in [Path(start)] if start else [Path.cwd(), Path(__file__).resolve()]:
        for p in [base, *base.parents]:
            cand = p / "verification" / "problems"
            if cand.is_dir():
                return cand
    return None


def find_problems(root=None) -> list:
    """Problem directories (those that hold a problem.json), sorted."""
    root = Path(root) if root else library_root()
    if root is None or not root.is_dir():
        return []
    return sorted(p.parent for p in root.rglob("problem.json"))


def load_problem(directory) -> ProblemSpec:
    d = Path(directory)
    errors = []
    pj, mj = d / "problem.json", d / "model.json"
    if not pj.is_file():
        raise ProblemError([f"{pj}: missing"])
    if not mj.is_file():
        errors.append(f"{mj}: missing")
    try:
        data = json.loads(pj.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        raise ProblemError([f"{pj}: not readable JSON ({e})"]) from e
    if not isinstance(data, dict):
        raise ProblemError([f"{pj}: must be a JSON object"])
    for key in ("id", "title", "tier", "source", "basis", "targets", "status"):
        if key not in data:
            errors.append(f"problem.json: '{key}' is required")
    if errors:
        raise ProblemError(errors)
    if data["id"] != d.name:
        errors.append(f"problem.json: id {data['id']!r} must equal the directory name {d.name!r}")
    if data["tier"] not in ("A", "B"):
        errors.append("problem.json: tier must be 'A' or 'B'")
    if data["basis"] not in BASES:
        errors.append(f"problem.json: basis must be one of {', '.join(BASES)}")
    if data["status"] not in ("active", "quarantined"):
        errors.append("problem.json: status must be 'active' or 'quarantined'")
    src = data["source"]
    if not isinstance(src, dict) or not src.get("work"):
        errors.append("problem.json: source must be an object with at least 'work'")
    notes = []
    if data.get("tier") == "B" and data.get("basis") == "derived-here":
        errors.append("problem.json: a tier-B problem reproduces a PUBLISHED problem; basis 'derived-here' is tier A")
    if data.get("tier") == "B" and isinstance(src, dict) and not src.get("location"):
        errors.append("problem.json: a published problem needs source.location (example / figure / equation number)")
    if data["status"] == "quarantined" and not (isinstance(data.get("quarantine"), dict)
                                                and str(data["quarantine"].get("diagnosis", "")).strip()):
        errors.append("problem.json: a quarantined problem needs quarantine.diagnosis")
    targets = data["targets"]
    if not isinstance(targets, list) or not targets:
        errors.append("problem.json: targets must be a non-empty list")
        targets = []
    for k, t in enumerate(targets):
        where = f"problem.json: targets[{k}]"
        if not isinstance(t, dict):
            errors.append(f"{where}: must be an object")
            continue
        if ("path" in t) == ("extract" in t):
            errors.append(f"{where}: exactly one of 'path' and 'extract' is required")
        if "extract" in t and _norm(str(t["extract"])) not in _EXTRACTOR_INDEX:
            errors.append(f"{where}: unknown extractor {t['extract']!r} (known: {', '.join(EXTRACTORS)})")
        v = t.get("value")
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
            errors.append(f"{where}: 'value' must be a finite number")
            continue
        problem = check_target_tolerance(t.get("tolerance"))
        if problem:
            errors.append(f"{where}: {problem}")
            continue
        if "rel" in t["tolerance"] and v == 0.0:
            errors.append(f"{where}: a relative tolerance on the value 0 can never be met; use {{\"abs\": y}}")
        if data.get("basis") == "derived-here" and "rel" in t["tolerance"] \
                and t["tolerance"]["rel"] > POLICY["targets"]["derivedHere"]:
            notes.append(f"targets[{k}]: relative tolerance {t['tolerance']['rel']:g} is looser than the policy for "
                         f"targets derived here ({POLICY['targets']['derivedHere']:g})")
        if data.get("basis") != "derived-here" and not str(t.get("precisionNote", "")).strip():
            notes.append(f"targets[{k}]: no precisionNote (the printed precision or the stated tolerance of the "
                         "source that the tolerance is derived from)")
    if errors:
        raise ProblemError(errors)
    return ProblemSpec(data["id"], d, data["title"], data["tier"], dict(src), data["basis"],
                       str(data.get("analysis", "")), targets, data["status"], data.get("quarantine"),
                       str(data.get("limitsNote", "")), notes)


# ====================================================================================== targets
def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def resolve_path(document, path: str):
    """Value at a JSON-pointer-like path ("/combinations/L/displacements/B/uz")."""
    if not path.startswith("/"):
        raise TargetError(f"path {path!r} must start with '/'")
    cur = document
    for raw in path[1:].split("/"):
        key = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(cur, dict):
            if key not in cur:
                near = ", ".join(list(cur)[:8])
                raise TargetError(f"path {path!r}: no key {key!r} (available: {near})")
            cur = cur[key]
        elif isinstance(cur, list):
            try:
                cur = cur[int(key)]
            except (ValueError, IndexError):
                raise TargetError(f"path {path!r}: {key!r} is not an index of a list of {len(cur)}") from None
        else:
            raise TargetError(f"path {path!r}: cannot descend into a {type(cur).__name__} at {key!r}")
    if isinstance(cur, bool) or not isinstance(cur, (int, float)):
        raise TargetError(f"path {path!r} is not a number")
    return float(cur)


def _combo(doc, t):
    cid = t.get("combination")
    if cid is None:
        if len(doc["combinations"]) == 1:
            cid = next(iter(doc["combinations"]))
        else:
            raise TargetError("'combination' is required (the results hold several)")
    if cid not in doc["combinations"]:
        raise TargetError(f"no combination {cid!r} in the results")
    return doc["combinations"][cid]


def _need(t, *names):
    for n in names:
        if n not in t:
            raise TargetError(f"extractor {t['extract']!r} needs '{n}'")


def _stations(doc, t):
    _need(t, "member", "quantity")
    c = _combo(doc, t)
    if t["member"] not in c["members"]:
        raise TargetError(f"no member {t['member']!r}")
    st = c["members"][t["member"]]["stations"]
    if t["quantity"] not in STATION_QUANTITIES:
        raise TargetError(f"quantity must be one of {', '.join(STATION_QUANTITIES)}")
    return st


def _x_buckling_factor(doc, design, t):
    """{"extract": "bucklingFactor", "combination": id, "mode": 1}: positive buckling factor (1 = lowest)."""
    cid = t.get("combination") or (next(iter(doc["buckling"])) if len(doc.get("buckling", {})) == 1 else None)
    if cid is None or cid not in doc.get("buckling", {}):
        raise TargetError(f"no buckling analysis of combination {cid!r} in the results")
    mode = int(t.get("mode", 1))
    f = doc["buckling"][cid]["factors"]
    if not 1 <= mode <= len(f):
        raise TargetError(f"mode {mode}: the results hold {len(f)} buckling factor(s)")
    return float(f[mode - 1])


def _x_member_station(doc, design, t):
    """{"extract": "memberStation", "combination", "member", "quantity": N|Vy|Vz|T|My|Mz|ux|uy|uz,
    "x": distance from the flexible start (file units) | "relative": 0..1, "side": "left" | "right"}:
    value AT A RESULT STATION (no interpolation: put a station there with analysis.resultStations or a
    load).  At a point load a station exists twice: "left" (default) is the value just before it."""
    st = _stations(doc, t)
    xs = st["x"]
    L = xs[-1]
    if ("x" in t) == ("relative" in t):
        raise TargetError("exactly one of 'x' and 'relative' is required")
    x = float(t["x"]) if "x" in t else float(t["relative"]) * L
    hits = [k for k, v in enumerate(xs) if abs(v - x) <= 1e-9 * max(L, 1e-300)]
    if not hits:
        raise TargetError(f"x = {x:g} is not a result station of member {t['member']!r} (stations: "
                          + ", ".join(f"{v:g}" for v in xs[:12]) + ("..." if len(xs) > 12 else "") + ")")
    side = t.get("side", "left")
    if side not in ("left", "right"):
        raise TargetError("side must be 'left' or 'right'")
    return float(st[t["quantity"]][hits[0] if side == "left" else hits[-1]])


def _x_member_extreme(doc, design, t):
    """{"extract": "memberExtreme", "combination", "member", "quantity", "kind": "max" | "min" | "absmax"}:
    extreme over the RESULT STATIONS of the member."""
    v = _stations(doc, t)[t["quantity"]]
    kind = t.get("kind", "absmax")
    if kind == "max":
        return float(max(v))
    if kind == "min":
        return float(min(v))
    if kind == "absmax":
        return float(max(abs(x) for x in v))
    raise TargetError("kind must be 'max', 'min' or 'absmax'")


def _x_node_displacement(doc, design, t):
    """{"extract": "nodeDisplacement", "combination", "node", "dof": ux|uy|uz|rx|ry|rz}."""
    _need(t, "node", "dof")
    return resolve_path(_combo(doc, t), f"/displacements/{t['node']}/{t['dof']}")


def _x_reaction(doc, design, t):
    """{"extract": "reaction", "combination", "node", "component": Fx|Fy|Fz|Mx|My|Mz}."""
    _need(t, "node", "component")
    return resolve_path(_combo(doc, t), f"/reactions/{t['node']}/{t['component']}")


def _x_sum_reactions(doc, design, t):
    """{"extract": "sumReactions", "combination", "component": Fx|Fy|Fz}: sum over the supports."""
    _need(t, "component")
    return float(sum(r[t["component"]] for r in _combo(doc, t)["reactions"].values()))


def _x_end_force(doc, design, t):
    """{"extract": "memberEndForce", "combination", "member", "end": "i" | "j", "component": Fx|...|Mz}:
    force by the joint ON the member, local axes."""
    _need(t, "member", "end", "component")
    return resolve_path(_combo(doc, t), f"/members/{t['member']}/endForces/{t['end']}/{t['component']}")


def _x_design_check(doc, design, t):
    """{"extract": "designCheck", "member": design member id, "check": name of the limit state (exact, or
    a unique part of it), "field": "ratio" | "capacity" | "nominal" | "demand", "combination": optional}:
    a number of a design check (the governing combination unless one is named)."""
    _need(t, "member", "check")
    if design is None:
        raise TargetError("no design results")
    if t["member"] not in design["members"]:
        raise TargetError(f"no design member {t['member']!r}")
    m = design["members"][t["member"]]
    checks = m["byCombination"].get(t["combination"], []) if "combination" in t else m["checks"]
    hit = [c for c in checks if c["name"] == t["check"]] or [c for c in checks if t["check"] in c["name"]]
    if len(hit) != 1:
        raise TargetError(f"check {t['check']!r}: {len(hit)} match(es) among " + "; ".join(c["name"] for c in checks))
    fld = t.get("field", "ratio")
    v = hit[0].get(fld)
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise TargetError(f"check {hit[0]['name']!r} has no number '{fld}' (status {hit[0].get('status')}: "
                          f"{hit[0].get('reason', '')})")
    return float(v)


EXTRACTORS = {"bucklingFactor": _x_buckling_factor, "memberStation": _x_member_station,
              "memberExtreme": _x_member_extreme, "nodeDisplacement": _x_node_displacement,
              "reaction": _x_reaction, "sumReactions": _x_sum_reactions, "memberEndForce": _x_end_force,
              "designCheck": _x_design_check}
_ALIASES = {"buckling": "bucklingFactor", "criticalfactor": "bucklingFactor", "criticalloadfactor": "bucklingFactor",
            "station": "memberStation", "stationvalue": "memberStation", "memberstationvalue": "memberStation",
            "membermax": "memberExtreme", "memberabsmax": "memberExtreme", "stationextreme": "memberExtreme",
            "displacement": "nodeDisplacement", "nodaldisplacement": "nodeDisplacement",
            "nodereaction": "reaction", "supportreaction": "reaction", "endforce": "memberEndForce",
            "memberendforces": "memberEndForce", "design": "designCheck"}
_EXTRACTOR_INDEX = {**{_norm(k): k for k in EXTRACTORS}, **{k: v for k, v in _ALIASES.items()}}


def _label(t: dict) -> str:
    if "path" in t:
        return ("design:" if t.get("document") == "design" else "") + t["path"]
    args = ", ".join(f"{k}={t[k]}" for k in t if k not in ("extract", "value", "unit", "tolerance", "precisionNote",
                                                            "sign", "note"))
    return f"{t['extract']}({args})"


def evaluate_targets(spec: ProblemSpec, results_document: dict, design_document: dict | None = None) -> list:
    """Compare every target with the engine's documents (file units).  Returns ``[TargetResult]``."""
    rows = []
    for k, t in enumerate(spec.targets):
        expected, tol = float(t["value"]), t["tolerance"]
        try:
            if "path" in t:
                doc = design_document if t.get("document") == "design" else results_document
                if doc is None:
                    raise TargetError("no design results")
                actual = resolve_path(doc, t["path"])
            else:
                actual = EXTRACTORS[_EXTRACTOR_INDEX[_norm(t["extract"])]](results_document, design_document, t)
        except TargetError as e:
            rows.append(TargetResult(k, _label(t), expected, None, None, None, tol, False, t.get("unit", ""), str(e)))
            continue
        except (KeyError, TypeError, ValueError) as e:
            rows.append(TargetResult(k, _label(t), expected, None, None, None, tol, False, t.get("unit", ""),
                                     f"cannot be evaluated: {type(e).__name__}: {e}"))
            continue
        diff = abs(actual - expected)
        rel = diff / abs(expected) if expected != 0.0 else None
        ok = diff <= tol["abs"] if "abs" in tol else (rel is not None and rel <= tol["rel"])
        rows.append(TargetResult(k, _label(t), expected, actual, diff, rel, tol, bool(ok), t.get("unit", "")))
    return rows


# ====================================================================================== verification
#: categories of a tier that is not passed which do NOT make an active problem fail (they are reported)
_NOT_A_FAILURE = ("not-representable", "not-covered", "not-run")


def verify_problem(directory, *, tiers="all", budget: float | None = None, cache=None,
                   foreign_timeout: float = 120.0, extra_evidence: bool = False) -> ProblemRecord:
    """Verify one problem of the library: its targets, and the tiers that apply to its model."""
    from .pipeline import verify_model
    t0 = time.perf_counter()
    d = Path(directory)
    try:
        spec = load_problem(d)
    except ProblemError as e:
        return ProblemRecord(d.name, str(d), "invalid", errors=e.errors, seconds=time.perf_counter() - t0)
    rec = ProblemRecord(spec.id, str(d), "failed", spec.status, spec.tier, spec.title, spec.source, spec.basis,
                        quarantine=spec.quarantine, notes=list(spec.notes), limits_note=spec.limits_note)
    model_record = verify_model(spec.model_path, tiers=tiers, budget=budget, cache=cache, targets=spec,
                                foreign_timeout=foreign_timeout, extra_evidence=extra_evidence,
                                name=f"{d.parent.name}/{spec.id}")
    rec.model_record = model_record
    failures = []
    if model_record.kind != "results" or model_record.status == "refused":
        why = "; ".join(p.get("message", "") for p in (model_record.refusal or {}).get("problems", [])[:2]) \
            or "the model is unstable (the engine returned a mechanism diagnosis)"
        failures.append(f"the engine produced no results for model.json: {why}")
    else:
        t = model_record.tiers[spec.tier]
        rec.targets = [TargetResult(**_target_kwargs(row)) for row in t.details.get("targets", [])]
        for r in rec.targets:
            if not r.passed:
                failures.append(f"target {r.label}: expected {r.expected:.10g}, engine "
                                + ("no value" if r.actual is None else f"{r.actual:.10g}")
                                + (f" ({r.error})" if r.error else
                                   f" (difference {r.difference:.3e}, tolerance {r.tolerance})"))
        for k, tr in model_record.tiers.items():
            if tr.applies and tr.passed is not True and tr.category not in _NOT_A_FAILURE:
                if k == spec.tier:                            # missed targets are listed above, one by one
                    auto = tr.details.get("failed") if k == "A" else None
                    if auto:
                        failures.append("tier A (ODE oracle against the engine): " + "; ".join(auto[:3]))
                    continue
                failures.append(f"tier {k}: {tr.summary}")
    rec.failures = failures
    rec.outcome = "quarantined" if spec.status == "quarantined" else ("passed" if not failures else "failed")
    rec.seconds = time.perf_counter() - t0
    return rec


def _target_kwargs(row: dict) -> dict:
    return {"index": row["index"], "label": row["target"], "expected": row["expected"], "actual": row["actual"],
            "difference": row["difference"], "relative": row["relative"], "tolerance": row["tolerance"],
            "passed": row["passed"], "unit": row.get("unit", ""), "error": row.get("error", "")}


def verify_library(root=None, *, tiers="all", budget: float | None = None, cache=None, only=None,
                   foreign_timeout: float = 120.0, progress=None) -> LibraryReport:
    """Verify every problem under ``root`` (default: ``verification/problems`` of the repository).

    ``budget``: seconds per problem.  ``only``: iterable of substrings; a problem is run when its
    "<source-slug>/<id>" contains one of them.  ``progress``: callable(ProblemRecord)."""
    from .pipeline import software_versions
    t0 = time.perf_counter()
    root_path = Path(root) if root else library_root()
    dirs = find_problems(root_path)
    if only:
        only = list(only)
        dirs = [d for d in dirs if any(o in f"{d.parent.name}/{d.name}" for o in only)]
    records = []
    for d in dirs:
        r = verify_problem(d, tiers=tiers, budget=budget, cache=cache, foreign_timeout=foreign_timeout)
        records.append(r)
        if progress:
            progress(r)
    return LibraryReport(str(root_path) if root_path else "", records, software_versions(), time.perf_counter() - t0)
