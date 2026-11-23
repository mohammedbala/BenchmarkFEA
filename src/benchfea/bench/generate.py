"""Seeded parametric model generators: ``generate(family, seed, difficulty) -> Generated``.

    from benchfea.bench import generate
    g = generate("portal_frame", seed=7, difficulty=2)
    g.document          # input document (dict) -- valid, stable, sized
    g.metadata          # parameters drawn, feature flags, measured difficulty, sizing, tier-C expectation
    g.json()            # canonical JSON text (byte-identical for the same family, seed, difficulty)

How a model is made (deterministic for a given ``(family, seed, difficulty)`` and options):

1. attempt k draws from ``Draw(stable_seed(family, seed, difficulty, k))`` (SHA-256 -> PCG64);
2. the family builds geometry, loads, combinations and design data in abstract axes with rounded
   engineering values, and names its sizing groups (beams, columns, braces, ...);
3. the measurable difficulty that needs no analysis is checked; a draw that cannot give the
   requested level is dropped before any analysis;
4. the SIZING LOOP analyses the model and runs the AISC 360-22 member checks, then moves each
   sizing group up or down the generic section catalogue until its largest design ratio is inside
   the target band (``RATIO_LOW`` .. ``RATIO_HIGH``) and the buckling factor is above
   ``BUCKLING_MIN``;
5. the model is accepted only if the analysis ran with all invariants passed, no design ratio is
   above ``RATIO_MAX``, the buckling factor is above ``BUCKLING_MIN`` and the MEASURED difficulty
   level equals the requested one; otherwise attempt k + 1 (at most ``max_attempts``).

The acceptance rules are the same for every family and are never relaxed for a model; a family /
level that cannot be produced raises :class:`GenerationFailed` with the reasons of every attempt.
"""
from __future__ import annotations

from typing import NamedTuple

from .. import __version__
from ..io import dumps
from ..model import build_model
from ..problems import ModelInvalid
from .difficulty import difficulty as measure_difficulty
from .difficulty import feature_flags
from .families import FAMILIES
from .families._common import UNIT_SYSTEMS, Draw, GenerationRejected, catalogue_size, stable_seed

__all__ = ["generate", "Generated", "GenerationFailed", "families", "levels", "tier_c_expected",
           "RATIO_LOW", "RATIO_HIGH", "RATIO_MAX", "BUCKLING_MIN", "GENERATOR_VERSION"]

GENERATOR_VERSION = "1"        # bump when a change alters the output of any (family, seed, difficulty)
RATIO_LOW = 0.50               # a group whose largest ratio is below this is tried one size smaller
RATIO_HIGH = 0.95              # a group whose largest ratio is above this is made one size larger
RATIO_MAX = 1.00               # no design ratio of an accepted model exceeds this
BUCKLING_MIN = 2.0             # smallest accepted elastic buckling factor of the checked combination
MAX_SIZING_ITERATIONS = 40


class GenerationFailed(RuntimeError):
    """No acceptable model within ``max_attempts``; ``.reasons`` lists why each attempt was dropped."""

    def __init__(self, message: str, reasons: list):
        super().__init__(message + "\n  " + "\n  ".join(reasons))
        self.reasons = reasons


class Generated(NamedTuple):
    document: dict
    metadata: dict

    def json(self) -> str:
        return dumps(self.document, "json")

    def yaml(self) -> str:
        return dumps(self.document, "yaml")


def families() -> list:
    """Names of the model families."""
    return list(FAMILIES)


def levels(family: str) -> tuple:
    """Difficulty levels the family can produce."""
    return tuple(FAMILIES[family].LEVELS)


# ----------------------------------------------------------------------------------------- tier C
def tier_c_expected(document: dict, flags: dict | None = None) -> dict:
    """What the two foreign solvers of ``benchfea.crosscheck`` are expected to make of the model,
    from the capability table (docs/CROSSCHECK.md section 4) -- an expectation derived from the
    features present, not a result.  Per combination: ``opensees`` / ``pynite`` in ``exact |
    approximate | not-eligible | unsupported`` and ``crossCheckable`` (both solvers eligible)."""
    flags = flags or feature_flags(document)
    notes = []
    base = {"opensees": "exact", "pynite": "exact"}
    rank = ["exact", "approximate", "not-eligible", "unsupported"]

    def worse(d, solver, level, note):
        if rank.index(level) > rank.index(d[solver]):
            d[solver] = level
        if note not in notes:
            notes.append(note)

    if flags["shearDeformation"]:
        worse(base, "pynite", "unsupported", "shear deformation: PyNite has none (only OpenSees eligible, first order)")
    if flags["tensionOnly"] or flags["compressionOnly"]:
        for s in base:
            worse(base, s, "unsupported", "tension-/compression-only members are not reproduced by either "
                                          "translator (a check is possible only on the frozen active set)")
    if flags["releasesPartialOther"]:
        worse(base, "pynite", "unsupported", "axial / shear / torsion connection springs: PyNite releases are boolean")
    if flags["releasesFullNonBending"]:
        worse(base, "opensees", "approximate", "full axial / shear / torsion releases: OpenSees exact only for "
                                               "members along global axes without offsets, else penalty springs")
    if flags["rigidOffsets"]:
        worse(base, "pynite", "approximate", "rigid offsets: stiff stub members in PyNite")
    if flags["rigidLink"] or flags["rigidDiaphragm"]:
        worse(base, "pynite", "approximate", "rigid links / diaphragms: stiff members in PyNite")
    if flags["buckling"]:
        notes.append("linear buckling is not reproduced by either foreign solver")
    analysis_type = document.get("analysis", {}).get("type", "firstOrder")
    by = {}
    for cid, c in document["combinations"].items():
        d = dict(base)
        if c.get("analysisType", analysis_type) == "secondOrder":
            for s in d:
                worse(d, s, "approximate", "second order: both foreign solvers are converged approximations")
            if flags["shearDeformation"]:
                worse(d, "opensees", "unsupported", "second order with shear deformation: no foreign solver")
            if flags["rigidOffsets"]:
                for s in d:
                    worse(d, s, "not-eligible", "second order with rigid offsets: approximation without an "
                                                "error estimate (not eligible)")
        if c.get("notionalLoad") is not None:
            for s in d:
                worse(d, s, "unsupported", "notional loads must be supplied to the translators by the caller")
        if c.get("stiffnessReduction"):
            notes.append("stiffness reduction: a computed tau_b must be supplied to the translators by the caller")
        d["crossCheckable"] = all(d[s] in ("exact", "approximate") for s in ("opensees", "pynite"))
        by[cid] = d
    worst = {s: max((by[c][s] for c in by), key=rank.index) for s in ("opensees", "pynite")}
    return {**worst, "crossCheckable": all(v["crossCheckable"] for v in by.values()),
            "byCombination": by, "notes": notes}


# ----------------------------------------------------------------------------------------- sizing
def _group_ratios(builder, model, design):
    """Largest non-advisory design ratio per sizing group and per design member."""
    by_group, by_member = {}, {}
    for dmid, r in design.members.items():
        ratios = [c.ratio for c in r.checks if c.status == "ok" and not c.advisory and c.ratio is not None]
        if not ratios:
            continue
        ratio = max(ratios)
        by_member[dmid] = ratio
        g = builder.member_group[model.design_members[dmid].members[0]]
        by_group[g] = max(by_group.get(g, 0.0), ratio)
    return by_group, by_member


def _steps(builder, model, design, locked, allow_down: bool) -> dict:
    """Catalogue steps wanted by the design ratios of each sizing group."""
    ratios, _ = _group_ratios(builder, model, design)
    step = {}
    for g, v in builder.groups.items():
        r = ratios.get(g)
        if r is None or v["fixed"]:
            continue
        if r > RATIO_HIGH:
            step[g] = 1 if r < 1.6 else 2 if r < 3.0 else 3
        elif allow_down and r < RATIO_LOW and v["index"] > 0 and g not in locked:
            step[g] = -1 if (r > 0.2 or v["index"] < 2) else -2
    return step


def _apply(builder, step: dict, history: list, locked: set) -> None:
    for g, s in step.items():
        v = builder.groups[g]
        new = v["index"] + s
        if new >= catalogue_size(v["kind"], builder.unit_system):
            raise GenerationRejected(f"group {g!r} needs a section beyond the catalogue")
        if s > 0 and any(h.get(g, 0) < 0 for h in history):
            locked.add(g)                                       # went down and came back up: stay
        v["index"] = max(new, 0)
    history.append(step)


def _analyse(doc):
    """-> (model, results or None, refusal codes)."""
    from ..solve import analyze
    model = build_model(doc)
    try:
        return model, analyze(model), []
    except ModelInvalid as exc:
        return model, None, sorted({p.code for p in exc.problems})


def _size(builder):
    """The sizing loop -> (document, model, results, design, first-order results, record).

    Phase 1 sizes on the FAST model (every combination first order, no buckling analysis): each
    group moves along the catalogue until its largest design ratio is in the band.  Phase 2
    analyses the complete document; groups whose ratio now exceeds ``RATIO_HIGH`` (second-order
    effects) go up, a combination beyond its critical load sends every group up, and a buckling
    factor below ``BUCKLING_MIN`` sends the column / brace groups up; phase 2 never goes down.
    """
    from .. import design as design_module

    groups = builder.groups
    locked, history = set(), []
    analyses = 0
    for _ in range(MAX_SIZING_ITERATIONS):
        doc = builder.render(fast=True)
        model, first, codes = _analyse(doc)
        analyses += 1
        if first is None:
            raise GenerationRejected("analysis refused: " + ", ".join(codes))
        design = design_module.check(model, first)
        step = _steps(builder, model, design, locked, True)
        if not step:
            break
        _apply(builder, step, history, locked)
    else:
        raise GenerationRejected("sizing did not settle")
    results = first
    if builder.nonlinear:
        for _ in range(MAX_SIZING_ITERATIONS):
            doc = builder.render()
            model, results, codes = _analyse(doc)
            analyses += 1
            if results is None:
                if not ({"stability.beyond-critical-load", "solve.not-converged"} & set(codes)):
                    raise GenerationRejected("analysis refused: " + ", ".join(codes))
                step = {g: 1 for g, v in groups.items() if not v["fixed"]}      # too slender
            else:
                design = design_module.check(model, results)
                step = _steps(builder, model, design, locked, False)
                factors = [float(b.factors[0]) for b in results.buckling.values() if len(b.factors)]
                if factors and min(factors) < BUCKLING_MIN:
                    up = [g for g, v in groups.items() if v["kind"] != "beam" and not v["fixed"]] or list(groups)
                    step.update({g: max(1, step.get(g, 0)) for g in up})
            if not step:
                break
            first = None                                        # sections changed: no longer the same model
            _apply(builder, step, history, locked)
        else:
            raise GenerationRejected("sizing did not settle")
    record = {"iterations": len(history) + 1, "analyses": analyses, "sections": builder.section_names()}
    return doc, model, results, design, first, record


# ----------------------------------------------------------------------------------------- generate
def generate(family: str, seed: int, difficulty: int, *, units: str | None = None,
             vertical_axis: str | None = None, features: dict | None = None,
             max_attempts: int = 60) -> Generated:
    """One model of ``family`` for ``(seed, difficulty)``.

    ``units`` ("kN-m" | "kN-mm" | "kip-in") and ``vertical_axis`` ("Y" | "Z") are drawn from the seed
    unless given.  ``features``: {flag: bool} forcing a feature of the family on or off (the flags of
    each family are ``FAMILIES[family].FEATURES``); a forced feature can make a level unreachable.
    Raises ``ValueError`` for an unknown family / level / flag and :class:`GenerationFailed` when no
    acceptable model is found.
    """
    if family not in FAMILIES:
        raise ValueError(f"unknown family {family!r}; families: {', '.join(FAMILIES)}")
    fam = FAMILIES[family]
    if difficulty not in fam.LEVELS:
        raise ValueError(f"family {family!r} produces difficulty levels {tuple(fam.LEVELS)}, not {difficulty!r}")
    if units is not None and units not in UNIT_SYSTEMS:
        raise ValueError(f"unknown unit system {units!r}; use one of {', '.join(UNIT_SYSTEMS)}")
    if vertical_axis not in (None, "Y", "Z"):
        raise ValueError("vertical_axis must be 'Y' or 'Z'")
    features = dict(features or {})
    unknown = sorted(set(features) - set(fam.FEATURES))
    if unknown:
        raise ValueError(f"family {family!r} has no feature flag {unknown}; flags: {', '.join(fam.FEATURES)}")
    seed = int(seed)
    reasons = []
    for attempt in range(max_attempts):
        draw = Draw(stable_seed(GENERATOR_VERSION, family, seed, difficulty, attempt))
        unit_system = units or draw.choice(sorted(UNIT_SYSTEMS))
        vertical = vertical_axis or draw.choice(("Y", "Z"))
        try:
            builder = fam.build(draw, difficulty, unit_system, vertical, features)
            builder.name = f"{family} seed {seed} level {difficulty}"
            pre = measure_difficulty(builder.render(), analyse=False)
            if not pre["levelRange"][0] <= difficulty <= pre["levelRange"][1]:
                raise GenerationRejected(f"draw gives level {pre['levelRange']} (score {pre['score']})")
            doc, model, results, design, first, sizing = _size(builder)
            measured = measure_difficulty(doc, results, first_order=first)
            if measured["level"] != difficulty:
                raise GenerationRejected(f"measured level {measured['level']} (score {measured['score']}, "
                                         f"amplification {measured['amplification']})")
            _, member_ratios = _group_ratios(builder, model, design)
            worst = max(member_ratios.values(), default=0.0)
            if worst > RATIO_MAX:
                raise GenerationRejected(f"design ratio {worst:.2f} above {RATIO_MAX}")
            if design.status == "fail":
                raise GenerationRejected("design status fail")
            factors = {cid: float(b.factors[0]) for cid, b in results.buckling.items() if len(b.factors)}
            if factors and min(factors.values()) < BUCKLING_MIN:
                raise GenerationRejected(f"buckling factor {min(factors.values()):.2f} below {BUCKLING_MIN}")
        except GenerationRejected as exc:
            reasons.append(f"attempt {attempt}: {exc}")
            continue
        flags = feature_flags(doc, model)
        ratios = sorted(member_ratios.values())
        in_band = sum(1 for r in ratios if 0.3 <= r <= 1.1)
        metadata = {
            "generator": {"version": GENERATOR_VERSION, "engine": __version__},
            "family": family, "seed": seed, "difficulty": difficulty, "attempt": attempt,
            "rejected": list(reasons),
            "units": unit_system, "verticalAxis": vertical,
            "parameters": builder.parameters,
            "requestedFeatures": features,
            "features": flags,
            "measured": measured,
            "sizing": sizing,
            "design": {"status": design.status, "nDesignMembers": len(design.members),
                       "ratioMin": round(ratios[0], 3) if ratios else None,
                       "ratioMax": round(ratios[-1], 3) if ratios else None,
                       "ratioMedian": round(ratios[len(ratios) // 2], 3) if ratios else None,
                       "fractionBetween0.3And1.1": round(in_band / len(ratios), 3) if ratios else None,
                       "notChecked": sorted({n.split(":")[0] for r in design.members.values() for n in r.not_checked})},
            "bucklingFactors": {k: round(v, 3) for k, v in factors.items()},
            "invariantsPassed": all(c.invariants.passed for c in results.combinations.values()),
            "tierC_expected": tier_c_expected(doc, flags),
        }
        return Generated(doc, metadata)
    raise GenerationFailed(f"no acceptable {family!r} model for seed {seed}, difficulty {difficulty} "
                           f"in {max_attempts} attempts", reasons)
