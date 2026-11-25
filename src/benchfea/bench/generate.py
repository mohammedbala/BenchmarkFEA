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
   requested level -- or would give it only for one side of the amplification band -- is dropped
   before any analysis;
4. the SIZING LOOP analyses the model and runs the AISC 360-22 member checks, then moves each
   sizing group up or down the generic section catalogue until its largest design ratio is inside
   the target band (``RATIO_LOW`` .. ``RATIO_HIGH``) and the buckling factor is above
   ``BUCKLING_MIN`` (stability check of the gravity combination, see ``_size``);
5. the model is accepted only if the analysis ran with all invariants passed, the design status
   is ``pass`` or ``conditional`` (no failed check and no applicable limit state left unchecked), no
   design ratio is above ``RATIO_MAX``, at least ``RATIO_BAND_FRACTION`` of the design members have their governing
   ratio inside ``RATIO_BAND``, the stability check passed, the second-order amplification is at
   most ``AMPLIFICATION_MAX`` and the MEASURED difficulty level equals the requested one; otherwise
   attempt k + 1 (at most ``max_attempts``).

The acceptance rules are the same for every family and are never relaxed for a model; a family /
level that cannot be produced raises :class:`GenerationFailed` with the reasons of every attempt.
"""
from __future__ import annotations

import copy
import os
import sys
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
           "RATIO_LOW", "RATIO_HIGH", "RATIO_MAX", "BUCKLING_MIN", "AMPLIFICATION_MAX", "RATIO_BAND",
           "RATIO_BAND_FRACTION", "GENERATOR_VERSION"]

GENERATOR_VERSION = "1"        # bump when a change alters the output of any (family, seed, difficulty)
RATIO_LOW = 0.50               # a group whose largest ratio is below this is tried one size smaller
RATIO_HIGH = 0.95              # a group whose largest ratio is above this is made one size larger
RATIO_MAX = 1.00               # no design ratio of an accepted model exceeds this
BUCKLING_MIN = 2.0             # smallest accepted elastic buckling factor of the checked combination
AMPLIFICATION_MAX = 1.5        # largest accepted second-order amplification of the largest translation
RATIO_BAND = (0.3, 1.1)        # "sensible" governing design ratio of a design member
RATIO_BAND_FRACTION = 0.5      # at least this fraction of the design members is inside RATIO_BAND
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
    approximate | not-eligible | unsupported`` and ``crossCheckable`` (both solvers eligible);
    ``buckling``: ``none`` (not asked for) | ``unsupported`` (asked for; no foreign solver has it).
    Overall ``crossCheckable``: every combination is, and the file asks for no buckling analysis."""
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
    buckling = "none"
    if flags["buckling"]:
        buckling = "unsupported"
        notes.append("linear buckling is not reproduced by either foreign solver: the buckling results of this "
                     "file cannot be cross-checked")
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
            if flags["rigidLink"] or flags["rigidDiaphragm"]:
                worse(d, "pynite", "not-eligible", "second order with rigid links / diaphragms: PyNite's stiff-member "
                                                   "emulation does not reproduce the second-order result (measured)")
        if c.get("notionalLoad") is not None:
            for s in d:
                worse(d, s, "unsupported", "notional loads must be supplied to the translators by the caller")
        if c.get("stiffnessReduction"):
            notes.append("stiffness reduction: a computed tau_b must be supplied to the translators by the caller")
        d["crossCheckable"] = all(d[s] in ("exact", "approximate") for s in ("opensees", "pynite"))
        by[cid] = d
    worst = {s: max((by[c][s] for c in by), key=rank.index) for s in ("opensees", "pynite")}
    return {**worst, "buckling": buckling,
            "crossCheckable": all(v["crossCheckable"] for v in by.values()) and buckling == "none",
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


def _stable_at(doc: dict, combination: str, factor: float) -> bool:
    """Is the model below its elastic critical load at ``factor`` times the loads of ``combination``?

    One second-order analysis of that combination with all load factors multiplied by ``factor``:
    the engine refuses it (``stability.beyond-critical-load``) exactly when the tangent stiffness of
    the loaded structure is no longer positive definite.  Used for models that do not ask for a
    buckling analysis in the file (the engine's buckling analysis is far slower)."""
    probe = copy.deepcopy(doc)
    c = probe["combinations"][combination]
    probe["combinations"] = {combination: {"factors": {k: v * factor for k, v in c["factors"].items()},
                                           "analysisType": "secondOrder"}}
    probe.get("analysis", {}).pop("buckling", None)
    probe.pop("design", None)
    used = set(c["factors"])
    probe["loadCases"] = {k: v for k, v in probe["loadCases"].items() if k in used}
    _, results, codes = _analyse(probe)
    if results is not None:
        return bool(results.combinations[combination].converged)
    if {"stability.beyond-critical-load", "solve.not-converged"} & set(codes):
        return False
    raise GenerationRejected("stability probe refused: " + ", ".join(codes))


def _size(builder):
    """The sizing loop -> (document, model, results, design, first-order results, record).

    Phase 1 sizes on the FAST model (every combination first order, no buckling analysis): each
    group moves along the catalogue until its largest design ratio is in the band.  Phase 2
    analyses the complete document; groups whose ratio now exceeds ``RATIO_HIGH`` (second-order
    effects) go up, and a combination beyond its critical load or a failed stability check sends
    the column / brace groups up (every group when the model has no column / brace group); phase 2
    never goes down.  Stability check of ``builder.stability_combination``: the buckling factor of
    the file's own buckling analysis if it asks for one, else :func:`_stable_at` ``BUCKLING_MIN``.
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
    stiffer = ([g for g, v in groups.items() if v["kind"] != "beam" and not v["fixed"]]
               or [g for g, v in groups.items() if not v["fixed"]])
    combo = builder.stability_combination
    stability = None
    fresh = True                                                # ``first`` belongs to the current sections
    for _ in range(MAX_SIZING_ITERATIONS):
        step = {}
        if builder.nonlinear or not fresh:
            doc = builder.render()
            model, results, codes = _analyse(doc)
            analyses += 1
            if results is None:
                if not ({"stability.beyond-critical-load", "solve.not-converged"} & set(codes)):
                    raise GenerationRejected("analysis refused: " + ", ".join(codes))
                step = {g: 1 for g in stiffer}                  # too slender
            else:
                design = design_module.check(model, results)
                step = _steps(builder, model, design, locked, False)
        if not step and combo is not None:
            if combo in results.buckling and len(results.buckling[combo].factors):
                factor = float(results.buckling[combo].factors[0])
                stability = {"combination": combo, "basis": "linear buckling analysis of the file",
                             "bucklingFactor": round(factor, 3), "minimum": BUCKLING_MIN}
                ok = factor >= BUCKLING_MIN
            else:
                analyses += 1
                ok = _stable_at(doc, combo, BUCKLING_MIN)
                stability = {"combination": combo, "bucklingFactor": None, "minimum": BUCKLING_MIN,
                             "basis": f"second-order analysis at {BUCKLING_MIN:g} x the combination is below the "
                                      "elastic critical load (tangent stiffness positive definite)"}
            if not ok:
                step = {g: 1 for g in stiffer}
        if not step:
            break
        first = None                                            # sections changed: no longer the same model
        fresh = False
        _apply(builder, step, history, locked)
    else:
        raise GenerationRejected("sizing did not settle")
    record = {"iterations": len(history) + 1, "analyses": analyses, "sections": builder.section_names(),
              "stability": stability}
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
        drawn_units, drawn_vertical = draw.choice(sorted(UNIT_SYSTEMS)), draw.choice(("Y", "Z"))
        # both draws are always consumed: forcing the units or the axis does not change the structure
        unit_system, vertical = units or drawn_units, vertical_axis or drawn_vertical
        try:
            builder = fam.build(draw, difficulty, unit_system, vertical, features)
            builder.name = f"{family} seed {seed} level {difficulty}"
            pre = measure_difficulty(builder.render(), analyse=False)
            if pre["levelRange"] != [difficulty, difficulty]:
                # also dropped: a draw whose level would depend on which side of the amplification
                # band the analysis lands (the level must not be a gamble on an expensive analysis)
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
            if design.status not in ("pass", "conditional"):
                # "incomplete": a limit state that applies is outside what the design module checks
                # (e.g. torsion above its screening threshold) -- no complete design reference
                blocked = sorted({n.split(":")[0] for r in design.members.values() for n in r.not_checked})
                raise GenerationRejected(f"design status {design.status}" + (f" ({', '.join(blocked)})" if blocked else ""))
            ratios = sorted(member_ratios.values())
            in_band = sum(1 for r in ratios if RATIO_BAND[0] <= r <= RATIO_BAND[1])
            if ratios and in_band < RATIO_BAND_FRACTION * len(ratios):
                raise GenerationRejected(f"only {in_band} of {len(ratios)} design members have a governing ratio "
                                         f"in {RATIO_BAND}")
            factors = {cid: float(b.factors[0]) for cid, b in results.buckling.items() if len(b.factors)}
            if factors and min(factors.values()) < BUCKLING_MIN:
                raise GenerationRejected(f"buckling factor {min(factors.values()):.2f} below {BUCKLING_MIN}")
            if measured["amplification"] is not None and measured["amplification"] > AMPLIFICATION_MAX:
                raise GenerationRejected(f"second-order amplification {measured['amplification']:.2f} above "
                                         f"{AMPLIFICATION_MAX}")
        except GenerationRejected as exc:
            reasons.append(f"attempt {attempt}: {exc}")
            if os.environ.get("BENCHFEA_GENERATE_TRACE"):
                print(f"[generate {family} seed {seed} level {difficulty}] {reasons[-1]}", file=sys.stderr, flush=True)
            continue
        flags = feature_flags(doc, model)
        metadata = {
            "generator": {"version": GENERATOR_VERSION, "engine": __version__},
            "family": family, "seed": seed, "difficulty": difficulty, "attempt": attempt,
            "rejected": list(reasons),
            "units": unit_system, "verticalAxis": vertical,
            "parameters": builder.parameters,
            "requestedFeatures": features,
            "features": flags,
            "measured": measured,
            "sizing": {k: v for k, v in sizing.items() if k != "stability"},
            "stability": sizing["stability"],
            "design": {"status": design.status, "nDesignMembers": len(design.members),
                       "ratioMin": round(float(ratios[0]), 3) if ratios else None,
                       "ratioMax": round(float(ratios[-1]), 3) if ratios else None,
                       "ratioMedian": round(float(ratios[len(ratios) // 2]), 3) if ratios else None,
                       "fractionBetween0.3And1.1": round(in_band / len(ratios), 3) if ratios else None,
                       "notChecked": sorted({n.split(":")[0] for r in design.members.values() for n in r.not_checked})},
            "bucklingFactors": {k: round(v, 3) for k, v in factors.items()},
            "invariantsPassed": all(c.invariants.passed for c in results.combinations.values()),
            "tierC_expected": tier_c_expected(doc, flags),
        }
        return Generated(doc, metadata)
    raise GenerationFailed(f"no acceptable {family!r} model for seed {seed}, difficulty {difficulty} "
                           f"in {max_attempts} attempts", reasons)
