"""The written tolerance policy of the verification pipeline (docs/VERIFICATION.md is the prose).

One set of numbers, applied uniformly.  Nothing in the pipeline changes a tolerance per model or per
problem; a case that misses is quarantined with a diagnosis.
"""
from __future__ import annotations

import math

__all__ = ["POLICY", "printed_precision_tolerance", "check_target_tolerance", "REL_FLOOR"]

REL_FLOOR = 1.0e-6

POLICY = {
    # tier A: engine against the ODE oracle (accurate to about 1e-8) / derived-here targets
    "A": {"relative": 1.0e-6,
          "scale": "largest value of the same kind (translation, rotation, force, moment) in the member line"},
    # tier B: published targets; tolerance from the printed precision or a tolerance stated by the source
    "B": {"relFloor": REL_FLOOR, "rule": "half a unit in the last printed digit, as a relative tolerance "
                                          "(floor 1e-6), or the acceptance tolerance stated by the source"},
    # tier C: benchfea.crosscheck defaults; the estimated error of a converged approximation widens the
    # tolerance of the comparisons it takes part in
    "C": {"relative": 1.0e-6, "maxApproximation": 5.0e-3, "requiredSolvers": 2, "estimateSafety": 2.0,
          "scale": "model-wide scales of benchfea.crosscheck.compare"},
    # tier D: the tolerances stored with every invariant by the solver, plus
    "D": {"superposition": 1.0e-9, "reciprocity": 1.0e-9},
    # tier E
    "E": {"limit": 1.0e-5,                  # extrapolated limit of a stiffness series against the engine
          "splitFirstOrder": 1.0e-8,        # member split, linear analysis
          "splitNonlinear": 1.0e-6,         # member split, iterated analysis (the iteration tolerance governs)
          "designVariables": 0.0,           # bitwise
          "seriesRatio": 10.0,              # stiffness ratio between successive members of a series
          "observedOrder": (0.5, 1.5)},     # accepted range of the observed order of convergence
    # tier F: no numeric tolerance (coverage registry)
    "F": {},
    # verification of a mechanism diagnosis
    "diagnosis": {"nullResidual": 1.0e-8, "singular": 1.0e-9},
}


def printed_precision_tolerance(printed: str) -> dict:
    """Tolerance implied by a number AS PRINTED in a source: half a unit in the last printed digit,
    expressed as a relative tolerance with the floor ``REL_FLOOR``.

    ``"0.3125"`` -> half unit 0.00005 -> rel 1.6e-4;  ``"1250"`` -> half unit 0.5 -> rel 4e-4;
    ``"1.20e3"`` -> half unit 5 -> rel 4.2e-3.  Trailing zeros of an integer without decimal point are
    taken as significant (the conservative reading: it gives the tighter tolerance).
    """
    s = printed.strip().lower().replace(",", "")
    mant, _, exp = s.partition("e")
    e = int(exp) if exp else 0
    digits = mant.lstrip("+-")
    decimals = len(digits.split(".")[1]) if "." in digits else 0
    value = abs(float(s))
    half = 0.5 * 10.0 ** (e - decimals)
    if value == 0.0:
        return {"abs": half}
    return {"rel": max(half / value, REL_FLOOR)}


def check_target_tolerance(tol: dict) -> str | None:
    """A target tolerance is ``{"rel": x}`` or ``{"abs": y}`` with a positive finite number; returns
    the problem as text or None."""
    if not isinstance(tol, dict) or len(tol) != 1 or next(iter(tol)) not in ("rel", "abs"):
        return 'tolerance must be {"rel": x} or {"abs": y}'
    v = next(iter(tol.values()))
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0.0:
        return "tolerance value must be a positive finite number"
    if "rel" in tol and v < REL_FLOOR:
        return f"relative tolerance below the policy floor {REL_FLOOR:g}"
    return None
