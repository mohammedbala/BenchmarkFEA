"""Train / held-out split of the benchmark (docs/BENCHMARK.md section 7).

The split is defined by the committed file ``splits.json`` next to this module and by one rule:

    is_held_out(family, curve_ball, seed)  ==  (family, curve_ball) is a held-out COMBINATION
                                               or  seed % modulus is a held-out RESIDUE

* held-out combinations: pairs (family, curve ball) that never appear in the train split, whatever
  the seed -- an agent trained on the train split has not seen that curve ball on that family;
* held-out seeds: for every other pair, the seeds with ``seed % 5 == 4`` are held out.

Which pairs are held-out combinations is not a choice made pair by pair: ``select_held_out`` applies
a rule that looks at nothing but names (docs/BENCHMARK.md section 7).  Per family, among the curve
balls that the applicability table (docs/CURVEBALLS.md section 5, copied into ``splits.json`` as
``heldOutSelection.applicable``) lists with an effect class other than ``invalid-model`` alone, the
``perFamily`` = 2 with the smallest SHA-256 of the text ``benchfea.bench|held-out-pair|<family>|<curve
ball>`` are held out.  Whether a pair ends up with verified instances plays no part.

The split also decides the base models: ``tasks.derived_seed`` puts the split into the generator seed
(the lowest bit is 1 for held-out, 0 for train), so the two splits never generate from the same seed.

``seed`` is the seed of the task instance as requested (the id of the instance); the derivation
chain that replaces a quarantined instance (``tasks.derived_seed``) does not change the split: the
split is a function of ``(family, curve_ball, seed)`` alone and ``is_held_out`` is its only definition.

    from benchfea.bench.splits import is_held_out, split_of, combinations
    is_held_out("portal_frame", "release_added", 9)      # True  (9 % 5 == 4)
    split_of("portal_frame", "release_added", 3)         # "train"
    combinations("train")                                # [(family, curve ball), ...] served by default

``servedCombinations`` in the file is the list of pairs for which the dataset build produced at least
one verified instance (measured, see data/tasks/manifest.json); the environment draws from it when
the caller does not name a family and a curve ball.  A pair outside the list can still be requested
explicitly; it then usually has no verified instance and the request is refused with the reason.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path

__all__ = ["is_held_out", "split_of", "combinations", "load_splits", "SPLITS", "SPLITS_FILE", "seeds",
           "held_out_rank", "select_held_out"]

SPLITS = ("train", "held-out")
SPLITS_FILE = Path(__file__).with_name("splits.json")


@lru_cache(maxsize=1)
def load_splits() -> dict:
    """The committed split definition (``splits.json``)."""
    d = json.loads(SPLITS_FILE.read_text(encoding="utf-8"))
    d["_heldOut"] = frozenset((f, c) for f, c in d["heldOutCombinations"])
    return d


def held_out_rank(family: str, curve_ball: str) -> str:
    """Sort key of the held-out selection: SHA-256 (hex) of ``benchfea.bench|held-out-pair|family|curve_ball``."""
    return hashlib.sha256(f"benchfea.bench|held-out-pair|{family}|{curve_ball}".encode("utf-8")).hexdigest()


def select_held_out(applicable: dict, per_family: int = 2) -> list:
    """The held-out combinations given ``applicable`` = {family: [curve balls eligible for that family]}:
    per family (in the order of ``applicable``) the ``per_family`` curve balls of smallest
    :func:`held_out_rank`.  ``splits.json`` must equal this for its own ``heldOutSelection`` (tested)."""
    out = []
    for family, names in applicable.items():
        ranked = sorted(set(names), key=lambda c: held_out_rank(family, c))
        out += [[family, c] for c in ranked[:per_family]]
    return out


def is_held_out(family: str, curve_ball: str, seed: int) -> bool:
    """Is the instance ``(family, curve_ball, seed)`` part of the held-out split?"""
    d = load_splits()
    if (family, curve_ball) in d["_heldOut"]:
        return True
    rule = d["seedRule"]
    return int(seed) % int(rule["modulus"]) in rule["heldOutResidues"]


def split_of(family: str, curve_ball: str, seed: int) -> str:
    return "held-out" if is_held_out(family, curve_ball, seed) else "train"


def combinations(split: str | None = None) -> list:
    """Served (family, curve ball) pairs.  ``split="train"``: the pairs that have train instances
    (the held-out combinations are left out); ``"held-out"``: every served pair (a held-out
    combination with any seed, any other pair with a held-out seed); ``None``: every served pair."""
    if split not in (None,) + SPLITS:
        raise ValueError(f"split must be one of {SPLITS}, not {split!r}")
    d = load_splits()
    pairs = [(f, c) for f, c in d["servedCombinations"]]
    if split == "train":
        pairs = [p for p in pairs if p not in d["_heldOut"]]
    return pairs


def seeds(split: str, family: str, curve_ball: str, count: int, start: int = 0) -> list:
    """The first ``count`` seeds >= ``start`` of the pair that belong to ``split``."""
    if split not in SPLITS:
        raise ValueError(f"split must be one of {SPLITS}, not {split!r}")
    out, s = [], int(start)
    held_combo = (family, curve_ball) in load_splits()["_heldOut"]
    if split == "train" and held_combo:
        return out
    while len(out) < count:
        if is_held_out(family, curve_ball, s) == (split == "held-out"):
            out.append(s)
        s += 1
    return out
