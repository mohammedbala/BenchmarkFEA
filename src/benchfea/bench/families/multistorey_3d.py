"""Multi-bay multi-storey space frames (true 3D): moment frames in both directions.

Rigid beam-column joints along x and y, fixed / pinned / partially fixed bases, one-way floors on
the beams along x, cladding on the end frames, wind along x and along y as storey forces on the
windward face.  Optionally second order and pinned secondary beams (the beams along y on interior
lines).  No rigid diaphragm (family ``diaphragm_frame`` has one).
"""
from __future__ import annotations

from ._common import BUCKLING_CHANCE, BUCKLING_FLAG, Builder, flag
from ._space import space_frame

NAME = "multistorey_3d"
LEVELS = (3, 4, 5)
FEATURES = {
    "buckling": BUCKLING_FLAG,
    "secondOrder": "strength combinations analysed second order",
    "pinnedSecondary": "beams along y on interior lines pinned at both ends (needs two bays along x)",
    "baseSprings": "partially fixed bases (rotational springs about both horizontal axes)",
}
_SHAPES = {3: ((1, 1, 1), (1, 1, 2), (2, 1, 1)), 4: ((2, 1, 1), (1, 1, 2), (2, 1, 2), (2, 2, 1)),
           5: ((2, 1, 2), (2, 2, 1), (2, 2, 2), (2, 1, 3))}


def build(draw, level: int, unit_system: str, vertical: str, features: dict) -> Builder:
    B = Builder(unit_system, vertical, NAME)
    shape = draw.choice(_SHAPES[level])
    second = flag(features, "secondOrder", draw.chance({3: 0.2, 4: 0.7, 5: 0.9}[level]))
    pinned = flag(features, "pinnedSecondary", draw.chance({3: 0.2, 4: 0.5, 5: 0.9}[level])) and shape[0] >= 2
    springs = flag(features, "baseSprings", draw.chance({3: 0.1, 4: 0.4, 5: 0.8}[level]))
    n_combos = {3: draw.integer(3, 5), 4: draw.integer(4, 7), 5: draw.integer(6, 8)}[level]
    buckling = flag(features, "buckling", draw.chance(BUCKLING_CHANCE[level]))
    B.parameters = space_frame(B, draw, shape, second=second, pinned_secondary=pinned, springs=springs,
                               diaphragm=False, braced_side=False, torsion=False, n_combos=n_combos, buckling=buckling)
    return B
