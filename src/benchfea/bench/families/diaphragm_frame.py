"""Space frames with a rigid diaphragm at every floor (3D, eccentric lateral load, torsion).

The frame of ``multistorey_3d`` with all nodes of a floor slaved to a master node at the mass
centre, which is eccentric to the grid centre (and never on a grid line).  The wind acts on the
master, optionally with an accidental torsional moment; an optional braced bay in one face makes
the lateral stiffness unsymmetric as well.  All slaves are in the plane of their master.
"""
from __future__ import annotations

from ._common import BUCKLING_CHANCE, BUCKLING_FLAG, Builder, flag
from ._space import space_frame

NAME = "diaphragm_frame"
LEVELS = (3, 4, 5)
FEATURES = {
    "buckling": BUCKLING_FLAG,
    "secondOrder": "strength combinations analysed second order (opt-in: PyNite cannot follow, see the module)",
    "accidentalTorsion": "torsional moment at the master with each wind force (5 % of the plan dimension)",
    "bracedSide": "diagonal bracing in the face y = 0 (unsymmetric lateral stiffness)",
    "pinnedSecondary": "beams along y on interior lines pinned at both ends (needs two bays along x)",
    "baseSprings": "partially fixed bases (rotational springs about both horizontal axes)",
}
_SHAPES = {3: ((1, 1, 1), (1, 1, 2), (2, 1, 1)), 4: ((1, 1, 2), (2, 1, 1), (2, 1, 2), (2, 2, 1)),
           5: ((2, 1, 2), (2, 2, 1), (1, 1, 3), (2, 2, 2))}


def build(draw, level: int, unit_system: str, vertical: str, features: dict) -> Builder:
    B = Builder(unit_system, vertical, NAME)
    shape = draw.choice(_SHAPES[level])
    # second order with a diaphragm is OPT-IN (never drawn; the draw is still consumed): PyNite's
    # stiff-member emulation of the diaphragm does not reproduce the second-order result (measured:
    # 0.5-0.7 % against the engine and OpenSees, which agree with each other to 1e-4)
    second = flag(features, "secondOrder", draw.chance(0.0))
    torsion = flag(features, "accidentalTorsion", draw.chance(0.6))
    braced = flag(features, "bracedSide", draw.chance({3: 0.3, 4: 0.5, 5: 0.6}[level]))
    pinned = flag(features, "pinnedSecondary", draw.chance({3: 0.1, 4: 0.5, 5: 0.9}[level])) and shape[0] >= 2
    springs = flag(features, "baseSprings", draw.chance({3: 0.1, 4: 0.3, 5: 0.7}[level]))
    n_combos = {3: draw.integer(3, 4), 4: draw.integer(4, 7), 5: draw.integer(6, 8)}[level]
    buckling = flag(features, "buckling", draw.chance(BUCKLING_CHANCE[level]))
    B.parameters = space_frame(B, draw, shape, second=second, pinned_secondary=pinned, springs=springs,
                               diaphragm=True, braced_side=braced, torsion=torsion, n_combos=n_combos, buckling=buckling)
    return B
