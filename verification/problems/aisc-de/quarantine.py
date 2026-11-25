"""Quarantined problems of verification/problems/aisc-de and their written diagnosis.

All ten problems are quarantined for the same reason: the source's model is not completely specified
(results of unnamed software; shear areas and the treatment of shear deformation are not given), and the
results do not reproduce at the printed precision with the model that CAN be built from what the source
gives.  Nothing was changed in a model, a target or a tolerance; in particular no shear area was taken
from a shape table or fitted into the models.  Numbers measured with ../check_local_v2b.py and with the
parameter study described in README.md.
"""

_COMMON = (
    "Cause: the source's model is not completely specified, and the missing part matters at the printed precision. "
    "The source's figures are output of unnamed analysis software; the source gives Ag and Ix of the W12x65 and Ix "
    "of the W18x40 but no shear areas, and does not say whether shear deformation is included. The model built here "
    "from the given data (no shear deformation) is systematically stiffer: first-order drift 1.8 % to 2.2 % smaller, "
    "base shears and column moments 0.4 % to 1.0 % larger than printed, vertical reactions reproduced (they follow "
    "from statics). Evidence that shear deformation in the source's software is the cause (parameter study, "
    "README.md; Example C.1B, first order, LRFD: printed drift 0.145 in, base shears 5.64 and 6.21 kips; without "
    "shear deformation 0.1422 in, 5.671 and 6.247 kips): with shear deformation and shear areas of about 4.5 in^2 "
    "(column) and 5 to 6 in^2 (beam) - values of the size of the web areas of such shapes, NOT taken from any table "
    "and NOT entered in the model - the same analysis gives 0.145 in, 5.64 and 6.21 kips. The beam area, which the "
    "source does not give either, is not the cause: between 5 in^2 and an axially rigid beam the average drift of "
    "grid lines B and C changes by less than 1e-5 in and the base shears by 0.15 %. The second-order to first-order "
    "drift ratios of the engine (C.1A LRFD 0.28935 / 0.17770 = 1.63, source 1.60; C.1B LRFD 1.45, source 1.41) are "
    "consistent with the same difference in first-order stiffness. Because the shear areas are not given by the "
    "source and no opened reference supplies them, the example cannot be completed without a guess: not usable as "
    "a target at the printed precision; kept, evaluated and reported, never counted as passed.")

_MISSED = {
    "c1a-first-order-lrfd": "4 of 7", "c1a-first-order-asd": "5 of 7", "c1a-second-order-lrfd": "4 of 7",
    "c1a-second-order-asd": "5 of 7", "c1b-first-order-lrfd": "4 of 7", "c1b-first-order-asd": "5 of 7",
    "c1b-second-order-lrfd": "4 of 7", "c1b-second-order-asd": "5 of 7", "c1c-first-order-lrfd": "4 of 7",
    "c1c-first-order-asd": "5 of 7",
}
_ASD = (" ASD: the leaning-column load is entered as printed in the source's figure (154 kips; 60.0 ft x 2.56 kip/ft "
        "is 153.6 kips).")

QUARANTINE = {pid: f"{n} targets missed (drift, base shears, and column moments where printed with three figures). "
              + _COMMON + (_ASD if pid.endswith("asd") and not pid.startswith("c1c") else "")
              for pid, n in _MISSED.items()}
