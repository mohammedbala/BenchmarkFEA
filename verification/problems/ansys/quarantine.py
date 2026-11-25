"""Quarantined problems of verification/problems/ansys and their written diagnoses.

"<id>-vendor-printed" problems hold the "Target" values PRINTED by the vendor page that do not reproduce
at their printed precision (half a unit in the last printed digit); the problem "<id>" with the
reference values recomputed from the same statement stays active.  Nothing was changed in a model, a
target or a tolerance.  The numbers quoted were measured with ../check_local_v2b.py and with
verification/derive/vendor_reference.py.
"""

QUARANTINE = {
    "vm2-beam-with-overhangs-vendor-printed":
        "Missed by 3.5e-4 relative: engine moment -6,000,000 lb in, i.e. a stress of -11,403.95 psi; printed target "
        "-11400.000 psi. The page prints a value rounded to three significant figures (evidently the textbook's "
        "rounded answer) with three decimals, so the printed-precision rule assigns it the floor tolerance 1e-6, which "
        "its real precision (5e-3 x 11,400 ... half a unit of the third figure) does not support. The exact value "
        "follows from statics: the moment in the middle portion is w a^2 / 2 = (10000/12) x 120^2 / 2 = 6.0e6 lb in, "
        "the stress 6.0e6 x 15 / 7892 = 11,403.95 psi (ODE oracle = engine to 1e-15). The page's deflection target "
        "(0.182 in) reproduces and is a target of the active problem.",
    "vm21-tie-rod-lateral-load-vendor-printed":
        "Both printed moments are missed at the boundary of the tolerance (1.015 and 1.0002 times half a unit of the "
        "last printed digit). (1) Stiffened: engine 4580.1507 lb in, printed 4580.1 - the page's value is the exact "
        "value cut off after the first decimal rather than rounded (4580.15 rounds to 4580.2); difference 0.0507 "
        "against the tolerance 0.05. (2) Unstiffened: the exact value is p L^2 / 8 = 1.79253 x 200^2 / 8 = 8962.65 lb "
        "in, exactly half a unit from the printed 8962.7; the comparison in floating point (8962.7 is not "
        "representable) falls on the failing side by 2e-11 of the value. The ODE oracle and the engine agree to "
        "1e-12 on 4580.1507 and 8962.65. The page's other four values (deflections and end slopes, five significant "
        "figures) reproduce and are targets of the active problem.",
}
