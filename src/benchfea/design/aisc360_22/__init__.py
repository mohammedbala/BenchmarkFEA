"""ANSI/AISC 360-22, Specification for Structural Steel Buildings (August 1, 2022).

One module per chapter; every equation is cited by number in the code and in
the ``clause`` of the checks it produces.

==================== ====================================================
module               content
==================== ====================================================
``b4_classification`` B4.1, Tables B4.1a / B4.1b
``d_tension``         D1 (advisory), D2, D3 (Ae = An U)
``e_compression``     E2 (advisory), E3, E4 (Eq. E4-2), E7
``g_shear``           G2.1, G4, G5, G6
``f_flexure``         F1 (Cb), F2, F3, F6, F7, F8
``h_combined``        H1.1, H1.2 (Eq. H1-1a/b); torsion screening (H3 is
                      not implemented); formulas of App. 8.1.2 (B1)
``serviceability``    deflection limits of the user (Chapter L has none)
==================== ====================================================

Design basis B3: LRFD (Eq. B3-1) or ASD (Eq. B3-2), applied by
:func:`benchfea.design.base.strength_check`.  The safety factors are the
printed values (1.67, 2.00, 1.50), not 1.5/phi.
"""
from ..registry import Edition, register
from . import b4_classification, d_tension, e_compression, f_flexure, g_shear, h_combined, serviceability

CODE = "AISC 360-22"

EDITION = register(Edition(
    code=CODE,
    title="ANSI/AISC 360-22 Specification for Structural Steel Buildings",
    limit_states=(d_tension.LIMIT_STATES + e_compression.LIMIT_STATES + g_shear.LIMIT_STATES
                  + f_flexure.LIMIT_STATES + h_combined.LIMIT_STATES),
    service=(serviceability.deflection_checks,)))

__all__ = ["CODE", "EDITION", "b4_classification", "d_tension", "e_compression", "f_flexure", "g_shear",
           "h_combined", "serviceability"]
