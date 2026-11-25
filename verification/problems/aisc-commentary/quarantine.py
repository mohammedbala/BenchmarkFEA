"""Quarantined problems of verification/problems/aisc-commentary and their written diagnoses.

Every problem listed here REPRODUCES its two targets within the Commentary's tolerance; it is quarantined
because tier C (cross-check with the foreign solvers) does not pass on its model, and a problem is not
published when a tier that applies does not pass (no majority vote).  Nothing was changed in a model, a
target or a tolerance.  The numbers were measured with ``benchfea.verify.verify_problem``.
"""

_C = ("Both targets are reproduced within the source's tolerance, and the engine agrees with the ODE oracle on this "
      "model (tier A, 1e-6 along the member). Tier C does not pass: {found}. Tier C allows a difference of 1e-6 plus "
      "the ESTIMATED discretisation error of the foreign solver (OpenSeesPy with a linearised P-Delta transformation "
      "and a finite number of elements per member), which is an estimate and not a bound (docs/VERIFICATION.md, "
      "section 7); here the difference to OpenSeesPy exceeds the estimate by less than 1 %{pyn}. The engine's result "
      "is the exact solution of the second-order member equations (it is the ODE oracle's), so the disagreement lies "
      "in the foreign solver's discretisation; but under the rule that nothing is published when a tier disagrees, "
      "this axial load level is quarantined. The levels below it, and all levels of the model with shear deformation "
      "(which the foreign solvers cannot represent, so that tier C does not apply), are active.")
_PYN = ("; PyNite is not eligible at this load level (its estimated approximation error 6.45e-3 exceeds the limit "
        "5e-3), so that only one foreign solver remains")

QUARANTINE = {
    "case-1-us-without-shear-p300": _C.format(found="engine and OpenSeesPy differ by 1.04e-3 in the rotation of the "
                                              "pinned end (allowed 1.03e-3)", pyn=""),
    "case-1-us-without-shear-p450": _C.format(found="engine and OpenSeesPy differ by 1.85e-3 in the rotation of the "
                                              "pinned end (allowed 1.84e-3)", pyn=_PYN),
    "case-1-si-without-shear-p1330": _C.format(found="engine and OpenSeesPy differ by 1.03e-3 in the rotation of the "
                                               "pinned end (allowed 1.03e-3; exceeded in the next digit)", pyn=""),
    "case-1-si-without-shear-p2000": _C.format(found="engine and OpenSeesPy differ by 1.85e-3 in the rotation of the "
                                               "pinned end (allowed 1.84e-3)", pyn=_PYN),
    "case-2-us-without-shear-p150": _C.format(found="engine and OpenSeesPy differ by 7.74e-4 in the rotation of the "
                                              "top (allowed 7.73e-4)", pyn=""),
    "case-2-us-without-shear-p200": _C.format(found="engine and OpenSeesPy differ by 1.51e-3 in the rotation of the "
                                              "top (allowed 1.51e-3; exceeded in the next digit)", pyn=""),
    "case-2-si-without-shear-p667": _C.format(found="engine and OpenSeesPy differ by 7.75e-4 in the rotation of the "
                                              "top (allowed 7.73e-4)", pyn=""),
    "case-2-si-without-shear-p890": _C.format(found="engine and OpenSeesPy differ by 1.52e-3 in the rotation of the "
                                              "top (allowed 1.51e-3)", pyn=""),
}
