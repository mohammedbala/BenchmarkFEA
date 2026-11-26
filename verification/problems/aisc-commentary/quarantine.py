"""Quarantined problems of verification/problems/aisc-commentary and their written diagnoses.

No problem of this source is quarantined at present: every problem reproduces its two targets within the
Commentary's tolerance and passes every tier that applies to its model.  (The models with shear
deformation cannot be cross-checked -- no foreign solver represents second order with shear deformation --
which their model records state; their targets are tier B and pass.)  A problem that fails a tier is
entered here with its diagnosis and is then not published.
"""

QUARANTINE = {}
