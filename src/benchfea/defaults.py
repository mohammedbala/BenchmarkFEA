"""Built-in defaults: the single source used by the schema generator
(``default`` annotations), the model builder and the "defaults used" report.

Clause references are to ANSI/AISC 360-22 (Specification for Structural Steel
Buildings, August 1, 2022); page numbers are the printed "16.1-nn" numbers.
"""
from __future__ import annotations

FORMAT_VERSION = "1.0.0"          # version written by this engine
SUPPORTED_MAJOR = 1
SUPPORTED_MINOR = 0               # highest minor version of major 1 understood

VERTICAL_DETECTION_SINE = 1.0e-6  # member is "vertical" when |x_local x vertical| <= this

CONVENTIONS = {
    "handedness": "right",
    "rotations": "right-hand-rule",
    "memberLocalX": "i-to-j",
    "defaultOrientation": "y-toward-vertical",
    "memberEndForces": "on-member-local",
    "internalForces": "positive-face",
    "reactions": "on-structure-global",
    "loads": "along-positive-axes",
}

ANALYSIS = {
    "type": "firstOrder",
    "shearDeformation": True,
    "resultStations": 11,
    "convergenceTolerance": 1.0e-6,
    "maxIterations": 50,
    "alpha": 1.0,
    "equilibriumTolerance": 1.0e-6,
    "pivotTolerance": 1.0e-10,
    "geometryTolerance": 1.0e-9,
}
ANALYSIS_REASONS = {
    "type": "linear elastic first-order analysis unless a combination or this option asks for second order",
    "shearDeformation": "shear deformation included (AISC 360-22 C2.1(a) requires flexural, shear and axial deformations)",
    "resultStations": "11 equally spaced stations per member (tenth points)",
    "convergenceTolerance": "relative displacement-increment norm for iterative analyses",
    "maxIterations": "iteration cap for second-order and tension/compression-only analyses",
    "alpha": "load multiplier for second-order analysis: 1.0 (LRFD); set 1.6 for ASD (AISC 360-22 C2.1(d))",
    "equilibriumTolerance": "relative out-of-balance force accepted by the invariants report",
    "pivotTolerance": "relative pivot threshold below which the structure is reported as a mechanism",
    "geometryTolerance": "lengths below this fraction of the model size are treated as zero",
}
STIFFNESS_REDUCTION = {"factor": 0.8, "tauB": "computed"}
STIFFNESS_REDUCTION_REASONS = {
    "factor": "0.80 on all stiffnesses (AISC 360-22 C2.3(a), p. 16.1-30)",
    "tauB": "tau_b computed from Eq. C2-2a/C2-2b with Pns = Fy*A (AISC 360-22 C2.3(b), p. 16.1-30)",
}
BUCKLING = {"modes": 3}
NOTIONAL = {"coefficient": 0.002}
NOTIONAL_REASON = "Ni = 0.002*alpha*Yi (AISC 360-22 Eq. C2-1, p. 16.1-29)"

COMBINATION = {"purpose": "strength", "stiffnessReduction": False}

#: Design variables with a fixed built-in default: name -> (value, reason).
#: Variables whose default is derived from the model (Lx, Ly, Lz, Lb, An, role)
#: are handled in designparams.py and documented in the schema descriptions.
DESIGN = {
    "code": ("AISC 360-22", "only supported design code and edition"),
    "method": ("LRFD", "engine default; AISC 360-22 B3 permits LRFD or ASD and prescribes no default"),
    "Kx": (1.0, "effective length = unbraced length (AISC 360-22 C3, p. 16.1-31)"),
    "Ky": (1.0, "effective length = unbraced length (AISC 360-22 C3, p. 16.1-31)"),
    "Kz": (1.0, "effective length = unbraced length (AISC 360-22 C3, p. 16.1-31)"),
    "Cb": ("computed", "Cb from the moment diagram of each unbraced segment (AISC 360-22 Eq. F1-1, p. 16.1-52)"),
    "Cmx": ("computed", "Eq. A-8-4 when no transverse load acts on the segment, otherwise 1.0 "
                        "(AISC 360-22 App. 8.1.2, p. 16.1-284)"),
    "Cmy": ("computed", "Eq. A-8-4 when no transverse load acts on the segment, otherwise 1.0 "
                        "(AISC 360-22 App. 8.1.2, p. 16.1-284)"),
    "sideswayInhibited": (False, "engine default: member ends assumed free to translate (not prescribed by AISC 360-22)"),
    "U": (1.0, "engine default: all cross-section elements connected, no shear lag (not prescribed by AISC 360-22)"),
    "maxSlendernessCompression": (200.0, "Lc/r preferably not above 200 (AISC 360-22 E2 User Note, p. 16.1-40); advisory"),
    "maxSlendernessTension": (300.0, "L/r preferably not above 300 (AISC 360-22 D1 User Note, p. 16.1-32); advisory"),
    "deflectionLimits": ([], "no deflection check: AISC 360-22 L2 (p. 16.1-174) sets no numeric limit"),
}

DOFS = ("ux", "uy", "uz", "rx", "ry", "rz")
