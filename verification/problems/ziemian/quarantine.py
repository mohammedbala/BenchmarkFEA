"""Quarantined problems of verification/problems/ziemian and their written diagnoses.

A problem is listed here when it does not reproduce within the written tolerance policy (half a unit in
the last digit the data set's spreadsheet displays) or when the engine cannot solve it.  Nothing was
changed in a model, a target or a tolerance to move a problem out of this list.  The numbers quoted were
measured with ../check_local_v2b.py and with verification/derive/ziemian_diagnosis.py (own plane-frame
program with the data set's discretisation; OpenSeesPy with the linearised and with the corotational
transformation); those are evidence for the diagnosis only, not a verification tier.
"""

_SECOND = (
    "Cause: different theory. The data set's WC values come from MASTAN2's incremental work-control analysis "
    "(step 0.001) in an updated Lagrangian formulation, 4 elements per member: a large-displacement analysis. The "
    "engine solves the classical small-displacement second-order theory (P-Delta and P-delta with exact stability "
    "functions; docs/THEORY.md 10 and 12), which neglects terms of second order in the kinematics such as the "
    "shortening of the chord of a bent member. Evidence (verification/derive/ziemian_diagnosis.py, case 'Right', "
    "joint displacements): OpenSeesPy with the linearised 'PDelta' transformation, 32 elements per member, agrees "
    "with the engine within {ep}; OpenSeesPy with the 'Corotational' (large-displacement) transformation agrees with "
    "the data set within {cw} (32 elements per member; {c1} with the data set's own 4 elements per member), while "
    "the engine differs from the data set by up to {ew}. The remaining difference between the corotational results "
    "and the data set (element formulation, increment size) was not resolved. The targets are outside what the "
    "engine's theory can reproduce at the printed precision; they are kept, evaluated and reported, and never "
    "counted as passed.")

_LOADS = (
    "Cause: precision of the input listing, not of the analysis. The data set's model file lists loads (and "
    "coordinates) with five significant figures, while its results were computed with the unrounded model and are "
    "displayed to 0.01 kip in, i.e. to 1e-6 ... 1e-5 of the moments; a load rounded in the fifth figure changes a "
    "moment by up to 5e-5 of its value. Evidence: the first-order joint displacements of the same model, displayed "
    "with three decimals, are reproduced; ")

_LBA = (
    "Cause: discretisation of the data set. Its ratio comes from an eigenvalue analysis with 4 elements per member "
    "and a polynomial geometric stiffness, which overestimates the critical load; the engine's members carry the "
    "exact stability functions. Evidence (verification/derive/ziemian_diagnosis.py): a plane-frame program written "
    "here with the data set's own elements (cubic elements, consistent geometric stiffness, axial forces of a "
    "linear analysis, perfect geometry, 0.8 E) gives {own}, the data set's stored value {ds} within {d:.0e}; the "
    "engine gives {eng}, lower by {lo:.1e} of the value, which is more than half a unit of the third decimal the "
    "data set displays ({shown}).")

QUARANTINE = {
    # ---------------------------------------------------------------- critical load ratios
    "frame-1-critical-load-ratio": _LBA.format(own="7.856995", ds="7.857009", d=1.7e-6, eng="7.856261", lo=9.5e-5,
                                               shown="7.857 against the engine's 7.856"),
    "frame-7-critical-load-ratio": _LBA.format(own="6.353614", ds="6.353606", d=1.2e-6, eng="6.353220", lo=6.1e-5,
                                               shown="6.354 against the engine's 6.353"),
    # ---------------------------------------------------------------- first-order moments
    "frame-1-first-order-moments-right":
        "One of three targets missed by 1.24 times the tolerance: member 2 (right column) engine 2131.4762 kip in, "
        "displayed 2131.47 (stored 2131.47251). " + _LOADS + "every first-order result of this frame (displacements "
        "and moments, both cases) exceeds the stored value by the same factor 1 + 1.73e-6, as a load factor would: "
        "the loads listed as 34.849 and 348.49 kip are not exactly the loads of the analysis.",
    "frame-2-first-order-moments-right":
        "2 of 4 targets missed, by up to 1.5e-5 relative (engine 778.888 and 3948.680 kip in, displayed 778.90 and "
        "3948.69; tolerance 0.005 kip in). " + _LOADS + "the differences of the four moments are -1.5e-5 ... "
        "-1.7e-6, not a common factor (nine listed loads, each rounded separately).",
    "frame-2-first-order-moments-left":
        "2 of 4 targets missed, by up to 9.5e-6 relative (engine 1013.750 and 4293.270 kip in, displayed 1013.76 and "
        "4293.28; tolerance 0.005 kip in). " + _LOADS + "see frame-2-first-order-moments-right.",
    "frame-9-first-order-moments-right":
        "7 of 10 targets missed, by up to 2.6e-5 relative (largest: member 2, engine 2974.823 kip in, displayed "
        "2974.90; tolerance 0.005 kip in). " + _LOADS + "the differences are 2e-6 ... 2.6e-5 with both signs "
        "(eighteen listed loads, each rounded separately).",
    "frame-9-first-order-moments-left":
        "8 of 10 targets missed, by up to 2.5e-5 relative (largest: member 3, engine 3037.076 kip in, displayed "
        "3037.15; tolerance 0.005 kip in). " + _LOADS + "see frame-9-first-order-moments-right.",
    # ---------------------------------------------------------------- first-order displacements
    "frame-9-first-order-displacements-right":
        "One of six targets missed by 1.002 times the tolerance, at a rounding boundary: joint 21, stored value "
        "-0.2375013 in, displayed -0.238; engine -0.2374990 in. The engine differs from the stored value by 2.3e-6 in "
        "(9.5e-6 relative, the size of the rounding of the listed loads, see frame-9-first-order-moments-right), "
        "which moves it to the other side of the rounding boundary -0.2375 of the displayed number; the policy "
        "compares with the number as displayed. The other five displacements are reproduced.",
    "frame-3-first-order-displacements-right":
        "One of four targets missed (joint 18: engine 0.57205 in, displayed 0.571, stored 0.571297); all four differ "
        "from the stored values by +3.1e-4, -2.7e-4, -1.4e-3, +1.3e-3 relative, two orders of magnitude more than for "
        "the other frames (1e-5). Cause NOT identified. What was established: the stiffness model is the data set's "
        "(the critical load ratio of the same model is reproduced to 6e-5, and to 6e-6 by the cubic-element program "
        "of verification/derive/ziemian_diagnosis.py relative to the engine); the differences alternate in sign "
        "between the left and the right column line, the pattern of a different distribution of the loads between "
        "the two column lines (the lateral displacements of this frame depend on the difference of the two gravity "
        "loads, listed as 235.25 and 239.88 kip, through the diagonals), not of a stiffness difference. Variants "
        "tried without success: E = 29,000 ksi for the pin-ended members, the two pin-ended sections exchanged. "
        "Open question for the data set, see README.md.",
    "frame-3-first-order-displacements-left":
        "One of four targets missed (joint 17: engine -0.56181 in, displayed -0.563, stored -0.562549). Same finding "
        "and same open question as frame-3-first-order-displacements-right: the first-order displacements of frame 3 "
        "differ from the data set by up to 1.4e-3 relative with alternating sign between the column lines; cause not "
        "identified.",
    # ---------------------------------------------------------------- second order
    "frame-1-second-order-right":
        "All 5 targets missed: displacements by 2.4e-3 and 3.3e-3 relative (engine 1.19383 and 1.16788 in, data set "
        "1.190821 and 1.163823), moments by 3e-5 ... 1e-3. " + _SECOND.format(ep="4.8e-5", cw="1.4e-3", c1="1.6e-3", ew="3.5e-3"),
    "frame-1-second-order-left":
        "All 5 targets missed: displacements by 2.7e-3 and 1.8e-3 relative, moments by 1.5e-4 ... 1.7e-3. "
        + _SECOND.format(ep="4.8e-5", cw="1.4e-3", c1="1.6e-3", ew="3.5e-3"),
    "frame-2-second-order-right":
        "All 7 targets missed: displacements by 1.6e-3 ... 1.8e-2 relative (eaves: engine 6.9466 and 9.5419 in, data "
        "set 7.0707 and 9.6204; ridge 12.1274 against 12.1078), moments by 2e-3 ... 1.6e-2. "
        + _SECOND.format(ep="3.0e-5", cw="1.5e-3", c1="1.2e-3", ew="1.8e-2"),
    "frame-2-second-order-left":
        "All 7 targets missed: displacements by 2.2e-3 ... 2.1e-2 relative, moments by 2.2e-3 ... 1.3e-2. "
        + _SECOND.format(ep="3.0e-5", cw="1.5e-3", c1="1.2e-3", ew="1.8e-2"),
    "frame-3-second-order-right":
        "All 4 targets missed by 1.1e-2 ... 1.2e-2 relative (engine 1.9003 ... 1.9683 in, data set 1.8783 ... 1.9467; "
        "critical load ratio of this frame 1.387, amplification about 3.6). "
        + _SECOND.format(ep="3.0e-4", cw="6.6e-3", c1="1.3e-2", ew="1.2e-2")
        + " In addition the first-order displacements of this frame are not reproduced (frame-3-first-order-"
          "displacements-right).",
    "frame-3-second-order-left":
        "All 4 targets missed by 1.05e-2 ... 1.2e-2 relative. "
        + _SECOND.format(ep="3.0e-4", cw="6.6e-3", c1="1.3e-2", ew="1.2e-2")
        + " In addition the first-order displacements of this frame are not reproduced (frame-3-first-order-"
          "displacements-left).",
    "frame-5-second-order-right":
        "All 4 targets missed by 2.2e-3 ... 1.6e-2 relative (largest: joint 25, engine 2.1776 in, data set 2.2139). "
        + _SECOND.format(ep="3.8e-5", cw="2.1e-3", c1="1.3e-3", ew="1.6e-2"),
    "frame-5-second-order-left":
        "All 4 targets missed by 4.6e-4 ... 1.3e-2 relative (largest: joint 25, engine -2.1120 in, data set -2.0854). "
        + _SECOND.format(ep="3.8e-5", cw="2.1e-3", c1="1.3e-3", ew="1.6e-2"),
    "frame-7-second-order-right":
        "The engine produces no results: it refuses the combination as beyond the elastic critical load (buckling "
        "factor about 0.85). Cause: in the case 'Right' the only diagonal brace (one pin-ended element, W4x13, "
        "375 in long) is in compression: first-order force 21.68 kip against its Euler load pi^2 (0.8 E) I / L^2 = "
        "18.45 kip (ratio 1.18; verification/derive/ziemian_diagnosis.py). The engine's pin-ended member is a "
        "beam-column that buckles at its Euler load; the data set's single pin-ended element has no member buckling "
        "(the data set reports displacements 2.340 and 2.342 in "
        "for this case). The description file of the data set states that only the tension brace is modelled, i.e. "
        "the frame is meant for the case 'Left'. Not usable as a target for this engine.",
    "frame-7-second-order-left":
        "3 of 4 targets missed: joints 14 and 18 by 5e-4 and 1.3e-3 relative (engine -2.4002 and -2.4059 in, data set "
        "-2.398659 and -2.408860), joint 29 by 3.0e-2 (engine -0.1891, data set -0.1947 in; a small displacement of "
        "the braced bay). " + _SECOND.format(ep="(case 'Right' not solvable by the engine, see frame-7-second-order-right; "
                                             "for the other frames 3e-5 ... 3e-4)", cw="1.1e-3", c1="1.3e-3",
                                             ew="3.0e-2 (case 'Left')"),
    "frame-9-second-order-right":
        "All 16 targets missed: displacements by 8e-3 ... 2.8e-2 relative (e.g. joint 13: engine -0.43277 in, data "
        "set -0.421399), moments by up to 1.1e-2. This frame carries gravity loads only; its small sway (0.4 ... 0.8 "
        "in) is sensitive to the horizontal movements caused by the shortening of the chords of the heavily loaded "
        "beams. " + _SECOND.format(ep="1.2e-4", cw="4.3e-3", c1="4.3e-3", ew="2.7e-2"),
    "frame-9-second-order-left":
        "All 16 targets missed: displacements by 4e-3 ... 1.4e-2 relative, moments by up to 2.6e-3. Gravity loads only, "
        "see frame-9-second-order-right. " + _SECOND.format(ep="1.2e-4", cw="4.3e-3", c1="4.3e-3", ew="2.7e-2"),
}
