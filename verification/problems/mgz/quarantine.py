"""Quarantined problems of verification/problems/mgz and their written diagnoses.

A problem is listed here when it does not reproduce within the written tolerance policy (half a unit in
the last printed digit of every target) or when another tier fails on its model.  Nothing was changed in
a model, a target or a tolerance to move a problem out of this list.  The numbers quoted were measured
with check_local.py, with ``benchfea verify`` and with small diagnosis scripts (solution of the source's
own printed equations; a plane-frame calculation with the source's polynomial geometric stiffness); the
diagnosis scripts are evidence for the diagnosis only, not a verification tier.
"""

_CONFIRMED = (" The engine's results for this model are confirmed by the other tiers that apply (tier C: OpenSeesPy "
              "and PyNite agree with the engine; tier D: equilibrium and the other invariants pass); the difference "
              "is between the source's printed digits and the exact solution of the stated problem.")
_HAND = ("Cause: rounding in the source's hand calculation - its stiffness coefficients are carried with three or four "
         "significant figures and later results are computed from rounded earlier ones, so the last printed digit of "
         "its results is not reliable, while the policy counts every printed digit as significant. ")

QUARANTINE = {
    "ex-2-2-truss-four-bar":
        "Displacement reproduced (2.5548 mm against the printed 2.55). The printed bar forces are missed by 1.05 and "
        "1.6 times the tolerance: engine 368.752 and -255.479 kN, printed 368.7 and -255.4 kN. " + _HAND +
        "Evidence: the source obtains the forces as products of intermediate values rounded to four digits "
        "(319.3 x 2/sqrt 3 and 180.6 x sqrt 2); the printed forces leave a residual of 0.20 kN in the horizontal "
        "equilibrium of joint a (2 x 368.7 cos 30 + 2 x 255.4 cos 45 = 999.80 kN against the load of 1000 kN), the "
        "engine's forces none." + _CONFIRMED,
    "ex-3-2-truss-three-bar":
        "6 of 9 targets missed, by up to 1.1e-3 relative (tolerances 1e-4 to 6e-4). Engine: u_a = 0.87018, "
        "v_a = 1.24297, u_b = -0.19298 mm; reactions -63.572, -383.000, -257.828 kN; bar forces 89.904, -63.572, "
        "515.657 kN. Printed: 0.871, 1.244, -0.193; -63.6, -383.4, -258.1; 90.0, -63.6, 516.2. " + _HAND +
        "Evidence: (1) the truss is statically determinate, so reactions and bar forces follow from statics alone: "
        "R_xc = -383.0 kN (the applied 383.0 kN; printed -383.4), R_yc = -257.83 kN (printed -258.1), "
        "F_ac = 2 x 257.83 = 515.66 kN (printed 516.2); the printed reactions leave residuals of 0.4 kN (x) and 0.3 kN "
        "(y). (2) The exact solution of the source's own printed 3 x 3 equations is 0.8702, 1.2431, -0.1931 mm - the "
        "engine's values, not the printed 0.871 and 1.244." + _CONFIRMED,
    "ex-3-3-truss-with-tie":
        "5 of 7 targets missed, by up to 4.2e-3 relative (tolerances 2e-4 to 1.3e-3). Engine: 0.38188, 1.22791, "
        "-0.43796 mm; reactions -144.273, -162.519, -177.127, -220.481 kN. Printed: 0.383, 1.228, -0.437; -144.3, "
        "-163.2, -177.4, -221.1. " + _HAND + "Evidence: the printed reactions leave residuals of 1.3 kN in x "
        "(383.0 - 163.2 - 221.1) and 0.3 kN in y; the exact solution of the source's own printed equations is 0.3819, "
        "1.2281, -0.4381 mm - the engine's values to the digits shown, not the printed 0.383 and -0.437." + _CONFIRMED,
    "ex-3-5-truss-five-bar":
        "All 13 targets missed, by 1e-4 to 3.1e-3 relative (tolerances 5e-5 to 1.2e-3). Engine: displacements -0.40824, "
        "9.81217, -2.23069, 10.92751, -7.80177 mm; reactions -282.840, -772.733, 1055.573 kN; bar forces -163.298, "
        "446.138, 892.276, -326.596, -772.733 kN. Printed: -0.407, 9.809, -2.232, 10.926, -7.801; -282.9, -772.0, "
        "1056.2; -162.8, 446.8, 891.4, -326.8, -773.2. " + _HAND + "Evidence: the truss is statically determinate "
        "(three reaction components, five bars, four joints), so reactions and bar forces follow from statics alone: "
        "R_xa = -282.84 kN, R_yb = 282.84 (8.660 + 10)/5 = 1055.57 kN, R_ya = -772.73 kN; the printed reactions leave a "
        "residual of 1.36 kN in the vertical direction (-772.0 + 1056.2 - 282.84). The exact solution of the source's "
        "printed (four-digit) equations is -0.4081, 9.8107, -2.2303, 10.9259, -7.8006 mm: within 3.2e-4 of the "
        "engine's values, and not the printed ones (-0.407, 9.809, -2.232)." + _CONFIRMED,
    "ex-3-6-truss-added-restraints":
        "u_d missed by 2.7 times the tolerance: engine 1.11534 mm, printed 1.114 mm (v_d reproduced: -1.90101 against "
        "-1.901). " + _HAND + "Evidence: the source multiplies the load by an inverse matrix printed with three digits "
        "((0.250 + 0.144) x 2.8284 = 1.114); the exact solution of its own printed stiffness equations is 1.1153 mm, "
        "the engine's value." + _CONFIRMED,
    "ex-4-9-beam-overhang":
        "The engine reproduces the printed displacements and the moment-diagram value at b (17.68 kN m; statics: "
        "5/sqrt 2 kN x 5 m = 17.678), but 3 of the 4 printed reactions are missed. Engine: R_xa = -3.5355 kN, "
        "R_ya = -3.3146 kN, R_mza = -8.8388 kN m, R_yb = 6.8501 kN. Printed: -3.60, -3.30, -8.80, 6.85. " + _HAND +
        "Evidence: the printed reactions do not satisfy equilibrium with the stated load - horizontal equilibrium "
        "requires R_xa = -5/sqrt 2 = -3.536 kN (printed -3.60, 1.8 % off), and moments about b require "
        "R_mza = 8 m x R_ya + 17.678 kN m, which the printed pair (-3.30, -8.80) misses by 0.08 kN m. The source obtains "
        "them by multiplying stiffness coefficients with displacements rounded to two or three significant figures "
        "(200 x 0.750 x 0.024 = 3.60)." + _CONFIRMED,
    "ex-4-10-beam-joint-moment":
        "The rotation and the five reactions are reproduced; the two member end moments at b are missed by 1.14 times "
        "the tolerance: engine 35714.3 and 14285.7 kN mm (exactly 5/7 and 2/7 of the applied 50 kN m), printed 35.72 "
        "and 14.28 kN m. " + _HAND + "Evidence: the source computes them from the rotation rounded to four digits "
        "(200 x 1.0e5 x 0.001786 = 35.72; the unrounded rotation is 1/560 = 0.0017857 rad, giving 35.714)."
        + _CONFIRMED,
    "ex-4-11-beam-settlement":
        "6 of 8 targets reproduced; R_ya and R_yb are missed by 1.2 and 2.4 times the tolerance: engine 15.1339 and "
        "-36.6482 kN, printed 15.14 and -36.66 kN. " + _HAND + "Evidence: the source's coefficients 0.00469 and "
        "0.00949 are 12EI/L^3/200 = 0.0046875 and 0.0094875 rounded to three digits; with them "
        "200 x (0.00949 x (-20) + 6.750 x 0.0009643) = -36.658 (printed -36.66), with the unrounded coefficients "
        "-36.648 (engine). The printed vertical reactions sum to -0.01 kN instead of zero." + _CONFIRMED,
    "ex-4-12-beam-bending-torsion":
        "6 of 9 targets missed. Engine: v_b = -0.54573 mm, theta_xb = 0.00904348 rad, R_ya = 0.41295 kN, "
        "R_mza = 1783.4 kN mm, R_yc = 0.58705 kN, R_mzc = -1415.0 kN mm. Printed: -0.545, 0.009044, 0.423, 1.781 kN m, "
        "0.586, -1.413 kN m. Two causes. (1) Suspected misprint: the printed reaction vector gives R_ya = 0.423 kN, "
        "but the source's own figure of the deflected structure labels this reaction .413 kN, and vertical equilibrium "
        "with the 1 kN load requires R_ya + R_yc = 1 (0.423 + 0.586 = 1.009; engine 0.41295 + 0.58705 = 1.00000). "
        "(2) " + _HAND + "The remaining misses are 6e-5 to 1.8e-3 relative (v_b from the three-digit coefficient "
        "0.00949; theta_xb = 40/(200 x 22.115) with the torsional stiffness rounded to five digits; reactions from the "
        "rounded displacements)." + _CONFIRMED,
    "ex-4-13-frame-two-member":
        "u_b and v_b missed by 1.4 and 1.8 times the tolerance: engine 0.44147 and -0.39989 mm, printed 0.4414 and "
        "-0.3998 mm (rotation reproduced). Cause: the printed values are truncated, not rounded. Evidence: the exact "
        "solution of the source's own printed 3 x 3 equations is 0.44147, -0.39988 mm, 0.0016943 rad - the engine's "
        "values; rounded to four digits they would print as 0.4415 and -0.3999. The policy (half a unit in the last "
        "printed digit) assumes rounding." + _CONFIRMED,
    "ex-4-15-beam-on-elastic-supports":
        "Both targets are reproduced by the model as written (support area entered as the source's defining expression "
        "22.5 x 15 / 29,000 in^2: deflection -0.23795 in, moment 547.22 in-kips; printed 0.238 and 547). The problem is "
        "nevertheless quarantined because the outcome depends on how a rounded input of the source is read: with the "
        "area as the source prints it, 0.0116 in^2 (22.43 kip/in instead of the stated 22.5 kip/in per support), the "
        "engine gives -0.23855 in and 547.65 in-kips, which miss the printed values by 1.1 and 1.3 times the "
        "tolerance. The source does not say which value its computer run used, and both targets are labels on figures. "
        "An input that decides the outcome within its own printed precision is recorded as an open question, not "
        "resolved by choosing the reading that passes.",
    "ex-5-3-gable-frame-stiffness":
        "16 of 18 stiffness coefficients reproduced. Missed: K(M_b, v_b) and K(M_c, v_c), by 1.4 times the tolerance: "
        "engine 200 x 17.3817 and 200 x (-7.7817), printed 17.381 and -7.781. Cause: a last-digit inconsistency in the "
        "source. Evidence: the coefficient is (6EI/L^2) cos(phi) = 18.75 x sqrt(55)/8 = 17.38171 (x 200); the source "
        "prints 17.381 in Example 5.3 and 17.382 for the same coefficient of the same member in Example 5.7 (p. 115)."
        + _CONFIRMED,
    "ex-5-4-space-truss":
        "8 of 15 targets missed, by up to 2.4e-3 relative for u_a (engine 0.17787 mm, printed 0.1783 mm) and up to "
        "7.5e-4 for the reactions; v_a reproduced (2.72196 against 2.722), w_a -0.48652 against -0.4863. " + _HAND +
        "Evidence: the printed displacements are the exact solution of the source's printed 3 x 3 matrix (0.17823, "
        "2.72215, -0.48631), whose coefficients are sums of products rounded to four digits (for example 2.5215 "
        "entered as 2.522); the small displacement u_a is sensitive to that rounding. The printed reactions leave a "
        "residual of 0.2 kN in the y direction (sum -599.8 kN against the load of 600 kN); the engine's none."
        + _CONFIRMED,
    "ex-5-6-beam-member-loads":
        "6 of 9 targets reproduced (reactions, support moment, end moments at b). Missed: theta_a (engine "
        "-5.6812e-4 rad, printed -5.684e-4), theta_b (0.69565e-4, printed 0.698e-4) and the moment under the 20 kN load "
        "(11297 kN mm, figure label 11.29 kN m), by 5.7, 4.7 and 1.5 times the tolerance. " + _HAND + "Evidence: the "
        "source solves with the fixed-end moments rounded to 10.67 and 3.73 kN m (exact 10.6667 and 3.7333): "
        "(-1.4 x 10.67 + 0.5 x 3.73)/1.15 x 5e-5 = -5.684e-4 and (0.5 x 10.67 - 3.73)/1.15 x 5e-5 = 0.698e-4; with the "
        "unrounded moments the same expressions give -5.6812e-4 and 0.6957e-4, the engine's values. The label 11.29 is "
        "6.87 x 3 - 9.32 computed from the rounded reactions." + _CONFIRMED,
    "ex-5-7-frame-inclined-member":
        "7 of 9 targets missed, by 0.1 % to 0.5 % (tolerances 5e-5 to 9e-4). Engine: u_b = 0.99500 mm, "
        "v_b = -4.98158 mm, theta_b = -5.3423e-4 rad; reactions 130.500 kN, 55.676 kN, 13374.6 kN mm, -149.250 kN, "
        "22.674 kN, -45356.6 kN mm. Printed: 0.9982 (0.09982 in the solution line, a misprint - 0.9982 is used in the "
        "reaction product), -4.996, -0.000534; 131.0, 55.4, 13.43e3, -149.7, 22.7, -45.41e3. " + _HAND + "Evidence: the "
        "printed displacements are the exact solution of the source's printed 3 x 3 matrix (1.394, 0.259, 7.031; 0.114, "
        "1.368; 2e5), but its coefficient 0.114 stands for (EA/L) sin^2 + (12EI/L^3) cos^2 + 12EI/L^3 = 0.10950 + "
        "0.00469 = 0.11419; the source rounds 0.1095 to 0.109 and the sum 0.11369 to 0.114, and v_b is sensitive to "
        "this coefficient. The printed vertical reactions leave a residual of 0.25 kN (55.4 + 22.7 - 46.35 - 32)."
        + _CONFIRMED,
    "ex-5-8-grid":
        "7 of 9 targets missed, by 1.6e-4 to 1.0e-3 relative (tolerances 8e-5 to 8e-4). Engine: w_b = -22.3845 mm, "
        "theta_xb = 4.19677e-3, theta_yb = 5.93197e-3 rad; R_za = 14.7524 kN, R_mxa = -6.4566 kN mm, "
        "R_mya = -36245.0 kN mm, R_zc = 5.2476 kN, R_mxc = -41974.2 kN mm, R_myc = -17.1114 kN mm. Printed: -22.37, "
        "4.195e-3, 5.931e-3; 14.74, -6.45, -36.21e3, 5.25, -41.94e3, -17.11. " + _HAND + "Evidence: the exact solution "
        "of the source's printed 3 x 3 equations (with 12EI/L^3 terms rounded to 0.00949 for 0.0094875) is "
        "w_b = -22.374 mm, the printed value; the printed vertical reactions sum to 19.99 kN against the load of "
        "20 kN." + _CONFIRMED,
    "ex-5-10-truss-lack-of-fit":
        "6 of 7 targets missed, by 1e-4 to 2.4e-3 relative (tolerances 6e-5 to 1e-3). Engine: -1.68343, -0.05192, "
        "-0.84458 mm; reactions R_yb = -278.223, R_xc = 760.120, R_yc = 278.223, R_xd = -760.120 kN. Printed: -1.683, "
        "-0.0518, -0.8448; -278.1, 760.0, 278.1, -760.2. " + _HAND + "Evidence: there is no external load, so "
        "horizontal equilibrium requires R_xc = -R_xd; the printed pair 760.0 and -760.2 cannot both be right, and the "
        "engine's 760.12 lies between them. The exact solution of the source's printed equations is -1.68347, "
        "-0.05183, -0.84473 mm." + _CONFIRMED,
    "ex-9-5-frame-critical-load":
        "Source value is a discretised approximation; no stated error bound. Engine (exact stability functions): "
        "6609.41 kN; printed 6630 kN (one element per member with the polynomial geometric stiffness of the source's "
        "Eq. 9.14); difference 0.31 %, tolerance from the printed precision 7.5e-5. Evidence that the difference is the "
        "source's discretisation error: a plane-frame calculation with the source's own element matrices (Eq. 4.34 and "
        "Eq. 9.14) gives 6629.6 kN with one element per member - the printed value - and 6610.8 and 6609.39 kN with two "
        "and eight elements per member, converging to the engine's value. (That calculation is a diagnosis aid, not a "
        "verification tier.) Buckling factors cannot be cross-checked by the foreign solvers of tier C.",
    "ex-9-5-frame-second-order":
        "All three targets missed by 3 % to 4.5 %: engine -971.71, 908.02, 280.31 kN m; figure labels of the source "
        "932, 869, 272 kN m. The source's values come from an incremental analysis with one element per member whose "
        "solution scheme and step size are not stated; no error bound is given. Evidence that the difference is not "
        "the source's discretisation: the source's own element matrices (Eq. 4.34 with the geometric stiffness of "
        "Eq. 9.14, one element per member), solved here to equilibrium with the member axial forces iterated, give "
        "971.5, 907.8 and 280.3 kN m - the engine's values within 0.03 % - so the difference of 4 % lies in the source's "
        "solution process (an incremental scheme without equilibrium iteration underestimates the displacements of a "
        "softening structure; the source discusses this drift in its Section 12.2) or in its geometry updating, which "
        "cannot be determined from the example. The other tiers that apply pass on this model (tier C: largest difference "
        "to OpenSeesPy 1.5e-4 against an allowed 4.4e-4, to PyNite 2.1e-5 against 3.2e-5; tiers D and E); the problem "
        "is quarantined for its published targets only. Read off a figure; approximate method in the source.",
    "ex-9-6-portal-second-order":
        "All three targets missed: engine -550.62, 9098.77, -2664.06 in-kips; figure labels of the source 519, 9078, "
        "2661 in-kips (6.1 %, 0.23 %, 0.12 %). The source's values come from an incremental analysis (two elements per "
        "column) whose solution scheme and step size are not stated; no error bound is given. Evidence that the "
        "difference is not the source's discretisation: the source's own element matrices (Eq. 4.34 and Eq. 9.14, two "
        "elements per column, one per girder segment), solved here to equilibrium, give 550.6, 9098.6 and 2664.5 "
        "in-kips - the engine's values within 0.02 %. The first-order moments of the same example (316, 8708, 2044) and "
        "its critical load factor (2.20) are reproduced (problems ex-9-6-portal-first-order and "
        "ex-9-6-portal-critical-load). The other tiers that apply pass on this model (tier C: largest difference to "
        "OpenSeesPy 1.9e-4 against an allowed 5.6e-4, to PyNite 1.9e-5 against 3.8e-5; tiers D and E); the problem is "
        "quarantined for its published targets only. Read off a figure; approximate method in the source.",
    "ex-9-8-column-braced-critical-load":
        "Source value is a discretised approximation; no stated error bound. Engine (exact stability functions): "
        "214.05 kips = pi^2 EI/L^2 with L = 156 in (buckling in two half waves, as the source's text states for a "
        "sufficiently stiff brace); label on the source's plot 215.7 kips (four elements with the polynomial geometric "
        "stiffness of Eq. 9.14); difference 0.77 %, tolerance from the printed precision 2.3e-4. Evidence: a plane-frame "
        "calculation with the source's own element matrices gives 215.66 kips with four elements - the printed value - "
        "and 214.06 kips with sixteen. Buckling factors cannot be cross-checked by the foreign solvers of tier C.",
    "ex-13-2-fixed-beam-three-loads":
        "5 of 6 targets reproduced (v1, theta1, v2, theta2, v3). Missed: theta3, engine 1.4375e-3 rad, printed "
        "1.48e-3 rad (2.9 % apart, tolerance 0.34 %). Suspected misprint or arithmetic slip in the source. Evidence: "
        "the source's own substructure equations give 1.4375e-3 - with its printed matrices, "
        "theta3 = -(1/4)(1.5e-3 x v2 + theta2) = -(1/4)(1.5e-3 x (-4.0) + 0.25e-3) = 1.4375e-3 rad for the printed "
        "v2 = -4.0 mm and theta2 = 0.25e-3 rad." + _CONFIRMED,
    "ex-13-3-three-storey-frame":
        "6 of 9 targets missed, by 2e-4 to 1.2e-3 relative (tolerances 9e-5 to 4e-4). Engine (member areas 1e8 mm^2 "
        "for the neglected axial deformation): u3 = 11.3718 mm, theta3 = -3.1419e-4, theta4 = -3.6513e-4, u1 = 14.0514 "
        "mm, theta1 = -0.5727e-4, theta2 = -1.5706e-4, u5 = 6.0426 mm, theta5 = -5.0750e-4, theta6 = -5.7658e-4 rad. "
        "Printed: 11.36, -3.141e-4, -3.652e-4, 14.05, -0.573e-4, -1.570e-4, 6.04, -5.070e-4, -5.759e-4. " + _HAND +
        "Evidence: the source's condensed stiffness coefficients carry three or four digits (0.00355, 4.961, "
        "15.257e4); the exact solution of its printed assembled equations is u3 = 11.3665 mm, which the source prints "
        "truncated as 11.36. The differences are not caused by the area used here for 'axial deformation neglected': "
        "the engine's values change by less than 1e-6 relative between areas of 1e8 and 1e9 mm^2." + _CONFIRMED,
    "ex-13-11-gable-frame":
        "7 of 9 targets missed, by 0.3 % to 0.6 % (tolerances 4e-4 to 4e-3). Engine: u_b = 4.8510, v_b = -0.02196 mm, "
        "theta_b = -1.3448e-3; u_c = 7.0575, v_c = -5.5590 mm, theta_c = 0.36612e-3; u_d = 9.2485, v_d = -0.04054 mm, "
        "theta_d = -0.12831e-3 rad. Printed: 4.88, -0.02, -1.351e-3; 7.09, -5.58, 0.368e-3; 9.28, -0.04, -0.129e-3. "
        "Cause: the source's hand solution of its own equations is not exact to the digits printed. Evidence: the exact "
        "solution of the printed symmetric-part equations is u_b = -2.198, v_c = -5.557 mm (the source prints -2.208 "
        "and -5.583) and of the printed antisymmetric-part equations u_b = 7.050, u_c = 7.058 mm (the source prints "
        "7.085 and 7.093); combined they give 4.852, 7.058, -5.557, 9.248 mm, which agree with the engine within "
        "0.04 %." + _CONFIRMED,
}
