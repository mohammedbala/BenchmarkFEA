"""Quarantined problems of verification/problems/staad and their written diagnoses.

"<id>-vendor-printed" problems hold the values PRINTED by the vendor page that do not reproduce at their
printed precision (half a unit in the last printed digit); the problem "<id>" with the reference values
recomputed from the same statement stays active.  A plain "<id>" listed here is quarantined as a whole
(the engine produces no results for it, or another tier fails on its model).  Nothing was changed in
a model, a target or a tolerance.  The numbers quoted were measured with ../check_local_v2b.py and with
verification/derive/vendor_reference.py.
"""

QUARANTINE = {
    "frame-support-settlement":
        "All 12 displacement targets are reproduced (largest relative difference 2.0e-4, inside the printed "
        "precision), but the model cannot pass tiers C, D and E, so the problem stays quarantined. Cause: the frame "
        "is statically determinate, so a support settlement produces a rigid-body motion with NO member forces and "
        "NO reactions. Every force of the solution is round-off of zero (largest member end force or moment 1.3e-8 "
        "in kN and mm, largest reaction 2e-12 kN), and the checks that compare forces relative to the largest "
        "force of the model then compare one rounding error with another: tier C reports the end moments of the "
        "engine and of PyNite as differing by 1.7 and 1.8 (relative) and the OpenSees route fails its own joint "
        "equilibrium check with a residual of order one; tier D reports station forces that differ from the end "
        "forces by 6e-4 to 4e-3; the member-split check of tier E is above its tolerance. None of these numbers "
        "says anything about the displacements, and none shows an error of the engine or of a foreign solver. The "
        "reference solution (plane frame program of verification/derive/vendor_reference.py) gives the rigid-body "
        "displacements the page prints (1.25 cm, 0.001667 rad, 0.4167 cm) with forces of 1e-12 kN. History: an "
        "earlier build of the engine refused this model ('invariant failed'), for the same reason inside its own "
        "invariants; a round-off floor on the force scales of the invariants removed that refusal "
        "(docs/THEORY.md 8.1). The verification tiers were not changed to let the problem through, and no load was "
        "added to the model to make its forces non-zero.",
    "portal-pinned-midspan-load-vendor-printed":
        "Missed by 2.3e-5 relative (13 times the tolerance): engine -2755.035 lb in, printed 2,754.97. The page's "
        "value comes from a closed-form expression that neglects axial deformation, evaluated on the page with "
        "rounded section properties (it prints I = 10.67 and 1.333 in^4 and '= 2,754 in-lb' next to the formula, "
        "2,754.97 in its table). The same expression with the exact properties (I_b / I_c = 8) gives 2755.10 lb in; "
        "the frame of the statement, with axial deformation, gives 2755.03 lb in (reference program of "
        "verification/derive/vendor_reference.py = engine to 1e-9). The printed digits beyond the fourth are not "
        "significant. The page's other two values (500 and 27.55 lb) reproduce and are targets of the active problem.",
    "portal-fixed-uniform-load-vendor-printed":
        "Missed by 8.1e-4 relative (7 times the tolerance): engine -44.436 kip ft, printed 44.40. The page's 'theory' "
        "value is a hand calculation (the page cites an example of a textbook on approximate / iterative frame "
        "analysis) carried to three figures. The exact value for members without axial deformation is w L^2 / 18 = "
        "44.444 kip ft (reference program with EA x 1e6); with the axial stiffness of the statement 44.436 kip ft "
        "(reference program = engine to 1e-9). The program result printed on the same page is 44.44.",
    "two-storey-frame-lateral-loads-vendor-printed":
        "Missed by 9.1e-3 relative: engine 178.007 kip ft, printed 176.40. The page's 'theory' value is a hand "
        "calculation that neglects axial deformation: the reference program with EA x 1e6 gives 176.364 kip ft, which "
        "rounds to the printed value; with the axial stiffness of the statement (A = 0.2 ft^2) it gives 178.007 kip "
        "ft = engine to 1e-9. The program result printed on the same page is 178.01 (the page itself reports a "
        "difference of 0.9 %).",
    "cantilever-column-p-delta-vendor-printed":
        "Missed by 6.4e-4 relative: engine 737.194 kip in, printed 736.72. The page's reference is a one-step "
        "approximation, M = H L + P x (FIRST-order tip deflection) = 5 x 144 + 3 x 5.57; the exact second-order "
        "solution of the stated beam-column (ODE oracle) is M = 737.194 kip in with a tip deflection of 5.731 in. "
        "The approximation omits the growth of the deflection under the axial load (P / P_e = 0.028). The page's "
        "first-order tip deflection 5.57 in reproduces and is a target of the active problem.",
    "cantilever-column-notional-load-vendor-printed":
        "Both values missed by 2.0 times the tolerance: moment engine 61.380 kip in, printed 61.37; deflection engine "
        "1.1460 in, printed 1.145. Cause: rounding in the page's hand evaluation of the closed form. The page carries "
        "alpha = 1.226, tan(alpha)/alpha = 2.273 and M0 = 2.25 ft kip and prints '5.115 ft k = 61.37 in k' although "
        "5.115 x 12 = 61.38; for the deflection it multiplies 0.451 by 2.537. In double precision the same closed "
        "forms give alpha = 1.22642, M = 0.09 x 300 x 2.273332 = 61.380 kip in and delta = 0.451228 x 2.539731 = "
        "1.14600 in, the ODE oracle's and the engine's values. The program results printed on the same page are "
        "61.36 and 1.14538.",
}
