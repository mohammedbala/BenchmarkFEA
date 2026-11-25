# Vendor verification problems: ANSYS Mechanical APDL Verification Manual

Source of the problem **statements**: ANSYS, Inc., *Mechanical APDL Verification Manual*, Release 2025 R2 (public
help site ansyshelp.ansys.com, pages `Hlp_V_VM<n>.html`, read 2026-11-22; (c) ANSYS, Inc., all rights reserved:
cited, not reproduced). Each page states a test case, cites a reference, and prints a "Target" value next to the
program's value. The page address is in every `problem.json`.

**The references the pages cite were not opened for this project**, and the problem sketches of the pages were not
readable. Therefore:

* a page is used only where its TEXT (test case, modelling notes, formulas) determines the supports and the load
  positions; how it does is written in the problem's `analysis`;
* basis of every problem: `recomputed-from-published-statement`. Reference values are computed here
  (`verification/derive/vendor_reference.py`: ODE oracle, plane frame / truss program, characteristic equation of
  the Euler column); tolerance 1e-6;
* the page's printed Target values are recorded as additional targets with the printed-precision tolerance; those
  that do not reproduce are moved to `<id>-vendor-printed`, quarantined with a diagnosis;
* the program's own results are never targets.

The problems are written by `build_problems.py` (helpers in `../vendor_common.py`).

## Result: 8 problems (6 statements): 6 active, 2 quarantined

| page | problem | what |
|---|---|---|
| VM2 | `vm2-beam-with-overhangs` | beam with loaded overhangs: moment, deflection |
| VM4 | `vm4-two-bar-truss` | two-bar truss: bar force, deflection |
| VM21 | `vm21-tie-rod-lateral-load` | tie rod: second-order analysis in TENSION, and first order |
| VM41 | `vm41-stiff-beam-on-rotational-spring` | rotational spring support |
| VM127 | `vm127-hinged-bar-buckling` | Euler load |
| VM217 | `vm217-portal-frame-uniform-load` | fixed-base portal (base fixity deduced from the page's two formulas) |

Quarantined, printed values that do not reproduce (the recomputed targets of the same statements are active):

* `vm2-beam-with-overhangs-vendor-printed`: stress target printed as -11400.000 psi - a three-figure value printed
  with three decimals; exact value of the statement 11,403.95 psi.
* `vm21-tie-rod-lateral-load-vendor-printed`: the two printed moments, 4580.1 (exact 4580.1507: cut off, not
  rounded) and 8962.7 (exact 8962.65: exactly half a unit away; the floating-point comparison fails by 2e-11).

Diagnoses: `quarantine.py`.

## Excluded pages

| page | reason |
|---|---|
| VM1 | load positions only in the sketch (the same textbook problem is in `../staad/axially-loaded-column`, where the page's input listing gives them) |
| VM10, VM12 | cross-section stress problems, no frame result beyond statics |
| VM14 | large-deflection analysis of an eccentrically loaded strut |
| VM18 | curved bar |
| VM136 | elastica (large rotations) |
| VM216 | lateral buckling of a right-angle frame of plate-like members |
| VM222 | warping torsion |
| VM50, VM59 and the other dynamic cases | dynamics |
| VM180 and the remaining pages of the manual | not examined |

## How to read the index

`targets` / `missed` and the differences are those of `../check_local_v2b.py` (engine against the targets of
`problem.json`); "largest difference / tolerance" above 1 means a target is missed. `status` is the declared
status. "Verified" in this project means only: passed the tiers listed in `docs/VERIFICATION.md` at the
tolerances stated there. The other tiers (cross-check with foreign solvers, invariants, sensitivity) are run by
`benchfea verify verification/problems/ansys/<id>`.

Tolerance policy (docs/VERIFICATION.md, section 5), applied uniformly and checked by `../check_local_v2b.py`:
a published number carries half a unit in its last printed digit (floor 1e-6); a number computed here
carries 1e-6. A problem that misses, or on whose model another tier does not pass, is quarantined with a written
diagnosis (`quarantine.py`); nothing is loosened or dropped.

## Index

| problem | source location | what it verifies | basis | targets | missed | largest rel. difference | largest difference / tolerance | status |
|---|---|---|---|---|---|---|---|---|
| `vm127-hinged-bar-buckling` | VM127, 'Buckling of a Bar with Hinged Ends (Line Elements)' | Slender bar with hinged ends: critical load | recomputed-from-published-statement | 2 | 0 | 3.7e-06 | 0.28 | active |
| `vm2-beam-with-overhangs` | VM2, 'Beam Stresses and Deflections' | Beam with two loaded overhangs: moment in the middle portion and midspan deflection | recomputed-from-published-statement | 5 | 0 | 2.5e-03 | 0.93 | active |
| `vm2-beam-with-overhangs-vendor-printed` | VM2, 'Beam Stresses and Deflections' | Beam with two loaded overhangs: moment in the middle portion and midspan deflection - values printed by the vendor that do not reproduce | published | 1 | 1 | 3.5e-04 | 346.79 | quarantined |
| `vm21-tie-rod-lateral-load` | VM21, 'Tie Rod with Lateral Loading' | Tie rod with uniform lateral load: deflection, end slope and moment, with and without the tensile force | recomputed-from-published-statement | 10 | 0 | 1.5e-05 | 0.72 | active |
| `vm21-tie-rod-lateral-load-vendor-printed` | VM21, 'Tie Rod with Lateral Loading' | Tie rod with uniform lateral load: deflection, end slope and moment, with and without the tensile force - values printed by the vendor that do not reproduce | published | 2 | 2 | 1.1e-05 | 1.01 | quarantined |
| `vm217-portal-frame-uniform-load` | VM217, 'Portal Frame Under Symmetric Loading' | Fixed-base portal frame with a stiff span under uniform load: joint rotation and largest moment | recomputed-from-published-statement | 7 | 0 | 2.0e-03 | 0.96 | active |
| `vm4-two-bar-truss` | VM4, 'Deflection of a Hinged Support' | Two hinged bars under a load at their common joint: bar force and deflection | recomputed-from-published-statement | 4 | 0 | 5.8e-16 | 0.00 | active |
| `vm41-stiff-beam-on-rotational-spring` | VM41, 'Small Deflection of a Rigid Beam' | Very stiff beam pinned on a weak rotational spring: small rotation | recomputed-from-published-statement | 6 | 0 | 1.7e-06 | 0.00 | active |
