# Vendor verification problems: STAAD.Pro verification examples

Source of the problem **statements**: Bentley Systems, *STAAD.Pro Help*, "Verification Examples" (public
documentation, docs.bentley.com; the pages of the STAAD.Pro 2023 sample files, read 2026-11-22; (c) Bentley
Systems, no reuse licence: cited, not reproduced). Each page states a problem, cites a textbook, prints a "theory"
value and the program's value, and lists the program input. The page title and address are in every
`problem.json`.

**The textbooks the pages cite were not opened for this project.** Therefore:

* basis of every problem: `recomputed-from-published-statement`. The statement (geometry, properties, loads,
  supports) is taken from the page - where its text and its input listing differ, the listing defines the problem
  and the difference is noted in the problem;
* the reference values are computed here, by means that share no code with the engine
  (`verification/derive/vendor_reference.py`: the independent ODE oracle for single members, a plane frame / truss
  and a space frame / truss direct-stiffness program written for this purpose, the characteristic equation for the
  Euler load); tolerance 1e-6;
* the page's printed "theory" values are recorded as additional targets with the printed-precision tolerance.
  Those that reproduce are targets of the same problem. Those that do not are moved to a separate problem
  `<id>-vendor-printed`, **quarantined with a diagnosis**, so that the recomputed targets stay active;
* STAAD.Pro's own results are never targets.

The problems are written by `build_problems.py` (helpers in `../vendor_common.py`).

## Result: 24 problems (19 statements): 16 active, 8 quarantined

Active: every recomputed target is reproduced to 1e-6 or better and the page's printed values at their printed
precision.

Quarantined, values printed by the vendor that do not reproduce (5; the recomputed targets of the same statements
are active):

| problem | printed | exact solution of the statement | cause |
|---|---|---|---|
| `portal-pinned-midspan-load-vendor-printed` | 2,754.97 lb in | 2755.03 | closed form without axial deformation, evaluated with rounded properties |
| `portal-fixed-uniform-load-vendor-printed` | 44.40 kip ft | 44.436 | three-figure hand calculation (44.444 without axial deformation) |
| `two-storey-frame-lateral-loads-vendor-printed` | 176.40 kip ft | 178.007 | hand calculation without axial deformation (176.36 with EA -> infinity) |
| `cantilever-column-p-delta-vendor-printed` | 736.72 kip in | 737.194 | one-step approximation M = H L + P x (first-order deflection) |
| `cantilever-column-notional-load-vendor-printed` | 61.37 kip in, 1.145 in | 61.380, 1.1460 | rounding in the page's evaluation of the closed form |

Quarantined as a whole (3):

* `beam-column-uniform-load`, `cantilever-column-notional-load`: all targets reproduced, engine = ODE oracle;
  tier C does not pass (the difference to OpenSeesPy exceeds, by less than 1 %, the tolerance derived from the
  estimated discretisation error of the foreign solver).
* `frame-support-settlement`: **the engine produces no results.** A statically determinate frame under a support
  settlement moves as a rigid body with no forces; the engine's invariant check measures equilibrium residuals
  relative to the forces present, which are round-off of zero, and refuses the solution ("invariant failed",
  residuals of order 0.3). This is a limitation of the engine's invariant check found by this problem (reported;
  no load was added to make it run).

Diagnoses with numbers: `quarantine.py`.

## Excluded pages

| page | reason |
|---|---|
| 2D Portal Reactions Sidesway 2 | section taken from the program's shape table (W12X26), properties not printed on the page |
| 2D Frame 2 Step P-Delta Displacement | sections from the shape table, self weight, a two-iteration P-Delta as "theory" |
| Column Pushover Displacement | inelastic |
| Axial Forces on a 3D Space Model | the labels of the result tables were lost in the saved page text (which value belongs to which member / reaction) |
| 3D Frame Max Forces: largest moments about local y and z | depend on the vendor's local-axis convention for an inclined member, not stated on the page (the largest axial force of the same page IS used: `space-frame-three-members`) |
| Beam on Elastic Foundation | elastic foundation; foundation modulus only in the figure |
| Bent Cantilever Deflection | geometry and load only in the figure |
| Thermal Loading on a Beam, Bent Beam Thermal Loading, Axial Forces in a Plane Frame 2, Stresses in a Cable due to Thermal Loading | thermal loading (outside the engine) |
| Curved Beam, Stresses in a Circular Beam | curved members |
| Stresses in a Tapered Cantilever, Twist in a Tapered Tube, End Moments in a Non Uniform Beam | non-prismatic members: not examined |
| 3x2 Plane Frame Moments | complete statement with 28 end moments from a hand method (differences to the program up to 2.2 % on the page); not transcribed in this pass |
| Forces on a Propped Cantilever 1, Tee Shaped Cantilever, Hanging Bar Axial Stress, Reactions in a 2D Truss Model 1 and 2, Axial Force in a 2D Plane Frame 1, Roof Truss Axial Forces, Axial Force on a Cable | not examined in this pass |
| V.04 plates and shells, V.05 solids, V.06 loading, V.08 dynamics, V.09 - V.11 design | outside the scope (frame statics, second order, buckling) |

## How to read the index

`targets` / `missed` and the differences are those of `../check_local_v2b.py` (engine against the targets of
`problem.json`); "largest difference / tolerance" above 1 means a target is missed. `status` is the declared
status. "Verified" in this project means only: passed the tiers listed in `docs/VERIFICATION.md` at the
tolerances stated there. The other tiers (cross-check with foreign solvers, invariants, sensitivity) are run by
`benchfea verify verification/problems/staad/<id>`.

Tolerance policy (docs/VERIFICATION.md, section 5), applied uniformly and checked by `../check_local_v2b.py`:
a published number carries half a unit in its last printed digit (floor 1e-6); a number computed here
carries 1e-6. A problem that misses, or on whose model another tier does not pass, is quarantined with a written
diagnosis (`quarantine.py`); nothing is loosened or dropped.

## Index

| problem | source location | what it verifies | basis | targets | missed | largest rel. difference | largest difference / tolerance | status |
|---|---|---|---|---|---|---|---|---|
| `axially-loaded-column` | 'V. Axially Loaded Column' | Bar held at both ends with two axial loads: support reactions | recomputed-from-published-statement | 6 | 0 | 1.9e-16 | 0.00 | active |
| `beam-column-uniform-load` | 'V. Direct Analysis of a Beam' | Simply supported beam-column with uniform load: second-order moment and deflection | recomputed-from-published-statement | 4 | 0 | 1.1e-04 | 0.19 | quarantined |
| `beam-partial-trapezoidal-load` | 'V. Deflection and Reactions in a Beam' | Pinned-fixed beam with a linearly varying load on part of the span | recomputed-from-published-statement | 9 | 0 | 1.8e-02 | 0.80 | active |
| `bent-cantilever-with-prop` | 'V. Support Reactions for a Simple Frame' | Bent cantilever with an intermediate lateral support: support reactions | recomputed-from-published-statement | 5 | 0 | 1.6e-14 | 0.00 | active |
| `cantilever-column-notional-load` | 'V. Direct Analysis of a Column' | Cantilever column with axial load and a small lateral load: second-order moment and deflection | recomputed-from-published-statement | 4 | 0 | 3.2e-11 | 0.00 | quarantined |
| `cantilever-column-notional-load-vendor-printed` | 'V. Direct Analysis of a Column' | Cantilever column with axial load and a small lateral load: second-order moment and deflection - values printed by the vendor that do not reproduce | published | 2 | 2 | 8.7e-04 | 2.00 | quarantined |
| `cantilever-column-p-delta` | 'V. Single Column P-Delta Analysis' | Cantilever column with axial and lateral load: second-order base moment | recomputed-from-published-statement | 5 | 0 | 3.0e-04 | 0.34 | active |
| `cantilever-column-p-delta-vendor-printed` | 'V. Single Column P-Delta Analysis' | Cantilever column with axial and lateral load: second-order base moment - values printed by the vendor that do not reproduce | published | 1 | 1 | 6.4e-04 | 94.85 | quarantined |
| `frame-support-settlement` | 'V. 2D Portal Reactions Sidesway 1' | Statically determinate frame: displacements caused by support settlements | recomputed-from-published-statement | 12 | - | - | - | quarantined (the engine refuses the model) |
| `pinned-column-buckling-factor` | 'V. Column Buckling Factor' | Pin-ended column: elastic buckling factor | recomputed-from-published-statement | 2 | 0 | 8.7e-05 | 0.61 | active |
| `portal-fixed-uniform-load` | 'V. 2D Portal Reactions 2' | Fixed-base square portal frame with a uniform load on the beam | recomputed-from-published-statement | 5 | 0 | 9.0e-16 | 0.00 | active |
| `portal-fixed-uniform-load-vendor-printed` | 'V. 2D Portal Reactions 2' | Fixed-base square portal frame with a uniform load on the beam - values printed by the vendor that do not reproduce | published | 1 | 1 | 8.1e-04 | 7.22 | quarantined |
| `portal-pinned-midspan-load` | 'V. 2D Portal Reactions 1' | Pinned-base portal frame with a load at midspan of the beam | recomputed-from-published-statement | 7 | 0 | 1.3e-05 | 0.07 | active |
| `portal-pinned-midspan-load-vendor-printed` | 'V. 2D Portal Reactions 1' | Pinned-base portal frame with a load at midspan of the beam - values printed by the vendor that do not reproduce | published | 1 | 1 | 2.3e-05 | 12.91 | quarantined |
| `propped-cantilever-self-weight` | 'V. Forces on a Propped Cantilever 2' | Propped cantilever under its own weight | recomputed-from-published-statement | 12 | 0 | 5.7e-04 | 0.80 | active |
| `space-frame-three-members` | 'V. 3D Frame Max Forces' | Space frame of three members under joint loads | recomputed-from-published-statement | 18 | 0 | 2.8e-04 | 0.08 | active |
| `space-truss-three-bars` | 'V. Axial Stress on a Truss Model' | Three-bar space truss: bar forces | recomputed-from-published-statement | 7 | 0 | 6.3e-04 | 0.39 | active |
| `stepped-shaft-torsion` | 'V. Torsion on a Stepped Cantilever' | Stepped cantilever shaft under three torques: angle of twist | recomputed-from-published-statement | 4 | 0 | 8.6e-04 | 0.74 | active |
| `truss-five-bars-thrust` | 'V. Reactions in a 2D Truss Model 3' | Statically indeterminate five-bar truss: horizontal reaction | recomputed-from-published-statement | 5 | 0 | 4.6e-04 | 0.81 | active |
| `truss-thirteen-bars-deflection` | 'V. Deflections in a 2D Truss Model' | Thirteen-bar truss with an overhang: joint deflection | recomputed-from-published-statement | 6 | 0 | 1.3e-04 | 0.07 | active |
| `two-bar-truss` | 'V. Stress in a 2D Truss Model' | Two-bar truss: bar force and joint deflection | recomputed-from-published-statement | 4 | 0 | 5.8e-16 | 0.00 | active |
| `two-bay-frame-beam-forces` | 'V. 2 Bay Frame Moments and Shear' | Two-bay frame with beams pinned to the outer columns: beam moment and shear | recomputed-from-published-statement | 11 | 0 | 6.8e-07 | 0.00 | active |
| `two-storey-frame-lateral-loads` | 'V. 1x2 Plane Frame Lateral Load' | One-bay two-storey frame under lateral joint loads | recomputed-from-published-statement | 6 | 0 | 6.1e-14 | 0.00 | active |
| `two-storey-frame-lateral-loads-vendor-printed` | 'V. 1x2 Plane Frame Lateral Load' | One-bay two-storey frame under lateral joint loads - values printed by the vendor that do not reproduce | published | 1 | 1 | 9.1e-03 | 321.42 | quarantined |
