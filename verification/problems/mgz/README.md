# Published verification problems: McGuire, Gallagher & Ziemian

Source of every problem in this directory: W. McGuire, R. H. Gallagher, R. D. Ziemian, *Matrix Structural
Analysis*, 2nd edition (copyright 2014 R. D. Ziemian). Nothing of the book is reproduced here. Each
`problem.json` cites the example or equation and the page, gives in this project's own words the input data
that define the problem, and lists the printed target numbers; `model.json` is a normal benchfea input file
in the book's units (kN and mm, or kip and in).

"Reproduced" / "passed" below means only: the engine's results agree with the listed targets within the
tolerances stated, and the other tiers that apply to the model did not disagree
([docs/VERIFICATION.md](../../../docs/VERIFICATION.md)). Nobody has reviewed or approved these files.

## Counts (current state)

| | number |
|---|---|
| problems encoded | 46 |
| active - every target within the tolerance policy, no tier disagreeing | 20 |
| quarantined - with a written diagnosis in `problem.json` | 26 |
| examples of the book not encoded (not applicable, excluded, nothing to compare) | listed below |

Of the 26 quarantined problems:

* 18 miss one or more targets only in the last printed digit(s) (differences of 6e-5 to 6e-3) because of
  rounding, truncation or a slip in the book's hand calculation; the evidence (equilibrium not satisfied by the
  printed values, the exact solution of the book's own printed equations, internal inconsistencies of the book)
  is in each diagnosis. Two of them also contain a suspected misprint (Example 4.12, reaction 0.423 against
  0.413 in the book's own figure; Example 13.2, rotation 1.48e-3 against 1.4375e-3 from the book's own
  equations).
* 2 compare the engine's exact-member critical load with a value the book computed with a stated number of
  elements and an approximate geometric stiffness (Examples 9.5 and 9.8): "source value is a discretised
  approximation; no stated error bound".
* 2 compare with second-order moments the book shows on figures, obtained by an incremental analysis whose
  scheme is not stated (Examples 9.5 and 9.6); differences of up to 6 %.
* 2 reproduce their (closed-form) targets to 1e-14 but fail tier C on the model, in the way described in
  docs/VERIFICATION.md under "Second-order cross-checks" (Examples 8.3 and 8.4).
* 1 reproduces its targets, but the outcome depends on how a rounded input of the book is read (Example 4.15).
* 1 misses two of eighteen stiffness coefficients in the last digit because the book prints the same coefficient
  as 17.381 and as 17.382 in two examples (Example 5.3).

Largest differences (engine against the targets, see the index): among the 20 active problems the largest
ratio of difference to tolerance is 0.83 and the largest relative difference 4.6e-3 (Example 2.1, a displacement
printed as 0.72 mm, tolerance 6.9e-3). Among the quarantined problems the largest relative difference of a
missed target is 6.1 % (Example 9.6, second-order moment at b) and the largest ratio of difference to tolerance
is 79 (Example 9.5, second-order moments); in the 18 problems quarantined for rounding in the book the missed
targets differ by 6e-5 to 6e-3, except the two suspected misprints (2.4 % and 2.9 %) and the reactions of
Example 4.9 (0.4 % and 1.8 %).

## Files

| file | purpose |
|---|---|
| `<problem-id>/model.json`, `<problem-id>/problem.json` | the problems (format: docs/VERIFICATION.md, section 3) |
| `build_problems.py`, `problems_ch2_5.py`, `problems_ch7_13.py` | the definitions the directories are written from; a published target is entered as printed and its tolerance computed by `benchfea.verify.printed_precision_tolerance` |
| `quarantine.py` | the register of quarantined problems with their diagnoses |
| `check_local.py` | self-check with the engine only: `python3 verification/problems/mgz/check_local.py [-v] [--table]` |

The full verification of a problem (all tiers) is `benchfea verify verification/problems/mgz/<problem-id>`.

## Rules followed

1. **Targets as printed.** `value` = the printed number (times a unit factor such as 1000 for kN m to kN mm, or
   the common factor 200 of a printed matrix); `precisionNote` records the printed form.
2. **Tolerance** = half a unit in the last printed digit, as a relative tolerance, floor 1e-6 - the one policy of
   docs/VERIFICATION.md. No tolerance was chosen per case; `check_local.py` recomputes every tolerance from the
   printed form and fails if one differs.
3. **Recomputed targets** (`basis: recomputed-from-published-statement`): where the book gives a closed-form
   expression or a symbolic matrix and no numbers, a numeric instance is chosen here, the printed expression is
   evaluated in double precision, and the tolerance is the floor 1e-6 (Examples 4.14, 8.3, 8.4; Equations 13.31
   to 13.33 and 13.38 to 13.40). Exact fractions printed by the book (Example 13.12) also carry the floor.
4. **A problem that misses is quarantined** with a diagnosis; no model, target or tolerance was changed to make
   a problem pass. The status of a problem is the status of all its targets: one missed target quarantines the
   problem.
5. **Expression before rounded display.** Where the book gives a quantity both as an expression and as a rounded
   number, the expression is the input: coordinates from the stated lengths and angles (4/tan 30 deg rather than
   6.928 m; sqrt 55 m rather than 7.416 m), load components 5/sqrt 2 kN, and so on. Where the book works with
   stated rounded load components (383.0 and 321.4 kN; 18.75 and -46.35 kN; 282.84 kN), those are the input.
6. **Discretised values of the book are not targets when the book also prints the theoretical value.** For
   Examples 9.3, 9.4 and 12.6-12.9 the target is the theoretical critical load the book prints; its own one- or
   two-element results (1,137; 746; 9651; 6884 kN) are recorded in the analysis text and not compared. Where the
   book prints only a discretised value, that value is the target, and the problem is quarantined unless the
   engine's exact-member result agrees within the printed precision (it does for Example 9.6, lambda = 2.20, and
   for the unbraced column of Example 9.8, 53.5 kips).

## Modelling conventions used in the models

* **Plane problems** (x-y plane of the book) are modelled in 3D with `verticalAxis: "Y"`: the book's x, y,
  theta_z are the model's X, Y, rz, counter-clockwise positive. With the default member orientation in-plane
  bending is about local z for every member, so the book's I is `Iz`. `Iy`, `J` (and `A` where the book gives
  none) are placeholders that do not influence the targets.
* **Space problems** (Examples 5.4, 5.8) use `verticalAxis: "Z"`, the book's axes unchanged.
* **Signs.** Reactions: the book's R = force by the support on the structure = the engine's reactions. Member
  end moments in the book's element convention = `endForces` of the engine. Bar forces: tension positive.
  Moment-diagram labels are magnitudes in the book; the sign in the engine's convention is derived by statics
  and written in the target's `sign` note.
* **Plane trusses**: members of type `truss`; the out-of-plane translation of every joint is held.
* **Axial deformation neglected by the book** (Examples 5.9, 13.3): member area 1e8 mm^2 for every member, one
  value for all such problems. Measured influence of the remaining axial flexibility: the targets change by less
  than 7e-7 relative between 1e8 and 1e9 mm^2 (with the areas of the real members, 4e3 mm^2, they would change by
  up to 1.8 %).
* **Plane stability problems** (Chapters 8, 9, 12): the out-of-plane DOFs of every node are held and
  Iy = 1000 Iz, so that no out-of-plane mode precedes the in-plane modes compared.
* **Axial-force member in a buckling analysis** (Example 9.1): the engine treats a compressed truss member as a
  pin-ended beam-column with the section's I. The compressed bar is given I = 1e10 mm^4 so that its own Euler
  load (1.2e6 kN) stays out of the comparison; with a placeholder I = 1 mm^4 the engine reports the member's
  Euler factor 1.2e-4 as the lowest mode.
* **Prescribed displacements and stiffness coefficients**: a prescribed joint displacement is a support with a
  `supportDisplacement` load (Examples 2.3, 4.11, 5.11); a stiffness coefficient is the reaction at a fully held
  node under an imposed unit displacement (Examples 4.14, 5.3; Equations 13.31, 13.38).
* **Lack of fit** (Example 5.10: bar 3 mm too short): the input format has no initial strain; the exactly
  equivalent loading is used (the support at the far end of the bar moved 3 mm), see the problem's analysis text.
* **Inclined roller** (Example 13.5): supports act along global axes only; the model is written in the book's
  joint axes of the roller (all data rotated by -30 deg). Only u_b'' is then directly comparable (-0.0942 mm
  against the printed -0.094). Rotated back, the engine gives u_a = 0.9282 mm and v_a = 1.1425 mm; the book
  prints 0.928 and 1.143 (not targets, because printed values would have to be transformed).

## Index

Differences as measured by `check_local.py --table` (engine against the targets; "largest difference /
tolerance" above 1 means a target is missed).

| problem | source location | pages | what it verifies | targets | missed | largest rel. difference | largest difference / tolerance | status |
|---|---|---|---|---|---|---|---|---|
| `eq-13-31-flexible-connections` | Equations 13.31, 13.32, 13.33 (Examples 13.6 and 13.7) | pp. 394-398 | Beam with flexible (semi-rigid) and pinned end connections: stiffness coefficients | 17 | 0 | 6.0e-16 | 0.00 | active |
| `eq-13-38-rigid-end-offsets` | Equations 13.37 to 13.40 (Example 13.8) | pp. 398-401 | Member with rigid end offsets (finite joint size): stiffness coefficients | 18 | 0 | 7.3e-16 | 0.00 | active |
| `ex-12-7-cantilever-column-critical-loads` | Examples 12.6, 12.7 and 12.9: theoretical values | pp. 362-363, 367 | Cantilever column: first and second elastic critical load | 2 | 0 | 3.0e-04 | 0.44 | active |
| `ex-13-10-three-span-beam` | Example 13.10 | pp. 406-407 | Three-span beam with fixed ends and two point loads | 4 | 0 | 1.7e-03 | 0.67 | active |
| `ex-13-11-gable-frame` | Example 13.11 | pp. 407-408 | Symmetric gable frame under lateral and ridge loads | 9 | 7 | 9.8e-02 | 12.35 | quarantined |
| `ex-13-12-spring-chain-changed` | Example 13.12, summary result (changed structure) | pp. 411-412 | Chain of axial springs: flexibility matrix of the changed structure | 9 | 0 | 6.1e-16 | 0.00 | active |
| `ex-13-12-spring-chain-original` | Example 13.12, Equation 13.1 of the example (flexibility of the original structure) | pp. 411-412 | Chain of axial springs: flexibility matrix of the original structure | 9 | 0 | 5.2e-16 | 0.00 | active |
| `ex-13-2-fixed-beam-three-loads` | Example 13.2 | pp. 380-381 | Fixed-ended beam with three point loads | 6 | 1 | 2.9e-02 | 8.50 | quarantined |
| `ex-13-3-three-storey-frame` | Example 13.3 | pp. 383-385 | Three-storey single-bay rigid frame, axial deformation neglected | 9 | 6 | 1.2e-03 | 13.64 | quarantined |
| `ex-13-5-truss-inclined-roller` | Example 13.5 (truss of Examples 3.1 and 3.2) | pp. 390-391 | Three-bar plane truss on an inclined roller | 1 | 0 | 1.9e-03 | 0.36 | active |
| `ex-2-1-truss-two-bar` | Example 2.1 | pp. 21-22 | Two-bar plane truss: joint displacement under a horizontal load | 4 | 0 | 4.6e-03 | 0.66 | active |
| `ex-2-2-truss-four-bar` | Example 2.2 | pp. 22-23 | Statically indeterminate four-bar plane truss | 5 | 4 | 1.9e-03 | 1.58 | quarantined |
| `ex-2-3-truss-prescribed-displacement` | Example 2.3 | pp. 23-24 | Force needed for a prescribed joint displacement of a two-bar truss | 2 | 0 | 1.1e-04 | 0.76 | active |
| `ex-3-2-truss-three-bar` | Example 3.2 (truss of Example 3.1) | pp. 35-38 | Three-bar plane truss: displacements, reactions, bar forces | 9 | 6 | 1.1e-03 | 10.87 | quarantined |
| `ex-3-3-truss-with-tie` | Example 3.3 | pp. 38-39 | Three-bar plane truss with an added horizontal tie | 7 | 5 | 4.2e-03 | 13.62 | quarantined |
| `ex-3-5-truss-five-bar` | Example 3.5 (truss of Example 3.4) | pp. 42-46 | Five-bar plane truss: displacements, reactions, bar forces | 13 | 13 | 3.1e-03 | 17.51 | quarantined |
| `ex-3-6-truss-added-restraints` | Example 3.6 | p. 46 | Five-bar plane truss with added horizontal restraints | 2 | 1 | 1.2e-03 | 2.69 | quarantined |
| `ex-4-10-beam-joint-moment` | Example 4.10 | pp. 81-82 | Two-span beam with fixed ends, moment applied at the interior support | 8 | 2 | 1.0e-03 | 1.14 | quarantined |
| `ex-4-11-beam-settlement` | Example 4.11 | pp. 82-83 | Two-span beam with fixed ends, settlement of the interior support | 8 | 2 | 4.0e-04 | 2.36 | quarantined |
| `ex-4-12-beam-bending-torsion` | Example 4.12 | pp. 83-84 | Fixed-ended beam with an eccentric load: bending and torsion | 9 | 6 | 2.4e-02 | 20.10 | quarantined |
| `ex-4-13-frame-two-member` | Example 4.13 | pp. 84-85 | Two-member rigid plane frame with a loaded joint | 3 | 2 | 2.6e-03 | 1.76 | quarantined |
| `ex-4-14-beam-internal-hinge` | Example 4.14 (stiffness matrix of the combined member) | pp. 86-87 | Two beam segments joined by a hinge: condensed stiffness coefficients | 8 | 0 | 1.4e-15 | 0.00 | active |
| `ex-4-15-beam-on-elastic-supports` | Example 4.15, part 2 (stiffness method) | pp. 87-88 | Beam on an elastic foundation modelled by 19 axial members | 2 | 0 | 4.1e-04 | 0.45 | quarantined |
| `ex-4-9-beam-overhang` | Example 4.9 | pp. 80-81 | Propped cantilever with overhang, inclined end load | 10 | 3 | 1.8e-02 | 12.89 | quarantined |
| `ex-5-10-truss-lack-of-fit` | Example 5.10 (truss of Example 3.3) | pp. 123-124 | Plane truss with a bar fabricated 3 mm too short | 7 | 6 | 2.4e-03 | 4.45 | quarantined |
| `ex-5-11-beam-settlement` | Example 5.11 | pp. 124-125 | Two-span beam: settlement of the interior support | 6 | 0 | 1.1e-03 | 0.83 | active |
| `ex-5-3-gable-frame-stiffness` | Example 5.3 (assembled global stiffness equation) | pp. 102-104 | Gable frame: assembled global stiffness coefficients (rotated members) | 18 | 2 | 9.2e-05 | 1.43 | quarantined |
| `ex-5-4-space-truss` | Example 5.4 | pp. 104-106 | Pin-jointed space truss: four bars to one loaded joint | 15 | 8 | 2.4e-03 | 8.66 | quarantined |
| `ex-5-6-beam-member-loads` | Example 5.6 | pp. 113-114 | Two-span beam with a distributed load and a point load between nodes | 9 | 3 | 3.4e-03 | 5.68 | quarantined |
| `ex-5-7-frame-inclined-member` | Example 5.7 | pp. 115-116 | Two-member rigid frame with an inclined member, joint load and distributed load | 9 | 7 | 5.0e-03 | 63.98 | quarantined |
| `ex-5-8-grid` | Example 5.8 | pp. 116-119 | Grid (plane frame loaded normal to its plane): bending and torsion | 9 | 7 | 1.0e-03 | 7.00 | quarantined |
| `ex-5-9-portal-axially-rigid` | Example 5.9 | pp. 119-120 | Portal frame, axial deformation neglected: vertical and horizontal load | 6 | 0 | 1.6e-03 | 0.81 | active |
| `ex-7-10-cantilever-shear-deformation` | Example 7.10 | p. 206 | End-loaded cantilever of rectangular section with shear deformation | 2 | 0 | 2.5e-04 | 0.50 | active |
| `ex-8-3-beam-column-end-moments` | Example 8.3, Equation f (the member of Example 9.3) | pp. 223-224 (data: p. 249) | Pin-ended beam-column with equal end moments: midheight deflection | 1 | 0 | 5.0e-15 | 0.00 | quarantined |
| `ex-8-4-cantilever-beam-column` | Example 8.4, Equations e and g (the member of Example 9.4) | pp. 225-226 (data: p. 250) | Cantilever with axial and lateral end load: compression and tension | 2 | 0 | 4.7e-15 | 0.00 | quarantined |
| `ex-9-1-truss-critical-load` | Example 9.1, case 1 (alpha = 0), Equation d | p. 247 | Two-bar truss: elastic critical load of a bar restrained by a weak tie | 1 | 0 | 1.1e-10 | 0.00 | active |
| `ex-9-3-column-pinned-critical-load` | Example 9.3, case 1 (alpha = 0): value 'from Example 8.3' | pp. 249-250 | Pin-ended column: elastic critical load | 1 | 0 | 1.5e-04 | 0.33 | active |
| `ex-9-4-column-cantilever-critical-load` | Example 9.4, case 1 (alpha = 0) | p. 250 | Cantilever column: elastic critical load | 1 | 0 | 1.5e-04 | 0.33 | active |
| `ex-9-5-frame-critical-load` | Example 9.5, case 1 (alpha = 0) | p. 251 | Two-storey column with beams: elastic critical load | 1 | 1 | 3.1e-03 | 41.18 | quarantined |
| `ex-9-5-frame-first-order` | Example 9.5, case 2 (alpha = 0.01, P = 3978 kN), moment diagram for Eq. 8.1 | p. 251 | Two-storey column with beams: first-order moments under vertical and lateral load | 3 | 0 | 8.7e-04 | 0.80 | active |
| `ex-9-5-frame-second-order` | Example 9.5, case 2 (alpha = 0.01, P = 3978 kN), moment diagram for Eq. 8.3 | p. 251 | Two-storey column with beams: second-order moments under vertical and lateral load | 3 | 3 | 4.5e-02 | 79.43 | quarantined |
| `ex-9-6-portal-critical-load` | Example 9.6, part 1 | p. 252 | Unsymmetrically loaded pin-based portal frame: critical load factor | 1 | 0 | 1.0e-03 | 0.46 | active |
| `ex-9-6-portal-first-order` | Example 9.6, part 2, moment diagram for Eq. 8.1 | p. 252 | Unsymmetrically loaded pin-based portal frame: first-order moments | 3 | 0 | 8.1e-04 | 0.51 | active |
| `ex-9-6-portal-second-order` | Example 9.6, part 2, moment diagram for Eq. 8.3 | p. 252 | Unsymmetrically loaded pin-based portal frame: second-order moments | 3 | 3 | 6.1e-02 | 63.23 | quarantined |
| `ex-9-8-column-braced-critical-load` | Example 9.8, part 1, sketch (a): plateau value; part 2 brace area | p. 254 | Column fully braced at midheight by an elastic strut: critical load | 1 | 1 | 7.6e-03 | 32.95 | quarantined |
| `ex-9-8-column-unbraced-critical-load` | Example 9.8, part 1, sketch (a): value at A_br = 0 | p. 254 | Column with a midheight strut of zero area: critical load | 1 | 0 | 2.4e-04 | 0.26 | active |

## Examples not encoded

**Not applicable: engine limit**

| example | reason |
|---|---|
| 2.4, 5.12, 5.13 | thermal loads: not in the input format |
| 2.5 | rigid-body kinematics of a bar on an inclined surface: no inclined supports, no result of this kind |
| 4.1 | strain energy (367.0 J): not a quantity of the results document |
| 4.15 part 1 | continuous elastic foundation (closed form 0.237 in, 554 in-kips): no foundation member |
| 5.5 | direction cosines of a rolled member: not a quantity of the results document |
| 7.2, 7.3, 7.9 | tapered member: prismatic members only |
| 7.4, 7.5, 7.6 | warping (nonuniform) torsion: the engine has Saint-Venant torsion only |
| 7.12 | circular arch by a curved element, axial deformation neglected: straight members only; a polygon of straight members would carry a discretisation error the book does not bound |
| 8.1, 8.2, 9.1 case 2, 9.2 | large-displacement response, limit points, snap-through: the engine's second-order analysis covers the effect of axial force on flexure (P-Delta, P-delta), not large displacements |
| 8.5, 8.6, 12.5, chapter 10 | material nonlinearity |
| 8.7 | nonlinear rotational spring |
| 9.9 | torsional buckling: the engine's geometric stiffness covers axial-force effects on flexure only |
| 9.10, 9.11, A.2 | lateral-torsional buckling: engine limit |
| 9.12, A.1, A.3 | torsional-flexural buckling: engine limit |
| 13.4 | general linear constraint equations (v2 = v3, theta2 = -theta3): the engine has rigid links and diaphragms only |

**Excluded: input incomplete or ambiguous in the book** (questions recorded, nothing guessed)

| example | reason |
|---|---|
| 9.7 | support type and member layout are given only by a small figure; the printed critical load (210 lb) is also a discretised value |
| 9.13 | section properties of the columns, the beam and the area of the diagonals are not printed (shape names only) |
| 9.14 | joint coordinates are not printed (isometric sketch) |
| 13.9 | section properties of the W-shapes are not printed |

**No result to compare, or covered otherwise**

| example | reason |
|---|---|
| 2.6, 4.2-4.7, 5.1, 5.2, 7.7, 7.8, 8.8, 8.9, 13.1, 13.6-13.8 | symbolic derivations (13.6-13.8 are covered by the numeric instances `eq-13-31-...` and `eq-13-38-...`; 4.14 by `ex-4-14-...`) |
| 3.1, 3.4, 4.8 | stiffness matrices of unsupported structures / element matrices; the same structures are solved in 3.2, 3.5, 4.9-4.12 |
| 7.10, 7.11 | 7.10 is encoded; 7.11 evaluates a shear-area formula (no analysis result) |
| 9.3 case 2, 9.4 cases 2 and 3, 9.5 curve, 9.6 curve, 9.8 part 2, 12.2, 12.3 | results shown only as curves; the closed forms behind 9.3 and 9.4 are covered by `ex-8-3-...` and `ex-8-4-...` |
| 12.1, 12.4, 12.6-12.11 (matrix iterations), 13.13 | illustrations of numerical procedures (the structural content of 12.6-12.9 is `ex-12-7-...`) |
| chapter 6, chapter 11 | chapter 6 (virtual work) was not reviewed for this library; chapter 11 treats equation solving |

## Open questions (for the project owner)

1. **Policy against hand-calculated textbook numbers.** The policy counts every printed digit of a target as
   significant. 18 problems of this book are quarantined only because its hand calculations carry three- or
   four-digit intermediate values, so that the last printed digit of a result is not reliable (documented case
   by case with equilibrium residuals and with the exact solution of the book's own printed equations). Whether a
   written, source-level rule for such sources is wanted (for example a tolerance from the number of digits the
   source carries in its working) is not for the problem author to decide; no tolerance was changed here.
2. **Example 4.15**: the support area is printed as "22.5 x 15 / 29,000 = 0.0116 in^2". With the expression the
   printed results are reproduced; with 0.0116 they are missed in the last digit. Which did the book's run use?
3. **Examples 9.5 and 9.6, second-order moments**: the book's incremental scheme (single-step or iterative, step
   size, geometry update) is not stated; its values are 3 % to 6 % below the equilibrium solution of its own
   element matrices.
4. **Example 9.6**: the label 1660 on the moment diagram is not identified by the book. (A diagnosis calculation
   with the book's discretisation gives 1668 in-kips for the second-order moment at midheight of column de.) Not
   used as a target.
5. **Example 9.3**: alpha is printed as 1.25e-4 in the statement and 1.125e-4 in the solution heading (it affects
   only the curve, which is not a target; `ex-8-3-...` uses 1.25e-4 as an instance of the closed form).
6. **Example 9.7**: are the two wall supports pinned (rotation free)? Excluded until answered.
7. **Example 13.3**: the figure shows the joint moments as curved arrows, the load vector as +2 x 10^3 kN mm;
   the load vector was followed (counter-clockwise).
8. **Examples 8.3 and 8.4** fail tier C for the reason already recorded for the tier-A beam-column problems
   (the estimated OpenSees discretisation error is not a bound).
