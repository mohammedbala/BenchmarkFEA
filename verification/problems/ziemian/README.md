# Published verification problems: Ziemian & Ziemian steel benchmark frames

## Attribution and licence

The frames in this directory are taken from

> C. W. Ziemian, R. D. Ziemian, "Steel benchmark frames for structural analysis and validation studies: Finite
> element models and numerical simulation data", *Data in Brief* 39 (2021) 107564,
> doi:[10.1016/j.dib.2021.107564](https://doi.org/10.1016/j.dib.2021.107564);
> data set: Mendeley Data, doi:[10.17632/39sjhchwtx.1](https://doi.org/10.17632/39sjhchwtx.1).

Licence of the article and of the data set: **Creative Commons Attribution 4.0 International (CC BY 4.0)**,
https://creativecommons.org/licenses/by/4.0/. The model data of six frames (node coordinates, connectivity, end
connections, section properties, loads) and the published elastic results used as targets are reproduced here
under that licence, **with changes**: converted to this project's input format; the data set's 4 elements per
member merged into one member with the quarter-point loads as point loads on the member; the plane model embedded
in 3D (out-of-plane DOFs held, Iy = 1000 Izz instead of the listed Iyy); 0.8 E entered as the modulus; results
entered as displayed by the data set's spreadsheet. The authors of the data set are not connected with this
project and have not reviewed this use. The same attribution is in every `problem.json` (`source.notes`).

The two spreadsheets of the data set are **not** in this repository. `build_from_dataset.py <directory>` reads
them ("FEA Model Details_All frames.xlsx" saved as `model.xlsx`, "Simulation Results_All frames.xlsx" as
`results.xlsx`) with `dataset.py` and writes the problem directories.

## What is used

Frames 1 (portal, fixed bases), 2 (gabled frame), 3 (two-storey frame with pin-ended beams and diagonals),
5 (two-bay frame with unequal column heights), 7 (braced two-bay frame with a lean-on column) and 9 (two-storey
two-bay frame under gravity loads only). Per frame, the **elastic** results only:

* first-order analysis ("LA"): lateral joint displacements, cases "Right" and "Left"; member design moments for
  frames 1, 2 and 9;
* second-order elastic analysis ("WC"): the same quantities;
* elastic critical load ratio ("LBA").

Not used: the approximate "SIPC" results and everything inelastic.

## Analysis conditions of the data set, and how each was verified

Article, Tables 2 and 3; the program is MASTAN2.

* **0.8 E** for all elastic analyses. The model file's element table points to its material 1 (E = 29,000 ksi);
  the published results are obtained with 0.8 E = 23,200 ksi. Verified for **each** frame used:

| frame | case | largest rel. difference of the first-order joint displacements with 0.8 E | same with E = 29,000 ksi |
|---|---|---|---|
| 1 | right / left | 1.7e-06 / 1.7e-06 | 0.20 |
| 2 | right / left | 4.3e-06 / 4.2e-06 | 0.20 |
| 3 | right / left | 1.4e-03 / 1.4e-03 | 0.20 |
| 5 | right / left | 1.3e-05 / 9.9e-06 | 0.20 |
| 7 | right / left | 1.9e-06 / 4.2e-07 | 0.20 |
| 9 | right / left | 9.5e-06 / 7.2e-06 | 0.20 |

  (frame 3: see "open questions").
* **No shear deformation** (the file's shear areas are "Inf").
* **Geometry**: LA and WC on the imperfect geometry, every node shifted horizontally by (its height)/500; LBA on
  the perfect geometry.
* **The cases "Right" and "Left"**: "Right" = nodes at x + y/500 with the lateral loads acting toward +x; "Left" =
  nodes at x - y/500 with the lateral loads acting toward -x; gravity loads unchanged. The model file lists the
  imperfect coordinates with the rightward lean for every frame and the lateral loads in the direction of the
  description file; for frames 5 and 7 that direction is leftward, so the file's own combination is neither
  published case. This definition reproduces the first-order results of both cases of every frame (table above);
  the other combinations differ by 2 % to 300 %.
* **Discretisation**: 4 elements per member, distributed loads lumped at the quarter points. The engine's member is
  exact for point loads, so one member per straight member line is the same structure.
* **Second-order analysis of the data set**: work-control incremental analysis, step 0.001, updated Lagrangian
  formulation with geometric stiffness matrices - a large-displacement analysis with a polynomial element.
* **Listing precision**: coordinates and loads are listed with five significant figures (the rafter quarter points
  of frame 2 at y = 339.88, 379.75, 419.62 for 339.875, 379.75, 419.625). Interior nodes within that rounding of
  the straight line between the member ends are taken as on the line.

## Printed precision

The results spreadsheet stores doubles and displays them with a fixed number of decimals: displacements "0.000",
moments "0.00", critical load ratios "0.000". The target is the number **as displayed**, the tolerance half a unit
of its last displayed digit; the stored value is recorded in each `precisionNote`. Entries stored as 0 (supports,
and values the data set filtered out) are not targets.

## Member numbering

The data set numbers its members without defining the numbering. `build_from_dataset.py --infer-members` identifies
a published member with the straight member line (a) whose Zzz x Fy equals the tab's "plastic section modulus x
Fy" and (b) whose largest first-order moment is closest to the published LA values of both cases, accepted only
if unique and closer than 1e-4. Frames 1, 2 and 9 are identified this way (the mapping is stored in `MEMBERS`);
frames 3, 5 and 7 list fewer members than they have straight member lines and get no moment targets. The
first-order moment targets are therefore **not independent evidence of the identification**; the second-order ones
use an identification fixed beforehand.

## Result: 36 problems, 14 active, 22 quarantined

* **First-order displacements** (12 problems): 9 active. Quarantined: frame 3 right and left (differences of up to
  1.4e-3, cause not identified, see below) and frame 9 right (one value on a rounding boundary of the displayed
  number).
* **First-order member moments** (6 problems): 1 active (frame 1 left). Five quarantined: the differences are
  1e-6 to 2.6e-5 of the moments, the size of the rounding of the listed loads (five significant figures), while the
  moments are displayed to 0.01 kip in (1e-6 to 1e-5 relative).
* **Critical load ratios** (6 problems): 4 active (frames 2, 3, 5, 9; differences 6.5e-6 to 1.5e-4 against a
  tolerance of 5e-4 absolute). Frames 1 and 7 quarantined: the engine is lower by 9.5e-5 and 6.1e-5 of the value,
  which crosses the third displayed decimal. `verification/derive/ziemian_diagnosis.py` reproduces the data set's
  stored values to 1e-6 ... 5e-6 with the data set's own discretisation (cubic elements with a consistent geometric
  stiffness), which overestimates the critical load; the engine's exact stability functions give the lower value.
* **Second-order analysis** (12 problems): **all quarantined**. The data set's WC analysis is a large-displacement
  analysis; the engine solves the classical small-displacement second-order theory. Measured (case "Right", joint
  displacements, `verification/derive/ziemian_diagnosis.py`):

  | frame | engine vs OpenSeesPy "PDelta", 32 el./member | OpenSeesPy "Corotational", 32 el./member, vs data set | engine vs data set |
  |---|---|---|---|
  | 1 | 4.8e-05 | 1.4e-03 | 3.5e-03 |
  | 2 | 3.0e-05 | 1.5e-03 | 1.8e-02 |
  | 3 | 3.0e-04 | 6.6e-03 | 1.2e-02 |
  | 5 | 3.8e-05 | 2.1e-03 | 1.6e-02 |
  | 7 | (engine refuses the case, see below) | 1.1e-03 | - |
  | 9 | 1.2e-04 | 4.3e-03 | 2.7e-02 |

  The engine reproduces the linearised (P-Delta) theory to 1e-4; a corotational analysis agrees with the data set
  to 1e-3 ... 7e-3; the two theories differ by 0.4 % to 2.7 % in the joint displacements of these frames (largest
  where the sway is small compared with the beam deflections: frame 9, gravity only). The remaining difference
  between the corotational results and the data set was not resolved. **This is a limit of the engine's theory
  relative to this data set, not a defect found in either.**
  Frame 7, case "Right": the engine refuses the analysis - the single diagonal brace is in compression beyond its
  Euler load (the description file states that only the tension brace is modelled, i.e. the case "Left").

Diagnoses with the measured numbers: `quarantine.py`.

## How to read the index

`targets` / `missed` and the differences are those of `../check_local_v2b.py` (engine against the targets of
`problem.json`); "largest difference / tolerance" above 1 means a target is missed. `status` is the declared
status. "Verified" in this project means only: passed the tiers listed in `docs/VERIFICATION.md` at the
tolerances stated there. The other tiers (cross-check with foreign solvers, invariants, sensitivity) are run by
`benchfea verify verification/problems/ziemian/<id>`.

Tolerance policy (docs/VERIFICATION.md, section 5), applied uniformly and checked by `../check_local_v2b.py`:
a published number carries half a unit in its last printed digit (floor 1e-6) - here: of the number as displayed by the spreadsheet; a number computed here
carries 1e-6. A problem that misses, or on whose model another tier does not pass, is quarantined with a written
diagnosis (`quarantine.py`); nothing is loosened or dropped.

## Index

| problem | source location | what it verifies | basis | targets | missed | largest rel. difference | largest difference / tolerance | status |
|---|---|---|---|---|---|---|---|---|
| `frame-1-critical-load-ratio` | Frame 1: 'FEA Model Details_All frames.xlsx', tab 'Frame 1'; 'Simulation Results_All frames.xlsx', tab 'SUMMARY', column acr | Ziemian frame 1 (single-storey portal frame, fixed bases): elastic critical load ratio | published | 1 | 1 | 9.4e-05 | 1.48 | quarantined |
| `frame-1-first-order-displacements-left` | Frame 1: 'FEA Model Details_All frames.xlsx', tab 'Frame 1'; 'Simulation Results_All frames.xlsx', tab 'Frame 1', columns LA / Left | Ziemian frame 1 (single-storey portal frame, fixed bases): first-order joint displacements, case 'left' | published | 2 | 0 | 4.3e-05 | 0.09 | active |
| `frame-1-first-order-displacements-right` | Frame 1: 'FEA Model Details_All frames.xlsx', tab 'Frame 1'; 'Simulation Results_All frames.xlsx', tab 'Frame 1', columns LA / Right | Ziemian frame 1 (single-storey portal frame, fixed bases): first-order joint displacements, case 'right' | published | 2 | 0 | 2.2e-05 | 0.04 | active |
| `frame-1-first-order-moments-left` | Frame 1: 'FEA Model Details_All frames.xlsx', tab 'Frame 1'; 'Simulation Results_All frames.xlsx', tab 'Frame 1', columns LA / Left | Ziemian frame 1 (single-storey portal frame, fixed bases): first-order member design moments, case 'left' | published | 3 | 0 | 2.9e-06 | 0.95 | active |
| `frame-1-first-order-moments-right` | Frame 1: 'FEA Model Details_All frames.xlsx', tab 'Frame 1'; 'Simulation Results_All frames.xlsx', tab 'Frame 1', columns LA / Right | Ziemian frame 1 (single-storey portal frame, fixed bases): first-order member design moments, case 'right' | published | 3 | 1 | 2.9e-06 | 1.24 | quarantined |
| `frame-1-second-order-left` | Frame 1: 'FEA Model Details_All frames.xlsx', tab 'Frame 1'; 'Simulation Results_All frames.xlsx', tab 'Frame 1', columns WC / Left | Ziemian frame 1 (single-storey portal frame, fixed bases): second-order joint displacements and member design moments, case 'left' | published | 5 | 5 | 2.7e-03 | 422.19 | quarantined |
| `frame-1-second-order-right` | Frame 1: 'FEA Model Details_All frames.xlsx', tab 'Frame 1'; 'Simulation Results_All frames.xlsx', tab 'Frame 1', columns WC / Right | Ziemian frame 1 (single-storey portal frame, fixed bases): second-order joint displacements and member design moments, case 'right' | published | 5 | 5 | 3.3e-03 | 355.01 | quarantined |
| `frame-2-critical-load-ratio` | Frame 2: 'FEA Model Details_All frames.xlsx', tab 'Frame 2'; 'Simulation Results_All frames.xlsx', tab 'SUMMARY', column acr | Ziemian frame 2 (single-storey gabled frame, pinned bases): elastic critical load ratio | published | 1 | 0 | 3.8e-05 | 0.45 | active |
| `frame-2-first-order-displacements-left` | Frame 2: 'FEA Model Details_All frames.xlsx', tab 'Frame 2'; 'Simulation Results_All frames.xlsx', tab 'Frame 2', columns LA / Left | Ziemian frame 2 (single-storey gabled frame, pinned bases): first-order joint displacements, case 'left' | published | 3 | 0 | 5.7e-05 | 0.59 | active |
| `frame-2-first-order-displacements-right` | Frame 2: 'FEA Model Details_All frames.xlsx', tab 'Frame 2'; 'Simulation Results_All frames.xlsx', tab 'Frame 2', columns LA / Right | Ziemian frame 2 (single-storey gabled frame, pinned bases): first-order joint displacements, case 'right' | published | 3 | 0 | 7.5e-05 | 0.82 | active |
| `frame-2-first-order-moments-left` | Frame 2: 'FEA Model Details_All frames.xlsx', tab 'Frame 2'; 'Simulation Results_All frames.xlsx', tab 'Frame 2', columns LA / Left | Ziemian frame 2 (single-storey gabled frame, pinned bases): first-order member design moments, case 'left' | published | 4 | 2 | 9.5e-06 | 1.96 | quarantined |
| `frame-2-first-order-moments-right` | Frame 2: 'FEA Model Details_All frames.xlsx', tab 'Frame 2'; 'Simulation Results_All frames.xlsx', tab 'Frame 2', columns LA / Right | Ziemian frame 2 (single-storey gabled frame, pinned bases): first-order member design moments, case 'right' | published | 4 | 2 | 1.5e-05 | 2.39 | quarantined |
| `frame-2-second-order-left` | Frame 2: 'FEA Model Details_All frames.xlsx', tab 'Frame 2'; 'Simulation Results_All frames.xlsx', tab 'Frame 2', columns WC / Left | Ziemian frame 2 (single-storey gabled frame, pinned bases): second-order joint displacements and member design moments, case 'left' | published | 7 | 7 | 2.1e-02 | 2258.88 | quarantined |
| `frame-2-second-order-right` | Frame 2: 'FEA Model Details_All frames.xlsx', tab 'Frame 2'; 'Simulation Results_All frames.xlsx', tab 'Frame 2', columns WC / Right | Ziemian frame 2 (single-storey gabled frame, pinned bases): second-order joint displacements and member design moments, case 'right' | published | 7 | 7 | 1.8e-02 | 1891.03 | quarantined |
| `frame-3-critical-load-ratio` | Frame 3: 'FEA Model Details_All frames.xlsx', tab 'Frame 3'; 'Simulation Results_All frames.xlsx', tab 'SUMMARY', column acr | Ziemian frame 3 (two-storey one-bay frame with pin-ended beams and diagonals): elastic critical load ratio | published | 1 | 0 | 6.5e-06 | 0.02 | active |
| `frame-3-first-order-displacements-left` | Frame 3: 'FEA Model Details_All frames.xlsx', tab 'Frame 3'; 'Simulation Results_All frames.xlsx', tab 'Frame 3', columns LA / Left | Ziemian frame 3 (two-storey one-bay frame with pin-ended beams and diagonals): first-order joint displacements, case 'left' | published | 4 | 1 | 2.1e-03 | 2.39 | quarantined |
| `frame-3-first-order-displacements-right` | Frame 3: 'FEA Model Details_All frames.xlsx', tab 'Frame 3'; 'Simulation Results_All frames.xlsx', tab 'Frame 3', columns LA / Right | Ziemian frame 3 (two-storey one-bay frame with pin-ended beams and diagonals): first-order joint displacements, case 'right' | published | 4 | 1 | 1.8e-03 | 2.10 | quarantined |
| `frame-3-second-order-left` | Frame 3: 'FEA Model Details_All frames.xlsx', tab 'Frame 3'; 'Simulation Results_All frames.xlsx', tab 'Frame 3', columns WC / Left | Ziemian frame 3 (two-storey one-bay frame with pin-ended beams and diagonals): second-order joint displacements, case 'left' | published | 4 | 4 | 1.2e-02 | 44.06 | quarantined |
| `frame-3-second-order-right` | Frame 3: 'FEA Model Details_All frames.xlsx', tab 'Frame 3'; 'Simulation Results_All frames.xlsx', tab 'Frame 3', columns WC / Right | Ziemian frame 3 (two-storey one-bay frame with pin-ended beams and diagonals): second-order joint displacements, case 'right' | published | 4 | 4 | 1.2e-02 | 44.52 | quarantined |
| `frame-5-critical-load-ratio` | Frame 5: 'FEA Model Details_All frames.xlsx', tab 'Frame 5'; 'Simulation Results_All frames.xlsx', tab 'SUMMARY', column acr | Ziemian frame 5 (single-storey two-bay frame with unequal column heights): elastic critical load ratio | published | 1 | 0 | 1.5e-04 | 0.99 | active |
| `frame-5-first-order-displacements-left` | Frame 5: 'FEA Model Details_All frames.xlsx', tab 'Frame 5'; 'Simulation Results_All frames.xlsx', tab 'Frame 5', columns LA / Left | Ziemian frame 5 (single-storey two-bay frame with unequal column heights): first-order joint displacements, case 'left' | published | 4 | 0 | 4.9e-04 | 0.89 | active |
| `frame-5-first-order-displacements-right` | Frame 5: 'FEA Model Details_All frames.xlsx', tab 'Frame 5'; 'Simulation Results_All frames.xlsx', tab 'Frame 5', columns LA / Right | Ziemian frame 5 (single-storey two-bay frame with unequal column heights): first-order joint displacements, case 'right' | published | 4 | 0 | 2.8e-04 | 0.95 | active |
| `frame-5-second-order-left` | Frame 5: 'FEA Model Details_All frames.xlsx', tab 'Frame 5'; 'Simulation Results_All frames.xlsx', tab 'Frame 5', columns WC / Left | Ziemian frame 5 (single-storey two-bay frame with unequal column heights): second-order joint displacements, case 'left' | published | 4 | 4 | 1.3e-02 | 53.98 | quarantined |
| `frame-5-second-order-right` | Frame 5: 'FEA Model Details_All frames.xlsx', tab 'Frame 5'; 'Simulation Results_All frames.xlsx', tab 'Frame 5', columns WC / Right | Ziemian frame 5 (single-storey two-bay frame with unequal column heights): second-order joint displacements, case 'right' | published | 4 | 4 | 1.6e-02 | 72.77 | quarantined |
| `frame-7-critical-load-ratio` | Frame 7: 'FEA Model Details_All frames.xlsx', tab 'Frame 7'; 'Simulation Results_All frames.xlsx', tab 'SUMMARY', column acr | Ziemian frame 7 (single-storey two-bay frame with a diagonal brace and a lean-on column): elastic critical load ratio | published | 1 | 1 | 1.2e-04 | 1.56 | quarantined |
| `frame-7-first-order-displacements-left` | Frame 7: 'FEA Model Details_All frames.xlsx', tab 'Frame 7'; 'Simulation Results_All frames.xlsx', tab 'Frame 7', columns LA / Left | Ziemian frame 7 (single-storey two-bay frame with a diagonal brace and a lean-on column): first-order joint displacements, case 'left' | published | 4 | 0 | 1.8e-03 | 0.67 | active |
| `frame-7-first-order-displacements-right` | Frame 7: 'FEA Model Details_All frames.xlsx', tab 'Frame 7'; 'Simulation Results_All frames.xlsx', tab 'Frame 7', columns LA / Right | Ziemian frame 7 (single-storey two-bay frame with a diagonal brace and a lean-on column): first-order joint displacements, case 'right' | published | 4 | 0 | 4.5e-03 | 0.93 | active |
| `frame-7-second-order-left` | Frame 7: 'FEA Model Details_All frames.xlsx', tab 'Frame 7'; 'Simulation Results_All frames.xlsx', tab 'Frame 7', columns WC / Left | Ziemian frame 7 (single-storey two-bay frame with a diagonal brace and a lean-on column): second-order joint displacements, case 'left' | published | 4 | 3 | 3.0e-02 | 11.80 | quarantined |
| `frame-7-second-order-right` | Frame 7: 'FEA Model Details_All frames.xlsx', tab 'Frame 7'; 'Simulation Results_All frames.xlsx', tab 'Frame 7', columns WC / Right | Ziemian frame 7 (single-storey two-bay frame with a diagonal brace and a lean-on column): second-order joint displacements, case 'right' | published | 2 | - | - | - | quarantined (the engine refuses the model) |
| `frame-9-critical-load-ratio` | Frame 9: 'FEA Model Details_All frames.xlsx', tab 'Frame 9'; 'Simulation Results_All frames.xlsx', tab 'SUMMARY', column acr | Ziemian frame 9 (two-storey two-bay frame, gravity loads only): elastic critical load ratio | published | 1 | 0 | 3.3e-05 | 0.17 | active |
| `frame-9-first-order-displacements-left` | Frame 9: 'FEA Model Details_All frames.xlsx', tab 'Frame 9'; 'Simulation Results_All frames.xlsx', tab 'Frame 9', columns LA / Left | Ziemian frame 9 (two-storey two-bay frame, gravity loads only): first-order joint displacements, case 'left' | published | 6 | 0 | 5.8e-04 | 0.96 | active |
| `frame-9-first-order-displacements-right` | Frame 9: 'FEA Model Details_All frames.xlsx', tab 'Frame 9'; 'Simulation Results_All frames.xlsx', tab 'Frame 9', columns LA / Right | Ziemian frame 9 (two-storey two-bay frame, gravity loads only): first-order joint displacements, case 'right' | published | 6 | 1 | 2.1e-03 | 1.00 | quarantined |
| `frame-9-first-order-moments-left` | Frame 9: 'FEA Model Details_All frames.xlsx', tab 'Frame 9'; 'Simulation Results_All frames.xlsx', tab 'Frame 9', columns LA / Left | Ziemian frame 9 (two-storey two-bay frame, gravity loads only): first-order member design moments, case 'left' | published | 10 | 8 | 4.3e-05 | 15.20 | quarantined |
| `frame-9-first-order-moments-right` | Frame 9: 'FEA Model Details_All frames.xlsx', tab 'Frame 9'; 'Simulation Results_All frames.xlsx', tab 'Frame 9', columns LA / Right | Ziemian frame 9 (two-storey two-bay frame, gravity loads only): first-order member design moments, case 'right' | published | 10 | 7 | 1.5e-04 | 15.49 | quarantined |
| `frame-9-second-order-left` | Frame 9: 'FEA Model Details_All frames.xlsx', tab 'Frame 9'; 'Simulation Results_All frames.xlsx', tab 'Frame 9', columns WC / Left | Ziemian frame 9 (two-storey two-bay frame, gravity loads only): second-order joint displacements and member design moments, case 'left' | published | 16 | 16 | 1.4e-02 | 1515.42 | quarantined |
| `frame-9-second-order-right` | Frame 9: 'FEA Model Details_All frames.xlsx', tab 'Frame 9'; 'Simulation Results_All frames.xlsx', tab 'Frame 9', columns WC / Right | Ziemian frame 9 (two-storey two-bay frame, gravity loads only): second-order joint displacements and member design moments, case 'right' | published | 16 | 16 | 2.8e-02 | 1396.97 | quarantined |

## Open questions

1. Frame 3: the first-order displacements differ from the data set by up to 1.4e-3 with alternating sign between
   the two column lines, although the critical load ratio of the same model is reproduced to 6.5e-6. Are the loads
   of the published LA / WC analyses of frame 3 the ones listed (235.25 and 239.88 kip, 1.544 kip twice)?
2. How are the members of the results file numbered (definition, and for frames 3, 5, 7 which elements form a
   "member")?
3. The first-order results of frame 1 are all smaller than the engine's by the same factor 1 - 1.73e-6: what are
   the unrounded loads?
4. Is there a tolerance the authors consider appropriate for comparing another program's second-order results
   with the WC values (the article reports the SIPC errors against WC, not an acceptance criterion)?
