# Published verification problems: AISC Design Examples v16.0, Examples C.1A, C.1B, C.1C (frame level)

Source: American Institute of Steel Construction, *Companion to the AISC Steel Construction Manual, Volume 1:
Design Examples*, Version 16.0, Examples C.1A (direct analysis method, pp. C-2 to C-6), C.1B (effective length
method, pp. C-7 to C-12) and C.1C (first-order analysis method, pp. C-13 to C-17). Nothing of the source is
reproduced here; the input data and the printed results are given with their location. The problems are written
by `build_problems.py`.

## The example

One-storey frame: two W12x65 columns (20 ft, pinned bases) and a W18x40 beam (30 ft, rigidly connected) stabilise
the gravity load of three further "leaning" columns, represented in the examples' analysis model by one pin-ended
leaning column linked to the frame (C.1A, C.1B; C.1C has no leaning column in the model). Loads: the maximum
gravity combination and a notional (C.1A, C.1B: 0.002 alpha Y_i) or additional (C.1C: 0.0042 Y_i) lateral load.
The examples print, in figures, the results of an analysis **by software**: drift, column moments at the beam,
base shears and vertical reactions; first order and second order; LRFD and ASD (ASD: analysed at 1.6 x the loads,
forces printed divided by 1.6, drift undivided).

Ten problems: C.1A first and second order (0.8 x stiffness), C.1B first and second order (full stiffness), C.1C
first order; each LRFD and ASD; 7 targets each.

## Status: all ten problems are QUARANTINED - the example is not usable at its printed precision

The source does not completely specify its analysis model:

| item | given by the source? |
|---|---|
| geometry, loads, notional loads, stiffness factor | yes |
| W12x65: Ag = 19.1 in^2, Ix = 533 in^4; W18x40: Ix = 612 in^4 | yes (C.1A-C) |
| E = 29,000 ksi | yes (used in the examples' evaluations) |
| area of the W18x40; properties of the leaning column and its link | no |
| whether the software includes axial and shear deformation; shear areas | no |
| the software | no ("analysis software") |

The model that can be built from the given data (columns with their Ag; beam, link and leaning column axially
rigid; **no shear deformation**) does not reproduce the printed results at their printed precision (three
significant figures): the first-order drift is 1.8 % to 2.2 % smaller, base shears and column moments are 0.4 %
to 1.0 % larger; the vertical reactions (statics) are reproduced.

Parameter study (Example C.1B, first order, LRFD; printed: drift 0.145 in, base shears 5.64 and 6.21 kips):

| model | drift (in) | base shears (kips) |
|---|---|---|
| as built here (no shear deformation, axially rigid beam) | 0.14216 | 5.671, 6.247 |
| beam area 5 in^2 instead of rigid (average of the drifts at B and C) | 0.14216 | 5.662, 6.238 |
| beam area 20 in^2 | 0.14216 | 5.669, 6.245 |
| with shear deformation, shear areas 4 in^2 (column) and 6 in^2 (beam) | 0.1451 | 5.635, 6.211 |
| with shear deformation, shear areas 5 in^2 (column) and 5 in^2 (beam) | 0.1451 | 5.642, 6.218 |

So the beam area, which is not given, does not matter (the average drift is unchanged to 1e-5 in); **shear
deformation with shear areas of the size of the web areas reproduces all three printed numbers**. The shear areas
in the last two rows are trial values of this study: they are not taken from any table, not entered in the models
and not used as evidence of anything except the diagnosis. Since no opened reference supplies the shear areas the
example's software used, the example cannot be completed without a guess; by the project rule ("where a reference
is silent, do not guess: exclude the affected case") the problems are kept, evaluated and reported, and never
counted as passed. Diagnosis: `quarantine.py`.

Remark (ASD, C.1A / C.1B): the leaning-column load is entered as printed in the figure, 154 kips (60.0 ft x 2.56
kip/ft is 153.6 kips).

## How to read the index

`targets` / `missed` and the differences are those of `../check_local_v2b.py` (engine against the targets of
`problem.json`); "largest difference / tolerance" above 1 means a target is missed. `status` is the declared
status. "Verified" in this project means only: passed the tiers listed in `docs/VERIFICATION.md` at the
tolerances stated there. The other tiers (cross-check with foreign solvers, invariants, sensitivity) are run by
`benchfea verify verification/problems/aisc-de/<id>`.

Tolerance policy (docs/VERIFICATION.md, section 5), applied uniformly and checked by `../check_local_v2b.py`:
a published number carries half a unit in its last printed digit (floor 1e-6); a number computed here
carries 1e-6. A problem that misses, or on whose model another tier does not pass, is quarantined with a written
diagnosis (`quarantine.py`); nothing is loosened or dropped.

## Index

| problem | source location | what it verifies | basis | targets | missed | largest rel. difference | largest difference / tolerance | status |
|---|---|---|---|---|---|---|---|---|
| `c1a-first-order-asd` | Example C.1A, Fig. C.1A-3, ASD column | AISC Design Example C1.A: first-order frame results, ASD | published | 7 | 5 | 1.9e-02 | 9.82 | quarantined |
| `c1a-first-order-lrfd` | Example C.1A, Fig. C.1A-3, LRFD column | AISC Design Example C1.A: first-order frame results, LRFD | published | 7 | 4 | 1.8e-02 | 7.40 | quarantined |
| `c1a-second-order-asd` | Example C.1A, Fig. C.1A-4, ASD column | AISC Design Example C1.A: second-order frame results, ASD | published | 7 | 5 | 7.4e-03 | 12.54 | quarantined |
| `c1a-second-order-lrfd` | Example C.1A, Fig. C.1A-4, LRFD column | AISC Design Example C1.A: second-order frame results, LRFD | published | 7 | 4 | 6.7e-03 | 8.38 | quarantined |
| `c1b-first-order-asd` | Example C.1B, Fig. C.1B-3, ASD column | AISC Design Example C1.B: first-order frame results, ASD | published | 7 | 5 | 2.2e-02 | 9.82 | quarantined |
| `c1b-first-order-lrfd` | Example C.1B, Fig. C.1B-3, LRFD column | AISC Design Example C1.B: first-order frame results, LRFD | published | 7 | 4 | 2.0e-02 | 7.40 | quarantined |
| `c1b-second-order-asd` | Example C.1B, Fig. C-1B-4 (as captioned), ASD column | AISC Design Example C1.B: second-order frame results, ASD | published | 7 | 5 | 1.4e-02 | 12.99 | quarantined |
| `c1b-second-order-lrfd` | Example C.1B, Fig. C-1B-4 (as captioned), LRFD column | AISC Design Example C1.B: second-order frame results, LRFD | published | 7 | 4 | 8.0e-03 | 8.71 | quarantined |
| `c1c-first-order-asd` | Example C.1C, Fig. C.1C-3, ASD column | AISC Design Example C1.C: first-order frame results, ASD | published | 7 | 5 | 2.0e-02 | 9.87 | quarantined |
| `c1c-first-order-lrfd` | Example C.1C, Fig. C.1C-3, LRFD column | AISC Design Example C1.C: first-order frame results, LRFD | published | 7 | 4 | 1.8e-02 | 10.73 | quarantined |

## Open questions

1. Which shear areas (and area of the W18x40) did the software of Examples C.1A-C use, and does it include shear
   deformation? With an opened section-property source (AISC Manual Part 1 / Shapes Database) the problems could be
   rebuilt and re-evaluated.
2. Leaning column and link of the examples' analysis model: axially rigid, or with which properties?
