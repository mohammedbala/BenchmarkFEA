# Published verification problems: AISC 360-22 Commentary C2.1 benchmark problems

Source: ANSI/AISC 360-22, *Specification for Structural Steel Buildings* (American Institute of Steel
Construction, August 1, 2022), Commentary to Section C2.1, benchmark problems for second-order analysis:
Case 1 (Fig. C-C2.2, simply supported beam-column, P-delta only) and Case 2 (Fig. C-C2.3, fixed-base cantilever,
P-Delta and P-delta), p. 16.1-337. Nothing of the source is reproduced here: the data that define the two problems
and the tabulated target numbers are given with their location, in this project's own words. The problems are
written by `build_problems.py`.

## What the source gives, and what it does not

* Given: W14x48 (W360x72), major-axis bending, E = 29,000 ksi (200 000 MPa), length 28.0 ft (8.53 m); Case 1:
  0.200 kip/ft (2.92 kN/m) transverse load; Case 2: 1.0 kip (4.45 kN) at the top; four axial load levels per case;
  moment and deflection for each level, with shear deformation and (in brackets) without, in a US customary and
  an SI table, each printed to three significant figures.
* Acceptance criterion **stated by the source** (same page): moments within 3 %, deflections within 5 %. This is
  the tolerance of every target (policy rule 2: a tolerance stated by the source itself).
* NOT given: I_x, A, the shear area of the W14x48, and G. **No section property is taken from memory or from a
  table that was not opened.** Exactly as in `tests/test_second_order_aisc_benchmarks.py`:
  * I_x is *calibrated* on the source's own first-order (P = 0) deflection without shear deformation of Case 2,
    I = H L^3 / (3 E Delta) = 483.9 in^4 (US table) / 4.021e8 mm^4 (SI table). That one number is therefore
    reproduced by construction and is not evidence (marked in its `precisionNote`).
  * the shear flexibility 1/(G A_s) is *calibrated* per case and per table on the difference between the P = 0
    deflections with and without shear deformation, which the source prints to ONE significant digit. The
    with-shear problems are a consistency test of the amplification, not a test of a shear area; their P = 0
    deflection is reproduced by construction.
    **They do not test shear deformation either**: the source's acceptance tolerance (3 % / 5 %) is wider than the
    whole effect of shear deformation in its tables (at most 1.3 % on a moment, 3.5 % on a deflection), so every
    with-shear target would also be met by an analysis without shear deformation (stated in each `limitsNote`).
  * A does not enter the results in this theory.
* Support idealisation of Case 1 (pin at the bottom, lateral roller at the loaded top) is read from the sketch.

The problems therefore verify **the second-order amplification** (what the benchmark is for), not section
properties.

## Problems

One problem per case, table, shear option and axial load level: 2 x 2 x 2 x 4 = 32 problems with two targets each
(64 tabulated numbers = every moment and deflection of both figures). US and SI tables are independently rounded
and are separate problems.

Result: **all 64 targets are reproduced** within the source's tolerance (6 of them by construction: the P = 0
deflections used for the calibration); the largest differences are 0.38 % on a
moment and 0.88 % on a deflection (both Case 1, SI table, P = 2000 kN, with shear deformation; the SI shear
calibration rests on one digit) - the size of the rounding of the printed values.

**Quarantined: none.** All 32 problems are active. The 16 models WITHOUT shear deformation also pass tier C
(OpenSeesPy and PyNite, within the error bound that tier C computes from three refinement levels of each foreign
solver; docs/CROSSCHECK.md 1.3). The models WITH shear deformation cannot be represented by the foreign solvers
in second order, so their model records are not cross-checked: their evidence is the targets (tier B), the
invariants and tier E only.

## How to read the index

`targets` / `missed` and the differences are those of `../check_local_v2b.py` (engine against the targets of
`problem.json`); "largest difference / tolerance" above 1 means a target is missed. `status` is the declared
status. "Verified" in this project means only: passed the tiers listed in `docs/VERIFICATION.md` at the
tolerances stated there. The other tiers (cross-check with foreign solvers, invariants, sensitivity) are run by
`benchfea verify verification/problems/aisc-commentary/<id>`.

Tolerance policy (docs/VERIFICATION.md, section 5), applied uniformly and checked by `../check_local_v2b.py`:
a published number carries half a unit in its last printed digit (floor 1e-6), or the acceptance tolerance the source itself states (here: 3 % / 5 %); a number computed here
carries 1e-6. A problem that misses, or on whose model another tier does not pass, is quarantined with a written
diagnosis (`quarantine.py`); nothing is loosened or dropped.

## Index

| problem | source location | what it verifies | basis | targets | missed | largest rel. difference | largest difference / tolerance | status |
|---|---|---|---|---|---|---|---|---|
| `case-1-si-with-shear-p0` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), SI table, column P = 0, unbracketed values | AISC benchmark Case 1 (P-delta only), SI table, with shear deformation, P = 0 | published | 2 | 0 | 1.6e-03 | 0.05 | active |
| `case-1-si-with-shear-p1330` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), SI table, column P = 1330, unbracketed values | AISC benchmark Case 1 (P-delta only), SI table, with shear deformation, P = 1330 | published | 2 | 0 | 4.4e-03 | 0.09 | active |
| `case-1-si-with-shear-p2000` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), SI table, column P = 2000, unbracketed values | AISC benchmark Case 1 (P-delta only), SI table, with shear deformation, P = 2000 | published | 2 | 0 | 8.8e-03 | 0.18 | active |
| `case-1-si-with-shear-p667` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), SI table, column P = 667, unbracketed values | AISC benchmark Case 1 (P-delta only), SI table, with shear deformation, P = 667 | published | 2 | 0 | 6.4e-03 | 0.13 | active |
| `case-1-si-without-shear-p0` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), SI table, column P = 0, values in brackets | AISC benchmark Case 1 (P-delta only), SI table, without shear deformation, P = 0 | published | 2 | 0 | 1.6e-03 | 0.05 | active |
| `case-1-si-without-shear-p1330` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), SI table, column P = 1330, values in brackets | AISC benchmark Case 1 (P-delta only), SI table, without shear deformation, P = 1330 | published | 2 | 0 | 7.8e-04 | 0.03 | active |
| `case-1-si-without-shear-p2000` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), SI table, column P = 2000, values in brackets | AISC benchmark Case 1 (P-delta only), SI table, without shear deformation, P = 2000 | published | 2 | 0 | 2.2e-03 | 0.04 | active |
| `case-1-si-without-shear-p667` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), SI table, column P = 667, values in brackets | AISC benchmark Case 1 (P-delta only), SI table, without shear deformation, P = 667 | published | 2 | 0 | 3.0e-03 | 0.06 | active |
| `case-1-us-with-shear-p0` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), US customary table, column P = 0, unbracketed values | AISC benchmark Case 1 (P-delta only), US table, with shear deformation, P = 0 | published | 2 | 0 | 8.5e-04 | 0.03 | active |
| `case-1-us-with-shear-p150` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), US customary table, column P = 150, unbracketed values | AISC benchmark Case 1 (P-delta only), US table, with shear deformation, P = 150 | published | 2 | 0 | 5.2e-03 | 0.10 | active |
| `case-1-us-with-shear-p300` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), US customary table, column P = 300, unbracketed values | AISC benchmark Case 1 (P-delta only), US table, with shear deformation, P = 300 | published | 2 | 0 | 3.9e-03 | 0.08 | active |
| `case-1-us-with-shear-p450` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), US customary table, column P = 450, unbracketed values | AISC benchmark Case 1 (P-delta only), US table, with shear deformation, P = 450 | published | 2 | 0 | 7.9e-03 | 0.16 | active |
| `case-1-us-without-shear-p0` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), US customary table, column P = 0, values in brackets | AISC benchmark Case 1 (P-delta only), US table, without shear deformation, P = 0 | published | 2 | 0 | 8.5e-04 | 0.03 | active |
| `case-1-us-without-shear-p150` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), US customary table, column P = 150, values in brackets | AISC benchmark Case 1 (P-delta only), US table, without shear deformation, P = 150 | published | 2 | 0 | 2.9e-03 | 0.06 | active |
| `case-1-us-without-shear-p300` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), US customary table, column P = 300, values in brackets | AISC benchmark Case 1 (P-delta only), US table, without shear deformation, P = 300 | published | 2 | 0 | 1.7e-03 | 0.06 | active |
| `case-1-us-without-shear-p450` | Commentary Section C2.1, Fig. C-C2.2 (benchmark problem Case 1), US customary table, column P = 450, values in brackets | AISC benchmark Case 1 (P-delta only), US table, without shear deformation, P = 450 | published | 2 | 0 | 2.2e-03 | 0.04 | active |
| `case-2-si-with-shear-p0` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), SI table, column P = 0, unbracketed values | AISC benchmark Case 2 (P-Delta and P-delta), SI table, with shear deformation, P = 0 | published | 2 | 0 | 1.1e-03 | 0.04 | active |
| `case-2-si-with-shear-p445` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), SI table, column P = 445, unbracketed values | AISC benchmark Case 2 (P-Delta and P-delta), SI table, with shear deformation, P = 445 | published | 2 | 0 | 1.4e-03 | 0.03 | active |
| `case-2-si-with-shear-p667` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), SI table, column P = 667, unbracketed values | AISC benchmark Case 2 (P-Delta and P-delta), SI table, with shear deformation, P = 667 | published | 2 | 0 | 2.9e-03 | 0.06 | active |
| `case-2-si-with-shear-p890` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), SI table, column P = 890, unbracketed values | AISC benchmark Case 2 (P-Delta and P-delta), SI table, with shear deformation, P = 890 | published | 2 | 0 | 9.0e-04 | 0.02 | active |
| `case-2-si-without-shear-p0` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), SI table, column P = 0, values in brackets | AISC benchmark Case 2 (P-Delta and P-delta), SI table, without shear deformation, P = 0 | published | 2 | 0 | 1.1e-03 | 0.04 | active |
| `case-2-si-without-shear-p445` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), SI table, column P = 445, values in brackets | AISC benchmark Case 2 (P-Delta and P-delta), SI table, without shear deformation, P = 445 | published | 2 | 0 | 1.4e-03 | 0.03 | active |
| `case-2-si-without-shear-p667` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), SI table, column P = 667, values in brackets | AISC benchmark Case 2 (P-Delta and P-delta), SI table, without shear deformation, P = 667 | published | 2 | 0 | 9.5e-04 | 0.03 | active |
| `case-2-si-without-shear-p890` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), SI table, column P = 890, values in brackets | AISC benchmark Case 2 (P-Delta and P-delta), SI table, without shear deformation, P = 890 | published | 2 | 0 | 5.6e-03 | 0.11 | active |
| `case-2-us-with-shear-p0` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), US customary table, column P = 0, unbracketed values | AISC benchmark Case 2 (P-Delta and P-delta), US table, with shear deformation, P = 0 | published | 2 | 0 | 3.7e-16 | 0.00 | active |
| `case-2-us-with-shear-p100` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), US customary table, column P = 100, unbracketed values | AISC benchmark Case 2 (P-Delta and P-delta), US table, with shear deformation, P = 100 | published | 2 | 0 | 2.2e-03 | 0.04 | active |
| `case-2-us-with-shear-p150` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), US customary table, column P = 150, unbracketed values | AISC benchmark Case 2 (P-Delta and P-delta), US table, with shear deformation, P = 150 | published | 2 | 0 | 1.2e-03 | 0.04 | active |
| `case-2-us-with-shear-p200` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), US customary table, column P = 200, unbracketed values | AISC benchmark Case 2 (P-Delta and P-delta), US table, with shear deformation, P = 200 | published | 2 | 0 | 3.1e-03 | 0.06 | active |
| `case-2-us-without-shear-p0` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), US customary table, column P = 0, values in brackets | AISC benchmark Case 2 (P-Delta and P-delta), US table, without shear deformation, P = 0 | published | 2 | 0 | 8.5e-16 | 0.00 | active |
| `case-2-us-without-shear-p100` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), US customary table, column P = 100, values in brackets | AISC benchmark Case 2 (P-Delta and P-delta), US table, without shear deformation, P = 100 | published | 2 | 0 | 7.5e-04 | 0.02 | active |
| `case-2-us-without-shear-p150` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), US customary table, column P = 150, values in brackets | AISC benchmark Case 2 (P-Delta and P-delta), US table, without shear deformation, P = 150 | published | 2 | 0 | 1.2e-03 | 0.04 | active |
| `case-2-us-without-shear-p200` | Commentary Section C2.1, Fig. C-C2.3 (benchmark problem Case 2), US customary table, column P = 200, values in brackets | AISC benchmark Case 2 (P-Delta and P-delta), US table, without shear deformation, P = 200 | published | 2 | 0 | 2.4e-03 | 0.05 | active |

## Open questions

1. Which I_x, A and shear area (and which definition of the shear area) did the Commentary use for the W14x48?
   Not stated; needs an opened section-property source (AISC Manual / Shapes Database were not available).
2. Is the Commentary's reference solution closed-form or computed, and with which shear formulation?
