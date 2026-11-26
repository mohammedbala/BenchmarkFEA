# References

This file records, for every reference this project names, whether it was **opened** (read from
the document itself during the project) or **not opened**, and what in the repository rests on it.

Rules that apply to all of them:

* Nothing in this repository is taken from a reference that was not opened, and no reference value
  is taken from memory.
* No text, table or figure of any reference is copied into the repository. A reference is cited
  (work, edition, section / equation / example number, page), and the value is then computed here
  or entered as a target number with its location.
* "Rests on" means: the code implements equations cited from that reference, or a verification
  problem compares the engine with numbers printed in it. It does not mean that the authors or
  publishers of a reference have anything to do with this project, or have reviewed it.
* "Verified" means only: passed the tiers listed in [VERIFICATION.md](VERIFICATION.md), at the
  tolerances stated there.

Statements followed by the HTML comment `<!-- RECHECK -->` (visible in the source of this file)
give counts that were changing while the file was written.

Problem counts below are directories under `verification/problems/`; the split into active and
quarantined is that of the `status` field of each `problem.json` at the time of writing.
<!-- RECHECK --> The current state is always in
[VERIFICATION_REPORT.md](VERIFICATION_REPORT.md).

---

## 1. Opened

### 1.1 Specifications, books and data

| # | Reference | How obtained | What was read | What rests on it | Licence / redistribution |
|---|---|---|---|---|---|
| R1 | American Institute of Steel Construction, *Specification for Structural Steel Buildings*, ANSI/AISC 360-22, dated August 1, 2022 (with Commentary; 780-page PDF; printed pages 16.1-N) | PDF supplied by the project owner. The AISC web site could not be opened by this project (it answered automated requests with an access challenge). | Read from the rendered pages, not from a text extraction. **Specification**: Chapter B, 16.1-15 to 16.1-25 (B3, B4.1 with Tables B4.1a/b, B4.2, B4.3); Chapter C, 16.1-26 to 16.1-31; Chapter D, 16.1-32 to 16.1-37; Chapter E, 16.1-38 to 16.1-49 (complete; E5, E6 noted at summary level); Chapter F, 16.1-50 to 16.1-74 (complete; F4, F5, F9–F13 noted at summary level); Chapter G, 16.1-75 to 16.1-81 (complete); Chapter H, 16.1-82 to 16.1-88; Chapter L, 16.1-174 to 16.1-175; Appendix 7, 16.1-279 to 16.1-282; Appendix 8, 16.1-283 to 16.1-286. **Commentary**: to B3 and B4 (16.1-313 to 16.1-330), C (16.1-333 to 16.1-344, including the benchmark problems of C2.1 on 16.1-337), D (16.1-345 to 16.1-352), E (16.1-353 to 16.1-369), F (16.1-370 to 16.1-395), G (16.1-396 to 16.1-403), H (16.1-404 to 16.1-417), L (16.1-542 to 16.1-547), Appendix 7 (16.1-648 to 16.1-658), Appendix 8 (16.1-659 to 16.1-666). Not read: Chapters A, I, J, K, M, N, Appendices 1–6 and their Commentary; the glossary was seen in a text extraction only. | **Code**: `src/benchfea/design/aisc360_22/` (all modules: classification, tension, compression, flexure, shear, combined forces, serviceability); `src/benchfea/solve/nonlinear.py` (C2.1(d) load multiplier, C2.3 stiffness reduction and Eq. C2-2a/b, C2.2b Eq. C2-1 notional loads); slenderness defaults in `src/benchfea/defaults.py`. **Tests**: the hand calculations `test_hand_*` in `tests/test_design_*.py` (tier F, side a); `tests/test_second_order_aisc_benchmarks.py`. **Problems**: `verification/problems/aisc-commentary/` — 32 problems, 64 targets (Commentary C2.1 Cases 1 and 2), all active. <!-- RECHECK --> | Copyright AISC; not redistributable. No AISC text, table or figure is in the repository; the user supplies their own copy to check a citation. |
| R2 | American Institute of Steel Construction, *Companion to the AISC Steel Construction Manual, Volume 1: Design Examples*, Version 16.0 (© AISC 2023; 1054-page PDF). Its preface names ANSI/AISC 360-22 and the 16th Edition Manual as its basis. No revision or errata date is printed on the front pages read. | PDF supplied by the project owner (same access problem as R1). | Read from the rendered pages: front matter (PDF pp. 1–6); Chapter C, C-1 to C-17 (Examples C.1A, C.1B, C.1C); Chapter D, D-1 to D-19 (D.1–D.6); Chapter E (E.1–E.4, E.9–E.12 read; E.5–E.8, E.13, E.14 from the index only); Chapter F, F-1 to F-44 (F.1–F.9) and F-65 to F-73 (F.14); Chapter G (G.1–G.7; G.8 on G-17, G-18); Chapter H, H-1 to H-33; Example A-6.5 (A6-19, A6-21); Part III, pages III-84 to III-86, III-90, III-91, III-93, III-94. | **Tests** (tier F, side b): `tests/design_examples/aisc_design_examples_v16*.json` — 43 case entries, 40 distinct examples or cases, numbers only with page references — used by `tests/test_design_*.py` and `tests/test_design_review.py` with one tolerance (1 % on a printed value). **Coverage registry**: `src/benchfea/verify/design_coverage.json`. **Problems**: `verification/problems/aisc-de/` — 10 problems, 70 targets (frame of Examples C.1A–C), **all quarantined**: the examples do not state the shear areas or whether their software includes shear deformation. | Copyright AISC; not redistributable. Only example identifiers, page numbers, input numbers and printed result numbers are in the repository. |
| R3 | W. McGuire, R. H. Gallagher, R. D. Ziemian, *Matrix Structural Analysis*, 2nd edition. The rights page of the electronic edition read says "Copyright © 2014 by Ronald D. Ziemian"; no publisher or ISBN is printed on the pages read. | PDF supplied by the project owner (file "Part 1"). The authors offer the electronic edition free of charge through the MASTAN2 site and a university repository; the automated download was blocked, so the owner supplied it. | Read as page images (PDF page = printed page + 21): front matter; Chapters 2–5; Chapter 7 (Sections 7.3–7.6 in part); Chapter 8; Chapter 9 including its problems; Chapter 11 Section 11.6; Chapter 12 Sections 12.1–12.5, 12.7, 12.8; Chapter 13; Appendix A. Not read: Chapter 6, Chapter 10, Section 12.6, the end-of-chapter problems of Chapters 2, 3 and 12, the contents and symbol lists. | **Problems**: `verification/problems/mgz/` — 46 problems, 277 targets; 41 with basis `published`, 5 `recomputed-from-published-statement` (a closed-form expression of the book evaluated here); 22 active, 24 quarantined (18 of them for rounding in the book's hand calculations). <!-- RECHECK --> **Code**: none. The engine's formulation is derived in [THEORY.md](THEORY.md) and does not cite the book. | All rights reserved by the rights holder; the rights page forbids reproduction. Only example numbers, page numbers, problem data restated in this project's words and printed result numbers are in the repository. |
| R4 | C. W. Ziemian, R. D. Ziemian, "Steel benchmark frames for structural analysis and validation studies: Finite element models and numerical simulation data", *Data in Brief* 39 (2021) 107564, doi:10.1016/j.dib.2021.107564. Data set: Mendeley Data, doi:10.17632/39sjhchwtx.1. | Article read on the open-access page of the publisher's archive (PMC8637320); data set files downloaded from data.mendeley.com (public). | Article: text and Tables 2 and 3 (analysis conditions). Data set: "Description Of Frames.pdf" (24 pages), "FEA Model Details_All frames.xlsx" and "Simulation Results_All frames.xlsx" — the tabs of frames 1, 2, 3, 5, 7, 9 and the summary tab. | **Problems**: `verification/problems/ziemian/` — 36 problems, 164 targets (first-order displacements and moments, second-order results, critical load ratios of six frames); 14 active, 22 quarantined (all 12 second-order problems: the data set's analysis is a large-displacement analysis). **Diagnosis script**: `verification/derive/ziemian_diagnosis.py`. | **CC BY 4.0.** The project owner approved embedding a few frames with attribution. The model data of six frames and the result numbers used as targets are reproduced, with changes, under that licence; the attribution, the licence link and the list of changes are in `verification/problems/ziemian/README.md` and in every `problem.json` of that directory (`source.notes`). Anyone redistributing those files must keep that attribution. The two spreadsheets themselves are not in the repository. |
| R5 | National Institute of Standards and Technology, Special Publication 811 (guide for the use of the SI) — Appendix B.8 and its footnotes. The edition and date of the page read were not recorded in the project notes. | Official public page, nist.gov | Appendix B.8 (conversion factors) and footnotes 22 and 23: inch and foot (exact), pound-force (exact value in newtons) | `src/benchfea/units.py` (exact conversion factors as fractions) | Work of the United States government; public |
| R6 | JSON Schema, draft 2020-12: Core (draft-bhutton-json-schema-01, 16 June 2022), Validation, and the meta-schema | Official public pages, json-schema.org/draft/2020-12/ | Core sections 4.3.1, 6.5, 8, 10, 11; Validation sections 5, 6, 7, 9; the meta-schema document | `src/benchfea/schema/` (the input schema and its generator), `src/benchfea/validate.py`. Checked with the installed `jsonschema` 4.26.0 | Public specification |
| R7 | Gymnasium documentation (Farama Foundation): API pages `api/env`, `api/registry`, `api/spaces`, `api/utils` | Official public pages, gymnasium.farama.org; signatures compared with the installed package `gymnasium` 1.4.0 | The four pages named | `src/benchfea/bench/env.py` (`BenchFEAEnv`), `tests/test_env.py` | Open-source package (optional dependency); documentation public |

### 1.2 Vendor verification manuals (problem statements only)

| # | Reference | How obtained | What was read | What rests on it | Licence / redistribution |
|---|---|---|---|---|---|
| R8 | Bentley Systems, *STAAD.Pro Help*, "Verification Examples" (the help build that refers to the STAAD.Pro 2023 sample files), read 2026-11-22 | Official public documentation site, docs.bentley.com (`LiveContent/web/STAAD.Pro Help-v21/en/…`). The newer help build could not be read. | Sections V.01 Beams, V.02 Trusses, V.03 Frames and V.07 Nonlinear Analysis: the pages listed in `verification/problems/staad/README.md` (19 statements used; the pages not used and the reason for each are listed there). V.04–V.06 and V.08–V.11 not read. Sketches were not readable; a page is used only where its text or input listing fixes the problem. | **Problems**: `verification/problems/staad/` — 24 problems, 138 targets, from 19 statements; 19 with basis `recomputed-from-published-statement`, 5 separate problems holding vendor-printed values that do not reproduce; 18 active, 6 quarantined. <!-- RECHECK --> The reference values of the recomputed targets are computed by `verification/derive/vendor_reference.py`. The program's own results are never targets. | © Bentley Systems; no reuse licence was seen. Cited, not reproduced: page title and address are in each `problem.json`. |
| R9 | ANSYS, Inc., *Mechanical APDL Verification Manual*, Release 2025 R2, read 2026-11-22 | Official public help site, ansyshelp.ansys.com (pages `Hlp_V_VM<n>.html`) | Test cases VM1, VM2, VM4, VM10, VM12, VM14, VM18, VM21, VM41, VM50, VM59, VM127, VM136, VM216, VM217, VM222 (text only; sketches not readable). Used: VM2, VM4, VM21, VM41, VM127, VM217. | **Problems**: `verification/problems/ansys/` — 8 problems, 37 targets, from 6 statements; 6 active, 2 quarantined (vendor-printed values that do not reproduce). Reference values by `verification/derive/vendor_reference.py`. | © ANSYS, Inc., all rights reserved. Cited, not reproduced. |

The textbooks that the vendor pages cite as the origin of their problems (named on each page and
in the corresponding `problem.json`) were **not opened**. This is why those problems have the basis
`recomputed-from-published-statement`: the statement comes from the vendor page, the target values
are computed in this repository, and the result is evidence of the same kind as tier A — two
computations by one project — not a reproduction of a textbook value.

### 1.3 Documentation of the independent solvers

| # | Reference | How obtained | What was read | What rests on it |
|---|---|---|---|---|
| R10 | OpenSeesPy documentation and OpenSees documentation | Official public pages: openseespydoc.readthedocs.io (`elasticBeamColumn`, `ElasticTimoshenkoBeam`, `LinearTransf`, `pdelta`, `corotational`, `geomTransf`, `eleload`, `rigidLink`, `rigidDiaphragm`, `equalDOF`, the constraint handler pages, `trussEle`, `ZeroLength`, `eigen`, `reactions`, `nodeReaction`, `eleForce`, `eleResponse`, `basicForce`); opensees.berkeley.edu wiki (`Elastic_Beam_Column_Element`, `Linear_Transformation`, `PDelta_Transformation`, `RigidDiaphragm_command`); opensees.github.io user manual (`elasticBeamColumn`, `PDelta`). Installed package `openseespy` 3.8.0.0 and its licence file. | The pages named. Where the documentation ("latest") and the installed version differ, the measured behaviour of version 3.8.0 is used and is pinned by `tests/test_crosscheck_measured.py`. | `src/benchfea/crosscheck/opensees.py`, `capabilities.py`; tier C |
| R11 | PyNite (PyNiteFEA) documentation and source | Documentation at pynite.readthedocs.io (`index`, `PDelta`, `member`, `stability`, `analysis`, `node`); source of the installed package `PyNiteFEA` 3.2.1 (`Member3D.py`, `FEModel3D.py`, `Analysis.py`, `PhysMember.py`, `Spring3D.py`: the parts on member stiffness, geometric stiffness, releases, loads, P-Delta procedure, reactions) | As named | `src/benchfea/crosscheck/pynite.py`, `capabilities.py`; tier C. No line of PyNite's source is in the repository. |
| R12 | Frame3DD: project page and manual (frame3dd.sourceforge.net; manual of the svn trunk, sections 7.2–7.13), parts of its source read as text | Official project site | Read for a formulation summary only. **The program was not run** (see 2.2). | One entry `not-run` in `src/benchfea/crosscheck/capabilities.py`. No result, formula or code of Frame3DD is in the repository. |

Full addresses of the solver documentation pages: `benchfea.crosscheck.capabilities.SOURCES`.

### 1.4 What rests on no reference

* The engine (`src/benchfea/solve/`): derived from the mechanics in [THEORY.md](THEORY.md); no
  textbook is cited and no code of another program was ported.
* The member-equation oracle (`src/benchfea/oracle/`), the fixed-end-force routine
  (`src/benchfea/crosscheck/fef.py`) and the "tableau" solver (`tests/independent/`): derived in
  their module docstrings.
* `verification/problems/derived/` — 17 tier-A problems, 101 targets, basis `derived-here`
  (computed by the oracle or by a numerically solved characteristic equation, script in
  `verification/derive/derive_problems.py`). They are never attributed to a book.
* Section property formulas in `src/benchfea/sections.py`: derived there; the shear areas are
  simple definitions, not the coefficients of any publication.
* The names "Euler–Bernoulli", "Timoshenko", "Engesser", "Haringx", "St. Venant", "Maxwell–Betti",
  "Richardson" are used as the common names of formulations and theorems. No publication under
  those names was opened or is cited.

---

## 2. Not opened

For each work below: **not opened; nothing in this repository is taken from it.** It must not be
cited as the source of any number, formula or statement here. The last column says what would
become possible if a user supplied it.

### 2.1 Books, papers, standards and data

| Work | Why it was not opened | What would become possible with it |
|---|---|---|
| A. Kassimali, *Matrix Analysis of Structures* | commercial textbook; no official open copy found | further published worked examples for tier B (plane and space frames, trusses, grids) |
| J. S. Przemieniecki, *Theory of Matrix Structural Analysis* | commercial (original and reprint); no official open copy found | a published statement of element stiffness and geometric stiffness matrices to compare the derivations of THEORY.md with |
| R. D. Cook, D. S. Malkus, M. E. Plesha, R. J. Witt, *Concepts and Applications of Finite Element Analysis* | commercial textbook | as above; published treatment of constraints and conditioning |
| C. A. Felippa, *Introduction to Finite Element Methods* (course notes) | the author's university pages no longer exist; the only copy found is on a third-party site, which is not an official location | as above |
| *Roark's Formulas for Stress and Strain* | commercial; electronic edition by site licence only | published closed-form beam cases as tier-B targets (today the closed forms are tier A, derived here); two vendor pages cite it |
| S. P. Timoshenko, J. M. Gere, *Theory of Elastic Stability* | commercial (original and reprint) | published buckling loads and beam-column solutions as tier-B targets for buckling, which today has only tier-A evidence; a published basis for the choice between the two shear formulations in second order |
| W. F. Chen, E. M. Lui, *Structural Stability* | commercial | published stability functions and frame examples as tier-B targets |
| G. R. Cowper, "The shear coefficient in Timoshenko's beam theory", *Journal of Applied Mechanics* 33 (1966) 335–340 | the publisher's page refused access; no open copy exists (bibliographic data confirmed from a citation index only) | published shear coefficients; the shear areas in `sections.py` could then be offered as that definition instead of the simple one |
| AISC *Steel Construction Manual*, 16th edition | commercial; AISC site not reachable | the tabulated section properties and strength tables that the Design Examples read; many Design Examples cases that are now partial reproductions (they take a tabulated value as input) could become full ones |
| AISC Shapes Database v16.0 | AISC site not reachable; no terms of use seen | (a) confirm the column names and units assumed by `src/benchfea/shapes.py`; (b) section properties of W14x48, W12x65, W18x40 and others, which would turn the calibrated AISC Commentary problems into real ones and allow the ten quarantined Design Examples frame problems to be rebuilt; (c) tier-F coverage of more limit states with real shapes |
| AISC Design Guide 28, *Stability Design of Steel Buildings* | AISC site not reachable | published worked frames for the direct analysis method (notional loads, stiffness reduction) as tier-B targets; guidance on questions recorded in OPEN_QUESTIONS.md about notional loads and τ_b |
| R. D. Ziemian (ed.), *Guide to Stability Design Criteria for Metal Structures* | commercial | background to the Specification's stability provisions |
| NAFEMS, *The Standard NAFEMS Benchmarks* (publication P18) | paid download; only the public catalogue page (titles, no data) was readable | independent published benchmarks |
| CSI verification manuals (SAP2000 / ETABS) | the official wiki returned no content to this project; the manuals ship with the software; only third-party copies were found, which were not used | a third family of vendor problem statements |
| ASCE/SEI 7 (referred to by the Commentary to Chapter L) | not supplied | load combinations and serviceability guidance; today limits and combinations are user input |
| ANSI/AISC 341 (seismic provisions) and the AISC Specification chapters not listed under R1 | not supplied / not read | nothing is implemented from them |
| The textbooks cited by the STAAD and ANSYS pages (strength of materials, structural analysis and statics textbooks named on each page) | not supplied | the vendor problems could change basis from `recomputed-from-published-statement` to `published`, and misprints or hand-method approximations in the vendor values could be traced |
| The original literature of each Ziemian frame (named in the data set's description file) | not looked for | the unrounded loads and the member numbering asked about in OPEN_QUESTIONS.md |

### 2.2 Programs that could not be run

| Program | Why | Consequence |
|---|---|---|
| Frame3DD | no binary and no installable package for the Python version in use; compiling and executing downloaded source was not permitted in the project environment. Its manual and parts of its source were read (R12). | not used as a solver; it would be a third independent solver with shear deformation and its own geometric stiffness, which would make shear-deformable models cross-checkable |
| MASTAN2 | needs MATLAB (or its runtime); none available. Not obtained. | not used; it is the program behind the Ziemian data set, so running it would resolve the remaining differences on those frames, and it would provide buckling factors from an outside program |

---

## 3. The independent solvers actually used

The project specification requires at least two independent solvers. **This requirement is met by
OpenSeesPy and PyNite.** The in-repository oracle, the tableau solver and the separate reference
programs of `verification/derive/` are by the same author as the engine; they are useful evidence
(they share no solver code with the engine) but they are **not counted** towards the two.

Both packages are optional (`pip install benchfea[crosscheck]`), are run in a subprocess, and are
never imported by the engine. No line of their source is in the repository.

| | OpenSeesPy | PyNite |
|---|---|---|
| Version used | `openseespy` 3.8.0.0 | `PyNiteFEA` 3.2.1 |
| Licence (as read in the installed package) | free for research, education and internal use; commercial redistribution needs a licence — optional test dependency only, never bundled | MIT |
| Element | `elasticBeamColumn`: prismatic three-dimensional Euler–Bernoulli element (no shear deformation); `ElasticTimoshenkoBeam` where shear deformation is asked for | 12 × 12 prismatic Euler–Bernoulli element |
| Shear deformation | first order only, through `ElasticTimoshenkoBeam`. Measured in 3.8.0: that element ignores the P-Delta transformation, drops the rigid arms of joint offsets, takes no release flags and no point loads | none (stated in its documentation) |
| Geometric stiffness (second order) | `PDelta` transformation: the chord-rotation term only (axial force / length on the end translations). The curvature of the member between its ends (P-δ) is not in one element; it is obtained by cutting each member into 32 elements (8, 16, 32 for the error bound). Axial forces iterated by Newton's method. The axial force of an element is that of its end j without the axial part of an element load, so the translator applies axial member loads at the element nodes | consistent cubic geometric stiffness in both bending planes, which captures P-δ approximately within one element; members are cut into 8 sub-members (2, 4, 8 for the error bound). PyNite's own matrix has two further terms (on the axial and on the torsional stiffness) and its own procedure is a single solve with first-order axial forces; the translator removes the two terms by adjusting section properties and iterates the axial forces to convergence, using PyNite's own member matrices, partitioning and solver |
| Constraints | `rigidLink('beam')`, `rigidDiaphragm`, `equalDOF`, enforced exactly by the `Transformation` handler (elimination). The same linear relation as the engine's, also in second order | none. Rigid links, diaphragms and offsets are emulated by very stiff members (error about 1/multiplier, with round-off growing with the multiplier): `approximate` in first order; `unsupported` in second order, because stiff members solve the rigid-body problem and not the linear constraint |
| Rigid offsets | `-jntOffset` on `elasticBeamColumn`; an added node with a rigid link otherwise. In second order: no error bound, not eligible | stiff stub member, `approximate` |
| Releases | element flags for full bending releases only. Rotational springs: an exact decomposition into parallel elements, computed by the translator. Axial, shear and torsion releases: a duplicate node with `equalDOF` (exact when the member axes are global, there is no offset and the joint is not a slave), otherwise penalty springs (`approximate`) | boolean release of any of the twelve end components (static condensation). Rotational springs through the same parallel-member decomposition. Partial axial, shear or torsion connections: `unsupported` |
| Member loads | native element loads only where measured to be exact (uniform, partial uniform and point forces on unreleased Euler–Bernoulli elements). Everything else as equivalent nodal loads computed by `benchfea.crosscheck.fef`, with the fixed-end forces added back | PyNite's own point and linearly varying member loads, on sub-members cut by the translator |
| Tension-/compression-only members | not reproduced. Tier C gives the solver the linear model of the engine's final set of active members and checks the sign conditions on the solver's results | same |
| Linear buckling | not implemented in the translator (feasible in principle from two tangent matrices) | no buckling analysis in the package |
| Linear equation solver | `UmfPack` (sparse). Measured: it does not report a singular matrix | sparse |
| Self-checks added by the translator | finite displacements; every constraint satisfied to 10⁻⁸; joint equilibrium to 10⁻⁶; reactions computed by the translator from joint equilibrium | member axes reproduced to 10⁻⁹; finite displacements |

Known differences from the engine, and how each is handled:

* **The engine's member is exact for the member differential equation** (first order, and second
  order through the transfer relation of the beam-column equation): one element per member, shear
  deformation optional, the P-δ effect included without subdivision. Both outside solvers
  discretise the second-order problem; their results carry an error bound computed from three
  refinement levels of that solver alone, and the comparison allows 10⁻⁶ plus the bounds of the
  two results compared. A solver whose bound exceeds 5·10⁻³ is not eligible. <!-- RECHECK -->
* **Shear deformation**: one eligible solver in first order, none in second order → not
  cross-checked.
* **Rigid links and diaphragms in second order**: OpenSees and the engine implement the same
  linear relation and agree; PyNite's emulation is the rigid body, which differs (7.6 % on one
  test model) → one eligible solver, not cross-checked. Which statement is the right one is an
  open question.
* **End shears in second order**: OpenSees reports shears normal to the displaced chord; the
  translator converts to components along the undeformed axes, which is what the engine reports.
* **Truss-type members in second order**: the engine treats them as pin-ended beam-columns; the
  translators subdivide such a member when it carries a load along its axis.
* **Shared code between the two routes**: the fixed-end-force and parallel-member routines
  (`benchfea.crosscheck.fef`) and the preparation of the model (`benchfea.crosscheck.prepare`,
  which takes loads and axes as resolved by `benchfea.model`). An error there would be common to
  both routes. See [../LIMITS.md](../LIMITS.md) section 4.5.

Details, measured behaviours and the capability table: [CROSSCHECK.md](CROSSCHECK.md) sections 3
to 5; machine-readable in `benchfea.crosscheck.CAPABILITIES`.

---

## 4. Other software the results depend on

Python 3.13, `numpy` 2.5.3, `scipy` 1.18.1 (dense linear algebra, `solve_bvp` and `solve_ivp` for
the oracle and the station checks), `jsonschema` 4.26.0, `pyyaml` 6.0.3; optional `gymnasium`
1.4.0 and `openpyxl`. These are the versions installed where the test suite was run; every
verification record stores the versions it was made with (`versions`). Their documentation is not
a reference for any engineering number.
