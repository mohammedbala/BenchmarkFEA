# Data card — BenchFEA tasks

Every number on this card stands between `generated` markers and is written by
`python3 tools/build_dataset_docs.py` from `data/tasks/manifest.json`, the task files and
`data/CHECKS.json`; `--check` (run by the test suite) fails when the card and the data disagree.
The text outside the markers is prose and holds no count.

## What it is

<!-- generated:headline -->
**748 verified tasks** from 263 verified cells of 528 (every family × every curve ball × seeds 0 and 4): 236 predict, 263 diagnose, 249 repair, 0 design; 348 train, 400 held-out.
<!-- /generated:headline -->

Benchmark tasks for agents that work on structural frame models, each with a hidden reference. A
task hands the agent a model of a steel frame (JSON input file of the `benchfea` engine, format
1.0.0) in which one thing was changed in a way that matters (a *curve ball*), and asks it to
**predict** named results, to **diagnose** what was changed and where (or, for an unstable model,
the mechanism), to **repair** the model, or to **design** its member sizes. Task types, answer
formats and scoring: [docs/BENCHMARK.md](../docs/BENCHMARK.md).

```
data/tasks/manifest.json                       counts, every task, every cell, the quarantine log, versions
data/tasks/{train,heldout}/<id>.task.json      public: statement, answer format, quantities / criteria
data/tasks/{train,heldout}/<id>.model.json     public: the input file of the task (one line of JSON)
data/answers/{train,heldout}/<id>.answers.json hidden: change record, reference solution, verification record
data/CHECKS.json                               summaries of the two checks run after the build
data/QUARANTINE_NOTES.md                       written diagnoses of the attempts logged as disagreement
data/user/                                     not part of the dataset (user-supplied files, never committed)
```

Hand the task and model files to an agent; keep `data/answers/` away from it. Task ids are opaque
(they do not name the curve ball). The manifest lists the curve ball of every cell next to its task
ids, so it belongs on the answers' side of that line too; the environment (`BenchFEA-v0`) does not
show it.

## How it was generated

<!-- generated:build -->
```
benchfea bench make-dataset --out data --seeds 0,4 --workers 4 --cache <dir> --clean
python3 tools/check_dataset_scores.py data --workers 4 --cache <dir>
python3 tools/audit_dataset.py data
python3 tools/build_dataset_docs.py
```

Versions of this build: benchfea 0.1.0, python 3.13.18, numpy 2.5.3, scipy 1.18.1, foreign:opensees 3.8.0.0, foreign:pynite 3.2.1, task 2.

Cells planned: 8 families × 33 curve balls × 2 seeds (0 and 4) = 528; with a verified instance: 263. Entries in the quarantine log of the manifest: 1138.
<!-- /generated:build -->

With the code of this repository; no file under `data/tasks` or `data/answers` is edited by hand.

1. **Model**: `benchfea.bench.generate(family, seed, difficulty)` — a seeded parametric generator
   of eight families with member sizes found by the engine's own sizing loop
   ([docs/DIFFICULTY.md](../docs/DIFFICULTY.md)). The models are synthetic: generic section
   catalogue, nominal steel properties, rounded loads. They are not real structures. The generator
   seed carries the split in its lowest bit, so the two splits never generate from the same seed.
2. **Curve ball**: `benchfea.bench.curveballs.apply(model, name, seed)` — one named change
   (supports, releases, offsets, constraints, sections, loads, units, second-order setting, design
   data), accepted only if its *measured* effect on the analysis or the design is material
   ([docs/CURVEBALLS.md](../docs/CURVEBALLS.md)). The exact change is recorded (JSON-pointer edits
   with before / after values).
3. **Verification**: the curve-balled model, and for repair tasks the model before the change, go
   through `benchfea.verify.verify_model` ([docs/VERIFICATION.md](../docs/VERIFICATION.md)).
   Only a reference that comes out `verified` is used. An instance that does not is skipped, the
   skip is logged, and the next seed of a documented chain is tried (at most three attempts); if
   none verifies, nothing is written for that cell. No tolerance was adjusted and no case was
   special-cased to make it pass.
4. **Tasks**: every task type that can be served on a verified instance becomes one task; the mode
   (`tool-allowed` / `no-solver`) is a fixed pseudo-random choice.

## Composition

<!-- generated:composition -->
| split | tasks | predict | diagnose | repair | design |
|---|---|---|---|---|---|
| train | 348 | 111 | 123 | 114 | 0 |
| held-out | 400 | 125 | 140 | 135 | 0 |
| all | 748 | 236 | 263 | 249 | 0 |

| mode | tasks | predict | diagnose | repair | design |
|---|---|---|---|---|---|
| `no-solver` | 364 | 106 | 128 | 130 | 0 |
| `tool-allowed` | 384 | 130 | 135 | 119 | 0 |

| family | tasks | train | held-out | cells | cells verified |
|---|---|---|---|---|---|
| `braced_frame` | 68 | 36 | 32 | 66 | 25 |
| `continuous_beam` | 81 | 36 | 45 | 66 | 28 |
| `diaphragm_frame` | 116 | 56 | 60 | 66 | 41 |
| `multistorey_2d` | 99 | 50 | 49 | 66 | 34 |
| `multistorey_3d` | 102 | 48 | 54 | 66 | 36 |
| `pipe_rack` | 95 | 43 | 52 | 66 | 34 |
| `portal_frame` | 94 | 39 | 55 | 66 | 33 |
| `stair_frame` | 93 | 40 | 53 | 66 | 32 |

| difficulty level | tasks | train | held-out |
|---|---|---|---|
| 1 | 63 | 27 | 36 |
| 2 | 134 | 69 | 65 |
| 3 | 240 | 101 | 139 |
| 4 | 196 | 104 | 92 |
| 5 | 115 | 47 | 68 |

| effect class | tasks | predict | diagnose | repair | design |
|---|---|---|---|---|---|
| `changes-analysis` | 695 | 236 | 236 | 223 | 0 |
| `creates-mechanism` | 53 | 0 | 27 | 26 | 0 |

| analysis order of the task's model | tasks | train | held-out | by family |
|---|---|---|---|---|
| first-order | 499 | 226 | 273 | `braced_frame` 27, `continuous_beam` 81, `diaphragm_frame` 116, `multistorey_2d` 35, `multistorey_3d` 49, `pipe_rack` 35, `portal_frame` 63, `stair_frame` 93 |
| second-order | 249 | 122 | 127 | `braced_frame` 41, `multistorey_2d` 64, `multistorey_3d` 53, `pipe_rack` 60, `portal_frame` 31 |

A model counts as second-order when at least one of its combinations is analysed second order.
<!-- /generated:composition -->

By curve ball (a cell is one family × curve ball × seed; it is verified when one attempt of its
chain gave a verified reference):

<!-- generated:curve-balls -->
| curve ball | cells | verified | tasks | train | held-out | cells not verified: stage:category of the last attempt |
|---|---|---|---|---|---|---|
| `constraint_added` | 16 | 2 | 6 | 6 | 0 | curve-ball:not-applicable 12, verify:not-representable 2 |
| `constraint_made_flexible` | 16 | 2 | 6 | 3 | 3 | curve-ball:not-applicable 14 |
| `constraint_removed` | 16 | 2 | 4 | 2 | 2 | curve-ball:not-applicable 14 |
| `constraint_slaves_dropped` | 16 | 2 | 6 | 3 | 3 | curve-ball:not-applicable 14 |
| `constraint_wrong_master` | 16 | 2 | 4 | 2 | 2 | curve-ball:not-applicable 14 |
| `constraint_wrong_slave` | 16 | 0 | 0 | 0 | 0 | curve-ball:not-applicable 14, verify:refused 2 |
| `design_bracing_points_changed` | 16 | 0 | 0 | 0 | 0 | verify:not-covered 8, curve-ball:not-applicable 7, verify:not-representable 1 |
| `design_k_changed` | 16 | 0 | 0 | 0 | 0 | verify:not-covered 13, curve-ball:not-applicable 3 |
| `design_lb_full_where_braced` | 16 | 0 | 0 | 0 | 0 | verify:not-covered 11, curve-ball:not-applicable 5 |
| `design_lb_short_where_unbraced` | 16 | 0 | 0 | 0 | 0 | curve-ball:not-applicable 14, verify:not-covered 2 |
| `design_unbraced_length_changed` | 16 | 0 | 0 | 0 | 0 | verify:not-covered 15, curve-ball:not-applicable 1 |
| `load_direction_changed` | 16 | 16 | 48 | 21 | 27 | - |
| `load_moved` | 16 | 16 | 48 | 15 | 33 | - |
| `load_sign_flipped` | 16 | 16 | 48 | 24 | 24 | - |
| `member_rolled_90` | 16 | 16 | 48 | 21 | 27 | - |
| `offset_added` | 16 | 12 | 31 | 17 | 14 | curve-ball:not-applicable 2, verify:not-representable 2 |
| `offset_removed` | 16 | 1 | 3 | 0 | 3 | curve-ball:not-applicable 15 |
| `offset_sign_reversed` | 16 | 0 | 0 | 0 | 0 | curve-ball:not-applicable 16 |
| `offset_wrong_axis` | 16 | 1 | 3 | 0 | 3 | curve-ball:not-applicable 15 |
| `release_added` | 16 | 15 | 45 | 21 | 24 | curve-ball:not-applicable 1 |
| `release_mechanism` | 16 | 16 | 31 | 14 | 17 | - |
| `release_pinned_frame` | 16 | 12 | 33 | 14 | 19 | curve-ball:not-applicable 4 |
| `release_removed` | 16 | 12 | 36 | 21 | 15 | curve-ball:not-applicable 4 |
| `second_order_off` | 16 | 4 | 9 | 9 | 0 | curve-ball:not-applicable 12 |
| `section_heavier` | 16 | 16 | 46 | 20 | 26 | - |
| `section_lighter` | 16 | 16 | 47 | 23 | 24 | - |
| `sections_swapped` | 16 | 14 | 42 | 18 | 24 | curve-ball:not-applicable 2 |
| `support_fixed_to_pinned` | 16 | 13 | 39 | 15 | 24 | curve-ball:not-applicable 3 |
| `support_pinned_to_fixed` | 16 | 12 | 36 | 18 | 18 | curve-ball:not-applicable 4 |
| `support_removed` | 16 | 15 | 40 | 19 | 21 | verify:not-representable 1 |
| `support_roller_direction` | 16 | 1 | 3 | 0 | 3 | curve-ball:not-applicable 15 |
| `support_spring_stiffness` | 16 | 13 | 39 | 18 | 21 | curve-ball:not-applicable 3 |
| `units_relabelled` | 16 | 16 | 47 | 24 | 23 | - |
<!-- /generated:curve-balls -->

## What "verified" means

**Only this: the reference passed every verification tier that applies to it, at the tolerances
stated in its record.** The record of each task (in its answers file) lists the tiers passed, the
tolerances and the largest differences found.

<!-- generated:verified -->
| tiers passed by the reference | tasks |
|---|---|
| `ACDE` | 66 |
| `CDE` | 629 |
| `NKCE` | 53 |

Tolerances (`benchfea.verify.POLICY`; relative to the scale stated in docs/VERIFICATION.md):

* tier A (ODE oracle, straight member lines): 1e-06;
* tier C (engine against OpenSeesPy and PyNite): 1e-06 plus the error bounds of the two results compared; a foreign result with a bound above 0.005 is not eligible; 2 eligible foreign solvers are required;
* tier D (invariants): superposition and reciprocity 1e-09, station values from statics 1e-09 (1e-06 in an iterated analysis), station values from an integration of the member equations 1e-06; equilibrium and energy at the tolerances the solver stores with each invariant;
* tier E (equivalent models): member split 1e-08 (1e-06 in an iterated analysis), limit of a stiffness series 1e-05, design variables bitwise;
* tier F (design results): every clause item of every design check used has a hand-calculation test and a reproduced strength of the AISC Design Examples v16.0 in the test suite (a coverage rule, no tolerance of its own);
* mechanism diagnoses (tiers N, K, C, E of the diagnosis route): null-space residual 2e-06, singular threshold 1e-09.
<!-- /generated:verified -->

What the tiers are: tier C — the engine, OpenSeesPy and PyNite agree, where the allowance for a
foreign solver's discretisation is an error bound computed from three refinement levels of that
solver ([docs/CROSSCHECK.md](../docs/CROSSCHECK.md)); tier D — equilibrium, energy and related
invariants, and the values at the result stations recomputed from statics and from an integration
of the member equations; tier E — equivalence checks of the modelling features; tier A — an
independent integration of the member differential equations, where the model is a straight member
line; N, K, C, E of the diagnosis route — null space of an independently assembled stiffness, the
engine's matrix, two foreign solvers under a probing load, invariance; tier F — coverage of every
design limit state used, required only where a task uses design results.

It does **not** mean that anyone reviewed, approved or certified a model, a result or a task. It
is a statement about agreement between computations made by this project's code and two
open-source solvers. The engine, the oracle and the translators to the foreign solvers share an
author and some code (fixed-end forces); an error common to all routes would not be detected.

## Checks run after the build

<!-- generated:checks -->
* `tools/check_dataset_scores.py` (normal scoring path, repairs and designs re-solved through the pipeline): the answer built from the hidden reference scores 1.0 on 748 of 748 tasks; the model as given, handed in as the repair, scores below the reference answer on 249 of 249 repair tasks (0 on 249; largest score 0.0). Findings: 0.
* `tools/audit_dataset.py` (not the pipeline that produced the references): 748 tasks on 86 base models; curve-balled models in both splits: 0; held-out tasks on a base model that train tasks also use: 0. 236 stable models solved again by the engine (1392 reference values recomputed, largest relative difference 0.0e+00) and by the independent tableau solver of `tests/independent`: 888 first-order combinations, largest relative difference 7.3e-12; 331 second-order combinations with the tableau refined (it takes the axial force constant over a piece of a member; 4, 16, 64, 256 pieces per member, raised while the difference to the engine is above 1e-3): 7 needed more than 16 pieces, the largest difference at the last count run is 1.2e-04, and between the last two counts run the difference falls by a factor of at least 6.7 (a finding below 4). Number of mechanisms of 27 unstable models from the null space of the tableau: 27 agree. Right and wrong answers through the scorer: 499 predict and diagnose tasks. First observation of 25 tasks probed for parts of the change record: 0 found. Findings: 0.

These are computations by this project's code. Nobody has reviewed, approved or certified a task, a model or a reference.
<!-- /generated:checks -->

## Split

`benchfea.bench.splits.is_held_out(family, curve_ball, seed)`, defined by the committed file
`src/benchfea/bench/splits.json`; how the held-out combinations are selected is in
docs/BENCHMARK.md section 7.

<!-- generated:split -->
| split | tasks | cells verified | seeds |
|---|---|---|---|
| train | 348 | 123 | 0 |
| held-out | 400 | 140 | 0 and 4 |

Held-out: 16 (family, curve ball) combinations with every seed, and of every other combination the seeds with `seed % 5` in [4].
<!-- /generated:split -->

The environment refuses held-out instances unless it is constructed with `split="held-out"`.
The split is by construction, not by secrecy: the generators are public code. Base models are not
shared between the splits (the audit line above counts the exceptions, which must be none).

## Known biases and gaps

<!-- generated:gaps -->
**Curve balls without any task (7 of 33)** — attempts of their cells, by stage, category and the tier named in the log:

| curve ball | cells | why (attempts of cells not served) |
|---|---|---|
| `constraint_wrong_slave` | 16 | curve-ball:not-applicable 42, verify:refused 6 |
| `design_bracing_points_changed` | 16 | verify:not-covered (tier F) 26, curve-ball:not-applicable 19, verify:not-representable (tier C) 3 |
| `design_k_changed` | 16 | verify:not-covered (tier F) 37, curve-ball:not-applicable 10, verify:not-representable (tier C) 1 |
| `design_lb_full_where_braced` | 16 | verify:not-covered (tier F) 33, curve-ball:not-applicable 13, verify:not-representable (tier C) 2 |
| `design_lb_short_where_unbraced` | 16 | curve-ball:not-applicable 43, verify:not-covered (tier F) 5 |
| `design_unbraced_length_changed` | 16 | verify:not-covered (tier F) 44, curve-ball:not-applicable 2, verify:not-representable (tier C) 2 |
| `offset_sign_reversed` | 16 | curve-ball:not-applicable 47, verify:not-representable (tier C) 1 |

**Curve balls with 1 to 3 tasks (3)**:

| curve ball | tasks | cells verified | why (attempts of cells not served) |
|---|---|---|---|
| `offset_removed` | 3 | 1 of 16 | curve-ball:not-applicable 45 |
| `offset_wrong_axis` | 3 | 1 of 16 | curve-ball:not-applicable 45 |
| `support_roller_direction` | 3 | 1 of 16 | curve-ball:not-applicable 45 |

**Design results.** Design tasks: 0. Predict tasks that ask for a design ratio: 0. Repair tasks with a design-ratio quantity: 0. Tasks of the 5 design-only curve balls (`design_*`): 0. Task types not served on verified cells because design results did not verify are counted in the quarantine summary.

**Second order.** 249 of 748 tasks are on a model with a second-order combination (33.3 %).

**Share of the applicable cells that are served, by family** (applicable = the last attempt was not "curve ball not applicable"):

| family | applicable cells | verified | share |
|---|---|---|---|
| `braced_frame` | 34 | 25 | 74 % |
| `continuous_beam` | 31 | 28 | 90 % |
| `diaphragm_frame` | 51 | 41 | 80 % |
| `multistorey_2d` | 40 | 34 | 85 % |
| `multistorey_3d` | 45 | 36 | 80 % |
| `pipe_rack` | 41 | 34 | 83 % |
| `portal_frame` | 41 | 33 | 80 % |
| `stair_frame` | 37 | 32 | 86 % |

**Held-out combinations** (16; 12 with a verified instance, 4 without — those are not served and were not replaced):

| family | curve ball | cells verified | tasks | cells not verified |
|---|---|---|---|---|
| `continuous_beam` | `release_mechanism` | 2 | 4 | - |
| `continuous_beam` | `release_added` | 2 | 6 | - |
| `portal_frame` | `sections_swapped` | 2 | 6 | - |
| `portal_frame` | `release_pinned_frame` | 2 | 5 | - |
| `multistorey_2d` | `design_unbraced_length_changed` | 0 | 0 | verify:not-covered 2 |
| `multistorey_2d` | `load_moved` | 2 | 6 | - |
| `multistorey_3d` | `section_heavier` | 2 | 5 | - |
| `multistorey_3d` | `support_fixed_to_pinned` | 2 | 6 | - |
| `braced_frame` | `load_moved` | 2 | 6 | - |
| `braced_frame` | `design_k_changed` | 0 | 0 | verify:not-covered 2 |
| `pipe_rack` | `member_rolled_90` | 2 | 6 | - |
| `pipe_rack` | `load_moved` | 2 | 6 | - |
| `stair_frame` | `support_roller_direction` | 0 | 0 | curve-ball:not-applicable 2 |
| `stair_frame` | `release_removed` | 2 | 6 | - |
| `diaphragm_frame` | `design_lb_full_where_braced` | 0 | 0 | verify:not-covered 2 |
| `diaphragm_frame` | `load_direction_changed` | 2 | 6 | - |
<!-- /generated:gaps -->

How to read this:

* **`curve-ball:not-applicable`** is most of every gap: a curve ball on constraints, offsets or
  roller supports has nothing to act on in most families.
* **Design-only curve balls and design tasks depend on tier F.** A design result is a verified
  reference only if every limit state it uses has a hand calculation and a reproduced strength of
  the AISC Design Examples in the test suite. Most generated models use at least one limit state
  without that (the registry `src/benchfea/verify/design_coverage.json` lists which), so the
  design-only curve balls, design-ratio quantities and design tasks are rare or absent. The code
  paths are tested with a device that counts every item as covered (`tests/bench_helpers.py`),
  which says nothing about any design value.
* **Tier C decides most of the rest.** A model that only one foreign solver can represent with an
  eligible error bound is not cross-checked and not served: second order together with rigid
  offsets or with rigid links or diaphragms; crossing truss braces under self weight in second
  order (PyNite fails on them, see QUARANTINE_NOTES.md); a node lying on another member's axis;
  any buckling analysis. The served tasks therefore lean towards models without those features.
* **No buckling factors** are asked for, and **`invalid-model` curve balls are not served** (there
  is no analysis to verify), nor are mechanisms that exist only through slack tension-only members.
* **Seeds.** With few seeds the train split has fewer instances than the held-out split (held-out
  combinations take every seed), and a (family, curve ball) pair has at most one instance per seed.
* **Synthetic models, one change per model.** Eight parametric families, generic sections, nominal
  steel, rounded loads.

## Quarantine log

<!-- generated:quarantine -->
1138 entries: 834 attempts that were skipped and 304 task types not served on a verified cell. Every entry carries its reason text.

| attempts skipped: stage:category | attempts | of those replaced by a later attempt | tier named in the reason | meaning |
|---|---|---|---|---|
| `curve-ball:not-applicable` | 645 | 25 | - | the curve ball does not fit the model (precondition not met, or no measured material effect); not a quarantine |
| `verify:not-covered` | 145 | 0 | tier F 145 | tier F: a design limit state used by the reference lacks coverage |
| `verify:not-representable` | 33 | 11 | tier C 33 | fewer than two foreign solvers can represent the model or give an eligible error bound, so it is not cross-checked |
| `verify:refused` | 6 | 0 | - | the curve-balled model is refused by validation (no analysis to verify) |
| `verify:disagreement` | 5 | 3 | tier C 2, tier K 2, tier D 1 | two computations of the same quantity differ by more than the tolerance of the tier |

| task type not served on a verified cell: reason | cells |
|---|---|
| design: the curve-balled model is a mechanism | 27 |
| design: the design results of the given model are not verified (with the design results, tier F) | 236 |
| predict: the curve-balled model is a mechanism | 27 |
| repair: the curve ball moves no significant result of the intended design beyond the zero-credit e | 10 |
| repair: the reference of the intended model is not verified | 4 |
<!-- /generated:quarantine -->

A quarantined attempt is left as it is: nothing is served from it, and no tolerance is changed to
admit it. The log (`quarantineLog` in the manifest) holds the reason text of every entry; written
diagnoses of the attempts logged as `disagreement` are in [QUARANTINE_NOTES.md](QUARANTINE_NOTES.md).

## Intended use

* Evaluating and training agents (language models with tools, scripted agents) on reading,
  checking and correcting structural analysis models; measuring whether an agent notices a
  modelling error that matters.
* The `train` split for development and training, `heldout` for reporting.

## Not intended for

* Design or assessment of real structures. The models are synthetic and the design checks cover a
  subset of AISC 360-22 member limit states.
* Claims about engineering competence in general: the tasks cover elastic first- and second-order
  frame analysis of eight parametric families and one change per model.
* A statement that any model or result is correct beyond what the tiers listed establish.
* Training on the held-out split, or tuning on it, and then reporting on it.

## Other limits of the tasks

* **Predict in mode tool-allowed is a tool-use task**: the values can be read off the analysis. The
  curve ball matters there only in that the model must be analysed as given. Mode no-solver asks for
  the reasoning; its tolerances are wide (5 % / 25 %) and a large frame is not a hand calculation.
* **Diagnose compares a description with a file.** The intended design is given completely, in
  words; finding the change is careful reading plus the judgement of what matters. The direction
  vocabulary is mechanical (docs/BENCHMARK.md section 2.2), not an engineering judgement.
* **Entities are matched by id, not by kind**, and recall is the best over the kinds of the change
  record: naming all members touched, or all nodes touched, is a complete answer.
* **Repair is judged on named results, not on the edit.** Any model that is stable, passes the
  pipeline and reproduces the named quantities of the intended design within tolerance scores
  full marks — including a model that differs from the intended one elsewhere. A correct repair
  that uses a feature the pipeline cannot verify scores 0 with that reason: an unverified result is
  never scored as correct.
* **The reference weight of a design task is the best found by the stated search, not a proven
  optimum.**
* **Scoring a repair or design runs the verification pipeline** (seconds to tens of seconds; budget
  300 s; what the pipeline could not run within the budget counts as not passed).
* **The held-out split is a convention.** Generators and curve balls are public code; the
  environment refuses held-out instances unless constructed for them, nothing more.

## How to rebuild

The commands are in "How it was generated". The build is deterministic and resumable (`--cache`);
with a fresh cache it takes hours on four cores. After a build run the two checks, then
`python3 tools/build_dataset_docs.py`, which rewrites this card, section 8 of docs/BENCHMARK.md and
`servedCombinations` in `src/benchfea/bench/splits.json`. A change of the generator, a curve ball,
the task code or the verification pipeline changes the dataset; `tests/test_dataset.py` rebuilds one
committed task from its key and fails when the data no longer match the code.

## Licence

MIT, the same as the code (the `LICENSE` file at the root of the repository). The dataset is
generated by the code of this repository and contains no third-party data: no table, text or figure
of AISC 360-22, of the AISC Design Examples or of any other reference is included, and no
proprietary section database was used (sections are a generic catalogue defined in the generator).
