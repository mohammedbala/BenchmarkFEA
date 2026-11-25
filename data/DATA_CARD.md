# Data card — BenchFEA tasks

## What it is

534 benchmark tasks for agents that work on structural frame models, each with a hidden
reference. A task hands the agent a model of a steel frame (JSON input file of the `benchfea`
engine, format 1.0.0) in which one thing was changed in a way that matters (a *curve ball*), and
asks it to **predict** named results, to **diagnose** what was changed and where (or, for an
unstable model, the mechanism), to **repair** the model, or (2 tasks) to **design** its member
sizes. Task types, answer formats and scoring:
[docs/BENCHMARK.md](../docs/BENCHMARK.md).

```
data/tasks/manifest.json                       counts, every task, every cell, the quarantine log, versions
data/tasks/{train,heldout}/<id>.task.json      public: statement, answer format, quantities / criteria
data/tasks/{train,heldout}/<id>.model.json     public: the input file of the task
data/answers/{train,heldout}/<id>.answers.json hidden: change record, reference solution, verification record
data/user/                                     not part of the dataset (user-supplied files, never committed)
```

Hand `data/tasks/` to an agent; keep `data/answers/` away from it. Task ids are opaque (they do not
name the curve ball).

## How it was generated

`benchfea bench make-dataset --out data --seeds 0,4 --workers 4 --cache <dir> --clean` with the code of this repository (benchfea 0.1.0, python 3.13.18, numpy 2.5.3, scipy 1.18.1, foreign:opensees 3.8.0.0, foreign:pynite 3.2.1, task 1). No hand editing.

1. **Model**: `benchfea.bench.generate(family, seed, difficulty)` — a seeded parametric generator
   of eight families (`braced_frame`, `continuous_beam`, `diaphragm_frame`, `multistorey_2d`, `multistorey_3d`, `pipe_rack`, `portal_frame`, `stair_frame`) with member sizes found by the engine's own sizing loop
   ([docs/DIFFICULTY.md](../docs/DIFFICULTY.md)). The models are synthetic: generic section
   catalogue, nominal steel properties, rounded loads. They are not real structures.
2. **Curve ball**: `benchfea.bench.curveballs.apply(model, name, seed)` — one of
   33 named changes (supports, releases, offsets, constraints, sections, loads, units,
   second-order setting, design data), accepted only if its *measured* effect on the analysis or
   the design is material ([docs/CURVEBALLS.md](../docs/CURVEBALLS.md)). The exact change is
   recorded (JSON-pointer edits with before / after values).
3. **Verification**: the curve-balled model, and for repair tasks the model before the change, go
   through `benchfea.verify.verify_model` ([docs/VERIFICATION.md](../docs/VERIFICATION.md)).
   Only a reference that comes out `verified` is used. An instance that does not is skipped, the
   skip is logged, and the next seed of a documented chain is tried (at most three attempts); if
   none verifies, nothing is written for that cell. No tolerance was adjusted and no case was
   special-cased to make it pass.
4. **Tasks**: every task type that can be served on a verified instance becomes one task; the mode
   (`tool-allowed` / `no-solver`) is a fixed pseudo-random choice.

Cells built: every family × every curve ball × seeds 0 and 4 = 528 cells, of which
192 gave a verified instance. Everything that was not served is in the
`quarantineLog` of the manifest with its reason (1286 entries); the counts by family and by
curve ball are in section 8 of docs/BENCHMARK.md.

Checked after the build: the answer built from the hidden reference scores 1.0 on all 534 tasks
through the normal scoring path, and the unrepaired model scores 0 on all 174 repair tasks.

Also run after the build, with means that are not the pipeline that produced the references
(`python3 tools/audit_dataset.py data`; the numbers below are its output on this dataset): every
stable model (166; all of their combinations are first-order analyses) solved by the independent
tableau solver of `tests/independent` — 823 combinations, largest relative difference to the engine
5.6·10⁻¹² on node displacements, reactions, end forces and station values; the number of mechanisms
of every unstable model (26) confirmed by the null space of the tableau; every task in the split
its (family, curve ball, seed) belongs to; a wrong answer scores below 1.0 on every predict and
diagnose task; nothing of the change record in the first observation of 25 probed tasks. The tool
is extra evidence, not a tier, and what it finds is computation, not review.

## What "verified" means

**Only this: the reference passed every verification tier that applies to it, at the tolerances
stated in its record.** The record of each task (in its answers file) lists the tiers passed, the
tolerances and the largest differences found. For the tasks here:

* analysis results — tier C (the engine, OpenSeesPy and PyNite agree: 10⁻⁶ relative plus the
  estimated discretisation error of an approximate foreign solve), tier D (equilibrium, energy and
  related invariants; the values at the result stations recomputed from statics and from an
  integration of the member equations), tier E (equivalence checks of the modelling features); tier A (an
  independent integration of the member differential equations) where the model is a straight
  member line;
* mechanism diagnoses — tiers N, K, C, E of the diagnosis route (null space of an independently
  assembled stiffness, the engine's matrix, two foreign solvers under a probing load, invariance);
* design results (only where a task uses them: 2 design tasks, 5 predict tasks that ask for a
  design ratio, 4 repair tasks with a design criterion) — in addition tier F: every design limit
  state used has a hand calculation and a reproduced value of the AISC Design Examples in the
  test suite.

It does **not** mean that anyone reviewed, approved or certified a model, a result or a task. It
is a statement about agreement between computations made by this project's code and two
open-source solvers. The engine, the oracle and the translators to the foreign solvers share an
author and some code (fixed-end forces); an error common to all routes would not be detected.

## Split

`benchfea.bench.splits.is_held_out(family, curve_ball, seed)`, defined by the committed file
`src/benchfea/bench/splits.json`:

* 16 (family, curve ball) combinations are held out entirely (two per family);
* of every other combination the seeds with `seed % 5 == 4` are held out.

| split | tasks |
|---|---|
| train | 248 |
| held-out | 286 |

The environment refuses held-out instances unless it is constructed with `split="held-out"`.
The split is by construction, not by secrecy: the generators are public code.

The split separates (family, curve ball, seed) instances, not base models: the 12 held-out instances
with seed 0 (held-out combinations; 34 tasks: 12 diagnose, 12 predict, 10 repair) are built on base
models that train tasks of the same family, seed and difficulty also use, and for the 10 repair
tasks among them the reference answer (the intended model) is the reference answer of train repair
tasks as well. No curve-balled model occurs in both splits. The 246 held-out tasks with seed 4
share no base model with the train split.

## Intended use

* Evaluating and training agents (language models with tools, scripted agents) on reading,
  checking and correcting structural analysis models; measuring whether an agent notices a
  modelling error that matters.
* The `train` split for development and training, `heldout` for reporting.

## Not intended for

* Design or assessment of real structures. The models are synthetic and the design checks cover a
  subset of AISC 360-22 member limit states.
* Claims about engineering competence in general: the tasks cover linear-elastic (and some
  second-order elastic) frame analysis of eight parametric families and one change per model.

## Limits

* **What "verified" covers.** Agreement of computations at the stated tolerances (analysis results:
  tiers C, D, E and, for straight member lines, A; mechanism diagnoses: tiers N, K, C, E). The
  engine, the ODE oracle and the translators to the two foreign solvers come from one project and
  share the fixed-end-force code; an error common to all routes would not be detected. Nobody has
  reviewed, approved or certified any task or reference.
* **Design results are almost never a verified reference.** Tier F requires, for every limit state
  a model uses, a hand calculation and a reproduced value of the AISC Design Examples in the test
  suite. Most generated models use limit states without full coverage (measured most often: Eq.
  F2-7, Eq. E4-1 / E4-2, Section H3 and the engine's torsion screening rule, Table B4.1a Case 8 /
  B4.1b Case 21). Consequences, all visible in the quarantine log: the five design-only curve balls
  (`design_*`) are not served at all; a `designRatio` quantity, a design criterion in a repair task
  and a design task exist only on the few models whose design results did verify (in this dataset:
  5 instances, all continuous beams — 2 design tasks, 5 predict tasks with a design ratio, 4 repair
  tasks with the design criterion); the governing design ratio is otherwise never asked for. The code paths are tested with a test device that counts every item as covered
  (`tests/bench_helpers.py`), which says nothing about any design value.
* **Second-order models are mostly quarantined.** Tier C compares the engine (exact stability
  functions) with foreign solvers whose error estimate is a leading-order value, not a bound
  (VERIFICATION.md, "Second-order cross-checks"); the difference often lands just above the allowed
  one and the case is quarantined as a `disagreement` (in this dataset: 261 of the 278 attempts
  skipped with category `disagreement` are tier-C differences of this kind, in the families with
  second-order combinations). Nothing was loosened. Of the cells whose last attempt was not
  "curve ball not applicable", 15 of 44 are served for `multistorey_2d`, 19 of 43 for
  `multistorey_3d`, 18 of 40 for `pipe_rack` and 19 of 32 for `braced_frame`, against 28 of 33 for
  `continuous_beam` and 30 of 38 for `stair_frame`: the served tasks lean towards first-order
  models.
* **Other quarantines seen in this build** (each left as it is, with its diagnosis in the log):
  the superposition check of tier D not passed (5 attempts): 1.58e-9 against 1.0e-9 once
  (`multistorey_2d`) and 0.64 to 2.33 four times (`stair_frame`, `support_removed`). Diagnosis of the
  latter (audit): the check takes its ratio array by array, and the column whose support was removed
  carries no force in a combination without self weight, so the ratio is one of two round-off numbers
  (1.0e-9 N against 4.6e-10 N where the largest member force is 1e7 N); the combination and the sum of
  its load cases agree to 6e-12 in displacements and reactions. A defect of the check on the
  conservative side: the attempts stay quarantined and the check was not changed; a mechanism made of a loaded node without members
  is not confirmed by the diagnosis tiers (12 attempts, `diaphragm_frame`, `constraint_removed` /
  `constraint_wrong_master`: never served); models that fewer than
  two foreign solvers can represent (8 attempts: second order with rigid offsets, a node lying on a
  member's axis). Two attempts (`continuous_beam`) that an earlier build of this dataset logged as a tier-A
  disagreement are not one: a regular result station one ulp away from a point load gives three stations
  at one point, and the comparison took the limit from the wrong side for the third; the comparison was
  corrected (it accepts either limit there, and tier D requires that the stations of the point show
  both), the analysis of that model verifies, and the two cells are now logged as `not-covered` (tier F). The generator raised an error for three `multistorey_2d` level-2 seeds of the
  chain (an empty design group; 24 attempts): logged as `generator-error`, the next seed is used.
* **No buckling factors.** A model that asks for a buckling analysis cannot pass tier C (no foreign
  solver reproduces it), so `bucklingFactor` is never asked for.
* **`invalid-model` curve balls are not served** (the model is refused by validation; there is no
  analysis the pipeline could verify), nor are mechanisms that the diagnosis tiers do not confirm
  (a loaded node without members) or that exist only through slack tension-only members.
* **Predict in mode tool-allowed is a tool-use task**: the values can be read off the analysis. The
  curve ball matters there only in that the model must be analysed as given. Mode no-solver asks for
  the reasoning; its tolerances are wide (5 % / 25 %) and a large frame is not a hand calculation.
* **Diagnose compares a description with a file.** The intended design is given completely, in
  words; finding the change is careful reading plus the judgement of what matters. One change per
  model; the direction vocabulary is mechanical (section 2.2), not an engineering judgement.
* **Entities are matched by id, not by kind**, and recall is the best over the kinds of the change
  record: naming all members touched, or all nodes touched, is a complete answer.
* **Repair is judged on named results, not on the edit.** Any model that is stable, passes the
  pipeline and reproduces the named quantities of the intended design within tolerance scores
  full marks — including a model that differs from the intended one elsewhere. A correct repair
  that uses a feature the pipeline cannot verify (fewer than two independent solvers can represent
  it) scores 0 with that reason: an unverified result is never scored as correct.
* **The reference weight of a design task is the best found by the stated search, not a proven
  optimum**; the search moves members group-wise and the catalogue is the generic one of the
  generators. A feasible submission that uses a limit state without tier-F coverage fails the gate.
* **Scoring a repair or design runs the verification pipeline** (seconds to tens of seconds; budget
  300 s; what the pipeline could not run within the budget counts as not passed).
* **The held-out split is a convention.** Generators and curve balls are public code; the
  environment refuses held-out instances unless constructed for them, nothing more.
* **Synthetic models.** Generic sections, nominal steel, rounded loads, eight parametric families;
  not real structures, and not a basis for design.

## Licence

The dataset is generated by the code of this repository and is distributed under the same terms as
the repository (see the repository's licence file; if the repository has none yet, no licence has
been granted for the data either). It contains no third-party data: no table, text or figure of
AISC 360-22, of the AISC Design Examples or of any other reference is included, and no proprietary
section database was used (sections are a generic catalogue defined in the generator).
