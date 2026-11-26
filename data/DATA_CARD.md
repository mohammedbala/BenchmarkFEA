# Data card — BenchFEA tasks

Every number on this card stands between `generated` markers and is written by
`python3 tools/build_dataset_docs.py` from `data/tasks/manifest.json`, the task files and
`data/CHECKS.json`; `--check` (run by the test suite) fails when the card and the data disagree.
The text outside the markers is prose and holds no count.

## What it is

<!-- generated:headline -->
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
data/user/                                     not part of the dataset (user-supplied files, never committed)
```

Hand the task and model files to an agent; keep `data/answers/` away from it. Task ids are opaque
(they do not name the curve ball). The manifest lists the curve ball of every cell next to its task
ids, so it belongs on the answers' side of that line too; the environment (`BenchFEA-v0`) does not
show it.

## How it was generated

<!-- generated:build -->
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
<!-- /generated:composition -->

By curve ball (a cell is one family × curve ball × seed; it is verified when one attempt of its
chain gave a verified reference):

<!-- generated:curve-balls -->
<!-- /generated:curve-balls -->

## What "verified" means

**Only this: the reference passed every verification tier that applies to it, at the tolerances
stated in its record.** The record of each task (in its answers file) lists the tiers passed, the
tolerances and the largest differences found.

<!-- generated:verified -->
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
<!-- /generated:checks -->

## Split

`benchfea.bench.splits.is_held_out(family, curve_ball, seed)`, defined by the committed file
`src/benchfea/bench/splits.json`; how the held-out combinations are selected is in
docs/BENCHMARK.md section 7.

<!-- generated:split -->
<!-- /generated:split -->

The environment refuses held-out instances unless it is constructed with `split="held-out"`.
The split is by construction, not by secrecy: the generators are public code. Base models are not
shared between the splits (the audit line above counts the exceptions, which must be none).

## Known biases and gaps

<!-- generated:gaps -->
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
  offsets, rigid links or diaphragms, or shear deformation; members close to buckling; partial
  releases in PyNite; tension-only members; any buckling analysis. The served tasks therefore lean
  towards models without those features.
* **No buckling factors** are asked for, and **`invalid-model` curve balls are not served** (there
  is no analysis to verify), nor are mechanisms that exist only through slack tension-only members.
* **Seeds.** With few seeds the train split has fewer instances than the held-out split (held-out
  combinations take every seed), and a (family, curve ball) pair has at most one instance per seed.
* **Synthetic models, one change per model.** Eight parametric families, generic sections, nominal
  steel, rounded loads.

## Quarantine log

<!-- generated:quarantine -->
<!-- /generated:quarantine -->

A quarantined attempt is left as it is: nothing is served from it, and no tolerance is changed to
admit it. The log (`quarantineLog` in the manifest) holds the reason text of every entry.

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
