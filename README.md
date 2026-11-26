# BenchmarkFEA

## 1. What it is

BenchmarkFEA is a program that analyses three-dimensional frames of beams, columns and braces
(first order, second order and linear buckling) and checks steel members against parts of
AISC 360-22. Around it is a benchmark for AI agents that work on structural models: an agent is
given a frame model in which one thing has been changed in a way that matters (a "curve ball"), and
its answer is scored against a reference solution. A reference solution is served only if it has
passed the verification tiers described in section 7.

## 2. What it does not do

Read [LIMITS.md](LIMITS.md) before you rely on any number. In short:

* It is not a design office tool. It has not been used on a real project, and nobody has checked,
  approved or certified it or any result it gives.
* Only straight prismatic members: no plates, shells, cables, tapered members, dynamics, plasticity
  or large displacements ([LIMITS.md, section 1](LIMITS.md#1-analysis)).
* The design checks cover part of AISC 360-22 for doubly symmetric sections only, and no
  connections ([LIMITS.md, section 3](LIMITS.md#3-design-checks)).
* Several kinds of result can be computed but cannot reach the status "verified", among them every
  buckling factor and nearly all design results
  ([LIMITS.md, section 4.3](LIMITS.md#43-what-cannot-currently-be-verified)).
* What it must not be used for: [LIMITS.md, section 7](LIMITS.md#7-do-not-use-this-software-for).

## 3. Install

You need Python 3.11 or newer and Git. Type the commands in a terminal (on Windows: PowerShell),
one line at a time.

```bash
git clone https://github.com/mohammedbala/BenchmarkFEA.git
cd BenchmarkFEA
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[env,crosscheck]"
```

* On Windows the fourth line is `.venv\Scripts\activate`, and `python3` may be called `python`.
* `source .venv/bin/activate` has to be repeated in every new terminal window before `benchfea`
  is found.
* `[env,crosscheck]` adds three optional packages: Gymnasium (for the benchmark environment) and
  OpenSeesPy and PyNite (two other analysis programs, used only to cross-check results, section 7).
  If they do not install on your computer, use `pip install -e .`: analysis, design checks, curve
  balls and scoring of predictions still work, but nothing can be verified.

Check the installation:

```bash
benchfea --version
benchfea validate examples/01_simple_beam.json
```

The second command must end with `valid (0 error(s), 0 warning(s))`.

## 4. Run one example and read the result

[`examples/01_simple_beam.json`](examples/01_simple_beam.json) is an 8 m simply supported floor
beam with a dead load of 12 kN/m and a live load of 18 kN/m. The file is plain text; open it in any
editor. [`examples/01_simple_beam.md`](examples/01_simple_beam.md) explains it line by line, and
[docs/INPUT_REFERENCE.md](docs/INPUT_REFERENCE.md) lists every field an input file may hold.

```bash
mkdir out
benchfea run examples/01_simple_beam.json -o out/beam.results.json
```

The program prints one line per load combination (three here), then the number of defaults it
used and the name of the results file. The line of the strength combination `ULS`
(1.2 × dead + 1.6 × live = 43.2 kN/m), wrapped here to fit the page:

```
ULS [firstOrder]: largest node displacement 0 m (ux at A); largest member-axis displacement
component 0.0344387 m (member B1, x = 4); sum of reactions (0, 0, 345.6) kN; invariants passed
(largest: rigid-body modes 3.80e-17)
```

| What is printed | What it means |
|---|---|
| `[firstOrder]` | linear elastic analysis; no second-order effects |
| `largest node displacement 0 m` | both nodes of this beam are supports, so no node moves |
| `largest member-axis displacement component 0.0344387 m (member B1, x = 4)` | the largest deflection along the member: 34.4 mm at 4 m from its start, which is midspan. The summary gives the size only; the sign is in the results file |
| `sum of reactions (0, 0, 345.6) kN` | the supports together push up with 345.6 kN = 43.2 kN/m × 8 m. The three numbers are the global X, Y and Z components; Z points up in this file. A reaction is the force of the support **on the structure** |
| `invariants passed` | the program's own checks of this solution (equilibrium of the whole structure, of every joint and of every member, and others) are satisfied |
| `(largest: rigid-body modes 3.80e-17)` | the name and size of the check that came closest to its tolerance. It is rounding error, and the last digits may differ on your computer |

All results are in `out/beam.results.json`, in the units of the input file (here kN and m). It is
a plain text file; [docs/OUTPUT_REFERENCE.md](docs/OUTPUT_REFERENCE.md) describes every key. Under
`combinations`, then `ULS`:

| Where in the file | Value | Meaning |
|---|---|---|
| `reactions` → `A` → `Fz` | 172.8 kN | vertical reaction at A, upward (= wL/2) |
| `members` → `B1` → `stations` → `Mz`, eleventh of 21 values | 345.6 kN·m | bending moment at midspan (= wL²/8). Positive `Mz` is sagging here: tension on the lower side |
| `members` → `B1` → `stations` → `Vy`, first value | −172.8 kN | shear force at the start of the member |
| `members` → `B1` → `stations` → `uy`, eleventh value | −0.0344 m | deflection at midspan. The member's local y axis points up, so a negative value is downward |

The sign rules for loads, reactions, member end forces and internal forces, and the directions of
the member axes, are in [docs/CONVENTIONS.md](docs/CONVENTIONS.md), with figures.

The steel member checks:

```bash
benchfea design examples/01_simple_beam.json -o out/beam.design.json
```

The program prints a table with one row per member (shortened here: `...` stands for the rest of
the line), followed by a warning that the combination `ULS` was analysed first order:

```
member  status         ratio  comb.     governing limit state [clause]
B1      conditional    0.717  SLS       deflection limit 2 (local y) [AISC 360-22 L1, L2 ...]
```

* `ratio 0.717` is demand divided by capacity for the check that governs. Here it is the
  deflection limit L/240 that the input file sets for the combination `SLS`: 23.9 mm of 33.3 mm.
  The limit is the user's; AISC 360-22 gives no number.
* `conditional` means every check passes, but the strength combination was analysed first order,
  and AISC 360-22 Chapter C asks for second-order effects in the required strengths. The program
  says so and does not amplify anything by itself.
* Every check, with the clause, the equations used and the values that went in, is in
  `out/beam.design.json` ([docs/OUTPUT_REFERENCE.md, section 2](docs/OUTPUT_REFERENCE.md#2-design-document)).

## 5. Apply one curve ball and see what changed

A curve ball is a named change to a model. `benchfea perturb --list` shows all 33;
[docs/CURVEBALLS.md](docs/CURVEBALLS.md) describes each. This one turns the pinned support at A
into a fixed support:

```bash
benchfea perturb examples/01_simple_beam.json support_pinned_to_fixed -o out/beam_fixed.json --record out/beam_fixed.record.json
benchfea run out/beam_fixed.json -o out/beam_fixed.results.json
```

The first command prints two lines (wrapped here). The first is:

```
support_pinned_to_fixed seed 0: support A: rotation(s) ry, rz fixed; effect class changes-analysis;
analysis change 100.0%, design ratio change 7.6%; 2 change(s): /supports/A/dofs, /supports/A/type
```

The second line says which verification tiers the changed model is expected to be able to pass.
It is an expectation, not the result of a verification run.

The change record `out/beam_fixed.record.json` says, in words:

* **What was changed.** Two entries of the input file: `/supports/A/type` went from `pinned` to
  `fixed`, and `/supports/A/dofs` (which held the twist of the beam at A) was removed because a
  fixed support holds every rotation already. The record keeps the old values, so the change can be
  undone exactly.
* **Why it matters.** Full fixity is rarely there unless the foundation is designed for it; a
  support modelled as fixed attracts a moment and relieves the rest of the structure.
* **What it did to the results, measured by running both models.** `analysis change 100.0%`: the
  largest change of any result, relative to the size of that kind of result, is 100 % — a moment
  where there was none. `design ratio change 7.6%`: the governing ratio went from 0.717 to 0.663.
  A change counts as a curve ball only when a measured effect is at least 5 %, or when it changes
  whether a member passes its checks; a change that makes no difference is refused.

Compare the two runs (combination `ULS`):

| | simply supported (section 4) | support A fixed |
|---|---|---|
| vertical reaction at A | 172.8 kN | 216.0 kN (= 5wL/8) |
| moment reaction at A (`My`) | 0 | −345.6 kN·m (size wL²/8) |
| largest deflection | 34.4 mm at 4.0 m | 14.3 mm at 4.8 m |

The beam has become a propped cantilever. The model still runs and its output looks as tidy as
before; that is the point of the exercise.

## 6. Score one agent answer

The benchmark's tasks are in [`data/tasks/`](data/DATA_CARD.md). A task of type "predict" names
some results and asks for their values. The file
[`examples/answers/te6bac428ea3e9b.answer.json`](examples/answers/te6bac428ea3e9b.answer.json) is
an answer to one task of the training set; five of its six numbers were read from an analysis of
the task's model and the sixth, `q6`, was given the wrong sign on purpose.

```bash
benchfea score examples/answers/te6bac428ea3e9b.answer.json --task data/tasks/train/te6bac428ea3e9b.task.json
```

It prints (`...` stands for the lines of `q2` to `q5`, which look like that of `q1`):

```
task te6bac428ea3e9b [predict, tool-allowed]: score 0.8333
  5 of 6 quantities within tolerance; score 0.833
  q1: credit 1.000 (relative error 4.528e-05; full credit up to 1%)
  ...
  q6: credit 0.000 (relative error 2.000e+00; full credit up to 1%)
```

The score is between 0 and 1. For this kind of task each quantity earns full credit up to 1 %
error and none from 5 %; the score is the mean. The reference values are in the hidden file
`data/answers/train/te6bac428ea3e9b.answers.json`, which an agent must not see. The other task
types (diagnose what was changed, repair the model) and their scoring are in
[docs/BENCHMARK.md](docs/BENCHMARK.md).

## 7. How verification works, and what "verified" means

**"Verified" means only: the reference solution passed the tiers listed in its record, at the
tolerances stated there. Nobody has checked, approved or certified any result.** It is a statement
that independent computations agree. It does not say that a model represents a real structure well
or that a design is adequate.

The tiers ([docs/VERIFICATION.md](docs/VERIFICATION.md)):

| Tier | What is compared |
|---|---|
| A | the engine against a separate program that integrates the differential equations of a straight member line directly, and against values derived in this repository |
| B | the engine against numbers printed in a publication, each cited |
| C | the engine against two other analysis programs (OpenSeesPy and PyNite) on the same model |
| D | the solution against itself: equilibrium, compatibility, energy, superposition, reciprocity, and the values between the member ends recomputed from statics |
| E | the model against an equivalent model that reaches the same answer by another route (members split, rigid parts replaced by stiff members, and so on) |
| F | for design results: every limit state used must have both a hand calculation from AISC 360-22 and a reproduced value of the AISC Design Examples in the test suite |

Try it (this needs OpenSeesPy and PyNite):

```bash
benchfea verify examples/01_simple_beam.json
```

**Tolerances.** There is one written policy for all models and problems (`benchfea.verify.POLICY`,
printed in full at the end of the report). Results of linear analyses must agree to 10⁻⁶. A
published number must be met to half a unit of its last printed digit, or to the tolerance the
source itself states. In second-order analyses the two other programs divide each member into
pieces, so each of their results carries an error bound computed from three levels of subdivision;
the comparison allows 10⁻⁶ plus these bounds, and a program whose bound is above 0.5 % does not
count. No tolerance is changed for a single case.

**Quarantine.** Every tier that applies must pass; there is no majority vote. If one tier is not
passed — two computations differ, fewer than two other programs can represent the model, or a check
could not be run — the result is *quarantined* with a written reason and is not served as a
reference. A published problem that the engine does not reproduce stays in the library as
quarantined, with a diagnosis; it is not dropped and its tolerance is not loosened.

**Current numbers**, from [docs/VERIFICATION_REPORT.md](docs/VERIFICATION_REPORT.md), which is
generated from actual runs:

<!-- from:VERIFICATION_REPORT -->
* Problem library: 173 problem(s) — 109 passed, 0 failed, 64 quarantined, 0 invalid.
* Of the 109 passed problem(s), 74 have a model whose record is `verified` (every applicable tier passed, including the cross-check by two foreign solvers).
* Model bank: 50 model(s) — 36 verified, 13 quarantined, 1 refused by the engine.
<!-- /from:VERIFICATION_REPORT -->

"Passed" is not always the reproduction of a printed number: for some problems only the statement
is published and the target values were computed in this repository by a separate program. The
report counts these separately, and [LIMITS.md, section 4.4](LIMITS.md#44-the-published-problems)
says what a sceptical reader should know about them.

## 8. Adding a model, a curve ball or a verification case

* **A model** is one input file. Start from an example, check it with `benchfea validate`, and see
  [CONTRIBUTING.md](CONTRIBUTING.md#a-an-example-model).
* **A curve ball** is a registered, seeded change whose effect is measured:
  [CONTRIBUTING.md](CONTRIBUTING.md#c-a-curve-ball).
* **A verification case** is a folder with a model and a `problem.json` that cites its source and
  lists the target values: [CONTRIBUTING.md](CONTRIBUTING.md#d-a-verification-problem).

The rules that apply to all three: no value from memory, no value from a source that was not
opened, no tolerance tuned to make a case pass. Every document of the repository is listed, with
one line each, in [docs/README.md](docs/README.md).

## What you must supply yourself

Nothing below is in this repository.

* **AISC 360-22** (*Specification for Structural Steel Buildings*, ANSI/AISC 360-22). The design
  checks cite its sections and equations; you need your own copy to read what is cited.
* **The AISC Shapes Database**, only if you want to name sections as standard shapes such as
  `W14X48`. Without it everything else works: give sections by their dimensions, as the
  examples do (a section given only by its properties can be analysed but not design-checked).
  * Where to put it: `data/user/aisc-shapes-database-v16.0.xlsx` (or `.csv`), or set the environment
    variable `BENCHFEA_SHAPES_DB` to the path of the file. An `.xlsx` file also needs
    `pip install -e ".[shapes]"`.
  * Expected format: the first row of the first worksheet (or of the CSV file) holds the column
    names `AISC_Manual_Label` or `EDI_Std_Nomenclature`, `Type`, `A`, `d`, `bf`, `tw`, `tf`, `Ht`,
    `B`, `tdes`, `OD`, `Ix`, `Zx`, `Sx`, `Iy`, `Zy`, `Sy`, `J`, `Cw` (and optionally `kdes`), with
    lengths in inches. This layout is an assumption: the database could not be opened while the
    reader was written. If your file differs, see [`data/user/README.md`](data/user/README.md).
  * The database is published by AISC. Check its licence terms yourself before you pass it on to
    anyone. Everything in `data/user/` except its README is ignored by Git, so the file is not
    committed by accident.
* **The AISC Steel Construction Manual and the AISC Design Examples** are not included. Some tests
  reproduce printed values of the Design Examples; the values are cited, the text is not copied.
* **OpenSeesPy and PyNite** (optional) are needed for tier C. They are installed by
  `pip install -e ".[crosscheck]"`. Without them no result can be verified.

## The Gymnasium environment

For those who program: the benchmark is also a [Gymnasium](https://gymnasium.farama.org/)
environment. Save the following as `try_env.py` in the repository folder and run
`python try_env.py`:

```python
import json
import gymnasium
import benchfea.bench.env  # registers "BenchFEA-v0"

env = gymnasium.make("BenchFEA-v0", split="train", dataset_dir="data")
obs, info = env.reset(seed=1, options={"task_type": "predict", "mode": "tool-allowed"})
print(obs["statement"])                                   # the task in words
obs, reward, terminated, truncated, info = env.step(json.dumps({"tool": "run_analysis", "args": {}}))
print(obs["tool_result"][:300])                           # what the analysis tool returned
obs, reward, terminated, truncated, info = env.step(json.dumps({"tool": "submit", "args": {"answer": {"values": {"q1": 0.0}}}}))
print("reward:", reward, "| episode ended:", terminated)  # the score of the submitted answer, 0 to 1
```

The tools an agent may call, the observation and the reward are described in
[docs/BENCHMARK.md, section 5](docs/BENCHMARK.md#5-the-environment); example agents are in
[`examples/agents/`](examples/agents/README.md).

## The dataset

<!-- from:DATA_CARD -->
**748 verified tasks** from 263 verified cells of 528 (every family × every curve ball × seeds 0 and 4): 236 predict, 263 diagnose, 249 repair, 0 design; 348 train, 400 held-out.
<!-- /from:DATA_CARD -->

The tasks are in `data/tasks/`, the hidden references in `data/answers/`. How the dataset was
built, what it contains, its known gaps (there is no design task) and how to rebuild it:
[data/DATA_CARD.md](data/DATA_CARD.md). The models are synthetic; they are not real structures.

## Licence

MIT; see [LICENSE](LICENSE). The licence covers the code and the generated dataset. It does not
cover the references named above, which are not in the repository.

## How to cite

> BenchmarkFEA contributors. *BenchmarkFEA: a 3D frame analysis engine and a benchmark environment
> for structural-modelling agents*, version 0.1.0, 2026.
> https://github.com/mohammedbala/BenchmarkFEA

If you use the verification problems under `verification/problems/ziemian/`, cite the data set
below as well.

## Acknowledgement

Six frames in `verification/problems/ziemian/`, and the published results used as their targets,
are taken from

> C. W. Ziemian, R. D. Ziemian, "Steel benchmark frames for structural analysis and validation
> studies: Finite element models and numerical simulation data", *Data in Brief* 39 (2021) 107564,
> [doi:10.1016/j.dib.2021.107564](https://doi.org/10.1016/j.dib.2021.107564); data set: Mendeley
> Data, [doi:10.17632/39sjhchwtx.1](https://doi.org/10.17632/39sjhchwtx.1),

licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
Changes were made: the data were converted to this project's input format, and the other changes
are listed in [`verification/problems/ziemian/README.md`](verification/problems/ziemian/README.md).
The authors of the data set are not connected with this project and have not reviewed this use.

All references that were opened, and what was taken from each, are listed in
[docs/REFERENCES.md](docs/REFERENCES.md).
