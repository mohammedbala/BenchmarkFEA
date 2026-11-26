# Contributing

The purpose of this repository is a set of reference solutions that can be trusted to the extent
stated in their verification records, and no further. A contribution that makes the engine faster
or larger while weakening that is not wanted. Please read the three rule sections at the end
(references, tolerances, wording) before you start; they are not negotiable.

## Set up

Python 3.11 or newer.

```bash
git clone https://github.com/mohammedbala/BenchmarkFEA.git
cd BenchmarkFEA
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,env,crosscheck,shapes]"
```

`dev` is pytest; `env` is Gymnasium; `crosscheck` is OpenSeesPy and PyNite, without which tier C
cannot run and the tests that need it fail; `shapes` is openpyxl, without which one test of the
`.xlsx` reader is skipped. Set one thread for the linear algebra library before running
anything long: the matrices are small, and a multi-threaded library was measured to be far slower.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
```

## Run the tests

The fast tests, while you work (about 2 minutes):

```bash
python3 -m pytest -q -m "not slow"
```

The full suite, before you propose a change (about 12 minutes; `-x` stops at the first failure):

```bash
python3 -m pytest -q -x -p no:cacheprovider
```

Nothing is skipped or weakened by the `slow` marker; it only names the slowest tests
(`tests/conftest.py`).

The dataset in `data/` was built with Python 3.13, numpy 2.5 and scipy 1.18. One test requires its
stored reference values to equal a fresh solve bit for bit; with the older numpy and scipy that are
installed on Python 3.11 it fails by a difference of about 4·10⁻¹² (measured). This is recorded as
an open question ([docs/OPEN_QUESTIONS.md](docs/OPEN_QUESTIONS.md), Q120); do not loosen the test.

## Run the verification tiers

One model, every tier that applies (add `--design` to include tier F). The exit code is 0 for
`verified` and 5 for `quarantined`; the output names the tier that was not passed and why:

```bash
benchfea verify examples/02_continuous_beam.json
```

The whole problem library. Exit code 0 means every *active* problem passes; quarantined problems
are evaluated and reported, never counted as passed:

```bash
benchfea verify --all
```

What the tiers are and what they prove: [docs/VERIFICATION.md](docs/VERIFICATION.md).

## Generated files

Several files are written by tools and must never be edited by hand. Each tool has a `--check`
mode that fails when the file on disk is not what the tool would write; continuous integration
runs them.

| File | Written by |
|---|---|
| `src/benchfea/schema/benchfea-input-1.0.schema.json` | `python3 -m benchfea.schema.build_schema` |
| `docs/INPUT_REFERENCE.md` | `python3 tools/build_input_reference.py` |
| `src/benchfea/verify/design_coverage.json` | `python3 tools/build_design_coverage.py` |
| `examples/expected.json`, the result tables in `examples/*.md`, `examples/01_simple_beam.yaml` | `python3 tools/build_examples.py` |
| `docs/CURVEBALLS.md` | `python3 tools/build_curveball_docs.py` (several minutes) |
| `src/benchfea/bench/splits.json` (the held-out combinations) | `python3 tools/apply_split_rule.py` |
| the marked blocks of `data/DATA_CARD.md` and `docs/BENCHMARK.md` | `python3 tools/build_dataset_docs.py` |
| `docs/VERIFICATION_REPORT.md` | `python3 tools/build_verification_report.py` (no `--check`: it is a record of a run) |
| `data/tasks/`, `data/answers/` | `benchfea bench make-dataset` (see [data/DATA_CARD.md](data/DATA_CARD.md), "How to rebuild") |

Links between the documents are tested by `python3 tools/check_links.py`.
`docs/OUTPUT_REFERENCE.md` and `docs/REQUIREMENTS_TRACE.md` are written by hand. A new key in an output
document must be added to the first (`tests/test_output_reference.py` fails otherwise); a change of
scope must be entered in the second.
`python3 tools/compare_dataset_sample.py` rebuilds a few cells of the dataset from nothing and
compares the files with the committed ones byte for byte (about 2 minutes for the default sample).

## How to add …

### (a) An example model

1. Write `examples/NN_name.json`. Use sections given by dimensions (`I`, `box`, `pipe`,
   `rectangle`), so that the file needs no shapes database. Every field is described in
   [docs/INPUT_REFERENCE.md](docs/INPUT_REFERENCE.md).
2. `benchfea validate examples/NN_name.json` must report no error; then run it and look at the
   results with the eyes of an engineer (reactions against the loads, deflected shape, signs).
3. Write `examples/NN_name.md` next to it: what the model is, what is worth noticing in the file,
   the commands, and the two headings `## Key results` and `## Design summary` (the tables under
   them are generated).
4. Name the quantities of the "Key results" table in `QUANTITIES` of `tests/test_examples.py`, and
   add every closed-form check you can derive (statics, a beam formula) to that test.
5. `python3 tools/build_examples.py`, add a row to the table in `examples/README.md`, run the
   tests.

The numbers of an example are what the engine computes. Do not present them as reference values of
a publication.

### (b) A generator family

A family is a module in `src/benchfea/bench/families/` that defines `NAME`, `LEVELS`, `FEATURES`
and `build(draw, level, unit_system, vertical, features)`; see the docstring of
`families/__init__.py` and an existing family.

1. Draw every random number from the `draw` object, so that the output is a function of
   (family, seed, difficulty) only.
2. Stay inside what both outside solvers can represent, or the models cannot be verified: no shear
   deformation, no rigid offsets together with second order, and the other conditions listed in
   `families/__init__.py`. A feature that leaves this envelope must be an opt-in flag.
3. Add the module name to the tuple in `families/__init__.py`.
4. Add the family to `tests/test_generators.py` (determinism, validity, the pinned output of one
   seed) and describe it in [docs/DIFFICULTY.md](docs/DIFFICULTY.md) with its measured difficulty
   properties.
5. Rebuild the curve-ball document and the split (`tools/build_curveball_docs.py`, then
   `tools/apply_split_rule.py`): the applicability table gains a row.
6. A change that alters the output of an existing (family, seed, difficulty) requires
   `GENERATOR_VERSION` in `src/benchfea/bench/generate.py` to be raised, and the dataset to be
   rebuilt.

### (c) A curve ball

A curve ball is a named, seeded, deterministic change to a model whose effect is *measured*, not
assumed ([docs/CURVEBALLS.md](docs/CURVEBALLS.md)).

1. **Registry.** Add a function decorated with `@register(...)` in the module of its category
   under `src/benchfea/bench/curveballs/` (see `supports.py` for a short one). The registration
   states `name`, `title`, `changes` (what it changes in the file), `why` (why it matters
   physically, in an engineer's words), `moves` (which results it moves), `effect_classes` and
   `precondition`. The function returns the candidate targets in the base model; it does not
   choose among them.
2. **Record.** You do not write the change record. It is the structural difference of the two
   documents, computed by `record.diff`, and it must apply and revert exactly; the tests check both
   directions for every curve ball.
3. **Measured materiality.** The registry analyses and designs the base and the changed model and
   accepts a candidate only when the measured effect reaches the thresholds of
   `measure.THRESHOLDS` and falls in a declared effect class. Do not add a curve ball that needs a
   lower threshold; a change that makes no difference is not a curve ball.
4. **Tests.** Add the curve ball to `tests/test_curveballs.py`: applicable and not applicable
   cases, determinism for a seed, and the direction of the effect you claim in `moves`.
5. **Documents regenerated.** `python3 tools/build_curveball_docs.py`, then
   `python3 tools/apply_split_rule.py`. If the second tool reports that the held-out combinations
   changed, the dataset has to be rebuilt.
6. A change to an existing curve ball that alters the output of any (model, name, seed) requires
   `CURVEBALL_VERSION` in `curveballs/registry.py` to be raised.

### (d) A verification problem

A problem is a directory `verification/problems/<source>/<problem-id>/` with two files
([docs/VERIFICATION.md, section 3](docs/VERIFICATION.md#3-problem-library)):

* `model.json` — a normal input file;
* `problem.json` — `id`, `title`, `tier` (`A` derived here, `B` published), `source` (`work`,
  `edition`, `location` = example, figure or equation number, `pages`, `notes`), `basis`,
  `analysis`, `targets`, `status` (`active` or `quarantined`), and for a quarantined problem
  `quarantine.diagnosis`.

**Basis values.** `published`: every target is a number printed in the source. 
`recomputed-from-published-statement`: only the problem statement is published; the targets were
computed in this repository by means that share no code with the engine. Such a problem is not a
reproduction of a published number, and the report counts it separately. `derived-here`: tier A;
the derivation script is kept under `verification/derive/`. A printed number the model was
calibrated on is marked "reproduced by construction, not evidence" in its `precisionNote`.

**Targets.** Each target is a path into the results document or a named extractor, a `value` in the
units of the model file and in the engine's sign convention (`sign` records the mapping from the
source's convention), and a `tolerance`.

**The single tolerance policy.** For a published number the tolerance is half a unit in the last
digit the source prints (`benchfea.verify.printed_precision_tolerance`), or the acceptance
tolerance the source itself states; `precisionNote` says which, and gives the number as printed.
For targets derived here it is 10⁻⁶. You may not choose a tolerance per case.

**The quarantine register.** A problem that does not reproduce, or whose model fails another tier,
gets `"status": "quarantined"` and a written diagnosis: what differs, by how much, the suspected
cause and the evidence for it. Sources with a builder script keep these diagnoses in one file
(for example `verification/problems/mgz/quarantine.py`). A quarantined problem stays in the
library, is still evaluated and reported, and is never counted as passed.

> **Never tune a tolerance. Never drop a case. Quarantine with a diagnosis.**

Do not change a model, a target or a tolerance to move a problem out of quarantine. If the
diagnosis is "transcription error in the problem file", fix the transcription and say so; if it is
anything else, the problem stays quarantined until the cause is removed.

Then: record the source in [docs/REFERENCES.md](docs/REFERENCES.md) (what was opened, which
pages), run `benchfea verify verification/problems/<source>/<problem-id>`, and rebuild the report
with `python3 tools/build_verification_report.py`.

### (e) A design check

1. **Implemented only from the text of the specification you have open.** Not from memory, not
   from another program, not from a summary or a textbook. If the text is silent or ambiguous, do
   not guess: return `notChecked` with the reason and record the question in
   [docs/OPEN_QUESTIONS.md](docs/OPEN_QUESTIONS.md).
2. **The section is cited at every equation**, in the code and in the `clause` of the check
   (section, equation number, printed page), so that a reader with the specification can follow
   each line. An engine rule that is not specification text is labelled as an engine rule.
3. **A hand-worked test**: `test_hand_*` in `tests/test_design_*.py`, with the arithmetic of the
   specification written out in the test body, reading the strength or the utilisation of the
   check.
4. **A Design Examples test**: a printed value of the AISC Design Examples reproduced through
   `design_examples.assert_printed` (one tolerance for all such tests), with the example number
   and page in the data file of the test. The printed value is cited; no text or table is copied.
5. `python3 tools/build_design_coverage.py`. **Without both tests the check stays outside tier F
   coverage**: it may be computed and reported, but no design result that uses it can be served as
   verified. That is the intended outcome, not a problem to work around.

## The reference rules

* **Nothing from memory.** No reference value, formula coefficient or table entry is typed from
  recollection.
* **Nothing from a source that was not opened.** If you did not read the page, you may not cite it
  and may not take anything from it. A vendor page that cites a textbook is evidence of the vendor
  page only.
* **Record every source in [docs/REFERENCES.md](docs/REFERENCES.md)**: what it is, how it was
  obtained, which pages were read, what uses it.
* **No copied text, tables or figures.** Cite, then compute. Problem data are restated in your own
  words; printed result numbers are quoted as numbers with their page.
* **No AISC documents and no shapes data in the repository.** Not the Specification, not the
  Manual, not the Design Examples, not the Shapes Database, not an extract of any of them. A user
  supplies their own copy ([`data/user/README.md`](data/user/README.md)).
* **No client or project data**, and no personal names, e-mail addresses or private file paths in
  any file.
* **The engine is written from the mechanics** ([docs/THEORY.md](docs/THEORY.md)). Do not port
  code from another analysis program.

## The wording rule

Never write that a result is checked, approved or certified by anyone — in code, output,
documents, commit messages or issue comments. "Verified" has one meaning here: *passed the tiers
listed, at the tolerances stated.* Where a check could not be made, say "not checked" and why.
