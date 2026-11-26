# Documents

One line for each document of the repository. Start with the first three.

## Start here

| Document | What it is |
|---|---|
| [../README.md](../README.md) | What the software is, how to install it, and one example each of running, changing, scoring and verifying a model |
| [../LIMITS.md](../LIMITS.md) | What the software does not do, what cannot be verified, and what it must not be used for |
| [CONVENTIONS.md](CONVENTIONS.md) | Axes, signs and units: the single place where they are defined, with figures |

## Using the engine

| Document | What it is |
|---|---|
| [INPUT_REFERENCE.md](INPUT_REFERENCE.md) | Every field of an input file: path, value, unit, default, description, example (generated from the schema) |
| [OUTPUT_REFERENCE.md](OUTPUT_REFERENCE.md) | Every key of the results file and of the design file the engine writes |
| [../examples/README.md](../examples/README.md) | Five worked input files, from a single beam to a braced frame, each with an explanation and its key results |
| [FORMAT_VERSIONING.md](FORMAT_VERSIONING.md) | What the `formatVersion` of a file means and which versions the engine accepts |
| [../data/user/README.md](../data/user/README.md) | Where to put your own copy of the AISC Shapes Database, and the file layout the reader expects |

## What the engine computes

| Document | What it is |
|---|---|
| [THEORY.md](THEORY.md) | The mechanics of the solver, derived from equilibrium and virtual work, with the test that checks each result |
| [ARCHITECTURE.md](ARCHITECTURE.md) | For developers: the modules, the data that passes between them, and the contracts of the solver and the design module |

## Verification

| Document | What it is |
|---|---|
| [VERIFICATION.md](VERIFICATION.md) | What "verified" means, the tiers A to F, the one tolerance policy, the quarantine rule, the format of a verification problem |
| [VERIFICATION_REPORT.md](VERIFICATION_REPORT.md) | The outcome of an actual run on every problem of the library and on the model bank, with every difference and tolerance (generated) |
| [CROSSCHECK.md](CROSSCHECK.md) | How models are translated to OpenSeesPy and PyNite, the error bound of their results, and the independent member-equation solver |
| [REFERENCES.md](REFERENCES.md) | Every reference the project names: whether it was opened, which pages were read, and what in the repository rests on it |
| [REQUIREMENTS_TRACE.md](REQUIREMENTS_TRACE.md) | Each requirement of the project specification, where it is met, and whether it is met, partly met or not met, with the reason |
| [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) | Every place where a reference or AISC 360-22 was silent or ambiguous, what the software does there now, and what an answer would change |
| [../verification/problems/aisc-commentary/README.md](../verification/problems/aisc-commentary/README.md) | The benchmark problems of the AISC 360-22 Commentary, Section C2.1 |
| [../verification/problems/aisc-de/README.md](../verification/problems/aisc-de/README.md) | Problems from the AISC Design Examples (all quarantined, with the diagnoses) |
| [../verification/problems/mgz/README.md](../verification/problems/mgz/README.md) | Examples from McGuire, Gallagher and Ziemian, *Matrix Structural Analysis*, and the quarantine diagnoses |
| [../verification/problems/ziemian/README.md](../verification/problems/ziemian/README.md) | Frames from the Ziemian and Ziemian data set (CC BY 4.0), with the attribution and the changes made |
| [../verification/problems/staad/README.md](../verification/problems/staad/README.md) | Problem statements from the public STAAD.Pro verification pages, with targets recomputed here |
| [../verification/problems/ansys/README.md](../verification/problems/ansys/README.md) | Problem statements from the public ANSYS verification manual, with targets recomputed here |

## The benchmark

| Document | What it is |
|---|---|
| [BENCHMARK.md](BENCHMARK.md) | Task types, answer formats, scoring, the Gymnasium environment, the train / held-out split, the dataset |
| [CURVEBALLS.md](CURVEBALLS.md) | Every curve ball: what it changes, why it matters, where it applies, with measured examples (generated) |
| [DIFFICULTY.md](DIFFICULTY.md) | The model generators, their families and the measured difficulty levels |
| [../data/DATA_CARD.md](../data/DATA_CARD.md) | The dataset in `data/`: how it was built, what it contains, known gaps, intended use, licence |
| [../data/QUARANTINE_NOTES.md](../data/QUARANTINE_NOTES.md) | Written diagnoses of the dataset build's attempts that were quarantined for a disagreement |
| [../examples/agents/README.md](../examples/agents/README.md) | Example agents for the environment |

## Contributing

| Document | What it is |
|---|---|
| [../CONTRIBUTING.md](../CONTRIBUTING.md) | Setting up, running tests and tiers, adding a model, a generator family, a curve ball, a verification problem or a design check; the reference, tolerance and wording rules |
| [../LICENSE](../LICENSE) | The MIT licence of the code and of the generated dataset |

Figures used by CONVENTIONS.md are in `figures/`; they are drawn by `figures/make_figures.py`.
