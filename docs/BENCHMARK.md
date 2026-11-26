# The benchmark: tasks, scoring, environment, dataset

An agent is handed a structural model in which something was changed in a way that matters (a
*curve ball*, [CURVEBALLS.md](CURVEBALLS.md)) and is asked to predict results, to say what is wrong,
or to repair it. Its answer is scored against a reference that passed the verification pipeline
([VERIFICATION.md](VERIFICATION.md)).

> **"Verified" means only: passed the tiers listed, at the tolerances stated.** It is a statement
> about agreement between computations. Nothing in this benchmark has been reviewed, approved or
> certified by anyone, and nothing here says so.

Code: `benchfea.bench.tasks` (cases and tasks), `benchfea.bench.scoring`, `benchfea.bench.env`
(Gymnasium environment `BenchFEA-v0`), `benchfea.bench.splits`, `benchfea.bench.dataset`.
Command line: `benchfea score`, `benchfea bench make-dataset`, `benchfea bench run-agent`.

## 1. Cases: only verified references are served

A **case** is one curve-balled model with its reference:

```python
from benchfea.bench.tasks import build_case, make_task
case = build_case("portal_frame", "support_fixed_to_pinned", seed=3, difficulty=2)
case["status"]        # "verified" | "unavailable"
case["skips"]         # every attempt that was dropped, with stage, category and reason
```

How a case is made (deterministic for `(family, curve ball, seed, difficulty)`):

1. the base model is `bench.generate(family, g, difficulty)`; its name is replaced by a neutral one
   (the generator's name states the seed) before anything is verified;
2. the curve ball is applied with the same seed `g` (`curveballs.apply`): the first target with a
   material, measured effect;
3. the curve-balled model is put through `benchfea.verify.verify_model` — the analysis results, or
   the mechanism diagnosis if the model is unstable — and, if that record is `verified` and the
   model has results, once more **with its design results** (`design=True`, tier F). The base model
   is verified in the same two ways;
4. the case is accepted if the reference it needs is `verified`:

| effect class of the curve ball | reference that must be verified |
|---|---|
| `changes-analysis` | analysis results of the curve-balled model |
| `creates-mechanism` | mechanism diagnosis of the curve-balled model (tiers N, K, C, E) |
| `design-only` | analysis results **and design results** of the curve-balled model |
| `invalid-model` | none exists: the model is refused by validation, the pipeline has nothing to verify; never served |

**Derivation chain.** The generator seed of attempt `k = 0, 1, 2` is a 31-bit number
(`tasks.derived_seed(family, curve_ball, seed, difficulty, k)`; at most `MAX_ATTEMPTS = 3`):

* its upper 30 bits are the first 30 bits of
  `SHA-256("benchfea.bench|task-chain|<TASK_VERSION>|<split>|family|seed|difficulty|k")`;
* its lowest bit is the split of the instance (section 7): 1 for held-out, 0 for train.

If the model cannot be generated, the curve ball does not apply, or the reference is not verified,
the attempt is **skipped and logged** and the next one is tried. The curve ball enters the chain
through the split only, so the instances of one split with the same `(family, seed, difficulty)`
share their base models, and the two splits never generate from the same seed. The instance seed is
never itself a generator seed. When no attempt
is left the case is `unavailable` and nothing is served for it. No tolerance is changed and no
quarantined reference is used; a skip entry reads
`{"attempt", "generatorSeed", "stage": "generate" | "curve-ball" | "verify", "category", "reason"}`.

What this means in numbers for the dataset of this repository is in section 8.

## 2. Task types

`make_task(case, task_type, mode)` returns two dicts:

* `task` — **public**: `id` (opaque), `taskType`, `mode`, `family`, `difficulty`, `split`,
  `statement` (text), `model` (the input file), `answerFormat`, and per type `quantities`,
  `criteria`, `intendedDesign`, `catalogue`; `verification` says that the reference is verified and
  what the word means — not which tiers (they would tell a mechanism from a stable model);
* `answers` — **hidden**: the case key (family, curve ball, seed, difficulty), the change record,
  the reference solution, the skip log and the verification record (status, tiers passed, per tier
  the summary, the tolerances and the largest differences found).

Modes:

* `tool-allowed` — the agent may run the analysis tool;
* `no-solver` — the run-analysis tool is disabled; the agent reasons from the file. The task is the
  same; a prediction is judged with wider tolerances (section 4).

### 2.1 predict

Give named results of the model **as given** (the curve ball is in it; the agent is not told).
The quantities are chosen by a seeded draw among results that are significant (at least 5 % of the
largest value of the same dimension in their combination) and that the curve ball moved by at
least 5 %, one per node / member: two nodal displacements, two reactions, one member end force, one
largest moment along a member (fewer if the model has fewer).

| kind | what | value in |
|---|---|---|
| `displacement` | `ux..rz` of a node in a combination, global axes | results `/combinations/<c>/displacements/<node>/<comp>` |
| `reaction` | `Fx..Mz` at a support, on the structure | `/combinations/<c>/reactions/<node>/<comp>` |
| `memberEndForce` | `Fx`, `My`, `Mz` at end `i` / `j`, local axes, on the member | `/combinations/<c>/members/<m>/endForces/<end>/<comp>` |
| `memberMaxMoment` | largest absolute `My` / `Mz` along a member (≥ 0) | max over `/combinations/<c>/members/<m>/stations/<comp>` |
| `bucklingFactor` | lowest positive buckling factor of a combination | `/buckling/<c>/factors/0` |
| `designRatio` | governing design ratio of a design member | design `/members/<m>/governing/ratio` |

`bucklingFactor` and `designRatio` are asked for only when the verified reference contains them.
A model with a buckling analysis cannot pass tier C, so `bucklingFactor` is never asked for with
the present pipeline. `designRatio` is added to a predict task when the design results of the
given model are verified (tier F) — the member whose governing ratio the curve ball moved most;
few generated models qualify (section 9; 5 predict tasks of the dataset).

Answer: `{"values": {"q1": 0.0123, "q2": -4.5, ...}}` in the units of the file.

Not available for a curve ball that creates a mechanism (there are no results).

### 2.2 diagnose

The statement gives the **intended design in words** (`tasks.describe` of the model before the
curve ball: every value, section by section — units, nodes, supports, materials, sections, members,
constraints, loads, combinations, analysis settings, design data) — not the base file. The agent
must name what was changed:

```json
{"category": "supports", "entities": {"nodes": ["N1"]}, "direction": "removed", "mechanism": null}
```

* `category` — one of the curve-ball categories: `constraints`, `offsets`, `releases`, `design`,
  `supports`, `members`, `units`, `loads`, `second_order`;
* `entities` — ids of what was changed, as named in the file (`{kind: [ids]}`, or a plain list);
* `direction` — derived mechanically from the edits of the change record
  (`tasks.change_direction`): `added` every edit adds something to the file; `removed` every edit
  removes something; `reversed` every edit replaces a number by its negative; `increased` /
  `decreased` every edit replaces a number by one of larger / smaller magnitude; `changed`
  anything else (a substituted name, a mixture);
* `mechanism` — `null`, or, **if the model as given is a mechanism**, its description
  `{"kind": "node" | "rigid-body" | "distributed", "motion": "translation" | "rotation" |
  "translation+rotation", "nodes": [ids of the nodes that move]}`. For a mechanism case this is the
  answer that is scored: it is matched to the verified diagnosis of the solver.

The statement is worded the same whether or not the model is a mechanism.

### 2.3 repair

The statement gives the intended design in words and the criteria; the agent submits a complete
model (in the environment: the working copy it edited).

| criterion | check |
|---|---|
| stable | the submitted model is accepted by the validator and analysed by the engine |
| verified | the reference of the **submitted** model passes `verify_model` (every applicable tier) |
| quantities | each named quantity is within 2 % of its value in the **verified intended design** (partial credit to 10 %) |
| design | every design check passes — **only** where the design results of the intended model are verified; otherwise the task states that this is not a criterion |

The named quantities are chosen like those of a predict task, from the results of the intended
design, **only** among those the curve ball moved by at least 15 % — beyond the zero-credit error,
so the model as given earns nothing on them and a repair can be told from no repair; if the curve
ball moved no significant result that far, no repair task is served on the case. For a mechanism
they are chosen among the significant results of the intended design (the model as given has
none). Their target values are not given. Values are compared physically: a repaired model written in other units is the same
structure.

A repair task needs the reference of the intended (base) model to be verified.

### 2.4 design

Choose the section of each member so that every AISC 360-22 member check passes at least weight.
Members that share a section in the given model form a group; a member may take any section of the
generic catalogue of the same shape (`bench.CATALOGUE`: I-shapes of the beam and column series,
boxes, pipes). Answer: `{"sections": {"B1": "I-450x200x14x9", ...}}` (members not named keep their
section).

* **Feasible**: analysed (strict), design status `pass` or `conditional`, and the reference of the
  submitted design passes `verify_model(design=True)`.
* **Score** = feasibility gate × min(1, reference weight / submitted weight); weight = Σ unit
  weight × area × flexible length.
* **Reference weight** = the *best found by this deterministic search* (`tasks.reference_design`):
  (1) every group at the heaviest section it may take — if that is not feasible there is no
  reference; (2) group by group, in file order, the lightest feasible section with the other groups
  unchanged; (3) repeat until a pass changes nothing (at most 4 passes). **It is not a proven
  optimum**; the search moves a group as a whole, so a submission that sizes members individually
  can be lighter (the score is capped at 1).

A design task is served only if the design results of the given model **and** of the reference
design are verified. With the present tier-F coverage few models qualify (section 9): the dataset
contains 2 design tasks (on continuous beams); the other code paths are tested with a test device
that counts every clause item as covered (`tests/bench_helpers.py`).

## 3. What each task type needs to be served

| task type | `changes-analysis` | `creates-mechanism` | `design-only` |
|---|---|---|---|
| predict | yes | no (no results) | design ratio only; needs verified design results |
| diagnose | yes | yes (mechanism matched to the verified diagnosis) | needs verified design results |
| repair | if the base reference is verified | if the base reference is verified | needs verified base design results |
| design | needs verified design results | no | needs verified design results |

`tasks.task_types_for(case)` gives the reason for each type that cannot be served.

## 4. Scoring

```python
from benchfea.bench.scoring import score_answer
out = score_answer(task, answers, answer)     # answer: dict or JSON text
out["score"], out["valid"], out["reasons"], out["summary"]
```

```
benchfea score answer.json --task data/tasks/train/<id>.task.json [--answers FILE] [--json]
```

The score is in [0, 1]; `reasons` is machine readable (per-quantity errors and credits, the criteria
that failed, the pipeline outcome). A malformed answer scores 0 with the reason and `valid: false`;
the scorer does not raise.

**Quantities** (predict, repair). Relative error
`e = |submitted − reference| / max(|reference|, floor)`, `floor` = 0.001 × the largest value of the
same dimension in the same combination. Credit: 1 for `e ≤ full`; 0 for `e ≥ zero`; in between
`1 − log(e / full) / log(zero / full)` (linear in log e; the geometric mean of the two bounds earns
one half). The score is the mean credit.

| task, mode | full credit up to | no credit from | quantity kinds |
|---|---|---|---|
| predict, tool-allowed | 1 % | 5 % | all |
| predict, no-solver | 5 % | 25 % | all |
| repair, either mode | 2 % | 10 % | all |

The full-credit tolerance is never tighter than twice the largest tolerance at which a reference
can have been verified (tier C admits a foreign solver up to an estimated error of 0.5 %).
A missing value, a text, NaN or infinity earns no credit for that quantity.

**Diagnosis of a change** (model as given is stable). Weights: entities 0.5, category 0.3,
direction 0.2.

* entities: F1 of *precision* = named ids that are entities of the change record (of any kind) /
  named ids, and *recall* = the best over the kinds of the record of (ids of that kind named) / (ids
  of that kind) — naming every member touched, or every node touched, is complete. The kind an
  answer files an id under is not compared;
* category: exact match; direction: exact match, credited only if an entity was found;
* a mechanism claimed for a stable model multiplies the score by 0.5.

**Diagnosis of a mechanism.** Against each mode of the verified diagnosis: nodes 0.6 (F1 of the
nodes named and the nodes of the mode), kind 0.25, motion 0.15 (kind and motion are credited only
if a node was found); the best mode counts. `mechanism: null` scores 0.

**Repair.** Gate (stable and verified, and the design criterion where it applies) × mean credit of
the named quantities. A submitted model that lacks a named node or member earns nothing on that
quantity.

**Design.** Section 2.4.

## 5. The environment

`benchfea.bench.env.BenchFEAEnv` is a `gymnasium.Env`, registered as `BenchFEA-v0` when the module
is imported. One episode is one task.

```python
import json, gymnasium
import benchfea.bench.env                               # registers "BenchFEA-v0"
env = gymnasium.make("BenchFEA-v0", split="train", max_steps=30, cache_dir="/tmp/benchfea-cache")
obs, info = env.reset(seed=7, options={"family": "portal_frame", "task_type": "predict"})
obs, reward, terminated, truncated, info = env.step(json.dumps({"tool": "run_analysis", "args": {}}))
```

Constructor: `split` (`"train"` | `"held-out"`), `max_steps` (step budget, default 30),
`cache_dir` (cache of generated and verified cases; default in memory), `max_draws`,
`score_budget`, `render_mode` (`None` | `"ansi"`), `families`, `curve_balls`, `difficulties`,
`seed_range` to restrict the random draw, and `dataset_dir`.

Two ways of serving tasks:

* **on the fly** (default) — the instance is generated, curve-balled and verified at `reset`
  (section 1). Unlimited instances; a reset costs from a second (a small beam) to minutes (a
  second-order 3D frame whose first attempts are quarantined). Use `cache_dir`, and restrict
  `families` / `difficulties` when the time matters;
* **from a stored dataset** — `BenchFEAEnv(dataset_dir="data")` draws uniformly among the stored
  tasks of the split that match the options; a reset costs no analysis. The stored task is the
  task the on-the-fly path builds for the same instance (tested).

**Spaces.**

* `action_space = Text(max_length=2 000 000, charset=printable ASCII)` — the JSON text of one tool
  call `{"tool": name, "args": {...}}`. `step` also accepts the same call as a dict.
* `observation_space = Dict(statement: Text, input_file: Text, tool_result: Text,
  steps_remaining: Discrete(max_steps + 1))` — `statement`: the task; `input_file`: the model as
  given (JSON text, the same in every observation of the episode); `tool_result`: JSON text of what
  the last tool returned; all ASCII.

**`reset(seed=None, options=None) -> (observation, info)`.** `options` may fix any of `family`,
`curve_ball`, `difficulty`, `task_type`, `mode`, `instance_seed`; an unknown key or value is a
`ValueError`. What is not fixed is drawn from `env.np_random` (seeded by `seed`): a served
(family, curve ball) pair of the split, a level of the family, an instance seed (moved up to the
next seed of the split), a task type that can be served on the case, a mode. A drawn instance
without a verified task is skipped (the reason goes to `env.skip_log`, harness side) and another is
drawn, at most `max_draws` times. If the options fix an instance that has no verified task, `reset`
raises `TaskUnavailable` with the log — it never serves an unverified reference.

**`step(action) -> (observation, reward, terminated, truncated, info)`.** Tools:

| tool | args | returns (`tool_result`) |
|---|---|---|
| `inspect` | `path`: JSON pointer (default: the whole model) | the value in the **working copy** |
| `edit` | `ops`: list of `{"op": "add" \| "remove" \| "replace", "path", "value"}` | applied together to the working copy, then validated; refused with the problems, changing nothing, if a path does not fit or the edited model is not valid. `add` sets an object member or inserts into an array (`-` appends) |
| `run_analysis` | `design`: bool | the engine on the working copy: a summary per combination (and per design member); a refusal (invalid, unstable with the mechanism diagnosis, beyond the critical load, …) comes back as `{"ok": false, "refusal", "problems", "diagnosis"}` — it never raises. **Disabled in mode `no-solver`** |
| `read_results` | `path`, `document`: `results` \| `design` | a value of the last analysis (results document, [CONVENTIONS.md](CONVENTIONS.md)) |
| `submit` | `answer` | ends the episode; the answer is scored. For a repair task without `answer.model` the working copy is submitted |

* **Reward**: 0 for every step except `submit`, whose reward is the score in [0, 1].
* **Termination**: `terminated` after `submit`. **Truncation**: `truncated` when `max_steps` tool
  calls were used without a `submit` (reward 0). Every call counts, a malformed or unknown one
  included (it returns an error observation). A value longer than 200 000 characters is not
  returned whole: the keys are, with a note to ask for a part.
* **`info`** before the end: `task_id`, `task_type`, `mode`, `family`, `difficulty`, `split`,
  `steps_used`, `max_steps`, `tools`, `answer_format`, `quantities` / `criteria` / `catalogue`,
  `verification` (`status: verified` and the meaning of the word), `instances_skipped`.
  **Nothing of the reference** — curve ball, change record, reference values, verification tiers,
  seeds — is in an observation or in `info` before the episode ends (tested). At the end
  `info["score"]` is the scoring record and `info["reference"]` the hidden part of the task.
* **Determinism**: the same `seed` and `options` give the same task and, for the same actions, the
  same observations, rewards and infos (tested with two environments and a scripted sequence).
* `env.unwrapped.hidden_answers()` is a harness-side accessor of the hidden part (used by the
  scripted oracle agent); it is not part of the agent's interface.

**`gymnasium.utils.env_checker.check_env`** passes on the environment, the render check included
(`tests/test_env.py`, with the random draw restricted to four small cases so that the test is
fast): spaces, reset / step signatures and return types, seeding determinism, step determinism,
fresh observation and info objects, render and close. Its only warning is that no render fps is
declared (the render mode is `"ansi"`: the statement and the last tool result as text). What
`check_env` cannot cover: it samples random text from the action space, which exercises only the
malformed-action path — the tools, the reward and the budget are tested separately.

## 6. Running an agent

An agent is a callable `agent(observation, info) -> action`.

```python
from benchfea.bench.env import run_episode
episode = run_episode(env, agent, seed=1, options={"task_type": "diagnose"})
episode["reward"], episode["steps"], episode["transcript"], episode["info"]["score"]
```

```
benchfea bench run-agent examples/agents/oracle_agent.py --dataset data --episodes 5 --seed 0
benchfea bench run-agent examples/agents/random_agent.py --dataset data --task-type predict
benchfea bench run-agent examples/agents/oracle_agent.py --family continuous_beam --episodes 3   # on the fly
benchfea bench run-agent my_agent.py:make_agent --split held-out --episodes 20 --json
```

The file must define `make_agent(env) -> agent`. In `examples/agents/`:

* `random_agent.py` — baseline: submits a no-op answer (zeros, a random category and member, the
  model as given). The floor.
* `oracle_agent.py` — scripted: reads the hidden reference through `hidden_answers()` and submits
  it. It is not solving anything; it shows that every served task is solvable and that the perfect
  answer scores 1.0 through the normal scoring path (repairs are re-solved through the pipeline).
* `llm_agent_template.py` — how to plug in a language model: `make_agent(env, llm=callable)`, where
  `llm(prompt: str) -> str` returns one tool call as JSON. No API key and no network call is in the
  repository; without `llm` a stub submits an empty answer.

Scoring a stored task without the environment: write the answer to a file and use `benchfea score`.

## 7. Train / held-out split

Defined by the committed file `src/benchfea/bench/splits.json` and one rule
(`benchfea.bench.splits.is_held_out(family, curve_ball, seed)`):

* **held-out combinations** — 16 (family, curve ball) pairs, two per family, that
  never appear in the train split whatever the seed;
* **held-out seeds** — for every other pair, the seeds with `seed % 5 == 4`.

`seed` is the seed of the instance as requested; the derivation chain of section 1 does not move an
instance to the other split. `is_held_out` is the only definition of the split; everything else
(the dataset directories, the environment, the generator seeds) calls it.

**How the held-out combinations are selected.** By a rule that looks at names only
(`splits.select_held_out`; the test suite checks that `splits.json` equals its result):

1. *Candidates of a family*: the curve balls whose cell for that family in the applicability table
   of [CURVEBALLS.md](CURVEBALLS.md) section 5 is not `-`, except cells that are `I`
   (`invalid-model`) alone — such a pair is never served, by the definition of section 1. The
   candidates are stored in `splits.json` (`heldOutSelection.applicable`).
2. *Order*: ascending SHA-256 (hexadecimal) of the UTF-8 text
   `benchfea.bench|held-out-pair|<family>|<curve ball>`.
3. *Selection*: the first two candidates of every family.

Whether a pair has verified instances plays no part, so a held-out combination can turn out to have
none; it is then simply not served (it is not replaced by another pair). When the applicability
table changes, the rule is applied again; nothing else changes the list.
An earlier list was partly chosen after it was known which pairs had verified instances (three
pairs had been replaced by served ones); the rule above replaces that list.

| family | held-out curve balls |
|---|---|
| `continuous_beam` | `release_mechanism`, `release_added` |
| `portal_frame` | `sections_swapped`, `release_pinned_frame` |
| `multistorey_2d` | `design_unbraced_length_changed`, `load_moved` |
| `multistorey_3d` | `section_heavier`, `support_fixed_to_pinned` |
| `braced_frame` | `load_moved`, `design_k_changed` |
| `pipe_rack` | `member_rolled_90`, `load_moved` |
| `stair_frame` | `support_roller_direction`, `release_removed` |
| `diaphragm_frame` | `design_lb_full_where_braced`, `load_direction_changed` |

`servedCombinations` in the file lists the pairs for which a dataset build produced a verified case;
the environment draws from it. It is a measurement and has to be refreshed by every rebuild of the
dataset (section 8).

**Policy.** `BenchFEAEnv(split="train")` — the default — refuses a held-out instance: `reset`
raises `ValueError` when the options name one, and the random draw never produces one. Held-out
instances are served only by an environment constructed with `split="held-out"`, which in turn
refuses train instances. Train on `train`; report on `held-out`; do not tune on `held-out`.
The split is by construction, not by secrecy: the generators are in the repository, so anyone can
generate held-out instances — the policy relies on using the environment as documented.

**Base models are not shared between the splits.** The generator seed of every attempt carries the
split in its lowest bit (section 1), so a base model of the train split and one of the held-out
split are never generated from the same `(family, seed, difficulty)` input — also not for a held-out
*combination* asked for with a train seed, whose base models used to be those of the train tasks of
the same family, seed and difficulty. What is checked, and how:

* *seeds* — for every planned cell of the dataset and every attempt of the chain, the sets of
  `(family, difficulty, generator seed)` of the two splits are disjoint (test suite, no analysis);
* *documents* — two different seeds could still give the same document. The test suite generates
  the base models of the planned cells of the cheap families (all attempts) and compares SHA-256
  hashes between the splits; `BENCHFEA_SPLIT_FULL=1` extends this to every planned cell (about an
  hour). On the dataset as built, `tools/audit_dataset.py` compares the hashes of the base
  models and of the curve-balled models between the splits; a shared one is a finding (section 8).

The split separates instances, not knowledge: families, curve balls and the generator are the same
code for both splits.

## 8. Dataset

A dataset of your own, in a directory of your choice (do not write into `data/`, which holds the
committed dataset; the commands that built it are listed further down in this section):

```
benchfea bench make-dataset --out my-dataset --seeds 0,1,4 --workers 4 --cache /tmp/benchfea-cache
```

The layout under the output directory, shown for the committed dataset in `data/`:

```
data/tasks/manifest.json                       counts, every task, every cell, the quarantine log, versions
data/tasks/{train,heldout}/<id>.task.json      public: statement, answer format, quantities / criteria
data/tasks/{train,heldout}/<id>.model.json     public: the input file of the task (one line of JSON)
data/answers/{train,heldout}/<id>.answers.json hidden: change record, reference solution, verification record
data/CHECKS.json                               summaries of the two checks run after a build (not written by the build)
data/QUARANTINE_NOTES.md                       written diagnoses of the attempts logged as disagreement
```

A **cell** is `(family, curve ball, seed, difficulty)`; `dataset.plan(seeds)` takes every family ×
every curve ball × the seeds, with the difficulty level rotating over the levels of the family
(`dataset.level_for`). Every task type that can be served on a verified cell becomes one task, in
the mode `dataset.mode_for` gives (a fixed pseudo-random choice per cell and task type). The build
is deterministic, parallel (`--workers`) and resumable (`--cache`: generated models, verification
records and foreign solves are cached; a second run re-reads them).

The reference solution in an answers file is what the task is scored on (the named values, the
change record with its entities and direction, the mechanism diagnosis, the intended model of a
repair) with the SHA-256 of the full results document; the full document is reproduced by
`benchfea run` on the model file. The verification record is stored without the per-check details.

**The dataset in this repository** ([data/DATA_CARD.md](../data/DATA_CARD.md)). Every number
below stands between `generated` markers and is written by `python3 tools/build_dataset_docs.py`
from `data/tasks/manifest.json` and `data/CHECKS.json`; the test suite runs its `--check`.

<!-- generated:headline -->
**748 verified tasks** from 263 verified cells of 528 (every family × every curve ball × seeds 0 and 4): 236 predict, 263 diagnose, 249 repair, 0 design; 348 train, 400 held-out.
<!-- /generated:headline -->

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

Cells by family (a cell is verified when one attempt of its chain gave a verified reference):

<!-- generated:families -->
| family | cells | verified | tasks | cells not verified: stage:category of the last attempt |
|---|---|---|---|---|
| `braced_frame` | 66 | 25 | 68 | curve-ball:not-applicable 32, verify:not-covered 7, verify:not-representable 2 |
| `continuous_beam` | 66 | 28 | 81 | curve-ball:not-applicable 35, verify:not-covered 3 |
| `diaphragm_frame` | 66 | 41 | 116 | curve-ball:not-applicable 15, verify:not-covered 8, verify:refused 2 |
| `multistorey_2d` | 66 | 34 | 99 | curve-ball:not-applicable 26, verify:not-covered 6 |
| `multistorey_3d` | 66 | 36 | 102 | curve-ball:not-applicable 21, verify:not-covered 8, verify:not-representable 1 |
| `pipe_rack` | 66 | 34 | 95 | curve-ball:not-applicable 25, verify:not-covered 5, verify:not-representable 2 |
| `portal_frame` | 66 | 33 | 94 | curve-ball:not-applicable 25, verify:not-covered 7, verify:not-representable 1 |
| `stair_frame` | 66 | 32 | 93 | curve-ball:not-applicable 29, verify:not-covered 5 |
<!-- /generated:families -->

Cells by curve ball:

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

Quarantine log (`quarantineLog` in the manifest):

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

Checks run after the build. The first puts every task through the normal scoring path; the second
uses means that are not the pipeline that produced the references (the split re-evaluated from
`splits.json`, base and curve-balled model hashes per split, every stable model solved again by the
engine and by the independent tableau solver of `tests/independent`, the number of mechanisms of
every unstable model from the null space of the tableau, right and wrong answers through the
scorer, a probe of the first observation for anything of the change record). The test suite repeats
the scoring for every predict, diagnose and design task and for a sample of the repair tasks, and
rebuilds one task from its key.

<!-- generated:checks -->
* `tools/check_dataset_scores.py` (normal scoring path, repairs and designs re-solved through the pipeline): the answer built from the hidden reference scores 1.0 on 748 of 748 tasks; the model as given, handed in as the repair, scores below the reference answer on 249 of 249 repair tasks (0 on 249; largest score 0.0). Findings: 0.
* `tools/audit_dataset.py` (not the pipeline that produced the references): 748 tasks on 86 base models; curve-balled models in both splits: 0; held-out tasks on a base model that train tasks also use: 0. 236 stable models solved again by the engine (1392 reference values recomputed, largest relative difference 0.0e+00) and by the independent tableau solver of `tests/independent`: 888 first-order combinations, largest relative difference 7.3e-12; 331 second-order combinations with the tableau refined (it takes the axial force constant over a piece of a member; 4, 16, 64, 256 pieces per member, raised while the difference to the engine is above 1e-3): 7 needed more than 16 pieces, the largest difference at the last count run is 1.2e-04, and between the last two counts run the difference falls by a factor of at least 6.7 (a finding below 4). Number of mechanisms of 27 unstable models from the null space of the tableau: 27 agree. Right and wrong answers through the scorer: 499 predict and diagnose tasks. First observation of 25 tasks probed for parts of the change record: 0 found. Findings: 0.

These are computations by this project's code. Nobody has reviewed, approved or certified a task, a model or a reference.
<!-- /generated:checks -->

Gaps of this build (curve balls with no or few tasks and why, design results, second order,
held-out combinations without a verified instance):

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

## 9. Limits

* **What "verified" covers.** Agreement of computations at the stated tolerances (analysis results:
  tiers C, D, E and, for straight member lines, A; mechanism diagnoses: tiers N, K, C, E). The
  engine, the ODE oracle and the translators to the two foreign solvers come from one project and
  share the fixed-end-force code; an error common to all routes would not be detected. Nobody has
  reviewed, approved or certified any task or reference.
* **Design results are rarely a verified reference.** Tier F requires, for every limit state a
  model uses, a hand calculation and a reproduced strength of the AISC Design Examples in the test
  suite (`src/benchfea/verify/design_coverage.json` lists the items with and without it). A model
  that uses one item without it has no verified design results: no design-only curve ball, no
  `designRatio` quantity, no design criterion in a repair task and no design task on it. How many
  there are in the dataset is in section 8 ("Design results"). The code paths are tested with a
  test device that counts every item as covered (`tests/bench_helpers.py`), which says nothing
  about any design value.
* **Second-order models are served only where two foreign solvers give an eligible error bound.**
  Tier C compares the engine (exact stability functions) with OpenSeesPy and PyNite, allowing each
  foreign result the error bound computed from three of its own refinement levels
  ([CROSSCHECK.md](CROSSCHECK.md) 1.3). Second order together with rigid offsets, rigid links or
  diaphragms, and crossing truss braces under self weight in second order, leave fewer than two
  eligible solvers; those models are quarantined as `not-representable`, not served. The share of
  second-order tasks is in section 8.
* **Every quarantine is left as it is**, with its reason in the log; the counts by stage, category
  and tier are in section 8, written diagnoses of the attempts logged as `disagreement` in
  `data/QUARANTINE_NOTES.md`. Nothing was loosened to admit a case.
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
