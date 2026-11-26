# Architecture (developer reference)

Concise map of the package and the Python APIs that modules share. Sign and
axis conventions are in [CONVENTIONS.md](CONVENTIONS.md); the file format rule
in [FORMAT_VERSIONING.md](FORMAT_VERSIONING.md).

## Pipeline

```
file (.json/.yaml) or dict
  └─ benchfea.io.load ──────────────► Document (data + line/col of every value)
       └─ benchfea.validate.validate ► list[Problem]  (schema layer + semantic layer, ALL problems)
            └─ benchfea.model.build_model ► Model  (internal units, defaults resolved, provenance)
                 ├─ benchfea.solve.analyze(model) ───────────► Results          [first-order core implemented]
                 │       └─ Results.to_document(model) ► "benchfea-results" 1.0.0, file units
                 └─ benchfea.design.check(model, results) ───► DesignResults    [D, E, F, G, H1, deflection limits]
                         └─ DesignResults.to_document(model) ► "benchfea-design" 1.0.0, file units
```

## Modules

| module | role |
|---|---|
| `benchfea.problems` | `Problem(path, message, hint, severity, code, line, col)`, `ModelInvalid`, `ModelUnstable`, `pointer()` |
| `benchfea.io` | `load(path \| dict) -> Document`, `loads(text, format)`, `dump(data, path)`, `dumps(data, format)`; position-tracking strict JSON parser and YAML loader |
| `benchfea.schema` | `benchfea-input-1.0.schema.json` (draft 2020-12), `load_schema()`, `validator()`; `build_schema.py` generates the JSON from `benchfea.defaults` |
| `benchfea.defaults` | every built-in default and its reason (single source for schema, builder and report) |
| `benchfea.units` | `UnitSystem(force, length)`: `factor(dim)` (exact `Fraction`), `to_internal`, `from_internal`, `label`; `INTERNAL` = N, mm |
| `benchfea.unitmap` | `unit_map(data)` (JSON pointer → dimension, from the schema's `x-unit`), `convert_document(data, target_units)` |
| `benchfea.validate` | `validate(source) -> list[Problem]`, `check(source)` (raises `ModelInvalid`), `validate_schema`, `validate_semantics` |
| `benchfea.geometry` | `local_triad(...)`, `member_geometry(...)`, `is_vertical`, `AXIS` |
| `benchfea.sections` | `SectionProps`, `i_shape`, `box`, `pipe`, `rectangle`, `rectangle_torsion_constant`, `from_definition`, `check_dimensions` |
| `benchfea.shapes` | reader of the user-supplied AISC Shapes Database (`lookup(name, length_unit) -> SectionProps`), `ShapesDatabaseError` with hint |
| `benchfea.designparams` | resolution of design variables (member > group > global > default) with provenance |
| `benchfea.model` | dataclasses of the resolved model and `build_model(source) -> Model` |
| `benchfea.results` | dataclasses `Results`, `CombinationResult`, `MemberResult`, `BucklingResult`, `Envelope`, `InvariantsReport`, `DesignResults`, `DesignMemberResult`, `DesignCheck` and the output-document writers |
| `benchfea.cli` | `benchfea validate \| run \| design \| perturb \| verify \| score` (all implemented; also `generate` and `bench make-dataset | run-agent`) |
| `benchfea.design` | steel member checks: `check(model, results)`, `check_member`, registry of editions, `aisc360_22` chapter modules — see "Contract for `benchfea.design`" |
| `benchfea.solve` | solver: `analyze(model) -> Results`; first-order core ("Solver core" below) and the nonlinear analyses — second order, linear buckling, tension/compression-only members, stiffness reduction, notional loads ("Nonlinear analyses" below); theory in [THEORY.md](THEORY.md) |
| `benchfea.oracle` | independent single-member / continuous-beam reference by direct integration of the member ODEs (`scipy.integrate.solve_bvp`): `Beam`, `Bar` — see [CROSSCHECK.md](CROSSCHECK.md) |
| `benchfea.crosscheck` | translators to OpenSeesPy and PyNite (optional packages, run in a subprocess), `run`, `compare`, `cross_check` (tier-C decision), fidelity records, `CAPABILITIES` — see [CROSSCHECK.md](CROSSCHECK.md) |
| `benchfea.verify` | verification pipeline: `verify_model(source, tiers=, design=, budget=) -> VerificationRecord` (`status` `verified` \| `quarantined` \| `refused`, per-tier `TierResult`), `verify_problem(dir)`, `verify_library()`, problem library loader (`verification/problems`), tier-F coverage registry, verification of a mechanism diagnosis, `verify.stations.station_statics / station_deflections` (the values at the result stations recomputed from the finished results without solver code, part of tier D); CLI `benchfea verify` — see [VERIFICATION.md](VERIFICATION.md) |
| `benchfea.bench` | seeded model generators and measured difficulty: `generate(family, seed, difficulty, units=, vertical_axis=, features=) -> Generated(document, metadata)` (families `continuous_beam`, `portal_frame`, `multistorey_2d`, `multistorey_3d`, `braced_frame`, `pipe_rack`, `stair_frame`, `diaphragm_frame`), `difficulty(document) -> measured properties + level`, `feature_flags`, `tier_c_expected`; CLI `benchfea generate`; worked examples in `examples/` — see [DIFFICULTY.md](DIFFICULTY.md) |
| `benchfea.bench.curveballs` | curve balls: `apply(document, name, seed, effect_class=None) -> Perturbed(document, record)` (named, seeded, deterministic changes whose effect is measured on analysis and design; `NotApplicable(reason)` when the model lacks the feature or no target is material), `available()`, `get(name) -> CurveBall`, `applicability(document)`, `ChangeRecord` (`changes` = JSON-pointer `Change(op, path, before, after, index)`, `entities`, `effect_class` `changes-analysis` \| `design-only` \| `creates-mechanism` \| `invalid-model`, measured `effect`, `reference`, `verification`; `apply(base)`, `revert(perturbed)`), `diff`, `apply_changes`, `invert_changes`, `evaluate`, `compare`, `THRESHOLDS`; CLI `benchfea perturb` (exit code 7: not applicable) — see [CURVEBALLS.md](CURVEBALLS.md) |
| `benchfea.bench.tasks` | verified cases and tasks: `build_case(family, curve_ball, seed, difficulty, cache_dir=) -> dict` (`status` `verified` \| `unavailable`, `skips`: the logged attempts of the derivation chain `derived_seed`), `make_task(case, task_type, mode) -> (task, answers)` (public / hidden part; task types `predict` \| `diagnose` \| `repair` \| `design`, modes `tool-allowed` \| `no-solver`), `TaskUnavailable(reason, log)`, `task_types_for`, `describe(document)` (the intended design in words), `extract(quantity, results, design)`, `select_quantities`, `change_direction`, `reference_design` (documented deterministic search), `catalogue_for`, `design_document`, `model_weight`, `design_feasible` — see [BENCHMARK.md](BENCHMARK.md) |
| `benchfea.bench.scoring` | `score_answer(task, answers, answer) -> {"score" in [0, 1], "valid", "reasons", "summary"}` (never raises for a bad answer; repairs and designs are re-solved through `verify_model`), `TOLERANCES`, `credit`, `relative_error`, `entity_match`, `reference_answer`; CLI `benchfea score ANSWER --task TASK [--answers FILE] [--json]` |
| `benchfea.bench.env` | Gymnasium environment `BenchFEAEnv` (id `BenchFEA-v0`; needs the optional `gymnasium`): `reset(seed=, options={family, curve_ball, difficulty, task_type, mode, instance_seed})`, `step(action)` with tool calls `inspect` \| `edit` \| `run_analysis` \| `read_results` \| `submit` as JSON text, reward on `submit` only, truncation on the step budget; `run_episode(env, agent)`, `apply_edit`; CLI `benchfea bench run-agent`; agents in `examples/agents/` |
| `benchfea.bench.splits` | train / held-out split: `is_held_out(family, curve_ball, seed)`, `split_of`, `combinations(split)`, `seeds(...)`; definition in the committed `splits.json` |
| `benchfea.bench.dataset` | `plan(seeds=) -> cells`, `build(root, cells, workers=, cache_dir=) -> manifest` (writes `tasks/{train,heldout}`, `answers/{train,heldout}`, `tasks/manifest.json` with counts and the quarantine log), `load_task`, `load_answers`, `iter_tasks`; CLI `benchfea bench make-dataset`; the built dataset and its card are in `data/` |

## Input format decisions

* **Collections are objects keyed by id** (`nodes`, `supports` [keyed by node
  id], `materials`, `sections`, `members`, `constraints`, `loadCases`,
  `combinations`, `design.groups`, `design.members`). JSON pointers such as
  `/members/B1/section` are stable under insertion/removal, duplicate ids are
  impossible in data and are caught at parse time in text (code
  `id.duplicate`). Ids match `^[A-Za-z0-9_][A-Za-z0-9_.\-]*$`, ≤ 64 chars.
  Only `loadCases/<id>/loads` and `deflectionLimits` are arrays.
* **JSON is canonical, YAML is equivalent**: a YAML file must mean exactly what
  its JSON transcription means. The YAML loader therefore refuses, with
  line/column, everything that has no JSON counterpart or that YAML 1.1 would
  silently reinterpret: non-string keys (`1:`, `no:`, `~:`), duplicate keys,
  merge keys, dates, NaN/Infinity, and numbers that are not plain decimals
  (`010` = octal 8, `0x10`, `1_000`, `1:30` = sexagesimal 90). A dict passed
  from Python is held to the same rules (string keys, finite numbers, JSON
  value types; tuples and numpy numbers are accepted as lists and numbers).
  Files are UTF-8 (a byte order mark is ignored).
* **Discriminated variants** (`sections`, `loads`, `constraints`) use a `type`
  property with `if/then` + `unevaluatedProperties: false`, which gives one
  precise error instead of `oneOf` noise.
* **Units**: every numeric field has an `x-unit` dimension; conversion is
  driven by the schema (`unitmap`), never by per-field code.
* **Analysis vs design data**: everything that changes analysis results lives
  in `members`, `sections`, …, `analysis`, and `combinations/*` (`analysisType`,
  `stiffnessReduction`, `notionalLoad`). The `design` block is design-only and
  is resolved into `Model.design_members`, which the solver must not read
  (tests build a model with and without a heavy `design` block — including
  overrides of E, G, Fy, Fu — and require every analysis field to be
  identical). Consequently the stiffness-reduction options of AISC 360-22 C2.3
  (factor, tau_b, member list) live in `analysis.stiffnessReduction` and are
  switched on per combination; they are not design variables, and a computed
  tau_b always uses the material's Fy, never a design override.
* **Prescribed displacements** are loads (`supportDisplacement`) on DOFs fixed
  by a support, so they are factored and combined like any other load.

## Internal conventions (binding for all modules)

* Internal units: **N, mm, rad** (stress MPa, moment N·mm).
* Node order = file order; node `index` = position; **DOF order per node:
  `ux, uy, uz, rx, ry, rz`**; global DOF = `6*index + k` (`Model.dof(node, dof)`).
* `Member.R` is 3×3 with **rows = local x, y, z in global components**:
  `v_local = R @ v_global`, `v_global = R.T @ v_local`.
* `Member.start/end/length` describe the **flexible length**;
  `offset_i/offset_j` are global vectors node → flexible end. The element DOFs
  at the flexible ends follow the node rigidly:
  `u_end = u_node + θ_node × offset`, `θ_end = θ_node`.
* `Member.k_end_i/k_end_j`: six connection stiffnesses per end in local DOF
  order; `math.inf` (`RIGID`) = no release, `0.0` = full release, `k` = partial
  (spring between the member end and the offset end/node). For truss-type
  members (`truss`, `tensionOnly`, `compressionOnly`) the solver must use an
  axial-only element and ignore `k_end_*` (stored as rotations released).
* Member loads are already in **member local axes**, positions are absolute
  distances from the flexible start (mm), projected loads already scaled:
  `MemberPointLoad(member, kind, vector(3), x)`,
  `MemberDistributedLoad(member, w1(3), w2(3), x1, x2)`.
  `SelfWeight(factor, direction(3) global unit, members)`: line load
  `factor * unit_weight * A` along `direction` on the flexible length.
  For truss-type members, member loads are carried to the end nodes as the
  reactions of a simply supported beam (station shears/moments of that beam).
* Member end forces vs station forces — see CONVENTIONS §4:
  `MemberResult.end_forces[0:6]` at i and `[6:12]` at j are forces **on the
  member**, local axes; `MemberResult.forces[k] = (N, Vy, Vz, T, My, Mz)` on
  the +x face; `forces[-1] == end_forces[6:]`, `forces[0] == -end_forces[:6]`.
* Reactions: forces by supports **on the structure**, global; shape
  `(n_nodes, 6)`.

## Model (`benchfea.model`)

```python
model = build_model("file.json")     # or a dict; raises ModelInvalid(problems)
model.nodes[id]         -> Node(id, index, xyz)
model.supports[node]    -> Support(node, fixed[6], spring[6])
model.materials[id]     -> Material(id, E, G, nu, unit_weight, Fy, Fu, description)
model.sections[id]      -> Section(id, kind, props: SectionProps, shape)   # .A .Iy .Iz .J .Asy .Asz .Cw .Sy .Sz .Zy .Zz .ry .rz .family .dims
model.members[id]       -> Member(id, type, i, j, section, material, start, end, length, R, offset_i, offset_j, k_end_i, k_end_j)
model.constraints[id]   -> Constraint(id, type, master, slaves, normal, dofs)
model.load_cases[id]    -> LoadCase(id, loads=(NodalLoad | MemberPointLoad | MemberDistributedLoad | SelfWeight | SupportDisplacement, ...))
model.combinations[id]  -> Combination(id, factors, analysis_type, purpose, stiffness_reduction, notional: NotionalLoad | None)
model.analysis          -> AnalysisOptions(type, shear_deformation, result_stations, convergence_tolerance, max_iterations, alpha,
                                           equilibrium_tolerance, pivot_tolerance, geometry_tolerance, buckling_combinations,
                                           buckling_modes, reduction_factor, tau_b, reduction_members)
model.design_members[id]-> DesignParams(id, members, length, section, material, values: {name: Resolved(value, source, path, reason)})
model.defaults_used     -> [DefaultUsed(path, value [file units], reason)]
model.warnings          -> [Problem]   (severity "warning")
model.units, model.vertical_axis, model.dof(node, dof), model.n_dof
```

`Asy`/`Asz` may be `None` (no shear deformation in that direction). Shear
deformation is used only when `analysis.shear_deformation` is true.

### Design variables (`DesignParams`)

Names (`benchfea.designparams.NAMES`): `code, method, role, E, G, Fy, Fu, Lx,
Ly, Lz, Lb, Kx, Ky, Kz, braceMajor, braceMinor, braceTorsion, braceTopFlange,
braceBottomFlange, Cb, Cmx, Cmy, sideswayInhibited, An, U,
maxSlendernessCompression, maxSlendernessTension, deflectionLimits`.
`p.Kx` is the value; `p.source("Kx")` is `"member" | "group:<id>" | "global" |
"default" | "bracing" | "material:<id>"`; `p.values["Kx"].path` is the JSON
pointer where the value is (or would be) set.

* Resolution order per variable: design member → groups containing it (highest
  `priority` wins; equal priority with different values is an error) →
  `design.defaults` → built-in default.
* An unbraced length and its bracing list (`Lx`↔`major`, `Ly`↔`minor`,
  `Lz`↔`torsion`, `Lb`↔`topFlange`/`bottomFlange`) describe the same thing:
  the higher level wins; both at the same level is an error. When brace points
  win, the length is `None` with source `"bracing"` and the design module works
  on the segments between brace points (fractions of `length`, ends implied).
* An unbraced length together with its bracing list in the same block, and
  relative brace positions outside (0, 1), are errors wherever they appear
  (also in a group that no design member uses).
* `Cb`, `Cmx`, `Cmy` are a number or `"computed"`.
* `E`, `G`, `Fy`, `Fu` are design variables like the others. Unless overridden
  at member, group or `design.defaults` level they are those of the first
  member's material (source `"material:<id>"`, G derived from E and nu when
  the material gives none; `Fy`/`Fu` are `None` when the material has none,
  with warning `design.missing-strength`). An override is used by the design
  checks only; the analysis keeps the material's values, and an `E` or `G`
  override that differs from the material raises the warning
  `design.stiffness-differs`.
* A design member is one analysis member or an ordered chain of collinear,
  identically oriented, prismatic members; with `autoMembers` (default) every
  member not in an explicit design member is its own design member.

### Defaults-used report

`Model.defaults_used` lists every default applied because the file was silent
(conventions, analysis options, stiffness-reduction settings and member list,
node `z`, derived G, section fallbacks, member type, the reference vector of
the default member orientation, offset axes, combination settings, load
`positionType`/`kind`, self-weight factor/members/gravity direction, diaphragm
normal, `design.autoMembers`, bracing `positionType`, and every design
variable), as `{path, value, reason}` with `value` in file units. Design
values taken from the material (`E, G, Fy, Fu`) are not defaults and are not
listed. Not reported: additive zero/empty/off defaults (zero load components,
no releases, no offsets, `roll = 0`, `w2 = w1`, load extents `start`/`end`,
`projected = false`, group `priority = 0`, no supports/constraints). Results
documents carry it as `defaultsUsed`.

## Contract for `benchfea.solve`

```python
def analyze(model: Model) -> Results
```

* One `CombinationResult` per combination, in internal units. Loads of the
  cases are summed with the factors; **every combination is solved on its own**
  (second-order and tension/compression-only analyses are nonlinear).
* `analysis_type == "secondOrder"`: geometric stiffness from the member axial
  forces, iterated until the relative displacement change is below
  `convergence_tolerance` (at most `max_iterations`); loads multiplied by
  `analysis.alpha` and results divided by it (AISC 360-22 C2.1(d)).
* `combination.stiffness_reduction`: multiply EA, EI, GJ, GAs of
  `reduction_members` (all if `None`) by `reduction_factor`, and the flexural
  stiffness additionally by τ_b (`tau_b`: `"computed"` per AISC 360-22
  Eq. C2-2a/b with P_ns = F_y A, `"none"`, or a number); report the τ_b used.
* `combination.notional`: add at every node a force
  `coefficient × (vertical component of the factored loads at that node,
  member loads lumped statically to the member ends)` along `direction`.
* `analysis.buckling_combinations`: linear buckling,
  `(K + λ Kg) φ = 0`, `buckling_modes` smallest positive λ.
* Constraints by exact transformation (master/slave elimination), not penalty.
* Mechanism / singular stiffness: raise `ModelUnstable` (a `ModelInvalid`) with
  problems of code `stability.mechanism`, path `/nodes/<id>`, naming the DOF.
* Fill `invariants` (global force and moment equilibrium of applied loads +
  reactions, joint residual, work = strain energy for linear cases) against
  `equilibrium_tolerance`; copy `model.defaults_used` into `Results`.

## Solver core (`benchfea.solve`) — as implemented

```python
from benchfea.solve import analyze, ModelRefused, ModelUnstable
results = analyze(model, strict=True)      # Results, internal units
```

| module | content |
|---|---|
| `solve.beam` | `SectionFlex`, `SpanLoads(points, dists)`, `integrate(L, flex, loads, s_start, y_start) -> Integration` (`.at(x, side)`): closed-form integration of the member equations, segment by segment |
| `solve.element` | `build_element(model, member, factors=None, active=True) -> MemberElement` (`k_local`, `G`, `k_global`, `dofs`, `fixed_end(loads)`, `equivalent_nodal_loads(loads)`, `recover(u_nodes, loads, stations) -> MemberRecovery`), `StiffnessFactors(axial, flexural, torsional, shear)`, `cantilever_flexibility`, `release_stiffness`, `rigid_transfer` |
| `solve.loads` | `LoadSet(nodal, prescribed, member)`, `case_loads`, `combine`, `combination_loads`, `applied_resultant`, `member_load_positions` |
| `solve.system` | `Structure(model, elements)`: `T`, `retained`, `free_r`, `fixed_r`, `inactive`, `assemble(matrices=None)`, `stiffness()`, `reduce(K)`, `factorize(K) -> Factor`, `solve(factor, F, prescribed) -> Solution(u, reactions, constraint_forces, F, residual)` |
| `solve.stability` | `diagnose(...) -> MechanismDiagnosis` |
| `solve.invariants` | `model_checks`, `solution_checks`, `superposition_check(model)`, `reciprocity_check(model, pairs=None)` |
| `solve.analysis` | `analyze`, `LinearAnalysis(model, elements=None)` (`factor`, `case`, `combination_loads`, `load_vector`, `solve_loads`, `recover`, `checks`, `result`), `ANALYSIS_TYPES`, `station_layout`, `unsupported_features` |
| `solve.errors` | `ModelRefused(ModelInvalid)`, `MechanismDiagnosis`, `MechanismMode`, `NodeParticipation` |

* **Dispatch**: `ANALYSIS_TYPES[combination.analysis_type](LinearAnalysis, combination) -> CombinationResult`;
  `"firstOrder"` and `"secondOrder"` are registered. A nonlinear module adds its entry and
  removes the corresponding refusal in `unsupported_features`; it reuses
  `Structure.assemble(matrices)` (e.g. geometric stiffness as per-member 12×12
  matrices in nodal global DOFs), `factorize`, `solve`, `LinearAnalysis.recover`
  and `build_element(..., factors, active)`.
* **Refusals**: `ModelRefused` with codes `analysis.not-implemented`,
  `member.release-mechanism`, `section.zero-stiffness`, `load.torque-on-truss`,
  `solve.invariant-failed` (`.results` holds the flagged results);
  `ModelUnstable` (code `stability.mechanism`) carries `.diagnosis`
  (`MechanismDiagnosis.to_dict()`: `classification` `"mechanism"` |
  `"ill-conditioned"`, `nMechanisms`, `modes[]` of kind `rigid-body` | `node` |
  `distributed` with `description, motion, cause, direction, axis, through,
  nodes[{node, share, dofs}], members`). Problem paths: `/supports` for
  rigid-body modes, `/nodes/<id>` otherwise.
  `cause` ∈ `no-restraint`, `negligible-springs` (rigid-body motion opposed only by
  support springs below the tolerance), `all-member-ends-released`, `only-truss-members`,
  `member-spin` (a member free to spin about its own axis), `no-stiffness`.
* **Independent review** (`tests/independent/`, THEORY 17): a reference solver that shares no
  code with the engine (`tableau.py`), direct OpenSeesPy / PyNite translations (`oracles.py`),
  `solve_bvp` single members and cut-and-balance statics; model bank in `models.py`, unstable
  models in `unstable.py`.
* **Stations**: `analysis.result_stations` equally spaced stations plus the
  ends of line loads plus two stations (left limit, right limit) at every
  interior point load of any load case; the same list in every combination.
  `MemberResult.x` therefore may contain an x twice.
* **Additions to `benchfea.results`**: `CombinationResult.constraint_forces`
  (n_nodes, 6) = forces by the rigid constraints on the nodes (document:
  `combinations/<id>/constraintForces/<node>` for masters and slaves);
  `Results.nodal_envelopes` = `{"displacements" | "reactions": NodalEnvelope(max,
  min, max_combination, min_combination)}` (document: `nodalEnvelopes/<kind>/<node>`).
* **Invariants** (names in `CombinationResult.invariants.checks`): `stiffness
  symmetry`, `rigid-body modes`, `positive definiteness`, `linear solve
  residual`, `global force equilibrium`, `global moment equilibrium`, `joint
  equilibrium`, `member equilibrium`, `member compatibility`, `work = strain
  energy` (definitions and tolerances: THEORY.md §8).
* **CLI**: `benchfea run FILE [-o OUT] [--json] [--allow-flagged]`; exit code 0
  ok, 1 model refused (invalid, unsupported, unstable — with the diagnosis in
  `--json`), 2 input unreadable / usage, 3 results written but flagged.
* **Inactive DOFs**: free DOFs of nodes without members (unless a master DOF
  that drives slaves) and rotations of nodes with only truss-type members are
  removed from the system and reported as 0 with a warning in `Results.warnings`.

## Nonlinear analyses (`benchfea.solve`) — as implemented

Theory, rules and limits: [THEORY.md](THEORY.md) sections 10–16.

| module | content |
|---|---|
| `solve.beamcolumn` | `BeamColumn.build(model, member, loads=None, n_start=0.0, factors=None, axial=None)`: the second-order piece (exact for a constant axial force; same attributes as `MemberElement` for assembly: `dofs, k_global, G, nodes, is_truss`; `stable`, `f_fixed`, `equivalent_nodal_loads()`, `recover(u_nodes, stations=(), checks=False) -> PieceRecovery`, `slopes(xi)`), `axial_state`, `consistent_geometric_stiffness(piece, axial)`, `expm` |
| `solve.refine` | `refine(model, divisions, unsplit=()) -> Mesh(model, base, pieces, bounds, n_nodes)` (derived model: original nodes first, internal nodes and pieces with ids `<member>@…`), `split_loads(mesh, loads)`, `piece_stations(bounds, stations)` |
| `solve.nonlinear` | `NonlinearSolver(la)`: `solve(combination) -> CombinationResult`, `loads`, `active_set`, `first_order`, `second_order`, `violations`; `CombinationRefused(problems, diagnosis, info)`, `notional_forces`, `tau_b`, `needs_nonlinear`, `second_order_checks` |
| `solve.buckling` | `buckling_analysis(nonlinear_solver, combination) -> BucklingResult`, `BUCKLING_TOL` |

* **Dispatch**: `ANALYSIS_TYPES = {"firstOrder", "secondOrder"}`. A first-order combination goes through
  `NonlinearSolver` when the model has tension/compression-only members or the combination has stiffness
  reduction or notional loads; otherwise the linear path is unchanged. Every combination is solved on its own.
* **Refusal of one combination** (`CombinationRefused`, collected by `analyze`): codes
  `stability.beyond-critical-load` (with `info["bucklingFactor"]`), `stability.slack-mechanism` (raised as
  `ModelUnstable` with `.diagnosis`), `solve.active-set-cycling`, `analysis.release-second-order`,
  `analysis.tau-b-zero`. `analyze` then raises with `.results` = the results of the other combinations and
  `.refused` = `{what: CombinationRefused}`; with `strict=False` it returns those results and a warning
  `NO RESULTS for …` per refused combination. A combination that does not converge is returned flagged
  (`converged = False`; strict: `solve.not-converged`).
* **Additions to `benchfea.results`**: `CombinationResult.solver` (plain data; document:
  `combinations/<id>/solver`): `secondOrder {iterations, converged, tolerance, maxIterations, history[{iteration,
  displacementChange, axialForceChange, tauBChange}], subdivisions{member: pieces}, largestPieceKh,
  subdivisionCheck{pieces, change, tolerance}}`, `activeSet {iterations, slack[], history[], strategy}`,
  `stiffnessReduction {clause, factor, tauB, members}`, `notionalLoads {clause, coefficient, direction, unit,
  level, nodes{node: [Fx, Fy, Fz]}}`, `tauBIterations`. `BucklingResult.negative_factors`, `.negative_modes`,
  `.member_shapes {member: (n_modes, n_stations, 3)}`, `.info` (document: `buckling/<id>/{factors, modes,
  negativeFactors, negativeModes, memberShapes, info}`; mode rotations are written per file length unit).
  `MemberResult.active` is `False` for a slack member. `CombinationResult.iterations` = second-order
  iterations (or τ_b iterations of a first-order combination), `alpha`, `tau_b` as defined before.
* **Station forces of a second-order result** are components along the undeformed member axes; the shear
  normal to the deformed axis is `Vy − N v′` and `dMz/dx = −Vy + N v′` (THEORY 10.8).
* **Invariants** of a second-order combination: the names of the linear checks plus `member equilibrium on
  the deformed axis`; `rigid-body modes` covers translations only; `global moment equilibrium` includes the
  couples of the member axial forces (THEORY 15).
* **CLI `run`**: per second-order combination the iteration record (converged, iterations, last changes,
  tolerance, alpha, pieces per member), slack members, τ_b range; per buckling combination the factors, the
  factors under reversed loads and the estimated relative error; the same under `secondOrder`, `slackMembers`,
  `tauB` and `buckling` in `--json`.

Limits of this part (for LIMITS.md) — THEORY 16: axial-force effects on flexure only (no lateral-torsional,
torsional or flexural-torsional buckling); small displacements (no bowing, no large rotation); Engesser shear
form; transverse translational releases refused in second-order / buckling analyses; rigid links and
diaphragms without geometric stiffness; truss-type members are pin-ended beam-columns (their own flexural
buckling counts); τ_b from the largest compression of each analysis member, P_ns = F_y A; notional loads from
the signed net downward load at each node; linear buckling on the first-order reference state with frozen
slack set and τ_b; dense matrices.

## Contract for `benchfea.design`

```python
from benchfea import design
design.check(model, results) -> DesignResults                 # every design member, strength combinations
design.check_member(ctx, demands, edition=None, extra_checks=(), conditions=(), incomplete=()) -> DesignMemberResult
design.member_demands(model, results, design_member_id) -> {combination id: Demand}   # maxima + stations in Demand.extra
design.analysis_basis(model, results) -> (basis, warnings, conditions, incomplete)    # how each combination was analysed
design.registry.get(code) -> Edition                          # unknown code: UnknownDesignCode (refusal)
```

Package layout (`src/benchfea/design/`):

| module | role |
|---|---|
| `registry` | `Edition(code, title, limit_states)`, `register`, `get`, `registered`. Only `"AISC 360-22"` is registered; any other `code` is refused (`UnknownDesignCode`; `check_member` then returns one `refused` check). A new edition = a new package with one module per chapter, registered from its `__init__` and imported in `design/__init__.py` |
| `base` | `DesignContext`, `Demand`, the check builders `strength_check`, `limit_check`, `not_checked`, `refused` |
| `section` | `DesignSection`: design view of a section (shape family, rolled / built-up, plate elements with width, thickness, ratio and provenance); `DesignSection.from_model(section)` and the constructors `i_shape`, `rectangular`, `round_hss`, `unsupported` |
| `demands` | `member_stations`, `demand_from_stations`, `member_demands`, `distance_to_zero` (L_v of G5) |
| `aisc360_22/` | `b4_classification` (B4.1, Tables B4.1a/b), `d_tension` (D1–D3), `e_compression` (E2, E3, E4-2, E7), `g_shear` (G2.1, G4, G5, G6), `f_flexure` (F1 Cb, F2, F3, F6, F7, F8), `h_combined` (H1.1, H1.2; torsion screening; App. 8.1.2 formulas), `serviceability` (deflection limits) |

**Limit-state function** (the unit every chapter module is made of; each module
lists its functions in `LIMIT_STATES`, an edition concatenates them):

```python
def limit_state(ctx: DesignContext, demand: Demand) -> list[DesignCheck]
```

* `ctx` — one design member: `member, section: DesignSection, method ("LRFD"|"ASD"),
  E, G, Fy, Fu, L, Lx, Ly, Lz, Kx, Ky, Kz, An, U, maxSlendernessCompression,
  maxSlendernessTension, sources {name: provenance}, code`.
  `DesignContext.from_model(model, id)` fills it from `Model.design_members`
  (AISC axis letters: x = major = local z, y = minor = local y; an unbraced
  length defined by brace points becomes the longest segment between them).
  `ctx.inputs("Fy", "Kx")` → `{name: (value, dimension, source)}`.
* `demand` — the required strengths of ONE combination (maxima over the stations
  of the design member's analysis members): `combination, Pt, x_Pt, Pc, x_Pc, Vy,
  x_Vy, Lv_y, Vz, x_Vz, Lv_z, V, x_V, Lv` (magnitudes; `x_*` from the design
  member start; `V` = resultant shear; `Lv*` = distance from maximum to zero
  shear, `None` if the shear never vanishes) and `extra` (dict) for the demands
  of later limit states (moments).
* It returns at least one check for every limit state it owns and never raises
  for something it cannot handle: outside the scope or data missing →
  `status = "notChecked"`; contradictory input → `status = "refused"`; both with
  `reason`. All equations used are unit-consistent (the engine calls with N, mm).
* `check_member` evaluates every limit state for every combination and keeps,
  per check `name`, the governing one (largest ratio; a check without a result
  for a nonzero demand is never hidden by another combination).

**`DesignCheck` record** (`benchfea.results`, internal units; `to_document`
converts): `name` (limit state), `clause` (every section and equation used,
e.g. `"AISC 360-22 E3, Eq. E3-2, E3-4; E2; B4.1, Table B4.1a; B3.1 Eq. B3-1"`),
`status` (`"ok" | "notChecked" | "refused"`) and `reason`, `method`, `inputs`
`{name: (value, dimension, source)}` with source `member | group:<id> | global |
default | bracing… | material:<id> | section:<id> | given | derived: …`,
`details` (intermediate values `{name: (value, dimension)}`), `nominal` (Rn),
`factor` + `factor_name` (`phi` or `Omega`), `capacity` (available strength,
B3-1 / B3-2), `demand`, `ratio`, `passed`, `combination`, `x`, `dimension`,
`advisory`, `notes`. When `status != "ok"`: `capacity`, `ratio`, `passed` are
`None` (the demand is kept).

`DesignMemberResult.passed` is true only if every non-advisory check passed
and no limit state with a nonzero demand is `notChecked`/`refused`;
`not_checked` lists every check without a result as `"name: status: reason"`;
`governing` = the non-advisory check with the largest ratio. Slenderness
recommendations (D1, E2 User Notes) are `advisory`. Strength checks use the
combinations with `purpose == "strength"`.

**Section → design view** (the input format has no field for how a shape is
made, so the mapping follows what the section data says and is reported in the
inputs/notes of every check): database W/M/S/HP → rolled I, h = d − 2·kdes
(optional `kdes` column of the user's database; absent → web-dependent checks
are `notChecked`); database HSS → rectangular / round HSS with `tdes`;
`type: I` with fillet radius r > 0 → rolled, h = d − 2tf − 2r; r = 0 →
built-up welded, h = d − 2tf; `type: box` → built-up **box section** (clear
widths, Table B4.1a Case 8), not an HSS; `type: pipe` → round HSS with t used
as the design wall thickness (the 0.93 factor of B4.2 is never applied by the
engine); everything else is out of scope (tension checks only, which need
only A). `DesignSection` constructors accept tabulated ratios (`bf/2tf`,
`h/tw`, `b/t`, `h/t`, `D/t`) directly; they then take precedence.

### Flexure, combined forces, serviceability, orchestration (as implemented)

* **Additions to the shared records** (all optional, defaults keep earlier code working):
  `DesignSection.Sx, Zx, Sy, Zy` (+ `dims["rts"]`, `dims["ho"]` when tabulated values are given);
  `DesignContext.Lb, braceTop, braceBottom, Cb, end_free, end_built_in`; `Edition.service`
  (serviceability functions `(model, results, design_member_id, ctx) -> [DesignCheck]`);
  `DesignMemberResult.status, by_combination`; `DesignResults.status, combinations`.
* **Moments in `Demand.extra`**: either the stations (`"x"`, `"forces"` (n, 6), and `"analysisType"`) — what
  `design.member_demands` supplies — or the numbers `"Mrx"`, `"Mry"`, `"Tr"`. A demand without either says
  nothing about flexure: no Chapter F / H check is produced and `check_member` reports it (status
  `"incomplete"`).
* **Chapter F** (`f_flexure`): `moments(demand)`, `Diagram` (cubic through station moments and slopes,
  dMz/dx = −Vy, dMy/dx = +Vz — exact between stations for the loads of the format in a first-order analysis),
  `cb_f1_1`, `segment_cb`, `ltb_segments(ctx, moments)`, `axis_strength(ctx, demand, "x"|"y"|"r") ->
  AxisResult(checks, Mc, blocked, governing)`, pure functions per equation. Major axis = `Mz` (top flange =
  +y side in compression when `Mz > 0`). Round sections use the resultant moment.
* **Chapter H** (`h_combined`): `h1_interaction`, `interaction` (Eq. H1-1a/b at every station AND at the
  stationary points of the interaction between stations — moments by the `Diagram` cubics, axial force
  linear, so the result does not depend on `resultStations` for the loads of the format — with the
  concurrent demands against the member's governing available strengths), `torsion` (H3 not implemented:
  engine screening rule τ_t ≤ 1 % of 0.6 Fy, else `notChecked`), `cm_a_8_4`, `pe1_a_8_5`, `b1_a_8_3`
  (formulas only — B1/B2 are never applied automatically).
* **Deflections** (`serviceability`): `deflection_checks(model, results, id)`: per entry of
  `deflectionLimits`, the largest deflection relative to the chord through the design-member ends over the
  listed (default: all serviceability) combinations against `min(L/spanRatio, absolute)`; enforced check.
* **Status** of a member / of the design (`design.PASS, CONDITIONAL, INCOMPLETE, FAIL`), first that applies:
  `fail` (a non-advisory ratio > 1) → `incomplete` (an applicable limit state `notChecked`/`refused`, a
  strength combination missing or not converged, demands without moments) → `conditional` (all pass, but a
  strength combination was analysed first-order: AISC 360-22 C1 requires second-order effects) → `pass`.
  `passed` keeps its meaning (no failed and no blocked strength check).
* **Analysis basis**: `DesignResults.combinations[cid]` = `{purpose, inResults, analysisType, secondOrder,
  converged, stiffnessReduction, notionalLoad, alpha, directAnalysisSettings}`; warnings for second-order
  without stiffness reduction / notional loads, ASD with alpha ≠ 1.6, and per member in compression for
  K = 1 without direct-analysis settings or K < 1.
* **Design never changes the analysis**: `design.check` only reads `model` and `results` (tests: results
  bitwise identical under any change of Lb, K, bracing, Cb, method, Fy, E in the `design` block; deep
  comparison of model and results before / after `check`).
* **Design document** additions: `status`, `combinations/<id>` (analysis basis),
  `members/<id>/{status, byCombination/<combination id>[checks]}`.
* **CLI**: `benchfea design FILE [-o OUT] [--results FILE] [--json] [--allow-flagged]`: analysis + design,
  one row per design member (status, ratio, combination, governing limit state [clause]); exit code 0
  `pass`/`conditional`, 4 `fail`/`incomplete`, 1 model refused, 2 input unreadable, 3 analysis flagged.

Limits of this part (for LIMITS.md):

* Not implemented → `notChecked`: I-shapes with noncompact / slender webs (F4, F5); rectangular sections
  with slender webs (Eq. F7-7); round HSS with D/t ≥ 0.45 E/Fy; F9–F13; H1.3, H2, H3 (torsion), H4;
  Appendix 8 B1/B2 as an automatic amplification; an I-shape whose local-z inertia is not the larger one;
  a rectangular section bent about local y when that is its strong axis (LTB).
* LTB bracing: only `bracing.topFlange` / `bottomFlange` (member ends are brace points); twist braces
  (`bracing.torsion`) are not counted. With `Lb` given and ≠ member length the brace positions are unknown:
  Cb = 1.0 unless given, with the largest moment of the member. Where the moment changes sign between the
  brace points of one flange, the segment between points where both flanges are braced is used.
* Free ends: F1's cantilever rule (Cb = 1.0, Lb = cantilever length) only for a member that is free at one
  end and built-in (six DOFs fixed) at the other, warping restraint assumed; other free-end segments are
  `notChecked` unless Cb is given. Load height (top-flange loading) is not considered.
* Se of slender HSS / box flanges: ineffective width removed from both flanges (no neutral-axis shift).
* H1: governing member strengths (lowest LTB strength of all segments) with concurrent station demands;
  the optional Cb increase under tension (Eq. H1-2) is not applied; maxima given as numbers are taken as
  concurrent.
* Torsion threshold (engine rule, not AISC): τ_t = T·max(tf, tw)/J (I) or T/(2 A0 t) (closed) ≤ 0.01 · 0.6 Fy.
* Deflection: maximum over the result stations only; members with a free end are `notChecked`;
  slenderness limits are advisory only (no input to enforce them).
* Moment interpolation between stations uses the station shears as slopes; with second-order results this
  is an approximation between stations.

### Known limits of the design checks (for LIMITS.md)

* Shapes: doubly symmetric I-shapes, rectangular HSS, box sections, round HSS /
  pipe. Channels, tees, angles, singly symmetric or unsymmetric shapes,
  `rectangle`, `properties` sections: compression and shear `notChecked`.
* Implemented: B3.1/B3.2; B4.1 (Table B4.1a Cases 1, 2, 5, 6, 8, 9; Table B4.1b
  Cases 10, 11, 13, 15, 17, 19, 20, 21); D1 (advisory), D2, D3-1; E2 (advisory),
  E3, E4-2, E7.1, E7.2; G2.1, G4, G5, G6; flexure, combined forces and deflection
  limits as listed above. Not implemented: D4–D6, E4-3/E4-4/E4-10..12, E5, E6,
  G2.2/G2.3 (tension field action), G2.4, G3, G7, Table D3.1 (U is an input),
  B4.3b net area (An is an input), block shear and all connection limit states.
* An = Ag and U = 1.0 by default (reported as defaults).
* Webs have no transverse stiffeners (kv = 5.34); no web openings; no holes.
* No file field for rolled/built-up or HSS/box (see the mapping above); HSS
  corner radius never known (flat width = outside − 3t).
* Compression with brace points: Lc = K × longest segment, with the maximum
  compression of the member (conservative when the force varies).
* E3 branch by Lc/r (E4 by Fy/Fe); Table E7.1 c2 (rounded) with be ≤ b;
  Eq. E7-7 capped at Ag; D/t ≥ 0.45E/Fy `notChecked`; torsional buckling (E4)
  evaluated where E4's first paragraph calls for it or leaves it open — not
  for a rolled I-shape or an HSS whose torsional unbraced and effective lengths
  do not exceed the lateral ones (`e_compression.e4_scope`; Cw = 0 for closed
  sections).
* A station force below 1e-9 of the member's force scale is zero for the design
  (`design/demands.py`, `ZERO_FORCE`): rounding errors of the analysis do not
  decide signs or "is there a force".
* G5: resultant shear; Eq. G5-2a not credited when the shear never vanishes
  along the member (L_v undefined). G6: 2 × Eq. G6-1 for the two flanges.
  No interaction between shear directions.
* Safety factors are the printed 1.67 / 2.00 / 1.50 (not 1.5/φ); no tolerance
  at classification boundaries (ratio = limit is the favourable class).

## Output documents

Key-by-key description for users: [OUTPUT_REFERENCE.md](OUTPUT_REFERENCE.md).

`Results.to_document(model)` → format `benchfea-results` 1.0.0;
`DesignResults.to_document(model)` → `benchfea-design` 1.0.0. Both are plain
JSON data in the **file's units** with `units`, `defaultsUsed` and `warnings`.
Results layout: `combinations/<id>/{displacements/<node>, reactions/<node>,
members/<id>/{endForces/{i,j}, stations/{x,N,Vy,Vz,T,My,Mz,ux,uy,uz}},
invariants, converged, iterations, alpha, tauB}`, `buckling/<id>/{factors,
modes}`, `envelopes/<member>/{x,max,min,maxCombination,minCombination}`.

## Validation problem codes

`parse.{json,yaml,io,extension,encoding,key-type,value}`, `id.duplicate`, `id.invalid`, `version.unsupported`, `schema.*`
(`required`, `unknown-property`, `enum`, `type`, `range`, `invalid`, …),
`ref.unknown-{node,section,material,member,case,combination,design-member}`,
`member.{same-node,zero-length,offset-zero-length,orientation-collinear,
release-mechanism,release-on-truss,near-vertical}`, `section.{geometry,
shapes-db}`, `constraint.{master-is-slave,slave-two-masters,chained,
slave-restrained,diaphragm-out-of-plane}`, `load.{position-range,
projected-local,zero-direction,displacement-dof-not-fixed,self-weight-zero,
empty}`, `combination.{notional-direction,reduction-first-order}`,
`analysis.tau-b-needs-fy`, `design.{member-reused,id-collision,not-connected,
not-collinear,orientation-differs,not-prismatic,group-conflict,
length-and-bracing,bracing-range,missing-strength,stiffness-differs,asd-alpha}`,
`node.{coincident,unused}`, `support.{none,no-restraint}`, `material.g-nu-inconsistent`, `case.unused`,
and (solver) `stability.mechanism`.

## Known limits of this layer (for LIMITS.md)

* Section properties: J of I-shapes (Σbt³/3) and boxes (Bredt), Cw of I-shapes
  (thin-walled), all shear areas (simple definitions, not Cowper coefficients)
  are approximations — see the table in `benchfea/sections.py`; `overrides`
  replaces any of them. Box sections have sharp corners.
* No distributed moments, no thermal loads, no member loads on rigid offsets,
  no tapered members, no warping torsion in analysis, no chained constraints,
  no inclined supports.
* Computed tau_b uses P_ns = F_y A (gross area). AISC 360-22 C2.3 defines
  P_ns = F_y A_e for slender-element sections; the effective area is not
  evaluated at this layer.
* Validation refuses a structure with no support at all (`support.none`);
  every other mechanism is found by the solver (`stability.mechanism`).
* `type: properties` sections carry no shape dimensions (design checks that
  need them must refuse).
* The AISC Shapes Database reader's column names and units are unverified
  assumptions (overridable); no AISC data is in the repository.
