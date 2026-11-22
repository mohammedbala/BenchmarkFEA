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
                 ├─ benchfea.solve.analyze(model) ───────────► Results          [to be implemented]
                 │       └─ Results.to_document(model) ► "benchfea-results" 1.0.0, file units
                 └─ benchfea.design.check(model, results) ───► DesignResults    [to be implemented]
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
| `benchfea.cli` | `benchfea validate \| run \| design \| perturb \| verify \| score` (only `validate` implemented so far) |
| `benchfea.solve`, `benchfea.design`, `benchfea.verification`, `benchfea.env` | **to be implemented** against the contracts below |

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

## Contract for `benchfea.design`

```python
def check(model: Model, results: Results) -> DesignResults
```

One `DesignMemberResult` per design member with a `DesignCheck` per limit
state: `name`, `clause` (e.g. `"AISC 360-22 F2.2, Eq. F2-2"`), `demand`,
`capacity`, `ratio`, `passed`, governing `combination` and location `x`,
`details` (intermediate values with dimensions) and `notes`. Strength checks
use combinations with `purpose == "strength"`; deflection limits use their own
combination lists. Limit states outside the implemented scope are listed in
`not_checked`, never silently skipped. Slenderness limits are `advisory`.

## Output documents

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
