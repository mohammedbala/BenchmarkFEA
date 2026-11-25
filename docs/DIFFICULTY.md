# Difficulty levels and model generators

`benchfea.bench` generates benchmark models and measures how hard a model is.

```python
from benchfea.bench import generate, difficulty, families, levels

g = generate("braced_frame", seed=3, difficulty=3)   # Generated(document, metadata)
g.document                    # the input document (dict): valid, stable, sized
g.metadata                    # parameters drawn, feature flags, measured difficulty, design ratios, ...
g.json()                      # canonical JSON text
difficulty(g.document)        # measured properties, score and level of ANY input document
```

## 1. What "difficulty" means here

A level is a statement about **measurable properties of the model**, not about how hard a
particular question on it is. `difficulty(document)` measures six properties, gives each a number
of points, and maps the sum (the *score*) to a level 1-5. Nothing is taken from the generator's
intentions: the same function measures a generated model, a hand-written one, or one that a curve
ball has changed.

| component | measured property | points |
|---|---|---|
| `size` | free degrees of freedom: not fixed by a support, not eliminated by a rigid constraint, not inactive; for a planar model only the three in-plane DOFs of each node count. ≤ 12 / ≤ 30 / ≤ 72 / ≤ 150 / more | 0 / 1 / 2 / 3 / 4 |
| `indeterminacy` | degree of static indeterminacy. 0 / 1-3 / 4-12 / 13-48 / more | 0 / 1 / 2 / 3 / 4 |
| `dimension` | the model is three-dimensional (not planar) | 0 / 1 |
| `features` | number of feature kinds present, at most 3 counted: member end releases (full or partial); rigid end offsets; rigid constraints (links, diaphragms); tension-/compression-only members; support springs or imposed support displacements | 0-3 |
| `secondOrder` | no second-order combination / some / some and the largest amplification ≥ 1.10 | 0 / 1 / 2 |
| `combinations` | number of load combinations: ≤ 2 / 3-5 / more | 0 / 1 / 2 |

| score | 0-2 | 3-5 | 6-8 | 9-11 | 12-16 |
|---|---|---|---|---|---|
| **level** | 1 | 2 | 3 | 4 | 5 |

Definitions (code: `benchfea/bench/difficulty.py`; the limits are the module constants
`SIZE_LIMITS`, `INDETERMINACY_LIMITS`, `LEVEL_LIMITS`, `AMPLIFICATION_BAND`):

* **Planar** means: all nodes lie in one plane normal to a global axis, every member has a
  principal axis along that normal, every load acts in the plane and every rigid offset lies in it.
  A plane frame that is modelled in 3D space with its out-of-plane DOFs restrained is planar.
* **Degree of static indeterminacy** = unknown forces − equilibrium equations. Unknowns: 6 per frame
  member (3 planar) minus its fully released end DOFs, 1 per truss-type member, 1 per support DOF
  that is fixed or sprung, 1 per DOF of a slave node that a rigid constraint eliminates. Equations: 6
  per node (3 planar); 3 (2 planar) at a node joined only by truss-type members; at a memberless
  master only the DOFs that drive slaves. For a stable model this is the dimension of the
  self-stress space. A tension-only member counts as a member (the model as posed, not the active
  set of one combination). Rigid constraints raise the count: each eliminated slave DOF brings a
  constraint force.
* **Amplification** of a second-order combination = largest *horizontal* nodal translation of the
  combination as analysed ÷ the same with the combination switched to first order (the sway
  amplification; vertical deflections would dilute it). A combination whose first-order horizontal
  translations are below 1 % of its largest translation (symmetric gravity loading) uses the largest
  translation of any direction instead. The component uses the largest value over the second-order
  combinations. It needs two analyses; `difficulty(document,
  analyse=False)` skips them and returns `level: None` with `levelRange` when the answer depends on
  it.
* The **features** component is capped at 3 so that a small model cannot reach a high level on
  features alone.

`difficulty()` returns every measured value (`freeDofs`, `staticIndeterminacy`, `planar`,
`featureKinds`, `nCombinations`, `nSecondOrderCombinations`, `amplification`, ...) together with
`components` (the points), `score` and `level`.

Measured levels of the worked examples (`examples/`): 01 simple beam — score 1, level 1; 02
continuous beam — 4, level 2; 03 portal frame — 4, level 2; 04 diaphragm frame — 7, level 3; 05
braced frame — 7, level 3.

What the level does **not** measure: the difficulty of the task posed on the model (the curve ball),
the amount of design input, numerical conditioning, or how hard the model is for the foreign solvers
of the cross-check.

## 2. Families

`generate(family, seed, difficulty, units=None, vertical_axis=None, features=None)`.

| family | what it is | levels | optional features (flags) |
|---|---|---|---|
| `continuous_beam` | beams over 1-7 spans, pattern live load | 1-3 | `hinge`, `partialFixity`, `springSupport`, `settlement`, `overhang`, `varyingLoad`, `pointLoads` |
| `portal_frame` | single-bay portal, flat or pitched | 1-3 | `pitched`, `secondOrder`, `baseSprings`, `semiRigid`, `rigidOffsets`, `craneLoad`, `buckling` |
| `multistorey_2d` | multi-bay multi-storey moment frame, plane model in 3D space with out-of-plane restraint | 2-4 | `secondOrder`, `pinnedBay`, `baseSprings`, `rigidOffsets`, `buckling` |
| `multistorey_3d` | space moment frame on a rectangular grid | 3-5 | `secondOrder`, `pinnedSecondary`, `baseSprings`, `buckling` |
| `braced_frame` | braced bay (X, chevron or single diagonal) in a pinned multi-storey frame, plane model | 2-4 | `tensionOnly`, `secondOrder`, `baseSprings`, `buckling` |
| `pipe_rack` | transverse moment bents, longitudinal struts and bracing, 1-2 tiers (3D) | 3-5 | `tensionOnly`, `secondOrder`, `baseSprings`, `buckling` |
| `stair_frame` | inclined stringers and landings: plane stair (levels 1-2), dog-leg stair in 3D (levels 3-4) | 1-4 | `posts`, `hingedFlights`, `springSupport`, `frontPosts`, `buckling` |
| `diaphragm_frame` | space frame with a rigid diaphragm per floor, master at an eccentric mass centre, lateral load with torsion | 3-5 | `secondOrder`, `accidentalTorsion`, `bracedSide`, `pinnedSecondary`, `baseSprings`, `buckling` |

`features={"flag": True | False}` forces a feature on or off; a flag left out is drawn from the seed
with a probability that grows with the level. Forcing a flag can make a level unreachable
(`GenerationFailed`). The unit system (`"kN-m"`, `"kN-mm"`, `"kip-in"`) and the vertical axis
(`"Y"`, `"Z"`) are drawn from the seed unless given; a family describes its structure in abstract
axes and the same structure is written in either convention.

Sections are generic parametric shapes from a small catalogue of round dimensions defined in
`families/_common.py` (`I-300x150x10x7`, `Box-150x150x6`, `Pipe-114x5`, ...). They are **not**
rolled shapes of any standard; no shape table is used. Loads are drawn from ranges chosen for
sensible proportions, and the combination factors are model data — they are not a statement of a
loading code.

## 3. How a model is made, and what the generator ensures

1. **Deterministic draws.** Attempt *k* of `(family, seed, difficulty)` draws from
   `numpy.random.Generator(PCG64(s))`, with *s* the first 128 bits of the SHA-256 of
   `"benchfea.bench|<generator version>|<family>|<seed>|<difficulty>|<k>"` — never Python's `hash()`.
   Only the raw 64-bit stream of PCG64 is used, turned into numbers by plain integer / IEEE
   arithmetic. Coordinates, loads and stiffnesses are rounded to engineering values (0.25 m / 0.5 ft
   grids, 0.5 kN, 0.5 kN/m, two significant digits for springs). The same call gives the same JSON
   text on every run; `GENERATOR_VERSION` is bumped when a change alters any output.
2. **Level before analysis.** The difficulty that needs no analysis is measured on the draw; a draw
   that would not give the requested level — or would give it only on one side of the amplification
   band — is dropped at once.
3. **Sizing.** Each sizing group (columns, beams, braces, ...) moves along the catalogue until its
   largest design ratio (AISC 360-22 member checks of `benchfea.design`, deflection limits included)
   is between `RATIO_LOW` = 0.50 and `RATIO_HIGH` = 0.95; first on the model with every combination
   first order, then on the complete model, where groups only go up.
4. **Acceptance** — the same rules for every family, never relaxed for a model:
   * the analysis ran and every invariant of every combination passed;
   * the design status is `pass` or `conditional`: no failed check, no design ratio above
     `RATIO_MAX` = 1.0, and no applicable limit state that the design module cannot check (for
     example torsion above its screening threshold);
   * at least half of the design members (`RATIO_BAND_FRACTION`) have their governing ratio inside
     `RATIO_BAND` = 0.3 … 1.1;
   * **stability** of the gravity combination `ULS1`: if the file asks for a linear buckling
     analysis (feature `buckling`), its factor is at least `BUCKLING_MIN` = 2.0; otherwise a
     second-order analysis of `ULS1` with all loads multiplied by 2.0 must be below the elastic
     critical load (the engine refuses it exactly when the tangent stiffness is no longer positive
     definite). `continuous_beam` has no axial force and no stability check;
   * second-order amplification at most `AMPLIFICATION_MAX` = 1.5;
   * the **measured** level equals the requested one.

   A draw that fails goes to attempt *k* + 1 (at most 60); the metadata lists the reasons of the
   rejected attempts (`rejected`).

Three features are **opt-in** and never drawn by default, because a model that has them cannot pass
the cross-check (tier C): `buckling` (a linear buckling analysis of `ULS1` requested in the file)
and `tensionOnly` (tension-only X-bracing in `braced_frame` and `pipe_rack`), which no foreign
solver reproduces, and `secondOrder` in `diaphragm_frame` (PyNite emulates a diaphragm with stiff
members and does not follow the second-order result: measured 0.5-0.7 % away from the engine and
OpenSees, which agree with each other to 1e-4). Ask for them with `features={"buckling": True}` /
`features={"tensionOnly": True}`; the rest of the model is the same as without the flag, apart from
what the flag itself changes (`buckling`: nothing else; `tensionOnly`: bracing pattern X, a lateral
load in every combination, first-order analysis and a buckling analysis in the file — section 5).

## 4. Metadata

`g.metadata` holds: `generator` (version), `family`, `seed`, `difficulty`, `attempt`, `rejected`,
`units`, `verticalAxis`, `parameters` (everything the family drew: geometry, loads, options),
`requestedFeatures`, `features` (what the document actually contains, measured by
`feature_flags(document)`), `measured` (the full output of `difficulty`), `sizing` (iterations,
analyses, section of each group), `stability` (which check was made, the buckling factor if one was
computed), `design` (status, number of design members, smallest / median / largest governing ratio,
fraction inside 0.3-1.1, limit states not checked), `bucklingFactors`, `invariantsPassed` and
`tierC_expected`.

### `tierC_expected`: can the model be cross-checked?

Verification tier C needs **both** foreign solvers (OpenSeesPy and PyNite) to represent a model
exactly or by a converged approximation (docs/CROSSCHECK.md section 4). The families stay inside
that envelope by construction:

* no shear deformation (`analysis.shearDeformation: false` — PyNite has none);
* no notional loads and no stiffness reduction (a minimum lateral load is an ordinary load case);
* rigid end offsets only in models without second-order combinations (second order with offsets is
  not eligible);
* connection springs only in bending (`rz` / `ry`); full releases only in bending;
* diaphragm slaves in the plane of their master, master off every grid line; no node on the axis of a
  member it does not belong to;
* linearly varying, partial and projected loads, point moments, support springs and settlements are
  used freely (both translators represent them exactly).

The exceptions are the opt-in features: **tension-only bracing** (no foreign solver reproduces
the active-set iteration), a **buckling analysis in the file** (no foreign solver computes linear
buckling) and **second order with a rigid diaphragm** (PyNite not eligible). A model with any of
them has `tierC_expected.crossCheckable = false`. `tierC_expected` is
derived from the capability table and the features present; it lists per combination and per solver
`exact | approximate | not-eligible | unsupported`, `buckling: none | unsupported`, and notes what
applies. It is an expectation, not a result: only running the cross-check decides — in particular
for second-order combinations, where both foreign solvers are converged approximations and the
cross-check passes only if their estimated errors are small enough.

## 5. Limits

* **Tension-only bracing (opt-in) is first order only**, and such a file always asks for a linear
  buckling analysis of `ULS1`, which serves as its stability check. Reason: the engine refuses the
  second-order analysis of a pinned frame with tension-only X-bracing under gravity plus a small
  lateral load ("beyond the elastic critical load", reported factor about 0.4 for
  `examples/05_braced_frame.json` at its own loads) although the linear buckling factor of the same
  combination on the working diagonals is 7.65. Until that is explained, second order and
  tension-only members are not combined by the generators.

* **Second-order models and the cross-check.** `tierC_expected` says "approximate" for a
  second-order combination: both foreign solvers converge to the beam-column problem. In the sample
  run when the generators were written, every first-order generated model that was put through
  `benchfea.verify.verify_model` passed tiers C, D and E, and every second-order one was quarantined
  by tier C by a small margin: the difference between the engine and OpenSees (1e-4 … 2e-4 relative)
  was a little larger than the error OpenSees was estimated to have, while a refinement study showed
  OpenSees converging towards the engine's value. Until the error estimate of the cross-check is
  settled, `features={"secondOrder": False}` gives models that pass.
* The level is a property of the model, not of the task (section 1).
* "Stable" is established for the gravity combination `ULS1` at twice its loads (or by its buckling
  factor); the other combinations are covered by the member checks and by the amplification limit,
  not by a buckling analysis of their own.
* Design ratios refer to the checks `benchfea.design` implements; a design status `conditional`
  (strength combinations analysed first order) is normal for models without second-order
  combinations.
* Generation time is dominated by the engine (analysis, member checks, buckling): well under a
  second for a level-1 beam, tens of seconds for mid-size frames, minutes for the largest 3D models
  or when a buckling analysis is requested.
* Combination factors, loads and proportions are generic model data chosen for plausible sizes;
  they do not implement a loading standard.
