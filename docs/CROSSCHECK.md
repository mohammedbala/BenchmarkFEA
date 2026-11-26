# Cross-checks against independent solvers, and the ODE oracle

Two packages give benchfea references that do not share its numerical path:

* `benchfea.oracle` — single members and continuous beams by **direct integration of the governing
  ODEs** (`scipy.integrate.solve_bvp`). No stiffness matrix, shape function or fixed-end force is
  involved. Used where a closed form exists (verification tier A) and to test the translators.
* `benchfea.crosscheck` — translators from the benchfea `Model` to **OpenSeesPy** and **PyNite**,
  a subprocess runner, a comparison in benchfea conventions and the automatic decision whether a
  model counts as cross-checked (verification tier C).

Both foreign packages are optional (`pip install benchfea[crosscheck]`); nothing in benchfea imports
them at import time, no line of their source is in this repository, and every test that needs them
skips when they are absent. Frame3DD and MASTAN2 could not be run in this project (no binary /
no MATLAB); Frame3DD appears in the capability table as `not-run`.

Statements below marked **measured** were obtained by numerical experiment with OpenSeesPy 3.8.0 and
PyNiteFEA 3.2.1 and are reproduced by `tests/test_crosscheck_measured.py`; if such a test fails after
an upgrade, the behaviour changed and the corresponding translator rule must be reviewed.
Statements marked **doc** come from the official documentation pages listed at the end.

## 1. Use

```python
from benchfea.model import build_model
from benchfea.solve import analyze
from benchfea.crosscheck import run, compare, cross_check, from_benchfea

model = build_model("frame.json")
r = run("opensees", model, "LC1")            # ForeignResult: displacements, reactions, end_forces, fidelity, info
r.fidelity.level                             # "exact" | "approximate" | "unsupported"
report = cross_check(model, "LC1", benchfea_result=analyze(model))
report.cross_checked, report.reasons, print(report.summary())
```

`ForeignResult` is in benchfea conventions and internal units (N, mm, rad): nodal arrays `(n_nodes, 6)`
in model node order, reactions = forces by the supports on the structure, `end_forces[member]` = 12
forces on the member in local axes at the ends of its flexible length. Results of a second-order
combination are already divided by `alpha`; **second-order end forces are components along the
undeformed local axes** (so that they are in equilibrium with the reactions).

Other entry points: `run_job(job)` (raw worker call), `convergence_study(solver, model, combo,
"subdivisions", [2, 4, 8, 64])`, `equilibrium(model, result)` (applied loads from their definitions +
reactions; no solver involved), `applied_resultant`, `CAPABILITIES` / `capability(solver, feature)`.

### 1.1 Subprocess

`run` writes the job (plain JSON: the list of library calls and the outputs to read) to a temporary
file and executes `python -m benchfea.crosscheck.worker job.json out.json` with a timeout
(default 300 s). A crash, a timeout, a missing package or a failed analysis come back as
`ForeignResult(ok=False, error=...)`; nothing raises into the caller. Reason: OpenSees terminates the
interpreter on some input errors (**measured**: a `vecxz` parallel to the member axis; a model with a
zero-length element). `in_process=True` runs the solver in the calling process (used for PyNite in the
tests, which is pure Python).

### 1.2 Fidelity record

Every translation returns, per model feature, one of

| level | meaning |
|---|---|
| `exact` | the foreign model is the same mathematical problem (differences are round-off) |
| `approximate` | the foreign model converges to the same problem; `error` = estimated relative error |
| `unsupported` | the solver cannot represent the feature or the translator refuses to guess; any numbers returned must not be used as a check |

The `error` of an approximate feature is computed by the runner from runs of the **foreign solver
alone** (`info["convergence"]`); the benchfea result never enters it. Items that share one refinement
parameter share one estimate and are counted once in the sum.

### 1.3 Error bound of an approximate feature

**Discretisation parameter** (`subdivisions`, elements per member, second order). The model is solved
at three levels n/4, n/2, n (OpenSees 8, 16, 32; PyNite 2, 4, 8); the finest is the delivered result.
With d₁ = |r(n/4) − r(n/2)| and d₂ = |r(n/2) − r(n)| (largest relative difference over all compared
quantities) and r(∞) the limit of the sequence,

```
|r(n) − r(∞)|  ≤  d₃ + d₄ + …  ≤  d₂ (ρ + ρ² + …)  =  d₂ ρ / (1 − ρ)
```

if every further level-to-level difference contracts at least by ρ. `runner.estimate_from_levels`
takes

* ρ = max(d₂/d₁, 2^−p_min): the contraction observed, but never better than that of the lowest order
  p_min the method is known to have (1 for both routes). An error a/n² + b/n contracts by almost 1/4 on
  coarse levels and by 1/2 in the limit; an observed contraction that is better than the guaranteed one
  must therefore not be extrapolated. For p_min = 1 the bound is **the last difference d₂ itself**
  whenever the observed order is between 1 and p_max;
* max(d₂, d₁ · 2^−p_max) in place of d₂: a contraction faster than the highest order the method can
  have (p_max = 2 for the chord-rotation geometric stiffness of OpenSees, 4 for the nodal displacements
  of PyNite's cubic element) is a cancellation between two levels, not convergence;
* **no estimate** when d₂ ≥ d₁: the three levels are not in the asymptotic range. A solver without an
  estimate is *not eligible* (reason "approximation without an error estimate … not-asymptotic") — never
  a pass;
* both differences below 10⁻⁹ (three decades under the comparison tolerance): round-off, the estimate is
  the larger of the two.

This replaces the former Richardson estimate |fine − coarse| / (ratio² − 1), which was the leading-order
value of a pure second-order error — not a bound, and too small as soon as a first-order term is present
(section 3.4). Against a pure second-order error the bound is 3× the true error, against a pure
first-order error 1×, against fourth-order convergence 15×.

**Penalty stiffness / stiff emulation member** (`penalty`, `stiffness`): unchanged — two levels with
ratio 10 and `2 × |fine − coarse|`; no error law is assumed because round-off grows with the stiffness
(measured).

**Acceptance** (`tests/crosscheck_acceptance.py`, run by `tests/test_crosscheck_estimate.py`). On models
with an answer that depends neither on the engine nor on a foreign solver — the ODE oracle for member
lines, the tableau solver of `tests/independent` for frames (run at two refinements of its own where the
axial force varies) — the estimated error of each foreign solver must not be below its true error.
Measured, `python3 tests/crosscheck_acceptance.py` (estimates and true errors relative to the model-wide
scales of section 1.4; "order" = observed order of the three levels):

| model | reference | OpenSees estimate | true error | est/true | order | PyNite estimate | true error | est/true | order |
|---|---|---|---|---|---|---|---|---|---|
| derived/beam-column-02 | oracle | 6.01e-04 | 2.00e-04 | 3.0 | 2.00 | 1.45e-04 | 9.41e-06 | 15.4 | 4.01 |
| derived/beam-column-05 | oracle | 2.39e-03 | 8.00e-04 | 3.0 | 1.98 | 5.38e-04 | 3.48e-05 | 15.5 | 3.99 |
| derived/beam-column-08 | oracle | 9.48e-03 (not eligible) | 3.19e-03 | 3.0 | 1.94 | 2.00e-03 | 1.32e-04 | 15.1 | 3.90 |
| derived/beam-column-tension-05 | oracle | 3.63e-04 | 1.19e-04 | 3.0 | 2.02 | 1.17e-04 | 7.52e-06 | 15.6 | 3.94 |
| aisc-commentary/case-1-si-without-shear-p1330 | oracle | 7.75e-04 | 2.58e-04 | 3.0 | 2.00 | 1.84e-04 | 1.20e-05 | 15.3 | 4.01 |
| aisc-commentary/case-1-si-without-shear-p2000 | oracle | 1.39e-03 | 4.63e-04 | 3.0 | 1.99 | 3.19e-04 | 2.08e-05 | 15.3 | 4.00 |
| aisc-commentary/case-1-si-without-shear-p667 | oracle | 3.35e-04 | 1.12e-04 | 3.0 | 2.00 | 8.21e-05 | 5.38e-06 | 15.3 | 4.01 |
| aisc-commentary/case-1-us-without-shear-p150 | oracle | 3.35e-04 | 1.12e-04 | 3.0 | 2.00 | 8.20e-05 | 5.37e-06 | 15.3 | 4.01 |
| aisc-commentary/case-1-us-without-shear-p300 | oracle | 7.77e-04 | 2.59e-04 | 3.0 | 2.00 | 1.84e-04 | 1.20e-05 | 15.3 | 4.01 |
| aisc-commentary/case-1-us-without-shear-p450 | oracle | 1.39e-03 | 4.63e-04 | 3.0 | 1.99 | 3.19e-04 | 2.08e-05 | 15.3 | 4.00 |
| aisc-commentary/case-2-si-without-shear-p445 | oracle | 2.94e-04 | 9.79e-05 | 3.0 | 2.00 | 6.85e-06 | 1.93e-06 | 3.5 | 3.99 |
| aisc-commentary/case-2-si-without-shear-p667 | oracle | 5.81e-04 | 1.94e-04 | 3.0 | 2.00 | 1.84e-05 | 3.84e-06 | 4.8 | 3.98 |
| aisc-commentary/case-2-si-without-shear-p890 | oracle | 1.14e-03 | 3.80e-04 | 3.0 | 1.99 | 4.57e-05 | 7.58e-06 | 6.0 | 3.98 |
| aisc-commentary/case-2-us-without-shear-p100 | oracle | 2.93e-04 | 9.78e-05 | 3.0 | 2.00 | 6.83e-06 | 1.93e-06 | 3.5 | 3.99 |
| aisc-commentary/case-2-us-without-shear-p150 | oracle | 5.80e-04 | 1.94e-04 | 3.0 | 2.00 | 1.84e-05 | 3.84e-06 | 4.8 | 3.98 |
| aisc-commentary/case-2-us-without-shear-p200 | oracle | 1.14e-03 | 3.79e-04 | 3.0 | 1.99 | 4.55e-05 | 7.55e-06 | 6.0 | 3.98 |
| bank/S01 | tableau | 3.72e-04 | 1.24e-04 | 3.0 | 1.99 | 1.38e-05 | 3.70e-06 | 3.7 | 3.18 |
| bank/S03 | tableau | 7.16e-05 | 2.39e-05 | 3.0 | 1.99 | 1.60e-06 | 9.54e-07 | 1.7 | 4.02 |
| bank/S04 | tableau | 6.18e-04 | 2.11e-04 | 2.9 | 1.98 | 4.62e-05 | 1.34e-05 | 3.4 | 3.87 |
| bank/S06 | tableau | 6.63e-05 | 2.21e-05 | 3.0 | 1.99 | 1.72e-03 | 1.53e-03 | 1.1 | 0.97 |
| bank/S07 | tableau | 7.69e-04 | 3.46e-04 | ≥ 1.8 | 1.97 | 1.75e-02 (not eligible) | 5.67e-03 | 3.1 | 2.08 |
| bank/S08 | tableau | 1.93e-04 | 6.44e-05 | 3.0 | 1.99 | unsupported (3.5) | 8.60e-02 |  |  |
| generated portal_frame-3-1 | tableau | 6.40e-05 | 2.11e-05 | 3.0 | 2.01 | 2.28e-04 | 7.86e-06 | 29.0 | 7.95 |
| generated portal_frame-3-2 | tableau | 2.57e-04 | 8.09e-05 | 3.2 | 1.52 | 7.62e-04 | 3.55e-04 | 2.1 | 2.48 |
| generated multistorey_2d-3-1 | tableau | 2.69e-04 | 8.99e-05 | 3.0 | 1.98 | 6.58e-06 | 8.18e-07 | 8.0 | 3.83 |
| generated multistorey_2d-4-2 | tableau | 5.97e-05 | 1.99e-05 | 3.0 | 1.99 | 3.80e-06 | 1.18e-06 | 3.2 | 4.18 |
| generated pipe_rack-3-1 | tableau | 2.20e-04 | 7.41e-05 | 3.0 | 1.99 | 3.02e-05 | 1.03e-05 | 2.9 | 1.34 |
| generated braced_frame-3-1 | tableau | 2.38e-05 | 7.96e-06 | 3.0 | 2.00 | 2.90e-06 | 8.47e-07 | 3.4 | 1.95 |
| generated multistorey_3d-4-1 | tableau | 1.29e-04 | 4.30e-05 | 3.0 | 2.00 | 2.04e-05 | 4.39e-06 | 4.6 | 3.01 |
| generated diaphragm_frame-4-1 | tableau | 5.47e-05 | 1.83e-05 | 3.0 | 1.98 | unsupported (3.5) | 8.09e-03 |  |  |

(Reference uncertainty: 10⁻⁸ for the oracle; for the tableau zero where the axial forces are constant,
at most 7·10⁻⁷ on the generated frames, 5·10⁻⁸ on S06 and 9·10⁻⁵ on S07, whose columns are close to
buckling under their own weight — its est/true is therefore given as a lower limit. S08 and
diaphragm_frame: PyNite's "true error" is its distance from the linear-constraint statement, section
3.5.)

No row has an estimate below the true error. est/true is 3 for OpenSees throughout (clean second order);
for PyNite it is between 1 and 16 except where the cancellation guard acts (portal_frame 3/1: observed
order 8, estimate d₁/16). The eligibility limit (`max_approximation` = 5·10⁻³) is unchanged: at 0.8 of
the Euler load (`beam-column-08`) the OpenSees bound is 9.5·10⁻³ and OpenSees is not eligible; bank
model S07 leaves PyNite not eligible (1.7·10⁻²).

**Not vacuous.** The engine is compared with every eligible solver, so the smaller of the two estimates
decides what a wrong engine result can hide behind. A relative error of 10⁻³ put on the largest
displacement of a second-order result is rejected by tier C on every model of the test
(`test_corrupted_engine_result_does_not_pass_tier_c`); in general it is rejected whenever
10⁻⁶ + estimate + true error < 10⁻³ for one eligible solver — which holds on every row of the table
that has two eligible solvers (largest value of the smaller estimate: 5.4·10⁻⁴ with a true error of
3.5·10⁻⁵, `beam-column-05`).

**Cost.** Tier C on generated level-4 frames, single BLAS thread, empty cache: `multistorey_2d` (35
members, 4 second-order + 2 first-order combinations) 7.1 s before → 19.0 s; `multistorey_3d` (21
members, 3 + 1 combinations) 3.4 s → 8.7 s. Per second-order combination: OpenSees 8/16/32 elements
1.25 s (0.67 s before), PyNite 2/4/8 sub-members with the iteration 3.4 s. Foreign results are cached
per (document, combination, solver, options, solver version, translator revision) within a verification
and, when a cache directory is given, across them.

### 1.4 Tier-C decision (`cross_check`)

A combination is **cross-checked** only if

1. at least two foreign solvers ran, have no `unsupported` feature, and every `approximate` feature
   has an estimate whose sum is ≤ `max_approximation` (default 5·10⁻³) — these are the *eligible* solvers;
2. all eligible solvers agree pairwise, and with the benchfea result if one is given, within
   `Tolerances` (default 10⁻⁶ per quantity) **plus the estimated errors of the two results compared**
   (bounds, section 1.3).

The report lists the reasons, every pairwise `Comparison` (largest relative difference per quantity
with location and both values, and the ten largest overall), the global equilibrium of each result,
and the fidelity records. Differences are relative to model-wide scales of the same physical kind:
largest translation; largest rotation (at least translation scale / model size); largest force among
reactions and end forces; largest moment (at least 10⁻³ × force scale × model size).

## 2. The ODE oracle (`benchfea.oracle`)

`Beam(length, EI, GAs=None, N=0.0, second_order=False, shear_model=None, plane="xy")` with
`.support(x, v=, r=, settlement=, rotation=)`, `.point(x, P=, M=)`, `.dist(x1, x2, w1, w2)`,
`.hinge(x, k=0)` and `.solve(tol)` → `v(x)`, `theta(x)`, `M(x, side)`, `Q(x, side)`, `V(x, side)`,
`reactions`. `Bar(length, EA)` is the axial / torsional analogue. `EI`, `GAs`, `N` may be piecewise
constant.

Equations (derived in the module docstring from the equilibrium of a displaced slice, plane sections
and the shear angle γ = v′ − θ; signs of CONVENTIONS §4):

```
v' = θ + V/GAs        θ' = M/EI        M' = −Q + N v'        Q' = −w
joints:  Q⁺ = Q⁻ − P + k v      M⁺ = M⁻ − Mo + k_r θ      hinge:  M⁻ = M⁺ = k_h (θ⁺ − θ⁻)
```

`Q` is the force along the undeformed transverse axis. First order: `V = Q`, `N` ignored. Second
order with shear deformation and axial force needs the definition of the shear force that enters the
shear strain — `"engesser"` (`V = Q − N v′`) or `"haringx"` (`V = Q − N θ`). The mechanics does not
single one out, so the oracle **refuses** such a problem unless `shear_model` is given.

The member is cut at every station, load end and property change; the segments are mapped to [0, 1]
and solved as one two-point BVP with the joint relations as boundary conditions, in scaled variables.
`tol` bounds the collocation residual, not the error: treat results as accurate to about 100·tol
(10⁻⁸ relative by default). A mechanism or a load at/above a buckling load makes the problem singular
and raises `OracleError`.

## 3. Formulation of each solver and how the translators use it

### 3.1 Summary

| | benchfea (target) | OpenSeesPy 3.8 | PyNite 3.2.1 |
|---|---|---|---|
| element | prismatic 3D frame, exact for the member ODE | `elasticBeamColumn` (Euler–Bernoulli) | 12×12 Euler–Bernoulli |
| shear deformation | optional (`analysis.shearDeformation`) | `ElasticTimoshenkoBeam`, first order only (**measured**: ignores P-Delta, mishandles `-jntOffset`, accepts no release flags, no `-beamPoint`) | none (**doc**: "not currently considered") → `unsupported` |
| local axes | `R` rows = x, y, z | `vecxz` = benchfea z; y = vecxz × x, z = x × y (**doc**, **measured**) | own Y-up rule + `rotation`; the worker computes the rotation from PyNite's own triad and verifies every member to 10⁻⁹ |
| geometric stiffness | (engine) | `PDelta`: chord-rotation term N/L on the end translations only (**doc**: no P-δ; **measured**: one element gives H/(3EI/L³ − P/L)) → subdivision, O(1/n²) when N is the mean axial force of the element (section 3.4) | consistent cubic Kg + axial term (EA+N)/L + torsional term N(Iy+Iz)/A (**measured**, source); the two extra terms are removed by the worker |
| axial force in Kg | (engine) | current force of element end j, without the axial part of an element load (**measured**); Newton iterations to convergence | PyNite: first-order force, one solve; the worker iterates to convergence with PyNite's own matrices (section 3.3) |
| constraints | exact master/slave (linear kinematic relation, also in second order) | `rigidLink('beam')`, `rigidDiaphragm`, `equalDOF` with the `Transformation` handler: exact; a constrained node must not be the retained node of another constraint (**doc**; **measured**: wrong results) | none → stiff members, `approximate` in first order, `unsupported` in second order (section 3.5) |
| rigid offsets | exact | `-jntOffset` (global vectors node → element end) on `elasticBeamColumn`; explicit node + `rigidLink('beam')` otherwise | none → stiff stub member, `approximate` |
| releases | per local DOF, full or spring | flags `-releasey/-releasez` (bending, full) only | boolean, any of the 12 end DOFs |
| member loads | point force/moment, linear line loads | `-beamUniform` (full/partial), `-beamPoint`, local axes; no point moment; **measured**: end intensities of a varying load are ignored without message; release flags condense fixed-end forces of full-span uniform loads only | point forces and moments, linearly varying line loads; **measured**: a varying load on a member that PyNite splits at ≥ 2 intermediate nodes violates statics |
| reactions | on structure | `nodeReaction` is zero at a fixed node retaining a constraint (**measured**) → computed by the translator | `RxnFX…` |
| linear solver | — | `UmfPack` (least round-off with penalties, but reports no singular matrix — **measured**) | sparse |

### 3.2 OpenSees translator (`benchfea.crosscheck.opensees`)

* **Members**: one `geomTransf` per element (`Linear` / `PDelta`), `elasticBeamColumn(A, E, G, J, Iy,
  Iz)` or `ElasticTimoshenkoBeam(E, G, A, J, Iy, Iz, Avy, Avz)` (argument order measured; `Avy` pairs
  with bending about z). A shear area that the section does not have is replaced by 10¹⁰·A.
* **Truss-type members**: one `elasticBeamColumn` with bending and torsion stiffness × 10⁻¹⁵. Nodes
  joined only by truss members get their rotations fixed (benchfea reports them as 0). In a second-order
  combination a truss member that carries a member load with a component **along its axis** is a
  pin-ended beam-column in benchfea (THEORY 10.6: its axial force varies, it bends, and the axial load
  acts on the bowed axis): it is then subdivided like a frame member, with its own bending stiffness,
  bending hinges at its two ends (parallel elements, below) and no torsion stiffness in its last element.
  With a constant axial force the single element is exact (mean axial force × chord rotation is the
  complete geometric effect of a pin-ended member; **measured** on bank model S06 with the truss loads
  turned transverse: tableau = engine to 2·10⁻¹³). A truss section without bending stiffness that would
  need this treatment is `unsupported`.
* **Bending connections** (hinge or rotational spring about local y/z at either end): *exact
  decomposition into parallel elements*. With joint i held, a member with connection stiffnesses
  k_i, k_j has the 2×2 stiffness K = F⁻¹, F = F₀ + (1/k_i)[[L², L],[L, 1]] + (1/k_j)[[0,0],[0,1]]
  (cantilever flexibility plus the two connection rotations). K is symmetric, so it is a combination
  of three matrices: K = a·K(rigid ends, EI·a, GA·a) + b·K(hinge at j, EI′) + c·K(hinge at i, EI′),
  EI′ = EI/(1 + 3EI/(GA L²)) — three linear equations for (a, b, c); for one spring on an
  Euler–Bernoulli member a = 1/(1 + 4EI/(kL)), b = 1 − a. The coefficients are non-negative
  (checked on random data). The translator emits a core element (stiffnesses × a) plus up to two
  partners per bending plane that carry only that plane's bending (other stiffnesses × 10⁻¹⁵) with a
  `-releasey/-releasez` flag. This works for any member orientation, with offsets, with shear
  deformation and in second order (only the core carries axial force, so the geometric stiffness is
  not duplicated). Verified against the ODE oracle (`tests/test_crosscheck_fef.py`,
  `test_crosscheck_translators.py`).
* **Axial / shear / torsion connections**: duplicate node at the flexible end; `equalDOF` for rigid
  directions, `zeroLength` spring with `-orient` = member axes for partial ones. Exact when the member
  axes are parallel to the global axes (or the whole triple is rigid), the end has no offset and the
  joint is not a slave; otherwise penalty springs (default multiplier 10⁶, coarse run 10⁵) →
  `approximate`.
* **Constraints**: `rigidLink('beam', master, slave)`, `rigidDiaphragm(normal, master, *slaves)`,
  `constraints('Transformation')`. A diaphragm slave outside the master's plane is `unsupported`
  (benchfea's rule — in-plane part of the offset only — is not a rigid-body constraint; benchfea itself
  then fails its moment-equilibrium invariant). Memberless master DOFs that drive nothing are fixed.
* **Supports**: `fix`; springs = `zeroLength` to a fixed twin node; imposed displacements = `sp` in
  the load pattern.
* **Member loads**: `eleLoad` only where OpenSees is exact (uniform full/partial and interior point
  forces on unreleased Euler–Bernoulli elements; `native_loads=False` switches it off). Everything
  else — linearly varying loads, point moments, loads on Timoshenko or released elements — is applied
  as **exact equivalent nodal loads**: `benchfea.crosscheck.fef` computes the end forces of the member
  with held joints by the force method (cantilever statics, unit-load theorem with curvature M/EI,
  shear angle V/GA and the connection rotations M/k, exact Gauss quadrature of the piecewise
  polynomials), the opposite forces are applied to the element nodes and the fixed-end forces are added
  back to the element end forces. The two routes agree with each other and with the oracle to 10⁻⁸.
* **Second order**: `PDelta`, every frame member split into `subdivisions` elements (default 32;
  levels 8, 16, 32 for the error bound), Newton with `NormDispIncr` 10⁻¹⁰. The **axial component of
  every member load goes to the element nodes** as exact equivalent nodal loads (`fef`), the transverse
  components stay native: with an axial `eleLoad` the transformation uses the wrong axial force
  (section 3.4). Second order with shear deformation is `unsupported` (the member is then modelled
  without shear deformation and the result must not be used).
* **Self-checks on every run** (a failure makes the result `ok = False`): finite and sane
  displacements; violation of every kinematic constraint recomputed from the displacements < 10⁻⁸;
  joint equilibrium residual at all free DOFs < 10⁻⁶ (reactions come from the same balance, transferred
  through the known constraints).

### 3.3 PyNite translator (`benchfea.crosscheck.pynite`)

* **Axes**: PyNite's rule is Y-up with special cases for vertical and horizontal members. No
  coordinate permutation is used: the worker reads `member.T()` of the unrotated member, sets the
  `rotation` that brings its y axis onto benchfea's and re-reads the triad; the run fails unless all
  triads match to 10⁻⁹. Proven on 16 directions / rolls / reference vectors with Z-up and Y-up models
  (`test_member_axes_reproduced_exactly`).
* **Releases**: boolean releases for full releases; rotational springs about local y/z through the
  same parallel-member decomposition as OpenSees, with the member loads as exact equivalent nodal
  loads; partial axial / shear / torsion connections are `unsupported`.
* **Rigid offsets, rigid links, diaphragms**: stiff members with section properties `stiffness` ×
  the stiffest member (scaled with the stub length so that translational and rotational flexibility
  stay below 1/`stiffness` of the reference). Default 10⁶, coarse run 10⁵. **Measured**: the emulation
  error falls like 1/multiplier while PyNite's round-off grows with it; short offset stubs are best
  near 10⁵, long diaphragm links still improve at 10⁷. A diaphragm link has local y = normal and
  releases out-of-plane shear, torsion and out-of-plane bending at the slave; slaves must be in the
  master's plane and not at the master's position.
* **Topology**: PyNite connects any node lying on a member's axis to that member; benchfea does not.
  Such a model is `unsupported` for PyNite.
* **Loads**: PyNite's own member loads in local axes on explicit sub-members (loads clipped by the
  translator, never PyNite's own splitting).
* **Second order**: explicit sub-members (default 8; levels 2, 4, 8 for the error bound). PyNite's
  `analyze_PDelta` is one solve with the geometric stiffness of the first-order axial forces, and its
  member matrix `kg(P)` has two terms that benchfea's theory does not have — an axial term that turns
  EA/L into (EA + P)/L and a torsional term that turns GJ/L into (GJ + P(Iy + Iz)/A)/L (**measured**:
  `test_pynite_geometric_stiffness_has_an_axial_and_a_torsional_term`) — with P taken as strain × EA.
  None of this is estimated; the worker removes it, using PyNite's own member matrices, partitioning
  and solver:
  * axial force of a member N = strain × E × (the given area);
  * section properties of every member (each has its own section): A′ = A − N/E, so that EA′ + N = EA,
    and J′ = J − N(Iy + Iz)/(A′G), so that GJ′ + N(Iy + Iz)/A′ = GJ;
  * K = Ke(A′, J′) + Σ PyNite member matrices `Kg(N)`, assembled by node number as PyNite does; solve;
    repeat until the displacements change by less than 10⁻¹⁰ (typically 3 steps).

  At convergence PyNite's tangent stiffness is the cubic element with the consistent geometric
  stiffness of the current axial forces — a discretisation of the same beam-column problem, which
  converges to the references (table of section 1.3: true errors 10⁻⁶…10⁻⁵ on frames, against
  10⁻⁴…10⁻³ from N/EA alone before). Besides the discretisation bound three small estimates are
  recorded from the run: the last change of the iteration; the residual of the two removed terms with
  the converged strains; and the effect of PyNite's own force recovery, which still takes
  P = strain × EA′ = N(1 − strain) in the geometric stiffness — `kg(N − P) d` of every member and its
  sum per node, relative to the force and moment scales (10⁻⁷…10⁻⁶ typically).
  The earlier estimates (3 × a global compliance measure for the missing iteration, max |N|/EA for the
  axial term) were not bounds: on a generated braced frame the true error was 6.6·10⁻⁴ against an
  estimate of 5.8·10⁻⁴, and with the axial term removed alone the missing iteration left errors up to
  3× the estimate (bank S03, S04). A truss member with an axial member load is subdivided as in the
  OpenSees route.
* **Rigid links, diaphragms, offsets in second order**: `unsupported` / no estimate (section 3.5).
* **Self-checks**: axes (above); finite and sane displacements.

### 3.4 What differs in second order, and how to read differences

| effect | OpenSees route | PyNite route |
|---|---|---|
| P-Δ (chord rotation) | exact per element | included |
| P-δ within a member | only through subdivision (error ∝ 1/n²; 32 elements: about 2·10⁻⁴ at half the buckling load) | consistent Kg: nodal displacements converge with the fourth power of the sub-member length |
| axial forces | iterated (Newton) | iterated by the worker with PyNite's matrices |
| axial stiffness | EA/L | EA/L after the removal of PyNite's axial term |
| released members | flag in the basic system | Kg condensed separately from Ke (vanishes with subdivision; observed order then about 1) |
| end shears | normal to the displaced chord (`localForce`); the translator reports the global force rotated to the undeformed axes instead | undeformed axes |
| reactions | equilibrated | PyNite's own recovery, with the axial force N(1 − strain) in the geometric stiffness (estimated) |

**Why the OpenSees discretisation error was not second order.** The chord-rotation geometric stiffness
of the `PDelta` transformation is N/L times the relative end translation. With N the *mean* axial force
of the element this is a second-order accurate discretisation (measured orders 1.97…2.02 in the table
of section 1.3). OpenSees takes for N the axial basic force of the element, which is the force of
**element end j without the axial part of an element load**: **measured**
(`test_opensees_pdelta_takes_the_axial_force_of_end_j_with_an_axial_element_load`) — a cantilever of one
element whose only axial load is an `eleLoad -beamUniform` with an axial component has the lateral
stiffness 3EI/L³, i.e. N = 0, although its mean axial force is wL/2. Per element the force is off by
half the axial load on that element, w h/2: an error of first order in the element length, in the
geometric stiffness of every element. Models with self weight on columns or with gravity line loads on
inclined members therefore converged like a/n² + b/n — measured on bank model S07 (columns under self
weight): error ratios 1.8…2.0 per halving, estimated error 2.3·10⁻² at 16 elements — and the Richardson
estimate with the order 2 understated the error (portal_frame 3/1: differences to the engine 4.2·10⁻⁴,
1.4·10⁻⁴, 5.2·10⁻⁵, 2.1·10⁻⁵, 9.5·10⁻⁶ at 8…128 elements). The translator now applies the axial component
of member loads at the element nodes (exact equivalent nodal loads); the axial force of an element is
then constant and equals the mean of the exact one for any axial load (lever rule), and the order is 2
on S07 as well. The bound of section 1.3 does not rely on this: it assumes first order.

### 3.5 Rigid links and diaphragms in second order

benchfea's constraint is the linear relation u_s = u_m + θ_m × r (THEORY 4); the second-order effect of
the forces the constraint transmits is not included, which THEORY 16 lists as a limit and the results
carry as a warning. OpenSees' `rigidLink` / `rigidDiaphragm` with the `Transformation` handler is the
same linear relation: same problem, fidelity `exact`.

The mechanics of a rigid body, to the order of the P-Δ terms: a body that turns by θ about the
diaphragm normal moves a slave at r by θ × r − ½θ² r. The in-plane force f_s that the body applies to
slave s does the second-order work −½θ² (f_s · r_s), so the body adds the rotational stiffness

```
K_θθ = − Σ_s f_s · r_s
```

at the master (for a rigid link in general ½(r fᵀ + f rᵀ) − (r·f) I on the three rotations; the same
derivation as the rigid end offset of THEORY 10.7, which the engine does include). It is of the same
order in the displacements as the N/L terms of the columns; it is positive (stabilising) when the links
are in tension.

Measured on bank model S08 (five columns under a rigid floor, twist 0.046 rad, displacements 3.25× the
first-order ones), largest displacement u_x of the master:

| statement of the problem | u_x [mm] |
|---|---|
| linear constraint: benchfea | 193.83 |
| linear constraint: OpenSees, 64 elements per member | 193.83 |
| linear constraint + K_θθ = −Σ f_s·r_s as a rotational spring at the master, iterated to its fixed point (278 kNm/rad), solved by benchfea | 179.11 |
| stiff members instead of the constraint: PyNite route (axial forces of the links iterated) | 177.67 |

So the stiff-member emulation reproduces the rigid-body statement (to 0.8 %; the remainder is PyNite's
separately condensed geometric stiffness of the released link directions), not the linear constraint,
and the two statements differ by 7.6 % on this model. Consequences:

* The PyNite route marks rigid links and diaphragms in a second-order combination **`unsupported`**
  (`constraints.<type>.secondOrder`): it solves a different problem, not an approximation of the same
  one, and no estimate of the difference is available. Such a combination has one eligible foreign
  solver and is **not cross-checked** (category `not-representable`). Before, the record said
  `approximate` with the emulation-stiffness estimate, and the model failed as a disagreement of 3.5 %
  (PyNite's one-pass result with first-order link forces).
* Which statement is right for a rigid diaphragm in second order: the mechanics says the rigid body
  carries the term; the linear constraint is the engine's documented simplification. That benchfea and
  OpenSees agree to 10⁻⁴ shows only that both implement the linear relation. The size of what is left
  out can be read from the engine's own results as |θ Σ f_s·r_s| against the applied torque (S08:
  13.6 kNm against 30 kNm); it is small when the floor forces are small against the torsional
  stiffness of the storey, and it is not small for leaning-column systems under torsion. Whether the
  engine should include the term (as it does for rigid offsets) is an open question for the engine
  owner (section 6); AISC 360-22 C2.1 does not address constraints.

## 4. Capabilities (`benchfea.crosscheck.CAPABILITIES`)

Machine-readable: `CAPABILITIES[solver][feature] = {"support", "how", "evidence"}` for the features in
`capabilities.FEATURES`; `FORMULATION[solver]` and `SOURCES[solver]` carry the formulation summary and
the documentation URLs. Summary:

| feature | OpenSees | PyNite |
|---|---|---|
| shear deformation | exact (first order) | unsupported |
| truss members, support springs, imposed displacements, nodal loads | exact | exact |
| full bending releases, rotational springs | exact | exact |
| axial / shear / torsion releases (full) | exact if member axes global, no offset, joint not a slave; else approximate | exact |
| axial / shear / torsion springs | as above | unsupported |
| rigid offsets | exact (first order) | approximate |
| rigid link, rigid diaphragm (coplanar slaves) | exact (linear relation, also in second order) | approximate; second order: unsupported (3.5) |
| point forces / moments, uniform, partial, linearly varying loads, self weight | exact | exact |
| second order | approximate (bound from three levels) | approximate (bound from three levels; axial forces iterated, extra Kg terms removed) |
| truss member with an axial member load, second order | approximate (subdivided pin-ended beam-column) | same |
| second order with shear deformation | unsupported | unsupported |
| second order with rigid offsets | approximate without estimate (not eligible) | same |
| linear buckling | not implemented (feasible from two tangent matrices) | unsupported |
| tension-/compression-only members | not implemented (tier C supplies the active-set model, with the loads of slack members as nodal loads) | same |
| "computed" τ_b, notional loads | values must be supplied by the caller | same |

## 5. Limits

* A model with shear deformation has only one eligible foreign solver (OpenSees, first order), so it
  cannot be cross-checked under the two-solver rule; second order with shear deformation has none.
* Second order: both routes are discretisations of the beam-column problem; the allowed difference to
  the engine is 10⁻⁶ plus the bound of the solver compared — typically 10⁻⁵…10⁻³, three times the true
  OpenSees error. A model whose bound exceeds 5·10⁻³ for one solver (members near their buckling load)
  has only one eligible solver. Rigid offsets in second order are not eligible; rigid links and
  diaphragms in second order have one eligible solver (3.5).
* The bound assumes that the level-to-level differences keep contracting at least like a first-order
  method beyond the finest level run. Three levels cannot prove that; the acceptance table is the
  evidence that it holds on the kinds of model listed there.
* The PyNite route drives PyNite's internals (`Analysis._partition`, `_store_displacements`,
  `_calc_reactions`, `Member3D.Kg`) of version 3.2.1; the measured tests fail if an upgrade changes them.
* PyNite emulations (offsets, links, diaphragms) reach about 10⁻⁵…10⁻⁴; their estimate widens the
  tolerance accordingly.
* Tension-/compression-only members, linear buckling, notional loads and a "computed" τ_b are not
  reproduced independently (see the table); warping torsion, tapered members and distributed moments
  do not exist in the model.
* The fixed-end forces and the parallel-member decomposition (`fef`) are shared by both translators for
  released / spring-connected / specially loaded members; they are independent of `benchfea.solve` and
  verified against the ODE oracle, but an error in them would be common to both foreign routes.
* The oracle covers one bending plane (plus the axial/torsional bar) of a straight member with
  piecewise constant properties; no biaxial coupling, no frames.
* Large models: the PyNite route is slow (dense per-member Python loops); OpenSees uses a sparse solver.

## 6. Open questions

1. Second-order member end shears: the translators report components along the undeformed local
   axes, as the engine does (THEORY 10.8).
2. Second order with shear deformation: Engesser or Haringx shear force? The oracle offers both and
   refuses to choose; no foreign solver available here covers the case.
3. Notional loads: sign for net uplift and the lumping of member loads to nodes are engine decisions;
   the translators take the nodal loads from the caller (`extra_nodal_loads`).
4. Geometric stiffness of rigid offset zones in second order: undefined in the contract; not eligible.
   Rigid links and diaphragms in second order: the engine's linear constraint leaves out the term
   −Σ f_s·r_s of section 3.5 (7.6 % on bank model S08). Should the engine include it, as it does for
   offsets? Until then such a combination is not cross-checked.
5. Rigid diaphragm with slaves off the master's plane: benchfea's rule fails its own moment-equilibrium
   invariant; unsupported in both translators.

## 7. Documentation pages used

OpenSees — openseespydoc.readthedocs.io/en/latest/src/: `elasticBeamColumn.html`,
`ElasticTimoshenkoBeam.html`, `LinearTransf.html`, `pdelta.html`, `eleload.html`, `rigidLink.html`,
`rigidDiaphragm.html`, `equalDOF.html`, `TransformationMethod.html`, `ZeroLength.html`,
`eleForce.html`; opensees.berkeley.edu/wiki/index.php/: `Linear_Transformation`,
`PDelta_Transformation`, `RigidDiaphragm_command`;
opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/elasticBeamColumn.html.
PyNite — pynite.readthedocs.io/en/latest/: `member.html`, `PDelta.html`, `analysis.html`,
`node.html`, `stability.html`.
Frame3DD (not run) — frame3dd.sourceforge.net, Frame3DD-manual.html in the project's svn trunk.
Full URLs: `benchfea.crosscheck.capabilities.SOURCES`.
