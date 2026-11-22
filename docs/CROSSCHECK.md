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

### Subprocess

`run` writes the job (plain JSON: the list of library calls and the outputs to read) to a temporary
file and executes `python -m benchfea.crosscheck.worker job.json out.json` with a timeout
(default 300 s). A crash, a timeout, a missing package or a failed analysis come back as
`ForeignResult(ok=False, error=...)`; nothing raises into the caller. Reason: OpenSees terminates the
interpreter on some input errors (**measured**: a `vecxz` parallel to the member axis; a model with a
zero-length element). `in_process=True` runs the solver in the calling process (used for PyNite in the
tests, which is pure Python).

### Fidelity record

Every translation returns, per model feature, one of

| level | meaning |
|---|---|
| `exact` | the foreign model is the same mathematical problem (differences are round-off) |
| `approximate` | the foreign model converges to the same problem; `error` = estimated relative error |
| `unsupported` | the solver cannot represent the feature or the translator refuses to guess; any numbers returned must not be used as a check |

The estimate of an approximate feature comes from a second run at a coarser value of its refinement
parameter (`info["convergence"]`):

| control | used for | estimate |
|---|---|---|
| `subdivisions` (elements per member) | second order | `\|fine − coarse\| / (ratio² − 1)` (error ∝ 1/n², measured for OpenSees; PyNite converges faster, so this is conservative there) |
| `penalty` (OpenSees), `stiffness` (PyNite) | penalty springs, stiff emulation members | `2 × \|fine − coarse\|` with ratio 10 — no error law is assumed because round-off grows with the stiffness (measured) |

### Tier-C decision (`cross_check`)

A combination is **cross-checked** only if

1. at least two foreign solvers ran, have no `unsupported` feature, and every `approximate` feature
   has an estimate whose sum is ≤ `max_approximation` (default 5·10⁻³) — these are the *eligible* solvers;
2. all eligible solvers agree pairwise, and with the benchfea result if one is given, within
   `Tolerances` (default 10⁻⁶ per quantity) **plus the estimated errors of the two results compared**.

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
| geometric stiffness | (engine) | `PDelta`: chord-rotation term N/L on the end translations only (**doc**: no P-δ; **measured**: one element gives H/(3EI/L³ − P/L)) → subdivision, O(1/n²) | consistent cubic Kg + axial term (EA+N)/L + torsional term; **one** solve with first-order N |
| axial force in Kg | (engine) | current force, Newton iterations to convergence | first-order force, not iterated; reactions/end forces post-processed with N from final displacements (not exactly equilibrated) |
| constraints | exact master/slave | `rigidLink('beam')`, `rigidDiaphragm`, `equalDOF` with the `Transformation` handler: exact; a constrained node must not be the retained node of another constraint (**doc**; **measured**: wrong results) | none → stiff members, `approximate` |
| rigid offsets | exact | `-jntOffset` (global vectors node → element end) on `elasticBeamColumn`; explicit node + `rigidLink('beam')` otherwise | none → stiff stub member, `approximate` |
| releases | per local DOF, full or spring | flags `-releasey/-releasez` (bending, full) only | boolean, any of the 12 end DOFs |
| member loads | point force/moment, linear line loads | `-beamUniform` (full/partial), `-beamPoint`, local axes; no point moment; **measured**: end intensities of a varying load are ignored without message; release flags condense fixed-end forces of full-span uniform loads only | point forces and moments, linearly varying line loads; **measured**: a varying load on a member that PyNite splits at ≥ 2 intermediate nodes violates statics |
| reactions | on structure | `nodeReaction` is zero at a fixed node retaining a constraint (**measured**) → computed by the translator | `RxnFX…` |
| linear solver | — | `UmfPack` (least round-off with penalties, but reports no singular matrix — **measured**) | sparse |

### 3.2 OpenSees translator (`benchfea.crosscheck.opensees`)

* **Members**: one `geomTransf` per element (`Linear` / `PDelta`), `elasticBeamColumn(A, E, G, J, Iy,
  Iz)` or `ElasticTimoshenkoBeam(E, G, A, J, Iy, Iz, Avy, Avz)` (argument order measured; `Avy` pairs
  with bending about z). A shear area that the section does not have is replaced by 10¹⁰·A.
* **Truss-type members**: one `elasticBeamColumn` with bending and torsion stiffness × 10⁻¹⁵, never
  subdivided. Nodes joined only by truss members get their rotations fixed (benchfea reports them as 0).
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
* **Second order**: `PDelta`, every frame member split into `subdivisions` elements (default 16,
  coarse 8), Newton with `NormDispIncr` 10⁻¹⁰. Second order with shear deformation is `unsupported`
  (the member is then modelled without shear deformation and the result must not be used).
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
* **Second order**: `analyze_PDelta` on explicit sub-members (default 4, coarse 2). Besides the
  discretisation estimate, two estimates come from the run itself:
  *axial-force iteration* — with K_T u = F, a change δKg changes the compliance by Fᵀδu = −uᵀ δKg u;
  δKg is linear in ΔN = N(final displacements) − N(first order) and uᵀ kg u of a piece per unit axial
  force is (chord rotation)² × length, so `3 × |Σ ΔN (chord rotation)² L| / |Fᵀu|` is recorded (a
  first-order estimate of a global measure, not a bound for every result);
  *axial term in Kg* — `max |N| / EA`.
* **Self-checks**: axes (above); finite and sane displacements.

### 3.4 What differs in second order, and how to read differences

| effect | OpenSees route | PyNite route |
|---|---|---|
| P-Δ (chord rotation) | exact per element | included |
| P-δ within a member | only through subdivision (error ∝ 1/n²; 16 elements: about 10⁻³ at half the buckling load) | consistent Kg: one element already within about 0.5 %, converges ∝ 1/n⁴ |
| axial forces | iterated | first order |
| axial stiffness | EA/L | (EA + N)/L |
| released members | flag in the basic system | Kg condensed separately from Ke (vanishes with subdivision) |
| end shears | normal to the displaced chord (`localForce`); the translator reports the global force rotated to the undeformed axes instead | undeformed axes |
| reactions | equilibrated | not exactly equilibrated (N recomputed after the solve) |

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
| rigid link, rigid diaphragm (coplanar slaves) | exact | approximate |
| point forces / moments, uniform, partial, linearly varying loads, self weight | exact | exact |
| second order | approximate (converged) | approximate (converged + formulation estimates) |
| second order with shear deformation | unsupported | unsupported |
| second order with rigid offsets | approximate without estimate (not eligible) | same |
| linear buckling | not implemented (feasible from two tangent matrices) | unsupported |
| tension-/compression-only members | not implemented | not implemented |
| "computed" τ_b, notional loads | values must be supplied by the caller | same |

## 5. Limits

* A model with shear deformation has only one eligible foreign solver (OpenSees, first order), so it
  cannot be cross-checked under the two-solver rule; second order with shear deformation has none.
* Second order: both routes are approximations of the beam-column problem; the widened tolerance is
  typically 10⁻⁴…10⁻³. Rigid offsets in second order are not eligible.
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
   axes. `benchfea.solve` has no second-order analysis yet; if it reports chord-normal shears the
   three-way test `test_second_order_three_way_when_the_engine_supports_it` will show it.
2. Second order with shear deformation: Engesser or Haringx shear force? The oracle offers both and
   refuses to choose; no foreign solver available here covers the case.
3. Notional loads: sign for net uplift and the lumping of member loads to nodes are engine decisions;
   the translators take the nodal loads from the caller (`extra_nodal_loads`).
4. Geometric stiffness of rigid offset zones in second order: undefined in the contract; not eligible.
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
