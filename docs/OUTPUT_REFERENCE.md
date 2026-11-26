# Output reference

The two files the engine writes, key by key. Both are plain JSON (or YAML, if the output file name
ends in `.yaml`), and **every number is in the units of the input file**. Axes and signs are those
of [CONVENTIONS.md](CONVENTIONS.md); they are not repeated here.

* `benchfea run FILE -o OUT` writes the **results document** (format `benchfea-results`, 1.0.0).
* `benchfea design FILE -o OUT` writes the **design document** (format `benchfea-design`, 1.0.0);
  with `--results FILE2` it also writes the results document it used.

This file is written by hand. `tests/test_output_reference.py` runs the five worked examples and
fails when an output document holds a key that is not named below; it does not check the
descriptions. The format versions and the rule for changing them are in
[FORMAT_VERSIONING.md](FORMAT_VERSIONING.md).

A path such as `combinations/<combination>/reactions/<node>/Fz` means: the key `combinations`,
then the id of a combination as written in the input file, then `reactions`, and so on.

When a model is refused (invalid input, a mechanism, a load beyond the elastic critical load, a
failed self-check without `--allow-flagged`), **no document is written**; the reason is printed
with a problem code, the place in the input file and a hint.

## 1. Results document

### 1.1 Top level

| Key | Content |
|---|---|
| `format`, `formatVersion` | `"benchfea-results"`, `"1.0.0"` |
| `engineVersion` | version of the `benchfea` package that wrote the file |
| `model` | the `name` of the input file |
| `units` | `{force, length}` as declared in the input file |
| `unitLabels` | the unit of a `length`, `rotation`, `force` and `moment` in this file, as text (for example `kN*mm`); rotations are always radians |
| `conventions` | `{verticalAxis}` of the input file |
| `combinations` | one entry per combination of the input file (1.2) |
| `buckling` | one entry per combination for which a buckling analysis was requested (1.3); empty otherwise |
| `envelopes` | per analysis member, the largest and smallest station forces over all combinations (1.4) |
| `nodalEnvelopes` | the largest and smallest node displacements and reactions over all combinations (1.4) |
| `defaultsUsed` | list of `{path, value, reason}`: every value the engine filled in because the input file did not give it |
| `warnings` | list of texts that concern the whole model |

### 1.2 `combinations/<combination>`

| Key | Content |
|---|---|
| `analysisType` | `"firstOrder"` or `"secondOrder"` |
| `converged`, `iterations` | whether the iteration of a second-order or one-sided-member analysis converged, and the number of iterations (first order without one-sided members: `true`, 1) |
| `alpha` | the factor of AISC 360-22 C2.1(d) used: loads were multiplied by it for the analysis and the results divided by it. The numbers in the file are already divided |
| `tauB` | `{<member>: value}`: the stiffness reduction factor used for each member; empty when no stiffness reduction was asked for |
| `displacements/<node>` | `ux, uy, uz` (length) and `rx, ry, rz` (radians) of every node, global axes |
| `reactions/<node>` | `Fx, Fy, Fz` (force) and `Mx, My, Mz` (moment) at every supported node, global axes: the force of the support **on the structure**. Spring supports are included |
| `constraintForces/<node>` | present when the model has rigid links or diaphragms: the forces the constraint exerts on its master and slave nodes, global axes, same six names |
| `members/<member>/active` | `false` for a tension-only or compression-only member that is slack in this combination |
| `members/<member>/endForces/i`, `…/j` | `Fx, Fy, Fz, Mx, My, Mz`: forces of the joints on the member at the two ends of its flexible length, member local axes |
| `members/<member>/stations` | lists of equal length, one value per result station: `x` (distance from the start of the flexible length), the internal forces `N, Vy, Vz, T, My, Mz`, and the displacements of the member axis `ux, uy, uz` in member local axes (total displacements, rigid-body motion included). At a point load a position occurs twice, with the values on either side |
| `invariants/passed` | `true` when every self-check of this solution is within its tolerance |
| `invariants/checks` | list of `{name, value, tolerance, passed, detail}`: each self-check (equilibrium of the whole structure, of the joints and of the members, symmetry, energy balance, and others), its measured relative error and its tolerance |
| `solver` | present for an iterated analysis: the record of the iteration. `secondOrder` holds `iterations`, `converged`, `tolerance`, `maxIterations`, `history` (the changes per iteration), `subdivisions` (pieces each member was cut into), `largestPieceKh` and `subdivisionCheck`; `activeSet` holds `iterations`, `slack` (the one-sided members found slack), `history` and `strategy` of the search; `stiffnessReduction` and `notionalLoads` record what was applied |
| `warnings` | list of texts that concern this combination |

### 1.3 `buckling/<combination>`

| Key | Content |
|---|---|
| `factors` | buckling factors in ascending order: multiples of the loads of the combination at which the elastic structure buckles |
| `modes` | one entry per factor: `{<node>: {ux, uy, uz, rx, ry, rz}}`. A mode is a shape without a scale; it is normalised so that the largest translation, over the nodes and along the members, is 1. Rotations are per unit of the file's length |
| `memberShapes/<member>` | one entry per factor: `{ux, uy, uz}`, lists over the result stations, member local axes, same normalisation |
| `negativeFactors`, `negativeModes` | present when the structure can also buckle under the **reversed** loads: the factors (negative, closest to zero first) and their shapes |
| `info` | how the factors were computed: `basis`, `piecesPerMember` and `subdivisions` (pieces per member), `freeDofs`, `estimatedRelativeError` and `estimatedRelativeErrorNegative` (an estimate from two subdivisions, not a bound), `tolerance`, `exactMemberStiffness`, `slack` (one-sided members taken as slack), `activeSet`, `tauB`, `stiffnessReduction`, `warnings` |

No buckling factor can reach the status "verified" ([LIMITS.md, section 4.3](../LIMITS.md#43-what-cannot-currently-be-verified)).

### 1.4 Envelopes

| Key | Content |
|---|---|
| `envelopes/<member>/x` | station positions, as in 1.2 |
| `envelopes/<member>/max`, `…/min` | `{N, Vy, Vz, T, My, Mz}`: lists over the stations with the largest and the smallest (algebraic) value over all combinations |
| `envelopes/<member>/maxCombination`, `…/minCombination` | same layout, holding the id of the combination that gives each value |
| `nodalEnvelopes/displacements/<node>`, `nodalEnvelopes/reactions/<node>` | `max`, `min` (six components each) and `maxCombination`, `minCombination` (the governing combination of each component; the first one when several tie). Reactions are listed for supported nodes only |

The envelopes run over **all** combinations of the file, strength and serviceability alike, and
over first-order and second-order combinations alike. There is no envelope of station
displacements.

## 2. Design document

### 2.1 Top level

| Key | Content |
|---|---|
| `format`, `formatVersion` | `"benchfea-design"`, `"1.0.0"` |
| `code` | the design code used, for example `"AISC 360-22"` |
| `model`, `units` | as in the results document |
| `status` | the worst member status: `pass`, `conditional`, `incomplete` or `fail` ([LIMITS.md, section 3.1](../LIMITS.md#31-the-rule-you-must-know)) |
| `combinations/<combination>` | the analysis basis of each combination as the design module sees it: `purpose` (`strength` or `serviceability`), `inResults`, `analysisType`, `secondOrder`, `converged`, `stiffnessReduction`, `notionalLoad`, `alpha`, `directAnalysisSettings` |
| `members/<design member>` | the checks of each design member (2.2) |
| `defaultsUsed` | as in the results document |
| `warnings` | list of texts that concern the whole design run |

### 2.2 `members/<design member>`

| Key | Content |
|---|---|
| `status` | `pass`, `conditional`, `incomplete` or `fail` |
| `passed` | `true` only when no checked limit state is exceeded and no applicable limit state was left unchecked (status `pass` or `conditional`) |
| `governing` | the check with the largest ratio among those that were made (one check, layout of 2.3); `null` if none was made |
| `checks` | list of checks: for each limit state, the combination that governs it |
| `byCombination/<combination>` | list of every check made for that combination |
| `notChecked` | list of texts: limit states that apply but were not checked, each with the reason |
| `parameters/<design variable>` | `{value, unit, source, path}`: every design variable as resolved for this member. `source` says where the value came from (`member`, `group:<id>`, `global`, `default`, `material:<id>`, `bracing`, or a derivation), `path` is the place in the input file |
| `warnings` | list of texts that concern this member |

### 2.3 One check

| Key | Content |
|---|---|
| `name` | the limit state, in words |
| `clause` | the citation: design code, section and equation numbers |
| `status` | `ok`; or `notChecked` / `refused`, in which case `reason` says why and `capacity`, `ratio`, `passed` are `null` |
| `reason` | text; empty when `status` is `ok` |
| `method` | `LRFD` or `ASD` |
| `demand` | required strength, or the measured quantity of a serviceability check |
| `nominal`, `factor`, `factorName` | nominal strength, and the resistance factor (`phi`) or safety factor (`Omega`) applied to it |
| `capacity` | available strength, or the limit of a serviceability check |
| `unit` | unit of `demand`, `nominal` and `capacity`, as text |
| `ratio` | `demand` / `capacity` |
| `passed` | `true`, `false`, or `null` when not checked |
| `advisory` | `true` for a check that never changes the member status (the slenderness notes) |
| `combination`, `x` | the combination, and the distance from the start of the design member, at which the ratio occurs |
| `inputs/<name>` | `{value, unit, source}`: every input the check used, with where it came from |
| `details/<name>` | `{value, unit}`: intermediate values of the calculation |
| `notes` | list of texts: assumptions, engine rules that are not text of the design code, and remarks |

Nothing in either document says that a result has been reviewed, approved or certified by anyone.
What "verified" means, and which results can reach that status, is in
[VERIFICATION.md](VERIFICATION.md).
