# Example 05 - braced three-storey frame with tension-only X-bracing

**File:** `05_braced_frame.json` (kN, m; vertical axis Z).

## What it is

Two bays of 6 m and three storeys (4.0 m, 3.5 m, 3.5 m) in the X-Z plane, all nodes held out of
plane. Columns are continuous and pinned at the base; beams are pinned to the columns
(`releases` of `rz` at both ends); bay A-B has X-bracing in every storey.

```
 A3 +-----+-----+ C3
    | \ / |     |
 A2 +--X--+-----+ C2       X = two diagonals, type "tensionOnly"
    | / \ |     |
 A1 +-----+-----+ C1
    | \ / |     |
 A0 ^--X--^-----^ C0
```

* **Tension-only members.** The diagonals have `"type": "tensionOnly"`: in each combination the engine
  finds which of them are stretched and switches the others off (a slack diagonal takes no force
  from the frame; it only carries its own weight to its two ends). Because
  that depends on the loads, every combination is solved on its own - results of load cases cannot be
  superposed.
* **Every combination has some lateral load.** Under gravity alone both diagonals of a panel shorten
  with the columns and go slack, and a frame with pinned beams and pinned bases would then be a
  mechanism (the engine refuses it and says why). The gravity combinations therefore include the load
  case `N`, a small lateral load of 0.2 % of the gravity load of each level.
* **Combinations** cover gravity, wind from both sides (`W` with factor +1 and -1 - the other set of
  diagonals works) and a minimum-gravity case.
* **Linear buckling** of `ULS1` (`analysis.buckling`), on the set of diagonals that work in `ULS1`.
* **Design groups** give all columns, floor beams, roof beams and braces their design variables in
  one place (`role`, K factors, flange bracing, deflection limits).

## Run it

```
benchfea validate examples/05_braced_frame.json
benchfea run      examples/05_braced_frame.json -o 05_braced_frame.results.json
benchfea design   examples/05_braced_frame.json -o 05_braced_frame.design.json
```

Numbers below are in the file's units (kN, m); signs follow docs/CONVENTIONS.md. They are what
the engine computes for this file (`tests/test_examples.py` recomputes every row); they are not
reference values from any publication.

## Key results

| quantity | value | unit |
|---|---|---|
| ULS2 (wind from the left): tension in brace X1a | 77.46 | kN |
| ULS3 (wind from the right): tension in brace X1b | 76.49 | kN |
| ULS2: number of slack braces | 3 | - |
| ULS1 (gravity + minimum lateral load): tension in brace X1a | 6.087 | kN |
| ULS1: axial force in column CB1 (compression negative) | -1013 | kN |
| SLS2: roof drift (A3 along X) | 0.008753 | m |
| ULS4: vertical reaction at A0 (windward column, least compression) | 153.5 | kN |
| buckling factor of ULS1, mode 1 | 7.653 | - |
| design: governing ratio, column CB1 | 0.8184 | - |
| design: governing ratio, floor beam BA1 | 0.9425 | - |
| design: governing ratio, brace X1a | 0.434 | - |

Under wind from the left the three `a` diagonals work and the three `b` diagonals are slack; under
wind from the right it is the other way round.

## Note on verification

No foreign solver of the cross-check (docs/CROSSCHECK.md) reproduces tension-only members or linear
buckling, so `benchfea verify` quarantines this file at tier C: its results are what the engine
computes, without an independent cross-check. Examples 01, 02 and 04 pass that tier.

## Design summary

| design member | status | ratio | combination | governing limit state |
|---|---|---|---|---|
| CA1 | conditional | 0.408 | ULS1 | compression: flexural buckling about the minor axis (y) |
| CA2 | conditional | 0.223 | ULS1 | combined flexure and compression |
| CA3 | conditional | 0.079 | ULS1 | compression: flexural buckling about the minor axis (y) |
| CB1 | conditional | 0.818 | ULS1 | compression: flexural buckling about the minor axis (y) |
| CB2 | conditional | 0.445 | ULS1 | combined flexure and compression |
| CB3 | conditional | 0.157 | ULS1 | compression: flexural buckling about the minor axis (y) |
| CC1 | conditional | 0.409 | ULS1 | compression: flexural buckling about the minor axis (y) |
| CC2 | conditional | 0.223 | ULS1 | combined flexure and compression |
| CC3 | conditional | 0.079 | ULS1 | compression: flexural buckling about the minor axis (y) |
| BA1 | conditional | 0.943 | SLS1 | deflection limit 1 (local y) |
| BB1 | conditional | 0.943 | SLS1 | deflection limit 1 (local y) |
| BA2 | conditional | 0.943 | SLS1 | deflection limit 1 (local y) |
| BB2 | conditional | 0.943 | SLS1 | deflection limit 1 (local y) |
| BA3 | conditional | 0.428 | SLS1 | deflection limit 1 (local y) |
| BB3 | conditional | 0.428 | SLS1 | deflection limit 1 (local y) |
| X1a | conditional | 0.434 | ULS2 | combined flexure and tension |
| X1b | conditional | 0.430 | ULS3 | combined flexure and tension |
| X2a | conditional | 0.290 | ULS2 | combined flexure and tension |
| X2b | conditional | 0.294 | ULS3 | combined flexure and tension |
| X3a | conditional | 0.131 | ULS2 | combined flexure and tension |
| X3b | conditional | 0.131 | ULS3 | combined flexure and tension |

Status `conditional` means: every check passes, but a strength combination was analysed first
order, and AISC 360-22 C1 asks for second-order effects in the required strengths (the engine says
so instead of amplifying silently).
