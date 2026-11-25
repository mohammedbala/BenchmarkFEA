# Example 01 - simply supported beam with design check

**File:** `01_simple_beam.json` (the same model in YAML: `01_simple_beam.yaml`).

## What it is

An 8 m floor beam on two supports, carrying a dead load of 12 kN/m and a live load of 18 kN/m. It is
the smallest complete model: two nodes, one member, two load cases, three combinations and a design
block.

Things worth noticing in the file:

* **Units and axes** are declared first: kN and m, vertical axis Z. Gravity loads are therefore
  negative `Z` values.
* **Supports in 3D.** The model lives in 3D space even though the beam is a plane problem. Node `A`
  is `pinned` (three translations) and also holds the twist `rx`; node `B` holds the vertical and the
  sideways translation and is free to move along the beam. With less than this the engine refuses the
  model and names the free motion.
* **Section.** A generic I-shape given by its dimensions (`type: "I"`): the engine computes the
  properties. With a fillet radius `r > 0` the design module treats it as a rolled shape.
* **`shearDeformation: false`** so that the deflection can be compared with the textbook formula
  (the default is `true`).
* **Combinations** carry a `purpose`: `strength` (default) feeds the strength checks, `serviceability`
  feeds the deflection limits.
* **Design block.** `bracing.topFlange` at the quarter points gives the unbraced length for
  lateral-torsional buckling (2 m); two deflection limits are checked, L/360 under live load and L/240
  under total load.

## Run it

```
benchfea validate examples/01_simple_beam.json
benchfea run      examples/01_simple_beam.json -o 01_simple_beam.results.json
benchfea design   examples/01_simple_beam.json -o 01_simple_beam.design.json
```

Numbers below are in the file's units (kN, m); signs follow docs/CONVENTIONS.md. They are what
the engine computes for this file (`tests/test_examples.py` recomputes every row); they are not
reference values from any publication.

## Key results

| quantity | value | unit |
|---|---|---|
| ULS: vertical reaction at A | 172.8 | kN |
| ULS: bending moment Mz at midspan | 345.6 | kN*m |
| ULS: shear force Vy at A | -172.8 | kN |
| SLS-L: midspan deflection (local y) | -0.01435 | m |
| SLS: midspan deflection (local y) | -0.02392 | m |
| design: flexure ratio (largest of the Chapter F checks) | 0.6629 | - |
| design: shear ratio | 0.2061 | - |
| design: deflection limit 1 ratio (L/360, live) | 0.6457 | - |
| design: deflection limit 2 ratio (L/240, total) | 0.7175 | - |
| design: governing ratio | 0.7175 | - |

Hand check (statics and the beam equation; derived here, checked in the test to 1e-10): with
w = 1.2 x 12 + 1.6 x 18 = 43.2 kN/m, the reaction is wL/2 = 172.8 kN and the midspan moment
wL^2/8 = 345.6 kN m; the live-load deflection is 5 w L^4 / (384 E Iz).

## Design summary

| design member | status | ratio | combination | governing limit state |
|---|---|---|---|---|
| B1 | conditional | 0.717 | SLS | deflection limit 2 (local y) |

Status `conditional` means: every check passes, but a strength combination was analysed first
order, and AISC 360-22 C1 asks for second-order effects in the required strengths (the engine says
so instead of amplifying silently).
