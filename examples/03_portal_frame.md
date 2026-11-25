# Example 03 - portal frame, second-order analysis, rigid end offsets

**File:** `03_portal_frame.json` (kN, mm; vertical axis Z).

## What it is

A 12 m x 6 m portal frame with pinned bases in the X-Z plane; nodes `N2` and `N3` are held out of
plane (`uy`, `rx`, `rz` fixed).

* **Rigid end offsets.** The beam is connected to the column centre lines, but it is flexible only
  between the column faces: `"offsets": {"axes": "global", "i": [150, 0, 0], "j": [-150, 0, 0]}`
  (half the 300 mm column depth). Beam loads, stations and end forces refer to the flexible length of
  11 700 mm, so "end j" of `B1` is the face of column `C2`.
* **Second-order analysis per combination.** `ULS1` and `ULS2` have
  `"analysisType": "secondOrder"`; the engine iterates on the axial forces with exact beam-column
  stability functions. `ULS2-first-order` repeats `ULS2` first order so that the two can be compared
  in one run.
* **Linear buckling** of `ULS1` is requested in `analysis.buckling`.
* **Design:** K = 1 for the columns because the strength combinations are second order. The engine
  warns that this is not a complete direct analysis (no stiffness reduction, no notional load) - both
  can be switched on per combination (`stiffnessReduction`, `notionalLoad`); they are left out here
  to keep the example small.

## Run it

```
benchfea validate examples/03_portal_frame.json
benchfea run      examples/03_portal_frame.json -o 03_portal_frame.results.json
benchfea design   examples/03_portal_frame.json -o 03_portal_frame.design.json
```

Numbers below are in the file's units (kN, mm); signs follow docs/CONVENTIONS.md. They are what
the engine computes for this file (`tests/test_examples.py` recomputes every row); they are not
reference values from any publication.

## Key results

| quantity | value | unit |
|---|---|---|
| ULS2 (second order): sway of N2 | 45.29 | mm |
| ULS2-first-order: sway of N2 | 41.77 | mm |
| sway amplification second order / first order | 1.084 | - |
| ULS2 (second order): beam moment Mz at the face of column C2 (end j) | -2.579e+05 | kN*mm |
| ULS2-first-order: beam moment Mz at the face of column C2 (end j) | -2.5e+05 | kN*mm |
| ULS2 (second order): horizontal reaction at N4 | -46.68 | kN |
| ULS2 (second order): vertical reaction at N4 | 220.8 | kN |
| ULS1 (second order): largest sagging moment in the beam | 6.525e+05 | kN*mm |
| buckling factor of ULS1, mode 1 | 8.69 | - |
| buckling factor of ULS1, mode 2 | 51.09 | - |
| SLS: sway of N2 | 25.09 | mm |
| design: governing ratio, column C1 | 0.5369 | - |
| design: governing ratio, column C2 | 0.6298 | - |
| design: governing ratio, beam B1 | 0.684 | - |

The sway under wind grows by the amplification shown when the destabilising effect of the gravity
loads is included, and the moment at the leeward beam end grows with it. The buckling factor says the
gravity combination `ULS1` is that many times below the elastic sway buckling load.

## Design summary

| design member | status | ratio | combination | governing limit state |
|---|---|---|---|---|
| C1 | conditional | 0.537 | ULS1 | combined flexure and compression |
| C2 | conditional | 0.630 | ULS2 | combined flexure and compression |
| B1 | conditional | 0.684 | SLS | deflection limit 1 (local y) |

`conditional` here comes from the combination `ULS2-first-order`, which is a strength combination
analysed first order.

## Note on verification

This example combines second-order analysis with rigid end offsets and asks for a buckling analysis.
Neither foreign solver of the cross-check (docs/CROSSCHECK.md) can reproduce these with an error
estimate, so `benchfea verify` quarantines this file at tier C: its results are what the engine
computes, without an independent cross-check. Examples 01, 02 and 04 pass that tier.
