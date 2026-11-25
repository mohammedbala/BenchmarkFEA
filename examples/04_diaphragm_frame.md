# Example 04 - 3D frame with a rigid diaphragm

**File:** `04_diaphragm_frame.json` (kip, in; vertical axis **Y**).

## What it is

A one-storey space frame, 240 in x 180 in in plan (X x Z) and 144 in high, with four fixed-base
columns and four beams.

* **Vertical axis Y.** Gravity acts along -Y; plan coordinates are X and Z. The default member
  orientation follows the vertical axis, so beams bend about their strong axis under gravity without
  any `orientation` entry, and all columns have their web along global X.
* **Rigid diaphragm.** `constraints.floor` ties the four floor nodes to the master `CM`: their two
  in-plane translations (`ux`, `uz`) and their rotation about the vertical axis (`ry`) follow the
  master as a rigid body; the vertical translation and the two other rotations of each node stay free.
  The diaphragm normal defaults to the vertical axis.
* **Master at the mass centre.** `CM` is a node without members at (132, 144, 108) - 12 in and 18 in
  away from the centre of the frame. The lateral loads `EX` and `EZ` are applied there, so the floor
  both translates and twists.
* **Results at a slave node** include `constraintForces`: the forces the diaphragm exerts on the node.

## Run it

```
benchfea validate examples/04_diaphragm_frame.json
benchfea run      examples/04_diaphragm_frame.json -o 04_diaphragm_frame.results.json
benchfea design   examples/04_diaphragm_frame.json -o 04_diaphragm_frame.design.json
```

Numbers below are in the file's units (kip, in); signs follow docs/CONVENTIONS.md. They are what
the engine computes for this file (`tests/test_examples.py` recomputes every row); they are not
reference values from any publication.

## Key results

| quantity | value | unit |
|---|---|---|
| SLS-X: displacement of the master CM along X | 0.5445 | in |
| SLS-X: displacement of the master CM along Z | -0.008063 | in |
| SLS-X: rotation of the floor about the vertical axis Y | 0.0006719 | rad |
| SLS-X: displacement of A1 along X (frame line z = 0) | 0.472 | in |
| SLS-X: displacement of D1 along X (frame line z = 180) | 0.5929 | in |
| SLS-X: base shear along X taken by columns CA + CB (line z = 0) | -8.863 | kip |
| SLS-X: base shear along X taken by columns CD + CC (line z = 180) | -11.14 | kip |
| ULS2: axial force in column CB (compression negative) | -18.09 | kip |
| design: governing ratio, column CB | 0.866 | - |
| design: governing ratio, beam BX2 | 0.4688 | - |

The two frame lines along X do not move equally (the floor rotates about Y), and the line nearer to
the load takes the larger share of the base shear. The test checks that every floor node moves
exactly as the rigid-body relation of docs/CONVENTIONS.md section 5 says, and that the base shears
add up to the applied 20 kip.

## Design summary

| design member | status | ratio | combination | governing limit state |
|---|---|---|---|---|
| CA | conditional | 0.788 | ULS3 | combined flexure and compression |
| CB | conditional | 0.866 | ULS3 | combined flexure and compression |
| CC | conditional | 0.861 | ULS3 | combined flexure and compression |
| CD | conditional | 0.849 | ULS3 | combined flexure and compression |
| BX1 | conditional | 0.435 | ULS1 | flexure, major axis: yielding |
| BX2 | conditional | 0.469 | ULS2 | combined flexure and compression |
| BZ1 | conditional | 0.201 | ULS3 | flexure, major axis: yielding |
| BZ2 | conditional | 0.219 | ULS3 | flexure, major axis: yielding |

Status `conditional` means: every check passes, but a strength combination was analysed first
order, and AISC 360-22 C1 asks for second-order effects in the required strengths (the engine says
so instead of amplifying silently).
