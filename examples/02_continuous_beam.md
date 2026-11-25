# Example 02 - two-span beam with a hinge and a partially fixed end

**File:** `02_continuous_beam.json` (kip, in; vertical axis Y).

## What it is

A beam over three supports with spans of 24 ft (A-B) and 30 ft (B-C):

```
   spring                          hinge
  A ~~~~~~~~~~~~~~~~~~~~ B ~~~~~~~ H o~~~~~~~~~~~~~~~~~~~~ C
  |<------ 288 in ------>|<- 90 ->|<------- 270 in ------->|
```

* **Partial fixity at A:** the support fixes the translations and holds the rotation about Z with a
  rotational spring, `"rz": {"spring": 400000}` (kip in / rad). The reaction moment at A is
  `-k x rotation`.
* **Internal hinge at H:** member `S2b` starts with `"releases": {"i": {"rz": true}}` - its bending
  moment about the local z axis is released at node H. The part H-C is then a simply supported
  "drop-in" span hanging on the cantilever B-H.
* **Pattern loading:** live load on span 1 (`L1`) and on span 2 (`L2`, with a 12 kip point load in
  the middle of `S2b`) are separate load cases, combined together and one at a time.
* **Settlement:** load case `S` imposes a displacement of 0.5 in downwards on support B. An imposed
  displacement is a load on a DOF that a support fixes, so it is factored and combined like any other
  load (combination `ULS4`).
* **Design member across the hinge:** span 2 is one design member made of `S2a` and `S2b`
  (`design.members.S2`); the group `beams` gives both spans the same bracing.

## Run it

```
benchfea validate examples/02_continuous_beam.json
benchfea run      examples/02_continuous_beam.json -o 02_continuous_beam.results.json
benchfea design   examples/02_continuous_beam.json -o 02_continuous_beam.design.json
```

Numbers below are in the file's units (kip, in); signs follow docs/CONVENTIONS.md. They are what
the engine computes for this file (`tests/test_examples.py` recomputes every row); they are not
reference values from any publication.

## Key results

| quantity | value | unit |
|---|---|---|
| ULS1: vertical reaction at A | 29.86 | kip |
| ULS1: vertical reaction at B | 169 | kip |
| ULS1: vertical reaction at C | 59.41 | kip |
| ULS1: moment taken by the rotational spring at A (reaction Mz) | 140.5 | kip*in |
| ULS1: rotation of the beam end at A | -0.0003513 | rad |
| ULS1: bending moment Mz over support B | -6841 | kip*in |
| ULS1: bending moment Mz at the hinge H | 0 | kip*in |
| ULS1: largest sagging moment in span 1 | 1063 | kip*in |
| ULS1: largest sagging moment in the drop-in part H-C | 4658 | kip*in |
| ULS4 (with 0.5 in settlement of B): bending moment Mz over support B | -6841 | kip*in |
| ULS4 (with 0.5 in settlement of B): moment taken by the spring at A | 593.9 | kip*in |
| ULS4 (with 0.5 in settlement of B): vertical reaction at B | 167.4 | kip |
| ULS1: deflection of the hinge H | -0.6109 | in |
| design: governing ratio, span 1 | 0.6404 | - |
| design: governing ratio, span 2 (S2a + S2b) | 0.6404 | - |

Hand checks (statics, derived here, in the test): the moment at the hinge is zero; the reaction at C
is that of the simply supported part H-C; the three vertical reactions add up to the total load.
B-H is a cantilever that carries the drop-in part at its tip, so the moment over B follows from
statics alone and the settlement of B does not change it (`ULS4` = `ULS1`). The settlement acts on
the statically indeterminate part: it changes the moment in the spring at A and the reactions of
span 1.

## Design summary

| design member | status | ratio | combination | governing limit state |
|---|---|---|---|---|
| S1 | conditional | 0.640 | ULS4 | flexure, major axis: yielding |
| S2 | conditional | 0.640 | ULS4 | flexure, major axis: yielding |

Status `conditional` means: every check passes, but a strength combination was analysed first
order, and AISC 360-22 C1 asks for second-order effects in the required strengths (the engine says
so instead of amplifying silently).
