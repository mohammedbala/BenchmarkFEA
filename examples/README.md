# Worked examples

Five input files, from a single beam to a braced multi-storey frame. Each has a short explanation
next to it (`<name>.md`) with the points worth noticing in the file, the commands to run it and a
table of key results. All sections are generic parametric shapes given by their dimensions, so no
shapes database is needed.

| file | what it shows | units, vertical axis |
|---|---|---|
| [01_simple_beam.json](01_simple_beam.md) (also `01_simple_beam.yaml`) | the smallest complete model; supports of a plane member in 3D space; design check with flange bracing and deflection limits | kN, m; Z |
| [02_continuous_beam.json](02_continuous_beam.md) | partial fixity (rotational spring support), an internal hinge (member end release), pattern loading, a support settlement, a design member made of two analysis members | kip, in; Y |
| [03_portal_frame.json](03_portal_frame.md) | rigid end offsets, second-order analysis per combination, linear buckling, K factors | kN, mm; Z |
| [04_diaphragm_frame.json](04_diaphragm_frame.md) | a 3D frame with vertical axis Y, a rigid diaphragm with an eccentric master node, torsion of the floor | kip, in; Y |
| [05_braced_frame.json](05_braced_frame.md) | tension-only X-bracing, combinations with reversing wind, linear buckling, design groups | kN, m; Z |

```
benchfea validate examples/01_simple_beam.json
benchfea run      examples/01_simple_beam.json -o beam.results.json
benchfea design   examples/01_simple_beam.json -o beam.design.json
```

`expected.json` holds the numbers of the "Key results" tables at full precision;
`tests/test_examples.py` recomputes them, checks that the tables quote them, and checks the
closed-form results that exist (statics, the beam deflection formula, the rigid-body motion of the
diaphragm). The numbers are what the engine computes for these files; they are not taken from any
publication.

More models of the same kinds, in any number, come from the generators:
`benchfea generate --list` (docs/DIFFICULTY.md).
