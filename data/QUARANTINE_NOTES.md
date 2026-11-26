# Notes on the quarantine log of this build

Written diagnoses of the attempts that the build logged with category `disagreement`, and of two
groups of `not-representable` attempts that share one cause. They belong to the manifest whose
SHA-256 begins with `c87dd884255cdd43`; the counts by category are generated on the data card
(`DATA_CARD.md`, "Quarantine log"). Every attempt listed here stays quarantined: nothing is served
from it, and no tolerance or check was changed because of it. The numbers were measured by running
the model of the attempt again (`benchfea.bench.tasks._base_model` with the generator seed of the
log entry, then `curveballs.apply` and `benchfea.verify.verify_model`); nobody has reviewed these
notes.

## Attempts logged as `disagreement` (5)

**1. `portal_frame` × `support_removed`, seed 0, level 3, attempt 1 (generator seed 1078338122) —
tier C. Cell not served.**
The engine and PyNite differ by 2.40e-4 on the shear of member C1 at end j in ULS1 (second order);
allowed 1.73e-4 (1e-6 plus the error bound PyNite's three refinement levels give). OpenSeesPy does
not converge on the second-order combinations (its Newton iteration stops at an increment norm of
4e-8 to 6e-7 against its test of 1e-10), so only one foreign solver is left in any case. Which side is
off: the independent tableau solver of `tests/independent` converges to the engine as its pieces
are refined — 2.8e-3, 1.7e-4, 1.1e-5, 6.6e-7 on the translations with 1, 4, 16, 64 pieces per
member, and its end forces agree with the engine to 2e-11 at every refinement. The difference is
therefore PyNite's discretisation error, and in this case it is 1.4 times the bound computed for
it: the bound of docs/CROSSCHECK.md 1.3 did not hold here. A bound that is too small makes tier C
stricter, not looser, so the outcome is a quarantine and not a wrong pass. The model is extreme:
with the support removed the frame hangs on one column, the first-order rotations reach 0.65 rad
and the elastic critical factor of ULS1 is 2.42.

**2. `braced_frame` × `load_direction_changed`, seed 4, level 3, attempt 1 (generator seed
1487383207) — tier D, station integration. Replaced by attempt 2.**
In ULS1 (second order) the deflection of brace X1 at midspan differs from the independent
integration of the member equations by 3.8e-6; tolerance 1e-6. X1 is a truss member, i.e. a
pin-ended strut carrying its self weight (docs/THEORY.md 10.6), and in this combination its axial
force is 357.3 kN in compression against an Euler load of 363.3 kN: 98 % of its own buckling load.
The deflection of such a strut is amplified by 1 / (1 − N/Pe) ≈ 60, and a relative change of N
changes it about 60 times as much; the second-order iteration stopped at a displacement change of
4.9e-8 and an axial-force change of 1.9e-8, which that factor turns into the order of the
difference found. In the other second-order combinations, where the brace is not near its Euler
load, the same check gives 9e-11 to 1e-9. Read as: the result is ill-conditioned at this load, not
that either computation is wrong; the estimate above explains the size but is not a proof. A model
with a brace at 98 % of its buckling load is not a reference to publish.

**3. `pipe_rack` × `offset_added`, seed 4, level 4, attempt 0 (generator seed 2031077189) —
tier E, rigid offsets; also tier C (not representable). Cell not served.**
Tier E replaces a rigid end offset by a stub member whose stiffness is raised tenfold per level and
compares the extrapolated limit with the engine: 1.17e-5 against a tolerance of 1e-5 (rotation of
node CL1, ULS5, second order). The differences to the engine fall about tenfold per level, as they must,
for the first five levels (ULS5: 5.8e-3, 5.9e-4, 5.9e-5, 5.9e-6, 9.1e-7) and then stop falling at
the sixth (9.7e-6); the seventh level does not converge at all in the second-order iteration. The
extrapolation uses the last three levels, so it takes in the sixth, where the stub is 1e6 times
stiffer than the member and the iterated solution of the stub model carries an error of 1e-6 to
1e-5. The first five levels approach the engine's result at the expected rate, to within 9e-7 to
1.6e-5 at the fifth level depending on the combination. This is a limit of the check in an iterated analysis with very
stiff stubs, not evidence against the engine's offsets; the attempt stays quarantined, and the
model could not be served anyway because no foreign solver gives an error bound for second order
with rigid offsets.

**4 and 5. `pipe_rack` × `release_mechanism`, seed 0 level 3 (generator seed 979462616) and seed 4
level 4 (generator seed 2031077189), attempt 0 — tier K. Both replaced by attempt 1.**
The diagnosis says: one mechanism, the rotation of a node about the global X axis, which no member
resists once the bending rotations of the three member ends there are released (the third member
has its torsion released at the other end). Tiers N, C and E pass: the independently assembled
stiffness has exactly a zero row and column for that rotation, both foreign solvers find no
equilibrium under a probing load, and the diagnosis is invariant. Tier K reports "K phi = 1.0 on
the engine's matrix". In the engine's matrix that row and column hold one entry, on the diagonal:
−1.9e-9 and −9.3e-10 in the two cases, against a largest entry of 1.6e10 and 1.7e10 — round-off of
the condensation of the released member ends, 1e-19 of the matrix scale. The engine counts the
mode correctly (one eigenvalue below the threshold). The residual of tier K is |K phi| divided by
|K| |phi| taken over the same column, so for a column that holds nothing but one round-off number
it is that number divided by itself. This is a defect of the check on the conservative side (in
the axis-aligned cases that pass, the entry is exactly zero and the residual is 0). The check was
not changed and the two attempts stay quarantined.

## `not-representable` attempts with a common cause

**PyNite stops with a division by zero on frames with crossing braces in second order (12
attempts with that error text, `braced_frame` and `pipe_rack`; one more with a failed refinement
level).** A truss member with a load along its length (its self
weight) is a pin-ended beam-column in the engine, and the translators subdivide it so that the
foreign solver can show its bending (docs/CROSSCHECK.md). The two diagonals of an X-brace are both
cut at mid-length, which puts two unconnected nodes at the crossing point. PyNite splits a member at
every node lying on it and so creates a piece of zero length between the two; its analysis stops
with `ZeroDivisionError`. The failure is recorded as "PyNite cannot represent the model", one
foreign solver is left, and the attempt is quarantined. Measured on one of them (`braced_frame` ×
`member_rolled_90`, seed 4, generator seed 1890047881: PyNite member `M10s3p0b` between nodes
`X81` and `X74`, both at the crossing). Consequence: second-order frames with crossing truss braces
under self weight are served only when a later attempt has no such brace, or not at all. Not
fixed in this build.

**OpenSeesPy does not converge (2 attempts, `portal_frame` × `support_removed`).** Other attempts
of the same curve ball on the same family as in note 1: a frame left hanging on one column.
