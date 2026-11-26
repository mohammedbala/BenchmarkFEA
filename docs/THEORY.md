# Theory of the linear solver core

Everything the solver uses is derived here from equilibrium, small-displacement
kinematics, linear elasticity and the principle of virtual work. No matrix or
formula is taken from a reference; each result names the test that verifies it
numerically against an independent solution of the governing differential
equations (`tests/solve_reference.py`: Chebyshev collocation of the first-order
system, cross-checked against `scipy.integrate.solve_bvp`).

Axis and sign conventions are those of [CONVENTIONS.md](CONVENTIONS.md); code
is in `src/benchfea/solve/`.

Sections: 1 element · 2 end releases · 3 rigid end offsets · 4 constraints ·
5 loads, fixed-end forces and recovery · 6 supports and solution ·
7 stability diagnosis · 8 invariants · 9 hooks for nonlinear analysis ·
10 second-order member · 11 linear buckling · 12 second-order analysis of a
combination · 13 tension-only and compression-only members · 14 stiffness
reduction and notional loads · 15 invariants of a second-order solution ·
16 scope and limits of the nonlinear analyses.

## 1. Element

### 1.1 Member equations

Local axis x along the member, y and z the section axes. Displacement of the
axis `u, v, w`, rotation of the section `θx, θy, θz` (right-hand rule). Plane
sections remain plane, so a point (y, z) of the section moves along x by

    u_x(y, z) = u − y θz + z θy,        ε = u′ − y θz′ + z θy′,     σ = E ε.

With the stress resultants of CONVENTIONS §4 (on the face with outward normal
+x): `N = ∫σ dA`, `Mz = −∫ y σ dA`, `My = +∫ z σ dA`, and y, z principal
centroidal axes (∫y dA = ∫z dA = ∫yz dA = 0):

    N = EA u′,      Mz = EIz θz′,      My = EIy θy′.

Shear strains of the axis: `γxy = v′ − θz`, `γxz = w′ + θy`; with the shear
areas `Asy, Asz`: `Vy = G Asy (v′ − θz)`, `Vz = G Asz (w′ + θy)`. Torsion
(uniform, St. Venant): `T = GJ θx′`. Hence the kinematic–elastic equations

    u′ = N / EA                 θx′ = T / GJ
    θz′ = Mz / EIz              v′ = θz + Vy / (G Asy)
    θy′ = My / EIy              w′ = −θy + Vz / (G Asz)                  (1)

Without shear deformation (option off, or no shear area) the terms `V / G As`
are dropped: `v′ = θz`, `w′ = −θy`, i.e. `Mz = EIz v″`, `My = −EIy w″`.

Equilibrium of a slice dx under line loads `(wx, wy, wz)` per unit length
(force balance; moment balance about the z and y axes through the slice):

    N′ = −wx     Vy′ = −wy     Vz′ = −wz     T′ = 0
    Mz′ = −Vy    My′ = +Vz                                                 (2)

A point force `P` at x = a makes the force resultants drop by `P` across a
(`s(a⁺) = s(a⁻) − P`), a point moment likewise for the moment resultants;
displacements and rotations are continuous.

### 1.2 Exact integration (`beam.py`)

Loads are point loads and piecewise linear line loads. Between two consecutive
load break points the right-hand sides of (2) are polynomials of degree ≤ 1,
so `N, V` have degree ≤ 2, `M` ≤ 3, rotations ≤ 4, deflections ≤ 5. `integrate`
marches from x = 0 with given start values, integrating (2) then (1) in closed
form on each segment (polynomial coefficients), applying the jumps at point
loads. There is no quadrature and no shape function; the only error is
round-off. *Verified:* `test_integrator_reproduces_the_ode_solution_*`.

### 1.3 Cantilever flexibility

Clamp end i, load end j by `Q = (Fx, Fy, Fz, Mx, My, Mz)` (forces on the
member). Equilibrium of the part between x and L gives the resultants

    N = Fx    Vy = Fy    Vz = Fz    T = Mx
    My(x) = My_j − (L − x) Fz        Mz(x) = Mz_j + (L − x) Fy            (3)

i.e. `s(x) = B(x) Q`. By the principle of complementary virtual work, for any
virtual end load δQ the displacement d of end j satisfies

    δQᵀ d = ∫₀ᴸ δsᵀ C s dx,   C = diag(1/EA, 1/GAsy, 1/GAsz, 1/GJ, 1/EIy, 1/EIz)

so `d = F_b Q` with `F_b = ∫ Bᵀ C B dx`. Carrying out the integrals with (3)
(order x, y, z, rx, ry, rz):

    F_b[x,x]   = L/EA                    F_b[rx,rx] = L/GJ
    F_b[y,y]   = L³/3EIz + L/GAsy        F_b[y,rz]  = L²/2EIz      F_b[rz,rz] = L/EIz
    F_b[z,z]   = L³/3EIy + L/GAsz        F_b[z,ry]  = −L²/2EIy     F_b[ry,ry] = L/EIy

(the sign of the z–ry coupling follows from the minus sign in `My` of (3)).
*Verified:* the closed form equals the exact integration of unit end loads and
the independent ODE solution
(`test_closed_form_cantilever_flexibility_equals_exact_integration_of_unit_end_loads`,
`test_cantilever_flexibility_against_the_independent_ode_solution`).

### 1.4 Stiffness of the 12 end DOFs

Let `H` be the rigid-body transfer from end i to end j,
`(u_j, θ_j)_rigid = H (u_i, θ_i)` with `u_j = u_i + θ_i × (L e_x)`. The
deformation of the cantilever is `d = u_j − H u_i = Γ u`, `Γ = [−H, I]`, and
`Q_j = K₀ d` with `K₀ = F_b⁻¹`. Equilibrium of the unloaded member reduced to
end i gives `Q_i = −Hᵀ Q_j` (force transfer is the transpose of the kinematic
transfer). Therefore

    f = [Q_i; Q_j] = Γᵀ K₀ Γ u = K u.

`K` is symmetric, positive semi-definite, and its null space is exactly the
six rigid-body motions (`Γ u = 0`). This is the Timoshenko-type member
stiffness (it contains the shear flexibility through `F_b`); it is never
written out entry by entry. *Verified:* every column against the ODE solution
with unit end displacements, with and without shear deformation
(`test_stiffness_columns_are_the_end_forces_of_the_ode_solution_with_unit_end_displacements`),
and the null space
(`test_stiffness_is_symmetric_positive_semidefinite_with_exactly_six_rigid_body_modes`).

### 1.5 Truss-type members

A truss member transmits axial force only: it is the member of 1.4 with both
bending rotations released at both ends and the twist released at one end. Its
stiffness between the end nodes is then `EA/L` on the axial DOFs and exactly
zero elsewhere, which is what the code stores (exact zeros, so that the
rotations of truss-only nodes can be recognised as having no stiffness).
Loads on the member are carried by it as a simply supported beam (section 5).

## 2. End releases (`element.release_stiffness`)

A release is a spring of stiffness `k ∈ [0, ∞]` between the node side and the
beam side of one end DOF (`k = ∞` rigid, `k = 0` fully released). The spring
transmits the same resultant as the section it sits at, and its deformation is
`e = resultant / k` (beam side minus node side at i; node side minus beam side
at j, so that a positive resultant stretches it — the same sign rule as (1)).

Two exact, algebraically equivalent ways to include a spring are used,
chosen by magnitude so that round-off never matters:

**Stiff springs** (`k ≥` the beam's own diagonal stiffness in that end DOF):
springs in series add flexibilities. A deformation `e` of spring a moves end j
of the cantilever by `s_a e`, where `s_a` is column a of `I` (spring at j) or
of `H` (spring at i: the whole beam follows rigidly); the spring force is
`s_aᵀ Q`. Hence `F₀ = F_b + Σ (1/k_a) s_a s_aᵀ`, `K₀ = F₀⁻¹`. `k = ∞` adds
nothing.

**Soft springs** (`k` below that, including 0): keep the beam-side DOFs `b`
as extra unknowns. The stiffness `Γᵀ K₀ Γ` then refers to the node-side DOFs
`r` of the rigid/stiff connections and to the beam-side DOFs `b`; each soft
spring adds `k (u_n − u_b)²/2` to the energy. Eliminating `u_b` (no external
force acts on it) gives

    K_bb* = K_bb + diag(k)
    K_rr ← K_rr − K_rb K_bb*⁻¹ K_br
    K_nr = diag(k) K_bb*⁻¹ K_br
    K_nn = diag(k) − diag(k) K_bb*⁻¹ diag(k)

For `k = 0` rows and columns `n` are exactly zero (a full release transmits
nothing), for small k the node sees `k` minus a correction of relative size
`k / K_bb`, so nothing large is subtracted. For `k → ∞` the first form is used
and `1/k → 0`. *Verified:* partial releases over ten decades against the ODE
with a spring boundary condition
(`test_partial_rotational_release_*`, `test_partial_translational_release_*`),
limits and monotonicity
(`test_release_limits_k_to_zero_and_k_to_infinity_*`,
`test_partial_release_tends_to_the_rigid_and_to_the_pinned_solution_and_is_monotonic`).

**Member-level mechanism.** `K_bb*` is singular iff there is a motion of the
beam-side DOFs with no strain energy in the beam and no energy in the soft
springs: a rigid-body motion of the member that the retained connections do
not prevent (twist released at both ends, the same translation at both ends,
both end moments plus a shear, …). Such a member is refused
(`member.release-mechanism`, naming the releases from the null vector); the
test is `λ_min < 1e-12 λ_max` on the unit-diagonal-scaled `K_bb*`, which also
refuses springs so soft that the member is a mechanism to round-off.

## 3. Rigid end offsets

The flexible length starts at `node i + o_i` and ends at `node j + o_j`. The
piece between node and flexible end is rigid, so for small displacements

    u_end = u_node + θ_node × o,      θ_end = θ_node,

i.e. `(u, θ)_end = E(o) (u, θ)_node` with `E = [[I, −[o]×], [0, I]]`. With the
rotation `R` to local axes, `u_local = G u_nodes`, `G = blockdiag(R E(o_i),
R E(o_j))` per end. Virtual work must be the same at the nodes and at the
ends, so nodal forces are `Gᵀ f` and `K_nodes = Gᵀ K G` (this contains the
lever arm `M_node = M_end + o × F_end`). Member loads act on the flexible
length only; the rigid pieces carry no load and no self weight.
*Verified:* `test_rigid_offsets_are_the_exact_small_displacement_kinematics`,
`test_offset_not_along_the_axis_eccentric_beam_statics_and_kinematics`,
and convergence of a stiff stub member to the offset
(`test_offset_converges_to_a_very_stiff_stub_member_with_an_added_node`).

## 4. Rigid links and diaphragms as exact constraints (`system.py`)

A rigid link makes the slave follow the master as a rigid body:
`u_s = u_m + θ_m × r`, `θ_s = θ_m`, `r = x_s − x_m`. A rigid diaphragm with
normal n imposes this for the two in-plane translations and the rotation about
n only, using the in-plane part of r and only the master's in-plane
translations and rotation about n.

These are linear homogeneous relations `u = T u_r` between all nodal DOFs and
the retained DOFs (all except the slave DOFs); imposed displacements enter
separately as prescribed values of retained DOFs, so the general form
`u = T u_r + g` is realised with `g` inside `u_r`. Substituting into the
potential energy `½uᵀKu − uᵀF` and making it stationary over `u_r`:

    Tᵀ K T u_r = Tᵀ F.

Slave DOFs are eliminated exactly; no stiff member and no penalty is used.
The force each constrained node receives from the rigid body is the residual
of its own equilibrium, `λ = K u − F − R_support` on the slave and master DOFs;
`Tᵀ λ = 0` on the free DOFs (the constraint forces do no work), and for a
rigid link or an in-plane diaphragm they are self-equilibrated.

A diaphragm slave that is not in the master's plane receives in-plane forces
at its own level, whose moment about in-plane axes nothing balances; the
global moment check then fails by exactly that amount and the result is
refused/flagged (section 8). *Verified:*
`test_rigid_link_is_the_exact_rigid_body_relation_and_reports_constraint_forces`,
`test_rigid_link_is_the_limit_of_a_very_stiff_member`,
`test_rigid_diaphragm_*`.

## 5. Loads, fixed-end forces, recovery

### 5.1 Member loads

The model layer hands over member loads in local axes on the flexible length
(global directions rotated by R; a projected load multiplied by
`|x_local × direction|`). Self weight is the line load
`factor · γ · A · (R · direction)` on the flexible length. Truss members
receive member loads too (not torques about their axis: refused).

### 5.2 Fixed-end forces

Hold both nodes. Take the same cantilever as in 1.3 with the span loads on it:
the resultants at its clamped end are the load resultant `P₀` (reduced to end
i), and exact integration (1.2) gives the displacement `u₀` of its free end
(plus `H_a P₀_a / k_a` for each stiff spring at i, which carries `P₀`).
Compatibility requires end j to return to zero: `Q_j = −K₀ u₀`, and equilibrium
gives `Q_i = −Hᵀ Q_j − P₀`. These are the fixed-end forces with the soft
beam-side DOFs still held; releasing them by the condensation of section 2,
`u_b = −K_bb*⁻¹ f_b`, `f ← f + K_·b u_b`, gives the fixed-end forces at the
node-side DOFs (`f_n = −k u_b`, exactly zero for a full release). The nodes
receive `−Gᵀ f_fixed`. The same derivation serves every load type, release
pattern and shear option, so no table of fixed-end formulas exists in the
code. *Verified:* `test_fixed_end_forces_against_the_ode_solution` (24
load/connection combinations), `test_equivalent_nodal_loads_balance_*`.

### 5.3 Recovery, stations

With the nodal displacements known: end displacements `G u`; beam-side
displacements of soft connections from the condensation; end forces
`f = K u + f_fixed` (forces on the member, local). The spring deformations
give the start values of the beam at x = 0; the member equations are then
integrated exactly from end i with `s(0⁻) = −Q_i` (1.2), which yields the
resultants and the displacements of the axis at any x, including both sides of
a point load. At x = L the integration must reproduce the end forces at j and
the displacement of end j — two closure checks that are reported as invariants
on every solve.

Stations: the requested number of equally spaced stations, the ends of line
loads, and two stations (left and right limit) at every interior point load of
any load case, so that every combination uses the same stations. The station
at 0 is the limit from outside (= minus the end forces at i), the station at L
likewise (= the end forces at j).

### 5.4 Combinations

All loads are linear, so a combination is the factored sum of the case load
sets; it is solved directly with the factorised stiffness. That this equals
the superposition of separately solved cases is checked by
`invariants.superposition_check`.

## 6. Supports and solution

Fixed DOFs are eliminated (their displacement is zero or the value imposed by
`supportDisplacement` loads); springs to ground add `k` to the diagonal.
With f = free and p = prescribed retained DOFs:

    K_ff u_f = F_f − K_fp u_p,       R_p = K_pf u_f + K_pp u_p − F_p.

Reactions are forces on the structure: `R_p` at fixed DOFs (for a master this
includes what its slaves transfer) and `−k u` at springs. `K_ff` is scaled to
unit diagonal and factorised by Cholesky, with two steps of iterative
refinement.

DOFs that have no stiffness at all and are not mechanisms of the structure are
taken out and reported as zero with a warning: all free DOFs of a node without
members (unless it is a master driving slaves), and the rotations of a node to
which only truss members without offsets connect. A load on such a DOF is
refused.

## 7. Stability diagnosis (`stability.py`)

The structure is refused when the smallest eigenvalue of the unit-diagonal
scaled `K_ff` is below `analysis.pivotTolerance` (or a Cholesky pivot is, or
the factorisation fails). With scaled eigenvalue ≤ 1e-12 the motion has zero
stiffness to round-off ("mechanism"), otherwise the matrix is
"ill-conditioned" (stiffness contrast ≈ 1/eigenvalue).

The null space is reported through basis-independent constructions:

1. *Rigid-body mechanisms*: for each connected part, the combinations of its 6
   rigid-body motions that move no fixed or sprung DOF and violate no
   constraint — the null space of a small matrix, put in reduced row echelon
   form over (TX, TY, TZ, RX, RY, RZ) about a supported node. Each is described
   as a translation, a rotation about an axis through a point, or both
   (cause `no-restraint`). Then the rigid-body motions that move no *fixed*
   DOF but do stretch support springs are tested against the numerical null
   space: those inside it are reported with the cause `negligible-springs`
   (the springs exist but are below the tolerance relative to the stiffness of
   the structure), instead of being left to item 3.
2. *Single-node mechanisms*: null vectors confined to one node (eigenvalue 1
   of `Q_aᵀQ_a` for an orthonormal null basis Q restricted to the node's
   rows), again in row echelon form, with the cause: all connected member ends
   released, only truss members, a member free to spin about its own axis
   (`member-spin`: every member at the node lies along the rotation axis, so
   the null vector is the twist of those members, held at neither end), or no
   stiffness.
3. *The rest*: the participation of every node and DOF = diagonal of the
   orthogonal projector onto the remaining null space (independent of the
   basis), with the members attached to the participating nodes.

Rotations are multiplied by the model size to compare them with translations.

## 8. Invariants on every solve (`invariants.py`)

| check | quantity | tolerance |
|---|---|---|
| stiffness symmetry | max\|K − Kᵀ\| / max\|K\| | 1e-10 |
| rigid-body modes | max\|K B\| / max(\|K\|\|B\|), B = 6 rigid-body motions of the unsupported, unconstrained structure | 1e-9 |
| positive definiteness | smallest eigenvalue of the scaled `K_ff` (and Cholesky pivots) | ≥ `pivotTolerance` |
| linear solve residual | max\|K_ff u_f − rhs\| / max(\|K_ff\|\|u_f\| + \|rhs\|) | 1e-9 |
| global force / moment equilibrium | applied loads (from their definitions, moments about the origin) + reactions | `equilibriumTolerance` |
| joint equilibrium | nodal loads + reactions + constraint forces − member end forces, every node | `equilibriumTolerance` |
| member equilibrium | integrated resultants at x = L − end forces at j | `equilibriumTolerance` |
| member compatibility | integrated displacements at x = L − displacement of end j | `equilibriumTolerance` |
| work = strain energy | uᵀK_members u − (loads + reactions)·u | `equilibriumTolerance` |

Force-type residuals are measured against the largest force-type quantity of
the check (at least moment scale / model size), moment-type residuals against
the largest moment-type quantity (at least force scale × model size).
A failed check makes `analyze` raise `ModelRefused`
(`solve.invariant-failed`); with `strict=False` the results are returned
flagged. `superposition_check` and `reciprocity_check` (Maxwell–Betti:
`δ_ab = δ_ba` for unit loads) are available to tests and to `verify`.

### 8.1 Round-off floor of the force scales

Every force of a solution is evaluated as a sum of products stiffness ×
displacement (minus fixed-end forces). When the true forces vanish — a
statically determinate structure under a support settlement moves as a rigid
body — the computed "forces" are the rounding errors of those sums, and a
residual divided by them is a ratio of two rounding errors. The quantity that
keeps a meaning is the size of the terms summed,

    S_f = max over the force rows of (|K| |u|),
    S_m = max( max over the moment rows of (|K| |u|),  S_f × model size ),
    S_e = |u|ᵀ |K| |u|                                   (work check),

with K the assembled member stiffness (the tangent stiffness in second
order). A sum of n products evaluated in floating point is in error by at
most about n ε Σ|terms| (ε = 2.2e-16; each product and each addition is
rounded once); a row of an element stiffness has n = 12 terms.
`ROUND_OFF = 16 ε = 3.6e-15` is taken. A residual is accepted when

    residual ≤ tolerance × (force scale of the check)   or   residual ≤ ROUND_OFF × S,

i.e. the reported value is `residual / max(force scale, ROUND_OFF × S / tolerance)`.
The floor never reduces a scale, so no check is looser than before wherever
the forces of the check exceed `ROUND_OFF × S / tolerance` (3.6e-9 S with the
default tolerance 1e-6). *Measured* on the model bank (27 + 8 + 3 models of
`tests/independent`, `examples/`, `verification/problems/`, valid fixtures;
218 models, 2885 residual ratios): 11 ratios changed, all of them in the
rigid-body settlement frame (`verification/problems/staad/frame-support-settlement`,
before 7e-4 … 1.0, after 9e-12 … 3.5e-8); every other ratio is bit-for-bit
the old one; the loaded case closest to the floor has a force scale 7.4 times
the floor (an "axially rigid" frame with A = 1e8 mm²). The largest rounding
error observed in the settlement frame is 0.6 ε S.

### 8.2 Superposition check

For every combination the factored sum of the separately solved cases is
compared with the combination solved directly, array by array: nodal
displacements; reaction forces; reaction moments; and per member the station
forces (N, Vy, Vz), the station moments (T, My, Mz) and the station
displacements. Forces and moments are separate arrays because their units
differ (in N, mm a moment array would otherwise hide every force). The value
of the check is the largest `max|sum − direct| / scale` over the arrays, with
`scale = max|direct|` of the array for displacements and

    scale = max( max|direct| of the array,  ROUND_OFF × S / tolerance )

for force arrays, where S is the stiffness × displacement scale of that array
— for a member `max(|k| |u_ends| + |fixed-end forces|)` over its 12 end DOFs
in local axes (moments: at least force scale × member length), for the
reactions S_f, S_m of 8.1 — summed over the direct solution and |factor| ×
each case. A member that carries no force in a combination while its nodes
move (or whose case forces cancel) has "forces" that are rounding errors in
both results; before this scale the check divided one by the other (values
0.6 … 2.3 against a tolerance of 1e-9). With the tolerance 1e-9 of tier D an
array is measured against its own size as soon as its forces exceed 3.6e-6 S;
below that, a difference is accepted only up to `ROUND_OFF × S`, the rounding
bound of the two evaluations compared. *Measured* (268 models: the bank of 8.1
plus two generated models per family and level, 17 805 arrays): 554 arrays
are measured against the floor, the largest difference among them is 5e-17 S;
3 old failures (all of the round-off kind), no new failure, largest new value
1.6e-11. Tests corrupt a result and show that a wrong moment of a loaded
member is caught at the same size as before (relative 3e-9 caught, 3e-10 not)
and that a spurious moment of 1e-12 S in a force-free member is caught.

## 9. Hooks for nonlinear analysis

`ANALYSIS_TYPES[type](analysis, combination)` solves one combination;
`firstOrder` (this core) and `secondOrder` (sections 10–15) are registered. A
nonlinear driver reuses
`Structure.assemble(matrices)` for the geometric stiffness (per-member 12×12
matrices in nodal global DOFs, `Gᵀ k_g G`), `Structure.factorize` (stability
diagnosis included), `Structure.solve`, `LinearAnalysis.recover`, and
`build_element(model, member, factors, active)` for reduced stiffness and for
members that are switched off.

## 10. Second-order member (`beamcolumn.py`)

As in sections 1–9, everything is derived here; the verification is against
independent integrations of the same differential equations
(`tests/second_order_reference.py`: `scipy.integrate.solve_bvp`, and the
classical fourth-order form of the equation as a second witness).

### 10.1 Assumptions

Small-displacement second-order theory: displacements and rotations are small
compared with the member length, but the equilibrium of a slice is written in
its displaced position. Of the internal forces only the **axial force N**
(tension positive) is allowed to act on the displacements; it is taken from the
previous iterate (section 12), keeps the direction of the undeformed member
axis, and may vary along the member (axial member loads). The axial and
torsional behaviour stay first order: `u′ = N/EA`, `N′ = −w_x`, `θx′ = T/GJ`.

### 10.2 Equations of one bending plane

Plane xy, deflection `v`, section rotation `θ = θz`, `M = Mz`, and **V = the
force component along the undeformed y axis** on the +x face. Slice between x
and x + dx in its displaced position; on its +x face `(N, V, M)` at the point
`(x + dx, v + dv)`, on its −x face `(−N, −V, −M)` at `(x, v)`, line load
`q = w_y`:

    force along y:              V′ = −q
    moment about the −x face:   dM + V dx − N dv = 0   ⇒   M′ = −V + N v′        (4)

(first order: `M′ = −V`, equation (2)). Elasticity: `θ′ = M/EI`.

*Shear deformation.* The second-order work of the axial force is that of N on
the shortening of the axis, `∫ v′²/2 dx`, so the potential energy of the plane
problem is

    Π = ∫ [ ½ EI θ′² + ½ G As (v′ − θ)² + ½ N v′² − q v ] dx .

Stationarity with respect to θ gives `M′ = −G As (v′ − θ)`, with respect to v
`(G As (v′ − θ) + N v′)′ = −q`. Comparing with (4): the shear strain is
produced by `Q = V − N v′`, the force **normal to the deformed axis**,

    v′ − θ = β (V − N v′),   β = 1/(G As)      ⇒      v′ = (θ + β V) / (1 + β N).

Together:

    v′ = (θ + βV)/(1 + βN)    θ′ = M/EI    M′ = (Nθ − V)/(1 + βN)    V′ = −q      (5)

`1 + βN > 0` is required (`N > −G As`); otherwise the piece is unstable in
shear and is treated as beyond its critical load. This is the classical
Engesser form of the shear-deformable beam-column. The other consistent choice
(shear strain produced by the force normal to the rotated *section*,
`V − Nθ`) differs by terms of relative order `N/(G As)` times the shear
deformation; see section 16. β = 0 gives `EI v⁗ − (N v′)′ = q`.

The xz plane obeys the same equations in the variables
`(w, −θy, −My, Vz)` with `q = w_z` (section 1.1: `w′ = −θy + …`,
`θy′ = My/EIy`, `My′ = +Vz` map onto `v′ = θ + …`, `θ′ = M/EI`, `M′ = −V`
under this substitution, and (4) is derived in the same way).
*Verified:* `test_second_order_in_the_other_bending_plane_and_for_a_rolled_member`.

### 10.3 Transfer matrix, stiffness, fixed-end forces

(5) is a linear system `z′ = A(N) z + g(x)`, `z = (v, θ, M, V)`. A piece
carries a line load `q = q_a + q_s x` over its whole length h and no load
inside (section 12 puts a piece boundary at every load break). With the two
extra states `(1, x)` the load becomes part of a homogeneous system
`y′ = Â y`, `y = (z, 1, x)`. For constant N the solution operator is the
matrix exponential,

    y(x) = exp(Â x) y(0),

which is exact for compression (trigonometric), tension (hyperbolic) and
N = 0 (polynomial) alike — the "stability functions" in numerical form, with
no special cases and no cancellation near N = 0. The state is scaled to
`(v/h, θ, M h/EI, V h²/EI)`; the exponential is evaluated by scaling and
squaring with a Taylor series (`beamcolumn.expm`, checked against SciPy).

Write the transfer over the piece as `d₁ = A d₀ + B s₀ + p_d`,
`s₁ = C d₀ + D s₀ + p_s` with `d = (v, θ)` and the work-conjugate
`s = (V, M)`. The forces on the piece are `−s₀` at the start and `+s₁` at the
end (CONVENTIONS §4). Solving the first relation for `s₀`:

    [F₀]   [ B⁻¹A          −B⁻¹  ] [d₀]   [ B⁻¹ p_d          ]
    [F₁] = [ C − D B⁻¹A     D B⁻¹ ] [d₁] + [ p_s − D B⁻¹ p_d  ]

— the second-order stiffness and the fixed-end forces of the piece, including
the effect of the axial force on the loads inside the piece (P-δ). `B` is
singular only at the critical load of the clamped–clamped piece; the pieces
are kept far shorter than that (section 12). The matrix is not symmetrised:
its symmetry is an invariant that is checked on every solve. *Verified:* for
N = 0 it reproduces the first-order element of sections 1–5 to round-off,
stiffness, fixed-end forces and stations
(`test_without_axial_force_the_piece_is_the_first_order_element`); with N,
every station value against `solve_bvp`
(`tests/test_second_order_beamcolumn.py`, P/Pe from 0.1 to 0.9, tension, with
and without shear deformation, point loads, moments, partial trapezoids).

### 10.4 Axial force varying along the piece

With axial member loads `N(x) = N₀ − w_a x − w_s x²/2` and `Â` depends on x.
The solution of `y′ = Â(x) y` over a step Δ is the time-ordered series
`Φ = I + ∫Â + ∫∫_{s₂<s₁} Â(s₁)Â(s₂) + …`. Writing `Φ = exp(Ω₁ + Ω₂ + …)`
and comparing orders: `Ω₁ = ∫Â ds` and, since `Ω₁²/2` contains the symmetric
part of the double integral, `Ω₂ = ½ ∫∫_{s₂<s₁} [Â(s₁), Â(s₂)]`. Expanding Â
about the midpoint, `Â(s) = Â_m + (s − m) Â′ + …`:

    Ω₁ = Δ Â_m + O(Δ³),     Ω₂ = −(Δ³/12) [Â_m, Â′_m] + O(Δ⁵),     Ω₃ = O(Δ⁵).

With the two Gauss points `s₁,₂ = m ∓ Δ√3/6`: `(Â₁ + Â₂)Δ/2 = ∫Â + O(Δ⁵)`
and `[Â₁, Â₂] = (√3Δ/3)[Â_m, Â′_m] + O(Δ³)`, hence

    Ω = (Δ/2)(Â₁ + Â₂) − (√3 Δ²/12) [Â₁, Â₂],      Φ = exp(Ω) + O(Δ⁵).

A piece with varying N is crossed in 8 such steps (`MAGNUS_STEPS`); the
global error is O(Δ⁴). For constant N the formula reduces to the exact
exponential. *Verified:*
`test_transfer_with_a_varying_axial_force_is_fourth_order_accurate`,
`test_axial_force_varying_along_the_member` (column under its own weight,
error below 1e-8), `test_axial_force_varying_along_the_member_self_weight_buckling`.

### 10.5 Axial and torsional part

Uniform bar, first order: stiffness `EA/h` (and `GJ/h`) between the ends. Bar
with both ends held under `w_x = w_a + w_s x`: by virtual work with the two
linear displacement fields `1 − x/h` and `x/h` (the exact homogeneous
solutions, so the result is exact) the end forces on the bar are
`−∫ w_x (1 − x/h) dx = −(w_a h/2 + w_s h²/6)` and
`−∫ w_x x/h dx = −(w_a h/2 + w_s h²/3)`. Point loads at the two ends of a
piece act on the held end DOFs directly.

### 10.6 End connections

A release or partial release is a spring `k` between the node side and the
beam side of an end DOF (section 2). With `c` = the spring-connected DOFs,
`r` the rigidly connected ones, `D = diag(k)`, `K` and `f₀` the stiffness and
fixed-end forces at the beam-side DOFs, the beam side carries the spring force
`D (u_n − u_b)`:

    (K_cc + D) u_b = D u_n − K_cr u_r − f₀_c ,          S = (K_cc + D)⁻¹
    K_rr ← K_rr − K_rc S K_cr        K_nr = D S K_cr        K_nn = D S K_cc
    f_r ← f₀_r − K_rc S f₀_c         f_n = D S f₀_c

(`K_nn = D − D S D = D S K_cc`: the second form has no cancellation for stiff
springs; for `k = 0` the rows are exactly zero.) `K_cc + D` is positive
definite as long as the piece itself has not buckled between its end
connections (e.g. a pin-ended member below its own critical load). When it is
not, the piece is flagged and the combination is refused (section 12): this
instability is invisible in the condensed matrix.

**Transverse translational releases** (`uy`, `uz` released or on a spring) are
refused in second-order and buckling analyses: across such a sliding
connection the axial force acts with an eccentricity equal to the slip, and
the moment this produces depends on how the connection is built (which side
carries the guide), which the model does not say.

**Truss-type members** (`truss`, `tensionOnly`, `compressionOnly`) are the
piece with both bending rotations released at the original member ends and no
torsional stiffness. They are therefore pin-ended beam-columns with the
`Iy`, `Iz` of their section: a transverse load on a compressed truss member is
amplified, and its own flexural buckling is part of the analysis. The twist of
the internal nodes of a subdivided truss member is not a degree of freedom of
the theory (no torque can reach it) and is held.

Consequences, intended and to be known when modelling:

* In a first-order analysis the bending stiffness of a truss-type member
  affects only the deflection of its own axis under its member loads; a
  placeholder `Iy`, `Iz` is harmless there.
* In a second-order or buckling analysis `Iy`, `Iz` are physical data. The
  own Euler load of the member is `π² E I_min / L²` (pin-ended, L = flexible
  length, stiffness reduction included). A linear buckling analysis returns
  the member's own factor `π² E I_min / (L² |N|)` when that is the lowest
  mode — with a placeholder `I = 1 mm⁴` this is a factor of order 1e-4, which
  is the correct answer for the model as given and not the buckling factor of
  the truss. A second-order analysis refuses a member beyond its own Euler
  load (`stability.beyond-critical-load`).
* The engine says so: when the lowest buckling factor is within 1 % of the
  own Euler factor of a compressed truss-type member, the buckling result
  carries a warning naming the member, its Euler load and its axial force
  (`buckling.truss_strut_warning`); the refusal of a second-order
  combination carries the same sentence. Results are not changed.
* To study the buckling of a truss as a system of bars only, give the bars an
  `I` large enough that their own factors lie above the modes of interest
  (the exact value then has no influence: in
  `verification/problems/mgz/ex-9-1-truss-critical-load` the factor is the
  same for `I` between 1e8 and 1e12 mm⁴).
* The input schema requires positive `Iy`, `Iz` for a `properties` section, so
  "no I" cannot occur; whether a given `I` is realistic cannot be decided by
  the engine (a flat-bar tie has a very small one), which is why there is a
  warning at the point where it governs and no validation rule.

### 10.7 Rigid end offsets

`u_end = u + θ × o` (section 3) is the first-order kinematics of the rigid
arm. To second order in the rotation vector θ the arm end moves by
`θ × o + ½ θ × (θ × o)`. The axial end force `f` that the arm applies to the
member end (`f = −N e_x` at end i, `+N e_x` at end j) does the second-order
work `f · ½ θ × (θ × o) = ½[(θ·o)(θ·f) − |θ|² (o·f)]`, i.e. the arm adds

    K_arm = ½ (o fᵀ + f oᵀ) − (o·f) I

to the rotational DOFs of the node (for an offset along the member in
compression P: `−P a` on the two transverse rotations — the inverted
pendulum). The nodal forces of the piece are `Gᵀ f_end + K_arm θ`.
*Verified:* `test_rigid_end_offset_carrying_axial_load` (boundary condition
`M = (H + Pθ) a` in the reference), `test_rigid_end_offset_carrying_the_load`
(buckling), `test_rigid_offset_equals_a_very_stiff_stub_member`.

### 10.8 Recovery; meaning of the reported forces

With the node displacements: end displacements, beam-side displacements from
10.6, end forces `F = K u_b + f₀`; the start state `(d₀, −F₀ − P₀)` is carried
along the piece with the transfer matrix of 10.3/10.4, so **station moments
and deflections are second-order accurate inside the member**. The axial
force at the stations is `N(x)` from statics with the current end force.

Reported `N, Vy, Vz` are components along the **undeformed** member axes
(CONVENTIONS §4 holds: station forces at the ends equal the end forces, and
end shears enter joint equilibrium directly). The shear force normal to the
deformed axis is `Vy − N v′`; in second-order results `dMz/dx = −Vy + N v′`
instead of `−Vy`.

## 11. Linear buckling (`buckling.py`)

### 11.1 Consistent geometric stiffness

For a load factor λ on a reference state with axial forces `N_ref(x)`, the
tangent stiffness of 10.3 is `K_T(λ)`. Its expansion about λ = 0,
`K_T = K + λ K_g + O(λ²)`, defines the linear eigenproblem
`(K + λ K_g) φ = 0`. Because the piece stiffness is the second derivative of
the stationary value of the energy Π of 10.2, and Π depends on N only through
`½ N v′²`, the derivative of the stationary value with respect to N is the
derivative of Π at the unperturbed (first-order) field:

    K_g[a, b] = ∫ N_ref(x) ( v_a′ v_b′ + w_a′ w_b′ ) dx  +  rigid-arm terms of 10.7,

with `v_a, w_a` the first-order displacement fields of the piece for a unit
end displacement a — the fields of the exact first-order element, releases
and shear deformation included (cubic in x, so a 4-point Gauss rule is exact
for `N_ref` up to quadratic). This is the *consistent* geometric stiffness; it
is linear in the reference forces. *Verified:* it equals the numerical
derivative of the exact stiffness
(`test_consistent_geometric_stiffness_is_the_derivative_of_the_exact_stiffness`).

### 11.2 Procedure

1. Reference state: **first-order** solution of the combination's factored
   loads (not multiplied by `analysis.alpha`), with its notional loads, its
   stiffness reduction (τ_b from the reference axial forces) and its set of
   slack one-sided members (section 13).
2. Every member is cut into n pieces (load break points are piece
   boundaries). `K` and `K_g` are assembled, constraints and supports
   eliminated as in the static analysis, and the symmetric generalised
   problem `(−K_g) φ = μ K φ`, `λ = 1/μ`, is solved (K positive definite).
   `μ > 0`: buckling under λ × the loads; `μ < 0`: buckling under the reversed
   loads, reported separately (`negative_factors`).
3. Each eigenpair is taken to the eigenpair of the exact problem
   `K_T(λ) φ = 0` on the same mesh by residual inverse iteration: λ from the
   Rayleigh functional `φᵀ K_T(λ) φ = 0` (secant), then
   `φ ← φ − K_T(σ)⁻¹ K_T(λ) φ` with σ just below the starting value. This
   removes the discretisation error of step 2, which decreases with the fourth
   power of the piece length (cubic fields); for a constant axial force per
   member the factors are then exact for the theory.
4. n = 4, 8, … (at most 32) until the positive factors change by less than
   `BUCKLING_TOL = 1e-7` between two subdivisions; the last change is reported
   as `estimatedRelativeError`. If an eigenpair cannot be polished, its mesh
   value is kept and this criterion decides.

**Stated accuracy:** the reported factors are converged to 1e-7 (relative)
with respect to subdivision, otherwise a warning says how far they still
move. *Verified:* the four classical end conditions, higher modes, spring
supports, releases, a leaning column on a rigid diaphragm, a loaded rigid
offset, a column under its own weight and the sway buckling of a portal frame
agree to 3e-7 with the roots of characteristic determinants built from the
general solution of `EI v⁗ + P v″ = 0` and the boundary conditions of each
case, solved numerically (`tests/test_buckling_columns.py`; no
effective-length factor is used anywhere); shear-deformable columns agree with
the eigenvalue of (5) found by `solve_bvp`.

Modes are normalised to a largest translation of 1 over all mesh nodes (so a
mode in which only members bow has small nodal values; `member_shapes` gives
the member axis at the result stations), with the sign that makes the first
translation component reaching that magnitude positive. Coincident factors
have no unique modes.

## 12. Second-order analysis of a combination (`nonlinear.py`)

Each combination is solved on its own from its factored loads.

1. Loads × `analysis.alpha` (AISC 360-22 C2.1(d)); results are divided by
   alpha at the end.
2. First-order solution (this also yields the mechanism diagnosis in terms of
   the original model). It gives the starting axial forces.
3. **Subdivision** (`refine.py`). Members are cut at the load break points.
   A-priori rule: no piece longer than `1.5 / k`, `k = √(2|N|/EI_min)` with
   the first-order N — not needed for accuracy with a constant N (the piece is
   exact) but it keeps the transfer matrices well conditioned in strong
   tension and far from the pole of 10.3 in compression; verified after
   convergence (`k h ≤ 3` with the final N, otherwise repeated with shorter
   pieces). Members whose axial force **varies** along their length are
   additionally cut into 2, 4, 8, … (at most 32) times as many pieces until
   the results (node displacements, station forces and moments, each relative
   to its largest value) change by less than `analysis.convergenceTolerance`
   between two levels; the change is recorded (`subdivisionCheck`). The
   internal nodes are ordinary DOFs of the solve, so the positive-definiteness
   test covers buckling of a member between its joints.
4. **Iteration.** `K_T(N⁽ᵏ⁻¹⁾) u⁽ᵏ⁾ = F(N⁽ᵏ⁻¹⁾)`: stiffness and fixed-end
   forces with the axial forces (and τ_b) of the previous iterate, then new
   axial forces from the end forces. Converged when the change of the
   displacements (2-norm over the original nodes, rotations × model size,
   relative), of the axial forces (largest change / largest force) and of τ_b
   are all below `analysis.convergenceTolerance`; at most
   `analysis.maxIterations`. Not converged: the result is flagged
   (`converged = False`, code `solve.not-converged`).
5. **Beyond the critical load.** If the tangent stiffness is not positive
   definite (or a piece has buckled between its end connections), the
   combination is refused (`stability.beyond-critical-load`) with a buckling
   factor: the largest multiple ≤ 1 of the axial forces of the last iterate at
   which the tangent stiffness is still positive definite, found by bisection
   with the exact pieces. *Verified:* the factor equals critical load / load
   (`test_load_beyond_the_elastic_critical_load_is_refused_with_the_buckling_factor`).
6. Results are mapped back to the stations of the original members; the
   iteration history, the tolerance and the number of pieces per member are
   recorded in `CombinationResult.solver`.

AISC 360-22 Commentary C2.1 benchmark problems (Cases 1 and 2), every
tabulated value, US and SI: largest difference 0.38 % (moments) and 0.88 %
(deflections) against the 3 % / 5 % acceptance criterion
(`tests/test_second_order_aisc_benchmarks.py`; the section properties are not
printed in the Commentary and are calibrated on its own P = 0 columns — see
the test and section 16).

## 13. Tension-only and compression-only members

**Rule.** A one-sided member is *slack* when the axial force it carries, or
would carry, has the wrong sign. The test quantity is its mean axial force
`EA × elongation / L`, with the elongation from the node displacements
(first-order kinematics, rigid offsets included) — defined whether the member
is active or not. A slack member has zero stiffness; its member loads (e.g.
self weight) still go to its end nodes as the reactions of a simply supported
beam. A member keeps its state while its strain is below 1e-9 of the largest
member strain.

**Iteration** (per combination, no superposition): start with all members
active; solve (first or second order); flip every member whose state
contradicts its strain; repeat until nothing changes. If a set of slack
members repeats, or a mechanism appears on the way, the iteration restarts and
flips one member at a time (the one with the largest wrong force). A second
cycle, or `analysis.maxIterations`, is refused (`solve.active-set-cycling`)
naming the members that keep changing.

**A trial set is not the structure.** In a trial set a tension-only member can
still be active while it carries compression (the first trial has every
member active). It cannot carry that force, so the second-order equilibrium
of the trial set — in particular its stability — says nothing about the
structure: slender tension-only diagonals in compression exceed their own
Euler load, or make the tangent stiffness of the trial set indefinite, long
before the structure with the diagonals that can work is near its critical
load. When the second-order solution of a trial set does not exist and the
first-order solution of the same set has one-sided members with the wrong
sign of force, the set is updated from those first-order violations and the
iteration continues (recorded in `solver.activeSet.history` as a note). The
combination is refused as `stability.beyond-critical-load` only for a set
whose first-order state has no wrong-signed member — an admissible set that is
unstable. The accepted result is always a converged second-order solution
whose own strains satisfy the rule above. *Verified*
(`tests/test_solve_roundoff_scales.py`): `examples/05_braced_frame.json`,
gravity combination in second order (formerly refused with "buckling factor
0.38"; linear buckling factor 7.65): agreement with the tableau reference of
`tests/independent` on the same active set to 2e-12, with the sign conditions
of all six diagonals evaluated from the reference's own results; tier C
(OpenSeesPy 1.3e-6 against its estimate 4.8e-6, PyNite within its estimate).
An enumeration of the 27 sets with at most one working diagonal per storey
finds no other admissible set.

**Slack members and stability.** If the converged set leaves a mechanism, the
combination is refused (`stability.slack-mechanism`, with the mechanism
diagnosis of section 7 and the slack members). The same refusal applies when
the result is not determined: some one-sided members end with exactly zero
force and the structure without them is a mechanism that can move without
loading any of them (a gravity-only load on a pin-jointed frame with
tension-only X-bracing: both diagonals shorten, and the frame can sway freely
between the two positions in which one of them becomes taut). If either sense
of that motion loads one of the idle members, the state is unique and
accepted (e.g. zero load).

## 14. Stiffness reduction and notional loads (AISC 360-22 C2.3, C2.2b)

* **C2.3(a)**: factor 0.80 (`analysis.stiffnessReduction.factor`) on EA, EI,
  GJ and G·As of the listed members (default: all).
* **C2.3(b)**: τ_b additionally on the flexural stiffness only (User Note of
  C2.3(b): `0.8 τ_b EI`, 0.8 × the other stiffnesses). `"computed"`:
  Eq. C2-2a (`αP_r/P_ns ≤ 0.5`: τ_b = 1.0) / Eq. C2-2b
  (`τ_b = 4 (αP_r/P_ns)(1 − αP_r/P_ns)`), `P_ns = F_y A` with the yield stress
  of the member's material. Engine rule where the Specification is silent:
  `αP_r` is the largest compression along the analysis member in the analysis
  itself (which runs at α × the loads), one τ_b per analysis member, updated
  in every iteration; tension or no axial force: 1.0; `αP_r ≥ P_ns` is refused
  (`analysis.tau-b-zero`).
* **C2.2b, Eq. C2-1**: `N_i = coefficient × α × Y_i` (default 0.002) at every
  node along the given horizontal direction; `Y_i` = the downward component of
  the loads applied at the node, member loads lumped to the two end nodes of
  their member by the lever rule. `Y_i` is signed (a net upward load gives a
  notional load against the direction). The notional loads are loads: the
  reactions contain them (User Note of C2.2b(a): fictitious base shear).

*Verified:* `test_stiffness_reduction_c2_3`, `test_notional_loads_c2_2b`,
`test_alpha_multiplies_the_loads_and_divides_the_results`.

## 15. Invariants of a second-order solution

Evaluated on the subdivided model with the terms of the theory:

| check | quantity |
|---|---|
| stiffness symmetry | of the assembled tangent stiffness |
| rigid-body modes | `K_T B = 0` for the three rigid **translations** only (a rigid rotation of a member that carries an axial force does produce forces) |
| positive definiteness | smallest scaled eigenvalue of the free-free tangent stiffness |
| linear solve residual | of the last solve |
| global force equilibrium | applied loads + reactions, components in the undeformed axes |
| global moment equilibrium | applied loads + reactions about the origin at the undeformed positions **minus** the couples `∫ N v′ dx` of every piece and the rigid-arm couples (from (4): the moment of all forces on a piece taken at the undeformed positions equals `∫ N v′ dx`) |
| joint equilibrium | nodal loads + reactions + constraint forces − second-order member end forces (arm terms included), internal nodes included |
| member equilibrium / compatibility | resultants and displacements carried over every piece reproduce its end forces and end displacements |
| member equilibrium on the deformed axis | end moments + end shears × length + load moments − `∫ N v′ dx` for every piece, the integral from the recovered deflection (`N (v_j − v_i)` for constant N, Gauss quadrature otherwise) |
| work = strain energy | `uᵀ K_T u` = (equivalent nodal loads + reactions)·u |

Tolerances and the round-off floor of the force scales as in section 8 (8.1,
with K the tangent stiffness and u including the internal nodes; the check on
the deformed axis uses the moment floor). *Verified:* in
`test_invariants_of_a_second_order_solution_are_evaluated_with_the_second_order_terms`
the moments of loads and reactions on the undeformed geometry differ by
exactly the couples of the member axial forces.

## 16. Scope and limits of the nonlinear analyses

* **Included:** the effect of the member axial force on flexure in both
  principal planes — joint translation (P-Δ) and member curvature (P-δ),
  compression and tension, with shear deformation — and on rigid end offsets.
* **Not included:** the effect of bending moments, shear forces and torque on
  the stiffness, and of the axial force on torsion (no `N r₀²` term, no
  warping). **Lateral-torsional, torsional and flexural-torsional buckling are
  therefore not captured** by the second-order or the buckling analysis.
* Small-displacement theory: no large rotations, no bowing (axial shortening
  due to bending), no change of the axial stiffness with N; the elongation
  used for one-sided members is first order.
* Shear-deformable members follow the Engesser form (10.2). The Haringx form
  would differ by terms of relative order `N/(G As)` × (shear deformation);
  the choice is not prescribed by AISC 360-22.
* Rigid links and diaphragms are linear kinematic constraints: the
  second-order effect of the forces they transmit is not included (a warning
  is attached to the results).
* Rigid offsets not along the member axis: only the axial end force acts on
  the rotating arm.
* The buckling analysis is linear: first-order reference state, axial forces
  proportional to the load factor, slack set and τ_b frozen at the reference
  loads.
* Dense matrices throughout: intended for frames of up to a few thousand
  degrees of freedom after subdivision.

## 17. Independent review of the solver (`tests/independent/`)

"Independent" here means independent code, not an independent person or body:
this review was made within the project that wrote the engine, and nobody
outside the project has reviewed the solver.

A second implementation, written without reading the engine's tests and sharing
no code with it (own reading of the input, own local axes, own load
resolution, file units), is kept as permanent tests:

* `tableau.py` — no stiffness matrix, no fixed-end forces, no elimination of
  slaves: node displacements, member end states, reactions and constraint
  multipliers are the unknowns of the equilibrium, compatibility and
  connection equations written down directly; the member is represented by
  `exp(A L)` of its full 12-state first-order system (second order: with the
  axial force in `A`); buckling factors are the roots of the determinant.
* `test_indep_bvp.py` — single members as multi-interval boundary value
  problems (`scipy.integrate.solve_bvp`); `test_indep_conventions.py` —
  station forces of a determinate tree by cutting (direct equilibrium).
* `oracles.py` — direct translations to OpenSeesPy and PyNite for the models
  each can represent exactly.

Largest relative differences measured (each quantity relative to the largest
value of its kind in the model; node displacements, reactions, constraint
forces, end forces, and forces and displacements at every station):

| analysis | models | tableau | OpenSees | PyNite |
|---|---|---|---|---|
| first order (rolls, partial/full releases, off-axis offsets, links, diaphragms, springs 1e-9…1e18, imposed displacements, trapezoid/partial/projected loads, self weight, shear deformation, Y-up/Z-up, five unit systems) | 27 | 4e-12 (1e-10 with 20 mm stubs next to 30 m spans) | 3e-12 (18 models; 1e-10 on the stub model) | 5e-13 (9 models; 2e-10 on the stub model) |
| second order, constant N per member | 5 (+ one-sided, edge cases) | 5e-13 | 2e-4 (32 elements/member, chord Kg) | 2e-4 (8 sub-members, one pass) |
| second order, N varying along members | 3 | converges to the engine at the reference's O(h²): 3e-8 … 4e-4 at 32 steps | — | — |
| linear buckling (2–3 modes, reversed loads, releases, truss members, diaphragm, offsets) | 7 | 7e-11 (1e-7 with varying N, 8 steps) | — | — |
| tension-/compression-only (first and second order) | 3 | 4e-11 (same slack sets) | — | 4e-14 (same slack sets, first order) |
| single member vs `solve_bvp` (first/second order, shear, spring end) | 6 | < 1e-7 (tolerance of the BVP) | | |

Findings of the review that are limits rather than errors:

* The invariants of section 8/15 close the member integration at x = L; a
  wrong value at an *interior* station would not be seen by them
  (`test_a_wrong_station_value_is_not_seen_by_the_invariants_but_by_the_reference`).
* An injected error must exceed `equilibriumTolerance` × the largest quantity
  of the check to be caught: 1e-4 on reactions, displacements, equivalent
  nodal loads or end forces is caught; a 1e-4 error in one small moment next
  to a large force × model size is not.
* `convergenceTolerance` below about 1e-9 cannot be reached when the axial
  force varies along members (round-off of the subdivision comparison): the
  result is flagged not converged after a long run.
* OpenSees `-jntOffset` with `PDelta` has no second-order term of the rigid
  arm (25 % low at an amplification of 3); the engine's term (10.7) agrees with
  stiff stub members in both programs.

