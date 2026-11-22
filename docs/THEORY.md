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
7 stability diagnosis · 8 invariants · 9 hooks for nonlinear analysis.

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
   as a translation, a rotation about an axis through a point, or both.
2. *Single-node mechanisms*: null vectors confined to one node (eigenvalue 1
   of `Q_aᵀQ_a` for an orthonormal null basis Q restricted to the node's
   rows), again in row echelon form, with the cause: all connected member ends
   released, only truss members, or no stiffness.
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

## 9. Hooks for nonlinear analysis

`ANALYSIS_TYPES[type](analysis, combination)` solves one combination; only
`firstOrder` is registered. A second-order driver can reuse
`Structure.assemble(matrices)` for the geometric stiffness (per-member 12×12
matrices in nodal global DOFs, `Gᵀ k_g G`), `Structure.factorize` (stability
diagnosis included), `Structure.solve`, `LinearAnalysis.recover`, and
`build_element(model, member, factors, active)` for reduced stiffness and for
members that are switched off.
