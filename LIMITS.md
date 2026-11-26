# Limits

This file says what the software does **not** do and what it must not be used for. It is written
for a practising structural engineer. Each limit is followed by the place where the behaviour is
defined (a source file or a document of this repository), so that it can be checked.

**Read this first.**

* The results are the output of software that has passed the tests listed in
  [docs/VERIFICATION.md](docs/VERIFICATION.md). They are not a substitute for an engineer's own
  checking of the model, the loads, the results and the design.
* Nothing here is certified, approved or reviewed by anyone. "Verified" in this repository has one
  meaning only: *the result passed the verification tiers listed in its record, at the tolerances
  stated in that record.* It is a statement that several computations agree. It says nothing about
  whether the model represents a real structure or whether a design is adequate.
* The software was written for a benchmark of structural-modelling agents. It is not a design
  office tool, it has no user base behind it, and it has not been used on a real project.
* When the engine cannot do something it is meant to refuse or to report "not checked", never to
  guess. If you find a case where it silently does something else, that is a defect.

Statements marked with an HTML comment `<!-- RECHECK -->` in the source of this file describe
behaviour that was being corrected while the file was written; they describe the code as found and
must be checked again.

Terms used below, in one line each:

* *First-order analysis*: equilibrium is written on the undeformed structure (ordinary linear
  frame analysis).
* *Second-order analysis*: the axial forces act on the deformed structure, so a compressed column
  that sways or bows gets additional moments (P-Δ = effect of joint sway, P-δ = effect of the
  member's own curvature).
* *Buckling factor*: the multiple of the applied loads at which the elastic structure loses
  stiffness.
* *Mechanism*: a structure, or part of one, that can move without any member resisting.
* *Release*: a member end connection that does not transmit one of the six force components
  (a hinge releases a moment).
* *Rigid offset*: a rigid arm between a node and the start of the flexible member.
* *Rigid link / rigid diaphragm*: a rule that makes nodes move together as a rigid body
  (all six movements, or only the movements in a floor plane).
* *Station*: a point along a member at which forces and displacements are reported.
* *Tier*: one kind of verification check (A to F; section 4).

---

## 1. Analysis

### 1.1 What is analysed

The engine analyses three-dimensional frames of **straight, prismatic, linear elastic** members
with six movements per node, under static loads. Everything else in this section is a limit of
that.

| Not done | What happens instead | Defined in |
|---|---|---|
| Material nonlinearity (yielding, plastic hinges, cracking, nonlinear springs) | not available; the material is elastic at any stress | `docs/THEORY.md` 1, 16 |
| Large displacements, large rotations, snap-through, cable behaviour | not available; see 1.2 | `docs/THEORY.md` 10.1, 16 |
| Dynamics: natural frequencies, response spectrum, time history, seismic analysis | not available | no such analysis type in `src/benchfea/schema/` |
| Plates, shells, walls, slabs, solids | not available; line members only | `docs/ARCHITECTURE.md` "Input format decisions" |
| Tapered, haunched or curved members | not available; a curved member must be approximated by straight pieces, and the error of that is yours to judge | `docs/ARCHITECTURE.md` "Known limits of this layer" |
| Warping torsion (the extra torsional stiffness of an I-section whose flanges cannot warp) | not available; torsion is uniform (St. Venant) torsion with the constant J only | `docs/THEORY.md` 1.1 |
| Temperature loads, shrinkage, lack of fit, prestress, any initial strain | not in the input format; a temperature load is refused by validation (`schema.enum`) | `src/benchfea/schema/` |
| Distributed moments along a member; loads applied on a rigid offset | not in the input format | `docs/ARCHITECTURE.md` "Known limits of this layer" |
| Supports inclined to the global axes | not available; supports and support springs act along the global X, Y, Z axes only | `docs/CONVENTIONS.md` 5 |
| General linear constraints (for example "these two rotations are equal and opposite") | not available; only rigid links and rigid diaphragms, and a slave of one constraint cannot be the master of another | `docs/THEORY.md` 4; codes `constraint.*` in `docs/ARCHITECTURE.md` |
| Elastic foundation along a member, soil | not available; only springs at nodes | `src/benchfea/schema/` |
| Moving loads, influence lines, pattern loading, automatic load combinations of a loading code | not available; you write the load cases and the factored combinations | `src/benchfea/schema/` |
| Construction stages, creep | not available | — |
| Self weight of rigid offsets | not included; self weight acts on the flexible length only | `docs/THEORY.md` 5.1 |

### 1.2 Second-order analysis is small-displacement theory

The second-order analysis includes the effect of member axial force on bending in both planes
(P-Δ and P-δ), in compression and in tension, with or without shear deformation. It is the theory
that AISC 360-22 Section C2.1 calls a second-order elastic analysis. It does **not** include:

* large rotations or the shortening of a member because it bends ("bowing");
* any effect of bending moments, shear forces or torque on the stiffness, or of axial force on
  torsion. **Lateral-torsional buckling, torsional buckling and flexural-torsional buckling are
  therefore not found by the second-order analysis or by the buckling analysis.** A beam that
  would buckle sideways is reported as stable by the analysis; that limit state exists only as a
  member check in the design module (section 3), and only for the shapes listed there.

Defined in: `docs/THEORY.md` 10.1 and 16.

Consequence measured in this project: published second-order results obtained with a
large-displacement program differ from this engine. On the six frames of the Ziemian & Ziemian
data set that were encoded, joint displacements differ by **0.4 % to 2.7 %** (largest on a frame
under gravity load only and on the gabled frame), while the engine agrees with another
small-displacement solver to about 10⁻⁴. Those twelve problems are quarantined, not passed.
Defined in: `verification/problems/ziemian/README.md`.

Other second-order limits:

* **Rigid links and rigid diaphragms carry no second-order effect.** The engine treats them as a
  linear rule between displacements. A real rigid floor that transmits forces while it twists adds
  a further term. On one test model (five columns under a rigid floor, strong twist) the two
  statements differ by 7.6 % in the largest displacement. The engine attaches a warning to such
  results. Do not rely on second-order results of a structure with leaning columns and significant
  floor torsion. Defined in: `docs/CROSSCHECK.md` 3.5; `docs/THEORY.md` 16.
* **Rigid offsets not along the member axis**: only the axial end force acts on the rotating arm;
  terms of the kind moment × rotation are left out. Defined in: `docs/THEORY.md` 10.7.
* **A release or spring in a transverse translation** (a sliding connection) is refused in a
  second-order or buckling analysis (`analysis.release-second-order`): the answer depends on how
  the connection is built, which the model does not say. Defined in:
  `src/benchfea/solve/nonlinear.py`.
* **Shear-deformable members** use one of two classical formulations (the "Engesser" form). The
  other ("Haringx") would give slightly different results; AISC 360-22 does not prescribe either.
  Defined in: `docs/THEORY.md` 10.2.
* **Loads beyond the elastic critical load are refused**, not solved
  (`stability.beyond-critical-load`); the message gives an approximate buckling factor. Defined
  in: `docs/THEORY.md` 12.
* **Convergence tolerance.** A `convergenceTolerance` below about 10⁻⁹ cannot be reached when the
  axial force varies along members (self weight of columns, axial line loads); the run is then slow
  and reported as not converged. Defined in: `docs/THEORY.md` 17.
* **Station forces in second order** (N, Vy, Vz) are components along the *undeformed* member axes.
  The shear normal to the deformed axis is a different number. The design checks use the values as
  reported. Defined in: `docs/THEORY.md` 10.8.

### 1.3 Buckling analysis

* It is a **linear** (eigenvalue) buckling analysis: the axial forces are those of a first-order
  analysis and are assumed to grow in proportion to the load factor. It is not a collapse load.
* It finds flexural (Euler-type) buckling of members and frames only — see 1.2 for what it cannot
  find.
* With tension-only or compression-only members, the set of slack members is frozen at the
  reference loads. With stiffness reduction, the reduced stiffness at the reference loads is used.
* Modes that share one buckling factor are not unique.

Defined in: `docs/THEORY.md` 11 and 16.

### 1.4 Truss-type members and pin-ended braces

A member of type `truss`, `tensionOnly` or `compressionOnly` carries axial force only.

* **In a second-order or buckling analysis such a member is a pin-ended strut with the Iy and Iz of
  its section.** If it is compressed beyond its own Euler load (π²EI/L²) the analysis is refused,
  and in a buckling analysis its own buckling appears as the lowest mode. A placeholder moment of
  inertia on a truss section therefore gives a refusal or a meaningless buckling factor. Example of
  the refusal, from an actual run of `benchfea run`:

  ```
  error [stability.beyond-critical-load] at /combinations/U: combination 'U': the load exceeds the
  elastic critical load ... the lowest buckling factor (0.0566077) is the own flexural buckling of
  truss-type member 'b': a truss-type member is analysed as a pin-ended strut with the Iy, Iz of
  its section (own Euler load pi^2 E I/L^2 = 3705.2 N, axial compression 65436 N ...)
  REFUSED, no results written
  ```

  A model that relies on a slender brace carrying compression beyond its Euler load (as some
  published models do, with a single brace element that cannot buckle) cannot be analysed.
  Defined in: `docs/THEORY.md` 10.6; `verification/problems/mgz/README.md` ("Axial-force member in
  a buckling analysis"); `verification/problems/ziemian/README.md` (frame 7).
* Loads applied along a truss-type member are carried to its ends as by a simply supported beam;
  a torque about its axis is refused (`load.torque-on-truss`). Defined in: `docs/THEORY.md` 1.5.
* A slack tension-only member still delivers its own weight to its end nodes, as a simply supported
  beam, without any second-order effect on that member. Defined in: `docs/THEORY.md` 13.
* A model in which idle one-sided members leave a mechanism (for example gravity load only on a
  pin-jointed panel with tension-only cross bracing) is refused (`stability.slack-mechanism`).
  When several such mechanisms exist the engine refuses rather than decide. Defined in:
  `docs/THEORY.md` 13.
* Tension-only members together with second-order analysis: an earlier build refused a braced
  frame at its own loads that has a buckling factor of 7.65; the build examined for this file
  analyses the same case. <!-- RECHECK --> Treat this combination with care and compare with a
  first-order run. Defined in: `src/benchfea/solve/nonlinear.py`.

### 1.5 Refusals that concern stable or nearly stable models

* **Badly conditioned models are refused as "ill-conditioned"**, even if they are stable. Measured:
  a member about 10⁸ times stiffer than its neighbours is refused; 10⁶ times is accepted. Very
  stiff "dummy" members are therefore not a way to model rigid parts; use rigid links, diaphragms
  or offsets. Stiffness contrasts beyond about 10¹⁰ cannot be told from a mechanism. Defined in:
  `docs/THEORY.md` 6 and 7 (`analysis.pivotTolerance`, default 10⁻¹⁰).
* **A rigid diaphragm whose slave nodes are not in the plane of the master node** is accepted by
  validation with a warning, but the result usually fails the engine's own moment-equilibrium check
  and is then refused. Keep the master in the plane of its slaves. Defined in: `docs/THEORY.md` 4;
  `docs/CROSSCHECK.md` 6.
* **Nodes without members, and the rotations of nodes joined only by truss-type members**, are
  taken out of the solution, reported as zero and listed in a warning. A load on one of them is
  refused as a mechanism. Defined in: `docs/THEORY.md` 6.
* **A load case in which every force is exactly zero** (a statically determinate frame under a
  support settlement moves as a rigid body) was refused by an earlier build because its
  equilibrium checks divided one rounding error by another. The build examined for this file
  analyses that case. <!-- RECHECK --> Defined in: `docs/THEORY.md` 8.1.
* Member ends whose rigid offsets overlap, so that the flexible length is reversed, are accepted
  without a warning. Defined in: `docs/THEORY.md` 17.

### 1.6 What the engine's own checks do and do not catch

Every solution is checked for equilibrium, symmetry and energy balance before it is written
("invariants"); a failure refuses the results or flags them. Limits of these checks:

* They are relative: an error must exceed the tolerance (default 10⁻⁶) times the largest force
  (or force × size of the model, for moments) in the check. A small error in one small quantity is
  not caught.
* The engine's checks close at the member ends. Values at stations between the ends are checked
  only by the verification layer (tier D), not by every solve.

Defined in: `docs/THEORY.md` 8, 15, 17; `src/benchfea/verify/stations.py`.

### 1.7 Size and speed

All matrices are dense and the member loop is plain Python. The engine is meant for frames of up
to a few thousand degrees of freedom. Measured for this file on one processor core
(`OMP_NUM_THREADS=1`), regular three-dimensional moment frames with fixed bases, one combination:

| free degrees of freedom | members | first order | second order |
|---|---|---|---|
| 384 | 160 | 0.1 s | 0.4 s |
| 1 296 | 576 | 1.3 s | 3.2 s |
| 3 072 | 1 408 | 6.5 s | 19 s |
| 6 000 | 2 800 | 33 s, about 2.4 GB of memory | not measured |

* Memory grows with the square and time with the cube of the number of degrees of freedom. There
  is no sparse solver. Ten thousand degrees of freedom are not practical.
* Second-order analysis cuts members into pieces internally where the axial force is large or
  varies; models with heavy columns, self weight and one-sided braces take much longer than the
  table suggests (an 18-member frame with second order, tension-only braces, self weight and three
  buckling modes took about 15 s).
* A buckling analysis adds roughly 10 to 60 s or more on frames of the size the generators make.
* Above 1 200 free degrees of freedom the stability test uses an estimate of the smallest
  eigenvalue instead of the exact value.
* **Set `OMP_NUM_THREADS=1`.** With a multi-threaded linear algebra library the small dense
  matrices were measured two orders of magnitude slower in this project's environment.

Defined in: `src/benchfea/solve/system.py` (`EXACT_EIGEN_LIMIT`); `docs/THEORY.md` 16;
`docs/VERIFICATION.md` 6.

### 1.8 Section properties computed from dimensions are approximations

For sections given by dimensions (`I`, `box`, `pipe`, `rectangle`) the engine computes the
properties itself. The torsion constant of I-shapes and boxes, the warping constant of I-shapes
and every shear area are simple thin-wall formulas. Box sections have sharp corners. Fillets of
rolled shapes are included in the area only if the radius `r` is given. If a tabulated value
matters, give it with `overrides` or use a `properties` section. Defined in:
`src/benchfea/sections.py`.

Shear deformation is included only if `analysis.shearDeformation` is true and the section has a
shear area for that direction; otherwise it is silently absent for that direction (reported in
the list of defaults). Defined in: `docs/THEORY.md` 1.3.

---

## 2. Input format

Defined in: `src/benchfea/schema/` (JSON Schema), `docs/CONVENTIONS.md`,
`docs/ARCHITECTURE.md` "Input format decisions", `docs/FORMAT_VERSIONING.md`.

* **One format, version 1.0.0**, as JSON or the equivalent YAML. There is no import from or export
  to any commercial program, no IFC, no drawing exchange. No graphical interface.
* **Units**: force N, kN, lbf or kip; length mm, m, in or ft. Nothing else (no cm, no tonne). One
  unit system per file; every number in the file is in it, including moduli and unit weights.
* **Vertical axis** Y or Z only.
* **Loads**: nodal forces and moments, point forces and point moments on members, linearly varying
  line forces on members (full or partial length; global, local or projected), self weight,
  imposed support displacements. Nothing else (see 1.1).
* **Members**: `frame`, `truss`, `tensionOnly`, `compressionOnly`; one section and one material
  per member; releases per end and per local direction, fully released or with a spring; rigid end
  offsets.
* **Sections**: `I` (doubly symmetric), `box`, `pipe`, `rectangle`, `properties` (numbers only),
  `standard` (a shape name looked up in a database you supply — see 3.6). No channels, tees,
  angles, built-up or composite sections by dimension. A singly symmetric or unsymmetric section
  can be entered only as `properties`, the local axes are then assumed to be its principal axes, no
  distance between shear centre and centroid is modelled, and the design module will not check it.
* **A YAML file must mean exactly what its JSON form means.** Dates, octal or hexadecimal numbers,
  `1_000`, merge keys, duplicate keys and non-text keys are refused with line and column, because
  YAML would otherwise reinterpret them silently.
* **The format cannot say** whether an I-shape is rolled or welded, whether a box is a
  manufactured hollow section (HSS) or a welded box, what the corner radius of a hollow section
  is, or whether the wall thickness given is nominal or design thickness. The design module
  decides these from the kind of section and says so in every check (section 3.3).
* **Design data never change analysis results**: effective length factors, unbraced lengths, brace
  points, the design method. This is intended and tested; it also means that a brace point given
  as design data does **not** restrain the member in the analysis.
* **Notional loads** have no concept of a storey or a level; see 3.5.
* The number of result stations per member is 2 to 1001 (default 11). Two stations (the ends only)
  are allowed and make several design checks impossible; they are then reported as not checked.
* Additions to the output documents were made without changing the output format version 1.0.0
  (new keys only). Defined in: `docs/FORMAT_VERSIONING.md`.

---

## 3. Design checks

The design module checks steel **members** against ANSI/AISC 360-22. It checks nothing else.

### 3.1 The rule you must know

Every check has a status. A check the software could not make is `notChecked` (with the reason) or
`refused`. **A member with any applicable check that is not checked is `incomplete` — never
`pass`.** The four member results are, from best to worst:

| status | meaning |
|---|---|
| `pass` | every applicable limit state was checked and is satisfied |
| `conditional` | as `pass`, but under a stated condition — most often that the forces come from a first-order analysis, so second-order effects are not in them |
| `incomplete` | at least one applicable limit state was not checked; **the member is not shown to be adequate** |
| `fail` | at least one checked limit state is exceeded |

The command `benchfea design` returns exit code 4 for `fail` and for `incomplete`, 0 for `pass`
and `conditional`. Example from an actual run (a solid rectangular section, which is outside the
scope):

```
member  status         ratio  comb.     governing limit state [clause]
c       incomplete     0.000  U         tension: yielding in the gross section [AISC 360-22 D2(a) ...]
          not checked -- compression: flexural buckling about the major axis (x): notChecked:
          section 'col' (section type 'rectangle') is outside the implemented scope (doubly
          symmetric I-shapes, rectangular HSS / box sections, round HSS / pipe)
```

A ratio shown next to `incomplete` is the largest ratio among the checks that *were* made; it is
not the utilisation of the member. Defined in: `src/benchfea/design/__init__.py`;
`docs/ARCHITECTURE.md` "Contract for `benchfea.design`".

### 3.2 Code and shapes

* **One code: AISC 360-22**, LRFD or ASD. Any other code name gives a single `refused` check. No
  earlier AISC edition, no Eurocode, no other standard.
* **Shapes checked**: doubly symmetric I-shapes, rectangular hollow sections and welded boxes,
  round hollow sections and pipes. **Not checked** (compression, flexure and shear are
  `notChecked`; tension yielding is still checked because it needs only the area): channels, tees,
  angles, double angles, singly symmetric or unsymmetric shapes, solid bars and rectangles,
  sections given by `properties` only.

Defined in: `docs/ARCHITECTURE.md` "Known limits of the design checks".

### 3.3 What is implemented, by section of AISC 360-22

| Chapter | Implemented | Not implemented → the check is `notChecked`, or the quantity is an input |
|---|---|---|
| B | B3.1 (LRFD), B3.2 (ASD) with the printed safety factors; B4.1 classification for the shapes above (Table B4.1a Cases 1, 2, 5, 6, 8, 9; Table B4.1b Cases 10, 11, 13, 15, 17, 19, 20, 21) | B4.2: the 0.93 factor on the wall thickness of hollow sections is never applied — the thickness in the file is used as the design thickness. B4.3 net area: not computed |
| D | D2 yielding and rupture; D1 slenderness as an **advisory** note | D3: the shear lag factor U and the net area An are **inputs**; the defaults are An = Ag and U = 1.0 (reported as defaults), which is the unconservative end. D4–D6 (built-up members, pins, eyebars) |
| E | E3 flexural buckling about both axes; E4 torsional buckling of doubly symmetric members (Eq. E4-2); E7 slender elements; E2 slenderness as an **advisory** note | E4 for singly symmetric and unsymmetric members; bracing offset from the shear centre; E5 single angles; E6 built-up members; buckling stress from a buckling analysis; round sections with D/t ≥ 0.45E/Fy |
| F | F2, F3 (I-shapes, compact web, major axis); F6 (I-shapes, minor axis); F7 (rectangular hollow sections and boxes: yielding, flange and web local buckling, lateral-torsional buckling); F8 (round); Cb from Eq. F1-1 | **F4 and F5: I-shapes with noncompact or slender webs (typical plate girders)**; slender webs of boxes (Eq. F7-7); F9–F12; F13 (holes, proportions); load height (loads applied at the top flange); any shape not listed |
| G | G2.1 (webs of I-shapes without stiffeners), G4, G5, G6 | Transverse stiffeners and tension field action (G2.2, G2.3), G2.4, G3, G7, web openings; no interaction of shear in two directions |
| H | H1.1 and H1.2 (axial force with bending), evaluated along the member | H1.3, H2, **H3 (torsion)**, H4. Torsion: see 3.4 |
| L | deflection of a member relative to the line through its ends, against limits **you give** | AISC 360-22 sets no deflection limit; no drift check; members with a free end are `notChecked` |
| C, App. 7, App. 8 | the analysis options of C2 (second order, stiffness reduction of C2.3, notional loads of Eq. C2-1) are in the analysis; warnings about inconsistent settings | The amplifiers B1 and B2 of Appendix 8 are **never applied**. The effective length method of Appendix 7 is not automated: K factors are inputs and are never computed or changed |
| J, K | — | **No connection design of any kind**: bolts, welds, block shear, bearing, gussets, base plates, column web panel zones, concentrated forces on flanges and webs |
| I | — | No composite members |
| App. 1–6 | — | No inelastic design, ponding, fatigue, fire, evaluation of existing structures, member stability bracing design |

Also absent: **all seismic provisions** (AISC 341 is not implemented and was not opened), member
design for any other material, foundations.

Defined in: `docs/ARCHITECTURE.md` "Contract for `benchfea.design`" and the two "limits"
subsections; module docstrings in `src/benchfea/design/aisc360_22/`.

How the kind of section is decided, since the file cannot say (reported in every check):

* a database W, M, S or HP shape is rolled; its web depth needs the database's design fillet
  dimension — without that column the web checks are `notChecked`;
* `type: I` with a fillet radius is rolled; with radius 0 it is treated as **welded built-up**,
  with a warning;
* `type: box` is a **welded box, not a hollow section**; rectangular hollow sections are reachable
  only through a database shape;
* `type: pipe` is a round hollow section.

Treating an intended hollow section as a welded box is slightly unconservative in shear (about 7 %
more shear area in one worked case). Defined in: `src/benchfea/design/section.py`.

### 3.4 Engine rules that are not AISC text

Where the Specification is silent the module either reports `notChecked` or applies a stated rule.
The rules are listed in [docs/OPEN_QUESTIONS.md](docs/OPEN_QUESTIONS.md); the ones that matter most:

* **Torsion.** H3 is not implemented. A member's torsion is ignored only if a simple estimate of
  its shear stress is at most 1 % of 0.6Fy — a threshold chosen for this engine, not taken from
  AISC. Above it the member is `incomplete`. Warping stresses are never considered.
* **First-order forces.** A design run on first-order results is at best `conditional`: the
  Specification requires second-order effects, and the module does not amplify. Run the strength
  combinations as second order.
* **Lateral-torsional buckling.** Brace points are the `bracing.topFlange` and
  `bracing.bottomFlange` positions only; a brace on the tension flange does not count; torsional
  braces do not count. If only a length `Lb` is given, Cb = 1.0. Cantilevers: only the case with a
  built-in support is checked (warping there is assumed prevented); overhangs are `notChecked`
  unless you give Cb. An earlier build chose the wrong unbraced segment when a moment at a pinned
  end was a rounding error instead of exactly zero; forces below 10⁻⁹ of the member's force scale
  are now treated as zero. <!-- RECHECK -->
* **Combined forces (H1)** use the forces at each point with the *lowest* strengths of the whole
  member. This is on the safe side and can overstate the ratio of a long member with several
  unbraced segments.
* **Compression with brace points** uses the longest segment with the largest compression of the
  member (safe side).
* **Torsional buckling (E4)** is evaluated where the first paragraph of E4 calls for it or leaves
  it open, and not for a rolled I-shape or hollow section whose torsional unbraced length does not
  exceed the lateral ones. An earlier build evaluated it for every member. <!-- RECHECK -->
* **Slenderness** (KL/r ≤ 200, L/r ≤ 300) is advisory, as in the Specification; exceeding it never
  fails a member.
* **Values between stations.** Moments are interpolated between stations (exact for first-order
  results with linearly varying loads, approximate in second order). Shear and axial force are
  taken at stations only. Deflections are taken at stations only. With few stations these are
  coarse; with two stations several checks are `notChecked`.

Defined in: `docs/ARCHITECTURE.md` (both "limits" subsections); `src/benchfea/design/demands.py`;
`src/benchfea/design/aisc360_22/f_flexure.py`, `h_combined.py`, `e_compression.py`.

### 3.5 Direct analysis method: what you must do yourself

The engine provides the ingredients of AISC 360-22 Chapter C, but does not decide for you:

* You choose second-order analysis, stiffness reduction and notional loads per combination. The
  design module warns when a strength combination lacks them; it does not add them.
* **Notional loads are computed per node** from the downward load applied at that node (member
  loads are shared between the member's two end nodes). There is no storey or level. A node with a
  net upward load gets a notional load in the opposite direction. Whether this is what the
  Specification means for your structure is for you to judge; the Specification defines notional
  loads by level. Defined in: `docs/THEORY.md` 14.
* **Initial imperfections are not modelled** except through notional loads. Member out-of-
  straightness is not modelled.
* The stiffness reduction factor τ_b, when computed, uses the gross area and one value per analysis
  member. Defined in: `docs/THEORY.md` 14.
* For ASD you set `analysis.alpha` to 1.6 yourself (a warning is given otherwise).

### 3.6 The shapes database is not included

Sections given by a standard shape name (`type: standard`) need the AISC Shapes Database. It is not
redistributable and **is not in this repository**; you supply your own copy. The reader's column
names and units are assumptions that could not be checked against the database (it could not be
opened for this project); they can be overridden. Until you have confirmed them against your copy,
compare the properties the engine reports with the printed table. Defined in:
`src/benchfea/shapes.py`.

No AISC document is in the repository either; the equations are implemented from the Specification
and cited by section and equation number.

---

## 4. Verification

### 4.1 What "verified" means

A reference solution is `verified` when it passed every tier that applies to it:

| tier | check |
|---|---|
| A | a single straight member line recomputed by direct integration of the member equations |
| B | targets of a published problem reproduced |
| C | the same model solved by two independent programs (OpenSeesPy and PyNite), which must agree with the engine and with each other |
| D | equilibrium, energy, superposition, reciprocity; station values recomputed from statics |
| E | the same model built a different way (members split, offsets replaced by stiff stubs, …) gives the same answer; design data change nothing in the analysis |
| F | for design results: every limit state used has both a hand calculation and a reproduced value of the AISC Design Examples in the test suite |

One tier not passed — for whatever reason, including "could not be run" — means `quarantined`.
There is no majority vote. Defined in: `docs/VERIFICATION.md` 1–2.

### 4.2 What "verified" does not mean

* Not reviewed, approved or certified by any person or body.
* Not that the model is a good model of anything.
* Not exact: linear results agree to about 10⁻⁶; second-order results agree only to the
  discretisation error of the other programs, 10⁻⁵ to 10⁻³. <!-- RECHECK -->
* Not that the capacity of an unused design check is right: tier F covers only limit states with a
  nonzero demand.
* A model with no loads verifies trivially.

### 4.3 What cannot currently be verified

The following results can be computed by the engine but **cannot reach the status `verified`**,
because fewer than two independent programs can represent them. They are tested in the engine's
own test suite, but that is not the same thing.

| Feature | Why | Defined in |
|---|---|---|
| Shear deformation | PyNite has none; OpenSees only in first order → one program | `docs/CROSSCHECK.md` 4, 5 |
| Partial releases in axial, shear or torsion (springs) | PyNite has none | same |
| **Any buckling factor** | neither translator computes buckling; buckling is checked only on single-member problems derived in this repository and on a few published problems | `docs/VERIFICATION.md` 7 |
| Rigid offsets in a second-order analysis | no error bound for either program | `docs/CROSSCHECK.md` 4 |
| Rigid links and diaphragms in a second-order analysis | the two statements of the problem differ (1.2) | `docs/CROSSCHECK.md` 3.5 |
| Second order with members close to their buckling load | the other programs' error bound exceeds the limit of 5·10⁻³ | `docs/VERIFICATION.md` 7 |
| Second order with shear deformation | no program | `docs/CROSSCHECK.md` 4 |
| A mechanism that exists only because tension-only members are slack | no independent computation of which members are slack | `docs/VERIFICATION.md` 4 |
| Most design results | tier F coverage is incomplete: at the time of writing 61 of 80 clause items are fully covered <!-- RECHECK -->; not covered are, among others, torsional buckling (E4), the computed r_ts of Eq. F2-7, box sections, slender flanges in flexure, round-section local buckling, web shear with reduced strength, deflection limits and the torsion screening rule. Most columns and unbraced I-beams use one of these | `docs/VERIFICATION.md` 2 (tier F); `src/benchfea/verify/design_coverage.json` |
| Anything, when OpenSeesPy and PyNite are not installed | tier C cannot run | `docs/VERIFICATION.md` 7 |

Second-order frames in general: with the rule used when the dataset was built, most second-order
models missed tier C by a small margin and were quarantined, so **no stable second-order result is
served as a reference in the dataset**. The rule has since been replaced by an error bound from
three refinement levels. <!-- RECHECK -->

### 4.4 The published problems

173 problems are in `verification/problems/`. At the time of writing 109 are active and 64 are
quarantined with a written diagnosis. <!-- RECHECK --> Current numbers are in
`docs/VERIFICATION_REPORT.md`. What a sceptical reader should know:

* **"Passed" is not always the reproduction of a printed number.** For the vendor problems (STAAD,
  ANSYS) and five textbook problems only the *statement* is published or usable; the target values
  were computed in this repository by separate small programs. That is two computations by one
  project, not an outside reference. The vendors' printed values that do not reproduce are kept as
  separate quarantined problems. Defined in: `docs/VERIFICATION.md` 2 (tier B).
* **The AISC Commentary benchmark problems** (32 problems) test the second-order amplification
  only. The Commentary gives no section properties, so the moment of inertia and the shear
  flexibility were calibrated on the Commentary's own zero-axial-load values; six targets are
  reproduced by construction. The sixteen "with shear deformation" problems do not test shear
  deformation: the Commentary's own tolerance (3 % and 5 %) is wider than the whole shear effect.
  Defined in: `verification/problems/aisc-commentary/README.md`.
* **The AISC Design Examples frame (10 problems) is quarantined**: the examples do not state the
  shear areas or whether the software behind their figures included shear deformation, and the
  model that can be built from the stated data is about 2 % stiffer. Defined in:
  `verification/problems/aisc-de/README.md`.
* **Eighteen textbook problems are quarantined only because of rounding in the book's hand
  calculations**; the engine agrees with the exact solution of the book's own equations. Defined
  in: `verification/problems/mgz/README.md`.
* **Published second-order results from large-displacement programs are not reproduced** (1.2).
* Sketches on the vendor pages could not be read; a page was used only where its text or input
  listing fixes the supports and loads. Two statements rest on that reading.
* Plane problems from the sources are modelled in three dimensions with out-of-plane movements
  held and a large out-of-plane stiffness, because the engine has no continuous out-of-plane
  restraint of a member.

### 4.5 Common-mode risks: where one mistake would pass every check

The checks are only as independent as their parts. These are shared:

* **One author.** The engine, the member-equation oracle of tier A, the translators to OpenSeesPy
  and PyNite, the separate reference programs for the vendor problems and the "tableau" solver
  used for extra evidence were all written within this one project. OpenSeesPy and PyNite
  themselves are the only outside computations.
* **Shared interpretation of the input.** Both translators, the station checks and the engine take
  the model as resolved by `benchfea.model` (local axes, load directions, projection of loads,
  unit conversion, self weight). A mistake there would be common to the engine, both outside
  solvers and the station checks, and no tier would see it. Only the tableau solver reads the raw
  file itself, and it is **not** a required tier. Defined in: `docs/VERIFICATION.md` 7;
  `src/benchfea/crosscheck/prepare.py`.
* **Shared fixed-end-force code.** For members with releases, spring connections or load types the
  outside programs lack, both translators use one routine of this project to turn member loads
  into node loads. It is tested against the oracle, but an error in it would be common to both
  outside routes. Defined in: `docs/CROSSCHECK.md` 5 (`src/benchfea/crosscheck/fef.py`).
* **Substitutions in tier C.** Notional loads, the computed τ_b and the set of slack one-sided
  members are taken from the engine (after an independent recomputation of the first two) and
  given to the outside programs. Each use is written in the record. Defined in:
  `docs/VERIFICATION.md` 2 (tier C).
* **The independent stiffness matrix** used to verify a mechanism diagnosis re-implements the
  engine's rule for degrees of freedom left out of the solution.
* **Tier F** rests on a naming convention of the tests and on the tests comparing with numbers
  worked by hand; the registry cannot prove that the hand numbers are right.
* **Design Examples inputs.** Where an example does not print a section property, the test
  supplies a stated placeholder or takes the example's tabulated strength as an input; those cases
  are partial reproductions.
* The PyNite route drives internal functions of PyNite 3.2.1 and removes two terms of its
  geometric stiffness so that it solves the same problem; it is PyNite's matrices but not PyNite's
  own procedure. Defined in: `docs/CROSSCHECK.md` 3.3, 5.

Frame3DD and MASTAN2, which the project specification named as further independent solvers, could
not be run. Defined in: `docs/REFERENCES.md`.

---

## 5. Benchmark

Defined in: `docs/BENCHMARK.md` 7–9, `data/DATA_CARD.md`, `docs/DIFFICULTY.md` 5,
`docs/CURVEBALLS.md`.

* **Only verified references are served, so the dataset is biased towards what can be verified.**
  At the time of the last build: 534 tasks from 192 verified cells of 528. <!-- RECHECK -->
  * All stable reference solutions are **first-order** analyses. Second-order models appear only
    as mechanism diagnoses.
  * No task asks for a buckling factor.
  * Models with shear deformation, partial releases and similar features (4.3) are absent.
  * Eight of the 33 curve balls have no served task at all: the five design-variable curve balls
    and three of the constraint curve balls. Several others have one or two served cells.
* **Design tasks are nearly absent**: 2 design tasks, 5 prediction tasks that ask for a design
  ratio and 4 repair tasks with a design criterion, all on continuous beams. A benchmark score
  says almost nothing about an agent's ability to design or check steel members.
* **The models are generated, not real.** Eight families of regular frames with a generic
  catalogue of round-dimension sections (no shape table), loads and combination factors that are
  plausible numbers and not a loading code. Stability is established for one strength combination
  only. The difficulty level is a property of the model, not of the question asked about it.
* **One change per model.** A curve ball changes one thing; real modelling errors come in groups.
* **Task types**: in tool-allowed mode a prediction task can be answered by running the analysis;
  it tests tool use. In no-solver mode the tolerances are wide (5 % full credit, 25 % none).
  A diagnosis task compares a complete description in words with the file. A repair is judged on
  named results, not on the edit, and a correct repair that the pipeline cannot verify scores 0.
  The reference of a design task is the best design found by a stated search, not a proven
  optimum.
* **The held-out split is a convention, not a secret.** Generators, curve balls and the split rule
  are public code; anyone can generate the held-out instances. The environment only refuses to
  serve them unless asked.
* **Held-out overlap.** In the stored dataset, 12 held-out instances (34 tasks) are built on base
  models that train tasks also use, and for 10 repair tasks the reference answer is also a train
  reference answer. The seed rule has since been changed so that base models are not shared, but
  the stored dataset is that of the earlier rule until it is rebuilt. <!-- RECHECK -->
* **Small sample.** Two seeds per (family, curve ball) cell; each train pair has one seed.
* **Tied to this code.** The dataset is reproducible only with the generator, curve-ball and
  verification code of the same version; a change in any of them changes the tasks.
* **No licence has been granted** for the code or the data: the repository has no licence file.
* The benchmark measures agreement with this engine's verified results. An agent that disagrees
  with the engine for a good reason (section 1.2, for example) is scored as wrong.

---

## 6. Things that look supported but are not

These exist in the input format, the output or the code, and do less than their name suggests.

| Item | What it actually does | Defined in |
|---|---|---|
| Design variables `role`, `Cmx`, `Cmy`, `sideswayInhibited` | accepted and reported; **used by no check** | `src/benchfea/designparams.py`; no use in `src/benchfea/design/` |
| `bracing.torsion` positions | used for the torsional buckling length of columns only; **not** counted as braces against lateral-torsional buckling | `docs/ARCHITECTURE.md` (flexure limits) |
| `maxSlendernessCompression`, `maxSlendernessTension` | set the limit of an **advisory** note; exceeding it never changes the member status | `src/benchfea/design/aisc360_22/d_tension.py`, `e_compression.py` |
| Appendix 8 functions (B1, Cm, Pe1) | exist as plain functions, tested against one example; **never applied** to any member | `src/benchfea/design/aisc360_22/h_combined.py` |
| `type: standard` sections | need a database that is not included, read with unconfirmed column names (3.6) | `src/benchfea/shapes.py` |
| `type: box` | checked as a welded box, **not** as a hollow structural section | `src/benchfea/design/section.py` |
| `type: rectangle` and `type: properties` | analysed; **no design check** beyond tension yielding | 3.2 |
| A brace point in the `design` block | changes design capacities only; **does not brace the analysis model** | `docs/ARCHITECTURE.md` "Input format decisions" |
| Buckling analysis | computed and reported; can never be `verified` (4.3); says nothing about lateral-torsional buckling (1.2) | `docs/THEORY.md` 11 |
| Rigid diaphragm in second order | computed with a warning; a term of the mechanics is left out (1.2) | `docs/CROSSCHECK.md` 3.5 |
| Diaphragm with slaves off the master's plane | accepted with a warning by validation, then usually refused after the solve | 1.5 |
| `resultStations: 2` | allowed; several design checks then cannot be made | 3.4 |
| Frame3DD in the capability table of the cross-check | listed as `not-run`: it was never executed | `src/benchfea/crosscheck/capabilities.py` |
| `--extra-evidence` (tableau solver) | recorded next to tier C; **never changes a status** | `docs/VERIFICATION.md` 2 |
| Tier B for an ordinary model | applies only to registered library problems; for your own model it is "not applicable", not "passed" | `docs/VERIFICATION.md` 2 |
| "Passed" problems with basis `recomputed-from-published-statement` | targets computed in this repository, not printed in the source (4.4) | `docs/VERIFICATION.md` 2 |
| `verification.route` in a curve-ball record | an expectation derived from the model's features, not the result of a verification run | `docs/CURVEBALLS.md` |
| `benchfea verify` on a design | tier F tells whether the limit states are *covered by tests*; it does not re-derive your member's capacity | `docs/VERIFICATION.md` 2 (tier F) |

---

## 7. Do not use this software for

* the design or assessment of a real structure without independent checking by a qualified
  engineer using other means;
* anything in section 1.1 (nonlinear material, dynamics and seismic analysis, plates, cables,
  temperature), by approximation or workaround;
* stability of beams against lateral-torsional buckling by analysis; plate girders; members with
  significant torsion; connections; composite construction; seismic design;
* slender structures where large displacements matter (cable-like systems, shallow arches,
  snap-through);
* models larger than a few thousand degrees of freedom;
* as an authority on what AISC 360-22 requires. Where the Specification is silent the module
  follows the readings recorded in [docs/OPEN_QUESTIONS.md](docs/OPEN_QUESTIONS.md), which are
  one project's readings and are waiting for decisions.
