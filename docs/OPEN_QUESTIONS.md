# Open questions

Every place where an opened reference or AISC 360-22 was silent or ambiguous, and every policy decision that
belongs to the project owner, collected from the build and review reports of the project, the reading notes of
the references and the READMEs of the verification problems, with duplicates merged.

The rule of the project is: **do not guess.** So for each question the software does one of three things, and
the entry says which:

* reports the affected check as *not checked* (the member is then `incomplete`), or refuses the model;
* quarantines or excludes the affected verification problem, with a written diagnosis;
* follows one stated reading, which is written in the output or in the document named — these are the entries
  to look at first, because a different answer changes results.

Nothing here has been decided by anyone outside the project, and none of the readings is endorsed by the
publisher of any reference.

Each entry gives: the question; **Where** it arises; what the software does **Now**; and what an answer would
**Change**. Source paths that begin with `design/`, `solve/`, `verify/`, `crosscheck/`, `bench/`, `oracle/` or
`schema/`, and bare module names such as `shapes.py`, are under `src/benchfea/`; document names are under `docs/`.
An HTML comment `<!-- RECHECK -->` in the source of this file marks a statement about behaviour that was being
changed while the file was written; it describes the code as found and must be checked again.

119 questions in 6 groups.

## Decisions needed

The most consequential first. Each line points to the full entry.

1. **Q63** — The second-order rule of tier C. It decides whether second-order references exist and how loose "verified" is for them. A new rule was put in while this file was written; it needs confirmation.
2. **Q64** — Printed-precision policy against hand-rounded textbook values: 18 McGuire/Gallagher/Ziemian problems are quarantined only for this.
3. **Q67** — Tier F granularity and coverage: as it stands almost no generated model can have a verified design reference, so the benchmark has almost no design content.
4. **Q69** — Whether problems recomputed from a published statement are tier A or tier B (30 problems).
5. **Q70** — Whether the AISC Commentary problems, whose section properties are calibrated on the Commentary's own values, count.
6. **Q71** — Whether the tableau solver becomes a required tier — the only protection against an input-interpretation error shared by all other routes.
7. **Q1** — Slenderness limits: advisory (as now) or pass/fail.
8. **Q2** — Notional loads without a storey concept: confirm or replace the per-node rule.
9. **Q43** — Rigid links and diaphragms in second order: include the missing term (7.6 % on one test model) or keep the linear constraint.
10. **Q6** — Design on first-order results: `conditional`, `incomplete`, or amplification by B1/B2.
11. **Q34** — Torsion: keep the engine's 1 % screening rule or implement H3 for hollow sections.
12. **Q46** — Truss-type members refusing the analysis at their own Euler load.
13. **Q87** — Large-displacement published results (Ziemian): accept the engine's small-displacement theory with these 12 problems quarantined, or ask for a large-displacement option.
14. **Q42** — The shear formulation in second order (Engesser, as implemented, or Haringx).
15. **Q72** — Buckling factors: publishable without a second independent solver, or not at all.
16. **Q101** — Rebuild the dataset under the new split rule.
17. **Q107** — A licence for the code and the data.

## A. AISC 360-22 interpretation

Places where the Specification is silent or can be read in more than one way for an automated member check. Section and equation numbers are those of ANSI/AISC 360-22.

**Q1. Slenderness limits: advisory or pass/fail?**
The limits L/r ≤ 300 (tension) and Lc/r ≤ 200 (compression) appear only in User Notes. Should exceeding them fail a member, warn, or be ignored — and with which length and which radius of gyration?
- **Where:** D1 and E2 User Notes; `design/aisc360_22/d_tension.py`, `e_compression.py`.
- **Now:** Reported as an advisory check with the limit from `maxSlendernessTension` / `maxSlendernessCompression` (defaults 300 / 200); it never changes the member status. Tension uses the design member length and the smaller radius; compression uses KL/r per axis.
- **Change:** With "pass/fail" slender bracing members would become `fail`; generated models and design tasks would change; tier F would need coverage of the limit.

**Q2. Notional loads without a storey concept.**
Eq. C2-1 defines the notional load at a *level* from the gravity load at that level. A general frame model has no levels. What is Y_i at a node, should non-gravity load cases be excluded, is a net uplift clamped to zero, in which direction(s) do the loads act in 3D, and are they applied at supported nodes?
- **Where:** C2.2b, Eq. C2-1; `solve/nonlinear.py` (`notional_forces`); THEORY.md 14.
- **Now:** Per node: coefficient × α × the signed net downward component of *all* loads of the combination applied at that node; member loads go to the member's two end nodes by the lever rule (rigid offsets ignored in the lever arms); an uplift node gets a load against the direction; one user-given horizontal direction per combination; applied at supported nodes too, so the reactions contain the notional base shear. The input format cannot mark a load case as gravity.
- **Change:** Any other reading changes every analysis result of a combination with notional loads, and the tier-C recomputation of the notional loads.

**Q3. τ_b: which axial force, and on which stiffness?**
C2.3(b) makes τ_b depend on the required axial strength, which depends on the analysis. Is P_r taken per analysis member, per design member or per piece; the maximum along the member; from the iterated second-order analysis or from a first-order one? The Commentary writes the reduced axial stiffness with τ_b, the Specification puts τ_b on the flexural stiffness only. For slender-element sections the Specification uses the effective area in P_ns.
- **Where:** C2.3(b), Eq. C2-2a/b and its Commentary; THEORY.md 14.
- **Now:** One τ_b per analysis member from its largest compression in the analysis itself (at α × loads), updated every iteration; τ_b on EI only (the Specification's wording), 0.8 on the other stiffnesses; P_ns = F_y × gross area with the material's F_y; refused when αP_r ≥ P_ns.
- **Change:** Results of every combination with `tau_b: computed`. Using the effective area needs the design module's section classification inside the analysis, which the architecture forbids today.

**Q4. Which members get the 0.8 stiffness reduction?**
C2.3(a) requires the reduction for members that contribute to the stability of the structure and permits it for all. Who decides?
- **Where:** C2.3(a); `analysis.stiffnessReduction.members`.
- **Now:** The user lists the members; the default is all members.
- **Change:** Nothing in the code; a benchmark rule would be needed if generated models are to use partial lists.

**Q5. Buckling analysis of a combination with stiffness reduction.**
Should a linear buckling analysis use the nominal stiffness always, or the reduced stiffness of the combination it refers to?
- **Where:** C2.3; `solve/buckling.py`; THEORY.md 11, 16.
- **Now:** Reduced stiffness (0.8 and τ_b) frozen at the reference loads; slack set of one-sided members frozen likewise.
- **Change:** Buckling factors of such combinations (by the factor of the reduction).

**Q6. Design on first-order results: `conditional`, `incomplete`, or amplified?**
The Specification requires second-order effects in the required strengths. When a strength combination was analysed first order, is `conditional` the right member status, or should members in compression be `incomplete`? Should the Appendix 8 amplifier B1 be offered as an opt-in?
- **Where:** C1, Appendix 8; `design/__init__.py` (`analysis_basis`).
- **Now:** Status at best `conditional`, with warnings; B1 and B2 are never applied; the design variables `Cmx`, `Cmy`, `sideswayInhibited`, `role` exist in the schema and are used by no check.
- **Change:** The status of every first-order design run; whether those four variables are implemented or removed from the format.

**Q7. Appendix 8 as a reference method.**
B1 needs C_m (undefined for general transverse loading), is undefined when αP_r ≥ P_e1, and applies to an unspecified moment along the member; B2 needs storeys, storey shear and the moment-frame share. Are Appendix 8 results acceptable as benchmark references at all?
- **Where:** Appendix 8, Eq. A-8-1 to A-8-8.
- **Now:** Not applied. Eq. A-8-3, A-8-4 and A-8-5 exist as plain functions, reproduced against one Design Example.
- **Change:** Whether B1/B2 are implemented; needs a storey concept in the input format (see the input-format group).

**Q8. Effective length method (Appendix 7).**
Appendix 7 requires K from a sidesway buckling analysis without prescribing a method; the Commentary's methods give different K. In Eq. A-7-3 the minimum additional lateral load is printed without α while the first term has it: for ASD is the minimum 0.0042 Y_i or 0.0042 α Y_i?
- **Where:** Appendix 7, 7.2 and 7.3, Eq. A-7-3; `verification/problems/aisc-de/` (Example C.1C).
- **Now:** K factors are user inputs and are never computed or changed; warnings are given for K < 1 and for K = 1 without the direct-analysis settings. The first-order method of 7.3 is not implemented; the Example C.1C problems enter the lateral load as the example prints it.
- **Change:** Whether K is ever computed; the loads of the (quarantined) C.1C problems under ASD.

**Q9. ASD: the 1.6 load multiplier.**
C2.1(d) analyses ASD combinations at 1.6 × the loads and divides the results by 1.6. Is the Commentary's alternative acceptable instead, and does α also multiply imposed support displacements?
- **Where:** C2.1(d); `analysis.alpha`; CROSSCHECK.md (prepare).
- **Now:** α multiplies every load including imposed displacements, results are divided by α; the user sets α; the design module warns for ASD with α ≠ 1.6.
- **Change:** Second-order ASD results with imposed displacements.

**Q10. Direct modelling of imperfections.**
C2.2a gives no magnitude or pattern for modelled imperfections (1/500 appears only in C2.2b and User Notes), and member out-of-straightness is not quantified.
- **Where:** C2.2a.
- **Now:** Not modelled. Only notional loads are available; a user can shift node coordinates by hand.
- **Change:** A new input feature; no present result.

**Q11. E4: when does torsional buckling apply to doubly symmetric members?**
E4 applies "when the torsional unbraced length exceeds the lateral unbraced length" without naming which lateral length, and names "certain" built-up members without saying which. The User Note table of E1 lists flexural buckling only for hollow sections, but the text of E4 does not exclude closed sections.
- **Where:** E4 first paragraph, Table User Note E1.1, E3 User Note; `design/aisc360_22/e_compression.py` (`e4_scope`).
- **Now:** Evaluated for built-up I-shapes and boxes always; for rolled I-shapes and hollow sections only when L_z exceeds the smaller of L_x, L_y or K_z L_z exceeds a lateral effective length; skipped otherwise, with the reason recorded. An earlier build evaluated it for every member, which made it govern by up to 0.35 % for some column proportions and blocked tier F for every I-shaped column. <!-- RECHECK -->
- **Change:** The governing compression limit state and P_c in H1 for short columns; how many generated models can have a verified design reference.

**Q12. E3: which branch test is normative?**
E3 states the inelastic/elastic limit both as Lc/r ≤ 4.71√(E/Fy) and as Fy/Fe ≤ 2.25; the two differ slightly (4.71²/π² = 2.2477).
- **Where:** E3, Eq. E3-2 / E3-3.
- **Now:** Lc/r for flexural buckling, Fy/Fe for torsional buckling; a note is added when the two tests disagree (effect below 0.05 %).
- **Change:** Nothing above 0.05 %.

**Q13. E7: effective area details.**
(a) Is the effective area computed per buckling mode with that mode's F_n, or once with the governing F_n? (b) Table E7.1 prints rounded c2 values that differ from Eq. E7-4. (c) No cap b_e ≤ b is stated. (d) Eq. E7-7 gives A_e > A_g just above D/t = 0.11E/Fy. (e) D/t ≥ 0.45E/Fy has no provision.
- **Where:** E7.1, Table E7.1, Eq. E7-3, E7-4, E7-7.
- **Now:** (a) per limit state with its own F_n; (b) the rounded table values, as the Design Examples use; (c) b_e capped at b; (d) A_e capped at A_g; (e) `notChecked`.
- **Change:** Capacities of slender-element columns by 0.1–1 %.

**Q14. Safety factors: printed 1.67 or 1.5/φ?**
The printed Ω values (1.67, 2.00, 1.50) differ from 1.5/φ by up to 0.2 %.
- **Where:** B3.2 and each chapter; Commentary to B3.
- **Now:** The printed values.
- **Change:** Every ASD capacity by up to 0.2 %.

**Q15. E and G: no value is enforced, and the US and SI values are not exact conversions.**
The Specification prints 29,000 ksi and 200 000 MPa (11,200 ksi and 77 200 MPa). Should the design module insist on them, and must results agree between unit systems?
- **Where:** Symbol lists of E3, E4, F2; `designparams.py`.
- **Now:** E and G come from the material (or from design variables that override them for the checks only); nothing is enforced. A warning exists when G and Poisson's ratio are inconsistent.
- **Change:** Whether a model with another E is `refused` by the design module.

**Q16. Classification exactly at a limit.**
No tolerance is given for a width-to-thickness ratio equal to λ_p or λ_r in floating point.
- **Where:** B4.1.
- **Now:** Equal to the limit is the favourable class; no tolerance.
- **Change:** The class of sections built exactly on a limit (the generators' catalogue avoids them only by chance).

**Q17. Web depth h of rolled I-shapes.**
h is the clear distance between flanges less the fillets; the Specification does not say which tabulated fillet dimension to use, nor what to do when only d, bf, tf, tw are known.
- **Where:** B4.1b; `design/section.py`.
- **Now:** Database shape: h = d − 2 k_des, and web checks are `notChecked` without that column. `type: I` with radius r: h = d − 2t_f − 2r. With r = 0 the shape is treated as welded (h = d − 2t_f) with a warning.
- **Change:** Web classification, shear and slender-web compression of rolled shapes.

**Q18. Hollow sections: design thickness, corner radius, flat width.**
The design wall thickness (0.93 t or t) depends on the production standard, which is not a geometric property; the corner radius is "if known"; tabulated b/t ratios may or may not take precedence. Round-section provisions are permitted for one pipe grade only by a User Note.
- **Where:** B4.1b(d), B4.2.
- **Now:** The thickness in the file is used as the design thickness (0.93 never applied, stated in the output); flat width = outside − 3t; tabulated ratios are used when a section is built through the Python interface; any `pipe` is a round hollow section.
- **Change:** Classification and capacities of hollow sections; needs a format field (input-format group).

**Q19. Net area and shear lag factor.**
A frame model has no connection data. Should rupture be checked with defaults, require input, or be skipped?
- **Where:** D2(b), D3, Table D3.1, B4.3.
- **Now:** A_n = A_g and U = 1.0 by default, reported as defaults (the unconservative end); Table D3.1 and B4.3b are not evaluated.
- **Change:** Tension rupture capacities; whether a member without A_n and U is `incomplete`.

**Q20. G2.1(b) for a rolled shape that fails G2.1(a).**
The web depth for G2.1(b) is worded for built-up sections only.
- **Where:** G2.1.
- **Now:** The h of G2.1(a) is reused, with a note.
- **Change:** Shear capacity of rolled shapes with slender webs (rare).

**Q21. Shear in the minor direction; shear in two directions.**
G4 and G6 imply but do not state the rule for the other wall pair / the two flanges; no interaction is given for simultaneous shear in two directions; for round sections the use of the resultant shear is natural but not stated.
- **Where:** G4, G5, G6.
- **Now:** G4 with the other wall pair; G6 as twice Eq. G6-1; each direction checked separately; round sections on the resultant shear.
- **Change:** Shear ratios under biaxial shear.

**Q22. G5: L_v when the shear never vanishes.**
L_v is the distance from maximum to zero shear; it is undefined for constant shear and ambiguous with several zeros.
- **Where:** G5, Eq. G5-2a.
- **Now:** Found by a piecewise-linear zero crossing of the station shears; when the shear never vanishes Eq. G5-2a is not credited (lower-bound capacity).
- **Change:** Shear capacity of round members with constant shear.

**Q23. Which brace points count against lateral-torsional buckling?**
L_b is defined by braces against lateral displacement of the compression flange *or* twist. Do `bracing.torsion` points count? Does a brace on the tension flange?
- **Where:** F2.2 definition of L_b; `design/aisc360_22/f_flexure.py` (`ltb_segments`).
- **Now:** Only `bracing.topFlange` / `bottomFlange` on the flange that is in compression; member ends are brace points; torsional braces and tension-flange braces are not counted.
- **Change:** LTB capacities of members with torsional bracing; the design-variable curve balls.

**Q24. C_b where the Specification gives no rule.**
(a) L_b given as a length without brace positions; (b) moment reversal inside a segment braced on one flange; (c) a segment with zero moment (Eq. F1-1 is 0/0); (d) free ends and overhangs other than the one cantilever case of F1.
- **Where:** F1, Eq. F1-1.
- **Now:** (a) C_b = 1.0 unless given, with the largest member moment; (b) fall back to the segment between points where both flanges are braced; (c) a segment without moment is not checked; (d) only a cantilever with a built-in support is checked (C_b = 1.0, warping restraint there *assumed*), any other free end is `notChecked` unless C_b is given. These are engine rules. An earlier build took a rounding-error moment at a pinned end as a sign change in (b); forces below 10⁻⁹ of the member's scale are now zero. <!-- RECHECK -->
- **Change:** Whether (a) and (b) should be `notChecked` instead; LTB ratios of braced beams.

**Q25. Load height.**
Loads applied above the shear centre (top-flange loading) reduce the LTB strength; the Specification has no equation, only a Commentary remark.
- **Where:** Commentary to F1.
- **Now:** Not considered anywhere; stated as a limit.
- **Change:** Needs an input field and a rule from an opened source.

**Q26. F7.2(c): effective section modulus of a slender flange.**
The Specification defines the effective width only; the method for S_e (neutral-axis shift or not) is not prescribed.
- **Where:** F7.2(c), Eq. F7-3 to F7-5.
- **Now:** Ineffective width removed from both flanges, no neutral-axis shift (the simpler, conservative alternative; reproduces Design Example F.8B).
- **Change:** Flexural capacity of slender-flange boxes by about 2 %.

**Q27. F7.3(c): slender webs of boxes.**
Eq. F7-7 refers to the reduction factor of F5, written for I-shapes with h_c/t_w; whether the cap a_w ≤ 10 applies to boxes is not stated.
- **Where:** F7.3(c), Eq. F7-7, F5-6.
- **Now:** `notChecked`.
- **Change:** Whether the check is implemented.

**Q28. F7.4: lateral-torsional buckling of rectangular sections about which axis?**
The text does not restrict LTB to major-axis bending; only a User Note says it does not occur for square sections or minor-axis bending.
- **Where:** F7.4, Eq. F7-8 to F7-11.
- **Now:** Checked only for bending about local z when I about z is the larger; a rectangular section bent about local y as its strong axis is `notChecked`.
- **Change:** Coverage of rolled-90° boxes.

**Q29. Eq. F2-7 and the warping constant.**
r_ts needs C_w; the User Note alternatives are advisory. Which is the reference when C_w is not tabulated? Eq. F2-2 and F2-3 meet only within about 0.1 % at L_b = L_r because of the printed coefficients.
- **Where:** F2.2, Eq. F2-4 to F2-7.
- **Now:** r_ts from Eq. F2-7 with the section's C_w (for a parametric I-shape the thin-wall formula), or the tabulated r_ts when given. No smoothing at L_r.
- **Change:** LTB capacities by a fraction of a percent; tier-F coverage of Eq. F2-7.

**Q30. F3.2 for rolled shapes.**
Eq. F3-2 uses k_c although the rolled-flange limit of Table B4.1b has none; the strength steps slightly at the slender limit.
- **Where:** F3.2, Eq. F3-2.
- **Now:** As printed.
- **Change:** Nothing practical (no rolled shape has a slender flange in the usual grades).

**Q31. H1: where is the interaction checked, and with which strengths?**
Are P_r, M_rx, M_ry the maxima along the member (possibly at different sections) or concurrent values at each section, and which M_c (governed by LTB of one segment) pairs with which section?
- **Where:** H1.1, Eq. H1-1a/b; `design/aisc360_22/h_combined.py`.
- **Now:** Concurrent forces at every station and at the stationary points between stations, against the member-governing strengths: P_c = lowest of all compression limit states (including torsional buckling), M_c = lowest flexural strength over all LTB segments. Conservative.
- **Change:** H1 ratios of long members with several unbraced segments (lower with segment-wise pairing).

**Q32. H1 for round sections.**
H1-1 sums the two moment ratios; whether the resultant moment may be used for a round section is not stated.
- **Where:** H1.1; F8.
- **Now:** The resultant moment (one flexural term).
- **Change:** H1 ratios of round members under biaxial bending (up to a factor √2 lower than the sum).

**Q33. H1.2 and H1.3 options.**
The C_b increase under axial tension (H1.2) is optional and P_r varies along L_b; H1.3 (separate in-plane and out-of-plane checks) has an unclear scope.
- **Where:** H1.2, H1.3.
- **Now:** Neither is applied.
- **Change:** Ratios of beam-columns in tension or under mainly major-axis bending (lower with the options).

**Q34. Torsion.**
H3 gives a neglect rule (T_r ≤ 20 % of T_c) for hollow sections only and no method for open sections, where a frame element has no warping stress. When may torsion be ignored?
- **Where:** H3.1, H3.2, H3.3; `h_combined.py` (`torsion`).
- **Now:** H3 is not implemented. Engine rule, not AISC: torsion is ignored only if the St. Venant shear stress estimate is ≤ 1 % of 0.6F_y; otherwise `notChecked`, and the member is `incomplete`.
- **Change:** Implementing H3.1/H3.2 for hollow sections would replace the 1 % rule there by the Specification's 20 %; many generated members would stop being `incomplete`; tier F needs coverage either way.

**Q35. Deflection checks.**
Chapter L gives no limits, combinations or reference lengths. Is deflection measured from the chord through the member ends; what is the reference for a cantilever; are combinations analysed with reduced stiffness admissible; is second-order analysis required?
- **Where:** L1, L2; `design/aisc360_22/serviceability.py`.
- **Now:** Limits and combinations are user input; deflection relative to the chord through the design-member ends, at stations only; members with a free end are `notChecked`; combinations with stiffness reduction are checked and flagged; no drift check.
- **Change:** Whether cantilevers can be checked; the meaning of every deflection ratio.

**Q36. Web classification under combined axial force and bending.**
Table B4.1a is for compression and B4.1b for flexure; nothing says which classifies the web of a beam-column.
- **Where:** B4.1.
- **Now:** B4.1a for the Chapter E checks, B4.1b for the Chapter F checks.
- **Change:** Capacities of beam-columns with webs between the two limits.

**Q37. Compression with intermediate brace points and a varying force.**
The Specification does not say how to infer unbraced lengths from an analysis model, nor which force pairs with which segment.
- **Where:** C3, E2.
- **Now:** L_c = K × the longest segment, with the largest compression anywhere in the member (conservative).
- **Change:** Compression ratios of braced columns with varying force.

**Q38. Elastic buckling stress from a buckling analysis.**
E3 permits F_e from an elastic buckling analysis instead of KL/r. Should reference solutions use it?
- **Where:** E3 (definition of F_e).
- **Now:** Not implemented; always Eq. E3-4 with the given K and L.
- **Change:** A link between the buckling analysis and the design module.

**Q39. Shapes that are not implemented.**
The reading notes list further ambiguities for channels (only F2 covers major-axis flexure; shear-centre loading), tees (F9 scope, stem depth), single and double angles (which leg in Table B4.1a, principal axes, how E4, E5 and E7 combine, H2 point checks), and unsymmetric shapes.
- **Where:** B4.1 Cases 3, 4; E4-3, E4-4, E5; F9–F12; G3; H2.
- **Now:** All such members are `notChecked` for compression, flexure and shear.
- **Change:** Only relevant if those shapes are to be implemented; each would need its own decisions.

**Q40. Sampling of demands.**
Shear, axial force and deflection are taken at result stations; `resultStations: 2` is allowed by the schema. Should the minimum be higher, or should a warning be given?
- **Where:** `design/demands.py`; schema `analysis.resultStations`.
- **Now:** Moments and the H1 interaction are evaluated between stations as well; with two stations the deflection check is `notChecked`; a shear extremum between stations (possible when a varying line load changes sign) is not captured.
- **Change:** Schema minimum (a format change), or a design warning.

**Q41. Demands without moment data.**
When the design module is called with forces that carry no moments (direct calls with numbers only), should flexure be skipped?
- **Where:** `design.check_member`.
- **Now:** No flexure or interaction check is produced and the member is `incomplete`, with a warning.
- **Change:** The contract of the Python interface only.

## B. Analysis formulation

Choices of the engine that no opened reference prescribes.

**Q42. Second order with shear deformation: Engesser or Haringx?**
The shear force that enters the shear strain of a compressed member can be defined normal to the deformed axis or normal to the rotated section. The two differ by terms of the order N/(G·A_s) times the shear deformation. AISC 360-22 does not prescribe one, and no opened reference decides.
- **Where:** THEORY.md 10.2; `solve/beamcolumn.py`; `oracle/beam.py`.
- **Now:** The engine uses the Engesser form (and requires N > −G·A_s). The oracle offers both and refuses to choose. No outside solver covers the case, so such results cannot be verified.
- **Change:** Second-order results with shear deformation, in the fourth digit or less for steel frames.

**Q43. Rigid links and diaphragms in second order.**
A rigid body that transmits forces while it rotates adds a rotational stiffness term (−Σ f·r at the master). The engine's constraint is the linear relation and leaves it out; for rigid end offsets the engine does include the analogous term. Should it be included for constraints?
- **Where:** THEORY.md 4, 16; CROSSCHECK.md 3.5.
- **Now:** Linear constraint, with a warning on the results. Measured difference between the two statements: 7.6 % in the largest displacement of one test model with leaning columns and torsion. Such combinations are not cross-checked (one eligible solver) and never served.
- **Change:** Second-order results of every model with constraints; they would become cross-checkable against the stiff-member emulation instead of against OpenSees.

**Q44. Rigid offsets not along the member axis in second order.**
Only the axial end force acts on the rotating arm (symmetric form); terms of the kind moment × rotation are neglected. Modelling the arm as a stiff member gives a different answer. Is this the intended convention?
- **Where:** THEORY.md 10.7.
- **Now:** As described; no outside solver has an error bound for offsets in second order, so the case cannot be verified.
- **Change:** Second-order results with transverse offsets.

**Q45. Sliding (transverse translational) releases in second order.**
The moment that the axial force produces across a sliding connection depends on how the connection is built.
- **Where:** `solve/nonlinear.py`, code `analysis.release-second-order`.
- **Now:** Refused for second-order and buckling analyses.
- **Change:** Needs a defined connection model.

**Q46. Own buckling of truss-type members.**
A truss-type member is a pin-ended strut with the I of its section in second-order and buckling analyses: compressed beyond its Euler load it refuses the analysis and it appears as buckling modes. Should member buckling of an axial-only member be left to Chapter E instead? This matters for `compressionOnly` members used as contact elements and for placeholder sections.
- **Where:** THEORY.md 10.6; refusal `stability.beyond-critical-load`.
- **Now:** Refused, with a message naming the member and its Euler load.
- **Change:** Whether braced frames with slender compression diagonals can be analysed in second order; Ziemian frame 7 "Right"; MGZ Example 9.1.

**Q47. Diaphragm with slaves outside the master's plane.**
The diaphragm rule uses the in-plane part of the offset only, which is not a rigid-body constraint when a slave is off the plane and then violates global moment equilibrium. Should validation refuse it, or should the constraint use the full offset?
- **Where:** THEORY.md 4; code `constraint.diaphragm-out-of-plane` (warning).
- **Now:** Accepted with a warning; the solver refuses the result when its moment check fails.
- **Change:** Either an earlier, clearer refusal or a different constraint.

**Q48. Default conditioning threshold.**
`pivotTolerance` 10⁻¹⁰ refuses stable models with stiffness contrasts of about 10⁸ as "ill-conditioned", and two almost-free release springs at one member as a member mechanism.
- **Where:** THEORY.md 6, 7; `defaults.py`.
- **Now:** As stated; the user can change the tolerance.
- **Change:** Which stiff-member models run; verification substitutes already use 10⁻¹⁵.

**Q49. Unused nodes and rotations of truss-only nodes.**
They are removed from the solution, reported as zero and listed in a warning. Should a benchmark model with an unused node be refused instead?
- **Where:** THEORY.md 6.
- **Now:** Removed with a warning; a load on such a degree of freedom is refused as a mechanism (causes `no-stiffness` / `only-truss-members`).
- **Change:** Validation of generated and submitted models.

**Q50. A loaded node without members.**
Is it a validation error with its own code, or a mechanism?
- **Where:** `solve/stability.py`; VERIFICATION.md 4.
- **Now:** A mechanism of kind "node", cause `no-stiffness`. When the dataset was built the diagnosis tiers could not confirm it (12 attempts quarantined, the two curve balls that produce it never served); the diagnosis tiers have since been extended to keep loaded degrees of freedom. <!-- RECHECK -->
- **Change:** Whether `constraint_removed` and `constraint_wrong_master` can be served.

**Q51. Mechanisms that are neither rigid-body nor single-node.**
Three hinges in line or a pin-jointed portal are reported as one "distributed" entry with a count and node participations, because no canonical basis exists. Is that enough as a reference answer?
- **Where:** THEORY.md 7.
- **Now:** As described; the verification compares the count and the participating nodes only.
- **Change:** Scoring of mechanism diagnoses.

**Q52. Refused combinations.**
A refused combination is delivered as an exception carrying the other combinations. Should the results document hold an explicit record of the refusal?
- **Where:** `solve/analysis.py`; results format.
- **Now:** Exception (strict) or omission with a `NO RESULTS` warning.
- **Change:** The results format (an addition).

**Q53. Buckling factor quoted in a "beyond critical load" refusal.**
It is obtained by bisection on scaled axial forces, not from the linear buckling analysis; the two agreed to better than 0.2 % where compared.
- **Where:** `solve/nonlinear.py`.
- **Now:** The message gives the value as approximate and states the basis it uses.
- **Change:** Wording only.

**Q54. Station shear in second order.**
Reported shears are components along the undeformed axes. Should the results also carry the shear normal to the deformed axis for the Chapter G checks and the moment interpolation?
- **Where:** THEORY.md 10.8; `design/demands.py`.
- **Now:** Only the undeformed-axis components; the design module uses them as they are.
- **Change:** Shear ratios and between-station moments in second order (small).

**Q55. A station-level statics check inside the engine.**
The engine's invariants close at the member ends; a wrong interior station value passes them. The verification layer now recomputes every station (tier D). Should the engine do so on every solve?
- **Where:** THEORY.md 8, 17; `verify/stations.py`.
- **Now:** Only in tier D.
- **Change:** Protection of ordinary runs outside the verification pipeline.

**Q56. Overlapping rigid offsets.**
Offsets that reverse the flexible length are accepted silently; local x flips.
- **Where:** THEORY.md 17.
- **Now:** Accepted; results are consistent with the independent reference.
- **Change:** A warning or a refusal.

**Q57. Subdivision tolerance tied to the iteration tolerance.**
With an axial force varying along members, the comparison between subdivision levels has a round-off floor near 10⁻⁹, so tighter `convergenceTolerance` values end as "not converged".
- **Where:** THEORY.md 10.4, 17.
- **Now:** As described; verification models use 10⁻⁹ at the tightest.
- **Change:** Usable range of the tolerance.

**Q58. A slack one-sided member with its own loads.**
It carries its loads to its ends as a simply supported beam in first order, including a small axial part, with no second-order effect on that member.
- **Where:** THEORY.md 13.
- **Now:** As described.
- **Change:** Results of braced frames with self weight on tension-only braces (small).

**Q59. Tension-only members with second-order analysis.**
A braced frame with a linear buckling factor of 7.65 was refused at its own loads as beyond critical (reported factor about 0.38) when its strength combination was switched to second order. Intended or a defect?
- **Where:** `solve/nonlinear.py`; `examples/05_braced_frame.json`.
- **Now:** The build examined for this file analyses that case (three iterations, three slack braces). The generators still never combine the two. <!-- RECHECK -->
- **Change:** Whether generated braced frames may use second order with tension-only bracing.

**Q60. States in which all forces vanish.**
A rigid-body response to a support settlement gives forces that are rounding errors; ratios of two rounding errors made the engine refuse a valid solution and made the superposition check quarantine correct models.
- **Where:** THEORY.md 8.1, 8.2.
- **Now:** A round-off floor on the force scales has since been added to the invariants and to the superposition check; the settlement frame now runs. <!-- RECHECK --> The station checks of tier D have their own floor.
- **Change:** Confirm the floor (16 machine epsilons × the stiffness-times-displacement scale) as policy.

**Q61. Several mechanisms with idle one-sided members.**
The uniqueness test of the active set decides only the case of one mechanism.
- **Where:** THEORY.md 13.
- **Now:** Several mechanisms with idle members are refused.
- **Change:** Some valid braced models are refused.

**Q62. Section properties from dimensions against tabulated values.**
Thin-wall J and C_w, simple shear areas and sharp box corners differ from tabulated properties. Which is the reference for a named shape?
- **Where:** `sections.py`.
- **Now:** Computed values unless `overrides` or a database shape gives others; the approximations are listed in the output.
- **Change:** Agreement with published examples that use tabulated properties.

## C. Verification policy

Rules of the verification pipeline that decide what may be called verified.

**Q63. Tier C in second order: how is the error of the outside solvers bounded?**
Both outside solvers discretise the beam-column problem, so the comparison must allow for their error. The first rule used a leading-order (Richardson) estimate, which equals the error itself when the engine is exact; correct models then missed by under 1 % to a factor 2, and most second-order models were quarantined (261 attempts in the dataset build, 10 AISC and vendor problems). Options were: keep it; a stated safety factor; compare with the extrapolated result; a different estimate.
- **Where:** VERIFICATION.md 2 (tier C); CROSSCHECK.md 1.3, 3.4; `crosscheck/runner.py`; `verify/policy.py`.
- **Now:** Replaced, while this file was written, by a bound from three refinement levels of each outside solver (the last level-to-level difference, enlarged when the contraction is poor; no estimate when the levels do not contract); `estimateSafety` stays 1.0; the translators were corrected for a first-order term. A verified second-order record now means agreement within about three times the OpenSees discretisation error. <!-- RECHECK --> The stored dataset predates this.
- **Change:** The owner should confirm the new rule, or choose another. It decides whether second-order references exist at all, and how loose "verified" is for them (10⁻⁵ to 10⁻³).

**Q64. Printed precision against hand-rounded textbook values.**
The one tolerance for a published target is half a unit of its last printed digit. Hand calculations in a textbook carry three or four digits in their working, so their last printed digit is often unreliable. Is a written, source-level rule wanted for such sources (for example a tolerance from the digits carried in the working)?
- **Where:** VERIFICATION.md 5; `verification/problems/mgz/`.
- **Now:** No tolerance was changed: 18 problems of McGuire/Gallagher/Ziemian are quarantined only for this reason (differences 6·10⁻⁵ to 6·10⁻³; the engine agrees with the exact solution of the book's own equations).
- **Change:** Up to 18 problems would become active.

**Q65. Precision of listed inputs.**
The Ziemian data set lists coordinates and loads to five significant figures but displays moments to 0.01 kip·in. Should the tolerance of a target account for the precision of the inputs?
- **Where:** VERIFICATION.md 5; `verification/problems/ziemian/`.
- **Now:** It does not: five first-order moment problems that differ by at most 2.6·10⁻⁵ are quarantined.
- **Change:** Those five problems.

**Q66. "Printed precision" of a spreadsheet.**
The data set stores doubles and displays three or two decimals. Is the displayed number the printed value?
- **Where:** `verification/problems/ziemian/README.md`.
- **Now:** Yes: the target is the number as displayed, the stored value is recorded in the note. With the stored doubles and the floor 10⁻⁶ nearly all first-order targets would be quarantined too.
- **Change:** The status of the 14 active Ziemian problems.

**Q67. Tier F: what must be covered, and at what granularity?**
(a) Should classification cases and advisory checks (slenderness, torsion screening) block a design reference? (b) Should a `notChecked` item with zero demand block, or only limit states with a demand ("used" = utilisation above 10⁻⁹, as implemented)? (c) A reproduced strength credits every clause item of its check, including classification cases reached with ratios the test invented. (d) Side (a) is recognised by a test naming convention. (e) Is coverage per equation the right unit, or per limit state and shape?
- **Where:** VERIFICATION.md 2 (tier F), 7; `verify/tier_f.py`; `tools/build_design_coverage.py`.
- **Now:** As in (b) and (c); at the time of writing 61 of 80 items are fully covered <!-- RECHECK -->; uncovered items include E4-1/E4-2, F2-7, F3-2, F6-2 to F6-4, F7-5, F7-6, F8-2 to F8-4, G2-3, G2-10, G2-11, Table B4.1a Case 8, Table B4.1b Cases 11 and 21, Sections L1/L2 and the torsion screening rule.
- **Change:** Whether most generated models can have a verified design reference; today almost none can, so design tasks are nearly absent from the benchmark.

**Q68. How are the tier-F gaps to be closed?**
Several limit states have no Design Example that prints a comparable strength (torsional buckling strength, slender flanges, box cases, round local buckling). Is a hand calculation alone enough for those, or must an opened source supply printed values?
- **Where:** `verify/design_coverage.json`; REFERENCES.md 2.1 (Manual, Shapes Database).
- **Now:** Both kinds are required; such items stay uncovered.
- **Change:** The same as the previous question.

**Q69. Problems recomputed from a published statement: tier A or tier B?**
For the vendor problems and five textbook problems the targets are computed in this repository; only the statement is published. They are labelled tier B with their own basis.
- **Where:** VERIFICATION.md 2 (tier B); `verification/problems/staad/`, `ansys/`, five in `mgz/`.
- **Now:** Tier B, basis `recomputed-from-published-statement`, counted separately in the report; the documents say they are evidence of the same kind as tier A.
- **Change:** Relabelling 30 problems to tier A; the headline count of "published problems reproduced".

**Q70. The calibrated AISC Commentary problems.**
The Commentary gives no section properties for its benchmark member, so I and the shear flexibility are calibrated on the Commentary's own zero-axial-load values. Should these 32 problems count at all, and do the 16 with-shear problems count given that the Commentary's tolerance is wider than the whole shear effect?
- **Where:** `verification/problems/aisc-commentary/README.md`.
- **Now:** Counted as active tier-B problems; six targets are marked "reproduced by construction"; each with-shear problem states that it does not test shear deformation.
- **Change:** 32 problems; with a sourced section table they would become real tests.

**Q71. Should the tableau solver be a required tier?**
It is the only solver that reads the raw input file itself; every other route takes the model as resolved by `benchfea.model`. It is by the same author as the engine.
- **Where:** VERIFICATION.md 2 (extra evidence); `tests/independent/tableau.py`; LIMITS.md 4.5.
- **Now:** Extra evidence only (`--extra-evidence`); never changes a status. The dataset audit ran it on all 166 stable models (largest difference 5.6·10⁻¹²).
- **Change:** Protection against an input-interpretation error common to all routes; cost per verification; it lives under `tests/`, not in the package.

**Q72. Buckling factors have no second independent solver.**
Neither translator computes buckling. May buckling references be published on tier A, D and E only, or should an eigenvalue route through OpenSees tangent matrices be added?
- **Where:** CROSSCHECK.md 4; VERIFICATION.md 7.
- **Now:** A model with buckling results can never pass tier C; no buckling factor is served.
- **Change:** Whether buckling tasks can exist.

**Q73. Features that only one outside solver represents.**
Shear deformation and partial axial/shear/torsion releases have one eligible solver (or none). Is one outside solver plus the in-repository oracle acceptable for them, or must a third solver be obtained?
- **Where:** CROSSCHECK.md 4, 5.
- **Now:** Not cross-checked, not publishable as verified.
- **Change:** Whether such models can be served.

**Q74. Slack members that carry member loads, in tier C.**
The outside solvers get the linear model of the engine's active set; loads on removed slack members must reach the nodes by the engine's rule.
- **Where:** VERIFICATION.md 2 (tier C substitutions).
- **Now:** When the examples were written this case was refused by tier C; the loads of slack members are now computed from the load definitions and passed as nodal loads. <!-- RECHECK -->
- **Change:** Whether tension-only models with self weight can be verified; confirm the lumping rule as part of the contract.

**Q75. Tolerance of the superposition and reciprocity checks.**
10⁻⁹ was missed at about 1.6·10⁻⁹ by one first-order model. Is the tolerance meant to be that tight?
- **Where:** `verify/policy.py` (tier D).
- **Now:** Unchanged; the attempt stays quarantined. <!-- RECHECK -->
- **Change:** A handful of quarantines.

**Q76. Design Examples: rounded chains and internal inconsistencies.**
Printed strengths come from chains of three-figure intermediates, from table look-ups, and are occasionally inconsistent (a printed 155 against 258/1.67 = 154.5; a printed 598 from an interpolated table against 594 from the equation). Match the printed value with a tolerance, or replicate the rounding?
- **Where:** `tests/design_examples/`.
- **Now:** One tolerance, 1 % on a printed value, never tuned per case; where the print is inconsistent the test asserts the equation's value and records the discrepancy. Largest deviation 0.95 %.
- **Change:** Nothing unless the tolerance is tightened.

**Q77. Design Examples that do not list their inputs.**
Several examples assert compactness or take strengths and lengths from Manual tables without listing the section properties. The tests supply stated placeholders or take the tabulated value as input.
- **Where:** `tests/test_design_*.py`, `tests/test_design_review.py`.
- **Now:** Such cases are partial reproductions and are marked so; a placeholder is never taken from memory. They still give tier-F credit to the check whose strength is reproduced.
- **Change:** Tier-F credit of a few items; resolved by the Manual or the Shapes Database.

**Q78. Mechanism diagnoses that cannot be verified.**
A mechanism through slack one-sided members, a "distributed" mode (counted, not tested as a vector) and an ill-conditioned classification near the tolerance.
- **Where:** VERIFICATION.md 4.
- **Now:** Slack mechanisms are always quarantined; the others are verified by count and participation only, or fail tier N.
- **Change:** Which mechanism tasks can be served.

**Q79. Quarantine is per problem, not per target.**
One missed target quarantines the whole problem, which is why vendor-printed values and load levels were split into separate problems.
- **Where:** `verify/problems.py`.
- **Now:** Per problem.
- **Change:** The number and naming of problems.

**Q80. Problems that need a transformed comparison.**
An inclined roller is modelled in rotated axes; values the source prints in other axes would need a propagated tolerance.
- **Where:** `verification/problems/mgz/` (Example 13.5).
- **Now:** Only the directly comparable value is a target.
- **Change:** Two further targets.

**Q81. Scope of a tolerance stated by the source.**
The Commentary's 3 % / 5 % is applied to all its tabulated values (with and without shear, US and SI tables). Is that the intended scope?
- **Where:** VERIFICATION.md 5, rule 2.
- **Now:** Applied to all 64 values; actual differences are at most 0.38 % and 0.88 %.
- **Change:** Nothing at the present differences.

## D. Published sources

Questions the opened sources leave open. Until answered, the affected problems are quarantined or excluded.

**Q82. AISC Commentary C2.1: section properties and method of the benchmark solution.**
Which I_x, area, shear area and shear modulus were used for the benchmark member; is the reference solution closed form or computed, and with which shear formulation? The support idealisation of Case 1 is read from a sketch.
- **Where:** Commentary C2.1, p. 16.1-337.
- **Now:** Properties calibrated on the Commentary's own values (Q70).
- **Change:** The 32 Commentary problems become independent tests.

**Q83. AISC Design Examples C.1A–C: the analysis model behind the figures.**
Which shear areas and beam area did the software use; does it include shear deformation; how are the leaning column and its link modelled?
- **Where:** Design Examples v16.0, pp. C-2 to C-17.
- **Now:** All ten problems quarantined: the model buildable from the stated data is 1.8–2.2 % stiffer in drift. A parameter study shows web-sized shear areas would reproduce the prints; those trial values are in no model.
- **Change:** Ten problems could be rebuilt.

**Q84. Ziemian frame 3: first-order mismatch.**
First-order displacements differ by up to 1.4·10⁻³ with alternating sign between column lines, although the critical load ratio reproduces to 6.5·10⁻⁶. Are the loads of the published runs the ones listed?
- **Where:** `verification/problems/ziemian/`.
- **Now:** Both first-order displacement problems of frame 3 quarantined; cause not identified.
- **Change:** Two problems; confidence in the frame-3 model.

**Q85. Ziemian: member numbering.**
The data set numbers members without defining the numbering; for frames 3, 5 and 7 it lists fewer members than there are member lines.
- **Where:** Results spreadsheet.
- **Now:** Inferred for frames 1, 2 and 9 from first-order moments and plastic moments (so the first-order moment targets are not independent evidence of the mapping); no moment targets for frames 3, 5, 7.
- **Change:** Moment targets for three more frames.

**Q86. Ziemian frame 1: unrounded loads.**
All first-order results are below the engine's by the same factor 1 − 1.73·10⁻⁶.
- **Where:** Model spreadsheet.
- **Now:** Within tolerance for displacements; contributes to the quarantined moment problems.
- **Change:** The moment problems of frame 1.

**Q87. Ziemian: comparing with the second-order (work-control) results.**
The data set's analysis is a large-displacement analysis; the engine's is small-displacement. The two theories differ by 0.4–2.7 % on these frames, and a corotational analysis still differs from the data set by 0.1–0.7 %, which is unresolved. Is there a tolerance the authors consider appropriate, and should the engine's theory be accepted as it is?
- **Where:** `verification/derive/ziemian_diagnosis.py`.
- **Now:** All 12 second-order problems quarantined as a difference of theory.
- **Change:** Twelve problems; whether a large-displacement option is wanted in the engine.

**Q88. Ziemian: critical load ratios of frames 1 and 7.**
The engine is lower by 9.5·10⁻⁵ and 6.1·10⁻⁵, which crosses the third displayed decimal; the data set's own discretisation (four cubic elements) reproduces its values.
- **Where:** Summary tab.
- **Now:** Two problems quarantined: the source value is a discretised approximation.
- **Change:** Two problems, under a rule for discretised source values.

**Q89. McGuire/Gallagher/Ziemian Example 4.15: which input?**
The support area is printed as an expression and as its rounded value; with the expression the printed results are met, with the rounded value they miss in the last digit.
- **Where:** Example 4.15.
- **Now:** Quarantined rather than choosing the reading that passes.
- **Change:** One problem.

**Q90. Examples 9.5 and 9.6: second-order moments.**
The book's incremental scheme, step size and geometry update are not stated; its values are 3–6 % below the equilibrium solution of its own element matrices. A label "1660" on the moment diagram of Example 9.6 is not identified.
- **Where:** Examples 9.5, 9.6.
- **Now:** Two problems quarantined; the label is not used.
- **Change:** Two problems.

**Q91. Source values computed with a stated discretisation.**
Examples 9.5 and 9.8 print critical loads from a few elements with an approximate geometric stiffness (0.3 % and 0.8 % above the exact value). Is a rule wanted for such targets?
- **Where:** Examples 9.5, 9.8.
- **Now:** Quarantined: "discretised approximation, no stated error bound". Where the book also prints the theoretical value, that is the target.
- **Change:** Two problems (and the two Ziemian ones above).

**Q92. Suspected misprints and ambiguities in the book.**
Example 4.12 (a reaction printed 0.423 against 0.413 in its own figure), Example 13.2 (a rotation 2.9 % off its own equations), Example 9.3 (a coefficient printed two ways), Example 13.3 (sense of joint moments drawn as arrows), the sign in the general statement of Eq. 13.29, Example 5.3 / 5.7 (one coefficient printed two ways).
- **Where:** As named.
- **Now:** Affected problems quarantined with the evidence; where a reading was needed it is stated in the problem.
- **Change:** Up to four problems.

**Q93. Examples excluded for missing data.**
Example 9.7 (support type only in a small figure), 9.13 and 13.9 (section properties not printed), 9.14 (coordinates not printed). Chapter 6 was not reviewed.
- **Where:** As named.
- **Now:** Excluded; listed in the README.
- **Change:** Further problems, with a shapes table or the missing data.

**Q94. STAAD 3D frame: local-axis convention.**
The page's largest moments about local y and z depend on the vendor's convention for an inclined member, which the page does not state.
- **Where:** STAAD.Pro Help, "3D Frame Max Forces".
- **Now:** Only the largest axial force is a target.
- **Change:** Two more targets.

**Q95. Vendor pages whose sketches could not be read.**
Supports and load positions of some pages are given only in a sketch. VM217's fixed bases are deduced from the page's two formulas; VM21's supports and VM4's load direction rest on the text. Where a page's text and input listing disagree, the listing is followed.
- **Where:** `verification/problems/staad/`, `ansys/`.
- **Now:** Used only where text or listing fix the problem; each reading is written in the problem's analysis text.
- **Change:** Confirmation of two statements; further pages.

**Q96. Design Examples: errata status and numbering.**
No revision date is printed beyond "Version 16.0"; one example cites an equation of F7 under a different number than the Specification prints.
- **Where:** Front matter; Example H.5C.
- **Now:** The Specification's numbering is cited.
- **Change:** Nothing unless an erratum changes a printed value.

**Q97. AISC Shapes Database: column names and units.**
The reader's header names, the rule for duplicated headers and the length unit are assumptions; the database could not be opened.
- **Where:** `shapes.py`.
- **Now:** Overridable assumptions, stated in the module and in LIMITS.md 3.6. Only the design fillet dimension is read for the web depth; tabulated width-to-thickness ratios are not read.
- **Change:** Whether `type: standard` sections can be trusted without a manual check.

**Q98. Parts of opened sources not yet used.**
STAAD: about nine pages not transcribed (one with 28 end moments of a hand method); ANSYS: the remaining pages; Design Examples: E.5–E.8, E.13, E.14, F.10–F.13 read from the index only.
- **Where:** READMEs under `verification/problems/`.
- **Now:** Listed as not examined.
- **Change:** More problems and tier-F coverage.

## E. Benchmark design

Policy choices of the task set, the scoring and the split.

**Q99. Design-only curve balls and tier F.**
A design-variable curve ball changes no analysis result. Its *diagnose* answer is the change record, not a design value. May it be served when only the analysis tiers pass?
- **Where:** BENCHMARK.md 1, 9; `bench/tasks.py`.
- **Now:** No: the reference is verified in design mode, tier F is not covered, and all five `design_*` curve balls are never served.
- **Change:** Five curve balls enter the dataset.

**Q100. Second order in generated models.**
Second order is drawn by default at levels 3–5, and with the earlier tier-C rule nearly all such models were quarantined. Should it be default-off until Q63 is settled?
- **Where:** `bench/generate.py`; DIFFICULTY.md.
- **Now:** Default on; the stored dataset therefore has no stable second-order reference. <!-- RECHECK -->
- **Change:** Composition of the dataset.

**Q101. Base models shared between train and held-out.**
In the stored dataset 12 held-out instances (34 tasks) use base models that train tasks also use; for 10 repair tasks the reference answer is also a train reference answer.
- **Where:** BENCHMARK.md 7; `bench/splits.py`.
- **Now:** The generator seed now carries the split, and the held-out pairs are chosen by a rule on names only; the stored dataset is that of the earlier rule until rebuilt. <!-- RECHECK -->
- **Change:** A rebuild (about two hours) and refreshed numbers.

**Q102. `invalid-model` curve balls.**
A curve ball that makes the model invalid has no analysis to verify. Should a diagnose task against the validator's problem list be served?
- **Where:** BENCHMARK.md 9.
- **Now:** Never served.
- **Change:** One curve ball (`constraint_wrong_slave`) in generated families.

**Q103. Curve balls whose reference cannot be verified.**
`constraint_removed`, `constraint_wrong_master` and `offset_sign_reversed` produce models whose reference is quarantined by construction. Keep them applicable with that recorded, or refuse them by default?
- **Where:** CURVEBALLS.md.
- **Now:** Applicable; the record says the reference is not expected to verify; with verification required they are refused.
- **Change:** Three curve balls.

**Q104. What counts as a change that matters?**
Is a change of a design ratio alone material? It is counted only for curve balls that change a member's resistance; for demand-only curve balls an analysis change of at least 5 % (or a pass/fail flip) is required.
- **Where:** CURVEBALLS.md; `bench/curveballs/measure.py`.
- **Now:** As described; thresholds are the same for every model.
- **Change:** Applicability of several curve balls.

**Q105. Should verification be the default when applying a curve ball?**
It is rigorous but takes seconds to minutes per candidate.
- **Where:** `bench/curveballs` (`apply(verify=)`).
- **Now:** Off by default in `apply`; always on in the dataset builder.
- **Change:** Convenience against the risk of an unverified reference outside the builder.

**Q106. Scoring rules introduced by the implementation.**
(a) A mechanism claimed for a stable model halves the diagnose score. (b) A correct repair that the pipeline cannot verify scores 0. (c) Tolerances: predict 1 % full / 5 % zero with tools, 5 % / 25 % without; repair 2 % / 10 %. (d) Diagnose weights. (e) Only quantities moved by at least 15 % are asked for in a repair.
- **Where:** BENCHMARK.md 4; `bench/scoring.py`.
- **Now:** As listed.
- **Change:** Scores; none of the references.

**Q107. Licence.**
The repository has no licence file, so no licence is granted for code or data; the Ziemian-derived problem files carry CC BY 4.0 obligations.
- **Where:** `data/DATA_CARD.md`; REFERENCES.md R4.
- **Now:** None.
- **Change:** Whether anyone else may use the benchmark.

**Q108. Dataset in the repository.**
24 MB of JSON, built from seeds 0 and 4 only, so each train pair has one seed. Keep it in the repository, compact it, or build more seeds (about an hour each)?
- **Where:** `data/`.
- **Now:** As built.
- **Change:** Repository size; statistical weight of scores.

**Q109. Features never drawn by default.**
Buckling requests, tension-only bracing and second order with diaphragms are opt-in because they cannot pass tier C; default models therefore contain none.
- **Where:** DIFFICULTY.md.
- **Now:** Opt-in.
- **Change:** Coverage of those features by the benchmark.

**Q110. Generator acceptance thresholds.**
Buckling factor at least 2, sway amplification at most 1.5, half the members with ratio 0.3–1.1, generic section catalogue, loads not from a loading code.
- **Where:** DIFFICULTY.md 3, 5.
- **Now:** As listed, never relaxed per family.
- **Change:** What the models look like.

**Q111. Generator defects found by downstream users.**
Three seeds of one family raised an invalid-model error (an empty design group) instead of rejecting the draw; generated documents shared one object between a member's two end releases.
- **Where:** `bench/families/`.
- **Now:** Logged as `generator-error` in the stored dataset (24 attempts). The build examined for this file generates the first of those seeds without error and without the shared object. <!-- RECHECK -->
- **Change:** A rebuild changes which base models are used.

**Q112. Numbers in the documents after a rebuild.**
The dataset tables of BENCHMARK.md and the data card were written by a one-off script and must be refreshed by hand.
- **Where:** BENCHMARK.md 8.
- **Now:** Stale after any rebuild.
- **Change:** A tool or a test that regenerates them.

## F. Input and output format

Questions whose answer is a change of the format (governed by FORMAT_VERSIONING.md).

**Q113. A field for fabrication and for hollow-section data.**
The format cannot say rolled or built-up, hollow section or welded box, fillet dimension, corner radius, nominal or design thickness. Should a field be added (a minor format change), and until then should a parametric `box` be `notChecked` rather than checked as a welded box?
- **Where:** `schema/`; `design/section.py`.
- **Now:** Decided from the section type and reported in every check; a `box` is checked as a welded box.
- **Change:** Classification and shear of boxes and of I-shapes without fillet radius.

**Q114. Storeys, levels and gravity cases.**
Notional loads, B2 and drift checks need levels and a marker for gravity load cases.
- **Where:** `schema/`.
- **Now:** Neither exists.
- **Change:** Q2 and Q7.

**Q115. Features met in the published problems that the format lacks.**
Initial strain or lack of fit, temperature loads, inclined supports, general linear constraints, a continuous elastic foundation, the unit cm.
- **Where:** READMEs under `verification/problems/`.
- **Now:** Absent; affected examples are excluded or modelled by an exactly equivalent loading or in rotated axes, stated in the problem.
- **Change:** Which published problems can be encoded.

**Q116. Design variables that no check uses.**
`role`, `Cmx`, `Cmy`, `sideswayInhibited`.
- **Where:** `designparams.py`.
- **Now:** Resolved and reported, unused.
- **Change:** Remove them (a major format change) or implement B1.

**Q117. Additions to the output documents without a version change.**
Keys were added to the results and design documents (`constraintForces`, `nodalEnvelopes`, `solver`, `status`, `byCombination`, check details) while the format version stayed 1.0.0.
- **Where:** FORMAT_VERSIONING.md "Output documents".
- **Now:** Treated as additive.
- **Change:** Whether output documents get their own minor versions.

**Q118. Exit codes of the command line.**
`run`: 1 refused, 2 unreadable, 3 written but flagged. `design`: 0 for pass and conditional, 4 for fail and incomplete. `verify`: 5 quarantined. Confirm the conventions, in particular 0 for `conditional`.
- **Where:** `cli.py`.
- **Now:** As listed.
- **Change:** Scripts that rely on them.

**Q119. What the shapes reader takes from the database.**
Should it also read tabulated width-to-thickness ratios, r_ts and h_o?
- **Where:** `shapes.py`.
- **Now:** Only the columns listed in the module, and the design fillet dimension if present.
- **Change:** Agreement with the Design Examples for database shapes.
