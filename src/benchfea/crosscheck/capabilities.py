"""Machine-readable record of what each independent solver can represent, and how.

Every entry is ``{"support": ..., "how": ..., "evidence": ...}`` with ``support`` one of

* ``"exact"``        the solver (through the translator) solves the same mathematical problem;
* ``"approximate"``  converges to it; the translator reports an estimated error;
* ``"unsupported"``  cannot be represented; the translator records it and the result is not a cross-check;
* ``"not-run"``      the solver could not be executed in this project (formulation from its documentation only).

``evidence`` is ``"measured"`` (numerical experiment with the installed version, reproduced by the
tests named in docs/CROSSCHECK.md), ``"documentation"`` (official documentation, URL in
``SOURCES``) or ``"source"`` (statement of the solver's own source/manual, not executed).
Versions the measurements refer to: OpenSeesPy 3.8.0, PyNiteFEA 3.2.1.

    from benchfea.crosscheck import capability
    capability("pynite", "shearDeformation")["support"]      # "unsupported"
"""
from __future__ import annotations

__all__ = ["CAPABILITIES", "FORMULATION", "SOURCES", "FEATURES", "capability"]

FEATURES = (
    "elementFormulation", "shearDeformation", "localAxes", "trussMembers", "rigidOffsets",
    "releasesFull", "releasesPartialRotational", "releasesPartialOther", "rigidLink", "rigidDiaphragm",
    "supportSprings", "prescribedDisplacements", "loadNodal", "loadMemberPointForce", "loadMemberPointMoment",
    "loadMemberUniform", "loadMemberPartialUniform", "loadMemberLinearlyVarying", "selfWeight",
    "secondOrder", "secondOrderWithShearDeformation", "geometricStiffness", "axialForceIteration",
    "linearBuckling", "tensionCompressionOnly", "stiffnessReduction", "notionalLoads",
)

SOURCES = {
    "opensees": [
        "https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/ElasticTimoshenkoBeam.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/LinearTransf.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/pdelta.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/eleload.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/rigidLink.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/rigidDiaphragm.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/equalDOF.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/TransformationMethod.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html",
        "https://openseespydoc.readthedocs.io/en/latest/src/eleForce.html",
        "https://opensees.berkeley.edu/wiki/index.php/Linear_Transformation",
        "https://opensees.berkeley.edu/wiki/index.php/PDelta_Transformation",
        "https://opensees.berkeley.edu/wiki/index.php/RigidDiaphragm_command",
        "https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/elasticBeamColumn.html",
    ],
    "pynite": [
        "https://pynite.readthedocs.io/en/latest/member.html",
        "https://pynite.readthedocs.io/en/latest/PDelta.html",
        "https://pynite.readthedocs.io/en/latest/analysis.html",
        "https://pynite.readthedocs.io/en/latest/node.html",
        "https://pynite.readthedocs.io/en/latest/stability.html",
    ],
    "frame3dd": [
        "https://frame3dd.sourceforge.net/",
        "https://svn.code.sourceforge.net/p/frame3dd/code/trunk/doc/Frame3DD-manual.html",
    ],
}

FORMULATION = {
    "opensees": {
        "package": "openseespy", "runnable": True,
        "element": "elasticBeamColumn: prismatic Euler-Bernoulli; ElasticTimoshenkoBeam: with shear deformation",
        "geometricStiffness": "geomTransf PDelta: chord-rotation ('string') term N/L on the end translations only; "
                              "no P-delta within an element; N = the axial force of element end j WITHOUT the "
                              "axial part of an element load (measured); axial force iterated by Newton",
        "constraints": "rigidLink('beam'), rigidDiaphragm, equalDOF enforced exactly by the Transformation "
                       "handler (retained node must not itself be constrained)",
        "licence": "see site-packages/openseespy/LICENSE.md (optional dependency, never vendored)",
    },
    "pynite": {
        "package": "PyNiteFEA", "runnable": True,
        "element": "prismatic Euler-Bernoulli 12x12, releases by static condensation",
        "geometricStiffness": "consistent (cubic) geometric stiffness in both bending planes plus an axial term and "
                              "a torsional term; analyze_PDelta is one solve with first-order axial forces.  The "
                              "worker removes the two extra terms through the section properties and iterates the "
                              "axial forces with PyNite's own member matrices and solver",
        "constraints": "none (no rigid links, diaphragms or offsets)",
        "licence": "MIT",
    },
    "frame3dd": {
        "package": None, "runnable": False,
        "element": "3D frame with optional shear deformation (manual / source, not executed)",
        "geometricStiffness": "consistent geometric stiffness including shear terms, axial terms zeroed, iterated "
                              "(manual / source, not executed)",
        "constraints": "none; approximate rigid node zones; no end releases",
        "licence": "GPL v3 or later (external program only)",
    },
}


def _e(support, how, evidence="measured"):
    return {"support": support, "how": how, "evidence": evidence}


CAPABILITIES = {
    "opensees": {
        "elementFormulation": _e("exact", "elasticBeamColumn (A, E, G, J, Iy, Iz)"),
        "shearDeformation": _e("exact", "ElasticTimoshenkoBeam (E, G, A, J, Iy, Iz, Avy, Avz), first order only"),
        "localAxes": _e("exact", "vecxz = benchfea local z; y = vecxz x x, z = x x y"),
        "trussMembers": _e("exact", "elasticBeamColumn with bending and torsion stiffness x 1e-15; second order "
                                    "with an axial member load: a subdivided pin-ended beam-column (approximate)"),
        "rigidOffsets": _e("exact", "-jntOffset on Euler-Bernoulli elements; explicit node + rigidLink('beam') for "
                                    "ElasticTimoshenkoBeam, which mishandles -jntOffset; second order: "
                                    "approximate without estimate"),
        "releasesFull": _e("exact", "bending: parallel elements with -releasey/-releasez; axial/shear/torsion: "
                                    "duplicate node + equalDOF (member axes parallel to global axes, no offset, "
                                    "joint not a constraint slave), otherwise penalty zeroLength (approximate)"),
        "releasesPartialRotational": _e("exact", "exact decomposition into parallel elements with release flags "
                                                 "(fef.bending_split)"),
        "releasesPartialOther": _e("exact", "zeroLength spring between duplicate nodes; rigid directions by "
                                            "equalDOF when the member axes are global, else penalty (approximate)"),
        "rigidLink": _e("exact", "rigidLink('beam') + constraints('Transformation'); second order: a linear "
                                 "kinematic relation, as in benchfea (no geometric effect of the transmitted forces)"),
        "rigidDiaphragm": _e("exact", "rigidDiaphragm + constraints('Transformation'); slaves must lie in the "
                                      "master's plane"),
        "supportSprings": _e("exact", "zeroLength to a fixed twin node"),
        "prescribedDisplacements": _e("exact", "sp in the load pattern"),
        "loadNodal": _e("exact", "load"),
        "loadMemberPointForce": _e("exact", "eleLoad -beamPoint on unreleased Euler-Bernoulli elements; otherwise "
                                            "equivalent nodal loads (fef)"),
        "loadMemberPointMoment": _e("exact", "no native load: equivalent nodal loads (fef)"),
        "loadMemberUniform": _e("exact", "eleLoad -beamUniform on unreleased Euler-Bernoulli elements; otherwise "
                                         "equivalent nodal loads (fef)"),
        "loadMemberPartialUniform": _e("exact", "eleLoad -beamUniform with aOverL, bOverL on unreleased "
                                                "Euler-Bernoulli elements; otherwise equivalent nodal loads (fef)"),
        "loadMemberLinearlyVarying": _e("exact", "OpenSees 3.8 ignores the end intensities in 3D: equivalent nodal "
                                                 "loads from exact fixed-end forces (fef)"),
        "selfWeight": _e("exact", "uniform line load computed by the translator"),
        "secondOrder": _e("approximate", "PDelta transformation with member subdivision, Newton; axial member "
                                         "loads at the nodes (O(1/n^2)); error bounded from three subdivision levels"),
        "secondOrderWithShearDeformation": _e("unsupported", "ElasticTimoshenkoBeam ignores the geometric "
                                                             "nonlinearity of the transformation"),
        "geometricStiffness": _e("approximate", "chord-rotation term only (see FORMULATION)"),
        "axialForceIteration": _e("exact", "Newton iterations on the current axial forces"),
        "linearBuckling": _e("unsupported", "no Kg eigen-solver; feasible from two tangent matrices (printA), not "
                                            "implemented in the translator"),
        "tensionCompressionOnly": _e("unsupported", "not implemented in the translator", "source"),
        "stiffnessReduction": _e("exact", "stiffnesses scaled by the translator; a 'computed' tau_b must be supplied"),
        "notionalLoads": _e("unsupported", "must be supplied by the caller as nodal loads", "source"),
    },
    "pynite": {
        "elementFormulation": _e("exact", "add_member with add_section(A, Iy, Iz, J)"),
        "shearDeformation": _e("unsupported", "no shear area input; transverse shear deformation not considered"),
        "localAxes": _e("exact", "rotation from PyNite's own triad (Y-up rule), verified per member in the worker"),
        "trussMembers": _e("exact", "moments released at both ends, torsion at end j; second order with an axial "
                                    "member load: subdivided like a frame member (approximate)"),
        "rigidOffsets": _e("approximate", "stiff member; measured optimum multiplier about 1e4 (round-off beyond)"),
        "releasesFull": _e("exact", "def_releases (boolean, any of the 12 end DOFs); second order: approximate "
                                    "(geometric stiffness condensed separately)"),
        "releasesPartialRotational": _e("exact", "parallel members (fef.bending_split) with equivalent nodal loads"),
        "releasesPartialOther": _e("unsupported", "releases are boolean"),
        "rigidLink": _e("approximate", "stiff member master -> slave; second order: unsupported (the stiff member "
                                       "carries geometric stiffness that the linear constraint does not have)"),
        "rigidDiaphragm": _e("approximate", "stiff members with out-of-plane releases at the slave; coplanar "
                                            "slaves only; second order: unsupported (as rigidLink)"),
        "supportSprings": _e("exact", "def_support_spring"),
        "prescribedDisplacements": _e("exact", "def_node_disp"),
        "loadNodal": _e("exact", "add_node_load"),
        "loadMemberPointForce": _e("exact", "add_member_pt_load, local axes, absolute position"),
        "loadMemberPointMoment": _e("exact", "add_member_pt_load with Mx / My / Mz"),
        "loadMemberUniform": _e("exact", "add_member_dist_load"),
        "loadMemberPartialUniform": _e("exact", "add_member_dist_load with x1, x2"),
        "loadMemberLinearlyVarying": _e("exact", "add_member_dist_load with w1 != w2 on a member WITHOUT "
                                                 "intermediate nodes (PyNite's own splitting at intermediate nodes "
                                                 "gives wrong results for varying loads; the translator never "
                                                 "relies on it)"),
        "selfWeight": _e("exact", "uniform line load computed by the translator"),
        "secondOrder": _e("approximate", "analyze_PDelta with explicit sub-members, then the axial forces iterated "
                                         "to convergence; error bounded from three subdivision levels"),
        "secondOrderWithShearDeformation": _e("unsupported", "no shear deformation"),
        "geometricStiffness": _e("approximate", "consistent cubic Kg; PyNite's axial term (EA + N)/L and torsional "
                                                "term N (Iy + Iz)/A are removed through A' and J' (residual measured)"),
        "axialForceIteration": _e("approximate", "PyNite itself: none; the worker repeats PyNite's second solution "
                                                 "step with the current axial forces until the displacements stop "
                                                 "changing (last change recorded)"),
        "linearBuckling": _e("unsupported", "no buckling eigen-analysis", "documentation"),
        "tensionCompressionOnly": _e("unsupported", "PyNite has tension/compression-only members, but the "
                                                    "translator does not use them", "source"),
        "stiffnessReduction": _e("exact", "stiffnesses scaled by the translator; a 'computed' tau_b must be supplied"),
        "notionalLoads": _e("unsupported", "must be supplied by the caller as nodal loads", "source"),
    },
    "frame3dd": {f: _e("not-run", "not executed in this project (no binary; GPL program)", "source")
                 for f in FEATURES},
}


def capability(solver: str, feature: str | None = None):
    """The record of one solver, or of one feature of it."""
    return CAPABILITIES[solver] if feature is None else CAPABILITIES[solver][feature]
