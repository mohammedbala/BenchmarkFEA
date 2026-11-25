"""benchfea model + combination  ->  OpenSeesPy model (docs/CROSSCHECK.md has the full account).

``translate`` runs in the calling process and imports nothing from OpenSees: it produces a *job*, a
plain JSON list of ``openseespy.opensees`` calls plus the outputs to read.  ``execute(job)`` performs
the calls; it is meant to run in a subprocess (``benchfea.crosscheck.runner``) because OpenSees
terminates the interpreter on some input errors.  The mapping of the raw output back to benchfea
conventions (``Translation.results``) again needs no OpenSees.

How each feature is represented (every statement about OpenSees 3.8 behaviour was measured and is
reproduced by tests/test_crosscheck_measured.py; see docs/CROSSCHECK.md):

* members: ``elasticBeamColumn`` (Euler-Bernoulli) or, with shear deformation in a first-order
  combination, ``ElasticTimoshenkoBeam``; one ``geomTransf`` per element with ``vecxz`` = benchfea's
  local z (OpenSees: y = vecxz x x, z = x x y, which returns benchfea's y and z).
* truss-type members: one ``elasticBeamColumn`` whose bending and torsion stiffness is multiplied by
  1e-15 (benchfea transmits only the axial force); never subdivided.
* bending connections (local ry, rz: hinge or rotational spring, either end, any orientation, with or
  without shear deformation): EXACT decomposition into elements in parallel between the same nodes --
  a rigidly connected core plus partners carrying the ``-releasey/-releasez`` flag at one end
  (``fef.bending_split``).  Released members never receive ``eleLoad`` (the flags condense the
  fixed-end forces of full-span uniform loads only).
* axial / shear / torsion connections: duplicate node at the flexible end, ``equalDOF`` for the rigid
  directions and a ``zeroLength`` spring (``-orient`` = member axes) for partial ones.  ``equalDOF``
  acts on global DOFs and its retained node must not itself be constrained, so this is exact when the
  member axes are parallel to the global axes (or the whole translation / rotation triple is rigid),
  the end has no rigid offset and the joint is not a constraint slave; otherwise the rigid directions
  get penalty springs (``penalty`` x member stiffness) and the feature is "approximate".
* rigid end offsets: ``-jntOffset`` on ``elasticBeamColumn``; ``ElasticTimoshenkoBeam`` mishandles
  ``-jntOffset``, so those elements (and node chains) use a node at the flexible end tied to the joint
  by ``rigidLink('beam')``.
* rigid links / diaphragms: ``rigidLink('beam')`` / ``rigidDiaphragm`` with the ``Transformation``
  constraint handler (exact); a diaphragm slave outside the master's plane is "unsupported".
* member loads: uniform (full or partial) line loads and interior point forces on unreleased
  Euler-Bernoulli elements use ``eleLoad`` (``native_loads=True``); linearly varying loads, point
  moments and every load on a Timoshenko or released element become exact equivalent nodal loads with
  the fixed-end forces of :mod:`benchfea.crosscheck.fef`, which are added back to the element end forces.
* second order: ``geomTransf('PDelta')`` (chord-rotation geometric stiffness only, no P-delta inside
  an element), Newton iterations on the current axial forces, every frame member divided into
  ``subdivisions`` elements.  The discretisation error is O(1/n^2) provided the axial force the
  transformation uses is the mean axial force of the element.  Measured: with an AXIAL component in
  an ``eleLoad`` the transformation takes the axial force of element end j (the basic force without
  the element-load part), which is off by half the load on the element and lowers the convergence
  to O(1/n).  In a second-order combination the axial component of every member load is therefore
  applied as equivalent nodal loads (the axial force of each element is then constant and equal to
  the mean of the exact one for a uniform load); the transverse components stay native.  The error
  is bounded by the runner from three subdivision levels (``runner.estimate_from_levels``), with no
  order assumed beyond "at least first".  Member end forces are reported along the UNDEFORMED local
  axes (forces from the global element forces, moments from the element ends).
* linear solver: ``UmfPack`` (measured: far less round-off with penalty springs than the dense
  solvers, but it does not report a singular matrix -- the self-check below does).
* reactions are NOT read from ``nodeReaction`` (measured: zero at a fixed node that retains a
  constraint); they are computed from joint equilibrium of the element end forces returned by
  OpenSees, the applied loads and the known constraint topology.  The residual at the free DOFs, the
  violation of the kinematic constraints (from the displacements) and the sanity of the displacements
  are a self-check: a run that fails it is reported as failed.
"""
from __future__ import annotations

import math

import numpy as np

from .fef import bending_split, clip_loads, fixed_end_forces
from .fidelity import APPROXIMATE, EXACT, UNSUPPORTED, Fidelity
from .foreign import ForeignResult, Translation, refinement_levels
from .prepare import Prepared, prepare

SOLVER = "opensees"
TINY = 1e-15                           # relative stiffness of a direction an element must not carry
NO_SHEAR_AREA_FACTOR = 1e10            # shear area standing for "no shear deformation in this direction"
DEFAULT_SUBDIVISIONS = (8, 16, 32)     # three levels, ratio 2; the last is the delivered run (second order)
SUBDIVISION_ORDERS = (1.0, 2.0)        # convergence orders the error bound allows for: at least first (assumed
                                       # for the future levels), at most second (chord-rotation geometric stiffness)
DEFAULT_PENALTY = (1e5, 1e6)           # coarse, fine

__all__ = ["translate", "execute", "SOLVER"]

_LINEAR_NOTE = (" (second order: a linear kinematic relation, like benchfea's -- the second-order effect of the "
                "forces the constraint transmits is in neither)")


def _aligned(R, tol=1e-12) -> bool:
    return bool(np.all((np.abs(R) < tol) | (np.abs(np.abs(R) - 1.0) < tol)))


class _Builder:
    def __init__(self, p: Prepared):
        self.p = p
        self.model_cmds: list = []
        self.n_tag = len(p.node_ids)                # structural node k has tag k + 1
        self.e_tag = 0
        self.mat: dict[float, int] = {}
        self.xyz: dict[int, np.ndarray] = {k + 1: p.xyz[k] for k in range(len(p.node_ids))}
        self.loads: dict[int, np.ndarray] = {}
        self.ele_loads: list = []
        self.elements: list = []                    # (tag, kind, node_i, node_j)
        self.constraints: list = []                 # dict(kind, r, c, dofs, arm, normal)
        self.fix: dict[int, np.ndarray] = {}

    def tag(self, nid: str) -> int:
        return self.p.index(nid) + 1

    def node(self, xyz) -> int:
        self.n_tag += 1
        self.xyz[self.n_tag] = np.asarray(xyz, float)
        return self.n_tag

    def element(self, kind, ni, nj) -> int:
        self.e_tag += 1
        self.elements.append((self.e_tag, kind, ni, nj))
        return self.e_tag

    def material(self, k: float) -> int:
        if k not in self.mat:
            self.mat[k] = len(self.mat) + 1
            self.model_cmds.append(["uniaxialMaterial", "Elastic", self.mat[k], float(k)])
        return self.mat[k]

    def add_load(self, tag, force, moment):
        v = self.loads.setdefault(tag, np.zeros(6))
        v[:3] += force
        v[3:] += moment

    def offset_node(self, N, pos, cmds):
        """Node at ``pos`` rigidly attached to node N (rigidLink 'beam').  The retained node of an
        OpenSees constraint must not itself be constrained, so for a slave of a rigid link the new node
        is attached to that link's master (the same rigid body).  Returns (tag, exact)."""
        t = self.node(pos)
        cmds.append(["node", t, *map(float, pos)])
        owner = [c for c in self.constraints if c["c"] == N]
        exact = True
        r = N
        if owner:
            if all(c["kind"] == "rigid" for c in owner):
                r = owner[0]["r"]
            else:
                exact = False
        cmds.append(["rigidLink", "beam", r, t])
        self.constrain("rigid", r, t, range(6))
        return t, exact

    def constrain(self, kind, r, c, dofs, normal=None):
        self.constraints.append({"kind": kind, "r": r, "c": c, "dofs": list(dofs),
                                 "arm": (self.xyz[c] - self.xyz[r]).tolist(), "normal": normal})


def translate(model, combination: str, *, subdivisions: int | None = None, penalty: float | None = None,
              native_loads: bool = True, tau_b: dict | None = None, extra_nodal_loads: dict | None = None,
              system: str | None = None, newton_tolerance: float = 1e-10, max_iterations: int = 50) -> Translation:
    """Build the OpenSees job of one combination.

    ``subdivisions``: elements per frame member (default 1 first order, 16 second order);
    ``penalty``: stiffness multiplier of penalty springs (default 1e6);
    ``native_loads``: use ``eleLoad`` where OpenSees is exact (otherwise every member load becomes
    equivalent nodal loads); ``tau_b`` / ``extra_nodal_loads``: see :func:`prepare`.
    """
    p = prepare(model, combination, tau_b=tau_b, extra_nodal_loads=extra_nodal_loads)
    fid = Fidelity(SOLVER)
    for feature, level, note in p.issues:
        fid.add(feature, level, note)
    b = _Builder(p)
    cmds = b.model_cmds
    n_sub = int(subdivisions) if subdivisions else (DEFAULT_SUBDIVISIONS[-1] if p.second_order else 1)
    pen = float(penalty) if penalty else DEFAULT_PENALTY[1]
    controls = {}
    transf_name = "PDelta" if p.second_order else "Linear"
    size = float(np.linalg.norm(np.ptp(p.xyz, axis=0))) or 1.0

    for k in range(len(p.node_ids)):
        cmds.append(["node", k + 1, *map(float, p.xyz[k])])

    # ------------------------------------------------------------------ supports
    for k, nid in enumerate(p.node_ids):
        flags = p.fixed[k].copy()
        flags[p.inactive.get(nid, [])] = True
        if flags.any():
            b.fix[k + 1] = flags
        if p.spring[k].any():
            fid.add("supports.spring", EXACT, "zeroLength spring to a fixed twin node")
            twin = b.node(p.xyz[k])
            cmds.append(["node", twin, *map(float, p.xyz[k])])
            b.fix[twin] = np.ones(6, bool)
            dirs = [d for d in range(6) if p.spring[k, d] != 0.0]
            e = b.element("support", k + 1, twin)
            cmds.append(["element", "zeroLength", e, k + 1, twin, "-mat", *[b.material(p.spring[k, d]) for d in dirs],
                         "-dir", *[d + 1 for d in dirs]])
    if p.prescribed.any():
        fid.add("loads.supportDisplacement", EXACT, "sp constraint in the load pattern")
    if p.load.any():
        fid.add("loads.nodal", EXACT)
    for k in range(len(p.node_ids)):
        if p.load[k].any():
            b.add_load(k + 1, p.load[k, :3], p.load[k, 3:])

    # ------------------------------------------------------------------ constraints
    for c in p.constraints:
        m_tag = b.tag(c.master)
        if c.type == "rigidLink":
            fid.add("constraints.rigidLink", EXACT, "rigidLink('beam'), Transformation handler" + _LINEAR_NOTE * p.second_order)
            for s in c.slaves:
                cmds.append(["rigidLink", "beam", m_tag, b.tag(s)])
                b.constrain("rigid", m_tag, b.tag(s), range(6))
        else:
            fid.add("constraints.rigidDiaphragm", EXACT,
                    "rigidDiaphragm, Transformation handler" + _LINEAR_NOTE * p.second_order)
            a = "XYZ".index(c.normal)
            dofs = [d for d in range(3) if d != a] + [3 + a]
            direct = []
            for s in c.slaves:
                s_tag = b.tag(s)
                if abs(float((b.xyz[s_tag] - b.xyz[m_tag])[a])) > 1e-9 * size:
                    fid.add("constraints.rigidDiaphragm", UNSUPPORTED,
                            f"slave {s!r} is not in the plane of the master: benchfea's rule (in-plane part of the "
                            "offset only) is not a rigid-body constraint and OpenSees' rigidDiaphragm is defined for "
                            "nodes in one plane; the slave is tied with its in-plane offset by OpenSees' own rule")
                direct.append(s_tag)
                b.constrain("diaphragm", m_tag, s_tag, dofs, a)
            cmds.append(["rigidDiaphragm", a + 1, m_tag, *direct])

    # ------------------------------------------------------------------ members
    members_out = {}
    for mid, m in p.members.items():
        timo = (m.GAy is not None or m.GAz is not None) and not m.truss
        if timo and p.second_order:
            fid.add("shearDeformation", UNSUPPORTED,
                    "OpenSees 3.8 ElasticTimoshenkoBeam ignores the geometric nonlinearity of the transformation "
                    "(measured); the member is modelled WITHOUT shear deformation")
            timo = False
        elif timo:
            fid.add("shearDeformation", EXACT, "ElasticTimoshenkoBeam (first order)")
        fid.add("members.truss" if m.truss else "members.frame", EXACT,
                "elasticBeamColumn with bending and torsion stiffness x 1e-12" if m.truss else "")
        if m.offset_i.any() or m.offset_j.any():
            if p.second_order:
                fid.add("offsets.secondOrder", APPROXIMATE,
                        "the geometric stiffness of the rigid offset zones is the solver's own (-jntOffset of the "
                        "PDelta transformation / linear rigidLink); no reference, no error estimate")
            fid.add("offsets", EXACT, "-jntOffset (Euler-Bernoulli elements) / rigidLink('beam') to a node at the "
                                      "flexible end (Timoshenko elements, node chains)")

        # ---- end connections.  Bending rotations (local ry, rz): parallel elements (below).
        #      Axial / shear / torsion connections: node chain at the flexible end.
        ends, bend = [], []
        for end, (nid, off, k_end, pos) in enumerate(((m.i, m.offset_i, m.k_i, m.start),
                                                      (m.j, m.offset_j, m.k_j, m.end))):
            N = b.tag(nid)
            if m.truss:
                ends.append((N, off))
                bend.append((0.0, 0.0))
                continue
            bend.append((k_end[4], k_end[5]))
            for k in k_end[4:]:
                if not math.isinf(k):
                    fid.add("releases.partial" if k != 0.0 else "releases.full", EXACT,
                            "bending connection: parallel elements with release flags (exact decomposition)")
            if all(math.isinf(k) for k in k_end[:4]):
                if timo and off.any():
                    # measured: ElasticTimoshenkoBeam with -jntOffset shortens the element but drops the rigid
                    # arms (equilibrium violated) in OpenSees 3.8
                    Pn, exact = b.offset_node(N, pos, cmds)
                    if not exact:
                        fid.add("offsets", UNSUPPORTED, f"member {mid!r}: rigid offset of a shear-deformable "
                                "member at a diaphragm slave node cannot be expressed with OpenSees constraints")
                    ends.append((Pn, np.zeros(3)))
                else:
                    ends.append((N, off))
                continue
            feature = "releases.partial" if any(k != 0.0 and not math.isinf(k) for k in k_end[:4]) \
                else "releases.full"
            if off.any():
                Pn, exact = b.offset_node(N, pos, cmds)
                if not exact:
                    fid.add("offsets", UNSUPPORTED, f"member {mid!r}: rigid offset at a diaphragm slave node "
                            "combined with a node chain cannot be expressed with OpenSees constraints")
            else:
                Pn = N
            D = b.node(pos)
            cmds.append(["node", D, *map(float, pos)])
            # equalDOF needs global DOF directions and a retained node that is not itself constrained
            # (measured: a constrained node used as retained node gives wrong results with the
            # Transformation handler); otherwise penalty springs
            direct_tie = Pn == N and not any(c["c"] == N for c in b.constraints)
            aligned = _aligned(m.R) and direct_tie
            kc = list(k_end[:4]) + [math.inf, math.inf]
            mats, dirs, tied, approx = [], [], [], False
            for triple in ((0, 1, 2), (3, 4, 5)):
                if direct_tie and all(math.isinf(kc[a]) for a in triple):
                    tied += list(triple)
                    continue
                for a in triple:
                    k = kc[a]
                    if math.isinf(k):
                        if aligned:
                            tied.append(3 * (a // 3) + int(np.argmax(np.abs(m.R[a % 3]))))
                            continue
                        approx = True
                        k = pen * (m.E * m.A / m.L if a < 3 else max(m.EIy, m.EIz) / m.L)
                    if k != 0.0:
                        mats.append(b.material(k))
                        dirs.append(a + 1)
            if tied:
                tied = sorted(tied)
                cmds.append(["equalDOF", Pn, D, *[d + 1 for d in tied]])
                b.constrain("equal", Pn, D, tied)
            if mats:
                e = b.element("spring", Pn, D)
                cmds.append(["element", "zeroLength", e, Pn, D, "-mat", *mats, "-dir", *dirs,
                             "-orient", *map(float, m.R[0]), *map(float, m.R[1])])
            if approx:
                fid.add(feature, APPROXIMATE, "axial / shear / torsion connection of a member whose axes are not "
                        "parallel to the global axes, or at a rigid offset / constrained node: the rigid directions "
                        "use penalty springs", control="penalty")
                controls["penalty"] = {"coarse": pen / 10.0, "fine": pen, "order": None}
            else:
                fid.add(feature, EXACT, "axial / shear / torsion connection: duplicate node with equalDOF / "
                                        "zeroLength spring at the flexible end")
            ends.append((D, np.zeros(3)))

        # ---- sub-elements along the flexible length
        bows = m.bows(p.second_order)
        if bows and not (m.EIy > 0.0 and m.EIz > 0.0):
            fid.add("members.truss", UNSUPPORTED,
                    f"second order: the truss-type member {mid!r} carries an axial member load (it is a pin-ended "
                    "beam-column in benchfea) but its section has no bending stiffness to model it with")
            bows = False
        n = n_sub if (bows or not m.truss) else 1
        if bows:
            fid.add("secondOrder", APPROXIMATE,
                    "truss-type member with an axial member load: pin-ended beam-column, subdivided like a frame "
                    "member (bending hinges at its two ends, torsion stiffness x 1e-15)", control="subdivisions")
            controls["subdivisions"] = refinement_levels(n_sub, SUBDIVISION_ORDERS,
                                                         DEFAULT_SUBDIVISIONS if subdivisions is None else None)
        elif p.second_order and not m.truss:
            fid.add("secondOrder", APPROXIMATE,
                    "PDelta transformation (chord-rotation geometric stiffness, no P-delta inside an element), "
                    "Newton iterations on the axial forces; error bounded from three subdivision levels",
                    control="subdivisions")
            controls["subdivisions"] = refinement_levels(n_sub, SUBDIVISION_ORDERS,
                                                         DEFAULT_SUBDIVISIONS if subdivisions is None else None)
        elif p.second_order:
            fid.add("secondOrder", EXACT, "truss-type member with a constant axial force: the chord-rotation "
                                          "geometric stiffness is the complete geometric stiffness of a pin-ended bar")
        chain = [ends[0][0]]
        for s in range(1, n):
            pos = m.start + (m.end - m.start) * (s / n)
            t = b.node(pos)
            cmds.append(["node", t, *map(float, pos)])
            chain.append(t)
        chain.append(ends[1][0])
        Ls = m.L / n
        zvec = [float(v) for v in m.R[2]]
        Iy_eff, Iz_eff = m.EIy / m.E, m.EIz / m.E
        if m.truss:                               # a truss section may come without bending properties
            Iy_eff = Iy_eff or m.A ** 2 / 12.0
            Iz_eff = Iz_eff or m.A ** 2 / 12.0
        J_full = m.J or m.A ** 2 / 6.0
        # (a pin-ended bar gets no bending stiffness from ``bending_split`` with both ends hinged; a bowing
        # truss member is a chain whose first and last element carry one hinge each)
        gay, gaz = (m.GAy, m.GAz) if timo else (None, None)
        groups, fefs = [], []
        for s in range(n):
            first, last = s == 0, s == n - 1
            oi = ends[0][1] if first else np.zeros(3)
            oj = ends[1][1] if last else np.zeros(3)
            ky = (bend[0][0] if first else math.inf, bend[1][0] if last else math.inf)
            kz = (bend[0][1] if first else math.inf, bend[1][1] if last else math.inf)
            released = not all(math.isinf(k) for k in ky + kz)
            # plane x-z (about y): EIy with shear along z;  plane x-y (about z): EIz with shear along y
            ay, by, cy = bending_split(Ls, m.E * Iy_eff, gaz, *ky)
            az, bz, cz = bending_split(Ls, m.E * Iz_eff, gay, *kz)
            ni, nj = chain[s], chain[s + 1]
            tags = []

            def transf():
                e = b.element("beam", ni, nj)
                t_args = ["geomTransf", transf_name, e, *zvec]
                if oi.any() or oj.any():
                    t_args += ["-jntOffset", *map(float, oi), *map(float, oj)]
                cmds.append(t_args)
                tags.append(e)
                return e

            e = transf()
            Iy_c, Iz_c = Iy_eff * max(ay, TINY), Iz_eff * max(az, TINY)
            # a truss-type member transmits no torque: twist released at end j (docs/THEORY.md 1.5), i.e. no
            # torsion stiffness in its last element; the other elements of a chain hold the twist of the
            # internal nodes, which no torque can reach
            J = J_full * (TINY if (m.truss and last) else 1.0)
            if timo:
                Avy = (gay / m.G) * (az if az > TINY else 1.0) if gay else NO_SHEAR_AREA_FACTOR * m.A
                Avz = (gaz / m.G) * (ay if ay > TINY else 1.0) if gaz else NO_SHEAR_AREA_FACTOR * m.A
                cmds.append(["element", "ElasticTimoshenkoBeam", e, ni, nj, m.E, m.G, m.A, J, Iy_c, Iz_c, Avy, Avz, e])
            else:
                cmds.append(["element", "elasticBeamColumn", e, ni, nj, m.A, m.E, m.G, J, Iy_c, Iz_c, e])
            # hinged partners (Euler-Bernoulli with the equivalent stiffness EI' of fef.bending_split)
            Iy_h = Iy_eff / (1.0 + 3.0 * m.E * Iy_eff / (gaz * Ls * Ls)) if gaz else Iy_eff
            Iz_h = Iz_eff / (1.0 + 3.0 * m.E * Iz_eff / (gay * Ls * Ls)) if gay else Iz_eff
            for coeff, flag, code, about_z in ((bz, "-releasez", 2, True), (cz, "-releasez", 1, True),
                                               (by, "-releasey", 2, False), (cy, "-releasey", 1, False)):
                if coeff > 1e-14:
                    e = transf()
                    iy = Iy_eff * TINY if about_z else Iy_h * coeff
                    iz = Iz_h * coeff if about_z else Iz_eff * TINY
                    cmds.append(["element", "elasticBeamColumn", e, ni, nj, m.A * TINY, m.E, m.G,
                                 (m.J or m.A ** 2 / 6.0) * TINY, iy, iz, e, flag, code])
            groups.append(tags)
            # ---- loads of this sub-element
            pts, ds = clip_loads(m.points, m.dists, s * Ls, (s + 1) * Ls, last, tol=1e-12 * m.L)
            native = native_loads and not timo and not released
            rest_p, rest_d = [], []
            # second order: the PDelta transformation takes the axial force of element end j, which does
            # not contain the axial part of an element load (measured) -> axial components go to the nodes
            split = p.second_order
            for kind, vec, x in pts:
                if native and kind == "force" and 1e-9 * Ls < x < Ls * (1 - 1e-9):
                    ax = float(vec[0]) if split else 0.0
                    if vec[1] != 0.0 or vec[2] != 0.0 or not ax:
                        b.ele_loads.append(["eleLoad", "-ele", tags[0], "-type", "-beamPoint", float(vec[1]),
                                            float(vec[2]), x / Ls, float(vec[0]) - ax])
                        fid.add("loads.memberPoint", EXACT, "eleLoad -beamPoint")
                    if ax:
                        rest_p.append((kind, np.array([ax, 0.0, 0.0]), x))
                        fid.add("loads.memberPoint", EXACT, "second order: axial component as equivalent nodal loads")
                else:
                    rest_p.append((kind, vec, x))
                    fid.add("loads.memberPoint", EXACT, "equivalent nodal loads from exact fixed-end forces")
            for w1, w2, x1, x2 in ds:
                if native and np.array_equal(w1, w2):
                    ax = float(w1[0]) if split else 0.0
                    if w1[1] != 0.0 or w1[2] != 0.0 or not ax:
                        args = ["eleLoad", "-ele", tags[0], "-type", "-beamUniform", float(w1[1]), float(w1[2]),
                                float(w1[0]) - ax]
                        if x1 > 1e-12 * Ls or x2 < Ls * (1 - 1e-12):
                            args += [x1 / Ls, x2 / Ls]
                        b.ele_loads.append(args)
                        fid.add("loads.memberDistributed", EXACT, "eleLoad -beamUniform")
                    if ax:
                        rest_d.append((np.array([ax, 0.0, 0.0]), np.array([ax, 0.0, 0.0]), x1, x2))
                        fid.add("loads.memberDistributed", EXACT,
                                "second order: axial component as equivalent nodal loads")
                else:
                    rest_d.append((w1, w2, x1, x2))
                    fid.add("loads.memberDistributed", EXACT, "equivalent nodal loads from exact fixed-end forces")
            f = np.zeros(12)
            if rest_p or rest_d:
                f = fixed_end_forces(Ls, m.E * Iy_eff, m.E * Iz_eff, gay, gaz, rest_p, rest_d, ky=ky, kz=kz)
                for blk, node, off in ((0, ni, oi), (6, nj, oj)):
                    F = -m.R.T @ f[blk:blk + 3]
                    M = -m.R.T @ f[blk + 3:blk + 6]
                    b.add_load(node, F, M + np.cross(off, F))
            fefs.append(f)
        members_out[mid] = {"first": groups[0], "last": groups[-1], "fef_i": fefs[0][:6].tolist(),
                            "fef_j": fefs[-1][6:].tolist()}

    for tag, flags in b.fix.items():
        cmds.append(["fix", tag, *[int(v) for v in flags]])

    pattern = [["load", tag, *map(float, v)] for tag, v in b.loads.items() if v.any()]
    pattern += b.ele_loads
    for k in range(len(p.node_ids)):
        for d in range(6):
            if p.fixed[k, d] and p.prescribed[k, d] != 0.0:
                pattern.append(["sp", k + 1, d + 1, float(p.prescribed[k, d])])
    job = {"solver": SOLVER, "model": cmds, "pattern": pattern, "secondOrder": p.second_order,
           "tolerance": newton_tolerance, "maxIterations": max_iterations,
           "system": system or "UmfPack",
           "nodes": list(range(1, b.n_tag + 1)), "elements": [e[0] for e in b.elements],
           "beams": [e[0] for e in b.elements if e[1] == "beam"]}
    options = {"subdivisions": n_sub, "penalty": pen, "native_loads": native_loads,
               "penalised": "penalty" in controls}

    def mapper(raw):
        return _map(p, b, members_out, fid, raw, options)

    return Translation(SOLVER, combination, job, fid, controls, options, mapper)


# ====================================================================================== post-processing
def _map(p: Prepared, b: _Builder, members_out, fid, raw, options) -> ForeignResult:
    disp = {int(k): np.array(v) for k, v in raw["disp"].items()}
    glob = {int(k): np.array(v) for k, v in raw["eleForce"].items()}
    loc = {int(k): np.array(v) for k, v in raw["localForce"].items()}
    n = len(p.node_ids)
    U = np.array([disp[k + 1] for k in range(n)])
    scale = float(np.abs(U).max()) or 1.0

    # kinematic constraints as OpenSees enforced them (self-check)
    violation = 0.0
    for c in b.constraints:
        ur, uc, arm = disp[c["r"]], disp[c["c"]], np.array(c["arm"])
        if c["kind"] == "rigid":
            expect = np.r_[ur[:3] + np.cross(ur[3:], arm), ur[3:]]
        elif c["kind"] == "equal":
            expect = ur
        else:
            a = c["normal"]
            arm_in = arm.copy()
            arm_in[a] = 0.0
            th = np.zeros(3)
            th[a] = ur[3 + a]
            expect = np.r_[ur[:3] + np.cross(th, arm_in), ur[3:]]
        violation = max(violation, float(np.abs((uc - expect)[c["dofs"]]).max()) / scale)

    # joint equilibrium: what each node needs from supports and constraints
    need = {t: -b.loads.get(t, np.zeros(6)).astype(float) for t in b.xyz}
    for tag, kind, ni, nj in b.elements:
        if kind == "support":
            continue
        need[ni] = need[ni] + glob[tag][:6]
        need[nj] = need[nj] + glob[tag][6:]
    pending = list(b.constraints)
    while pending:
        retained = {c["r"] for c in pending}
        leaves = [c for c in pending if c["c"] not in retained]
        if not leaves:                                    # pragma: no cover - translator never builds cycles
            raise RuntimeError("cyclic constraints")
        for c in leaves:
            d, r = c["dofs"], c["r"]
            v = np.zeros(6)
            v[d] = need[c["c"]][d]
            need[c["c"]][d] = 0.0
            arm = np.array(c["arm"])
            if c["kind"] == "equal":
                need[r] += v
            elif c["kind"] == "rigid":
                need[r][:3] += v[:3]
                need[r][3:] += v[3:] + np.cross(arm, v[:3])
            else:                                         # diaphragm: in-plane forces, moment about the normal
                a = c["normal"]
                arm[a] = 0.0
                need[r][:3] += v[:3]
                need[r][3 + a] += v[3 + a] + np.cross(arm, v[:3])[a]
        pending = [c for c in pending if c not in leaves]
    reactions = np.zeros((n, 6))
    # scales of the residual: largest element end force; largest end moment, but at least
    # 1e-3 x force scale x model size (a model without bending has only round-off moments)
    size = float(np.linalg.norm(np.ptp(p.xyz, axis=0))) or 1.0
    ef = np.array([glob[tag] for tag, kind, _, _ in b.elements if kind != "support"]).reshape(-1, 12)
    fmax = max(float(np.abs(ef[:, [0, 1, 2, 6, 7, 8]]).max()), 1e-300)
    mmax = max(float(np.abs(ef[:, [3, 4, 5, 9, 10, 11]]).max()), 1e-3 * fmax * size)
    residual = 0.0
    for t, v in need.items():
        mask = b.fix.get(t, np.zeros(6, bool)).copy()
        if t <= n:
            mask |= p.spring[t - 1] != 0.0
            reactions[t - 1] = np.where(mask, v, 0.0)
            reactions[t - 1, p.inactive.get(p.node_ids[t - 1], [])] = 0.0
        if (~mask).any():
            rel = np.abs(v) / np.r_[[fmax] * 3, [mmax] * 3]
            residual = max(residual, float(rel[~mask].max()))

    end_forces = {}
    for mid, out in members_out.items():
        m = p.members[mid]
        fi = sum(loc[t][:6] for t in out["first"])
        fj = sum(loc[t][6:] for t in out["last"])
        if p.second_order:
            # force components along the UNDEFORMED local axes (OpenSees localForce shears are normal
            # to the displaced chord); moments are those of the element ends
            fi[:3] = m.R @ sum(glob[t][:3] for t in out["first"])
            fj[:3] = m.R @ sum(glob[t][6:9] for t in out["last"])
        end_forces[mid] = np.r_[fi + np.array(out["fef_i"]), fj + np.array(out["fef_j"])] / p.alpha

    info = {"converged": raw.get("converged", True), "iterations": raw.get("iterations"),
            "constraintViolation": violation, "jointResidual": residual, "alpha": p.alpha,
            "version": raw.get("version"), **options}
    limit = max(1e-6, 1e-11 * options["penalty"]) if options.get("penalised") else 1e-6
    # UmfPack does not report a singular matrix (measured): a mechanism shows as non-finite or absurd
    # displacements and a joint residual of order one
    sane = bool(np.all(np.isfinite(U))) and float(np.abs(U[:, :3]).max()) < 1e6 * size
    ok = bool(raw.get("converged", True)) and sane and violation < 1e-8 and residual < limit
    error = "" if ok else ("singular or ill-conditioned system (mechanism in the foreign model?): " * (not sane)
                           + f"self-check failed: constraint violation {violation:.2e}, joint equilibrium "
                           f"residual {residual:.2e}, converged {raw.get('converged', True)}")
    if not ok:
        return ForeignResult(SOLVER, p.combination, False, fid, error=error, info=info)
    return ForeignResult(SOLVER, p.combination, True, fid, U / p.alpha, reactions / p.alpha, end_forces, "", info)


# ====================================================================================== worker
def execute(job: dict) -> dict:
    """Run a job with OpenSeesPy (call this in a subprocess: see module docstring)."""
    import openseespy.opensees as ops

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for c in job["model"]:
        getattr(ops, c[0])(*c[1:])
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for c in job["pattern"]:
        getattr(ops, c[0])(*c[1:])
    ops.system(job["system"])
    ops.numberer("RCM")
    ops.constraints("Transformation")
    if job["secondOrder"]:
        ops.test("NormDispIncr", job["tolerance"], job["maxIterations"])
        ops.algorithm("Newton")
    else:
        ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    rc = ops.analyze(1)
    out = {"ok": rc == 0, "converged": rc == 0, "version": str(ops.version())}
    if rc != 0:
        out["error"] = f"OpenSees analyze() returned {rc} (singular system or no convergence; see stderr)"
        return out
    if job["secondOrder"]:
        try:
            out["iterations"] = int(ops.testIter())
        except Exception:                                   # pragma: no cover - optional diagnostic
            out["iterations"] = None
    out["disp"] = {str(t): list(ops.nodeDisp(t)) for t in job["nodes"]}
    out["eleForce"] = {str(t): list(ops.eleForce(t)) for t in job["elements"]}
    out["localForce"] = {str(t): list(ops.eleResponse(t, "localForce")) for t in job["beams"]}
    if job.get("nodeReaction"):                             # diagnostic only: the translator never uses it
        ops.reactions()
        out["nodeReaction"] = {str(t): list(ops.nodeReaction(t)) for t in job["nodeReaction"]}
    return out
