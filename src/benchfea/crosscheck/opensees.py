"""benchfea model + combination  ->  OpenSeesPy model (docs/CROSSCHECK.md has the full account).

``translate`` runs in the calling process and imports nothing from OpenSees: it produces a *job*, a
plain JSON list of ``openseespy.opensees`` calls plus the outputs to read.  ``execute(job)`` performs
the calls; it is meant to run in a subprocess (``benchfea.crosscheck.runner``) because OpenSees
terminates the interpreter on some input errors.  The mapping of the raw output back to benchfea
conventions (``Translation.results``) again needs no OpenSees.

How each feature is represented (every statement about OpenSees 3.8 behaviour was measured, see
docs/CROSSCHECK.md):

* members: ``elasticBeamColumn`` (Euler-Bernoulli) or, with shear deformation in a first-order
  combination, ``ElasticTimoshenkoBeam``; one ``geomTransf`` per element with ``vecxz`` = benchfea's
  local z (OpenSees: y = vecxz x x, z = x x y, which returns benchfea's y and z).
* truss-type members: ``elasticBeamColumn`` with both bending moments released at both ends and the
  torsion constant multiplied by 1e-12 (benchfea transmits no torque through a truss member).
* rigid end offsets: ``-jntOffset`` of the transformation; where the end also needs a node chain (see
  releases) an explicit node at the flexible end tied to the joint by ``rigidLink('beam')``.
* releases: full release of My / Mz only -> element flags ``-releasey/-releasez`` (exact, any
  orientation); anything else (partial springs, axial/shear/torsion releases, any release on a
  Timoshenko element) -> duplicate node at the flexible end connected by ``equalDOF`` for rigid
  directions and a ``zeroLength`` spring (``-orient`` = member axes) for partial ones.  ``equalDOF``
  works on global DOFs, so this is exact when the member axes are parallel to the global axes or the
  whole translation / rotation triple is rigid; otherwise the rigid directions of a partly released
  triple get a penalty spring (``penalty`` x member stiffness) and the feature is "approximate".
* rigid links / diaphragms: ``rigidLink('beam')`` / ``rigidDiaphragm`` with the ``Transformation``
  constraint handler (exact); a diaphragm slave outside the master's plane is "unsupported".
* member loads: uniform (full or partial) line loads and interior point forces on Euler-Bernoulli
  elements use ``eleLoad`` (``native_loads=True``); linearly varying loads, point moments and every
  load on a Timoshenko element are turned into exact equivalent nodal loads with the fixed-end forces
  of :mod:`benchfea.crosscheck.fef`, which are added back to the element end forces.
* second order: ``geomTransf('PDelta')`` (chord-rotation geometric stiffness only, no P-delta inside
  an element), Newton iterations on the current axial forces, every frame member divided into
  ``subdivisions`` elements; the discretisation error is O(1/n^2) and is estimated by the runner from
  two subdivision levels.
* reactions are NOT read from ``nodeReaction`` (measured: zero at a fixed node that retains a
  constraint); they are computed from joint equilibrium of the element end forces returned by
  OpenSees, the applied loads and the known constraint topology, and the residual at the free DOFs
  is reported as a self-check together with the violation of the kinematic constraints.
"""
from __future__ import annotations

import math

import numpy as np

from .fef import bending_split, clip_loads, fixed_end_forces
from .fidelity import APPROXIMATE, EXACT, UNSUPPORTED, Fidelity
from .foreign import ForeignResult, Translation
from .prepare import Prepared, prepare

SOLVER = "opensees"
TINY = 1e-15                           # relative stiffness of a direction an element must not carry
NO_SHEAR_AREA_FACTOR = 1e10            # shear area standing for "no shear deformation in this direction"
DEFAULT_SUBDIVISIONS = (8, 16)         # coarse, fine (second order)
DEFAULT_PENALTY = (1e5, 1e6)           # coarse, fine

__all__ = ["translate", "execute", "SOLVER"]


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
    n_sub = int(subdivisions) if subdivisions else (DEFAULT_SUBDIVISIONS[1] if p.second_order else 1)
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
            fid.add("constraints.rigidLink", EXACT, "rigidLink('beam'), Transformation handler")
            for s in c.slaves:
                cmds.append(["rigidLink", "beam", m_tag, b.tag(s)])
                b.constrain("rigid", m_tag, b.tag(s), range(6))
        else:
            fid.add("constraints.rigidDiaphragm", EXACT, "rigidDiaphragm, Transformation handler")
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
                fid.add("offsets", APPROXIMATE, "second order: no geometric stiffness is attributed to the rigid "
                                                "offset zones by this translation; not compared with a reference")
            else:
                fid.add("offsets", EXACT, "-jntOffset (Euler-Bernoulli elements) / rigidLink('beam') to a node at "
                                          "the flexible end (Timoshenko elements, node chains)")

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
                    # measured: ElasticTimoshenkoBeam does not honour -jntOffset in OpenSees 3.8
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
        n = 1 if m.truss else n_sub
        if p.second_order and not m.truss:
            fid.add("secondOrder", APPROXIMATE,
                    "PDelta transformation (chord-rotation geometric stiffness, no P-delta inside an element), "
                    "Newton iterations on the axial forces; error O(1/n^2) in the number of elements per member",
                    control="subdivisions")
            controls["subdivisions"] = {"coarse": DEFAULT_SUBDIVISIONS[0] if subdivisions is None
                                        else max(1, n_sub // 2), "fine": n_sub, "order": 2}
        elif p.second_order:
            fid.add("secondOrder", EXACT, "truss-type member: chord-rotation geometric stiffness is the complete "
                                          "geometric stiffness of a pin-ended bar")
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
        J = (m.J or m.A ** 2 / 6.0) * (TINY if m.truss else 1.0)
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
            for kind, vec, x in pts:
                if native and kind == "force" and 1e-9 * Ls < x < Ls * (1 - 1e-9):
                    b.ele_loads.append(["eleLoad", "-ele", tags[0], "-type", "-beamPoint", float(vec[1]),
                                        float(vec[2]), x / Ls, float(vec[0])])
                    fid.add("loads.memberPoint", EXACT, "eleLoad -beamPoint")
                else:
                    rest_p.append((kind, vec, x))
                    fid.add("loads.memberPoint", EXACT, "equivalent nodal loads from exact fixed-end forces")
            for w1, w2, x1, x2 in ds:
                if native and np.array_equal(w1, w2):
                    args = ["eleLoad", "-ele", tags[0], "-type", "-beamUniform", float(w1[1]), float(w1[2]),
                            float(w1[0])]
                    if x1 > 1e-12 * Ls or x2 < Ls * (1 - 1e-12):
                        args += [x1 / Ls, x2 / Ls]
                    b.ele_loads.append(args)
                    fid.add("loads.memberDistributed", EXACT, "eleLoad -beamUniform")
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
    fmax = max(max(float(np.abs(v[:3]).max()) for v in need.values()), 1e-300)
    mmax = max(max(float(np.abs(v[3:]).max()) for v in need.values()), 1e-300)
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
    ok = bool(raw.get("converged", True)) and violation < 1e-8 and residual < limit
    error = "" if ok else (f"self-check failed: constraint violation {violation:.2e}, joint equilibrium "
                           f"residual {residual:.2e}, converged {raw.get('converged', True)}")
    return ForeignResult(SOLVER, p.combination, ok, fid, U / p.alpha, reactions / p.alpha, end_forces, error, info)


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
    return out
