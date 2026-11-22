"""benchfea model + combination  ->  PyNite (PyNiteFEA 3.x) model (docs/CROSSCHECK.md has the full account).

``translate`` imports nothing from PyNite and produces a JSON *job* (calls on ``FEModel3D`` plus the
member orientations and the outputs to read); ``execute(job)`` runs it (in a subprocess through
``benchfea.crosscheck.runner``); ``Translation.results`` maps the raw output back.

How each feature is represented:

* local axes: PyNite derives member axes from a Y-up rule of its own.  The worker reads the triad
  PyNite gives the unrotated member (``member.T()``), computes the ``rotation`` that turns its y axis
  into benchfea's, sets it and re-reads the triad: the run fails unless every member's (and
  sub-member's) triad equals benchfea's ``R`` to 1e-9.  No coordinate permutation is used.
* members: Euler-Bernoulli only -> shear deformation is "unsupported" (the model is then built
  without it and must not be used as a cross-check).
* truss-type members: both bending moments released at both ends, torsion released at end j.
* releases: full releases -> PyNite's boolean releases (any of the 12 end DOFs, exact first order).
  Partial ROTATIONAL springs about local y / z -> exact parallel members (``fef.bending_split``) with
  the member loads applied as exact equivalent nodal loads (``fef``).  Partial axial / shear / torsion
  springs -> "unsupported" (PyNite has no partial release).
* rigid end offsets, rigid links, rigid diaphragms: PyNite has no constraint and no offset facility.
  They are emulated by stiff members (``stiffness`` x the stiffest member of the model, uniform in
  the stub length) -> "approximate"; the runner solves at two stiffness levels and reports the
  difference as estimated error.  A diaphragm link keeps axial force, in-plane shear and the moment
  about the normal and releases the rest at the slave; slaves outside the master's plane are
  "unsupported".
* member loads: PyNite's own point loads (forces and moments) and linearly varying line loads in
  local axes, absolute positions from the flexible start.
* second order: ``analyze_PDelta`` = ONE solve with the consistent geometric stiffness built from the
  FIRST-ORDER axial forces (no iteration on N), members subdivided into explicit sub-members (PyNite's
  own splitting of a member at intermediate nodes mishandles linearly varying loads, measured).  Besides the
  discretisation error (two subdivision levels) the fidelity record carries two estimates computed
  from the run itself: the effect of not iterating N and of PyNite's axial term in Kg.
"""
from __future__ import annotations

import math

import numpy as np

from .fef import bending_split, clip_loads, fixed_end_forces
from .fidelity import APPROXIMATE, EXACT, UNSUPPORTED, Fidelity
from .foreign import ForeignResult, Translation
from .prepare import Prepared, prepare

SOLVER = "pynite"
TINY = 1e-15
DEFAULT_SUBDIVISIONS = (2, 4)          # coarse, fine (second order)
DEFAULT_STIFFNESS = (1e3, 1e4)         # coarse, fine multiplier of stiff emulation members (measured optimum:
                                       # the round-off error of PyNite's solve grows like 1e-9 x multiplier)
DOFS = ("DX", "DY", "DZ", "RX", "RY", "RZ")
FORCES = ("FX", "FY", "FZ", "MX", "MY", "MZ")
LOCAL = ("Fx", "Fy", "Fz")
LOCAL_M = ("Mx", "My", "Mz")

__all__ = ["translate", "execute", "finalise", "SOLVER"]


def translate(model, combination: str, *, subdivisions: int | None = None, stiffness: float | None = None,
              tau_b: dict | None = None, extra_nodal_loads: dict | None = None) -> Translation:
    """Build the PyNite job of one combination.  ``subdivisions``: sub-members per frame member
    (default 1 first order, 4 second order); ``stiffness``: multiplier of the stiff members that stand
    for rigid offsets / links / diaphragms (default 1e4)."""
    p = prepare(model, combination, tau_b=tau_b, extra_nodal_loads=extra_nodal_loads)
    fid = Fidelity(SOLVER)
    for feature, level, note in p.issues:
        fid.add(feature, level, note)
    fid.add("localAxes", EXACT, "rotation computed from PyNite's own triad and verified in the worker")
    n_sub = int(subdivisions) if subdivisions else (DEFAULT_SUBDIVISIONS[1] if p.second_order else 1)
    stiff = float(stiffness) if stiffness else DEFAULT_STIFFNESS[1]
    controls = {}
    cmds: list = []
    orient: list = []
    n = len(p.node_ids)
    xyz: dict[str, np.ndarray] = {f"N{k}": p.xyz[k] for k in range(n)}
    size = float(np.linalg.norm(np.ptp(p.xyz, axis=0))) or 1.0
    loads = {name: np.zeros(6) for name in xyz}
    segments: list = []                    # (member name, node_i, node_j)
    counter = [0]

    def name_of(nid):
        return f"N{p.index(nid)}"

    def new_node(pos):
        counter[0] += 1
        name = f"X{counter[0]}"
        xyz[name] = np.asarray(pos, float)
        loads[name] = np.zeros(6)
        return name

    # reference stiffness of the emulation members
    ms = list(p.members.values())
    L_ref = min(m.L for m in ms)
    E_ref, G_ref = max(m.E for m in ms), max(m.G for m in ms)
    A_ref = max(m.A for m in ms)
    I_ref = max(max(m.EIy, m.EIz) / m.E for m in ms)
    J_ref = max(max(m.J for m in ms), I_ref)

    def stiff_member(name, ni, nj, releases=None, y_axis=None):
        length = float(np.linalg.norm(xyz[nj] - xyz[ni]))
        r = length / L_ref
        cmds.append(["add_material", name, E_ref, G_ref, 0.3, 0.0])
        cmds.append(["add_section", name, stiff * A_ref * r, stiff * I_ref * max(r, r ** 3),
                     stiff * I_ref * max(r, r ** 3), stiff * J_ref * r])
        cmds.append(["add_member", name, ni, nj, name, name])
        if releases:
            cmds.append(["def_releases", name, *releases])
        if y_axis is not None:
            orient.append([name, [float(v) for v in y_axis]])
        segments.append((name, ni, nj))
        controls["stiffness"] = {"coarse": stiff / 10.0, "fine": stiff, "order": None}

    # ------------------------------------------------------------------ constraints
    for c in p.constraints:
        mn = name_of(c.master)
        for s in c.slaves:
            sn = name_of(s)
            d = xyz[sn] - xyz[mn]
            if float(np.linalg.norm(d)) <= 1e-9 * size:
                fid.add("constraints." + c.type, UNSUPPORTED, f"slave {s!r} coincides with the master: no stiff "
                                                               "member can be placed between them")
                continue
            if c.type == "rigidLink":
                fid.add("constraints.rigidLink", APPROXIMATE, "stiff member master -> slave", control="stiffness")
                stiff_member(f"L{len(segments)}", mn, sn)
            else:
                a = "XYZ".index(c.normal)
                if abs(float(d[a])) > 1e-9 * size:
                    fid.add("constraints.rigidDiaphragm", UNSUPPORTED,
                            f"slave {s!r} is not in the plane of the master (benchfea uses the in-plane part of "
                            "the offset only; a stiff member cannot do that)")
                    continue
                fid.add("constraints.rigidDiaphragm", APPROXIMATE,
                        "stiff member master -> slave, local y = diaphragm normal, out-of-plane shear, torsion and "
                        "out-of-plane bending released at the slave", control="stiffness")
                normal = np.eye(3)[a]
                #            Dxi   Dyi   Dzi   Rxi   Ryi   Rzi   Dxj   Dyj   Dzj   Rxj   Ryj   Rzj
                rel = [False] * 6 + [False, True, False, True, False, True]
                stiff_member(f"L{len(segments)}", mn, sn, rel, normal)

    # ------------------------------------------------------------------ members
    members_out = {}
    all_names: list = []
    for k_m, (mid, m) in enumerate(p.members.items()):
        base = f"M{k_m}"
        if (m.GAy is not None or m.GAz is not None) and not m.truss:
            fid.add("shearDeformation", UNSUPPORTED,
                    "PyNite members are Euler-Bernoulli (no shear area input); the member is modelled WITHOUT "
                    "shear deformation")
        fid.add("members.truss" if m.truss else "members.frame", EXACT,
                "moments released at both ends, torsion at end j" if m.truss else "")
        ends = []
        for nid, off, pos in ((m.i, m.offset_i, m.start), (m.j, m.offset_j, m.end)):
            N = name_of(nid)
            if off.any():
                fid.add("offsets", APPROXIMATE, "stiff member from the joint to the flexible end",
                        control="stiffness")
                if p.second_order:
                    fid.add("offsets", APPROXIMATE, "second order: the stiff stub takes part in PyNite's geometric "
                                                    "stiffness; not compared with a reference")
                P = new_node(pos)
                stiff_member(f"O{len(segments)}", N, P)
                ends.append(P)
            else:
                ends.append(N)
        k_i, k_j = (m.k_i, m.k_j) if not m.truss else ((math.inf,) * 3 + (math.inf, 0.0, 0.0),
                                                       (math.inf,) * 3 + (0.0, 0.0, 0.0))
        partial = [a for kk in (k_i, k_j) for a, k in enumerate(kk) if k != 0.0 and not math.isinf(k)]
        if any(a < 4 for a in partial):
            fid.add("releases.partial", UNSUPPORTED, "PyNite releases are boolean: no partial axial / shear / "
                                                     "torsion connection (treated as rigid here)")
        partial_bend = any(a >= 4 for a in partial)
        if not m.truss:
            for kk in (k_i, k_j):
                if any(k == 0.0 for k in kk):
                    fid.add("releases.full", EXACT, "def_releases")
        if partial_bend:
            fid.add("releases.partial", EXACT, "rotational springs about local y / z: parallel members (exact "
                    "decomposition) with member loads as exact equivalent nodal loads")
        # sub-division (second order): explicit sub-members between intermediate nodes
        nsub = 1 if (m.truss or not p.second_order) else n_sub
        if p.second_order and not m.truss:
            fid.add("secondOrder", APPROXIMATE,
                    "analyze_PDelta: consistent geometric stiffness, one solve with first-order axial forces; "
                    "discretisation error estimated from two subdivision levels", control="subdivisions")
            controls["subdivisions"] = {"coarse": DEFAULT_SUBDIVISIONS[0] if subdivisions is None
                                        else max(1, n_sub // 2), "fine": n_sub, "order": 2}
            if any(k != math.inf for k in k_i + k_j):
                fid.add("releases.full", APPROXIMATE, "second order: PyNite condenses the geometric stiffness of a "
                        "released member separately from the elastic stiffness (difference vanishes with "
                        "subdivision)", control="subdivisions")
        elif p.second_order:
            fid.add("secondOrder", APPROXIMATE, "truss-type member in analyze_PDelta (first-order axial force)",
                    control="subdivisions")
            controls.setdefault("subdivisions", {"coarse": DEFAULT_SUBDIVISIONS[0] if subdivisions is None
                                                 else max(1, n_sub // 2), "fine": n_sub, "order": 2})
        chain = [ends[0]] + [new_node(m.start + (m.end - m.start) * (s / nsub)) for s in range(1, nsub)] + [ends[1]]
        Ls = m.L / nsub
        Iy, Iz = m.EIy / m.E, m.EIz / m.E
        if m.truss:
            Iy, Iz = Iy or m.A ** 2 / 12.0, Iz or m.A ** 2 / 12.0
        J = m.J or m.A ** 2 / 6.0
        groups, fefs = [], []
        for s in range(nsub):
            first, last = s == 0, s == nsub - 1
            ni, nj = chain[s], chain[s + 1]
            ki = k_i if first else (math.inf,) * 6
            kj = k_j if last else (math.inf,) * 6
            names = []

            def add(A, iy, iz, j, rel):
                name = f"{base}s{s}p{len(names)}"
                names.append(name)
                all_names.append(name)
                cmds.append(["add_material", name, m.E, m.G, 0.3, 0.0])
                cmds.append(["add_section", name, A, iy, iz, j])
                cmds.append(["add_member", name, ni, nj, name, name])
                if any(rel):
                    cmds.append(["def_releases", name, *rel])
                orient.append([name, [float(v) for v in m.R[1]]])
                segments.append((name, ni, nj))

            boolean = [k == 0.0 for k in ki] + [k == 0.0 for k in kj]
            springs = any(k != 0.0 and not math.isinf(k) for k in ki[4:] + kj[4:])
            if not springs:
                add(m.A, Iy, Iz, J, boolean)
            else:
                # bending connections through parallel members; other (full) releases stay on the core
                ay, by, cy = bending_split(Ls, m.E * Iy, None, ki[4], kj[4])
                az, bz, cz = bending_split(Ls, m.E * Iz, None, ki[5], kj[5])
                core_rel = list(boolean)
                for idx in (4, 5, 10, 11):
                    core_rel[idx] = False
                add(m.A, Iy * max(ay, TINY), Iz * max(az, TINY), J, core_rel)
                for coeff, idx, about_z in ((bz, 11, True), (cz, 5, True), (by, 10, False), (cy, 4, False)):
                    if coeff > 1e-14:
                        rel = [False] * 12
                        rel[idx] = True
                        add(m.A * TINY, Iy * TINY if about_z else Iy * coeff, Iz * coeff if about_z else Iz * TINY,
                            J * TINY, rel)
            groups.append(names)
            # ---- loads of this sub-member (clipped here: PyNite's own splitting of a physical member
            #      mishandles linearly varying loads, measured)
            pts, ds = clip_loads(m.points, m.dists, s * Ls, (s + 1) * Ls, last, tol=1e-12 * m.L)
            fef = np.zeros(12)
            if springs and (pts or ds):
                fef = fixed_end_forces(Ls, m.E * Iy, m.E * Iz, None, None, pts, ds,
                                       ky=(ki[4], kj[4]), kz=(ki[5], kj[5]))
                for blk, node in ((0, ni), (6, nj)):
                    loads[node][:3] -= m.R.T @ fef[blk:blk + 3]
                    loads[node][3:] -= m.R.T @ fef[blk + 3:blk + 6]
                fid.add("loads.memberPoint" if pts else "loads.memberDistributed", EXACT,
                        "equivalent nodal loads from exact fixed-end forces (member with rotational springs)")
            else:
                for kind, vec, x in pts:
                    fid.add("loads.memberPoint", EXACT, "add_member_pt_load (local axes)")
                    for a in range(3):
                        if vec[a] != 0.0:
                            cmds.append(["add_member_pt_load", names[0], (LOCAL if kind == "force" else LOCAL_M)[a],
                                         float(vec[a]), float(x)])
                for w1, w2, x1, x2 in ds:
                    fid.add("loads.memberDistributed", EXACT, "add_member_dist_load (local axes, linearly varying)")
                    for a in range(3):
                        if w1[a] != 0.0 or w2[a] != 0.0:
                            cmds.append(["add_member_dist_load", names[0], LOCAL[a], float(w1[a]), float(w2[a]),
                                         float(x1), float(x2)])
            fefs.append(fef)
        members_out[mid] = {"first": groups[0], "last": groups[-1], "fef_i": fefs[0][:6].tolist(),
                            "fef_j": fefs[-1][6:].tolist()}

    # ------------------------------------------------------------------ nodes, supports, nodal loads
    head = [["add_node", name, *map(float, pos)] for name, pos in xyz.items()]
    for k, nid in enumerate(p.node_ids):
        flags = p.fixed[k].copy()
        flags[p.inactive.get(nid, [])] = True
        if flags.any():
            cmds.append(["def_support", f"N{k}", *[bool(v) for v in flags]])
        for d in range(6):
            if p.spring[k, d] != 0.0:
                fid.add("supports.spring", EXACT, "def_support_spring")
                cmds.append(["def_support_spring", f"N{k}", DOFS[d], float(p.spring[k, d])])
            if p.fixed[k, d] and p.prescribed[k, d] != 0.0:
                fid.add("loads.supportDisplacement", EXACT, "def_node_disp")
                cmds.append(["def_node_disp", f"N{k}", DOFS[d], float(p.prescribed[k, d])])
        if p.load[k].any():
            fid.add("loads.nodal", EXACT)
            loads[f"N{k}"] += p.load[k]
    for name, v in loads.items():
        for d in range(6):
            if v[d] != 0.0:
                cmds.append(["add_node_load", name, FORCES[d], float(v[d])])

    # PyNite connects every node that lies on a member's axis to that member
    tol = 1e-9 * size
    for name, ni, nj in segments:
        a, b = xyz[ni], xyz[nj]
        ab = b - a
        for other, pos in xyz.items():
            if other in (ni, nj):
                continue
            t = float((pos - a) @ ab) / float(ab @ ab)
            if 0.0 < t < 1.0 and float(np.linalg.norm(pos - (a + t * ab))) <= tol:
                fid.add("topology", UNSUPPORTED, f"node {other} lies on member {name} without being one of its "
                        "ends: PyNite would connect it to the member (benchfea does not)")

    job = {"solver": SOLVER, "model": head + cmds, "orient": orient, "secondOrder": p.second_order,
           "nodes": list(xyz), "members": all_names}
    options = {"subdivisions": n_sub, "stiffness": stiff}

    def mapper(raw):
        return _map(p, members_out, fid, raw, options)

    return Translation(SOLVER, combination, job, fid, controls, options, mapper)


def _map(p: Prepared, members_out, fid, raw, options) -> ForeignResult:
    n = len(p.node_ids)
    U = np.array([raw["disp"][f"N{k}"] for k in range(n)])
    Rx = np.array([raw["rxn"][f"N{k}"] for k in range(n)])
    mask = p.fixed | (p.spring != 0.0)
    for nid, dofs in p.inactive.items():
        mask[p.index(nid), dofs] = False
    reactions = np.where(mask, Rx, 0.0)
    end_forces = {}
    for mid, out in members_out.items():
        fi = sum(np.array(raw["f"][name])[:6] for name in out["first"]) + np.array(out["fef_i"])
        fj = sum(np.array(raw["f"][name])[6:] for name in out["last"]) + np.array(out["fef_j"])
        end_forces[mid] = np.r_[fi, fj] / p.alpha
    info = {"alpha": p.alpha, "version": raw.get("version"), "axesError": raw.get("axesError"), **options}
    if p.second_order:
        info["secondOrder"] = dict(raw.get("secondOrderInfo") or {})
        info["secondOrder"]["strain"] = max(abs(float(end_forces[mid][6])) * p.alpha / (m.E * m.A)
                                            for mid, m in p.members.items())
    return ForeignResult(SOLVER, p.combination, True, fid, U / p.alpha, reactions / p.alpha, end_forces, "", info)


def finalise(res: ForeignResult) -> None:
    """Second-order runs: add the two formulation-error estimates that come from the run itself.

    * PyNite builds Kg from the FIRST-ORDER axial forces.  With delta = largest change of a member
      axial force between the first-order run and the P-Delta run (relative to the largest axial
      force), kappa = relative size of the second-order increment of the displacements and the
      amplification 1/(1 - kappa), a first-order sensitivity estimate of the effect is
      ``delta * kappa / (1 - kappa)`` (heuristic: second-order increments scale with N).
    * PyNite's Kg has an axial term that changes the axial stiffness from EA/L to (EA + N)/L: relative
      effect at most ``max |N| / EA``.
    """
    so = res.info.get("secondOrder")
    if not so:
        return
    kappa = min(so["kappa"], 0.9)
    res.fidelity.add("secondOrder.axialForceIteration", APPROXIMATE,
                     "Kg from first-order axial forces, not iterated (estimate: delta * kappa / (1 - kappa))",
                     error=so["delta"] * kappa / (1.0 - kappa))
    res.fidelity.add("secondOrder.axialTermInKg", APPROXIMATE,
                     "axial stiffness (EA + N)/L in PyNite's geometric stiffness (estimate: max |N|/EA)",
                     error=so["strain"])


# ====================================================================================== worker
def execute(job: dict) -> dict:
    """Run a job with PyNite (normally in a subprocess)."""
    import Pynite
    from Pynite import FEModel3D

    model = FEModel3D()
    for c in job["model"]:
        getattr(model, c[0])(*c[1:])
    target = {}
    for name, y in job["orient"]:
        mem = model.members[name]
        y = np.array(y, float)
        T0 = np.array(mem.T())[:3, :3]
        angle = math.degrees(math.atan2(float(y @ T0[2]), float(y @ T0[1])))
        for trial in (angle, -angle):
            mem.rotation = trial
            if np.abs(np.array(mem.T())[1, :3] - y).max() < 1e-9:
                break
        else:
            return {"ok": False, "error": f"could not reproduce the local axes of member {name} in PyNite"}
        target[name] = y
    combo = "Combo 1"

    def collect():
        disp = {nm: [float(getattr(model.nodes[nm], d)[combo]) for d in DOFS] for nm in job["nodes"]}
        rxn = {nm: [float(getattr(model.nodes[nm], "Rxn" + d)[combo]) for d in FORCES] for nm in job["nodes"]}
        f, axes, axial = {}, 0.0, {}
        for nm in job["members"]:
            mem = model.members[nm]
            subs = list(mem.sub_members.values()) if getattr(mem, "sub_members", None) else [mem]
            # order the sub-members from end i to end j of the physical member
            xi = np.array([mem.i_node.X, mem.i_node.Y, mem.i_node.Z], float)
            subs.sort(key=lambda s: float(np.linalg.norm(np.array([s.i_node.X, s.i_node.Y, s.i_node.Z]) - xi)))
            fi = np.array(subs[0].f(combo), float).ravel()
            fj = np.array(subs[-1].f(combo), float).ravel()
            f[nm] = [*fi[:6], *fj[6:]]
            for s in subs:
                axes = max(axes, float(np.abs(np.array(s.T())[1, :3] - target[nm]).max()))
                axial[s.name] = float(np.array(s.f(combo)).ravel()[6])
        return disp, rxn, f, axes, axial

    out = {"ok": True, "version": getattr(Pynite, "__version__", None)}
    try:
        if job["secondOrder"]:
            model.analyze_linear(check_stability=False)
            d1, _, _, _, n1 = collect()
            model.analyze_PDelta(check_stability=False)
            disp, rxn, f, axes, n2 = collect()
            u1 = np.array([v[:3] for v in d1.values()])
            u2 = np.array([v[:3] for v in disp.values()])
            nmax = max((abs(v) for v in n2.values()), default=0.0) or 1.0
            out["secondOrderInfo"] = {
                "delta": max((abs(n2[k] - n1.get(k, 0.0)) for k in n2), default=0.0) / nmax,
                "kappa": float(np.abs(u2 - u1).max() / (np.abs(u2).max() or 1.0))}
        else:
            model.analyze_linear(check_stability=False)
            disp, rxn, f, axes, _ = collect()
    except Exception as e:                                   # noqa: BLE001 - library failure is a result
        return {"ok": False, "error": f"PyNite analysis failed: {type(e).__name__}: {e}"}
    if axes > 1e-9:
        return {"ok": False, "error": f"PyNite member axes differ from benchfea's by {axes:.2e}"}
    if not all(math.isfinite(v) for row in disp.values() for v in row):
        return {"ok": False, "error": "PyNite returned non-finite displacements (singular stiffness matrix)"}
    out.update({"disp": disp, "rxn": rxn, "f": f, "axesError": axes})
    return out
