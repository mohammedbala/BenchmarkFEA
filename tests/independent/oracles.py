"""Direct translations of a raw benchfea input dictionary to OpenSeesPy and PyNite (independent review).

These translators were written for the review and do not use ``benchfea.crosscheck`` or any other
benchfea code: geometry, local axes and local load components come from ``tableau.Tableau`` (the
reviewer's own reading of docs/CONVENTIONS.md), everything is in the units of the file.

Each ``run_*`` returns ``None`` when the solver cannot represent the model EXACTLY (the caller then
skips that oracle for that model) or a dict ``{"displacements" (n, 6), "reactions" (n, 6),
"end_forces" {member: (12,)}}`` in benchfea conventions.  Facts about the two packages used here
(element argument order, ``vecxz``, ``-jntOffset``, ``-releasez/-releasey``, what PDelta contains, what
PyNite's P-Delta does) were measured in this project (docs/CROSSCHECK.md).
"""
from __future__ import annotations

import math

import numpy as np

from tableau import Tableau

TRUSS = ("truss", "tensionOnly", "compressionOnly")


def _axis_aligned(R):
    return all(abs(abs(row).max() - 1.0) < 1e-12 for row in R)


def _truss_only_nodes(t):
    att = {n: [] for n in t.node_ids}
    for m in t.members.values():
        att[m.ni].append(m)
        att[m.nj].append(m)
    return [n for n, ms in att.items() if ms and all(m.type in TRUSS and not (m.oi.any() or m.oj.any()) for m in ms)]


def have(name):
    try:
        __import__({"opensees": "openseespy.opensees", "pynite": "Pynite"}[name])
        return True
    except Exception:
        return False


# ============================================================================== OpenSees
def run_opensees(raw, cid, second_order=False, nsub=1):
    import openseespy.opensees as ops
    t = Tableau(raw)
    t.set_combination(cid)
    if any(m.type in ("tensionOnly", "compressionOnly") for m in t.members.values()):
        return None
    slaves = {s for c in raw.get("constraints", {}).values() for s in c["slaves"]}
    if any(s in raw.get("supports", {}) for s in slaves):
        return None
    for m in t.members.values():
        if any(a < m.L and f[3:].any() for a, f in m.points) or any(f[3:].any() for a, f in m.points):
            return None                                   # member point moments: no element load for them
        if any(not np.allclose(w1, w2, rtol=0, atol=0) for _, _, w1, w2 in m.dists):
            return None                                   # trapezoids are silently wrong in OpenSees 3.8 (3D)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    tag = {n: k + 1 for k, n in enumerate(t.node_ids)}
    for n in t.node_ids:
        ops.node(tag[n], *map(float, t.xyz[n]))
    nxt = [len(t.node_ids) + 1]

    def new_node(xyz):
        ops.node(nxt[0], *map(float, xyz))
        nxt[0] += 1
        return nxt[0] - 1
    fixed = t.fixed.reshape(-1, 6).copy()
    for n in _truss_only_nodes(t):
        k = t.idx[n]
        for a in (3, 4, 5):
            if not t.spring[6 * k + a] and n not in slaves:
                fixed[k, a] = True
    etag, mtag = [1], [1]
    spring_nodes = {}
    held_nodes = []
    for k, n in enumerate(t.node_ids):
        if fixed[k].any():
            ops.fix(tag[n], *[int(v) for v in fixed[k]])
        ks = t.spring[6 * k:6 * k + 6]
        if ks.any():
            g = new_node(t.xyz[n])
            ops.fix(g, 1, 1, 1, 1, 1, 1)
            mats, dirs = [], []
            for a in range(6):
                if ks[a]:
                    ops.uniaxialMaterial("Elastic", mtag[0], float(ks[a]))
                    mats.append(mtag[0])
                    dirs.append(a + 1)
                    mtag[0] += 1
            ops.element("zeroLength", 100000 + etag[0], g, tag[n], "-mat", *mats, "-dir", *dirs)
            etag[0] += 1
            spring_nodes[n] = g
    for c in raw.get("constraints", {}).values():
        if c["type"] == "rigidLink":
            for s in c["slaves"]:
                ops.rigidLink("beam", tag[c["master"]], tag[s])
        else:
            perp = {"X": 1, "Y": 2, "Z": 3}[c.get("normal", t.vertical)]
            ops.rigidDiaphragm(perp, tag[c["master"]], *[tag[s] for s in c["slaves"]])
    ttag = [1]
    first_last = {}
    sub_elems = {}
    for m in t.members.values():
        d = raw["members"][m.id]
        mat, sec = raw["materials"][d["material"]], raw["sections"][d["section"]]
        E = mat["E"]
        G = mat["G"] if "G" in mat else E / (2 * (1 + mat["nu"]))
        A, J, Iy, Iz = sec["A"], sec.get("J", 1.0), sec.get("Iy", 1.0), sec.get("Iz", 1.0)
        shear = bool(m.by or m.bz)
        ends = [tag[m.ni], tag[m.nj]]
        flags = {4: 0, 5: 0}                             # local ry, rz -> release code
        if m.type in TRUSS:
            flags = {4: 3, 5: 3}
            J = J * 1e-12
            if shear:
                return None
        else:
            for e in (0, 1):
                ke = m.k[6 * e:6 * e + 6]
                if np.all(ke == np.inf):
                    continue
                off = m.oi if e == 0 else m.oj
                # a duplicate node tied to a constraint slave would be a chained constraint, which the
                # Transformation handler does not resolve (measured: the element loads are lost)
                if off.any() or not _axis_aligned(m.R) or (m.ni if e == 0 else m.nj) in slaves:
                    simple = all(ke[a] == np.inf for a in (0, 1, 2, 3)) and all(ke[a] in (0.0, np.inf) for a in (4, 5))
                    if not simple or shear:
                        return None
                    for a in (4, 5):
                        if ke[a] == 0.0:
                            flags[a] |= (1 if e == 0 else 2)
                    continue
                # member axes parallel to the global axes, no offset: duplicate node, exact for every load
                node = m.ni if e == 0 else m.nj
                dup = new_node(t.xyz[node])
                mats, dirs, same, held = [], [], [], []
                for a in range(6):
                    gdir = int(np.abs(m.R[a % 3]).argmax()) + (3 if a >= 3 else 0)
                    if ke[a] == np.inf:
                        g = 6 * t.idx[node] + gdir
                        if t.fixed[g] and t.pres[g]:
                            return None
                        (held if t.fixed[g] else same).append(gdir + 1)   # a fixed retained DOF: fix the copy
                    elif ke[a] > 0:
                        ops.uniaxialMaterial("Elastic", mtag[0], float(ke[a]))
                        mats.append(mtag[0])
                        dirs.append(gdir + 1)
                        mtag[0] += 1
                if held:
                    ops.fix(dup, *[int(a + 1 in held) for a in range(6)])
                    held_nodes.append((node, dup, [a - 1 for a in held]))
                if same:
                    ops.equalDOF(tag[node], dup, *sorted(same))      # one constraint per node pair
                if mats:
                    ops.element("zeroLength", 100000 + etag[0], tag[node], dup, "-mat", *mats, "-dir", *dirs)
                    etag[0] += 1
                ends[e] = dup
        full = all(x1 == 0.0 and x2 == m.L for x1, x2, _, _ in m.dists) and not m.points
        if (flags[4] or flags[5]) and not full:
            return None         # measured: the release flags condense the fixed-end forces of full-span uniform loads only
        if shear and (not full or m.oi.any() or m.oj.any()):
            return None         # measured: ElasticTimoshenkoBeam takes no -beamPoint / partial load, mishandles -jntOffset
        xs = [m.L * k / nsub for k in range(nsub + 1)]
        chain = [ends[0]] + [new_node(m.start + m.R[0] * x) for x in xs[1:-1]] + [ends[1]]
        tags = []
        for k in range(nsub):
            oi = m.oi if k == 0 else np.zeros(3)
            oj = m.oj if k == nsub - 1 else np.zeros(3)
            ops.geomTransf("PDelta" if second_order else "Linear", ttag[0], *map(float, m.R[2]),
                           "-jntOffset", *map(float, oi), *map(float, oj))
            if shear:
                Avy = 1.0 / (G * m.by) if m.by else 1e9 * A
                Avz = 1.0 / (G * m.bz) if m.bz else 1e9 * A
                ops.element("ElasticTimoshenkoBeam", etag[0], chain[k], chain[k + 1], E, G, A, J, Iy, Iz, Avy, Avz, ttag[0])
            else:
                args = []
                for a, name in ((5, "-releasez"), (4, "-releasey")):
                    code = (flags[a] & 1 if k == 0 else 0) | (flags[a] & 2 if k == nsub - 1 else 0)
                    if code:
                        args += [name, code]
                ops.element("elasticBeamColumn", etag[0], chain[k], chain[k + 1], A, E, G, J, Iy, Iz, ttag[0], *args)
            tags.append(etag[0])
            etag[0] += 1
            ttag[0] += 1
        sub_elems[m.id] = (tags, xs)
        first_last[m.id] = (tags[0], tags[-1])
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for k, n in enumerate(t.node_ids):
        f = t.F[6 * k:6 * k + 6]
        if f.any():
            ops.load(tag[n], *map(float, f))
        for a in range(6):
            if t.fixed[6 * k + a] and t.pres[6 * k + a]:
                ops.sp(tag[n], a + 1, float(t.pres[6 * k + a]))
    for m in t.members.values():
        tags, xs = sub_elems[m.id]
        for x1, x2, w1, _ in m.dists:
            for k, e in enumerate(tags):
                a, b = max(x1, xs[k]), min(x2, xs[k + 1])
                if b - a > 1e-12 * m.L:
                    h = xs[k + 1] - xs[k]
                    part = [] if (a == xs[k] and b == xs[k + 1]) else [float((a - xs[k]) / h), float((b - xs[k]) / h)]
                    ops.eleLoad("-ele", e, "-type", "-beamUniform", float(w1[1]), float(w1[2]), float(w1[0]), *part)
        for a, f in m.points:
            k = min(int(np.searchsorted(xs, a, side="right")) - 1, len(tags) - 1)
            h = xs[k + 1] - xs[k]
            ops.eleLoad("-ele", tags[k], "-type", "-beamPoint", float(f[1]), float(f[2]), float((a - xs[k]) / h), float(f[0]))
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("FixedNumIter", 12 if second_order else 3)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("OpenSees analysis failed")
    ops.reactions()
    n = len(t.node_ids)
    U = np.array([ops.nodeDisp(tag[nid]) for nid in t.node_ids])
    R = np.array([ops.nodeReaction(tag[nid]) for nid in t.node_ids]) * t.fixed.reshape(-1, 6)
    for nid, g in spring_nodes.items():
        R[t.idx[nid]] += np.array(ops.nodeReaction(g)) * (t.spring[6 * t.idx[nid]:6 * t.idx[nid] + 6] != 0)
    for nid, g, dofs in held_nodes:
        R[t.idx[nid], dofs] += np.array(ops.nodeReaction(g))[dofs]
    for c in raw.get("constraints", {}).values():
        R[t.idx[c["master"]]] = np.nan       # OpenSees does not include the constraint forces in the reaction of a retained node
    ef = {}
    for mid, (a, b) in first_last.items():
        fa, fb = ops.eleResponse(a, "localForce"), ops.eleResponse(b, "localForce")
        ef[mid] = np.array(list(fa[:6]) + list(fb[6:]))
    ops.wipe()
    return {"displacements": U, "reactions": R, "end_forces": ef}


# ============================================================================== PyNite
def run_pynite(raw, cid, second_order=False, nsub=1):
    from Pynite import FEModel3D
    t = Tableau(raw)
    t.set_combination(cid)
    if raw.get("constraints"):
        return None
    for m in t.members.values():
        if m.by or m.bz or m.oi.any() or m.oj.any():
            return None
        if m.type == "frame" and any(k not in (0.0, np.inf) for k in m.k):
            return None
    fe = FEModel3D()
    for n in t.node_ids:
        fe.add_node(n, *map(float, t.xyz[n]))
    fixed = t.fixed.reshape(-1, 6).copy()
    for n in _truss_only_nodes(t):
        fixed[t.idx[n], 3:] |= (t.spring[6 * t.idx[n] + 3:6 * t.idx[n] + 6] == 0)
    names = ("DX", "DY", "DZ", "RX", "RY", "RZ")
    for k, n in enumerate(t.node_ids):
        if fixed[k].any():
            fe.def_support(n, *[bool(v) for v in fixed[k]])
        for a in range(6):
            if t.spring[6 * k + a]:
                fe.def_support_spring(n, names[a], float(t.spring[6 * k + a]))
            if t.fixed[6 * k + a] and t.pres[6 * k + a]:
                fe.def_node_disp(n, names[a], float(t.pres[6 * k + a]))
    ends = {}
    for m in t.members.values():
        d = raw["members"][m.id]
        mat, sec = raw["materials"][d["material"]], raw["sections"][d["section"]]
        E = mat["E"]
        G = mat["G"] if "G" in mat else E / (2 * (1 + mat["nu"]))
        mname, sname = "M_" + m.id, "S_" + m.id
        fe.add_material(mname, E, G, 0.3, 0.0)
        fe.add_section(sname, sec["A"], sec.get("Iy", 1.0), sec.get("Iz", 1.0), sec.get("J", 1.0))
        chain = [m.ni]
        for k in range(1, nsub):
            nm = f"{m.id}__n{k}"
            fe.add_node(nm, *map(float, m.start + m.R[0] * (m.L * k / nsub)))
            chain.append(nm)
        chain.append(m.nj)
        subs = []
        for k in range(nsub):
            nm = m.id if nsub == 1 else f"{m.id}__{k}"
            fe.add_member(nm, chain[k], chain[k + 1], mname, sname, 0.0,
                          tension_only=m.type == "tensionOnly", comp_only=m.type == "compressionOnly")
            T = np.array(fe.members[nm].T())[:3, :3]             # PyNite's own default triad (rows x, y, z)
            ang = math.degrees(math.atan2(m.R[1] @ T[2], m.R[1] @ T[1]))
            fe.members[nm].rotation = ang
            T = np.array(fe.members[nm].T())[:3, :3]
            if not np.allclose(T, m.R, atol=1e-9):
                raise RuntimeError(f"PyNite triad of {m.id} does not reproduce the reference triad")
            subs.append(nm)
        rel_i = [bool(k == 0.0) for k in m.k[:6]]
        rel_j = [bool(k == 0.0) for k in m.k[6:]]
        if nsub == 1:
            if any(rel_i + rel_j):
                fe.def_releases(m.id, *rel_i, *rel_j)
        else:
            if any(rel_i):
                fe.def_releases(subs[0], *rel_i, *([False] * 6))
            if any(rel_j):
                fe.def_releases(subs[-1], *([False] * 6), *rel_j)
        xs = [m.L * k / nsub for k in range(nsub + 1)]
        for x1, x2, w1, w2 in m.dists:
            for k, nm in enumerate(subs):
                a, b = max(x1, xs[k]), min(x2, xs[k + 1])
                if b - a > 1e-12 * m.L:
                    wa = w1 + (w2 - w1) * (a - x1) / (x2 - x1)
                    wb = w1 + (w2 - w1) * (b - x1) / (x2 - x1)
                    for c, dn in enumerate(("Fx", "Fy", "Fz")):
                        if wa[c] or wb[c]:
                            fe.add_member_dist_load(nm, dn, float(wa[c]), float(wb[c]), float(a - xs[k]), float(b - xs[k]), case="C")
        for a, f in m.points:
            k = min(int(np.searchsorted(xs, a, side="right")) - 1, nsub - 1)
            for c, dn in enumerate(("Fx", "Fy", "Fz", "Mx", "My", "Mz")):
                if f[c]:
                    fe.add_member_pt_load(subs[k], dn, float(f[c]), float(a - xs[k]), case="C")
        ends[m.id] = (subs[0], subs[-1])
    for k, n in enumerate(t.node_ids):
        for a, dn in enumerate(("FX", "FY", "FZ", "MX", "MY", "MZ")):
            if t.F[6 * k + a]:
                fe.add_node_load(n, dn, float(t.F[6 * k + a]), case="C")
    fe.add_load_combo("K", {"C": 1.0})
    if second_order:
        fe.analyze_PDelta(check_stability=False)
    else:
        fe.analyze(check_stability=False, check_statics=False)
    U = np.array([[getattr(fe.nodes[n], a)["K"] for a in names] for n in t.node_ids])
    rn = ("RxnFX", "RxnFY", "RxnFZ", "RxnMX", "RxnMY", "RxnMZ")
    R = np.array([[getattr(fe.nodes[n], a)["K"] for a in rn] for n in t.node_ids])
    R = R * (t.fixed.reshape(-1, 6) | (t.spring.reshape(-1, 6) != 0))
    ef = {}
    for mid, (a, b) in ends.items():
        fa = np.asarray(fe.members[a].f("K")).ravel()
        fb = np.asarray(fe.members[b].f("K")).ravel()
        ef[mid] = np.concatenate([fa[:6], fb[6:]])
    slack = {mid for mid, (a, b) in ends.items() if not fe.members[a].active["K"]}
    return {"displacements": U, "reactions": R, "end_forces": ef, "slack": slack}
