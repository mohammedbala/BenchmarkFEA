"""Diagnosis of a singular (or numerically singular) stiffness matrix.

Input: the eigen-decomposition of the free-free stiffness scaled to unit diagonal.  Eigenvectors with
eigenvalue below the tolerance span the mechanisms.  A basis of a null space is not unique, so the
report is built only from basis-independent quantities:

1. rigid-body mechanisms: for every connected part of the structure the rigid-body motions that
   violate no support, spring or constraint (null space of a small matrix, reduced to row echelon
   form over TX, TY, TZ, RX, RY, RZ so that the basis is canonical);
2. single-node mechanisms: null vectors that live on one node only (row echelon form over its DOFs);
3. whatever remains ("distributed"): the participation of every node, i.e. the diagonal of the
   orthogonal projector onto the remaining null space (basis independent).

Rotations are multiplied by the model size to compare them with translations.
"""
from __future__ import annotations

import numpy as np

from ..model import DOFS
from .errors import MechanismDiagnosis, MechanismMode, NodeParticipation

AXES = ("X", "Y", "Z")
TOL = 1.0e-7            # relative tolerance of the geometric tests on unit-normalised null vectors
MAX_NODES = 12          # nodes listed per mode
TRUE_MECHANISM = 1.0e-12   # scaled eigenvalue at or below this: zero stiffness to round-off = a mechanism


def _rref(M: np.ndarray, tol: float = TOL) -> np.ndarray:
    """Reduced row echelon form of the row space of M (rows = vectors); zero rows dropped."""
    M = np.array(M, float)
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        if r == rows:
            break
        scale = max(np.abs(M[r:]).max(), 1e-300)
        p = r + int(np.argmax(np.abs(M[r:, c])))
        if abs(M[p, c]) <= tol * scale:
            M[r:, c] = 0.0
            continue
        M[[r, p]] = M[[p, r]]
        M[r] /= M[r, c]
        for k in range(rows):
            if k != r:
                M[k] -= M[k, c] * M[r]
        r += 1
    M = M[:r]
    M[np.abs(M) <= tol] = 0.0
    return M


def _null(M: np.ndarray, ncols: int, tol: float = TOL) -> np.ndarray:
    """Orthonormal basis (columns) of the null space of M (rows x ncols)."""
    if M.size == 0:
        return np.eye(ncols)
    _, s, vt = np.linalg.svd(np.atleast_2d(M), full_matrices=True)
    rank = int(np.sum(s > tol * max(s.max(), 1e-300))) if s.size else 0
    if s.size and s.max() <= 1e-300:
        rank = 0
    return vt[rank:].T


def _unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def _axis_name(v) -> str | None:
    v = _unit(v)
    for k, a in enumerate(AXES):
        if abs(abs(v[k]) - 1.0) <= 1e-6:
            return a
    return None


def _fmt(v, nd=4) -> str:
    return "(" + ", ".join(f"{x:.{nd}g}" for x in (float(round(c, 9)) + 0.0 for c in v)) + ")"


def _components(structure) -> list[list[int]]:
    """Connected parts (node indices) through members and constraints; nodes without members are left out
    unless a constraint ties them in."""
    model = structure.model
    parent = list(range(len(structure.node_ids)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    used = set()
    for e in structure.elements.values():
        a, b = (model.nodes[n].index for n in e.nodes)
        parent[find(a)] = find(b)
        used.update((a, b))
    for c in model.constraints.values():
        a = model.nodes[c.master].index
        for s in c.slaves:
            b = model.nodes[s].index
            parent[find(a)] = find(b)
            used.update((a, b))
    groups: dict[int, list[int]] = {}
    for a in sorted(used):
        groups.setdefault(find(a), []).append(a)
    return sorted(groups.values(), key=lambda g: g[0])


def _participation(structure, Wvecs: np.ndarray, members=True):
    """Node participations from weighted vectors (n, m) spanning a subspace with orthonormal columns."""
    p = (Wvecs ** 2).sum(axis=1).reshape(-1, 6)
    total = p.sum()
    if total <= 0:
        return (), ()
    per_node = p.sum(axis=1) / total
    order = sorted((k for k in range(len(per_node)) if per_node[k] > 1e-9), key=lambda k: (-round(per_node[k], 9), k))
    nodes = tuple(NodeParticipation(structure.node_ids[k], float(per_node[k]),
                                    {DOFS[c]: float(p[k, c] / p[k].sum()) for c in range(6)})
                  for k in order[:MAX_NODES])
    involved = {structure.node_ids[k] for k in order}
    mem = tuple(mid for mid, e in structure.elements.items() if e.nodes[0] in involved or e.nodes[1] in involved) \
        if members else ()
    return nodes, mem


def diagnose(structure, lam: np.ndarray, vec: np.ndarray, d: np.ndarray, tol: float) -> MechanismDiagnosis:
    n = structure.n
    null = lam < tol
    m = int(null.sum())
    numerics = {"smallestEigenvalue": float(lam[0]), "largestEigenvalue": float(lam[-1]),
                "eigenvaluesBelowTolerance": [float(x) for x in lam[null][:20]]}
    Lc = structure.size if structure.size > 0 else 1.0
    w = np.tile([1.0, 1.0, 1.0, Lc, Lc, Lc], n // 6)
    if m == 0:      # not positive definite to tolerance by the pivot test only
        null = np.zeros(len(lam), bool)
        null[0] = True
        m = 1
    ur = np.zeros((len(structure.retained), m))
    ur[structure.free_r] = vec[:, null] / d[:, None]
    N = structure.T @ ur                                    # (n, m) nodal null vectors
    Q, _ = np.linalg.qr(w[:, None] * N)                     # orthonormal in the weighted metric
    modes: list[MechanismMode] = []
    found = []                                              # weighted vectors of the classified modes

    def in_null(v_w):
        v_w = v_w / max(np.linalg.norm(v_w), 1e-300)
        return np.linalg.norm(v_w - Q @ (Q.T @ v_w)) <= 1e-6

    # ---- 1. rigid-body mechanisms of each connected part
    comps = _components(structure)
    whole = len(comps) == 1
    for comp in comps:
        ref = next((k for k in comp if structure.fixed[6 * k:6 * k + 3].any()), comp[0])
        x0 = structure.xyz[ref]
        B = np.zeros((n, 6))
        for k in comp:
            rel = structure.xyz[k] - x0
            for a in range(3):
                e = np.eye(3)[a]
                B[6 * k + a, a] = 1.0                       # translation along a
                B[6 * k:6 * k + 3, 3 + a] = np.cross(e, rel)    # rotation about a through x0
                B[6 * k + 3 + a, 3 + a] = 1.0
        dofs = np.concatenate([6 * k + np.arange(6) for k in comp])
        live = np.zeros(n, bool)
        live[dofs] = True
        live &= ~structure.inactive
        cs = np.array([1.0, 1.0, 1.0, 1.0 / Lc, 1.0 / Lc, 1.0 / Lc])      # coefficient scaling (rotations * Lc)
        Bw = (w[:, None] * B) * cs
        rows = [Bw[live & (structure.fixed | (structure.spring > 0.0))]]
        viol = Bw - (w[:, None] * (structure.T @ B[structure.retained])) * cs    # constraint violation
        rows.append(viol[live & structure.slave])
        Z = _null(np.vstack(rows), 6)
        if Z.shape[1] == 0:
            continue
        coef = _rref((Z * cs[:, None]).T)                   # rows = (t, r) in physical units, canonical
        for c in coef:
            v = B @ c
            full = np.linalg.norm(w * v)
            v[structure.inactive] = 0.0
            vw = w * v
            if np.linalg.norm(vw) <= 1e-9 * full or not in_null(vw):
                continue                 # moves only DOFs that are not part of the solution, or is not a null vector
            t, r = c[:3], c[3:]
            part = "" if whole else f" of the part containing node {structure.node_ids[comp[0]]}"
            if not r.any():
                ax = _axis_name(t)
                desc = (f"rigid-body translation in {ax}{part}: no restraint" if ax else
                        f"rigid-body translation along {_fmt(_unit(t))}{part}: no restraint")
                mode = dict(motion="translation", direction=tuple(_unit(t)))
            else:
                axis = _unit(r)
                pitch = float(t @ axis) / np.linalg.norm(r)
                point = x0 + np.cross(r, t) / float(r @ r)
                hit = next((structure.node_ids[k] for k in comp          # a node of this part on the axis
                            if np.linalg.norm(np.cross(structure.xyz[k] - point, axis)) <= 1e-7 * Lc), None)
                ax = _axis_name(axis)
                where = f"through node {hit}" if hit else f"through {_fmt(point)} mm"
                line = f"an axis parallel to {ax}" if ax else f"the axis {_fmt(axis)}"
                if abs(pitch) <= 1e-7 * Lc:
                    desc = f"rigid-body rotation about {line} {where}{part}: no restraint"
                    mode = dict(motion="rotation", axis=tuple(axis), through=tuple(point))
                else:
                    desc = (f"combined rigid-body rotation about {line} {where} and translation along "
                            f"{_fmt(_unit(t))}{part}: no restraint")
                    mode = dict(motion="translation+rotation", axis=tuple(axis), through=tuple(point),
                                direction=tuple(_unit(t)))
            q = vw / np.linalg.norm(vw)
            nodes, _ = _participation(structure, q[:, None], members=False)
            modes.append(MechanismMode("rigid-body", desc, cause="no-restraint", nodes=nodes, **mode))
            found.append(q)

    # ---- 2. mechanisms confined to a single node
    for k, nid in enumerate(structure.node_ids):
        Qa = Q[6 * k:6 * k + 6]
        if not Qa.any():
            continue
        s, vecs = np.linalg.eigh(Qa.T @ Qa)
        sel = s > 1.0 - 1e-7
        if not sel.any():
            continue
        local = (Qa @ vecs[:, sel]).T / w[6 * k:6 * k + 6]          # rows: physical 6-vectors at the node
        local = local / np.abs(local * w[:6]).max(axis=1, keepdims=True)
        for g in _rref(local):
            v = np.zeros(n)
            v[6 * k:6 * k + 6] = g
            vw = w * v
            vw /= np.linalg.norm(vw)
            if found and np.linalg.matrix_rank(np.column_stack(found + [vw]), tol=1e-6) == len(found):
                continue                                             # already reported (e.g. a one-node part)
            t, r = g[:3], g[3:]
            attached = structure.attached[nid]
            if r.any() and not t.any():
                ax = _axis_name(r)
                what = f"rotation about {ax}" if ax else f"rotation about the axis {_fmt(_unit(r))}"
                kind = dict(motion="rotation", axis=tuple(_unit(r)))
            elif t.any() and not r.any():
                ax = _axis_name(t)
                what = f"translation in {ax}" if ax else f"translation along {_fmt(_unit(t))}"
                kind = dict(motion="translation", direction=tuple(_unit(t)))
            else:
                what = f"combined translation {_fmt(t)} and rotation {_fmt(r)}"
                kind = dict(motion="translation+rotation")
            if attached and all(e.is_truss for e in attached):
                cause, why = "only-truss-members", ("only truss-type members connect to the node" if not t.any() else
                                                    "the connected truss-type members do not restrain it")
            elif attached and not t.any() and _all_ends_released(structure, nid, _unit(r)):
                cause, why = "all-member-ends-released", "all connected member ends released"
            else:
                cause, why = "no-stiffness", "no member, support or constraint resists it"
            share = (w[:6] * g) ** 2
            modes.append(MechanismMode(
                "node", f"node {nid} {what} unrestrained: {why}", cause=cause,
                nodes=(NodeParticipation(nid, 1.0, {DOFS[c]: float(share[c] / share.sum()) for c in range(6)}),),
                members=tuple(e.id for e in attached), **kind))
            found.append(vw)

    # ---- 3. the rest: participation of the remaining null space
    rank = np.linalg.matrix_rank(np.column_stack(found), tol=1e-6) if found else 0
    if rank < m:
        if found:
            C = Q.T @ np.column_stack(found)                 # (m, p)
            Zc = _null(C.T, m, tol=1e-6)                     # coordinates orthogonal to everything found
        else:
            Zc = np.eye(m)
        rest = Q @ Zc
        count = rest.shape[1]
        nodes, mem = _participation(structure, rest)
        listed = ", ".join(f"{p.node} ({'/'.join(k for k, v in p.dofs.items() if v > 1e-6)})" for p in nodes)
        t_share = float(sum(p.share * sum(v for k, v in p.dofs.items() if k.startswith("u")) for p in nodes))
        motion = "translation" if t_share > 1 - 1e-9 else ("rotation" if t_share < 1e-9 else "mixed")
        modes.append(MechanismMode(
            "distributed", f"{count} mechanism(s) involving several nodes; largest participation: {listed}",
            motion=motion, cause="no-stiffness", nodes=nodes, members=mem, count=count))
    kind = "mechanism" if lam[0] <= TRUE_MECHANISM else "ill-conditioned"
    return MechanismDiagnosis(m, tuple(modes), tol, len(structure.free_r), numerics, kind)


def _all_ends_released(structure, nid: str, axis: np.ndarray) -> bool:
    """True when every member end at the node has a full release (k = 0) for each local rotation that
    has a component along ``axis`` (frame members only; truss-type members never hold a rotation)."""
    for e in structure.attached[nid]:
        if e.is_truss:
            continue
        m = structure.model.members[e.id]
        comp = np.asarray(m.R, float) @ axis
        for end, node in (("i", e.nodes[0]), ("j", e.nodes[1])):
            if node != nid:
                continue
            k = e.k_end_i if end == "i" else e.k_end_j
            if any(abs(comp[a]) > 1e-9 and k[3 + a] != 0.0 for a in range(3)):
                return False
    return True
