"""Internal subdivision of members for the second-order and buckling analyses.

``refine(model, divisions)`` returns a :class:`Mesh`: a derived model in which every member is cut
into pieces at the load break points of all load cases (point loads, ends of line loads) and each
part between two break points is cut into ``divisions[member]`` equal pieces.  The original nodes
keep their index (and DOF numbers); the internal nodes are appended.  End releases and rigid end
offsets stay at the original member ends, internal joints are rigid.  Members listed in ``unsplit``
(slack tension/compression-only members) stay whole.

Loads are split with the pieces (:func:`split_loads`): a point load goes to the end of the piece that
finishes at it (to the start of the first piece when it sits at x = 0), a line load is cut at the piece
boundaries.  Results are mapped back to the stations of the original member (:func:`piece_stations`).
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np

from ..model import RIGID, Member, Node
from .beam import SpanLoads
from .loads import LoadSet, member_load_positions

__all__ = ["Mesh", "refine", "split_loads", "piece_stations"]

RIGID6 = (RIGID,) * 6
SEP = "@"          # not allowed in ids of the input format, so derived ids cannot collide


@dataclass
class Mesh:
    model: object                         # the refined model (same supports, constraints, analysis options)
    base: object                          # the original model
    pieces: dict                          # original member id -> list of piece ids (start to end)
    bounds: dict                          # original member id -> np.ndarray of piece boundaries (x from the start)
    n_nodes: int                          # number of original nodes

    @property
    def n_dof_base(self) -> int:
        return 6 * self.n_nodes


def refine(model, divisions: dict | int = 1, unsplit=()) -> Mesh:
    pos = member_load_positions(model)
    nodes = dict(model.nodes)
    members, pieces, bounds = {}, {}, {}
    for mid, m in model.members.items():
        L = float(m.length)
        if mid in unsplit:
            members[mid], pieces[mid], bounds[mid] = m, [mid], np.array([0.0, L])
            continue
        n = int(divisions if isinstance(divisions, int) else divisions.get(mid, 1))
        brk = sorted({0.0, L} | {x for group in pos[mid] for x in group if 0.0 < x < L})
        b = [0.0]
        for a, c in zip(brk[:-1], brk[1:]):
            b += [a + (c - a) * k / n for k in range(1, n)] + [c]
        b = np.array(b)
        ex = np.asarray(m.R, float)[0]
        ids = [f"{mid}{SEP}{k}" for k in range(len(b) - 1)]
        nid = [m.i] + [f"{mid}{SEP}n{k}" for k in range(1, len(b) - 1)] + [m.j]
        for k in range(1, len(b) - 1):
            nodes[nid[k]] = Node(nid[k], len(nodes), m.start + b[k] * ex)
        zero = np.zeros(3)
        for k, pid in enumerate(ids):
            first, last = k == 0, k == len(ids) - 1
            members[pid] = Member(
                pid, m.type, nid[k], nid[k + 1], m.section, m.material,
                m.start + b[k] * ex, m.start + b[k + 1] * ex, float(b[k + 1] - b[k]), m.R,
                m.offset_i if first else zero, m.offset_j if last else zero,
                tuple(m.k_end_i) if first else RIGID6, tuple(m.k_end_j) if last else RIGID6)
        pieces[mid], bounds[mid] = ids, b
    refined = dataclasses.replace(model, nodes=nodes, members=members, design_members={}, defaults_used=[],
                                  warnings=[])
    return Mesh(refined, model, pieces, bounds, len(model.nodes))


def split_loads(mesh: Mesh, loads: LoadSet) -> LoadSet:
    """The load set of the original model expressed on the refined model."""
    nodal = np.zeros(mesh.model.n_dof)
    nodal[:mesh.n_dof_base] = loads.nodal
    out = LoadSet(nodal, dict(loads.prescribed))
    for mid, span in loads.member.items():
        ids, b = mesh.pieces[mid], mesh.bounds[mid]
        if len(ids) == 1 and ids[0] == mid:
            out.member[mid] = span
            continue
        L = float(b[-1])
        pts = {pid: [] for pid in ids}
        dst = {pid: [] for pid in ids}
        for x, f in span.points:
            x = min(max(float(x), 0.0), L)
            k = int(np.argmin(np.abs(b - x)))                # load positions are piece boundaries
            if k == 0:
                pts[ids[0]].append((0.0, f))
            else:
                pts[ids[k - 1]].append((float(b[k] - b[k - 1]), f))
        for x1, x2, w1, w2 in span.dists:
            if not x2 > x1:
                continue
            w1, w2 = np.asarray(w1, float), np.asarray(w2, float)
            for k, pid in enumerate(ids):
                a, c = b[k], b[k + 1]
                if x1 <= a and c <= x2:
                    dst[pid].append((0.0, float(c - a), w1 + (w2 - w1) * ((a - x1) / (x2 - x1)),
                                     w1 + (w2 - w1) * ((c - x1) / (x2 - x1))))
        for pid in ids:
            if pts[pid] or dst[pid]:
                out.member[pid] = SpanLoads(tuple(pts[pid]), tuple(dst[pid]))
    return out


def piece_stations(bounds: np.ndarray, stations) -> list[list]:
    """For the stations (x, side) of an original member: per piece, the list of
    (index in ``stations``, local x, side)."""
    m = len(bounds) - 1
    out = [[] for _ in range(m)]
    for n, (x, side) in enumerate(stations):
        hit = np.flatnonzero(bounds == x)
        if len(hit):
            k = int(hit[0])
            if k == 0:
                out[0].append((n, 0.0, side))
            elif k == m:
                out[m - 1].append((n, float(bounds[m] - bounds[m - 1]), side))
            elif side < 0:
                out[k - 1].append((n, float(bounds[k] - bounds[k - 1]), -1))
            else:
                out[k].append((n, 0.0, +1))
        else:
            k = min(max(int(np.searchsorted(bounds, x)) - 1, 0), m - 1)
            out[k].append((n, float(x - bounds[k]), side))
    return out
