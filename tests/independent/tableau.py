"""Independent reference solver for the review of ``benchfea.solve`` ("tableau" formulation).

Written for the independent review; it shares NO code with ``benchfea`` (not even the loader, the unit
conversion, the local-axis routine or the section routines): it reads the raw input dictionary, works
in the units of the file and implements docs/CONVENTIONS.md on its own.  It is deliberately a different
numerical method from the engine:

* **No stiffness matrix, no fixed-end forces, no flexibility matrix, no master/slave elimination.**
  The unknowns are the node displacements, for every member the 12 beam-side end displacements and the
  6 stress resultants at its start, one reaction per fixed DOF and one multiplier per constraint
  equation.  The equations are written down directly:

  1. member: the state z = (u, v, w, tx, ty, tz, N, Vy, Vz, T, My, Mz) obeys the first-order system
     z' = A z + g(x) (kinematics + elasticity + equilibrium of a slice, written out in ``_state_matrix``);
     its solution over the member is z(L) = Phi z(0) + p with Phi = exp(A L) (``scipy.linalg.expm`` of
     the full 12x12 system) and p by Gauss-Legendre quadrature of the variation-of-constants integral.
     The first six rows are 6 equations per member (compatibility).
  2. end connections, per end DOF: rigid -> beam side = node side (through the rigid offset);
     spring -> end force = k (node side - beam side); released -> end force = 0.
  3. node equilibrium: member end forces (carried through the offsets), nodal loads, support reactions,
     support springs and constraint multipliers sum to zero.
  4. supports: u = imposed value;  5. constraints: slave = rigid-body image of the master.

  The resulting unsymmetric system is equilibrated and solved by LU with partial pivoting.

* Second order (docs/THEORY.md 10.2 as the *statement of the problem*, not its code): the slice is in
  equilibrium in its displaced position under the axial force N (tension +) acting along the undeformed
  axis: Mz' = -Vy + N v', My' = Vz - N w', shear strain driven by V - N v' (Engesser form).  A changes
  with N; exp(A L) then contains the stability functions.  N is iterated.
* Buckling: the load factor at which the determinant of the tableau (with N = factor x first-order
  axial forces) changes sign (scan + Brent) -- the transcendental problem itself, no geometric
  stiffness matrix.
* Tension/compression-only members: plain active-set loop.

Sign conventions of the outputs are those of docs/CONVENTIONS.md (reactions on the structure, end
forces on the member in local axes, station forces on the +x face).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from scipy.optimize import brentq

DOF = ("ux", "uy", "uz", "rx", "ry", "rz")
AX = {"X": 0, "Y": 1, "Z": 2}
GP, GW = np.polynomial.legendre.leggauss(10)


def cross_matrix(v):
    x, y, z = v
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0.0]])


def triad(p, q, vertical, ref=None, roll=0.0):
    """Rows = local x, y, z (CONVENTIONS section 2)."""
    x = (q - p) / np.linalg.norm(q - p)
    up = np.eye(3)[AX[vertical]]
    if ref is None:
        v = up if np.linalg.norm(np.cross(x, up)) > 1e-6 else np.array([1.0, 0, 0])
    else:
        v = np.asarray(ref, float)
    z = np.cross(x, v)
    z /= np.linalg.norm(z)
    y = np.cross(z, x)
    r = math.radians(roll)
    return np.array([x, y * math.cos(r) + z * math.sin(r), -y * math.sin(r) + z * math.cos(r)])


class Member:
    def __init__(self, mid, d, data):
        self.id = mid
        self.type = d.get("type", "frame")
        self.ni, self.nj = d["i"], d["j"]
        Xi, Xj = data.xyz[self.ni], data.xyz[self.nj]
        ori = d.get("orientation", {})
        ref, roll = ori.get("referenceVector"), ori.get("roll", 0.0)
        off = d.get("offsets", {})
        oi, oj = np.array(off.get("i", [0, 0, 0]), float), np.array(off.get("j", [0, 0, 0]), float)
        if off.get("axes", "global") == "local":
            Rn = triad(Xi, Xj, data.vertical, ref, roll)
            oi, oj = Rn.T @ oi, Rn.T @ oj
        self.oi, self.oj = oi, oj
        self.start, self.end = Xi + oi, Xj + oj
        self.L = float(np.linalg.norm(self.end - self.start))
        self.R = triad(self.start, self.end, data.vertical, ref, roll)
        mat, sec = data.raw["materials"][d["material"]], data.raw["sections"][d["section"]]
        E = mat["E"]
        G = mat["G"] if "G" in mat else E / (2 * (1 + mat["nu"]))
        self.gammaA = mat.get("unitWeight", 0.0) * sec["A"]
        self.EA, self.GJ = E * sec["A"], G * sec.get("J", 1.0)
        self.EIy, self.EIz = E * sec.get("Iy", 1.0), E * sec.get("Iz", 1.0)
        shear = data.raw.get("analysis", {}).get("shearDeformation", True)
        self.by = 1.0 / (G * sec["Asy"]) if shear and sec.get("Asy") else 0.0
        self.bz = 1.0 / (G * sec["Asz"]) if shear and sec.get("Asz") else 0.0
        # end connections: inf rigid, 0 released, k spring
        k = np.full(12, np.inf)
        if self.type != "frame":
            k[[4, 5, 9, 10, 11]] = 0.0
        else:
            for e, base in (("i", 0), ("j", 6)):
                for name, v in d.get("releases", {}).get(e, {}).items():
                    a = base + DOF.index(name)
                    k[a] = 0.0 if v is True else (np.inf if v is False else float(v["spring"]))
        self.k = k
        self.active = True
        self.points, self.dists = [], []        # loads of the current combination, local axes
        self.Nstart = None                      # None = first order; else axial resultant at x = 0-

    # ------------------------------------------------------------------ member ODE
    def axial(self, x):
        """Axial resultant just left of x from statics, given N(0-) = Nstart."""
        n = self.Nstart
        for a, f in self.points:
            if a < x:
                n -= f[0]
        for x1, x2, w1, w2 in self.dists:
            b = min(x, x2)
            if b > x1:
                s = (w2[0] - w1[0]) / (x2 - x1)
                n -= w1[0] * (b - x1) + s * (b - x1) ** 2 / 2
        return n

    def _state_matrix(self, N):
        A = np.zeros((12, 12))
        dy, dz = 1 + self.by * N, 1 + self.bz * N
        A[0, 6] = 1 / self.EA                       # u' = N / EA
        A[1, 5], A[1, 7] = 1 / dy, self.by / dy     # v' = (tz + Vy/GAsy) / (1 + N/GAsy)
        A[2, 4], A[2, 8] = -1 / dz, self.bz / dz    # w' = (-ty + Vz/GAsz) / (1 + N/GAsz)
        A[3, 9] = 1 / self.GJ
        A[4, 10] = 1 / self.EIy
        A[5, 11] = 1 / self.EIz
        A[10, 8], A[10, 4] = 1 / dz, N / dz         # My' = Vz - N w'
        A[11, 7], A[11, 5] = -1 / dy, N / dy        # Mz' = -Vy + N v'
        return A

    def _phi(self, h, N):
        L, EI = self.L, max(self.EIy, self.EIz)
        s = np.array([1 / L] * 3 + [1.0] * 3 + [L * L / EI] * 3 + [L / EI] * 3)
        return (expm((self._state_matrix(N) * h) * s[:, None] / s[None, :])) / s[:, None] * s[None, :]

    def _w(self, x):
        w = np.zeros(3)
        for x1, x2, w1, w2 in self.dists:
            if x1 <= x <= x2:
                w += w1 + (w2 - w1) * (x - x1) / (x2 - x1)
        return w

    def _breaks(self, extra=(), nsub=1):
        xs = {0.0, self.L} | {float(a) for a, _ in self.points} | {float(v) for d in self.dists for v in d[:2]}
        xs |= {float(x) for x in extra} | {self.L * k / nsub for k in range(1, nsub)}
        return sorted(x for x in xs if 0.0 <= x <= self.L)

    def march(self, stations=(), nsub=1):
        """Affine map of the start state: returns [(x, side, Phi, p)] with z(x, side) = Phi z(0-) + p for the
        requested stations, and the same for (L, +1) as last entry."""
        want = {}
        for x, side in stations:
            want.setdefault(float(x), []).append(side)
        bs = self._breaks(want, nsub)
        Phi, p = np.eye(12), np.zeros(12)
        out = []

        def jump(x, p):
            p = p.copy()
            for a, f in self.points:
                if a == x:
                    p[6:] -= f
            return p
        for n, x in enumerate(bs):
            if -1 in want.get(x, []):
                out.append((x, -1, Phi.copy(), p.copy()))
            p = jump(x, p)
            if +1 in want.get(x, []):
                out.append((x, +1, Phi.copy(), p.copy()))
            if n + 1 < len(bs):
                h = bs[n + 1] - x
                N = 0.0 if self.Nstart is None else self.axial(x + h / 2)
                seg = self._phi(h, N)
                q = np.zeros(12)
                if self.dists:
                    for t, wt in zip(GP, GW):
                        tk = h / 2 * (1 + t)
                        w = self._w(x + tk)
                        if w.any():
                            g = np.zeros(12)
                            g[6:9] = -w
                            q += wt * h / 2 * (self._phi(h - tk, N) @ g)
                Phi, p = seg @ Phi, seg @ p + q
        out.append((self.L, +1, Phi, p))
        return out


class Tableau:
    def __init__(self, raw):
        self.raw = raw
        self.vertical = raw["conventions"]["verticalAxis"]
        self.node_ids = list(raw["nodes"])
        self.idx = {n: k for k, n in enumerate(self.node_ids)}
        self.xyz = {n: np.array([v["x"], v["y"], v.get("z", 0.0)], float) for n, v in raw["nodes"].items()}
        self.members = {mid: Member(mid, d, self) for mid, d in raw["members"].items()}
        nn = len(self.node_ids)
        self.fixed = np.zeros(6 * nn, bool)
        self.spring = np.zeros(6 * nn)
        for nid, s in raw.get("supports", {}).items():
            st = {"fixed": [1] * 6, "pinned": [1, 1, 1, 0, 0, 0], None: [0] * 6,
                  "roller": [int(a == AX[self.vertical]) for a in range(3)] + [0, 0, 0]}[s.get("type")]
            for a in range(6):
                self.fixed[6 * self.idx[nid] + a] = bool(st[a])
            for name, v in s.get("dofs", {}).items():
                g = 6 * self.idx[nid] + DOF.index(name)
                self.fixed[g] = v == "fixed"
                self.spring[g] = v["spring"] if isinstance(v, dict) else 0.0
        # constraint rows g . U = 0
        self.crows = []
        for c in raw.get("constraints", {}).values():
            m = self.idx[c["master"]]
            for sid in c["slaves"]:
                s = self.idx[sid]
                r = self.xyz[sid] - self.xyz[c["master"]]
                if c["type"] == "rigidLink":
                    for a in range(3):          # u_s - u_m - theta_m x r = 0 ;  theta x r = -[r]x theta
                        row = np.zeros(6 * nn)
                        row[6 * s + a], row[6 * m + a] = 1.0, -1.0
                        row[6 * m + 3:6 * m + 6] = cross_matrix(r)[a]
                        self.crows.append(row)
                    for a in range(3, 6):
                        row = np.zeros(6 * nn)
                        row[6 * s + a], row[6 * m + a] = 1.0, -1.0
                        self.crows.append(row)
                else:
                    n = AX[c.get("normal", self.vertical)]
                    e = np.eye(3)[n]
                    r = r - (r @ e) * e
                    lever = np.cross(e, r)       # displacement of the slave per unit rotation about n
                    for a in range(3):
                        if a == n:
                            continue
                        row = np.zeros(6 * nn)
                        row[6 * s + a], row[6 * m + a], row[6 * m + 3 + n] = 1.0, -1.0, -lever[a]
                        self.crows.append(row)
                    row = np.zeros(6 * nn)
                    row[6 * s + 3 + n], row[6 * m + 3 + n] = 1.0, -1.0
                    self.crows.append(row)
        self.second_order = False
        self.arm = True
        self.nsub = 1

    # ------------------------------------------------------------------ loads
    def set_combination(self, cid):
        nn = len(self.node_ids)
        self.F = np.zeros(6 * nn)
        self.pres = np.zeros(6 * nn)
        for m in self.members.values():
            m.points, m.dists = [], []
        for case, fac in self.raw["combinations"][cid]["factors"].items():
            for ld in self.raw["loadCases"][case].get("loads", []):
                t = ld["type"]
                if t == "nodal":
                    for a, name in enumerate(("Fx", "Fy", "Fz", "Mx", "My", "Mz")):
                        self.F[6 * self.idx[ld["node"]] + a] += fac * ld.get(name, 0.0)
                elif t == "supportDisplacement":
                    for a, name in enumerate(DOF):
                        self.pres[6 * self.idx[ld["node"]] + a] += fac * ld.get(name, 0.0)
                elif t == "selfWeight":
                    g = np.array(ld["direction"], float) if "direction" in ld else -np.eye(3)[AX[self.vertical]]
                    g = g / np.linalg.norm(g)
                    for mid in ld.get("members", list(self.members)):
                        m = self.members[mid]
                        w = fac * ld.get("factor", 1.0) * m.gammaA * (m.R @ g)
                        if w.any():
                            m.dists.append((0.0, m.L, w, w))
                else:
                    m = self.members[ld["member"]]
                    dn = ld["direction"]
                    e = np.eye(3)[AX[dn.upper()]]
                    loc = e if dn.islower() else m.R @ e
                    scale = 1.0 if ld.get("positionType", "relative") == "absolute" else m.L
                    if t == "memberPoint":
                        f = np.zeros(6)
                        f[3 if ld.get("kind", "force") == "moment" else 0:][:3] = fac * ld["value"] * loc
                        m.points.append((ld["position"] * scale, f))
                    else:
                        x1 = ld["start"] * scale if "start" in ld else 0.0
                        x2 = ld["end"] * scale if "end" in ld else m.L
                        pr = 1.0
                        if ld.get("projected", False):
                            pr = np.linalg.norm(np.cross(m.R[0], e))
                        w1 = fac * pr * ld["w1"] * loc
                        w2 = fac * pr * ld.get("w2", ld["w1"]) * loc
                        m.dists.append((x1, x2, w1, w2))

    # ------------------------------------------------------------------ system
    def _G(self, m, end):
        """(6, 6): local displacements of the member end from the 6 global node DOFs."""
        o = m.oi if end == 0 else m.oj
        E = np.eye(6)
        E[:3, 3:] = -cross_matrix(o)
        T = np.zeros((6, 6))
        T[:3, :3] = T[3:, 3:] = m.R
        return T @ E

    def assemble(self, n_scale=None):
        """Matrix and right-hand side.  ``n_scale``: buckling factor on the stored member axial forces."""
        nn = len(self.node_ids)
        mem = [m for m in self.members.values() if m.active]
        nU = 6 * nn
        col_m = {m.id: nU + 18 * k for k, m in enumerate(mem)}        # d0(6) d1(6) s0(6)
        fixed = np.flatnonzero(self.fixed)
        nR, nC = len(fixed), len(self.crows)
        cR = nU + 18 * len(mem)
        cC = cR + nR
        n = cC + nC
        M, b = np.zeros((n, n)), np.zeros(n)
        b[:nU] = -self.F
        M[np.arange(nU), np.arange(nU)] -= self.spring          # node equilibrium rows 0..nU-1
        for k, g in enumerate(fixed):
            M[g, cR + k] = 1.0                                   # reaction on the node
            M[cR + k, g] = 1.0                                   # u = imposed
            b[cR + k] = self.pres[g]
        for k, row in enumerate(self.crows):
            M[:nU, cC + k] += row
            M[cC + k, :nU] = row
        row = nU
        self._maps = {}
        for m in mem:
            c = col_m[m.id]
            if n_scale is not None:
                keep = m.Nstart
                m.Nstart = (keep or 0.0) * n_scale
            *_, (_, _, Phi, p) = m.march(nsub=self.nsub)
            Nax = (m.axial(0.0), m.axial(m.L)) if m.Nstart is not None else (0.0, 0.0)
            if n_scale is not None:
                m.Nstart = keep
            self._maps[m.id] = (Phi, p)
            # 1. compatibility: d1 - Phi_dd d0 - Phi_ds s0 = p_d
            M[row:row + 6, c + 6:c + 12] = np.eye(6)
            M[row:row + 6, c:c + 6] = -Phi[:6, :6]
            M[row:row + 6, c + 12:c + 18] = -Phi[:6, 6:]
            b[row:row + 6] = p[:6]
            row += 6
            for end, node in ((0, m.ni), (1, m.nj)):
                gn = 6 * self.idx[node]
                G = self._G(m, end)
                # end force on the member Q = Qd d0 + Qs s0 + q
                if end == 0:
                    Qd, Qs, q = np.zeros((6, 6)), -np.eye(6), np.zeros(6)
                else:
                    Qd, Qs, q = Phi[6:, :6], Phi[6:, 6:], p[6:]
                # 3. node equilibrium: the node receives -G^T Q
                M[gn:gn + 6, c:c + 6] -= G.T @ Qd
                M[gn:gn + 6, c + 12:c + 18] -= G.T @ Qs
                b[gn:gn + 6] += G.T @ q
                if self.second_order and self.arm:
                    # second-order rotation of the rigid arm under the axial end force f (THEORY 10.7):
                    # the node carries the extra stiffness 1/2 (o f^T + f o^T) - (o.f) I on its rotations
                    o = m.oi if end == 0 else m.oj
                    f = m.R[0] * (-Nax[0] if end == 0 else Nax[1])
                    M[gn + 3:gn + 6, gn + 3:gn + 6] -= 0.5 * (np.outer(o, f) + np.outer(f, o)) - (o @ f) * np.eye(3)
                # 2. connections
                for a in range(6):
                    k = m.k[6 * end + a]
                    dcol = c + 6 * end + a
                    if k == np.inf:
                        M[row, dcol] = 1.0
                        M[row, gn:gn + 6] = -G[a]
                    else:
                        M[row, c:c + 6] = Qd[a]
                        M[row, c + 12:c + 18] = Qs[a]
                        b[row] = -q[a]
                        if k > 0:
                            M[row, gn:gn + 6] -= k * G[a]
                            M[row, dcol] += k
                    row += 1
        # node DOFs that nothing touches, and the rotations of nodes with truss-type members only (a truss
        # member has no rotational end connection; its twist is held at end i only to fix its own spin)
        att = {nid: [] for nid in self.node_ids}
        for m in self.members.values():
            att[m.ni].append(m)
            att[m.nj].append(m)
        slave = np.zeros(nU, bool)
        for r_ in self.crows:
            slave |= r_ != 0
        for g in range(nU):
            nid = self.node_ids[g // 6]
            truss_rot = (g % 6 >= 3 and att[nid] and not self.spring[g] and not slave[g] and
                         all(m.type != "frame" and not (m.oi.any() or m.oj.any()) for m in att[nid]))
            if not self.fixed[g] and (truss_rot or (not np.any(M[nU:, g]) and not np.any(M[:nU, g]))):
                M[g, :] = 0.0
                M[g, g] = 1.0
                b[g] = 0.0
        self._layout = (nU, col_m, cR, fixed, mem)
        return M, b

    @staticmethod
    def _equilibrate(M):
        r, c = np.ones(len(M)), np.ones(len(M))
        A = M.copy()
        for _ in range(12):
            rr = 1 / np.sqrt(np.maximum(np.abs(A).max(axis=1), 1e-300))
            A *= rr[:, None]
            cc = 1 / np.sqrt(np.maximum(np.abs(A).max(axis=0), 1e-300))
            A *= cc[None, :]
            r, c = r * rr, c * cc
        return A, r, c

    def _solve_linear(self):
        M, b = self.assemble()
        A, r, c = self._equilibrate(M)
        lu = lu_factor(A)
        y = lu_solve(lu, r * b)
        for _ in range(3):
            y += lu_solve(lu, r * b - A @ y)
        x = c * y
        self.residual = np.abs(A @ y - r * b).max() / max(np.abs(r * b).max(), 1e-300)
        self.rcond = 1.0 / np.linalg.cond(A)
        return x

    def solve(self, cid, second_order=False, stations=None, tol=1e-11, max_iter=200):
        """Solve one combination.  ``stations``: {member: [(x, side), ...]} (file length unit)."""
        self.set_combination(cid)
        self.second_order = second_order
        for m in self.members.values():
            m.active, m.Nstart = True, None
        slack_hist = []
        while True:
            for m in self.members.values():
                m.Nstart = None
            x = self._solve_linear()
            if second_order:
                for it in range(max_iter):
                    nU, col_m, *_ = self._layout
                    new = {mid: x[c + 12] for mid, c in col_m.items()}
                    old = {mid: self.members[mid].Nstart or 0.0 for mid in new}
                    for mid, v in new.items():
                        self.members[mid].Nstart = v
                    scale = max(max(abs(v) for v in new.values()), 1e-300)
                    if it > 0 and max(abs(new[k] - old[k]) for k in new) < tol * scale:
                        break
                    x = self._solve_linear()
                else:
                    raise RuntimeError("second-order iteration of the reference did not converge")
                self.iterations = it
            nU, col_m, cR, fixed, mem = self._layout
            U = x[:nU].reshape(-1, 6)
            # one-sided members
            flip = []
            for m in self.members.values():
                if m.type not in ("tensionOnly", "compressionOnly"):
                    continue
                sign = 1.0 if m.type == "tensionOnly" else -1.0
                if m.active:
                    if sign * x[col_m[m.id] + 12] < 0:
                        flip.append(m)
                else:
                    ui = self._G(m, 0) @ U[self.idx[m.ni]]
                    uj = self._G(m, 1) @ U[self.idx[m.nj]]
                    if sign * (uj[0] - ui[0]) > 0:
                        flip.append(m)
            if not flip:
                break
            for m in flip:
                m.active = not m.active
            state = tuple(m.active for m in self.members.values())
            if state in slack_hist:
                raise RuntimeError("active set of the reference cycles")
            slack_hist.append(state)
        R = np.zeros(nU)
        R[fixed] = x[cR:cR + len(fixed)]
        R -= self.spring * x[:nU]
        out = {"displacements": U, "reactions": R.reshape(-1, 6), "members": {}}
        for m in mem:
            c = col_m[m.id]
            z0 = np.concatenate([x[c:c + 6], x[c + 12:c + 18]])
            st = (stations or {}).get(m.id, [])
            states = [(xs, side, Phi @ z0 + p) for xs, side, Phi, p in m.march(st, nsub=self.nsub)]
            zL = states[-1][2]
            # station list in the requested order
            lookup = {(xs, side): z for xs, side, z in states[:-1]}
            ordered = np.array([lookup[(float(xs), side)] for xs, side in st]).reshape(-1, 12)
            out["members"][m.id] = {"end_forces": np.concatenate([-z0[6:], zL[6:]]),
                                    "forces": ordered[:, 6:], "displacements": ordered[:, :3],
                                    "state": ordered}
        for m in self.members.values():
            if not m.active:
                out["members"][m.id] = {"end_forces": np.zeros(12), "slack": True}
        return out

    # ------------------------------------------------------------------ buckling
    def buckling(self, cid, n_modes=1, lam_max=50.0, steps=400, negative=False):
        """Smallest load factors at which the tableau (N = factor x first-order N) is singular."""
        self.solve(cid, second_order=False)
        x = self._solve_linear()
        nU, col_m, *_ = self._layout
        for mid, c in col_m.items():
            self.members[mid].Nstart = x[c + 12]
        self.second_order = True
        M0, _ = self.assemble(0.0)
        _, r, c = self._equilibrate(M0)

        def det(lam):
            M, _ = self.assemble(lam)
            s, logdet = np.linalg.slogdet(M * r[:, None] * c[None, :])
            return s, logdet
        s0, _ = det(0.0)
        f = lambda lam: det(lam)[0] * s0
        roots, prev, a = [], 1.0, 0.0
        sgn = -1.0 if negative else 1.0
        for k in range(1, steps + 1):
            bq = sgn * lam_max * k / steps
            cur = f(bq)
            if cur * prev < 0:
                roots.append(brentq(f, a, bq, xtol=1e-12, rtol=1e-12))
                if len(roots) == n_modes:
                    break
            prev, a = cur, bq
        self.second_order = False
        return roots
