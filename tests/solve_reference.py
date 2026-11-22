"""Helpers for the solver tests.

* ``make(...)``: a small model document in N and mm (file units = internal units).
* ``beam_bvp(...)``: an INDEPENDENT reference for one bending plane of a straight member: the governing
  first-order system of docs/CONVENTIONS.md section 4 is solved as a two-point boundary value problem
  by Chebyshev collocation (one sub-interval per load segment, jump conditions at point loads).  The
  exact solutions are piecewise polynomials of degree <= 5 and the collocation polynomials have
  degree 9, so the reference is accurate to round-off.  It shares no code with ``benchfea.solve``.
* ``beam_bvp_scipy(...)``: the same problem handed to ``scipy.integrate.solve_bvp`` (adaptive cubic
  collocation, residual tolerance 1e-9): a second, slower witness used to check the first.
* ``bar_bvp(...)``: the axial / torsional first-order system (Chebyshev collocation).
"""
import numpy as np
from scipy.integrate import solve_bvp

E, G = 200000.0, 77000.0
SEC = {"A": 5000.0, "Iy": 8.0e6, "Iz": 4.0e7, "J": 3.0e5, "Asy": 2500.0, "Asz": 3000.0}
EA, EIy, EIz, GJ = E * SEC["A"], E * SEC["Iy"], E * SEC["Iz"], G * SEC["J"]
GAy, GAz = G * SEC["Asy"], G * SEC["Asz"]


def make(nodes, members, supports, cases, combos=None, analysis=None, constraints=None, sections=None,
         materials=None, vertical="Z", design=None):
    """Model document. nodes: {id: (x, y, z)}; members: {id: (i, j) | dict}; supports: {id: "fixed" | dict};
    cases: {id: [load, ...]}; combos default: one combination per case, factor 1, same id."""
    doc = {
        "formatVersion": "1.0.0", "units": {"force": "N", "length": "mm"},
        "conventions": {"verticalAxis": vertical},
        "nodes": {k: {"x": float(v[0]), "y": float(v[1]), "z": float(v[2])} for k, v in nodes.items()},
        "supports": {k: ({"type": v} if isinstance(v, str) else v) for k, v in supports.items()},
        "materials": materials or {"M": {"E": E, "G": G, "unitWeight": 7.7e-5}},
        "sections": sections or {"S": {"type": "properties", **SEC}},
        "members": {}, "loadCases": {k: {"loads": list(v)} for k, v in cases.items()},
        "combinations": combos or {k: {"factors": {k: 1.0}} for k in cases},
    }
    for k, v in members.items():
        m = {"i": v[0], "j": v[1]} if isinstance(v, (tuple, list)) else dict(v)
        m.setdefault("section", next(iter(doc["sections"])))
        m.setdefault("material", next(iter(doc["materials"])))
        doc["members"][k] = m
    if analysis:
        doc["analysis"] = analysis
    if constraints:
        doc["constraints"] = constraints
    if design:
        doc["design"] = design
    return doc


N_CHEB = 9


def _cheb(n=N_CHEB):
    """Chebyshev-Lobatto nodes on [0, 1], differentiation matrix and barycentric weights."""
    j = np.arange(n + 1)
    t = 0.5 * (1.0 - np.cos(np.pi * j / n))
    wts = (-1.0) ** j
    wts[0] *= 0.5
    wts[-1] *= 0.5
    D = np.zeros((n + 1, n + 1))
    for i in range(n + 1):
        for k in range(n + 1):
            if i != k:
                D[i, k] = (wts[k] / wts[i]) / (t[i] - t[k])
        D[i, i] = -D[i].sum()
    return t, D, wts


def _interp(t, wts, values, x):
    d = x - t
    if np.any(d == 0.0):
        return values[int(np.argmin(np.abs(d)))]
    c = wts / d
    return (c @ values) / c.sum()


def _affine(fun, phys, n):
    """Boundary conditions are affine in the state: residual = r0 + J z; rows scaled to unit size."""
    r0 = np.asarray(fun(phys(np.zeros(n))), float)
    J = np.array([np.asarray(fun(phys(e)), float) - r0 for e in np.eye(n)]).T
    sc = np.abs(J).max(axis=1)
    return r0 / sc, J / sc[:, None]


def _collocate(xs, A, g, bc0, bcL, jump, phys):
    """Solve z' = A z + g(k, x) on the segments xs[k]..xs[k+1] (z continuous except for the jumps
    ``z_right - z_left = -jump[x]`` at interior break points) with two-point boundary conditions.
    Returns f(x, side) -> physical state."""
    ns = A.shape[0]
    m = len(xs) - 1
    t, D, wts = _cheb()
    n1 = len(t)
    size = ns * n1 * m
    M = np.zeros((size, size))
    rhs = np.zeros(size)

    def col(k, i, c):                     # unknown: state c at node i of segment k
        return (k * n1 + i) * ns + c
    row = 0
    for k in range(m):
        h = xs[k + 1] - xs[k]
        for i in range(1, n1):            # collocation at every node but the first
            gi = h * np.asarray(g(k, xs[k] + t[i] * h), float)
            for c in range(ns):
                for jn in range(n1):
                    M[row, col(k, jn, c)] += D[i, jn]
                for c2 in range(ns):
                    M[row, col(k, i, c2)] -= h * A[c, c2]
                rhs[row] = gi[c]
                row += 1
    half = ns // 2
    r0, J = _affine(bc0, phys, ns)
    for q in range(half):
        M[row, [col(0, 0, c) for c in range(ns)]] = J[q]
        rhs[row] = -r0[q]
        row += 1
    r0, J = _affine(bcL, phys, ns)
    for q in range(half):
        M[row, [col(m - 1, n1 - 1, c) for c in range(ns)]] = J[q]
        rhs[row] = -r0[q]
        row += 1
    for k in range(m - 1):
        jv = jump(xs[k + 1])
        for c in range(ns):
            M[row, col(k + 1, 0, c)] = 1.0
            M[row, col(k, n1 - 1, c)] = -1.0
            rhs[row] = -jv[c]
            row += 1
    assert row == size
    z = np.linalg.solve(M, rhs).reshape(m, n1, ns)

    def f(x, side=-1):
        x = float(x)
        k = int(np.searchsorted(xs, x, side="left" if side < 0 else "right")) - 1
        k = min(max(k, 0), m - 1)
        return phys(_interp(t, wts, z[k], (x - xs[k]) / (xs[k + 1] - xs[k])))
    return f


def _breaks(L, dist, points):
    return sorted({0.0, float(L)} | {float(p[0]) for p in points if 0 < p[0] < L}
                  | {float(b) for d in dist for b in d[:2] if 0 < b < L})


def beam_bvp(L, EI, GAs, bc0, bcL, dist=(), points=(), plane="xy"):
    """Bending in one plane.  State y = (d, th, M, V):

        plane "xy":  v' = th_z + Vy/GAs,   th_z' = Mz/EIz,  Mz' = -Vy,  Vy' = -w_y
        plane "xz":  w' = -th_y + Vz/GAs,  th_y' = My/EIy,  My' = +Vz,  Vz' = -w_z

    ``dist`` = [(x1, x2, w1, w2)]: line loads varying linearly from w1 at x1 to w2 at x2.
    ``points`` = [(x, P, Mo)]: at a point force P the
    shear drops by P, at a point moment Mo the moment drops by Mo (resultants on the +x face).
    ``bc0(y)``, ``bcL(y)`` return two residuals each (affine in y); they see the state just inside the
    member, i.e. after the loads at x = 0 and before the loads at x = L -- loads exactly at the ends
    must be put in the boundary conditions by the caller.  ``GAs = None``: no shear deformation.
    Returns ``f(x, side=-1) -> (d, th, M, V)``.
    """
    s = 1.0 if plane == "xy" else -1.0
    xs = _breaks(L, dist, points)
    w = _line_load(xs, dist)
    # scaled state (all lengths): z = (d, th L, M L^2 / EI, V L^3 / EI)
    cM, cV = L ** 2 / EI, L ** 3 / EI
    sh = EI / (L ** 3 * GAs) if GAs else 0.0
    A = np.array([[0.0, s / L, 0.0, sh], [0.0, 0.0, 1.0 / L, 0.0], [0.0, 0.0, 0.0, -s / L], [0.0, 0.0, 0.0, 0.0]])
    jumps = {}
    for x, P, Mo in points:
        jv = jumps.setdefault(float(x), np.zeros(4))
        jv += [0.0, 0.0, Mo * cM, P * cV]
    return _collocate(xs, A, lambda k, x: [0.0, 0.0, 0.0, -w(k, x) * cV], bc0, bcL,
                      lambda x: jumps.get(x, np.zeros(4)), lambda z: np.array([z[0], z[1] / L, z[2] / cM, z[3] / cV]))


def bar_bvp(L, EAx, bc0, bcL, dist=(), points=()):
    """Axial (or torsional) member: u' = N/EA, N' = -w_x; at a point load P the force drops by P.
    Returns ``f(x, side=-1) -> (u, N)``."""
    xs = _breaks(L, dist, points)
    w = _line_load(xs, dist)
    A = np.array([[0.0, 1.0 / L], [0.0, 0.0]])          # scaled state z = (u, N L / EA)
    jumps = {}
    for x, P in points:
        jumps[float(x)] = jumps.get(float(x), np.zeros(2)) + [0.0, P * L / EAx]
    return _collocate(xs, A, lambda k, x: [0.0, -w(k, x) * L / EAx], bc0, bcL,
                      lambda x: jumps.get(x, np.zeros(2)), lambda z: np.array([z[0], z[1] * EAx / L]))


def beam_bvp_scipy(L, EI, GAs, bc0, bcL, dist=(), points=(), plane="xy", tol=1e-9):
    """Same problem and interface as :func:`beam_bvp`, solved by ``scipy.integrate.solve_bvp``."""
    s = 1.0 if plane == "xy" else -1.0
    xs = _breaks(L, dist, points)
    w = _line_load(xs, dist)
    m = len(xs) - 1
    h = np.diff(xs)
    jump = {}
    for x, P, Mo in points:
        jv = jump.setdefault(float(x), [0.0, 0.0])
        jv[0] += P
        jv[1] += Mo
    cM, cV = L ** 2 / EI, L ** 3 / EI
    sh = EI / (L ** 3 * GAs) if GAs else 0.0

    def rhs(t, z):
        out = np.zeros_like(z)
        for k in range(m):
            d, th, M, V = z[4 * k:4 * k + 4]
            x = xs[k] + t * h[k]
            out[4 * k] = h[k] * (s * th / L + V * sh)
            out[4 * k + 1] = h[k] * M / L
            out[4 * k + 2] = h[k] * (-s * V) / L
            out[4 * k + 3] = h[k] * (-w(k, x)) * cV
        return out

    def phys(z):
        return np.array([z[0], z[1] / L, z[2] / cM, z[3] / cV])

    r00, J0 = _affine(bc0, phys, 4)
    r0L, JL = _affine(bcL, phys, 4)

    def bc(za, zb):
        res = list(r00 + J0 @ za[:4]) + list(r0L + JL @ zb[4 * (m - 1):])
        for k in range(m - 1):
            P, Mo = jump.get(xs[k + 1], (0.0, 0.0))
            left, right = zb[4 * k:4 * k + 4], za[4 * k + 4:4 * k + 8]
            res += [right[0] - left[0], right[1] - left[1], right[2] - left[2] + Mo * cM, right[3] - left[3] + P * cV]
        return np.array(res)

    t = np.linspace(0.0, 1.0, 41)
    sol = solve_bvp(rhs, bc, t, np.zeros((4 * m, t.size)), tol=tol, max_nodes=50000)
    assert sol.success, sol.message

    def f(x, side=-1):
        x = float(x)
        k = int(np.searchsorted(xs, x, side="left" if side < 0 else "right")) - 1
        k = min(max(k, 0), m - 1)
        return phys(sol.sol((x - xs[k]) / h[k])[4 * k:4 * k + 4])
    return f


def _line_load(xs, dist):
    """w(k, x) on segment k = [xs[k], xs[k+1]]: sum of the line loads that cover the whole segment."""
    def w(k, x):
        a, b = xs[k], xs[k + 1]
        out = 0.0 * x
        for x1, x2, w1, w2 in dist:
            if x1 <= a and b <= x2:
                out = out + w1 + (w2 - w1) * (x - x1) / (x2 - x1)
        return out
    return w


# boundary conditions for beam_bvp (state = (d, th, M, V) just inside the member)
def clamped(y):
    return [y[0], y[1]]


def pinned(y):
    return [y[0], y[2]]


def free(y):
    return [y[2], y[3]]


def guided(y):
    return [y[1], y[3]]
