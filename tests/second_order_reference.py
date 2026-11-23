"""INDEPENDENT references for the second-order and buckling tests.  Nothing here imports
``benchfea.solve``.

* ``beam_column_bvp``: one bending plane of a straight member under a given axial force, solved as a
  two-point boundary value problem with ``scipy.integrate.solve_bvp`` (adaptive collocation, residual
  tolerance 1e-10).  Governing equations (docs/THEORY.md section 10; N tension positive, V = force
  along the UNDEFORMED y axis on the +x face, beta = 1 / (G As)):

      v' = (theta + beta V) / (1 + beta N)     theta' = M / EI
      M' = -V + N v'                           V' = -q(x)

  The member may consist of several segments (point loads / supports in between): the segments are
  stacked on a common unit interval and coupled by interface conditions.
* ``euler_bernoulli_bvp``: the same problem without shear deformation in its classical fourth-order
  form EI v'''' - N v'' = q, i.e. a differently written system (state v, v', v'', v'''), used to check
  the first reference.
* ``critical_k``: smallest root of a characteristic determinant.  The general solution of
  EI v'''' + P v'' = 0 is v = A + B x + C sin(kx) + D cos(kx), k^2 = P / EI (check: substitute).  The
  four boundary conditions give a 4x4 homogeneous system in (A, B, C, D); the critical load is the
  smallest k for which its determinant vanishes, found numerically (sign change + brentq).  No
  effective-length factor is quoted anywhere.
"""
import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq


def beam_column_bvp(lengths, EI, GAs, N, q, bc, tol=1e-10):
    """Segments of lengths ``lengths`` (one number or a list).  ``N`` and ``q`` are functions of the
    distance x from the member start, or numbers, or a list with one function / number per segment
    (so that they may jump at a segment boundary).  ``bc(start, end)``
    receives the physical states (v, theta, M, V) of every segment at its start and its end, as two
    arrays of shape (n_segments, 4), and returns 4 * n_segments residuals.
    Returns f(x) -> (v, theta, M, V), with ``f(x, right=True)`` selecting the segment that starts at x."""
    lengths = [float(lengths)] if np.isscalar(lengths) else [float(x) for x in lengths]
    m = len(lengths)
    starts = np.concatenate([[0.0], np.cumsum(lengths)])
    L = float(starts[-1])

    def per_segment(f):
        fs = list(f) if isinstance(f, (list, tuple)) else [f] * m
        return [g if callable(g) else (lambda x, c=float(g): c + 0.0 * x) for g in fs]
    Nf, qf = per_segment(N), per_segment(q)
    beta = 1.0 / GAs if GAs else 0.0
    S = np.array([L, 1.0, EI / L, EI / L ** 2])              # physical = S * scaled

    def rhs(s, y):
        out = np.zeros_like(y)
        for k in range(m):
            x = starts[k] + lengths[k] * s
            v, th, M, V = (y[4 * k + c] * S[c] for c in range(4))
            n = Nf[k](x)
            dv = (th + beta * V) / (1.0 + beta * n)
            d = [dv, M / EI, -V + n * dv, -qf[k](x)]
            for c in range(4):
                out[4 * k + c] = lengths[k] * d[c] / S[c]
        return out

    def res(ya, yb):
        a = ya.reshape(m, 4) * S
        b = yb.reshape(m, 4) * S
        r = np.asarray(bc(a, b), float)
        return r

    # scale the residuals so that the Newton iteration is well conditioned
    probe = np.ones((m, 4)) * S
    unit = np.zeros(4 * m)
    for k in range(4 * m):
        e = np.zeros(4 * m)
        e[k] = 1.0
        d = np.abs(np.asarray(bc(e.reshape(m, 4) * S, np.zeros((m, 4))), float)
                   - np.asarray(bc(np.zeros((m, 4)), np.zeros((m, 4))), float))
        d2 = np.abs(np.asarray(bc(np.zeros((m, 4)), e.reshape(m, 4) * S), float)
                    - np.asarray(bc(np.zeros((m, 4)), np.zeros((m, 4))), float))
        unit = np.maximum(unit, np.maximum(d, d2))
    del probe
    unit[unit == 0.0] = 1.0

    s0 = np.linspace(0.0, 1.0, 41)
    sol = solve_bvp(rhs, lambda ya, yb: res(ya, yb) / unit, s0, np.zeros((4 * m, len(s0))), tol=tol, max_nodes=200000)
    assert sol.success, sol.message

    def f(x, right=False):
        x = float(x)
        k = int(np.searchsorted(starts, x, side="right" if right else "left")) - 1
        k = min(max(k, 0), m - 1)
        return sol.sol((x - starts[k]) / lengths[k])[4 * k:4 * k + 4] * S
    return f


def euler_bernoulli_bvp(L, EI, N, q, bc, tol=1e-10):
    """EI v'''' - N v'' = q (constant N, no shear deformation).  ``bc(a, b)`` receives (v, v', v'', v''')
    at the two ends and returns 4 residuals.  Returns f(x) -> (v, v', v'', v''')."""
    qf = q if callable(q) else (lambda x, c=float(q): c + 0.0 * x)
    S = np.array([L, 1.0, 1.0 / L, 1.0 / L ** 2])

    def rhs(s, y):
        return np.vstack([y[1], y[2], y[3], (N * y[2] * S[2] + qf(L * s)) / EI * L / S[3]])

    zero = np.zeros(4)
    unit = np.zeros(4)
    for k in range(4):
        e = np.zeros(4)
        e[k] = S[k]
        base = np.asarray(bc(zero, zero), float)
        unit = np.maximum(unit, np.maximum(np.abs(np.asarray(bc(e, zero), float) - base),
                                           np.abs(np.asarray(bc(zero, e), float) - base)))
    unit[unit == 0.0] = 1.0
    s0 = np.linspace(0.0, 1.0, 41)
    sol = solve_bvp(rhs, lambda ya, yb: np.asarray(bc(ya * S, yb * S), float) / unit, s0, np.zeros((4, len(s0))),
                    tol=tol, max_nodes=200000)
    assert sol.success, sol.message
    return lambda x: sol.sol(float(x) / L) * S


def column_shape(k, x):
    """Rows: v, v', v'', v''' of the basis (1, x, sin kx, cos kx) at x."""
    s, c = np.sin(k * x), np.cos(k * x)
    return np.array([[1.0, x, s, c],
                     [0.0, 1.0, k * c, -k * s],
                     [0.0, 0.0, -k * k * s, -k * k * c],
                     [0.0, 0.0, -k ** 3 * c, k ** 3 * s]])


def critical_k(rows, k_max, n_scan=4000, count=1):
    """Smallest k in (0, k_max) at which det([rows(k)]) changes sign (``count`` > 1: the first ``count``
    such k as a list).  ``rows(k)`` returns the 4x4 matrix of the boundary conditions applied to the
    basis (1, x, sin kx, cos kx)."""
    def det(k):
        M = np.asarray(rows(k), float)
        M = M / np.abs(M).max(axis=1, keepdims=True)
        return np.linalg.det(M)
    ks = np.linspace(k_max / n_scan, k_max, n_scan)
    vals = [det(k) for k in ks]
    found = []
    for a, b, fa, fb in zip(ks[:-1], ks[1:], vals[:-1], vals[1:]):
        if fa * fb < 0.0:
            found.append(brentq(det, a, b, xtol=1e-15, rtol=1e-14))
            if len(found) == count:
                return found[0] if count == 1 else found
    raise AssertionError("not enough roots of the characteristic determinant below k_max")


def critical_load_bvp(L, EI, GAs, bc, P_guess, v_guess, norm, tol=1e-10):
    """Critical load of a (shear-deformable) column as an eigenvalue of the homogeneous beam-column
    equations with N = -P: ``solve_bvp`` with P as the unknown parameter.  ``bc(start, end)`` -> 4
    residuals on the physical states (v, theta, M, V); ``norm`` = (``"start"`` | ``"end"``, component)
    names the scaled state component (v / L, theta, M L / EI, V L^2 / EI) that is set to 1 to fix the
    amplitude of the mode; ``v_guess(s)`` is a first guess of the mode shape on s = x / L (it selects
    the mode) and ``P_guess`` of the load."""
    beta = 1.0 / GAs if GAs else 0.0
    S = np.array([L, 1.0, EI / L, EI / L ** 2])

    def rhs(s, y, p):
        P = p[0] * EI / L ** 2
        th, M, V = y[1] * S[1], y[2] * S[2], y[3] * S[3]
        dv = (th + beta * V) / (1.0 - beta * P)
        d = [dv, M / EI, -V - P * dv, np.zeros_like(dv)]
        return np.vstack([L * d[c] / S[c] for c in range(4)])

    zero = np.zeros(4)
    unit = np.zeros(4)
    for k in range(4):
        e = np.zeros(4)
        e[k] = S[k]
        unit = np.maximum(unit, np.maximum(np.abs(np.asarray(bc(e, zero), float)), np.abs(np.asarray(bc(zero, e), float))))
    unit[unit == 0.0] = 1.0
    end, comp = norm
    s0 = np.linspace(0.0, 1.0, 81)
    v0 = np.asarray(v_guess(s0), float)
    th0 = np.gradient(v0, s0)
    m0 = np.gradient(th0, s0)
    y0 = np.vstack([v0, th0, m0, -np.gradient(m0, s0) - P_guess * L ** 2 / EI * th0])
    y0 = y0 / (y0[comp, 0] if end == "start" else y0[comp, -1])

    def res(ya, yb, p):
        r = list(np.asarray(bc(ya * S, yb * S), float) / unit)
        r.append((ya if end == "start" else yb)[comp] - 1.0)
        return np.array(r)
    sol = solve_bvp(rhs, res, s0, y0, p=[P_guess * L ** 2 / EI], tol=tol, max_nodes=100000)
    assert sol.success, sol.message
    return float(sol.p[0] * EI / L ** 2)


# Rows of the characteristic matrix for the classical end conditions (N = -P, k^2 = P / EI):
# held translation: v = 0; held rotation: v' = 0; free rotation: M = EI v'' = 0;
# free translation: V = -EI v''' - P v' = 0, i.e. v''' + k^2 v' = 0.
HELD_ROWS = {
    "held": lambda k, x: column_shape(k, x)[0],
    "no rotation": lambda k, x: column_shape(k, x)[1],
    "no moment": lambda k, x: column_shape(k, x)[2],
    "no shear": lambda k, x: column_shape(k, x)[3] + k * k * column_shape(k, x)[1],
}
