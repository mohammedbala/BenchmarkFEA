"""Independent review: single members by direct numerical integration of the governing ODE as a
two-point boundary value problem (``scipy.integrate.solve_bvp``) -- no transfer matrix, no stiffness.

One bending plane of a prismatic member, state (v, theta, M, V) with the conventions of
docs/CONVENTIONS.md section 4 (derived in ``tableau.py``): with the axial force N (tension +) acting on
the displaced slice and shear flexibility beta = 1/(G As),

    v' = (theta + beta V) / (1 + beta N)     theta' = M / EI
    M' = (N theta - V) / (1 + beta N)        V' = -q(x)

(N = 0: first order).  The engine and the tableau reference are both compared with it.
"""
import numpy as np
import pytest
from scipy.integrate import solve_bvp

from harness import case, compare, member, nodes, run_benchfea, worst
from models import start

L, EI, GAS = 5.0, 2.0e8 * 1.5e-4, 7.7e7 * 0.0012


def bvp(N, beta, loads, end):
    """Clamped at x = 0; ``end`` = 'pinned' (v = 0, M = 0) or 'spring' (v = 0, M = -k theta, k = 9000).

    The load is piecewise linear, so the member is cut at the load break points; every smooth piece is
    mapped onto t in [0, 1] and all pieces are solved together, with continuity of (v, theta, M, V) at
    the cuts as boundary conditions.  States are scaled to (v/L, theta, M L/EI, V L^2/EI).
    Returns a function x -> (v, theta, M, V)."""
    breaks = sorted({0.0, L} | {v for ld in loads for v in ld[:2]})
    segs = list(zip(breaks[:-1], breaks[1:]))
    ns = len(segs)
    S = np.array([1.0 / L, 1.0, L / EI, L * L / EI])

    def q(x, a, b):                                       # line load on the open piece (a, b)
        out = np.zeros_like(x)
        for x1, x2, w1, w2 in loads:
            if x1 <= a and b <= x2:
                out += w1 + (w2 - w1) * (x - x1) / (x2 - x1)
        return out

    def rhs(t, y):
        out = np.zeros_like(y)
        for k, (a, b) in enumerate(segs):
            v, th, M, V = y[4 * k:4 * k + 4] / S[:, None]
            x = a + t * (b - a)
            f = np.vstack([(th + beta * V) / (1 + beta * N), M / EI, (N * th - V) / (1 + beta * N), -q(x, a, b)])
            out[4 * k:4 * k + 4] = (b - a) * f * S[:, None]
        return out

    def bc(ya, yb):
        res = [ya[0], ya[1]]
        for k in range(ns - 1):
            res += list(yb[4 * k:4 * k + 4] - ya[4 * k + 4:4 * k + 8])
        last = yb[4 * (ns - 1):] / S
        res += [yb[4 * (ns - 1)], (last[2] + (9000.0 * last[1] if end == "spring" else 0.0)) * S[2]]
        return np.array(res)
    t = np.linspace(0.0, 1.0, 41)
    sol = solve_bvp(rhs, bc, t, np.zeros((4 * ns, len(t))), tol=1e-11, max_nodes=20000)
    assert sol.success and sol.rms_residuals.max() < 1e-9, sol.message

    def at(x):
        out = np.zeros((4, len(x)))
        for n, xv in enumerate(x):
            k = min(max(int(np.searchsorted(breaks, xv, side="right")) - 1, 0), ns - 1)
            a, b = segs[k]
            out[:, n] = sol.sol((xv - a) / (b - a))[4 * k:4 * k + 4] / S
        return out
    return at


@pytest.mark.parametrize("P, shear, end", [(0.0, False, "pinned"), (0.0, True, "pinned"), (2500.0, False, "pinned"),
                                           (2500.0, True, "spring"), (-4000.0, True, "pinned"), (6000.0, True, "spring")])
def test_propped_member_with_partial_trapezoid_against_solve_bvp(P, shear, end):
    """P > 0 compression (Euler load of the pinned-pinned member 11 800 kN), P < 0 tension."""
    loads = [(0.6, 3.1, -12.0, -40.0), (3.1, 5.0, 25.0, 25.0)]
    d = start(shear=shear, convergenceTolerance=1.0e-12)
    d["sections"]["S"].update({"Asy": 0.0012})
    d["sections"]["S"].pop("Asz", None)
    nodes(d, A=(0, 0, 0), B=(L, 0, 0))
    d["supports"] = {"A": {"type": "fixed"}, "B": {"dofs": {"uy": "fixed", "uz": "fixed", "rx": "fixed", "rz": "fixed",
                                                              **({"ry": "fixed"} if end == "spring" else {})}}}
    rel = {"releases": {"j": {"rz": {"spring": 9000.0}}}} if end == "spring" else {}
    member(d, "m", "A", "B", **rel)
    case(d, "L", {"type": "nodal", "node": "B", "Fx": -P},
         *[{"type": "memberDistributed", "member": "m", "direction": "y", "w1": w1, "w2": w2, "start": x1, "end": x2,
            "positionType": "absolute"} for x1, x2, w1, w2 in loads],
         analysisType="secondOrder" if P else "firstOrder")
    _, res, doc = run_benchfea(d)
    m = doc["combinations"]["L"]["members"]["m"]
    x = np.array(m["stations"]["x"])
    sol = bvp(-P, 1.0 / GAS if shear else 0.0, loads, end)
    v, th, M, V = sol(np.minimum(x, L))
    scale_m, scale_v, scale_d = np.abs(M).max(), np.abs(V).max(), np.abs(v).max()
    assert np.abs(np.array(m["stations"]["Mz"]) - M).max() <= 1e-7 * scale_m
    assert np.abs(np.array(m["stations"]["Vy"]) - V).max() <= 1e-7 * scale_v
    assert np.abs(np.array(m["stations"]["uy"]) - v).max() <= 1e-7 * scale_d
    assert np.allclose(m["stations"]["N"], -P, atol=1e-9 * max(abs(P), 1.0))
    # and the tableau reference reproduces the same member
    assert worst(compare(d, "L", second_order=bool(P), doc=doc)) <= 1e-9
