"""Station statics: the internal forces at every result station recomputed from statics.

Tier C compares node displacements, reactions and member END forces with the foreign solvers; the
forces at the result stations BETWEEN the member ends (the values a largest span moment is read from)
are not part of that comparison.  This check closes the gap without any code of ``benchfea.solve``:

    for every member and every station, the free body between the start of the flexible length and
    the station is in equilibrium under  (a) minus the station forces at x = 0,  (b) the member loads
    of the combination between 0 and the station, taken from the load DEFINITIONS
    (``benchfea.crosscheck.prepare``),  (c) the station forces at the station.

It also requires the station forces at x = 0 and x = L to equal -(end forces at i) and the end forces
at j (docs/CONVENTIONS.md 4), which ties the stations to the end forces that tier C cross-checks.

A station AT a point load holds the limit from the left or the limit from the right.  Besides the two
stations the solver puts there, a regular station can fall on the same point within round-off, and
whether it shows the left or the right limit is decided by an ulp of its position; such a station is
compared with both limits and the nearer one counts.  The stations of an interior point load together
must show BOTH limits (a jump that is missing in the results fails).

First order: the free body is taken on the undeformed axis; the identity is exact (round-off).
Second order (docs/THEORY.md 10.2, 10.8): the axial force components act with the lever of the
deflection of the axis, ``dMz/dx = -Vy + N v'`` and ``dMy/dx = Vz - N w'``; the deflections are the
station displacements of the results (multiplied by ``alpha``: stored results are divided by it).  The
check then states that the stored moments, shears, axial forces AND deflections of a member are
consistent with equilibrium on the deformed axis at every station; it does not recompute the
deflections (``station_deflections`` below does).  With an axial LINE load on a member in a
second-order combination the lever integral needs the deflection between the stations: it is taken by
the trapezoid rule, and unless a bound of that integral is negligible (below a tenth of the tolerance)
the bending moments of that member are left to ``station_deflections`` (listed under ``deferred``),
which integrates the member equations and compares the moments at every station.
"""
from __future__ import annotations

import numpy as np

from ..crosscheck.prepare import prepare

__all__ = ["station_statics", "station_deflections"]

_EPS_X = 1.0e-9                      # two positions closer than this (relative to L) are the same point
#: Scales: the largest station force and the largest station moment of the combination.  Only so that a
#: quantity that is zero in the whole model (the moments of a truss) does not turn round-off into a ratio,
#: each scale is at least this fraction of the other one times / divided by the size of the model.
_FLOOR = 1.0e-6


def _line_load(w1, w2, x1, x2, x):
    """Resultant force of a linearly varying line load on [x1, min(x, x2)] and its first moment
    ``int (s - x) w(s) ds`` (vectors of 3, local axes)."""
    b = min(x, x2)
    if b <= x1:
        return np.zeros(3), np.zeros(3), 0.0
    T, c = b - x1, x1 - x
    g = (w2 - w1) / (x2 - x1)
    force = w1 * T + g * (T * T / 2.0)
    first = w1 * (T * T / 2.0 + c * T) + g * (T ** 3 / 3.0 + c * T * T / 2.0)
    return force, first, T


def station_statics(model, results, tol_linear: float, tol_iterated: float) -> dict:
    """Returns ``{"passed", "failed": [text], "largest": {...}, "combinations": {cid: {...}},
    "stations": n, "notEvaluated": [text], "deferred": [text]}`` (``deferred``: members whose bending
    moments are compared by :func:`station_deflections` instead).  Differences are relative to the
    largest station force (forces) / the largest station moment (moments) of the combination
    (``_FLOOR``: only a quantity that is zero in the whole model gets a floor)."""
    xyz = np.array([nd.xyz for nd in model.nodes.values()], float)
    size = float(np.linalg.norm(xyz.max(axis=0) - xyz.min(axis=0))) if len(xyz) else 1.0
    size = size if size > 0.0 else 1.0
    failed, not_evaluated, deferred, per_combo = [], [], [], {}
    worst_all = {"value": 0.0}
    n_stations = 0
    for cid, cres in results.combinations.items():
        combo = model.combinations[cid]
        second = combo.analysis_type == "secondOrder"
        iterated = second or bool(combo.stiffness_reduction and model.analysis.tau_b == "computed") \
            or any(not mr.active for mr in cres.members.values()) \
            or any(m.type in ("tensionOnly", "compressionOnly") for m in model.members.values())
        tol = tol_iterated if iterated else tol_linear
        alpha = float(getattr(cres, "alpha", 1.0) or 1.0)
        p = prepare(model, cid, tau_b=dict(getattr(cres, "tau_b", {}) or {}), extra_nodal_loads={})
        fs = max((float(np.abs(mr.forces[:, :3]).max(initial=0.0)) for mr in cres.members.values()), default=0.0)
        ms = max((float(np.abs(mr.forces[:, 3:]).max(initial=0.0)) for mr in cres.members.values()), default=0.0)
        fs, ms = max(fs, _FLOOR * ms / size), max(ms, _FLOOR * fs * size)
        scale = np.array([max(fs, 1e-300)] * 3 + [max(ms, 1e-300)] * 3)
        worst = {"value": 0.0}

        def note(value, where, quantity):
            if value > worst["value"]:
                worst.update(value=float(value), where=where, quantity=quantity, tolerance=tol)

        for mid, mr in cres.members.items():
            pm = p.members[mid]
            x = np.asarray(mr.x, float)
            F = np.asarray(mr.forces, float)
            n = len(x)
            n_stations += n
            if not mr.active:
                v = float(np.abs(F).max(initial=0.0) / scale.min()) if F.size else 0.0
                v = max(v, float(np.abs(np.asarray(mr.end_forces, float)).max(initial=0.0) / scale.min()))
                note(v, f"inactive member {mid}", "station forces")
                if v > tol:
                    failed.append(f"{cid}: inactive member {mid} carries forces")
                continue
            # ---- the ends of the stations are the end forces (CONVENTIONS 4)
            ef = np.asarray(mr.end_forces, float)
            d0 = float((np.abs(F[0] + ef[:6]) / scale).max())
            dL = float((np.abs(F[-1] - ef[6:]) / scale).max())
            note(max(d0, dL), f"member {mid} station forces at the ends against the end forces", "end")
            if max(d0, dL) > tol:
                failed.append(f"{cid}: member {mid}: station forces at the ends differ from the end forces by "
                              f"{max(d0, dL):.2e}")
            L = float(x[-1] - x[0]) if n > 1 else 0.0
            eps = _EPS_X * max(L, 1e-300)
            loads_scale = 1.0 / alpha                    # stored results are divided by alpha
            pts = [(kind, np.asarray(vec, float) * loads_scale, float(xs) + float(x[0])) for kind, vec, xs in pm.points]
            dsts = [(np.asarray(w1, float) * loads_scale, np.asarray(w2, float) * loads_scale,
                     float(a) + float(x[0]), float(b) + float(x[0])) for w1, w2, a, b in pm.dists]
            U = np.asarray(mr.displacements, float)
            moments_ok = True
            if second and any(abs(w1[0]) + abs(w2[0]) > 0.0 for w1, w2, _, _ in dsts):
                bound = sum(0.5 * (abs(w1[0]) + abs(w2[0])) * (b - a) for w1, w2, a, b in dsts) \
                    * alpha * 2.0 * float(np.abs(U[:, 1:3]).max(initial=0.0))
                if bound / scale[3] > 0.1 * tol:
                    moments_ok = False
                    deferred.append(f"{cid}: member {mid} (axial line load in a second-order combination; lever "
                                    f"bound {bound / scale[3]:.1e} of the moment scale)")

            def deflection_at(s):
                k = int(np.argmin(np.abs(x - s)))
                return U[k, 1:3] if abs(x[k] - s) <= eps else None

            sides = {}                                     # station -> "left" | "right" | "both" at a point load
            for k in range(1, n):
                xk = float(x[k])
                Fv = F[0, :3].copy()
                # (r_0 - r_x) x F(0)
                Mv = F[0, 3:] + np.cross([float(x[0]) - xk, 0.0, 0.0], F[0, :3])
                at = np.zeros(6)                               # the point loads AT this station (within eps)
                lever_missing = False
                if second:
                    d = alpha * (U[0, 1:3] - U[k, 1:3])              # (dv, dw) of r_0 - r_x
                    Mv += np.array([0.0, d[1] * F[0, 0], -d[0] * F[0, 0]])
                for kind, vec, xs in pts:
                    if abs(xs - xk) <= eps:
                        at -= np.r_[vec, 0.0, 0.0, 0.0] if kind == "force" else np.r_[0.0, 0.0, 0.0, vec]
                    elif xs < xk:
                        if kind == "force":
                            Fv -= vec
                            Mv -= np.cross([xs - xk, 0.0, 0.0], vec)
                            if second and vec[0] != 0.0:
                                ui = deflection_at(xs)
                                if ui is None:
                                    lever_missing = True
                                else:
                                    d = alpha * (ui - U[k, 1:3])
                                    Mv -= np.array([0.0, d[1] * vec[0], -d[0] * vec[0]])
                        else:
                            Mv -= vec
                for w1, w2, a, b in dsts:
                    force, first, T = _line_load(w1, w2, a, b, xk)
                    if T > 0.0:
                        Fv -= force
                        Mv -= np.cross([1.0, 0.0, 0.0], first)
                        if second and (w1[0] != 0.0 or w2[0] != 0.0) and moments_ok:
                            # lever of the axial line load: trapezoid over the stations (negligible, see above)
                            idx = [j for j in range(n) if a - eps <= x[j] <= min(xk, b) + eps]
                            if len(idx) > 1:
                                wx = np.array([w1[0] + (w2[0] - w1[0]) * (x[j] - a) / (b - a) for j in idx])
                                dd = alpha * (U[idx, 1:3] - U[k, 1:3])
                                iv = np.trapezoid(dd[:, 0] * wx, x[idx])
                                iw = np.trapezoid(dd[:, 1] * wx, x[idx])
                                Mv -= np.array([0.0, iw, -iv])
                skip_m = not moments_ok or lever_missing
                if lever_missing and moments_ok:
                    moments_ok = False
                    not_evaluated.append(f"{cid}: member {mid}: bending moments not evaluated (no station at an "
                                         "axial point load)")

                def difference(extra):
                    dd = np.abs(np.r_[Fv, Mv] + extra - F[k]) / scale
                    if skip_m:
                        dd[4:] = 0.0
                    return dd
                if not at.any():
                    diff = difference(0.0)
                elif k == n - 1:
                    diff = difference(at)                      # the last station is the limit from outside
                else:
                    # a station AT a point load (within the round-off of its position) holds the limit from the
                    # left or from the right; which of the two is decided by an ulp of x, so either is accepted
                    # here and the stations of the point together must show both (checked below)
                    dl, dr = difference(0.0), difference(at)
                    ok_l, ok_r = bool(dl.max() <= tol), bool(dr.max() <= tol)
                    sides[k] = "both" if ok_l and ok_r else ("left" if ok_l else ("right" if ok_r else "none"))
                    diff = dl if dl.max() <= dr.max() else dr
                j = int(np.argmax(diff))
                note(float(diff[j]), f"member {mid} at x = {xk:.6g} mm", ("N", "Vy", "Vz", "T", "My", "Mz")[j])
                if diff[j] > tol:
                    failed.append(f"{cid}: member {mid} station {k} (x = {xk:.6g} mm): "
                                  f"{('N', 'Vy', 'Vz', 'T', 'My', 'Mz')[j]} differs from statics by {diff[j]:.2e} "
                                  f"(tolerance {tol:.0e})")
            # every interior point load must be seen from both sides (docs/THEORY.md 5.3: two stations)
            for xs in sorted({xs for _, _, xs in pts if eps < xs - float(x[0]) and xs < float(x[-1]) - eps}):
                here = [sides[k] for k in range(1, n - 1) if abs(float(x[k]) - xs) <= eps and k in sides]
                if "none" in here:
                    continue                                   # reported above
                if not (any(v in ("left", "both") for v in here) and any(v in ("right", "both") for v in here)):
                    failed.append(f"{cid}: member {mid}: the stations at the point load at x = {xs:.6g} mm do not "
                                  "show both the limit from the left and the limit from the right")
        per_combo[cid] = {**worst, "tolerance": tol, "secondOrder": second}
        if worst["value"] > worst_all["value"]:
            worst_all = {**worst, "combination": cid}
    return {"passed": not failed and not not_evaluated, "failed": failed, "notEvaluated": not_evaluated,
            "deferred": deferred, "largest": worst_all, "combinations": per_combo, "stations": n_stations}


# ====================================================================================== deflections
def _plane_solution(x, brk, EI, beta, N_of, q_of, jumps, y0, second, alpha, atol):
    """Integrate one bending plane from the first station (docs/THEORY.md 1.1 and 10.2, eq. (5)):

        v' = (th + beta V) / (1 + beta N)    th' = M / EI    M' = (N th - V) / (1 + beta N)    V' = -q

    (first order: the same with N = 0) as an initial value problem with ``scipy.integrate.solve_ivp``.
    The start rotation is not part of the results: the particular solution (th0 = 0) and the
    homogeneous one (th0 = 1) are integrated together and returned at the stations ``x``.
    ``jumps``: [(position, dM, dV)] of point loads; ``brk``: positions where the loading changes.
    Returns ``(left, right)``: the state just before and just after every station position (they differ
    only at a point load), or None if the integration failed."""
    from scipy.integrate import solve_ivp

    def rhs(s, y):
        N = alpha * N_of(s) if second else 0.0
        den = 1.0 + beta * N
        q = q_of(s)
        out = np.empty(8)
        for o, qq in ((0, q), (4, 0.0)):
            v, th, M, V = y[o:o + 4]
            out[o] = (th + beta * V) / den
            out[o + 1] = M / EI
            out[o + 2] = (N * th - V) / den
            out[o + 3] = -qq
        return out

    y = np.r_[y0, 0.0, 1.0, 0.0, 0.0]
    cuts = sorted({float(x[0]), float(x[-1])} | {float(b) for b in brk if x[0] < b < x[-1]})
    L = float(x[-1] - x[0])
    eps = _EPS_X * max(L, 1e-300)
    n = len(x)
    left, right = np.zeros((n, 8)), np.zeros((n, 8))        # state just before / just after position x[k]
    done = np.zeros(n, bool)
    for a, b in zip(cuts[:-1], cuts[1:]):
        before = y
        for pos, dM, dV in jumps:                              # point loads at the start of this piece
            if abs(pos - a) <= eps:
                y = y + np.r_[0.0, 0.0, dM, dV, 0.0, 0.0, 0.0, 0.0]
        for k in range(n):                                     # stations AT the start of this piece
            if abs(x[k] - a) <= eps:
                left[k], right[k], done[k] = before, y, True
        idx = [k for k in range(n) if a + eps < x[k] < b - eps]
        te = sorted({float(x[k]) for k in idx} | {b})
        r = solve_ivp(rhs, (a, b), y, method="DOP853", t_eval=te, rtol=1e-12, atol=atol)
        if not r.success:
            return None
        for k in idx:
            j = int(np.argmin(np.abs(np.asarray(te) - float(x[k]))))
            left[k] = right[k] = r.y[:, j]
            done[k] = True
        y = r.y[:, -1]
    before = y
    for pos, dM, dV in jumps:                                  # a point load at the very end
        if abs(pos - cuts[-1]) <= eps:
            y = y + np.r_[0.0, 0.0, dM, dV, 0.0, 0.0, 0.0, 0.0]
    for k in range(n):
        if abs(x[k] - cuts[-1]) <= eps:
            left[k], right[k], done[k] = before, y, True
    left[0] = right[0] = np.r_[y0, 0.0, 1.0, 0.0, 0.0]         # the first station is the limit from outside
    left[-1] = right[-1]                                       # ... and so is the last one
    return (left, right) if done.all() else None


def station_deflections(model, results, tol: float) -> dict:
    """The displacements of the member axis and (again) the bending moments and shears at the result
    stations, recomputed by integrating the member differential equations from the first station with
    the stored start values (no solver code; initial value problems, ``scipy.integrate.solve_ivp``).

    Per member and bending plane: start deflection, moment and shear from the first station; the start
    rotation, which the results do not hold, from the deflection of the LAST station; compared at
    every interior station: deflection, moment, shear.  Axial displacement: ``u' = N / EA`` with N from
    statics.  The displacement of the member axis at each end must equal that of the node (carried over
    the rigid offset), unless a translational connection of that end is released.  Differences
    relative to the largest member-axis translation / station moment / station force of the
    combination.  Returns the same structure as :func:`station_statics`."""
    xyz = np.array([nd.xyz for nd in model.nodes.values()], float)
    size = float(np.linalg.norm(xyz.max(axis=0) - xyz.min(axis=0))) if len(xyz) else 1.0
    size = size if size > 0.0 else 1.0
    failed, not_evaluated, per_combo = [], [], {}
    worst_all = {"value": 0.0}
    n_members = 0
    for cid, cres in results.combinations.items():
        combo = model.combinations[cid]
        second = combo.analysis_type == "secondOrder"
        alpha = float(getattr(cres, "alpha", 1.0) or 1.0)
        p = prepare(model, cid, tau_b=dict(getattr(cres, "tau_b", {}) or {}), extra_nodal_loads={})
        fs = max((float(np.abs(mr.forces[:, :3]).max(initial=0.0)) for mr in cres.members.values()), default=0.0)
        ms = max((float(np.abs(mr.forces[:, 3:]).max(initial=0.0)) for mr in cres.members.values()), default=0.0)
        fs, ms = max(fs, _FLOOR * ms / size, 1e-300), max(ms, _FLOOR * fs * size, 1e-300)
        ds = max(max((float(np.abs(mr.displacements).max(initial=0.0)) for mr in cres.members.values()), default=0.0),
                 1e-300)
        worst = {"value": 0.0}

        def note(value, where, quantity):
            if value > worst["value"]:
                worst.update(value=float(value), where=where, quantity=quantity, tolerance=tol)

        for mid, mr in cres.members.items():
            if not mr.active:
                continue
            pm = p.members[mid]
            x = np.asarray(mr.x, float)
            if len(x) < 3:
                continue
            n_members += 1
            F, U = np.asarray(mr.forces, float), np.asarray(mr.displacements, float)
            x0 = float(x[0])
            pts = [(kind, np.asarray(vec, float) / alpha, float(xs) + x0) for kind, vec, xs in pm.points]
            dsts = [(np.asarray(w1, float) / alpha, np.asarray(w2, float) / alpha, float(a) + x0, float(b) + x0)
                    for w1, w2, a, b in pm.dists]
            brk = [xs for _, _, xs in pts] + [a for *_, a, _ in dsts] + [b for *_, b in dsts]

            def line(comp):
                def q(s):
                    return sum(float(w1[comp] + (w2[comp] - w1[comp]) * (s - a) / (b - a))
                               for w1, w2, a, b in dsts if a <= s <= b)
                return q

            def axial(s):                                      # N(s) from statics, right limit
                N = F[0, 0]
                for kind, vec, xs in pts:
                    if kind == "force" and xs <= s:
                        N -= vec[0]
                for w1, w2, a, b in dsts:
                    N -= _line_load(w1, w2, a, b, s)[0][0]
                return float(N)

            # ---- the ends of the member axis move with the nodes (rigid offset kinematics), unless a
            #      translational connection of that end is released or a spring
            mm = model.members[mid]
            Rm = np.asarray(mm.R, float)
            for end, node, off, kk, row in (("i", mm.i, mm.offset_i, mm.k_end_i, 0), ("j", mm.j, mm.offset_j, mm.k_end_j, -1)):
                if any(k != np.inf for k in kk[:3]):
                    continue
                un = np.asarray(cres.displacements[model.nodes[node].index], float)
                expected = Rm @ (un[:3] + np.cross(un[3:], np.asarray(off, float)))
                d = float(np.abs(expected - U[row]).max()) / ds
                note(d, f"member {mid} end {end} against node {node}", "displacement of the member end")
                if d > tol:
                    failed.append(f"{cid}: member {mid}: the displacement of the member axis at end {end} differs "
                                  f"from the displacement of node {node} by {d:.2e} (tolerance {tol:.0e})")
            bad = False
            # plane xy: (v, theta_z, Mz, Vy), q = w_y;  plane xz: (w, -theta_y, -My, Vz), q = w_z
            for plane, (ic, im, sgn, EI, GA, qc) in {"xy": (1, 5, 1.0, pm.EIz, pm.GAy, 1),
                                                     "xz": (2, 4, -1.0, pm.EIy, pm.GAz, 2)}.items():
                beta = 0.0 if not GA else 1.0 / float(GA)
                jumps = []
                for kind, vec, xs in pts:
                    if kind == "force" and vec[qc] != 0.0:
                        jumps.append((xs, 0.0, -float(vec[qc])))
                    elif kind == "moment" and vec[2 if plane == "xy" else 1] != 0.0:
                        jumps.append((xs, -sgn * float(vec[2 if plane == "xy" else 1]), 0.0))
                y0 = np.array([U[0, ic], 0.0, sgn * F[0, im], F[0, ic]])
                if second and beta and (1.0 + beta * alpha * min(axial(float(s)) for s in x)) <= 0.0:
                    not_evaluated.append(f"{cid}: member {mid}: compression beyond the shear stiffness")
                    bad = True
                    break
                Lm = float(x[-1] - x0)
                atol = 1e-13 * np.array([ds, ds / Lm, ms, fs, Lm, 1.0, max(fs * Lm, float(EI) / Lm), fs])
                sol = _plane_solution(x, brk, float(EI), beta, axial, line(qc), jumps, y0, second, alpha, atol)
                if sol is None or abs(sol[1][-1, 4]) < 1e-12 * max(float(x[-1] - x0), 1e-300):
                    not_evaluated.append(f"{cid}: member {mid}, plane {plane}: the integration of the member "
                                         "equations did not succeed")
                    bad = True
                    break
                th0 = (U[-1, ic] - sol[1][-1, 0]) / sol[1][-1, 4]
                # at a point load a station holds the limit from the left or from the right (decided by an
                # ulp of its position): the nearer of the two integrated states is compared; that the
                # stations of the point show both limits is checked by ``station_statics``
                both = [side[:, :4] + th0 * side[:, 4:] for side in sol]
                pick = []
                for k in range(len(x)):
                    cand = [max(abs(st[k, 0] - U[k, ic]) / ds, abs(sgn * st[k, 2] - F[k, im]) / ms,
                                abs(st[k, 3] - F[k, ic]) / fs) for st in both]
                    pick.append(both[0][k] if cand[0] <= cand[1] else both[1][k])
                state = np.array(pick)
                dv = np.abs(state[1:-1, 0] - U[1:-1, ic]) / ds
                dm = np.abs(sgn * state[1:, 2] - F[1:, im]) / ms
                dq = np.abs(state[1:, 3] - F[1:, ic]) / fs
                for arr, off, what in ((dv, 1, "deflection " + ("uy" if plane == "xy" else "uz")),
                                       (dm, 1, "Mz" if plane == "xy" else "My"), (dq, 1, "Vy" if plane == "xy" else "Vz")):
                    if arr.size:
                        k = int(np.argmax(arr))
                        note(float(arr[k]), f"member {mid} at x = {x[k + off]:.6g} mm", what)
                        if arr[k] > tol:
                            failed.append(f"{cid}: member {mid} at x = {x[k + off]:.6g} mm: {what} differs from the "
                                          f"integrated member equations by {arr[k]:.2e} (tolerance {tol:.0e})")
            if bad:
                continue
            # axial displacement: u' = N / EA (first order in both analyses)
            EA = float(pm.E * pm.A)
            u = U[0, 0]
            eps = _EPS_X * max(float(x[-1] - x0), 1e-300)
            gp, gw = np.polynomial.legendre.leggauss(3)        # exact: N is at most quadratic between loads
            for k in range(1, len(x)):
                a, b = float(x[k - 1]), float(x[k])
                if b - a > eps:
                    s = 0.5 * (a + b) + 0.5 * (b - a) * gp
                    u += 0.5 * (b - a) * float(sum(w * axial(si) for w, si in zip(gw, s))) / EA
                d = abs(u - U[k, 0]) / ds
                note(d, f"member {mid} at x = {b:.6g} mm", "axial displacement ux")
                if d > tol:
                    failed.append(f"{cid}: member {mid} at x = {b:.6g} mm: axial displacement differs from the "
                                  f"integrated N / EA by {d:.2e} (tolerance {tol:.0e})")
        per_combo[cid] = {**worst, "tolerance": tol, "secondOrder": second}
        if worst["value"] > worst_all["value"]:
            worst_all = {**worst, "combination": cid}
    return {"passed": not failed and not not_evaluated, "failed": failed, "notEvaluated": not_evaluated,
            "largest": worst_all, "combinations": per_combo, "members": n_members}
