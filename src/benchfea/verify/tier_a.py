"""Tier A (automatic part): the engine against the independent ODE oracle (``benchfea.oracle``).

Applies to models the oracle can represent exactly: ONE STRAIGHT LINE of frame members (a single
member, a continuous beam, a beam-column), with

* supports and support springs at the nodes, imposed support displacements,
* every load of the format (nodal forces and moments, member point forces and moments, linearly varying
  line loads, self weight), in both bending planes, axially and in torsion,
* bending releases (full or spring) that are an internal hinge of the line or a released end at a
  rotation-fixed end node,
* first-order analysis with or without shear deformation, and second-order analysis when the axial
  force is piecewise constant (no axial line load) -- without shear deformation, because the oracle
  refuses to choose between the Engesser and Haringx shear forces (docs/CROSSCHECK.md 6.2),
* stiffness reduction with a given tau_b.

The line is turned into two ``oracle.Beam`` problems (local x-y and x-z planes) and two ``oracle.Bar``
problems (axial, torsion) from the load definitions (``crosscheck.prepare``: no solver code), solved by
``scipy.integrate.solve_bvp`` and compared with the engine at every node and result station.  The
oracle shares no stiffness matrix, shape function or fixed-end force with the engine.  Targets produced
this way are "derived here"; they are not attributed to any book.

Everything else (frames, constraints, offsets, trusses, notional loads, computed tau_b, buckling) is
outside this tier: ``applies = False`` with the reason.
"""
from __future__ import annotations

import time

import numpy as np

from ..crosscheck.prepare import prepare
from ..model import RIGID
from ..oracle import Bar, Beam, OracleError
from .policy import POLICY
from .record import TierResult

__all__ = ["tier_a", "line_representation", "NotRepresentable", "solve_line"]

FIXED = "fixed"


class NotRepresentable(Exception):
    """The oracle cannot represent the model; the message is the reason."""


# ====================================================================================== representation
def line_representation(model) -> dict:
    """Geometry of the member line or raise :class:`NotRepresentable`.

    Returns ``{"R", "sign", "gdof", "nodes": [ids along the line], "s": {node: coordinate},
    "members": [ids along the line], "L"}``; ``gdof[a]`` / ``sign[a]`` map local DOF a (0..5) to the
    global DOF and its sign."""
    if not model.members:
        raise NotRepresentable("no members")
    if model.constraints:
        raise NotRepresentable("rigid links / diaphragms")
    R0 = None
    for mid, m in model.members.items():
        if m.type != "frame":
            raise NotRepresentable(f"member {mid!r} is a {m.type} member")
        if m.has_offsets:
            raise NotRepresentable(f"member {mid!r} has rigid offsets")
        R = np.asarray(m.R, float)
        if R0 is None:
            R0 = R
        elif not np.allclose(R, R0, atol=1e-12):
            raise NotRepresentable("the members are not collinear with identical local axes")
        for end, k in (("i", m.k_end_i), ("j", m.k_end_j)):
            if any(k[a] != RIGID for a in (0, 1, 2, 3)):
                raise NotRepresentable(f"member {mid!r}: axial, shear or torsion release at end {end}")
    if not np.allclose(np.sort(np.abs(R0), axis=1)[:, :2], 0.0, atol=1e-12):
        raise NotRepresentable("the member axes are not parallel to the global axes (supports could not be "
                               "mapped one-to-one onto the bending planes)")
    gdof = [int(np.argmax(np.abs(R0[a]))) for a in range(3)]
    sign = [float(np.sign(R0[a, gdof[a]])) for a in range(3)]
    gdof, sign = gdof + [g + 3 for g in gdof], sign + sign
    ex = R0[0]
    used = {n for m in model.members.values() for n in (m.i, m.j)}
    if set(model.nodes) != used:
        raise NotRepresentable("a node without members")
    x0 = min((np.asarray(model.nodes[n].xyz, float) for n in used), key=lambda p: float(p @ ex))
    s = {}
    for n in used:
        d = np.asarray(model.nodes[n].xyz, float) - x0
        if np.linalg.norm(d - (d @ ex) * ex) > 1e-9 * max(1.0, np.linalg.norm(d)):
            raise NotRepresentable("the nodes are not on one straight line")
        s[n] = float(d @ ex)
    nodes = sorted(used, key=lambda n: s[n])
    by_start = {}
    for mid, m in model.members.items():
        if m.i in by_start:
            raise NotRepresentable("two members start at the same node")
        by_start[m.i] = mid
    members = []
    for a, b in zip(nodes[:-1], nodes[1:]):
        mid = by_start.get(a)
        if mid is None or model.members[mid].j != b:
            raise NotRepresentable("the members do not form one chain along the line")
        members.append(mid)
    if len(members) != len(model.members):
        raise NotRepresentable("the members do not form one chain along the line")
    return {"R": R0, "sign": sign, "gdof": gdof, "nodes": nodes, "s": s, "members": members,
            "L": s[nodes[-1]] - s[nodes[0]]}


def _combination_problem(model, line, cid):
    """Everything the oracle needs for one combination, or raise NotRepresentable."""
    combo = model.combinations[cid]
    if combo.notional is not None:
        raise NotRepresentable("notional loads")
    if combo.stiffness_reduction and model.analysis.tau_b == "computed":
        raise NotRepresentable("computed tau_b (it depends on the solution being checked)")
    p = prepare(model, cid)
    second = p.second_order
    nodes, s, sign, gdof = line["nodes"], line["s"], line["sign"], line["gdof"]
    idx = {n: model.nodes[n].index for n in nodes}
    # ---- hinges and released ends, per plane (local rotation DOF 5 = plane xy, 4 = plane xz)
    hinges = {4: {}, 5: {}}            # x -> stiffness
    end_release = {4: {}, 5: {}}       # end node -> stiffness (member end released at a rotation-fixed end node)
    follow = {4: {}, 5: {}}            # node -> side (+1/-1) of the member the node rotation belongs to
    for a in (4, 5):
        g = gdof[a]
        for k, n in enumerate(nodes):
            left = model.members[line["members"][k - 1]].k_end_j[a] if k > 0 else None
            right = model.members[line["members"][k]].k_end_i[a] if k < len(nodes) - 1 else None
            rel = [v for v in (left, right) if v is not None and v != RIGID]
            if not rel:
                continue
            sup = model.supports.get(n)
            rot_fixed = bool(sup is not None and sup.fixed[g])
            rot_spring = float(sup.spring[g]) if sup is not None else 0.0
            moment = any(abs(pp.load[idx[n], g]) > 0.0 for pp in (p,)) or _any_nodal_moment(model, n, g)
            if k in (0, len(nodes) - 1):
                if not rot_fixed or moment:
                    raise NotRepresentable(f"released member end at end node {n!r} whose rotation is not fixed "
                                           "(or that carries a nodal moment)")
                end_release[a][n] = float(rel[0])
            else:
                if len(rel) == 2:
                    raise NotRepresentable(f"both member ends released at node {n!r}")
                if rot_fixed or rot_spring or moment:
                    raise NotRepresentable(f"hinge at node {n!r} together with a rotational support or a nodal "
                                           "moment (the side of the hinge matters)")
                hinges[a][s[n]] = float(rel[0])
                follow[a][n] = +1 if left is not None and left != RIGID else -1
    # ---- loads
    axial_line_load = any(abs(w[0]) > 0.0 for m in p.members.values() for w1, w2, _, _ in m.dists for w in (w1, w2))
    if second:
        if axial_line_load:
            raise NotRepresentable("second-order analysis with an axial line load (the axial force is not "
                                   "piecewise constant)")
    return {"p": p, "second": second, "hinges": hinges, "end_release": end_release, "follow": follow, "idx": idx}


def _any_nodal_moment(model, node, g) -> bool:
    from ..model import NodalLoad
    return any(isinstance(ld, NodalLoad) and ld.node == node and abs(ld.values[g]) > 0.0
               for case in model.load_cases.values() for ld in case.loads)


# ====================================================================================== oracle solution
def solve_line(model, line, cid, tol=1e-10):
    """Solve one combination with the oracle.  Returns a dict of callables/arrays in LOCAL axes with
    results already divided by alpha: ``u, v, w, tx, ty, tz`` (x, side), ``N, Vy, Vz, T, My, Mz``
    (x, side) and ``reactions[node] = (6,) local``."""
    prob = _combination_problem(model, line, cid)
    p, second = prob["p"], prob["second"]
    nodes, s, sign, gdof, L = line["nodes"], line["s"], line["sign"], line["gdof"], line["L"]
    idx = prob["idx"]
    spans = [(s[model.members[mid].i], s[model.members[mid].j], p.members[mid]) for mid in line["members"]]
    pw = lambda f: [(a, b, f(m)) for a, b, m in spans]

    def local(arr, n, a):                       # local component a of a nodal array (n_nodes, 6) global
        return sign[a] * float(arr[idx[n], gdof[a]])

    def add_loads(line_obj, comp, moment_comp=None, node_dof=None, node_moment_dof=None):
        for n in nodes:
            f = local(p.load, n, node_dof)
            mo = local(p.load, n, node_moment_dof) if node_moment_dof is not None else 0.0
            if moment_comp is None:
                if f:
                    line_obj.point(s[n], f)
            elif f or mo:
                line_obj.point(s[n], P=f, M=mo)
        for a, b, m in spans:
            for kind, vec, x in m.points:
                if kind == "force" and vec[comp]:
                    line_obj.point(a + x, vec[comp]) if moment_comp is None else line_obj.point(a + x, P=vec[comp])
                if kind == "moment" and moment_comp is not None and vec[moment_comp]:
                    line_obj.point(a + x, M=vec[moment_comp])
            for w1, w2, x1, x2 in m.dists:
                if w1[comp] or w2[comp]:
                    line_obj.dist(a + x1, a + x2, w1[comp], w2[comp])

    def support_value(n, a):
        sup = model.supports.get(n)
        if sup is None:
            return None
        g = gdof[a]
        if sup.fixed[g]:
            return FIXED
        return float(sup.spring[g]) if sup.spring[g] else None

    # ---- axial and torsion
    def bar(a, stiff, torque):
        b = Bar(L, pw(stiff))
        any_support = False
        for n in nodes:
            v = support_value(n, a)
            if v is not None:
                b.support(s[n], u=v, settlement=local(p.prescribed, n, a))
                any_support = True
            f = local(p.load, n, a)
            if f:
                b.point(s[n], f)
        loaded = any(local(p.load, n, a) for n in nodes)
        for x0, _, m in spans:
            for kind, vec, x in m.points:
                if kind == ("moment" if torque else "force") and vec[0]:
                    b.point(x0 + x, vec[0])
                    loaded = True
            if not torque:
                for w1, w2, x1, x2 in m.dists:
                    if w1[0] or w2[0]:
                        b.dist(x0 + x1, x0 + x2, w1[0], w2[0])
                        loaded = True
        if not any_support:
            if loaded:
                raise NotRepresentable("axial / torsional load without a support in that direction")
            return None
        return b.solve(tol)

    axial = bar(0, lambda m: m.E * m.A, False)
    torsion = bar(3, lambda m: m.G * m.J, True)
    # ---- axial force for the second-order bending problems
    n_spec = 0.0
    if second and axial is not None:
        br = axial.breaks
        n_spec = [(float(a), float(b), axial.N(0.5 * (a + b), -1)) for a, b in zip(br[:-1], br[1:])]
        if not any(abs(v) > 0.0 for *_, v in n_spec):
            n_spec = 0.0
    shear = bool(model.analysis.shear_deformation)
    planes = {}
    for plane, a_t, a_r, comp, mcomp, ei, ga in (
            ("xy", 1, 5, 1, 2, lambda m: m.EIz, lambda m: m.GAy),
            ("xz", 2, 4, 2, 1, lambda m: m.EIy, lambda m: m.GAz)):
        gas = None
        if shear and any(ga(m) for _, _, m in spans):
            if not all(ga(m) for _, _, m in spans):
                raise NotRepresentable("shear deformation on some members of the line only")
            gas = pw(ga)
            if second and n_spec != 0.0:
                raise NotRepresentable("second-order analysis with shear deformation and axial force: the shear "
                                       "force definition (Engesser / Haringx) is an open question")
        bm = Beam(L, pw(ei), GAs=gas, N=n_spec, second_order=bool(second), plane=plane)
        for n in nodes:
            v, r = support_value(n, a_t), support_value(n, a_r)
            if n in prob["end_release"][a_r]:
                k = prob["end_release"][a_r][n]
                r = None if k == 0.0 else k
            if v is not None or r is not None:
                bm.support(s[n], v=v, r=r, settlement=local(p.prescribed, n, a_t),
                           rotation=local(p.prescribed, n, a_r))
        for x, k in prob["hinges"][a_r].items():
            bm.hinge(x, k)
        add_loads(bm, comp, mcomp, a_t, a_r)
        planes[plane] = bm.solve(tol)
    al = p.alpha
    xy, xz = planes["xy"], planes["xz"]
    zero = lambda x, side=-1: 0.0
    out = {
        "u": (lambda x, side=-1: axial.u(x, side) / al) if axial else zero,
        "v": lambda x, side=-1: xy.v(x, side) / al, "w": lambda x, side=-1: xz.v(x, side) / al,
        "tx": (lambda x, side=-1: torsion.u(x, side) / al) if torsion else zero,
        "ty": lambda x, side=-1: xz.theta(x, side) / al, "tz": lambda x, side=-1: xy.theta(x, side) / al,
        "N": (lambda x, side=-1: axial.N(x, side) / al) if axial else zero,
        "Vy": lambda x, side=-1: xy.Q(x, side) / al, "Vz": lambda x, side=-1: xz.Q(x, side) / al,
        "T": (lambda x, side=-1: torsion.N(x, side) / al) if torsion else zero,
        "My": lambda x, side=-1: xz.M(x, side) / al, "Mz": lambda x, side=-1: xy.M(x, side) / al,
        "residuals": {"xy": xy.rms_residual, "xz": xz.rms_residual,
                      "axial": axial.rms_residual if axial else None,
                      "torsion": torsion.rms_residual if torsion else None},
    }
    reactions = {}
    for n in nodes:
        if n not in model.supports:
            continue
        x = s[n]
        r = np.zeros(6)

        def get(sol, key=x):
            for kx, val in sol.reactions.items():
                if abs(kx - key) <= 1e-9 * L:
                    return val
            return None
        if axial and get(axial) is not None:
            r[0] = get(axial)
        if torsion and get(torsion) is not None:
            r[3] = get(torsion)
        if get(xy) is not None:
            r[1], r[5] = get(xy)
        if get(xz) is not None:
            r[2], r[4] = get(xz)
        reactions[n] = r / al
    out["reactions"] = reactions
    out["problem"] = prob
    return out


# ====================================================================================== comparison
def _sides(xs, L):
    """Side (-1 left limit, +1 right limit) of every station of a member (a repeated x is left, right)."""
    out = []
    for k, x in enumerate(xs):
        if k > 0 and xs[k - 1] == x:
            out.append(+1)
        elif k == 0:
            out.append(+1)
        else:
            out.append(-1)
    return out


def _compare(model, line, cid, cres, sol) -> dict:
    s, nodes, sign, gdof, L = line["s"], line["nodes"], line["sign"], line["gdof"], line["L"]
    prob = sol["problem"]
    rows = []                       # (kind, where, engine, oracle)
    names_d = ("u", "v", "w", "tx", "ty", "tz")
    for k, n in enumerate(nodes):
        x = s[n]
        side = -1 if k == len(nodes) - 1 else (+1 if k == 0 else -1)
        for a, name in enumerate(names_d):
            sd = side
            if a in (4, 5):
                if n in prob["end_release"][a]:
                    continue                                # node rotation = imposed value, not the member's
                sd = prob["follow"][a].get(n, side)
            eng = sign[a] * float(cres.displacements[model.nodes[n].index, gdof[a]])
            rows.append(("translation" if a < 3 else "rotation", f"node {n} local {name}", eng, sol[name](x, sd)))
        if n in sol["reactions"]:
            for a in range(6):
                if a in (4, 5) and n in prob["end_release"][a]:
                    continue
                eng = sign[a] * float(cres.reactions[model.nodes[n].index, gdof[a]])
                rows.append(("force" if a < 3 else "moment", f"reaction {n} local {'Fx Fy Fz Mx My Mz'.split()[a]}",
                             eng, float(sol["reactions"][n][a])))
    names_f = ("N", "Vy", "Vz", "T", "My", "Mz")
    for mid in line["members"]:
        m = model.members[mid]
        mr = cres.members[mid]
        x0 = s[m.i]
        Lm = float(m.length)
        n_st = len(mr.x)
        for k, (x, side) in enumerate(zip(mr.x, _sides(list(mr.x), Lm))):
            xx = min(max(x0 + float(x), 0.0), L)
            if abs(float(x) - Lm) <= 1e-12 * Lm:
                side = -1
            # An interior station AT a point load holds the limit from the left or from the right; which
            # one is decided by an ulp of its position (a regular station that falls on a point load gives
            # three stations there: left, left, right or left, right, right).  The oracle value compared is
            # therefore the nearer of its two limits; that the stations of the point show BOTH limits is
            # checked by tier D (benchfea.verify.stations).
            interior = 0 < k < n_st - 1

            def oracle(name, engine):
                o = sol[name](xx, side)
                if interior:
                    other = sol[name](xx, -side)
                    if other != o and abs(engine - other) < abs(engine - o):
                        return other
                return o
            for a, name in enumerate(names_f):
                e = float(mr.forces[k, a])
                rows.append(("force" if a < 3 else "moment", f"member {mid} station {k} (x = {float(x):.6g} mm) {name}",
                             e, oracle(name, e)))
            for a, name in enumerate(("u", "v", "w")):
                e = float(mr.displacements[k, a])
                rows.append(("translation", f"member {mid} station {k} (x = {float(x):.6g} mm) u{'xyz'[a]}",
                             e, oracle(name, e)))
    scale = {}
    for kind in ("translation", "rotation", "force", "moment"):
        vals = [abs(v) for kd, _, e, o in rows if kd == kind for v in (e, o)]
        scale[kind] = max(vals, default=0.0)
    scale["rotation"] = max(scale["rotation"], scale["translation"] / L)
    scale["moment"] = max(scale["moment"], 1e-3 * scale["force"] * L)
    worst = {}
    for kind, where, e, o in rows:
        sc = scale[kind]
        if sc == 0.0:
            continue
        d = abs(e - o) / sc
        if kind not in worst or d > worst[kind]["value"]:
            worst[kind] = {"quantity": kind, "where": where, "engine": e, "oracle": o, "value": d}
    return {"worst": worst, "scales": scale, "compared": len(rows), "residuals": sol["residuals"]}


def tier_a(model, results) -> TierResult:
    t0 = time.perf_counter()
    tol = POLICY["A"]["relative"]
    try:
        line = line_representation(model)
    except NotRepresentable as e:
        return TierResult("A", False, None, f"not a model the ODE oracle can represent: {e}", [], {},
                          {"reason": str(e)}, time.perf_counter() - t0)
    per_combo, failed, largest, not_rep = {}, [], [], {}
    for cid, cres in results.combinations.items():
        try:
            sol = solve_line(model, line, cid)
        except NotRepresentable as e:
            not_rep[cid] = str(e)
            continue
        except OracleError as e:
            per_combo[cid] = {"passed": False, "error": f"the oracle could not solve the problem: {e}"}
            failed.append(f"{cid}: oracle failure ({e})")
            continue
        cmp = _compare(model, line, cid, cres, sol)
        ok = all(w["value"] <= tol for w in cmp["worst"].values())
        per_combo[cid] = {"passed": ok, **cmp}
        for w in cmp["worst"].values():
            largest.append({"what": w["quantity"], "combination": cid, "value": w["value"], "tolerance": tol,
                            "where": w["where"], "engine": w["engine"], "oracle": w["oracle"]})
        if not ok:
            w = max(cmp["worst"].values(), key=lambda v: v["value"])
            failed.append(f"{cid}: {w['quantity']} differs by {w['value']:.2e} at {w['where']}")
    details = {"line": {"nodes": line["nodes"], "members": line["members"]}, "combinations": per_combo,
               "notRepresentable": not_rep, "failed": failed,
               "notCovered": [f"buckling of {c}" for c in results.buckling]}
    if not per_combo:
        why = "; ".join(f"{c}: {r}" for c, r in not_rep.items())
        return TierResult("A", False, None, f"no combination the ODE oracle can represent ({why})", [], {},
                          details, time.perf_counter() - t0)
    largest.sort(key=lambda d: -d["value"])
    # a combination the oracle cannot represent is left to the other tiers; the tier passes on the
    # combinations it covers and says which ones it does not
    passed = not failed
    summary = (f"{len(per_combo)} combination(s) against the ODE oracle: " + ("agree" if passed else "DIFFER: "
               + "; ".join(failed[:3])))
    if not_rep:
        summary += f"; {len(not_rep)} combination(s) not representable (left to the other tiers)"
    if results.buckling:
        summary += "; buckling factors are not covered by the automatic oracle check"
    return TierResult("A", True, passed, summary, largest[:6], dict(POLICY["A"]), details, time.perf_counter() - t0)
