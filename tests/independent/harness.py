"""Shared helpers of the independent review tests: model builder, benchfea runner, comparison."""
from __future__ import annotations

import copy

import numpy as np

from benchfea.model import build_model
from benchfea.solve import analyze

from tableau import DOF, Tableau

FORCE = {"N": 1.0, "kN": 1.0e3, "lbf": 4.4482216152605, "kip": 4448.2216152605}
LENGTH = {"mm": 1.0, "m": 1.0e3, "in": 25.4, "ft": 304.8}


def base(units=("kN", "m"), vertical="Z", **analysis):
    d = {"formatVersion": "1.0.0", "name": "independent review model",
         "units": {"force": units[0], "length": units[1]}, "conventions": {"verticalAxis": vertical},
         "nodes": {}, "supports": {}, "materials": {}, "sections": {}, "members": {}, "constraints": {},
         "loadCases": {}, "combinations": {}}
    if analysis:
        d["analysis"] = analysis
    return d


def nodes(d, **xyz):
    for k, v in xyz.items():
        d["nodes"][k] = {"x": v[0], "y": v[1], "z": v[2]}


def member(d, mid, i, j, sec="S", mat="M", **kw):
    d["members"][mid] = {"i": i, "j": j, "section": sec, "material": mat, **kw}


def case(d, cid, *loads, factor=1.0, combo=None, **combo_kw):
    d["loadCases"][cid] = {"loads": list(loads)}
    d["combinations"][combo or cid] = {"factors": {cid: factor}, **combo_kw}


def run_benchfea(raw, strict=True):
    model = build_model(copy.deepcopy(raw))
    res = analyze(model, strict=strict)
    return model, res, res.to_document(model)


def rel(a, b, scale=None):
    a, b = np.asarray(a, float), np.asarray(b, float)
    s = scale if scale is not None else max(np.abs(a).max(initial=0.0), np.abs(b).max(initial=0.0))
    return 0.0 if s == 0.0 else float(np.abs(a - b).max(initial=0.0) / s)


def stations_of(doc_combo):
    """{member: [(x, side)]} from a benchfea output combination: a repeated x is (left, right); x = 0 is
    the limit from outside (-1), x = L the limit from outside (+1)."""
    out = {}
    for mid, m in doc_combo["members"].items():
        xs = m["stations"]["x"]
        st = []
        for k, x in enumerate(xs):
            dup_prev = k > 0 and xs[k - 1] == x
            st.append((x, +1 if (dup_prev or k == len(xs) - 1) else -1))
        out[mid] = st
    return out


def compare(raw, cid, second_order=False, doc=None, nsub=1, skip_members=(), arm=True):
    """Solve with benchfea and with the tableau reference; return the largest relative differences per
    quantity (each relative to the largest value of that quantity kind in the model)."""
    if doc is None:
        _, _, doc = run_benchfea(raw)
    c = doc["combinations"][cid]
    t = Tableau(raw)
    t.nsub, t.arm = nsub, arm
    lunit = LENGTH[raw["units"]["length"]]
    st = stations_of(c)
    # station x was converted mm -> file unit by benchfea; snap the ends exactly
    for mid, lst in st.items():
        L = t.members[mid].L
        st[mid] = [(0.0 if k == 0 else (L if k == len(lst) - 1 else min(x, L)), s) for k, (x, s) in enumerate(lst)]
    # ... and snap stations onto the reference's own load positions (they differ by round-off of the
    # unit conversion; a station must sit exactly on a point load to see both sides of the jump)
    t.set_combination(cid)
    for mid, lst in st.items():
        m = t.members[mid]
        bs = np.array(m._breaks())
        snapped = []
        for x, s_ in lst:
            k = int(np.abs(bs - x).argmin())
            snapped.append((float(bs[k]) if abs(bs[k] - x) < 1e-9 * m.L else x, s_))
        st[mid] = snapped
    ref = t.solve(cid, second_order=second_order, stations=st)
    ids = list(raw["nodes"])
    U = np.array([[c["displacements"][n][d] for d in DOF] for n in ids])
    F6 = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")
    R = np.array([[c["reactions"].get(n, {}).get(f, 0.0) for f in F6] for n in ids])
    size = max(np.ptp(np.array([list(t.xyz[n]) for n in ids]), axis=0).max(), 1e-300)
    out = {}
    out["translations"] = rel(U[:, :3], ref["displacements"][:, :3])
    out["rotations"] = rel(U[:, 3:], ref["displacements"][:, 3:])
    fs = max(np.abs(R[:, :3]).max(), np.abs(ref["reactions"][:, :3]).max(),
             max((np.abs(m["end_forces"][[0, 1, 2, 6, 7, 8]]).max() for m in ref["members"].values()), default=0.0))
    ms = max(np.abs(R[:, 3:]).max(), np.abs(ref["reactions"][:, 3:]).max(),
             max((np.abs(m["end_forces"][[3, 4, 5, 9, 10, 11]]).max() for m in ref["members"].values()), default=0.0),
             fs * size * 1e-3)
    out["reaction forces"] = rel(R[:, :3], ref["reactions"][:, :3], fs)
    out["reaction moments"] = rel(R[:, 3:], ref["reactions"][:, 3:], ms)
    ef = em = sf = sm = sd = 0.0
    dscale = max(np.abs(U[:, :3]).max(), np.abs(ref["displacements"][:, :3]).max(), 1e-300)
    for mid, m in c["members"].items():
        if mid in skip_members:
            continue
        r = ref["members"][mid]
        e = np.array([m["endForces"][end][f] for end in "ij" for f in F6])
        ef = max(ef, rel(e[[0, 1, 2, 6, 7, 8]], r["end_forces"][[0, 1, 2, 6, 7, 8]], fs))
        em = max(em, rel(e[[3, 4, 5, 9, 10, 11]], r["end_forces"][[3, 4, 5, 9, 10, 11]], ms))
        if r.get("slack"):
            assert m["active"] is False, f"{mid}: reference slack, benchfea active"
            continue
        assert m["active"] is True, f"{mid}: benchfea slack, reference active"
        s = m["stations"]
        S = np.array([s[k] for k in ("N", "Vy", "Vz", "T", "My", "Mz")]).T
        D = np.array([s[k] for k in ("ux", "uy", "uz")]).T
        sf = max(sf, rel(S[:, :3], r["forces"][:, :3], fs))
        sm = max(sm, rel(S[:, 3:], r["forces"][:, 3:], ms))
        sd = max(sd, rel(D, r["displacements"], dscale))
    out.update({"end forces": ef, "end moments": em, "station forces": sf, "station moments": sm,
                "station displacements": sd})
    out["_ref"], out["_doc"], out["_tableau"] = ref, c, t
    return out


def worst(out):
    return max(v for k, v in out.items() if not k.startswith("_"))


def doc_arrays(raw, c):
    ids = list(raw["nodes"])
    F6 = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")
    U = np.array([[c["displacements"][n][d] for d in DOF] for n in ids])
    R = np.array([[c["reactions"].get(n, {}).get(k, 0.0) for k in F6] for n in ids])
    E = {mid: np.array([m["endForces"][e][k] for e in "ij" for k in F6]) for mid, m in c["members"].items()}
    return U, R, E


def compare_foreign(raw, c, f, shears=True):
    """Largest relative differences between a benchfea output combination ``c`` and a foreign result ``f``
    (reactions that the foreign solver does not report are NaN and skipped)."""
    U, R, E = doc_arrays(raw, c)
    mask = ~np.isnan(f["reactions"])
    Rf = np.where(mask, f["reactions"], 0.0)
    Rb = np.where(mask, R, 0.0)
    ids = list(E)
    Eb, Ef = np.array([E[m] for m in ids]), np.array([f["end_forces"][m] for m in ids])
    fi, mi = ([0, 1, 2, 6, 7, 8] if shears else [0, 6]), [3, 4, 5, 9, 10, 11]
    return {"translations": rel(U[:, :3], f["displacements"][:, :3]), "rotations": rel(U[:, 3:], f["displacements"][:, 3:]),
            "reaction forces": rel(Rb[:, :3], Rf[:, :3]), "reaction moments": rel(Rb[:, 3:], Rf[:, 3:]),
            "end forces": rel(Eb[:, fi], Ef[:, fi]), "end moments": rel(Eb[:, mi], Ef[:, mi])}
