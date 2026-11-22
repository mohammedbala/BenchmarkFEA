"""Required strengths of a design member from the analysis results.

A design member is an ordered chain of collinear analysis members with the
same local axes, so their station forces ``(N, Vy, Vz, T, My, Mz)`` (stress
resultants on the +x face, docs/CONVENTIONS.md section 4; N > 0 = tension) can
be read along one axis: x = distance from the flexible start of the first
member.
"""
from __future__ import annotations

import numpy as np

from .base import Demand

N, VY, VZ = 0, 1, 2


def member_stations(model, combination_result, design_member_id: str):
    """Stations of a design member in one combination -> (x (n,), forces (n, 6)), x from the
    design member start, members concatenated in order (junction stations appear twice)."""
    p = model.design_members[design_member_id]
    first = model.members[p.members[0]]
    xs, fs = [], []
    for mid in p.members:
        m, r = model.members[mid], combination_result.members[mid]
        xs.append(float(np.linalg.norm(m.start - first.start)) + np.asarray(r.x, float))
        fs.append(np.asarray(r.forces, float))
    return np.concatenate(xs), np.vstack(fs)


def distance_to_zero(x: np.ndarray, s: np.ndarray, k: int) -> float | None:
    """Distance from station ``k`` to the nearest point where the (piecewise linear)
    signed quantity ``s`` is zero; ``None`` when it does not vanish (or ``s[k]`` is zero)."""
    if s[k] == 0.0:
        return None
    best = None
    for a in range(len(x) - 1):
        sa, sb = s[a], s[a + 1]
        if sa == 0.0:
            x0 = x[a]
        elif sb == 0.0:
            x0 = x[a + 1]
        elif sa * sb < 0.0:
            x0 = x[a] + (x[a + 1] - x[a]) * sa / (sa - sb)
        else:
            continue
        dist = abs(float(x0) - float(x[k]))
        best = dist if best is None else min(best, dist)
    return best


def demand_from_stations(x, forces, combination: str = "") -> Demand:
    """Maximum tension, compression and shears over the stations of one combination."""
    x, f = np.asarray(x, float), np.asarray(forces, float)
    n, vy, vz = f[:, N], f[:, VY], f[:, VZ]
    kt, kc = int(np.argmax(n)), int(np.argmin(n))
    ky, kz = int(np.argmax(np.abs(vy))), int(np.argmax(np.abs(vz)))
    v = np.hypot(vy, vz)
    kv = int(np.argmax(v))
    # resultant: zero of the shear vector projected on the direction of the maximum shear
    proj = (vy * vy[kv] + vz * vz[kv]) / v[kv] if v[kv] > 0.0 else np.zeros_like(v)

    def at(k, value):
        return float(x[k]) if value > 0.0 else None

    Pt, Pc = max(float(n[kt]), 0.0), max(float(-n[kc]), 0.0)
    Vy, Vz, V = float(abs(vy[ky])), float(abs(vz[kz])), float(v[kv])
    return Demand(combination=combination, Pt=Pt, x_Pt=at(kt, Pt), Pc=Pc, x_Pc=at(kc, Pc),
                  Vy=Vy, x_Vy=at(ky, Vy), Lv_y=distance_to_zero(x, vy, ky),
                  Vz=Vz, x_Vz=at(kz, Vz), Lv_z=distance_to_zero(x, vz, kz),
                  V=V, x_V=at(kv, V), Lv=distance_to_zero(x, proj, kv))


def member_demands(model, results, design_member_id: str, purpose: str | None = "strength") -> dict:
    """{combination id: Demand} of a design member, for the combinations of ``purpose``
    (``None``: all) that are present in ``results``."""
    out = {}
    for cid, combo in model.combinations.items():
        if (purpose is None or combo.purpose == purpose) and cid in results.combinations:
            x, f = member_stations(model, results.combinations[cid], design_member_id)
            out[cid] = demand_from_stations(x, f, cid)
    return out
