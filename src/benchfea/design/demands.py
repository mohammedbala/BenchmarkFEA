"""Required strengths of a design member from the analysis results.

A design member is an ordered chain of collinear analysis members with the
same local axes, so their station forces ``(N, Vy, Vz, T, My, Mz)`` (stress
resultants on the +x face, docs/CONVENTIONS.md section 4; N > 0 = tension) can
be read along one axis: x = distance from the flexible start of the first
member.

Zero of a computed force (:data:`ZERO_FORCE`, :func:`clean_stations`).  The design checks take discrete
decisions from the station forces: which flange is in compression and whether the moment changes sign
in an unbraced segment (Chapter F), whether a member is in tension or in compression and whether it
bends about an axis at all (Chapters D, E, H), where the shear vanishes (G).  A station force that is
zero in exact arithmetic is returned by the analysis as a rounding error of either sign (a moment of
-2e-7 N mm at a released beam end next to a span moment of 2e8 N mm), and a decision taken on its sign
or on ``!= 0`` is then arbitrary.  The decisions therefore see the stations AFTER this rule:

    a force component (N, Vy, Vz) with |value| <= ZERO_FORCE x F_ref is zero,
    a moment component (T, My, Mz) with |value| <= ZERO_FORCE x M_ref is zero,

    F_ref = max(largest |N|, |Vy|, |Vz| of the design member in the combination, M_max / L),
    M_ref = max(largest |T|, |My|, |Mz| of the design member in the combination, F_max x L),

L = length of the design member.  The cross terms are the force that would produce the largest moment
over the member length and the moment the largest force would produce over it, so a member with
moments only (or forces only) still has a scale of the other kind.

ZERO_FORCE = 1e-9.  The station forces are sums of products stiffness x displacement evaluated in double
precision; their rounding error is bounded by ROUND_OFF = 16 x 2.2e-16 = 3.6e-15 times the size S of the
terms summed (``benchfea.solve.invariants``; measured: 0.6 x 2.2e-16 x S).  Relative to the member's own
forces this is 3.6e-15 x S / F: 1e-15 in the case that showed the defect, 4e-11 for S / F = 1e4 (an
ordinary frame, axial stiffness against bending), and it reaches 1e-9 only for S / F = 3e5 (bound) to
7e6 (measured error), i.e. for members whose forces are six orders below stiffness x displacement.  The
threshold is therefore above the rounding error except in such extremely stiff models, where a rounding
error can survive the rule (the decision is then as arbitrary as it was before, never worse).  On the
other side, a real force of 1e-9 of the largest one changes no design ratio by more than about that
fraction, so the rule cannot hide anything that matters.  Nothing larger than the threshold is altered.
"""
from __future__ import annotations

import numpy as np

from .base import Demand

N, VY, VZ = 0, 1, 2
#: a station force / moment below this fraction of the member's force / moment scale is zero (module docstring)
ZERO_FORCE = 1.0e-9


def force_scales(x, forces) -> tuple[float, float]:
    """(F_ref, M_ref) of the stations of one design member in one combination (module docstring)."""
    x, f = np.asarray(x, float), np.asarray(forces, float)
    if not f.size:
        return 0.0, 0.0
    L = float(x[-1] - x[0]) if len(x) else 0.0
    fmax, mmax = float(np.abs(f[:, :3]).max()), float(np.abs(f[:, 3:]).max())
    if L > 0.0:
        return max(fmax, mmax / L), max(mmax, fmax * L)
    return fmax, mmax


def clean_stations(x, forces) -> np.ndarray:
    """Copy of the station forces (n, 6) with the components that are zero by the rule of the module
    docstring set to exactly 0.0."""
    f = np.array(forces, float)
    fref, mref = force_scales(x, f)
    f[:, :3][np.abs(f[:, :3]) <= ZERO_FORCE * fref] = 0.0
    f[:, 3:][np.abs(f[:, 3:]) <= ZERO_FORCE * mref] = 0.0
    return f


def member_stations(model, combination_result, design_member_id: str):
    """Stations of a design member in one combination -> (x (n,), forces (n, 6)), x from the
    design member start, members concatenated in order (junction stations appear twice).  The forces
    are returned after :func:`clean_stations` (rounding errors of zero set to zero)."""
    p = model.design_members[design_member_id]
    first = model.members[p.members[0]]
    xs, fs = [], []
    for mid in p.members:
        m, r = model.members[mid], combination_result.members[mid]
        xs.append(float(np.linalg.norm(m.start - first.start)) + np.asarray(r.x, float))
        fs.append(np.asarray(r.forces, float))
    x = np.concatenate(xs)
    return x, clean_stations(x, np.vstack(fs))


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
