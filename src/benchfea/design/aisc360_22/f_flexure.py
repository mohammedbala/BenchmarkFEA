"""AISC 360-22 Chapter F -- members in flexure (pp. 16.1-50 to 16.1-65).

F1(a): phi_b = 0.90 (LRFD), Omega_b = 1.67 (ASD) for every limit state.

Implemented (one check per limit state; the largest ratio governs, which is
"the lower/lowest value" of the nominal strengths of each section):

* F2  doubly symmetric I-shapes, compact web and flanges, major axis:
  yielding Mn = Mp = Fy Zx (Eq. F2-1); lateral-torsional buckling (F2.2):
  Lb <= Lp does not apply; Lp < Lb <= Lr Eq. F2-2; Lb > Lr Eq. F2-3 with Fcr of
  Eq. F2-4; Lp Eq. F2-5, Lr Eq. F2-6, rts Eq. F2-7, c = 1 (Eq. F2-8a).
* F3  same with noncompact / slender flanges: LTB as F2.2; compression flange
  local buckling Eq. F3-1 (noncompact), Eq. F3-2 with kc (slender).
* F6  I-shapes, minor axis: yielding Mp = Fy Zy <= 1.6 Fy Sy (Eq. F6-1); flange
  local buckling Eq. F6-2 (noncompact), F6-3 / F6-4 (slender).
* F7  rectangular HSS and box sections, either axis: yielding (Eq. F7-1);
  flange local buckling Eq. F7-2 (noncompact), Eq. F7-3 with the effective
  width of Eq. F7-4 (HSS) / F7-5 (box) (slender); web local buckling Eq. F7-6
  (noncompact); LTB Eq. F7-8 to F7-11.
* F8  round HSS with D/t < 0.45 E/Fy: yielding (Eq. F8-1), local buckling
  Eq. F8-2 (noncompact), F8-3 / F8-4 (slender).
* F1  Cb by Eq. F1-1 from the moment diagram of each unbraced segment, or Cb
  given (design variable ``Cb``).

NOT implemented -> ``notChecked`` with the reason: I-shapes with noncompact
webs (F4) or slender webs (F5); rectangular sections with slender webs
(Eq. F7-7: the use of Eq. F5-6 for a box is not fully defined); round HSS with
D/t >= 0.45 E/Fy; every other shape (F9 to F12); F13 (holes, proportions).

Axes.  AISC x = major = bending about local z (station moment ``Mz``), AISC y
= minor = bending about local y (``My``); the mapping is fixed
(docs/CONVENTIONS.md section 3).  ``Mz > 0`` puts the +y ("top") flange in
compression (sigma = N/A - Mz y/Iz + My z/Iy).  An I-shape whose inertia about
local z is not the larger one is ``notChecked`` (F2/F3/F6 presume x = major).

Required flexural strength.  ``Demand.extra`` carries either the stations of
the combination (``"x"``, ``"forces"`` (n, 6) = N, Vy, Vz, T, My, Mz) or the
numbers ``"Mrx"``, ``"Mry"`` (and ``"Tr"``).  Without either, the demand says
nothing about flexure and NO Chapter F check is produced (the orchestration
then reports that flexure was not evaluated).  Yielding and local buckling
are checked with the largest moment of the member; LTB with the largest
moment of each unbraced segment (Commentary F1: Mn is compared with the
maximum moment in the unbraced segment).

Moment diagram between stations (derived here, no reference used).  For a
straight member, moment equilibrium of a slice dx about local z gives
dMz/dx = -Vy and about local y gives dMy/dx = +Vz (docs/CONVENTIONS.md
section 4; first-order equilibrium), and dV/dx = -w.  Between two stations
the moment is therefore a polynomial whose end values AND end slopes are
known from the station results; the cubic through these four numbers (cubic
Hermite interpolation) is the exact moment for any load that varies linearly
between the stations (w linear -> V quadratic -> M cubic).  Stations exist at
every point load and at the ends of every line load, so the interpolation is
exact for the loads of the input format in a first-order analysis.  In a
second-order analysis dM/dx also contains N dv/dx, so the interpolation is
then an approximation between stations (the station values themselves are
used unchanged).  Tests compare it with closed-form beam moments.

Unbraced segments for LTB (F2.2: Lb = length between points braced against
lateral displacement of the compression flange or against twist):

* ``Lb`` given and equal to the member length (the default): one segment, the
  whole member.
* ``Lb`` given and different from the member length: the positions of the
  brace points are not known, so Cb = 1.0 unless given, with the largest
  moment of the member (Commentary F1: Cb = 1.0 is conservative, except for
  some unbraced overhangs and top-flange loading -- see LIMITS).
* brace points given (``bracing.topFlange`` / ``bottomFlange``, the member
  ends being brace points by definition of the input format): per flange, the
  segments between ITS brace points.  A segment is checked for that flange
  when the flange is in compression over the whole segment.  When the moment
  changes sign inside a segment of a flange, that flange's brace points do
  not brace the compression flange everywhere; the segment between the
  nearest points where BOTH flanges are braced is then used (engine rule:
  only restraint that is certain is counted).  Eq. F1-1 itself covers reverse
  curvature of doubly symmetric members (F1(c)).
  Twist brace points (``bracing.torsion``) are NOT counted as LTB braces.
* F1-1 requires both ends of the segment braced.  A segment that ends at a
  FREE end of the member (no support, constraint or other member at the node)
  is not such a segment.  F1 gives one rule for it: cantilevers with warping
  prevented at the support and the free end unbraced, Cb = 1.0.  This is
  applied when the segment is the whole member and its other end is a
  built-in support (all six DOFs fixed; the frame model has no warping DOF,
  so "warping prevented" is ASSUMED for a built-in support and said so in the
  check), with Lb = the cantilever length.  Any other free-end segment
  (overhangs, a free end beyond a brace point) is ``notChecked`` unless both
  ``Cb`` and ``Lb`` are given.
* a segment without moment is not checked (Eq. F1-1 is 0/0 there).

Choices where the Specification leaves room (reported in each check):

* rts: Eq. F2-7 (needs Iy, Cw) unless a tabulated ``rts`` is given; ho = d - tf
  (distance between the centroids of two equal rectangular flanges) unless
  given.
* F7.2(c) effective section modulus Se: the Specification defines only be.  Se
  is computed by removing the ineffective width b - be from BOTH flanges (no
  shift of the neutral axis): Ieff = I - 2 [(b - be) t^3/12 + (b - be) t
  ((H - t)/2)^2], Se = Ieff / (H/2), H = outside dimension in the plane of
  bending.  The Commentary to F7 calls this alternative slightly conservative.
* F7.4 LTB is evaluated only for bending about local z when Ix > Iy (the
  bracing data refers to the +y / -y flanges); bending about the axis of
  smaller or equal inertia has no LTB (F7.4 User Note).  A rectangular section
  bent about local y with Iy > Ix is ``notChecked`` for LTB.
* round HSS: the RESULTANT moment sqrt(Mx^2 + My^2) is used (every axis of a
  circular section is a principal axis, so two moment components are one
  uniaxial moment about the axis of their resultant).
* the optional Cb increase for concurrent tension (H1.2, Eq. H1-2) is not applied.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..base import LENGTH, NONE, STRESS, TEXT, DesignContext, Demand, missing, not_checked, positive, \
    strength_check
from ..section import BOX, I_SHAPE, RECT_HSS, ROUND_HSS
from . import b4_classification as b4

PHI_B, OMEGA_B = 0.90, 1.67                 # F1(a)
MOMENT = "moment"
MODULUS = "length^3"

NAME_Y = {"x": "flexure, major axis: yielding", "y": "flexure, minor axis: yielding"}
NAME_FLB = {"x": "flexure, major axis: compression flange local buckling",
            "y": "flexure, minor axis: flange local buckling"}
NAME_WLB = {"x": "flexure, major axis: web local buckling", "y": "flexure, minor axis: web local buckling"}
NAME_LTB = "flexure, major axis: lateral-torsional buckling"
NAME_LTB_MINOR = "flexure, minor axis: lateral-torsional buckling"
NAME_Y_ROUND = "flexure (resultant moment, round section): yielding"
NAME_LB_ROUND = "flexure (resultant moment, round section): local buckling"
CLAUSE = "AISC 360-22 Chapter F"

_AXIS_TEXT = {"x": "major axis (bending about local z)", "y": "minor axis (bending about local y)"}


# ============================================================================== moment diagrams
class Diagram:
    """One bending moment along the design member: station values ``m`` and slopes ``dm`` = dM/dx
    (from the shear forces), interpolated by the cubic through values and slopes (module docstring)."""

    def __init__(self, x, m, dm):
        self.x, self.m, self.dm = (np.asarray(v, float) for v in (x, m, dm))

    def _pieces(self):
        for k in range(len(self.x) - 1):
            h = self.x[k + 1] - self.x[k]
            if h > 0.0:
                yield k, h

    @staticmethod
    def _coefficients(m0, d0, m1, d1, h):
        """Cubic p(s), s in [0, 1], with p(0) = m0, p(1) = m1, p'(0) = h d0, p'(1) = h d1."""
        return (2.0 * m0 + h * d0 - 2.0 * m1 + h * d1, -3.0 * m0 - 2.0 * h * d0 + 3.0 * m1 - h * d1, h * d0, m0)

    def at(self, x: float) -> float:
        """Moment at ``x``; where two stations share ``x`` (a jump) the value of larger magnitude."""
        best = None
        for k, h in self._pieces():
            if self.x[k] <= x <= self.x[k + 1]:
                s = (x - self.x[k]) / h
                a3, a2, a1, a0 = self._coefficients(self.m[k], self.dm[k], self.m[k + 1], self.dm[k + 1], h)
                v = ((a3 * s + a2) * s + a1) * s + a0
                if s == 0.0:
                    v = self.m[k]
                elif s == 1.0:
                    v = self.m[k + 1]
                best = v if best is None or abs(v) > abs(best) else best
        if best is None:
            raise ValueError(f"x = {x} is outside the stations")
        return float(best)

    def points(self, a: float, b: float) -> list:
        """(x, M) at ``a``, ``b``, every station in between and every interior extremum of the cubics."""
        out = [(a, self.at(a)), (b, self.at(b))]
        out += [(float(x), float(m)) for x, m in zip(self.x, self.m) if a <= x <= b]
        for k, h in self._pieces():
            lo, hi = max(a, self.x[k]), min(b, self.x[k + 1])
            if hi <= lo:
                continue
            a3, a2, a1, a0 = self._coefficients(self.m[k], self.dm[k], self.m[k + 1], self.dm[k + 1], h)
            for s in _quadratic_roots(3.0 * a3, 2.0 * a2, a1):
                xs = self.x[k] + s * h
                if 0.0 < s < 1.0 and lo < xs < hi:
                    out.append((float(xs), float(((a3 * s + a2) * s + a1) * s + a0)))
        return out

    def extremes(self, a: float, b: float):
        """-> (largest |M|, its x, smallest signed M, largest signed M) in [a, b]."""
        pts = self.points(a, b)
        x_abs, m_abs = max(pts, key=lambda p: abs(p[1]))
        return abs(m_abs), x_abs, min(p[1] for p in pts), max(p[1] for p in pts)

    def cubics(self):
        """(k, x_k, h, (a3, a2, a1, a0)) of every piece between two stations with h > 0:
        M(x_k + s h) = ((a3 s + a2) s + a1) s + a0, s in [0, 1]."""
        for k, h in self._pieces():
            yield k, float(self.x[k]), float(h), self._coefficients(self.m[k], self.dm[k], self.m[k + 1],
                                                                    self.dm[k + 1], h)


def resultant_stationary_points(Mx: "Diagram", My: "Diagram") -> list:
    """Points BETWEEN the stations where the resultant moment sqrt(Mx^2 + My^2) is stationary ->
    [(x, Mx, My)].  On a piece both moments are cubics in s, so q(s) = Mx(s)^2 + My(s)^2 is a polynomial of
    degree 6 and q'(s) one of degree 5; its real roots in (0, 1) are found numerically (companion matrix,
    ``numpy.roots``).  The two diagrams have the same stations."""
    out = []
    for (_, x0, h, cx), (_, _, _, cy) in zip(Mx.cubics(), My.cubics()):
        for s in resultant_roots(cx, cy):
            out.append((x0 + s * h, float(np.polyval(cx, s)), float(np.polyval(cy, s))))
    return out


def resultant_roots(cx, cy) -> list:
    """Real roots s in (0, 1) of d/ds [px(s)^2 + py(s)^2] for two cubics given by their coefficients
    (highest power first)."""
    px, py = np.poly1d(cx), np.poly1d(cy)
    dq = (px * px + py * py).deriv()
    if not np.any(dq.coeffs != 0.0):
        return []
    return [float(r.real) for r in np.roots(dq.coeffs)
            if abs(r.imag) <= 1.0e-9 * max(1.0, abs(r.real)) and 0.0 < r.real < 1.0]


def _quadratic_roots(a: float, b: float, c: float) -> list:
    if a == 0.0:
        return [] if b == 0.0 else [-c / b]
    disc = b * b - 4.0 * a * c
    if disc < 0.0:
        return []
    r = math.sqrt(disc)
    q = -0.5 * (b + math.copysign(r, b))
    return [q / a] + ([c / q] if q != 0.0 else [])


@dataclass(frozen=True)
class Moments:
    """Flexural demand of one combination (magnitudes; ``x_*`` from the design member start)."""
    Mrx: float
    x_Mrx: float | None
    Mry: float
    x_Mry: float | None
    Mr: float                         # largest resultant sqrt(Mx^2 + My^2) (stations and between them)
    x_Mr: float | None
    Tr: float
    x_Tr: float | None
    Mx: Diagram | None = None         # signed diagrams (None when only numbers were given)
    My: Diagram | None = None
    x: np.ndarray | None = None       # stations
    N: np.ndarray | None = None
    T: np.ndarray | None = None


def moments(demand: Demand) -> Moments | None:
    """Flexural demand from ``Demand.extra``; ``None`` when the demand carries no moment data."""
    e = demand.extra or {}
    if "forces" in e and "x" in e:
        x, f = np.asarray(e["x"], float), np.asarray(e["forces"], float)
        Mx = Diagram(x, f[:, 5], -f[:, 1])          # dMz/dx = -Vy
        My = Diagram(x, f[:, 4], f[:, 2])           # dMy/dx = +Vz
        L0, L1 = float(x[0]), float(x[-1])
        mx, x_mx, _, _ = Mx.extremes(L0, L1)
        my, x_my, _, _ = My.extremes(L0, L1)
        res = np.hypot(f[:, 5], f[:, 4])
        kr, kt = int(np.argmax(res)), int(np.argmax(np.abs(f[:, 3])))
        Tr = float(abs(f[kt, 3]))
        mr, x_mr = float(res[kr]), float(x[kr])
        # the largest resultant may lie between two stations (same interpolation as for Mrx, Mry)
        for xs, vx, vy in resultant_stationary_points(Mx, My):
            if math.hypot(vx, vy) > mr:
                mr, x_mr = math.hypot(vx, vy), xs
        return Moments(mx, x_mx if mx > 0.0 else None, my, x_my if my > 0.0 else None, mr,
                       x_mr if mr > 0.0 else None, Tr, float(x[kt]) if Tr > 0.0 else None,
                       Mx, My, x, f[:, 0], f[:, 3])
    if any(k in e for k in ("Mrx", "Mry", "Tr")):
        mx, my = abs(float(e.get("Mrx", 0.0))), abs(float(e.get("Mry", 0.0)))
        return Moments(mx, None, my, None, math.hypot(mx, my), None, abs(float(e.get("Tr", 0.0))), None)
    return None


# ============================================================================== F1: Cb
def cb_f1_1(Mmax: float, MA: float, MB: float, MC: float) -> float:
    """Eq. F1-1: Cb = 12.5 Mmax / (2.5 Mmax + 3 MA + 4 MB + 3 MC) (absolute values)."""
    Mmax, MA, MB, MC = abs(Mmax), abs(MA), abs(MB), abs(MC)
    return 12.5 * Mmax / (2.5 * Mmax + 3.0 * MA + 4.0 * MB + 3.0 * MC)


def segment_cb(diagram: Diagram, a: float, b: float):
    """Eq. F1-1 for the segment [a, b] -> (Cb, Mmax, MA, MB, MC, x of Mmax); Cb is None when Mmax = 0."""
    Mmax, x_max, _, _ = diagram.extremes(a, b)
    MA, MB, MC = (abs(diagram.at(a + q * (b - a))) for q in (0.25, 0.5, 0.75))
    return (cb_f1_1(Mmax, MA, MB, MC) if Mmax > 0.0 else None), Mmax, MA, MB, MC, x_max


@dataclass(frozen=True)
class Segment:
    """One unbraced segment to be checked for LTB in one combination."""
    a: float | None                   # from / to (None: location not known)
    b: float | None
    Lb: float
    Mr: float                         # required strength = largest |M| in the segment
    x: float | None
    Cb: float | None                  # None with `blocked`
    cb_text: str
    blocked: str = ""                 # reason why LTB cannot be checked for this segment
    note: str = ""


def _cuts(fractions, L):
    return [0.0, *(f * L for f in fractions), L]


def ltb_segments(ctx: DesignContext, mo: Moments) -> list:
    """The unbraced segments of the design member with a moment, for one combination (rules in the
    module docstring)."""
    L = ctx.L
    given = ctx.Cb if positive(ctx.Cb) else None
    braces = bool(ctx.braceTop or ctx.braceBottom)
    if mo.Mx is None:                                               # numbers only: no diagram
        if ctx.Lb is not None:
            Lb = ctx.Lb
        elif braces:
            Lb = max(max(q - p for p, q in zip(c, c[1:])) for c in (_cuts(ctx.braceTop, L),
                                                                     _cuts(ctx.braceBottom, L)))
        else:
            Lb = L
        if given is not None:
            return [Segment(None, None, Lb, mo.Mrx, None, given, "Cb given")]
        return [Segment(None, None, Lb, mo.Mrx, None, 1.0, "Cb = 1.0: no moment diagram was given (Eq. F1-1 "
                        "cannot be evaluated)")]
    D = mo.Mx
    x0, span = float(D.x[0]), float(D.x[-1] - D.x[0])

    def at(a):
        """Station coordinate of the design-member position ``a`` (0 ... L); the stations span the member."""
        return x0 + min(max(a / L, 0.0), 1.0) * span

    def make(a, b, note=""):
        """Segment [a, b] (design member coordinates) with its Cb, or blocked."""
        Cb, Mmax, MA, MB, MC, x_max = segment_cb(D, at(a), at(b))
        if Mmax == 0.0:
            return None
        free = (a <= 0.0 and ctx.end_free[0]) or (b >= L and ctx.end_free[1])
        kw = dict(a=a, b=b, Lb=b - a, Mr=Mmax, x=x_max - x0, note=note)
        if given is not None:
            return Segment(Cb=given, cb_text="Cb given", **kw)
        if free:
            whole = a <= 0.0 and b >= L
            root_built_in = ctx.end_built_in[1] if ctx.end_free[0] else ctx.end_built_in[0]
            if whole and root_built_in:
                return Segment(Cb=1.0, cb_text="Cb = 1.0 (F1: cantilever with the free end unbraced; warping "
                               "ASSUMED prevented at the built-in support -- the frame model has no warping "
                               "degree of freedom; Lb = cantilever length)", **kw)
            return Segment(Cb=None, cb_text="", blocked="the unbraced segment ends at a free end of the member, "
                           "which is not a brace point: Eq. F1-1 needs both ends of the segment braced and F1 "
                           "gives Cb only for a cantilever with warping prevented at the support; give Cb and "
                           "Lb explicitly", **kw)
        return Segment(Cb=Cb, cb_text=f"Cb by Eq. F1-1 (Mmax = {Mmax:.6g}, MA = {MA:.6g}, MB = {MB:.6g}, "
                       f"MC = {MC:.6g})", **kw)

    if ctx.Lb is not None:
        if abs(ctx.Lb - L) <= 1.0e-9 * L:
            seg = make(0.0, L)
            return [seg] if seg else []
        Mmax, x_max, _, _ = D.extremes(x0, float(D.x[-1]))
        if Mmax == 0.0:
            return []
        if given is not None:
            return [Segment(None, None, ctx.Lb, Mmax, x_max - x0, given, "Cb given")]
        return [Segment(None, None, ctx.Lb, Mmax, x_max - x0, 1.0, "Cb = 1.0: Lb is given as a length that "
                        "differs from the member length, so the positions of the brace points (and the moment "
                        "diagram of each segment) are not known; give bracing.topFlange / bottomFlange for "
                        "Eq. F1-1")]
    # --- brace points: per flange, the segments between its brace points
    top, bottom = _cuts(ctx.braceTop, L), _cuts(ctx.braceBottom, L)
    both = sorted(set(top) & set(bottom))
    out, seen = [], set()

    def add(a, b, note=""):
        if (a, b) not in seen:
            seen.add((a, b))
            seg = make(a, b, note)
            if seg:
                out.append(seg)

    for cuts, sign, name in ((top, 1.0, "top"), (bottom, -1.0, "bottom")):
        for a, b in zip(cuts, cuts[1:]):
            _, _, lo, hi = D.extremes(at(a), at(b))
            lo, hi = (lo, hi) if sign > 0.0 else (-hi, -lo)         # as seen by this flange: > 0 = compression
            if hi <= 0.0:
                continue                                            # this flange is never in compression here
            if lo >= 0.0:
                add(a, b)                                           # in compression throughout: braced at both ends
                continue
            a2 = max(p for p in both if p <= a)
            b2 = min(p for p in both if p >= b)
            add(a2, b2, "" if (a2, b2) == (a, b) else
                f"the moment changes sign between the {name}-flange brace points at {a:.6g} and {b:.6g}: "
                f"the segment between points where both flanges are braced ({a2:.6g} to {b2:.6g}) is "
                "used (engine rule)")
    return out


# ============================================================================== F2 / F3 pure functions
def lp_i(ry: float, E: float, Fy: float) -> float:
    """Eq. F2-5: Lp = 1.76 ry sqrt(E/Fy)."""
    return 1.76 * ry * math.sqrt(E / Fy)


def lr_i(rts: float, E: float, Fy: float, J: float, c: float, Sx: float, ho: float) -> float:
    """Eq. F2-6: Lr = 1.95 rts E/(0.7 Fy) sqrt(Jc/(Sx ho) + sqrt((Jc/(Sx ho))^2 + 6.76 (0.7 Fy/E)^2))."""
    X = J * c / (Sx * ho)
    return 1.95 * rts * E / (0.7 * Fy) * math.sqrt(X + math.sqrt(X * X + 6.76 * (0.7 * Fy / E) ** 2))


def rts_f2_7(Iy: float, Cw: float, Sx: float) -> float:
    """Eq. F2-7: rts^2 = sqrt(Iy Cw) / Sx."""
    return math.sqrt(math.sqrt(Iy * Cw) / Sx)


def fcr_f2_4(Cb: float, E: float, Lb: float, rts: float, J: float, c: float, Sx: float, ho: float) -> float:
    """Eq. F2-4: Fcr = Cb pi^2 E / (Lb/rts)^2 sqrt(1 + 0.078 Jc/(Sx ho) (Lb/rts)^2)."""
    lam2 = (Lb / rts) ** 2
    return Cb * math.pi ** 2 * E / lam2 * math.sqrt(1.0 + 0.078 * J * c / (Sx * ho) * lam2)


def inelastic_ltb(Cb: float, Mp: float, Fy: float, Sx: float, Lb: float, Lp: float, Lr: float) -> float:
    """Eq. F2-2 (also Eq. F7-8): Mn = Cb [Mp - (Mp - 0.7 Fy Sx) (Lb - Lp)/(Lr - Lp)] <= Mp."""
    return min(Mp, Cb * (Mp - (Mp - 0.7 * Fy * Sx) * (Lb - Lp) / (Lr - Lp)))


def ltb_i_shape(E, Fy, Mp, Sx, ry, rts, J, c, ho, Lb, Cb):
    """F2.2 -> (Mn or None when LTB does not apply, equation text, Lp, Lr, Fcr or None)."""
    Lp, Lr = lp_i(ry, E, Fy), lr_i(rts, E, Fy, J, c, Sx, ho)
    if Lb <= Lp:                                                    # F2.2(a)
        return None, "F2.2(a)", Lp, Lr, None
    if Lb <= Lr:                                                    # F2.2(b)
        return inelastic_ltb(Cb, Mp, Fy, Sx, Lb, Lp, Lr), "F2-2", Lp, Lr, None
    Fcr = fcr_f2_4(Cb, E, Lb, rts, J, c, Sx, ho)                    # F2.2(c)
    return min(Fcr * Sx, Mp), "F2-3, F2-4", Lp, Lr, Fcr


def flb_i_major(E, Fy, Mp, Sx, lam, lam_p, lam_r, kc):
    """F3.2 -> (Mn, equation): Eq. F3-1 (noncompact flange) or Eq. F3-2 (slender flange)."""
    if lam <= lam_r:
        return Mp - (Mp - 0.7 * Fy * Sx) * (lam - lam_p) / (lam_r - lam_p), "F3-1"
    return 0.9 * E * kc * Sx / lam ** 2, "F3-2"


# ============================================================================== F6 pure functions
def mp_minor(Fy: float, Zy: float, Sy: float) -> float:
    """Eq. F6-1: Mn = Mp = Fy Zy <= 1.6 Fy Sy."""
    return min(Fy * Zy, 1.6 * Fy * Sy)


def flb_i_minor(E, Fy, Mp, Sy, lam, lam_p, lam_r):
    """F6.2 -> (Mn, equation): Eq. F6-2 (noncompact) or Eq. F6-3 with Fcr of Eq. F6-4 (slender)."""
    if lam <= lam_r:
        return Mp - (Mp - 0.70 * Fy * Sy) * (lam - lam_p) / (lam_r - lam_p), "F6-2"
    return 0.70 * E / lam ** 2 * Sy, "F6-3, F6-4"


# ============================================================================== F7 pure functions
def local_buckling_f7(Mp, Fy, S, lam, lam_p, lam_r) -> float:
    """Eq. F7-2 (flange) / Eq. F7-6 (web): Mn = Mp - (Mp - Fy S) (lam - lam_p)/(lam_r - lam_p) <= Mp."""
    return min(Mp, Mp - (Mp - Fy * S) * (lam - lam_p) / (lam_r - lam_p))


def effective_width_f7(b: float, t: float, E: float, Fy: float, hss: bool) -> float:
    """Eq. F7-4 (HSS, 0.38) / Eq. F7-5 (box, 0.34): be = 1.92 t sqrt(E/Fy) [1 - c/(b/t) sqrt(E/Fy)] <= b."""
    root = math.sqrt(E / Fy)
    return min(b, 1.92 * t * root * (1.0 - (0.38 if hss else 0.34) / (b / t) * root))


def effective_modulus_both_flanges(I: float, H: float, t: float, b: float, be: float) -> float:
    """Se with the ineffective width removed from both flanges (module docstring):
    Ieff = I - 2 [(b - be) t^3/12 + (b - be) t ((H - t)/2)^2], Se = Ieff / (H/2)."""
    lost = b - be
    return (I - 2.0 * (lost * t ** 3 / 12.0 + lost * t * ((H - t) / 2.0) ** 2)) / (H / 2.0)


def ltb_rectangular(E, Fy, Mp, Sx, ry, J, Ag, Lb, Cb):
    """F7.4 -> (Mn or None when LTB does not apply, equation, Lp, Lr)."""
    root = math.sqrt(J * Ag)
    Lp = 0.13 * E * ry * root / Mp                                  # Eq. F7-10
    Lr = 2.0 * E * ry * root / (0.7 * Fy * Sx)                      # Eq. F7-11
    if Lb <= Lp:                                                    # F7.4(a)
        return None, "F7.4(a)", Lp, Lr
    if Lb <= Lr:                                                    # Eq. F7-8
        return inelastic_ltb(Cb, Mp, Fy, Sx, Lb, Lp, Lr), "F7-8", Lp, Lr
    return min(Mp, 2.0 * E * Cb * root / (Lb / ry)), "F7-9", Lp, Lr     # Eq. F7-9


# ============================================================================== F8 pure functions
def round_local_buckling(E, Fy, S, D_over_t, slender: bool):
    """F8.2 -> (Mn, equation): Eq. F8-2 (noncompact) or Eq. F8-3 with Fcr of Eq. F8-4 (slender walls)."""
    if not slender:
        return (0.021 * E / D_over_t + Fy) * S, "F8-2"
    return 0.33 * E / D_over_t * S, "F8-3, F8-4"


# ============================================================================== checks
@dataclass
class AxisResult:
    """Chapter F checks of one bending axis in one combination."""
    checks: list
    Mc: float | None = None           # smallest available strength over the limit states (and LTB segments)
    blocked: str = ""                 # why Mc is not available
    governing: str = ""               # name of the limit state that gives Mc


def _finish(checks, extra_capacities=()) -> AxisResult:
    bad = [c for c in checks if c.status != "ok"]
    if bad:
        return AxisResult(checks, None, f"{bad[0].name}: {bad[0].reason}")
    caps = [(c.capacity, c.name) for c in checks] + list(extra_capacities)
    if not caps:
        return AxisResult(checks, None, "no flexural limit state was evaluated")
    Mc, name = min(caps, key=lambda p: p[0])
    return AxisResult(checks, Mc, "", name)


def _strength(name, clause, ctx, Mn, Mr, x, combination, inputs, details, notes=()):
    return strength_check(name, f"AISC 360-22 {clause}; F1(a)", ctx, nominal=Mn, phi=PHI_B, omega=OMEGA_B,
                          demand=Mr, dimension=MOMENT, combination=combination, x=x, inputs=inputs,
                          details=details, notes=list(notes))


def _basic(ctx: DesignContext, *moduli: str):
    sec = ctx.section
    if sec.shape is None:
        return sec.reason
    miss = missing(ctx, "E", "Fy")
    if miss:
        return f"no {' and '.join(miss)} (material or design variable)"
    lacking = [n for n in moduli if not positive(getattr(sec, n))]
    if lacking:
        return f"section {sec.id!r} gives no {', '.join(lacking)} (section moduli needed by Chapter F)"
    return None


def _ltb_checks(ctx, mo, combination, capacity_of, clause_of, inputs, details, notes):
    """LTB over the unbraced segments -> (check or None, [(capacity, name)] of all segments, blocked reason).

    ``capacity_of(Lb, Cb) -> (Mn or None, equation, {details})``."""
    segments = ltb_segments(ctx, mo)
    rows, caps = [], []
    blocked = [s for s in segments if s.blocked]
    if blocked:
        s = max(blocked, key=lambda s: s.Mr)
        return not_checked(NAME_LTB, f"AISC 360-22 {clause_of('')}; F1", s.blocked, ctx, demand=s.Mr,
                           dimension=MOMENT, combination=combination, x=s.x, inputs=inputs,
                           details={**details, "segment": (_seg_text(s), TEXT)}, notes=list(notes)), [], s.blocked
    for s in segments:
        Mn, eq, more = capacity_of(s.Lb, s.Cb)
        rows.append((s, Mn, eq, more))
    summary = "; ".join(f"{_seg_text(s)}: Lb = {s.Lb:.6g}, Cb = {s.Cb:.4g}, Mr = {s.Mr:.6g}, "
                        + ("LTB does not apply" if Mn is None else f"Mn = {Mn:.6g} (Eq. {eq})")
                        for s, Mn, eq, more in rows)
    applies = [r for r in rows if r[1] is not None]
    if not applies:
        return None, [], ""
    factor = PHI_B if ctx.method == "LRFD" else 1.0 / OMEGA_B
    caps = [(factor * Mn, NAME_LTB) for _, Mn, _, _ in applies]
    s, Mn, eq, more = max(applies, key=lambda r: r[0].Mr / r[1])
    seg_notes = [n for n in (s.note, s.cb_text) if n]
    check = _strength(NAME_LTB, f"{clause_of(eq)}; F1" + (", Eq. F1-1" if "F1-1" in s.cb_text else ""), ctx, Mn,
                      s.Mr, s.x, combination,
                      {**inputs, "Lb": (s.Lb, LENGTH, _lb_source(ctx, s)), "Cb": (s.Cb, NONE, _cb_source(ctx, s))},
                      {**details, **more, "governing segment": (_seg_text(s), TEXT),
                       "all segments": (summary, TEXT)}, [*notes, *seg_notes])
    return check, caps, ""


def _seg_text(s: Segment) -> str:
    return "location not known" if s.a is None else f"from {s.a:.6g} to {s.b:.6g}"


def _lb_source(ctx, s):
    if ctx.Lb is not None:
        return ctx.source("Lb")
    if ctx.braceTop or ctx.braceBottom:
        return "bracing: segment between brace points"
    return "derived: design member length"


def _cb_source(ctx, s):
    if positive(ctx.Cb):
        return ctx.source("Cb")
    return "derived: Eq. F1-1" if "F1-1" in s.cb_text else "derived: Cb = 1.0 (see notes)"


def _no_ltb_detail(ctx, mo) -> dict:
    """Detail entry stating why no LTB check appears."""
    if not ltb_segments(ctx, mo):
        return {"LTB": ("no moment about the major axis in any unbraced segment", TEXT)}
    return {"LTB": ("does not apply: Lb <= Lp in every unbraced segment with a moment", TEXT)}


# ------------------------------------------------------------------------------ I-shapes
def _i_major(ctx: DesignContext, mo: Moments, combination: str) -> AxisResult:
    sec, E, Fy = ctx.section, ctx.E, ctx.Fy
    kw = dict(demand=mo.Mrx, dimension=MOMENT, combination=combination, x=mo.x_Mrx)
    reason = _basic(ctx, "Zx")
    if not reason and positive(sec.Ix) and positive(sec.Iy) and not sec.Ix > sec.Iy:
        reason = (f"section {sec.id!r}: the inertia about local z (Ix = {sec.Ix:g}) is not larger than about "
                  f"local y (Iy = {sec.Iy:g}); AISC 360-22 F2, F3 and F6 presume x = major axis")
    cls = None
    if not reason:
        cls = b4.classify_flexure(sec, E, Fy, b4.MAJOR)
        reason = cls.reason if cls.status != "ok" else None
    if reason:
        return _finish([not_checked(NAME_Y["x"], "AISC 360-22 F2, F3", reason, ctx, **kw)])
    flange, web = cls.elements
    if web.cls != b4.COMPACT:
        section = "F4 (noncompact web)" if web.cls == b4.NONCOMPACT else "F5 (slender web)"
        return _finish([not_checked(
            NAME_Y["x"], f"AISC 360-22 {section.split()[0]}", f"the web of section {sec.id!r} is {web.cls} in "
            f"flexure (h/tw = {web.ratio:.4g} > lambda_pw = {web.lambda_p:.4g}, Table B4.1b Case 15): Section "
            f"{section} is not implemented", ctx, details=cls.details(), **kw)])
    Zx, Sx = sec.Zx, sec.Sx
    Mp = Fy * Zx                                                    # Eq. F2-1
    inputs = {**ctx.inputs("E", "Fy"), **sec.inputs("Zx", "shape", "fabrication"), **cls.inputs()}
    if positive(Sx):
        inputs.update(sec.inputs("Sx"))
    details = {**cls.details(), "Mp = Fy Zx (Eq. F2-1)": (Mp, MOMENT)}
    compact = flange.cls == b4.COMPACT
    base = "F2" if compact else "F3"
    checks = []

    # --- lateral-torsional buckling (F2.2; F3.1 refers to F2.2)
    ry = sec.ry
    rts, ho = sec.dims.get("rts"), sec.dims.get("ho")
    ltb_inputs, ltb_details = dict(inputs), dict(details)
    if positive(ry):
        ltb_inputs["ry"] = (ry, LENGTH, sec.source("A"))
    if rts is not None:
        ltb_inputs["rts"] = (rts, LENGTH, sec.source("rts"))
    elif positive(sec.Iy) and positive(sec.Cw) and positive(Sx):
        rts = rts_f2_7(sec.Iy, sec.Cw, Sx)                          # Eq. F2-7
        ltb_details["rts (Eq. F2-7)"] = (rts, LENGTH)
        ltb_inputs.update(sec.inputs("Iy", "Cw"))
    if ho is not None:
        ltb_inputs["ho"] = (ho, LENGTH, sec.source("ho"))
    elif "d" in sec.dims and "tf" in sec.dims:
        ho = sec.dims["d"] - sec.dims["tf"]
        ltb_details["ho = d - tf (distance between flange centroids)"] = (ho, LENGTH)
    if positive(sec.J):
        ltb_inputs.update(sec.inputs("J"))
    full = positive(rts) and positive(ho) and positive(sec.J) and positive(Sx)

    def capacity_of(Lb, Cb):
        if Lb <= lp_i(ry, E, Fy):
            return None, "F2.2(a)", {}
        Mn, eq, Lp, Lr, Fcr = ltb_i_shape(E, Fy, Mp, Sx, ry, rts, sec.J, 1.0, ho, Lb, Cb)
        more = {"Lp (Eq. F2-5)": (Lp, LENGTH), "Lr (Eq. F2-6)": (Lr, LENGTH), "c (Eq. F2-8a)": (1.0, NONE)}
        if Fcr is not None:
            more["Fcr (Eq. F2-4)"] = (Fcr, STRESS)
        return Mn, eq, more

    def clause_of(eq):
        if not eq:
            return f"{base}{'.1' if base == 'F3' else '.2'}"
        return (("F3.1; " if base == "F3" else "") + f"F2.2, Eq. {eq}, F2-5, F2-6"
                + (", F2-7" if "rts (Eq. F2-7)" in ltb_details else "") + ", F2-8a; F2.1, Eq. F2-1; B4.1, Table B4.1b")

    ltb_caps, ltb_detail = [], {}
    segments = ltb_segments(ctx, mo)
    Lp = lp_i(ry, E, Fy) if positive(ry) else None                  # Eq. F2-5
    # ry is needed for any unbraced length > 0; Sx, rts, ho, J only for a segment longer than Lp
    beyond = [s for s in segments if s.blocked or s.Lb > (Lp or 0.0)]
    ltb_reason = None
    if beyond and Lp is None:
        ltb_reason = f"section {sec.id!r} gives no ry (Eq. F2-5)"
    elif beyond and not full:
        lack = [n for n, v in (("Sx", Sx), ("rts (or Iy and Cw for Eq. F2-7)", rts), ("ho (or d and tf)", ho),
                               ("J", sec.J)) if not positive(v)]
        ltb_reason = (f"section {sec.id!r} gives no {', '.join(lack)}: needed by Eq. F2-4 / F2-6 for an "
                      f"unbraced length beyond Lp = {Lp:.6g}")
    if ltb_reason:
        worst = max(beyond, key=lambda s: s.Mr)
        checks.append(not_checked(NAME_LTB, f"AISC 360-22 {clause_of('')}", ltb_reason, ctx, demand=worst.Mr,
                                  dimension=MOMENT, combination=combination, x=worst.x, inputs=ltb_inputs,
                                  details=ltb_details))
    elif beyond:
        c, ltb_caps, _ = _ltb_checks(ctx, mo, combination, capacity_of, clause_of, ltb_inputs, ltb_details,
                                     sec.notes)
        if c is not None:
            checks.append(c)
    else:
        ltb_detail = _no_ltb_detail(ctx, mo)

    # --- yielding (F2.1) or compression flange local buckling (F3.2)
    if compact:
        checks.insert(0, _strength(NAME_Y["x"], "F2.1, Eq. F2-1; B4.1, Table B4.1b", ctx, Mp, mo.Mrx, mo.x_Mrx,
                                   combination, inputs, {**details, **ltb_detail}, sec.notes))
    elif not positive(Sx):
        checks.insert(0, not_checked(NAME_FLB["x"], "AISC 360-22 F3.2", f"section {sec.id!r} gives no Sx (Eq. "
                                     "F3-1, F3-2)", ctx, inputs=inputs, details=details, **kw))
    else:
        kc = b4.kc_factor(sec.elements["web"].ratio)                # F3.2: kc = 4 / sqrt(h/tw), 0.35 ... 0.76
        Mn, eq = flb_i_major(E, Fy, Mp, Sx, flange.ratio, flange.lambda_p, flange.lambda_r, kc)
        more = {"kc (F3.2)": (kc, NONE)} if eq == "F3-2" else {}
        checks.insert(0, _strength(NAME_FLB["x"], f"F3.2, Eq. {eq}; F2.1, Eq. F2-1; B4.1, Table B4.1b", ctx, Mn,
                                   mo.Mrx, mo.x_Mrx, combination, inputs, {**details, **more, **ltb_detail},
                                   [*sec.notes, f"{flange.cls} flange: Section F3"]))
    return _finish(checks, ltb_caps)


def _i_minor(ctx: DesignContext, mo: Moments, combination: str) -> AxisResult:
    sec, E, Fy = ctx.section, ctx.E, ctx.Fy
    kw = dict(demand=mo.Mry, dimension=MOMENT, combination=combination, x=mo.x_Mry)
    reason = _basic(ctx, "Zy", "Sy")
    if not reason and positive(sec.Ix) and positive(sec.Iy) and not sec.Ix > sec.Iy:
        reason = (f"section {sec.id!r}: the inertia about local z is not larger than about local y; AISC 360-22 "
                  "F6 presumes y = minor axis")
    cls = None
    if not reason:
        cls = b4.classify_flexure(sec, E, Fy, b4.MINOR)
        reason = cls.reason if cls.status != "ok" else None
    if reason:
        return _finish([not_checked(NAME_Y["y"], "AISC 360-22 F6", reason, ctx, **kw)])
    flange = cls.elements[0]
    Mp = mp_minor(Fy, sec.Zy, sec.Sy)                               # Eq. F6-1
    inputs = {**ctx.inputs("E", "Fy"), **sec.inputs("Zy", "Sy", "shape"), **cls.inputs()}
    details = {**cls.details(), "Fy Zy": (Fy * sec.Zy, MOMENT), "1.6 Fy Sy": (1.6 * Fy * sec.Sy, MOMENT),
               "Mp (Eq. F6-1)": (Mp, MOMENT)}
    checks = [_strength(NAME_Y["y"], "F6.1, Eq. F6-1; B4.1, Table B4.1b", ctx, Mp, mo.Mry, mo.x_Mry, combination,
                        inputs, details, sec.notes)]
    if flange.cls != b4.COMPACT:
        Mn, eq = flb_i_minor(E, Fy, Mp, sec.Sy, flange.ratio, flange.lambda_p, flange.lambda_r)
        checks.append(_strength(NAME_FLB["y"], f"F6.2, Eq. {eq}; F6.1, Eq. F6-1; B4.1, Table B4.1b", ctx, Mn,
                                mo.Mry, mo.x_Mry, combination, inputs, details, sec.notes))
    return _finish(checks)


# ------------------------------------------------------------------------------ rectangular HSS / box
def _rect(ctx: DesignContext, mo: Moments, combination: str, axis: str) -> AxisResult:
    sec, E, Fy = ctx.section, ctx.E, ctx.Fy
    major = axis == "x"
    Mr, x = (mo.Mrx, mo.x_Mrx) if major else (mo.Mry, mo.x_Mry)
    kw = dict(demand=Mr, dimension=MOMENT, combination=combination, x=x)
    zname, sname = ("Zx", "Sx") if major else ("Zy", "Sy")
    reason = _basic(ctx, zname)
    cls = None
    if not reason:
        cls = b4.classify_flexure(sec, E, Fy, b4.MAJOR if major else b4.MINOR)
        reason = cls.reason if cls.status != "ok" else None
    if reason:
        return _finish([not_checked(NAME_Y[axis], "AISC 360-22 F7", reason, ctx, **kw)])
    Z, S = getattr(sec, zname), getattr(sec, sname)
    flange, web = cls.elements
    hss = sec.shape == RECT_HSS
    Mp = Fy * Z                                                     # Eq. F7-1
    inputs = {**ctx.inputs("E", "Fy"), **sec.inputs(zname, "shape"), **cls.inputs()}
    if positive(S):
        inputs.update(sec.inputs(sname))
    no_s = f"section {sec.id!r} gives no {sname} (elastic section modulus)"
    details = {**cls.details(), "Mp = Fy Z (Eq. F7-1)": (Mp, MOMENT)}
    notes = [*sec.notes, f"bending about the {_AXIS_TEXT[axis]}"]
    tail = "; F7.1, Eq. F7-1; B4.1, Table B4.1b"
    checks = [_strength(NAME_Y[axis], "F7.1, Eq. F7-1; B4.1, Table B4.1b", ctx, Mp, Mr, x, combination, inputs,
                        details, notes)]

    # --- F7.2 flange local buckling
    if flange.cls == b4.NONCOMPACT and not positive(S):
        checks.append(not_checked(NAME_FLB[axis], "AISC 360-22 F7.2(b)", no_s, ctx, inputs=inputs, details=details,
                                  **kw))
    elif flange.cls == b4.NONCOMPACT:
        Mn = local_buckling_f7(Mp, Fy, S, flange.ratio, flange.lambda_p, flange.lambda_r)       # Eq. F7-2
        checks.append(_strength(NAME_FLB[axis], "F7.2(b), Eq. F7-2" + tail, ctx, Mn, Mr, x, combination, inputs,
                                details, notes))
    elif flange.cls == b4.SLENDER:
        el = flange.element
        t = el.thickness
        Hout = sec.dims.get("H" if major else "B")
        I = sec.Ix if major else sec.Iy
        lack = [n for n, v in (("the wall thickness t", t), ("the outside dimension in the plane of bending", Hout),
                               ("the moment of inertia about the axis of bending", I)) if not positive(v)]
        if lack:
            checks.append(not_checked(NAME_FLB[axis], "AISC 360-22 F7.2(c)", f"slender flange (lambda = "
                                      f"{flange.ratio:.4g} > {flange.lambda_r:.4g}): the effective section modulus "
                                      f"of Eq. F7-3 needs {', '.join(lack)} of section {sec.id!r}", ctx,
                                      inputs=inputs, details=details, **kw))
        else:
            b = el.width if el.width is not None else flange.ratio * t
            be = effective_width_f7(b, t, E, Fy, hss)               # Eq. F7-4 / F7-5
            Se = effective_modulus_both_flanges(I, Hout, t, b, be)
            Mn = Fy * Se                                            # Eq. F7-3
            eq = "F7-4" if hss else "F7-5"
            more = {"b": (b, LENGTH), f"be (Eq. {eq})": (be, LENGTH),
                    "Se (ineffective width removed from both flanges)": (Se, MODULUS)}
            checks.append(_strength(
                NAME_FLB[axis], f"F7.2(c), Eq. F7-3, {eq}" + tail, ctx, Mn, Mr, x, combination,
                {**inputs, "t": (t, LENGTH, sec.source("t")),
                 "outside dimension in the plane of bending": (Hout, LENGTH, sec.source("H"))},
                {**details, **more},
                [*notes, "Se: the Specification defines only be; the ineffective width b - be is removed from "
                         "both flanges without shifting the neutral axis (Commentary F7: slightly conservative)"]))

    # --- F7.3 web local buckling
    if web.cls == b4.NONCOMPACT and not positive(S):
        checks.append(not_checked(NAME_WLB[axis], "AISC 360-22 F7.3(b)", no_s, ctx, inputs=inputs, details=details,
                                  **kw))
    elif web.cls == b4.NONCOMPACT:
        Mn = local_buckling_f7(Mp, Fy, S, web.ratio, web.lambda_p, web.lambda_r)                # Eq. F7-6
        checks.append(_strength(NAME_WLB[axis], "F7.3(b), Eq. F7-6" + tail, ctx, Mn, Mr, x, combination, inputs,
                                details, notes))
    elif web.cls == b4.SLENDER:
        checks.append(not_checked(
            NAME_WLB[axis], "AISC 360-22 F7.3(c)", f"slender web (h/t = {web.ratio:.4g} > {web.lambda_r:.4g}): "
            "Eq. F7-7 (Rpg of Eq. F5-6 with aw = 2 h tw/(b tf)) is not implemented -- Eq. F5-6 is written with "
            "hc/tw and an upper limit on aw whose application to a box is not stated; box sections with slender "
            "webs and slender flanges are not covered by the Specification", ctx, inputs=inputs, details=details,
            **kw))

    # --- F7.4 lateral-torsional buckling
    ltb_caps = []
    Ix, Iy = sec.Ix, sec.Iy
    if not (positive(Ix) and positive(Iy)):
        checks[0].notes.append("LTB (F7.4): the moments of inertia are not given, so the stronger axis is not known")
        if major and Mr > 0.0:
            checks.append(not_checked(NAME_LTB, "AISC 360-22 F7.4", f"section {sec.id!r} gives no Ix and Iy: it "
                                      "cannot be decided whether lateral-torsional buckling applies to bending "
                                      "about local z", ctx, inputs=inputs, details=details, **kw))
    elif major and Ix > Iy:
        ry = sec.ry
        if not (positive(sec.J) and positive(sec.A) and positive(ry)):
            checks.append(not_checked(NAME_LTB, "AISC 360-22 F7.4", f"section {sec.id!r} gives no J, Ag or ry "
                                      "(Eq. F7-10)", ctx, inputs=inputs, details=details, **kw))
        else:
            ltb_inputs = {**inputs, **sec.inputs("J", "A"), "ry": (ry, LENGTH, sec.source("A"))}
            Lp = 0.13 * E * ry * math.sqrt(sec.J * sec.A) / Mp       # Eq. F7-10 (needs no elastic modulus)
            # Sx is needed (Eq. F7-8, F7-11) only for a segment longer than Lp; F7.4(a): no LTB when Lb <= Lp
            segments = ltb_segments(ctx, mo)
            beyond = [] if any(s.blocked for s in segments) else [s for s in segments if s.Lb > Lp]

            def capacity_of(Lb, Cb):
                if Lb <= Lp:                                        # F7.4(a)
                    return None, "F7.4(a)", {"Lp (Eq. F7-10)": (Lp, LENGTH)}
                Mn, eq, _, Lr = ltb_rectangular(E, Fy, Mp, S, ry, sec.J, sec.A, Lb, Cb)
                return Mn, eq, {"Lp (Eq. F7-10)": (Lp, LENGTH), "Lr (Eq. F7-11)": (Lr, LENGTH)}

            def clause_of(eq):
                return "F7.4" + (f", Eq. {eq}, F7-10, F7-11; F7.1, Eq. F7-1" if eq else "")

            if beyond and not positive(S):
                worst = max(beyond, key=lambda s: s.Mr)
                checks.append(not_checked(
                    NAME_LTB, "AISC 360-22 F7.4", f"section {sec.id!r} gives no {sname}: needed by Eq. F7-8 / "
                    f"F7-11 for an unbraced length beyond Lp = {Lp:.6g} (Eq. F7-10)", ctx, demand=worst.Mr,
                    dimension=MOMENT, combination=combination, x=worst.x, inputs=ltb_inputs,
                    details={**details, "Lp (Eq. F7-10)": (Lp, LENGTH)}))
            else:
                c, ltb_caps, _ = _ltb_checks(ctx, mo, combination, capacity_of, clause_of, ltb_inputs, details,
                                             notes)
                if c is not None:
                    checks.append(c)
                else:
                    checks[0].details.update(_no_ltb_detail(ctx, mo))
                    checks[0].details["Lp (Eq. F7-10)"] = (Lp, LENGTH)
    elif not major and Iy > Ix:
        if Mr > 0.0:
            checks.append(not_checked(
                NAME_LTB_MINOR, "AISC 360-22 F7.4", f"section {sec.id!r} is stronger about local y (Iy = {Iy:g} > "
                f"Ix = {Ix:g}): bending about local y is its major-axis bending, for which the model has no LTB "
                "bracing data (bracing.topFlange / bottomFlange refer to the +y / -y flanges); orient the "
                "section with its strong axis as local z", ctx, inputs=inputs, details=details, **kw))
    else:
        checks[0].details["LTB"] = ("does not apply: bending about the axis of smaller or equal inertia "
                                    "(F7.4 User Note)", TEXT)
    return _finish(checks, ltb_caps)


# ------------------------------------------------------------------------------ round HSS
def _round(ctx: DesignContext, mo: Moments, combination: str) -> AxisResult:
    sec, E, Fy = ctx.section, ctx.E, ctx.Fy
    kw = dict(demand=mo.Mr, dimension=MOMENT, combination=combination, x=mo.x_Mr)
    reason = _basic(ctx, "Zx")
    cls = None
    if not reason:
        cls = b4.classify_flexure(sec, E, Fy)
        reason = cls.reason if cls.status != "ok" else None
    if reason:
        return _finish([not_checked(NAME_Y_ROUND, "AISC 360-22 F8", reason, ctx, **kw)])
    wall = cls.elements[0]
    limit = 0.45 * E / Fy                                           # F8: applies for D/t < 0.45 E/Fy
    details = {**cls.details(), "0.45 E/Fy (F8)": (limit, NONE)}
    if not wall.ratio < limit:
        return _finish([not_checked(NAME_Y_ROUND, "AISC 360-22 F8", f"D/t = {wall.ratio:g} is not below 0.45 E/Fy "
                                    f"= {limit:g}: Section F8 does not apply to such round HSS", ctx,
                                    details=details, **kw)])
    Z, S = sec.Zx, sec.Sx
    Mp = Fy * Z                                                     # Eq. F8-1
    inputs = {**ctx.inputs("E", "Fy"), "Z": (Z, MODULUS, sec.source("Zx")), **sec.inputs("shape"), **cls.inputs()}
    if positive(S):
        inputs["S"] = (S, MODULUS, sec.source("Sx"))
    details["Mp = Fy Z (Eq. F8-1)"] = (Mp, MOMENT)
    notes = [*sec.notes, "resultant moment sqrt(Mx^2 + My^2): every axis of a circular section is a principal axis"]
    checks = [_strength(NAME_Y_ROUND, "F8.1, Eq. F8-1; B4.1, Table B4.1b", ctx, Mp, mo.Mr, mo.x_Mr, combination,
                        inputs, details, notes)]
    if wall.cls != b4.COMPACT and not positive(S):
        checks.append(not_checked(NAME_LB_ROUND, "AISC 360-22 F8.2", f"section {sec.id!r} gives no elastic section "
                                  "modulus S (Eq. F8-2, F8-3)", ctx, inputs=inputs, details=details, **kw))
    elif wall.cls != b4.COMPACT:
        Mn, eq = round_local_buckling(E, Fy, S, wall.ratio, wall.cls == b4.SLENDER)
        checks.append(_strength(NAME_LB_ROUND, f"F8.2, Eq. {eq}; B4.1, Table B4.1b", ctx, Mn, mo.Mr, mo.x_Mr,
                                combination, inputs, details, notes))
    return _finish(checks)


# ============================================================================== public
def axis_strength(ctx: DesignContext, demand: Demand, axis: str) -> AxisResult | None:
    """Chapter F checks and the governing available strength for one axis: ``"x"`` (major), ``"y"`` (minor)
    or ``"r"`` (resultant, round sections).  ``None`` when the demand carries no moment data."""
    mo = moments(demand)
    if mo is None:
        return None
    shape, c = ctx.section.shape, demand.combination
    if shape == ROUND_HSS:
        return _round(ctx, mo, c) if axis == "r" else AxisResult([], None, "round section: resultant moment")
    if axis == "r":
        return AxisResult([], None, "not a round section")
    if shape == I_SHAPE:
        return _i_major(ctx, mo, c) if axis == "x" else _i_minor(ctx, mo, c)
    if shape in (RECT_HSS, BOX):
        return _rect(ctx, mo, c, axis)
    Mr, x = (mo.Mrx, mo.x_Mrx) if axis == "x" else (mo.Mry, mo.x_Mry)
    return _finish([not_checked(NAME_Y[axis], CLAUSE, ctx.section.reason or "section outside the implemented "
                                "scope", ctx, demand=Mr, dimension=MOMENT, combination=c, x=x)])


def flexure_major(ctx: DesignContext, demand: Demand) -> list:
    """Major-axis flexure (F2, F3, F7); for a round section the checks of the resultant moment (F8)."""
    res = axis_strength(ctx, demand, "r" if ctx.section.shape == ROUND_HSS else "x")
    return [] if res is None else res.checks


def flexure_minor(ctx: DesignContext, demand: Demand) -> list:
    """Minor-axis flexure (F6, F7); nothing for a round section (covered by the resultant moment)."""
    if ctx.section.shape == ROUND_HSS:
        return []
    res = axis_strength(ctx, demand, "y")
    return [] if res is None else res.checks


LIMIT_STATES = (flexure_major, flexure_minor)


def check(ctx: DesignContext, demand: Demand) -> list:
    """All Chapter F checks of one combination."""
    return [c for fn in LIMIT_STATES for c in fn(ctx, demand)]
