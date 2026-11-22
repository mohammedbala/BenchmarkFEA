"""Generate the SVG figures of docs/CONVENTIONS.md.

    python docs/figures/make_figures.py

Pure standard library; the output is deterministic and committed next to this
script.  The member triads are computed with benchfea.geometry.local_triad, so
the figures cannot drift from the engine.
"""
from __future__ import annotations

import math
from pathlib import Path

from benchfea.geometry import local_triad

OUT = Path(__file__).parent
INK, MUTED, RED, GREEN, BLUE, ORANGE = "#1f2933", "#6b7785", "#c0392b", "#1e8449", "#1f5fa8", "#b9770e"
AXIS_COLOR = {"x": RED, "y": GREEN, "z": BLUE, "X": RED, "Y": GREEN, "Z": BLUE}


class Svg:
    def __init__(self, w, h, title):
        self.w, self.h = w, h
        self.parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            f'font-family="Helvetica, Arial, sans-serif" font-size="13" role="img" aria-label="{title}">',
            f"<title>{title}</title>",
            f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
        ]

    def line(self, p, q, color=INK, width=1.5, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.parts.append(f'<line x1="{p[0]:.1f}" y1="{p[1]:.1f}" x2="{q[0]:.1f}" y2="{q[1]:.1f}" '
                          f'stroke="{color}" stroke-width="{width}"{d}/>')

    def poly(self, pts, fill="#e4e9ef", stroke=INK, width=1.2):
        s = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        self.parts.append(f'<polygon points="{s}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')

    def head(self, p, q, color, size=9):
        a = math.atan2(q[1] - p[1], q[0] - p[0])
        pts = [q, (q[0] - size * math.cos(a - 0.35), q[1] - size * math.sin(a - 0.35)),
               (q[0] - size * math.cos(a + 0.35), q[1] - size * math.sin(a + 0.35))]
        self.poly(pts, fill=color, stroke=color, width=0.5)

    def arrow(self, p, q, color=INK, width=2.0, double=False):
        """Arrow from p to q; double = two heads at q (moment / rotation vector, right-hand rule)."""
        self.line(p, q, color, width)
        self.head(p, q, color)
        if double:
            a = math.atan2(q[1] - p[1], q[0] - p[0])
            self.head(p, (q[0] - 8 * math.cos(a), q[1] - 8 * math.sin(a)), color)

    def arc(self, c, r, a0, a1, color=INK, width=2.0):
        """Circular arrow about c from angle a0 to a1 (degrees, counter-clockwise positive on screen)."""
        def pt(a):
            return c[0] + r * math.cos(math.radians(a)), c[1] - r * math.sin(math.radians(a))
        p, q = pt(a0), pt(a1)
        sweep = 0 if a1 > a0 else 1
        large = 1 if abs(a1 - a0) > 180 else 0
        self.parts.append(f'<path d="M{p[0]:.1f},{p[1]:.1f} A{r},{r} 0 {large} {sweep} {q[0]:.1f},{q[1]:.1f}" '
                          f'fill="none" stroke="{color}" stroke-width="{width}"/>')
        back = pt(a1 - math.copysign(12, a1 - a0))
        self.head(back, q, color)

    def text(self, p, s, color=INK, size=13, anchor="start", bold=False, italic=False):
        style = (' font-weight="bold"' if bold else "") + (' font-style="italic"' if italic else "")
        self.parts.append(f'<text x="{p[0]:.1f}" y="{p[1]:.1f}" fill="{color}" font-size="{size}" '
                          f'text-anchor="{anchor}"{style}>{s}</text>')

    def dot(self, p, r=4, color=INK):
        self.parts.append(f'<circle cx="{p[0]:.1f}" cy="{p[1]:.1f}" r="{r}" fill="{color}"/>')

    def save(self, name):
        (OUT / name).write_text("\n".join(self.parts + ["</svg>"]) + "\n", encoding="utf-8")
        print("wrote", OUT / name)


def projector(origin, scale, vertical="Z"):
    """Oblique projection: global X to the right, the vertical axis up, the third axis receding up-right."""
    ca, sa = 0.5 * math.cos(math.radians(35)), 0.5 * math.sin(math.radians(35))

    def proj(p):
        x, y, z = p
        if vertical == "Y":                     # Y up, Z towards the viewer (right-handed)
            x, y, z = x, -z, y
        return origin[0] + scale * (x + ca * y), origin[1] - scale * (z + sa * y)
    return proj


def add(p, v, s=1.0):
    return tuple(a + s * b for a, b in zip(p, v))


def global_axes():
    g = Svg(880, 330, "Global axes, loads, displacements and reactions")
    for k, (vertical, ox) in enumerate((("Z", 90), ("Y", 520))):
        pr = projector((ox, 250), 120, vertical)
        g.text((ox - 80, 28), f'verticalAxis: "{vertical}"', bold=True, size=14)
        for name, v in (("X", (1, 0, 0)), ("Y", (0, 1, 0)), ("Z", (0, 0, 1))):
            g.arrow(pr((0, 0, 0)), pr(v), AXIS_COLOR[name], 2.5)
            lab = pr(add((0, 0, 0), v, 1.12))
            g.text((lab[0] - 4, lab[1] + 4), name + (" (up)" if name == vertical else ""), AXIS_COLOR[name], bold=True)
        g.dot(pr((0, 0, 0)))
        # a node with its positive force and moment components
        n = (1.25, 0.0, 0.9) if vertical == "Z" else (1.25, 0.9, 0.0)
        g.dot(pr(n), 5)
        g.text((pr(n)[0] - 44, pr(n)[1] - 8), "node", MUTED)
        for name, v in (("X", (1, 0, 0)), ("Y", (0, 1, 0)), ("Z", (0, 0, 1))):
            c = AXIS_COLOR[name]
            g.arrow(pr(n), pr(add(n, v, 0.42)), c)
            g.arrow(pr(add(n, v, 0.5)), pr(add(n, v, 0.8)), c, 1.6, double=True)
            lab = pr(add(n, v, 0.92))
            g.text((lab[0] + 2, lab[1] + 4), f"F{name.lower()}, u{name.lower()} / M{name.lower()}, r{name.lower()}", c, 11)
        grav = pr(add(n, (0, 0, -1) if vertical == "Z" else (0, -1, 0), 0.55))
        g.arrow(pr(n), grav, ORANGE, 2, )
        g.text((grav[0] + 6, grav[1]), "gravity (default)", ORANGE, 11)
    g.text((20, 312), "Single arrow: force / translation, positive along the axis.  Double arrow: moment / rotation "
                      "vector, positive by the right-hand rule.", MUTED, 12)
    g.save("global_axes.svg")


def member_axes():
    g = Svg(880, 450, "Member local axes: default orientation rule")
    pr = projector((60, 350), 85, "Z")
    for name, v in (("X", (1, 0, 0)), ("Y", (0, 1, 0)), ("Z", (0, 0, 1))):
        g.arrow(pr((0, 0, 0)), pr(add((0, 0, 0), v, 0.45)), AXIS_COLOR[name], 1.8)
        lab = pr(add((0, 0, 0), v, 0.56))
        g.text((lab[0] - 4, lab[1] + 4), name, AXIS_COLOR[name], bold=True)
    g.text((20, 392), "global axes (Z vertical)", MUTED, 12)
    cases = [("beam", (1.3, 0, 2.4), (3.9, 0, 2.4), 0.0, "horizontal member: y = +vertical", (170, 62)),
             ("column", (5.0, 0, 0.5), (5.0, 0, 2.9), 0.0, "vertical member: y = global +X", (420, 62)),
             ("inclined", (1.5, 0, 0.3), (3.7, 0, 1.5), 0.0, "inclined member: y in the vertical plane, upward",
              (170, 372)),
             ("rolled", (6.7, 0, 2.4), (8.6, 0, 2.4), 30.0, "roll = +30 deg about +x", (640, 62))]
    for name, a, b, roll, note, note_at in cases:
        R = local_triad(a, b, "Z", roll_deg=roll)
        g.line(pr(a), pr(b), INK, 5)
        g.dot(pr(a), 5)
        g.dot(pr(b), 5)
        g.text((pr(a)[0] - 14, pr(a)[1] + 20), "i", bold=True)
        g.text((pr(b)[0] + 8, pr(b)[1] + 20), "j", bold=True)
        mid = add(a, tuple(q - p for p, q in zip(a, b)), 0.45)
        for k, ax in enumerate("xyz"):
            tip = add(mid, tuple(R[k]), 0.62)
            g.arrow(pr(mid), pr(tip), AXIS_COLOR[ax], 2.2)
            lab = pr(add(mid, tuple(R[k]), 0.74))
            g.text((lab[0] - 3, lab[1] + 4), ax, AXIS_COLOR[ax], bold=True, italic=True)
        g.text(note_at, note, MUTED, 12)
    g.text((20, 420), "local x: from node i to node j.   z = unit(x cross v),  y = z cross x,   with v = the vertical "
                      "axis, or global +X for a vertical member.", MUTED, 12)
    g.text((20, 438), "orientation.referenceVector replaces v;  orientation.roll then rotates y and z about +x.",
           MUTED, 12)
    g.save("member_local_axes.svg")


def section_axes():
    g = Svg(880, 300, "Cross-section axes seen from the member end j looking towards i")
    # view from +x looking back (x points to the viewer): y up, z to the LEFT
    def shape(cx, cy, d, bf, tf, tw, roll=0.0):
        pts = [(-bf / 2, d / 2), (bf / 2, d / 2), (bf / 2, d / 2 - tf), (tw / 2, d / 2 - tf), (tw / 2, -d / 2 + tf),
               (bf / 2, -d / 2 + tf), (bf / 2, -d / 2), (-bf / 2, -d / 2), (-bf / 2, -d / 2 + tf),
               (-tw / 2, -d / 2 + tf), (-tw / 2, d / 2 - tf), (-bf / 2, d / 2 - tf)]
        c, s = math.cos(math.radians(roll)), math.sin(math.radians(roll))
        # screen: horizontal = -z, vertical = y; a positive roll about +x (towards the viewer) is counter-clockwise
        g.poly([(cx + (h * c - v * s), cy - (h * s + v * c)) for h, v in pts])
        return c, s

    for cx, roll, label in ((150, 0.0, "roll = 0"), (430, 30.0, "roll = +30 deg")):
        c, s = shape(cx, 150, 150, 100, 12, 8, roll)
        g.arrow((cx, 150), (cx - 95 * s, 150 - 95 * c), GREEN, 2.2)
        g.text((cx - 95 * s - 6, 150 - 95 * c - 8), "y", GREEN, bold=True, italic=True)
        g.arrow((cx, 150), (cx - 95 * c, 150 + 95 * s), BLUE, 2.2)
        g.text((cx - 95 * c - 14, 150 + 95 * s + 4), "z", BLUE, bold=True, italic=True)
        g.dot((cx, 150), 4, RED)
        g.text((cx + 8, 150 + 16), "x (towards viewer)", RED, 11, italic=True)
        g.text((cx - 40, 280), label, bold=True)
        if roll:
            g.arc((cx, 150), 62, 90, 90 + roll, ORANGE)
            g.text((cx - 10, 78), "roll", ORANGE, 12)
    for k, s in enumerate(["y: depth direction (web), d / H / h",
                           "z: width direction (flanges), bf / B / b",
                           "Iz, Sz, Zz: bending about z = AISC x-x (major)",
                           "Iy, Sy, Zy: bending about y = AISC y-y (minor)",
                           "top flange = flange on the +y side",
                           "Asy: shear along y (web);  Asz: shear along z"]):
        g.text((560, 70 + 24 * k), s, INK if k < 2 else MUTED, 12)
    g.save("section_axes.svg")


def member_forces():
    g = Svg(880, 580, "Member end forces and internal forces")
    for row, (plane, v, m, mcw, bend) in enumerate((
            ("local x-y plane (z towards the viewer)", "y", "Mz", False, "Mz &gt; 0: tension on the -y side"),
            ("local x-z plane (y away from the viewer)", "z", "My", True, "My &gt; 0: tension on the +z side"))):
        y0 = 150 + row * 270
        col = AXIS_COLOR[v]
        g.text((20, y0 - 112), plane, bold=True, size=14)
        # axes
        g.arrow((30, y0 + 70), (75, y0 + 70), RED, 1.8)
        g.text((80, y0 + 74), "x", RED, bold=True, italic=True)
        g.arrow((30, y0 + 70), (30, y0 + 25), col, 1.8)
        g.text((24, y0 + 18), v, col, bold=True, italic=True)
        # (a) member end forces: on the member, along the positive local axes at BOTH ends
        xa, xb = 130, 330
        g.poly([(xa, y0 - 9), (xb, y0 - 9), (xb, y0 + 9), (xa, y0 + 9)])
        g.text((xa - 6, y0 + 34), "i", bold=True)
        g.text((xb + 2, y0 + 34), "j", bold=True)
        for x, sgn, end in ((xa, -1, "i"), (xb, 1, "j")):
            x1 = x + sgn * 8
            g.arrow((x1 - 46 if sgn < 0 else x1, y0), (x1 if sgn < 0 else x1 + 46, y0), RED)
            g.text(((x1 - 62) if sgn < 0 else x1 + 50, y0 - 6 if sgn < 0 else y0 + 4), "Fx", RED, 12,
                   anchor="middle" if sgn < 0 else "start")
            g.arrow((x, y0 - 12), (x, y0 - 62), col)
            g.text((x + 6, y0 - 56), f"F{v}", col, 12)
            ccw = not mcw                 # positive moment appears counter-clockwise unless the axis recedes
            if sgn > 0:
                g.arc((x + 4, y0), 30, -60 if ccw else 60, 60 if ccw else -60, col)
            else:
                g.arc((x - 4, y0), 30, 120 if ccw else 240, 240 if ccw else 120, col)
            g.text((x + sgn * 40 - (10 if sgn < 0 else 0), y0 + 44), m, col, 12, anchor="end" if sgn < 0 else "start")
        g.text((xa - 40, y0 + 78), "(a) member end forces: joint on member, +local axes at both ends", MUTED, 12)
        # (b) internal forces on a cut segment: positive on the +x face
        xc, xd = 560, 670
        g.poly([(xc, y0 - 9), (xd, y0 - 9), (xd, y0 + 9), (xc, y0 + 9)])
        g.arrow((xd + 8, y0), (xd + 54, y0), RED)
        g.text((xd + 58, y0 + 4), "N", RED, 12)
        g.arrow((xd, y0 - 12), (xd, y0 - 62), col)
        g.text((xd + 6, y0 - 56), f"V{v}", col, 12)
        g.arrow((xc - 8, y0), (xc - 54, y0), RED)
        g.text((xc - 68, y0 + 4), "N", RED, 12)
        g.arrow((xc, y0 + 12), (xc, y0 + 62), col)
        g.text((xc - 22, y0 + 60), f"V{v}", col, 12)
        if mcw:
            g.arc((xd + 4, y0), 30, 60, -60, col)
            g.arc((xc - 4, y0), 30, 120, 240, col)
        else:
            g.arc((xd + 4, y0), 30, -60, 60, col)
            g.arc((xc - 4, y0), 30, 240, 120, col)
        g.text((xd + 40, y0 + 44), m, col, 12)
        g.text((xc - 60, y0 + 44 if not mcw else y0 - 30), m, col, 12)
        g.text((xd - 12, y0 + 30), "+x face", MUTED, 11)
        g.text((xc - 70, y0 + 92), "(b) internal forces on a cut piece, positive as drawn", MUTED, 12)
        g.text((xc - 70, y0 + 110), bend, MUTED, 12)
    g.text((20, 566), "Station forces at x = L equal the end forces at j; at x = 0 they are minus the end forces at i.  "
                      "T and Mx follow the right-hand rule about +x.", MUTED, 12)
    g.save("member_forces.svg")


if __name__ == "__main__":
    global_axes()
    member_axes()
    section_axes()
    member_forces()
