"""Unstable models of the independent review: each builder returns a raw input that must be refused
with a diagnosis naming what moves."""
from __future__ import annotations

from harness import case, member, nodes
from models import start

UNSTABLE = {}


def unstable(f):
    UNSTABLE[f.__name__] = f
    return f


@unstable
def U01_supports_leave_a_rigid_translation():
    d = start()
    nodes(d, A=(0, 0, 0), B=(5, 0, 0), C=(5, 4, 0))
    d["supports"] = {n: {"type": "roller"} for n in "ABC"}          # vertical only: slides in X, Y, spins about Z
    member(d, "m1", "A", "B")
    member(d, "m2", "B", "C")
    case(d, "L", {"type": "nodal", "node": "B", "Fz": -10})
    return d


@unstable
def U02_plane_frame_on_two_pins_rotates_about_the_support_line():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 4), C=(6, 0, 4), D=(6, 0, 0))
    d["supports"] = {"A": {"type": "pinned"}, "D": {"type": "pinned"}}
    member(d, "c1", "A", "B")
    member(d, "g", "B", "C")
    member(d, "c2", "D", "C")
    case(d, "L", {"type": "nodal", "node": "B", "Fx": 10})
    return d


@unstable
def U03_joint_with_all_member_ends_released():
    d = start()
    nodes(d, A=(0, 0, 0), B=(4, 0, 0), C=(8, 0, 0), D=(4, 0, -3))
    d["supports"] = {"A": {"type": "fixed"}, "C": {"type": "fixed"}, "D": {"type": "fixed"}}
    member(d, "m1", "A", "B", releases={"j": {"ry": True, "rz": True, "rx": True}})
    member(d, "m2", "B", "C", releases={"i": {"ry": True, "rz": True, "rx": True}})
    member(d, "m3", "D", "B", releases={"j": {"ry": True, "rz": True, "rx": True}})
    case(d, "L", {"type": "nodal", "node": "B", "Fz": -10})
    return d


@unstable
def U04_collinear_pin_ended_members():
    d = start()
    nodes(d, A=(0, 0, 0), B=(3, 1.5, 0.6), C=(6, 3, 1.2))
    d["supports"] = {"A": {"type": "pinned"}, "C": {"type": "pinned"}}
    member(d, "t1", "A", "B", type="truss")
    member(d, "t2", "B", "C", type="truss")
    case(d, "L", {"type": "nodal", "node": "B", "Fz": -10})
    return d


@unstable
def U05_member_spins_about_its_own_axis():
    d = start()
    nodes(d, A=(0, 0, 0), B=(3, 2, 4), C=(7, 2, 4), D=(7, 2, 0))
    d["supports"] = {"A": {"type": "pinned"}, "D": {"type": "fixed"}}
    member(d, "strut", "A", "B", releases={"j": {"rx": True}})          # ball joint below, twist released above
    member(d, "g", "B", "C")
    member(d, "c", "D", "C")
    case(d, "L", {"type": "nodal", "node": "B", "Fz": -10})
    return d


@unstable
def U06_storey_held_by_the_diaphragm_only():
    d = start()
    pts = [(0, 0), (6, 0), (6, 5), (0, 5)]
    for k, (x, y) in enumerate(pts):
        nodes(d, **{f"g{k}": (x, y, 0), f"u{k}": (x, y, 3.5)})
        d["supports"][f"g{k}"] = {"type": "pinned", "dofs": {"rz": "fixed"}}
        member(d, f"c{k}", f"g{k}", f"u{k}", releases={"j": {"ry": True, "rz": True}})
    for k in range(4):
        member(d, f"b{k}", f"u{k}", f"u{(k + 1) % 4}")
    d["constraints"] = {"floor": {"type": "rigidDiaphragm", "master": "u0", "slaves": ["u1", "u2", "u3"]}}
    case(d, "L", {"type": "nodal", "node": "u0", "Fx": 10})
    return d


@unstable
def U07_pin_jointed_portal_sways():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 4), C=(6, 0, 4), D=(6, 0, 0))
    d["supports"] = {"A": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}},
                     "D": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}}}
    member(d, "c1", "A", "B")
    member(d, "g", "B", "C", releases={"i": {"rz": True}, "j": {"rz": True}})
    member(d, "c2", "D", "C")
    case(d, "L", {"type": "nodal", "node": "B", "Fx": 10})
    return d


@unstable
def U08_second_part_not_connected_to_the_supported_part():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 4), P=(5, 0, 0), Q=(5, 0, 3), R=(8, 1, 3))
    d["supports"] = {"A": {"type": "fixed"}}
    member(d, "c", "A", "B")
    member(d, "f1", "P", "Q")
    member(d, "f2", "Q", "R")
    case(d, "L", {"type": "nodal", "node": "B", "Fx": 10})
    return d


@unstable
def U09_tension_only_bracing_goes_slack():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 4), C=(5, 0, 4), D=(5, 0, 0))
    d["supports"] = {"A": {"type": "pinned"}, "D": {"type": "pinned"}}
    for n in "BC":
        d["supports"][n] = {"dofs": {"uy": "fixed"}}
    member(d, "c1", "A", "B", type="truss")
    member(d, "g", "B", "C", type="truss")
    member(d, "c2", "D", "C", type="truss")
    member(d, "x1", "A", "C", type="tensionOnly")                      # only one diagonal: works for +X only
    case(d, "L", {"type": "nodal", "node": "B", "Fx": -10})
    return d


@unstable
def U10_gravity_only_on_tension_only_x_bracing_is_not_determined():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 4), C=(5, 0, 4), D=(5, 0, 0))
    d["supports"] = {"A": {"type": "pinned"}, "D": {"type": "pinned"}}
    for n in "BC":
        d["supports"][n] = {"dofs": {"uy": "fixed"}}
    member(d, "c1", "A", "B", type="truss")
    member(d, "g", "B", "C", type="truss")
    member(d, "c2", "D", "C", type="truss")
    member(d, "x1", "A", "C", type="tensionOnly")
    member(d, "x2", "D", "B", type="tensionOnly")
    case(d, "L", {"type": "nodal", "node": "B", "Fz": -10}, {"type": "nodal", "node": "C", "Fz": -10})
    return d


@unstable
def U11_beam_free_to_slide_along_its_axis():
    d = start()
    nodes(d, A=(0, 0, 0), B=(4, 0, 0), C=(9, 0, 0))
    d["supports"] = {"A": {"dofs": {"uy": "fixed", "uz": "fixed", "rx": "fixed"}}, "B": {"dofs": {"uy": "fixed", "uz": "fixed"}},
                     "C": {"dofs": {"uy": "fixed", "uz": "fixed"}}}
    member(d, "s1", "A", "B")
    member(d, "s2", "B", "C")
    case(d, "L", {"type": "nodal", "node": "B", "Fz": -10})
    return d


@unstable
def U12_member_release_mechanism():
    d = start()
    nodes(d, A=(0, 0, 0), B=(4, 0, 0))
    d["supports"] = {"A": {"type": "fixed"}, "B": {"type": "fixed"}}
    member(d, "m", "A", "B", releases={"i": {"rz": True}, "j": {"rz": True, "uy": True}})
    case(d, "L", {"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -3})
    return d
